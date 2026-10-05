#!/usr/bin/env python3
"""On-demand physical-source, energy and failure checks for curved refinement.

Use an existing Slurm allocation. Retains inputs, residual diagnostics and
receivers; --compare-to replays per-component comparisons against another run.
Not registered with pytest or CTest.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np
from compare_waveforms import stations, compare
from validation import check_solver_log


def physical(h=100, duration=.8, extra=False, width=1, attenuation=False,
             gradient=False, cap=None, moment='1e15'):
    text=f'grid h={h} x=4000 y=4000 z=4000 lat=37 lon=-122 az=0 extrapolate={width}\n'
    text+=f'time t={duration}\nsupergrid gp=8\nfileio path=output verbose=3\n'
    text+='developer cfl='+('0.8' if gradient else '0.3')
    if cap is not None:text+=f' curvimaxit={cap}'
    text+='\ntopography input=gaussian zmax=2000 order=4 gaussianAmp=100 gaussianXc=2000 gaussianYc=2000 gaussianLx=1500 gaussianLy=1500\n'
    # Keep at least twelve vertical points in each curved layer for the SBP
    # closures, including the coarse layer below the deeper curved interface.
    text+='refinement '+('zmax=300 zmax=900' if extra else 'zmax=800')+'\n'
    if attenuation:text+='attenuation nmech=1\n'
    text+='block vp=6000 vs=3464 rho=2700'
    if gradient:text+=' vpgrad=0.1 vsgrad=0.05 rhograd=0.02'
    if attenuation:text+=' qp=200 qs=100'
    text+='\nsource x=1637 y=1643 z=1000 mxy='+moment
    text+=' type=Gaussian freq='+('20' if gradient else '8')+' t0=0.05\n'
    return text+'rechdf5 infile=stations.h5 outfile=receivers.h5 writeEvery=3\n'


def fixtures(repo):
    cases={'physical':physical(), 'long-duration':physical(duration=1.6),
           'resolution80':physical(h=80),
           'resolution66':physical(h='66.66666666666667'),
           'multiple-levels':physical(extra=True)}
    for width in range(1,6):
        for attenuation in (False,True):
            cases[f'gradient-extrap{width}-att{int(attenuation)}']=physical(
                width=width,attenuation=attenuation,gradient=True)
    energy=(repo/'pytest/reference/curvimeshrefine/energy-1.in').read_text()
    cases['energy']=energy.replace('verbose=1 path=energy-1','verbose=3 path=output')
    cases.update({f'cap{limit}':physical(cap=limit) for limit in (1,30)})
    for label,value in (('zero','0'),('negative','-1'),('fraction','1.5'),
                        ('text','abc'),('overflow','999999999999999999999999')):
        cases[f'invalid-{label}']=physical(cap=value)
    # Inject a nonfinite source only in this deliberately invalid fixture. It
    # must fail collectively at the interface, not be lost by a MAX reduction.
    cases['nonfinite']=physical(moment='nan')
    return cases


def residuals(log, energy=False):
    text=log.read_text()
    rows=re.findall(r'Curvilinear interface: grid=(\d+) time=(\S+) iterations=(\d+) '
                    r'relative_residual=(\S+) scaled_residual=(\S+)',text)
    if not rows:raise ValueError(f'No interface residual evidence: {log}')
    values=np.asarray([[float(t),int(it),float(relative),float(scaled)]
                       for _,t,it,relative,scaled in rows])
    if not np.all(np.isfinite(values)) or np.any(values[:,1]>500):
        raise ValueError(f'Invalid interface residuals: {log}')
    if np.any((values[:,2]>1e-6)&(values[:,3]>1e-6)):
        raise ValueError(f'Unconverged interface solve: {log}')
    if np.max(values[:,0])<(99.9 if energy else .79):
        raise ValueError(f'Interface run ends before the required arrival: {log}')
    return {'solves':len(rows),'grids':sorted({int(row[0]) for row in rows}),
            'max_iterations':int(np.max(values[:,1])),
            'max_relative_residual':float(np.max(values[:,2])),
            'max_scaled_residual':float(np.max(values[:,3]))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sw4',type=Path,required=True)
    parser.add_argument('--backend',choices=('OPENMP','CUDA','HIP'),default='OPENMP')
    parser.add_argument('--tasks',type=int,default=4)
    parser.add_argument('--nodes',type=int,default=1)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--case',action='append')
    parser.add_argument('--compare-to',type=Path)
    args=parser.parse_args()
    if args.tasks<1 or args.nodes<1:parser.error('Positive tasks/nodes required')
    repo=Path(__file__).resolve().parents[2]
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_') or k=='SLURM_JOB_ID'}
    env.update(OMP_NUM_THREADS='1',OMP_PLACES='threads',OMP_PROC_BIND='close',
               MPICH_GPU_SUPPORT_ENABLED='0',DISABLE_LUSTRE_STRIPE='1')
    command=['srun','--exclusive','--exact','--cpu-bind=threads','-N',str(args.nodes),
             '-n',str(args.tasks),'-c','1']
    command+=['--gres=none'] if args.backend=='OPENMP' else ['--gpus-per-task=1','--gpu-bind=single:1']
    command+=[str(args.sw4.resolve(strict=True)),'case.in']
    cases=fixtures(repo)
    if args.case and set(args.case)-cases.keys():parser.error('Unknown case')
    report={}
    for name,text in cases.items():
        if args.case and name not in args.case:continue
        work=root/name;work.mkdir(exist_ok=False)
        stations(work/'stations.h5', [('station',(1800,1800,600)),('edge',(2800,2600,1400))])
        (work/'case.in').write_text(text)
        log=work/'run.log'
        with log.open('w') as stream:
            result=subprocess.run(command,cwd=work,env=env,stdout=stream,
                                  stderr=subprocess.STDOUT,timeout=1800)
        output=log.read_text()
        if name.startswith('invalid-'):
            if not result.returncode or 'curvimaxit must be a positive integer' not in output or 'Begin time stepping' in output:
                raise ValueError(f'Iteration input not rejected before stepping: {log}')
            record={'rejected_before_stepping':True}
        elif name.startswith('cap') or name=='nonfinite':
            failure=re.search(r'ERROR: no convergence in curvilinear interface at time \S+ after (\d+) iterations; residual=(\S+)',output)
            if not result.returncode or not failure or 'program sw4 finished!' in output:
                raise ValueError(f'Invalid interface was not fatal: {log}')
            if name.startswith('cap') and int(failure[1])!=int(name[3:]):
                raise ValueError(f'Incorrect iteration limit: {log}')
            if name=='nonfinite' and np.isfinite(float(failure[2])):
                raise ValueError(f'Nonfinite residual not detected: {log}')
            record={'failed_collectively':True,'iterations':int(failure[1]),
                    'failure_residual':failure[2]}
        else:
            if result.returncode or 'program sw4 finished!' not in output:
                raise ValueError(f'Curvilinear simulation failed: {log}')
            check_solver_log(log)
            record=residuals(log,name=='energy')
            if name=='energy':
                values=np.loadtxt(work/'output/energy.log')
                if values.ndim!=1 or len(values)<500 or not np.all(np.isfinite(values)) or values[2]<=0:
                    raise ValueError(f'Invalid energy history: {work}')
                drift=float(np.max(np.abs(values[2:]-values[2]))/values[2])
                # The existing SW4 energy test bounds change per timestep.
                # Check both signs and retain cumulative drift separately: the
                # unchanged 1e-6 interface tolerance permits accumulated error.
                step_change=float(np.max(np.abs(np.diff(values[2:])))/values[2])
                if step_change>1e-10:
                    raise ValueError(f'Curved-grid energy change per step: {step_change}')
                record.update(energy_steps=len(values),relative_energy_drift=drift,
                              max_relative_step_change=step_change)
            else:
                path=work/'output/receivers.h5'
                record['waveforms']=compare(path,path,2e-5,1e-12,arrival_after=.5)
                if args.compare_to:
                    record['reference_comparison']=compare(args.compare_to/name/'output/receivers.h5',
                                                           path,2e-5,1e-12,arrival_after=.5)
            if name=='physical' and record['max_iterations']<=31:
                raise ValueError('Physical fixture did not exercise the former iteration cap')
            if name=='multiple-levels' and len(record['grids'])<2:
                raise ValueError('Multiple-level fixture did not solve two curved interfaces')
        report[name]=record
        (root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        print('PASS:',name,record,flush=True)


if __name__=='__main__':main()
