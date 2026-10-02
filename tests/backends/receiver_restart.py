#!/usr/bin/env python3
"""On-demand receiver/PDE checkpoint equivalence, all quantities and three formats."""
import argparse
import json
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np
from receiver_formats import read_sac
from validation import check_solver_log

MODES=('displacement','velocity','div','curl','strains','displacementgradient')


def receivers():
    lines=[]
    for mode in MODES:
        for geo in ((0,1) if mode in ('displacement','velocity') else (0,)):
            for kind in ('sac','text','hdf1','hdf3'):
                name=f'{mode}-geo{geo}-{kind}'
                options={'sac':'sacformat=1 usgsformat=1',
                         'text':'sacformat=0 usgsformat=1',
                         'hdf1':f'sacformat=0 usgsformat=0 hdf5format=1 hdf5file={name}.h5 downSample=1',
                         'hdf3':f'sacformat=0 usgsformat=0 hdf5format=1 hdf5file={name}.h5 downSample=3'}[kind]
                lines.append(f'rec x=1800 y=1800 z=1600 file={name} sta=station variables={mode} nsew={geo} writeEvery=7 {options}\n')
    return ''.join(lines)


def capture(root):
    traces={}
    for path in root.iterdir():
        if path.suffix=='.h5':
            with h5py.File(path) as stream:
                group=stream['station'];n=int(group['NPTS'][0])
                for key in group:
                    if group[key].shape and group[key].shape[0]>=n and not key.endswith(('CMPAZ','CMPINC')):
                        traces[f'{path.name}/{key}']=group[key][:n].astype(float)
        elif path.suffix=='.txt':
            data=np.loadtxt(path)
            for column in range(data.shape[1]): traces[f'{path.name}/{column}']=data[:,column]
        elif path.name.split('-')[0] in MODES and path.suffix not in ('.bak','.in','.sw4checkpoint','.log'):
            _,data,header,utc=read_sac(path)
            traces[path.name]=data
            traces[path.name+'/times']=float(header[5])+np.arange(len(data))*float(header[0])
    return traces


def state(path):
    arrays={}
    with h5py.File(path) as stream:
        def read(name,obj):
            if isinstance(obj,h5py.Dataset): arrays[name]=np.asarray(obj[()])
        stream.visititems(read)
    if not arrays: raise ValueError('No checkpoint state to compare')
    return arrays


def compare(reference,actual,label,storage):
    if set(reference)!=set(actual):raise ValueError(f'{label}: datasets/components differ')
    error=0.
    signal=0.
    for name,want in reference.items():
        got=actual[name]
        if want.shape!=got.shape or not np.all(np.isfinite(want)) or not np.all(np.isfinite(got)):
            raise ValueError(f'{label}: invalid {name}')
        peak=float(np.max(np.abs(want))) if want.size else 0.
        delta=float(np.max(np.abs(got-want))) if want.size else 0.
        tolerance=(1e-12+2e-6*peak) if storage else (1e-12+2e-5*peak)
        if delta>tolerance:raise ValueError(f'{label}: {name} error={delta} limit={tolerance}')
        error=max(error,delta);signal=max(signal,peak)
    if not signal:raise ValueError(f'{label}: no nonzero data')
    return error


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sw4',type=Path,required=True)
    parser.add_argument('--backend',choices=('OPENMP','CUDA','HIP'),default='OPENMP')
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--tasks',type=int,default=2)
    parser.add_argument('--case',action='append')
    args=parser.parse_args()
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_') or k=='SLURM_JOB_ID'}
    env.update(OMP_NUM_THREADS='1',MPICH_GPU_SUPPORT_ENABLED='0')
    command=['srun','--exclusive','--exact','-N','1','-n',str(args.tasks),'-c','1']
    command+=['--gres=none'] if args.backend=='OPENMP' else ['--gpus-per-task=1','--gpu-bind=single:1']
    command += [str(args.sw4.resolve(strict=True))]
    report={}
    for az,att,ref in ((0,False,False),(27,False,False),(27,True,False),(27,False,True)):
        name=f'az{az}-att{int(att)}-ref{int(ref)}'
        if args.case and name not in args.case:continue
        work=root/name;work.mkdir(exist_ok=True)
        common=f'grid h=100 x=4000 y=4000 z=4000 lat=37 lon=-122 az={az}\n'
        common+='time steps=48 utcstart=10/02/2026:01:02:03.234567\nfileio path=.\nsupergrid gp=8\n'
        if ref:common+='refinement zmax=2000\ndeveloper ctol=1e-10 cmaxit=400 crelax=0.2\n'
        if att:common+='attenuation nmech=3\n'
        common+='block vp=6000 vs=3464 rho=2700'+(' qp=200 qs=100' if att else '')+'\n'
        common+='source x=1637 y=1643 z=1000 mxy=1e15 mxz=2e15 myy=5e14 type=Gaussian freq=10 t0=0.3\n'
        common+=receivers()+'checkpoint cycleInterval=24 restartpath=. file=restart hdf5=yes'
        def launch(label,restart=None,steps=48):
            path=work/f'{label}.in'
            path.write_text(common.replace('time steps=48',f'time steps={steps}')+(f' restartfile={restart}' if restart else '')+'\n')
            with (work/f'{label}.log').open('w') as log:
                result=subprocess.run(command+[path.name],cwd=work,env=env,stdout=log,
                    stderr=subprocess.STDOUT,timeout=600)
            if result.returncode:raise RuntimeError(f'Solver failed: {work}/{label}.log')
            check_solver_log(work/f'{label}.log')
        launch('uninterrupted')
        baseline=capture(work);pde=state(work/'restart.cycle=48.sw4checkpoint')
        for repetition,cycle in enumerate((24,24),1):
            launch(f'restart{cycle}-{repetition}',f'restart.cycle={cycle:02d}.sw4checkpoint')
            wave_error=compare(baseline,capture(work),name,True)
            state_error=compare(pde,state(work/'restart.cycle=48.sw4checkpoint'),name,False)
            report[f'{name}/cycle{cycle}/repeat{repetition}']={'wave_max_difference':wave_error,'state_max_difference':state_error}
            (root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
            print('PASS:',name,'restart',cycle,report[f'{name}/cycle{cycle}/repeat{repetition}'],flush=True)


if __name__=='__main__':main()
