#!/usr/bin/env python3
"""On-demand input-deck observation, parallel-event and endian regressions.

Run inside a CPU allocation. Uses real forward and native inversion executables,
including the observation parser, rather than constructing receivers directly.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

import h5py
import numpy as np
from receiver_formats import read_sac
from validation import check_solver_log


def launch(executable, work, text, tasks, failure=None, label="run"):
    work.mkdir(parents=True, exist_ok=True)
    (work/f'{label}.in').write_text(text)
    env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_') or k=='SLURM_JOB_ID'}
    env.update(OMP_NUM_THREADS='2', MPICH_GPU_SUPPORT_ENABLED='0', DISABLE_LUSTRE_STRIPE='1')
    command=['srun','--exclusive','--exact','--cpu-bind=none','-N','1','-n',str(tasks),'-c','2']
    command+=['--gres=none']
    command+=[str(executable),f'{label}.in']
    with (work/f'{label}.log').open('w') as stream:
        result=subprocess.run(command,cwd=work,env=env,stdout=stream,
                              stderr=subprocess.STDOUT,timeout=600)
    log=(work/f'{label}.log').read_text()
    if failure:
        if result.returncode==0 or failure not in log:
            raise ValueError(f'Expected {failure!r} failure in {work}/run.log')
        return log
    if result.returncode:raise ValueError(f'Solver failed: {work}/run.log')
    check_solver_log(work/f'{label}.log')
    return log


def close(want,got,label,rtol=2e-6):
    want=np.asarray(want);got=np.asarray(got)
    if want.shape!=got.shape or not np.all(np.isfinite(want)) or not np.all(np.isfinite(got)):
        raise ValueError(f'{label}: invalid shapes/data')
    peak=float(np.max(np.abs(want)))
    error=float(np.max(np.abs(want-got)))
    if error>1e-12+rtol*peak:raise ValueError(f'{label}: error {error}, peak {peak}')
    return error


def observations(work):
    values=[]
    for event in (0,1):
        files=list((work/f'event{event}').glob('*_obs.txt'))
        if len(files)!=1:raise ValueError(f'Expected one observation for event {event}: {files}')
        values.append(np.loadtxt(files[0]))
    return values


def close_traces(want,got,label):
    if want.ndim!=2 or got.shape!=want.shape or want.shape[1]!=4:
        raise ValueError(f'{label}: invalid trace shape')
    # Do not let time-axis magnitude hide errors in a weak waveform component.
    close(want[:,0],got[:,0],label+' time')
    return max(close(want[:,c],got[:,c],label+f' component {c}') for c in (1,2,3))


def gradient(work):
    with h5py.File(work/'event0/dfm.h5') as stream:
        value=stream['dfm'][()]
    if not np.any(value):raise ValueError('Zero inversion gradient')
    # Scale derivatives with the declared parameter scales before comparison;
    # density and stiffness derivatives have very different physical units.
    if value.shape[-1]%3:raise ValueError('Invalid gradient component shape')
    return value*np.tile([2700.,3.24e10,3.24e10],value.shape[-1]//3)


def close_gradients(want,got,label):
    if want.shape!=got.shape:raise ValueError(f'{label}: gradient shapes differ')
    return max(close(want[...,c::3],got[...,c::3],label+f' parameter {c}',2e-5)
               for c in range(3))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sw4',type=Path,required=True)
    parser.add_argument('--mopt',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--case',action='append')
    args=parser.parse_args()
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    sw4=args.sw4.resolve(strict=True);mopt=args.mopt.resolve(strict=True)
    base='grid h=100 x=4000 y=4000 z=4000 lat=37 lon=-122 az=27\n'
    base+='time t=0.5 utcstart=10/02/2026:01:02:03.234567\nsupergrid gp=8\ndeveloper cfl=0.8\n'
    base+='block vp=6000 vs=3464 rho=2700\n'
    source='source x=1600 y=1600 z=600 mxy=1e15 t0=0.05 freq=20 type=Gaussian'
    fixture=root/'fixture'
    launch(sw4,fixture,base+'fileio path=.\n'+source+'\n'+
           'rec x=1800 y=1800 z=200 file=station sta=station nsew=0 sacformat=1 usgsformat=1\n',2)
    original=np.loadtxt(fixture/'station.txt')
    files={suffix:(fixture/f'station.{suffix}').read_bytes() for suffix in ('x','y','z')}
    report={}

    def record(name,result):
        report[name]=result
        (root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        print('PASS:',name,result,flush=True)

    def selected(name):return not args.case or name in args.case

    def inversion_input(work,kind,parallel=False,absolute=False,coords=True,basename='station'):
        optimizer=work/'optimizer';optimizer.mkdir(parents=True,exist_ok=True)
        text=base+f'fileio path={optimizer}\n'
        for event in (0,1):
            output=work/f'event{event}';output.mkdir(parents=True,exist_ok=True)
            text+=f'event name=event{event} path={output} obspath={root}/observed{event} parallel={"yes" if parallel else "no"}\n'
            text+=f'time t=0.5 utcstart=10/02/2026:01:02:03.234567 event=event{event}\n'
            text+=source+f' event=event{event}\n'
            position='x=1800 y=1800 ' if coords else ''
            prefix=str(root/f'observed{event}'/basename) if absolute else basename
            options=f'file={prefix}' if kind=='usgs' else f'sacfile1={prefix}.x sacfile2={prefix}.y sacfile3={prefix}.z'
            text+=f'observation {position}z=200 {options} event=event{event} sta=station\n'
        # No optimization update: evaluate the actual initial objective/gradient.
        text+='mparcart nx=4 ny=4 nz=4 init=0\n'
        text+='mrun task=minvert win_mode=none tsoutput=off pseudohessian=off writedfm=on\n'
        text+='nlcg maxit=1 tolerance=1e100\n'
        text+='mscalefactors rho=2700 mu=3.24e10 lambda=3.24e10 misfit=1\n'
        return text

    def write_observed(endian='<',invalid=None):
        for event in (0,1):
            target=root/f'observed{event}';target.mkdir(exist_ok=True)
            factor=1.25+event
            lines=(fixture/'station.txt').read_text().splitlines()
            data=original.copy();data[:,1:]*=factor
            header='\n'.join(line for line in lines if line.startswith('#'))+'\n'
            with (target/'station.txt').open('w') as stream:
                stream.write(header);np.savetxt(stream,data,fmt='%.17g')
            for suffix,raw in files.items():
                floats=np.frombuffer(raw,dtype='<f4',count=70).copy()
                ints=np.frombuffer(raw,dtype='<i4',count=40,offset=280).copy()
                samples=np.frombuffer(raw,dtype='<f4',offset=632).copy()*factor
                if invalid=='undefined':floats[31:33]=-12345
                if invalid=='nan':floats[31]=np.nan
                if invalid=='range':floats[31]=100
                binary=floats.astype(endian+'f4').tobytes()+ints.astype(endian+'i4').tobytes()+raw[440:632]+samples.astype(endian+'f4').tobytes()
                if invalid=='truncated':binary=binary[:200]
                (target/f'station.{suffix}').write_bytes(binary)

    for kind in ('usgs','sac'):
        name=f'events-{kind}'
        if not selected(name):continue
        write_observed()
        serial=root/f'{name}-serial';parallel=root/f'{name}-parallel';absolute=root/f'{name}-absolute'
        logs=[]
        for work,split,full in ((serial,False,False),(parallel,True,False),(absolute,True,True)):
            deck=inversion_input(work,kind,split,full)
            tasks=4 if split else 2
            logs.append(launch(mopt,work,deck,tasks))
            launch(mopt,work,deck.replace('task=minvert','task=computegrad').replace('nlcg maxit=1 tolerance=1e100','lbfgs maxit=0'),tasks,label='gradient')
        serial_observations=observations(serial)
        for event,data in enumerate(serial_observations):
            expected=original.copy();expected[:,1:]*=1.25+event
            close_traces(expected,data,f'{name} event {event}')
        errors={}
        for work,log in zip((parallel,absolute),logs[1:]):
            for event,data in enumerate(observations(work)):
                close_traces(serial_observations[event],data,f'{work.name} observations {event}')
            errors[work.name]=close_gradients(gradient(serial),gradient(work),f'{name} full gradient')
            reference=float(re.findall(r'Initial misfit=\s*(\S+)',logs[0])[0])
            values=[float(x) for x in re.findall(r'Initial misfit=\s*(\S+)',log)]
            if not values or reference<=0:raise ValueError('Missing/nonpositive inversion misfit')
            for value in values:close(reference,value,f'{name} misfit',2e-5)
        record(name,{'gradient_max_errors':errors,'misfit':reference})

    if selected('sac-endian'):
        outputs=[]
        for endian,label in (('<','little'),('>','big')):
            write_observed(endian)
            work=root/f'endian-{label}'
            log=launch(mopt,work,inversion_input(work,'sac',coords=False),2)
            deck=inversion_input(work,'sac',coords=False)
            launch(mopt,work,deck.replace('task=minvert','task=computegrad').replace('nlcg maxit=1 tolerance=1e100','lbfgs maxit=0'),2,label='gradient')
            values=observations(work)
            outputs.append((values,gradient(work),float(re.findall(r'Initial misfit=\s*(\S+)',log)[0])))
        for event in (0,1):close_traces(outputs[0][0][event],outputs[1][0][event],f'endian event {event}')
        error=close_gradients(outputs[0][1],outputs[1][1],'endian full gradient')
        close(outputs[0][2],outputs[1][2],'endian misfit',2e-5)
        # Explicit coordinates must remain usable with undefined header location.
        write_observed(invalid='undefined')
        work=root/'endian-explicit'
        launch(mopt,work,inversion_input(work,'sac',coords=True),2)
        close_traces(outputs[0][0][0],observations(work)[0],'explicit coordinate samples')
        write_observed()
        for event in (0,1):
            for suffix in ('x','y','z'):
                (root/f'observed{event}/s.{suffix}').write_bytes((root/f'observed{event}/station.{suffix}').read_bytes())
        work=root/'endian-short-name'
        launch(mopt,work,inversion_input(work,'sac',coords=False,basename='s'),2)
        for event,data in enumerate(observations(work)):
            close_traces(outputs[0][0][event],data,f'short SAC name event {event}')
        record('sac-endian',{'gradient_max_error':error})

    for invalid in ('undefined','nan','range','truncated'):
        name=f'sac-{invalid}'
        if not selected(name):continue
        write_observed(invalid=invalid)
        work=root/name
        failure='invalid station coordinates' if invalid!='truncated' else 'invalid SAC observation'
        launch(mopt,work,inversion_input(work,'sac',coords=False),2,failure)
        record(name,{'rejected':True})

    for kind in ('usgs','sac'):
        name=f'events-{kind}-missing'
        if not selected(name):continue
        write_observed()
        missing=root/'observed1'/('station.txt' if kind=='usgs' else 'station.x')
        missing.rename(missing.with_suffix(missing.suffix+'.missing'))
        work=root/name
        launch(mopt,work,inversion_input(work,kind,True),4,
               'Could not open USGS receiver' if kind=='usgs' else 'invalid SAC observation')
        record(name,{'rejected':True})


if __name__=='__main__':main()
