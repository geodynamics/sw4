#!/usr/bin/env python3
"""On-demand checkpoint tests for shortened and misaligned receiver histories."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import h5py
import numpy as np
from receiver_restart import capture, compare, state
from validation import check_solver_log


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sw4',type=Path,required=True)
    parser.add_argument('--backend',choices=('OPENMP','CUDA','HIP'),default='OPENMP')
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--tasks',type=int,default=2)
    parser.add_argument('--case',action='append')
    args=parser.parse_args()
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    sw4=args.sw4.resolve(strict=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_') or k=='SLURM_JOB_ID'}
    env.update(OMP_NUM_THREADS='1',MPICH_GPU_SUPPORT_ENABLED='0',DISABLE_LUSTRE_STRIPE='1')
    command=['srun','--exclusive','--exact','--cpu-bind=none','-N','1','-n',str(args.tasks),'-c','1']
    command+=['--gres=none'] if args.backend=='OPENMP' else ['--gpus-per-task=1','--gpu-bind=single:1']
    command+=[str(sw4),'case.in']
    common='grid h=100 x=4000 y=4000 z=4000 lat=37 lon=-122 az=27\n'
    common+='time steps=48 utcstart=10/02/2026:01:02:03.234567\nfileio path=.\nsupergrid gp=8\n'
    common+='block vp=6000 vs=3464 rho=2700\nsource x=1637 y=1643 z=1000 mxy=1e15 mxz=2e15 myy=5e14 type=Gaussian freq=10 t0=0.3\n'
    checkpoint='checkpoint cycleInterval=24 restartpath=. file=restart hdf5=yes'
    report={}

    def launch(work,text,label):
        (work/'case.in').write_text(text)
        path=work/f'{label}.log'
        with path.open('w') as stream:
            result=subprocess.run(command,cwd=work,env=env,stdout=stream,
                                  stderr=subprocess.STDOUT,timeout=600)
        return result,path

    def record(name,result):
        report[name]=result
        (root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        print('PASS:',name,result,flush=True)

    for kind in ('text','sac','hdf1','hdf3'):
        mutations=('short','sufficient','tail','origin','utc')+ (('nonuniform',) if kind=='text' else ())
        if args.case and not any(f'{kind}-{mutation}' in args.case for mutation in mutations):continue
        options={'text':'sacformat=0 usgsformat=1','sac':'sacformat=1 usgsformat=0',
                 'hdf1':'sacformat=0 usgsformat=0 hdf5format=1 hdf5file=history.h5 downSample=1',
                 'hdf3':'sacformat=0 usgsformat=0 hdf5format=1 hdf5file=history.h5 downSample=3'}[kind]
        text=common+f'rec x=1800 y=1800 z=1600 file=displacement-history sta=station variables=displacement nsew=1 {options}\n'+checkpoint+'\n'
        baseline=root/f'{kind}-baseline';baseline.mkdir(exist_ok=True)
        result,path=launch(baseline,text,'uninterrupted')
        if result.returncode:raise ValueError(f'Baseline failed: {path}')
        check_solver_log(path)
        baseline_wave=capture(baseline);baseline_pde=state(baseline/'restart.cycle=48.sw4checkpoint')
        # Persist independent reference evidence before any restart overwrites it.
        np.savez_compressed(baseline/'reference-waveforms.npz',**baseline_wave)
        np.savez_compressed(baseline/'reference-state.npz',**baseline_pde)
        for mutation in mutations:
            name=f'{kind}-{mutation}'
            if args.case and name not in args.case:continue
            work=root/name;shutil.copytree(baseline,work)
            histories=[work/'displacement-history.txt'] if kind=='text' else (
                [work/f'displacement-history.{suffix}' for suffix in ('e','n','u')] if kind=='sac' else [work/'history.h5'])
            npts=2 if mutation=='short' else (24 if kind!='hdf3' else 8)
            for history in histories:
                if kind=='text':
                    lines=history.read_text().splitlines()
                    header='\n'.join(line for line in lines if line.startswith('#'))+'\n'
                    values=np.loadtxt(history)
                    if mutation in ('short','sufficient'):values=values[:npts]
                    if mutation=='origin':values[:,0]+=.1
                    if mutation=='nonuniform':values[2,0]+=.001
                    if mutation=='utc':header=header.replace(':03.234567',':04.234567')
                    with history.open('w') as stream:
                        stream.write(header);np.savetxt(stream,values,fmt='%.17g')
                elif kind=='sac':
                    raw=history.read_bytes()
                    floats=np.frombuffer(raw,dtype='<f4',count=70).copy()
                    ints=np.frombuffer(raw,dtype='<i4',count=40,offset=280).copy()
                    samples=np.frombuffer(raw,dtype='<f4',offset=632).copy()
                    if mutation in ('short','sufficient'):
                        samples=samples[:npts];ints[9]=len(samples)
                        floats[6]=floats[5]+(len(samples)-1)*floats[0]
                    if mutation=='origin':floats[5:7]+=.1
                    if mutation=='utc':ints[4]+=1
                    history.write_bytes(floats.tobytes()+ints.tobytes()+raw[440:632]+samples.tobytes())
                else:
                    with h5py.File(history,'r+') as stream:
                        if mutation in ('short','sufficient'):
                            group=stream['station'];group['NPTS'][...]=npts
                            for component in ('EW','NS','UP'):
                                values=group[component][:npts];del group[component]
                                # Preserve the originally planned storage capacity.
                                group.create_dataset(component,data=np.pad(values,(0,49-len(values))))
                        if mutation=='origin':stream['STARTTIME'][...]=.1
                        if mutation=='utc':
                            stamp=stream.attrs['DATETIME']
                            if isinstance(stamp,bytes):stamp=stamp.decode()
                            stream.attrs.modify('DATETIME',np.bytes_(stamp.replace(':03.',':04.')))
            before={str(p.relative_to(work)):digest(p) for p in histories}
            result,path=launch(work,text.rstrip()+' restartfile=restart.cycle=24.sw4checkpoint\n','restart')
            if mutation in ('sufficient','tail'):
                if result.returncode:raise ValueError(f'Valid restart failed: {path}')
                check_solver_log(path)
                wave=compare(baseline_wave,capture(work),name,True)
                pde=compare(baseline_pde,state(work/'restart.cycle=48.sw4checkpoint'),name,False)
                record(name,{'wave_max_difference':wave,'state_max_difference':pde,
                             'reference_wave_sha256':digest(baseline/'reference-waveforms.npz'),
                             'reference_state_sha256':digest(baseline/'reference-state.npz')})
            else:
                log=path.read_text()
                messages={'short':'ends before checkpoint cycle','origin':'time origin',
                          'utc':'time origin','nonuniform':'Nonuniform USGS'}
                expected=messages[mutation]
                if kind.startswith('hdf') and mutation in ('origin','utc'):expected='Could not restore receiver history'
                if result.returncode==0 or expected not in log or 'Begin time stepping' in log:
                    raise ValueError(f'Expected rejection before stepping: {path}, diagnostic={expected}')
                after={str(p.relative_to(work)):digest(p) for p in histories}
                if before!=after:raise ValueError(f'Failed restart rewrote receiver history: {work}')
                record(name,{'rejected_before_stepping':True,'histories_preserved':True})


if __name__=='__main__':main()
