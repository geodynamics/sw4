#!/usr/bin/env python3
"""Run production halo sentinels explicitly in a Slurm allocation (no pytest/CTest)."""
import argparse
import os
from pathlib import Path
import subprocess


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--backend',choices=('OPENMP','CUDA','HIP'),default='OPENMP')
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--nodes',type=int,default=1)
    parser.add_argument('--tasks',type=int,default=4)
    args=parser.parse_args()
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_') or k=='SLURM_JOB_ID'}
    env.update(OMP_NUM_THREADS='1')
    probe=str(args.probe.resolve(strict=True))
    for name,x,y in (('square-uneven',4600,5100),('x-strip',14000,2000),('y-strip',2000,14000)):
        work=root/name;work.mkdir(exist_ok=True)
        (work/'case.in').write_text(f'grid h=100 x={x} y={y} z=4000\n'
            'time steps=2\nfileio path=output\nsupergrid gp=8\n'
            'refinement zmax=2000\nblock vp=6000 vs=3464 rho=2700\n')
        command=['srun','--exclusive','--exact','-N',str(args.nodes),'-n',str(args.tasks),'-c','1']
        command+=['--gres=none'] if args.backend=='OPENMP' else ['--gpus-per-task=1','--gpu-bind=single:1']
        with (work/'run.log').open('w') as log:
            result=subprocess.run(command+[probe,'case.in'],cwd=work,env=env,
                stdout=log,stderr=subprocess.STDOUT,timeout=300)
        output=(work/'run.log').read_text()
        if result.returncode or 'PASS: 1/3/4/21-component halo sentinels' not in output:
            raise RuntimeError(f'Halo oracle failed: {work}/run.log\n{output[-4000:]}')
        expected=(args.tasks,1) if name=='x-strip' else (1,args.tasks) if name=='y-strip' else (2,2) if args.tasks==4 else None
        if expected and f'dimensions={expected[0]}x{expected[1]}' not in output:
            raise RuntimeError(f'Unexpected MPI decomposition: {work}/run.log')
        print('PASS:',name,'ranks=',args.tasks,'nodes=',args.nodes,flush=True)


if __name__=='__main__':main()
