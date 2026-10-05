#!/usr/bin/env python3
"""On-demand single-node OpenMP affinity and Cartesian refinement regression."""
import argparse
import json
import os
from pathlib import Path
import subprocess

import numpy as np
from validation import check_solver_log


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sw4',type=Path,required=True)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--tasks',type=int,default=8)
    args=parser.parse_args()
    if args.tasks<2:parser.error('Use at least two MPI ranks')
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_') or k=='SLURM_JOB_ID'}
    env.update(OMP_NUM_THREADS='1',OMP_PLACES='threads',OMP_PROC_BIND='close',
               MPICH_GPU_SUPPORT_ENABLED='0')
    launcher=['srun','--exclusive','--exact','--gres=none','--cpu-bind=threads',
              '-N','1','-n',str(args.tasks),'-c','1']
    with (root/'mpi.log').open('w') as log:
        result=subprocess.run(launcher+[str(args.probe.resolve(strict=True))],env=env,
                              stdout=log,stderr=subprocess.STDOUT,timeout=120)
    if result.returncode or 'PASS: distinct logical CPUs and 10000' not in (root/'mpi.log').read_text():
        raise RuntimeError(f'Affinity/communication probe failed: {root}/mpi.log')
    fixture=Path(__file__).resolve().parents[2]/'pytest/reference/energy/energy-mr-4th-1.in'
    (root/'case.in').write_bytes(fixture.read_bytes())
    with (root/'run.log').open('w') as log:
        result=subprocess.run(launcher+[str(args.sw4.resolve(strict=True)),'case.in'],
                              cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    check_solver_log(root/'run.log')
    if result.returncode or 'program sw4 finished!' not in (root/'run.log').read_text():
        raise RuntimeError(f'Cartesian refinement failed: {root}/run.log')
    energy=np.loadtxt(root/'energy-mr-4th-1/energy.log')
    if len(energy)!=535 or not np.all(np.isfinite(energy)) or energy[2]<=0:
        raise ValueError('Incomplete or invalid 535-step energy history')
    drift=float(np.max(np.abs(energy[2:]-energy[2]))/energy[2])
    if drift>1e-10:raise ValueError(f'Energy conservation error: {drift}')
    (root/'result.json').write_text(json.dumps({'ranks':args.tasks,'threads':1,
        'steps':535,'relative_energy_drift':drift},indent=2)+'\n')
    print('PASS: single-thread MPI affinity, halo/reductions, and 535-step conserved energy',flush=True)


if __name__=='__main__':main()
