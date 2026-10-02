#!/usr/bin/env python3
"""Run existing HDF5 regressions on demand with CPU or GPU launch settings.

Run from a login node with SLURM_JOB_ID set to an existing allocation and
--work-dir on a filesystem shared with its compute nodes. This script does
not submit jobs or register tests with pytest/CTest.
"""
import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sw4', type=Path, required=True)
    parser.add_argument('--backend', choices=('OPENMP', 'CUDA', 'HIP'), default='OPENMP')
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--tasks', type=int, default=4)
    args = parser.parse_args()
    if args.tasks < 1:
        parser.error('--tasks must be positive')
    sw4 = args.sw4.resolve(strict=True)
    real_srun = shutil.which('srun')
    if not real_srun:
        parser.error('srun is required; run inside or against a Slurm allocation')
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[2]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', MPICH_GPU_SUPPORT_ENABLED='0')
    with tempfile.TemporaryDirectory(prefix='launch-', dir=work) as name:
        launch = Path(name)
        wrapper = launch / 'srun'
        extra = ['--exclusive', '--exact']
        if args.backend != 'OPENMP':
            extra += ['--gpus-per-task=1', '--gpu-bind=closest']
        wrapper.write_text('#!/bin/sh\nexec ' + shlex.join([real_srun, *extra]) + ' "$@"\n')
        wrapper.chmod(0o755)
        env['PATH'] = str(launch) + os.pathsep + env.get('PATH', '')
        env['TMPDIR'] = str(work)
        cases = [
            ('receiver_hdf5_metadata.py', ['--tasks', str(args.tasks)]),
            ('receiver_hdf5_restart_failure.py', []),
            ('restart_hdf5_output_offsets.py', ['--tasks', str(args.tasks)]),
            ('ssi_zfp_dataset_lifetime.py', ['--launcher', 'srun', '--tasks', str(args.tasks)]),
        ]
        for script, options in cases:
            log = work / (Path(script).stem + '.log')
            command = [sys.executable, str(repo / 'pytest' / script), '--sw4', str(sw4), *options]
            with log.open('w') as stream:
                result = subprocess.run(command, cwd=work, env=env, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=1200)
            output = log.read_text()
            if result.returncode or 'PASS:' not in output or 'SKIP:' in output:
                print(output, file=sys.stderr)
                raise RuntimeError(f'{script} failed or skipped; see {log}')
            print(f'PASS: {script} ({log})', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
