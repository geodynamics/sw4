#!/usr/bin/env python3
"""On-demand full-rate USGS plus HDF5 restart and inconsistent-history checks.

Run inside a Slurm allocation. Each successful restart is compared with retained
uninterrupted waveforms and PDE state; rejected inputs must leave both histories
untouched. No tests are registered in the existing pytest suites.
"""
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
from receiver_restart import MODES, capture, compare, state
from validation import check_solver_log


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sw4', type=Path, required=True)
    parser.add_argument('--backend', choices=('OPENMP', 'CUDA', 'HIP'), default='OPENMP')
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--tasks', type=int, default=2)
    parser.add_argument('--case', action='append', help='Select baseline name, e.g. displacement-geo1-ds3-cycle24')
    parser.add_argument('--implicit-utc', action='store_true',
                        help='Use the default wall clock and test differing stored epochs')
    args = parser.parse_args()
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith('SLURM_') or k == 'SLURM_JOB_ID'}
    env.update(OMP_NUM_THREADS='1', MPICH_GPU_SUPPORT_ENABLED='0', DISABLE_LUSTRE_STRIPE='1')
    command = ['srun', '--exclusive', '--exact', '--cpu-bind=none', '-N', '1', '-n', str(args.tasks), '-c', '1']
    command += ['--gres=none'] if args.backend == 'OPENMP' else ['--gpus-per-task=1', '--gpu-bind=single:1']
    command += [str(args.sw4.resolve(strict=True)), 'case.in']
    report = {}

    def launch(work, text, label):
        (work / 'case.in').write_text(text)
        path = work / f'{label}.log'
        with path.open('w') as stream:
            result = subprocess.run(command, cwd=work, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=600)
        return result, path

    cases = [(mode, geo, ds, cycle) for mode in MODES
             for geo in ((0, 1) if mode in ('displacement', 'velocity') else (0,))
             for ds in (1, 3) for cycle in (24, 25)]
    for mode, geo, ds, cycle in cases:
        name = f'{mode}-geo{geo}-ds{ds}-cycle{cycle}'
        if args.case and name not in args.case:
            continue
        steps = 2 * cycle
        final_checkpoint = f'restart.cycle={steps}.sw4checkpoint'
        text = 'grid h=100 x=4000 y=4000 z=4000 lat=37 lon=-122 az=27\n'
        text += f'time steps={steps}' + ('' if args.implicit_utc else ' utcstart=10/02/2026:01:02:03.234567')
        text += '\nfileio path=.\nsupergrid gp=8\n'
        text += 'block vp=6000 vs=3464 rho=2700\n'
        text += 'source x=1637 y=1643 z=1000 mxy=1e15 mxz=2e15 myy=5e14 type=Gaussian freq=10 t0=0.3\n'
        text += f'rec x=1800 y=1800 z=1600 file={mode}-history sta=station variables={mode} nsew={geo} '
        text += f'sacformat=0 usgsformat=1 hdf5format=1 hdf5file=history.h5 downSample={ds} writeEvery=7\n'
        text += f'checkpoint cycleInterval={cycle} restartpath=. file=restart hdf5=yes\n'
        baseline = root / f'{name}-baseline'
        baseline.mkdir()
        result, path = launch(baseline, text, 'uninterrupted')
        if result.returncode:
            raise ValueError(f'Baseline failed: {path}')
        check_solver_log(path)
        reference_wave = capture(baseline)
        reference_state = state(baseline / final_checkpoint)
        np.savez_compressed(baseline / 'reference-waveforms.npz', **reference_wave)
        np.savez_compressed(baseline / 'reference-state.npz', **reference_state)
        reference_hashes = {file: digest(baseline / file) for file in
                            ('reference-waveforms.npz', 'reference-state.npz')}
        mutations = ['tail', 'sufficient']
        if name == 'displacement-geo1-ds3-cycle24':
            mutations += ['changed-tail', 'text-short', 'text-missing', 'text-origin', 'text-utc',
                          'text-nonuniform', 'disagree', 'hdf-short', 'hdf-origin',
                          'hdf-utc', 'hdf-downsample', 'hdf-component']
            if args.implicit_utc:
                mutations = [m for m in mutations if m not in ('text-utc', 'hdf-utc')] + ['epoch']
        for mutation in mutations:
            label = f'{name}-{mutation}'
            work = root / label
            shutil.copytree(baseline, work)
            usgs = work / f'{mode}-history.txt'
            hdf = work / 'history.h5'
            if mutation.startswith('text-') or mutation in ('sufficient', 'disagree', 'changed-tail', 'epoch'):
                lines = usgs.read_text().splitlines()
                header = '\n'.join(line for line in lines if line.startswith('#')) + '\n'
                values = np.loadtxt(usgs)
                if mutation == 'sufficient':
                    values = values[:cycle]
                elif mutation == 'text-short':
                    values = values[:cycle - 1]
                elif mutation == 'text-origin':
                    values[:, 0] += .1
                elif mutation == 'text-utc':
                    header = header.replace(':03.234567', ':04.234567')
                elif mutation == 'text-nonuniform':
                    values[2, 0] += .001
                elif mutation == 'disagree':
                    values[18, 1] += .01
                elif mutation == 'changed-tail':
                    values[cycle:, 1:] += .01
                elif mutation == 'epoch':
                    header = re.sub(r'/\d{4}:', '/2000:', header)
                if mutation == 'text-missing':
                    usgs.unlink()
                else:
                    with usgs.open('w') as stream:
                        stream.write(header)
                        np.savetxt(stream, values, fmt='%.17g')
            if mutation.startswith('hdf-') or mutation in ('sufficient', 'changed-tail', 'epoch'):
                with h5py.File(hdf, 'r+') as stream:
                    group = stream['station']
                    components = [key for key in group
                                  if key + 'CMPAZ' in group and key + 'CMPINC' in group]
                    if not components:
                        raise ValueError(f'No HDF5 waveform components: {hdf}')
                    if mutation in ('sufficient', 'hdf-short'):
                        npts = (cycle - 1) // ds + 1 if mutation == 'sufficient' else 2
                        group['NPTS'][...] = npts
                        for component in components:
                            group[component][npts:] = 0
                    elif mutation == 'changed-tail':
                        for component in components:
                            group[component][(cycle + ds - 1) // ds:] = .02
                    elif mutation == 'hdf-origin':
                        stream['STARTTIME'][...] = .1
                    elif mutation == 'hdf-utc':
                        stamp = stream.attrs['DATETIME']
                        if isinstance(stamp, bytes):
                            stamp = stamp.decode()
                        stream.attrs.modify('DATETIME', np.bytes_(stamp.replace(':03.', ':04.')))
                    elif mutation == 'hdf-downsample':
                        stream['DOWNSAMPLE'][...] = ds + 1
                    elif mutation == 'hdf-component':
                        del group[components[0]]
                    elif mutation == 'epoch':
                        stamp = stream.attrs['DATETIME']
                        if isinstance(stamp, bytes):
                            stamp = stamp.decode()
                        stream.attrs.modify('DATETIME', np.bytes_(re.sub(r'^\d{4}-', '2004-', stamp)))
            histories = [p for p in (usgs, hdf) if p.exists()]
            before = {p.name: digest(p) for p in histories}
            restart = text.rstrip() + f' restartfile=restart.cycle={cycle}.sw4checkpoint\n'
            result, path = launch(work, restart, 'restart')
            if mutation in ('tail', 'sufficient', 'changed-tail', 'epoch'):
                if result.returncode:
                    raise ValueError(f'Valid mixed restart failed: {path}')
                check_solver_log(path)
                record = {'wave_max_difference': compare(reference_wave, capture(work), label, True),
                          'state_max_difference': compare(reference_state, state(work / final_checkpoint), label, False),
                          'reference_directory': baseline.name, 'reference_sha256': reference_hashes,
                          'final_checkpoint': final_checkpoint}
            else:
                expected = {'text-short': 'ends before checkpoint cycle',
                            'text-missing': 'Could not open USGS receiver',
                            'text-origin': 'time origin', 'text-utc': 'time origin',
                            'text-nonuniform': 'Nonuniform USGS',
                            'disagree': 'histories disagree before checkpoint',
                            'hdf-short': 'ends before checkpoint cycle',
                            'hdf-origin': 'Could not restore receiver history',
                            'hdf-utc': 'Could not restore receiver history',
                            'hdf-downsample': 'downsample mismatch',
                            'hdf-component': 'Could not restore receiver history'}[mutation]
                log = path.read_text()
                if result.returncode == 0 or expected not in log or 'Begin time stepping' in log:
                    raise ValueError(f'Expected rejection before stepping: {path}, diagnostic={expected}')
                after = {p.name: digest(p) for p in (usgs, hdf) if p.exists()}
                if before != after:
                    raise ValueError(f'Failed restart rewrote receiver history: {work}')
                record = {'rejected_before_stepping': True, 'histories_preserved': True}
            report[label] = record
            (root / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
            print('PASS:', label, record, flush=True)
    if not report:
        raise ValueError('No selected test cases')


if __name__ == '__main__':
    main()
