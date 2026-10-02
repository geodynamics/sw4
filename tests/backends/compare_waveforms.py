#!/usr/bin/env python3
"""Compare small native/GPU simulations on demand in an existing Slurm allocation.

Checks inline and Sfile materials, Cartesian station metadata over topography,
attenuation, and text/HDF5 SRFs. Inputs, output files and launch logs are kept.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np

GRID = 'grid h=200 x=4000 y=4000 z=4000 lat=37 lon=-122 az=0\n'
BLOCK = 'block vp=6000 vs=3464 rho=2700\nblock vp=4000 vs=2000 rho=2600 z2=800\n'
SOURCE = 'source x=1600 y=1600 z=600 mxy=1e15 t0=0.05 freq=20 type=Gaussian\n'
COMMON = 'time steps=16\nsupergrid gp=4\nfileio path=output\n'
RECEIVERS = 'rechdf5 infile=stations.h5 outfile=receivers.h5 writeEvery=3\n'


def stations(path):
    with h5py.File(path, 'w') as f:
        for name, xyz in [('station', (1800, 1800, 200)), ('edge', (3200, 2800, 200))]:
            g = f.create_group(name)
            g['STX,STY,STZ'] = np.asarray(xyz, dtype=np.float64)
            g['USEZVALUE'] = np.asarray([1], dtype=np.int32)
            g['ISNSEW'] = np.asarray([0], dtype=np.int32)


def srf(path):
    # A smooth pulse with nonzero integral plus a zero-slip point exercises
    # normalization and source omission in both text and HDF5 readers.
    pulse = np.sin(np.linspace(0, np.pi, 32)) ** 2 * 100
    lines = ['2.0', 'PLANE 1', '-121.984375 37.015625 2 1 0.4 0.2',
             '0 90 0.6 0 0', 'POINTS 2']
    for zero in (False, True):
        values = np.zeros_like(pulse) if zero else pulse
        lines += ['-121.984375 37.015625 0.6 0 90 4e8 0.0 0.01 3.464e5 2.7',
                  f'0 {0 if zero else 15.5} 32 0 0 0 0']
        lines += [' '.join(f'{x:.9g}' for x in values[i:i+8]) for i in range(0, 32, 8)]
    path.write_text('\n'.join(lines) + '\n')


def run(exe, directory, text, tasks, gpu):
    directory.mkdir(parents=True, exist_ok=True)
    stations(directory / 'stations.h5')
    input_file = directory / 'case.in'
    input_file.write_text(text)
    command = ['srun', '--exclusive', '--exact', '-N', '1', '-n', str(tasks), '-c', '1']
    if gpu:
        command += ['--gpus-per-task=1', '--gpu-bind=single:1']
    command += [str(exe), str(input_file)]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', MPICH_GPU_SUPPORT_ENABLED='0')
    with (directory / 'run.log').open('w') as log:
        result = subprocess.run(command, cwd=directory, env=env, stdout=log,
                                stderr=subprocess.STDOUT, timeout=300)
    if result.returncode:
        raise RuntimeError(f'SW4 failed: {directory / "run.log"}')
    return directory / 'output' / 'receivers.h5'


def compare(a, b, rtol, atol):
    errors = {}
    with h5py.File(a) as fa, h5py.File(b) as fb:
        np.testing.assert_allclose(fa['DELTA'][:], fb['DELTA'][:], rtol=1e-12, atol=1e-12)
        for station in ('station', 'edge'):
            ga, gb = fa[station], fb[station]
            for name in ('NPTS', 'ISNSEW', 'STX,STY,STZ', 'ACTUALSTX,STY,STZ'):
                np.testing.assert_allclose(ga[name][:], gb[name][:], rtol=1e-12, atol=1e-12)
            signal = 0.0
            for component in ('X', 'Y', 'Z'):
                x, y = ga[component][:], gb[component][:]
                assert x.shape == y.shape and np.all(np.isfinite(x)) and np.all(np.isfinite(y))
                scale = float(np.max(np.abs(x)))
                delta = float(np.max(np.abs(x-y)))
                assert delta <= atol + rtol * scale, (station, component, delta, scale)
                signal = max(signal, scale)
                errors[f'{station}/{component}'] = {'max_difference': delta, 'reference_peak': scale}
            assert signal > 0, f'No signal at {station}'
    return errors


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu', type=Path, required=True)
    p.add_argument('--gpu', type=Path, required=True)
    p.add_argument('--work-dir', type=Path, required=True)
    p.add_argument('--case', action='append', choices=('inline', 'topography', 'attenuation', 'mesh-refinement', 'srf-text', 'srf-hdf5', 'sfile'), help='Run only selected cases; default: all')
    p.add_argument('--tasks', type=int, default=2)
    p.add_argument('--rtol', type=float, default=2e-5)
    p.add_argument('--atol', type=float, default=1e-12)
    args = p.parse_args()
    if args.tasks < 1:
        p.error('--tasks must be positive')
    cpu, gpu = args.cpu.resolve(strict=True), args.gpu.resolve(strict=True)
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[2]
    rupture = work / 'rupture.srf'
    srf(rupture)
    conversion = subprocess.run([sys.executable, str(repo/'tools/srf2hdf5.py'),
                                 str(rupture), str(work/'rupture.h5')],
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (work/'convert-srf.log').write_text(conversion.stdout)
    if conversion.returncode or not (work/'rupture.h5').exists():
        raise RuntimeError('SRF conversion failed')
    cases = {
        'inline': GRID + COMMON + BLOCK + SOURCE + RECEIVERS,
        'topography': GRID + COMMON + 'topography input=gaussian zmax=2000 order=4 '
            'gaussianAmp=100 gaussianXc=2000 gaussianYc=2000 gaussianLx=1500 gaussianLy=1500\n'
            + BLOCK + SOURCE + RECEIVERS,
        'mesh-refinement': GRID + COMMON + 'refinement zmax=2000\n' + BLOCK + SOURCE + RECEIVERS,
        'attenuation': GRID + COMMON + 'attenuation nmech=1\n'
            + BLOCK.replace('rho=2700', 'rho=2700 qp=200 qs=100').replace('rho=2600', 'rho=2600 qp=200 qs=100')
            + SOURCE + RECEIVERS,
        'srf-text': GRID + COMMON + BLOCK + f'rupture file={rupture}\n' + RECEIVERS,
        'srf-hdf5': GRID + COMMON + BLOCK + f'rupturehdf5 file={work / "rupture.h5"}\n' + RECEIVERS,
    }
    model_dir = work / 'generate-sfile'
    run(cpu, model_dir, cases['mesh-refinement'] + 'sfileoutput file=model sampleFactorH=2 sampleFactorV=1\n', args.tasks, False)
    model = model_dir / 'output/model.sfile'
    assert model.exists(), model
    cases['sfile'] = GRID + COMMON + f'sfile filename={model.name} directory={model.parent}\n' + SOURCE + RECEIVERS
    report = {}
    for name, text in cases.items():
        if args.case and name not in args.case:
            continue
        a = run(cpu, work/name/'cpu', text, args.tasks, False)
        b = run(gpu, work/name/'gpu', text, args.tasks, True)
        report[name] = compare(a, b, args.rtol, args.atol)
        (work/'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
        print(f'PASS: {name}', flush=True)
    # Readers of the same rupture representation should agree on each backend.
    for backend in ('cpu', 'gpu') if not args.case or {'srf-text', 'srf-hdf5'} <= set(args.case) else ():
        compare(work/'srf-text'/backend/'output/receivers.h5',
                work/'srf-hdf5'/backend/'output/receivers.h5', args.rtol, args.atol)
        print(f'PASS: text/HDF5 SRF consistency on {backend}', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
