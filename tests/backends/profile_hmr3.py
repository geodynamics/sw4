#!/usr/bin/env python3
"""On-demand four-node hmr3 comparison in an existing Slurm allocation.

Only model/output paths are relocated. Each full simulation retains its input,
log, waveform and executable identity. Not registered with pytest or CTest.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time

import numpy as np


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def relocate(text, model):
    result = []
    counts = {'rfile': 0, 'topography': 0, 'fileio': 0}
    for line in text.splitlines():
        command = line.lstrip().split(maxsplit=1)[0] if line.strip() else ''
        if command == 'rfile':
            line = re.sub(r'filename=\S+', f'filename={model.name}', line)
            line = re.sub(r'directory=\S+', f'directory={model.parent}', line)
        elif command == 'topography':
            line = re.sub(r'\bfile=\S+', f'file={model}', line)
        elif command == 'fileio':
            line = re.sub(r'path=\S+', 'path=output', line)
        if command in counts:
            counts[command] += 1
        result.append(line)
    if any(v != 1 for v in counts.values()):
        raise ValueError(f'Expected the original hmr3 model/output commands: {counts}')
    return '\n'.join(result) + '\n'


def compare(a, b, rtol, atol):
    x, y = np.loadtxt(a), np.loadtxt(b)
    if x.shape != y.shape or x.ndim != 2 or x.shape[1] != 4:
        raise ValueError(f'Waveform dimensions differ: {x.shape}, {y.shape}')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError('Nonfinite station waveform')
    np.testing.assert_allclose(x[:, 0], y[:, 0], rtol=1e-10, atol=1e-10)
    peaks = np.max(np.abs(x[:, 1:]), axis=0)
    delta = np.max(np.abs(x[:, 1:] - y[:, 1:]), axis=0)
    if not np.any(peaks > 0) or np.any(delta > atol + rtol * peaks):
        raise ValueError(f'Waveform mismatch: delta={delta}, reference peaks={peaks}')
    return {'max_difference': delta.tolist(), 'reference_peak': peaks.tolist()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--input', type=Path, required=True, help='Exact performance/large/hmr3.in from raja')
    p.add_argument('--model', type=Path, required=True)
    p.add_argument('--work-dir', type=Path, required=True)
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--only', choices=('baseline', 'candidate'), help='Run one variant; summarize with both variants later')
    p.add_argument('--resume', action='store_true', help='Reuse completed runs with matching identities')
    p.add_argument('--max-slowdown', type=float, default=0.05)
    p.add_argument('--rtol', type=float, default=2e-5)
    p.add_argument('--atol', type=float, default=1e-12)
    args = p.parse_args()
    if args.repeats < 1 or args.max_slowdown < 0 or args.rtol < 0 or args.atol < 0:
        p.error('Invalid repeats or acceptance threshold')
    if not os.environ.get('SLURM_JOB_ID'):
        p.error('An existing four-node Slurm allocation is required')
    executables = {key: getattr(args, key).resolve(strict=True) for key in ('baseline', 'candidate')}
    model, source = args.model.resolve(strict=True), args.input.resolve(strict=True)
    if any(c.isspace() for c in str(model)):
        p.error('SW4 model paths must not contain whitespace')
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    text = relocate(source.read_text(), model)
    common = {'source_input_sha256': digest(source), 'relocated_input_sha256': hashlib.sha256(text.encode()).hexdigest(),
              'model': str(model), 'model_size': model.stat().st_size, 'model_mtime_ns': model.stat().st_mtime_ns,
              'nodes': 4, 'tasks': 16, 'tasks_per_node': 4, 'cpus_per_task': 4, 'gpus_per_task': 1,
              'omp_threads': 1, 'mpi_gpu_support': 0}
    records = {}
    for i in range(args.repeats):
        order = ('baseline', 'candidate') if i % 2 == 0 else ('candidate', 'baseline')
        for key in order:
            if args.only and key != args.only:
                continue
            exe = executables[key]
            identity = dict(common, executable=str(exe), executable_sha256=digest(exe))
            directory = work / f'{key}-{i}'
            record_file = directory / 'result.json'
            if args.resume and record_file.exists():
                record = json.loads(record_file.read_text())
                if record['identity'] != identity:
                    raise ValueError(f'Cannot resume a different simulation: {record_file}')
            else:
                directory.mkdir(exist_ok=False)
                input_file = directory / 'hmr3.in'
                input_file.write_text(text)
                command = ['srun', '--exclusive', '--exact', '-N', '4', '-n', '16', '--ntasks-per-node=4',
                           '-c', '4', '--gpus-per-task=1', '--gpu-bind=single:1', str(exe), str(input_file)]
                # A coordinator step can export one-node/one-CPU limits. Launch
                # against the allocation, without inheriting those step limits.
                env = {k: v for k, v in os.environ.items() if not k.startswith('SLURM_')}
                env.update(SLURM_JOB_ID=os.environ['SLURM_JOB_ID'], OMP_NUM_THREADS='1',
                           MPICH_GPU_SUPPORT_ENABLED='0')
                started = time.monotonic()
                with (directory / 'run.log').open('w') as log:
                    result = subprocess.run(command, cwd=directory, env=env, stdout=log,
                                            stderr=subprocess.STDOUT, timeout=2400)
                wall = time.monotonic() - started
                log_text = (directory / 'run.log').read_text()
                match = re.search(r'Execution time, solver phase\s+(?:(\d+) hours?\s+)?(?:(\d+) minutes?\s+)?([\d.]+) seconds', log_text)
                if result.returncode or not match or 'program sw4 finished!' not in log_text:
                    raise RuntimeError(f'Incomplete SW4 run: {directory / "run.log"}')
                waveform = directory / 'output/sta1.txt'
                compare(waveform, waveform, args.rtol, args.atol)
                record = {'identity': identity, 'job_id': os.environ['SLURM_JOB_ID'],
                          'node_list': os.environ.get('SLURM_JOB_NODELIST'), 'command': command,
                          'wall_seconds': wall, 'solver_seconds': 3600 * int(match.group(1) or 0) + 60 * int(match.group(2) or 0) + float(match.group(3)),
                          'waveform_sha256': digest(waveform)}
                record_file.write_text(json.dumps(record, indent=2) + '\n')
            records[f'{key}-{i}'] = record
            print(f'PASS {key}-{i}: solver {record["solver_seconds"]:.2f}s, wall {record["wall_seconds"]:.2f}s', flush=True)
    if args.only:
        return
    comparisons = [compare(work/f'baseline-{i}/output/sta1.txt', work/f'candidate-{i}/output/sta1.txt',
                           args.rtol, args.atol) for i in range(args.repeats)]
    medians = {key: statistics.median(records[f'{key}-{i}']['solver_seconds'] for i in range(args.repeats))
               for key in executables}
    slowdown = medians['candidate'] / medians['baseline'] - 1
    passed = slowdown <= args.max_slowdown
    report = {'runs': records, 'waveform_comparisons': comparisons, 'rtol': args.rtol, 'atol': args.atol,
              'median_solver_seconds': medians, 'slowdown': slowdown, 'max_slowdown': args.max_slowdown,
              'performance_passed': passed}
    (work/'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'Waveforms PASS; candidate slowdown {slowdown:.2%}; performance {"PASS" if passed else "FAIL"}', flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
