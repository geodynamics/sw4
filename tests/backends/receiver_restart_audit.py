#!/usr/bin/env python3
"""Replay retained receiver/PDE references against final restarted outputs on demand.

Run in a compute allocation after receiver_restart[ _negative].py. This reads
artifacts without rerunning or rewriting simulations. Earlier restart repetitions
are recorded by the runner; this independent audit checks the retained final one.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from receiver_restart import capture, compare, state


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-dir', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit = []
    for root in args.evidence_dir:
        results = json.loads((root / 'comparison.json').read_text())
        cases = {}
        for name, result in results.items():
            if result.get('rejected_before_stepping'):
                continue
            folder = name.split('/')[0]
            if 'reference_sha256' in result:
                reference = root / folder
                hashes = result['reference_sha256']
            else:
                reference = root / (folder.split('-')[0] + '-baseline')
                hashes = {'reference-waveforms.npz': result['reference_wave_sha256'],
                          'reference-state.npz': result['reference_state_sha256']}
            cases[folder] = (reference, hashes)
        if not cases:
            raise ValueError(f'No successful restarted outputs: {root}')
        for folder, (reference, hashes) in cases.items():
            for filename, expected in hashes.items():
                if digest(reference / filename) != expected:
                    raise ValueError(f'Reference identity changed: {reference / filename}')
            with np.load(reference / 'reference-waveforms.npz', allow_pickle=False) as stream:
                waveform = {name: stream[name] for name in stream.files}
            with np.load(reference / 'reference-state.npz', allow_pickle=False) as stream:
                pde = {name: stream[name] for name in stream.files}
            work = root / folder
            restored = capture(work)
            wave_error = compare(waveform, restored, str(work), True)
            state_error = compare(pde, state(work / 'restart.cycle=48.sw4checkpoint'), str(work), False)
            files = {name.split('/')[0] for name in restored}
            files.add('restart.cycle=48.sw4checkpoint')
            restored_hashes = {name: digest(work / name) for name in sorted(files)}
            audit.append({'root': str(root), 'case': folder, 'reference': str(reference),
                          'reference_sha256': hashes, 'restored_sha256': restored_hashes,
                          'wave_max_difference': wave_error,
                          'state_max_difference': state_error})
            print('PASS:', work, 'retained references match final restored outputs', flush=True)
    args.output.write_text(json.dumps(audit, indent=2) + '\n')


if __name__ == '__main__':
    main()
