#!/usr/bin/env python3
"""On-demand negative checks: failed solves and noise-only traces cannot pass."""
from pathlib import Path
import tempfile

import h5py
import numpy as np

from compare_waveforms import compare
from profile_hmr3 import compare as compare_text
from validation import check_solver_log


def rejected(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError('Invalid scientific evidence was accepted')


def main():
    with tempfile.TemporaryDirectory() as name:
        root = Path(name)
        log = root/'run.log'
        log.write_text('program sw4 finished!\n')
        check_solver_log(log)
        for failure in ('WARNING: no convergence in interface iteration',
                        'HDF5-DIAG: Error detected in HDF5', 'Error norm = nan',
                        'Error norm = +inf'):
            log.write_text(failure+'\nprogram sw4 finished!\n')
            rejected(lambda: check_solver_log(log))
        wave = root/'receiver.h5'
        for amplitude in (0., 1e-14, 1e-4):
            with h5py.File(wave,'w') as stream:
                stream['DELTA']=[.01]
                for station in ('station','edge'):
                    group=stream.create_group(station)
                    group['NPTS']=[40]; group['ISNSEW']=[0]
                    group['STX,STY,STZ']=[1.,2.,3.]
                    group['ACTUALSTX,STY,STZ']=[1.,2.,3.]
                    for component in ('X','Y','Z'):
                        group[component]=np.full(40,amplitude if component=='X' else 0.)
            action=lambda: compare(wave,wave,2e-5,1e-12,arrival_after=.15)
            text=root/'receiver.txt'
            np.savetxt(text,np.column_stack((np.arange(40)*.01,np.full(40,amplitude),np.zeros((40,2)))))
            if amplitude<1e-12:
                rejected(action); rejected(lambda: compare_text(text,text,2e-5,1e-12))
            else:
                action(); compare_text(text,text,2e-5,1e-12)
        # Good amplitudes with an insufficient trace length still fail.
        rejected(lambda: compare(wave,wave,2e-5,1e-12,arrival_after=.5))
    print('PASS: acceptance rejects failed solves, noise-only traces and missing arrival windows')


if __name__=='__main__': main()
