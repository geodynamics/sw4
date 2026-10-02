#!/usr/bin/env python3
"""Check that a damaged receiver HDF5 file stops an SW4 restart.

Run inside a one-node CPU allocation with an HDF5-enabled SW4 executable.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import h5py


INPUT = """grid h=200 x=30000 y=30000 z=17000 extrapolate=1
time steps=8
fileio path=.
supergrid gp=30
block vp=6000 vs=3464 rho=2700
block vp=4000 vs=2000 rho=2600 z2=1000
source x=15000 y=15000 z=2000 mxy=1e18 t0=0.36 freq=16.6667 type=Gaussian
rec x=15600 y=15800 z=0 file=sta01 usgsformat=0 sacformat=0 hdf5format=1 hdf5file=receiver.h5 writeEvery=2
checkpoint cycleInterval=4 restartpath=. file=restart hdf5=yes{restartfile}
"""


def run(sw4, input_file, run_dir):
    result = subprocess.run(
        ["srun", "-N", "1", "-n", "1", "-c", "1", str(sw4), str(input_file)],
        cwd=str(run_dir), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True,
    )
    return result.returncode, result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sw4", required=True, type=Path)
    args = parser.parse_args()
    sw4 = args.sw4.resolve(strict=True)

    with tempfile.TemporaryDirectory(prefix="sw4-receiver-restart-") as name:
        run_dir = Path(name)
        base_input = run_dir / "base.in"
        restart_input = run_dir / "restart.in"
        base_input.write_text(INPUT.format(restartfile=""))
        restart_input.write_text(INPUT.format(
            restartfile=" restartfile=restart.cycle=4.sw4checkpoint"
        ))

        code, output = run(sw4, base_input, run_dir)
        receiver = run_dir / "receiver.h5"
        checkpoint = run_dir / "restart.cycle=4.sw4checkpoint"
        if code != 0 or not receiver.is_file() or not checkpoint.is_file():
            print("FAIL: base case did not create receiver and checkpoint")
            print(output[-8000:])
            return 1

        original = run_dir / "receiver-original.h5"
        shutil.copy2(receiver, original)
        with h5py.File(receiver, "r+") as stream:
            stream["sta01/X"][5] = 12345678.0
        code, output = run(sw4, restart_input, run_dir)
        if code != 0 or "Time step" not in output:
            print("FAIL: valid receiver restart did not advance")
            print(output[-8000:])
            return 1
        with h5py.File(receiver, "r") as stream:
            if float(stream["sta01/X"][5]) == 12345678.0:
                print("FAIL: restart did not replace samples beyond checkpoint")
                return 1

        shutil.copy2(original, receiver)
        os.truncate(receiver, 128)
        code, output = run(sw4, restart_input, run_dir)
        if code == 0 or "Could not restore receiver history" not in output:
            print("FAIL: truncated receiver restart did not fail closed")
            print(output[-8000:])
            return 1
        if "Time step" in output:
            print("FAIL: truncated receiver restart entered time stepping")
            print(output[-8000:])
            return 1

        shutil.copy2(original, receiver)
        with h5py.File(receiver, "r+") as stream:
            stream["sta01/NPTS"][...] = 2
        with h5py.File(receiver, "r") as stream:
            short_npts = int(stream["sta01/NPTS"][0])
            downsample = int(stream["DOWNSAMPLE"][0])
        if short_npts != 2 or (short_npts - 1) * downsample + 1 >= 4:
            print("FAIL: test fixture is not shorter than checkpoint:",
                  short_npts, downsample)
            return 1
        code, output = run(sw4, restart_input, run_dir)
        if code == 0 or "ends before checkpoint cycle" not in output:
            print("FAIL: short receiver history did not stop the restart")
            print("NPTS:", short_npts, "DOWNSAMPLE:", downsample)
            print(output[-8000:])
            return 1
        if "Time step" in output:
            print("FAIL: short receiver history entered time stepping")
            print(output[-8000:])
            return 1

        shutil.copy2(original, receiver)
        with h5py.File(receiver, "r+") as stream:
            stream["DOWNSAMPLE"][...] = 2
        code, output = run(sw4, restart_input, run_dir)
        if code == 0 or "receiver downsample mismatch" not in output:
            print("FAIL: mismatched receiver downSample was accepted")
            print(output[-8000:])
            return 1
        if "Time step" in output:
            print("FAIL: mismatched receiver downSample entered time stepping")
            print(output[-8000:])
            return 1

    print("PASS: valid receiver restart advanced; truncated, short, and mismatched histories stopped before time stepping")
    return 0


if __name__ == "__main__":
    sys.exit(main())
