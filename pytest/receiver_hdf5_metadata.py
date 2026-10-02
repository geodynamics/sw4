#!/usr/bin/env python3
"""Verify CPU receiver metadata and separate HDF5 output files.

Run in an allocation with an HDF5-enabled CPU SW4 executable. The case includes
an ASCII receiver before two rechdf5 commands, repeated station names across
output files, Cartesian coordinates without ISNSEW, and mixed orientations.
"""

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


INPUT = """grid h=200 x=4000 y=4000 z=3000 extrapolate=1
time steps=4
fileio path=.
supergrid gp=4
block vp=6000 vs=3464 rho=2700
source x=1600 y=1600 z=200 mxy=1e18 t0=0.01 freq=16.6667 type=Gaussian
rec x=1800 y=1800 z=0 file=ascii usgsformat=1 sacformat=0
rechdf5 infile=stations.h5 outfile=first.h5 writeEvery=2
rechdf5 infile=stations.h5 outfile=second.h5 writeEvery=2
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sw4", required=True, type=Path)
    parser.add_argument("--tasks", type=int, default=2)
    parser.add_argument("--launcher", default="srun")
    args = parser.parse_args()
    sw4 = args.sw4.resolve(strict=True)

    with tempfile.TemporaryDirectory(prefix="sw4-receiver-metadata-") as name:
        root = Path(name)
        with h5py.File(root / "stations.h5", "w") as stream:
            for station, coords, orientation in (
                ("sta01", (1800, 1800, 0), None),
                ("sta02", (2000, 1800, 0), 0),
            ):
                group = stream.create_group(station)
                group["STX,STY,STZ"] = np.asarray(coords, dtype=np.float64)
                if orientation is not None:
                    group["ISNSEW"] = np.asarray([orientation], dtype=np.int32)
                group["USEZVALUE"] = np.asarray([1], dtype=np.int32)
                group["WINDOWS"] = np.asarray([0, 1, 2, 3], dtype=np.float64)
        input_file = root / "case.in"
        # Also exercise an input line beyond the former 256-byte buffer.
        input_file.write_text(INPUT.replace(
            "file=ascii usgsformat=1 sacformat=0",
            "file=ascii usgsformat=1 sacformat=0 #" + "x" * 280,
        ))
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        result = subprocess.run(
            [args.launcher, "-n", str(args.tasks), str(sw4), str(input_file)],
            cwd=root, env=env, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=120,
        )
        if result.returncode != 0 or "HDF5-DIAG" in result.stdout:
            raise AssertionError("SW4 failed:\n" + result.stdout)
        for filename in ("first.h5", "second.h5"):
            with h5py.File(root / filename, "r") as stream:
                for station, orientation, components in (
                    ("sta01", 1, ("EW", "NS", "UP")),
                    ("sta02", 0, ("X", "Y", "Z")),
                ):
                    group = stream[station]
                    assert int(group["ISNSEW"][0]) == orientation
                    assert int(group["NPTS"][0]) == 5
                    for component in components:
                        values = group[component][:]
                        assert values.shape == (5,)
                        assert np.all(np.isfinite(values))
                    assert np.max(np.abs(group[components[0]][:])) > 0
        # Unreadable coordinates must stop parsing before the solver starts.
        with h5py.File(root / "stations.h5", "r+") as stream:
            group = stream.create_group("sta03")
            group["STX,STY,STZ"] = np.asarray([b"x", b"y", b"z"], dtype="S1")
        failed = subprocess.run(
            [args.launcher, "-n", str(args.tasks), str(sw4), str(input_file)],
            cwd=root, env=env, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=120,
        )
        assert failed.returncode != 0, failed.stdout
        assert "Failed reading station metadata" in failed.stdout, failed.stdout
        assert "Time step" not in failed.stdout, failed.stdout
        print("PASS: long input, Cartesian station metadata, mixed orientations, "
              "two receiver output files, and malformed metadata rejection")


if __name__ == "__main__":
    main()
