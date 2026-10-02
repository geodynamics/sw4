#!/usr/bin/env python3
"""Optional MPI regression for repeated ZFP-compressed SSI writes.

Exit status 77 means the test was skipped because SW4, HDF5, H5Z-ZFP, or
h5py is unavailable. Run this test inside an allocation on systems where MPI
launches are not permitted on login nodes.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


SKIP = 77
ZFP_FILTER_ID = 32013


def skip(message):
    print("SKIP:", message)
    return SKIP


def linked_optional_libraries(sw4):
    ldd = shutil.which("ldd")
    if ldd is None:
        return False, "ldd is unavailable; cannot verify optional libraries"
    result = subprocess.run(
        [ldd, str(sw4)], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )
    output = result.stdout.lower()
    if result.returncode != 0:
        return False, "ldd could not inspect the SW4 executable"
    if "not found" in output:
        return False, "the SW4 executable has unresolved shared libraries"
    if "hdf5" not in output:
        return False, "the SW4 executable is not linked with HDF5"
    if "h5zzfp" not in output:
        # CPU Make/CMake builds can link H5Z-ZFP statically.
        nm = shutil.which("nm")
        if nm is None:
            return False, "nm is unavailable; cannot verify static H5Z-ZFP"
        symbols = subprocess.run(
            [nm, "-g", str(sw4)], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL
        )
        if symbols.returncode != 0 or "H5Z_zfp_initialize" not in symbols.stdout:
            return False, "the SW4 executable is not linked with H5Z-ZFP"
    return True, ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sw4", required=True, type=Path,
                        help="path to an SW4 executable built with HDF5 and ZFP")
    parser.add_argument("--launcher", default="mpirun",
                        help="MPI launcher (default: mpirun)")
    parser.add_argument("--tasks", default=4, type=int,
                        help="number of MPI tasks (default: 4)")
    args = parser.parse_args()

    if not args.sw4.is_file() or not os.access(args.sw4, os.X_OK):
        return skip("SW4 executable is unavailable")
    if shutil.which(args.launcher) is None:
        return skip("MPI launcher '{}' is unavailable".format(args.launcher))
    available, reason = linked_optional_libraries(args.sw4)
    if not available:
        return skip(reason)

    try:
        import h5py
    except ImportError:
        return skip("h5py is unavailable")

    input_path = (Path(__file__).resolve().parent / "reference" / "hdf5" /
                  "ssi-zfp-dataset-lifetime.in")
    if not input_path.is_file():
        print("ERROR: missing regression input", input_path)
        return 1

    with tempfile.TemporaryDirectory(prefix="sw4-ssi-zfp-") as tmp:
        run_dir = Path(tmp)
        command = [args.launcher, "-n", str(args.tasks),
                   str(args.sw4.resolve()), str(input_path)]
        result = subprocess.run(
            command, cwd=str(run_dir), text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        )
        print(result.stdout)
        if result.returncode != 0:
            print("ERROR: SW4 returned", result.returncode)
            return 1
        if "HDF5-DIAG" in result.stdout:
            print("ERROR: HDF5 reported a diagnostic")
            return 1

        output = run_dir / "ssi-zfp.ssi"
        if not output.is_file():
            print("ERROR: missing SSI output", output)
            return 1

        with h5py.File(output, "r") as h5:
            if int(h5["lastsw4timestep"][0]) != 6:
                print("ERROR: incomplete SSI timestep progress")
                return 1
            if int(h5["lastoutputindex"][0]) != 5:
                print("ERROR: incomplete SSI output-index progress")
                return 1
            for component in range(3):
                dataset = h5["vel_{} ijk layout".format(component)]
                filters = [
                    dataset.id.get_create_plist().get_filter(index)[0]
                    for index in range(
                        dataset.id.get_create_plist().get_nfilters()
                    )
                ]
                if ZFP_FILTER_ID not in filters:
                    print("ERROR: velocity component", component,
                          "does not use the H5Z-ZFP filter")
                    return 1

    print("PASS: repeated ZFP-compressed SSI writes completed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
