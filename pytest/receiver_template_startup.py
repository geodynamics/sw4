#!/usr/bin/env python3
"""Compare native SW4 receiver outputs and retain startup timing evidence.

Run with EQSIM_PYTHON inside an authorized allocation. Both binaries must have
HDF5 enabled. This is a bounded regression, not a representative scale benchmark.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time
import threading

import h5py
import numpy as np


INPUT = """grid h=200 x=4000 y=4000 z=3000 extrapolate=1
time steps=6
fileio path=.
supergrid gp=4
block vp=6000 vs=3464 rho=2700
source x=1600 y=1600 z=200 mxy=1e18 t0=0.01 freq=16.6667 type=Gaussian
rec x=1800 y=1800 z=0 file=ascii usgsformat=1 sacformat=0
rechdf5 infile={template} outfile=receiver.h5 writeEvery=2 downSample=2
"""


def template(path, count, malformed=False):
    with h5py.File(path, "w") as stream:
        stream.create_dataset("header", data=[1])
        for index in range(count):
            group = stream.create_group("s{:06d}".format(index))
            coords = (800 + 200 * (index % 13),
                      800 + 200 * ((index // 13) % 13), 0)
            group["STX,STY,STZ"] = np.asarray(coords, dtype=np.float64)
            if index % 3:
                group["ISNSEW"] = np.asarray([index % 2], dtype=np.int32)
            if index % 7 == 0:
                # Both coordinate representations must retain Cartesian priority.
                group["STLA,STLO,STDP"] = np.asarray([31, -117, 0], dtype=np.float64)
            group["USEZVALUE"] = np.asarray([1], dtype=np.int32)
        if malformed:
            group = stream.create_group("z_bad")
            group["STX,STY,STZ"] = np.asarray([b"x", b"y", b"z"], dtype="S1")


def run(args, binary, fixture, destination, expect_failure=False, second=False):
    destination.mkdir()
    deck = INPUT.format(template=fixture)
    if second:
        deck += "rechdf5 infile={} outfile=second.h5 writeEvery=2 downSample=2\n".format(fixture)
    (destination / "case.in").write_text(deck)
    command = [args.launcher, "--overlap", "--exact", "--nodes=1", "--ntasks=" + str(args.tasks),
               "--cpus-per-task=" + str(args.cpus_per_task), "--cpu-bind=cores"]
    if args.gpu:
        command.extend(["--gpus-per-task=1", "--gpu-bind=closest"])
    command.extend([str(binary), str(destination / "case.in")])
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    start = time.monotonic()
    first_step = None
    with (destination / "sw4.log").open("w") as log:
        process = subprocess.Popen(command, cwd=destination, env=environment,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True)
        timed_out = threading.Event()
        def terminate():
            timed_out.set()
            process.terminate()
        timer = threading.Timer(120, terminate)
        timer.daemon = True
        timer.start()
        try:
            # SW4 output can be buffered; this observed timestamp includes that.
            for line in process.stdout:
                log.write(line)
                if first_step is None and re.search(r"Time step\s+(?:=\s*)?1\b", line):
                    first_step = time.monotonic() - start
            code = process.wait(timeout=120)
        except BaseException:
            process.kill()
            process.wait()
            raise
        finally:
            timer.cancel()
        if timed_out.is_set():
            raise TimeoutError("SW4 regression exceeded 120 seconds")
    output = (destination / "sw4.log").read_text()
    stepped = re.search(r"(?i)(?:\btime step\s*(?:=\s*)?\d+|\bbegin time stepping|\btime stepping finished)", output)
    if expect_failure:
        assert code != 0 and not stepped, output[-8000:]
    else:
        assert code == 0 and "HDF5-DIAG" not in output, output[-8000:]
        assert stepped, output[-8000:]
    return {"command": command, "exit_code": code,
            "elapsed": time.monotonic() - start,
            "first_step_observed": first_step,
            "receiver_setup": re.findall(r"Receiver HDF5 setup:.*", output)}


def compare(original, changed, count):
    peak = 0.0
    with h5py.File(original, "r") as baseline, h5py.File(changed, "r") as candidate:
        assert set(baseline) == set(candidate)
        stations = [name for name in baseline if name.startswith("s")]
        assert len(stations) == count
        for name in stations:
            left, right = baseline[name], candidate[name]
            assert set(left) == set(right), name
            for field in left:
                if field == "DATETIME":
                    continue  # Wall-clock publication metadata is not solver motion.
                a, b = left[field][...], right[field][...]
                if np.issubdtype(a.dtype, np.number):
                    np.testing.assert_allclose(a, b, rtol=2e-12, atol=1e-15,
                                               equal_nan=True, err_msg=name+"/"+field)
                else:
                    np.testing.assert_array_equal(a, b, err_msg=name+"/"+field)
            assert int(right["NPTS"][0]) == 4
            components = ("EW", "NS", "UP") if int(right["ISNSEW"][0]) else ("X", "Y", "Z")
            for component in components:
                values = right[component][...]
                assert values.shape == (4,) and np.all(np.isfinite(values))
                peak = max(peak, float(np.max(np.abs(values))))
    assert peak > 0, "Regression fixture produced no nonzero receiver motion"
    return peak


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-sw4", type=Path, required=True)
    parser.add_argument("--candidate-sw4", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--tasks", type=int, default=2)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--cpus-per-task", type=int, default=4)
    parser.add_argument("--launcher", default="srun")
    args = parser.parse_args()
    args.work_dir = args.work_dir.resolve()
    args.work_dir.mkdir(parents=True, exist_ok=False)
    args.baseline_sw4 = args.baseline_sw4.resolve(strict=True)
    args.candidate_sw4 = args.candidate_sw4.resolve(strict=True)
    results = {}
    try:
        for count in (8, 2053):
            fixture = args.work_dir / ("template-{}.h5".format(count))
            template(fixture, count)
            for label, binary in (("baseline", args.baseline_sw4),
                                  ("candidate", args.candidate_sw4)):
                directory = args.work_dir / (label + "-" + str(count))
                results[directory.name] = run(args, binary, fixture, directory,
                                              second=count == 8)
            for filename in (("receiver.h5", "second.h5") if count == 8 else ("receiver.h5",)):
                peak = compare(args.work_dir / ("baseline-" + str(count)) / filename,
                               args.work_dir / ("candidate-" + str(count)) / filename, count)
                print("PASS:", count, filename, "metadata and motion match; peak", peak)
        for name, count, malformed in (("missing", 0, False), ("bad", 0, True),
                                       ("bad-after-batch", 1024, True)):
            fixture = args.work_dir / (name + ".h5")
            if name != "missing":
                template(fixture, count, malformed)
            results[name] = run(args, args.candidate_sw4, fixture,
                                args.work_dir / name, expect_failure=True)
            print("PASS:", name, "fails before time stepping")
    finally:
        (args.work_dir / "results.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
