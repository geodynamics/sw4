#!/usr/bin/env python3
"""Check SW4 receiver startup and depth with native HDF5-enabled binaries.

Run in an authorized allocation. The startup mode compares matched baseline
and candidate builds. The depth mode compares template and inline receivers
over Gaussian topography. Startup timings are smoke-test measurements.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time
import signal

import h5py
import numpy as np


STARTUP_INPUT = """grid h=200 x=4000 y=4000 z=3000 extrapolate=1
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
            if index % 11 == 0:
                group["WINDOWS"] = np.asarray([0.1, 0.2, 0.3, 0.4])
        if malformed:
            group = stream.create_group("z_bad")
            group["STX,STY,STZ"] = np.asarray([b"x", b"y", b"z"], dtype="S1")


def launch(command, destination, environment, timeout=120, grace=5):
    """Retain output without blocking on a pipe held by a surviving child."""
    start = time.monotonic()
    deadline = start + timeout
    first_step = None
    path = destination / "sw4.log"
    with path.open("w") as log, path.open("r") as observed:
        process = subprocess.Popen(command, cwd=destination, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            while True:
                if first_step is None:
                    observed.seek(0)
                    if re.search(r"Time step\s+(?:=\s*)?1\b", observed.read()):
                        first_step = time.monotonic() - start
                try:
                    code = process.wait(timeout=max(0, min(0.1, deadline - time.monotonic())))
                    if first_step is None:
                        observed.seek(0)
                        if re.search(r"Time step\s+(?:=\s*)?1\b", observed.read()):
                            first_step = time.monotonic() - start
                    return code, first_step
                except subprocess.TimeoutExpired:
                    if time.monotonic() >= deadline:
                        raise TimeoutError("SW4 regression exceeded {} seconds".format(timeout))
        except BaseException:
            # The launcher and local descendants have a task-owned process group.
            for sig in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(process.pid, sig)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=grace)
                except subprocess.TimeoutExpired:
                    pass
            raise


def run_solver(args, binary, destination, deck, *, layout=()):
    destination.mkdir()
    deck_path = destination / "case.in"
    deck_path.write_text(deck)
    command = [args.launcher, "--overlap", "--exact", "--nodes=1",
               f"--ntasks={args.tasks}", f"--cpus-per-task={args.cpus_per_task}", *layout]
    if args.job_id:
        command.append(f"--jobid={args.job_id}")
    if args.gpu:
        command.extend(["--gpus-per-task=1", "--gpu-bind=closest"])
    command.extend([str(binary), str(deck_path)])
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    start = time.monotonic()
    code, first_step = launch(command, destination, environment)
    return {
        "command": command,
        "exit_code": code,
        "elapsed": time.monotonic() - start,
        "first_step_observed": first_step,
    }, (destination / "sw4.log").read_text()


def run_startup_case(args, binary, fixture, destination, expected_error=None, second=False):
    deck = STARTUP_INPUT.format(template=fixture)
    if second:
        deck += f"rechdf5 infile={fixture} outfile=second.h5 writeEvery=2 downSample=2\n"
    result, output = run_solver(args, binary, destination, deck, layout=("--cpu-bind=cores",))
    code = result["exit_code"]
    stepped = re.search(r"(?i)(?:\btime step\s*(?:=\s*)?\d+|\bbegin time stepping|\btime stepping finished)", output)
    if expected_error is not None:
        assert code != 0 and not stepped and expected_error in output, (
            f"Expected solver error {expected_error!r} before stepping:\n" + output[-8000:])
    else:
        assert code == 0 and "HDF5-DIAG" not in output, output[-8000:]
        assert stepped, output[-8000:]
    result["receiver_setup"] = re.findall(r"Receiver HDF5 setup:.*", output)
    return result


def compare_station(left, right, *, ignore_datetime=False):
    assert set(left) == set(right), left.name
    for field in left:
        if ignore_datetime and field == "DATETIME":
            continue  # Wall-clock publication metadata is not solver motion.
        a, b = left[field][...], right[field][...]
        if np.issubdtype(a.dtype, np.number):
            np.testing.assert_allclose(a, b, rtol=2e-12, atol=1e-15,
                                       equal_nan=True, err_msg=left.name + "/" + field)
        else:
            np.testing.assert_array_equal(a, b, err_msg=left.name + "/" + field)


def compare(original, changed, count):
    peak = 0.0
    with h5py.File(original, "r") as baseline, h5py.File(changed, "r") as candidate:
        assert set(baseline) == set(candidate)
        stations = [name for name in baseline if name.startswith("s")]
        assert len(stations) == count
        for name in stations:
            left, right = baseline[name], candidate[name]
            compare_station(left, right, ignore_datetime=True)
            assert int(right["NPTS"][0]) == 4
            components = ("EW", "NS", "UP") if int(right["ISNSEW"][0]) else ("X", "Y", "Z")
            for component in components:
                values = right[component][...]
                assert values.shape == (4,) and np.all(np.isfinite(values))
                peak = max(peak, float(np.max(np.abs(values))))
    assert peak > 0, "Regression fixture produced no nonzero receiver motion"
    return peak


DEPTH_INPUT = """grid h=200 x=4000 y=4000 z=3000 extrapolate=1
time steps=6
fileio path=.
supergrid gp=4
topography input=gaussian zmax=2000 order=4 gaussianAmp=100 gaussianXc=2000 gaussianYc=2000 gaussianLx=1500 gaussianLy=1500
block vp=6000 vs=3464 rho=2700
source x=1600 y=1600 z=200 mxy=1e18 t0=0.01 freq=16.6667 type=Gaussian
rec x=1800 y=1800 z=300 file=s_absolute usgsformat=0 sacformat=0 hdf5format=1 hdf5file=inline.h5 writeEvery=2 downSample=2 nsew=0
rec x=1800 y=1800 depth=300 file=s_relative usgsformat=0 sacformat=0 hdf5format=1 hdf5file=inline.h5 writeEvery=2 downSample=2 nsew=0
rechdf5 infile={template} outfile=receiver.h5 writeEvery=2 downSample=2
"""


def run_depth_case(args, binary, root, template, input_text=DEPTH_INPUT,
                   expect_topography_failure=False):
    result, output = run_solver(args, binary, root, input_text.format(template=template),
                                layout=("--ntasks-per-node=2",))
    code = result["exit_code"]
    if expect_topography_failure and code != 0:
        assert 'GaussianHill::inverse__grid_mapping: poor convergence' in output, output[-8000:]
        assert '-1.00000000e+38' in output and 'Time step ' not in output, output[-8000:]
        return None
    assert code==0 and 'HDF5-DIAG' not in output,output[-8000:]
    return root


def check_depth(root, unsupported=False):
    with h5py.File(root/'inline.h5','r') as inline,h5py.File(root/'receiver.h5','r') as template:
        compare_station(inline['s_relative'],template['s_relative'])
        assert any(not np.allclose(inline['s_absolute'][component][...],
                                   inline['s_relative'][component][...],
                                   rtol=2e-12,atol=1e-15) for component in 'XYZ'), \
            'Depth fixture does not distinguish absolute and relative station placement'
        if unsupported:
            differences=[]
            for component in 'XYZ':
                a=inline['s_absolute'][component][...]
                b=template['s_absolute'][component][...]
                differences.append(not np.allclose(a,b,rtol=2e-12,atol=1e-15))
            assert any(differences),'Baseline did not reproduce the missing absolute-depth behavior'
        else:
            compare_station(inline['s_absolute'],template['s_absolute'])
            peak=max(float(np.max(np.abs(template[name][component][...])))
                     for name in ('s_absolute','s_relative') for component in 'XYZ')
            assert peak>0


def elevated_absolute(args, binary):
    """Both representations must retain a receiver below terrain at negative z."""
    template = args.work_dir / 'elevated-stations.h5'
    with h5py.File(template, 'w') as stream:
        for name, z, use_z in [('s_absolute', -50, 1), ('s_relative', 50, 0)]:
            group = stream.create_group(name)
            group['STX,STY,STZ'] = np.asarray([2000, 2000, z], dtype=np.float64)
            group['ISNSEW'] = np.asarray([0], dtype=np.int32)
            group['USEZVALUE'] = np.asarray([use_z], dtype=np.int32)
    deck = DEPTH_INPUT.replace('x=1800 y=1800 z=300', 'x=2000 y=2000 z=-50')
    deck = deck.replace('x=1800 y=1800 depth=300', 'x=2000 y=2000 depth=50')
    root = run_depth_case(args, binary, args.work_dir/'elevated', template, input_text=deck)
    with h5py.File(root/'inline.h5', 'r') as inline, h5py.File(root/'receiver.h5', 'r') as receivers:
        for name in ('s_absolute', 's_relative'):
            assert name in inline and name in receivers, 'Valid elevated receiver omitted: ' + name
            compare_station(inline[name], receivers[name])
        for component in 'XYZ':
            absolute = receivers['s_absolute'][component][...]
            relative = receivers['s_relative'][component][...]
            assert np.all(np.isfinite(absolute))
            np.testing.assert_allclose(absolute, relative, rtol=2e-12, atol=1e-15)
        assert any(np.max(np.abs(receivers['s_absolute'][component][...])) > 0 for component in 'XYZ')
    print('PASS: absolute z=-50 and relative depth=50 retain identical nonzero motion below terrain z=-100')


def run_startup(args):
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
                results[directory.name] = run_startup_case(args, binary, fixture, directory,
                                              second=count == 8)
            for filename in (("receiver.h5", "second.h5") if count == 8 else ("receiver.h5",)):
                peak = compare(args.work_dir / ("baseline-" + str(count)) / filename,
                               args.work_dir / ("candidate-" + str(count)) / filename, count)
                print("PASS:", count, filename, "metadata and motion match; peak", peak)
        for name, count, malformed in (("missing", 0, False), ("bad", 0, True),
                                       ("bad-after-batch", 1024, True),
                                       ("bad-size", 0, False), ("bad-size-after-batch", 1024, False)):
            fixture = args.work_dir / (name + ".h5")
            if name != "missing":
                template(fixture, count, malformed)
                if name.startswith("bad-size"):
                    with h5py.File(fixture, "a") as stream:
                        group = stream.create_group("z_bad")
                        group["STX,STY,STZ"] = np.asarray([800, 800, 0], dtype=np.float64)
                        group["WINDOWS"] = np.ones(5)
            if name == "missing":
                expected_error = f"Error opening file [{fixture}]"
            elif malformed:
                expected_error = "Error reading rechdf5 station z_bad dataset STX,STY,STZ"
            else:
                expected_error = "Invalid rechdf5 station z_bad dataset WINDOWS: expected 4 elements, found 5"
            results[name] = run_startup_case(args, args.candidate_sw4, fixture,
                                args.work_dir / name, expected_error=expected_error)
            print("PASS:", name, "fails before time stepping")
    finally:
        (args.work_dir / "results.json").write_text(json.dumps(results, indent=2) + "\n")


def run_depth(args):
    args.work_dir=args.work_dir.resolve()
    args.work_dir.mkdir(parents=True,exist_ok=False)
    template=args.work_dir/'stations.h5'
    with h5py.File(template,'w') as stream:
        for name,use_z in [('s_absolute',1),('s_relative',0)]:
            group=stream.create_group(name)
            group['STX,STY,STZ']=np.asarray([1800,1800,300],dtype=np.float64)
            group['ISNSEW']=np.asarray([0],dtype=np.int32)
            group['USEZVALUE']=np.asarray([use_z],dtype=np.int32)
    if args.unsupported_baseline:
        baseline=args.unsupported_baseline.resolve(strict=True)
        inline_only="\n".join(line for line in DEPTH_INPUT.splitlines()
                               if 'file=s_relative ' not in line and not line.startswith(('rechdf5 ','topography ')))+"\n"
        control=run_depth_case(args,baseline,args.work_dir/'baseline-flat-control',template,input_text=inline_only)
        with h5py.File(control/'inline.h5','r') as stream:
            assert any(np.max(np.abs(stream['s_absolute'][component][...]))>0 for component in 'XYZ')
        failed=run_depth_case(args,baseline,args.work_dir/'baseline',template,expect_topography_failure=True)
        if failed is not None:
            check_depth(failed,unsupported=True)
        print('PASS: baseline reproduces unsupported depth behavior; flat absolute inline control works')
    check_depth(run_depth_case(args,args.sw4.resolve(strict=True),args.work_dir/'candidate',template))
    print('PASS: template absolute/relative depths match inline receivers over Gaussian topography')
    elevated_absolute(args, args.sw4.resolve(strict=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    launcher = argparse.ArgumentParser(add_help=False)
    launcher.add_argument("--work-dir", type=Path, required=True,
                          help="New shared directory for run artifacts")
    launcher.add_argument("--launcher", default="srun")
    launcher.add_argument("--job-id", help="Existing Slurm allocation ID")
    launcher.add_argument("--gpu", action="store_true", help="Request one GPU per rank")

    startup = modes.add_parser("startup", parents=[launcher],
                               help="Compare waveforms and verify metadata failures")
    startup.add_argument("--baseline-sw4", type=Path, required=True)
    startup.add_argument("--candidate-sw4", type=Path, required=True)
    startup.add_argument("--tasks", type=int, default=2)
    startup.add_argument("--cpus-per-task", type=int, default=4)
    startup.set_defaults(check=run_startup)

    depth = modes.add_parser("depth", parents=[launcher],
                             help="Check absolute and relative receiver depths with two ranks")
    depth.add_argument("--sw4", type=Path, required=True)
    depth.add_argument("--unsupported-baseline", type=Path,
                       help="Optional older build expected to lack depth support")
    depth.set_defaults(check=run_depth, tasks=2, cpus_per_task=4)

    args = parser.parse_args()
    args.check(args)


if __name__ == "__main__":
    main()
