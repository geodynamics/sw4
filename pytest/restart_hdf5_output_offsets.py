#!/usr/bin/env python3
"""Compare receiver and SSI output after two restarts with an uninterrupted run.

Run inside a one-node GPU allocation with a GPU/HDF5-enabled SW4 executable.
"""

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

import h5py
import numpy as np


INPUT = """grid h=200 x=30000 y=30000 z=17000 extrapolate=1
time steps={steps}
fileio path=.
supergrid gp=30
block vp=6000 vs=3464 rho=2700
block vp=4000 vs=2000 rho=2600 z2=1000
source x=11000 y=11000 z=200 mxy=1e18 t0=0.01 freq=16.6667 type=Gaussian
rec x=11000 y=11000 z=0 file=sta01 usgsformat=0 sacformat=0 hdf5format=1 hdf5file=receiver.h5 writeEvery=2 downSample={downsample}
rec x=11200 y=11000 z=0 file=sta02 usgsformat=0 sacformat=0 hdf5format=1 hdf5file=receiver.h5 writeEvery=2 downSample={downsample}
rec x=11000 y=11000 z=0 file=vel_xyz usgsformat=0 sacformat=0 hdf5format=1 hdf5file=velocity.h5 writeEvery=2 downSample={downsample} variables=velocity nsew=0
rec x=11200 y=11000 z=0 file=vel_enu usgsformat=0 sacformat=0 hdf5format=1 hdf5file=velocity.h5 writeEvery=2 downSample={downsample} variables=velocity nsew=1
ssioutput file=ssi-restart xmin=10000 xmax=12000 ymin=10000 ymax=12000 depth=0 precision=4 bufferInterval={buffer_interval} dumpInterval={dump_interval}{ssi_options}
checkpoint cycleInterval={checkpoint_interval} restartpath=. file=restart hdf5=yes{restartfile}
"""

RECEIVER_COMPONENTS = {
    "receiver.h5": {"sta01": ("X", "Y", "Z"), "sta02": ("X", "Y", "Z")},
    "velocity.h5": {
        "vel_xyz": ("Vx", "Vy", "Vz"),
        "vel_enu": ("Vew", "Vns", "Vup"),
    },
}


def run(sw4, input_file, run_dir, tasks):
    command = ["srun", "-N", "1", "-n", str(tasks), "-c", "1",
               "--gpus-per-task=1", "--gpu-bind=closest", str(sw4),
               str(input_file)]
    result = subprocess.run(command, cwd=str(run_dir), text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=180)
    if result.returncode != 0 or "HDF5-DIAG" in result.stdout:
        raise AssertionError("SW4 failed for {} (exit {}):\n{}".format(
            input_file.name, result.returncode, result.stdout[-8000:]))


def capture(run_dir, steps, downsample, dump_interval, read_ssi_data=True):
    receiver = {}
    extent = steps // downsample + 1
    for filename, stations in RECEIVER_COMPONENTS.items():
        # The centered velocity derivative is only saved during time stepping.
        samples = steps if filename == "velocity.h5" else steps + 1
        expected_npts = (samples - 1) // downsample + 1
        with h5py.File(run_dir / filename, "r") as stream:
            for station, components in stations.items():
                group = stream[station]
                npts = int(group["NPTS"][0])
                if npts != expected_npts:
                    raise AssertionError("{} NPTS {} != {}".format(
                        station, npts, expected_npts))
                for component in components:
                    if group[component].shape != (extent,):
                        raise AssertionError("{} {} shape {} != {}".format(
                            station, component, group[component].shape, (extent,)))
                receiver[station] = {
                    "NPTS": npts,
                    "components": {
                        component: group[component][:npts].copy()
                        for component in components
                    },
                }

    expected_cycles = list(range(dump_interval, steps + 1, dump_interval))
    if expected_cycles[-1] != steps:
        expected_cycles.append(steps)
    with h5py.File(run_dir / "ssi-restart.ssi", "r") as stream:
        progress = (int(stream["lastsw4timestep"][0]),
                    int(stream["lastoutputindex"][0]))
        expected_progress = (steps, len(expected_cycles) - 1)
        if progress != expected_progress:
            raise AssertionError("SSI progress {} != {} for steps {}".format(
                progress, expected_progress, steps))
        ssi = {}
        for component in range(3):
            name = "vel_{} ijk layout".format(component)
            dataset = stream[name]
            if dataset.shape[0] != len(expected_cycles):
                raise AssertionError("SSI {} has {} slots, expected {}".format(
                    name, dataset.shape[0], len(expected_cycles)))
            if dataset.id.get_storage_size() == 0:
                raise AssertionError("SSI {} has no stored data".format(name))
            if read_ssi_data:
                ssi[name] = dataset[:].copy()
    return receiver, ssi


def damage_uncheckpointed_tail(run_dir, checkpoint, downsample,
                               dump_interval, poison_ssi=True):
    for filename, stations in RECEIVER_COMPONENTS.items():
        samples = checkpoint if filename == "velocity.h5" else checkpoint + 1
        first = (samples - 1) // downsample + 1
        with h5py.File(run_dir / filename, "r+") as stream:
            for station, components in stations.items():
                group = stream[station]
                npts = int(group["NPTS"][0])
                if first >= npts:
                    raise AssertionError("No receiver tail after checkpoint")
                for component in components:
                    group[component][first:npts] = 1234567.0
                group["NPTS"][...] = first

    with h5py.File(run_dir / "ssi-restart.ssi", "r+") as stream:
        first = checkpoint // dump_interval
        for component in range(3):
            dataset = stream["vel_{} ijk layout".format(component)]
            if first >= dataset.shape[0]:
                raise AssertionError("No SSI tail after checkpoint")
            if poison_ssi:
                dataset[first:] = 1234567.0
        stream["lastsw4timestep"][...] = first * dump_interval
        stream["lastoutputindex"][...] = first - 1


def compare(actual, expected, label, ssi_atol):
    actual_receiver, actual_ssi = actual
    expected_receiver, expected_ssi = expected
    for station in expected_receiver:
        got = actual_receiver[station]
        want = expected_receiver[station]
        if got["NPTS"] != want["NPTS"]:
            raise AssertionError("{} {} NPTS {} != {}".format(
                label, station, got["NPTS"], want["NPTS"]))
        for component in want["components"]:
            np.testing.assert_allclose(
                got["components"][component], want["components"][component],
                rtol=2e-5, atol=1e-8,
                err_msg="{} {} {}".format(label, station, component))
    for name in expected_ssi:
        np.testing.assert_allclose(actual_ssi[name], expected_ssi[name],
                                   rtol=2e-5, atol=ssi_atol,
                                   err_msg="{} {}".format(label, name))


def check_case(sw4, tasks, steps, downsample, dump_interval,
               buffer_interval, checkpoints, root, ssi_options=""):
    run_dir = root / "steps{}-downsample{}-dump{}-buffer{}".format(
        steps, downsample, dump_interval, buffer_interval)
    run_dir.mkdir()
    base = run_dir / "base.in"
    base.write_text(INPUT.format(steps=steps, downsample=downsample,
                                 dump_interval=dump_interval,
                                 buffer_interval=buffer_interval,
                                 ssi_options=ssi_options,
                                 checkpoint_interval=checkpoints[0],
                                 restartfile=""))
    run(sw4, base, run_dir, tasks)
    read_ssi_data = not ssi_options
    expected = capture(run_dir, steps, downsample, dump_interval, read_ssi_data)
    for name, station in expected[0].items():
        traces = np.stack(list(station["components"].values()))
        if not (np.isfinite(traces).all() and np.any(traces)):
            raise AssertionError("Receiver {} has no finite nonzero motion".format(name))
    ssi_peak = (max(np.max(np.abs(values)) for values in expected[1].values())
                if read_ssi_data else None)
    if read_ssi_data and not (np.isfinite(ssi_peak) and ssi_peak > 0):
        raise AssertionError("SSI fixture has no finite nonzero motion")

    for checkpoint in checkpoints:
        restart = run_dir / "restart-{}.in".format(checkpoint)
        checkpoint_name = "restart.cycle={:0{}d}.sw4checkpoint".format(
            checkpoint, len(str(steps)))
        restart.write_text(INPUT.format(
            steps=steps, downsample=downsample,
            dump_interval=dump_interval, buffer_interval=buffer_interval,
            ssi_options=ssi_options,
            checkpoint_interval=checkpoints[0],
            restartfile=" restartfile={}".format(checkpoint_name)))
        if not (run_dir / checkpoint_name).is_file():
            raise AssertionError("Missing checkpoint {}; files: {}".format(
                checkpoint, sorted(p.name for p in run_dir.iterdir())))
        damage_uncheckpointed_tail(run_dir, checkpoint, downsample,
                                   dump_interval, read_ssi_data)
        run(sw4, restart, run_dir, tasks)
        compare(capture(run_dir, steps, downsample, dump_interval, read_ssi_data), expected,
                "steps {} downsample {} restart {}".format(
                    steps, downsample, checkpoint),
                0.02 if ssi_options else 1e-8)
    print("PASS: steps={}, receiver downSample={}, SSI dumpInterval={}, "
          "bufferInterval={}, restarts {} and {}{}{}".format(
              steps, downsample, dump_interval, buffer_interval,
              *checkpoints, ssi_options,
              " (SSI metadata only)" if ssi_options else ""))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sw4", required=True, type=Path)
    parser.add_argument("--tasks", type=int, default=2)
    args = parser.parse_args()
    sw4 = args.sw4.resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="sw4-hdf5-restarts-") as name:
        root = Path(name)
        check_case(sw4, args.tasks, 11, 2, 3, 2, (4, 8), root)
        check_case(sw4, args.tasks, 12, 1, 3, 2, (5, 10), root)
        check_case(sw4, args.tasks, 19, 3, 4, 3, (7, 14), root,
                   " zfp-accuracy=1e-2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
