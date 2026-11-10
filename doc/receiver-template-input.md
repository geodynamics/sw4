# Receiver-template input

This checkout targets `raja`. The shared input contract below applies to all three companion changes: [native CPU (#251)](https://github.com/geodynamics/sw4/pull/251), [standalone RAJA (#252)](https://github.com/geodynamics/sw4/pull/252), and [unified CPU/GPU (#253)](https://github.com/geodynamics/sw4/pull/253). Branch-specific sections describe those companion implementations; they do not add another backend to this checkout.

`rechdf5` and `sachdf5` metadata are read once by rank zero of each participating simulation MPI group. SW4 distributes
ordered station records in batches of up to 1,024, then performs the existing
coordinate conversion, domain filtering, and receiver construction on every
rank. Temporary distribution memory is bounded by a batch rather than the full
template. No input schema or command-line option changes are required.

The packed stream contains names, coordinates, coordinate representation,
orientation, and effective observation-window settings. Absolute depth (`USEZVALUE`) and four-value
`WINDOWS` retain the native CPU interpretation on every backend. Dataset reads and station-group traversal occur only on that group root.
Distribution uses the simulation communicator, so parallel events cannot mix
their receiver metadata or wait for unrelated event ranks.
Read or close failures stop all ranks before time stepping, including failures
after a completed batch. A missing template is an input error.

Receiver output grouping, geometry, sampling, output locations, and owner-rank
observation waveform reads retain their existing behavior. This does not move
receiver outputs or checkpoint histories into node-local storage.

Node-local `sbcast` staging is optional for this reader. One shared-file read is
the normal input path. The reader does not require independent ranks to open
separate replicas through a global MPI-IO communicator.

Rank zero reports metadata-read, distribution, and receiver-construction times
separately. The read time includes file open/close and HDF5 traversal, excluding
the timed distribution and construction phases. These are root wall-clock
measurements; they are not summed CPU time or a scale-performance guarantee.
The existing per-receiver topography/ownership MPI collectives remain in place.

## Bounded metadata regression

Use the repository's supported Python environment with h5py and NumPy. This
harness compiles the actual reader against real HDF5 with a replayed MPI stream
and a receiver-construction stub. It compares station records against a supplied
baseline revision, checks that only root opens/traverses the candidate template,
and checks failure propagation. It does not validate real MPI progress or SW4
numerical outputs.

```sh
: "${EQSIM_PYTHON:?select the supported scientific Python}"
"$EQSIM_PYTHON" pytest/receiver_template_metadata.py \
  --hdf5-prefix "$HDF5_ROOT" --mpi-include "$MPI_INCLUDE" \
  --baseline-ref <base-commit>
```

## Native regression in an allocation

Build HDF5-enabled baseline and candidate binaries with the same backend,
compiler options, and dependencies. Run only inside an authorized allocation;
select a shared scratch output path that does not already exist.

```sh
"$EQSIM_PYTHON" pytest/receiver_template_startup.py \
  --baseline-sw4 /path/to/baseline/sw4 --candidate-sw4 /path/to/candidate/sw4 \
  --work-dir "$SCRATCH/receiver-template-regression" --tasks 2
```

Add `--gpu` for a GPU build and allocate at least one GPU per MPI rank. The
regression preserves logs and a JSON result file. It compares receiver metadata
and finite, nonzero waveforms for 8 and 2,053 stations, exercises repeated station
names in separate output files, crosses batch boundaries, and checks missing or
malformed inputs before and after a full batch. Tests use the existing step-zero
and downsampling contract. The numerical comparison tolerance is `rtol=2e-12`,
`atol=1e-15` for matched builds; it is not a cross-backend acceptance threshold.

Startup timestamps observed from solver output include output buffering. These
small regression cases do not replace repeated, representative-scale comparisons
with matched receiver geometry, rank layout, inputs, and cache conditions.

## Observation-window parity

Observation windows use the native CPU five-sample tanh taper. Four-value
`WINDOWS` selects two intervals and takes their maximum weight where they
overlap. The second interval follows the existing native CPU activation rule
(at least one bound is positive). `WindowL`/`WindowR` selects one interval;
the second interval is disabled initially and when switching back from two
intervals. This fixes the formerly unbounded default second window on native
CPU, which could override a requested single window.

The standalone RAJA and unified GPU changes use this weighting instead of their previous cosine tapers. Its retained weighted-product helper uses the same window
weight as its misfit and cross-correlation paths. These are deliberate changes to windowed
misfit calculations; forward waveform comparisons alone do not establish
inversion equivalence. The production setters and misfit weighting blocks are
exercised with single intervals, disjoint/overlapping pairs, and mode changes:

```sh
"$EQSIM_PYTHON" pytest/receiver_windows.py
```

For the standalone RAJA metadata regression, use a native CPU baseline revision to check
new `USEZVALUE` and `WINDOWS` support. Its previous RAJA baseline intentionally
lacks these fields; native CPU and unified already supported them. The native startup regression compares matched backend
builds and exercises input, output, sampling, and failure behavior separately.

## Absolute and surface-relative depth control

`pytest/receiver_template_depth.py` compares template receivers with inline
absolute and surface-relative receivers over Gaussian topography. Its positive-depth fixture
requires the two placements to produce different waveforms, so a test cannot
pass merely because both requested depths sample the same grid point.

```sh
"$EQSIM_PYTHON" pytest/receiver_template_depth.py --sw4 /path/to/candidate/sw4 \
  --work-dir "$SCRATCH/receiver-depth-regression"
```

Use `--gpu` for CUDA execution. The standalone RAJA change also corrects the Gaussian interpolator's
success status: computing a valid analytical height must not report failure.
The optional `--unsupported-baseline /path/to/old/sw4` control records the old
RAJA failure or absolute-depth mismatch while retaining a valid flat-terrain inline control.

The metadata harness also accepts `--single-precision` to exercise float receiver
state against a double-precision reference. Scalar HDF5 window bounds are read
into double storage before conversion, so an eight-byte HDF5 read cannot overwrite
a neighboring float field.

The startup regression defaults to four CPUs per MPI rank, with one OpenMP
thread. Allocate those CPUs explicitly; record the resource shape when comparing reader and
receiver-construction costs, which can differ substantially between builds.

## Independent simulation groups

`pytest/receiver_template_groups.py` runs the production reader with real MPI
and HDF5 on four ranks split into two groups, using different receiver files
and batch counts. Receiver construction is a stub that performs a real group
collective. The check verifies group isolation and one input read per group;
it does not replace a multi-event inversion regression.

```sh
"$EQSIM_PYTHON" pytest/receiver_template_groups.py \
  --hdf5-prefix "$HDF5_ROOT" --mpi-prefix "$MPI_ROOT" \
  --work-dir "$SCRATCH/receiver-group-regression"
```

On the unified branch, also run `receiver_windows.py --raja`. The helper first
preprocesses the selected implementation (including inline header methods), so
the GPU checks cannot accidentally exercise the native CPU window code.
Smoke startup timestamps begin at `srun` and exclude module setup/allocation
wait; they do not establish representative-scale performance.

Metadata dataset element counts are checked before fixed-buffer reads: one for
scalar flags/bounds, three for coordinates, and four for `WINDOWS`. Scalar and
one-element-array representations remain accepted. Invalid counts use the same
metadata-failure path, including after a complete batch.

The local runner controls need no allocation:

```sh
"$EQSIM_PYTHON" pytest/test_receiver_template_startup_timeout.py
```

They cover normal exit, a launcher ignoring SIGTERM, and a surviving descendant
holding output open. Native runs retain output in a file and terminate the
launcher's process group with a bounded grace period on timeout. First-step
observation polls at 0.1-second intervals and includes solver output buffering.

The depth regression also requires an absolute receiver at
`z=-50 m` and its equivalent `depth=50 m` receiver beneath the Gaussian summit
at `z=-100 m` to exist and produce identical nonzero motion. Absolute receivers
use the actual terrain for validity checks without adding it to their z value.
