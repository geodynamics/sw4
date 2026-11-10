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
and checks failure propagation. Scalar window bounds are checked against the
station and command inputs. The harness does not validate real MPI progress or SW4
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
and downsampling contract. Negative cases require the expected SW4 metadata error
before time stepping. The numerical comparison tolerance is `rtol=2e-12`,
`atol=1e-15` for matched builds; it is not a cross-backend acceptance threshold.

Startup timestamps observed from solver output include output buffering. These
small regression cases do not replace repeated, representative-scale comparisons
with matched receiver geometry, rank layout, inputs, and cache conditions.

## Observation-window parity

Observation windows use the native CPU five-sample tanh taper. Four-value
`WINDOWS` selects two intervals and takes their maximum weight where they
overlap. Scalar `WindowL` and `WindowR` values override the command bounds only for
the station that supplies them. A station's four-value `WINDOWS` takes
precedence over its scalar bounds.
The second interval follows the existing native CPU activation rule
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

## Receiver-objective numerical checks

Use the CPU-only implementation as the behavioral reference. CPU corrections in
this change have independent evidence: a closed-form linear-waveform case for
normalized L2 derivatives, an exactly known Y-only correlation with X excluded,
and an AddressSanitizer endpoint overread. The L2 objective definition and CPU
window taper are retained. Its shift derivatives now differentiate the returned
normalized score, including the shift-dependent observation-energy denominator.
The positive derivative-squared curvature term `dd1shift` has the same normalized
units. Observation component masks apply independently, and interpolation edge
stencils fit zero-based sample storage. GPU paths retain the CPU observation
support rule, iteration limit and zero-curvature adjoint-source fallback. If a Newton
step has zero curvature or a nonfinite result, the search keeps its best
sampled shift and identifies the station in the log. Correlation trials
outside observation support contribute zero.

```sh
: "${EQSIM_PYTHON:?select the supported scientific Python}"
"$EQSIM_PYTHON" pytest/receiver_objective_numerics.py \
  --work-dir /new/path/objective-check --asan \
  --cpu-reference /path/to/cpu-checkout
```

Add `--raja` to select the unified GPU preprocessor path. The harness compiles the
complete production misfit, cross-correlation, maximum search, adjoint-source
and interpolation methods with controlled construction/storage. It checks both
normalization modes and sample representations, single/two/overlapping windows,
component exclusion, finite-difference first/second derivatives, and the dot
product of the returned data-adjoint source with a synthetic-data perturbation.
Analytic controls and exact-size endpoint buffers provide independent checks.

For native validation, first build the selected HDF5-enabled backend normally.
Within an authorized allocation, with that build's compiler environment, create
a separate executable linked against its actual production objects:

```sh
"$EQSIM_PYTHON" pytest/receiver_objective_native_build.py \
  --build-dir /path/to/existing-cmake-build --work-dir /new/path/native-check
```

For standalone RAJA's Make build, use `--make-optdir optimize` instead of
`--build-dir` (or the actual output-directory name). The builder prints compiler
and link commands and generates a main that runs the receiver tests after normal
SW4 setup. Launch `objectives-native` with the normal backend MPI/GPU resource
shape and a small, flat uniform-material input with (1800, 1800, 0) inside the grid.
The harness requires all 20 cases to run. Two MPI ranks were tested.
Native checks use real EW/TimeSeries objects, double and HDF5-float sample storage,
window/component cases, unequal sampling, endpoint data and zero observations.

These are receiver-objective and data-adjoint directional checks. Travel-time
finite differences require an isolated interior correlation maximum; zero-data
controls check finite CPU-compatible fallback behavior. Full wave-equation
adjoint propagation, complete model inversion, HIP execution and representative
startup-performance benchmarking are separate validation scopes.
