# Production review corrections

The review in `review-unified-backends-62825bc7/REPORT.md` applies to commit
`62825bc7`. Corrections below are undergoing fresh validation; the earlier
validation record does not establish acceptance of these changes.

## Source preparation and interface injection

Source preparation now retains raw discrete histories until optional filtering,
then performs one spline conversion. Prepared copies retain their representation.
This corrects the inherited native discrete-force double conversion and the GPU
force constructor's uninitialized readiness flag and mismatched integer storage.
All three force histories and all six moment histories are filtered separately
using the existing scalar filter before spline conversion. GPU forcing evaluates
all those histories. Cartesian interface source partitioning follows the native
reference's one-third/two-thirds weights.

The user approved the review fixes and the additional multicomponent filtering
correction. These changes can alter affected simulations; compare earlier results
with that distinction in mind.

Build the on-demand native probe with:

```sh
cmake --build build --target sw4_source_check
srun -n 2 build/bin/sw4_source_check case.in
```

Use a Cartesian, unrefined input that contains `(1800, 1800, 1000)` inside its
physical domain. The probe links the real solver. It checks scalar forces,
three-component forces, and six-component moments, with and without filtering,
including copies before/after preparation and repeated preparation. Volume sums
check forces; all nine first spatial moments check the six independent moment
histories against the raw-sample/filter oracle.

Fresh two-rank OpenMP execution on Perlmutter passed all six combinations, with
maximum absolute error `2.7e-15` (limit `2e-10`). Evidence:
`/pscratch/sd/h/houhun/sw4-review-fixes.l5lzya/source-check/run-v4.log`.
The combined CUDA/native-inversion configuration also compiles. CUDA propagation,
interface sweeps, broader regressions and performance are still pending.

These additional checks are excluded from default builds and are not registered
with pytest or CTest. The existing pytest suites are unchanged.

## Receiver format consistency

SAC UTC stores milliseconds; SW4 internally stores microseconds. SAC now writes
the millisecond portion into `NZMSEC` and includes the remaining microseconds in
`B`, `E`, and `O`. Reading SAC converts milliseconds back to microseconds. USGS
text writes a six-digit fractional second, including leading zeroes. Samples and
the simulation's time integration are unchanged by these metadata corrections.

Rechdf5 component inclination and azimuth now have the correct names and match
SAC's Cartesian or geographic component conventions, including derivative
output modes. Geographic EW/NS azimuths are 90/0 degrees. The user approved both
metadata corrections because they can affect downstream timing/rotation.

```sh
cmake --build build --target sw4_receiver_check
python tests/backends/receiver_formats.py --backend OPENMP \
  --sw4 build/bin/sw4 --reader build/bin/sw4_receiver_check \
  --work-dir "$SCRATCH/sw4-formats"
```

The on-demand driver checks displacement, velocity, divergence, curl, strains,
and displacement gradient; Cartesian/geographic components; flat/Gaussian
topography; HDF5 downsampling by one/three; partial writes and the final sample;
UTC fractions `.234567`, `.000123`, and `.999999`; origin/sample timing and
component metadata. It rejects nonfinite samples and missing signals.
SAC/rechdf5 float32 samples are compared with USGS solver-precision samples using
`1e-12 + 2e-6 * component_peak`. Displacement/velocity use metres/metres per
second, matching USGS headings and HDF5's `UNIT` attribute.

Downsampling is an HDF5-only option in the existing documented contract. A
downsampled HDF5 trace is checked against the corresponding SAC/USGS samples;
full-rate outputs have matching counts and samples. SAC `B` must be interpreted
together with its reference UTC, rather than compared alone with USGS times.

Fresh OpenMP execution passed all 22 format cases. Four displacement fixtures
also passed read-back through SW4's real SAC, USGS and rechdf5 readers, with
start-time/sample-interval errors below `1e-8` seconds and matching vertical
samples within the same floating-point tolerance. Evidence:
`/pscratch/sd/h/houhun/sw4-review-fixes.l5lzya/formats-cpu-v6/`.
CUDA format execution remains pending allocation.

## Cartesian anisotropic residuals

The inherited native Cartesian C operator computed but did not store the three
residual components in its interior and bottom boundary closure. Those rows now
store the existing residuals with the same `1/h²` scaling as the top closure.
The user explicitly approved this scientific correction. Affected anisotropic
runs can differ substantially from earlier results.

```sh
cmake --build build --target sw4_anisotropic_check
srun -n 2 build/bin/sw4_anisotropic_check case.in
```

This on-demand probe calls the real operator with each quadratic monomial in
all three displacement components. It checks every physical vertical row,
including both closures, against the independent continuum tensor contraction
`C_(i,a,l,b) * d_a d_b u_l`. Materials include an isotropic tensor and a coupled,
positive definite anisotropic tensor. The unfixed operator failed the initial
isotropic test with maximum error 14. The corrected operator passed the broader
check with maximum error `8.5e-12` (limit `2e-10`). Evidence:
`/pscratch/sd/h/houhun/sw4-review-fixes.l5lzya/aniso-operator-general.log`.

The native propagation check also passed the isotropic limit, including a remote
receiver. It uses a 100 m grid and eight supergrid points. Coarser initial inputs
produced growing amplitudes and are retained as failed evidence, not acceptance.
The production timestep/default grid policies were not changed. Evidence:
`/pscratch/sd/h/houhun/sw4-review-fixes.l5lzya/review-cpu-fine/`.
CUDA propagation and anisotropic topography remain additional acceptance gates.

## Remaining review corrections and acceptance gates

Sfile setup now communicates Vp/Vs and Q arrays through host packing while those
arrays reside in host storage. Curvilinear material extrapolation uses the native
loop, including attenuation fields. GPU anisotropic setup skips undefined
isotropic refinement/material arrays and waits for the GPU stream before host
operators read solution data. Both material checks use the native Vp/Vs lower
limit `sqrt(4/3)`. SAC source parsing reads lengths before allocating histories.
Both backends complete all four requests in the curvilinear setup exchange.
Event defaults follow native behavior; unsafe GPU parallel events are rejected.

SZ checkpoint compression now uses each grid's actual one-dimensional dataset
length and configured precision, releases filter metadata, and checks encoding
availability. CMake requires both SZ and the HDF5 filter and verifies linking.
SZ runtime validation is unavailable on this system because those dependencies
are absent; configuration rejects that missing dependency explicitly.

Validation drivers reject nonconvergence, nonfinite output and noise-only
references. Interface comparisons require meaningful late arrivals on both
grids. Performance records include the material-model SHA-256; resumed results
require valid logs and unchanged receiver checksums. The negative acceptance
checks passed (`python tests/backends/check_acceptance.py`).

Fresh OpenMP results include the existing solver suite (20 passes, five existing
skips), four existing HDF5 regressions, and 22 source/interface/material/event
propagation cases, including the anisotropic isotropic limit. The default
standalone OpenMP build succeeded with RAJA and Umpire discovery disabled.
Fresh CUDA with native inversion compiled successfully. These compile results
do not substitute for CUDA runtime acceptance. Evidence is retained under
`/pscratch/sd/h/houhun/sw4-review-fixes.l5lzya/`.

Additional curvilinear physical-source fixtures currently fail their convergence
gate on CPU. Their failures are retained in `curvi-controls/`, `curvi-thicker/`
and `curvi-constant/`. They are not counted as passing extrapolation validation.
GPU runtime, four-node performance, the current inversion rerun, velocity
receiver read-back, and ROCm/HIP execution remain open gates.

## Extended receiver checks and queued execution

The forward format fixtures now use the same stable 100 m/eight-point supergrid
configuration used for the source checks, with `cfl=0.8`. All 22 OpenMP forward
format cases passed. SAC and rechdf5 stored samples agree exactly for these
fixtures; comparisons permit floating-point differences and compare USGS text
using the documented component-relative tolerance. Evidence:
`/pscratch/sd/h/houhun/sw4-review-fixes.l5lzya/formats-cpu-fine-forward/`.

The optional real-reader probe now checks all three components and selects
displacement/velocity and Cartesian/geographic filenames. Geographic
displacement round-trip passed, including horizontal components. Two newly
exposed inherited failures await scientific approval rather than being hidden
by a looser threshold:

- Rechdf5 velocity reading searches for displacement datasets and fails.
  `velocity-read-before.log` retains the reproduced failure. GPU SAC restart
  also selects Cartesian displacement filenames regardless of mode/orientation.
- Cartesian SAC observation loading applies a geographic conversion to grid
  components; a tested horizontal sample changed by about `1e-6 m` from
  `2e-3 m`. USGS/rechdf5 load the original grid components.
  `displacement-xyz-read-difference.log` retains the reproduced failure.

The on-demand probe therefore intentionally fails affected Cartesian/velocity
read-back checks until their corrections are approved and validated. Forward
output equivalence is distinct from observation/restart equivalence. Native
SAC restart's raw copy of geographic components and rechdf5's handling of
`ignore_utc` also need follow-up inspection/validation.

The fixed source regression inputs additionally reject displacement above
1000 m, a deliberately generous sanity ceiling for their fixed source strengths.
This prevents agreement between two growing solutions from being counted as a
pass; it is not a general amplitude limit for SW4 simulations.

Queued CUDA jobs use frozen executable copies with SHA-256 checks:
`59219272` runs the unchanged solver suite, forward format checks and the
converged source/material comparisons. `59219270` runs three alternating pairs
of the exact RAJA hmr3 case on four nodes, followed by two-node HDF5 regressions.
Pending state is not runtime or performance acceptance. The earlier two queued
GPU jobs were replaced before execution to separate these workloads.
