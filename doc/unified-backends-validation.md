# Perlmutter validation, 2026-10-02

This is a historical record of the initial unification. The later review fixes,
supported receiver contracts and current qualification status are recorded in
[the updated review corrections](unified-backends-updated-review.md).

**Review correction:** the seven-case waveform comparator reported a pass for
the small Cartesian refinement fixture even though both solvers printed severe
interface nonconvergence. That fixture is not a valid correctness acceptance
result. See `review-unified-backends-62825bc7/REPORT.md` finding F01. The independent
existing solver regressions and full hmr3 measurements remain separate evidence.
Review fixes and new acceptance runs are in progress.

The native OpenMP and CUDA configurations pass the checks below. HIP compilation,
correctness and performance remain to be tested on Frontier. The existing
`pytest/` and `pytest-sw4mopt/` sources are unchanged from `fix-sfile-srf`.

The [machine-readable record](validation/perlmutter-2026-10-02.json) contains
executable/input identities, allocation nodes, commands, timings, tolerances and
waveform differences. Logs and outputs are retained under
`/pscratch/sd/h/houhun/sw4-unified.4MtlHr`.

## Correctness and build checks

| Check | Result | Retained evidence |
| --- | --- | --- |
| Native level-zero suite | 20 passed, 5 optional cases skipped | `suite-cpu-v4.log` |
| CUDA level-zero suite | 20 passed, 5 optional cases skipped | `suite-cuda-release.log` |
| Native inversion | Gradient, Hessian, misfit curve and one-point inversion passed | `suite-mopt-v4.log`, `mopt-final/batch.log` |
| Native HDF5 regressions | 4 passed | `hdf5-cpu-final.log` |
| CUDA HDF5 regressions | 4 passed | `hdf5-cuda-release.log` |
| CPU/CUDA waveform comparisons | All 7 cases passed | `comparison-release.log`, `comparison-release/comparison.json` |
| Text/HDF5 SRF consistency | Passed on both backends | `comparison-release.log` |
| Default native build | Built and ran with OpenMP; no GPU libraries or package discovery | `build-default-final.log`, `energy-default-final/run.log` |
| CUDA compatibility probe | Passed | `cuda-compat-release.log` |
| Fresh CUDA configuration | Host compiler selected automatically; optional transport link generated | `configure-transport-v2.log` |
| Cray CUDA MPI transport | GPU-aware MPI refinement passed and agreed with staged execution | `energy-transport-probe/run.log`, `transport-comparison.json` |

The seven waveform cases cover inline material, topography, Cartesian mesh
refinement, attenuation, text SRF, HDF5 SRF and Sfile material. They require finite,
nonzero traces and check receiver metadata. Waveform acceptance uses
`max_difference <= 1e-12 + 2e-5 * reference_peak` for each component; byte equality
is not required. Existing solver checks retain their reference tolerances.
The inversion one-point case was rerun successfully in CPU job `59203804` after
an earlier allocation expired; the four-case result combines the successful
individual checks.

The release CUDA executable was compiled from solver sources at `5055a1a6`.
`e55076b4` changes configuration and profiling, with no changes under `src/`.
Its fresh configuration and optional transport link were checked separately;
the transport executable reuses the release objects with the generated CMake
link command. Native inversion remains an OpenMP executable.

## Four-node performance

The baseline is an executable built from exact `raja` commit `52b51c92`, rather
than a pre-existing binary of uncertain provenance. Both executables ran the
complete nine-second `raja:performance/large/hmr3.in`, with only model and output
paths relocated. The checked-in input is byte-identical to that branch's input.
The model is the same `USGSBayAreaVM-08.3.0-corder.rfile` for all runs.

Allocation `59206992` used `nid[001288-001289,001292-001293]`: four nodes,
16 A100 GPUs, 16 MPI ranks, four CPU cores per rank and one OpenMP thread.
Three runs per executable alternated in order. MPI buffers were staged and
`MPICH_GPU_SUPPORT_ENABLED=0`.

| Executable | Solver times, seconds | Median, seconds |
| --- | --- | --- |
| Exact RAJA baseline | 82.47, 82.16, 82.13 | 82.16 |
| Unified CUDA release | 79.36, 79.65, 79.65 | 79.65 |

The unified median is **3.06% faster**, passing the requirement of at most 5%
slowdown. All three paired station-waveform comparisons pass the stated
tolerance. These are measured results for this input, allocation and toolchain.

Both builds use CUDA 13.2, GNU 13.2.1, Cray MPICH 9.1, RAJA 2026.07.0 and
Umpire 2026.07.1, with PROJ, parallel HDF5, ZFP and MPI FFTW enabled. Both omit
host OpenMP for CUDA. The legacy Make baseline enables CUDA fast math and pool
prefetch; the unified build uses ordinary floating-point compilation and leaves
pool prefetch off. Per-array prefetch is off for both.

The apparent Cartesian-refinement stall was reproduced with host OpenMP enabled,
`OMP_PROC_BIND=spread`, `OMP_PLACES=threads` and restrictive one-core binding.
The same executable passed with `OMP_PROC_BIND=false`, and the release build
passed that restrictive launch after host OpenMP was disabled. CUDA now defaults
host OpenMP to off, matching the working Make workflow. Native OpenMP and the
HIP Make threading default are retained.

See [build and on-demand test instructions](unified-backends.md) for reproduction
and the Frontier HIP handoff.
