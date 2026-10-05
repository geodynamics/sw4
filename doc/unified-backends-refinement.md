# Refinement validation and MPI affinity

The eight-rank, one-thread Cartesian energy timeout was reproduced outside SW4
with an MPI/OpenMP program containing only halo exchanges and reductions.
On Perlmutter, `srun -n8 -c1 --cpu-bind=cores` combined with inherited
`OMP_PLACES=threads OMP_PROC_BIND=spread` placed pairs of ranks on the same logical
CPU after OpenMP initialized. Their stacks stopped in Cray MPI's shared-memory
collective barrier. Changing MPI thread initialization or disabling optimized
Allreduce did not resolve this launch configuration.

With `--cpu-bind=threads`, the isolated program completed 10,000 checked
halo/reduction iterations. Allocating two logical CPUs per rank or disabling
OpenMP binding also completed the reproducer. The unchanged production CPU solver
then completed the original 535-step `energy-mr-4th-1.in` test with eight ranks and
one thread; relative energy drift after initialization was `7.39282463e-15`.
These results identify the launch affinity as the cause of this timeout and do
not justify changing the Cartesian interface equations or adding MPI barriers.

The affected on-demand drivers now explicitly select `--cpu-bind=threads` for
their one-thread ranks. The existing `pytest/` and `pytest-sw4mopt/` files are
unchanged. When launching their CPU solver tests on Perlmutter, replace their
legacy `--cpu_bind=cores -c1` launch binding with `--cpu-bind=threads`, or allocate
a full physical core per rank. Check actual affinity after OpenMP initializes.
See [NERSC's affinity guide](https://docs.nersc.gov/jobs/affinity/).

An additional regression is available explicitly in a single-node CPU allocation:

```sh
cmake --build <cpu-build> --target sw4_mpi_check
/global/cfs/cdirs/m3354/perl/conda_eqsim/bin/python \
  tests/backends/single_thread_refinement.py \
  --sw4 <cpu-build>/bin/sw4 --probe <cpu-build>/bin/sw4_mpi_check \
  --tasks 8 --work-dir <scratch-directory>
```

The probe checks distinct logical CPUs after an OpenMP region and verifies every
Allreduce result across 10,000 halo exchange sequences. The solver regression
requires a completed 535-step history and relative energy drift at most `1e-10`.
The probe target is excluded from ordinary builds; no tests are registered with
CTest or the existing pytest suites.

Diagnosis is retained under
`/pscratch/sd/h/houhun/sw4-interface-fixes-20261005`, using CPU allocation
`59379208` and the archived `c28cc06d` production core. Failed launches, affinity
snapshots, diagnostic stacks and private instrumented builds are retained
separately from accepted runs. The private curvilinear iteration experiments in
that directory are diagnostic evidence, not a released production build.

## Fresh validation of the affinity correction

Commit `48407f2ac2daf5f3d41d9c27033a8cdb995d6e3b` was archived and built cleanly
for standalone OpenMP and CUDA. It changes validation launchers and adds an
explicitly built MPI probe; production solver sources are unchanged from
`cf06e954`. The [validation record](validation/perlmutter-refinement-2026-10-05.json)
contains build, source archive, executable, driver and evidence identities.

| Check | OpenMP | CUDA |
| --- | --- | --- |
| Existing solver suite, eight ranks and one thread | 20 passed, 5 disabled optional cases | Not repeated; production solver unchanged |
| Independent MPI affinity and 10,000 reductions | Passed | CPU-only runtime probe |
| Original 535-step Cartesian refinement energy | Drift `7.39e-15` | Not repeated |
| Halo sentinels, square/X-strip/Y-strip | 3 shapes, 8 ranks | 3 shapes, 4 nodes, both array layouts |
| Refined receiver/PDE restart | 2 restorations, 8 ranks | 2 restorations, 2 GPUs |
| Short/misaligned/sufficient receiver history | 21 cases, 8 ranks | 21 cases, 2 GPUs |
| Mixed-output receiver/PDE restart | 16 cases, 8 ranks | 16 cases, 2 GPUs |
| Independent retained mixed-reference replay | 5 restored cases | 5 restored cases |
| Malformed receiver input and missing station | 23 native cases | Not repeated |

Mixed tests cover geographic displacement with decimation 3 at checkpoint 24,
Cartesian velocity with decimation 1 at checkpoint 25, and the associated
negative/tail variants. Restarts compare every component/dataset using its own
peak, with receiver tolerance `1e-12 + 2e-6 * peak` and PDE tolerance
`1e-12 + 2e-5 * peak`. No floating-point byte identity is required. The CUDA
allocation was `59381003`, containing three nodes of A100-SXM4-40GB devices and
one node of A100-SXM4-80GB devices. The halo oracle uses one GPU per node.
The `48407f2a` milestone changes no production hot path; hmr3 performance was not
remeasured for these test-launch changes.

## Curvilinear convergence correction

Before `a4504307`, the physical-source Gaussian-topography case exhausted the
production solver's 31-iteration curvilinear block Jacobi allowance and copied
an unconverged interface into the simulation. Cartesian `cmaxit`
does not configure this loop. A private diagnostic build retaining the equations,
relaxation and `1e-6` tolerances converged all 832 interface solves of the retained
0.8-second constant-material case, using at most 373 iterations. Its first solve
requiring more than 30 iterations is near simulation time 0.098 seconds.
One-, four- and eight-rank runs also converged and their receiver comparisons
passed the backend tolerance. Ten gradient-material cases, covering attenuation
on/off and extrapolation widths 1–5, needed at most 365 iterations.

An independent basis-vector diagnostic checked the actual interface operator
against its assembled diagonal blocks at interior and boundary-adjacent points;
the largest relative discrepancy was `6.26e-9`, below the solve tolerance. These
checks support an inadequate iteration cap rather than changing the equations.

The prepared patch raises the default cap to 500, adds a checked positive-integer
`developer curvimaxit=...` setting, enforces the exact limit, and aborts instead
of copying an unconverged interface. Verbosity 3 reports the iteration count and
relative/scaled residuals. A private clean CPU build of that proposal completed
all 832 solves, stopped after exactly 1/30 iterations when those limits were
requested, and rejected zero, negative, fractional, nonnumeric and overflowing
limits during parsing. The patch and private results are retained in the evidence
root as `curvi-convergence-proposal.patch` and `proposal-*-results.json`.

The user approved the correction after review `3296c3cd`. Commit `a4504307`
applies it to both OpenMP and RAJA backends, preserving the operator, relaxation
and `1e-6` relative/absolute tolerances. `developer curvimaxit=<positive integer>`
sets an exact maximum number of Jacobi updates; the default is 500. This is a
limit, not a convergence guarantee. If neither existing residual criterion is
satisfied at the limit, the run fails collectively before copying interface
planes or attenuation state into the evolving solution. A nonfinite residual
also causes collective failure. Local NaNs are mapped to infinity before device
and MPI maximum reductions so another finite residual cannot hide them.

Use `fileio verbose=3` to retain each solve's grid, time, update count and residual
ratios. Commit `9e3bb035` preserves 17 significant digits in these diagnostics
even when GPU initialization changes global output formatting. Equations,
scientific precision and stopping tolerances are unchanged. Runs that previously
continued after nonconvergence can now converge with additional iterations or
terminate with an error; partial receiver output from a failed run is not an
accepted simulation.

The explicitly invoked driver `tests/backends/curvilinear_refinement.py` covers
physical interface arrivals, gradient materials with attenuation on/off and
extrapolation widths 1–5, several resolutions, multiple curved interfaces,
longer propagation, curved-grid energy, exact small-limit failures, invalid
limits and localized nonfinite source input. `--compare-to <reference run>`
compares every receiver component with `1e-12 + 2e-5 * reference_peak`; successful
solves must satisfy the original residual criteria. It is not registered with
pytest or CTest. Existing suites remain unchanged.

The energy check uses the existing SW4 bound of `1e-10` relative change per
timestep, strengthened to check both signs, and records cumulative drift
separately. In the four-rank curved fixture, a pre-fix `48407f2a` run reproduces
the corrected CPU energy history exactly: maximum timestep change `7.41e-11`,
cumulative drift `7.80e-9` over 768 samples. The initial cumulative `1e-10`
validation bound was inappropriate for this unchanged interface tolerance;
the failed check and baseline comparison are retained. This is evidence about
that fixture, not a long-duration energy guarantee or permission to loosen the
solver's scientific tolerances. HIP execution remains deferred to Frontier.

The follow-up audit also found a separate inherited overflow in the residual
scale: `(global_nx-1)*(global_ny-1)` multiplies 32-bit integers before conversion
to `float_sw4`. The extracted production expression triggers signed-overflow
sanitization for a 50,001-by-50,001 interface. A negative wrapped scale can falsely
satisfy the absolute stopping criterion. Evaluating the same formula in floating
point avoids this overflow; that additional correction is prepared privately and
awaits separate scientific approval. Interfaces whose horizontal interval product
exceeds `INT_MAX` remain outside the qualified scope. The proposed patch and
arithmetic checks are retained as `scale-correction-proposal.*` in the new evidence
root `/pscratch/sd/h/houhun/sw4-curvi-review-20261005`.

## Validation of the approved correction

Clean CPU and CUDA builds pin production source `9e3bb035`; later commits change
only on-demand tests and documentation. The source/CMake file manifest verifies
that the tested production files match the branch. Qualification drivers pin
`6eb6f15e`, with the earlier driver identities retained for individual runs.

| Check | OpenMP | CUDA |
| --- | --- | --- |
| Physical-source constant/gradient material, attenuation, extrapolation, resolution and curved-level configurations | 15 passed | 15 passed against converged CPU references |
| Curved-grid energy, 768 samples | Timestep change `7.41e-11`; total drift `7.80e-9`, unchanged from pre-fix CPU | Passed timestep criterion |
| Manufactured elastic/attenuating solutions, three resolutions each | 6 passed; orders `3.93–4.03` | 6 passed; same error norms and convergence orders |
| Exact caps 1/30, five malformed limits, localized NaN | 8 expected failures | 8 expected failures |
| Physical source across decompositions | 1/4/8 ranks; multi-interface and heterogeneous attenuation also checked at 8 | 1/4/8 ranks; 8-rank physical, multi-interface and heterogeneous attenuation run on 4 nodes |
| Existing enabled solver tests | 20 passed; 5 optional cases disabled | 20 passed; 5 optional cases disabled |
| Production halo sentinels, square/X-strip/Y-strip | 3 shapes, 8 ranks | 3 shapes, 4 nodes, both layouts |

The physical-source case converges all 832 solves with at most 373 updates.
The two-curved-interface case converges 3,320 solves with at most 380 updates;
the 1.6-second case converges 1,660 solves. Every accepted solve's finite residuals
meet at least one of the unchanged relative/scaled absolute criteria. Component
waveforms require meaningful late arrivals and backend error at most
`1e-12 + 2e-5 * reference_peak`. The manufactured solutions also reproduce the
repository's independent exact-solution error references; halving spacing gives
approximately sixteenfold error reduction, including attenuation state.

Fresh `performance/large/hmr3.in` runs alternate three RAJA/unified pairs on the
same four nodes and sixteen GPUs, after compilation completes. RAJA solver times
are 82.58, 82.64 and 82.18 seconds; unified times are 79.43, 79.09 and 79.06 seconds.
The medians, 82.58 and 79.09 seconds, give 4.23% lower solver time for unified.
All three waveform comparisons pass. The allocation contains two nodes of 40 GB
A100s and two nodes of 80 GB A100s. This is a scoped workload result, not a general
speedup or proof that previously unconverged cases cost no additional iterations.

The [qualification record](validation/perlmutter-curvilinear-2026-10-05.json)
retains builds, caches, driver/source archives, every input and residual log,
receiver and energy artifacts, replay results, hardware and performance
identities. Failed preparation trials and the stopped low-precision CUDA trial
are excluded and retained separately. Additional checks remain on demand and
existing pytest suites are unchanged.

| Commit | Change | Retained validation |
| --- | --- | --- |
| `9aa87956` | Actual HDF5 ISO epoch mutation and explicit-epoch rejection | 62 CPU/CUDA cases and 24 restored-output replays |
| `a4504307` | Approved finite/converged interface gate and configurable exact cap | Physical arrivals, exact-limit failures, NaN rejection; final build below |
| `4c09211b` | Epoch evidence and initial curved regression driver | Epoch qualification; later driver corrections below |
| `9e3bb035` | Full-precision convergence diagnostics | Clean CPU/CUDA builds, 22 configurations each, rank/halo/energy/performance replay |
| `9d5f18fa` | Valid multi-interface layers and exact residual gates | Two interfaces, 3,320 converged solves on CPU/CUDA and multiple decompositions |
| `a64c40ce` | Established energy timestep criterion, both signs | Identical pre-fix CPU history and passing CPU/CUDA energy checks |
| `6eb6f15e` | Three-resolution manufactured-solution regression | Elastic/displacement/attenuation error orders `3.93–4.03` on both backends |
