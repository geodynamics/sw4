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
The current milestone changes no production hot path; hmr3 performance was not
remeasured for these test-launch changes.

## Curvilinear convergence proposal awaiting approval

The physical-source Gaussian-topography case still fails in the production
solver: its independent curvilinear block Jacobi loop stops after 31 iterations
and copies an unconverged interface into the simulation. Cartesian `cmaxit`
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

This solver change awaits the user's required scientific-correctness approval.
It can change runs that previously continued unconverged. It has not been applied
to the production branch, committed as solver code, or qualified on CUDA/HIP.
CUDA numerical and hmr3 performance checks remain required after approval.
