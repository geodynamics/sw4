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
