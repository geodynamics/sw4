# Updated production review corrections

The six findings in `review-unified-backends-dec83d29/REPORT.md` were confirmed
against the current implementation. The user authorized their root-cause fixes.
This record separates completed checks from remaining production qualification.
The existing `pytest/` and `pytest-sw4mopt/` suites remain unchanged.

## GPU halos and interface source coverage (R1, R6)

GPU exchange now accepts stiffness arrays with 21 components. Both component-major
and interleaved layouts construct the matching packing descriptors and buffers;
the latter previously constructed only MPI datatypes. Buffer sizing includes the
21-component counts, and existing producer/MPI/unpack completion remains intact.
CMake provides `sw4_gpu` as a reusable core so explicit probes use the same
transport and initialization as the forward executable.

```sh
cmake --build build --target sw4_halo_check sw4_interface_check
python tests/backends/halo_transport.py --probe build/bin/sw4_halo_check \
  --backend CUDA --nodes 4 --tasks 4 --work-dir "$SCRATCH/sw4-halos"
srun -n 4 build/bin/sw4_interface_check refined-case.in
```

The sentinel fills rank/component/index-specific values with a pending GPU
producer, checks every array entry after exchange, and repeats with new values.
It covers 1/3/4/21 components, physical boundaries, corners, two Cartesian grids,
uneven extents and square/X-strip/Y-strip decompositions. Native checks use
component-major layout; GPU checks exercise both layouts. Four-rank CPU and
four-node/four-rank CUDA staged-buffer/Umpire checks passed all three shapes.
Evidence is under `/pscratch/sd/h/houhun/sw4-review-r1-r6/halo-cpu` and
`halo-cuda-four`. Managed-buffer transport and HIP remain separate gates.

The source oracle generates positions from actual grid spacing and dimensions,
checks all six special interface stencils, far controls and the exact interface
plus/minus epsilon, and integrates independent forces and all nine first moments
with the SBP volume weights. It also checks sources spanning MPI ownership.
The propagation driver generates the corresponding depth sweep and verifies its
mapping using each executed solver's printed grid dimensions. The earlier
1747/1847/... sweep did not cover half the special stencils.

An immutable source archive can be configured with
`-DSW4_SOURCE_REVISION=<full-40-digit-commit>` to retain its revision in executable
version metadata. Retain the archive/source, configure commands, dependency
versions and binary/script hashes alongside execution logs; the version string
alone does not establish provenance.
