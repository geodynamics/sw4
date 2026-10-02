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

## Receiver read-back and restoration (R2–R5)

SAC allocation and sample basis are now separate decisions. Geographic histories
are always inverted to internal grid components, including restart; Cartesian
histories retain their grid samples. All quantities use one shared component and
filename mapping for output and restart, including velocity, curl, divergence,
strains and displacement gradients. Missing/truncated/nonfinite/incompatible
histories fail before their data are used. USGS ingestion now handles the
quantity's complete column set instead of assuming three displacement columns.

New vector SAC output records `KUSER0=SW4XYZ` or `SW4ENU`. Cartesian files also
store the grid-to-geographic horizontal matrix in `USER0`–`USER3`. Matching-grid
observations retain their original samples; another grid uses the stored matrix
and the receiving grid's inverse. Unmarked external observations continue to use
`CMPAZ/CMPINC` for instrument rotation. For legacy unmarked SW4 Cartesian SAC
observations, specify `observation ... sacbasis=grid`; the default is `auto`.
Explicit grid input requires matching X/Y/Z (or velocity) component names and
rejects an incompatible geographic marker. Legacy restart uses its configured
output basis and validated component names.

Receiver HDF5 selects the requested quantity's exact datasets and units,
validates all component extents/counts/finiteness and downsampling, and keeps
restart allocation/timing intact. Its metadata readers validate scalar/string
sizes; waveform reads use an explicit bounded memory dataspace. `ignore_utc`
now deliberately preserves the simulation's reference UTC. New files record
relative `STARTTIME` and dimensionless quantities use `UNIT=1`; legacy derivative
files without UNIT remain readable through their exact component names.
The remaining GPU UTC-difference millisecond conversion now uses SW4's internal
microsecond scale, consistent with the previously approved timestamp fix.

```sh
cmake --build build --target sw4_receiver_check
python tests/backends/receiver_formats.py --sw4 build/bin/sw4 \
  --reader build/bin/sw4_receiver_check --work-dir "$SCRATCH/sw4-formats"
python tests/backends/receiver_restart.py --sw4 build/bin/sw4 \
  --work-dir "$SCRATCH/sw4-restarts"
python tests/backends/receiver_negative.py --reader build/bin/sw4_receiver_check \
  --fixture "$SCRATCH/sw4-formats/displacement-nsew0-ds1-topo0-utc234567" \
  --work-dir "$SCRATCH/sw4-invalid-receivers"
```

CPU forward-format checks passed 22 cases. A full checkpoint matrix passed eight
comparisons: two successive restorations per elastic azimuth-0/27, attenuating
azimuth-27 and refined azimuth-27 case. Every receiver quantity was exercised
through SAC, USGS and HDF5 with downsample 1/3. Receiver comparisons use
`1e-12 + 2e-6 * component_peak`, independently for each component; final PDE state
uses `1e-12 + 2e-5 * component_peak`. Observed PDE state differences were zero;
receiver differences were at float32 storage precision. Nineteen real-reader
failure/basis cases passed, including missing/truncated files, malformed scalar
metadata, wrong quantities/units/basis/downsampling, nonfinite samples,
big-endian SAC and externally oriented orthogonal channels at 13/45 degrees.

These are developmental-build results under
`/pscratch/sd/h/houhun/sw4-review-r1-r6`; immutable-source CPU/CUDA reruns and
refreshed performance acceptance remain necessary. Restart validation keeps the
planned simulation duration unchanged. Existing fixed-size HDF5 receiver output
cannot yet be extended from a shorter originally planned duration.
