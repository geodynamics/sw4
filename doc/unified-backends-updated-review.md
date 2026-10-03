# Updated production review corrections

This is the historical validation record for `788f96a0`. The subsequent review
identified observation workflow and restart-validation defects, plus a coverage
overstatement: its 32 named propagation comparisons represent 30 unique physical
setups and omit exact-interface force/moment propagation. See
[the new review corrections](unified-backends-reader-workflows.md) for the fixes
and replacement qualification record. The old numerical comparisons remain
retained evidence; their case names alone do not establish scientific coverage.

All six findings in `review-unified-backends-dec83d29/REPORT.md` were confirmed
and their authorized root causes corrected. The existing `pytest/` and
`pytest-sw4mopt/` suites remain unchanged. Additional regression drivers and probe
targets run explicitly on demand; they are not added to the default test run.

The final compiled solver is `bf983ac72bac1f45ebeee328cc415287f936d17f`.
The later `53ff7c7c` commit corrects test case naming only. An immutable Git archive,
clean CPU and CUDA builds, configuration scripts, binary hashes and run logs are
retained under `/pscratch/sd/h/houhun/sw4-review-r1-r6/release-bf983ac7`.
The archive SHA-256 is
`4779ab24f61c63117379affa4b9aa8c6b981c23226dbef8fc9641bee99ba6623`.
The [machine-readable validation record](validation/perlmutter-updated-review-2026-10-02.json)
identifies each build and check, including earlier evidence used for unchanged
source/material operators. The result qualifies the configurations and fixtures
below; it does not establish unrestricted production acceptance.

## Root causes and corrections

| Finding | Correction |
| --- | --- |
| R1: GPU anisotropic setup rejected 21-component halos | Complete 21-component dispatch, packing descriptors and buffer sizing for component-major and interleaved arrays, preserving producer/MPI/unpack completion. |
| R2: geographic SAC restart changed the component basis | Restore geographic histories through the inverse output transformation; keep internal recording arrays in grid components and separate basis conversion from restart allocation. |
| R3: GPU restart selected displacement filenames for every mode | Share quantity/component/filename definitions between output and readers, including velocity, curl, divergence, strains and displacement gradients. |
| R4: receiver HDF5 could not read its own velocity output | Select the requested quantity's exact datasets and units, validate all extents and metadata, and preserve restart timing/allocation. |
| R5: Cartesian SAC observation read-back rotated samples incorrectly | Mark new files with their basis and store the horizontal transformation; preserve same-grid samples while retaining external instrument rotation and explicit legacy grid input. |
| R6: propagation fixtures missed special interface source stencils | Derive positions from actual grid dimensions/spacing, independently verify forces and all nine first moments, and retain distinct exact-interface/epsilon propagation cases. |

### GPU transport and source placement

CMake now exposes a reusable `sw4_gpu` core so the GPU probes use the production
transport and initialization. The halo sentinel starts a pending GPU producer,
fills rank/component/index-specific values, checks every array entry after
exchange, and repeats with new values. It exercises 1/3/4/21 components, physical
boundaries, corners, two Cartesian grids, uneven extents and square/X-strip/Y-strip
MPI decompositions. Native checks use component-major layout; CUDA checks use both
layouts. Four-rank native and four-node/four-rank CUDA checks passed all three
shapes with staged buffers, Umpire and `MPICH_GPU_SUPPORT_ENABLED=0`.

The source oracle covers all six special interface positions, far controls, the
exact interface and both epsilon offsets, using SBP volume weights and sources
spanning MPI ownership. Its independent expectations cover forces and all nine
first moments. All 22 cases passed on CPU and four-node CUDA.

The propagation driver verifies its positions against each executed solver's
printed grids. Ten-digit formatting collapsed exact and epsilon-sided depths into
one filename; the old `z2000` decks actually used positive epsilon. The first
sweep completed 28 unique setups on `47390e78`. Commit `53ff7c7c` retained full
depth precision, but the four subsequent named runs on `bf983ac7` added only two
new physical setups: two positive-epsilon runs duplicated earlier inputs. Thus
this historical set contains 32 named comparisons, 30 unique physical setups,
and no exact-interface force/moment propagation. The new review record replaces
this acceptance inventory using input-derived scientific identities.

### Receiver basis, quantities and checked ingestion

New vector SAC output records `KUSER0=SW4XYZ` or `SW4ENU`. Cartesian files also
store the grid-to-geographic horizontal matrix in `USER0`–`USER3`. Matching-grid
observations retain their samples; another grid uses the stored matrix and the
receiving grid's inverse. Unmarked external observations continue to use
`CMPAZ/CMPINC` for instrument rotation. For legacy unmarked SW4 Cartesian SAC
observations, specify `observation ... sacbasis=grid`; the default is `auto`.
Explicit grid input requires matching component names and rejects an incompatible
geographic marker. Legacy restart uses its configured basis and validated names.

SAC reads check payload bounds before allocation, endian conversion, time-series
shape, counts, timing and finite values. HDF5 validates scalar/string metadata,
component extents, quantities/units, downsampling, UTC and samples, and uses bounded
memory dataspaces. New HDF5 output records relative `STARTTIME` and dimensionless
`UNIT=1`; legacy derivative files without UNIT remain readable through their exact
component names. `ignore_utc` deliberately preserves the simulation reference UTC.
USGS ingestion validates the complete quantity-specific column set, units, uniform
times and finite samples, and restores geographic vectors to grid components.

Additional issues found while exercising these paths were corrected:

- GPU UTC differences and input-date parsing used milliseconds where SW4 stores
  microseconds. Both now use the CPU convention. Fractional seconds are parsed as
  double, avoiding float rounding of `59.234567` and `59.999999`; tests compare
  output UTC with the original input, not only with other formats.
- Missing station HDF5 files previously printed an error and silently omitted
  receivers. Input now fails before time stepping and closes its file-access
  property list; CPU and CUDA failure paths were verified.
- Corrupt HDF5 `NPTS` could trigger history allocation before component extents
  were checked. Validation now precedes allocation. A billion-sample declaration
  with a short dataset reproduced `std::bad_alloc` under a 2 GiB address-space
  limit before the fix; the final reader reports a checked failure under the same
  limit. Negative-reader fixtures now include valid station initialization so
  they actually exercise the intended reader failure.
- The checked reader retains the existing downsample-mismatch diagnostic, allowing
  the unchanged HDF5 restart regression to verify the failure contract.

## Final validation

Fresh builds used GNU 13.2.1, Cray MPI 9.1.0, parallel HDF5 1.14.3.9, CUDA 13.2
with architecture 80, RAJA/Umpire 2026.07 and double precision. The native inversion
executable was run from both the standalone CPU and combined CUDA build trees.
Build scripts and cache/binary hashes are recorded with the evidence.

| Check | CPU | CUDA / combined build |
| --- | --- | --- |
| Existing solver suite | 20 passed, 5 optional skips | 20 passed, same 5 skips |
| Existing native inversion, 32 MPI ranks × 2 OpenMP threads | 4 passed | 4 passed from CUDA tree's native executable |
| SAC/USGS/rechdf5 output plus actual reader, fractional UTC | 24 cases passed | 24 cases passed |
| Negative reader, basis and missing-station checks | 23 passed | 22 reader cases plus separate missing-station failure passed |
| Elastic/attenuating/refined receiver checkpoint matrix | 8 comparisons passed | 8 comparisons passed |
| Fully coupled anisotropic flat/topographic checkpoint matrix, 4 ranks | 4 comparisons passed | 4 comparisons passed |
| Existing HDF5 regression helpers | 4 passed | 4 passed |
| Halo sentinel, square and both strip decompositions | 3 shapes passed | 3 shapes passed on 4 nodes, both layouts |
| Independent interface force/moment oracle | 22 cases passed | 22 cases passed on 4 nodes |
| General CPU/CUDA waveforms, including text/HDF5 SRFs | 7 comparisons passed | Same comparisons |
| Source/material/event propagation, with completed epsilon coverage | 32 named comparisons passed (30 unique physical setups) | Same comparisons; build split described above |
| Final anisotropic isotropic-limit/reference checks | Passed | Passed |

The native source-history oracle passed six filtered/unfiltered scalar, force and
moment combinations with maximum error `2.66454e-15`. The independent native
anisotropic polynomial oracle covered all Cartesian interior/top/bottom rows,
isotropic and coupled tensors, with maximum error `8.41283e-12`. The final coupled
anisotropic CPU/CUDA comparisons differed by at most `8.14e-15` in waveforms and
`8.57e-14` in checkpoint fields across the flat/topographic fixtures.

Waveform acceptance uses `1e-12 + 2e-5 * reference_peak` separately for each
component. Float32 SAC/HDF5 storage comparisons use `1e-12 + 2e-6 * component_peak`;
restart PDE comparisons use `1e-12 + 2e-5 * dataset_peak`. Sample axes, component
names, units, counts, orientation and finite values are also checked. Same-backend
restart PDE differences were zero in the retained fixtures; receiver differences
were at float32 storage precision. Byte identity is not required for floating
point acceptance.

### Four-node performance acceptance

The exact `raja:performance/large/hmr3.in` input and retained RAJA baseline binary
were compared against the final solver on the same four Perlmutter nodes, using
16 A100-SXM4-40GB GPUs, staged MPI and three alternating baseline/candidate pairs.
Input, executable and model hashes, launch commands and hardware identity are in
the validation record. All three waveform comparisons passed the component-wise
threshold. The median solver time was **82.01 s for RAJA and 79.26 s for the final
solver**, a **3.35% reduction**, passing the 5% maximum-slowdown gate.
These timings qualify this workload and allocation, not every SW4 workload.

## Running additional checks on demand

```sh
cmake --build build --target sw4_halo_check sw4_interface_check sw4_receiver_check
python tests/backends/halo_transport.py --probe build/bin/sw4_halo_check \
  --backend CUDA --nodes 4 --tasks 4 --work-dir "$SCRATCH/sw4-halos"
srun -n 4 build/bin/sw4_interface_check refined-case.in
python tests/backends/receiver_formats.py --sw4 build/bin/sw4 \
  --reader build/bin/sw4_receiver_check --work-dir "$SCRATCH/sw4-formats"
python tests/backends/receiver_restart.py --sw4 build/bin/sw4 \
  --work-dir "$SCRATCH/sw4-restarts"
python tests/backends/receiver_restart.py --sw4 build/bin/sw4 --anisotropic \
  --work-dir "$SCRATCH/sw4-anisotropic-restarts"
python tests/backends/receiver_negative.py --reader build/bin/sw4_receiver_check \
  --fixture "$SCRATCH/sw4-formats/displacement-nsew0-ds1-topo0-utc234567" \
  --work-dir "$SCRATCH/sw4-invalid-receivers"
```

Run MPI/scientific workloads inside an allocation. Use the retained launch scripts
for the exact validated rank/node layouts. Probe targets are `EXCLUDE_FROM_ALL`
and are not registered in CTest. For an immutable archive build, set
`-DSW4_SOURCE_REVISION=<full-40-digit-commit>`; retain the source archive and hashes
because the executable version string alone does not prove build provenance.

## Remaining qualification limits

- HIP compilation, device linking, correctness and performance require Frontier
  validation; Perlmutter provides no positive HIP evidence.
- Physical-source curvilinear refinement fixtures still fail the inherited CPU
  interface convergence cap (31 iterations). Those failed fixtures are excluded;
  the numerical iteration policy was not changed. Passed existing curvilinear and
  topography checks do not resolve that separate convergence issue.
- Displacement and velocity cover grid/geographic output; divergence, curl,
  strains and gradients cover grid output. Geographic curl USGS output is an
  inherited unsupported case awaiting approval for correction. Geographic strain
  and gradient text output are also unsupported.
- The observation parser still selects displacement. Forward velocity read-back
  and restart checks do not qualify velocity adjoint/inversion.
- Restart tests preserve the originally planned duration. Fixed-size HDF5 receiver
  output cannot yet extend a shorter originally planned run.
- Anisotropic refinement and attenuation remain unsupported by the existing
  solver. Variable-coefficient anisotropic convergence, nonunit supergrid stretch,
  long-time stability and broader physical validation require further checks.
- Managed-buffer MPI, allocator/transport variants, optional SZ, broader scaling
  and long production duration remain separate acceptance gates. No GPU memory
  sanitizer or exhaustive supported-combination qualification is claimed.
