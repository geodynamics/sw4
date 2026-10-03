# Observation and restart workflow corrections

The findings N1–N4 in `review-unified-backends-788f96a0/REPORT.md` were checked
against current code and actual retained inputs. The user authorized fixes and
on-demand regressions. Existing pytest suites remain unchanged.

The solver corrections are committed in `62f1ca12`, `5c236501` and `95443f0e`.
Final qualification uses clean CPU/CUDA builds from an immutable archive of
`95443f0e`, under `/pscratch/sd/h/houhun/sw4-review-n1-n4/release-95443f0e`.
Later commits only add on-demand tests and documentation.

- N1: a shared observation path resolver maps local event indices to global
  observation directories and preserves absolute input paths. Absolute input
  paths use a local diagnostic basename after ingestion, preventing invalid
  output paths. Short SAC filenames no longer underflow the suffix check.
- N2: USGS history must reach the checkpoint before its last-valid index is reset.
  SAC, USGS and HDF5 also validate the reference UTC plus relative start against
  the allocated simulation axis, with storage-aware timestamp tolerances. Absolute
  UTC matching applies when `time utcstart=...` supplies an explicit epoch.
  Without it, the default wall clock changes between invocations and cannot be
  used as an absolute restart identity; relative-axis validation still applies.
- N3: station discovery uses the checked SAC decoder before placement, handles
  either endian, and validates inferred coordinates. Explicit coordinate overrides
  remain available with undefined header locations. Truncated SAC headers are
  rejected before accessing uninitialized integer fields.
- N4: propagation coverage comes from executed source coordinates and actual grid
  dimensions. The audit rejects labels that disagree with inputs and duplicated
  scientific setups, and requires all six special stencils plus exact and both
  epsilon-sided force/moment positions.

Pre-fix real-solver tests reproduced wrong-event ingestion, acceptance of short
USGS history and omitted big-endian observations. The existing HDF5 restart test
also caught the initially overstrict absolute-UTC check for inputs without an
explicit epoch; `95443f0e` corrects that regression. Additional tests deliberately
replace stored default epochs with the year 2000, require successful restoration,
and still reject a shifted relative start.

The approved N1–N4 validation gates pass for the configurations below.
An additional mixed-output restart defect remains open pending scientific approval. The machine-readable
[validation record](validation/perlmutter-reader-workflows-2026-10-02.json) retains
source/tree/archive, executable, build-cache, driver, reference and restored-output
identities. CPU build/validation used allocation `59245457`; four-node CUDA used
`59247613`; the combined tree's native inversion suite used CPU allocation
`59249195`. Earlier failed staging attempts and the overstrict UTC intermediate
build are retained separately and excluded from acceptance.

| Check | Standalone CPU tree | Combined CUDA tree |
| --- | --- | --- |
| Existing solver suite, level 0 | 20 passed, 5 optional cases disabled | 20 passed, 5 optional cases disabled |
| Existing native inversion, 32 MPI ranks x 2 OpenMP threads | 4 passed | 4 passed, native solver |
| Real observation parser, serial/parallel paths, objective and scaled gradient | 9 groups passed | 9 groups passed, native solver |
| SAC/USGS/rechdf5 samples and metadata | 24 cases passed | 24 cases passed |
| Short, sufficient, extra-tail and misaligned receiver histories | 21 cases passed | 21 cases passed |
| Default-epoch restoration and shifted-relative-origin rejection | 8 cases passed | 8 cases passed |
| Elastic/attenuating/refined restart comparisons | 8 passed | 8 passed |
| Coupled anisotropic flat/topographic restart comparisons | 4 passed | 4 passed |
| Existing HDF5 helpers, including SSI/ZFP lifetime | 4 passed | 4 passed |
| Halo sentinel, both strip decompositions and square | 3 shapes passed | 3 shapes passed on 4 nodes, both layouts |
| Independent interface force/nine-moment oracle | 22 cases passed | 22 cases passed on 4 nodes |
| Source/material/event propagation | 32 fresh comparisons; 22 unique interface setups | Same comparisons, GPU on 4 nodes |
| General propagation including Sfile and text/HDF5 SRFs | 7 comparisons passed | Same comparisons |

The native negative-reader matrix passed 23 cases. The source preparation oracle
passed six scalar/force/moment filtered/unfiltered cases, including copies and
repeated preparation; maximum error was `2.66454e-15`. The independent anisotropic
polynomial oracle passed all interior and both closure rows, with maximum error
`8.41283e-12`. Both backends also passed the propagated anisotropic isotropic limit.
The coverage audit rejects the old mislabeled inputs and a deliberately repeated
input collection. The independent restart replay passed 36 retained final cases
and recorded SHA-256 identities of both references and restored artifacts.

The existing solver suite's five skips are the explicitly disabled legacy
HDF5/Geodyn cases (`-u0 -g0`), not numerical failures. Separate HDF5, receiver,
restart, Sfile and SRF drivers cover the changed paths. Existing pytest files are
unchanged from `788f96a0`.

| Commit | Milestone | Final validation gate |
| --- | --- | --- |
| `62f1ca12` | Event paths, checked SAC placement, restart history/time-axis validation | Real inversion observations and negative/positive restarts |
| `5c236501` | Absolute observation diagnostics and short SAC names | Relative/absolute paths and short-name parser cases |
| `8f3f4fb8` | On-demand observation/restart/input-coverage regressions | Both native inversion trees and CPU/CUDA forward runs |
| `95443f0e` | Explicit versus wall-clock-default UTC handling | Existing HDF5 restart and default-epoch cases |
| `564d3378` | Deterministic implicit-UTC regressions | All four receiver storage modes |
| `9df9de3e` | Independent retained-reference replay | Final restarted waveforms and PDE checkpoints |
| `587906d4` | Restored artifact identities | SHA-256 of the independently audited outputs |


The final four-node, sixteen-GPU `performance/large/hmr3.in` comparison ran three
alternating baseline/candidate pairs. Baseline solver times were 81.79, 81.43 and
81.85 seconds; candidate times were 78.79, 78.66 and 78.51 seconds. The medians
are 81.79 versus 78.66 seconds, a 3.83% improvement, within the 5% allowed slowdown
limit. All three waveform comparisons passed. This allocation had one node of
A100-SXM4-40GB devices and three nodes of A100-SXM4-80GB devices; both executables
ran on the same allocation with four MPI ranks per node, one GPU per rank,
`OMP_NUM_THREADS=1`, staged buffers and `MPICH_GPU_SUPPORT_ENABLED=0`. Per-run hardware and simulation identities are retained with these paired measurements.

The on-demand drivers live in `tests/backends/` and are not registered with the
existing pytest suites or CTest. In an allocated CPU job, use the EQSIM Python to
run `receiver_observations.py --sw4 <CPU sw4> --mopt <sw4mopt> --work-dir <scratch>`
and `receiver_restart_negative.py --sw4 <sw4> --work-dir <scratch>`.
Add `--implicit-utc` for default-epoch coverage and `--backend CUDA` for GPU
restart checks. CUDA-tree `sw4mopt` uses the native CPU solver; a GPU forward run
does not process inversion observation commands and cannot qualify that parser.
`receiver_restart_audit.py --evidence-dir <completed restart dir> --output <json>`
independently replays retained waveform and PDE NPZ references against the final
restored outputs. Earlier repeated restarts are checked during execution; their
intermediate rewritten outputs are not independently retained.

`review_regressions.py --gpu-nodes 4` runs CPU/GPU propagation comparisons;
`interface_coverage.py --evidence-dir <completed comparison dir> --output <json>`
requires 22 unique interface inputs and rejects duplicated or mislabeled setups.
`profile_hmr3.py` profiles the original raja input with alternating baseline and
candidate runs on four nodes and sixteen GPUs.

Waveform comparisons use each component's own peak: backend tolerance is
`1e-12 + 2e-5 * peak`, storage/restart receiver tolerance is
`1e-12 + 2e-6 * peak`. PDE checkpoint comparisons use each dataset's own peak
with `1e-12 + 2e-5 * peak`. Inversion gradients are scaled by density/stiffness
parameter scales and checked separately for each parameter type at the backend
tolerance. Floating-point byte identity is not required.


Qualification remains limited to the tested configurations and short regressions.
HIP requires Frontier compilation, device linking, execution and performance
checks. Physical-source curvilinear refinement still encounters the inherited
iteration cap in the excluded trial fixtures. Geographic curl USGS and geographic
tensor text output remain unsupported; velocity adjoint/inversion is not qualified.
Receiver HDF5 restart retains the originally planned storage duration and cannot
extend it. Parallel-event observation routing is tested, but multi-event
checkpoint naming/restoration is not qualified by these tests. Managed-buffer
MPI, optional SZ, heterogeneous material convergence, wider scaling and long
production runs remain separate acceptance gates.


An additional inherited receiver defect was reproduced after completing the main
matrix. With `sacformat=0 usgsformat=1 hdf5format=1 downSample=3`, `doRestart`
reads HDF5 and replaces the valid full-rate text prefix with cubic interpolation.
A complete HDF5 history changes text by up to 2.73% of a component peak. An exactly
sufficient HDF5 history leaves full-rate rows 22 and 23 unfilled at checkpoint 24,
with error up to 97.1% of a component peak. Both simulations finish and their PDE
checkpoint arrays match the uninterrupted baseline exactly. Thus PDE restart
parity does not establish receiver-output correctness for this combination.

Evidence is retained in `/pscratch/sd/h/houhun/sw4-review-n1-n4/mixed-restart-probe`,
using the archived `95443f0e` CPU executable and allocation `59250344`.
The proposed correction is to restore the available full-rate text history while
validating the requested HDF5 history and their common stored samples. This
changes delivered receiver waveforms and awaits the user's scientific-correctness
approval. The acceptance matrix above qualifies separate full-rate SAC/text and
HDF5 receivers; it does not qualify this mixed receiver configuration.
