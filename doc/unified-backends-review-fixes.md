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
