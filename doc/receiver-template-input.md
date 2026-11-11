# Receiver-template input

Rank zero of each simulation MPI group reads the HDF5 station metadata and
broadcasts ordered batches to its group. The input schema is unchanged.
Receiver construction, sampling, output files and owner-rank waveform reads
follow the existing paths. Missing or invalid station metadata stops setup
before time stepping. The setup log separates read, distribution and receiver
construction times.

`USEZVALUE` selects an absolute z coordinate. SW4 checks the receiver against
the local terrain without adding the surface z to the supplied coordinate.
Scalar `WindowL` and `WindowR` override the command bounds for their station.
Four-value `WINDOWS` takes precedence and uses the maximum of the two
five-sample tanh tapers where they overlap. The second interval is active
when at least one of its bounds is positive.

The normalized L2 time-shift derivatives include the shift dependence of the
observation energy. A travel-time measurement requires a resolved interior
correlation peak. Unresolved windows contribute zero score and zero adjoint
forcing. SW4 logs the reason and sampled lag. Native CPU inversion fixes
the accepted set at the first evaluation and rejects trial models that lose
an accepted window. Restart with reviewed windows if no trial is admissible.

Run the receiver checks with the repository's scientific Python. Each mode
lists its inputs in `--help`.

| Test entry point | Modes |
| --- | --- |
| `pytest/receiver_template_metadata.py` | `replay` checks station records and malformed inputs locally. `groups` checks two independent MPI groups. |
| `pytest/receiver_template_startup.py` | `startup` compares matched HDF5-enabled builds. `depth` checks receiver placement; `geometry` checks source z and depth-gradient materials on flat/Gaussian grids. |
| `pytest/receiver_objective_numerics.py` | `windows` and `numerics` check receiver methods locally. `native-build` links a driver against an existing SW4 build. |

Use a native CPU baseline for metadata parity. Run `groups`, `startup`,
`depth`, and the native objective executable inside an authorized allocation.
For native objectives, use a flat uniform-material deck with (1800, 1800, 0)
inside the grid. The Taylor mode checks the material gradient with the
native CPU sw4 and sw4mopt executables in an allocation.

Geometry checks compare native rho/Vp/Vs planes with analytic depth gradients,
match relative and absolute source positions, and verify above-surface source
exclusion on one Gaussian hill and vertical material plane. The Taylor check
covers one Vp direction in a small uniform model.
GPU model-level adjoints, representative-recording peak selection, HIP and
representative-scale startup performance remain unverified.
