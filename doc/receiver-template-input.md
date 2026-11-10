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
observation energy. When Newton refinement fails, SW4 logs the station and
uses the best sampled shift. Travel-time adjoint checks apply to isolated
correlation peaks. Review sampled fallbacks before using their sources in
inversion.

Run the receiver checks with the repository's scientific Python. Each mode
lists its inputs in `--help`.

| Test entry point | Modes |
| --- | --- |
| `pytest/receiver_template_metadata.py` | `replay` checks station records and malformed inputs locally. `groups` checks two independent MPI groups. |
| `pytest/receiver_template_startup.py` | `startup` compares matched HDF5-enabled builds. `depth` checks template and inline receivers over topography. |
| `pytest/receiver_objective_numerics.py` | `windows` and `numerics` check receiver methods locally. `native-build` links a driver against an existing SW4 build. |

Use a native CPU baseline for metadata parity. Run `groups`, `startup`,
`depth`, and the native objective executable inside an authorized allocation.
For native objectives, use a flat uniform-material deck with (1800, 1800, 0)
inside the grid.
