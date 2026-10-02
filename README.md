SW4 - Seismic Waves, 4th order accuracy
===========================================================
[![License GPL2+:](https://img.shields.io/badge/License-GPL%202%2B-red)](https://github.com/geodynamics/sw4/blob/master/LICENSE.txt)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8322590.svg)](https://doi.org/10.5281/zenodo.8322590)
[![pdf manual](https://img.shields.io/badge/get-PDF-green.svg)](https://github.com/geodynamics/sw4/blob/master/doc/SW4_UsersGuide.pdf)

More information
----------------
Please see the [Users Guide](https://github.com/geodynamics/sw4/blob/master/doc/SW4_UsersGuide.pdf) for more information regarding the use of SW4.

Additional resources:
- [SW4](https://geodynamics.org/resources/sw4/about)
- [Forum](https://community.geodynamics.org/c/sw4/32)

Cite Code As
------------
Petersson, N. Anders, Sjogreen, Bjorn, Tang, Houjun, & Pankajakshan, Ramesh. (2023, September 6). geodynamics/sw4: SW4, version 3.0. doi:10.5281/zenodo.8322590, url: https://doi.org/10.5281/zenodo.8322590

Primary References
------------------
Zhang, L.;  Wang, S.; and Petersson, N.A. (2021), Elastic Wave Propagation in Curvilinear Coordinates with Mesh Refinement Interfaces by a Fourth Order Finite Difference Method, SIAM J. Sci. Comp.  43(2) pp. A1472-A1496. doi:10.1137/20M1339702

Petersson, N.A.; Sjögreen, B. (2015), Wave propagation in anisotropic elastic materials and curvilinear coordinates using a summation-by-parts finite difference method, Journal of Computational Physics, 299, 820-841, doi: 10.1016/j.jcp.2015.07.023, url: http://linkinghub.elsevier.com/retrieve/pii/S0021999115004684

Petersson, N.A.; Sjögreen, B. (2012), Stable and efficient modeling of anelastic attenuation in seismic wave propagation, Communications in Computational Physics, 12 (01), 193-225. doi:10.4208/cicp.201010.090611a

Sjögreen, B.; Petersson, N.A. (2012), A Fourth Order Accurate Finite Difference Scheme for the Elastic Wave Equation in Second Order Formulation, Journal of Scientific Computing, 52 (1), 17-48, doi: 10.1007/s10915-011-9531-1, url: http://link.springer.com/10.1007/s10915-011-9531-1



User's Guide
------------
Petersson, N.A.; Sjögreen, B.; Tang, H. (2023), User's Guide to SW4, version 3.0, LLNL-SM-741439

Files
-----

- LICENCE.txt -    GNU General Public Licence version 2
- INSTALL.txt  - Information on how to build SW4
- README.txt  - This file!
- wave.txt    - Text file containing the "SW4 Lives" banner
- Makefile     - Main makefile

Directories
-----------
- configs/    -  Directory containing configuration files for "make"
- src/        -  C++ and Fortran source code for SW4
- tools/       - Matlab/Octave scripts for post processing and analysis of results
- examples/   - Sample SW4 input files
- optimize/    - Directory for object files and the optimized SW4 executable
- debug/       - Directory for object files and a SW4 executable with debug symbols


CPU HDF5 output
---------------

An event can write separate receiver HDF5 files with repeated station names.
Receiver datasets include step zero when downsampled; restarts overwrite samples
after the checkpoint and stop if the saved history cannot be restored.

SSI files expose the solver `timestep`, effective `outputtimestep`, and durable
`lastsw4timestep` / `lastoutputindex` progress markers. Checkpoints flush partial
SSI buffers and receiver output before saving the restart state.

CPU allocation regressions are provided in
`pytest/receiver_hdf5_metadata.py`, `pytest/restart_hdf5_output_offsets.py`,
`pytest/receiver_hdf5_restart_failure.py`, and
`pytest/ssi_zfp_dataset_lifetime.py`. Run them with an HDF5-enabled CPU executable;
compressed cases require ZFP. Set `TMPDIR` to shared storage for their `srun` steps.


License
-------

SW4 is published under [GPL v2 or newer](LICENSE.txt).



## Native CPU and optional GPU builds

CMake defaults to standalone native OpenMP. Select `SW4_BACKEND=CUDA` or
`SW4_BACKEND=HIP` explicitly to build the retained RAJA GPU implementation.
RAJA/Umpire and GPU toolchains are unnecessary for the default CPU build.
See [build options and backend tests](doc/unified-backends.md).
