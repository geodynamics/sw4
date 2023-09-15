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

License
-------

SW4 is published under [GPL v2 or newer](LICENSE.txt).


Release Notes
-------------
**v3.0** [2023-08-30]

Various bug fixes and new features:

- Curvilinear mesh refinement.
- Read material properties in the sfile and GeoModelGrids (HDF5) formats.
- Output material properties data in the sfile format.
- Output cross-section image files in the HDF5 format.
- Output time-history data at different locations in the HDF5 format.
- Output near-surface sub-volume data (optionally with ZFP compression) in the HDF5 format.
- Some support for full-waveform inversion of the material model

