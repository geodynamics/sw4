  etree = no
  proj = yes
  SW4ROOT = /usr/apps/wpp
# assumes module load intel
# FC = mpiifort
  FC = mpif90
  CXX = mpicxx -std=c++11
  MKL_PATH = /usr/tce/packages/mkl/mkl-2019.0/lib
  HDF5ROOT   = ${HDF5}
  FFTWHOME = /usr/tce/packages/fftw/fftw-3.3.8-mvapich2-2.3-intel-19.0.3

EXTRA_CXX_FLAGS = -std=c++11
EXTRA_LINK_FLAGS = -Wl,-rpath=$(FFTWHOME)/lib -Wl,-rpath=$(SW4ROOT)/lib -Wl,-rpath=${MKL_PATH} -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -lifcore

    
#  EXTRA_CXX_FLAGS = -O
#  EXTRA_FORT_FLAGS = -O
# must load appropriate modules for hdf5-parallel library support
# ( module add hdf5-parallel )
#

  debugdir := debug
  optdir := optimize
