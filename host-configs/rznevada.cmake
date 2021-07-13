

set (LAPACK_DIR $ENV{CRAY_LIBSCI_PREFIX_DIR} CACHE STRING "")
set (PROJ4_DIR "/usr/workspace/wsrzd/ramesh/Project6/2021/HIP/RZNevada/proj-5.0.0/install" CACHE STRING "")
set (ENABLE_MPI ON CACHE BOOL "" FORCE)
set (ENABLE_HIP ON CACHE BOOL "" FORCE)
set (RNABLE_OPENMP OFF CACHE BOOL "" FORCE)
set (ENABLE_RAJA ON CACHE BOOL "" FORCE)


set (MPI_HOME "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.5.6-aocc-4.1.0/")

set(MPI_CXX_COMPILER "${MPI_HOME}/bin/mpic++" CACHE PATH "")

set(HIPCC_VERSION "rocm-4.1.0" CACHE STRING "")
set(HIPCC_HOME $ENV{ROCM_PATH} )

set(CMAKE_C_COMPILER   "${HIPCC_HOME}/bin/hipcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "${HIPCC_HOME}/bin/hipcc" CACHE PATH "")
set(BLT_CXX_STD "c++11" CACHE STRING "")

