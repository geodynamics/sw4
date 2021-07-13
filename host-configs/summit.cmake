set (LAPACK_DIR $ENV{OLCF_NETLIB_LAPACK_ROOT} CACHE STRING "")
set (PROJ4_DIR "/gpfs/alpine/geo130/proj-shared/AUGUST_CAMPAIGN/proj-5.0.0/install" CACHE STRING "")
set (ENABLE_MPI ON CACHE BOOL "" FORCE)
set (ENABLE_CUDA ON CACHE BOOL "" FORCE)
set (ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
set (ENABLE_RAJA ON CACHE BOOL "" FORCE)


set (CUDA_ARCH "sm_70" CACHE STRING "" FORCE)
set (MPI_HOME $ENV{OLCF_SPECTRUM_MPI_ROOT})

set(MPI_CXX_COMPILER "${MPI_HOME}/bin/mpic++" CACHE PATH "")

set(GCC_VERSION "gcc-8.3.1" CACHE STRING "")
set(GCC_HOME $ENV{OLCF_GCC_ROOT} )

set(CMAKE_C_COMPILER   "${GCC_HOME}/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "${GCC_HOME}/bin/g++" CACHE PATH "")
set(BLT_CXX_STD "c++11" CACHE STRING "")

#------------------------------------------------------------------------------
# CUDA support
#------------------------------------------------------------------------------


set(CUDA_TOOLKIT_ROOT_DIR "$ENV{CUDA_TOOLKIT_ROOT_DIR}" CACHE PATH "")
set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" CACHE PATH "")
set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE PATH "")

#set(_cuda_arch "sm_70")
set(CMAKE_CUDA_FLAGS "-lineinfo -restrict --expt-extended-lambda " CACHE STRING "")

set(CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "" FORCE)

# nvcc does not like gtest's 'pthreads' flag
set(gtest_disable_pthreads ON CACHE BOOL "")
