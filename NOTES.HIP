Runs on AMD devices can be slow due to tasks using all cores when -c and OMP_NUM_THREADS are not specified.
pip install h5py might be required for pass HDF5 tests in pytest


// Using MPI with Device memory buffers:

0. ml craype-accel-amd-gfx90a
1. Build SW4 with -DSW4_MANAGED_MPI_BUFFERS=1 ( see make.raja.crusher.GPU_MPI)
2. Add -L${CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa to the link line ( see make.raja.crusher.GPU_MPI)
3. Add export MPICH_GPU_SUPPORT_ENABLED=1 to run scrit/env

