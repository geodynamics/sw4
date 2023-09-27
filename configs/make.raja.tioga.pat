FC = gfortran
LINKER = CC
CXX = $(PREP) CC
GIT_VERSION := "$(shell git describe --abbrev=4 --dirty --always --tags)"

proj = yes
ckernel=yes 
openmp=yes
raja_cuda=yes
umpire = yes
fftw = no
hdf5= yes
zfp=no
spill_warns = no
caliper=yes
raja_only=no
dbg=no
register_count=no

SW4_GPU_ARCH = --offload-arch=gfx90a

RAJA_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/RAJA-v0.14.1/install_crusher_510
RAJA_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/RAJA-v2022.03.0/install_tioga_560
RAJA_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/RAJA-v2022.03.0/install_tioga_560_amdclang

PROJ_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/proj-5.0.0/install_500
PROJ_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/proj-9.0.0/install_550

UMPIRE_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/umpire-6.0.0/install_510
UMPIRE_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/umpire-2022.03.1/install_tioga_560
UMPIRE_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/umpire-2022.03.1/install_tioga_560_amdclang


HDF5_HOME = $(OLCF_HDF5_ROOT)
HDF5_HOME = $(HDF5_DIR)

ZFP_HOME=/gpfs/alpine/csc300/world-shared/spock/zfp
H5Z_HOME=/gpfs/alpine/csc300/world-shared/spock/H5Z-ZFP/install

CALIPER_HOME = /gpfs/alpine/geo130/proj-shared/ramesh/2021/SPOCK/Caliper/install
CALIPER_HOME = /usr/workspace/wsrzd/ramesh/Project6/2021/HIP/NEW/DECEMBER_2021/Caliper/install_release452
CALIPER_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/Caliper/install_release510
CALIPER_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/CALIPER/Caliper-2.9.1/install
CALIPER_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/CALIPER/Caliper-2.9.1/install_mpi
CALIPER_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/CALIPER/Caliper-2.9.1/install_no_gotcha
CALIPER_HOME = /usr/workspace/wsb/ramesh/CORAL2/2022/EAS3_ACCEPTANCE_SW4/CALIPER/Caliper/install

HIP_ROOT_DIR = $(HIP_PATH)
HSA_ROOT_DIR = $(ROCM_PATH)/hsa

LAPACK_DIR = $(CRAY_LIBSCI_PREFIX_DIR)
LAPACK_LIB = ${LAPACK_DIR}/lib
LAPACK_INC = ${LAPACK_DIR}/include 

MPI_HOME = /opt/cray/pe/cray-mvapich2_nogpu/2.3.5/infiniband/cray/10.0
MPI_HOME = ${MPICH_DIR}

EXTRA_FORT_FLAGS = 

MORE_FLAGS = -DENABLE_HIP=1 -DENABLE_MPI_TIMING_BARRIER=1 -DSW4_STAGED_MPI_BUFFERS=1 -DUSE_DIRECT_INVERSE=1 -I ${LAPACK_INC} -I$(MPI_HOME)/include -I${ROCM_PATH}/include -D__HIP_ARCH_GFX90a__=1 -O3 -x hip -Wno-inconsistent-missing-override $(SW4_GPU_ARCH) --rocm-path=$(ROCM_PATH) -I $(ROCM_PATH)/include/roctracer -DCAMP_USE_PLATFORM_DEFAULT_STREAM=1 -DSW4_GIT_VERSION=\"$(GIT_VERSION)\" -DSW4_GHCOF_NO_GP_IS_ZERO=1  -ggdb -std=c++17 -Wno-unused-result -DSOURCE_INVERSION=1 -fgpu-inline-threshold=7400

MORE_LINK_FLAGS =  -O3  -std=c++17 -L $(LAPACK_LIB)  -Wl,-rpath=$(LAPACK_LIB)  -lsci_amd_mp  -L $(RAJA_HOME)/lib -lRAJA -L $(MPI_HOME)/lib -lmpi -L $(UMPIRE_HOME)/lib -lumpire $(GCC_LINK_LINE) $(SW4_GPU_ARCH)  -L ${ROCM_PATH}/roctracer/lib -lroctracer64 -L${CRAY_MPICH_ROOTDIR}/gtl/lib -Wl,-rpath=${CRAY_MPICH_ROOTDIR}/gtl/lib/ -lmpi_gtl_hsa -lcamp -lflang -lpgmath -lflangrti -ldl -L /opt/cray/pe/gcc/12.1.0/snos/lib64/ -Wl,-rpath=/opt/cray/pe/gcc/12.1.0/snos/lib64/ -lgfortran --hip-link  -L ${ROCM_PATH}/lib -lamdhip64

ifeq ($(umpire),yes)
MORE_FLAGS+= -DSW4_USE_UMPIRE=1 -I $(UMPIRE_HOME)/include
MORE_LINK_FLAGS+=  -L$(UMPIRE_HOME)/lib -lumpire
endif 

ifeq ($(fftw),yes)
MORE_FLAGS+= -DENABLE_FFTW=1 -I $(FFTW_DIR)/include
MORE_LINK_FLAGS+=-L $(FFTW_DIR)/lib -Wl,-rpath=$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3
endif 

ifeq ($(hdf5),yes)
MORE_FLAGS+= -DUSE_HDF5=1 -I $(HDF5_HOME)/include
MORE_LINK_FLAGS+= -L $(HDF5_HOME)/lib -Wl,-rpath=$(HDF5_HOME)/lib -lhdf5_hl -lhdf5
endif 
 
ifeq ($(zfp),yes)
MORE_FLAGS+= -DUSE_ZFP=1 -I $(ZFP_HOME)/include -I $(H5Z_HOME)/include
MORE_LINK_FLAGS+= -Wl,-rpath=$(ZFP_HOME)/lib-Wl,-rpath=$(H5Z_HOME)/lib -L $(H5Z_HOME)/lib -L $(ZFP_HOME)/lib -lh5zzfp -lzfp -lu -L/opt/cray/pe/cce/12.0.0/cce/x86_64/lib/
endif

ifeq ($(openmp),yes)
MORE_FLAGS+=  -fopenmp
MORE_LINK_FLAGS+= -fopenmp
else
MORE_FLAGS+= 
endif 

ifeq ($(spill_warns),yes)
MORE_FLAGS+=  -Xptxas -v,--warn-on-spills
endif 

ifeq ($(raja_only),yes)
MORE_FLAGS+=  -DRAJA_ONLY=1
endif 

ifeq ($(caliper),yes)
MORE_FLAGS+= -DENABLE_CALIPER=1 -I$(CALIPER_HOME)/include
MORE_LINK_FLAGS+= -Wl,-rpath=$(CALIPER_HOME)/lib64 -L $(CALIPER_HOME)/lib64 -lcaliper
endif 

ifeq ($(proj),yes)
MORE_FLAGS+= -I$(PROJ_HOME)/include -DENABLE_PROJ
MORE_LINK_FLAGS+=-Wl,-rpath=$(PROJ_HOME)/lib64  -L$(PROJ_HOME)/lib64 -lproj
endif 

ifeq ($(dbg),yes)
MORE_FLAGS+=  -ggdb -DPEEKS_GALORE=1 
endif 

ifeq ($(register_count),yes)
MORE_FLAGS+=  --save-temps -g -ggdb -DNO_DEVICE_FUNCTION_POINTERS=1
else
MORE_FLAGS+= -fgpu-rdc 
MORE_LINK_FLAGS+= -fgpu-rdc
endif 


LINKFLAGS = 

EXTRA_CXX_FLAGS =  -O3 $(MORE_FLAGS) -I$(RAJA_HOME)/include -DRAJA_USE_RESTRICT_PTR -DCUDA_CODE

EXTRA_LINK_FLAGS =  $(MORE_LINK_FLAGS) 

