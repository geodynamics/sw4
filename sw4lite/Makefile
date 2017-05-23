#-----------------------------------------------------------------------
# Usage:
# make sw4lite [ckernel=yes/no single=yes/no] 
# (default is ckernel=no, single=no, i.e. use the Fortran kernel and double precision)
#
# To define the compilers on LC:
# use ic-16.0.210
# use mvapich2-intel-2.1
#

ifneq ($(ckernel),yes)
   ckernel := "no"
endif

ifneq ($(single),yes)
   single := "no"
endif

HOSTNAME := $(shell hostname)

# default optimization level (may be modified below)
OPT = -O3

# Fourier is Anders' mac powerbook
ifeq ($(findstring fourier,$(HOSTNAME)),fourier)
  FC = mpif90
  CC = gcc
  CXX = mpicxx
  OMPOPT = -fopenmp
  EXTRA_LINK_FLAGS = -framework Accelerate -L/opt/local/lib/gcc5 -lgfortran
  openmp = yes
# LC quartz is a large cluster of Intel Haswell nodes
  computername := fourier
else ifeq ($(findstring quartz,$(HOSTNAME)),quartz)
  FC = mpifort
  CXX = mpicxx
  OMPOPT = -fopenmp -std=c++11 
  MKL_PATH = /usr/tce/packages/mkl/mkl-11.3.3/lib
  EXTRA_CXX_FLAGS = -xCORE-AVX2 -I /usr/workspace/wsb/ramesh/Quartz/R2/RAJA/install/include
  EXTRA_FORT_FLAGS = -xCORE-AVX2
  EXTRA_LINK_FLAGS = -fopenmp -Wl,-rpath=$(SW4ROOT)/lib -Wl,-rpath=${MKL_PATH} -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -lifcore -L /usr/workspace/wsb/ramesh/Quartz/R2/RAJA/install/lib -lRAJA
  openmp = yes
  LINKER = mpicxx
  computername := quartz
 RAJA_LOCATION=/usr/workspace/wsb/ramesh/Quartz/RAJA_APRIL_10/RAJA/install
else ifeq ($(findstring ray,$(HOSTNAME)),ray)
  FC = mpif90
  CXX = nvcc
  OPT = -O2 -ccbin mpixlC -std=c++11 --expt-extended-lambda -restrict -arch=sm_60 -I $(CUDA_INCLUDES) -I /g/g0/ramesh/Project6/Ray/RAJA/install/include  --x cu -DUSE_NVTX -DRAJA_USE_CUDA -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -DCUDA_CODE 
  OMPOPT = 
  EXTRA_LINK_FLAGS = -L /usr/tcetmp/packages/xl/xl-beta-2017.03.28/xlf/16.1.0/lib/ -L /usr/tcetmp/packages/blas/blas-3.6.0-xlf-15.1.5/lib -L /usr/tcetmp/packages/lapack/lapack-3.6.0-xlf-15.1.5/lib/ -lxlf90 -llapack -lblas  -L /usr/local/cuda/lib64 -lcudart -lnvToolsExt -lcuda
#EXTRA_LINK_FLAGS = -L /usr/tcetmp/packages/blas/blas-3.6.0-gfortran-4.8.5/lib/ -L /usr/tcetmp/packages/lapack/lapack-3.6.0-gfortran-4.8.5/lib/ -llapack -lblas -lm -lgfortran  -L /usr/local/cuda/lib64 -lcudart -lnvToolsExt -lcuda
 RAJA_LOCATION =  /g/g0/ramesh/Project6/Ray/RAJA/install/
  LINKER = mpixlC
 LINKFLAGS = 
  computername := ray
# LC cab is a large cluster of Intel Xeon nodes
else ifeq ($(findstring cab,$(HOSTNAME)),cab)
# assumes: use ic_16.0.210, use mvapich2-intel-2.1
  FC = mpif90
  CXX = mpic++
  OMPOPT = -fopenmp
  EXTRA_LINK_FLAGS = -lblas -llapack
  computername := cab
# NERSC carl is a single node Intel KNL machine
else ifeq ($(findstring carl,$(HOSTNAME)),carl)
  FC = mpiifort
  CXX = mpiicpc
  OMPOPT = -qopenmp
  MKL_PATH = /usr/common/software/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64
#  EXTRA_CXX_FLAGS = -xmic-avx512
  EXTRA_CXX_FLAGS = -xmic-avx2
#  EXTRA_CXX_FLAGS = -no-vec -no-simd
  EXTRA_FORT_FLAGS = -xmic-avx512,core-avx2
  EXTRA_LINK_FLAGS = -Wl,-rpath=$(SW4ROOT)/lib -Wl,-rpath=${MKL_PATH} -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -lifcore
  openmp = yes
  computername := carl
#
# LBL quadknl is Hans' single node Intel KNL machine
#
else ifeq ($(findstring quadknl,$(HOSTNAME)),quadknl)
  FC = mpiifort
  CXX = mpiicpc
  CC = icc
  OMPOPT = -qopenmp
  MKL_PATH = /opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64
  OTHER_PATH = /opt/intel/compilers_and_libraries_2017/linux/lib/intel64
  BASIC_PATH = /opt/intel/lib/intel64
  EXTRA_CXX_FLAGS = -xmic-avx512
  EXTRA_C_FLAGS = -xmic-avx512
  EXTRA_FORT_FLAGS = -xmic-avx512 
  EXTRA_LINK_FLAGS = -Wl,-rpath=${MKL_PATH} -Wl,-rpath=${BASIC_PATH} -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -Wl,-rpath=${OTHER_PATH} -L${OTHER_PATH} -lifcore -L${BASIC_PATH} -limf -lsvml -lintlc
# for building testil
#  EXTRA_LINK_FLAGS = -Wl,-rpath=${OTHER_PATH} -L${OTHER_PATH} -lifcore -L${BASIC_PATH} -limf -lsvml -lintlc -lm -ldl 
  openmp = yes
  computername := quadknl
#
# Cori
else ifeq ($(findstring cori,$(HOSTNAME)),cori)
# cray compiler wrappers
  FC = ftn
  CXX = CC
# Cray compiler:  (requires modules: PrgEnv-cray and cray-mpich)
# rhs4sgcurv_rev.C gives segfault with -O, try -O1 instead
# the ftn compiler does not accept -O. Need to specify -O2
#  OPT = -O2
#  EXTRA_LINK_FLAGS =  # only needs module load cray-mpich
# Intel compiler:
  OPT = -O3 -qoverride-limits
  OMPOPT = -qopenmp
  EXTRA_LINK_FLAGS =  -lpthread -lm -ldl -lifcore
  openmp = yes
  computername := cori
else ifeq ($(findstring valhall,$(HOSTNAME)),valhall)
  FC = mpif90_ifort
  CXX = mpicc_icc
  OMPOPT = -qopenmp
  EXTRA_LINK_FLAGS = -lifcore -lblas -llapack
  openmp = yes
  computername := valhall
else
  FC  = mpif90
  CXX = mpic++
  OMPOPT = -fopenmp
  EXTRA_LINK_FLAGS = -lgfortran -lblas -llapack
  openmp = no
  debug = no
endif

ifeq ($(debug),yes)
   optlevel = DEBUG
else
   debug := "no"
   optlevel = OPTIMIZE
endif

ifeq ($(optlevel),DEBUG)
   FFLAGS    = -g 
   CXXFLAGS  = -g -I../src -DBZ_DEBUG
   CFLAGS    = -g
else
   FFLAGS   = $(OPT)
   CXXFLAGS = $(OPT) -I../src
   CFLAGS   = $(OPT)
endif

ifeq ($(openmp),yes)
  FFLAGS   += $(OMPOPT)
  CXXFLAGS += $(OMPOPT) -DSW4_OPENMP
  CFLAGS   += $(OMPOPT) -DSW4_OPENMP
endif

ifeq ($(ckernel),yes)
   CXXFLAGS += -DSW4_CROUTINES
endif

ifeq ($(ckernel),yes)
   ifeq ($(openmp),yes)
	ifeq ($(raja),yes)
		debugdir:= debug_mp_c_raja
		optdir:= optimize_mp_c_raja
		CXXFLAGS += -std=c++11 -I $(RAJA_LOCATION)/include  -DSW4_CROUTINES 
		EXTRA_LINK_FLAGS += -qopenmp -L $(RAJA_LOCATION)/lib -lRAJA 
	else	
		debugdir := debug_mp_c
		optdir   := optimize_mp_c
	endif
   else
	ifeq ($(raja),yes)
		debugdir:= debug_c_raja
		optdir:= optimize_c_raja
		CXXFLAGS += -I $(RAJA_LOCATION)/include  -DSW4_CROUTINES 
		EXTRA_LINK_FLAGS += -L $(RAJA_LOCATION)/lib -lRAJA 
	else	
		debugdir := debug_mp_c
		optdir   := optimize_mp_c
	endif
    endif
else
   ifeq ($(openmp),yes)
      debugdir := debug_mp
      optdir   := optimize_mp
   else
      debugdir := debug
      optdir   := optimize
   endif
endif

ifeq ($(single),yes)
   debugdir := $(debugdir)_sp
   optdir   := $(optdir)_sp
   CXXFLAGS += -I../src/float
else
   CXXFLAGS += -I../src/double
endif

ifneq ($(origin computername),undefined)
   debugdir := $(debugdir)_$(computername)
   optdir   := $(optdir)_$(computername)
endif

ifdef EXTRA_CXX_FLAGS
   CXXFLAGS += $(EXTRA_CXX_FLAGS)
endif

ifdef EXTRA_C_FLAGS
   CFLAGS += $(EXTRA_C_FLAGS)
endif

ifdef EXTRA_FORT_FLAGS
   FFLAGS += $(EXTRA_FORT_FLAGS)
endif

ifdef EXTRA_LINK_FLAGS
   linklibs += $(EXTRA_LINK_FLAGS)
endif

ifeq ($(optlevel),DEBUG)
   builddir = $(debugdir)
else
   builddir = $(optdir)
endif

OBJ  = main.o EW.o Sarray.o Source.o SuperGrid.o GridPointSource.o time_functions.o EW_cuda.o ew-cfromfort.o rhs4sg.o rhs4sg_rev.o EWCuda.o CheckPoint.o Parallel_IO.o EW-dg.o MaterialData.o MaterialBlock.o Polynomial.o SecondOrderSection.o Filter.o TimeSeries.o sacsubc.o curvilinear-c.o rhs4sgcurv.o rhs4sgcurv_rev.o

FOROBJ = rhs4th3fort.o boundaryOp.o addsgd.o solerr3.o bcfort.o type_defs.o lglnodes.o dgmodule.o metric.o curvilinear4sg.o freesurfcurvisg.o

ifneq ($(ckernel),yes)
  OBJ += $(FOROBJ)
endif
FOBJ = $(addprefix $(builddir)/,$(OBJ))

sw4lite: $(FOBJ)
	@echo "********* User configuration variables **************"
	@echo "ckernel=" $(ckernel) " debug=" $(debug) " proj=" $(proj) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "CXX=" $(CXX) "EXTRA_CXX_FLAGS"= $(EXTRA_CXX_FLAGS)
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(LINKER) $(LINKFLAGS) -o $@ $(OBJ) $(linklibs)
	@cat wave.txt
	@echo "*** Build directory: " $(builddir) " ***"

OBJ2 =  testil.o rhs4sg.o rhs4sg_rev.o gettimec.o rhs4th3fort.o boundaryOp.o
FOBJ2 = $(addprefix $(builddir)/,$(OBJ2))

testil : $(FOBJ2)
	$(CXX) $(CXXFLAGS) -o $@ $(FOBJ2) $(linklibs)

$(builddir)/%.o:src/%.C
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(CXX) $(CXXFLAGS) -c ../$< 

$(builddir)/%.o:src/%.c
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(CC) $(CFLAGS) -c ../$< 

$(builddir)/%.o:src/%.f
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$< 

$(builddir)/%.o:src/%.f90
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$< 

clean:
	/bin/mkdir -p $(optdir)
	/bin/mkdir -p $(debugdir)
	/bin/mkdir -p optimize_mp_c
	cd $(optdir); /bin/rm -f sw4lite *.o; cd ../$(debugdir); /bin/rm -f sw4lite *.o; cd ../optimize_mp_c; /bin/rm -f sw4lite *.o


