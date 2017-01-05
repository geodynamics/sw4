# To define the compilers on LC:
# use ic-16.0.210
# use mvapich2-intel-2.1
#
HOSTNAME := $(shell hostname)

ifeq ($(findstring fourier,$(HOSTNAME)),fourier)
  FC = mpif90
  CXX = mpicxx
  OMPOPT = -fopenmp
  EXTRA_LINK_FLAGS = -framework Accelerate -L/opt/local/lib/gcc48 -lgfortran
  openmp = yes
else ifeq ($(findstring quartz,$(HOSTNAME)),quartz)
  FC = mpifort
  CXX = mpicxx
  OMPOPT = -qopenmp
  MKL_PATH = /usr/tce/packages/mkl/mkl-11.3.3/lib
  EXTRA_CXX_FLAGS = -xCORE-AVX2
  EXTRA_FORT_FLAGS = -xCORE-AVX2
  EXTRA_LINK_FLAGS = -Wl,-rpath=$(SW4ROOT)/lib -Wl,-rpath=${MKL_PATH} -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -lifcore
  openmp = yes
else ifeq ($(findstring cab,$(HOSTNAME)),cab)
# assumes: use ic_16.0.210, use mvapich2-intel-2.1
  FC = mpif90
  CXX = mpic++
  OMPOPT = -fopenmp
  EXTRA_LINK_FLAGS = -lblas -llapack
else ifeq ($(findstring carl,$(HOSTNAME)),carl)
  FC = mpiifort
  CXX = mpiicpc
  OMPOPT = -qopenmp
  MKL_PATH = /usr/common/software/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64
  EXTRA_CXX_FLAGS = -xCORE-AVX2
  EXTRA_FORT_FLAGS = -xCORE-AVX2
  EXTRA_LINK_FLAGS = -Wl,-rpath=$(SW4ROOT)/lib -Wl,-rpath=${MKL_PATH} -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -lifcore
  openmp = yes
else ifeq ($(findstring valhall,$(HOSTNAME)),valhall)
  FC = mpif90_ifort
  CXX = mpicc_icc
  OMPOPT = -qopenmp
  EXTRA_LINK_FLAGS = -lifcore -lblas -llapack
  openmp = yes
  debug = no
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
   FFLAGS   = -O3 -g
   CXXFLAGS = -O2 -I../src -g
   CFLAGS   = -O2 -g
endif

ifeq ($(openmp),yes)
  FFLAGS   += $(OMPOPT)
  CXXFLAGS += $(OMPOPT) -DSW4_OPENMP
  CFLAGS   += $(OMPOPT) -DSW4_OPENMP
endif

ifeq ($(conly),yes)
   CXXFLAGS += -DSW4_CROUTINES
   ifeq ($(openmp),yes)
	CXXFLAGS += $(OMPOPT) -DSW4_OPENMP
   endif
endif

fullpath := $(shell pwd)

HOSTNAME := $(shell hostname)
UNAME := $(shell uname)

ifeq ($(conly),yes)
   ifeq ($(openmp),yes)
      debugdir := debug_mp_c
      optdir   := optimize_mp_c
   else
      debugdir := debug_c
      optdir   := optimize_c
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

ifdef EXTRA_CXX_FLAGS
   CXXFLAGS += $(EXTRA_CXX_FLAGS)
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

ifneq ($(conly),yes)
  OBJ += $(FOROBJ)
endif
FOBJ = $(addprefix $(builddir)/,$(OBJ))


sw4lite: $(FOBJ)
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " proj=" $(proj) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "CXX=" $(CXX) "EXTRA_CXX_FLAGS"= $(EXTRA_CXX_FLAGS)
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(linklibs)
	@cat wave.txt
	@echo "*** Build directory: " $(builddir) " ***"

$(builddir)/%.o:src/%.C
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(CXX) $(CXXFLAGS) -c ../$< 

$(builddir)/%.o:src/%.f
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$< 

$(builddir)/%.o:src/%.f90
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$< 

clean:
	/bin/mkdir -p $(optdir)
	/bin/mkdir -p $(debugdir)
	cd $(optdir); /bin/rm -f sw4lite *.o; cd ../$(debugdir); /bin/rm -f sw4lite *.o


