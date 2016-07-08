# To define the compilers on LC:
# use ic-16.0.210
# use mvapich2-intel-2.1
#
FC  = mpif90
CXX = mpic++
OMPOPT = -qopenmp
EXTRA_LINK_FLAGS = -lifcore
openmp = yes

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

ifeq ($(nofort),yes)
   CXXFLAGS += -DSW4_CROUTINES
   ifeq ($(openmp),yes)
	CXXFLAGS += $(OMPOPT) -DSW4_OPENMP
   endif
endif

fullpath := $(shell pwd)

HOSTNAME := $(shell hostname)
UNAME := $(shell uname)

ifeq ($(nofort),yes)
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

OBJ  = main.o EW.o Sarray.o Source.o SuperGrid.o GridPointSource.o time_functions.o EW_cuda.o ew-cfromfort.o rhs4sg.o EWCuda.o CheckPoint.o Parallel_IO.o

FOROBJ = rhs4th3fort.o boundaryOp.o addsgd.o solerr3.o bcfort.o

ifneq ($(nofort),yes)
  OBJ += $(FOROBJ)
endif
FOBJ = $(addprefix $(builddir)/,$(OBJ))


sw4g: $(FOBJ)
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


