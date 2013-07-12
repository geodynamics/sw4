#-----------------------------------------------------------------------
# Usage:
# make sw4 [debug=0/1] [etree=yes/no]
#
# This Makefile asumes that the following environmental variables have been assigned:
# CXX = C++ compiler
# CC  = C-compiler
# FC  = Fortran-77 compiler
# SW4ROOT = path to third party libraries (for etree=yes). 
#
# Note: third party libraries should have include files in $(SW4ROOT)/include, libraries in $(SW4ROOT)/lib
#
# The following environmental variables are optional:
# EXTRA_CXX_FLAGS  = additional c++ compiler flags
# EXTRA_FORT_FLAGS = additional fortran compiler flags
# EXTRA_LINK_FLAGS = additional arguments to the linker
#
# There a three ways of assigning the environmental variables:
# 1) Set them in your .cshrc (or similar) file
# 2) Set them in the configs/make.inc file
# 3) Set them on the command line before running make
#
#-----------------------------------------------------------------------
# Do not make changes below this line (don't blame us if you do!)
#-----------------------------------------------------------------------

# change this to yes/no
ifndef debug
   optlevel = OPTIMIZE
else
   optlevel = DEBUG
endif

ifeq ($(optlevel),DEBUG)
   FFLAGS    = -g
   CXXFLAGS  = -g 
   CFLAGS    = -g
else
   FFLAGS   = -O3 
   CXXFLAGS = -O 
   CFLAGS   = -O 
endif

debugdir := debug_v1.0
optdir := optimize_v1.0

HOSTNAME := $(shell hostname)
UNAME := $(shell uname)

SW4INC    = $(SW4ROOT)/include
SW4LIB    = $(SW4ROOT)/lib

emptystring := ""
foundincfile := $(emptystring)

# how can we check if the file configs/make.inc exists?
ifeq ($(UNAME),Darwin)
  include configs/make.osx
  foundincfile := "configs/make.osx"
endif

# put the variables in the configs/make.xyz file
ifeq ($(UNAME),Linux)
  FC = mpiifort
  CC = mpiicc
  CXX = mpiicpc
  linklibs =  -llapack -lblas -lifcore
# for Bjorn's tux box
  ifeq ($(findstring tux270,$(HOSTNAME)),tux270)
    linklibs = -L/home/sjogreen2/apps/lapack-3.2 -llapack -lblas -lgfortran
  endif
# for Anders' tux box
  ifeq ($(findstring tux355,$(HOSTNAME)),tux355)
    linklibs = -L/home/andersp/lib -llapack -lblas -lgfortran
    FC = gfortran
    CC = mpicc
    CXX = mpicxx
  endif
# For Edison at NERSC
  ifeq ($(findstring edison,$(HOSTNAME)),edison)
    FC  = ftn
    CC  = cc
    CXX = CC
    linklibs = -L/global/u2/b/bjorns/local/lapack-3.4.0 -llapack -lrefblas -lifcore
  endif
  ifeq ($(etree),yes)
     CXXFLAGS += -DENABLE_ETREE
     CXXFLAGS += -I$(SW4INC)
     linklibs += -Wl,-rpath=$(SW4LIB) -L$(SW4LIB) -lcencalvm -lproj
  endif
endif

ifdef EXTRA_CXX_FLAGS
   CXXFLAGS += $(EXTRA_CXX_FLAGS)
endif

ifdef EXTRA_FORT_FLAGS
   FFLAGS += $(EXTRA_FORT_FLAGS)
endif

ifeq ($(etree),yes)
   CXXFLAGS += -DENABLE_ETREE -I$(SW4INC)
   linklibs += -lcencalvm -lproj
endif

ifdef EXTRA_LINK_FLAGS
   linklibs += $(EXTRA_LINK_FLAGS)
endif

ifeq ($(optlevel),DEBUG)
   builddir = $(debugdir)
else
   builddir = $(optdir)
endif

QUADPACK = dqags.o dqagse.o  dqaws.o  dqawse.o  dqc25s.o \
           dqcheb.o  dqelg.o  dqk15w.o  dqk21.o  dqmomo.o \
           dqpsrt.o  dqwgts.o  qaws.o  qawse.o  qc25s.o \
           qcheb.o  qk15w.o  qmomo.o  qpsrt.o  qwgts.o xerror.o d1mach.o r1mach.o

OBJ  = main.o EW.o Sarray.o version.o parseInputFile.o ForcingTwilight.o \
       curvilinearGrid.o boundaryOp.o bcfort.o twilightfort.o rhs4th3fort.o \
       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
       solve.o solerr3.o Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o testsrc.o \
       TimeSeries.o sacsubc.o SuperGrid.o addsgd.o velsum.o rayleighfort.o energy4.o TestRayleighWave.o \
       MaterialPfile.o Filter.o Polynomial.o SecondOrderSection.o time_functions.o Qspline.o \
       lamb_exact_numquad.o twilightsgfort.o EtreeFile.o MaterialIfile.o GeographicProjection.o \
       rhs4curvilinear.o rhs4curvilinearsg.o gradients.o Image3D.o MaterialVolimagefile.o MaterialInvtest.o \
       invtestmtrl.o projectmtrl.o

# prefix object files with build directory
FOBJ = $(addprefix $(builddir)/,$(OBJ)) $(addprefix $(builddir)/,$(QUADPACK))

ifeq ($(foundincfile), $(emptystring))
sw4:: $(FOBJ) $(builddir)
	@echo "ERROR: configuarion file not found"
else
sw4: $(FOBJ) $(builddir)
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(QUADPACK) $(linklibs)
	@cat wave.txt
	@echo "*** Configuration file: " $(foundincfile) "***"
	@echo "*** Build directory: " $(builddir) " ***"
endif

$(builddir):
	/bin/mkdir -p $(builddir)

$(builddir)/version.o:src/version.C $(builddir) .FORCE
	cd $(builddir); $(CXX) $(CXXFLAGS) -DEW_MADEBY=\"$(USER)\"  -DEW_OPT_LEVEL=\"$(optlevel)\" -DEW_COMPILER=\""$(shell which $(CXX))"\" -DEW_LIBDIR=\"${SW4LIB}\" -DEW_INCDIR=\"${SW4INC}\" -DEW_HOSTNAME=\""$(shell hostname)"\" -DEW_WHEN=\""$(shell date)"\" -c ../$<

# having version.o depend on .FORCE has the effect of always building version.o
.FORCE:

$(builddir)/%.o:src/%.f
	cd $(builddir); $(FC) $(FFLAGS) -c ../$<

$(builddir)/%.o:src/quadpack/%.f
	cd $(builddir); $(FC) $(FFLAGS) -c ../$<

$(builddir)/%.o:src/%.C
	 cd $(builddir); $(CXX) $(CXXFLAGS) -c ../$< 

clean:
	/bin/mkdir -p $(optdir)
	/bin/mkdir -p $(debugdir)
	cd $(optdir); /bin/rm -f sw4 $(OBJ); cd ../$(debugdir); /bin/rm -f sw4 $(OBJ)
