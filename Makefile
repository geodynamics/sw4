#-----------------------------------------------------------------------
# Usage:
#   Default is: debug=no prec=double openmp=yes hdf5=no fftw=no
# make sw4     [debug=yes/no] [prec=single/double] [openmp=yes/no] [hdf5=yes/no] [fftw=yes/no]
# make sw4mopt [debug=yes/no] [prec=single/double] [openmp=yes/no] [hdf5=yes/no] [fftw=yes/no]
#
# This Makefile asumes that the following environmental variables have been assigned,
# see note below.
# etree = [yes/no]
# proj = [yes/no]
# CXX = C++ compiler
# FC  = Fortran-77 compiler
#
# SW4ROOT = path to third party libraries (used when etree=yes and proj=yes). 
# HDF5ROOT = path to hdf5 library and include files (used when hdf5=yes).
# FFTWROOT = path to fftw library and include files (used when fftw=yes).
# Note: third party libraries should have include files in $(SW4ROOT)/include, libraries in $(SW4ROOT)/lib
# Note: HDF5ROOT and FFTWROOT can be left undefined if these libraries are 
#       available by some other mechanism, such as 'module add'.
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

ifeq ($(debug),yes)
   optlevel = DEBUG
else
   debug := "no"
   optlevel = OPTIMIZE
endif

ifeq ($(optlevel),DEBUG)
   FFLAGS    = -g
   CXXFLAGS  = -g -I../src  -DBZ_DEBUG
   CFLAGS    = -g
else
   FFLAGS   = -O3 
   CXXFLAGS = -O3 -I../src 
   CFLAGS   = -O3 
endif

fullpath := $(shell pwd)

HOSTNAME := $(shell hostname)
UNAME := $(shell uname)

debugdir := debug
optdir := optimize

SW4INC    = $(SW4ROOT)/include
SW4LIB    = $(SW4ROOT)/lib
#Default, override with configs/make.name. Preferably, FFTW is installed under SW4ROOT
FFTWHOME  = $(SW4ROOT)


emptystring := ""
foundincfile := $(emptystring)

# check if the file configs/make.inc exists?
USERMAKE := $(shell if test -r configs/make.inc; then echo "configs/make.inc"; fi)

ifeq ($(USERMAKE),configs/make.inc)
  include configs/make.inc
  foundincfile := "configs/make.inc"
else

# if configs/make.inc does not exist
  ifeq ($(UNAME),Darwin)
  # for Anders' old laptop
    ifeq ($(findstring chebyshev,$(HOSTNAME)),chebyshev)
      include configs/make.chebyshev
      foundincfile := "configs/make.chebyshev"
  # for Anders' new laptop
    else ifeq ($(findstring fourier,$(HOSTNAME)),fourier)
      include configs/make.fourier
      foundincfile := "configs/make.fourier"
   # for any other MacOS system
    else
       include configs/make.osx
       foundincfile := "configs/make.osx"
    endif
  endif

# put the variables in the configs/make.xyz file
  ifeq ($(UNAME),Linux)
# For Quartz at LC
    ifeq ($(findstring quartz,$(HOSTNAME)),quartz)
      include configs/make.quartz
      foundincfile := "configs/make.quartz"
# object code goes in machine specific directory on LC
      debugdir := debug_quartz
      optdir   := optimize_quartz
    else ifeq ($(findstring lassen,$(HOSTNAME)),lassen)
# Lassen @ LC, cpu only
      include configs/make.lassen
      foundincfile := "configs/make.lassen"
      debugdir := debug_lassen
      optdir   := optimize_lassen
# Cori @ NERSC
    else ifeq ($(findstring cori,$(HOSTNAME)),cori)
      include configs/make.cori
      foundincfile := "configs/make.cori"
  # for Bjorn's tux box
    else ifeq ($(findstring tux405,$(HOSTNAME)),tux405)
      include configs/make.tux405
      foundincfile := "configs/make.tux405"
  # for Anders' tux box
    else ifeq ($(findstring tux355,$(HOSTNAME)),tux355)
      include configs/make.tux355
      foundincfile := "configs/make.tux355"
  # For Edison at NERSC
    else ifeq ($(findstring edison,$(HOSTNAME)),edison)
      include configs/make.edison
      foundincfile := "configs/make.edison"
  # For Ray at LC, running on CPUs only
    else ifeq ($(findstring ray,$(HOSTNAME)),ray)
      include configs/make.ray
      foundincfile := "configs/make.ray"
# object code goes in machine specific directory on LC
      debugdir := debug_raycpu
      optdir := optimize_raycpu
    endif
  endif

endif

# This needs to be added before the OMP flags to work on a mac with the Apple clang compiler
ifdef EXTRA_CXX_FLAGS
   CXXFLAGS += $(EXTRA_CXX_FLAGS)
endif

# openmp=yes is default
ifeq ($(openmp),no)
   CXXFLAGS += -DSW4_NOOMP
else
   debugdir := $(debugdir)_mp
   optdir   := $(optdir)_mp
   CXXFLAGS += -fopenmp
   FFLAGS   += -fopenmp
endif

ifdef EXTRA_FORT_FLAGS
   FFLAGS += $(EXTRA_FORT_FLAGS)
endif

ifeq ($(etree),yes)
   CXXFLAGS += -DENABLE_ETREE -DENABLE_PROJ4 -I$(SW4INC)
   linklibs += -L$(SW4LIB) -lcencalvm -lproj
else ifeq ($(proj),yes)
   CXXFLAGS += -DENABLE_PROJ4 -I$(SW4INC)
   linklibs += -L$(SW4LIB) -lproj 
   etree := "no"
else
   etree := "no"
   proj  := "no"
endif

# FFTW needed for random material. If FFTWHOME undefined, it is assumed that 
#   fftw has been defined by adding a module (or similar) from the OS.
ifeq ($(fftw),yes)
   ifdef FFTWHOME
      CXXFLAGS += -DENABLE_FFTW -I$(FFTWHOME)/include 
   else
      CXXFLAGS += -DENABLE_FFTW 
   endif
   linklibs += -L$(FFTWHOME)/lib -lfftw3_mpi -lfftw3 
endif

ifeq ($(prec),single)
   debugdir := $(debugdir)_sp
   optdir   := $(optdir)_sp
   CXXFLAGS += -I../src/float
else
   CXXFLAGS += -I../src/double
endif

# hdf5=no is the default
ifeq ($(hdf5),yes)
   # PROVIDE HDF5ROOT in configs/make.xyz, e.g.
   CXXFLAGS  += -I$(HDF5ROOT)/include -DUSE_HDF5
   EXTRA_LINK_FLAGS += -L$(HDF5ROOT)/lib -lhdf5_hl -lhdf5
endif

ifdef EXTRA_LINK_FLAGS
   linklibs += $(EXTRA_LINK_FLAGS)
endif

ifeq ($(optlevel),DEBUG)
   builddir = $(debugdir)
else
   builddir = $(optdir)
endif

# routines from quadpack for numerical quadrature
QUADPACK = dqags.o dqagse.o  dqaws.o  dqawse.o  dqc25s.o \
           dqcheb.o  dqelg.o  dqk15w.o  dqk21.o  dqmomo.o \
           dqpsrt.o  dqwgts.o  qaws.o  qawse.o  qc25s.o \
           qcheb.o  qk15w.o  qmomo.o  qpsrt.o  qwgts.o xerror.o d1mach.o r1mach.o

# sw4 main program
OBJSW4 = main.o

# basic sw4 classes and functions
# The code includes MaterialInvTest, invtestmtrl and projectmtrl to build sw4 with support for sw4mopt

OBJ  = EW.o Sarray.o version.o parseInputFile.o ForcingTwilight.o \
       curvilinearGrid.o   \
       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
       solve.o  Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o  \
       TimeSeries.o sacsubc.o SuperGrid.o  TestRayleighWave.o \
       MaterialPfile.o Filter.o Polynomial.o SecondOrderSection.o time_functions.o Qspline.o \
       EtreeFile.o MaterialIfile.o GeographicProjection.o \
       Image3D.o ESSI3D.o ESSI3DHDF5.o MaterialVolimagefile.o MaterialRfile.o MaterialSfile.o  \
       AnisotropicMaterialBlock.o  sacutils.o  DataPatches.o addmemvarforcing2.o \
       consintp.o  oddIoddJinterp.o evenIoddJinterp.o MaterialInvtest.o \
       oddIevenJinterp.o evenIevenJinterp.o CheckPoint.o geodyn.o AllDims.o Patch.o RandomizedMaterial.o  \
       sw4-prof.o

# Fortran routines (lamb_exact_numquad needs QUADPACK)
 OBJ +=  rayleighfort.o lamb_exact_numquad.o 

# new C-routines converted from fortran
 OBJ += addsgdc.o bcfortc.o bcfortanisgc.o bcfreesurfcurvanic.o boundaryOpc.o energy4c.o checkanisomtrlc.o \
        computedtanisoc.o curvilinear4sgc.o gradientsc.o randomfield3dc.o innerloop-ani-sgstr-vcc.o ilanisocurvc.o \
        rhs4curvilinearc.o rhs4curvilinearsgc.o rhs4th3fortc.o solerr3c.o testsrcc.o rhs4th3windc.o \
        tw_aniso_forcec.o tw_aniso_force_ttc.o velsumc.o twilightfortc.o twilightsgfortc.o tw_ani_stiffc.o \
        anisomtrltocurvilinearc.o scalar_prodc.o updatememvarc.o addsg4windc.o bndryOpNoGhostc.o rhs4th3windc2.o


# OpenMP & C-version of the F-77 routine curvilinear4sg() is in rhs4sgcurv.o
# Source optimization
#OBJOPT = optmain.o linsolvelu.o solve-backward.o ConvParOutput.o 

# Material optimization
MOBJOPT  = moptmain.o solve-backward-allpars.o lbfgs.o nlcg.o ProjectMtrl.o \
           MaterialParameterization.o Mopt.o MaterialParCartesian.o InterpolateMaterial.o \
	   MaterialParCartesianVels.o MaterialParCartesianVp.o MParGridFile.o MaterialParCartesianVsVp.o \
           MaterialParAllpts.o

# prefix object files with build directory
FSW4 = $(addprefix $(builddir)/,$(OBJSW4))

FOBJ = $(addprefix $(builddir)/,$(OBJ)) $(addprefix $(builddir)/,$(QUADPACK))

# Source optimization
#FOBJOPT = $(addprefix $(builddir)/,$(OBJOPT)) $(addprefix $(builddir)/,$(OBJ)) $(addprefix $(builddir)/,$(QUADPACK))

# Material optimization
FMOBJOPT = $(addprefix $(builddir)/,$(MOBJOPT)) $(addprefix $(builddir)/,$(QUADPACK))


# prefix
sw4: $(FSW4) $(FOBJ)
	@echo "*** Configuration file: '" $(foundincfile) "' ***"
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " proj=" $(proj) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "CXX=" $(CXX) "EXTRA_CXX_FLAGS"= $(EXTRA_CXX_FLAGS)
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ main.o $(OBJ) $(QUADPACK) $(linklibs)
# test: linking with openmp for the routine rhs4sgcurv.o
#	cd $(builddir); $(CXX) $(CXXFLAGS) -qopenmp -o $@ main.o $(OBJ) $(QUADPACK) $(linklibs)
	@cat wave.txt
	@echo "*** Build directory: " $(builddir) " ***"

sw4mopt: $(FOBJ) $(FMOBJOPT) 
	@echo "*** Configuration file: '" $(foundincfile) "' ***"
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " proj=" $(proj) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "CXX=" $(CXX) "EXTRA_CXX_FLAGS"= $(EXTRA_CXX_FLAGS)
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ $(MOBJOPT) $(OBJ) $(QUADPACK) $(linklibs)
	@echo " "
	@echo "******* sw4mopt was built successfully *******" 
	@echo " "
	@echo "*** Build directory: " $(builddir) " ***"

# tarball
sw4-v1.1.tgz:  $(FSW4) $(FOBJ)
	rm -rf sw4-v1.1
	mkdir sw4-v1.1
	cp -r src configs tools examples doc Makefile wave.txt CMakeLists.txt INSTALL.txt LICENSE.txt README.txt sw4-v1.1
	tar czf $@ sw4-v1.1
	rm -rf sw4-v1.1 

# test
$(builddir)/rhs4sgcurv.o:src/rhs4sgcurv.C
	cd $(builddir); $(CXX) $(CXXFLAGS) -c ../$<
#	cd $(builddir); $(CXX) $(CXXFLAGS) -qopenmp -c ../$<

$(builddir)/version.o:src/version.C .FORCE
	cd $(builddir); $(CXX) $(CXXFLAGS) -DEW_MADEBY=\"$(USER)\"  -DEW_OPT_LEVEL=\"$(optlevel)\" -DEW_COMPILER=\""$(shell which $(CXX))"\" -DEW_LIBDIR=\"${SW4LIB}\" -DEW_INCDIR=\"${SW4INC}\" -DEW_HOSTNAME=\""$(shell hostname)"\" -DEW_WHEN=\""$(shell date)"\" -c ../$<

# having version.o depend on .FORCE has the effect of always building version.o
.FORCE:

$(builddir)/%.o:src/%.f
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FC_FIXED_FORMAT) $(FFLAGS) -c ../$<

$(builddir)/%.o:src/%.f90
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$<

$(builddir)/%.o:src/quadpack/%.f
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FC_FIXED_FORMAT) $(FFLAGS) -c ../$<

$(builddir)/%.o:src/%.C
	/bin/mkdir -p $(builddir)
	 cd $(builddir); $(CXX) $(CXXFLAGS) -c ../$< 

$(builddir)/RandomizedMaterial.o:src/RandomizedMaterial.C
	/bin/mkdir -p $(builddir)
	 cd $(builddir); $(CXX) $(CXXFLAGS) -std=c++11 -c ../$< 

clean:
	/bin/mkdir -p $(optdir)
	/bin/mkdir -p $(debugdir)
	cd $(optdir); /bin/rm -f sw4 sw4mopt *.o; cd ../$(debugdir); /bin/rm -f sw4 sw4mopt *.o

# Special rule for the target test
test:
	echo "Running tests..."
	/opt/local/bin/ctest --force-new-ctest-process $(ARGS)
