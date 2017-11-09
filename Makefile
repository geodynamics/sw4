#-----------------------------------------------------------------------
# Usage:
# make sw4 [debug=yes/no] [prec=single/double] [fortran=yes/no] [openmp=yes/no]
#   Default is: debug=no prec=double fortran=yes openmp=yes
# This Makefile asumes that the following environmental variables have been assigned:
# etree = [yes/no]
# proj = [yes/no]
# CXX = C++ compiler
# FC  = Fortran-77 compiler
# SW4ROOT = path to third party libraries (used when etree=yes). 
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
   FFLAGS   = -O3
# AP (160419) Note that cmake uses -O3 instead of -O for CXX and C
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
    ifeq ($(findstring yorkville,$(HOSTNAME)),yorkville)
      include configs/make.yorkville
      foundincfile := "configs/make.yorkville"
  # for Anders' new laptop
    else ifeq ($(findstring fourier,$(HOSTNAME)),fourier)
      include configs/make.fourier
      foundincfile := "configs/make.fourier"
    endif
  endif
  
  # put the variables in the configs/make.xyz file
  ifeq ($(UNAME),Linux)
  # For Cab at LC
    ifeq ($(findstring cab,$(HOSTNAME)),cab)
      include configs/make.cab
      foundincfile := "configs/make.cab"
# object code goes in machine specific directory on LC
      debugdir := debug_cab
      optdir := optimize_cab
  # for Sierra (old Sierra) on LC, almost same machine as Cab
    else ifeq ($(findstring sierra,$(HOSTNAME)),sierra)
      include configs/make.sierra
      foundincfile := "configs/make.sierra"
      debugdir := debug_sierra
      optdir := optimize_sierra
# For Quartz at LC (why doesn't this work when HOSTNAME is quartz770 ?
    else ifeq ($(findstring quartz,$(HOSTNAME)),quartz)
      include configs/make.quartz
      foundincfile := "configs/make.quartz"
    else ifeq ($(findstring cori,$(HOSTNAME)),cori)
      include configs/make.cori
      foundincfile := "configs/make.cori"
  # for Bjorn's tux box
    else ifeq ($(findstring tux337,$(HOSTNAME)),tux337)
      include configs/make.tux337
      foundincfile := "configs/make.tux337"
  # for Anders' tux box
    else ifeq ($(findstring tux355,$(HOSTNAME)),tux355)
      include configs/make.tux355
      foundincfile := "configs/make.tux355"
  # For Edison at NERSC
    else ifeq ($(findstring edison,$(HOSTNAME)),edison)
      include configs/make.edison
      foundincfile := "configs/make.edison"
  # For Vulcan at LC
    else ifeq ($(findstring vulcan,$(HOSTNAME)),vulcan)
      include configs/make.bgq
      foundincfile := "configs/make.bgq"
# object code goes in machine specific directory on LC
      debugdir := debug_vulcan
      optdir := optimize_vulcan
    endif
  endif

endif

ifdef EXTRA_CXX_FLAGS
   CXXFLAGS += $(EXTRA_CXX_FLAGS)
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

# openmp=yes is default
ifeq ($(openmp),no)
   CXXFLAGS += -DSW4_NOOMP
else
   debugdir := $(debugdir)_mp
   optdir   := $(optdir)_mp
   CXXFLAGS += -fopenmp
   FFLAGS   += -fopenmp
endif

ifeq ($(fortran),no)
else
   debugdir := $(debugdir)_fort
   optdir   := $(optdir)_fort
   CXXFLAGS += -DSW4_NOC
endif

ifeq ($(prec),single)
   debugdir := $(debugdir)_sp
   optdir   := $(optdir)_sp
   CXXFLAGS += -I../src/float
else
   CXXFLAGS += -I../src/double
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

# sw4 main program (kept separate)
OBJSW4 = main.o

OBJ  = EW.o Sarray.o version.o parseInputFile.o ForcingTwilight.o \
       curvilinearGrid.o boundaryOp.o bndryOpNoGhost.o  bcfort.o twilightfort.o rhs4th3fort.o \
       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
       solve.o solerr3.o Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o testsrc.o \
       TimeSeries.o sacsubc.o SuperGrid.o addsgd.o velsum.o rayleighfort.o energy4.o TestRayleighWave.o \
       MaterialPfile.o Filter.o Polynomial.o SecondOrderSection.o time_functions.o Qspline.o \
       lamb_exact_numquad.o twilightsgfort.o EtreeFile.o MaterialIfile.o GeographicProjection.o \
       rhs4curvilinear.o curvilinear4.o rhs4curvilinearsg.o curvilinear4sg.o gradients.o Image3D.o \
       MaterialVolimagefile.o MaterialRfile.o randomfield3d.o innerloop-ani-sgstr-vc.o bcfortanisg.o \
       AnisotropicMaterialBlock.o checkanisomtrl.o computedtaniso.o sacutils.o ilanisocurv.o \
       anisomtrltocurvilinear.o bcfreesurfcurvani.o tw_ani_stiff.o tw_aniso_force.o tw_aniso_force_tt.o \
       rhs4th3fortwind.o updatememvar.o addmemvarforcing2.o addsg4wind.o consintp.o scalar_prod.o

# new C-routines converted from fortran
 OBJ += addsgdc.o bcfortc.o bcfortanisgc.o bcfreesurfcurvanic.o boundaryOpc.o energy4c.o checkanisomtrlc.o \
        computedtanisoc.o curvilinear4sgc.o gradientsc.o randomfield3dc.o innerloop-ani-sgstr-vcc.o ilanisocurvc.o \
        rhs4curvilinearc.o rhs4curvilinearsgc.o rhs4th3fortc.o rhs4th3fortwindc.o solerr3c.o testsrcc.o \
        tw_aniso_forcec.o tw_aniso_force_ttc.o velsumc.o twilightfortc.o twilightsgfortc.o tw_ani_stiffc.o \
        anisomtrltocurvilinearc.o scalar_prodc.o updatememvarc.o addsg4windc.o


# OpenMP & C-version of the F-77 routine curvilinear4sg() is in rhs4sgcurv.o

# prefix object files with build directory
FSW4 = $(addprefix $(builddir)/,$(OBJSW4))
FOBJ = $(addprefix $(builddir)/,$(OBJ)) $(addprefix $(builddir)/,$(QUADPACK))

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
	cd $(builddir); $(FC) $(FFLAGS) -c ../$<

$(builddir)/%.o:src/%.f90
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$<

$(builddir)/%.o:src/quadpack/%.f
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$<

$(builddir)/%.o:src/%.C
	/bin/mkdir -p $(builddir)
	 cd $(builddir); $(CXX) $(CXXFLAGS) -c ../$< 

clean:
	/bin/mkdir -p $(optdir)
	/bin/mkdir -p $(debugdir)
	cd $(optdir); /bin/rm -f sw4 *.o; cd ../$(debugdir); /bin/rm -f sw4 *.o

# Special rule for the target test
test:
	echo "Running tests..."
	/opt/local/bin/ctest --force-new-ctest-process $(ARGS)
