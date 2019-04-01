#-----------------------------------------------------------------------
# Usage:
# make sw4mopt [debug=yes/no]
# make sw4     [debug=yes/no]
# make sw4opt  [debug=yes/no]
#
# This Makefile asumes that the following environmental variables have been assigned:
# etree = [yes/no]
# proj = [yes/no]
# CXX = C++ compiler
# FC  = Fortran-77 compiler
# SW4ROOT = path to third party libraries (used when etree=yes and proj=yes). 
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
   CXXFLAGS  = -g -I../src-sw4 -I../src-sw4/double
   CFLAGS    = -g
else
   FFLAGS   = -O3 
   CXXFLAGS = -O -I../src-sw4 -I../src-sw4/double
   CFLAGS   = -O 
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
# For Quartz at LC
    ifeq ($(findstring quartz,$(HOSTNAME)),quartz)
      include configs/make.quartz
      foundincfile := "configs/make.quartz"
# object code goes in machine specific directory on LC
      debugdir := debug_quartz
      optdir   := optimize_quartz
    else ifeq ($(findstring lassen,$(HOSTNAME)),lassen)
# Lassen @ LC
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

# openmp=yes is default
ifeq ($(openmp),no)
   CXXFLAGS += -DSW4_NOOMP
else
   debugdir := $(debugdir)_mp
   optdir   := $(optdir)_mp
   CXXFLAGS += -fopenmp
   FFLAGS   += -fopenmp
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

# sw4 main program (kept separate)
OBJSW4 = main.o

# basic sw4 classes and functions
# Note that include MaterialInvTest, invtestmtrl and projectmtrl to allow sw4 to be built with -DENABLE_OPT
#OBJ  = EW.o Sarray.o version.o parseInputFile.o ForcingTwilight.o \
#       curvilinearGrid.o boundaryOp.o bcfort.o twilightfort.o rhs4th3fort.o \
#       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
#       solve.o solerr3.o Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o testsrc.o \
#       TimeSeries.o sacsubc.o SuperGrid.o addsgd.o velsum.o rayleighfort.o energy4.o TestRayleighWave.o \
#       MaterialPfile.o Filter.o Polynomial.o SecondOrderSection.o time_functions.o Qspline.o \
#       lamb_exact_numquad.o twilightsgfort.o EtreeFile.o MaterialIfile.o GeographicProjection.o \
#       rhs4curvilinear.o curvilinear4.o rhs4curvilinearsg.o curvilinear4sg.o gradients.o Image3D.o \
#       MaterialVolimagefile.o MaterialRfile.o randomfield3d.o innerloop-ani-sgstr-vc.o bcfortanisg.o \
#       AnisotropicMaterialBlock.o checkanisomtrl.o computedtaniso.o \
#       MaterialInvtest.o invtestmtrl.o projectmtrl.o 

OBJ  = EW.o Sarray.o version.o parseInputFile.o ForcingTwilight.o \
       curvilinearGrid.o   \
       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
       solve.o  Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o  \
       TimeSeries.o sacsubc.o SuperGrid.o  rayleighfort.o TestRayleighWave.o \
       MaterialPfile.o Filter.o Polynomial.o SecondOrderSection.o time_functions.o Qspline.o \
       lamb_exact_numquad.o  EtreeFile.o MaterialIfile.o GeographicProjection.o \
       Image3D.o ESSI3D.o ESSI3DHDF5.o MaterialVolimagefile.o MaterialRfile.o   \
       AnisotropicMaterialBlock.o  sacutils.o  DataPatches.o addmemvarforcing2.o \
       consintp.o  oddIoddJinterp.o evenIoddJinterp.o MaterialInvtest.o \
       oddIevenJinterp.o evenIevenJinterp.o CheckPoint.o geodyn.o AllDims.o Patch.o RandomizedMaterial.o

#removed: boundaryOp.o bndryOpNoGhost.o bcfort.o twilightfort.o rhs4th3fort.o solerr3.o testsrc.o \
#         addsgd.o velsum.o energy4.o twilightsgfort.o rhs4curvilinear.o curvilinear4.o rhs4curvilinearsg.o \
#         curvilinear4sg.o gradients.o randomfield3d.o innerloop-ani-sgstr-vc.o bcfortanisg.o checkanisomtrl.o \
#         computedtaniso.o ilanisocurv.o \
#         anisomtrltocurvilinear.o bcfreesurfcurvani.o tw_ani_stiff.o tw_aniso_force.o tw_aniso_force_tt.o \
#         updatememvar.o addsg4wind.o scalar_prod.o
# new C-routines converted from fortran
 OBJ += addsgdc.o bcfortc.o bcfortanisgc.o bcfreesurfcurvanic.o boundaryOpc.o energy4c.o checkanisomtrlc.o \
        computedtanisoc.o curvilinear4sgc.o gradientsc.o randomfield3dc.o innerloop-ani-sgstr-vcc.o ilanisocurvc.o \
        rhs4curvilinearc.o rhs4curvilinearsgc.o rhs4th3fortc.o solerr3c.o testsrcc.o rhs4th3windc.o \
        tw_aniso_forcec.o tw_aniso_force_ttc.o velsumc.o twilightfortc.o twilightsgfortc.o tw_ani_stiffc.o \
        anisomtrltocurvilinearc.o scalar_prodc.o updatememvarc.o addsg4windc.o bndryOpNoGhostc.o rhs4th3windc2.o

# source optimization
OBJOPT = optmain.o linsolvelu.o solve-backward.o ConvParOutput.o 

# material optimization
MOBJOPT  = moptmain.o solve-backward-allpars.o lbfgs.o nlcg.o ProjectMtrl.o \
           MaterialParameterization.o Mopt.o MaterialParCartesian.o InterpolateMaterial.o \
	   MaterialParCartesianVels.o MaterialParCartesianVp.o
           

OBJL  = lamb_one_point.o

OBJP = test_proj.o

OBJCONV = convert_etree.o

OBJCONV2 = convert_small.o

# prefix object files with build directory
FSW4 = $(addprefix $(builddir)/,$(OBJSW4))

FOBJ = $(addprefix $(builddir)/,$(OBJ)) $(addprefix $(builddir)/,$(QUADPACK))

FOBJOPT = $(addprefix $(builddir)/,$(OBJOPT)) $(addprefix $(builddir)/,$(OBJ)) $(addprefix $(builddir)/,$(QUADPACK))

FMOBJOPT = $(addprefix $(builddir)/,$(MOBJOPT)) $(addprefix $(builddir)/,$(QUADPACK))

FOBJL = $(addprefix $(builddir)/,$(OBJL)) $(addprefix $(builddir)/,$(QUADPACK))

FOBJP = $(addprefix $(builddir)/,$(OBJP))

FOBJCONV = $(addprefix $(builddir)/,$(OBJCONV))

FOBJCONV2 = $(addprefix $(builddir)/,$(OBJCONV2))

CXXFLAGS += -I../src-opt


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

sw4opt: $(FOBJ) $(FOBJOPT)
	@echo "*** Configuration file: '" $(foundincfile) "' ***"
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " proj=" $(proj) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "CXX=" $(CXX) "EXTRA_CXX_FLAGS"= $(EXTRA_CXX_FLAGS)
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ $(OBJOPT) $(OBJ) $(QUADPACK) $(linklibs)
	@echo " "
	@echo "******* sw4opt was built successfully *******" 
	@echo " "
	@echo "*** Build directory: " $(builddir) " ***"

sw4: $(FSW4) $(FOBJ)
	@echo "*** Configuration file: '" $(foundincfile) "' ***"
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " proj=" $(proj) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "CXX=" $(CXX) "EXTRA_CXX_FLAGS"= $(EXTRA_CXX_FLAGS)
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ main.o $(OBJ) $(QUADPACK) $(linklibs)
	@cat wave.txt
	@echo "*** Build directory: " $(builddir) " ***"

lamb1: $(FOBJL)
	@echo "*** Configuration file: '" $(foundincfile) "' ***"
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(FC) $(FFLAGS) -o $@ $(OBJL) $(QUADPACK) $(linklibs)
#	@cat "Done building lamb1 executable"
#	@echo "*** Build directory: " $(builddir) " ***"

test_proj: $(FOBJP)
	@echo "*** Configuration file: '" $(foundincfile) "' ***"
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ $(OBJP) $(linklibs)

convert_etree: $(FOBJCONV)
	EXTRA_CXX_FLAGS += -I../src-opt
	@echo "*** Configuration file: '" $(foundincfile) "' ***"
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ $(OBJCONV) $(linklibs)

convert_small: $(FOBJCONV2)
	EXTRA_CXX_FLAGS += -I../src-opt
	@echo "*** Configuration file: '" $(foundincfile) "' ***"
	@echo "********* User configuration variables **************"
	@echo "debug=" $(debug) " etree=" $(etree) " SW4ROOT"= $(SW4ROOT) 
	@echo "FC=" $(FC) " EXTRA_FORT_FLAGS=" $(EXTRA_FORT_FLAGS)
	@echo "EXTRA_LINK_FLAGS"= $(EXTRA_LINK_FLAGS)
	@echo "******************************************************"
	cd $(builddir); $(CXX) $(CXXFLAGS) -o $@ $(OBJCONV2) $(linklibs)

$(builddir)/version.o:src-sw4/version.C .FORCE
	cd $(builddir); $(CXX) $(CXXFLAGS) -DEW_MADEBY=\"$(USER)\"  -DEW_OPT_LEVEL=\"$(optlevel)\" -DEW_COMPILER=\""$(shell which $(CXX))"\" -DEW_LIBDIR=\"${SW4LIB}\" -DEW_INCDIR=\"${SW4INC}\" -DEW_HOSTNAME=\""$(shell hostname)"\" -DEW_WHEN=\""$(shell date)"\" -c ../$<

# having version.o depend on .FORCE has the effect of always building version.o
.FORCE:

$(builddir)/%.o:src-sw4/quadpack/%.f
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$<

$(builddir)/%.o:src-sw4/%.f
	/bin/mkdir -p $(builddir)
	cd $(builddir); $(FC) $(FFLAGS) -c ../$<

$(builddir)/%.o:src-sw4/%.C
	/bin/mkdir -p $(builddir)
	 cd $(builddir); $(CXX) $(CXXFLAGS) -c ../$< 

$(builddir)/%.o:src-opt/%.C
	/bin/mkdir -p $(builddir)
	 cd $(builddir); $(CXX) $(CXXFLAGS) -c ../$< 

$(builddir)/RandomizedMaterial.o:src-sw4/RandomizedMaterial.C
	/bin/mkdir -p $(builddir)
	 cd $(builddir); $(CXX) $(CXXFLAGS) -std=c++11 -c ../$< 

clean:
	/bin/mkdir -p $(optdir)
	/bin/mkdir -p $(debugdir)
	cd $(optdir); /bin/rm -f sw4opt sw4mopt lamb1 *.o; cd ../$(debugdir); /bin/rm -f sw4opt sw4mopt lambd1 *.o


