ifndef debug
   optlevel = OPTIMIZE
else
   optlevel = DEBUG
endif

UNAME := $(shell uname)
HOSTNAME := $(shell hostname)

SW4INC    = $(SW4ROOT)/include
SW4LIB    = $(SW4ROOT)/lib

ifeq ($(optlevel),DEBUG)
   FFLAGS    = -g
   CXXFLAGS  = -g 
   CFLAGS    = -g
else
   FFLAGS   = -O3 
   CXXFLAGS = -O 
   CFLAGS   = -O 
endif

# for mac
ifeq ($(UNAME),Darwin)
  FC = openmpif77
  CC = openmpicc
  CXX = openmpicxx
  linklibs = -framework vecLib -L/opt/local/lib/gcc45 -lgfortran

  ifeq ($(optlevel),DEBUG)
    FFLAGS = -g -fno-underscoring
  else
    FFLAGS = -O3 -fno-underscoring
  endif
endif

# for linux
ifeq ($(UNAME),Linux)
  FC = mpiifort
  CC = mpiicc
  CXX = mpiicpc
  linklibs =  -llapack -lblas -lifcore
# for Bjorn's tux box
  ifeq ($(findstring tux,$(HOSTNAME)),tux)
     FC = gfortran
     CC = mpicc
     CXX = mpicxx
    linklibs = -L/home/sjogreen2/apps/lapack-3.2 -llapack -lblas -lgfortran
  endif
endif

ifeq ($(etree),yes)
   CXXFLAGS += -DENABLE_ETREE
   CXXFLAGS += -I$(SW4INC)
   linklibs += -Wl,-rpath=$(SW4LIB) -L$(SW4LIB) -lcencalvm -lproj
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
       lamb_exact_numquad.o twilightsgfort.o EtreeFile.o MaterialIfile.o

OBJOPT  = optmain.o EW.o Sarray.o version.o parseInputFile.o ForcingTwilight.o \
       curvilinearGrid.o boundaryOp.o bcfort.o twilightfort.o rhs4th3fort.o \
       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
       solve.o solerr3.o Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o testsrc.o \
       TimeSeries.o sacsubc.o SuperGrid.o addsgd.o linsolvelu.o velsum.o solve-backward.o \
       rayleighfort.o energy4.o TestRayleighWave.o MaterialPfile.o \
       Filter.o Polynomial.o SecondOrderSection.o time_functions.o Qspline.o \
       lamb_exact.o twilightsgfort.o EtreeFile.o MaterialIfile.o

OBJ2 = test.o boundaryOp.o boundaryDisc.o gettimec.o 
OBJ3 = test3.o boundaryOp.o gettimec.o
OBJ4 = filterTest.o Filter.o SecondOrderSection.o Polynomial.o
sw4: $(OBJ) $(QUADPACK)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(QUADPACK) $(linklibs)

sw4opt: $(OBJOPT) $(QUADPACK)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJOPT) $(QUADPACK) $(linklibs)

version.o:version.C
	$(CXX) $(CXXFLAGS) -DEW_MADEBY=\"$(USER)\"  -DEW_OPT_LEVEL=\"$(optlevel)\" -DEW_COMPILER=\""$(shell which $(CXX))"\" -DEW_BASEDIR=\"NotUsed\" -DEW_HOSTNAME=\""$(shell hostname)"\" -DEW_WHEN=\""$(shell date)"\" -c $<

test3: $(OBJ3)
	$(FC) -O3 -Wall -o $@ $(OBJ3)

test: $(OBJ2)
	$(FC) -O3 -Wall -o $@ $(OBJ2)

filterTest: $(OBJ4)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ4)

%.o:%.f
	$(FC) $(FFLAGS) -c $<

%.o:%.c
	$(CC) $(CFLAGS) -c $<

%.o:%.C
	$(CXX) $(CXXFLAGS) -c $<

%.o:quadpack/%.f
	$(FC) $(FFLAGS) -c $<

clean:
	/bin/rm -f sw4 sw4opt $(OBJ) $(OBJOPT) $(QUADPACK)
