ifndef debug
   optlevel = OPTIMIZE
else
   optlevel = DEBUG
endif

UNAME := $(shell uname)
HOSTNAME := $(shell hostname)

ifeq ($(optlevel),DEBUG)
   FFLAGS    = -g
   CXXFLAGS  = -g 
   CFLAGS = -g
else
   FFLAGS = -O3
   CXXFLAGS = -O
   CFLAGS = -O
endif

# for mac
ifeq ($(UNAME),Darwin)
  FC = gfortran
  CC = mpicc
  CXX = mpicxx
  linklibs = -framework vecLib -lgfortran

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
  linklibs = -llapack -lblas -lifcore
# for Bjorn's tux box
  ifeq ($(findstring tux,$(HOSTNAME)),tux)
     FC = gfortran
     CC = mpicc
     CXX = mpicxx
    linklibs = -L/home/sjogreen2/apps/lapack-3.2 -llapack -lblas -lgfortran
  endif
endif


OBJ  = main.o EW.o Sarray.o version.o parseInputFile.o ForcingTwilight.o \
       curvilinearGrid.o boundaryOp.o bcfort.o twilightfort.o rhs4th3fort.o \
       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
       solve.o solerr3.o Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o testsrc.o \
       TimeSeries.o sacsubc.o SuperGrid.o addsgd.o velsum.o rayleighfort.o TestRayleighWave.o

OBJOPT  = optmain.o EW.o Sarray.o version.o parseInputFile.o ForcingTwilight.o \
       curvilinearGrid.o boundaryOp.o bcfort.o twilightfort.o rhs4th3fort.o \
       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
       solve.o solerr3.o Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o testsrc.o \
       TimeSeries.o sacsubc.o SuperGrid.o addsgd.o linsolvelu.o velsum.o solve-backward.o  rayleighfort.o

OBJ2 = test.o boundaryOp.o boundaryDisc.o gettimec.o 
OBJ3 = test3.o boundaryOp.o gettimec.o

s4f: $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $(OBJ) $(linklibs)

sbp4opt: $(OBJOPT)
	$(CXX) $(CFLAGS) -o $@ $(OBJOPT) $(linklibs)

version.o:version.C
	$(CXX) $(CXXFLAGS) -DEW_MADEBY=\"$(USER)\"  -DEW_OPT_LEVEL=\"$(optlevel)\" -DEW_COMPILER=\""$(shell which $(CXX))"\" -DEW_BASEDIR=\"NotUsed\" -DEW_HOSTNAME=\""$(shell hostname)"\" -DEW_WHEN=\""$(shell date)"\" -c $<

test3: $(OBJ3)
	gfortran -O3 -Wall -o $@ $(OBJ3)

test: $(OBJ2)
	gfortran -O3 -Wall -o $@ $(OBJ2)

%.o:%.f
	$(FC) $(FFLAGS) -c $<

%.o:%.c
	$(CC) $(CFLAGS) -c $<

%.o:%.C
	$(CXX) $(CFLAGS) -c $<

clean:
	/bin/rm -f s4f sbp4opt $(OBJ)
