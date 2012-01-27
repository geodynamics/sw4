ifndef debug
   optlevel = OPTIMIZE
else
   optlevel = DEBUG
endif

# for mac
FC = gfortran
CC = mpicc
CXX = mpicxx
linklibs = -L/home/sjogreen2/apps/lapack-3.2 -llapack -lblas -lgfortran

# for atlas
#FC = mpiifort
#CC = mpiicc
#CXX = mpiicpc
#linklibs = -llapack -lblas -lifcore


ifeq ($(optlevel),DEBUG)
#   FFLAGS    = -g -fno-underscoring
   FFLAGS    = -g
   CXXFLAGS  = -g 
   CFLAGS = -g
else
#   FFLAGS = -O3 -fno-underscoring
   FFLAGS = -O3
   CXXFLAGS = -O
   CFLAGS = -O
endif


OBJ  = main.o EW.o Sarray.o GeographicCoord.o version.o parseInputFile.o ForcingTwilight.o \
       curvilinearGrid.o boundaryOp.o bcfort.o twilightfort.o rhs4th3fort.o \
       parallelStuff.o Source.o MaterialProperty.o MaterialData.o material.o setupRun.o \
       solve.o solerr3.o Parallel_IO.o Image.o GridPointSource.o MaterialBlock.o

OBJ2 = test.o boundaryOp.o boundaryDisc.o gettimec.o 
OBJ3 = test3.o boundaryOp.o gettimec.o

sbp4f: $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $(OBJ) $(linklibs)

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
	/bin/rm -f ewave $(OBJ)
