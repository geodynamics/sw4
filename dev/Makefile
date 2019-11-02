.SUFFIXES: .f .o
.SUFFIXES: .f90 .o
.SUFFIXES: .f90 .mod
.SUFFIXES: .C .o
BINDIR = ./
FC = gfortran
CXX = g++
#FC = gfortran-mp-7
#FC = mpif90

FFLAGS = -O1 -g -fbacktrace -fcheck=bounds
#FINC   = -I/Users/siyangwang/libs/silo/include
LDFLAGS= -llapack -lblas -lgfortran
#-L/Users/siyangwang/libs/silo/lib -lsilo
CXXFLAGS = -std=c++11 -g -O3
LD = $(CXX)

# Object files

OBJS = main.o fortran_functions_current.o Sarray.o cfunctions.o

#
.PHONY:copy

$(BINDIR)/x: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(LDFLAGS)


.f.o:
	$(FC) -c $(FFLAGS) $(FINC) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FINC) $<
.f90.mod:
	$(FC) -c $(FFLAGS) $(FINC) $<
.C.o:
	$(CXX) -c $(CXXFLAGS) $(CINC) $<

%.o: %.mod

# The dependencies between modules

x: Farray.h Sarray.h cfunctions.h
main.o: Farray.h Sarray.h cfunctions.h
cfunctions.o: Farray.h Sarray.h cfunctions.h
clean:
	rm -rf *.x *.txt *.o *.a *.mod s0* f0* smultivar.0* fmultivar.0* *__genmod.f90 *~ *#

format:
	clang-format -style Google -i *.C
	clang-format -style Google -i *.h
