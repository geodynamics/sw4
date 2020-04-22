.SUFFIXES: .f .o
.SUFFIXES: .f90 .o
.SUFFIXES: .f90 .mod

BINDIR = ./
FC = gfortran
#FC = gfortran-mp-7
#FC = mpif90

FFLAGS = -O#1 -g -fbacktrace -fcheck=bounds
#FINC   = -I/Users/siyangwang/libs/silo/include
LDFLAGS= -llapack -lblas
#-L/Users/siyangwang/libs/silo/lib -lsilo

LD = $(FC)

# Object files

OBJS = energy_4th_20200108.o problemsetup_new_3d.o SBP_operator.o
#OBJS = Elastic_3D_Interface_nonconforming_it.o problemsetup_new_3d.o SBP_operator.o Iterative_solver.o

#
.PHONY:copy

$(BINDIR)/energy.x: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(LDFLAGS)


.f.o:
	$(FC) -c $(FFLAGS) $(FINC) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FINC) $<
.f90.mod:
	$(FC) -c $(FFLAGS) $(FINC) $<


# The dependencies between modules

#Elastic_3D_Interface_nonconforming_it.f90: problemsetup_new_3d.mod SBP_operator.mod Iterative_solver.mod
energy_4th_20200108.f90: problemsetup_new_3d.mod SBP_operator.mod

clean:
	rm -rf *.x *.o *.a *.mod *.txt s0* f0* smultivar.0* fmultivar.0* *__genmod.f90 *~ *#
