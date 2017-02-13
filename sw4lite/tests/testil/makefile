#clang, use clang
#use clang-3.9.0
#FC = clang
#CC = clang
#CXX = clang++
#FFLAGS = -O2
#CFLAGS = -O2
#CXXFLAGS = -O2
#pgi, use pgi-15.7
#pgi, use pgi-16.7
# 
#FC = pgf77
#CC = pgcc
#CXX =  pgc++
#FFLAGS = -O3 
#CFLAGS = -O3 
#CXXFLAGS = -O3
#OMPOPT = -mp
#intel,  use ic-16.0.150
# newer compiler, ic-16.0.210
FC = ifort
CC = icc
CXX= icpc
FFLAGS =  -O2 -warn all 
CFLAGS =  -O2 -Wall 
CXXFLAGS =-O2 -Wall 
OMPOPT = -qopenmp
#FFLAGS =  -g -warn all
#CFLAGS =  -g -warn all
#CXXFLAGS =-g -warn all

#gnu, use gcc-4.9.3p
#FC = gfortran 
#CC = gcc 
#CXX= g++
#OMPOPT = -fopenmp
#FFLAGS = -g -Wall
#CFLAGS = -g -Wall
#FFLAGS = -O2 -ftree-vectorize -ftree-vectorizer-verbose=6 -Wall
#CXXFLAGS = -O2 -ftree-vectorize -ftree-vectorizer-verbose=6 -Wall
#CFLAGS   = -O2 -ftree-vectorize -ftree-vectorizer-verbose=6 -Wall

#IBM
#FC = xlf
#CC = xlc
#CXX= xlC
#FFLAGS = -O2 
#CFLAGS = -O2 
#CXXFLAGS = -O3 

ifeq ($(openmp),yes)
   FFLAGS += $(OMPOPT)
   CXXFLAGS += $(OMPOPT)
   CFLAGS += $(OMPOPT)
endif

COBJ = rhs4sg.o rhs4sg_rev.o gettimec.o grid-utilities.o rhs4sgcurv.o rhs4sgcurv_rev.o metricC.o \
       addsgd4fort.o addsgd4fort_rev.o
FOBJ = rhs4th3fortsgstr.o boundaryOp.o curvilinear4sg.o metric.o addsgd4.o

main : $(FOBJ) $(COBJ) testil.o
	$(CXX) $(CXXFLAGS) -o $@ testil.o $(FOBJ) $(COBJ)

