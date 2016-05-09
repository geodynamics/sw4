#clang, use clang
#FC = clang
#CC = clang
#CXX = clang++
#FFLAGS = -O2
#CFLAGS = -O2

#pgi, use pgi-15.7
# 
#FC = pgf77
#CC = pgcc
#CXX =  pgCC
#FFLAGS = -O3
#CFLAGS = -O3

#intel,  use ic-16.0.150
#
FC = ifort
CC = icc
CXX= icpc
FFLAGS = -O2 -g
CFLAGS = -O2 -g 
CXXFLAGS = -O2 -g

#gnu, use gcc-4.9.3p
#FC = gfortran 
#CC = gcc 
#CXX= g++
#FFLAGS = -O2 -ftree-vectorize -ftree-vectorizer-verbose=6
#CFLAGS = -O3 -ftree-vectorize -ftree-vectorizer-verbose=6

COBJ = rhs4sg.o rhs4sg_rev.o gettimec.o
FOBJ = rhs4th3fortsgstr.o boundaryOp.o

main : $(FOBJ) $(COBJ) testil.o
	$(CXX) -O2 -o $@ testil.o $(FOBJ) $(COBJ)

