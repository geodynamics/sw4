ifeq ($(debug),yes)
   CXXFLAGS += -g
else
   CXXFLAGS += -O2 
endif

ifneq ($(hdf5),no)
#   HDF5ROOT   = /opt/hdf5-1.8/usr
#   HDF5ROOT   = /g/g91/bjorn/local
   HDF5ROOT   = /usr/local/Cellar/hdf5/1.10.2_1
   CXXFLAGS  += -I$(HDF5ROOT)/include -DIOTEST_HDF5
   LINKFLAGS += -L$(HDF5ROOT)/lib
   linklibs  = -lhdf5_hl -lhdf5
endif

CXX = mpicxx

IDIR = $(CURDIR)/

OBJ = iotest.o Parallel_IO-sw4.o 

iotest: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LINKFLAGS) $(linklibs)

gentestfile: gentestfile.o
	$(CXX) $(CXXFLAGS) -o $@ gentestfile.o $(LINKFLAGS) $(linklibs)

%.o:%.C
	$(CXX) $(CXXFLAGS) -c $<

clean:
	/bin/rm -f *.o; /bin/rm -f iotest;/bin/rm -f gentestfile
