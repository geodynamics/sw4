proj = yes
FC = ifort
LINKER = mpicxx
CXX = mpicxx
RAJA_LOCATION = /usr/workspace/wsb/deg/ramesh/RAJA_053/RAJA/install_ray053
RAJA_LOCATION = /usr/workspace/wsb/deg/ramesh/RAJA_01_24_2018/RAJA/install
RAJA_LOCATION = /usr/workspace/wsb/deg/ramesh/RAJA_053/RAJA/install_053NOMP
RAJA_LOCATION = /usr/workspace/wsb/deg/ramesh/RAJA_01_24_2018/RAJA/install_nomp
RAJA_LOCATION = /usr/workspace/wsb/ramesh/RAJA/2018/RAJA/install-ray-nomp
RAJA_LOCATION = /usr/workspace/wsb/ramesh/RAJA/2018/RAJA/install-quartz
PROJ_HOME = /usr/workspace/wsb/ramesh/Project6/Proj/proj-5.0.0/install-quartz
CALIPER_LOCATION = /usr/workspace/wsb/ramesh/Caliper/Caliper/install_quartz
EXTRA_FORT_FLAGS = 
LINKFLAGS = 
EXTRA_CXX_FLAGS =  -O3 -std=c++11 -I$(RAJA_LOCATION)/include  -I$(CALIPER_LOCATION)/include -DENABLE_CALIPER -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -I$(PROJ_HOME)/include
EXTRA_CXX_FLAGS2 =  -O3 -std=c++11 -I$(RAJA_LOCATION)/include  -I$(CALIPER_LOCATION)/include -DENABLE_CALIPER -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -I$(PROJ_HOME)/include

MKL_PATH = /usr/tce/packages/mkl/mkl-2018.0/lib
EXTRA_LINK_FLAGS = -fopenmp -Wl,-rpath=$(SW4ROOT)/lib -Wl,-rpath=${MKL_PATH} -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -lifcore -Wl,-rpath=$(CALIPER_LOCATION)/lib64 -L $(CALIPER_LOCATION)/lib64 -lcaliper -L $(PROJ_HOME)/lib64 -lproj -Wl,-rpath=$(PROJ_HOME)/lib64

