proj = yes
FC = xlf
LINKER = mpixlC
CXX = nvcc
RAJA_LOCATION = /usr/workspace/wsb/deg/ramesh/RAJA_053/RAJA/install_ray053
RAJA_LOCATION = /usr/workspace/wsb/deg/ramesh/RAJA_01_24_2018/RAJA/install
RAJA_LOCATION = /usr/workspace/wsb/deg/ramesh/RAJA_053/RAJA/install_053NOMP
RAJA_LOCATION = /usr/workspace/wsb/deg/ramesh/RAJA_01_24_2018/RAJA/install_nomp
RAJA_LOCATION = /usr/workspace/wsb/ramesh/RAJA/2018/RAJA/install-sierra-nomp-9.2
PROJ_HOME = /usr/workspace/wsb/ramesh/Project6/Proj/proj-5.0.0/install/
CALIPER_LOCATION = /usr/workspace/wsb/ramesh/Caliper/Caliper/install_ray_9.2
EXTRA_FORT_FLAGS = 
LINKFLAGS =  -arch=sm_70
EXTRA_CXX_FLAGS =  -O3 -use_fast_math -ccbin mpixlC -Xcompiler="-qmaxmem=-1 -qstrict -qsmp=omp" -std=c++11 --expt-extended-lambda -restrict -arch=sm_70 -DENABLE_CALIPER -I $(CUDA_HOME)/include -I$(RAJA_LOCATION)/include  -I$(CALIPER_LOCATION)/include --x cu -DUSE_NVTX -DRAJA_USE_CUDA -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -DCUDA_CODE -DENABLE_CUDA -dc -I$(PROJ_HOME)/include

EXTRA_CXX_FLAGS2 =  -O3 -use_fast_math -ccbin mpixlC -Xcompiler="-qmaxmem=-1 -qstrict -qsmp=omp " -std=c++11 --expt-extended-lambda -restrict -arch=sm_70 -DENABLE_CALIPER -I $(CUDA_HOME)/include -I$(RAJA_LOCATION)/include  --x cu  -DUSE_NVTX -DRAJA_USE_CUDA -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -DCUDA_CODE -DENABLE_CUDA 

EXTRA_LINK_FLAGS = -qsmp=omp -dc -L /usr/tcetmp/packages/xl/xl-beta-2017.03.28/xlf/16.1.0/lib/ -L /usr/tcetmp/packages/blas/blas-3.6.0-xlf-15.1.5/lib -L /usr/tcetmp/packages/lapack/lapack-3.6.0-xlf-15.1.5/lib/ -lxlf90 -llapack -lblas  -L $(CUDA_HOME)/lib64 -lcudart -lnvToolsExt -lcuda -L $(RAJA_LOCATION)/lib -lRAJA  -lcudadevrt -Wl,-rpath=$(CALIPER_LOCATION)/lib64 -L $(CALIPER_LOCATION)/lib64 -lcaliper -L $(PROJ_HOME)/lib -lproj
