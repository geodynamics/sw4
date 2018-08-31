proj = yes
ckernel=yes 
openmp=no
raja_cuda=yes
FC = xlf
LINKER = mpixlC
CXX = nvcc
RAJA_LOCATION = /usr/workspace/wsrzd/ramesh/Project6/SW4/RAJA/install_cuda9.2.148
PROJ_HOME = /usr/workspace/wsrzd/ramesh/Project6/SW4/proj-5.0.0/install170818
CALIPER_LOCATION = /usr/workspace/wsrzd/ramesh/Project6/SW4/Caliper/install170818
UMPIRE_HOME = /usr/workspace/wsrzd/ramesh/Project6/SW4/Umpire/Umpire/install
EXTRA_FORT_FLAGS = -qstrict
LINKFLAGS =  -arch=sm_70
EXTRA_CXX_FLAGS =  -O3 -use_fast_math -DDISABLE_PREFETCH  -DSW4_USE_UMPIRE2 -ccbin mpixlC -Xcompiler="-qmaxmem=-1 -qstrict " -std=c++11 --expt-extended-lambda -restrict -arch=sm_70 -DENABLE_CALIPER2 -I $(CUDA_HOME)/include -I$(RAJA_LOCATION)/include  -I$(CALIPER_LOCATION)/include --x cu -DUSE_NVTX -DRAJA_USE_CUDA -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -DCUDA_CODE -DENABLE_CUDA -dc -I$(PROJ_HOME)/include -I$(UMPIRE_HOME)/include 

EXTRA_CXX_FLAGS2 =  -O3 -use_fast_math -DDISABLE_PREFETCH -DSW4_USE_UMPIRE2 -ccbin mpixlC -Xcompiler="-qmaxmem=-1 -qstrict " -std=c++11 --expt-extended-lambda -restrict -arch=sm_70 -DENABLE_CALIPER2 -I $(CUDA_HOME)/include -I$(RAJA_LOCATION)/include  --x cu  -DUSE_NVTX -DRAJA_USE_CUDA -DSW4_CROUTINES -DRAJA_USE_RESTRICT_PTR -DCUDA_CODE -DENABLE_CUDA -I$(UMPIRE_HOME)/include 

EXTRA_LINK_FLAGS =  -dc -L /usr/tcetmp/packages/xl/xl-beta-2018.08.08/xlf/16.1.1/lib -lxlf90 -lxlfmath -L /usr/tcetmp/packages/blas/blas-3.6.0-xlf-15.1.5/lib -L /usr/tcetmp/packages/lapack/lapack-3.6.0-xlf-15.1.5/lib/ -llapack -lblas  -L $(CUDA_HOME)/lib64 -lcudart -lnvToolsExt -lcuda -L $(RAJA_LOCATION)/lib -lRAJA  -lcudadevrt -Wl,-rpath=$(CALIPER_LOCATION)/lib64 -L $(CALIPER_LOCATION)/lib64 -lcaliper -L $(PROJ_HOME)/lib64 -lproj -L /usr/lib64/nvidia -lnvidia-ml  -L$(UMPIRE_HOME)/lib -lumpire -lumpire_util -lumpire_tpl_judy -lumpire_resource -lumpire_strategy -lumpire_resource -lumpire_op
