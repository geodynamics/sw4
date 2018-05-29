#ifndef __POLICIES_H__
#define __POLICIES_H__
#include "RAJA/RAJA.hpp"
#ifdef ENABLE_CUDA
using REDUCTION_POLICY = RAJA::cuda_reduce<1024>;

typedef RAJA::cuda_exec<1024> DEFAULT_LOOP1;
#define SW4_FORCEINLINE __forceinline__
using XRHS_POL = 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernel<
       RAJA::statement::For<0, RAJA::cuda_block_exec, 
			    RAJA::statement::For<1, RAJA::cuda_block_exec, 
						 RAJA::statement::For<2, RAJA::cuda_thread_exec,
								      RAJA::statement::Lambda<0> >>>>>;

using DEFAULT_LOOP2 = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
  RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>, 
		       RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
					    RAJA::statement::Lambda<0> >>>>;

using DEFAULT_LOOP2X = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
  RAJA::statement::For<1, RAJA::cuda_block_exec, 
  RAJA::statement::For<0, RAJA::cuda_thread_exec,
  RAJA::statement::Lambda<0> >>>>;


using DEFAULT_LOOP3 = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16 >,
								   RAJA::statement::Lambda<0> >>>>>;
using DEFAULT_LOOP4 =
  RAJA::KernelPolicy<
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<1>,
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<1>,
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<1024>,
  RAJA::statement::For<3, RAJA::seq_exec,
								   RAJA::statement::Lambda<0> >>>>>>;
using RHS4_EXEC_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16 >,
								   RAJA::statement::Lambda<0> >>>>>;
			 
			 
using ICSTRESS_EXEC_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
  RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>, 
		       RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
					    RAJA::statement::Lambda<0> >>>>;


using CONSINTP_EXEC_POL1 = ICSTRESS_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<32>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<32>>> >;


using CONSINTP_EXEC_POL3 = RHS4_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;


using CONSINTP_EXEC_POL4  = ICSTRESS_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>>> >;


using CONSINTP_EXEC_POL5 =  ICSTRESS_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<32>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<32>>> >;


using PRELIM_CORR_EXEC_POL1 =  DEFAULT_LOOP2X;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>>> >;

using PRELIM_PRED_EXEC_POL1 =  ICSTRESS_EXEC_POL;
			 
using ENFORCEBC_CORR_EXEC_POL1 =  ICSTRESS_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>>> >;

// Policy used in EW::get_exact_point_source
using GEPS_EXEC_POL = RHS4_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;

typedef RAJA::cuda_exec<1024> SARRAY_LOOP_POL1 ;

using SARRAY_LOOP_POL2=
  RAJA::KernelPolicy<
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<1>,
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<1>,
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<1024>,
								   RAJA::statement::For<3, RAJA::seq_exec,
								   RAJA::statement::Lambda<0> >>>>>>;

typedef RAJA::cuda_exec<1024> PREDFORT_LOOP_POL;

typedef RAJA::cuda_exec<1024> CORRFORT_LOOP_POL;

typedef RAJA::cuda_exec<1024> DPDMTFORT_LOOP_POL;

using DPDMT_WIND_LOOP_POL=  SARRAY_LOOP_POL2;
  // RAJA::KernelPolicy<
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<1>>,
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<1>>,
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<1024>>>,
  // RAJA::statement::For<3, RAJA::cuda_loop_exec> >;

//using COPY_KPLANE_EXEC_POL =  RHS4_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;
using COPY_KPLANE_EXEC_POL  = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;



using ENERGY4CI_EXEC_POL =  RHS4_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;

// Increasing the 16,16 block size below can cause failures with the error message:
//CUDAassert: too many resources requested for launch
// /usr/workspace/wsb/deg/ramesh/RAJA_01_24_2018/RAJA/install_nomp/include/RAJA/policy/cuda/MemUtils_CUDA.hpp 216


using ODDIODDJ_EXEC_POL1 = ICSTRESS_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>> >;

using ODDIODDJ_EXEC_POL2 = RHS4_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;

using ODDIEVENJ_EXEC_POL1 =  ICSTRESS_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>> >;

using ODDIEVENJ_EXEC_POL2 =  RHS4_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;

using BCFORT_EXEC_POL1 =  RHS4_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;


using BCFORT_EXEC_POL2 = ICSTRESS_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>>> >;


using EVENIODDJ_EXEC_POL = 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernel<
  RAJA::statement::For<1, RAJA::cuda_block_exec, 
  RAJA::statement::For<0, RAJA::cuda_thread_exec,
  RAJA::statement::Lambda<0> >>>>;

#define SYNC_DEVICE SW4_CheckDeviceError(cudaDeviceSynchronize())
#define SYNC_STREAM SW4_CheckDeviceError(cudaStreamSynchronize(0))




//*****************************************************************************************
//****************************** OMP POLICIES *******************************************






#else

#define SW4_FORCEINLINE


using REDUCTION_POLICY = RAJA::omp_reduce;

using DEFAULT_LOOP1 = RAJA::omp_parallel_for_exec;


using DEFAULT_LOOP2 = 
    RAJA::KernelPolicy<
      RAJA::statement::For<1, RAJA::omp_parallel_for_exec,
  RAJA::statement::For<0, RAJA::omp_parallel_for_exec,           
          RAJA::statement::Lambda<0>
        > 
      > 
    >;


using DEFAULT_LOOP3 = 
    RAJA::KernelPolicy<
  RAJA::statement::For<2, RAJA::omp_parallel_for_exec, 
      RAJA::statement::For<1, RAJA::omp_parallel_for_exec, 
        RAJA::statement::For<0, RAJA::seq_exec,           
          RAJA::statement::Lambda<0>
  >
        > 
      > 
    >;

using DEFAULT_LOOP4 = 
    RAJA::KernelPolicy<
  RAJA::statement::For<3, RAJA::omp_parallel_for_exec, 
  RAJA::statement::For<2, RAJA::omp_parallel_for_exec, 
      RAJA::statement::For<1, RAJA::omp_parallel_for_exec, 
        RAJA::statement::For<0, RAJA::seq_exec,           
          RAJA::statement::Lambda<0>
  >
  >
        > 
      > 
    >;


using XRHS_POL = DEFAULT_LOOP3;

using RHS4_EXEC_POL = DEFAULT_LOOP3;


using ICSTRESS_EXEC_POL = DEFAULT_LOOP2;
 


using CONSINTP_EXEC_POL1 = DEFAULT_LOOP2;


using CONSINTP_EXEC_POL3 = DEFAULT_LOOP3;
  


using CONSINTP_EXEC_POL4 = DEFAULT_LOOP2;


using PRELIM_CORR_EXEC_POL1 = DEFAULT_LOOP2;

using PRELIM_PRED_EXEC_POL1 = DEFAULT_LOOP2;

using ENFORCEBC_CORR_EXEC_POL1 = DEFAULT_LOOP2;


// Policy used in EW::get_exact_point_source
using GEPS_EXEC_POL = DEFAULT_LOOP3;
 

typedef RAJA::omp_parallel_for_exec SARRAY_LOOP_POL1;

using SARRAY_LOOP_POL2 = DEFAULT_LOOP4;
 

typedef RAJA::omp_parallel_for_exec PREDFORT_LOOP_POL;

typedef RAJA::omp_parallel_for_exec CORRFORT_LOOP_POL;

typedef RAJA::omp_parallel_for_exec DPDMTFORT_LOOP_POL;


using DPDMT_WIND_LOOP_POL = DEFAULT_LOOP4;

using COPY_KPLANE_EXEC_POL = DEFAULT_LOOP3;
 
using ENERGY4CI_EXEC_POL = DEFAULT_LOOP3;


using ODDIODDJ_EXEC_POL1 = DEFAULT_LOOP2;
 

using ODDIODDJ_EXEC_POL2 = DEFAULT_LOOP3;
 

using ODDIEVENJ_EXEC_POL1 = DEFAULT_LOOP2;
 
using ODDIEVENJ_EXEC_POL2 = DEFAULT_LOOP3;


using BCFORT_EXEC_POL1 = DEFAULT_LOOP3;


using BCFORT_EXEC_POL2 = DEFAULT_LOOP2;

using EVENIODDJ_EXEC_POL = DEFAULT_LOOP2;

#define SYNC_DEVICE
#define SYNC_STREAM

#endif

#endif
