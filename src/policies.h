#ifndef __POLICIES_H__
#define __POLICIES_H__
#include "RAJA/RAJA.hpp"
#ifdef ENABLE_CUDA

typedef RAJA::cuda_exec<1024> DEFAULT_LOOP1;


using DEFAULT_LOOP2 = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
  RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>, 
		       RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
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


using PRELIM_CORR_EXEC_POL1 =  ICSTRESS_EXEC_POL;
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

using COPY_KPLANE_EXEC_POL =  RHS4_EXEC_POL;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;

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

#define SYNC_DEVICE cudaDeviceSynchronize()

//****************************** OMP POLICIES *******************************************
#else

using RHS4_EXEC_POL =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1>,         
  RAJA::statement::For<2> > > ;

using ICSTRESS_EXEC_POL = RAJA::KernelPolicy< 
  RAJA::statement::For<0, RAJA::parallel_exec>,
  RAJA::statement::For<1, RAJA::simd_exec> >;


using CONSINTP_EXEC_POL1 = RAJA::KernelPolicy< 
  RAJA::statement::For<0, RAJA::parallel_exec>,
  RAJA::statement::For<1, RAJA::simd_exec> >;


using CONSINTP_EXEC_POL3 =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1>,         
  RAJA::statement::For<2> > > ;


using CONSINTP_EXEC_POL4 = RAJA::KernelPolicy< 
  RAJA::statement::For<0, RAJA::parallel_exec>,
  RAJA::statement::For<1, RAJA::simd_exec> >;

using PRELIM_CORR_EXEC_POL1 = RAJA::KernelPolicy< 
  RAJA::statement::For<0, RAJA::parallel_exec>,
  RAJA::statement::For<1, RAJA::simd_exec> >;

using ENFORCEBC_CORR_EXEC_POL1 = RAJA::KernelPolicy< 
  RAJA::statement::For<0, RAJA::parallel_exec>,
  RAJA::statement::For<1, RAJA::parallel_exec> >;

// Policy used in EW::get_exact_point_source
using GEPS_EXEC_POL =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1>,         
  RAJA::statement::For<2> > > ;

typedef RAJA::omp_parallel_for_exec SARRAY_LOOP_POL1;

using SARRAY_LOOP_POL2 =
  RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::omp_parallel_for_exec>,
  RAJA::statement::For<1, RAJA::omp_parallel_for_exec>,
  RAJA::statement::For<2, RAJA::simd_exec>,
  RAJA::statement::For<3, RAJA::seq_exec>
  >;

typedef RAJA::omp_parallel_for_exec PREDFORT_LOOP_POL;

typedef RAJA::omp_parallel_for_exec CORRFORT_LOOP_POL;

typedef RAJA::omp_parallel_for_exec DPDMTFORT_LOOP_POL;


using DPDMT_WIND_LOOP_POL =
  RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::omp_parallel_for_exec>,
  RAJA::statement::For<1, RAJA::omp_parallel_for_exec>,
  RAJA::statement::For<2, RAJA::simd_exec>,
  RAJA::statement::For<3, RAJA::seq_exec>
  >;

using COPY_KPLANE_EXEC_POL =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1>,         
  RAJA::statement::For<2> > > ;

using ENERGY4CI_EXEC_POL = RAJA::KernelPolicy< 
  RAJA::statement::For<0, RAJA::parallel_exec>,
  RAJA::statement::For<1, RAJA::parallel_exec>,
  RAJA::statement::For<2, RAJA::simd_exec> >;

using ODDIODDJ_EXEC_POL1 =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1> > > ;

using ODDIODDJ_EXEC_POL2 =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1>,         
  RAJA::statement::For<2> > > ;

using ODDIEVENJ_EXEC_POL1 =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1> > > ;

using ODDIEVENJ_EXEC_POL2 =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1>,         
  RAJA::statement::For<2> > > ;

using BCFORT_EXEC_POL1 =
  RAJA::KernelPolicy< 
  RAJA::statement::OmpParallelCollapse< 
  RAJA::statement::For<0>,         
  RAJA::statement::For<1>,         
  RAJA::statement::For<2> > > ;

using BCFORT_EXEC_POL2 = RAJA::KernelPolicy< 
  RAJA::statement::For<0, RAJA::parallel_exec>,
  RAJA::statement::For<1, RAJA::simd_exec> >;

#define SYNC_DEVICE

#endif

#endif
