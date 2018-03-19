#ifndef __POLICIES_H__
#define __POLICIES_H__
#include "RAJA/RAJA.hpp"
#ifdef ENABLE_CUDA
using RHS4_EXEC_POL = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<4>>, 
  RAJA::nested::For<2, RAJA::cuda_threadblock_z_exec<16>>> >;


using ICSTRESS_EXEC_POL = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<4>>> >;


using CONSINTP_EXEC_POL1 = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<32>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<32>>> >;


using CONSINTP_EXEC_POL3 = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<4>>, 
  RAJA::nested::For<2, RAJA::cuda_threadblock_z_exec<16>>> >;


using CONSINTP_EXEC_POL4 = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<4>>> >;


using CONSINTP_EXEC_POL5 = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<32>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<32>>> >;


using PRELIM_CORR_EXEC_POL1 = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<4>>> >;


using ENFORCEBC_CORR_EXEC_POL1 = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<4>>> >;

// Policy used in EW::get_exact_point_source
using GEPS_EXEC_POL = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<4>>, 
  RAJA::nested::For<2, RAJA::cuda_threadblock_z_exec<16>>> >;

typedef RAJA::cuda_exec<1024> SARRAY_LOOP_POL1 ;

using SARRAY_LOOP_POL2=
  RAJA::nested::Policy<
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_z_exec<1>>,
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<1>>,
  RAJA::nested::For<2, RAJA::cuda_threadblock_x_exec<1024>>>,
  RAJA::nested::For<3, RAJA::cuda_loop_exec> >;

typedef RAJA::cuda_exec<1024> PREDFORT_LOOP_POL;

typedef RAJA::cuda_exec<1024> CORRFORT_LOOP_POL;

typedef RAJA::cuda_exec<1024> DPDMTFORT_LOOP_POL;

using DPDMT_WIND_LOOP_POL=
  RAJA::nested::Policy<
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_z_exec<1>>,
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<1>>,
  RAJA::nested::For<2, RAJA::cuda_threadblock_x_exec<1024>>>,
  RAJA::nested::For<3, RAJA::cuda_loop_exec> >;

using COPY_KPLANE_EXEC_POL = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<16>>, 
  RAJA::nested::For<2, RAJA::cuda_threadblock_z_exec<16>>> >;

using ENERGY4CI_EXEC_POL = 
  RAJA::nested::Policy< 
  RAJA::nested::CudaCollapse<
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<16>>, 
  RAJA::nested::For<2, RAJA::cuda_threadblock_z_exec<16>>> >;


#define SYNC_DEVICE cudaDeviceSynchronize()

//****************************** OMP POLICIES *******************************************
#else

using RHS4_EXEC_POL =
  RAJA::nested::Policy< 
  RAJA::nested::OmpParallelCollapse< 
  RAJA::nested::For<0>,         
  RAJA::nested::For<1>,         
  RAJA::nested::For<2> > > ;

using ICSTRESS_EXEC_POL = RAJA::nested::Policy< 
  RAJA::nested::For<0, RAJA::parallel_exec>,
  RAJA::nested::For<1, RAJA::simd_exec> >;


using CONSINTP_EXEC_POL1 = RAJA::nested::Policy< 
  RAJA::nested::For<0, RAJA::parallel_exec>,
  RAJA::nested::For<1, RAJA::simd_exec> >;


using CONSINTP_EXEC_POL3 =
  RAJA::nested::Policy< 
  RAJA::nested::OmpParallelCollapse< 
  RAJA::nested::For<0>,         
  RAJA::nested::For<1>,         
  RAJA::nested::For<2> > > ;


using CONSINTP_EXEC_POL4 = RAJA::nested::Policy< 
  RAJA::nested::For<0, RAJA::parallel_exec>,
  RAJA::nested::For<1, RAJA::simd_exec> >;

using PRELIM_CORR_EXEC_POL1 = RAJA::nested::Policy< 
  RAJA::nested::For<0, RAJA::parallel_exec>,
  RAJA::nested::For<1, RAJA::simd_exec> >;

using ENFORCEBC_CORR_EXEC_POL1 = RAJA::nested::Policy< 
  RAJA::nested::For<0, RAJA::parallel_exec>,
  RAJA::nested::For<1, RAJA::parallel_exec> >;

// Policy used in EW::get_exact_point_source
using GEPS_EXEC_POL =
  RAJA::nested::Policy< 
  RAJA::nested::OmpParallelCollapse< 
  RAJA::nested::For<0>,         
  RAJA::nested::For<1>,         
  RAJA::nested::For<2> > > ;

typedef RAJA::omp_parallel_for_exec SARRAY_LOOP_POL1;

using SARRAY_LOOP_POL2 =
  RAJA::nested::Policy<
  RAJA::nested::For<0, RAJA::omp_parallel_for_exec>,
  RAJA::nested::For<1, RAJA::omp_parallel_for_exec>,
  RAJA::nested::For<2, RAJA::simd_exec>,
  RAJA::nested::For<3, RAJA::seq_exec>
  >;

typedef RAJA::omp_parallel_for_exec PREDFORT_LOOP_POL;

typedef RAJA::omp_parallel_for_exec CORRFORT_LOOP_POL;

typedef RAJA::omp_parallel_for_exec DPDMTFORT_LOOP_POL;


using DPDMT_WIND_LOOP_POL =
  RAJA::nested::Policy<
  RAJA::nested::For<0, RAJA::omp_parallel_for_exec>,
  RAJA::nested::For<1, RAJA::omp_parallel_for_exec>,
  RAJA::nested::For<2, RAJA::simd_exec>,
  RAJA::nested::For<3, RAJA::seq_exec>
  >;

using COPY_KPLANE_EXEC_POL =
  RAJA::nested::Policy< 
  RAJA::nested::OmpParallelCollapse< 
  RAJA::nested::For<0>,         
  RAJA::nested::For<1>,         
  RAJA::nested::For<2> > > ;

using ENERGY4CI_EXEC_POL = RAJA::nested::Policy< 
  RAJA::nested::For<0, RAJA::parallel_exec>,
  RAJA::nested::For<1, RAJA::parallel_exec>,
  RAJA::nested::For<2, RAJA::simd_exec> >;

#define SYNC_DEVICE

#endif

#endif
