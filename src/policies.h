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
  RAJA::nested::For<0, RAJA::cuda_threadblock_x_exec<4>>, 
  RAJA::nested::For<1, RAJA::cuda_threadblock_y_exec<4>>> >;


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


#define SYNC_DEVICE cudaDeviceSynchronize()
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

#define SYNC_DEVICE

#endif

#endif
