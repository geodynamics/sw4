#ifndef __POLICIES_H__
#define __POLICIES_H__
#include "RAJA/RAJA.hpp"
#ifdef ENABLE_CUDA

#define SW4_FORCEINLINE __forceinline__
#define SYNC_DEVICE SW4_CheckDeviceError(cudaDeviceSynchronize())
#define SYNC_STREAM SW4_CheckDeviceError(cudaStreamSynchronize(0))
#define SW4_PEEK SW4_CheckDeviceError(cudaPeekAtLastError());
//   SW4_CheckDeviceError(cudaStreamSynchronize(0));
typedef RAJA::cuda_exec<1024> DEFAULT_LOOP1;
typedef RAJA::cuda_exec<1024,true> DEFAULT_LOOP1_ASYNC;
using REDUCTION_POLICY = RAJA::cuda_reduce;

typedef RAJA::cuda_exec<1024> PREDFORT_LOOP_POL;
typedef RAJA::cuda_exec<512,true> PREDFORT_LOOP_POL_ASYNC;

typedef RAJA::cuda_exec<1024> CORRFORT_LOOP_POL;
typedef RAJA::cuda_exec<512,true> CORRFORT_LOOP_POL_ASYNC;

typedef RAJA::cuda_exec<1024> DPDMTFORT_LOOP_POL;


typedef RAJA::cuda_exec<256,true> DPDMTFORT_LOOP_POL_ASYNC;

typedef RAJA::cuda_exec<1024> SARRAY_LOOP_POL1 ;

#if SW4_RAJA_VERSION==6




using XRHS_POL = 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernel<
       RAJA::statement::For<0, RAJA::cuda_block_exec, 
			    RAJA::statement::For<1, RAJA::cuda_block_exec, 
						 RAJA::statement::For<2, RAJA::cuda_thread_exec,
								      RAJA::statement::Lambda<0> >>>>>;

using XRHS_POL_ASYNC = 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernelAsync<
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

using DEFAULT_LOOP2X_ASYNC = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
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

using RHS4_EXEC_POL_ASYNC = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
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

using ICSTRESS_EXEC_POL_ASYNC = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
  RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>, 
  RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
					    RAJA::statement::Lambda<0> >>>>;


using CONSINTP_EXEC_POL1 = 
RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
  RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>, 
  RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
					    RAJA::statement::Lambda<0> >>>>;
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
using PRELIM_CORR_EXEC_POL1_ASYNC =  DEFAULT_LOOP2X_ASYNC;

  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>>> >;

using PRELIM_PRED_EXEC_POL1 =  ICSTRESS_EXEC_POL;
using PRELIM_PRED_EXEC_POL1_ASYNC =  ICSTRESS_EXEC_POL_ASYNC;
			 
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




using SARRAY_LOOP_POL2=
  RAJA::KernelPolicy<
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<1>,
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<1>,
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<1024>,
								   RAJA::statement::For<3, RAJA::seq_exec,
								   RAJA::statement::Lambda<0> >>>>>>;








using DPDMT_WIND_LOOP_POL=  SARRAY_LOOP_POL2;
using DPDMT_WIND_LOOP_POL_ASYNC = 
RAJA::KernelPolicy<
  RAJA::statement::CudaKernelAsync<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<1>,
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<1>,
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<1024>,
								   RAJA::statement::For<3, RAJA::seq_exec,
								   RAJA::statement::Lambda<0> >>>>>>;

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
// using COPY_KPLANE_EXEC_POL  =
//   RAJA::KernelPolicy<
//   RAJA::statement::CudaKernel<
//     RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>,
// 			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>,
// 					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
// 								   RAJA::statement::Lambda<0> >>>>>;

using COPY_KPLANE_EXEC_POL =
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_block_exec, 
			 RAJA::statement::For<1, RAJA::cuda_block_exec, 
					      RAJA::statement::For<2, RAJA::cuda_thread_exec,
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


using ODDIODDJ_EXEC_POL1_ASYNC = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
  RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>, 
  RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
					    RAJA::statement::Lambda<0> >>>>;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>> >;

using ODDIODDJ_EXEC_POL2 = RHS4_EXEC_POL_ASYNC;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16>>> >;

using ODDIEVENJ_EXEC_POL1 =  ICSTRESS_EXEC_POL_ASYNC;
  // RAJA::KernelPolicy< 
  // RAJA::statement::CudaCollapse<
  // RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>>, 
  // RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>>> >;

using ODDIEVENJ_EXEC_POL2 =  RHS4_EXEC_POL_ASYNC;
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


using EVENIODDJ_EXEC_POL_ASYNC= 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernelAsync<
  RAJA::statement::For<1, RAJA::cuda_block_exec, 
  RAJA::statement::For<0, RAJA::cuda_thread_exec,
  RAJA::statement::Lambda<0> >>>>;

using EVENIEVENJ_EXEC_POL =  EVENIODDJ_EXEC_POL_ASYNC;



using TWILIGHTSG_POL = RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;

#elif SW4_RAJA_VERSION==7



using DEFAULT_LOOP2X = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
  RAJA::statement::For<1, RAJA::cuda_block_x_loop, 
  RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
  RAJA::statement::Lambda<0> >>>>;

using DEFAULT_LOOP2X_ASYNC = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
  RAJA::statement::For<1, RAJA::cuda_block_x_loop, 
  RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
  RAJA::statement::Lambda<0> >>>>;


using DEFAULT_LOOP3 = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
			      RAJA::statement::Tile<2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
          RAJA::statement::For<0, RAJA::cuda_thread_z_direct,
            RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
				 RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
						      RAJA::statement::Lambda<0> >>>>>>>>;

using SARRAY_LOOP_POL2 = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<1>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<1, RAJA::statement::tile_fixed<1>, RAJA::cuda_block_x_loop,
			      RAJA::statement::Tile<2, RAJA::statement::tile_fixed<1024>, RAJA::cuda_block_z_loop,
          RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
            RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
				 RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
						      RAJA::statement::For<3, RAJA::seq_exec,
									   RAJA::statement::Lambda<0> >>>>>>>>>;



// using RHS4_EXEC_POL = 
//   RAJA::KernelPolicy< 
//   RAJA::statement::CudaKernel<
//     RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
// 			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
// 					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16 >,
// 								   RAJA::statement::Lambda<0> >>>>>;

using ICSTRESS_EXEC_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::Tile<0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
			  RAJA::statement::Tile<1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
						RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
								     RAJA::statement::For<1, RAJA::cuda_thread_y_loop,		 
											  RAJA::statement::Lambda<0> >>>>>>;

using ICSTRESS_EXEC_POL_ASYNC = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
    RAJA::statement::Tile<0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
			  RAJA::statement::Tile<1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
						RAJA::statement::For<0, RAJA::cuda_thread_x_direct,
								     RAJA::statement::For<1, RAJA::cuda_thread_y_direct,		 
											  RAJA::statement::Lambda<0> >>>>>>;

// using RHS4_EXEC_POL_ASYNC = 
//   RAJA::KernelPolicy< 
//   RAJA::statement::CudaKernelAsync<
//     RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
// 			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
// 					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<16 >,
// 								   RAJA::statement::Lambda<0> >>>>>;

using RHS4_EXEC_POL = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
			      RAJA::statement::Tile<2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_z_loop,
          RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
            RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
				 RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
						      RAJA::statement::Lambda<0> >>>>>>>>;


using RHS4_EXEC_POL_ASYNC_OLDE = RAJA::KernelPolicy<
    RAJA::statement::CudaKernelAsync<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
			      RAJA::statement::Tile<2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_z_loop,
          RAJA::statement::For<0, RAJA::cuda_thread_y_loop,
            RAJA::statement::For<1, RAJA::cuda_thread_x_loop,
				 RAJA::statement::For<2, RAJA::cuda_thread_z_loop,
						      RAJA::statement::Lambda<0> >>>>>>>>;

using RHS4_EXEC_POL_ASYNC = RAJA::KernelPolicy<
	  RAJA::statement::CudaKernelFixedAsync<256,
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
			      RAJA::statement::Tile<2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
          RAJA::statement::For<0, RAJA::cuda_thread_z_direct,
            RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
				 RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
						      RAJA::statement::Lambda<0> >>>>>>>>;


using CONSINTP_EXEC_POL1 = RAJA::KernelPolicy<
  RAJA::statement::CudaKernel<
    RAJA::statement::Tile<0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
			  RAJA::statement::Tile<1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
						RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
								     RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
											  RAJA::statement::Lambda<0> >>>>>>;


using ODDIODDJ_EXEC_POL1_ASYNC = RAJA::KernelPolicy<
  RAJA::statement::CudaKernelFixedAsync<256,
    RAJA::statement::Tile<0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
			  RAJA::statement::Tile<1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
						RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
								     RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
											  RAJA::statement::Lambda<0> >>>>>>;

using ODDIODDJ_EXEC_POL2_ASYNC = RHS4_EXEC_POL_ASYNC;



using EVENIODDJ_EXEC_POL_ASYNC= 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernelAsync<
  RAJA::statement::For<1, RAJA::cuda_block_x_loop, 
  RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
  RAJA::statement::Lambda<0> >>>>;

using EVENIEVENJ_EXEC_POL_ASYNC =  EVENIODDJ_EXEC_POL_ASYNC;


using ODDIEVENJ_EXEC_POL1_ASYNC=  ICSTRESS_EXEC_POL_ASYNC;
 

using ODDIEVENJ_EXEC_POL2_ASYNC  =  RHS4_EXEC_POL_ASYNC;



using XRHS_POL = 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernel<
       RAJA::statement::For<0, RAJA::cuda_block_exec, 
			    RAJA::statement::For<1, RAJA::cuda_block_exec, 
						 RAJA::statement::For<2, RAJA::cuda_thread_exec,
								      RAJA::statement::Lambda<0> >>>>>;

using XRHS_POL_ASYNC = 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernelAsync<
       RAJA::statement::For<0, RAJA::cuda_block_x_loop, 
			    RAJA::statement::For<1, RAJA::cuda_block_y_loop, 
						 RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
								      RAJA::statement::Lambda<0> >>>>>;




using TWILIGHTSG_POL = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
			      RAJA::statement::Tile<2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_x_loop,
          RAJA::statement::For<0, RAJA::cuda_thread_z_direct,
            RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
				 RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
						      RAJA::statement::Lambda<0> >>>>>>>>;

using CONSINTP_EXEC_POL4  = ICSTRESS_EXEC_POL;



using CONSINTP_EXEC_POL5 =  ICSTRESS_EXEC_POL;
 


using PRELIM_CORR_EXEC_POL1 =  DEFAULT_LOOP2X;
using PRELIM_CORR_EXEC_POL1_ASYNC =  DEFAULT_LOOP2X_ASYNC;

 

using PRELIM_PRED_EXEC_POL1 =  ICSTRESS_EXEC_POL;
using PRELIM_PRED_EXEC_POL1_ASYNC =  ICSTRESS_EXEC_POL_ASYNC;
			 
//using ENFORCEBC_CORR_EXEC_POL1 =  ICSTRESS_EXEC_POL;

using ENFORCEBC_CORR_EXEC_POL1 = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelFixed<256,
    RAJA::statement::Tile<1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
			  RAJA::statement::Tile<0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
						RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
								     RAJA::statement::For<0, RAJA::cuda_thread_y_direct,		 
											  RAJA::statement::Lambda<0> >>>>>>;

using BCFORT_EXEC_POL1 =  RHS4_EXEC_POL;
using BCFORT_EXEC_POL2 = ICSTRESS_EXEC_POL;

using ENERGY4CI_EXEC_POL =  RHS4_EXEC_POL;

#endif


//*****************************************************************************************
//****************************** OMP POLICIES *******************************************






#else

#define SW4_FORCEINLINE


using REDUCTION_POLICY = RAJA::omp_reduce;

using DEFAULT_LOOP1 = RAJA::omp_parallel_for_exec;
using DEFAULT_LOOP1_ASYNC = RAJA::omp_parallel_for_exec;

using DEFAULT_LOOP2 = 
    RAJA::KernelPolicy<
      RAJA::statement::For<1, RAJA::omp_parallel_for_exec,
  RAJA::statement::For<0, RAJA::omp_parallel_for_exec,           
          RAJA::statement::Lambda<0>
        > 
      > 
    >;

using DEFAULT_LOOP2X = DEFAULT_LOOP2;
using DEFAULT_LOOP2X_ASYNC = DEFAULT_LOOP2;

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
using XRHS_POL_ASYNC = DEFAULT_LOOP3;

using RHS4_EXEC_POL = DEFAULT_LOOP3;
using RHS4_EXEC_POL_ASYNC = DEFAULT_LOOP3;


using ICSTRESS_EXEC_POL = DEFAULT_LOOP2;
using ICSTRESS_EXEC_POL_ASYNC = DEFAULT_LOOP2;
 


using CONSINTP_EXEC_POL1 = DEFAULT_LOOP2;


using CONSINTP_EXEC_POL3 = DEFAULT_LOOP3;
  


using CONSINTP_EXEC_POL4 = DEFAULT_LOOP2;


using PRELIM_CORR_EXEC_POL1 = DEFAULT_LOOP2;
using PRELIM_CORR_EXEC_POL1_ASYNC = DEFAULT_LOOP2;

using PRELIM_PRED_EXEC_POL1 = DEFAULT_LOOP2;
using PRELIM_PRED_EXEC_POL1_ASYNC = DEFAULT_LOOP2;

using ENFORCEBC_CORR_EXEC_POL1 = DEFAULT_LOOP2;


// Policy used in EW::get_exact_point_source
using GEPS_EXEC_POL = DEFAULT_LOOP3;
 

typedef RAJA::omp_parallel_for_exec SARRAY_LOOP_POL1;

using SARRAY_LOOP_POL2 = DEFAULT_LOOP4;
 

typedef RAJA::omp_parallel_for_exec PREDFORT_LOOP_POL;
typedef RAJA::omp_parallel_for_exec PREDFORT_LOOP_POL_ASYNC;

typedef RAJA::omp_parallel_for_exec CORRFORT_LOOP_POL;
typedef RAJA::omp_parallel_for_exec CORRFORT_LOOP_POL_ASYNC;

typedef RAJA::omp_parallel_for_exec DPDMTFORT_LOOP_POL;
typedef RAJA::omp_parallel_for_exec DPDMTFORT_LOOP_POL_ASYNC;


using DPDMT_WIND_LOOP_POL = DEFAULT_LOOP4;
using DPDMT_WIND_LOOP_POL_ASYNC  = DEFAULT_LOOP4;

using COPY_KPLANE_EXEC_POL = DEFAULT_LOOP3;
 
using ENERGY4CI_EXEC_POL = DEFAULT_LOOP3;


using ODDIODDJ_EXEC_POL1_ASYNC = DEFAULT_LOOP2;
 

using ODDIODDJ_EXEC_POL2_ASYNC = DEFAULT_LOOP3;
 

using ODDIEVENJ_EXEC_POL1_ASYNC = DEFAULT_LOOP2;
 
using ODDIEVENJ_EXEC_POL2 = DEFAULT_LOOP3;


using BCFORT_EXEC_POL1 = DEFAULT_LOOP3;


using BCFORT_EXEC_POL2 = DEFAULT_LOOP2;

using EVENIODDJ_EXEC_POL_ASYNC = DEFAULT_LOOP2;

using EVENIEVENJ_EXEC_POL_ASYNC = DEFAULT_LOOP2;


using TWILIGHTSG_POL = DEFAULT_LOOP3;
#define SYNC_DEVICE
#define SYNC_STREAM

#endif

#endif
