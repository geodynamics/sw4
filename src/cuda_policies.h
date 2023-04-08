#ifndef __CUDA_POLICIES_H__
#define __CUDA_POLICIES_H__
#include "RAJA/RAJA.hpp"

#define SW4_FORCEINLINE __forceinline__
#define SYNC_DEVICE SW4_CheckDeviceError(cudaDeviceSynchronize())
#define SYNC_STREAM SW4_CheckDeviceError(cudaStreamSynchronize(0))
#define SW4_PEEK SW4_CheckDeviceError(cudaPeekAtLastError());

#define SW4_MALLOC_MANAGED(addr, size) (cudaMallocManaged(addr, size))
#define SW4_MALLOC_DEVICE(addr, size) (cudaMalloc(addr, size))
#define SW4_MALLOC_PINNED(addr, size) \
  (cudaHostAlloc(addr, size, cudaHostAllocMapped))

#define SW4_FREE_MANAGED(addr) (cudaFree(addr))
#define SW4_FREE_DEVICE(addr) (cudaFree(addr))
#define SW4_FREE_PINNED(addr) (cudaFreeHost(addr))

#define SW4_DEVICE_SUCCESS cudaSuccess
//   SW4_CheckDeviceError(cudaStreamSynchronize(0));
typedef RAJA::cuda_exec<1024> DEFAULT_LOOP1;
typedef RAJA::cuda_exec<1024, true> DEFAULT_LOOP1_ASYNC;
using REDUCTION_POLICY = RAJA::cuda_reduce;

typedef RAJA::cuda_exec<1024> PREDFORT_LOOP_POL;
typedef RAJA::cuda_exec<512, true> PREDFORT_LOOP_POL_ASYNC;

typedef RAJA::cuda_exec<1024> CORRFORT_LOOP_POL;
typedef RAJA::cuda_exec<512, true> CORRFORT_LOOP_POL_ASYNC;

typedef RAJA::cuda_exec<1024> DPDMTFORT_LOOP_POL;

typedef RAJA::cuda_exec<256, true> DPDMTFORT_LOOP_POL_ASYNC;

typedef RAJA::cuda_exec<1024> SARRAY_LOOP_POL1;

using DEFAULT_LOOP2X = RAJA::KernelPolicy<RAJA::statement::CudaKernel<
    RAJA::statement::For<1, RAJA::cuda_block_x_loop,
                         RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                                              RAJA::statement::Lambda<0>>>>>;

using DEFAULT_LOOP2X_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::For<
        1, RAJA::cuda_block_x_loop,
        RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                             RAJA::statement::Lambda<0>>>>>;

using DEFAULT_LOOP3 = RAJA::KernelPolicy<RAJA::statement::CudaKernelFixed<
    256,
    RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_z_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using DEFAULT_LOOP3_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
        256,
        RAJA::statement::Tile<
            0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
            RAJA::statement::Tile<
                1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
                RAJA::statement::Tile<
                    2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                    RAJA::statement::For<
                        0, RAJA::cuda_thread_z_direct,
                        RAJA::statement::For<
                            1, RAJA::cuda_thread_y_direct,
                            RAJA::statement::For<
                                2, RAJA::cuda_thread_x_direct,
                                RAJA::statement::Lambda<0>>>>>>>>>;

using SARRAY_LOOP_POL2 =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<1>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<1>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<1024>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_x_direct,
                        RAJA::statement::For<
                            2, RAJA::cuda_thread_z_direct,
                            RAJA::statement::For<
                                3, RAJA::seq_exec,
                                RAJA::statement::Lambda<0>>>>>>>>>>;

// using RHS4_EXEC_POL =
//   RAJA::KernelPolicy<
//   RAJA::statement::CudaKernel<
//     RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>,
// 			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>,
// 					      RAJA::statement::For<2,
// RAJA::cuda_threadblock_exec<16 >,
// RAJA::statement::Lambda<0> >>>>>;

using ICSTRESS_EXEC_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::For<
                0, RAJA::cuda_thread_x_loop,
                RAJA::statement::For<1, RAJA::cuda_thread_y_loop,
                                     RAJA::statement::Lambda<0>>>>>>>;

using ICSTRESS_EXEC_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
        256,
        RAJA::statement::Tile<
            0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_x_direct,
                    RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                                         RAJA::statement::Lambda<0>>>>>>>;

// using RHS4_EXEC_POL_ASYNC =
//   RAJA::KernelPolicy<
//   RAJA::statement::CudaKernelAsync<
//     RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>,
// 			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>,
// 					      RAJA::statement::For<2,
// RAJA::cuda_threadblock_exec<16 >,
// RAJA::statement::Lambda<0> >>>>>;

using RHS4_EXEC_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_x_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using RHS4_EXEC_POL_ASYNC_OLDE =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_loop,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_x_loop,
                        RAJA::statement::For<2, RAJA::cuda_thread_z_loop,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using RHS4_EXEC_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
        256,
        RAJA::statement::Tile<
            0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
            RAJA::statement::Tile<
                1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
                RAJA::statement::Tile<
                    2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                    RAJA::statement::For<
                        0, RAJA::cuda_thread_z_direct,
                        RAJA::statement::For<
                            1, RAJA::cuda_thread_y_direct,
                            RAJA::statement::For<
                                2, RAJA::cuda_thread_x_direct,
                                RAJA::statement::Lambda<0>>>>>>>>>;

using CONSSW4_TYPEP_EXEC_POL1 =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
            RAJA::statement::For<
                0, RAJA::cuda_thread_y_direct,
                RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
                                     RAJA::statement::Lambda<0>>>>>>>;

using ODDIODDJ_EXEC_POL1_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
        256,
        RAJA::statement::Tile<
            0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_direct,
                    RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
                                         RAJA::statement::Lambda<0>>>>>>>;

using ODDIODDJ_EXEC_POL2_ASYNC = RHS4_EXEC_POL_ASYNC;

// Check this
using EVENIODDJ_EXEC_POL_ASYNC_ORG =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::For<
        1, RAJA::cuda_block_x_loop,
        RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                             RAJA::statement::Lambda<0>>>>>;

using EVENIODDJ_EXEC_POL_ASYNC = ODDIODDJ_EXEC_POL1_ASYNC;

using EVENIEVENJ_EXEC_POL_ASYNC = EVENIODDJ_EXEC_POL_ASYNC;

using ODDIEVENJ_EXEC_POL1_ASYNC = ICSTRESS_EXEC_POL_ASYNC;

using ODDIEVENJ_EXEC_POL2_ASYNC = RHS4_EXEC_POL_ASYNC;

using XRHS_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::For<
        0, RAJA::cuda_block_exec,
        RAJA::statement::For<
            1, RAJA::cuda_block_exec,
            RAJA::statement::For<2, RAJA::cuda_thread_exec,
                                 RAJA::statement::Lambda<0>>>>>>;

using XRHS_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::For<
        0, RAJA::cuda_block_x_loop,
        RAJA::statement::For<
            1, RAJA::cuda_block_y_loop,
            RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
                                 RAJA::statement::Lambda<0>>>>>>;

using TWILIGHTSG_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_x_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_z_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using CONSSW4_TYPEP_EXEC_POL4 = ICSTRESS_EXEC_POL;

using CONSSW4_TYPEP_EXEC_POL5 = ICSTRESS_EXEC_POL;

using PRELIM_CORR_EXEC_POL1 = DEFAULT_LOOP2X;
using PRELIM_CORR_EXEC_POL1_ASYNC = DEFAULT_LOOP2X_ASYNC;

using PRELIM_PRED_EXEC_POL1 = ICSTRESS_EXEC_POL;
using PRELIM_PRED_EXEC_POL1_ASYNC = ICSTRESS_EXEC_POL_ASYNC;

// using ENFORCEBC_CORR_EXEC_POL1 =  ICSTRESS_EXEC_POL;

using ENFORCEBC_CORR_EXEC_POL1 =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixed<
        256,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
                RAJA::statement::For<
                    1, RAJA::cuda_thread_x_direct,
                    RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
                                         RAJA::statement::Lambda<0>>>>>>>;

using BCFORT_EXEC_POL1 = RHS4_EXEC_POL;
using BCFORT_EXEC_POL2 = ICSTRESS_EXEC_POL;

using ENERGY4CI_EXEC_POL = RHS4_EXEC_POL;

// Next 4 in solve.C
using DHI_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::For<
        0, RAJA::cuda_block_x_loop,
        RAJA::statement::For<
            1, RAJA::cuda_block_y_loop,
            RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
                                 RAJA::statement::Lambda<0>>>>>>;

using GIG_POL_ASYNC = RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<
    RAJA::statement::For<0, RAJA::cuda_block_x_loop,
                         RAJA::statement::For<1, RAJA::cuda_thread_x_loop,
                                              RAJA::statement::Lambda<0>>>>>;

using AVS_POL_ASYNC = RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
    256, RAJA::statement::Tile<
             0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
             RAJA::statement::Tile<
                 1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                 RAJA::statement::For<
                     0, RAJA::cuda_thread_y_direct,
                     RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
                                          RAJA::statement::Lambda<0>>>>>>>;

using EBFA_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::For<
                0, RAJA::cuda_thread_x_direct,
                RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                                     RAJA::statement::Lambda<0>>>>>>>;

// void EW::get_exact_point_source in EW.C
using GEPS_POL = RAJA::KernelPolicy<
    // RAJA::statement::CudaKernelExt<RAJA::cuda_explicit_launch<false, 0,
    // 256>,
    RAJA::statement::CudaKernelFixed<
        256,
        // RAJA::statement::CudaKernel<
        // RAJA::statement::CudaKernelOcc<
        RAJA::statement::Tile<
            0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
                RAJA::statement::Tile<
                    2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_z_loop,
                    RAJA::statement::For<
                        0, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<
                            1, RAJA::cuda_thread_x_direct,
                            RAJA::statement::For<
                                2, RAJA::cuda_thread_z_direct,
                                RAJA::statement::Lambda<0>>>>>>>>>;

// CurvilinearSw4_Typeerface2::bnd_zero
using BZ_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_z_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

// CurvilinearSw4_Typeerface2::injection
using INJ_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_z_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using INJ_POL2_ASYNC = RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<
    RAJA::statement::For<1, RAJA::cuda_block_x_loop,
                         RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                                              RAJA::statement::Lambda<0>>>>>;

// CurvilinearSw4_Typeerface2::communicate_array

using CA_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_z_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

// Sarray::assign

using SAA_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            3, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    1, RAJA::cuda_thread_y_direct,
                    RAJA::statement::For<
                        3, RAJA::cuda_thread_x_direct,
                        RAJA::statement::For<
                            2, RAJA::cuda_thread_z_direct,
                            RAJA::statement::For<
                                0, RAJA::seq_exec,
                                RAJA::statement::Lambda<0>>>>>>>>>>;

// Sarray::insert_sw4_typeersection(
using SII_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_z_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

// TestEcons::get_ubnd(
using TGU_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_z_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_x_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

// in addmemvarforcing2.C
using AMVPCa_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_loop,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_x_loop,
                        RAJA::statement::For<2, RAJA::cuda_thread_z_loop,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using AMVPCu_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_loop,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_x_loop,
                        RAJA::statement::For<2, RAJA::cuda_thread_z_loop,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using AMVC2Ca_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_loop,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_x_loop,
                        RAJA::statement::For<2, RAJA::cuda_thread_z_loop,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using AMVC2Cu_POL =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_loop,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_x_loop,
                        RAJA::statement::For<2, RAJA::cuda_thread_z_loop,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

// In addsg4windc.C
using ASG4WC_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::For<
        2, RAJA::cuda_block_z_loop,
        RAJA::statement::For<
            1, RAJA::cuda_block_y_loop,
            RAJA::statement::For<
                0, RAJA::cuda_thread_x_direct,
                RAJA::statement::For<3, RAJA::seq_exec,
                                     RAJA::statement::Lambda<0>>>>>>>;

// In addsgdc.C
using ADDSGD_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
        256,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                3, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_x_loop,
                RAJA::statement::Tile<
                    2, RAJA::statement::tile_fixed<1>, RAJA::cuda_block_z_loop,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_loop,
                        RAJA::statement::For<
                            3, RAJA::cuda_thread_x_loop,
                            RAJA::statement::For<
                                2, RAJA::cuda_thread_z_loop,
                                RAJA::statement::For<
                                    0, RAJA::seq_exec,
                                    RAJA::statement::Lambda<0>>>>>>>>>>;

using ADDSGD_POL2_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
        256,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<
                3, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                RAJA::statement::Tile<
                    2, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_y_direct,
                        RAJA::statement::For<
                            3, RAJA::cuda_thread_x_direct,
                            RAJA::statement::For<
                                2, RAJA::cuda_thread_z_direct,
                                RAJA::statement::For<
                                    0, RAJA::seq_exec,
                                    RAJA::statement::Lambda<0>>>>>>>>>>;

// in bcforce.C
using BCFORT_EXEC_POL2_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_y_direct,
                    RAJA::statement::For<
                        1, RAJA::cuda_thread_x_direct,
                        RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
                                             RAJA::statement::Lambda<0>>>>>>>>>;

using BCFORT_EXEC_POL3_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::For<
                0, RAJA::cuda_thread_x_direct,
                RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                                     RAJA::statement::Lambda<0>>>>>>>;

// in curvilinear4sgc.C
using CURV_POL_ORG =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::For<
        0, RAJA::cuda_block_exec,
        RAJA::statement::For<
            1, RAJA::cuda_block_exec,
            RAJA::statement::For<2, RAJA::cuda_thread_exec,
                                 RAJA::statement::Lambda<0>>>>>>;
using CURV_POL = DEFAULT_LOOP3;
// in parallelStuff.C
using BUFFER_POL = RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<
    RAJA::statement::For<1, RAJA::cuda_block_x_loop,
                         RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                                              RAJA::statement::Lambda<0>>>>>;

// in rhs3cuvilinearsgc.C
using RHS4CU_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::For<
                0, RAJA::cuda_thread_x_direct,
                RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                                     RAJA::statement::Lambda<0>>>>>>>;

// in rhs4th3fortc.C
using XRHS_POL2 =
    RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::For<
        0, RAJA::cuda_block_x_loop,
        RAJA::statement::For<
            1, RAJA::cuda_block_y_loop,
            RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
                                 RAJA::statement::Lambda<0>>>>>>;
using RHS4TH3_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
        256,
        RAJA::statement::Tile<
            0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
            RAJA::statement::Tile<
                1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
                RAJA::statement::Tile<
                    2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                    RAJA::statement::For<
                        0, RAJA::cuda_thread_z_direct,
                        RAJA::statement::For<
                            1, RAJA::cuda_thread_y_direct,
                            RAJA::statement::For<
                                2, RAJA::cuda_thread_x_direct,
                                RAJA::statement::Lambda<0>>>>>>>>>;

using RHS4TH3_POL2_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
        256,
        RAJA::statement::Tile<
            0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
            RAJA::statement::Tile<
                1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
                RAJA::statement::Tile<
                    2, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
                    RAJA::statement::For<
                        0, RAJA::cuda_thread_z_direct,
                        RAJA::statement::For<
                            1, RAJA::cuda_thread_y_direct,
                            RAJA::statement::For<
                                2, RAJA::cuda_thread_x_direct,
                                RAJA::statement::Lambda<0>>>>>>>>>;

/* using RHS4TH3_POL2_ASYNC = RAJA::statement::CudaKernelFixedAsync< */
/*   256, */
/*   RAJA::statement::Tile< */
/*   0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop, */
/*   RAJA::statement::Tile< */
/*   1, RAJA::statement::tile_fixed<4>, */
/*   RAJA::cuda_block_y_loop, */
/*   RAJA::statement::Tile< */
/*   2, RAJA::statement::tile_fixed<16>, */
/*   RAJA::cuda_block_x_loop, */
/*   RAJA::statement::For< */
/*   0, RAJA::cuda_thread_z_direct, */
/*   RAJA::statement::For< */
/*   1, RAJA::cuda_thread_y_direct, */
/*   RAJA::statement::For< */
/*   2, RAJA::cuda_thread_x_direct, */
/*   RAJA::statement::Lambda<0>>>>>>>>>; */

using VBSC_POL = DEFAULT_LOOP2X;

using AFCC_POL_ASYNC =
    RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
        0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
        RAJA::statement::Tile<
            1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::For<
                0, RAJA::cuda_thread_x_direct,
                RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                                     RAJA::statement::Lambda<0>>>>>>>;

// In updatememvarc.C
/* using MPFC_POL_ASYNC = */
/*     RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::For<
 */
/*         0, RAJA::cuda_threadblock_exec<1>, */
/*         RAJA::statement::For< */
/*             1, RAJA::cuda_threadblock_exec<1>, */
/*             RAJA::statement::For<2, RAJA::cuda_threadblock_exec<1024>, */
/*                                  RAJA::statement::Lambda<0>>>>>>; */

using MPFC_POL_ASYNC = DEFAULT_LOOP3_ASYNC;

// IN EW.C
using FORCE_LOOP_ASYNC = RAJA::cuda_exec<32, true>;
using FORCETT_LOOP_ASYNC = RAJA::cuda_exec<32, true>;

using COPY_KPLANE_EXEC_POL = DEFAULT_LOOP3;

using DPDMT_WIND_LOOP_POL_ASYNC = DEFAULT_LOOP3;
#endif
