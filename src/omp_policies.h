#ifndef __OMP_POLICIES_H__
#define __OMP_POLICIES_H__
#include "RAJA/RAJA.hpp"

#define SW4_FORCEINLINE
#define SW4_PEEK

using REDUCTION_POLICY = RAJA::omp_reduce;

using DEFAULT_LOOP1 = RAJA::omp_parallel_for_exec;
using DEFAULT_LOOP1_ASYNC = RAJA::omp_parallel_for_exec;

using DEFAULT_LOOP2 = RAJA::KernelPolicy<RAJA::statement::For<
    1, RAJA::omp_parallel_for_exec,
    RAJA::statement::For<0, RAJA::seq_exec, RAJA::statement::Lambda<0>>>>;

using DEFAULT_LOOP2X = DEFAULT_LOOP2;
using DEFAULT_LOOP2X_ASYNC = DEFAULT_LOOP2;

using DEFAULT_LOOP3 = RAJA::KernelPolicy<RAJA::statement::For<
    2, RAJA::omp_parallel_for_exec,
    RAJA::statement::For<
        1, RAJA::seq_exec,
        RAJA::statement::For<0, RAJA::simd_exec, RAJA::statement::Lambda<0>>>>>;

using DEFAULT_LOOP4 = RAJA::KernelPolicy<RAJA::statement::For<
    3, RAJA::omp_parallel_for_exec,
    RAJA::statement::For<
        2, RAJA::seq_exec,
        RAJA::statement::For<
            1, RAJA::seq_exec,
            RAJA::statement::For<0, RAJA::simd_exec,
                                 RAJA::statement::Lambda<0>>>>>>;

using ADDSGD_POL_ASYNC = RAJA::KernelPolicy<RAJA::statement::For<
    0, RAJA::seq_exec,
    RAJA::statement::For<
        1, RAJA::omp_parallel_for_exec,
        RAJA::statement::For<
            2, RAJA::seq_exec,
            RAJA::statement::For<3, RAJA::simd_exec,
                                 RAJA::statement::Lambda<0>>>>>>;

using XRHS_POL = RAJA::KernelPolicy<RAJA::statement::For<
    0, RAJA::omp_parallel_for_exec,
    RAJA::statement::For<
        1, RAJA::seq_exec,
        RAJA::statement::For<2, RAJA::simd_exec, RAJA::statement::Lambda<0>>>>>;
using XRHS_POL_ASYNC = XRHS_POL;

using RHS4_EXEC_POL = XRHS_POL;
using RHS4_EXEC_POL_ASYNC = XRHS_POL;

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
using DPDMT_WIND_LOOP_POL_ASYNC = DEFAULT_LOOP4;

using COPY_KPLANE_EXEC_POL = DEFAULT_LOOP3;

// using ENERGY4CI_EXEC_POL = DEFAULT_LOOP3;

using ENERGY4CI_EXEC_POL = RAJA::KernelPolicy<RAJA::statement::For<
    0, RAJA::omp_parallel_for_exec,
    RAJA::statement::For<
        1, RAJA::seq_exec,
        RAJA::statement::For<2, RAJA::seq_exec,  // Changing this to simd_exec
                                                 // caused the 4 energy tests in
                                                 // pytest to fail Ramesh
                                                 // Pankajakshan April 26 2019
                             RAJA::statement::Lambda<0>>>>>;

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

// Next 4 in solve.C
using DHI_POL_ASYNC = DEFAULT_LOOP3;

using GIG_POL = DEFAULT_LOOP2;

using AVS_POL_ASYNC = DEFAULT_LOOP2;

using EBFA_POL = DEFAULT_LOOP2;

// void EW::get_exact_point_source in EW.C
using GEPS_POL = DEFAULT_LOOP3;

// CurvilinearInterface2.C

using BZ_POL_ASYNC = DEFAULT_LOOP3;

// Injection in CurvilinearInterface2.C

using INJ_POL_ASYNC = DEFAULT_LOOP3;

using INJ_POL2_ASYNC = DEFAULT_LOOP2;
// CurvilinearInterface2::communicate_array

using CA_POL = DEFAULT_LOOP3;

// Sarray::assign

using SAA_POL = DEFAULT_LOOP4;

// Sarray::insert_intersection(
using SII_POL = DEFAULT_LOOP3;

// TestEcons::get_ubnd(
using TGU_POL_ASYNC = DEFAULT_LOOP3;

// in addmemvarforcing2.C
using AMVPCa_POL = DEFAULT_LOOP3;
using AMVPCu_POL = DEFAULT_LOOP3;
using AMVC2Ca_POL_ASYNC = DEFAULT_LOOP3;
using AMVC2Cu_POL = DEFAULT_LOOP3;

// In addsg4windc.C
using ASG4WC_POL_ASYNC = DEFAULT_LOOP4;

// In addsgdc.C
// using ADDSGD_POL_ASYNC = DEFAULT_LOOP4;

using ADDSGD_POL2_ASYNC = DEFAULT_LOOP4;

// in bcforce.C
using BCFORT_EXEC_POL2_ASYNC = DEFAULT_LOOP3;
using BCFORT_EXEC_POL3_ASYNC = DEFAULT_LOOP2;
// in curvilinear4sgc.C
using CURV_POL = XRHS_POL;

// in parallelStuff.C
using BUFFER_POL = DEFAULT_LOOP2;

// in rhs3cuvilinearsgc.C
using RHS4CU_POL_ASYNC = DEFAULT_LOOP2;

// in rhs4th3fortc.C
using XRHS_POL2 = XRHS_POL;
using RHS4TH3_POL_ASYNC = XRHS_POL;
using RHS4TH3_POL2_ASYNC = XRHS_POL;

using VBSC_POL = DEFAULT_LOOP2;

using AFCC_POL = DEFAULT_LOOP2;

using MPFC_POL_ASYNC = DEFAULT_LOOP3;

// in EW.C
using FORCE_LOOP_ASYNC = RAJA::omp_parallel_for_exec;
using FORCETT_LOOP_ASYNC = RAJA::omp_parallel_for_exec;

using GIG_POL_ASYNC = DEFAULT_LOOP2;
using AFCC_POL_ASYNC = DEFAULT_LOOP2;

#endif
