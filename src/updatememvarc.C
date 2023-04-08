//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// # Produced at the Lawrence Livermore National Laboratory.
// #
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// #
// # LLNL-CODE-643337
// #
// # All rights reserved.
// #
// # This file is part of SW4, Version: 1.0
// #
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
// #
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991.
// #
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details.
// #
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

#include <sys/types.h>

#include "Mspace.h"
#include "caliper.h"
#include "foralls.h"
#include "policies.h"
#include "sw4.h"
// 90 GB/s on V100 !! with 16 4 4
// 500 GBs on V100 with 16 16 4. But this is slower fro long runs !!
void memvar_pred_fort_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                         sw4_type kfirst, sw4_type klast, float_sw4* __restrict__ alp,
                         float_sw4* __restrict__ alm, float_sw4* __restrict__ u,
                         float_sw4 omega, float_sw4 dt, sw4_type domain) {
  SW4_MARK_FUNCTION;
  //***********************************************************************
  //***
  //*** domain = 0 --> Entire domain
  //*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  //*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  //***
  // AP: this routine implementes the 2nd order predictor step for evolving the
  // memory variables
  //***********************************************************************

  float_sw4 dto = dt * omega;
  float_sw4 icp = 1 / (1.0 / 2 + 1 / (2 * dto));
  float_sw4 cm = 0.5 - 1 / (2 * dto);
  sw4_type k1, k2;
  if (domain == 0) {
    k1 = kfirst;
    k2 = klast;
  } else if (domain == 1) {
    k1 = kfirst;
    k2 = kfirst + 2;
  } else if (domain == 2) {
    k1 = klast - 2;
    k2 = klast;
  }
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const sw4_type base = -ifirst - ni * jfirst - nij * kfirst;
  //    for( sw4_type c=0 ; c < 3 ;c++)
  // #pragma omp parallel for
  //       for( sw4_type k=k1 ; k <= k2 ; k++)
  // 	 for( sw4_type j=jfirst ; j<= jlast; j++ )
  // 	    for( sw4_type i=ifirst ; i<= ilast; i++ )
  // 	    {

  ASSERT_MANAGED(alp);
  ASSERT_MANAGED(alm);
  ASSERT_MANAGED(u);
  // PREFETCH(alp);
  // PREFETCH(alm);
  // PREFETCH(u);
  RAJA::RangeSegment i_range(ifirst, ilast + 1);
  RAJA::RangeSegment j_range(jfirst, jlast + 1);
  RAJA::RangeSegment k_range(k1, k2 + 1);
  //  RAJA::RangeSegment c_range(0, 3);

  // Note both POL3 and POL3X_ take about the same time on Hayward with 16 ranks

#if !defined(RAJA_ONLY) && defined(ENABLE_GPU)
  Range<16> I(ifirst, ilast + 1);
  Range<4> J(jfirst, jlast + 1);
  Range<4> K(k1, k2 + 1);
  forall3async(I, J, K, [=] RAJA_DEVICE(sw4_type i, sw4_type j, sw4_type k) {
#else
  RAJA::kernel<MPFC_POL_ASYNC>(
      RAJA::make_tuple(k_range, j_range, i_range),
      [=] RAJA_DEVICE(sw4_type k, sw4_type j, sw4_type i) {
#endif
    size_t ind = base + i + ni * j + nij * k;
#pragma unroll 3
    for (sw4_type c = 0; c < 3; c++)
      alp[ind + c * nijk] =
          icp * (-cm * alm[ind + c * nijk] + u[ind + c * nijk]);
  });  // SYNC_STREAM;

  return;
}

//-----------------------------------------------------------------------
// 618 Gb/s on V100
void memvar_corr_fort_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                         sw4_type kfirst, sw4_type klast, float_sw4* __restrict__ alp,
                         float_sw4* __restrict__ alm,
                         float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                         float_sw4* __restrict__ um, float_sw4 omega,
                         float_sw4 dt, sw4_type domain) {
  SW4_MARK_FUNCTION;
  //***********************************************************************
  //***
  //*** domain = 0 --> Entire domain
  //*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  //*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  //
  // AP Apr. 3, 2017: corrector step for updating memory variables
  // AP June 14, 2017: make corrector step independent of predictor step to
  // simplify the mesh refinement algorithm
  //***********************************************************************

  const float_sw4 i6 = 1.0 / 6;
  float_sw4 dto = dt * omega;

  float_sw4 icp = 1 / (0.5 + 1 / (2 * dto) + dto / 4 + dto * dto / 12);
  float_sw4 cm = 1 / (2 * dto) + dto / 4 - 0.5 - dto * dto / 12;

  sw4_type k1, k2;
  if (domain == 0) {
    k1 = kfirst;
    k2 = klast;
  } else if (domain == 1) {
    k1 = kfirst;
    k2 = kfirst + 2;
  } else if (domain == 2) {
    k1 = klast - 2;
    k2 = klast;
  }
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const sw4_type base = -ifirst - ni * jfirst - nij * kfirst;
  //    for( sw4_type c=0 ; c < 3 ;c++)
  // #pragma omp parallel for
  //       for( sw4_type k=k1 ; k <= k2 ; k++)
  // 	 for( sw4_type j=jfirst ; j<= jlast; j++ )
  // 	    for( sw4_type i=ifirst ; i<= ilast; i++ )
  // 	    {
  RAJA::RangeSegment i_range(ifirst, ilast + 1);
  RAJA::RangeSegment j_range(jfirst, jlast + 1);
  RAJA::RangeSegment k_range(k1, k2 + 1);
  //  RAJA::RangeSegment c_range(0, 3);

#if !defined(RAJA_ONLY) && defined(ENABLE_GPU)
#ifdef ENABLE_CUDA
  Range<16> I(ifirst, ilast + 1);
  Range<4> J(jfirst, jlast + 1);
  Range<4> K(k1, k2 + 1);
#else
  // 64,2,2 is 196ms, 256,2,2 is 182ms, 512,2,1 is 181ms, 1024,1,1 is 169ms
  Range<1024> I(ifirst, ilast + 1);
  Range<1> J(jfirst, jlast + 1);
  Range<1> K(k1, k2 + 1);
#endif
  forall3async(I, J, K, [=] RAJA_DEVICE(sw4_type i, sw4_type j, sw4_type k) {
#else
  RAJA::kernel<MPFC_POL_ASYNC>(
      RAJA::make_tuple(k_range, j_range, i_range),
      [=] RAJA_DEVICE(sw4_type k, sw4_type j, sw4_type i) {
#endif
    size_t ind = base + i + ni * j + nij * k;
    // Note that alp is ASSIGNED by this formula
#pragma unroll 3
    for (sw4_type c = 0; c < 3; c++)
      alp[ind + c * nijk] =
          icp * (cm * alm[ind + c * nijk] + u[ind + c * nijk] +
                 i6 * (dto * dto * u[ind + c * nijk] +
                       dto * (up[ind + c * nijk] - um[ind + c * nijk]) +
                       (up[ind + c * nijk] - 2 * u[ind + c * nijk] +
                        um[ind + c * nijk])));
  });  // SYNC_STREAM;

  return;
}

//-----------------------------------------------------------------------
// 565 GB/s on V100
void memvar_corr_fort_wind_ci(
    sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst, sw4_type klast,
    float_sw4* __restrict__ alp, sw4_type d1b, sw4_type d1e, sw4_type d2b, sw4_type d2e, sw4_type d3b,
    sw4_type d3e, float_sw4* __restrict__ alm, float_sw4* __restrict__ up,
    float_sw4* __restrict__ u, float_sw4* __restrict__ um, float_sw4 omega,
    float_sw4 dt, sw4_type domain) {
  SW4_MARK_FUNCTION;
  //***********************************************************************
  //***
  //*** domain = 0 --> Entire domain
  //*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  //*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  //***
  // AP Apr. 3, 2017: corrector step for updating memory variables
  // AP June 14, 2017: make corrector step independent of predictor step to
  // simplify the mesh refinement algorithm
  //***********************************************************************

  const float_sw4 i6 = 1.0 / 6;
  float_sw4 dto = dt * omega;
  float_sw4 icp = 1 / (0.5 + 1 / (2 * dto) + dto / 4 + dto * dto / 12);
  float_sw4 cm = 1 / (2 * dto) + dto / 4 - 0.5 - dto * dto / 12;

  sw4_type k1, k2;
  if (domain == 0) {
    k1 = kfirst;
    k2 = klast;
  } else if (domain == 1) {
    k1 = kfirst;
    k2 = kfirst + 2;
  } else if (domain == 2) {
    k1 = klast - 2;
    k2 = klast;
  }
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const sw4_type base = -ifirst - ni * jfirst - nij * kfirst;

  const size_t dni = d1e - d1b + 1;
  const size_t dnij = dni * (d2e - d2b + 1);
  const size_t dnijk = dnij * (d3e - d3b + 1);
  const sw4_type dbase = -d1b - dni * d2b - dnij * d3b;

  //  real*8 up(3,d1b:d1e, d2b:d2e, d3b:d3e)
  //  real*8  u(3,d1b:d1e, d2b:d2e, d3b:d3e);
  //  real*8 um(3,d1b:d1e, d2b:d2e, d3b:d3e);
  //  real*8 alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast) // different sizes
  //  here real*8 alm(3,d1b:d1e, d2b:d2e, d3b:d3e)

  //   for( sw4_type c=0 ; c < 3 ;c++)
  // #pragma omp parallel for
  //       for( sw4_type k=k1 ; k <= k2 ; k++)
  // 	 for( sw4_type j=jfirst ; j<= jlast; j++ )
  // 	    for( sw4_type i=ifirst ; i<= ilast; i++ )

  RAJA::RangeSegment i_range(ifirst, ilast + 1);
  RAJA::RangeSegment j_range(jfirst, jlast + 1);
  RAJA::RangeSegment k_range(k1, k2 + 1);
  //  RAJA::RangeSegment c_range(0, 3);

  RAJA::kernel<XRHS_POL_ASYNC>(
      RAJA::make_tuple(k_range, j_range, i_range),
      [=] RAJA_DEVICE(sw4_type k, sw4_type j, sw4_type i) {
        size_t ind = base + i + ni * j + nij * k;
        size_t dind = dbase + i + dni * j + dnij * k;
#pragma unroll
        for (sw4_type c = 0; c < 3; c++)
          alp[ind + c * nijk] =
              icp * (cm * alm[dind + c * dnijk] + u[dind + c * dnijk] +
                     i6 * (dto * dto * u[dind + c * dnijk] +
                           dto * (up[dind + c * dnijk] - um[dind + c * dnijk]) +
                           (up[dind + c * dnijk] - 2 * u[dind + c * dnijk] +
                            um[dind + c * dnijk])));
      });
  // SYNC_STREAM;
}
