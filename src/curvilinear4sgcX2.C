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

#include "Mspace.h"
#include "caliper.h"
#include "foralls.h"
#include "policies.h"
#include "sw4.h"

void curvilinear4sgX1_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_mu,
    float_sw4* __restrict__ a_lambda, float_sw4* __restrict__ a_met,
    float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_lu, int* onesided,
    float_sw4* __restrict__ a_acof, float_sw4* __restrict__ a_bope,
    float_sw4* __restrict__ a_ghcof, float_sw4* __restrict__ a_acof_no_gp,
    float_sw4* __restrict__ a_ghcof_no_gp, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry, int nk, char op) {
  SW4_MARK_FUNCTION;
  //      subroutine CURVILINEAR4SG( ifirst, ilast, jfirst, jlast, kfirst,
  //     *                         klast, u, mu, la, met, jac, lu,
  //     *                         onesided, acof, bope, ghcof, strx, stry,
  //     *                         op )

  // Routine with supergrid stretchings strx and stry. No stretching
  // in z, since top is always topography, and bottom always interface
  // to a deeper Cartesian grid.
  // opcount:
  //      Interior (k>6), 2126 arithmetic ops.
  //      Boundary discretization (1<=k<=6 ), 6049 arithmetic ops.

  //   const float_sw4 a1 =0;
  float_sw4 a1 = 0;
  float_sw4 sgn = 1;
  if (op == '=') {
    a1 = 0;
    sgn = 1;
  } else if (op == '+') {
    a1 = 1;
    sgn = 1;
  } else if (op == '-') {
    a1 = 1;
    sgn = -1;
  }
#ifdef ENABLE_CUDA
#define NO_COLLAPSE 1
#endif

  const float_sw4 i6 = 1.0 / 6;
  const float_sw4 tf = 0.75;
  const float_sw4 c1 = 2.0 / 3;
  const float_sw4 c2 = -1.0 / 12;

  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;
  const int base4 = base - nijk;
  const int ifirst0 = ifirst;
  const int jfirst0 = jfirst;

  // Direct reuse of fortran code by these macro definitions:
  // Direct reuse of fortran code by these macro definitions:
#define mu(i, j, k) a_mu[base + (i) + ni * (j) + nij * (k)]
#define la(i, j, k) a_lambda[base + (i) + ni * (j) + nij * (k)]
#define jac(i, j, k) a_jac[base + (i) + ni * (j) + nij * (k)]
#define u(c, i, j, k) a_u[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define lu(c, i, j, k) a_lu[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define met(c, i, j, k) a_met[base4 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define strx(i) a_strx[i - ifirst0]
#define stry(j) a_stry[j - jfirst0]
#define acof(i, j, k) a_acof[(i - 1) + 6 * (j - 1) + 48 * (k - 1)]
#define bope(i, j) a_bope[i - 1 + 6 * (j - 1)]
#define ghcof(i) a_ghcof[i - 1]
#define acof_no_gp(i, j, k) a_acof_no_gp[(i - 1) + 6 * (j - 1) + 48 * (k - 1)]
#define ghcof_no_gp(i) a_ghcof_no_gp[i - 1]

#ifdef PEEKS_GALORE
  SW4_PEEK;
  SYNC_DEVICE;
#endif
  SW4_MARK_BEGIN("CURVI::cuvilinear4sgc");
  // SYNC_STREAM; // CURVI_CPU
  /// CURVIMR ADDITION
  if (onesided[5] == 1) {
// #pragma omp for
//     for (int k = nk - 5; k <= nk; k++)
//       for (int j = jfirst + 2; j <= jlast - 2; j++)
// #pragma omp simd
// #pragma ivdep
//         for (int i = ifirst + 2; i <= ilast - 2; i++) {
#if defined(NO_COLLAPSE)
    // LOOP -1
    // 32,4,2 is 4% slower. 32 4 4 does not fit
    Range<16> II(ifirst + 2, ilast - 1);
    Range<4> JJ(jfirst + 2, jlast - 1);
    Range<6> KK(nk - 5, nk + 1);
    // Register count goes upto 254. Runtime goes up by factor of 2.8X
    //     Range<16> JJ2(jfirst + 2, jlast - 1);
    //     forall2async(II, JJ2,[=] RAJA_DEVICE(int i, int j) {
    // #pragma unroll
    // 	for (int kk=-5;kk<1;kk++){
    // 	  int k=nk+kk;
    forall3async(II, JJ, KK, [=] RAJA_DEVICE(int i, int j, int k) {
    // forall3X results in a 2.5X slowdown even though registers drop from
    // 168 to 130
    // forall3X<256>(ifirst + 2, ilast - 1,jfirst + 2, jlast - 1,nk-5,nk+1,
    //    [=] RAJA_DEVICE(int i, int j, int k) {
#else
    RAJA::RangeSegment kk_range(nk - 5, nk + 1);
    RAJA::RangeSegment jj_range(jfirst + 2, jlast - 1);
    RAJA::RangeSegment ii_range(ifirst + 2, ilast - 1);
    RAJA::kernel<
        CURV_POL>(RAJA::make_tuple(kk_range, jj_range, ii_range), [=] RAJA_DEVICE(
                                                                      int k,
                                                                      int j,
                                                                      int i) {
#endif
      // 5 ops
      const int N = 4;
      float_sw4 ijac = strx(i) * stry(j) / jac(i, j, k);
      float_sw4 istry = 1 / (stry(j));
      float_sw4 istrx = 1 / (strx(i));
      float_sw4 istrxy = istry * istrx;

      float_sw4 r1 = 0, r2 = 0, r3 = 0;
      float_sw4 cof1, cof2, cof3, cof4, cof5;
      float_sw4 mux1, mux2, mux3, mux4;

      // pp derivative (u) (u-eq)
      // 53 ops, tot=58
      if (N == 1) {
        cof1 = (2 * mu(i - 2, j, k) + la(i - 2, j, k)) * met(1, i - 2, j, k) *
               met(1, i - 2, j, k) * strx(i - 2);
        cof2 = (2 * mu(i - 1, j, k) + la(i - 1, j, k)) * met(1, i - 1, j, k) *
               met(1, i - 1, j, k) * strx(i - 1);
        cof3 = (2 * mu(i, j, k) + la(i, j, k)) * met(1, i, j, k) *
               met(1, i, j, k) * strx(i);
        cof4 = (2 * mu(i + 1, j, k) + la(i + 1, j, k)) * met(1, i + 1, j, k) *
               met(1, i + 1, j, k) * strx(i + 1);
        cof5 = (2 * mu(i + 2, j, k) + la(i + 2, j, k)) * met(1, i + 2, j, k) *
               met(1, i + 2, j, k) * strx(i + 2);

        mux1 = cof2 - tf * (cof3 + cof1);
        mux2 = cof1 + cof4 + 3 * (cof3 + cof2);
        mux3 = cof2 + cof5 + 3 * (cof4 + cof3);
        mux4 = cof4 - tf * (cof3 + cof5);

        r1 = r1 + i6 *
                      (mux1 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                       mux2 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                       mux3 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                       mux4 * (u(1, i + 2, j, k) - u(1, i, j, k))) *
                      istry;

        // qq derivative (u) (u-eq)
        // 43 ops, tot=101
        cof1 = (mu(i, j - 2, k)) * met(1, i, j - 2, k) * met(1, i, j - 2, k) *
               stry(j - 2);
        cof2 = (mu(i, j - 1, k)) * met(1, i, j - 1, k) * met(1, i, j - 1, k) *
               stry(j - 1);
        cof3 = (mu(i, j, k)) * met(1, i, j, k) * met(1, i, j, k) * stry(j);
        cof4 = (mu(i, j + 1, k)) * met(1, i, j + 1, k) * met(1, i, j + 1, k) *
               stry(j + 1);
        cof5 = (mu(i, j + 2, k)) * met(1, i, j + 2, k) * met(1, i, j + 2, k) *
               stry(j + 2);

        mux1 = cof2 - tf * (cof3 + cof1);
        mux2 = cof1 + cof4 + 3 * (cof3 + cof2);
        mux3 = cof2 + cof5 + 3 * (cof4 + cof3);
        mux4 = cof4 - tf * (cof3 + cof5);

        r1 = r1 + i6 *
                      (mux1 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                       mux2 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                       mux3 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                       mux4 * (u(1, i, j + 2, k) - u(1, i, j, k))) *
                      istrx;
      }
      // pp derivative (v) (v-eq)
      // 43 ops, tot=144
      if (N == 2) {
        cof1 = (mu(i - 2, j, k)) * met(1, i - 2, j, k) * met(1, i - 2, j, k) *
               strx(i - 2);
        cof2 = (mu(i - 1, j, k)) * met(1, i - 1, j, k) * met(1, i - 1, j, k) *
               strx(i - 1);
        cof3 = (mu(i, j, k)) * met(1, i, j, k) * met(1, i, j, k) * strx(i);
        cof4 = (mu(i + 1, j, k)) * met(1, i + 1, j, k) * met(1, i + 1, j, k) *
               strx(i + 1);
        cof5 = (mu(i + 2, j, k)) * met(1, i + 2, j, k) * met(1, i + 2, j, k) *
               strx(i + 2);

        mux1 = cof2 - tf * (cof3 + cof1);
        mux2 = cof1 + cof4 + 3 * (cof3 + cof2);
        mux3 = cof2 + cof5 + 3 * (cof4 + cof3);
        mux4 = cof4 - tf * (cof3 + cof5);

        r2 = r2 + i6 *
                      (mux1 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                       mux2 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                       mux3 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                       mux4 * (u(2, i + 2, j, k) - u(2, i, j, k))) *
                      istry;

        // qq derivative (v) (v-eq)
        // 53 ops, tot=197
        cof1 = (2 * mu(i, j - 2, k) + la(i, j - 2, k)) * met(1, i, j - 2, k) *
               met(1, i, j - 2, k) * stry(j - 2);
        cof2 = (2 * mu(i, j - 1, k) + la(i, j - 1, k)) * met(1, i, j - 1, k) *
               met(1, i, j - 1, k) * stry(j - 1);
        cof3 = (2 * mu(i, j, k) + la(i, j, k)) * met(1, i, j, k) *
               met(1, i, j, k) * stry(j);
        cof4 = (2 * mu(i, j + 1, k) + la(i, j + 1, k)) * met(1, i, j + 1, k) *
               met(1, i, j + 1, k) * stry(j + 1);
        cof5 = (2 * mu(i, j + 2, k) + la(i, j + 2, k)) * met(1, i, j + 2, k) *
               met(1, i, j + 2, k) * stry(j + 2);
        mux1 = cof2 - tf * (cof3 + cof1);
        mux2 = cof1 + cof4 + 3 * (cof3 + cof2);
        mux3 = cof2 + cof5 + 3 * (cof4 + cof3);
        mux4 = cof4 - tf * (cof3 + cof5);

        r2 = r2 + i6 *
                      (mux1 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                       mux2 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                       mux3 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                       mux4 * (u(2, i, j + 2, k) - u(2, i, j, k))) *
                      istrx;
      }
      if (N == 3) {
        // pp derivative (w) (w-eq)
        // 43 ops, tot=240
        cof1 = (mu(i - 2, j, k)) * met(1, i - 2, j, k) * met(1, i - 2, j, k) *
               strx(i - 2);
        cof2 = (mu(i - 1, j, k)) * met(1, i - 1, j, k) * met(1, i - 1, j, k) *
               strx(i - 1);
        cof3 = (mu(i, j, k)) * met(1, i, j, k) * met(1, i, j, k) * strx(i);
        cof4 = (mu(i + 1, j, k)) * met(1, i + 1, j, k) * met(1, i + 1, j, k) *
               strx(i + 1);
        cof5 = (mu(i + 2, j, k)) * met(1, i + 2, j, k) * met(1, i + 2, j, k) *
               strx(i + 2);

        mux1 = cof2 - tf * (cof3 + cof1);
        mux2 = cof1 + cof4 + 3 * (cof3 + cof2);
        mux3 = cof2 + cof5 + 3 * (cof4 + cof3);
        mux4 = cof4 - tf * (cof3 + cof5);

        r3 = r3 + i6 *
                      (mux1 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                       mux2 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                       mux3 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                       mux4 * (u(3, i + 2, j, k) - u(3, i, j, k))) *
                      istry;

        // qq derivative (w) (w-eq)
        // 43 ops, tot=283
        cof1 = (mu(i, j - 2, k)) * met(1, i, j - 2, k) * met(1, i, j - 2, k) *
               stry(j - 2);
        cof2 = (mu(i, j - 1, k)) * met(1, i, j - 1, k) * met(1, i, j - 1, k) *
               stry(j - 1);
        cof3 = (mu(i, j, k)) * met(1, i, j, k) * met(1, i, j, k) * stry(j);
        cof4 = (mu(i, j + 1, k)) * met(1, i, j + 1, k) * met(1, i, j + 1, k) *
               stry(j + 1);
        cof5 = (mu(i, j + 2, k)) * met(1, i, j + 2, k) * met(1, i, j + 2, k) *
               stry(j + 2);
        mux1 = cof2 - tf * (cof3 + cof1);
        mux2 = cof1 + cof4 + 3 * (cof3 + cof2);
        mux3 = cof2 + cof5 + 3 * (cof4 + cof3);
        mux4 = cof4 - tf * (cof3 + cof5);

        r3 = r3 + i6 *
                      (mux1 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                       mux2 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                       mux3 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                       mux4 * (u(3, i, j + 2, k) - u(3, i, j, k))) *
                      istrx;
      }
      // All rr-derivatives at once
      // averaging the coefficient
      // 54*8*8+25*8 = 3656 ops, tot=3939
      float_sw4 mucofu2, mucofuv, mucofuw, mucofvw, mucofv2, mucofw2;
      for (int q = nk - 7; q <= nk; q++) {
        mucofu2 = 0;
        mucofuv = 0;
        mucofuw = 0;
        mucofvw = 0;
        mucofv2 = 0;
        mucofw2 = 0;
        for (int m = nk - 7; m <= nk; m++) {
          mucofu2 += acof_no_gp(nk - k + 1, nk - q + 1, nk - m + 1) *
                     ((2 * mu(i, j, m) + la(i, j, m)) * met(2, i, j, m) *
                          strx(i) * met(2, i, j, m) * strx(i) +
                      mu(i, j, m) * (met(3, i, j, m) * stry(j) *
                                         met(3, i, j, m) * stry(j) +
                                     met(4, i, j, m) * met(4, i, j, m)));
          mucofv2 += acof_no_gp(nk - k + 1, nk - q + 1, nk - m + 1) *
                     ((2 * mu(i, j, m) + la(i, j, m)) * met(3, i, j, m) *
                          stry(j) * met(3, i, j, m) * stry(j) +
                      mu(i, j, m) * (met(2, i, j, m) * strx(i) *
                                         met(2, i, j, m) * strx(i) +
                                     met(4, i, j, m) * met(4, i, j, m)));
          mucofw2 +=
              acof_no_gp(nk - k + 1, nk - q + 1, nk - m + 1) *
              ((2 * mu(i, j, m) + la(i, j, m)) * met(4, i, j, m) *
                   met(4, i, j, m) +
               mu(i, j, m) *
                   (met(2, i, j, m) * strx(i) * met(2, i, j, m) * strx(i) +
                    met(3, i, j, m) * stry(j) * met(3, i, j, m) * stry(j)));
          mucofuv += acof_no_gp(nk - k + 1, nk - q + 1, nk - m + 1) *
                     (mu(i, j, m) + la(i, j, m)) * met(2, i, j, m) *
                     met(3, i, j, m);
          mucofuw += acof_no_gp(nk - k + 1, nk - q + 1, nk - m + 1) *
                     (mu(i, j, m) + la(i, j, m)) * met(2, i, j, m) *
                     met(4, i, j, m);
          mucofvw += acof_no_gp(nk - k + 1, nk - q + 1, nk - m + 1) *
                     (mu(i, j, m) + la(i, j, m)) * met(3, i, j, m) *
                     met(4, i, j, m);
        }

        // Computing the second derivative,
        if (N == 1)
          r1 += istrxy * mucofu2 * u(1, i, j, q) + mucofuv * u(2, i, j, q) +
                istry * mucofuw * u(3, i, j, q);
        if (N == 2)
          r2 += mucofuv * u(1, i, j, q) + istrxy * mucofv2 * u(2, i, j, q) +
                istrx * mucofvw * u(3, i, j, q);
        if (N == 3)
          r3 += istry * mucofuw * u(1, i, j, q) +
                istrx * mucofvw * u(2, i, j, q) +
                istrxy * mucofw2 * u(3, i, j, q);
      }

      // Ghost point values, only nonzero for k=nk.
      // 72 ops., tot=4011
      mucofu2 = ghcof_no_gp(nk - k + 1) *
                ((2 * mu(i, j, nk) + la(i, j, nk)) * met(2, i, j, nk) *
                     strx(i) * met(2, i, j, nk) * strx(i) +
                 mu(i, j, nk) *
                     (met(3, i, j, nk) * stry(j) * met(3, i, j, nk) * stry(j) +
                      met(4, i, j, nk) * met(4, i, j, nk)));
      mucofv2 = ghcof_no_gp(nk - k + 1) *
                ((2 * mu(i, j, nk) + la(i, j, nk)) * met(3, i, j, nk) *
                     stry(j) * met(3, i, j, nk) * stry(j) +
                 mu(i, j, nk) *
                     (met(2, i, j, nk) * strx(i) * met(2, i, j, nk) * strx(i) +
                      met(4, i, j, nk) * met(4, i, j, nk)));
      mucofw2 = ghcof_no_gp(nk - k + 1) *
                ((2 * mu(i, j, nk) + la(i, j, nk)) * met(4, i, j, nk) *
                     met(4, i, j, nk) +
                 mu(i, j, nk) *
                     (met(2, i, j, nk) * strx(i) * met(2, i, j, nk) * strx(i) +
                      met(3, i, j, nk) * stry(j) * met(3, i, j, nk) * stry(j)));
      mucofuv = ghcof_no_gp(nk - k + 1) * (mu(i, j, nk) + la(i, j, nk)) *
                met(2, i, j, nk) * met(3, i, j, nk);
      mucofuw = ghcof_no_gp(nk - k + 1) * (mu(i, j, nk) + la(i, j, nk)) *
                met(2, i, j, nk) * met(4, i, j, nk);
      mucofvw = ghcof_no_gp(nk - k + 1) * (mu(i, j, nk) + la(i, j, nk)) *
                met(3, i, j, nk) * met(4, i, j, nk);
      if (N == 1)
        r1 += istrxy * mucofu2 * u(1, i, j, nk + 1) +
              mucofuv * u(2, i, j, nk + 1) +
              istry * mucofuw * u(3, i, j, nk + 1);
      if (N == 2)
        r2 += mucofuv * u(1, i, j, nk + 1) +
              istrxy * mucofv2 * u(2, i, j, nk + 1) +
              istrx * mucofvw * u(3, i, j, nk + 1);
      if (N == 3)
        r3 += istry * mucofuw * u(1, i, j, nk + 1) +
              istrx * mucofvw * u(2, i, j, nk + 1) +
              istrxy * mucofw2 * u(3, i, j, nk + 1);

      // pq-derivatives (u-eq)
      // 38 ops., tot=4049
      if (N == 1) {
        r1 +=
            c2 * (mu(i, j + 2, k) * met(1, i, j + 2, k) * met(1, i, j + 2, k) *
                      (c2 * (u(2, i + 2, j + 2, k) - u(2, i - 2, j + 2, k)) +
                       c1 * (u(2, i + 1, j + 2, k) - u(2, i - 1, j + 2, k))) -
                  mu(i, j - 2, k) * met(1, i, j - 2, k) * met(1, i, j - 2, k) *
                      (c2 * (u(2, i + 2, j - 2, k) - u(2, i - 2, j - 2, k)) +
                       c1 * (u(2, i + 1, j - 2, k) - u(2, i - 1, j - 2, k)))) +
            c1 * (mu(i, j + 1, k) * met(1, i, j + 1, k) * met(1, i, j + 1, k) *
                      (c2 * (u(2, i + 2, j + 1, k) - u(2, i - 2, j + 1, k)) +
                       c1 * (u(2, i + 1, j + 1, k) - u(2, i - 1, j + 1, k))) -
                  mu(i, j - 1, k) * met(1, i, j - 1, k) * met(1, i, j - 1, k) *
                      (c2 * (u(2, i + 2, j - 1, k) - u(2, i - 2, j - 1, k)) +
                       c1 * (u(2, i + 1, j - 1, k) - u(2, i - 1, j - 1, k))));

        // qp-derivatives (u-eq)
        // 38 ops. tot=4087
        r1 +=
            c2 * (la(i + 2, j, k) * met(1, i + 2, j, k) * met(1, i + 2, j, k) *
                      (c2 * (u(2, i + 2, j + 2, k) - u(2, i + 2, j - 2, k)) +
                       c1 * (u(2, i + 2, j + 1, k) - u(2, i + 2, j - 1, k))) -
                  la(i - 2, j, k) * met(1, i - 2, j, k) * met(1, i - 2, j, k) *
                      (c2 * (u(2, i - 2, j + 2, k) - u(2, i - 2, j - 2, k)) +
                       c1 * (u(2, i - 2, j + 1, k) - u(2, i - 2, j - 1, k)))) +
            c1 * (la(i + 1, j, k) * met(1, i + 1, j, k) * met(1, i + 1, j, k) *
                      (c2 * (u(2, i + 1, j + 2, k) - u(2, i + 1, j - 2, k)) +
                       c1 * (u(2, i + 1, j + 1, k) - u(2, i + 1, j - 1, k))) -
                  la(i - 1, j, k) * met(1, i - 1, j, k) * met(1, i - 1, j, k) *
                      (c2 * (u(2, i - 1, j + 2, k) - u(2, i - 1, j - 2, k)) +
                       c1 * (u(2, i - 1, j + 1, k) - u(2, i - 1, j - 1, k))));
      }
      // pq-derivatives (v-eq)
      // 38 ops. , tot=4125
      if (N == 2) {
        r2 +=
            c2 * (la(i, j + 2, k) * met(1, i, j + 2, k) * met(1, i, j + 2, k) *
                      (c2 * (u(1, i + 2, j + 2, k) - u(1, i - 2, j + 2, k)) +
                       c1 * (u(1, i + 1, j + 2, k) - u(1, i - 1, j + 2, k))) -
                  la(i, j - 2, k) * met(1, i, j - 2, k) * met(1, i, j - 2, k) *
                      (c2 * (u(1, i + 2, j - 2, k) - u(1, i - 2, j - 2, k)) +
                       c1 * (u(1, i + 1, j - 2, k) - u(1, i - 1, j - 2, k)))) +
            c1 * (la(i, j + 1, k) * met(1, i, j + 1, k) * met(1, i, j + 1, k) *
                      (c2 * (u(1, i + 2, j + 1, k) - u(1, i - 2, j + 1, k)) +
                       c1 * (u(1, i + 1, j + 1, k) - u(1, i - 1, j + 1, k))) -
                  la(i, j - 1, k) * met(1, i, j - 1, k) * met(1, i, j - 1, k) *
                      (c2 * (u(1, i + 2, j - 1, k) - u(1, i - 2, j - 1, k)) +
                       c1 * (u(1, i + 1, j - 1, k) - u(1, i - 1, j - 1, k))));

        //* qp-derivatives (v-eq)
        // 38 ops., tot=4163
        r2 +=
            c2 * (mu(i + 2, j, k) * met(1, i + 2, j, k) * met(1, i + 2, j, k) *
                      (c2 * (u(1, i + 2, j + 2, k) - u(1, i + 2, j - 2, k)) +
                       c1 * (u(1, i + 2, j + 1, k) - u(1, i + 2, j - 1, k))) -
                  mu(i - 2, j, k) * met(1, i - 2, j, k) * met(1, i - 2, j, k) *
                      (c2 * (u(1, i - 2, j + 2, k) - u(1, i - 2, j - 2, k)) +
                       c1 * (u(1, i - 2, j + 1, k) - u(1, i - 2, j - 1, k)))) +
            c1 * (mu(i + 1, j, k) * met(1, i + 1, j, k) * met(1, i + 1, j, k) *
                      (c2 * (u(1, i + 1, j + 2, k) - u(1, i + 1, j - 2, k)) +
                       c1 * (u(1, i + 1, j + 1, k) - u(1, i + 1, j - 1, k))) -
                  mu(i - 1, j, k) * met(1, i - 1, j, k) * met(1, i - 1, j, k) *
                      (c2 * (u(1, i - 1, j + 2, k) - u(1, i - 1, j - 2, k)) +
                       c1 * (u(1, i - 1, j + 1, k) - u(1, i - 1, j - 1, k))));
      }
      // rp - derivatives
      // 24*8 = 192 ops, tot=4355
      float_sw4 dudrm2 = 0, dudrm1 = 0, dudrp1 = 0, dudrp2 = 0;
      float_sw4 dvdrm2 = 0, dvdrm1 = 0, dvdrp1 = 0, dvdrp2 = 0;
      float_sw4 dwdrm2 = 0, dwdrm1 = 0, dwdrp1 = 0, dwdrp2 = 0;
      for (int q = nk - 7; q <= nk; q++) {
        dudrm2 -= bope(nk - k + 1, nk - q + 1) * u(1, i - 2, j, q);
        dvdrm2 -= bope(nk - k + 1, nk - q + 1) * u(2, i - 2, j, q);
        dwdrm2 -= bope(nk - k + 1, nk - q + 1) * u(3, i - 2, j, q);
        dudrm1 -= bope(nk - k + 1, nk - q + 1) * u(1, i - 1, j, q);
        dvdrm1 -= bope(nk - k + 1, nk - q + 1) * u(2, i - 1, j, q);
        dwdrm1 -= bope(nk - k + 1, nk - q + 1) * u(3, i - 1, j, q);
        dudrp2 -= bope(nk - k + 1, nk - q + 1) * u(1, i + 2, j, q);
        dvdrp2 -= bope(nk - k + 1, nk - q + 1) * u(2, i + 2, j, q);
        dwdrp2 -= bope(nk - k + 1, nk - q + 1) * u(3, i + 2, j, q);
        dudrp1 -= bope(nk - k + 1, nk - q + 1) * u(1, i + 1, j, q);
        dvdrp1 -= bope(nk - k + 1, nk - q + 1) * u(2, i + 1, j, q);
        dwdrp1 -= bope(nk - k + 1, nk - q + 1) * u(3, i + 1, j, q);
      }

      // rp derivatives (u-eq)
      // 67 ops, tot=4422
      if (N == 1)
        r1 += (c2 * ((2 * mu(i + 2, j, k) + la(i + 2, j, k)) *
                         met(2, i + 2, j, k) * met(1, i + 2, j, k) *
                         strx(i + 2) * dudrp2 +
                     la(i + 2, j, k) * met(3, i + 2, j, k) *
                         met(1, i + 2, j, k) * dvdrp2 * stry(j) +
                     la(i + 2, j, k) * met(4, i + 2, j, k) *
                         met(1, i + 2, j, k) * dwdrp2 -
                     ((2 * mu(i - 2, j, k) + la(i - 2, j, k)) *
                          met(2, i - 2, j, k) * met(1, i - 2, j, k) *
                          strx(i - 2) * dudrm2 +
                      la(i - 2, j, k) * met(3, i - 2, j, k) *
                          met(1, i - 2, j, k) * dvdrm2 * stry(j) +
                      la(i - 2, j, k) * met(4, i - 2, j, k) *
                          met(1, i - 2, j, k) * dwdrm2)) +
               c1 * ((2 * mu(i + 1, j, k) + la(i + 1, j, k)) *
                         met(2, i + 1, j, k) * met(1, i + 1, j, k) *
                         strx(i + 1) * dudrp1 +
                     la(i + 1, j, k) * met(3, i + 1, j, k) *
                         met(1, i + 1, j, k) * dvdrp1 * stry(j) +
                     la(i + 1, j, k) * met(4, i + 1, j, k) *
                         met(1, i + 1, j, k) * dwdrp1 -
                     ((2 * mu(i - 1, j, k) + la(i - 1, j, k)) *
                          met(2, i - 1, j, k) * met(1, i - 1, j, k) *
                          strx(i - 1) * dudrm1 +
                      la(i - 1, j, k) * met(3, i - 1, j, k) *
                          met(1, i - 1, j, k) * dvdrm1 * stry(j) +
                      la(i - 1, j, k) * met(4, i - 1, j, k) *
                          met(1, i - 1, j, k) * dwdrm1))) *
              istry;

      // rp derivatives (v-eq)
      // 42 ops, tot=4464
      if (N == 2)
        r2 += c2 * (mu(i + 2, j, k) * met(3, i + 2, j, k) *
                        met(1, i + 2, j, k) * dudrp2 +
                    mu(i + 2, j, k) * met(2, i + 2, j, k) *
                        met(1, i + 2, j, k) * dvdrp2 * strx(i + 2) * istry -
                    (mu(i - 2, j, k) * met(3, i - 2, j, k) *
                         met(1, i - 2, j, k) * dudrm2 +
                     mu(i - 2, j, k) * met(2, i - 2, j, k) *
                         met(1, i - 2, j, k) * dvdrm2 * strx(i - 2) * istry)) +
              c1 * (mu(i + 1, j, k) * met(3, i + 1, j, k) *
                        met(1, i + 1, j, k) * dudrp1 +
                    mu(i + 1, j, k) * met(2, i + 1, j, k) *
                        met(1, i + 1, j, k) * dvdrp1 * strx(i + 1) * istry -
                    (mu(i - 1, j, k) * met(3, i - 1, j, k) *
                         met(1, i - 1, j, k) * dudrm1 +
                     mu(i - 1, j, k) * met(2, i - 1, j, k) *
                         met(1, i - 1, j, k) * dvdrm1 * strx(i - 1) * istry));

      // rp derivatives (w-eq)
      // 38 ops, tot=4502
      if (N == 3)
        r3 += istry * (c2 * (mu(i + 2, j, k) * met(4, i + 2, j, k) *
                                 met(1, i + 2, j, k) * dudrp2 +
                             mu(i + 2, j, k) * met(2, i + 2, j, k) *
                                 met(1, i + 2, j, k) * dwdrp2 * strx(i + 2) -
                             (mu(i - 2, j, k) * met(4, i - 2, j, k) *
                                  met(1, i - 2, j, k) * dudrm2 +
                              mu(i - 2, j, k) * met(2, i - 2, j, k) *
                                  met(1, i - 2, j, k) * dwdrm2 * strx(i - 2))) +
                       c1 * (mu(i + 1, j, k) * met(4, i + 1, j, k) *
                                 met(1, i + 1, j, k) * dudrp1 +
                             mu(i + 1, j, k) * met(2, i + 1, j, k) *
                                 met(1, i + 1, j, k) * dwdrp1 * strx(i + 1) -
                             (mu(i - 1, j, k) * met(4, i - 1, j, k) *
                                  met(1, i - 1, j, k) * dudrm1 +
                              mu(i - 1, j, k) * met(2, i - 1, j, k) *
                                  met(1, i - 1, j, k) * dwdrm1 * strx(i - 1))));

      // rq - derivatives
      // 24*8 = 192 ops , tot=4694

      dudrm2 = 0;
      dudrm1 = 0;
      dudrp1 = 0;
      dudrp2 = 0;
      dvdrm2 = 0;
      dvdrm1 = 0;
      dvdrp1 = 0;
      dvdrp2 = 0;
      dwdrm2 = 0;
      dwdrm1 = 0;
      dwdrp1 = 0;
      dwdrp2 = 0;
      for (int q = nk - 7; q <= nk; q++) {
        dudrm2 -= bope(nk - k + 1, nk - q + 1) * u(1, i, j - 2, q);
        dvdrm2 -= bope(nk - k + 1, nk - q + 1) * u(2, i, j - 2, q);
        dwdrm2 -= bope(nk - k + 1, nk - q + 1) * u(3, i, j - 2, q);
        dudrm1 -= bope(nk - k + 1, nk - q + 1) * u(1, i, j - 1, q);
        dvdrm1 -= bope(nk - k + 1, nk - q + 1) * u(2, i, j - 1, q);
        dwdrm1 -= bope(nk - k + 1, nk - q + 1) * u(3, i, j - 1, q);
        dudrp2 -= bope(nk - k + 1, nk - q + 1) * u(1, i, j + 2, q);
        dvdrp2 -= bope(nk - k + 1, nk - q + 1) * u(2, i, j + 2, q);
        dwdrp2 -= bope(nk - k + 1, nk - q + 1) * u(3, i, j + 2, q);
        dudrp1 -= bope(nk - k + 1, nk - q + 1) * u(1, i, j + 1, q);
        dvdrp1 -= bope(nk - k + 1, nk - q + 1) * u(2, i, j + 1, q);
        dwdrp1 -= bope(nk - k + 1, nk - q + 1) * u(3, i, j + 1, q);
      }

      // rq derivatives (u-eq)
      // 42 ops, tot=4736
      if (N == 1)
        r1 += c2 * (mu(i, j + 2, k) * met(3, i, j + 2, k) *
                        met(1, i, j + 2, k) * dudrp2 * stry(j + 2) * istrx +
                    mu(i, j + 2, k) * met(2, i, j + 2, k) *
                        met(1, i, j + 2, k) * dvdrp2 -
                    (mu(i, j - 2, k) * met(3, i, j - 2, k) *
                         met(1, i, j - 2, k) * dudrm2 * stry(j - 2) * istrx +
                     mu(i, j - 2, k) * met(2, i, j - 2, k) *
                         met(1, i, j - 2, k) * dvdrm2)) +
              c1 * (mu(i, j + 1, k) * met(3, i, j + 1, k) *
                        met(1, i, j + 1, k) * dudrp1 * stry(j + 1) * istrx +
                    mu(i, j + 1, k) * met(2, i, j + 1, k) *
                        met(1, i, j + 1, k) * dvdrp1 -
                    (mu(i, j - 1, k) * met(3, i, j - 1, k) *
                         met(1, i, j - 1, k) * dudrm1 * stry(j - 1) * istrx +
                     mu(i, j - 1, k) * met(2, i, j - 1, k) *
                         met(1, i, j - 1, k) * dvdrm1));

      // rq derivatives (v-eq)
      // 70 ops, tot=4806
      if (N == 2)
        r2 += c2 * (la(i, j + 2, k) * met(2, i, j + 2, k) *
                        met(1, i, j + 2, k) * dudrp2 +
                    (2 * mu(i, j + 2, k) + la(i, j + 2, k)) *
                        met(3, i, j + 2, k) * met(1, i, j + 2, k) * dvdrp2 *
                        stry(j + 2) * istrx +
                    la(i, j + 2, k) * met(4, i, j + 2, k) *
                        met(1, i, j + 2, k) * dwdrp2 * istrx -
                    (la(i, j - 2, k) * met(2, i, j - 2, k) *
                         met(1, i, j - 2, k) * dudrm2 +
                     (2 * mu(i, j - 2, k) + la(i, j - 2, k)) *
                         met(3, i, j - 2, k) * met(1, i, j - 2, k) * dvdrm2 *
                         stry(j - 2) * istrx +
                     la(i, j - 2, k) * met(4, i, j - 2, k) *
                         met(1, i, j - 2, k) * dwdrm2 * istrx)) +
              c1 * (la(i, j + 1, k) * met(2, i, j + 1, k) *
                        met(1, i, j + 1, k) * dudrp1 +
                    (2 * mu(i, j + 1, k) + la(i, j + 1, k)) *
                        met(3, i, j + 1, k) * met(1, i, j + 1, k) * dvdrp1 *
                        stry(j + 1) * istrx +
                    la(i, j + 1, k) * met(4, i, j + 1, k) *
                        met(1, i, j + 1, k) * dwdrp1 * istrx -
                    (la(i, j - 1, k) * met(2, i, j - 1, k) *
                         met(1, i, j - 1, k) * dudrm1 +
                     (2 * mu(i, j - 1, k) + la(i, j - 1, k)) *
                         met(3, i, j - 1, k) * met(1, i, j - 1, k) * dvdrm1 *
                         stry(j - 1) * istrx +
                     la(i, j - 1, k) * met(4, i, j - 1, k) *
                         met(1, i, j - 1, k) * dwdrm1 * istrx));

      // rq derivatives (w-eq)
      // 39 ops, tot=4845
      if (N == 3)
        r3 += (c2 * (mu(i, j + 2, k) * met(3, i, j + 2, k) *
                         met(1, i, j + 2, k) * dwdrp2 * stry(j + 2) +
                     mu(i, j + 2, k) * met(4, i, j + 2, k) *
                         met(1, i, j + 2, k) * dvdrp2 -
                     (mu(i, j - 2, k) * met(3, i, j - 2, k) *
                          met(1, i, j - 2, k) * dwdrm2 * stry(j - 2) +
                      mu(i, j - 2, k) * met(4, i, j - 2, k) *
                          met(1, i, j - 2, k) * dvdrm2)) +
               c1 * (mu(i, j + 1, k) * met(3, i, j + 1, k) *
                         met(1, i, j + 1, k) * dwdrp1 * stry(j + 1) +
                     mu(i, j + 1, k) * met(4, i, j + 1, k) *
                         met(1, i, j + 1, k) * dvdrp1 -
                     (mu(i, j - 1, k) * met(3, i, j - 1, k) *
                          met(1, i, j - 1, k) * dwdrm1 * stry(j - 1) +
                      mu(i, j - 1, k) * met(4, i, j - 1, k) *
                          met(1, i, j - 1, k) * dvdrm1))) *
              istrx;

      // pr and qr derivatives at once
      // in loop: 8*(53+53+43) = 1192 ops, tot=6037
      for (int q = nk - 7; q <= nk; q++) {
        // (u-eq)
        // 53 ops
        if (N == 1)
          r1 -= bope(nk - k + 1, nk - q + 1) *
                (
                    // pr
                    (2 * mu(i, j, q) + la(i, j, q)) * met(2, i, j, q) *
                        met(1, i, j, q) *
                        (c2 * (u(1, i + 2, j, q) - u(1, i - 2, j, q)) +
                         c1 * (u(1, i + 1, j, q) - u(1, i - 1, j, q))) *
                        strx(i) * istry +
                    mu(i, j, q) * met(3, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(2, i + 2, j, q) - u(2, i - 2, j, q)) +
                         c1 * (u(2, i + 1, j, q) - u(2, i - 1, j, q))) +
                    mu(i, j, q) * met(4, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(3, i + 2, j, q) - u(3, i - 2, j, q)) +
                         c1 * (u(3, i + 1, j, q) - u(3, i - 1, j, q))) *
                        istry
                    // qr
                    + mu(i, j, q) * met(3, i, j, q) * met(1, i, j, q) *
                          (c2 * (u(1, i, j + 2, q) - u(1, i, j - 2, q)) +
                           c1 * (u(1, i, j + 1, q) - u(1, i, j - 1, q))) *
                          stry(j) * istrx +
                    la(i, j, q) * met(2, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(2, i, j + 2, q) - u(2, i, j - 2, q)) +
                         c1 * (u(2, i, j + 1, q) - u(2, i, j - 1, q))));

        // (v-eq)
        // 53 ops
        if (N == 2)
          r2 -= bope(nk - k + 1, nk - q + 1) *
                (
                    // pr
                    la(i, j, q) * met(3, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(1, i + 2, j, q) - u(1, i - 2, j, q)) +
                         c1 * (u(1, i + 1, j, q) - u(1, i - 1, j, q))) +
                    mu(i, j, q) * met(2, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(2, i + 2, j, q) - u(2, i - 2, j, q)) +
                         c1 * (u(2, i + 1, j, q) - u(2, i - 1, j, q))) *
                        strx(i) * istry
                    // qr
                    + mu(i, j, q) * met(2, i, j, q) * met(1, i, j, q) *
                          (c2 * (u(1, i, j + 2, q) - u(1, i, j - 2, q)) +
                           c1 * (u(1, i, j + 1, q) - u(1, i, j - 1, q))) +
                    (2 * mu(i, j, q) + la(i, j, q)) * met(3, i, j, q) *
                        met(1, i, j, q) *
                        (c2 * (u(2, i, j + 2, q) - u(2, i, j - 2, q)) +
                         c1 * (u(2, i, j + 1, q) - u(2, i, j - 1, q))) *
                        stry(j) * istrx +
                    mu(i, j, q) * met(4, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(3, i, j + 2, q) - u(3, i, j - 2, q)) +
                         c1 * (u(3, i, j + 1, q) - u(3, i, j - 1, q))) *
                        istrx);

        // (w-eq)
        // 43 ops
        if (N == 3)
          r3 -= bope(nk - k + 1, nk - q + 1) *
                (
                    // pr
                    la(i, j, q) * met(4, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(1, i + 2, j, q) - u(1, i - 2, j, q)) +
                         c1 * (u(1, i + 1, j, q) - u(1, i - 1, j, q))) *
                        istry +
                    mu(i, j, q) * met(2, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(3, i + 2, j, q) - u(3, i - 2, j, q)) +
                         c1 * (u(3, i + 1, j, q) - u(3, i - 1, j, q))) *
                        strx(i) * istry
                    // qr
                    + mu(i, j, q) * met(3, i, j, q) * met(1, i, j, q) *
                          (c2 * (u(3, i, j + 2, q) - u(3, i, j - 2, q)) +
                           c1 * (u(3, i, j + 1, q) - u(3, i, j - 1, q))) *
                          stry(j) * istrx +
                    la(i, j, q) * met(4, i, j, q) * met(1, i, j, q) *
                        (c2 * (u(2, i, j + 2, q) - u(2, i, j - 2, q)) +
                         c1 * (u(2, i, j + 1, q) - u(2, i, j - 1, q))) *
                        istrx);
      }

      // 12 ops, tot=6049
      if (N == 1) lu(1, i, j, k) = a1 * lu(1, i, j, k) + sgn * r1 * ijac;
      if (N == 2) lu(2, i, j, k) = a1 * lu(2, i, j, k) + sgn * r2 * ijac;
      if (N == 3) lu(3, i, j, k) = a1 * lu(3, i, j, k) + sgn * r3 * ijac;
    });
  }
  SW4_MARK_END("CURVI::cuvilinear4sgc");
#ifdef PEEKS_GALORE
  SW4_PEEK;
  SYNC_DEVICE;
#endif

  // SYNC_STREAM; // NOW BEING DONE at the end of evalRHS
#undef mu
#undef la
#undef jac
#undef u
#undef lu
#undef met
#undef strx
#undef stry
#undef acof
#undef bope
#undef ghcof
#undef acof_no_gp
#undef ghcof_no_gp
}
