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

#include "EW.h"
#include "sw4.h"

void EW::randomfield3d_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                          sw4_type kfirst, sw4_type klast, sw4_type nig, sw4_type njg, sw4_type nkg,
                          sw4_type gh, float_sw4* __restrict__ a_w,
                          float_sw4* __restrict__ a_wgh, float_sw4 dist,
                          float_sw4 distz, float_sw4 h, sw4_type* randw,
                          float_sw4* __restrict__ a_saverands, sw4_type p, sw4_type pz) {
  sw4_type iseed1 = randw[0];
  sw4_type iseed2 = randw[1];
  sw4_type iseed3 = randw[2];
  float_sw4* _b = new float_sw4[2 * p + 1];
  float_sw4* _bz = new float_sw4[2 * pz + 1];
  const sw4_type ni = ilast - ifirst + 1;
  const sw4_type nij = ni * (jlast - jfirst + 1);
  const sw4_type base = -ifirst - ni * jfirst - nij * kfirst;
  const sw4_type nisr = ni + 2 * p;
  const sw4_type nijsr = nisr * (jlast - jfirst + 1 + 2 * p);
  const sw4_type basesr = -ifirst + p - nisr * (jfirst - p) - nijsr * (1 - pz);
#define b(i) _b[(i) + p]
#define bz(i) _bz[(i) + pz]
#define w(i, j, k) a_w[base + (i) + ni * (j) + nij * (k)]
#define wgh(i, j, k) a_wgh[base + (i) + ni * (j) + nij * (k)]
#define saverands(i, j, k) a_saverands[basesr + (i) + nisr * (j) + nijsr * (k)]

  for (sw4_type k = -p; k <= p; k++) b(k) = exp(-k * k * h * h / (2 * dist * dist));
  for (sw4_type k = -pz; k <= pz; k++)
    bz(k) = exp(-k * k * h * h / (2 * distz * distz));
  size_t npts = static_cast<size_t>(ilast - ifirst + 1) * (jlast - jfirst + 1) *
                (klast - kfirst + 1);
#pragma omp parallel for
  for (sw4_type i = 0; i < npts; i++) a_w[i] = a_wgh[i] = 0;

  // wgh is masw4_typeained because of boundary effects. The stencil is cut
  // at boundaries and will use fewer points there.

  // Can not use OpenMP for this loop, since it would unsync the random number
  // generator
  for (sw4_type k = 1 - pz; k <= nkg + gh; k++)
    for (sw4_type j = 1 - gh; j <= njg + gh; j++)
#pragma ivdep
#pragma simd
      for (sw4_type i = 1 - gh; i <= nig + gh; i++) {
        // Random number generator, expanded sw4_typeo loop
        sw4_type krand = iseed1 / 206;
        iseed1 = 157 * (iseed1 - krand * 206) - krand * 21;
        if (iseed1 < 0) iseed1 = iseed1 + 32363;
        krand = iseed2 / 217;
        iseed2 = 146 * (iseed2 - krand * 217) - krand * 45;
        if (iseed2 < 0) iseed2 = iseed2 + 31727;
        krand = iseed3 / 222;
        iseed3 = 142 * (iseed3 - krand * 222) - krand * 133;
        if (iseed3 < 0) iseed3 = iseed3 + 31657;
        sw4_type iz = iseed1 - iseed2;
        if (iz > 706) iz = iz - 32362;
        iz = iz + iseed3;
        if (iz < 1) iz = iz + 32362;
        // First test if loop over stencil is necessary at all
        if (i + p >= ifirst && i - p <= ilast && j + p >= jfirst &&
            j - p <= jlast && k + pz >= kfirst && k - pz <= klast) {
          // Compute the random number in [-1,1], and loop over stencil
          float randno = 2 * static_cast<float>(iz) * 3.0899e-5 - 1;
          if (k <= 1 + pz) saverands(i, j, k) = randno;
          for (sw4_type kk = -pz; kk <= pz; kk++)
            if (k - kk >= kfirst && k - kk <= klast) {
              for (sw4_type jj = -p; jj <= p; jj++)
                if (j - jj >= jfirst && j - jj <= jlast) {
                  for (sw4_type ii = -p; ii <= p; ii++)
                    if (i - ii >= ifirst && i - ii <= ilast) {
                      w(i - ii, j - jj, k - kk) +=
                          b(ii) * b(jj) * bz(kk) * randno;
                      wgh(i - ii, j - jj, k - kk) +=
                          b(ii) * b(jj) * bz(kk) * b(ii) * b(jj) * bz(kk);
                    }
                }
            }
        }
      }
  // Save state of random number generator
  randw[0] = iseed1;
  randw[1] = iseed2;
  randw[2] = iseed3;

  // Divide weights
#pragma omp parallel for
  for (sw4_type i = 0; i < npts; i++) a_w[i] /= sqrt(a_wgh[i]);

  //   for( sw4_type k=kfirst ; k <= klast ; k++ )
  //      for( sw4_type j=jfirst ; j <= jlast ; j++ )
  //	 for( sw4_type i=ifirst ; i <= ilast ; i++ )
  //	    w(i,j,k) = w(i,j,k)/sqrt(wgh(i,j,k));

  delete[] _b;
  delete[] _bz;
#undef b
#undef bz
#undef w
#undef wz
#undef saverands
}

//-----------------------------------------------------------------------
void EW::randomfield3dc_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                           sw4_type kfirst, sw4_type klast, sw4_type nig, sw4_type njg, sw4_type nkg,
                           sw4_type gh, float_sw4* __restrict__ a_w,
                           float_sw4* __restrict__ a_wgh, float_sw4 dist,
                           float_sw4 distz, float_sw4 h,
                           float_sw4* __restrict__ a_z, sw4_type* randw,
                           float_sw4* __restrict__ a_saverands, sw4_type p, sw4_type pz) {
  sw4_type iseed1 = randw[0];
  sw4_type iseed2 = randw[1];
  sw4_type iseed3 = randw[2];
  float_sw4* _b = new float_sw4[2 * p + 1];
  const sw4_type ni = ilast - ifirst + 1;
  const sw4_type nij = ni * (jlast - jfirst + 1);
  const sw4_type base = -ifirst - ni * jfirst - nij * kfirst;
  const sw4_type nisr = ni + 2 * p;
  const sw4_type nijsr = nisr * (jlast - jfirst + 1 + 2 * p);
  const sw4_type basesr = -ifirst + p - nisr * (jfirst - p) - nijsr * (nkg - pz);
#define b(i) _b[(i) + p]
#define w(i, j, k) a_w[base + (i) + ni * (j) + nij * (k)]
#define z(i, j, k) a_z[base + (i) + ni * (j) + nij * (k)]
#define wgh(i, j, k) a_wgh[base + (i) + ni * (j) + nij * (k)]
#define saverands(i, j, k) a_saverands[basesr + (i) + nisr * (j) + nijsr * (k)]

  for (sw4_type k = -p; k <= p; k++) b(k) = exp(-k * k * h * h / (2 * dist * dist));
  float_sw4 i2dz = 1 / (2 * distz * distz);
  size_t npts = static_cast<size_t>(ilast - ifirst + 1) * (jlast - jfirst + 1) *
                (klast - kfirst + 1);

#pragma omp parallel for
  for (sw4_type i = 0; i < npts; i++) a_w[i] = a_wgh[i] = 0;

  // wgh is masw4_typeained because of boundary effects. The stencil is cut
  // at boundaries and will use fewer points there.

  // Can not use OpenMP for this loop, since it would unsync the random number
  // generator
  for (sw4_type k = 1 - pz; k <= nkg + gh; k++)
    for (sw4_type j = 1 - gh; j <= njg + gh; j++)
#pragma ivdep
#pragma simd
      for (sw4_type i = 1 - gh; i <= nig + gh; i++) {
        // Random number generator, expanded sw4_typeo loop
        sw4_type krand = iseed1 / 206;
        iseed1 = 157 * (iseed1 - krand * 206) - krand * 21;
        if (iseed1 < 0) iseed1 = iseed1 + 32363;
        krand = iseed2 / 217;
        iseed2 = 146 * (iseed2 - krand * 217) - krand * 45;
        if (iseed2 < 0) iseed2 = iseed2 + 31727;
        krand = iseed3 / 222;
        iseed3 = 142 * (iseed3 - krand * 222) - krand * 133;
        if (iseed3 < 0) iseed3 = iseed3 + 31657;
        sw4_type iz = iseed1 - iseed2;
        if (iz > 706) iz = iz - 32362;
        iz = iz + iseed3;
        if (iz < 1) iz = iz + 32362;
        // First test if loop over stencil is necessary at all
        if (i + p >= ifirst && i - p <= ilast && j + p >= jfirst &&
            j - p <= jlast && k + pz >= kfirst && k - pz <= klast) {
          // Compute the random number in [-1,1], and loop over stencil
          float randno;
          if (k <= nkg - pz)
            randno = 2 * static_cast<float>(iz) * 3.0899e-5 - 1;
          else
            randno = saverands(i, j, k);
          for (sw4_type kk = -pz; kk <= pz; kk++)
            if (k - kk >= kfirst && k - kk <= klast) {
              for (sw4_type jj = -p; jj <= p; jj++)
                if (j - jj >= jfirst && j - jj <= jlast) {
                  for (sw4_type ii = -p; ii <= p; ii++)
                    if (i - ii >= ifirst && i - ii <= ilast) {
                      float_sw4 bz = exp(
                          -(z(i - ii, j - jj, k - kk) - z(i - ii, j - jj, k)) *
                          (z(i - ii, j - jj, k - kk) - z(i - ii, j - jj, k)) *
                          i2dz);
                      w(i - ii, j - jj, k - kk) += b(ii) * b(jj) * bz * randno;
                      wgh(i - ii, j - jj, k - kk) +=
                          b(ii) * b(ii) * b(jj) * b(jj) * bz * bz;
                    }
                }
            }
        }
      }
  // Save and return state of random number generator
  randw[0] = iseed1;
  randw[1] = iseed2;
  randw[2] = iseed3;

// Divide weights
#pragma omp parallel for
  for (sw4_type i = 0; i < npts; i++) a_w[i] /= sqrt(a_wgh[i]);
  delete[] _b;
}

//-----------------------------------------------------------------------
void EW::perturbvelocity_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                            sw4_type kfirst, sw4_type klast, float_sw4* __restrict__ a_vs,
                            float_sw4* __restrict__ a_vp,
                            float_sw4* __restrict__ a_per, float_sw4 amp,
                            float_sw4 grad, float_sw4 zmin, float_sw4 h,
                            float_sw4 plimit) {
  size_t ind = 0;
#pragma omp parallel for
  for (sw4_type k = kfirst; k <= klast; k++) {
    float_sw4 A = amp + grad * (zmin + (k - 1) * h);
    for (sw4_type j = jfirst; j <= jlast; j++)
#pragma ivdep
#pragma simd
      for (sw4_type i = ifirst; i <= ilast; i++) {
        float_sw4 perijk = a_per[ind] > -plimit ? a_per[ind] : -plimit;
        perijk = perijk < plimit ? perijk : plimit;
        a_vs[ind] *= 1 + A * perijk;
        a_vp[ind] *= 1 + A * perijk;
        ind++;
      }
  }
}

//-----------------------------------------------------------------------
void EW::perturbvelocityc_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                             sw4_type kfirst, sw4_type klast,
                             float_sw4* __restrict__ a_vs,
                             float_sw4* __restrict__ a_vp,
                             float_sw4* __restrict__ a_per, float_sw4 amp,
                             float_sw4 grad, float_sw4* a_z, float_sw4 plimit) {
  size_t ind = 0;
#pragma omp parallel for
  for (sw4_type k = kfirst; k <= klast; k++)
    for (sw4_type j = jfirst; j <= jlast; j++)
#pragma ivdep
#pragma simd
      for (sw4_type i = ifirst; i <= ilast; i++) {
        float_sw4 A = amp + grad * a_z[ind];
        float_sw4 perijk = a_per[ind] > -plimit ? a_per[ind] : -plimit;
        perijk = perijk < plimit ? perijk : plimit;
        a_vs[ind] *= 1 + A * perijk;
        a_vp[ind] *= 1 + A * perijk;
        ind++;
      }
}
