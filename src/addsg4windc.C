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
//-----------------------------------------------------------------------
// Adds 4th order artificial disssipation for super-grid damping layers
//
// Windowed version used with mesh refinement
//
//-----------------------------------------------------------------------
#include <sys/types.h>
#include "sw4.h"

// extern "C" {

void addsg4wind_ci(float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
                   float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
                   float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
                   float_sw4* __restrict__ a_dcz,
                   float_sw4* __restrict__ a_strx,
                   float_sw4* __restrict__ a_stry,
                   float_sw4* __restrict__ a_strz,
                   float_sw4* __restrict__ a_cox, float_sw4* __restrict__ a_coy,
                   float_sw4* __restrict__ a_coz, int ifirst, int ilast,
                   int jfirst, int jlast, int kfirst, int klast, float_sw4 beta,
                   int kupb, int kupe, int kwindb, int kwinde) {
  //***********************************************************************
  //*** Version with correct density scaling and supergrid stretching.
  //*** cox, coy, coz are corner factors that reduce the damping near edges and
  //corners
  //***
  //***********************************************************************

  //  real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  //  real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  //// up can have different size in k-index
  //  real*8 up(3,ifirst:ilast,jfirst:jlast,kupb:kupe)
  //  real*8  rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
  //  real*8 dcx(ifirst:ilast), strx(ifirst:ilast), cox(ifirst:ilast)
  //  real*8 dcy(jfirst:jlast), stry(jfirst:jlast), coy(jfirst:jlast)
  //  real*8 dcz(kfirst:klast), strz(kfirst:klast), coz(kfirst:klast)

  // this routine uses un-divided differences in x and t
  if (beta == 0) return;

  float_sw4 coeff = beta;
  // beta is the supergrid damping coefficient as entered in the input file
  //
  // add in the SG damping
  //
  // There are enough ghost points to always use the interior formula
  // But only in (i,j) because the k-window may be near a refinement bndry
  //
  // the corner tapering is applied by replacing
  // strx -> strx*coy(j)*coz(k)
  // stry -> stry*cox(i)*coz(k)
  // strz -> strz*cox(i)*coy(j)
  //
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t nijkp = nij * (kupe - kupb + 1);
  const int base = -ifirst - ni * jfirst - nij * kfirst;
  //   const int basep= -ifirst-ni*jfirst-nij*kupb;
  const int base3 = -ifirst - ni * jfirst - nij * kfirst - nijk;
  const int base3p = -ifirst - ni * jfirst - nij * kupb - nijkp;

#define rho(i, j, k) a_rho[base + (i) + ni * (j) + nij * (k)]
#define u(c, i, j, k) a_u[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define um(c, i, j, k) a_um[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define up(c, i, j, k) a_up[base3p + (i) + ni * (j) + nij * (k) + nijkp * (c)]
#define dcx(i) a_dcx[i - ifirst]
#define strx(i) a_strx[i - ifirst]
#define cox(i) a_cox[i - ifirst]
#define dcy(j) a_dcy[j - jfirst]
#define stry(j) a_stry[j - jfirst]
#define coy(j) a_coy[j - jfirst]
#define dcz(k) a_dcz[k - kfirst]
#define strz(k) a_strz[k - kfirst]
#define coz(k) a_coz[k - kfirst]

  for (int c = 1; c <= 3; c++)
    for (int k = kwindb; k <= kwinde; k++)
#pragma omp parallel for
      for (int j = jfirst + 2; j <= jlast - 2; j++)
#pragma ivdep
        //#pragma simd
        for (int i = ifirst + 2; i <= ilast - 2; i++) {
          up(c, i, j, k) -=
              coeff *
              (
                  // x-differences
                  strx(i) * coy(j) * coz(k) *
                      (rho(i + 1, j, k) * dcx(i + 1) *
                           (u(c, i + 2, j, k) - 2 * u(c, i + 1, j, k) +
                            u(c, i, j, k)) -
                       2 * rho(i, j, k) * dcx(i) *
                           (u(c, i + 1, j, k) - 2 * u(c, i, j, k) +
                            u(c, i - 1, j, k)) +
                       rho(i - 1, j, k) * dcx(i - 1) *
                           (u(c, i, j, k) - 2 * u(c, i - 1, j, k) +
                            u(c, i - 2, j, k)) -
                       rho(i + 1, j, k) * dcx(i + 1) *
                           (um(c, i + 2, j, k) - 2 * um(c, i + 1, j, k) +
                            um(c, i, j, k)) +
                       2 * rho(i, j, k) * dcx(i) *
                           (um(c, i + 1, j, k) - 2 * um(c, i, j, k) +
                            um(c, i - 1, j, k)) -
                       rho(i - 1, j, k) * dcx(i - 1) *
                           (um(c, i, j, k) - 2 * um(c, i - 1, j, k) +
                            um(c, i - 2, j, k)))
                  // y-differences
                  + stry(j) * cox(i) * coz(k) *
                        (rho(i, j + 1, k) * dcy(j + 1) *
                             (u(c, i, j + 2, k) - 2 * u(c, i, j + 1, k) +
                              u(c, i, j, k)) -
                         2 * rho(i, j, k) * dcy(j) *
                             (u(c, i, j + 1, k) - 2 * u(c, i, j, k) +
                              u(c, i, j - 1, k)) +
                         rho(i, j - 1, k) * dcy(j - 1) *
                             (u(c, i, j, k) - 2 * u(c, i, j - 1, k) +
                              u(c, i, j - 2, k)) -
                         rho(i, j + 1, k) * dcy(j + 1) *
                             (um(c, i, j + 2, k) - 2 * um(c, i, j + 1, k) +
                              um(c, i, j, k)) +
                         2 * rho(i, j, k) * dcy(j) *
                             (um(c, i, j + 1, k) - 2 * um(c, i, j, k) +
                              um(c, i, j - 1, k)) -
                         rho(i, j - 1, k) * dcy(j - 1) *
                             (um(c, i, j, k) - 2 * um(c, i, j - 1, k) +
                              um(c, i, j - 2, k)))) /
              rho(i, j, k);
        }
}

//}
