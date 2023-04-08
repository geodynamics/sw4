#include <cmath>

#include "sw4.h"
void ilanisocurv_ci(
    sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst, sw4_type klast, sw4_type nk,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_c,
    float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_lu, sw4_type* onesided,
    float_sw4* __restrict__ a_acof, float_sw4* __restrict__ a_bope,
    float_sw4* __restrict__ a_ghcof, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry, float_sw4* __restrict__ a_strz) {
  //  Assumes that metric terms have been merged sw4_typeo the material tensor a_c
  //  before calling this routine

  const float_sw4 i6 = 1.0 / 6;
  const float_sw4 a1 = 2.0 / 3;
  const float_sw4 a2 = -1.0 / 12;

  const sw4_type ni = ilast - ifirst + 1;
  const sw4_type nij = ni * (jlast - jfirst + 1);
  const sw4_type nijk = nij * (klast - kfirst + 1);
  const sw4_type base = -(ifirst + ni * jfirst + nij * kfirst);
  const sw4_type base3 = base - nijk;
  const sw4_type ifirst0 = ifirst;
  const sw4_type jfirst0 = jfirst;
  const sw4_type kfirst0 = kfirst;

  // Direct reuse of fortran code by these macro definitions:
#define c(m, i, j, k) a_c[base3 + (i) + ni * (j) + nij * (k) + nijk * (m)]
#define jac(i, j, k) a_jac[base + (i) + ni * (j) + nij * (k)]
#define u(m, i, j, k) a_u[base3 + (i) + ni * (j) + nij * (k) + nijk * (m)]
#define lu(m, i, j, k) a_lu[base3 + (i) + ni * (j) + nij * (k) + nijk * (m)]
#define strx(i) a_strx[i - ifirst0]
#define stry(j) a_stry[j - jfirst0]
#define strz(k) a_strz[k - kfirst0]
#define acof(i, j, k) a_acof[(i - 1) + 6 * (j - 1) + 48 * (k - 1)]
#define bop(i, j) a_bope[i - 1 + 6 * (j - 1)]
#define ghcof(i) a_ghcof[i - 1]

#pragma omp parallel
  {
    sw4_type kb = kfirst + 2;
    sw4_type ke = klast - 2;
    if (onesided[4] == 1) {
      kb = 7;
#pragma omp for
      for (sw4_type k = 1; k <= 6; k++)
        for (sw4_type j = jfirst + 2; j <= jlast - 2; j++)
#pragma simd
#pragma ivdep
          for (sw4_type i = ifirst + 2; i <= ilast - 2; i++) {
            float_sw4 r1 = 0, r2 = 0, r3 = 0;
            float_sw4 ac1, ac2, ac3, ac4, ac5, ac6;
            float_sw4 dum2, dum1, du, dup1, dup2;
            float_sw4 ijac = 1 / jac(i, j, k);

            float_sw4 cm2 = c(1, i - 1, j, k) * strx(i - 1) -
                            0.75 * (c(1, i, j, k) * strx(i) +
                                    c(1, i - 2, j, k) * strx(i - 2));
            float_sw4 cm1 =
                c(1, i - 2, j, k) * strx(i - 2) +
                c(1, i + 1, j, k) * strx(i + 1) +
                3 * (c(1, i, j, k) * strx(i) + c(1, i - 1, j, k) * strx(i - 1));
            float_sw4 cp1 =
                c(1, i - 1, j, k) * strx(i - 1) +
                c(1, i + 2, j, k) * strx(i + 2) +
                3 * (c(1, i, j, k) * strx(i) + c(1, i + 1, j, k) * strx(i + 1));
            float_sw4 cp2 = c(1, i + 1, j, k) * strx(i + 1) -
                            0.75 * (c(1, i, j, k) * strx(i) +
                                    c(1, i + 2, j, k) * strx(i + 2));

            r1 = r1 + i6 * strx(i) *
                          (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
            cm2 = c(2, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(2, i, j, k) * strx(i) +
                          c(2, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(2, i - 2, j, k) * strx(i - 2) +
                c(2, i + 1, j, k) * strx(i + 1) +
                3 * (c(2, i, j, k) * strx(i) + c(2, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(2, i - 1, j, k) * strx(i - 1) +
                c(2, i + 2, j, k) * strx(i + 2) +
                3 * (c(2, i, j, k) * strx(i) + c(2, i + 1, j, k) * strx(i + 1));
            cp2 = c(2, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(2, i, j, k) * strx(i) +
                          c(2, i + 2, j, k) * strx(i + 2));

            r1 = r1 + i6 * strx(i) *
                          (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
            r2 = r2 + i6 * strx(i) *
                          (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
            cm2 = c(3, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(3, i, j, k) * strx(i) +
                          c(3, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(3, i - 2, j, k) * strx(i - 2) +
                c(3, i + 1, j, k) * strx(i + 1) +
                3 * (c(3, i, j, k) * strx(i) + c(3, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(3, i - 1, j, k) * strx(i - 1) +
                c(3, i + 2, j, k) * strx(i + 2) +
                3 * (c(3, i, j, k) * strx(i) + c(3, i + 1, j, k) * strx(i + 1));
            cp2 = c(3, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(3, i, j, k) * strx(i) +
                          c(3, i + 2, j, k) * strx(i + 2));

            r1 = r1 + i6 * strx(i) *
                          (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
            r3 = r3 + i6 * strx(i) *
                          (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
            cm2 = c(4, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(4, i, j, k) * strx(i) +
                          c(4, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(4, i - 2, j, k) * strx(i - 2) +
                c(4, i + 1, j, k) * strx(i + 1) +
                3 * (c(4, i, j, k) * strx(i) + c(4, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(4, i - 1, j, k) * strx(i - 1) +
                c(4, i + 2, j, k) * strx(i + 2) +
                3 * (c(4, i, j, k) * strx(i) + c(4, i + 1, j, k) * strx(i + 1));
            cp2 = c(4, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(4, i, j, k) * strx(i) +
                          c(4, i + 2, j, k) * strx(i + 2));

            r2 = r2 + i6 * strx(i) *
                          (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
            cm2 = c(5, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(5, i, j, k) * strx(i) +
                          c(5, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(5, i - 2, j, k) * strx(i - 2) +
                c(5, i + 1, j, k) * strx(i + 1) +
                3 * (c(5, i, j, k) * strx(i) + c(5, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(5, i - 1, j, k) * strx(i - 1) +
                c(5, i + 2, j, k) * strx(i + 2) +
                3 * (c(5, i, j, k) * strx(i) + c(5, i + 1, j, k) * strx(i + 1));
            cp2 = c(5, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(5, i, j, k) * strx(i) +
                          c(5, i + 2, j, k) * strx(i + 2));

            r2 = r2 + i6 * strx(i) *
                          (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
            r3 = r3 + i6 * strx(i) *
                          (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
            cm2 = c(6, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(6, i, j, k) * strx(i) +
                          c(6, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(6, i - 2, j, k) * strx(i - 2) +
                c(6, i + 1, j, k) * strx(i + 1) +
                3 * (c(6, i, j, k) * strx(i) + c(6, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(6, i - 1, j, k) * strx(i - 1) +
                c(6, i + 2, j, k) * strx(i + 2) +
                3 * (c(6, i, j, k) * strx(i) + c(6, i + 1, j, k) * strx(i + 1));
            cp2 = c(6, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(6, i, j, k) * strx(i) +
                          c(6, i + 2, j, k) * strx(i + 2));

            r3 = r3 + i6 * strx(i) *
                          (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
            cm2 = c(7, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(7, i, j, k) * stry(j) +
                          c(7, i, j - 2, k) * stry(j - 2));
            cm1 =
                c(7, i, j - 2, k) * stry(j - 2) +
                c(7, i, j + 1, k) * stry(j + 1) +
                3 * (c(7, i, j, k) * stry(j) + c(7, i, j - 1, k) * stry(j - 1));
            cp1 =
                c(7, i, j - 1, k) * stry(j - 1) +
                c(7, i, j + 2, k) * stry(j + 2) +
                3 * (c(7, i, j, k) * stry(j) + c(7, i, j + 1, k) * stry(j + 1));
            cp2 = c(7, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(7, i, j, k) * stry(j) +
                          c(7, i, j + 2, k) * stry(j + 2));

            r1 = r1 + i6 * stry(j) *
                          (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
            cm2 = c(8, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(8, i, j, k) * stry(j) +
                          c(8, i, j - 2, k) * stry(j - 2));
            cm1 =
                c(8, i, j - 2, k) * stry(j - 2) +
                c(8, i, j + 1, k) * stry(j + 1) +
                3 * (c(8, i, j, k) * stry(j) + c(8, i, j - 1, k) * stry(j - 1));
            cp1 =
                c(8, i, j - 1, k) * stry(j - 1) +
                c(8, i, j + 2, k) * stry(j + 2) +
                3 * (c(8, i, j, k) * stry(j) + c(8, i, j + 1, k) * stry(j + 1));
            cp2 = c(8, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(8, i, j, k) * stry(j) +
                          c(8, i, j + 2, k) * stry(j + 2));

            r1 = r1 + i6 * stry(j) *
                          (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
            r2 = r2 + i6 * stry(j) *
                          (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
            cm2 = c(9, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(9, i, j, k) * stry(j) +
                          c(9, i, j - 2, k) * stry(j - 2));
            cm1 =
                c(9, i, j - 2, k) * stry(j - 2) +
                c(9, i, j + 1, k) * stry(j + 1) +
                3 * (c(9, i, j, k) * stry(j) + c(9, i, j - 1, k) * stry(j - 1));
            cp1 =
                c(9, i, j - 1, k) * stry(j - 1) +
                c(9, i, j + 2, k) * stry(j + 2) +
                3 * (c(9, i, j, k) * stry(j) + c(9, i, j + 1, k) * stry(j + 1));
            cp2 = c(9, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(9, i, j, k) * stry(j) +
                          c(9, i, j + 2, k) * stry(j + 2));

            r1 = r1 + i6 * stry(j) *
                          (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
            r3 = r3 + i6 * stry(j) *
                          (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
            cm2 = c(10, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(10, i, j, k) * stry(j) +
                          c(10, i, j - 2, k) * stry(j - 2));
            cm1 = c(10, i, j - 2, k) * stry(j - 2) +
                  c(10, i, j + 1, k) * stry(j + 1) +
                  3 * (c(10, i, j, k) * stry(j) +
                       c(10, i, j - 1, k) * stry(j - 1));
            cp1 = c(10, i, j - 1, k) * stry(j - 1) +
                  c(10, i, j + 2, k) * stry(j + 2) +
                  3 * (c(10, i, j, k) * stry(j) +
                       c(10, i, j + 1, k) * stry(j + 1));
            cp2 = c(10, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(10, i, j, k) * stry(j) +
                          c(10, i, j + 2, k) * stry(j + 2));

            r2 = r2 + i6 * stry(j) *
                          (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
            cm2 = c(11, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(11, i, j, k) * stry(j) +
                          c(11, i, j - 2, k) * stry(j - 2));
            cm1 = c(11, i, j - 2, k) * stry(j - 2) +
                  c(11, i, j + 1, k) * stry(j + 1) +
                  3 * (c(11, i, j, k) * stry(j) +
                       c(11, i, j - 1, k) * stry(j - 1));
            cp1 = c(11, i, j - 1, k) * stry(j - 1) +
                  c(11, i, j + 2, k) * stry(j + 2) +
                  3 * (c(11, i, j, k) * stry(j) +
                       c(11, i, j + 1, k) * stry(j + 1));
            cp2 = c(11, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(11, i, j, k) * stry(j) +
                          c(11, i, j + 2, k) * stry(j + 2));

            r2 = r2 + i6 * stry(j) *
                          (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
            r3 = r3 + i6 * stry(j) *
                          (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
            cm2 = c(12, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(12, i, j, k) * stry(j) +
                          c(12, i, j - 2, k) * stry(j - 2));
            cm1 = c(12, i, j - 2, k) * stry(j - 2) +
                  c(12, i, j + 1, k) * stry(j + 1) +
                  3 * (c(12, i, j, k) * stry(j) +
                       c(12, i, j - 1, k) * stry(j - 1));
            cp1 = c(12, i, j - 1, k) * stry(j - 1) +
                  c(12, i, j + 2, k) * stry(j + 2) +
                  3 * (c(12, i, j, k) * stry(j) +
                       c(12, i, j + 1, k) * stry(j + 1));
            cp2 = c(12, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(12, i, j, k) * stry(j) +
                          c(12, i, j + 2, k) * stry(j + 2));

            r3 = r3 + i6 * stry(j) *
                          (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
            r1 = r1 + ghcof(k) * c(13, i, j, 1) * u(1, i, j, 0) +
                 ghcof(k) * c(14, i, j, 1) * u(2, i, j, 0) +
                 ghcof(k) * c(15, i, j, 1) * u(3, i, j, 0);
            r2 = r2 + ghcof(k) * c(14, i, j, 1) * u(1, i, j, 0) +
                 ghcof(k) * c(16, i, j, 1) * u(2, i, j, 0) +
                 ghcof(k) * c(17, i, j, 1) * u(3, i, j, 0);
            r3 = r3 + ghcof(k) * c(15, i, j, 1) * u(1, i, j, 0) +
                 ghcof(k) * c(17, i, j, 1) * u(2, i, j, 0) +
                 ghcof(k) * c(18, i, j, 1) * u(3, i, j, 0);
            for (sw4_type q = 1; q <= 8; q++) {
              ac1 = 0;
              ac2 = 0;
              ac3 = 0;
              ac4 = 0;
              ac5 = 0;
              ac6 = 0;
              for (sw4_type m = 1; m <= 8; m++) {
                ac1 = ac1 + acof(k, q, m) * c(13, i, j, m);
                ac2 = ac2 + acof(k, q, m) * c(14, i, j, m);
                ac3 = ac3 + acof(k, q, m) * c(15, i, j, m);
                ac4 = ac4 + acof(k, q, m) * c(16, i, j, m);
                ac5 = ac5 + acof(k, q, m) * c(17, i, j, m);
                ac6 = ac6 + acof(k, q, m) * c(18, i, j, m);
              }
              r1 = r1 + ac1 * u(1, i, j, q) + ac2 * u(2, i, j, q) +
                   ac3 * u(3, i, j, q);
              r2 = r2 + ac2 * u(1, i, j, q) + ac4 * u(2, i, j, q) +
                   ac5 * u(3, i, j, q);
              r3 = r3 + ac3 * u(1, i, j, q) + ac5 * u(2, i, j, q) +
                   ac6 * u(3, i, j, q);
            }

            dup1 = a2 * (u(1, i + 1, j + 2, k) - u(1, i + 1, j - 2, k)) +
                   a1 * (u(1, i + 1, j + 1, k) - u(1, i + 1, j - 1, k));
            dum1 = a2 * (u(1, i - 1, j + 2, k) - u(1, i - 1, j - 2, k)) +
                   a1 * (u(1, i - 1, j + 1, k) - u(1, i - 1, j - 1, k));
            dup2 = a2 * (u(1, i + 2, j + 2, k) - u(1, i + 2, j - 2, k)) +
                   a1 * (u(1, i + 2, j + 1, k) - u(1, i + 2, j - 1, k));
            dum2 = a2 * (u(1, i - 2, j + 2, k) - u(1, i - 2, j - 2, k)) +
                   a1 * (u(1, i - 2, j + 1, k) - u(1, i - 2, j - 1, k));
            r1 = r1 + stry(j) * strx(i) *
                          (a1 * (c(19, i + 1, j, k) * dup1 -
                                 c(19, i - 1, j, k) * dum1) +
                           a2 * (c(19, i + 2, j, k) * dup2 -
                                 c(19, i - 2, j, k) * dum2));
            r2 = r2 + stry(j) * strx(i) *
                          (a1 * (c(22, i + 1, j, k) * dup1 -
                                 c(22, i - 1, j, k) * dum1) +
                           a2 * (c(22, i + 2, j, k) * dup2 -
                                 c(22, i - 2, j, k) * dum2));
            r3 = r3 + stry(j) * strx(i) *
                          (a1 * (c(25, i + 1, j, k) * dup1 -
                                 c(25, i - 1, j, k) * dum1) +
                           a2 * (c(25, i + 2, j, k) * dup2 -
                                 c(25, i - 2, j, k) * dum2));

            dup1 = a2 * (u(2, i + 1, j + 2, k) - u(2, i + 1, j - 2, k)) +
                   a1 * (u(2, i + 1, j + 1, k) - u(2, i + 1, j - 1, k));
            dum1 = a2 * (u(2, i - 1, j + 2, k) - u(2, i - 1, j - 2, k)) +
                   a1 * (u(2, i - 1, j + 1, k) - u(2, i - 1, j - 1, k));
            dup2 = a2 * (u(2, i + 2, j + 2, k) - u(2, i + 2, j - 2, k)) +
                   a1 * (u(2, i + 2, j + 1, k) - u(2, i + 2, j - 1, k));
            dum2 = a2 * (u(2, i - 2, j + 2, k) - u(2, i - 2, j - 2, k)) +
                   a1 * (u(2, i - 2, j + 1, k) - u(2, i - 2, j - 1, k));
            r1 = r1 + stry(j) * strx(i) *
                          (a1 * (c(20, i + 1, j, k) * dup1 -
                                 c(20, i - 1, j, k) * dum1) +
                           a2 * (c(20, i + 2, j, k) * dup2 -
                                 c(20, i - 2, j, k) * dum2));
            r2 = r2 + stry(j) * strx(i) *
                          (a1 * (c(23, i + 1, j, k) * dup1 -
                                 c(23, i - 1, j, k) * dum1) +
                           a2 * (c(23, i + 2, j, k) * dup2 -
                                 c(23, i - 2, j, k) * dum2));
            r3 = r3 + stry(j) * strx(i) *
                          (a1 * (c(26, i + 1, j, k) * dup1 -
                                 c(26, i - 1, j, k) * dum1) +
                           a2 * (c(26, i + 2, j, k) * dup2 -
                                 c(26, i - 2, j, k) * dum2));

            dup1 = a2 * (u(3, i + 1, j + 2, k) - u(3, i + 1, j - 2, k)) +
                   a1 * (u(3, i + 1, j + 1, k) - u(3, i + 1, j - 1, k));
            dum1 = a2 * (u(3, i - 1, j + 2, k) - u(3, i - 1, j - 2, k)) +
                   a1 * (u(3, i - 1, j + 1, k) - u(3, i - 1, j - 1, k));
            dup2 = a2 * (u(3, i + 2, j + 2, k) - u(3, i + 2, j - 2, k)) +
                   a1 * (u(3, i + 2, j + 1, k) - u(3, i + 2, j - 1, k));
            dum2 = a2 * (u(3, i - 2, j + 2, k) - u(3, i - 2, j - 2, k)) +
                   a1 * (u(3, i - 2, j + 1, k) - u(3, i - 2, j - 1, k));
            r1 = r1 + stry(j) * strx(i) *
                          (a1 * (c(21, i + 1, j, k) * dup1 -
                                 c(21, i - 1, j, k) * dum1) +
                           a2 * (c(21, i + 2, j, k) * dup2 -
                                 c(21, i - 2, j, k) * dum2));
            r2 = r2 + stry(j) * strx(i) *
                          (a1 * (c(24, i + 1, j, k) * dup1 -
                                 c(24, i - 1, j, k) * dum1) +
                           a2 * (c(24, i + 2, j, k) * dup2 -
                                 c(24, i - 2, j, k) * dum2));
            r3 = r3 + stry(j) * strx(i) *
                          (a1 * (c(27, i + 1, j, k) * dup1 -
                                 c(27, i - 1, j, k) * dum1) +
                           a2 * (c(27, i + 2, j, k) * dup2 -
                                 c(27, i - 2, j, k) * dum2));

            dup1 = a2 * (u(1, i + 2, j + 1, k) - u(1, i - 2, j + 1, k)) +
                   a1 * (u(1, i + 1, j + 1, k) - u(1, i - 1, j + 1, k));
            dum1 = a2 * (u(1, i + 2, j - 1, k) - u(1, i - 2, j - 1, k)) +
                   a1 * (u(1, i + 1, j - 1, k) - u(1, i - 1, j - 1, k));
            dup2 = a2 * (u(1, i + 2, j + 2, k) - u(1, i - 2, j + 2, k)) +
                   a1 * (u(1, i + 1, j + 2, k) - u(1, i - 1, j + 2, k));
            dum2 = a2 * (u(1, i + 2, j - 2, k) - u(1, i - 2, j - 2, k)) +
                   a1 * (u(1, i + 1, j - 2, k) - u(1, i - 1, j - 2, k));
            r1 = r1 + strx(i) * stry(j) *
                          (a1 * (c(19, i, j + 1, k) * dup1 -
                                 c(19, i, j - 1, k) * dum1) +
                           a2 * (c(19, i, j + 2, k) * dup2 -
                                 c(19, i, j - 2, k) * dum2));
            r2 = r2 + strx(i) * stry(j) *
                          (a1 * (c(20, i, j + 1, k) * dup1 -
                                 c(20, i, j - 1, k) * dum1) +
                           a2 * (c(20, i, j + 2, k) * dup2 -
                                 c(20, i, j - 2, k) * dum2));
            r3 = r3 + strx(i) * stry(j) *
                          (a1 * (c(21, i, j + 1, k) * dup1 -
                                 c(21, i, j - 1, k) * dum1) +
                           a2 * (c(21, i, j + 2, k) * dup2 -
                                 c(21, i, j - 2, k) * dum2));

            dup1 = a2 * (u(2, i + 2, j + 1, k) - u(2, i - 2, j + 1, k)) +
                   a1 * (u(2, i + 1, j + 1, k) - u(2, i - 1, j + 1, k));
            dum1 = a2 * (u(2, i + 2, j - 1, k) - u(2, i - 2, j - 1, k)) +
                   a1 * (u(2, i + 1, j - 1, k) - u(2, i - 1, j - 1, k));
            dup2 = a2 * (u(2, i + 2, j + 2, k) - u(2, i - 2, j + 2, k)) +
                   a1 * (u(2, i + 1, j + 2, k) - u(2, i - 1, j + 2, k));
            dum2 = a2 * (u(2, i + 2, j - 2, k) - u(2, i - 2, j - 2, k)) +
                   a1 * (u(2, i + 1, j - 2, k) - u(2, i - 1, j - 2, k));
            r1 = r1 + strx(i) * stry(j) *
                          (a1 * (c(22, i, j + 1, k) * dup1 -
                                 c(22, i, j - 1, k) * dum1) +
                           a2 * (c(22, i, j + 2, k) * dup2 -
                                 c(22, i, j - 2, k) * dum2));
            r2 = r2 + strx(i) * stry(j) *
                          (a1 * (c(23, i, j + 1, k) * dup1 -
                                 c(23, i, j - 1, k) * dum1) +
                           a2 * (c(23, i, j + 2, k) * dup2 -
                                 c(23, i, j - 2, k) * dum2));
            r3 = r3 + strx(i) * stry(j) *
                          (a1 * (c(24, i, j + 1, k) * dup1 -
                                 c(24, i, j - 1, k) * dum1) +
                           a2 * (c(24, i, j + 2, k) * dup2 -
                                 c(24, i, j - 2, k) * dum2));

            dup1 = a2 * (u(3, i + 2, j + 1, k) - u(3, i - 2, j + 1, k)) +
                   a1 * (u(3, i + 1, j + 1, k) - u(3, i - 1, j + 1, k));
            dum1 = a2 * (u(3, i + 2, j - 1, k) - u(3, i - 2, j - 1, k)) +
                   a1 * (u(3, i + 1, j - 1, k) - u(3, i - 1, j - 1, k));
            dup2 = a2 * (u(3, i + 2, j + 2, k) - u(3, i - 2, j + 2, k)) +
                   a1 * (u(3, i + 1, j + 2, k) - u(3, i - 1, j + 2, k));
            dum2 = a2 * (u(3, i + 2, j - 2, k) - u(3, i - 2, j - 2, k)) +
                   a1 * (u(3, i + 1, j - 2, k) - u(3, i - 1, j - 2, k));
            r1 = r1 + strx(i) * stry(j) *
                          (a1 * (c(25, i, j + 1, k) * dup1 -
                                 c(25, i, j - 1, k) * dum1) +
                           a2 * (c(25, i, j + 2, k) * dup2 -
                                 c(25, i, j - 2, k) * dum2));
            r2 = r2 + strx(i) * stry(j) *
                          (a1 * (c(26, i, j + 1, k) * dup1 -
                                 c(26, i, j - 1, k) * dum1) +
                           a2 * (c(26, i, j + 2, k) * dup2 -
                                 c(26, i, j - 2, k) * dum2));
            r3 = r3 + strx(i) * stry(j) *
                          (a1 * (c(27, i, j + 1, k) * dup1 -
                                 c(27, i, j - 1, k) * dum1) +
                           a2 * (c(27, i, j + 2, k) * dup2 -
                                 c(27, i, j - 2, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              dup2 = dup2 + bop(k, q) * u(1, i + 2, j, q);
              dup1 = dup1 + bop(k, q) * u(1, i + 1, j, q);
              dum1 = dum1 + bop(k, q) * u(1, i - 1, j, q);
              dum2 = dum2 + bop(k, q) * u(1, i - 2, j, q);
            }
            r1 = r1 + strx(i) * (a1 * (c(28, i + 1, j, k) * dup1 -
                                       c(28, i - 1, j, k) * dum1) +
                                 a2 * (c(28, i + 2, j, k) * dup2 -
                                       c(28, i - 2, j, k) * dum2));
            r2 = r2 + strx(i) * (a1 * (c(31, i + 1, j, k) * dup1 -
                                       c(31, i - 1, j, k) * dum1) +
                                 a2 * (c(31, i + 2, j, k) * dup2 -
                                       c(31, i - 2, j, k) * dum2));
            r3 = r3 + strx(i) * (a1 * (c(34, i + 1, j, k) * dup1 -
                                       c(34, i - 1, j, k) * dum1) +
                                 a2 * (c(34, i + 2, j, k) * dup2 -
                                       c(34, i - 2, j, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              dup2 = dup2 + bop(k, q) * u(2, i + 2, j, q);
              dup1 = dup1 + bop(k, q) * u(2, i + 1, j, q);
              dum1 = dum1 + bop(k, q) * u(2, i - 1, j, q);
              dum2 = dum2 + bop(k, q) * u(2, i - 2, j, q);
            }
            r1 = r1 + strx(i) * (a1 * (c(29, i + 1, j, k) * dup1 -
                                       c(29, i - 1, j, k) * dum1) +
                                 a2 * (c(29, i + 2, j, k) * dup2 -
                                       c(29, i - 2, j, k) * dum2));
            r2 = r2 + strx(i) * (a1 * (c(32, i + 1, j, k) * dup1 -
                                       c(32, i - 1, j, k) * dum1) +
                                 a2 * (c(32, i + 2, j, k) * dup2 -
                                       c(32, i - 2, j, k) * dum2));
            r3 = r3 + strx(i) * (a1 * (c(35, i + 1, j, k) * dup1 -
                                       c(35, i - 1, j, k) * dum1) +
                                 a2 * (c(35, i + 2, j, k) * dup2 -
                                       c(35, i - 2, j, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              dup2 = dup2 + bop(k, q) * u(3, i + 2, j, q);
              dup1 = dup1 + bop(k, q) * u(3, i + 1, j, q);
              dum1 = dum1 + bop(k, q) * u(3, i - 1, j, q);
              dum2 = dum2 + bop(k, q) * u(3, i - 2, j, q);
            }
            r1 = r1 + strx(i) * (a1 * (c(30, i + 1, j, k) * dup1 -
                                       c(30, i - 1, j, k) * dum1) +
                                 a2 * (c(30, i + 2, j, k) * dup2 -
                                       c(30, i - 2, j, k) * dum2));
            r2 = r2 + strx(i) * (a1 * (c(33, i + 1, j, k) * dup1 -
                                       c(33, i - 1, j, k) * dum1) +
                                 a2 * (c(33, i + 2, j, k) * dup2 -
                                       c(33, i - 2, j, k) * dum2));
            r3 = r3 + strx(i) * (a1 * (c(36, i + 1, j, k) * dup1 -
                                       c(36, i - 1, j, k) * dum1) +
                                 a2 * (c(36, i + 2, j, k) * dup2 -
                                       c(36, i - 2, j, k) * dum2));

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              du = a2 * (u(1, i + 2, j, q) - u(1, i - 2, j, q)) +
                   a1 * (u(1, i + 1, j, q) - u(1, i - 1, j, q));
              ;
              ac1 = ac1 + bop(k, q) * c(28, i, j, q) * du;
              ac2 = ac2 + bop(k, q) * c(29, i, j, q) * du;
              ac3 = ac3 + bop(k, q) * c(30, i, j, q) * du;
            }
            r1 = r1 + strx(i) * ac1;
            r2 = r2 + strx(i) * ac2;
            r3 = r3 + strx(i) * ac3;

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              du = a2 * (u(2, i + 2, j, q) - u(2, i - 2, j, q)) +
                   a1 * (u(2, i + 1, j, q) - u(2, i - 1, j, q));
              ;
              ac1 = ac1 + bop(k, q) * c(31, i, j, q) * du;
              ac2 = ac2 + bop(k, q) * c(32, i, j, q) * du;
              ac3 = ac3 + bop(k, q) * c(33, i, j, q) * du;
            }
            r1 = r1 + strx(i) * ac1;
            r2 = r2 + strx(i) * ac2;
            r3 = r3 + strx(i) * ac3;

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              du = a2 * (u(3, i + 2, j, q) - u(3, i - 2, j, q)) +
                   a1 * (u(3, i + 1, j, q) - u(3, i - 1, j, q));
              ;
              ac1 = ac1 + bop(k, q) * c(34, i, j, q) * du;
              ac2 = ac2 + bop(k, q) * c(35, i, j, q) * du;
              ac3 = ac3 + bop(k, q) * c(36, i, j, q) * du;
            }
            r1 = r1 + strx(i) * ac1;
            r2 = r2 + strx(i) * ac2;
            r3 = r3 + strx(i) * ac3;

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              dup2 = dup2 + bop(k, q) * u(1, i, j + 2, q);
              dup1 = dup1 + bop(k, q) * u(1, i, j + 1, q);
              dum1 = dum1 + bop(k, q) * u(1, i, j - 1, q);
              dum2 = dum2 + bop(k, q) * u(1, i, j - 2, q);
            }
            r1 = r1 + stry(j) * (a1 * (c(37, i, j + 1, k) * dup1 -
                                       c(37, i, j - 1, k) * dum1) +
                                 a2 * (c(37, i, j + 2, k) * dup2 -
                                       c(37, i, j - 2, k) * dum2));
            r2 = r2 + stry(j) * (a1 * (c(40, i, j + 1, k) * dup1 -
                                       c(40, i, j - 1, k) * dum1) +
                                 a2 * (c(40, i, j + 2, k) * dup2 -
                                       c(40, i, j - 2, k) * dum2));
            r3 = r3 + stry(j) * (a1 * (c(43, i, j + 1, k) * dup1 -
                                       c(43, i, j - 1, k) * dum1) +
                                 a2 * (c(43, i, j + 2, k) * dup2 -
                                       c(43, i, j - 2, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              dup2 = dup2 + bop(k, q) * u(2, i, j + 2, q);
              dup1 = dup1 + bop(k, q) * u(2, i, j + 1, q);
              dum1 = dum1 + bop(k, q) * u(2, i, j - 1, q);
              dum2 = dum2 + bop(k, q) * u(2, i, j - 2, q);
            }
            r1 = r1 + stry(j) * (a1 * (c(38, i, j + 1, k) * dup1 -
                                       c(38, i, j - 1, k) * dum1) +
                                 a2 * (c(38, i, j + 2, k) * dup2 -
                                       c(38, i, j - 2, k) * dum2));
            r2 = r2 + stry(j) * (a1 * (c(41, i, j + 1, k) * dup1 -
                                       c(41, i, j - 1, k) * dum1) +
                                 a2 * (c(41, i, j + 2, k) * dup2 -
                                       c(41, i, j - 2, k) * dum2));
            r3 = r3 + stry(j) * (a1 * (c(44, i, j + 1, k) * dup1 -
                                       c(44, i, j - 1, k) * dum1) +
                                 a2 * (c(44, i, j + 2, k) * dup2 -
                                       c(44, i, j - 2, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              dup2 = dup2 + bop(k, q) * u(3, i, j + 2, q);
              dup1 = dup1 + bop(k, q) * u(3, i, j + 1, q);
              dum1 = dum1 + bop(k, q) * u(3, i, j - 1, q);
              dum2 = dum2 + bop(k, q) * u(3, i, j - 2, q);
            }
            r1 = r1 + stry(j) * (a1 * (c(39, i, j + 1, k) * dup1 -
                                       c(39, i, j - 1, k) * dum1) +
                                 a2 * (c(39, i, j + 2, k) * dup2 -
                                       c(39, i, j - 2, k) * dum2));
            r2 = r2 + stry(j) * (a1 * (c(42, i, j + 1, k) * dup1 -
                                       c(42, i, j - 1, k) * dum1) +
                                 a2 * (c(42, i, j + 2, k) * dup2 -
                                       c(42, i, j - 2, k) * dum2));
            r3 = r3 + stry(j) * (a1 * (c(45, i, j + 1, k) * dup1 -
                                       c(45, i, j - 1, k) * dum1) +
                                 a2 * (c(45, i, j + 2, k) * dup2 -
                                       c(45, i, j - 2, k) * dum2));

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              du = a2 * (u(1, i, j + 2, q) - u(1, i, j - 2, q)) +
                   a1 * (u(1, i, j + 1, q) - u(1, i, j - 1, q));
              ;
              ac1 = ac1 + bop(k, q) * c(37, i, j, q) * du;
              ac2 = ac2 + bop(k, q) * c(38, i, j, q) * du;
              ac3 = ac3 + bop(k, q) * c(39, i, j, q) * du;
            }
            r1 = r1 + stry(j) * ac1;
            r2 = r2 + stry(j) * ac2;
            r3 = r3 + stry(j) * ac3;

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              du = a2 * (u(2, i, j + 2, q) - u(2, i, j - 2, q)) +
                   a1 * (u(2, i, j + 1, q) - u(2, i, j - 1, q));
              ;
              ac1 = ac1 + bop(k, q) * c(40, i, j, q) * du;
              ac2 = ac2 + bop(k, q) * c(41, i, j, q) * du;
              ac3 = ac3 + bop(k, q) * c(42, i, j, q) * du;
            }
            r1 = r1 + stry(j) * ac1;
            r2 = r2 + stry(j) * ac2;
            r3 = r3 + stry(j) * ac3;

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = 1; q <= 8; q++) {
              du = a2 * (u(3, i, j + 2, q) - u(3, i, j - 2, q)) +
                   a1 * (u(3, i, j + 1, q) - u(3, i, j - 1, q));
              ;
              ac1 = ac1 + bop(k, q) * c(43, i, j, q) * du;
              ac2 = ac2 + bop(k, q) * c(44, i, j, q) * du;
              ac3 = ac3 + bop(k, q) * c(45, i, j, q) * du;
            }
            r1 = r1 + stry(j) * ac1;
            r2 = r2 + stry(j) * ac2;
            r3 = r3 + stry(j) * ac3;
            lu(1, i, j, k) = r1 * ijac;
            lu(2, i, j, k) = r2 * ijac;
            lu(3, i, j, k) = r3 * ijac;
          }
    }
    if (onesided[5] == 1) {
      ke = nk - 6;
#pragma omp for
      for (sw4_type k = nk - 5; k <= nk; k++)
        for (sw4_type j = jfirst + 2; j <= jlast - 2; j++)
#pragma simd
#pragma ivdep
          for (sw4_type i = ifirst + 2; i <= ilast - 2; i++) {
            float_sw4 r1 = 0, r2 = 0, r3 = 0;
            float_sw4 ac1, ac2, ac3, ac4, ac5, ac6;
            float_sw4 dum2, dum1, du, dup1, dup2;
            float_sw4 ijac = 1 / jac(i, j, k);
            float_sw4 cm2 = c(1, i - 1, j, k) * strx(i - 1) -
                            0.75 * (c(1, i, j, k) * strx(i) +
                                    c(1, i - 2, j, k) * strx(i - 2));
            float_sw4 cm1 =
                c(1, i - 2, j, k) * strx(i - 2) +
                c(1, i + 1, j, k) * strx(i + 1) +
                3 * (c(1, i, j, k) * strx(i) + c(1, i - 1, j, k) * strx(i - 1));
            float_sw4 cp1 =
                c(1, i - 1, j, k) * strx(i - 1) +
                c(1, i + 2, j, k) * strx(i + 2) +
                3 * (c(1, i, j, k) * strx(i) + c(1, i + 1, j, k) * strx(i + 1));
            float_sw4 cp2 = c(1, i + 1, j, k) * strx(i + 1) -
                            0.75 * (c(1, i, j, k) * strx(i) +
                                    c(1, i + 2, j, k) * strx(i + 2));

            r1 = r1 + i6 * strx(i) *
                          (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
            cm2 = c(2, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(2, i, j, k) * strx(i) +
                          c(2, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(2, i - 2, j, k) * strx(i - 2) +
                c(2, i + 1, j, k) * strx(i + 1) +
                3 * (c(2, i, j, k) * strx(i) + c(2, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(2, i - 1, j, k) * strx(i - 1) +
                c(2, i + 2, j, k) * strx(i + 2) +
                3 * (c(2, i, j, k) * strx(i) + c(2, i + 1, j, k) * strx(i + 1));
            cp2 = c(2, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(2, i, j, k) * strx(i) +
                          c(2, i + 2, j, k) * strx(i + 2));

            r1 = r1 + i6 * strx(i) *
                          (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
            r2 = r2 + i6 * strx(i) *
                          (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
            cm2 = c(3, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(3, i, j, k) * strx(i) +
                          c(3, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(3, i - 2, j, k) * strx(i - 2) +
                c(3, i + 1, j, k) * strx(i + 1) +
                3 * (c(3, i, j, k) * strx(i) + c(3, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(3, i - 1, j, k) * strx(i - 1) +
                c(3, i + 2, j, k) * strx(i + 2) +
                3 * (c(3, i, j, k) * strx(i) + c(3, i + 1, j, k) * strx(i + 1));
            cp2 = c(3, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(3, i, j, k) * strx(i) +
                          c(3, i + 2, j, k) * strx(i + 2));

            r1 = r1 + i6 * strx(i) *
                          (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
            r3 = r3 + i6 * strx(i) *
                          (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
            cm2 = c(4, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(4, i, j, k) * strx(i) +
                          c(4, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(4, i - 2, j, k) * strx(i - 2) +
                c(4, i + 1, j, k) * strx(i + 1) +
                3 * (c(4, i, j, k) * strx(i) + c(4, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(4, i - 1, j, k) * strx(i - 1) +
                c(4, i + 2, j, k) * strx(i + 2) +
                3 * (c(4, i, j, k) * strx(i) + c(4, i + 1, j, k) * strx(i + 1));
            cp2 = c(4, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(4, i, j, k) * strx(i) +
                          c(4, i + 2, j, k) * strx(i + 2));

            r2 = r2 + i6 * strx(i) *
                          (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
            cm2 = c(5, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(5, i, j, k) * strx(i) +
                          c(5, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(5, i - 2, j, k) * strx(i - 2) +
                c(5, i + 1, j, k) * strx(i + 1) +
                3 * (c(5, i, j, k) * strx(i) + c(5, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(5, i - 1, j, k) * strx(i - 1) +
                c(5, i + 2, j, k) * strx(i + 2) +
                3 * (c(5, i, j, k) * strx(i) + c(5, i + 1, j, k) * strx(i + 1));
            cp2 = c(5, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(5, i, j, k) * strx(i) +
                          c(5, i + 2, j, k) * strx(i + 2));

            r2 = r2 + i6 * strx(i) *
                          (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
            r3 = r3 + i6 * strx(i) *
                          (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
            cm2 = c(6, i - 1, j, k) * strx(i - 1) -
                  0.75 * (c(6, i, j, k) * strx(i) +
                          c(6, i - 2, j, k) * strx(i - 2));
            cm1 =
                c(6, i - 2, j, k) * strx(i - 2) +
                c(6, i + 1, j, k) * strx(i + 1) +
                3 * (c(6, i, j, k) * strx(i) + c(6, i - 1, j, k) * strx(i - 1));
            cp1 =
                c(6, i - 1, j, k) * strx(i - 1) +
                c(6, i + 2, j, k) * strx(i + 2) +
                3 * (c(6, i, j, k) * strx(i) + c(6, i + 1, j, k) * strx(i + 1));
            cp2 = c(6, i + 1, j, k) * strx(i + 1) -
                  0.75 * (c(6, i, j, k) * strx(i) +
                          c(6, i + 2, j, k) * strx(i + 2));

            r3 = r3 + i6 * strx(i) *
                          (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
            cm2 = c(7, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(7, i, j, k) * stry(j) +
                          c(7, i, j - 2, k) * stry(j - 2));
            cm1 =
                c(7, i, j - 2, k) * stry(j - 2) +
                c(7, i, j + 1, k) * stry(j + 1) +
                3 * (c(7, i, j, k) * stry(j) + c(7, i, j - 1, k) * stry(j - 1));
            cp1 =
                c(7, i, j - 1, k) * stry(j - 1) +
                c(7, i, j + 2, k) * stry(j + 2) +
                3 * (c(7, i, j, k) * stry(j) + c(7, i, j + 1, k) * stry(j + 1));
            cp2 = c(7, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(7, i, j, k) * stry(j) +
                          c(7, i, j + 2, k) * stry(j + 2));

            r1 = r1 + i6 * stry(j) *
                          (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
            cm2 = c(8, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(8, i, j, k) * stry(j) +
                          c(8, i, j - 2, k) * stry(j - 2));
            cm1 =
                c(8, i, j - 2, k) * stry(j - 2) +
                c(8, i, j + 1, k) * stry(j + 1) +
                3 * (c(8, i, j, k) * stry(j) + c(8, i, j - 1, k) * stry(j - 1));
            cp1 =
                c(8, i, j - 1, k) * stry(j - 1) +
                c(8, i, j + 2, k) * stry(j + 2) +
                3 * (c(8, i, j, k) * stry(j) + c(8, i, j + 1, k) * stry(j + 1));
            cp2 = c(8, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(8, i, j, k) * stry(j) +
                          c(8, i, j + 2, k) * stry(j + 2));

            r1 = r1 + i6 * stry(j) *
                          (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
            r2 = r2 + i6 * stry(j) *
                          (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
            cm2 = c(9, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(9, i, j, k) * stry(j) +
                          c(9, i, j - 2, k) * stry(j - 2));
            cm1 =
                c(9, i, j - 2, k) * stry(j - 2) +
                c(9, i, j + 1, k) * stry(j + 1) +
                3 * (c(9, i, j, k) * stry(j) + c(9, i, j - 1, k) * stry(j - 1));
            cp1 =
                c(9, i, j - 1, k) * stry(j - 1) +
                c(9, i, j + 2, k) * stry(j + 2) +
                3 * (c(9, i, j, k) * stry(j) + c(9, i, j + 1, k) * stry(j + 1));
            cp2 = c(9, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(9, i, j, k) * stry(j) +
                          c(9, i, j + 2, k) * stry(j + 2));

            r1 = r1 + i6 * stry(j) *
                          (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
            r3 = r3 + i6 * stry(j) *
                          (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                           cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                           cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                           cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
            cm2 = c(10, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(10, i, j, k) * stry(j) +
                          c(10, i, j - 2, k) * stry(j - 2));
            cm1 = c(10, i, j - 2, k) * stry(j - 2) +
                  c(10, i, j + 1, k) * stry(j + 1) +
                  3 * (c(10, i, j, k) * stry(j) +
                       c(10, i, j - 1, k) * stry(j - 1));
            cp1 = c(10, i, j - 1, k) * stry(j - 1) +
                  c(10, i, j + 2, k) * stry(j + 2) +
                  3 * (c(10, i, j, k) * stry(j) +
                       c(10, i, j + 1, k) * stry(j + 1));
            cp2 = c(10, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(10, i, j, k) * stry(j) +
                          c(10, i, j + 2, k) * stry(j + 2));

            r2 = r2 + i6 * stry(j) *
                          (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
            cm2 = c(11, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(11, i, j, k) * stry(j) +
                          c(11, i, j - 2, k) * stry(j - 2));
            cm1 = c(11, i, j - 2, k) * stry(j - 2) +
                  c(11, i, j + 1, k) * stry(j + 1) +
                  3 * (c(11, i, j, k) * stry(j) +
                       c(11, i, j - 1, k) * stry(j - 1));
            cp1 = c(11, i, j - 1, k) * stry(j - 1) +
                  c(11, i, j + 2, k) * stry(j + 2) +
                  3 * (c(11, i, j, k) * stry(j) +
                       c(11, i, j + 1, k) * stry(j + 1));
            cp2 = c(11, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(11, i, j, k) * stry(j) +
                          c(11, i, j + 2, k) * stry(j + 2));

            r2 = r2 + i6 * stry(j) *
                          (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
            r3 = r3 + i6 * stry(j) *
                          (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                           cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                           cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                           cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
            cm2 = c(12, i, j - 1, k) * stry(j - 1) -
                  0.75 * (c(12, i, j, k) * stry(j) +
                          c(12, i, j - 2, k) * stry(j - 2));
            cm1 = c(12, i, j - 2, k) * stry(j - 2) +
                  c(12, i, j + 1, k) * stry(j + 1) +
                  3 * (c(12, i, j, k) * stry(j) +
                       c(12, i, j - 1, k) * stry(j - 1));
            cp1 = c(12, i, j - 1, k) * stry(j - 1) +
                  c(12, i, j + 2, k) * stry(j + 2) +
                  3 * (c(12, i, j, k) * stry(j) +
                       c(12, i, j + 1, k) * stry(j + 1));
            cp2 = c(12, i, j + 1, k) * stry(j + 1) -
                  0.75 * (c(12, i, j, k) * stry(j) +
                          c(12, i, j + 2, k) * stry(j + 2));

            r3 = r3 + i6 * stry(j) *
                          (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                           cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                           cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                           cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
            r1 = r1 + ghcof(nk - k + 1) * c(13, i, j, nk) * u(1, i, j, nk + 1) +
                 ghcof(nk - k + 1) * c(14, i, j, nk) * u(2, i, j, nk + 1) +
                 ghcof(nk - k + 1) * c(15, i, j, nk) * u(3, i, j, nk + 1);
            r2 = r2 + ghcof(nk - k + 1) * c(14, i, j, nk) * u(1, i, j, nk + 1) +
                 ghcof(nk - k + 1) * c(16, i, j, nk) * u(2, i, j, nk + 1) +
                 ghcof(nk - k + 1) * c(17, i, j, nk) * u(3, i, j, nk + 1);
            r3 = r3 + ghcof(nk - k + 1) * c(15, i, j, nk) * u(1, i, j, nk + 1) +
                 ghcof(nk - k + 1) * c(17, i, j, nk) * u(2, i, j, nk + 1) +
                 ghcof(nk - k + 1) * c(18, i, j, nk) * u(3, i, j, nk + 1);
            for (sw4_type q = nk - 7; q <= nk; q++) {
              ac1 = 0;
              ac2 = 0;
              ac3 = 0;
              ac4 = 0;
              ac5 = 0;
              ac6 = 0;
              for (sw4_type m = nk - 7; m <= nk; m++) {
                ac1 = ac1 +
                      acof(nk - k + 1, nk - q + 1, nk - m + 1) * c(13, i, j, m);
                ac2 = ac2 +
                      acof(nk - k + 1, nk - q + 1, nk - m + 1) * c(14, i, j, m);
                ac3 = ac3 +
                      acof(nk - k + 1, nk - q + 1, nk - m + 1) * c(15, i, j, m);
                ac4 = ac4 +
                      acof(nk - k + 1, nk - q + 1, nk - m + 1) * c(16, i, j, m);
                ac5 = ac5 +
                      acof(nk - k + 1, nk - q + 1, nk - m + 1) * c(17, i, j, m);
                ac6 = ac6 +
                      acof(nk - k + 1, nk - q + 1, nk - m + 1) * c(18, i, j, m);
              }
              r1 = r1 + ac1 * u(1, i, j, q) + ac2 * u(2, i, j, q) +
                   ac3 * u(3, i, j, q);
              r2 = r2 + ac2 * u(1, i, j, q) + ac4 * u(2, i, j, q) +
                   ac5 * u(3, i, j, q);
              r3 = r3 + ac3 * u(1, i, j, q) + ac5 * u(2, i, j, q) +
                   ac6 * u(3, i, j, q);
            }

            dup1 = a2 * (u(1, i + 1, j + 2, k) - u(1, i + 1, j - 2, k)) +
                   a1 * (u(1, i + 1, j + 1, k) - u(1, i + 1, j - 1, k));
            dum1 = a2 * (u(1, i - 1, j + 2, k) - u(1, i - 1, j - 2, k)) +
                   a1 * (u(1, i - 1, j + 1, k) - u(1, i - 1, j - 1, k));
            dup2 = a2 * (u(1, i + 2, j + 2, k) - u(1, i + 2, j - 2, k)) +
                   a1 * (u(1, i + 2, j + 1, k) - u(1, i + 2, j - 1, k));
            dum2 = a2 * (u(1, i - 2, j + 2, k) - u(1, i - 2, j - 2, k)) +
                   a1 * (u(1, i - 2, j + 1, k) - u(1, i - 2, j - 1, k));
            r1 = r1 + stry(j) * strx(i) *
                          (a1 * (c(19, i + 1, j, k) * dup1 -
                                 c(19, i - 1, j, k) * dum1) +
                           a2 * (c(19, i + 2, j, k) * dup2 -
                                 c(19, i - 2, j, k) * dum2));
            r2 = r2 + stry(j) * strx(i) *
                          (a1 * (c(22, i + 1, j, k) * dup1 -
                                 c(22, i - 1, j, k) * dum1) +
                           a2 * (c(22, i + 2, j, k) * dup2 -
                                 c(22, i - 2, j, k) * dum2));
            r3 = r3 + stry(j) * strx(i) *
                          (a1 * (c(25, i + 1, j, k) * dup1 -
                                 c(25, i - 1, j, k) * dum1) +
                           a2 * (c(25, i + 2, j, k) * dup2 -
                                 c(25, i - 2, j, k) * dum2));

            dup1 = a2 * (u(2, i + 1, j + 2, k) - u(2, i + 1, j - 2, k)) +
                   a1 * (u(2, i + 1, j + 1, k) - u(2, i + 1, j - 1, k));
            dum1 = a2 * (u(2, i - 1, j + 2, k) - u(2, i - 1, j - 2, k)) +
                   a1 * (u(2, i - 1, j + 1, k) - u(2, i - 1, j - 1, k));
            dup2 = a2 * (u(2, i + 2, j + 2, k) - u(2, i + 2, j - 2, k)) +
                   a1 * (u(2, i + 2, j + 1, k) - u(2, i + 2, j - 1, k));
            dum2 = a2 * (u(2, i - 2, j + 2, k) - u(2, i - 2, j - 2, k)) +
                   a1 * (u(2, i - 2, j + 1, k) - u(2, i - 2, j - 1, k));
            r1 = r1 + stry(j) * strx(i) *
                          (a1 * (c(20, i + 1, j, k) * dup1 -
                                 c(20, i - 1, j, k) * dum1) +
                           a2 * (c(20, i + 2, j, k) * dup2 -
                                 c(20, i - 2, j, k) * dum2));
            r2 = r2 + stry(j) * strx(i) *
                          (a1 * (c(23, i + 1, j, k) * dup1 -
                                 c(23, i - 1, j, k) * dum1) +
                           a2 * (c(23, i + 2, j, k) * dup2 -
                                 c(23, i - 2, j, k) * dum2));
            r3 = r3 + stry(j) * strx(i) *
                          (a1 * (c(26, i + 1, j, k) * dup1 -
                                 c(26, i - 1, j, k) * dum1) +
                           a2 * (c(26, i + 2, j, k) * dup2 -
                                 c(26, i - 2, j, k) * dum2));

            dup1 = a2 * (u(3, i + 1, j + 2, k) - u(3, i + 1, j - 2, k)) +
                   a1 * (u(3, i + 1, j + 1, k) - u(3, i + 1, j - 1, k));
            dum1 = a2 * (u(3, i - 1, j + 2, k) - u(3, i - 1, j - 2, k)) +
                   a1 * (u(3, i - 1, j + 1, k) - u(3, i - 1, j - 1, k));
            dup2 = a2 * (u(3, i + 2, j + 2, k) - u(3, i + 2, j - 2, k)) +
                   a1 * (u(3, i + 2, j + 1, k) - u(3, i + 2, j - 1, k));
            dum2 = a2 * (u(3, i - 2, j + 2, k) - u(3, i - 2, j - 2, k)) +
                   a1 * (u(3, i - 2, j + 1, k) - u(3, i - 2, j - 1, k));
            r1 = r1 + stry(j) * strx(i) *
                          (a1 * (c(21, i + 1, j, k) * dup1 -
                                 c(21, i - 1, j, k) * dum1) +
                           a2 * (c(21, i + 2, j, k) * dup2 -
                                 c(21, i - 2, j, k) * dum2));
            r2 = r2 + stry(j) * strx(i) *
                          (a1 * (c(24, i + 1, j, k) * dup1 -
                                 c(24, i - 1, j, k) * dum1) +
                           a2 * (c(24, i + 2, j, k) * dup2 -
                                 c(24, i - 2, j, k) * dum2));
            r3 = r3 + stry(j) * strx(i) *
                          (a1 * (c(27, i + 1, j, k) * dup1 -
                                 c(27, i - 1, j, k) * dum1) +
                           a2 * (c(27, i + 2, j, k) * dup2 -
                                 c(27, i - 2, j, k) * dum2));

            dup1 = a2 * (u(1, i + 2, j + 1, k) - u(1, i - 2, j + 1, k)) +
                   a1 * (u(1, i + 1, j + 1, k) - u(1, i - 1, j + 1, k));
            dum1 = a2 * (u(1, i + 2, j - 1, k) - u(1, i - 2, j - 1, k)) +
                   a1 * (u(1, i + 1, j - 1, k) - u(1, i - 1, j - 1, k));
            dup2 = a2 * (u(1, i + 2, j + 2, k) - u(1, i - 2, j + 2, k)) +
                   a1 * (u(1, i + 1, j + 2, k) - u(1, i - 1, j + 2, k));
            dum2 = a2 * (u(1, i + 2, j - 2, k) - u(1, i - 2, j - 2, k)) +
                   a1 * (u(1, i + 1, j - 2, k) - u(1, i - 1, j - 2, k));
            r1 = r1 + strx(i) * stry(j) *
                          (a1 * (c(19, i, j + 1, k) * dup1 -
                                 c(19, i, j - 1, k) * dum1) +
                           a2 * (c(19, i, j + 2, k) * dup2 -
                                 c(19, i, j - 2, k) * dum2));
            r2 = r2 + strx(i) * stry(j) *
                          (a1 * (c(20, i, j + 1, k) * dup1 -
                                 c(20, i, j - 1, k) * dum1) +
                           a2 * (c(20, i, j + 2, k) * dup2 -
                                 c(20, i, j - 2, k) * dum2));
            r3 = r3 + strx(i) * stry(j) *
                          (a1 * (c(21, i, j + 1, k) * dup1 -
                                 c(21, i, j - 1, k) * dum1) +
                           a2 * (c(21, i, j + 2, k) * dup2 -
                                 c(21, i, j - 2, k) * dum2));

            dup1 = a2 * (u(2, i + 2, j + 1, k) - u(2, i - 2, j + 1, k)) +
                   a1 * (u(2, i + 1, j + 1, k) - u(2, i - 1, j + 1, k));
            dum1 = a2 * (u(2, i + 2, j - 1, k) - u(2, i - 2, j - 1, k)) +
                   a1 * (u(2, i + 1, j - 1, k) - u(2, i - 1, j - 1, k));
            dup2 = a2 * (u(2, i + 2, j + 2, k) - u(2, i - 2, j + 2, k)) +
                   a1 * (u(2, i + 1, j + 2, k) - u(2, i - 1, j + 2, k));
            dum2 = a2 * (u(2, i + 2, j - 2, k) - u(2, i - 2, j - 2, k)) +
                   a1 * (u(2, i + 1, j - 2, k) - u(2, i - 1, j - 2, k));
            r1 = r1 + strx(i) * stry(j) *
                          (a1 * (c(22, i, j + 1, k) * dup1 -
                                 c(22, i, j - 1, k) * dum1) +
                           a2 * (c(22, i, j + 2, k) * dup2 -
                                 c(22, i, j - 2, k) * dum2));
            r2 = r2 + strx(i) * stry(j) *
                          (a1 * (c(23, i, j + 1, k) * dup1 -
                                 c(23, i, j - 1, k) * dum1) +
                           a2 * (c(23, i, j + 2, k) * dup2 -
                                 c(23, i, j - 2, k) * dum2));
            r3 = r3 + strx(i) * stry(j) *
                          (a1 * (c(24, i, j + 1, k) * dup1 -
                                 c(24, i, j - 1, k) * dum1) +
                           a2 * (c(24, i, j + 2, k) * dup2 -
                                 c(24, i, j - 2, k) * dum2));

            dup1 = a2 * (u(3, i + 2, j + 1, k) - u(3, i - 2, j + 1, k)) +
                   a1 * (u(3, i + 1, j + 1, k) - u(3, i - 1, j + 1, k));
            dum1 = a2 * (u(3, i + 2, j - 1, k) - u(3, i - 2, j - 1, k)) +
                   a1 * (u(3, i + 1, j - 1, k) - u(3, i - 1, j - 1, k));
            dup2 = a2 * (u(3, i + 2, j + 2, k) - u(3, i - 2, j + 2, k)) +
                   a1 * (u(3, i + 1, j + 2, k) - u(3, i - 1, j + 2, k));
            dum2 = a2 * (u(3, i + 2, j - 2, k) - u(3, i - 2, j - 2, k)) +
                   a1 * (u(3, i + 1, j - 2, k) - u(3, i - 1, j - 2, k));
            r1 = r1 + strx(i) * stry(j) *
                          (a1 * (c(25, i, j + 1, k) * dup1 -
                                 c(25, i, j - 1, k) * dum1) +
                           a2 * (c(25, i, j + 2, k) * dup2 -
                                 c(25, i, j - 2, k) * dum2));
            r2 = r2 + strx(i) * stry(j) *
                          (a1 * (c(26, i, j + 1, k) * dup1 -
                                 c(26, i, j - 1, k) * dum1) +
                           a2 * (c(26, i, j + 2, k) * dup2 -
                                 c(26, i, j - 2, k) * dum2));
            r3 = r3 + strx(i) * stry(j) *
                          (a1 * (c(27, i, j + 1, k) * dup1 -
                                 c(27, i, j - 1, k) * dum1) +
                           a2 * (c(27, i, j + 2, k) * dup2 -
                                 c(27, i, j - 2, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              dup2 = dup2 - bop(nk - k + 1, nk - q + 1) * u(1, i + 2, j, q);
              dup1 = dup1 - bop(nk - k + 1, nk - q + 1) * u(1, i + 1, j, q);
              dum1 = dum1 - bop(nk - k + 1, nk - q + 1) * u(1, i - 1, j, q);
              dum2 = dum2 - bop(nk - k + 1, nk - q + 1) * u(1, i - 2, j, q);
            }
            r1 = r1 + strx(i) * (a1 * (c(28, i + 1, j, k) * dup1 -
                                       c(28, i - 1, j, k) * dum1) +
                                 a2 * (c(28, i + 2, j, k) * dup2 -
                                       c(28, i - 2, j, k) * dum2));
            r2 = r2 + strx(i) * (a1 * (c(31, i + 1, j, k) * dup1 -
                                       c(31, i - 1, j, k) * dum1) +
                                 a2 * (c(31, i + 2, j, k) * dup2 -
                                       c(31, i - 2, j, k) * dum2));
            r3 = r3 + strx(i) * (a1 * (c(34, i + 1, j, k) * dup1 -
                                       c(34, i - 1, j, k) * dum1) +
                                 a2 * (c(34, i + 2, j, k) * dup2 -
                                       c(34, i - 2, j, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              dup2 = dup2 - bop(nk - k + 1, nk - q + 1) * u(2, i + 2, j, q);
              dup1 = dup1 - bop(nk - k + 1, nk - q + 1) * u(2, i + 1, j, q);
              dum1 = dum1 - bop(nk - k + 1, nk - q + 1) * u(2, i - 1, j, q);
              dum2 = dum2 - bop(nk - k + 1, nk - q + 1) * u(2, i - 2, j, q);
            }
            r1 = r1 + strx(i) * (a1 * (c(29, i + 1, j, k) * dup1 -
                                       c(29, i - 1, j, k) * dum1) +
                                 a2 * (c(29, i + 2, j, k) * dup2 -
                                       c(29, i - 2, j, k) * dum2));
            r2 = r2 + strx(i) * (a1 * (c(32, i + 1, j, k) * dup1 -
                                       c(32, i - 1, j, k) * dum1) +
                                 a2 * (c(32, i + 2, j, k) * dup2 -
                                       c(32, i - 2, j, k) * dum2));
            r3 = r3 + strx(i) * (a1 * (c(35, i + 1, j, k) * dup1 -
                                       c(35, i - 1, j, k) * dum1) +
                                 a2 * (c(35, i + 2, j, k) * dup2 -
                                       c(35, i - 2, j, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              dup2 = dup2 - bop(nk - k + 1, nk - q + 1) * u(3, i + 2, j, q);
              dup1 = dup1 - bop(nk - k + 1, nk - q + 1) * u(3, i + 1, j, q);
              dum1 = dum1 - bop(nk - k + 1, nk - q + 1) * u(3, i - 1, j, q);
              dum2 = dum2 - bop(nk - k + 1, nk - q + 1) * u(3, i - 2, j, q);
            }
            r1 = r1 + strx(i) * (a1 * (c(30, i + 1, j, k) * dup1 -
                                       c(30, i - 1, j, k) * dum1) +
                                 a2 * (c(30, i + 2, j, k) * dup2 -
                                       c(30, i - 2, j, k) * dum2));
            r2 = r2 + strx(i) * (a1 * (c(33, i + 1, j, k) * dup1 -
                                       c(33, i - 1, j, k) * dum1) +
                                 a2 * (c(33, i + 2, j, k) * dup2 -
                                       c(33, i - 2, j, k) * dum2));
            r3 = r3 + strx(i) * (a1 * (c(36, i + 1, j, k) * dup1 -
                                       c(36, i - 1, j, k) * dum1) +
                                 a2 * (c(36, i + 2, j, k) * dup2 -
                                       c(36, i - 2, j, k) * dum2));

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              du = a2 * (u(1, i + 2, j, q) - u(1, i - 2, j, q)) +
                   a1 * (u(1, i + 1, j, q) - u(1, i - 1, j, q));
              ;
              ac1 = ac1 - bop(nk - k + 1, nk - q + 1) * c(28, i, j, q) * du;
              ac2 = ac2 - bop(nk - k + 1, nk - q + 1) * c(29, i, j, q) * du;
              ac3 = ac3 - bop(nk - k + 1, nk - q + 1) * c(30, i, j, q) * du;
            }
            r1 = r1 + strx(i) * ac1;
            r2 = r2 + strx(i) * ac2;
            r3 = r3 + strx(i) * ac3;

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              du = a2 * (u(2, i + 2, j, q) - u(2, i - 2, j, q)) +
                   a1 * (u(2, i + 1, j, q) - u(2, i - 1, j, q));
              ;
              ac1 = ac1 - bop(nk - k + 1, nk - q + 1) * c(31, i, j, q) * du;
              ac2 = ac2 - bop(nk - k + 1, nk - q + 1) * c(32, i, j, q) * du;
              ac3 = ac3 - bop(nk - k + 1, nk - q + 1) * c(33, i, j, q) * du;
            }
            r1 = r1 + strx(i) * ac1;
            r2 = r2 + strx(i) * ac2;
            r3 = r3 + strx(i) * ac3;

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              du = a2 * (u(3, i + 2, j, q) - u(3, i - 2, j, q)) +
                   a1 * (u(3, i + 1, j, q) - u(3, i - 1, j, q));
              ;
              ac1 = ac1 - bop(nk - k + 1, nk - q + 1) * c(34, i, j, q) * du;
              ac2 = ac2 - bop(nk - k + 1, nk - q + 1) * c(35, i, j, q) * du;
              ac3 = ac3 - bop(nk - k + 1, nk - q + 1) * c(36, i, j, q) * du;
            }
            r1 = r1 + strx(i) * ac1;
            r2 = r2 + strx(i) * ac2;
            r3 = r3 + strx(i) * ac3;

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              dup2 = dup2 - bop(nk - k + 1, nk - q + 1) * u(1, i, j + 2, q);
              dup1 = dup1 - bop(nk - k + 1, nk - q + 1) * u(1, i, j + 1, q);
              dum1 = dum1 - bop(nk - k + 1, nk - q + 1) * u(1, i, j - 1, q);
              dum2 = dum2 - bop(nk - k + 1, nk - q + 1) * u(1, i, j - 2, q);
            }
            r1 = r1 + stry(j) * (a1 * (c(37, i, j + 1, k) * dup1 -
                                       c(37, i, j - 1, k) * dum1) +
                                 a2 * (c(37, i, j + 2, k) * dup2 -
                                       c(37, i, j - 2, k) * dum2));
            r2 = r2 + stry(j) * (a1 * (c(40, i, j + 1, k) * dup1 -
                                       c(40, i, j - 1, k) * dum1) +
                                 a2 * (c(40, i, j + 2, k) * dup2 -
                                       c(40, i, j - 2, k) * dum2));
            r3 = r3 + stry(j) * (a1 * (c(43, i, j + 1, k) * dup1 -
                                       c(43, i, j - 1, k) * dum1) +
                                 a2 * (c(43, i, j + 2, k) * dup2 -
                                       c(43, i, j - 2, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              dup2 = dup2 - bop(nk - k + 1, nk - q + 1) * u(2, i, j + 2, q);
              dup1 = dup1 - bop(nk - k + 1, nk - q + 1) * u(2, i, j + 1, q);
              dum1 = dum1 - bop(nk - k + 1, nk - q + 1) * u(2, i, j - 1, q);
              dum2 = dum2 - bop(nk - k + 1, nk - q + 1) * u(2, i, j - 2, q);
            }
            r1 = r1 + stry(j) * (a1 * (c(38, i, j + 1, k) * dup1 -
                                       c(38, i, j - 1, k) * dum1) +
                                 a2 * (c(38, i, j + 2, k) * dup2 -
                                       c(38, i, j - 2, k) * dum2));
            r2 = r2 + stry(j) * (a1 * (c(41, i, j + 1, k) * dup1 -
                                       c(41, i, j - 1, k) * dum1) +
                                 a2 * (c(41, i, j + 2, k) * dup2 -
                                       c(41, i, j - 2, k) * dum2));
            r3 = r3 + stry(j) * (a1 * (c(44, i, j + 1, k) * dup1 -
                                       c(44, i, j - 1, k) * dum1) +
                                 a2 * (c(44, i, j + 2, k) * dup2 -
                                       c(44, i, j - 2, k) * dum2));

            dum2 = 0;
            dum1 = 0;
            dup1 = 0;
            dup2 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              dup2 = dup2 - bop(nk - k + 1, nk - q + 1) * u(3, i, j + 2, q);
              dup1 = dup1 - bop(nk - k + 1, nk - q + 1) * u(3, i, j + 1, q);
              dum1 = dum1 - bop(nk - k + 1, nk - q + 1) * u(3, i, j - 1, q);
              dum2 = dum2 - bop(nk - k + 1, nk - q + 1) * u(3, i, j - 2, q);
            }
            r1 = r1 + stry(j) * (a1 * (c(39, i, j + 1, k) * dup1 -
                                       c(39, i, j - 1, k) * dum1) +
                                 a2 * (c(39, i, j + 2, k) * dup2 -
                                       c(39, i, j - 2, k) * dum2));
            r2 = r2 + stry(j) * (a1 * (c(42, i, j + 1, k) * dup1 -
                                       c(42, i, j - 1, k) * dum1) +
                                 a2 * (c(42, i, j + 2, k) * dup2 -
                                       c(42, i, j - 2, k) * dum2));
            r3 = r3 + stry(j) * (a1 * (c(45, i, j + 1, k) * dup1 -
                                       c(45, i, j - 1, k) * dum1) +
                                 a2 * (c(45, i, j + 2, k) * dup2 -
                                       c(45, i, j - 2, k) * dum2));

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              du = a2 * (u(1, i, j + 2, q) - u(1, i, j - 2, q)) +
                   a1 * (u(1, i, j + 1, q) - u(1, i, j - 1, q));
              ;
              ac1 = ac1 - bop(nk - k + 1, nk - q + 1) * c(37, i, j, q) * du;
              ac2 = ac2 - bop(nk - k + 1, nk - q + 1) * c(38, i, j, q) * du;
              ac3 = ac3 - bop(nk - k + 1, nk - q + 1) * c(39, i, j, q) * du;
            }
            r1 = r1 + stry(j) * ac1;
            r2 = r2 + stry(j) * ac2;
            r3 = r3 + stry(j) * ac3;

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              du = a2 * (u(2, i, j + 2, q) - u(2, i, j - 2, q)) +
                   a1 * (u(2, i, j + 1, q) - u(2, i, j - 1, q));
              ;
              ac1 = ac1 - bop(nk - k + 1, nk - q + 1) * c(40, i, j, q) * du;
              ac2 = ac2 - bop(nk - k + 1, nk - q + 1) * c(41, i, j, q) * du;
              ac3 = ac3 - bop(nk - k + 1, nk - q + 1) * c(42, i, j, q) * du;
            }
            r1 = r1 + stry(j) * ac1;
            r2 = r2 + stry(j) * ac2;
            r3 = r3 + stry(j) * ac3;

            ac1 = 0;
            ac2 = 0;
            ac3 = 0;
            for (sw4_type q = nk - 7; q <= nk; q++) {
              du = a2 * (u(3, i, j + 2, q) - u(3, i, j - 2, q)) +
                   a1 * (u(3, i, j + 1, q) - u(3, i, j - 1, q));
              ;
              ac1 = ac1 - bop(nk - k + 1, nk - q + 1) * c(43, i, j, q) * du;
              ac2 = ac2 - bop(nk - k + 1, nk - q + 1) * c(44, i, j, q) * du;
              ac3 = ac3 - bop(nk - k + 1, nk - q + 1) * c(45, i, j, q) * du;
            }
            r1 = r1 + stry(j) * ac1;
            r2 = r2 + stry(j) * ac2;
            r3 = r3 + stry(j) * ac3;
            lu(1, i, j, k) = r1 * ijac;
            lu(2, i, j, k) = r2 * ijac;
            lu(3, i, j, k) = r3 * ijac;
          }
    }
#pragma omp for
    for (sw4_type k = kb; k <= ke; k++)
      for (sw4_type j = jfirst + 2; j <= jlast - 2; j++)
#pragma simd
#pragma ivdep
        for (sw4_type i = ifirst + 2; i <= ilast - 2; i++) {
          float_sw4 r1 = 0, r2 = 0, r3 = 0;
          float_sw4 dum2, dum1, dup1, dup2;
          float_sw4 ijac = 1 / jac(i, j, k);
          float_sw4 cm2 = c(1, i - 1, j, k) * strx(i - 1) -
                          0.75 * (c(1, i, j, k) * strx(i) +
                                  c(1, i - 2, j, k) * strx(i - 2));
          float_sw4 cm1 =
              c(1, i - 2, j, k) * strx(i - 2) +
              c(1, i + 1, j, k) * strx(i + 1) +
              3 * (c(1, i, j, k) * strx(i) + c(1, i - 1, j, k) * strx(i - 1));
          float_sw4 cp1 =
              c(1, i - 1, j, k) * strx(i - 1) +
              c(1, i + 2, j, k) * strx(i + 2) +
              3 * (c(1, i, j, k) * strx(i) + c(1, i + 1, j, k) * strx(i + 1));
          float_sw4 cp2 = c(1, i + 1, j, k) * strx(i + 1) -
                          0.75 * (c(1, i, j, k) * strx(i) +
                                  c(1, i + 2, j, k) * strx(i + 2));

          r1 = r1 + i6 * strx(i) *
                        (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                         cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                         cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                         cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
          cm2 = c(2, i - 1, j, k) * strx(i - 1) -
                0.75 *
                    (c(2, i, j, k) * strx(i) + c(2, i - 2, j, k) * strx(i - 2));
          cm1 = c(2, i - 2, j, k) * strx(i - 2) +
                c(2, i + 1, j, k) * strx(i + 1) +
                3 * (c(2, i, j, k) * strx(i) + c(2, i - 1, j, k) * strx(i - 1));
          cp1 = c(2, i - 1, j, k) * strx(i - 1) +
                c(2, i + 2, j, k) * strx(i + 2) +
                3 * (c(2, i, j, k) * strx(i) + c(2, i + 1, j, k) * strx(i + 1));
          cp2 = c(2, i + 1, j, k) * strx(i + 1) -
                0.75 *
                    (c(2, i, j, k) * strx(i) + c(2, i + 2, j, k) * strx(i + 2));

          r1 = r1 + i6 * strx(i) *
                        (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                         cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                         cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                         cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
          r2 = r2 + i6 * strx(i) *
                        (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                         cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                         cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                         cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
          cm2 = c(3, i - 1, j, k) * strx(i - 1) -
                0.75 *
                    (c(3, i, j, k) * strx(i) + c(3, i - 2, j, k) * strx(i - 2));
          cm1 = c(3, i - 2, j, k) * strx(i - 2) +
                c(3, i + 1, j, k) * strx(i + 1) +
                3 * (c(3, i, j, k) * strx(i) + c(3, i - 1, j, k) * strx(i - 1));
          cp1 = c(3, i - 1, j, k) * strx(i - 1) +
                c(3, i + 2, j, k) * strx(i + 2) +
                3 * (c(3, i, j, k) * strx(i) + c(3, i + 1, j, k) * strx(i + 1));
          cp2 = c(3, i + 1, j, k) * strx(i + 1) -
                0.75 *
                    (c(3, i, j, k) * strx(i) + c(3, i + 2, j, k) * strx(i + 2));

          r1 = r1 + i6 * strx(i) *
                        (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                         cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                         cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                         cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
          r3 = r3 + i6 * strx(i) *
                        (cm2 * (u(1, i - 2, j, k) - u(1, i, j, k)) +
                         cm1 * (u(1, i - 1, j, k) - u(1, i, j, k)) +
                         cp1 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                         cp2 * (u(1, i + 2, j, k) - u(1, i, j, k)));
          cm2 = c(4, i - 1, j, k) * strx(i - 1) -
                0.75 *
                    (c(4, i, j, k) * strx(i) + c(4, i - 2, j, k) * strx(i - 2));
          cm1 = c(4, i - 2, j, k) * strx(i - 2) +
                c(4, i + 1, j, k) * strx(i + 1) +
                3 * (c(4, i, j, k) * strx(i) + c(4, i - 1, j, k) * strx(i - 1));
          cp1 = c(4, i - 1, j, k) * strx(i - 1) +
                c(4, i + 2, j, k) * strx(i + 2) +
                3 * (c(4, i, j, k) * strx(i) + c(4, i + 1, j, k) * strx(i + 1));
          cp2 = c(4, i + 1, j, k) * strx(i + 1) -
                0.75 *
                    (c(4, i, j, k) * strx(i) + c(4, i + 2, j, k) * strx(i + 2));

          r2 = r2 + i6 * strx(i) *
                        (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                         cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                         cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                         cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
          cm2 = c(5, i - 1, j, k) * strx(i - 1) -
                0.75 *
                    (c(5, i, j, k) * strx(i) + c(5, i - 2, j, k) * strx(i - 2));
          cm1 = c(5, i - 2, j, k) * strx(i - 2) +
                c(5, i + 1, j, k) * strx(i + 1) +
                3 * (c(5, i, j, k) * strx(i) + c(5, i - 1, j, k) * strx(i - 1));
          cp1 = c(5, i - 1, j, k) * strx(i - 1) +
                c(5, i + 2, j, k) * strx(i + 2) +
                3 * (c(5, i, j, k) * strx(i) + c(5, i + 1, j, k) * strx(i + 1));
          cp2 = c(5, i + 1, j, k) * strx(i + 1) -
                0.75 *
                    (c(5, i, j, k) * strx(i) + c(5, i + 2, j, k) * strx(i + 2));

          r2 = r2 + i6 * strx(i) *
                        (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                         cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                         cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                         cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
          r3 = r3 + i6 * strx(i) *
                        (cm2 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                         cm1 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                         cp1 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                         cp2 * (u(2, i + 2, j, k) - u(2, i, j, k)));
          cm2 = c(6, i - 1, j, k) * strx(i - 1) -
                0.75 *
                    (c(6, i, j, k) * strx(i) + c(6, i - 2, j, k) * strx(i - 2));
          cm1 = c(6, i - 2, j, k) * strx(i - 2) +
                c(6, i + 1, j, k) * strx(i + 1) +
                3 * (c(6, i, j, k) * strx(i) + c(6, i - 1, j, k) * strx(i - 1));
          cp1 = c(6, i - 1, j, k) * strx(i - 1) +
                c(6, i + 2, j, k) * strx(i + 2) +
                3 * (c(6, i, j, k) * strx(i) + c(6, i + 1, j, k) * strx(i + 1));
          cp2 = c(6, i + 1, j, k) * strx(i + 1) -
                0.75 *
                    (c(6, i, j, k) * strx(i) + c(6, i + 2, j, k) * strx(i + 2));

          r3 = r3 + i6 * strx(i) *
                        (cm2 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                         cm1 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                         cp1 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                         cp2 * (u(3, i + 2, j, k) - u(3, i, j, k)));
          cm2 = c(7, i, j - 1, k) * stry(j - 1) -
                0.75 *
                    (c(7, i, j, k) * stry(j) + c(7, i, j - 2, k) * stry(j - 2));
          cm1 = c(7, i, j - 2, k) * stry(j - 2) +
                c(7, i, j + 1, k) * stry(j + 1) +
                3 * (c(7, i, j, k) * stry(j) + c(7, i, j - 1, k) * stry(j - 1));
          cp1 = c(7, i, j - 1, k) * stry(j - 1) +
                c(7, i, j + 2, k) * stry(j + 2) +
                3 * (c(7, i, j, k) * stry(j) + c(7, i, j + 1, k) * stry(j + 1));
          cp2 = c(7, i, j + 1, k) * stry(j + 1) -
                0.75 *
                    (c(7, i, j, k) * stry(j) + c(7, i, j + 2, k) * stry(j + 2));

          r1 = r1 + i6 * stry(j) *
                        (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                         cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                         cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                         cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
          cm2 = c(8, i, j - 1, k) * stry(j - 1) -
                0.75 *
                    (c(8, i, j, k) * stry(j) + c(8, i, j - 2, k) * stry(j - 2));
          cm1 = c(8, i, j - 2, k) * stry(j - 2) +
                c(8, i, j + 1, k) * stry(j + 1) +
                3 * (c(8, i, j, k) * stry(j) + c(8, i, j - 1, k) * stry(j - 1));
          cp1 = c(8, i, j - 1, k) * stry(j - 1) +
                c(8, i, j + 2, k) * stry(j + 2) +
                3 * (c(8, i, j, k) * stry(j) + c(8, i, j + 1, k) * stry(j + 1));
          cp2 = c(8, i, j + 1, k) * stry(j + 1) -
                0.75 *
                    (c(8, i, j, k) * stry(j) + c(8, i, j + 2, k) * stry(j + 2));

          r1 = r1 + i6 * stry(j) *
                        (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                         cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                         cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                         cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
          r2 = r2 + i6 * stry(j) *
                        (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                         cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                         cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                         cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
          cm2 = c(9, i, j - 1, k) * stry(j - 1) -
                0.75 *
                    (c(9, i, j, k) * stry(j) + c(9, i, j - 2, k) * stry(j - 2));
          cm1 = c(9, i, j - 2, k) * stry(j - 2) +
                c(9, i, j + 1, k) * stry(j + 1) +
                3 * (c(9, i, j, k) * stry(j) + c(9, i, j - 1, k) * stry(j - 1));
          cp1 = c(9, i, j - 1, k) * stry(j - 1) +
                c(9, i, j + 2, k) * stry(j + 2) +
                3 * (c(9, i, j, k) * stry(j) + c(9, i, j + 1, k) * stry(j + 1));
          cp2 = c(9, i, j + 1, k) * stry(j + 1) -
                0.75 *
                    (c(9, i, j, k) * stry(j) + c(9, i, j + 2, k) * stry(j + 2));

          r1 = r1 + i6 * stry(j) *
                        (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                         cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                         cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                         cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
          r3 = r3 + i6 * stry(j) *
                        (cm2 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                         cm1 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                         cp1 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                         cp2 * (u(1, i, j + 2, k) - u(1, i, j, k)));
          cm2 = c(10, i, j - 1, k) * stry(j - 1) -
                0.75 * (c(10, i, j, k) * stry(j) +
                        c(10, i, j - 2, k) * stry(j - 2));
          cm1 =
              c(10, i, j - 2, k) * stry(j - 2) +
              c(10, i, j + 1, k) * stry(j + 1) +
              3 * (c(10, i, j, k) * stry(j) + c(10, i, j - 1, k) * stry(j - 1));
          cp1 =
              c(10, i, j - 1, k) * stry(j - 1) +
              c(10, i, j + 2, k) * stry(j + 2) +
              3 * (c(10, i, j, k) * stry(j) + c(10, i, j + 1, k) * stry(j + 1));
          cp2 = c(10, i, j + 1, k) * stry(j + 1) -
                0.75 * (c(10, i, j, k) * stry(j) +
                        c(10, i, j + 2, k) * stry(j + 2));

          r2 = r2 + i6 * stry(j) *
                        (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                         cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                         cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                         cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
          cm2 = c(11, i, j - 1, k) * stry(j - 1) -
                0.75 * (c(11, i, j, k) * stry(j) +
                        c(11, i, j - 2, k) * stry(j - 2));
          cm1 =
              c(11, i, j - 2, k) * stry(j - 2) +
              c(11, i, j + 1, k) * stry(j + 1) +
              3 * (c(11, i, j, k) * stry(j) + c(11, i, j - 1, k) * stry(j - 1));
          cp1 =
              c(11, i, j - 1, k) * stry(j - 1) +
              c(11, i, j + 2, k) * stry(j + 2) +
              3 * (c(11, i, j, k) * stry(j) + c(11, i, j + 1, k) * stry(j + 1));
          cp2 = c(11, i, j + 1, k) * stry(j + 1) -
                0.75 * (c(11, i, j, k) * stry(j) +
                        c(11, i, j + 2, k) * stry(j + 2));

          r2 = r2 + i6 * stry(j) *
                        (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                         cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                         cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                         cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
          r3 = r3 + i6 * stry(j) *
                        (cm2 * (u(2, i, j - 2, k) - u(2, i, j, k)) +
                         cm1 * (u(2, i, j - 1, k) - u(2, i, j, k)) +
                         cp1 * (u(2, i, j + 1, k) - u(2, i, j, k)) +
                         cp2 * (u(2, i, j + 2, k) - u(2, i, j, k)));
          cm2 = c(12, i, j - 1, k) * stry(j - 1) -
                0.75 * (c(12, i, j, k) * stry(j) +
                        c(12, i, j - 2, k) * stry(j - 2));
          cm1 =
              c(12, i, j - 2, k) * stry(j - 2) +
              c(12, i, j + 1, k) * stry(j + 1) +
              3 * (c(12, i, j, k) * stry(j) + c(12, i, j - 1, k) * stry(j - 1));
          cp1 =
              c(12, i, j - 1, k) * stry(j - 1) +
              c(12, i, j + 2, k) * stry(j + 2) +
              3 * (c(12, i, j, k) * stry(j) + c(12, i, j + 1, k) * stry(j + 1));
          cp2 = c(12, i, j + 1, k) * stry(j + 1) -
                0.75 * (c(12, i, j, k) * stry(j) +
                        c(12, i, j + 2, k) * stry(j + 2));

          r3 = r3 + i6 * stry(j) *
                        (cm2 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                         cm1 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                         cp1 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                         cp2 * (u(3, i, j + 2, k) - u(3, i, j, k)));
          cm2 = c(13, i, j, k - 1) * strz(k - 1) -
                0.75 * (c(13, i, j, k) * strz(k) +
                        c(13, i, j, k - 2) * strz(k - 2));
          cm1 =
              c(13, i, j, k - 2) * strz(k - 2) +
              c(13, i, j, k + 1) * strz(k + 1) +
              3 * (c(13, i, j, k) * strz(k) + c(13, i, j, k - 1) * strz(k - 1));
          cp1 =
              c(13, i, j, k - 1) * strz(k - 1) +
              c(13, i, j, k + 2) * strz(k + 2) +
              3 * (c(13, i, j, k) * strz(k) + c(13, i, j, k + 1) * strz(k + 1));
          cp2 = c(13, i, j, k + 1) * strz(k + 1) -
                0.75 * (c(13, i, j, k) * strz(k) +
                        c(13, i, j, k + 2) * strz(k + 2));

          r1 = r1 + i6 * strz(k) *
                        (cm2 * (u(1, i, j, k - 2) - u(1, i, j, k)) +
                         cm1 * (u(1, i, j, k - 1) - u(1, i, j, k)) +
                         cp1 * (u(1, i, j, k + 1) - u(1, i, j, k)) +
                         cp2 * (u(1, i, j, k + 2) - u(1, i, j, k)));
          cm2 = c(14, i, j, k - 1) * strz(k - 1) -
                0.75 * (c(14, i, j, k) * strz(k) +
                        c(14, i, j, k - 2) * strz(k - 2));
          cm1 =
              c(14, i, j, k - 2) * strz(k - 2) +
              c(14, i, j, k + 1) * strz(k + 1) +
              3 * (c(14, i, j, k) * strz(k) + c(14, i, j, k - 1) * strz(k - 1));
          cp1 =
              c(14, i, j, k - 1) * strz(k - 1) +
              c(14, i, j, k + 2) * strz(k + 2) +
              3 * (c(14, i, j, k) * strz(k) + c(14, i, j, k + 1) * strz(k + 1));
          cp2 = c(14, i, j, k + 1) * strz(k + 1) -
                0.75 * (c(14, i, j, k) * strz(k) +
                        c(14, i, j, k + 2) * strz(k + 2));

          r1 = r1 + i6 * strz(k) *
                        (cm2 * (u(2, i, j, k - 2) - u(2, i, j, k)) +
                         cm1 * (u(2, i, j, k - 1) - u(2, i, j, k)) +
                         cp1 * (u(2, i, j, k + 1) - u(2, i, j, k)) +
                         cp2 * (u(2, i, j, k + 2) - u(2, i, j, k)));
          r2 = r2 + i6 * strz(k) *
                        (cm2 * (u(1, i, j, k - 2) - u(1, i, j, k)) +
                         cm1 * (u(1, i, j, k - 1) - u(1, i, j, k)) +
                         cp1 * (u(1, i, j, k + 1) - u(1, i, j, k)) +
                         cp2 * (u(1, i, j, k + 2) - u(1, i, j, k)));
          cm2 = c(15, i, j, k - 1) * strz(k - 1) -
                0.75 * (c(15, i, j, k) * strz(k) +
                        c(15, i, j, k - 2) * strz(k - 2));
          cm1 =
              c(15, i, j, k - 2) * strz(k - 2) +
              c(15, i, j, k + 1) * strz(k + 1) +
              3 * (c(15, i, j, k) * strz(k) + c(15, i, j, k - 1) * strz(k - 1));
          cp1 =
              c(15, i, j, k - 1) * strz(k - 1) +
              c(15, i, j, k + 2) * strz(k + 2) +
              3 * (c(15, i, j, k) * strz(k) + c(15, i, j, k + 1) * strz(k + 1));
          cp2 = c(15, i, j, k + 1) * strz(k + 1) -
                0.75 * (c(15, i, j, k) * strz(k) +
                        c(15, i, j, k + 2) * strz(k + 2));

          r1 = r1 + i6 * strz(k) *
                        (cm2 * (u(3, i, j, k - 2) - u(3, i, j, k)) +
                         cm1 * (u(3, i, j, k - 1) - u(3, i, j, k)) +
                         cp1 * (u(3, i, j, k + 1) - u(3, i, j, k)) +
                         cp2 * (u(3, i, j, k + 2) - u(3, i, j, k)));
          r3 = r3 + i6 * strz(k) *
                        (cm2 * (u(1, i, j, k - 2) - u(1, i, j, k)) +
                         cm1 * (u(1, i, j, k - 1) - u(1, i, j, k)) +
                         cp1 * (u(1, i, j, k + 1) - u(1, i, j, k)) +
                         cp2 * (u(1, i, j, k + 2) - u(1, i, j, k)));
          cm2 = c(16, i, j, k - 1) * strz(k - 1) -
                0.75 * (c(16, i, j, k) * strz(k) +
                        c(16, i, j, k - 2) * strz(k - 2));
          cm1 =
              c(16, i, j, k - 2) * strz(k - 2) +
              c(16, i, j, k + 1) * strz(k + 1) +
              3 * (c(16, i, j, k) * strz(k) + c(16, i, j, k - 1) * strz(k - 1));
          cp1 =
              c(16, i, j, k - 1) * strz(k - 1) +
              c(16, i, j, k + 2) * strz(k + 2) +
              3 * (c(16, i, j, k) * strz(k) + c(16, i, j, k + 1) * strz(k + 1));
          cp2 = c(16, i, j, k + 1) * strz(k + 1) -
                0.75 * (c(16, i, j, k) * strz(k) +
                        c(16, i, j, k + 2) * strz(k + 2));

          r2 = r2 + i6 * strz(k) *
                        (cm2 * (u(2, i, j, k - 2) - u(2, i, j, k)) +
                         cm1 * (u(2, i, j, k - 1) - u(2, i, j, k)) +
                         cp1 * (u(2, i, j, k + 1) - u(2, i, j, k)) +
                         cp2 * (u(2, i, j, k + 2) - u(2, i, j, k)));
          cm2 = c(17, i, j, k - 1) * strz(k - 1) -
                0.75 * (c(17, i, j, k) * strz(k) +
                        c(17, i, j, k - 2) * strz(k - 2));
          cm1 =
              c(17, i, j, k - 2) * strz(k - 2) +
              c(17, i, j, k + 1) * strz(k + 1) +
              3 * (c(17, i, j, k) * strz(k) + c(17, i, j, k - 1) * strz(k - 1));
          cp1 =
              c(17, i, j, k - 1) * strz(k - 1) +
              c(17, i, j, k + 2) * strz(k + 2) +
              3 * (c(17, i, j, k) * strz(k) + c(17, i, j, k + 1) * strz(k + 1));
          cp2 = c(17, i, j, k + 1) * strz(k + 1) -
                0.75 * (c(17, i, j, k) * strz(k) +
                        c(17, i, j, k + 2) * strz(k + 2));

          r2 = r2 + i6 * strz(k) *
                        (cm2 * (u(3, i, j, k - 2) - u(3, i, j, k)) +
                         cm1 * (u(3, i, j, k - 1) - u(3, i, j, k)) +
                         cp1 * (u(3, i, j, k + 1) - u(3, i, j, k)) +
                         cp2 * (u(3, i, j, k + 2) - u(3, i, j, k)));
          r3 = r3 + i6 * strz(k) *
                        (cm2 * (u(2, i, j, k - 2) - u(2, i, j, k)) +
                         cm1 * (u(2, i, j, k - 1) - u(2, i, j, k)) +
                         cp1 * (u(2, i, j, k + 1) - u(2, i, j, k)) +
                         cp2 * (u(2, i, j, k + 2) - u(2, i, j, k)));
          cm2 = c(18, i, j, k - 1) * strz(k - 1) -
                0.75 * (c(18, i, j, k) * strz(k) +
                        c(18, i, j, k - 2) * strz(k - 2));
          cm1 =
              c(18, i, j, k - 2) * strz(k - 2) +
              c(18, i, j, k + 1) * strz(k + 1) +
              3 * (c(18, i, j, k) * strz(k) + c(18, i, j, k - 1) * strz(k - 1));
          cp1 =
              c(18, i, j, k - 1) * strz(k - 1) +
              c(18, i, j, k + 2) * strz(k + 2) +
              3 * (c(18, i, j, k) * strz(k) + c(18, i, j, k + 1) * strz(k + 1));
          cp2 = c(18, i, j, k + 1) * strz(k + 1) -
                0.75 * (c(18, i, j, k) * strz(k) +
                        c(18, i, j, k + 2) * strz(k + 2));

          r3 = r3 + i6 * strz(k) *
                        (cm2 * (u(3, i, j, k - 2) - u(3, i, j, k)) +
                         cm1 * (u(3, i, j, k - 1) - u(3, i, j, k)) +
                         cp1 * (u(3, i, j, k + 1) - u(3, i, j, k)) +
                         cp2 * (u(3, i, j, k + 2) - u(3, i, j, k)));

          dup1 = a2 * (u(1, i + 1, j + 2, k) - u(1, i + 1, j - 2, k)) +
                 a1 * (u(1, i + 1, j + 1, k) - u(1, i + 1, j - 1, k));
          dum1 = a2 * (u(1, i - 1, j + 2, k) - u(1, i - 1, j - 2, k)) +
                 a1 * (u(1, i - 1, j + 1, k) - u(1, i - 1, j - 1, k));
          dup2 = a2 * (u(1, i + 2, j + 2, k) - u(1, i + 2, j - 2, k)) +
                 a1 * (u(1, i + 2, j + 1, k) - u(1, i + 2, j - 1, k));
          dum2 = a2 * (u(1, i - 2, j + 2, k) - u(1, i - 2, j - 2, k)) +
                 a1 * (u(1, i - 2, j + 1, k) - u(1, i - 2, j - 1, k));
          r1 = r1 + stry(j) * strx(i) *
                        (a1 * (c(19, i + 1, j, k) * dup1 -
                               c(19, i - 1, j, k) * dum1) +
                         a2 * (c(19, i + 2, j, k) * dup2 -
                               c(19, i - 2, j, k) * dum2));
          r2 = r2 + stry(j) * strx(i) *
                        (a1 * (c(22, i + 1, j, k) * dup1 -
                               c(22, i - 1, j, k) * dum1) +
                         a2 * (c(22, i + 2, j, k) * dup2 -
                               c(22, i - 2, j, k) * dum2));
          r3 = r3 + stry(j) * strx(i) *
                        (a1 * (c(25, i + 1, j, k) * dup1 -
                               c(25, i - 1, j, k) * dum1) +
                         a2 * (c(25, i + 2, j, k) * dup2 -
                               c(25, i - 2, j, k) * dum2));

          dup1 = a2 * (u(2, i + 1, j + 2, k) - u(2, i + 1, j - 2, k)) +
                 a1 * (u(2, i + 1, j + 1, k) - u(2, i + 1, j - 1, k));
          dum1 = a2 * (u(2, i - 1, j + 2, k) - u(2, i - 1, j - 2, k)) +
                 a1 * (u(2, i - 1, j + 1, k) - u(2, i - 1, j - 1, k));
          dup2 = a2 * (u(2, i + 2, j + 2, k) - u(2, i + 2, j - 2, k)) +
                 a1 * (u(2, i + 2, j + 1, k) - u(2, i + 2, j - 1, k));
          dum2 = a2 * (u(2, i - 2, j + 2, k) - u(2, i - 2, j - 2, k)) +
                 a1 * (u(2, i - 2, j + 1, k) - u(2, i - 2, j - 1, k));
          r1 = r1 + stry(j) * strx(i) *
                        (a1 * (c(20, i + 1, j, k) * dup1 -
                               c(20, i - 1, j, k) * dum1) +
                         a2 * (c(20, i + 2, j, k) * dup2 -
                               c(20, i - 2, j, k) * dum2));
          r2 = r2 + stry(j) * strx(i) *
                        (a1 * (c(23, i + 1, j, k) * dup1 -
                               c(23, i - 1, j, k) * dum1) +
                         a2 * (c(23, i + 2, j, k) * dup2 -
                               c(23, i - 2, j, k) * dum2));
          r3 = r3 + stry(j) * strx(i) *
                        (a1 * (c(26, i + 1, j, k) * dup1 -
                               c(26, i - 1, j, k) * dum1) +
                         a2 * (c(26, i + 2, j, k) * dup2 -
                               c(26, i - 2, j, k) * dum2));

          dup1 = a2 * (u(3, i + 1, j + 2, k) - u(3, i + 1, j - 2, k)) +
                 a1 * (u(3, i + 1, j + 1, k) - u(3, i + 1, j - 1, k));
          dum1 = a2 * (u(3, i - 1, j + 2, k) - u(3, i - 1, j - 2, k)) +
                 a1 * (u(3, i - 1, j + 1, k) - u(3, i - 1, j - 1, k));
          dup2 = a2 * (u(3, i + 2, j + 2, k) - u(3, i + 2, j - 2, k)) +
                 a1 * (u(3, i + 2, j + 1, k) - u(3, i + 2, j - 1, k));
          dum2 = a2 * (u(3, i - 2, j + 2, k) - u(3, i - 2, j - 2, k)) +
                 a1 * (u(3, i - 2, j + 1, k) - u(3, i - 2, j - 1, k));
          r1 = r1 + stry(j) * strx(i) *
                        (a1 * (c(21, i + 1, j, k) * dup1 -
                               c(21, i - 1, j, k) * dum1) +
                         a2 * (c(21, i + 2, j, k) * dup2 -
                               c(21, i - 2, j, k) * dum2));
          r2 = r2 + stry(j) * strx(i) *
                        (a1 * (c(24, i + 1, j, k) * dup1 -
                               c(24, i - 1, j, k) * dum1) +
                         a2 * (c(24, i + 2, j, k) * dup2 -
                               c(24, i - 2, j, k) * dum2));
          r3 = r3 + stry(j) * strx(i) *
                        (a1 * (c(27, i + 1, j, k) * dup1 -
                               c(27, i - 1, j, k) * dum1) +
                         a2 * (c(27, i + 2, j, k) * dup2 -
                               c(27, i - 2, j, k) * dum2));

          dup1 = a2 * (u(1, i + 2, j + 1, k) - u(1, i - 2, j + 1, k)) +
                 a1 * (u(1, i + 1, j + 1, k) - u(1, i - 1, j + 1, k));
          dum1 = a2 * (u(1, i + 2, j - 1, k) - u(1, i - 2, j - 1, k)) +
                 a1 * (u(1, i + 1, j - 1, k) - u(1, i - 1, j - 1, k));
          dup2 = a2 * (u(1, i + 2, j + 2, k) - u(1, i - 2, j + 2, k)) +
                 a1 * (u(1, i + 1, j + 2, k) - u(1, i - 1, j + 2, k));
          dum2 = a2 * (u(1, i + 2, j - 2, k) - u(1, i - 2, j - 2, k)) +
                 a1 * (u(1, i + 1, j - 2, k) - u(1, i - 1, j - 2, k));
          r1 = r1 + strx(i) * stry(j) *
                        (a1 * (c(19, i, j + 1, k) * dup1 -
                               c(19, i, j - 1, k) * dum1) +
                         a2 * (c(19, i, j + 2, k) * dup2 -
                               c(19, i, j - 2, k) * dum2));
          r2 = r2 + strx(i) * stry(j) *
                        (a1 * (c(20, i, j + 1, k) * dup1 -
                               c(20, i, j - 1, k) * dum1) +
                         a2 * (c(20, i, j + 2, k) * dup2 -
                               c(20, i, j - 2, k) * dum2));
          r3 = r3 + strx(i) * stry(j) *
                        (a1 * (c(21, i, j + 1, k) * dup1 -
                               c(21, i, j - 1, k) * dum1) +
                         a2 * (c(21, i, j + 2, k) * dup2 -
                               c(21, i, j - 2, k) * dum2));

          dup1 = a2 * (u(2, i + 2, j + 1, k) - u(2, i - 2, j + 1, k)) +
                 a1 * (u(2, i + 1, j + 1, k) - u(2, i - 1, j + 1, k));
          dum1 = a2 * (u(2, i + 2, j - 1, k) - u(2, i - 2, j - 1, k)) +
                 a1 * (u(2, i + 1, j - 1, k) - u(2, i - 1, j - 1, k));
          dup2 = a2 * (u(2, i + 2, j + 2, k) - u(2, i - 2, j + 2, k)) +
                 a1 * (u(2, i + 1, j + 2, k) - u(2, i - 1, j + 2, k));
          dum2 = a2 * (u(2, i + 2, j - 2, k) - u(2, i - 2, j - 2, k)) +
                 a1 * (u(2, i + 1, j - 2, k) - u(2, i - 1, j - 2, k));
          r1 = r1 + strx(i) * stry(j) *
                        (a1 * (c(22, i, j + 1, k) * dup1 -
                               c(22, i, j - 1, k) * dum1) +
                         a2 * (c(22, i, j + 2, k) * dup2 -
                               c(22, i, j - 2, k) * dum2));
          r2 = r2 + strx(i) * stry(j) *
                        (a1 * (c(23, i, j + 1, k) * dup1 -
                               c(23, i, j - 1, k) * dum1) +
                         a2 * (c(23, i, j + 2, k) * dup2 -
                               c(23, i, j - 2, k) * dum2));
          r3 = r3 + strx(i) * stry(j) *
                        (a1 * (c(24, i, j + 1, k) * dup1 -
                               c(24, i, j - 1, k) * dum1) +
                         a2 * (c(24, i, j + 2, k) * dup2 -
                               c(24, i, j - 2, k) * dum2));

          dup1 = a2 * (u(3, i + 2, j + 1, k) - u(3, i - 2, j + 1, k)) +
                 a1 * (u(3, i + 1, j + 1, k) - u(3, i - 1, j + 1, k));
          dum1 = a2 * (u(3, i + 2, j - 1, k) - u(3, i - 2, j - 1, k)) +
                 a1 * (u(3, i + 1, j - 1, k) - u(3, i - 1, j - 1, k));
          dup2 = a2 * (u(3, i + 2, j + 2, k) - u(3, i - 2, j + 2, k)) +
                 a1 * (u(3, i + 1, j + 2, k) - u(3, i - 1, j + 2, k));
          dum2 = a2 * (u(3, i + 2, j - 2, k) - u(3, i - 2, j - 2, k)) +
                 a1 * (u(3, i + 1, j - 2, k) - u(3, i - 1, j - 2, k));
          r1 = r1 + strx(i) * stry(j) *
                        (a1 * (c(25, i, j + 1, k) * dup1 -
                               c(25, i, j - 1, k) * dum1) +
                         a2 * (c(25, i, j + 2, k) * dup2 -
                               c(25, i, j - 2, k) * dum2));
          r2 = r2 + strx(i) * stry(j) *
                        (a1 * (c(26, i, j + 1, k) * dup1 -
                               c(26, i, j - 1, k) * dum1) +
                         a2 * (c(26, i, j + 2, k) * dup2 -
                               c(26, i, j - 2, k) * dum2));
          r3 = r3 + strx(i) * stry(j) *
                        (a1 * (c(27, i, j + 1, k) * dup1 -
                               c(27, i, j - 1, k) * dum1) +
                         a2 * (c(27, i, j + 2, k) * dup2 -
                               c(27, i, j - 2, k) * dum2));

          dup1 = a2 * (u(1, i + 1, j, k + 2) - u(1, i + 1, j, k - 2)) +
                 a1 * (u(1, i + 1, j, k + 1) - u(1, i + 1, j, k - 1));
          dum1 = a2 * (u(1, i - 1, j, k + 2) - u(1, i - 1, j, k - 2)) +
                 a1 * (u(1, i - 1, j, k + 1) - u(1, i - 1, j, k - 1));
          dup2 = a2 * (u(1, i + 2, j, k + 2) - u(1, i + 2, j, k - 2)) +
                 a1 * (u(1, i + 2, j, k + 1) - u(1, i + 2, j, k - 1));
          dum2 = a2 * (u(1, i - 2, j, k + 2) - u(1, i - 2, j, k - 2)) +
                 a1 * (u(1, i - 2, j, k + 1) - u(1, i - 2, j, k - 1));
          r1 = r1 + strz(k) * strx(i) *
                        (a1 * (c(28, i + 1, j, k) * dup1 -
                               c(28, i - 1, j, k) * dum1) +
                         a2 * (c(28, i + 2, j, k) * dup2 -
                               c(28, i - 2, j, k) * dum2));
          r2 = r2 + strz(k) * strx(i) *
                        (a1 * (c(31, i + 1, j, k) * dup1 -
                               c(31, i - 1, j, k) * dum1) +
                         a2 * (c(31, i + 2, j, k) * dup2 -
                               c(31, i - 2, j, k) * dum2));
          r3 = r3 + strz(k) * strx(i) *
                        (a1 * (c(34, i + 1, j, k) * dup1 -
                               c(34, i - 1, j, k) * dum1) +
                         a2 * (c(34, i + 2, j, k) * dup2 -
                               c(34, i - 2, j, k) * dum2));

          dup1 = a2 * (u(2, i + 1, j, k + 2) - u(2, i + 1, j, k - 2)) +
                 a1 * (u(2, i + 1, j, k + 1) - u(2, i + 1, j, k - 1));
          dum1 = a2 * (u(2, i - 1, j, k + 2) - u(2, i - 1, j, k - 2)) +
                 a1 * (u(2, i - 1, j, k + 1) - u(2, i - 1, j, k - 1));
          dup2 = a2 * (u(2, i + 2, j, k + 2) - u(2, i + 2, j, k - 2)) +
                 a1 * (u(2, i + 2, j, k + 1) - u(2, i + 2, j, k - 1));
          dum2 = a2 * (u(2, i - 2, j, k + 2) - u(2, i - 2, j, k - 2)) +
                 a1 * (u(2, i - 2, j, k + 1) - u(2, i - 2, j, k - 1));
          r1 = r1 + strz(k) * strx(i) *
                        (a1 * (c(29, i + 1, j, k) * dup1 -
                               c(29, i - 1, j, k) * dum1) +
                         a2 * (c(29, i + 2, j, k) * dup2 -
                               c(29, i - 2, j, k) * dum2));
          r2 = r2 + strz(k) * strx(i) *
                        (a1 * (c(32, i + 1, j, k) * dup1 -
                               c(32, i - 1, j, k) * dum1) +
                         a2 * (c(32, i + 2, j, k) * dup2 -
                               c(32, i - 2, j, k) * dum2));
          r3 = r3 + strz(k) * strx(i) *
                        (a1 * (c(35, i + 1, j, k) * dup1 -
                               c(35, i - 1, j, k) * dum1) +
                         a2 * (c(35, i + 2, j, k) * dup2 -
                               c(35, i - 2, j, k) * dum2));

          dup1 = a2 * (u(3, i + 1, j, k + 2) - u(3, i + 1, j, k - 2)) +
                 a1 * (u(3, i + 1, j, k + 1) - u(3, i + 1, j, k - 1));
          dum1 = a2 * (u(3, i - 1, j, k + 2) - u(3, i - 1, j, k - 2)) +
                 a1 * (u(3, i - 1, j, k + 1) - u(3, i - 1, j, k - 1));
          dup2 = a2 * (u(3, i + 2, j, k + 2) - u(3, i + 2, j, k - 2)) +
                 a1 * (u(3, i + 2, j, k + 1) - u(3, i + 2, j, k - 1));
          dum2 = a2 * (u(3, i - 2, j, k + 2) - u(3, i - 2, j, k - 2)) +
                 a1 * (u(3, i - 2, j, k + 1) - u(3, i - 2, j, k - 1));
          r1 = r1 + strz(k) * strx(i) *
                        (a1 * (c(30, i + 1, j, k) * dup1 -
                               c(30, i - 1, j, k) * dum1) +
                         a2 * (c(30, i + 2, j, k) * dup2 -
                               c(30, i - 2, j, k) * dum2));
          r2 = r2 + strz(k) * strx(i) *
                        (a1 * (c(33, i + 1, j, k) * dup1 -
                               c(33, i - 1, j, k) * dum1) +
                         a2 * (c(33, i + 2, j, k) * dup2 -
                               c(33, i - 2, j, k) * dum2));
          r3 = r3 + strz(k) * strx(i) *
                        (a1 * (c(36, i + 1, j, k) * dup1 -
                               c(36, i - 1, j, k) * dum1) +
                         a2 * (c(36, i + 2, j, k) * dup2 -
                               c(36, i - 2, j, k) * dum2));

          dup1 = a2 * (u(1, i + 2, j, k + 1) - u(1, i - 2, j, k + 1)) +
                 a1 * (u(1, i + 1, j, k + 1) - u(1, i - 1, j, k + 1));
          dum1 = a2 * (u(1, i + 2, j, k - 1) - u(1, i - 2, j, k - 1)) +
                 a1 * (u(1, i + 1, j, k - 1) - u(1, i - 1, j, k - 1));
          dup2 = a2 * (u(1, i + 2, j, k + 2) - u(1, i - 2, j, k + 2)) +
                 a1 * (u(1, i + 1, j, k + 2) - u(1, i - 1, j, k + 2));
          dum2 = a2 * (u(1, i + 2, j, k - 2) - u(1, i - 2, j, k - 2)) +
                 a1 * (u(1, i + 1, j, k - 2) - u(1, i - 1, j, k - 2));
          r1 = r1 + strx(i) * strz(k) *
                        (a1 * (c(28, i, j, k + 1) * dup1 -
                               c(28, i, j, k - 1) * dum1) +
                         a2 * (c(28, i, j, k + 2) * dup2 -
                               c(28, i, j, k - 2) * dum2));
          r2 = r2 + strx(i) * strz(k) *
                        (a1 * (c(29, i, j, k + 1) * dup1 -
                               c(29, i, j, k - 1) * dum1) +
                         a2 * (c(29, i, j, k + 2) * dup2 -
                               c(29, i, j, k - 2) * dum2));
          r3 = r3 + strx(i) * strz(k) *
                        (a1 * (c(30, i, j, k + 1) * dup1 -
                               c(30, i, j, k - 1) * dum1) +
                         a2 * (c(30, i, j, k + 2) * dup2 -
                               c(30, i, j, k - 2) * dum2));

          dup1 = a2 * (u(2, i + 2, j, k + 1) - u(2, i - 2, j, k + 1)) +
                 a1 * (u(2, i + 1, j, k + 1) - u(2, i - 1, j, k + 1));
          dum1 = a2 * (u(2, i + 2, j, k - 1) - u(2, i - 2, j, k - 1)) +
                 a1 * (u(2, i + 1, j, k - 1) - u(2, i - 1, j, k - 1));
          dup2 = a2 * (u(2, i + 2, j, k + 2) - u(2, i - 2, j, k + 2)) +
                 a1 * (u(2, i + 1, j, k + 2) - u(2, i - 1, j, k + 2));
          dum2 = a2 * (u(2, i + 2, j, k - 2) - u(2, i - 2, j, k - 2)) +
                 a1 * (u(2, i + 1, j, k - 2) - u(2, i - 1, j, k - 2));
          r1 = r1 + strx(i) * strz(k) *
                        (a1 * (c(31, i, j, k + 1) * dup1 -
                               c(31, i, j, k - 1) * dum1) +
                         a2 * (c(31, i, j, k + 2) * dup2 -
                               c(31, i, j, k - 2) * dum2));
          r2 = r2 + strx(i) * strz(k) *
                        (a1 * (c(32, i, j, k + 1) * dup1 -
                               c(32, i, j, k - 1) * dum1) +
                         a2 * (c(32, i, j, k + 2) * dup2 -
                               c(32, i, j, k - 2) * dum2));
          r3 = r3 + strx(i) * strz(k) *
                        (a1 * (c(33, i, j, k + 1) * dup1 -
                               c(33, i, j, k - 1) * dum1) +
                         a2 * (c(33, i, j, k + 2) * dup2 -
                               c(33, i, j, k - 2) * dum2));

          dup1 = a2 * (u(3, i + 2, j, k + 1) - u(3, i - 2, j, k + 1)) +
                 a1 * (u(3, i + 1, j, k + 1) - u(3, i - 1, j, k + 1));
          dum1 = a2 * (u(3, i + 2, j, k - 1) - u(3, i - 2, j, k - 1)) +
                 a1 * (u(3, i + 1, j, k - 1) - u(3, i - 1, j, k - 1));
          dup2 = a2 * (u(3, i + 2, j, k + 2) - u(3, i - 2, j, k + 2)) +
                 a1 * (u(3, i + 1, j, k + 2) - u(3, i - 1, j, k + 2));
          dum2 = a2 * (u(3, i + 2, j, k - 2) - u(3, i - 2, j, k - 2)) +
                 a1 * (u(3, i + 1, j, k - 2) - u(3, i - 1, j, k - 2));
          r1 = r1 + strx(i) * strz(k) *
                        (a1 * (c(34, i, j, k + 1) * dup1 -
                               c(34, i, j, k - 1) * dum1) +
                         a2 * (c(34, i, j, k + 2) * dup2 -
                               c(34, i, j, k - 2) * dum2));
          r2 = r2 + strx(i) * strz(k) *
                        (a1 * (c(35, i, j, k + 1) * dup1 -
                               c(35, i, j, k - 1) * dum1) +
                         a2 * (c(35, i, j, k + 2) * dup2 -
                               c(35, i, j, k - 2) * dum2));
          r3 = r3 + strx(i) * strz(k) *
                        (a1 * (c(36, i, j, k + 1) * dup1 -
                               c(36, i, j, k - 1) * dum1) +
                         a2 * (c(36, i, j, k + 2) * dup2 -
                               c(36, i, j, k - 2) * dum2));

          dup1 = a2 * (u(1, i, j + 1, k + 2) - u(1, i, j + 1, k - 2)) +
                 a1 * (u(1, i, j + 1, k + 1) - u(1, i, j + 1, k - 1));
          dum1 = a2 * (u(1, i, j - 1, k + 2) - u(1, i, j - 1, k - 2)) +
                 a1 * (u(1, i, j - 1, k + 1) - u(1, i, j - 1, k - 1));
          dup2 = a2 * (u(1, i, j + 2, k + 2) - u(1, i, j + 2, k - 2)) +
                 a1 * (u(1, i, j + 2, k + 1) - u(1, i, j + 2, k - 1));
          dum2 = a2 * (u(1, i, j - 2, k + 2) - u(1, i, j - 2, k - 2)) +
                 a1 * (u(1, i, j - 2, k + 1) - u(1, i, j - 2, k - 1));
          r1 = r1 + strz(k) * stry(j) *
                        (a1 * (c(37, i, j + 1, k) * dup1 -
                               c(37, i, j - 1, k) * dum1) +
                         a2 * (c(37, i, j + 2, k) * dup2 -
                               c(37, i, j - 2, k) * dum2));
          r2 = r2 + strz(k) * stry(j) *
                        (a1 * (c(40, i, j + 1, k) * dup1 -
                               c(40, i, j - 1, k) * dum1) +
                         a2 * (c(40, i, j + 2, k) * dup2 -
                               c(40, i, j - 2, k) * dum2));
          r3 = r3 + strz(k) * stry(j) *
                        (a1 * (c(43, i, j + 1, k) * dup1 -
                               c(43, i, j - 1, k) * dum1) +
                         a2 * (c(43, i, j + 2, k) * dup2 -
                               c(43, i, j - 2, k) * dum2));

          dup1 = a2 * (u(2, i, j + 1, k + 2) - u(2, i, j + 1, k - 2)) +
                 a1 * (u(2, i, j + 1, k + 1) - u(2, i, j + 1, k - 1));
          dum1 = a2 * (u(2, i, j - 1, k + 2) - u(2, i, j - 1, k - 2)) +
                 a1 * (u(2, i, j - 1, k + 1) - u(2, i, j - 1, k - 1));
          dup2 = a2 * (u(2, i, j + 2, k + 2) - u(2, i, j + 2, k - 2)) +
                 a1 * (u(2, i, j + 2, k + 1) - u(2, i, j + 2, k - 1));
          dum2 = a2 * (u(2, i, j - 2, k + 2) - u(2, i, j - 2, k - 2)) +
                 a1 * (u(2, i, j - 2, k + 1) - u(2, i, j - 2, k - 1));
          r1 = r1 + strz(k) * stry(j) *
                        (a1 * (c(38, i, j + 1, k) * dup1 -
                               c(38, i, j - 1, k) * dum1) +
                         a2 * (c(38, i, j + 2, k) * dup2 -
                               c(38, i, j - 2, k) * dum2));
          r2 = r2 + strz(k) * stry(j) *
                        (a1 * (c(41, i, j + 1, k) * dup1 -
                               c(41, i, j - 1, k) * dum1) +
                         a2 * (c(41, i, j + 2, k) * dup2 -
                               c(41, i, j - 2, k) * dum2));
          r3 = r3 + strz(k) * stry(j) *
                        (a1 * (c(44, i, j + 1, k) * dup1 -
                               c(44, i, j - 1, k) * dum1) +
                         a2 * (c(44, i, j + 2, k) * dup2 -
                               c(44, i, j - 2, k) * dum2));

          dup1 = a2 * (u(3, i, j + 1, k + 2) - u(3, i, j + 1, k - 2)) +
                 a1 * (u(3, i, j + 1, k + 1) - u(3, i, j + 1, k - 1));
          dum1 = a2 * (u(3, i, j - 1, k + 2) - u(3, i, j - 1, k - 2)) +
                 a1 * (u(3, i, j - 1, k + 1) - u(3, i, j - 1, k - 1));
          dup2 = a2 * (u(3, i, j + 2, k + 2) - u(3, i, j + 2, k - 2)) +
                 a1 * (u(3, i, j + 2, k + 1) - u(3, i, j + 2, k - 1));
          dum2 = a2 * (u(3, i, j - 2, k + 2) - u(3, i, j - 2, k - 2)) +
                 a1 * (u(3, i, j - 2, k + 1) - u(3, i, j - 2, k - 1));
          r1 = r1 + strz(k) * stry(j) *
                        (a1 * (c(39, i, j + 1, k) * dup1 -
                               c(39, i, j - 1, k) * dum1) +
                         a2 * (c(39, i, j + 2, k) * dup2 -
                               c(39, i, j - 2, k) * dum2));
          r2 = r2 + strz(k) * stry(j) *
                        (a1 * (c(42, i, j + 1, k) * dup1 -
                               c(42, i, j - 1, k) * dum1) +
                         a2 * (c(42, i, j + 2, k) * dup2 -
                               c(42, i, j - 2, k) * dum2));
          r3 = r3 + strz(k) * stry(j) *
                        (a1 * (c(45, i, j + 1, k) * dup1 -
                               c(45, i, j - 1, k) * dum1) +
                         a2 * (c(45, i, j + 2, k) * dup2 -
                               c(45, i, j - 2, k) * dum2));

          dup1 = a2 * (u(1, i, j + 2, k + 1) - u(1, i, j - 2, k + 1)) +
                 a1 * (u(1, i, j + 1, k + 1) - u(1, i, j - 1, k + 1));
          dum1 = a2 * (u(1, i, j + 2, k - 1) - u(1, i, j - 2, k - 1)) +
                 a1 * (u(1, i, j + 1, k - 1) - u(1, i, j - 1, k - 1));
          dup2 = a2 * (u(1, i, j + 2, k + 2) - u(1, i, j - 2, k + 2)) +
                 a1 * (u(1, i, j + 1, k + 2) - u(1, i, j - 1, k + 2));
          dum2 = a2 * (u(1, i, j + 2, k - 2) - u(1, i, j - 2, k - 2)) +
                 a1 * (u(1, i, j + 1, k - 2) - u(1, i, j - 1, k - 2));
          r1 = r1 + stry(j) * strz(k) *
                        (a1 * (c(37, i, j, k + 1) * dup1 -
                               c(37, i, j, k - 1) * dum1) +
                         a2 * (c(37, i, j, k + 2) * dup2 -
                               c(37, i, j, k - 2) * dum2));
          r2 = r2 + stry(j) * strz(k) *
                        (a1 * (c(38, i, j, k + 1) * dup1 -
                               c(38, i, j, k - 1) * dum1) +
                         a2 * (c(38, i, j, k + 2) * dup2 -
                               c(38, i, j, k - 2) * dum2));
          r3 = r3 + stry(j) * strz(k) *
                        (a1 * (c(39, i, j, k + 1) * dup1 -
                               c(39, i, j, k - 1) * dum1) +
                         a2 * (c(39, i, j, k + 2) * dup2 -
                               c(39, i, j, k - 2) * dum2));

          dup1 = a2 * (u(2, i, j + 2, k + 1) - u(2, i, j - 2, k + 1)) +
                 a1 * (u(2, i, j + 1, k + 1) - u(2, i, j - 1, k + 1));
          dum1 = a2 * (u(2, i, j + 2, k - 1) - u(2, i, j - 2, k - 1)) +
                 a1 * (u(2, i, j + 1, k - 1) - u(2, i, j - 1, k - 1));
          dup2 = a2 * (u(2, i, j + 2, k + 2) - u(2, i, j - 2, k + 2)) +
                 a1 * (u(2, i, j + 1, k + 2) - u(2, i, j - 1, k + 2));
          dum2 = a2 * (u(2, i, j + 2, k - 2) - u(2, i, j - 2, k - 2)) +
                 a1 * (u(2, i, j + 1, k - 2) - u(2, i, j - 1, k - 2));
          r1 = r1 + stry(j) * strz(k) *
                        (a1 * (c(40, i, j, k + 1) * dup1 -
                               c(40, i, j, k - 1) * dum1) +
                         a2 * (c(40, i, j, k + 2) * dup2 -
                               c(40, i, j, k - 2) * dum2));
          r2 = r2 + stry(j) * strz(k) *
                        (a1 * (c(41, i, j, k + 1) * dup1 -
                               c(41, i, j, k - 1) * dum1) +
                         a2 * (c(41, i, j, k + 2) * dup2 -
                               c(41, i, j, k - 2) * dum2));
          r3 = r3 + stry(j) * strz(k) *
                        (a1 * (c(42, i, j, k + 1) * dup1 -
                               c(42, i, j, k - 1) * dum1) +
                         a2 * (c(42, i, j, k + 2) * dup2 -
                               c(42, i, j, k - 2) * dum2));

          dup1 = a2 * (u(3, i, j + 2, k + 1) - u(3, i, j - 2, k + 1)) +
                 a1 * (u(3, i, j + 1, k + 1) - u(3, i, j - 1, k + 1));
          dum1 = a2 * (u(3, i, j + 2, k - 1) - u(3, i, j - 2, k - 1)) +
                 a1 * (u(3, i, j + 1, k - 1) - u(3, i, j - 1, k - 1));
          dup2 = a2 * (u(3, i, j + 2, k + 2) - u(3, i, j - 2, k + 2)) +
                 a1 * (u(3, i, j + 1, k + 2) - u(3, i, j - 1, k + 2));
          dum2 = a2 * (u(3, i, j + 2, k - 2) - u(3, i, j - 2, k - 2)) +
                 a1 * (u(3, i, j + 1, k - 2) - u(3, i, j - 1, k - 2));
          r1 = r1 + stry(j) * strz(k) *
                        (a1 * (c(43, i, j, k + 1) * dup1 -
                               c(43, i, j, k - 1) * dum1) +
                         a2 * (c(43, i, j, k + 2) * dup2 -
                               c(43, i, j, k - 2) * dum2));
          r2 = r2 + stry(j) * strz(k) *
                        (a1 * (c(44, i, j, k + 1) * dup1 -
                               c(44, i, j, k - 1) * dum1) +
                         a2 * (c(44, i, j, k + 2) * dup2 -
                               c(44, i, j, k - 2) * dum2));
          r3 = r3 + stry(j) * strz(k) *
                        (a1 * (c(45, i, j, k + 1) * dup1 -
                               c(45, i, j, k - 1) * dum1) +
                         a2 * (c(45, i, j, k + 2) * dup2 -
                               c(45, i, j, k - 2) * dum2));

          lu(1, i, j, k) = r1 * ijac;
          lu(2, i, j, k) = r2 * ijac;
          lu(3, i, j, k) = r3 * ijac;
        }
  }
}
