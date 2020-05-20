#include <cmath>
#include <cstdlib>

#include "Mspace.h"
#include "caliper.h"
#include "foralls.h"
#include "forallsdecl.h"
#include "policies.h"
#include "sw4.h"
//#ifdef ENABLE_GPU
// j#include <cuda_runtime.h>
//#include "optimizedcuda.h"
//#endif
//#include "tests.h"
// extern "C" {

void rhs4th3fort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                    int klast, int nk, int* __restrict__ onesided,
                    float_sw4* __restrict__ a_acof,
                    float_sw4* __restrict__ a_bope,
                    float_sw4* __restrict__ a_ghcof,
                    float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
                    float_sw4* __restrict__ a_mu,
                    float_sw4* __restrict__ a_lambda, float_sw4 h, char op) {
  SW4_MARK_FUNCTION;
  // Direct reuse of fortran code by these macro definitions:
#define mu(i, j, k) a_mu[base + i + ni * (j) + nij * (k)]
#define la(i, j, k) a_lambda[base + i + ni * (j) + nij * (k)]
  // Reversed indexation
#define u(c, i, j, k) a_u[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define lu(c, i, j, k) a_lu[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define strx(i) 1
#define stry(j) 1
#define strz(k) 1
#define acof(i, j, k) a_acof[(i - 1) + 6 * (j - 1) + 48 * (k - 1)]
#define bope(i, j) a_bope[i - 1 + 6 * (j - 1)]
#define ghcof(i) a_ghcof[i - 1]

  const float_sw4 i6 = 1.0 / 6;
  const float_sw4 i12 = 1.0 / 12;
  const float_sw4 i144 = 1.0 / 144;
  const float_sw4 tf = 0.75;

  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;

  int k1, k2;

  float_sw4 cof = 1.0 / (h * h);
  int a1;

  if (op == '=')
    a1 = 0;
  else if (op == '+')
    a1 = 1;
  else if (op == '-') {
    a1 = 1;
    cof = -cof;
  }

  k1 = kfirst + 2;
  if (onesided[4] == 1) k1 = 7;
  k2 = klast - 2;
  if (onesided[5] == 1) k2 = nk - 6;

  ASSERT_MANAGED(a_mu);
  ASSERT_MANAGED(a_lambda);
  ASSERT_MANAGED(a_u);
  ASSERT_MANAGED(a_lu);
  ASSERT_MANAGED(a_acof);
  ASSERT_MANAGED(a_bope);
  ASSERT_MANAGED(a_ghcof);
  {
    RAJA::RangeSegment k_range(k1, k2 + 1);
    RAJA::RangeSegment j_range(jfirst + 2, jlast - 1);
    RAJA::RangeSegment i_range(ifirst + 2, ilast - 1);
    RAJA::kernel<RHS4_EXEC_POL>(
        RAJA::make_tuple(k_range, j_range, i_range),
        [=] RAJA_DEVICE(int k, int j, int i) {
          float_sw4 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4, muz1, muz2,
              muz3, muz4;
          float_sw4 r1, r2, r3;
          /* from inner_loop_4a, 28x3 = 84 ops */
          mux1 = mu(i - 1, j, k) * strx(i - 1) -
                 tf * (mu(i, j, k) * strx(i) + mu(i - 2, j, k) * strx(i - 2));
          mux2 = mu(i - 2, j, k) * strx(i - 2) + mu(i + 1, j, k) * strx(i + 1) +
                 3 * (mu(i, j, k) * strx(i) + mu(i - 1, j, k) * strx(i - 1));
          mux3 = mu(i - 1, j, k) * strx(i - 1) + mu(i + 2, j, k) * strx(i + 2) +
                 3 * (mu(i + 1, j, k) * strx(i + 1) + mu(i, j, k) * strx(i));
          mux4 = mu(i + 1, j, k) * strx(i + 1) -
                 tf * (mu(i, j, k) * strx(i) + mu(i + 2, j, k) * strx(i + 2));

          muy1 = mu(i, j - 1, k) * stry(j - 1) -
                 tf * (mu(i, j, k) * stry(j) + mu(i, j - 2, k) * stry(j - 2));
          muy2 = mu(i, j - 2, k) * stry(j - 2) + mu(i, j + 1, k) * stry(j + 1) +
                 3 * (mu(i, j, k) * stry(j) + mu(i, j - 1, k) * stry(j - 1));
          muy3 = mu(i, j - 1, k) * stry(j - 1) + mu(i, j + 2, k) * stry(j + 2) +
                 3 * (mu(i, j + 1, k) * stry(j + 1) + mu(i, j, k) * stry(j));
          muy4 = mu(i, j + 1, k) * stry(j + 1) -
                 tf * (mu(i, j, k) * stry(j) + mu(i, j + 2, k) * stry(j + 2));

          muz1 = mu(i, j, k - 1) * strz(k - 1) -
                 tf * (mu(i, j, k) * strz(k) + mu(i, j, k - 2) * strz(k - 2));
          muz2 = mu(i, j, k - 2) * strz(k - 2) + mu(i, j, k + 1) * strz(k + 1) +
                 3 * (mu(i, j, k) * strz(k) + mu(i, j, k - 1) * strz(k - 1));
          muz3 = mu(i, j, k - 1) * strz(k - 1) + mu(i, j, k + 2) * strz(k + 2) +
                 3 * (mu(i, j, k + 1) * strz(k + 1) + mu(i, j, k) * strz(k));
          muz4 = mu(i, j, k + 1) * strz(k + 1) -
                 tf * (mu(i, j, k) * strz(k) + mu(i, j, k + 2) * strz(k + 2));
          /* xx, yy, and zz derivatives:*/
          /* 75 ops */
          r1 = i6 * (strx(i) * ((2 * mux1 + la(i - 1, j, k) * strx(i - 1) -
                                 tf * (la(i, j, k) * strx(i) +
                                       la(i - 2, j, k) * strx(i - 2))) *
                                    (u(1, i - 2, j, k) - u(1, i, j, k)) +
                                (2 * mux2 + la(i - 2, j, k) * strx(i - 2) +
                                 la(i + 1, j, k) * strx(i + 1) +
                                 3 * (la(i, j, k) * strx(i) +
                                      la(i - 1, j, k) * strx(i - 1))) *
                                    (u(1, i - 1, j, k) - u(1, i, j, k)) +
                                (2 * mux3 + la(i - 1, j, k) * strx(i - 1) +
                                 la(i + 2, j, k) * strx(i + 2) +
                                 3 * (la(i + 1, j, k) * strx(i + 1) +
                                      la(i, j, k) * strx(i))) *
                                    (u(1, i + 1, j, k) - u(1, i, j, k)) +
                                (2 * mux4 + la(i + 1, j, k) * strx(i + 1) -
                                 tf * (la(i, j, k) * strx(i) +
                                       la(i + 2, j, k) * strx(i + 2))) *
                                    (u(1, i + 2, j, k) - u(1, i, j, k))) +
                     stry(j) * (muy1 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                                muy2 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                                muy3 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                                muy4 * (u(1, i, j + 2, k) - u(1, i, j, k))) +
                     strz(k) * (muz1 * (u(1, i, j, k - 2) - u(1, i, j, k)) +
                                muz2 * (u(1, i, j, k - 1) - u(1, i, j, k)) +
                                muz3 * (u(1, i, j, k + 1) - u(1, i, j, k)) +
                                muz4 * (u(1, i, j, k + 2) - u(1, i, j, k))));

          /* 75 ops */
          r2 = i6 * (strx(i) * (mux1 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                                mux2 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                                mux3 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                                mux4 * (u(2, i + 2, j, k) - u(2, i, j, k))) +
                     stry(j) * ((2 * muy1 + la(i, j - 1, k) * stry(j - 1) -
                                 tf * (la(i, j, k) * stry(j) +
                                       la(i, j - 2, k) * stry(j - 2))) *
                                    (u(2, i, j - 2, k) - u(2, i, j, k)) +
                                (2 * muy2 + la(i, j - 2, k) * stry(j - 2) +
                                 la(i, j + 1, k) * stry(j + 1) +
                                 3 * (la(i, j, k) * stry(j) +
                                      la(i, j - 1, k) * stry(j - 1))) *
                                    (u(2, i, j - 1, k) - u(2, i, j, k)) +
                                (2 * muy3 + la(i, j - 1, k) * stry(j - 1) +
                                 la(i, j + 2, k) * stry(j + 2) +
                                 3 * (la(i, j + 1, k) * stry(j + 1) +
                                      la(i, j, k) * stry(j))) *
                                    (u(2, i, j + 1, k) - u(2, i, j, k)) +
                                (2 * muy4 + la(i, j + 1, k) * stry(j + 1) -
                                 tf * (la(i, j, k) * stry(j) +
                                       la(i, j + 2, k) * stry(j + 2))) *
                                    (u(2, i, j + 2, k) - u(2, i, j, k))) +
                     strz(k) * (muz1 * (u(2, i, j, k - 2) - u(2, i, j, k)) +
                                muz2 * (u(2, i, j, k - 1) - u(2, i, j, k)) +
                                muz3 * (u(2, i, j, k + 1) - u(2, i, j, k)) +
                                muz4 * (u(2, i, j, k + 2) - u(2, i, j, k))));

          /* 75 ops */
          r3 = i6 * (strx(i) * (mux1 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                                mux2 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                                mux3 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                                mux4 * (u(3, i + 2, j, k) - u(3, i, j, k))) +
                     stry(j) * (muy1 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                                muy2 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                                muy3 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                                muy4 * (u(3, i, j + 2, k) - u(3, i, j, k))) +
                     strz(k) * ((2 * muz1 + la(i, j, k - 1) * strz(k - 1) -
                                 tf * (la(i, j, k) * strz(k) +
                                       la(i, j, k - 2) * strz(k - 2))) *
                                    (u(3, i, j, k - 2) - u(3, i, j, k)) +
                                (2 * muz2 + la(i, j, k - 2) * strz(k - 2) +
                                 la(i, j, k + 1) * strz(k + 1) +
                                 3 * (la(i, j, k) * strz(k) +
                                      la(i, j, k - 1) * strz(k - 1))) *
                                    (u(3, i, j, k - 1) - u(3, i, j, k)) +
                                (2 * muz3 + la(i, j, k - 1) * strz(k - 1) +
                                 la(i, j, k + 2) * strz(k + 2) +
                                 3 * (la(i, j, k + 1) * strz(k + 1) +
                                      la(i, j, k) * strz(k))) *
                                    (u(3, i, j, k + 1) - u(3, i, j, k)) +
                                (2 * muz4 + la(i, j, k + 1) * strz(k + 1) -
                                 tf * (la(i, j, k) * strz(k) +
                                       la(i, j, k + 2) * strz(k + 2))) *
                                    (u(3, i, j, k + 2) - u(3, i, j, k))));

          /* Mixed derivatives: */
          /* 29ops /mixed derivative */
          /* 116 ops for r1 */
          /*   (la*v_y)_x */
          r1 = r1 +
               strx(i) * stry(j) * i144 *
                   (la(i - 2, j, k) *
                        (u(2, i - 2, j - 2, k) - u(2, i - 2, j + 2, k) +
                         8 * (-u(2, i - 2, j - 1, k) + u(2, i - 2, j + 1, k))) -
                    8 * (la(i - 1, j, k) *
                         (u(2, i - 1, j - 2, k) - u(2, i - 1, j + 2, k) +
                          8 * (-u(2, i - 1, j - 1, k) +
                               u(2, i - 1, j + 1, k)))) +
                    8 * (la(i + 1, j, k) *
                         (u(2, i + 1, j - 2, k) - u(2, i + 1, j + 2, k) +
                          8 * (-u(2, i + 1, j - 1, k) +
                               u(2, i + 1, j + 1, k)))) -
                    (la(i + 2, j, k) *
                     (u(2, i + 2, j - 2, k) - u(2, i + 2, j + 2, k) +
                      8 * (-u(2, i + 2, j - 1, k) + u(2, i + 2, j + 1, k)))))
               /*   (la*w_z)_x */
               +
               strx(i) * strz(k) * i144 *
                   (la(i - 2, j, k) *
                        (u(3, i - 2, j, k - 2) - u(3, i - 2, j, k + 2) +
                         8 * (-u(3, i - 2, j, k - 1) + u(3, i - 2, j, k + 1))) -
                    8 * (la(i - 1, j, k) *
                         (u(3, i - 1, j, k - 2) - u(3, i - 1, j, k + 2) +
                          8 * (-u(3, i - 1, j, k - 1) +
                               u(3, i - 1, j, k + 1)))) +
                    8 * (la(i + 1, j, k) *
                         (u(3, i + 1, j, k - 2) - u(3, i + 1, j, k + 2) +
                          8 * (-u(3, i + 1, j, k - 1) +
                               u(3, i + 1, j, k + 1)))) -
                    (la(i + 2, j, k) *
                     (u(3, i + 2, j, k - 2) - u(3, i + 2, j, k + 2) +
                      8 * (-u(3, i + 2, j, k - 1) + u(3, i + 2, j, k + 1)))))
               /*   (mu*v_x)_y */
               +
               strx(i) * stry(j) * i144 *
                   (mu(i, j - 2, k) *
                        (u(2, i - 2, j - 2, k) - u(2, i + 2, j - 2, k) +
                         8 * (-u(2, i - 1, j - 2, k) + u(2, i + 1, j - 2, k))) -
                    8 * (mu(i, j - 1, k) *
                         (u(2, i - 2, j - 1, k) - u(2, i + 2, j - 1, k) +
                          8 * (-u(2, i - 1, j - 1, k) +
                               u(2, i + 1, j - 1, k)))) +
                    8 * (mu(i, j + 1, k) *
                         (u(2, i - 2, j + 1, k) - u(2, i + 2, j + 1, k) +
                          8 * (-u(2, i - 1, j + 1, k) +
                               u(2, i + 1, j + 1, k)))) -
                    (mu(i, j + 2, k) *
                     (u(2, i - 2, j + 2, k) - u(2, i + 2, j + 2, k) +
                      8 * (-u(2, i - 1, j + 2, k) + u(2, i + 1, j + 2, k)))))
               /*   (mu*w_x)_z */
               +
               strx(i) * strz(k) * i144 *
                   (mu(i, j, k - 2) *
                        (u(3, i - 2, j, k - 2) - u(3, i + 2, j, k - 2) +
                         8 * (-u(3, i - 1, j, k - 2) + u(3, i + 1, j, k - 2))) -
                    8 * (mu(i, j, k - 1) *
                         (u(3, i - 2, j, k - 1) - u(3, i + 2, j, k - 1) +
                          8 * (-u(3, i - 1, j, k - 1) +
                               u(3, i + 1, j, k - 1)))) +
                    8 * (mu(i, j, k + 1) *
                         (u(3, i - 2, j, k + 1) - u(3, i + 2, j, k + 1) +
                          8 * (-u(3, i - 1, j, k + 1) +
                               u(3, i + 1, j, k + 1)))) -
                    (mu(i, j, k + 2) *
                     (u(3, i - 2, j, k + 2) - u(3, i + 2, j, k + 2) +
                      8 * (-u(3, i - 1, j, k + 2) + u(3, i + 1, j, k + 2)))));

          /* 116 ops for r2 */
          /*   (mu*u_y)_x */
          r2 = r2 +
               strx(i) * stry(j) * i144 *
                   (mu(i - 2, j, k) *
                        (u(1, i - 2, j - 2, k) - u(1, i - 2, j + 2, k) +
                         8 * (-u(1, i - 2, j - 1, k) + u(1, i - 2, j + 1, k))) -
                    8 * (mu(i - 1, j, k) *
                         (u(1, i - 1, j - 2, k) - u(1, i - 1, j + 2, k) +
                          8 * (-u(1, i - 1, j - 1, k) +
                               u(1, i - 1, j + 1, k)))) +
                    8 * (mu(i + 1, j, k) *
                         (u(1, i + 1, j - 2, k) - u(1, i + 1, j + 2, k) +
                          8 * (-u(1, i + 1, j - 1, k) +
                               u(1, i + 1, j + 1, k)))) -
                    (mu(i + 2, j, k) *
                     (u(1, i + 2, j - 2, k) - u(1, i + 2, j + 2, k) +
                      8 * (-u(1, i + 2, j - 1, k) + u(1, i + 2, j + 1, k)))))
               /* (la*u_x)_y */
               +
               strx(i) * stry(j) * i144 *
                   (la(i, j - 2, k) *
                        (u(1, i - 2, j - 2, k) - u(1, i + 2, j - 2, k) +
                         8 * (-u(1, i - 1, j - 2, k) + u(1, i + 1, j - 2, k))) -
                    8 * (la(i, j - 1, k) *
                         (u(1, i - 2, j - 1, k) - u(1, i + 2, j - 1, k) +
                          8 * (-u(1, i - 1, j - 1, k) +
                               u(1, i + 1, j - 1, k)))) +
                    8 * (la(i, j + 1, k) *
                         (u(1, i - 2, j + 1, k) - u(1, i + 2, j + 1, k) +
                          8 * (-u(1, i - 1, j + 1, k) +
                               u(1, i + 1, j + 1, k)))) -
                    (la(i, j + 2, k) *
                     (u(1, i - 2, j + 2, k) - u(1, i + 2, j + 2, k) +
                      8 * (-u(1, i - 1, j + 2, k) + u(1, i + 1, j + 2, k)))))
               /* (la*w_z)_y */
               +
               stry(j) * strz(k) * i144 *
                   (la(i, j - 2, k) *
                        (u(3, i, j - 2, k - 2) - u(3, i, j - 2, k + 2) +
                         8 * (-u(3, i, j - 2, k - 1) + u(3, i, j - 2, k + 1))) -
                    8 * (la(i, j - 1, k) *
                         (u(3, i, j - 1, k - 2) - u(3, i, j - 1, k + 2) +
                          8 * (-u(3, i, j - 1, k - 1) +
                               u(3, i, j - 1, k + 1)))) +
                    8 * (la(i, j + 1, k) *
                         (u(3, i, j + 1, k - 2) - u(3, i, j + 1, k + 2) +
                          8 * (-u(3, i, j + 1, k - 1) +
                               u(3, i, j + 1, k + 1)))) -
                    (la(i, j + 2, k) *
                     (u(3, i, j + 2, k - 2) - u(3, i, j + 2, k + 2) +
                      8 * (-u(3, i, j + 2, k - 1) + u(3, i, j + 2, k + 1)))))
               /* (mu*w_y)_z */
               +
               stry(j) * strz(k) * i144 *
                   (mu(i, j, k - 2) *
                        (u(3, i, j - 2, k - 2) - u(3, i, j + 2, k - 2) +
                         8 * (-u(3, i, j - 1, k - 2) + u(3, i, j + 1, k - 2))) -
                    8 * (mu(i, j, k - 1) *
                         (u(3, i, j - 2, k - 1) - u(3, i, j + 2, k - 1) +
                          8 * (-u(3, i, j - 1, k - 1) +
                               u(3, i, j + 1, k - 1)))) +
                    8 * (mu(i, j, k + 1) *
                         (u(3, i, j - 2, k + 1) - u(3, i, j + 2, k + 1) +
                          8 * (-u(3, i, j - 1, k + 1) +
                               u(3, i, j + 1, k + 1)))) -
                    (mu(i, j, k + 2) *
                     (u(3, i, j - 2, k + 2) - u(3, i, j + 2, k + 2) +
                      8 * (-u(3, i, j - 1, k + 2) + u(3, i, j + 1, k + 2)))));
          /* 116 ops for r3 */
          /*  (mu*u_z)_x */
          r3 = r3 +
               strx(i) * strz(k) * i144 *
                   (mu(i - 2, j, k) *
                        (u(1, i - 2, j, k - 2) - u(1, i - 2, j, k + 2) +
                         8 * (-u(1, i - 2, j, k - 1) + u(1, i - 2, j, k + 1))) -
                    8 * (mu(i - 1, j, k) *
                         (u(1, i - 1, j, k - 2) - u(1, i - 1, j, k + 2) +
                          8 * (-u(1, i - 1, j, k - 1) +
                               u(1, i - 1, j, k + 1)))) +
                    8 * (mu(i + 1, j, k) *
                         (u(1, i + 1, j, k - 2) - u(1, i + 1, j, k + 2) +
                          8 * (-u(1, i + 1, j, k - 1) +
                               u(1, i + 1, j, k + 1)))) -
                    (mu(i + 2, j, k) *
                     (u(1, i + 2, j, k - 2) - u(1, i + 2, j, k + 2) +
                      8 * (-u(1, i + 2, j, k - 1) + u(1, i + 2, j, k + 1)))))
               /* (mu*v_z)_y */
               +
               stry(j) * strz(k) * i144 *
                   (mu(i, j - 2, k) *
                        (u(2, i, j - 2, k - 2) - u(2, i, j - 2, k + 2) +
                         8 * (-u(2, i, j - 2, k - 1) + u(2, i, j - 2, k + 1))) -
                    8 * (mu(i, j - 1, k) *
                         (u(2, i, j - 1, k - 2) - u(2, i, j - 1, k + 2) +
                          8 * (-u(2, i, j - 1, k - 1) +
                               u(2, i, j - 1, k + 1)))) +
                    8 * (mu(i, j + 1, k) *
                         (u(2, i, j + 1, k - 2) - u(2, i, j + 1, k + 2) +
                          8 * (-u(2, i, j + 1, k - 1) +
                               u(2, i, j + 1, k + 1)))) -
                    (mu(i, j + 2, k) *
                     (u(2, i, j + 2, k - 2) - u(2, i, j + 2, k + 2) +
                      8 * (-u(2, i, j + 2, k - 1) + u(2, i, j + 2, k + 1)))))
               /*   (la*u_x)_z */
               +
               strx(i) * strz(k) * i144 *
                   (la(i, j, k - 2) *
                        (u(1, i - 2, j, k - 2) - u(1, i + 2, j, k - 2) +
                         8 * (-u(1, i - 1, j, k - 2) + u(1, i + 1, j, k - 2))) -
                    8 * (la(i, j, k - 1) *
                         (u(1, i - 2, j, k - 1) - u(1, i + 2, j, k - 1) +
                          8 * (-u(1, i - 1, j, k - 1) +
                               u(1, i + 1, j, k - 1)))) +
                    8 * (la(i, j, k + 1) *
                         (u(1, i - 2, j, k + 1) - u(1, i + 2, j, k + 1) +
                          8 * (-u(1, i - 1, j, k + 1) +
                               u(1, i + 1, j, k + 1)))) -
                    (la(i, j, k + 2) *
                     (u(1, i - 2, j, k + 2) - u(1, i + 2, j, k + 2) +
                      8 * (-u(1, i - 1, j, k + 2) + u(1, i + 1, j, k + 2)))))
               /* (la*v_y)_z */
               +
               stry(j) * strz(k) * i144 *
                   (la(i, j, k - 2) *
                        (u(2, i, j - 2, k - 2) - u(2, i, j + 2, k - 2) +
                         8 * (-u(2, i, j - 1, k - 2) + u(2, i, j + 1, k - 2))) -
                    8 * (la(i, j, k - 1) *
                         (u(2, i, j - 2, k - 1) - u(2, i, j + 2, k - 1) +
                          8 * (-u(2, i, j - 1, k - 1) +
                               u(2, i, j + 1, k - 1)))) +
                    8 * (la(i, j, k + 1) *
                         (u(2, i, j - 2, k + 1) - u(2, i, j + 2, k + 1) +
                          8 * (-u(2, i, j - 1, k + 1) +
                               u(2, i, j + 1, k + 1)))) -
                    (la(i, j, k + 2) *
                     (u(2, i, j - 2, k + 2) - u(2, i, j + 2, k + 2) +
                      8 * (-u(2, i, j - 1, k + 2) + u(2, i, j + 1, k + 2)))));

          /* 9 ops */
          lu(1, i, j, k) = a1 * lu(1, i, j, k) + cof * r1;
          lu(2, i, j, k) = a1 * lu(2, i, j, k) + cof * r2;
          lu(3, i, j, k) = a1 * lu(3, i, j, k) + cof * r3;
        });  // End of rhs4th3fort_ci LOOP 1
    SYNC_STREAM;
    if (onesided[4] == 1) {
      RAJA::RangeSegment k_range(1, 6 + 1);
      RAJA::RangeSegment j_range(jfirst + 2, jlast - 1);
      RAJA::RangeSegment i_range(ifirst + 2, ilast - 1);
      RAJA::kernel<RHS4_EXEC_POL>(
          RAJA::make_tuple(k_range, j_range, i_range),
          [=] RAJA_DEVICE(int k, int j, int i) {
            float_sw4 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4;
            float_sw4 r1, r2, r3, mu1zz, mu2zz, mu3zz, lap2mu, mucof;
            /* from inner_loop_4a */
            mux1 = mu(i - 1, j, k) * strx(i - 1) -
                   tf * (mu(i, j, k) * strx(i) + mu(i - 2, j, k) * strx(i - 2));
            mux2 = mu(i - 2, j, k) * strx(i - 2) +
                   mu(i + 1, j, k) * strx(i + 1) +
                   3 * (mu(i, j, k) * strx(i) + mu(i - 1, j, k) * strx(i - 1));
            mux3 = mu(i - 1, j, k) * strx(i - 1) +
                   mu(i + 2, j, k) * strx(i + 2) +
                   3 * (mu(i + 1, j, k) * strx(i + 1) + mu(i, j, k) * strx(i));
            mux4 = mu(i + 1, j, k) * strx(i + 1) -
                   tf * (mu(i, j, k) * strx(i) + mu(i + 2, j, k) * strx(i + 2));

            muy1 = mu(i, j - 1, k) * stry(j - 1) -
                   tf * (mu(i, j, k) * stry(j) + mu(i, j - 2, k) * stry(j - 2));
            muy2 = mu(i, j - 2, k) * stry(j - 2) +
                   mu(i, j + 1, k) * stry(j + 1) +
                   3 * (mu(i, j, k) * stry(j) + mu(i, j - 1, k) * stry(j - 1));
            muy3 = mu(i, j - 1, k) * stry(j - 1) +
                   mu(i, j + 2, k) * stry(j + 2) +
                   3 * (mu(i, j + 1, k) * stry(j + 1) + mu(i, j, k) * stry(j));
            muy4 = mu(i, j + 1, k) * stry(j + 1) -
                   tf * (mu(i, j, k) * stry(j) + mu(i, j + 2, k) * stry(j + 2));

            r1 = i6 * (strx(i) * ((2 * mux1 + la(i - 1, j, k) * strx(i - 1) -
                                   tf * (la(i, j, k) * strx(i) +
                                         la(i - 2, j, k) * strx(i - 2))) *
                                      (u(1, i - 2, j, k) - u(1, i, j, k)) +
                                  (2 * mux2 + la(i - 2, j, k) * strx(i - 2) +
                                   la(i + 1, j, k) * strx(i + 1) +
                                   3 * (la(i, j, k) * strx(i) +
                                        la(i - 1, j, k) * strx(i - 1))) *
                                      (u(1, i - 1, j, k) - u(1, i, j, k)) +
                                  (2 * mux3 + la(i - 1, j, k) * strx(i - 1) +
                                   la(i + 2, j, k) * strx(i + 2) +
                                   3 * (la(i + 1, j, k) * strx(i + 1) +
                                        la(i, j, k) * strx(i))) *
                                      (u(1, i + 1, j, k) - u(1, i, j, k)) +
                                  (2 * mux4 + la(i + 1, j, k) * strx(i + 1) -
                                   tf * (la(i, j, k) * strx(i) +
                                         la(i + 2, j, k) * strx(i + 2))) *
                                      (u(1, i + 2, j, k) - u(1, i, j, k))) +
                       stry(j) * (+muy1 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                                  muy2 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                                  muy3 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                                  muy4 * (u(1, i, j + 2, k) - u(1, i, j, k))));

            /* (mu*uz)_z can not be centered */
            /* second derivative (mu*u_z)_z at grid point z_k */
            /* averaging the coefficient, */
            /* leave out the z-supergrid stretching strz, since it will */
            /* never be used together with the sbp-boundary operator */
            mu1zz = 0;
            mu2zz = 0;
            mu3zz = 0;
            for (int q = 1; q <= 8; q++) {
              //		     lap2mu= 0;
              //		     mucof = 0;
              //		     for( m=1 ; m<=8; m++ )
              //		     {
              //			mucof  += acof(k,q,m)*mu(i,j,m);
              //			lap2mu +=
              // acof(k,q,m)*(la(i,j,m)+2*mu(i,j,m));
              //		     }
              lap2mu = acof(k, q, 1) * (la(i, j, 1) + 2 * mu(i, j, 1)) +
                       acof(k, q, 2) * (la(i, j, 2) + 2 * mu(i, j, 2)) +
                       acof(k, q, 3) * (la(i, j, 3) + 2 * mu(i, j, 3)) +
                       acof(k, q, 4) * (la(i, j, 4) + 2 * mu(i, j, 4)) +
                       acof(k, q, 5) * (la(i, j, 5) + 2 * mu(i, j, 5)) +
                       acof(k, q, 6) * (la(i, j, 6) + 2 * mu(i, j, 6)) +
                       acof(k, q, 7) * (la(i, j, 7) + 2 * mu(i, j, 7)) +
                       acof(k, q, 8) * (la(i, j, 8) + 2 * mu(i, j, 8));
              mucof =
                  acof(k, q, 1) * mu(i, j, 1) + acof(k, q, 2) * mu(i, j, 2) +
                  acof(k, q, 3) * mu(i, j, 3) + acof(k, q, 4) * mu(i, j, 4) +
                  acof(k, q, 5) * mu(i, j, 5) + acof(k, q, 6) * mu(i, j, 6) +
                  acof(k, q, 7) * mu(i, j, 7) + acof(k, q, 8) * mu(i, j, 8);
              mu1zz += mucof * u(1, i, j, q);
              mu2zz += mucof * u(2, i, j, q);
              mu3zz += lap2mu * u(3, i, j, q);
            }

            /* ghost point only influences the first point (k=1) because
             * ghcof(k)=0 for k>=2*/
            r1 = r1 + (mu1zz + ghcof(k) * mu(i, j, 1) * u(1, i, j, 0));

            r2 = i6 * (strx(i) * (mux1 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                                  mux2 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                                  mux3 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                                  mux4 * (u(2, i + 2, j, k) - u(2, i, j, k))) +
                       stry(j) * ((2 * muy1 + la(i, j - 1, k) * stry(j - 1) -
                                   tf * (la(i, j, k) * stry(j) +
                                         la(i, j - 2, k) * stry(j - 2))) *
                                      (u(2, i, j - 2, k) - u(2, i, j, k)) +
                                  (2 * muy2 + la(i, j - 2, k) * stry(j - 2) +
                                   la(i, j + 1, k) * stry(j + 1) +
                                   3 * (la(i, j, k) * stry(j) +
                                        la(i, j - 1, k) * stry(j - 1))) *
                                      (u(2, i, j - 1, k) - u(2, i, j, k)) +
                                  (2 * muy3 + la(i, j - 1, k) * stry(j - 1) +
                                   la(i, j + 2, k) * stry(j + 2) +
                                   3 * (la(i, j + 1, k) * stry(j + 1) +
                                        la(i, j, k) * stry(j))) *
                                      (u(2, i, j + 1, k) - u(2, i, j, k)) +
                                  (2 * muy4 + la(i, j + 1, k) * stry(j + 1) -
                                   tf * (la(i, j, k) * stry(j) +
                                         la(i, j + 2, k) * stry(j + 2))) *
                                      (u(2, i, j + 2, k) - u(2, i, j, k))));

            /* ghost point only influences the first point (k=1) because
             * ghcof(k)=0 for k>=2 */
            r2 = r2 + (mu2zz + ghcof(k) * mu(i, j, 1) * u(2, i, j, 0));

            r3 = i6 * (strx(i) * (mux1 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                                  mux2 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                                  mux3 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                                  mux4 * (u(3, i + 2, j, k) - u(3, i, j, k))) +
                       stry(j) * (muy1 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                                  muy2 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                                  muy3 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                                  muy4 * (u(3, i, j + 2, k) - u(3, i, j, k))));
            /* ghost point only influences the first point (k=1) because
             * ghcof(k)=0 for k>=2 */
            r3 = r3 + (mu3zz + ghcof(k) * (la(i, j, 1) + 2 * mu(i, j, 1)) *
                                   u(3, i, j, 0));

            /* cross-terms in first component of rhs */
            /*   (la*v_y)_x */
            r1 = r1 +
                 strx(i) * stry(j) *
                     (i144 *
                          (la(i - 2, j, k) *
                               (u(2, i - 2, j - 2, k) - u(2, i - 2, j + 2, k) +
                                8 * (-u(2, i - 2, j - 1, k) +
                                     u(2, i - 2, j + 1, k))) -
                           8 * (la(i - 1, j, k) *
                                (u(2, i - 1, j - 2, k) - u(2, i - 1, j + 2, k) +
                                 8 * (-u(2, i - 1, j - 1, k) +
                                      u(2, i - 1, j + 1, k)))) +
                           8 * (la(i + 1, j, k) *
                                (u(2, i + 1, j - 2, k) - u(2, i + 1, j + 2, k) +
                                 8 * (-u(2, i + 1, j - 1, k) +
                                      u(2, i + 1, j + 1, k)))) -
                           (la(i + 2, j, k) *
                            (u(2, i + 2, j - 2, k) - u(2, i + 2, j + 2, k) +
                             8 * (-u(2, i + 2, j - 1, k) +
                                  u(2, i + 2, j + 1, k)))))
                      /*   (mu*v_x)_y */
                      +
                      i144 *
                          (mu(i, j - 2, k) *
                               (u(2, i - 2, j - 2, k) - u(2, i + 2, j - 2, k) +
                                8 * (-u(2, i - 1, j - 2, k) +
                                     u(2, i + 1, j - 2, k))) -
                           8 * (mu(i, j - 1, k) *
                                (u(2, i - 2, j - 1, k) - u(2, i + 2, j - 1, k) +
                                 8 * (-u(2, i - 1, j - 1, k) +
                                      u(2, i + 1, j - 1, k)))) +
                           8 * (mu(i, j + 1, k) *
                                (u(2, i - 2, j + 1, k) - u(2, i + 2, j + 1, k) +
                                 8 * (-u(2, i - 1, j + 1, k) +
                                      u(2, i + 1, j + 1, k)))) -
                           (mu(i, j + 2, k) *
                            (u(2, i - 2, j + 2, k) - u(2, i + 2, j + 2, k) +
                             8 * (-u(2, i - 1, j + 2, k) +
                                  u(2, i + 1, j + 2, k))))));
            /*   (la*w_z)_x: NOT CENTERED */
            float_sw4 u3zip2 = 0;
            float_sw4 u3zip1 = 0;
            float_sw4 u3zim1 = 0;
            float_sw4 u3zim2 = 0;
            for (int q = 1; q <= 8; q++) {
              u3zip2 += bope(k, q) * u(3, i + 2, j, q);
              u3zip1 += bope(k, q) * u(3, i + 1, j, q);
              u3zim1 += bope(k, q) * u(3, i - 1, j, q);
              u3zim2 += bope(k, q) * u(3, i - 2, j, q);
            }
            float_sw4 lau3zx =
                i12 *
                (-la(i + 2, j, k) * u3zip2 + 8 * la(i + 1, j, k) * u3zip1 -
                 8 * la(i - 1, j, k) * u3zim1 + la(i - 2, j, k) * u3zim2);
            r1 = r1 + strx(i) * lau3zx;
            /*   (mu*w_x)_z: NOT CENTERED */
            float_sw4 mu3xz = 0;
            for (int q = 1; q <= 8; q++)
              mu3xz +=
                  bope(k, q) * (mu(i, j, q) * i12 *
                                (-u(3, i + 2, j, q) + 8 * u(3, i + 1, j, q) -
                                 8 * u(3, i - 1, j, q) + u(3, i - 2, j, q)));
            r1 = r1 + strx(i) * mu3xz;

            /* cross-terms in second component of rhs */
            /*   (mu*u_y)_x */
            r2 = r2 +
                 strx(i) * stry(j) *
                     (i144 *
                          (mu(i - 2, j, k) *
                               (u(1, i - 2, j - 2, k) - u(1, i - 2, j + 2, k) +
                                8 * (-u(1, i - 2, j - 1, k) +
                                     u(1, i - 2, j + 1, k))) -
                           8 * (mu(i - 1, j, k) *
                                (u(1, i - 1, j - 2, k) - u(1, i - 1, j + 2, k) +
                                 8 * (-u(1, i - 1, j - 1, k) +
                                      u(1, i - 1, j + 1, k)))) +
                           8 * (mu(i + 1, j, k) *
                                (u(1, i + 1, j - 2, k) - u(1, i + 1, j + 2, k) +
                                 8 * (-u(1, i + 1, j - 1, k) +
                                      u(1, i + 1, j + 1, k)))) -
                           (mu(i + 2, j, k) *
                            (u(1, i + 2, j - 2, k) - u(1, i + 2, j + 2, k) +
                             8 * (-u(1, i + 2, j - 1, k) +
                                  u(1, i + 2, j + 1, k)))))
                      /* (la*u_x)_y  */
                      +
                      i144 *
                          (la(i, j - 2, k) *
                               (u(1, i - 2, j - 2, k) - u(1, i + 2, j - 2, k) +
                                8 * (-u(1, i - 1, j - 2, k) +
                                     u(1, i + 1, j - 2, k))) -
                           8 * (la(i, j - 1, k) *
                                (u(1, i - 2, j - 1, k) - u(1, i + 2, j - 1, k) +
                                 8 * (-u(1, i - 1, j - 1, k) +
                                      u(1, i + 1, j - 1, k)))) +
                           8 * (la(i, j + 1, k) *
                                (u(1, i - 2, j + 1, k) - u(1, i + 2, j + 1, k) +
                                 8 * (-u(1, i - 1, j + 1, k) +
                                      u(1, i + 1, j + 1, k)))) -
                           (la(i, j + 2, k) *
                            (u(1, i - 2, j + 2, k) - u(1, i + 2, j + 2, k) +
                             8 * (-u(1, i - 1, j + 2, k) +
                                  u(1, i + 1, j + 2, k))))));
            /* (la*w_z)_y : NOT CENTERED */
            float_sw4 u3zjp2 = 0;
            float_sw4 u3zjp1 = 0;
            float_sw4 u3zjm1 = 0;
            float_sw4 u3zjm2 = 0;
            for (int q = 1; q <= 8; q++) {
              u3zjp2 += bope(k, q) * u(3, i, j + 2, q);
              u3zjp1 += bope(k, q) * u(3, i, j + 1, q);
              u3zjm1 += bope(k, q) * u(3, i, j - 1, q);
              u3zjm2 += bope(k, q) * u(3, i, j - 2, q);
            }
            float_sw4 lau3zy =
                i12 *
                (-la(i, j + 2, k) * u3zjp2 + 8 * la(i, j + 1, k) * u3zjp1 -
                 8 * la(i, j - 1, k) * u3zjm1 + la(i, j - 2, k) * u3zjm2);

            r2 = r2 + stry(j) * lau3zy;

            /* (mu*w_y)_z: NOT CENTERED */
            float_sw4 mu3yz = 0;
            for (int q = 1; q <= 8; q++)
              mu3yz +=
                  bope(k, q) * (mu(i, j, q) * i12 *
                                (-u(3, i, j + 2, q) + 8 * u(3, i, j + 1, q) -
                                 8 * u(3, i, j - 1, q) + u(3, i, j - 2, q)));

            r2 = r2 + stry(j) * mu3yz;

            /* No centered cross terms in r3 */
            /*  (mu*u_z)_x: NOT CENTERED */
            float_sw4 u1zip2 = 0;
            float_sw4 u1zip1 = 0;
            float_sw4 u1zim1 = 0;
            float_sw4 u1zim2 = 0;
            for (int q = 1; q <= 8; q++) {
              u1zip2 += bope(k, q) * u(1, i + 2, j, q);
              u1zip1 += bope(k, q) * u(1, i + 1, j, q);
              u1zim1 += bope(k, q) * u(1, i - 1, j, q);
              u1zim2 += bope(k, q) * u(1, i - 2, j, q);
            }
            float_sw4 mu1zx =
                i12 *
                (-mu(i + 2, j, k) * u1zip2 + 8 * mu(i + 1, j, k) * u1zip1 -
                 8 * mu(i - 1, j, k) * u1zim1 + mu(i - 2, j, k) * u1zim2);
            r3 = r3 + strx(i) * mu1zx;

            /* (mu*v_z)_y: NOT CENTERED */
            float_sw4 u2zjp2 = 0;
            float_sw4 u2zjp1 = 0;
            float_sw4 u2zjm1 = 0;
            float_sw4 u2zjm2 = 0;
            for (int q = 1; q <= 8; q++) {
              u2zjp2 += bope(k, q) * u(2, i, j + 2, q);
              u2zjp1 += bope(k, q) * u(2, i, j + 1, q);
              u2zjm1 += bope(k, q) * u(2, i, j - 1, q);
              u2zjm2 += bope(k, q) * u(2, i, j - 2, q);
            }
            float_sw4 mu2zy =
                i12 *
                (-mu(i, j + 2, k) * u2zjp2 + 8 * mu(i, j + 1, k) * u2zjp1 -
                 8 * mu(i, j - 1, k) * u2zjm1 + mu(i, j - 2, k) * u2zjm2);
            r3 = r3 + stry(j) * mu2zy;

            /*   (la*u_x)_z: NOT CENTERED */
            float_sw4 lau1xz = 0;
            for (int q = 1; q <= 8; q++)
              lau1xz +=
                  bope(k, q) * (la(i, j, q) * i12 *
                                (-u(1, i + 2, j, q) + 8 * u(1, i + 1, j, q) -
                                 8 * u(1, i - 1, j, q) + u(1, i - 2, j, q)));
            r3 = r3 + strx(i) * lau1xz;

            /* (la*v_y)_z: NOT CENTERED */
            float_sw4 lau2yz = 0;
            for (int q = 1; q <= 8; q++)
              lau2yz +=
                  bope(k, q) * (la(i, j, q) * i12 *
                                (-u(2, i, j + 2, q) + 8 * u(2, i, j + 1, q) -
                                 8 * u(2, i, j - 1, q) + u(2, i, j - 2, q)));
            r3 = r3 + stry(j) * lau2yz;

            lu(1, i, j, k) = a1 * lu(1, i, j, k) + cof * r1;
            lu(2, i, j, k) = a1 * lu(2, i, j, k) + cof * r2;
            lu(3, i, j, k) = a1 * lu(3, i, j, k) + cof * r3;
          });  // End of rhs4th3fort_ci LOOP 2
      SYNC_STREAM;
    }
    if (onesided[5] == 1) {
      RAJA::RangeSegment k_range(nk - 5, nk + 1);
      RAJA::RangeSegment j_range(jfirst + 2, jlast - 1);
      RAJA::RangeSegment i_range(ifirst + 2, ilast - 1);
      RAJA::kernel<RHS4_EXEC_POL>(
          RAJA::make_tuple(k_range, j_range, i_range),
          [=] RAJA_DEVICE(int k, int j, int i) {
            float_sw4 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4;
            float_sw4 r1, r2, r3;
            // #pragma omp for
            // 	 for(  k = nk-5 ; k <= nk ; k++ )
            // 	    for(  j=jfirst+2; j<=jlast-2; j++ )
            // #pragma simd
            // #pragma ivdep
            // 	       for(  i=ifirst+2; i<=ilast-2; i++ )
            // 	       {
            /* from inner_loop_4a */
            mux1 = mu(i - 1, j, k) * strx(i - 1) -
                   tf * (mu(i, j, k) * strx(i) + mu(i - 2, j, k) * strx(i - 2));
            mux2 = mu(i - 2, j, k) * strx(i - 2) +
                   mu(i + 1, j, k) * strx(i + 1) +
                   3 * (mu(i, j, k) * strx(i) + mu(i - 1, j, k) * strx(i - 1));
            mux3 = mu(i - 1, j, k) * strx(i - 1) +
                   mu(i + 2, j, k) * strx(i + 2) +
                   3 * (mu(i + 1, j, k) * strx(i + 1) + mu(i, j, k) * strx(i));
            mux4 = mu(i + 1, j, k) * strx(i + 1) -
                   tf * (mu(i, j, k) * strx(i) + mu(i + 2, j, k) * strx(i + 2));

            muy1 = mu(i, j - 1, k) * stry(j - 1) -
                   tf * (mu(i, j, k) * stry(j) + mu(i, j - 2, k) * stry(j - 2));
            muy2 = mu(i, j - 2, k) * stry(j - 2) +
                   mu(i, j + 1, k) * stry(j + 1) +
                   3 * (mu(i, j, k) * stry(j) + mu(i, j - 1, k) * stry(j - 1));
            muy3 = mu(i, j - 1, k) * stry(j - 1) +
                   mu(i, j + 2, k) * stry(j + 2) +
                   3 * (mu(i, j + 1, k) * stry(j + 1) + mu(i, j, k) * stry(j));
            muy4 = mu(i, j + 1, k) * stry(j + 1) -
                   tf * (mu(i, j, k) * stry(j) + mu(i, j + 2, k) * stry(j + 2));

            /* xx, yy, and zz derivatives: */
            /* note that we could have introduced intermediate variables for the
             * average of lambda  */
            /* in the same way as we did for mu */
            r1 = i6 * (strx(i) * ((2 * mux1 + la(i - 1, j, k) * strx(i - 1) -
                                   tf * (la(i, j, k) * strx(i) +
                                         la(i - 2, j, k) * strx(i - 2))) *
                                      (u(1, i - 2, j, k) - u(1, i, j, k)) +
                                  (2 * mux2 + la(i - 2, j, k) * strx(i - 2) +
                                   la(i + 1, j, k) * strx(i + 1) +
                                   3 * (la(i, j, k) * strx(i) +
                                        la(i - 1, j, k) * strx(i - 1))) *
                                      (u(1, i - 1, j, k) - u(1, i, j, k)) +
                                  (2 * mux3 + la(i - 1, j, k) * strx(i - 1) +
                                   la(i + 2, j, k) * strx(i + 2) +
                                   3 * (la(i + 1, j, k) * strx(i + 1) +
                                        la(i, j, k) * strx(i))) *
                                      (u(1, i + 1, j, k) - u(1, i, j, k)) +
                                  (2 * mux4 + la(i + 1, j, k) * strx(i + 1) -
                                   tf * (la(i, j, k) * strx(i) +
                                         la(i + 2, j, k) * strx(i + 2))) *
                                      (u(1, i + 2, j, k) - u(1, i, j, k))) +
                       stry(j) * (+muy1 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                                  muy2 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                                  muy3 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                                  muy4 * (u(1, i, j + 2, k) - u(1, i, j, k))));

            /* all indices ending with 'b' are indices relative to the boundary,
             * going into the domain (1,2,3,...)*/
            int kb = nk - k + 1;
            /* all coefficient arrays (acof, bope, ghcof) should be indexed with
             * these indices */
            /* all solution and material property arrays should be indexed with
             * (i,j,k) */

            /* (mu*uz)_z can not be centered */
            /* second derivative (mu*u_z)_z at grid point z_k */
            /* averaging the coefficient */
            float_sw4 mu1zz = 0;
            float_sw4 mu2zz = 0;
            float_sw4 mu3zz = 0;
            for (int qb = 1; qb <= 8; qb++) {
              float_sw4 mucof = 0;
              float_sw4 lap2mu = 0;
              for (int mb = 1; mb <= 8; mb++) {
                mucof += acof(kb, qb, mb) * mu(i, j, nk - mb + 1);
                lap2mu += acof(kb, qb, mb) *
                          (2 * mu(i, j, nk - mb + 1) + la(i, j, nk - mb + 1));
              }
              mu1zz += mucof * u(1, i, j, nk - qb + 1);
              mu2zz += mucof * u(2, i, j, nk - qb + 1);
              mu3zz += lap2mu * u(3, i, j, nk - qb + 1);
            }
            /* computing the second derivative */
            /* ghost point only influences the first point (k=1) because
             * ghcof(k)=0 for k>=2*/
            r1 = r1 + (mu1zz + ghcof(kb) * mu(i, j, nk) * u(1, i, j, nk + 1));

            r2 = i6 * (strx(i) * (mux1 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                                  mux2 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                                  mux3 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                                  mux4 * (u(2, i + 2, j, k) - u(2, i, j, k))) +
                       stry(j) * ((2 * muy1 + la(i, j - 1, k) * stry(j - 1) -
                                   tf * (la(i, j, k) * stry(j) +
                                         la(i, j - 2, k) * stry(j - 2))) *
                                      (u(2, i, j - 2, k) - u(2, i, j, k)) +
                                  (2 * muy2 + la(i, j - 2, k) * stry(j - 2) +
                                   la(i, j + 1, k) * stry(j + 1) +
                                   3 * (la(i, j, k) * stry(j) +
                                        la(i, j - 1, k) * stry(j - 1))) *
                                      (u(2, i, j - 1, k) - u(2, i, j, k)) +
                                  (2 * muy3 + la(i, j - 1, k) * stry(j - 1) +
                                   la(i, j + 2, k) * stry(j + 2) +
                                   3 * (la(i, j + 1, k) * stry(j + 1) +
                                        la(i, j, k) * stry(j))) *
                                      (u(2, i, j + 1, k) - u(2, i, j, k)) +
                                  (2 * muy4 + la(i, j + 1, k) * stry(j + 1) -
                                   tf * (la(i, j, k) * stry(j) +
                                         la(i, j + 2, k) * stry(j + 2))) *
                                      (u(2, i, j + 2, k) - u(2, i, j, k))));

            /* (mu*vz)_z can not be centered */
            /* second derivative (mu*v_z)_z at grid point z_k */
            /* averaging the coefficient: already done above */
            r2 = r2 + (mu2zz + ghcof(kb) * mu(i, j, nk) * u(2, i, j, nk + 1));

            r3 = i6 * (strx(i) * (mux1 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                                  mux2 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                                  mux3 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                                  mux4 * (u(3, i + 2, j, k) - u(3, i, j, k))) +
                       stry(j) * (muy1 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                                  muy2 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                                  muy3 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                                  muy4 * (u(3, i, j + 2, k) - u(3, i, j, k))));
            r3 = r3 + (mu3zz + ghcof(kb) * (la(i, j, nk) + 2 * mu(i, j, nk)) *
                                   u(3, i, j, nk + 1));

            /* cross-terms in first component of rhs */
            /*   (la*v_y)_x */
            r1 = r1 +
                 strx(i) * stry(j) *
                     (i144 *
                          (la(i - 2, j, k) *
                               (u(2, i - 2, j - 2, k) - u(2, i - 2, j + 2, k) +
                                8 * (-u(2, i - 2, j - 1, k) +
                                     u(2, i - 2, j + 1, k))) -
                           8 * (la(i - 1, j, k) *
                                (u(2, i - 1, j - 2, k) - u(2, i - 1, j + 2, k) +
                                 8 * (-u(2, i - 1, j - 1, k) +
                                      u(2, i - 1, j + 1, k)))) +
                           8 * (la(i + 1, j, k) *
                                (u(2, i + 1, j - 2, k) - u(2, i + 1, j + 2, k) +
                                 8 * (-u(2, i + 1, j - 1, k) +
                                      u(2, i + 1, j + 1, k)))) -
                           (la(i + 2, j, k) *
                            (u(2, i + 2, j - 2, k) - u(2, i + 2, j + 2, k) +
                             8 * (-u(2, i + 2, j - 1, k) +
                                  u(2, i + 2, j + 1, k)))))
                      /*   (mu*v_x)_y */
                      +
                      i144 *
                          (mu(i, j - 2, k) *
                               (u(2, i - 2, j - 2, k) - u(2, i + 2, j - 2, k) +
                                8 * (-u(2, i - 1, j - 2, k) +
                                     u(2, i + 1, j - 2, k))) -
                           8 * (mu(i, j - 1, k) *
                                (u(2, i - 2, j - 1, k) - u(2, i + 2, j - 1, k) +
                                 8 * (-u(2, i - 1, j - 1, k) +
                                      u(2, i + 1, j - 1, k)))) +
                           8 * (mu(i, j + 1, k) *
                                (u(2, i - 2, j + 1, k) - u(2, i + 2, j + 1, k) +
                                 8 * (-u(2, i - 1, j + 1, k) +
                                      u(2, i + 1, j + 1, k)))) -
                           (mu(i, j + 2, k) *
                            (u(2, i - 2, j + 2, k) - u(2, i + 2, j + 2, k) +
                             8 * (-u(2, i - 1, j + 2, k) +
                                  u(2, i + 1, j + 2, k))))));
            /*   (la*w_z)_x: NOT CENTERED */
            float_sw4 u3zip2 = 0;
            float_sw4 u3zip1 = 0;
            float_sw4 u3zim1 = 0;
            float_sw4 u3zim2 = 0;
            for (int qb = 1; qb <= 8; qb++) {
              u3zip2 -= bope(kb, qb) * u(3, i + 2, j, nk - qb + 1);
              u3zip1 -= bope(kb, qb) * u(3, i + 1, j, nk - qb + 1);
              u3zim1 -= bope(kb, qb) * u(3, i - 1, j, nk - qb + 1);
              u3zim2 -= bope(kb, qb) * u(3, i - 2, j, nk - qb + 1);
            }
            float_sw4 lau3zx =
                i12 *
                (-la(i + 2, j, k) * u3zip2 + 8 * la(i + 1, j, k) * u3zip1 -
                 8 * la(i - 1, j, k) * u3zim1 + la(i - 2, j, k) * u3zim2);
            r1 = r1 + strx(i) * lau3zx;

            /*   (mu*w_x)_z: NOT CENTERED */
            float_sw4 mu3xz = 0;
            for (int qb = 1; qb <= 8; qb++)
              mu3xz -= bope(kb, qb) * (mu(i, j, nk - qb + 1) * i12 *
                                       (-u(3, i + 2, j, nk - qb + 1) +
                                        8 * u(3, i + 1, j, nk - qb + 1) -
                                        8 * u(3, i - 1, j, nk - qb + 1) +
                                        u(3, i - 2, j, nk - qb + 1)));

            r1 = r1 + strx(i) * mu3xz;

            /* cross-terms in second component of rhs */
            /*   (mu*u_y)_x */
            r2 = r2 +
                 strx(i) * stry(j) *
                     (i144 *
                          (mu(i - 2, j, k) *
                               (u(1, i - 2, j - 2, k) - u(1, i - 2, j + 2, k) +
                                8 * (-u(1, i - 2, j - 1, k) +
                                     u(1, i - 2, j + 1, k))) -
                           8 * (mu(i - 1, j, k) *
                                (u(1, i - 1, j - 2, k) - u(1, i - 1, j + 2, k) +
                                 8 * (-u(1, i - 1, j - 1, k) +
                                      u(1, i - 1, j + 1, k)))) +
                           8 * (mu(i + 1, j, k) *
                                (u(1, i + 1, j - 2, k) - u(1, i + 1, j + 2, k) +
                                 8 * (-u(1, i + 1, j - 1, k) +
                                      u(1, i + 1, j + 1, k)))) -
                           (mu(i + 2, j, k) *
                            (u(1, i + 2, j - 2, k) - u(1, i + 2, j + 2, k) +
                             8 * (-u(1, i + 2, j - 1, k) +
                                  u(1, i + 2, j + 1, k)))))
                      /* (la*u_x)_y */
                      +
                      i144 *
                          (la(i, j - 2, k) *
                               (u(1, i - 2, j - 2, k) - u(1, i + 2, j - 2, k) +
                                8 * (-u(1, i - 1, j - 2, k) +
                                     u(1, i + 1, j - 2, k))) -
                           8 * (la(i, j - 1, k) *
                                (u(1, i - 2, j - 1, k) - u(1, i + 2, j - 1, k) +
                                 8 * (-u(1, i - 1, j - 1, k) +
                                      u(1, i + 1, j - 1, k)))) +
                           8 * (la(i, j + 1, k) *
                                (u(1, i - 2, j + 1, k) - u(1, i + 2, j + 1, k) +
                                 8 * (-u(1, i - 1, j + 1, k) +
                                      u(1, i + 1, j + 1, k)))) -
                           (la(i, j + 2, k) *
                            (u(1, i - 2, j + 2, k) - u(1, i + 2, j + 2, k) +
                             8 * (-u(1, i - 1, j + 2, k) +
                                  u(1, i + 1, j + 2, k))))));
            /* (la*w_z)_y : NOT CENTERED */
            float_sw4 u3zjp2 = 0;
            float_sw4 u3zjp1 = 0;
            float_sw4 u3zjm1 = 0;
            float_sw4 u3zjm2 = 0;
            for (int qb = 1; qb <= 8; qb++) {
              u3zjp2 -= bope(kb, qb) * u(3, i, j + 2, nk - qb + 1);
              u3zjp1 -= bope(kb, qb) * u(3, i, j + 1, nk - qb + 1);
              u3zjm1 -= bope(kb, qb) * u(3, i, j - 1, nk - qb + 1);
              u3zjm2 -= bope(kb, qb) * u(3, i, j - 2, nk - qb + 1);
            }
            float_sw4 lau3zy =
                i12 *
                (-la(i, j + 2, k) * u3zjp2 + 8 * la(i, j + 1, k) * u3zjp1 -
                 8 * la(i, j - 1, k) * u3zjm1 + la(i, j - 2, k) * u3zjm2);
            r2 = r2 + stry(j) * lau3zy;

            /* (mu*w_y)_z: NOT CENTERED */
            float_sw4 mu3yz = 0;
            for (int qb = 1; qb <= 8; qb++)
              mu3yz -= bope(kb, qb) * (mu(i, j, nk - qb + 1) * i12 *
                                       (-u(3, i, j + 2, nk - qb + 1) +
                                        8 * u(3, i, j + 1, nk - qb + 1) -
                                        8 * u(3, i, j - 1, nk - qb + 1) +
                                        u(3, i, j - 2, nk - qb + 1)));
            r2 = r2 + stry(j) * mu3yz;

            /* No centered cross terms in r3 */
            /*  (mu*u_z)_x: NOT CENTERED */
            float_sw4 u1zip2 = 0;
            float_sw4 u1zip1 = 0;
            float_sw4 u1zim1 = 0;
            float_sw4 u1zim2 = 0;
            for (int qb = 1; qb <= 8; qb++) {
              u1zip2 -= bope(kb, qb) * u(1, i + 2, j, nk - qb + 1);
              u1zip1 -= bope(kb, qb) * u(1, i + 1, j, nk - qb + 1);
              u1zim1 -= bope(kb, qb) * u(1, i - 1, j, nk - qb + 1);
              u1zim2 -= bope(kb, qb) * u(1, i - 2, j, nk - qb + 1);
            }
            float_sw4 mu1zx =
                i12 *
                (-mu(i + 2, j, k) * u1zip2 + 8 * mu(i + 1, j, k) * u1zip1 -
                 8 * mu(i - 1, j, k) * u1zim1 + mu(i - 2, j, k) * u1zim2);
            r3 = r3 + strx(i) * mu1zx;

            /* (mu*v_z)_y: NOT CENTERED */
            float_sw4 u2zjp2 = 0;
            float_sw4 u2zjp1 = 0;
            float_sw4 u2zjm1 = 0;
            float_sw4 u2zjm2 = 0;
            for (int qb = 1; qb <= 8; qb++) {
              u2zjp2 -= bope(kb, qb) * u(2, i, j + 2, nk - qb + 1);
              u2zjp1 -= bope(kb, qb) * u(2, i, j + 1, nk - qb + 1);
              u2zjm1 -= bope(kb, qb) * u(2, i, j - 1, nk - qb + 1);
              u2zjm2 -= bope(kb, qb) * u(2, i, j - 2, nk - qb + 1);
            }
            float_sw4 mu2zy =
                i12 *
                (-mu(i, j + 2, k) * u2zjp2 + 8 * mu(i, j + 1, k) * u2zjp1 -
                 8 * mu(i, j - 1, k) * u2zjm1 + mu(i, j - 2, k) * u2zjm2);
            r3 = r3 + stry(j) * mu2zy;

            /*   (la*u_x)_z: NOT CENTERED */
            float_sw4 lau1xz = 0;
            for (int qb = 1; qb <= 8; qb++)
              lau1xz -= bope(kb, qb) * (la(i, j, nk - qb + 1) * i12 *
                                        (-u(1, i + 2, j, nk - qb + 1) +
                                         8 * u(1, i + 1, j, nk - qb + 1) -
                                         8 * u(1, i - 1, j, nk - qb + 1) +
                                         u(1, i - 2, j, nk - qb + 1)));
            r3 = r3 + strx(i) * lau1xz;

            /* (la*v_y)_z: NOT CENTERED */
            float_sw4 lau2yz = 0;
            for (int qb = 1; qb <= 8; qb++) {
              lau2yz -= bope(kb, qb) * (la(i, j, nk - qb + 1) * i12 *
                                        (-u(2, i, j + 2, nk - qb + 1) +
                                         8 * u(2, i, j + 1, nk - qb + 1) -
                                         8 * u(2, i, j - 1, nk - qb + 1) +
                                         u(2, i, j - 2, nk - qb + 1)));
            }
            r3 = r3 + stry(j) * lau2yz;

            lu(1, i, j, k) = a1 * lu(1, i, j, k) + cof * r1;
            lu(2, i, j, k) = a1 * lu(2, i, j, k) + cof * r2;
            lu(3, i, j, k) = a1 * lu(3, i, j, k) + cof * r3;
          });  // End of rhs4th3fort_ci LOOP 3
      SYNC_STREAM;
    }
  }
#undef mu
#undef la
#undef u
#undef lu
#undef strx
#undef stry
#undef strz
#undef acof
#undef bope
#undef ghcof
}

//-----------------------------------------------------------------------
void rhs4th3fortsgstr_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nk,
    int* __restrict__ onesided, float_sw4* __restrict__ a_acof,
    float_sw4* __restrict__ a_bope, float_sw4* __restrict__ a_ghcof,
    float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, float_sw4 h,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
    float_sw4* __restrict__ a_strz, char op) {
  SW4_MARK_FUNCTION;
  // using XRHS_POL =
  //   RAJA::KernelPolicy<
  //   RAJA::statement::CudaKernel<
  //     RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>,
  // 			    RAJA::statement::For<1,
  // RAJA::cuda_threadblock_exec<4>,
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
  // RAJA::statement::Lambda<0> >>>>>;
#ifdef ENABLE_CUDA
#if SW4_RAJA_VERSION == 6
  using XRHS_POL2 =
      RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::For<
          0, RAJA::cuda_block_exec,
          RAJA::statement::For<
              1, RAJA::cuda_block_exec,
              RAJA::statement::For<2, RAJA::cuda_thread_exec,
                                   RAJA::statement::Lambda<0>>>>>>;
#elif SW4_RAJA_VERSION == 7

  using XRHS_POL2 =
      RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::For<
          0, RAJA::cuda_block_x_loop,
          RAJA::statement::For<
              1, RAJA::cuda_block_y_loop,
              RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
                                   RAJA::statement::Lambda<0>>>>>>;

#endif
#else
  using XRHS_POL2 = XRHS_POL;
#endif
  // Xrhs4th3fortsgstr_ci<XRHS_POL2>(ifirst,ilast,jfirst,jlast,kfirst,klast,nk,onesided,a_acof,a_bope,a_ghcof,a_lu,a_u,a_mu,a_lambda,h,a_strx,a_stry,a_strz,op);
  // return;
  // This would work to create multi-dimensional C arrays:
  //   float_sw4** b_ar=(float_sw4*)malloc(ni*nj*sizeof(float_sw4*));
  //   for( int j=0;j<nj;j++)
  //      b_ar[j] = &a_lu[j-1+ni*(1-1)];
  //#define ar(i,j) b_ar[j][i];

  // Direct reuse of fortran code by these macro definitions:
#define mu(i, j, k) a_mu[base + i + ni * (j) + nij * (k)]
#define la(i, j, k) a_lambda[base + i + ni * (j) + nij * (k)]
  // Reversed indexation
#define u(c, i, j, k) a_u[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define lu(c, i, j, k) a_lu[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define strx(i) a_strx[i - ifirst0]
#define stry(j) a_stry[j - jfirst0]
#define strz(k) a_strz[k - kfirst0]
#define acof(i, j, k) a_acof[(i - 1) + 6 * (j - 1) + 48 * (k - 1)]
#define bope(i, j) a_bope[i - 1 + 6 * (j - 1)]
#define ghcof(i) a_ghcof[i - 1]
#define ltu(c, i, j, k) tu[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
  //   const float_sw4 a1   = 0;
  const float_sw4 i6 = 1.0 / 6;
  const float_sw4 i12 = 1.0 / 12;
  const float_sw4 i144 = 1.0 / 144;
  const float_sw4 tf = 0.75;

  const int ni = ilast - ifirst + 1;
  //   const int nj= jlast-jfirst+1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;
  // const int nic  = 3*ni;
  // const int nijc = 3*nij;
  const int ifirst0 = ifirst;
  const int jfirst0 = jfirst;
  const int kfirst0 = kfirst;

  int k1, k2;
  // int i, j, k, m,a1;
  // float_sw4 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4, muz1, muz2, muz3,
  // muz4; float_sw4 r1,r2,r3,mucof, mu1zz, mu2zz, mu3zz; float_sw4 lap2mu,
  // u3zip2, u3zip1, u3zim1, u3zim2, lau3zx, mu3xz, u3zjp2, u3zjp1, u3zjm1,
  // u3zjm2; float_sw4 lau3zy, mu3yz, mu1zx, mu2zy, u1zip2, u1zip1, u1zim1,
  // u1zim2; float_sw4 u2zjp2, u2zjp1, u2zjm1, u2zjm2, lau1xz, lau2yz;

  float_sw4 cof = 1.0 / (h * h);
  int a1;

  if (op == '=')
    a1 = 0;
  else if (op == '+')
    a1 = 1;
  else if (op == '-') {
    a1 = 1;
    cof = -cof;
  }

  k1 = kfirst + 2;
  if (onesided[4] == 1) k1 = 7;
  k2 = klast - 2;
  if (onesided[5] == 1) k2 = nk - 6;

  ASSERT_MANAGED(a_mu);
  ASSERT_MANAGED(a_lambda);
  ASSERT_MANAGED(a_u);
  ASSERT_MANAGED(a_lu);
  ASSERT_MANAGED(a_acof);
  ASSERT_MANAGED(a_bope);
  ASSERT_MANAGED(a_ghcof);
  ASSERT_MANAGED(a_strx);
  ASSERT_MANAGED(a_stry);
  ASSERT_MANAGED(a_strz);

  PREFETCH(a_mu);      // Needed
  PREFETCH(a_lambda);  // Needed
                       //   int nkk=klast-kfirst+1;
#ifdef ENABLE_CUDA
                       // rhs4th3fortsgstr_ciopt(2,ni-2,2,nj-2,2,nkk-2,
  //    			  ni,nj,nkk,
  // 			  a_lu,a_u,
  // 			  a_mu,a_lambda,
  // 			  h,a_strx,a_stry,a_strz,op);

  // if ( (onesided[4] == 1 ) ||  (onesided[5] == 1 )) {
  //   std::cout<<"ERROR:: ONE SIDED !! THIS NEEDS TO BE FIXED\n";
  // } else return;
#endif

  // PREFETCH(a_u);
  // PREFETCH(a_lu);
  // PREFETCH(a_strx);
  // PREFETCH(a_stry);
  // PREFETCH(a_strz);
  // using RHS_POL =
  //   RAJA::KernelPolicy<
  //   RAJA::statement::CudaKernel<
  //     RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>,
  // 			    RAJA::statement::For<1,
  // RAJA::cuda_threadblock_exec<4>,
  // RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
  // RAJA::statement::Lambda<0> >>>>>;

  // using RHS_POL = XRHS_POL;
  using RHS_POL = XRHS_POL_ASYNC;
  {
    // #pragma omp parallel private(k,i,j,mux1,mux2,mux3,mux4,muy1,muy2,muy3,muy4,\
//               r1,r2,r3,mucof,mu1zz,mu2zz,mu3zz,lap2mu,q,u3zip2,u3zip1,\
//               u3zim1,u3zim2,lau3zx,mu3xz,u3zjp2,u3zjp1,u3zjm1,u3zjm2,lau3zy,\
//               mu3yz,mu1zx,u1zip2,u1zip1,u1zim1,u1zim2,\
// 	      u2zjp2,u2zjp1,u2zjm1,u2zjm2,mu2zy,lau1xz,lau2yz,kb,qb,mb,muz1,muz2,muz3,muz4)
    //    {
    // #pragma omp for
    //    for( k= k1; k <= k2 ; k++ )
    //       for( j=jfirst+2; j <= jlast-2 ; j++ )
    // #pragma simd
    // #pragma ivdep
    // 	 for( i=ifirst+2; i <= ilast-2 ; i++ )
    // 	 {
    RAJA::RangeSegment k_range(k1, k2 + 1);
    RAJA::RangeSegment j_range(jfirst + 2, jlast - 1);
    RAJA::RangeSegment i_range(ifirst + 2, ilast - 1);
    SW4_MARK_BEGIN("rhs4th3fortsgstr_ci::LOOP1");
#ifdef ENABLE_CUDA
#define NO_COLLAPSE 1
#endif
#if defined(NO_COLLAPSE)

#ifdef SW4_AUTOTUNE
    RangeAT<384, __LINE__, 1> IJK_AT(ifirst + 2, ilast - 1, jfirst + 2,
                                     jlast - 1, k1, k2 + 1);
    forall3asyncAT<__LINE__>(IJK_AT, [=] RAJA_DEVICE(int i, int j, int k) {
#else
    Range<16> I(ifirst + 2, ilast - 1);
    Range<4> J(jfirst + 2, jlast - 1);
    Range<4> K(k1, k2 + 1);
    forall3async(I, J, K, [=] RAJA_DEVICE(int i, int j, int k) {
#endif
#else
    RAJA::kernel<
        XRHS_POL_ASYNC>(RAJA::make_tuple(k_range, j_range, i_range), [=] RAJA_DEVICE(
                                                                         int k,
                                                                         int j,
                                                                         int i) {
#endif
      float_sw4 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4, muz1, muz2,
          muz3, muz4;
      float_sw4 r1, r2, r3;

      /* from inner_loop_4a, 28x3 = 84 ops */
      mux1 = mu(i - 1, j, k) * strx(i - 1) -
             tf * (mu(i, j, k) * strx(i) + mu(i - 2, j, k) * strx(i - 2));
      mux2 = mu(i - 2, j, k) * strx(i - 2) + mu(i + 1, j, k) * strx(i + 1) +
             3 * (mu(i, j, k) * strx(i) + mu(i - 1, j, k) * strx(i - 1));
      mux3 = mu(i - 1, j, k) * strx(i - 1) + mu(i + 2, j, k) * strx(i + 2) +
             3 * (mu(i + 1, j, k) * strx(i + 1) + mu(i, j, k) * strx(i));
      mux4 = mu(i + 1, j, k) * strx(i + 1) -
             tf * (mu(i, j, k) * strx(i) + mu(i + 2, j, k) * strx(i + 2));

      muy1 = mu(i, j - 1, k) * stry(j - 1) -
             tf * (mu(i, j, k) * stry(j) + mu(i, j - 2, k) * stry(j - 2));
      muy2 = mu(i, j - 2, k) * stry(j - 2) + mu(i, j + 1, k) * stry(j + 1) +
             3 * (mu(i, j, k) * stry(j) + mu(i, j - 1, k) * stry(j - 1));
      muy3 = mu(i, j - 1, k) * stry(j - 1) + mu(i, j + 2, k) * stry(j + 2) +
             3 * (mu(i, j + 1, k) * stry(j + 1) + mu(i, j, k) * stry(j));
      muy4 = mu(i, j + 1, k) * stry(j + 1) -
             tf * (mu(i, j, k) * stry(j) + mu(i, j + 2, k) * stry(j + 2));

      muz1 = mu(i, j, k - 1) * strz(k - 1) -
             tf * (mu(i, j, k) * strz(k) + mu(i, j, k - 2) * strz(k - 2));
      muz2 = mu(i, j, k - 2) * strz(k - 2) + mu(i, j, k + 1) * strz(k + 1) +
             3 * (mu(i, j, k) * strz(k) + mu(i, j, k - 1) * strz(k - 1));
      muz3 = mu(i, j, k - 1) * strz(k - 1) + mu(i, j, k + 2) * strz(k + 2) +
             3 * (mu(i, j, k + 1) * strz(k + 1) + mu(i, j, k) * strz(k));
      muz4 = mu(i, j, k + 1) * strz(k + 1) -
             tf * (mu(i, j, k) * strz(k) + mu(i, j, k + 2) * strz(k + 2));
      /* xx, yy, and zz derivatives:*/
      /* 75 ops */
      r1 =
          i6 *
          (strx(i) *
               ((2 * mux1 + la(i - 1, j, k) * strx(i - 1) -
                 tf * (la(i, j, k) * strx(i) + la(i - 2, j, k) * strx(i - 2))) *
                    (u(1, i - 2, j, k) - u(1, i, j, k)) +
                (2 * mux2 + la(i - 2, j, k) * strx(i - 2) +
                 la(i + 1, j, k) * strx(i + 1) +
                 3 * (la(i, j, k) * strx(i) + la(i - 1, j, k) * strx(i - 1))) *
                    (u(1, i - 1, j, k) - u(1, i, j, k)) +
                (2 * mux3 + la(i - 1, j, k) * strx(i - 1) +
                 la(i + 2, j, k) * strx(i + 2) +
                 3 * (la(i + 1, j, k) * strx(i + 1) + la(i, j, k) * strx(i))) *
                    (u(1, i + 1, j, k) - u(1, i, j, k)) +
                (2 * mux4 + la(i + 1, j, k) * strx(i + 1) -
                 tf * (la(i, j, k) * strx(i) + la(i + 2, j, k) * strx(i + 2))) *
                    (u(1, i + 2, j, k) - u(1, i, j, k))) +
           stry(j) * (muy1 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                      muy2 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                      muy3 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                      muy4 * (u(1, i, j + 2, k) - u(1, i, j, k))) +
           strz(k) * (muz1 * (u(1, i, j, k - 2) - u(1, i, j, k)) +
                      muz2 * (u(1, i, j, k - 1) - u(1, i, j, k)) +
                      muz3 * (u(1, i, j, k + 1) - u(1, i, j, k)) +
                      muz4 * (u(1, i, j, k + 2) - u(1, i, j, k))));

      /* 75 ops */
      r2 =
          i6 *
          (strx(i) * (mux1 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                      mux2 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                      mux3 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                      mux4 * (u(2, i + 2, j, k) - u(2, i, j, k))) +
           stry(j) *
               ((2 * muy1 + la(i, j - 1, k) * stry(j - 1) -
                 tf * (la(i, j, k) * stry(j) + la(i, j - 2, k) * stry(j - 2))) *
                    (u(2, i, j - 2, k) - u(2, i, j, k)) +
                (2 * muy2 + la(i, j - 2, k) * stry(j - 2) +
                 la(i, j + 1, k) * stry(j + 1) +
                 3 * (la(i, j, k) * stry(j) + la(i, j - 1, k) * stry(j - 1))) *
                    (u(2, i, j - 1, k) - u(2, i, j, k)) +
                (2 * muy3 + la(i, j - 1, k) * stry(j - 1) +
                 la(i, j + 2, k) * stry(j + 2) +
                 3 * (la(i, j + 1, k) * stry(j + 1) + la(i, j, k) * stry(j))) *
                    (u(2, i, j + 1, k) - u(2, i, j, k)) +
                (2 * muy4 + la(i, j + 1, k) * stry(j + 1) -
                 tf * (la(i, j, k) * stry(j) + la(i, j + 2, k) * stry(j + 2))) *
                    (u(2, i, j + 2, k) - u(2, i, j, k))) +
           strz(k) * (muz1 * (u(2, i, j, k - 2) - u(2, i, j, k)) +
                      muz2 * (u(2, i, j, k - 1) - u(2, i, j, k)) +
                      muz3 * (u(2, i, j, k + 1) - u(2, i, j, k)) +
                      muz4 * (u(2, i, j, k + 2) - u(2, i, j, k))));

      /* 75 ops */
      r3 =
          i6 *
          (strx(i) * (mux1 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                      mux2 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                      mux3 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                      mux4 * (u(3, i + 2, j, k) - u(3, i, j, k))) +
           stry(j) * (muy1 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                      muy2 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                      muy3 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                      muy4 * (u(3, i, j + 2, k) - u(3, i, j, k))) +
           strz(k) *
               ((2 * muz1 + la(i, j, k - 1) * strz(k - 1) -
                 tf * (la(i, j, k) * strz(k) + la(i, j, k - 2) * strz(k - 2))) *
                    (u(3, i, j, k - 2) - u(3, i, j, k)) +
                (2 * muz2 + la(i, j, k - 2) * strz(k - 2) +
                 la(i, j, k + 1) * strz(k + 1) +
                 3 * (la(i, j, k) * strz(k) + la(i, j, k - 1) * strz(k - 1))) *
                    (u(3, i, j, k - 1) - u(3, i, j, k)) +
                (2 * muz3 + la(i, j, k - 1) * strz(k - 1) +
                 la(i, j, k + 2) * strz(k + 2) +
                 3 * (la(i, j, k + 1) * strz(k + 1) + la(i, j, k) * strz(k))) *
                    (u(3, i, j, k + 1) - u(3, i, j, k)) +
                (2 * muz4 + la(i, j, k + 1) * strz(k + 1) -
                 tf * (la(i, j, k) * strz(k) + la(i, j, k + 2) * strz(k + 2))) *
                    (u(3, i, j, k + 2) - u(3, i, j, k))));

      /* Mixed derivatives: */
      /* 29ops /mixed derivative */
      /* 116 ops for r1 */
      /*   (la*v_y)_x */
      r1 = r1 +
           strx(i) * stry(j) * i144 *
               (la(i - 2, j, k) *
                    (u(2, i - 2, j - 2, k) - u(2, i - 2, j + 2, k) +
                     8 * (-u(2, i - 2, j - 1, k) + u(2, i - 2, j + 1, k))) -
                8 * (la(i - 1, j, k) *
                     (u(2, i - 1, j - 2, k) - u(2, i - 1, j + 2, k) +
                      8 * (-u(2, i - 1, j - 1, k) + u(2, i - 1, j + 1, k)))) +
                8 * (la(i + 1, j, k) *
                     (u(2, i + 1, j - 2, k) - u(2, i + 1, j + 2, k) +
                      8 * (-u(2, i + 1, j - 1, k) + u(2, i + 1, j + 1, k)))) -
                (la(i + 2, j, k) *
                 (u(2, i + 2, j - 2, k) - u(2, i + 2, j + 2, k) +
                  8 * (-u(2, i + 2, j - 1, k) + u(2, i + 2, j + 1, k)))))
           /*   (la*w_z)_x */
           + strx(i) * strz(k) * i144 *
                 (la(i - 2, j, k) *
                      (u(3, i - 2, j, k - 2) - u(3, i - 2, j, k + 2) +
                       8 * (-u(3, i - 2, j, k - 1) + u(3, i - 2, j, k + 1))) -
                  8 * (la(i - 1, j, k) *
                       (u(3, i - 1, j, k - 2) - u(3, i - 1, j, k + 2) +
                        8 * (-u(3, i - 1, j, k - 1) + u(3, i - 1, j, k + 1)))) +
                  8 * (la(i + 1, j, k) *
                       (u(3, i + 1, j, k - 2) - u(3, i + 1, j, k + 2) +
                        8 * (-u(3, i + 1, j, k - 1) + u(3, i + 1, j, k + 1)))) -
                  (la(i + 2, j, k) *
                   (u(3, i + 2, j, k - 2) - u(3, i + 2, j, k + 2) +
                    8 * (-u(3, i + 2, j, k - 1) + u(3, i + 2, j, k + 1)))))
           /*   (mu*v_x)_y */
           + strx(i) * stry(j) * i144 *
                 (mu(i, j - 2, k) *
                      (u(2, i - 2, j - 2, k) - u(2, i + 2, j - 2, k) +
                       8 * (-u(2, i - 1, j - 2, k) + u(2, i + 1, j - 2, k))) -
                  8 * (mu(i, j - 1, k) *
                       (u(2, i - 2, j - 1, k) - u(2, i + 2, j - 1, k) +
                        8 * (-u(2, i - 1, j - 1, k) + u(2, i + 1, j - 1, k)))) +
                  8 * (mu(i, j + 1, k) *
                       (u(2, i - 2, j + 1, k) - u(2, i + 2, j + 1, k) +
                        8 * (-u(2, i - 1, j + 1, k) + u(2, i + 1, j + 1, k)))) -
                  (mu(i, j + 2, k) *
                   (u(2, i - 2, j + 2, k) - u(2, i + 2, j + 2, k) +
                    8 * (-u(2, i - 1, j + 2, k) + u(2, i + 1, j + 2, k)))))
           /*   (mu*w_x)_z */
           + strx(i) * strz(k) * i144 *
                 (mu(i, j, k - 2) *
                      (u(3, i - 2, j, k - 2) - u(3, i + 2, j, k - 2) +
                       8 * (-u(3, i - 1, j, k - 2) + u(3, i + 1, j, k - 2))) -
                  8 * (mu(i, j, k - 1) *
                       (u(3, i - 2, j, k - 1) - u(3, i + 2, j, k - 1) +
                        8 * (-u(3, i - 1, j, k - 1) + u(3, i + 1, j, k - 1)))) +
                  8 * (mu(i, j, k + 1) *
                       (u(3, i - 2, j, k + 1) - u(3, i + 2, j, k + 1) +
                        8 * (-u(3, i - 1, j, k + 1) + u(3, i + 1, j, k + 1)))) -
                  (mu(i, j, k + 2) *
                   (u(3, i - 2, j, k + 2) - u(3, i + 2, j, k + 2) +
                    8 * (-u(3, i - 1, j, k + 2) + u(3, i + 1, j, k + 2)))));

      /* 116 ops for r2 */
      /*   (mu*u_y)_x */
      r2 = r2 +
           strx(i) * stry(j) * i144 *
               (mu(i - 2, j, k) *
                    (u(1, i - 2, j - 2, k) - u(1, i - 2, j + 2, k) +
                     8 * (-u(1, i - 2, j - 1, k) + u(1, i - 2, j + 1, k))) -
                8 * (mu(i - 1, j, k) *
                     (u(1, i - 1, j - 2, k) - u(1, i - 1, j + 2, k) +
                      8 * (-u(1, i - 1, j - 1, k) + u(1, i - 1, j + 1, k)))) +
                8 * (mu(i + 1, j, k) *
                     (u(1, i + 1, j - 2, k) - u(1, i + 1, j + 2, k) +
                      8 * (-u(1, i + 1, j - 1, k) + u(1, i + 1, j + 1, k)))) -
                (mu(i + 2, j, k) *
                 (u(1, i + 2, j - 2, k) - u(1, i + 2, j + 2, k) +
                  8 * (-u(1, i + 2, j - 1, k) + u(1, i + 2, j + 1, k)))))
           /* (la*u_x)_y */
           + strx(i) * stry(j) * i144 *
                 (la(i, j - 2, k) *
                      (u(1, i - 2, j - 2, k) - u(1, i + 2, j - 2, k) +
                       8 * (-u(1, i - 1, j - 2, k) + u(1, i + 1, j - 2, k))) -
                  8 * (la(i, j - 1, k) *
                       (u(1, i - 2, j - 1, k) - u(1, i + 2, j - 1, k) +
                        8 * (-u(1, i - 1, j - 1, k) + u(1, i + 1, j - 1, k)))) +
                  8 * (la(i, j + 1, k) *
                       (u(1, i - 2, j + 1, k) - u(1, i + 2, j + 1, k) +
                        8 * (-u(1, i - 1, j + 1, k) + u(1, i + 1, j + 1, k)))) -
                  (la(i, j + 2, k) *
                   (u(1, i - 2, j + 2, k) - u(1, i + 2, j + 2, k) +
                    8 * (-u(1, i - 1, j + 2, k) + u(1, i + 1, j + 2, k)))))
           /* (la*w_z)_y */
           + stry(j) * strz(k) * i144 *
                 (la(i, j - 2, k) *
                      (u(3, i, j - 2, k - 2) - u(3, i, j - 2, k + 2) +
                       8 * (-u(3, i, j - 2, k - 1) + u(3, i, j - 2, k + 1))) -
                  8 * (la(i, j - 1, k) *
                       (u(3, i, j - 1, k - 2) - u(3, i, j - 1, k + 2) +
                        8 * (-u(3, i, j - 1, k - 1) + u(3, i, j - 1, k + 1)))) +
                  8 * (la(i, j + 1, k) *
                       (u(3, i, j + 1, k - 2) - u(3, i, j + 1, k + 2) +
                        8 * (-u(3, i, j + 1, k - 1) + u(3, i, j + 1, k + 1)))) -
                  (la(i, j + 2, k) *
                   (u(3, i, j + 2, k - 2) - u(3, i, j + 2, k + 2) +
                    8 * (-u(3, i, j + 2, k - 1) + u(3, i, j + 2, k + 1)))))
           /* (mu*w_y)_z */
           + stry(j) * strz(k) * i144 *
                 (mu(i, j, k - 2) *
                      (u(3, i, j - 2, k - 2) - u(3, i, j + 2, k - 2) +
                       8 * (-u(3, i, j - 1, k - 2) + u(3, i, j + 1, k - 2))) -
                  8 * (mu(i, j, k - 1) *
                       (u(3, i, j - 2, k - 1) - u(3, i, j + 2, k - 1) +
                        8 * (-u(3, i, j - 1, k - 1) + u(3, i, j + 1, k - 1)))) +
                  8 * (mu(i, j, k + 1) *
                       (u(3, i, j - 2, k + 1) - u(3, i, j + 2, k + 1) +
                        8 * (-u(3, i, j - 1, k + 1) + u(3, i, j + 1, k + 1)))) -
                  (mu(i, j, k + 2) *
                   (u(3, i, j - 2, k + 2) - u(3, i, j + 2, k + 2) +
                    8 * (-u(3, i, j - 1, k + 2) + u(3, i, j + 1, k + 2)))));
      /* 116 ops for r3 */
      /*  (mu*u_z)_x */
      r3 = r3 +
           strx(i) * strz(k) * i144 *
               (mu(i - 2, j, k) *
                    (u(1, i - 2, j, k - 2) - u(1, i - 2, j, k + 2) +
                     8 * (-u(1, i - 2, j, k - 1) + u(1, i - 2, j, k + 1))) -
                8 * (mu(i - 1, j, k) *
                     (u(1, i - 1, j, k - 2) - u(1, i - 1, j, k + 2) +
                      8 * (-u(1, i - 1, j, k - 1) + u(1, i - 1, j, k + 1)))) +
                8 * (mu(i + 1, j, k) *
                     (u(1, i + 1, j, k - 2) - u(1, i + 1, j, k + 2) +
                      8 * (-u(1, i + 1, j, k - 1) + u(1, i + 1, j, k + 1)))) -
                (mu(i + 2, j, k) *
                 (u(1, i + 2, j, k - 2) - u(1, i + 2, j, k + 2) +
                  8 * (-u(1, i + 2, j, k - 1) + u(1, i + 2, j, k + 1)))))
           /* (mu*v_z)_y */
           + stry(j) * strz(k) * i144 *
                 (mu(i, j - 2, k) *
                      (u(2, i, j - 2, k - 2) - u(2, i, j - 2, k + 2) +
                       8 * (-u(2, i, j - 2, k - 1) + u(2, i, j - 2, k + 1))) -
                  8 * (mu(i, j - 1, k) *
                       (u(2, i, j - 1, k - 2) - u(2, i, j - 1, k + 2) +
                        8 * (-u(2, i, j - 1, k - 1) + u(2, i, j - 1, k + 1)))) +
                  8 * (mu(i, j + 1, k) *
                       (u(2, i, j + 1, k - 2) - u(2, i, j + 1, k + 2) +
                        8 * (-u(2, i, j + 1, k - 1) + u(2, i, j + 1, k + 1)))) -
                  (mu(i, j + 2, k) *
                   (u(2, i, j + 2, k - 2) - u(2, i, j + 2, k + 2) +
                    8 * (-u(2, i, j + 2, k - 1) + u(2, i, j + 2, k + 1)))))
           /*   (la*u_x)_z */
           + strx(i) * strz(k) * i144 *
                 (la(i, j, k - 2) *
                      (u(1, i - 2, j, k - 2) - u(1, i + 2, j, k - 2) +
                       8 * (-u(1, i - 1, j, k - 2) + u(1, i + 1, j, k - 2))) -
                  8 * (la(i, j, k - 1) *
                       (u(1, i - 2, j, k - 1) - u(1, i + 2, j, k - 1) +
                        8 * (-u(1, i - 1, j, k - 1) + u(1, i + 1, j, k - 1)))) +
                  8 * (la(i, j, k + 1) *
                       (u(1, i - 2, j, k + 1) - u(1, i + 2, j, k + 1) +
                        8 * (-u(1, i - 1, j, k + 1) + u(1, i + 1, j, k + 1)))) -
                  (la(i, j, k + 2) *
                   (u(1, i - 2, j, k + 2) - u(1, i + 2, j, k + 2) +
                    8 * (-u(1, i - 1, j, k + 2) + u(1, i + 1, j, k + 2)))))
           /* (la*v_y)_z */
           + stry(j) * strz(k) * i144 *
                 (la(i, j, k - 2) *
                      (u(2, i, j - 2, k - 2) - u(2, i, j + 2, k - 2) +
                       8 * (-u(2, i, j - 1, k - 2) + u(2, i, j + 1, k - 2))) -
                  8 * (la(i, j, k - 1) *
                       (u(2, i, j - 2, k - 1) - u(2, i, j + 2, k - 1) +
                        8 * (-u(2, i, j - 1, k - 1) + u(2, i, j + 1, k - 1)))) +
                  8 * (la(i, j, k + 1) *
                       (u(2, i, j - 2, k + 1) - u(2, i, j + 2, k + 1) +
                        8 * (-u(2, i, j - 1, k + 1) + u(2, i, j + 1, k + 1)))) -
                  (la(i, j, k + 2) *
                   (u(2, i, j - 2, k + 2) - u(2, i, j + 2, k + 2) +
                    8 * (-u(2, i, j - 1, k + 2) + u(2, i, j + 1, k + 2)))));

      /* 9 ops */
      lu(1, i, j, k) = a1 * lu(1, i, j, k) + cof * r1;
      lu(2, i, j, k) = a1 * lu(2, i, j, k) + cof * r2;
      lu(3, i, j, k) = a1 * lu(3, i, j, k) + cof * r3;
    });  // END OF rhs4th3fortsgstr_ci LOOP 1
    // SYNC_STREAM;
    // printf("END LOOP1\n");
    SW4_MARK_END("rhs4th3fortsgstr_ci::LOOP1");
    if (onesided[4] == 1) {
      RAJA::RangeSegment k_range(1, 6 + 1);
      RAJA::RangeSegment j_range(jfirst + 2, jlast - 1);
      RAJA::RangeSegment i_range(ifirst + 2, ilast - 1);

      SW4_MARK_BEGIN("rhs4th3fortsgstr_ci::LOOP2");
      // printf("START LOOP2 \n");

      // using LOCAL_POL_ASYNC_OLDE = RAJA::KernelPolicy<
      //   RAJA::statement::CudaKernelAsync<
      //     RAJA::statement::Tile<0, RAJA::statement::tile_fixed<6>,
      //     RAJA::cuda_block_z_loop,
      // RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>,
      // RAJA::cuda_block_y_loop, 		      RAJA::statement::Tile<2,
      // RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
      //   RAJA::statement::For<0, RAJA::cuda_thread_z_loop,
      //     RAJA::statement::For<1, RAJA::cuda_thread_y_loop,
      // 			 RAJA::statement::For<2,
      // RAJA::cuda_thread_x_loop,
      // RAJA::statement::Lambda<0> >>>>>>>>;

#ifdef ENABLE_CUDA
      using LOCAL_POL_ASYNC =
          RAJA::KernelPolicy<RAJA::statement::CudaKernelFixedAsync<
              256,
              RAJA::statement::Tile<
                  0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
                  RAJA::statement::Tile<
                      1, RAJA::statement::tile_fixed<4>,
                      RAJA::cuda_block_y_loop,
                      RAJA::statement::Tile<
                          2, RAJA::statement::tile_fixed<16>,
                          RAJA::cuda_block_x_loop,
                          RAJA::statement::For<
                              0, RAJA::cuda_thread_z_direct,
                              RAJA::statement::For<
                                  1, RAJA::cuda_thread_y_direct,
                                  RAJA::statement::For<
                                      2, RAJA::cuda_thread_x_direct,
                                      RAJA::statement::Lambda<0>>>>>>>>>;
#else
      using LOCAL_POL_ASYNC = XRHS_POL;
#endif

      RAJA::kernel<LOCAL_POL_ASYNC>(
          RAJA::make_tuple(k_range, j_range, i_range),
          [=] RAJA_DEVICE(int k, int j, int i) {
            float_sw4 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4;
            float_sw4 r1, r2, r3;
            // #pragma omp for
            // 	 for( k=1 ; k<= 6 ; k++ )
            // /* the centered stencil can be used in the x- and y-directions */
            // 	    for( j=jfirst+2; j<=jlast-2; j++ )
            // #pragma simd
            // #pragma ivdep
            // 	       for( i=ifirst+2; i<=ilast-2; i++ )
            // 	       {
            /* from inner_loop_4a */
            mux1 = mu(i - 1, j, k) * strx(i - 1) -
                   tf * (mu(i, j, k) * strx(i) + mu(i - 2, j, k) * strx(i - 2));
            mux2 = mu(i - 2, j, k) * strx(i - 2) +
                   mu(i + 1, j, k) * strx(i + 1) +
                   3 * (mu(i, j, k) * strx(i) + mu(i - 1, j, k) * strx(i - 1));
            mux3 = mu(i - 1, j, k) * strx(i - 1) +
                   mu(i + 2, j, k) * strx(i + 2) +
                   3 * (mu(i + 1, j, k) * strx(i + 1) + mu(i, j, k) * strx(i));
            mux4 = mu(i + 1, j, k) * strx(i + 1) -
                   tf * (mu(i, j, k) * strx(i) + mu(i + 2, j, k) * strx(i + 2));

            muy1 = mu(i, j - 1, k) * stry(j - 1) -
                   tf * (mu(i, j, k) * stry(j) + mu(i, j - 2, k) * stry(j - 2));
            muy2 = mu(i, j - 2, k) * stry(j - 2) +
                   mu(i, j + 1, k) * stry(j + 1) +
                   3 * (mu(i, j, k) * stry(j) + mu(i, j - 1, k) * stry(j - 1));
            muy3 = mu(i, j - 1, k) * stry(j - 1) +
                   mu(i, j + 2, k) * stry(j + 2) +
                   3 * (mu(i, j + 1, k) * stry(j + 1) + mu(i, j, k) * stry(j));
            muy4 = mu(i, j + 1, k) * stry(j + 1) -
                   tf * (mu(i, j, k) * stry(j) + mu(i, j + 2, k) * stry(j + 2));

            r1 = i6 * (strx(i) * ((2 * mux1 + la(i - 1, j, k) * strx(i - 1) -
                                   tf * (la(i, j, k) * strx(i) +
                                         la(i - 2, j, k) * strx(i - 2))) *
                                      (u(1, i - 2, j, k) - u(1, i, j, k)) +
                                  (2 * mux2 + la(i - 2, j, k) * strx(i - 2) +
                                   la(i + 1, j, k) * strx(i + 1) +
                                   3 * (la(i, j, k) * strx(i) +
                                        la(i - 1, j, k) * strx(i - 1))) *
                                      (u(1, i - 1, j, k) - u(1, i, j, k)) +
                                  (2 * mux3 + la(i - 1, j, k) * strx(i - 1) +
                                   la(i + 2, j, k) * strx(i + 2) +
                                   3 * (la(i + 1, j, k) * strx(i + 1) +
                                        la(i, j, k) * strx(i))) *
                                      (u(1, i + 1, j, k) - u(1, i, j, k)) +
                                  (2 * mux4 + la(i + 1, j, k) * strx(i + 1) -
                                   tf * (la(i, j, k) * strx(i) +
                                         la(i + 2, j, k) * strx(i + 2))) *
                                      (u(1, i + 2, j, k) - u(1, i, j, k))) +
                       stry(j) * (+muy1 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                                  muy2 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                                  muy3 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                                  muy4 * (u(1, i, j + 2, k) - u(1, i, j, k))));

            /* (mu*uz)_z can not be centered */
            /* second derivative (mu*u_z)_z at grid point z_k */
            /* averaging the coefficient, */
            /* leave out the z-supergrid stretching strz, since it will */
            /* never be used together with the sbp-boundary operator */
            float_sw4 mu1zz = 0;
            float_sw4 mu2zz = 0;
            float_sw4 mu3zz = 0;
            for (int q = 1; q <= 8; q++) {
              //		     lap2mu= 0;
              //		     mucof = 0;
              //		     for( m=1 ; m<=8; m++ )
              //		     {
              //			mucof  += acof(k,q,m)*mu(i,j,m);
              //			lap2mu +=
              // acof(k,q,m)*(la(i,j,m)+2*mu(i,j,m));
              //		     }
              float_sw4 lap2mu =
                  acof(k, q, 1) * (la(i, j, 1) + 2 * mu(i, j, 1)) +
                  acof(k, q, 2) * (la(i, j, 2) + 2 * mu(i, j, 2)) +
                  acof(k, q, 3) * (la(i, j, 3) + 2 * mu(i, j, 3)) +
                  acof(k, q, 4) * (la(i, j, 4) + 2 * mu(i, j, 4)) +
                  acof(k, q, 5) * (la(i, j, 5) + 2 * mu(i, j, 5)) +
                  acof(k, q, 6) * (la(i, j, 6) + 2 * mu(i, j, 6)) +
                  acof(k, q, 7) * (la(i, j, 7) + 2 * mu(i, j, 7)) +
                  acof(k, q, 8) * (la(i, j, 8) + 2 * mu(i, j, 8));
              float_sw4 mucof =
                  acof(k, q, 1) * mu(i, j, 1) + acof(k, q, 2) * mu(i, j, 2) +
                  acof(k, q, 3) * mu(i, j, 3) + acof(k, q, 4) * mu(i, j, 4) +
                  acof(k, q, 5) * mu(i, j, 5) + acof(k, q, 6) * mu(i, j, 6) +
                  acof(k, q, 7) * mu(i, j, 7) + acof(k, q, 8) * mu(i, j, 8);
              mu1zz += mucof * u(1, i, j, q);
              mu2zz += mucof * u(2, i, j, q);
              mu3zz += lap2mu * u(3, i, j, q);
            }

            /* ghost point only influences the first point (k=1) because
             * ghcof(k)=0 for k>=2*/
            r1 = r1 + (mu1zz + ghcof(k) * mu(i, j, 1) * u(1, i, j, 0));

            r2 = i6 * (strx(i) * (mux1 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                                  mux2 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                                  mux3 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                                  mux4 * (u(2, i + 2, j, k) - u(2, i, j, k))) +
                       stry(j) * ((2 * muy1 + la(i, j - 1, k) * stry(j - 1) -
                                   tf * (la(i, j, k) * stry(j) +
                                         la(i, j - 2, k) * stry(j - 2))) *
                                      (u(2, i, j - 2, k) - u(2, i, j, k)) +
                                  (2 * muy2 + la(i, j - 2, k) * stry(j - 2) +
                                   la(i, j + 1, k) * stry(j + 1) +
                                   3 * (la(i, j, k) * stry(j) +
                                        la(i, j - 1, k) * stry(j - 1))) *
                                      (u(2, i, j - 1, k) - u(2, i, j, k)) +
                                  (2 * muy3 + la(i, j - 1, k) * stry(j - 1) +
                                   la(i, j + 2, k) * stry(j + 2) +
                                   3 * (la(i, j + 1, k) * stry(j + 1) +
                                        la(i, j, k) * stry(j))) *
                                      (u(2, i, j + 1, k) - u(2, i, j, k)) +
                                  (2 * muy4 + la(i, j + 1, k) * stry(j + 1) -
                                   tf * (la(i, j, k) * stry(j) +
                                         la(i, j + 2, k) * stry(j + 2))) *
                                      (u(2, i, j + 2, k) - u(2, i, j, k))));

            /* ghost point only influences the first point (k=1) because
             * ghcof(k)=0 for k>=2 */
            r2 = r2 + (mu2zz + ghcof(k) * mu(i, j, 1) * u(2, i, j, 0));

            r3 = i6 * (strx(i) * (mux1 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                                  mux2 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                                  mux3 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                                  mux4 * (u(3, i + 2, j, k) - u(3, i, j, k))) +
                       stry(j) * (muy1 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                                  muy2 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                                  muy3 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                                  muy4 * (u(3, i, j + 2, k) - u(3, i, j, k))));
            /* ghost point only influences the first point (k=1) because
             * ghcof(k)=0 for k>=2 */
            r3 = r3 + (mu3zz + ghcof(k) * (la(i, j, 1) + 2 * mu(i, j, 1)) *
                                   u(3, i, j, 0));

            /* cross-terms in first component of rhs */
            /*   (la*v_y)_x */
            r1 = r1 +
                 strx(i) * stry(j) *
                     (i144 *
                          (la(i - 2, j, k) *
                               (u(2, i - 2, j - 2, k) - u(2, i - 2, j + 2, k) +
                                8 * (-u(2, i - 2, j - 1, k) +
                                     u(2, i - 2, j + 1, k))) -
                           8 * (la(i - 1, j, k) *
                                (u(2, i - 1, j - 2, k) - u(2, i - 1, j + 2, k) +
                                 8 * (-u(2, i - 1, j - 1, k) +
                                      u(2, i - 1, j + 1, k)))) +
                           8 * (la(i + 1, j, k) *
                                (u(2, i + 1, j - 2, k) - u(2, i + 1, j + 2, k) +
                                 8 * (-u(2, i + 1, j - 1, k) +
                                      u(2, i + 1, j + 1, k)))) -
                           (la(i + 2, j, k) *
                            (u(2, i + 2, j - 2, k) - u(2, i + 2, j + 2, k) +
                             8 * (-u(2, i + 2, j - 1, k) +
                                  u(2, i + 2, j + 1, k)))))
                      /*   (mu*v_x)_y */
                      +
                      i144 *
                          (mu(i, j - 2, k) *
                               (u(2, i - 2, j - 2, k) - u(2, i + 2, j - 2, k) +
                                8 * (-u(2, i - 1, j - 2, k) +
                                     u(2, i + 1, j - 2, k))) -
                           8 * (mu(i, j - 1, k) *
                                (u(2, i - 2, j - 1, k) - u(2, i + 2, j - 1, k) +
                                 8 * (-u(2, i - 1, j - 1, k) +
                                      u(2, i + 1, j - 1, k)))) +
                           8 * (mu(i, j + 1, k) *
                                (u(2, i - 2, j + 1, k) - u(2, i + 2, j + 1, k) +
                                 8 * (-u(2, i - 1, j + 1, k) +
                                      u(2, i + 1, j + 1, k)))) -
                           (mu(i, j + 2, k) *
                            (u(2, i - 2, j + 2, k) - u(2, i + 2, j + 2, k) +
                             8 * (-u(2, i - 1, j + 2, k) +
                                  u(2, i + 1, j + 2, k))))));
            /*   (la*w_z)_x: NOT CENTERED */
            float_sw4 u3zip2 = 0;
            float_sw4 u3zip1 = 0;
            float_sw4 u3zim1 = 0;
            float_sw4 u3zim2 = 0;
            for (int q = 1; q <= 8; q++) {
              u3zip2 += bope(k, q) * u(3, i + 2, j, q);
              u3zip1 += bope(k, q) * u(3, i + 1, j, q);
              u3zim1 += bope(k, q) * u(3, i - 1, j, q);
              u3zim2 += bope(k, q) * u(3, i - 2, j, q);
            }
            float_sw4 lau3zx =
                i12 *
                (-la(i + 2, j, k) * u3zip2 + 8 * la(i + 1, j, k) * u3zip1 -
                 8 * la(i - 1, j, k) * u3zim1 + la(i - 2, j, k) * u3zim2);
            r1 = r1 + strx(i) * lau3zx;
            /*   (mu*w_x)_z: NOT CENTERED */
            float_sw4 mu3xz = 0;
            for (int q = 1; q <= 8; q++)
              mu3xz +=
                  bope(k, q) * (mu(i, j, q) * i12 *
                                (-u(3, i + 2, j, q) + 8 * u(3, i + 1, j, q) -
                                 8 * u(3, i - 1, j, q) + u(3, i - 2, j, q)));
            r1 = r1 + strx(i) * mu3xz;

            /* cross-terms in second component of rhs */
            /*   (mu*u_y)_x */
            r2 = r2 +
                 strx(i) * stry(j) *
                     (i144 *
                          (mu(i - 2, j, k) *
                               (u(1, i - 2, j - 2, k) - u(1, i - 2, j + 2, k) +
                                8 * (-u(1, i - 2, j - 1, k) +
                                     u(1, i - 2, j + 1, k))) -
                           8 * (mu(i - 1, j, k) *
                                (u(1, i - 1, j - 2, k) - u(1, i - 1, j + 2, k) +
                                 8 * (-u(1, i - 1, j - 1, k) +
                                      u(1, i - 1, j + 1, k)))) +
                           8 * (mu(i + 1, j, k) *
                                (u(1, i + 1, j - 2, k) - u(1, i + 1, j + 2, k) +
                                 8 * (-u(1, i + 1, j - 1, k) +
                                      u(1, i + 1, j + 1, k)))) -
                           (mu(i + 2, j, k) *
                            (u(1, i + 2, j - 2, k) - u(1, i + 2, j + 2, k) +
                             8 * (-u(1, i + 2, j - 1, k) +
                                  u(1, i + 2, j + 1, k)))))
                      /* (la*u_x)_y  */
                      +
                      i144 *
                          (la(i, j - 2, k) *
                               (u(1, i - 2, j - 2, k) - u(1, i + 2, j - 2, k) +
                                8 * (-u(1, i - 1, j - 2, k) +
                                     u(1, i + 1, j - 2, k))) -
                           8 * (la(i, j - 1, k) *
                                (u(1, i - 2, j - 1, k) - u(1, i + 2, j - 1, k) +
                                 8 * (-u(1, i - 1, j - 1, k) +
                                      u(1, i + 1, j - 1, k)))) +
                           8 * (la(i, j + 1, k) *
                                (u(1, i - 2, j + 1, k) - u(1, i + 2, j + 1, k) +
                                 8 * (-u(1, i - 1, j + 1, k) +
                                      u(1, i + 1, j + 1, k)))) -
                           (la(i, j + 2, k) *
                            (u(1, i - 2, j + 2, k) - u(1, i + 2, j + 2, k) +
                             8 * (-u(1, i - 1, j + 2, k) +
                                  u(1, i + 1, j + 2, k))))));
            /* (la*w_z)_y : NOT CENTERED */
            float_sw4 u3zjp2 = 0;
            float_sw4 u3zjp1 = 0;
            float_sw4 u3zjm1 = 0;
            float_sw4 u3zjm2 = 0;
            for (int q = 1; q <= 8; q++) {
              u3zjp2 += bope(k, q) * u(3, i, j + 2, q);
              u3zjp1 += bope(k, q) * u(3, i, j + 1, q);
              u3zjm1 += bope(k, q) * u(3, i, j - 1, q);
              u3zjm2 += bope(k, q) * u(3, i, j - 2, q);
            }
            float_sw4 lau3zy =
                i12 *
                (-la(i, j + 2, k) * u3zjp2 + 8 * la(i, j + 1, k) * u3zjp1 -
                 8 * la(i, j - 1, k) * u3zjm1 + la(i, j - 2, k) * u3zjm2);

            r2 = r2 + stry(j) * lau3zy;

            /* (mu*w_y)_z: NOT CENTERED */
            float_sw4 mu3yz = 0;
            for (int q = 1; q <= 8; q++)
              mu3yz +=
                  bope(k, q) * (mu(i, j, q) * i12 *
                                (-u(3, i, j + 2, q) + 8 * u(3, i, j + 1, q) -
                                 8 * u(3, i, j - 1, q) + u(3, i, j - 2, q)));

            r2 = r2 + stry(j) * mu3yz;

            /* No centered cross terms in r3 */
            /*  (mu*u_z)_x: NOT CENTERED */
            float_sw4 u1zip2 = 0;
            float_sw4 u1zip1 = 0;
            float_sw4 u1zim1 = 0;
            float_sw4 u1zim2 = 0;
            for (int q = 1; q <= 8; q++) {
              u1zip2 += bope(k, q) * u(1, i + 2, j, q);
              u1zip1 += bope(k, q) * u(1, i + 1, j, q);
              u1zim1 += bope(k, q) * u(1, i - 1, j, q);
              u1zim2 += bope(k, q) * u(1, i - 2, j, q);
            }
            float_sw4 mu1zx =
                i12 *
                (-mu(i + 2, j, k) * u1zip2 + 8 * mu(i + 1, j, k) * u1zip1 -
                 8 * mu(i - 1, j, k) * u1zim1 + mu(i - 2, j, k) * u1zim2);
            r3 = r3 + strx(i) * mu1zx;

            /* (mu*v_z)_y: NOT CENTERED */
            float_sw4 u2zjp2 = 0;
            float_sw4 u2zjp1 = 0;
            float_sw4 u2zjm1 = 0;
            float_sw4 u2zjm2 = 0;
            for (int q = 1; q <= 8; q++) {
              u2zjp2 += bope(k, q) * u(2, i, j + 2, q);
              u2zjp1 += bope(k, q) * u(2, i, j + 1, q);
              u2zjm1 += bope(k, q) * u(2, i, j - 1, q);
              u2zjm2 += bope(k, q) * u(2, i, j - 2, q);
            }
            float_sw4 mu2zy =
                i12 *
                (-mu(i, j + 2, k) * u2zjp2 + 8 * mu(i, j + 1, k) * u2zjp1 -
                 8 * mu(i, j - 1, k) * u2zjm1 + mu(i, j - 2, k) * u2zjm2);
            r3 = r3 + stry(j) * mu2zy;

            /*   (la*u_x)_z: NOT CENTERED */
            float_sw4 lau1xz = 0;
            for (int q = 1; q <= 8; q++)
              lau1xz +=
                  bope(k, q) * (la(i, j, q) * i12 *
                                (-u(1, i + 2, j, q) + 8 * u(1, i + 1, j, q) -
                                 8 * u(1, i - 1, j, q) + u(1, i - 2, j, q)));
            r3 = r3 + strx(i) * lau1xz;

            /* (la*v_y)_z: NOT CENTERED */
            float_sw4 lau2yz = 0;
            for (int q = 1; q <= 8; q++)
              lau2yz +=
                  bope(k, q) * (la(i, j, q) * i12 *
                                (-u(2, i, j + 2, q) + 8 * u(2, i, j + 1, q) -
                                 8 * u(2, i, j - 1, q) + u(2, i, j - 2, q)));
            r3 = r3 + stry(j) * lau2yz;

            lu(1, i, j, k) = a1 * lu(1, i, j, k) + cof * r1;
            lu(2, i, j, k) = a1 * lu(2, i, j, k) + cof * r2;
            lu(3, i, j, k) = a1 * lu(3, i, j, k) + cof * r3;
          });  // End of rhs4th3fortsgstr_ci LOOP 2
      // SYNC_STREAM;
      SW4_MARK_END("rhs4th3fortsgstr_ci::LOOP2");
      // printf("END LOOP2 \n");
    }
    if (onesided[5] == 1) {
      // printf("START LOOP3 \n");
      RAJA::RangeSegment k_range(nk - 5, nk + 1);
      RAJA::RangeSegment j_range(jfirst + 2, jlast - 1);
      RAJA::RangeSegment i_range(ifirst + 2, ilast - 1);

      SW4_MARK_BEGIN("rhs4th3fortsgstr_ci::LOOP3");
#ifdef ENABLE_CUDA
      using LOCAL_POL_ASYNC1 = RAJA::KernelPolicy<
          // RAJA::statement::CudaKernelExt<RAJA::cuda_explicit_launch<true, 0,
          // 256>,
          RAJA::statement::CudaKernelFixedAsync<
              256,
              RAJA::statement::Tile<
                  0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_z_loop,
                  RAJA::statement::Tile<
                      1, RAJA::statement::tile_fixed<4>,
                      RAJA::cuda_block_y_loop,
                      RAJA::statement::Tile<
                          2, RAJA::statement::tile_fixed<16>,
                          RAJA::cuda_block_x_loop,
                          RAJA::statement::For<
                              0, RAJA::cuda_thread_z_direct,
                              RAJA::statement::For<
                                  1, RAJA::cuda_thread_y_direct,
                                  RAJA::statement::For<
                                      2, RAJA::cuda_thread_x_direct,
                                      RAJA::statement::Lambda<0>>>>>>>>>;

      // 	using LOCAL_POL_ASYNC2 =
      // RAJA::KernelPolicy<
      // 	  RAJA::statement::CudaKernelFixedAsync<256,
      //   RAJA::statement::For<0, RAJA::cuda_block_y_loop,
      // 			    RAJA::statement::For<1,
      // RAJA::cuda_block_x_loop,
      // RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
      // RAJA::statement::Lambda<0> >>>>>;
#else
      using LOCAL_POL_ASYNC1 = XRHS_POL;
#endif

      RAJA::kernel<LOCAL_POL_ASYNC1>(
          RAJA::make_tuple(k_range, j_range, i_range),
          [=] RAJA_DEVICE(int k, int j, int i) {
            float_sw4 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4;
            float_sw4 r1, r2, r3;
            // #pragma omp for
            // 	 for(  k = nk-5 ; k <= nk ; k++ )
            // 	    for(  j=jfirst+2; j<=jlast-2; j++ )
            // #pragma simd
            // #pragma ivdep
            // 	       for(  i=ifirst+2; i<=ilast-2; i++ )
            // 	       {
            /* from inner_loop_4a */
            mux1 = mu(i - 1, j, k) * strx(i - 1) -
                   tf * (mu(i, j, k) * strx(i) + mu(i - 2, j, k) * strx(i - 2));
            mux2 = mu(i - 2, j, k) * strx(i - 2) +
                   mu(i + 1, j, k) * strx(i + 1) +
                   3 * (mu(i, j, k) * strx(i) + mu(i - 1, j, k) * strx(i - 1));
            mux3 = mu(i - 1, j, k) * strx(i - 1) +
                   mu(i + 2, j, k) * strx(i + 2) +
                   3 * (mu(i + 1, j, k) * strx(i + 1) + mu(i, j, k) * strx(i));
            mux4 = mu(i + 1, j, k) * strx(i + 1) -
                   tf * (mu(i, j, k) * strx(i) + mu(i + 2, j, k) * strx(i + 2));

            muy1 = mu(i, j - 1, k) * stry(j - 1) -
                   tf * (mu(i, j, k) * stry(j) + mu(i, j - 2, k) * stry(j - 2));
            muy2 = mu(i, j - 2, k) * stry(j - 2) +
                   mu(i, j + 1, k) * stry(j + 1) +
                   3 * (mu(i, j, k) * stry(j) + mu(i, j - 1, k) * stry(j - 1));
            muy3 = mu(i, j - 1, k) * stry(j - 1) +
                   mu(i, j + 2, k) * stry(j + 2) +
                   3 * (mu(i, j + 1, k) * stry(j + 1) + mu(i, j, k) * stry(j));
            muy4 = mu(i, j + 1, k) * stry(j + 1) -
                   tf * (mu(i, j, k) * stry(j) + mu(i, j + 2, k) * stry(j + 2));

            /* xx, yy, and zz derivatives: */
            /* note that we could have introduced intermediate variables for the
             * average of lambda  */
            /* in the same way as we did for mu */
            r1 = i6 * (strx(i) * ((2 * mux1 + la(i - 1, j, k) * strx(i - 1) -
                                   tf * (la(i, j, k) * strx(i) +
                                         la(i - 2, j, k) * strx(i - 2))) *
                                      (u(1, i - 2, j, k) - u(1, i, j, k)) +
                                  (2 * mux2 + la(i - 2, j, k) * strx(i - 2) +
                                   la(i + 1, j, k) * strx(i + 1) +
                                   3 * (la(i, j, k) * strx(i) +
                                        la(i - 1, j, k) * strx(i - 1))) *
                                      (u(1, i - 1, j, k) - u(1, i, j, k)) +
                                  (2 * mux3 + la(i - 1, j, k) * strx(i - 1) +
                                   la(i + 2, j, k) * strx(i + 2) +
                                   3 * (la(i + 1, j, k) * strx(i + 1) +
                                        la(i, j, k) * strx(i))) *
                                      (u(1, i + 1, j, k) - u(1, i, j, k)) +
                                  (2 * mux4 + la(i + 1, j, k) * strx(i + 1) -
                                   tf * (la(i, j, k) * strx(i) +
                                         la(i + 2, j, k) * strx(i + 2))) *
                                      (u(1, i + 2, j, k) - u(1, i, j, k))) +
                       stry(j) * (+muy1 * (u(1, i, j - 2, k) - u(1, i, j, k)) +
                                  muy2 * (u(1, i, j - 1, k) - u(1, i, j, k)) +
                                  muy3 * (u(1, i, j + 1, k) - u(1, i, j, k)) +
                                  muy4 * (u(1, i, j + 2, k) - u(1, i, j, k))));

            /* all indices ending with 'b' are indices relative to the boundary,
             * going into the domain (1,2,3,...)*/
            int kb = nk - k + 1;
            /* all coefficient arrays (acof, bope, ghcof) should be indexed with
             * these indices */
            /* all solution and material property arrays should be indexed with
             * (i,j,k) */

            /* (mu*uz)_z can not be centered */
            /* second derivative (mu*u_z)_z at grid point z_k */
            /* averaging the coefficient */
            float_sw4 mu1zz = 0;
            float_sw4 mu2zz = 0;
            float_sw4 mu3zz = 0;
            for (int qb = 1; qb <= 8; qb++) {
              float_sw4 mucof = 0;
              float_sw4 lap2mu = 0;
              for (int mb = 1; mb <= 8; mb++) {
                mucof += acof(kb, qb, mb) * mu(i, j, nk - mb + 1);
                lap2mu += acof(kb, qb, mb) *
                          (2 * mu(i, j, nk - mb + 1) + la(i, j, nk - mb + 1));
              }
              mu1zz += mucof * u(1, i, j, nk - qb + 1);
              mu2zz += mucof * u(2, i, j, nk - qb + 1);
              mu3zz += lap2mu * u(3, i, j, nk - qb + 1);
            }
            /* computing the second derivative */
            /* ghost point only influences the first point (k=1) because
             * ghcof(k)=0 for k>=2*/
            r1 = r1 + (mu1zz + ghcof(kb) * mu(i, j, nk) * u(1, i, j, nk + 1));

            r2 = i6 * (strx(i) * (mux1 * (u(2, i - 2, j, k) - u(2, i, j, k)) +
                                  mux2 * (u(2, i - 1, j, k) - u(2, i, j, k)) +
                                  mux3 * (u(2, i + 1, j, k) - u(2, i, j, k)) +
                                  mux4 * (u(2, i + 2, j, k) - u(2, i, j, k))) +
                       stry(j) * ((2 * muy1 + la(i, j - 1, k) * stry(j - 1) -
                                   tf * (la(i, j, k) * stry(j) +
                                         la(i, j - 2, k) * stry(j - 2))) *
                                      (u(2, i, j - 2, k) - u(2, i, j, k)) +
                                  (2 * muy2 + la(i, j - 2, k) * stry(j - 2) +
                                   la(i, j + 1, k) * stry(j + 1) +
                                   3 * (la(i, j, k) * stry(j) +
                                        la(i, j - 1, k) * stry(j - 1))) *
                                      (u(2, i, j - 1, k) - u(2, i, j, k)) +
                                  (2 * muy3 + la(i, j - 1, k) * stry(j - 1) +
                                   la(i, j + 2, k) * stry(j + 2) +
                                   3 * (la(i, j + 1, k) * stry(j + 1) +
                                        la(i, j, k) * stry(j))) *
                                      (u(2, i, j + 1, k) - u(2, i, j, k)) +
                                  (2 * muy4 + la(i, j + 1, k) * stry(j + 1) -
                                   tf * (la(i, j, k) * stry(j) +
                                         la(i, j + 2, k) * stry(j + 2))) *
                                      (u(2, i, j + 2, k) - u(2, i, j, k))));

            /* (mu*vz)_z can not be centered */
            /* second derivative (mu*v_z)_z at grid point z_k */
            /* averaging the coefficient: already done above */
            r2 = r2 + (mu2zz + ghcof(kb) * mu(i, j, nk) * u(2, i, j, nk + 1));

            r3 = i6 * (strx(i) * (mux1 * (u(3, i - 2, j, k) - u(3, i, j, k)) +
                                  mux2 * (u(3, i - 1, j, k) - u(3, i, j, k)) +
                                  mux3 * (u(3, i + 1, j, k) - u(3, i, j, k)) +
                                  mux4 * (u(3, i + 2, j, k) - u(3, i, j, k))) +
                       stry(j) * (muy1 * (u(3, i, j - 2, k) - u(3, i, j, k)) +
                                  muy2 * (u(3, i, j - 1, k) - u(3, i, j, k)) +
                                  muy3 * (u(3, i, j + 1, k) - u(3, i, j, k)) +
                                  muy4 * (u(3, i, j + 2, k) - u(3, i, j, k))));
            r3 = r3 + (mu3zz + ghcof(kb) * (la(i, j, nk) + 2 * mu(i, j, nk)) *
                                   u(3, i, j, nk + 1));

            /* cross-terms in first component of rhs */
            /*   (la*v_y)_x */
            r1 = r1 +
                 strx(i) * stry(j) *
                     (i144 *
                          (la(i - 2, j, k) *
                               (u(2, i - 2, j - 2, k) - u(2, i - 2, j + 2, k) +
                                8 * (-u(2, i - 2, j - 1, k) +
                                     u(2, i - 2, j + 1, k))) -
                           8 * (la(i - 1, j, k) *
                                (u(2, i - 1, j - 2, k) - u(2, i - 1, j + 2, k) +
                                 8 * (-u(2, i - 1, j - 1, k) +
                                      u(2, i - 1, j + 1, k)))) +
                           8 * (la(i + 1, j, k) *
                                (u(2, i + 1, j - 2, k) - u(2, i + 1, j + 2, k) +
                                 8 * (-u(2, i + 1, j - 1, k) +
                                      u(2, i + 1, j + 1, k)))) -
                           (la(i + 2, j, k) *
                            (u(2, i + 2, j - 2, k) - u(2, i + 2, j + 2, k) +
                             8 * (-u(2, i + 2, j - 1, k) +
                                  u(2, i + 2, j + 1, k)))))
                      /*   (mu*v_x)_y */
                      +
                      i144 *
                          (mu(i, j - 2, k) *
                               (u(2, i - 2, j - 2, k) - u(2, i + 2, j - 2, k) +
                                8 * (-u(2, i - 1, j - 2, k) +
                                     u(2, i + 1, j - 2, k))) -
                           8 * (mu(i, j - 1, k) *
                                (u(2, i - 2, j - 1, k) - u(2, i + 2, j - 1, k) +
                                 8 * (-u(2, i - 1, j - 1, k) +
                                      u(2, i + 1, j - 1, k)))) +
                           8 * (mu(i, j + 1, k) *
                                (u(2, i - 2, j + 1, k) - u(2, i + 2, j + 1, k) +
                                 8 * (-u(2, i - 1, j + 1, k) +
                                      u(2, i + 1, j + 1, k)))) -
                           (mu(i, j + 2, k) *
                            (u(2, i - 2, j + 2, k) - u(2, i + 2, j + 2, k) +
                             8 * (-u(2, i - 1, j + 2, k) +
                                  u(2, i + 1, j + 2, k))))));
            /*   (la*w_z)_x: NOT CENTERED */
            float_sw4 u3zip2 = 0;
            float_sw4 u3zip1 = 0;
            float_sw4 u3zim1 = 0;
            float_sw4 u3zim2 = 0;
            for (int qb = 1; qb <= 8; qb++) {
              u3zip2 -= bope(kb, qb) * u(3, i + 2, j, nk - qb + 1);
              u3zip1 -= bope(kb, qb) * u(3, i + 1, j, nk - qb + 1);
              u3zim1 -= bope(kb, qb) * u(3, i - 1, j, nk - qb + 1);
              u3zim2 -= bope(kb, qb) * u(3, i - 2, j, nk - qb + 1);
            }
            float_sw4 lau3zx =
                i12 *
                (-la(i + 2, j, k) * u3zip2 + 8 * la(i + 1, j, k) * u3zip1 -
                 8 * la(i - 1, j, k) * u3zim1 + la(i - 2, j, k) * u3zim2);
            r1 = r1 + strx(i) * lau3zx;

            /*   (mu*w_x)_z: NOT CENTERED */
            float_sw4 mu3xz = 0;
            for (int qb = 1; qb <= 8; qb++)
              mu3xz -= bope(kb, qb) * (mu(i, j, nk - qb + 1) * i12 *
                                       (-u(3, i + 2, j, nk - qb + 1) +
                                        8 * u(3, i + 1, j, nk - qb + 1) -
                                        8 * u(3, i - 1, j, nk - qb + 1) +
                                        u(3, i - 2, j, nk - qb + 1)));

            r1 = r1 + strx(i) * mu3xz;

            /* cross-terms in second component of rhs */
            /*   (mu*u_y)_x */
            r2 = r2 +
                 strx(i) * stry(j) *
                     (i144 *
                          (mu(i - 2, j, k) *
                               (u(1, i - 2, j - 2, k) - u(1, i - 2, j + 2, k) +
                                8 * (-u(1, i - 2, j - 1, k) +
                                     u(1, i - 2, j + 1, k))) -
                           8 * (mu(i - 1, j, k) *
                                (u(1, i - 1, j - 2, k) - u(1, i - 1, j + 2, k) +
                                 8 * (-u(1, i - 1, j - 1, k) +
                                      u(1, i - 1, j + 1, k)))) +
                           8 * (mu(i + 1, j, k) *
                                (u(1, i + 1, j - 2, k) - u(1, i + 1, j + 2, k) +
                                 8 * (-u(1, i + 1, j - 1, k) +
                                      u(1, i + 1, j + 1, k)))) -
                           (mu(i + 2, j, k) *
                            (u(1, i + 2, j - 2, k) - u(1, i + 2, j + 2, k) +
                             8 * (-u(1, i + 2, j - 1, k) +
                                  u(1, i + 2, j + 1, k)))))
                      /* (la*u_x)_y */
                      +
                      i144 *
                          (la(i, j - 2, k) *
                               (u(1, i - 2, j - 2, k) - u(1, i + 2, j - 2, k) +
                                8 * (-u(1, i - 1, j - 2, k) +
                                     u(1, i + 1, j - 2, k))) -
                           8 * (la(i, j - 1, k) *
                                (u(1, i - 2, j - 1, k) - u(1, i + 2, j - 1, k) +
                                 8 * (-u(1, i - 1, j - 1, k) +
                                      u(1, i + 1, j - 1, k)))) +
                           8 * (la(i, j + 1, k) *
                                (u(1, i - 2, j + 1, k) - u(1, i + 2, j + 1, k) +
                                 8 * (-u(1, i - 1, j + 1, k) +
                                      u(1, i + 1, j + 1, k)))) -
                           (la(i, j + 2, k) *
                            (u(1, i - 2, j + 2, k) - u(1, i + 2, j + 2, k) +
                             8 * (-u(1, i - 1, j + 2, k) +
                                  u(1, i + 1, j + 2, k))))));
            /* (la*w_z)_y : NOT CENTERED */
            float_sw4 u3zjp2 = 0;
            float_sw4 u3zjp1 = 0;
            float_sw4 u3zjm1 = 0;
            float_sw4 u3zjm2 = 0;
            for (int qb = 1; qb <= 8; qb++) {
              u3zjp2 -= bope(kb, qb) * u(3, i, j + 2, nk - qb + 1);
              u3zjp1 -= bope(kb, qb) * u(3, i, j + 1, nk - qb + 1);
              u3zjm1 -= bope(kb, qb) * u(3, i, j - 1, nk - qb + 1);
              u3zjm2 -= bope(kb, qb) * u(3, i, j - 2, nk - qb + 1);
            }
            float_sw4 lau3zy =
                i12 *
                (-la(i, j + 2, k) * u3zjp2 + 8 * la(i, j + 1, k) * u3zjp1 -
                 8 * la(i, j - 1, k) * u3zjm1 + la(i, j - 2, k) * u3zjm2);
            r2 = r2 + stry(j) * lau3zy;

            /* (mu*w_y)_z: NOT CENTERED */
            float_sw4 mu3yz = 0;
            for (int qb = 1; qb <= 8; qb++)
              mu3yz -= bope(kb, qb) * (mu(i, j, nk - qb + 1) * i12 *
                                       (-u(3, i, j + 2, nk - qb + 1) +
                                        8 * u(3, i, j + 1, nk - qb + 1) -
                                        8 * u(3, i, j - 1, nk - qb + 1) +
                                        u(3, i, j - 2, nk - qb + 1)));
            r2 = r2 + stry(j) * mu3yz;

            /* No centered cross terms in r3 */
            /*  (mu*u_z)_x: NOT CENTERED */
            float_sw4 u1zip2 = 0;
            float_sw4 u1zip1 = 0;
            float_sw4 u1zim1 = 0;
            float_sw4 u1zim2 = 0;
            for (int qb = 1; qb <= 8; qb++) {
              u1zip2 -= bope(kb, qb) * u(1, i + 2, j, nk - qb + 1);
              u1zip1 -= bope(kb, qb) * u(1, i + 1, j, nk - qb + 1);
              u1zim1 -= bope(kb, qb) * u(1, i - 1, j, nk - qb + 1);
              u1zim2 -= bope(kb, qb) * u(1, i - 2, j, nk - qb + 1);
            }
            float_sw4 mu1zx =
                i12 *
                (-mu(i + 2, j, k) * u1zip2 + 8 * mu(i + 1, j, k) * u1zip1 -
                 8 * mu(i - 1, j, k) * u1zim1 + mu(i - 2, j, k) * u1zim2);
            r3 = r3 + strx(i) * mu1zx;

            /* (mu*v_z)_y: NOT CENTERED */
            float_sw4 u2zjp2 = 0;
            float_sw4 u2zjp1 = 0;
            float_sw4 u2zjm1 = 0;
            float_sw4 u2zjm2 = 0;
            for (int qb = 1; qb <= 8; qb++) {
              u2zjp2 -= bope(kb, qb) * u(2, i, j + 2, nk - qb + 1);
              u2zjp1 -= bope(kb, qb) * u(2, i, j + 1, nk - qb + 1);
              u2zjm1 -= bope(kb, qb) * u(2, i, j - 1, nk - qb + 1);
              u2zjm2 -= bope(kb, qb) * u(2, i, j - 2, nk - qb + 1);
            }
            float_sw4 mu2zy =
                i12 *
                (-mu(i, j + 2, k) * u2zjp2 + 8 * mu(i, j + 1, k) * u2zjp1 -
                 8 * mu(i, j - 1, k) * u2zjm1 + mu(i, j - 2, k) * u2zjm2);
            r3 = r3 + stry(j) * mu2zy;

            /*   (la*u_x)_z: NOT CENTERED */
            float_sw4 lau1xz = 0;
            for (int qb = 1; qb <= 8; qb++)
              lau1xz -= bope(kb, qb) * (la(i, j, nk - qb + 1) * i12 *
                                        (-u(1, i + 2, j, nk - qb + 1) +
                                         8 * u(1, i + 1, j, nk - qb + 1) -
                                         8 * u(1, i - 1, j, nk - qb + 1) +
                                         u(1, i - 2, j, nk - qb + 1)));
            r3 = r3 + strx(i) * lau1xz;

            /* (la*v_y)_z: NOT CENTERED */
            float_sw4 lau2yz = 0;
            for (int qb = 1; qb <= 8; qb++) {
              lau2yz -= bope(kb, qb) * (la(i, j, nk - qb + 1) * i12 *
                                        (-u(2, i, j + 2, nk - qb + 1) +
                                         8 * u(2, i, j + 1, nk - qb + 1) -
                                         8 * u(2, i, j - 1, nk - qb + 1) +
                                         u(2, i, j - 2, nk - qb + 1)));
            }
            r3 = r3 + stry(j) * lau2yz;

            lu(1, i, j, k) = a1 * lu(1, i, j, k) + cof * r1;
            lu(2, i, j, k) = a1 * lu(2, i, j, k) + cof * r2;
            lu(3, i, j, k) = a1 * lu(3, i, j, k) + cof * r3;
          });  // End of rhs4th3fortsgstr_ci LOOP 3
      // SYNC_STREAM;
      SW4_MARK_END("rhs4th3fortsgstr_ci::LOOP3");
      // printf("END LOOP3 \n");
    }
    // SYNC_STREAM; // BEING DONE AT THE END OF evalRHS
  }
#undef mu
#undef la
#undef u
#undef lu
#undef strx
#undef stry
#undef strz
#undef acof
#undef bope
#undef ghcof
}

//-----------------------------------------------------------------------
void rhserrfort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, int nz, float_sw4 h, float_sw4* __restrict__ a_fo,
                   float_sw4* __restrict__ a_u2, float_sw4 lowZ[3],
                   float_sw4 interZ[3], float_sw4 highZ[3]) {
  SW4_MARK_FUNCTION;
#define fo(c, i, j, k) a_fo[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define u2(c, i, j, k) a_u2[base3 + i + ni * (j) + nij * (k) + nijk * (c)]

  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;

  for (int c = 1; c <= 3; c++) {
    float_sw4 lowz = 0;
#pragma omp parallel
#pragma omp for reduction(max : lowz)
    for (int k = 1; k <= 6; k++)
      for (int j = jfirst + 2; j <= jlast - 2; j++)
        for (int i = ifirst + 2; i <= ilast - 2; i++) {
          float_sw4 er = abs(fo(c, i, j, k) - u2(c, i, j, k));
          if (lowz < er) lowz = er;
        }
    lowZ[c - 1] = lowz;
  }

  for (int c = 1; c <= 3; c++) {
    float_sw4 interz = 0;
#pragma omp parallel
#pragma omp for reduction(max : interz)
    for (int k = 7; k <= nz - 6; k++)
      for (int j = jfirst + 2; j <= jlast - 2; j++)
        for (int i = ifirst + 2; i <= ilast - 2; i++) {
          float_sw4 er = abs(fo(c, i, j, k) - u2(c, i, j, k));
          if (interz < er) interz = er;
        }
    interZ[c - 1] = interz;
  }

  for (int c = 1; c <= 3; c++) {
    float_sw4 highz = 0;
#pragma omp parallel
#pragma omp for reduction(max : highz)
    for (int k = nz - 5; k <= nz; k++)
      for (int j = jfirst + 2; j <= jlast - 2; j++)
        for (int i = ifirst + 2; i <= ilast - 2; i++) {
          float_sw4 er = abs(fo(c, i, j, k) - u2(c, i, j, k));
          if (highz < er) highz = er;
        }
    highZ[c - 1] = highz;
  }
#undef fo
#undef u2
}

//-----------------------------------------------------------------------
void rhouttlumf_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, int nz, float_sw4* __restrict__ a_uacc,
                   float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_fo,
                   float_sw4* __restrict__ a_rho, float_sw4 lowZ[3],
                   float_sw4 interZ[3], float_sw4 highZ[3]) {
  SW4_MARK_FUNCTION;
#define fo(c, i, j, k) a_fo[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define uacc(c, i, j, k) a_uacc[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define lu(c, i, j, k) a_lu[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define rho(i, j, k) a_rho[base + i + ni * (j) + nij * (k)]

  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;

  for (int c = 1; c <= 3; c++) {
    float_sw4 lowz = 0;
#pragma omp parallel
#pragma omp for reduction(max : lowz)
    for (int k = 1; k <= 6; k++)
      for (int j = jfirst + 2; j <= jlast - 2; j++)
        for (int i = ifirst + 2; i <= ilast - 2; i++) {
          float_sw4 er = abs(rho(i, j, k) * uacc(c, i, j, k) - lu(c, i, j, k) -
                             fo(c, i, j, k));
          if (lowz < er) lowz = er;
        }
    lowZ[c - 1] = lowz;
  }

  for (int c = 1; c <= 3; c++) {
    float_sw4 interz = 0;
#pragma omp parallel
#pragma omp for reduction(max : interz)
    for (int k = 7; k <= nz - 6; k++)
      for (int j = jfirst + 2; j <= jlast - 2; j++)
        for (int i = ifirst + 2; i <= ilast - 2; i++) {
          float_sw4 er = abs(rho(i, j, k) * uacc(c, i, j, k) - lu(c, i, j, k) -
                             fo(c, i, j, k));
          if (interz < er) interz = er;
        }
    interZ[c - 1] = interz;
  }

  for (int c = 1; c <= 3; c++) {
    float_sw4 highz = 0;
#pragma omp parallel
#pragma omp for reduction(max : highz)
    for (int k = nz - 5; k <= nz; k++)
      for (int j = jfirst + 2; j <= jlast - 2; j++)
        for (int i = ifirst + 2; i <= ilast - 2; i++) {
          float_sw4 er = abs(rho(i, j, k) * uacc(c, i, j, k) - lu(c, i, j, k) -
                             fo(c, i, j, k));
          if (highz < er) highz = er;
        }
    highZ[c - 1] = highz;
  }
#undef fo
#undef uacc
#undef lu
#undef rho
}

//-----------------------------------------------------------------------
void predfort_ci(int ib, int ie, int jb, int je, int kb, int ke,
                 float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                 float_sw4* __restrict__ um, float_sw4* __restrict__ lu,
                 float_sw4* __restrict__ fo, float_sw4* __restrict__ rho,
                 float_sw4 dt2) {
  SW4_MARK_FUNCTION;
  const size_t npts =
      static_cast<size_t>((ie - ib + 1)) * (je - jb + 1) * (ke - kb + 1);
#// pragma omp parallel for                  \
    // #pragma ivdep                         \
    // #pragma simd                          \
    //    for( size_t i=0 ; i < npts ; i++ ) \
    //    {
  ASSERT_MANAGED(up);
  ASSERT_MANAGED(um);
  ASSERT_MANAGED(fo);
  ASSERT_MANAGED(u);
  ASSERT_MANAGED(lu);
  ASSERT_MANAGED(rho);
  std::cout<<"RUNNING ON CPU\n";
  for( size_t i=0 ; i < npts ; i++ ) {
float_sw4 dt2orh = dt2 / rho[i];
        up[i] = 2 * u[i] - um[i] + dt2orh * (lu[i] + fo[i]);
        up[i + npts] = 2 * u[i + npts] - um[i + npts] +
                       dt2orh * (lu[i + npts] + fo[i + npts]);
        up[i + 2 * npts] = 2 * u[i + 2 * npts] - um[i + 2 * npts] +
                           dt2orh * (lu[i + 2 * npts] + fo[i + 2 * npts]);
	std::cout<<"PRED "<<i<<" "<<up[i+npts]<<" "<<u[i+npts]<<" "<<um[i+npts]<<" "<<lu[i+npts]<<" "<<fo[i+npts]<<"\n";
  }
  return;
    
  RAJA::forall<PREDFORT_LOOP_POL_ASYNC>(
					//RAJA::forall<RAJA::seq_exec>(
      RAJA::RangeSegment(0, npts), [=] RAJA_DEVICE(size_t i) {
        float_sw4 dt2orh = dt2 / rho[i];
        up[i] = 2 * u[i] - um[i] + dt2orh * (lu[i] + fo[i]);
        up[i + npts] = 2 * u[i + npts] - um[i + npts] +
                       dt2orh * (lu[i + npts] + fo[i + npts]);
        up[i + 2 * npts] = 2 * u[i + 2 * npts] - um[i + 2 * npts] +
                           dt2orh * (lu[i + 2 * npts] + fo[i + 2 * npts]);
      });  // SYNC_STREAM;
}

//-----------------------------------------------------------------------
void corrfort_ci(int ib, int ie, int jb, int je, int kb, int ke,
                 float_sw4* __restrict__ up, float_sw4* __restrict__ lu,
                 float_sw4* __restrict__ fo, float_sw4* __restrict__ rho,
                 float_sw4 dt4) {
  SW4_MARK_FUNCTION;
  const float_sw4 dt4i12 = dt4 / 12;
  const size_t npts =
      static_cast<size_t>((ie - ib + 1)) * (je - jb + 1) * (ke - kb + 1);
  // #pragma omp parallel for
  // #pragma ivdep
  // #pragma simd
  //    for( size_t i=0 ; i < npts ; i++ )
  //    {
  RAJA::forall<CORRFORT_LOOP_POL_ASYNC>(
      RAJA::RangeSegment(0, npts), [=] RAJA_DEVICE(size_t i) {
        float_sw4 dt4i12orh = dt4i12 / rho[i];
        up[i] += dt4i12orh * (lu[i] + fo[i]);
        up[i + npts] += dt4i12orh * (lu[i + npts] + fo[i + npts]);
        up[i + 2 * npts] += dt4i12orh * (lu[i + 2 * npts] + fo[i + 2 * npts]);
      });  // SYNC_STREAM;
}

//-----------------------------------------------------------------------
void dpdmtfort_ci(int ib, int ie, int jb, int je, int kb, int ke,
                  float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                  float_sw4* __restrict__ um, float_sw4* __restrict__ u2,
                  float_sw4 dt2i) {
  SW4_MARK_FUNCTION;
  const size_t npts =
      static_cast<size_t>((ie - ib + 1)) * (je - jb + 1) * (ke - kb + 1);
  // #pragma omp parallel for
  // #pragma ivdep
  // #pragma simd
  //   for( size_t i = 0 ; i < 3*npts ; i++ )
  RAJA::forall<DPDMTFORT_LOOP_POL>(RAJA::RangeSegment(0, 3 * npts),
                                   [=] RAJA_DEVICE(size_t i) {
                                     u2[i] = dt2i * (up[i] - 2 * u[i] + um[i]);
                                   });  // SYNC_STREAM;
}

//-----------------------------------------------------------------------
void dpdmtfortatt_ci(int ib, int ie, int jb, int je, int kb, int ke,
                     float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                     float_sw4* __restrict__ um, float_sw4 dt2i) {
  SW4_MARK_FUNCTION;
  const size_t npts =
      static_cast<size_t>((ie - ib + 1)) * (je - jb + 1) * (ke - kb + 1);
  // #pragma omp parallel for
  // #pragma ivdep
  // #pragma simd
  //   for( size_t i = 0 ; i < 3*npts ; i++ )
  RAJA::forall<DPDMTFORT_LOOP_POL_ASYNC>(
      RAJA::RangeSegment(0, 3 * npts), [=] RAJA_DEVICE(size_t i) {
        um[i] = dt2i * (up[i] - 2 * u[i] + um[i]);
      });  // SYNC_STREAM;

  // Forallasync is the same as RAJA:;forall
  // forallasync(0,3*npts,[=] RAJA_DEVICE(size_t i){
  //     um[i] = dt2i*(up[i]-2*u[i]+um[i]);}); //SYNC_STREAM;
}

// //-----------------------------------------------------------------------
// void EW::updatememvar_ci( int ifirst, int ilast, int jfirst, int jlast, int
// kfirst, int klast, 			  float_sw4* __restrict__ a_alp,
// float_sw4*
// __restrict__ a_alm, 			  float_sw4* __restrict__ a_up,
// float_sw4* __restrict__ a_u, 			  float_sw4*
// __restrict__ a_um, 			  float_sw4 omega, float_sw4 dt, int
// domain
// )
// {
//    const float_sw4 i6 = 1.0/6;
//    const int ni    = ilast-ifirst+1;
//    const int nij   = ni*(jlast-jfirst+1);
//    const int nijk  = nij*(klast-kfirst+1);
//    const int base  = -(ifirst+ni*jfirst+nij*kfirst);
//    const int base3 = base-nijk;

// #define alp(c,i,j,k) a_alp[base3+i+ni*(j)+nij*(k)+nijk*(c)]
// #define alm(c,i,j,k) a_alm[base3+i+ni*(j)+nij*(k)+nijk*(c)]
// #define up(c,i,j,k)   a_up[base3+i+ni*(j)+nij*(k)+nijk*(c)]
// #define u(c,i,j,k)     a_u[base3+i+ni*(j)+nij*(k)+nijk*(c)]
// #define um(c,i,j,k)   a_um[base3+i+ni*(j)+nij*(k)+nijk*(c)]

//    float_sw4 dto = dt*omega;
//    float_sw4 icp = 1.0/( 1.0/2 + 1/(2*dto) + dto/4 + dto*dto/12 );
//    float_sw4 cm  =     1.0/2 - 1/(2*dto) - dto/4 + dto*dto/12;

//    int k1, k2;
//    if( domain == 0 )
//    {
//       k1 = kfirst;
//       k2 = klast;
//    }
//    else if( domain==1 )
//    {
//       k1 = kfirst;
//       k2 = kfirst+2;
//    }
//    else if( domain == 2 )
//    {
//       k1 = klast-2;
//       k2 = klast;
//    }

//    for( int c=1 ; c<=3 ; c++ )
// #pragma omp parallel
// #pragma omp for
//       for( int k=k1 ; k<=k2 ; k++ )
// 	 for( int j=jfirst ; j<=jlast ; j++ )
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst ; i<=ilast ; i++ )
// 	    {
// 	       alp(c,i,j,k) = icp*(-cm*alm(c,i,j,k) + u(c,i,j,k) +
// 			  i6*( dto*dto*u(c,i,j,k) +
// dto*(up(c,i,j,k)-um(c,i,j,k))
// + 			       (up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k))  ) );
// 	    }
// #undef alp
// #undef alm
// #undef up
// #undef u
// #undef um
// }

//-----------------------------------------------------------------------
void satt_ci(float_sw4* __restrict__ up, float_sw4* __restrict__ qs,
             float_sw4 dt, float_sw4 cfreq, int ifirst, int ilast, int jfirst,
             int jlast, int kfirst, int klast) {
  SW4_MARK_FUNCTION;
  const size_t npts = static_cast<size_t>((ilast - ifirst + 1)) *
                      (jlast - jfirst + 1) * (klast - kfirst + 1);
  const float_sw4 efact = M_PI * cfreq * dt;
#pragma omp parallel for
#pragma ivdep
#pragma simd
  for (size_t i = 0; i < npts; i++) {
    float_sw4 fact = exp(-efact / qs[i]);
    up[i] *= fact;
    up[i + npts] *= fact;
    up[i + 2 * npts] *= fact;
  }
}

//-----------------------------------------------------------------------
void solveattfreeac_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* __restrict__ a_alpha,
                       float_sw4 cof, float_sw4* __restrict__ a_up) {
  SW4_MARK_FUNCTION;
#define alpha(c, i, j, k) a_alpha[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define up(c, i, j, k) a_up[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;

#pragma omp parallel
  {
    int k = 0;
#pragma omp for
    for (int j = jfirst + 2; j <= jlast - 2; j++)
#pragma ivdep
#pragma simd
      for (int i = ifirst + 2; i <= ilast - 2; i++)
        alpha(1, i, j, k) += cof * up(1, i, j, k);
#pragma omp for
    for (int j = jfirst + 2; j <= jlast - 2; j++)
#pragma ivdep
#pragma simd
      for (int i = ifirst + 2; i <= ilast - 2; i++)
        alpha(2, i, j, k) += cof * up(2, i, j, k);
#pragma omp for
    for (int j = jfirst + 2; j <= jlast - 2; j++)
#pragma ivdep
#pragma simd
      for (int i = ifirst + 2; i <= ilast - 2; i++)
        alpha(3, i, j, k) += cof * up(3, i, j, k);
  }
#undef alpha
#undef up
}

//-----------------------------------------------------------------------
void solveattfreec_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_mu,
    float_sw4* __restrict__ a_la, float_sw4* __restrict__ a_muve,
    float_sw4* __restrict__ a_lave, float_sw4* __restrict__ a_bforcerhs,
    float_sw4* __restrict__ a_met, float_sw4 s[5], int usesg,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry) {
  SW4_MARK_FUNCTION;
#define u(c, i, j, k) a_u[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define met(c, i, j, k) a_met[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define mu(i, j, k) a_mu[base + i + ni * (j) + nij * (k)]
#define la(i, j, k) a_la[base + i + ni * (j) + nij * (k)]
#define muve(i, j) a_muve[base0 + i + ni * (j)]
#define lambdave(i, j) a_lave[base0 + i + ni * (j)]
#define sgstrx(i) a_strx[i - ifirst]
#define sgstry(j) a_stry[j - jfirst]
#define bforcerhs(c, i, j) a_bforcerhs[base03 + (i) + ni * (j) + nij * (c)]

  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base0 = -(ifirst + ni * jfirst);
  const int base03 = base0 - nij;
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;

#pragma omp parallel
  {
    float_sw4 sgx = 1, sgy = 1, isgx = 1, isgy = 1, s0i = 1 / s[0];
    int k = 1, kl = 1;
#pragma omp for
    for (int j = jfirst + 2; j <= jlast - 2; j++)
#pragma ivdep
#pragma simd
      for (int i = ifirst + 2; i <= ilast - 2; i++) {
        float_sw4 mupt = mu(i, j, k) - muve(i, j);
        float_sw4 lapt = la(i, j, k) - lambdave(i, j);
        if (usesg == 1) {
          sgx = sgstrx(i);
          sgy = sgstry(j);
          isgy = 1 / sgy;
          isgx = 1 / sgx;
        }
        float_sw4 m2sg = sqrt(sgx * isgy);
        float_sw4 m3sg = 1 / m2sg;
        float_sw4 m4sg = isgx * m2sg;

        float_sw4 ac = sgx * isgy * met(2, i, j, k) * met(2, i, j, k) +
                       isgx * sgy * met(3, i, j, k) * met(3, i, j, k) +
                       isgx * isgy * met(4, i, j, k) * met(4, i, j, k);
        float_sw4 bc = 1 / (mupt * ac);
        float_sw4 cc = (mupt + lapt) / (2 * mupt + lapt) * bc / ac;
        float_sw4 dc = cc * (m2sg * met(2, i, j, k) * bforcerhs(1, i, j) +
                             m3sg * met(3, i, j, k) * bforcerhs(2, i, j) +
                             m4sg * met(4, i, j, k) * bforcerhs(3, i, j));
        u(1, i, j, k - kl) =
            s0i * (bc * bforcerhs(1, i, j) - dc * met(2, i, j, k) * m2sg);
        u(2, i, j, k - kl) =
            s0i * (bc * bforcerhs(2, i, j) - dc * met(3, i, j, k) * m3sg);
        u(3, i, j, k - kl) =
            s0i * (bc * bforcerhs(3, i, j) - dc * met(4, i, j, k) * m4sg);
      }
  }
#undef u
#undef mu
#undef la
#undef met
#undef muve
#undef lambdave
#undef sgstrx
#undef sgstry
#undef bforcerhs
}

//-----------------------------------------------------------------------
void addbstresswresc_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nz,
    float_sw4* __restrict__ a_alphap, float_sw4* __restrict__ a_alpham,
    float_sw4* __restrict__ a_muve, float_sw4* __restrict__ a_lave,
    float_sw4* __restrict__ a_bforcerhs, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_met, int side,
    float_sw4 dt, float_sw4 omegave, float_sw4* __restrict__ a_memforce,
    float_sw4* __restrict__ a_muvebnd, float_sw4* __restrict__ a_lambdavebnd,
    float_sw4 s[5], float_sw4& cof, int usesg, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry) {
  SW4_MARK_FUNCTION;
#define u(c, i, j, k) a_u[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define um(c, i, j, k) a_um[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define met(c, i, j, k) a_met[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define alphap(c, i, j, k) \
  a_alphap[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define alpham(c, i, j, k) \
  a_alpham[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define muve(i, j, k) a_muve[base + i + ni * (j) + nij * (k)]
#define lave(i, j, k) a_lave[base + i + ni * (j) + nij * (k)]
#define muvebnd(i, j) a_muvebnd[base0 + i + ni * (j)]
#define lambdavebnd(i, j) a_lambdavebnd[base0 + i + ni * (j)]
#define sgstrx(i) a_strx[i - ifirst]
#define sgstry(j) a_stry[j - jfirst]
  //#define bforcerhs(c,i,j) a_bforcerhs[base0+c-1+3*(i+ni*(j))]
  //#define memforce(c,i,j)   a_memforce[base0+c-1+3*(i+ni*(j))]
#define bforcerhs(c, i, j) a_bforcerhs[base03 + (i) + ni * (j) + nij * (c)]
#define memforce(c, i, j) a_memforce[base03 + (i) + ni * (j) + nij * (c)]

  const float_sw4 i6 = 1.0 / 6;
  const float_sw4 c1 = 2.0 / 3;
  const float_sw4 c2 = -1.0 / 12;
  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  //   const int base0 = -(ifirst+ni*jfirst);
  const int base0 = -(ifirst + ni * jfirst);
  const int base03 = base0 - nij;
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;
  int k, kl;
  if (side == 5) {
    k = 1;
    kl = 1;
  }
  if (side == 6) {
    k = nz;
    kl = -1;
  }
  float_sw4 omdt = omegave * dt;
  float_sw4 cp = 0.5 + 1 / (2 * omdt) + omdt / 4 + omdt * omdt / 12;
  float_sw4 cm = 0.5 - 1 / (2 * omdt) - omdt / 4 + omdt * omdt / 12;
  cof = (omdt + 1) / (6 * cp);
#pragma omp parallel
  {
    float_sw4 sgx = 1, sgy = 1, isgx = 1, isgy = 1;
#pragma omp for
    for (int j = jfirst + 2; j <= jlast - 2; j++)
#pragma ivdep
#pragma simd
      for (int i = ifirst + 2; i <= ilast - 2; i++) {
        float_sw4 r1 =
            (-cm * alpham(1, i, j, k - kl) +
             (4 + omdt * omdt) * i6 * u(1, i, j, k - kl) +
             i6 * (1 - omdt) * um(1, i, j, k - kl) + memforce(1, i, j)) /
            cp;
        float_sw4 r2 =
            (-cm * alpham(2, i, j, k - kl) +
             (4 + omdt * omdt) * i6 * u(2, i, j, k - kl) +
             i6 * (1 - omdt) * um(2, i, j, k - kl) + memforce(2, i, j)) /
            cp;
        float_sw4 r3 =
            (-cm * alpham(3, i, j, k - kl) +
             (4 + omdt * omdt) * i6 * u(3, i, j, k - kl) +
             i6 * (1 - omdt) * um(3, i, j, k - kl) + memforce(3, i, j)) /
            cp;
        alphap(1, i, j, k - kl) = r1;
        alphap(2, i, j, k - kl) = r2;
        alphap(3, i, j, k - kl) = r3;
        muvebnd(i, j) += cof * muve(i, j, k);
        lambdavebnd(i, j) += cof * lave(i, j, k);

        //	    sgx = usesg*sgstrx(i)+(1-usesg)*1;
        //	    sgx = usesg*sgstry(j)+(1-usesg)*1;
        //	    isgx =1/sgx;
        //	    isgy =1/sgy;
        if (usesg == 1) {
          sgx = sgstrx(i);
          sgy = sgstry(j);
          isgy = 1 / sgy;
          isgx = 1 / sgx;
        }

        //   tangential derivatives
        float_sw4 rhs1 =
            // pr
            (2 * muve(i, j, k) + lave(i, j, k)) * met(2, i, j, k) *
                met(1, i, j, k) *
                (c2 * (alphap(1, i + 2, j, k) - alphap(1, i - 2, j, k)) +
                 c1 * (alphap(1, i + 1, j, k) - alphap(1, i - 1, j, k))) *
                sgx * isgy +
            muve(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(2, i + 2, j, k) - alphap(2, i - 2, j, k)) +
                 c1 * (alphap(2, i + 1, j, k) - alphap(2, i - 1, j, k))) +
            muve(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(3, i + 2, j, k) - alphap(3, i - 2, j, k)) +
                 c1 * (alphap(3, i + 1, j, k) - alphap(3, i - 1, j, k))) *
                isgy
            // qr
            + muve(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                  (c2 * (alphap(1, i, j + 2, k) - alphap(1, i, j - 2, k)) +
                   c1 * (alphap(1, i, j + 1, k) - alphap(1, i, j - 1, k))) *
                  isgx * sgy +
            lave(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(2, i, j + 2, k) - alphap(2, i, j - 2, k)) +
                 c1 * (alphap(2, i, j + 1, k) - alphap(2, i, j - 1, k)));

        // (v-eq)
        float_sw4 rhs2 =
            // pr
            lave(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(1, i + 2, j, k) - alphap(1, i - 2, j, k)) +
                 c1 * (alphap(1, i + 1, j, k) - alphap(1, i - 1, j, k))) +
            muve(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(2, i + 2, j, k) - alphap(2, i - 2, j, k)) +
                 c1 * (alphap(2, i + 1, j, k) - alphap(2, i - 1, j, k))) *
                sgx * isgy
            // qr
            + muve(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                  (c2 * (alphap(1, i, j + 2, k) - alphap(1, i, j - 2, k)) +
                   c1 * (alphap(1, i, j + 1, k) - alphap(1, i, j - 1, k))) +
            (2 * muve(i, j, k) + lave(i, j, k)) * met(3, i, j, k) *
                met(1, i, j, k) *
                (c2 * (alphap(2, i, j + 2, k) - alphap(2, i, j - 2, k)) +
                 c1 * (alphap(2, i, j + 1, k) - alphap(2, i, j - 1, k))) *
                sgy * isgx +
            muve(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(3, i, j + 2, k) - alphap(3, i, j - 2, k)) +
                 c1 * (alphap(3, i, j + 1, k) - alphap(3, i, j - 1, k))) *
                isgx;

        // (w-eq)
        float_sw4 rhs3 =
            // pr
            lave(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(1, i + 2, j, k) - alphap(1, i - 2, j, k)) +
                 c1 * (alphap(1, i + 1, j, k) - alphap(1, i - 1, j, k))) *
                isgy +
            muve(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(3, i + 2, j, k) - alphap(3, i - 2, j, k)) +
                 c1 * (alphap(3, i + 1, j, k) - alphap(3, i - 1, j, k))) *
                sgx * isgy
            // qr
            + muve(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                  (c2 * (alphap(3, i, j + 2, k) - alphap(3, i, j - 2, k)) +
                   c1 * (alphap(3, i, j + 1, k) - alphap(3, i, j - 1, k))) *
                  sgy * isgx +
            lave(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(2, i, j + 2, k) - alphap(2, i, j - 2, k)) +
                 c1 * (alphap(2, i, j + 1, k) - alphap(2, i, j - 1, k))) *
                isgx;

        // normal derivatives
        float_sw4 un1 = s[1] * alphap(1, i, j, k) +
                        s[2] * alphap(1, i, j, k + kl) +
                        s[3] * alphap(1, i, j, k + 2 * kl) +
                        s[4] * alphap(1, i, j, k + 3 * kl) + s[0] * r1;
        float_sw4 vn1 = s[1] * alphap(2, i, j, k) +
                        s[2] * alphap(2, i, j, k + kl) +
                        s[3] * alphap(2, i, j, k + 2 * kl) +
                        s[4] * alphap(2, i, j, k + 3 * kl) + s[0] * r2;
        float_sw4 wn1 = s[1] * alphap(3, i, j, k) +
                        s[2] * alphap(3, i, j, k + kl) +
                        s[3] * alphap(3, i, j, k + 2 * kl) +
                        s[4] * alphap(3, i, j, k + 3 * kl) + s[0] * r3;

        float_sw4 m2sg = sqrt(sgx * isgy);
        float_sw4 m3sg = 1 / m2sg;
        float_sw4 m4sg = isgx * m2sg;

        float_sw4 rtu = un1 * m2sg * met(2, i, j, k) +
                        vn1 * m3sg * met(3, i, j, k) +
                        wn1 * m4sg * met(4, i, j, k);
        float_sw4 ac = sgx * isgy * met(2, i, j, k) * met(2, i, j, k) +
                       sgy * isgx * met(3, i, j, k) * met(3, i, j, k) +
                       isgx * isgy * met(4, i, j, k) * met(4, i, j, k);
        rhs1 = rhs1 +
               (muve(i, j, k) + lave(i, j, k)) * rtu * m2sg * met(2, i, j, k) +
               muve(i, j, k) * ac * un1;
        rhs2 = rhs2 +
               (muve(i, j, k) + lave(i, j, k)) * rtu * m3sg * met(3, i, j, k) +
               muve(i, j, k) * ac * vn1;
        rhs3 = rhs3 +
               (muve(i, j, k) + lave(i, j, k)) * rtu * m4sg * met(4, i, j, k) +
               muve(i, j, k) * ac * wn1;
        bforcerhs(1, i, j) += rhs1;
        bforcerhs(2, i, j) += rhs2;
        bforcerhs(3, i, j) += rhs3;
      }
  }
#undef u
#undef um
#undef met
#undef alphap
#undef alpham
#undef muve
#undef lave
#undef muvebnd
#undef lambdavebnd
#undef sgstrx
#undef sgstry
#undef bforcerhs
#undef memforce
}

//-----------------------------------------------------------------------
void ve_bndry_stress_curvi_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nz,
    float_sw4* __restrict__ a_alphap, float_sw4* __restrict__ a_muve,
    float_sw4* __restrict__ a_lave, float_sw4* __restrict__ a_bforcerhs,
    float_sw4* __restrict__ a_met, int side, float_sw4* __restrict__ sbop,
    int usesg, float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry) {
  SW4_MARK_FUNCTION;
  SW4_MARK_BEGIN("HOST CODE");
#define alphap(c, i, j, k) \
  a_alphap[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define muve(i, j, k) a_muve[base + i + ni * (j) + nij * (k)]
#define lave(i, j, k) a_lave[base + i + ni * (j) + nij * (k)]
#define met(c, i, j, k) a_met[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define sgstrx(i) a_strx[i - ifirst]
#define sgstry(j) a_stry[j - jfirst]
  //#define bforcerhs(c,i,j) a_bforcerhs[base0+c-1+3*(i+ni*(j))]
  //#define memforce(c,i,j)   a_memforce[base0+c-1+3*(i+ni*(j))]
#define bforcerhs(c, i, j) a_bforcerhs[base03 + (i) + ni * (j) + nij * (c)]

  // const float_sw4 i6 = 1.0/6;
  const float_sw4 c1 = 2.0 / 3;
  const float_sw4 c2 = -1.0 / 12;
  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  //   const int base0 = -(ifirst+ni*jfirst);
  const int base0 = -(ifirst + ni * jfirst);
  const int base03 = base0 - nij;
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;
  int k, kl;
  if (side == 5) {
    k = 1;
    kl = 1;
  }
  if (side == 6) {
    k = nz;
    kl = -1;
  }
  //#pragma omp parallel
  {
    //
// #pragma omp for
//       for( int j=jfirst+2 ; j<=jlast-2 ; j++ )
// #pragma ivdep
// 	 for( int i=ifirst+2 ; i<=ilast-2 ; i++ )
// 	 {
#ifdef ENABLE_CUDA
    using LOCAL_POL_ORIG =
        RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::For<
            0, RAJA::cuda_threadblock_exec<16>,
            RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
                                 RAJA::statement::Lambda<0>>>>>;

    using LOCAL_POL =
        RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::For<
            0, RAJA::cuda_block_exec,
            RAJA::statement::For<1, RAJA::cuda_block_exec,
                                 RAJA::statement::Lambda<0>>>>>;
#else
    using LOCAL_POL = DEFAULT_LOOP2;
#endif
    RAJA::RangeSegment i_range(ifirst + 2, ilast - 1);
    RAJA::RangeSegment j_range(jfirst + 2, jlast - 1);
    SW4_MARK_END("HOST CODE");
#ifdef ENABLE_CUDA
#define NO_COLLAPSE 1
#endif
#ifdef NO_COLLAPSE
    Range<16> I(ifirst + 2, ilast - 1);
    Range<4> J(jfirst + 2, jlast - 1);
    forall2async(I, J, [=] RAJA_DEVICE(int i, int j) {
#else
    RAJA::kernel<LOCAL_POL>(RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(
                                                                    int j,
                                                                    int i) {
#endif
      float_sw4 sgx = 1, sgy = 1, isgx = 1, isgy = 1;
      if (usesg == 1) {
        sgx = sgstrx(i);
        sgy = sgstry(j);
        isgy = 1 / sgy;
        isgx = 1 / sgx;
      }

      //   tangential derivatives
      float_sw4 rhs1 =
          // pr
          (2 * muve(i, j, k) + lave(i, j, k)) * met(2, i, j, k) *
              met(1, i, j, k) *
              (c2 * (alphap(1, i + 2, j, k) - alphap(1, i - 2, j, k)) +
               c1 * (alphap(1, i + 1, j, k) - alphap(1, i - 1, j, k))) *
              sgx * isgy +
          muve(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(2, i + 2, j, k) - alphap(2, i - 2, j, k)) +
               c1 * (alphap(2, i + 1, j, k) - alphap(2, i - 1, j, k))) +
          muve(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(3, i + 2, j, k) - alphap(3, i - 2, j, k)) +
               c1 * (alphap(3, i + 1, j, k) - alphap(3, i - 1, j, k))) *
              isgy
          // qr
          + muve(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(1, i, j + 2, k) - alphap(1, i, j - 2, k)) +
                 c1 * (alphap(1, i, j + 1, k) - alphap(1, i, j - 1, k))) *
                isgx * sgy +
          lave(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(2, i, j + 2, k) - alphap(2, i, j - 2, k)) +
               c1 * (alphap(2, i, j + 1, k) - alphap(2, i, j - 1, k)));

      // (v-eq)
      float_sw4 rhs2 =
          // pr
          lave(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(1, i + 2, j, k) - alphap(1, i - 2, j, k)) +
               c1 * (alphap(1, i + 1, j, k) - alphap(1, i - 1, j, k))) +
          muve(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(2, i + 2, j, k) - alphap(2, i - 2, j, k)) +
               c1 * (alphap(2, i + 1, j, k) - alphap(2, i - 1, j, k))) *
              sgx * isgy
          // qr
          + muve(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(1, i, j + 2, k) - alphap(1, i, j - 2, k)) +
                 c1 * (alphap(1, i, j + 1, k) - alphap(1, i, j - 1, k))) +
          (2 * muve(i, j, k) + lave(i, j, k)) * met(3, i, j, k) *
              met(1, i, j, k) *
              (c2 * (alphap(2, i, j + 2, k) - alphap(2, i, j - 2, k)) +
               c1 * (alphap(2, i, j + 1, k) - alphap(2, i, j - 1, k))) *
              sgy * isgx +
          muve(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(3, i, j + 2, k) - alphap(3, i, j - 2, k)) +
               c1 * (alphap(3, i, j + 1, k) - alphap(3, i, j - 1, k))) *
              isgx;

      // (w-eq)
      float_sw4 rhs3 =
          // pr
          lave(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(1, i + 2, j, k) - alphap(1, i - 2, j, k)) +
               c1 * (alphap(1, i + 1, j, k) - alphap(1, i - 1, j, k))) *
              isgy +
          muve(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(3, i + 2, j, k) - alphap(3, i - 2, j, k)) +
               c1 * (alphap(3, i + 1, j, k) - alphap(3, i - 1, j, k))) *
              sgx * isgy
          // qr
          + muve(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                (c2 * (alphap(3, i, j + 2, k) - alphap(3, i, j - 2, k)) +
                 c1 * (alphap(3, i, j + 1, k) - alphap(3, i, j - 1, k))) *
                sgy * isgx +
          lave(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
              (c2 * (alphap(2, i, j + 2, k) - alphap(2, i, j - 2, k)) +
               c1 * (alphap(2, i, j + 1, k) - alphap(2, i, j - 1, k))) *
              isgx;

      // normal derivatives
      float_sw4 un1 = sbop[0] * alphap(1, i, j, k - kl) +
                      sbop[1] * alphap(1, i, j, k) +
                      sbop[2] * alphap(1, i, j, k + kl) +
                      sbop[3] * alphap(1, i, j, k + 2 * kl) +
                      sbop[4] * alphap(1, i, j, k + 3 * kl) +
                      sbop[5] * alphap(1, i, j, k + 4 * kl);
      float_sw4 vn1 = sbop[0] * alphap(2, i, j, k - kl) +
                      sbop[1] * alphap(2, i, j, k) +
                      sbop[2] * alphap(2, i, j, k + kl) +
                      sbop[3] * alphap(2, i, j, k + 2 * kl) +
                      sbop[4] * alphap(2, i, j, k + 3 * kl) +
                      sbop[5] * alphap(2, i, j, k + 4 * kl);
      float_sw4 wn1 = sbop[0] * alphap(3, i, j, k - kl) +
                      sbop[1] * alphap(3, i, j, k) +
                      sbop[2] * alphap(3, i, j, k + kl) +
                      sbop[3] * alphap(3, i, j, k + 2 * kl) +
                      sbop[4] * alphap(3, i, j, k + 3 * kl) +
                      sbop[5] * alphap(3, i, j, k + 4 * kl);
      float_sw4 m2sg = sqrt(sgx * isgy);
      float_sw4 m3sg = 1 / m2sg;
      float_sw4 m4sg = isgx * m2sg;
      float_sw4 rtu = un1 * m2sg * met(2, i, j, k) +
                      vn1 * m3sg * met(3, i, j, k) +
                      wn1 * m4sg * met(4, i, j, k);
      float_sw4 ac = sgx * isgy * met(2, i, j, k) * met(2, i, j, k) +
                     sgy * isgx * met(3, i, j, k) * met(3, i, j, k) +
                     isgx * isgy * met(4, i, j, k) * met(4, i, j, k);
      rhs1 = rhs1 +
             (muve(i, j, k) + lave(i, j, k)) * rtu * m2sg * met(2, i, j, k) +
             muve(i, j, k) * ac * un1;
      rhs2 = rhs2 +
             (muve(i, j, k) + lave(i, j, k)) * rtu * m3sg * met(3, i, j, k) +
             muve(i, j, k) * ac * vn1;
      rhs3 = rhs3 +
             (muve(i, j, k) + lave(i, j, k)) * rtu * m4sg * met(4, i, j, k) +
             muve(i, j, k) * ac * wn1;
      bforcerhs(1, i, j) = rhs1 + bforcerhs(1, i, j);
      bforcerhs(2, i, j) = rhs2 + bforcerhs(2, i, j);
      bforcerhs(3, i, j) = rhs3 + bforcerhs(3, i, j);
    });  // SYNC_STREAM;
  }
}
#undef alphap
#undef muve
#undef lave
#undef met
#undef sgstrx
#undef sgstry
#undef bforcerhs

//-----------------------------------------------------------------------
void att_free_curvi_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_mu,
    float_sw4* __restrict__ a_lambda, float_sw4* __restrict__ a_bforcerhs,
    float_sw4* __restrict__ a_met, float_sw4* __restrict__ sbop, int usesg,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry) {
  SW4_MARK_FUNCTION;
#define u(c, i, j, k) a_u[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define mu(i, j, k) a_mu[base + i + ni * (j) + nij * (k)]
#define la(i, j, k) a_lambda[base + i + ni * (j) + nij * (k)]
#define met(c, i, j, k) a_met[base3 + i + ni * (j) + nij * (k) + nijk * (c)]
#define sgstrx(i) a_strx[i - ifirst]
#define sgstry(j) a_stry[j - jfirst]
#define bforcerhs(c, i, j) a_bforcerhs[base03 + (i) + ni * (j) + nij * (c)]

  const float_sw4 c1 = 2.0 / 3;
  const float_sw4 c2 = -1.0 / 12;
  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base0 = -(ifirst + ni * jfirst);
  const int base03 = base0 - nij;
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;

  // Hardcoded for the k=1 surface
  int k = 1, kl = 1;
  //#pragma omp parallel
  {
    //     float_sw4 sgx = 1, sgy = 1, isgx = 1, isgy = 1;
    float_sw4 s0i = 1 / sbop[0];
// #pragma omp for
//       for( int j=jfirst+2 ; j<=jlast-2 ; j++ )
// #pragma ivdep
// 	 for( int i=ifirst+2 ; i<=ilast-2 ; i++ )
// 	 {
#ifdef ENABLE_CUDA

#if SW4_RAJA_VERSION == 6
    using LOCAL_POL_ASYNC = RAJA::KernelPolicy<
        RAJA::statement::CudaKernelAsync<RAJA::statement::For<
            0, RAJA::cuda_threadblock_exec<16>,
            RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
                                 RAJA::statement::Lambda<0>>>>>;
#elif SW4_RAJA_VERSION == 7

    using LOCAL_POL_ASYNC = RAJA::KernelPolicy<
        RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
            0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
            RAJA::statement::Tile<
                1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
                RAJA::statement::For<
                    0, RAJA::cuda_thread_x_direct,
                    RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                                         RAJA::statement::Lambda<0>>>>>>>;

#endif
#else

    using LOCAL_POL_ASYNC = DEFAULT_LOOP2;

#endif
    RAJA::RangeSegment i_range(ifirst + 2, ilast - 1);
    RAJA::RangeSegment j_range(jfirst + 2, jlast - 1);
    RAJA::kernel<LOCAL_POL_ASYNC>(
        RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
          float_sw4 sgx = 1, sgy = 1, isgx = 1, isgy = 1;
          if (usesg == 1) {
            sgx = sgstrx(i);
            sgy = sgstry(j);
            isgx = 1 / sgx;
            isgy = 1 / sgy;
          }

          // First tangential derivatives
          float_sw4 rhs1 =
              // pr
              (2 * mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) *
                  met(1, i, j, k) *
                  (c2 * (u(1, i + 2, j, k) - u(1, i - 2, j, k)) +
                   c1 * (u(1, i + 1, j, k) - u(1, i - 1, j, k))) *
                  sgx * isgy +
              mu(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(2, i + 2, j, k) - u(2, i - 2, j, k)) +
                   c1 * (u(2, i + 1, j, k) - u(2, i - 1, j, k))) +
              mu(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(3, i + 2, j, k) - u(3, i - 2, j, k)) +
                   c1 * (u(3, i + 1, j, k) - u(3, i - 1, j, k))) *
                  isgy
              // qr
              + mu(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                    (c2 * (u(1, i, j + 2, k) - u(1, i, j - 2, k)) +
                     c1 * (u(1, i, j + 1, k) - u(1, i, j - 1, k))) *
                    isgx * sgy +
              la(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(2, i, j + 2, k) - u(2, i, j - 2, k)) +
                   c1 * (u(2, i, j + 1, k) - u(2, i, j - 1, k))) -
              bforcerhs(1, i, j);

          // (v-eq)
          float_sw4 rhs2 =
              // pr
              la(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(1, i + 2, j, k) - u(1, i - 2, j, k)) +
                   c1 * (u(1, i + 1, j, k) - u(1, i - 1, j, k))) +
              mu(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(2, i + 2, j, k) - u(2, i - 2, j, k)) +
                   c1 * (u(2, i + 1, j, k) - u(2, i - 1, j, k))) *
                  sgx * isgy
              // qr
              + mu(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                    (c2 * (u(1, i, j + 2, k) - u(1, i, j - 2, k)) +
                     c1 * (u(1, i, j + 1, k) - u(1, i, j - 1, k))) +
              (2 * mu(i, j, k) + la(i, j, k)) * met(3, i, j, k) *
                  met(1, i, j, k) *
                  (c2 * (u(2, i, j + 2, k) - u(2, i, j - 2, k)) +
                   c1 * (u(2, i, j + 1, k) - u(2, i, j - 1, k))) *
                  sgy * isgx +
              mu(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(3, i, j + 2, k) - u(3, i, j - 2, k)) +
                   c1 * (u(3, i, j + 1, k) - u(3, i, j - 1, k))) *
                  isgx -
              bforcerhs(2, i, j);

          // (w-eq)
          float_sw4 rhs3 =
              // pr
              la(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(1, i + 2, j, k) - u(1, i - 2, j, k)) +
                   c1 * (u(1, i + 1, j, k) - u(1, i - 1, j, k))) *
                  isgy +
              mu(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(3, i + 2, j, k) - u(3, i - 2, j, k)) +
                   c1 * (u(3, i + 1, j, k) - u(3, i - 1, j, k))) *
                  sgx * isgy
              // qr
              + mu(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                    (c2 * (u(3, i, j + 2, k) - u(3, i, j - 2, k)) +
                     c1 * (u(3, i, j + 1, k) - u(3, i, j - 1, k))) *
                    sgy * isgx +
              la(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(2, i, j + 2, k) - u(2, i, j - 2, k)) +
                   c1 * (u(2, i, j + 1, k) - u(2, i, j - 1, k))) *
                  isgx -
              bforcerhs(3, i, j);

          // Normal derivatives
          float_sw4 ac = sgx * isgy * met(2, i, j, k) * met(2, i, j, k) +
                         sgy * isgx * met(3, i, j, k) * met(3, i, j, k) +
                         met(4, i, j, k) * met(4, i, j, k) * isgy * isgx;
          float_sw4 bc = 1 / (mu(i, j, k) * ac);
          float_sw4 cc = (mu(i, j, k) + la(i, j, k)) /
                         (2 * mu(i, j, k) + la(i, j, k)) * bc / ac;

          float_sw4 xoysqrt = sqrt(sgx * isgy);
          float_sw4 yoxsqrt = 1 / xoysqrt;
          float_sw4 isqrtxy = isgx * xoysqrt;
          float_sw4 dc = cc * (xoysqrt * met(2, i, j, k) * rhs1 +
                               yoxsqrt * met(3, i, j, k) * rhs2 +
                               isqrtxy * met(4, i, j, k) * rhs3);

          u(1, i, j, k - kl) =
              -s0i * (sbop[1] * u(1, i, j, k) + sbop[2] * u(1, i, j, k + kl) +
                      sbop[3] * u(1, i, j, k + 2 * kl) +
                      sbop[4] * u(1, i, j, k + 3 * kl) + bc * rhs1 -
                      dc * met(2, i, j, k) * xoysqrt);
          u(2, i, j, k - kl) =
              -s0i * (sbop[1] * u(2, i, j, k) + sbop[2] * u(2, i, j, k + kl) +
                      sbop[3] * u(2, i, j, k + 2 * kl) +
                      sbop[4] * u(2, i, j, k + 3 * kl) + bc * rhs2 -
                      dc * met(3, i, j, k) * yoxsqrt);
          u(3, i, j, k - kl) =
              -s0i * (sbop[1] * u(3, i, j, k) + sbop[2] * u(3, i, j, k + kl) +
                      sbop[3] * u(3, i, j, k + 2 * kl) +
                      sbop[4] * u(3, i, j, k + 3 * kl) + bc * rhs3 -
                      dc * met(4, i, j, k) * isqrtxy);
        });  // SYNC_STREAM;
  }
#undef u
#undef mu
#undef la
#undef met
#undef sgstrx
#undef sgstry
#undef bforcerhs
}

//}
// #ifdef ENABLE_CUDA
// void rhs4th3fortsgstr_ciopt( int ifirst, int ilast, int jfirst, int jlast,
// int kfirst, int klast, 			     int ni,int nj,int nk,
// float_sw4* a_lu, float_sw4* a_u, 			     float_sw4* a_mu,
// float_sw4* a_lambda, 			     float_sw4 h, float_sw4*
// a_strx, float_sw4* a_stry, 			     float_sw4* a_strz, char op
// ){
//   SW4_MARK_FUNCTION;
//   int njcomp = jlast - jfirst + 1;
//   dim3 blocks = dim3((ni+BX-1)/BX, (njcomp+BY-1)/BY, 1);
//   dim3 threads = dim3(BX, BY, 1);

//   int a1;
//   float_sw4 cof = 1.0/(h*h);
//    if( op == '=' )
//       a1 = 0;
//    else if( op == '+' )
//       a1 = 1;
//    else if( op == '-' )
//    {
//       a1 = 1;
//       cof = -cof;
//    }
//    //std::cout<<"Calling optimized rhs4_v2\n";
//   rhs4_v2<1,0><<<blocks,threads,0,0>>>
//     (ifirst, ilast, jfirst, jlast, kfirst, klast,
//      ni, nj, nk,
//      a_lu, a_u,
//      a_mu, a_lambda,
//      a_strx, a_stry, a_strz, h,a1,cof);
//   //SW4_CheckDeviceError(cudaPeekAtLastError());
//   SW4_CheckDeviceError(cudaStreamSynchronize(0));
//   //std::cout<<"Done optimized rhs4_v2\n";
// }

// #endif
