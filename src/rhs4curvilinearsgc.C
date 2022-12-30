#include <cmath>

#include "EW.h"
#include "Mspace.h"
#include "caliper.h"
#include "sw4.h"

//-----------------------------------------------------------------------
void EW::freesurfcurvisg_ci(int ib, int ie, int jb, int je, int kb, int ke,
                            int nz, int side, float_sw4* __restrict__ a_u,
                            float_sw4* __restrict__ a_mu,
                            float_sw4* __restrict__ a_la,
                            float_sw4* __restrict__ a_met, float_sw4* s,
                            float_sw4* __restrict__ a_forcing,
                            float_sw4* __restrict__ a_strx,
                            float_sw4* __restrict__ a_stry) {
  SW4_MARK_FUNCTION;
  const float_sw4 c1 = 2.0 / 3, c2 = -1.0 / 12;

  const int ni = ie - ib + 1;
  const int nij = ni * (je - jb + 1);
  const int nijk = ni * (je - jb + 1) * (ke - kb + 1);
  const int base = -(ib + ni * jb + nij * kb);
  const int basef = -(ib + ni * jb);
  const int base4 = base - nijk;
  const int base3 = base - nijk;
  const int nic3 = 3 * ni;
#define mu(i, j, k) a_mu[base + i + ni * (j) + nij * (k)]
#define la(i, j, k) a_la[base + i + ni * (j) + nij * (k)]
#define met(c, i, j, k) a_met[base4 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define u(c, i, j, k) a_u[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define forcing(c, i, j) a_forcing[3 * basef - 1 + (c) + 3 * (i) + nic3 * (j)]
#define strx(i) a_strx[i - ib]
#define stry(j) a_stry[j - jb]

  int k, kl;
  if (side == 5) {
    k = 1;
    kl = 1;
  } else if (side == 6) {
    k = nz;
    kl = -1;
  }

  float_sw4 s0i = 1 / s[0];
  // #pragma omp parallel for
  //    for( int j= jb+2; j<=je-2 ; j++ )
  //    {

  // #pragma ivdep
  // #pragma simd
  //       for( int i= ib+2; i<=ie-2 ; i++ )
  //       {

  RAJA::RangeSegment i_range(ib + 2, ie - 1);
  RAJA::RangeSegment j_range(jb + 2, je - 1);
  RAJA::kernel<RHS4CU_POL_ASYNC>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
        float_sw4 istrx = 1 / strx(i);
        float_sw4 istry = 1 / stry(j);
        // First tangential derivatives
        float_sw4 rhs1 =
            // pr
            (2 * mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) *
                met(1, i, j, k) *
                (c2 * (u(1, i + 2, j, k) - u(1, i - 2, j, k)) +
                 c1 * (u(1, i + 1, j, k) - u(1, i - 1, j, k))) *
                strx(i) * istry +
            mu(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                (c2 * (u(2, i + 2, j, k) - u(2, i - 2, j, k)) +
                 c1 * (u(2, i + 1, j, k) - u(2, i - 1, j, k))) +
            mu(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                (c2 * (u(3, i + 2, j, k) - u(3, i - 2, j, k)) +
                 c1 * (u(3, i + 1, j, k) - u(3, i - 1, j, k))) *
                istry
            // qr
            + mu(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(1, i, j + 2, k) - u(1, i, j - 2, k)) +
                   c1 * (u(1, i, j + 1, k) - u(1, i, j - 1, k))) *
                  istrx * stry(j) +
            la(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                (c2 * (u(2, i, j + 2, k) - u(2, i, j - 2, k)) +
                 c1 * (u(2, i, j + 1, k) - u(2, i, j - 1, k))) -
            forcing(1, i, j);

        //(v-eq)
        float_sw4 rhs2 =
            // pr
            la(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                (c2 * (u(1, i + 2, j, k) - u(1, i - 2, j, k)) +
                 c1 * (u(1, i + 1, j, k) - u(1, i - 1, j, k))) +
            mu(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                (c2 * (u(2, i + 2, j, k) - u(2, i - 2, j, k)) +
                 c1 * (u(2, i + 1, j, k) - u(2, i - 1, j, k))) *
                strx(i) * istry
            // qr
            + mu(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(1, i, j + 2, k) - u(1, i, j - 2, k)) +
                   c1 * (u(1, i, j + 1, k) - u(1, i, j - 1, k))) +
            (2 * mu(i, j, k) + la(i, j, k)) * met(3, i, j, k) *
                met(1, i, j, k) *
                (c2 * (u(2, i, j + 2, k) - u(2, i, j - 2, k)) +
                 c1 * (u(2, i, j + 1, k) - u(2, i, j - 1, k))) *
                stry(j) * istrx +
            mu(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                (c2 * (u(3, i, j + 2, k) - u(3, i, j - 2, k)) +
                 c1 * (u(3, i, j + 1, k) - u(3, i, j - 1, k))) *
                istrx -
            forcing(2, i, j);

        // (w-eq)
        float_sw4 rhs3 =
            // pr
            la(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                (c2 * (u(1, i + 2, j, k) - u(1, i - 2, j, k)) +
                 c1 * (u(1, i + 1, j, k) - u(1, i - 1, j, k))) *
                istry +
            mu(i, j, k) * met(2, i, j, k) * met(1, i, j, k) *
                (c2 * (u(3, i + 2, j, k) - u(3, i - 2, j, k)) +
                 c1 * (u(3, i + 1, j, k) - u(3, i - 1, j, k))) *
                strx(i) * istry
            // qr
            + mu(i, j, k) * met(3, i, j, k) * met(1, i, j, k) *
                  (c2 * (u(3, i, j + 2, k) - u(3, i, j - 2, k)) +
                   c1 * (u(3, i, j + 1, k) - u(3, i, j - 1, k))) *
                  stry(j) * istrx +
            la(i, j, k) * met(4, i, j, k) * met(1, i, j, k) *
                (c2 * (u(2, i, j + 2, k) - u(2, i, j - 2, k)) +
                 c1 * (u(2, i, j + 1, k) - u(2, i, j - 1, k))) *
                istrx -
            forcing(3, i, j);

        // Normal derivatives
        float_sw4 ac = strx(i) * istry * met(2, i, j, k) * met(2, i, j, k) +
                       stry(j) * istrx * met(3, i, j, k) * met(3, i, j, k) +
                       met(4, i, j, k) * met(4, i, j, k) * istry * istrx;
        float_sw4 bc = 1 / (mu(i, j, k) * ac);
        float_sw4 cc = (mu(i, j, k) + la(i, j, k)) /
                       (2 * mu(i, j, k) + la(i, j, k)) * bc / ac;

        float_sw4 xoysqrt = sqrt(strx(i) * istry);
        float_sw4 yoxsqrt = 1 / xoysqrt;
        float_sw4 isqrtxy = istrx * xoysqrt;
        float_sw4 dc = cc * (xoysqrt * met(2, i, j, k) * rhs1 +
                             yoxsqrt * met(3, i, j, k) * rhs2 +
                             isqrtxy * met(4, i, j, k) * rhs3);

      u(1, i, j, k - kl) =
          -s0i *
          (s[1] * u(1, i, j, k) + s[2] * u(1, i, j, k + kl) +
           s[3] * u(1, i, j, k + 2 * kl) + s[4] * u(1, i, j, k + 3 * kl) +
           kl * bc * rhs1 - kl * dc * met(2, i, j, k) * xoysqrt);
      u(2, i, j, k - kl) =
          -s0i *
          (s[1] * u(2, i, j, k) + s[2] * u(2, i, j, k + kl) +
           s[3] * u(2, i, j, k + 2 * kl) + s[4] * u(2, i, j, k + 3 * kl) +
           kl * bc * rhs2 - kl * dc * met(3, i, j, k) * yoxsqrt);
      u(3, i, j, k - kl) =
          -s0i *
          (s[1] * u(3, i, j, k) + s[2] * u(3, i, j, k + kl) +
           s[3] * u(3, i, j, k + 2 * kl) + s[4] * u(3, i, j, k + 3 * kl) +
           kl * bc * rhs3 - kl * dc * met(4, i, j, k) * isqrtxy);
    }
  }
#undef mu
#undef la
#undef met
#undef u
#undef forcing
#undef strx
#undef stry
}

//-----------------------------------------------------------------------
void EW::getsurfforcingsg_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int k,
    float_sw4* __restrict__ a_met, float_sw4* __restrict__ a_jac,
    float_sw4* __restrict__ a_tau, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry, float_sw4* __restrict__ a_forcing) {
  SW4_MARK_FUNCTION;
  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = ni * (jlast - jfirst + 1) * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int basef = -(ifirst + ni * jfirst);
  const int base3 = base - nijk;
  // const int basef3= basef-nij;
  const int nic3 = 3 * ni;
  // const int nic6  = 6*ni;

#define met(c, i, j, k) a_met[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define jac(i, j, k) a_jac[base + (i) + ni * (j) + nij * (k)]
#define forcing(c, i, j) a_forcing[3 * basef - 1 + (c) + 3 * (i) + nic3 * (j)]
  //#define tau(c,i,j)         a_tau[6*basef-1+(c)+6*(i)+nic6*(j)]
#define tau(c, i, j) a_tau[basef - nij + (i) + ni * (j) + nij * (c)]
#define strx(i) a_strx[i - ifirst]
#define stry(j) a_stry[j - jfirst]

#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++) {
    float_sw4 istry = 1 / stry(j);
#pragma ivdep
    //#pragma simd
    for (int i = ifirst; i <= ilast; i++) {
      float_sw4 istrx = 1 / strx(i);
      float_sw4 sqjac = sqrt(jac(i, j, k));
      forcing(1, i, j) =
          sqjac * (istry * met(2, i, j, k) * tau(1, i, j) +
                   istrx * met(3, i, j, k) * tau(2, i, j) +
                   istrx * istry * met(4, i, j, k) * tau(3, i, j));
      forcing(2, i, j) =
          sqjac * (istry * met(2, i, j, k) * tau(2, i, j) +
                   istrx * met(3, i, j, k) * tau(4, i, j) +
                   istrx * istry * met(4, i, j, k) * tau(5, i, j));
      forcing(3, i, j) =
          sqjac * (istry * met(2, i, j, k) * tau(3, i, j) +
                   istrx * met(3, i, j, k) * tau(5, i, j) +
                   istrx * istry * met(4, i, j, k) * tau(6, i, j));
    }
  }
#undef met
#undef jac
#undef forcing
#undef tau
#undef strx
#undef stry
}

//-----------------------------------------------------------------------
void EW::subsurfforcingsg_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int k,
    float_sw4* __restrict__ a_met, float_sw4* __restrict__ a_jac,
    float_sw4* __restrict__ a_tau, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry, float_sw4* __restrict__ a_forcing) {
  SW4_MARK_FUNCTION;
  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = ni * (jlast - jfirst + 1) * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int basef = -(ifirst + ni * jfirst);
  const int base3 = base - nijk;
  // const int basef3= basef-nij;
  const int nic3 = 3 * ni;
  // const int nic6  = 6*ni;

#define met(c, i, j, k) a_met[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define jac(i, j, k) a_jac[base + (i) + ni * (j) + nij * (k)]
#define forcing(c, i, j) a_forcing[3 * basef - 1 + (c) + 3 * (i) + nic3 * (j)]
  //#define tau(c,i,j)         a_tau[6*basef-1+(c)+6*(i)+nic6*(j)]
#define tau(c, i, j) a_tau[basef - nij + (i) + ni * (j) + nij * (c)]
#define strx(i) a_strx[i - ifirst]
#define stry(j) a_stry[j - jfirst]

#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++) {
    float_sw4 istry = 1 / stry(j);
#pragma ivdep
    //#pragma simd
    for (int i = ifirst; i <= ilast; i++) {
      float_sw4 istrx = 1 / strx(i);
      float_sw4 sqjac = sqrt(jac(i, j, k));
      forcing(1, i, j) -=
          sqjac * (istry * met(2, i, j, k) * tau(1, i, j) +
                   istrx * met(3, i, j, k) * tau(2, i, j) +
                   istrx * istry * met(4, i, j, k) * tau(3, i, j));
      forcing(2, i, j) -=
          sqjac * (istry * met(2, i, j, k) * tau(2, i, j) +
                   istrx * met(3, i, j, k) * tau(4, i, j) +
                   istrx * istry * met(4, i, j, k) * tau(5, i, j));
      forcing(3, i, j) -=
          sqjac * (istry * met(2, i, j, k) * tau(3, i, j) +
                   istrx * met(3, i, j, k) * tau(5, i, j) +
                   istrx * istry * met(4, i, j, k) * tau(6, i, j));
    }
  }
#undef met
#undef jac
#undef forcing
#undef tau
#undef strx
#undef stry
}
