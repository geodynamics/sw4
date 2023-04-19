#include "sw4.h"
#include "caliper.h"
//-----------------------------------------------------------------------
void add_pseudohessian_terms2(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nk,
    int ifirstact, int ilastact, int jfirstact, int jlastact, int kfirstact,
    int klastact, float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_rho,
    float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, float_sw4 h,
    float_sw4 dt, int onesided[6], int varcase, float_sw4* __restrict__ a_bope,
    float_sw4* __restrict__ a_acof, float_sw4* a_ghcof,
    float_sw4* __restrict__ a_ph) {
  SW4_MARK_FUNCTION;
  const float_sw4 ih2 = 1.0 / (h * h);
  const float_sw4 idt2 = 1 / (dt * dt);
  const float_sw4 d4a = 2.0 / 3;
  const float_sw4 d4b = -1.0 / 12;
  const int ni = ilast - ifirst + 1;
  const int nij = ni * (jlast - jfirst + 1);
  const int nijk = nij * (klast - kfirst + 1);
  const int base = -(ifirst + ni * jfirst + nij * kfirst);
  const int base3 = base - nijk;
  const float_sw4 o6 = 1.0 / 6.0;
  const float_sw4 o24 = 1.0 / 24.0;

#define rho(i, j, k) a_rho[base + (i) + ni * (j) + nij * (k)]
#define mu(i, j, k) a_mu[base + (i) + ni * (j) + nij * (k)]
#define lambda(i, j, k) a_lambda[base + (i) + ni * (j) + nij * (k)]
#define um(c, i, j, k) a_um[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define u(c, i, j, k) a_u[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define up(c, i, j, k) a_up[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define ph(c, i, j, k) a_ph[base3 + (i) + ni * (j) + nij * (k) + nijk * (c)]
#define bope(k, q) a_bope[k - 1 + 6 * (q - 1)]
#define acof(k, p, q) a_acof[k - 1 + 6 * (p - 1) + 48 * (q - 1)]
#define ghcof(k) a_ghcof[k - 1]
#define SQR(x) ((x) * (x))

#define srcla(c, i, j, k) \
  srcla_[(i + 2) + 5 * (j + 2) + 25 * (k - 1) + (c - 1) * 250]
#define srcmu(c, i, j, k) \
  srcmu_[(i + 2) + 5 * (j + 2) + 25 * (k - 1) + (c - 1) * 250]

  float_sw4* srcmu_ = new float_sw4[750];
  float_sw4* srcla_ = new float_sw4[750];

  int nb = 6, wb = 8;

  int kstart = kfirstact;
  if (onesided[4] == 1 && kfirstact <= nb) {
    // SBP boundary operators
    kstart = nb + 1;
    for (int n = kfirstact; n <= nb; n++)

      // #pragma omp parallel for
      for (int m = jfirstact; m <= jlastact; m++) /* #pragma ivdep */
        for (int l = ifirstact; l <= ilastact; l++) {
          //	 srcla_ = srcpla[s];
          //	 srcmu_ = srcpmu[s];
          for (int c = 0; c < 3 * 5 * 5 * 10; c++) srcla_[c] = srcmu_[c] = 0;

          float_sw4 du, dv, dw;
          //         float_sw4 srcla1, srcla2, srcla3, srcmu1=0, srcmu2=0,
          //         srcmu3=0;

          // d L^u/d(lambda):
          srcla(1, +2, 0, n) = -0.125 * u(1, l, m, n) + o6 * u(1, l + 1, m, n) -
                               o24 * u(1, l + 2, m, n);
          srcla(1, +1, 0, n) = o6 * u(1, l - 1, m, n) + 0.5 * u(1, l, m, n) -
                               5 * o6 * u(1, l + 1, m, n) +
                               o6 * u(1, l + 2, m, n);
          srcla(1, 0, 0, n) = -0.125 * (u(1, l - 2, m, n) + u(1, l + 2, m, n)) +
                              0.5 * (u(1, l - 1, m, n) + u(1, l + 1, m, n)) -
                              0.75 * u(1, l, m, n);
          srcla(1, -1, 0, n) = o6 * u(1, l - 2, m, n) -
                               5 * o6 * u(1, l - 1, m, n) +
                               0.5 * u(1, l, m, n) + o6 * u(1, l + 1, m, n);
          srcla(1, -2, 0, n) = -o24 * u(1, l - 2, m, n) +
                               o6 * u(1, l - 1, m, n) - 0.125 * u(1, l, m, n);

          srcmu(1, +2, 0, n) = 2 * srcla(1, +2, 0, n);
          srcmu(1, +1, 0, n) = 2 * srcla(1, +1, 0, n);
          srcmu(1, 0, 0, n) = 2 * srcla(1, 0, 0, n);
          srcmu(1, -1, 0, n) = 2 * srcla(1, -1, 0, n);
          srcmu(1, -2, 0, n) = 2 * srcla(1, -2, 0, n);

          dv = d4b * (u(2, l, m + 2, n) - u(2, l, m - 2, n)) +
               d4a * (u(2, l, m + 1, n) - u(2, l, m - 1, n));
          srcla(1, +2, 0, n) += -d4b * dv;
          srcla(1, +1, 0, n) += -d4a * dv;
          srcla(1, -1, 0, n) += d4a * dv;
          srcla(1, -2, 0, n) += d4b * dv;

          if (n > nb)
            dw = d4b * (u(3, l, m, n + 2) - u(3, l, m, n - 2)) +
                 d4a * (u(3, l, m, n + 1) - u(3, l, m, n - 1));
          else {
            dw = 0;
            for (int q = 1; q <= wb; q++) dw += bope(n, q) * u(3, l, m, q);
          }

          srcla(1, +2, 0, n) += -d4b * dw;
          srcla(1, +1, 0, n) += -d4a * dw;
          srcla(1, -1, 0, n) += d4a * dw;
          srcla(1, -2, 0, n) += d4b * dw;

          // d L^u/d(mu):
          srcmu(1, 0, +2, n) += -0.125 * u(1, l, m, n) +
                                o6 * u(1, l, m + 1, n) -
                                o24 * u(1, l, m + 2, n);
          srcmu(1, 0, +1, n) += o6 * u(1, l, m - 1, n) + 0.5 * u(1, l, m, n) -
                                5 * o6 * u(1, l, m + 1, n) +
                                o6 * u(1, l, m + 2, n);
          srcmu(1, 0, 0, n) +=
              -0.125 * (u(1, l, m - 2, n) + u(1, l, m + 2, n)) +
              0.5 * (u(1, l, m - 1, n) + u(1, l, m + 1, n)) -
              0.75 * u(1, l, m, n);
          srcmu(1, 0, -1, n) += o6 * u(1, l, m - 2, n) -
                                5 * o6 * u(1, l, m - 1, n) +
                                0.5 * u(1, l, m, n) + o6 * u(1, l, m + 1, n);
          srcmu(1, 0, -2, n) += -o24 * u(1, l, m - 2, n) +
                                o6 * u(1, l, m - 1, n) - 0.125 * u(1, l, m, n);

          float_sw4 ddu;
          for (int k = 1; k <= nb; k++) {
            ddu = 0;
            for (int q = 1; q <= wb; q++) ddu += acof(k, q, n) * u(1, l, m, q);
            srcmu(1, 0, 0, k) += ddu;
          }
          if (n == 1) srcmu(1, 0, 0, n) += ghcof(1) * u(1, l, m, 0);

          if (n > nb - 2)
            srcmu(1, 0, 0, n + 2) += -0.125 * u(1, l, m, n) +
                                     o6 * u(1, l, m, n + 1) -
                                     o24 * u(1, l, m, n + 2);
          if (n > nb - 1)
            srcmu(1, 0, 0, n + 1) +=
                o6 * u(1, l, m, n - 1) + 0.5 * u(1, l, m, n) -
                5 * o6 * u(1, l, m, n + 1) + o6 * u(1, l, m, n + 2);
          if (n > nb)
            srcmu(1, 0, 0, n) +=
                -0.125 * (u(1, l, m, n - 2) + u(1, l, m, n + 2)) +
                0.5 * (u(1, l, m, n - 1) + u(1, l, m, n + 1)) -
                0.75 * u(1, l, m, n);
          if (n > nb + 1)
            srcmu(1, 0, 0, n - 1) +=
                o6 * u(1, l, m, n - 2) - 5 * o6 * u(1, l, m, n - 1) +
                0.5 * u(1, l, m, n) + o6 * u(1, l, m, n + 1);

          dv = d4b * (u(2, l + 2, m, n) - u(2, l - 2, m, n)) +
               d4a * (u(2, l + 1, m, n) - u(2, l - 1, m, n));
          srcmu(1, 0, +2, n) += -d4b * dv;
          srcmu(1, 0, +1, n) += -d4a * dv;
          srcmu(1, 0, -1, n) += d4a * dv;
          srcmu(1, 0, -2, n) += d4b * dv;

          dw = d4b * (u(3, l + 2, m, n) - u(3, l - 2, m, n)) +
               d4a * (u(3, l + 1, m, n) - u(3, l - 1, m, n));
          for (int k = 1; k <= nb; k++) srcmu(1, 0, 0, k) += bope(k, n) * dw;
          if (n > nb - 2) srcmu(1, 0, 0, n + 2) += -d4b * dw;
          if (n > nb - 1) srcmu(1, 0, 0, n + 1) += -d4a * dw;
          if (n > nb + 1) srcmu(1, 0, 0, n - 1) += d4a * dw;

          // Transform and square
          if (varcase >= 2) {
            for (int k = 1; k <= nb + 2; k++)
              srcmu(1, 0, 0, k) = srcmu(1, 0, 0, k) - 2 * srcla(1, 0, 0, k);
            srcmu(1, -2, 0, n) = srcmu(1, -2, 0, n) - 2 * srcla(1, -2, 0, n);
            srcmu(1, -1, 0, n) = srcmu(1, -1, 0, n) - 2 * srcla(1, -1, 0, n);
            srcmu(1, 1, 0, n) = srcmu(1, 1, 0, n) - 2 * srcla(1, 1, 0, n);
            srcmu(1, 2, 0, n) = srcmu(1, 2, 0, n) - 2 * srcla(1, 2, 0, n);
            srcmu(1, 0, -2, n) = srcmu(1, 0, -2, n) - 2 * srcla(1, 0, -2, n);
            srcmu(1, 0, -1, n) = srcmu(1, 0, -1, n) - 2 * srcla(1, 0, -1, n);
            srcmu(1, 0, 1, n) = srcmu(1, 0, 1, n) - 2 * srcla(1, 0, 1, n);
            srcmu(1, 0, 2, n) = srcmu(1, 0, 2, n) - 2 * srcla(1, 0, 2, n);
            //            srcrho = srcrho + mu*srcmu
          }
          float_sw4 phtermmu = 0;
          for (int k = 1; k <= nb + 2; k++)
            phtermmu += srcmu(1, 0, 0, k) * srcmu(1, 0, 0, k);
          phtermmu += srcmu(1, -2, 0, n) * srcmu(1, -2, 0, n) +
                      srcmu(1, -1, 0, n) * srcmu(1, -1, 0, n) +
                      srcmu(1, 1, 0, n) * srcmu(1, 1, 0, n) +
                      srcmu(1, 2, 0, n) * srcmu(1, 2, 0, n) +
                      srcmu(1, 0, -2, n) * srcmu(1, 0, -2, n) +
                      srcmu(1, 0, -1, n) * srcmu(1, 0, -1, n) +
                      srcmu(1, 0, 1, n) * srcmu(1, 0, 1, n) +
                      srcmu(1, 0, 2, n) * srcmu(1, 0, 2, n);
          float_sw4 phtermla = 0;
          for (int k = 1; k <= nb + 2; k++)
            phtermla += srcla(1, 0, 0, k) * srcla(1, 0, 0, k);
          phtermla += srcla(1, -2, 0, n) * srcla(1, -2, 0, n) +
                      srcla(1, -1, 0, n) * srcla(1, -1, 0, n) +
                      srcla(1, 1, 0, n) * srcla(1, 1, 0, n) +
                      srcla(1, 2, 0, n) * srcla(1, 2, 0, n) +
                      srcla(1, 0, -2, n) * srcla(1, 0, -2, n) +
                      srcla(1, 0, -1, n) * srcla(1, 0, -1, n) +
                      srcla(1, 0, 1, n) * srcla(1, 0, 1, n) +
                      srcla(1, 0, 2, n) * srcla(1, 0, 2, n);
          float_sw4 phtermrho =
              SQR(idt2 * (up(1, l, m, n) - 2 * u(1, l, m, n) + um(1, l, m, n)));
          if (varcase == 1) {
            ph(1, l, m, n) += phtermrho;
            ph(2, l, m, n) += phtermmu;
            ph(3, l, m, n) += phtermla;
          } else if (varcase == 2) {
            float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
            float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
            ph(1, l, m, n) +=
                phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else if (varcase == 3) {
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else {
            float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
            ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
            //            ph(3,l,m,n) +=
            //            4*rho(l,m,n)*mu(l,m,n)*(irat*irat*phtermmu+phtermla);
          }

          // d L^v/d(lambda):
          srcla(2, 0, +2, n) = -0.125 * u(2, l, m, n) + o6 * u(2, l, m + 1, n) -
                               o24 * u(2, l, m + 2, n);
          srcla(2, 0, +1, n) = o6 * u(2, l, m - 1, n) + 0.5 * u(2, l, m, n) -
                               5 * o6 * u(2, l, m + 1, n) +
                               o6 * u(2, l, m + 2, n);
          srcla(2, 0, 0, n) = -0.125 * (u(2, l, m - 2, n) + u(2, l, m + 2, n)) +
                              0.5 * (u(2, l, m - 1, n) + u(2, l, m + 1, n)) -
                              0.75 * u(2, l, m, n);
          srcla(2, 0, -1, n) = o6 * u(2, l, m - 2, n) -
                               5 * o6 * u(2, l, m - 1, n) +
                               0.5 * u(2, l, m, n) + o6 * u(2, l, m + 1, n);
          srcla(2, 0, -2, n) = -o24 * u(2, l, m - 2, n) +
                               o6 * u(2, l, m - 1, n) - 0.125 * u(2, l, m, n);

          srcmu(2, 0, +2, n) = 2 * srcla(2, 0, +2, n);
          srcmu(2, 0, +1, n) = 2 * srcla(2, 0, +1, n);
          srcmu(2, 0, 0, n) = 2 * srcla(2, 0, 0, n);
          srcmu(2, 0, -1, n) = 2 * srcla(2, 0, -1, n);
          srcmu(2, 0, -2, n) = 2 * srcla(2, 0, -2, n);

          du = d4b * (u(1, l + 2, m, n) - u(1, l - 2, m, n)) +
               d4a * (u(1, l + 1, m, n) - u(1, l - 1, m, n));
          srcla(2, 0, +2, n) += -d4b * du;
          srcla(2, 0, +1, n) += -d4a * du;
          srcla(2, 0, -1, n) += d4a * du;
          srcla(2, 0, -2, n) += d4b * du;

          if (n > nb)
            dw = d4b * (u(3, l, m, n + 2) - u(3, l, m, n - 2)) +
                 d4a * (u(3, l, m, n + 1) - u(3, l, m, n - 1));
          else {
            dw = 0;
            for (int q = 1; q <= wb; q++) dw += bope(n, q) * u(3, l, m, q);
          }
          srcla(2, 0, +2, n) += -d4b * dw;
          srcla(2, 0, +1, n) += -d4a * dw;
          srcla(2, 0, -1, n) += d4a * dw;
          srcla(2, 0, -2, n) += d4b * dw;

          // dL^v/d(mu)
          srcmu(2, +2, 0, n) += -0.125 * u(2, l, m, n) +
                                o6 * u(2, l + 1, m, n) -
                                o24 * u(2, l + 2, m, n);
          srcmu(2, +1, 0, n) += o6 * u(2, l - 1, m, n) + 0.5 * u(2, l, m, n) -
                                5 * o6 * u(2, l + 1, m, n) +
                                o6 * u(2, l + 2, m, n);
          srcmu(2, 0, 0, n) +=
              -0.125 * (u(2, l - 2, m, n) + u(2, l + 2, m, n)) +
              0.5 * (u(2, l - 1, m, n) + u(2, l + 1, m, n)) -
              0.75 * u(2, l, m, n);
          srcmu(2, -1, 0, n) += o6 * u(2, l - 2, m, n) -
                                5 * o6 * u(2, l - 1, m, n) +
                                0.5 * u(2, l, m, n) + o6 * u(2, l + 1, m, n);
          srcmu(2, -2, 0, n) += -o24 * u(2, l - 2, m, n) +
                                o6 * u(2, l - 1, m, n) - 0.125 * u(2, l, m, n);

          for (int k = 1; k <= nb; k++) {
            ddu = 0;
            for (int q = 1; q <= wb; q++) ddu += acof(k, q, n) * u(2, l, m, q);
            srcmu(2, 0, 0, k) += ddu;
          }
          if (n == 1) srcmu(2, 0, 0, n) += ghcof(1) * u(2, l, m, 0);

          if (n > nb - 2)
            srcmu(2, 0, 0, n + 2) += -0.125 * u(2, l, m, n) +
                                     o6 * u(2, l, m, n + 1) -
                                     o24 * u(2, l, m, n + 2);
          if (n > nb - 1)
            srcmu(2, 0, 0, n + 1) +=
                o6 * u(2, l, m, n - 1) + 0.5 * u(2, l, m, n) -
                5 * o6 * u(2, l, m, n + 1) + o6 * u(2, l, m, n + 2);
          if (n > nb)
            srcmu(2, 0, 0, n) +=
                -0.125 * (u(2, l, m, n - 2) + u(2, l, m, n + 2)) +
                0.5 * (u(2, l, m, n - 1) + u(2, l, m, n + 1)) -
                0.75 * u(2, l, m, n);
          if (n > nb + 1)
            srcmu(2, 0, 0, n - 1) +=
                o6 * u(2, l, m, n - 2) - 5 * o6 * u(2, l, m, n - 1) +
                0.5 * u(2, l, m, n) + o6 * u(2, l, m, n + 1);

          du = d4b * (u(1, l, m + 2, n) - u(1, l, m - 2, n)) +
               d4a * (u(1, l, m + 1, n) - u(1, l, m - 1, n));
          srcmu(2, +2, 0, n) += -d4b * du;
          srcmu(2, +1, 0, n) += -d4a * du;
          srcmu(2, -1, 0, n) += d4a * du;
          srcmu(2, -2, 0, n) += d4b * du;

          dw = d4b * (u(3, l, m + 2, n) - u(3, l, m - 2, n)) +
               d4a * (u(3, l, m + 1, n) - u(3, l, m - 1, n));
          for (int k = 1; k <= nb; k++) srcmu(2, 0, 0, k) += bope(k, n) * dw;
          if (n > nb - 2) srcmu(2, 0, 0, n + 2) += -d4b * dw;
          if (n > nb - 1) srcmu(2, 0, 0, n + 1) += -d4a * dw;
          if (n > nb + 1) srcmu(2, 0, 0, n - 1) += d4a * dw;

          // Transform and square
          if (varcase >= 2) {
            for (int k = 1; k <= nb + 2; k++)
              srcmu(2, 0, 0, k) = srcmu(2, 0, 0, k) - 2 * srcla(2, 0, 0, k);
            srcmu(2, -2, 0, n) = srcmu(2, -2, 0, n) - 2 * srcla(2, -2, 0, n);
            srcmu(2, -1, 0, n) = srcmu(2, -1, 0, n) - 2 * srcla(2, -1, 0, n);
            srcmu(2, 1, 0, n) = srcmu(2, 1, 0, n) - 2 * srcla(2, 1, 0, n);
            srcmu(2, 2, 0, n) = srcmu(2, 2, 0, n) - 2 * srcla(2, 2, 0, n);
            srcmu(2, 0, -2, n) = srcmu(2, 0, -2, n) - 2 * srcla(2, 0, -2, n);
            srcmu(2, 0, -1, n) = srcmu(2, 0, -1, n) - 2 * srcla(2, 0, -1, n);
            srcmu(2, 0, 1, n) = srcmu(2, 0, 1, n) - 2 * srcla(2, 0, 1, n);
            srcmu(2, 0, 2, n) = srcmu(2, 0, 2, n) - 2 * srcla(2, 0, 2, n);
          }
          phtermmu = 0;
          for (int k = 1; k <= nb + 2; k++)
            phtermmu += srcmu(2, 0, 0, k) * srcmu(2, 0, 0, k);
          phtermmu += srcmu(2, -2, 0, n) * srcmu(2, -2, 0, n) +
                      srcmu(2, -1, 0, n) * srcmu(2, -1, 0, n) +
                      srcmu(2, 1, 0, n) * srcmu(2, 1, 0, n) +
                      srcmu(2, 2, 0, n) * srcmu(2, 2, 0, n) +
                      srcmu(2, 0, -2, n) * srcmu(2, 0, -2, n) +
                      srcmu(2, 0, -1, n) * srcmu(2, 0, -1, n) +
                      srcmu(2, 0, 1, n) * srcmu(2, 0, 1, n) +
                      srcmu(2, 0, 2, n) * srcmu(2, 0, 2, n);
          phtermla = 0;
          for (int k = 1; k <= nb + 2; k++)
            phtermla += srcla(2, 0, 0, k) * srcla(2, 0, 0, k);
          phtermla += srcla(2, -2, 0, n) * srcla(2, -2, 0, n) +
                      srcla(2, -1, 0, n) * srcla(2, -1, 0, n) +
                      srcla(2, 1, 0, n) * srcla(2, 1, 0, n) +
                      srcla(2, 2, 0, n) * srcla(2, 2, 0, n) +
                      srcla(2, 0, -2, n) * srcla(2, 0, -2, n) +
                      srcla(2, 0, -1, n) * srcla(2, 0, -1, n) +
                      srcla(2, 0, 1, n) * srcla(2, 0, 1, n) +
                      srcla(2, 0, 2, n) * srcla(2, 0, 2, n);
          phtermrho =
              SQR(idt2 * (up(2, l, m, n) - 2 * u(2, l, m, n) + um(2, l, m, n)));
          if (varcase == 1) {
            ph(1, l, m, n) += phtermrho;
            ph(2, l, m, n) += phtermmu;
            ph(3, l, m, n) += phtermla;
          } else if (varcase == 2) {
            float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
            float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
            ph(1, l, m, n) +=
                phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else if (varcase == 3) {
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else {
            float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
            //            ph(3,l,m,n) +=
            //            4*rho(l,m,n)*mu(l,m,n)*(irat*irat*phtermmu+phtermla);
            ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          }

          // dL^w/d(lambda)
          for (int k = 1; k <= nb; k++) {
            ddu = 0;
            for (int q = 1; q <= wb; q++) ddu += acof(k, q, n) * u(3, l, m, q);
            srcla(3, 0, 0, k) += ddu;
          }
          if (n == 1) srcla(3, 0, 0, n) += ghcof(1) * u(3, l, m, 0);

          if (n > nb - 2)
            srcla(3, 0, 0, n + 2) += -0.125 * u(3, l, m, n) +
                                     o6 * u(3, l, m, n + 1) -
                                     o24 * u(3, l, m, n + 2);
          if (n > nb - 1)
            srcla(3, 0, 0, n + 1) +=
                o6 * u(3, l, m, n - 1) + 0.5 * u(3, l, m, n) -
                5 * o6 * u(3, l, m, n + 1) + o6 * u(3, l, m, n + 2);
          if (n > nb)
            srcla(3, 0, 0, n) +=
                -0.125 * (u(3, l, m, n - 2) + u(3, l, m, n + 2)) +
                0.5 * (u(3, l, m, n - 1) + u(3, l, m, n + 1)) -
                0.75 * u(3, l, m, n);
          if (n > nb + 1)
            srcla(3, 0, 0, n - 1) +=
                o6 * u(3, l, m, n - 2) - 5 * o6 * u(3, l, m, n - 1) +
                0.5 * u(3, l, m, n) + o6 * u(3, l, m, n + 1);

          //      srcla(3,0,0,+2) =
          //      -0.125*u(3,l,m,n)+o6*u(3,l,m,n+1)-o24*u(3,l,m,n+2);
          //      srcla(3,0,0,+1) =
          //      o6*u(3,l,m,n-1)+0.5*u(3,l,m,n)-5*o6*u(3,l,m,n+1)+o6*u(3,l,m,n+2);
          //      srcla(3,0,0, 0) =
          //      -0.125*(u(3,l,m,n-2)+u(3,l,m,n+2))+0.5*(u(3,l,m,n-1)+u(3,l,m,n+1))-0.75*u(3,l,m,n);
          //      srcla(3,0,0,-1) =
          //      o6*u(3,l,m,n-2)-5*o6*u(3,l,m,n-1)+0.5*u(3,l,m,n)+o6*u(3,l,m,n+1);
          //      srcla(3,0,0,-2) =
          //      -o24*u(3,l,m,n-2)+o6*u(3,l,m,n-1)-0.125*u(3,l,m,n);

          for (int q = 1; q <= wb + 2; q++) {
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
          }

          dv = d4b * (u(2, l, m + 2, n) - u(2, l, m - 2, n)) +
               d4a * (u(2, l, m + 1, n) - u(2, l, m - 1, n));
          du = d4b * (u(1, l + 2, m, n) - u(1, l - 2, m, n)) +
               d4a * (u(1, l + 1, m, n) - u(1, l - 1, m, n));
          for (int k = 1; k <= nb; k++)
            srcla(3, 0, 0, k) += bope(k, n) * (du + dv);
          if (n > nb - 2) srcla(3, 0, 0, n + 2) += -d4b * (du + dv);
          if (n > nb - 1) srcla(3, 0, 0, n + 1) += -d4a * (du + dv);
          if (n > nb + 1) srcla(3, 0, 0, n - 1) += d4a * (du + dv);

          // dL^w/d(mu)
          srcmu(3, +2, 0, n) += -0.125 * u(3, l, m, n) +
                                o6 * u(3, l + 1, m, n) -
                                o24 * u(3, l + 2, m, n);
          srcmu(3, +1, 0, n) += o6 * u(3, l - 1, m, n) + 0.5 * u(3, l, m, n) -
                                5 * o6 * u(3, l + 1, m, n) +
                                o6 * u(3, l + 2, m, n);
          srcmu(3, 0, 0, n) +=
              -0.125 * (u(3, l - 2, m, n) + u(3, l + 2, m, n)) +
              0.5 * (u(3, l - 1, m, n) + u(3, l + 1, m, n)) -
              0.75 * u(3, l, m, n);
          srcmu(3, -1, 0, n) += o6 * u(3, l - 2, m, n) -
                                5 * o6 * u(3, l - 1, m, n) +
                                0.5 * u(3, l, m, n) + o6 * u(3, l + 1, m, n);
          srcmu(3, -2, 0, n) += -o24 * u(3, l - 2, m, n) +
                                o6 * u(3, l - 1, m, n) - 0.125 * u(3, l, m, n);

          srcmu(3, 0, +2, n) += -0.125 * u(3, l, m, n) +
                                o6 * u(3, l, m + 1, n) -
                                o24 * u(3, l, m + 2, n);
          srcmu(3, 0, +1, n) += o6 * u(3, l, m - 1, n) + 0.5 * u(3, l, m, n) -
                                5 * o6 * u(3, l, m + 1, n) +
                                o6 * u(3, l, m + 2, n);
          srcmu(3, 0, 0, n) +=
              -0.125 * (u(3, l, m - 2, n) + u(3, l, m + 2, n)) +
              0.5 * (u(3, l, m - 1, n) + u(3, l, m + 1, n)) -
              0.75 * u(3, l, m, n);
          srcmu(3, 0, -1, n) += o6 * u(3, l, m - 2, n) -
                                5 * o6 * u(3, l, m - 1, n) +
                                0.5 * u(3, l, m, n) + o6 * u(3, l, m + 1, n);
          srcmu(3, 0, -2, n) += -o24 * u(3, l, m - 2, n) +
                                o6 * u(3, l, m - 1, n) - 0.125 * u(3, l, m, n);

          if (n > nb) {
            du = d4b * (u(1, l, m, n + 2) - u(1, l, m, n - 2)) +
                 d4a * (u(1, l, m, n + 1) - u(1, l, m, n - 1));
            dv = d4b * (u(2, l, m, n + 2) - u(2, l, m, n - 2)) +
                 d4a * (u(2, l, m, n + 1) - u(2, l, m, n - 1));
          } else {
            du = dv = 0;
            for (int q = 1; q <= wb; q++) {
              du += bope(n, q) * u(1, l, m, q);
              dv += bope(n, q) * u(2, l, m, q);
            }
          }

          srcmu(3, +2, 0, n) += -d4b * du;
          srcmu(3, +1, 0, n) += -d4a * du;
          srcmu(3, -1, 0, n) += d4a * du;
          srcmu(3, -2, 0, n) += d4b * du;

          srcmu(3, 0, +2, n) += -d4b * dv;
          srcmu(3, 0, +1, n) += -d4a * dv;
          srcmu(3, 0, -1, n) += d4a * dv;
          srcmu(3, 0, -2, n) += d4b * dv;

          // Square and transform
          if (varcase >= 2) {
            for (int k = 1; k <= nb + 2; k++)
              srcmu(3, 0, 0, k) = srcmu(3, 0, 0, k) - 2 * srcla(3, 0, 0, k);
            srcmu(3, -2, 0, n) = srcmu(3, -2, 0, n) - 2 * srcla(3, -2, 0, n);
            srcmu(3, -1, 0, n) = srcmu(3, -1, 0, n) - 2 * srcla(3, -1, 0, n);
            srcmu(3, 1, 0, n) = srcmu(3, 1, 0, n) - 2 * srcla(3, 1, 0, n);
            srcmu(3, 2, 0, n) = srcmu(3, 2, 0, n) - 2 * srcla(3, 2, 0, n);
            srcmu(3, 0, -2, n) = srcmu(3, 0, -2, n) - 2 * srcla(3, 0, -2, n);
            srcmu(3, 0, -1, n) = srcmu(3, 0, -1, n) - 2 * srcla(3, 0, -1, n);
            srcmu(3, 0, 1, n) = srcmu(3, 0, 1, n) - 2 * srcla(3, 0, 1, n);
            srcmu(3, 0, 2, n) = srcmu(3, 0, 2, n) - 2 * srcla(3, 0, 2, n);
          }
          phtermmu = 0;
          for (int k = 1; k <= nb + 2; k++)
            phtermmu += srcmu(3, 0, 0, k) * srcmu(3, 0, 0, k);
          phtermmu += srcmu(3, -2, 0, n) * srcmu(3, -2, 0, n) +
                      srcmu(3, -1, 0, n) * srcmu(3, -1, 0, n) +
                      srcmu(3, 1, 0, n) * srcmu(3, 1, 0, n) +
                      srcmu(3, 2, 0, n) * srcmu(3, 2, 0, n) +
                      srcmu(3, 0, -2, n) * srcmu(3, 0, -2, n) +
                      srcmu(3, 0, -1, n) * srcmu(3, 0, -1, n) +
                      srcmu(3, 0, 1, n) * srcmu(3, 0, 1, n) +
                      srcmu(3, 0, 2, n) * srcmu(3, 0, 2, n);
          phtermla = 0;
          for (int k = 1; k <= nb + 2; k++)
            phtermla += srcla(3, 0, 0, k) * srcla(3, 0, 0, k);
          phtermla += srcla(3, -2, 0, n) * srcla(3, -2, 0, n) +
                      srcla(3, -1, 0, n) * srcla(3, -1, 0, n) +
                      srcla(3, 1, 0, n) * srcla(3, 1, 0, n) +
                      srcla(3, 2, 0, n) * srcla(3, 2, 0, n) +
                      srcla(3, 0, -2, n) * srcla(3, 0, -2, n) +
                      srcla(3, 0, -1, n) * srcla(3, 0, -1, n) +
                      srcla(3, 0, 1, n) * srcla(3, 0, 1, n) +
                      srcla(3, 0, 2, n) * srcla(3, 0, 2, n);
          phtermrho =
              SQR(idt2 * (up(3, l, m, n) - 2 * u(3, l, m, n) + um(3, l, m, n)));

          if (varcase == 1) {
            ph(1, l, m, n) += phtermrho;
            ph(2, l, m, n) += phtermmu;
            ph(3, l, m, n) += phtermla;
          } else if (varcase == 2) {
            float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
            float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
            ph(1, l, m, n) +=
                phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else if (varcase == 3) {
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else {
            float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
            ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          }
        }
  }

#undef srcla
#undef srcmu

#define srcla(c, i, j, k) \
  srcla_[(i + 2) + 5 * (j + 2) + 25 * (k - nk + 9) + (c - 1) * 250]
#define srcmu(c, i, j, k) \
  srcmu_[(i + 2) + 5 * (j + 2) + 25 * (k - nk + 9) + (c - 1) * 250]

  int kend = klastact;
  if (onesided[5] == 1 && klastact >= nk - nb + 1) {
    // SBP boundary operators at lower bndry
    kend = nk - nb;
    for (int n = nk - nb + 1; n <= klastact; n++)
      for (int m = jfirstact; m <= jlastact; m++)
        for (int l = ifirstact; l <= ilastact; l++) {
          for (int c = 0; c < 3 * 5 * 5 * 10; c++) srcla_[c] = srcmu_[c] = 0;

          float_sw4 du, dv, dw;

          // d L^u/d(lambda):
          srcla(1, +2, 0, n) = -0.125 * u(1, l, m, n) + o6 * u(1, l + 1, m, n) -
                               o24 * u(1, l + 2, m, n);
          srcla(1, +1, 0, n) = o6 * u(1, l - 1, m, n) + 0.5 * u(1, l, m, n) -
                               5 * o6 * u(1, l + 1, m, n) +
                               o6 * u(1, l + 2, m, n);
          srcla(1, 0, 0, n) = -0.125 * (u(1, l - 2, m, n) + u(1, l + 2, m, n)) +
                              0.5 * (u(1, l - 1, m, n) + u(1, l + 1, m, n)) -
                              0.75 * u(1, l, m, n);
          srcla(1, -1, 0, n) = o6 * u(1, l - 2, m, n) -
                               5 * o6 * u(1, l - 1, m, n) +
                               0.5 * u(1, l, m, n) + o6 * u(1, l + 1, m, n);
          srcla(1, -2, 0, n) = -o24 * u(1, l - 2, m, n) +
                               o6 * u(1, l - 1, m, n) - 0.125 * u(1, l, m, n);

          srcmu(1, +2, 0, n) = 2 * srcla(1, +2, 0, n);
          srcmu(1, +1, 0, n) = 2 * srcla(1, +1, 0, n);
          srcmu(1, 0, 0, n) = 2 * srcla(1, 0, 0, n);
          srcmu(1, -1, 0, n) = 2 * srcla(1, -1, 0, n);
          srcmu(1, -2, 0, n) = 2 * srcla(1, -2, 0, n);

          dv = d4b * (u(2, l, m + 2, n) - u(2, l, m - 2, n)) +
               d4a * (u(2, l, m + 1, n) - u(2, l, m - 1, n));
          srcla(1, +2, 0, n) += -d4b * dv;
          srcla(1, +1, 0, n) += -d4a * dv;
          srcla(1, -1, 0, n) += d4a * dv;
          srcla(1, -2, 0, n) += d4b * dv;

          if (n < nk - nb + 1)
            dw = d4b * (u(3, l, m, n + 2) - u(3, l, m, n - 2)) +
                 d4a * (u(3, l, m, n + 1) - u(3, l, m, n - 1));
          else {
            dw = 0;
            for (int q = 1; q <= wb; q++)
              dw -= bope(nk - n + 1, q) * u(3, l, m, nk - q + 1);
          }

          srcla(1, +2, 0, n) += -d4b * dw;
          srcla(1, +1, 0, n) += -d4a * dw;
          srcla(1, -1, 0, n) += d4a * dw;
          srcla(1, -2, 0, n) += d4b * dw;

          // d L^u/d(mu):
          srcmu(1, 0, +2, n) += -0.125 * u(1, l, m, n) +
                                o6 * u(1, l, m + 1, n) -
                                o24 * u(1, l, m + 2, n);
          srcmu(1, 0, +1, n) += o6 * u(1, l, m - 1, n) + 0.5 * u(1, l, m, n) -
                                5 * o6 * u(1, l, m + 1, n) +
                                o6 * u(1, l, m + 2, n);
          srcmu(1, 0, 0, n) +=
              -0.125 * (u(1, l, m - 2, n) + u(1, l, m + 2, n)) +
              0.5 * (u(1, l, m - 1, n) + u(1, l, m + 1, n)) -
              0.75 * u(1, l, m, n);
          srcmu(1, 0, -1, n) += o6 * u(1, l, m - 2, n) -
                                5 * o6 * u(1, l, m - 1, n) +
                                0.5 * u(1, l, m, n) + o6 * u(1, l, m + 1, n);
          srcmu(1, 0, -2, n) += -o24 * u(1, l, m - 2, n) +
                                o6 * u(1, l, m - 1, n) - 0.125 * u(1, l, m, n);

          float_sw4 ddu;
          for (int k = 1; k <= nb; k++) {
            ddu = 0;
            for (int q = 1; q <= wb; q++)
              ddu += acof(k, q, nk - n + 1) * u(1, l, m, nk - q + 1);
            srcmu(1, 0, 0, nk - k + 1) += ddu;
          }
          if (n == nk) srcmu(1, 0, 0, n) += ghcof(1) * u(1, l, m, nk + 1);

          //         if( n > nb-2 )
          if (n < nk - nb + 3)
            srcmu(1, 0, 0, n - 2) += -0.125 * u(1, l, m, n) +
                                     o6 * u(1, l, m, n - 1) -
                                     o24 * u(1, l, m, n - 2);
          //	 if( n> nb-1 )
          if (n < nk - nb + 2)
            srcmu(1, 0, 0, n - 1) +=
                o6 * u(1, l, m, n - 2) - 5 * o6 * u(1, l, m, n - 1) +
                0.5 * u(1, l, m, n) + o6 * u(1, l, m, n + 1);
          //	 if( n > nb )
          if (n < nk - nb + 1)
            srcmu(1, 0, 0, n) +=
                -0.125 * (u(1, l, m, n - 2) + u(1, l, m, n + 2)) +
                0.5 * (u(1, l, m, n - 1) + u(1, l, m, n + 1)) -
                0.75 * u(1, l, m, n);
          //	 if( n > nb+1 )
          if (n < nk - nb)
            srcmu(1, 0, 0, n + 1) +=
                o6 * u(1, l, m, n - 1) + 0.5 * u(1, l, m, n) -
                5 * o6 * u(1, l, m, n + 1) + o6 * u(1, l, m, n + 2);
          // //	 srcmu(1,0,0,-2) +=
          // -o24*u(1,l,m,n-2)+o6*u(1,l,m,n-1)-0.125*u(1,l,m,n);

          dv = d4b * (u(2, l + 2, m, n) - u(2, l - 2, m, n)) +
               d4a * (u(2, l + 1, m, n) - u(2, l - 1, m, n));
          srcmu(1, 0, +2, n) += -d4b * dv;
          srcmu(1, 0, +1, n) += -d4a * dv;
          srcmu(1, 0, -1, n) += d4a * dv;
          srcmu(1, 0, -2, n) += d4b * dv;

          dw = d4b * (u(3, l + 2, m, n) - u(3, l - 2, m, n)) +
               d4a * (u(3, l + 1, m, n) - u(3, l - 1, m, n));
          for (int k = 1; k <= nb; k++)
            srcmu(1, 0, 0, nk - k + 1) -= bope(k, nk - n + 1) * dw;
          if (n < nk - nb + 3) srcmu(1, 0, 0, n - 2) += d4b * dw;
          if (n < nk - nb + 2) srcmu(1, 0, 0, n - 1) += d4a * dw;
          if (n < nk - nb) srcmu(1, 0, 0, n + 1) += -d4a * dw;

          // Transform and square
          if (varcase >= 2) {
            for (int k = nk - nb - 1; k <= nk; k++)
              srcmu(1, 0, 0, k) = srcmu(1, 0, 0, k) - 2 * srcla(1, 0, 0, k);
            srcmu(1, -2, 0, n) = srcmu(1, -2, 0, n) - 2 * srcla(1, -2, 0, n);
            srcmu(1, -1, 0, n) = srcmu(1, -1, 0, n) - 2 * srcla(1, -1, 0, n);
            srcmu(1, 1, 0, n) = srcmu(1, 1, 0, n) - 2 * srcla(1, 1, 0, n);
            srcmu(1, 2, 0, n) = srcmu(1, 2, 0, n) - 2 * srcla(1, 2, 0, n);
            srcmu(1, 0, -2, n) = srcmu(1, 0, -2, n) - 2 * srcla(1, 0, -2, n);
            srcmu(1, 0, -1, n) = srcmu(1, 0, -1, n) - 2 * srcla(1, 0, -1, n);
            srcmu(1, 0, 1, n) = srcmu(1, 0, 1, n) - 2 * srcla(1, 0, 1, n);
            srcmu(1, 0, 2, n) = srcmu(1, 0, 2, n) - 2 * srcla(1, 0, 2, n);
            //            srcrho = srcrho + mu*srcmu
          }
          float_sw4 phtermmu = 0;
          for (int k = nk - nb - 1; k <= nk; k++)
            phtermmu += srcmu(1, 0, 0, k) * srcmu(1, 0, 0, k);
          phtermmu += srcmu(1, -2, 0, n) * srcmu(1, -2, 0, n) +
                      srcmu(1, -1, 0, n) * srcmu(1, -1, 0, n) +
                      srcmu(1, 1, 0, n) * srcmu(1, 1, 0, n) +
                      srcmu(1, 2, 0, n) * srcmu(1, 2, 0, n) +
                      srcmu(1, 0, -2, n) * srcmu(1, 0, -2, n) +
                      srcmu(1, 0, -1, n) * srcmu(1, 0, -1, n) +
                      srcmu(1, 0, 1, n) * srcmu(1, 0, 1, n) +
                      srcmu(1, 0, 2, n) * srcmu(1, 0, 2, n);
          float_sw4 phtermla = 0;
          for (int k = nk - nb - 1; k <= nk; k++)
            phtermla += srcla(1, 0, 0, k) * srcla(1, 0, 0, k);
          phtermla += srcla(1, -2, 0, n) * srcla(1, -2, 0, n) +
                      srcla(1, -1, 0, n) * srcla(1, -1, 0, n) +
                      srcla(1, 1, 0, n) * srcla(1, 1, 0, n) +
                      srcla(1, 2, 0, n) * srcla(1, 2, 0, n) +
                      srcla(1, 0, -2, n) * srcla(1, 0, -2, n) +
                      srcla(1, 0, -1, n) * srcla(1, 0, -1, n) +
                      srcla(1, 0, 1, n) * srcla(1, 0, 1, n) +
                      srcla(1, 0, 2, n) * srcla(1, 0, 2, n);
          float_sw4 phtermrho =
              SQR(idt2 * (up(1, l, m, n) - 2 * u(1, l, m, n) + um(1, l, m, n)));
          if (varcase == 1) {
            ph(1, l, m, n) += phtermrho;
            ph(2, l, m, n) += phtermmu;
            ph(3, l, m, n) += phtermla;
          } else if (varcase == 2) {
            float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
            float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
            ph(1, l, m, n) +=
                phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else if (varcase == 3) {
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else {
            float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
            ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
            //            ph(3,l,m,n) +=
            //            4*rho(l,m,n)*mu(l,m,n)*(irat*irat*phtermmu+phtermla);
          }

          // d L^v/d(lambda):
          srcla(2, 0, +2, n) = -0.125 * u(2, l, m, n) + o6 * u(2, l, m + 1, n) -
                               o24 * u(2, l, m + 2, n);
          srcla(2, 0, +1, n) = o6 * u(2, l, m - 1, n) + 0.5 * u(2, l, m, n) -
                               5 * o6 * u(2, l, m + 1, n) +
                               o6 * u(2, l, m + 2, n);
          srcla(2, 0, 0, n) = -0.125 * (u(2, l, m - 2, n) + u(2, l, m + 2, n)) +
                              0.5 * (u(2, l, m - 1, n) + u(2, l, m + 1, n)) -
                              0.75 * u(2, l, m, n);
          srcla(2, 0, -1, n) = o6 * u(2, l, m - 2, n) -
                               5 * o6 * u(2, l, m - 1, n) +
                               0.5 * u(2, l, m, n) + o6 * u(2, l, m + 1, n);
          srcla(2, 0, -2, n) = -o24 * u(2, l, m - 2, n) +
                               o6 * u(2, l, m - 1, n) - 0.125 * u(2, l, m, n);

          srcmu(2, 0, +2, n) = 2 * srcla(2, 0, +2, n);
          srcmu(2, 0, +1, n) = 2 * srcla(2, 0, +1, n);
          srcmu(2, 0, 0, n) = 2 * srcla(2, 0, 0, n);
          srcmu(2, 0, -1, n) = 2 * srcla(2, 0, -1, n);
          srcmu(2, 0, -2, n) = 2 * srcla(2, 0, -2, n);

          du = d4b * (u(1, l + 2, m, n) - u(1, l - 2, m, n)) +
               d4a * (u(1, l + 1, m, n) - u(1, l - 1, m, n));
          srcla(2, 0, +2, n) += -d4b * du;
          srcla(2, 0, +1, n) += -d4a * du;
          srcla(2, 0, -1, n) += d4a * du;
          srcla(2, 0, -2, n) += d4b * du;

          if (n < nk - nb + 1)
            dw = d4b * (u(3, l, m, n + 2) - u(3, l, m, n - 2)) +
                 d4a * (u(3, l, m, n + 1) - u(3, l, m, n - 1));
          else {
            dw = 0;
            for (int q = 1; q <= wb; q++)
              dw -= bope(nk - n + 1, q) * u(3, l, m, nk - q + 1);
          }
          srcla(2, 0, +2, n) += -d4b * dw;
          srcla(2, 0, +1, n) += -d4a * dw;
          srcla(2, 0, -1, n) += d4a * dw;
          srcla(2, 0, -2, n) += d4b * dw;

          // dL^v/d(mu)
          srcmu(2, +2, 0, n) += -0.125 * u(2, l, m, n) +
                                o6 * u(2, l + 1, m, n) -
                                o24 * u(2, l + 2, m, n);
          srcmu(2, +1, 0, n) += o6 * u(2, l - 1, m, n) + 0.5 * u(2, l, m, n) -
                                5 * o6 * u(2, l + 1, m, n) +
                                o6 * u(2, l + 2, m, n);
          srcmu(2, 0, 0, n) +=
              -0.125 * (u(2, l - 2, m, n) + u(2, l + 2, m, n)) +
              0.5 * (u(2, l - 1, m, n) + u(2, l + 1, m, n)) -
              0.75 * u(2, l, m, n);
          srcmu(2, -1, 0, n) += o6 * u(2, l - 2, m, n) -
                                5 * o6 * u(2, l - 1, m, n) +
                                0.5 * u(2, l, m, n) + o6 * u(2, l + 1, m, n);
          srcmu(2, -2, 0, n) += -o24 * u(2, l - 2, m, n) +
                                o6 * u(2, l - 1, m, n) - 0.125 * u(2, l, m, n);

          for (int k = 1; k <= nb; k++) {
            ddu = 0;
            for (int q = 1; q <= wb; q++)
              ddu += acof(k, q, nk - n + 1) * u(2, l, m, nk - q + 1);
            srcmu(2, 0, 0, nk - k + 1) += ddu;
          }
          if (n == nk) srcmu(2, 0, 0, n) += ghcof(1) * u(2, l, m, nk + 1);

          if (n < nk - nb + 3)
            srcmu(2, 0, 0, n - 2) += -0.125 * u(2, l, m, n) +
                                     o6 * u(2, l, m, n - 1) -
                                     o24 * u(2, l, m, n - 2);
          if (n < nk - nb + 2)
            srcmu(2, 0, 0, n - 1) +=
                o6 * u(2, l, m, n - 2) + 0.5 * u(2, l, m, n) -
                5 * o6 * u(2, l, m, n - 1) + o6 * u(2, l, m, n + 1);
          if (n < nk - nb + 1)
            srcmu(2, 0, 0, n) +=
                -0.125 * (u(2, l, m, n - 2) + u(2, l, m, n + 2)) +
                0.5 * (u(2, l, m, n - 1) + u(2, l, m, n + 1)) -
                0.75 * u(2, l, m, n);
          if (n < nk - nb)
            srcmu(2, 0, 0, n + 1) +=
                o6 * u(2, l, m, n - 1) - 5 * o6 * u(2, l, m, n + 1) +
                0.5 * u(2, l, m, n) + o6 * u(2, l, m, n + 2);

          du = d4b * (u(1, l, m + 2, n) - u(1, l, m - 2, n)) +
               d4a * (u(1, l, m + 1, n) - u(1, l, m - 1, n));
          srcmu(2, +2, 0, n) += -d4b * du;
          srcmu(2, +1, 0, n) += -d4a * du;
          srcmu(2, -1, 0, n) += d4a * du;
          srcmu(2, -2, 0, n) += d4b * du;

          dw = d4b * (u(3, l, m + 2, n) - u(3, l, m - 2, n)) +
               d4a * (u(3, l, m + 1, n) - u(3, l, m - 1, n));
          for (int k = 1; k <= nb; k++)
            srcmu(2, 0, 0, nk - k + 1) -= bope(k, nk - n + 1) * dw;
          if (n < nk - nb + 3) srcmu(2, 0, 0, n - 2) += d4b * dw;
          if (n < nk - nb + 2) srcmu(2, 0, 0, n - 1) += d4a * dw;
          if (n < nk - nb) srcmu(2, 0, 0, n + 1) += -d4a * dw;

          // Transform and square
          if (varcase >= 2) {
            for (int k = nk - nb - 1; k <= nk; k++)
              srcmu(2, 0, 0, k) = srcmu(2, 0, 0, k) - 2 * srcla(2, 0, 0, k);
            srcmu(2, -2, 0, n) = srcmu(2, -2, 0, n) - 2 * srcla(2, -2, 0, n);
            srcmu(2, -1, 0, n) = srcmu(2, -1, 0, n) - 2 * srcla(2, -1, 0, n);
            srcmu(2, 1, 0, n) = srcmu(2, 1, 0, n) - 2 * srcla(2, 1, 0, n);
            srcmu(2, 2, 0, n) = srcmu(2, 2, 0, n) - 2 * srcla(2, 2, 0, n);
            srcmu(2, 0, -2, n) = srcmu(2, 0, -2, n) - 2 * srcla(2, 0, -2, n);
            srcmu(2, 0, -1, n) = srcmu(2, 0, -1, n) - 2 * srcla(2, 0, -1, n);
            srcmu(2, 0, 1, n) = srcmu(2, 0, 1, n) - 2 * srcla(2, 0, 1, n);
            srcmu(2, 0, 2, n) = srcmu(2, 0, 2, n) - 2 * srcla(2, 0, 2, n);
          }
          phtermmu = 0;
          for (int k = nk - nb - 1; k <= nk; k++)
            phtermmu += srcmu(2, 0, 0, k) * srcmu(2, 0, 0, k);
          phtermmu += srcmu(2, -2, 0, n) * srcmu(2, -2, 0, n) +
                      srcmu(2, -1, 0, n) * srcmu(2, -1, 0, n) +
                      srcmu(2, 1, 0, n) * srcmu(2, 1, 0, n) +
                      srcmu(2, 2, 0, n) * srcmu(2, 2, 0, n) +
                      srcmu(2, 0, -2, n) * srcmu(2, 0, -2, n) +
                      srcmu(2, 0, -1, n) * srcmu(2, 0, -1, n) +
                      srcmu(2, 0, 1, n) * srcmu(2, 0, 1, n) +
                      srcmu(2, 0, 2, n) * srcmu(2, 0, 2, n);
          phtermla = 0;
          for (int k = nk - nb - 1; k <= nk; k++)
            phtermla += srcla(2, 0, 0, k) * srcla(2, 0, 0, k);
          phtermla += srcla(2, -2, 0, n) * srcla(2, -2, 0, n) +
                      srcla(2, -1, 0, n) * srcla(2, -1, 0, n) +
                      srcla(2, 1, 0, n) * srcla(2, 1, 0, n) +
                      srcla(2, 2, 0, n) * srcla(2, 2, 0, n) +
                      srcla(2, 0, -2, n) * srcla(2, 0, -2, n) +
                      srcla(2, 0, -1, n) * srcla(2, 0, -1, n) +
                      srcla(2, 0, 1, n) * srcla(2, 0, 1, n) +
                      srcla(2, 0, 2, n) * srcla(2, 0, 2, n);
          phtermrho =
              SQR(idt2 * (up(2, l, m, n) - 2 * u(2, l, m, n) + um(2, l, m, n)));
          if (varcase == 1) {
            ph(1, l, m, n) += phtermrho;
            ph(2, l, m, n) += phtermmu;
            ph(3, l, m, n) += phtermla;
          } else if (varcase == 2) {
            float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
            float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
            ph(1, l, m, n) +=
                phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else if (varcase == 3) {
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else {
            float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
            //            ph(3,l,m,n) +=
            //            4*rho(l,m,n)*mu(l,m,n)*(irat*irat*phtermmu+phtermla);
            ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          }

          // dL^w/d(lambda)
          for (int k = 1; k <= nb; k++) {
            ddu = 0;
            for (int q = 1; q <= wb; q++)
              ddu += acof(k, q, nk - n + 1) * u(3, l, m, nk - q + 1);
            srcla(3, 0, 0, nk - k + 1) += ddu;
          }
          if (n == nk) srcla(3, 0, 0, n) += ghcof(1) * u(3, l, m, nk + 1);

          if (n < nk - nb + 3)
            srcla(3, 0, 0, n - 2) += -0.125 * u(3, l, m, n) +
                                     o6 * u(3, l, m, n - 1) -
                                     o24 * u(3, l, m, n - 2);
          if (n < nk - nb + 2)
            srcla(3, 0, 0, n - 1) +=
                o6 * u(3, l, m, n - 2) + 0.5 * u(3, l, m, n) -
                5 * o6 * u(3, l, m, n - 1) + o6 * u(3, l, m, n + 1);
          if (n < nk - nb + 1)
            srcla(3, 0, 0, n) +=
                -0.125 * (u(3, l, m, n - 2) + u(3, l, m, n + 2)) +
                0.5 * (u(3, l, m, n - 1) + u(3, l, m, n + 1)) -
                0.75 * u(3, l, m, n);
          if (n < nk - nb)
            srcla(3, 0, 0, n + 1) +=
                o6 * u(3, l, m, n - 1) - 5 * o6 * u(3, l, m, n + 1) +
                0.5 * u(3, l, m, n) + o6 * u(3, l, m, n + 2);

          for (int q = nk - 9; q <= nk; q++) {
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
            srcmu(3, 0, 0, q) = 2 * srcla(3, 0, 0, q);
          }

          dv = d4b * (u(2, l, m + 2, n) - u(2, l, m - 2, n)) +
               d4a * (u(2, l, m + 1, n) - u(2, l, m - 1, n));
          du = d4b * (u(1, l + 2, m, n) - u(1, l - 2, m, n)) +
               d4a * (u(1, l + 1, m, n) - u(1, l - 1, m, n));
          for (int k = 1; k <= nb; k++)
            srcla(3, 0, 0, nk - k + 1) -= bope(k, nk - n + 1) * (du + dv);
          if (n < nk - nb + 3) srcla(3, 0, 0, n - 2) += d4b * (du + dv);
          if (n < nk - nb + 2) srcla(3, 0, 0, n - 1) += d4a * (du + dv);
          if (n < nk - nb) srcla(3, 0, 0, n + 1) += -d4a * (du + dv);

          // dL^w/d(mu)
          srcmu(3, +2, 0, n) += -0.125 * u(3, l, m, n) +
                                o6 * u(3, l + 1, m, n) -
                                o24 * u(3, l + 2, m, n);
          srcmu(3, +1, 0, n) += o6 * u(3, l - 1, m, n) + 0.5 * u(3, l, m, n) -
                                5 * o6 * u(3, l + 1, m, n) +
                                o6 * u(3, l + 2, m, n);
          srcmu(3, 0, 0, n) +=
              -0.125 * (u(3, l - 2, m, n) + u(3, l + 2, m, n)) +
              0.5 * (u(3, l - 1, m, n) + u(3, l + 1, m, n)) -
              0.75 * u(3, l, m, n);
          srcmu(3, -1, 0, n) += o6 * u(3, l - 2, m, n) -
                                5 * o6 * u(3, l - 1, m, n) +
                                0.5 * u(3, l, m, n) + o6 * u(3, l + 1, m, n);
          srcmu(3, -2, 0, n) += -o24 * u(3, l - 2, m, n) +
                                o6 * u(3, l - 1, m, n) - 0.125 * u(3, l, m, n);

          srcmu(3, 0, +2, n) += -0.125 * u(3, l, m, n) +
                                o6 * u(3, l, m + 1, n) -
                                o24 * u(3, l, m + 2, n);
          srcmu(3, 0, +1, n) += o6 * u(3, l, m - 1, n) + 0.5 * u(3, l, m, n) -
                                5 * o6 * u(3, l, m + 1, n) +
                                o6 * u(3, l, m + 2, n);
          srcmu(3, 0, 0, n) +=
              -0.125 * (u(3, l, m - 2, n) + u(3, l, m + 2, n)) +
              0.5 * (u(3, l, m - 1, n) + u(3, l, m + 1, n)) -
              0.75 * u(3, l, m, n);
          srcmu(3, 0, -1, n) += o6 * u(3, l, m - 2, n) -
                                5 * o6 * u(3, l, m - 1, n) +
                                0.5 * u(3, l, m, n) + o6 * u(3, l, m + 1, n);
          srcmu(3, 0, -2, n) += -o24 * u(3, l, m - 2, n) +
                                o6 * u(3, l, m - 1, n) - 0.125 * u(3, l, m, n);

          if (n < nk - nb + 1) {
            du = d4b * (u(1, l, m, n + 2) - u(1, l, m, n - 2)) +
                 d4a * (u(1, l, m, n + 1) - u(1, l, m, n - 1));
            dv = d4b * (u(2, l, m, n + 2) - u(2, l, m, n - 2)) +
                 d4a * (u(2, l, m, n + 1) - u(2, l, m, n - 1));
          } else {
            du = dv = 0;
            for (int q = 1; q <= wb; q++) {
              du -= bope(nk - n + 1, q) * u(1, l, m, nk - q + 1);
              dv -= bope(nk - n + 1, q) * u(2, l, m, nk - q + 1);
            }
          }

          srcmu(3, +2, 0, n) += -d4b * du;
          srcmu(3, +1, 0, n) += -d4a * du;
          srcmu(3, -1, 0, n) += d4a * du;
          srcmu(3, -2, 0, n) += d4b * du;

          srcmu(3, 0, +2, n) += -d4b * dv;
          srcmu(3, 0, +1, n) += -d4a * dv;
          srcmu(3, 0, -1, n) += d4a * dv;
          srcmu(3, 0, -2, n) += d4b * dv;

          // Square and transform
          if (varcase >= 2) {
            for (int k = nk - nb - 1; k <= nk; k++)
              srcmu(3, 0, 0, k) = srcmu(3, 0, 0, k) - 2 * srcla(3, 0, 0, k);
            srcmu(3, -2, 0, n) = srcmu(3, -2, 0, n) - 2 * srcla(3, -2, 0, n);
            srcmu(3, -1, 0, n) = srcmu(3, -1, 0, n) - 2 * srcla(3, -1, 0, n);
            srcmu(3, 1, 0, n) = srcmu(3, 1, 0, n) - 2 * srcla(3, 1, 0, n);
            srcmu(3, 2, 0, n) = srcmu(3, 2, 0, n) - 2 * srcla(3, 2, 0, n);
            srcmu(3, 0, -2, n) = srcmu(3, 0, -2, n) - 2 * srcla(3, 0, -2, n);
            srcmu(3, 0, -1, n) = srcmu(3, 0, -1, n) - 2 * srcla(3, 0, -1, n);
            srcmu(3, 0, 1, n) = srcmu(3, 0, 1, n) - 2 * srcla(3, 0, 1, n);
            srcmu(3, 0, 2, n) = srcmu(3, 0, 2, n) - 2 * srcla(3, 0, 2, n);
          }
          phtermmu = 0;
          for (int k = nk - nb - 1; k <= nk; k++)
            phtermmu += srcmu(3, 0, 0, k) * srcmu(3, 0, 0, k);
          phtermmu += srcmu(3, -2, 0, n) * srcmu(3, -2, 0, n) +
                      srcmu(3, -1, 0, n) * srcmu(3, -1, 0, n) +
                      srcmu(3, 1, 0, n) * srcmu(3, 1, 0, n) +
                      srcmu(3, 2, 0, n) * srcmu(3, 2, 0, n) +
                      srcmu(3, 0, -2, n) * srcmu(3, 0, -2, n) +
                      srcmu(3, 0, -1, n) * srcmu(3, 0, -1, n) +
                      srcmu(3, 0, 1, n) * srcmu(3, 0, 1, n) +
                      srcmu(3, 0, 2, n) * srcmu(3, 0, 2, n);
          phtermla = 0;
          for (int k = nk - nb - 1; k <= nk; k++)
            phtermla += srcla(3, 0, 0, k) * srcla(3, 0, 0, k);
          phtermla += srcla(3, -2, 0, n) * srcla(3, -2, 0, n) +
                      srcla(3, -1, 0, n) * srcla(3, -1, 0, n) +
                      srcla(3, 1, 0, n) * srcla(3, 1, 0, n) +
                      srcla(3, 2, 0, n) * srcla(3, 2, 0, n) +
                      srcla(3, 0, -2, n) * srcla(3, 0, -2, n) +
                      srcla(3, 0, -1, n) * srcla(3, 0, -1, n) +
                      srcla(3, 0, 1, n) * srcla(3, 0, 1, n) +
                      srcla(3, 0, 2, n) * srcla(3, 0, 2, n);
          phtermrho =
              SQR(idt2 * (up(3, l, m, n) - 2 * u(3, l, m, n) + um(3, l, m, n)));

          if (varcase == 1) {
            ph(1, l, m, n) += phtermrho;
            ph(2, l, m, n) += phtermmu;
            ph(3, l, m, n) += phtermla;
          } else if (varcase == 2) {
            float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
            float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
            ph(1, l, m, n) +=
                phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else if (varcase == 3) {
            ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
            ph(3, l, m, n) += 4 * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          } else {
            float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
            ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                              (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          }
        }
  }

#undef srcla
#undef srcmu

#define srcla(c, i, j, k) \
  srcla_[(i + 2) + 5 * (j + 2) + 25 * (k + 2) + (c - 1) * 125]
#define srcmu(c, i, j, k) \
  srcmu_[(i + 2) + 5 * (j + 2) + 25 * (k + 2) + (c - 1) * 125]

  // Interior operators
  // #pragma omp parallel for
  for (int n = kstart; n <= kend; n++)
    for (int m = jfirstact; m <= jlastact; m++) /* #pragma ivdep */
      for (int l = ifirstact; l <= ilastact; l++) {
        for (int c = 0; c < 3 * 5 * 5 * 5; c++) srcla_[c] = srcmu_[c] = 0;

        // d L^u/d(lambda):
        srcla(1, +2, 0, 0) = -0.125 * u(1, l, m, n) + o6 * u(1, l + 1, m, n) -
                             o24 * u(1, l + 2, m, n);
        srcla(1, +1, 0, 0) = o6 * u(1, l - 1, m, n) + 0.5 * u(1, l, m, n) -
                             5 * o6 * u(1, l + 1, m, n) +
                             o6 * u(1, l + 2, m, n);
        srcla(1, 0, 0, 0) = -0.125 * (u(1, l - 2, m, n) + u(1, l + 2, m, n)) +
                            0.5 * (u(1, l - 1, m, n) + u(1, l + 1, m, n)) -
                            0.75 * u(1, l, m, n);
        srcla(1, -1, 0, 0) = o6 * u(1, l - 2, m, n) -
                             5 * o6 * u(1, l - 1, m, n) + 0.5 * u(1, l, m, n) +
                             o6 * u(1, l + 1, m, n);
        srcla(1, -2, 0, 0) = -o24 * u(1, l - 2, m, n) + o6 * u(1, l - 1, m, n) -
                             0.125 * u(1, l, m, n);

        srcmu(1, +2, 0, 0) = 2 * srcla(1, +2, 0, 0);
        srcmu(1, +1, 0, 0) = 2 * srcla(1, +1, 0, 0);
        srcmu(1, 0, 0, 0) = 2 * srcla(1, 0, 0, 0);
        srcmu(1, -1, 0, 0) = 2 * srcla(1, -1, 0, 0);
        srcmu(1, -2, 0, 0) = 2 * srcla(1, -2, 0, 0);

        float_sw4 dv = d4b * (u(2, l, m + 2, n) - u(2, l, m - 2, n)) +
                       d4a * (u(2, l, m + 1, n) - u(2, l, m - 1, n));
        srcla(1, +2, 0, 0) += -d4b * dv;
        srcla(1, +1, 0, 0) += -d4a * dv;
        srcla(1, -1, 0, 0) += d4a * dv;
        srcla(1, -2, 0, 0) += d4b * dv;

        float_sw4 dw = d4b * (u(3, l, m, n + 2) - u(3, l, m, n - 2)) +
                       d4a * (u(3, l, m, n + 1) - u(3, l, m, n - 1));
        srcla(1, +2, 0, 0) += -d4b * dw;
        srcla(1, +1, 0, 0) += -d4a * dw;
        srcla(1, -1, 0, 0) += d4a * dw;
        srcla(1, -2, 0, 0) += d4b * dw;

        // dL^u/d(mu)
        srcmu(1, 0, +2, 0) += -0.125 * u(1, l, m, n) + o6 * u(1, l, m + 1, n) -
                              o24 * u(1, l, m + 2, n);
        srcmu(1, 0, +1, 0) += o6 * u(1, l, m - 1, n) + 0.5 * u(1, l, m, n) -
                              5 * o6 * u(1, l, m + 1, n) +
                              o6 * u(1, l, m + 2, n);
        srcmu(1, 0, 0, 0) += -0.125 * (u(1, l, m - 2, n) + u(1, l, m + 2, n)) +
                             0.5 * (u(1, l, m - 1, n) + u(1, l, m + 1, n)) -
                             0.75 * u(1, l, m, n);
        srcmu(1, 0, -1, 0) += o6 * u(1, l, m - 2, n) -
                              5 * o6 * u(1, l, m - 1, n) + 0.5 * u(1, l, m, n) +
                              o6 * u(1, l, m + 1, n);
        srcmu(1, 0, -2, 0) += -o24 * u(1, l, m - 2, n) +
                              o6 * u(1, l, m - 1, n) - 0.125 * u(1, l, m, n);

        srcmu(1, 0, 0, +2) += -0.125 * u(1, l, m, n) + o6 * u(1, l, m, n + 1) -
                              o24 * u(1, l, m, n + 2);
        srcmu(1, 0, 0, +1) += o6 * u(1, l, m, n - 1) + 0.5 * u(1, l, m, n) -
                              5 * o6 * u(1, l, m, n + 1) +
                              o6 * u(1, l, m, n + 2);
        srcmu(1, 0, 0, 0) += -0.125 * (u(1, l, m, n - 2) + u(1, l, m, n + 2)) +
                             0.5 * (u(1, l, m, n - 1) + u(1, l, m, n + 1)) -
                             0.75 * u(1, l, m, n);
        srcmu(1, 0, 0, -1) += o6 * u(1, l, m, n - 2) -
                              5 * o6 * u(1, l, m, n - 1) + 0.5 * u(1, l, m, n) +
                              o6 * u(1, l, m, n + 1);
        srcmu(1, 0, 0, -2) += -o24 * u(1, l, m, n - 2) +
                              o6 * u(1, l, m, n - 1) - 0.125 * u(1, l, m, n);

        dv = d4b * (u(2, l + 2, m, n) - u(2, l - 2, m, n)) +
             d4a * (u(2, l + 1, m, n) - u(2, l - 1, m, n));
        srcmu(1, 0, +2, 0) += -d4b * dv;
        srcmu(1, 0, +1, 0) += -d4a * dv;
        srcmu(1, 0, -1, 0) += d4a * dv;
        srcmu(1, 0, -2, 0) += d4b * dv;

        dw = d4b * (u(3, l + 2, m, n) - u(3, l - 2, m, n)) +
             d4a * (u(3, l + 1, m, n) - u(3, l - 1, m, n));
        srcmu(1, 0, 0, +2) += -d4b * dw;
        srcmu(1, 0, 0, +1) += -d4a * dw;
        srcmu(1, 0, 0, -1) += d4a * dw;
        srcmu(1, 0, 0, -2) += d4b * dw;

        // Transform and square
        if (varcase >= 2) {
          srcmu(1, 0, 0, -2) = srcmu(1, 0, 0, -2) - 2 * srcla(1, 0, 0, -2);
          srcmu(1, 0, 0, -1) = srcmu(1, 0, 0, -1) - 2 * srcla(1, 0, 0, -1);
          srcmu(1, 0, 0, 0) = srcmu(1, 0, 0, 0) - 2 * srcla(1, 0, 0, 0);
          srcmu(1, 0, 0, 1) = srcmu(1, 0, 0, 1) - 2 * srcla(1, 0, 0, 1);
          srcmu(1, 0, 0, 2) = srcmu(1, 0, 0, 2) - 2 * srcla(1, 0, 0, 2);
          srcmu(1, -2, 0, 0) = srcmu(1, -2, 0, 0) - 2 * srcla(1, -2, 0, 0);
          srcmu(1, -1, 0, 0) = srcmu(1, -1, 0, 0) - 2 * srcla(1, -1, 0, 0);
          srcmu(1, 1, 0, 0) = srcmu(1, 1, 0, 0) - 2 * srcla(1, 1, 0, 0);
          srcmu(1, 2, 0, 0) = srcmu(1, 2, 0, 0) - 2 * srcla(1, 2, 0, 0);
          srcmu(1, 0, -2, 0) = srcmu(1, 0, -2, 0) - 2 * srcla(1, 0, -2, 0);
          srcmu(1, 0, -1, 0) = srcmu(1, 0, -1, 0) - 2 * srcla(1, 0, -1, 0);
          srcmu(1, 0, 1, 0) = srcmu(1, 0, 1, 0) - 2 * srcla(1, 0, 1, 0);
          srcmu(1, 0, 2, 0) = srcmu(1, 0, 2, 0) - 2 * srcla(1, 0, 2, 0);
        }
        float_sw4 phtermmu = 0;
        phtermmu += srcmu(1, 0, 0, -2) * srcmu(1, 0, 0, -2) +
                    srcmu(1, 0, 0, -1) * srcmu(1, 0, 0, -1) +
                    srcmu(1, 0, 0, 0) * srcmu(1, 0, 0, 0) +
                    srcmu(1, 0, 0, 1) * srcmu(1, 0, 0, 1) +
                    srcmu(1, 0, 0, 2) * srcmu(1, 0, 0, 2);
        phtermmu += srcmu(1, -2, 0, 0) * srcmu(1, -2, 0, 0) +
                    srcmu(1, -1, 0, 0) * srcmu(1, -1, 0, 0) +
                    srcmu(1, 1, 0, 0) * srcmu(1, 1, 0, 0) +
                    srcmu(1, 2, 0, 0) * srcmu(1, 2, 0, 0) +
                    srcmu(1, 0, -2, 0) * srcmu(1, 0, -2, 0) +
                    srcmu(1, 0, -1, 0) * srcmu(1, 0, -1, 0) +
                    srcmu(1, 0, 1, 0) * srcmu(1, 0, 1, 0) +
                    srcmu(1, 0, 2, 0) * srcmu(1, 0, 2, 0);
        float_sw4 phtermla = 0;
        phtermla += srcla(1, 0, 0, -2) * srcla(1, 0, 0, -2) +
                    srcla(1, 0, 0, -1) * srcla(1, 0, 0, -1) +
                    srcla(1, 0, 0, 0) * srcla(1, 0, 0, 0) +
                    srcla(1, 0, 0, 1) * srcla(1, 0, 0, 1) +
                    srcla(1, 0, 0, 2) * srcla(1, 0, 0, 2);
        phtermla += srcla(1, -2, 0, 0) * srcla(1, -2, 0, 0) +
                    srcla(1, -1, 0, 0) * srcla(1, -1, 0, 0) +
                    srcla(1, 1, 0, 0) * srcla(1, 1, 0, 0) +
                    srcla(1, 2, 0, 0) * srcla(1, 2, 0, 0) +
                    srcla(1, 0, -2, 0) * srcla(1, 0, -2, 0) +
                    srcla(1, 0, -1, 0) * srcla(1, 0, -1, 0) +
                    srcla(1, 0, 1, 0) * srcla(1, 0, 1, 0) +
                    srcla(1, 0, 2, 0) * srcla(1, 0, 2, 0);
        float_sw4 phtermrho =
            SQR(idt2 * (up(1, l, m, n) - 2 * u(1, l, m, n) + um(1, l, m, n)));
        if (varcase == 1) {
          ph(1, l, m, n) += phtermrho;
          ph(2, l, m, n) += phtermmu;
          ph(3, l, m, n) += phtermla;
        } else if (varcase == 2) {
          float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
          float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
          ph(1, l, m, n) +=
              phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
          ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
          ph(3, l, m, n) +=
              4 * rho(l, m, n) * (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
        } else if (varcase == 3) {
          ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
          ph(3, l, m, n) +=
              4 * rho(l, m, n) * (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
        } else {
          float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
          //            ph(3,l,m,n) +=
          //            4*rho(l,m,n)*mu(l,m,n)*(irat*irat*phtermmu+phtermla);
          ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                            (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
        }

        // d L^v/d(lambda):
        srcla(2, 0, +2, 0) = -0.125 * u(2, l, m, n) + o6 * u(2, l, m + 1, n) -
                             o24 * u(2, l, m + 2, n);
        srcla(2, 0, +1, 0) = o6 * u(2, l, m - 1, n) + 0.5 * u(2, l, m, n) -
                             5 * o6 * u(2, l, m + 1, n) +
                             o6 * u(2, l, m + 2, n);
        srcla(2, 0, 0, 0) = -0.125 * (u(2, l, m - 2, n) + u(2, l, m + 2, n)) +
                            0.5 * (u(2, l, m - 1, n) + u(2, l, m + 1, n)) -
                            0.75 * u(2, l, m, n);
        srcla(2, 0, -1, 0) = o6 * u(2, l, m - 2, n) -
                             5 * o6 * u(2, l, m - 1, n) + 0.5 * u(2, l, m, n) +
                             o6 * u(2, l, m + 1, n);
        srcla(2, 0, -2, 0) = -o24 * u(2, l, m - 2, n) + o6 * u(2, l, m - 1, n) -
                             0.125 * u(2, l, m, n);

        srcmu(2, 0, +2, 0) = 2 * srcla(2, 0, +2, 0);
        srcmu(2, 0, +1, 0) = 2 * srcla(2, 0, +1, 0);
        srcmu(2, 0, 0, 0) = 2 * srcla(2, 0, 0, 0);
        srcmu(2, 0, -1, 0) = 2 * srcla(2, 0, -1, 0);
        srcmu(2, 0, -2, 0) = 2 * srcla(2, 0, -2, 0);

        float_sw4 du = d4b * (u(1, l + 2, m, n) - u(1, l - 2, m, n)) +
                       d4a * (u(1, l + 1, m, n) - u(1, l - 1, m, n));
        srcla(2, 0, +2, 0) += -d4b * du;
        srcla(2, 0, +1, 0) += -d4a * du;
        srcla(2, 0, -1, 0) += d4a * du;
        srcla(2, 0, -2, 0) += d4b * du;

        dw = d4b * (u(3, l, m, n + 2) - u(3, l, m, n - 2)) +
             d4a * (u(3, l, m, n + 1) - u(3, l, m, n - 1));
        srcla(2, 0, +2, 0) += -d4b * dw;
        srcla(2, 0, +1, 0) += -d4a * dw;
        srcla(2, 0, -1, 0) += d4a * dw;
        srcla(2, 0, -2, 0) += d4b * dw;

        // dL^v/d(mu)
        srcmu(2, +2, 0, 0) += -0.125 * u(2, l, m, n) + o6 * u(2, l + 1, m, n) -
                              o24 * u(2, l + 2, m, n);
        srcmu(2, +1, 0, 0) += o6 * u(2, l - 1, m, n) + 0.5 * u(2, l, m, n) -
                              5 * o6 * u(2, l + 1, m, n) +
                              o6 * u(2, l + 2, m, n);
        srcmu(2, 0, 0, 0) += -0.125 * (u(2, l - 2, m, n) + u(2, l + 2, m, n)) +
                             0.5 * (u(2, l - 1, m, n) + u(2, l + 1, m, n)) -
                             0.75 * u(2, l, m, n);
        srcmu(2, -1, 0, 0) += o6 * u(2, l - 2, m, n) -
                              5 * o6 * u(2, l - 1, m, n) + 0.5 * u(2, l, m, n) +
                              o6 * u(2, l + 1, m, n);
        srcmu(2, -2, 0, 0) += -o24 * u(2, l - 2, m, n) +
                              o6 * u(2, l - 1, m, n) - 0.125 * u(2, l, m, n);

        srcmu(2, 0, 0, +2) += -0.125 * u(2, l, m, n) + o6 * u(2, l, m, n + 1) -
                              o24 * u(2, l, m, n + 2);
        srcmu(2, 0, 0, +1) += o6 * u(2, l, m, n - 1) + 0.5 * u(2, l, m, n) -
                              5 * o6 * u(2, l, m, n + 1) +
                              o6 * u(2, l, m, n + 2);
        srcmu(2, 0, 0, 0) += -0.125 * (u(2, l, m, n - 2) + u(2, l, m, n + 2)) +
                             0.5 * (u(2, l, m, n - 1) + u(2, l, m, n + 1)) -
                             0.75 * u(2, l, m, n);
        srcmu(2, 0, 0, -1) += o6 * u(2, l, m, n - 2) -
                              5 * o6 * u(2, l, m, n - 1) + 0.5 * u(2, l, m, n) +
                              o6 * u(2, l, m, n + 1);
        srcmu(2, 0, 0, -2) += -o24 * u(2, l, m, n - 2) +
                              o6 * u(2, l, m, n - 1) - 0.125 * u(2, l, m, n);

        du = d4b * (u(1, l, m + 2, n) - u(1, l, m - 2, n)) +
             d4a * (u(1, l, m + 1, n) - u(1, l, m - 1, n));
        srcmu(2, +2, 0, 0) += -d4b * du;
        srcmu(2, +1, 0, 0) += -d4a * du;
        srcmu(2, -1, 0, 0) += d4a * du;
        srcmu(2, -2, 0, 0) += d4b * du;

        dw = d4b * (u(3, l, m + 2, n) - u(3, l, m - 2, n)) +
             d4a * (u(3, l, m + 1, n) - u(3, l, m - 1, n));
        srcmu(2, 0, 0, +2) += -d4b * dw;
        srcmu(2, 0, 0, +1) += -d4a * dw;
        srcmu(2, 0, 0, -1) += d4a * dw;
        srcmu(2, 0, 0, -2) += d4b * dw;

        // Transform and square
        if (varcase >= 2) {
          srcmu(2, 0, 0, -2) = srcmu(2, 0, 0, -2) - 2 * srcla(2, 0, 0, -2);
          srcmu(2, 0, 0, -1) = srcmu(2, 0, 0, -1) - 2 * srcla(2, 0, 0, -1);
          srcmu(2, 0, 0, 0) = srcmu(2, 0, 0, 0) - 2 * srcla(2, 0, 0, 0);
          srcmu(2, 0, 0, 1) = srcmu(2, 0, 0, 1) - 2 * srcla(2, 0, 0, 1);
          srcmu(2, 0, 0, 2) = srcmu(2, 0, 0, 2) - 2 * srcla(2, 0, 0, 2);
          srcmu(2, -2, 0, 0) = srcmu(2, -2, 0, 0) - 2 * srcla(2, -2, 0, 0);
          srcmu(2, -1, 0, 0) = srcmu(2, -1, 0, 0) - 2 * srcla(2, -1, 0, 0);
          srcmu(2, 1, 0, 0) = srcmu(2, 1, 0, 0) - 2 * srcla(2, 1, 0, 0);
          srcmu(2, 2, 0, 0) = srcmu(2, 2, 0, 0) - 2 * srcla(2, 2, 0, 0);
          srcmu(2, 0, -2, 0) = srcmu(2, 0, -2, 0) - 2 * srcla(2, 0, -2, 0);
          srcmu(2, 0, -1, 0) = srcmu(2, 0, -1, 0) - 2 * srcla(2, 0, -1, 0);
          srcmu(2, 0, 1, 0) = srcmu(2, 0, 1, 0) - 2 * srcla(2, 0, 1, 0);
          srcmu(2, 0, 2, 0) = srcmu(2, 0, 2, 0) - 2 * srcla(2, 0, 2, 0);
        }
        phtermmu = 0;
        phtermmu += srcmu(2, 0, 0, -2) * srcmu(2, 0, 0, -2) +
                    srcmu(2, 0, 0, -1) * srcmu(2, 0, 0, -1) +
                    srcmu(2, 0, 0, 0) * srcmu(2, 0, 0, 0) +
                    srcmu(2, 0, 0, 1) * srcmu(2, 0, 0, 1) +
                    srcmu(2, 0, 0, 2) * srcmu(2, 0, 0, 2);
        phtermmu += srcmu(2, -2, 0, 0) * srcmu(2, -2, 0, 0) +
                    srcmu(2, -1, 0, 0) * srcmu(2, -1, 0, 0) +
                    srcmu(2, 1, 0, 0) * srcmu(2, 1, 0, 0) +
                    srcmu(2, 2, 0, 0) * srcmu(2, 2, 0, 0) +
                    srcmu(2, 0, -2, 0) * srcmu(2, 0, -2, 0) +
                    srcmu(2, 0, -1, 0) * srcmu(2, 0, -1, 0) +
                    srcmu(2, 0, 1, 0) * srcmu(2, 0, 1, 0) +
                    srcmu(2, 0, 2, 0) * srcmu(2, 0, 2, 0);
        phtermla = 0;
        phtermla += srcla(2, 0, 0, -2) * srcla(2, 0, 0, -2) +
                    srcla(2, 0, 0, -1) * srcla(2, 0, 0, -1) +
                    srcla(2, 0, 0, 0) * srcla(2, 0, 0, 0) +
                    srcla(2, 0, 0, 1) * srcla(2, 0, 0, 1) +
                    srcla(2, 0, 0, 2) * srcla(2, 0, 0, 2);
        phtermla += srcla(2, -2, 0, 0) * srcla(2, -2, 0, 0) +
                    srcla(2, -1, 0, 0) * srcla(2, -1, 0, 0) +
                    srcla(2, 1, 0, 0) * srcla(2, 1, 0, 0) +
                    srcla(2, 2, 0, 0) * srcla(2, 2, 0, 0) +
                    srcla(2, 0, -2, 0) * srcla(2, 0, -2, 0) +
                    srcla(2, 0, -1, 0) * srcla(2, 0, -1, 0) +
                    srcla(2, 0, 1, 0) * srcla(2, 0, 1, 0) +
                    srcla(2, 0, 2, 0) * srcla(2, 0, 2, 0);
        phtermrho =
            SQR(idt2 * (up(2, l, m, n) - 2 * u(2, l, m, n) + um(2, l, m, n)));
        if (varcase == 1) {
          ph(1, l, m, n) += phtermrho;
          ph(2, l, m, n) += phtermmu;
          ph(3, l, m, n) += phtermla;
        } else if (varcase == 2) {
          float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
          float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
          ph(1, l, m, n) +=
              phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
          ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
          ph(3, l, m, n) +=
              4 * rho(l, m, n) * (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
        } else if (varcase == 3) {
          ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
          ph(3, l, m, n) +=
              4 * rho(l, m, n) * (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
        } else {
          float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
          ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                            (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          //            ph(3,l,m,n) +=
          //            4*rho(l,m,n)*mu(l,m,n)*(irat*irat*phtermmu+phtermla);
        }

        // dL^w/d(lambda)
        srcla(3, 0, 0, +2) = -0.125 * u(3, l, m, n) + o6 * u(3, l, m, n + 1) -
                             o24 * u(3, l, m, n + 2);
        srcla(3, 0, 0, +1) = o6 * u(3, l, m, n - 1) + 0.5 * u(3, l, m, n) -
                             5 * o6 * u(3, l, m, n + 1) +
                             o6 * u(3, l, m, n + 2);
        srcla(3, 0, 0, 0) = -0.125 * (u(3, l, m, n - 2) + u(3, l, m, n + 2)) +
                            0.5 * (u(3, l, m, n - 1) + u(3, l, m, n + 1)) -
                            0.75 * u(3, l, m, n);
        srcla(3, 0, 0, -1) = o6 * u(3, l, m, n - 2) -
                             5 * o6 * u(3, l, m, n - 1) + 0.5 * u(3, l, m, n) +
                             o6 * u(3, l, m, n + 1);
        srcla(3, 0, 0, -2) = -o24 * u(3, l, m, n - 2) + o6 * u(3, l, m, n - 1) -
                             0.125 * u(3, l, m, n);

        srcmu(3, 0, 0, +2) = 2 * srcla(3, 0, 0, +2);
        srcmu(3, 0, 0, +1) = 2 * srcla(3, 0, 0, +1);
        srcmu(3, 0, 0, 0) = 2 * srcla(3, 0, 0, 0);
        srcmu(3, 0, 0, -1) = 2 * srcla(3, 0, 0, -1);
        srcmu(3, 0, 0, -2) = 2 * srcla(3, 0, 0, -2);

        du = d4b * (u(1, l + 2, m, n) - u(1, l - 2, m, n)) +
             d4a * (u(1, l + 1, m, n) - u(1, l - 1, m, n));
        srcla(3, 0, 0, +2) += -d4b * du;
        srcla(3, 0, 0, +1) += -d4a * du;
        srcla(3, 0, 0, -1) += d4a * du;
        srcla(3, 0, 0, -2) += d4b * du;

        dv = d4b * (u(2, l, m + 2, n) - u(2, l, m - 2, n)) +
             d4a * (u(2, l, m + 1, n) - u(2, l, m - 1, n));
        srcla(3, 0, 0, +2) += -d4b * dv;
        srcla(3, 0, 0, +1) += -d4a * dv;
        srcla(3, 0, 0, -1) += d4a * dv;
        srcla(3, 0, 0, -2) += d4b * dv;

        // dL^w/d(mu)
        srcmu(3, +2, 0, 0) += -0.125 * u(3, l, m, n) + o6 * u(3, l + 1, m, n) -
                              o24 * u(3, l + 2, m, n);
        srcmu(3, +1, 0, 0) += o6 * u(3, l - 1, m, n) + 0.5 * u(3, l, m, n) -
                              5 * o6 * u(3, l + 1, m, n) +
                              o6 * u(3, l + 2, m, n);
        srcmu(3, 0, 0, 0) += -0.125 * (u(3, l - 2, m, n) + u(3, l + 2, m, n)) +
                             0.5 * (u(3, l - 1, m, n) + u(3, l + 1, m, n)) -
                             0.75 * u(3, l, m, n);
        srcmu(3, -1, 0, 0) += o6 * u(3, l - 2, m, n) -
                              5 * o6 * u(3, l - 1, m, n) + 0.5 * u(3, l, m, n) +
                              o6 * u(3, l + 1, m, n);
        srcmu(3, -2, 0, 0) += -o24 * u(3, l - 2, m, n) +
                              o6 * u(3, l - 1, m, n) - 0.125 * u(3, l, m, n);

        srcmu(3, 0, +2, 0) += -0.125 * u(3, l, m, n) + o6 * u(3, l, m + 1, n) -
                              o24 * u(3, l, m + 2, n);
        srcmu(3, 0, +1, 0) += o6 * u(3, l, m - 1, n) + 0.5 * u(3, l, m, n) -
                              5 * o6 * u(3, l, m + 1, n) +
                              o6 * u(3, l, m + 2, n);
        srcmu(3, 0, 0, 0) += -0.125 * (u(3, l, m - 2, n) + u(3, l, m + 2, n)) +
                             0.5 * (u(3, l, m - 1, n) + u(3, l, m + 1, n)) -
                             0.75 * u(3, l, m, n);
        srcmu(3, 0, -1, 0) += o6 * u(3, l, m - 2, n) -
                              5 * o6 * u(3, l, m - 1, n) + 0.5 * u(3, l, m, n) +
                              o6 * u(3, l, m + 1, n);
        srcmu(3, 0, -2, 0) += -o24 * u(3, l, m - 2, n) +
                              o6 * u(3, l, m - 1, n) - 0.125 * u(3, l, m, n);

        du = d4b * (u(1, l, m, n + 2) - u(1, l, m, n - 2)) +
             d4a * (u(1, l, m, n + 1) - u(1, l, m, n - 1));
        srcmu(3, +2, 0, 0) += -d4b * du;
        srcmu(3, +1, 0, 0) += -d4a * du;
        srcmu(3, -1, 0, 0) += d4a * du;
        srcmu(3, -2, 0, 0) += d4b * du;

        dv = d4b * (u(2, l, m, n + 2) - u(2, l, m, n - 2)) +
             d4a * (u(2, l, m, n + 1) - u(2, l, m, n - 1));
        srcmu(3, 0, +2, 0) += -d4b * dv;
        srcmu(3, 0, +1, 0) += -d4a * dv;
        srcmu(3, 0, -1, 0) += d4a * dv;
        srcmu(3, 0, -2, 0) += d4b * dv;

        // Transform and square
        if (varcase >= 2) {
          srcmu(3, 0, 0, -2) = srcmu(3, 0, 0, -2) - 2 * srcla(3, 0, 0, -2);
          srcmu(3, 0, 0, -1) = srcmu(3, 0, 0, -1) - 2 * srcla(3, 0, 0, -1);
          srcmu(3, 0, 0, 0) = srcmu(3, 0, 0, 0) - 2 * srcla(3, 0, 0, 0);
          srcmu(3, 0, 0, 1) = srcmu(3, 0, 0, 1) - 2 * srcla(3, 0, 0, 1);
          srcmu(3, 0, 0, 2) = srcmu(3, 0, 0, 2) - 2 * srcla(3, 0, 0, 2);
          srcmu(3, -2, 0, 0) = srcmu(3, -2, 0, 0) - 2 * srcla(3, -2, 0, 0);
          srcmu(3, -1, 0, 0) = srcmu(3, -1, 0, 0) - 2 * srcla(3, -1, 0, 0);
          srcmu(3, 1, 0, 0) = srcmu(3, 1, 0, 0) - 2 * srcla(3, 1, 0, 0);
          srcmu(3, 2, 0, 0) = srcmu(3, 2, 0, 0) - 2 * srcla(3, 2, 0, 0);
          srcmu(3, 0, -2, 0) = srcmu(3, 0, -2, 0) - 2 * srcla(3, 0, -2, 0);
          srcmu(3, 0, -1, 0) = srcmu(3, 0, -1, 0) - 2 * srcla(3, 0, -1, 0);
          srcmu(3, 0, 1, 0) = srcmu(3, 0, 1, 0) - 2 * srcla(3, 0, 1, 0);
          srcmu(3, 0, 2, 0) = srcmu(3, 0, 2, 0) - 2 * srcla(3, 0, 2, 0);
        }
        phtermmu = 0;
        phtermmu += srcmu(3, 0, 0, -2) * srcmu(3, 0, 0, -2) +
                    srcmu(3, 0, 0, -1) * srcmu(3, 0, 0, -1) +
                    srcmu(3, 0, 0, 0) * srcmu(3, 0, 0, 0) +
                    srcmu(3, 0, 0, 1) * srcmu(3, 0, 0, 1) +
                    srcmu(3, 0, 0, 2) * srcmu(3, 0, 0, 2);
        phtermmu += srcmu(3, -2, 0, 0) * srcmu(3, -2, 0, 0) +
                    srcmu(3, -1, 0, 0) * srcmu(3, -1, 0, 0) +
                    srcmu(3, 1, 0, 0) * srcmu(3, 1, 0, 0) +
                    srcmu(3, 2, 0, 0) * srcmu(3, 2, 0, 0) +
                    srcmu(3, 0, -2, 0) * srcmu(3, 0, -2, 0) +
                    srcmu(3, 0, -1, 0) * srcmu(3, 0, -1, 0) +
                    srcmu(3, 0, 1, 0) * srcmu(3, 0, 1, 0) +
                    srcmu(3, 0, 2, 0) * srcmu(3, 0, 2, 0);
        phtermla = 0;
        phtermla += srcla(3, 0, 0, -2) * srcla(3, 0, 0, -2) +
                    srcla(3, 0, 0, -1) * srcla(3, 0, 0, -1) +
                    srcla(3, 0, 0, 0) * srcla(3, 0, 0, 0) +
                    srcla(3, 0, 0, 1) * srcla(3, 0, 0, 1) +
                    srcla(3, 0, 0, 2) * srcla(3, 0, 0, 2);
        phtermla += srcla(3, -2, 0, 0) * srcla(3, -2, 0, 0) +
                    srcla(3, -1, 0, 0) * srcla(3, -1, 0, 0) +
                    srcla(3, 1, 0, 0) * srcla(3, 1, 0, 0) +
                    srcla(3, 2, 0, 0) * srcla(3, 2, 0, 0) +
                    srcla(3, 0, -2, 0) * srcla(3, 0, -2, 0) +
                    srcla(3, 0, -1, 0) * srcla(3, 0, -1, 0) +
                    srcla(3, 0, 1, 0) * srcla(3, 0, 1, 0) +
                    srcla(3, 0, 2, 0) * srcla(3, 0, 2, 0);
        phtermrho =
            SQR(idt2 * (up(3, l, m, n) - 2 * u(3, l, m, n) + um(3, l, m, n)));
        if (varcase == 1) {
          ph(1, l, m, n) += phtermrho;
          ph(2, l, m, n) += phtermmu;
          ph(3, l, m, n) += phtermla;
        } else if (varcase == 2) {
          float_sw4 cp2 = (2 * mu(l, m, n) + lambda(l, m, n)) / rho(l, m, n);
          float_sw4 cs2 = mu(l, m, n) / rho(l, m, n);
          ph(1, l, m, n) +=
              phtermrho + cs2 * cs2 * phtermmu + cp2 * cp2 * phtermla;
          ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
          ph(3, l, m, n) +=
              4 * rho(l, m, n) * (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
        } else if (varcase == 3) {
          ph(2, l, m, n) += 4 * rho(l, m, n) * mu(l, m, n) * phtermmu;
          ph(3, l, m, n) +=
              4 * rho(l, m, n) * (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
        } else {
          float_sw4 irat = mu(l, m, n) / (2 * mu(l, m, n) + lambda(l, m, n));
          ph(3, l, m, n) += 16 * SQR(1 - 2 * irat * irat) * rho(l, m, n) *
                            (2 * mu(l, m, n) + lambda(l, m, n)) * phtermla;
          //            ph(3,l,m,n) +=
          //            4*rho(l,m,n)*mu(l,m,n)*(irat*irat*phtermmu+phtermla);
        }
      }
#undef srcla
#undef srcmu
  delete[] srcmu_;
  delete[] srcla_;
}
