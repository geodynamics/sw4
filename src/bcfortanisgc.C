#include "EW.h"
#include "F77_FUNC.h"

extern "C" {
void F77_FUNC(dgesv, DGESV)(sw4_type*, sw4_type*, double*, sw4_type*, sw4_type*, double*, sw4_type*,
                            sw4_type*);
}

void EW::bcfortanisg_ci(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                        int wind[36], sw4_type nx, sw4_type ny, sw4_type nz, float_sw4* u,
                        float_sw4 h, boundaryConditionType bccnd[6],
                        float_sw4 sbop[5], float_sw4* c, float_sw4* bforce1,
                        float_sw4* bforce2, float_sw4* bforce3,
                        float_sw4* bforce4, float_sw4* bforce5,
                        float_sw4* bforce6, float_sw4* strx, float_sw4* stry) {
  const float_sw4 d4a = 2.0 / 3.0;
  const float_sw4 d4b = -1.0 / 12.0;
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  const size_t npts =
      static_cast<size_t>((ie - ib + 1)) * (je - jb + 1) * (ke - kb + 1);
  for (sw4_type s = 0; s < 6; s++) {
    if (bccnd[s] == bDirichlet || bccnd[s] == bSuperGrid) {
      size_t idel = 1 + wind[1 + 6 * s] - wind[6 * s];
      size_t ijdel = idel * (1 + wind[3 + 6 * s] - wind[2 + 6 * s]);
      if (s == 0) {
        //            #pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce1[3 * qq];
              u[ind + npts] = bforce1[1 + 3 * qq];
              u[ind + 2 * npts] = bforce1[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 1) {
        //            #pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce2[3 * qq];
              u[ind + npts] = bforce2[1 + 3 * qq];
              u[ind + 2 * npts] = bforce2[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 2) {
        //            #pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce3[3 * qq];
              u[ind + npts] = bforce3[1 + 3 * qq];
              u[ind + 2 * npts] = bforce3[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 3) {
        //            #pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce4[3 * qq];
              u[ind + npts] = bforce4[1 + 3 * qq];
              u[ind + 2 * npts] = bforce4[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 4) {
        //            #pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce5[3 * qq];
              u[ind + npts] = bforce5[1 + 3 * qq];
              u[ind + 2 * npts] = bforce5[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 5) {
        //            #pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce6[3 * qq];
              u[ind + npts] = bforce6[1 + 3 * qq];
              u[ind + 2 * npts] = bforce6[2 + 3 * qq];
              qq++;
            }
          }
        }
      }
    } else if (bccnd[s] == bPeriodic) {
      if (s == 0) {
#pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + nx;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 1) {
#pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - nx;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 2) {
#pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + ni * ny;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 3) {
#pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - ni * ny;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 4) {
#pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + nij * nz;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 5) {
#pragma omp parallel for
        for (sw4_type k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (sw4_type j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (sw4_type i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - nij * nz;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      }
    } else if (bccnd[s] == bStressFree) {
      REQUIRE2(s == 4 || s == 5,
               "EW::bcfortanisg_ci,  ERROR: Free surface condition"
                   << " not implemented for side " << s << endl);
      double s0i = 1 / sbop[0];
      if (s == 4) {
        sw4_type k = 1;
#pragma omp parallel for
        for (sw4_type j = jb + 2; j <= je - 2; j++)
          for (sw4_type i = ib + 2; i <= ie - 2; i++) {
            size_t qq = i - ib + ni * (j - jb);
            size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
            float_sw4 du = strx[i - ib] * (d4a * (u[ind + 1] - u[ind - 1]) +
                                           d4b * (u[ind + 2] - u[ind - 2]));
            float_sw4 dv =
                strx[i - ib] * (d4a * (u[npts + ind + 1] - u[npts + ind - 1]) +
                                d4b * (u[npts + ind + 2] - u[npts + ind - 2]));
            float_sw4 dw =
                strx[i - ib] *
                (d4a * (u[2 * npts + ind + 1] - u[2 * npts + ind - 1]) +
                 d4b * (u[2 * npts + ind + 2] - u[2 * npts + ind - 2]));

            float_sw4 rhs1 = c[ind + 2 * npts] * du + c[ind + 7 * npts] * dv +
                             c[ind + 11 * npts] * dw;
            float_sw4 rhs2 = c[ind + 4 * npts] * du + c[ind + 9 * npts] * dv +
                             c[ind + 13 * npts] * dw;
            float_sw4 rhs3 = c[ind + 5 * npts] * du + c[ind + 10 * npts] * dv +
                             c[ind + 14 * npts] * dw;

            du = stry[j - jb] * (d4a * (u[ind + ni] - u[ind - ni]) +
                                 d4b * (u[ind + 2 * ni] - u[ind - 2 * ni]));
            dv = stry[j - jb] *
                 (d4a * (u[npts + ind + ni] - u[npts + ind - ni]) +
                  d4b * (u[npts + ind + 2 * ni] - u[npts + ind - 2 * ni]));
            dw = stry[j - jb] *
                 (d4a * (u[2 * npts + ind + ni] - u[2 * npts + ind - ni]) +
                  d4b * (u[2 * npts + ind + 2 * ni] -
                         u[2 * npts + ind - 2 * ni]));
            rhs1 += c[ind + 7 * npts] * du + c[ind + 12 * npts] * dv +
                    c[ind + 13 * npts] * dw;
            rhs2 += c[ind + 9 * npts] * du + c[ind + 16 * npts] * dv +
                    c[ind + 18 * npts] * dw;
            rhs3 += c[ind + 10 * npts] * du + c[ind + 17 * npts] * dv +
                    c[ind + 19 * npts] * dw;

            du = dv = dw = 0;
            for (sw4_type w = 1; w <= 4; w++) {
              du += sbop[w] * u[ind + nij * (w - 1)];
              dv += sbop[w] * u[npts + ind + nij * (w - 1)];
              dw += sbop[w] * u[2 * npts + ind + nij * (w - 1)];
            }
            rhs1 += c[ind + 11 * npts] * du + c[ind + 13 * npts] * dv +
                    c[ind + 14 * npts] * dw - h * bforce5[3 * qq];
            rhs2 += c[ind + 13 * npts] * du + c[ind + 18 * npts] * dv +
                    c[ind + 19 * npts] * dw - h * bforce5[1 + 3 * qq];
            rhs3 += c[ind + 14 * npts] * du + c[ind + 19 * npts] * dv +
                    c[ind + 20 * npts] * dw - h * bforce5[2 + 3 * qq];
            // Solve system for ghost point values
            float_sw4 x[3] = {rhs1, rhs2, rhs3};
            float_sw4 a[9];
            a[0] = c[ind + 11 * npts];
            a[1] = c[ind + 13 * npts];
            a[2] = c[ind + 14 * npts];
            a[3] = c[ind + 13 * npts];
            a[4] = c[ind + 18 * npts];
            a[5] = c[ind + 19 * npts];
            a[6] = c[ind + 14 * npts];
            a[7] = c[ind + 19 * npts];
            a[8] = c[ind + 20 * npts];
            sw4_type dim = 3, one = 1, info = 0, ipiv[3];
            F77_FUNC(dgesv, DGESV)(&dim, &one, a, &dim, ipiv, x, &dim, &info);
            if (info != 0)
              cout << "ERROR in bcfortanisosg_ci, call to DGESV returned info "
                   << info << endl;
            u[ind - nij] = -s0i * x[0];
            u[npts + ind - nij] = -s0i * x[1];
            u[2 * npts + ind - nij] = -s0i * x[2];
          }
      } else {
        sw4_type k = nz;
#pragma omp parallel for
        for (sw4_type j = jb + 2; j <= je - 2; j++)
          for (sw4_type i = ib + 2; i <= ie - 2; i++) {
            size_t qq = i - ib + ni * (j - jb);
            size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
            float_sw4 du = strx[i - ib] * (d4a * (u[ind + 1] - u[ind - 1]) +
                                           d4b * (u[ind + 2] - u[ind - 2]));
            float_sw4 dv =
                strx[i - ib] * (d4a * (u[npts + ind + 1] - u[npts + ind - 1]) +
                                d4b * (u[npts + ind + 2] - u[npts + ind - 2]));
            float_sw4 dw =
                strx[i - ib] *
                (d4a * (u[2 * npts + ind + 1] - u[2 * npts + ind - 1]) +
                 d4b * (u[2 * npts + ind + 2] - u[2 * npts + ind - 2]));

            float_sw4 rhs1 = c[ind + 2 * npts] * du + c[ind + 7 * npts] * dv +
                             c[ind + 11 * npts] * dw;
            float_sw4 rhs2 = c[ind + 4 * npts] * du + c[ind + 9 * npts] * dv +
                             c[ind + 13 * npts] * dw;
            float_sw4 rhs3 = c[ind + 5 * npts] * du + c[ind + 10 * npts] * dv +
                             c[ind + 14 * npts] * dw;

            du = stry[j - jb] * (d4a * (u[ind + ni] - u[ind - ni]) +
                                 d4b * (u[ind + 2 * ni] - u[ind - 2 * ni]));
            dv = stry[j - jb] *
                 (d4a * (u[npts + ind + ni] - u[npts + ind - ni]) +
                  d4b * (u[npts + ind + 2 * ni] - u[npts + ind - 2 * ni]));
            dw = stry[j - jb] *
                 (d4a * (u[2 * npts + ind + ni] - u[2 * npts + ind - ni]) +
                  d4b * (u[2 * npts + ind + 2 * ni] -
                         u[2 * npts + ind - 2 * ni]));
            rhs1 += c[ind + 7 * npts] * du + c[ind + 12 * npts] * dv +
                    c[ind + 13 * npts] * dw;
            rhs2 += c[ind + 9 * npts] * du + c[ind + 16 * npts] * dv +
                    c[ind + 18 * npts] * dw;
            rhs3 += c[ind + 10 * npts] * du + c[ind + 17 * npts] * dv +
                    c[ind + 19 * npts] * dw;

            du = dv = dw = 0;
            for (sw4_type w = 1; w <= 4; w++) {
              du -= sbop[w] * u[ind - nij * (w - 1)];
              dv -= sbop[w] * u[npts + ind - nij * (w - 1)];
              dw -= sbop[w] * u[2 * npts + ind - nij * (w - 1)];
            }
            rhs1 += c[ind + 11 * npts] * du + c[ind + 13 * npts] * dv +
                    c[ind + 14 * npts] * dw - h * bforce6[3 * qq];
            rhs2 += c[ind + 13 * npts] * du + c[ind + 18 * npts] * dv +
                    c[ind + 19 * npts] * dw - h * bforce6[1 + 3 * qq];
            rhs3 += c[ind + 14 * npts] * du + c[ind + 19 * npts] * dv +
                    c[ind + 20 * npts] * dw - h * bforce6[2 + 3 * qq];

            // Solve system for ghost point values
            float_sw4 x[3] = {rhs1, rhs2, rhs3};
            float_sw4 a[9];
            a[0] = c[ind + 11 * npts];
            a[1] = c[ind + 13 * npts];
            a[2] = c[ind + 14 * npts];
            a[3] = c[ind + 13 * npts];
            a[4] = c[ind + 18 * npts];
            a[5] = c[ind + 19 * npts];
            a[6] = c[ind + 14 * npts];
            a[7] = c[ind + 19 * npts];
            a[8] = c[ind + 20 * npts];
            sw4_type dim = 3, one = 1, info = 0, ipiv[3];
            F77_FUNC(dgesv, DGESV)(&dim, &one, a, &dim, ipiv, x, &dim, &info);
            if (info != 0)
              cout << "ERROR in bcfortanisosg_ci, call to DGESV returned info "
                   << info << endl;
            u[ind + nij] = s0i * x[0];
            u[npts + ind + nij] = s0i * x[1];
            u[2 * npts + ind + nij] = s0i * x[2];
          }
      }
    }
  }
}
