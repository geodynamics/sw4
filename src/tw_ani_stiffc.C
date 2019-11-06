#include "EW.h"
#include "sw4.h"

//-----------------------------------------------------------------------
void EW::tw_ani_stiff_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4 h, float_sw4 zmin,
                         float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                         float_sw4* __restrict__ a_rho, float_sw4 a_phc[21],
                         float_sw4* __restrict__ a_cm)

{
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);

#define phc(c) a_phc[c - 1]
#define rho(i, j, k) a_rho[ind]
#define cm(c, i, j, k) a_cm[ind + (c - 1) * nijk]

#pragma omp parallel
#pragma omp for
  for (int k = kfirst; k <= klast; k++)
    for (int j = jfirst; j <= jlast; j++)
#pragma ivdep
#pragma simd
      for (int i = ifirst; i <= ilast; i++) {
        size_t ind = base + i + ni * j + nij * k;
        float_sw4 x = (i - 1) * h;
        float_sw4 y = (j - 1) * h;
        float_sw4 z = zmin + (k - 1) * h;
        // note that the constant in the diagonal elements is 10 instead of 2
        rho(i, j, k) = amprho * (2 + sin(omm * x + phm) * cos(omm * y + phm) *
                                         sin(omm * z + phm));
        cm(1, i, j, k) = 10 + sin(omm * x + phc(1)) * cos(omm * y + phc(1)) *
                                  sin(omm * z + phc(1));
        cm(2, i, j, k) = 2 + sin(omm * x + phc(2)) * cos(omm * y + phc(1)) *
                                 sin(omm * z + phc(2));
        cm(3, i, j, k) = 2 + sin(omm * x + phc(3)) * cos(omm * y + phc(2)) *
                                 sin(omm * z + phc(3));
        cm(4, i, j, k) = 2 + sin(omm * x + phc(4)) * cos(omm * y + phc(2)) *
                                 sin(omm * z + phc(4));
        cm(5, i, j, k) = 2 + sin(omm * x + phc(5)) * cos(omm * y + phc(3)) *
                                 sin(omm * z + phc(5));
        cm(6, i, j, k) = 2 + sin(omm * x + phc(6)) * cos(omm * y + phc(3)) *
                                 sin(omm * z + phc(6));
        cm(7, i, j, k) = 10 + sin(omm * x + phc(7)) * cos(omm * y + phc(4)) *
                                  sin(omm * z + phc(7));
        cm(8, i, j, k) = 2 + sin(omm * x + phc(8)) * cos(omm * y + phc(4)) *
                                 sin(omm * z + phc(8));
        cm(9, i, j, k) = 2 + sin(omm * x + phc(9)) * cos(omm * y + phc(5)) *
                                 sin(omm * z + phc(9));
        cm(10, i, j, k) = 2 + sin(omm * x + phc(10)) * cos(omm * y + phc(5)) *
                                  sin(omm * z + phc(1));
        cm(11, i, j, k) = 2 + sin(omm * x + phc(11)) * cos(omm * y + phc(6)) *
                                  sin(omm * z + phc(2));
        cm(12, i, j, k) = 10 + sin(omm * x + phc(12)) * cos(omm * y + phc(6)) *
                                   sin(omm * z + phc(3));
        cm(13, i, j, k) = 2 + sin(omm * x + phc(13)) * cos(omm * y + phc(7)) *
                                  sin(omm * z + phc(4));
        cm(14, i, j, k) = 2 + sin(omm * x + phc(14)) * cos(omm * y + phc(7)) *
                                  sin(omm * z + phc(5));
        cm(15, i, j, k) = 2 + sin(omm * x + phc(15)) * cos(omm * y + phc(8)) *
                                  sin(omm * z + phc(6));
        cm(16, i, j, k) = 10 + sin(omm * x + phc(16)) * cos(omm * y + phc(8)) *
                                   sin(omm * z + phc(7));
        cm(17, i, j, k) = 2 + sin(omm * x + phc(17)) * cos(omm * y + phc(9)) *
                                  sin(omm * z + phc(8));
        cm(18, i, j, k) = 2 + sin(omm * x + phc(18)) * cos(omm * y + phc(9)) *
                                  sin(omm * z + phc(9));
        cm(19, i, j, k) = 10 + sin(omm * x + phc(19)) * cos(omm * y + phc(10)) *
                                   sin(omm * z + phc(10));
        cm(20, i, j, k) = 2 + sin(omm * x + phc(20)) * cos(omm * y + phc(10)) *
                                  sin(omm * z + phc(11));
        cm(21, i, j, k) = 10 + sin(omm * x + phc(21)) * cos(omm * y + phc(10)) *
                                   sin(omm * z + phc(12));
      }
#undef rho
#undef cm
#undef phc
}

//-----------------------------------------------------------------------
void EW::tw_ani_curvi_stiff_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
    float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
    float_sw4* __restrict__ zz, float_sw4 omm, float_sw4 phm, float_sw4 amprho,
    float_sw4* __restrict__ a_rho, float_sw4 a_phc[21],
    float_sw4* __restrict__ a_cm)

{
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);

#define phc(c) a_phc[c - 1]
#define rho(i, j, k) a_rho[ind]
#define cm(c, i, j, k) a_cm[ind + (c - 1) * nijk]

#pragma omp parallel
#pragma omp for
  for (int k = kfirst; k <= klast; k++)
    for (int j = jfirst; j <= jlast; j++)
#pragma ivdep
#pragma simd
      for (int i = ifirst; i <= ilast; i++) {
        size_t ind = base + i + ni * j + nij * k;
        float_sw4 x = xx[ind];
        float_sw4 y = yy[ind];
        float_sw4 z = zz[ind];
        // note that the constant in the diagonal elements is 10 instead of 2
        rho(i, j, k) = amprho * (2 + sin(omm * x + phm) * cos(omm * y + phm) *
                                         sin(omm * z + phm));
        cm(1, i, j, k) = 10 + sin(omm * x + phc(1)) * cos(omm * y + phc(1)) *
                                  sin(omm * z + phc(1));
        cm(2, i, j, k) = 2 + sin(omm * x + phc(2)) * cos(omm * y + phc(1)) *
                                 sin(omm * z + phc(2));
        cm(3, i, j, k) = 2 + sin(omm * x + phc(3)) * cos(omm * y + phc(2)) *
                                 sin(omm * z + phc(3));
        cm(4, i, j, k) = 2 + sin(omm * x + phc(4)) * cos(omm * y + phc(2)) *
                                 sin(omm * z + phc(4));
        cm(5, i, j, k) = 2 + sin(omm * x + phc(5)) * cos(omm * y + phc(3)) *
                                 sin(omm * z + phc(5));
        cm(6, i, j, k) = 2 + sin(omm * x + phc(6)) * cos(omm * y + phc(3)) *
                                 sin(omm * z + phc(6));
        cm(7, i, j, k) = 10 + sin(omm * x + phc(7)) * cos(omm * y + phc(4)) *
                                  sin(omm * z + phc(7));
        cm(8, i, j, k) = 2 + sin(omm * x + phc(8)) * cos(omm * y + phc(4)) *
                                 sin(omm * z + phc(8));
        cm(9, i, j, k) = 2 + sin(omm * x + phc(9)) * cos(omm * y + phc(5)) *
                                 sin(omm * z + phc(9));
        cm(10, i, j, k) = 2 + sin(omm * x + phc(10)) * cos(omm * y + phc(5)) *
                                  sin(omm * z + phc(1));
        cm(11, i, j, k) = 2 + sin(omm * x + phc(11)) * cos(omm * y + phc(6)) *
                                  sin(omm * z + phc(2));
        cm(12, i, j, k) = 10 + sin(omm * x + phc(12)) * cos(omm * y + phc(6)) *
                                   sin(omm * z + phc(3));
        cm(13, i, j, k) = 2 + sin(omm * x + phc(13)) * cos(omm * y + phc(7)) *
                                  sin(omm * z + phc(4));
        cm(14, i, j, k) = 2 + sin(omm * x + phc(14)) * cos(omm * y + phc(7)) *
                                  sin(omm * z + phc(5));
        cm(15, i, j, k) = 2 + sin(omm * x + phc(15)) * cos(omm * y + phc(8)) *
                                  sin(omm * z + phc(6));
        cm(16, i, j, k) = 10 + sin(omm * x + phc(16)) * cos(omm * y + phc(8)) *
                                   sin(omm * z + phc(7));
        cm(17, i, j, k) = 2 + sin(omm * x + phc(17)) * cos(omm * y + phc(9)) *
                                  sin(omm * z + phc(8));
        cm(18, i, j, k) = 2 + sin(omm * x + phc(18)) * cos(omm * y + phc(9)) *
                                  sin(omm * z + phc(9));
        cm(19, i, j, k) = 10 + sin(omm * x + phc(19)) * cos(omm * y + phc(10)) *
                                   sin(omm * z + phc(10));
        cm(20, i, j, k) = 2 + sin(omm * x + phc(20)) * cos(omm * y + phc(10)) *
                                  sin(omm * z + phc(11));
        cm(21, i, j, k) = 10 + sin(omm * x + phc(21)) * cos(omm * y + phc(10)) *
                                   sin(omm * z + phc(12));
      }
#undef rho
#undef cm
#undef phc
}
