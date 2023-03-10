
#include "EW.h"
#include "GridGenerator.h"

//-----------------------------------------------------------------------
int GridGenerator::metric_ci(int ib, int ie, int jb, int je, int kb, int ke,
                             float_sw4* __restrict__ a_x,
                             float_sw4* __restrict__ a_y,
                             float_sw4* __restrict__ a_z,
                             float_sw4* __restrict__ a_met,
                             float_sw4* __restrict__ a_jac) {
  const float_sw4 c1 = 2.0 / 3, c2 = -1.0 / 12;
  const float_sw4 fs = 5.0 / 6, ot = 1.0 / 12, ft = 4.0 / 3, os = 1.0 / 6,
                  d3 = 14.0 / 3;
  const int ni = ie - ib + 1;
  const int nij = ni * (je - jb + 1);
  const int nijk = nij * (ke - kb + 1);
  const int base = -(ib + ni * jb + nij * kb);
  const int base4 = base - nijk;
  int ecode = 0;
  //   const int nic  = 4*ni;
  //   const int nijc = 4*nij;
#define x(i, j, k) a_x[base + i + ni * (j) + nij * (k)]
#define y(i, j, k) a_y[base + i + ni * (j) + nij * (k)]
#define z(i, j, k) a_z[base + i + ni * (j) + nij * (k)]
#define jac(i, j, k) a_jac[base + i + ni * (j) + nij * (k)]
#define met(c, i, j, k) a_met[base4 + (i) + ni * (j) + nij * (k) + nijk * (c)]

  double h = x(ib + 1, jb, kb) - x(ib, jb, kb);

#pragma omp parallel for reduction(+ : ecode)
  for (int k = kb; k <= ke; k++)
    for (int j = jb; j <= je; j++)
      //#pragma ivdep
      //#pragma omp simd
      for (int i = ib; i <= ie; i++) {
        // k-derivatives
        double zr, zp, zq, sqzr;
        if (k >= kb + 2 && k <= ke - 2)
          zr = c2 * (z(i, j, k + 2) - z(i, j, k - 2)) +
               c1 * (z(i, j, k + 1) - z(i, j, k - 1));
        else if (k == kb) {
          zr = -2.25 * z(i, j, k) + (4 + fs) * z(i, j, k + 1) -
               d3 * z(i, j, k + 2) + 3 * z(i, j, k + 3) -
               (1 + ot) * z(i, j, k + 4) + os * z(i, j, k + 5);
        } else if (k == kb + 1) {
          zr = -os * z(i, j, k - 1) - 1.25 * z(i, j, k) +
               (1 + ft) * z(i, j, k + 1) - ft * z(i, j, k + 2) +
               0.5 * z(i, j, k + 3) - ot * z(i, j, k + 4);
        } else if (k == ke - 1) {
          zr = os * z(i, j, k + 1) + 1.25 * z(i, j, k) -
               (1 + ft) * z(i, j, k - 1) + ft * z(i, j, k - 2) -
               0.5 * z(i, j, k - 3) + ot * z(i, j, k - 4);
        } else if (k == ke) {
          zr = 2.25 * z(i, j, k) - (4 + fs) * z(i, j, k - 1) +
               d3 * z(i, j, k - 2) - 3 * z(i, j, k - 3) +
               (1 + ot) * z(i, j, k - 4) - os * z(i, j, k - 5);
        }

        // j-derivatives
        if (j >= jb + 2 && j <= je - 2) {
          zq = c2 * (z(i, j + 2, k) - z(i, j - 2, k)) +
               c1 * (z(i, j + 1, k) - z(i, j - 1, k));
        } else if (j == jb) {
          zq = -2.25 * z(i, j, k) + (4 + fs) * z(i, j + 1, k) -
               d3 * z(i, j + 2, k) + 3 * z(i, j + 3, k) -
               (1 + ot) * z(i, j + 4, k) + os * z(i, j + 5, k);
        } else if (j == jb + 1) {
          zq = -os * z(i, j - 1, k) - 1.25 * z(i, j, k) +
               (1 + ft) * z(i, j + 1, k) - ft * z(i, j + 2, k) +
               0.5 * z(i, j + 3, k) - ot * z(i, j + 4, k);
        } else if (j == je - 1) {
          zq = os * z(i, j + 1, k) + 1.25 * z(i, j, k) -
               (1 + ft) * z(i, j - 1, k) + ft * z(i, j - 2, k) -
               0.5 * z(i, j - 3, k) + ot * z(i, j - 4, k);
        } else if (j == je) {
          zq = 2.25 * z(i, j, k) - (4 + fs) * z(i, j - 1, k) +
               d3 * z(i, j - 2, k) - 3 * z(i, j - 3, k) +
               (1 + ot) * z(i, j - 4, k) - os * z(i, j - 5, k);
        }

        // i-derivatives
        if (i >= ib + 2 && i <= ie - 2) {
          zp = c2 * (z(i + 2, j, k) - z(i - 2, j, k)) +
               c1 * (z(i + 1, j, k) - z(i - 1, j, k));
        } else if (i == ib) {
          zp = -2.25 * z(i, j, k) + (4 + fs) * z(i + 1, j, k) -
               d3 * z(i + 2, j, k) + 3 * z(i + 3, j, k) -
               (1 + ot) * z(i + 4, j, k) + os * z(i + 5, j, k);
        } else if (i == ib + 1) {
          zp = -os * z(i - 1, j, k) - 1.25 * z(i, j, k) +
               (1 + ft) * z(i + 1, j, k) - ft * z(i + 2, j, k) +
               0.5 * z(i + 3, j, k) - ot * z(i + 4, j, k);
        } else if (i == ie - 1) {
          zp = os * z(i + 1, j, k) + 1.25 * z(i, j, k) -
               (1 + ft) * z(i - 1, j, k) + ft * z(i - 2, j, k) -
               0.5 * z(i - 3, j, k) + ot * z(i - 4, j, k);
        } else if (i == ie) {
          zp = 2.25 * z(i, j, k) - (4 + fs) * z(i - 1, j, k) +
               d3 * z(i - 2, j, k) - 3 * z(i - 3, j, k) +
               (1 + ot) * z(i - 4, j, k) - os * z(i - 5, j, k);
        }

        // Compute the metric
        if (zr <= 0) {
          ecode = -1;
          // tmp
          cout << "zr = " << zr << " at " << i << " " << j << " " << k << endl;
          cout << "x,y,z = " << x(i, j, k) << " " << y(i, j, k) << " "
               << z(i, j, k) << endl;
          //               printf("ib=%d, ie=%d, jb=%d, je=%d, kb=%d, ke=%d,
          //               proc=%d\n", ib, ie, jb, je, kb, ke, m_myRank);
          printf("ib=%d, ie=%d, jb=%d, je=%d, kb=%d, ke=%d \n", ib, ie, jb, je,
                 kb, ke);
          // end tmp
        }
        sqzr = sqrt(zr);
        jac(i, j, k) = h * h * zr;
        met(1, i, j, k) = sqzr;
        met(2, i, j, k) = -zp / sqzr;
        met(3, i, j, k) = -zq / sqzr;
        met(4, i, j, k) = h / sqzr;
      }
  return ecode;
#undef x
#undef y
#undef z
#undef jac
#undef met
}

//-----------------------------------------------------------------------
bool GridGenerator::interpolate_topography(EW* a_ew, float_sw4 x, float_sw4 y,
                                           float_sw4& z, Sarray& topo) {
  // Interpolate the topography
  //
  // if (q,r) is on this processor then
  // Return true and assign z corresponding to (q,r)
  //
  // Returns false if
  // 1) (q,r) is outside the global parameter domain (expanded by ghost points)
  // 2) There is no curvilinear grid.
  // 3) (q,r) is not on this processor
  //
  // The parameters are normalized such that 1 <= q <= Nx is the full domain
  // (without ghost points),
  //  1 <= r <= Ny.

  if (!a_ew->topographyExists()) return false;

  int gTop = a_ew->mNumberOfGrids - 1;
  float_sw4 h = a_ew->mGridSize[gTop];
  float_sw4 q = x / h + 1.0;
  float_sw4 r = y / h + 1.0;
  // Find topography at (q,r), tau=tau(q,r)
  // Nearest grid point:
  int iNear = static_cast<int>(round(q));
  int jNear = static_cast<int>(round(r));
  float_sw4 tau;

  if (fabs(iNear - q) < 1.e-9 && fabs(jNear - r) < 1.e-9) {
    // At a grid point, evaluate topography at that point
    if (topo.in_range(1, iNear, jNear, 1))
      tau = topo(iNear, jNear, 1);
    else
      return false;
  } else {
    // Not at a grid  point, interpolate the topography
    // Nearest lower grid point
    iNear = static_cast<int>(floor(q));
    jNear = static_cast<int>(floor(r));
    if (topo.in_range(1, iNear - 3, jNear - 3, 1) &&
        topo.in_range(1, iNear + 4, jNear + 4, 1)) {
      float_sw4 a6cofi[8], a6cofj[8];
      gettopowgh(q - iNear, a6cofi);
      gettopowgh(r - jNear, a6cofj);
      tau = 0;
      for (int l = -3; l <= 4; l++)
        for (int k = -3; k <= 4; k++)
          tau += a6cofi[k + 3] * a6cofj[l + 3] * topo(k + iNear, l + jNear, 1);
    } else {
      return false;
    }
  }
  z = -tau;
  return true;
}

//-----------------------------------------------------------------------
void GridGenerator::gettopowgh(float_sw4 ai, float_sw4 wgh[8]) const {
  float_sw4 pol = ai * ai * ai * ai * ai * ai * ai *
                  (-251 + 135 * ai + 25 * ai * ai - 33 * ai * ai * ai +
                   6 * ai * ai * ai * ai) /
                  720;
  wgh[0] = -1.0 / 60 * ai + 1.0 / 180 * ai * ai + 1.0 / 48 * ai * ai * ai +
           23.0 / 144 * ai * ai * ai * ai -
           (17.0 * ai + 223.0) * ai * ai * ai * ai * ai / 720 - pol;
  wgh[1] = 3.0 / 20 * ai - 3.0 / 40 * ai * ai - 1.0 / 6 * ai * ai * ai -
           13.0 / 12 * ai * ai * ai * ai + 97.0 / 45 * ai * ai * ai * ai * ai +
           1.0 / 6 * ai * ai * ai * ai * ai * ai + 7 * pol;
  wgh[2] = -0.75 * ai + 0.75 * ai * ai + (13.0 + 155 * ai) * ai * ai * ai / 48 -
           103.0 / 16 * ai * ai * ai * ai * ai -
           121.0 / 240 * ai * ai * ai * ai * ai * ai - 21 * pol;
  wgh[3] = 1 - 49.0 / 36 * ai * ai - 49.0 / 9 * ai * ai * ai * ai +
           385.0 / 36 * ai * ai * ai * ai * ai +
           61.0 / 72 * ai * ai * ai * ai * ai * ai + 35 * pol;
  wgh[4] = 0.75 * ai + 0.75 * ai * ai - 13.0 / 48 * ai * ai * ai +
           89.0 / 16 * ai * ai * ai * ai -
           1537.0 / 144 * ai * ai * ai * ai * ai -
           41.0 / 48 * ai * ai * ai * ai * ai * ai - 35 * pol;
  wgh[5] = -3.0 / 20 * ai - 3.0 / 40 * ai * ai + 1.0 / 6 * ai * ai * ai -
           41.0 / 12 * ai * ai * ai * ai + 6.4 * ai * ai * ai * ai * ai +
           31.0 / 60 * ai * ai * ai * ai * ai * ai + 21 * pol;
  wgh[6] = 1.0 / 60 * ai + 1.0 / 180 * ai * ai - 1.0 / 48 * ai * ai * ai +
           167.0 / 144 * ai * ai * ai * ai -
           1537.0 / 720 * ai * ai * ai * ai * ai -
           25.0 / 144 * ai * ai * ai * ai * ai * ai - 7 * pol;
  wgh[7] = -1.0 / 6 * ai * ai * ai * ai + 11.0 / 36 * ai * ai * ai * ai * ai +
           1.0 / 40 * ai * ai * ai * ai * ai * ai + pol;
}

//-----------------------------------------------------------------------
bool GridGenerator::exact_metric(EW* a_ew, int g, Sarray& a_jac,
                                 Sarray& a_met) {
  std::cout << "GridGenerator: Exact metric is not available \n" << std::endl;
  return true;
}

//-----------------------------------------------------------------------
void GridGenerator::get_gridgen_info(int& order, float_sw4& zetaBreak) const {
  order = m_grid_interpolation_order;
  zetaBreak = m_zetaBreak;
}

//-----------------------------------------------------------------------
float_sw4 GridGenerator::curvilinear_interface_parameter(EW* a_ew, int gcurv) {
  if (gcurv < 0)
    return 0;
  else
    return (m_topo_zmax - a_ew->m_curviRefLev[gcurv]) / m_topo_zmax;
}

//-----------------------------------------------------------------------
void GridGenerator::fill_topo(Sarray& topo, float_sw4 h) {
  std::cout << "GridGenerator: Exact topography is not available \n"
            << std::endl;
}
//-----------------------------------------------------------------------
void GridGenerator::grid_mapping_diff(EW* a_ew, float_sw4 q, float_sw4 r,
                                      float_sw4 s, int g, int ic, int jc,
                                      int kc, float_sw4& zq, float_sw4& zr,
                                      float_sw4& zs, float_sw4& zqq,
                                      float_sw4& zqr, float_sw4& zqs,
                                      float_sw4& zrr, float_sw4& zrs,
                                      float_sw4& zss) {
  std::cout << "GridGenerator: grid_mapping_diff NYI " << std::endl;
}
