
#include "EW.h"
#include "GridGeneratorGaussianHill.h"

GridGeneratorGaussianHill::GridGeneratorGaussianHill(
    float_sw4 topo_zmax, bool always_new, bool analytic_metric,
    int grid_interpolation_order, float_sw4 zetaBreak, float_sw4 amp,
    float_sw4 xc, float_sw4 yc, float_sw4 lx, float_sw4 ly)
    : GridGenerator(topo_zmax, always_new, grid_interpolation_order, zetaBreak),
      m_analytic_metric(analytic_metric),
      m_amp(amp),
      m_xc(xc),
      m_yc(yc),
      m_lx(lx),
      m_ly(ly) {
  m_ixl2 = 1 / (m_lx * m_lx);
  m_iyl2 = 1 / (m_ly * m_ly);
}

//-----------------------------------------------------------------------
bool GridGeneratorGaussianHill::grid_mapping(EW* a_ew, float_sw4 q, float_sw4 r,
                                             float_sw4 s, int g, float_sw4& x,
                                             float_sw4& y, float_sw4& z) {
  float_sw4 h = a_ew->mGridSize[g];
  x = (q - 1) * h;
  y = (r - 1) * h;
  float_sw4 tau = top(x, y);
  if (m_always_new ||
      a_ew->mNumberOfGrids - a_ew->mNumberOfCartesianGrids > 1) {
    // new
    float_sw4 s1 = curvilinear_interface_parameter(
        a_ew, g - a_ew->mNumberOfCartesianGrids);
    float_sw4 s0 = curvilinear_interface_parameter(
        a_ew, g - a_ew->mNumberOfCartesianGrids - 1);
    float_sw4 Ztop = s1 * (-tau) + (1 - s1) * m_topo_zmax;
    float_sw4 Zbot = s0 * (-tau) + (1 - s0) * m_topo_zmax;

    float_sw4 Nz_real =
        static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
    float_sw4 inzr = 1.0 / Nz_real;
    float_sw4 zeta = static_cast<float_sw4>((s - a_ew->m_kStartInt[g]) * inzr);
    z = (1.0 - zeta) * Ztop + zeta * Zbot;
    return true;
  } else {
    // old
    int nz = a_ew->m_global_nz[g];
    float_sw4 zu1 = m_topo_zmax - (nz - 1) * h;
    float_sw4 izb = 1.0 / (m_zetaBreak * (nz - 1));
    z = m_topo_zmax - (nz - s) * h;
    float_sw4 zeta = (s - 1) * izb;
    if (zeta < 1) {
      float_sw4 omsm = (1 - zeta);
      for (int l = 2; l <= m_grid_interpolation_order; l++) omsm *= (1 - zeta);
      z -= omsm * (zu1 + tau);
    }
    return true;
  }
}

//-----------------------------------------------------------------------
bool GridGeneratorGaussianHill::inverse_grid_mapping(EW* a_ew, float_sw4 x,
                                                     float_sw4 y, float_sw4 z,
                                                     int g, float_sw4& q,
                                                     float_sw4& r,
                                                     float_sw4& s) {
  float_sw4 h = a_ew->mGridSize[g];
  q = x / h + 1;
  r = y / h + 1;
  int i = static_cast<int>(round(q));
  int j = static_cast<int>(round(r));
  if (a_ew->interior_point_in_proc(i, j, g)) {
    float_sw4 tau = top(x, y);
    if (m_always_new ||
        a_ew->mNumberOfGrids - a_ew->mNumberOfCartesianGrids > 1) {
      // new
      float_sw4 s1 = curvilinear_interface_parameter(
          a_ew, g - a_ew->mNumberOfCartesianGrids);
      float_sw4 s0 = curvilinear_interface_parameter(
          a_ew, g - a_ew->mNumberOfCartesianGrids - 1);
      float_sw4 Ztop = s1 * (-tau) + (1 - s1) * m_topo_zmax;
      float_sw4 Zbot = s0 * (-tau) + (1 - s0) * m_topo_zmax;

      if (Ztop <= z && z <= Zbot || (g == a_ew->mNumberOfGrids - 1 &&
                                     (Ztop - h * 0.5 <= z && z <= Zbot))) {
        // Point is found on grid:
        float_sw4 Nz_real =
            static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
        float_sw4 zeta = (z - Ztop) / (Zbot - Ztop);
        s = (Nz_real)*zeta + a_ew->m_kStartInt[g];
        //         s = (z-ztop)/(zbot-ztop)*(Nz-1)+1;
        //         if( s < 1 )
        //            s = 1;
        return true;
      } else
        return false;
    } else {
      // old
      int nz = a_ew->m_global_nz[g];
      //   If z is in the Cartesian part of grid, this is the s value:
      s = (z - m_topo_zmax) / h + nz;
      float_sw4 zlim = m_topo_zmax - (nz - 1) * (1 - m_zetaBreak) * h;
      if (z >= zlim)
        return true;
      else {
        // z is in the curvilinear part, solve non-linear equation
        float_sw4 z0 = m_topo_zmax - (nz - 1) * h + tau;
        float_sw4 izb = 1.0 / (m_zetaBreak * (nz - 1));
        float_sw4 tol = 1e-12;
        float_sw4 er = tol + 1;
        int maxit = 10;
        int it = 0;
        while (er > tol && it < maxit) {
          float_sw4 omra = 1 - (s - 1) * izb;
          float_sw4 omsm = omra;
          for (int l = 2; l <= m_grid_interpolation_order - 1; l++)
            omsm *= omra;
          float_sw4 dfcn = h + izb * m_grid_interpolation_order * omsm * z0;
          omsm *= omra;
          float_sw4 fcn = m_topo_zmax - (nz - s) * h - omsm * z0 - z;
          float_sw4 sp = s - fcn / dfcn;
          er = abs(sp - s);
          s = sp;
          it++;
        }
        if (er > tol) {
          cout << "GaussianHill::inverse__grid_mapping: poor convergence for "
                  "X0, Y0, Z0 = "
               << x << ", " << y << ", " << z << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
      return true;
    }
  }
}

//-----------------------------------------------------------------------
void GridGeneratorGaussianHill::grid_mapping_diff(
    EW* a_ew, float_sw4 q, float_sw4 r, float_sw4 s, int g, int ic, int jc,
    int kc, float_sw4& zq, float_sw4& zr, float_sw4& zs, float_sw4& zqq,
    float_sw4& zqr, float_sw4& zqs, float_sw4& zrr, float_sw4& zrs,
    float_sw4& zss) {
  float_sw4 h = a_ew->mGridSize[g];
  //   x = (q-1)*h;
  //   y = (r-1)*h;
  //   float_sw4 tau = top(x,y);
  if (m_always_new ||
      a_ew->mNumberOfGrids - a_ew->mNumberOfCartesianGrids > 1) {
    // new
    float_sw4 s1 = curvilinear_interface_parameter(
        a_ew, g - a_ew->mNumberOfCartesianGrids);
    float_sw4 s0 = curvilinear_interface_parameter(
        a_ew, g - a_ew->mNumberOfCartesianGrids - 1);

    // Topography
    float_sw4 tau, tauq, taur, tauqq, tauqr, taurr;
    topder((q - 1) * h, (r - 1) * h, tau, tauq, taur);
    topder2nd((q - 1) * h, (r - 1) * h, tauqq, tauqr, taurr);
    tauq *= h;
    taur *= h;
    tauqq *= h * h;
    tauqr *= h * h;
    taurr *= h * h;

    // Upper and lower interfaces for this grid
    float_sw4 Ztop = s1 * (-tau) + (1 - s1) * m_topo_zmax;
    float_sw4 Ztopq = s1 * (-tauq);
    float_sw4 Ztopr = s1 * (-taur);
    float_sw4 Ztopqq = s1 * (-tauqq);
    float_sw4 Ztopqr = s1 * (-tauqr);
    float_sw4 Ztoprr = s1 * (-taurr);

    float_sw4 Zbot = s0 * (-tau) + (1 - s0) * m_topo_zmax;
    float_sw4 Zbotq = s0 * (-tauq);
    float_sw4 Zbotr = s0 * (-taur);
    float_sw4 Zbotqq = s0 * (-tauqq);
    float_sw4 Zbotqr = s0 * (-tauqr);
    float_sw4 Zbotrr = s0 * (-taurr);

    float_sw4 Nz_real =
        static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
    float_sw4 inzr = 1.0 / Nz_real;
    float_sw4 zeta = static_cast<float_sw4>((s - a_ew->m_kStartInt[g]) * inzr);

    zq = (1.0 - zeta) * Ztopq + zeta * Zbotq;
    zr = (1.0 - zeta) * Ztopr + zeta * Zbotr;
    zs = (Zbot - Ztop) * inzr;
    zqq = (1.0 - zeta) * Ztopqq + zeta * Zbotqq;
    zqr = (1.0 - zeta) * Ztopqr + zeta * Zbotqr;
    zrr = (1.0 - zeta) * Ztoprr + zeta * Zbotrr;
    zqs = (Zbotq - Ztopq) * inzr;
    zrs = (Zbotr - Ztopr) * inzr;
    zss = 0;
  } else {
    // old
    int nz = a_ew->m_global_nz[g];
    float_sw4 zu1 = m_topo_zmax - (nz - 1) * h;
    float_sw4 izb = 1.0 / (m_zetaBreak * (nz - 1));
    //      z = m_topo_zmax - (nz-s)*h;
    float_sw4 zeta = (s - 1) * izb;
    zq = zr = 0;
    zs = h * (nz - 1);
    if (zeta < 1) {
      int m = m_grid_interpolation_order;  // shorter name
      float_sw4 tau, tauq, taur, tauqq, tauqr, taurr;
      topder((q - 1) * h, (r - 1) * h, tau, tauq, taur);
      topder2nd((q - 1) * h, (r - 1) * h, tauqq, tauqr, taurr);
      float_sw4 inzm1 = 1.0 / (nz - 1);
      float_sw4 omsm = (1 - zeta);
      for (int l = 2; l <= m - 2; l++) omsm *= (1 - zeta);
      zss -= (m / m_zetaBreak) * ((m - 1) / m_zetaBreak) * omsm * (zu1 + tau);
      zss *= inzm1 * inzm1;
      omsm *= (1 - zeta);
      zs += (m / m_zetaBreak) * omsm * (zu1 + tau);
      zs *= inzm1;
      omsm *= (1 - zeta);

      zq = omsm * (-tauq) * h;
      zr = omsm * (-taur) * h;
      zqq = omsm * (-tauqq) * h * h;
      zqr = omsm * (-tauqr) * h * h;
      zrr = omsm * (-taurr) * h * h;
      //         z -=  omsm*(zu1 + tau);
    }
  }
}

//-----------------------------------------------------------------------
void GridGeneratorGaussianHill::generate_grid_and_met(
    EW* a_ew, int g, Sarray& a_x, Sarray& a_y, Sarray& a_z, Sarray& a_jac,
    Sarray& a_met, bool a_comm) {
  int ncurv = a_ew->mNumberOfGrids - a_ew->mNumberOfCartesianGrids;
  if (m_always_new || ncurv > 1)
    generate_grid_and_met_new_gh(a_ew, g, a_x, a_y, a_z, a_jac, a_met);
  else
    generate_grid_and_met_old_gh(a_ew, a_x, a_y, a_z, a_jac, a_met);
  if (!m_analytic_metric) {
    int ierr = a_ew->metric_ci(a_x.m_ib, a_x.m_ie, a_x.m_jb, a_x.m_je, a_x.m_kb,
                               a_x.m_ke, a_x.c_ptr(), a_y.c_ptr(), a_z.c_ptr(),
                               a_met.c_ptr(), a_jac.c_ptr());
    if (a_comm) {
      a_ew->communicate_array(a_jac, g);
      a_ew->communicate_array(a_met, g);
    }
  }
}

//-----------------------------------------------------------------------
void GridGeneratorGaussianHill::generate_grid_and_met_new_gh(
    EW* a_ew, int g, Sarray& a_x, Sarray& a_y, Sarray& a_z, Sarray& a_jac,
    Sarray& a_met) {
  int iSurfTop = g - a_ew->mNumberOfCartesianGrids;
  if (0 <= iSurfTop &&
      iSurfTop <= a_ew->mNumberOfGrids - a_ew->mNumberOfCartesianGrids - 1) {
    float_sw4 h = a_ew->mGridSize[g];
    float_sw4 s1 = curvilinear_interface_parameter(a_ew, iSurfTop);
    float_sw4 s0 = curvilinear_interface_parameter(a_ew, iSurfTop - 1);
    //      std::cout << " grid " << g << " interface parameters s0,s1 " << s0
    //      << " " << s1 <<
    //         " kint = " << a_ew->m_kStartInt[g] << " " << a_ew->m_kEndInt[g]
    //         << std::endl;
    //      std::cout << " array dims " << a_x.m_ib << " " << a_x.m_ie << " " <<
    //      a_x.m_jb << " " << a_x.m_je << " "
    //                << a_x.m_kb << " " << a_x.m_ke << std::endl;
    for (int j = a_x.m_jb; j <= a_x.m_je; j++)
      for (int i = a_x.m_ib; i <= a_x.m_ie; i++) {
        float_sw4 X0 = (i - 1) * h;
        float_sw4 Y0 = (j - 1) * h;

        // Topography
        float_sw4 tau, taup, tauq;
        topder(X0, Y0, tau, taup, tauq);

        // Upper and lower interfaces for this grid
        float_sw4 Ztop = s1 * (-tau) + (1 - s1) * m_topo_zmax;
        float_sw4 Ztopp = s1 * (-taup);
        float_sw4 Ztopq = s1 * (-tauq);

        float_sw4 Zbot = s0 * (-tau) + (1 - s0) * m_topo_zmax;
        float_sw4 Zbotp = s0 * (-taup);
        float_sw4 Zbotq = s0 * (-tauq);

        // Linear interpolation in the vertical direction
        float_sw4 Nz_real =
            static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
        float_sw4 inzr = 1.0 / Nz_real;
#pragma omp parallel for
        for (int k = a_x.m_kb; k <= a_x.m_ke; k++) {
          float_sw4 zeta =
              static_cast<float_sw4>((k - a_ew->m_kStartInt[g]) * inzr);
          a_x(i, j, k) = X0;
          a_y(i, j, k) = Y0;
          a_z(i, j, k) = (1.0 - zeta) * Ztop + zeta * Zbot;

          float_sw4 zp = (1.0 - zeta) * Ztopp + zeta * Zbotp;
          float_sw4 zq = (1.0 - zeta) * Ztopq + zeta * Zbotq;
          float_sw4 zr = Zbot - Ztop;
          zp = zp * h;
          zq = zq * h;
          zr = zr * inzr;
          float_sw4 sqzr = sqrt(zr);
          a_jac(i, j, k) = h * h * zr;
          a_met(1, i, j, k) = sqzr;
          a_met(2, i, j, k) = -zp / sqzr;
          a_met(3, i, j, k) = -zq / sqzr;
          a_met(4, i, j, k) = h / sqzr;
        }
      }
  }
}

//-----------------------------------------------------------------------
void GridGeneratorGaussianHill::generate_grid_and_met_old_gh(
    EW* a_ew, Sarray& a_x, Sarray& a_y, Sarray& a_z, Sarray& a_jac,
    Sarray& a_met) {
  int g = a_ew->mNumberOfGrids - 1;
  float_sw4 h = a_ew->mGridSize[g];
  int nz = a_ew->m_global_nz[g];
  float_sw4 zu1 = m_topo_zmax - (nz - 1) * h;
  float_sw4 inzm1 = 1.0 / (nz - 1);
  float_sw4 izb = 1.0 / (m_zetaBreak * (nz - 1));
  int m = m_grid_interpolation_order;  // shorter name
  for (int k = a_x.m_kb; k <= a_x.m_ke; k++)
    for (int j = a_x.m_jb; j <= a_x.m_je; j++)
      for (int i = a_x.m_ib; i <= a_x.m_ie; i++) {
        // Topography
        float_sw4 tau, taup, tauq;
        topder((i - 1) * h, (j - 1) * h, tau, taup, tauq);

        // 1. Generate grid and its derivatives
        a_x(i, j, k) = (i - 1) * h;
        a_y(i, j, k) = (j - 1) * h;
        a_z(i, j, k) = m_topo_zmax - (nz - k) * h;
        float_sw4 s = (k - 1) * izb;
        float_sw4 zp = 0, zq = 0, zr = h * (nz - 1);
        if (s < 1) {
          float_sw4 omsm = (1 - s);
          for (int l = 2; l <= m - 1; l++) omsm *= (1 - s);
          zr += (m / m_zetaBreak) * omsm * (zu1 + tau);
          omsm *= (1 - s);
          a_z(i, j, k) -= omsm * (zu1 + tau);
          zp = omsm * (-taup);
          zq = omsm * (-tauq);
        }
        zr = zr * inzm1;

        // 2. Check if ok
        //	    REQUIRE2( zr >=0 , "Error, zr = " << zr << " at " << i << "
        //"
        //		      <<  j << " " <<  k << std::endl);
        if (zr < 0) {
          std::cout << "Error, zr = " << zr << " at " << i << " " << j << " "
                    << k << std::endl;
          std::cout << " s= " << s << " zu1, tau, m, m_zetaBreak" << zu1 << " "
                    << tau << " " << m << " " << m_zetaBreak << std::endl;
          exit(0);
          return;
        }

        // 3. Compute metric and Jacobian
        float_sw4 sqzr = sqrt(zr);
        a_jac(i, j, k) = h * h * zr;
        a_met(1, i, j, k) = sqzr;
        a_met(2, i, j, k) = -zp / sqzr;
        a_met(3, i, j, k) = -zq / sqzr;
        a_met(4, i, j, k) = h / sqzr;
      }
}

//-----------------------------------------------------------------------
bool GridGeneratorGaussianHill::interpolate_topography(EW* a_ew, float_sw4 x,
                                                       float_sw4 y,
                                                       float_sw4& z,
                                                       Sarray& topo) {
  //   float_sw4 h = a_ew->mGridSize[a_ew->mNumberOfGrids-1];
  //   float_sw4 x = (q-1)*h;
  //   float_sw4 y = (r-1)*h;
  z = -m_amp *
      exp(-(x - m_xc) * (x - m_xc) * m_ixl2 - (y - m_yc) * (y - m_yc) * m_iyl2);
  return true;
}

//-----------------------------------------------------------------------
bool GridGeneratorGaussianHill::exact_metric(EW* a_ew, int g, Sarray& a_jac,
                                             Sarray& a_met) {
  Sarray x(a_jac), y(a_jac), z(a_jac);
  int ncurv = a_ew->mNumberOfGrids - a_ew->mNumberOfCartesianGrids;
  if (m_always_new || ncurv > 1)
    generate_grid_and_met_new_gh(a_ew, g, x, y, z, a_jac, a_met);
  else
    generate_grid_and_met_old_gh(a_ew, x, y, z, a_jac, a_met);
  return true;
}

//-----------------------------------------------------------------------
void GridGeneratorGaussianHill::fill_topo(Sarray& topo, float_sw4 h) {
  for (int i = topo.m_ib; i <= topo.m_ie; ++i)
    for (int j = topo.m_jb; j <= topo.m_je; ++j) {
      float_sw4 x = (i - 1) * h;
      float_sw4 y = (j - 1) * h;
      // positive topography  is up (negative z)
      topo(i, j, 1) = m_amp * exp(-(x - m_xc) * (x - m_xc) * m_ixl2 -
                                  (y - m_yc) * (y - m_yc) * m_iyl2);
    }
}

//-----------------------------------------------------------------------
void GridGeneratorGaussianHill::generate_z_and_j(EW* a_ew, int g, Sarray& z,
                                                 Sarray& J) {
  int iSurfTop = g - a_ew->mNumberOfCartesianGrids;
  if (0 <= iSurfTop &&
      iSurfTop <= a_ew->mNumberOfGrids - a_ew->mNumberOfCartesianGrids - 1) {
    float_sw4 h = a_ew->mGridSize[g];
    float_sw4 s1 = curvilinear_interface_parameter(a_ew, iSurfTop);
    float_sw4 s0 = curvilinear_interface_parameter(a_ew, iSurfTop - 1);
    for (int j = z.m_jb; j <= z.m_je; j++)
      for (int i = z.m_ib; i <= z.m_ie; i++) {
        float_sw4 X0 = (i - 1) * h;
        float_sw4 Y0 = (j - 1) * h;

        // Topography
        float_sw4 tau, taup, tauq;
        topder(X0, Y0, tau, taup, tauq);

        // Upper and lower interfaces for this grid
        float_sw4 Ztop = s1 * (-tau) + (1 - s1) * m_topo_zmax;
        float_sw4 Ztopp = s1 * (-taup);
        float_sw4 Ztopq = s1 * (-tauq);

        float_sw4 Zbot = s0 * (-tau) + (1 - s0) * m_topo_zmax;
        float_sw4 Zbotp = s0 * (-taup);
        float_sw4 Zbotq = s0 * (-tauq);

        // Linear interpolation in the vertical direction
        float_sw4 Nz_real =
            static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
        float_sw4 inzr = 1.0 / Nz_real;
#pragma omp parallel for
        for (int k = z.m_kb; k <= z.m_ke; k++) {
          float_sw4 zeta =
              static_cast<float_sw4>((k - a_ew->m_kStartInt[g]) * inzr);
          z(i, j, k) = (1.0 - zeta) * Ztop + zeta * Zbot;
          J(i, j, k) = h * h * (Zbot - Ztop) * inzr;
        }
      }
  }
}
