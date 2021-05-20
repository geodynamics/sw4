
#include <iostream>
#include "EW.h"
#include "TestGrid.h"

TestGrid::TestGrid(float_sw4 topo_zmax, float_sw4 amp, float_sw4 xc,
                   float_sw4 yc, float_sw4 lx, float_sw4 ly)
    : m_topo_zmax(topo_zmax),
      m_amp(amp),
      m_xc(xc),
      m_yc(yc),
      m_lx(lx),
      m_ly(ly) {
  m_ixl2 = 1 / (m_lx * m_lx);
  m_iyl2 = 1 / (m_ly * m_ly);
}

void TestGrid::generate_grid_and_met(EW* a_ew, int g, Sarray& a_x, Sarray& a_y,
                                     Sarray& a_z, Sarray& a_jac,
                                     Sarray& a_met) {
  int iSurfTop = g - a_ew->mNumberOfCartesianGrids;
  if (0 <= iSurfTop &&
      iSurfTop <= a_ew->mNumberOfGrids - a_ew->mNumberOfCartesianGrids - 1) {
    float_sw4 h = a_ew->mGridSize[g];
    float_sw4 s1 = a_ew->curvilinear_interface_parameter(iSurfTop);
    float_sw4 s0 = a_ew->curvilinear_interface_parameter(iSurfTop - 1);
    std::cout << " grid " << g << " interface parameters s0,s1 " << s0 << " "
              << s1 << " kint = " << a_ew->m_kStartInt[g] << " "
              << a_ew->m_kEndInt[g] << std::endl;
    std::cout << " array dims " << a_x.m_ib << " " << a_x.m_ie << " "
              << a_x.m_jb << " " << a_x.m_je << " " << a_x.m_kb << " "
              << a_x.m_ke << std::endl;
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
  //   else if( 0 <= g && g < a_ew->mNumberOfCartesianGrids )
  //   {
  //      float_sw4 h  = a_ew->mGridSize[g];
  //      for (int k= a_x.m_kb; k<= a_x.m_ke; k++)
  //         for (int j=a_x.m_jb; j<=a_x.m_je; j++)
  //            for (int i=a_x.m_ib; i<=a_x.m_ie; i++)
  //            {
  //               a_x(i,j,k) = (i-1)*h;
  //               a_y(i,j,k) = (j-1)*h;
  //               a_z(i,j,k) = (k-1)*h+zmin;
  //            }
  //   }
}
