#ifndef SW4_TESTGRID_H
#define SW4_TESTGRID_H

#include <cmath>
#include "sw4.h"

class EW;
class Sarray;

class TestGrid {
  // Gaussian hill grid for testing purpose
  float_sw4 m_amp, m_xc, m_yc, m_lx, m_ly, m_ixl2, m_iyl2, m_topo_zmax;

  // Describing the topography
  inline float_sw4 top(float_sw4 x, float_sw4 y) {
    return m_amp * exp(-(x - m_xc) * (x - m_xc) * m_ixl2 -
                       (y - m_yc) * (y - m_yc) * m_iyl2);
  }
  inline void topder(float_sw4 x, float_sw4 y, float_sw4& tau, float_sw4& taup,
                     float_sw4& tauq) {
    tau = m_amp * exp(-(x - m_xc) * (x - m_xc) * m_ixl2 -
                      (y - m_yc) * (y - m_yc) * m_iyl2);
    taup = -2 * (x - m_xc) * m_ixl2 * tau;
    tauq = -2 * (y - m_yc) * m_iyl2 * tau;
  }

 public:
  TestGrid(float_sw4 topo_zmax, float_sw4 amp, float_sw4 xc, float_sw4 yc,
           float_sw4 lx, float_sw4 ly);
  void generate_grid_and_met(EW* a_ew, int g, Sarray& a_x, Sarray& a_y,
                             Sarray& a_z, Sarray& a_jac, Sarray& a_met);
};
#endif
