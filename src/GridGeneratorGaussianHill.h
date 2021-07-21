#ifndef GRIDGENERATORGAUSSIANHILL_H
#define GRIDGENERATORGAUSSIANHILL_H

#include "GridGenerator.h"

class GridGeneratorGaussianHill : public GridGenerator {
  bool m_analytic_metric;
  float_sw4 m_amp, m_xc, m_yc, m_lx, m_ly, m_ixl2, m_iyl2;

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
  inline void topder2nd(float_sw4 x, float_sw4 y, float_sw4& taupp,
                        float_sw4& taupq, float_sw4& tauqq) {
    float_sw4 tau = m_amp * exp(-(x - m_xc) * (x - m_xc) * m_ixl2 -
                                (y - m_yc) * (y - m_yc) * m_iyl2);
    taupp =
        -2 * m_ixl2 * tau + 4 * (x - m_xc) * (x - m_xc) * m_ixl2 * m_ixl2 * tau;
    taupq = 4 * (x - m_xc) * (y - m_yc) * m_ixl2 * m_iyl2 * tau;
    tauqq =
        -2 * m_iyl2 * tau + 4 * (y - m_yc) * (y - m_yc) * m_iyl2 * m_iyl2 * tau;
  }
  void generate_grid_and_met_new_gh(EW* a_ew, int g, Sarray& a_x, Sarray& a_y,
                                    Sarray& a_z, Sarray& a_jac, Sarray& a_met);
  void generate_grid_and_met_old_gh(EW* a_ew, Sarray& a_x, Sarray& a_y,
                                    Sarray& a_z, Sarray& a_jac, Sarray& a_met);

 public:
  GridGeneratorGaussianHill(float_sw4 topo_zmax, bool always_new,
                            bool analytic_metric, int grid_interpolation_order,
                            float_sw4 zetaBreak, float_sw4 amp, float_sw4 xc,
                            float_sw4 yc, float_sw4 lx, float_sw4 ly);

  virtual void assignInterfaceSurfaces(EW* a_ew, Sarray& TopoGridExt){};

  virtual void generate_grid_and_met(EW* a_ew, int g, Sarray& a_x, Sarray& a_y,
                                     Sarray& a_z, Sarray& a_jac, Sarray& a_met,
                                     bool a_comm = true);
  virtual bool grid_mapping(EW* a_ew, float_sw4 p, float_sw4 q, float_sw4 r,
                            int g, float_sw4& x, float_sw4& y, float_sw4& z);
  virtual bool inverse_grid_mapping(EW* a_ew, float_sw4 x, float_sw4 y,
                                    float_sw4 z, int g, float_sw4& p,
                                    float_sw4& q, float_sw4& r);
  virtual void grid_mapping_diff(EW* a_ew, float_sw4 q, float_sw4 r,
                                 float_sw4 s, int g, int ic, int jc, int kc,
                                 float_sw4& zq, float_sw4& zr, float_sw4& zs,
                                 float_sw4& zqq, float_sw4& zqr, float_sw4& zqs,
                                 float_sw4& zrr, float_sw4& zrs,
                                 float_sw4& zss);
  virtual int interpolate_topography(EW* a_ew, float_sw4 q, float_sw4 r,
                                      float_sw4& z, Sarray& topo);
  virtual bool exact_metric(EW* a_ew, int g, Sarray& a_jac, Sarray& a_met);
  virtual void fill_topo(Sarray& topo, float_sw4 h);
  virtual void generate_z_and_j(EW* a_ew, int g, Sarray& z, Sarray& J);
};
#endif
