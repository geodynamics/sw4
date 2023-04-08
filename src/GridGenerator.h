#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H
class Sarray;
class EW;

#include "sw4.h"

class GridGenerator {
 protected:
  float_sw4 m_topo_zmax, m_zetaBreak;
  bool m_always_new;
  sw4_type m_grid_sw4_typeerpolation_order;
  sw4_type metric_ci(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke, float_sw4* a_x,
                float_sw4* a_y, float_sw4* a_z, float_sw4* a_met,
                float_sw4* a_jac);
  void gettopowgh(float_sw4 ai, float_sw4 wgh[8]) const;
  float_sw4 curvilinear_sw4_typeerface_parameter(EW* a_ew, sw4_type gcurv);

 public:
  GridGenerator(float_sw4 topo_zmax, bool always_new,
                sw4_type grid_sw4_typeerpolation_order, float_sw4 zetaBreak)
      : m_topo_zmax(topo_zmax),
        m_zetaBreak(zetaBreak),
        m_always_new(always_new),
        m_grid_sw4_typeerpolation_order(grid_sw4_typeerpolation_order){};
  bool curviCartIsSmooth(sw4_type ncurv) { return ncurv == 1 && !m_always_new; }
  float_sw4 get_topo_zmax() { return m_topo_zmax; };
  void get_gridgen_info(sw4_type& order, float_sw4& zetaBreak) const;

  virtual void assignSw4_TypeerfaceSurfaces(EW* a_ew, Sarray& TopoGridExt) = 0;
  virtual void generate_z_and_j(EW* a_ew, sw4_type g, Sarray& z, Sarray& J) = 0;
  virtual void generate_grid_and_met(EW* a_ew, sw4_type g, Sarray& a_x, Sarray& a_y,
                                     Sarray& a_z, Sarray& a_jac, Sarray& a_met,
                                     bool a_comm = true) = 0;
  virtual bool grid_mapping(EW* a_ew, float_sw4 p, float_sw4 q, float_sw4 r,
                            sw4_type g, float_sw4& x, float_sw4& y,
                            float_sw4& z) = 0;
  virtual bool inverse_grid_mapping(EW* a_ew, float_sw4 x, float_sw4 y,
                                    float_sw4 z, sw4_type g, float_sw4& p,
                                    float_sw4& q, float_sw4& r) = 0;
  virtual void grid_mapping_diff(EW* a_ew, float_sw4 q, float_sw4 r,
                                 float_sw4 s, sw4_type g, sw4_type ic, sw4_type jc, sw4_type kc,
                                 float_sw4& zq, float_sw4& zr, float_sw4& zs,
                                 float_sw4& zqq, float_sw4& zqr, float_sw4& zqs,
                                 float_sw4& zrr, float_sw4& zrs,
                                 float_sw4& zss);
  virtual sw4_type sw4_typeerpolate_topography(EW* a_ew, float_sw4 x, float_sw4 y,
                                     float_sw4& z, Sarray& topo);
  virtual bool exact_metric(EW* a_ew, sw4_type g, Sarray& a_jac, Sarray& a_met);
  virtual void fill_topo(Sarray& topo, float_sw4 h);
};
#endif
