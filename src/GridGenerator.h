#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H
class Sarray;
class EW;

#include "sw4.h"

class GridGenerator {
 protected:
  float_sw4 m_topo_zmax, m_zetaBreak;
  bool m_always_new;
  int m_grid_interpolation_order;
  int metric_ci(int ib, int ie, int jb, int je, int kb, int ke, float_sw4* a_x,
                float_sw4* a_y, float_sw4* a_z, float_sw4* a_met,
                float_sw4* a_jac);
  void gettopowgh(float_sw4 ai, float_sw4 wgh[8]) const;
  float_sw4 curvilinear_interface_parameter(EW* a_ew, int gcurv);

 public:
  GridGenerator(float_sw4 topo_zmax, bool always_new,
                int grid_interpolation_order, float_sw4 zetaBreak)
      : m_topo_zmax(topo_zmax),
        m_zetaBreak(zetaBreak),
        m_always_new(always_new),
        m_grid_interpolation_order(grid_interpolation_order){};
  bool curviCartIsSmooth(int ncurv) { return ncurv == 1 && !m_always_new; }
  float_sw4 get_topo_zmax() { return m_topo_zmax; };
  void get_gridgen_info(int& order, float_sw4& zetaBreak) const;

  virtual void assignInterfaceSurfaces(EW* a_ew, Sarray& TopoGridExt) = 0;
  virtual void generate_z_and_j(EW* a_ew, int g, Sarray& z, Sarray& J) = 0;
  virtual void generate_grid_and_met(EW* a_ew, int g, Sarray& a_x, Sarray& a_y,
                                     Sarray& a_z, Sarray& a_jac, Sarray& a_met,
                                     bool a_comm = true) = 0;
  virtual bool grid_mapping(EW* a_ew, float_sw4 p, float_sw4 q, float_sw4 r,
                            int g, float_sw4& x, float_sw4& y,
                            float_sw4& z) = 0;
  virtual bool inverse_grid_mapping( EW* a_ew, float_sw4 x, float_sw4 y, float_sw4 z, int g,
                                      float_sw4& p, float_sw4& q, float_sw4& r, bool interior=true )=0;
  virtual void grid_mapping_diff(EW* a_ew, float_sw4 q, float_sw4 r,
                                 float_sw4 s, int g, int ic, int jc, int kc,
                                 float_sw4& zq, float_sw4& zr, float_sw4& zs,
                                 float_sw4& zqq, float_sw4& zqr, float_sw4& zqs,
                                 float_sw4& zrr, float_sw4& zrs,
                                 float_sw4& zss);
  virtual int interpolate_topography(EW* a_ew, float_sw4 x, float_sw4 y,
                                     float_sw4& z, Sarray& topo);
  virtual bool exact_metric(EW* a_ew, int g, Sarray& a_jac, Sarray& a_met);
  virtual void fill_topo(Sarray& topo, float_sw4 h);
};
#endif
