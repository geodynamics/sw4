#ifndef GRIDGENERATORGENERAL_H
#define GRIDGENERATORGENERAL_H

#include <vector>

#include "GridGenerator.h"

class GridGeneratorGeneral : public GridGenerator {
  std::vector<Sarray> m_curviInterface;

  void generate_grid_and_met_new(EW* a_ew, int g, Sarray& a_x, Sarray& a_y,
                                 Sarray& a_z, Sarray& a_jac, Sarray& a_met);
  void generate_grid_and_met_old(EW* a_ew, Sarray& a_x, Sarray& a_y,
                                 Sarray& a_z, Sarray& a_jac, Sarray& a_met);
  bool grid_mapping_old(float_sw4 p, float_sw4 q, float_sw4 r, int g,
                        float_sw4& x, float_sw4& y, float_sw4& z,
                        Sarray& TopoGridExt, float_sw4 h, int Nz);
  bool inverse_grid_mapping_old(EW* a_ew, float_sw4 p, float_sw4 q, float_sw4 r,
                                int g, float_sw4& x, float_sw4& y, float_sw4& z,
                                Sarray& TopoGridExt, float_sw4 h, int Nz);
  bool inverse_grid_mapping_new(EW* a_ew, float_sw4 x, float_sw4 y, float_sw4 z,
                                int g, float_sw4& q, float_sw4& r, float_sw4& s,
                                float_sw4 h, int Nz);

 public:
  GridGeneratorGeneral(float_sw4 topo_zmax, bool always_new,
                       int grid_interpolation_order, float_sw4 zetaBreak);
  void assignInterfaceSurfaces(EW* a_ew, Sarray& TopoGridExt);
  virtual void generate_grid_and_met(EW* a_ew, int g, Sarray& a_x, Sarray& a_y,
                                     Sarray& a_z, Sarray& a_jac, Sarray& a_met,
                                     bool a_comm = true);
  virtual bool grid_mapping(EW* a_ew, float_sw4 p, float_sw4 q, float_sw4 r,
                            int g, float_sw4& x, float_sw4& y, float_sw4& z);
  virtual bool inverse_grid_mapping(EW* a_ew, float_sw4 x, float_sw4 y,
                                    float_sw4 z, int g, float_sw4& p,
                                    float_sw4& q, float_sw4& r);
};
#endif
