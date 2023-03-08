#ifndef MPARGRIDFILE_H
#define MPARGRIDFILE_H

#include <string>

#include "sw4.h"

// Class representing a file holding a material parameter grid.
// Variable `vars' has the following meaning:
//  vars =1 (rho,mu,lambda) are stored on file
//        2 (rho,cs,cp) stored on file
//        3 (cp) stored on file
//        4 (cs,cp) stored on file
//  The parameters are in xms, addressed as
//  xms[c+nc*(i-1)+nc*ni*(j-1)+nc*ni*nj*(k-1)]; where nc=1,2, or 3, depending on
//  value of `vars'.

class MParGridFile {
 public:
  MParGridFile(std::string fname);
  MParGridFile(float_sw4* xms, int vars, int nx, int ny, int nz, float_sw4 hx,
               float_sw4 hy, float_sw4 hz, float_sw4 xmin, float_sw4 ymin,
               float_sw4 zmin, int update = 1);
  ~MParGridFile();
  void interpolate_to_other(float_sw4* xms, int vars, int nx, int ny, int nz,
                            float_sw4 hx, float_sw4 hy, float_sw4 hz,
                            float_sw4 xmin, float_sw4 ymin, float_sw4 zmin);
  void write_mparcartfile(std::string fname);
  bool is_update() { return m_update == 1; }

 private:
  void read_mparcartfile(std::string fname);
  int m_vars, m_nxr, m_nyr, m_nzr, m_update;
  float_sw4 m_xrmin, m_yrmin, m_zrmin, m_hxr, m_hyr, m_hzr;
  float_sw4* m_xms;
};
#endif
