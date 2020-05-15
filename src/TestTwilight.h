#ifndef SW4_TESTTWILIGHT_H
#define SW4_TESTTWILIGHT_H

#include "Sarray.h"
#include "sw4.h"

class TestTwilight {
  float_sw4 m_omega, m_c, m_phase, m_momega, m_mphase, m_amprho, m_ampmu,
      m_amplambda;
  bool m_sw4twilight;

 public:
  TestTwilight(float_sw4 omega, float_sw4 c, float_sw4 phase, float_sw4 momega,
               float_sw4 mphase, float_sw4 amprho, float_sw4 ampmu,
               float_sw4 amplambda);
  void get_rho(Sarray& rho, Sarray& x, Sarray& y, Sarray& z);
  void get_mula(Sarray& mu, Sarray& lambda, Sarray& x, Sarray& y, Sarray& z);
  void get_ubnd(Sarray& u, Sarray& x, Sarray& y, Sarray& z, float_sw4 t,
                int npts, int sides[6]);
  void get_ubnd(Sarray& u, float_sw4 h, float_sw4 zmin, float_sw4 t, int npts,
                int sides[6]);
  void get_mula_att(Sarray& muve, Sarray& lambdave, Sarray& x, Sarray& y,
                    Sarray& z);
  void get_bnd_att(Sarray& AlphaVE, Sarray& x, Sarray& y, Sarray& z,
                   float_sw4 t, int npts, int sides[6]);
};
#endif
