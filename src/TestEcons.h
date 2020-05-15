#ifndef TESTECONS_H
#define TESTECONS_H

#include "Sarray.h"

class TestEcons {
  float_sw4 m_amp, m_cpocs;

 public:
  TestEcons(float_sw4 amp, float_sw4 cpocs);
  void get_rhobnd(Sarray& rho, int npts, int sides[6]);
  void get_mulabnd(Sarray& mu, Sarray& labmda, int npts, int sides[6]);
  void get_ubnd(Sarray& u, int npts, int sides[6]);
};

#endif
