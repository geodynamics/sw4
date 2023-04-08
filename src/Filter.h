// -*-c++-*-
//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// # Produced at the Lawrence Livermore National Laboratory.
// #
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// #
// # LLNL-CODE-643337
// #
// # All rights reserved.
// #
// # This file is part of SW4, Version: 1.0
// #
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
// #
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991.
// #
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details.
// #
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#ifndef FILTER_H
#define FILTER_H

#include <vector>

#include "SecondOrderSection.h"
#include "sw4.h"

using namespace std;
// support for derived quantities of the time derivative are not yet implemented
enum FilterType { lowPass, bandPass };

class Filter {
  friend std::ostream &operator<<(std::ostream &output, const Filter &s);

 public:
  Filter(FilterType type, unsigned sw4_type numberOfPoles,
         unsigned sw4_type numberOfPasses, float_sw4 f1, float_sw4 f2);
  ~Filter();
  void evaluate(sw4_type N, float_sw4 *u, float_sw4 *mf);
  void computeSOS(float_sw4 dt);
  float_sw4 estimatePrecursor();

  FilterType get_type() { return m_type; }
  sw4_type get_order() { return m_poles; }
  sw4_type get_passes() { return m_passes; }
  float_sw4 get_corner_freq1() { return m_f1; }
  float_sw4 get_corner_freq2() { return m_f2; }

 private:
  Filter();
  float_sw4 realPoleBP(float_sw4 f1, float_sw4 f2, float_sw4 dt,
                       SecondOrderSection *&sos_ptr);
  float_sw4 complexConjugatedPolesBP(float_sw4 f1, float_sw4 f2, float_sw4 dt,
                                     float_sw4 alpha,
                                     SecondOrderSection *&sos1_ptr,
                                     SecondOrderSection *&sos2_ptr);
  float_sw4 realPoleLP(float_sw4 fc, float_sw4 dt,
                       SecondOrderSection *&sos_ptr);
  float_sw4 complexConjugatedPolesLP(float_sw4 fc, float_sw4 dt,
                                     float_sw4 alpha,
                                     SecondOrderSection *&sos_ptr);
  void a2d(float_sw4 n[3], float_sw4 d[3], Polynomial &b, Polynomial &a);

  FilterType m_type;
  float_sw4 m_dt, m_f1, m_f2;
  unsigned sw4_type m_passes, m_poles;
  vector<SecondOrderSection *> m_SOSp;
  sw4_type m_numberOfSOS;
  sw4_type m_real_poles, m_complex_pairs;
  bool m_initialized;
  float_sw4 m_pole_min_re;
};

#endif
