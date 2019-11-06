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
#include <math.h>

#include "TestRayleighWave.h"

extern "C" {
double rvel(double *lambda, double *mu);
}

TestRayleighWave::TestRayleighWave(float_sw4 rho, float_sw4 cs, float_sw4 cp,
                                   int nwl, float_sw4 xmax)
    : m_rho(rho),
      m_cs(cs),
      m_cp(cp),
      m_alpha(0.0)
// m_omega must be assigned before the exact solution can be evaluated
{
  // calculate wave length and wave number
  float_sw4 Lwave = xmax / (cos(m_alpha) * nwl);
  m_omega = 2 * M_PI / Lwave;

  // calculate the phase velocity
  m_mu = static_cast<double>(m_cs * m_cs * m_rho);
  m_lambda = static_cast<double>(m_cp * m_cp * m_rho - 2 * m_mu);
  double xi = rvel(&m_lambda, &m_mu);
  m_cr = xi * m_cs;
}
