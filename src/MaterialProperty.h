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
#ifndef WPP_MATERIALPROPERTY_H
#define WPP_MATERIALPROPERTY_H

#include <string>
#include "sw4.h"

class MaterialProperty {
 public:
  // -----------------------------------------------------------------
  // class MaterialProperty
  //
  // Defines Vp, Vs, Rho and, optionally, Qp and Qs.
  // ------------------------------------------------------------------
  MaterialProperty(int id, float_sw4 vp0, float_sw4 vp1, float_sw4 vp2,
                   float_sw4 vs0, float_sw4 vs1, float_sw4 vs2, float_sw4 rho0,
                   float_sw4 rho1, float_sw4 rho2, float_sw4 qp, float_sw4 qs);

  void setSqrtCoefficients(float_sw4 vp1o2, float_sw4 vs1o2, float_sw4 rho1o2);

  int m_materialID;
  float_sw4 m_vp0, m_vp1, m_vp2, m_vp1o2, m_vs0, m_vs1, m_vs2, m_vs1o2, m_rho0,
      m_rho1, m_rho2, m_rho1o2, m_qp, m_qs;

 private:
  MaterialProperty();
};
#endif
