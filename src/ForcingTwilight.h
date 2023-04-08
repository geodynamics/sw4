//-*-c++-*-
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
#ifndef FORCING_TWILIGHT_H
#define FORCING_TWILIGHT_H

#include <string>

#include "Sarray.h"
#include "boundaryConditionTypes.h"
#include "sw4.h"

class ForcingTwilight {
 public:
  ForcingTwilight(float_sw4 omega, float_sw4 c, float_sw4 phase,
                  float_sw4 momega, float_sw4 mphase, float_sw4 amprho,
                  float_sw4 ampmu, float_sw4 amplambda);

  // void default_bcs( boundaryConditionType bcs[6] );

  float_sw4 m_omega, m_c, m_phase, m_momega, m_mphase, m_amprho, m_ampmu,
      m_amplambda;
  bool m_use_attenuation;
  sw4_type m_number_mechanisms;

 private:
  ForcingTwilight(const ForcingTwilight&);
  ForcingTwilight& operator=(const ForcingTwilight&);
};

#endif
