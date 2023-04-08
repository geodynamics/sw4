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
#ifndef MATERIALINVTEST_H
#define MATERIALINVTEST_H

#include "MaterialData.h"

class EW;

class MaterialInvtest : public MaterialData {
 public:
  MaterialInvtest(EW* a_ew, sw4_type nr);
  void set_material_properties(std::vector<Sarray>& rho,
                               std::vector<Sarray>& cs, std::vector<Sarray>& cp,
                               std::vector<Sarray>& xis,
                               std::vector<Sarray>& xip);

 private:
  sw4_type m_nr;
  EW* mEW;
  void invtestmtrl(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                   float_sw4* rho, float_sw4* cs, float_sw4* cp, float_sw4 h,
                   float_sw4 zmin, sw4_type nr);
  void invtestmtrlc(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                    float_sw4* rho, float_sw4* cs, float_sw4* cp, float_sw4* xx,
                    float_sw4* yy, float_sw4* zz, sw4_type nr);
};

#endif
