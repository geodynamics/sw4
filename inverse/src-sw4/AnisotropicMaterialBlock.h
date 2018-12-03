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
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
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
#ifndef ANISOTROPIC_MATERIALBLOCK_H
#define ANISOTROPIC_MATERIALBLOCK_H

#include "AnisotropicMaterial.h"

class EW;

class AnisotropicMaterialBlock : public AnisotropicMaterial
{
 public:
   AnisotropicMaterialBlock( EW * a_ew, float_sw4 rho, float_sw4 c[21], float_sw4 xmin, float_sw4 xmax, float_sw4 ymin,
			     float_sw4 ymax, float_sw4 zmin, float_sw4 zmax ); 
   virtual void set_material_properties( std::vector<Sarray> &rho, std::vector<Sarray> &c );

   void set_gradients( float_sw4 rhograd, float_sw4 cgrad[21] );
   void set_absoluteDepth( bool absDepth );

private:
  bool inside_block( float_sw4 x, float_sw4 y, float_sw4 z );
  float_sw4 m_rho, m_rhograd, m_c[21], m_cgrad[21];
  float_sw4 m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
  float_sw4 m_tol;
  bool m_absoluteDepth;
  EW *mEW;
};

#endif
