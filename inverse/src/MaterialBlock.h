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
#ifndef MATERIALBLOCK_H
#define MATERIALBLOCK_H

#include "MaterialData.h"

class EW;
class MaterialBlock : public MaterialData
{
 public:
   MaterialBlock( EW * a_ew,double rho, double vs, double vp, double xmin, double xmax, double ymin,
		  double ymax, double zmin, double zmax, double qs=-1, double qp=-1,
		  double freq=1 );

   virtual void set_material_properties( std::vector<Sarray> &rho, std::vector<Sarray> &cs,
					 std::vector<Sarray> &cp,
					 std::vector<Sarray>& xis, std::vector<Sarray>& xip);

   void set_gradients( double rhograd, double vsgrad, double vpgrad );
   void set_absoluteDepth( bool absDepth );

private:
  bool inside_block( double x, double y, double z );
  double m_rho, m_vp, m_vs, m_qp, m_qs, m_freq;
  double m_vpgrad, m_vsgrad, m_rhograd;
  double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
  double m_tol;
  bool m_absoluteDepth;
  EW *mEW; // where is this pointer needed?
};

#endif
