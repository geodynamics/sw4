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
#include "EW.h"
#include "F77_FUNC.h"
#include "MaterialInvtest.h"

extern "C" {
   void F77_FUNC(invtestmtrl,INVTESTMTRL)(  int*, int*, int*, int*, int*, int*, double*, double*,
				       double*, double*, double*, int* );
   void F77_FUNC(invtestmtrlc,INVTESTMTRLC)( int*, int*, int*, int*, int*, int*, double*, double*,
					double*, double*, double*, double*, int* );
}

//-----------------------------------------------------------------------
MaterialInvtest::MaterialInvtest( EW* a_ew, int nr ) :
   mEW(a_ew),
   m_nr(nr)
{
   mCoversAllPoints = true;
}

//-----------------------------------------------------------------------
void MaterialInvtest::set_material_properties( std::vector<Sarray> & rho, std::vector<Sarray> & cs,
					       std::vector<Sarray> & cp, std::vector<Sarray>& xis,
					       std::vector<Sarray>& xip)
{
   int ifirst, ilast, jfirst, jlast, kfirst, klast;
   double* rho_ptr, *cs_ptr, *cp_ptr;

   for( int g = 0 ; g < mEW->mNumberOfCartesianGrids; g++) 
   {
      rho_ptr   = rho[g].c_ptr();
      cs_ptr    = cs[g].c_ptr();
      cp_ptr    = cp[g].c_ptr();
      ifirst = mEW->m_iStart[g];
      ilast  = mEW->m_iEnd[g];
      jfirst = mEW->m_jStart[g];
      jlast  = mEW->m_jEnd[g];
      kfirst = mEW->m_kStart[g];
      klast  = mEW->m_kEnd[g];
      double h = mEW->mGridSize[g]; 
      double zmin = mEW->m_zmin[g];
      F77_FUNC(invtestmtrl,INVTESTMTRL)( &ifirst, &ilast, &jfirst, &jlast,
					 &kfirst, &klast, rho_ptr, cs_ptr,
					 cp_ptr, &h, &zmin, &m_nr);
   }
   if( mEW->topographyExists() )
   {
      int g = mEW->mNumberOfGrids-1;
      rho_ptr   = rho[g].c_ptr();
      cs_ptr    = cs[g].c_ptr();
      cp_ptr    = cp[g].c_ptr();
      ifirst = mEW->m_iStart[g];
      ilast  = mEW->m_iEnd[g];
      jfirst = mEW->m_jStart[g];
      jlast  = mEW->m_jEnd[g];
      kfirst = mEW->m_kStart[g];
      klast  = mEW->m_kEnd[g];
      double* x_ptr = mEW->mX.c_ptr();
      double* y_ptr = mEW->mY.c_ptr();
      double* z_ptr = mEW->mZ.c_ptr();
      F77_FUNC(invtestmtrlc,INVTESTMTRLC)( &ifirst, &ilast, &jfirst, &jlast,
					   &kfirst, &klast, rho_ptr, cs_ptr,
					   cp_ptr, x_ptr, y_ptr, z_ptr, &m_nr);
   }
}
