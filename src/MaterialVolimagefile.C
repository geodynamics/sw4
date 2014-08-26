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
#include "MaterialVolimagefile.h"

//-----------------------------------------------------------------------
MaterialVolimagefile::MaterialVolimagefile( EW* a_ew, bool rhomula, std::string path, std::string rhofile, 
	std::string mufile, std::string lambdafile, std::string qpfile, std::string qsfile ) :
   mEW(a_ew)
{
   mCoversAllPoints = true;
   m_rhomula = rhomula;
   m_path = path;
   m_rho = rhofile;
   m_mu = mufile;
   m_lambda = lambdafile;
   m_qp = qpfile;
   m_qs = qsfile;
}

//-----------------------------------------------------------------------
void MaterialVolimagefile::set_material_properties( std::vector<Sarray> & rho, std::vector<Sarray> & cs,
				std::vector<Sarray> & cp, std::vector<Sarray>& xis, std::vector<Sarray>& xip)
{
   mEW->read_volimage( m_path, m_rho, rho );
   mEW->read_volimage( m_path, m_mu, cs );
   mEW->read_volimage( m_path, m_lambda, cp );
   //   mEW->check_min_max_int( rho );

   mEW->communicate_arrays( rho );
   mEW->communicate_arrays( cs );
   mEW->communicate_arrays( cp );
   mEW->update_curvilinear_cartesian_interface( rho );
   mEW->update_curvilinear_cartesian_interface( cs );
   mEW->update_curvilinear_cartesian_interface( cp );

   if( xis[0].is_defined() )
   {
      mEW->read_volimage( m_path, m_qs, xis );
      mEW->communicate_arrays( xis );
      mEW->update_curvilinear_cartesian_interface( xis );
// tmp
      // if (mEW->proc_zero())
      // 	cout << "qs[0] array is defined in MatVolimgfile()" << endl;
   }
//    else
//    {
// // tmp
//       if (mEW->proc_zero())
// 	cout << "qs[0] array is NOT defined in MatVolimgfile()" << endl;
//    }
   
   if( xip[0].is_defined() )
   {
      mEW->read_volimage( m_path, m_qp, xip );
      mEW->communicate_arrays( xip );
      mEW->update_curvilinear_cartesian_interface( xip );
   }


   if( m_rhomula )
   {
      for( int g = 0 ; g < mEW->mNumberOfGrids; g++) 
      {
	 // Convert to cp,cs,rho, here cs is mu, cp is lambda
	 double* rhoptr = rho[g].c_ptr();
	 double* csptr  = cs[g].c_ptr();
	 double* cpptr  = cp[g].c_ptr();
	 for( size_t i = 0 ; i < rho[g].npts() ; i++ )
	 {
            if( rhoptr[i] != -1 )
	    {
	       csptr[i] = csptr[i]/rhoptr[i];
	       cpptr[i] = 2*csptr[i] + cpptr[i]/rhoptr[i];
	       csptr[i] = sqrt(csptr[i]);
	       cpptr[i] = sqrt(cpptr[i]);
	    }
	 }
      }
   }
}
