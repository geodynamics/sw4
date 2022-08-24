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
#include "MaterialInvtest.h"


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
   float_sw4* rho_ptr, *cs_ptr, *cp_ptr;

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
      float_sw4 h = mEW->mGridSize[g]; 
      float_sw4 zmin = mEW->m_zmin[g];
      invtestmtrl( ifirst, ilast, jfirst, jlast, kfirst, klast,
		   rho_ptr, cs_ptr, cp_ptr, h, zmin, m_nr );
   }
   for( int g=mEW->mNumberOfCartesianGrids; g < mEW->mNumberOfGrids; g++)
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
      float_sw4* x_ptr = mEW->mX[g].c_ptr();
      float_sw4* y_ptr = mEW->mY[g].c_ptr();
      float_sw4* z_ptr = mEW->mZ[g].c_ptr();
      invtestmtrlc( ifirst, ilast, jfirst, jlast, kfirst, klast,
		    rho_ptr, cs_ptr, cp_ptr, x_ptr, y_ptr, z_ptr, m_nr );
   }
}

//-----------------------------------------------------------------------
void MaterialInvtest::invtestmtrl( int ib, int ie, int jb, int je, int kb, int ke, 
				   float_sw4* rho, float_sw4* cs, float_sw4* cp,
				   float_sw4 h, float_sw4 zmin, int nr )
{
   const float_sw4 ep = 0.01;
   const float_sw4 pi2 = atan(1.0)*8;  // 2*pi
   const float_sw4 om = pi2;
   const float_sw4 xminbox = 0.4;
   const float_sw4 xmaxbox = 0.6;
   const float_sw4 yminbox = 0.4;
   const float_sw4 ymaxbox = 0.6;
   const float_sw4 zminbox = 0.1;
   const float_sw4 zmaxbox = 0.3;
   float_sw4 omx, omy, omz, phip, phir, phis;
   float_sw4 cp1, cp2, cp3, cp4, cs1, cs2, cs3, cs4;
   float_sw4 rho1, rho2, rho3, rho4, z1, z2, z3, Lz=200;
   if( nr == 3 )
   {
      omx = pi2*2.0/30000;
      omy = pi2*2.0/30000;
      omz = pi2*1.0/17000;
      phip = 0.3;
      phir = 0.17;
      phis = 0.08;
   }
   else if( nr == 4 )
   {
      cp1=4000;
      cp2=6000;
      cp3=4000;
      cp4=6000;

      cs1=2000;
      cs2=3464;
      cs3=2000;
      cs4=3464;

      rho1 = 2600;
      rho2 = 2700;
      rho3 = 2600;
      rho4 = 2700;

      z1=1000;
      z2=3000;
      z3=5000;
      Lz=200;
   }
   size_t ind=0;
   for( int k=kb ; k<= ke ; k++ )
   {
      float_sw4 z = zmin + (k-1)*h;
      for( int j=jb ; j<= je ; j++ )
      {
	 float_sw4 y = (j-1)*h;
	 for( int i=ib ; i<= ie ; i++ )
	 {
	    float_sw4 x = (i-1)*h;
	    if( nr == 1 )
	    {
	       rho[ind] = 1+ep*sin(om*x+0.13)*sin(om*y)*sin(om*z);
	       cs[ind]  = 2+ep*cos(om*x)*sin(om*y)*cos(om*z+0.01);
	       cp[ind]  = 4+ep*sin(om*x+0.4)*sin(om*y)*cos(om*z+0.1);
	    }
	    else if( nr == 2 )
	    {
	       rho[ind] = 1;
	       cs[ind]  = 2;
	       cp[ind]  = 4;
	       if( x >= xminbox && x <= xmaxbox &&
		   y >= yminbox && y <= ymaxbox && 
		   z >= zminbox && z <= zmaxbox   )
	       {
		  rho[ind] = 1.3;
		  cs[ind]  = 0.5;
		  cp[ind]  = 2;
	       }
	    }
	    else if( nr == 3 )
	    {
                  rho[ind]=2650+
		     50*sin(omx*x)*sin(omy*y)*cos(omz*z+phir);
                  cp[ind]=5000+
                     1000*sin(omx*x+phip)*sin(omy*y)*cos(omz*z);
                  cs[ind] = cp[ind]/( 1.8660 + 
		     0.1339*cos(omx*x)*sin(omy*y+phis)*cos(omz*z) );
	    }
	    else if( nr == 4 )
	    {
	       rho[ind] = rho1 
		  + 0.5*(rho2-rho1)*(1.0 + tanh((z-z1)/Lz))
		  + 0.5*(rho3-rho2)*(1.0 + tanh((z-z2)/Lz))
		  + 0.5*(rho4-rho3)*(1.0 + tanh((z-z3)/Lz));
	       cp[ind] = cp1 
		  + 0.5*(cp2-cp1)*(1.0 + tanh((z-z1)/Lz))
		  + 0.5*(cp3-cp2)*(1.0 + tanh((z-z2)/Lz))
		  + 0.5*(cp4-cp3)*(1.0 + tanh((z-z3)/Lz));
	       cs[ind] = cs1
		  + 0.5*(cs2-cs1)*(1.0 + tanh((z-z1)/Lz))
		  + 0.5*(cs3-cs2)*(1.0 + tanh((z-z2)/Lz))
		  + 0.5*(cs4-cs3)*(1.0 + tanh((z-z3)/Lz));
	    }
	    ind++;
	 }
      }
   }
}

//-----------------------------------------------------------------------
void MaterialInvtest::invtestmtrlc( int ib, int ie, int jb, int je, int kb, int ke, 
				    float_sw4* rho, float_sw4* cs, float_sw4* cp,
				    float_sw4* xx, float_sw4* yy, float_sw4* zz, int nr )
{
   const float_sw4 ep = 0.01;
   const float_sw4 pi2 = atan(1.0)*8;  // 2*pi
   const float_sw4 om = pi2;
   const float_sw4 xminbox = 0.4;
   const float_sw4 xmaxbox = 0.6;
   const float_sw4 yminbox = 0.4;
   const float_sw4 ymaxbox = 0.6;
   const float_sw4 zminbox = 0.1;
   const float_sw4 zmaxbox = 0.3;
   float_sw4 omx, omy, omz, phip, phir, phis;
   float_sw4 cp1, cp2, cp3, cp4, cs1, cs2, cs3, cs4;
   float_sw4 rho1, rho2, rho3, rho4, z1, z2, z3, Lz=200;
   if( nr == 3 )
   {
      omx = pi2*2.0/30000;
      omy = pi2*2.0/30000;
      omz = pi2*1.0/17000;
      phip = 0.3;
      phir = 0.17;
      phis = 0.08;
   }
   else if( nr == 4 )
   {
      cp1=4000;
      cp2=6000;
      cp3=4000;
      cp4=6000;

      cs1=2000;
      cs2=3464;
      cs3=2000;
      cs4=3464;

      rho1 = 2600;
      rho2 = 2700;
      rho3 = 2600;
      rho4 = 2700;

      z1=1000;
      z2=3000;
      z3=5000;
      Lz=200;
   }
   size_t ind=0;
   for( int k=kb ; k<= ke ; k++ )
      for( int j=jb ; j<= je ; j++ )
	 for( int i=ib ; i<= ie ; i++ )
	 {
	    float_sw4 x = xx[ind];
            float_sw4 y = yy[ind];
            float_sw4 z = zz[ind];
	    if( nr == 1 )
	    {
	       rho[ind] = 1+ep*sin(om*x+0.13)*sin(om*y)*sin(om*z);
	       cs[ind]  = 2+ep*cos(om*x)*sin(om*y)*cos(om*z+0.01);
	       cp[ind]  = 4+ep*sin(om*x+0.4)*sin(om*y)*cos(om*z+0.1);
	    }
	    else if( nr == 2 )
	    {
	       rho[ind] = 1;
	       cs[ind]  = 2;
	       cp[ind]  = 4;
	       if( x >= xminbox && x <= xmaxbox &&
		   y >= yminbox && y <= ymaxbox && 
		   z >= zminbox && z <= zmaxbox   )
	       {
		  rho[ind] = 1.3;
		  cs[ind]  = 0.5;
		  cp[ind]  = 2;
	       }
	    }
	    else if( nr == 3 )
	    {
	       rho[ind]=2650 +
		     50*sin(omx*x)*sin(omy*y)*cos(omz*z+phir);
	       cp[ind]=5000 +
                    1000*sin(omx*x+phip)*sin(omy*y)*cos(omz*z);
               cs[ind] = cp[ind]/( 1.8660 + 
		     0.1339*cos(omx*x)*sin(omy*y+phis)*cos(omz*z) );
                  //	       cs[ind] =2732 +
                  //		     732*cos(omx*x)*sin(omy*y+phis)*cos(omz*z);
	    }
	    else if( nr == 4 )
	    {
	       rho[ind] = rho1 
                     + 0.5*(rho2-rho1)*(1.0 + tanh((z-z1)/Lz))
                     + 0.5*(rho3-rho2)*(1.0 + tanh((z-z2)/Lz))
		     + 0.5*(rho4-rho3)*(1.0 + tanh((z-z3)/Lz));
	       cp[ind] = cp1 
                     + 0.5*(cp2-cp1)*(1.0 + tanh((z-z1)/Lz))
                     + 0.5*(cp3-cp2)*(1.0 + tanh((z-z2)/Lz))
		     + 0.5*(cp4-cp3)*(1.0 + tanh((z-z3)/Lz));
	       cs[ind] = cs1
                     + 0.5*(cs2-cs1)*(1.0 + tanh((z-z1)/Lz))
                     + 0.5*(cs3-cs2)*(1.0 + tanh((z-z2)/Lz))
	 	     + 0.5*(cs4-cs3)*(1.0 + tanh((z-z3)/Lz));
	    }
	    ind++;
	 }
}
