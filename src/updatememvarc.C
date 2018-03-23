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

#include "sw4.h"
#include <sys/types.h>
#include "policies.h"
#include "caliper.h"

void memvar_pred_fort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			  float_sw4* __restrict__ alp, float_sw4* __restrict__ alm,
			  float_sw4* __restrict__ u, float_sw4 omega, float_sw4 dt, int domain )
{
SW4_MARK_FUNCTION;
//***********************************************************************
//*** 
//*** domain = 0 --> Entire domain
//*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
//*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
//***
// AP: this routine implementes the 2nd order predictor step for evolving the memory variables
//***********************************************************************      

   float_sw4 dto = dt*omega;
   float_sw4 icp = 1/( 1.0/2 + 1/(2*dto) );
   float_sw4 cm = 0.5 - 1/(2*dto);

   int k1,k2;
   if( domain==0 )
   {
      k1 = kfirst;
      k2 = klast;
   }
   else if( domain==1 )
   {
      k1 = kfirst;
      k2 = kfirst+2;
   }
   else if( domain == 2 )
   {
      k1 = klast-2;
      k2 = klast;
   }
   const size_t ni  = ilast-ifirst+1;
   const size_t nij = ni*(jlast-jfirst+1);
   const size_t nijk= nij*(klast-kfirst+1);
   const int base = -ifirst-ni*jfirst-nij*kfirst;
//    for( int c=0 ; c < 3 ;c++)
// #pragma omp parallel for
//       for( int k=k1 ; k <= k2 ; k++)
// 	 for( int j=jfirst ; j<= jlast; j++ )
// 	    for( int i=ifirst ; i<= ilast; i++ )
// 	    {
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment c_range(0,3);
   RAJA::kernel<DEFAULT_LOOP4>(
				  RAJA::make_tuple(i_range,j_range,k_range,c_range),
				  [=]RAJA_DEVICE (int i,int j, int k,int c) {
	       size_t ind = base+i+ni*j+nij*k;
	       alp[ind+c*nijk] = icp*(-cm*alm[ind+c*nijk] + u[ind+c*nijk] );
				  });
}

//-----------------------------------------------------------------------
void memvar_corr_fort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			  float_sw4* __restrict__ alp, float_sw4* __restrict__ alm,
			  float_sw4* __restrict__ up, float_sw4* __restrict__ u, 
			  float_sw4* __restrict__ um, float_sw4 omega, float_sw4 dt, int domain )
{
  SW4_MARK_FUNCTION;
  //***********************************************************************
  //*** 
  //*** domain = 0 --> Entire domain
  //*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  //*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  //
  // AP Apr. 3, 2017: corrector step for updating memory variables
  // AP June 14, 2017: make corrector step independent of predictor step to simplify
  // the mesh refinement algorithm
  //***********************************************************************      

   const float_sw4 i6=1.0/6;
   float_sw4 dto = dt*omega;

   float_sw4 icp = 1/( 0.5 + 1/(2*dto) + dto/4 + dto*dto/12 );
   float_sw4 cm = 1/(2*dto) + dto/4  - 0.5 - dto*dto/12;

   int k1,k2;
   if( domain==0 )
   {
      k1 = kfirst;
      k2 = klast;
   }
   else if( domain==1 )
   {
      k1 = kfirst;
      k2 = kfirst+2;
   }
   else if( domain == 2 )
   {
      k1 = klast-2;
      k2 = klast;
   }
   const size_t ni  = ilast-ifirst+1;
   const size_t nij = ni*(jlast-jfirst+1);
   const size_t nijk= nij*(klast-kfirst+1);
   const int base = -ifirst-ni*jfirst-nij*kfirst;
//    for( int c=0 ; c < 3 ;c++)
// #pragma omp parallel for
//       for( int k=k1 ; k <= k2 ; k++)
// 	 for( int j=jfirst ; j<= jlast; j++ )
// 	    for( int i=ifirst ; i<= ilast; i++ )
// 	    {
RAJA::RangeSegment i_range(ifirst,ilast+1);
RAJA::RangeSegment j_range(jfirst,jlast+1);
RAJA::RangeSegment k_range(k1,k2+1);
RAJA::RangeSegment c_range(0,3);
RAJA::kernel<DEFAULT_LOOP4>(RAJA::make_tuple(i_range,j_range,k_range,c_range),
			    [=]RAJA_DEVICE (int i,int j, int k,int c) {
			    size_t ind = base+i+ni*j+nij*k;
			    // Note that alp is ASSIGNED by this formula
			    alp[ind+c*nijk] = icp*( cm*alm[ind+c*nijk] + u[ind+c*nijk] + i6* ( dto*dto*u[ind+c*nijk] + 
											       dto*(up[ind+c*nijk]-um[ind+c*nijk]) + (up[ind+c*nijk]-2*u[ind+c*nijk]+um[ind+c*nijk]) ) );
			    });
}

//-----------------------------------------------------------------------
void memvar_corr_fort_wind_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			       float_sw4* __restrict__  alp,
			       int d1b, int d1e, int d2b, int d2e, int d3b, int d3e,
                               float_sw4* __restrict__ alm, float_sw4* __restrict__ up,
			       float_sw4* __restrict__ u, float_sw4* __restrict__ um,
			       float_sw4 omega, float_sw4 dt, int domain )
{
  SW4_MARK_FUNCTION;
  //***********************************************************************
  //*** 
  //*** domain = 0 --> Entire domain
  //*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  //*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  //***
  // AP Apr. 3, 2017: corrector step for updating memory variables
  // AP June 14, 2017: make corrector step independent of predictor step to simplify
  // the mesh refinement algorithm
  //***********************************************************************      

   const float_sw4 i6=1.0/6;
   float_sw4 dto = dt*omega;
   float_sw4 icp = 1/( 0.5 + 1/(2*dto) + dto/4 + dto*dto/12 );
   float_sw4 cm = 1/(2*dto) + dto/4  - 0.5 - dto*dto/12;

   int k1,k2;
   if( domain==0 )
   {
      k1 = kfirst;
      k2 = klast;
   }
   else if( domain==1 )
   {
      k1 = kfirst;
      k2 = kfirst+2;
   }
   else if( domain == 2 )
   {
      k1 = klast-2;
      k2 = klast;
   }
   const size_t ni  = ilast-ifirst+1;
   const size_t nij = ni*(jlast-jfirst+1);
   const size_t nijk= nij*(klast-kfirst+1);
   const int base = -ifirst-ni*jfirst-nij*kfirst;

   const size_t dni  = d1e-d1b+1;
   const size_t dnij = dni*(d2e-d2b+1);
   const size_t dnijk= dnij*(d3e-d3b+1);
   const int dbase = -d1b-dni*d2b-dnij*d3b;


   //  real*8 up(3,d1b:d1e, d2b:d2e, d3b:d3e)
   //  real*8  u(3,d1b:d1e, d2b:d2e, d3b:d3e);
   //  real*8 um(3,d1b:d1e, d2b:d2e, d3b:d3e);
   //  real*8 alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast) // different sizes here
   //  real*8 alm(3,d1b:d1e, d2b:d2e, d3b:d3e)

   for( int c=0 ; c < 3 ;c++)
#pragma omp parallel for
      for( int k=k1 ; k <= k2 ; k++)
	 for( int j=jfirst ; j<= jlast; j++ )
	    for( int i=ifirst ; i<= ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
	       size_t dind = dbase+i+dni*j+dnij*k;
	       alp[ind+c*nijk] = icp*( cm*alm[dind+c*dnijk] + u[dind+c*dnijk] + i6* ( dto*dto*u[dind+c*dnijk] + 
	      dto*(up[dind+c*dnijk]-um[dind+c*dnijk]) + (up[dind+c*dnijk]-2*u[dind+c*dnijk]+um[dind+c*dnijk]) ) );
	    }
}

