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
//-----------------------------------------------------------------------
//  Adds 4th order artificial disssipation for super-grid damping layers
//
//-----------------------------------------------------------------------

#include "EW.h"
#include "policies.h"
#include "caliper.h"
//-----------------------------------------------------------------------
void EW::addsgd4_ci( int ifirst, int ilast, int jfirst, int jlast,
		     int kfirst, int klast,
		     float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		     float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		     float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
		     float_sw4* __restrict__ a_dcz, float_sw4* __restrict__ a_strx, 
		     float_sw4* __restrict__ a_stry, float_sw4* __restrict__ a_strz,
		     float_sw4* __restrict__ a_cox,  float_sw4* __restrict__ a_coy,
		     float_sw4* __restrict__ a_coz,
		     float_sw4 beta )
{
  SW4_MARK_FUNCTION;
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
#define strz(k) a_strz[(k-kfirst)]
#define dcz(k) a_dcz[(k-kfirst)]
#define coz(k) a_coz[(k-kfirst)]

      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);

// AP: The for c loop could be inside the for i loop. The simd, ivdep pragmas should be outside the inner-most loop
// #pragma omp parallel
//       {
//       for( int c=0 ; c < 3 ; c++ )
// #pragma omp for
//       for( int k=kfirst+2; k <= klast-2 ; k++ )
// 	 for( int j=jfirst+2; j <= jlast-2 ; j++ )
// #pragma simd
// #pragma ivdep
// 	    for( int i=ifirst+2; i <= ilast-2 ; i++ )
// 	    {
      ASSERT_MANAGED(a_dcx);
      ASSERT_MANAGED(a_cox);
      ASSERT_MANAGED(a_strx);
using ADDSGD_POL  =
       RAJA::KernelPolicy<
       RAJA::statement::CudaKernel<
	 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>,
			      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<4>,
						   RAJA::statement::For<3, RAJA::cuda_threadblock_exec<32>,
									RAJA::statement::For<0, RAJA::seq_exec,
											     RAJA::statement::Lambda<0> >>>>>>;
RAJA::RangeSegment i_range(ifirst+2,ilast-1);
RAJA::RangeSegment j_range(jfirst+2,jlast-1);
RAJA::RangeSegment k_range(kfirst+2,klast-1);
RAJA::RangeSegment c_range(0,3);
RAJA::kernel<ADDSGD_POL>(
			    RAJA::make_tuple(c_range,k_range,j_range,i_range),
			    [=]RAJA_DEVICE (int c,int k, int j,int i) {
			      float_sw4 birho=beta/rho(i,j,k);
			      {
		  up(c,i,j,k) -= birho*( 
		  // x-differences
		   strx(i)*coy(j)*coz(k)*(
       rho(i+1,j,k)*dcx(i+1)*
                   ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
      -2*rho(i,j,k)*dcx(i)  *
                   ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
      +rho(i-1,j,k)*dcx(i-1)*
                   ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
      -rho(i+1,j,k)*dcx(i+1)*
                   (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
      +2*rho(i,j,k)*dcx(i)  *
                   (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
      -rho(i-1,j,k)*dcx(i-1)*
                   (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
// y-differences
      stry(j)*cox(i)*coz(k)*(
      +rho(i,j+1,k)*dcy(j+1)*
                   ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
      -2*rho(i,j,k)*dcy(j)  *
                   ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
      +rho(i,j-1,k)*dcy(j-1)*
                   ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
      -rho(i,j+1,k)*dcy(j+1)*
                   (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
      +2*rho(i,j,k)*dcy(j)  *
                   (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
      -rho(i,j-1,k)*dcy(j-1)*
                   (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) +
       strz(k)*cox(i)*coy(j)*(
// z-differences
      +rho(i,j,k+1)*dcz(k+1)* 
                 ( u(c,i,j,k+2) -2*u(c,i,j,k+1)+ u(c,i,j,k  )) 
      -2*rho(i,j,k)*dcz(k)  *
                 ( u(c,i,j,k+1) -2*u(c,i,j,k  )+ u(c,i,j,k-1))
      +rho(i,j,k-1)*dcz(k-1)*
                 ( u(c,i,j,k  ) -2*u(c,i,j,k-1)+ u(c,i,j,k-2)) 
      -rho(i,j,k+1)*dcz(k+1)*
                 (um(c,i,j,k+2)-2*um(c,i,j,k+1)+um(c,i,j,k  )) 
      +2*rho(i,j,k)*dcz(k)  *
                 (um(c,i,j,k+1)-2*um(c,i,j,k  )+um(c,i,j,k-1)) 
      -rho(i,j,k-1)*dcz(k-1)*
                 (um(c,i,j,k  )-2*um(c,i,j,k-1)+um(c,i,j,k-2)) ) 
					 );

	       }
			    }); SYNC_STREAM;
//   }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef strz
#undef dcz
#undef coz
   }
}

//-----------------------------------------------------------------------
void EW::addsgd6_ci( int ifirst, int ilast, int jfirst, int jlast,
		     int kfirst, int klast,
		     float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		     float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		     float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
		     float_sw4* __restrict__ a_dcz, float_sw4* __restrict__ a_strx,
		     float_sw4* __restrict__ a_stry, float_sw4* __restrict__ a_strz,
		     float_sw4* __restrict__ a_cox,  float_sw4* __restrict__ a_coy,
		     float_sw4* __restrict__ a_coz, float_sw4 beta )
{
  SW4_MARK_FUNCTION;
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
#define strz(k) a_strz[(k-kfirst)]
#define dcz(k) a_dcz[(k-kfirst)]
#define coz(k) a_coz[(k-kfirst)]
      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
#pragma omp parallel
      {
      for( int c=0 ; c < 3 ; c++ )
#pragma omp for
      for( int k=kfirst+3; k <= klast-3 ; k++ )
	 for( int j=jfirst+3; j <= jlast-3 ; j++ )
#pragma simd
#pragma ivdep
	    for( int i=ifirst+3; i <= ilast-3 ; i++ )
	    {
	       float_sw4 birho=0.5*beta/rho(i,j,k);
	       {
		 up(c,i,j,k) += birho*( 
       strx(i)*coy(j)*coz(k)*(
// x-differences
         (rho(i+2,j,k)*dcx(i+2)+rho(i+1,j,k)*dcx(i+1))*(
         u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
      -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
      -3*(rho(i+1,j,k)*dcx(i+1)+rho(i,j,k)*dcx(i))*(
         u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
      -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
      +3*(rho(i,j,k)*dcx(i)+rho(i-1,j,k)*dcx(i-1))*(
         u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
      -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
       - (rho(i-1,j,k)*dcx(i-1)+rho(i-2,j,k)*dcx(i-2))*(
         u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
      -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
                 ) +  stry(j)*cox(i)*coz(k)*(
// y-differences
         (rho(i,j+2,k)*dcy(j+2)+rho(i,j+1,k)*dcy(j+1))*(
         u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
      -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
      -3*(rho(i,j+1,k)*dcy(j+1)+rho(i,j,k)*dcy(j))*(
         u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
      -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
      +3*(rho(i,j,k)*dcy(j)+rho(i,j-1,k)*dcy(j-1))*(
         u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
      -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
       - (rho(i,j-1,k)*dcy(j-1)+rho(i,j-2,k)*dcy(j-2))*(
         u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
      -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
                 ) +  strz(k)*cox(i)*coy(j)*(
// z-differences
         ( rho(i,j,k+2)*dcz(k+2) + rho(i,j,k+1)*dcz(k+1) )*(
         u(c,i,j,k+3)- 3*u(c,i,j,k+2)+ 3*u(c,i,j,k+1)- u(c,i,  j,k) 
      -(um(c,i,j,k+3)-3*um(c,i,j,k+2)+3*um(c,i,j,k+1)-um(c,i,  j,k)) )
      -3*(rho(i,j,k+1)*dcz(k+1)+rho(i,j,k)*dcz(k))*(
         u(c,i,j,k+2) -3*u(c,i,j,k+1)+ 3*u(c,i,  j,k)- u(c,i,j,k-1) 
      -(um(c,i,j,k+2)-3*um(c,i,j,k+1)+3*um(c,i,  j,k)-um(c,i,j,k-1)) )
      +3*(rho(i,j,k)*dcz(k)+rho(i,j,k-1)*dcz(k-1))*(
         u(c,i,j,k+1)- 3*u(c,i,  j,k)+ 3*u(c,i,j,k-1)-u(c,i,j,k-2) 
      -(um(c,i,j,k+1)-3*um(c,i,  j,k)+3*um(c,i,j,k-1)-um(c,i,j,k-2)) )
       - (rho(i,j,k-1)*dcz(k-1)+rho(i,j,k-2)*dcz(k-2))*(
         u(c,i,  j,k) -3*u(c,i,j,k-1)+ 3*u(c,i,j,k-2)- u(c,i,j,k-3)
      -(um(c,i,  j,k)-3*um(c,i,j,k-1)+3*um(c,i,j,k-2)-um(c,i,j,k-3)) )
					     )  );
	       }
	    }
      }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef strz
#undef dcz
#undef coz
   }
}

//-----------------------------------------------------------------------
void EW::addsgd4c_ci( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		      float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		      float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
		      float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
		      float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_cox,
		      float_sw4* __restrict__ a_coy, float_sw4 beta )
{
  SW4_MARK_FUNCTION;
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define jac(i,j,k) a_jac[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]

      const size_t ni   =     (ilast-ifirst+1);
      const size_t nij  =  ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
      
// #pragma omp parallel
//       {
//       for( int c=0 ; c < 3 ; c++ )
// #pragma omp for
//       for( int k=kfirst+2; k <= klast-2 ; k++ )
// 	 for( int j=jfirst+2; j <= jlast-2 ; j++ )
// #pragma simd
// #pragma ivdep
// 	    for( int i=ifirst+2; i <= ilast-2 ; i++ )
// 	    {
using ADDSGD_POL2  =
       RAJA::KernelPolicy<
       RAJA::statement::CudaKernel<
	 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>,
			      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<4>,
						   RAJA::statement::For<3, RAJA::cuda_threadblock_exec<32>,
									RAJA::statement::For<0, RAJA::seq_exec,
											     RAJA::statement::Lambda<0> >>>>>>;
RAJA::RangeSegment i_range(ifirst+2,ilast-1);
RAJA::RangeSegment j_range(jfirst+2,jlast-1);
RAJA::RangeSegment k_range(kfirst+2,klast-1);
RAJA::RangeSegment c_range(0,3);
RAJA::kernel<ADDSGD_POL2>(
			    RAJA::make_tuple(c_range,k_range,j_range,i_range),
			    [=]RAJA_DEVICE (int c,int k, int j,int i) {
			      
	       float_sw4 irhoj=beta/(rho(i,j,k)*jac(i,j,k));
	       {
		  up(c,i,j,k) -= irhoj*( 
		  // x-differences
		   strx(i)*coy(j)*(
		   rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*
                   ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
		   -2*rho(i,j,k)*dcx(i)*jac(i,j,k)*
                   ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
		   +rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*
                   ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
		   -rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*
                   (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
		   +2*rho(i,j,k)*dcx(i)*jac(i,j,k)*
                   (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
		   -rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*
                   (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
// y-differences
		   stry(j)*cox(i)*(
		    +rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*
                   ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
		    -2*rho(i,j,k)*dcy(j)*jac(i,j,k)*
                   ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
		    +rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*
                   ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
		    -rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*
                   (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
		    +2*rho(i,j,k)*dcy(j)*jac(i,j,k)*
                   (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
		    -rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*
		    (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) );
	       }
			    } ); SYNC_STREAM;
//}
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef jac
   }
}

//-----------------------------------------------------------------------
void EW::addsgd6c_ci(  int ifirst, int ilast, int jfirst, int jlast,
		       int kfirst, int klast,
		       float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		       float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		       float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
		       float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
		       float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_cox,
		       float_sw4* __restrict__ a_coy, float_sw4 beta )
{
  SW4_MARK_FUNCTION;
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define jac(i,j,k) a_jac[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
#pragma omp parallel
      {
      for( int c=0 ; c < 3 ; c++ )
#pragma omp for
      for( int k=kfirst+3; k <= klast-3 ; k++ )
	 for( int j=jfirst+3; j <= jlast-3 ; j++ )
#pragma simd
#pragma ivdep
	    for( int i=ifirst+3; i <= ilast-3 ; i++ )
	    {
	       float_sw4 birho=0.5*beta/(rho(i,j,k)*jac(i,j,k));
	       {
		 up(c,i,j,k) += birho*( 
       strx(i)*coy(j)*(
// x-differences
      (rho(i+2,j,k)*dcx(i+2)*jac(i+2,j,k)+rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k))*(
         u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
      -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
      -3*(rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)+rho(i,j,k)*dcx(i)*jac(i,j,k))*(
         u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
      -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
      +3*(rho(i,j,k)*dcx(i)*jac(i,j,k)+rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k))*(
         u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
      -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
      - (rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)+rho(i-2,j,k)*dcx(i-2)*jac(i-2,j,k))*(
         u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
      -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
		       ) +  stry(j)*cox(i)*( 
// y-differences
     (rho(i,j+2,k)*dcy(j+2)*jac(i,j+2,k)+rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k))*(
         u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
      -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
     -3*(rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)+rho(i,j,k)*dcy(j)*jac(i,j,k))*(
         u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
      -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
     +3*(rho(i,j,k)*dcy(j)*jac(i,j,k)+rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k))*(
         u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
      -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
     - (rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)+rho(i,j-2,k)*dcy(j-2)*jac(i,j-2,k))*(
         u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
      -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
					     )  );
	       }
	    }
      }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef jac
   }
}
