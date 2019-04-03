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

void energy4_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		 int i1, int i2, int j1, int j2, int k1, int k2, int* onesided,
		 float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_up,
		 float_sw4* __restrict__ a_rho, float_sw4 h, float_sw4* a_strx, float_sw4* a_stry,
		 float_sw4* a_strz, float_sw4& a_energy )
{
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = nij*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int base3 = base-nijk;
#define rho(i,j,k)  a_rho[base +(i)+ni*(j)+nij*(k)]   
#define um(c,i,j,k)  a_um[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define u(c,i,j,k)    a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define up(c,i,j,k)  a_up[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define strx(i) a_strx[i-ifirst]
#define stry(j) a_stry[j-jfirst]
#define strz(k) a_strz[k-kfirst]

   const float_sw4 normwgh[4]={17.0/48,59.0/48,43.0/48,49.0/48};
   float_sw4 energy=0;
#pragma omp parallel for  reduction(+:energy)
   for( int k = k1; k <= k2 ; k++ )
	 for( int j = j1; j <= j2 ; j++ )
	    //#pragma simd
#pragma ivdep	 
	    for( int i = i1; i <= i2 ; i++ )
	    {
               float_sw4 term =(
                 (up(1,i,j,k)-u(1,i,j,k))*(up(1,i,j,k)-u(1,i,j,k)) + 
                 (up(2,i,j,k)-u(2,i,j,k))*(up(2,i,j,k)-u(2,i,j,k)) +
                 (up(3,i,j,k)-u(3,i,j,k))*(up(3,i,j,k)-u(3,i,j,k)) -
                   up(1,i,j,k)*(up(1,i,j,k)-2*u(1,i,j,k)+um(1,i,j,k)) -
                   up(2,i,j,k)*(up(2,i,j,k)-2*u(2,i,j,k)+um(2,i,j,k)) -
                   up(3,i,j,k)*(up(3,i,j,k)-2*u(3,i,j,k)+um(3,i,j,k))
		 //                   )*rho(i,j,k);
  	       )*rho(i,j,k)/(strx(i)*stry(j)*strz(k));
	       float_sw4 normfact = 1;
               if( k <= 4 && onesided[4] == 1 )
                  normfact = normwgh[k-1];
               if( k >= k2-3 && onesided[5] == 1 )
                  normfact = normwgh[k2-k];
               energy += normfact*h*h*h*term;
	    }
   a_energy = energy;
}
#undef u
#undef up
#undef um
#undef rho
#undef strx
#undef stry
#undef strz
//-----------------------------------------------------------------------
void energy4c_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		 int i1, int i2, int j1, int j2, int k1, int k2, int* onesided,
		 float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_up,
		 float_sw4* __restrict__ a_rho, float_sw4* __restrict__ a_jac, float_sw4& a_energy )
{
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = nij*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int base3 = base-nijk;
#define rho(i,j,k)  a_rho[base +(i)+ni*(j)+nij*(k)]   
#define jac(i,j,k)  a_jac[base +(i)+ni*(j)+nij*(k)]   
#define um(c,i,j,k)  a_um[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define u(c,i,j,k)    a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define up(c,i,j,k)  a_up[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
   const float_sw4 normwgh[4]={17.0/48,59.0/48,43.0/48,49.0/48};
   float_sw4 energy=0;
#pragma omp parallel for reduction(+:energy)
   for( int k = k1; k <= k2 ; k++ )
	 for( int j = j1; j <= j2 ; j++ )
	    //#pragma simd
#pragma ivdep	 
	    for( int i = i1; i <= i2 ; i++ )
	    {
               float_sw4 term =(
                 (up(1,i,j,k)-u(1,i,j,k))*(up(1,i,j,k)-u(1,i,j,k)) + 
                 (up(2,i,j,k)-u(2,i,j,k))*(up(2,i,j,k)-u(2,i,j,k)) +
                 (up(3,i,j,k)-u(3,i,j,k))*(up(3,i,j,k)-u(3,i,j,k)) -
                   up(1,i,j,k)*(up(1,i,j,k)-2*u(1,i,j,k)+um(1,i,j,k)) -
                   up(2,i,j,k)*(up(2,i,j,k)-2*u(2,i,j,k)+um(2,i,j,k)) -
                   up(3,i,j,k)*(up(3,i,j,k)-2*u(3,i,j,k)+um(3,i,j,k))
				)*rho(i,j,k)*jac(i,j,k);
	       float_sw4 normfact = 1;
               if( k <= 4 && onesided[4] == 1 )
                  normfact = normwgh[k-1];
               if( k >= k2-3 && onesided[5] == 1 )
                  normfact = normwgh[k2-k];
               energy += normfact*term;
	    }
   a_energy = energy;
}
