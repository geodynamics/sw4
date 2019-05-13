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

void scalar_prod_ci( int is, int ie, int js, int je, int ks, int ke,
		     int i1, int i2, int j1, int j2, int k1, int k2,
		     int onesided[6], float_sw4* __restrict__ a_u, 
		     float_sw4* __restrict__ a_v,
		     float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
		     float_sw4* __restrict__ a_strz, float_sw4& sc_prod )
{
   // is,ie,js,je,ks,ke declared size
   // i1,i2,j1,j2,k1,k2 domain over which scalar product is computed
   const int ni    = ie-is+1;
   const int nij   = ni*(je-js+1);
   const int nijk  = nij*(ke-ks+1);
   const int base  = -(is+ni*js+nij*ks);
   //   const int base3 = base-nijk;
   //#define u(c,i,j,k) a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
   //#define v(c,i,j,k) a_v[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
   //#define strx(i) a_strx[i-is]
   //#define stry(j) a_stry[j-js]
   //#define strz(k) a_strz[k-ks]
   const float_sw4 normwgh[4]={17.0/48,59.0/48,43.0/48,49.0/48};

   float_sw4 scprod_loc=0;
#pragma omp parallel for reduction(+:scprod_loc)
   for( int k=k1 ;k <= k2 ;k++ )
      for( int j=j1 ;j <= j2 ;j++ )
	 for( int i=i1 ;i <= i2 ; i++ )
	 {
	    size_t ind=base+i+ni*j+nij*k;
// NOTE: the scalar product is scaled by the stretching
//	    float_sw4 term =(u(1,i,j,k)*v(1,i,j,k) + u(2,i,j,k)*v(2,i,j,k) + u(3,i,j,k)*v(3,i,j,k))/
//	       (strx(i)*stry(j)*strz(k));
	    float_sw4 term =(a_u[ind]*a_v[ind]+a_u[ind+nijk]*a_v[ind+nijk] + a_u[ind+2*nijk]*a_v[ind+2*nijk])/
	       (a_strx[i-is]*a_stry[j-js]*a_strz[k-ks]);
	    float_sw4 normfact = 1;
           if( k <= 4 && onesided[4] == 1 )
              normfact = normwgh[k-1];
           if( k >= k2-3 && onesided[5] == 1 )
              normfact = normwgh[k2-k];
           scprod_loc += normfact*term;
	 }
   sc_prod = scprod_loc;
}
