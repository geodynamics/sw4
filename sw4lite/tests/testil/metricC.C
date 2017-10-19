#include "sw4.h"
#include <cmath>

//-----------------------------------------------------------------------
int metric( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* a_x,
		 float_sw4* a_y, float_sw4* a_z, float_sw4* a_met, float_sw4* a_jac )
{
   const float_sw4 c1=2.0/3, c2=-1.0/12;
   const float_sw4 fs= 5.0/6, ot=1.0/12, ft=4.0/3, os=1.0/6, d3=14.0/3;
   const int ni    = ie-ib+1;
   const int nij   = ni*(je-jb+1);
   const int base  = -(ib+ni*jb+nij*kb);
   const int base4 = 4*base-1;
   const int nic  = 4*ni;
   const int nijc = 4*nij;
#define x(i,j,k)     a_x[base+i+ni*(j)+nij*(k)]
#define y(i,j,k)     a_y[base+i+ni*(j)+nij*(k)]
#define z(i,j,k)     a_z[base+i+ni*(j)+nij*(k)]
#define jac(i,j,k)   a_jac[base+i+ni*(j)+nij*(k)]
#define met(c,i,j,k) a_met[base4+c+4*(i)+nic*(j)+nijc*(k)]

   double h = x(ib+1,jb,kb)-x(ib,jb,kb);
   for( int k = kb; k <= ke ; k++ )
      for( int j = jb; j <= je ; j++ )
	 for( int i = ib; i <= ie ; i++ )
	 {
    // k-derivatives
	    double zr, zp, zq, sqzr;
	    if( k >= kb+2 && k <= ke-2 )
                  zr = c2*(z(i,j,k+2)-z(i,j,k-2)) +
		     c1*(z(i,j,k+1)-z(i,j,k-1));
	    else if( k == kb )
	    {
	       zr=-2.25*z(i,j,k)+(4+fs)*z(i,j,k+1)-d3*z(i,j,k+2)+
		  3*z(i,j,k+3)-(1+ot)*z(i,j,k+4) +os*z(i,j,k+5);
	    }
	    else if( k == kb+1 )
	    {
	       zr = -os*z(i,j,k-1) -1.25*z(i,j,k)+(1+ft)*z(i,j,k+1)
		  - ft*z(i,j,k+2) + 0.5*z(i,j,k+3) -ot*z(i,j,k+4);
	    }
	    else if( k == ke-1 )
	    {
	       zr =  os*z(i,j,k+1) +1.25*z(i,j,k)-(1+ft)*z(i,j,k-1)
		  + ft*z(i,j,k-2) - 0.5*z(i,j,k-3) + ot*z(i,j,k-4);
	    }
	    else if( k == ke )
	    {
                  zr= 2.25*z(i,j,k)-(4+fs)*z(i,j,k-1)+d3*z(i,j,k-2)-
		           3*z(i,j,k-3)+(1+ot)*z(i,j,k-4) -os*z(i,j,k-5);
	    }
               
// j-derivatives
	    if( j >= jb+2 && j <= je-2 )
	    {
                  zq = c2*(z(i,j+2,k)-z(i,j-2,k)) + 
		     c1*(z(i,j+1,k)-z(i,j-1,k));
	    }
	    else if( j == jb )
	    {
                  zq=-2.25*z(i,j,k)+(4+fs)*z(i,j+1,k)-d3*z(i,j+2,k)+
		  3*z(i,j+3,k)-(1+ot)*z(i,j+4,k) +os*z(i,j+5,k);
	    }
	    else if( j == jb+1 )
	    {
                  zq = -os*z(i,j-1,k) -1.25*z(i,j,k)+(1+ft)*z(i,j+1,k)
		  - ft*z(i,j+2,k) + 0.5*z(i,j+3,k) -ot*z(i,j+4,k);
	    }
	    else if( j == je-1 )
	    {
                  zq = os*z(i,j+1,k) +1.25*z(i,j,k)-(1+ft)*z(i,j-1,k)
		  + ft*z(i,j-2,k) - 0.5*z(i,j-3,k) + ot*z(i,j-4,k);
	    }
	    else if( j == je )
	    {
                  zq= 2.25*z(i,j,k)-(4+fs)*z(i,j-1,k)+d3*z(i,j-2,k)-
		  3*z(i,j-3,k)+(1+ot)*z(i,j-4,k) -os*z(i,j-5,k);
	    }

// i-derivatives
	    if( i >= ib+2 && i <= ie-2 )
	    {
	       zp= c2*(z(i+2,j,k)-z(i-2,j,k)) + 
                     c1*(z(i+1,j,k)-z(i-1,j,k));
	    }
	    else if( i == ib )
	    {
	       zp=-2.25*z(i,j,k)+(4+fs)*z(i+1,j,k)-d3*z(i+2,j,k)+
		  3*z(i+3,j,k)-(1+ot)*z(i+4,j,k) +os*z(i+5,j,k);
	    }
	    else if( i == ib+1 )
	    {
	       zp = -os*z(i-1,j,k) -1.25*z(i,j,k)+(1+ft)*z(i+1,j,k)
		  - ft*z(i+2,j,k) + 0.5*z(i+3,j,k) - ot*z(i+4,j,k);
	    }
	    else if( i == ie-1)
	    {
	       zp =  os*z(i+1,j,k) +1.25*z(i,j,k)-(1+ft)*z(i-1,j,k)
		  + ft*z(i-2,j,k) - 0.5*z(i-3,j,k) + ot*z(i-4,j,k);
	    }
	    else if( i == ie)
	    {
	       zp= 2.25*z(i,j,k)-(4+fs)*z(i-1,j,k)+d3*z(i-2,j,k)-
		  3*z(i-3,j,k)+(1+ot)*z(i-4,j,k) -os*z(i-5,j,k);
	    }

// Compute the metric
	    if( zr <= 0 )
	       return -1;

	    sqzr = sqrt(zr);
	    jac(i,j,k) = h*h*zr;
	    met(1,i,j,k) = sqzr;
	    met(2,i,j,k) = -zp/sqzr;
	    met(3,i,j,k) = -zq/sqzr;
	    met(4,i,j,k) = h/sqzr;
	 }
   return 0;
#undef x
#undef y
#undef z
#undef jac
#undef met
}

//-----------------------------------------------------------------------
int metric_rev( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* a_x,
		float_sw4* a_y, float_sw4* a_z, float_sw4* a_met, float_sw4* a_jac )
{
   const float_sw4 c1=2.0/3, c2=-1.0/12;
   const float_sw4 fs= 5.0/6, ot=1.0/12, ft=4.0/3, os=1.0/6, d3=14.0/3;
   const int ni    = ie-ib+1;
   const int nij   = ni*(je-jb+1);
   const int nijk  = nij*(ke-kb+1);
   const int base  = -(ib+ni*jb+nij*kb);
   const int base4 = base-nijk;
   //   const int nic  = 4*ni;
   //   const int nijc = 4*nij;
#define x(i,j,k)     a_x[base+i+ni*(j)+nij*(k)]
#define y(i,j,k)     a_y[base+i+ni*(j)+nij*(k)]
#define z(i,j,k)     a_z[base+i+ni*(j)+nij*(k)]
#define jac(i,j,k)   a_jac[base+i+ni*(j)+nij*(k)]
#define met(c,i,j,k) a_met[base4+(i)+ni*(j)+nij*(k)+nijk*(c)]

   double h = x(ib+1,jb,kb)-x(ib,jb,kb);
   for( int k = kb; k <= ke ; k++ )
      for( int j = jb; j <= je ; j++ )
	 for( int i = ib; i <= ie ; i++ )
	 {
    // k-derivatives
	    double zr, zp, zq, sqzr;
	    if( k >= kb+2 && k <= ke-2 )
                  zr = c2*(z(i,j,k+2)-z(i,j,k-2)) +
		     c1*(z(i,j,k+1)-z(i,j,k-1));
	    else if( k == kb )
	    {
	       zr=-2.25*z(i,j,k)+(4+fs)*z(i,j,k+1)-d3*z(i,j,k+2)+
		  3*z(i,j,k+3)-(1+ot)*z(i,j,k+4) +os*z(i,j,k+5);
	    }
	    else if( k == kb+1 )
	    {
	       zr = -os*z(i,j,k-1) -1.25*z(i,j,k)+(1+ft)*z(i,j,k+1)
		  - ft*z(i,j,k+2) + 0.5*z(i,j,k+3) -ot*z(i,j,k+4);
	    }
	    else if( k == ke-1 )
	    {
	       zr =  os*z(i,j,k+1) +1.25*z(i,j,k)-(1+ft)*z(i,j,k-1)
		  + ft*z(i,j,k-2) - 0.5*z(i,j,k-3) + ot*z(i,j,k-4);
	    }
	    else if( k == ke )
	    {
                  zr= 2.25*z(i,j,k)-(4+fs)*z(i,j,k-1)+d3*z(i,j,k-2)-
		           3*z(i,j,k-3)+(1+ot)*z(i,j,k-4) -os*z(i,j,k-5);
	    }
               
// j-derivatives
	    if( j >= jb+2 && j <= je-2 )
	    {
                  zq = c2*(z(i,j+2,k)-z(i,j-2,k)) + 
		     c1*(z(i,j+1,k)-z(i,j-1,k));
	    }
	    else if( j == jb )
	    {
                  zq=-2.25*z(i,j,k)+(4+fs)*z(i,j+1,k)-d3*z(i,j+2,k)+
		  3*z(i,j+3,k)-(1+ot)*z(i,j+4,k) +os*z(i,j+5,k);
	    }
	    else if( j == jb+1 )
	    {
                  zq = -os*z(i,j-1,k) -1.25*z(i,j,k)+(1+ft)*z(i,j+1,k)
		  - ft*z(i,j+2,k) + 0.5*z(i,j+3,k) -ot*z(i,j+4,k);
	    }
	    else if( j == je-1 )
	    {
                  zq = os*z(i,j+1,k) +1.25*z(i,j,k)-(1+ft)*z(i,j-1,k)
		  + ft*z(i,j-2,k) - 0.5*z(i,j-3,k) + ot*z(i,j-4,k);
	    }
	    else if( j == je )
	    {
                  zq= 2.25*z(i,j,k)-(4+fs)*z(i,j-1,k)+d3*z(i,j-2,k)-
		  3*z(i,j-3,k)+(1+ot)*z(i,j-4,k) -os*z(i,j-5,k);
	    }

// i-derivatives
	    if( i >= ib+2 && i <= ie-2 )
	    {
	       zp= c2*(z(i+2,j,k)-z(i-2,j,k)) + 
                     c1*(z(i+1,j,k)-z(i-1,j,k));
	    }
	    else if( i == ib )
	    {
	       zp=-2.25*z(i,j,k)+(4+fs)*z(i+1,j,k)-d3*z(i+2,j,k)+
		  3*z(i+3,j,k)-(1+ot)*z(i+4,j,k) +os*z(i+5,j,k);
	    }
	    else if( i == ib+1 )
	    {
	       zp = -os*z(i-1,j,k) -1.25*z(i,j,k)+(1+ft)*z(i+1,j,k)
		  - ft*z(i+2,j,k) + 0.5*z(i+3,j,k) - ot*z(i+4,j,k);
	    }
	    else if( i == ie-1)
	    {
	       zp =  os*z(i+1,j,k) +1.25*z(i,j,k)-(1+ft)*z(i-1,j,k)
		  + ft*z(i-2,j,k) - 0.5*z(i-3,j,k) + ot*z(i-4,j,k);
	    }
	    else if( i == ie)
	    {
	       zp= 2.25*z(i,j,k)-(4+fs)*z(i-1,j,k)+d3*z(i-2,j,k)-
		  3*z(i-3,j,k)+(1+ot)*z(i-4,j,k) -os*z(i-5,j,k);
	    }

// Compute the metric
	    if( zr <= 0 )
	       return -1;

	    sqzr = sqrt(zr);
	    jac(i,j,k) = h*h*zr;
	    met(1,i,j,k) = sqzr;
	    met(2,i,j,k) = -zp/sqzr;
	    met(3,i,j,k) = -zq/sqzr;
	    met(4,i,j,k) = h/sqzr;
	 }
   return 0;
#undef x
#undef y
#undef z
#undef jac
#undef met
}
