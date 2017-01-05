#include <cmath>
#include "sw4.h"
#include "EW.h"

//-----------------------------------------------------------------------
void EW::gridinfo( int ib, int ie, int jb, int je, int kb, int ke,
		   float_sw4* met, float_sw4* jac, float_sw4&  minj,
		   float_sw4& maxj )
{
   // met not used, might be in the future
   maxj = -1e30;
   minj =  1e30;
   size_t npts = (static_cast<size_t>(ie-ib+1))*(je-jb+1)*(ke-kb+1);
   for( int i= 0 ; i < npts ; i++ )
   {
      maxj = jac[i]>maxj ? jac[i] : maxj;
      minj = jac[i]<minj ? jac[i] : minj;
   }
}

//-----------------------------------------------------------------------
int EW::metric( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* a_x,
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
int EW::metric_rev( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* a_x,
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

//-----------------------------------------------------------------------
void EW::metricexgh( int ib, int ie, int jb, int je, int kb, int ke,
		     int nz, float_sw4* a_x, float_sw4* a_y, float_sw4* a_z, 
		     float_sw4* a_met, float_sw4* a_jac, int order,
                     float_sw4 sb, float_sw4 zmax, float_sw4 amp, float_sw4 xc,
		     float_sw4 yc, float_sw4 xl, float_sw4 yl )
{
// Exact metric derivatives for the Gaussian hill topography

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
   double ixl2 = 1/(xl*xl);
   double iyl2 = 1/(yl*yl);
   for( int k = kb; k <= ke ; k++ )
      for( int j = jb; j <= je ; j++ )
	 for( int i = ib; i <= ie ; i++ )
	 {
	    double zp, zq, zr, zz;
	    double s = (k-1.0)/(nz-1.0);
	    if( s < sb )
	    {
	       double sdb = s/sb;
	       double tau  = amp*exp( - (x(i,j,1)-xc)*(x(i,j,1)-xc)*ixl2 
				      - (y(i,j,1)-yc)*(y(i,j,1)-yc)*iyl2 );
               double taup = -2*(x(i,j,1)-xc)*ixl2*tau;
               double tauq = -2*(y(i,j,1)-yc)*iyl2*tau;
               double p1 = 1-sdb;
               double p2 = 1;
	       double powvar = 1-sdb;
	       for( int l=2; l <= order-1; l++ )
	       {
		  //		  p1 = p1 + (1-sdb)**l;
		  //		  p2 = p2 + l*(1-sdb)**(l-1);
		  p2 += l*powvar;
		  powvar *= (1-sdb);
		  p1 += powvar;
	       }
               zp = taup*( -(1-sdb)+sdb*p1 );
               zq = tauq*( -(1-sdb)+sdb*p1);
               zr = (tau+zmax+(zmax+tau-h*sb*(nz-1))*p1 -
                            sdb*(zmax+tau-h*sb*(nz-1))*p2 )/sb;
               zz = (1-sdb)*(-tau) + 
		  sdb*(zmax+(zmax+tau-h*sb*(nz-1))*p1);
	    }
	    else
	    {
	       zp = 0;
	       zq = 0;
	       zr = h*(nz-1);
	       zz = zmax + (s-sb)*h*(nz-1);
	    }

 // Convert to 'undivided differences'
	    zp = zp*h;
	    zq = zq*h;
	    zr = zr/(nz-1);
                  
// Formulas from metric evaluation numerically
	    float_sw4 sqzr = sqrt(zr);
	    jac(i,j,k)   = h*h*zr;
	    met(1,i,j,k) = sqzr;
	    met(2,i,j,k) = -zp/sqzr;
	    met(3,i,j,k) = -zq/sqzr;
	    met(4,i,j,k) = h/sqzr;
	 }
#undef x
#undef y
#undef z
#undef jac
#undef met
}

//-----------------------------------------------------------------------
void EW::metricexgh_rev( int ib, int ie, int jb, int je, int kb, int ke,
			 int nz, float_sw4* a_x, float_sw4* a_y, float_sw4* a_z, 
			 float_sw4* a_met, float_sw4* a_jac, int order,
			 float_sw4 sb, float_sw4 zmax, float_sw4 amp, float_sw4 xc,
			 float_sw4 yc, float_sw4 xl, float_sw4 yl )
{
// Exact metric derivatives for the Gaussian hill topography

   const int ni    = ie-ib+1;
   const int nij   = ni*(je-jb+1);
   const int nijk  = nij*(ke-kb+1);
   const int base  = -(ib+ni*jb+nij*kb);
   const int base4 = base-nijk;
   const int nic  = 4*ni;
   const int nijc = 4*nij;
#define x(i,j,k)     a_x[base+i+ni*(j)+nij*(k)]
#define y(i,j,k)     a_y[base+i+ni*(j)+nij*(k)]
#define z(i,j,k)     a_z[base+i+ni*(j)+nij*(k)]
#define jac(i,j,k)   a_jac[base+i+ni*(j)+nij*(k)]
#define met(c,i,j,k) a_met[base4+(i)+ni*(j)+nij*(k)+nijk*(c)]

   double h = x(ib+1,jb,kb)-x(ib,jb,kb);
   double ixl2 = 1/(xl*xl);
   double iyl2 = 1/(yl*yl);
   for( int k = kb; k <= ke ; k++ )
      for( int j = jb; j <= je ; j++ )
	 for( int i = ib; i <= ie ; i++ )
	 {
	    double zp, zq, zr, zz;
	    double s = (k-1.0)/(nz-1.0);
	    if( s < sb )
	    {
	       double sdb = s/sb;
	       double tau  = amp*exp( - (x(i,j,1)-xc)*(x(i,j,1)-xc)*ixl2 
				      - (y(i,j,1)-yc)*(y(i,j,1)-yc)*iyl2 );
               double taup = -2*(x(i,j,1)-xc)*ixl2*tau;
               double tauq = -2*(y(i,j,1)-yc)*iyl2*tau;
               double p1 = 1-sdb;
               double p2 = 1;
	       double powvar = 1-sdb;
	       for( int l=2; l <= order-1; l++ )
	       {
		  //		  p1 = p1 + (1-sdb)**l;
		  //		  p2 = p2 + l*(1-sdb)**(l-1);
		  p2 += l*powvar;
		  powvar *= (1-sdb);
		  p1 += powvar;
	       }
               zp = taup*( -(1-sdb)+sdb*p1 );
               zq = tauq*( -(1-sdb)+sdb*p1);
               zr = (tau+zmax+(zmax+tau-h*sb*(nz-1))*p1 -
                            sdb*(zmax+tau-h*sb*(nz-1))*p2 )/sb;
               zz = (1-sdb)*(-tau) + 
		  sdb*(zmax+(zmax+tau-h*sb*(nz-1))*p1);
	    }
	    else
	    {
	       zp = 0;
	       zq = 0;
	       zr = h*(nz-1);
	       zz = zmax + (s-sb)*h*(nz-1);
	    }

 // Convert to 'undivided differences'
	    zp = zp*h;
	    zq = zq*h;
	    zr = zr/(nz-1);
                  
// Formulas from metric evaluation numerically
	    float_sw4 sqzr = sqrt(zr);
	    jac(i,j,k)   = h*h*zr;
	    met(1,i,j,k) = sqzr;
	    met(2,i,j,k) = -zp/sqzr;
	    met(3,i,j,k) = -zq/sqzr;
	    met(4,i,j,k) = h/sqzr;
	 }
#undef x
#undef y
#undef z
#undef jac
#undef met
}

//-----------------------------------------------------------------------
void EW::freesurfcurvisg( int ib, int ie, int jb, int je, int kb, int ke,
			  int nz, int side, float_sw4* a_u, float_sw4* a_mu,
			  float_sw4* a_la, float_sw4* a_met, float_sw4* s,
			  float_sw4* a_forcing, float_sw4* a_strx, float_sw4* a_stry )
{
   const float_sw4 c1=2.0/3, c2=-1.0/12;

   //      integer ifirst, ilast, jfirst, jlast, kfirst, klast
   //      integer i, j, k, kl, nz, side
   //      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
   //      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
   //      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
   //      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
   //      real*8 forcing(3,ifirst:ilast,jfirst:jlast)
   //      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
   //      real*8 s(0:4), rhs1, rhs2, rhs3, s0i, ac, bc, cc, dc
   //      real*8 istry, istrx, xoysqrt, yoxsqrt, isqrtxy

   const int ni    = ie-ib+1;
   const int nij   = ni*(je-jb+1);
   const int base  = -(ib+ni*jb+nij*kb);
   const int base4 = 4*base-1;
   const int base3 = 3*base-1;
   const int nic4  = 4*ni;
   const int nijc4 = 4*nij;
   const int nic3  = 3*ni;
   const int nijc3 = 3*nij;
   const int basef = -3*(ib+ni*jb)-1;
#define x(i,j,k)     a_x[base+i+ni*(j)+nij*(k)]
#define y(i,j,k)     a_y[base+i+ni*(j)+nij*(k)]
#define z(i,j,k)     a_z[base+i+ni*(j)+nij*(k)]
#define mu(i,j,k)    a_mu[base+i+ni*(j)+nij*(k)]
#define la(i,j,k)    a_la[base+i+ni*(j)+nij*(k)]
#define met(c,i,j,k) a_met[base4+c+4*(i)+nic4*(j)+nijc4*(k)]
#define u(c,i,j,k)   a_u[base3+c+3*(i)+nic3*(j)+nijc3*(k)]
#define forcing(c,i,j)   a_forcing[basef+(c)+3*(i)+nic3*(j)]
#define strx(i) a_strx[(i-ib)]
#define stry(j) a_stry[(j-jb)]

   int k, kl;
   if( side == 5 )
   {
      k = 1;
      kl= 1;
   }
   else if( side == 6 )
   {
      k = nz;
      kl= -1;
   }

   float_sw4 s0i = 1/s[0];
   for( int j= jb+2; j<=je-2 ; j++ )
   {
      float_sw4 istry = 1/stry(j);
      for( int i= ib+2; i<=ie-2 ; i++ )
      {
	 float_sw4 istrx = 1/strx(i);

    // First tangential derivatives
            float_sw4 rhs1 = 
// pr
        (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
               c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
               c1*(u(1,i+1,j,k)-u(1,i-1,j,k))  )*strx(i)*istry 
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
             c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
       + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
             c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*istry   
// qr
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
             c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )*istrx*stry(j) 
       + la(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )  -
	       forcing(1,i,j);

// (v-eq)
            float_sw4 rhs2 = 
// pr
         la(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
             c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   ) 
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
             c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  )*strx(i)*istry 
// qr
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
             c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   ) 
      + (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*stry(j)*istrx 
       + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
             c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*istrx -
	       forcing(2,i,j);

// (w-eq)
            float_sw4 rhs3 = 
// pr
         la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
             c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )*istry 
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
             c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*strx(i)*istry
// qr 
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
             c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*stry(j)*istrx
       + la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*istrx -
	       forcing(3,i,j);

// Normal derivatives
            float_sw4 ac = strx(i)*istry*met(2,i,j,k)*met(2,i,j,k)+
	       stry(j)*istrx*met(3,i,j,k)*met(3,i,j,k)+
	       met(4,i,j,k)*met(4,i,j,k)*istry*istrx;
            float_sw4 bc = 1/(mu(i,j,k)*ac);
            float_sw4 cc = (mu(i,j,k)+la(i,j,k))/(2*mu(i,j,k)+la(i,j,k))*bc/ac;

            float_sw4 xoysqrt = sqrt(strx(i)*istry);
            float_sw4 yoxsqrt = 1/xoysqrt;
            float_sw4 isqrtxy = istrx*xoysqrt;
            float_sw4 dc = cc*( xoysqrt*met(2,i,j,k)*rhs1 + 
		      yoxsqrt*met(3,i,j,k)*rhs2 + isqrtxy*met(4,i,j,k)*rhs3);

            u(1,i,j,k-kl) = -s0i*(  s[1]*u(1,i,j,k)+s[2]*u(1,i,j,k+kl)+
                s[3]*u(1,i,j,k+2*kl)+s[4]*u(1,i,j,k+3*kl) + bc*rhs1 - 
				    dc*met(2,i,j,k)*xoysqrt );
            u(2,i,j,k-kl) = -s0i*(  s[1]*u(2,i,j,k)+s[2]*u(2,i,j,k+kl)+
                s[3]*u(2,i,j,k+2*kl)+s[4]*u(2,i,j,k+3*kl) + bc*rhs2 - 
				    dc*met(3,i,j,k)*yoxsqrt );
            u(3,i,j,k-kl) = -s0i*(  s[1]*u(3,i,j,k)+s[2]*u(3,i,j,k+kl)+
                s[3]*u(3,i,j,k+2*kl)+s[4]*u(3,i,j,k+3*kl) + bc*rhs3 - 
				    dc*met(4,i,j,k)*isqrtxy );
	       }
   }
#undef x
#undef y
#undef z
#undef mu
#undef la
#undef met
#undef u
#undef forcing
#undef strx
#undef stry
}

//-----------------------------------------------------------------------
void EW::freesurfcurvisg_rev( int ib, int ie, int jb, int je, int kb, int ke,
			  int nz, int side, float_sw4* a_u, float_sw4* a_mu,
			  float_sw4* a_la, float_sw4* a_met, float_sw4* s,
			  float_sw4* a_forcing, float_sw4* a_strx, float_sw4* a_stry )
{
   const float_sw4 c1=2.0/3, c2=-1.0/12;
   //      integer ifirst, ilast, jfirst, jlast, kfirst, klast
   //      integer i, j, k, kl, nz, side
   //      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
   //      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
   //      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
   //      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
   //      real*8 forcing(3,ifirst:ilast,jfirst:jlast)
   //      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
   //      real*8 s(0:4), rhs1, rhs2, rhs3, s0i, ac, bc, cc, dc
   //      real*8 istry, istrx, xoysqrt, yoxsqrt, isqrtxy

   const int ni    = ie-ib+1;
   const int nij   = ni*(je-jb+1);
   const int nijk  = ni*(je-jb+1)*(ke-kb+1);
   const int base  = -(ib+ni*jb+nij*kb);
   const int basef = -(ib+ni*jb);
   const int base4 = base-nijk;
   const int base3 = base-nijk;
   //   const int nic4  = 4*ni;
   //   const int nijc4 = 4*nij;
   const int nic3  = 3*ni;
   //   const int nijc3 = 3*nij;
#define x(i,j,k)     a_x[base+i+ni*(j)+nij*(k)]
#define y(i,j,k)     a_y[base+i+ni*(j)+nij*(k)]
#define z(i,j,k)     a_z[base+i+ni*(j)+nij*(k)]
#define mu(i,j,k)    a_mu[base+i+ni*(j)+nij*(k)]
#define la(i,j,k)    a_la[base+i+ni*(j)+nij*(k)]
#define met(c,i,j,k) a_met[base4+(i)+ni*(j)+nij*(k)+nijk*(c)]
#define u(c,i,j,k)   a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]
#define forcing(c,i,j)   a_forcing[3*basef-1+(c)+3*(i)+nic3*(j)]
#define strx(i) a_strx[(i-ib)]
#define stry(j) a_stry[(j-jb)]

   int k, kl;
   if( side == 5 )
   {
      k = 1;
      kl= 1;
   }
   else if( side == 6 )
   {
      k = nz;
      kl= -1;
   }

   float_sw4 s0i = 1/s[0];
   for( int j= jb+2; j<=je-2 ; j++ )
   {
      float_sw4 istry = 1/stry(j);
      for( int i= ib+2; i<=ie-2 ; i++ )
      {
	 float_sw4 istrx = 1/strx(i);

    // First tangential derivatives
            float_sw4 rhs1 = 
// pr
        (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
               c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
               c1*(u(1,i+1,j,k)-u(1,i-1,j,k))  )*strx(i)*istry 
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
             c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
       + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
             c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*istry   
// qr
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
             c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )*istrx*stry(j) 
       + la(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )  -
	       forcing(1,i,j);

// (v-eq)
            float_sw4 rhs2 = 
// pr
         la(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
             c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   ) 
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
             c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  )*strx(i)*istry 
// qr
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
             c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   ) 
      + (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*stry(j)*istrx 
       + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
             c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*istrx -
	       forcing(2,i,j);

// (w-eq)
            float_sw4 rhs3 = 
// pr
         la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
             c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )*istry 
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
             c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*strx(i)*istry
// qr 
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
             c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*stry(j)*istrx
       + la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*istrx -
	       forcing(3,i,j);

// Normal derivatives
            float_sw4 ac = strx(i)*istry*met(2,i,j,k)*met(2,i,j,k)+
	       stry(j)*istrx*met(3,i,j,k)*met(3,i,j,k)+
	       met(4,i,j,k)*met(4,i,j,k)*istry*istrx;
            float_sw4 bc = 1/(mu(i,j,k)*ac);
            float_sw4 cc = (mu(i,j,k)+la(i,j,k))/(2*mu(i,j,k)+la(i,j,k))*bc/ac;

            float_sw4 xoysqrt = sqrt(strx(i)*istry);
            float_sw4 yoxsqrt = 1/xoysqrt;
            float_sw4 isqrtxy = istrx*xoysqrt;
            float_sw4 dc = cc*( xoysqrt*met(2,i,j,k)*rhs1 + 
		      yoxsqrt*met(3,i,j,k)*rhs2 + isqrtxy*met(4,i,j,k)*rhs3);

            u(1,i,j,k-kl) = -s0i*(  s[1]*u(1,i,j,k)+s[2]*u(1,i,j,k+kl)+
                s[3]*u(1,i,j,k+2*kl)+s[4]*u(1,i,j,k+3*kl) + bc*rhs1 - 
				    dc*met(2,i,j,k)*xoysqrt );
            u(2,i,j,k-kl) = -s0i*(  s[1]*u(2,i,j,k)+s[2]*u(2,i,j,k+kl)+
                s[3]*u(2,i,j,k+2*kl)+s[4]*u(2,i,j,k+3*kl) + bc*rhs2 - 
				    dc*met(3,i,j,k)*yoxsqrt );
            u(3,i,j,k-kl) = -s0i*(  s[1]*u(3,i,j,k)+s[2]*u(3,i,j,k+kl)+
                s[3]*u(3,i,j,k+2*kl)+s[4]*u(3,i,j,k+3*kl) + bc*rhs3 - 
				    dc*met(4,i,j,k)*isqrtxy );
	       }
   }
}
