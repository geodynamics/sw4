#include <cmath>
#include "sw4.h"
#include "EW.h"

//-----------------------------------------------------------------------
void EW::gridinfo_ci( int ib, int ie, int jb, int je, int kb, int ke,
		   float_sw4* __restrict__ met, float_sw4* __restrict__ jac, 
		   float_sw4&  minj, float_sw4& maxj )
{
   // met not used, might be in the future
   maxj = -1e30;
   minj =  1e30;
   size_t npts = (static_cast<size_t>(ie-ib+1))*(je-jb+1)*(ke-kb+1);
#pragma omp parallel for reduction(max:maxj) reduction(min:minj)
   for( int i= 0 ; i < npts ; i++ )
   {
      maxj = jac[i]>maxj ? jac[i] : maxj;
      minj = jac[i]<minj ? jac[i] : minj;
   }
}

//-----------------------------------------------------------------------
int EW::metric_ci( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* __restrict__ a_x,
		   float_sw4* __restrict__ a_y, float_sw4* __restrict__ a_z,
		   float_sw4* __restrict__ a_met, float_sw4* __restrict__ a_jac )
{
   const float_sw4 c1=2.0/3, c2=-1.0/12;
   const float_sw4 fs= 5.0/6, ot=1.0/12, ft=4.0/3, os=1.0/6, d3=14.0/3;
   const int ni    = ie-ib+1;
   const int nij   = ni*(je-jb+1);
   const int nijk  = nij*(ke-kb+1);
   const int base  = -(ib+ni*jb+nij*kb);
   const int base4 = base-nijk;
   int ecode = 0;
   //   const int nic  = 4*ni;
   //   const int nijc = 4*nij;
#define x(i,j,k)     a_x[base+i+ni*(j)+nij*(k)]
#define y(i,j,k)     a_y[base+i+ni*(j)+nij*(k)]
#define z(i,j,k)     a_z[base+i+ni*(j)+nij*(k)]
#define jac(i,j,k)   a_jac[base+i+ni*(j)+nij*(k)]
#define met(c,i,j,k) a_met[base4+(i)+ni*(j)+nij*(k)+nijk*(c)]

   double h = x(ib+1,jb,kb)-x(ib,jb,kb);

#pragma omp parallel for reduction(+:ecode)
   for( int k = kb; k <= ke ; k++ )
      for( int j = jb; j <= je ; j++ )
#pragma ivdep
#pragma simd
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
	    {
	       ecode = -1;
	       //	       cout << "zr = " << zr << " at " << i << " " << j << " " << k << endl;
	       //	       cout << "x,y,z = " << x(i,j,k) << " " << y(i,j,k) << " " << z(i,j,k) << endl;
	    }
	    sqzr = sqrt(zr);
	    jac(i,j,k) = h*h*zr;
	    met(1,i,j,k) = sqzr;
	    met(2,i,j,k) = -zp/sqzr;
	    met(3,i,j,k) = -zq/sqzr;
	    met(4,i,j,k) = h/sqzr;
	 }
   return ecode;
#undef x
#undef y
#undef z
#undef jac
#undef met
}

//-----------------------------------------------------------------------
void EW::metricexgh_ci( int ib, int ie, int jb, int je, int kb, int ke,
			int nz, float_sw4* __restrict__ a_x, float_sw4* __restrict__ a_y, 
                        float_sw4* __restrict__ a_z, float_sw4* __restrict__ a_met, 
			float_sw4* __restrict__ a_jac, int order,
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
#pragma omp parallel for
   for( int k = kb; k <= ke ; k++ )
      for( int j = jb; j <= je ; j++ )
#pragma ivdep
#pragma simd
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
void EW::freesurfcurvi_ci( int ib, int ie, int jb, int je, int kb, int ke,
			   int nz, int side, float_sw4* __restrict__ a_u, 
			   float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_la, 
			   float_sw4* __restrict__ a_met, float_sw4* s,
			   float_sw4* __restrict__ a_forcing )
{
   const float_sw4 c1=2.0/3, c2=-1.0/12;

   const int ni    = ie-ib+1;
   const int nij   = ni*(je-jb+1);
   const int nijk  = ni*(je-jb+1)*(ke-kb+1);
   const int base  = -(ib+ni*jb+nij*kb);
   const int basef = -(ib+ni*jb);
   const int base4 = base-nijk;
   const int base3 = base-nijk;
   const int nic3  = 3*ni;
   //#define x(i,j,k)     a_x[base+i+ni*(j)+nij*(k)]
   //#define y(i,j,k)     a_y[base+i+ni*(j)+nij*(k)]
   //#define z(i,j,k)     a_z[base+i+ni*(j)+nij*(k)]
#define mu(i,j,k)    a_mu[base+i+ni*(j)+nij*(k)]
#define la(i,j,k)    a_la[base+i+ni*(j)+nij*(k)]
#define met(c,i,j,k) a_met[base4+(i)+ni*(j)+nij*(k)+nijk*(c)]
#define u(c,i,j,k)   a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]
#define forcing(c,i,j)   a_forcing[3*basef-1+(c)+3*(i)+nic3*(j)]

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
#pragma omp parallel for
   for( int j= jb+2; j<=je-2 ; j++ )
   {
#pragma ivdep
#pragma simd
      for( int i= ib+2; i<=ie-2 ; i++ )
      {
    // First tangential derivatives
            float_sw4 rhs1 = 
// pr
        (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
               c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
               c1*(u(1,i+1,j,k)-u(1,i-1,j,k))  )
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
             c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
       + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
             c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )
// qr
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
             c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )
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
             c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  )
// qr
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
             c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   ) 
      + (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )
       + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
             c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   ) -
	       forcing(2,i,j);

// (w-eq)
            float_sw4 rhs3 = 
// pr
         la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
             c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
             c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )
// qr 
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
             c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )
       + la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  ) -
	       forcing(3,i,j);

// Normal derivatives
            float_sw4 ac = met(2,i,j,k)*met(2,i,j,k)+ met(3,i,j,k)*met(3,i,j,k)+
	         met(4,i,j,k)*met(4,i,j,k);
            float_sw4 bc = 1/(mu(i,j,k)*ac);
            float_sw4 cc = (mu(i,j,k)+la(i,j,k))/(2*mu(i,j,k)+la(i,j,k))*bc/ac;
            float_sw4 dc = cc*( met(2,i,j,k)*rhs1 + met(3,i,j,k)*rhs2 
				+ met(4,i,j,k)*rhs3 );

            u(1,i,j,k-kl) = -s0i*(  s[1]*u(1,i,j,k)+s[2]*u(1,i,j,k+kl)+
                s[3]*u(1,i,j,k+2*kl)+s[4]*u(1,i,j,k+3*kl) + bc*rhs1 - 
				    dc*met(2,i,j,k) );
            u(2,i,j,k-kl) = -s0i*(  s[1]*u(2,i,j,k)+s[2]*u(2,i,j,k+kl)+
                s[3]*u(2,i,j,k+2*kl)+s[4]*u(2,i,j,k+3*kl) + bc*rhs2 - 
				    dc*met(3,i,j,k) );
            u(3,i,j,k-kl) = -s0i*(  s[1]*u(3,i,j,k)+s[2]*u(3,i,j,k+kl)+
                s[3]*u(3,i,j,k+2*kl)+s[4]*u(3,i,j,k+3*kl) + bc*rhs3 - 
				    dc*met(4,i,j,k) );
      }
   }
   //#undef x
   //#undef y
   //#undef z
#undef mu
#undef la
#undef met
#undef u
#undef forcing
}

//-----------------------------------------------------------------------
void EW::getsurfforcing_ci( int ifirst, int ilast, int jfirst, int jlast,
			    int kfirst, int klast, int k, float_sw4* __restrict__ a_met,
			    float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_tau,
			    float_sw4* __restrict__ a_forcing )
{
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = ni*(jlast-jfirst+1)*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int basef = -(ifirst+ni*jfirst);
   const int base3 = base-nijk;
   const int basef3= basef-nij;
   const int nic3  = 3*ni;
   const int nic6  = 6*ni;

#define met(c,i,j,k)   a_met[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]
#define jac(i,j,k)     a_jac[base+(i)+ni*(j)+nij*(k)]
#define forcing(c,i,j) a_forcing[3*basef-1+(c)+3*(i)+nic3*(j)]
#define tau(c,i,j)     a_tau[basef-nij +(i)+ni*(j)+nij*(c)]
   //#define tau(c,i,j)         a_tau[6*basef-1+(c)+6*(i)+nic6*(j)]

#pragma omp parallel for
   for( int j=jfirst ; j <= jlast ; j++ )
#pragma ivdep
#pragma simd
      for( int i=ifirst ; i <=ilast ; i++ )
      {
	 float_sw4 sqjac = sqrt(jac(i,j,k));
	 forcing(1,i,j) =  sqjac*( met(2,i,j,k)*tau(1,i,j)+
	   met(3,i,j,k)*tau(2,i,j)+met(4,i,j,k)*tau(3,i,j) );
	 forcing(2,i,j) =  sqjac*( met(2,i,j,k)*tau(2,i,j)+
	   met(3,i,j,k)*tau(4,i,j)+met(4,i,j,k)*tau(5,i,j) );
	 forcing(3,i,j) =  sqjac*( met(2,i,j,k)*tau(3,i,j)+
           met(3,i,j,k)*tau(5,i,j)+met(4,i,j,k)*tau(6,i,j) );
      }
#undef met
#undef jac
#undef forcing
#undef tau
}

//-----------------------------------------------------------------------
void EW::getsurfforcinggh_ci( int ifirst, int ilast, int jfirst, int jlast,
			      int kfirst, int klast, int k, float_sw4 h, 
			      float_sw4* __restrict__ a_tau,
			      float_sw4* __restrict__ a_forcing, float_sw4 amp,
			      float_sw4 xc, float_sw4 yc, float_sw4 xl, float_sw4 yl )
{
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = ni*(jlast-jfirst+1)*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int basef = -(ifirst+ni*jfirst);
   const int base3 = base-nijk;
   const int basef3= basef-nij;
   const int nic3  = 3*ni;
   const int nic6  = 6*ni;

#define forcing(c,i,j) a_forcing[3*basef-1+(c)+3*(i)+nic3*(j)]
#define tau(c,i,j)     a_tau[basef-nij +(i)+ni*(j)+nij*(c)]
//#define tau(c,i,j)         a_tau[6*basef-1+(c)+6*(i)+nic6*(j)]
   float_sw4 ixl2 = 1/(xl*xl);
   float_sw4 iyl2 = 1/(yl*yl);
#pragma omp parallel for
   for( int j=jfirst ; j <= jlast ; j++ )
   {
      float_sw4 y = (j-1)*h;
      for( int i=ifirst ; i <=ilast ; i++ )
      {
	 float_sw4 x = (i-1)*h;
         float_sw4 efact = amp*exp(-(x-xc)*(x-xc)*ixl2 - (y-yc)*(y-yc)*iyl2 );
         float_sw4 zp = 2*(x-xc)*ixl2*efact;
         float_sw4 zq = 2*(y-yc)*iyl2*efact;
	 forcing(1,i,j) =  h*h*( -zp*tau(1,i,j)-zq*tau(2,i,j)+tau(3,i,j) );
	 forcing(2,i,j) =  h*h*( -zp*tau(2,i,j)-zq*tau(4,i,j)+tau(5,i,j) );
	 forcing(3,i,j) =  h*h*( -zp*tau(3,i,j)-zq*tau(5,i,j)+tau(6,i,j) );
      }
#undef forcing
#undef tau
   }
}

//-----------------------------------------------------------------------
void EW::subsurfforcing_ci( int ifirst, int ilast, int jfirst, int jlast,
			    int kfirst, int klast, int k, float_sw4* __restrict__ a_met,
			    float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_tau,
			    float_sw4* __restrict__ a_forcing )
{
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = ni*(jlast-jfirst+1)*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int basef = -(ifirst+ni*jfirst);
   const int base3 = base-nijk;
   const int basef3= basef-nij;
   const int nic3  = 3*ni;
   const int nic6  = 6*ni;

#define met(c,i,j,k)   a_met[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]
#define jac(i,j,k)     a_jac[base+(i)+ni*(j)+nij*(k)]
#define forcing(c,i,j) a_forcing[3*basef-1+(c)+3*(i)+nic3*(j)]
#define tau(c,i,j)     a_tau[basef-nij +(i)+ni*(j)+nij*(c)]
   //#define tau(c,i,j)         a_tau[6*basef-1+(c)+6*(i)+nic6*(j)]

#pragma omp parallel for
   for( int j=jfirst ; j <= jlast ; j++ )
      for( int i=ifirst ; i <=ilast ; i++ )
      {
	 float_sw4 sqjac = sqrt(jac(i,j,k));
	 forcing(1,i,j) -=  sqjac*( met(2,i,j,k)*tau(1,i,j)+
	   met(3,i,j,k)*tau(2,i,j)+met(4,i,j,k)*tau(3,i,j) );
	 forcing(2,i,j) -=  sqjac*( met(2,i,j,k)*tau(2,i,j)+
	   met(3,i,j,k)*tau(4,i,j)+met(4,i,j,k)*tau(5,i,j) );
	 forcing(3,i,j) -=  sqjac*( met(2,i,j,k)*tau(3,i,j)+
           met(3,i,j,k)*tau(5,i,j)+met(4,i,j,k)*tau(6,i,j) );
      }

#undef met
#undef jac
#undef forcing
#undef tau
}

//-----------------------------------------------------------------------
void EW::addbstressc_ci( int ifirst, int ilast, int jfirst, int jlast,
			 int kfirst, int klast, int nz, float_sw4* __restrict__ a_u, 
			 float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_la,
			 float_sw4* __restrict__ a_bs, float_sw4* __restrict__ a_met,
			 int side, float_sw4* s, char op, int ghterm, int usesg,
			 float_sw4* __restrict__ a_sgstrx, float_sw4* __restrict__ a_sgstry )
{
   const float_sw4 c1=2.0/3, c2=-1.0/12;

   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = ni*(jlast-jfirst+1)*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int basef = -(ifirst+ni*jfirst);
   const int base3 = base-nijk;
   const int basef3= basef-nij;

#define mu(i,j,k)    a_mu[base+i+ni*(j)+nij*(k)]
#define la(i,j,k)    a_la[base+i+ni*(j)+nij*(k)]
#define met(c,i,j,k) a_met[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]
#define u(c,i,j,k)   a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]
#define bs(c,i,j)   a_bs[basef3+(i)+ni*(j) + nij*(c)]
#define sgstrx(i) a_sgstrx[i-ifirst]
#define sgstry(j) a_sgstry[j-jfirst]

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
   int a1, a2;
   if( op == '=' )
   {
      a1=0; a2=1;
   }
   else if( op == '+' )
   {
      a1=1; a2=1;
   }
   else if( op == '-' )
   {
      a1=1; a2=-1;
   }

#pragma omp parallel for
   for( int j= jfirst+2; j <= jlast-2 ; j++ )
   {
      float_sw4 sgy  = usesg ? sgstry(j) : 1 ;
      float_sw4 isgy = 1/sgy;
#pragma ivdep
#pragma simd
      for( int i= ifirst+2; i<=ilast-2 ; i++ )
      {
	 float_sw4 sgx = usesg ? sgstrx(i) : 1;
	 float_sw4 isgx = 1/sgx;
// First, tangential derivatives
         float_sw4 rhs1 = 
 // pr
           (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
               c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
               c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )*sgx*isgy 
       + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
             c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
       + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
             c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*isgy   
// qr
       +    mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
             c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )*isgx*sgy 
       + la(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  );

// (v-eq)
	 float_sw4 rhs2 = 
 // pr
          la(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
             c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   ) 
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
             c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  )*sgx*isgy 
// qr
       +    mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
             c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )
      + (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*sgy*isgx 
       + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
             c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*isgx;
// (w-eq)
	 float_sw4 rhs3 = 
// pr
           la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
             c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )*isgy 
       + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
             c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*sgx*isgy 
// qr 
       +    mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
             c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
             c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*sgy*isgx 
       + la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
             c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
             c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*isgx;
     
// then, normal derivatives
           float_sw4 un1 = s[1]*u(1,i,j,k)+s[2]*u(1,i,j,k+kl)+s[3]*u(1,i,j,k+2*kl)
	      +s[4]*u(1,i,j,k+3*kl);
           float_sw4 vn1 = s[1]*u(2,i,j,k)+s[2]*u(2,i,j,k+kl)+s[3]*u(2,i,j,k+2*kl)
	      +s[4]*u(2,i,j,k+3*kl);
           float_sw4 wn1 = s[1]*u(3,i,j,k)+s[2]*u(3,i,j,k+kl)+s[3]*u(3,i,j,k+2*kl)
	      +s[4]*u(3,i,j,k+3*kl);
           if( ghterm == 1 )
	   {
              un1 = un1 + s[0]*u(1,i,j,k-kl);
              vn1 = vn1 + s[0]*u(2,i,j,k-kl);
              wn1 = wn1 + s[0]*u(3,i,j,k-kl);
           }
           float_sw4 m2sg = sqrt(sgx*isgy);
           float_sw4 m3sg = 1/m2sg;
           float_sw4 m4sg = isgx*m2sg;
           float_sw4 rtu = un1*m2sg*met(2,i,j,k) + vn1*m3sg*met(3,i,j,k) +
	                   wn1*m4sg*met(4,i,j,k);
           float_sw4 ac  = sgx*isgy*met(2,i,j,k)*met(2,i,j,k) 
	               +   sgy*isgx*met(3,i,j,k)*met(3,i,j,k) 
	               +  isgx*isgy*met(4,i,j,k)*met(4,i,j,k);
           rhs1 = rhs1 + (mu(i,j,k)+la(i,j,k))*rtu*m2sg*met(2,i,j,k) +
	      mu(i,j,k)*ac*un1;
           rhs2 = rhs2 + (mu(i,j,k)+la(i,j,k))*rtu*m3sg*met(3,i,j,k) +
	      mu(i,j,k)*ac*vn1;
           rhs3 = rhs3 + (mu(i,j,k)+la(i,j,k))*rtu*m4sg*met(4,i,j,k) +
	      mu(i,j,k)*ac*wn1;
           bs(1,i,j) = a1*bs(1,i,j) + a2*rhs1;
           bs(2,i,j) = a1*bs(2,i,j) + a2*rhs2;
           bs(3,i,j) = a1*bs(3,i,j) + a2*rhs3;
      }
#undef mu
#undef la
#undef met
#undef u
#undef bs
#undef sgstrx
#undef sgstry
   }
}
