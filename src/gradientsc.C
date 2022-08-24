#include "sw4.h"

void addgradrho_ci(  int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		     int ifirstact, int ilastact, int jfirstact, int jlastact, 
		     int kfirstact, int klastact, int nk,
		     float_sw4* __restrict__ a_kap, float_sw4* __restrict__ a_kapacc, 
		     float_sw4* __restrict__ a_um,  float_sw4* __restrict__ a_u,
		     float_sw4* __restrict__ a_up,  float_sw4* __restrict__ a_uacc,
		     float_sw4* __restrict__ a_grho,
		     float_sw4 dt, float_sw4 h, int onesided[6])
{
   const float_sw4 idt= 1.0/dt;
   const float_sw4 dt2o12 = dt*dt/12;
   const float_sw4 h3 = h*h*h;
   const float_sw4 normwgh[4]={17.0/48,59.0/48,43.0/48,49.0/48};   
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = nij*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int base3 = base-nijk;
#define grho(i,j,k) a_grho[base +(i)+ni*(j)+nij*(k)]   
#define um(c,i,j,k)   a_um[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define u(c,i,j,k)     a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define up(c,i,j,k)   a_up[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define kap(c,i,j,k) a_kap[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define kapacc(c,i,j,k) a_kapacc[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define uacc(c,i,j,k)     a_uacc[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#pragma omp parallel for
   for( int k = kfirstact; k <= klastact ; k++ )
      for( int j = jfirstact; j <= jlastact ; j++ )
#pragma ivdep	 
	 for( int i = ifirstact; i <= ilastact ; i++ )
	 {
	    float_sw4 normfact = h3;
	    if( k <= 4 && onesided[4] == 1 )
	       normfact *= normwgh[k-1];
            if( k >= nk-3 && onesided[5] == 1 )
               normfact *= normwgh[nk-k];
	    grho(i,j,k) +=  (
          kap(1,i,j,k)*(up(1,i,j,k)-2*u(1,i,j,k)+um(1,i,j,k))*idt*idt
       + dt2o12*kapacc(1,i,j,k)*uacc(1,i,j,k) +
          kap(2,i,j,k)*(up(2,i,j,k)-2*u(2,i,j,k)+um(2,i,j,k))*idt*idt
       + dt2o12*kapacc(2,i,j,k)*uacc(2,i,j,k) +
          kap(3,i,j,k)*(up(3,i,j,k)-2*u(3,i,j,k)+um(3,i,j,k))*idt*idt
	  + dt2o12*kapacc(3,i,j,k)*uacc(3,i,j,k) )*normfact;
	 }
}

//-----------------------------------------------------------------------
void addgradrhoc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		     int ifirstact, int ilastact, int jfirstact, int jlastact, 
		     int kfirstact, int klastact, int nk,
		     float_sw4* __restrict__ a_kap, float_sw4* __restrict__ a_kapacc, 
		     float_sw4* __restrict__ a_um,  float_sw4* __restrict__ a_u,
		     float_sw4* __restrict__ a_up,  float_sw4* __restrict__ a_uacc,
		     float_sw4* __restrict__ a_grho,
		     float_sw4 dt, float_sw4* __restrict__ a_jac, int onesided[6])
{
   const float_sw4 idt= 1.0/dt;
   const float_sw4 dt2o12 = dt*dt/12;
   const float_sw4 normwgh[4]={17.0/48,59.0/48,43.0/48,49.0/48};   
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = nij*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int base3 = base-nijk;
#define grho(i,j,k) a_grho[base +(i)+ni*(j)+nij*(k)]   
#define jac(i,j,k)   a_jac[base +(i)+ni*(j)+nij*(k)]   
#define um(c,i,j,k)   a_um[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define u(c,i,j,k)     a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define up(c,i,j,k)   a_up[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define kap(c,i,j,k) a_kap[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define kapacc(c,i,j,k) a_kapacc[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define uacc(c,i,j,k)     a_uacc[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#pragma omp parallel for
   for( int k = kfirstact; k <= klastact ; k++ )
      for( int j = jfirstact; j <= jlastact ; j++ )
#pragma ivdep	 
	 for( int i = ifirstact; i <= ilastact ; i++ )
	 {
	    float_sw4 normfact = jac(i,j,k);
	    if( k <= 4 && onesided[4] == 1 )
	       normfact *= normwgh[k-1];
            if( k >= nk-3 && onesided[5] == 1 )
               normfact *= normwgh[nk-k];
	    grho(i,j,k) +=  (
          kap(1,i,j,k)*(up(1,i,j,k)-2*u(1,i,j,k)+um(1,i,j,k))*idt*idt
       + dt2o12*kapacc(1,i,j,k)*uacc(1,i,j,k) +
          kap(2,i,j,k)*(up(2,i,j,k)-2*u(2,i,j,k)+um(2,i,j,k))*idt*idt
       + dt2o12*kapacc(2,i,j,k)*uacc(2,i,j,k) +
          kap(3,i,j,k)*(up(3,i,j,k)-2*u(3,i,j,k)+um(3,i,j,k))*idt*idt
	  + dt2o12*kapacc(3,i,j,k)*uacc(3,i,j,k) )*normfact;
	 }
}

//-----------------------------------------------------------------------
void addgradmula_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		     int ifirstact, int ilastact, int jfirstact, int jlastact, 
		     int kfirstact, int klastact, int nk,
		     float_sw4* __restrict__ a_kap, float_sw4* __restrict__ a_kapacc, 
		     float_sw4* __restrict__ a_u,   float_sw4* __restrict__ a_uacc,
		     float_sw4* __restrict__ a_gmu, float_sw4* __restrict__ a_glambda,
		     float_sw4 dt, float_sw4 h, int onesided[6],
		     int nb, int wb, float_sw4* __restrict__ a_bop )
{
   // nk is number of points, not counting ghost points.
   const float_sw4 h3 = h*h*h;
   const float_sw4 ih2= 1.0/(h*h);
   const float_sw4 dt2o12 = dt*dt/12;
   const float_sw4 wgh[4]={17.0/48,59.0/48,43.0/48,49.0/48};   
   const float_sw4 d4a= 2.0/3;
   const float_sw4 d4b=-1.0/12;
   const float_sw4 c6 = 1.0/18;
   const float_sw4 c8 = 1.0/144;
   const float_sw4 al1 = 181507.0/1719312;
   const float_sw4 al2 = -1441.0/39984;
   const float_sw4 al3 = -2593.0/151704;
   const float_sw4 al4 = 11.0/3528;
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = nij*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int base3 = base-nijk;
#define gmu(i,j,k) a_gmu[base +(i)+ni*(j)+nij*(k)]   
#define glambda(i,j,k) a_glambda[base +(i)+ni*(j)+nij*(k)]   
#define u(c,i,j,k)     a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define uacc(c,i,j,k)     a_uacc[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define kap(c,i,j,k) a_kap[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define kapacc(c,i,j,k) a_kapacc[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define bop(q,k) a_bop[q-1+4*(k-1)]

   //#pragma omp parallel
   {
      int kstart = kfirstact;
      int kend   = klastact;
      if( klastact >= nk-3 && onesided[5] == 1 )
         kend = nk-4;
      if( kfirstact <= 4 && onesided[4] == 1 )
      {
         kstart = 5;
	 float_sw4 w8[4]={0,0,1,1};
         float_sw4 w6m[4]={0,0,al1,al1+al2};
	 float_sw4 w6p[4]={0,al1,al1+al2,al1+al2+al3};
	 for( int k=kfirstact; k <= 4; k++ )
#pragma omp parallel for
            for( int j=jfirstact; j <= jlastact; j++ )
#pragma ivdep
               for( int i=ifirstact; i <= ilastact; i++ )
	       {
                  float_sw4 normfact = h3*wgh[k-1];
// Diagonal terms
		  float_sw4 dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
		     d4a*(u(1,i+1,j,k)-u(1,i-1,j,k));
		  float_sw4 dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
		     d4a*(u(2,i,j+1,k)-u(2,i,j-1,k));
		  float_sw4 dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
		     d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k));
		  float_sw4 dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
		     d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k));

		  float_sw4 duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
		     d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k));
		  float_sw4 dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
		     d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k));
		  float_sw4 dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
		     d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k));
		  float_sw4 dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
		     d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k));
		  float_sw4 dwz=0, dmz=0, dwaz=0, dmaz=0;
		  for(int m=1; m<=wb ; m++ )
		  {
		     dwz += bop(k,m)*u(3,i,j,m);
		     dmz += bop(k,m)*kap(3,i,j,m);
		     dwaz+= bop(k,m)*uacc(3,i,j,m);
		     dmaz+= bop(k,m)*kapacc(3,i,j,m);
		  }

		  glambda(i,j,k) = glambda(i,j,k) + (
		     (dux+dvy+dwz)*(dkx+dly+dmz) +
		     dt2o12*( (duax+dvay+dwaz)*(dkx+dly+dmz) +
			      (dux+dvy+dwz)*(dkax+dlay+dmaz) )
						  )*normfact*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + 2*(
		       dux*dkx+dvy*dly+dwz*dmz +
		       dt2o12*( (duax*dkx+dvay*dly+dwaz*dmz)+
				(dux*dkax+dvy*dlay+dwz*dmaz) )
					    )*normfact*ih2;
// Off diagonal stresses
//  xy
		  float_sw4 stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
		     d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
		     d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
		     d4a*(u(1,i,j+1,k)-u(1,i,j-1,k));
		  float_sw4 stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
                      d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
                      d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
		     d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k));
		  float_sw4 stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
                      d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
                      d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
		     d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k));
		  float_sw4 stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
                      d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
                      d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
		     d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k));
		  gmu(i,j,k)= gmu(i,j,k) + (stuxy*stkxy + 
		       dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )*normfact*ih2;

//  xz
		  float_sw4 stuxz = d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
		     d4a*(u(3,i+1,j,k)-u(3,i-1,j,k));
		  float_sw4 stkxz = d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
		     d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k));
		  float_sw4 stuaxz = d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
		     d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k));
		  float_sw4 stkaxz = d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
		     d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k));
		  for(int m=1; m<=wb ; m++ )
		  {
		     stuxz += bop(k,m)*u(1,i,j,m);
		     stkxz += bop(k,m)*kap(1,i,j,m);
		     stuaxz+= bop(k,m)*uacc(1,i,j,m);
		     stkaxz+= bop(k,m)*kapacc(1,i,j,m);
		  }
		  gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
		    dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )*normfact*ih2;
//  yz
		  float_sw4 stuyz = d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
		     d4a*(u(3,i,j+1,k)-u(3,i,j-1,k));
		  float_sw4 stkyz = d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
		     d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k));
		  float_sw4 stuayz = d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
		     d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k));
		  float_sw4 stkayz = d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
		     d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k));
		  for(int m=1; m<=wb ; m++ )
		  {
		     stuyz += bop(k,m)*u(2,i,j,m);
		     stkyz += bop(k,m)*kap(2,i,j,m);
		     stuayz+= bop(k,m)*uacc(2,i,j,m);
		     stkayz+= bop(k,m)*kapacc(2,i,j,m);
		  }
		  gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
			dt2o12*( stuyz*stkayz+stuayz*stkyz) )*normfact*ih2;

// Pos. def extra terms
// x-direction
		  float_sw4 d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
		  3*u(1,i,j,k) -   u(1,i-1,j,k);
		  float_sw4 d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
                  3*u(1,i-1,j,k)-  u(1,i-2,j,k);
		  float_sw4 d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
		  3*kap(1,i,j,k) -   kap(1,i-1,j,k);
		  float_sw4 d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
		  3*kap(1,i-1,j,k) - kap(1,i-2,j,k);
		  float_sw4 d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
		  3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k);
		  float_sw4 d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
		  3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k);
		  float_sw4 d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
		  3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k);
		  float_sw4 d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
		  3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k);
		  float_sw4 pd = (c6*0.5*( d3up*d3kp+d3um*d3km  +
           dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2*wgh[k-1];
		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k)     = gmu(i,j,k) + 2*pd;
		  d3up = u(2,i+2,j,k)-3*u(2,i+1,j,k)+3*u(2,i,j,k)-
		     u(2,i-1,j,k);
		  d3um =u(2,i+1,j,k)-3*u(2,i,j,k)+3*u(2,i-1,j,k)-
		     u(2,i-2,j,k);
		  d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
		     3*kap(2,i,j,k)-kap(2,i-1,j,k);
		  d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
		     3*kap(2,i-1,j,k)-kap(2,i-2,j,k);
		  d3uap = uacc(2,i+2,j,k)-3*uacc(2,i+1,j,k)+
		     3*uacc(2,i,j,k)- uacc(2,i-1,j,k);
		  d3uam =uacc(2,i+1,j,k)-3*uacc(2,i,j,k)+
		     3*uacc(2,i-1,j,k)- uacc(2,i-2,j,k);
		  d3kap = kapacc(2,i+2,j,k)-3*kapacc(2,i+1,j,k)+
		     3*kapacc(2,i,j,k)-kapacc(2,i-1,j,k);
		  d3kam = kapacc(2,i+1,j,k)-3*kapacc(2,i,j,k)+
		     3*kapacc(2,i-1,j,k)-kapacc(2,i-2,j,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) ))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[k-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+3*u(3,i,j,k)-
		                 u(3,i-1,j,k);
		  d3um =u(3,i+1,j,k)-3*u(3,i,j,k)+3*u(3,i-1,j,k)-
		     u(3,i-2,j,k);
		     d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
			3*kap(3,i,j,k)-kap(3,i-1,j,k);
		     d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
			3*kap(3,i-1,j,k)-kap(3,i-2,j,k);
		     d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
			3*uacc(3,i,j,k)- uacc(3,i-1,j,k);
		     d3uam =uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
			3*uacc(3,i-1,j,k)- uacc(3,i-2,j,k);
		     d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
			3*kapacc(3,i,j,k)-kapacc(3,i-1,j,k);
		     d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
			3*kapacc(3,i-1,j,k)-kapacc(3,i-2,j,k);
		     pd = (c6*( 0.5*( d3up*d3kp + d3um*d3km +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[k-1];
		     gmu(i,j,k) = gmu(i,j,k) + pd;

// y-direction
		     d3up = u(1,i,j+2,k)-3*u(1,i,j+1,k)+3*u(1,i,j,k)-
			               u(1,i,j-1,k);
		     d3um =u(1,i,j+1,k)-3*u(1,i,j,k)+3*u(1,i,j-1,k)-
			u(1,i,j-2,k);
		     d3kp = kap(1,i,j+2,k)-3*kap(1,i,j+1,k)+
			3*kap(1,i,j,k)-kap(1,i,j-1,k);
		     d3km = kap(1,i,j+1,k)-3*kap(1,i,j,k)+
			3*kap(1,i,j-1,k)-kap(1,i,j-2,k);
		     d3uap = uacc(1,i,j+2,k)-3*uacc(1,i,j+1,k)+
			3*uacc(1,i,j,k)- uacc(1,i,j-1,k);
		     d3uam =uacc(1,i,j+1,k)-3*uacc(1,i,j,k)+
			3*uacc(1,i,j-1,k)- uacc(1,i,j-2,k);
		     d3kap = kapacc(1,i,j+2,k)-3*kapacc(1,i,j+1,k)+
			3*kapacc(1,i,j,k)-kapacc(1,i,j-1,k);
		     d3kam = kapacc(1,i,j+1,k)-3*kapacc(1,i,j,k)+
			3*kapacc(1,i,j-1,k)-kapacc(1,i,j-2,k);
		     pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[k-1];
		     gmu(i,j,k) = gmu(i,j,k) + pd;

		     d3up = u(2,i,j+2,k)-3*u(2,i,j+1,k)+3*u(2,i,j,k)-
			u(2,i,j-1,k);
		     d3um =u(2,i,j+1,k)-3*u(2,i,j,k)+3*u(2,i,j-1,k)-
			u(2,i,j-2,k);
		     d3kp = kap(2,i,j+2,k)-3*kap(2,i,j+1,k)+
			3*kap(2,i,j,k)-kap(2,i,j-1,k);
		     d3km = kap(2,i,j+1,k)-3*kap(2,i,j,k)+
			3*kap(2,i,j-1,k)-kap(2,i,j-2,k);
		     d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
			3*uacc(2,i,j,k) -   uacc(2,i,j-1,k);
		     d3uam =uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
			3*uacc(2,i,j-1,k)-  uacc(2,i,j-2,k);
		     d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
			3*kapacc(2,i,j,k)-    kapacc(2,i,j-1,k);
		     d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
			3*kapacc(2,i,j-1,k)-  kapacc(2,i,j-2,k);

		     pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[k-1];
		     glambda(i,j,k) = glambda(i,j,k) + pd;
		     gmu(i,j,k) = gmu(i,j,k) + 2*pd;

		     d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+3*u(3,i,j,k)-
			u(3,i,j-1,k);
		     d3um =u(3,i,j+1,k)-3*u(3,i,j,k)+3*u(3,i,j-1,k)-
			u(3,i,j-2,k);
		     d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
			3*kap(3,i,j,k)-kap(3,i,j-1,k);
		     d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
			3*kap(3,i,j-1,k)-kap(3,i,j-2,k);
		     d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
			3*uacc(3,i,j,k)- uacc(3,i,j-1,k);
		     d3uam =uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
			3*uacc(3,i,j-1,k)-uacc(3,i,j-2,k);
		     d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
			3*kapacc(3,i,j,k)-kapacc(3,i,j-1,k);
		     d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
			3*kapacc(3,i,j-1,k)-kapacc(3,i,j-2,k);
		     pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[k-1];
		     gmu(i,j,k) = gmu(i,j,k) + pd;

// z-direction
		     if( k >= 2 )
		     {
			d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
			   u(1,i,j,k-1);
			d3um = u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
			   u(1,i,j,k-2);
			d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
			   3*kap(1,i,j,k)-kap(1,i,j,k-1);
			d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
			   3*kap(1,i,j,k-1)-kap(1,i,j,k-2);
			d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
			   3*uacc(1,i,j,k)-uacc(1,i,j,k-1);
			d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
			   3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2);
			d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
			   3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1);
			d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
			   3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2);
			pd = ( ( 0.5*( w6p[k-1]*d3up*d3kp+w6m[k-1]*d3um*d3km  +
             dt2o12*(w6p[k-1]*d3up*d3kap + w6m[k-1]*d3um*d3kam+ 
                     w6p[k-1]*d3uap*d3kp + w6m[k-1]*d3uam*d3km)))
                 + w8[k-1]*c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2;
			gmu(i,j,k) = gmu(i,j,k) + pd;
			if( k==4 )
			{
			   pd = 0.5*( -al4*d3up*d3kp+
			dt2o12*(-al4*d3up*d3kap -al4*d3uap*d3kp ) )*h3*ih2;
			   gmu(i,j,k+1) = gmu(i,j,k+1) + pd;
			}
			d3up = u(2,i,j,k+2)-3*u(2,i,j,k+1)+
			   3*u(2,i,j,k)- u(2,i,j,k-1);
			d3um =u(2,i,j,k+1)-3*u(2,i,j,k)+
			   3*u(2,i,j,k-1)- u(2,i,j,k-2);
			d3kp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
			   3*kap(2,i,j,k)-kap(2,i,j,k-1);
			d3km = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
			   3*kap(2,i,j,k-1)-kap(2,i,j,k-2);
			d3uap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
			   3*uacc(2,i,j,k)- uacc(2,i,j,k-1);
			d3uam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
			   3*uacc(2,i,j,k-1)- uacc(2,i,j,k-2);
			d3kap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
			   3*kapacc(2,i,j,k)-kapacc(2,i,j,k-1);
			d3kam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
			   3*kapacc(2,i,j,k-1)-kapacc(2,i,j,k-2);
			pd = (( 0.5*( w6p[k-1]*d3up*d3kp+w6m[k-1]*d3um*d3km  +
             dt2o12*(w6p[k-1]*d3up*d3kap + w6m[k-1]*d3um*d3kam+
                     w6p[k-1]*d3uap*d3kp + w6m[k-1]*d3uam*d3km)))
                 + w8[k-1]*c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2;
			gmu(i,j,k) = gmu(i,j,k) + pd;
			if( k == 4 )
			{
			   pd = 0.5*( -al4*d3up*d3kp+
			    dt2o12*(-al4*d3up*d3kap -al4*d3uap*d3kp ) )*h3*ih2;
			   gmu(i,j,k+1) = gmu(i,j,k+1) + pd;
			}
			d3up = u(3,i,j,k+2) - 3*u(3,i,j,k+1) +
			   3*u(3,i,j,k)   -   u(3,i,j,k-1);
			d3um = u(3,i,j,k+1) - 3*u(3,i,j,k) +
			   3*u(3,i,j,k-1) -   u(3,i,j,k-2);
			d3kp = kap(3,i,j,k+2) - 3*kap(3,i,j,k+1) +
			   3*kap(3,i,j,k)   -   kap(3,i,j,k-1);
			d3km = kap(3,i,j,k+1) - 3*kap(3,i,j,k  ) +
			   3*kap(3,i,j,k-1) -   kap(3,i,j,k-2);

			d3uap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
			   3*uacc(3,i,j,k)  -  uacc(3,i,j,k-1);
			d3uam = uacc(3,i,j,k+1)-3*uacc(3,i,j,k) +
			   3*uacc(3,i,j,k-1)-  uacc(3,i,j,k-2);

			d3kap = kapacc(3,i,j,k+2)- 3*kapacc(3,i,j,k+1)+
			   3*kapacc(3,i,j,k)  -   kapacc(3,i,j,k-1);
			d3kam = kapacc(3,i,j,k+1) -3*kapacc(3,i,j,k)+
			   3*kapacc(3,i,j,k-1) -  kapacc(3,i,j,k-2);

			pd = (0.5*( w6p[k-1]*d3up*d3kp  + w6m[k-1]*d3um*d3km  +
		      dt2o12*(    w6p[k-1]*d3up*d3kap + w6m[k-1]*d3um*d3kam + 
                               w6p[k-1]*d3uap*d3kp + w6m[k-1]*d3uam*d3km) )
                 + w8[k-1]*c8*( (d3up-d3um)*(d3kp-d3km) +
                       dt2o12*(  (d3uap-d3uam)*(d3kp-d3km) +
                                 (d3up-d3um)*(d3kap-d3kam)   )))*h3*ih2;
			glambda(i,j,k) = glambda(i,j,k) + pd;
			gmu(i,j,k) = gmu(i,j,k) + 2*pd;
			if( k == 4 )
			{
			   pd = 0.5*( -al4*d3up*d3kp+
			dt2o12*(-al4*d3up*d3kap -al4*d3uap*d3kp ) )*h3*ih2;
			   glambda(i,j,k+1) = glambda(i,j,k+1) + pd;
			   gmu(i,j,k+1)     = gmu(i,j,k+1)   +  2*pd;
			}
		     }
	       }
      }
      
#pragma omp parallel for
	 for( int k=kstart; k <= kend; k++ )
            for( int j=jfirstact; j <= jlastact; j++ )
#pragma ivdep
               for( int i=ifirstact; i <= ilastact; i++ )
	       {
		  float_sw4 normfact = h3;
// Diagonal terms
		  float_sw4 dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
		  d4a*(u(1,i+1,j,k)-u(1,i-1,j,k));
		  float_sw4 dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
		  d4a*(u(2,i,j+1,k)-u(2,i,j-1,k));
		  float_sw4 dwz = d4b*(u(3,i,j,k+2)-u(3,i,j,k-2))+
		  d4a*(u(3,i,j,k+1)-u(3,i,j,k-1));
		  float_sw4 dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
		  d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k));
		  float_sw4 dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
		  d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k));
		  float_sw4 dmz = d4b*(kap(3,i,j,k+2)-kap(3,i,j,k-2))+
		  d4a*(kap(3,i,j,k+1)-kap(3,i,j,k-1));

		  float_sw4 duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
		  d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k));
		  float_sw4 dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
		  d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k));
		  float_sw4 dwaz = d4b*(uacc(3,i,j,k+2)-uacc(3,i,j,k-2))+
		  d4a*(uacc(3,i,j,k+1)-uacc(3,i,j,k-1));
		  float_sw4 dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
		  d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k));
		  float_sw4 dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
		  d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k));
		  float_sw4 dmaz = d4b*(kapacc(3,i,j,k+2)-kapacc(3,i,j,k-2))+
		  d4a*(kapacc(3,i,j,k+1)-kapacc(3,i,j,k-1));

		  glambda(i,j,k) = glambda(i,j,k) + (
                (dux+dvy+dwz)*(dkx+dly+dmz) +
                dt2o12*( (duax+dvay+dwaz)*(dkx+dly+dmz) +
                         (dux+dvy+dwz)*(dkax+dlay+dmaz) )
						     )*normfact*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + 2*(
               dux*dkx+dvy*dly+dwz*dmz +
                dt2o12*( (duax*dkx+dvay*dly+dwaz*dmz)+
                         (dux*dkax+dvy*dlay+dwz*dmaz) )
					       )*normfact*ih2;

// Off diagonal stresses
//  xy
		  float_sw4 stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
                      d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
                      d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
		  d4a*(u(1,i,j+1,k)-u(1,i,j-1,k));
		  float_sw4 stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
                      d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
                      d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
		  d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k));
		  float_sw4 stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
                       d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
                       d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
		  d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k));
		  float_sw4 stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
                       d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
                       d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
		  d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k));

		  gmu(i,j,k)= gmu(i,j,k) + (stuxy*stkxy + 
			   dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )*normfact*ih2;

//  xz
		  float_sw4 stuxz = d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
                      d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
                      d4b*(u(1,i,j,k+2)-u(1,i,j,k-2))+
		  d4a*(u(1,i,j,k+1)-u(1,i,j,k-1));
		  float_sw4 stkxz = d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
                      d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k))+
                      d4b*(kap(1,i,j,k+2)-kap(1,i,j,k-2))+
		  d4a*(kap(1,i,j,k+1)-kap(1,i,j,k-1));
		  float_sw4 stuaxz = d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
                       d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k))+
                       d4b*(uacc(1,i,j,k+2)-uacc(1,i,j,k-2))+
		  d4a*(uacc(1,i,j,k+1)-uacc(1,i,j,k-1));
		  float_sw4 stkaxz = d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
                       d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k))+
                       d4b*(kapacc(1,i,j,k+2)-kapacc(1,i,j,k-2))+
		  d4a*(kapacc(1,i,j,k+1)-kapacc(1,i,j,k-1));
		  gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
			  dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )*normfact*ih2;

//  yz
		  float_sw4 stuyz = d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
                      d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
                      d4b*(u(2,i,j,k+2)-u(2,i,j,k-2))+
		  d4a*(u(2,i,j,k+1)-u(2,i,j,k-1));
		  float_sw4 stkyz = d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
                      d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k))+
                      d4b*(kap(2,i,j,k+2)-kap(2,i,j,k-2))+
		  d4a*(kap(2,i,j,k+1)-kap(2,i,j,k-1));
		  float_sw4 stuayz = d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
                       d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k))+
                       d4b*(uacc(2,i,j,k+2)-uacc(2,i,j,k-2))+
		  d4a*(uacc(2,i,j,k+1)-uacc(2,i,j,k-1));
		  float_sw4 stkayz = d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
                       d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k))+
                       d4b*(kapacc(2,i,j,k+2)-kapacc(2,i,j,k-2))+
		  d4a*(kapacc(2,i,j,k+1)-kapacc(2,i,j,k-1));
		  gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
	           dt2o12*( stuyz*stkayz+stuayz*stkyz) )*normfact*ih2;
// Pos. def extra terms
// x-direction
		  float_sw4 d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
		  3*u(1,i,j,k) -   u(1,i-1,j,k);
		  float_sw4 d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
                  3*u(1,i-1,j,k)-  u(1,i-2,j,k);
		  float_sw4 d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
		  3*kap(1,i,j,k) -   kap(1,i-1,j,k);
		  float_sw4 d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
		  3*kap(1,i-1,j,k) - kap(1,i-2,j,k);
		  float_sw4 d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
		  3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k);
		  float_sw4 d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
		  3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k);
		  float_sw4 d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
		  3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k);
		  float_sw4 d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
		  3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k);
		  float_sw4 pd = (c6*0.5*( d3up*d3kp+d3um*d3km  +
           dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2;
		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k)     = gmu(i,j,k) + 2*pd;

		  d3up = u(2,i+2,j,k)- 3*u(2,i+1,j,k) +
		  3*u(2,i,j,k)  -   u(2,i-1,j,k);
		  d3um = u(2,i+1,j,k)- 3*u(2,i,j,k) +
		  3*u(2,i-1,j,k)-   u(2,i-2,j,k);
		  d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
		  3*kap(2,i,j,k)  -  kap(2,i-1,j,k);
		  d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
		  3*kap(2,i-1,j,k)-  kap(2,i-2,j,k);
		  d3uap = uacc(2,i+2,j,k) - 3*uacc(2,i+1,j,k)+
		  3*uacc(2,i,j,k)   -   uacc(2,i-1,j,k);
		  d3uam = uacc(2,i+1,j,k) - 3*uacc(2,i,j,k)+
		  3*uacc(2,i-1,j,k) -   uacc(2,i-2,j,k);
		  d3kap = kapacc(2,i+2,j,k) - 3*kapacc(2,i+1,j,k)+
		  3*kapacc(2,i,j,k)   -   kapacc(2,i-1,j,k);
		  d3kam = kapacc(2,i+1,j,k) - 3*kapacc(2,i,j,k)+
		  3*kapacc(2,i-1,j,k) -   kapacc(2,i-2,j,k);
		  pd = (c6*0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+
		  3*u(3,i,j,k)  -  u(3,i-1,j,k);
		  d3um = u(3,i+1,j,k)-3*u(3,i,j,k)+
		  3*u(3,i-1,j,k)-  u(3,i-2,j,k);
		  d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
		  3*kap(3,i,j,k)  -  kap(3,i-1,j,k);
		  d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
		  3*kap(3,i-1,j,k)-  kap(3,i-2,j,k);
		  d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
		  3*uacc(3,i,j,k) -   uacc(3,i-1,j,k);
		  d3uam = uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
		  3*uacc(3,i-1,j,k)-  uacc(3,i-2,j,k);
		  d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
		  3*kapacc(3,i,j,k) -   kapacc(3,i-1,j,k);
		  d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
		  3*kapacc(3,i-1,j,k)-  kapacc(3,i-2,j,k);

		  pd = ( c6*0.5*( d3up*d3kp + d3um*d3km +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + pd;

// y-direction
		  d3up = u(1,i,j+2,k) - 3*u(1,i,j+1,k) +
		  3*u(1,i,j,k)   -   u(1,i,j-1,k);
		  d3um = u(1,i,j+1,k) - 3*u(1,i,j,k) +
		  3*u(1,i,j-1,k) -   u(1,i,j-2,k);
		  d3kp = kap(1,i,j+2,k) - 3*kap(1,i,j+1,k) +
		  3*kap(1,i,j,k)   -   kap(1,i,j-1,k);
		  d3km = kap(1,i,j+1,k) - 3*kap(1,i,j,k)+
		  3*kap(1,i,j-1,k) -   kap(1,i,j-2,k);
		  d3uap = uacc(1,i,j+2,k)- 3*uacc(1,i,j+1,k)+
		  3*uacc(1,i,j,k) -    uacc(1,i,j-1,k);
		  d3uam =uacc(1,i,j+1,k) - 3*uacc(1,i,j,k)+
		  3*uacc(1,i,j-1,k) -   uacc(1,i,j-2,k);
		  d3kap = kapacc(1,i,j+2,k) - 3*kapacc(1,i,j+1,k)+
		  3*kapacc(1,i,j,k)   -   kapacc(1,i,j-1,k);
		  d3kam = kapacc(1,i,j+1,k) - 3*kapacc(1,i,j,k)+
		  3*kapacc(1,i,j-1,k) -   kapacc(1,i,j-2,k);

		  pd = ( c6*0.5*( d3up*d3kp + d3um*d3km  +
           dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam) ) ) )*h3*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(2,i,j+2,k) - 3*u(2,i,j+1,k)+
		  3*u(2,i,j,k)   -   u(2,i,j-1,k);
		  d3um =u(2,i,j+1,k)  - 3*u(2,i,j,k)+
                  3*u(2,i,j-1,k) -    u(2,i,j-2,k);
		  d3kp = kap(2,i,j+2,k) - 3*kap(2,i,j+1,k)+
		  3*kap(2,i,j,k)   -   kap(2,i,j-1,k);
		  d3km = kap(2,i,j+1,k) - 3*kap(2,i,j,k)+
		  3*kap(2,i,j-1,k) -   kap(2,i,j-2,k);
		  d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
		  3*uacc(2,i,j,k)  -  uacc(2,i,j-1,k);
		  d3uam = uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
		  3*uacc(2,i,j-1,k)-  uacc(2,i,j-2,k);
		  d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
		  3*kapacc(2,i,j,k)  -  kapacc(2,i,j-1,k);
		  d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
		  3*kapacc(2,i,j-1,k)-  kapacc(2,i,j-2,k);

		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2;

		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k)     = gmu(i,j,k) + 2*pd;

		  d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+
		  3*u(3,i,j,k)  -  u(3,i,j-1,k);
		  d3um = u(3,i,j+1,k)-3*u(3,i,j,k)+
		  3*u(3,i,j-1,k)-  u(3,i,j-2,k);
		  d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
		  3*kap(3,i,j,k)  -  kap(3,i,j-1,k);
		  d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
		  3*kap(3,i,j-1,k)-  kap(3,i,j-2,k);
		  d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
		  3*uacc(3,i,j,k)  -  uacc(3,i,j-1,k);
		  d3uam = uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
		  3*uacc(3,i,j-1,k)-  uacc(3,i,j-2,k);
		  d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
		  3*kapacc(3,i,j,k)  -  kapacc(3,i,j-1,k);
		  d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
		  3*kapacc(3,i,j-1,k)-  kapacc(3,i,j-2,k);
		  pd = ( c6*0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam) ) ) )*h3*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + pd;

// z-direction
		  d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
		  u(1,i,j,k-1);
		  d3um =u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
		  u(1,i,j,k-2);
		  d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
		  3*kap(1,i,j,k)-kap(1,i,j,k-1);
		  d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
		  3*kap(1,i,j,k-1)-kap(1,i,j,k-2);
		  d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
		  3*uacc(1,i,j,k)-uacc(1,i,j,k-1);
		  d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
		  3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2);
		  d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
		  3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1);
		  d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
		  3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2);
		  pd = ( c6*0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(2,i,j,k+2)-3*u(2,i,j,k+1)+
		  3*u(2,i,j,k)  -  u(2,i,j,k-1);
		  d3um =u(2,i,j,k+1)-3*u(2,i,j,k)+
                  3*u(2,i,j,k-1)-  u(2,i,j,k-2);
		  d3kp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
		  3*kap(2,i,j,k)  -  kap(2,i,j,k-1);
		  d3km = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
		  3*kap(2,i,j,k-1)-  kap(2,i,j,k-2);
		  d3uap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
		  3*uacc(2,i,j,k)  -  uacc(2,i,j,k-1);
		  d3uam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
		  3*uacc(2,i,j,k-1)-  uacc(2,i,j,k-2);
		  d3kap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
		  3*kapacc(2,i,j,k)  -  kapacc(2,i,j,k-1);
		  d3kam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
		  3*kapacc(2,i,j,k-1) - kapacc(2,i,j,k-2);

		  pd = ( c6*0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(3,i,j,k+2)-3*u(3,i,j,k+1)+3*u(3,i,j,k)-
		  u(3,i,j,k-1);
		  d3um =u(3,i,j,k+1)-3*u(3,i,j,k)+3*u(3,i,j,k-1)-
		  u(3,i,j,k-2);
		  d3kp = kap(3,i,j,k+2)-3*kap(3,i,j,k+1)+
		  3*kap(3,i,j,k)-kap(3,i,j,k-1);
		  d3km = kap(3,i,j,k+1)-3*kap(3,i,j,k)+
		  3*kap(3,i,j,k-1)-kap(3,i,j,k-2);
		  d3uap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
		  3*uacc(3,i,j,k)-uacc(3,i,j,k-1);
		  d3uam =uacc(3,i,j,k+1)-3*uacc(3,i,j,k)+
		  3*uacc(3,i,j,k-1)- uacc(3,i,j,k-2);
		  d3kap = kapacc(3,i,j,k+2)-3*kapacc(3,i,j,k+1)+
		  3*kapacc(3,i,j,k)-kapacc(3,i,j,k-1);
		  d3kam = kapacc(3,i,j,k+1)-3*kapacc(3,i,j,k)+
		  3*kapacc(3,i,j,k-1)-kapacc(3,i,j,k-2);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2;
		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k)     = gmu(i,j,k)    + 2*pd;
	       }
	 }
      if( klastact >= nk-3 && onesided[5] == 1 )
      {
	 float_sw4 w8[4]={0,0,1,1};
         float_sw4 w6p[4]={0,0,al1,al1+al2};
	 float_sw4 w6m[4]={0,al1,al1+al2,al1+al2+al3};
	 for( int k=nk-3; k <= klastact; k++ )
#pragma omp parallel for
            for( int j=jfirstact; j <= jlastact; j++ )
#pragma ivdep
               for( int i=ifirstact; i <= ilastact; i++ )
	       {
                  int kk=nk-k+1;
                  float_sw4 normfact = h3*wgh[kk-1];
// Diagonal terms
		  float_sw4 dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
		     d4a*(u(1,i+1,j,k)-u(1,i-1,j,k));
		  float_sw4 dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
		     d4a*(u(2,i,j+1,k)-u(2,i,j-1,k));
		  float_sw4 dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
		     d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k));
		  float_sw4 dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
		     d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k));

		  float_sw4 duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
		     d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k));
		  float_sw4 dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
		     d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k));
		  float_sw4 dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
		     d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k));
		  float_sw4 dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
		     d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k));
		  float_sw4 dwz=0, dmz=0, dwaz=0, dmaz=0;
		  for(int m=1; m<=wb ; m++ )
		  {
		     dwz -= bop(kk,m)*u(3,i,j,nk-m+1);
		     dmz -= bop(kk,m)*kap(3,i,j,nk-m+1);
		     dwaz-= bop(kk,m)*uacc(3,i,j,nk-m+1);
		     dmaz-= bop(kk,m)*kapacc(3,i,j,nk-m+1);
		  }

		  glambda(i,j,k) = glambda(i,j,k) + (
		     (dux+dvy+dwz)*(dkx+dly+dmz) +
		     dt2o12*( (duax+dvay+dwaz)*(dkx+dly+dmz) +
			      (dux+dvy+dwz)*(dkax+dlay+dmaz) )
						  )*normfact*ih2;
		  gmu(i,j,k) = gmu(i,j,k) + 2*(
		       dux*dkx+dvy*dly+dwz*dmz +
		       dt2o12*( (duax*dkx+dvay*dly+dwaz*dmz)+
				(dux*dkax+dvy*dlay+dwz*dmaz) )
					    )*normfact*ih2;
// Off diagonal stresses
//  xy
		  float_sw4 stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
		     d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
		     d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
		     d4a*(u(1,i,j+1,k)-u(1,i,j-1,k));
		  float_sw4 stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
                      d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
                      d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
		     d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k));
		  float_sw4 stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
                      d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
                      d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
		     d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k));
		  float_sw4 stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
                      d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
                      d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
		     d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k));
		  gmu(i,j,k)= gmu(i,j,k) + (stuxy*stkxy + 
		       dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )*normfact*ih2;

//  xz
		  float_sw4 stuxz = d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
		     d4a*(u(3,i+1,j,k)-u(3,i-1,j,k));
		  float_sw4 stkxz = d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
		     d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k));
		  float_sw4 stuaxz = d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
		     d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k));
		  float_sw4 stkaxz = d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
		     d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k));
		  for(int m=1; m<=wb ; m++ )
		  {
		     stuxz -= bop(kk,m)*u(1,i,j,nk-m+1);
		     stkxz -= bop(kk,m)*kap(1,i,j,nk-m+1);
		     stuaxz-= bop(kk,m)*uacc(1,i,j,nk-m+1);
		     stkaxz-= bop(kk,m)*kapacc(1,i,j,nk-m+1);
		  }
		  gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
		    dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )*normfact*ih2;
//  yz
		  float_sw4 stuyz = d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
		     d4a*(u(3,i,j+1,k)-u(3,i,j-1,k));
		  float_sw4 stkyz = d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
		     d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k));
		  float_sw4 stuayz = d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
		     d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k));
		  float_sw4 stkayz = d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
		     d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k));
		  for(int m=1; m<=wb ; m++ )
		  {
		     stuyz -= bop(kk,m)*u(2,i,j,nk-m+1);
		     stkyz -= bop(kk,m)*kap(2,i,j,nk-m+1);
		     stuayz-= bop(kk,m)*uacc(2,i,j,nk-m+1);
		     stkayz-= bop(kk,m)*kapacc(2,i,j,nk-m+1);
		  }
		  gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
			dt2o12*( stuyz*stkayz+stuayz*stkyz) )*normfact*ih2;

// Pos. def extra terms
// x-direction
		  float_sw4 d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
		  3*u(1,i,j,k) -   u(1,i-1,j,k);
		  float_sw4 d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
                  3*u(1,i-1,j,k)-  u(1,i-2,j,k);
		  float_sw4 d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
		  3*kap(1,i,j,k) -   kap(1,i-1,j,k);
		  float_sw4 d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
		  3*kap(1,i-1,j,k) - kap(1,i-2,j,k);
		  float_sw4 d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
		  3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k);
		  float_sw4 d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
		  3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k);
		  float_sw4 d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
		  3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k);
		  float_sw4 d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
		  3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k);
		  float_sw4 pd = (c6*0.5*( d3up*d3kp+d3um*d3km  +
           dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2*wgh[kk-1];
		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k)     = gmu(i,j,k) + 2*pd;
		  d3up = u(2,i+2,j,k)-3*u(2,i+1,j,k)+3*u(2,i,j,k)-
		     u(2,i-1,j,k);
		  d3um =u(2,i+1,j,k)-3*u(2,i,j,k)+3*u(2,i-1,j,k)-
		     u(2,i-2,j,k);
		  d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
		     3*kap(2,i,j,k)-kap(2,i-1,j,k);
		  d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
		     3*kap(2,i-1,j,k)-kap(2,i-2,j,k);
		  d3uap = uacc(2,i+2,j,k)-3*uacc(2,i+1,j,k)+
		     3*uacc(2,i,j,k)- uacc(2,i-1,j,k);
		  d3uam =uacc(2,i+1,j,k)-3*uacc(2,i,j,k)+
		     3*uacc(2,i-1,j,k)- uacc(2,i-2,j,k);
		  d3kap = kapacc(2,i+2,j,k)-3*kapacc(2,i+1,j,k)+
		     3*kapacc(2,i,j,k)-kapacc(2,i-1,j,k);
		  d3kam = kapacc(2,i+1,j,k)-3*kapacc(2,i,j,k)+
		     3*kapacc(2,i-1,j,k)-kapacc(2,i-2,j,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) ))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[kk-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+3*u(3,i,j,k)-
		                 u(3,i-1,j,k);
		  d3um =u(3,i+1,j,k)-3*u(3,i,j,k)+3*u(3,i-1,j,k)-
		     u(3,i-2,j,k);
		     d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
			3*kap(3,i,j,k)-kap(3,i-1,j,k);
		     d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
			3*kap(3,i-1,j,k)-kap(3,i-2,j,k);
		     d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
			3*uacc(3,i,j,k)- uacc(3,i-1,j,k);
		     d3uam =uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
			3*uacc(3,i-1,j,k)- uacc(3,i-2,j,k);
		     d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
			3*kapacc(3,i,j,k)-kapacc(3,i-1,j,k);
		     d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
			3*kapacc(3,i-1,j,k)-kapacc(3,i-2,j,k);
		     pd = (c6*( 0.5*( d3up*d3kp + d3um*d3km +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[kk-1];
		     gmu(i,j,k) = gmu(i,j,k) + pd;

// y-direction
		     d3up = u(1,i,j+2,k)-3*u(1,i,j+1,k)+3*u(1,i,j,k)-
			               u(1,i,j-1,k);
		     d3um =u(1,i,j+1,k)-3*u(1,i,j,k)+3*u(1,i,j-1,k)-
			u(1,i,j-2,k);
		     d3kp = kap(1,i,j+2,k)-3*kap(1,i,j+1,k)+
			3*kap(1,i,j,k)-kap(1,i,j-1,k);
		     d3km = kap(1,i,j+1,k)-3*kap(1,i,j,k)+
			3*kap(1,i,j-1,k)-kap(1,i,j-2,k);
		     d3uap = uacc(1,i,j+2,k)-3*uacc(1,i,j+1,k)+
			3*uacc(1,i,j,k)- uacc(1,i,j-1,k);
		     d3uam =uacc(1,i,j+1,k)-3*uacc(1,i,j,k)+
			3*uacc(1,i,j-1,k)- uacc(1,i,j-2,k);
		     d3kap = kapacc(1,i,j+2,k)-3*kapacc(1,i,j+1,k)+
			3*kapacc(1,i,j,k)-kapacc(1,i,j-1,k);
		     d3kam = kapacc(1,i,j+1,k)-3*kapacc(1,i,j,k)+
			3*kapacc(1,i,j-1,k)-kapacc(1,i,j-2,k);
		     pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[kk-1];
		     gmu(i,j,k) = gmu(i,j,k) + pd;

		     d3up = u(2,i,j+2,k)-3*u(2,i,j+1,k)+3*u(2,i,j,k)-
			u(2,i,j-1,k);
		     d3um =u(2,i,j+1,k)-3*u(2,i,j,k)+3*u(2,i,j-1,k)-
			u(2,i,j-2,k);
		     d3kp = kap(2,i,j+2,k)-3*kap(2,i,j+1,k)+
			3*kap(2,i,j,k)-kap(2,i,j-1,k);
		     d3km = kap(2,i,j+1,k)-3*kap(2,i,j,k)+
			3*kap(2,i,j-1,k)-kap(2,i,j-2,k);
		     d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
			3*uacc(2,i,j,k) -   uacc(2,i,j-1,k);
		     d3uam =uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
			3*uacc(2,i,j-1,k)-  uacc(2,i,j-2,k);
		     d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
			3*kapacc(2,i,j,k)-    kapacc(2,i,j-1,k);
		     d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
			3*kapacc(2,i,j-1,k)-  kapacc(2,i,j-2,k);

		     pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[kk-1];
		     glambda(i,j,k) = glambda(i,j,k) + pd;
		     gmu(i,j,k) = gmu(i,j,k) + 2*pd;

		     d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+3*u(3,i,j,k)-
			u(3,i,j-1,k);
		     d3um =u(3,i,j+1,k)-3*u(3,i,j,k)+3*u(3,i,j-1,k)-
			u(3,i,j-2,k);
		     d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
			3*kap(3,i,j,k)-kap(3,i,j-1,k);
		     d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
			3*kap(3,i,j-1,k)-kap(3,i,j-2,k);
		     d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
			3*uacc(3,i,j,k)- uacc(3,i,j-1,k);
		     d3uam =uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
			3*uacc(3,i,j-1,k)-uacc(3,i,j-2,k);
		     d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
			3*kapacc(3,i,j,k)-kapacc(3,i,j-1,k);
		     d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
			3*kapacc(3,i,j-1,k)-kapacc(3,i,j-2,k);
		     pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh[kk-1];
		     gmu(i,j,k) = gmu(i,j,k) + pd;

// z-direction
		     if( k <= nk-1 )
		     {
			d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
			   u(1,i,j,k-1);
			d3um = u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
			   u(1,i,j,k-2);
			d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
			   3*kap(1,i,j,k)-kap(1,i,j,k-1);
			d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
			   3*kap(1,i,j,k-1)-kap(1,i,j,k-2);
			d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
			   3*uacc(1,i,j,k)-uacc(1,i,j,k-1);
			d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
			   3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2);
			d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
			   3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1);
			d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
			   3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2);
			pd = ( ( 0.5*( w6p[kk-1]*d3up*d3kp+w6m[kk-1]*d3um*d3km  +
             dt2o12*(w6p[kk-1]*d3up*d3kap + w6m[kk-1]*d3um*d3kam+ 
                     w6p[kk-1]*d3uap*d3kp + w6m[kk-1]*d3uam*d3km)))
                 + w8[kk-1]*c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2;
			gmu(i,j,k) = gmu(i,j,k) + pd;
			if( k==nk-3 )
			{
			   pd = 0.5*( -al4*d3um*d3km+
			dt2o12*(-al4*d3um*d3kam -al4*d3uam*d3km ) )*h3*ih2;
			   gmu(i,j,k-1) = gmu(i,j,k-1) + pd;
			}
			d3up = u(2,i,j,k+2)-3*u(2,i,j,k+1)+
			   3*u(2,i,j,k)- u(2,i,j,k-1);
			d3um =u(2,i,j,k+1)-3*u(2,i,j,k)+
			   3*u(2,i,j,k-1)- u(2,i,j,k-2);
			d3kp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
			   3*kap(2,i,j,k)-kap(2,i,j,k-1);
			d3km = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
			   3*kap(2,i,j,k-1)-kap(2,i,j,k-2);
			d3uap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
			   3*uacc(2,i,j,k)- uacc(2,i,j,k-1);
			d3uam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
			   3*uacc(2,i,j,k-1)- uacc(2,i,j,k-2);
			d3kap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
			   3*kapacc(2,i,j,k)-kapacc(2,i,j,k-1);
			d3kam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
			   3*kapacc(2,i,j,k-1)-kapacc(2,i,j,k-2);
			pd = (( 0.5*( w6p[kk-1]*d3up*d3kp+w6m[kk-1]*d3um*d3km  +
             dt2o12*(w6p[kk-1]*d3up*d3kap + w6m[kk-1]*d3um*d3kam+
                     w6p[kk-1]*d3uap*d3kp + w6m[kk-1]*d3uam*d3km)))
                 + w8[kk-1]*c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2;
			gmu(i,j,k) = gmu(i,j,k) + pd;
			if( k == nk-3 )
			{
			   pd = 0.5*( -al4*d3um*d3km+
			    dt2o12*(-al4*d3um*d3kam -al4*d3uam*d3km ) )*h3*ih2;
			   gmu(i,j,k-1) = gmu(i,j,k-1) + pd;
			}
			d3up = u(3,i,j,k+2) - 3*u(3,i,j,k+1) +
			   3*u(3,i,j,k)   -   u(3,i,j,k-1);
			d3um = u(3,i,j,k+1) - 3*u(3,i,j,k) +
			   3*u(3,i,j,k-1) -   u(3,i,j,k-2);
			d3kp = kap(3,i,j,k+2) - 3*kap(3,i,j,k+1) +
			   3*kap(3,i,j,k)   -   kap(3,i,j,k-1);
			d3km = kap(3,i,j,k+1) - 3*kap(3,i,j,k  ) +
			   3*kap(3,i,j,k-1) -   kap(3,i,j,k-2);

			d3uap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
			   3*uacc(3,i,j,k)  -  uacc(3,i,j,k-1);
			d3uam = uacc(3,i,j,k+1)-3*uacc(3,i,j,k) +
			   3*uacc(3,i,j,k-1)-  uacc(3,i,j,k-2);

			d3kap = kapacc(3,i,j,k+2)- 3*kapacc(3,i,j,k+1)+
			   3*kapacc(3,i,j,k)  -   kapacc(3,i,j,k-1);
			d3kam = kapacc(3,i,j,k+1) -3*kapacc(3,i,j,k)+
			   3*kapacc(3,i,j,k-1) -  kapacc(3,i,j,k-2);

			pd = (0.5*( w6p[kk-1]*d3up*d3kp  + w6m[kk-1]*d3um*d3km  +
		      dt2o12*(    w6p[kk-1]*d3up*d3kap + w6m[kk-1]*d3um*d3kam + 
                               w6p[kk-1]*d3uap*d3kp + w6m[kk-1]*d3uam*d3km) )
                 + w8[kk-1]*c8*( (d3up-d3um)*(d3kp-d3km) +
                       dt2o12*(  (d3uap-d3uam)*(d3kp-d3km) +
                                 (d3up-d3um)*(d3kap-d3kam)   )))*h3*ih2;
			glambda(i,j,k) = glambda(i,j,k) + pd;
			gmu(i,j,k) = gmu(i,j,k) + 2*pd;
			if( k == nk-3 )
			{
			   pd = 0.5*( -al4*d3um*d3km+
			dt2o12*(-al4*d3um*d3kam -al4*d3uam*d3km ) )*h3*ih2;
			   glambda(i,j,k-1) = glambda(i,j,k-1) + pd;
			   gmu(i,j,k-1)     = gmu(i,j,k-1)   +  2*pd;
			}
		     }
	       }
      }
}

//-----------------------------------------------------------------------
void addgradmulac_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		     int ifirstact, int ilastact, int jfirstact, int jlastact, 
		     int kfirstact, int klastact, int nk,
		     float_sw4* __restrict__ a_kap, float_sw4* __restrict__ a_kapacc, 
		     float_sw4* __restrict__ a_u,   float_sw4* __restrict__ a_uacc,
		     float_sw4* __restrict__ a_gmu, float_sw4* __restrict__ a_glambda,
		     float_sw4 dt, float_sw4 h, float_sw4* __restrict__ a_met,
		     float_sw4* __restrict__ a_jac, int onesided[6],
		     int nb, int wb, float_sw4* __restrict__ a_bop )
{
   const float_sw4 dt2o12 = dt*dt/12;
   const float_sw4 wgh[4]={17.0/48,59.0/48,43.0/48,49.0/48};   
   const float_sw4 d4a= 2.0/3;
   const float_sw4 d4b=-1.0/12;
   const float_sw4 c6 = 1.0/18;
   const float_sw4 c8 = 1.0/144;
   const float_sw4 al1 = 181507.0/1719312;
   const float_sw4 al2 = -1441.0/39984;
   const float_sw4 al3 = -2593.0/151704;
   const float_sw4 al4 = 11.0/3528;
   const int ni    = ilast-ifirst+1;
   const int nij   = ni*(jlast-jfirst+1);
   const int nijk  = nij*(klast-kfirst+1);
   const int base  = -(ifirst+ni*jfirst+nij*kfirst);
   const int base3 = base-nijk;
#define jac(i,j,k) a_jac[base +(i)+ni*(j)+nij*(k)]   
#define gmu(i,j,k) a_gmu[base +(i)+ni*(j)+nij*(k)]   
#define glambda(i,j,k) a_glambda[base +(i)+ni*(j)+nij*(k)]   
#define u(c,i,j,k)     a_u[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define uacc(c,i,j,k)     a_uacc[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define kap(c,i,j,k) a_kap[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define kapacc(c,i,j,k) a_kapacc[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define met(c,i,j,k) a_met[base3+(i)+ni*(j)+nij*(k)+nijk*(c)]   
#define bop(q,k) a_bop[q-1+4*(k-1)]

   int kstart = kfirstact;
   int kend   = klastact;
   if( klastact >= nk-3 && onesided[5] == 1 )
      kend = nk-4;
   //#pragma omp parallel
	 {
      if( kfirstact <= 4 && onesided[4] == 1 )
      {
         kstart = 5;
	 float_sw4 w8[4] ={0,0,1,1};
         float_sw4 w6m[4]={0,0,al1,al1+al2};
	 float_sw4 w6p[4]={0,al1,al1+al2,al1+al2+al3};
         //#pragma omp for
	 for( int k=kfirstact; k <= 4; k++ )
            for( int j=jfirstact; j <= jlastact; j++ )
               //#pragma ivdep
               for( int i=ifirstact; i <= ilastact; i++ )
	       {
		  float_sw4 normfact = wgh[k-1];
// Diagonal terms
		  float_sw4 dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
		     d4a*(u(1,i+1,j,k)-u(1,i-1,j,k));
		  float_sw4 dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
		     d4a*(u(2,i,j+1,k)-u(2,i,j-1,k));
		  float_sw4 dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
		     d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k));
		  float_sw4 dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
		     d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k));

		  float_sw4 duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
		     d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k));
		  float_sw4 dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
		     d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k));
		  float_sw4 dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
		     d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k));
		  float_sw4 dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
		     d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k));
		  float_sw4 duz = 0, dvz = 0, dwz = 0, dkz = 0, dlz = 0, dmz = 0;
		  float_sw4 duaz=0, dvaz=0, dwaz=0, dkaz=0, dlaz=0, dmaz=0;
		  for( int m=1; m <= wb ;m++ )
		  {
		     duz += bop(k,m)*u(1,i,j,m);
		     dvz += bop(k,m)*u(2,i,j,m);
		     dwz += bop(k,m)*u(3,i,j,m);

		     dkz += bop(k,m)*kap(1,i,j,m);
		     dlz += bop(k,m)*kap(2,i,j,m);
		     dmz += bop(k,m)*kap(3,i,j,m);

		     duaz += bop(k,m)*uacc(1,i,j,m);
		     dvaz += bop(k,m)*uacc(2,i,j,m);
		     dwaz += bop(k,m)*uacc(3,i,j,m);

		     dkaz += bop(k,m)*kapacc(1,i,j,m);
		     dlaz += bop(k,m)*kapacc(2,i,j,m);
		     dmaz += bop(k,m)*kapacc(3,i,j,m);
		  }
		  dux = met(1,i,j,k)*dux + met(2,i,j,k)*duz;
		  dvy = met(1,i,j,k)*dvy + met(3,i,j,k)*dvz;
		  dkx = met(1,i,j,k)*dkx + met(2,i,j,k)*dkz;
		  dly = met(1,i,j,k)*dly + met(3,i,j,k)*dlz;
		  duax = met(1,i,j,k)*duax + met(2,i,j,k)*duaz;
		  dvay = met(1,i,j,k)*dvay + met(3,i,j,k)*dvaz;
		  dkax = met(1,i,j,k)*dkax + met(2,i,j,k)*dkaz;
		  dlay = met(1,i,j,k)*dlay + met(3,i,j,k)*dlaz;

		  float_sw4 dwzm  = dwz*met(4,i,j,k);
		  float_sw4 dmzm  = dmz*met(4,i,j,k);
		  float_sw4 dwazm = dwaz*met(4,i,j,k);
		  float_sw4 dmazm = dmaz*met(4,i,j,k);

		  glambda(i,j,k) = glambda(i,j,k) + (
                (dux+dvy+dwzm)*(dkx+dly+dmzm) +
                dt2o12*( (duax+dvay+dwazm)*(dkx+dly+dmzm) +
                         (dux+dvy+dwzm)*(dkax+dlay+dmazm) )
						     )*normfact;
		  gmu(i,j,k) = gmu(i,j,k) + 2*(
               dux*dkx+dvy*dly+dwzm*dmzm +
                dt2o12*( (duax*dkx+dvay*dly+dwazm*dmzm)+
                         (dux*dkax+dvy*dlay+dwzm*dmazm) )
					       )*normfact;
// Off diagonal stresses
//  xy
		  float_sw4 stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
                      d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
                      d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
		     d4a*(u(1,i,j+1,k)-u(1,i,j-1,k));
		  float_sw4 stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
                      d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
                      d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
		     d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k));
		  float_sw4 stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
                      d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
                      d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
		     d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k));
		  float_sw4 stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
                      d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
                      d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
		     d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k));

		  stuxy = met(1,i,j,k)*stuxy+met(2,i,j,k)*dvz +
		     met(3,i,j,k)*duz;
		  stkxy = met(1,i,j,k)*stkxy+met(2,i,j,k)*dlz +
		     met(3,i,j,k)*dkz;
		  stuaxy= met(1,i,j,k)*stuaxy+met(2,i,j,k)*dvaz +
		     met(3,i,j,k)*duaz;
		  stkaxy= met(1,i,j,k)*stkaxy+met(2,i,j,k)*dlaz +
		     met(3,i,j,k)*dkaz;

		  gmu(i,j,k)= gmu(i,j,k) + (
           stuxy*stkxy + dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )*
		     normfact;
//  xz
		  float_sw4 stuxz = met(1,i,j,k)*(
		      d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
		      d4a*(u(3,i+1,j,k)-u(3,i-1,j,k)) );
		  float_sw4 stkxz = met(1,i,j,k)*(
		      d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
                      d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k)) );
		  float_sw4 stuaxz = met(1,i,j,k)*(
                      d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
                      d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k)) );
		  float_sw4 stkaxz = met(1,i,j,k)*(
                      d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
                      d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k)) );

		  stuxz = stuxz + met(2,i,j,k)*dwz+met(4,i,j,k)*duz;
		  stkxz = stkxz + met(2,i,j,k)*dmz+met(4,i,j,k)*dkz;
		  stuaxz= stuaxz+ met(2,i,j,k)*dwaz+met(4,i,j,k)*duaz;
		  stkaxz= stkaxz+ met(2,i,j,k)*dmaz+met(4,i,j,k)*dkaz;

		  gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
		       dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )*normfact;
//  yz
		  float_sw4 stuyz = met(1,i,j,k)*(d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
						  d4a*(u(3,i,j+1,k)-u(3,i,j-1,k)));
		  float_sw4 stkyz = met(1,i,j,k)*(
                      d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
                      d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k)) );
		  float_sw4 stuayz = met(1,i,j,k)*(
                      d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
                      d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k)) );
		  float_sw4 stkayz = met(1,i,j,k)*(
                      d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
                      d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k)) );

		  stuyz = stuyz + met(3,i,j,k)*dwz+met(4,i,j,k)*dvz;
		  stkyz = stkyz + met(3,i,j,k)*dmz+met(4,i,j,k)*dlz;
		  stuayz= stuayz+ met(3,i,j,k)*dwaz+met(4,i,j,k)*dvaz;
		  stkayz= stkayz+ met(3,i,j,k)*dmaz+met(4,i,j,k)*dlaz;

		  gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
 	               dt2o12*( stuyz*stkayz+stuayz*stkyz) )*normfact;
// Pos. def extra terms
// x-direction
//   metric is constant for the x and y directions, just use as a norm weight
		  float_sw4 m1sq = met(1,i,j,k)*met(1,i,j,k);
		  float_sw4 d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
		     3*u(1,i,j,k) -   u(1,i-1,j,k);
		  float_sw4 d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
		     3*u(1,i-1,j,k)-  u(1,i-2,j,k);
		  float_sw4 d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
		     3*kap(1,i,j,k) -   kap(1,i-1,j,k);
		  float_sw4 d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
		     3*kap(1,i-1,j,k) - kap(1,i-2,j,k);
		  float_sw4 d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
		     3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k);
		  float_sw4 d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
		     3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k);
		  float_sw4 d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
		     3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k);
		  float_sw4 d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
		     3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k);
		  float_sw4 pd = (c6*0.5*( d3up*d3kp+d3um*d3km  +
           dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*m1sq*wgh[k-1];
		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k)     = gmu(i,j,k) + 2*pd;
		  d3up = u(2,i+2,j,k)-3*u(2,i+1,j,k)+3*u(2,i,j,k)-
		     u(2,i-1,j,k);
		  d3um =u(2,i+1,j,k)-3*u(2,i,j,k)+3*u(2,i-1,j,k)-
		     u(2,i-2,j,k);
		  d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
		     3*kap(2,i,j,k)-kap(2,i-1,j,k);
		  d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
		     3*kap(2,i-1,j,k)-kap(2,i-2,j,k);
		  d3uap = uacc(2,i+2,j,k)-3*uacc(2,i+1,j,k)+
		     3*uacc(2,i,j,k)- uacc(2,i-1,j,k);
		  d3uam =uacc(2,i+1,j,k)-3*uacc(2,i,j,k)+
		     3*uacc(2,i-1,j,k)- uacc(2,i-2,j,k);
		  d3kap = kapacc(2,i+2,j,k)-3*kapacc(2,i+1,j,k)+
		     3*kapacc(2,i,j,k)-kapacc(2,i-1,j,k);
		  d3kam = kapacc(2,i+1,j,k)-3*kapacc(2,i,j,k)+
		     3*kapacc(2,i-1,j,k)-kapacc(2,i-2,j,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) ))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[k-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+3*u(3,i,j,k)-
		     u(3,i-1,j,k);
		  d3um =u(3,i+1,j,k)-3*u(3,i,j,k)+3*u(3,i-1,j,k)-
		     u(3,i-2,j,k);
		  d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
		     3*kap(3,i,j,k)-kap(3,i-1,j,k);
		  d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
		     3*kap(3,i-1,j,k)-kap(3,i-2,j,k);
		  d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
		     3*uacc(3,i,j,k)- uacc(3,i-1,j,k);
		  d3uam =uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
		     3*uacc(3,i-1,j,k)- uacc(3,i-2,j,k);
		  d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
		     3*kapacc(3,i,j,k)-kapacc(3,i-1,j,k);
		  d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
		     3*kapacc(3,i-1,j,k)-kapacc(3,i-2,j,k);
		  pd = (c6*( 0.5*( d3up*d3kp + d3um*d3km +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[k-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;

// y-direction
		  d3up = u(1,i,j+2,k)-3*u(1,i,j+1,k)+3*u(1,i,j,k)-
		                   u(1,i,j-1,k);
		  d3um =u(1,i,j+1,k)-3*u(1,i,j,k)+3*u(1,i,j-1,k)-
		     u(1,i,j-2,k);
		  d3kp = kap(1,i,j+2,k)-3*kap(1,i,j+1,k)+
		     3*kap(1,i,j,k)-kap(1,i,j-1,k);
		  d3km = kap(1,i,j+1,k)-3*kap(1,i,j,k)+
		     3*kap(1,i,j-1,k)-kap(1,i,j-2,k);
		  d3uap = uacc(1,i,j+2,k)-3*uacc(1,i,j+1,k)+
		     3*uacc(1,i,j,k)- uacc(1,i,j-1,k);
		  d3uam =uacc(1,i,j+1,k)-3*uacc(1,i,j,k)+
		     3*uacc(1,i,j-1,k)- uacc(1,i,j-2,k);
		  d3kap = kapacc(1,i,j+2,k)-3*kapacc(1,i,j+1,k)+
		     3*kapacc(1,i,j,k)-kapacc(1,i,j-1,k);
		  d3kam = kapacc(1,i,j+1,k)-3*kapacc(1,i,j,k)+
		     3*kapacc(1,i,j-1,k)-kapacc(1,i,j-2,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[k-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(2,i,j+2,k)-3*u(2,i,j+1,k)+3*u(2,i,j,k)-
		     u(2,i,j-1,k);
		  d3um =u(2,i,j+1,k)-3*u(2,i,j,k)+3*u(2,i,j-1,k)-
		     u(2,i,j-2,k);
		  d3kp = kap(2,i,j+2,k)-3*kap(2,i,j+1,k)+
		     3*kap(2,i,j,k)-kap(2,i,j-1,k);
		  d3km = kap(2,i,j+1,k)-3*kap(2,i,j,k)+
		     3*kap(2,i,j-1,k)-kap(2,i,j-2,k);
		  d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
		     3*uacc(2,i,j,k)-uacc(2,i,j-1,k);
		  d3uam =uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
		     3*uacc(2,i,j-1,k)-uacc(2,i,j-2,k);
		  d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
		     3*kapacc(2,i,j,k)-kapacc(2,i,j-1,k);
		  d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
		     3*kapacc(2,i,j-1,k)-kapacc(2,i,j-2,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[k-1];
		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k) = gmu(i,j,k) + 2*pd;

		  d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+3*u(3,i,j,k)-
                     u(3,i,j-1,k);
		  d3um =u(3,i,j+1,k)-3*u(3,i,j,k)+3*u(3,i,j-1,k)-
		     u(3,i,j-2,k);
		  d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
		     3*kap(3,i,j,k)-kap(3,i,j-1,k);
		  d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
		     3*kap(3,i,j-1,k)-kap(3,i,j-2,k);
		  d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
		     3*uacc(3,i,j,k)- uacc(3,i,j-1,k);
		  d3uam =uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
		     3*uacc(3,i,j-1,k)-uacc(3,i,j-2,k);
		  d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
		     3*kapacc(3,i,j,k)-kapacc(3,i,j-1,k);
		  d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
		     3*kapacc(3,i,j-1,k)-kapacc(3,i,j-2,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[k-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;
// z-direction
		  if( k>=2 )
		  {
// All derivatives are needed
		  d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
                     u(1,i,j,k-1);
 	          d3um = u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
                     u(1,i,j,k-2);
                  d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
		     3*kap(1,i,j,k)-kap(1,i,j,k-1);
                  d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
		     3*kap(1,i,j,k-1)-kap(1,i,j,k-2);
                  d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
		     3*uacc(1,i,j,k)-uacc(1,i,j,k-1);
                  d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
		     3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2);
                  d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
		     3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1);
                  d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
		     3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2);

                  float_sw4 d3vp = u(2,i,j,k+2)-3*u(2,i,j,k+1)+3*u(2,i,j,k)-
                     u(2,i,j,k-1);
                  float_sw4 d3vm = u(2,i,j,k+1)-3*u(2,i,j,k)+3*u(2,i,j,k-1)-
                     u(2,i,j,k-2);
                  float_sw4 d3lp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
		     3*kap(2,i,j,k)-kap(2,i,j,k-1);
                  float_sw4 d3lm = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
		     3*kap(2,i,j,k-1)-kap(2,i,j,k-2);
                  float_sw4 d3vap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
		     3*uacc(2,i,j,k)-uacc(2,i,j,k-1);
                  float_sw4 d3vam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
		     3*uacc(2,i,j,k-1)-uacc(2,i,j,k-2);
                  float_sw4 d3lap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
		     3*kapacc(2,i,j,k)-kapacc(2,i,j,k-1);
                  float_sw4 d3lam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
		     3*kapacc(2,i,j,k-1)-kapacc(2,i,j,k-2);

                  float_sw4 d3wp = u(3,i,j,k+2)-3*u(3,i,j,k+1)+3*u(3,i,j,k)-
                     u(3,i,j,k-1);
                  float_sw4 d3wm = u(3,i,j,k+1)-3*u(3,i,j,k)+3*u(3,i,j,k-1)-
                     u(3,i,j,k-2);
                  float_sw4 d3mp = kap(3,i,j,k+2)-3*kap(3,i,j,k+1)+
		     3*kap(3,i,j,k)-kap(3,i,j,k-1);
                  float_sw4 d3mm = kap(3,i,j,k+1)-3*kap(3,i,j,k)+
		     3*kap(3,i,j,k-1)-kap(3,i,j,k-2);
                  float_sw4 d3wap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
		     3*uacc(3,i,j,k)-uacc(3,i,j,k-1);
                  float_sw4 d3wam =uacc(3,i,j,k+1)-3*uacc(3,i,j,k)+
		     3*uacc(3,i,j,k-1)-uacc(3,i,j,k-2);
                  float_sw4 d3map = kapacc(3,i,j,k+2)-3*kapacc(3,i,j,k+1)+
		     3*kapacc(3,i,j,k)-kapacc(3,i,j,k-1);
                  float_sw4 d3mam = kapacc(3,i,j,k+1)-3*kapacc(3,i,j,k)+
		     3*kapacc(3,i,j,k-1)-kapacc(3,i,j,k-2);

#define SQR(x) (x)*(x)
                  float_sw4 mucof =SQR(met(2,i,j,k))+SQR(met(3,i,j,k))+SQR(met(4,i,j,k));
                  float_sw4 mucofp=SQR(met(2,i,j,k+1))+SQR(met(3,i,j,k+1))+SQR(met(4,i,j,k+1));
// u-u
                  pd = ( ( 0.5*( w6p[k-1]*d3up*d3kp+w6m[k-1]*d3um*d3km  +
             dt2o12*(w6p[k-1]*d3up*d3kap + w6m[k-1]*d3um*d3kam+ 
                     w6p[k-1]*d3uap*d3kp + w6m[k-1]*d3uam*d3km)))
                 + w8[k-1]*c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) );
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(2,i,j,k)));
                  glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(2,i,j,k));
                  if( k==4 )
		  {
                     pd = 0.5*( -al4*d3up*d3kp+
				  dt2o12*(-al4*d3up*d3kap -al4*d3uap*d3kp ) );
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  
			pd*(mucofp+SQR(met(2,i,j,k+1)));
                     glambda(i,j,k+1) = glambda(i,j,k+1) + 
			pd*SQR(met(2,i,j,k+1));
                  }
// u-v
                  pd = ( ( 0.5*( w6p[k-1]*d3vp*d3kp+w6m[k-1]*d3vm*d3km  +
             dt2o12*(w6p[k-1]*d3vp*d3kap + w6m[k-1]*d3vm*d3kam+ 
                     w6p[k-1]*d3vap*d3kp + w6m[k-1]*d3vam*d3km)))
                 + w8[k-1]*c8*( (d3vp-d3vm)*(d3kp-d3km) +
             dt2o12*(  (d3vap-d3vam)*(d3kp-d3km)+
                                 (d3vp-d3vm)*(d3kap-d3kam))) )
		     *met(2,i,j,k)*met(3,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==4 )
		  {
                     pd = 0.5*( -al4*d3vp*d3kp+
                    dt2o12*(-al4*d3vp*d3kap -al4*d3vap*d3kp ) )*
			met(2,i,j,k+1)*met(3,i,j,k+1);
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd;
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd;
                  }
// u-w 
                  pd = ( ( 0.5*( w6p[k-1]*d3wp*d3kp+w6m[k-1]*d3wm*d3km  +
             dt2o12*(w6p[k-1]*d3wp*d3kap + w6m[k-1]*d3wm*d3kam+ 
                     w6p[k-1]*d3wap*d3kp + w6m[k-1]*d3wam*d3km)))
                 + w8[k-1]*c8*( (d3wp-d3wm)*(d3kp-d3km) +
             dt2o12*(  (d3wap-d3wam)*(d3kp-d3km)+
                                 (d3wp-d3wm)*(d3kap-d3kam))) )
		     *met(2,i,j,k)*met(4,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==4 )
		  {
                     pd = 0.5*( -al4*d3wp*d3kp+
                        dt2o12*(-al4*d3wp*d3kap -al4*d3wap*d3kp ) )*
			           met(2,i,j,k+1)*met(4,i,j,k+1);
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd;
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd;
                  }
// v-u
                  pd = ( ( 0.5*( w6p[k-1]*d3up*d3lp+w6m[k-1]*d3um*d3lm  +
             dt2o12*(w6p[k-1]*d3up*d3lap + w6m[k-1]*d3um*d3lam+ 
                     w6p[k-1]*d3uap*d3lp + w6m[k-1]*d3uam*d3lm)))
                 + w8[k-1]*c8*( (d3up-d3um)*(d3lp-d3lm) +
             dt2o12*(  (d3uap-d3uam)*(d3lp-d3lm)+
                                 (d3up-d3um)*(d3lap-d3lam))) )
		     *met(2,i,j,k)*met(3,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==4 )
		  {
                     pd = 0.5*( -al4*d3up*d3lp+
                    dt2o12*(-al4*d3up*d3lap -al4*d3uap*d3lp ) )*
			met(2,i,j,k+1)*met(3,i,j,k+1);
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd;
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd;
                  }

// v-v
                  pd = ( ( 0.5*( w6p[k-1]*d3vp*d3lp+w6m[k-1]*d3vm*d3lm  +
             dt2o12*(w6p[k-1]*d3vp*d3lap + w6m[k-1]*d3vm*d3lam+ 
                     w6p[k-1]*d3vap*d3lp + w6m[k-1]*d3vam*d3lm)))
                 + w8[k-1]*c8*( (d3vp-d3vm)*(d3lp-d3lm) +
             dt2o12*(  (d3vap-d3vam)*(d3lp-d3lm)+
		       (d3vp-d3vm)*(d3lap-d3lam))) );
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(3,i,j,k)));
                  glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(3,i,j,k));
                  if( k==4 )
		  {
                     pd = 0.5*( -al4*d3vp*d3lp+
				  dt2o12*(-al4*d3vp*d3lap -al4*d3vap*d3lp ) );
                     glambda(i,j,k+1) = glambda(i,j,k+1) + 
			pd*SQR(met(3,i,j,k+1));
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  
			pd*(mucofp+SQR(met(3,i,j,k+1)));
                  }

// v-w
                  pd = ( ( 0.5*( w6p[k-1]*d3wp*d3lp+w6m[k-1]*d3wm*d3lm  +
             dt2o12*(w6p[k-1]*d3wp*d3lap + w6m[k-1]*d3wm*d3lam+ 
                     w6p[k-1]*d3wap*d3lp + w6m[k-1]*d3wam*d3lm)))
                + w8[k-1]*c8*( (d3wp-d3wm)*(d3lp-d3lm) +
             dt2o12*(  (d3wap-d3wam)*(d3lp-d3lm)+
                                 (d3wp-d3wm)*(d3lap-d3lam))) )
		     *met(3,i,j,k)*met(4,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==4 )
		  {
                     pd = 0.5*( -al4*d3wp*d3lp+
                    dt2o12*(-al4*d3wp*d3lap -al4*d3wap*d3lp ) )*
			met(3,i,j,k+1)*met(4,i,j,k+1);
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd;
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd;
                  }

// w-u
                  pd = ( ( 0.5*( w6p[k-1]*d3up*d3mp+w6m[k-1]*d3um*d3mm  +
             dt2o12*(w6p[k-1]*d3up*d3map + w6m[k-1]*d3um*d3mam+ 
                     w6p[k-1]*d3uap*d3mp + w6m[k-1]*d3uam*d3mm)))
                 + w8[k-1]*c8*( (d3up-d3um)*(d3mp-d3mm) +
             dt2o12*(  (d3uap-d3uam)*(d3mp-d3mm)+
                                 (d3up-d3um)*(d3map-d3mam))) )
		     *met(2,i,j,k)*met(4,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k == 4 )
		  {
                     pd = 0.5*( -al4*d3up*d3mp+
                    dt2o12*(-al4*d3up*d3map -al4*d3uap*d3mp ) )*
			met(2,i,j,k+1)*met(4,i,j,k+1);
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd;
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd;
                  }

// w-v
                  pd = ( ( 0.5*( w6p[k-1]*d3vp*d3mp+w6m[k-1]*d3vm*d3mm  +
             dt2o12*(w6p[k-1]*d3vp*d3map + w6m[k-1]*d3vm*d3mam+ 
                     w6p[k-1]*d3vap*d3mp + w6m[k-1]*d3vam*d3mm)))
                 + w8[k-1]*c8*( (d3vp-d3vm)*(d3mp-d3mm) +
             dt2o12*(  (d3vap-d3vam)*(d3mp-d3mm)+
                                 (d3vp-d3vm)*(d3map-d3mam))) )
		     *met(3,i,j,k)*met(4,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==4 )
		  {
                     pd = 0.5*( -al4*d3vp*d3mp+
                    dt2o12*(-al4*d3vp*d3map -al4*d3vap*d3mp ) )*
			met(3,i,j,k+1)*met(4,i,j,k+1);
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd;
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd;
                  }
// w-w
                  pd = ( ( 0.5*( w6p[k-1]*d3wp*d3mp+w6m[k-1]*d3wm*d3mm  +
             dt2o12*(w6p[k-1]*d3wp*d3map + w6m[k-1]*d3wm*d3mam+ 
                     w6p[k-1]*d3wap*d3mp + w6m[k-1]*d3wam*d3mm)))
                 + w8[k-1]*c8*( (d3wp-d3wm)*(d3mp-d3mm) +
             dt2o12*(  (d3wap-d3wam)*(d3mp-d3mm)+
		       (d3wp-d3wm)*(d3map-d3mam))) );
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(4,i,j,k)));
                  glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(4,i,j,k));

		  if( k == 4 )
		  {
                     pd = 0.5*( -al4*d3wp*d3mp+
				  dt2o12*(-al4*d3wp*d3map -al4*d3wap*d3mp ) );
                     glambda(i,j,k+1) = glambda(i,j,k+1) + 
			pd*SQR(met(4,i,j,k+1));
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  
			pd*(mucofp+SQR(met(4,i,j,k+1)));
		  }
		  }
	       }
      }


      //#pragma omp parallel for
      for( int k=kstart; k <= kend; k++ )
	 for( int j=jfirstact; j <= jlastact; j++ )
            //#pragma ivdep
	    for( int i=ifirstact; i <= ilastact; i++ )
	    {
// Diagonal terms
               float_sw4 dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
		  d4a*(u(1,i+1,j,k)-u(1,i-1,j,k));
               float_sw4 dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
		  d4a*(u(2,i,j+1,k)-u(2,i,j-1,k));
               float_sw4 dwz = d4b*(u(3,i,j,k+2)-u(3,i,j,k-2))+
		  d4a*(u(3,i,j,k+1)-u(3,i,j,k-1));
               float_sw4 dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
		  d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k));
               float_sw4 dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
		  d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k));
               float_sw4 dmz = d4b*(kap(3,i,j,k+2)-kap(3,i,j,k-2))+
		  d4a*(kap(3,i,j,k+1)-kap(3,i,j,k-1));

               float_sw4 duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
		  d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k));
               float_sw4 dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
		  d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k));
               float_sw4 dwaz = d4b*(uacc(3,i,j,k+2)-uacc(3,i,j,k-2))+
		  d4a*(uacc(3,i,j,k+1)-uacc(3,i,j,k-1));
               float_sw4 dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
		  d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k));
               float_sw4 dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
		  d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k));
               float_sw4 dmaz = d4b*(kapacc(3,i,j,k+2)-kapacc(3,i,j,k-2))+
		  d4a*(kapacc(3,i,j,k+1)-kapacc(3,i,j,k-1));

               float_sw4 duz = d4b*(u(1,i,j,k+2)-u(1,i,j,k-2))+
		  d4a*(u(1,i,j,k+1)-u(1,i,j,k-1));
               float_sw4 dvz = d4b*(u(2,i,j,k+2)-u(2,i,j,k-2))+
		  d4a*(u(2,i,j,k+1)-u(2,i,j,k-1));

               float_sw4 dkz = d4b*(kap(1,i,j,k+2)-kap(1,i,j,k-2))+
		  d4a*(kap(1,i,j,k+1)-kap(1,i,j,k-1));
               float_sw4 dlz = d4b*(kap(2,i,j,k+2)-kap(2,i,j,k-2))+
		  d4a*(kap(2,i,j,k+1)-kap(2,i,j,k-1));

               float_sw4 duaz= d4b*(uacc(1,i,j,k+2)-uacc(1,i,j,k-2))+
		  d4a*(uacc(1,i,j,k+1)-uacc(1,i,j,k-1));
               float_sw4 dvaz= d4b*(uacc(2,i,j,k+2)-uacc(2,i,j,k-2))+
		  d4a*(uacc(2,i,j,k+1)-uacc(2,i,j,k-1));

               float_sw4 dkaz= d4b*(kapacc(1,i,j,k+2)-kapacc(1,i,j,k-2))+
		  d4a*(kapacc(1,i,j,k+1)-kapacc(1,i,j,k-1));
               float_sw4 dlaz= d4b*(kapacc(2,i,j,k+2)-kapacc(2,i,j,k-2))+
		  d4a*(kapacc(2,i,j,k+1)-kapacc(2,i,j,k-1));

               dux = met(1,i,j,k)*dux + met(2,i,j,k)*duz;
               dvy = met(1,i,j,k)*dvy + met(3,i,j,k)*dvz;
               dkx = met(1,i,j,k)*dkx + met(2,i,j,k)*dkz;
               dly = met(1,i,j,k)*dly + met(3,i,j,k)*dlz;
               duax = met(1,i,j,k)*duax + met(2,i,j,k)*duaz;
               dvay = met(1,i,j,k)*dvay + met(3,i,j,k)*dvaz;
               dkax = met(1,i,j,k)*dkax + met(2,i,j,k)*dkaz;
               dlay = met(1,i,j,k)*dlay + met(3,i,j,k)*dlaz;

               float_sw4 dwzm  = dwz*met(4,i,j,k);    
               float_sw4 dmzm  = dmz*met(4,i,j,k);
               float_sw4 dwazm = dwaz*met(4,i,j,k);
               float_sw4 dmazm = dmaz*met(4,i,j,k);

               glambda(i,j,k) = glambda(i,j,k) + (
                (dux+dvy+dwzm)*(dkx+dly+dmzm) +
                dt2o12*( (duax+dvay+dwazm)*(dkx+dly+dmzm) +
                         (dux+dvy+dwzm)*(dkax+dlay+dmazm) )
						  );
               gmu(i,j,k) = gmu(i,j,k) + 2*(
               dux*dkx+dvy*dly+dwzm*dmzm +
                dt2o12*( (duax*dkx+dvay*dly+dwazm*dmzm)+
                         (dux*dkax+dvy*dlay+dwzm*dmazm) )
					    );
// Off diagonal stresses
//  xy
               float_sw4 stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
                      d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
                      d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
		  d4a*(u(1,i,j+1,k)-u(1,i,j-1,k));
               float_sw4 stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
                      d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
                      d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
		  d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k));
               float_sw4 stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
                      d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
                      d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
		  d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k));
               float_sw4 stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
                      d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
                      d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
		  d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k));

               stuxy = met(1,i,j,k)*stuxy+met(2,i,j,k)*dvz +
		  met(3,i,j,k)*duz;
               stkxy = met(1,i,j,k)*stkxy+met(2,i,j,k)*dlz +
		  met(3,i,j,k)*dkz;
               stuaxy= met(1,i,j,k)*stuaxy+met(2,i,j,k)*dvaz +
		  met(3,i,j,k)*duaz;
               stkaxy= met(1,i,j,k)*stkaxy+met(2,i,j,k)*dlaz +
		  met(3,i,j,k)*dkaz;

               gmu(i,j,k)= gmu(i,j,k) + (
		 stuxy*stkxy + dt2o12*(stuxy*stkaxy+stuaxy*stkxy) );

//  xz
               float_sw4 stuxz = met(1,i,j,k)*(d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
					       d4a*(u(3,i+1,j,k)-u(3,i-1,j,k)) );
               float_sw4 stkxz = met(1,i,j,k)*(
                      d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
                      d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k)) );
               float_sw4 stuaxz = met(1,i,j,k)*(
                      d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
                      d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k)) );
               float_sw4 stkaxz = met(1,i,j,k)*(
                      d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
                      d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k)) );

               stuxz = stuxz + met(2,i,j,k)*dwz+met(4,i,j,k)*duz;
               stkxz = stkxz + met(2,i,j,k)*dmz+met(4,i,j,k)*dkz;
               stuaxz= stuaxz+ met(2,i,j,k)*dwaz+met(4,i,j,k)*duaz;
               stkaxz= stkaxz+ met(2,i,j,k)*dmaz+met(4,i,j,k)*dkaz;

               gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
					  dt2o12*(stuxz*stkaxz + stuaxz*stkxz) );
//  yz
               float_sw4 stuyz = met(1,i,j,k)*(d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
					       d4a*(u(3,i,j+1,k)-u(3,i,j-1,k)));
               float_sw4 stkyz = met(1,i,j,k)*(
                      d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
                      d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k)) );
               float_sw4 stuayz = met(1,i,j,k)*(
                      d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
                      d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k)) );
               float_sw4 stkayz = met(1,i,j,k)*(
                      d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
                      d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k)) );

               stuyz = stuyz + met(3,i,j,k)*dwz+met(4,i,j,k)*dvz;
               stkyz = stkyz + met(3,i,j,k)*dmz+met(4,i,j,k)*dlz;
               stuayz= stuayz+ met(3,i,j,k)*dwaz+met(4,i,j,k)*dvaz;
               stkayz= stkayz+ met(3,i,j,k)*dmaz+met(4,i,j,k)*dlaz;

               gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
				  dt2o12*( stuyz*stkayz+stuayz*stkyz) );

// Pos. def extra terms
// x-direction
//   metric is constant for the x and y directions, just use as a norm weight
               float_sw4 m1sq = met(1,i,j,k)*met(1,i,j,k);

               float_sw4 d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
		  3*u(1,i,j,k) -   u(1,i-1,j,k);
               float_sw4 d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
                  3*u(1,i-1,j,k)-  u(1,i-2,j,k);
               float_sw4 d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
		  3*kap(1,i,j,k) -   kap(1,i-1,j,k);
               float_sw4 d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
		  3*kap(1,i-1,j,k) - kap(1,i-2,j,k);
               float_sw4 d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
		  3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k);
               float_sw4 d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
		  3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k);
               float_sw4 d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
		  3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k);
               float_sw4 d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
		  3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k);
               float_sw4 pd = (c6*0.5*( d3up*d3kp+d3um*d3km  +
           dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*m1sq;
               glambda(i,j,k) = glambda(i,j,k) + pd;
               gmu(i,j,k)     = gmu(i,j,k) + 2*pd;
               d3up = u(2,i+2,j,k)-3*u(2,i+1,j,k)+3*u(2,i,j,k)-
		  u(2,i-1,j,k);
               d3um =u(2,i+1,j,k)-3*u(2,i,j,k)+3*u(2,i-1,j,k)-
		  u(2,i-2,j,k);
               d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
		  3*kap(2,i,j,k)-kap(2,i-1,j,k);
               d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
		  3*kap(2,i-1,j,k)-kap(2,i-2,j,k);
               d3uap = uacc(2,i+2,j,k)-3*uacc(2,i+1,j,k)+
		  3*uacc(2,i,j,k)- uacc(2,i-1,j,k);
               d3uam =uacc(2,i+1,j,k)-3*uacc(2,i,j,k)+
		  3*uacc(2,i-1,j,k)- uacc(2,i-2,j,k);
               d3kap = kapacc(2,i+2,j,k)-3*kapacc(2,i+1,j,k)+
		  3*kapacc(2,i,j,k)-kapacc(2,i-1,j,k);
               d3kam = kapacc(2,i+1,j,k)-3*kapacc(2,i,j,k)+
		  3*kapacc(2,i-1,j,k)-kapacc(2,i-2,j,k);
               pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) ))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq;
               gmu(i,j,k) = gmu(i,j,k) + pd;

               d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+3*u(3,i,j,k)-
		  u(3,i-1,j,k);
               d3um =u(3,i+1,j,k)-3*u(3,i,j,k)+3*u(3,i-1,j,k)-
		  u(3,i-2,j,k);
               d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
		  3*kap(3,i,j,k)-kap(3,i-1,j,k);
               d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
		  3*kap(3,i-1,j,k)-kap(3,i-2,j,k);
               d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
		  3*uacc(3,i,j,k)- uacc(3,i-1,j,k);
               d3uam =uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
                  3*uacc(3,i-1,j,k)- uacc(3,i-2,j,k);
               d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
		  3*kapacc(3,i,j,k)-kapacc(3,i-1,j,k);
               d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
		  3*kapacc(3,i-1,j,k)-kapacc(3,i-2,j,k);
               pd = (c6*( 0.5*( d3up*d3kp + d3um*d3km +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq;
               gmu(i,j,k) = gmu(i,j,k) + pd;

// y-direction
               d3up = u(1,i,j+2,k)-3*u(1,i,j+1,k)+3*u(1,i,j,k)-
		  u(1,i,j-1,k);
               d3um =u(1,i,j+1,k)-3*u(1,i,j,k)+3*u(1,i,j-1,k)-
		  u(1,i,j-2,k);
               d3kp = kap(1,i,j+2,k)-3*kap(1,i,j+1,k)+
		  3*kap(1,i,j,k)-kap(1,i,j-1,k);
               d3km = kap(1,i,j+1,k)-3*kap(1,i,j,k)+
		  3*kap(1,i,j-1,k)-kap(1,i,j-2,k);
               d3uap = uacc(1,i,j+2,k)-3*uacc(1,i,j+1,k)+
                  3*uacc(1,i,j,k)- uacc(1,i,j-1,k);
               d3uam =uacc(1,i,j+1,k)-3*uacc(1,i,j,k)+
                  3*uacc(1,i,j-1,k)- uacc(1,i,j-2,k);
               d3kap = kapacc(1,i,j+2,k)-3*kapacc(1,i,j+1,k)+
		  3*kapacc(1,i,j,k)-kapacc(1,i,j-1,k);
               d3kam = kapacc(1,i,j+1,k)-3*kapacc(1,i,j,k)+
		  3*kapacc(1,i,j-1,k)-kapacc(1,i,j-2,k);
               pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq;
               gmu(i,j,k) = gmu(i,j,k) + pd;

               d3up = u(2,i,j+2,k)-3*u(2,i,j+1,k)+3*u(2,i,j,k)-
		  u(2,i,j-1,k);
               d3um =u(2,i,j+1,k)-3*u(2,i,j,k)+3*u(2,i,j-1,k)-
		  u(2,i,j-2,k);
               d3kp = kap(2,i,j+2,k)-3*kap(2,i,j+1,k)+
		  3*kap(2,i,j,k)-kap(2,i,j-1,k);
               d3km = kap(2,i,j+1,k)-3*kap(2,i,j,k)+
		  3*kap(2,i,j-1,k)-kap(2,i,j-2,k);
               d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
		  3*uacc(2,i,j,k)-uacc(2,i,j-1,k);
               d3uam =uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
		  3*uacc(2,i,j-1,k)-uacc(2,i,j-2,k);
               d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
		  3*kapacc(2,i,j,k)-kapacc(2,i,j-1,k);
               d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
		  3*kapacc(2,i,j-1,k)-kapacc(2,i,j-2,k);
               pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq;
               glambda(i,j,k) = glambda(i,j,k) + pd;
               gmu(i,j,k) = gmu(i,j,k) + 2*pd;

               d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+3*u(3,i,j,k)-
		  u(3,i,j-1,k);
               d3um =u(3,i,j+1,k)-3*u(3,i,j,k)+3*u(3,i,j-1,k)-
		  u(3,i,j-2,k);
               d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
		  3*kap(3,i,j,k)-kap(3,i,j-1,k);
               d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
		  3*kap(3,i,j-1,k)-kap(3,i,j-2,k);
               d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
		  3*uacc(3,i,j,k)- uacc(3,i,j-1,k);
               d3uam =uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
		  3*uacc(3,i,j-1,k)-uacc(3,i,j-2,k);
               d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
		  3*kapacc(3,i,j,k)-kapacc(3,i,j-1,k);
               d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
		  3*kapacc(3,i,j-1,k)-kapacc(3,i,j-2,k);
               pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq;
               gmu(i,j,k) = gmu(i,j,k) + pd;

// z-direction
// All derivatives are needed
               d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
		  u(1,i,j,k-1);
               d3um = u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
		  u(1,i,j,k-2);
               d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
		  3*kap(1,i,j,k)-kap(1,i,j,k-1);
               d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
		  3*kap(1,i,j,k-1)-kap(1,i,j,k-2);
               d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
		  3*uacc(1,i,j,k)-uacc(1,i,j,k-1);
               d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
		  3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2);
               d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
		  3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1);
               d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
		  3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2);

               float_sw4 d3vp = u(2,i,j,k+2)-3*u(2,i,j,k+1)+3*u(2,i,j,k)-
		  u(2,i,j,k-1);
               float_sw4 d3vm = u(2,i,j,k+1)-3*u(2,i,j,k)+3*u(2,i,j,k-1)-
		  u(2,i,j,k-2);
               float_sw4 d3lp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
		  3*kap(2,i,j,k)-kap(2,i,j,k-1);
               float_sw4 d3lm = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
		  3*kap(2,i,j,k-1)-kap(2,i,j,k-2);
               float_sw4 d3vap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
		  3*uacc(2,i,j,k)-uacc(2,i,j,k-1);
               float_sw4 d3vam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
		  3*uacc(2,i,j,k-1)-uacc(2,i,j,k-2);
               float_sw4 d3lap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
		  3*kapacc(2,i,j,k)-kapacc(2,i,j,k-1);
               float_sw4 d3lam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
		  3*kapacc(2,i,j,k-1)-kapacc(2,i,j,k-2);

               float_sw4 d3wp = u(3,i,j,k+2)-3*u(3,i,j,k+1)+3*u(3,i,j,k)-
		  u(3,i,j,k-1);
               float_sw4 d3wm = u(3,i,j,k+1)-3*u(3,i,j,k)+3*u(3,i,j,k-1)-
		  u(3,i,j,k-2);
               float_sw4 d3mp = kap(3,i,j,k+2)-3*kap(3,i,j,k+1)+
		  3*kap(3,i,j,k)-kap(3,i,j,k-1);
               float_sw4 d3mm = kap(3,i,j,k+1)-3*kap(3,i,j,k)+
		  3*kap(3,i,j,k-1)-kap(3,i,j,k-2);
               float_sw4 d3wap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
		  3*uacc(3,i,j,k)-uacc(3,i,j,k-1);
               float_sw4 d3wam =uacc(3,i,j,k+1)-3*uacc(3,i,j,k)+
		  3*uacc(3,i,j,k-1)-uacc(3,i,j,k-2);
               float_sw4 d3map = kapacc(3,i,j,k+2)-3*kapacc(3,i,j,k+1)+
		  3*kapacc(3,i,j,k)-kapacc(3,i,j,k-1);
               float_sw4 d3mam = kapacc(3,i,j,k+1)-3*kapacc(3,i,j,k)+
		  3*kapacc(3,i,j,k-1)-kapacc(3,i,j,k-2);

               float_sw4 mucof=SQR(met(2,i,j,k))+SQR(met(3,i,j,k))+SQR(met(4,i,j,k));
// u-u
               pd = ( c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap + d3um*d3kam+ 
                     d3uap*d3kp + d3uam*d3km)))
                 +  c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) );
               gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(2,i,j,k)));
               glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(2,i,j,k));

// u-v
               pd = ( c6*( 0.5*( d3vp*d3kp+d3vm*d3km  +
             dt2o12*(d3vp*d3kap + d3vm*d3kam+ 
                     d3vap*d3kp + d3vam*d3km)))
                 + c8*( (d3vp-d3vm)*(d3kp-d3km) +
             dt2o12*(  (d3vap-d3vam)*(d3kp-d3km)+
                                 (d3vp-d3vm)*(d3kap-d3kam))) )
		  *met(2,i,j,k)*met(3,i,j,k);
               gmu(i,j,k) = gmu(i,j,k) + pd;
               glambda(i,j,k) = glambda(i,j,k) + pd;

// u-w 
               pd = ( c6*( 0.5*( d3wp*d3kp+d3wm*d3km  +
             dt2o12*(d3wp*d3kap + d3wm*d3kam+ 
                     d3wap*d3kp + d3wam*d3km)))
                 + c8*( (d3wp-d3wm)*(d3kp-d3km) +
             dt2o12*(  (d3wap-d3wam)*(d3kp-d3km)+
                                 (d3wp-d3wm)*(d3kap-d3kam))) )
		  *met(2,i,j,k)*met(4,i,j,k);
               gmu(i,j,k) = gmu(i,j,k) + pd;
               glambda(i,j,k) = glambda(i,j,k) + pd;

// v-u
               pd = ( c6*( 0.5*( d3up*d3lp+d3um*d3lm  +
             dt2o12*(d3up*d3lap + d3um*d3lam+ 
                     d3uap*d3lp + d3uam*d3lm)))
                 + c8*( (d3up-d3um)*(d3lp-d3lm) +
             dt2o12*(  (d3uap-d3uam)*(d3lp-d3lm)+
                                 (d3up-d3um)*(d3lap-d3lam))) )
		  *met(2,i,j,k)*met(3,i,j,k);
               gmu(i,j,k) = gmu(i,j,k) + pd;
               glambda(i,j,k) = glambda(i,j,k) + pd;

// v-v
               pd = ( c6*( 0.5*( d3vp*d3lp+d3vm*d3lm  +
             dt2o12*(d3vp*d3lap + d3vm*d3lam+ 
                     d3vap*d3lp + d3vam*d3lm)))
                 + c8*( (d3vp-d3vm)*(d3lp-d3lm) +
             dt2o12*(  (d3vap-d3vam)*(d3lp-d3lm)+
		       (d3vp-d3vm)*(d3lap-d3lam))) );
               gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(3,i,j,k)));
               glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(3,i,j,k));

// v-w
               pd = ( c6*( 0.5*( d3wp*d3lp+d3wm*d3lm  +
             dt2o12*(d3wp*d3lap + d3wm*d3lam+ 
                     d3wap*d3lp + d3wam*d3lm)))
                 + c8*( (d3wp-d3wm)*(d3lp-d3lm) +
             dt2o12*(  (d3wap-d3wam)*(d3lp-d3lm)+
                                 (d3wp-d3wm)*(d3lap-d3lam))) )
		  *met(3,i,j,k)*met(4,i,j,k);
               gmu(i,j,k) = gmu(i,j,k) + pd;
               glambda(i,j,k) = glambda(i,j,k) + pd;

// w-u
               pd = ( c6*( 0.5*( d3up*d3mp+d3um*d3mm  +
             dt2o12*(d3up*d3map + d3um*d3mam+ 
                     d3uap*d3mp + d3uam*d3mm)))
                 + c8*( (d3up-d3um)*(d3mp-d3mm) +
             dt2o12*(  (d3uap-d3uam)*(d3mp-d3mm)+
                                 (d3up-d3um)*(d3map-d3mam))) )
		  *met(2,i,j,k)*met(4,i,j,k);
               gmu(i,j,k) = gmu(i,j,k) + pd;
               glambda(i,j,k) = glambda(i,j,k) + pd;

// w-v
               pd = ( c6*( 0.5*( d3vp*d3mp+d3vm*d3mm  +
             dt2o12*(d3vp*d3map + d3vm*d3mam+ 
                     d3vap*d3mp + d3vam*d3mm)))
                 + c8*( (d3vp-d3vm)*(d3mp-d3mm) +
             dt2o12*(  (d3vap-d3vam)*(d3mp-d3mm)+
                                 (d3vp-d3vm)*(d3map-d3mam))) )
		  *met(3,i,j,k)*met(4,i,j,k);
               gmu(i,j,k) = gmu(i,j,k) + pd;
               glambda(i,j,k) = glambda(i,j,k) + pd;

// w-w
               pd = ( c6*( 0.5*( d3wp*d3mp+d3wm*d3mm  +
             dt2o12*(d3wp*d3map + d3wm*d3mam+ 
                     d3wap*d3mp + d3wam*d3mm)))
                 + c8*( (d3wp-d3wm)*(d3mp-d3mm) +
             dt2o12*(  (d3wap-d3wam)*(d3mp-d3mm)+
		       (d3wp-d3wm)*(d3map-d3mam))) );
               gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(4,i,j,k)));
               glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(4,i,j,k));
	    }
	 }
      if( klastact >= nk-3 && onesided[5] == 1 )
      {
	 float_sw4 w8[4]={0,0,1,1};
         float_sw4 w6p[4]={0,0,al1,al1+al2};
	 float_sw4 w6m[4]={0,al1,al1+al2,al1+al2+al3};
	 for( int k=nk-3; k <= klastact; k++ )
#pragma omp parallel for
            for( int j=jfirstact; j <= jlastact; j++ )
#pragma ivdep
               for( int i=ifirstact; i <= ilastact; i++ )
	       {
                  int kk=nk-k+1;
		  float_sw4 normfact = wgh[kk-1];
// Diagonal terms
		  float_sw4 dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
		     d4a*(u(1,i+1,j,k)-u(1,i-1,j,k));
		  float_sw4 dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
		     d4a*(u(2,i,j+1,k)-u(2,i,j-1,k));
		  float_sw4 dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
		     d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k));
		  float_sw4 dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
		     d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k));

		  float_sw4 duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
		     d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k));
		  float_sw4 dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
		     d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k));
		  float_sw4 dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
		     d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k));
		  float_sw4 dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
		     d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k));
		  float_sw4 duz = 0, dvz = 0, dwz = 0, dkz = 0, dlz = 0, dmz = 0;
		  float_sw4 duaz=0, dvaz=0, dwaz=0, dkaz=0, dlaz=0, dmaz=0;
		  for( int m=1; m <= wb ;m++ )
		  {
		     duz -= bop(kk,m)*u(1,i,j,nk-m+1);
		     dvz -= bop(kk,m)*u(2,i,j,nk-m+1);
		     dwz -= bop(kk,m)*u(3,i,j,nk-m+1);

		     dkz -= bop(kk,m)*kap(1,i,j,nk-m+1);
		     dlz -= bop(kk,m)*kap(2,i,j,nk-m+1);
		     dmz -= bop(kk,m)*kap(3,i,j,nk-m+1);

		     duaz -= bop(kk,m)*uacc(1,i,j,nk-m+1);
		     dvaz -= bop(kk,m)*uacc(2,i,j,nk-m+1);
		     dwaz -= bop(kk,m)*uacc(3,i,j,nk-m+1);

		     dkaz -= bop(kk,m)*kapacc(1,i,j,nk-m+1);
		     dlaz -= bop(kk,m)*kapacc(2,i,j,nk-m+1);
		     dmaz -= bop(kk,m)*kapacc(3,i,j,nk-m+1);
		  }
		  dux = met(1,i,j,k)*dux + met(2,i,j,k)*duz;
		  dvy = met(1,i,j,k)*dvy + met(3,i,j,k)*dvz;
		  dkx = met(1,i,j,k)*dkx + met(2,i,j,k)*dkz;
		  dly = met(1,i,j,k)*dly + met(3,i,j,k)*dlz;
		  duax = met(1,i,j,k)*duax + met(2,i,j,k)*duaz;
		  dvay = met(1,i,j,k)*dvay + met(3,i,j,k)*dvaz;
		  dkax = met(1,i,j,k)*dkax + met(2,i,j,k)*dkaz;
		  dlay = met(1,i,j,k)*dlay + met(3,i,j,k)*dlaz;

		  float_sw4 dwzm  = dwz*met(4,i,j,k);
		  float_sw4 dmzm  = dmz*met(4,i,j,k);
		  float_sw4 dwazm = dwaz*met(4,i,j,k);
		  float_sw4 dmazm = dmaz*met(4,i,j,k);

		  glambda(i,j,k) = glambda(i,j,k) + (
                (dux+dvy+dwzm)*(dkx+dly+dmzm) +
                dt2o12*( (duax+dvay+dwazm)*(dkx+dly+dmzm) +
                         (dux+dvy+dwzm)*(dkax+dlay+dmazm) )
						     )*normfact;
		  gmu(i,j,k) = gmu(i,j,k) + 2*(
               dux*dkx+dvy*dly+dwzm*dmzm +
                dt2o12*( (duax*dkx+dvay*dly+dwazm*dmzm)+
                         (dux*dkax+dvy*dlay+dwzm*dmazm) )
					       )*normfact;
// Off diagonal stresses
//  xy
		  float_sw4 stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
                      d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
                      d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
		     d4a*(u(1,i,j+1,k)-u(1,i,j-1,k));
		  float_sw4 stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
                      d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
                      d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
		     d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k));
		  float_sw4 stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
                      d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
                      d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
		     d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k));
		  float_sw4 stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
                      d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
                      d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
		     d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k));

		  stuxy = met(1,i,j,k)*stuxy+met(2,i,j,k)*dvz +
		     met(3,i,j,k)*duz;
		  stkxy = met(1,i,j,k)*stkxy+met(2,i,j,k)*dlz +
		     met(3,i,j,k)*dkz;
		  stuaxy= met(1,i,j,k)*stuaxy+met(2,i,j,k)*dvaz +
		     met(3,i,j,k)*duaz;
		  stkaxy= met(1,i,j,k)*stkaxy+met(2,i,j,k)*dlaz +
		     met(3,i,j,k)*dkaz;

		  gmu(i,j,k)= gmu(i,j,k) + (
           stuxy*stkxy + dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )*
		     normfact;
//  xz
		  float_sw4 stuxz = met(1,i,j,k)*(
		      d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
		      d4a*(u(3,i+1,j,k)-u(3,i-1,j,k)) );
		  float_sw4 stkxz = met(1,i,j,k)*(
		      d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
                      d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k)) );
		  float_sw4 stuaxz = met(1,i,j,k)*(
                      d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
                      d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k)) );
		  float_sw4 stkaxz = met(1,i,j,k)*(
                      d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
                      d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k)) );

		  stuxz = stuxz + met(2,i,j,k)*dwz+met(4,i,j,k)*duz;
		  stkxz = stkxz + met(2,i,j,k)*dmz+met(4,i,j,k)*dkz;
		  stuaxz= stuaxz+ met(2,i,j,k)*dwaz+met(4,i,j,k)*duaz;
		  stkaxz= stkaxz+ met(2,i,j,k)*dmaz+met(4,i,j,k)*dkaz;

		  gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
		       dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )*normfact;
//  yz
		  float_sw4 stuyz = met(1,i,j,k)*(d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
						  d4a*(u(3,i,j+1,k)-u(3,i,j-1,k)));
		  float_sw4 stkyz = met(1,i,j,k)*(
                      d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
                      d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k)) );
		  float_sw4 stuayz = met(1,i,j,k)*(
                      d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
                      d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k)) );
		  float_sw4 stkayz = met(1,i,j,k)*(
                      d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
                      d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k)) );

		  stuyz = stuyz + met(3,i,j,k)*dwz+met(4,i,j,k)*dvz;
		  stkyz = stkyz + met(3,i,j,k)*dmz+met(4,i,j,k)*dlz;
		  stuayz= stuayz+ met(3,i,j,k)*dwaz+met(4,i,j,k)*dvaz;
		  stkayz= stkayz+ met(3,i,j,k)*dmaz+met(4,i,j,k)*dlaz;

		  gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
 	               dt2o12*( stuyz*stkayz+stuayz*stkyz) )*normfact;
// Pos. def extra terms
// x-direction
//   metric is constant for the x and y directions, just use as a norm weight
		  float_sw4 m1sq = met(1,i,j,k)*met(1,i,j,k);
		  float_sw4 d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
		     3*u(1,i,j,k) -   u(1,i-1,j,k);
		  float_sw4 d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
		     3*u(1,i-1,j,k)-  u(1,i-2,j,k);
		  float_sw4 d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
		     3*kap(1,i,j,k) -   kap(1,i-1,j,k);
		  float_sw4 d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
		     3*kap(1,i-1,j,k) - kap(1,i-2,j,k);
		  float_sw4 d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
		     3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k);
		  float_sw4 d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
		     3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k);
		  float_sw4 d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
		     3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k);
		  float_sw4 d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
		     3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k);
		  float_sw4 pd = (c6*0.5*( d3up*d3kp+d3um*d3km  +
           dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam)) ) )*m1sq*wgh[kk-1];
		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k)     = gmu(i,j,k) + 2*pd;
		  d3up = u(2,i+2,j,k)-3*u(2,i+1,j,k)+3*u(2,i,j,k)-
		     u(2,i-1,j,k);
		  d3um =u(2,i+1,j,k)-3*u(2,i,j,k)+3*u(2,i-1,j,k)-
		     u(2,i-2,j,k);
		  d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
		     3*kap(2,i,j,k)-kap(2,i-1,j,k);
		  d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
		     3*kap(2,i-1,j,k)-kap(2,i-2,j,k);
		  d3uap = uacc(2,i+2,j,k)-3*uacc(2,i+1,j,k)+
		     3*uacc(2,i,j,k)- uacc(2,i-1,j,k);
		  d3uam =uacc(2,i+1,j,k)-3*uacc(2,i,j,k)+
		     3*uacc(2,i-1,j,k)- uacc(2,i-2,j,k);
		  d3kap = kapacc(2,i+2,j,k)-3*kapacc(2,i+1,j,k)+
		     3*kapacc(2,i,j,k)-kapacc(2,i-1,j,k);
		  d3kam = kapacc(2,i+1,j,k)-3*kapacc(2,i,j,k)+
		     3*kapacc(2,i-1,j,k)-kapacc(2,i-2,j,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) ))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[kk-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+3*u(3,i,j,k)-
		     u(3,i-1,j,k);
		  d3um =u(3,i+1,j,k)-3*u(3,i,j,k)+3*u(3,i-1,j,k)-
		     u(3,i-2,j,k);
		  d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
		     3*kap(3,i,j,k)-kap(3,i-1,j,k);
		  d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
		     3*kap(3,i-1,j,k)-kap(3,i-2,j,k);
		  d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
		     3*uacc(3,i,j,k)- uacc(3,i-1,j,k);
		  d3uam =uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
		     3*uacc(3,i-1,j,k)- uacc(3,i-2,j,k);
		  d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
		     3*kapacc(3,i,j,k)-kapacc(3,i-1,j,k);
		  d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
		     3*kapacc(3,i-1,j,k)-kapacc(3,i-2,j,k);
		  pd = (c6*( 0.5*( d3up*d3kp + d3um*d3km +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[kk-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;

// y-direction
		  d3up = u(1,i,j+2,k)-3*u(1,i,j+1,k)+3*u(1,i,j,k)-
		                   u(1,i,j-1,k);
		  d3um =u(1,i,j+1,k)-3*u(1,i,j,k)+3*u(1,i,j-1,k)-
		     u(1,i,j-2,k);
		  d3kp = kap(1,i,j+2,k)-3*kap(1,i,j+1,k)+
		     3*kap(1,i,j,k)-kap(1,i,j-1,k);
		  d3km = kap(1,i,j+1,k)-3*kap(1,i,j,k)+
		     3*kap(1,i,j-1,k)-kap(1,i,j-2,k);
		  d3uap = uacc(1,i,j+2,k)-3*uacc(1,i,j+1,k)+
		     3*uacc(1,i,j,k)- uacc(1,i,j-1,k);
		  d3uam =uacc(1,i,j+1,k)-3*uacc(1,i,j,k)+
		     3*uacc(1,i,j-1,k)- uacc(1,i,j-2,k);
		  d3kap = kapacc(1,i,j+2,k)-3*kapacc(1,i,j+1,k)+
		     3*kapacc(1,i,j,k)-kapacc(1,i,j-1,k);
		  d3kam = kapacc(1,i,j+1,k)-3*kapacc(1,i,j,k)+
		     3*kapacc(1,i,j-1,k)-kapacc(1,i,j-2,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)) )
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[kk-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;

		  d3up = u(2,i,j+2,k)-3*u(2,i,j+1,k)+3*u(2,i,j,k)-
		     u(2,i,j-1,k);
		  d3um =u(2,i,j+1,k)-3*u(2,i,j,k)+3*u(2,i,j-1,k)-
		     u(2,i,j-2,k);
		  d3kp = kap(2,i,j+2,k)-3*kap(2,i,j+1,k)+
		     3*kap(2,i,j,k)-kap(2,i,j-1,k);
		  d3km = kap(2,i,j+1,k)-3*kap(2,i,j,k)+
		     3*kap(2,i,j-1,k)-kap(2,i,j-2,k);
		  d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
		     3*uacc(2,i,j,k)-uacc(2,i,j-1,k);
		  d3uam =uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
		     3*uacc(2,i,j-1,k)-uacc(2,i,j-2,k);
		  d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
		     3*kapacc(2,i,j,k)-kapacc(2,i,j-1,k);
		  d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
		     3*kapacc(2,i,j-1,k)-kapacc(2,i,j-2,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[kk-1];
		  glambda(i,j,k) = glambda(i,j,k) + pd;
		  gmu(i,j,k) = gmu(i,j,k) + 2*pd;

		  d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+3*u(3,i,j,k)-
                     u(3,i,j-1,k);
		  d3um =u(3,i,j+1,k)-3*u(3,i,j,k)+3*u(3,i,j-1,k)-
		     u(3,i,j-2,k);
		  d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
		     3*kap(3,i,j,k)-kap(3,i,j-1,k);
		  d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
		     3*kap(3,i,j-1,k)-kap(3,i,j-2,k);
		  d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
		     3*uacc(3,i,j,k)- uacc(3,i,j-1,k);
		  d3uam =uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
		     3*uacc(3,i,j-1,k)-uacc(3,i,j-2,k);
		  d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
		     3*kapacc(3,i,j,k)-kapacc(3,i,j-1,k);
		  d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
		     3*kapacc(3,i,j-1,k)-kapacc(3,i,j-2,k);
		  pd = (c6*( 0.5*( d3up*d3kp+d3um*d3km  +
             dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
                 + c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh[kk-1];
		  gmu(i,j,k) = gmu(i,j,k) + pd;
// z-direction
		  if( k<=nk-1 )
		  {
// All derivatives are needed
		  d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
                     u(1,i,j,k-1);
 	          d3um = u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
                     u(1,i,j,k-2);
                  d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
		     3*kap(1,i,j,k)-kap(1,i,j,k-1);
                  d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
		     3*kap(1,i,j,k-1)-kap(1,i,j,k-2);
                  d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
		     3*uacc(1,i,j,k)-uacc(1,i,j,k-1);
                  d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
		     3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2);
                  d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
		     3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1);
                  d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
		     3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2);

                  float_sw4 d3vp = u(2,i,j,k+2)-3*u(2,i,j,k+1)+3*u(2,i,j,k)-
                     u(2,i,j,k-1);
                  float_sw4 d3vm = u(2,i,j,k+1)-3*u(2,i,j,k)+3*u(2,i,j,k-1)-
                     u(2,i,j,k-2);
                  float_sw4 d3lp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
		     3*kap(2,i,j,k)-kap(2,i,j,k-1);
                  float_sw4 d3lm = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
		     3*kap(2,i,j,k-1)-kap(2,i,j,k-2);
                  float_sw4 d3vap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
		     3*uacc(2,i,j,k)-uacc(2,i,j,k-1);
                  float_sw4 d3vam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
		     3*uacc(2,i,j,k-1)-uacc(2,i,j,k-2);
                  float_sw4 d3lap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
		     3*kapacc(2,i,j,k)-kapacc(2,i,j,k-1);
                  float_sw4 d3lam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
		     3*kapacc(2,i,j,k-1)-kapacc(2,i,j,k-2);

                  float_sw4 d3wp = u(3,i,j,k+2)-3*u(3,i,j,k+1)+3*u(3,i,j,k)-
                     u(3,i,j,k-1);
                  float_sw4 d3wm = u(3,i,j,k+1)-3*u(3,i,j,k)+3*u(3,i,j,k-1)-
                     u(3,i,j,k-2);
                  float_sw4 d3mp = kap(3,i,j,k+2)-3*kap(3,i,j,k+1)+
		     3*kap(3,i,j,k)-kap(3,i,j,k-1);
                  float_sw4 d3mm = kap(3,i,j,k+1)-3*kap(3,i,j,k)+
		     3*kap(3,i,j,k-1)-kap(3,i,j,k-2);
                  float_sw4 d3wap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
		     3*uacc(3,i,j,k)-uacc(3,i,j,k-1);
                  float_sw4 d3wam =uacc(3,i,j,k+1)-3*uacc(3,i,j,k)+
		     3*uacc(3,i,j,k-1)-uacc(3,i,j,k-2);
                  float_sw4 d3map = kapacc(3,i,j,k+2)-3*kapacc(3,i,j,k+1)+
		     3*kapacc(3,i,j,k)-kapacc(3,i,j,k-1);
                  float_sw4 d3mam = kapacc(3,i,j,k+1)-3*kapacc(3,i,j,k)+
		     3*kapacc(3,i,j,k-1)-kapacc(3,i,j,k-2);

#define SQR(x) (x)*(x)
                  float_sw4 mucof =SQR(met(2,i,j,k))+SQR(met(3,i,j,k))+SQR(met(4,i,j,k));
                  float_sw4 mucofm=SQR(met(2,i,j,k-1))+SQR(met(3,i,j,k-1))+SQR(met(4,i,j,k-1));
// u-u
                  pd = ( ( 0.5*( w6p[kk-1]*d3up*d3kp+w6m[kk-1]*d3um*d3km  +
             dt2o12*(w6p[kk-1]*d3up*d3kap + w6m[kk-1]*d3um*d3kam+ 
                     w6p[kk-1]*d3uap*d3kp + w6m[kk-1]*d3uam*d3km)))
                 + w8[kk-1]*c8*( (d3up-d3um)*(d3kp-d3km) +
             dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
		       (d3up-d3um)*(d3kap-d3kam))) );
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(2,i,j,k)));
                  glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(2,i,j,k));
                  if( k==nk-3 )
		  {
                     pd = 0.5*( -al4*d3um*d3km+
				  dt2o12*(-al4*d3um*d3kam -al4*d3uam*d3km ) );
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  
			pd*(mucofm+SQR(met(2,i,j,k-1)));
                     glambda(i,j,k-1) = glambda(i,j,k-1) + 
			pd*SQR(met(2,i,j,k-1));
                  }
// u-v
                  pd = ( ( 0.5*( w6p[kk-1]*d3vp*d3kp+w6m[kk-1]*d3vm*d3km  +
             dt2o12*(w6p[kk-1]*d3vp*d3kap + w6m[kk-1]*d3vm*d3kam+ 
                     w6p[kk-1]*d3vap*d3kp + w6m[kk-1]*d3vam*d3km)))
                 + w8[kk-1]*c8*( (d3vp-d3vm)*(d3kp-d3km) +
             dt2o12*(  (d3vap-d3vam)*(d3kp-d3km)+
                                 (d3vp-d3vm)*(d3kap-d3kam))) )
		     *met(2,i,j,k)*met(3,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==nk-3 )
		  {
                     pd = 0.5*( -al4*d3vm*d3km+
                    dt2o12*(-al4*d3vm*d3kam -al4*d3vam*d3km ) )*
			met(2,i,j,k-1)*met(3,i,j,k-1);
                     glambda(i,j,k-1) = glambda(i,j,k-1) + pd;
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  pd;
                  }
// u-w 
                  pd = ( ( 0.5*( w6p[kk-1]*d3wp*d3kp+w6m[kk-1]*d3wm*d3km  +
             dt2o12*(w6p[kk-1]*d3wp*d3kap + w6m[kk-1]*d3wm*d3kam+ 
                     w6p[kk-1]*d3wap*d3kp + w6m[kk-1]*d3wam*d3km)))
                 + w8[kk-1]*c8*( (d3wp-d3wm)*(d3kp-d3km) +
             dt2o12*(  (d3wap-d3wam)*(d3kp-d3km)+
                                 (d3wp-d3wm)*(d3kap-d3kam))) )
		     *met(2,i,j,k)*met(4,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==nk-3 )
		  {
                     pd = 0.5*( -al4*d3wm*d3km+
                        dt2o12*(-al4*d3wm*d3kam -al4*d3wam*d3km ) )*
			           met(2,i,j,k-1)*met(4,i,j,k-1);
                     glambda(i,j,k-1) = glambda(i,j,k-1) + pd;
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  pd;
                  }
// v-u
                  pd = ( ( 0.5*( w6p[kk-1]*d3up*d3lp+w6m[kk-1]*d3um*d3lm  +
             dt2o12*(w6p[kk-1]*d3up*d3lap + w6m[kk-1]*d3um*d3lam+ 
                     w6p[kk-1]*d3uap*d3lp + w6m[kk-1]*d3uam*d3lm)))
                 + w8[kk-1]*c8*( (d3up-d3um)*(d3lp-d3lm) +
             dt2o12*(  (d3uap-d3uam)*(d3lp-d3lm)+
                                 (d3up-d3um)*(d3lap-d3lam))) )
		     *met(2,i,j,k)*met(3,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==nk-3 )
		  {
                     pd = 0.5*( -al4*d3um*d3lm+
                    dt2o12*(-al4*d3um*d3lam -al4*d3uam*d3lm ) )*
			met(2,i,j,k-1)*met(3,i,j,k-1);
                     glambda(i,j,k-1) = glambda(i,j,k-1) + pd;
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  pd;
                  }

// v-v
                  pd = ( ( 0.5*( w6p[kk-1]*d3vp*d3lp+w6m[kk-1]*d3vm*d3lm  +
             dt2o12*(w6p[kk-1]*d3vp*d3lap + w6m[kk-1]*d3vm*d3lam+ 
                     w6p[kk-1]*d3vap*d3lp + w6m[kk-1]*d3vam*d3lm)))
                 + w8[kk-1]*c8*( (d3vp-d3vm)*(d3lp-d3lm) +
             dt2o12*(  (d3vap-d3vam)*(d3lp-d3lm)+
		       (d3vp-d3vm)*(d3lap-d3lam))) );
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(3,i,j,k)));
                  glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(3,i,j,k));
                  if( k==nk-3 )
		  {
                     pd = 0.5*( -al4*d3vm*d3lm+
				  dt2o12*(-al4*d3vm*d3lam -al4*d3vam*d3lm ) );
                     glambda(i,j,k-1) = glambda(i,j,k-1) + 
			pd*SQR(met(3,i,j,k-1));
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  
			pd*(mucofm+SQR(met(3,i,j,k-1)));
                  }

// v-w
                  pd = ( ( 0.5*( w6p[kk-1]*d3wp*d3lp+w6m[kk-1]*d3wm*d3lm  +
             dt2o12*(w6p[kk-1]*d3wp*d3lap + w6m[kk-1]*d3wm*d3lam+ 
                     w6p[kk-1]*d3wap*d3lp + w6m[kk-1]*d3wam*d3lm)))
                + w8[kk-1]*c8*( (d3wp-d3wm)*(d3lp-d3lm) +
             dt2o12*(  (d3wap-d3wam)*(d3lp-d3lm)+
                                 (d3wp-d3wm)*(d3lap-d3lam))) )
		     *met(3,i,j,k)*met(4,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==nk-3 )
		  {
                     pd = 0.5*( -al4*d3wm*d3lm+
                    dt2o12*(-al4*d3wm*d3lam -al4*d3wam*d3lm ) )*
			met(3,i,j,k-1)*met(4,i,j,k-1);
                     glambda(i,j,k-1) = glambda(i,j,k-1) + pd;
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  pd;
                  }

// w-u
                  pd = ( ( 0.5*( w6p[kk-1]*d3up*d3mp+w6m[kk-1]*d3um*d3mm  +
             dt2o12*(w6p[kk-1]*d3up*d3map + w6m[kk-1]*d3um*d3mam+ 
                     w6p[kk-1]*d3uap*d3mp + w6m[kk-1]*d3uam*d3mm)))
                 + w8[kk-1]*c8*( (d3up-d3um)*(d3mp-d3mm) +
             dt2o12*(  (d3uap-d3uam)*(d3mp-d3mm)+
                                 (d3up-d3um)*(d3map-d3mam))) )
		     *met(2,i,j,k)*met(4,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k == nk-3 )
		  {
                     pd = 0.5*( -al4*d3um*d3mm+
                    dt2o12*(-al4*d3um*d3mam -al4*d3uam*d3mm ) )*
			met(2,i,j,k-1)*met(4,i,j,k-1);
                     glambda(i,j,k-1) = glambda(i,j,k-1) + pd;
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  pd;
                  }

// w-v
                  pd = ( ( 0.5*( w6p[kk-1]*d3vp*d3mp+w6m[kk-1]*d3vm*d3mm  +
             dt2o12*(w6p[kk-1]*d3vp*d3map + w6m[kk-1]*d3vm*d3mam+ 
                     w6p[kk-1]*d3vap*d3mp + w6m[kk-1]*d3vam*d3mm)))
                 + w8[kk-1]*c8*( (d3vp-d3vm)*(d3mp-d3mm) +
             dt2o12*(  (d3vap-d3vam)*(d3mp-d3mm)+
                                 (d3vp-d3vm)*(d3map-d3mam))) )
		     *met(3,i,j,k)*met(4,i,j,k);
                  gmu(i,j,k) = gmu(i,j,k) + pd;
                  glambda(i,j,k) = glambda(i,j,k) + pd;
                  if( k==nk-3 )
		  {
                     pd = 0.5*( -al4*d3vm*d3mm+
                    dt2o12*(-al4*d3vm*d3mam -al4*d3vam*d3mm ) )*
			met(3,i,j,k-1)*met(4,i,j,k-1);
                     glambda(i,j,k-1) = glambda(i,j,k-1) + pd;
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  pd;
                  }
// w-w
                  pd = ( ( 0.5*( w6p[kk-1]*d3wp*d3mp+w6m[kk-1]*d3wm*d3mm  +
             dt2o12*(w6p[kk-1]*d3wp*d3map + w6m[kk-1]*d3wm*d3mam+ 
                     w6p[kk-1]*d3wap*d3mp + w6m[kk-1]*d3wam*d3mm)))
                 + w8[kk-1]*c8*( (d3wp-d3wm)*(d3mp-d3mm) +
             dt2o12*(  (d3wap-d3wam)*(d3mp-d3mm)+
		       (d3wp-d3wm)*(d3map-d3mam))) );
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+SQR(met(4,i,j,k)));
                  glambda(i,j,k) = glambda(i,j,k) + pd*SQR(met(4,i,j,k));

		  if( k == nk-3 )
		  {
                     pd = 0.5*( -al4*d3wm*d3mm+
				  dt2o12*(-al4*d3wm*d3mam -al4*d3wam*d3mm ) );
                     glambda(i,j,k-1) = glambda(i,j,k-1) + 
			pd*SQR(met(4,i,j,k-1));
                     gmu(i,j,k-1)     = gmu(i,j,k-1)   +  
			pd*(mucofm+SQR(met(4,i,j,k-1)));
		  }
		  }
	       }
      }


}

