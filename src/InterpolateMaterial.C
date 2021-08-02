#include "EW.h"
#ifndef SW4_NOOMP
#include <omp.h>
#endif
#include <cfloat>

//-----------------------------------------------------------------------
void interpolatemtrl( int nx, int  ny, int nz, float_sw4 xmin, float_sw4 ymin,
		      float_sw4 zmin, float_sw4 hx, float_sw4 hy, float_sw4 hz,
		      float_sw4* __restrict__ a_rho, float_sw4* __restrict__ a_mu,
		      float_sw4*  __restrict__ a_lambda, int ib,
		      int ie, int jb, int je, int kb, int ke, int ibact, int ieact,
		      int jbact, int jeact, int kbact, int keact,
		      float_sw4* __restrict__ a_rhogrid, float_sw4* __restrict__ a_mugrid,
		      float_sw4* __restrict__ a_lambdagrid, float_sw4 hf, float_sw4 zmingrid )
{

#define rhogrid(i,j,k)       a_rhogrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define mugrid(i,j,k)         a_mugrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define lambdagrid(i,j,k) a_lambdagrid[i-ib+nig*(j-jb)+nijg*(k-kb)]

#define rho(i,j,k)       a_rho[i+nx2*(j)+nxy2*(k)]
#define mu(i,j,k)         a_mu[i+nx2*(j)+nxy2*(k)]
#define lambda(i,j,k) a_lambda[i+nx2*(j)+nxy2*(k)]

   const int nig = ie-ib+1;
   const int nijg= nig*(je-jb+1);
   const int nx2 = nx+2;
   const int nxy2= nx2*(ny+2);
   
   float_sw4 wgz, wgy, wgx;   
#pragma omp parallel private(wgx,wgy,wgz)
   {   
#pragma omp for
   for( int k=kbact ; k <= keact ; k++ )
   {
      float_sw4 z = (k-1)*hf + zmingrid;
      int kc = (z-zmin)/hz;
      if( kc<0 )
      {
	 kc  = 0;
	 wgz = 0;
      }
      else if( kc>nz )
      {
	 kc  = nz;
	 wgz = 1;
      }
      else
	 wgz = (z-zmin)/hz-kc;
      for( int j=jbact ; j<= jeact ;j++ )
      {
	 float_sw4 y = (j-1)*hf;
	 int jc = (y-ymin)/hy;
	 if( jc<0 )
	 {
	    jc  = 0;
	    wgy = 0;
	 }
	 else if( jc>ny )
	 {
	    jc  = ny;
	    wgy = 1;
	 }
	 else
	    wgy = (y-ymin)/hy-jc;
	 for( int i=ibact; i <= ieact ; i++ )
	 {
	    float_sw4 x = (i-1)*hf;
            int ic = (x-xmin)/hx;
	    if( ic<0 )
	    {
	       ic  = 0;
	       wgx = 0;
	    }
	    else if( ic>nx )
	    {
	       ic  = nx;
	       wgx = 1;
	    }
	    else
	       wgx = (x-xmin)/hx-ic;

	    rhogrid(i,j,k) = rhogrid(i,j,k) +
           (1-wgz)*(
              (1-wgy)*( (1-wgx)*rho(ic,jc,  kc) + wgx*rho(ic+1,jc,  kc) ) +
	        wgy  *( (1-wgx)*rho(ic,jc+1,kc) + wgx*rho(ic+1,jc+1,kc) ) ) 
            + wgz*(
              (1-wgy)*( (1-wgx)*rho(ic,jc,  kc+1) + wgx*rho(ic+1,jc,  kc+1) ) +
	        wgy  *( (1-wgx)*rho(ic,jc+1,kc+1) + wgx*rho(ic+1,jc+1,kc+1) ) );
	    mugrid(i,j,k) = mugrid(i,j,k) +
           (1-wgz)*(
             (1-wgy)*( (1-wgx)*mu(ic,jc,  kc) + wgx*mu(ic+1,jc,  kc) ) +
                 wgy*( (1-wgx)*mu(ic,jc+1,kc) + wgx*mu(ic+1,jc+1,kc) ) ) 
           + wgz*(
             (1-wgy)*( (1-wgx)*mu(ic,jc,  kc+1) + wgx*mu(ic+1,jc,  kc+1) ) +
	         wgy*( (1-wgx)*mu(ic,jc+1,kc+1) + wgx*mu(ic+1,jc+1,kc+1) ) );
	    lambdagrid(i,j,k) = lambdagrid(i,j,k) + 
           (1-wgz)*(
        (1-wgy)*( (1-wgx)*lambda(ic,jc,  kc) + wgx*lambda(ic+1,jc,  kc) ) +
            wgy*( (1-wgx)*lambda(ic,jc+1,kc) + wgx*lambda(ic+1,jc+1,kc) ) ) 
           + wgz*(
        (1-wgy)*( (1-wgx)*lambda(ic,jc,  kc+1) + wgx*lambda(ic+1,jc,  kc+1) )+
            wgy*( (1-wgx)*lambda(ic,jc+1,kc+1) + wgx*lambda(ic+1,jc+1,kc+1) ) );
	 }
      }
   }
   }
#undef rhogrid
#undef mugrid
#undef lambdagrid
#undef rho
#undef mu
#undef lambda
}

//-----------------------------------------------------------------------
void interpolatemtrlc( int nx, int  ny, int nz, float_sw4 xmin, float_sw4 ymin,
		      float_sw4 zmin, float_sw4 hx, float_sw4 hy, float_sw4 hz,
		      float_sw4* __restrict__ a_rho, float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, int ib,
		      int ie, int jb, int je, int kb, int ke, int ibact, int ieact,
		      int jbact, int jeact, int kbact, int keact, float_sw4* __restrict__ a_rhogrid,
		      float_sw4* __restrict__ a_mugrid, float_sw4* __restrict__ a_lambdagrid, 
		      float_sw4 hf, float_sw4* __restrict__ a_zgrid )
{

#define zgrid(i,j,k)           a_zgrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define rhogrid(i,j,k)       a_rhogrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define mugrid(i,j,k)         a_mugrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define lambdagrid(i,j,k) a_lambdagrid[i-ib+nig*(j-jb)+nijg*(k-kb)]

#define rho(i,j,k)       a_rho[i+nx2*(j)+nxy2*(k)]
#define mu(i,j,k)         a_mu[i+nx2*(j)+nxy2*(k)]
#define lambda(i,j,k) a_lambda[i+nx2*(j)+nxy2*(k)]

   const int nig = ie-ib+1;
   const int nijg= nig*(je-jb+1);
   const int nx2 = nx+2;
   const int nxy2= nx2*(ny+2);
   
   float_sw4 wgz, wgy, wgx;
#pragma omp parallel private(wgx,wgy,wgz)
   {   
#pragma omp for
   for( int k=kbact ; k <= keact ; k++ )
      for( int j=jbact ; j<= jeact ;j++ )
      {
	 float_sw4 y = (j-1)*hf;
	 int jc = (y-ymin)/hy;
	 if( jc<0 )
	 {
	    jc  = 0;
	    wgy = 0;
	 }
	 else if( jc>ny )
	 {
	    jc  = ny;
	    wgy = 1;
	 }
	 else
	    wgy = (y-ymin)/hy-jc;
	 for( int i=ibact; i <= ieact ; i++ )
	 {
	    int kc = (zgrid(i,j,k)-zmin)/hz;
	    if( kc < 0 )
	    {
	       kc  = 0;
	       wgz = 0;
	    }
	    else if( kc>nz )
	    {
	       kc  = nz;
	       wgz = 1;
	    }
	    else
	       wgz = (zgrid(i,j,k)-zmin)/hz-kc;
	    float_sw4 x = (i-1)*hf;
            int ic = (x-xmin)/hx;
	    if( ic<0 )
	    {
	       ic  = 0;
	       wgx = 0;
	    }
	    else if( ic>nx )
	    {
	       ic  = nx;
	       wgx = 1;
	    }
	    else
	       wgx = (x-xmin)/hx-ic;

	    rhogrid(i,j,k) = rhogrid(i,j,k) +
           (1-wgz)*(
             (1-wgy)*( (1-wgx)*rho(ic,jc,kc) + wgx*rho(ic+1,jc,kc) )+
              wgy*( (1-wgx)*rho(ic,jc+1,kc) + wgx*rho(ic+1,jc+1,kc) ) ) 
           + wgz*(
         (1-wgy)*( (1-wgx)*rho(ic,jc,kc+1) + wgx*rho(ic+1,jc,kc+1) )+
	 wgy*( (1-wgx)*rho(ic,jc+1,kc+1) + wgx*rho(ic+1,jc+1,kc+1) ) );
	    mugrid(i,j,k) = mugrid(i,j,k) + 
           (1-wgz)*(
             (1-wgy)*( (1-wgx)*mu(ic,jc,kc) + wgx*mu(ic+1,jc,kc) )+
              wgy*( (1-wgx)*mu(ic,jc+1,kc) + wgx*mu(ic+1,jc+1,kc) ) ) 
           + wgz*(
        (1-wgy)*( (1-wgx)*mu(ic,jc,kc+1) + wgx*mu(ic+1,jc,kc+1) )+
	wgy*( (1-wgx)*mu(ic,jc+1,kc+1) + wgx*mu(ic+1,jc+1,kc+1) ) );
	    lambdagrid(i,j,k) = lambdagrid(i,j,k) + 
           (1-wgz)*(
        (1-wgy)*( (1-wgx)*lambda(ic,jc,kc) + wgx*lambda(ic+1,jc,kc) )+
        wgy*( (1-wgx)*lambda(ic,jc+1,kc) + wgx*lambda(ic+1,jc+1,kc) ) ) 
           + wgz*(
      (1-wgy)*( (1-wgx)*lambda(ic,jc,kc+1) + wgx*lambda(ic+1,jc,kc+1))+
      wgy*((1-wgx)*lambda(ic,jc+1,kc+1) + wgx*lambda(ic+1,jc+1,kc+1)));
	 }
      }
   }
#undef zgrid
#undef rhogrid
#undef mugrid
#undef lambdagrid
#undef rho
#undef mu
#undef lambda
}

//-----------------------------------------------------------------------
void gradients( int nx, int  ny, int nz, float_sw4 xmin, float_sw4 ymin,
		float_sw4 zmin, float_sw4 hx, float_sw4 hy, float_sw4 hz,
		float_sw4* __restrict__ a_gradrho, float_sw4* __restrict__ a_gradmu, float_sw4* __restrict__ a_gradlambda,
		int ib, int ie, int jb, int je, int kb, int ke, int ibact, int ieact,
		int jbact, int jeact, int kbact, int keact, float_sw4* __restrict__ a_gradrhogrid,
		float_sw4* __restrict__ a_gradmugrid, float_sw4* __restrict__ a_gradlambdagrid, 
		float_sw4 hf, float_sw4 zmingrid )
{
//* Chain rule of trilinear interpolation

#define gradrhogrid(i,j,k)       a_gradrhogrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define gradmugrid(i,j,k)         a_gradmugrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define gradlambdagrid(i,j,k) a_gradlambdagrid[i-ib+nig*(j-jb)+nijg*(k-kb)]

#define gradrho(i,j,k)       a_gradrho[i+nx2*(j)+nxy2*(k)]
#define gradmu(i,j,k)         a_gradmu[i+nx2*(j)+nxy2*(k)]
#define gradlambda(i,j,k) a_gradlambda[i+nx2*(j)+nxy2*(k)]

   const int nig = ie-ib+1;
   const int nijg= nig*(je-jb+1);
   const int nx2 = nx+2;
   const int nxy2 = nx2 *(ny+2);
   const int nxyz2= nxy2*(nz+2);

   int nth = 1;
#ifndef SW4_NOOMP
#pragma omp parallel
   if( omp_get_thread_num() == 0 )
      nth=omp_get_num_threads();
#endif

   float_sw4* _gradrho = new float_sw4[nxyz2*nth];
   float_sw4* _gradmu = new float_sw4[nxyz2*nth];
   float_sw4* _gradlambda = new float_sw4[nxyz2*nth];
#define gradrhol(i,j,k,t) _gradrho[i+nx2*(j)+nxy2*(k)+nxyz2*(t)]
#define gradmul(i,j,k,t) _gradmu[i+nx2*(j)+nxy2*(k)+nxyz2*(t)]
#define gradlambdal(i,j,k,t) _gradlambda[i+nx2*(j)+nxy2*(k)+nxyz2*(t)]

#pragma omp parallel for
   for( size_t ind=0 ; ind < nxyz2*nth ; ind++ )
      _gradrho[ind]=_gradmu[ind]=_gradlambda[ind]=0;

#pragma omp parallel for 
   for( int k=kbact ; k <= keact ; k++ )
   {
      float_sw4 z = (k-1)*hf + zmingrid;
      int kc = (z-zmin)/hz;
      float_sw4 wgx, wgy, wgz;
      int t=0;
#ifndef SW4_NOOMP
      t=omp_get_thread_num();
#endif
      if( kc<0 )
      {
	 kc  = 0;
	 wgz = 0;
      }
      else if( kc>nz )
      {
	 kc  = nz;
	 wgz = 1;
      }
      else
	 wgz = (z-zmin)/hz-kc;
      for( int j=jbact ; j<= jeact ;j++ )
      {
	 float_sw4 y = (j-1)*hf;
	 int jc = (y-ymin)/hy;
	 if( jc<0 )
	 {
	    jc  = 0;
	    wgy = 0;
	 }
	 else if( jc>ny )
	 {
	    jc  = ny;
	    wgy = 1;
	 }
	 else
	    wgy = (y-ymin)/hy-jc;
	 for( int i=ibact; i <= ieact ; i++ )
	 {
	    float_sw4 x = (i-1)*hf;
            int ic = (x-xmin)/hx;
	    if( ic<0 )
	    {
	       ic  = 0;
	       wgx = 0;
	    }
	    else if( ic>nx )
	    {
	       ic  = nx;
	       wgx = 1;
	    }
	    else
	       wgx = (x-xmin)/hx-ic;

	    //  for( int l3=0 ; l3 <= 1 ; l3++ )
	    //    for( int l2=0 ; l2 <= 1 ; l2++ )
	    //	    for( int l1=0 ; l1 <= 1 ; l1++ )
	    //	      float_sw4 wgh3 = ((1-l3)*(1-wgz)+l3*wgz)*((1-l2)*(1-wgy)+l2*wgy)*((1-l1)*(1-wgx)+l1*wgx);
	    //	      gradrho(ic+l1,jc+l2,kc+l3)    += wgh3*gradrhogrid(i,j,k);
	    //	      gradmu(ic+l1,jc+l2,kc+l3)     += wgh3*gradmugrid(i,j,k);
	    //	      gradlambda(ic+l1,jc+l2,kc+l3) += wgh3*gradmugrid(i,j,k);	     

	    gradrhol(ic,jc,kc,t)    += (1-wgz)*(1-wgy)*(1-wgx)*gradrhogrid(i,j,k);
	    gradrhol(ic+1,jc,kc,t)  += (1-wgz)*(1-wgy)*  wgx  *gradrhogrid(i,j,k);
	    gradrhol(ic,jc+1,kc,t)  += (1-wgz)*  wgy*  (1-wgx)*gradrhogrid(i,j,k);
	    gradrhol(ic+1,jc+1,kc,t)+= (1-wgz)*  wgy*    wgx  *gradrhogrid(i,j,k);

	    gradrhol(ic,jc,kc+1,t)    += wgz*(1-wgy)*(1-wgx)*gradrhogrid(i,j,k);
	    gradrhol(ic+1,jc,kc+1,t)  += wgz*(1-wgy)*wgx*gradrhogrid(i,j,k);
	    gradrhol(ic,jc+1,kc+1,t)  += wgz*wgy*(1-wgx)*gradrhogrid(i,j,k);
	    gradrhol(ic+1,jc+1,kc+1,t)+= wgz*wgy*wgx*gradrhogrid(i,j,k);

	    gradmul(ic,jc,kc,t)     += (1-wgz)*(1-wgy)*(1-wgx)*gradmugrid(i,j,k);
	    gradmul(ic+1,jc,kc,t)   += (1-wgz)*(1-wgy)*wgx*gradmugrid(i,j,k);
	    gradmul(ic,jc+1,kc,t)   += (1-wgz)*wgy*(1-wgx)*gradmugrid(i,j,k);
	    gradmul(ic+1,jc+1,kc,t) += (1-wgz)*wgy*wgx*gradmugrid(i,j,k);

	    gradmul(ic,jc,kc+1,t)    += wgz*(1-wgy)*(1-wgx)*gradmugrid(i,j,k);
	    gradmul(ic+1,jc,kc+1,t)  += wgz*(1-wgy)*wgx*gradmugrid(i,j,k);
	    gradmul(ic,jc+1,kc+1,t)  +=  wgz*wgy*(1-wgx)*gradmugrid(i,j,k);
	    gradmul(ic+1,jc+1,kc+1,t)+= wgz*wgy*wgx*gradmugrid(i,j,k);

	    gradlambdal(ic,jc,kc,t) += (1-wgz)*(1-wgy)*(1-wgx)*gradlambdagrid(i,j,k);
	    gradlambdal(ic+1,jc,kc,t) += (1-wgz)*(1-wgy)*wgx*gradlambdagrid(i,j,k);
	    gradlambdal(ic,jc+1,kc,t) += (1-wgz)*wgy*(1-wgx)*gradlambdagrid(i,j,k);
	    gradlambdal(ic+1,jc+1,kc,t) += (1-wgz)*wgy*wgx*gradlambdagrid(i,j,k);

	    gradlambdal(ic,jc,kc+1,t) += wgz*(1-wgy)*(1-wgx)*gradlambdagrid(i,j,k);
	    gradlambdal(ic+1,jc,kc+1,t) += wgz*(1-wgy)*wgx*gradlambdagrid(i,j,k);
	    gradlambdal(ic,jc+1,kc+1,t) += wgz*wgy*(1-wgx)*gradlambdagrid(i,j,k);
	    gradlambdal(ic+1,jc+1,kc+1,t) += wgz*wgy*wgx*gradlambdagrid(i,j,k);
	 }
      }
   }

#pragma omp parallel for
   for( int kc=0 ; kc <= nz+1 ; kc++ )
      for( int jc=0 ; jc <= ny+1 ; jc++ )
	 for( int ic=0 ; ic <= nx+1 ; ic++ )
	    for( int t=0 ; t < nth ; t++ )
	    {
	       gradrho(ic,jc,kc)    += gradrhol(ic,jc,kc,t);
	       gradmu(ic,jc,kc)     += gradmul(ic,jc,kc,t);
	       gradlambda(ic,jc,kc) += gradlambdal(ic,jc,kc,t);
	    }

#undef gradrhogrid
#undef gradmugrid
#undef gradlambdagrid
#undef gradrho
#undef gradmu
#undef gradlambda
#undef gradrhol
#undef gradmul
#undef gradlambdal
   delete[] _gradrho;
   delete[] _gradmu;
   delete[] _gradlambda;
}

//-----------------------------------------------------------------------
void gradientsc( int nx, int  ny, int nz, float_sw4 xmin, float_sw4 ymin,
		 float_sw4 zmin, float_sw4 hx, float_sw4 hy, float_sw4 hz,
		 float_sw4* __restrict__ a_gradrho, float_sw4* __restrict__ a_gradmu, float_sw4* __restrict__ a_gradlambda,
		 int ib, int ie, int jb, int je, int kb, int ke, int ibact, int ieact,
		 int jbact, int jeact, int kbact, int keact, float_sw4* __restrict__ a_gradrhogrid,
		 float_sw4* __restrict__ a_gradmugrid, float_sw4* __restrict__ a_gradlambdagrid, 
		 float_sw4 hf, float_sw4* __restrict__ a_zgrid )
{
#define zgrid(i,j,k)                   a_zgrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define gradrhogrid(i,j,k)       a_gradrhogrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define gradmugrid(i,j,k)         a_gradmugrid[i-ib+nig*(j-jb)+nijg*(k-kb)]
#define gradlambdagrid(i,j,k) a_gradlambdagrid[i-ib+nig*(j-jb)+nijg*(k-kb)]

#define gradrho(i,j,k)       a_gradrho[i+nx2*(j)+nxy2*(k)]
#define gradmu(i,j,k)         a_gradmu[i+nx2*(j)+nxy2*(k)]
#define gradlambda(i,j,k) a_gradlambda[i+nx2*(j)+nxy2*(k)]

   const int nig = ie-ib+1;
   const int nijg= nig*(je-jb+1);
   const int nx2 = nx+2;
   const int nxy2= nx2*(ny+2);
   
   float_sw4 wgz, wgy, wgx;
#pragma omp parallel private(wgx,wgy,wgz)
   {   
#pragma omp for
   for( int k=kbact ; k <= keact ; k++ )
      for( int j=jbact ; j <= jeact ;j++ )
      {
	 float_sw4 y = (j-1)*hf;
	 int jc = (y-ymin)/hy;
	 if( jc<0 )
	 {
	    jc  = 0;
	    wgy = 0;
	 }
	 else if( jc>ny )
	 {
	    jc  = ny;
	    wgy = 1;
	 }
	 else
	    wgy = (y-ymin)/hy-jc;
	 for( int i=ibact; i <= ieact ; i++ )
	 {
	    int kc = (zgrid(i,j,k)-zmin)/hz;
	    if( kc <0 )
	    {
	       kc  = 0;
	       wgz = 0;
	    }
	    else if( kc>nz )
	    {
	       kc  = nz;
	       wgz = 1;
	    }
	    else
	       wgz = (zgrid(i,j,k)-zmin)/hz-kc;
	    float_sw4 x = (i-1)*hf;
            int ic = (x-xmin)/hx;
	    if( ic<0 )
	    {
	       ic  = 0;
	       wgx = 0;
	    }
	    else if( ic>nx )
	    {
	       ic  = nx;
	       wgx = 1;
	    }
	    else
	       wgx = (x-xmin)/hx-ic;
	    gradrho(ic,jc,kc) = gradrho(ic,jc,kc) +
	       (1-wgz)*(1-wgy)*(1-wgx)*gradrhogrid(i,j,k);
	    gradrho(ic+1,jc,kc) = gradrho(ic+1,jc,kc) +
	       (1-wgz)*(1-wgy)*wgx*gradrhogrid(i,j,k);
	    gradrho(ic,jc+1,kc) = gradrho(ic,jc+1,kc) +
	       (1-wgz)*wgy*(1-wgx)*gradrhogrid(i,j,k);
	    gradrho(ic+1,jc+1,kc) = gradrho(ic+1,jc+1,kc) +
	       (1-wgz)*wgy*wgx*gradrhogrid(i,j,k);
	    gradrho(ic,jc,kc+1) = gradrho(ic,jc,kc+1) +
	       wgz*(1-wgy)*(1-wgx)*gradrhogrid(i,j,k);
	    gradrho(ic+1,jc,kc+1) = gradrho(ic+1,jc,kc+1) +
	       wgz*(1-wgy)*wgx*gradrhogrid(i,j,k);
	    gradrho(ic,jc+1,kc+1) = gradrho(ic,jc+1,kc+1) +
	       wgz*wgy*(1-wgx)*gradrhogrid(i,j,k);
	    gradrho(ic+1,jc+1,kc+1) = gradrho(ic+1,jc+1,kc+1) +
	       wgz*wgy*wgx*gradrhogrid(i,j,k);
	    gradmu(ic,jc,kc) = gradmu(ic,jc,kc) +
	       (1-wgz)*(1-wgy)*(1-wgx)*gradmugrid(i,j,k);
	    gradmu(ic+1,jc,kc) = gradmu(ic+1,jc,kc) +
	       (1-wgz)*(1-wgy)*wgx*gradmugrid(i,j,k);
	    gradmu(ic,jc+1,kc) = gradmu(ic,jc+1,kc) +
	       (1-wgz)*wgy*(1-wgx)*gradmugrid(i,j,k);
	    gradmu(ic+1,jc+1,kc) = gradmu(ic+1,jc+1,kc) +
	       (1-wgz)*wgy*wgx*gradmugrid(i,j,k);
	    gradmu(ic,jc,kc+1) = gradmu(ic,jc,kc+1) +
	       wgz*(1-wgy)*(1-wgx)*gradmugrid(i,j,k);
	    gradmu(ic+1,jc,kc+1) = gradmu(ic+1,jc,kc+1) +
	       wgz*(1-wgy)*wgx*gradmugrid(i,j,k);
	    gradmu(ic,jc+1,kc+1) = gradmu(ic,jc+1,kc+1) +
	       wgz*wgy*(1-wgx)*gradmugrid(i,j,k);
	    gradmu(ic+1,jc+1,kc+1) = gradmu(ic+1,jc+1,kc+1) +
	       wgz*wgy*wgx*gradmugrid(i,j,k);
	    gradlambda(ic,jc,kc) = gradlambda(ic,jc,kc) +
	       (1-wgz)*(1-wgy)*(1-wgx)*gradlambdagrid(i,j,k);
	    gradlambda(ic+1,jc,kc) = gradlambda(ic+1,jc,kc) +
	       (1-wgz)*(1-wgy)*wgx*gradlambdagrid(i,j,k);
	    gradlambda(ic,jc+1,kc) = gradlambda(ic,jc+1,kc) +
	       (1-wgz)*wgy*(1-wgx)*gradlambdagrid(i,j,k);
	    gradlambda(ic+1,jc+1,kc) = gradlambda(ic+1,jc+1,kc) +
	       (1-wgz)*wgy*wgx*gradlambdagrid(i,j,k);
	    gradlambda(ic,jc,kc+1) = gradlambda(ic,jc,kc+1) +
	       wgz*(1-wgy)*(1-wgx)*gradlambdagrid(i,j,k);
	    gradlambda(ic+1,jc,kc+1) = gradlambda(ic+1,jc,kc+1) +
	       wgz*(1-wgy)*wgx*gradlambdagrid(i,j,k);
	    gradlambda(ic,jc+1,kc+1) = gradlambda(ic,jc+1,kc+1) +
	       wgz*wgy*(1-wgx)*gradlambdagrid(i,j,k);
	    gradlambda(ic+1,jc+1,kc+1) = gradlambda(ic+1,jc+1,kc+1) +
	       wgz*wgy*wgx*gradlambdagrid(i,j,k);
	 }
      }
      }
#undef zgrid
#undef gradrhogrid
#undef gradmugrid
#undef gradlambdagrid
#undef gradrho
#undef gradmu
#undef gradlambda
}

//-----------------------------------------------------------------------
void EW::interpolate( int nx, int ny, int nz, double xmin, double ymin, double zmin, double hx,
		      double hy, double hz, Sarray& rho, Sarray& mu, Sarray& lambda,
		      int grid, Sarray& rhogrid, Sarray& mugrid, Sarray& lambdagrid, bool update )
//
// Computes the material on computational grid from a given material perturbation on a parameter grid.
//
// Computes: rhogrid := I(rho)+mRho[grid]  (similar for mu, lambda)
//           where I(rho) interpolates input perturbation `rho' onto the computational grid.
//
//   Input: nx, ny, nz       - Dimensions of parameter grid
//          xmin, ymin, zmin - Origin of parameter grid.
//          hx, hy, hz       - Spacing of parameter grid.
//          rho, mu, lambda  - Material perturbation on parameter grid
//
//   Output: rhogrid, mugrid, lambdagrid - Material (perturbation+base material) on computational grid: g=`grid'
//
{
   int ifirst=m_iStart[grid];
   int ilast=m_iEnd[grid];
   int jfirst=m_jStart[grid];
   int jlast=m_jEnd[grid];
   int kfirst=m_kStart[grid];
   int klast=m_kEnd[grid];
   int ifirstact=m_iStartAct[grid];
   int ilastact=m_iEndAct[grid];
   int jfirstact=m_jStartAct[grid];
   int jlastact=m_jEndAct[grid];
   int kfirstact=m_kStartAct[grid];
   int klastact=m_kEndAct[grid];

   if( update )
   {
 // Start with a copy of the reference material
      rhogrid.copy( mRho[grid] );
      mugrid.copy( mMu[grid] );
      lambdagrid.copy( mLambda[grid] );
   }
   else
   {
 // Start from zero
      rhogrid.define( mRho[grid] );
      mugrid.define( mMu[grid] );
      lambdagrid.define( mLambda[grid] );
      rhogrid.set_to_zero();
      mugrid.set_to_zero();
      lambdagrid.set_to_zero();
   }
   
   double* rhop = rho.c_ptr();
   double* mup = mu.c_ptr();
   double* lambdap = lambda.c_ptr();
   double* rhogp = rhogrid.c_ptr();
   double* mugp = mugrid.c_ptr();
   double* lambdagp = lambdagrid.c_ptr();
   if( topographyExists() && grid >= mNumberOfCartesianGrids )
   {
      interpolatemtrlc(nx, ny, nz, xmin, ymin, zmin, 
		       hx, hy, hz, rhop, mup, lambdap,
		       ifirst, ilast, jfirst, jlast, kfirst, klast,
		       ifirstact, ilastact, jfirstact, jlastact, kfirstact, klastact, 
		       rhogp, mugp, lambdagp, mGridSize[grid], mZ[grid].c_ptr() );
   }
   else
   {
      interpolatemtrl(nx, ny, nz, xmin, ymin, zmin, 
		      hx, hy, hz, rhop, mup, lambdap,
		      ifirst, ilast, jfirst, jlast, kfirst, klast,
		      ifirstact, ilastact, jfirstact, jlastact, kfirstact, klastact, 
		      rhogp, mugp, lambdagp, mGridSize[grid], m_zmin[grid] );
   }
   communicate_array(rhogrid,grid);
   communicate_array(mugrid,grid);
   communicate_array(lambdagrid,grid);
}

//-----------------------------------------------------------------------
void EW::interpolate_to_coarse( int nx, int ny, int nz, double xmin, double ymin,
				double zmin, double hx, double hy, double hz,
				Sarray& rho, Sarray& mu, Sarray& lambda,
				vector<Sarray>& rhogrid, vector<Sarray>& mugrid,
				vector<Sarray>& lambdagrid, bool update )
{
// Compute material perturbation on parameter grid from a material that is given on the computational grid.
//
// Computes:  rho := I(rhogrid-mRho)   (and similar for mu and lambda)
//            where I() interpolates from computational grid onto the parameter grid.
//   
//   Input: nx, ny, nz       - Dimensions of parameter grid
//          xmin, ymin, zmin - Origin of parameter grid.
//          hx, hy, hz       - Spacing of parameter grid.
//          rhogrid, mugrid, lambdagrid - The material on the computational grids.
//
//   Output: rho, mu, lambda  - Material perturbation on parameter grid.
//

   int a1=1;
   if( !update )
      a1 = 0;

   int ig, jg, kg, g;
   rho.set_to_zero();
   mu.set_to_zero();
   lambda.set_to_zero();
   for( int k=1 ; k <= nz ; k++ )
      for( int j=1 ; j <= ny ; j++ )
	 for( int i=1 ; i <= nx ; i++ )
	 {
            double x = xmin + i*hx;
	    double y = ymin + j*hy;
	    double z = zmin + k*hz;
            computeNearestLowGridPoint( ig, jg, kg, g, x, y, z );
            if( interior_point_in_proc( ig, jg, g) )
	    {
	       double h = mGridSize[g];
	       double wghx = x/h-ig+1;
	       double wghy = y/h-jg+1;
	       double wghz = z/h-kg+1;
               rho(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*rhogrid[g](ig,jg,kg)+wghx*rhogrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*rhogrid[g](ig,jg+1,kg)+wghx*rhogrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*rhogrid[g](ig,jg,kg+1)+wghx*rhogrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*rhogrid[g](ig,jg+1,kg+1)+wghx*rhogrid[g](ig+1,jg+1,kg+1)) 
		  - a1*((1-wghy)*(1-wghz)*(
			           (1-wghx)*mRho[g](ig,jg,kg)+wghx*mRho[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mRho[g](ig,jg+1,kg)+wghx*mRho[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mRho[g](ig,jg,kg+1)+wghx*mRho[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mRho[g](ig,jg+1,kg+1)+wghx*mRho[g](ig+1,jg+1,kg+1))  );
               mu(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*mugrid[g](ig,jg,kg)+wghx*mugrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mugrid[g](ig,jg+1,kg)+wghx*mugrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mugrid[g](ig,jg,kg+1)+wghx*mugrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mugrid[g](ig,jg+1,kg+1)+wghx*mugrid[g](ig+1,jg+1,kg+1))
		  - a1*( (1-wghy)*(1-wghz)*(
			           (1-wghx)*mMu[g](ig,jg,kg)+wghx*mMu[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mMu[g](ig,jg+1,kg)+wghx*mMu[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mMu[g](ig,jg,kg+1)+wghx*mMu[g](ig+1,jg,kg+1))+
		      (wghy)*(wghz)*(  (1-wghx)*mMu[g](ig,jg+1,kg+1)+wghx*mMu[g](ig+1,jg+1,kg+1)) );
               lambda(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*lambdagrid[g](ig,jg,kg)+wghx*lambdagrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*lambdagrid[g](ig,jg+1,kg)+wghx*lambdagrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*lambdagrid[g](ig,jg,kg+1)+wghx*lambdagrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*lambdagrid[g](ig,jg+1,kg+1)+wghx*lambdagrid[g](ig+1,jg+1,kg+1))
		  -a1*((1-wghy)*(1-wghz)*(
			           (1-wghx)*mLambda[g](ig,jg,kg)+wghx*mLambda[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mLambda[g](ig,jg+1,kg)+wghx*mLambda[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mLambda[g](ig,jg,kg+1)+wghx*mLambda[g](ig+1,jg,kg+1))+
		    (wghy)*(wghz)*(  (1-wghx)*mLambda[g](ig,jg+1,kg+1)+wghx*mLambda[g](ig+1,jg+1,kg+1)));
      // Could do trilinear intp through ig,ig+1,jg,jg+1,kg,kg+1 instead
	       //               rho(i,j,k)    =    rhogrid[g](ig,jg,kg)-mRho[g](ig,jg,kg);
	       //               mu(i,j,k)     =     mugrid[g](ig,jg,kg)-mMu[g](ig,jg,kg);
	       //               lambda(i,j,k) = lambdagrid[g](ig,jg,kg)-mLambda[g](ig,jg,kg);
	    }
	 }
   Sarray tmp;
   tmp.copy(rho);
   MPI_Allreduce(tmp.c_ptr(),rho.c_ptr(),rho.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   tmp.copy(mu);
   MPI_Allreduce(tmp.c_ptr(),mu.c_ptr(),mu.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   tmp.copy(lambda);
   MPI_Allreduce(tmp.c_ptr(),lambda.c_ptr(),lambda.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

//-----------------------------------------------------------------------
void EW::interpolate_to_coarse_vel( int nx, int ny, int nz, double xmin, double ymin,
				    double zmin, double hx, double hy, double hz,
				    Sarray& rho, Sarray& cs, Sarray& cp,
				    vector<Sarray>& rhogrid, vector<Sarray>& mugrid,
				    vector<Sarray>& lambdagrid )
{
// Compute material perturbation on parameter grid from a material that is given on the computational grid.
//
// Computes:  rho := I(rhogrid-mRho)   (and similar for mu and lambda)
//            cs  := I(csgrid-mCs)
//            cp  := I(cpgrid-mCp)
//            where I() interpolates from computational grid onto the parameter grid.
//            csgrid and cpgrid are obtained by transforming (rhogrid,mugrid,lambdagrid)
//            mCs and mCp are obtained by transforming (mRho, mMu, mLambda)
//   
//   Input: nx, ny, nz       - Dimensions of parameter grid
//          xmin, ymin, zmin - Origin of parameter grid.
//          hx, hy, hz       - Spacing of parameter grid.
//          rhogrid, mugrid, lambdagrid - The material (rho,mu,lambda) on the computational grids.
//
//   Output: rho, cs, cp  - Material velocity perturbation on parameter grid.
//

   int ig, jg, kg, g;
   rho.set_to_zero();
   cs.set_to_zero();
   cp.set_to_zero();
   vector<bool> done(mNumberOfGrids);
   for( int g=0 ; g < mNumberOfGrids ; g++ )
      done[g] = false;
   vector<Sarray> csdiff(mNumberOfGrids);
   vector<Sarray> cpdiff(mNumberOfGrids);
   for( int k=1 ; k <= nz ; k++ )
      for( int j=1 ; j <= ny ; j++ )
	 for( int i=1 ; i <= nx ; i++ )
	 {
            double x = xmin + i*hx;
	    double y = ymin + j*hy;
	    double z = zmin + k*hz;
            computeNearestLowGridPoint( ig, jg, kg, g, x, y, z );
            if( interior_point_in_proc( ig, jg, g) )
	    {
	       if( !done[g] )
	       {
		  csdiff[g].define(mugrid[g]);
		  cpdiff[g].define(mugrid[g]);
		  for( int k1=rhogrid[g].m_kb ; k1 <= rhogrid[g].m_ke ; k1++)
		     for( int j1=rhogrid[g].m_jb ; j1 <= rhogrid[g].m_je ; j1++)
			for( int i1=rhogrid[g].m_ib ; i1 <= rhogrid[g].m_ie ; i1++)
			{
			   csdiff[g](i1,j1,k1)=sqrt(mugrid[g](i1,j1,k1)/rhogrid[g](i1,j1,k1))-
			                       sqrt(mMu[g](i1,j1,k1)/mRho[g](i1,j1,k1));
			   cpdiff[g](i1,j1,k1)=
			     sqrt((2*mugrid[g](i1,j1,k1)+lambdagrid[g](i1,j1,k1))/rhogrid[g](i1,j1,k1)) -
		             sqrt((2*mMu[g](i1,j1,k1)   +mLambda[g](i1,j1,k1)   )/mRho[g](i1,j1,k1));
			}
		  done[g] = true;
	       }
	       double h = mGridSize[g];
	       double wghx = x/h-ig+1;
	       double wghy = y/h-jg+1;
	       double wghz = z/h-kg+1;
               rho(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*rhogrid[g](ig,jg,kg)+wghx*rhogrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*rhogrid[g](ig,jg+1,kg)+wghx*rhogrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*rhogrid[g](ig,jg,kg+1)+wghx*rhogrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*rhogrid[g](ig,jg+1,kg+1)+wghx*rhogrid[g](ig+1,jg+1,kg+1)) 
		  -      ((1-wghy)*(1-wghz)*(
			           (1-wghx)*mRho[g](ig,jg,kg)+wghx*mRho[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mRho[g](ig,jg+1,kg)+wghx*mRho[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mRho[g](ig,jg,kg+1)+wghx*mRho[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mRho[g](ig,jg+1,kg+1)+wghx*mRho[g](ig+1,jg+1,kg+1))  );
               cs(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*csdiff[g](ig,jg,kg)+wghx*csdiff[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*csdiff[g](ig,jg+1,kg)+wghx*csdiff[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*csdiff[g](ig,jg,kg+1)+wghx*csdiff[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*csdiff[g](ig,jg+1,kg+1)+wghx*csdiff[g](ig+1,jg+1,kg+1));
               cp(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*cpdiff[g](ig,jg,kg)+wghx*cpdiff[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*cpdiff[g](ig,jg+1,kg)+wghx*cpdiff[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*cpdiff[g](ig,jg,kg+1)+wghx*cpdiff[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*cpdiff[g](ig,jg+1,kg+1)+wghx*cpdiff[g](ig+1,jg+1,kg+1));
	    }
	 }
   Sarray tmp;
   tmp.copy(rho);
   MPI_Allreduce(tmp.c_ptr(),rho.c_ptr(),rho.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   tmp.copy(cs);
   MPI_Allreduce(tmp.c_ptr(),cs.c_ptr(),cs.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   tmp.copy(cp);
   MPI_Allreduce(tmp.c_ptr(),cp.c_ptr(),cp.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

//-----------------------------------------------------------------------
void EW::interpolate_base_to_coarse( int nx, int ny, int nz, double xmin, double ymin,
				     double zmin, double hx, double hy, double hz,
				     Sarray& rho, Sarray& mu, Sarray& lambda )
{
// Compute material perturbation on parameter grid from the background material.
//
// Computes:  rho := I(mRho)
//            mu   := I(mMu)
//            lambda  := I(mLambda)
//            where I() interpolates from computational grid onto the parameter grid.
//            csgrid and cpgrid are obtained by transforming (rhogrid,mugrid,lambdagrid)
//            mCs and mCp are obtained by transforming (mRho, mMu, mLambda)
//   
//   Input: nx, ny, nz       - Dimensions of parameter grid
//          xmin, ymin, zmin - Origin of parameter grid.
//          hx, hy, hz       - Spacing of parameter grid.
//
//   Output: rho, mu, lambda  - Material velocity perturbation on parameter grid.
//
   int ig, jg, kg, g;
   rho.set_to_zero();
   mu.set_to_zero();
   lambda.set_to_zero();
   for( int k=1 ; k <= nz ; k++ )
      for( int j=1 ; j <= ny ; j++ )
	 for( int i=1 ; i <= nx ; i++ )
	 {
            double x = xmin + i*hx;
	    double y = ymin + j*hy;
	    double z = zmin + k*hz;
            computeNearestLowGridPoint( ig, jg, kg, g, x, y, z );
            if( interior_point_in_proc( ig, jg, g) )
	    {
	       double h = mGridSize[g];
	       double wghx = x/h-ig+1;
	       double wghy = y/h-jg+1;
	       double wghz = z/h-kg+1;
               rho(i,j,k) = ((1-wghy)*(1-wghz)*(
			           (1-wghx)*mRho[g](ig,jg,kg)+wghx*mRho[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mRho[g](ig,jg+1,kg)+wghx*mRho[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mRho[g](ig,jg,kg+1)+wghx*mRho[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mRho[g](ig,jg+1,kg+1)+wghx*mRho[g](ig+1,jg+1,kg+1))  );
               mu(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*mMu[g](ig,jg,kg)+wghx*mMu[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mMu[g](ig,jg+1,kg)+wghx*mMu[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mMu[g](ig,jg,kg+1)+wghx*mMu[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mMu[g](ig,jg+1,kg+1)+wghx*mMu[g](ig+1,jg+1,kg+1));
               lambda(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*mLambda[g](ig,jg,kg)+wghx*mLambda[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mLambda[g](ig,jg+1,kg)+wghx*mLambda[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mLambda[g](ig,jg,kg+1)+wghx*mLambda[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mLambda[g](ig,jg+1,kg+1)+wghx*mLambda[g](ig+1,jg+1,kg+1));
	    }
	 }
   Sarray tmp;
   tmp.copy(rho);
   MPI_Allreduce(tmp.c_ptr(),rho.c_ptr(),rho.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   tmp.copy(mu);
   MPI_Allreduce(tmp.c_ptr(),mu.c_ptr(),mu.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   tmp.copy(lambda);
   MPI_Allreduce(tmp.c_ptr(),lambda.c_ptr(),lambda.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

//-----------------------------------------------------------------------
void EW::interpolate_base_to_coarse_vel( int nx, int ny, int nz, double xmin, double ymin,
					 double zmin, double hx, double hy, double hz,
					 Sarray& rho, Sarray& cs, Sarray& cp )
{
// Compute material perturbation on parameter grid from the background material.
//
// Computes:  rho := I(mRho)
//            cs  := I(mCs)
//            cp  := I(mCp)
//            where I() interpolates from computational grid onto the parameter grid.
//            csgrid and cpgrid are obtained by transforming (rhogrid,mugrid,lambdagrid)
//            mCs and mCp are obtained by transforming (mRho, mMu, mLambda)
//   
//   Input: nx, ny, nz       - Dimensions of parameter grid
//          xmin, ymin, zmin - Origin of parameter grid.
//          hx, hy, hz       - Spacing of parameter grid.
//
//   Output: rho, cs, cp  - Material velocity perturbation on parameter grid.
//
   int ig, jg, kg, g;
   rho.set_to_zero();
   cs.set_to_zero();
   cp.set_to_zero();
   vector<bool> done(mNumberOfGrids);
   for( int g=0 ; g < mNumberOfGrids ; g++ )
      done[g] = false;
   vector<Sarray> csgrid(mNumberOfGrids);
   vector<Sarray> cpgrid(mNumberOfGrids);
   for( int k=1 ; k <= nz ; k++ )
      for( int j=1 ; j <= ny ; j++ )
	 for( int i=1 ; i <= nx ; i++ )
	 {
            double x = xmin + i*hx;
	    double y = ymin + j*hy;
	    double z = zmin + k*hz;
            computeNearestLowGridPoint( ig, jg, kg, g, x, y, z );
            if( interior_point_in_proc( ig, jg, g) )
	    {
	       if( !done[g] )
	       {
		  csgrid[g].define(mMu[g]);
		  cpgrid[g].define(mMu[g]);
		  for( int k1=mMu[g].m_kb ; k1 <= mMu[g].m_ke ; k1++)
		     for( int j1=mMu[g].m_jb ; j1 <= mMu[g].m_je ; j1++)
			for( int i1=mMu[g].m_ib ; i1 <= mMu[g].m_ie ; i1++)
			{
			   csgrid[g](i1,j1,k1)= sqrt(mMu[g](i1,j1,k1)/mRho[g](i1,j1,k1));
			   cpgrid[g](i1,j1,k1)= sqrt((2*mMu[g](i1,j1,k1)+mLambda[g](i1,j1,k1))/mRho[g](i1,j1,k1));
			}
		  done[g] = true;
	       }
	       double h = mGridSize[g];
	       double wghx = x/h-ig+1;
	       double wghy = y/h-jg+1;
	       double wghz = z/h-kg+1;
               rho(i,j,k) = ((1-wghy)*(1-wghz)*(
			           (1-wghx)*mRho[g](ig,jg,kg)+wghx*mRho[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mRho[g](ig,jg+1,kg)+wghx*mRho[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mRho[g](ig,jg,kg+1)+wghx*mRho[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mRho[g](ig,jg+1,kg+1)+wghx*mRho[g](ig+1,jg+1,kg+1))  );
               cs(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*csgrid[g](ig,jg,kg)+wghx*csgrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*csgrid[g](ig,jg+1,kg)+wghx*csgrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*csgrid[g](ig,jg,kg+1)+wghx*csgrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*csgrid[g](ig,jg+1,kg+1)+wghx*csgrid[g](ig+1,jg+1,kg+1));
               cp(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*cpgrid[g](ig,jg,kg)+wghx*cpgrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*cpgrid[g](ig,jg+1,kg)+wghx*cpgrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*cpgrid[g](ig,jg,kg+1)+wghx*cpgrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*cpgrid[g](ig,jg+1,kg+1)+wghx*cpgrid[g](ig+1,jg+1,kg+1));
	    }
	 }
   Sarray tmp;
   tmp.copy(rho);
   MPI_Allreduce(tmp.c_ptr(),rho.c_ptr(),rho.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   tmp.copy(cs);
   MPI_Allreduce(tmp.c_ptr(),cs.c_ptr(),cs.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   tmp.copy(cp);
   MPI_Allreduce(tmp.c_ptr(),cp.c_ptr(),cp.npts(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

//-----------------------------------------------------------------------
void EW::interpolation_gradient( int nx, int ny, int nz, double xmin, double ymin, double zmin, double hx,
		      double hy, double hz, Sarray& gradrho, Sarray& gradmu, Sarray& gradlambda,
		      int grid, Sarray& gradrhogrid, Sarray& gradmugrid, Sarray& gradlambdagrid )
{
   int ifirst= m_iStart[grid];
   int ilast = m_iEnd[grid];
   int jfirst= m_jStart[grid];
   int jlast = m_jEnd[grid];
   int kfirst= m_kStart[grid];
   int klast = m_kEnd[grid];

   int ifirstact = m_iStartAct[grid];
   int ilastact  = m_iEndAct[grid];
   int jfirstact = m_jStartAct[grid];
   int jlastact  = m_jEndAct[grid];
   int kfirstact = m_kStartAct[grid];
   int klastact  = m_kEndAct[grid];
   if( ifirstact < m_iStartInt[grid] )
      ifirstact = m_iStartInt[grid];
   if( ilastact > m_iEndInt[grid] )
      ilastact = m_iEndInt[grid];
   if( jfirstact < m_jStartInt[grid] )
      jfirstact = m_jStartInt[grid];
   if( jlastact > m_jEndInt[grid] )
      jlastact = m_jEndInt[grid];
   
   double* grhop = gradrho.c_ptr();
   double* gmup = gradmu.c_ptr();
   double* glambdap = gradlambda.c_ptr();
   double* grhogp = gradrhogrid.c_ptr();
   double* gmugp = gradmugrid.c_ptr();
   double* glambdagp = gradlambdagrid.c_ptr();
   if( topographyExists() && grid >= mNumberOfCartesianGrids )
   {
      gradientsc(nx, ny, nz, xmin, ymin, zmin, 
		 hx, hy, hz, grhop, gmup, glambdap,
		 ifirst, ilast, jfirst, jlast, kfirst, klast,
		 ifirstact, ilastact, jfirstact, jlastact, kfirstact, klastact, 
		 grhogp, gmugp, glambdagp, mGridSize[grid], mZ[grid].c_ptr() );
   }
   else
   {
      gradients( nx, ny, nz, xmin, ymin, zmin, 
		 hx, hy, hz, grhop, gmup, glambdap,
		 ifirst, ilast, jfirst, jlast, kfirst, klast,
		 ifirstact, ilastact, jfirstact, jlastact, kfirstact, klastact, 
		 grhogp, gmugp, glambdagp, mGridSize[grid], m_zmin[grid] );
   }
   //   communicate_array(rhogrid);
   //   communicate_array(mugrid);
   //   communicate_array(lambdagrid);
}

//-----------------------------------------------------------------------
void EW::update_and_transform_material( int g, Sarray& rho, Sarray& mu, Sarray& lambda,
float_sw4 vp_min, float_sw4 vp_max, float_sw4 vs_min, float_sw4 vs_max, int wave_mode)
{
   // Input (rho,mu,lambda) on grid g contains (rho,cs,cp) update from base material.
   // This routine returns (rho,mu,lambda) computed by adding the input update variables
   // to the base material (rhoB,csB,cpB) and transforming back to Lame` parameters.
   float_sw4* rhop=rho.c_ptr();
   float_sw4* mup=mu.c_ptr();
   float_sw4* lap=lambda.c_ptr();
   float_sw4* m_rho=mRho[g].c_ptr();
   float_sw4* m_mu=mMu[g].c_ptr();
   float_sw4* m_la=mLambda[g].c_ptr();

   float_sw4 cs_min, cs_max;
      cs_min=FLT_MAX;
      cs_max=FLT_MIN;

   for( size_t ind=0 ; ind < rho.m_npts ; ind++ )
   {
      // Add base material to update in velocity variables
	  float_sw4 cs, cp;
	  if(wave_mode>0) { // S-wave or both
	     //cs = sqrt(m_mu[ind]/m_rho[ind])              + mup[ind];
		 // impose a max percentage of 15%
		 cs = sqrt(m_mu[ind]/m_rho[ind]);
		 if(fabs(mup[ind]/cs)>0.15) cs = cs + mup[ind]/fabs(mup[ind])*cs*0.15;
		 else cs = cs + mup[ind];
	  }
	  else
         cs = sqrt(m_mu[ind]/m_rho[ind]);     
	  
	  
	  if(vs_min>0 && cs<vs_min) cs=vs_min;  // global velocity constraints
	  if(vs_max>0 && cs>vs_max) cs=vs_max;

      if(wave_mode==0 || wave_mode==2)  // P-wave or both
         cp = sqrt((2*m_mu[ind]+m_la[ind])/m_rho[ind])+ lap[ind];
	  else
	     cp = sqrt((2*m_mu[ind]+m_la[ind])/m_rho[ind]);

	  if(vp_min>0 && cp<vp_min) cp=vp_min;
	  if(vp_max>0 && cp>vp_max) cp=vp_max;

      float_sw4 rho= m_rho[ind]                             + rhop[ind];

	    if(cs<cs_min) cs_min=cs;
        if(cs>cs_max) cs_max=cs;

      // return total material as Lame parameters
      rhop[ind]= rho;
      mup[ind] = cs*cs*rho;      // update vs if S mode
      lap[ind] = (cp*cp-2*cs*cs)*rho;
	 
   }
    cout<< ">>>>>>>>>>>>>>>> get_material: vs_min=" << vs_min << " vs_max=" << vs_max << " cs_min=" << cs_min << " cs_max=" << cs_max << endl;

}

//-----------------------------------------------------------------------
void EW::transform_gradient( Sarray& rho, Sarray& mu, Sarray& lambda, 
			     Sarray& grho, Sarray& gmu, Sarray& glambda )
{
   // (grho,gmu,glambda) is gradient on grid g with respect to (rho,mu,lambda) at the
   //  grid points.
   // This routine returns the gradient with respect to (rho,cs,cp) at the grid points.
   // NOTE: input (grho,gmu,glambda) are overwritten with the transformed gradient.

   float_sw4* rhop=rho.c_ptr();
   float_sw4* mup=mu.c_ptr();
   float_sw4* lap=lambda.c_ptr();
   float_sw4* grhop=grho.c_ptr();
   float_sw4* gmup=gmu.c_ptr();
   float_sw4* glambdap=glambda.c_ptr();
   for( size_t ind=0 ; ind < rho.m_npts ; ind++ )
   {
      // Chain rule tranformation
      float_sw4 cs=sqrt(mup[ind]/rhop[ind]);
      float_sw4 cp=sqrt((2*mup[ind]+lap[ind])/rhop[ind]);
      grhop[ind] = cs*cs*gmup[ind]+(cp*cp-2*cs*cs)*glambdap[ind]+grhop[ind];
      gmup[ind]     = 2*rhop[ind]*cs*gmup[ind]-4*rhop[ind]*cs*glambdap[ind];
      glambdap[ind] = 2*rhop[ind]*cp*glambdap[ind];
   }
}
