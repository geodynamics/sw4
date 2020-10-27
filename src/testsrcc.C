#include "sw4.h"
#include "EW.h"
void EW::testsrc_ci( float_sw4* __restrict__ f, int ib, int ie, int jb, int je, int kb, int ke,
		     int nk, int wind[6], float_sw4 zmin, float_sw4 h, int kx[3], 
		     int ky[3], int kz[3], float_sw4 mom[3] ) 
{
   const size_t ni = ie-ib+1;
   const size_t nij = ni*(je-jb+1);
   const size_t nijk = nij*(ke-kb+1);
   const size_t base  = -(ib+ni*jb+nij*kb);
   const float_sw4 h3 = h*h*h;
   float_sw4 wgh[4] ={ 17.0/48, 59.0/48, 43.0/48, 49.0/48 };
   //   float_sw4 mom1=mom[0],mom2=mom[1],mom3=mom[2];
   float_sw4 mom1=0, mom2=0, mom3=0;
#pragma omp parallel
#pragma omp for reduction(+:mom1,mom2,mom3)
   for( int k=wind[4]; k <= wind[5]; k++ )
   {
      float_sw4 z = zmin + (k-1)*h;
      for( int j=wind[2]; j <= wind[3]; j++ )
      {
	 float_sw4 y = (j-1)*h;
	 for( int i=wind[0]; i <= wind[1]; i++ )
	 {
	    float_sw4 normfact = h3;
            float_sw4 x = (i-1)*h;
	    if( k <= 4 )
	       normfact *= wgh[k-1];
	    if( k>= nk-3 && k <= nk )
	       normfact *= wgh[nk-k];
	    float_sw4 x1, x2, x3, y1, y2, y3, z1, z2, z3;
	    if( kx[0]== 0 )
	       x1 = 1;
	    else
	       x1 = pow(x,kx[0]);
	    if( kx[1]== 0 )
	       x2 = 1;
	    else
	       x2 = pow(x,kx[1]);
	    if( kx[2]== 0 )
	       x3 = 1;
	    else
	       x3 = pow(x,kx[2]);
	    if( ky[0]== 0 )
	       y1 = 1;
	    else
	       y1 = pow(y,ky[0]);
	    if( ky[1]== 0 )
	       y2 = 1;
	    else
	       y2 = pow(y,ky[1]);
	    if( ky[2]== 0 )
	       y3 = 1;
	    else
	       y3 = pow(y,ky[2]);
	    if( kz[0]== 0 )
	       z1 = 1;
	    else
	       z1 = pow(z,kz[0]);
	    if( kz[1]== 0 )
	       z2 = 1;
	    else
	       z2 = pow(z,kz[1]);
	    if( kz[2]== 0 )
	       z3 = 1;
	    else
	       z3 = pow(z,kz[2]);
	    size_t ind = base+i+ni*j+nij*k;
	    mom1 += x1*y1*z1*f[ind]*normfact;
	    mom2 += x2*y2*z2*f[ind+nijk]*normfact;
	    mom3 += x3*y3*z3*f[ind+2*nijk]*normfact;
	 }
      }
   }
   mom[0] += mom1;
   mom[1] += mom2;
   mom[2] += mom3;
}

//-----------------------------------------------------------------------
void EW::testsrcc_ci( float_sw4* __restrict__ f, int ib, int ie, int jb, int je, 
                      int kb, int ke, int nk, int g, int wind[6], 
		      int kx[3], int ky[3], int kz[3], float_sw4 mom[3] )
{
   const size_t ni = ie-ib+1;
   const size_t nij = ni*(je-jb+1);
   const size_t nijk = nij*(ke-kb+1);
   const size_t base  = -(ib+ni*jb+nij*kb);
   //   const float_sw4 h3 = h*h*h;
   float_sw4 wgh[4] ={ 17.0/48, 59.0/48, 43.0/48, 49.0/48 };
   //   float_sw4 mom1=mom[0],mom2=mom[1],mom3=mom[2];
   float_sw4 mom1=0, mom2=0, mom3=0;
#pragma omp parallel
#pragma omp for reduction(+:mom1,mom2,mom3)
   for( int k=wind[4]; k <= wind[5]; k++ )
   {
      for( int j=wind[2]; j <= wind[3]; j++ )
      {
	 for( int i=wind[0]; i <= wind[1]; i++ )
	 {
            float_sw4 x = mX[g](i,j,k);
            float_sw4 y = mY[g](i,j,k);
            float_sw4 z = mZ[g](i,j,k);            
	    float_sw4 normfact = mJ[g](i,j,k);

	    if( k <= 4 )
	       normfact *= wgh[k-1];
	    if( k>= nk-3 && k <= nk )
	       normfact *= wgh[nk-k];
	    float_sw4 x1, x2, x3, y1, y2, y3, z1, z2, z3;
	    if( kx[0]== 0 )
	       x1 = 1;
	    else
	       x1 = pow(x,kx[0]);
	    if( kx[1]== 0 )
	       x2 = 1;
	    else
	       x2 = pow(x,kx[1]);
	    if( kx[2]== 0 )
	       x3 = 1;
	    else
	       x3 = pow(x,kx[2]);
	    if( ky[0]== 0 )
	       y1 = 1;
	    else
	       y1 = pow(y,ky[0]);
	    if( ky[1]== 0 )
	       y2 = 1;
	    else
	       y2 = pow(y,ky[1]);
	    if( ky[2]== 0 )
	       y3 = 1;
	    else
	       y3 = pow(y,ky[2]);
	    if( kz[0]== 0 )
	       z1 = 1;
	    else
	       z1 = pow(z,kz[0]);
	    if( kz[1]== 0 )
	       z2 = 1;
	    else
	       z2 = pow(z,kz[1]);
	    if( kz[2]== 0 )
	       z3 = 1;
	    else
	       z3 = pow(z,kz[2]);
	    size_t ind = base+i+ni*j+nij*k;
	    mom1 += x1*y1*z1*f[ind]*normfact;
	    mom2 += x2*y2*z2*f[ind+nijk]*normfact;
	    mom3 += x3*y3*z3*f[ind+2*nijk]*normfact;
            //            if( f[ind] != 0 )
            //            std::cout << " TS " << g << " (i,j,k)" << i << " " << j << " " << k << " " << x1*y1*z1*f[ind] << std::endl;
	 }
      }
   }
   mom[0] += mom1;
   mom[1] += mom2;
   mom[2] += mom3;
   //   std::cout << "mom 1 " << mom1 << std::endl;
}
