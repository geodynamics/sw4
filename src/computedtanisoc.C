#include "sw4.h"
#include "EW.h"
#include "F77_FUNC.h"
extern "C" {
void F77_FUNC(dspev,DSPEV)(char & JOBZ, char & UPLO, int & N, double *AP, double *W, double *Z, int & LDZ, double *WORK, int & INFO);
}

void EW::maxwave( float_sw4 c[21], float_sw4 rho, float_sw4& eigestimate )
{
// Compute the maximum of the sum of eigenvalues of the matrix 
//   \sum_{i,j} ki* kj* Aij
//   where \sum_{i} ki*ki = 1, scaled by the density.
// This is equal to (4*mu+lambda)/rho in the isotropic case 
   int info, three=3, one=1;
   char n='N', l='L';
   double eg[3], a[6], work[9], z;
   a[0] = c[0]+c[6]+c[11];
   a[1] = c[1]+c[8]+c[13];
   a[2] = c[2]+c[9]+c[14];
   a[3] = c[6]+c[15]+c[18];
   a[4] = c[7]+c[16]+c[19];
   a[5] = c[11]+c[18]+c[20];
   F77_FUNC(dspev,DSPEV)(n, l, three, a, eg, &z, one, work, info );
   if( info != 0 )
      cout << "ERROR in EW::maxwave:" <<
	 " could not compute eigenvalues. info from DSPEV = " <<
	 info << endl;
   float_sw4 mxeig = eg[0]>eg[1]?eg[0]:eg[1];
   if( eg[2] > mxeig )
      mxeig = eg[2];
   eigestimate = mxeig/rho;
   //   eigestimate = MAX(eg[0],eg[1],eg[2])/rho;
}

void EW::maxwavecurv( float_sw4 c[45], float_sw4 rho, float_sw4 jac, float_sw4& eigestimate )
{
 // Traces of matrices, in order xx,xy,xz,yy,yz,zz
   int info, three=3, one=1;
   char n='N', l='L';
   double eg[3], a[6], work[9], z;
   a[0] = c[0] +c[3] +c[5];
   a[1] = c[18]+c[22]+c[26];
   a[2] = c[27]+c[31]+c[35];
   a[3] = c[6] +c[9] +c[11];
   a[4] = c[36]+c[40]+c[44];
   a[5] = c[12]+c[15]+c[17];
   F77_FUNC(dspev,DSPEV)(n, l, three, a, eg, &z, one, work, info );
   if( info != 0 )
      cout << "ERROR in EW::maxwavecurv:" <<
	 " could not compute eigenvalues. info from DSPEV = " <<
	 info << endl;
   float_sw4 mxeig = eg[0]>eg[1]?eg[0]:eg[1];
   if( eg[2] > mxeig )
      mxeig = eg[2];
   eigestimate = mxeig/(jac*rho);
}

void EW::computedtaniso2_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			    float_sw4* __restrict__ rho, float_sw4* __restrict__ c, float_sw4 cfl,
			     float_sw4 dx, float_sw4& a_dtloc )
{
   size_t npts = static_cast<size_t>(ilast-ifirst+1)*(jlast-jfirst+1)*(klast-kfirst+1);
   float_sw4 dtloc=1.0e19;
   float_sw4 cfldx2 = cfl*cfl*dx*dx;
#pragma omp parallel   
   {
#pragma omp for reduction(min:dtloc)
      for( size_t ind = 0 ; ind < npts ; ind++ ) 
      {
	 float_sw4 cpt[21], eigestimate;
	 for( int m=0 ; m < 21 ; m++ )
	    cpt[m] = c[ind+m*npts];
	 maxwave( cpt, rho[ind], eigestimate );
	 float_sw4 dtgp2 = cfldx2/eigestimate;
	 if( dtgp2 < dtloc*dtloc )
	    dtloc = sqrt(dtgp2);
      }
   }
   a_dtloc = dtloc;
}

void EW::computedtaniso2curv_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
				 float_sw4* __restrict__ rho, float_sw4* __restrict__ c, float_sw4* jac,
			         float_sw4 cfl, float_sw4& a_dtloc )
{
   size_t npts = static_cast<size_t>(ilast-ifirst+1)*(jlast-jfirst+1)*(klast-kfirst+1);
   float_sw4 dtloc=1.0e19;
   float_sw4 cfl2 = cfl*cfl;
#pragma omp parallel   
   {
#pragma omp for reduction(min:dtloc)
      for( size_t ind = 0 ; ind < npts ; ind++ ) 
      {
	 float_sw4 cpt[45], eigestimate;
	 for( int m=0 ; m < 45 ; m++ )
	    cpt[m] = c[ind+m*npts];
	 maxwavecurv( cpt, rho[ind], jac[ind], eigestimate );
	 float_sw4 dtgp2 = cfl2/eigestimate;
	 if( dtgp2 < dtloc*dtloc )
	    dtloc = sqrt(dtgp2);
      }
   }
   a_dtloc = dtloc;
}

