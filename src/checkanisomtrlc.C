#include "sw4.h"
#include "EW.h"
#include "F77_FUNC.h"
extern "C" {
void F77_FUNC(dspev,DSPEV)(char & JOBZ, char & UPLO, int & N, double *AP, double *W, double *Z, int & LDZ, double *WORK, int & INFO);
}
void EW::checkanisomtrl_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			    float_sw4* __restrict__ rho, float_sw4* __restrict__ c, float_sw4& a_rhomin, 
			    float_sw4& a_rhomax, float_sw4& a_eigmin, float_sw4& a_eigmax )
			 
{
   float_sw4 rhomin =  1.0e38;
   float_sw4 rhomax = -1.0e38;
   float_sw4 eigmin =  1.0e38;
   float_sw4 eigmax = -1.0e38;
   size_t npts = static_cast<size_t>(ilast-ifirst+1)*(jlast-jfirst+1)*(klast-kfirst+1);
#pragma omp parallel   
   {
#pragma omp for reduction(max:rhomax,eigmax) reduction(min:rhomin,eigmin)
      for( size_t ind = 0 ; ind < npts ; ind++ ) 
      {
	 int info=0, six=6, one=1;
	 char n='N', l='L';
	 double a[21], eig[6], work[18], z;
	 if( rho[ind] < rhomin )
	    rhomin = rho[ind];
	 if( rho[ind] > rhomax )
	    rhomax = rho[ind];
	 for( int m=0 ; m < 21 ; m++ )
	    a[m] = c[ind+m*npts];
	 F77_FUNC(dspev,DSPEV)(n, l, six, a, eig, &z, one, work, info );
	 if( info != 0 )
	    cout << "ERROR in check_anisotropic_materialc:"<<
	       " could not compute eigenvalues. info from DSPEV = "<<
	       info<< endl;
	 if( eig[0] < eigmin )
	    eigmin = eig[0];
	 if( eig[5] > eigmax )
	    eigmax = eig[5];
      }
   }
   a_rhomin = rhomin;
   a_rhomax = rhomax;
   a_eigmin = eigmin;
   a_eigmax = eigmax;
}
