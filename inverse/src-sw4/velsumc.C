#include "sw4.h"
#include "EW.h"

void EW::velsum_ci( int is, int ie, int js, int je, int ks, int ke,
		    int i1, int i2, int j1, int j2, int k1, int k2,
		    float_sw4* __restrict__ mu, float_sw4* __restrict__ lambda,
		    float_sw4* __restrict__ rho, float_sw4& cp, float_sw4& cs,
		    size_t& npts ) 
{
   const size_t ni    = ie-is+1;
   const size_t nij   = ni*(je-js+1);
   const size_t nijk  = nij*(ke-ks+1);
   const size_t base  = -(is+ni*js+nij*ks);
   cp = cs = 0;
   npts = static_cast<size_t>(i2-i1+1)*static_cast<size_t>(j2-j1+1)*(k2-k1+1);
#pragma omp parallel
#pragma omp for reduction(+:cp,cs)
   for( int k=k1; k<=k2; k++ )
      for( int j=j1; j<=j2; j++ )
	 for( int i=i1; i<=i2; i++ )
	 {
	    size_t ind = base + i + ni*j + nij*k;
	    cp += sqrt( (2*mu[ind]+lambda[ind])/rho[ind] );
	    cs += sqrt( mu[ind]/rho[ind] );
	 }
}
