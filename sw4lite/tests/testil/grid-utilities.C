#include <cstdlib>
#include <cmath>
void generate_ghgrid( int ib, int ie, int jb, int je, int kb, int ke,
		      double h, double* x, double* y, double* z, double topo_zmax )
{
   // Default in sw4
   int grid_interpolation_order = 4;
   double zetaBreak = 0.95;

   // Hardcoded for domain of size 1:
   double GaussianLx =0.15, GaussianLy=0.10, GaussianXc=0.6, GaussianYc=0.4, GaussianAmp=0.1;
   double igx2 = 1/(GaussianLx*GaussianLx);
   double igy2 = 1/(GaussianLy*GaussianLy);

   size_t ni=ie-ib+1, nj=je-jb+1, nk=ke-kb+1;
   for (int k=kb; k <= ke ; k++ )
      for (int j=jb; j <= je ; j++ )
	 for (int i=ib; i <= ie ; i++ )
	 {
	    size_t ind = i-ib+ni*(j-jb)+ni*nj*(k-kb);
	    x[ind] = (i-1)*h;
	    y[ind] = (j-1)*h;
	    double ztopo = - GaussianAmp*exp(-(x[ind]-GaussianXc)*(x[ind]-GaussianXc)*igx2
			      -(y[ind]-GaussianYc)*(y[ind]-GaussianYc)*igy2 );
	    double izb = 1.0/(zetaBreak*(nk-1));
	    double sa  = (k-1)*izb;
	    double omsm = (1-sa);
	    for( int l=2 ; l <= grid_interpolation_order ; l++ )
	       omsm *= (1-sa);
	    if( sa >= 1 )
	       z[ind] = topo_zmax - (nk-k)*h;
	    else
	       z[ind] = topo_zmax - (nk-k)*h - omsm*(topo_zmax-(nk-1)*h-ztopo);
	 }
}

