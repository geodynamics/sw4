#include "EW.h"
#include "F77_FUNC.h"

extern "C" {
   void F77_FUNC(interpolatemtrl,INTERPOLATEMTRL)( int*, int*, int*, double*, double*, double*,
						double*, double*, double*, double*, double*, double*,
						int*, int*, int*, int*, int*, int*,
						int*, int*, int*, int*, int*, int*,
						double*, double*, double*, double*, double* );
   void F77_FUNC(interpolatemtrlc,INTERPOLATEMTRLC)( int*, int*, int*, double*, double*, double*,
						double*, double*, double*, double*, double*, double*,
						int*, int*, int*, int*, int*, int*,
						int*, int*, int*, int*, int*, int*,
						double*, double*, double*, double*, double* );
   void F77_FUNC(gradients,GRADIENTS)( int*, int*, int*, double*, double*, double*,
						double*, double*, double*, double*, double*, double*,
						int*, int*, int*, int*, int*, int*,
						int*, int*, int*, int*, int*, int*,
						double*, double*, double*, double*, double* );
   void F77_FUNC(gradientsc,GRADIENTSC)( int*, int*, int*, double*, double*, double*,
						double*, double*, double*, double*, double*, double*,
						int*, int*, int*, int*, int*, int*,
						int*, int*, int*, int*, int*, int*,
						double*, double*, double*, double*, double* );
}


//-----------------------------------------------------------------------
void EW::interpolate( int nx, int ny, int nz, double xmin, double ymin, double zmin, double hx,
		      double hy, double hz, Sarray& rho, Sarray& mu, Sarray& lambda,
		      int grid, Sarray& rhogrid, Sarray& mugrid, Sarray& lambdagrid )
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

 // Start with a copy of the reference material
   rhogrid.copy( mRho[grid] );
   mugrid.copy( mMu[grid] );
   lambdagrid.copy( mLambda[grid] );
   
   double* rhop = rho.c_ptr();
   double* mup = mu.c_ptr();
   double* lambdap = lambda.c_ptr();
   double* rhogp = rhogrid.c_ptr();
   double* mugp = mugrid.c_ptr();
   double* lambdagp = lambdagrid.c_ptr();
   if( topographyExists() && grid == mNumberOfGrids -1 )
   {
      F77_FUNC(interpolatemtrlc,INTERPOLATEMTRLC)(&nx, &ny, &nz, &xmin, &ymin, &zmin, 
						&hx, &hy, &hz, rhop, mup, lambdap,
						&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						&ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact, &klastact, 
						rhogp, mugp, lambdagp, &mGridSize[grid], mZ.c_ptr() );
   }
   else
   {
      F77_FUNC(interpolatemtrl,INTERPOLATEMTRL)(&nx, &ny, &nz, &xmin, &ymin, &zmin, 
						&hx, &hy, &hz, rhop, mup, lambdap,
						&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						&ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact, &klastact, 
						rhogp, mugp, lambdagp, &mGridSize[grid], &m_zmin[grid] );
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
				vector<Sarray>& lambdagrid )
{
// Simple and not so accurate...
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
		  - ((1-wghy)*(1-wghz)*(
			           (1-wghx)*mRho[g](ig,jg,kg)+wghx*mRho[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mRho[g](ig,jg+1,kg)+wghx*mRho[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mRho[g](ig,jg,kg+1)+wghx*mRho[g](ig+1,jg,kg+1))+
		     (wghy)*(wghz)*(  (1-wghx)*mRho[g](ig,jg+1,kg+1)+wghx*mRho[g](ig+1,jg+1,kg+1))  );
               mu(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*mugrid[g](ig,jg,kg)+wghx*mugrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mugrid[g](ig,jg+1,kg)+wghx*mugrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mugrid[g](ig,jg,kg+1)+wghx*mugrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mugrid[g](ig,jg+1,kg+1)+wghx*mugrid[g](ig+1,jg+1,kg+1))
		  - ( (1-wghy)*(1-wghz)*(
			           (1-wghx)*mMu[g](ig,jg,kg)+wghx*mMu[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mMu[g](ig,jg+1,kg)+wghx*mMu[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mMu[g](ig,jg,kg+1)+wghx*mMu[g](ig+1,jg,kg+1))+
		      (wghy)*(wghz)*(  (1-wghx)*mMu[g](ig,jg+1,kg+1)+wghx*mMu[g](ig+1,jg+1,kg+1)) );
               lambda(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*lambdagrid[g](ig,jg,kg)+wghx*lambdagrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*lambdagrid[g](ig,jg+1,kg)+wghx*lambdagrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*lambdagrid[g](ig,jg,kg+1)+wghx*lambdagrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*lambdagrid[g](ig,jg+1,kg+1)+wghx*lambdagrid[g](ig+1,jg+1,kg+1))
		  -((1-wghy)*(1-wghz)*(
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
   if( topographyExists() && grid == mNumberOfGrids -1 )
   {
      F77_FUNC(gradientsc,GRADIENTSC)(&nx, &ny, &nz, &xmin, &ymin, &zmin, 
						&hx, &hy, &hz, grhop, gmup, glambdap,
						&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						&ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact, &klastact, 
						grhogp, gmugp, glambdagp, &mGridSize[grid], mZ.c_ptr() );
   }
   else
   {
      F77_FUNC(gradients,GRADIENTS)( &nx, &ny, &nz, &xmin, &ymin, &zmin, 
				     &hx, &hy, &hz, grhop, gmup, glambdap,
				     &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
				     &ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact, &klastact, 
				     grhogp, gmugp, glambdagp, &mGridSize[grid], &m_zmin[grid] );
   }
   //   communicate_array(rhogrid);
   //   communicate_array(mugrid);
   //   communicate_array(lambdagrid);
}
