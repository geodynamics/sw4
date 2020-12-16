#ifdef OPTTEST_MODE
#include "dummy-classes.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <mpi.h>
#else
#include "EW.h"
#include "Mopt.h"
#include "MaterialParameterization.h"
#include "compute_f.h"
#endif

#ifdef USE_HDF5
#include <sachdf5.h>
#endif

using namespace std;

//-----------------------------------------------------------------------
void wolfecondition( EW& simulation, vector<vector<Source*> >& GlobalSources,
		     vector<vector<TimeSeries*> >& GlobalTimeSeries, 
		     vector<vector<TimeSeries*> >& GlobalObservations,
		     int nspar, int nmpars, double* xs, int nm, double* xm, double* ps, double* pm,
		     double* xsnew, double* xmnew, double f, double* dfsnew, double* dfmnew,
		     double& lambda, double maxstep, double minlambda, double cglen, double alpha,
		     double initslope, int myRank, double& fnew, bool& maxtaken, int& retcode,
		     Mopt* mopt )

//-----------------------------------------------------------------------
//  Algorithm A6.3.1-mod from book by Dennis and Schnabel. This routine only
//  implements the mod-part of Algorithm A6.3.1. Call it from 'linesearch'.
//
//-----------------------------------------------------------------------   

{
   double beta=0.9, lambdaprev, fpprev, maxlambda;
   int ns = nspar+nmpars;
   compute_f_and_df( simulation, nspar, nmpars, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, fnew, dfsnew, dfmnew, myRank, mopt );
   double slope=0;
   if( nm >  0 )
   {
      double slopetmp=0;
      for( int i=0 ; i < nm ; i++ )
	 slopetmp += pm[i]*dfmnew[i];
      MPI_Allreduce( &slopetmp, &slope, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   }
   for( int i=0 ; i < ns ; i++ )
      slope += ps[i]*dfsnew[i];
   if( slope < beta*initslope )
   {
      if( lambda == 1 && cglen < maxstep )
      {
	 maxlambda = maxstep/cglen;
	 do
	 {
	    lambdaprev = lambda;
	    fpprev = fnew;
	    lambda = 2*lambda<maxlambda ? 2*lambda:maxlambda;
	    for( int i=0 ; i < nm ; i++ )
	       xmnew[i] = xm[i]+lambda*pm[i];
	    for( int i=0 ; i < ns ; i++ )
	       xsnew[i] = xs[i]+lambda*ps[i];
	    compute_f( simulation, nspar, nmpars, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries, 
		       GlobalObservations, fnew, mopt );
	    if( fnew <= f + alpha*lambda*initslope )
	    {
	       compute_f_and_df( simulation, nspar, nmpars, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries,
				 GlobalObservations, fnew, dfsnew, dfmnew, myRank, mopt );
               slope = 0;
	       if( nm > 0 )
	       {
		  double slopetmp = 0;
		  for( int i=0 ; i < nm ; i++ )
		     slopetmp += pm[i]*dfmnew[i];
		  MPI_Allreduce( &slopetmp, &slope, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	       }
	       for( int i=0 ; i < ns ; i++ )
		  slope += ps[i]*dfsnew[i];
	    }
	 }
	 while( (fnew <= f + alpha*lambda*initslope) && (slope < beta*initslope) && (lambda < maxlambda) );
      }
      if( (lambda<1) || ( lambda>1 && fnew > f + alpha*lambda*initslope ) )
      {
         double lambdalo = lambda < lambdaprev?lambda:lambdaprev;
	 double lambdadiff = fabs(lambda-lambdaprev);
         double flo, fhi;
	 if( lambda < lambdaprev )
	 {
            flo = fnew;
	    fhi = fpprev;
	 }
	 else
	 {
	    flo = fpprev;
	    fhi = fnew;
	 }
	 do
	 {
	    double lambdaincr = -slope*lambdadiff*lambdadiff/(2*(fhi-(flo+slope*lambdadiff)));
	    if( lambdaincr < 0.2*lambdadiff )
	       lambdaincr = 0.2*lambdadiff;
	    lambda = lambdalo + lambdaincr;
	    for( int i=0 ; i < nm ; i++ )
	       xmnew[i] = xm[i]+lambda*pm[i];
	    for( int i=0 ; i < ns ; i++ )
	       xsnew[i] = xs[i]+lambda*ps[i];
	    compute_f( simulation, nspar, nmpars, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries,
		        GlobalObservations, fnew, mopt );
	    if( fnew > f + alpha*lambda*initslope )
	    {
	       lambdadiff = lambdaincr;
	       fhi = fnew;
	    }
	    else
	    {
	       compute_f_and_df( simulation, nspar, nmpars, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries,
				 GlobalObservations, fnew, dfsnew, dfmnew, myRank, mopt );
               slope = 0;
	       if( nm > 0 )
	       {
		  double slopetmp = 0;
		  for( int i=0 ; i < nm ; i++ )
		     slopetmp += pm[i]*dfmnew[i];
		  MPI_Allreduce( &slopetmp, &slope, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	       }
	       for( int i=0 ; i < ns ; i++ )
		  slope += ps[i]*dfsnew[i];
	       if( slope < beta*initslope )
	       {
		  lambdalo = lambda;
		  lambdadiff = lambdadiff - lambdaincr;
		  flo = fnew;
	       }
	    }
	 }
	 while( (slope < beta*initslope) && (lambdadiff >= minlambda) );
	 if( slope < beta*initslope )
	 {
	    fnew = flo;
	    for( int i=0 ; i < nm ; i++ )
	       xmnew[i] = xm[i]+lambdalo*pm[i];
	    for( int i=0 ; i < ns ; i++ )
	       xsnew[i] = xs[i]+lambdalo*ps[i];
	 }
      }
   }
   retcode = 0;
   maxtaken = false;
   if( lambda*cglen > 0.99*maxstep )
      maxtaken = true;
}

//-----------------------------------------------------------------------
void linesearch( EW& simulation, vector<vector<Source*> >& GlobalSources,
		 vector<vector<TimeSeries*> >& GlobalTimeSeries, 
		 vector<vector<TimeSeries*> >& GlobalObservations,
		 int nspar, int nmpars, double* xs, int nmpard_global, int nmpard, 
		 double* xm, double f, double* dfs,
		 double* dfm, double* ps, double* pm, double cgstep, double maxstep, double steptol,
		 double* xsnew, double* xmnew, double& fnew, double* sfs, double* sfm,
		 int myRank, int& retcode, int& nstep_reductions, bool testing, double* dfsnew,
		 double* dfmnew, Mopt* mopt )

//-----------------------------------------------------------------------
// Line seach by backtracking for CG.
// Algorithm A6.3.1 in book by Dennis and Schnabel, adapted for CG.
// 
// CG update without scaling is x^{k+1} = x^k + alpha_k p^k
//
// Scaled update: x^{k+1} = x^k + lambda*lambda_{max}*p^k/|| D p^k ||
//  where D is the diagonal matrix of scale factors.
// lambda is a normalized step length, where lambda = 1 corresponds 
// to the maximum allowed step. 
//
//  The maximum step is restricted by lambda_{max}. In unscaled units 
//  the max step is the minimum of alpha_k and lambda_{max}/|| D p^k ||. 
//
//  The scaled update gives:
//   || Dx^{k+1} - Dx^k || = lambda*lambda_{max} 
//  Thus lambda_{max} is the maximum change of x in the scaled units (lambda=1).
//
//  The minimum step length is set to approximately lambda = steptol/lambda_{max}
//  so that 
//   || Dx^{k+1} - Dx^k || = steptol
//  in maximum norm, i.e., steptol is the smallest allowable change in scaled units.
//
// Input: simulation - Simulation object
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        x          - Vector of source parameters, current iterate
//        f          - Misfit at x.
//        df         - Gradient at x.
//        p          - Conjugate gradient vector = search direction
//        cgstep     - Original estimate of step length, alpha_k above.
//        maxstep    - Maximum allowed step lambda_{max} above. 
//        steptol    - Minimum allowed step, for the scaled variables.
//        sf         - Scale factors.
//        myRank     - ID of this processor.
//
// Output: xnew    - New iterate.
//         fnew    - Misfit at new iterate.
//         retcode - Return code, indicating 0 = success, >0 = failure
//
//-----------------------------------------------------------------------

{
   bool maxtaken = false;
   int ns;
   double alpha = 1e-6, cglen=0, ang=0;
   vector<double> lambdasave, fcnsave;

   ns = nspar + nmpars;
   for( int i=0 ; i < ns ; i++ )
      xsnew[i] = xs[i];
   for( int i=0 ; i < nmpard ; i++ )
      xmnew[i] = xm[i];

   retcode = 2;

   cglen = ang = 0;
   if( nmpard_global > 0 )
   {
      for( int i=0 ; i < nmpard ; i++ )
      {
	 cglen += pm[i]*pm[i]/(sfm[i]*sfm[i]);
	 ang   += pm[i]*dfm[i];
      }
      double tmpv[2], tmp[2];
      tmpv[0] = cglen;
      tmpv[1] = ang;
      MPI_Allreduce( tmpv, tmp, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      cglen = tmp[0];
      ang = tmp[1];
   }
   for( int i=0 ; i < ns ; i++ )
   {
      cglen += ps[i]*ps[i]/(sfs[i]*sfs[i]);
      ang   += ps[i]*dfs[i];
   }
   if( ang >= 0 )
   {
      if( myRank == 0 )
      {
	 cout << "LINESEARCH: Warning, direction is not a descent direction" << endl;
         cout << "   switching to gradient direction search" << endl;
      }

      cglen = 0;
      if( nmpard_global > 0 )
      {
	 double cglenloc = 0;
	 for( int i=0 ; i < nmpard ; i++ )
	 {
	    pm[i] = -dfm[i]/(sfm[i]*sfm[i]);
	    cglenloc += pm[i]*pm[i]/(sfm[i]*sfm[i]);
	 }
	 MPI_Allreduce(&cglenloc,&cglen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
      for( int i=0 ; i < ns ; i++ )
      {
	 ps[i] = -dfs[i]/(sfs[i]*sfs[i]);
	 cglen += ps[i]*ps[i]/(sfs[i]*sfs[i]);
      }
   }
   cglen = sqrt(cglen);

   // p is scaled such that lambda=1 corresponds to the 
   // cg step x^{k+1} = x^k + cgstep*p^k unless maxstep is set lower.
   // i.e., maxstep is the maximum allowable relative change in x, ||x^{k+1}-x^k||/typx < maxstep


   if( cglen*cgstep > maxstep )
   {
      for( int i=0 ; i < nmpard ; i++ )
	 pm[i] = pm[i]*(maxstep/cglen);
      for( int i=0 ; i < ns ; i++ )
	 ps[i] = ps[i]*(maxstep/cglen);
      cglen = maxstep;
   }
   else
   {
      for( int i=0; i < nmpard ; i++ )
	 pm[i] = cgstep*pm[i];
      for( int i=0; i < ns ; i++ )
	 ps[i] = cgstep*ps[i];
      cglen = cgstep*cglen;
   }
      // minlambda is scaled such that x^{k+1}-x^k = lambda*p^k with lambda=minlambda implies
      //            || x^{k+1}-x^k ||_\infty < minlambda (termination criterium)
   double initslope=0;
   if( nmpard_global > 0 )
   {
      double initslopetmp = 0;
      for( int i=0; i < nmpard ; i++ )
         initslopetmp += dfm[i]*pm[i];
      MPI_Allreduce(&initslopetmp, &initslope,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   }
   for( int i=0; i < ns ; i++ )
      initslope += dfs[i]*ps[i];

   double rellength = 0;
   if( nmpard_global > 0 )
   {
      for( int i=0; i < nmpard ; i++ )
      {
	 double rlocal;
         if( fabs(xm[i]) > sfm[i] )
	    rlocal = fabs(pm[i])/fabs(xm[i]);
	 else
	    rlocal = fabs(pm[i])/sfm[i];
	 if( rlocal > rellength )
	    rellength = rlocal;
      }
      double rellengthtmp=rellength;
      MPI_Allreduce(&rellengthtmp,&rellength,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   }
   for( int i=0; i < ns ; i++ )
   {
      double rlocal;
      if( fabs(xs[i]) > sfs[i] )
	 rlocal = fabs(ps[i])/fabs(xs[i]);
      else
	 rlocal = fabs(ps[i])/sfs[i];
      if( rlocal > rellength )
	 rellength = rlocal;
   }
   double minlambda = steptol/rellength;
      cout << "steptol = " << steptol << " rellength " << rellength << " ns = " << ns << " minlambda=" << minlambda << endl; 
   double lambda = 1, fnewprev, lambdaprev;
   nstep_reductions = 0;
   
   // Restrict step size to inside of domain (physical range of material parameters)
   int ok=0;
   int ntries = 0;
   while( !ok && ntries < 30 )
   {

       for( int i=0; i < ns ; i++ )
	     xsnew[i] = xs[i] + lambda*ps[i];

      for( int i=0; i < nmpard ; i++ )
	 xmnew[i] = xm[i] + lambda*pm[i];

      int ng = simulation.mNumberOfGrids;
      vector<Sarray> rho(ng), mu(ng), la(ng);

	  //checkMinMax(nmpars, &xs[nspar], "linesearch: xs");
	  //checkMinMax(nmpars, &ps[nspar], "linesearch: ps");
	  std::cout << "scaling lambda=" << lambda << std::endl;
	  //checkMinMax(nmpars, &xsnew[nspar], "linesearch: xsnew");
      mopt->m_mp->get_material( nmpard, xmnew, nmpars, &xsnew[nspar], rho, mu, la, 
	    mopt->get_vp_min(), mopt->get_vp_max(), mopt->get_vs_min(), mopt->get_vs_max(), mopt->get_wave_mode());
      int ret_code = simulation.check_material( rho, mu, la, ok );

	  
      if( !ok )
      {
         if (myRank == 0)
         {
            printf("Material is not valid at relative step length (lambda) = %e, error code = %d, reducing lambda...\n", 
                   lambda, ret_code);
         }
         
	 lambda *= 0.7;
      }
      
      ntries++;
   }
   if( myRank == 0 && lambda != 1 )
      cout << "After constraining the material, lambda (scalefactor for step length) = " << lambda << " number of reductions = " << ntries << endl;

   if( !ok )
   {
      // Error, could not constrain material
      retcode = 3;
      // Take a small step. AP: calling compute_f when the material is out of bounds will lead to MPI-Abort()
      for( int i=0 ; i < ns ; i++ )
	 xsnew[i] = xs[i]; /* + lambda*ps[i]; */
      for( int i=0 ; i < nmpard ; i++ )
	 xmnew[i] = xm[i]; /* + lambda*pm[i];*/

	  if( myRank == 0 ) std::cout << "retcode=3" << std::endl;
      compute_f( simulation, nspar, nmpars, xsnew, nmpard, xmnew, GlobalSources, GlobalTimeSeries,
		 GlobalObservations, fnew, mopt );
      return;
   }
   
   bool lambdaprevdefined = false;
   while( retcode == 2 )
   {
	  if( myRank == 0 ) std::cout << "retcode=2" << std::endl;

      for( int i=0; i < ns ; i++ )
	     xsnew[i] = xs[i] + lambda*ps[i];

      for( int i=0; i < nmpard ; i++ )
	     xmnew[i] = xm[i] + lambda*pm[i];

    //if( myRank == 0 ) {  // Wei debugging
	//  save_array_to_disk(ns, ps, "ps.bin"); // multi-component *ns
	//  save_array_to_disk(ns, xsnew, "xsnew.bin"); // multi-component *ns
    //}

	  MPI_Barrier(MPI_COMM_WORLD); // Wei added for debugging

    //if( myRank == 0 ) {  // Wei debugging
	//  save_array_to_disk(ns, ps, "ps.bin"); // multi-component *ns
	//  save_array_to_disk(ns, xsnew, "xsnew.bin"); // multi-component *ns
    //}

	  MPI_Barrier(MPI_COMM_WORLD); // Wei added for debugging

      compute_f( simulation, nspar, nmpars, xsnew, nmpard, xmnew, GlobalSources, GlobalTimeSeries,
		 GlobalObservations, fnew, mopt );

            if( myRank == 0 )
            	 cout << "Evaluate at lambda = " << lambda << " gives f = " << fnew << " (initslope = "
            	      << initslope << ") " << endl;

      lambdasave.push_back(lambda);
      fcnsave.push_back(fnew);

      if( fnew < f + alpha*lambda*initslope )
      {
         if( mopt->m_wolfe )
	    wolfecondition( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nspar, nmpars, xs, nmpard, xm,
	           ps, pm, xsnew, xmnew, f, dfsnew, dfmnew, lambda, maxstep, minlambda, cglen, alpha, initslope, myRank,
			 fnew, maxtaken, retcode, mopt );
	 else
	 {
	    // Found satisfactory step
	    retcode = 0;
	    if( lambda == 1 && (cglen > 0.99*maxstep) )
	       maxtaken = true;
	    return;
	 }
      }
      else if( lambda<minlambda || ( std::isnan(lambda) || std::isnan(fnew)) )
      {
	    // Could not find satisfactory step
	 ang = 0;
         if( nmpard_global > 0 )
	 {
	    for( int i=0 ; i < nmpard ; i++ )
	       ang += pm[i]*dfm[i];
            double angtmp=ang;
            MPI_Allreduce(&angtmp,&ang,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	 }
	 for( int i=0 ; i < ns ; i++ )
	    ang += ps[i]*dfs[i];
  
    //if( myRank == 0 ) {  // Wei debugging
	//  save_array_to_disk(ns, ps, "ps.bin"); // multi-component *ns
    //}


	 if( myRank == 0 )
	 {
	    cout << "LINESEARCH: no satsifactory step found \n";
	    cout << "cg-alpha = " << cgstep << endl;
	    cout << "scprod = " << ang << endl;
	       //	       cout << "search direction = " << endl;
	       //	       for( int i=0 ; i < n ; i++ )
	       //		  cout << " " << p[i] << endl;
	       //	       cout << "gradient direction = " << endl;
	       //	       for( int i=0 ; i < n ; i++ )
	       //		  cout << " " << df[i] << endl;
	    cout << "maxstep   = " << maxstep << endl;
	    cout << "minlambda = " << minlambda << endl;
	    cout << "initslope = " << initslope << endl;
	    cout << " f = " << f << " fnew = " << fnew <<  " lambda = " << lambda ;
	    cout << " fnew-(f+alpha*lambda*initslope) = " << fnew-(f+alpha*lambda*initslope) << endl;
	    cout << "lambda and f(lambda) tried = " << endl;
	    for( int i=0 ; i < lambdasave.size() ; i++ )
	       cout << lambdasave[i] << " " << fcnsave[i] << endl;
	       //	       cout << "x = " << endl;
	       //	       for( int i=0 ; i < n ; i++ )
	       //		  cout << " " << x[i] << endl;
	 }
	 retcode = 1;
	    // Take a full step
	    for( int i=0 ; i < ns ; i++ )
	       xsnew[i] = xs[i] + ps[i];

	    for( int i=0 ; i < nmpard ; i++ )
	       xmnew[i] = xm[i] + pm[i];

  
    //if( myRank == 0 ) {  // Wei debugging
	//  save_array_to_disk(ns, xsnew, "xsnew.bin"); // multi-component *ns
    //}

  
    //if( myRank == 0 ) {  // Wei debugging
	//  save_array_to_disk(ns, xsnew, "xsnew.bin"); // multi-component *ns
    //}

	 // Compute return value for fnew
	    //            if( !testing && (xnew[2] < 0 && p[2] != 0) )
	    //	    {
	    //               lambda = -x[2]/p[2];
	    //	       for( int i=0 ; i < n ; i++ )
	    //		  xnew[i] = x[i] + lambda*p[i];
	    //	    }

	 std::cout << "retcode=1" << std::endl;
	 compute_f( simulation, nspar, nmpars, xsnew, nmpard, xmnew, GlobalSources,
		    GlobalTimeSeries, GlobalObservations, fnew, mopt );
	 return;
      } // end if ( Could not find satisfactory step )      
      else
      {
	    // Reduce step size
         if( myRank == 0 )
	    cout << " step failed for lambda = " << lambda;

	 double ltemp;
         if( !lambdaprevdefined )
	    //	 if( lambda == 1 )
	 {
	    ltemp = -initslope/(2*(fnew-f-initslope));
            if( ltemp > 0.5*lambda )
	       ltemp = 0.5*lambda;
	 }
	 else
	 {
	    double r1 =(fnew-f-lambda*initslope)/(lambda-lambdaprev);
	    double r2 =(fnewprev-f-lambdaprev*initslope)/(lambda-lambdaprev);
	    double a = 1/(lambda*lambda)*r1 - 1/(lambdaprev*lambdaprev)*r2;
	    double b = -lambdaprev/(lambda*lambda)*r1 +
		                       lambda/(lambdaprev*lambdaprev)*r2;
	    double disc = b*b-3*a*initslope;
	    if( a == 0 )
	       ltemp = -initslope/(2*b);
	    else
	       ltemp = (-b+sqrt(disc))/(3*a);
	    if( ltemp> 0.5*lambda )
	       ltemp = 0.5*lambda;
	 }
	 lambdaprev = lambda;
         lambdaprevdefined = true;
	 fnewprev = fnew;
	 if( ltemp < 0.1*lambda )
	    lambda = 0.1*lambda;
	 else
	    lambda = ltemp;

         if( myRank == 0 )
	    cout << " ... trying new lambda = " << lambda << endl;

	 nstep_reductions++;
      } // end else (do reduce the step size)
      
   }

   cout << " Line search done" << endl;
}

//-----------------------------------------------------------------------
void lbfgs( EW& simulation, int nspar, int nmpars, double* xs, 
            int nmpard, double* xm, 
	    vector<vector<Source*> >& GlobalSources,
	    vector<vector<TimeSeries*> >& GlobalTimeSeries,
	    vector<vector<TimeSeries*> >& GlobalObservations,
	    int myRank, Mopt* mopt )

//            MaterialParameterization* mp, int maxit, double tolerance,
//	    bool dolinesearch, int m, int ihess, bool use_wolfe, bool mcheck,
//	    bool output_ts, vector<Image*>& images )

//-----------------------------------------------------------------------
// l-BFGS minimization of misfit function.
// 
// Input: simulation - simulation object
//        nspar      - Number of parameters kept global, ie., exist in each processor.
//        xs         - The global parameters.
//        sf         - Scale factors for xs
//        nmpar      - Number of parameters for distributed arrays, ie, each processor 
//                     holds only a part of the array. nmpar is the size of the array in this processor.
//        xm         - This processor's part of the distributed parameter vector, used as initial guess.
//        sfm        - Scale factors for xm.
//        GlobalSource - Object representing the unique source.
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        m           - Number of previous vectors to save.
//        myRank      - ID of this processor.
//
// Output: xs, xm - The minimum is returned.
//
//
//-----------------------------------------------------------------------
{
   int j, k, ns;
   bool converged = false;
   //   bool dolinesearch = true;
   //   bool fletcher_reeves=true;
   //   bool use_wolfe = false;
   //   int varcase=0, stepselection=0;

   // Do not allow source to rise closer than this to the surface:
   double rnorm, f;
   //   double dfs[11], ds[11], da[11], xa[11], dfps[11], dx[11]
   double* dfs, *ds, *da, *xa, *dfps, *dx;
   double* dfm;
   bool testing=false;
   int nreductions = 0;
   int verbose = 0;
   bool hscale = (mopt->m_ihess_guess == 1);
   int maxit = mopt->m_maxit, m=mopt->m_nbfgs_vectors;
   bool dolinesearch=mopt->m_dolinesearch;
   double tolerance = mopt->m_tolerance;
   double* sf =mopt->m_sfs;
   double* sfm=mopt->m_sfm;

   // used variables: maxit, tolerance, dolinesearch
   //   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
   //				dolinesearch, varcase, testing );


   ns = nspar + nmpars;
   if(myRank==0) cout << "nspar=" << nspar << " nmpars=" << nmpars << " ns=" << ns << endl;

   if( maxit == 0 )
      return;
   
   int nmpard_global=0;
   MPI_Allreduce( &nmpard, &nmpard_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

   FILE *fd;
   FILE *fdx;
   const string parfile = mopt->m_path + "parameters.bin";
   if( myRank == 0 )
   {
      const string convfile = mopt->m_path + "convergence.log";
      fd = fopen(convfile.c_str(),"w");
      fprintf(fd, "it  max-nrm-gradient  max-nrm-update   misfit       step-length-reductions\n");
      fflush(fd);  //Wei added

      const string parafile = mopt->m_path + "parameters.log";
      if( nspar > 0 )
	 fdx=fopen(parafile.c_str(),"w");
   }


   if( nmpard > 0 )
      dfm = new double[nmpard];
   if( ns > 0 )
   {
      dfs  = new double[ns];
      ds   = new double[ns];
      da   = new double[ns];
      xa   = new double[ns];
      dfps = new double[ns];
      dx   = new double[ns];
   }
   
   if( myRank == 0 ) {
      cout << "Begin L-BFGS iteration by evaluating initial misfit and gradient..." << endl;
      checkMinMax(ns, xs, "xs");
	  checkMinMax(nmpard, xm, "xm");
    }

   compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, f, dfs, dfm, myRank, mopt, 0 );
    

	if( myRank == 0 ) {
      cout << "initial misfit f=" << f << endl;
	  cout << "ns=" << ns << endl;
      checkMinMax(ns, dfs, "gradient to update model perturbation (most critical): dfs");
	  checkMinMax(nmpard, dfm, "dfm");
	  //save_array_to_disk(ns, dfs, "dfs.bin"); // multi-component *ns
    }

   if( mopt->m_output_ts )
   {
     for( int e=0 ; e < GlobalTimeSeries.size() ; e++ )
     {
#ifdef USE_HDF5
       // Tang: need to create a HDF5 file before writing
       if (GlobalTimeSeries[e].size() > 0 && GlobalTimeSeries[e][0]->getUseHDF5()) {
         for (int tsi = 0; tsi < GlobalTimeSeries[e].size(); tsi++) 
           GlobalTimeSeries[e][tsi]->resetHDF5file();
         if(myRank == 0) 
           createTimeSeriesHDF5File(GlobalTimeSeries[e], GlobalTimeSeries[e][0]->getNsteps(), GlobalTimeSeries[e][0]->getDt(), "_ini");
         MPI_Barrier(MPI_COMM_WORLD);
       }
#endif

       for( int m = 0; m < GlobalTimeSeries[e].size(); m++ )
		    GlobalTimeSeries[e][m]->writeFile( "_ini" );
     }
   }

   if( myRank == 0 )
   {
      cout << "Initial misfit= "  << f << endl;
      if( nspar > 0 )
      {
	 cout << " scaled source gradient = ";
	 for( int i=0 ; i < nspar ; i++ )
	 {
	    cout << dfs[i]*sf[i] << " ";
	    if( i==5 )
	       cout << endl << "      " ;
	 }
      }
   }
   //   return;
   rnorm = 0;
   if( nmpard_global > 0 )
   {
      for( int i=0 ; i < nmpard ; i++ )
         rnorm = rnorm > fabs(dfm[i])*sfm[i] ? rnorm : fabs(dfm[i])*sfm[i];
      double rnormtmp = rnorm;
      MPI_Allreduce(&rnormtmp, &rnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   }
   for( int i=0 ; i < ns ; i++ )
      rnorm = rnorm > fabs(dfs[i])*sf[i] ? rnorm : fabs(dfs[i])*sf[i];
   if( myRank == 0 )
      cout << "Max norm of scaled total gradient = " << rnorm << endl;

   //   cout << endl;
   //   fprintf(fd, "%i %15.7g %15.7g %15.7g %i\n", 0, rnorm, 0.0, f, 0 );
   //   fprintf(fdx, "%i %i %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g\n",
   //	   0,0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10] );

   // s and y stores the m vectors

   double* s  = NULL;
   double* y  = NULL;
   if( ns > 0 )
   {
      s  = new double[m*ns]; 
      y  = new double[m*ns]; 

   }
   double* sm = NULL;
   double* ym = NULL;
   if( nmpard > 0 )
   {
      sm = new double[m*nmpard]; 
      ym = new double[m*nmpard]; 
   }

   double* rho = new double[m];
   double* al  = new double[m];

   double* dftemp=NULL;
   if( ns > 0 )
      dftemp  = new double[ns];
   
   
   double* dfmtemp = NULL;
   double* dm = NULL;
   double* dam = NULL;
   double* xam = NULL;
   double* dfpm = NULL;
   double* dmsave=NULL;
   double* dssave= NULL;
   if( ns > 0 )
      dssave  = new double[ns];

   if( nmpard > 0 )
   {
      dfmtemp = new double[nmpard];
      dm   = new double[nmpard];
      dam  = new double[nmpard];
      xam  = new double[nmpard];
      dfpm = new double[nmpard];
      dmsave = new double[nmpard];
   }

   
   for( int i=0 ; i < ns ; i++ )
   {
      ds[i] = da[i] = xa[i] = dx[i] = dfps[i] =0;
   }
   // kf points to the first vector
   int kf = 0;
   int it = 1;

   while( it <= maxit && !converged )
   {
      // perform the two-loop recursion (Alg 7.4) to compute search direction d=-H*df
      int nv = it-1 < m ? it-1 : m ;
      if( nv == 0 )
      {
	 for( int i=0 ; i < ns ; i++ )
	    ds[i] = -sf[i]*sf[i]*dfs[i];
	 for( int i=0 ; i < nmpard ; i++ )
	    dm[i] = -sfm[i]*sfm[i]*dfm[i];
      }
      else
      {
	 for( int j=0 ; j < ns ; j++ )
	    dftemp[j] = dfs[j];
	 for( int j=0 ; j < nmpard ; j++ )
	    dfmtemp[j] = dfm[j];

         int vi = kf+nv-1;
	 if( vi > m-1 )
	    vi = vi-m;
	 for( int i=nv ; i >= 1 ; i-- )
	 {
	    double scprods[2]={0,0};
	    if( nmpard_global > 0 )
	    {
	       double scprodsloc[2]={0,0};
	       for( int j=0 ; j < nmpard ; j++ )
	       {
		  scprodsloc[0] += sm[j+nmpard*vi]*ym[j+nmpard*vi];
		  scprodsloc[1] += sm[j+nmpard*vi]*dfmtemp[j];
	       }
	       MPI_Allreduce( scprodsloc, scprods, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	    }
	    for( int j=0 ; j < ns ; j++ )
	    {
	       scprods[0] += s[j+ns*vi]*y[j+ns*vi];
	       scprods[1] += s[j+ns*vi]*dftemp[j];
	    }
	    rho[vi] = 1/scprods[0];
	    al[vi] = scprods[1]*rho[vi];
	    for( int j=0 ; j < ns ; j++ )
	       dftemp[j] -= al[vi]*y[j+ns*vi];
	    for( int j=0 ; j < nmpard ; j++ )
	       dfmtemp[j] -= al[vi]*ym[j+nmpard*vi];
	    vi--;
	    if( vi == -1 )
	       vi = nv-1;
	 }

	 // H0*df
         if( hscale || nv == 0 )
	 {
	    // Use scale factors for diagonal H0
	    for( int  i=0 ; i < ns ; i++ )
	       ds[i] = sf[i]*sf[i]*dftemp[i];
	    for( int  i=0 ; i < nmpard ; i++ )
	       dm[i] = sfm[i]*sfm[i]*dfmtemp[i];
	 }
	 else
	 {
	    // Use diagonal H0 given by formula (7.20) 
	    //	    double gamn = 0, gamd=0;
            vi = kf+nv-1;
	    if( vi > m-1 )
	       vi = vi-m;

            double gams[2]={0,0};
	    if( nmpard_global > 0 )
	    {
               double gamsloc[2]={0,0};
	       for( int i=0 ; i < nmpard ; i++ )
	       {
		  gamsloc[0] += sm[i+nmpard*vi]*ym[i+nmpard*vi];
		  gamsloc[1] += ym[i+nmpard*vi]*ym[i+nmpard*vi];
	       }
	       MPI_Allreduce( gamsloc, gams, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	    }
	    for( int i=0 ; i < ns ; i++ )
	    {
	       gams[0] += s[i+ns*vi]*y[i+ns*vi];
	       gams[1] += y[i+ns*vi]*y[i+ns*vi];
	    }
	    gams[0] = gams[0]/gams[1];

	    for( int  i=0 ; i < ns ; i++ )
	       ds[i] = gams[0]*dftemp[i];
	    for( int  i=0 ; i < nmpard ; i++ )
	       dm[i] = gams[0]*dfmtemp[i];
	 }
         vi = kf;
	 for( int i=1 ; i <= nv ; i++ )
	 {
	    double scprod = 0;
	    if( nmpard_global > 0 )
	    {
               double scprodloc=0;
	       for( int j=0 ; j < nmpard ; j++ )
		  scprodloc += ym[j+nmpard*vi]*dm[j];
	       MPI_Allreduce( &scprodloc, &scprod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	    }
	    for( int j=0 ; j < ns ; j++ )
	       scprod += y[j+ns*vi]*ds[j];

	    for( int j=0 ; j < nmpard ; j++ )
	       dm[j] += (al[vi]-rho[vi]*scprod)*sm[j+nmpard*vi];
	    for( int j=0 ; j < ns ; j++ )
	       ds[j] += (al[vi]-rho[vi]*scprod)*s[j+ns*vi];
            vi++;

	    if( vi > nv-1 )
	       vi = 0;
	 }
	 for( int i=0 ; i < ns ; i++ )
	    ds[i] = -ds[i];
	 for( int i=0 ; i < nmpard ; i++ )
	    dm[i] = -dm[i];
      }
      // Done with two-loop recursion (Alg 7.4) 

      // Next do line search, or update with constant step length
      double alpha = 1;
      if( dolinesearch )
      {
	 for( int i=0 ; i < ns ; i++ )
	    da[i] = ds[i];
	 for( int i=0 ; i < nmpard ; i++ )
	    dam[i] = dm[i];
	 int retcode;
         double fp;
	 if( myRank == 0 )
	    cout << "Line search...  iteration= " << it << endl;

	 linesearch( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations,
		     nspar, nmpars, xs, nmpard_global, nmpard, xm, f, dfs, dfm, da, dam,
		     fabs(alpha), 10.0, tolerance, xa, xam, fp, sf, sfm, myRank,
		     retcode, nreductions, testing, dfps, dfpm, mopt );

	if( myRank == 0 )  cout << "Line search done...  retcode= " << retcode << endl;

         if( retcode == 3 )
	 	{
            if( myRank == 0 )
	       cout << "ERROR exit, could not find a steplength that gives a valid material model" << endl;
	    	it = maxit+1; // should terminate the iteration
	 	}

	 if( myRank == 0 )
         {
	    printf(".. return code %d misfit changed from %e to %e\n",  retcode,  f, fp);
	    cout << " .. return code "  << retcode << " misfit changed from " << f << " to " << fp << endl;
		cout << "line search alpha=" << alpha << endl;
		checkMinMax(ns, ds, "ds");
         }
         
      }
      else   // no line-search
      {
		cout << " update model w/o linesearch" << endl;

		for( int i=0 ; i < ns ; i++ )
			xa[i] = xs[i] + alpha*ds[i];   // updated model
		for( int i=0 ; i < nmpard ; i++ )
			xam[i] = xm[i] + alpha*dm[i];   
      }



      // xa is now the new iterate
      // store x^{k+1}-x^k in d to save memory
      double dxnorm = 0;
      if( nmpard_global > 0 )
      {
			for( int i=0 ; i < nmpard ; i++ ) // distributed parameters
			{
				double locnorm = fabs(xm[i]-xam[i]);
		// why scaling?
				// if( fabs(xam[i])> sfm[i] )
				//    locnorm /= fabs(xam[i]);
				// else
				//    locnorm /= sfm[i];
				if( locnorm > dxnorm )
				dxnorm = locnorm;
				dmsave[i] = dm[i];
				dm[i] = xam[i] - xm[i];
				xm[i]  = xam[i];
			}
				double dxnormloc=dxnorm;
				MPI_Allreduce(&dxnormloc,&dxnorm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      } // end if nmpard > 0 (distributed parameters)      
      

      for( int i=0 ; i < ns ; i++ ) // src and shared material parameters
      {
			double locnorm = fabs(xs[i]-xa[i]);
		// why scaling?
			// if( fabs(xa[i])> sf[i] )
			//    locnorm /= fabs(xa[i]);
			// else
			//    locnorm /= sf[i];
			if( locnorm > dxnorm )
				dxnorm = locnorm;
				dssave[i] = ds[i];
			ds[i] = xa[i] - xs[i];
			xs[i]  = xa[i];
      }
      //      if( myRank == 0 )
      //	 for( int i=0 ; i < ns ; i++ )
      //	    cout << " i="  << i << " " << ds[i] << " " << xs[i] << endl;


    //if( myRank == 0 ) {
	//  save_array_to_disk(ns, xs, "xs.bin");
      //checkMinMax(ns, xa, "updated model");
	//}

      compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
			GlobalObservations, f, dfps, dfpm, myRank, mopt, it );


      // Check Wolfe condition:
      double sc[2]={0,0};
      if( nmpard_global > 0 )
      {
			double sctmp[2]={0,0};
		for( int i=0 ; i < nmpard ; i++ )
		{
			sctmp[0] += dmsave[i]*dfm[i];
			sctmp[1] += dmsave[i]*dfpm[i];
		}
			MPI_Allreduce(sctmp,sc,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
      for( int i=0 ; i < ns ; i++ )
      {
         sc[0] += dssave[i]*dfs[i];
         sc[1] += dssave[i]*dfps[i];
      }
      if( verbose >= 2 &&  myRank == 0 )
	   cout << "Wolfe condition " <<  sc[1] << " " << sc[0] << " quotient " << sc[1]/sc[0] << " should be >= beta " << endl;

      rnorm = 0;
      if( nmpard_global > 0 )
      {
         double rnormloc = 0;
	 for( int i=0 ; i < nmpard ; i++ )
	    if( fabs(dfpm[i])*sfm[i] > rnormloc )
	       rnormloc = fabs(dfpm[i])*sfm[i];
         MPI_Allreduce(&rnormloc,&rnorm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      }

      for( int i=0 ; i < ns ; i++ )
	 if( fabs(dfps[i])*sf[i] > rnorm )
	    rnorm = fabs(dfps[i])*sf[i];

// Save the time series after each sub-iteration.
//      for( int ts=0 ; ts < GlobalTimeSeries.size() ; ts++ )
//	 GlobalTimeSeries[ts]->writeFile();

// Save shared material parameters, for restart.
      mopt->m_mp->write_parameters(parfile.c_str(),nmpars,xs);
	

// Check that wave speeds do not become too high or too low.
//      simulation.material_correction( nmpard, xm );
//      mp->constrain_material( nmpard, xm, nmpars, &xs[nspar] );

      if( myRank == 0 )
      {
	 cout << "-----------------------------------------------------------------------" << endl;
	 cout << " it=" << it << " max-norm scaled gradient= " << rnorm << " max-norm model update= " << dxnorm << endl;
	 cout << " Misfit= "  << f << endl;

	 fprintf(fd, "%i %15.7e %15.7e %15.7e %i\n", it, rnorm, dxnorm, f, nreductions );
	 fflush(fd);

         if( nspar>0 )
	 {
	    cout << " new x = " ;
	    for( int i=0 ; i < nspar ; i++ )
	    {
	       cout << xs[i] << " ";
	       if( i==5 )
		  cout << endl << "      " ;
	    }
	    cout << endl;
	    cout << " scaled source gradient = " ;
	    for( int i=0 ; i < nspar ; i++ )
	    {
	       cout << dfps[i]*sf[i] << " ";
	       if( i==5 )
		  cout << endl << "      " ;
	    }
	    cout << endl;

	    fprintf(fdx, "%d ", it );
            for( int p=0 ; p < nspar ; p++ )
	       fprintf(fdx, "%15.7g ",xs[p] );
            fprintf(fdx,"\n");
	    fflush(fdx);
	 }

      }
      converged = rnorm < tolerance;

      // Update s and y vectors
      if( it > m )
      {
	 for( int i=0 ; i < ns ; i++ )
	 {
	    s[i+ns*kf] = ds[i];
	    y[i+ns*kf] = dfps[i]-dfs[i];
	 }
	 for( int i=0 ; i < nmpard ; i++ )
	 {
	    sm[i+nmpard*kf] = dm[i];
	    ym[i+nmpard*kf] = dfpm[i]-dfm[i];
	 }
	 kf++;
	 if( kf > m-1 )
	    kf = 0;
      }
      else
      {
	 for( int i=0 ; i < ns ; i++ )
	 {
	    s[i+ns*(it-1)] = ds[i];
	    y[i+ns*(it-1)] = dfps[i]-dfs[i];
	 }
	 for( int i=0 ; i < nmpard ; i++ )
	 {
	    sm[i+nmpard*(it-1)] = dm[i];
	    ym[i+nmpard*(it-1)] = dfpm[i]-dfm[i];
	 }
      }

      // Advance one iteration
      it++;
      for( int i=0 ; i < ns ; i++ )
	     dfs[i] = dfps[i];
      for( int i=0 ; i < nmpard ; i++ )
	    dfm[i] = dfpm[i];
   } // end while it<maxit && !converged
   
   	 // if(myRank==0) save_array_to_disk(ns, dfs, "dfs_end.bin"); // multi-component *ns

   if( myRank == 0 )
   {
      fclose(fd);
      if( nspar > 0 )
	 fclose(fdx);
   }
   delete[] al;
   delete[] rho;
   if( ns > 0 )
   {
      delete[] s;
      delete[] y;
      delete[] dftemp;
      delete[] dssave;
      delete[] dfs;
      delete[] ds;
      delete[] da;
      delete[] xa;
      delete[] dfps;
      delete[] dx;
   }
   if( nmpard > 0 )
   {
      delete[] dfm;
      delete[] sm;
      delete[] ym;
      delete[] dfmtemp;
      delete[] dm;
      delete[] dam;
      delete[] xam;
      delete[] dfpm;
      delete[] dmsave;
   }
}
