#include "EW.h"
using namespace std;

void compute_f( EW& simulation, int ns, double xs[11], int nm, double* xm,
		vector<Source*>& GlobalSources,
		vector<TimeSeries*>& GlobalTimeSeries,
		vector<TimeSeries*>& GlobalObservations,
		double& mf );

void compute_f_and_df( EW& simulation, int ns, double xs[11], int nmpar, double* xm,
		       vector<Source*>& GlobalSources,
		       vector<TimeSeries*>& GlobalTimeSeries,
		       vector<TimeSeries*>& GlobalObservations, 
		       double& f, double dfs[11], double* dfm, int myrank );

//-----------------------------------------------------------------------
void wolfecondition( EW& simulation, vector<Source*>& GlobalSources,
		     vector<TimeSeries*>& GlobalTimeSeries, vector<TimeSeries*>& GlobalObservations,
		     int ns, double xs[11], int nm, double* xm, double ps[11], double* pm,
		     double xsnew[11], double* xmnew, double f, double dfsnew[11], double* dfmnew,
		     double& lambda, double maxstep, double minlambda, double cglen, double alpha,
		     double initslope, int myRank, double& fnew, bool& maxtaken, int& retcode )

//-----------------------------------------------------------------------
//  Algorithm A6.3.1-mod from book by Dennis and Schnabel. This routine only
//  implements the mod-part of Algorithm A6.3.1. Call it from 'linesearch'.
//
//-----------------------------------------------------------------------   

{
   double beta=0.9, lambdaprev, fpprev, maxlambda;
   compute_f_and_df( simulation, ns, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, fnew, dfsnew, dfmnew, myRank );
   double slope=0, slopetmp=0;
   for( int i=0 ; i < nm ; i++ )
      slopetmp += pm[i]*dfmnew[i];
   MPI_Allreduce( &slopetmp, &slope, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
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
	    compute_f( simulation, ns, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries, GlobalObservations,
		       fnew );
	    if( fnew <= f + alpha*lambda*initslope )
	    {
	       compute_f_and_df( simulation, ns, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries,
				 GlobalObservations, fnew, dfsnew, dfmnew, myRank );
               slopetmp = 0;
	       for( int i=0 ; i < nm ; i++ )
		  slopetmp += pm[i]*dfmnew[i];
	       MPI_Allreduce( &slopetmp, &slope, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
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
	    compute_f( simulation, ns, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries, GlobalObservations,
		       fnew );
	    if( fnew > f + alpha*lambda*initslope )
	    {
	       lambdadiff = lambdaincr;
	       fhi = fnew;
	    }
	    else
	    {
	       compute_f_and_df( simulation, ns, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries,
				 GlobalObservations, fnew, dfsnew, dfmnew, myRank );
               slopetmp = 0;
	       for( int i=0 ; i < nm ; i++ )
		  slopetmp += pm[i]*dfmnew[i];
	       MPI_Allreduce( &slopetmp, &slope, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
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
void linesearch( EW& simulation, vector<Source*>& GlobalSources,
		 vector<TimeSeries*>& GlobalTimeSeries, vector<TimeSeries*>& GlobalObservations,
		 int ns, int nm, double xs[11], double* xm, double f, double dfs[11], double* dfm,
		 double ps[11], double* pm, double cgstep, double maxstep, double steptol,
		 double xsnew[11], double* xmnew, double& fnew, double sfs[11], double* sfm,
		 int myRank, int& retcode, int& nstep_reductions, bool testing, double dfsnew[11],
		 double* dfmnew )

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
//         retcode - Return code, indicating 0 - success, 1 - failure
//
//-----------------------------------------------------------------------

{
   bool maxtaken = false;
   double alpha = 1e-6, cglen=0, ang=0;
   vector<double> lambdasave, fcnsave;

   for( int i=0 ; i < 11 ; i++ )
      xsnew[i] = xs[i];
   for( int i=0 ; i < nm ; i++ )
      xmnew[i] = xm[i];

   retcode = 2;

   for( int i=0 ; i < nm ; i++ )
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
      double cglenloc = 0;
      for( int i=0 ; i < nm ; i++ )
      {
	 pm[i] = -dfm[i]/(sfm[i]*sfm[i]);
	 cglenloc += pm[i]*pm[i]/(sfm[i]*sfm[i]);
      }
      MPI_Allreduce(&cglenloc,&cglen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
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
      for( int i=0 ; i < nm ; i++ )
	 pm[i] = pm[i]*(maxstep/cglen);
      for( int i=0 ; i < ns ; i++ )
	 ps[i] = ps[i]*(maxstep/cglen);
      cglen = maxstep;
   }
   else
   {
      for( int i=0; i < nm ; i++ )
	 pm[i] = cgstep*pm[i];
      for( int i=0; i < ns ; i++ )
	 ps[i] = cgstep*ps[i];
      cglen = cgstep*cglen;
   }
      // minlambda is scaled such that x^{k+1}-x^k = lambda*p^k with lambda=minlambda implies
      //            || x^{k+1}-x^k ||_\infty < minlambda (termination criterium)
      double initslopetmp = 0;
      for( int i=0; i < nm ; i++ )
         initslopetmp += dfm[i]*pm[i];
      double initslope;
      MPI_Allreduce(&initslopetmp, &initslope,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      for( int i=0; i < ns ; i++ )
         initslope += dfs[i]*ps[i];

      double rellength = 0;
      for( int i=0; i < nm ; i++ )
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
      double lambda = 1, fnewprev, lambdaprev;
      nstep_reductions = 0;
      while( retcode == 2 )
      {
	 for( int i=0; i < ns ; i++ )
	    xsnew[i] = xs[i] + lambda*ps[i];
	 for( int i=0; i < nm ; i++ )
	    xmnew[i] = xm[i] + lambda*pm[i];

	 compute_f( simulation, ns, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries, GlobalObservations,
		    fnew );
         lambdasave.push_back(lambda);
         fcnsave.push_back(fnew);
         if( fnew < f + alpha*lambda*initslope )
	 {
	    // Found satisfactory step
            wolfecondition( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, ns, xs, nm, xm, ps, pm,
			    xsnew, xmnew, f, dfsnew, dfmnew, lambda, maxstep, minlambda, cglen, alpha, initslope, myRank,
			    fnew, maxtaken, retcode );
	       //	    retcode = 0;
	       //	    if( lambda == 1 && (cglen > 0.99*maxstep) )
	       //	       maxtaken = true;
	    return;
	 }
	 else if( lambda<minlambda )
	 {
	    // Could not find satisfactory step
	    ang = 0;
	    for( int i=0 ; i < nm ; i++ )
	       ang += pm[i]*dfm[i];
            double angtmp=ang;
            MPI_Allreduce(&angtmp,&ang,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    for( int i=0 ; i < ns ; i++ )
	       ang += ps[i]*dfs[i];

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
	    for( int i=0 ; i < nm ; i++ )
	       xmnew[i] = xm[i] + pm[i];

	    // Compute return value for fnew
	    //            if( !testing && (xnew[2] < 0 && p[2] != 0) )
	    //	    {
	    //               lambda = -x[2]/p[2];
	    //	       for( int i=0 ; i < n ; i++ )
	    //		  xnew[i] = x[i] + lambda*p[i];
	    //	    }
	    compute_f( simulation, ns, xsnew, nm, xmnew, GlobalSources, GlobalTimeSeries, GlobalObservations,
		       fnew );
	    return;
	 }
	 else
	 {
	    // Reduce step size
            double ltemp;
	    if( lambda == 1 )
	       ltemp = -initslope/(2*(fnew-f-initslope));
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
	    fnewprev = fnew;
	    if( ltemp < 0.1*lambda )
	       lambda = 0.1*lambda;
	    else
	       lambda = ltemp;
	    nstep_reductions++;
	 }
      }
}

//-----------------------------------------------------------------------
void lbfgs( EW& simulation, int ns, double xs[11], double sf[11], int nmpar,
	    double* xm, double* sfm,
	    vector<Source*>& GlobalSources,
	    vector<TimeSeries*>& GlobalTimeSeries,
	    vector<TimeSeries*> & GlobalObservations,
	    int m, int myRank )

//-----------------------------------------------------------------------
// l-BFGS minimization of misfit function.
// 
// Input: simulation - simulation object
//        ns         - Number of parameters kept global, ie., exist in each processor.
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
// Note: Configuration of the algorithm is done by a call 
//       to get_cgparameters of the simulation object.
//
//-----------------------------------------------------------------------
{
   //   const int n = 11;
   int nvar, j, k;
   bool done = false;
   bool dolinesearch = true;
   bool fletcher_reeves=true;
   int maxit, maxrestart, varcase=0, stepselection=0;

   // Do not allow source to rise closer than this to the surface:
   double tolerance;
   double f, dfs[11], ds[11], d2f[121], da[11], xa[11], dfps[11], dx[11], rnorm;
   double* dfm;
   bool testing=false, hscale=false;
   int nreductions = 0;
   
   // used variables: maxrestart, tolerance, dolinesearch
   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
				dolinesearch, varcase, testing );

   if( maxrestart == 0 )
      return;
   
   FILE *fd;
   FILE *fdx;
   if( myRank == 0 )
   {
      const string convfile = simulation.getOutputPath() + "convergence.log";
      fd = fopen(convfile.c_str(),"w");
      const string parafile = simulation.getOutputPath() + "parameters.log";
      fdx=fopen(parafile.c_str(),"w");
   }

   if( nmpar > 0 )
      dfm = new double[nmpar];
   
   if( myRank == 0 )
      cout << "Begin L-BFGS iteration by evaluating initial misfit and gradient..." << endl;

   compute_f_and_df( simulation, ns, xs, nmpar, xm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, f, dfs, dfm, myRank );
   if( myRank == 0 )
   {
      cout << "Initial misfit= "  << f << endl;
      if( ns > 0 )
      {
	 cout << " scaled source gradient = ";
	 for( int i=0 ; i < ns ; i++ )
	 {
	    cout << dfs[i]*sf[i] << " ";
	    if( i==5 )
	       cout << endl << "      " ;
	 }
      }
   }

   rnorm = 0;
   for( int i=0 ; i < nmpar ; i++ )
      rnorm = rnorm > fabs(dfm[i])*sfm[i] ? rnorm : fabs(dfm[i])*sfm[i];
   double rnormtmp = rnorm;
   MPI_Allreduce(&rnormtmp, &rnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
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
   if( nmpar > 0 )
   {
      sm = new double[m*nmpar]; 
      ym = new double[m*nmpar]; 
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
   double dssave[11];

   if( nmpar > 0 )
   {
      dfmtemp = new double[nmpar];
      dm   = new double[nmpar];
      dam  = new double[nmpar];
      xam  = new double[nmpar];
      dfpm = new double[nmpar];
      dmsave = new double[nmpar];
   }

   
   for( int i=0 ; i < 11 ; i++ )
   {
      ds[i] = da[i] = xa[i] = dx[i] = dfs[i]=dfps[i] =0;
   }
   // kf points to the first vector
   int kf = 0;
   int it = 1;

   while( it <= maxrestart && !done )
   {
      // perform the two-loop recursion (Alg 7.4) to compute search direction d=-H*df
      int nv = it-1 < m ? it-1 : m ;
      if( nv == 0 )
      {
	 for( int i=0 ; i < ns ; i++ )
	    ds[i] = -sf[i]*sf[i]*dfs[i];
	 for( int i=0 ; i < nmpar ; i++ )
	    dm[i] = -sfm[i]*sfm[i]*dfm[i];
      }
      else
      {
	 for( int j=0 ; j < ns ; j++ )
	    dftemp[j] = dfs[j];
	 for( int j=0 ; j < nmpar ; j++ )
	    dfmtemp[j] = dfm[j];
	 
         int vi = kf+nv-1;
	 if( vi > m-1 )
	    vi = vi-m;
	 for( int i=nv ; i >= 1 ; i-- )
	 {
            double scprodsloc[2]={0,0};
	    for( int j=0 ; j < nmpar ; j++ )
	    {
	       scprodsloc[0] += sm[j+nmpar*vi]*ym[j+nmpar*vi];
	       scprodsloc[1] += sm[j+nmpar*vi]*dfmtemp[j];
	    }
            double scprods[2];            
            MPI_Allreduce( scprodsloc, scprods, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	    for( int j=0 ; j < ns ; j++ )
	    {
	       scprods[0] += s[j+ns*vi]*y[j+ns*vi];
	       scprods[1] += s[j+ns*vi]*dftemp[j];
	    }
	    rho[vi] = 1/scprods[0];
	    al[vi] = scprods[1]*rho[vi];
	    for( int j=0 ; j < ns ; j++ )
	       dftemp[j] -= al[vi]*y[j+ns*vi];
	    for( int j=0 ; j < nmpar ; j++ )
	       dfmtemp[j] -= al[vi]*ym[j+nmpar*vi];
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
	    for( int  i=0 ; i < nmpar ; i++ )
	       dm[i] = sfm[i]*sfm[i]*dfmtemp[i];
	 }
	 else
	 {
	    // Use diagonal H0 given by formula (7.20) 
	    //	    double gamn = 0, gamd=0;
            vi = kf+nv-1;
	    if( vi > m-1 )
	       vi = vi-m;

            double gamsloc[2]={0,0};
	    for( int i=0 ; i < nmpar ; i++ )
	    {
	       gamsloc[0] += sm[i+nmpar*vi]*ym[i+nmpar*vi];
	       gamsloc[1] += ym[i+nmpar*vi]*ym[i+nmpar*vi];
	    }
            double gams[2];
            MPI_Allreduce( gamsloc, gams, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	    for( int i=0 ; i < ns ; i++ )
	    {
	       gams[0] += s[i+ns*vi]*y[i+ns*vi];
	       gams[1] += y[i+ns*vi]*y[i+ns*vi];
	    }
	    gams[0] = gams[0]/gams[1];

	    for( int  i=0 ; i < ns ; i++ )
	       ds[i] = gams[0]*dftemp[i];
	    for( int  i=0 ; i < nmpar ; i++ )
	       dm[i] = gams[0]*dfmtemp[i];
	 }
         vi = kf;
	 for( int i=1 ; i <= nv ; i++ )
	 {
            double scprodloc = 0;
	    for( int j=0 ; j < nmpar ; j++ )
	       scprodloc += ym[j+nmpar*vi]*dm[j];
            double scprod;
            MPI_Allreduce( &scprodloc, &scprod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	    for( int j=0 ; j < ns ; j++ )
	       scprod += y[j+ns*vi]*ds[j];

	    for( int j=0 ; j < nmpar ; j++ )
	       dm[j] += (al[vi]-rho[vi]*scprod)*sm[j+nmpar*vi];
	    for( int j=0 ; j < ns ; j++ )
	       ds[j] += (al[vi]-rho[vi]*scprod)*s[j+ns*vi];
            vi++;

	    if( vi > nv-1 )
	       vi = 0;
	 }
	 for( int i=0 ; i < ns ; i++ )
	    ds[i] = -ds[i];
	 for( int i=0 ; i < nmpar ; i++ )
	    dm[i] = -dm[i];
      }
   // Done with two-loop recursion (Alg 7.4) 

   // Next do line search, or update with constant step length
      double alpha = 1.0;
      if( dolinesearch )
      {
	 for( int i=0 ; i < ns ; i++ )
	    da[i] = ds[i];
	 for( int i=0 ; i < nmpar ; i++ )
	    dam[i] = dm[i];
	 int retcode;
         double fp;
	 if( myRank == 0 )
	    cout << "Line search.. " << endl;

	 linesearch( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations,
		     ns, nmpar, xs, xm, f, dfs, dfm, da, dam, fabs(alpha), 10.0, tolerance*0.01, xa, xam,
		     fp, sf, sfm, myRank,
		     retcode, nreductions, testing, dfps, dfpm );

	 if( myRank == 0 )
	    cout << " .. return code "  << retcode << " misfit changed from " << f << " to " << fp << endl;
      }
      else
      {
	 for( int i=0 ; i < ns ; i++ )
	    xa[i] = xs[i] + alpha*ds[i];
	 for( int i=0 ; i < nmpar ; i++ )
	    xam[i] = xm[i] + alpha*dm[i];
      }

   // xa is now the new iterate
   // store x^{k+1}-x^k in d to save memory
      double dxnorm = 0;
      for( int i=0 ; i < nmpar ; i++ )
      {
	 double locnorm = fabs(xm[i]-xam[i]);
	 if( fabs(xam[i])> sfm[i] )
	    locnorm /= fabs(xam[i]);
	 else
	    locnorm /= sfm[i];
	 if( locnorm > dxnorm )
	    dxnorm = locnorm;
         dmsave[i] = dm[i];
	 dm[i] = xam[i] - xm[i];
	 xm[i]  = xam[i];
      }
      double dxnormloc=dxnorm;
      MPI_Allreduce(&dxnormloc,&dxnorm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      for( int i=0 ; i < ns ; i++ )
      {
	 double locnorm = fabs(xs[i]-xa[i]);
	 if( fabs(xa[i])> sf[i] )
	    locnorm /= fabs(xa[i]);
	 else
	    locnorm /= sf[i];
	 if( locnorm > dxnorm )
	    dxnorm = locnorm;
         dssave[i] = ds[i];
	 ds[i] = xa[i] - xs[i];
	 xs[i]  = xa[i];
      }

      compute_f_and_df( simulation, ns, xs, nmpar, xm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, f, dfps, dfpm, myRank );
 // Check Wolfe condition:
      double sctmp[2]={0,0};
      double sc[2]={0,0};
      for( int i=0 ; i < nmpar ; i++ )
      {
         sctmp[0] += dmsave[i]*dfm[i];
	 sctmp[1] += dmsave[i]*dfpm[i];
      }
      MPI_Allreduce(sctmp,sc,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      for( int i=0 ; i < ns ; i++ )
      {
         sc[0] += dssave[i]*dfs[i];
         sc[1] += dssave[i]*dfps[i];
      }
      if( myRank == 0 )
	 cout << "Wolfe condition " <<  sc[1] << " " << sc[0] << " quotient " << sc[1]/sc[0] << " should be >= beta " << endl;

      double rnormloc = 0;
      for( int i=0 ; i < nmpar ; i++ )
	 if( fabs(dfpm[i])*sfm[i] > rnormloc )
	    rnormloc = fabs(dfpm[i])*sfm[i];
      MPI_Allreduce(&rnormloc,&rnorm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      for( int i=0 ; i < ns ; i++ )
	 if( fabs(dfps[i])*sf[i] > rnorm )
	    rnorm = fabs(dfps[i])*sf[i];

// Save the time series after each sub-iteration.
      for( int ts=0 ; ts < GlobalTimeSeries.size() ; ts++ )
	 GlobalTimeSeries[ts]->writeFile();

// Check that wave speeds do not become too high or too low.
      simulation.material_correction( nmpar, xm );

      if( myRank == 0 )
      {
	 cout << "-----------------------------------------------------------------------" << endl;
	 cout << " it=" << it << " dfnorm= " << rnorm << " dxnorm= " << dxnorm << endl;
         if( ns>0 )
	 {
	    cout << " new x = " ;
	    for( int i=0 ; i < ns ; i++ )
	    {
	       cout << xs[i] << " ";
	       if( i==5 )
		  cout << endl << "      " ;
	    }
	    cout << endl;
	 }
	 cout << " Misfit= "  << f << endl;
         if( ns > 0 )
	 {
	    cout << " scaled gradient = " ;
	    for( int i=0 ; i < ns ; i++ )
	    {
	       cout << dfps[i]*sf[i] << " ";
	       if( i==5 )
		  cout << endl << "      " ;
	    }
	    cout << endl;
	 }
	 fprintf(fd, "%i %15.7g %15.7g %15.7g %i\n", it, rnorm, dxnorm, f, nreductions );
         if( ns>0)
	 {
	    fprintf(fdx, "%i %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g\n",
		 it, xs[0], xs[1], xs[2], xs[3], xs[4], xs[5], xs[6], xs[7], xs[8], xs[9], xs[10] );
	    fflush(fdx);
	 }
	 fflush(fd);
      }
      done = rnorm < tolerance;

      // Update s and y vectors
      if( it > m )
      {
	 for( int i=0 ; i < ns ; i++ )
	 {
	    s[i+ns*kf] = ds[i];
	    y[i+ns*kf] = dfps[i]-dfs[i];
	 }
	 for( int i=0 ; i < nmpar ; i++ )
	 {
	    sm[i+nmpar*kf] = dm[i];
	    ym[i+nmpar*kf] = dfpm[i]-dfm[i];
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
	 for( int i=0 ; i < nmpar ; i++ )
	 {
	    sm[i+nmpar*(it-1)] = dm[i];
	    ym[i+nmpar*(it-1)] = dfpm[i]-dfm[i];
	 }
      }

    // Advance one iteration
      it++;
      for( int i=0 ; i < ns ; i++ )
	 dfs[i] = dfps[i];
      for( int i=0 ; i < nmpar ; i++ )
	 dfm[i] = dfpm[i];
   }
   if( myRank == 0 )
   {
      fclose(fd);
      fclose(fdx);
   }
   delete[] al;
   delete[] rho;
   if( ns > 0 )
   {
      delete[] s;
      delete[] y;
      delete[] dftemp;
   }
   if( nmpar > 0 )
   {
      delete[] dfm;
      delete[] sm;
      delete[] ym;
      delete[] dfmtemp;
      delete[] dm;
      delete[] dam;
      delete[] xam;
      delete[] dfpm;
   }
}
