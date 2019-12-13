#ifdef OPTTEST_MODE
#include "dummy-classes.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <mpi.h>
#else
#include "EW.h"
#include "MaterialParameterization.h"
#include "Mopt.h"
#include "compute_f.h"
#endif

#ifdef USE_HDF5
#include <sachdf5.h>
#endif

using namespace std;

//-----------------------------------------------------------------------
void nlcg( EW& simulation, int nspar, int nmpars, double* xs, 
	   int nmpard, double* xm, 
	   vector<vector<Source*> >& GlobalSources,
	   vector<vector<TimeSeries*> >& GlobalTimeSeries,
	   vector<vector<TimeSeries*> >& GlobalObservations,
	   int myRank, Mopt* mopt )

{
   int ns, verbose = -1, nreductions = 0;
   double rnorm, f;
   bool testing=false;
   int maxrestart = mopt->m_maxit;
   int maxit      = mopt->m_maxsubit;
   double tolerance = mopt->m_tolerance;
   bool fletcher_reeves = mopt->m_fletcher_reeves, dolinesearch=mopt->m_dolinesearch;
   double* sfs = mopt->m_sfs;
   double* sfm = mopt->m_sfm;

   ns = nspar + nmpars;

   if( maxrestart == 0 )
      return;

   int nmpard_global=0;
   MPI_Allreduce( &nmpard, &nmpard_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

   if( maxit <= 0 && nmpard == 0 )
      maxit = ns + nmpard;
   else if( maxit <= 0 )
      maxit = 10;

   FILE *fd;
   FILE *fdx;
   const string parfile = mopt->m_path + "parameters.bin";
   if( myRank == 0 )
   {
      const string convfile = mopt->m_path + "convergence.log";

      fd = fopen(convfile.c_str(),"w");
      fprintf(fd, "it  sub-it  max-nrm-gradient  max-nrm-model-update  misfit\n");
      
      const string parafile = mopt->m_path + "parameters.log";
      if( nspar > 0 )
	 fdx=fopen(parafile.c_str(),"w");
   }
   if( myRank == 0 )
      cout << "Begin NLCG iteration by evaluating initial misfit and gradient..." << endl;

   double* dfs, *dfm, *ds, *dm, *das, *dam, *xas, *xam, *dfps, *dfpm;
   if( ns > 0 )
   {
      dfs = new double[ns];
      ds  = new double[ns];
      das = new double[ns];
      xas = new double[ns];
      dfps= new double[ns];
   }
   if( nmpard > 0 )
   {
      dfm  = new double[nmpard];
      dm   = new double[nmpard];
      dam  = new double[nmpard];
      xam  = new double[nmpard];
      dfpm = new double[nmpard];
   }

   compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, f, dfs, dfm, myRank, mopt );
   if( mopt->m_output_ts )
   {
     for( int e=0 ; e < GlobalTimeSeries.size() ; e++ )
     {
#ifdef USE_HDF5
       // Tang: need to create a HDF5 file before writing
       if (GlobalTimeSeries[e].size() > 0 && GlobalTimeSeries[e][0]->getUseHDF5()) {
         if(myRank == 0) 
           createTimeSeriesHDF5File(GlobalTimeSeries[e], GlobalTimeSeries[e][0]->getNsteps(), GlobalTimeSeries[e][0]->getDt(), "_ini");
         for (int tsi = 0; tsi < GlobalTimeSeries[e].size(); tsi++) 
           GlobalTimeSeries[e][tsi]->resetHDF5file();
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
	    cout << dfs[i]*sfs[i] << " ";
	    if( i==5 )
	       cout << endl << "      " ;
	 }
      }
   }
   rnorm = 0;
   if( nmpard_global > 0 )
   {
      for( int i=0 ; i < nmpard ; i++ )
	 rnorm = rnorm > fabs(dfm[i])*sfm[i] ? rnorm : fabs(dfm[i])*sfm[i];
      double rnormtmp = rnorm;
      MPI_Allreduce(&rnormtmp, &rnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   }
   for( int i=0 ; i < ns ; i++ )
      rnorm = rnorm > fabs(dfs[i])*sfs[i] ? rnorm : fabs(dfs[i])*sfs[i];
   if( myRank == 0 )
      cout << "Max norm of scaled total gradient = " << rnorm << endl;

// Start CG iterations
   int j = 1;
   bool done = false;
   int it= 0;
   while( j <= maxrestart && !done )
   {
      for( int i=0 ; i < ns ; i++ )
         ds[i] = -dfs[i]*sfs[i]*sfs[i];
      for( int i=0 ; i < nmpard ; i++ )
         dm[i] = -dfm[i]*sfm[i]*sfm[i];
      int k = 1;
      while( k <= maxit && rnorm > tolerance )
      {  
	 double alpha, fp, dtHd, h=1e-6;
         double hi=1/h;

	 // Compute step length
	 if( myRank == 0 && verbose > 2 )
	    cout << "Step length computation " << endl;

         // Perturb x in the direction of d
	 double normd=0;
	 if( nmpard_global > 0 )
	 {
	    double normdlocal=0;
	    for( int i=0 ; i < nmpard ; i++ )
	       normdlocal += dm[i]*dm[i];
	    MPI_Allreduce( &normdlocal, &normd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 }
	 for( int i=0 ; i < ns ; i++ )
	    normd += ds[i]*ds[i];
	 normd = sqrt(normd);
	 h = 1e-6/normd;
	 hi = 1/h;

	 for( int i=0 ; i < ns ; i++ )
	    xas[i] = xs[i] + h*ds[i];
	 for( int i=0 ; i < nmpard ; i++ )
	    xam[i] = xm[i] + h*dm[i];

	 compute_f_and_df( simulation, nspar, nmpars, xas, nmpard, xam, GlobalSources, GlobalTimeSeries,
			   GlobalObservations, fp, dfps, dfpm, myRank, mopt );
         dtHd  = 0;
	 alpha = 0;
	 if( nmpard_global > 0 )
	 {
            double dtHdlocal=0, alphaloc=0;
            for( int i=0 ; i < nmpard ; i++ )
	    {
               dtHdlocal += dm[i]*(dfpm[i]-dfm[i])*hi;
               alphaloc  += dfm[i]*dm[i];
	    }
	    MPI_Allreduce( &dtHdlocal, &dtHd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	    MPI_Allreduce( &alphaloc, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 }
	 for( int i=0 ; i < ns ; i++ )
	 {
	    dtHd += ds[i]*(dfps[i]-dfs[i])*hi;
	    alpha += ds[i]*dfs[i];
	 }
	 alpha = -alpha/dtHd;
         if( dolinesearch )
         {
	    for( int i=0 ; i < ns ; i++ )
	       das[i] = ds[i];
	    for( int i=0 ; i < nmpard ; i++ )
	       dam[i] = dm[i];
            int retcode;
            if( myRank == 0 && verbose > 2 )
	       cout << "Line search.. " << endl;
	    linesearch( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations,
			nspar, nmpars, xs, nmpard_global, nmpard, xm, f, dfs, dfm, das, dam,
			fabs(alpha), 1.0, tolerance, xas, xam, fp, sfs, sfm, myRank,
			retcode, nreductions, testing, dfps, dfpm, mopt );
            if( myRank == 0 && verbose > 2 )
	       cout << " .. return code "  << retcode << " misfit changed from " << f << " to " << fp << endl;
	 }
         else
	 {
	    for( int i=0 ; i < ns ; i++ )
	       xas[i] = xs[i] + alpha*ds[i];
	    for( int i=0 ; i < nmpard ; i++ )
	       xam[i] = xm[i] + alpha*dm[i];
	 }         

   // xa is now the new iterate
   // store x^{k+1}-x^k in d to save memory
	 double dxnorm = 0;
	 if( nmpard_global > 0 )
	 {
	    for( int i=0 ; i < nmpard ; i++ )
	    {
	       double locnorm = fabs(xm[i]-xam[i]);
	       if( fabs(xam[i])> sfm[i] )
		  locnorm /= fabs(xam[i]);
	       else
		  locnorm /= sfm[i];
	       if( locnorm > dxnorm )
		  dxnorm = locnorm;
	       xm[i]  = xam[i];
	    }
	    double dxnormloc=dxnorm;
	    MPI_Allreduce(&dxnormloc,&dxnorm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	 }
	 for( int i=0 ; i < ns ; i++ )
	 {
	    double locnorm = fabs(xs[i]-xas[i]);
	    if( fabs(xas[i])> sfs[i] )
	       locnorm /= fabs(xas[i]);
	    else
	       locnorm /= sfs[i];
	    if( locnorm > dxnorm )
	       dxnorm = locnorm;
	    xs[i]  = xas[i];
	 }
	 compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
			   GlobalObservations, f, dfps, dfpm, myRank, mopt, it );

// Compute norm of new gradient
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
	    if( fabs(dfps[i])*sfs[i] > rnorm )
	       rnorm = fabs(dfps[i])*sfs[i];

// Save the time series after each sub-iteration.
//         for( int ts=0 ; ts < GlobalTimeSeries.size() ; ts++ )
//	    GlobalTimeSeries[ts]->writeFile();

	 // Save material parameters, for restart.
	 mopt->m_mp->write_parameters(parfile.c_str(),nmpars,xs);

// Print some information on standard output
	 if( myRank == 0 )
	 {
	    cout << "-----------------------------------------------------------------------" << endl;
	    cout << " it=" << j << " " << k << " max-norm scaled gradient= " << rnorm << " max-norm model update= " << dxnorm << endl;
	    cout << " Misfit= "  << f << endl;
	 }

// Save convergence and source parameters to file
         if( myRank == 0 )
	 {
	    fprintf(fd, "%i %i %15.7g %15.7g %15.7g %i\n", j, k, rnorm, dxnorm, f, nreductions );
	    fflush(fd);

	    if( nspar>0)
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
		  cout << dfps[i]*sfs[i] << " ";
		  if( i==5 )
		     cout << endl << "      " ;
	       }
	       cout << endl;
	       fprintf(fdx, "%d  %d ", j, k );
	       for( int p=0 ; p < nspar ; p++ )
		  fprintf(fdx, "%15.7g ",xs[p] );
	       fprintf(fdx,"\n");
	       fflush(fdx);
	    }

	 }

         double beta;
         if( k < maxit )
	 {
            double num=0, mix=0, den=0;
            if( nmpard_global > 0 )
	    {
               double varloc[3]={0,0,0};
	       double vars[3]={0,0,0};
	       for( int i=0 ; i < nmpard ; i++ )
	       {
		  varloc[0] +=dfpm[i]*dfpm[i]*sfm[i]*sfm[i];
		  varloc[1] += dfm[i]*dfm[i]*sfm[i]*sfm[i];
		  varloc[2] += dfm[i]*dfpm[i]*sfm[i]*sfm[i];
	       }
	       MPI_Allreduce(varloc,vars,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	       num=vars[0];
	       den=vars[1];
	       mix=vars[2];
	    }
	    for( int i=0 ; i < ns ; i++ )
	    {
               num +=dfps[i]*dfps[i]*sfs[i]*sfs[i];
	       den +=  dfs[i]*dfs[i]*sfs[i]*sfs[i];
	       mix += dfs[i]*dfps[i]*sfs[i]*sfs[i];
	    }

	    if( fletcher_reeves )
	       beta = num/den;
	    else
	    {
	       beta = (num-mix)/den;
               if( beta <= 0 )
	       {
                  k = maxit;
		  beta = 0;
                  if( myRank == 0 && verbose > 2 )
		     cout << "P-R restart because beta <= 0 " << endl;
	       }
	    }
	 }
	 else
	    beta = 0;
         for( int i=0 ; i < nmpard ; i++ )
	 {
	    dm[i]  = -sfm[i]*sfm[i]*dfpm[i] + beta*dm[i];
            dfm[i] = dfpm[i];
	 }
         for( int i=0 ; i < ns ; i++ )
	 {
	    ds[i]  = -sfs[i]*sfs[i]*dfps[i] + beta*ds[i];
            dfs[i] = dfps[i];
	 }
         k++;
	 it++;
      }
      j++;
      done = rnorm < tolerance;
   }
   mopt->m_mp->write_parameters(parfile.c_str(),nmpars,xs);

   if( myRank == 0 )
   {
      fclose(fd);
      if( nspar > 0 )
	 fclose(fdx);
   }

   if( ns > 0 )
   {
      delete[] dfs ;
      delete[] ds  ;
      delete[] das ;
      delete[] xas ;
      delete[] dfps;
   }
   if( nmpard > 0 )
   {
      delete[] dfm;
      delete[] dm;
      delete[] dam;
      delete[] xam;
      delete[] dfpm;
   }
}
