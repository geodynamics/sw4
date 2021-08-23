#include "EW.h"
#include "DataPatches.h"
#include <cstring>
#include "version.h"
//#include "MaterialParameterization.h"
#include "MaterialParCartesian.h"
#include "Mopt.h"
#include "compute_f.h"
#include "sw4-prof.h"

#ifndef SW4_NOOMP
#include <omp.h>
#endif

#include <fcntl.h>
#include <unistd.h>
#include <fstream>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifdef USE_HDF5
#include <sachdf5.h>
#endif


double t0;

float_sw4 checkTimeSeriesLimits(vector<TimeSeries*>& v)  ;

void usage(string thereason)
{
  cout << endl
       << "sbp4mopt - Summation by parts 4th order inverse seismic wave solver"  << endl << endl
       << "Usage: sbp4mopt [-v] file.in" << endl
       << "\t -v:      prints out the version info" << endl
       << "\t file.in: an input file" << endl << endl
       << "Reason for message: " << thereason << endl;
}

//-----------------------------------------------------------------------
void set_source_pars( int nspar, double srcpars[11], double* xs )
{
   // nspar =11, all parameters=(x0,y0,z0,m_{xx}-m_{zz},t0,freq)
   // nspar =10, no freq., parameters=(x0,y0,z0,m_{xx}-m_{zz},t0)
   // nspar = 9 , no freq or t0, parameters=(x0,y0,z0,m_{xx}-m_{zz})
   // nspar = 0, source not part of inversion.
   if( nspar == 11 )
      for( int i=0; i<11 ;i++ )
	 srcpars[i] = xs[i];
   else if( nspar == 10 )
      for( int i=0; i<10 ;i++ )
	 srcpars[i] = xs[i];
   else if( nspar == 9 )
      for( int i=0; i<9 ;i++ )
	 srcpars[i] = xs[i];
   else if( nspar == 6 )
      for( int i=0; i<6 ;i++ )
	 srcpars[i+3] = xs[i];      
   else if( nspar !=0 )
      cout << "Error in set_source_pars, nspar = " << nspar
	   << " undefined case"<< endl;
}

//-----------------------------------------------------------------------
void get_source_pars( int nspar, double srcpars[11], double* xs )
{
   // nspar =11, all parameters=(x0,y0,z0,m_{xx}-m_{zz},t0,freq)
   // nspar =10, no freq., parameters=(x0,y0,z0,m_{xx}-m_{zz},t0)
   // nspar = 9 , no freq or t0, parameters=(x0,y0,z0,m_{xx}-m_{zz})
   // nspar = 0, source not part of inversion.
   if( nspar == 11 )
      for( int i=0; i<11 ;i++ )
	 xs[i] = srcpars[i];
   else if( nspar == 10 )
      for( int i=0; i<10 ;i++ )
	 xs[i] = srcpars[i];
   else if( nspar == 9 )
      for( int i=0; i<9 ;i++ )
	 xs[i] = srcpars[i];
   else if( nspar == 6 )
      for( int i=0 ; i< 6 ;i++)
	 xs[i] = srcpars[i+3];
   else if( nspar !=0 )
      cout << "Error in get_source_pars, nspar = " << nspar
	   << " undefined case"<< endl;
}

//-----------------------------------------------------------------------
void  normalize_gradient_ph( vector<Sarray>& pseudo_hessian, 
                                    vector<Sarray>& gRho,
                                    vector<Sarray>& gMu, 
                                    vector<Sarray>& gLambda, 
                                    float_sw4 eps, int phcase )
{
   int g=0;
   int ib=pseudo_hessian[g].m_ib, ie=pseudo_hessian[g].m_ie;
   int jb=pseudo_hessian[g].m_jb, je=pseudo_hessian[g].m_je;
   int kb=pseudo_hessian[g].m_kb, ke=pseudo_hessian[g].m_ke;   
   // q&d for single grid VsVp case
   float_sw4 maxnorm[3]={0,0,0};

   for( int k=kb ; k <= ke ;k++)
      for( int j=jb ; j <= je ;j++)
         for( int i=ib ; i <= ie ;i++)
            for( int c=0 ; c < 3; c++ )
         {
            if( pseudo_hessian[g](c+1,i,j,k) > maxnorm[c] )
               maxnorm[c] = pseudo_hessian[g](c+1,i,j,k);
         }
   
   float_sw4 mxnormloc[3]={maxnorm[0],maxnorm[1],maxnorm[2]};
   MPI_Allreduce( mxnormloc,maxnorm,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   for( int k=kb ; k <= ke ;k++)
      for( int j=jb ; j <= je ;j++)
         for( int i=ib ; i <= ie ;i++)
         {
            gRho[g](i,j,k)    /= pseudo_hessian[g](1,i,j,k)/maxnorm[0]+eps;
            gMu[g](i,j,k)     /= pseudo_hessian[g](2,i,j,k)/maxnorm[1]+eps;
            gLambda[g](i,j,k) /= pseudo_hessian[g](3,i,j,k)/maxnorm[2]+eps;            
         }
}

//-----------------------------------------------------------------------
void normalize_pseudohessian( int nmpars, float_sw4* phs, int nmpard, 
                              float_sw4* phd, float_sw4 eps, int phcase )
{
   float_sw4 mxnorm[3]={0,0,0};
   int ncomp;
   if( phcase==1 || phcase==2 )
      ncomp=3;
   else if( phcase==3 )
      ncomp=2;
   else
      ncomp=1;
   for( int m=0; m < nmpars ; m++ )
      if( std::isnan(phs[m]) )
         std::cout << "ph is nan at m= " << m << std::endl;

   int npts=nmpars/ncomp;
   for( int m=0; m < npts ; m++ )
      for( int c=0; c < ncomp ; c++ )
      {
         if( mxnorm[c]< phs[m*ncomp+c] )
            mxnorm[c]=phs[m*ncomp+c];
      }
   npts = nmpard/ncomp;
   for( int m=0; m < npts ; m++ )
      for( int c=0; c < ncomp ; c++ )
      {
         if( mxnorm[c]< phd[m*ncomp+c] )
            mxnorm[c]=phd[m*ncomp+c];
      }
   if( nmpard > 0 )
   {
      float_sw4 mxnormloc[3]={mxnorm[0],mxnorm[1],mxnorm[2]};
      MPI_Allreduce( mxnormloc,mxnorm,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   }
   //   std::cout << "mxnorm = " << mxnorm[0] << " " << mxnorm[1] << " " << mxnorm[2] << std::endl;
   npts=nmpars/ncomp;
   for( int m=0; m < npts ; m++ )
      for( int c=0; c < ncomp ; c++ )
         phs[m*ncomp+c] = phs[m*ncomp+c]/mxnorm[c]+eps;
   npts = nmpard/ncomp;
   for( int m=0; m < npts ; m++ )
      for( int c=0; c < ncomp ; c++ )
         phd[m*ncomp+c] = phd[m*ncomp+c]/mxnorm[c]+eps;
}

//-----------------------------------------------------------------------
void compute_f( EW& simulation, int nspar, int nmpars, double* xs,
		int nmpard, double* xm,
		vector<vector<Source*> >& GlobalSources,
		vector<vector<TimeSeries*> >& GlobalTimeSeries,
		vector<vector<TimeSeries*> >& GlobalObservations,
		double& mf, Mopt *mopt )
//-----------------------------------------------------------------------
// Compute misfit.
//
// Input: simulation - Simulation object
//        nspar  - Number of source parameters
//        nmpars - Number of shared material parameters
//        xs     - Vector of shared unknown, xs[0..nspar-1] are the source
//                 parameters, xs[nspar..nmpars+nspar-1] are the material parameters.
//        nmpard - Number of distributed material parameters
//        xm     - Vector of distributed unknown, of size nmpard.
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        mopt - Pointer to the class Mopt object
//
// Output: GlobalTimeSeries - The solution of the forward problem at the stations.
//         mf               - The misfit.
//-----------------------------------------------------------------------
{
   // Source optimization
   //   vector<Source*> src(1);
   //   src[0] = GlobalSources[0][0]->copy(" ");
   //   // fetch all 11 source parameters, and merge in the unknowns
   //   double srcpars[11];
   //   src[0]->get_parameters( srcpars );
   //   set_source_pars( nspar, srcpars, xs );
   //   src[0]->set_parameters( srcpars );


   sw4_profile->time_stamp("Enter compute_f");
// Translate one-dimensional parameter vector xm to material data (rho,mu,lambda)
   int ng = simulation.mNumberOfGrids;
   vector<Sarray> rho(ng), mu(ng), lambda(ng);

//New
   int nms, nmd, nmpard_global;
   mopt->m_mp->get_nr_of_parameters( nms, nmd, nmpard_global );
   if( nms != nmpars || nmd != nmpard )
      cout << "compute_f: WARNING, inconsistent number of material parameters" << endl;

//solveTT
   double *coarse = new double[nmpars];
   mopt->m_mp->get_base_parameters(nmpard,xm,nmpars,coarse,simulation.mRho,simulation.mMu,simulation.mLambda );
   checkMinMax(nmpars/2, coarse, "coarse:");
   float_sw4 freq;

   for( int e=0 ; e < simulation.getNumberOfEvents() ; e++ )
   {
     freq= mopt->get_freq_gradsmooth()>0.? mopt->get_freq_gradsmooth() : GlobalSources[e][0]->getFrequency();
     simulation.solveTT(GlobalSources[e][0], GlobalTimeSeries[e], coarse, nmpars, mopt->m_mp, 
     mopt->get_wave_mode(), mopt->get_twin_shift(), mopt->get_twin_scale(), freq, e, simulation.getRank());
   }
   delete[] coarse;
   MPI_Barrier(MPI_COMM_WORLD);
   
   //rho, mu lambda are local volumes, xs assigned to base mu/lambda
   mopt->m_mp->get_material( nmpard, xm, nmpars, &xs[nspar], rho, mu, lambda, 
       mopt->get_vp_min(), mopt->get_vp_max(), mopt->get_vs_min(), mopt->get_vs_max(), mopt->get_wave_mode()); 
   int ok=1;
   if( mopt->m_mcheck )
      simulation.check_material( rho, mu, lambda, ok );
   VERIFY2( ok, "ERROR: compute_f Material check failed\n" );


// Old
 //   simulation.parameters_to_material( nm, xm, rho, mu, lambda );

// Run forward problem with guessed source, upred_saved,ucorr_saved are allocated
// inside solve_allpars. U and Um are final time solutions, to be used as 'initial' data
// when reconstructing U backwards.
   vector<DataPatches*> upred_saved(ng), ucorr_saved(ng);
   vector<Sarray> U(ng), Um(ng), ph(ng);
   mf = 0;
   if( !mopt->m_test_regularizer ){
   for( int e=0 ; e < simulation.getNumberOfEvents() ; e++ )
   {
//	 simulation.solve( src, GlobalTimeSeries[e], mu, lambda, rho, U, Um, upred_saved, ucorr_saved, false, e );
      std::cout << "compute_f forward solve" << " time from t0=" << MPI_Wtime()-t0 << std::endl;
      sw4_profile->time_stamp("forward solve" );
      simulation.solve( GlobalSources[e], GlobalTimeSeries[e], mu, lambda, rho, U, Um, upred_saved, ucorr_saved, 
      false, e, mopt->m_nsteps_in_memory, 0, ph );
      sw4_profile->time_stamp("done forward solve" ); 
      std::cout << "done compute_f forward solve" << " time from t0=" << MPI_Wtime()-t0 << std::endl;

     // Wei added sync here
     MPI_Barrier(MPI_COMM_WORLD);

//        Compute misfit
      if( mopt->m_misfit == Mopt::L2 )
      {
	 double dshift, ddshift, dd1shift;

	 for( int m = 0 ; m < GlobalTimeSeries[e].size() ; m++ )
	    mf += GlobalTimeSeries[e][m]->misfit( *GlobalObservations[e][m], NULL, dshift, ddshift, dd1shift );
      }
      else if( mopt->m_misfit == Mopt::CROSSCORR )
      {
	 for( int m = 0 ; m < GlobalTimeSeries[e].size() ; m++ )
	 {
	    mf += GlobalTimeSeries[e][m]->misfit2( *GlobalObservations[e][m], NULL );
	    //	    if( e==0 && m== 0 )
	    //	       exit(0);
	 }
      }
      // Give back memory
      for( unsigned int g=0 ; g < ng ; g++ )
      {
         delete upred_saved[g];
         delete ucorr_saved[g];
      }
   }
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   }
// add in a Tikhonov regularizing term:
   bool tikhonovreg=false;
   if( tikhonovreg )
   {
      double tcoff = (1.0/(nmpars+nmpard_global))*(mopt->m_reg_coeff);
      if( tcoff != 0 )
      {
// Shared parameters
	 double tikhonov=0;
	 for (int q=nspar; q<nspar+nmpars; q++)
	    tikhonov += SQR( (xs[q] - mopt->m_xs0[q])/mopt->m_sfs[q]);
	 mf += tcoff*tikhonov;

// Distributed parameters
	 double tikhonovd = 0;
	 for (int q=0; q<nmpard; q++)
	    tikhonovd += SQR( (xm[q] - mopt->m_xm0[q])/mopt->m_sfm[q]);
	 MPI_Allreduce( &tikhonovd, &tikhonov, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 mf += tcoff*tikhonov;
      }
   }
   else
   {
      double mf_reg;
      double* dmfs, *dmfd;
      mopt->m_mp->get_regularizer( nmpard, xm, nmpars, xs, mopt->m_xm0, mopt->m_xs0,
				   mopt->m_reg_coeff, rho, mu, lambda, mf_reg, mopt->m_sfm,
				   mopt->m_sfs, false, dmfd, dmfs);
      if( mopt->m_test_regularizer )
	 mf = mf_reg;
      else
	 mf += mf_reg;
   }
   sw4_profile->time_stamp("Exit compute_f");
}

//-----------------------------------------------------------------------
void compute_f_and_df( EW& simulation, int nspar, int nmpars, double* xs,
		       int nmpard, double* xm,
		       vector<vector<Source*> >& GlobalSources,
		       vector<vector<TimeSeries*> >& GlobalTimeSeries,
		       vector<vector<TimeSeries*> >& GlobalObservations, 
		       double& f, double* dfs, double* dfm, int myrank,
                       Mopt* mopt, int it )

//		       MaterialParameterization* mp, bool mcheck, bool output_ts,
//		       vector<Image*>& images )

//-----------------------------------------------------------------------
// Compute misfit and its gradient.
//
// Input: simulation - Simulation object
//        nspar  - Number of source parameters (shared among processors).
//        nmpars - Number of shared material parameters.
//        xs     - Vector of shared unknown, xs[0..nspar-1] are the source
//                 parameters, xs[nspar..nmpars+nspar-1] are the material parameters.
//        nmpard - Number of distributed material parameters
//        xm     - Vector of distributed unknown, of size nmpard.
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        mp - Pointer to object describing the parameterization of the material.
//
// Output: GlobalTimeSeries - The solution of the forward problem at the stations.
//         f                - The misfit.
//         dfs              - Gradient wrt to the source of misfit.
//         dfm              - Gradient wrt to the material of misfit.
//-----------------------------------------------------------------------
{
   int verbose = 0;
   sw4_profile->time_stamp("Enter compute_f_and_df");
   // source optimization
   //   vector<Source*> src(1);
   //   src[0] = GlobalSources[0]->copy(" ");

   //   // fetch all 11 source parameters, and merge in the unknowns
   //   double srcpars[11];
   //   src[0]->get_parameters( srcpars );
   //   set_source_pars( nspar, srcpars, xs );
   //   src[0]->set_parameters( srcpars );

// Translate one-dimensional parameter vector (xm,xs) to material data (rho,mu,lambda)
   int ng = simulation.mNumberOfGrids;
   vector<Sarray> rho(ng), mu(ng), lambda(ng);

   int nms, nmd, nmpard_global;
   int ok=1;

     //for( int m = 0 ; m < GlobalObservations[0].size() ; m++ ) 
        //if(GlobalObservations[0][m]->myPoint()) cout << "obs check station=" << GlobalObservations[0][m]->getStationName() << " start time=" << GlobalObservations[0][m]->getStartTime() << " shift=" <<  GlobalObservations[0][m]->getTimeShift() << endl;

   mopt->m_mp->get_nr_of_parameters( nms, nmd, nmpard_global );
   if( nms != nmpars || nmd != nmpard )
      cout << "compute_f_and_df: WARNING, inconsistent number of material parameters" << endl;


   // solveTT
   double *coarse = new double[nmpars]; // globally shared base model for computing TT
   mopt->m_mp->get_base_parameters(nmpard,xm,nmpars,coarse,simulation.mRho,simulation.mMu,simulation.mLambda );
   checkMinMax(nmpars/2, coarse, "coarse2:");
   float_sw4 freq;
   for( int e=0 ; e < simulation.getNumberOfEvents() ; e++ )
   {
     freq= mopt->get_freq_gradsmooth()>0.? mopt->get_freq_gradsmooth() : GlobalSources[e][0]->getFrequency();
     simulation.solveTT(GlobalSources[e][0], GlobalTimeSeries[e], coarse, nmpars, mopt->m_mp, 
     mopt->get_wave_mode(), mopt->get_twin_shift(), mopt->get_twin_scale(), freq, e, myrank);
   }
    MPI_Barrier(MPI_COMM_WORLD);
   delete[] coarse;

   checkMinMax(nmpars, &xs[nspar], "compute_f_and_df: xs");
   //rho,mu,lambda local volumes get updated
   mopt->m_mp->get_material( nmpard, xm, nmpars, &xs[nspar], rho, mu, lambda,
     mopt->get_vp_min(), mopt->get_vp_max(), mopt->get_vs_min(), mopt->get_vs_max(), mopt->get_wave_mode());

   if( mopt->m_mcheck )
      simulation.check_material( rho, mu, lambda, ok, 2 );
   VERIFY2( ok, "ERROR: Material check failed\n" );

// Run forward problem with guessed source, upred_saved,ucorr_saved are allocated
// inside solve_allpars. U and Um are final time solutions, to be used as 'initial' data
// when reconstructing U backwards.
   vector<DataPatches*> upred_saved(ng), ucorr_saved(ng);
   vector<Sarray> U(ng), Um(ng);
   vector<Sarray> gRho(ng), gMu(ng), gLambda(ng);  // gradients of model parms
   vector<Sarray> pseudo_hessian(ng);
   f = 0;
   
   float_sw4 *dfsevent, *dfmevent;
   if( nmpars > 0 )
      dfsevent = new float_sw4[nmpars];
   if( nmpard > 0 )
      dfmevent = new float_sw4[nmpard];

   for( int m=0 ; m < nmpars ; m++ )
      dfs[m+nspar] = 0;
   for( int m=0 ; m < nmpard ; m++ )
      dfm[m] = 0;

  if( !mopt->m_test_regularizer )
  {
        mopt->init_pseudohessian( pseudo_hessian );
        int phcase = mopt->get_pseudo_hessian_case();

   for( int e=0 ; e < simulation.getNumberOfEvents() ; e++ )
   {
      std::cout << "compute_f_df forward solve"  << " time from t0=" << MPI_Wtime()-t0 << std::endl;
      sw4_profile->time_stamp("forward solve" );
       simulation.solve( GlobalSources[e], GlobalTimeSeries[e], mu, lambda, rho, U, Um, upred_saved, ucorr_saved, true, e, 
     mopt->m_nsteps_in_memory, phcase, pseudo_hessian );

      // check if U has nan
      U[0].checknan("forward U");
      sw4_profile->time_stamp("done forward solve" );
      std::cout << "done compute_f_df forward solve" << " time from t0=" << MPI_Wtime()-t0 << std::endl;
      
      //check simulation limits
     cout << "forward modeling GlobalTimeSeries max=" << checkTimeSeriesLimits(GlobalTimeSeries[e]) << endl;
     //Wei added sync
     MPI_Barrier(MPI_COMM_WORLD);
     
     
// Compute misfit, 'diffs' will hold the source for the adjoint problem

// 1. Copy computed time series into diffs[m]
     
      vector<TimeSeries*> diffs; // station# * npts

      for( int m=0 ; m < GlobalTimeSeries[e].size() ; m++ )
      {
	    if( mopt->m_output_ts && it >= 0 )
	       GlobalTimeSeries[e][m]->writeFile();
       TimeSeries *elem = GlobalTimeSeries[e][m]->copy( &simulation, "diffsrc", false);  // Wei add true to append filename substring
	    diffs.push_back(elem);
      }

      //for( int m=0 ; m < GlobalTimeSeries[e].size() ; m++ )
            //if(diffs[m]->myPoint()) std::cout << "rank=" << myrank << " m=" << m << " max copy of syn=" << diffs[m]->getMaxValue(0) << std::endl;

            if(myrank == 0) 
            {
              if( mopt->m_misfit == Mopt::L2 )
                  createTimeSeriesHDF5File(diffs, GlobalTimeSeries[e][0]->getNsteps(), simulation.getTimeStep(), "_adj_l2");
              if(mopt->m_misfit == Mopt::CROSSCORR)
                  createTimeSeriesHDF5File(diffs, GlobalTimeSeries[e][0]->getNsteps(), simulation.getTimeStep(), "_adj_cross");
            }
            MPI_Barrier(MPI_COMM_WORLD); 
      
      cout << "copy of GlobalTimeSeries max=" << checkTimeSeriesLimits(diffs) << endl;

// 2. misfit function also updates diffs := this - observed
      if( mopt->m_misfit == Mopt::L2 )
      {
	 double dshift, ddshift, dd1shift;

     for( int m = 0 ; m < GlobalTimeSeries[e].size() ; m++ ) 
     {

	     f += GlobalTimeSeries[e][m]->misfit( *GlobalObservations[e][m], diffs[m], dshift, ddshift, dd1shift );
      // QC adj source by wei
         #if USE_HDF5
         // Allocate HDF5 fid for later file write
         if(m == 0) diffs[0]->allocFid();
         else       diffs[m]->setFidPtr(diffs[0]->getFidPtr());
                     
         diffs[m]->setTS0Ptr(diffs[0]);
         diffs[m]->syncSolFloats();
         diffs[m]->writeFile("_adj_l2");
         if(diffs[m]->myPoint()) std::cout << "rank=" << myrank << " m=" << m << " max obs=" << GlobalObservations[e][m]->getMaxValue(0) 
               << " syn=" << GlobalTimeSeries[e][m]->getMaxValue(0) << " adj=" << diffs[m]->getMaxValue(0) << std::endl;
         #endif
       } // loop over m

        std::cout << ">>>>>>>>>>>>>>>>>>> adj source written to files" << std::endl;
      }
      else if( mopt->m_misfit == Mopt::CROSSCORR )
         {
            for( int m = 0 ; m < GlobalTimeSeries[e].size() ; m++ )
            {
        //cout << "obs start time=" << GlobalObservations[e][m]->getStartTime() << " shift=" <<  GlobalObservations[e][m]->getTimeShift() << endl;
        //cout << "syn start time=" << GlobalTimeSeries[e][m]->getStartTime() << " shift=" <<  GlobalTimeSeries[e][m]->getTimeShift() << endl;

               f += GlobalTimeSeries[e][m]->misfit2( *GlobalObservations[e][m], diffs[m] );

               #if USE_HDF5
               // Allocate HDF5 fid for later file write
               if(m == 0) diffs[0]->allocFid();
               else       diffs[m]->setFidPtr(diffs[0]->getFidPtr());
                           
               diffs[m]->setTS0Ptr(diffs[0]);
               diffs[m]->syncSolFloats();
               diffs[m]->writeFile("_adj_cross");  // myPoint checked internally
               if(diffs[m]->myPoint()) {
                  std::cout << "rank=" << myrank << " m=" << m << " max obs=" << GlobalObservations[e][m]->getMaxValue(0) 
                     << " syn=" << GlobalTimeSeries[e][m]->getMaxValue(0) << " adj=" << diffs[m]->getMaxValue(0) << std::endl;
               }
               #endif
            } // loop over m
         } // end of Mopt

      double dfsrc[11];
      get_source_pars( nspar, dfsrc, dfs );   
      std::cout << "backward+adjoint solve" << " time from t0=" << MPI_Wtime()-t0 << std::endl;

      // debug break for windowing 
      //MPI_Barrier(MPI_COMM_WORLD);
      //MPI_Abort(MPI_COMM_WORLD, 1);

      sw4_profile->time_stamp("backward+adjoint solve" );
      
        for( int m = 0 ; m < diffs.size() ; m++ ) {
         if(diffs[m]->myPoint() && myrank==0) {
            diffs[m]->getMaxValue(0);
            //std::cout << "rank=" << myrank << " m=" << m << " max obs=" << GlobalObservations[e][m]->getMaxValue(0) 
            //<< " syn=" << GlobalTimeSeries[e][m]->getMaxValue(0) << " adj=" << diffs[m]->getMaxValue(0) << std::endl;
            }
        }

      simulation.solve_backward_allpars( GlobalSources[e], rho, mu, lambda,  diffs, U, Um, upred_saved, ucorr_saved, dfsrc, gRho, gMu, gLambda, e );
      sw4_profile->time_stamp("done backward+adjoint solve" );
      cout << "done adjoint solve:" << " time from t0=" << MPI_Wtime()-t0 << " gLambda[0] npts=" << gLambda[0].npts() << " min=" << gLambda[0].minimum() << " max=" << gLambda[0].maximum() << endl;

      //      mopt->m_mp->get_gradient( nmpard, xm, nmpars, &xs[nspar], &dfs[nspar], dfm, gRho, gMu, gLambda );
      //      mopt->m_mp->gradient_transformation( rho, mu, lambda, gRho, gMu, gLambda );
      
      //gMu[0].save_to_disk("gMu_ini.say");
      //gLambda[0].save_to_disk("gLambda_ini.say"); // local size with halos
      //std::cout << "gMu_ini min=" << gMu[0].minimum() << " max=" << gMu[0].maximum() << std::endl;
      //std::cout << "gLambda_ini min=" << gLambda[0].minimum() << " max=" << gLambda[0].maximum() << std::endl;

      //gMu[0].save_to_disk("gMu.say");
      //gLambda[0].save_to_disk("gLambda.say"); // local size with halos
      //std::cout << "gMu min=" << gMu[0].minimum() << " max=" << gMu[0].maximum() << std::endl;
      //std::cout << "gLambda min=" << gLambda[0].minimum() << " max=" << gLambda[0].maximum() << std::endl;

      mopt->m_mp->get_gradient( nmpard, xm, nmpars, &xs[nspar], dfsevent, dfmevent, 
      rho, mu, lambda, gRho, gMu, gLambda);

      //gMu[0].save_to_disk("gMu_xform.say");
      //gLambda[0].save_to_disk("gLambda_xform.say"); // local size with halos

      //std::cout << "gMu_xform min=" << gMu[0].minimum() << " max=" << gMu[0].maximum() << std::endl;
      //std::cout << "gLambda_xform min=" << gLambda[0].minimum() << " max=" << gLambda[0].maximum() << std::endl;
    
      //save_array_to_disk(nmpars, dfsevent, "dfsevent.bin");

      //mopt->m_mp->smooth_gradient(dfsevent); // needs to act on dfsevent instead of dfs
    if( phcase > 0)
      {
// Interpolate pseudo-hessian to parameter grid
         float_sw4* phs=0, *phm=0;
         if( nmpars > 0 )
            phs = new float_sw4[nmpars];
         if(  nmpard > 0 )
            phm = new float_sw4[nmpard];
         mopt->m_mp->interpolate_pseudohessian(nmpars, phs, nmpard, phm, pseudo_hessian);
         float_sw4 eps=1e-3;
         normalize_pseudohessian( nmpars, phs, nmpard, phm, eps, phcase );

         std::cout << "nmpars=" << nmpars << std::endl;
         
         //mopt->m_mp->smooth_gradient(phs); // smooth phs
         //save_array_to_disk(nmpars, phs, "phs.bin");

// ..scale the gradient
         //  float_sw4* sfs=mopt->m_sfs;
         for( int m=0 ; m < nmpars ; m++ )
            dfsevent[m+nspar] *= 1.0/phs[m];
         for( int m=0 ; m < nmpard ; m++ )
            dfmevent[m] *= 1.0/phm[m];

         //save_array_to_disk(nmpars, dfsevent, "dfsevent_scaled.bin");
// ..and give back memory
         if( nmpars> 0 )
            delete[] phs;
         if( nmpard> 0 )
            delete[] phm;

         // For plotting purpose:
      //normalize_gradient_ph( pseudo_hessian, gRho, gMu, gLambda, eps, phcase );
      }  // if pH scaling

   
    float_sw4 freq= mopt->get_freq_gradsmooth()>0.? mopt->get_freq_gradsmooth() : GlobalSources[e][0]->getFrequency();
    mopt->m_mp->smooth_gradient(dfsevent, rho, mu, lambda, freq, GlobalSources[e][0]->getZ0()); // needs to act on dfsevent instead of dfs
    //save_array_to_disk(nmpars, dfsevent, "dfsevent_smoothed.bin");

      for( int m=0 ; m < nmpars ; m++ )
	      dfs[m+nspar] += dfsevent[m];
      for( int m=0 ; m < nmpard ; m++ )
	      dfm[m] += dfmevent[m];

      
// 3. Give back memory
      for( unsigned int m = 0 ; m < GlobalTimeSeries[e].size() ; m++ )
	        delete diffs[m];
      diffs.clear();

   for( unsigned int g=0 ; g < ng ; g++ )
   {
      delete upred_saved[g];
      delete ucorr_saved[g];
   }

   }  // over all events

   double mftmp = f;
   MPI_Allreduce(&mftmp,&f,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   
   /*
      if( phcase > 0)
      {
// Interpolate pseudo-hessian to parameter grid
         float_sw4* phs=0, *phm=0;
         if( nmpars > 0 )
            phs = new float_sw4[nmpars];
         if(  nmpard > 0 )
            phm = new float_sw4[nmpard];
         mopt->m_mp->interpolate_pseudohessian(nmpars, phs, nmpard, phm, pseudo_hessian);
         float_sw4 eps=1e-3;
         normalize_pseudohessian( nmpars, phs, nmpard, phm, eps, phcase );

// ..scale the gradient
         //         float_sw4* sfs=mopt->m_sfs;
         for( int m=0 ; m < nmpars ; m++ )
            dfs[m+nspar] *= 1.0/phs[m];
         //         float_sw4* sfm=mopt->m_sfm;
         for( int m=0 ; m < nmpard ; m++ )
            dfm[m] *= 1.0/phm[m];

// precondition gradients further with gaussian smoothing
         //mopt->m_mp->smooth_gradient(dfs);

// ..and give back memory
         if( nmpars> 0 )
            delete[] phs;
         if( nmpard> 0 )
            delete[] phm;

         // For plotting purpose:
         normalize_gradient_ph( pseudo_hessian, gRho, gMu, gLambda, eps, phcase );
      }  // if pH scaling
 */

  } // test_regularizer
  
// add in a Tikhonov regularizing term:
   bool tikhonovreg=false;
   if( tikhonovreg )
   {
      double tcoff = (1.0/(nmpars+nmpard_global))*(mopt->m_reg_coeff);
      if( tcoff != 0 )
      {
	 double tikhonov=0;
	 for (int q=nspar; q<nspar+nmpars; q++)
	 {
	    tikhonov +=  SQR( (xs[q] - mopt->m_xs0[q])/mopt->m_sfs[q]);
	    dfs[q] += 2*tcoff*(xs[q] - mopt->m_xs0[q])/SQR(mopt->m_sfs[q]);
	 }
	 f += tcoff*tikhonov;

	 double tikhonovd = 0;
	 for (int q=0; q<nmpard; q++)
	 {
	    tikhonovd += SQR( (xm[q] - mopt->m_xm0[q])/mopt->m_sfm[q]);
	    dfm[q] += 2*tcoff*(xm[q] - mopt->m_xm0[q])/SQR(mopt->m_sfm[q]);
	 }
	 MPI_Allreduce( &tikhonovd, &tikhonov, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	 f += tcoff*tikhonov;
      }
   }
   else
   {
      double mf_reg;
      mopt->m_mp->get_regularizer( nmpard, xm, nmpars, xs, mopt->m_xm0, mopt->m_xs0,
				   mopt->m_reg_coeff, rho, mu, lambda, mf_reg, mopt->m_sfm,
				   mopt->m_sfs, true, dfmevent, dfsevent);
      if( mopt->m_test_regularizer )
      {
	 f = mf_reg;
	 for( int m=0 ; m < nmpars ; m++ )
	    dfs[m+nspar] = dfsevent[m];
	 for( int m=0 ; m < nmpard ; m++ )
	    dfm[m] = dfmevent[m];
      }
      else
      {
	 f += mf_reg;
	 for( int m=0 ; m < nmpars ; m++ )
	    dfs[m+nspar] += dfsevent[m];
	 for( int m=0 ; m < nmpard ; m++ )
	    dfm[m] += dfmevent[m];
      }
   }

   if( myrank == 0 && verbose >= 1 )
   {
      cout.precision(16);  
      cout << " Misfit (objective functional) is f = " << f << endl;
   }

// Get gradient by solving the adjoint problem:
//   double dfsrc[11];
   //   vector<Sarray> gRho(ng), gMu(ng), gLambda(ng);
   //   simulation.solve_backward_allpars( src, rho, mu, lambda,  diffs, U, Um, upred_saved, ucorr_saved, dfsrc, gRho, gMu, gLambda );
   //   get_source_pars( nspar, dfsrc, dfs );   
   //   mopt->m_mp->get_gradient( nmpard, xm, nmpars, &xs[nspar], &dfs[nspar], dfm, gRho, gMu, gLambda );
   
// Give back memory
   if( nmpars > 0 )
      delete[] dfsevent;
   if( nmpard > 0 )
      delete[] dfmevent;




   sw4_profile->time_stamp("Save images");   
   if( it >= 0 )
   {
// 2D images     
     for( int im=0 ; im < mopt->m_image_files.size() ; im++ )
     {
       int ng=simulation.mNumberOfGrids;
       Image* image = mopt->m_image_files[im];
       if( image->timeToWrite( it ) )
       {
	 if(image->mMode == Image::RHO )
	   image->computeImageQuantity(rho, 1);
	 else if(image->mMode == Image::MU )
	   image->computeImageQuantity(mu, 1);
	 else if(image->mMode == Image::LAMBDA )
	   image->computeImageQuantity(lambda, 1);
	 else if(image->mMode == Image::P )
	   image->computeImagePvel(mu, lambda, rho);
	 else if(image->mMode == Image::S )
	   image->computeImageSvel(mu, rho);
	 else if(image->mMode == Image::GRADRHO )
	   image->computeImageQuantity( gRho, 1 );
	 else if(image->mMode == Image::GRADMU )
	   image->computeImageQuantity( gMu, 1 );
	 else if(image->mMode == Image::GRADLAMBDA )
	   image->computeImageQuantity( gLambda, 1 );
	 else if(image->mMode == Image::GRADP )
	   image->compute_image_gradp( gLambda, mu, lambda, rho );
	 else if(image->mMode == Image::GRADS )
	   image->compute_image_grads( gMu, gLambda, mu, rho );
	 //	    string path = simulation.getOutputPath();
	 //	    image->writeImagePlane_2( it, path, 0 );
	 image->writeImagePlane_2( it, mopt->m_path, 0 );
       } // end if time to write
     } // end for all images
     
     // 3D images
     EW *ew_ptr = mopt->get_EWptr();
     
     for( int i3=0 ; i3<mopt->m_3dimage_files.size() ; i3++ )
// rho, mu occur twice because we don't save Up, Qp and Qs.
       mopt->m_3dimage_files[i3]->update_image( it, 0.0, 1.0, rho, rho, mu, lambda,
						gRho, gMu, gLambda, rho, mu, mopt->m_path,
						ew_ptr->mZ ); 

   } // end if it>=0
   sw4_profile->time_stamp("Exit compute_f_and_df");   
}


//-----------------------------------------------------------------------
void restrict( int active[6], int wind[6], double* xm, double* xmi )
{
   int ni = (wind[1]-wind[0]+1);
   int nj = (wind[3]-wind[2]+1);
   int nia= (active[1]-active[0]+1);
   int nja= (active[3]-active[2]+1);

   for( int k=wind[4] ; k<=wind[5] ; k++ )
      for( int j=wind[2] ; j<=wind[3] ; j++ )
	 for( int i=wind[0] ; i<=wind[1] ; i++ )
	 {
            size_t indi = (i-wind[0]) + ni*(j-wind[2]) + ni*nj*(k-wind[4]);
	    size_t ind  = (i-active[0]) + nia*(j-active[2]) + nia*nja*(k-active[4]);
	    xmi[3*indi]   = xm[3*ind];
	    xmi[3*indi+1] = xm[3*ind+1];
	    xmi[3*indi+2] = xm[3*ind+2];
	 }
}

//-----------------------------------------------------------------------
void gradient_test( EW& simulation, vector<vector<Source*> >& GlobalSources, 
		    vector<vector<TimeSeries*> >& GlobalTimeSeries,
		    vector<vector<TimeSeries*> >& GlobalObservations, 
		    int nspar, int nmpars, double* xs, int nmpard, double* xm,
		    int myRank, Mopt* mopt )
{
   // nspar:  Number of parameters in source description, when solving for the source
   // nmpars: Number of parameters in material description, non-distributed
   // nmpard: Number of parameters in material description, distributed over the mpi tasks
   //
   int ns = nspar+nmpars;
   double* dfs;
   if( ns > 0 )
      dfs = new double[ns];
   double* dfm;
   if( nmpard > 0 )
      dfm = new double[nmpard];

   //   int sharedpars = 1;

   double f, fp;
      
   vector<Image*> im;
   compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, f, dfs, dfm, myRank, mopt ); //mp, false, false, im );
   if( myRank == 0 )
   {
//      cout << "Initial f = " << f << " " << endl;
     printf("Unperturbed objective function f=%e\n", f);
   }
   
   double h=1e-6;

   int nms, nmd, nmpard_global;
   mopt->m_mp->get_nr_of_parameters( nms, nmd, nmpard_global ) ;
   std::ofstream dftest;
   if( (ns>0 || nmpard_global > 0) && myRank == 0 )
   {
      //      string fname = simulation.getOutputPath()+"GradientTest.txt";
      string fname = mopt->m_path+"GradientTest.txt";
      dftest.open(fname.c_str());
   }
   if( ns>0 )
   {
     double* sf = mopt->m_sfs;
      if( myRank == 0 )
      {
	printf("Gradient testing shared parameters :\n");
	printf("Param#  f-perturbed     step-size     f'-adjoint      f'-div-diff    error \n");
      }
      
      for( int ind=0 ; ind < ns ; ind++ )
      {
	 h = 3e-8*sf[ind]; // multiply step length by scale factor
	 xs[ind] += h;
	 compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		 GlobalObservations, fp, mopt );
	 double dfnum = (fp-f)/h;
	 double dfan  = dfs[ind];
         double relerr = fabs(dfan-dfnum)/(fabs(dfan)+1e-10);
	 if( myRank == 0/* && relerr > 1e-6*/ )
	 {
	   printf("%4d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n", ind, fp, h, dfan, dfnum, dfan-dfnum);
	   
            // cout << " ind = " << ind << "f = " << fp << " h= " << h << " dfan = " << dfan
	    // 	 << " dfnum = " << dfnum << " err = " << dfan-dfnum << endl;
	 }
	 if( myRank == 0 )
	 {
	   dftest << ind << " " << fp << " " << h << " " << dfan << " " << dfnum << " " << dfan-dfnum << endl;
	 }

	 xs[ind] -= h; // reset parameter #ind
	 if( ind > 0 && (ind % 100 == 0) && myRank == 0 )
	    cout << "Done ind = " << ind << endl;
      }
   }
   if( nmpard_global > 0 )
   {
      double* sf = mopt->m_sfm;
      if( myRank == 0 )
	 cout << "Gradient testing distributed parameters :" << endl;
      for( size_t indg = 0 ; indg < nmpard_global ; indg++ )
      {
         ssize_t ind = mopt->m_mp->local_index(indg);
	 if( ind >=0 )
         {
 	    h = 3e-8*sf[ind];
	    xm[ind] += h;
	 }
	 compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		 GlobalObservations, fp, mopt );
	 double dfnum = (fp-f)/h;
	 double dfan;
	 if( ind >=0 )
	    dfan = dfm[ind];

	 //	 if( myRank == 0 )
	 if( ind >= 0 )
	 {	    
	   cout << "( " << myRank << "), f = " << fp << " h= " << h << " dfan = " << dfan << " dfnum = " << dfnum << " err = " << dfan-dfnum << " relerr= " << (dfan-dfnum)/dfnum << endl;
	    dftest << ind << " " << fp << " " << h << " " << dfan << " " << dfnum << " " << dfan-dfnum << endl;
	 }
         if( ind >= 0 )
	    xm[ind] -= h;
      }
   }
   if( dftest.is_open() )
      dftest.close();

   if( ns > 0 )
      delete[] dfs;
   if( nmpard > 0 )
      delete[] dfm;
}

//-----------------------------------------------------------------------
void hessian_test( EW& simulation, vector<vector<Source*> >& GlobalSources, 
		   vector<vector<TimeSeries*> >& GlobalTimeSeries,
		   vector<vector<TimeSeries*> >& GlobalObservations, 
		   int nspar, int nmpars, double* xs, int nmpard, double* xm,
		   //		   int myRank, MaterialParameterization* mp, double* sf, double* sfm )
		   int myRank, Mopt* mopt )
{
	      // Hessian test
   int ns = nspar+nmpars;
   double f;
   double* dfs;
   if( ns > 0 )
      dfs = new double[ns];
   double* dfm;
   if( nmpard > 0 )
      dfm = new double[nmpard];

   vector<Image*> im;
   compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, f, dfs, dfm, myRank, mopt ); //mp, false, false, im );

   double* dfsp;
   if( ns > 0 )
      dfsp= new double[ns]; 

   double* dfmp;
   if( nmpard > 0 )
      dfmp= new double[nmpard]; 

   int sharedpars = 1;
   bool ascii_output = true;
   double h =1e-6;
   int grid = 0;

   double* sf= mopt->m_sfs;
   if( sharedpars == 1 )
   {
      double* hess = new double[ns*ns];
      double fp;
      if( myRank == 0 )
	 cout << "Hessian computation of shared parameters :" << endl;
      
      for( int ind=0 ; ind < ns ; ind++ )
      {
         h = 3e-8*sf[ind];
	 xs[ind] += h;
	 compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
			   GlobalObservations, fp, dfsp, dfmp, myRank, mopt );// mp, false, false, im );
	 for( int p= 0 ; p < ns ; p++ )
	    hess[p+ns*ind] = (dfsp[p]-dfs[p])/h;
	 xs[ind] -= h;
         if( myRank == 0 )
	    cout << " done "  << ind << endl;
      }
      if( myRank == 0 )
      {
	 //	 string fname = simulation.getOutputPath()+"hessian.bin";
	 string fname = mopt->m_path+"hessian.bin";
	 int fid = open( fname.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
	 if (fid == -1 )
	 {
	    VERIFY2(0, "ERROR: error opening file hessians.bin for writing header");
	    exit(-1);
	 }
	 int prec = 8;
	 size_t nr=write(fid,&prec,sizeof(int));
	 int npatch=1;
	 nr=write(fid,&npatch,sizeof(int));
	 int nc = 1;
	 nr=write(fid,&nc,sizeof(int));
	 double gz=simulation.mGridSize[0];
	 nr=write(fid,&gz,sizeof(double));
	 int dims[6]={1,ns,1,ns,1,1};
	 nr=write(fid,dims,6*sizeof(int));
	 nr=write(fid,hess,ns*ns*sizeof(double));
	 close(fid);
	 if( ascii_output )
	 {
	    //	    fname = simulation.getOutputPath()+"Hessian.txt";
	    fname = mopt->m_path+"Hessian.txt";
	    ofstream htest(fname.c_str());
	    for( int i=0 ; i < ns ; i++ )
	    {
	       for( int j=0 ; j < ns ; j++ )
		  htest << hess[j+ns*i] << " ";
	       htest << endl;
	    }
	    htest.close();
	 }
      }
      
   }
   else
   {
      int active[6], activeg[6];
      activeg[0] = simulation.m_iStartActGlobal[grid];
      activeg[1] = simulation.m_iEndActGlobal[grid];
      activeg[2] = simulation.m_jStartActGlobal[grid];
      activeg[3] = simulation.m_jEndActGlobal[grid];
      activeg[4] = simulation.m_kStartActGlobal[grid];
      activeg[5] = simulation.m_kEndActGlobal[grid];
      active[0]  = simulation.m_iStartAct[grid];
      active[1]  = simulation.m_iEndAct[grid];
      active[2]  = simulation.m_jStartAct[grid];
      active[3]  = simulation.m_jEndAct[grid];
      active[4]  = simulation.m_kStartAct[grid];
      active[5]  = simulation.m_kEndAct[grid];
      int wind[6];
	      // wind is intersection of 'active' and 'interior of proc.'
      for( int i=0 ; i < 6  ; i++ )
	 wind[i] = active[i];
      if( active[1] >= active[0] )
      {
	 if( wind[0] < simulation.m_iStartInt[grid] )
	    wind[0] = simulation.m_iStartInt[grid];
	 if( wind[1] > simulation.m_iEndInt[grid] )
	    wind[1] = simulation.m_iEndInt[grid];
	 if( wind[0] > wind[1] )
	    wind[1] = wind[0]-1;
      }
      if( active[3] >= active[2] )
      {
	 if( wind[2] < simulation.m_jStartInt[grid] )
	    wind[2] = simulation.m_jStartInt[grid];
	 if( wind[3] > simulation.m_jEndInt[grid] )
	    wind[3] = simulation.m_jEndInt[grid];
	 if( wind[2] > wind[3] )
	    wind[3] = wind[2]-1;
      }
      if( active[5] >= active[4] )
      {
	 if( wind[4] < simulation.m_kStartInt[grid] )
	    wind[4] = simulation.m_kStartInt[grid];
	 if( wind[5] > simulation.m_kEndInt[grid] )
	    wind[5] = simulation.m_kEndInt[grid];
	 if( wind[4] > wind[5] )
	    wind[5] = wind[4]-1;
      }
      int start[3]   ={wind[0]-activeg[0],wind[2]-activeg[2],wind[4]-activeg[4]};
      int locsize[3] ={wind[1]-wind[0]+1,wind[3]-wind[2]+1,wind[5]-wind[4]+1};
      int globsize[3]={activeg[1]-activeg[0]+1,activeg[3]-activeg[2]+1,activeg[5]-activeg[4]+1};
              
	      //              cout << myRank << " " << activeg[2] << " " << activeg[3] << " " <<  active[2] << " " << active[3] <<  " " << wind[2] << " " << wind[3] << endl;
     //              cout << myRank << " " << simulation .m_jStart[grid] << " " << simulation.m_jEnd[grid] << endl;
	      //              int start[3]   ={active[0]-activeg[0],active[2]-activeg[2],active[4]-activeg[4]};
	      //	      int locsize[3] ={active[1]-active[0]+1,active[3]-active[2]+1,active[5]-active[4]+1};
	      //	      int globsize[3]={activeg[1]-activeg[0]+1,activeg[3]-activeg[2]+1,activeg[5]-activeg[4]+1};
      int nptsbuf = 1000000;
      int iwrite = myRank == 0 || ( myRank % 8 == 0 );

      Parallel_IO* pio = new Parallel_IO(1,1,globsize,locsize,start,nptsbuf,0);
      int fid;
      //      string fname = simulation.getOutputPath()+"hessian.bin";
      string fname = mopt->m_path+"hessian.bin";
      if( myRank == 0 )
      {
		 // create file, write header
	 //	 int fid = open( "hessdir/hessians.bin", O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
	 //	 fid = open( "hessdir/hessian.bin", O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
	 int fid = open( fname.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0660 );

	 if (fid == -1 )
	 {
	    VERIFY2(0, "ERROR: error opening file hessian.bin for writing header");
	    exit(-1);
	 }
	 int prec = 8;
	 size_t nr=write(fid,&prec,sizeof(int));
	 int npatch=1;
	 nr=write(fid,&npatch,sizeof(int));
	 int nc = 3;
	 nr=write(fid,&nc,sizeof(int));
	 double gz=simulation.mGridSize[0];
	 nr=write(fid,&gz,sizeof(double));
	 int dims[6]={activeg[0],globsize[0],activeg[2],globsize[1],activeg[4],globsize[2]};
	 nr=write(fid,dims,6*sizeof(int));
      }
      else
	 fid = open( fname.c_str(), O_WRONLY );

      size_t offset = sizeof(double)+9*sizeof(int);
      size_t colsize = globsize[0]*((size_t)globsize[1])*globsize[2]*3;

	      // Verify IO
      int npts = (wind[1]-wind[0]+1)*(wind[3]-wind[2]+1)*(wind[5]-wind[4]+1);
      double* xmi;
      if( npts > 0 )
	 xmi = new double[3*npts];
	      //	          simulation.get_material_parameter( nmpar, xm );
	      //		  if( npts > 0 )
	      //		     restrict( active, wind, xm, xmi );
	      //			  pio->write_array( &fid, 3, xmi, offset, "double");
	      //			  close(fid);
	      //			  exit(0);
      for( int kper = activeg[4] ; kper <= activeg[5] ; kper++ )
	 for( int jper = activeg[2] ; jper <= activeg[3] ; jper++ )
	    for( int iper = activeg[0] ; iper <= activeg[1] ; iper++ )
	       for( int var = 0 ; var < 3 ; var++ )
	       {
	       // perturb material 
	       //	       simulation.perturb_mtrl(iper,jper,kper,h,grid,var);
	       //	       simulation.get_material_parameter( nmpard, xm );
		  ssize_t pind = mopt->m_mp->parameter_index(iper,jper,kper,grid,var);
		  if( pind >= 0 )
		     xm[pind] += h;
		  //	       mp->perturb_material(iper,jper,kper,grid,var,h,xs,xm);
		  compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
				    GlobalObservations, f, dfs, dfmp, myRank, mopt ); //mp, false, false, im );
		  for( int p= 0 ; p < nmpard ; p++ )
		     dfmp[p] = (dfmp[p]-dfm[p])/h;
			  // Save Hessian column
		  restrict( active, wind, dfmp, xmi );
		  pio->write_array( &fid, 3, xmi, offset, "double");
		  offset += colsize*sizeof(double);
		       // restore material 
	       //	       simulation.perturb_mtrl(iper,jper,kper,-h,grid,var);
	       //	       mp->perturb_material(iper,jper,kper,grid,var,-h,xs,xm);
		  if( pind >= 0 )
		     xm[pind] -= h;
		   
	       }
      close(fid);
      if( npts > 0 )
	 delete[] xmi;
   }
   if( nmpard > 0 )
   {
      delete[] dfm;
      delete[] dfmp;
   }
   if( ns > 0 )
   {
      delete[] dfs;
      delete[] dfsp;
   }
}


//-----------------------------------------------------------------------
void misfit_curve( int i, int j, int k, int var, double pmin, double pmax,
		   int npts, EW& simulation, MaterialParameterization* mp,
		   int nspar, int nmpars, double* xs, int nmpard, double* xm,
		   vector<vector<Source*> >& GlobalSources, 
		   vector<vector<TimeSeries*> >& GlobalTimeSeries,
		   vector<vector<TimeSeries*> >& GlobalObservations, int myRank,
		   Mopt* mopt )
{
   double* fcn = new double[npts];
   //   ssize_t ind=mp->parameter_index(i,j,k,0,var);
   int ind=mp->parameter_index(i,j,k,0,var);
   double xoriginal;
   bool ascii_output = true;
   double p;
   string lfname = mopt->m_path+"mflocal.txt";
   ofstream localmfout;

   int tmp=ind;
   MPI_Allreduce(&tmp,&ind,1,MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if( ind == -1 )
   {
      if( myRank == 0 )
	 cout << "misfit_curve: ERROR: index out of range " << endl;
      return;
   }

   xoriginal = xs[ind];
// Output materials
   if( mopt->m_misfit1d_images )
   {
      for( int m=0 ; m < npts ; m++ )
      {
	 if( npts==1 )
	    p = pmin;
	 else
	    p = pmin + static_cast<double>(m)/(npts-1)*(pmax-pmin);
	 xs[ind] = xoriginal+p;

	 int ng=simulation.mNumberOfGrids;
	 vector<Sarray> rho(ng), mu(ng), lambda(ng);
	 mopt->m_mp->get_material( nmpard, xm, nmpars, &xs[nspar], rho, mu, lambda,
      mopt->get_vp_min(), mopt->get_vp_max(), mopt->get_vs_min(), mopt->get_vs_max(), mopt->get_wave_mode());

	 for( int im=0 ; im < mopt->m_image_files.size() ; im++ )
	 {
	    Image* image = mopt->m_image_files[im];
	    if(image->mMode == Image::RHO )
	       image->computeImageQuantity(rho, 1);
	    else if(image->mMode == Image::MU )
	       image->computeImageQuantity(mu, 1);
	    else if(image->mMode == Image::LAMBDA )
	       image->computeImageQuantity(lambda, 1);
	    else if(image->mMode == Image::P )
	       image->computeImagePvel(mu, lambda, rho);
	    else if(image->mMode == Image::S )
	       image->computeImageSvel(mu, rho);
	    image->writeImagePlane_2( m, mopt->m_path, 0 );
	 }
      }
   }
   else
   {

      if( myRank == 0 )
	 localmfout.open(lfname.c_str());
      for( int m=0 ; m < npts ; m++ )
      {
	 if( npts==1 )
	    p = pmin;
	 else
	    p = pmin + static_cast<double>(m)/(npts-1)*(pmax-pmin);

	 xs[ind] = xoriginal+p;
	 double f;
	 compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		    GlobalObservations, f, mopt );
	 fcn[m] = f;
	 if( mopt->m_output_ts )
	 {
	    if( myRank == 0 )
	       cout << "Total misfit at p= " << p << " is " << f << " decomposition into stations and events: " << endl;

	    for( int e=0 ; e < GlobalTimeSeries.size(); e++ )
	       for( int s=0 ; s < GlobalTimeSeries[e].size() ; s++ )
	       {
		  GlobalTimeSeries[e][s]->writeFile();
		  double dshift, ddshift, dd1shift;
		  double mf; 
		  if( mopt->m_misfit == Mopt::L2 )
		     mf = GlobalTimeSeries[e][s]->misfit( *GlobalObservations[e][s], NULL, dshift, ddshift, dd1shift );
		  else if(  mopt->m_misfit == Mopt::CROSSCORR )
		     mf = GlobalTimeSeries[e][s]->misfit2( *GlobalObservations[e][s], NULL );
		  double mftmp = mf;
		  MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		  if( myRank == 0 )
		  {
		     cout << "     Event " << e+1 << " station " << s+1 << " misfit is " << mf << endl;
		     localmfout << " " << mf;
		  }
	       }
	    localmfout << endl;
	 }
      }
      if( myRank == 0 )
      {
	 localmfout.close();
      //      string fname = simulation.getOutputPath()+"fsurf.bin";
	 string fname = mopt->m_path+"fsurf.bin";
	 int fd=open(fname.c_str(), O_CREAT|O_TRUNC|O_WRONLY, 0660 );
	 int dims=1;
	 size_t nr = write(fd,&dims,sizeof(int));
	 nr = write(fd,&npts,sizeof(int));
	 nr = write(fd,&pmin,sizeof(double));
	 nr = write(fd,&pmax,sizeof(double));
	 nr = write(fd,fcn,npts*sizeof(double));
	 close(fd);
	 if( ascii_output )
	 {
	 //	 fname = simulation.getOutputPath()+"Misfit1d.txt";
	    fname = mopt->m_path+"Misfit1d.txt";
	    ofstream fcurve(fname.c_str());
	    for( int m=0 ; m < npts ; m++ )
	    {
	       if( npts > 1 )
		  fcurve << pmin + static_cast<double>(m)/(npts-1)*(pmax-pmin) << " " << fcn[m] << " " << endl;
	       else
		  fcurve << pmin  << " " << fcn[m] << " " << endl;		  
	    }
	    fcurve.close();
	 }
      }
   }
}

//-----------------------------------------------------------------------
void misfit_surface( int ix1, int jx1, int kx1, int ix2, int jx2, int kx2,
		     int varx1, int varx2, double pmin1, double pmax1,
		     double pmin2, double pmax2, int npts1, int npts2, EW& simulation,
		     MaterialParameterization* mp, int nspar, int nmpars,
		     double* xs, int nmpard, double* xm,
		     vector<vector<Source*> >& GlobalSources, 
		     vector<vector<TimeSeries*> >& GlobalTimeSeries,
		     vector<vector<TimeSeries*> >& GlobalObservations, int myRank,
		     Mopt* mopt )
{
   double* fcn = new double[npts1*npts2];
   ssize_t ind1=mp->parameter_index(ix1,jx1,kx1,0,varx1);
   ssize_t ind2=mp->parameter_index(ix2,jx2,kx2,0,varx2);
   double xoriginal1 = xs[ind1];
   double xoriginal2 = xs[ind2];
   for( int m1=0 ; m1 < npts1 ; m1++ )
   {
      for( int m2=0 ; m2 < npts2 ; m2++ )
      {
	 double p1 = pmin1 + static_cast<double>(m1)/(npts1-1)*(pmax1-pmin1);
	 xs[ind1] = xoriginal1+p1;
	 double p2 = pmin2 + static_cast<double>(m2)/(npts2-1)*(pmax2-pmin2);
	 xs[ind2] = xoriginal2+p2;
	 double f;
	 compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		    GlobalObservations, f, mopt );
	 fcn[m1+npts1*m2] = f;
      }
      if( myRank == 0 )
	 cout << "Done m = " << m1 << endl;
   }
   if( myRank == 0 )
   {
      string fname = mopt->m_path+"fsurf.bin";
      //      string fname = simulation.getOutputPath()+"fsurf.bin";
      int fd=open(fname.c_str(), O_CREAT|O_TRUNC|O_WRONLY, 0660 );
      int dims=2;
      size_t nr = write(fd,&dims,sizeof(int));
      nr = write(fd,&npts1,sizeof(int));
      nr = write(fd,&npts2,sizeof(int));
      nr = write(fd,&pmin1,sizeof(double));
      nr = write(fd,&pmax1,sizeof(double));
      nr = write(fd,&pmin2,sizeof(double));
      nr = write(fd,&pmax2,sizeof(double));
      nr = write(fd,fcn,npts1*npts2*sizeof(double));
      close(fd);
   }
}

//-----------------------------------------------------------------------
int start_minv( int argc, char **argv, string& input_file, int& myRank,
		int& nProcs )
{
  stringstream reason;

  // Initialize MPI...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // mpi2 adds on four more args...  [-p4pg, dir, -p4wd, dir]
  int mpi2args = 4;
  if (argc != 2 && argc != 3 && argc != (2+mpi2args) && argc != (3+mpi2args) )
    {
      reason << "Wrong number of args (1-2), not: " << argc-1 << endl; 
      
      if (myRank == 0)
      {
	for (int i = 0; i < argc; ++i)
	  cout << "Argv[" << i << "] = " << argv[i] << endl;

	usage(reason.str());
      }
// Stop MPI
      MPI_Finalize();
      return 1;
    }
  else if (strcmp(argv[1],"-v") == 0 )
  {
     if (myRank == 0)
        cout << ewversion::getVersionInfo() << endl;
// Stop MPI
     MPI_Finalize();
     return 2;
  }
  else
     if (argc == 1) 
     {
        reason  << "ERROR: ****No input file specified!" << endl;
        for (int i = 0; i < argc; ++i)
           reason << "Argv[" << i << "] = " << argv[i] << endl;
        
        if (myRank == 0) usage(reason.str());
// Stop MPI
	MPI_Finalize();
        return 1;
     }

  else
    input_file = argv[1];

  if (myRank == 0) 
  {
    cout << ewversion::getVersionInfo() << endl;
    cout << "Input file: " << input_file << endl;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  return 0;
}

//-----------------------------------------------------------------------
int main(int argc, char **argv)
{
  string fileName;
  int myRank, nProcs;
  int status = start_minv( argc, argv, fileName, myRank, nProcs );
  
  t0 = MPI_Wtime();

  if( status == 0 )
  {
// Save the source description here
     vector<vector<Source*> > GlobalSources; 
// Save the time series here
     vector<vector<TimeSeries*> > GlobalTimeSeries;
     vector<vector<TimeSeries*> > GlobalObservations;

// make a new simulation object by reading the input file 'fileName'
     EW simulation(fileName, GlobalSources, GlobalObservations, true );

     if (!simulation.wasParsingSuccessful())
     {
	if (myRank == 0)
	   cout << "Error: there were problems parsing the input file" << endl;
	status = 1;
     }
     else
     {
// get the simulation object ready for time-stepping
	simulation.setupRun( GlobalSources );
	if (!simulation.isInitialized())
	{ 
	   if (myRank == 0)
	      cout << "Error: simulation object not ready for time stepping" << endl;
	   status=1;
	}
	//	else if( GlobalSources[0].size() != 1 )
	//	{
	//	   if (myRank == 0)
	//	      cout << "Source optmization only implemented for a single source" << endl;
	//	}
	else
	{
// Successful initialization


	   simulation.setQuiet(true);
	   //	   simulation.setQuiet(false);
// Make observations aware of the utc reference time, if set.
// Filter observed data if required
	   GlobalTimeSeries.resize(GlobalObservations.size());

      std::cout << "Global Observations size max of e=" << GlobalObservations.size() << std::endl;

	   for( int e=0 ; e < GlobalObservations.size() ; e++ )
	   {
         for( int m = 0; m < GlobalObservations[e].size(); m++ )
         {
          if(myRank==0 && GlobalObservations[e][m]->myPoint()) 
          {
            std::cout << "obs nsamples=" << GlobalObservations[e][m]->getLastTimeStep() << " sample dt=" << GlobalObservations[e][m]->getDt() << std::endl;
            createTimeSeriesHDF5File(GlobalObservations[e], GlobalObservations[e][m]->getLastTimeStep()+1, GlobalObservations[e][m]->getDt(), "_filtered");
            break;
          }
         }
         MPI_Barrier(MPI_COMM_WORLD); // wait for file creation

	      for( int m = 0; m < GlobalObservations[e].size(); m++ )
	      {
	      //simulation.set_utcref( *GlobalObservations[m] );

            if(myRank==0 && GlobalObservations[e][m]->myPoint()) 
            std::cout << "e=" << e << " m=" << m << " GlobalObservations[e][m] Nsteps=" <<  
            GlobalObservations[e][m]->getNsteps() << " dt=" << 
            GlobalObservations[e][m]->getDt() << " zTopo=" << GlobalObservations[e][m]->getZtopo() << std::endl;
         
            //GlobalObservations[e][m]->writeFileUSGS("_obs");
            if(simulation.m_prefilter_sources && simulation.m_filter_observations )
            {
            GlobalObservations[e][m]->filter_data(simulation.m_filterobs_ptr);
            //GlobalObservations[e][m]->writeFile( "_filtered" );

               #if USE_HDF5
               // Allocate HDF5 fid for later file write
               if(m == 0) GlobalObservations[e][0]->allocFid();
               else       GlobalObservations[e][m]->setFidPtr(GlobalObservations[e][0]->getFidPtr());
                        
               GlobalObservations[e][m]->setTS0Ptr(GlobalObservations[e][0]);
               GlobalObservations[e][m]->syncSolFloats();
               GlobalObservations[e][m]->writeFile("_filtered");  // myPoint checked internally
               #endif
		      }
         } // m


MPI_Barrier(MPI_COMM_WORLD);  // added for debugging

//  First copy observations to GlobalTimeSeries, later, solve will insert 
//  the simulation time step and start time into GlobalTimeSeries.

      std::cout << "Copy observations to GlobalTimeSeries" << std::endl;

	      for( int m = 0; m < GlobalObservations[e].size(); m++ )
	      {
		 TimeSeries *elem = GlobalObservations[e][m]->copy( &simulation, "_out", true );
		 GlobalTimeSeries[e].push_back(elem);
#if USE_HDF5
                 // Allocate HDF5 fid for later file write
                 if (elem->getUseHDF5()) {
                   if(m == 0) { 
                     elem->allocFid();
                     elem->setTS0Ptr(elem);
                   }
                   else {
                     elem->setFidPtr(GlobalTimeSeries[e][0]->getFidPtr());
                     elem->setTS0Ptr(GlobalTimeSeries[e][0]);
                   }
                 }
#endif
	      } // m
	   
   }  // e


MPI_Barrier(MPI_COMM_WORLD); // added for debugging


// Configure optimizer 
   Mopt* mopt = new Mopt( &simulation );
	   mopt->parseInputFileOpt( fileName );
	   if (myRank == 0)
	   {
	      int nth=1;
#ifndef SW4_NOOMP
#pragma omp parallel
	      {
		 if( omp_get_thread_num() == 0 )
		    nth=omp_get_num_threads();
	      }
#endif
	      if( nth == 1 )
	      {
		 if( nProcs > 1 )
		    cout << "Running sw4mopt on " << nProcs << " processors..." << endl;
		 else
		    cout << "Running sw4mopt on " << nProcs << " processor..."  << endl;
	      }
	      else
	      {
		 if( nProcs > 1 )
	       // Assume same number of threads for each MPI-task.
		    cout << "Running sw4mopt on " <<  nProcs << " processors, using " << nth << " threads/processor..." << endl;
		 else
		    cout << "Running sw4mopt on " <<  nProcs << " processor, using " << nth << " threads..." << endl;
	      }
	      cout << "Writing output to directory: " << mopt->getPath() << endl;
	   }
// Create a profiling object
	   if( mopt->m_do_profiling )
	      sw4_profile = new SW4Prof(mopt->getPath());
	   else
	      sw4_profile = new SW4Prof0();

// Select material parameterization
           MaterialParameterization* mp = mopt->m_mp;

// figure out how many parameters we need.
//	Guess: nmpars - number of non-distributed (=shared) material parameters, exist copies in each proc.
//             nmpard - Number of distributed material parameters, size of part in my processor.
//	       nmpard_global - number of distributed parameters, total number over all processors
           int nmpars, nmpard, nmpard_global;
	   mp->get_nr_of_parameters( nmpars, nmpard, nmpard_global );

	   double* xm=NULL;
           if( nmpard > 0 )
	      xm = new double[nmpard];

// nspar - Number of parameters in source description. These are always non-distributed (=shared)
	   int nspar=mopt->m_nspar;
// ns - Total number of non-distributed (=shared) parameters.
           int ns = nmpars + nspar;
	   double *xs = new double[ns];

// Default initial guess, the input source, stored in GlobalSources[0], will do nothing if nspar=0.
           double xspar[11];
	   GlobalSources[0][0]->get_parameters(xspar);
	   get_source_pars( nspar, xspar, xs );

      std::cout << "nmpars=" << nmpars << std::endl;


// Initialize the material parameters
      mp->get_parameters(nmpard,xm,nmpars,&xs[nspar],simulation.mRho,simulation.mMu,simulation.mLambda );

   //           string parname = simulation.getOutputPath() + "mtrlpar-init.bin";
           string parname = mopt->m_path + "mtrlpar-init.bin";
	   mp->write_parameters(parname.c_str(),nmpars,&xs[nspar]);

// store initial material parameters in mp->m_xs0 (needed for Tikhonov regularization)
           mopt->set_baseMat(xs,xm);

// Scale factors
//	   double* sf  = NULL;
//           double* sfm = NULL;
           // if( ns > 0 )
	   //    sf  = new double[ns];
	   //           if( nmpard > 0 )
	   //              sfm = new double[nmpard];
// both the source scale factors and the shared material scale factors are in the array 'sf'
// both types of scale factors are now set in set_sscalefactors()

// Shared scale factors (source and material)
           mopt->set_sscalefactors();
// Distributed material scale factors 
	   mopt->set_dscalefactors();
// for backwards compatibility
//	   sf = mopt->m_sfs;
// tmp
	   // if (myRank == 0)
	   // {
	   //   printf("TEST: moptmain: nspar=%d, nstot=%d\n", mopt->m_nspar, mopt->m_nstot);
	   //   for (int q=0; q<mopt->m_nstot; q++)
	   //   {
	   //     printf("m_sfs[%d]=%e\n", q, mopt->m_sfs[q]);
	   //   }	     
	   // }
// end tmp
	   
// Typical sizes, initialize as scale factors
//           double* typxs=NULL;
//	   double* typxd=NULL;
//	   if( ns > 0 )
//	   {
//	      typxs = new double[ns];
//              for( int i=0; i < ns ; i++ )
//		 typxs[i] = sf[i];
//	   }
//	   if( nmpard > 0 )
//	   {
//	      typxd = new double[nmpard];
//              for( int i=0; i < nmpard ; i++ )
//		 typxd[i] = sfm[i];
//	   }

	   // Possibility to override default: typx = scale factors
	   // Commented out, not working properly atm, will fix later
	   //	   mopt->set_typx( nmpars, &sf[nspar], &typxs[nspar] );
	   //	   mopt->set_typx( nmpard, sfm, typxd );

// Output source initial guess
	   if( myRank == 0 && nspar > 0 )
	   {
	      cout << "Initial source guess : \n";
              char* names[11] = {" x0 = "," y0 = ", " z0 = "," mxx = ", " mxy = ", " mxz = ",
			   " myy = ", " myz = ", " mzz = ", " t0 = ", " freq = " };
              for( int i=0 ; i < nspar ; i++ )
	      {
		 cout << names[i] << xs[i] ;
                 if ( (i+1) % 3 == 0 || i == nspar-1 )
		    cout << endl;
	      }
	   }

           if( mopt->m_opttest == 2 )
	   {
// Gradient_test compares the computed gradient with the gradient obtained by numerical differentiation.
              gradient_test( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nspar, nmpars, xs,
			     nmpard, xm, myRank, mopt );
	   }
	   else if( mopt->m_opttest == 3 )
	   {
// Hessian_test outputs the Hessian computed by numerical differentiation.
              hessian_test( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nspar, nmpars, xs,
			    nmpard, xm, myRank, mopt );
	   }
	   else if( mopt->m_opttest == 4 )
	   {
// Compute and save a one dimensional cut through the objective function
              int npts = mopt->m_nsurfpts;
	      int ix=mopt->m_itest, jx=mopt->m_jtest, kx=mopt->m_ktest, varx=mopt->m_var;
              double pmin = mopt->m_pmin, pmax = mopt->m_pmax;
	      //              double pmin=-300,pmax=300; // rho
	      //	      double pmin=-1.57e9, pmax=1.57e9; //mu
	      //                         double pmin=-2.533e9, pmax=2.533e9; //lambda
	      //	      double pmin=-2,pmax=2;
              misfit_curve( ix, jx, kx, varx, pmin, pmax, npts, simulation, mp, nspar, nmpars, xs,
			    nmpard, xm, GlobalSources, GlobalTimeSeries, GlobalObservations, myRank, mopt );
	   }
	   else if( mopt->m_opttest == 5 )
	   {
// Compute and save a two dimensional cut through the objective function
              int npts1 = mopt->m_nsurfpts;
	      int ix1=mopt->m_itest, jx1=mopt->m_jtest, kx1=mopt->m_ktest, varx1=mopt->m_var;
              double pmin1 = mopt->m_pmin, pmax1 = mopt->m_pmax;

              int npts2 = mopt->m_nsurfpts2;
	      int ix2=mopt->m_itest2, jx2=mopt->m_jtest2, kx2=mopt->m_ktest2, varx2=mopt->m_var2;
              double pmin2 = mopt->m_pmin2, pmax2 = mopt->m_pmax2;

              misfit_surface( ix1, jx1, kx1, ix2, jx2, kx2, varx1, varx2, pmin1, pmax1,
			      pmin2, pmax2, npts1, npts2, simulation, mp, nspar, nmpars,
			      xs, nmpard, xm, GlobalSources, GlobalTimeSeries, GlobalObservations,
			      myRank, mopt );
	   }
           else if( mopt->m_opttest == 6 )
	   {
// Solve forward problem to generate synthetic data
              double f;
	      compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
			 GlobalObservations, f, mopt );
	      for( int e=0 ; e < simulation.getNumberOfEvents() ; e++ ) 
              {
#ifdef USE_HDF5
                 // Tang: need to create a HDF5 file before writing
                 if (GlobalTimeSeries[e].size() > 0 && GlobalTimeSeries[e][0]->getUseHDF5()) {
                   for (int tsi = 0; tsi < GlobalTimeSeries[e].size(); tsi++) 
                     GlobalTimeSeries[e][tsi]->resetHDF5file();
                   if(myRank == 0) 
                     createTimeSeriesHDF5File(GlobalTimeSeries[e], GlobalTimeSeries[e][0]->getNsteps(), simulation.getTimeStep(), "");
                   MPI_Barrier(MPI_COMM_WORLD);
                 }
#endif
		 for( int m=0 ; m < GlobalTimeSeries[e].size() ; m++ )
		    GlobalTimeSeries[e][m]->writeFile( );
              }
	   }
	   else if( mopt->m_opttest == 7 )
	   {
      // Project material onto a Cartesian material parameterization grid
	      CHECK_INPUT( mopt->m_mpcart0 != NULL, "ERROR, there is no Cartesian material parameterization defined\n");
	      mopt->m_mpcart0->project_and_write( simulation.mRho, simulation.mMu, simulation.mLambda,
						  "projmtrl.mpc");
	   }
           else if( mopt->m_opttest == 1 )
	   {
// Run optimizer (default)
	      sw4_profile->time_stamp("Start optimizer");

         cout << " start lbfgs time from t0=" << MPI_Wtime()-t0 << endl;
	      if( mopt->m_optmethod == 1 )
		 lbfgs( simulation, nspar, nmpars, xs, nmpard, xm, 
			GlobalSources, GlobalTimeSeries,
			GlobalObservations, myRank, mopt );


         //if(myRank==0) save_array_to_disk(nmpars, xs, "xs_lbfgs.bin"); 
         //if(myRank==0) save_array_to_disk(nmpard, xm, "xm_lbfgs.bin"); 

	      else if( mopt->m_optmethod == 2 )
		  nlcg( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		       GlobalObservations, myRank, mopt );
	      sw4_profile->time_stamp("Done optimizer");
	      sw4_profile->flush();

         cout << " done optimizer time from t0=" << MPI_Wtime()-t0 << endl;

         //Wei added to assure memory release prior to next event
         MPI_Barrier(MPI_COMM_WORLD);
	   }
           else
	      if( myRank == 0 )
		 cout << "ERROR: m_opttest = " << mopt->m_opttest << " is not a valid choice" << endl;

           {
              int ng = simulation.mNumberOfGrids;
              vector<Sarray> rho(ng), mu(ng), lambda(ng);
              mopt->m_mp->get_material( nmpard, xm, nmpars, &xs[nspar], rho, mu, lambda,
                mopt->get_vp_min(), mopt->get_vp_max(), mopt->get_vs_min(), mopt->get_vs_max(), mopt->get_wave_mode());

              for( int i3=0 ; i3<mopt->m_sfiles.size() ; i3++ )
                mopt->m_sfiles[i3]->force_write_image( 0, 0, rho, rho, mu, lambda, rho, mu, lambda, rho, lambda, simulation.getOutputPath(), simulation.mZ ); 
           }

	   if( myRank == 0 )
	   {
	      cout << "============================================================" << endl
		   << " sw4mopt ( Material/Source estimation solver) finished! " << endl
		   << "============================================================" << endl;
         cout << " time from t0=" << MPI_Wtime()-t0 << endl;
	   }
	}
     }

  }
  else if( status == 1 )
     cout  << "============================================================" << endl
	   << "The execution on proc " << myRank << " was unsuccessful." << endl
	   << "============================================================" << endl;
  if( status == 2 )
     status = 0;
// Stop MPI
  MPI_Finalize();
  return 0;
  // Note: Always return 0, to avoid having one error message per process from LC slurmd.
  //  return status;
} 


float_sw4 checkTimeSeriesLimits(vector<TimeSeries* >& v)  
{
float_sw4 value;
float_sw4 fmax=-1e20;   
   for(vector<TimeSeries*>::iterator it = v.begin(); it != v.end(); ++it)
   {
      if((*it)->myPoint()) {
         value =(*it)->getMaxValue(0);
         if(value>fmax) fmax=value;
         }
   }
return fmax;
}
