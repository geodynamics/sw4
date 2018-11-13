#include "EW.h"
#include "F77_FUNC.h"

extern "C" {
   void F77_FUNC(hdirichlet5,HDIRICHLET5)( int*, int*, int*, int*, int*, int*, int*, int*,
					   int*, int*, int*, int*, double* );
}

void EW::solve_backward_allpars( vector<Source*> & a_Sources,
				 vector<Sarray>& a_Rho, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
				 vector<TimeSeries*> & a_TimeSeries, vector<Sarray>& Up, vector<Sarray>& U,
				 vector<DataPatches*>& Upred_saved, vector<DataPatches*>& Ucorr_saved,
				 double gradientsrc[11], vector<Sarray>& gRho, vector<Sarray>& gMu, vector<Sarray>& gLambda )
{
// solution arrays
   vector<Sarray> F, Lk, Kacc, Kp, Km, K, Um, Uacc;
   //   vector<Sarray> gRho, gMu, gLambda;
   vector<Sarray*> AlphaVE, AlphaVEm, AlphaVEp;
   vector<double **> BCForcing;

   F.resize(mNumberOfGrids);
   Lk.resize(mNumberOfGrids);
   Kacc.resize(mNumberOfGrids);
   Uacc.resize(mNumberOfGrids);
   Kp.resize(mNumberOfGrids);
   Km.resize(mNumberOfGrids);
   K.resize(mNumberOfGrids);
   Um.resize(mNumberOfGrids);
   //   gRho.resize(mNumberOfGrids);
   //   gMu.resize(mNumberOfGrids);
   //   gLambda.resize(mNumberOfGrids);

// Allocate pointers, even if attenuation not used, for avoid segfault in parameter list with mMuVE[g], etc...
   AlphaVE.resize(mNumberOfGrids);
   AlphaVEm.resize(mNumberOfGrids);
   AlphaVEp.resize(mNumberOfGrids);

   BCForcing.resize(mNumberOfGrids);

   int ifirst, ilast, jfirst, jlast, kfirst, klast;
   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      BCForcing[g] = new double *[6];
      for(int side=0; side < 6; side++)
      {
	 BCForcing[g][side]=NULL;
         if (m_bcType[g][side] == bStressFree || m_bcType[g][side] == bDirichlet ||
	     m_bcType[g][side] == bSuperGrid)
         {
    	    BCForcing[g][side] = new double[3*m_NumberOfBCPoints[g][side]];
         }
      }
      ifirst = m_iStart[g];
      ilast  = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast  = m_jEnd[g];
      kfirst = m_kStart[g];
      klast  = m_kEnd[g];

      F[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Lk[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Kacc[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Kp[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Km[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      K[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Um[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      Uacc[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);

      gRho[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
      gMu[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
      gLambda[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
      gRho[g].set_to_zero();
      gMu[g].set_to_zero();
      gLambda[g].set_to_zero();
   }


// will accumulate the gradient of the misfit in these arrays
   for( int s=0 ; s < 11 ; s++ )
      gradientsrc[s] = 0;
   //   for( int s=0 ; s < nmpar ; s++ )
   //      gradientm[s] = 0;

// will accumulate the Hessian of the misfit in this array
//   for( int s=0 ; s < 121 ; s++ )
//      hessian[s] = 0;
   
// the Source objects get discretized into GridPointSource objects
   vector<GridPointSource*> point_sources;

// Transfer source terms to each individual grid as point sources at grid points.
   for( unsigned int i=0 ; i < a_Sources.size() ; i++ )
	 a_Sources[i]->set_grid_point_sources4( this, point_sources );

   //   if (!m_testing && m_prefilter_sources)
   //   {
   ////  Replace the time function by a filtered one, represented by a (long) vector holding values at each time step   
   //      for( int s=0; s < point_sources.size(); s++ ) 
   //	 point_sources[s]->discretizeTimeFuncAndFilter(mTstart, mDt, mNumberOfTimeSteps, m_filter_ptr);
   //   } 

// Sort sources wrt spatial location, needed for thread parallel computing
  vector<int> identsources;
  sort_grid_point_sources( point_sources, identsources );

   // Initial data
   for(int g=0 ; g < mNumberOfGrids ; g++ )
   {
      Kp[g].set_to_zero();
      K[g].set_to_zero();
   }
   double t = mDt*(mNumberOfTimeSteps-1);
   int beginCycle = 1;

   double time_measure[8];
   double time_sum[8]={0,0,0,0,0,0,0,0};
   double time_start_solve = MPI_Wtime();

// Save all time series
//	for (int ts=0; ts<a_TimeSeries.size(); ts++)
//	   a_TimeSeries[ts]->writeFile();

   // Backward time stepping loop
   for( int currentTimeStep = mNumberOfTimeSteps ; currentTimeStep >= beginCycle; currentTimeStep-- )
   {    
      time_measure[0] = MPI_Wtime();
      evalRHS( K, a_Mu, a_Lambda, Lk, AlphaVE );
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         F[g].set_to_zero();
      evalPredictor( Km, K, Kp, a_Rho, Lk, F );

      time_measure[1] = MPI_Wtime();

     // Boundary conditions on predictor
      //      enforceDirichlet5( Km );
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( Km[g], g );
      cartesian_bc_forcing( t-mDt, BCForcing, a_Sources );
      enforceBC( Km, a_Mu, a_Lambda, t-mDt, BCForcing );

      time_measure[2] = MPI_Wtime();

    // Corrector
      for( int s= 0 ; s < a_TimeSeries.size() ; s++ )
	 a_TimeSeries[s]->use_as_forcing( currentTimeStep-1, F, mGridSize, mDt, mJ, topographyExists() );

      evalDpDmInTime( Kp, K, Km, Kacc ); 
      evalRHS( Kacc, a_Mu, a_Lambda, Lk, AlphaVEm );
      evalCorrector( Km, a_Rho, Lk, F );

      time_measure[3] = MPI_Wtime();

    // Add in super-grid damping terms
      if (usingSupergrid())
	 addSuperGridDamping( Km, K, Kp, a_Rho );

      time_measure[4] = MPI_Wtime();

     // Boundary conditions on corrector
      //      enforceDirichlet5( Km );
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( Km[g], g );
      //      cartesian_bc_forcing( t-mDt, BCForcing );
      enforceBC( Km, a_Mu, a_Lambda, t-mDt, BCForcing );

      // U-backward solution, predictor
      evalRHS( U, a_Mu, a_Lambda, Lk, AlphaVE );
      Force( t, F, point_sources, identsources );
      evalPredictor( Um, U, Up, a_Rho, Lk, F );
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( Um[g], g );
      cartesian_bc_forcing( t-mDt, BCForcing, a_Sources );
      enforceBC( Um, a_Mu, a_Lambda, t-mDt, BCForcing );
      
     // U-backward solution, corrector
      evalDpDmInTime( Up, U, Um, Uacc ); 

      // set boundary data on Uacc, from forward solver
      for( int g=0 ; g < mNumberOfGrids ; g++ )
	 Upred_saved[g]->pop( Uacc[g], currentTimeStep );
      enforceBC( Uacc, a_Mu, a_Lambda, t, BCForcing );

      evalRHS( Uacc, a_Mu, a_Lambda, Lk, AlphaVEm );
      Force_tt( t, F, point_sources, identsources );
      evalCorrector( Um, a_Rho, Lk, F );

      for( int g=0 ; g < mNumberOfGrids ; g++ )
	 Ucorr_saved[g]->pop( Um[g], currentTimeStep-2 );

      // set boundary data on U
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( Um[g], g );
      cartesian_bc_forcing( t-mDt, BCForcing, a_Sources );
      enforceBC( Um, a_Mu, a_Lambda, t-mDt, BCForcing );

      time_measure[5] = MPI_Wtime();

      // Accumulate the gradient
      for( int s=0 ; s < point_sources.size() ; s++ )
      {
	 point_sources[s]->add_to_gradient( K, Kacc, t, mDt, gradientsrc, mGridSize, mJ, topographyExists() );
	 //	 point_sources[s]->add_to_hessian(  K, Kacc, t, mDt, hessian,  mGridSize );
      }
      add_to_grad( K, Kacc, Um, U, Up, Uacc, gRho, gMu, gLambda );      

      time_measure[6] = MPI_Wtime();
      t -= mDt;
      cycleSolutionArrays( Kp, K, Km );
      cycleSolutionArrays( Up, U, Um );
      time_measure[7] = MPI_Wtime();
      time_sum[0] += time_measure[1]-time_measure[0]; // Predictor
      time_sum[1] += time_measure[2]-time_measure[1]; // Predictor, bc
      time_sum[2] += time_measure[3]-time_measure[2]; // Corrector
      time_sum[3] += time_measure[4]-time_measure[3]; // Super grid damping
      time_sum[4] += time_measure[5]-time_measure[4]; // Corrector, bc
      time_sum[5] += time_measure[6]-time_measure[5]; // Gradient accumulation
      time_sum[6] += time_measure[7]-time_measure[6]; // Cycle arrays
   }
   time_sum[7] = MPI_Wtime() - time_start_solve; // Total solver time
   //   cout << "Final t = " << t << endl;

   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      double upmx[3]={0,0,0}, umx[3]={0,0,0}, tmp[3];
      for( int k=m_kStartAct[g] ; k <= m_kEndAct[g] ; k++ )
	 for( int j=m_jStartAct[g] ; j <= m_jEndAct[g] ; j++ )
	    for( int i=m_iStartAct[g] ; i <= m_iEndAct[g] ; i++ )
	       for( int c=0 ; c < 3 ; c++ )
	       {
		  if( fabs(Up[g](c+1,i,j,k) )>upmx[c] )
		     upmx[c] = fabs(Up[g](c+1,i,j,k));
		  if( fabs( U[g](c+1,i,j,k) )>umx[c] )
		     umx[c] = fabs(U[g](c+1,i,j,k));
	       }
      tmp[0] = upmx[0];
      tmp[1] = upmx[1];
      tmp[2] = upmx[2];
      MPI_Allreduce( tmp, upmx, 3, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
      tmp[0] = umx[0];
      tmp[1] = umx[1];
      tmp[2] = umx[2];
      MPI_Allreduce( tmp, umx, 3, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
      if( m_myRank == 0 )
      {
         cout << "Grid nr. " << g << ": " << endl;
	 cout << "   Max norm of backed out U  = " << upmx[0] <<  " " << upmx[1] <<  " " << upmx[2] << endl;
	 cout << "   Max norm of backed out Um = " << umx[0]  <<  " " << umx[1]  <<  " " << umx[2]  << endl;
      }
   }
   communicate_arrays( gRho );
   communicate_arrays( gMu );
   communicate_arrays( gLambda );
   

   //   gRho[0].save_to_disk("grho.bin");
   //   gMu[0].save_to_disk("gmu.bin");
   //   gLambda[0].save_to_disk("glambda.bin");

    for( int i3 = 0 ; i3 < mImage3DFiles.size() ; i3++ )
       mImage3DFiles[i3]->force_write_image( t, 0, Up, a_Rho, a_Mu, a_Lambda, gRho, gMu, gLambda, mQp, mQs, mPath, mZ );

    for( int i2 = 0 ; i2 < mImageFiles.size() ; i2++ )
    {
       if( mImageFiles[i2]->needs_mgrad() )
	  mImageFiles[i2]->output_image( 0, t, mDt, Up, U, Um, a_Rho, a_Mu, a_Lambda, gRho, gMu, gLambda, a_Sources, 0 );
    }
   //   Up[0].save_to_disk("ubackedout.bin");
   //   U[0].save_to_disk("umbackedout.bin");
   // 
    //    material_to_parameters( nmpar, gradientm, gRho, gMu, gLambda );

   // Sum source gradient contributions from all processors
   double gradtmp[11];
   for( int s=0 ; s < 11 ; s++ )
      gradtmp[s] = gradientsrc[s];
   MPI_Allreduce( gradtmp, gradientsrc, 11, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

   //   // Sum Hessian contributions from all processors
   //   double hesstmp[121];
   //   for( int s=0 ; s < 121 ; s++ )
   //      hesstmp[s] = hessian[s];
   //   MPI_Allreduce( hesstmp, hessian, 121, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

   //   // Symmetry gives the lower half of matrix:
   //   for( int m= 0 ; m < 11 ; m++ )
   //      for( int j=0 ; j<m ; j++ )
   //	 hessian[m+11*j] = hessian[j+11*m];

   // Give back memory
   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      for(int side=0; side < 6; side++)
	 if( BCForcing[g][side] != NULL )
	    delete[] BCForcing[g][side];
      delete[] BCForcing[g];
   }
   for( int s=0 ; s < point_sources.size() ; s++ )
      delete point_sources[s];
}

//-----------------------------------------------------------------------
void EW::enforceDirichlet5( vector<Sarray> & a_U )
{
  for(int g=0 ; g<mNumberOfGrids; g++ )
  {
    double* u_ptr = a_U[g].c_ptr();
    int ifirst, ilast, jfirst, jlast, kfirst, klast;
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];

    int iafirst, ialast, jafirst, jalast, kafirst, kalast;
    iafirst = m_iStartAct[g];
    ialast  = m_iEndAct[g];
    jafirst = m_jStartAct[g];
    jalast  = m_jEndAct[g];
    kafirst = m_kStartAct[g];
    kalast  = m_kEndAct[g];
    for( int k=kalast+1 ; k <= klast ; k++ )
       for( int j=jfirst ; j <= jlast ; j++ )
	  for( int i=ifirst ; i <= ilast ; i++ )
	     a_U[g](1,i,j,k)=a_U[g](2,i,j,k)=a_U[g](3,i,j,k)=0;
    for( int k=kfirst ; k <= klast ; k++ )
    {
       for( int j=jfirst ; j <= jafirst-1 ; j++ )
	  for( int i=ifirst ; i <= ilast ; i++ )
	     a_U[g](1,i,j,k)=a_U[g](2,i,j,k)=a_U[g](3,i,j,k)=0;
       for( int j=jalast+1 ; j <= jlast ; j++ )
	  for( int i=ifirst ; i <= ilast ; i++ )
	     a_U[g](1,i,j,k)=a_U[g](2,i,j,k)=a_U[g](3,i,j,k)=0;
       for( int j=jfirst ; j <= jlast ; j++ )
	  for( int i=ifirst ; i <= iafirst-1 ; i++ )
	     a_U[g](1,i,j,k)=a_U[g](2,i,j,k)=a_U[g](3,i,j,k)=0;
       for( int j=jfirst ; j <= jlast ; j++ )
	  for( int i=ialast+1 ; i <= ilast ; i++ )
	     a_U[g](1,i,j,k)=a_U[g](2,i,j,k)=a_U[g](3,i,j,k)=0;
    }
    //    F77_FUNC(hdirichlet5,HDIRICHLET5)( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
    //				       &klast, &iafirst, &ialast, &jafirst, &jalast,
    //				       &kafirst, &kalast, u_ptr );
  }
}
