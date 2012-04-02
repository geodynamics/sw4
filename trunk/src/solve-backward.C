#include "EW.h"

void EW::solve_backward( vector<Source*> & a_Sources, vector<TimeSeries*> & a_TimeSeries, double gradient[11],
			 double hessian[121] )
{
// solution arrays
   vector<Sarray> F, Lk, Kacc, Kp, Km, K;
   vector<double **> BCForcing;
 
   F.resize(mNumberOfGrids);
   Lk.resize(mNumberOfGrids);
   Kacc.resize(mNumberOfGrids);
   Kp.resize(mNumberOfGrids);
   Km.resize(mNumberOfGrids);
   K.resize(mNumberOfGrids);
   BCForcing.resize(mNumberOfGrids);

   int ifirst, ilast, jfirst, jlast, kfirst, klast;
   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      BCForcing[g] = new double *[6];
      for(int side=0; side < 6; side++)
      {
	 BCForcing[g][side]=NULL;
         if (m_bcType[g][side] == bStressFree || m_bcType[g][side] == bDirichlet || m_bcType[g][side] == bSuperGrid)
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
   }


// will accumulate the gradient of the misfit in this array
   for( int s=0 ; s < 11 ; s++ )
      gradient[s] = 0;

// will accumulate the Hessian of the misfit in this array
   for( int s=0 ; s < 121 ; s++ )
      hessian[s] = 0;
   
// the Source objects get discretized into GridPointSource objects
   vector<GridPointSource*> point_sources;

// Transfer source terms to each individual grid as point sources at grid points.
   for( unsigned int i=0 ; i < a_Sources.size() ; i++ )
      if (!a_Sources[i]->ignore())
	 a_Sources[i]->set_grid_point_sources4( this, point_sources );


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

   // Backward time stepping loop
   for( int currentTimeStep = mNumberOfTimeSteps ; currentTimeStep >= beginCycle; currentTimeStep-- )
   {    
      time_measure[0] = MPI_Wtime();
      evalRHS( K, Lk );
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         F[g].set_to_zero();
      evalPredictor( Km, K, Kp, Lk, F );

      time_measure[1] = MPI_Wtime();

     // Boundary conditions on predictor
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( Km[g], g );
      cartesian_bc_forcing( t-mDt, BCForcing );
      enforceBC( Km, t-mDt, BCForcing );

      time_measure[2] = MPI_Wtime();

    // Corrector
      for( int s= 0 ; s < a_TimeSeries.size() ; s++ )
	 a_TimeSeries[s]->use_as_forcing( currentTimeStep-1, F, mGridSize, mDt );

      evalDpDmInTime( Kp, K, Km, Kacc ); 
      evalRHS( Kacc, Lk );
      evalCorrector( Km, Lk, F );

      time_measure[3] = MPI_Wtime();
    // Add in super-grid damping terms
      if (usingSupergrid())
	 addSuperGridDamping( Km, K, Kp );

      time_measure[4] = MPI_Wtime();
     // Boundary conditions on corrector
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( Km[g], g );
      //      cartesian_bc_forcing( t-mDt, BCForcing );
      enforceBC( Km, t-mDt, BCForcing );

      time_measure[5] = MPI_Wtime();
      // Accumulate the gradient
      for( int s=0 ; s < point_sources.size() ; s++ )
      {
	 point_sources[s]->add_to_gradient( K, Kacc, t, mDt, gradient, mGridSize );
	 point_sources[s]->add_to_hessian( K, Kacc, t, mDt, hessian, mGridSize );
      }
      
      time_measure[6] = MPI_Wtime();
      t -= mDt;
      cycleSolutionArrays( Kp, K, Km );
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

   // 
   // Sum gradient contributions from all processors
   double gradtmp[11];
   for( int s=0 ; s < 11 ; s++ )
      gradtmp[s] = gradient[s];
   MPI_Allreduce( gradtmp, gradient, 11, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

   // Sum Hessian contributions from all processors
   double hesstmp[121];
   for( int s=0 ; s < 121 ; s++ )
      hesstmp[s] = hessian[s];
   MPI_Allreduce( hesstmp, hessian, 121, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

   // Symmetry gives the lower half of matrix:
   for( int m= 0 ; m < 11 ; m++ )
      for( int j=0 ; j<m ; j++ )
	 hessian[m+11*j] = hessian[j+11*m];
}

//------------------------------------------------------------------------
void EW::cycleSolutionArrays(vector<Sarray> & a_Um, vector<Sarray> & a_U, vector<Sarray> & a_Up ) 
{
  for (int g=0; g<mNumberOfGrids; g++)
  {
    double *tmp = a_Um[g].c_ptr();
    a_Um[g].reference(a_U[g].c_ptr());
    a_U[g].reference(a_Up[g].c_ptr());
    a_Up[g].reference(tmp);
  }
}
