#include "EW.h"

void EW::solve_dudp( vector<Source*>& a_Sources, vector<Sarray>& a_Rho, 
                     vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                     vector<TimeSeries*> & a_TimeSeries, 
                     vector<TimeSeries*> & GlobalObservations,
                     vector<Sarray>& Um, vector<Sarray>& U,
                     vector<Sarray>& dUm, vector<Sarray>& dU,
                     double& misfit, double& dmisfitdp, 
                     int di, int dj, int dk, int dgrid, int event )
{
   const float_sw4 dt2   = mDt*mDt;
   const float_sw4 todt4 = 12/(mDt*mDt*mDt*mDt);
   int ng=mNumberOfGrids;
   vector<Sarray> F(ng), Lu(ng), Uacc(ng), Up(ng), dUp(ng), G(ng);
   vector<Sarray*> AlphaVE(ng), AlphaVEm(ng), AlphaVEp(ng);
   vector<float_sw4**> BCForcing(ng);
   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      F[g].define(U[g]);
      Lu[g].define(U[g]);
      Uacc[g].define(U[g]);
      Up[g].define(U[g]);
      dUp[g].define(U[g]);
      G[g].define(U[g]);
      F[g].set_to_zero();
      G[g].set_to_zero();
      U[g].set_to_zero();
      Um[g].set_to_zero();
      dU[g].set_to_zero();
      dUm[g].set_to_zero();
   }
   bool duforcing=interior_point_in_proc(di,dj,dgrid);

   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      BCForcing[g] = new float_sw4*[6];
      for (int side=0; side < 6; side++)
      {
	 BCForcing[g][side]=NULL;
	 if (m_bcType[g][side] == bStressFree || 
             m_bcType[g][side] == bDirichlet || 
             m_bcType[g][side] == bSuperGrid)
	 {
	    BCForcing[g][side] = new float_sw4[3*m_NumberOfBCPoints[g][side]];
	 }
      }
   }

   vector<TimeSeries*> dudpTimeSeries;
   for(int ts=0; ts<a_TimeSeries.size(); ts++)
      dudpTimeSeries.push_back( a_TimeSeries[ts]->copy( this, "dudp" ) );
   for(int ts=0; ts<a_TimeSeries.size(); ts++)
   {
     a_TimeSeries[ts]->allocateRecordingArrays( mNumberOfTimeSteps[event]+1, mTstart, mDt);
     dudpTimeSeries[ts]->allocateRecordingArrays( mNumberOfTimeSteps[event]+1, mTstart, mDt);
   }

   vector<GridPointSource*> point_sources;
   for( unsigned int i=0 ; i < a_Sources.size() ; i++ )
      a_Sources[i]->set_grid_point_sources4( this, point_sources );
   vector<int> identsources;
   sort_grid_point_sources( point_sources, identsources );

   double t = mTstart;
   for (int ts=0; ts<a_TimeSeries.size(); ts++)
   {
      if (a_TimeSeries[ts]->getMode() != TimeSeries::Velocity && a_TimeSeries[ts]->myPoint())
      {
         vector<float_sw4> uRec;
         int i0 = a_TimeSeries[ts]->m_i0;
         int j0 = a_TimeSeries[ts]->m_j0;
         int k0 = a_TimeSeries[ts]->m_k0;
         int grid0 = a_TimeSeries[ts]->m_grid0;
         extractRecordData(a_TimeSeries[ts]->getMode(), i0, j0, k0, grid0, 
                           uRec, Um, U); 
         a_TimeSeries[ts]->recordData(uRec);
         extractRecordData(dudpTimeSeries[ts]->getMode(), i0, j0, k0, grid0, 
                           uRec, dUm, dU); 
         dudpTimeSeries[ts]->recordData(uRec);
      }
   }

   Force( t, F, point_sources, identsources );

   for( int n = 0; n <= mNumberOfTimeSteps[event]; n++)
   {    
// evaluate spatial discretization
      evalRHS( U, a_Mu, a_Lambda, Lu, AlphaVE ); // save Lu in composite grid 'Lu'

// take predictor step, store in Up
      evalPredictor( Up, U, Um, a_Rho, Lu, F );

// communicate across processor boundaries
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( Up[g], g );

// calculate boundary forcing at time t+mDt
      cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );

// update ghost points in Up
      enforceBC( Up, a_Mu, a_Lambda, AlphaVEp, t+mDt, BCForcing );

// precompute F_tt(t)
      Force_tt( t, F, point_sources, identsources );

// *** 4th order in TIME interface conditions for the predictor
      enforceIC( Up, U, Um, AlphaVEp, AlphaVE, AlphaVEm, t, true, F, point_sources );


// corrector step for
      evalDpDmInTime( Up, U, Um, Uacc ); // store result in Uacc

      evalRHS( Uacc, a_Mu, a_Lambda, Lu, AlphaVEm );
      evalCorrector( Up, a_Rho, Lu, F );

// add in super-grid damping terms
      if (usingSupergrid())
         addSuperGridDamping( Up, U, Um, a_Rho );

// communicate across processor boundaries
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( Up[g], g );

// calculate boundary forcing at time t+mDt 
      cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );

// update ghost points in Up
      enforceBC( Up, a_Mu, a_Lambda, AlphaVEp, t+mDt, BCForcing );

// compute forcing for next time step here so it can be used in enforceIC()
      Force( t+mDt, F, point_sources, identsources );

      enforceIC( Up, U, Um, AlphaVEp, AlphaVE, AlphaVEm, t, false, F, point_sources );


// Derivative computation:
      evalRHS( dU, a_Mu, a_Lambda, Lu, AlphaVE ); // save Lu in composite grid 'Lu'

      if( duforcing )
      {
         G[dgrid](1,di,dj,dk)=-Uacc[dgrid](1,di,dj,dk);
         G[dgrid](2,di,dj,dk)=-Uacc[dgrid](2,di,dj,dk);
         G[dgrid](3,di,dj,dk)=-Uacc[dgrid](3,di,dj,dk);
      }
// take predictor step, store in Up
      evalPredictor( dUp, dU, dUm, a_Rho, Lu, G );

// communicate across processor boundaries
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( dUp[g], g );

// calculate boundary forcing at time t+mDt
      cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );

// update ghost points in Up
      enforceBC( dUp, a_Mu, a_Lambda, AlphaVEp, t+mDt, BCForcing );

// precompute F_tt(t)
//      Force_tt( t, F, point_sources, identsources );
      if( duforcing )
      {
         G[dgrid](1,di,dj,dk)=-todt4*(Up[dgrid](1,di,dj,dk)-dt2*Uacc[dgrid](1,di,dj,dk)-
                                  2*U[dgrid](1,di,dj,dk)+Um[dgrid](1,di,dj,dk));
         G[dgrid](2,di,dj,dk)=-todt4*(Up[dgrid](2,di,dj,dk)-dt2*Uacc[dgrid](2,di,dj,dk)-
                                  2*U[dgrid](2,di,dj,dk)+Um[dgrid](2,di,dj,dk));
         G[dgrid](3,di,dj,dk)=-todt4*(Up[dgrid](3,di,dj,dk)-dt2*Uacc[dgrid](3,di,dj,dk)-
                                  2*U[dgrid](3,di,dj,dk)+Um[dgrid](3,di,dj,dk));
      }
// *** 4th order in TIME interface conditions for the predictor
      enforceIC( dUp, dU, dUm, AlphaVEp, AlphaVE, AlphaVEm, t, true, G, point_sources );

// corrector step for
      evalDpDmInTime( dUp, dU, dUm, Uacc ); // store result in Uacc

      evalRHS( Uacc, a_Mu, a_Lambda, Lu, AlphaVEm );
      evalCorrector( dUp, a_Rho, Lu, G );
      if (usingSupergrid())
         addSuperGridDamping( dUp, dU, dUm, a_Rho );

// communicate across processor boundaries
      for(int g=0 ; g < mNumberOfGrids ; g++ )
         communicate_array( dUp[g], g );

// calculate boundary forcing at time t+mDt 
      cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );

// update ghost points in Up
      enforceBC( dUp, a_Mu, a_Lambda, AlphaVEp, t+mDt, BCForcing );

// compute forcing for next time step here so it can be used in enforceIC()
//      Force( t+mDt, F, point_sources, identsources );
      if( duforcing )
      {
         evalLupt( Up, a_Mu, a_Lambda, Lu, dgrid, di, dj, dk );
         G[dgrid](1,di,dj,dk) = -(Lu[dgrid](1,di,dj,dk)+F[dgrid](1,di,dj,dk))/a_Rho[dgrid](di,dj,dk);
         G[dgrid](2,di,dj,dk) = -(Lu[dgrid](2,di,dj,dk)+F[dgrid](2,di,dj,dk))/a_Rho[dgrid](di,dj,dk);
         G[dgrid](3,di,dj,dk) = -(Lu[dgrid](3,di,dj,dk)+F[dgrid](3,di,dj,dk))/a_Rho[dgrid](di,dj,dk);
      }

      enforceIC( dUp, dU, dUm, AlphaVEp, AlphaVE, AlphaVEm, t, false, G, point_sources );
       
// increment time
       t += mDt;

// periodically, print time stepping info to stdout
       printTime( n, t, n == mNumberOfTimeSteps[event] ); 

// save the current solution on receiver records (time-derivative require Up and Um 
// for a 2nd order approximation, so do this before cycling the arrays).
       for (int ts=0; ts<a_TimeSeries.size(); ts++)
       {
          if (a_TimeSeries[ts]->myPoint())
          {
             vector<float_sw4> uRec;
             int i0 = a_TimeSeries[ts]->m_i0;
             int j0 = a_TimeSeries[ts]->m_j0;
             int k0 = a_TimeSeries[ts]->m_k0;
             int grid0 = a_TimeSeries[ts]->m_grid0;

// note that the solution on the new time step is in Up
// also note that all quantities related to velocities lag by one time step; they are not
// saved before the time stepping loop started
             extractRecordData(a_TimeSeries[ts]->getMode(), i0, j0, k0, grid0, 
                               uRec, Um, Up);
             a_TimeSeries[ts]->recordData(uRec);

             extractRecordData(dudpTimeSeries[ts]->getMode(), i0, j0, k0, grid0, 
                               uRec, dUm, dUp);
             dudpTimeSeries[ts]->recordData(uRec);
          }
       }
       cycleSolutionArrays( Um,  U,  Up, AlphaVEm, AlphaVE, AlphaVEp);
       cycleSolutionArrays(dUm, dU, dUp, AlphaVEm, AlphaVE, AlphaVEp);
   }
   // End of time stepping.
   // Now compute misfit and its derivative
   misfit  = dmisfitdp = 0;
   float_sw4 misfitlocal, dmisfitlocal;
   for (int ts=0; ts<a_TimeSeries.size(); ts++)
   {
      a_TimeSeries[ts]->misfitanddudp( GlobalObservations[ts], dudpTimeSeries[ts],
                                       misfitlocal, dmisfitlocal );
      misfit    += misfitlocal;
      dmisfitdp += dmisfitlocal;
   }
   misfitlocal = misfit;
   MPI_Allreduce( &misfitlocal, &misfit, 1, MPI_DOUBLE, MPI_SUM, m_1d_communicator );
   dmisfitlocal = dmisfitdp;
   MPI_Allreduce( &dmisfitlocal, &dmisfitdp, 1, MPI_DOUBLE, MPI_SUM, m_1d_communicator );

   for(int ts=0; ts<dudpTimeSeries.size(); ts++)
      delete dudpTimeSeries[ts];

}
