//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 

#include "EW.h"
#include "impose_cartesian_bc.h"
#include "cf_interface.h"
#include "f_interface.h"
#include "Mspace.h"
#include "policies.h"
#include "caliper.h"
#ifdef ENABLE_CUDA
#include "cuda_profiler_api.h"
#endif

#define SQR(x) ((x)*(x))

//--------------------------------------------------------------------
void EW::solve( vector<Source*> & a_Sources, vector<TimeSeries*> & a_TimeSeries )
{
  SW4_MARK_FUNCTION;
// solution arrays
  vector<Sarray> F, Lu, Uacc, Up, Um, U;
  vector<Sarray*> AlphaVE, AlphaVEm, AlphaVEp;
// vectors of pointers to hold boundary forcing arrays in each grid
  vector<float_sw4 **> BCForcing;

  BCForcing.resize(mNumberOfGrids);
  F.resize(mNumberOfGrids);
  Lu.resize(mNumberOfGrids);
  Uacc.resize(mNumberOfGrids);
  Up.resize(mNumberOfGrids);
  Um.resize(mNumberOfGrids);
  U.resize(mNumberOfGrids);

// Allocate pointers, even if attenuation not used, to avoid segfault in parameter list with mMuVE[g], etc...
  AlphaVE.resize(mNumberOfGrids);
  AlphaVEm.resize(mNumberOfGrids);
  AlphaVEp.resize(mNumberOfGrids);
  if (m_use_attenuation && m_number_mechanisms > 0)
  {
    for( int g = 0; g <mNumberOfGrids; g++ )
    {
      AlphaVE[g]  = new Sarray[m_number_mechanisms];
      AlphaVEp[g] = new Sarray[m_number_mechanisms];
      AlphaVEm[g] = new Sarray[m_number_mechanisms];
    }
  }

  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  for( int g = 0; g <mNumberOfGrids; g++ )
  {
    BCForcing[g] = new float_sw4 *[6];
    for (int side=0; side < 6; side++)
    {
      BCForcing[g][side]=NULL;
      if (m_bcType[g][side] == bStressFree || m_bcType[g][side] == bDirichlet || m_bcType[g][side] == bSuperGrid)
      {
    	BCForcing[g][side] = SW4_NEW(Managed,float_sw4[3*m_NumberOfBCPoints[g][side]]);
      }
      
    }

    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];

    kfirst = m_kStart[g];
    klast = m_kEnd[g];

    F[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
    Lu[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
    Uacc[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
    Up[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
    Um[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
    U[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
    if (m_use_attenuation && m_number_mechanisms > 0)
    {
      for (int a=0; a<m_number_mechanisms; a++)
      {
	AlphaVE[g][a].define( 3,ifirst,ilast,jfirst,jlast,kfirst,klast);
	AlphaVEp[g][a].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
	AlphaVEm[g][a].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      }
    }
  }
// done allocating solution arrays

// Set the number of time steps, allocate the recording arrays, and set reference time in all time series objects  
#pragma omp parallel for
  for (int ts=0; ts<a_TimeSeries.size(); ts++)
  {
     a_TimeSeries[ts]->allocateRecordingArrays( mNumberOfTimeSteps+1, mTstart, mDt); // AP: added one to mNumber...
     // In forward solve, the output receivers will use the same UTC as the
     // global reference utc0, therefore, set station utc equal reference utc.
     //     if( m_utc0set )
     //	a_TimeSeries[ts]->set_station_utc( m_utc0 );
  }
  if( !mQuiet && mVerbose >=3 && proc_zero() )
    printf("***  Allocated all receiver time series\n");

// Reset image time to zero, in case we are rerunning the solver
#pragma omp parallel for
  for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
     mImageFiles[fIndex]->initializeTime();
   
// the Source objects get discretized into GridPointSource objects
  vector<GridPointSource*> point_sources;

// Transfer source terms to each individual grid as point sources at grid points.
  SW4_MARK_BEGIN("set_grid_point_sources4");
  for( unsigned int i=0 ; i < a_Sources.size() ; i++ )
      a_Sources[i]->set_grid_point_sources4( this, point_sources );
   SW4_MARK_END("set_grid_point_sources4");
 // Debug
  // for (int proc = 0; proc<m_nProcs; proc++)
  //    if (proc == m_myRank)
  //    {
  //       int nSources=0;
  //       for( unsigned int i=0 ; i < a_Sources.size() ; i++ )
  //       {
  //          if (a_Sources[i]->m_timeFuncIsReady) nSources++;
  //       }
  //       printf("\n**** MPI-task #%d needs %d source terms  ********\n\n", proc, nSources);     
  // }
// end debug
  
// modification of time functions by prefiltering is currently done in preprocessSources()
  // only reported here
  if (!m_testing && m_prefilter_sources)
  {
    if (!mQuiet && proc_zero() )
    {
      if (m_filter_ptr->get_type()==lowPass)
	printf("Lowpass filtering all source time functions to corner frequency fc2=%e\n", 
	       m_filter_ptr->get_corner_freq2());
      else if (m_filter_ptr->get_type()==bandPass)
	printf("Bandpass filtering all source time functions to corner frequencies fc1=%e and fc2=%e\n", 
	       m_filter_ptr->get_corner_freq1(), m_filter_ptr->get_corner_freq2());
    }
    
// tmp
    // if (proc_zero() && point_sources.size()>0)
    // {
    //   printf("Saving one un-filtered original time function\n");
	 
    //   FILE *tf=fopen("g0.dat","w");
    //   float_sw4 t;
    //   float_sw4 gt;
    //   for (int i=0; i<=mNumberOfTimeSteps; i++)
    //   {
    // 	t = mTstart + i*mDt;
    // 	gt = point_sources[0]->getTimeFunc(t);
    // 	fprintf(tf, "%e %.18e\n", t, gt);
    //   }
    //   fclose(tf);
    // }

// 3. Replace the time function by a filtered one, represented by a (long) vector holding values at each time step   
//    for( int s=0; s < point_sources.size(); s++ ) 
//      point_sources[s]->discretizeTimeFuncAndFilter(mTstart, mDt, mNumberOfTimeSteps, m_filter_ptr);

// tmp
//    if (proc_zero() && point_sources.size()>0)
//    {
//      printf("Saving one filtered discretized time function\n");
//	 
//      FILE *tf=fopen("g1.dat","w");
//      float_sw4 t;
//      float_sw4 gt, gt1, gt2;
//      for (int i=0; i<=mNumberOfTimeSteps; i++)
//      {
//    	t = mTstart + i*mDt;
//    	gt = point_sources[0]->getTimeFunc(t);
//    	gt1 = point_sources[0]->evalTimeFunc_t(t);
//    	gt2 = point_sources[0]->evalTimeFunc_tt(t);
//    	fprintf(tf, "%e  %.18e  %.18e  %.18e\n", t, gt, gt1, gt2);
//      }
//      fclose(tf);
//    }
       
  } // end if prefiltering

// AP changed to false
  bool output_timefunc = false;
  if( output_timefunc )
  {
     int has_source_id=-1, has_source_max;
     if( point_sources.size() > 0 )
	has_source_id = m_myRank;

     MPI_Allreduce( &has_source_id, &has_source_max, 1, MPI_INT, MPI_MAX, m_cartesian_communicator );
     if( m_myRank == has_source_max )
     {
       if (!mQuiet && mVerbose >=1 )
	 printf("*** Saving one discretized time function ***\n");

// tmp
       // printf("mTstart = %e, mDt = %e\n", mTstart, mDt);
       // printf("GridPointSource::mT0 = %e\n", point_sources[0]->mT0);
// end tmp
       
//building the file name...
       string filename;
       if( mPath != "." )
	 filename += mPath;
       filename += "g1.dat";	 

       FILE *tf=fopen(filename.c_str(),"w");
       float_sw4 t;
       float_sw4 gt, gt1, gt2;
       for (int i=0; i<=mNumberOfTimeSteps; i++)
       {
	 //           for( int sb=0 ; sb < 10 ; sb++ )
	 //	   {
	 //	      t = mTstart + i*mDt + 0.1*sb*mDt;
	 t = mTstart + i*mDt;
	 gt = point_sources[0]->getTimeFunc(t);
	 gt1 = point_sources[0]->evalTimeFunc_t(t);
	 gt2 = point_sources[0]->evalTimeFunc_tt(t);
	 fprintf(tf, "%.18e  %.18e  %.18e  %.18e\n", t, gt, gt1, gt2);
	 //	   }
       }
       fclose(tf);
     }
  }

  if( !mQuiet && mVerbose && proc_zero() )
  {
    cout << endl << "***  Starting solve ***" << endl;
  }
  printPreamble(a_Sources);

  
// Set up timers
  double time_start_solve = MPI_Wtime();
  double time_measure[20];
  double time_sum[9]={0,0,0,0,0,0,0,0,0};
  double bc_time_measure[5]={0,0,0,0,0};

// Sort sources wrt spatial location, needed for thread parallel computing
  vector<int> identsources;
  sort_grid_point_sources( point_sources, identsources );

// Assign initial data
  int beginCycle; 
  float_sw4 t;
  if( m_check_point->do_restart() )
  {
     double timeRestartBegin = MPI_Wtime();
     m_check_point->read_checkpoint( t, beginCycle, Um, U,
				     AlphaVEm, AlphaVE );
// tmp
     if (proc_zero())
        printf("After reading checkpoint data: beginCycle=%d, t=%e\n", beginCycle, t);
// end tmp     

     // Make sure the TimeSeries output has the correct time shift,
     // and know's it's a restart
     double timeSeriesRestartBegin = MPI_Wtime();
     for (int ts=0; ts<a_TimeSeries.size(); ts++)
     {
       a_TimeSeries[ts]->doRestart(this, false, t, beginCycle);
     }
     double timeSeriesRestart = MPI_Wtime() - timeSeriesRestartBegin;
     if( proc_zero() && m_output_detailed_timing )
     {
        cout << "Wallclock time to read checkpoint file: " << timeSeriesRestartBegin-timeRestartBegin << " seconds " << endl;
        cout << "Wallclock time to read " << a_TimeSeries.size() << " sets of station files: " << timeSeriesRestart << " seconds " << endl;
     }

// Reset image time to the time corresponding to restart
#pragma omp parallel for
  for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
     mImageFiles[fIndex]->initializeTime(t);
   


     // Restart data is defined at ghost point outside physical boundaries, still
     // need to communicate solution arrays to define it a parallel overlap points
     beginCycle++; // needs to be one step ahead of 't', see comment 5 lines below
  }
  else
  {
// NOTE: time stepping loop starts at currentTimeStep = beginCycle; ends at currentTimeStep <= mNumberOfTimeSteps
// However, the time variable 't' is incremented at the end of the time stepping loop. Thus the time step index is one step
// ahead of 't' at the start.
     beginCycle = 1; 
     t = mTstart;
     initialData(mTstart, U, AlphaVE);
     initialData(mTstart-mDt, Um, AlphaVEm );
  }
  
  if ( !mQuiet && mVerbose && proc_zero() )
    cout << "  Initial data has been assigned" << endl;

// // Assign Up to make it possible to output velDiv and velCurl images of the initial data
//   for( g=0 ; g<mNumberOfCartesianGrids; g++ )
//   {
//     m_forcing->get_initial_data_Cartesian( m_zmin[g], mGridSize[g], mTstart+mDt, Up[g] );
//   }
  
//   if ( topographyExists() )
//   {
//     g = mNumberOfGrids-1;
//     m_forcing->get_initial_data_Curvilinear( mX, mY, mZ, mTstart+mDt, Up[g] );
//   }
  
// moved below, after enforcing BC
// // save any images for cycle = 0 (initial data) ?
//   update_images( 0, t, U, Um, Up, mRho, mMu, mLambda, a_Sources, 1 );
//   for( int i3 = 0 ; i3 < mImage3DFiles.size() ; i3++ )
//     mImage3DFiles[i3]->update_image( t, 0, mDt, U, mRho, mMu, mLambda, mRho, mMu, mLambda, mQp, mQs, mPath, mZ );

// do some testing...
  if (m_twilight_forcing && getVerbosity() >= 3) // only do these tests if verbose>=3
  {
    if ( !mQuiet && proc_zero() )
      cout << "***Twilight Testing..." << endl;

// output some internal flags        
    for(int g=0; g<mNumberOfGrids; g++)
    {
      printf("proc=%i, Onesided[grid=%i]:", m_myRank, g);
      for (int q=0; q<6; q++)
	printf(" os[%i]=%i", q, m_onesided[g][q]);
      printf("\n");
      printf("proc=%i, bcType[grid=%i]:", m_myRank, g);
      for (int q=0; q<6; q++)
	printf(" bc[%i]=%i", q, m_bcType[g][q]);
      printf("\n");
    }

// test accuracy of spatial approximation
    if ( proc_zero() )
      printf("\n Testing the accuracy of the spatial difference approximation\n");
    exactRhsTwilight(t, F);
    evalRHS( U, mMu, mLambda, Up, AlphaVE ); // save Lu in composite grid 'Up'


// evaluate and print errors
    float_sw4 * lowZ = new float_sw4[3*mNumberOfGrids];
    float_sw4 * interiorZ = new float_sw4[3*mNumberOfGrids];
    float_sw4 * highZ = new float_sw4[3*mNumberOfGrids];
    //    float_sw4 lowZ[3], interiorZ[3], highZ[3];
    bndryInteriorDifference( F, Up, lowZ, interiorZ, highZ );

    float_sw4* tmp= new float_sw4[3*mNumberOfGrids];
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = lowZ[i];
    MPI_Reduce( tmp, lowZ, 3*mNumberOfGrids, m_mpifloat, MPI_MAX, 0, m_cartesian_communicator );
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = interiorZ[i];
    MPI_Reduce( tmp, interiorZ, 3*mNumberOfGrids, m_mpifloat, MPI_MAX, 0, m_cartesian_communicator );
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = highZ[i];
    MPI_Reduce( tmp, highZ, 3*mNumberOfGrids, m_mpifloat, MPI_MAX, 0, m_cartesian_communicator );

    if ( proc_zero() )
    {
      for( int g=0 ; g < mNumberOfGrids ; g++ )
      {
	printf("Grid nr: %3i \n", g );
	printf("Max errors low-k boundary RHS:  %15.7e  %15.7e  %15.7e\n", lowZ[3*g], lowZ[3*g+1], lowZ[3*g+2]);
	printf("Max errors interior RHS:        %15.7e  %15.7e  %15.7e\n", interiorZ[3*g], interiorZ[3*g+1], interiorZ[3*g+2]);
	printf("Max errors high-k boundary RHS: %15.7e  %15.7e  %15.7e\n", highZ[3*g], highZ[3*g+1], highZ[3*g+2]);
      }
    }
  
// test accuracy of forcing
    evalRHS( U, mMu, mLambda, Lu, AlphaVE ); // save Lu in composite grid 'Lu'
    Force( t, F, point_sources, identsources );
    exactAccTwilight( t, Uacc ); // save Utt in Uacc
    test_RhoUtt_Lu( Uacc, Lu, F, lowZ, interiorZ, highZ );

    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = lowZ[i];
    MPI_Reduce( tmp, lowZ, 3*mNumberOfGrids, m_mpifloat, MPI_MAX, 0, m_cartesian_communicator );
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = interiorZ[i];
    MPI_Reduce( tmp, interiorZ, 3*mNumberOfGrids, m_mpifloat, MPI_MAX, 0, m_cartesian_communicator );
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = highZ[i];
    MPI_Reduce( tmp, highZ, 3*mNumberOfGrids, m_mpifloat, MPI_MAX, 0, m_cartesian_communicator );

    if ( proc_zero() )
    {
      printf("Testing accuracy of rho*utt - L(u) = F\n");
      for( int g=0 ; g < mNumberOfGrids ; g++ )
      {
	printf("Grid nr: %3i \n", g );
	printf("Max errors low-k boundary RHS:  %15.7e  %15.7e  %15.7e\n",lowZ[3*g],lowZ[3*g+1],lowZ[3*g+2]);
	printf("Max errors interior RHS:        %15.7e  %15.7e  %15.7e\n",interiorZ[3*g],interiorZ[3*g+1],interiorZ[3*g+2]);
	printf("Max errors high-k boundary RHS: %15.7e  %15.7e  %15.7e\n",highZ[3*g],highZ[3*g+1],highZ[3*g+2]);
      }
    }
    delete[] tmp;
    delete[] lowZ;
    delete[] interiorZ;
    delete[] highZ;
  } // end m_twilight_forcing    

// after checkpoint restart, we must communicate the memory variables
  if(  m_check_point->do_restart() && m_use_attenuation && (m_number_mechanisms > 0) )
  {
// AlphaVE
// communicate across processor boundaries
     for(int g=0 ; g < mNumberOfGrids ; g++ )
     {
        for(int m=0 ; m < m_number_mechanisms; m++ )
           communicate_array( AlphaVE[g][m], g );
     }
// AlphaVEm
// communicate across processor boundaries
     for(int g=0 ; g < mNumberOfGrids ; g++ )
     {
        for(int m=0 ; m < m_number_mechanisms; m++ )
           communicate_array( AlphaVEm[g][m], g );
     }
  } // end if checkpoint restarting
  

// enforce bc on initial data
// U
// communicate across processor boundaries
  for(int g=0 ; g < mNumberOfGrids ; g++ )
    communicate_array( U[g], g );

  //    U[0].save_to_disk("u-dbg0.bin");
  //    U[1].save_to_disk("u-dbg1.bin");

// boundary forcing
  cartesian_bc_forcing( t, BCForcing, a_Sources );
// OLD
//  if( m_use_attenuation && m_number_mechanisms > 0 )
//     addAttToFreeBcForcing( AlphaVE, BCForcing, m_sbop );
// enforce boundary condition
  if( m_anisotropic )
     enforceBCanisotropic( U, mC, t, BCForcing );
  else
     enforceBC( U, mMu, mLambda, t, BCForcing );   
// Impose un-coupled free surface boundary condition with visco-elastic terms for 'Up'
  if( m_use_attenuation && (m_number_mechanisms > 0) )
  {
     enforceBCfreeAtt2( U, mMu, mLambda, AlphaVE, BCForcing );
  }

  //    U[0].save_to_disk("u-dbg0-bc.bin");
  //    U[1].save_to_disk("u-dbg1-bc.bin");

// Um
// communicate across processor boundaries
  for(int g=0 ; g < mNumberOfGrids ; g++ )
    communicate_array( Um[g], g );

  //    Um[0].save_to_disk("um-dbg0.bin");
  //    Um[1].save_to_disk("um-dbg1.bin");

// boundary forcing
  cartesian_bc_forcing( t-mDt, BCForcing, a_Sources );
// OLD
// if( m_use_attenuation && m_number_mechanisms > 0 )
//    addAttToFreeBcForcing( AlphaVEm, BCForcing, m_sbop );

// enforce boundary condition
  if( m_anisotropic )
     enforceBCanisotropic( Um, mC, t-mDt, BCForcing );
  else
     enforceBC( Um, mMu, mLambda, t-mDt, BCForcing );
// Impose un-coupled free surface boundary condition with visco-elastic terms for 'Up'
  if( m_use_attenuation && (m_number_mechanisms > 0) )
  {
     enforceBCfreeAtt2( Um, mMu, mLambda, AlphaVEm, BCForcing );
  }


// more testing
  if (m_twilight_forcing && m_check_point->do_restart() && getVerbosity()>=3)
  {
    if ( proc_zero() )
      printf("Checking the accuracy of the checkpoint data\n");

// check the accuracy of the initial data, store exact solution in Up, ignore AlphaVE
    float_sw4 errInf=0, errL2=0, solInf=0, solL2=0;
    exactSol( t, Up, AlphaVEp, a_Sources );

    normOfDifference( Up, U, errInf, errL2, solInf, a_Sources );

    if ( proc_zero() )
       printf("\n Checkpoint errors in U: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);

    if( m_use_attenuation )
    {
       vector<Sarray> Aex(mNumberOfGrids), A(mNumberOfGrids);
       for( int g=0 ; g < mNumberOfGrids ; g++ )
       {
          Aex[g].copy(AlphaVEp[g][0]); // only checking mechanism m=0
          A[g].copy(AlphaVE[g][0]);
       }
       normOfDifference( Aex, A, errInf, errL2, solInf, a_Sources );
       if ( proc_zero() )
          printf(" Checkpoint solution errors, attenuation at t: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);
    }
    
// Now check Um and AlpphaVEm
    exactSol( t-mDt, Up, AlphaVEp, a_Sources );

    normOfDifference( Up, Um, errInf, errL2, solInf, a_Sources );

    if ( proc_zero() )
       printf("\n Checkpoint errors in Um: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);

    if( m_use_attenuation )
    {
       vector<Sarray> Aex(mNumberOfGrids), A(mNumberOfGrids);
       for( int g=0 ; g < mNumberOfGrids ; g++ )
       {
          Aex[g].copy(AlphaVEp[g][0]); // only checking mechanism m=0
          A[g].copy(AlphaVEm[g][0]);
       }
       normOfDifference( Aex, A, errInf, errL2, solInf, a_Sources );
       if ( proc_zero() )
          printf(" Checkpoint solution errors, attenuation at t-dt: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);
    }
  } // end if twilight testing
  

// test if the spatial operator is self-adjoint (only works without mesh refinement)
  if (m_energy_test && getVerbosity() >= 1 && getNumberOfGrids() == 1)
  {
    if ( proc_zero() )
    {
      printf("Using the intial data to check if the spatial operator is self-adjoint\n");
    }
     
// compute Uacc = L(U) and Vacc=L(V); V=Um
    evalRHS( U, mMu, mLambda, Lu, AlphaVE ); // save Lu in composite grid 'Lu'
    evalRHS( Um, mMu, mLambda, Uacc, AlphaVE ); // save Lu in composite grid 'Lu'
// should not be necessary to communicate across processor boundaries to make ghost points agree
  
// evaluate (V, Uacc) and (U, Vacc) and compare!

// NOTE: scalalarProd() is not implemented for curvilinear grids
    float_sw4 sp_vLu = scalarProduct( Um,Lu );
    float_sw4 sp_uLv = scalarProduct( U,Uacc );
    
    if ( proc_zero() )
    {
       printf("Scalar products (Um, L(U)) = %e and (U, L(Um)) = %e, diff=%e\n", sp_vLu, sp_uLv, sp_vLu-sp_uLv);
    }

// tmp
// save U
  // int g=0;
  // char fname[100];
  // sprintf(fname,"ux-%i.dat",m_myRank);
  // FILE *fp=fopen(fname,"w");
  // printf("Saving tmp file=%s, g=%i, m_jStart=%i, m_jEnd=%i\n", fname, g, m_jStart[g], m_jEnd[g]);
  // for ( int j = m_jStart[g]; j<=m_jEnd[g]; j++ )
  //    fprintf(fp,"%d %e\n", j, U[g](1,35,j,35));
  // fclose(fp);  


// tmp: save L(Um) in U
    // evalRHS( Um, mMu, mLambda, U, AlphaVE );
    
  } // end m_energy_test ...

  if( m_moment_test )
     test_sources( point_sources, a_Sources, F, identsources );

// save initial data on receiver records
  vector<float_sw4> uRec;
  for (int ts=0; ts<a_TimeSeries.size(); ts++)
  {
// can't compute a 2nd order accurate time derivative at this point
// therefore, don't record anything related to velocities for the initial data
    if (a_TimeSeries[ts]->getMode() != TimeSeries::Velocity && a_TimeSeries[ts]->myPoint())
    {
      int i0 = a_TimeSeries[ts]->m_i0;
      int j0 = a_TimeSeries[ts]->m_j0;
      int k0 = a_TimeSeries[ts]->m_k0;
      int grid0 = a_TimeSeries[ts]->m_grid0;
      extractRecordData(a_TimeSeries[ts]->getMode(), i0, j0, k0, grid0, 
			uRec, Um, U); 
      a_TimeSeries[ts]->recordData(uRec);
    }
  }

// save any images for cycle = 0 (initial data), or beginCycle-1 (checkpoint restart)
  update_images( beginCycle-1, t, U, Um, Up, mRho, mMu, mLambda, a_Sources, 1 );
  for( int i3 = 0 ; i3 < mImage3DFiles.size() ; i3++ )
    mImage3DFiles[i3]->update_image( beginCycle-1, t, mDt, U, mRho, mMu, mLambda, mRho, mMu, mLambda, mQp, mQs, mPath, mZ );

  FILE *lf=NULL;
// open file for saving norm of error
  if ( (m_lamb_test || m_point_source_test || m_rayleigh_wave_test || m_error_log) && proc_zero() )
  {
    string path=getPath();

    stringstream fileName;
    if( path != "." )
      fileName << path;
    
    if (m_error_log)
      fileName << m_error_log_file;
    else if (m_lamb_test)
      fileName << "LambErr.txt";
    else if (m_point_source_test)
      fileName << "PointSourceErr.txt";
    else
      fileName << "RayleighErr.txt";
    lf = fopen(fileName.str().c_str(),"w");
  }
  // DEBUG
  //     for( int s = 0 ; s < point_sources.size() ; s++ )
  //        point_sources[s]->print_info();
    
// output flags and settings that affect the run
   if( proc_zero() && mVerbose >= 1 )
   {
      printf("\nReporting SW4 internal flags and settings:\n");
      printf("m_testing=%s, twilight=%s, point_source=%s, moment_test=%s, energy_test=%s, " 
             "lamb_test=%s, rayleigh_test=%s\n",
             m_testing?"yes":"no",
             m_twilight_forcing?"yes":"no",
             m_point_source_test?"yes":"no",
             m_moment_test?"yes":"no",
             m_energy_test?"yes":"no",
             m_lamb_test?"yes":"no",
             m_rayleigh_wave_test?"yes":"no");
      printf("m_use_supergrid=%s\n", usingSupergrid()?"yes":"no");
      printf("End report of internal flags and settings\n\n");
   }
   
  if ( !mQuiet && proc_zero() )
    cout << "  Begin time stepping..." << endl;

  // Prefetch Sarrays before starting time stepping
  SarrayVectorPrefetch(Up);
  SarrayVectorPrefetch(Um);
  SarrayVectorPrefetch(U);
  SarrayVectorPrefetch(Uacc);
  SarrayVectorPrefetch(F);
  SarrayVectorPrefetch(mMu);
  SarrayVectorPrefetch(mLambda);
  SarrayVectorPrefetch(Lu);
  SarrayVectorPrefetch(AlphaVE,m_number_mechanisms);
  SarrayVectorPrefetch(AlphaVEm,m_number_mechanisms);
  SarrayVectorPrefetch(AlphaVEp,m_number_mechanisms);
  // End prefetch


// Begin time stepping loop
  for( int g=0 ; g < mNumberOfGrids ; g++ )
     Up[g].set_to_zero();


// test: compute forcing for the first time step before the loop to get started
       Force( t, F, point_sources, identsources );
// end test
       std::chrono::high_resolution_clock::time_point t1,t2;
// BEGIN TIME STEPPING LOOP
  PROFILER_START;
  SW4_MARK_BEGIN("TIME_STEPPING");
  for( int currentTimeStep = beginCycle; currentTimeStep <= mNumberOfTimeSteps; currentTimeStep++)
  {    
    time_measure[0] = MPI_Wtime();
    if (currentTimeStep==mNumberOfTimeSteps) t1 = std::chrono::high_resolution_clock::now();
// all types of forcing...
    bool trace =false;
    int dbgproc = 1;

    // if( trace && m_myRank == dbgproc )
    //    cout <<" before Forcing" << endl;
    // Force( t, F, point_sources, identsources );

    // if( m_output_detailed_timing )
    //    time_measure[1] = MPI_Wtime();

    // if( trace && m_myRank == dbgproc )
    //    cout <<" after Forcing" << endl;

    if( m_checkfornan )
    {
       check_for_nan( F, 1, "F" );
       check_for_nan( U, 1, "U" );
    }

// evaluate right hand side
    if( m_anisotropic )
       evalRHSanisotropic( U, mC, Lu );
    else
       evalRHS( U, mMu, mLambda, Lu, AlphaVE ); // save Lu in composite grid 'Lu'

    if( m_output_detailed_timing )
       time_measure[1] = MPI_Wtime();

    if( trace && m_myRank == dbgproc )
       cout <<" after evalRHS" << endl;

    if( m_checkfornan )
       check_for_nan( Lu, 1, "Lu pred. " );

// take predictor step, store in Up
    evalPredictor( Up, U, Um, mRho, Lu, F );    
    SW4_MARK_BEGIN("COMM_WINDOW");
    if( m_output_detailed_timing )
       time_measure[2] = MPI_Wtime();

    if( trace &&  m_myRank == dbgproc )
       cout <<" after evalPredictor" << endl;

    SW4_MARK_BEGIN("COMM_ACTUAL");
// communicate across processor boundaries
    for(int g=0 ; g < mNumberOfGrids ; g++ )
       communicate_array( Up[g], g );
    SW4_MARK_END("COMM_ACTUAL");
    if( m_output_detailed_timing )
       time_measure[3] = MPI_Wtime();

    if( trace && m_myRank == dbgproc )
       cout <<" after communicate_array " << endl;

// calculate boundary forcing at time t+mDt
    cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );

    if( m_output_detailed_timing )
       time_measure[4] = MPI_Wtime();

// NEW (Apr. 3, 2017) PC-time stepping for the memory variable
    if( m_use_attenuation && m_number_mechanisms > 0 )
       updateMemVarPred( AlphaVEp, AlphaVEm, U, t );

    if( m_output_detailed_timing )
       time_measure[5] = MPI_Wtime();
    SW4_MARK_END("COMM_WINDOW");
// update ghost points in Up
    if( m_anisotropic )
       enforceBCanisotropic( Up, mC, t+mDt, BCForcing );
    else
       enforceBC( Up, mMu, mLambda, t+mDt, BCForcing );

// NEW
// Impose un-coupled free surface boundary condition with visco-elastic terms
    if( m_use_attenuation && m_number_mechanisms > 0 )
       enforceBCfreeAtt2( Up, mMu, mLambda, AlphaVEp, BCForcing );
    
    if( m_output_detailed_timing )
       time_measure[6] = MPI_Wtime();

    if( trace && m_myRank == dbgproc )
       cout <<" after enforceBC" << endl;

    if( m_checkfornan )
       check_for_nan( Up, 1, "U pred. " );

// Grid refinement interface conditions:
// *** 2nd order in TIME
    if (mOrder == 2)
    {
      SW4_MARK_BEGIN("mOrder=2");
// add super-grid damping terms before enforcing interface conditions
// (otherwise, Up doesn't have the correct values on the interface)
       if (usingSupergrid())
       {
	  addSuperGridDamping( Up, U, Um, mRho );
       }
// Also add Arben's simplified attenuation
       if (m_use_attenuation && m_number_mechanisms == 0)
       {
	 simpleAttenuation( Up );
       }
       if( m_output_detailed_timing )
          time_measure[7] = MPI_Wtime();

// test: compute forcing for next time step here so it can be used in enforceIC()
       Force( t+mDt, F, point_sources, identsources );

       if( m_output_detailed_timing )
          time_measure[8] = MPI_Wtime();
// end test

// interface conditions for 2nd order in time
// NOTE: this routine calls preliminary_predictor for t+dt, which needs F(t+dt). It is computed at the top of next time step
       enforceIC2( Up, U, Um, AlphaVEp, t, F, point_sources );

       if( m_output_detailed_timing )
          time_measure[17] = MPI_Wtime();
       SW4_MARK_END("mOrder=2");
    }
    else // 4th order time stepping
    {
       SW4_MARK_BEGIN("MPI_WTIME");
       if( m_output_detailed_timing )
          time_measure[7] = MPI_Wtime();
       SW4_MARK_END("MPI_WTIME");
// test: precompute F_tt(t)
       Force_tt( t, F, point_sources, identsources );

       if( m_output_detailed_timing )
          time_measure[8] = MPI_Wtime();
// end test 

// *** 4th order in TIME interface conditions for the predictor
// June 14, 2017: adding AlphaVE & AlphaVEm
// NOTE: true means call preliminary_corrector, which needs F_tt(t) & is computed 5 lines down
       enforceIC( Up, U, Um, AlphaVEp, AlphaVE, AlphaVEm, t, true, F, point_sources );

       if( m_output_detailed_timing )
          time_measure[9] = MPI_Wtime();
    }
//
// corrector step for
// *** 4th order in time ***
//
    if (mOrder == 4)
    {
       // Force_tt( t, F, point_sources, identsources );

       // if( m_output_detailed_timing )
       //    time_measure[10] = MPI_Wtime();

       evalDpDmInTime( Up, U, Um, Uacc ); // store result in Uacc
       if( trace && m_myRank == dbgproc )
          cout <<" after evalDpDmInTime" << endl;

       if( m_checkfornan )
	  check_for_nan( Uacc, 1, "uacc " );

// July 22,  4th order update for memory variables
       if( m_use_attenuation && m_number_mechanisms > 0 )
         updateMemVarCorr( AlphaVEp, AlphaVEm, Up, U, Um, t );

       if( m_use_attenuation && m_number_mechanisms > 0 )
          evalDpDmInTimeAtt( AlphaVEp, AlphaVE, AlphaVEm ); // store AlphaVEacc in AlphaVEm
       if( trace && m_myRank == dbgproc )
	  cout <<" after evalDpDmInTimeAtt" << endl;

       if( m_output_detailed_timing )
          time_measure[10] = MPI_Wtime();

       if( m_anisotropic )
	  evalRHSanisotropic( Uacc, mC, Lu );
       else
	  evalRHS( Uacc, mMu, mLambda, Lu, AlphaVEm );

       if( m_output_detailed_timing )
          time_measure[11] = MPI_Wtime();

       if( trace && m_myRank == dbgproc )
	  cout <<" after evalRHS" << endl;
       
       if( m_checkfornan )
	  check_for_nan( Lu, 1, "L(uacc) " );

       evalCorrector( Up, mRho, Lu, F );

       if( m_output_detailed_timing )
          time_measure[12] = MPI_Wtime();

// add in super-grid damping terms
       if (usingSupergrid())
       {
	  addSuperGridDamping( Up, U, Um, mRho );
       }

// Arben's simplified attenuation
       if (m_use_attenuation && m_number_mechanisms == 0)
       {
	 simpleAttenuation( Up );
       }

       if( m_output_detailed_timing )
          time_measure[13] = MPI_Wtime();

// communicate across processor boundaries
       for(int g=0 ; g < mNumberOfGrids ; g++ )
	  communicate_array( Up[g], g );

       if( m_output_detailed_timing )
          time_measure[14] = MPI_Wtime();

// calculate boundary forcing at time t+mDt (do we really need to call this fcn again???)
       cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );

// update ghost points in Up

       if( m_anisotropic )
	  enforceBCanisotropic( Up, mC, t+mDt, BCForcing );
       else
	  enforceBC( Up, mMu, mLambda, t+mDt, BCForcing );

// NEW (Apr. 4, 2017)
// Impose un-coupled free surface boundary condition with visco-elastic terms for 'Up'
       if( m_use_attenuation && (m_number_mechanisms > 0) )
       {
          enforceBCfreeAtt2( Up, mMu, mLambda, AlphaVEp, BCForcing );
       }

       if( m_output_detailed_timing )
          time_measure[15] = MPI_Wtime();
    
       if( trace && m_myRank == dbgproc )
          cout <<" before Forcing" << endl;

// test: compute forcing for next time step here so it can be used in enforceIC()
       Force( t+mDt, F, point_sources, identsources );

       if( m_output_detailed_timing )
          time_measure[16] = MPI_Wtime();

       if( trace && m_myRank == dbgproc )
          cout <<" after Forcing" << endl;
// end test

// interface conditions for the corrector
// June 14, 2017: adding AlphaVE & AlphaVEm
// NOTE: false means call preliminary_predictor for t+dt, which needs F(t+dt). It is computed at the top of next time step
       enforceIC( Up, U, Um, AlphaVEp, AlphaVE, AlphaVEm, t, false, F, point_sources );

       if( m_output_detailed_timing )
          time_measure[17] = MPI_Wtime();
       
    }// end if mOrder == 4
    
    if( m_checkfornan )
       check_for_nan( Up, 1, "Up" );

    //    Um[0].save_to_disk("um-dbg0.bin");
    //    Um[1].save_to_disk("um-dbg1.bin");
    //    U[0].save_to_disk("u-dbg0.bin");
    //    U[1].save_to_disk("u-dbg1.bin");
    //    Up[0].save_to_disk("up-dbg0.bin");
    //    Up[1].save_to_disk("up-dbg1.bin");
    //    mRho[0].save_to_disk("rho-dbg0.bin");
    //    mRho[1].save_to_disk("rho-dbg1.bin");
    //    mMu[0].save_to_disk("mu-dbg0.bin");
    //    mMu[1].save_to_disk("mu-dbg1.bin");
    //    mLambda[0].save_to_disk("lambda-dbg0.bin");
    //    mLambda[1].save_to_disk("lambda-dbg1.bin");
    //    exit(0);
// increment time
    t += mDt;

// periodically, print time stepping info to stdout
    printTime( currentTimeStep, t, currentTimeStep == mNumberOfTimeSteps ); 
    //    printTime( currentTimeStep, t, true ); 

// Images have to be written before the solution arrays are cycled, because both Up and Um are needed
// to compute a centered time derivative
//
// AP: Note to self: Any quantity related to velocities will be lagged by one time step
//
    update_images( currentTimeStep, t, Up, U, Um, mRho, mMu, mLambda, a_Sources, currentTimeStep == mNumberOfTimeSteps );
    for( int i3 = 0 ; i3 < mImage3DFiles.size() ; i3++ )
      mImage3DFiles[i3]->update_image( currentTimeStep, t, mDt, Up, mRho, mMu, mLambda, mRho, mMu, mLambda, 
				       mQp, mQs, mPath, mZ ); // mRho, mMu, mLambda occur twice because we don't use gradRho etc.
    
// save the current solution on receiver records (time-derivative require Up and Um for a 2nd order
// approximation, so do this before cycling the arrays)
    for (int ts=0; ts<a_TimeSeries.size(); ts++)
    {
      if (a_TimeSeries[ts]->myPoint())
      {
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
      }
    }


// Write check point, if requested (timeToWrite returns false if checkpointing is not used)
    if( m_check_point->timeToWrite( t, currentTimeStep, mDt ) )
    {
       double time_chkpt=MPI_Wtime();
       m_check_point->write_checkpoint( t, currentTimeStep, U, Up, AlphaVE, AlphaVEp );
       double time_chkpt_tmp =MPI_Wtime()-time_chkpt;
       if( mVerbose >= 2 )

       {
	  MPI_Allreduce( &time_chkpt_tmp, &time_chkpt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	  if( m_myRank == 0 )
	     cout << "Wallclock time to write check point file " << time_chkpt << " seconds " << endl;
       }
       // Force write all the TimeSeries files for restart
       double time_chkpt_timeseries=MPI_Wtime();
       for (int ts=0; ts<a_TimeSeries.size(); ts++)
       {
         a_TimeSeries[ts]->writeFile();
       }
	     double time_chkpt_timeseries_tmp=MPI_Wtime()-time_chkpt_timeseries;
       if( m_output_detailed_timing )
       {
	        MPI_Allreduce( &time_chkpt_timeseries_tmp, &time_chkpt_timeseries, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	        if( m_myRank == 0 )
	          cout << "Wallclock time to write all checkpoint time series files "
              << time_chkpt_timeseries << " seconds " << endl;
       }
   }

// Energy evaluation, requires all three time levels present, do before cycle arrays.
    if( m_output_detailed_timing )
       time_measure[18] = MPI_Wtime();

// // Energy evaluation, requires all three time levels present, do before cycle arrays.
    if( m_energy_test )
       compute_energy( mDt, currentTimeStep == mNumberOfTimeSteps, Um, U, Up, currentTimeStep  );

// cycle the solution arrays
    cycleSolutionArrays(Um, U, Up, AlphaVEm, AlphaVE, AlphaVEp);

    if( m_output_detailed_timing )
       time_measure[19] = MPI_Wtime();
	  
// evaluate error for some test cases
    if (m_lamb_test || m_point_source_test || m_rayleigh_wave_test )
    {
      float_sw4 errInf=0, errL2=0, solInf=0, solL2=0;
      exactSol( t, Up, AlphaVE, a_Sources ); // store exact solution in Up

      if (m_lamb_test)
	normOfSurfaceDifference( Up, U, errInf, errL2, solInf, solL2, a_Sources);
      else if (m_point_source_test || m_rayleigh_wave_test)
	normOfDifference( Up, U, errInf, errL2, solInf, a_Sources );

      if ( proc_zero() )
// output time, Linf-err, Linf-sol-err
	fprintf(lf, "%e %15.7e %15.7e %15.7e\n", t, errInf, errL2, solInf);
    }

// // See if it is time to write a restart file
// //      if (mRestartDumpInterval > 0 &&  currentTimeStep % mRestartDumpInterval == 0)
// //        serialize(currentTimeStep, U, Um);  

    if( m_output_detailed_timing )
    {
       if (mOrder == 4)
       {
          time_sum[0] += time_measure[19]-time_measure[0]; // total
          time_sum[1] += time_measure[1]-time_measure[0] + time_measure[11]-time_measure[10]; // div-stress
          time_sum[2] += time_measure[8]-time_measure[7] + time_measure[16]-time_measure[15]; // forcing
          time_sum[3] += time_measure[4]-time_measure[3]+time_measure[6]-time_measure[5]+time_measure[7]-time_measure[6] +
             time_measure[15]-time_measure[14];//  bc
          time_sum[4] += time_measure[13]-time_measure[12]; // super-grid
          time_sum[5] += time_measure[3]-time_measure[2] + time_measure[14]-time_measure[13]; // communicate
          time_sum[6] += time_measure[9]-time_measure[8] + time_measure[17]-time_measure[16]; // mesh ref
          time_sum[7] += time_measure[18]-time_measure[17]; // images + time-series
//          time_sum[8] += 0;
          time_sum[8] += time_measure[2]-time_measure[1] + time_measure[5]-time_measure[4] + 
             time_measure[10]-time_measure[9] + time_measure[12]-time_measure[11] +
             time_measure[19]-time_measure[18]; // updates
          
       }
       else
       { // 2nd order in time algorithm
          time_sum[0] += time_measure[19]-time_measure[0]; // total
          time_sum[1] = 0; // update later
          time_sum[2] = 0;
          time_sum[3] = 0;
          time_sum[4] = 0;
          time_sum[5] = 0;
          time_sum[6] = 0;
          time_sum[7] = 0;
          time_sum[8] = 0;
       }
    }
    if (currentTimeStep==mNumberOfTimeSteps) {
      t2 = std::chrono::high_resolution_clock::now();
      if (proc_zero()){
	std::cout<<" Time for the last time step is "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<" ms \n";
      }
    }
  } // end time stepping loop
SW4_MARK_END("TIME_STEPPING");
       //cudaProfilerStop();
  if ( !mQuiet && proc_zero() )
    cout << "  Time stepping finished..." << endl;

//   delete[] wk;

//   if( ind != 0 )
//      delete[] ind;
  
   double time_end_solve = MPI_Wtime();
   print_execution_time( time_start_solve, time_end_solve, "solver phase" );

   if( m_output_detailed_timing )
     print_execution_times( time_sum );

// check the accuracy of the final solution, store exact solution in Up, ignore AlphaVE
   if( exactSol( t, Up, AlphaVEp, a_Sources ) )
   {
     float_sw4 errInf=0, errL2=0, solInf=0, solL2=0;

// tmp: output exact sol for Lamb's prolem 
//      cout << *mGlobalUniqueSources[0] << endl;
//       Image* im = new Image( this, 0, 1, 0, 1, "exact", 1 , Image::UZ, Image::Z, 0.0, true );
//       im->computeGridPtIndex();
//       im->allocatePlane();
//       im->computeImageQuantity(Up, 3); // z-component
//       string path=".";
//       im->writeImagePlane_2(1,path);

// depending on the test case, we should compare in the interior, or only on the surface
      if (m_lamb_test)
	normOfSurfaceDifference( Up, U, errInf, errL2, solInf, solL2, a_Sources);
      else
	normOfDifference( Up, U, errInf, errL2, solInf, a_Sources );

      if ( proc_zero() )
      {         
	 printf("\n Final solution errors: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);

// output time, Linf-err, Linf-sol-err
         if ( m_error_log )
         {
            fprintf(lf, "Final time\n");
            fprintf(lf, "%e\n", t);
            fprintf(lf, "Displacement variables (errInf, errL2, solInf)\n");            
            fprintf(lf, "%15.7e %15.7e %15.7e\n", errInf, errL2, solInf);
         }
         
      }
      
      if( m_twilight_forcing && m_use_attenuation )
      {
         vector<Sarray> Aex(mNumberOfGrids), A(mNumberOfGrids);
         for( int g=0 ; g < mNumberOfGrids ; g++ )
	 {
	    Aex[g].copy(AlphaVEp[g][0]);
            A[g].copy(AlphaVE[g][0]);
	 }
	 normOfDifference( Aex, A, errInf, errL2, solInf, a_Sources );
	 if ( proc_zero() )
         {
	    printf("\n Final solution errors, attenuation: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);
            if ( m_error_log )
            {
               fprintf(lf, "Attennuation variables (errInf, errL2, solInf)\n");
               fprintf(lf, "%15.7e %15.7e %15.7e\n", errInf, errL2, solInf);
            }
            
         }
         
      }
// test
//      int g=mNumberOfCartesianGrids - 1;
//      Up[g].set_to_minusOne();
//      U[g].set_to_zero();
//      normOfSurfaceDifference( Up, U, errInf, errL2, solInf, solL2, a_Sources);
//      if ( proc_zero() )
//	 printf("\n Surface norm of 1: Inf = %15.7e, L2 = %15.7e\n", errInf, errL2);
   } // end if exactSol
   
// close error log file for testing
  if ((m_lamb_test || m_point_source_test || m_rayleigh_wave_test || m_error_log) && proc_zero() )
  {
    fclose(lf);
    printf("**** Closing file with solution errors for testing\n");
  }

   finalizeIO();
   cout.flush(); cerr.flush();

   // Give back memory
   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      for(int side=0; side < 6; side++)
	 if( BCForcing[g][side] != NULL )
	   ::operator delete[] (BCForcing[g][side],Managed);
      delete[] BCForcing[g];
   }
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   //std::cerr<<"Deleting point_sources "<<myRank<"\n";

   //std::cerr<<"WARNING :: delete of point_sources turned off for speed in lines 1112-1113 of solve.C \n";
   for( int s = 0 ; s < point_sources.size(); s++ )
     delete point_sources[s];
   //std::cerr<<"Done "<<myRank<<"\n";
// why is this barrier needed???
   MPI_Barrier(MPI_COMM_WORLD);

//   if( m_forcing->knows_exact() )
//      computeSolutionError(U, mTime, AlphaVE ); // note that final solution ends up in U after the call to cycleSolutionArrays()

} // end EW::solve()

//------------------------------------------------------------------------
void EW::cycleSolutionArrays(vector<Sarray> & a_Um, vector<Sarray> & a_U, vector<Sarray> & a_Up, 
			     vector<Sarray*> & a_AlphaVEm, vector<Sarray*> & a_AlphaVE, vector<Sarray*> & a_AlphaVEp)
 {
   SW4_MARK_FUNCTION;
  for (int g=0; g<mNumberOfGrids; g++)
  {
    float_sw4 *tmp = a_Um[g].c_ptr();
    a_Um[g].reference(a_U[g].c_ptr());
    a_U[g].reference(a_Up[g].c_ptr());
    a_Up[g].reference(tmp);
    for( int a = 0 ; a < m_number_mechanisms ; a++ )
    {
       float_sw4 *tmp = a_AlphaVEm[g][a].c_ptr();
       a_AlphaVEm[g][a].reference(a_AlphaVE[g][a].c_ptr());
       a_AlphaVE[g][a].reference(a_AlphaVEp[g][a].c_ptr());
       a_AlphaVEp[g][a].reference(tmp);
    }
  }
}

//---------------------------------------------------------------------------
void EW::enforceBC( vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		    float_sw4 t, vector<float_sw4 **> & a_BCForcing )
{
  SW4_MARK_FUNCTION;
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz;
  float_sw4 *u_ptr, *mu_ptr, *la_ptr, h;
  boundaryConditionType *bcType_ptr;
  float_sw4 *bforce_side0_ptr, *bforce_side1_ptr, *bforce_side2_ptr, *bforce_side3_ptr, *bforce_side4_ptr, *bforce_side5_ptr;
  int *wind_ptr;
  float_sw4 om=0, ph=0, cv=0;
    
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    u_ptr    = a_U[g].c_ptr();
    mu_ptr    = a_Mu[g].c_ptr();
    la_ptr    = a_Lambda[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    nx = m_global_nx[g];
    ny = m_global_ny[g];
    nz = m_global_nz[g];
    
    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    bcType_ptr = m_bcType[g]; // get a pointer to the boundary conditions for grid 'g'
    
    wind_ptr = m_BndryWindow[g];// get a pointer to the boundary window array for grid 'g'
    //    cout << "Grid: " << g << endl;
    //    for( int s=0 ; s < 6 ; s++ )
    //       cout << " side " << s << " wind = " << wind_ptr[6*s] << " " << wind_ptr[6*s+1] << " " << wind_ptr[6*s+2] << " " 
    //	    << wind_ptr[6*s+3] << " " << wind_ptr[6*s+4] << " " << wind_ptr[6*s+5] << endl;
    int topo=topographyExists() && g == mNumberOfGrids-1;
    
// THESE ARRAYS MUST BE FILLED IN BEFORE CALLING THIS ROUTINE
// for periodic bc, a_BCForcing[g][s] == NULL, so you better not access theses arrays in that case
    bforce_side0_ptr = a_BCForcing[g][0]; // low-i bndry forcing array pointer
    bforce_side1_ptr = a_BCForcing[g][1]; // high-i bndry forcing array pointer
    bforce_side2_ptr = a_BCForcing[g][2]; // low-j bndry forcing array pointer
    bforce_side3_ptr = a_BCForcing[g][3]; // high-j bndry forcing array pointer
    bforce_side4_ptr = a_BCForcing[g][4]; // low-k bndry forcing array pointer
    bforce_side5_ptr = a_BCForcing[g][5]; // high-k bndry forcing array pointer
    
    if( usingSupergrid() )
    {
       if( m_croutines )
	  bcfortsg_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, 
		    wind_ptr, nx, ny, nz,
		    u_ptr, h, bcType_ptr, m_sbop, mu_ptr, la_ptr, t,
		    bforce_side0_ptr, bforce_side1_ptr, 
		    bforce_side2_ptr, bforce_side3_ptr, 
		    bforce_side4_ptr, bforce_side5_ptr, 
		    om, ph, cv, m_sg_str_x[g], m_sg_str_y[g] );
       else
	  bcfortsg( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
		    wind_ptr, &nx, &ny, &nz,
		    u_ptr, &h, bcType_ptr, m_sbop, mu_ptr, la_ptr, &t,
		    bforce_side0_ptr, bforce_side1_ptr, 
		    bforce_side2_ptr, bforce_side3_ptr, 
		    bforce_side4_ptr, bforce_side5_ptr, 
		    &om, &ph, &cv, m_sg_str_x[g], m_sg_str_y[g] );
       int side;
       if( topo == 1 && m_bcType[g][4] == bStressFree )
       {
	  side = 5;
	  if( m_croutines )
	     freesurfcurvisg_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
				 nz, side, u_ptr, mu_ptr, la_ptr, mMetric.c_ptr(),
				 m_sbop, bforce_side4_ptr, m_sg_str_x[g], m_sg_str_y[g] );
	  else
	     freesurfcurvisg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			     &nz, &side, u_ptr, mu_ptr, la_ptr, mMetric.c_ptr(),
			     m_sbop, bforce_side4_ptr, m_sg_str_x[g], m_sg_str_y[g] );
       }
    }
    else
    {
       if( m_croutines )
	  bcfort_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, 
		     wind_ptr, nx, ny, nz,
		     u_ptr, h, bcType_ptr, m_sbop, mu_ptr, la_ptr, t,
		     bforce_side0_ptr, bforce_side1_ptr, 
		     bforce_side2_ptr, bforce_side3_ptr, 
		     bforce_side4_ptr, bforce_side5_ptr, 
		     om, ph, cv, topo );
       else
	  bcfort( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
		  wind_ptr, &nx, &ny, &nz,
		  u_ptr, &h, bcType_ptr, m_sbop, mu_ptr, la_ptr, &t,
		  bforce_side0_ptr, bforce_side1_ptr, 
		  bforce_side2_ptr, bforce_side3_ptr, 
		  bforce_side4_ptr, bforce_side5_ptr, 
		  &om, &ph, &cv, &topo );
       int side;
       if( topo == 1 && m_bcType[g][4] == bStressFree )
       {
	  side = 5;
	  if( m_croutines )
	     freesurfcurvi_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
			       side, u_ptr, mu_ptr, la_ptr, mMetric.c_ptr(), m_sbop, bforce_side4_ptr );
	  else
	     freesurfcurvi(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
			   &side, u_ptr, mu_ptr, la_ptr, mMetric.c_ptr(), m_sbop, bforce_side4_ptr );
       }
       if( topo == 1 && m_bcType[g][5] == bStressFree )
       {
	  side = 6;
	  if( m_croutines )
	     freesurfcurvi_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
			      side, u_ptr, mu_ptr, la_ptr, mMetric.c_ptr(), m_sbop, bforce_side5_ptr );
	  else
	     freesurfcurvi(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
			   &side, u_ptr, mu_ptr, la_ptr, mMetric.c_ptr(), m_sbop, bforce_side5_ptr );
       }
    }
  }
// interface between curvilinear and top Cartesian grid
   if (topographyExists())
   {
      int nc = 3;
      int g = mNumberOfCartesianGrids-1;
      int gc = mNumberOfGrids-1;
      int mgp = getNumberOfGhostPoints();
      //      float_sw4 nrm[3]={0,0,0};
      //      int q, i, j;
// inject solution values between lower boundary of gc and upper boundary of g

      SView &a_UgV = a_U[g].getview();
      SView &a_UgcV = a_U[gc].getview();
      ASSERT_MANAGED(a_U[g].c_ptr());
      ASSERT_MANAGED(a_U[gc].c_ptr());
      int kstartg = m_kStart[g];
      int kendgc = m_kEnd[gc];
      RAJA::RangeSegment j_range(m_jStart[g],m_jEnd[g]+1);
      RAJA::RangeSegment i_range(m_iStart[g],m_iEnd[g]+1);
      RAJA::kernel<ENFORCEBC_CORR_EXEC_POL1>(
			   RAJA::make_tuple(j_range,i_range),
			   [=]RAJA_DEVICE (int j,int i) {
// #pragma omp parallel for
//       for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
// 	 for( int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
// 	 {
// assign ghost points in the Cartesian grid
	    for (int q = 0; q < mgp; q++) // only once when m_ghost_points==1
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_UgV(c,i,j,kstartg + q) = a_UgcV(c,i,j,kendgc-2*mgp + q);
	    }
// assign ghost points in the Curvilinear grid
	    for (int q = 0; q <= mgp; q++) // twice when m_ghost_points==1 (overwrites solution on the common grid line)
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_UgcV(c,i,j,kendgc-q) = a_UgV(c,i,j,kstartg+2*mgp - q);
	    }
	    //	    // Verify that the grid lines are the same
	    //            if( i>= 1 && i <= m_global_nx[g] && j >= 1 && j <= m_global_ny[g] )
	    //            for( int c=1; c <= nc ; c++ )
	    //	    {
	    //	       float_sw4 ndiff=fabs(a_U[gc](c,i,j,m_kEnd[gc]-m_ghost_points)-a_U[g](c,i,j,m_kStart[g]+m_ghost_points));
	    //	       if( ndiff > nrm[c-1] )
	    //		  nrm[c-1] = ndiff;
	    //	    }
			   }); // End of RAJA loop
      //      cout << "Difference of curvilinear and cartesian common grid line = "<< nrm[0] << " " << nrm[1] << " " << nrm[2] << endl;
   }
}

//---------------------------------------------------------------------------
void EW::enforceBCanisotropic( vector<Sarray> & a_U, vector<Sarray>& a_C,
		    float_sw4 t, vector<float_sw4 **> & a_BCForcing )
{
  SW4_MARK_FUNCTION;
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz;
  float_sw4 *u_ptr, *c_ptr, h;
  boundaryConditionType *bcType_ptr;
  float_sw4 *bforce_side0_ptr, *bforce_side1_ptr, *bforce_side2_ptr, *bforce_side3_ptr, *bforce_side4_ptr, *bforce_side5_ptr;
  int *wind_ptr;
  //float_sw4 om=0, ph=0, cv=0;
    
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    u_ptr    = a_U[g].c_ptr();
    c_ptr    = a_C[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    nx = m_global_nx[g];
    ny = m_global_ny[g];
    nz = m_global_nz[g];
    
    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    bcType_ptr = m_bcType[g]; // get a pointer to the boundary conditions for grid 'g'
    
    wind_ptr = m_BndryWindow[g];// get a pointer to the boundary window array for grid 'g'
    //    cout << "Grid: " << g << endl;
    //    for( int s=0 ; s < 6 ; s++ )
    //       cout << " side " << s << " wind = " << wind_ptr[6*s] << " " << wind_ptr[6*s+1] << " " << wind_ptr[6*s+2] << " " 
    //	    << wind_ptr[6*s+3] << " " << wind_ptr[6*s+4] << " " << wind_ptr[6*s+5] << endl;
    int topo=topographyExists() && g == mNumberOfGrids-1;
    
// THESE ARRAYS MUST BE FILLED IN BEFORE CALLING THIS ROUTINE
// for periodic bc, a_BCForcing[g][s] == NULL, so you better not access
// theses arrays in that case
    bforce_side0_ptr = a_BCForcing[g][0]; // low-i bndry forcing array pointer
    bforce_side1_ptr = a_BCForcing[g][1]; // high-i bndry forcing array pointer
    bforce_side2_ptr = a_BCForcing[g][2]; // low-j bndry forcing array pointer
    bforce_side3_ptr = a_BCForcing[g][3]; // high-j bndry forcing array pointer
    bforce_side4_ptr = a_BCForcing[g][4]; // low-k bndry forcing array pointer
    bforce_side5_ptr = a_BCForcing[g][5]; // high-k bndry forcing array pointer
    
    if( m_croutines )
       bcfortanisg_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, 
		    wind_ptr, nx, ny, nz,
		    u_ptr, h, bcType_ptr, m_sbop, c_ptr, 
		    bforce_side0_ptr, bforce_side1_ptr, 
		    bforce_side2_ptr, bforce_side3_ptr, 
		    bforce_side4_ptr, bforce_side5_ptr, 
		    m_sg_str_x[g], m_sg_str_y[g] );
    else
       bcfortanisg( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
		    wind_ptr, &nx, &ny, &nz,
		    u_ptr, &h, bcType_ptr, m_sbop, c_ptr, 
		    bforce_side0_ptr, bforce_side1_ptr, 
		    bforce_side2_ptr, bforce_side3_ptr, 
		    bforce_side4_ptr, bforce_side5_ptr, 
		    m_sg_str_x[g], m_sg_str_y[g] );
       if( topographyExists() && g == mNumberOfGrids-1 && m_bcType[g][4] == bStressFree )
       {
	  int fside = 5;
	  float_sw4* cc_ptr = mCcurv.c_ptr();
	  if( m_croutines )
	     bcfreesurfcurvani_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
			       nz, u_ptr, cc_ptr, fside, m_sbop, bforce_side4_ptr,
			       bforce_side5_ptr, m_sg_str_x[g], m_sg_str_y[g] );
	  else
	     bcfreesurfcurvani(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			       &nz, u_ptr, cc_ptr, &fside, m_sbop, bforce_side4_ptr,
			       bforce_side5_ptr, m_sg_str_x[g], m_sg_str_y[g] );
       }
  }
  update_curvilinear_cartesian_interface( a_U );
}

//-----------------------------------------------------------------------
void EW::update_curvilinear_cartesian_interface( vector<Sarray>& a_U )
{
   if (topographyExists())
   {
      int g  = mNumberOfCartesianGrids-1;
      int gc = mNumberOfGrids-1;
      int nc = a_U[0].m_nc;
// inject solution values between lower boundary of gc and upper boundary of g
#pragma omp parallel for
      for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	 for( int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
	 {
// assign ghost points in the Cartesian grid
	    for (int q = 0; q < m_ghost_points; q++) 
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_U[g](c,i,j,m_kStart[g] + q) = a_U[gc](c,i,j,m_kEnd[gc]-2*m_ghost_points + q);
	    }
// assign ghost points in the Curvilinear grid
	    for (int q = 0; q <= m_ghost_points; q++) // (overwrites solution on the common grid line)
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_U[gc](c,i,j,m_kEnd[gc]-q) = a_U[g](c,i,j,m_kStart[g]+2*m_ghost_points - q);
	    }
	 }
   }
}
//-----------------Mesh refinement interface condition for 4th order predictor-corrector scheme------------------------
void EW::enforceIC( vector<Sarray>& a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
                    vector<Sarray*>& a_AlphaVEp, vector<Sarray*>& a_AlphaVE, vector<Sarray*>& a_AlphaVEm, 
		    float_sw4 time, bool predictor, vector<Sarray> &F, vector<GridPointSource*> & point_sources )
{
SW4_MARK_FUNCTION;
   for( int g = 0 ; g < mNumberOfCartesianGrids-1 ; g++ )
   {
      // Interpolate between g and g+1, assume factor 2 refinement with at least three ghost points
      VERIFY2( m_ghost_points >= 3,
		  "enforceIC Error: "<<
                  "Number of ghost points must be three or more, not " << m_ghost_points );

      //Sarray Unextf, Unextc, Bf, Bc, Uf_tt, Uc_tt;
      
      
      int ibf=m_iStart[g+1], ief=m_iEnd[g+1], jbf=m_jStart[g+1], jef=m_jEnd[g+1];
      int kf = m_global_nz[g+1];
      int ibc=m_iStart[g], iec=m_iEnd[g], jbc=m_jStart[g], jec=m_jEnd[g];
      int kc = 1;
      SW4_MARK_BEGIN("enforceIC::Allocs");
  // fine side
      Sarray Unextf(3,ibf,ief,jbf,jef,kf,kf,__FILE__,__LINE__); // only needs k=kf (on the interface)
      Sarray Bf(3,ibf,ief,jbf,jef,kf,kf,__FILE__,__LINE__);
// coarse side
      Sarray Unextc(3,ibc,iec,jbc,jec,kc,kc,__FILE__,__LINE__); // only needs k=kc (on the interface)
      Sarray Bc(3,ibc,iec,jbc,jec,kc,kc,__FILE__,__LINE__);
      Unextf.set_to_zero_async();
      Bf.set_to_zero_async();
      //std::cout<<"BF ARRAY "<<Bf.c_ptr()<<"\n";
      Unextc.set_to_zero_async();
      Bc.set_to_zero_async();
// to compute the corrector we need the acceleration in the vicinity of the interface
      Sarray Uf_tt(3,ibf,ief,jbf,jef,kf-7,kf+1,__FILE__,__LINE__);
      Sarray Uc_tt(3,ibc,iec,jbc,jec,kc-1,kc+7,__FILE__,__LINE__);
// reuse Utt to hold the acceleration of the memory variables
      SW4_MARK_END("enforceIC::Allocs");
    // Set to zero the ghost point values that are unknowns when solving the interface condition. Assume that Dirichlet data
    // is already set on ghost points on the (supergrid) sides, which are not treated as unknown variables.
      dirichlet_hom_ic( a_Up[g+1], g+1, kf+1, true ); // inside=true
      dirichlet_hom_ic( a_Up[g], g, kc-1, true );
      if( m_doubly_periodic )
      {
	 dirichlet_hom_ic( a_Up[g+1], g+1, kf+1, false ); // inside=false
	 dirichlet_hom_ic( a_Up[g], g, kc-1, false );
      }

      if( predictor ) // In the predictor step, (Unextc, Unextf) represent the displacement after the corrector step
      {
//  REMARK: June 15, 2017: if predictor == true, the memory variable a_alphaVEp holds the predicted
// (2nd order) values on entry. However, the interior contribution to the displacement on the
// interface depends on the corrected memory variable. The memory variable is updated within
// compute_preliminary_corrector. 

// get the interior contribution to the displacements on the interface for the corrector (depends on the corrector value of AlphaVEp)
	 compute_preliminary_corrector( a_Up[g+1], a_U[g+1], a_Um[g+1],
                                        a_AlphaVEp[g+1], a_AlphaVE[g+1], a_AlphaVEm[g+1],
                                        Uf_tt, Unextf, g+1, kf, time, F[g+1], point_sources );
         compute_preliminary_corrector( a_Up[g], a_U[g], a_Um[g],
                                        a_AlphaVEp[g], a_AlphaVE[g], a_AlphaVEm[g],
                                        Uc_tt, Unextc, g, kc, time, F[g], point_sources );
	 if( !m_doubly_periodic )
	 {
// dirichlet conditions for Unextc in super-grid layer at time t+dt
            dirichlet_LRic( Unextc, g, kc, time+mDt, 1 ); 
	 }
      }
      else // In the corrector step, (Unextc, Unextf) represent the displacement after next predictor step
      {
	 compute_preliminary_predictor( a_Up[g+1], a_U[g+1], a_AlphaVEp[g+1], Unextf, g+1, kf,  time+mDt,
                                        F[g+1], point_sources );
	 compute_preliminary_predictor( a_Up[g], a_U[g], a_AlphaVEp[g], Unextc, g, kc, time+mDt, 
                                        F[g], point_sources );

	 if( !m_doubly_periodic )
	 {
// dirichlet conditions for Unextc in super-grid layer at time t+2*dt
            dirichlet_LRic( Unextc, g, kc, time+2*mDt, 1 ); 
	 }
      }

      compute_icstresses( a_Up[g+1], Bf, g+1, kf, m_sg_str_x[g+1], m_sg_str_y[g+1]);
      compute_icstresses( a_Up[g], Bc, g, kc, m_sg_str_x[g], m_sg_str_y[g] );

// NEW June 13, 2017: add in the visco-elastic boundary traction
      if( m_use_attenuation && m_number_mechanisms > 0 )
      {
         for( int a=0 ; a < m_number_mechanisms ; a++ )
         {
// the visco-elastic stresses depend on the predictor values of AlphaVEp (if predictor == true)
            add_ve_stresses( a_AlphaVEp[g+1][a], Bf, g+1, kf, a, m_sg_str_x[g+1], m_sg_str_y[g+1]); 
            add_ve_stresses( a_AlphaVEp[g][a],     Bc, g,    kc, a, m_sg_str_x[g],     m_sg_str_y[g]   );
         }
      }

// from enforceIC2()
      if( !m_doubly_periodic )
      {
//  dirichlet condition for Bf in the super-grid layer at time t+dt (also works with twilight)
         dirichlet_LRstress( Bf, g+1, kf, time+mDt, 1 ); 
      }

      SW4_MARK_BEGIN("enforceIC::MPI2DCOMM");
      communicate_array_2d( Unextf, g+1, kf );
      communicate_array_2d( Unextc, g, kc );
      communicate_array_2d( Bf, g+1, kf );
      communicate_array_2d( Bc, g, kc );
      SW4_MARK_END("enforceIC::MPI2DCOMM");
      // Up to here, interface stresses and displacement (Bc,Bf) and (Unextc, Unextf) were computed with 
      // correct ghost point values in the corners (supergrid layers).
      // In the following iteration, we use centered formulas for interpolation and restriction, all the 
      // way up to the Dirichlet boundaries.
      // We must therefore set the corner ghost point values to zero.
      // This is needed in the call to consintp() because Up[g+1] (fine grid) is used in a 7-point restriction
      // stencil, and Up[g] (coarse grid) is used in a 4-point interpolation stencil (ic-1,jc-1) -> (ic+2, jc+2)
      //
      // Note: the inside flag is now false -> only zero out the ghost points in the corners, i.e., above (or below) 
      // the sides where dirichlet boundary conditions are imposed. 
      if( !m_doubly_periodic )
      {
	 dirichlet_hom_ic( a_Up[g+1], g+1, kf+1, false );
	 dirichlet_hom_ic( a_Up[g], g, kc-1, false );
      }
      // Initial guesses for grid interface iteration
      gridref_initial_guess( a_Up[g+1], g+1, false );
      gridref_initial_guess( a_Up[g], g, true );

      float_sw4 cof = predictor ? 12 : 1;
      int is_periodic[2] ={0,0};
      if( m_doubly_periodic )
	 is_periodic[0] = is_periodic[1] = 1;
      
      consintp( a_Up[g+1], Unextf, Bf, mMu[g+1], mLambda[g+1], mRho[g+1], mGridSize[g+1],
		a_Up[g],   Unextc, Bc, mMu[g],   mLambda[g],   mRho[g],   mGridSize[g], 
		cof, g, g+1, is_periodic );
      //      CHECK_INPUT(false," controlled termination");

      // Finally, restore the ghost point values on the sides of the domain.
      // Note: these ghost point values might never be used ?
      if( !m_doubly_periodic )
      {
         dirichlet_LRic( a_Up[g+1], g+1, kf+1, time+mDt, 1 );
         dirichlet_LRic( a_Up[g], g, kc-1, time+mDt, 1 );
      }
       
   } // end for g...
} // enforceIC

//-----------------------Special case for 2nd order time stepper----------------------------------------------------
void EW::enforceIC2( vector<Sarray>& a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
                     vector<Sarray*>& a_AlphaVEp,
                     float_sw4 time, vector<Sarray> &F, vector<GridPointSource*> &point_sources )
{
SW4_MARK_FUNCTION;
   bool predictor = false;   // or true???
   for( int g = 0 ; g < mNumberOfCartesianGrids-1 ; g++ )
   {
      // Interpolate between g and g+1, assume factor 2 refinement with at least three ghost points
      VERIFY2( m_ghost_points >= 3,
		  "enforceIC2 Error: "<<
                  "Number of ghost points must be three or more, not " << m_ghost_points );

      Sarray Unextf, Unextc, Bf, Bc;
      int ibf=m_iStart[g+1], ief=m_iEnd[g+1], jbf=m_jStart[g+1], jef=m_jEnd[g+1];
      int kf = m_global_nz[g+1];
      int ibc=m_iStart[g], iec=m_iEnd[g], jbc=m_jStart[g], jec=m_jEnd[g];
      int kc = 1;
// fine side
      Unextf.define(3,ibf,ief,jbf,jef,kf,kf); // only needs k=kf (on the interface)
      Bf.define(3,ibf,ief,jbf,jef,kf,kf);
// coarse side
      Unextc.define(3,ibc,iec,jbc,jec,kc,kc); // only needs k=kc (on the interface)
      Bc.define(3,ibc,iec,jbc,jec,kc,kc);

// test: check that the condition Up[g+1](kf) = P Up[g](1) is satisfied on the interface
      if (mVerbose >= 3)
         check_displacement_continuity( a_Up[g+1], a_Up[g], g+1, g);
      
    //  Zero out the ghost point values that are unknowns when solving the interface condition. Assume that Dirichlet data
    // are already set on ghost points on the other (supergrid) sides, which are not treated as unknown variables.
      dirichlet_hom_ic( a_Up[g+1], g+1, kf+1, true ); // inside=true
      dirichlet_hom_ic( a_Up[g], g, kc-1, true );

      if( m_doubly_periodic )
      {
	 dirichlet_hom_ic( a_Up[g+1], g+1, kf+1, false ); // inside=false
	 dirichlet_hom_ic( a_Up[g], g, kc-1, false );
      }

//  compute contribution to the displacements at the next time level (Unextc, Unextf) from the interior grid points in Up
// note: t+dt refers to the time level for the forcing. Unextf lives on time level t+2*dt
//
      // Check super-grid terms !
      //
      compute_preliminary_predictor( a_Up[g+1], a_U[g+1], a_AlphaVEp[g+1], Unextf, g+1, kf, time+mDt, F[g+1], point_sources );
      compute_preliminary_predictor( a_Up[g], a_U[g], a_AlphaVEp[g], Unextc, g, kc, time+mDt, F[g], point_sources );

      if( !m_doubly_periodic )
      {
// dirichlet conditions for Unextc in super-grid layer at time t+2*dt
         dirichlet_LRic( Unextc, g, kc, time+2*mDt, 1 ); 
      }

// test: assign exact (twilight) solution to a_Up
//      initialData(time+mDt, a_Up, a_AlphaVEp);
// end test
//  compute contribution to the normal stresses (Bc, Bf) from the interior grid points in Up
      compute_icstresses( a_Up[g+1], Bf, g+1, kf, m_sg_str_x[g+1], m_sg_str_y[g+1]);
      compute_icstresses( a_Up[g], Bc, g, kc, m_sg_str_x[g], m_sg_str_y[g]);

// add in the visco-elastic boundary traction
      if( m_use_attenuation && m_number_mechanisms > 0 )
      {
         for( int a=0 ; a < m_number_mechanisms ; a++ )
         {
            add_ve_stresses( a_AlphaVEp[g+1][a], Bf, g+1, kf, a, m_sg_str_x[g+1], m_sg_str_y[g+1]);
            add_ve_stresses( a_AlphaVEp[g][a],     Bc, g,    kc, a, m_sg_str_x[g],     m_sg_str_y[g]   );
         }
      }

      if( !m_doubly_periodic )
      {
//  dirichlet condition for Bf in the super-grid layer at time t+dt (also works with Dirichlet condition and twilight)
         dirichlet_LRstress( Bf, g+1, kf, time+mDt, 1 ); 
      }
      
      communicate_array_2d( Unextf, g+1, kf );
      communicate_array_2d( Unextc, g, kc );
      communicate_array_2d( Bf, g+1, kf );
      communicate_array_2d( Bc, g, kc );

      // Up to here, interface stresses and displacement (Bc,Bf) and (Unextc, Unextf) were computed with correct ghost point values
      // in the corners.
      // In the following iteration, we use centered formulas for interpolation and restriction, all the way up to the Dirichlet boundaries.
      // We must therefore set the corner ghost point values to zero.
      // This is needed in the call to consintp() because Up[g+1] (fine grid) is used in a 7-point restriction
      // stencil, and Up[g] (coarse grid) is used in a 4-point interpolation stencil (ic-1,jc-1) -> (ic+2, jc+2)
      //
      // Note: the inside flag is now false -> only zero out the ghost points in the corners, i.e., above (or below) the sides where
      // dirichlet boundary conditions are imposed. 
      if( !m_doubly_periodic ) 
      {
         dirichlet_hom_ic( a_Up[g+1], g+1, kf+1, false );
         dirichlet_hom_ic( a_Up[g], g, kc-1, false );
      }
      // Initial guesses for grid interface iteration
      gridref_initial_guess( a_Up[g+1], g+1, false );
      gridref_initial_guess( a_Up[g], g, true );
      float_sw4 cof = predictor ? 12 : 1;
      int is_periodic[2] ={0,0};
      if( m_doubly_periodic )
	 is_periodic[0] = is_periodic[1] = 1;

// // testing (twilight)
//       initialData(time+mDt, a_Up, a_AlphaVEp);
//       // Iteratively determine the ghost point values in Up to satisfy the jump conditions
//       checkintp( a_Up[g+1], Unextf, Bf, mMu[g+1], mLambda[g+1], mRho[g+1], mGridSize[g+1],
//                  a_Up[g],   Unextc, Bc, mMu[g],   mLambda[g],   mRho[g],   mGridSize[g], 
//                  cof, g, g+1, is_periodic, time+mDt);
// // end test
      
// Iteratively determine the ghost point values in Up to satisfy the jump conditions
      consintp( a_Up[g+1], Unextf, Bf, mMu[g+1], mLambda[g+1], mRho[g+1], mGridSize[g+1],
		a_Up[g],   Unextc, Bc, mMu[g],   mLambda[g],   mRho[g],   mGridSize[g], 
		cof, g, g+1, is_periodic);
      //      CHECK_INPUT(false," controlled termination");

      // Finally, restore the corner ghost point values (above and below) the Dirichlet sides of the domain.
      // Note: these ghost point values might never be used
      if( !m_doubly_periodic )
      {
         dirichlet_LRic( a_Up[g+1], g+1, kf+1, time+mDt, 1 );
         dirichlet_LRic( a_Up[g], g, kc-1, time+mDt, 1 );
      }
   }
}// end enforceIC2


//-----------------------------------------------------------------------
void EW::check_corrector( Sarray& Uf, Sarray& Uc, Sarray& Unextf, Sarray& Unextc, int kf, int kc )
{
   int ic =2, jc=4, c=1;
   int i = 2*ic-1, j=2*jc-1;
   //   cout <<"check " << Uf(c,i,j,kf) << " " << Uc(c,ic,jc,kc) << " " << Uf(c,i,j,kf)-Uc(c,ic,jc,kc) << endl;
   //   cout << "  next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << endl;
   i = 2*ic;
   //   j = 2*jc;
   float_sw4 pci   = (9*(Uc(c,ic,jc,kc)+Uc(c,ic+1,jc,kc))-(Uc(c,ic-1,jc,kc)+Uc(c,ic+2,jc,kc)))/16;
   float_sw4 pcim  = (9*(Uc(c,ic,jc-1,kc)+Uc(c,ic+1,jc-1,kc))-(Uc(c,ic-1,jc-1,kc)+Uc(c,ic+2,jc-1,kc)))/16;
   float_sw4 pcip  = (9*(Uc(c,ic,jc+1,kc)+Uc(c,ic+1,jc+1,kc))-(Uc(c,ic-1,jc+1,kc)+Uc(c,ic+2,jc+1,kc)))/16;
   float_sw4 pcipp = (9*(Uc(c,ic,jc+2,kc)+Uc(c,ic+1,jc+2,kc))-(Uc(c,ic-1,jc+2,kc)+Uc(c,ic+2,jc+2,kc)))/16;

   //float_sw4 pc = ( 9*(pci+pcip)-(pcim+pcipp))/16;

   float_sw4 pcj = (9*(Uc(c,ic,jc,kc)+Uc(c,ic,jc+1,kc))-(Uc(c,ic,jc-1,kc)+Uc(c,ic,jc+2,kc)))/16;
   cout <<"check " << Uf(c,i,j,kf) << " " << Uc(c,ic,jc,kc) << " " << pci << " " << Uf(c,i,j,kf)-pci << endl;
   //   cout << "  next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << endl;
   float_sw4 pcni   = (9*(Unextc(c,ic,jc,kc)+Unextc(c,ic+1,jc,kc))-(Unextc(c,ic-1,jc,kc)+Unextc(c,ic+2,jc,kc)))/16;
     cout << "  check next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << " " << pcni << " " << Unextf(c,i,j,kf)-pcni << endl;
   //   cout << " check " << Uc(c,ic-2,jc,kc) << " " << Unextc(c,ic-2,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic-1,jc,kc) << " " << Unextc(c,ic-1,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic,jc,kc)  << " " << Unextc(c,ic,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic+1,jc,kc) << " " << Unextc(c,ic+1,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic+2,jc,kc) << " " << Unextc(c,ic+2,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic+3,jc,kc) << " " << Unextc(c,ic+3,jc,kc) << endl;

}

//-----------------------------------------------------------------------
void EW::check_displacement_continuity( Sarray& Uf, Sarray& Uc, int gf, int gc )
{
   int icb, ifb, ice, ife, jcb, jfb, jce, jfe, nkf;
   
   icb = m_iStartInt[gc];
   ifb = m_iStartInt[gf];

   ice = m_iEndInt[gc];
   ife = m_iEndInt[gf];
   
   jcb = m_jStartInt[gc];
   jfb = m_jStartInt[gf];

   jce = m_jEndInt[gc];
   jfe = m_jEndInt[gf];

   nkf = m_global_nz[gf];

   float_sw4 l2err=0, l2err_global=0;
   
// for i=2*ic-1 and j=2*jc-1: Enforce continuity of displacements and normal stresses along the interface
   for( int jc= jcb ; jc <= jce ; jc++ )
      for( int ic= icb ; ic <= ice ; ic++ )
      {
// i odd, j odd
         int i=2*ic-1, j=2*jc-1;
         
         for (int c=1; c<=3; c++)
         {
            l2err += (Uc(c,ic,jc,1)-Uf(c,i,j,nkf))*(Uc(c,ic,jc,1)-Uf(c,i,j,nkf));
         }
      }
   MPI_Allreduce( &l2err, &l2err_global, 1, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

   l2err_global = sqrt(l2err_global);

   if (proc_zero())
      cout << "Fine-coarse displacement missmatch = " << l2err_global << endl;
               
   // //   cout <<"check " << Uf(c,i,j,kf) << " " << Uc(c,ic,jc,kc) << " " << Uf(c,i,j,kf)-Uc(c,ic,jc,kc) << endl;
   // //   cout << "  next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << endl;
   // i = 2*ic;
   // //   j = 2*jc;
   // double pci   = (9*(Uc(c,ic,jc,kc)+Uc(c,ic+1,jc,kc))-(Uc(c,ic-1,jc,kc)+Uc(c,ic+2,jc,kc)))/16;
   // double pcim  = (9*(Uc(c,ic,jc-1,kc)+Uc(c,ic+1,jc-1,kc))-(Uc(c,ic-1,jc-1,kc)+Uc(c,ic+2,jc-1,kc)))/16;
   // double pcip  = (9*(Uc(c,ic,jc+1,kc)+Uc(c,ic+1,jc+1,kc))-(Uc(c,ic-1,jc+1,kc)+Uc(c,ic+2,jc+1,kc)))/16;
   // double pcipp = (9*(Uc(c,ic,jc+2,kc)+Uc(c,ic+1,jc+2,kc))-(Uc(c,ic-1,jc+2,kc)+Uc(c,ic+2,jc+2,kc)))/16;
   // double pc = ( 9*(pci+pcip)-(pcim+pcipp))/16;

   // double pcj = (9*(Uc(c,ic,jc,kc)+Uc(c,ic,jc+1,kc))-(Uc(c,ic,jc-1,kc)+Uc(c,ic,jc+2,kc)))/16;
   // cout <<"check " << Uf(c,i,j,kf) << " " << Uc(c,ic,jc,kc) << " " << pci << " " << Uf(c,i,j,kf)-pci << endl;
   // //   cout << "  next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << endl;
   // double pcni   = (9*(Unextc(c,ic,jc,kc)+Unextc(c,ic+1,jc,kc))-(Unextc(c,ic-1,jc,kc)+Unextc(c,ic+2,jc,kc)))/16;
   //   cout << "  check next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << " " << pcni << " " << Unextf(c,i,j,kf)-pcni << endl;
   //   cout << " check " << Uc(c,ic-2,jc,kc) << " " << Unextc(c,ic-2,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic-1,jc,kc) << " " << Unextc(c,ic-1,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic,jc,kc)  << " " << Unextc(c,ic,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic+1,jc,kc) << " " << Unextc(c,ic+1,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic+2,jc,kc) << " " << Unextc(c,ic+2,jc,kc) << endl;
   //   cout << "       " << Uc(c,ic+3,jc,kc) << " " << Unextc(c,ic+3,jc,kc) << endl;

}

//-----------------------------------------------------------------------
void EW::dirichlet_hom_ic( Sarray& U, int g, int k, bool inner )
{
SW4_MARK_FUNCTION;
 
 using LOCAL_POL = 
   RAJA::KernelPolicy< 
   RAJA::statement::CudaKernelAsync<
     RAJA::statement::For<0, RAJA::cuda_block_exec, 
			  RAJA::statement::For<1, RAJA::cuda_block_exec, 
					       RAJA::statement::For<2, RAJA::cuda_thread_exec,
								    RAJA::statement::Lambda<0> >>>>>;
 RAJA::RangeSegment c_range(1,U.m_nc+1);
 SView &UV = U.getview();
   // zero out all ghost points
   if( !inner )
   {
     RAJA::RangeSegment jall(m_jStart[g],m_jEnd[g]+1),jzero(m_jStart[g],1);
     RAJA::RangeSegment iall(m_iStart[g],m_iEnd[g]+1),izero(m_iStart[g],1);
 
     
     // Outer layer of non-unknown ghost points
      if( m_iStartInt[g] == 1 )
      {
      // low i-side
// #pragma omp parallel for
// 	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
// 	    for( int i=m_iStart[g] ; i <= 0 ; i++ )
// 	       for( int c=1 ; c <= U.m_nc ; c++ )
		 RAJA::kernel<LOCAL_POL>(
					 RAJA::make_tuple(c_range,jall,izero),
					 [=]RAJA_DEVICE (int c,int j,int i) 
					 {
					   UV(c,i,j,k) = 0;});
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
// #pragma omp parallel for
// 	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
// 	    for( int i=m_iEndInt[g]+1 ; i <= m_iEnd[g] ; i++ )
// 	       for( int c=1 ; c <= U.m_nc ; c++ )
// 		  U(c,i,j,k) = 0;

	RAJA::RangeSegment iend(m_iEndInt[g]+1,m_iEnd[g]+1);
	RAJA::kernel<LOCAL_POL>(
				RAJA::make_tuple(c_range,jall,iend),
					 [=]RAJA_DEVICE (int c,int j,int i) 
					 {
					   UV(c,i,j,k) = 0;});
	
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
// #pragma omp parallel for
// 	 for( int j=m_jStart[g] ; j <= 0 ; j++ )
// 	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
// 	       for( int c=1 ; c <= U.m_nc ; c++ )
// 		  U(c,i,j,k) = 0;
	RAJA::kernel<LOCAL_POL>(
				RAJA::make_tuple(c_range,jzero,iall),
					 [=]RAJA_DEVICE (int c,int j,int i) 
					 {
					   UV(c,i,j,k) = 0;});


      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
// #pragma omp parallel for
// 	 for( int j=m_jEndInt[g]+1 ; j <= m_jEnd[g] ; j++ )
// 	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
// 	       for( int c=1 ; c <= U.m_nc ; c++ )
// 		  U(c,i,j,k) = 0;
	RAJA::RangeSegment jend(m_jEndInt[g]+1,m_jEnd[g]+1);
	RAJA::kernel<LOCAL_POL>(
				RAJA::make_tuple(c_range,jend,iall),
				[=]RAJA_DEVICE (int c,int j,int i) 
				{
				  UV(c,i,j,k) = 0;});
      }
   }
   else
   {
      // Interior, unknown ghost points.
      int ib, ie, jb, je;
      if( m_iStartInt[g] == 1 )
         ib = 1;
      else
	 ib = m_iStart[g];
      if( m_iEndInt[g] == m_global_nx[g] )
	 ie = m_global_nx[g];
      else
	 ie = m_iEnd[g];
      if( m_jStartInt[g] == 1 )
	 jb = 1;
      else
	 jb = m_jStart[g];
      if( m_jEndInt[g] == m_global_ny[g] )
	 je = m_global_ny[g];
      else
	 je = m_jEnd[g];
      RAJA::RangeSegment j_range(jb,je+1),i_range(ib,ie+1);
      
// #pragma omp parallel for
//       for( int j=jb ; j <= je ; j++ )
// 	 for( int i=ib ; i <= ie ; i++ )
// 	    for( int c=1 ; c <= U.m_nc ; c++ )
// 	       U(c,i,j,k) = 0;
//    }
      RAJA::kernel<LOCAL_POL>(
			      RAJA::make_tuple(c_range,j_range,i_range),
			      [=]RAJA_DEVICE (int c,int j,int i) 
			      {
				UV(c,i,j,k) = 0;}); //SYNC_STREAM;
   }
}

//-----------------------------------------------------------------------
void EW::dirichlet_twilight_ic( Sarray& U, int g, int kic, float_sw4 t )
{
SW4_MARK_FUNCTION;
   // assign exact solution at all ghost points with a given 'k'-index
   if (m_twilight_forcing)
   {
      int i1 = m_iStart[g], i2=m_iEnd[g];
      int j1 = m_jStart[g], j2=m_jEnd[g];
      int kdb=U.m_kb, kde=U.m_ke;
      float_sw4 om = m_twilight_forcing->m_omega;
      float_sw4 ph = m_twilight_forcing->m_phase;
      float_sw4 cv = m_twilight_forcing->m_c;
      float_sw4 h  = mGridSize[g];
      float_sw4* u_ptr = U.c_ptr();
      if( m_croutines )
	 twilightfortwind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
			      kdb, kde, u_ptr, t, om, cv, ph, 
			      h, m_zmin[g],
			      i1, i2, j1, j2, kic, kic );      
      else
	 twilightfortwind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
			   &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
			   &h, &m_zmin[g],
			   &i1, &i2, &j1, &j2, &kic, &kic );      
   }
   
}

//-----------------------------------------------------------------------
void EW::dirichlet_LRic( Sarray& U, int g, int kic, float_sw4 t, int adj )
{
SW4_MARK_FUNCTION;
   // Put back exact solution at the ghost points that don't participate in the interface conditions, i.e. at supergrid points
   //   int k = upper ? 0 : m_global_nz[g]+1;
   //
   // set adj= 0 for ghost pts + boundary pt
   //          1 for only ghost pts.

   int kdb=U.m_kb, kde=U.m_ke;
   SView &Uv= U.getview();
   if( !m_twilight_forcing )
   {
      if( m_iStartInt[g] == 1 )
      {
	 // low i-side
	RAJA::RangeSegment j_range(m_jStart[g] ,m_jEnd[g]+1);
	RAJA::RangeSegment i_range(m_iStart[g] ,1-adj+1);
	RAJA::RangeSegment c_range(1,U.m_nc+1);
	RAJA::kernel<RHS4_EXEC_POL_ASYNC>(
				    RAJA::make_tuple(j_range, i_range,c_range),
				    [=]RAJA_DEVICE (int j,int i,int c) {
				      Uv(c,i,j,kic)=0;});
// #pragma omp parallel for
// 	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
// 	    for( int i=m_iStart[g] ; i <= 1-adj ; i++ )
// 	       for( int c=1 ; c <= U.m_nc ; c++ )
// 		  U(c,i,j,kic) = 0;
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	RAJA::RangeSegment j_range(m_jStart[g] ,m_jEnd[g]+1);
	RAJA::RangeSegment i_range(m_iEndInt[g]+adj ,m_iEnd[g]+1);
	RAJA::RangeSegment c_range(1,U.m_nc+1);
	RAJA::kernel<RHS4_EXEC_POL_ASYNC>(
				    RAJA::make_tuple(j_range, i_range,c_range),
				    [=]RAJA_DEVICE (int j,int i,int c) {
				      Uv(c,i,j,kic)=0;});
// #pragma omp parallel for
// 	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
// 	    for( int i=m_iEndInt[g]+adj ; i <= m_iEnd[g] ; i++ )
// 	       for( int c=1 ; c <= U.m_nc ; c++ )
// 		  U(c,i,j,kic) = 0;
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	RAJA::RangeSegment j_range(m_jStart[g] ,1-adj +1);
	RAJA::RangeSegment i_range(m_iStart[g] ,m_iEnd[g] +1 );
	RAJA::RangeSegment c_range(1,U.m_nc+1);
	RAJA::kernel<RHS4_EXEC_POL_ASYNC>(
				    RAJA::make_tuple(j_range, i_range,c_range),
				    [=]RAJA_DEVICE (int j,int i,int c) {
				      Uv(c,i,j,kic)=0;});
// #pragma omp parallel for
// 	 for( int j=m_jStart[g] ; j <= 1-adj ; j++ )
// 	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
// 	       for( int c=1 ; c <= U.m_nc ; c++ )
// 		  U(c,i,j,kic) = 0;
      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	RAJA::RangeSegment j_range(m_jEndInt[g]+adj ,m_jEnd[g]+1);
	RAJA::RangeSegment i_range(m_iStart[g] ,m_iEnd[g] +1 );
	RAJA::RangeSegment c_range(1,U.m_nc+1);
	RAJA::kernel<RHS4_EXEC_POL_ASYNC>(
				    RAJA::make_tuple(j_range, i_range,c_range),
				    [=]RAJA_DEVICE (int j,int i,int c) {
				      Uv(c,i,j,kic)=0;});
// #pragma omp parallel for
// 	 for( int j=m_jEndInt[g]+adj ; j <= m_jEnd[g] ; j++ )
// 	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
// 	       for( int c=1 ; c <= U.m_nc ; c++ )
// 		  U(c,i,j,kic) = 0;
      }
   }
   else
   {
      float_sw4 om = m_twilight_forcing->m_omega;
      float_sw4 ph = m_twilight_forcing->m_phase;
      float_sw4 cv = m_twilight_forcing->m_c;
      float_sw4 h  = mGridSize[g];
      float_sw4* u_ptr = U.c_ptr();
      if( m_iStartInt[g] == 1 )
      {
	 // low i-side
	 int i1 = m_iStart[g], i2=m_iStartInt[g]-adj;
	 int j1 = m_jStart[g], j2=m_jEnd[g];
	 if( m_croutines )
	    twilightfortwind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
			      kdb, kde, u_ptr, t, om, cv, ph, 
			      h, m_zmin[g], i1, i2, j1, j2, kic, kic );
	 else
	    twilightfortwind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
			      &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
			      &h, &m_zmin[g],
			      &i1, &i2, &j1, &j2, &kic, &kic );
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	 int i1 = m_iEndInt[g]+adj, i2=m_iEnd[g];
	 int j1 = m_jStart[g], j2=m_jEnd[g];
	 if( m_croutines )
	    twilightfortwind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
				 kdb, kde, u_ptr, t, om, cv, ph, 
				 h, m_zmin[g], i1, i2, j1, j2, kic, kic );
	 else
	    twilightfortwind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
			      &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
			      &h, &m_zmin[g],
			      &i1, &i2, &j1, &j2, &kic, &kic );
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	 int i1 = m_iStart[g], i2=m_iEnd[g];
	 int j1 = m_jStart[g], j2=m_jStartInt[g]-adj;
	 if( m_croutines )
	    twilightfortwind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
				 kdb, kde, u_ptr, t, om, cv, ph, 
				 h, m_zmin[g], i1, i2, j1, j2, kic, kic );
	 else
	    twilightfortwind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
			      &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
			      &h, &m_zmin[g],
			      &i1, &i2, &j1, &j2, &kic, &kic );
      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	 int i1 = m_iStart[g], i2=m_iEnd[g];
	 int j1 = m_jEndInt[g]+adj, j2=m_jEnd[g];
	 if( m_croutines )
	    twilightfortwind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
				 kdb, kde, u_ptr, t, om, cv, ph, 
				 h, m_zmin[g],i1, i2, j1, j2, kic, kic );
	 else
	    twilightfortwind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
			      &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
			      &h, &m_zmin[g],
			      &i1, &i2, &j1, &j2, &kic, &kic );
      }
   }      
}

//-----------------------------------------------------------------------
void EW::dirichlet_LRstress( Sarray& B, int g, int kic, float_sw4 t, int adj )
{
SW4_MARK_FUNCTION;
   // Exact stresses at the ghost points that don't participate in the interface conditions,
   // i.e. at supergrid (Dirichlet) points
   //   int k = upper ? 0 : m_global_nz[g]+1;
   //
   // set adj= 0 for ghost pts + boundary pt
   //          1 for only ghost pts.


   //int kdb=B.m_kb, kde=B.m_ke;

//   printf("dirichlet_LRstress> kdb=%d, kde=%d\n", kdb, kde);
 SView &Bv= B.getview();
 //SYNC_STREAM; // Since this is running on the host
   if( !m_twilight_forcing )
   {
      if( m_iStartInt[g] == 1 )
      {
	 // low i-side
	RAJA::RangeSegment j_range(m_jStart[g] ,m_jEnd[g]+1);
	RAJA::RangeSegment i_range(m_iStart[g] ,1-adj+1);
	RAJA::RangeSegment c_range(1,B.m_nc+1);
	RAJA::kernel<RHS4_EXEC_POL_ASYNC>(
				    RAJA::make_tuple(j_range, i_range,c_range),
				    [=]RAJA_DEVICE (int j,int i,int c) {
				      Bv(c,i,j,kic)=0;});
#// pragma omp parallel for
 // 	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
 // 	    for( int i=m_iStart[g] ; i <= 1-adj ; i++ )
 // 	       for( int c=1 ; c <= B.m_nc ; c++ )
 // 		  B(c,i,j,kic) = 0;
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	RAJA::RangeSegment j_range(m_jStart[g] ,m_jEnd[g]+1);
	RAJA::RangeSegment i_range(m_iEndInt[g]+adj ,m_iEnd[g]+1);
	RAJA::RangeSegment c_range(1,B.m_nc+1);
	RAJA::kernel<RHS4_EXEC_POL_ASYNC>(
				    RAJA::make_tuple(j_range, i_range,c_range),
				    [=]RAJA_DEVICE (int j,int i,int c) {
				      Bv(c,i,j,kic)=0;});
#// pragma omp parallel for
 // 	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
 // 	    for( int i=m_iEndInt[g]+adj ; i <= m_iEnd[g] ; i++ )
 // 	       for( int c=1 ; c <= B.m_nc ; c++ )
 // 		  B(c,i,j,kic) = 0;
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	RAJA::RangeSegment j_range(m_jStart[g] ,1-adj +1);
	RAJA::RangeSegment i_range(m_iStart[g] ,m_iEnd[g] +1 );
	RAJA::RangeSegment c_range(1,B.m_nc+1);
	RAJA::kernel<RHS4_EXEC_POL_ASYNC>(
				    RAJA::make_tuple(j_range, i_range,c_range),
				    [=]RAJA_DEVICE (int j,int i,int c) {
				      Bv(c,i,j,kic)=0;});
// 	 for( int j=m_jStart[g] ; j <= 1-adj ; j++ )
// #pragma omp parallel for
// 	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
// 	       for( int c=1 ; c <= B.m_nc ; c++ )
// 		  B(c,i,j,kic) = 0;
      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	RAJA::RangeSegment j_range(m_jEndInt[g]+adj ,m_jEnd[g]+1);
	RAJA::RangeSegment i_range(m_iStart[g] ,m_iEnd[g] +1 );
	RAJA::RangeSegment c_range(1,B.m_nc+1);
	RAJA::kernel<RHS4_EXEC_POL_ASYNC>(
				    RAJA::make_tuple(j_range, i_range,c_range),
				    [=]RAJA_DEVICE (int j,int i,int c) {
				      Bv(c,i,j,kic)=0;});
// 	 for( int j=m_jEndInt[g]+adj ; j <= m_jEnd[g] ; j++ )
// #pragma omp parallel for
// 	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
// 	       for( int c=1 ; c <= B.m_nc ; c++ )
// 		  B(c,i,j,kic) = 0;
      }
   }
   else // twilight forcing below
   {
// dbg
      if (false && g==0)
      {
         int i=0, j=25, k=1;
         float_sw4 x=(i-1)*mGridSize[g];
         float_sw4 y=(j-1)*mGridSize[g];
         float_sw4 z=(k-1)*mGridSize[g] + m_zmin[g];

         printf("3: x=%e, y=%e, z=%e, mu=%e\n", x, y, z, mMu[g](i,j,k));
      }
// get array pointers for fortran
      float_sw4* mu_ptr    = mMu[g].c_ptr();
      float_sw4* la_ptr    = mLambda[g].c_ptr();
      float_sw4 om = m_twilight_forcing->m_omega;
      float_sw4 ph = m_twilight_forcing->m_phase;
      float_sw4 cv = m_twilight_forcing->m_c;
      float_sw4 h  = mGridSize[g];
      float_sw4* b_ptr = B.c_ptr();
      float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
      float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();

      float_sw4* mua_ptr = NULL;
      float_sw4* laa_ptr   = NULL;
      if (m_use_attenuation)
      {
         mua_ptr    = mMuVE[g][0].c_ptr();
         laa_ptr    = mLambdaVE[g][0].c_ptr();
      }

      if( m_iStartInt[g] == 1 )
      {
	 // low i-side
	 int i1 = m_iStart[g], i2=m_iStartInt[g]-adj;
	 int j1 = m_jStart[g], j2=m_jEnd[g];
         if( usingSupergrid() )
         {
// assigns B
	    if( m_croutines )
	       twfrsurfzsg_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
				    h, kic, t, om, cv, ph, omstrx, omstry,
				    b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
	    else
	       twfrsurfzsg_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
				 h, kic, t, om, cv, ph, omstrx, omstry,
				 b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
            if (m_use_attenuation) // only 1 mechanism with twilight forcing
            {
// adds attenuation to B
	       if( m_croutines )
		  twfrsurfzsg_att_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                     h, kic, t, om, cv, ph, omstrx, omstry,
                                     b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
	       else
		  twfrsurfzsg_att_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                     h, kic, t, om, cv, ph, omstrx, omstry,
                                     b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
            }
         }
         else
         {
// assigns B
	    if( m_croutines )
	       twfrsurfz_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
			       h, kic, t, om, cv, ph, b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
	    else
	       twfrsurfz_wind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
			       &h, &kic, &t, &om, &cv, &ph, b_ptr, mu_ptr, la_ptr, &m_zmin[g], &i1, &i2, &j1, &j2 );
            if (m_use_attenuation) // only 1 mechanism with twilight forcing
            {
// adds attenuation to B
	       if( m_croutines )
		  twfrsurfz_att_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                   h, kic, t, om, cv, ph, b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
	       else
		  twfrsurfz_att_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                   h, kic, t, om, cv, ph, b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
            }
         } // else (not supergrid)
      } //  if( m_iStartInt[g] == 1 )
      
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	 int i1 = m_iEndInt[g]+adj, i2=m_iEnd[g];
	 int j1 = m_jStart[g], j2=m_jEnd[g];
         if( usingSupergrid() )
         {
// assigns B
	    if( m_croutines )
	       twfrsurfzsg_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
	    else
	       twfrsurfzsg_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
            if (m_use_attenuation) // only 1 mechanism with twilight forcing
            {
// adds attenuation to B
	       if( m_croutines )
		  twfrsurfzsg_att_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                     h, kic, t, om, cv, ph, omstrx, omstry,
                                     b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
	       else
		  twfrsurfzsg_att_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                     h, kic, t, om, cv, ph, omstrx, omstry,
                                     b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
            }
         }
         else
         {
// assigns B
	    if( m_croutines )
	       twfrsurfz_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                            h, kic, t, om, cv, ph, b_ptr, 
                            mu_ptr, la_ptr, m_zmin[g],
                            i1, i2, j1, j2 );
	    else
	       twfrsurfz_wind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
                            &h, &kic, &t, &om, &cv, &ph, b_ptr, 
                            mu_ptr, la_ptr, &m_zmin[g],
                            &i1, &i2, &j1, &j2 );
            if (m_use_attenuation) // only 1 mechanism with twilight forcing
            {
// adds attenuation to B
	       if( m_croutines )
		  twfrsurfz_att_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                   h, kic, t, om, cv, ph, b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
	       else
		  twfrsurfz_att_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                   h, kic, t, om, cv, ph, b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
            }
         } // else (not supergrid)
      } // end if( m_iEndInt[g] == m_global_nx[g] )
      
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	 int i1 = m_iStart[g], i2=m_iEnd[g];
	 int j1 = m_jStart[g], j2=m_jStartInt[g]-adj;
         if( usingSupergrid() )
         {
// assigns B
	    if( m_croutines )
	       twfrsurfzsg_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
	    else
	       twfrsurfzsg_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
            if (m_use_attenuation) // only 1 mechanism with twilight forcing
            {
// adds attenuation to B
	       if( m_croutines )
		  twfrsurfzsg_att_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                     h, kic, t, om, cv, ph, omstrx, omstry,
                                     b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
	       else
		  twfrsurfzsg_att_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                     h, kic, t, om, cv, ph, omstrx, omstry,
                                     b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
            }
         }
         else
         {
// assigns B
	    if( m_croutines )
	       twfrsurfz_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                            h, kic, t, om, cv, ph, b_ptr, 
                            mu_ptr, la_ptr, m_zmin[g],
                            i1, i2, j1, j2 );
	    else
	       twfrsurfz_wind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
                            &h, &kic, &t, &om, &cv, &ph, b_ptr, 
                            mu_ptr, la_ptr, &m_zmin[g],
                            &i1, &i2, &j1, &j2 );
            if (m_use_attenuation) // only 1 mechanism with twilight forcing
            {
// adds attenuation to B
	       if( m_croutines )
		  twfrsurfz_att_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                   h, kic, t, om, cv, ph, b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
	       else
		  twfrsurfz_att_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                   h, kic, t, om, cv, ph, b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
            }
         } // else (not supergrid)
      } // end if( m_jStartInt[g] == 1 )
      
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	 int i1 = m_iStart[g], i2=m_iEnd[g];
	 int j1 = m_jEndInt[g]+adj, j2=m_jEnd[g];
         if( usingSupergrid() )
         {
// assigns B
	    if( m_croutines )
	       twfrsurfzsg_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
	    else
	       twfrsurfzsg_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
            if (m_use_attenuation) // only 1 mechanism with twilight forcing
            {
// adds attenuation to B
	       if( m_croutines )
		  twfrsurfzsg_att_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                     h, kic, t, om, cv, ph, omstrx, omstry,
                                     b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
	       else
		  twfrsurfzsg_att_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                     h, kic, t, om, cv, ph, omstrx, omstry,
                                     b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
            }
         }
         else
         {
// assigns B
	    if( m_croutines )
	       twfrsurfz_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                            h, kic, t, om, cv, ph, b_ptr, 
                            mu_ptr, la_ptr, m_zmin[g],
                            i1, i2, j1, j2 );
	    else
	       twfrsurfz_wind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
                            &h, &kic, &t, &om, &cv, &ph, b_ptr, 
                            mu_ptr, la_ptr, &m_zmin[g],
                            &i1, &i2, &j1, &j2 );
            if (m_use_attenuation) // only 1 mechanism with twilight forcing
            {
// adds attenuation to B
	       if( m_croutines )
		  twfrsurfz_att_wind_ci( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                   h, kic, t, om, cv, ph, b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
	       else
		  twfrsurfz_att_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                   h, kic, t, om, cv, ph, b_ptr, mua_ptr, laa_ptr, m_zmin[g], i1, i2, j1, j2 );
            }
         }
      } // end if( m_jEndInt[g] == m_global_ny[g] )
      
   } // else twilight
   
}

//-----------------------------------------------------------------------
void EW::gridref_initial_guess( Sarray& u, int g, bool upper )
{
SW4_MARK_FUNCTION;
// Extrapolate the initial guess from neighboring point.
   int k, s;
   if( upper )
   {
      k = 0;
      s = 1;
   }
   else
   {
      k = m_kEndInt[g]+1;
      s = -1;
   }
   int ib, ie, jb, je;
   if( m_iStartInt[g] == 1 )
      ib = 1;
   else
      ib = m_iStart[g];
   if( m_iEndInt[g] == m_global_nx[g] )
      ie = m_global_nx[g];
   else
      ie = m_iEnd[g];
   if( m_jStartInt[g] == 1 )
      jb = 1;
   else
      jb = m_jStart[g];
   if( m_jEndInt[g] == m_global_ny[g] )
      je = m_global_ny[g];
   else
      je = m_jEnd[g];
   
   SView &uV = u.getview();
   
   RAJA::RangeSegment j_range(jb,je+1),i_range(ib,ie+1);
   
   // for( int j=jb ; j <= je ; j++ )
   //    for( int i=ib ; i <= ie ; i++ )
   using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
  RAJA::statement::For<0, RAJA::cuda_block_exec, 
  RAJA::statement::For<1, RAJA::cuda_thread_exec,
  RAJA::statement::Lambda<0> >>>>;

   RAJA::kernel<LOCAL_POL>(
			   RAJA::make_tuple(j_range,i_range),
			   [=]RAJA_DEVICE (int j,int i) 
			   {
			     uV(1,i,j,k) = uV(1,i,j,k+s);
			     uV(2,i,j,k) = uV(2,i,j,k+s);
			     uV(3,i,j,k) = uV(3,i,j,k+s);
			   }); //SYNC_STREAM;
}

//-----------------------------------------------------------------------
void EW::compute_preliminary_corrector( Sarray& a_Up, Sarray& a_U, Sarray& a_Um,
                                        Sarray* a_AlphaVEp, Sarray* a_AlphaVE, Sarray* a_AlphaVEm,
					Sarray& Utt, Sarray& Unext, int g, int kic, float_sw4 t, 
					Sarray &Ftt, vector<GridPointSource*> &point_sources )
{
SW4_MARK_FUNCTION;
   //
   // NOTE: This routine is called by enforceIC() after the predictor stage to calculate the interior contribution to
   // the corrector on the interface at time step.
   // It is NOT called by enforceIC2(), which handles 2nd order time stepping.
   // Super-grid dissipation is added because it is part of the interior contribution to the corrected displacement.
   //
   float_sw4 idt2 = 1/(mDt*mDt);
   // to evaluate L(Up_tt) for k=kic, we need Up(k) in the vicinity of the interface
   // Note: Utt is needed at all points (interior + ghost) to evaluate L(Utt) in all interior points

   Utt.prefetch();
   a_U.prefetch();
   a_Up.prefetch();
   a_Um.prefetch();
   if (m_croutines) // optimized C-version for reversed index ordering
     dpdmt_wind( Utt.m_ib, Utt.m_ie, Utt.m_jb, Utt.m_je, Utt.m_kb, Utt.m_ke, a_U.m_kb, a_U.m_ke,
		 a_Up.c_ptr(), a_U.c_ptr(), a_Um.c_ptr(), Utt.c_ptr(), idt2 );
   else
   {
     for( int k=Utt.m_kb ; k <= Utt.m_ke ; k++ )
       for( int j=Utt.m_jb ; j <= Utt.m_je ; j++ )
   	 for( int i=Utt.m_ib ; i <= Utt.m_ie ; i++ )
   	 {
   	    Utt(1,i,j,k) = idt2*(a_Up(1,i,j,k)-2*a_U(1,i,j,k)+a_Um(1,i,j,k));
   	    Utt(2,i,j,k) = idt2*(a_Up(2,i,j,k)-2*a_U(2,i,j,k)+a_Um(2,i,j,k));
   	    Utt(3,i,j,k) = idt2*(a_Up(3,i,j,k)-2*a_U(3,i,j,k)+a_Um(3,i,j,k));
   	 }
   }
   

// all points (for mMu, mLambda

   int ib=m_iStart[g], jb=m_jStart[g], kb=m_kStart[g];
   int ie=m_iEnd[g], je=m_jEnd[g], ke=m_kEnd[g];
// k-indices for Utt
   int kbu = Utt.m_kb, keu= Utt.m_ke;

   // Compute L(Utt) at k=kic.
   char op='=';
   int nz = m_global_nz[g];
   Sarray Lutt(3,ib,ie,jb,je,kic,kic,__FILE__,__LINE__);
   Lutt.set_to_zero_async(); // Keep memory checker happy
// Note: 6 first arguments of the function call:
// (ib,ie), (jb,je), (kb,ke) is the declared size of mMu and mLambda in the (i,j,k)-directions, respectively

// there are C and Fortran versions of this routine that are selected by the Makefile
   rhs4th3wind( ib, ie, jb, je, kb, ke, nz, m_onesided[g], m_acof,
		       m_bope, m_ghcof, Lutt.c_ptr(), Utt.c_ptr(), mMu[g].c_ptr(),
		       mLambda[g].c_ptr(), mGridSize[g], m_sg_str_x[g], m_sg_str_y[g],
		       m_sg_str_z[g], op, kbu, keu, kic, kic );

// Plan: 1) loop over all mechanisms, 2) precompute d^2 alpha/dt^2 and store in Utt, 3) accumulate visco-elastic stresses
   if( m_use_attenuation && m_number_mechanisms > 0 )
   {
      op = '-'; // Subtract Lu := Lu - L_a(alpha)
      for( int a=0 ; a < m_number_mechanisms ; a++ )
      {
// compute corrected memory variable for mechanism 'a', near the interface, store in 'Utt'
         updateMemVarCorrNearInterface( Utt, a_AlphaVEm[a], a_Up, a_U, a_Um, t, a, g );
         
// assume a_U and a_AlphaVE have the same dimensions for all mechanisms
	 SW4_MARK_BEGIN("PRED_DEVICE_LOOP");
	 SView &UttV = Utt.getview();
	 SView &a_AlphaVEV = a_AlphaVE[a].getview();
	 SView &a_AlphaVEmV = a_AlphaVEm[a].getview();
	 RAJA::RangeSegment k_range(Utt.m_kb, Utt.m_ke +1);
	 RAJA::RangeSegment j_range(Utt.m_jb, Utt.m_je+1);
	 RAJA::RangeSegment i_range(Utt.m_ib, Utt.m_ie+1 );
         // for( int k=Utt.m_kb ; k <= Utt.m_ke ; k++ )
         //    for( int j=Utt.m_jb ; j <= Utt.m_je ; j++ )
         //       for( int i=Utt.m_ib ; i <= Utt.m_ie ; i++ )
	  RAJA::kernel<XRHS_POL_ASYNC>(
				 RAJA::make_tuple(k_range, j_range,i_range),
				 [=]RAJA_DEVICE (int k,int j,int i)
               {
// corrector value of AlphaVEp in variable 'Utt'
                  UttV(1,i,j,k) = idt2*(UttV(1,i,j,k)-2*a_AlphaVEV(1,i,j,k)+a_AlphaVEmV(1,i,j,k));
                  UttV(2,i,j,k) = idt2*(UttV(2,i,j,k)-2*a_AlphaVEV(2,i,j,k)+a_AlphaVEmV(2,i,j,k));
                  UttV(3,i,j,k) = idt2*(UttV(3,i,j,k)-2*a_AlphaVEV(3,i,j,k)+a_AlphaVEmV(3,i,j,k));
               }); //SYNC_STREAM;
	 SW4_MARK_END("PRED_DEVICE_LOOP");
// NEW June 13, 2017: add in visco-elastic terms
         float_sw4* mua_ptr = mMuVE[g][a].c_ptr();
         float_sw4* lambdaa_ptr = mLambdaVE[g][a].c_ptr();
	 rhs4th3wind( ib, ie, jb, je, kb, ke, nz, m_onesided[g], m_acof_no_gp, 
			     m_bope, m_ghcof_no_gp, Lutt.c_ptr(), Utt.c_ptr(), mua_ptr, //use stencil WITHOUT ghost points
			     lambdaa_ptr, mGridSize[g], m_sg_str_x[g], m_sg_str_y[g],
			     m_sg_str_z[g], op, kbu, keu, kic, kic );
      } // end for a      
   } // end if using attenuation
   
   // Compute forcing_{tt} at k=kic
//    Sarray force(3,ib,ie,jb,je,kic,kic);
//    if( m_twilight_forcing )
//    {
//       float_sw4 om = m_twilight_forcing->m_omega;
//       float_sw4 ph = m_twilight_forcing->m_phase;
//       float_sw4 cv = m_twilight_forcing->m_c;
//       float_sw4 omm= m_twilight_forcing->m_momega;
//       float_sw4 phm= m_twilight_forcing->m_mphase;
//       float_sw4 amprho   = m_twilight_forcing->m_amprho;
//       float_sw4 ampmu    = m_twilight_forcing->m_ampmu;
//       float_sw4 amplambda= m_twilight_forcing->m_amplambda;
//       if( m_croutines )
// 	 forcingttfort_ci( ib, ie, jb, je, kic, kic, force.c_ptr(), t, om, cv, ph, omm, phm,
// 			   amprho, ampmu, amplambda, mGridSize[g], m_zmin[g] );
//       else
// 	 forcingttfort( &ib, &ie, &jb, &je, &kic, &kic, force.c_ptr(), &t, &om, &cv, &ph, &omm, &phm,
// 			&amprho, &ampmu, &amplambda, &mGridSize[g], &m_zmin[g] );
//    }
//    else if( m_rayleigh_wave_test || m_energy_test )
//       force.set_to_zero();
//    else
//    {
//      // Default: m_point_source_test, m_lamb_test or full seismic case
//       force.set_to_zero();
// // NOTE: this routine needs to be reworked!
// // AP: Can we do omp for around this loop?
//       for( int s = 0 ; s < point_sources.size() ; s++ )
//       {
// 	 if( point_sources[s]->m_grid == g && point_sources[s]->m_k0 == kic )
// 	 {
// 	    float_sw4 fxyz[3];
// 	    point_sources[s]->getFxyztt(t,fxyz);
// 	    force(1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
// 	    force(2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
// 	    force(3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
// 	 }
//       }
//    }

   float_sw4 cof = mDt*mDt*mDt*mDt/12.0;

   // if (m_croutines)
   //   update_unext( ib, ie, jb, je, kb, ke, Unext.c_ptr(), a_Up.c_ptr(), Lutt.c_ptr(), force.c_ptr(), 
   //      	   mRho[g].c_ptr(), cof, kic);
   // else
   {
// #pragma omp parallel for
//      for( int j=Unext.m_jb+2 ; j <= Unext.m_je-2 ; j++ )
// #pragma omp simd
//        for( int i=Unext.m_ib+2 ; i <= Unext.m_ie-2 ; i++ )
//        {
     RAJA::RangeSegment j_range(Unext.m_jb+2,Unext.m_je-1);
     RAJA::RangeSegment i_range(Unext.m_ib+2,Unext.m_ie-1 );
     SView &mRhogV = mRho[g].getview();
     SView &UnextV =  Unext.getview();
     SView &a_UpV = a_Up.getview();
     SView &LuttV = Lutt.getview();
     SView &FttV = Ftt.getview();
     RAJA::kernel<PRELIM_CORR_EXEC_POL1_ASYNC>(
					 RAJA::make_tuple(j_range,i_range),
					 [=]RAJA_DEVICE (int j,int i) {
					   float_sw4 irho=cof/mRhogV(i,j,kic);
					   UnextV(1,i,j,kic) = a_UpV(1,i,j,kic) + irho*(LuttV(1,i,j,kic)+FttV(1,i,j,kic)); // +force(1,i,j,kic));
					   UnextV(2,i,j,kic) = a_UpV(2,i,j,kic) + irho*(LuttV(2,i,j,kic)+FttV(2,i,j,kic)); // +force(2,i,j,kic));
					   UnextV(3,i,j,kic) = a_UpV(3,i,j,kic) + irho*(LuttV(3,i,j,kic)+FttV(3,i,j,kic));//+force(3,i,j,kic));
					 }); SYNC_STREAM;
   }

// add in super-grid damping terms (does it make a difference?)
   if (usingSupergrid()) // Assume 4th order AD, Cartesian grid
   {
// assign array pointers on the fly
      // what time levels of U should be used here? NOTE: t_n and t_{n-1} ?  Up is the predictor for U(t_{n+1})
// July 22: Changed time levels for (Up, U) to (U, Um)
//      addsg4wind( &mDt, &mGridSize[g], Unext.c_ptr(), a_Up.c_ptr(), a_U.c_ptr(), mRho[g].c_ptr(),
      if( m_croutines )
	 addsg4wind_ci( Unext.c_ptr(), a_U.c_ptr(), a_Um.c_ptr(), mRho[g].c_ptr(),
			m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g],
			m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
			ib, ie, jb, je, kb, ke, m_supergrid_damping_coefficient, Unext.m_kb, Unext.m_ke, kic, kic );
      else
	 addsg4wind( &mDt, &mGridSize[g], Unext.c_ptr(), a_U.c_ptr(), a_Um.c_ptr(), mRho[g].c_ptr(),
		     m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g],
		     m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
		     ib, ie, jb, je, kb, ke, m_supergrid_damping_coefficient, Unext.m_kb, Unext.m_ke, kic, kic );
      // Note: the last four arguments define the declared size of Unext, followed by the lower and upper boundaries of the k-window
   }
} //end compute_preliminary_corrector


//-----------------------------------------------------------------------
void EW::compute_preliminary_predictor( Sarray& a_Up, Sarray& a_U, Sarray* a_AlphaVEp, Sarray& Unext,
					int g, int kic, float_sw4 t, Sarray &F, vector<GridPointSource*> &point_sources )
{
  SW4_MARK_FUNCTION;
   //
   // NOTE: This routine is called by enforceIC() after the corrector stage to calculate the interior contribution to
   // the predictor on the interface at next time step.
   // It is also called by enforceIC2(), which handles 2nd order time stepping. Only in that case should super-grid dissipation be added.
   //
   int ib=m_iStart[g], jb=m_jStart[g], kb=m_kStart[g];
   int ie=m_iEnd[g], je=m_jEnd[g], ke=m_kEnd[g];

   // Compute L(Up) at k=kic.
   Sarray Lu(3,ib,ie,jb,je,kic,kic,__FILE__,__LINE__);
   Lu.set_to_zero_async();  // Keep memory checker happy
   char op='=';
   int nz = m_global_nz[g];
// Note: 6 first arguments of the function call:
// (ib,ie), (jb,je), (kb,ke) is the declared size of mMu and mLambda in the (i,j,k)-directions, respectively
   rhs4th3wind( ib, ie, jb, je, kb, ke, nz, m_onesided[g], m_acof, // ghost point operators for elastic part
		       m_bope, m_ghcof, Lu.c_ptr(), a_Up.c_ptr(), mMu[g].c_ptr(),
		       mLambda[g].c_ptr(), mGridSize[g], m_sg_str_x[g], m_sg_str_y[g],
		       m_sg_str_z[g], op, kb, ke, kic, kic ); 
// Note: 4 last arguments of the above function call:
// (kb,ke) is the declared size of Up in the k-direction
// (kic,kic) is the declared size of Lu in the k-direction

   // add in visco-elastic terms
   if( m_use_attenuation && m_number_mechanisms > 0 )
   {
      op = '-'; // Subtract Lu := Lu - L_a(alpha)
      for( int a=0 ; a < m_number_mechanisms ; a++ )
      {
         float_sw4* alpha_ptr = a_AlphaVEp[a].c_ptr();
         float_sw4* mua_ptr = mMuVE[g][a].c_ptr();
         float_sw4* lambdaa_ptr = mLambdaVE[g][a].c_ptr();
	 rhs4th3wind( ib, ie, jb, je, kb, ke, nz, m_onesided[g], m_acof_no_gp, // NO ghost points for attenuation
			     m_bope, m_ghcof_no_gp, Lu.c_ptr(), alpha_ptr, mua_ptr,
			     lambdaa_ptr, mGridSize[g], m_sg_str_x[g], m_sg_str_y[g],
			     m_sg_str_z[g], op, kb, ke, kic, kic );
      }
   }
   
   // Compute forcing at k=kic
   // Sarray f(3,ib,ie,jb,je,kic,kic);
   // if( m_twilight_forcing )
   // {
   //    float_sw4 om = m_twilight_forcing->m_omega;
   //    float_sw4 ph = m_twilight_forcing->m_phase;
   //    float_sw4 cv = m_twilight_forcing->m_c;
   //    float_sw4 omm= m_twilight_forcing->m_momega;
   //    float_sw4 phm= m_twilight_forcing->m_mphase;
   //    float_sw4 amprho=m_twilight_forcing->m_amprho;
   //    float_sw4 ampmu=m_twilight_forcing->m_ampmu;
   //    float_sw4 amplambda=m_twilight_forcing->m_amplambda;
   //    if( usingSupergrid() )
   //    {
   //       float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
   //       float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
   //       float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
   //       if( m_croutines )
   //          forcingfortsg_ci( ib, ie, jb, je, kic, kic, f.c_ptr(), t, om, cv, ph, omm, phm,
   //      		      amprho, ampmu, amplambda, mGridSize[g], m_zmin[g],
   //      		      omstrx, omstry, omstrz );
   //       else
   //          forcingfortsg(  &ib, &ie, &jb, &je, &kic, &kic, f.c_ptr(), &t, &om, &cv, &ph, &omm, &phm,
   //      		    &amprho, &ampmu, &amplambda, &mGridSize[g], &m_zmin[g],
   //      		    &omstrx, &omstry, &omstrz );
   //       if( m_use_attenuation ) // NOTE: forcingfortsgatt only adds in the visco-elastic terms to 'f'
   //       {
   //          if( m_croutines )
   //             forcingfortsgatt_ci( ib, ie, jb, je, kic, kic, f.c_ptr(), t, om, cv, ph, omm, phm,
   //      			    amprho, ampmu, amplambda, mGridSize[g], m_zmin[g],
   //      			    omstrx, omstry, omstrz );
   //          else
   //             forcingfortsgatt(  &ib, &ie, &jb, &je, &kic, &kic, f.c_ptr(), &t, &om, &cv, &ph, &omm, &phm,
   //      			  &amprho, &ampmu, &amplambda, &mGridSize[g], &m_zmin[g],
   //      			  &omstrx, &omstry, &omstrz );
   //       }
   //    }
   //    else
   //    {
   //       if( m_croutines )
   //          forcingfort_ci( ib, ie, jb, je, kic, kic, f.c_ptr(), t, om, cv, ph, omm, phm,
   //      		    amprho, ampmu, amplambda, mGridSize[g], m_zmin[g] );
   //       else
   //          forcingfort( &ib, &ie, &jb, &je, &kic, &kic, f.c_ptr(), &t, &om, &cv, &ph, &omm, &phm,
   //      		 &amprho, &ampmu, &amplambda, &mGridSize[g], &m_zmin[g] );
   //    }
   // } // end twilight
   // else if( m_rayleigh_wave_test || m_energy_test )
   //    f.set_to_zero();
   // else
   // {
   //   // Default: m_point_source_test, m_lamb_test or full seismic case
   //    f.set_to_zero();
   //    for( int s = 0 ; s < point_sources.size() ; s++ )
   //    {
   //       if( point_sources[s]->m_grid == g && point_sources[s]->m_k0 == kic )
   //       {
   //          float_sw4 fxyz[3];
   //          point_sources[s]->getFxyz(t,fxyz);
   //          f(1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
   //          f(2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
   //          f(3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
   //       }
   //    }
   // } // end else... (not twilight, rayleigh_test, or energy_test)
   
   const float_sw4 cof = mDt*mDt;
// initialize
   Unext.set_to_zero_async();
   SView &UnextV = Unext.getview();
   SView &a_UpV = a_Up.getview();
   SView &a_UV = a_U.getview();
   SView &LuV = Lu.getview();
   SView &FV = F.getview();
   SView &mRhogV = mRho[g].getview();

   // SView &UnextV = *new SView(Unext);
   // SView &a_UpV = *new SView(a_Up);
   // SView &a_UV = *new SView(a_U);
   // SView &LuV = *new SView(Lu);
   // SView &FV = *new SView(F);
   // SView &mRhogV = *new SView(mRho[g]);
   RAJA::RangeSegment j_range(jb+2,je-1);
   RAJA::RangeSegment i_range(ib+2,ie-1);
   RAJA::kernel<PRELIM_PRED_EXEC_POL1_ASYNC>(
			RAJA::make_tuple(j_range,i_range),
			[=]RAJA_DEVICE (int j,int i) {
#// pragma omp parallel for
//    for( int j=jb+2 ; j <= je-2 ; j++ )
// #pragma omp simd
//       for( int i=ib+2 ; i <= ie-2 ; i++ )
//       {
	 float_sw4 irho=cof/mRhogV(i,j,kic);
	 UnextV(1,i,j,kic) = 2*a_UpV(1,i,j,kic) - a_UV(1,i,j,kic) + irho*(LuV(1,i,j,kic)+FV(1,i,j,kic)); //+f(1,i,j,kic));
	 UnextV(2,i,j,kic) = 2*a_UpV(2,i,j,kic) - a_UV(2,i,j,kic) + irho*(LuV(2,i,j,kic)+FV(2,i,j,kic)); //+f(2,i,j,kic));
	 UnextV(3,i,j,kic) = 2*a_UpV(3,i,j,kic) - a_UV(3,i,j,kic) + irho*(LuV(3,i,j,kic)+FV(3,i,j,kic)); //+f(3,i,j,kic));
			}); //SYNC_STREAM;
// add in super-grid damping terms
   if (mOrder==2 && usingSupergrid()) // only needed for 2nd order time-stepping. Assume 4th order AD, Cartesian grid
   {
// assign array pointers on the fly
      if( m_croutines )
	 addsg4wind_ci( Unext.c_ptr(), a_Up.c_ptr(), a_U.c_ptr(), mRho[g].c_ptr(),
			m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g],
			m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
			ib, ie, jb, je, kb, ke, m_supergrid_damping_coefficient, Unext.m_kb, Unext.m_ke, kic, kic );
      else
	 addsg4wind( &mDt, &mGridSize[g], Unext.c_ptr(), a_Up.c_ptr(), a_U.c_ptr(), mRho[g].c_ptr(),
		     m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g],
		     m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
		     ib, ie, jb, je, kb, ke, m_supergrid_damping_coefficient, Unext.m_kb, Unext.m_ke, kic, kic );
      // Note: the last four arguments define the declared size of Unext, followed by the lower and upper boundaries of the k-window
   }
}

//-----------------------------------------------------------------------
void EW::compute_icstresses( Sarray& a_Up, Sarray& B, int g, int kic,
			     float_sw4* a_str_x, float_sw4* a_str_y )
{
SW4_MARK_FUNCTION;
   const float_sw4 a1=2.0/3, a2=-1.0/12;
   bool upper = (kic == 1);
   int k=kic;
   float_sw4 ih = 1/mGridSize[g];
   int ifirst = a_Up.m_ib;
   int jfirst = a_Up.m_jb;
#define str_x(i) a_str_x[(i-ifirst)]   
#define str_y(j) a_str_y[(j-jfirst)]
   //Sarray &mMug = mMu[g];
   //Sarray &mLambdag = mLambda[g];
   //float_sw4 *Up_data = a_Up.c_ptr();
   //size_t m_base,m_offc,m_offi,m_offj,m_offk;
   //m_base = a_Up.m_base;
   //m_offc = a_Up.m_offc;
   //m_offi = a_Up.m_offi;
   //m_offj = a_Up.m_offj;
   // m_offk = a_Up.m_offk;
   
   ASSERT_MANAGED(a_Up.c_ptr());
   ASSERT_MANAGED(B.c_ptr());
   ASSERT_MANAGED( mMu[g].c_ptr());
   ASSERT_MANAGED( mLambda[g].c_ptr());
   ASSERT_MANAGED(m_sbop);
   ASSERT_MANAGED(a_str_x);
   ASSERT_MANAGED(a_str_y);
   SView &UpV = a_Up.getview();
   SView &BV = B.getview();
   SView &mMuV = mMu[g].getview();
   SView &mLambdaV = mLambda[g].getview();

   int istart,jstart;
   istart = B.m_ib+2;
   jstart = B.m_jb+2;
   
   a_Up.prefetch();
   B.prefetch();
   mMu[g].prefetch();
   mLambda[g].prefetch();
   float_sw4* lm_sbop = m_sbop;
   RAJA::RangeSegment j_range(B.m_jb+2,B.m_je-1);
   RAJA::RangeSegment i_range(B.m_ib+2,B.m_ie-1);
   
   RAJA::kernel<ICSTRESS_EXEC_POL_ASYNC>(
			RAJA::make_tuple(j_range,i_range),
			[=]RAJA_DEVICE (int j,int i) {

			  //#pragma omp parallel for
			  //   for( int j=B.m_jb+2 ; j <= B.m_je-2 ; j++ )
			  //#pragma omp simd
			  //     for( int i=B.m_ib+2 ; i <= B.m_ie-2 ; i++ )
			  //      {
			  //UpV.print(true);
			  //BV(1,i,j,k) = 0.0;
			  //printf("HERE %d %d %lf %lf \n",i,j,UpV(1,i,j,k-1),BV(1,i,j,k));
			  float_sw4 uz, vz, wz;	 
			  uz = vz = wz = 0.0;
			  //printf("UZ %lf %lf %lf\n",uz,vz,wz);
			  if( upper )
			    {
			      for( int m=0 ; m <= 4 ; m++ )
				{
				  uz += lm_sbop[m]*UpV(1,i,j,k+m-1);
				  vz += lm_sbop[m]*UpV(2,i,j,k+m-1);
				  wz += lm_sbop[m]*UpV(3,i,j,k+m-1);
				}
			    }
			  else
			    {
			      for( int m=0 ; m <= 4 ; m++ )
				{
				  uz -= lm_sbop[m]*UpV(1,i,j,k+1-m);
				  vz -= lm_sbop[m]*UpV(2,i,j,k+1-m);
				  wz -= lm_sbop[m]*UpV(3,i,j,k+1-m);
				}
			    }
			  
         BV(1,i,j,k) = ih*mMuV(i,j,k)*(
	 			       str_x(i)*( a2*(UpV(3,i+2,j,k)-UpV(3,i-2,j,k))+
	 					  a1*(UpV(3,i+1,j,k)-UpV(3,i-1,j,k))) + (uz)  );
         BV(2,i,j,k) = ih*mMuV(i,j,k)*(
            str_y(j)*( a2*(UpV(3,i,j+2,k)-UpV(3,i,j-2,k))+
                       a1*(UpV(3,i,j+1,k)-UpV(3,i,j-1,k))) +
            (vz)  );
         BV(3,i,j,k) = ih*((2*mMuV(i,j,k)+mLambdaV(i,j,k))*(wz) + mLambdaV(i,j,k)*(
                             str_x(i)*( a2*(UpV(1,i+2,j,k)-UpV(1,i-2,j,k))+a1*(UpV(1,i+1,j,k)-UpV(1,i-1,j,k))) +
                             str_y(j)*( a2*(UpV(2,i,j+2,k)-UpV(2,i,j-2,k))+a1*(UpV(2,i,j+1,k)-UpV(2,i,j-1,k))) ) );
         
			}); //SYNC_STREAM;
   //std::cout<<"And we are DONE WITH void EW::compute_icstresses\n";
#undef str_x
#undef str_y
}

//-----------------------------------------------------------------------
void EW::add_ve_stresses( Sarray& a_Up, Sarray& B, int g, int kic, int a_mech, float_sw4* a_str_x, float_sw4* a_str_y)
{
SW4_MARK_FUNCTION;
   const float_sw4 a1=2.0/3, a2=-1.0/12;
   bool upper = (kic == 1);
   int k=kic;
   float_sw4 ih = 1/mGridSize[g];

   int ifirst = a_Up.m_ib;
   int jfirst = a_Up.m_jb;
#define str_x(i) a_str_x[(i-ifirst)]   
#define str_y(j) a_str_y[(j-jfirst)]   


   SView &a_UpV=a_Up.getview();
   SView &BV = B.getview();
   SView &mMuVEV = mMuVE[g][a_mech].getview();
   SView &mLambdaVEV = mLambdaVE[g][a_mech].getview();
   float_sw4 *m_sbop_no_gpV = m_sbop_no_gp;

   RAJA::RangeSegment j_range(B.m_jb+2,B.m_je-2+1);
   RAJA::RangeSegment i_range(B.m_ib+2,B.m_ie-2+1);
// NEW July 21: use new operators WITHOUT ghost points (m_sbop -> m_sbop_no_gp)
// #pragma omp parallel for
//    for( int j=B.m_jb+2 ; j <= B.m_je-2 ; j++ )
// #pragma omp simd
//       for( int i=B.m_ib+2 ; i <= B.m_ie-2 ; i++ )
RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(
				    RAJA::make_tuple(j_range,i_range),
				    [=]RAJA_DEVICE (int j,int i) 
      {
	 float_sw4 uz=0, vz=0, wz=0;
	 if( upper )
	 {
	    for( int m=0 ; m <= 5 ; m++ )
	    {
	       uz += m_sbop_no_gpV[m]*a_UpV(1,i,j,k+m-1);
	       vz += m_sbop_no_gpV[m]*a_UpV(2,i,j,k+m-1);
	       wz += m_sbop_no_gpV[m]*a_UpV(3,i,j,k+m-1);
	    }
	 }
	 else
	 {
	    for( int m=0 ; m <= 5 ; m++ )
	    {
	       uz -= m_sbop_no_gpV[m]*a_UpV(1,i,j,k+1-m);
	       vz -= m_sbop_no_gpV[m]*a_UpV(2,i,j,k+1-m);
	       wz -= m_sbop_no_gpV[m]*a_UpV(3,i,j,k+1-m);
	    }
	 }
// subtract the visco-elastic contribution from mechanism 'a_mech'
         BV(1,i,j,k) = BV(1,i,j,k) - ih*mMuVEV(i,j,k)*(
            str_x(i)*( a2*(a_UpV(3,i+2,j,k)-a_UpV(3,i-2,j,k))+
                       a1*(a_UpV(3,i+1,j,k)-a_UpV(3,i-1,j,k))) + (uz)  );
         BV(2,i,j,k) = BV(2,i,j,k) - ih*mMuVEV(i,j,k)*(
            str_y(j)*( a2*(a_UpV(3,i,j+2,k)-a_UpV(3,i,j-2,k))+
                       a1*(a_UpV(3,i,j+1,k)-a_UpV(3,i,j-1,k))) + (vz)  );
         BV(3,i,j,k) = BV(3,i,j,k) - ih*((2*mMuVEV(i,j,k)+mLambdaVEV(i,j,k))*(wz) + mLambdaVEV(i,j,k)*(
                                          str_x(i)*( a2*(a_UpV(1,i+2,j,k)-a_UpV(1,i-2,j,k))+a1*(a_UpV(1,i+1,j,k)-a_UpV(1,i-1,j,k))) +
                                          str_y(j)*( a2*(a_UpV(2,i,j+2,k)-a_UpV(2,i,j-2,k))+a1*(a_UpV(2,i,j+1,k)-a_UpV(2,i,j-1,k))) ) );
         
      }); //SYNC_STREAM;
#undef str_x
#undef str_y
}

//------------------------------------------------------------------------------
void EW::cartesian_bc_forcing(float_sw4 t, vector<float_sw4 **> & a_BCForcing,
			      vector<Source*>& a_sources )
// assign the boundary forcing arrays a_BCForcing[g][side]
{
  SW4_MARK_FUNCTION;
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz;
  float_sw4 *mu_ptr, *la_ptr, h, zmin;
  boundaryConditionType *bcType_ptr;
  float_sw4 *bforce_side0_ptr, *bforce_side1_ptr, *bforce_side2_ptr, *bforce_side3_ptr, *bforce_side4_ptr, *bforce_side5_ptr;
  int *wind_ptr;
  float_sw4 om=0, ph=0, cv=0, omm;
    
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    SW4_MARK_BEGIN("cart_bc_forcing_iniial");
    mu_ptr    = mMu[g].c_ptr();
    la_ptr    = mLambda[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    nx = m_global_nx[g];
    ny = m_global_ny[g];
    nz = m_global_nz[g];
    
    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    bcType_ptr = m_bcType[g]; // pointer to the local bc array
    zmin = m_zmin[g];
    int curvilinear = topographyExists() && g == mNumberOfGrids-1;
    
    wind_ptr = m_BndryWindow[g];
    
// pointers to the six sides of the cube
    bforce_side0_ptr = a_BCForcing[g][0]; // low-i bndry forcing array pointer
    bforce_side1_ptr = a_BCForcing[g][1]; // high-i bndry forcing array pointer
    bforce_side2_ptr = a_BCForcing[g][2]; // low-j bndry forcing array pointer
    bforce_side3_ptr = a_BCForcing[g][3]; // high-j bndry forcing array pointer
    bforce_side4_ptr = a_BCForcing[g][4]; // low-k bndry forcing array pointer
    bforce_side5_ptr = a_BCForcing[g][5]; // high-k bndry forcing array pointer
    //if (!m_croutines) std::cerr<<"HERE IS THE PROBLME\n";
    SW4_MARK_END("cart_bc_forcing_iniial");
    if (m_twilight_forcing)
    {
       float_sw4 phc[21]; // move these angles to the EW class
       om = m_twilight_forcing->m_omega;
       ph = m_twilight_forcing->m_phase;
       cv = m_twilight_forcing->m_c;
       omm = m_twilight_forcing->m_momega;

      // need to store all the phase angle constants somewhere
      for (int i=0; i<21; i++)
         phc[i] = i*10*M_PI/180;

// the following code can probably be improved by introducing a loop over all sides,
// but bStressFree is only implemented for side=4 and 5, so there must be some special cases
      int k = 1;
      if (m_bcType[g][0] == bDirichlet || m_bcType[g][0] == bSuperGrid )
      {
	SW4_MARK_BEGIN("LOOP1");
         if( !curvilinear )
	 {
	    if( m_croutines )
	       twdirbdry_ci( &wind_ptr[0], h, t, om, cv, ph, bforce_side0_ptr, m_zmin[g] );
	    else
	       twdirbdry( &wind_ptr[0], &h, &t, &om, &cv, &ph, bforce_side0_ptr, &m_zmin[g] );
	 }
	 else
	 {
	    if( m_croutines )
	       twdirbdryc_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
			      &wind_ptr[0], t, om, cv, ph, bforce_side0_ptr,
			      mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	    else
	       twdirbdryc( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                                             &wind_ptr[0], &t, &om, &cv, &ph, bforce_side0_ptr,
					     mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	 }
	 SW4_MARK_END("LOOP1");
      }

      if (m_bcType[g][1] == bDirichlet || m_bcType[g][1] == bSuperGrid )
      {
         if( !curvilinear )
	 {
	    if( m_croutines )
	       twdirbdry_ci( &wind_ptr[6], h, t, om, cv, ph, bforce_side1_ptr, m_zmin[g] );
	    else
	       twdirbdry( &wind_ptr[6], &h, &t, &om, &cv, &ph, bforce_side1_ptr, &m_zmin[g] );
	 }
	 else
	 {
	    if( m_croutines )
	       twdirbdryc_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
		       &wind_ptr[6], t, om, cv, ph, bforce_side1_ptr,
		       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	    else
	       twdirbdryc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			  &wind_ptr[6], &t, &om, &cv, &ph, bforce_side1_ptr,
			  mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	 }
      }
      	SW4_MARK_BEGIN("LOOP2");
      if (m_bcType[g][2] == bDirichlet || m_bcType[g][2] == bSuperGrid)
      {
	 if( !curvilinear )
	 {
	    if( m_croutines )
	       twdirbdry_ci( &wind_ptr[6*2], h, t, om, cv, ph, bforce_side2_ptr, m_zmin[g] );
	    else
	       twdirbdry( &wind_ptr[6*2], &h, &t, &om, &cv, &ph, bforce_side2_ptr, &m_zmin[g] );
	 }
         else
	 {
	    if( m_croutines )
	       twdirbdryc_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
		       &wind_ptr[6*2], t, om, cv, ph, bforce_side2_ptr,
		       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	    else
	       twdirbdryc( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			   &wind_ptr[6*2], &t, &om, &cv, &ph, bforce_side2_ptr,
			   mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	 }
      }

      if (m_bcType[g][3] == bDirichlet || m_bcType[g][3] == bSuperGrid)
      {
         if( !curvilinear )
	 {
	    if( m_croutines )
	       twdirbdry_ci( &wind_ptr[6*3], h, t, om, cv, ph, bforce_side3_ptr, m_zmin[g] );
	    else
	       twdirbdry( &wind_ptr[6*3], &h, &t, &om, &cv, &ph, bforce_side3_ptr, &m_zmin[g] );
	 }
         else
	 {
	    if( m_croutines )
	       twdirbdryc_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
		       &wind_ptr[6*3], t, om, cv, ph, bforce_side3_ptr,
		       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	    else
	       twdirbdryc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			  &wind_ptr[6*3], &t, &om, &cv, &ph, bforce_side3_ptr,
			  mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	 }
      }
      SW4_MARK_END("LOOP2");
      SW4_MARK_BEGIN("LOOP3");
      if (m_bcType[g][4] == bDirichlet || m_bcType[g][4] == bSuperGrid)
      {
	 if( !curvilinear )
	 {
	    if( m_croutines )
	       twdirbdry_ci( &wind_ptr[6*4], h, t, om, cv, ph, bforce_side4_ptr, m_zmin[g] );
	    else
	       twdirbdry( &wind_ptr[6*4], &h, &t, &om, &cv, &ph, bforce_side4_ptr, &m_zmin[g] );
	 }
         else
	 {
	    if( m_croutines )
	       twdirbdryc_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
		       &wind_ptr[6*4], t, om, cv, ph, bforce_side4_ptr,
		       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	    else
	       twdirbdryc( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			   &wind_ptr[6*4], &t, &om, &cv, &ph, bforce_side4_ptr,
			   mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	 }
      }
      else if (m_bcType[g][4] == bStressFree)
      {
	 k = 1;
         if( m_anisotropic )
         {
// curvilinear anisotropic case is not yet implemented
            CHECK_INPUT (!curvilinear, "cartesian_bc_forcing> bStressFree not implemented for anisotropic materials and curvilinear grids" <<endl);

	    if( m_croutines )
	       tw_aniso_free_surf_z_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, k, t, om,
					cv, ph, omm, phc, bforce_side4_ptr, h, m_zmin[g] );            
	    else
	       tw_aniso_free_surf_z( ifirst, ilast, jfirst, jlast, kfirst, klast, k, t, om,
				     cv, ph, omm, phc, bforce_side4_ptr, h, m_zmin[g] );            
         }
         else
         { //isotropic stuff
            
            if( usingSupergrid() && !curvilinear )
            {
               float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
               float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
	       if( m_croutines )
		  twfrsurfzsgstr_ci( ifirst, ilast, jfirst, jlast, kfirst, 
				  klast, h, k, t, om, cv, ph, omstrx, omstry,
				  bforce_side4_ptr, mu_ptr, la_ptr, m_zmin[g] );
	       else
		  twfrsurfzsgstr( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				  &klast, &h, &k, &t, &om, &cv, &ph, &omstrx, &omstry,
				  bforce_side4_ptr, mu_ptr, la_ptr, &m_zmin[g] );
               if( m_use_attenuation )
               {
                  float_sw4* mua_ptr    = mMuVE[g][0].c_ptr();
                  float_sw4* laa_ptr    = mLambdaVE[g][0].c_ptr();
		  if( m_croutines )
		     twfrsurfzsgstratt_ci( ifirst, ilast, jfirst, jlast, kfirst, 
					   klast, h, k, t, om, cv, ph, omstrx, omstry,
					   bforce_side4_ptr, mua_ptr, laa_ptr, m_zmin[g] );
		  else
		     twfrsurfzsgstratt( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					&klast, &h, &k, &t, &om, &cv, &ph, &omstrx, &omstry,
					bforce_side4_ptr, mua_ptr, laa_ptr, &m_zmin[g] );
               }
            }
            else if( !usingSupergrid() && curvilinear )
            {
               // Stress tensor on boundary
               Sarray tau(6,ifirst,ilast,jfirst,jlast,1,1);
               // Get twilight stress tensor, tau.
	       if( m_croutines )
		  twstensor_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
				k, t, om, cv, ph,
				mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(), mu_ptr, la_ptr );
	       else
		  twstensor( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			     &k, &t, &om, &cv, &ph,
			     mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(), mu_ptr, la_ptr );
               // Compute boundary forcing for given stress tensor, tau.

	       if( m_croutines )
		  getsurfforcing_ci( ifirst, ilast, jfirst, jlast, kfirst,
				     klast, k, mMetric.c_ptr(), mJ.c_ptr(),
				     tau.c_ptr(), bforce_side4_ptr );
	       else
		  getsurfforcing( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
				  &klast, &k, mMetric.c_ptr(), mJ.c_ptr(),
				  tau.c_ptr(), bforce_side4_ptr );

               if( m_use_attenuation )
               {
                  float_sw4* mua_ptr    = mMuVE[g][0].c_ptr();
                  float_sw4* laa_ptr    = mLambdaVE[g][0].c_ptr();
		  if( m_croutines )
		     twstensoratt_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
				      k, t, om, cv, ph,
				      mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(), mua_ptr, laa_ptr );
		  else
		     twstensoratt( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
				   &k, &t, &om, &cv, &ph,
				   mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(), mua_ptr, laa_ptr );
		  if( m_croutines )
		     subsurfforcing_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, k,
					mMetric.c_ptr(), mJ.c_ptr(), tau.c_ptr(), bforce_side4_ptr );
		  else
		     subsurfforcing( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
				     &klast, &k, mMetric.c_ptr(), mJ.c_ptr(),
				     tau.c_ptr(), bforce_side4_ptr );
               }
            }
            else if( !usingSupergrid() && !curvilinear )
            {
	       if( m_croutines )
		  twfrsurfz_ci( ifirst, ilast, jfirst, jlast, kfirst, 
				klast, h, k, t, om, cv, ph,
				bforce_side4_ptr, mu_ptr, la_ptr, m_zmin[g] );
	       else
		  twfrsurfz( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
			     &klast, &h, &k, &t, &om, &cv, &ph,
			     bforce_side4_ptr, mu_ptr, la_ptr, &m_zmin[g] );
               if( m_use_attenuation )
               {
                  float_sw4* mua_ptr    = mMuVE[g][0].c_ptr();
                  float_sw4* laa_ptr    = mLambdaVE[g][0].c_ptr();
		  if( m_croutines )
		     twfrsurfzatt_ci( ifirst, ilast, jfirst, jlast, kfirst, 
				      klast, h, k, t, om, cv, ph,
				      bforce_side4_ptr, mua_ptr, laa_ptr, m_zmin[g] );
		  else
		     twfrsurfzatt( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				   &klast, &h, &k, &t, &om, &cv, &ph,
				   bforce_side4_ptr, mua_ptr, laa_ptr, &m_zmin[g] );
	       
               }
            }
            else if( usingSupergrid() && curvilinear )
            {
               float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
               float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();

               // Stress tensor on boundary
               Sarray tau(6,ifirst,ilast,jfirst,jlast,1,1);
               // Get twilight stress tensor, tau.
	       if( m_croutines )
		  twstensorsg_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
				  k, t, om, cv, ph,
				  mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(),
				  mu_ptr, la_ptr, omstrx, omstry );
	       else
		  twstensorsg( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			       &k, &t, &om, &cv, &ph,
			       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(),
			       mu_ptr, la_ptr, &omstrx, &omstry );
               // Compute boundary forcing for given stress tensor, tau.
	       if( m_croutines )
		  getsurfforcingsg_ci( ifirst, ilast, jfirst, jlast, kfirst,
				       klast, k, mMetric.c_ptr(), mJ.c_ptr(),
				       tau.c_ptr(), m_sg_str_x[g], m_sg_str_y[g], bforce_side4_ptr );
	       else
		  getsurfforcingsg( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
				    &klast, &k, mMetric.c_ptr(), mJ.c_ptr(),
				    tau.c_ptr(), m_sg_str_x[g], m_sg_str_y[g], bforce_side4_ptr );

               if( m_use_attenuation )
               {
                  float_sw4* mua_ptr    = mMuVE[g][0].c_ptr();
                  float_sw4* laa_ptr    = mLambdaVE[g][0].c_ptr();
		  if( m_croutines )
		  {
		     twstensorsgatt_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
					k, t, om, cv, ph,
					mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(),
					mua_ptr, laa_ptr, omstrx, omstry );
		     subsurfforcingsg_ci( ifirst, ilast, jfirst, jlast, kfirst,
					  klast, k, mMetric.c_ptr(), mJ.c_ptr(),
					  tau.c_ptr(), m_sg_str_x[g], m_sg_str_y[g],
					  bforce_side4_ptr );
		  }
		  else
		  {
		     twstensorsgatt( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
				     &k, &t, &om, &cv, &ph,
				     mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(),
				     mua_ptr, laa_ptr, &omstrx, &omstry );
		     subsurfforcingsg( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
				       &klast, &k, mMetric.c_ptr(), mJ.c_ptr(),
				       tau.c_ptr(), m_sg_str_x[g], m_sg_str_y[g],
				       bforce_side4_ptr );
		  }
               }
            } // end supergrid && curvilinear
         
         } // end isotropic case
         
      } // end side==4 is bStressFree
      SW4_MARK_END("LOOP3");
  
      if (m_bcType[g][5] == bDirichlet || m_bcType[g][5] == bSuperGrid)
      {
	SW4_MARK_BEGIN("LOOP4");
	 if( !curvilinear )
	 {
	    if( m_croutines )
	       twdirbdry_ci( &wind_ptr[6*5], h, t, om, cv, ph, bforce_side5_ptr, m_zmin[g] );
	    else
	       twdirbdry( &wind_ptr[6*5], &h, &t, &om, &cv, &ph, bforce_side5_ptr, &m_zmin[g] );
	 }
	 else
	 {
	    if( m_croutines )
	       twdirbdryc_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
		       &wind_ptr[6*5], t, om, cv, ph, bforce_side5_ptr,
		       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	    else
	       twdirbdryc( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			   &wind_ptr[6*5], &t, &om, &cv, &ph, bforce_side5_ptr,
			   mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
	 }
	 SW4_MARK_END("LOOP4");
      }
      else if (m_bcType[g][5] == bStressFree)
      {
	SW4_MARK_BEGIN("LOOP5");
	 k = nz;
         if( m_anisotropic )
         {
// curvilinear anisotropic case is not yet implemented
            CHECK_INPUT (!curvilinear, "cartesian_bc_forcing> bStressFree not implemented for anisotropic materials and curvilinear grids" <<endl);

	    if( m_croutines )
	       tw_aniso_free_surf_z_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, k, t, om,
					cv, ph, omm, phc, bforce_side5_ptr, h, m_zmin[g] );
	    else
	       tw_aniso_free_surf_z( ifirst, ilast, jfirst, jlast, kfirst, klast, k, t, om,
				     cv, ph, omm, phc, bforce_side5_ptr, h, m_zmin[g] );            
         }
         else
         { //isotropic stuff

            if( usingSupergrid() )
            {
               float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
               float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
               twfrsurfzsgstr( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
			       &klast, &h, &k, &t, &om, &cv, &ph, &omstrx, &omstry,
			       bforce_side5_ptr, mu_ptr, la_ptr, &m_zmin[g] );
               if( m_use_attenuation )
               {
                  float_sw4* mua_ptr    = mMuVE[g][0].c_ptr();
                  float_sw4* laa_ptr    = mLambdaVE[g][0].c_ptr();
                  twfrsurfzsgstratt( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				     &klast, &h, &k, &t, &om, &cv, &ph, &omstrx, &omstry,
				     bforce_side5_ptr, mua_ptr, laa_ptr, &m_zmin[g] );
	       
               }
            }
            else
            {
               twfrsurfz( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
			  &klast, &h, &k, &t, &om, &cv, &ph,
			  bforce_side5_ptr, mu_ptr, la_ptr, &m_zmin[g] );
               if( m_use_attenuation )
               {
                  float_sw4* mua_ptr    = mMuVE[g][0].c_ptr();
                  float_sw4* laa_ptr    = mLambdaVE[g][0].c_ptr();
                  twfrsurfzatt( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				&klast, &h, &k, &t, &om, &cv, &ph,
				bforce_side5_ptr, mua_ptr, laa_ptr, &m_zmin[g] );
               }
            } // end ! supergrid
            
         } // end isotropic case

	 SW4_MARK_END("LOOP5");

      } // end bStressFree on side 5

    }
    else if (m_rayleigh_wave_test)
    {

      int q;
      float_sw4 lambda, mu, rho, cr, omega, alpha;
      
      lambda = m_rayleigh_wave_test->m_lambda;
      mu = m_rayleigh_wave_test->m_mu;
      rho = m_rayleigh_wave_test->m_rho;
      cr = m_rayleigh_wave_test->m_cr;
      omega = m_rayleigh_wave_test->m_omega;
      alpha = m_rayleigh_wave_test->m_alpha;

// homogneous free surface bc (low-z)
      for (q=0; q<3*m_NumberOfBCPoints[g][4]; q++)
	bforce_side4_ptr[q] = 0.;

// assign exact solution on bottom (high-z)

      if (m_bcType[g][5] == bDirichlet)
      {
	SW4_MARK_BEGIN("raydirbdry");
	raydirbdry( bforce_side5_ptr, &wind_ptr[6*5], &t, &lambda, &mu, &rho, &cr, 
		    &omega, &alpha, &h, &zmin );
	SW4_MARK_END("raydirbdry");
      }

     //  subroutine RAYDIRBDRY( bforce, wind, t, lambda, mu, rho, cr, 
     // +     omega, alpha, h, zmin )
    }
    else if( m_point_source_test )
    {
       for( int side=0 ; side < 6 ; side++ )
	  if( m_bcType[g][side] == bDirichlet )
	     get_exact_point_source( a_BCForcing[g][side], t, g, *a_sources[0], &wind_ptr[6*side] );
	  else
	     for (int q=0; q<3*m_NumberOfBCPoints[g][side]; q++)
		a_BCForcing[g][side][q] = 0.;

    }
    else
    {
// no boundary forcing
// we can do the same loop for all types of bc. For bParallel boundaries, numberOfBCPoints=0
    SW4_MARK_BEGIN("LOOP6");
    
      //for (q=0; q<3*m_NumberOfBCPoints[g][0]; q++)
	RAJA::forall<DEFAULT_LOOP1> (RAJA::RangeSegment(0,3*m_NumberOfBCPoints[g][0]),[=] RAJA_DEVICE(int q){
	    bforce_side0_ptr[q] = 0.;});
	//for (q=0; q<3*m_NumberOfBCPoints[g][1]; q++)
	RAJA::forall<DEFAULT_LOOP1> (RAJA::RangeSegment(0,3*m_NumberOfBCPoints[g][1]),[=] RAJA_DEVICE(int q){
	    bforce_side1_ptr[q] = 0.;});
	//for (q=0; q<3*m_NumberOfBCPoints[g][2]; q++)
	RAJA::forall<DEFAULT_LOOP1> (RAJA::RangeSegment(0,3*m_NumberOfBCPoints[g][2]),[=] RAJA_DEVICE(int q){
	    bforce_side2_ptr[q] = 0.;});
	//for (q=0; q<3*m_NumberOfBCPoints[g][3]; q++)
	RAJA::forall<DEFAULT_LOOP1> (RAJA::RangeSegment(0,3*m_NumberOfBCPoints[g][3]),[=] RAJA_DEVICE(int q){
	    bforce_side3_ptr[q] = 0.;});
	//for (q=0; q<3*m_NumberOfBCPoints[g][4]; q++)
	RAJA::forall<DEFAULT_LOOP1> (RAJA::RangeSegment(0,3*m_NumberOfBCPoints[g][4]),[=] RAJA_DEVICE(int q){
	    bforce_side4_ptr[q] = 0.;});
	//for (q=0; q<3*m_NumberOfBCPoints[g][5]; q++)
	RAJA::forall<DEFAULT_LOOP1> (RAJA::RangeSegment(0,3*m_NumberOfBCPoints[g][5]),[=] RAJA_DEVICE(int q){
	    bforce_side5_ptr[q] = 0.;});
	SYNC_STREAM;
      SW4_MARK_END("LOOP6");
    }
  }
}

// //---------------------------------------------------------------------------
// void eval_curvilinear_bc_stress(Sarray & a_u, double ** bcForcing, Sarray & a_x, Sarray & a_y, Sarray & a_z,
// 				Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J)
// {
// // 4D macros swap the last and first indices to compensate for different conventions between how 
// // the Sarrays were allocated and how this routine was originally written
// #define u(i,j,k,c) u_[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)]
// #define q(i,j,k,c) q_[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)]
// #define r(i,j,k,c) r_[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)]
// #define s(i,j,k,c) s_[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)]
// // 3D array macros are special cases of the 4D macros with c=1 and nc=1
// #define x(i,j,k) x_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
// #define y(i,j,k) y_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
// #define z(i,j,k) z_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
// #define J(i,j,k) J_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
// #define mu(i,j,k) mu_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
// #define lam(i,j,k) lam_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
// // not necessary to store lambda + 2*mu in separate array
// #define lam2mu(i,j,k) (lam(i,j,k) + 2.*mu(i,j,k))
  
// // extract pointers for the macros
// // 4D arrays
//   double * u_=a_u.c_ptr();
//   double * q_=a_q.c_ptr();
//   double * r_=a_r.c_ptr();
//   double * s_=a_s.c_ptr();
// // 3D arrays
//   double * x_=a_x.c_ptr();
//   double * y_=a_y.c_ptr();
//   double * z_=a_z.c_ptr();
//   double * mu_=a_mu.c_ptr();
//   double * lam_=a_lam.c_ptr();
//   double * J_=a_J.c_ptr();
  
// // all 3D/4D Sarrays must have the same number of grid points and the same starting/ending indices
//   int m_nc = a_q.m_nc;
//   int m_ni = a_q.m_ni;
//   int m_nj = a_q.m_nj;
//   int m_nk = a_q.m_nk;
// // to mimic the original coding:
// // setting starting indices to one 
// // setting ending indices to equal the number of points in each dimension
//   int m_ib = 1;
//   int m_jb = 1;
//   int m_kb = 1;
//   int Nx = a_q.m_ni;
//   int Ny = a_q.m_nj;
//   int Nz = a_q.m_nk;

//   int i, j, k, q, side, ind;
   
// // only implemented for the low-k boundary
//   side=4;
  
//   {
    
//     ind = 0;
// #define E1(i,j,t1,t2,i2,t3,i3,t4) (t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2))
// #define E12(i,j,t1,t2,i2,t3,i3,t4) (0.5*(t1(i,j,1)*t2(i,j,1,i2)*t3(i,j,1,i3)*t4(i,j,1)\
// 				 +t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2)))
// #define E32(i,j,t1,t2,i2,t3,i3,t4) (0.5*(t1(i,j,3)*t2(i,j,3,i2)*t3(i,j,3,i3)*t4(i,j,3)\
// 				 +t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2)))
//     double Dr0[3],Dq0[3],Dsp[3],ugp[3],b[3],A[9],x[3],Ax[3];
//     int info,N1,N2,N3 ,ipv[3]; 
//     N1=3;N2=1;N3=3;
//     k=2; ///NOTE!!! 
//     for(int j=2; j<=Ny-1; j++)
//       for(int i=2; i<=Nx-1; i++)
//       {
// 	for(int c=1; c<=3;c++){
// 	  Dr0[c-1]=0.5*(u(i,j+1,2,c)-u(i,j-1,2,c));
// 	  Dq0[c-1]=0.5*(u(i+1,j,2,c)-u(i-1,j,2,c));
// 	  Dsp[c-1]=u(i,j,3,c)-u(i,j,2,c);
// 	}
// 	b[0]=-(E1(i,j,J,s,1,q,1,lam2mu)*Dq0[0]
// 	       +E1(i,j,J,s,1,r,1,lam2mu)*Dr0[0]	
// 	       +E1(i,j,J,s,1,q,2,lam)*Dq0[1]
// 	       +E1(i,j,J,s,1,r,2,lam)*Dr0[1]	 
// 	       +E1(i,j,J,s,1,q,3,lam)*Dq0[2]
// 	       +E1(i,j,J,s,1,r,3,lam)*Dr0[2]	 
// 	       +E1(i,j,J,s,2,q,1,mu)*Dq0[1]
// 	       +E1(i,j,J,s,2,r,1,mu)*Dr0[1]	 
// 	       +E1(i,j,J,s,2,q,2,mu)*Dq0[0]
// 	       +E1(i,j,J,s,2,r,2,mu)*Dr0[0]	 
// 	       +E1(i,j,J,s,3,q,1,mu)*Dq0[2]
// 	       +E1(i,j,J,s,3,r,1,mu)*Dr0[2]	 
// 	       +E1(i,j,J,s,3,q,3,mu)*Dq0[0]
// 	       +E1(i,j,J,s,3,r,3,mu)*Dr0[0]	 
// 	       +0.5*(
// 		 E32(i,j,J,s,1,s,1,lam2mu)*Dsp[0]
// 		 +E32(i,j,J,s,1,s,2,lam)*Dsp[1]
// 		 +E32(i,j,J,s,1,s,3,lam)*Dsp[2]
// 		 +E32(i,j,J,s,2,s,1,mu)*Dsp[1]
// 		 +E32(i,j,J,s,2,s,2,mu)*Dsp[0]
// 		 +E32(i,j,J,s,3,s,1,mu)*Dsp[2]
// 		 +E32(i,j,J,s,3,s,3,mu)*Dsp[0]));
// 	b[1]=-(E1(i,j,J,s,1,q,1,mu)*Dq0[1]
// 	       +E1(i,j,J,s,1,r,1,mu)*Dr0[1]	
// 	       +E1(i,j,J,s,1,q,2,mu)*Dq0[0]
// 	       +E1(i,j,J,s,1,r,2,mu)*Dr0[0]	 
// 	       +E1(i,j,J,s,2,q,2,lam2mu)*Dq0[1]
// 	       +E1(i,j,J,s,2,r,2,lam2mu)*Dr0[1]
// 	       +E1(i,j,J,s,2,r,1,lam)*Dr0[0]	 
// 	       +E1(i,j,J,s,2,q,1,lam)*Dq0[0]
// 	       +E1(i,j,J,s,2,r,3,lam)*Dr0[2]	 
// 	       +E1(i,j,J,s,2,q,3,lam)*Dq0[2]
// 	       +E1(i,j,J,s,3,q,2,mu)*Dq0[2]
// 	       +E1(i,j,J,s,3,r,2,mu)*Dr0[2]	 
// 	       +E1(i,j,J,s,3,q,3,mu)*Dq0[1]
// 	       +E1(i,j,J,s,3,r,3,mu)*Dr0[1]	 
// 	       +0.5*(E32(i,j,J,s,2,s,2,lam2mu)*Dsp[1]
// 		     +E32(i,j,J,s,2,s,1,lam)*Dsp[0]
// 		     +E32(i,j,J,s,2,s,3,lam)*Dsp[2]
// 		     +E32(i,j,J,s,1,s,1,mu)*Dsp[1]
// 		     +E32(i,j,J,s,1,s,2,mu)*Dsp[0]
// 		     +E32(i,j,J,s,3,s,2,mu)*Dsp[2]
// 		     +E32(i,j,J,s,3,s,3,mu)*Dsp[1]));
	    
// 	b[2]=-(E1(i,j,J,s,1,q,1,mu)*Dq0[2]
// 	       +E1(i,j,J,s,1,r,1,mu)*Dr0[2]	
// 	       +E1(i,j,J,s,1,q,3,mu)*Dq0[0]
// 	       +E1(i,j,J,s,1,r,3,mu)*Dr0[0]	 
// 	       +E1(i,j,J,s,3,q,3,lam2mu)*Dq0[2]
// 	       +E1(i,j,J,s,3,r,3,lam2mu)*Dr0[2]
// 	       +E1(i,j,J,s,3,r,1,lam)*Dr0[0]	 
// 	       +E1(i,j,J,s,3,q,1,lam)*Dq0[0]
// 	       +E1(i,j,J,s,3,r,2,lam)*Dr0[1]	 
// 	       +E1(i,j,J,s,3,q,2,lam)*Dq0[1]
// 	       +E1(i,j,J,s,2,q,2,mu)*Dq0[2]
// 	       +E1(i,j,J,s,2,r,2,mu)*Dr0[2]	 
// 	       +E1(i,j,J,s,2,q,3,mu)*Dq0[1]
// 	       +E1(i,j,J,s,2,r,3,mu)*Dr0[1]	 
// 	       +0.5*(E32(i,j,J,s,3,s,3,lam2mu)*Dsp[2]
// 		     +E32(i,j,J,s,3,s,1,lam)*Dsp[0]
// 		     +E32(i,j,J,s,3,s,2,lam)*Dsp[1]
// 		     +E32(i,j,J,s,1,s,1,mu)*Dsp[2]
// 		     +E32(i,j,J,s,1,s,3,mu)*Dsp[0]
// 		     +E32(i,j,J,s,2,s,2,mu)*Dsp[2]
// 		     +E32(i,j,J,s,2,s,3,mu)*Dsp[1]));


// 	A[0]=(
// 	  E12(i,j,J,s,1,s,1,lam2mu)
// 	  +E12(i,j,J,s,2,s,2,mu)
// 	  +E12(i,j,J,s,3,s,3,mu));
// 	A[3]=(E12(i,j,J,s,1,s,2,lam)
// 	      +E12(i,j,J,s,2,s,1,mu));
// 	A[6]=(E12(i,j,J,s,1,s,3,lam)
// 	      +E12(i,j,J,s,3,s,1,mu));
// // u, v, w in v eq.
// 	A[1]=(E12(i,j,J,s,2,s,1,lam)
// 	      +E12(i,j,J,s,1,s,2,mu));
// 	A[4]=(E12(i,j,J,s,3,s,3,mu)
// 	      +E12(i,j,J,s,1,s,1,mu)
// 	      +E12(i,j,J,s,2,s,2,lam2mu));
// 	A[7]=(E12(i,j,J,s,2,s,3,lam)
// 	      +E12(i,j,J,s,3,s,2,mu));
// 	// u, v, w in w eq.
// 	A[2]=(E12(i,j,J,s,3,s,1,lam)
// 	      +E12(i,j,J,s,1,s,3,mu));
// 	A[5]=(E12(i,j,J,s,2,s,3,mu)
// 	      +E12(i,j,J,s,3,s,2,lam));
// 	A[8]=(E12(i,j,J,s,3,s,3,lam2mu)
// 	      +E12(i,j,J,s,1,s,1,mu)
// 	      +E12(i,j,J,s,2,s,2,mu));
// 	for(int c=0; c<9; c++)
// 	  A[c]*=0.5;

// 	for (int c=1; c<=3; c++)
// 	  x[c-1] = u(i,j,2,c) - u(i,j,1,c);
	    
// 	Ax[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
// 	Ax[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
// 	Ax[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];

// 	for(int c=1; c<=3;c++)
// 	{
// // this routine is used to accumulate viscoelastic boundary stresses from each mechanism
// 	  bcForcing[side][ind+c-1] += Ax[c-1] - b[c-1];
// 	}
// 	ind += 3;
	    
//       }
    
//   }
  
  
// #undef u
// #undef x
// #undef y
// #undef z
// #undef q
// #undef r
// #undef s
// #undef J
// #undef mu
// #undef lam
// #undef lam2mu
  
// }

//-----------------------------------------------------------------------
void EW::test_sources( vector<GridPointSource*>& a_point_sources,
		       vector<Source*>& a_global_unique_sources, vector<Sarray>& a_F,
		       vector<int>& identsources )
{
  SW4_MARK_FUNCTION;
// Check the source discretization
  int kx[3] = {0,0,0};
  int ky[3] = {0,0,0};
  int kz[3] = {0,0,0};
  float_sw4 moments[3], momexact[3];
  int nsourcesloc = a_point_sources.size();
  int nsources;
  MPI_Allreduce( &nsourcesloc, &nsources, 1, MPI_INT, MPI_SUM, m_cartesian_communicator );

  if( proc_zero() )
  {
     cout << "Source test " << endl;
     cout << "source size = " << a_global_unique_sources.size() << endl;
     cout << "grid point source size = " << nsources << endl;
  }
  for( int c=0; c <= 7 ; c++ )
  {
     kx[0] = c;
     ky[1] = c;
     kz[2] = c;
     testSourceDiscretization( kx, ky, kz, moments, a_point_sources, a_F, identsources );
     a_global_unique_sources[0]->exact_testmoments( kx, ky, kz, momexact );
     if( proc_zero() )
     {
        for( int comp = 0 ; comp < 3 ; comp++ )
           cout << kx[comp] << " " << ky[comp] << " " << kz[comp] << " computed " << moments[comp] <<
              " exact " << momexact[comp] << " difference = " << moments[comp]-momexact[comp] << endl;
     }
  }  
  kx[0] = 1;
  ky[0] = 1;
  kz[0] = 1;
  kx[1] = 2;
  ky[1] = 1;
  kz[1] = 1;
  kx[2] = 1;
  ky[2] = 2;
  kz[2] = 1;
  testSourceDiscretization( kx, ky, kz, moments, a_point_sources, a_F, identsources );
  a_global_unique_sources[0]->exact_testmoments( kx, ky, kz, momexact );
  if( proc_zero() )
  {
     for( int comp = 0 ; comp < 3 ; comp++ )
        cout << kx[comp] << " " << ky[comp] << " " << kz[comp] << " computed " << moments[comp] <<
           " exact " << momexact[comp] << " difference = " << moments[comp]-momexact[comp] << endl;
  }
  kx[0] = 3;
  ky[0] = 2;
  kz[0] = 2;
  kx[1] = 2;
  ky[1] = 3;
  kz[1] = 2;
  kx[2] = 2;
  ky[2] = 2;
  kz[2] = 3;
  testSourceDiscretization( kx, ky, kz, moments, a_point_sources, a_F, identsources );
  a_global_unique_sources[0]->exact_testmoments( kx, ky, kz, momexact );
  if( proc_zero() )
  {
     for( int comp = 0 ; comp < 3 ; comp++ )
        cout << kx[comp] << " " << ky[comp] << " " << kz[comp] << " computed " << moments[comp] <<
           " exact " << momexact[comp] << " difference = " << moments[comp]-momexact[comp] << endl;
  }
  kx[0] = 4;
  ky[0] = 3;
  kz[0] = 3;
  kx[1] = 3;
  ky[1] = 4;
  kz[1] = 3;
  kx[2] = 3;
  ky[2] = 3;
  kz[2] = 4;
  testSourceDiscretization( kx, ky, kz, moments, a_point_sources, a_F, identsources );
  a_global_unique_sources[0]->exact_testmoments( kx, ky, kz, momexact );
  if( proc_zero() )
  {
     for( int comp = 0 ; comp < 3 ; comp++ )
        cout << kx[comp] << " " << ky[comp] << " " << kz[comp] << " computed " << moments[comp] <<
           " exact " << momexact[comp] << " difference = " << moments[comp]-momexact[comp] << endl;
  }
}

//-----------------------------------------------------------------------
void EW::testSourceDiscretization( int kx[3], int ky[3], int kz[3],
				   float_sw4 moments[3],
				   vector<GridPointSource*>& point_sources,
				   vector<Sarray>& F, vector<int>& identsources )
{
   // Evaluate sources at a large time (assume that the time function is=1 at t=infinity)
   // Compute moments, integrals of the source times polynomials of degree (kx,ky,kz).
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 h;

//tmp  
  if (proc_zero())
    printf("Inside testSourceDiscretization\n");
  
  // Impose source
  for( int g=0 ; g<mNumberOfGrids; g++ )
     F[g].set_to_zero();
#pragma omp parallel for
  for( int r=0 ; r < identsources.size()-1 ; r++ )
  {
     int s0=identsources[r];
     int g= point_sources[s0]->m_grid;	
     int i= point_sources[s0]->m_i0;
     int j= point_sources[s0]->m_j0;
     int k= point_sources[s0]->m_k0;
     float_sw4 f1=0, f2=0, f3=0;
     for( int s = identsources[r] ; s < identsources[r+1] ; s++ )
     {
	float_sw4 fxyz[3];
	point_sources[s]->getFxyz_notime(fxyz);
	f1 += fxyz[0];
	f2 += fxyz[1];
	f3 += fxyz[2];
     }
     F[g](1,i,j,k) += f1;
     F[g](2,i,j,k) += f2;
     F[g](3,i,j,k) += f3;
  }
  //  for( int s= 0 ; s < point_sources.size() ; s++ ) 
  //  {
  //     float_sw4 fxyz[3];
  //     point_sources[s]->getFxyz_notime( fxyz );
  //     int g = point_sources[s]->m_grid;
  //     F[g](1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
  //     F[g](2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
  //     F[g](3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
  //  }
  float_sw4 momgrid[3]={0,0,0};
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];  
    h      = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    float_sw4* f_ptr = F[g].c_ptr();
    int wind[6];
    wind[0] = m_iStartInt[g];
    wind[1] = m_iEndInt[g];
    wind[2] = m_jStartInt[g];
    wind[3] = m_jEndInt[g];
    wind[4] = m_kStartInt[g];
    wind[5] = m_kEndInt[g];
    int nz = m_global_nz[g];
    if( m_croutines )
       testsrc_ci( f_ptr, ifirst, ilast, jfirst, jlast, kfirst, klast,
		   nz, wind, m_zmin[g], h, kx, ky, kz, momgrid );
    else
       testsrc( f_ptr, &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
		&nz, wind, &m_zmin[g], &h, kx, ky, kz, momgrid );
  }
  MPI_Allreduce( momgrid, moments, 3, m_mpifloat, MPI_SUM, m_cartesian_communicator );
}

//-----------------------------------------------------------------------
void EW::extractRecordData(TimeSeries::receiverMode mode, int i0, int j0, int k0, int g0, 
			   vector<float_sw4> &uRec, vector<Sarray> &Um2, vector<Sarray> &U)
{
SW4_MARK_FUNCTION;
  if (mode == TimeSeries::Displacement)
  {
    uRec.resize(3);
    uRec[0] = U[g0](1, i0, j0, k0);
    uRec[1] = U[g0](2, i0, j0, k0);
    uRec[2] = U[g0](3, i0, j0, k0);
  }
  else if (mode == TimeSeries::Velocity)
  {
    uRec.resize(3);
    uRec[0] = (U[g0](1, i0, j0, k0) - Um2[g0](1, i0, j0, k0))/(2*mDt);
    uRec[1] = (U[g0](2, i0, j0, k0) - Um2[g0](2, i0, j0, k0))/(2*mDt);
    uRec[2] = (U[g0](3, i0, j0, k0) - Um2[g0](3, i0, j0, k0))/(2*mDt);
  }
  else if(mode == TimeSeries::Div)
  {
    uRec.resize(1);
    if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
    {
//      int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      float_sw4 factor = 1.0/(2*mGridSize[g0]);
      uRec[0] = ((U[g0](1,i0+1, j0, k0) - U[g0](1,i0-1, j0, k0)+
		  U[g0](2,i0, j0+1, k0) - U[g0](2,i0, j0-1, k0)+
		  U[g0](3,i0, j0, k0+1) - U[g0](3,i0, j0, k0-1))*factor);
    }
    else // must be curvilinear
    {
//      int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
       float_sw4 factor = 0.5/sqrt(mJ(i0,j0,k0));
       uRec[0] = ( ( mMetric(1,i0,j0,k0)*(U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))+
		     mMetric(1,i0,j0,k0)*(U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))+
		     mMetric(2,i0,j0,k0)*(U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))+
		     mMetric(3,i0,j0,k0)*(U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))+
		     mMetric(4,i0,j0,k0)*(U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1))  )*factor);
    }
  } // end div
  else if(mode == TimeSeries::Curl)
  {
    uRec.resize(3);
    if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      float_sw4 factor = 1.0/(2*mGridSize[g0]);
      float_sw4 duydx = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0))*factor;
      float_sw4 duzdx = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0))*factor;
      float_sw4 duxdy = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0))*factor;
      float_sw4 duzdy = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0))*factor;
      float_sw4 duxdz = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))*factor;
      float_sw4 duydz = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))*factor;
//       if( m_xycomponent )
//       {
      uRec[0] = ( duzdy-duydz );
      uRec[1] = ( duxdz-duzdx );
      uRec[2] = ( duydx-duxdy );
//       }
//       else
//       {
// 	 float_sw4 uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
// 	 float_sw4 uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
// 	 mRecordedUX.push_back( uew );
// 	 mRecordedUY.push_back( uns );
// 	 mRecordedUZ.push_back( -(duydx-duxdy) );
//       }
    }
    else // must be curvilinear
    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      float_sw4 factor = 0.5/sqrt(mJ(i0,j0,k0));
      float_sw4 duxdq = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0));
      float_sw4 duydq = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0));
      float_sw4 duzdq = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0));
      float_sw4 duxdr = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0));
      float_sw4 duydr = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0));
      float_sw4 duzdr = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0));
      float_sw4 duxds = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1));
      float_sw4 duyds = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1));
      float_sw4 duzds = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1));
      float_sw4 duzdy = mMetric(1,i0,j0,k0)*duzdr+mMetric(3,i0,j0,k0)*duzds;
      float_sw4 duydz = mMetric(4,i0,j0,k0)*duyds;
      float_sw4 duxdz = mMetric(4,i0,j0,k0)*duxds;
      float_sw4 duzdx = mMetric(1,i0,j0,k0)*duzdq+mMetric(2,i0,j0,k0)*duzds;
      float_sw4 duydx = mMetric(1,i0,j0,k0)*duydq+mMetric(2,i0,j0,k0)*duyds;
      float_sw4 duxdy = mMetric(1,i0,j0,k0)*duxdr+mMetric(3,i0,j0,k0)*duxds;
//       if( m_xycomponent )
//       {
      uRec[0] = (duzdy-duydz)*factor;
      uRec[1] = (duxdz-duzdx)*factor;
      uRec[2] = (duydx-duxdy)*factor;
//       }
//       else
//       {
// 	 float_sw4 uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
// 	 float_sw4 uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
// 	 mRecordedUX.push_back( uew*factor );
// 	 mRecordedUY.push_back( uns*factor );
// 	 mRecordedUZ.push_back( -(duydx-duxdy)*factor );
//       }
    }
  } // end Curl
  else if(mode == TimeSeries::Strains )
  {
     uRec.resize(6);
    if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      float_sw4 factor = 1.0/(2*mGridSize[g0]);
      float_sw4 duydx = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0))*factor;
      float_sw4 duzdx = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0))*factor;
      float_sw4 duxdy = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0))*factor;
      float_sw4 duzdy = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0))*factor;
      float_sw4 duxdz = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))*factor;
      float_sw4 duydz = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))*factor;
      float_sw4 duxdx = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))*factor;
      float_sw4 duydy = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))*factor;
      float_sw4 duzdz = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1))*factor;
      uRec[0] = ( duxdx );
      uRec[1] = ( duydy );
      uRec[2] = ( duzdz );
      uRec[3] = ( 0.5*(duydx+duxdy) );
      uRec[4] = ( 0.5*(duzdx+duxdz) );
      uRec[5] = ( 0.5*(duydz+duzdy) );
   }
    else // must be curvilinear
   {
//       int i=m_i0, j=m_j0, k0=m_k00, g0=m_grid0;
      float_sw4 factor = 0.5/sqrt(mJ(i0,j0,k0));
      float_sw4 duxdq = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0));
      float_sw4 duydq = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0));
      float_sw4 duzdq = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0));
      float_sw4 duxdr = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0));
      float_sw4 duydr = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0));
      float_sw4 duzdr = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0));
      float_sw4 duxds = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1));
      float_sw4 duyds = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1));
      float_sw4 duzds = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1));
      float_sw4 duzdy = (mMetric(1,i0,j0,k0)*duzdr+mMetric(3,i0,j0,k0)*duzds)*factor;
      float_sw4 duydz = (mMetric(4,i0,j0,k0)*duyds)*factor;
      float_sw4 duxdz = (mMetric(4,i0,j0,k0)*duxds)*factor;
      float_sw4 duzdx = (mMetric(1,i0,j0,k0)*duzdq+mMetric(2,i0,j0,k0)*duzds)*factor;
      float_sw4 duydx = (mMetric(1,i0,j0,k0)*duydq+mMetric(2,i0,j0,k0)*duyds)*factor;
      float_sw4 duxdy = (mMetric(1,i0,j0,k0)*duxdr+mMetric(3,i0,j0,k0)*duxds)*factor;
      float_sw4 duxdx = ( mMetric(1,i0,j0,k0)*(U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))+
		       mMetric(2,i0,j0,k0)*(U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1)) )*factor;
      float_sw4 duydy = ( mMetric(1,i0,j0,k0)*(U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))+
		       mMetric(3,i0,j0,k0)*(U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1)) )*factor;
      float_sw4 duzdz = ( mMetric(4,i0,j0,k0)*(U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1)) )*factor;
      uRec[0] = ( duxdx );
      uRec[1] = ( duydy );
      uRec[2] = ( duzdz );
      uRec[3] = ( 0.5*(duydx+duxdy) );
      uRec[4] = ( 0.5*(duzdx+duxdz) );
      uRec[5] = ( 0.5*(duydz+duzdy) );
   }
  } // end Strains
  else if(mode == TimeSeries::DisplacementGradient )
  {
     uRec.resize(9);
     if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
     {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
	float_sw4 factor = 1.0/(2*mGridSize[g0]);
	float_sw4 duydx = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0))*factor;
	float_sw4 duzdx = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0))*factor;
	float_sw4 duxdy = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0))*factor;
	float_sw4 duzdy = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0))*factor;
	float_sw4 duxdz = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))*factor;
	float_sw4 duydz = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))*factor;
	float_sw4 duxdx = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))*factor;
	float_sw4 duydy = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))*factor;
	float_sw4 duzdz = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1))*factor;
	uRec[0] =  duxdx;
	uRec[1] =  duxdy;
	uRec[2] =  duxdz;
	uRec[3] =  duydx;
	uRec[4] =  duydy;
	uRec[5] =  duydz;
	uRec[6] =  duzdx;
	uRec[7] =  duzdy;
	uRec[8] =  duzdz;
     }
     else // must be curvilinear
     {
//       int i=m_i0, j=m_j0, k0=m_k00, g0=m_grid0;
	float_sw4 factor = 0.5/sqrt(mJ(i0,j0,k0));
	float_sw4 duxdq = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0));
	float_sw4 duydq = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0));
	float_sw4 duzdq = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0));
	float_sw4 duxdr = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0));
	float_sw4 duydr = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0));
	float_sw4 duzdr = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0));
	float_sw4 duxds = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1));
	float_sw4 duyds = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1));
	float_sw4 duzds = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1));
	float_sw4 duzdy = (mMetric(1,i0,j0,k0)*duzdr+mMetric(3,i0,j0,k0)*duzds)*factor;
	float_sw4 duydz = (mMetric(4,i0,j0,k0)*duyds)*factor;
	float_sw4 duxdz = (mMetric(4,i0,j0,k0)*duxds)*factor;
	float_sw4 duzdx = (mMetric(1,i0,j0,k0)*duzdq+mMetric(2,i0,j0,k0)*duzds)*factor;
	float_sw4 duydx = (mMetric(1,i0,j0,k0)*duydq+mMetric(2,i0,j0,k0)*duyds)*factor;
	float_sw4 duxdy = (mMetric(1,i0,j0,k0)*duxdr+mMetric(3,i0,j0,k0)*duxds)*factor;
	float_sw4 duxdx = ( mMetric(1,i0,j0,k0)*(U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))+
		       mMetric(2,i0,j0,k0)*(U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1)) )*factor;
	float_sw4 duydy = ( mMetric(1,i0,j0,k0)*(U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))+
		       mMetric(3,i0,j0,k0)*(U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1)) )*factor;
	float_sw4 duzdz = ( mMetric(4,i0,j0,k0)*(U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1)) )*factor;
	uRec[0] =  duxdx;
	uRec[1] =  duxdy;
	uRec[2] =  duxdz;
	uRec[3] =  duydx;
	uRec[4] =  duydy;
	uRec[5] =  duydz;
	uRec[6] =  duzdx;
	uRec[7] =  duzdy;
	uRec[8] =  duzdz;
     }

  } // end DisplacementGradient
  return;
}

//---------------------------------------------------------------------------
void EW::addSuperGridDamping(vector<Sarray> & a_Up, vector<Sarray> & a_U,
			     vector<Sarray> & a_Um, vector<Sarray> & a_Rho )
{
SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *up_ptr, *u_ptr, *um_ptr;
  
  int g;
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    u_ptr   = a_U[g].c_ptr();
    um_ptr  = a_Um[g].c_ptr();
    float_sw4* rho_ptr = a_Rho[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    if( m_sg_damping_order == 4 )
    {
       if( topographyExists() && g == mNumberOfGrids-1 )
       {
	  if( m_croutines )
	     addsgd4c_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, up_ptr, u_ptr, um_ptr,
			  rho_ptr, m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], m_sg_str_y[g],
			  mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g], m_supergrid_damping_coefficient );
          else
	     addsgd4c( &mDt, up_ptr, u_ptr, um_ptr, rho_ptr, m_sg_dc_x[g], m_sg_dc_y[g],
		        m_sg_str_x[g], m_sg_str_y[g], mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g],
		       &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
       }
       else
       {
	  if( m_croutines )
	     addsgd4_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, up_ptr, u_ptr, um_ptr, rho_ptr,
			 m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], 
			 m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
			 m_supergrid_damping_coefficient );
	  else
	     addsgd4( &mDt, &mGridSize[g], up_ptr, u_ptr, um_ptr, rho_ptr, m_sg_dc_x[g], m_sg_dc_y[g],
		      m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g], m_sg_corner_x[g],
		      m_sg_corner_y[g], m_sg_corner_z[g],
		      &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
       }
    }
    else if(  m_sg_damping_order == 6 )
    {
       if( topographyExists() && g == mNumberOfGrids-1 )
       {
	  if( m_croutines )
	     addsgd6c_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, up_ptr, u_ptr, um_ptr,
			  rho_ptr, m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], m_sg_str_y[g],
			  mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g], m_supergrid_damping_coefficient );
	  else
	     addsgd6c( &mDt, up_ptr, u_ptr, um_ptr, rho_ptr, m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g],
		       m_sg_str_y[g], mJ.c_ptr(), m_sg_corner_x[g], m_sg_corner_y[g],
		       &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
       }
       else
       {
	  if( m_croutines )
	     addsgd6_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, up_ptr, u_ptr, um_ptr, rho_ptr,
			 m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], 
			 m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
			 m_supergrid_damping_coefficient );
	  else
	     addsgd6( &mDt, &mGridSize[g], up_ptr, u_ptr, um_ptr, rho_ptr, m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g],
		      m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g], m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
		      &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
       }
    }
  }
  SYNC_STREAM;
}

//---------------------------------------------------------------------------
void EW::simpleAttenuation( vector<Sarray> & a_Up )
{
SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *up_ptr, cfreq, dt;
// Qs is saved in the EW object as mQs
// center frequency is called m_att_max_frecuency
// time step is called mDt  
  dt = mDt;
  cfreq = m_att_max_frequency;
  
  int g;
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    float_sw4* qs_ptr =mQs[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];

    if( m_croutines )
       satt_ci( up_ptr, qs_ptr, dt, cfreq,
		ifirst, ilast, jfirst, jlast, kfirst, klast );
    else
       satt( up_ptr, qs_ptr, &dt, &cfreq,
	     &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast );
  }
}

//-----------------------------------------------------------------------
void EW::enforceBCfreeAtt2( vector<Sarray>& a_Up, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                            vector<Sarray*>& a_AlphaVEp, vector<float_sw4 **>& a_BCForcing ) 
{
SW4_MARK_FUNCTION;
// AP: Apr. 3, 2017: Decoupled enforcement of the free surface bc with PC time stepping for memory variables
   int sg = usingSupergrid();
#ifdef ENABLE_CUDA
   using LOCAL_POL = 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernel<
       RAJA::statement::For<0, RAJA::cuda_threadblock_exec<16>, 
			    RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
						 RAJA::statement::Lambda<0> >>>>;
#else
   using LOCAL_POL = DEFAULT_LOOP2;
#endif
   SView *viewArray = viewArrayActual;
   SView *viewArray2 = viewArrayActual +  m_number_mechanisms;
   SView *viewArray3 = viewArrayActual +  2*m_number_mechanisms;
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
      int ifirst = m_iStart[g];
      int ilast  = m_iEnd[g];
      int jfirst = m_jStart[g];
      int jlast  = m_jEnd[g];
      int kfirst = m_kStart[g];
      int klast  = m_kEnd[g];
      float_sw4 h   = mGridSize[g];
      int topo = topographyExists() && g == mNumberOfGrids-1;

      SView &a_UpgV = a_Up[g].getview();
      SView &a_MugV = a_Mu[g].getview();
      SView &a_LambdagV = a_Lambda[g].getview();
     // SView &a_AlphaVEpgV = a_AlphaVEp[g]->getview(); // PBUGS
      for(int m=0;m <  m_number_mechanisms; m++){
	viewArray[m] = a_AlphaVEp[g][m].getview();
	viewArray2[m] = mMuVE[g][m].getview();
	viewArray3[m] = mLambdaVE[g][m].getview();
      }
      float_sw4 *m_sg_str_xg = m_sg_str_x[g];
      float_sw4 *m_sg_str_yg = m_sg_str_y[g];
      float_sw4 *m_sg_str_zg = m_sg_str_z[g];
      PREFETCH(m_sg_str_x[g]);
      PREFETCH(m_sg_str_y[g]);
      //      PREFETCH(viewArray);
      if( m_bcType[g][4] == bStressFree && !topo ) // Cartesian case
      {
	SW4_MARK_BEGIN("enforceBCfreeAtt2::SET1");
	 //	    if( m_croutines )
	 //	       memvarforcesurf_ci( ifirst, ilast, jfirst, jlast, k, mf, a_t, om,
	 //				   cv, ph, mOmegaVE[0], mDt, h, m_zmin[g] );
	 //	    else
	 //	       memvarforcesurf( &ifirst, &ilast, &jfirst, &jlast, &k, mf, &a_t, &om,
	 //				&cv, &ph, &mOmegaVE[0], &mDt, &h, &m_zmin[g] );
	 //	 }
	 //	 else
	 //	    memforce.set_value(0.0);
	 //const float_sw4 i6  = 1.0/6;
	 const float_sw4 d4a = 2.0/3;
	 const float_sw4 d4b =-1.0/12;
	 float_sw4* forcing = a_BCForcing[g][4];
         ASSERT_MANAGED(forcing);
         ASSERT_MANAGED(m_sbop);
         ASSERT_MANAGED(m_sbop_no_gp);
	 int ni = (ilast-ifirst+1);
         int m_number_mechanisms_local = m_number_mechanisms; // Because the lambda does not capture member arrays and variables.
	 float_sw4* lm_sbop = m_sbop; // Because the lambda does not capture member arrays and variables.
	 float_sw4* lm_sbop_no_gp = m_sbop_no_gp; // Because the lambda does not capture member arrays and variables.
	 //#pragma omp parallel
	 {
	    //	 float_sw4* r1 = new float_sw4[m_number_mechanisms];
	    //	 float_sw4* r2 = new float_sw4[m_number_mechanisms];
	    //	 float_sw4* r3 = new float_sw4[m_number_mechanisms];
	    //	 float_sw4* cof = new float_sw4[m_number_mechanisms];
	   //#pragma omp for
	   //for( int j=jfirst+2 ; j<=jlast-2 ; j++ )
	   //for( int i=ifirst+2 ; i<=ilast-2 ; i++ )
	   //{
	      RAJA::RangeSegment j_range(jfirst+2,jlast-1);
	      RAJA::RangeSegment i_range(ifirst+2,ilast-1 );
	      RAJA::kernel<LOCAL_POL>(
					 RAJA::make_tuple(j_range,i_range),
					 [=]RAJA_DEVICE (int j,int i) {
	       float_sw4 g1, g2, g3, acof, bcof;
	       int ind = i-ifirst + ni*(j-jfirst);
	       float_sw4 a4ci, b4ci, a4cj, b4cj;
	       a4ci = a4cj = d4a;
	       b4ci = b4cj = d4b;
	       if( sg == 1 )
	       {
		  a4ci = d4a*m_sg_str_xg[i-ifirst];
		  b4ci = d4b*m_sg_str_xg[i-ifirst];
		  a4cj = d4a*m_sg_str_yg[j-jfirst];
		  b4cj = d4b*m_sg_str_yg[j-jfirst];
	       }
// this would be more efficient if done in Fortran
// first add interior elastic terms (use ghost point stencils)
	       g1 =  h*forcing[3*ind]
                  - a_MugV(i,j,1)*
                  (lm_sbop[1]*a_UpgV(1,i,j,1) + lm_sbop[2]*a_UpgV(1,i,j,2)
                   + lm_sbop[3]*a_UpgV(1,i,j,3) + lm_sbop[4]*a_UpgV(1,i,j,4) + lm_sbop[5]*a_UpgV(1,i,j,5)
                   + a4ci*(a_UpgV(3,i+1,j,1)-a_UpgV(3,i-1,j,1))
                   + b4ci*(a_UpgV(3,i+2,j,1)-a_UpgV(3,i-2,j,1)) );

	       g2 =  h*forcing[3*ind+1]
                  - a_MugV(i,j,1)*
                  ( lm_sbop[1]*a_UpgV(2,i,j,1) + lm_sbop[2]*a_UpgV(2,i,j,2)
                    + lm_sbop[3]*a_UpgV(2,i,j,3) + lm_sbop[4]*a_UpgV(2,i,j,4) + lm_sbop[5]*a_UpgV(2,i,j,5)
                    + a4cj*(a_UpgV(3,i,j+1,1)-a_UpgV(3,i,j-1,1))
                    + b4cj*(a_UpgV(3,i,j+2,1)-a_UpgV(3,i,j-2,1)) );

	       g3 =  h*forcing[3*ind+2]
                  - (2*a_MugV(i,j,1)+a_LambdagV(i,j,1))*
                  ( lm_sbop[1]*a_UpgV(3,i,j,1) + lm_sbop[2]*a_UpgV(3,i,j,2)
                    + lm_sbop[3]*a_UpgV(3,i,j,3) + lm_sbop[4]*a_UpgV(3,i,j,4)  + lm_sbop[5]*a_UpgV(3,i,j,5) )
                  - a_LambdagV(i,j,1)*
                  ( a4ci*(a_UpgV(1,i+1,j,1)-a_UpgV(1,i-1,j,1))
                    + b4ci*(a_UpgV(1,i+2,j,1)-a_UpgV(1,i-2,j,1)) 
                    + a4cj*(a_UpgV(2,i,j+1,1)-a_UpgV(2,i,j-1,1))
                    + b4cj*(a_UpgV(2,i,j+2,1)-a_UpgV(2,i,j-2,1)) );
               
	       acof = a_MugV(i,j,1);
	       bcof = 2*a_MugV(i,j,1)+a_LambdagV(i,j,1);
	       for( int a=0 ; a < m_number_mechanisms_local ; a++ )
	       {
// this would be more efficient if done in Fortran
// Add in visco-elastic contributions (NOT using ghost points)
// mu*( a1_z + a3_x )
		  g1 = g1 + viewArray2[a](i,j,1)*
                     ( lm_sbop_no_gp[0]*viewArray[a](1,i,j,0) + lm_sbop_no_gp[1]*viewArray[a](1,i,j,1)
                       + lm_sbop_no_gp[2]*viewArray[a](1,i,j,2)
                       + lm_sbop_no_gp[3]*viewArray[a](1,i,j,3)
                       + lm_sbop_no_gp[4]*viewArray[a](1,i,j,4) + lm_sbop_no_gp[5]*viewArray[a](1,i,j,5)
                       + a4ci*(viewArray[a](3,i+1,j,1)-viewArray[a](3,i-1,j,1))
                       + b4ci*(viewArray[a](3,i+2,j,1)-viewArray[a](3,i-2,j,1)) );
// mu*( a2_z + a3_y )
		  g2 = g2 + viewArray2[a](i,j,1)*
                     ( lm_sbop_no_gp[0]*viewArray[a](2,i,j,0) + lm_sbop_no_gp[1]*viewArray[a](2,i,j,1)
                       + lm_sbop_no_gp[2]*viewArray[a](2,i,j,2)
                       + lm_sbop_no_gp[3]*viewArray[a](2,i,j,3)
                       + lm_sbop_no_gp[4]*viewArray[a](2,i,j,4) + lm_sbop_no_gp[5]*viewArray[a](2,i,j,5)
                       + a4cj*(viewArray[a](3,i,j+1,1)-viewArray[a](3,i,j-1,1))
                       + b4cj*(viewArray[a](3,i,j+2,1)-viewArray[a](3,i,j-2,1)) );
// (2*mu + lambda)*( a3_z ) + lambda*( a1_x + a2_y )
		  g3 = g3 + (2*viewArray2[a](i,j,1)+viewArray3[a](i,j,1))*
                     (lm_sbop_no_gp[0]*viewArray[a](3,i,j,0)
                      + lm_sbop_no_gp[1]*viewArray[a](3,i,j,1) + lm_sbop_no_gp[2]*viewArray[a](3,i,j,2)
                      + lm_sbop_no_gp[3]*viewArray[a](3,i,j,3)
                      + lm_sbop_no_gp[4]*viewArray[a](3,i,j,4) + lm_sbop_no_gp[5]*viewArray[a](3,i,j,5) )
                     + viewArray3[a](i,j,1)*
                     ( a4ci*(viewArray[a](1,i+1,j,1)-viewArray[a](1,i-1,j,1))
                       + b4ci*(viewArray[a](1,i+2,j,1)-viewArray[a](1,i-2,j,1)) 
                       + a4cj*(viewArray[a](2,i,j+1,1)-viewArray[a](2,i,j-1,1))
                       + b4cj*(viewArray[a](2,i,j+2,1)-viewArray[a](2,i,j-2,1)) );
	       } // end for all mechanisms
// solve for the ghost point value of Up (stencil uses ghost points for the elastic variable)
	       a_UpgV(1,i,j,0) = g1/(acof*lm_sbop[0]);
	       a_UpgV(2,i,j,0) = g2/(acof*lm_sbop[0]);
	       a_UpgV(3,i,j,0) = g3/(bcof*lm_sbop[0]);
					 }); //SYNC_STREAM;

	 }
	 SW4_MARK_END("enforceBCfreeAtt2::SET1");
      }    // end if bcType[g][4] == bStressFree   
      if( m_bcType[g][5] == bStressFree  )
      {
	SW4_MARK_BEGIN("enforceBCfreeAtt2::SET 2");
         int nk=m_global_nz[g];
	 //const float_sw4 i6  = 1.0/6;
	 const float_sw4 d4a = 2.0/3;
	 const float_sw4 d4b =-1.0/12;
	 float_sw4* forcing = a_BCForcing[g][5];
	 int ni = (ilast-ifirst+1);
#pragma omp parallel
	 {
#pragma omp for
	 for( int j=jfirst+2 ; j<=jlast-2 ; j++ )
	    for( int i=ifirst+2 ; i<=ilast-2 ; i++ )
	    {
	       float_sw4 g1, g2, g3, acof, bcof;
	       int ind = i-ifirst + ni*(j-jfirst);
	       float_sw4 a4ci, b4ci, a4cj, b4cj;
	       a4ci = a4cj = d4a;
	       b4ci = b4cj = d4b;
	       if( sg == 1 )
	       {
		  a4ci = d4a*m_sg_str_x[g][i-ifirst];
		  b4ci = d4b*m_sg_str_x[g][i-ifirst];
		  a4cj = d4a*m_sg_str_y[g][j-jfirst];
		  b4cj = d4b*m_sg_str_y[g][j-jfirst];
	       }
// add in contributions from elastic terms
	       g1 = h*forcing[3*ind] - a_Mu[g](i,j,nk)*
                  (-m_sbop[1]*a_Up[g](1,i,j,nk) - m_sbop[2]*a_Up[g](1,i,j,nk-1) 
                   -m_sbop[3]*a_Up[g](1,i,j,nk-2) - m_sbop[4]*a_Up[g](1,i,j,nk-3)
                   -m_sbop[5]*a_Up[g](1,i,j,nk-4)
                   + a4ci*(a_Up[g](3,i+1,j,nk)-a_Up[g](3,i-1,j,nk))
                   + b4ci*(a_Up[g](3,i+2,j,nk)-a_Up[g](3,i-2,j,nk)) );

	       g2 = h*forcing[3*ind+1] - a_Mu[g](i,j,nk)*
                  (-m_sbop[1]*a_Up[g](2,i,j,nk) - m_sbop[2]*a_Up[g](2,i,j,nk-1)
                   -m_sbop[3]*a_Up[g](2,i,j,nk-2) - m_sbop[4]*a_Up[g](2,i,j,nk-3)
                   -m_sbop[5]*a_Up[g](2,i,j,nk-4)
                   + a4cj*(a_Up[g](3,i,j+1,nk)-a_Up[g](3,i,j-1,nk))
                   + b4cj*(a_Up[g](3,i,j+2,nk)-a_Up[g](3,i,j-2,nk)) );

	       g3 =  h*forcing[3*ind+2] - (2*a_Mu[g](i,j,nk)+a_Lambda[g](i,j,nk))*
                  ( -m_sbop[1]*a_Up[g](3,i,j,nk) - m_sbop[2]*a_Up[g](3,i,j,nk-1) 
                    -m_sbop[3]*a_Up[g](3,i,j,nk-2) - m_sbop[4]*a_Up[g](3,i,j,nk-3)
                    -m_sbop[5]*a_Up[g](3,i,j,nk-4) )
                  - a_Lambda[g](i,j,nk)*
                  ( a4ci*(a_Up[g](1,i+1,j,nk)-a_Up[g](1,i-1,j,nk))
                    + b4ci*(a_Up[g](1,i+2,j,nk)-a_Up[g](1,i-2,j,nk)) 
                    + a4cj*(a_Up[g](2,i,j+1,nk)-a_Up[g](2,i,j-1,nk))
                    + b4cj*(a_Up[g](2,i,j+2,nk)-a_Up[g](2,i,j-2,nk)) );
               
	       acof = a_Mu[g](i,j,nk);
	       bcof = 2*a_Mu[g](i,j,nk)+a_Lambda[g](i,j,nk);
	       for( int a=0 ; a < m_number_mechanisms ; a++ )
	       {
// visco-elastic terms (NOT using ghost points)
		  g1 = g1 + mMuVE[g][a](i,j,nk)*
                     (-m_sbop_no_gp[0]*a_AlphaVEp[g][a](1,i,j,nk+1) - m_sbop_no_gp[1]*a_AlphaVEp[g][a](1,i,j,nk)
                      -m_sbop_no_gp[2]*a_AlphaVEp[g][a](1,i,j,nk-1)
                      -m_sbop_no_gp[3]*a_AlphaVEp[g][a](1,i,j,nk-2) - m_sbop_no_gp[4]*a_AlphaVEp[g][a](1,i,j,nk-3)
                      -m_sbop_no_gp[5]*a_AlphaVEp[g][a](1,i,j,nk-4)
                      +a4ci*(a_AlphaVEp[g][a](3,i+1,j,nk)-a_AlphaVEp[g][a](3,i-1,j,nk))
                      +b4ci*(a_AlphaVEp[g][a](3,i+2,j,nk)-a_AlphaVEp[g][a](3,i-2,j,nk)));

		  g2 = g2 + mMuVE[g][a](i,j,nk)*
                     ( - m_sbop_no_gp[0]*a_AlphaVEp[g][a](2,i,j,nk+1) -m_sbop_no_gp[1]*a_AlphaVEp[g][a](2,i,j,nk)
                       - m_sbop_no_gp[2]*a_AlphaVEp[g][a](2,i,j,nk-1)
                       - m_sbop_no_gp[3]*a_AlphaVEp[g][a](2,i,j,nk-2) -m_sbop_no_gp[4]*a_AlphaVEp[g][a](2,i,j,nk-3)
                       - m_sbop_no_gp[5]*a_AlphaVEp[g][a](2,i,j,nk-4)
                       + a4cj*(a_AlphaVEp[g][a](3,i,j+1,nk)-a_AlphaVEp[g][a](3,i,j-1,nk))
                       + b4cj*(a_AlphaVEp[g][a](3,i,j+2,nk)-a_AlphaVEp[g][a](3,i,j-2,nk)) );
                                                  
		  g3 = g3 + (2*mMuVE[g][a](i,j,nk)+mLambdaVE[g][a](i,j,nk))*
                     (- m_sbop_no_gp[0]*a_AlphaVEp[g][a](3,i,j,nk+1) - m_sbop_no_gp[1]*a_AlphaVEp[g][a](3,i,j,nk)
                      - m_sbop_no_gp[2]*a_AlphaVEp[g][a](3,i,j,nk-1)
                      - m_sbop_no_gp[3]*a_AlphaVEp[g][a](3,i,j,nk-2) - m_sbop_no_gp[4]*a_AlphaVEp[g][a](3,i,j,nk-3)
                      - m_sbop_no_gp[5]*a_AlphaVEp[g][a](3,i,j,nk-4))
                     + mLambdaVE[g][a](i,j,nk)*
                     (a4ci*(a_AlphaVEp[g][a](1,i+1,j,nk)-a_AlphaVEp[g][a](1,i-1,j,nk))
                      + b4ci*(a_AlphaVEp[g][a](1,i+2,j,nk)-a_AlphaVEp[g][a](1,i-2,j,nk)) 
                      + a4cj*(a_AlphaVEp[g][a](2,i,j+1,nk)-a_AlphaVEp[g][a](2,i,j-1,nk))
                      + b4cj*(a_AlphaVEp[g][a](2,i,j+2,nk)-a_AlphaVEp[g][a](2,i,j-2,nk)) );
	       }
               // solve for the ghost point value of Up (using the ghost point stencil)
	       a_Up[g](1,i,j,nk+1) = g1/(-m_sbop[0]*acof);
	       a_Up[g](2,i,j,nk+1) = g2/(-m_sbop[0]*acof);
	       a_Up[g](3,i,j,nk+1) = g3/(-m_sbop[0]*bcof);
	    }
	 }
	 SW4_MARK_END("enforceBCfreeAtt2::SET 2");
      }// end if bcType[g][5] == bStressFree
   
// all the curvilinear code needs to be overhauled
      if( m_bcType[g][4] == bStressFree && topo && g == mNumberOfGrids-1 )
      {
	SW4_MARK_BEGIN("enforceBCfreeAtt2::SET 3");
         float_sw4* mu_p = a_Mu[g].c_ptr();
         float_sw4* la_p = a_Lambda[g].c_ptr();
         float_sw4* up_p = a_Up[g].c_ptr();
         int side = 5;
         int nz = m_global_nz[g];
         //int ghno = 0;
         //char op = '-';
	 float_sw4* forcing = a_BCForcing[g][4];
	 int usesg = usingSupergrid() ? 1 : 0;

// make a local copy of the boundary forcing array to simplify access
         Sarray bforcerhs(3,ifirst,ilast,jfirst,jlast,1,1,__FILE__,__LINE__);
         bforcerhs.assign(forcing,0);

	 //	 if( m_croutines )
	 //	    addbstressc_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
	 //			    nz, up_p, mu_p, la_p, bforcerhs.c_ptr(), mMetric.c_ptr(), 
	 //			    side, m_sbop, op, ghno, usesg, m_sg_str_x[g], m_sg_str_y[g] );
	 //	 else
	 //	    addbstressc( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
	 //			 &nz, up_p, mu_p, la_p, bforcerhs.c_ptr(), mMetric.c_ptr(), 
	 //			 &side, m_sbop, &op, &ghno, &usesg, m_sg_str_x[g], m_sg_str_y[g] );

         for( int a = 0 ; a < m_number_mechanisms ; a++ )
         {
	   SW4_MARK_BEGIN("CPTR");
            float_sw4* mu_ve_p   = mMuVE[g][a].c_ptr();
            float_sw4* lave_p   = mLambdaVE[g][a].c_ptr();
            float_sw4* alphap_p = a_AlphaVEp[g][a].c_ptr();
	    SW4_MARK_END("CPTR");
	    //std::cout<<"Pointers "<<mu_ve_p<<" "<<lave_p<<" "<<alphap_p<<"\n";
            // This function adds the visco-elastic boundary stresses to bforcerhs
	    if(  m_croutines )
	       ve_bndry_stress_curvi_ci( ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
                                  alphap_p, mu_ve_p, lave_p, bforcerhs.c_ptr(), mMetric.c_ptr(), side,
                                  m_sbop_no_gp, usesg, m_sg_str_x[g], m_sg_str_y[g] ); // no ghost points here
	    else
	       ve_bndry_stress_curvi(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
				     alphap_p, mu_ve_p, lave_p, bforcerhs.c_ptr(), mMetric.c_ptr(), &side,
				     m_sbop_no_gp, &usesg, m_sg_str_x[g], m_sg_str_y[g] ); // no ghost points here
         } // end for a...
         
// update GHOST POINT VALUES OF UP
	 if( m_croutines )
	    att_free_curvi_ci( ifirst, ilast, jfirst, jlast, kfirst, klast,
                         up_p, mu_p, la_p, bforcerhs.c_ptr(), mMetric.c_ptr(), m_sbop, // use ghost points
                         usesg, m_sg_str_x[g], m_sg_str_y[g] );
	 else
	    att_free_curvi (&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			    up_p, mu_p, la_p, bforcerhs.c_ptr(), mMetric.c_ptr(), m_sbop, // use ghost points
			    &usesg, m_sg_str_x[g], m_sg_str_y[g] );
	 SW4_MARK_END("enforceBCfreeAtt2::SET 3");
	 SYNC_STREAM;
      } // end if bcType[g][4] == bStressFree && topography
     
   }  // end for g=0,.
   //::operator delete[](viewArray,Managed);
}
