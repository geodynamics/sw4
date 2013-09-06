#include "EW.h"
#include "impose_cartesian_bc.h"
#include "DataPatches.h"

//--------------------------------------------------------------------
void EW::solve_allpars( vector<Source*> & a_Sources, vector<Sarray> &a_Rho,
			vector<Sarray> &a_Mu, vector<Sarray> &a_Lambda,
			vector<TimeSeries*> & a_TimeSeries, 
                        vector<Sarray>& U, vector<Sarray>& Um,
			vector<DataPatches*>& Upred_saved_sides,
			vector<DataPatches*>& Ucorr_saved_sides )
{
// solution arrays
  vector<Sarray> F, Lu, Uacc, Up;
  vector<Sarray*> AlphaVE, AlphaVEm, AlphaVEp;
// vectors of pointers to hold boundary forcing arrays in each grid
  vector<double **> BCForcing;

  BCForcing.resize(mNumberOfGrids);
  F.resize(mNumberOfGrids);
  Lu.resize(mNumberOfGrids);
  Uacc.resize(mNumberOfGrids);
  Up.resize(mNumberOfGrids);
  Um.resize(mNumberOfGrids);
  U.resize(mNumberOfGrids);
// Allocate pointers, even if attenuation not used, for avoid segfault in parameter list with mMuVE[g], etc...
  AlphaVE.resize(mNumberOfGrids);
  AlphaVEm.resize(mNumberOfGrids);
  AlphaVEp.resize(mNumberOfGrids);
  if (m_use_attenuation)
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
    BCForcing[g] = new double *[6];
    for (int side=0; side < 6; side++)
    {
      BCForcing[g][side]=NULL;
      if (m_bcType[g][side] == bStressFree || m_bcType[g][side] == bDirichlet || m_bcType[g][side] == bSuperGrid)
      {
    	BCForcing[g][side] = new double[3*m_NumberOfBCPoints[g][side]];
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

    if (m_use_attenuation)
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

// Allocate boundary sides
  for( int g=0 ; g < mNumberOfGrids ; g++ )
  {
     stringstream procno;
     procno << m_myRank << "." << g ; 
     //     string logname(getlogin());

// Local disks on LC seem to be setup with directory /tmp/username when user username starts a job
     string upred_name = mTempPath + "upred" + procno.str() + ".bin";
     string ucorr_name = mTempPath + "ucorr" + procno.str() + ".bin";
      //     string upred_name = "/tmp/" + logname + "/upred" + procno.str() + ".bin";
      //     string ucorr_name = "/tmp/" + logname + "/ucorr" + procno.str() + ".bin";
     int imin, imax, jmin, jmax, kmax;
     if( m_iStartAct[g] <= m_iEndAct[g] && m_iStartAct[g] <= m_iEndAct[g] && m_iStartAct[g] <= m_iEndAct[g]  )
     {
	imin = m_iStartAct[g]-1;
	imax = m_iEndAct[g]+1;
	jmin = m_jStartAct[g]-1;
	jmax = m_jEndAct[g]+1;
	kmax = m_kEndAct[g]+1;
     }
     else
     {
	// empty active domain
        imin =  0;
	imax = -1;
	jmin =  0;
	jmax = -1;
	kmax = -1;
     }
     Upred_saved_sides[g] = new DataPatches( upred_name.c_str() ,U[g],imin,imax,jmin,jmax,kmax,2,20,mDt );
     Ucorr_saved_sides[g] = new DataPatches( ucorr_name.c_str() ,U[g],imin,imax,jmin,jmax,kmax,2,20,mDt );
     //     cout << "sides saved for i=[" << imin << " , " << imax << "] j=[" << jmin << " , " << jmax << "] k=[" << 1 << " , " << kmax << "]"<< endl;
     size_t maxsizeloc = Upred_saved_sides[g]->get_noofpoints();
     size_t maxsize;
     int mpisizelong, mpisizelonglong, mpisizeint;
     MPI_Type_size(MPI_LONG,&mpisizelong );
     MPI_Type_size(MPI_LONG_LONG,&mpisizelonglong );
     MPI_Type_size(MPI_INT,&mpisizeint );
     if( sizeof(size_t) == mpisizelong )
	MPI_Allreduce( &maxsizeloc, &maxsize, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD );
     else if( sizeof(size_t) == mpisizelonglong )
	MPI_Allreduce( &maxsizeloc, &maxsize, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD );
     else if( sizeof(size_t) == mpisizeint )
	MPI_Allreduce( &maxsizeloc, &maxsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

     if( proc_zero() )
	cout << "Maximum temporary file size on grid " << g << " is " << maxsize << " doubles for each time step "<<endl;

     if( !mQuiet && mVerbose >= 5 && proc_zero() )
	cout << "  Temporary files " << upred_name << " and " << ucorr_name << " will hold " <<
	   Upred_saved_sides[g]->get_noofpoints() << " values each, for each time step";
  }
// Set the number of time steps, allocate the recording arrays, and set reference time in all time series objects  
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
  for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
     mImageFiles[fIndex]->initializeTime();
   
// the Source objects get discretized into GridPointSource objects
  vector<GridPointSource*> point_sources;

// Transfer source terms to each individual grid as point sources at grid points.
  for( unsigned int i=0 ; i < a_Sources.size() ; i++ )
      a_Sources[i]->set_grid_point_sources4( this, point_sources );

 // Debug
  //    for( int i=0 ; i < point_sources.size() ; i++ )
  //       cout << *point_sources[i] << endl;

// modify the time functions if prefiltering is enabled
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
  } // end if prefiltering

  bool output_timefunc = true;
  if( output_timefunc )
  {
     int has_source_id=-1, has_source_max;
     if( point_sources.size() > 0 )
	has_source_id = m_myRank;

     MPI_Allreduce( &has_source_id, &has_source_max, 1, MPI_INT, MPI_MAX, m_cartesian_communicator );
     if( m_myRank == has_source_max )
     {
       if (!mQuiet)
	 printf("Saving one discretized time function\n");

//building the file name...
       string filename;
       if( mPath != "." )
	 filename += mPath;
       filename += "g1.dat";	 

       FILE *tf=fopen(filename.c_str(),"w");
       double t;
       double gt, gt1, gt2;
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

// this is the time level
  double t=mTstart;
  
// Set up timers
  double time_start_solve = MPI_Wtime();
  double time_measure[7];
  double time_sum[7]={0,0,0,0,0,0,0};
  double bc_time_measure[5]={0,0,0,0,0};

  int beginCycle = 1; // also set in setupRun(), perhaps make member variable?

// Assign initial data
//  cout << getRank() << "before initial data" << endl;
  initialData(mTstart, U, AlphaVE);
  //  cout << getRank() << "U given initial data" << endl;
  initialData(mTstart-mDt, Um, AlphaVEm );
  //  cout << getRank() << "Um given initial data" << endl;
  //  U[0].save_to_disk("ustart.bin");

  if ( !mQuiet && mVerbose && proc_zero() )
     cout << "  Initial data has been assigned" << endl;

// save any images for cycle = 0 (initial data) ?
  update_images( 0, t, U, Um, Up, a_Rho, a_Mu, a_Lambda, a_Sources, 1 );

  //  cout << getRank() << "images updated" << endl;

// enforce bc on initial data
// communicate across processor boundaries
  for(int g=0 ; g < mNumberOfGrids ; g++ )
     communicate_array( U[g], g );
  //  cout << getRank() << "done communicate array" << endl;
// boundary forcing
  cartesian_bc_forcing( t, BCForcing, a_Sources );
  //  cout << getRank() << "done bc_forcing" << endl;
// enforce boundary condition
  enforceBC( U, a_Mu, a_Lambda, t, BCForcing );   

  //  cout << getRank() << "done bc U" << endl;
// Um
// communicate across processor boundaries
  for(int g=0 ; g < mNumberOfGrids ; g++ )
     communicate_array( Um[g], g );
// boundary forcing
  cartesian_bc_forcing( t-mDt, BCForcing, a_Sources );
// enforce boundary condition
  enforceBC( Um, a_Mu, a_Lambda, t-mDt, BCForcing );

  //  cout << getRank() << "Done bc Um" << endl;

  if (m_twilight_forcing)
  {
// more testing
    if ( proc_zero() )
    {
      printf("About to call exactSolTwilight\n");
    }
// check the accuracy of the initial data, store exact solution in Up, ignore AlphaVE
    exactSol( t, Up, AlphaVE, a_Sources );
    double errInf, errL2;
    normOfDifferenceGhostPoints( Up, U, errInf, errL2 );
    if ( proc_zero() )
      printf("\n Ghost point errors: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);
  }

  if( m_moment_test )
    test_sources( point_sources, a_Sources, F );

  if ( !mQuiet && proc_zero() )
    cout << "  Begin time stepping..." << endl;

// save initial data on receiver records
  vector<double> uRec;
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
  //  cout << getRank() << "Done receiver save" << endl;

  FILE *lf=NULL;
// open file for saving norm of error
  if ( (m_lamb_test || m_point_source_test || m_rayleigh_wave_test ) && proc_zero() )
  {
    string path=getOutputPath();

    stringstream fileName;
    if( path != "." )
      fileName << path;
    
    if (m_lamb_test)
      fileName << "LambErr.txt";
    else if (m_point_source_test)
      fileName << "PointSourceErr.txt";
    else
      fileName << "RayleighErr.txt";
    lf = fopen(fileName.str().c_str(),"w");
  }
  //  // DEBUG
  //     for( int s = 0 ; s < point_sources.size() ; s++ )
  //        point_sources[s]->print_info();
  //  cout << getRank() << "pushing initial data " << endl;    

  for( int g=0 ; g < mNumberOfGrids ; g++ )
  {
     Upred_saved_sides[g]->push( Um[g], -1 );
     Upred_saved_sides[g]->push( U[g], 0 );
     Ucorr_saved_sides[g]->push( Um[g], -1 );
     Ucorr_saved_sides[g]->push( U[g], 0 );
  }
// Begin time stepping loop
  for( int currentTimeStep = beginCycle; currentTimeStep <= mNumberOfTimeSteps; currentTimeStep++)
  {    
     //     cout << getRank() << "Starting time step " << currentTimeStep << endl;
    time_measure[0] = MPI_Wtime();

// all types of forcing...
    Force( t, F, point_sources );
      
    if( m_checkfornan )
    {
       check_for_nan( F, 1, "F" );
       check_for_nan( U, 1, "U" );
    }

// evaluate right hand side
    evalRHS( U, a_Mu, a_Lambda, Lu, AlphaVE ); // save Lu in composite grid 'Lu'

// take predictor step, store in Up
    evalPredictor(Up, U, Um, a_Rho, Lu, F );    

    time_measure[1] = MPI_Wtime();
    time_measure[2] = MPI_Wtime();

// communicate across processor boundaries
    for(int g=0 ; g < mNumberOfGrids ; g++ )
      communicate_array( Up[g], g );
// calculate boundary forcing at time t+mDt
    cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );
// update ghost points in Up
    enforceBC( Up, a_Mu, a_Lambda, t+mDt, BCForcing );
// interpolate across mesh refinement boundaries (?)
 //    check_consintp( Up[0], Up[1], AlphaVEp[0], AlphaVEp[1] );
    time_measure[3] = time_measure[4] = MPI_Wtime();

// get 4th order in time
    if (mOrder == 4)
    {
       Force_tt( t, F, point_sources );
      
       evalDpDmInTime( Up, U, Um, Uacc ); // store result in Uacc
       for( int g=0 ; g < mNumberOfGrids ; g++ )
	  Upred_saved_sides[g]->push( Uacc[g], currentTimeStep );
      
       evalRHS( Uacc, a_Mu, a_Lambda, Lu, AlphaVEm );

       evalCorrector( Up, a_Rho, Lu, F );
      
// add in super-grid damping terms
       if (usingSupergrid())
       {
	  addSuperGridDamping( Up, U, Um, a_Rho );
       }

       time_measure[4] = MPI_Wtime();

// also check out EW::update_all_boundaries 
// communicate across processor boundaries
       for(int g=0 ; g < mNumberOfGrids ; g++ )
	  communicate_array( Up[g], g );
// calculate boundary forcing at time t+mDt (do we really need to call this fcn again???)
       cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );
// update ghost points in Up
       enforceBC( Up, a_Mu, a_Lambda, t+mDt, BCForcing );
       for( int g=0 ; g < mNumberOfGrids ; g++ )
	  Ucorr_saved_sides[g]->push( Up[g], currentTimeStep );
    }
    
    if( m_checkfornan )
       check_for_nan( F, 1, "Up" );

// increment time
    t += mDt;

    time_measure[5] = MPI_Wtime();	  

// periodically, print time stepping info to stdout
    printTime( currentTimeStep, t, currentTimeStep == mNumberOfTimeSteps ); 

// Images have to be written before the solution arrays are cycled, because both Up and Um are needed
// to compute a centered time derivative
//
// AP: Note to self: Any quantity related to velocities will be lagged by one time step
//
    update_images( currentTimeStep, t, Up, U, Um, a_Rho, a_Mu, a_Lambda,
		   a_Sources, currentTimeStep == mNumberOfTimeSteps );
    
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

// // Energy evaluation, requires all three time levels present, do before cycle arrays.
    if( m_energy_test )
       compute_energy( mDt, currentTimeStep == mNumberOfTimeSteps, Um, U, Up, currentTimeStep  );

// cycle the solution arrays
    cycleSolutionArrays(Um, U, Up, AlphaVEm, AlphaVE, AlphaVEp);

// evaluate error for some test cases
    if (m_lamb_test || m_point_source_test || m_rayleigh_wave_test )
    {
      double errInf=0, errL2=0, solInf=0, solL2=0;
      exactSol( t, Up, AlphaVE, a_Sources ); // store exact solution in Up

      if (m_lamb_test)
	normOfSurfaceDifference( Up, U, errInf, errL2, solInf, solL2, a_Sources);
      else if (m_point_source_test || m_rayleigh_wave_test)
	normOfDifference( Up, U, errInf, errL2, solInf, a_Sources );

      if ( proc_zero() )
// output time, Linf-err, Linf-sol-err
	fprintf(lf, "%e %15.7e %15.7e %15.7e\n", t, errInf, errL2, solInf);
    }

    time_measure[6] = MPI_Wtime();	  	
    time_sum[0] += time_measure[1]-time_measure[0] + time_measure[4]-time_measure[3]; // step
    time_sum[1] += time_measure[2]-time_measure[1] + time_measure[5]-time_measure[4]; // bcs
    time_sum[2] += time_measure[6]-time_measure[5]; // image & sac
    time_sum[3] += time_measure[6]-time_measure[0];//  total
  } // end time stepping loop

  if ( !mQuiet && proc_zero() )
    cout << "  Time stepping finished..." << endl;

// close error file for Lamb's test
  if ((m_lamb_test || m_point_source_test || m_rayleigh_wave_test) && proc_zero() )
  {
     fclose(lf);
     printf("**** Closed file with solution errors for testing\n");
  }
  double time_end_solve = MPI_Wtime();
  print_execution_time( time_start_solve, time_end_solve, "solver phase" );

   if( m_output_detailed_timing )
     print_execution_times( time_sum );

// check the accuracy of the final solution, store exact solution in Up, ignore AlphaVE
   if( exactSol( t, Up, AlphaVE, a_Sources ) )
   {
      double errInf, errL2, solInf, solL2;
// depending on the test case, we should compare in the interior, or only on the surface
      if (m_lamb_test)
	 normOfSurfaceDifference( Up, U, errInf, errL2, solInf, solL2, a_Sources);
      else
	 normOfDifference( Up, U, errInf, errL2, solInf, a_Sources );

      if ( proc_zero() )
	 printf("\n Final solution errors: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);
   }
   finalizeIO();
   cout.flush(); cerr.flush();
   // Give back memory
   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      for(int side=0; side < 6; side++)
	 if( BCForcing[g][side] != NULL )
	    delete[] BCForcing[g][side];
      delete[] BCForcing[g];
   }
   for( int s = 0 ; s < point_sources.size(); s++ )
      delete point_sources[s];

// why is this barrier needed???
   MPI_Barrier(MPI_COMM_WORLD);
} // end EW::solve2()


