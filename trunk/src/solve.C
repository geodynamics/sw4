#include "EW.h"
#include "impose_cartesian_bc.h"
//#include "impose_curvilinear_bc.h"
#include "F77_FUNC.h"

extern "C" {
void F77_FUNC(dgesv,DGESV)(int & N1, int & N2, double *A, int & LDA, int *IPIV, double * B, int &N3, int &INFO);
void F77_FUNC(factorizeinterfacematrices,FACTORIZEINTERFACEMATRICES)( int*, int*, int*, int*, int*,
								      double*, double*, double*,
								      int*, int*, int*, int*, int*, double*,
								      double*, double*, double*, double*, 
								      double*, int*, double*, int* );
void F77_FUNC(acc_bc_free_i,ACC_BC_FREE_I)
  ( int*, int*, int*, double*, double*, double*, double*, double*, int*, int*, double* );

void F77_FUNC(acc_bc_free_j,ACC_BC_FREE_J)
  ( int*, int*, int*, double*, double*, double*, double*, double*, int*, int*, double* );

void F77_FUNC(acc_bc_free_k,ACC_BC_FREE_K)
  ( int*, int*, int*, double*, double*, double*, double*, double*, int*, int*, double* );

void F77_FUNC(bcfort, BCFORT)( int*, int*, int*, int*, int*, int*, 
			       int *, int*, int*, int*,
			       double*, double*, boundaryConditionType*, double *, double*, double*, double*,
			       double* bf0_p, double* bf1_p, 
			       double* bf2_p, double*bf3_p, 
			       double*bf4_p, double*bf5_p, 
			       double*, double*, double* );
void F77_FUNC(twfrsurfz, TWFRSURFZ)(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p, int * kfirst_p, 
				  int * klast_p, int * nx_p, int * ny_p, int * nz_p, double* h_p, int * k_p,
				  double* t_p, double* om_p, double* cv_p, double* ph_p,
				  double* bforce_side5_ptr, double* mu_ptr, double* la_ptr );
void F77_FUNC(twdirbdry,TWDIRBDRY)( int *wind_ptr, double *h_p, double *t_p, double *om_p, double * cv_p, 
				    double *ph_p,  double * bforce_side_ptr );
}

//--------------------------------------------------------------------
void EW::solve( vector<Source*> & a_GlobalUniqueSources )
{
// solution arrays
  vector<Sarray> F, Lu, Uacc, Up, Um, U;
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
   
// the Source objects get discretized into GridPointSource objects
  vector<GridPointSource*> point_sources;

// Transfer source terms to each individual grid as point sources at grid points.
     for( unsigned int i=0 ; i < a_GlobalUniqueSources.size() ; i++ )
       if (!a_GlobalUniqueSources[i]->ignore())
	 a_GlobalUniqueSources[i]->set_grid_point_sources( this, point_sources );

  if( mVerbose && proc_zero() )
  {
    cout << endl << "***  Starting solve ***" << endl;
  }
  printPreamble();

// this is the time level
  double t=mTstart;
  
// Set up timers
  double time_start_solve = MPI_Wtime();
  double time_measure[7];
  double time_sum[7]={0,0,0,0,0,0,0};
  double bc_time_measure[5]={0,0,0,0,0};

  int beginCycle = 1; // also set in setupRun(), perhaps make member variable?

// Assign initial data
  int g; 
  if (m_twilight_forcing)
  {
    exactSolTwilight(mTstart, U, AlphaVE);
    exactSolTwilight(mTstart-mDt, Um, AlphaVEm);
  }
  else
  {
// homogeneous initial data is the default
    for( g=0 ; g<mNumberOfGrids; g++ )
    {
      U[g].set_to_zero();
      Um[g].set_to_zero();
      for( int a=0 ; a < m_number_mechanisms ; a++ )
      {
	AlphaVEm[g][a].set_to_zero();
	AlphaVE[g][a].set_to_zero();
      }
    }
  }
  
  if ( mVerbose && proc_zero() )
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
  
// // save any images for cycle = 0 (initial data) ?
//   mTime += mDt; // need to increase time by dt to get the initial velocity computation right (uses Up and Um)
//   update_images(beginCycle-1);
//   mTime -= mDt; // reset mTime

// do some testing...
  // if ( proc_zero() )
  //   cout << "***Testing..." << endl;
// tmp
//   for(int g=0; g<mNumberOfGrids; g++)
//   {
//     printf("proc=%i, Onesided[grid=%i]:", m_myRank, g);
//     for (int q=0; q<6; q++)
//       printf(" os[%i]=%i", q, m_onesided[g][q]);
//     printf("\n");
//     printf("proc=%i, bcType[grid=%i]:", m_myRank, g);
//     for (int q=0; q<6; q++)
//       printf(" bc[%i]=%i", q, m_bcType[g][q]);
//     printf("\n");
//   }

// test accuracy of spatial approximation
//   if ( proc_zero() )
//     printf("\n Testing the accuracy of the spatial difference approximation\n");
//   exactRhsTwilight(t, F);
//   evalRHS( U, Up ); // save Lu in composite grid 'Up'
// // evaluate and print errors
//   double lowZ[3], interiorZ[3], highZ[3];
//   bndryInteriorDifference( F, Up, lowZ, interiorZ, highZ );
//   if ( proc_zero() )
//   {
//     printf("Max errors low-k boundary RHS:  %15.7e  %15.7e  %15.7e\n", lowZ[0], lowZ[1], lowZ[2]);
//     printf("Max errors interior RHS:        %15.7e  %15.7e  %15.7e\n", interiorZ[0], interiorZ[1], interiorZ[2]);
//     printf("Max errors high-k boundary RHS: %15.7e  %15.7e  %15.7e\n", highZ[0], highZ[1], highZ[2]);
//   }
  
// // c test accuracy of forcing
//   evalRHS( U, Lu ); // save Lu in composite grid 'Lu'
//   exactForceTwilight( t, F );
//   exactAccTwilight( t, Uacc ); // save Utt in Uacc
//   test_RhoUtt_Lu( Uacc, Lu, F, lowZ, interiorZ, highZ );
//   if ( proc_zero() )
//   {
//     printf("Testing accuracy of rho*utt - L(u) = F\n");
//     printf("Max errors low-k boundary RHS:  %15.7e  %15.7e  %15.7e\n", lowZ[0], lowZ[1], lowZ[2]);
//     printf("Max errors interior RHS:        %15.7e  %15.7e  %15.7e\n", interiorZ[0], interiorZ[1], interiorZ[2]);
//     printf("Max errors high-k boundary RHS: %15.7e  %15.7e  %15.7e\n", highZ[0], highZ[1], highZ[2]);
//   }

  if ( proc_zero() )
  {
    printf("About to call enforceBC on intial data\n");
  }
  
// enforce bc on initial data
// U
// communicate across processor boundaries
  for(int g=0 ; g < mNumberOfGrids ; g++ )
    communicate_array( U[g], g );
// boundary forcing
  cartesian_bc_forcing( t, BCForcing );
// enforce boundary condition
  enforceBC( U, t, BCForcing );   

// Um
// communicate across processor boundaries
  for(int g=0 ; g < mNumberOfGrids ; g++ )
    communicate_array( Um[g], g );
// boundary forcing
  cartesian_bc_forcing( t-mDt, BCForcing );
// enforce boundary condition
  enforceBC( Um, t-mDt, BCForcing );

  if ( proc_zero() )
  {
    printf("About to call exactSolTwilight\n");
  }

// more testing
// check the accuracy of the initial data, store exact solution in Up, ignore AlphaVE
  // exactSolTwilight( t, Up, AlphaVE);

  // double errInf, errL2;
  // normOfDifferenceGhostPoints( Up, U, errInf, errL2 );
  // if ( proc_zero() )
  //   printf("\n Ghost point errors: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);

  if ( proc_zero() )
    cout << "  Begin time stepping..." << endl;

// // save initial data on sac files
//   update_SACs( beginCycle-1 );

// // Begin time stepping loop
  for( int currentTimeStep = beginCycle; currentTimeStep <= mNumberOfTimeSteps; currentTimeStep++)
  {    
    time_measure[0] = MPI_Wtime();

//     if( m_forcing->use_input_sources() )
//     {
// // tmp
// //      printf("There are %i point sources on proc # %i\n", point_sources.size(), m_myRank);

// // assign interior forcing
//         for( int g=0 ; g<mNumberOfGrids; g++ )
//           F[g].set_to_zero();

//        for( int s= 0 ; s < point_sources.size() ; s++ )
//        {
// 	  //	  double rhoi = 1/mRho(sN,sM,sL);
//           double fxyz[3];
// 	  point_sources[s]->getFxyz(mTime,fxyz);
//           int g = point_sources[s]->m_grid;
//           double rhoi = 1.0/(mRho[g].c_ptr()[ind[s]]);
//           double* fp = F[g].c_ptr();
// // tmp
// //          printf("Source # %i proc # %i has index %i %i %i on grid %i, value = (%e, %e, %e)\n", s, m_myRank,
// //                    point_sources[s]->m_i0, point_sources[s]->m_j0, point_sources[s]->m_k0, g,
// //                    fxyz[0]*rhoi, fxyz[1]*rhoi, fxyz[2]*rhoi );

// // the treatment of rho used to be inconsistent between the Cartesian and curvilinear inner loop!
// // in the cartesian inner loop, the forcing is not divided by rho, while 
// // it is (was) divided by rho in the curvilinear inner loop
//           fp[3*ind[s]  ] += fxyz[0]*rhoi;
//           fp[3*ind[s]+1] += fxyz[1]*rhoi;
//           fp[3*ind[s]+2] += fxyz[2]*rhoi;
//        }
//     }
//     else
//     {

// twilight forcing...
    if (m_twilight_forcing)
    {
      exactForceTwilight( t, F );
    }
    
// evaluate right hand side
    evalRHS( U, Lu ); // save Lu in composite grid 'Lu'

// take predictor step, store in Up
    evalPredictor(Up, U, Um, Lu, F );    

    time_measure[1] = MPI_Wtime();
    time_measure[2] = MPI_Wtime();

// communicate across processor boundaries
    for(int g=0 ; g < mNumberOfGrids ; g++ )
      communicate_array( Up[g], g );
// calculate boundary forcing at time t+mDt
    cartesian_bc_forcing( t+mDt, BCForcing );
// update ghost points in Up
    enforceBC( Up, t+mDt, BCForcing );

//     //    check_consintp( Up[0], Up[1], AlphaVEp[0], AlphaVEp[1] );
    time_measure[3] = time_measure[4] = MPI_Wtime();

// get 4th order in time
    if (mOrder == 4)
    {
      if (m_twilight_forcing)
      {
	exactForce_ttTwilight( t, F );
      }
      
      evalDpDmInTime( Up, U, Um, Uacc ); // store result in Uacc
      
      evalRHS( Uacc, Lu );

      evalCorrector( Up, Lu, F );
      
      time_measure[4] = MPI_Wtime();

// also check out EW::update_all_boundaries 
// communicate across processor boundaries
      for(int g=0 ; g < mNumberOfGrids ; g++ )
	communicate_array( Up[g], g );
// calculate boundary forcing at time t+mDt
      cartesian_bc_forcing( t+mDt, BCForcing );
// update ghost points in Up
      enforceBC( Up, t+mDt, BCForcing );
    }
    
// increment time
    t += mDt;

    time_measure[5] = MPI_Wtime();	  

// periodically, print time stepping info to stdout
    printTime( currentTimeStep, t, currentTimeStep == mNumberOfTimeSteps ); 

// // Images have to be written before the solution arrays are cycled, because both Up and Um are needed
// // to compute a centered time derivative
    update_images( currentTimeStep, t, Up );
    
// // Energy evaluation, requires all three time levels present, do before cycle arrays.
//      if( m_energy_log || m_energy_print )
//         compute_energy( mDt, currentTimeStep == mNumberOfTimeSteps );

// cycle the solution arrays
    cycleSolutionArrays(Um, U, Up, AlphaVEm, AlphaVE, AlphaVEp);

// // sac files save the value of U at certain grid points, so need to be called after the solution
// // arrays have been cycled 
//     update_SACs( currentTimeStep );

    time_measure[6] = MPI_Wtime();	  	
	  
// // See if it is time to write a restart file
// //      if (mRestartDumpInterval > 0 &&  currentTimeStep % mRestartDumpInterval == 0)
// //        serialize(currentTimeStep, U, Um);  

    time_sum[0] += time_measure[1]-time_measure[0] + time_measure[4]-time_measure[3]; // step
    time_sum[1] += time_measure[2]-time_measure[1] + time_measure[5]-time_measure[4]; // bcs
    time_sum[2] += time_measure[6]-time_measure[5]; // image & sac
    time_sum[3] += time_measure[6]-time_measure[0];//  total
// where does bc_time_measure come from?
    // time_sum[4] += bc_time_measure[1]-bc_time_measure[0]; // comm. + overlap grids.
    // time_sum[5] += bc_time_measure[2]-bc_time_measure[1]+bc_time_measure[4]-bc_time_measure[3]; // update proc boundary
    // time_sum[6] += bc_time_measure[3]-bc_time_measure[2]; // impose bc.

  } // end time stepping loop

  if ( proc_zero() )
    cout << "  Time stepping finished..." << endl;

//   delete[] wk;

//   if( ind != 0 )
//      delete[] ind;

   double time_end_solve = MPI_Wtime();
   print_execution_time( time_start_solve, time_end_solve, "solver phase" );

   if( m_output_detailed_timing )
     print_execution_times( time_sum );

   if (m_twilight_forcing)
   {
// check the accuracy of the final solution, store exact solution in Up, ignore AlphaVE
     exactSolTwilight( t, Up, AlphaVE);

     double errInf, errL2;
     normOfDifference( Up, U, errInf, errL2 );
     if ( proc_zero() )
       printf("\n Final solution errors: Linf = %15.7e, L2 = %15.7e\n", errInf, errL2);
   }

    finalizeIO();
    cout.flush(); cerr.flush();

// why is this barrier needed???
    MPI_Barrier(MPI_COMM_WORLD);

//   if( m_forcing->knows_exact() )
//      computeSolutionError(U, mTime, AlphaVE ); // note that final solution ends up in U after the call to cycleSolutionArrays()

} // end EW::solve()

//------------------------------------------------------------------------
void EW::cycleSolutionArrays(vector<Sarray> & a_Um, vector<Sarray> & a_U, vector<Sarray> & a_Up, 
			     vector<Sarray*> & a_AlphaVEm, vector<Sarray*> & a_AlphaVE, vector<Sarray*> & a_AlphaVEp)
{
  for (int g=0; g<mNumberOfGrids; g++)
  {
    double *tmp = a_Um[g].c_ptr();
    a_Um[g].reference(a_U[g].c_ptr());
    a_U[g].reference(a_Up[g].c_ptr());
    a_Up[g].reference(tmp);
    for( int a = 0 ; a < m_number_mechanisms ; a++ )
    {
       double *tmp = a_AlphaVEm[g][a].c_ptr();
       a_AlphaVEm[g][a].reference(a_AlphaVE[g][a].c_ptr());
       a_AlphaVE[g][a].reference(a_AlphaVEp[g][a].c_ptr());
       a_AlphaVEp[g][a].reference(tmp);
    }
  }
}

//---------------------------------------------------------------------------
void EW::enforceBC( vector<Sarray> & a_U, double t, vector<double **> & a_BCForcing )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz;
  double *u_ptr, *mu_ptr, *la_ptr, h;
  boundaryConditionType *bcType_ptr;
  double *bforce_side0_ptr, *bforce_side1_ptr, *bforce_side2_ptr, *bforce_side3_ptr, *bforce_side4_ptr, *bforce_side5_ptr;
  int *wind_ptr;
  double om=0, ph=0, cv=0;
    
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    u_ptr    = a_U[g].c_ptr();
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
    bcType_ptr = m_bcType[g]; // is this right???
    
    wind_ptr = m_BndryWindow[g];
    
// THESE ARRAYS MUST BE FILLED IN BEFORE CALLING THIS ROUTINE
    bforce_side0_ptr = a_BCForcing[g][0]; // low-i bndry forcing array pointer
    bforce_side1_ptr = a_BCForcing[g][1]; // high-i bndry forcing array pointer
    bforce_side2_ptr = a_BCForcing[g][2]; // low-j bndry forcing array pointer
    bforce_side3_ptr = a_BCForcing[g][3]; // high-j bndry forcing array pointer
    bforce_side4_ptr = a_BCForcing[g][4]; // low-k bndry forcing array pointer
    bforce_side5_ptr = a_BCForcing[g][5]; // high-k bndry forcing array pointer
    
    F77_FUNC(bcfort, BCFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
			      wind_ptr, &nx, &ny, &nz,
			      u_ptr, &h, bcType_ptr, m_sbop, mu_ptr, la_ptr, &t,
			      bforce_side0_ptr, bforce_side1_ptr, 
			      bforce_side2_ptr, bforce_side3_ptr, 
			      bforce_side4_ptr, bforce_side5_ptr, 
			      &om, &ph, &cv );
    

      
  }
}

//-----------------------------------------------------------------------
// void EW::update_all_boundaries(vector<Sarray> &U, vector<Sarray> &UM, double t, vector<Sarray*> &AlphaVE )
// {
//   int g;
// boundary forcing on Cartesian grids
//   for (g=0; g<mNumberOfCartesianGrids; g++)
//   {
//     cartesian_bc_forcing(mMu[g], m_BCForcing[g], m_NumberOfBCPoints[g], m_bcType[g], m_forcing, m_zmin[g], 
// 			 t, mGridSize[g], m_onesided[g], m_curlcoeff, m_number_mechanisms, AlphaVE[g],
// 			 mMuVE[g], mLambdaVE[g] );
//   }
// // boundary forcing for curvilinear grid
//   if (topographyExists())
//   {
//     g = mNumberOfGrids - 1;
//     curvilinear_bc_forcing(m_BCForcing[g], m_NumberOfBCPoints[g],
// 			   mX, mY, mZ,
// 			   mMu[g], mLambda[g], mQ, mR, mS, mJ,
// 			   mTime, m_bcType[g], m_forcing, m_number_mechanisms, AlphaVE[g],
// 			   mMuVE[g], mLambdaVE[g] );
// tmp
// save the bc forcing on side=4 (assume one proc)
//
//     if (m_bcType[g][4] == bStressFree)
//     {
//       FILE *tf=fopen("bstress.txt","w");
      
//       int ind=0;
//       fprintf(tf,"%i %i\n", mX.m_ni-2, mX.m_nj-2);
//       for (int j=2; j<=mX.m_nj-1; j++)
// 	for (int i=2; i<=mX.m_ni-1; i++)
// 	{
// 	  fprintf(tf,"%e %e %e\n", m_BCForcing[g][4][ind], m_BCForcing[g][4][ind+1], m_BCForcing[g][4][ind+2]);
// 	  ind+=3;
// 	}

//       fclose(tf);
//     }
// end tmp
//  } // end if topographyExists()
  
// boundary forcing (NEW style)
//   cartesian_bc_forcing( t );
    
//   if( mNumberOfGrids > 1 )
//   {
//     for( g=0 ; g < mNumberOfGrids ; g++ )
//        communicate_array( U[g], g );
// // Why is it necessary to impose the physical bc twice?   
//     enforceBC( U, t );
// //    impose_physical_bc( U, UM, t, AlphaVE );
// // interpolate in vertical direction
// // not yet activated
// //    interpolate_between_grids( U, UM, t, AlphaVE );
//   }

// // communicate across processor boundaries
//   for( g=0 ; g < mNumberOfGrids ; g++ )
//     communicate_array( U[g], g );
// // physical bc's on the outer boundaries
//   enforceBC( U, t );
// //  impose_physical_bc( U, UM, t, AlphaVE );
// // communicate across processor boundaries
//   for( g=0 ; g < mNumberOfGrids ; g++ )
//     communicate_array( U[g], g );

// tmp
// testing eval_curvilinear_bc_stress
//   if (topographyExists())
//   {
//     g = mNumberOfGrids - 1;
//     eval_curvilinear_bc_stress(U[g], m_BCForcing[g], mX, mY, mZ, mMu[g], mLambda[g],
// 			       mQ, mR, mS, mJ);

// // save the bc forcing on side=4 (assume one proc)
// //
//     if (m_bcType[g][4] == bStressFree)
//     {
//       FILE *tf=fopen("estress.txt","w");
      
//       int ind=0;
//       fprintf(tf,"%i %i\n", mX.m_ni-2, mX.m_nj-2);
//       for (int j=2; j<=mX.m_nj-1; j++)
// 	for (int i=2; i<=mX.m_ni-1; i++)
// 	{
// 	  fprintf(tf,"%e %e %e\n", m_BCForcing[g][4][ind], m_BCForcing[g][4][ind+1], m_BCForcing[g][4][ind+2]);
// 	  ind+=3;
// 	}

//       fclose(tf);
//       MPI_Abort(MPI_COMM_WORLD, 1);
//     }
//   }
// end tmp  
//}

//-----------------------------------------------------------------------
// void EW::impose_physical_bc(vector<Sarray> &U, vector<Sarray> &UM, double t, 
// 			    vector<Sarray*> &AlphaVE )
// {
//   int g;

//   // if( m_do_geodynbc )
//   //    impose_geodyn_ibcdata( U, UM, t );

//   for( g = 0 ; g < mNumberOfCartesianGrids ; g++ ) // only Cartesian grids
// //  for( g = 0 ; g < mNumberOfGrids ; g++ ) // tmp: all grids
//   {
//     double h= mGridSize[g];
// // Find out boundary condition types, and set some straightforward bc:s
//     for( int side= 0 ; side < 6 ; side++ )
//     {
// // new approach: assign all sides with dirichlet and free surface bc in one go.
//       impose_cartesian_bc( U[g], mMu[g], mLambda[g], m_BCForcing[g], m_bcType[g], // note that U points to the same array as a_Up
// 			   m_onesided[g], m_curlcoeff, h ); 

//     } // end for all Cartesian grids

//     if (topographyExists())
//     {
// //     if (proc_zero())
// //       cout << "*** Calling impose_curvilinear_bc ***" << endl;
    
//       g = mNumberOfGrids-1;
//       impose_curvilinear_bc(U[g], m_BCForcing[g],
// 			    mX, mY, mZ, mMu[g], mLambda[g],
// 			    mQ, mR, mS, mJ, t, m_bcType[g] ); 

//     }
  
//   }
// }



// //-----------------------------------------------------------------------
// void EW::bc_dirichlet( Sarray& u, int grid, double t, int side,
// 			 Forcing* forcing, double h )
// {
// // this routine sets Dirichlet conditions on the boundary points for Cartesian grids (and not the ghost points!)
//    int wind[6];
//    u.side_plane( side, wind );
// // tmp	
// //   printf("bc_dirichlet side=%i, wind=(%i %i %i %i %i %i)\n", side, wind[0], wind[1], wind[2], wind[3], wind[4], wind[5]);

//    int sidecorr[3] = {0,0,0};
//    if( side == 0 )
//       sidecorr[0] = 1;
//    else if( side == 1 )
//       sidecorr[0] = -1;
//    else if( side == 2 )
//       sidecorr[1] = 1;
//    else if( side == 3 )
//       sidecorr[1] =-1;
//    else if( side == 4 )
//       sidecorr[2] = 1;
//    else if( side == 5 )
//       sidecorr[2] =-1;

//    for( int k = wind[4] ; k <= wind[5] ; k++ )
//       for( int j = wind[2] ; j <= wind[3] ; j++ )
// 	 for( int i = wind[0] ; i <= wind[1] ; i++ )
// 	 {
//             double x, y, z;
//             if( grid < mNumberOfCartesianGrids )
// 	    {
// 	       x = (i+sidecorr[0]-1)*h;
//                y = (j+sidecorr[1]-1)*h;
// 	       z = (k+sidecorr[2]-1)*h + m_zmin[grid];
// //                x = (i-1)*h;
// //                y = (j-1)*h;
// // 	          z = (k-1)*h + m_zmin[grid];
// 	    }
// 	    else
// 	    {
// 	       x = mX(i+sidecorr[0], j+sidecorr[1], k+sidecorr[2]);
//                y = mY(i+sidecorr[0], j+sidecorr[1], k+sidecorr[2]);
// 	       z = mZ(i+sidecorr[0], j+sidecorr[1], k+sidecorr[2]);
// 	    }
// // assign the boundary value straight into the solution array
// 	    forcing->get_exact( x, y, z, t, &(u(1,i+sidecorr[0], j+sidecorr[1], k+sidecorr[2])), h );
// 	 }
// }

//------------------------------------------------------------------------------
void EW::cartesian_bc_forcing(double t, vector<double **> & a_BCForcing )
// assign the boundary forcing arrays a_BCForcing[g][side]
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz;
  double *u_ptr, *mu_ptr, *la_ptr, h;
  boundaryConditionType *bcType_ptr;
  double *bforce_side0_ptr, *bforce_side1_ptr, *bforce_side2_ptr, *bforce_side3_ptr, *bforce_side4_ptr, *bforce_side5_ptr;
  int *wind_ptr;
  double om=0, ph=0, cv=0;
    
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
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
    
    wind_ptr = m_BndryWindow[g];
    
// pointers to the six sides of the cube
    bforce_side0_ptr = a_BCForcing[g][0]; // low-i bndry forcing array pointer
    bforce_side1_ptr = a_BCForcing[g][1]; // high-i bndry forcing array pointer
    bforce_side2_ptr = a_BCForcing[g][2]; // low-j bndry forcing array pointer
    bforce_side3_ptr = a_BCForcing[g][3]; // high-j bndry forcing array pointer
    bforce_side4_ptr = a_BCForcing[g][4]; // low-k bndry forcing array pointer
    bforce_side5_ptr = a_BCForcing[g][5]; // high-k bndry forcing array pointer

    if (m_twilight_forcing)
    {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;

// the following code can probably be improved by introducing a loop over all sides,
// but bStressFree is only implemented for side=4 and 5, so there must be some special cases
      int k = 1;
      if (m_bcType[g][0] == bDirichlet)
      {
	F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[0], &h, &t, &om, &cv, &ph, bforce_side0_ptr );
//              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce6 )
      }

      if (m_bcType[g][1] == bDirichlet)
      {
	F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6], &h, &t, &om, &cv, &ph, bforce_side1_ptr );
      }

      if (m_bcType[g][2] == bDirichlet)
      {
	F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6*2], &h, &t, &om, &cv, &ph, bforce_side2_ptr );
      }

      if (m_bcType[g][3] == bDirichlet)
      {
	F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6*3], &h, &t, &om, &cv, &ph, bforce_side3_ptr );
      }

      if (m_bcType[g][4] == bDirichlet)
      {
	F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6*4], &h, &t, &om, &cv, &ph, bforce_side4_ptr );
      }
      else if (m_bcType[g][4] == bStressFree)
      {
	 k = 1;
	 F77_FUNC(twfrsurfz, TWFRSURFZ)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				       &klast, &nx, &ny, &nz, &h, &k, &t, &om, &cv, &ph,
				       bforce_side4_ptr, mu_ptr, la_ptr );
      }

      if (m_bcType[g][5] == bDirichlet)
      {
	F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6*5], &h, &t, &om, &cv, &ph, bforce_side5_ptr );
      }
      else if (m_bcType[g][5] == bStressFree)
      {
	 k = nz;
	 F77_FUNC(twfrsurfz, TWFRSURFZ)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				       &klast, &nx, &ny, &nz, &h, &k, &t, &om, &cv, &ph,
				       bforce_side5_ptr, mu_ptr, la_ptr );
      }
      
//               call TWFRSURFZ( ifirst, ilast, jfirst, jlast, kfirst, 
//      +             klast, nx, ny, nz, h, k, t, om, cv, ph, bforce5, mu, 
//      +             la )
    }
    else
    {
// no boundary forcing
// we can do the same loop for all types of bc. For bParallel boundaries, numberOfBCPoints=0
      int q;
      for (q=0; q<m_NumberOfBCPoints[g][0]; q++)
	bforce_side0_ptr[q] = 0.;
      for (q=0; q<m_NumberOfBCPoints[g][1]; q++)
	bforce_side1_ptr[q] = 0.;
      for (q=0; q<m_NumberOfBCPoints[g][2]; q++)
	bforce_side2_ptr[q] = 0.;
      for (q=0; q<m_NumberOfBCPoints[g][3]; q++)
	bforce_side3_ptr[q] = 0.;
      for (q=0; q<m_NumberOfBCPoints[g][4]; q++)
	bforce_side4_ptr[q] = 0.;
      for (q=0; q<m_NumberOfBCPoints[g][5]; q++)
	bforce_side5_ptr[q] = 0.;
    }
  }
}



//------------------------------------------------------------------------------
// void
// impose_curvilinear_bc(Sarray & a_u, double ** bcForcing, Sarray & a_x, Sarray & a_y, Sarray & a_z,
// 		      Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J,
// 		      double t, boundaryConditionType bcType[] )

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
   
// //
// // This routine only knows about dirichlet and free surface conditions
// //

//   side = 2;
//   ind = 0;
//   if (bcType[2]==bDirichlet || bcType[2]==bSuperGrid)
//   {
// // set the dirichlet condition on the boundary point
//     int j=2;
//     for(int k=1; k<=Nz; k++)
//       for(int i=1; i<=Nx; i++)
//       {
// 	for(int c=1; c<=3;c++)
// 	  u(i,j,k,c) = bcForcing[side][ind+c-1];
// 	ind += 3;
//       }
//   }

//   side = 3;
//   ind = 0;
//   if (bcType[3]==bDirichlet || bcType[3]==bSuperGrid)
//   {
// // set the dirichlet condition on the boundary point
//     j=Ny-1;
//     for(int k=1; k<=Nz; k++)
//       for(int i=1; i<=Nx; i++)
//       {
// 	for(int c=1; c<=3;c++)
// 	  u(i,j,k,c) = bcForcing[side][ind+c-1];
// 	ind += 3;
//       }
//   }

//   side = 0;
//   ind = 0;
//   if (bcType[0]==bDirichlet || bcType[0]==bSuperGrid)
//   {
// // set the dirichlet condition on the boundary point
//     i=2;
//     for(int k=1; k<=Nz; k++)
//       for(int j=1; j<=Ny; j++)
//       {
// 	for(int c=1; c<=3;c++)
// 	  u(i,j,k,c) = bcForcing[side][ind+c-1];
// 	ind += 3;
//       }
//   }
  
//   side = 1;
//   ind = 0;
//   if (bcType[1]==bDirichlet || bcType[1]==bSuperGrid)
//   {
// // set the dirichlet condition on the boundary point
//     i=Nx-1;
//     for(int k=1; k<=Nz; k++)
//       for(int j=1; j<=Ny; j++)
//       {
// 	for(int c=1; c<=3;c++)
// 	  u(i,j,k,c) = bcForcing[side][ind+c-1];
// 	ind += 3;
//       }
//   }

//   side = 5;
//   ind = 0;
//   if (bcType[5]==bDirichlet || bcType[5]==bSuperGrid)
//   {
// // set the dirichlet condition on the boundary point
//     k=Nz-1;
//     for(int j=1; j<=Ny; j++)
//       for(int i=1; i<=Nx; i++)
//       {
// 	for(int c=1; c<=3;c++)
// 	  u(i,j,k,c) = bcForcing[side][ind+c-1];
// 	ind += 3;
//       }
//   }
//   else if (bcType[5]==bStressFree)
//   {
// #define E1(i,j,t1,t2,i2,t3,i3,t4) (t1(i,j,Nz-1)*t2(i,j,Nz-1,i2)*t3(i,j,Nz-1,i3)*t4(i,j,Nz-1))
// #define E12(i,j,t1,t2,i2,t3,i3,t4) (0.5*(t1(i,j,Nz)*t2(i,j,Nz,i2)*t3(i,j,Nz,i3)*t4(i,j,Nz)\
// 				 +t1(i,j,Nz-1)*t2(i,j,Nz-1,i2)*t3(i,j,Nz-1,i3)*t4(i,j,Nz-1)))
// #define E32(i,j,t1,t2,i2,t3,i3,t4) (0.5*(t1(i,j,Nz-2)*t2(i,j,Nz-2,i2)*t3(i,j,Nz-2,i3)*t4(i,j,Nz-2)\
// 				 +t1(i,j,Nz-1)*t2(i,j,Nz-1,i2)*t3(i,j,Nz-1,i3)*t4(i,j,Nz-1)))
//     double Dr0[3],Dq0[3],Dsp[3],ugp[3],b[3],A[9];
//     int info,N1,N2,N3 ,ipv[3]; 
//     N1=3;N2=1;N3=3;
//     for(int j=2; j<=Ny-1; j++)
//       for(int i=2; i<=Nx-1; i++)
//       {
// // note that the free surface conditions are imposed on the boundary itself (as opposed to the ghost point)
// 	for(int c=1; c<=3;c++){
// 	  Dr0[c-1]=0.5*(u(i,j+1,Nz-1,c)-u(i,j-1,Nz-1,c));
// 	  Dq0[c-1]=0.5*(u(i+1,j,Nz-1,c)-u(i-1,j,Nz-1,c));
// 	  Dsp[c-1]=u(i,j,Nz-1,c)-u(i,j,Nz-2,c);
// 	}
// 	b[0]=(E1(i,j,J,s,1,q,1,lam2mu)*Dq0[0]
// 	      +E1(i,j,J,s,1,r,1,lam2mu)*Dr0[0]	
// 	      +E1(i,j,J,s,1,q,2,lam)*Dq0[1]
// 	      +E1(i,j,J,s,1,r,2,lam)*Dr0[1]	 
// 	      +E1(i,j,J,s,1,q,3,lam)*Dq0[2]
// 	      +E1(i,j,J,s,1,r,3,lam)*Dr0[2]	 
// 	      +E1(i,j,J,s,2,q,1,mu)*Dq0[1]
// 	      +E1(i,j,J,s,2,r,1,mu)*Dr0[1]	 
// 	      +E1(i,j,J,s,2,q,2,mu)*Dq0[0]
// 	      +E1(i,j,J,s,2,r,2,mu)*Dr0[0]	 
// 	      +E1(i,j,J,s,3,q,1,mu)*Dq0[2]
// 	      +E1(i,j,J,s,3,r,1,mu)*Dr0[2]	 
// 	      +E1(i,j,J,s,3,q,3,mu)*Dq0[0]
// 	      +E1(i,j,J,s,3,r,3,mu)*Dr0[0]	 
// 	      +0.5*(
// 		E32(i,j,J,s,1,s,1,lam2mu)*Dsp[0]
// 		+E32(i,j,J,s,1,s,2,lam)*Dsp[1]
// 		+E32(i,j,J,s,1,s,3,lam)*Dsp[2]
// 		+E32(i,j,J,s,2,s,1,mu)*Dsp[1]
// 		+E32(i,j,J,s,2,s,2,mu)*Dsp[0]
// 		+E32(i,j,J,s,3,s,1,mu)*Dsp[2]
// 		+E32(i,j,J,s,3,s,3,mu)*Dsp[0])); 
// 	b[1]=(E1(i,j,J,s,1,q,1,mu)*Dq0[1]
// 	      +E1(i,j,J,s,1,r,1,mu)*Dr0[1]	
// 	      +E1(i,j,J,s,1,q,2,mu)*Dq0[0]
// 	      +E1(i,j,J,s,1,r,2,mu)*Dr0[0]	 
// 	      +E1(i,j,J,s,2,q,2,lam2mu)*Dq0[1]
// 	      +E1(i,j,J,s,2,r,2,lam2mu)*Dr0[1]
// 	      +E1(i,j,J,s,2,r,1,lam)*Dr0[0]	 
// 	      +E1(i,j,J,s,2,q,1,lam)*Dq0[0]
// 	      +E1(i,j,J,s,2,r,3,lam)*Dr0[2]	 
// 	      +E1(i,j,J,s,2,q,3,lam)*Dq0[2]
// 	      +E1(i,j,J,s,3,q,2,mu)*Dq0[2]
// 	      +E1(i,j,J,s,3,r,2,mu)*Dr0[2]	 
// 	      +E1(i,j,J,s,3,q,3,mu)*Dq0[1]
// 	      +E1(i,j,J,s,3,r,3,mu)*Dr0[1]	 
// 	      +0.5*(E32(i,j,J,s,2,s,2,lam2mu)*Dsp[1]
// 		    +E32(i,j,J,s,2,s,1,lam)*Dsp[0]
// 		    +E32(i,j,J,s,2,s,3,lam)*Dsp[2]
// 		    +E32(i,j,J,s,1,s,1,mu)*Dsp[1]
// 		    +E32(i,j,J,s,1,s,2,mu)*Dsp[0]
// 		    +E32(i,j,J,s,3,s,2,mu)*Dsp[2]
// 		    +E32(i,j,J,s,3,s,3,mu)*Dsp[1]));
	    
// 	b[2]=(E1(i,j,J,s,1,q,1,mu)*Dq0[2]
// 	      +E1(i,j,J,s,1,r,1,mu)*Dr0[2]	
// 	      +E1(i,j,J,s,1,q,3,mu)*Dq0[0]
// 	      +E1(i,j,J,s,1,r,3,mu)*Dr0[0]	 
// 	      +E1(i,j,J,s,3,q,3,lam2mu)*Dq0[2]
// 	      +E1(i,j,J,s,3,r,3,lam2mu)*Dr0[2]
// 	      +E1(i,j,J,s,3,r,1,lam)*Dr0[0]	 
// 	      +E1(i,j,J,s,3,q,1,lam)*Dq0[0]
// 	      +E1(i,j,J,s,3,r,2,lam)*Dr0[1]	 
// 	      +E1(i,j,J,s,3,q,2,lam)*Dq0[1]
// 	      +E1(i,j,J,s,2,q,2,mu)*Dq0[2]
// 	      +E1(i,j,J,s,2,r,2,mu)*Dr0[2]	 
// 	      +E1(i,j,J,s,2,q,3,mu)*Dq0[1]
// 	      +E1(i,j,J,s,2,r,3,mu)*Dr0[1]	 
// 	      +0.5*(E32(i,j,J,s,3,s,3,lam2mu)*Dsp[2]
// 		    +E32(i,j,J,s,3,s,1,lam)*Dsp[0]
// 		    +E32(i,j,J,s,3,s,2,lam)*Dsp[1]
// 		    +E32(i,j,J,s,1,s,1,mu)*Dsp[2]
// 		    +E32(i,j,J,s,1,s,3,mu)*Dsp[0]
// 		    +E32(i,j,J,s,2,s,2,mu)*Dsp[2]
// 		    +E32(i,j,J,s,2,s,3,mu)*Dsp[1]));

// // add in boundary forcing
// 	for(int c=1; c<=3;c++)
// 	{
// 	  b[c-1]+= bcForcing[side][ind+c-1];
// 	}
// 	ind += 3;
	    
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
// // u, v, w in w eq.
// 	A[2]=(+E12(i,j,J,s,3,s,1,lam)
// 	      +E12(i,j,J,s,1,s,3,mu));
// 	A[5]=(E12(i,j,J,s,2,s,3,mu)
// 	      +E12(i,j,J,s,3,s,2,lam));
// 	A[8]=(+E12(i,j,J,s,3,s,3,lam2mu)
// 	      +E12(i,j,J,s,1,s,1,mu)
// 	      +E12(i,j,J,s,2,s,2,mu));
// 	for(int c=0; c<9; c++)
// 	  A[c]*=0.5;
	    
// 	F77_FUNC(dgesv,DGESV)(N1, N2, &A[0], N1, &ipv[0], &b[0], N3, info);
// 	for(int c=1; c<=3;c++)
// 	  u(i,j,Nz,c)=u(i,j,Nz-1,c)-b[c-1];
//       }
// #undef E1
// #undef E12
// #undef E32
      
//   }
  
//   side = 4;
//   ind = 0;
//   if (bcType[4]==bDirichlet || bcType[4]==bSuperGrid)
//   {
// // set the dirichlet condition on the boundary point
//     k=2;
//     for(int j=1; j<=Ny; j++)
//       for(int i=1; i<=Nx; i++)
//       {
// 	for(int c=1; c<=3;c++)
// 	  u(i,j,k,c) = bcForcing[side][ind+c-1];
// 	ind += 3;
//       }
//   }

//   if (bcType[4]==bStressFree)
//   {
// #define E1(i,j,t1,t2,i2,t3,i3,t4) (t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2))
// #define E12(i,j,t1,t2,i2,t3,i3,t4) (0.5*(t1(i,j,1)*t2(i,j,1,i2)*t3(i,j,1,i3)*t4(i,j,1)\
// 				 +t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2)))
// #define E32(i,j,t1,t2,i2,t3,i3,t4) (0.5*(t1(i,j,3)*t2(i,j,3,i2)*t3(i,j,3,i3)*t4(i,j,3)\
// 				 +t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2)))
//     double Dr0[3],Dq0[3],Dsp[3],ugp[3],b[3],A[9];
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

// 	for(int c=1; c<=3;c++)
// 	{
// 	  b[c-1]+= bcForcing[side][ind+c-1];
// 	}
// 	ind += 3;

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
	    
// 	F77_FUNC(dgesv,DGESV)(N1, N2 ,&A[0] ,N1 ,&ipv[0] , &b[0] ,N3 ,info);
	    
// // update ghost point
// 	for(int c=1; c<=3;c++)
// 	  u(i,j,1,c)=(u(i,j,2,c)-b[c-1]);
	    
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

//---------------------------------------------------------------
// void curvilinear_bc_forcing(double ** bcForcing,
// 		       int * numberOfBCPoints,
// 		       Sarray & a_x,
// 		       Sarray & a_y,
// 		       Sarray & a_z,
// 		       Sarray & a_mu,
// 		       Sarray & a_lam,
// 		       Sarray & a_q,
// 		       Sarray & a_r,
// 		       Sarray & a_s,
// 		       Sarray & a_J,
// 		       double t,
// 		       boundaryConditionType bcType[],
// 		       Forcing * f,
// 		       int noof_mechanisms, Sarray* alphaVE, Sarray* muVE, Sarray* lambdaVE )
// {
// // 4D macros swap the last and first indices to compensate for different conventions between how 
// // the Sarrays were allocated and how this routine was originally written
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

// // bcType[0-5] are low-q, high-q, low-r, high-r, low-s, high-s
// // this routine only knows about Dirichlet and Stress-free boundary conditions

// // exact solution must be known everywhere in the domain
//   if (f->knows_exact() && !f->exact_only_surface() )
//   {
//     int i, j, k, q, ind, side;
//     double tau[6], uEx[3], dum=0.;
  
//     side = 4;
//     ind = 0;
//  // Dirichlet. k=1 is the ghost point
//     if (bcType[4]==bDirichlet || bcType[4]==bSuperGrid) 
//     {
// // set the dirichlet condition on the boundary point
//       k=2;
//       for(int j=1; j<=Ny; j++)
// 	for(int i=1; i<=Nx; i++)
// 	{
// 	  f->get_exact( x(i,j,k), y(i,j,k), z(i,j,k), t, uEx, dum);	  
// 	  for(int c=1; c<=3;c++)
// 	  {
// 	    bcForcing[side][ind+c-1] = uEx[c-1];
// 	  }
// 	  ind += 3;
// 	}
//     }
// // stress-free at k=2
//     else if (bcType[4]==bStressFree)
//     {
//       for(int j=2; j<=Ny-1; j++)
// 	for(int i=2; i<=Nx-1; i++)
//         {
// 	  f->stress_tensor( x(i,j,2), y(i,j,2), z(i,j,2), t, tau);	  

// 	  bcForcing[side][ind] = J(i,j,2)
// 	    *(s(i,j,2,1)*tau[0]
// 	      +s(i,j,2,2)*tau[1]
// 	      +s(i,j,2,3)*tau[2]
// 	      );
// 	  bcForcing[side][ind+1] = J(i,j,2)
// 	    *(s(i,j,2,2)*tau[3]
// 	      +s(i,j,2,1)*tau[1]
// 	      +s(i,j,2,3)*tau[4]
// 	      );
// 	  bcForcing[side][ind+2] = J(i,j,2)
// 	    *(s(i,j,2,3)*tau[5]
// 	      +s(i,j,2,1)*tau[2]
// 	      +s(i,j,2,2)*tau[4]
// 	      );

// // tmp
// // 	  if (i==50 && j==30)
// // 	  {
// // 	    printf("curvilinear_bc_forcing: t=%e x=%e y=%e z=%e\n"
// // 		   "J=%e, sx=%e, sy=%e, sz=%e\n"
// // 		   "tau[0]=%e, tau[1]=%e, tau[2]=%e, tau[3]=%e, tau[4]=%e, tau[5]=%e\n", t,  x(i,j,2), y(i,j,2), z(i,j,2), 
// // //		   bcForcing[side][ind], bcForcing[side][ind+1], bcForcing[side][ind+2],
// // 		   J(i,j,2), s(i,j,2,1), s(i,j,2,2), s(i,j,2,3),
// // 		   tau[0], tau[1], tau[2], tau[3], tau[4], tau[5]);
// // 	  }

// 	  ind += 3;
// 	}
//     } // end bStressFree
  
//     side = 5;
//     ind = 0;
// // Dirichlet (k=Nz is the ghost point)
//     if (bcType[5]==bDirichlet || bcType[5]==bSuperGrid)
//     {
// // set the dirichlet condition on the boundary point
//       k=Nz-1;
//       for(int j=1; j<=Ny; j++)
// 	for(int i=1; i<=Nx; i++)
// 	{
// 	  f->get_exact( x(i,j,k), y(i,j,k), z(i,j,k), t, uEx, dum);	  
// 	  for(int c=1; c<=3;c++)
// 	  {
// 	    bcForcing[side][ind+c-1] = uEx[c-1];
// 	  }
// 	  ind += 3;
// 	}
//     }

// // Stress-free at k=Nz-1
//     else if (bcType[5]==bStressFree)
//     {
//       for(int j=2; j<=Ny-1; j++)
// 	for(int i=2; i<=Nx-1; i++)
//         {
// 	  f->stress_tensor( x(i,j,Nz-1), y(i,j,Nz-1), z(i,j,Nz-1), t, tau);	  
	  
// 	  bcForcing[side][ind]  = -J(i,j,Nz-1)
// 	    *(s(i,j,Nz-1,1)*tau[0]
// 	      +s(i,j,Nz-1,2)*tau[1]
// 	      +s(i,j,Nz-1,3)*tau[2]
// 	      );
// 	  bcForcing[side][ind+1] = -J(i,j,Nz-1)
// 	    *(s(i,j,Nz-1,2)*tau[3]
// 	      +s(i,j,Nz-1,1)*tau[1]
// 	      +s(i,j,Nz-1,3)*tau[4]
// 	      );
// 	  bcForcing[side][ind+2] = -J(i,j,Nz-1)
// 	    *(s(i,j,Nz-1,3)*tau[5]
// 	      +s(i,j,Nz-1,1)*tau[2]
// 	      +s(i,j,Nz-1,2)*tau[4]
// 	      );
// 	  ind += 3;
// 	}
//     } // end if stress free at k=Nz
  
//     side = 2;
//     ind = 0;
// // Dirichlet (j=1 is the ghost point)
//     if (bcType[2]==bDirichlet || bcType[2]==bSuperGrid)
//     {
// // set the dirichlet condition on the boundary point
//       int j=2;
//       for(int k=1; k<=Nz; k++)
// 	for(int i=1; i<=Nx; i++)
// 	{
// 	  f->get_exact( x(i,j,k), y(i,j,k), z(i,j,k), t, uEx, dum);	  
// 	  for(int c=1; c<=3;c++)
// 	  {
// 	    bcForcing[side][ind+c-1] = uEx[c-1];
// 	  }
// 	  ind += 3;
// 	}
//     }

//     side = 3;
//     ind = 0;
// // Dirichlet (j=Ny is the ghost point)
//     if (bcType[3]==bDirichlet || bcType[3]==bSuperGrid)
//     { 
// // set the dirichlet condition on the boundary point
//       j=Ny-1;
//       for(int k=1; k<=Nz; k++)
// 	for(int i=1; i<=Nx; i++)
// 	{
// 	  f->get_exact( x(i,j,k), y(i,j,k), z(i,j,k), t, uEx, dum);	  
// 	  for(int c=1; c<=3;c++)
// 	  {
// 	    bcForcing[side][ind+c-1] = uEx[c-1];
// 	  }
// 	  ind += 3;
// 	}
//     }

//     side = 0;
//     ind = 0;
// // Dirichlet (i=1 is the ghost point)
//     if (bcType[0]==bDirichlet || bcType[0]==bSuperGrid)
//     { 
// // set the dirichlet condition on the boundary point
//       i=2;
//       for(int k=1; k<=Nz; k++)
// 	for(int j=1; j<=Ny; j++)
// 	{
// 	  f->get_exact( x(i,j,k), y(i,j,k), z(i,j,k), t, uEx, dum);	  
// 	  for(int c=1; c<=3;c++)
// 	  {
// 	    bcForcing[side][ind+c-1] =uEx[c-1];
// 	  }
// 	  ind += 3;
// 	}
//     }

//     side = 1;
//     ind = 0;
// // Dirichlet on i=Nx
//     if (bcType[1]==bDirichlet || bcType[1]==bSuperGrid)
//     {
// // set the dirichlet condition on the boundary point
//       i=Nx-1;
//       for(int k=1; k<=Nz; k++)
// 	for(int j=1; j<=Ny; j++)
// 	{
// 	  f->get_exact( x(i,j,k), y(i,j,k), z(i,j,k), t, uEx, dum);	  
// 	  for(int c=1; c<=3;c++)
// 	  {
// 	    bcForcing[side][ind+c-1] = uEx[c-1];
// 	  }
// 	  ind += 3;
// 	}
    
//     }
//   } // end forcing knows exact solution
//   else // just assign zeros to the forcing array
//   {
//     for (int side = 0; side < 6; side++)
//     {
//       if (bcType[side]==bDirichlet || bcType[side]==bSuperGrid || bcType[side]==bStressFree)
//       {
// 	for (int i=0; i<3*numberOfBCPoints[side]; i++)
// 	  bcForcing[side][i]=0.;
//       }
//     }
//   } // end homogeneous forcing

// // Add attenuation contribution to free surface bcForcing.
//   if (bcType[4]==bStressFree)
//   {
// // Add attenuation contribution to free surface b.c:s
//       for( int a = 0 ; a < noof_mechanisms ; a++ )
// 	eval_curvilinear_bc_stress(alphaVE[a], bcForcing, a_x, a_y, a_z, muVE[a], lambdaVE[a], a_q, a_r, a_s, a_J);
//   }
  
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


//---------------------------------------------------------------------------
void eval_curvilinear_bc_stress(Sarray & a_u, double ** bcForcing, Sarray & a_x, Sarray & a_y, Sarray & a_z,
				Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J)
{
// 4D macros swap the last and first indices to compensate for different conventions between how 
// the Sarrays were allocated and how this routine was originally written
#define u(i,j,k,c) u_[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)]
#define q(i,j,k,c) q_[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)]
#define r(i,j,k,c) r_[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)]
#define s(i,j,k,c) s_[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)]
// 3D array macros are special cases of the 4D macros with c=1 and nc=1
#define x(i,j,k) x_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
#define y(i,j,k) y_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
#define z(i,j,k) z_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
#define J(i,j,k) J_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
#define mu(i,j,k) mu_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
#define lam(i,j,k) lam_[(i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb)]
// not necessary to store lambda + 2*mu in separate array
#define lam2mu(i,j,k) (lam(i,j,k) + 2.*mu(i,j,k))
  
// extract pointers for the macros
// 4D arrays
  double * u_=a_u.c_ptr();
  double * q_=a_q.c_ptr();
  double * r_=a_r.c_ptr();
  double * s_=a_s.c_ptr();
// 3D arrays
  double * x_=a_x.c_ptr();
  double * y_=a_y.c_ptr();
  double * z_=a_z.c_ptr();
  double * mu_=a_mu.c_ptr();
  double * lam_=a_lam.c_ptr();
  double * J_=a_J.c_ptr();
  
// all 3D/4D Sarrays must have the same number of grid points and the same starting/ending indices
  int m_nc = a_q.m_nc;
  int m_ni = a_q.m_ni;
  int m_nj = a_q.m_nj;
  int m_nk = a_q.m_nk;
// to mimic the original coding:
// setting starting indices to one 
// setting ending indices to equal the number of points in each dimension
  int m_ib = 1;
  int m_jb = 1;
  int m_kb = 1;
  int Nx = a_q.m_ni;
  int Ny = a_q.m_nj;
  int Nz = a_q.m_nk;

  int i, j, k, q, side, ind;
   
// only implemented for the low-k boundary
  side=4;
  
  {
    
    ind = 0;
#define E1(i,j,t1,t2,i2,t3,i3,t4) (t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2))
#define E12(i,j,t1,t2,i2,t3,i3,t4) (0.5*(t1(i,j,1)*t2(i,j,1,i2)*t3(i,j,1,i3)*t4(i,j,1)\
				 +t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2)))
#define E32(i,j,t1,t2,i2,t3,i3,t4) (0.5*(t1(i,j,3)*t2(i,j,3,i2)*t3(i,j,3,i3)*t4(i,j,3)\
				 +t1(i,j,2)*t2(i,j,2,i2)*t3(i,j,2,i3)*t4(i,j,2)))
    double Dr0[3],Dq0[3],Dsp[3],ugp[3],b[3],A[9],x[3],Ax[3];
    int info,N1,N2,N3 ,ipv[3]; 
    N1=3;N2=1;N3=3;
    k=2; ///NOTE!!! 
    for(int j=2; j<=Ny-1; j++)
      for(int i=2; i<=Nx-1; i++)
      {
	for(int c=1; c<=3;c++){
	  Dr0[c-1]=0.5*(u(i,j+1,2,c)-u(i,j-1,2,c));
	  Dq0[c-1]=0.5*(u(i+1,j,2,c)-u(i-1,j,2,c));
	  Dsp[c-1]=u(i,j,3,c)-u(i,j,2,c);
	}
	b[0]=-(E1(i,j,J,s,1,q,1,lam2mu)*Dq0[0]
	       +E1(i,j,J,s,1,r,1,lam2mu)*Dr0[0]	
	       +E1(i,j,J,s,1,q,2,lam)*Dq0[1]
	       +E1(i,j,J,s,1,r,2,lam)*Dr0[1]	 
	       +E1(i,j,J,s,1,q,3,lam)*Dq0[2]
	       +E1(i,j,J,s,1,r,3,lam)*Dr0[2]	 
	       +E1(i,j,J,s,2,q,1,mu)*Dq0[1]
	       +E1(i,j,J,s,2,r,1,mu)*Dr0[1]	 
	       +E1(i,j,J,s,2,q,2,mu)*Dq0[0]
	       +E1(i,j,J,s,2,r,2,mu)*Dr0[0]	 
	       +E1(i,j,J,s,3,q,1,mu)*Dq0[2]
	       +E1(i,j,J,s,3,r,1,mu)*Dr0[2]	 
	       +E1(i,j,J,s,3,q,3,mu)*Dq0[0]
	       +E1(i,j,J,s,3,r,3,mu)*Dr0[0]	 
	       +0.5*(
		 E32(i,j,J,s,1,s,1,lam2mu)*Dsp[0]
		 +E32(i,j,J,s,1,s,2,lam)*Dsp[1]
		 +E32(i,j,J,s,1,s,3,lam)*Dsp[2]
		 +E32(i,j,J,s,2,s,1,mu)*Dsp[1]
		 +E32(i,j,J,s,2,s,2,mu)*Dsp[0]
		 +E32(i,j,J,s,3,s,1,mu)*Dsp[2]
		 +E32(i,j,J,s,3,s,3,mu)*Dsp[0]));
	b[1]=-(E1(i,j,J,s,1,q,1,mu)*Dq0[1]
	       +E1(i,j,J,s,1,r,1,mu)*Dr0[1]	
	       +E1(i,j,J,s,1,q,2,mu)*Dq0[0]
	       +E1(i,j,J,s,1,r,2,mu)*Dr0[0]	 
	       +E1(i,j,J,s,2,q,2,lam2mu)*Dq0[1]
	       +E1(i,j,J,s,2,r,2,lam2mu)*Dr0[1]
	       +E1(i,j,J,s,2,r,1,lam)*Dr0[0]	 
	       +E1(i,j,J,s,2,q,1,lam)*Dq0[0]
	       +E1(i,j,J,s,2,r,3,lam)*Dr0[2]	 
	       +E1(i,j,J,s,2,q,3,lam)*Dq0[2]
	       +E1(i,j,J,s,3,q,2,mu)*Dq0[2]
	       +E1(i,j,J,s,3,r,2,mu)*Dr0[2]	 
	       +E1(i,j,J,s,3,q,3,mu)*Dq0[1]
	       +E1(i,j,J,s,3,r,3,mu)*Dr0[1]	 
	       +0.5*(E32(i,j,J,s,2,s,2,lam2mu)*Dsp[1]
		     +E32(i,j,J,s,2,s,1,lam)*Dsp[0]
		     +E32(i,j,J,s,2,s,3,lam)*Dsp[2]
		     +E32(i,j,J,s,1,s,1,mu)*Dsp[1]
		     +E32(i,j,J,s,1,s,2,mu)*Dsp[0]
		     +E32(i,j,J,s,3,s,2,mu)*Dsp[2]
		     +E32(i,j,J,s,3,s,3,mu)*Dsp[1]));
	    
	b[2]=-(E1(i,j,J,s,1,q,1,mu)*Dq0[2]
	       +E1(i,j,J,s,1,r,1,mu)*Dr0[2]	
	       +E1(i,j,J,s,1,q,3,mu)*Dq0[0]
	       +E1(i,j,J,s,1,r,3,mu)*Dr0[0]	 
	       +E1(i,j,J,s,3,q,3,lam2mu)*Dq0[2]
	       +E1(i,j,J,s,3,r,3,lam2mu)*Dr0[2]
	       +E1(i,j,J,s,3,r,1,lam)*Dr0[0]	 
	       +E1(i,j,J,s,3,q,1,lam)*Dq0[0]
	       +E1(i,j,J,s,3,r,2,lam)*Dr0[1]	 
	       +E1(i,j,J,s,3,q,2,lam)*Dq0[1]
	       +E1(i,j,J,s,2,q,2,mu)*Dq0[2]
	       +E1(i,j,J,s,2,r,2,mu)*Dr0[2]	 
	       +E1(i,j,J,s,2,q,3,mu)*Dq0[1]
	       +E1(i,j,J,s,2,r,3,mu)*Dr0[1]	 
	       +0.5*(E32(i,j,J,s,3,s,3,lam2mu)*Dsp[2]
		     +E32(i,j,J,s,3,s,1,lam)*Dsp[0]
		     +E32(i,j,J,s,3,s,2,lam)*Dsp[1]
		     +E32(i,j,J,s,1,s,1,mu)*Dsp[2]
		     +E32(i,j,J,s,1,s,3,mu)*Dsp[0]
		     +E32(i,j,J,s,2,s,2,mu)*Dsp[2]
		     +E32(i,j,J,s,2,s,3,mu)*Dsp[1]));


	A[0]=(
	  E12(i,j,J,s,1,s,1,lam2mu)
	  +E12(i,j,J,s,2,s,2,mu)
	  +E12(i,j,J,s,3,s,3,mu));
	A[3]=(E12(i,j,J,s,1,s,2,lam)
	      +E12(i,j,J,s,2,s,1,mu));
	A[6]=(E12(i,j,J,s,1,s,3,lam)
	      +E12(i,j,J,s,3,s,1,mu));
// u, v, w in v eq.
	A[1]=(E12(i,j,J,s,2,s,1,lam)
	      +E12(i,j,J,s,1,s,2,mu));
	A[4]=(E12(i,j,J,s,3,s,3,mu)
	      +E12(i,j,J,s,1,s,1,mu)
	      +E12(i,j,J,s,2,s,2,lam2mu));
	A[7]=(E12(i,j,J,s,2,s,3,lam)
	      +E12(i,j,J,s,3,s,2,mu));
	// u, v, w in w eq.
	A[2]=(E12(i,j,J,s,3,s,1,lam)
	      +E12(i,j,J,s,1,s,3,mu));
	A[5]=(E12(i,j,J,s,2,s,3,mu)
	      +E12(i,j,J,s,3,s,2,lam));
	A[8]=(E12(i,j,J,s,3,s,3,lam2mu)
	      +E12(i,j,J,s,1,s,1,mu)
	      +E12(i,j,J,s,2,s,2,mu));
	for(int c=0; c<9; c++)
	  A[c]*=0.5;

	for (int c=1; c<=3; c++)
	  x[c-1] = u(i,j,2,c) - u(i,j,1,c);
	    
	Ax[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
	Ax[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
	Ax[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];

	for(int c=1; c<=3;c++)
	{
// this routine is used to accumulate viscoelastic boundary stresses from each mechanism
	  bcForcing[side][ind+c-1] += Ax[c-1] - b[c-1];
	}
	ind += 3;
	    
      }
    
  }
  
  
#undef u
#undef x
#undef y
#undef z
#undef q
#undef r
#undef s
#undef J
#undef mu
#undef lam
#undef lam2mu
  
}


