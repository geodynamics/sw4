#include "EW.h"
#include <cstring>

// making directories
#include <errno.h>
//extern int errno;
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include <list>

#include "F77_FUNC.h"
extern "C" {
void F77_FUNC(exactmatfort,EXACTMATFORT)(int*, int*, int*, int*, int*, int*, double*, double*, double*, 
					 double*, double*, double*, double*, double*, double*, double*);
void F77_FUNC(wavepropbop_4, WAVEPROPBOP_4)(double *, double *, double *, double *, double *, double *, double *);
void F77_FUNC(varcoeffs4,VARCOEFFS4)(double *, double *);
void F77_FUNC(bopext4th,BOPEXT4TH)(double *, double *);

void F77_FUNC(dspev,DSPEV)(char & JOBZ, char & UPLO, int & N, double *AP, double *W, double *Z, int & LDZ, double *WORK, int & INFO);
void F77_FUNC(dgels,DGELS)(char & TRANS, int & M, int & N, int & NRHS, double *A, int & LDA, double *B, int & LDB, double *WORK, 
			   int & LWORK, int & INFO);
}

#define SQR(x) ((x)*(x))

//----------------------------------------------
void EW::setupRun( )
{
   if( mIsInitialized && proc_zero() )
      cout << " WARNING, calling setupRun twice " << endl;

  double time_start = MPI_Wtime();

// m_testing == true is one of the pointers to testing modes is assigned
// add other pointers to the list of testing modes as they get implemented
  m_testing = (m_twilight_forcing || m_point_source_test || m_lamb_test || m_rayleigh_wave_test ); 

// tmp
  if (mVerbose && proc_zero() )
    cout << " *** Testing = " << m_testing << endl;
    
  if( mVerbose && proc_zero() )
  {
    cout << "  Using Bjorn's (fast) parallel" << " IO library" << endl;
    if (m_pfs)
    {
      cout << "Assuming a PARALLEL file system" << endl;
      cout << "Writing images from (up to) " << m_nwriters << " procs" << endl;
    }
    else
    {
      cout << "Assuming a SERIAL file system." << endl;
    }
  }
  
// setup coefficients for SBP operators
  setupSBPCoeff();
  
  if( proc_zero() && mVerbose >=3)
  {
    int g;
    for( g=mNumberOfCartesianGrids-1; g>=0; g-- )
    {
      cout << "Info: Grid #" << g <<" min z-coordinate: " << m_zmin[g] << endl;
    }
    for( mNumberOfGrids-1; g>=0; g-- )
    {
      cout << "Info: mGridSize[" << g << "]=" << mGridSize[g] << endl;
    }
    cout << endl;
  }

// Set boundary conditions, if not done in input file.
  if( !mbcsSet )
    default_bcs( );

// check if there are any supergrid boundaries. If so, setup the supergrid tapering functions using
// the parameters m_sg_thickness, m_sg_transition, m_supergrid_damping_coefficient
  setup_supergrid( );

// assign m_bcType and m_onesided based on the global boundary conditions and parallel overlap boundaries. 
  assign_local_bcs(); 
  initializePaddingCells();

// for all sides with dirichlet, free surface, supergrid, or periodic boundary conditions,
// save the extent of the multi-D boundary window, and allocate arrays to hold the boundary forcing
  int wind[6], npts;
// with 2 ghost points, we need twice as many elements in the m_BCForcing arrays!!!

  for (int g=0; g<mNumberOfGrids; g++ )
  {
// tmp
//    printf("Allocating boundary forcing arrays for grid g=%i, # ghost points=%i\n", g, m_ghost_points);
    
    for(int side=0; side<6 ; side++ )
    {
      if (m_bcType[g][side] == bStressFree || m_bcType[g][side] == bDirichlet || 
	  m_bcType[g][side] == bSuperGrid  || m_bcType[g][side] == bPeriodic)
      {

// modify the window for stress free bc to only hold one plane
	if (m_bcType[g][side] == bStressFree)
	{
	  side_plane( g, side, wind, 1 );
// when calling side_plane with nGhost=1, you get the outermost grid plane
// for Free surface conditions, we apply the forcing on the boundary itself, i.e., just 
// inside the ghost points
// add/subtract the ghost point offset
	  if( side == 0 )
	  {
	    wind[0] += m_ghost_points;   wind[1] = wind[0];
	  }
	  else if( side == 1 )
	  {
	    wind[0] -= m_ghost_points;   wind[1] = wind[0];
	  }
	  else if( side == 2 )
	  {
	    wind[2] += m_ghost_points; wind[3] = wind[2];
	  }
	  else if( side == 3 )
	  {
	    wind[2]  -= m_ghost_points; wind[3] = wind[2];
	  }
	  else if( side == 4 )
	  {
	    wind[4] += m_ghost_points;
	    wind[5] = wind[4];
	  }
	  else
	  {
	    wind[4] -= m_ghost_points;
	    wind[5] = wind[4];
	  }
	}
	else // for Dirichlet, super grid, and periodic conditions, we
	     // apply the forcing directly on the ghost points
	{
	  side_plane( g, side, wind, m_ghost_points );
	}
	
	npts = (wind[5]-wind[4]+1)*
	  (wind[3]-wind[2]+1)*
	  (wind[1]-wind[0]+1);
// copy wind into m_BndryWindow
	for (int qq=0; qq<6; qq++)
	  m_BndryWindow[g][qq+side*6]=wind[qq];
	
// periodic conditions don't need any bc forcing array
	if (m_bcType[g][side] != bPeriodic)
	  m_NumberOfBCPoints[g][side] = npts;
	else
	  m_NumberOfBCPoints[g][side] = 0;
	
      } // end if

// tmp
      // printf("proc=%i side=%i, bc=%i, npts=%i, m_BndryWindow[g][side]=[%i %i %i %i %i %i]\n", m_myRank, side, 
      // 	     m_bcType[g][side], m_NumberOfBCPoints[g][side],
      // 	     m_BndryWindow[g][0 + side*6], m_BndryWindow[g][1 + side*6], m_BndryWindow[g][2 + side*6], 
      // 	     m_BndryWindow[g][3 + side*6], m_BndryWindow[g][4 + side*6], m_BndryWindow[g][5 + side*6]);
    } // end for side...
  } // end for grid...
  
// communicate the metric across processor boundaries (one-sided differences were used in setup_metric())
  if (topographyExists())
  {
    communicate_array( mJ, mNumberOfGrids-1 );
    communicate_array( mQ, mNumberOfGrids-1 );
    communicate_array( mR, mNumberOfGrids-1 );
    communicate_array( mS, mNumberOfGrids-1 );

// test interpolation of metric
//     double dist=0., X0, Y0, Z0, q0, r0, s0, qX[3], rX[3], sX[3];
//     int i, j, k, gTop = mNumberOfGrids-1;
//     for (k=m_kStart[gTop]; k<=m_kEnd[gTop]; k++)
//       for (j=m_jStart[gTop]; j<=m_jEnd[gTop]; j++)
// 	for (i=m_iStart[gTop]; i<=m_iEnd[gTop]; i++)
// 	{
// // evaluate the jacobian
// 	  find_curvilinear_derivatives_at_point( (double) i, (double) j, (double) k, qX, rX, sX );
// 	  dist += SQR(qX[0] - mQ(1,i,j,k)) + SQR(qX[1] - mQ(2,i,j,k)) + SQR(qX[2] - mQ(3,i,j,k))
// 	    + SQR(rX[0] - mR(1,i,j,k)) + SQR(rX[1] - mR(2,i,j,k)) + SQR(rX[2] - mR(3,i,j,k))
// 	    + SQR(sX[0] - mS(1,i,j,k)) + SQR(sX[1] - mS(2,i,j,k)) + SQR(sX[2] - mS(3,i,j,k));
// 	}
//     double totalDist;
//     MPI_Allreduce( &dist, &totalDist, 1, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );
//     if (m_mxyRank == 0)
//       printf("L2-error in metric of curvilinear mapping = %e\n", sqrt(totalDist));

  }

  MPI_Barrier(MPI_COMM_WORLD);

  string cachePath = mPath;
   
// tmp
  if ( mVerbose >= 3 )
   {
     int top = mNumberOfCartesianGrids-1;
     printf("=================Processor #%i index bounds====================\n"
	    "m_iStart=%i, m_iEnd=%i, m_global_nx=%i, m_jStart=%i, m_jEnd=%i, m_global_ny=%i\n",
	    m_myRank, 
	    m_iStart[top], m_iEnd[top], m_global_nx[top], 
	    m_jStart[top], m_jEnd[top], m_global_ny[top]);
     printf("=================Processor #%i interior index bounds====================\n"
	    "m_iStartInt=%i, m_iEndInt=%i, m_jStartInt=%i, m_jEndInt=%i\n",
	    m_myRank, 
	    m_iStartInt[top], m_iEndInt[top], 
	    m_jStartInt[top], m_jEndInt[top]);
     printf("=================Processor #%i Boundary Conditions in top grid====================\n"
	    "bc[0]=%i, bc[1]=%i, bc[2]=%i, bc[3]=%i, bc[4]=%i, bc[5]=%i\n",
	    m_myRank, 
	    m_bcType[top][0], m_bcType[top][1], m_bcType[top][2], m_bcType[top][3], m_bcType[top][4], m_bcType[top][5]);
  
   }

//tmp : make sure that all bc print statements are done. Then abort
//   MPI_Barrier(MPI_COMM_WORLD);
//   ASSERT(0)
   
//    if( m_output_load )
//       print_loadbalance_info();

  string saved_path = mPath;

  int beginCycle = 1;

// Initialize IO
  create_output_directory( );

// set material properties
  set_materials();

// evaluate resolution
  // double minvsoh;
  // compute_minvsoverh( minvsoh );
  // if (proc_zero())
  // {
  //   printf("\n***** PPW = minVs/h/maxFrequency ********\n");
  //   for (int g=0; g<mNumberOfCartesianGrids; g++)
  //   {
  //     printf("g=%i, h=%e, minVs/h=%g (Cartesian)\n", g, mGridSize[g], mMinVsOverH[g]);
  //   }
  //   if (topographyExists())
  //   {
  //     int g = mNumberOfGrids-1;
  //     printf("g=%i, h=%e, minVs/h=%g (curvilinear)\n", g, mGridSize[g], mMinVsOverH[g]);
  //   }
  //   printf("\n");
  // }
  
  if (usingSupergrid())
  {
// taper mu, lambda to near zero while making rho large at outflow (Dirichlet) boundaries
    supergrid_taper_material();
// 1-D damping coefficients in the (x,y,z) directions
    assign_supergrid_damping_arrays();
  }

// convert Qp and Qs to muVE, lambdaVE, and compute unrelaxed lambda, mu
  if (usingAttenuation())
  {
    if (m_myRank == 0) cout << "Viscoelastic not yet functional" << endl;
    return;
//     setup_viscoelastic( minvsoh );
  }
  
  if (proc_zero())
  {
    double lat[4],lon[4];
    computeGeographicCoord(0.0,           0.0,           lon[0], lat[0]);
    computeGeographicCoord(m_global_xmax, 0.0,           lon[1], lat[1]);
    computeGeographicCoord(m_global_xmax, m_global_ymax, lon[2], lat[2]);
    computeGeographicCoord(0.0,           m_global_ymax, lon[3], lat[3]);
    printf("Geographic coordinates of the corners of the computational grid:\n");
    for (int q=0; q<4; q++)
      printf("%i: Lon= %e, Lat=%e\n", q, lon[q], lat[q]);
    printf("\n");
  }
  
  if( mVerbose && proc_zero() )
    cout << "  Assigned material properties" << endl;

// compute time-step and number of time steps. 
// Note: SW4 always ends the simulation at mTmax, whether prefilter is enabled or not.
// This behavior is different from WPP
  computeDT( );

// should we initialize all images after the prefilter time offset stuff?

// Initialize image files: set time, tell images about grid hierarchy.
  initialize_image_files();
  if( mVerbose && proc_zero() )
    cout << "*** Initialized Images" << endl;

// // is the curvilinear grid ok?
//   if (topographyExists() && m_minJacobian <=0. && proc_zero()) // m_minJacobian is the global minimum and should be the same on all processes
//   {
//     printf("FATAL ERROR: min(Jacobian) = %e < 0 in the curvilinear grid.\n", m_minJacobian);
//     printf("You can try to increase zmax or decrease order in the topography command\n");
//     MPI_Abort(MPI_COMM_WORLD, 1);
//   }

// // Allocate work space for one sided operators in time stepping loop
//   int wksize=0;
//   for( int g=0; g < mNumberOfGrids ; g++ )
//   {
//      for( int side=0 ; side <= 1 ; side++ )
// 	if( m_onesided[g][side] == 1 )
// 	   wksize = wksize < mU[g].m_nj*mU[g].m_nk ? mU[g].m_nj*mU[g].m_nk : wksize; 
//      for( int side=2 ; side <= 3 ; side++ )
// 	if( m_onesided[g][side] == 1 )
// 	   wksize = wksize < mU[g].m_ni*mU[g].m_nk ? mU[g].m_ni*mU[g].m_nk : wksize;
//      for( int side=4 ; side <= 5 ; side++ )
// 	if( m_onesided[g][side] == 1 )
// 	   wksize = wksize < mU[g].m_nj*mU[g].m_ni ? mU[g].m_nj*mU[g].m_ni : wksize; 
//   }
//   double* wk = new double[4*wksize];

// if we got this far, the object should be ready for time-stepping
  mIsInitialized = true;

  double time_start_solve = MPI_Wtime();
  print_execution_time( time_start, time_start_solve, "start up phase" );
}

//-----------------------------------------------------------------------
void EW::preprocessSources( vector<Source*> & a_GlobalUniqueSources )
{
// make sure that the material model is in place
  if (!mIsInitialized)
  {
    if (proc_zero())
      printf("ERROR: Calling preprocessSources before the material model is initialized\n");
    return;
  }
  
// setup when there are point moment tensor sources of point forces,
// i.e., not twilight, energy conservation, or rayleigh surface wave test
  if( !m_twilight_forcing && !m_energy_test && !m_rayleigh_wave_test ) 
  {

// sanity checks for various testing modes
    if (m_testing)
    {
      bool sources_ok=true;
      
      if (mVerbose && proc_zero())
	cout << "Sanity testing of the source" << endl;

      if( (m_point_source_test || m_lamb_test) && a_GlobalUniqueSources.size() != 1 )
      {
	if (proc_zero())
	  cout << "Error: Point Source Test and Lamb Test must have one single source" << endl
	       << "       The number of defined sources is " << a_GlobalUniqueSources.size() << endl; 
	sources_ok=false;
      }

      if( m_point_source_test && !(a_GlobalUniqueSources[0]->getName() == "VerySmoothBump" ||
				   a_GlobalUniqueSources[0]->getName() == "SmoothWave" ||
				   a_GlobalUniqueSources[0]->getName() == "Gaussian") )
      {
	if (proc_zero())
	  cout << "Error: Point Source Test can only have source types" 
	       << " VerySmoothBump, SmoothWave, or Gaussian" << endl
	       << "  Input name is " << a_GlobalUniqueSources[0]->getName() << endl;
	sources_ok=false;
      }

      if( m_lamb_test && a_GlobalUniqueSources[0]->isMomentSource() )
      {
	if (proc_zero())
	  cout << "Error: Lamb's Test must have one point force" << endl
	       << "       The defined source is a moment tensor" << endl;
	sources_ok=false;
      }

      if( m_lamb_test )
      {
	double fx, fy, fz, z0, freq;
	a_GlobalUniqueSources[0]->getForces( fx, fy, fz );
	z0 = a_GlobalUniqueSources[0]->getZ0();
	freq = a_GlobalUniqueSources[0]->getFrequency();
	
	if ( !(a_GlobalUniqueSources[0]->getName() == "VerySmoothBump" &&
	       freq == 1.0 && z0 == 0.0 && fx == 0.0 && fy == 0.0) )
	{
	  if (proc_zero())
	    cout << "Error: Lamb Test assumes a 'VerySmoothBump' time function with freq=1" 
		 << " on z=0 with fx=0 and fy=0" << endl
		 << " The specified source has a '" << a_GlobalUniqueSources[0]->getName() 
		 << "' time function with freq=" << freq 
		 << " on z=" << z0 << " with fx=" << fx << " and fy=" << fy << endl;
	  sources_ok=false;
	}
      } // end if m_lamb_test

      if (!sources_ok) // can't use this setup
      {
	return; // not setting mSourcesOK to true should stop the calling program
      }
      
    } // end if m_testing
    else
    {
// Set up 'normal' sources for point_source_test, lamb_test, or standard seismic case.
// Correct source location for discrepancy between raw and smoothed topography
      for( unsigned int i=0 ; i < a_GlobalUniqueSources.size() ; i++ )
	a_GlobalUniqueSources[i]->correct_Z_level(); // also sets the ignore flag for sources that are above the topography

// limit max freq parameter (right now the raw freq parameter in the time function) Either rad/s or Hz depending on the time fcn
      // if (m_limit_source_freq)
      // {
      // 	if (mVerbose && proc_zero() )
      // 	  printf(" Limiting the freq parameter in all source time functions to the max value %e\n", m_source_freq_max);

      // 	for( int s=0; s < a_GlobalUniqueSources.size(); s++ ) 
      // 	  a_GlobalUniqueSources[s]->setMaxFrequency( m_source_freq_max );
       
      // } // end limit_source_freq
     
// check how deep the sources go
      double zSource, zMax=m_global_zmin, zMaxGlobal, zMin=m_global_zmax, zMinGlobal;
      for( int s=0; s < a_GlobalUniqueSources.size(); s++ ) 
      {
	zSource = a_GlobalUniqueSources[s]->getZ0( );
	if (zSource > zMax)
	  zMax = zSource;
	if (zSource < zMin)
	  zMin = zSource;
      }
// compute global max over all processors
      MPI_Allreduce( &zMax, &zMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator);
      MPI_Allreduce( &zMin, &zMinGlobal, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
      if (mVerbose && proc_zero() )
	printf(" Min source z-level: %e, max source z-level: %e\n", zMinGlobal, zMaxGlobal);

// Need to set the frequency to 1/dt for Dirac source
      for( int s=0 ; s  < a_GlobalUniqueSources.size(); s++ )
	 if( a_GlobalUniqueSources[s]->getTfunc() == iDirac )
	    a_GlobalUniqueSources[s]->setFrequency( 1.0/mDt );

// Modify the time functions if prefiltering is enabled
      if (m_prefilter_sources)
      {
// tell the filter about the time step and compute the second order sections
	m_filter_ptr->computeSOS( mDt );

// output details about the filter
	if (mVerbose>=3 && proc_zero() )
	  cout << *m_filter_ptr;
	
// 1. Make sure the smallest time offset is at least t0_min + (timeFcn dependent offset for centered fcn's)
	double dt0 = 0;
	double dt0loc, dt0max, t0_min;
	t0_min = m_filter_ptr->estimatePrecursor();
// tmp
	if ( proc_zero() )
	  printf("Filter precursor = %e\n", t0_min);
	
// old estimate for 2-pole low-pass Butterworth
//	t0_min = 4./m_filter_ptr->get_corner_freq2();
	
	for( int s=0; s < a_GlobalUniqueSources.size(); s++ ) 
	{
	  dt0loc = a_GlobalUniqueSources[s]->compute_t0_increase( t0_min );
	  
	  if( dt0loc > dt0 )
	    dt0 = dt0loc;
	}
// compute global max over all processors
	MPI_Allreduce( &dt0, &dt0max, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator);

// dt0max is the maxima over all dt0loc in all processors. 
// If it is positive, the t0 field in all source commands should be incremented 
// by at least this amount. Otherwise, there might be significant artifacts from 
// a sudden start of some source.
	if (dt0max > 0.)
	{
// Don't mess with t0.
// Instead, warn the user of potential transients due to unsmooth start
//	  for( int s=0; s < a_GlobalUniqueSources.size(); s++ ) 
//	    a_GlobalUniqueSources[s]->adjust_t0( dt0max );
	  if ( proc_zero() )
	    printf("\n*** WARNING: the 2 pass prefilter has an estimated precursor of length %e s\n"
		   "*** To avoid artifacts due to sudden startup, increase t0 in all source commands by at least %e\n\n",
		   t0_min, dt0max);
	}
// need to remember the time shift so we can compensate for it when writing sac and image files
//	m_t0Shift = dt0max;

// Do the filtering
	for( int s=0; s < a_GlobalUniqueSources.size(); s++ ) 
           a_GlobalUniqueSources[s]->filter_timefunc( m_filter_ptr, mTstart, mDt, mNumberOfTimeSteps );
      }

// TODO: check that t0 is large enough even when prefilter is NOT used

      if (proc_zero())
	getGMTOutput( a_GlobalUniqueSources );

// // Precompute grid point index for point sources, for better efficiency
//   int* ind=0;
//   if( m_forcing->use_input_sources() && m_point_sources.size() > 0 )
//   {
//     ind = new int[m_point_sources.size()];
//     for( int s= 0 ; s < m_point_sources.size() ; s++ ) 
//     {
// // tmp
// //      cout << *m_point_sources[s] << endl;
      
//       int g = m_point_sources[s]->m_grid;
//       ind[s] = mF[g].index(m_point_sources[s]->m_i0,m_point_sources[s]->m_j0,m_point_sources[s]->m_k0);
//     }
//   }

    } // end normal seismic setup

  } // end if ( !m_twilight_forcing && !m_energy_test && !m_rayleigh_wave_test ) 
   
  mSourcesOK = true;

} // end preprocessSources


//-----------------------------------------------------------------------
void EW::setupSBPCoeff()
{
  double gh2; // this coefficient is also stored in m_ghcof[0]
  if (mVerbose >=1 && m_myRank == 0)
    cout << "Setting up SBP boundary stencils" << endl;
  
// get coefficients for difference approximation of 2nd derivative with variable coefficients
//      call VARCOEFFS4( acof, ghcof )
  F77_FUNC(varcoeffs4,VARCOEFFS4)(m_acof, m_ghcof);
// get coefficients for difference approximation of 1st derivative
//      call WAVEPROPBOP_4( iop, iop2, bop, bop2, gh2, hnorm, sbop )
  F77_FUNC(wavepropbop_4,WAVEPROPBOP_4)(m_iop, m_iop2, m_bop, m_bop2, &gh2, m_hnorm, m_sbop);
// extend the definition of the 1st derivative tothe first 6 points
//      call BOPEXT4TH( bop, bope )
  F77_FUNC(bopext4th,BOPEXT4TH)(m_bop, m_bope);

}


//-----------------------------------------------------------------------
void EW::set_materials()
// Fill in material properties on all grids
// The material objects set velocities and stored Vs in mMu, Vp in mLambda. 
// Converted to (mu,lambda) by the call to convert_material_to_mulambda().
//
// After materials are set, call check_materials to make sure (mu,lambda,rho) values make sense
{  
  int g;
  
  if( !m_testing )
  {
// make sure all Nmat=m_number_material_surfaces material properties are defined
    
    if (m_number_material_surfaces > 0) // ifile stuff
    {
      if (m_materials.size() < m_number_material_surfaces)
      {
	cerr << "set_materials: ERROR: There are "<< m_number_material_surfaces << " material surfaces but only "
	     << m_materials.size() << " defined materials" << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
      }
// sort the m_materials vector in increasing id order...
      if (m_materials.size() > 1) // make sure there is more than 1 element
      {
	MaterialProperty *mp;
	for (int iStart=0; iStart < m_materials.size()-1; iStart++)
	{
	  int idmin=m_materials[iStart]->m_materialID;
	  for (int q=iStart+1; q<m_materials.size(); q++)
	  {
	    if (m_materials[q]->m_materialID < idmin)
	    {
// record new min value and swap elements q and iStart i m_materials vector
	      idmin = m_materials[q]->m_materialID;
	      mp = m_materials[q];
	      m_materials[q] = m_materials[iStart];
	      m_materials[iStart] = mp;
	    } // end if
	  } // end for q
	} // end for iStart
	       
      } // end if size > 1

// output a list of material id's
      if (proc_zero() && mVerbose>=3)
      {
	cout << "**** Material ID's: ********" << endl;
	for (int q=0; q<m_materials.size(); q++)
	  cout << "Material[" << q << "]->ID=" << m_materials[q]->m_materialID << endl;
      }
      
    }
    

// figure out if we may ignore some blocks (i.e. an efile command which is followed 
// by a block command covering all grid points). 
    int lastAllCoveringBlock=0;
    for( unsigned int b = 0 ; b < m_mtrlblocks.size() ; b++ )
    {
      if (m_mtrlblocks[b]->coversAllPoints())
        lastAllCoveringBlock=b;
    }
// tmp
    if (proc_zero())
    {
      if (lastAllCoveringBlock == 0)
	cout << "Considering all material blocks" << endl;
      else
	cout << "Only considering material blocks with index >= " << lastAllCoveringBlock << endl;
    } // end if proc_zero()
    
    for( unsigned int b = lastAllCoveringBlock ; b < m_mtrlblocks.size() ; b++ )
    {
      m_mtrlblocks[b]->set_material_properties(mRho, mMu, mLambda, mQs, mQp); // this is where all the assignments are done

    } // end for
   
// extrapolate to define material properties above the free surface (topography)
   g = mNumberOfGrids-1;
   
   bool linearExtrapolation=false;
// note that material thresholding for vs and vp happens further down in this procedure
   extrapolateInZ(mRho[g], false, 0., linearExtrapolation); 
   extrapolateInZ(mLambda[g], false, 0., linearExtrapolation); 
   extrapolateInZ(mMu[g], false, 0., linearExtrapolation);

   if( m_use_attenuation )
   {
     extrapolateInZ(mQs[g], false, 0., linearExtrapolation); 
     extrapolateInZ(mQp[g], false, 0., linearExtrapolation); 
   }

// extrapolate material properties to mesh refinement boundaries (e.g. for doing the LOH cases more accurately)
    if (proc_zero() && mVerbose>=3)
    {
      printf("setMaterials> mMaterialExtrapolate = %i, mNumberOfCartesianGrids=%i\n", mMaterialExtrapolate, mNumberOfCartesianGrids);
    }
    
    if (mMaterialExtrapolate > 0 && mNumberOfCartesianGrids > 1)
    {
      int kFrom;
      for (g=0; g<mNumberOfCartesianGrids; g++)
      {
	if (g < mNumberOfCartesianGrids-1) // extrapolate to top
	{
	  kFrom = m_kStart[g]+mMaterialExtrapolate;

	  if (proc_zero() && mVerbose>=3)
	    printf("setMaterials> top extrapol, g=%i, kFrom=%i, kStart=%i\n", g, kFrom, m_kStart[g]);

	  for (int k = m_kStart[g]; k < kFrom; ++k)
	    for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
	      for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
	      {
		mRho[g](i,j,k) = mRho[g](i,j,kFrom);
		mMu[g](i,j,k)  = mMu[g](i,j,kFrom);
		mLambda[g](i,j,k) = mLambda[g](i,j,kFrom);
	      }

	  if( m_use_attenuation )
	  {
	    for (int k = m_kStart[g]; k < kFrom; ++k)
	      for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
		for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
		{
		  mQs[g](i,j,k) = mQs[g](i,j,kFrom);
		  mQp[g](i,j,k) = mQp[g](i,j,kFrom);
		}
	  }
	  
	} // end extrapolat to top

	if (g > 0) // extrapolate to bottom
	{
	  kFrom = m_kEnd[g]-mMaterialExtrapolate;

	  if (proc_zero() && mVerbose>=3)
	    printf("setMaterials> bottom extrapol, g=%i, kFrom=%i, kEnd=%i\n", g, kFrom, m_kEnd[g]);

	  for (int k = kFrom+1; k <= m_kEnd[g]; ++k)
	    for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
	      for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
	      {
		mRho[g](i,j,k) = mRho[g](i,j,kFrom);
		mMu[g](i,j,k)  = mMu[g](i,j,kFrom);
		mLambda[g](i,j,k) = mLambda[g](i,j,kFrom);
	      }

	  if( m_use_attenuation )
	  {
	    for (int k = kFrom+1; k <= m_kEnd[g]; ++k)
	      for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
		for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
		{
		  mQs[g](i,j,k) = mQs[g](i,j,kFrom);
		  mQp[g](i,j,k) = mQp[g](i,j,kFrom);
		}
	  }
	  
	} // end extrapolate to bottom
	
      } // end for g      
    } // end if mMaterialExtrapolate > 0 ...
    
// tmp
//    printf("\n useVelocityThresholds=%i vpMin=%e vsMin=%e\n\n", m_useVelocityThresholds, m_vpMin, m_vsMin);
// threshold material velocities
    if (m_useVelocityThresholds)
    {
      for (g=0; g<mNumberOfGrids; g++)
	for (int k = m_kStart[g]; k <= m_kEnd[g]; k++)
	    for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
	      for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
	      {
		if (mMu[g](i,j,k) < m_vsMin) mMu[g](i,j,k) = m_vsMin;
		if (mLambda[g](i,j,k) < m_vpMin) mLambda[g](i,j,k) = m_vpMin;
	      }
    }
    
    convert_material_to_mulambda( );

// do the viscoelastic materials later (after estimating the resolution)
  } // end if !m_testing, i.e., not Twilight, point source or Lamb's test
  else if (m_twilight_forcing) 
  {
// tmp
    if (proc_zero())
      cout << "******************************" << endl
	   << " ASSIGNING TWILIGHT MATERIALS " << endl
	   << "******************************" << endl;

// For some forcings (such as twilight forcing) the material is set here.
      double xP, yP, zP;
      
      int ifirst, ilast, jfirst, jlast, kfirst, klast;
      double *rho_ptr, *mu_ptr, *la_ptr, h, zmin, omm, phm, amprho, ampmu, ampla;
	
      for (g=0; g<mNumberOfCartesianGrids; g++)
      {
	rho_ptr = mRho[g].c_ptr();
	mu_ptr  = mMu[g].c_ptr();
	la_ptr  = mLambda[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	h = mGridSize[g];
	zmin = m_zmin[g];
	omm = m_twilight_forcing->m_momega;
	phm = m_twilight_forcing->m_mphase;
	amprho = m_twilight_forcing->m_amprho;
	ampmu = m_twilight_forcing->m_ampmu;
	ampla = m_twilight_forcing->m_amplambda;
	
     //  subroutine exactmatfort( ifirst, ilast, jfirst, jlast, kfirst, 
     // +     klast, rho, mu, la, omm, phm, amprho, ampmu, amplambda, h, 
     // +     zmin )
	F77_FUNC(exactmatfort,EXACTMATFORT)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					    &klast, rho_ptr, mu_ptr, la_ptr, &omm, &phm, 
					    &amprho, &ampmu, &ampla, &h, &zmin );

// Need to communicate across material boundaries
	  communicate_array( mRho[g], g );
	  communicate_array( mMu[g], g );
	  communicate_array( mLambda[g], g );
      }
//       if (topographyExists())
//       {
// 	g = mNumberOfGrids-1;
// 	for (int k=m_kStart[g]; k<=m_kEnd[g]; k++)
// 	  for (int j=m_jStart[g]; j<=m_jEnd[g]; j++)
// 	    for (int i=m_iStart[g]; i<=m_iEnd[g]; i++)
// 	    {
// 	      xP = mX(i,j,k);
// 	      yP = mY(i,j,k);
// 	      zP = mZ(i,j,k);
// // what about Qp, Qs???
// 	      m_twilight_forcing->get_mtrl( xP, yP, zP, mRho[g](i,j,k), mMu[g](i,j,k), mLambda[g](i,j,k) );
// 	    }
//   // Need this for Energy testing, random material will not agree on processor boundaries.
// 	  communicate_array( mRho[g], g );
// 	  communicate_array( mMu[g], g );
// 	  communicate_array( mLambda[g], g );
//       }
      
  } // end material set by forcing mode (for testing)
  else if ( m_point_source_test )
  {
     for (g=0; g<mNumberOfCartesianGrids; g++)
     {
	mRho[g].set_value( m_point_source_test->m_rho );
	mMu[g].set_value( m_point_source_test->m_mu );
	mLambda[g].set_value( m_point_source_test->m_lambda );
     }
  }
  else if ( m_lamb_test )
  {
     for (g=0; g<mNumberOfCartesianGrids; g++)
     {
	mRho[g].set_value( m_lamb_test->m_rho );
	mMu[g].set_value( m_lamb_test->m_mu );
	mLambda[g].set_value( m_lamb_test->m_lambda );
     }
  }
  else if ( m_rayleigh_wave_test )
  {
     for (g=0; g<mNumberOfCartesianGrids; g++)
     {
	mRho[g].set_value( m_rayleigh_wave_test->m_rho );
	mMu[g].set_value( m_rayleigh_wave_test->m_mu );
	mLambda[g].set_value( m_rayleigh_wave_test->m_lambda );
     }
  }
  else if ( m_energy_test )
  {
     double cpocs = m_energy_test->m_cpcsratio;
     for (g=0; g<mNumberOfCartesianGrids; g++)
     {
	double* rho_ptr    = mRho[g].c_ptr();
	double* mu_ptr     = mMu[g].c_ptr();
	double* lambda_ptr = mLambda[g].c_ptr();
	for( int i=0 ; i < (m_iEnd[g]-m_iStart[g]+1)*(m_jEnd[g]-m_jStart[g]+1)*(m_kEnd[g]-m_kStart[g]+1); i++ )
	{
	   rho_ptr[i]    = drand48()+2;
	   mu_ptr[i]     = drand48()+2;
           lambda_ptr[i] = mu_ptr[i]*(cpocs*cpocs-2)+drand48();
	}
	communicate_array( mRho[g], g );
	communicate_array( mMu[g], g );
	communicate_array( mLambda[g], g );
     }
  }
  check_materials( );

}

//-----------------------------------------------------------------------
void EW::create_output_directory( )
{
   if (proc_zero()) 
   {

     cout << "----------------------------------------------------" << endl
	  << " Making Output Directory: " << mPath << endl
	  << "\t\t" << endl;

      // Create directory where all these files will be written.
      int err = mkdirs(mPath);

      if (err == 0)
	cout << "... Done!" << endl
	     << "----------------------------------------------------" << endl;
      else
      {
// fatal error
	cerr << endl << "******** Failed to create the output directory *******" << endl << endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

// check that we have write permission on the directory
      if (access(mPath.c_str(),W_OK)!=0)
      {
// fatal error
	cerr << endl << "Error: No write permission on output directory: " << mPath << endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
   }
  // Let processor 0 finish first!
   cout.flush();  cerr.flush();
   MPI_Barrier(MPI_COMM_WORLD);

// Check that the mPath directory exists from all processes
   struct stat statBuf;
   int statErr = stat(mPath.c_str(), &statBuf);
   CHECK_INPUT(statErr == 0 && S_ISDIR(statBuf.st_mode), "Error: " << mPath << " is not a directory" << endl);
   
// check that all processes have write permission on the directory
   CHECK_INPUT(access(mPath.c_str(),W_OK)==0,
	   "Error: No write permission on output directory: " << mPath << endl);
}

//-----------------------------------------------------------------------
void EW::computeDT()
{
  if (mVerbose >= 1 && proc_zero())
  {
    printf("*** computing the time step ***\n");
  }
  
   double dtloc=1.e10, dtGP;
   mDt = dtloc;

   double factor = 4;
   double loceig;
   
// temporarily changing the CFL number to suit 2nd order time stepping (and agree with test3 fortran code)
//   mCFL = 0.9;

   int g;
   for (int g=0; g<mNumberOfCartesianGrids; g++)
   {
     for (int k=m_kStart[g]; k<=m_kEnd[g]; k++)
       for (int j=m_jStart[g]; j<=m_jEnd[g]; j++)
	 for (int i=m_iStart[g]; i<=m_iEnd[g]; i++)
	 {
	    loceig = factor*mMu[g](i,j,k) + mLambda[g](i,j,k);
	    for( int a=0 ; a < m_number_mechanisms ; a++ )
	       loceig += factor*mMuVE[g][a](i,j,k) + mLambdaVE[g][a](i,j,k);
            loceig /= mRho[g](i,j,k);
	    //	   loceig = (factor*mMu[g](i,j,k) + mLambda[g](i,j,k) )/mRho[g](i,j,k);
	   dtGP = mCFL*mGridSize[g]/sqrt( loceig );
	   //	   dtloc = min(dtloc,dtGP);
	   dtloc = dtloc < dtGP ? dtloc : dtGP;
	 }
   }

// look at the curvilinear grid
   double dtCurv=1.e18, dtCurvGP;
   if (topographyExists())
   {
// tmp
//     FILE *fp=fopen("eigen.ext", "w");
     
     g = mNumberOfGrids-1;
     double Amat[6], la, mu, la2mu;
     int N=3, LDZ=1, INFO;
     char JOBZ='N', UPLO='L';
     double W[3], Z[1], WORK[9];
// do consider ghost points (especially the ghost line above the topography might be important)
     for (int k=m_kStart[g]; k<=m_kEnd[g]; k++)
       for (int j=m_jStart[g]; j<=m_jEnd[g]; j++)
	 for (int i=m_iStart[g]; i<=m_iEnd[g]; i++)
	 {
	   la = mLambda[g](i,j,k);
	   mu = mMu[g](i,j,k);
	   for( int a = 0 ; a < m_number_mechanisms ; a++ )
	   {
	      la += mLambdaVE[g][a](i,j,k);
              mu += mMuVE[g][a](i,j,k);
	   }
	   la2mu = la + 2.*mu;
	   
// A11
	   Amat[0] = -4.*(SQR(mQ(1,i,j,k))*la2mu + SQR(mQ(2,i,j,k))*mu + SQR(mQ(3,i,j,k))*mu 
			  + SQR(mR(1,i,j,k))*la2mu + SQR(mR(2,i,j,k))*mu + SQR(mR(3,i,j,k))*mu
			  + SQR(mS(1,i,j,k))*la2mu + SQR(mS(2,i,j,k))*mu + SQR(mS(3,i,j,k))*mu);
// A21 = A12
	   Amat[1] = -4.*(mQ(1,i,j,k)*mQ(2,i,j,k) + mR(1,i,j,k)*mR(2,i,j,k) + mS(1,i,j,k)*mS(2,i,j,k))*(mu+la);
// A31 = A13	   
	   Amat[2] = -4.*(mQ(1,i,j,k)*mQ(3,i,j,k) + mR(1,i,j,k)*mR(3,i,j,k) + mS(1,i,j,k)*mS(3,i,j,k))*(mu+la);
// A22	   
	   Amat[3] = -4.*(SQR(mQ(1,i,j,k))*mu + SQR(mQ(2,i,j,k))*la2mu + SQR(mQ(3,i,j,k))*mu 
			  + SQR(mR(1,i,j,k))*mu + SQR(mR(2,i,j,k))*la2mu + SQR(mR(3,i,j,k))*mu
			  + SQR(mS(1,i,j,k))*mu + SQR(mS(2,i,j,k))*la2mu + SQR(mS(3,i,j,k))*mu);
// A32 = A23
	   Amat[4] = -4.*(mQ(2,i,j,k)*mQ(3,i,j,k) + mR(2,i,j,k)*mR(3,i,j,k) + mS(2,i,j,k)*mS(3,i,j,k))*(mu+la);
// A33
	   Amat[5] = -4.*(SQR(mQ(1,i,j,k))*mu + SQR(mQ(2,i,j,k))*mu + SQR(mQ(3,i,j,k))*la2mu 
			  + SQR(mR(1,i,j,k))*mu + SQR(mR(2,i,j,k))*mu + SQR(mR(3,i,j,k))*la2mu
			  + SQR(mS(1,i,j,k))*mu + SQR(mS(2,i,j,k))*mu + SQR(mS(3,i,j,k))*la2mu);
// calculate eigenvalues of symmetric matrix
	   F77_FUNC(dspev,DSPEV)(JOBZ, UPLO, N, Amat, W, Z, LDZ, WORK, INFO);
	   if (INFO != 0)
	   {
	     printf("ERROR: computeDT: dspev returned INFO = %i for grid point (%i, %i, %i)\n", INFO, i, j, k);
	     printf("lambda = %e, mu = %e\n", la, mu);
	     MPI_Abort(MPI_COMM_WORLD, 1);
	   }
	   
// tmp
//	   fprintf(fp,"%i %i %i %e %e %e\n", i, j, k, W[0], W[1], W[2]);

// eigenvalues in ascending order: W[0] < W[1] < W[2]
	   if (W[0] >= 0.)
	   {
	     printf("ERROR: computeDT: determining eigenvalue is non-negative; W[0] = %e at curvilinear grid point (%i, %i, %i)\n", W[0], i, j, k);
	     MPI_Abort(MPI_COMM_WORLD, 1);
	   }
// local time step
	   dtCurvGP = dtGP = mCFL*sqrt(4.*mRho[g](i,j,k)/(-W[0]));
	   dtloc = min(dtloc, dtGP);
	   dtCurv = min(dtCurv, dtCurvGP);
	 }
// tmp
//     fclose(fp);
   } // end if topographyExists()
   
   mDt = dtloc;

// compute the global minima
    MPI_Allreduce( &dtloc, &mDt, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);

// global minima for curvilinear grid
    if (topographyExists())
    {
      double dtCurvGlobal;
      MPI_Allreduce( &dtCurv, &dtCurvGlobal, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
      if (mVerbose >= 1 && proc_zero())
      {
	printf("INFO: Smallest stable time step for curvilinear grid only: %e\n", dtCurvGlobal);
      }
    }

    if (mVerbose >= 1 && proc_zero())
    {
      cout << "order of accuracy=" << mOrder << " CFL=" << mCFL << " prel. time step=" << mDt << endl;
    }
    
    if (mTimeIsSet)
    {
// constrain the dt based on the goal time
//      VERIFY2(mTmax > mTstart,"*** ERROR: Tstart is greater than Tmax! ***");  
      mNumberOfTimeSteps = static_cast<int> ((mTmax - mTstart) / mDt + 0.5); 
      mNumberOfTimeSteps = (mNumberOfTimeSteps==0)? 1: mNumberOfTimeSteps;
// the resulting mDt could be slightly too large, because the numberOfTimeSteps is rounded to the nearest int
      mDt = (mTmax - mTstart) / mNumberOfTimeSteps;
    }
}

//-----------------------------------------------------------------------
int EW::mkdirs(const string& path)
{

   string pathTemp(path.begin(), path.end()); 
   //-----------------------------------------------------------------
   // Recursively call stat and then mkdir on each sub-directory in 'path'
   //-----------------------------------------------------------------
   string sep = "/";
   char* token = strtok(const_cast<char*>(pathTemp.c_str()), sep.c_str());

   stringstream pathsofar;

// for checking the status:
   struct stat statBuf;
   int statErr;
   
   // If there's a leading slash, put it back on...
   if (strncmp(pathTemp.c_str(), sep.c_str(), 1) == 0) pathsofar << sep;

   while (token != NULL)
   {
      pathsofar << token << sep;

// test: check the status of the path so far...
//      cout << "Calling stat() on path: " << pathsofar.str() << endl;
      statErr = stat(pathsofar.str().c_str(), &statBuf);
      if (statErr == 0)
      {
//	cout << "stat() returned successfully." << endl;
	if ( S_ISDIR(statBuf.st_mode) )
	{
//	  cout << "stat() says: '" << pathsofar.str() << "' is a directory." << endl;
// it already exists, this is okay, let's get the next directory in the string and skip to the while statement
	  token = strtok(NULL, sep.c_str());
	  continue;
	}
	else
	{
	  cerr << "stat() says: '" << pathsofar.str() << "' is not a directory." << endl;
	// real error, let's bail...
	  return -1;
	}
	
      }
      else
      {
//	cerr << "stat() returned an error code." << endl;
	if (errno == EACCES)
	{
	  cerr << "Error: **Search permission is denied for one of the directories in the path prefix of " << pathsofar.str() << endl;
	  return -1;
	}
	else if (errno == ENOTDIR)
	{
	  cerr << "Error: **A component of the path '" <<  pathsofar.str() << "' is not a directory. " << endl;
	  return -1;
	}
 	else if (errno == ENOENT)
 	{
// this means that we need to call mkdir to create the directory
	  if (mVerbose >=2) 
	    cout << "Info: **stat returned ENOENT (the path does not exist, or the path " << endl
		 << "      is an empty string) " << pathsofar.str() << endl;
 	}
	else
	{
	  if (mVerbose >=2) 
	    cout << "Info: **stat returned other error code for path: " << pathsofar.str() << endl;
	}
      }

// if we got this far, then 'pathsofar' does not exists

// tmp
      if (mVerbose >=2) cout << "Calling mkdir() on path: " << pathsofar.str() << endl;
// old code for recursively making the output directory
       if (mkdir(pathsofar.str().c_str(), 
                S_IWUSR | S_IXUSR | S_IRUSR | S_IRGRP | S_IXGRP ) // why do we need group permissions?
          == -1)
      {
	if (mVerbose >=2) cout << "mkdir() returned an error code." << endl;
         // check error conditions
	if (errno == EEXIST)
	{
// can this ever happen since we called stat(), which said that the directory did not exist ???
	  if (mVerbose >=2) cout << "Info: ** The directory already exists:" << pathsofar.str() << endl;
	  
	  // it already exists, this is okay!
	  token = strtok(NULL, sep.c_str());
	  continue;
	}
	else if (errno == EACCES)
	  cerr << "Error: **Write permission is denied for the parent directory in which the new directory is to be added." << pathsofar.str() << endl;
	else if (errno == EMLINK)
	  cerr << "Error: **The parent directory has too many links (entries)." << 
	    pathsofar.str() << endl;
	else if (errno == ENOSPC)
	  cerr << "Error: **The file system doesn't have enough room to create the new directory." <<
	    pathsofar.str() << endl;
	else if (errno == EROFS)
	  cerr << "Error: **  The parent directory of the directory being created is on a read-only file system and cannot be modified." << pathsofar.str() << endl;
	else if (errno == ENOSPC)
	  cerr << "Error: ** The new directory cannot be created because the user's disk quota is exhausted." << pathsofar.str() << endl;
	// real error, let's bail...
	return -1;
      }
      else
      {
	if (mVerbose >=2) cout << "mkdir() returned successfully." << endl;

// are there more directories to be made?
	token = strtok(NULL, sep.c_str());
      }
   }
   return 0;
}

//-----------------------------------------------------------------------
void EW::setup_supergrid( )
{
  if (mVerbose >= 3 && proc_zero())
    cout << "*** Inside setup_supergrid ***" << endl;
  
// check to see if there are any supergrid boundary conditions
  for (int b=0; b<6; b++)
  {
    if (mbcGlobalType[b] == bSuperGrid)
    {
      m_use_supergrid=true;
    }
  }

  if (mVerbose && proc_zero() && m_use_supergrid)
    cout << "Detected at least one boundary with supergrid conditions" << endl;
  
  int gTop = mNumberOfCartesianGrids-1;
  int gBot = 0;

  m_supergrid_taper_x.define_taper( (mbcGlobalType[0] == bSuperGrid), 0.0, (mbcGlobalType[1] == bSuperGrid), m_global_xmax, 
				    m_sg_gp_thickness*mGridSize[gTop], m_sg_gp_transition*mGridSize[gTop] );
  m_supergrid_taper_y.define_taper( (mbcGlobalType[2] == bSuperGrid), 0.0, (mbcGlobalType[3] == bSuperGrid), m_global_ymax, 
				    m_sg_gp_thickness*mGridSize[gTop], m_sg_gp_transition*mGridSize[gTop] );

// Note that we use the grid size on the bottom grid to define the transition width in the z-direction 
  if (topographyExists()) // the taper function in z is currently not defined for a non-planar top surface
    m_supergrid_taper_z.define_taper( false, 0.0, (mbcGlobalType[5] == bSuperGrid), m_global_zmax, 
				      m_sg_gp_thickness*mGridSize[gBot], m_sg_gp_transition*mGridSize[gBot] );
  else
    m_supergrid_taper_z.define_taper( (mbcGlobalType[4] == bSuperGrid), 0.0, (mbcGlobalType[5] == bSuperGrid), m_global_zmax, 
				      m_sg_gp_thickness*mGridSize[gBot], m_sg_gp_transition*mGridSize[gBot] );
// tmp
//   if (mVerbose >= 2 && proc_zero())
//   {
//     printf("********** Super-grid parameters (x, y, z)-directions:\n");
//     m_supergrid_taper_x.print_parameters();
//     m_supergrid_taper_y.print_parameters();
//     m_supergrid_taper_z.print_parameters();
//   }
  
}

//-----------------------------------------------------------------------
void EW::assign_supergrid_damping_arrays()
{
  int g, i, j, k, topCartesian;
  double x, y, z;
  
// resize the vectors for the pointers
  m_sg_dc_x.resize(mNumberOfGrids);
  m_sg_dc_y.resize(mNumberOfGrids);
  m_sg_dc_z.resize(mNumberOfGrids);
  
// allocate storage for 1-D damping coefficients on each grid
  for( g=0 ; g<mNumberOfGrids; g++) 
  {
    m_sg_dc_x[g] = new double[m_iEnd[g]-m_iStart[g]+1];
    m_sg_dc_y[g] = new double[m_jEnd[g]-m_jStart[g]+1];
    m_sg_dc_z[g] = new double[m_kEnd[g]-m_kStart[g]+1];
  }
#define dcx(i,g) (m_sg_dc_x[g])[i-m_iStart[g]]
#define dcy(j,g) (m_sg_dc_y[g])[j-m_jStart[g]]
#define dcz(k,g) (m_sg_dc_z[g])[k-m_kStart[g]]
  topCartesian = mNumberOfCartesianGrids-1;

// Note: compared to WPP2, we don't need to center the damping coefficients on the half-point anymore,
// because the damping term is now 4th order: D+D-( a(x) D+D- ut(x) )

  for( g=0 ; g<mNumberOfGrids; g++)  
  {
    for( i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
    {
      x = (i-1)*mGridSize[g];
      dcx(i,g) = m_supergrid_taper_x.dampingCoeff(x);
    }
    for( j = m_jStart[g] ; j <= m_jEnd[g] ; j++ )
    {
      y = (j-1)*mGridSize[g];
      dcy(j,g) = m_supergrid_taper_y.dampingCoeff(y);
    }
    if (g > topCartesian) // must be the curvilinear grid
    {
// supergrid damping in the vertical (k-) direction on a curvilinear grid is not defined
      for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
      {
	dcz(k,g) = 0.;
      }
    }
    else
    {
      for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
      {
	z = m_zmin[g] + (k-1)*mGridSize[g];
	dcz(k,g) = m_supergrid_taper_z.dampingCoeff(z);
      }
    }
    
  } // end for g...
  
// tmp: save the damping coefficients in matlab format
  // FILE *fp;
  // char fName[80];
 
  // g = 0;
  
  // sprintf(fName,"dcx-p%i.ext", m_myRank);
  // fp = fopen(fName,"w");
  // for( i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
  // {
  //   x = (i-1)*mGridSize[g];
  //   fprintf(fp,"%e %e\n", x, dcx(i,g));
  // }
  // fclose(fp);

  // sprintf(fName,"dcy-p%i.ext", m_myRank);
  // fp = fopen(fName,"w");
  // for( j = m_jStart[g] ; j <= m_jEnd[g] ; j++ )
  // {
  //   y = (j-1)*mGridSize[g];
  //   fprintf(fp,"%e %e\n", y, dcy(j,g));
  // }
  // fclose(fp);

  // sprintf(fName,"dcz-p%i.ext", m_myRank);
  // fp = fopen(fName,"w");
  // for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
  // {
  //   z = m_zmin[g] + (k-1)*mGridSize[g];
  //   fprintf(fp,"%e %e\n", z, dcz(k,g));
  // }
  // fclose(fp);

#undef dcx
#undef dcy
#undef dcz
}

//-----------------------------------------------------------------------
void EW::supergrid_taper_material( )
{
  double x, y, z;
  int g;
  if (mVerbose >= 3 && proc_zero())
    cout << "*** Inside supergrid_taper_material ***" << endl;
  
  for( g=0 ; g<mNumberOfCartesianGrids; g++) // do Cartesian grids first
  {
    for( int k = m_kStart[g] ; k <= m_kEnd[g]; k++ )
      for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	for( int i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	{
	  x = (i-1)*mGridSize[g];
	  y = (j-1)*mGridSize[g];
	  z = m_zmin[g] + (k-1)*mGridSize[g];
	  
	  mMu[g](i,j,k) *= m_supergrid_taper_x.velocityCoeff(x) * m_supergrid_taper_y.velocityCoeff(y) * m_supergrid_taper_z.velocityCoeff(z);
	  mLambda[g](i,j,k) *= m_supergrid_taper_x.velocityCoeff(x) * m_supergrid_taper_y.velocityCoeff(y) * m_supergrid_taper_z.velocityCoeff(z);
	  mRho[g](i,j,k) /=  m_supergrid_taper_x.velocityCoeff(x) * m_supergrid_taper_y.velocityCoeff(y) * m_supergrid_taper_z.velocityCoeff(z);
	} // end for i,j,k
  } // end for g

  if (topographyExists())
  {
    g = mNumberOfGrids - 1;
    for( int k = m_kStart[g] ; k <= m_kEnd[g]; k++ )
      for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	for( int i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	{
	  x = (i-1)*mGridSize[g];
	  y = (j-1)*mGridSize[g];
// NOTE: no tapering in z because there is (usually) a free surface boundary on top	  
	  mMu[g](i,j,k) *= m_supergrid_taper_x.velocityCoeff(x) * m_supergrid_taper_y.velocityCoeff(y);
	  mLambda[g](i,j,k) *= m_supergrid_taper_x.velocityCoeff(x) * m_supergrid_taper_y.velocityCoeff(y);
	  mRho[g](i,j,k) /=  m_supergrid_taper_x.velocityCoeff(x) * m_supergrid_taper_y.velocityCoeff(y);
	  
	} // end for i,j,k
  } // end if topographyExists()
  
}

