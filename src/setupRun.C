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
#include <cstring>

// making directories
#include <errno.h>
//extern int errno;
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include <list>
#include <unistd.h>

#include "F77_FUNC.h"
extern "C" {
   void tw_ani_stiff(int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, double h, double zmin,
                     double omm, double phm, double amprho, double *rho, double *phc, double *cm);

   void tw_ani_curvi_stiff(int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, double *xx, double *yy, double *zz,
                     double omm, double phm, double amprho, double *rho, double *phc, double *cm);
   
   void F77_FUNC(exactmatfort,EXACTMATFORT)(int*, int*, int*, int*, int*, int*, double*, double*, double*, 
					 double*, double*, double*, double*, double*, double*, double*);
   void F77_FUNC(exactmatfortc,EXACTMATFORTC)( int*, int*, int*, int*, int*, int*, double*, double*, double*,
                                            double*, double*, double*, double*, double*, double*, double*, double*); 
void F77_FUNC(wavepropbop_4, WAVEPROPBOP_4)(double *, double *, double *, double *, double *, double *, double *);
void F77_FUNC(varcoeffs4,VARCOEFFS4)(double *, double *);
void F77_FUNC(bopext4th,BOPEXT4TH)(double *, double *);

void F77_FUNC(dspev,DSPEV)(char & JOBZ, char & UPLO, int & N, double *AP, double *W, double *Z, int & LDZ, double *WORK, int & INFO);
void F77_FUNC(dgels,DGELS)(char & TRANS, int & M, int & N, int & NRHS, double *A, int & LDA, double *B, int & LDB, double *WORK, 
			   int & LWORK, int & INFO);
void F77_FUNC(randomfield3d,RANDOMFIELD3D)( int *, int *, int *, int *, int *, int *, int*, int*, int*, 
					    int*, double*, double*, double*, double*, double*, int*, double*, int*, int* );
void F77_FUNC(randomfield3dc,RANDOMFIELD3DC)( int *, int *, int *, int *, int *, int *, int*, int*, int*, 
                                              int*, double*, double*, double*, double*, double*, double*, int*, double*, int*, int* );

void F77_FUNC(perturbvelocity,PERTURBVELOCITY)( int *, int *, int *, int *, int *, int *, double*, 
						double*, double*, double*, double*, double*, double* );
void F77_FUNC(perturbvelocityc,PERTURBVELOCITYC)( int *, int *, int *, int *, int *, int *, double*, 
						  double*, double*, double*, double*, double* );
void F77_FUNC(checkanisomtrl,CHECKANISOMTRL)( int *, int *, int *, int *, int *, int *, double*, 
					      double*, double*, double*, double*, double* );
void F77_FUNC(computedtaniso,COMPUTEDTANISO)( int *, int *, int *, int *, int *, int *, double*,
					      double*, double*, double*, double* );
void F77_FUNC(computedtaniso2,COMPUTEDTANISO2)( int *, int *, int *, int *, int *, int *, double*,
					      double*, double*, double*, double* );
void F77_FUNC(computedtaniso2curv,COMPUTEDTANISO2CURV)( int *, int *, int *, int *, int *, int *, double*,
					      double*, double*, double*, double* );
  void anisomtrltocurvilinear( int*, int*, int*, int*, int*, int*, double*, double*, double* );
}

#define SQR(x) ((x)*(x))

//----------------------------------------------
void EW::setupRun( vector<Source*> & a_GlobalUniqueSources )
{
   if( mIsInitialized && proc_zero() )
      cout << " WARNING, calling setupRun twice " << endl;

  double time_start = MPI_Wtime();

// m_testing == true is one of the pointers to testing modes is assigned
// add other pointers to the list of testing modes as they get implemented
  m_testing = (m_twilight_forcing || m_point_source_test || m_lamb_test || m_rayleigh_wave_test || m_energy_test ); 

// tmp
  if (mVerbose && proc_zero() )
    cout << " *** Testing = " << m_testing << endl;
    
  if( mVerbose && proc_zero() )
  {
    cout << "  Using Bjorn's fast (parallel)" << " IO library" << endl;
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
    for( g=mNumberOfGrids-1; g>=0; g-- )
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

// Check that f.d. operators fit inside the domains.
  check_dimensions();

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
	else if( m_bcType[g][side] == bPeriodic )
	  side_plane( g, side, wind, m_ghost_points );
	else // for Dirichlet, super grid, and periodic conditions, we
	     // apply the forcing directly on the ghost points
	{
	   if( m_mesh_refinements && side < 4 )
	      side_plane( g, side, wind, m_ghost_points+1 );
	   else
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
//  if (topographyExists())
//  {
//    communicate_array( mJ, mNumberOfGrids-1 );
//    communicate_array( mMetric, mNumberOfGrids-1 );
    //    communicate_array( mQ, mNumberOfGrids-1 );
    //    communicate_array( mR, mNumberOfGrids-1 );
    //    communicate_array( mS, mNumberOfGrids-1 );

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

//  }

  MPI_Barrier(MPI_COMM_WORLD);

  //  string cachePath = mPath;
   
// tmp
  if ( mVerbose >= 3 )
   {
     int top = mNumberOfGrids-1;
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

//  string saved_path = mPath;

  int beginCycle = 1;

// Initialize IO
  create_output_directory( );

  if (proc_zero())
  {
    double lat[4],lon[4];
    computeGeographicCoord(0.0,           0.0,           lon[0], lat[0]);
    computeGeographicCoord(m_global_xmax, 0.0,           lon[1], lat[1]);
    computeGeographicCoord(m_global_xmax, m_global_ymax, lon[2], lat[2]);
    computeGeographicCoord(0.0,           m_global_ymax, lon[3], lat[3]);
// test inverse mapping too
    double xc[4], yc[4];
    for (int q=0; q<4; q++)
       computeCartesianCoord(xc[q], yc[q], lon[q], lat[q]);
    
    printf("Geographic and Cartesian coordinates of the corners of the computational grid:\n");
    for (int q=0; q<4; q++)
    {
//       printf("%i: Lon= %e, Lat=%e\n", q, lon[q], lat[q]);
       printf("%i: Lon= %e, Lat=%e, x=%e, y=%e\n", q, lon[q], lat[q], xc[q], yc[q]);
    }
    
    printf("\n");
  }
  
// set material properties
  if( m_anisotropic )
     set_anisotropic_materials();
  else
     set_materials();

// evaluate resolution
  double minvsoh;
  if( !m_anisotropic )
  {
     compute_minvsoverh( minvsoh );
     if (proc_zero())
     {
	printf("\n***** PPW = minVs/h/maxFrequency ********\n");
	for (int g=0; g<mNumberOfCartesianGrids; g++)
	{
	   printf("g=%i, h=%e, minVs/h=%g (Cartesian)\n", g, mGridSize[g], mMinVsOverH[g]);
	}
	if (topographyExists())
	{
	   int g = mNumberOfGrids-1;
	   printf("g=%i, h=%e, minVs/h=%g (curvilinear)\n", g, mGridSize[g], mMinVsOverH[g]);
	}
	printf("\n");
     }
  }
  
  assign_supergrid_damping_arrays();

// convert Qp and Qs to muVE, lambdaVE, and compute unrelaxed lambda, mu
  if( m_use_attenuation )
  {
     setup_attenuation_relaxation( minvsoh );
     if(  m_twilight_forcing )
	setup_viscoelastic_tw();
     else
	setup_viscoelastic();
  }
  
  if( mVerbose && proc_zero() )
    cout << "  Assigned material properties" << endl;

// compute time-step and number of time steps. 
// Note: SW4 always ends the simulation at mTmax, whether prefilter is enabled or not.
// This behavior is different from WPP
  if( m_anisotropic )
     computeDTanisotropic();
  else
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

// do source specific initialization such as prefiltering, scaling by mu, compute epicenter, etc (see below)
  preprocessSources( a_GlobalUniqueSources );
  
  double time_start_solve = MPI_Wtime();
  print_execution_time( time_start, time_start_solve, "start up phase" );
}

//-----------------------------------------------------------------------
void EW::preprocessSources( vector<Source*> & a_GlobalUniqueSources )
{
// This routine should be called once, after setupRun (can we include it in setupRun?)

// For non-testing modes, this routine performs the following tasks:
// 1: computes the epicenter, i.e. the location of the source with the earliest start time
// 2: if any source returns get_ShearModulusFactor() == true, it multiplies Mij by the shear modulus
// 3: evaluate the min and max z-coordinate of all sources
// 4: sets freq=1/dt in all sources with iDirac time-function
// 5: modifies the time function if pre-filtering is enabled *** MOVED TO Source:prepareTimeFunc() ***
// 6: saves a GMT file if requested

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
				   a_GlobalUniqueSources[0]->getName() == "C6SmoothBump" ||
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
	
	if ( !( (a_GlobalUniqueSources[0]->getName() == "VerySmoothBump" ||
		 a_GlobalUniqueSources[0]->getName() == "C6SmoothBump" ) &&
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
// normal source case starts here...

// find the epicenter, i.e., the location of the source with the lowest value of t0

// get the epicenter
      compute_epicenter( a_GlobalUniqueSources );
      
// Set up 'normal' sources for point_source_test, lamb_test, or standard seismic case.

// if the sources were defined by a rupture file, we need to multiply the Mij coefficients by mu (shear modulus)
      bool need_mu_corr=false;
      
      for( int i=0 ; i < a_GlobalUniqueSources.size() ; i++ )
	need_mu_corr = (need_mu_corr || a_GlobalUniqueSources[i]->get_CorrectForMu( ));

// must communicate need_mu_corr
      int mu_corr_global=0, mu_corr_loc= need_mu_corr? 1:0;
      
// take max over all procs to communicate 
      MPI_Allreduce( &mu_corr_loc, &mu_corr_global, 1, MPI_INT, MPI_MAX, m_cartesian_communicator);
      need_mu_corr=(bool) mu_corr_global;

// tmp
//      printf(" Proc #%i, sources needs correction for shear modulus: %s\n", getRank(), need_mu_corr? "TRUE":"FALSE");

      if (!mQuiet && mVerbose >= 3 && proc_zero() )
	printf(" Some sources needs correction for shear modulus: %s\n", need_mu_corr? "TRUE":"FALSE");

      if (need_mu_corr)
      {
	int nSources=a_GlobalUniqueSources.size();
	// if (proc_zero())
	//   printf("Number of sources: %i\n", nSources);

// allocate a temp array for the mu value at all source locations
	double *mu_source_loc = new double[nSources];
	double *mu_source_global = new double[nSources];
// initialize
	int s;
	for (s=0; s<nSources; s++)
	{
	  mu_source_loc[s]=-1.0;
	  mu_source_global[s]=-1.0;
	}
	
// fill in the values that are known to this processor
	for (s=0; s<nSources; s++)
	  if (a_GlobalUniqueSources[s]->myPoint())
	  {
	    int is=a_GlobalUniqueSources[s]->m_i0;
	    int js=a_GlobalUniqueSources[s]->m_j0;
	    int ks=a_GlobalUniqueSources[s]->m_k0;
	    int gs=a_GlobalUniqueSources[s]->m_grid;
	    
// tmp
	    mu_source_loc[s] = mMu[gs](is,js,ks); 
// printf("Proc #%i, source#%i, i=%i, j=%i, k=%i, g=%i, mu=%e\n", getRank(), s, is, js, ks, gs, mu_source_loc[s]);
	  }
// take max over all procs: communicate 
	MPI_Allreduce( mu_source_loc, mu_source_global, nSources, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator);

	if (!mQuiet && mVerbose >= 3 && proc_zero() )
	  printf(" Done communicating shear modulus to all procs\n");
// tmp
	// for (s=0; s<nSources; s++)
	//   printf("Proc #%i, source#%i, mu=%e\n", getRank(), s, mu_source_global[s]);

// scale all moments components
	double mu, mxx, mxy, mxz, myy, myz, mzz;
	for (s=0; s<nSources; s++)
	  if (a_GlobalUniqueSources[s]->get_CorrectForMu())
	  {
	    mu = mu_source_global[s];
	    a_GlobalUniqueSources[s]->getMoments( mxx, mxy, mxz, myy, myz, mzz);
	    a_GlobalUniqueSources[s]->setMoments( mu*mxx, mu*mxy, mu*mxz, mu*myy, mu*myz, mu*mzz);
// lower the flag
	    a_GlobalUniqueSources[s]->set_CorrectForMu(false);
	  }
	
// cleanup
	delete[] mu_source_loc;
	delete[] mu_source_global;
      } // end if need_mu_corr
      
      
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
      if (!mQuiet && mVerbose >= 1 && proc_zero() )
	printf(" Min source z-level: %e, max source z-level: %e\n", zMinGlobal, zMaxGlobal);

// Need to set the frequency to 1/dt for Dirac source
      for( int s=0 ; s  < a_GlobalUniqueSources.size(); s++ )
	 if( a_GlobalUniqueSources[s]->getTfunc() == iDirac )
	    a_GlobalUniqueSources[s]->setFrequency( 1.0/mDt );

      if (m_prefilter_sources)
      {
// tell the filter about the time step and compute the second order sections
	m_filter_ptr->computeSOS( mDt );

// output details about the filter
	if (mVerbose>=3 && proc_zero() )
           cout << *m_filter_ptr;

	if (!getQuiet() && proc_zero() )
        {
           printf("Filter precursor = %e\n", m_filter_ptr->estimatePrecursor() );
        }
        
      }
      
// // Modify the time functions if prefiltering is enabled
	
// // 1. Make sure the smallest time offset is at least t0_min + (timeFcn dependent offset for centered fcn's)
// 	double dt0 = 0;
// 	double dt0loc, dt0max, t0_min;
// 	t0_min = m_filter_ptr->estimatePrecursor();
// // tmp
// 	if (!mQuiet && proc_zero() )
// 	  printf("Filter precursor = %e\n", t0_min);
	
// // old estimate for 2-pole low-pass Butterworth
// //	t0_min = 4./m_filter_ptr->get_corner_freq2();
	
// 	for( int s=0; s < a_GlobalUniqueSources.size(); s++ ) 
// 	{
// 	  dt0loc = a_GlobalUniqueSources[s]->compute_t0_increase( t0_min );
	  
// 	  if( dt0loc > dt0 )
// 	    dt0 = dt0loc;
// 	}
// // compute global max over all processors
// 	MPI_Allreduce( &dt0, &dt0max, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator);

// // dt0max is the maxima over all dt0loc in all processors. 
// // If it is positive, the t0 field in all source commands should be incremented 
// // by at least this amount. Otherwise, there might be significant artifacts from 
// // a sudden start of some source.
// 	if (dt0max > 0.)
// 	{
// // Don't mess with t0.
// // Instead, warn the user of potential transients due to unsmooth start
// //	  for( int s=0; s < a_GlobalUniqueSources.size(); s++ ) 
// //	    a_GlobalUniqueSources[s]->adjust_t0( dt0max );
// 	  if ( !mQuiet && proc_zero() )
// 	    printf("\n*** WARNING: the 2 pass prefilter has an estimated precursor of length %e s\n"
// 		   "*** To avoid artifacts due to sudden startup, increase t0 in all source commands by at least %e\n\n",
// 		   t0_min, dt0max);
// 	}

// // Do the filtering
// 	for( int s=0; s < a_GlobalUniqueSources.size(); s++ ) 
//            a_GlobalUniqueSources[s]->filter_timefunc( m_filter_ptr, mTstart, mDt, mNumberOfTimeSteps );
//       }

// // TODO: check that t0 is large enough even when prefilter is NOT used

      if (proc_zero())
	saveGMTFile( a_GlobalUniqueSources );

    } // end normal seismic setup

  } // end if ( !m_twilight_forcing && !m_energy_test && !m_rayleigh_wave_test ) 
   
  mSourcesOK = true;


} // end preprocessSources

//-----------------------------------------------------------------------
void EW::compute_epicenter( vector<Source*> & a_GlobalUniqueSources ) 
{
  // To find out which event goes first, we need to query all sources
  double earliestTime=0.;
  double epiLat = 0.0, epiLon = 0.0, epiDepth = 0.0;
      
  if (a_GlobalUniqueSources.size() == 0)
  {
    // Put a location for a source at the origin
    computeGeographicCoord(0., 0., epiLat, epiLon);
  }
  else
  {
    // Find the source with the earliest time...
    const Source* firstSource = a_GlobalUniqueSources[0];
    earliestTime = firstSource->getOffset();
      
    for (unsigned int i = 1; i < a_GlobalUniqueSources.size(); ++i)
    {
      if (a_GlobalUniqueSources[i]->getOffset() < earliestTime)
      {
	firstSource = a_GlobalUniqueSources[i];
	earliestTime = firstSource->getOffset();
      }
    }
    
    computeGeographicCoord(firstSource->getX0(), firstSource->getY0(), epiLon, epiLat );
    epiDepth = firstSource->getDepth(); // corrected for topography!
  }

  set_epicenter(epiLat, epiLon, epiDepth, earliestTime);
  
}

//-----------------------------------------------------------------------
void EW::setupSBPCoeff()
{
  double gh2; // this coefficient is also stored in m_ghcof[0]
  if (mVerbose >=1 && m_myRank == 0)
    cout << "Setting up SBP boundary stencils" << endl;
  if( m_croutines )
     GetStencilCoefficients( m_acof, m_ghcof, m_bope, m_sbop );
  else
  {
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
// If material surfaces (Ifiles) are defined, sort them wrt their IDs
     if(m_materials.size() > 0)
     {
// Sort the m_materials vector in increasing id order...but only if there are more than one element.
	if (m_materials.size() > 1) 
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
       if (m_mtrlblocks[b]->coversAllPoints())
	  lastAllCoveringBlock=b;
// tmp
    if (proc_zero())
    {
      if (lastAllCoveringBlock == 0)
	cout << "Considering all material blocks" << endl;
      else
	cout << "Only considering material blocks with index >= " << lastAllCoveringBlock << endl;
    } // end if proc_zero()
    for( unsigned int b = lastAllCoveringBlock ; b < m_mtrlblocks.size() ; b++ )
       m_mtrlblocks[b]->set_material_properties(mRho, mMu, mLambda, mQs, mQp); 

//   bool linearExtrapolation=false;
// note that material thresholding for vs and vp happens further down in this procedure
//   extrapolateInZ( mRho[g],    false, 0., linearExtrapolation ); 
//   extrapolateInZ( mLambda[g], false, 0., linearExtrapolation ); 
//   extrapolateInZ( mMu[g],     false, 0., linearExtrapolation );
//   if( m_use_attenuation )
//   {
//     extrapolateInZ(mQs[g], false, 0., linearExtrapolation); 
//     extrapolateInZ(mQp[g], false, 0., linearExtrapolation); 
//   }

// extrapolate to define material properties above the free surface (topography)
   g = mNumberOfGrids-1;
   extrapolateInZ( g, mRho[g],    true, false ); 
   extrapolateInZ( g, mLambda[g], true, false ); 
   extrapolateInZ( g, mMu[g],     true, false );
   if( m_use_attenuation )
   {
      extrapolateInZ(g, mQs[g], true, false );
      extrapolateInZ(g, mQp[g], true, false );
   }

// Extrapolate to ghost points at bottom of domain, if necessary
   g = 0;
   extrapolateInZ( g, mRho[g],    false, true ); 
   extrapolateInZ( g, mLambda[g], false, true ); 
   extrapolateInZ( g, mMu[g],     false, true );
   if( m_use_attenuation )
   {
      extrapolateInZ(g, mQs[g], false, true );
      extrapolateInZ(g, mQp[g], false, true );
   }


// extrapolate material properties to mesh refinement boundaries (e.g. for doing the LOH cases more accurately)
    if (!mQuiet && proc_zero() && mVerbose>=3)
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

	  if (!mQuiet && proc_zero() && mVerbose>=3)
	    printf("setMaterials> top extrapol, g=%i, kFrom=%i, kStart=%i\n", g, kFrom, m_kStart[g]);

	  for (int k = m_kStart[g]; k < kFrom; ++k)
	    for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
	      for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
	      {
		mRho[g](i,j,k)    = mRho[g](i,j,kFrom);
		mMu[g](i,j,k)     = mMu[g](i,j,kFrom);
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

	  if (!mQuiet && proc_zero() && mVerbose>=3)
	    printf("setMaterials> bottom extrapol, g=%i, kFrom=%i, kEnd=%i\n", g, kFrom, m_kEnd[g]);

	  for (int k = kFrom+1; k <= m_kEnd[g]; ++k)
	    for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
	      for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
	      {
		mRho[g](i,j,k)    = mRho[g](i,j,kFrom);
		mMu[g](i,j,k)     = mMu[g](i,j,kFrom);
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
// Extrapolate to ghost points in x and y, if they were not set by the previous routines.
//    cout << "min rho before " << mRho[0].minimum() << endl;
    extrapolateInXY( mRho );
    extrapolateInXY( mMu );
    extrapolateInXY( mLambda );
    if( m_use_attenuation )
    {
      extrapolateInXY(mQs);
      extrapolateInXY(mQp);
    }
    //    cout << "min rho after " << mRho[0].minimum() << endl;
    if( m_use_attenuation && m_qmultiplier != 1 )
    {
       for (int g=0; g<mNumberOfGrids; g++)
	  for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
	     for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
		for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
		{
		   mQs[g](i,j,k) = m_qmultiplier*mQs[g](i,j,k);
		   mQp[g](i,j,k) = m_qmultiplier*mQp[g](i,j,k);
		}
    }

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
    
// add random perturbation
    if( m_randomize )
       perturb_velocities( mMu, mLambda );

    convert_material_to_mulambda( );
    
    check_for_nan( mMu, 1,"mu ");       
    check_for_nan( mLambda, 1,"lambda ");       
    check_for_nan( mRho, 1,"rho ");       

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
// Need to communicate across material boundaries, why ?
//	  communicate_array( mRho[g], g );
//	  communicate_array( mMu[g], g );
//	  communicate_array( mLambda[g], g );
      }
      if( topographyExists() )
      {
 	  g = mNumberOfGrids-1;
	  rho_ptr = mRho[g].c_ptr();
	  mu_ptr  = mMu[g].c_ptr();
	  la_ptr  = mLambda[g].c_ptr();
	  ifirst = m_iStart[g];
	  ilast  = m_iEnd[g];
	  jfirst = m_jStart[g];
	  jlast  = m_jEnd[g];
	  kfirst = m_kStart[g];
	  klast  = m_kEnd[g];
          double* x_ptr = mX.c_ptr();
          double* y_ptr = mY.c_ptr();
          double* z_ptr = mZ.c_ptr();
	  omm = m_twilight_forcing->m_momega;
	  phm = m_twilight_forcing->m_mphase;
	  amprho = m_twilight_forcing->m_amprho;
	  ampmu = m_twilight_forcing->m_ampmu;
	  ampla = m_twilight_forcing->m_amplambda;
	  F77_FUNC(exactmatfortc,EXACTMATFORTC)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					      &klast, rho_ptr, mu_ptr, la_ptr, &omm, &phm, 
					&amprho, &ampmu, &ampla, x_ptr, y_ptr, z_ptr );
//   // Need this for Energy testing, random material will not agree on processor boundaries.
// 	  communicate_array( mRho[g], g );
// 	  communicate_array( mMu[g], g );
// 	  communicate_array( mLambda[g], g );
      }
      
  } // end material set by forcing mode (for testing)
  else if ( m_point_source_test )
  {
     for (g=0; g<mNumberOfGrids; g++)
     {
	mRho[g].set_value( m_point_source_test->m_rho );
	mMu[g].set_value( m_point_source_test->m_mu );
	mLambda[g].set_value( m_point_source_test->m_lambda );
     }
  }
  else if ( m_lamb_test )
  {
     for (g=0; g<mNumberOfGrids; g++)
     {
	mRho[g].set_value( m_lamb_test->m_rho );
	mMu[g].set_value( m_lamb_test->m_mu );
	mLambda[g].set_value( m_lamb_test->m_lambda );
     }
  }
  else if ( m_rayleigh_wave_test )
  {
     for (g=0; g<mNumberOfGrids; g++)
     {
	mRho[g].set_value( m_rayleigh_wave_test->m_rho );
	mMu[g].set_value( m_rayleigh_wave_test->m_mu );
	mLambda[g].set_value( m_rayleigh_wave_test->m_lambda );
     }
  }
  else if ( m_energy_test )
  {
     double cpocs = m_energy_test->m_cpcsratio;
     for (g=0; g<mNumberOfGrids; g++)
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
	// Communication needed here. Because of random material data,
	// processor overlap points will otherwise have different values
	// in different processors.
        material_ic( mRho );
        material_ic( mMu );
	material_ic( mLambda );
	communicate_array( mRho[g], g );
	communicate_array( mMu[g], g );
	communicate_array( mLambda[g], g );
     }
  }
  check_materials( );

} // end EW::set_materials()


//-----------------------------------------------------------------------
void EW::set_anisotropic_materials()
{
   if (!m_testing)
   {      
      int lastAllCoveringBlock=0;
      for( unsigned int b = 0 ; b < m_anisotropic_mtrlblocks.size() ; b++ )
         if (m_anisotropic_mtrlblocks[b]->coversAllPoints())
            lastAllCoveringBlock=b;

      for( unsigned int b = lastAllCoveringBlock ; b < m_anisotropic_mtrlblocks.size() ; b++ )
         m_anisotropic_mtrlblocks[b]->set_material_properties( mRho, mC );

      if (proc_zero())
      {
         if (lastAllCoveringBlock == 0)
            cout << "Considering all material blocks" << endl;
         else
            cout << "Only considering material blocks with index >= " << lastAllCoveringBlock << endl;
      } // end if proc_zero()

      int g = mNumberOfGrids-1;
      extrapolateInZ( g, mRho[g],    true, false ); 
      extrapolateInZvector( g, mC[g],    true, false ); 

      g = 0;
      extrapolateInZ( g, mRho[g],    false, true ); 
      extrapolateInZvector( g, mC[g],    false, true ); 
      if (mMaterialExtrapolate > 0 && mNumberOfCartesianGrids > 1)
      {
         int kFrom;
         for (g=0; g<mNumberOfCartesianGrids; g++)
         {
            if (g < mNumberOfCartesianGrids-1) // extrapolate to top
            {
               kFrom = m_kStart[g]+mMaterialExtrapolate;
               if (!mQuiet && proc_zero() && mVerbose>=3)
                  printf("setMaterials> top extrapol, g=%i, kFrom=%i, kStart=%i\n", g, kFrom, m_kStart[g]);
               for (int k = m_kStart[g]; k < kFrom; ++k)
                  for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
                     for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
                     {
                        mRho[g](i,j,k)    = mRho[g](i,j,kFrom);
                        for( int m=1 ; m <= 21 ; m++ )
                           mC[g](m,i,j,k)    = mC[g](m,i,j,kFrom);
                     }
            } // end extrapolat to top
            if (g > 0) // extrapolate to bottom
            {
               kFrom = m_kEnd[g]-mMaterialExtrapolate;
               if (!mQuiet && proc_zero() && mVerbose>=3)
                  printf("setMaterials> bottom extrapol, g=%i, kFrom=%i, kEnd=%i\n", g, kFrom, m_kEnd[g]);
               for (int k = kFrom+1; k <= m_kEnd[g]; ++k)
                  for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
                     for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
                     {
                        mRho[g](i,j,k)    = mRho[g](i,j,kFrom);
                        for( int m=1 ; m <= 21 ; m++ )
                           mC[g](m,i,j,k)    = mC[g](m,i,j,kFrom);
                     }
            } 
         }
      } 
      extrapolateInXY( mRho );
      extrapolateInXYvector( mC );
      check_anisotropic_material( mRho, mC );
      if( m_energy_test )
      {
         material_ic( mRho );
         material_ic( mC );
         communicate_array( mRho[g], g );
         communicate_array( mC[g], g );
      }
      if( topographyExists() )
      {
         int g=mNumberOfGrids-1;
         anisomtrltocurvilinear( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
				 mMetric.c_ptr(), mC[g].c_ptr(), mCcurv.c_ptr() );
      }
   }// end if !m_testing, i.e., not Twilight
   else if (m_twilight_forcing) 
   {
// tmp
      if (proc_zero())
         cout << "******************************" << endl
              << " ASSIGNING ANISOTROPIC TWILIGHT MATERIALS " << endl
              << "******************************" << endl;

// For some forcings (such as twilight forcing) the material is set here.
      double xP, yP, zP;
      
      int ifirst, ilast, jfirst, jlast, kfirst, klast, g;
      double *cm_ptr, *rho_ptr, h, zmin, omm, phm, amprho, ampmu, ampla;
      double phc[21]; // move these angles to the EW class

      // need to store all the phase angle constants somewhere
     phc[0]=0;
      for (int i=0; i<21; i++)
         phc[i] = i*10*M_PI/180;
	
      for (g=0; g<mNumberOfCartesianGrids; g++)
      {
         rho_ptr = mRho[g].c_ptr();
         cm_ptr = mC[g].c_ptr();

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

// setup density (rho)
// setup rho and stiffness matrix         
        tw_ani_stiff(ifirst, ilast, jfirst, jlast, kfirst, klast, h, zmin,
                     omm, phm, amprho, rho_ptr, phc, cm_ptr);
         
// also need rho
      }
      if( topographyExists() )
      {
         g = mNumberOfGrids-1;
         rho_ptr = mRho[g].c_ptr();
         cm_ptr = mC[g].c_ptr();

         double* x_ptr = mX.c_ptr();
         double* y_ptr = mY.c_ptr();
         double* z_ptr = mZ.c_ptr();

         ifirst = m_iStart[g];
         ilast  = m_iEnd[g];
         jfirst = m_jStart[g];
         jlast  = m_jEnd[g];
         kfirst = m_kStart[g];
         klast  = m_kEnd[g];

         omm = m_twilight_forcing->m_momega;
         phm = m_twilight_forcing->m_mphase;
         amprho = m_twilight_forcing->m_amprho;

         if (proc_zero() )
            printf("set_anisotropic_mat> before tw_ani_curvi_stiff\n");
         tw_ani_curvi_stiff(ifirst, ilast, jfirst, jlast, kfirst, klast, x_ptr, y_ptr, z_ptr,
                            omm, phm, amprho, rho_ptr, phc, cm_ptr);
         if (proc_zero() )
            printf("set_anisotropic_mat> after tw_ani_curvi_stiff\n");

      }

// make sure the stiffness matrix is positive definite
      check_anisotropic_material( mRho, mC );
      if( topographyExists() )
      {
         int g=mNumberOfGrids-1;
         anisomtrltocurvilinear( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
				 mMetric.c_ptr(), mC[g].c_ptr(), mCcurv.c_ptr() );
      }
      
   } // end if m_twilight
   else
   {
      CHECK_INPUT(false, "Error:  the only test for an anisotropic material is twilight" << endl);
   }
   
} // end set_anisotropic_materials()


//-----------------------------------------------------------------------
void EW::check_anisotropic_material( vector<Sarray>& rho, vector<Sarray>& c )
{
   double rhomin=1e38, rhomax=-1e38;
   double eigmin=1e38, eigmax=-1e38;
   double rhominloc, rhomaxloc, eigminloc, eigmaxloc;
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      double* rho_ptr = rho[g].c_ptr();
      double* c_ptr = c[g].c_ptr();
      F77_FUNC(checkanisomtrl,CHECKANISOMTRL)(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
					      &m_kStart[g], &m_kEnd[g],
					      rho_ptr, c_ptr,
					      &rhominloc, &rhomaxloc, &eigminloc, &eigmaxloc );
      if( rhominloc < rhomin )
	 rhomin = rhominloc;
      if( rhomaxloc > rhomax )
	 rhomax = rhomaxloc;
      if( eigminloc < eigmin )
	 eigmin = eigminloc;
      if( eigmaxloc > eigmax )
	 eigmax = eigmaxloc;
   }
   rhominloc = rhomin;
   rhomaxloc = rhomax;
   eigminloc = eigmin;
   eigmaxloc = eigmax;
   MPI_Allreduce( &rhominloc, &rhomin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &rhomaxloc, &rhomax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &eigminloc, &eigmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &eigmaxloc, &eigmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

   if( proc_zero() )
   {
      cout << " Material properties " << endl;
      cout <<  rhomin <<  " <=  Density <= " << rhomax  << endl;
      cout <<  eigmin <<  " <=  Eig(c)  <= " << eigmax << endl;
      CHECK_INPUT(eigmin > 0, "ERROR: tensor of elastic coefficients is not positive definite. The problem is not well posed" << endl);
   }
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
  if (!mQuiet && mVerbose >= 1 && proc_zero())
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
           double jinv = 1/mJ(i,j,k);
// A11
//	   Amat[0] = -4.*(SQR(mQ(1,i,j,k))*la2mu + SQR(mQ(2,i,j,k))*mu + SQR(mQ(3,i,j,k))*mu 
//		        + SQR(mR(1,i,j,k))*la2mu + SQR(mR(2,i,j,k))*mu + SQR(mR(3,i,j,k))*mu
//			 + SQR(mS(1,i,j,k))*la2mu + SQR(mS(2,i,j,k))*mu + SQR(mS(3,i,j,k))*mu);
           Amat[0] = -4*(SQR(mMetric(1,i,j,k))*la2mu + SQR(mMetric(1,i,j,k))*mu + 
			 SQR(mMetric(2,i,j,k))*la2mu + SQR(mMetric(3,i,j,k))*mu + SQR(mMetric(4,i,j,k))*mu)*jinv;
// A21 = A12
//	   Amat[1] = -4.*(mQ(1,i,j,k)*mQ(2,i,j,k) + mR(1,i,j,k)*mR(2,i,j,k) + mS(1,i,j,k)*mS(2,i,j,k))*(mu+la);
	   Amat[1] = -4.*mMetric(2,i,j,k)*mMetric(3,i,j,k)*(mu+la)*jinv;
// A31 = A13	   
//	   Amat[2] = -4.*(mQ(1,i,j,k)*mQ(3,i,j,k) + mR(1,i,j,k)*mR(3,i,j,k) + mS(1,i,j,k)*mS(3,i,j,k))*(mu+la);
	   Amat[2] = -4.*mMetric(2,i,j,k)*mMetric(4,i,j,k)*(mu+la)*jinv;
// A22	   
	   Amat[3] = -4.*(SQR(mMetric(1,i,j,k))*mu + SQR(mMetric(1,i,j,k))*la2mu +
		        + SQR(mMetric(2,i,j,k))*mu + SQR(mMetric(3,i,j,k))*la2mu + SQR(mMetric(4,i,j,k))*mu)*jinv;
// A32 = A23
//	   Amat[4] = -4.*(mQ(2,i,j,k)*mQ(3,i,j,k) + mR(2,i,j,k)*mR(3,i,j,k) + mS(2,i,j,k)*mS(3,i,j,k))*(mu+la);
	   Amat[4] = -4.*mMetric(3,i,j,k)*mMetric(4,i,j,k)*(mu+la)*jinv;
// A33
	   Amat[5] = -4.*(SQR(mMetric(1,i,j,k))*mu + SQR(mMetric(1,i,j,k))*mu
			+ SQR(mMetric(2,i,j,k))*mu + SQR(mMetric(3,i,j,k))*mu + SQR(mMetric(4,i,j,k))*la2mu)*jinv;
// calculate eigenvalues of symmetric matrix
	   F77_FUNC(dspev,DSPEV)(JOBZ, UPLO, N, Amat, W, Z, LDZ, WORK, INFO);

	   if (INFO != 0)
	   {
	     printf("ERROR: computeDT: dspev returned INFO = %i for grid point (%i, %i, %i)\n", INFO, i, j, k);
	     printf("lambda = %e, mu = %e\n", la, mu);
             printf("Jacobian = %15.7g \n",mJ(i,j,k));
	     printf("Matrix = \n");
	     printf(" %15.7g  %15.7g %15.7g \n",Amat[0],Amat[1],Amat[2]);
	     printf(" %15.7g  %15.7g %15.7g \n",Amat[1],Amat[3],Amat[4]);
	     printf(" %15.7g  %15.7g %15.7g \n",Amat[2],Amat[4],Amat[5]);
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
//           double asym = SQR(mMetric(2,i,j,k))+SQR(mMetric(3,i,j,k))+SQR(mMetric(4,i,j,k));
//           double m1 = SQR(mMetric(1,i,j,k));
//           double eg1 = m1*(3*mu+la)+mu*asym;
//	   double eg2 = (mu+la)*(0.5*(m1+asym) +0.5*sqrt( (m1+asym)*(m1+asym)-4*SQR(mMetric(4,i,j,k))*m1)) + mu*(2*m1+asym);
//           eg1 = eg1*4*jinv;
//           eg2 = eg2*4*jinv;
//           cout << "eig solver " << W[0] << " analytic values " << eg1 << " and " << eg2 << endl;

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

    //   cout << "cfl = " << mCFL << endl;
    //   cout << "dtloc = " << mDt << endl;

// global minima for curvilinear grid
    if (topographyExists())
    {
      double dtCurvGlobal;
      MPI_Allreduce( &dtCurv, &dtCurvGlobal, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
      if (!mQuiet && mVerbose >= 1 && proc_zero())
      {
	printf("INFO: Smallest stable time step for curvilinear grid only: %e\n", dtCurvGlobal);
      }
    }

    if (!mQuiet && mVerbose >= 1 && proc_zero())
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
void EW::computeDTanisotropic()
{
   if (!mQuiet && mVerbose >= 1 && proc_zero())
   {
      printf("*** computing the time step ***\n");
   }
   double dtproc=1.e38;
   for (int g=0; g<mNumberOfCartesianGrids; g++)
   {
      double* rho_ptr = mRho[g].c_ptr();
      double* c_ptr = mC[g].c_ptr();
      double dtgrid;
      F77_FUNC(computedtaniso2,COMPUTEDTANISO2)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
						 &m_kStart[g], &m_kEnd[g],
						 rho_ptr, c_ptr, &mCFL, &mGridSize[g], &dtgrid );
      if( dtgrid < dtproc )
	 dtproc = dtgrid;
   }
   if( topographyExists() )
   {
      int g=mNumberOfGrids-1;
      double* rho_ptr = mRho[g].c_ptr();
      double* c_ptr = mCcurv.c_ptr();
      double* jac_ptr = mJ.c_ptr();
      double dtgrid;
      F77_FUNC(computedtaniso2curv,COMPUTEDTANISO2CURV)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
							 &m_kStart[g], &m_kEnd[g],
							 rho_ptr, c_ptr, jac_ptr, &mCFL, &dtgrid );
      if( dtgrid < dtproc )
	 dtproc = dtgrid;
   }
// compute the global minima
    MPI_Allreduce( &dtproc, &mDt, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
    if (!mQuiet && mVerbose >= 1 && proc_zero())
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
   
   //   string pathTemp(path.begin(), path.end()); 
   string pathTemp = path;
   //-----------------------------------------------------------------
   // Recursively call stat and then mkdir on each sub-directory in 'path'
   //-----------------------------------------------------------------
   string sep = "/";

   char * pathtemparg = new char[pathTemp.length()+1];
   strcpy(pathtemparg,pathTemp.c_str());
   char* token = strtok( pathtemparg, sep.c_str() );
//   char* token = strtok(const_cast<char*>(pathTemp.c_str()), sep.c_str());

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
	  delete[] pathtemparg;
	  return -1;
	}
	
      }
      else
      {
//	cerr << "stat() returned an error code." << endl;
	if (errno == EACCES)
	{
	  cerr << "Error: **Search permission is denied for one of the directories in the path prefix of " << pathsofar.str() << endl;
	  delete[] pathtemparg;
	  return -1;
	}
	else if (errno == ENOTDIR)
	{
	  cerr << "Error: **A component of the path '" <<  pathsofar.str() << "' is not a directory. " << endl;
	  delete[] pathtemparg;
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
	delete[] pathtemparg;
	return -1;
      }
      else
      {
	if (mVerbose >=2) cout << "mkdir() returned successfully." << endl;

// are there more directories to be made?
	token = strtok(NULL, sep.c_str());
      }
   }
   delete[] pathtemparg;
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
  vector<double> sg_width(mNumberOfCartesianGrids);
  for( int g=0 ; g < mNumberOfCartesianGrids ; g++ )
  {
     if( m_use_sg_width )
	sg_width[g] = m_supergrid_width;
     else
	sg_width[g] = m_sg_gp_thickness*mGridSize[g];

     m_supergrid_taper_x[g].define_taper( (mbcGlobalType[0] == bSuperGrid), 0.0,
					  (mbcGlobalType[1] == bSuperGrid), m_global_xmax, 
					   sg_width[g] );
     m_supergrid_taper_y[g].define_taper( (mbcGlobalType[2] == bSuperGrid), 0.0,
					  (mbcGlobalType[3] == bSuperGrid), m_global_ymax, 
					   sg_width[g] );
  }
  if( topographyExists() )
  {
     int g=mNumberOfGrids-1;
     m_supergrid_taper_x[g].define_taper( (mbcGlobalType[0] == bSuperGrid), 0.0,
					  (mbcGlobalType[1] == bSuperGrid), m_global_xmax, 
					  sg_width[gTop] );
     m_supergrid_taper_y[g].define_taper( (mbcGlobalType[2] == bSuperGrid), 0.0,
					  (mbcGlobalType[3] == bSuperGrid), m_global_ymax, 
					  sg_width[gTop] );
  }
  if( mNumberOfGrids == 1 )
  {
     m_supergrid_taper_z[0].define_taper( !topographyExists() && (mbcGlobalType[4] == bSuperGrid), 0.0,
					  (mbcGlobalType[5] == bSuperGrid), m_global_zmax,
					  sg_width[0] );
  }
  else
  {
     m_supergrid_taper_z[mNumberOfGrids-1].define_taper( !topographyExists() && (mbcGlobalType[4] == bSuperGrid), 0.0,
							 false, m_global_zmax, sg_width[gTop] );
     m_supergrid_taper_z[0].define_taper( false, 0.0, mbcGlobalType[5]==bSuperGrid, m_global_zmax,
					  sg_width[0] );
     for( int g=1 ; g < mNumberOfGrids-1 ; g++ )
	m_supergrid_taper_z[g].define_taper( false, 0.0, false, 0.0, sg_width[g] );
  }
  
  for( int g=0 ; g < mNumberOfGrids ; g++ )
  {
   // Add one to thickness to allow two layers of internal ghost points
     int sgpts = m_sg_gp_thickness;
     if( m_use_sg_width )
	sgpts = m_supergrid_width/mGridSize[g];
     int imin = 1+sgpts, imax = m_global_nx[g]-sgpts, jmin=1+sgpts, jmax=m_global_ny[g]-sgpts;
     int kmax=m_global_nz[g]-sgpts;

     // Only grid 0 has super grid boundary at the bottom
     if( g > 0 )
	kmax = m_global_nz[g];

     //     cout << "Active region for backward solver: " << imin+1 << " " << imax-1 << " " << jmin+1 << " " << jmax-1
     //	  << " " << 1+1 << " " << kmax-1 << endl;

     m_iStartActGlobal[g] = m_iStartAct[g] = imin+1;
     m_iEndActGlobal[g]   = m_iEndAct[g]   = imax-1;
     m_jStartActGlobal[g] = m_jStartAct[g] = jmin+1;
     m_jEndActGlobal[g]   = m_jEndAct[g]   = jmax-1;
     m_kStartActGlobal[g] = m_kStartAct[g] = 1;
     m_kEndActGlobal[g]   = m_kEndAct[g]   = kmax-1;

     if( m_iStartAct[g] < m_iStart[g] )
	m_iStartAct[g] = m_iStart[g];
     if( m_jStartAct[g] < m_jStart[g] )
	m_jStartAct[g] = m_jStart[g];
     if( m_iEndAct[g] > m_iEnd[g] )
	m_iEndAct[g] = m_iEnd[g];
     if( m_jEndAct[g] > m_jEnd[g] )
	m_jEndAct[g] = m_jEnd[g];

     // If empty, set dimensions so that imax-imin+1=0, to avoid negative element count.
     if( m_iStartAct[g] > m_iEndAct[g] )
	m_iStartAct[g] = m_iEndAct[g]+1;
     if( m_jStartAct[g] > m_jEndAct[g] )
	m_jStartAct[g] = m_jEndAct[g]+1;
     if( m_kStartAct[g] > m_kEndAct[g] )
	m_kStartAct[g] = m_kEndAct[g]+1;
  }
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

  m_sg_str_x.resize(mNumberOfGrids);
  m_sg_str_y.resize(mNumberOfGrids);
  m_sg_str_z.resize(mNumberOfGrids);

// new corner taper functions to reduce strength of damping near the edges and corners
  m_sg_corner_x.resize(mNumberOfGrids);
  m_sg_corner_y.resize(mNumberOfGrids);
  m_sg_corner_z.resize(mNumberOfGrids);
  
// allocate storage for 1-D damping coefficients on each grid
  for( g=0 ; g<mNumberOfGrids; g++) 
  {
    m_sg_dc_x[g]  = new double[m_iEnd[g]-m_iStart[g]+1];
    m_sg_dc_y[g]  = new double[m_jEnd[g]-m_jStart[g]+1];
    m_sg_dc_z[g]  = new double[m_kEnd[g]-m_kStart[g]+1];

    m_sg_str_x[g] = new double[m_iEnd[g]-m_iStart[g]+1];
    m_sg_str_y[g] = new double[m_jEnd[g]-m_jStart[g]+1];
    m_sg_str_z[g] = new double[m_kEnd[g]-m_kStart[g]+1];

// new corner taper functions to reduce strength of damping near the edges and corners
    m_sg_corner_x[g] = new double[m_iEnd[g]-m_iStart[g]+1];
    m_sg_corner_y[g] = new double[m_jEnd[g]-m_jStart[g]+1];
    m_sg_corner_z[g] = new double[m_kEnd[g]-m_kStart[g]+1];
  }

#define dcx(i,g) (m_sg_dc_x[g])[i-m_iStart[g]]
#define dcy(j,g) (m_sg_dc_y[g])[j-m_jStart[g]]
#define dcz(k,g) (m_sg_dc_z[g])[k-m_kStart[g]]

#define strx(i,g) (m_sg_str_x[g])[i-m_iStart[g]]
#define stry(j,g) (m_sg_str_y[g])[j-m_jStart[g]]
#define strz(k,g) (m_sg_str_z[g])[k-m_kStart[g]]

#define cornerx(i,g) (m_sg_corner_x[g])[i-m_iStart[g]]
#define cornery(j,g) (m_sg_corner_y[g])[j-m_jStart[g]]
#define cornerz(k,g) (m_sg_corner_z[g])[k-m_kStart[g]]

  //  topCartesian = mNumberOfCartesianGrids-1;
// Note: compared to WPP2, we don't need to center the damping coefficients on the half-point anymore,
// because the damping term is now 4th order: D+D-( a(x) D+D- ut(x) )

  topCartesian = mNumberOfCartesianGrids-1;
  if( m_use_supergrid )
  {
// tmp
//    printf("SG: using supergrid!\n");
    
     if( m_twilight_forcing )
     {
// tmp
//       printf("SG: twilight setup!\n");

	for( g=0 ; g<mNumberOfGrids; g++)  
	{
	   for( i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	   {
	      x = (i-1)*mGridSize[g];
	      dcx(i,g)  = 0;
	      cornerx(i,g)  = 1;
	      strx(i,g) = m_supergrid_taper_x[g].tw_stretching(x);
	   }
	   for( j = m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	   {
	      y = (j-1)*mGridSize[g];
	      dcy(j,g)  = 0;
	      cornery(j,g)  = 1;
	      stry(j,g) = m_supergrid_taper_y[g].tw_stretching(y);
	   }
	   if ( g > topCartesian || (0 < g && g < mNumberOfGrids-1) ) // curvilinear grid or refinement grid.
	   {
	      for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	      {
		 dcz(k,g)  = 0.;
		 cornerz(k,g)  = 1.;
		 strz(k,g) = 1;
	      }
	   }
	   else
	   {
	      for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	      {
		 z = m_zmin[g] + (k-1)*mGridSize[g];
		 dcz(k,g)  = 0;
		 cornerz(k,g) = 1.;
		 strz(k,g) = m_supergrid_taper_z[g].tw_stretching(z);
	      }
	   }
	}
     } // end if (m_twilight_forcing) ...
     else
     { // standard case starts here
// tmp
//       printf("SG: standard case!\n");
	for( g=0 ; g<mNumberOfGrids; g++)  
	{
	   for( i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	   {
	      x = (i-1)*mGridSize[g];
	      dcx(i,g)  = m_supergrid_taper_x[g].dampingCoeff(x);
	      strx(i,g) = m_supergrid_taper_x[g].stretching(x);
	      cornerx(i,g)  = m_supergrid_taper_x[g].cornerTaper(x);
	   }
	   for( j = m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	   {
	      y = (j-1)*mGridSize[g];
	      dcy(j,g)  = m_supergrid_taper_y[g].dampingCoeff(y);
	      stry(j,g) = m_supergrid_taper_y[g].stretching(y);
	      cornery(j,g)  = m_supergrid_taper_y[g].cornerTaper(y);
	   }
	   if (g > topCartesian || (0 < g && g < mNumberOfGrids-1)  ) // Curvilinear or refinement grid
	   {
// No supergrid damping in the vertical (k-) direction on a curvilinear or refinement grid.
	      for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	      {
		 dcz(k,g) = 0.;
		 strz(k,g) = 1;
		 cornerz(k,g) = 1.;
	      }
	   }
	   else
	   {
	      for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	      {
		z = m_zmin[g] + (k-1)*mGridSize[g];
		dcz(k,g)  = m_supergrid_taper_z[g].dampingCoeff(z);
		strz(k,g) = m_supergrid_taper_z[g].stretching(z);
		cornerz(k,g) = m_supergrid_taper_z[g].cornerTaper(z);
	      }
	   }
	} // end for g...
     }
  } // end if m_use_supergrid  
  else //
  {
// tmp
//       printf("SG: supergrid not used!\n");

// Supergrid not used, but define arrays to simplify coding in some places.
     for( int g=0 ; g < mNumberOfGrids ; g++ )
     {
	for( i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	{
	   dcx(i,g)  = 0;
	   strx(i,g) = 1;
	   cornerx(i,g) = 1.;
	}
	for( j = m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	{
	   dcy(j,g)  = 0;
	   stry(j,g) = 1;
	   cornery(j,g) = 1.;
	}
	for( k = m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	{
	   dcz(k,g)  = 0.;
	   strz(k,g) = 1;
	   cornerz(k,g) = 1.;
	}
     }
  }
  
// tmp
// save the 1-D arrays cornerx, cornery, cornerz
// grid g=0 (for now)
  // g=0;
  // char fname[100];
  // sprintf(fname,"cx-%i.dat",m_myRank);
  // FILE *fp=fopen(fname,"w");
  // printf("Saving tmp file=%s, g=%i, m_iStart=%i, m_iEnd=%i\n", fname, g, m_iStart[g], m_iEnd[g]);
  // for ( i = m_iStart[g]; i<=m_iEnd[g]; i++ )
  //   fprintf(fp,"%i %e %e %e\n", i, cornerx(i,g), strx(i,g), strx(i,g)*dcx(i,g));
  // fclose(fp);  

  // sprintf(fname,"cy-%i.dat",m_myRank);
  // fp=fopen(fname,"w");
  // printf("Saving tmp file=%s, g=%i, m_jStart=%i, m_jEnd=%i\n", fname, g, m_jStart[g], m_jEnd[g]);
  // for ( j = m_jStart[g]; j<=m_jEnd[g]; j++ )
  //   fprintf(fp,"%i %e %e %e\n", j, cornery(j,g), stry(j,g), stry(j,g)*dcy(j,g));
  // fclose(fp);  

  // sprintf(fname,"cz-%i.dat",m_myRank);
  // fp=fopen(fname,"w");
  // printf("Saving tmp file=%s, g=%i, m_kStart=%i, m_kEnd=%i\n", fname, g, m_kStart[g], m_kEnd[g]);
  // for ( k = m_kStart[g]; k<=m_kEnd[g]; k++ )
  //   fprintf(fp,"%i %e %e %e\n", k, cornerz(k,g), strz(k,g), strz(k,g)*dcz(k,g));
  // fclose(fp);  
// end tmp

#undef dcx
#undef dcy
#undef dcz
#undef strx
#undef stry
#undef strz
#undef cornerx
#undef cornery
#undef cornerz
}

//-----------------------------------------------------------------------
void EW::material_ic( vector<Sarray>& a_mtrl )
{
// interface between curvilinear and top Cartesian grid
   if (topographyExists())
   {
      //      int nc=1;
      int g  = mNumberOfCartesianGrids-1;
      int gc = mNumberOfGrids-1;
      int nc = a_mtrl[g].ncomp();
      int q, i, j;
// inject values between lower boundary of gc and upper boundary of g
      for( j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	 for( i = m_iStart[g]; i <= m_iEnd[g]; i++ )
	 {
// assign ghost points in the Cartesian grid
	    for (q = 0; q < m_ghost_points; q++) // only once when m_ghost_points==1
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_mtrl[g](c,i,j,m_kStart[g] + q) = a_mtrl[gc](c,i,j,m_kEnd[gc]-2*m_ghost_points + q);
	    }
// assign ghost points in the Curvilinear grid
	    for (q = 0; q <= m_ghost_points; q++) // twice when m_ghost_points==1 (overwrites solution on the common grid line)
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_mtrl[gc](c,i,j,m_kEnd[gc]-q) = a_mtrl[g](c,i,j,m_kStart[g]+2*m_ghost_points - q);
	    }
	 }
   }
}

//-----------------------------------------------------------------------
void EW::perturb_velocities( vector<Sarray>& a_vs, vector<Sarray>& a_vp )
{
   int g = mNumberOfGrids-1;
   int p = m_random_dist/mGridSize[g]+1;
   int pz= m_random_distz/mGridSize[g]+1;
   //   int p = 2*m_random_dist/mGridSize[g]+1;
   //   int pz= 2*m_random_distz/mGridSize[g]+1;
   int ifirst = m_iStart[g];
   int ilast  = m_iEnd[g];
   int jfirst = m_jStart[g];
   int jlast  = m_jEnd[g];
   int kfirst = m_kStart[g];
   int klast  = m_kEnd[g];
   // Need to save random numbers on the overlap between the topograpy-grid 
   // and the uppermost Cartesian grid, in order to make the perturbation
   // continuous across the grid/grid interface. The overlap is always saved,
   // even if there is only a Cartesian grid, to simplify the code.
   //
   // NOTE: This might not work if grid refinement is put into the code later.
   //
   Sarray saverand(ifirst-p,ilast+p,jfirst-p,jlast+p,kfirst-pz,kfirst+pz);
   double* saverand_ptr = saverand.c_ptr();

   for ( g=0; g<mNumberOfGrids; g++)
   {
      double* vs_ptr  = a_vs[g].c_ptr();
      double* vp_ptr  = a_vp[g].c_ptr();
      ifirst = m_iStart[g];
      ilast  = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast  = m_jEnd[g];
      kfirst = m_kStart[g];
      klast  = m_kEnd[g];
      int nx = m_global_nx[g];
      int ny = m_global_ny[g];
      int nz = m_global_nz[g];
      int ghost = m_ghost_points;
      double h = mGridSize[g];
      Sarray pert(ifirst,ilast,jfirst,jlast,kfirst,klast);
      Sarray wgh(ifirst,ilast,jfirst,jlast,kfirst,klast);
      double* pert_ptr = pert.c_ptr();
      double* wgh_ptr = wgh.c_ptr();

      if( g == mNumberOfGrids-1 && topographyExists() )
      {
         F77_FUNC(randomfield3dc,RANDOMFIELD3DC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
		     &nx, &ny, &nz, &ghost, pert_ptr, wgh_ptr, &m_random_dist,
			      &m_random_distz, &h, mZ.c_ptr(), m_random_seed, saverand_ptr, &p, &pz );
         F77_FUNC(perturbvelocityc,PERTURBVELOCITYC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						      vs_ptr, vp_ptr, pert_ptr, &m_random_amp, &m_random_amp_grad,
						      mZ.c_ptr() );
      }
      else
      {
         F77_FUNC(randomfield3d,RANDOMFIELD3D)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
		     &nx, &ny, &nz, &ghost, pert_ptr, wgh_ptr, &m_random_dist,
						&m_random_distz, &h, m_random_seed, saverand_ptr, &p, &pz );
         F77_FUNC(perturbvelocity,PERTURBVELOCITY)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						    vs_ptr, vp_ptr, pert_ptr, &m_random_amp, &m_random_amp_grad,
						    &m_zmin[g], &h );
      }
   }
}
