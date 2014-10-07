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
//#include "impose_curvilinear_bc.h"
#include "F77_FUNC.h"

extern "C" {
void F77_FUNC(satt,SATT)(double *up, double *qs, double *dt, double *cfreq, int *ifirst, int *ilast, 
			 int *jfirst, int *jlast, int *kfirst, int *klast);


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
			       double*, double*, double*, int* );
void F77_FUNC(freesurfcurvi,FREESURFCURVI)(int*, int*, int*, int*, int*, int*, int*, int*,
					  double*, double*, double*, double*, double*, double* );
void F77_FUNC(freesurfcurvisg,FREESURFCURVISG)(int*, int*, int*, int*, int*, int*, int*, int*,
					       double*, double*, double*, double*, double*,
					       double*, double*, double* );
void F77_FUNC(bcfortsg, BCFORTSG)( int*, int*, int*, int*, int*, int*, 
			       int *, int*, int*, int*,
			       double*, double*, boundaryConditionType*, double *, double*, double*, double*,
			       double* bf0_p, double* bf1_p, 
			       double* bf2_p, double*bf3_p, 
			       double*bf4_p, double*bf5_p, 
			       double*, double*, double*, double*, double* );
   void F77_FUNC(bcfortanisg, BCFORTANISG)( int*, int*, int*, int*, int*, int*,  int*, int*, int*, int*,
					 double*, double*, boundaryConditionType*, double*, double*,
					 double*, double*, double*, double*, double*, double*,
					 double*, double* );

void F77_FUNC(twfrsurfz, TWFRSURFZ)(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p, int * kfirst_p, 
				  int * klast_p, double* h_p, int * k_p,
				  double* t_p, double* om_p, double* cv_p, double* ph_p,
				    double* bforce_side5_ptr, double* mu_ptr, double* la_ptr, double* zmin );
void F77_FUNC(twfrsurfzatt, TWFRSURFZATT)(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p, int * kfirst_p, 
				  int * klast_p, double* h_p, int * k_p,
				  double* t_p, double* om_p, double* cv_p, double* ph_p,
				    double* bforce_side5_ptr, double* mu_ptr, double* la_ptr, double* zmin );
void F77_FUNC(twfrsurfzsgstr, TWFRSURFZSGSTR)(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p, int * kfirst_p, 
				  int * klast_p, double* h_p, int * k_p,
				  double* t_p, double* om_p, double* cv_p, double* ph_p,
					      double* omstrx_p, double* omstry_p,
					      double* bforce_side5_ptr, double* mu_ptr, double* la_ptr, double* zmin );
void F77_FUNC(twfrsurfzsgstratt, TWFRSURFZSGSTRATT)(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p,
						    int * kfirst_p, int * klast_p, double* h_p, int * k_p,
				  double* t_p, double* om_p, double* cv_p, double* ph_p,
					      double* omstrx_p, double* omstry_p,
					      double* bforce_side5_ptr, double* mu_ptr, double* la_ptr, double* zmin );
void F77_FUNC(memvarforcesurf,MEMVARFORCESURF)( int*, int*, int*, int*, int*, double*, double*, double*, 
						   double*, double*, double*, double*, double*, double* );
   void F77_FUNC(memvarforcesurfc,MEMVARFORCESURFC)( int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, 
						  double*, double*, double*, double*, double*, double*, double* );
void F77_FUNC(twdirbdry,TWDIRBDRY)( int *wind_ptr, double *h_p, double *t_p, double *om_p, double * cv_p, 
				    double *ph_p,  double * bforce_side_ptr, double* zmin );
   void F77_FUNC(twdirbdryc,TWDIRBDRYC)( int* ifirst, int* ilast, int* jfirst, int* jlast, int* kfirst, int* klast,
                                         int *wind_ptr, double *t_p, double *om_p, double * cv_p, 
				      double *ph_p,  double * bforce_side_ptr, double* x, double* y, double* z );
void F77_FUNC(testsrc, TESTSRC )( double* f_ptr, int* ifirst, int* ilast, int* jfirst, int* jlast, int* kfirst,
				   int* klast, int* nz, int* wind, double* m_zmin, double* h, int* kx, int* ky, int* kz,
				   double* momgrid );
void F77_FUNC(addsgd,ADDSGD) (double* dt, double *h, double *a_U, double*a_Um, double*a_Up, 
			      double *sg_dc_x, double* sg_dc_y, double* sg_dc_z,
			      int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, double* damping_coefficient );
void F77_FUNC(addsgd4,ADDSGD4) (double* dt, double *h, double *a_Up, double*a_U, double*a_Um, double* Rho,
				double *sg_dc_x, double* sg_dc_y, double* sg_dc_z, double* sg_str_x, double* sg_str_y, double* sg_str_z,
				double* sg_corner_x, double* sg_corner_y, double* sg_corner_z,
				int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, double* damping_coefficient );
void F77_FUNC(addsgd6,ADDSGD6) (double* dt, double *h, double *a_Up, double*a_U, double*a_Um, double* Rho,
				double *sg_dc_x, double* sg_dc_y, double* sg_dc_z, double* sg_str_x, double* sg_str_y, double* sg_str_z,
				double* sg_corner_x, double* sg_corner_y, double* sg_corner_z,
				int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, double* damping_coefficient );
void F77_FUNC(addsgd4c,ADDSGD4C) (double* dt, double *a_Up, double*a_U, double*a_Um, double* Rho,
				  double *sg_dc_x, double* sg_dc_y, double* sg_str_x, double* sg_str_y, double* jac,
				  double* sg_corner_x, double* sg_corner_y, 
				  int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, double* damping_coefficient );
void F77_FUNC(addsgd6c,ADDSGD6C) (double* dt, double *a_Up, double*a_U, double*a_Um, double* Rho,
				  double *sg_dc_x, double* sg_dc_y, double* sg_str_x, double* sg_str_y, double* jac,
				  double* sg_corner_x, double* sg_corner_y, 
				  int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, double* damping_coefficient );
//  subroutine RAYDIRBDRY( bforce, wind, t, lambda, mu, rho, cr, 
// +     omega, alpha, h, zmin )
   void F77_FUNC(raydirbdry,RAYDIRBDRY)( double *bforce_side_ptr, int *wind_ptr, double *t, double *lambda,
					 double *mu, double *rho,
				         double *cr, double *omega, double *alpha, double *h, double *zmin );
   void F77_FUNC(twstensor,TWSTENSOR)( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* kz,
				       double* t, double* om, double* c, double* ph, double* xx, double* yy, double* zz,
				       double* tau, double* mu, double* lambda );
   void F77_FUNC(twstensoratt,TWSTENSORATT)( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* kz,
				       double* t, double* om, double* c, double* ph, double* xx, double* yy, double* zz,
				       double* tau, double* mu, double* lambda );
   void F77_FUNC(twstensorsg,TWSTENSORSG)( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* kz,
				       double* t, double* om, double* c, double* ph, double* xx, double* yy, double* zz,
					   double* tau, double* mu, double* lambda, double* omstrx, double* omstry );
   void F77_FUNC(twstensorsgatt,TWSTENSORSGATT)( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, 
                                          int* kz,
				       double* t, double* om, double* c, double* ph, double* xx, double* yy, double* zz,
					   double* tau, double* mu, double* lambda, double* omstrx, double* omstry );
   void F77_FUNC(getsurfforcing,GETSURFFORCING)( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast,
						 int* k, double* met, double* jac, double* tau, double* forcing );
   void F77_FUNC(subsurfforcing,SUBSURFFORCING)( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast,
						 int* k, double* met, double* jac, double* tau, double* forcing );
   void F77_FUNC(getsurfforcinggh,GETSURFFORCINGGH)( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast,
						    int* k, double* h, double* tau, double* forcing, double* amp, double* xc,
						    double* yc, double* xl, double* yl );
   void F77_FUNC(getsurfforcingsg,GETSURFFORCINGSG)( 
	   int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* k,
           double* met, double* jac, double* tau, double* strx, double* stry, double* forcing );
   void F77_FUNC(subsurfforcingsg,SUBSURFFORCINGSG)( 
	   int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* k,
           double* met, double* jac, double* tau, double* strx, double* stry, double* forcing );
   void F77_FUNC(addbstressc,ADDBSTRESSC)( int*, int*, int*, int*, int*, int*, int*, double*,  
					   double*,  double*,  double*, double*, int*,  double*, char*,
					   int*, int*,  double*,  double* );
   void F77_FUNC(addbstresswresc,ADDBSTRESSWRESC)( int*, int*, int*, int*, int*, int*, int*, double*, double*,    
						   double*, double*, double*, double*, double*, double*, int*, 
						   double*, double*, double*, double*, double*, double*, double*,
						   int*, double*, double* );
   void F77_FUNC(solveattfreec,SOLVEATTFREEC)( int*, int*, int*, int*, int*, int*, double*, double*,    
					       double*, double*, double*, double*, double*, double*, int*, double*, double* );
   void F77_FUNC(solveattfreeac,SOLVEATTFREEAC)( int*, int*, int*, int*, int*, int*, double*, double*, double*);
}


//--------------------------------------------------------------------
void EW::solve( vector<Source*> & a_Sources, vector<TimeSeries*> & a_TimeSeries )
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
    
// tmp
    // if (proc_zero() && point_sources.size()>0)
    // {
    //   printf("Saving one un-filtered original time function\n");
	 
    //   FILE *tf=fopen("g0.dat","w");
    //   double t;
    //   double gt;
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
//      double t;
//      double gt, gt1, gt2;
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
  bool output_timefunc = true;
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
  initialData(mTstart, U, AlphaVE);
  initialData(mTstart-mDt, Um, AlphaVEm );
  
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
  
// save any images for cycle = 0 (initial data) ?
  update_images( 0, t, U, Um, Up, mRho, mMu, mLambda, a_Sources, 1 );
  for( int i3 = 0 ; i3 < mImage3DFiles.size() ; i3++ )
    mImage3DFiles[i3]->update_image( t, 0, mDt, U, mRho, mMu, mLambda, mRho, mMu, mLambda, mQp, mQs, mPath, mZ );

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
    double * lowZ = new double[3*mNumberOfGrids];
    double * interiorZ = new double[3*mNumberOfGrids];
    double * highZ = new double[3*mNumberOfGrids];
    //    double lowZ[3], interiorZ[3], highZ[3];
    bndryInteriorDifference( F, Up, lowZ, interiorZ, highZ );

    double* tmp= new double[3*mNumberOfGrids];
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = lowZ[i];
    MPI_Reduce( tmp, lowZ, 3*mNumberOfGrids, MPI_DOUBLE, MPI_MAX, 0, m_cartesian_communicator );
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = interiorZ[i];
    MPI_Reduce( tmp, interiorZ, 3*mNumberOfGrids, MPI_DOUBLE, MPI_MAX, 0, m_cartesian_communicator );
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = highZ[i];
    MPI_Reduce( tmp, highZ, 3*mNumberOfGrids, MPI_DOUBLE, MPI_MAX, 0, m_cartesian_communicator );

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
  
// c test accuracy of forcing
    evalRHS( U, mMu, mLambda, Lu, AlphaVE ); // save Lu in composite grid 'Lu'
    Force( t, F, point_sources );
    exactAccTwilight( t, Uacc ); // save Utt in Uacc
    test_RhoUtt_Lu( Uacc, Lu, F, lowZ, interiorZ, highZ );

    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = lowZ[i];
    MPI_Reduce( tmp, lowZ, 3*mNumberOfGrids, MPI_DOUBLE, MPI_MAX, 0, m_cartesian_communicator );
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = interiorZ[i];
    MPI_Reduce( tmp, interiorZ, 3*mNumberOfGrids, MPI_DOUBLE, MPI_MAX, 0, m_cartesian_communicator );
    for( int i=0 ; i < 3*mNumberOfGrids ; i++ )
      tmp[i] = highZ[i];
    MPI_Reduce( tmp, highZ, 3*mNumberOfGrids, MPI_DOUBLE, MPI_MAX, 0, m_cartesian_communicator );

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

// enforce bc on initial data
// U
// communicate across processor boundaries
  for(int g=0 ; g < mNumberOfGrids ; g++ )
    communicate_array( U[g], g );
// boundary forcing
  cartesian_bc_forcing( t, BCForcing, a_Sources );
  if( m_use_attenuation && m_number_mechanisms > 0 )
     addAttToFreeBcForcing( AlphaVE, BCForcing, m_sbop );
// enforce boundary condition
  if( m_anisotropic )
     enforceBCanisotropic( U, mC, t, BCForcing );
  else
     enforceBC( U, mMu, mLambda, t, BCForcing );   

// Um
// communicate across processor boundaries
  for(int g=0 ; g < mNumberOfGrids ; g++ )
    communicate_array( Um[g], g );
// boundary forcing
  cartesian_bc_forcing( t-mDt, BCForcing, a_Sources );
  if( m_use_attenuation && m_number_mechanisms > 0 )
     addAttToFreeBcForcing( AlphaVEm, BCForcing, m_sbop );

// enforce boundary condition
  if( m_anisotropic )
     enforceBCanisotropic( Um, mC, t-mDt, BCForcing );
  else
     enforceBC( Um, mMu, mLambda, t-mDt, BCForcing );

  if (m_twilight_forcing && getVerbosity()>=3)
  {
// more testing
    if ( proc_zero() )
    {
      printf("Checking the accuracy of the initial data\n");
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
    
// Begin time stepping loop
  for( int currentTimeStep = beginCycle; currentTimeStep <= mNumberOfTimeSteps; currentTimeStep++)
  {    
    time_measure[0] = MPI_Wtime();

// all types of forcing...
    Force( t, F, point_sources );

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

    if( m_checkfornan )
       check_for_nan( Lu, 1, "Lu pred. " );

// take predictor step, store in Up
    evalPredictor( Up, U, Um, mRho, Lu, F );    

    time_measure[1] = MPI_Wtime();
    time_measure[2] = MPI_Wtime();

// communicate across processor boundaries
    for(int g=0 ; g < mNumberOfGrids ; g++ )
       communicate_array( Up[g], g );

// calculate boundary forcing at time t+mDt
    cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );

// update ghost points in Up
    if( m_anisotropic )
       enforceBCanisotropic( Up, mC, t+mDt, BCForcing );
    else
       enforceBC( Up, mMu, mLambda, t+mDt, BCForcing );

    if( m_checkfornan )
       check_for_nan( Up, 1, "U pred. " );

    if( m_use_attenuation && m_number_mechanisms > 0 )
    {
// Update memory variables
       updateMemoryVariables( AlphaVEp, AlphaVEm, Up, U, Um, t );
// Impose coupled free surface boundary condition
       enforceBCfreeAtt( Up, U, Um, mMu, mLambda, AlphaVEp, AlphaVEm, BCForcing, m_sbop, t );
    }
    time_measure[3] = time_measure[4] = MPI_Wtime();

// get 4th order in time
    if (mOrder == 4)
    {
       Force_tt( t, F, point_sources );
       evalDpDmInTime( Up, U, Um, Uacc ); // store result in Uacc

       if( m_checkfornan )
	  check_for_nan( Uacc, 1, "uacc " );

       if( m_use_attenuation && m_number_mechanisms > 0 )
          evalDpDmInTimeAtt( AlphaVEp, AlphaVE, AlphaVEm ); // store AlphaVEacc in AlphaVEm

       if( m_anisotropic )
	  evalRHSanisotropic( Uacc, mC, Lu );
       else
	  evalRHS( Uacc, mMu, mLambda, Lu, AlphaVEm );
       
       if( m_checkfornan )
	  check_for_nan( Lu, 1, "L(uacc) " );

       evalCorrector( Up, mRho, Lu, F );
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
       time_measure[4] = MPI_Wtime();

// also check out EW::update_all_boundaries 
// communicate across processor boundaries
       for(int g=0 ; g < mNumberOfGrids ; g++ )
	  communicate_array( Up[g], g );

// calculate boundary forcing at time t+mDt (do we really need to call this fcn again???)
       cartesian_bc_forcing( t+mDt, BCForcing, a_Sources );
       if( m_use_attenuation && (m_number_mechanisms > 0) )
       {
	  addAttToFreeBcForcing( AlphaVEp, BCForcing, m_sbop );
       }
// update ghost points in Up
       if( m_anisotropic )
	  enforceBCanisotropic( Up, mC, t+mDt, BCForcing );
       else
	  enforceBC( Up, mMu, mLambda, t+mDt, BCForcing );
    }
    if( m_checkfornan )
       check_for_nan( Up, 1, "Up" );

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
     double errInf, errL2, solInf, solL2;

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
	    delete[] BCForcing[g][side];
      delete[] BCForcing[g];
   }
   for( int s = 0 ; s < point_sources.size(); s++ )
      delete point_sources[s];

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
void EW::enforceBC( vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		    double t, vector<double **> & a_BCForcing )
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
// for periodic bc, a_BCForcing[g][s] == NULL, so you better not access the
// theses arrays in that case
    bforce_side0_ptr = a_BCForcing[g][0]; // low-i bndry forcing array pointer
    bforce_side1_ptr = a_BCForcing[g][1]; // high-i bndry forcing array pointer
    bforce_side2_ptr = a_BCForcing[g][2]; // low-j bndry forcing array pointer
    bforce_side3_ptr = a_BCForcing[g][3]; // high-j bndry forcing array pointer
    bforce_side4_ptr = a_BCForcing[g][4]; // low-k bndry forcing array pointer
    bforce_side5_ptr = a_BCForcing[g][5]; // high-k bndry forcing array pointer
    
    if( usingSupergrid() )
    {
       F77_FUNC(bcfortsg, BCFORTSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
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
          F77_FUNC(freesurfcurvisg,FREESURFCURVISG)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						    &nz, &side, u_ptr, mu_ptr, la_ptr, mMetric.c_ptr(),
						    m_sbop, bforce_side4_ptr, m_sg_str_x[g],
						    m_sg_str_y[g] );
       }
    }
    else
    {
       F77_FUNC(bcfort, BCFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
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
          F77_FUNC(freesurfcurvi,FREESURFCURVI)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
			       &side, u_ptr, mu_ptr, la_ptr, mMetric.c_ptr(), m_sbop, bforce_side4_ptr );
       }
       if( topo == 1 && m_bcType[g][5] == bStressFree )
       {
	  side = 6;
          F77_FUNC(freesurfcurvi,FREESURFCURVI)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
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

      //      double nrm[3]={0,0,0};
      int q, i, j;
// inject solution values between lower boundary of gc and upper boundary of g
      for( j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	 for( i = m_iStart[g]; i <= m_iEnd[g]; i++ )
	 {
// assign ghost points in the Cartesian grid
	    for (q = 0; q < m_ghost_points; q++) // only once when m_ghost_points==1
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_U[g](c,i,j,m_kStart[g] + q) = a_U[gc](c,i,j,m_kEnd[gc]-2*m_ghost_points + q);
	    }
// assign ghost points in the Curvilinear grid
	    for (q = 0; q <= m_ghost_points; q++) // twice when m_ghost_points==1 (overwrites solution on the common grid line)
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_U[gc](c,i,j,m_kEnd[gc]-q) = a_U[g](c,i,j,m_kStart[g]+2*m_ghost_points - q);
	    }
	    //	    // Verify that the grid lines are the same
	    //            if( i>= 1 && i <= m_global_nx[g] && j >= 1 && j <= m_global_ny[g] )
	    //            for( int c=1; c <= nc ; c++ )
	    //	    {
	    //	       double ndiff=fabs(a_U[gc](c,i,j,m_kEnd[gc]-m_ghost_points)-a_U[g](c,i,j,m_kStart[g]+m_ghost_points));
	    //	       if( ndiff > nrm[c-1] )
	    //		  nrm[c-1] = ndiff;
	    //	    }
	 }
      //      cout << "Difference of curvilinear and cartesian common grid line = "<< nrm[0] << " " << nrm[1] << " " << nrm[2] << endl;
   }
}

//---------------------------------------------------------------------------
void EW::enforceBCanisotropic( vector<Sarray> & a_U, vector<Sarray>& a_C,
		    double t, vector<double **> & a_BCForcing )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz;
  double *u_ptr, *c_ptr, h;
  boundaryConditionType *bcType_ptr;
  double *bforce_side0_ptr, *bforce_side1_ptr, *bforce_side2_ptr, *bforce_side3_ptr, *bforce_side4_ptr, *bforce_side5_ptr;
  int *wind_ptr;
  double om=0, ph=0, cv=0;
    
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
// for periodic bc, a_BCForcing[g][s] == NULL, so you better not access the
// theses arrays in that case
    bforce_side0_ptr = a_BCForcing[g][0]; // low-i bndry forcing array pointer
    bforce_side1_ptr = a_BCForcing[g][1]; // high-i bndry forcing array pointer
    bforce_side2_ptr = a_BCForcing[g][2]; // low-j bndry forcing array pointer
    bforce_side3_ptr = a_BCForcing[g][3]; // high-j bndry forcing array pointer
    bforce_side4_ptr = a_BCForcing[g][4]; // low-k bndry forcing array pointer
    bforce_side5_ptr = a_BCForcing[g][5]; // high-k bndry forcing array pointer
    
    //    if( usingSupergrid() )
    //    {
       F77_FUNC(bcfortanisg, BCFORTANISG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
				       wind_ptr, &nx, &ny, &nz,
				       u_ptr, &h, bcType_ptr, m_sbop, c_ptr, 
				       bforce_side0_ptr, bforce_side1_ptr, 
				       bforce_side2_ptr, bforce_side3_ptr, 
				       bforce_side4_ptr, bforce_side5_ptr, 
				       m_sg_str_x[g], m_sg_str_y[g] );
       //    }
    //    else
    //    {
       //       F77_FUNC(bcfort, BCFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
       //				 wind_ptr, &nx, &ny, &nz,
       //				 u_ptr, &h, bcType_ptr, m_sbop, c_ptr,
       //				 bforce_side0_ptr, bforce_side1_ptr, 
       //				 bforce_side2_ptr, bforce_side3_ptr, 
       //				 bforce_side4_ptr, bforce_side5_ptr )
       //				
    //  }
  }
}

//-----------------------------------------------------------------------
void EW::update_curvilinear_cartesian_interface( vector<Sarray>& a_U )
{
   if (topographyExists())
   {
      int g  = mNumberOfCartesianGrids-1;
      int gc = mNumberOfGrids-1;
      int q, i, j;
      int nc = a_U[0].m_nc;
// inject solution values between lower boundary of gc and upper boundary of g
      for( j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	 for( i = m_iStart[g]; i <= m_iEnd[g]; i++ )
	 {
// assign ghost points in the Cartesian grid
	    for (q = 0; q < m_ghost_points; q++) 
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_U[g](c,i,j,m_kStart[g] + q) = a_U[gc](c,i,j,m_kEnd[gc]-2*m_ghost_points + q);
	    }
// assign ghost points in the Curvilinear grid
	    for (q = 0; q <= m_ghost_points; q++) // (overwrites solution on the common grid line)
	    {
	       for( int c = 1; c <= nc ; c++ )
		  a_U[gc](c,i,j,m_kEnd[gc]-q) = a_U[g](c,i,j,m_kStart[g]+2*m_ghost_points - q);
	    }
	 }
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
void EW::cartesian_bc_forcing(double t, vector<double **> & a_BCForcing,
			      vector<Source*>& a_sources )
// assign the boundary forcing arrays a_BCForcing[g][side]
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz;
  double *u_ptr, *mu_ptr, *la_ptr, h, zmin;
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

    if (m_twilight_forcing)
    {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;

// the following code can probably be improved by introducing a loop over all sides,
// but bStressFree is only implemented for side=4 and 5, so there must be some special cases
      int k = 1;
      if (m_bcType[g][0] == bDirichlet || m_bcType[g][0] == bSuperGrid )
      {
         if( !curvilinear )
	    F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[0], &h, &t, &om, &cv, &ph, bforce_side0_ptr, &m_zmin[g] );
	 else
	    F77_FUNC(twdirbdryc,TWDIRBDRYC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                                             &wind_ptr[0], &t, &om, &cv, &ph, bforce_side0_ptr,
					     mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
      }

      if (m_bcType[g][1] == bDirichlet || m_bcType[g][1] == bSuperGrid )
      {
         if( !curvilinear )
	    F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6], &h, &t, &om, &cv, &ph, bforce_side1_ptr, &m_zmin[g] );
	 else
	    F77_FUNC(twdirbdryc,TWDIRBDRYC)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					    &wind_ptr[6], &t, &om, &cv, &ph, bforce_side1_ptr,
					    mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
      }

      if (m_bcType[g][2] == bDirichlet || m_bcType[g][2] == bSuperGrid)
      {
	 if( !curvilinear )
	    F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6*2], &h, &t, &om, &cv, &ph, bforce_side2_ptr, &m_zmin[g] );
         else
	    F77_FUNC(twdirbdryc,TWDIRBDRYC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					     &wind_ptr[6*2], &t, &om, &cv, &ph, bforce_side2_ptr,
					     mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
      }

      if (m_bcType[g][3] == bDirichlet || m_bcType[g][3] == bSuperGrid)
      {
         if( !curvilinear )
	    F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6*3], &h, &t, &om, &cv, &ph, bforce_side3_ptr, &m_zmin[g] );
         else
	    F77_FUNC(twdirbdryc,TWDIRBDRYC)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					    &wind_ptr[6*3], &t, &om, &cv, &ph, bforce_side3_ptr,
					    mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
      }

      if (m_bcType[g][4] == bDirichlet || m_bcType[g][4] == bSuperGrid)
      {
	 if( !curvilinear )
	    F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6*4], &h, &t, &om, &cv, &ph, bforce_side4_ptr, &m_zmin[g] );
         else
	    F77_FUNC(twdirbdryc,TWDIRBDRYC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					     &wind_ptr[6*4], &t, &om, &cv, &ph, bforce_side4_ptr,
					     mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
      }
      else if (m_bcType[g][4] == bStressFree)
      {
	 k = 1;
         if( usingSupergrid() && !curvilinear )
	 {
            double omstrx = m_supergrid_taper_x.get_tw_omega();
            double omstry = m_supergrid_taper_y.get_tw_omega();
	    F77_FUNC(twfrsurfzsgstr, TWFRSURFZSGSTR)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
						      &klast, &h, &k, &t, &om, &cv, &ph, &omstrx, &omstry,
						      bforce_side4_ptr, mu_ptr, la_ptr, &m_zmin[g] );
            if( m_use_attenuation )
	    {
	       double* mua_ptr    = mMuVE[g][0].c_ptr();
	       double* laa_ptr    = mLambdaVE[g][0].c_ptr();
	       F77_FUNC(twfrsurfzsgstratt, TWFRSURFZSGSTRATT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
							       &klast, &h, &k, &t, &om, &cv, &ph, &omstrx, &omstry,
							       bforce_side4_ptr, mua_ptr, laa_ptr, &m_zmin[g] );
	       
	    }
	 }
         else if( !usingSupergrid() && curvilinear )
	 {
	    // Stress tensor on boundary
            Sarray tau(6,ifirst,ilast,jfirst,jlast,1,1);
	    // Get twilight stress tensor, tau.
            F77_FUNC(twstensor,TWSTENSOR)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					   &k, &t, &om, &cv, &ph,
					   mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(), mu_ptr, la_ptr );
	    // Compute boundary forcing for given stress tensor, tau.

	    F77_FUNC(getsurfforcing,GETSURFFORCING)( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
						     &klast, &k, mMetric.c_ptr(), mJ.c_ptr(),
						     tau.c_ptr(), bforce_side4_ptr );
	    //	    F77_FUNC(getsurfforcinggh,GETSURFFORCINGGH)( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
	    //							 &klast, &k, &h, tau.c_ptr(), bforce_side4_ptr,
	    //							 &m_GaussianAmp, &m_GaussianXc, &m_GaussianYc,
	    //							 &m_GaussianLx, &m_GaussianLy );
            if( m_use_attenuation )
	    {
	       double* mua_ptr    = mMuVE[g][0].c_ptr();
	       double* laa_ptr    = mLambdaVE[g][0].c_ptr();
	       F77_FUNC(twstensoratt,TWSTENSORATT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						    &k, &t, &om, &cv, &ph,
						    mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(), mua_ptr, laa_ptr );
	       F77_FUNC(subsurfforcing,SUBSURFFORCING)( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
							&klast, &k, mMetric.c_ptr(), mJ.c_ptr(),
							tau.c_ptr(), bforce_side4_ptr );
	    }
	 }
	 else if( !usingSupergrid() && !curvilinear )
	 {
	    F77_FUNC(twfrsurfz, TWFRSURFZ)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					    &klast, &h, &k, &t, &om, &cv, &ph,
					    bforce_side4_ptr, mu_ptr, la_ptr, &m_zmin[g] );
            if( m_use_attenuation )
	    {
	       double* mua_ptr    = mMuVE[g][0].c_ptr();
	       double* laa_ptr    = mLambdaVE[g][0].c_ptr();
	       F77_FUNC(twfrsurfzatt, TWFRSURFZATT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					    &klast, &h, &k, &t, &om, &cv, &ph,
					    bforce_side4_ptr, mua_ptr, laa_ptr, &m_zmin[g] );
	       
	    }
	 }
	 else if( usingSupergrid() && curvilinear )
	 {
            double omstrx = m_supergrid_taper_x.get_tw_omega();
            double omstry = m_supergrid_taper_y.get_tw_omega();

	    // Stress tensor on boundary
            Sarray tau(6,ifirst,ilast,jfirst,jlast,1,1);
	    // Get twilight stress tensor, tau.
            F77_FUNC(twstensorsg,TWSTENSORSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					   &k, &t, &om, &cv, &ph,
					   mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(),
					       mu_ptr, la_ptr, &omstrx, &omstry );
	    // Compute boundary forcing for given stress tensor, tau.
	    F77_FUNC(getsurfforcingsg,GETSURFFORCINGSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
			         &klast, &k, mMetric.c_ptr(), mJ.c_ptr(),
				 tau.c_ptr(), m_sg_str_x[g], m_sg_str_y[g], bforce_side4_ptr );

            if( m_use_attenuation )
	    {
	       double* mua_ptr    = mMuVE[g][0].c_ptr();
	       double* laa_ptr    = mLambdaVE[g][0].c_ptr();
	       F77_FUNC(twstensorsgatt,TWSTENSORSGATT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
							&k, &t, &om, &cv, &ph,
							mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), tau.c_ptr(),
							mua_ptr, laa_ptr, &omstrx, &omstry );
	       F77_FUNC(subsurfforcingsg,SUBSURFFORCINGSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst,
							    &klast, &k, mMetric.c_ptr(), mJ.c_ptr(),
							    tau.c_ptr(), m_sg_str_x[g], m_sg_str_y[g],
							    bforce_side4_ptr );
	    }
      	 }
      }

      if (m_bcType[g][5] == bDirichlet || m_bcType[g][5] == bSuperGrid)
      {
	 if( !curvilinear )
	    F77_FUNC(twdirbdry,TWDIRBDRY)( &wind_ptr[6*5], &h, &t, &om, &cv, &ph, bforce_side5_ptr, &m_zmin[g] );
	 else
	    F77_FUNC(twdirbdryc,TWDIRBDRYC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					     &wind_ptr[6*5], &t, &om, &cv, &ph, bforce_side5_ptr,
					     mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
      }
      else if (m_bcType[g][5] == bStressFree)
      {
	 k = nz;
         if( usingSupergrid() )
	 {
            double omstrx = m_supergrid_taper_x.get_tw_omega();
            double omstry = m_supergrid_taper_y.get_tw_omega();
	    F77_FUNC(twfrsurfzsgstr, TWFRSURFZSGSTR)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
						      &klast, &h, &k, &t, &om, &cv, &ph, &omstrx, &omstry,
						      bforce_side5_ptr, mu_ptr, la_ptr, &m_zmin[g] );
            if( m_use_attenuation )
	    {
	       double* mua_ptr    = mMuVE[g][0].c_ptr();
	       double* laa_ptr    = mLambdaVE[g][0].c_ptr();
	       F77_FUNC(twfrsurfzsgstratt, TWFRSURFZSGSTRATT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
							       &klast, &h, &k, &t, &om, &cv, &ph, &omstrx, &omstry,
							       bforce_side5_ptr, mua_ptr, laa_ptr, &m_zmin[g] );
	       
	    }
	 }
	 else
	 {
	    F77_FUNC(twfrsurfz, TWFRSURFZ)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					    &klast, &h, &k, &t, &om, &cv, &ph,
					    bforce_side5_ptr, mu_ptr, la_ptr, &m_zmin[g] );
            if( m_use_attenuation )
	    {
	       double* mua_ptr    = mMuVE[g][0].c_ptr();
	       double* laa_ptr    = mLambdaVE[g][0].c_ptr();
	       F77_FUNC(twfrsurfzatt, TWFRSURFZATT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
						     &klast, &h, &k, &t, &om, &cv, &ph,
						     bforce_side5_ptr, mua_ptr, laa_ptr, &m_zmin[g] );
	    }
	 }
      }
    }
    else if (m_rayleigh_wave_test)
    {
      int q;
      double lambda, mu, rho, cr, omega, alpha;
      
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
	F77_FUNC(raydirbdry,RAYDIRBDRY)( bforce_side5_ptr, &wind_ptr[6*5], &t, &lambda, &mu, &rho, &cr, 
					 &omega, &alpha, &h, &zmin );
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
      int q;
      for (q=0; q<3*m_NumberOfBCPoints[g][0]; q++)
	bforce_side0_ptr[q] = 0.;
      for (q=0; q<3*m_NumberOfBCPoints[g][1]; q++)
	bforce_side1_ptr[q] = 0.;
      for (q=0; q<3*m_NumberOfBCPoints[g][2]; q++)
	bforce_side2_ptr[q] = 0.;
      for (q=0; q<3*m_NumberOfBCPoints[g][3]; q++)
	bforce_side3_ptr[q] = 0.;
      for (q=0; q<3*m_NumberOfBCPoints[g][4]; q++)
	bforce_side4_ptr[q] = 0.;
      for (q=0; q<3*m_NumberOfBCPoints[g][5]; q++)
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

//-----------------------------------------------------------------------
void EW::test_sources( vector<GridPointSource*>& a_point_sources,
		       vector<Source*>& a_global_unique_sources, vector<Sarray>& a_F )
{
// Check the source discretization
  int kx[3] = {0,0,0};
  int ky[3] = {0,0,0};
  int kz[3] = {0,0,0};
  double moments[3], momexact[3];
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
     testSourceDiscretization( kx, ky, kz, moments, a_point_sources, a_F );
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
  testSourceDiscretization( kx, ky, kz, moments, a_point_sources, a_F );
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
  testSourceDiscretization( kx, ky, kz, moments, a_point_sources, a_F );
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
  testSourceDiscretization( kx, ky, kz, moments, a_point_sources, a_F );
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
				   double moments[3],
				   vector<GridPointSource*>& point_sources,
				   vector<Sarray>& F )
{
   // Evaluate sources at a large time (assume that the time function is=1 at t=infinity)
   // Compute moments, integrals of the source times polynomials of degree (kx,ky,kz).
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast;
  double h;

//tmp  
  if (proc_zero())
    printf("Inside testSourceDiscretization\n");
  
  // Impose source
  for( int g=0 ; g<mNumberOfGrids; g++ )
     F[g].set_to_zero();
  for( int s= 0 ; s < point_sources.size() ; s++ ) 
  {
     double fxyz[3];
     point_sources[s]->getFxyz_notime( fxyz );
     int g = point_sources[s]->m_grid;
     F[g](1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
     F[g](2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
     F[g](3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
  }
  double momgrid[3]={0,0,0};
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];  
    h      = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    double* f_ptr = F[g].c_ptr();
    int wind[6];
    wind[0] = m_iStartInt[g];
    wind[1] = m_iEndInt[g];
    wind[2] = m_jStartInt[g];
    wind[3] = m_jEndInt[g];
    wind[4] = m_kStartInt[g];
    wind[5] = m_kEndInt[g];
    int nz = m_global_nz[g];
    F77_FUNC( testsrc, TESTSRC )( f_ptr, &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
				  &nz, wind, &m_zmin[g], &h, kx, ky, kz, momgrid );
  }
  MPI_Allreduce( momgrid, moments, 3, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );
}

//-----------------------------------------------------------------------
void EW::extractRecordData(TimeSeries::receiverMode mode, int i0, int j0, int k0, int g0, 
			   vector<double> &uRec, vector<Sarray> &Um2, vector<Sarray> &U)
{
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
      double factor = 1.0/(2*mGridSize[g0]);
      uRec[0] = ((U[g0](1,i0+1, j0, k0) - U[g0](1,i0-1, j0, k0)+
		  U[g0](2,i0, j0+1, k0) - U[g0](2,i0, j0-1, k0)+
		  U[g0](3,i0, j0, k0+1) - U[g0](3,i0, j0, k0-1))*factor);
    }
    else // must be curvilinear
    {
//      int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
       double factor = 0.5/sqrt(mJ(i0,j0,k0));
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
      double factor = 1.0/(2*mGridSize[g0]);
      double duydx = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0))*factor;
      double duzdx = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0))*factor;
      double duxdy = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0))*factor;
      double duzdy = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0))*factor;
      double duxdz = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))*factor;
      double duydz = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))*factor;
//       if( m_xycomponent )
//       {
      uRec[0] = ( duzdy-duydz );
      uRec[1] = ( duxdz-duzdx );
      uRec[2] = ( duydx-duxdy );
//       }
//       else
//       {
// 	 double uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
// 	 double uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
// 	 mRecordedUX.push_back( uew );
// 	 mRecordedUY.push_back( uns );
// 	 mRecordedUZ.push_back( -(duydx-duxdy) );
//       }
    }
    else // must be curvilinear
    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      double factor = 0.5/sqrt(mJ(i0,j0,k0));
      double duxdq = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0));
      double duydq = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0));
      double duzdq = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0));
      double duxdr = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0));
      double duydr = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0));
      double duzdr = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0));
      double duxds = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1));
      double duyds = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1));
      double duzds = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1));
      double duzdy = mMetric(1,i0,j0,k0)*duzdr+mMetric(3,i0,j0,k0)*duzds;
      double duydz = mMetric(4,i0,j0,k0)*duyds;
      double duxdz = mMetric(4,i0,j0,k0)*duxds;
      double duzdx = mMetric(1,i0,j0,k0)*duzdq+mMetric(2,i0,j0,k0)*duzds;
      double duydx = mMetric(1,i0,j0,k0)*duydq+mMetric(2,i0,j0,k0)*duyds;
      double duxdy = mMetric(1,i0,j0,k0)*duxdr+mMetric(3,i0,j0,k0)*duxds;
//       if( m_xycomponent )
//       {
      uRec[0] = (duzdy-duydz)*factor;
      uRec[1] = (duxdz-duzdx)*factor;
      uRec[2] = (duydx-duxdy)*factor;
//       }
//       else
//       {
// 	 double uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
// 	 double uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
// 	 mRecordedUX.push_back( uew*factor );
// 	 mRecordedUY.push_back( uns*factor );
// 	 mRecordedUZ.push_back( -(duydx-duxdy)*factor );
//       }
    }
  } // end Curl
  else if(mode == TimeSeries::Strains)
  {
    uRec.resize(6);
    if (g0 < mNumberOfCartesianGrids) // must be a Cartesian grid
    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
      double factor = 1.0/(2*mGridSize[g0]);
      double duydx = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0))*factor;
      double duzdx = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0))*factor;
      double duxdy = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0))*factor;
      double duzdy = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0))*factor;
      double duxdz = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1))*factor;
      double duydz = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1))*factor;
      double duxdx = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))*factor;
      double duydy = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))*factor;
      double duzdz = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1))*factor;
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
      double factor = 0.5/sqrt(mJ(i0,j0,k0));
      double duxdq = (U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0));
      double duydq = (U[g0](2,i0+1,j0,k0) - U[g0](2,i0-1,j0,k0));
      double duzdq = (U[g0](3,i0+1,j0,k0) - U[g0](3,i0-1,j0,k0));
      double duxdr = (U[g0](1,i0,j0+1,k0) - U[g0](1,i0,j0-1,k0));
      double duydr = (U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0));
      double duzdr = (U[g0](3,i0,j0+1,k0) - U[g0](3,i0,j0-1,k0));
      double duxds = (U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1));
      double duyds = (U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1));
      double duzds = (U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1));
      double duzdy = (mMetric(1,i0,j0,k0)*duzdr+mMetric(3,i0,j0,k0)*duzds)*factor;
      double duydz = (mMetric(4,i0,j0,k0)*duyds)*factor;
      double duxdz = (mMetric(4,i0,j0,k0)*duxds)*factor;
      double duzdx = (mMetric(1,i0,j0,k0)*duzdq+mMetric(2,i0,j0,k0)*duzds)*factor;
      double duydx = (mMetric(1,i0,j0,k0)*duydq+mMetric(2,i0,j0,k0)*duyds)*factor;
      double duxdy = (mMetric(1,i0,j0,k0)*duxdr+mMetric(3,i0,j0,k0)*duxds)*factor;
      double duxdx = ( mMetric(1,i0,j0,k0)*(U[g0](1,i0+1,j0,k0) - U[g0](1,i0-1,j0,k0))+
		       mMetric(2,i0,j0,k0)*(U[g0](1,i0,j0,k0+1) - U[g0](1,i0,j0,k0-1)) )*factor;
      double duydy = ( mMetric(1,i0,j0,k0)*(U[g0](2,i0,j0+1,k0) - U[g0](2,i0,j0-1,k0))+
		       mMetric(3,i0,j0,k0)*(U[g0](2,i0,j0,k0+1) - U[g0](2,i0,j0,k0-1)) )*factor;
      double duzdz = ( mMetric(4,i0,j0,k0)*(U[g0](3,i0,j0,k0+1) - U[g0](3,i0,j0,k0-1)) )*factor;
      uRec[0] = ( duxdx );
      uRec[1] = ( duydy );
      uRec[2] = ( duzdz );
      uRec[3] = ( 0.5*(duydx+duxdy) );
      uRec[4] = ( 0.5*(duzdx+duxdz) );
      uRec[5] = ( 0.5*(duydz+duzdy) );
   }
  } // end Strains
  
  return;
}

//---------------------------------------------------------------------------
void EW::addSuperGridDamping(vector<Sarray> & a_Up, vector<Sarray> & a_U,
			     vector<Sarray> & a_Um, vector<Sarray> & a_Rho )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *up_ptr, *u_ptr, *um_ptr, dt2i;
  
  int g;
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    u_ptr   = a_U[g].c_ptr();
    um_ptr  = a_Um[g].c_ptr();
    double* rho_ptr = a_Rho[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    //    F77_FUNC(addsgd,ADDSGD) ( &mDt, &mGridSize[g], up_ptr, u_ptr, um_ptr, 
    //			      m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g],
    //			      &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
    if( m_ghost_points == 2 )
    {
       if( topographyExists() && g == mNumberOfGrids-1 )
	  F77_FUNC(addsgd4c,ADDSGD4C) ( &mDt, up_ptr, u_ptr, um_ptr, rho_ptr,
					m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], m_sg_str_y[g], mJ.c_ptr(),
					m_sg_corner_x[g], m_sg_corner_y[g],
					&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
       else
	  F77_FUNC(addsgd4,ADDSGD4) ( &mDt, &mGridSize[g], up_ptr, u_ptr, um_ptr, rho_ptr,
				      m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g],
				      m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
				      &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
    }
    else if(  m_ghost_points == 3 )
    {
       if( topographyExists() && g == mNumberOfGrids-1 )
	  F77_FUNC(addsgd6c,ADDSGD6C) ( &mDt, up_ptr, u_ptr, um_ptr, rho_ptr,
					m_sg_dc_x[g], m_sg_dc_y[g], m_sg_str_x[g], m_sg_str_y[g], mJ.c_ptr(),
					m_sg_corner_x[g], m_sg_corner_y[g],
					&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
       else
	  F77_FUNC(addsgd6,ADDSGD6) ( &mDt, &mGridSize[g], up_ptr, u_ptr, um_ptr, rho_ptr,
				      m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g],
				      m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
				      &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &m_supergrid_damping_coefficient );
    }
  }
}

//---------------------------------------------------------------------------
void EW::simpleAttenuation( vector<Sarray> & a_Up )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *up_ptr, cfreq, dt;
// Qs is saved in the EW object as mQs
// center frequency is called m_att_max_frecuency
// time step is called mDt  
  dt = mDt;
  cfreq = m_att_max_frequency;
  
  int g;
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    double* qs_ptr =mQs[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    F77_FUNC(satt,SATT) ( up_ptr, qs_ptr, &dt, &cfreq,
			  &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast );
  }
}

//-----------------------------------------------------------------------
void EW::enforceBCfreeAtt( vector<Sarray>& a_Up, vector<Sarray>& a_U, vector<Sarray>& a_Um, 
			   vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
			   vector<Sarray*>& a_AlphaVEp, vector<Sarray*>& a_AlphaVEm,
                           vector<double **>& a_BCForcing, double bop[5], double a_t )
{
   int sg = usingSupergrid();
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
      int ifirst = m_iStart[g];
      int ilast  = m_iEnd[g];
      int jfirst = m_jStart[g];
      int jlast  = m_jEnd[g];
      int kfirst = m_kStart[g];
      int klast  = m_kEnd[g];
      double h   = mGridSize[g];
      int topo = topographyExists() && g == mNumberOfGrids-1;
      if( m_bcType[g][4] == bStressFree && !topo )
      {
	 // Note: Only one memforce, because twilight assumes nmech=1.
	 Sarray memforce(3,ifirst,ilast,jfirst,jlast,1,1);
	 if( m_twilight_forcing )
	 {
            double om = m_twilight_forcing->m_omega;
	    double cv = m_twilight_forcing->m_c;
	    double ph = m_twilight_forcing->m_phase;
	    double* mf = memforce.c_ptr();
            int k=0; // Ghost point
	    F77_FUNC(memvarforcesurf,MEMVARFORCESURF)( &ifirst, &ilast, &jfirst, &jlast, &k, mf, &a_t, &om,
						      &cv, &ph, &mOmegaVE[0], &mDt, &h, &m_zmin[g] );
	 }
	 else
	    memforce.set_value(0.0);

	 double g1, g2, g3, r1[8], r2[8], r3[8], cof[8], acof, bcof, a4ci, b4ci, a4cj, b4cj;
	 const double i6  = 1.0/6;
	 const double d4a = 2.0/3;
	 const double d4b =-1.0/12;
	 double* forcing = a_BCForcing[g][4];
	 int ni = (ilast-ifirst+1);
	 for( int j=jfirst+2 ; j<=jlast-2 ; j++ )
	    for( int i=ifirst+2 ; i<=ilast-2 ; i++ )
	    {
	       int ind = i-ifirst + ni*(j-jfirst);
	       //               if( i==23 && j==18 )
	       //		  cout << "bforce rhs " << forcing[3*ind] << endl;
	       a4ci = a4cj = d4a;
	       b4ci = b4cj = d4b;
	       if( sg == 1 )
	       {
		  a4ci = d4a*m_sg_str_x[g][i-ifirst];
		  b4ci = d4b*m_sg_str_x[g][i-ifirst];
		  a4cj = d4a*m_sg_str_y[g][j-jfirst];
		  b4cj = d4b*m_sg_str_y[g][j-jfirst];
	       }
	       g1 = -a_Mu[g](i,j,1)*(bop[1]*a_Up[g](1,i,j,1) + bop[2]*a_Up[g](1,i,j,2) +
				     bop[3]*a_Up[g](1,i,j,3) + bop[4]*a_Up[g](1,i,j,4) +
				     a4ci*(a_Up[g](3,i+1,j,1)-a_Up[g](3,i-1,j,1))
			           + b4ci*(a_Up[g](3,i+2,j,1)-a_Up[g](3,i-2,j,1)) ) + h*forcing[3*ind];

	       g2 = -a_Mu[g](i,j,1)*( bop[1]*a_Up[g](2,i,j,1) + bop[2]*a_Up[g](2,i,j,2) +
				      bop[3]*a_Up[g](2,i,j,3) + bop[4]*a_Up[g](2,i,j,4) +
                                       a4cj*(a_Up[g](3,i,j+1,1)-a_Up[g](3,i,j-1,1))
				     + b4cj*(a_Up[g](3,i,j+2,1)-a_Up[g](3,i,j-2,1)) ) + h*forcing[3*ind+1];

	       g3 = -(2*a_Mu[g](i,j,1)+a_Lambda[g](i,j,1))*(
                                     bop[1]*a_Up[g](3,i,j,1) + bop[2]*a_Up[g](3,i,j,2)+
				     bop[3]*a_Up[g](3,i,j,3) + bop[4]*a_Up[g](3,i,j,4) ) -
		  a_Lambda[g](i,j,1)*( a4ci*(a_Up[g](1,i+1,j,1)-a_Up[g](1,i-1,j,1))
				     + b4ci*(a_Up[g](1,i+2,j,1)-a_Up[g](1,i-2,j,1)) 
				     + a4cj*(a_Up[g](2,i,j+1,1)-a_Up[g](2,i,j-1,1))
			             + b4cj*(a_Up[g](2,i,j+2,1)-a_Up[g](2,i,j-2,1)) ) + h*forcing[3*ind+2];
	       acof = a_Mu[g](i,j,1);
	       bcof = 2*a_Mu[g](i,j,1)+a_Lambda[g](i,j,1);
	       for( int a=0 ; a < m_number_mechanisms ; a++ )
	       {
		  double omdt = mOmegaVE[a]*mDt;
		  double cp = 0.5 + 1/(2*omdt) + omdt/4 + omdt*omdt/12;
		  double cm = 0.5 - 1/(2*omdt) - omdt/4 + omdt*omdt/12;

		  cof[a]= (omdt+1)/(6*cp);
		  r1[a] = (-cm*a_AlphaVEm[g][a](1,i,j,0)+(4+omdt*omdt)*i6*a_U[g](1,i,j,0)+i6*(1-omdt)*a_Um[g](1,i,j,0)+
			   memforce(1,i,j,1))/cp;
		  r2[a] = (-cm*a_AlphaVEm[g][a](2,i,j,0)+(4+omdt*omdt)*i6*a_U[g](2,i,j,0)+i6*(1-omdt)*a_Um[g](2,i,j,0)+
			   memforce(2,i,j,1))/cp;
		  r3[a] = (-cm*a_AlphaVEm[g][a](3,i,j,0)+(4+omdt*omdt)*i6*a_U[g](3,i,j,0)+i6*(1-omdt)*a_Um[g](3,i,j,0)+
			   memforce(3,i,j,1))/cp;

		  g1 = g1 + mMuVE[g][a](i,j,1)*( bop[1]*a_AlphaVEp[g][a](1,i,j,1)+bop[2]*a_AlphaVEp[g][a](1,i,j,2)+
						 bop[3]*a_AlphaVEp[g][a](1,i,j,3)+bop[4]*a_AlphaVEp[g][a](1,i,j,4)+
							    a4ci*(a_AlphaVEp[g][a](3,i+1,j,1)-a_AlphaVEp[g][a](3,i-1,j,1))
							  + b4ci*(a_AlphaVEp[g][a](3,i+2,j,1)-a_AlphaVEp[g][a](3,i-2,j,1)) +
							    bop[0]*r1[a] );
		  g2 = g2 + mMuVE[g][a](i,j,1)*( bop[1]*a_AlphaVEp[g][a](2,i,j,1)+bop[2]*a_AlphaVEp[g][a](2,i,j,2)+
						 bop[3]*a_AlphaVEp[g][a](2,i,j,3)+bop[4]*a_AlphaVEp[g][a](2,i,j,4)+
							    a4cj*(a_AlphaVEp[g][a](3,i,j+1,1)-a_AlphaVEp[g][a](3,i,j-1,1))
							  + b4cj*(a_AlphaVEp[g][a](3,i,j+2,1)-a_AlphaVEp[g][a](3,i,j-2,1)) +
							    bop[0]*r2[a] );
		  g3 = g3 + (2*mMuVE[g][a](i,j,1)+mLambdaVE[g][a](i,j,1))*(
				    bop[1]*a_AlphaVEp[g][a](3,i,j,1)+bop[2]*a_AlphaVEp[g][a](3,i,j,2)+
				    bop[3]*a_AlphaVEp[g][a](3,i,j,3)+bop[4]*a_AlphaVEp[g][a](3,i,j,4)+
				    bop[0]*r3[a]) +
				    mLambdaVE[g][a](i,j,1)*(
					  a4ci*(a_AlphaVEp[g][a](1,i+1,j,1)-a_AlphaVEp[g][a](1,i-1,j,1))
					+ b4ci*(a_AlphaVEp[g][a](1,i+2,j,1)-a_AlphaVEp[g][a](1,i-2,j,1)) 
					+ a4cj*(a_AlphaVEp[g][a](2,i,j+1,1)-a_AlphaVEp[g][a](2,i,j-1,1))
				        + b4cj*(a_AlphaVEp[g][a](2,i,j+2,1)-a_AlphaVEp[g][a](2,i,j-2,1)) );
		  acof -=    mMuVE[g][a](i,j,1)*cof[a];
		  bcof -= (2*mMuVE[g][a](i,j,1)+mLambdaVE[g][a](i,j,1))*cof[a];
	       }
	       a_Up[g](1,i,j,0) = g1/(bop[0]*acof);
	       a_Up[g](2,i,j,0) = g2/(bop[0]*acof);
	       a_Up[g](3,i,j,0) = g3/(bop[0]*bcof);
	       for( int a=0 ; a < m_number_mechanisms ; a++ )
	       {
		  a_AlphaVEp[g][a](1,i,j,0) = cof[a]*a_Up[g](1,i,j,0)+ r1[a];
		  a_AlphaVEp[g][a](2,i,j,0) = cof[a]*a_Up[g](2,i,j,0)+ r2[a];
		  a_AlphaVEp[g][a](3,i,j,0) = cof[a]*a_Up[g](3,i,j,0)+ r3[a];
	       }
	    }
      }
      if( m_bcType[g][5] == bStressFree  )
      {
         int nk=m_global_nz[g];
	 // Note: Only one memforce, because twilight assumes nmech=1.
	 Sarray memforce(3,ifirst,ilast,jfirst,jlast,1,1);
	 if( m_twilight_forcing )
	 {
            double om = m_twilight_forcing->m_omega;
	    double cv = m_twilight_forcing->m_c;
	    double ph = m_twilight_forcing->m_phase;
	    double* mf = memforce.c_ptr();
            int k=nk+1; // Ghost point
	    F77_FUNC(memvarforcesurf,MEMVARFORCESURF)( &ifirst, &ilast, &jfirst, &jlast, &k, mf, &a_t, &om,
						      &cv, &ph, &mOmegaVE[0], &mDt, &h, &m_zmin[g] );
	 }
	 else
	    memforce.set_value(0.0);

	 double g1, g2, g3, r1[8], r2[8], r3[8], cof[8], acof, bcof, a4ci, b4ci, a4cj, b4cj;
	 const double i6  = 1.0/6;
	 const double d4a = 2.0/3;
	 const double d4b =-1.0/12;
	 double* forcing = a_BCForcing[g][5];
	 int ni = (ilast-ifirst+1);
	 for( int j=jfirst+2 ; j<=jlast-2 ; j++ )
	    for( int i=ifirst+2 ; i<=ilast-2 ; i++ )
	    {
	       int ind = i-ifirst + ni*(j-jfirst);
	       a4ci = a4cj = d4a;
	       b4ci = b4cj = d4b;
	       if( sg == 1 )
	       {
		  a4ci = d4a*m_sg_str_x[g][i-ifirst];
		  b4ci = d4b*m_sg_str_x[g][i-ifirst];
		  a4cj = d4a*m_sg_str_y[g][j-jfirst];
		  b4cj = d4b*m_sg_str_y[g][j-jfirst];
	       }
	       g1 = -a_Mu[g](i,j,nk)*(-bop[1]*a_Up[g](1,i,j,nk) - bop[2]*a_Up[g](1,i,j,nk-1) 
				     -bop[3]*a_Up[g](1,i,j,nk-2) - bop[4]*a_Up[g](1,i,j,nk-3) +
				     a4ci*(a_Up[g](3,i+1,j,nk)-a_Up[g](3,i-1,j,nk))
			           + b4ci*(a_Up[g](3,i+2,j,nk)-a_Up[g](3,i-2,j,nk)) ) + h*forcing[3*ind];

	       g2 = -a_Mu[g](i,j,nk)*(-bop[1]*a_Up[g](2,i,j,nk) - bop[2]*a_Up[g](2,i,j,nk-1) -
				      bop[3]*a_Up[g](2,i,j,nk-2) - bop[4]*a_Up[g](2,i,j,nk-3) +
                                       a4cj*(a_Up[g](3,i,j+1,nk)-a_Up[g](3,i,j-1,nk))
				     + b4cj*(a_Up[g](3,i,j+2,nk)-a_Up[g](3,i,j-2,nk)) ) + h*forcing[3*ind+1];

	       g3 = -(2*a_Mu[g](i,j,nk)+a_Lambda[g](i,j,nk))*(
                                     -bop[1]*a_Up[g](3,i,j,nk) - bop[2]*a_Up[g](3,i,j,nk-1) 
				     -bop[3]*a_Up[g](3,i,j,nk-2) - bop[4]*a_Up[g](3,i,j,nk-3) ) -
		  a_Lambda[g](i,j,nk)*( a4ci*(a_Up[g](1,i+1,j,nk)-a_Up[g](1,i-1,j,nk))
				     + b4ci*(a_Up[g](1,i+2,j,nk)-a_Up[g](1,i-2,j,nk)) 
				     + a4cj*(a_Up[g](2,i,j+1,nk)-a_Up[g](2,i,j-1,nk))
			             + b4cj*(a_Up[g](2,i,j+2,nk)-a_Up[g](2,i,j-2,nk)) ) + h*forcing[3*ind+2];
	       acof = a_Mu[g](i,j,nk);
	       bcof = 2*a_Mu[g](i,j,nk)+a_Lambda[g](i,j,nk);
	       for( int a=0 ; a < m_number_mechanisms ; a++ )
	       {
		  double omdt = mOmegaVE[a]*mDt;
		  double cp = 0.5 + 1/(2*omdt) + omdt/4 + omdt*omdt/12;
		  double cm = 0.5 - 1/(2*omdt) - omdt/4 + omdt*omdt/12;

		  cof[a]= (omdt+1)/(6*cp);
		  r1[a] = (-cm*a_AlphaVEm[g][a](1,i,j,nk+1)+(4+omdt*omdt)*i6*a_U[g](1,i,j,nk+1)+i6*(1-omdt)*a_Um[g](1,i,j,nk+1)+
			   memforce(1,i,j,1))/cp;
		  r2[a] = (-cm*a_AlphaVEm[g][a](2,i,j,nk+1)+(4+omdt*omdt)*i6*a_U[g](2,i,j,nk+1)+i6*(1-omdt)*a_Um[g](2,i,j,nk+1)+
			   memforce(2,i,j,1))/cp;
		  r3[a] = (-cm*a_AlphaVEm[g][a](3,i,j,nk+1)+(4+omdt*omdt)*i6*a_U[g](3,i,j,nk+1)+i6*(1-omdt)*a_Um[g](3,i,j,nk+1)+
			   memforce(3,i,j,1))/cp;

		  g1 = g1 + mMuVE[g][a](i,j,nk)*( -bop[1]*a_AlphaVEp[g][a](1,i,j,nk)-bop[2]*a_AlphaVEp[g][a](1,i,j,nk-1)-
						   bop[3]*a_AlphaVEp[g][a](1,i,j,nk-2)-bop[4]*a_AlphaVEp[g][a](1,i,j,nk-3)+
							    a4ci*(a_AlphaVEp[g][a](3,i+1,j,nk)-a_AlphaVEp[g][a](3,i-1,j,nk))
							  + b4ci*(a_AlphaVEp[g][a](3,i+2,j,nk)-a_AlphaVEp[g][a](3,i-2,j,nk)) 
							   - bop[0]*r1[a] );
		  g2 = g2 + mMuVE[g][a](i,j,nk)*( -bop[1]*a_AlphaVEp[g][a](2,i,j,nk)-bop[2]*a_AlphaVEp[g][a](2,i,j,nk-1)-
						   bop[3]*a_AlphaVEp[g][a](2,i,j,nk-2)-bop[4]*a_AlphaVEp[g][a](2,i,j,nk-3)+
							    a4cj*(a_AlphaVEp[g][a](3,i,j+1,nk)-a_AlphaVEp[g][a](3,i,j-1,nk))
							  + b4cj*(a_AlphaVEp[g][a](3,i,j+2,nk)-a_AlphaVEp[g][a](3,i,j-2,nk)) -
							    bop[0]*r2[a] );
		  g3 = g3 + (2*mMuVE[g][a](i,j,nk)+mLambdaVE[g][a](i,j,nk))*(
				   -bop[1]*a_AlphaVEp[g][a](3,i,j,nk)-bop[2]*a_AlphaVEp[g][a](3,i,j,nk-1)
				   -bop[3]*a_AlphaVEp[g][a](3,i,j,nk-2)-bop[4]*a_AlphaVEp[g][a](3,i,j,nk-3)-
				    bop[0]*r3[a]) +
				    mLambdaVE[g][a](i,j,nk)*(
					  a4ci*(a_AlphaVEp[g][a](1,i+1,j,nk)-a_AlphaVEp[g][a](1,i-1,j,nk))
					+ b4ci*(a_AlphaVEp[g][a](1,i+2,j,nk)-a_AlphaVEp[g][a](1,i-2,j,nk)) 
					+ a4cj*(a_AlphaVEp[g][a](2,i,j+1,nk)-a_AlphaVEp[g][a](2,i,j-1,nk))
				        + b4cj*(a_AlphaVEp[g][a](2,i,j+2,nk)-a_AlphaVEp[g][a](2,i,j-2,nk)) );
		  acof -=    mMuVE[g][a](i,j,nk)*cof[a];
		  bcof -= (2*mMuVE[g][a](i,j,nk)+mLambdaVE[g][a](i,j,nk))*cof[a];
	       }
	       a_Up[g](1,i,j,nk+1) = g1/(-bop[0]*acof);
	       a_Up[g](2,i,j,nk+1) = g2/(-bop[0]*acof);
	       a_Up[g](3,i,j,nk+1) = g3/(-bop[0]*bcof);
	       for( int a=0 ; a < m_number_mechanisms ; a++ )
	       {
		  a_AlphaVEp[g][a](1,i,j,nk+1) = cof[a]*a_Up[g](1,i,j,nk+1)+ r1[a];
		  a_AlphaVEp[g][a](2,i,j,nk+1) = cof[a]*a_Up[g](2,i,j,nk+1)+ r2[a];
		  a_AlphaVEp[g][a](3,i,j,nk+1) = cof[a]*a_Up[g](3,i,j,nk+1)+ r3[a];
	       }
	    }
      }
      if( m_bcType[g][4] == bStressFree && topo && g == mNumberOfGrids-1 )
      {
	 // Note: Only one memforce, because twilight assumes nmech=1.
	 Sarray memforce(3,ifirst,ilast,jfirst,jlast,1,1);
	 double* mf = memforce.c_ptr();
	 if( m_twilight_forcing )
	 {
            double om = m_twilight_forcing->m_omega;
	    double cv = m_twilight_forcing->m_c;
	    double ph = m_twilight_forcing->m_phase;
            int k=0; // Ghost point
	    F77_FUNC(memvarforcesurfc,MEMVARFORCESURFC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &k,
							 mf, &a_t, &om, &cv, &ph, &mOmegaVE[0], &mDt, mX.c_ptr(),
							 mY.c_ptr(), mZ.c_ptr() );
	 }
	 else
	    memforce.set_value(0.0);

	 double* mu_p = a_Mu[g].c_ptr();
	 double* la_p = a_Lambda[g].c_ptr();
	 double* up_p = a_Up[g].c_ptr();
         int side = 5;
	 int nz = m_global_nz[g];
         int ghno = 0;
         char op = '-';
	 double* forcing = a_BCForcing[g][4];
	 int usesg = usingSupergrid() ? 1 : 0;

         Sarray bforcerhs(3,ifirst,ilast,jfirst,jlast,1,1);
         bforcerhs.assign(forcing);

	 F77_FUNC(addbstressc,ADDBSTRESSC)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					   &nz, up_p, mu_p, la_p, bforcerhs.c_ptr(), mMetric.c_ptr(), 
					   &side, m_sbop, &op, &ghno, &usesg, m_sg_str_x[g], m_sg_str_y[g] );
         double cof[8];
	 double* u_p      = a_U[g].c_ptr();
	 double* um_p     = a_Um[g].c_ptr();

         Sarray mubnd(ifirst,ilast,jfirst,jlast,1,1);
         Sarray lambdabnd(ifirst,ilast,jfirst,jlast,1,1);
         mubnd.set_to_zero();
	 lambdabnd.set_to_zero();
         for( int a = 0 ; a < m_number_mechanisms ; a++ )
	 {
	    double* muve_p   = mMuVE[g][a].c_ptr();
	    double* lave_p   = mLambdaVE[g][a].c_ptr();
	    double* alphap_p = a_AlphaVEp[g][a].c_ptr();
	    double* alpham_p = a_AlphaVEm[g][a].c_ptr();
	    // This function will 1.update bforcerhs, 2.update alphap, 3.accumulate mubnd and lambdabnd,
	    //  4.compute cof.
	    F77_FUNC(addbstresswresc,ADDBSTRESSWRESC)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						      &nz, alphap_p, alpham_p, muve_p, lave_p, bforcerhs.c_ptr(), 
						      u_p, um_p, mMetric.c_ptr(), &side, &mDt, &mOmegaVE[a], mf, 
						      mubnd.c_ptr(), lambdabnd.c_ptr(), m_sbop, &cof[a], &usesg,
						      m_sg_str_x[g], m_sg_str_y[g] );
	 }
         F77_FUNC(solveattfreec,SOLVEATTFREEC)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					       up_p, mu_p, la_p, mubnd.c_ptr(), lambdabnd.c_ptr(),
					       bforcerhs.c_ptr(), mMetric.c_ptr(), m_sbop, &usesg,
					       m_sg_str_x[g], m_sg_str_y[g] );
         for( int a = 0 ; a < m_number_mechanisms ; a++ )
	 {
	    double* alphap_p = a_AlphaVEp[g][a].c_ptr();
	    F77_FUNC(solveattfreeac,SOLVEATTFREEAC)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
						    alphap_p, &cof[a], up_p );
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::addAttToFreeBcForcing( vector<Sarray*>& a_AlphaVEp,
				vector<double**>& a_BCForcing, double bop[5] )
{
   int sg = usingSupergrid();
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
      int ifirst = m_iStart[g];
      int ilast  = m_iEnd[g];
      int jfirst = m_jStart[g];
      int jlast  = m_jEnd[g];
      int kfirst = m_kStart[g];
      int klast  = m_kEnd[g];
      double ih = 1.0/mGridSize[g];
      int topo=topographyExists() && g == mNumberOfGrids-1;
      if( (m_bcType[g][4] == bStressFree) && !topo )
      {
	 double a4ci, b4ci, a4cj, b4cj;
	 const double d4a = 2.0/3;
	 const double d4b =-1.0/12;
	 double* forcing = a_BCForcing[g][4];
	 int ni = (ilast-ifirst+1);
	 for( int j=jfirst+2 ; j<=jlast-2 ; j++ )
	    for( int i=ifirst+2 ; i<=ilast-2 ; i++ )
	    {
	       int ind = i-ifirst + ni*(j-jfirst);
	       a4ci = a4cj = d4a;
	       b4ci = b4cj = d4b;
	       if( sg == 1 )
	       {
		  a4ci = d4a*m_sg_str_x[g][i-ifirst];
		  b4ci = d4b*m_sg_str_x[g][i-ifirst];
		  a4cj = d4a*m_sg_str_y[g][j-jfirst];
		  b4cj = d4b*m_sg_str_y[g][j-jfirst];
	       }

               for( int a=0 ; a < m_number_mechanisms ; a++ )
	       {
		  forcing[3*ind] +=  ih*mMuVE[g][a](i,j,1)*(
			      bop[0]*a_AlphaVEp[g][a](1,i,j,0)+bop[1]*a_AlphaVEp[g][a](1,i,j,1)+
			      bop[2]*a_AlphaVEp[g][a](1,i,j,2)+bop[3]*a_AlphaVEp[g][a](1,i,j,3)+
			      bop[4]*a_AlphaVEp[g][a](1,i,j,4) +
 			        a4ci*(a_AlphaVEp[g][a](3,i+1,j,1)-a_AlphaVEp[g][a](3,i-1,j,1))
			      + b4ci*(a_AlphaVEp[g][a](3,i+2,j,1)-a_AlphaVEp[g][a](3,i-2,j,1)) );

		  forcing[3*ind+1] += ih*mMuVE[g][a](i,j,1)*(
			      bop[0]*a_AlphaVEp[g][a](2,i,j,0)+bop[1]*a_AlphaVEp[g][a](2,i,j,1)+
			      bop[2]*a_AlphaVEp[g][a](2,i,j,2)+bop[3]*a_AlphaVEp[g][a](2,i,j,3)+
			      bop[4]*a_AlphaVEp[g][a](2,i,j,4) +
 				      a4cj*(a_AlphaVEp[g][a](3,i,j+1,1)-a_AlphaVEp[g][a](3,i,j-1,1))
				    + b4cj*(a_AlphaVEp[g][a](3,i,j+2,1)-a_AlphaVEp[g][a](3,i,j-2,1)) );

		  forcing[3*ind+2] += ih*(2*mMuVE[g][a](i,j,1)+mLambdaVE[g][a](i,j,1))*(
                              bop[0]*a_AlphaVEp[g][a](3,i,j,0)+bop[1]*a_AlphaVEp[g][a](3,i,j,1)+
			      bop[2]*a_AlphaVEp[g][a](3,i,j,2)+bop[3]*a_AlphaVEp[g][a](3,i,j,3)+
			      bop[4]*a_AlphaVEp[g][a](3,i,j,4) ) +
				    ih*mLambdaVE[g][a](i,j,1)*(
					  a4ci*(a_AlphaVEp[g][a](1,i+1,j,1)-a_AlphaVEp[g][a](1,i-1,j,1))
					+ b4ci*(a_AlphaVEp[g][a](1,i+2,j,1)-a_AlphaVEp[g][a](1,i-2,j,1)) 
					+ a4cj*(a_AlphaVEp[g][a](2,i,j+1,1)-a_AlphaVEp[g][a](2,i,j-1,1))
					+ b4cj*(a_AlphaVEp[g][a](2,i,j+2,1)-a_AlphaVEp[g][a](2,i,j-2,1)) );
	       }
	    }
      }
      if( m_bcType[g][5] == bStressFree )
      {
	 double a4ci, b4ci, a4cj, b4cj;
	 const double d4a = 2.0/3;
	 const double d4b =-1.0/12;
	 double* forcing = a_BCForcing[g][5];
	 int ni = (ilast-ifirst+1);
	 for( int j=jfirst+2 ; j<=jlast-2 ; j++ )
	    for( int i=ifirst+2 ; i<=ilast-2 ; i++ )
	    {
	       int ind = i-ifirst + ni*(j-jfirst);
	       a4ci = a4cj = d4a;
	       b4ci = b4cj = d4b;
	       if( sg == 1 )
	       {
		  a4ci = d4a*m_sg_str_x[g][i-ifirst];
		  b4ci = d4b*m_sg_str_x[g][i-ifirst];
		  a4cj = d4a*m_sg_str_y[g][j-jfirst];
		  b4cj = d4b*m_sg_str_y[g][j-jfirst];
	       }
               int nk=m_global_nz[g];
               for( int a=0 ; a < m_number_mechanisms ; a++ )
	       {
		  forcing[3*ind] +=  ih*mMuVE[g][a](i,j,nk)*(
 		          -(  bop[0]*a_AlphaVEp[g][a](1,i,j,nk+1)+bop[1]*a_AlphaVEp[g][a](1,i,j,nk)+
			      bop[2]*a_AlphaVEp[g][a](1,i,j,nk-1)+bop[3]*a_AlphaVEp[g][a](1,i,j,nk-2)+
			      bop[4]*a_AlphaVEp[g][a](1,i,j,nk-3) ) +
 			        a4ci*(a_AlphaVEp[g][a](3,i+1,j,nk)-a_AlphaVEp[g][a](3,i-1,j,nk))
			      + b4ci*(a_AlphaVEp[g][a](3,i+2,j,nk)-a_AlphaVEp[g][a](3,i-2,j,nk)) );

		  forcing[3*ind+1] += ih*mMuVE[g][a](i,j,nk)*(
		          -( bop[0]*a_AlphaVEp[g][a](2,i,j,nk+1)+bop[1]*a_AlphaVEp[g][a](2,i,j,nk)+
			     bop[2]*a_AlphaVEp[g][a](2,i,j,nk-1)+bop[3]*a_AlphaVEp[g][a](2,i,j,nk-2)+
			     bop[4]*a_AlphaVEp[g][a](2,i,j,nk-3) ) +
 				      a4cj*(a_AlphaVEp[g][a](3,i,j+1,nk)-a_AlphaVEp[g][a](3,i,j-1,nk))
				    + b4cj*(a_AlphaVEp[g][a](3,i,j+2,nk)-a_AlphaVEp[g][a](3,i,j-2,nk)) );

		  forcing[3*ind+2] += ih*(2*mMuVE[g][a](i,j,nk)+mLambdaVE[g][a](i,j,nk))*(
			  -( bop[0]*a_AlphaVEp[g][a](3,i,j,nk+1)+bop[1]*a_AlphaVEp[g][a](3,i,j,nk)+
			     bop[2]*a_AlphaVEp[g][a](3,i,j,nk-1)+bop[3]*a_AlphaVEp[g][a](3,i,j,nk-2)+
			     bop[4]*a_AlphaVEp[g][a](3,i,j,nk-3)) ) +
				    ih*mLambdaVE[g][a](i,j,nk)*(
					  a4ci*(a_AlphaVEp[g][a](1,i+1,j,nk)-a_AlphaVEp[g][a](1,i-1,j,nk))
					+ b4ci*(a_AlphaVEp[g][a](1,i+2,j,nk)-a_AlphaVEp[g][a](1,i-2,j,nk)) 
					+ a4cj*(a_AlphaVEp[g][a](2,i,j+1,nk)-a_AlphaVEp[g][a](2,i,j-1,nk))
					+ b4cj*(a_AlphaVEp[g][a](2,i,j+2,nk)-a_AlphaVEp[g][a](2,i,j-2,nk)) );
	       }
	    }
      }
      if( (m_bcType[g][4] == bStressFree) && topo )
      {
         int nz = m_global_nz[g];
         char op = '+';
	 int side = 5, ghyes=1;
         int usesg = usingSupergrid() ? 1:0;
	 double* forcing = a_BCForcing[g][4];
         for( int a=0 ; a < m_number_mechanisms ; a++ )
	 {
	    double* muve_p  = mMuVE[g][a].c_ptr();
	    double* lave_p  = mLambdaVE[g][a].c_ptr();
	    double* alpha_p = a_AlphaVEp[g][a].c_ptr();
	    F77_FUNC(addbstressc,ADDBSTRESSC)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
					      &nz, alpha_p, muve_p, lave_p, forcing, mMetric.c_ptr(),
					      &side, bop, &op, &ghyes, &usesg, m_sg_str_x[g], m_sg_str_y[g] );
	 }
      }
   }
}
