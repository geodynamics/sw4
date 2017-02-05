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
   void tw_aniso_free_surf_z(int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                             int kz, double t, double om, double cv, double ph, double omm, double* phc,
                             double* bforce, double h, double zmin );

   void twfrsurfz_wind( int *ifirst, int *ilast, int *jfirst, int *jlast, int *kfirst, int *klast,
                        double *h, int* kz, double *t, double *omega, double *c, double *phase, double *bforce,
                        double *mu, double *lambda, double *zmin,
                        int *i1, int *i2, int *j1, int *j2 );

   void twfrsurfzsg_wind( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                          double h, int kz, double t, double omega, double c, double phase, double omstrx, double omstry,
                          double *bforce, double *mu, double *lambda, double zmin,
                          int i1, int i2, int j1, int j2 );


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
   void addsg4wind( double *dt, double *h, double *up, double *u, double *um, double *rho,
                    double *dcx, double *dcy, double*  dcz, double *strx, double *stry, double *strz,
                    double *cox, double *coy, double *coz, 
                    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, double beta,
                    int kupb, int kupe, int kwindb, int kwinde );
   
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
   void F77_FUNC(bcfreesurfcurvani,BCFREESURFCURVANI)( int*, int*, int*, int*, int*, int*, int*, double*, double*, int*,
						       double*, double*, double*, double*, double* );
   void F77_FUNC(twilightfortwind,TWILIGHTFORTWIND)( int*, int*, int*, int*, int*, int*, double*, double*,
						     double*, double*, double*, double*, double*, int*, int*,
						     int*, int*, int*, int* );
   void rhs4th3fortwind( int*, int*, int*, int*, int*, int*, int*, int*, double*,
                         double*, double*,double*, double*, double*, double*, double*,
                         double*, double*, double*, char*, int*, int*, int*, int* );
   void forcingfort( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, 
                double*, double*, double*, double*, double*, double*, double*, double* );
   void forcingfortsg(int*, int*, int*, int*, int*, 
                      int*, double*, double*, double*, double*, double*, double*, double*, 
                      double*, double*, double*, double*, double*,double*,double*,double* );
   void F77_FUNC(forcingttfort,FORCINGTTFORT)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, 
					       double*, double*, double*, double*, double*, double*, double*, double* );
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
    
// output flags and settings that affect the run
   if( proc_zero() && mVerbose >= 3 )
   {
      printf("\nReporting SW4 internal flags and settings:\n");
      printf("m_testing=%s, twilight=%s, point_source=%s, moment_test=%s, energy_test=%s," 
             "rayleigh_test=%s\n",
             m_testing?"yes":"no",
             m_twilight_forcing?"yes":"no",
             m_point_source_test?"yes":"no",
             m_moment_test?"yes":"no",
             m_energy_test?"yes":"no",
             m_lamb_test?"yes":"no",
             m_rayleigh_wave_test?"yes":"no");
      printf("m_use_supergrid=%s\n", m_use_supergrid?"yes":"no");
      printf("End report settings\n\n");
   }
   
  if ( !mQuiet && proc_zero() )
    cout << "  Begin time stepping..." << endl;

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

// Grid refinement interface conditions:
    if (mOrder == 2)
       enforceIC2( Up, U, Um, t, point_sources );
    else
       enforceIC( Up, U, Um, t, true, point_sources );

    
    time_measure[3] = time_measure[4] = MPI_Wtime();

//
// *** 2nd order in TIME
//
    if (mOrder == 2)
    {
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
       
    } // end mOrder == 2
//
// *** 4th order in time ***
//
    else if (mOrder == 4)
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
       enforceIC( Up, U, Um, t, false, point_sources );

    }// end if mOrder == 4
    
    
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
// for periodic bc, a_BCForcing[g][s] == NULL, so you better not access
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
       if( topographyExists() && g == mNumberOfGrids-1 && m_bcType[g][4] == bStressFree )
       {
	  int fside = 5;
	  double* cc_ptr = mCcurv.c_ptr();
          F77_FUNC(bcfreesurfcurvani,BCFREESURFCURVANI)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
							&nz, u_ptr, cc_ptr, &fside, m_sbop, bforce_side4_ptr,
                                                        bforce_side5_ptr, m_sg_str_x[g], m_sg_str_y[g] );
       }

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
  update_curvilinear_cartesian_interface( a_U );
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
//--------------------Mesh refinement interface condition for 4th order predictor-corrector scheme----------------------------
void EW::enforceIC( vector<Sarray>& a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		    double time, bool predictor, vector<GridPointSource*> point_sources )
{
   for( int g = 0 ; g < mNumberOfCartesianGrids-1 ; g++ )
   {
      // Interpolate between g and g+1, assume factor 2 refinement with at least three ghost points
      VERIFY2( m_ghost_points >= 3,
		  "enforceIC Error: "<<
                  "Number of ghost points must be three or more, not " << m_ghost_points );

      Sarray Unextf, Unextc, Bf, Bc, Uf_tt, Uc_tt;
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
// to compute the corrector we need the acceleration in the vicinity of the interface
      Uf_tt.define(3,ibf,ief,jbf,jef,kf-7,kf+1);
      Uc_tt.define(3,ibc,iec,jbc,jec,kc-1,kc+7);

    // Set ghost point values that are unknowns when solving the interface condition to zero. Assume that Dirichlet data
    // are already set on ghost points on the (supergrid) sides, which are not treated as unknown variables.
      dirichlet_hom_ic( a_Up[g+1], g+1, kf+1, true ); // inside=true
      dirichlet_hom_ic( a_Up[g], g, kc-1, true );
      if( m_doubly_periodic )
      {
	 dirichlet_hom_ic( a_Up[g+1], g+1, kf+1, false ); // inside=false
	 dirichlet_hom_ic( a_Up[g], g, kc-1, false );
      }

      if( predictor ) // In the predictor step, (Unextc, Unextf) represent the displacement after the corrector step
      {
// TEST: compute_preliminary_corrector by first assigning exact ghost point values to Up; inspect Unextf & Unextc
// alphave is only used in visco-elastic mode
// source is not used with twilight mode
         // vector<Sarray*> alphave;
         // vector<Source*> sources;
         // exactSol( time+mDt, a_Up, alphave, sources );

	 compute_preliminary_corrector( a_Up[g+1], a_U[g+1], a_Um[g+1], Uf_tt, Unextf, g+1, kf, time, point_sources );
         compute_preliminary_corrector( a_Up[g], a_U[g], a_Um[g], Uc_tt, Unextc, g, kc, time, point_sources );
	 if( !m_doubly_periodic )
	 {
// dirichlet conditions for Unextc in super-grid layer at time t+dt
            dirichlet_LRic( Unextc, g, kc, time+mDt, 1 ); 
	    // dirichlet_LRic( Unextc, g, kc, t+mDt, 0 );
// only enable for testing
//	    dirichlet_LRic( Unextf, g+1, kf, time+mDt, 1 );
	 }
// // save Unextf on file (single proc)
//          FILE *fp;
//          char fname[100];
//          int j;
//          double x;
         
//          j = 0.5*(m_jStart[g+1] + m_jEnd[g+1]);
//          sprintf(fname,"ufi-j=%d.txt", j);
//          fp = fopen(fname,"w");
//          for (int i=m_iStart[g+1]; i<=m_iEnd[g+1]; i++)
//          {
//             x = (i-1)*mGridSize[g+1];
//             fprintf(fp,"%e %e %e %e\n", x, Unextf(1,i,j,kf), Unextf(2,i,j,kf), Unextf(3,i,j,kf));
//          }
//          fclose(fp);
//          CHECK_INPUT(false," controlled termination");
// // end test         
      }
      else // In the corrector step, (Unextc, Unextf) represent the displacement after next predictor step
      {
	 compute_preliminary_predictor( a_Up[g+1], a_U[g+1], Unextf, g+1, kf, time+mDt, point_sources );
	 compute_preliminary_predictor( a_Up[g], a_U[g], Unextc, g, kc, time+mDt, point_sources );

	 if( !m_doubly_periodic )
	 {
// dirichlet conditions for Unextc in super-grid layer at time t+2*dt
            dirichlet_LRic( Unextc, g, kc, time+2*mDt, 1 ); 
	    // dirichlet_LRic( Unextc, g, kc, t+2*mDt, 0 );
	    // dirichlet_LRic( Unextf, g+1, kf, t+2*mDt, 0 );
	 }
      }
      compute_icstresses( a_Up[g+1], Bf, g+1, kf, m_sg_str_x[g+1], m_sg_str_y[g+1] );
      compute_icstresses( a_Up[g], Bc, g, kc, m_sg_str_x[g], m_sg_str_y[g] );

// from enforceIC2()
      if( !m_doubly_periodic )
      {
//  dirichlet condition for Bf in the super-grid layer at time t+dt (also works with twilight)
         dirichlet_LRstress( Bf, g+1, kf, time+mDt, 1 ); 
      }
      
      communicate_array_2d( Unextf, g+1, kf );
      communicate_array_2d( Unextc, g, kc );
      communicate_array_2d( Bf, g+1, kf );
      communicate_array_2d( Bc, g, kc );

      // Up to here, interface stresses and displacement (Bc,Bf) and (Unextc, Unextf) were computed with correct ghost point values
      // in the corners.
      // In the following iteration we use centered formulas for interpolation and restriction, all the way up to the Dirichlet boundaries.
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

      double cof = predictor ? 12 : 1;
      int is_periodic[2] ={0,0};
      if( m_doubly_periodic )
	 is_periodic[0] = is_periodic[1] = 1;
      
      consintp( a_Up[g+1], Unextf, Bf, mMu[g+1], mLambda[g+1], mRho[g+1], mGridSize[g+1],
		a_Up[g],   Unextc, Bc, mMu[g],   mLambda[g],   mRho[g],   mGridSize[g], 
		cof, g, g+1, is_periodic );
      //      CHECK_INPUT(false," controlled termination");

      // Finally, restore the ghost point values on the sides of the domain.
      //
      // Note: these ghost point values might never be used
      //
      if( !m_doubly_periodic )
      {
	 // dirichlet_LRic( a_Up[g+1], g+1, kf+1, t+mDt, 0 );
	 // dirichlet_LRic( a_Up[g], g, kc-1, t+mDt, 0 );
         dirichlet_LRic( a_Up[g+1], g+1, kf+1, time+mDt, 1 );
         dirichlet_LRic( a_Up[g], g, kc-1, time+mDt, 1 );
      }
   }
}

//-----------------------Special case for 2nd order time stepper----------------------------------------------------
void EW::enforceIC2( vector<Sarray>& a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		    double time, vector<GridPointSource*> point_sources )
{
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
      compute_preliminary_predictor( a_Up[g+1], a_U[g+1], Unextf, g+1, kf, time+mDt, point_sources );
      compute_preliminary_predictor( a_Up[g], a_U[g], Unextc, g, kc, time+mDt, point_sources );

      if( !m_doubly_periodic )
      {
// dirichlet conditions for Unextc in super-grid layer at time t+2*dt
         dirichlet_LRic( Unextc, g, kc, time+2*mDt, 1 ); 
      }

//  compute contribution to the normal stresses (Bc, Bf) from the interior grid points in Up
      compute_icstresses( a_Up[g+1], Bf, g+1, kf, m_sg_str_x[g+1], m_sg_str_y[g+1] );
      compute_icstresses( a_Up[g], Bc, g, kc, m_sg_str_x[g], m_sg_str_y[g] );

      if( !m_doubly_periodic )
      {
//  dirichlet condition for Bf in the super-grid layer at time t+dt (also works with twilight)
         dirichlet_LRstress( Bf, g+1, kf, time+mDt, 1 ); 
      }
      
      communicate_array_2d( Unextf, g+1, kf );
      communicate_array_2d( Unextc, g, kc );
      communicate_array_2d( Bf, g+1, kf );
      communicate_array_2d( Bc, g, kc );

      // Up to here, interface stresses and displacement (Bc,Bf) and (Unextc, Unextf) were computed with correct ghost point values
      // in the corners.
      // In the following iteration we use centered formulas for interpolation and restriction, all the way up to the Dirichlet boundaries.
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

      double cof = predictor ? 12 : 1;
      int is_periodic[2] ={0,0};
      if( m_doubly_periodic )
	 is_periodic[0] = is_periodic[1] = 1;

// Iteratively determine the ghost point values in Up to satisfy the jump conditions
      consintp( a_Up[g+1], Unextf, Bf, mMu[g+1], mLambda[g+1], mRho[g+1], mGridSize[g+1],
		a_Up[g],   Unextc, Bc, mMu[g],   mLambda[g],   mRho[g],   mGridSize[g], 
		cof, g, g+1, is_periodic);
      //      CHECK_INPUT(false," controlled termination");

      // Finally, restore the corner ghost point values (above and below) the Dirichlet sides of the domain.
      //
      // Note: these ghost point values might never be used
      //
      if( !m_doubly_periodic )
      {
         dirichlet_LRic( a_Up[g+1], g+1, kf+1, time+mDt, 1 );
         dirichlet_LRic( a_Up[g], g, kc-1, time+mDt, 1 );
      }
            
// FROM WPP
   // // add -(1/rho)*L(alpha) to forcing. Note: forcing array will be modified
   // // however this is just a temporary array in the calling routine.
   // for( int a = 0 ; a < m_number_mechanisms ; a++ )
   // {
   //    double* alpha_p  = alpha_a[a].c_ptr();
   //    double* alphaf_p = alphaf_a[a].c_ptr();
   //    double* muve_p      = mMuVE[g][a].c_ptr();
   //    double* lambdave_p  = mLambdaVE[g][a].c_ptr();
   //    double* muvef_p     = mMuVE[g+1][a].c_ptr();
   //    double* lambdavef_p = mLambdaVE[g+1][a].c_ptr();
   //    // f := f - L_a(alpha)/rho
   //    F77_FUNC(cikplaneatt,CIKPLANEATT)( &ni, &nj, &nk, alpha_p, f_p, &side_low,
   //      			         &hc, muve_p, lambdave_p, rho_p );
   //    F77_FUNC(cikplaneatt,CIKPLANEATT)( &nif, &njf, &nkf, alphaf_p, ff_p, &side_high,
   //      			         &hf, muvef_p, lambdavef_p, rhof_p );
   // }

      //      if( predictor )
      //      {
      //	 compute_preliminary_corrector( a_Up[g+1], a_U[g+1], a_Um[g+1], Unextf, g+1, kf, t, point_sources );
      //         compute_preliminary_corrector( a_Up[g], a_U[g], a_Um[g], Unextc, g, kc, t, point_sources );
      //	 dirichlet_LRic( Unextc, g, kc, t+mDt, 0 );
      //      }
      //     else
      //      {
      //	 compute_preliminary_predictor( a_Up[g+1], a_U[g+1], Unextf, g+1, kf, t+mDt, point_sources );
      //	 compute_preliminary_predictor( a_Up[g], a_U[g], Unextc, g, kc, t+mDt, point_sources );
      //	 dirichlet_LRic( Unextc, g, kc, t+2*mDt, 0 );
      //      }
      //      check_corrector( a_Up[g+1], a_Up[g], Unextf, Unextc, kf, kc );
   }
}


//-----------------------------------------------------------------------
void EW::check_corrector( Sarray& Uf, Sarray& Uc, Sarray& Unextf, Sarray& Unextc, int kf, int kc )
{
   int ic =2, jc=4, c=1;
   int i = 2*ic-1, j=2*jc-1;
   //   cout <<"check " << Uf(c,i,j,kf) << " " << Uc(c,ic,jc,kc) << " " << Uf(c,i,j,kf)-Uc(c,ic,jc,kc) << endl;
   //   cout << "  next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << endl;
   i = 2*ic;
   //   j = 2*jc;
   double pci   = (9*(Uc(c,ic,jc,kc)+Uc(c,ic+1,jc,kc))-(Uc(c,ic-1,jc,kc)+Uc(c,ic+2,jc,kc)))/16;
   double pcim  = (9*(Uc(c,ic,jc-1,kc)+Uc(c,ic+1,jc-1,kc))-(Uc(c,ic-1,jc-1,kc)+Uc(c,ic+2,jc-1,kc)))/16;
   double pcip  = (9*(Uc(c,ic,jc+1,kc)+Uc(c,ic+1,jc+1,kc))-(Uc(c,ic-1,jc+1,kc)+Uc(c,ic+2,jc+1,kc)))/16;
   double pcipp = (9*(Uc(c,ic,jc+2,kc)+Uc(c,ic+1,jc+2,kc))-(Uc(c,ic-1,jc+2,kc)+Uc(c,ic+2,jc+2,kc)))/16;
   double pc = ( 9*(pci+pcip)-(pcim+pcipp))/16;

   double pcj = (9*(Uc(c,ic,jc,kc)+Uc(c,ic,jc+1,kc))-(Uc(c,ic,jc-1,kc)+Uc(c,ic,jc+2,kc)))/16;
   cout <<"check " << Uf(c,i,j,kf) << " " << Uc(c,ic,jc,kc) << " " << pci << " " << Uf(c,i,j,kf)-pci << endl;
   //   cout << "  next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << endl;
   double pcni   = (9*(Unextc(c,ic,jc,kc)+Unextc(c,ic+1,jc,kc))-(Unextc(c,ic-1,jc,kc)+Unextc(c,ic+2,jc,kc)))/16;
     cout << "  check next " << Unextf(c,i,j,kf) << " " << Unextc(c,ic,jc,kc) << " " << pcni << " " << Unextf(c,i,j,kf)-pcni << endl;
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
   // zero out all ghost points
   if( !inner )
   {
      // Outer layer of non-unknown ghost points
      if( m_iStartInt[g] == 1 )
      {
      // low i-side
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= 0 ; i++ )
	       for( int c=1 ; c <= U.m_nc ; c++ )
		  U(c,i,j,k) = 0;
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iEndInt[g]+1 ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= U.m_nc ; c++ )
		  U(c,i,j,k) = 0;
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	 for( int j=m_jStart[g] ; j <= 0 ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= U.m_nc ; c++ )
		  U(c,i,j,k) = 0;
      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	 for( int j=m_jEndInt[g]+1 ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= U.m_nc ; c++ )
		  U(c,i,j,k) = 0;
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
      for( int j=jb ; j <= je ; j++ )
	 for( int i=ib ; i <= ie ; i++ )
	    for( int c=1 ; c <= U.m_nc ; c++ )
	       U(c,i,j,k) = 0;
   }
}

//-----------------------------------------------------------------------
void EW::dirichlet_twilight_ic( Sarray& U, int g, int kic, double t )
{
   // assign exact solution at all ghost points with a given 'k'-index
   if (m_twilight_forcing)
   {
      int i1 = m_iStart[g], i2=m_iEnd[g];
      int j1 = m_jStart[g], j2=m_jEnd[g];
      int kdb=U.m_kb, kde=U.m_ke;
      double om = m_twilight_forcing->m_omega;
      double ph = m_twilight_forcing->m_phase;
      double cv = m_twilight_forcing->m_c;
      double h  = mGridSize[g];
      double* u_ptr = U.c_ptr();
      F77_FUNC(twilightfortwind,TWILIGHTFORTWIND)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
                                                   &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
                                                   &h, &m_zmin[g],
                                                   &i1, &i2, &j1, &j2, &kic, &kic );      
   }
   
}

//-----------------------------------------------------------------------
void EW::dirichlet_LRic( Sarray& U, int g, int kic, double t, int adj )
{
   // Put back exact solution at the ghost points that don't participate in the interface conditions, i.e. at supergrid points
   //   int k = upper ? 0 : m_global_nz[g]+1;
   //
   // set adj= 0 for ghost pts + boundary pt
   //          1 for only ghost pts.

   int kdb=U.m_kb, kde=U.m_ke;
   if( !m_twilight_forcing )
   {
      if( m_iStartInt[g] == 1 )
      {
	 // low i-side
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= 1-adj ; i++ )
	       for( int c=1 ; c <= U.m_nc ; c++ )
		  U(c,i,j,kic) = 0;
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iEndInt[g]+adj ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= U.m_nc ; c++ )
		  U(c,i,j,kic) = 0;
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	 for( int j=m_jStart[g] ; j <= 1-adj ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= U.m_nc ; c++ )
		  U(c,i,j,kic) = 0;
      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	 for( int j=m_jEndInt[g]+adj ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= U.m_nc ; c++ )
		  U(c,i,j,kic) = 0;
      }
   }
   else
   {
      double om = m_twilight_forcing->m_omega;
      double ph = m_twilight_forcing->m_phase;
      double cv = m_twilight_forcing->m_c;
      double h  = mGridSize[g];
      double* u_ptr = U.c_ptr();
      if( m_iStartInt[g] == 1 )
      {
	 // low i-side
	 int i1 = m_iStart[g], i2=m_iStartInt[g]-adj;
	 int j1 = m_jStart[g], j2=m_jEnd[g];
	 F77_FUNC(twilightfortwind,TWILIGHTFORTWIND)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
						      &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
						      &h, &m_zmin[g],
				                      &i1, &i2, &j1, &j2, &kic, &kic );
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	 int i1 = m_iEndInt[g]+adj, i2=m_iEnd[g];
	 int j1 = m_jStart[g], j2=m_jEnd[g];
	 F77_FUNC(twilightfortwind,TWILIGHTFORTWIND)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
						      &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
						      &h, &m_zmin[g],
						      &i1, &i2, &j1, &j2, &kic, &kic );
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	 int i1 = m_iStart[g], i2=m_iEnd[g];
	 int j1 = m_jStart[g], j2=m_jStartInt[g]-adj;
	 F77_FUNC(twilightfortwind,TWILIGHTFORTWIND)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
						      &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
						      &h, &m_zmin[g],
						      &i1, &i2, &j1, &j2, &kic, &kic );
      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	 int i1 = m_iStart[g], i2=m_iEnd[g];
	 int j1 = m_jEndInt[g]+adj, j2=m_jEnd[g];
	 F77_FUNC(twilightfortwind,TWILIGHTFORTWIND)( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
						      &kdb, &kde, u_ptr, &t, &om, &cv, &ph, 
						      &h, &m_zmin[g],
						      &i1, &i2, &j1, &j2, &kic, &kic );
      }
   }      
}

//-----------------------------------------------------------------------
void EW::dirichlet_LRstress( Sarray& B, int g, int kic, double t, int adj )
{
   // Exact stresses at the ghost points that don't participate in the interface conditions,
   // i.e. at supergrid points
   //   int k = upper ? 0 : m_global_nz[g]+1;
   //
   // set adj= 0 for ghost pts + boundary pt
   //          1 for only ghost pts.

   int kdb=B.m_kb, kde=B.m_ke;
//   printf("dirichlet_LRstress> kdb=%d, kde=%d\n", kdb, kde);
   
   if( !m_twilight_forcing )
   {
      if( m_iStartInt[g] == 1 )
      {
	 // low i-side
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= 1-adj ; i++ )
	       for( int c=1 ; c <= B.m_nc ; c++ )
		  B(c,i,j,kic) = 0;
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iEndInt[g]+adj ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= B.m_nc ; c++ )
		  B(c,i,j,kic) = 0;
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	 for( int j=m_jStart[g] ; j <= 1-adj ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= B.m_nc ; c++ )
		  B(c,i,j,kic) = 0;
      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	 for( int j=m_jEndInt[g]+adj ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       for( int c=1 ; c <= B.m_nc ; c++ )
		  B(c,i,j,kic) = 0;
      }
   }
   else
   {
// dbg
      if (false && g==0)
      {
         int i=0, j=25, k=1;
         double x=(i-1)*mGridSize[g];
         double y=(j-1)*mGridSize[g];
         double z=(k-1)*mGridSize[g] + m_zmin[g];

         printf("3: x=%e, y=%e, z=%e, mu=%e\n", x, y, z, mMu[g](i,j,k));
      }
// get array pointers for fortran
      double* mu_ptr    = mMu[g].c_ptr();
      double* la_ptr    = mLambda[g].c_ptr();
      double om = m_twilight_forcing->m_omega;
      double ph = m_twilight_forcing->m_phase;
      double cv = m_twilight_forcing->m_c;
      double h  = mGridSize[g];
      double* b_ptr = B.c_ptr();
      double omstrx = m_supergrid_taper_x[g].get_tw_omega();
      double omstry = m_supergrid_taper_y[g].get_tw_omega();

      if( m_iStartInt[g] == 1 )
      {
	 // low i-side
	 int i1 = m_iStart[g], i2=m_iStartInt[g]-adj;
	 int j1 = m_jStart[g], j2=m_jEnd[g];
         if( usingSupergrid() )
         {
            twfrsurfzsg_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
         }
         else
         {
            twfrsurfz_wind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
                            &h, &kic, &t, &om, &cv, &ph, b_ptr, mu_ptr, la_ptr, &m_zmin[g], &i1, &i2, &j1, &j2 );
         }
         
      }
      if( m_iEndInt[g] == m_global_nx[g] )
      {
	 // high i-side
	 int i1 = m_iEndInt[g]+adj, i2=m_iEnd[g];
	 int j1 = m_jStart[g], j2=m_jEnd[g];
         if( usingSupergrid() )
         {
            twfrsurfzsg_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
         }
         else
         {
            twfrsurfz_wind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
                            &h, &kic, &t, &om, &cv, &ph, b_ptr, 
                            mu_ptr, la_ptr, &m_zmin[g],
                            &i1, &i2, &j1, &j2 );
         }
      }
      if( m_jStartInt[g] == 1 )
      {
	 // low j-side
	 int i1 = m_iStart[g], i2=m_iEnd[g];
	 int j1 = m_jStart[g], j2=m_jStartInt[g]-adj;
         if( usingSupergrid() )
         {
            twfrsurfzsg_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
         }
         else
         {
            twfrsurfz_wind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
                            &h, &kic, &t, &om, &cv, &ph, b_ptr, 
                            mu_ptr, la_ptr, &m_zmin[g],
                            &i1, &i2, &j1, &j2 );
         }
      }
      if( m_jEndInt[g] == m_global_ny[g] )
      {
	 // high j-side
	 int i1 = m_iStart[g], i2=m_iEnd[g];
	 int j1 = m_jEndInt[g]+adj, j2=m_jEnd[g];
         if( usingSupergrid() )
         {
            twfrsurfzsg_wind( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                              h, kic, t, om, cv, ph, omstrx, omstry,
                              b_ptr, mu_ptr, la_ptr, m_zmin[g], i1, i2, j1, j2 );
         }
         else
         {
            twfrsurfz_wind( &m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
                            &h, &kic, &t, &om, &cv, &ph, b_ptr, 
                            mu_ptr, la_ptr, &m_zmin[g],
                            &i1, &i2, &j1, &j2 );
         }
      }
   }      
}

//-----------------------------------------------------------------------
void EW::gridref_initial_guess( Sarray& u, int g, bool upper )
{
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

   for( int j=jb ; j <= je ; j++ )
      for( int i=ib ; i <= ie ; i++ )
      {
	 u(1,i,j,k) = u(1,i,j,k+s);
	 u(2,i,j,k) = u(2,i,j,k+s);
	 u(3,i,j,k) = u(3,i,j,k+s);
      }
}

//-----------------------------------------------------------------------
void EW::compute_preliminary_corrector( Sarray& a_Up, Sarray& a_U, Sarray& a_Um, Sarray& Utt, Sarray& Unext,
					int g, int kic, double t, vector<GridPointSource*> point_sources )
{
   double idt2 = 1/(mDt*mDt);
   // to evaluate L(Up_tt) for k=kic, we need Up(k) in the vicinity of the interface
   // Note: Utt is needed at all points (interior + ghost) to evaluate L(Utt) in all interior points
   for( int k=Utt.m_kb ; k <= Utt.m_ke ; k++ )
      for( int j=Utt.m_jb ; j <= Utt.m_je ; j++ )
	 for( int i=Utt.m_ib ; i <= Utt.m_ie ; i++ )
	 {
	    Utt(1,i,j,k) = idt2*(a_Up(1,i,j,k)-2*a_U(1,i,j,k)+a_Um(1,i,j,k));
	    Utt(2,i,j,k) = idt2*(a_Up(2,i,j,k)-2*a_U(2,i,j,k)+a_Um(2,i,j,k));
	    Utt(3,i,j,k) = idt2*(a_Up(3,i,j,k)-2*a_U(3,i,j,k)+a_Um(3,i,j,k));
	 }
// all points (for mMu, mLambda
   int ib=m_iStart[g], jb=m_jStart[g], kb=m_kStart[g];
   int ie=m_iEnd[g], je=m_jEnd[g], ke=m_kEnd[g];
// k-indices for Utt
   int kbu = Utt.m_kb, keu= Utt.m_ke;

   // Compute L(Utt) at k=kic.
   char op='=';
   int nz = m_global_nz[g];
   Sarray Lutt(3,ib,ie,jb,je,kic,kic);
// Note: 6 first arguments of the function call:
// (ib,ie), (jb,je), (kb,ke) is the declared size of mMu and mLambda in the (i,j,k)-directions, respectively
   rhs4th3fortwind( &ib, &ie, &jb, &je, &kb, &ke, &nz, m_onesided[g], m_acof,
					      m_bope, m_ghcof, Lutt.c_ptr(), Utt.c_ptr(), mMu[g].c_ptr(),
					      mLambda[g].c_ptr(), &mGridSize[g], m_sg_str_x[g], m_sg_str_y[g],
					      m_sg_str_z[g], &op, &kbu, &keu, &kic, &kic );
// Note: 4 last arguments of the above function call:
// (kbu,keu) is the declared size of Utt in the k-direction
// (kic,kic) is the declared size of Lutt in the k-direction

   // Compute forcing_{tt} at k=kic
   Sarray force(3,ib,ie,jb,je,kic,kic);
   if( m_twilight_forcing )
   {
      double om = m_twilight_forcing->m_omega;
      double ph = m_twilight_forcing->m_phase;
      double cv = m_twilight_forcing->m_c;
      double omm= m_twilight_forcing->m_momega;
      double phm= m_twilight_forcing->m_mphase;
      double amprho   = m_twilight_forcing->m_amprho;
      double ampmu    = m_twilight_forcing->m_ampmu;
      double amplambda= m_twilight_forcing->m_amplambda;
      F77_FUNC(forcingttfort,FORCINGTTFORT)( &ib, &ie, &jb, &je, &kic, &kic, force.c_ptr(), &t, &om, &cv, &ph, &omm, &phm,
		 &amprho, &ampmu, &amplambda, &mGridSize[g], &m_zmin[g] );
   }
   else if( m_rayleigh_wave_test || m_energy_test )
      force.set_to_zero();
   else
   {
     // Default: m_point_source_test, m_lamb_test or full seismic case
      force.set_to_zero();
      for( int s = 0 ; s < point_sources.size() ; s++ )
      {
	 if( point_sources[s]->m_grid == g && point_sources[s]->m_k0 == kic )
	 {
	    double fxyz[3];
	    point_sources[s]->getFxyztt(t,fxyz);
	    force(1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
	    force(2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
	    force(3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
	 }
      }
   }

   double cof = mDt*mDt*mDt*mDt/12.0;
   // for( int j=Unext.m_jb ; j <= Unext.m_je ; j++ )
   //    for( int i=Unext.m_ib ; i <= Unext.m_ie ; i++ )
   for( int j=Unext.m_jb+2 ; j <= Unext.m_je-2 ; j++ )
      for( int i=Unext.m_ib+2 ; i <= Unext.m_ie-2 ; i++ )
      {
	 double irho=cof/mRho[g](i,j,kic);
	 Unext(1,i,j,kic) = a_Up(1,i,j,kic) + irho*(Lutt(1,i,j,kic)+force(1,i,j,kic));
	 Unext(2,i,j,kic) = a_Up(2,i,j,kic) + irho*(Lutt(2,i,j,kic)+force(2,i,j,kic));
	 Unext(3,i,j,kic) = a_Up(3,i,j,kic) + irho*(Lutt(3,i,j,kic)+force(3,i,j,kic));
      }
}

//-----------------------------------------------------------------------
void EW::compute_preliminary_predictor( Sarray& a_Up, Sarray& a_U, Sarray& Unext,
					int g, int kic, double t, vector<GridPointSource*> point_sources )
{
   int ib=m_iStart[g], jb=m_jStart[g], kb=m_kStart[g];
   int ie=m_iEnd[g], je=m_jEnd[g], ke=m_kEnd[g];

   // Compute L(Up) at k=kic.
   Sarray Lu(3,ib,ie,jb,je,kic,kic);
   char op='=';
   int nz = m_global_nz[g];
// Note: 6 first arguments of the function call:
// (ib,ie), (jb,je), (kb,ke) is the declared size of mMu and mLambda in the (i,j,k)-directions, respectively
   rhs4th3fortwind( &ib, &ie, &jb, &je, &kb, &ke, &nz, m_onesided[g], m_acof,
                    m_bope, m_ghcof, Lu.c_ptr(), a_Up.c_ptr(), mMu[g].c_ptr(),
                    mLambda[g].c_ptr(), &mGridSize[g], m_sg_str_x[g], m_sg_str_y[g],
                    m_sg_str_z[g], &op, &kb, &ke, &kic, &kic ); 
// Note: 4 last arguments of the above function call:
// (kb,ke) is the declared size of Up in the k-direction
// (kic,kic) is the declared size of Lu in the k-direction
   
   // Compute forcing at k=kic
   Sarray f(3,ib,ie,jb,je,kic,kic);
   if( m_twilight_forcing )
   {
      double om = m_twilight_forcing->m_omega;
      double ph = m_twilight_forcing->m_phase;
      double cv = m_twilight_forcing->m_c;
      double omm= m_twilight_forcing->m_momega;
      double phm= m_twilight_forcing->m_mphase;
      double amprho=m_twilight_forcing->m_amprho;
      double ampmu=m_twilight_forcing->m_ampmu;
      double amplambda=m_twilight_forcing->m_amplambda;
      if( usingSupergrid() )
      {
         double omstrx = m_supergrid_taper_x[g].get_tw_omega();
         double omstry = m_supergrid_taper_y[g].get_tw_omega();
         double omstrz = m_supergrid_taper_z[g].get_tw_omega();
         forcingfortsg(  &ib, &ie, &jb, &je, &kic, &kic, f.c_ptr(), &t, &om, &cv, &ph, &omm, &phm,
                         &amprho, &ampmu, &amplambda, &mGridSize[g], &m_zmin[g],
                         &omstrx, &omstry, &omstrz );
      }
      else
      {
         forcingfort( &ib, &ie, &jb, &je, &kic, &kic, f.c_ptr(), &t, &om, &cv, &ph, &omm, &phm,
                      &amprho, &ampmu, &amplambda, &mGridSize[g], &m_zmin[g] );
      }
           
   }
   else if( m_rayleigh_wave_test || m_energy_test )
      f.set_to_zero();
   else
   {
     // Default: m_point_source_test, m_lamb_test or full seismic case
      f.set_to_zero();
      for( int s = 0 ; s < point_sources.size() ; s++ )
      {
	 if( point_sources[s]->m_grid == g && point_sources[s]->m_k0 == kic )
	 {
	    double fxyz[3];
	    point_sources[s]->getFxyz(t,fxyz);
	    f(1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
	    f(2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
	    f(3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
	 }
      }
   }
   double cof = mDt*mDt;
// initialize
   Unext.set_to_zero();
   for( int j=jb+2 ; j <= je-2 ; j++ )
      for( int i=ib+2 ; i <= ie-2 ; i++ )
      {
	 double irho=cof/mRho[g](i,j,kic);
	 Unext(1,i,j,kic) = 2*a_Up(1,i,j,kic) - a_U(1,i,j,kic) + irho*(Lu(1,i,j,kic)+f(1,i,j,kic));
	 Unext(2,i,j,kic) = 2*a_Up(2,i,j,kic) - a_U(2,i,j,kic) + irho*(Lu(2,i,j,kic)+f(2,i,j,kic));
	 Unext(3,i,j,kic) = 2*a_Up(3,i,j,kic) - a_U(3,i,j,kic) + irho*(Lu(3,i,j,kic)+f(3,i,j,kic));
      }
// add in super-grid damping terms
   if (usingSupergrid()) // assume 4th order AD, Cartesian grid
   {
// assign array pointers on the fly
      addsg4wind( &mDt, &mGridSize[g], Unext.c_ptr(), a_Up.c_ptr(), a_U.c_ptr(), mRho[g].c_ptr(),
                  m_sg_dc_x[g], m_sg_dc_y[g], m_sg_dc_z[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g],
                  m_sg_corner_x[g], m_sg_corner_y[g], m_sg_corner_z[g],
                  ib, ie, jb, je, kb, ke, m_supergrid_damping_coefficient, Unext.m_kb, Unext.m_ke, kic, kic );
      // Note: the last four arguments define the declared size of Unext, followed by the lower and upper boundaries of the k-window
   }
}

//-----------------------------------------------------------------------
void EW::compute_icstresses( Sarray& a_Up, Sarray& B, int g, int kic,
			     double* a_str_x, double* a_str_y )
{
   const double a1=2.0/3, a2=-1.0/12;
   bool upper = (kic == 1);
   int k=kic;
   double ih = 1/mGridSize[g];
   double uz, vz, wz;
   int ifirst = a_Up.m_ib;
   int jfirst = a_Up.m_jb;
#define str_x(i) a_str_x[(i-ifirst)]   
#define str_y(j) a_str_y[(j-jfirst)]   

   for( int j=B.m_jb+2 ; j <= B.m_je-2 ; j++ )
      for( int i=B.m_ib+2 ; i <= B.m_ie-2 ; i++ )
      {
	 uz = vz = wz = 0;
	 if( upper )
	 {
	    for( int m=0 ; m <= 4 ; m++ )
	    {
	       uz += m_sbop[m]*a_Up(1,i,j,k+m-1);
	       vz += m_sbop[m]*a_Up(2,i,j,k+m-1);
	       wz += m_sbop[m]*a_Up(3,i,j,k+m-1);
	    }
	 }
	 else
	 {
	    for( int m=0 ; m <= 4 ; m++ )
	    {
	       uz -= m_sbop[m]*a_Up(1,i,j,k+1-m);
	       vz -= m_sbop[m]*a_Up(2,i,j,k+1-m);
	       wz -= m_sbop[m]*a_Up(3,i,j,k+1-m);
	    }
	 }
	 B(1,i,j,k) = ih*mMu[g](i,j,k)*(
	    str_x(i)*( a2*(a_Up(3,i+2,j,k)-a_Up(3,i-2,j,k))+
		       a1*(a_Up(3,i+1,j,k)-a_Up(3,i-1,j,k))) + (uz)  );
	 B(2,i,j,k) = ih*mMu[g](i,j,k)*(
	          str_y(j)*( a2*(a_Up(3,i,j+2,k)-a_Up(3,i,j-2,k))+
			     a1*(a_Up(3,i,j+1,k)-a_Up(3,i,j-1,k))) +
	      (vz)  );
	 B(3,i,j,k) = ih*((2*mMu[g](i,j,k)+mLambda[g](i,j,k))*(wz) + mLambda[g](i,j,k)*(
  	   str_x(i)*( a2*(a_Up(1,i+2,j,k)-a_Up(1,i-2,j,k))+a1*(a_Up(1,i+1,j,k)-a_Up(1,i-1,j,k))) +
	   str_y(j)*( a2*(a_Up(2,i,j+2,k)-a_Up(2,i,j-2,k))+a1*(a_Up(2,i,j+1,k)-a_Up(2,i,j-1,k))) ) );
      }
#undef str_x
#undef str_y
}

//-----------------------------------------------------------------------
void EW::consintp( Sarray& Uf, Sarray& Unextf, Sarray& Bf, Sarray& Muf, Sarray& Lambdaf, Sarray& Rhof, double hf,
		   Sarray& Uc, Sarray& Unextc, Sarray& Bc, Sarray& Muc, Sarray& Lambdac, Sarray& Rhoc, double hc,
		   double cof, int gc, int gf, int is_periodic[2])
{
   // At boundaries to the left and right, at least three ghost points are required
   // e.g., domain in i-direction:   i=-2,-1,0,1,2,...,Ni,Ni+1,Ni+2,Ni+3
   // we solve for ghost points at i=2,..,Ni-1, assuming Dirichlet conditions given on i=1,i=Ni and
   // at the ghost points i=-2,-1,0,Ni+1,Ni+2,Ni+3. 
   //
   // The arrays Bf,Unextf, etc are assumed to be computed with all z-ghost points zero,
   // i.e., Uf(i,Nkf+1) was zero for i=-2,-1,..,Ni+3 when Bf and Unextf were evaluated from Uf 
   // (similarly for Unextc and Bc).
   //
   // Before this routine is called, correct boundary values for
   // Uf(-2,Nkf+1),..Uf(1,Nkf+1) (and similarly at the upper i-boundary, and for Uc) must be imposed.
   // In this way the restriction and prolongation stencils can be computed without any special
   // treatment near the (i,j)-boundaries. 

   double *a_strc_x = m_sg_str_x[gc];
   double *a_strc_y = m_sg_str_y[gc];
   double *a_strf_x  = m_sg_str_x[gf];
   double *a_strf_y  = m_sg_str_y[gf];

// stretching on the coarse side
#define strc_x(i) a_strc_x[(i-m_iStart[gc])]   
#define strc_y(j) a_strc_y[(j-m_jStart[gc])]   

// stretching on the fine side
#define strf_x(i) a_strf_x[(i-m_iStart[gf])]   
#define strf_y(j) a_strf_y[(j-m_jStart[gf])]   

      
   const double i16 = 1.0/16;
   const double i256 = 1.0/256;
   const double i1024 = 1.0/1024;

   int jcb, jce, icb, ice, jfb, jfe, ifb, ife, nkf;
   double nuf = mDt*mDt/(cof*hf*hf); // cof=12 for the predictor, cof=1 for the corrector (argument to this routine)
   double nuc = mDt*mDt/(cof*hc*hc);
   double ihc = 1/hc, ihf=1/hf;
   double jacerr = m_citol+1,jacerr0;
   double a11, a12, a21, a22, b1, b2, r1, r2, r3, deti, relax;
   int it = 0;
   relax = m_cirelfact;
 
   icb = m_iStartInt[gc];
   ifb = m_iStartInt[gf];

   ice = m_iEndInt[gc];
   ife = m_iEndInt[gf];
   
   jcb = m_jStartInt[gc];
   jfb = m_jStartInt[gf];

   jce = m_jEndInt[gc];
   jfe = m_jEndInt[gf];

   nkf = m_global_nz[gf];
   Sarray Mlf(m_iStart[gf],m_iEnd[gf],m_jStart[gf],m_jEnd[gf],nkf,nkf);
   for( int j=m_jStart[gf] ; j<=m_jEnd[gf] ; j++ )
      for( int i=m_iStart[gf] ; i<=m_iEnd[gf] ; i++ )
	 Mlf(i,j,nkf) = 2*Muf(i,j,nkf)+Lambdaf(i,j,nkf); // 2*mu + lambda on the fine grid
   Sarray Morc(m_iStart[gc],m_iEnd[gc],m_jStart[gc],m_jEnd[gc],1,1);
   Sarray Mlrc(m_iStart[gc],m_iEnd[gc],m_jStart[gc],m_jEnd[gc],1,1);
   for( int jc=m_jStart[gc] ; jc<=m_jEnd[gc] ; jc++ )
      for( int ic=m_iStart[gc] ; ic<=m_iEnd[gc] ; ic++ )
      {
	 double irho=1/Rhoc(ic,jc,1);
//	 Morc(ic,jc,1) = Muc(ic,jc,1)*irho; // mu/rho on the coarse grid (no stretching)
//	 Mlrc(ic,jc,1) = (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*irho; // (2*mu+lambda)/rho on the coarse grid
// new: include stretching in Morc and Mlrc
	 Morc(ic,jc,1) = Muc(ic,jc,1)*irho/(strc_x(ic)*strc_y(jc)); // mu/rho/(strc_x*strc_y) on the coarse grid
	 Mlrc(ic,jc,1) = (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*irho/(strc_x(ic)*strc_y(jc)); // (2*mu+lambda)/rho/(strc_x*strc_y)
// new: include stretching terms in Unextc
         for (int c=1; c<=3; c++)
            Unextc(c,ic,jc,1) = Unextc(c,ic,jc,1)/(strc_x(ic)*strc_y(jc));
      }

// new: include stretching terms in Unextf
   for( int j=m_jStart[gf] ; j<=m_jEnd[gf] ; j++ )
      for( int i=m_iStart[gf] ; i<=m_iEnd[gf] ; i++ )
         for (int c=1; c<=3; c++)
         {
            Unextf(c,i,j,nkf) = Unextf(c,i,j,nkf)/(strf_x(i)*strf_y(j));
         }
      
// start iteration
   while( jacerr > m_citol && it < m_cimaxiter )
   {
      double rmax[6]={0,0,0,0,0,0};
//
// REMARK: check jump condition in the presence of stretching function;
// stretching function may be different in the fine and coarse grids!
//
// for i=2*ic-1 and j=2*jc-1: Enforce continuity of displacements and normal stresses along the interface
      for( int jc= jcb ; jc <= jce ; jc++ )
	 for( int ic= icb ; ic <= ice ; ic++ )
	 {
// i odd, j odd
	    int i=2*ic-1, j=2*jc-1;
            // setup 2x2 system matrix
            // unknowns: (Uf, Uc)
// eqn 1: continuity of normal stress
	    a11 = 0.25*Muf(i,j,nkf)*m_sbop[0]*ihf; // ihf = 1/h on the fine grid
	    a12 =      Muc(ic,jc,1)*m_sbop[0]*ihc;  // ihc = 1/h on the coarse grid
// eqn 2: continuity of displacement
// nuf = dt^2/(cof * hf^2), nuc = dt^2/(cof * hc^2)
            // NEED stretching terms in the matrix elements!!!
	    // a21 = nuf/Rhof(i,j,nkf)*Muf(i,j,nkf)*m_ghcof[0]; 
	    // a22 =-nuc/Rhoc(ic,jc,1)*Muc(ic,jc,1)*m_ghcof[0];
	    a21 = nuf/Rhof(i,j,nkf)*Muf(i,j,nkf)*m_ghcof[0]/(strf_x(i)*strf_y(j)); 
	    a22 =-nuc/Rhoc(ic,jc,1)*Muc(ic,jc,1)*m_ghcof[0]/(strc_x(ic)*strc_y(jc)); 
	    for( int c=1 ; c <= 2 ; c++ ) //  the 2 tangential components ? 
	    {
// apply the restriction operator to the normal stress on the interface (Bf is on the fine grid)
	       b1  = i1024*( 
                  Bf(c,i-3,j-3,nkf)-9*Bf(c,i-3,j-1,nkf)-16*Bf(c,i-3,j,nkf)-9*Bf(c,i-3,j+1,nkf)+Bf(c,i-3,j+3,nkf) +
                  9*(-Bf(c,i-1,j-3,nkf)+9*Bf(c,i-1,j-1,nkf)+16*Bf(c,i-1,j,nkf)+9*Bf(c,i-1,j+1,nkf)-Bf(c,i-1,j+3,nkf)) +
                  16*(-Bf(c,i,  j-3,nkf)+9*Bf(c,i,  j-1,nkf)+16*Bf(c,i,  j,nkf)+9*Bf(c,i,  j+1,nkf)-Bf(c,i,  j+3,nkf)) + 
                  9*(-Bf(c,i+1,j-3,nkf)+9*Bf(c,i+1,j-1,nkf)+16*Bf(c,i+1,j,nkf)+9*Bf(c,i+1,j+1,nkf)-Bf(c,i+1,j+3,nkf)) +
                  Bf(c,i+3,j-3,nkf)-9*Bf(c,i+3,j-1,nkf)-16*Bf(c,i+3,j,nkf)-9*Bf(c,i+3,j+1,nkf)+Bf(c,i+3,j+3,nkf) );

	       b1 = b1 - i1024*m_sbop[0]*ihf*(
                  Muf(i-3,j-3,nkf)*Uf(c,i-3,j-3,nkf+1) - 9*Muf(i-3,j-1,nkf)*Uf(c,i-3,j-1,nkf+1)
                  -16*Muf(i-3,j,  nkf)*Uf(c,i-3,  j,nkf+1) - 9*Muf(i-3,j+1,nkf)*Uf(c,i-3,j+1,nkf+1)
                  +Muf(i-3,j+3,nkf)*Uf(c,i-3,j+3,nkf+1) +					    
                  9*(  -Muf(i-1,j-3,nkf)*Uf(c,i-1,j-3,nkf+1) + 9*Muf(i-1,j-1,nkf)*Uf(c,i-1,j-1,nkf+1) 
                       +16*Muf(i-1,j,  nkf)*Uf(c,i-1,j,  nkf+1) + 9*Muf(i-1,j+1,nkf)*Uf(c,i-1,j+1,nkf+1) 
                       -Muf(i-1,j+3,nkf)*Uf(c,i-1,j+3,nkf+1) ) +
                  16*(  -Muf(i,  j-3,nkf)*Uf(c,i,  j-3,nkf+1) + 9*Muf(i,  j-1,nkf)*Uf(c,i,  j-1,nkf+1) 
                        + 9*Muf(i,  j+1,nkf)*Uf(c,i,  j+1,nkf+1)
                        -Muf(i,  j+3,nkf)*Uf(c,i,  j+3,nkf+1) ) + 
                  9*(  -Muf(i+1,j-3,nkf)*Uf(c,i+1,j-3,nkf+1) + 9*Muf(i+1,j-1,nkf)*Uf(c,i+1,j-1,nkf+1)
                       +16*Muf(i+1,j,  nkf)*Uf(c,i+1,j,  nkf+1) + 9*Muf(i+1,j+1,nkf)*Uf(c,i+1,j+1,nkf+1)
                       -Muf(i+1,j+3,nkf)*Uf(c,i+1,j+3,nkf+1) ) +
                  Muf(i+3,j-3,nkf)*Uf(c,i+3,j-3,nkf+1) - 9*Muf(i+3,j-1,nkf)*Uf(c,i+3,j-1,nkf+1)
                  -16*Muf(i+3,j,  nkf)*Uf(c,i+3,j,  nkf+1) - 9*Muf(i+3,j+1,nkf)*Uf(c,i+3,j+1,nkf+1)
                  +Muf(i+3,j+3,nkf)*Uf(c,i+3,j+3,nkf+1) );

	       b1 = b1 - Bc(c,ic,jc,1);
            // NEED stretching terms in b2!!!
	       b2 = Unextc(c,ic,jc,1)-Unextf(c,i,j,nkf); // add stretching functions
//	       b2 = Unextc(c,ic,jc,1)/(strc_x(ic)*strc_y(jc)) - Unextf(c,i,j,nkf)/(strf_x(i)*strf_y(j)); // stretching functions added
	       deti=1/(a11*a22-a12*a21);
	       r1 = Uf(c,i,j,nkf+1);
	       r2 = Uc(c,ic,jc,0);
// solve the linear 2x2 system
	       Uf(c,i,j,nkf+1) = deti*( a22*b1-a12*b2);
	       Uc(c,ic,jc,0)   = deti*(-a21*b1+a11*b2);
// damp the update of the ghost point values (r1, r2) hold previous values
	       Uf(c,i,j,nkf+1) = relax*Uf(c,i,j,nkf+1) + (1-relax)*r1;
	       Uc(c,ic,jc,0)   = relax*Uc(c,ic,jc,0)   + (1-relax)*r2;
// change in solution
	       r1 = r1-Uf(c,i,j,nkf+1);
	       r2 = r2-Uc(c,ic,jc,0);
	       rmax[c-1] = rmax[c-1] > fabs(r1) ? rmax[c-1] : fabs(r1);
	       rmax[c-1] = rmax[c-1] > fabs(r2) ? rmax[c-1] : fabs(r2);
	       //	       if( c == 2 && ic == 12 && jc == 13 )
	       //	       {
	       //	         cout << "i,j " << i << " " << j << " " << b1 << " " << b2 << " " << r1 << " " << r2 << endl;
	       //		 cout << "   " << Uf(c,i,j,nkf+1) << " " << Uc(c,ic,jc,0) << " " << 
	       //		    a21*Uf(c,i,j,nkf+1)+Unextf(c,i,j,nkf) << " " << -a22*Uc(c,ic,jc,0)+Unextc(c,ic,jc,1) << endl;
	       //	       }
	    } // end for c=1,2
            
// setup the matrix for the 3rd component of the normal stress (different coefficients)
	    a11 = 0.25*Mlf(i,j,nkf)*m_sbop[0]*ihf;
	    a12 =    (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*m_sbop[0]*ihc;

            // NEED stretching terms in the matrix elements a21 & a22!!!
	    // a21 = nuf/Rhof(i,j,nkf)*Mlf(i,j,nkf)*m_ghcof[0];
	    // a22 =-nuc/Rhoc(ic,jc,1)*(2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*m_ghcof[0];
	    a21 = nuf/Rhof(i,j,nkf)*Mlf(i,j,nkf)*m_ghcof[0]/(strf_x(i)*strf_y(j)); 
	    a22 =-nuc/Rhoc(ic,jc,1)*(2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*m_ghcof[0]/(strc_x(ic)*strc_y(jc)); 
// apply the restriction operator to the fine grid normal stress grid function (Bf)
	    b1  = i1024*( 
                    Bf(3,i-3,j-3,nkf)-9*Bf(3,i-3,j-1,nkf)-16*Bf(3,i-3,j,nkf)-9*Bf(3,i-3,j+1,nkf)+Bf(3,i-3,j+3,nkf) +
	        9*(-Bf(3,i-1,j-3,nkf)+9*Bf(3,i-1,j-1,nkf)+16*Bf(3,i-1,j,nkf)+9*Bf(3,i-1,j+1,nkf)-Bf(3,i-1,j+3,nkf)) +
	       16*(-Bf(3,i,  j-3,nkf)+9*Bf(3,i,  j-1,nkf)+16*Bf(3,i,  j,nkf)+9*Bf(3,i,  j+1,nkf)-Bf(3,i,  j+3,nkf)) + 
	        9*(-Bf(3,i+1,j-3,nkf)+9*Bf(3,i+1,j-1,nkf)+16*Bf(3,i+1,j,nkf)+9*Bf(3,i+1,j+1,nkf)-Bf(3,i+1,j+3,nkf)) +
	  	    Bf(3,i+3,j-3,nkf)-9*Bf(3,i+3,j-1,nkf)-16*Bf(3,i+3,j,nkf)-9*Bf(3,i+3,j+1,nkf)+Bf(3,i+3,j+3,nkf) );
	    b1 = b1 - i1024*m_sbop[0]*ihf*(
                Mlf(i-3,j-3,nkf)*Uf(3,i-3,j-3,nkf+1) - 9*Mlf(i-3,j-1,nkf)*Uf(3,i-3,j-1,nkf+1)
   	    -16*Mlf(i-3,j,  nkf)*Uf(3,i-3,  j,nkf+1) - 9*Mlf(i-3,j+1,nkf)*Uf(3,i-3,j+1,nkf+1)
               +Mlf(i-3,j+3,nkf)*Uf(3,i-3,j+3,nkf+1) +					    
	  9*(  -Mlf(i-1,j-3,nkf)*Uf(3,i-1,j-3,nkf+1) + 9*Mlf(i-1,j-1,nkf)*Uf(3,i-1,j-1,nkf+1) 
	    +16*Mlf(i-1,j,  nkf)*Uf(3,i-1,j,  nkf+1) + 9*Mlf(i-1,j+1,nkf)*Uf(3,i-1,j+1,nkf+1) 
	       -Mlf(i-1,j+3,nkf)*Uf(3,i-1,j+3,nkf+1) ) +
	 16*(  -Mlf(i,  j-3,nkf)*Uf(3,i,  j-3,nkf+1) + 9*Mlf(i,  j-1,nkf)*Uf(3,i,  j-1,nkf+1) 
	                                             + 9*Mlf(i,  j+1,nkf)*Uf(3,i,  j+1,nkf+1)
               -Mlf(i,  j+3,nkf)*Uf(3,i,  j+3,nkf+1) ) + 
	  9*(  -Mlf(i+1,j-3,nkf)*Uf(3,i+1,j-3,nkf+1) + 9*Mlf(i+1,j-1,nkf)*Uf(3,i+1,j-1,nkf+1)
            +16*Mlf(i+1,j,  nkf)*Uf(3,i+1,j,  nkf+1) + 9*Mlf(i+1,j+1,nkf)*Uf(3,i+1,j+1,nkf+1)
	       -Mlf(i+1,j+3,nkf)*Uf(3,i+1,j+3,nkf+1) ) +
	        Mlf(i+3,j-3,nkf)*Uf(3,i+3,j-3,nkf+1) - 9*Mlf(i+3,j-1,nkf)*Uf(3,i+3,j-1,nkf+1)
	    -16*Mlf(i+3,j,  nkf)*Uf(3,i+3,j,  nkf+1) - 9*Mlf(i+3,j+1,nkf)*Uf(3,i+3,j+1,nkf+1)
	       +Mlf(i+3,j+3,nkf)*Uf(3,i+3,j+3,nkf+1) );

// setup the RHS
	       b1 = b1 - Bc(3,ic,jc,1);
	       b2 = Unextc(3,ic,jc,1)-Unextf(3,i,j,nkf); // need stretching terms in b2!!!
//	       b2 = Unextc(3,ic,jc,1)/(strc_x(ic)*strc_y(jc)) - Unextf(3,i,j,nkf)/(strf_x(i)*strf_y(j)); // stretching functions added
	       deti=1/(a11*a22-a12*a21);
// previous values
	       r1 = Uf(3,i,j,nkf+1);
	       r2 = Uc(3,ic,jc,0);
// solve the 2x2 system for component 3 of Uf and Uc
	       Uf(3,i,j,nkf+1) = deti*( a22*b1-a12*b2);
	       Uc(3,ic,jc,0)   = deti*(-a21*b1+a11*b2);
// relax the updated value
	       Uf(3,i,j,nkf+1) = relax*Uf(3,i,j,nkf+1) + (1-relax)*r1;
	       Uc(3,ic,jc,0)   = relax*Uc(3,ic,jc,0)   + (1-relax)*r2;
// change in ghost point values
	       r1 = r1-Uf(3,i,j,nkf+1);
	       r2 = r2-Uc(3,ic,jc,0);
	       rmax[2] = rmax[2] > fabs(r1) ? rmax[2] : fabs(r1);
	       rmax[2] = rmax[2] > fabs(r2) ? rmax[2] : fabs(r2);
	       //	       int c=3;
	       //	       if( c == 3 && ic == 12 && jc == 13 )
	       //	       {
	       //	         cout << "i,j " << i << " " << j << " " << b1 << " " << b2 << " " << r1 << " " << r2 << endl;
	       //		 cout << "   " << Uf(c,i,j,nkf+1) << " " << Uc(c,ic,jc,0) << " " << 
	       //		    a21*Uf(c,i,j,nkf+1)+Unextf(c,i,j,nkf) << " " << -a22*Uc(c,ic,jc,0)+Unextc(c,ic,jc,1) << endl;
	       //	       }
	 }
//      
// Enforce continuity of displacements along the interface (for fine ghost points in between coarse points)
//
// TODO: insert coarse and fine stretching functions below
//
      int ic, jc;
      for( int j=jfb ; j <= jfe ; j++ )
	 for( int i=ifb ; i <= ife ; i++ )
	 {
	    if( !( (i % 2 == 1 && j % 2 == 1 ) ) ) // not both i and j are odd (handled above)
	    {
// updated components 1,2 of the ghost point value of Uf
	       for( int c=1 ; c <= 2 ; c++ )
	       {
		  if( (j % 2 == 0) && (i % 2 == 1) ) // j is even, i is odd
		  {
		     ic = (i+1)/2;
		     jc = j/2;
// All Unextc/str_coarse
		     b1 = i16*(-Unextc(c,ic,jc-1,1)+9*(Unextc(c,ic,jc,1)+Unextc(c,ic,jc+1,1))-Unextc(c,ic,jc+2,1));
// All Uc/str_coarse
		     b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(c,ic,jc-1,0)*Morc(ic,jc-1,1) + 
						      9*Uc(c,ic,jc  ,0)*Morc(ic,jc  ,1) + 
						      9*Uc(c,ic,jc+1,0)*Morc(ic,jc+1,1)
						      -Uc(c,ic,jc+2,0)*Morc(ic,jc+2,1) );
		  }
		  if( (j % 2 == 1) && (i % 2 == 0) ) // j is odd, i is even
		  {
		     ic = i/2;
		     jc = (j+1)/2;
// All Unextc/str_coarse
		     b1 = i16*(-Unextc(c,ic-1,jc,1)+9*(Unextc(c,ic,jc,1)+Unextc(c,ic+1,jc,1))-Unextc(c,ic+2,jc,1));
// All Uc/str_coarse
		     b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)+ 
						      9*Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+ 
						      9*Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1)
						      -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1));
		  }
		  if( (j % 2 == 0) && (i % 2 == 0) ) // i is even, j is even
		  {
		     ic = i/2;
		     jc = j/2;
// All Unextc/str_coarse
		     b1 = i256*
                            ( Unextc(c,ic-1,jc-1,1)-9*(Unextc(c,ic,jc-1,1)+Unextc(c,ic+1,jc-1,1))+Unextc(c,ic+2,jc-1,1)
	 	        + 9*(-Unextc(c,ic-1,jc,  1)+9*(Unextc(c,ic,jc,  1)+Unextc(c,ic+1,jc,  1))-Unextc(c,ic+2,jc,  1)  
			     -Unextc(c,ic-1,jc+1,1)+9*(Unextc(c,ic,jc+1,1)+Unextc(c,ic+1,jc+1,1))-Unextc(c,ic+2,jc+1,1))
		             +Unextc(c,ic-1,jc+2,1)-9*(Unextc(c,ic,jc+2,1)+Unextc(c,ic+1,jc+2,1))+Unextc(c,ic+2,jc+2,1) );

// All Uc/str_coarse
		     b1 = b1 + nuc*m_ghcof[0]*i256*(
                        Uc(c,ic-1,jc-1,0)*Morc(ic-1,jc-1,1)-9*(Uc(c,ic,  jc-1,0)*Morc(ic,  jc-1,1)+Uc(c,ic+1,jc-1,0)*Morc(ic+1,jc-1,1)) +
                        Uc(c,ic+2,jc-1,0)*Morc(ic+2,jc-1,1)
                        + 9*(
                           -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)
                           +9*(Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1))
                           -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1)  
                           -Uc(c,ic-1,jc+1,0)*Morc(ic-1,jc+1,1)+
                           9*(Uc(c,ic,  jc+1,0)*Morc(ic,  jc+1,1)+Uc(c,ic+1,jc+1,0)*Morc(ic+1,jc+1,1))
                           -Uc(c,ic+2,jc+1,0)*Morc(ic+2,jc+1,1)
                           )
                        + Uc(c,ic-1,jc+2,0)*Morc(ic-1,jc+2,1)
                        -9*(Uc(c,ic,  jc+2,0)*Morc(ic,  jc+2,1) +Uc(c,ic+1,jc+2,0)*Morc(ic+1,jc+2,1))
                        +Uc(c,ic+2,jc+2,0)*Morc(ic+2,jc+2,1)
                        );
		  }
		  b1 = b1 - Unextf(c,i,j,nkf); // b1*str_fine on rhs
//		  b1 = b1 - Unextf(c,i,j,nkf)/(strf_x(i)*strf_y(j)); // b1*str_fine on rhs
//                  a11 = nuf*m_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
                  a11 = nuf*m_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf))/(strf_x(i)*strf_y(j));
		  r3 = Uf(c,i,j,nkf+1); // save old value for relaxation
// update ghost point value Uf(c,i,j,nkf+1)
//		  Uf(c,i,j,nkf+1) = b1*Rhof(i,j,nkf)/(nuf*m_ghcof[0]*Muf(i,j,nkf)); // without stretching
                  Uf(c,i,j,nkf+1) = b1/a11; // with stretching
		  Uf(c,i,j,nkf+1) = relax*Uf(c,i,j,nkf+1)+(1-relax)*r3;
		  //		  if( i == 4 && j == 7 && c == 1)
		  //		     cout << "in loop " << -a11*Uf(c,i,j,nkf+1) + b1  << endl;
// change in ghost point value
		  r3 = r3 - Uf(c,i,j,nkf+1);
		  rmax[c-1+3] = rmax[c-1+3] > fabs(r3) ? rmax[c-1+3] : fabs(r3);
	       } // end for c=1,2
               
// work on componet 3 of the ghost point value of Uf
	       if( (j % 2 == 0) && (i % 2 == 1) ) // j even, i odd
	       {
		  ic = (i+1)/2;
		  jc = j/2;
// All Unextc/str_coarse
		  b1 = i16*(-Unextc(3,ic,jc-1,1)+9*(Unextc(3,ic,jc,1)+Unextc(3,ic,jc+1,1))-Unextc(3,ic,jc+2,1));
// All Uc/str_coarse
// stretching coeff in Mlrc
		  b1 = b1 + nuc*m_ghcof[0]*i16*(   - Uc(3,ic,jc-1,0)*Mlrc(ic,jc-1,1) + 
						   9*Uc(3,ic,jc  ,0)*Mlrc(ic,jc  ,1) + 
						   9*Uc(3,ic,jc+1,0)*Mlrc(ic,jc+1,1)
						    -Uc(3,ic,jc+2,0)*Mlrc(ic,jc+2,1) );
	       }
	       if( (j % 2 == 1) && (i % 2 == 0) ) // j odd, i even
	       {
		  ic = i/2;
		  jc = (j+1)/2;
// All Unextc/str_coarse
		  b1 = i16*(-Unextc(3,ic-1,jc,1)+9*(Unextc(3,ic,jc,1)+Unextc(3,ic+1,jc,1))-Unextc(3,ic+2,jc,1));
// All Uc/str_coarse
		  b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(3,ic-1,jc,0)*Mlrc(ic-1,jc,1)+ 
						  9*Uc(3,ic,  jc,0)*Mlrc(ic,  jc,1)+ 
						  9*Uc(3,ic+1,jc,0)*Mlrc(ic+1,jc,1)
						   -Uc(3,ic+2,jc,0)*Mlrc(ic+2,jc,1));
	       }
	       if( (j % 2 == 0) && (i % 2 == 0) ) // j even, i even
	       {
		  ic = i/2;
		  jc = j/2;
// All Unextc/str_coarse
// stretching coeff in Unextc
		  b1 = i256*
                     ( Unextc(3,ic-1,jc-1,1)-9*(Unextc(3,ic,jc-1,1)+Unextc(3,ic+1,jc-1,1))+Unextc(3,ic+2,jc-1,1)
                       + 9*(-Unextc(3,ic-1,jc,  1)+9*(Unextc(3,ic,jc,  1)+Unextc(3,ic+1,jc,  1))-Unextc(3,ic+2,jc,  1)  
                            -Unextc(3,ic-1,jc+1,1)+9*(Unextc(3,ic,jc+1,1)+Unextc(3,ic+1,jc+1,1))-Unextc(3,ic+2,jc+1,1))
                       +Unextc(3,ic-1,jc+2,1)-9*(Unextc(3,ic,jc+2,1)+Unextc(3,ic+1,jc+2,1))+Unextc(3,ic+2,jc+2,1) );

// All Uc/str_coarse
// stretching coeff in Mlrc
		  b1 = b1 + nuc*m_ghcof[0]*i256*(
                     Uc(3,ic-1,jc-1,0)*Mlrc(ic-1,jc-1,1)
                     -9*(Uc(3,ic,  jc-1,0)*Mlrc(ic,  jc-1,1)+
                         Uc(3,ic+1,jc-1,0)*Mlrc(ic+1,jc-1,1))+
                     Uc(3,ic+2,jc-1,0)*Mlrc(ic+2,jc-1,1)
                     + 9*(
                        -Uc(3,ic-1,jc,0)*Mlrc(ic-1,jc,1)+
                        9*( Uc(3,ic,  jc,0)*Mlrc(ic,  jc,1)+
                            Uc(3,ic+1,jc,0)*Mlrc(ic+1,jc,1))
                        -Uc(3,ic+2,jc,0)*Mlrc(ic+2,jc,1)  
                        -Uc(3,ic-1,jc+1,0)*Mlrc(ic-1,jc+1,1)+
                        9*( Uc(3,ic,  jc+1,0)*Mlrc(ic,  jc+1,1)+
                            Uc(3,ic+1,jc+1,0)*Mlrc(ic+1,jc+1,1))
                        -Uc(3,ic+2,jc+1,0)*Mlrc(ic+2,jc+1,1)  )
                     + Uc(3,ic-1,jc+2,0)*Mlrc(ic-1,jc+2,1)
                     -9*(Uc(3,ic,  jc+2,0)*Mlrc(ic,  jc+2,1) +
                         Uc(3,ic+1,jc+2,0)*Mlrc(ic+1,jc+2,1)) +
                     Uc(3,ic+2,jc+2,0)*Mlrc(ic+2,jc+2,1)   );
	       } // end  j even, i even
// right hand side is mismatch in displacement                
	       b1 = b1 - Unextf(3,i,j,nkf); // b1*str_fine on rhs
		  //	       a11 = nuf*m_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
               a11 = nuf*m_ghcof[0]*Muf(i,j,nkf)/Rhof(i,j,nkf)/(strf_x(i)*strf_y(j));
	       r3 = Uf(3,i,j,nkf+1); // save previous value for relaxation below
// solve for the ghost point value Uf(3,i,j,nkf+1)
//	       Uf(3,i,j,nkf+1) = b1*Rhof(i,j,nkf)/(nuf*m_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))); // no stretching
	       Uf(3,i,j,nkf+1) = b1/a11;
	       Uf(3,i,j,nkf+1) = relax*Uf(3,i,j,nkf+1)+(1-relax)*r3;
	       r3 = r3 - Uf(3,i,j,nkf+1);
	       rmax[2+3] = rmax[2+3] > fabs(r3) ? rmax[2+3] : fabs(r3);

	    } // end if not ( i%2=1 &&  j%2=1), i.e, not both i and j are odd

            // (i,j) both odd is handled by the first iteration
            
	 } // end for all fine grid points on the interface
      
      //   skipthis:
      communicate_array_2d( Uf, gf, nkf+1 );
      communicate_array_2d( Uc, gc, 0 );
      double jacerrtmp = 0;
      for (int q=0; q<6; q++)
         jacerrtmp += rmax[q];
      
      MPI_Allreduce( &jacerrtmp, &jacerr, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
      if( it == 0 )
	 jacerr0 = jacerr;
      if( jacerr0 > 0 )
	 jacerr = jacerr/jacerr0;
      it++;

   } // end while jacerr > eps
   
   if( jacerr > m_citol && proc_zero() )
      cout << "EW::consintp, Warning, no convergence. err = " << jacerr << " tol= " << m_citol << endl;
      
   if( proc_zero() && mVerbose >= 4 )
      cout << "EW::consintp, no of iterations= " << it << " Jac iteration error= " << jacerr << endl;
#undef strc_x
#undef strc_y
#undef strf_x
#undef strf_y
}


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
  double om=0, ph=0, cv=0, omm;
    
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
       double phc[21]; // move these angles to the EW class
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
         if( m_anisotropic )
         {
// curvilinear anisotropic case is not yet implemented
            CHECK_INPUT (!curvilinear, "cartesian_bc_forcing> bStressFree not implemented for anisotropic materials and curvilinear grids" <<endl);

            tw_aniso_free_surf_z( ifirst, ilast, jfirst, jlast, kfirst, klast, k, t, om, cv, ph, omm, phc, bforce_side4_ptr, h, m_zmin[g] );            
         }
         else
         { //isotropic stuff
            
            if( usingSupergrid() && !curvilinear )
            {
               double omstrx = m_supergrid_taper_x[g].get_tw_omega();
               double omstry = m_supergrid_taper_y[g].get_tw_omega();
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
               double omstrx = m_supergrid_taper_x[g].get_tw_omega();
               double omstry = m_supergrid_taper_y[g].get_tw_omega();

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
            } // end supergrid && curvilinear
         
         } // end isotropic case
         
      } // end side==4 is bStressFree

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
         if( m_anisotropic )
         {
// curvilinear anisotropic case is not yet implemented
            CHECK_INPUT (!curvilinear, "cartesian_bc_forcing> bStressFree not implemented for anisotropic materials and curvilinear grids" <<endl);

            tw_aniso_free_surf_z( ifirst, ilast, jfirst, jlast, kfirst, klast, k, t, om, cv, ph, omm, phc, bforce_side5_ptr, h, m_zmin[g] );            
         }
         else
         { //isotropic stuff

            if( usingSupergrid() )
            {
               double omstrx = m_supergrid_taper_x[g].get_tw_omega();
               double omstry = m_supergrid_taper_y[g].get_tw_omega();
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
            } // end ! supergrid
            
         } // end isotropic case
         
      } // end bStressFree on side 5
      
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
  else if(mode == TimeSeries::Strains )
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
  else if(mode == TimeSeries::DisplacementGradient )
  {
     uRec.resize(9);
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
    if( m_sg_damping_order == 4 )
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
    else if(  m_sg_damping_order == 6 )
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
