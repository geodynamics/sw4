#include "mpi.h"

#include "EW.h"

#include <cstring>

#include <iostream>
#include <sstream>
#include <unistd.h>
#include <algorithm>

#include "startEnd.h"
#include "version.h"
#include "F77_FUNC.h"

extern "C" {
void F77_FUNC(corrfort,CORRFORT)(int*, int*, int*, int*, int*, int*, 
				 double*, double*, double*, double*, double* );
void F77_FUNC(dpdmtfort,DPDMTFORT)(int*, int* , int*, int*, int*, int*, 
				  double*, double*, double*, double*, double* );    
void F77_FUNC(predfort,PREDFORT)(int*, int*, int*, int*, int*, int*, 
				 double*, double*, double*, double*, double*, double*, double*);    
void F77_FUNC(rhouttlumf, RHOUTTLUMF)(int*, int*, int*, int*, int*, int*, 
				      int*, double*, double*, double*, double*, double*, double*, double*);
void F77_FUNC(forcingfort,FORCINGFORT)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
				       double*, double*, double*, double*, double* );
void F77_FUNC(forcingfortc,FORCINGFORTC)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
					 double*, double*, double*, double*, double*, double* );
void F77_FUNC(forcingfortsg,FORCINGFORTSG)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
   	 			       double*, double*, double*, double*, double*,double*,double*,double* );
void F77_FUNC(forcingfortcsg,FORCINGFORTCSG)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
			     double*, double*, double*, double*, double*,double*,double*,double*, double* );
void F77_FUNC(forcingttfortsg,FORCINGTTFORTSG)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
					       double*, double*, double*, double*, double*, double*, double*, double* );
void F77_FUNC(forcingttfort,FORCINGTTFORT)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
				       double*, double*, double*, double*, double* );
void F77_FUNC(forcingttfortc,FORCINGTTFORTC)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
					     double*, double*, double*, double*, double*, double* );
void F77_FUNC(forcingttfortcsg,FORCINGTTFORTCSG)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
			     double*, double*, double*, double*, double*,double*,double*,double*, double* );
void F77_FUNC(exactaccfort,EXACTACCFORT)(int*, int*, int*, int*, int*, int*, double*, double*, double*, 
					 double*, double*, double*, double* );
void F77_FUNC(exactaccfortc,EXACTACCFORTC)(int*, int*, int*, int*, int*, int*, double*, double*, double*, 
					   double*, double*, double*, double*, double* );
void F77_FUNC(rhserrfort, RHSERRFORT)(int*, int*, int*, int*, int*, int*, int*, double*,
				      double*, double*, double*, double*, double*);
void F77_FUNC(rhs4th3fort,RHS4TH3FORT)(int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*,
				       double*, double*, double*, double*, double* );
void F77_FUNC(rhs4th3fortsgstr,RHS4TH3FORTSGSTR)(int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*,
						 double*, double*, double*, double*, double*, double*,double*,double* );
void F77_FUNC(exactrhsfort,EXACTRHSFORT)( int*, int*, int*, int*, int*, int*, double*, double*, 
					  double*, double*, double*, double*, double*, double*, double*, double*,
					  double*, double* );
void F77_FUNC(exactrhsfortc,EXACTRHSFORTC)( int*, int*, int*, int*, int*, int*, double*, double*, 
					  double*, double*, double*, double*, double*, double*, double*, double*,
					    double*, double*, double* );
void F77_FUNC(exactrhsfortsg,EXACTRHSFORTSG)( int*, int*, int*, int*, int*, int*, double*, double*,
					      double*, double*, double*, double*, double*,
					      double*, double*, double*, double*, double*,
					      double*, double*, double* );
void F77_FUNC(exactrhsfortsgc,EXACTRHSFORTSGC)( int*, int*, int*, int*, int*, int*, double*, double*,
						double*, double*, double*, double*, double*,
						double*, double*, double*, double*, double*,
						double*, double*, double*, double* );
void F77_FUNC(solerr3, SOLERR3)(int*, int*, int*, int*, int*, int*, double *h, double *uex, double *u, double *li,
				double *l2, double *xli, double *zmin, double *x0, double *y0, double *z0, double *radius,
				int *imin, int *imax, int *jmin, int *jmax, int *kmin, int *kmax);
   void F77_FUNC(solerr3c, SOLERR3c)(int*, int*, int*, int*, int*, int*, double *uex, double *u, double* x, double* y,
                                     double* z, double* jac, double *li, double *l2, double *xli, 
				     double *x0, double *y0, double *z0, double *radius,
				     int *imin, int *imax, int *jmin, int *jmax, int *kmin, int *kmax);
void F77_FUNC(solerrgp, SOLERRGP)(int*, int*, int*, int*, int*, int*, double*, double*, double*, double *li,
				double *l2 );
void F77_FUNC(twilightfort,TWILIGHTFORT)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, 
					  double*, double*, double* );
void F77_FUNC(twilightfortc,TWILIGHTFORTC)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, 
					    double*, double*, double*, double* );
//  subroutine rayleighfort( ifirst, ilast, jfirst, jlast, kfirst, klast,
// +     u, t, lambda, mu, rho, cr, omega, alpha, h, zmin )
void F77_FUNC(rayleighfort,RAYLEIGHFORT)( int*ifirst, int*ilast, int*jfirst, int*jlast, int*kfirst, int*klast, 
					  double*u, double*t, double*lambda, double*mu, 
					  double*rho, double*cr, double*omega, double *alpha, double *h, double *zmin);
void F77_FUNC(velsum,VELSUM)( int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
			      double*, double*, double*, double*, double*, double* );
void F77_FUNC(energy4,ENERGY4)( int*, int*, int*, int*, int*, int*,  int*, int*, int*, int*, int*, int*, int*,
                                   double*, double*, double*, double*, double*, double* );
void F77_FUNC(energy4c,ENERGY4C)( int*, int*, int*, int*, int*, int*,  int*, int*, int*, int*, int*, int*, int*,
                                   double*, double*, double*, double*, double*, double* );
void F77_FUNC(lambexact,LAMBEXACT)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*,
				    double*, double*, double*, int* );
void F77_FUNC(curvilinear4,CURVILINEAR4)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*,
					double*, int*, double*, double*, double* );
void F77_FUNC(curvilinear4sg,CURVILINEAR4SG)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*,
					      double*, double*, int*, double*, double*, double*, double*, double* );

void F77_FUNC(addgradrho,ADDGRADRHO)( int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
				       double*, double*, double*, double*, double*, double*, double*,
				       double*, double*, int* );
void F77_FUNC(addgradrhoc,ADDGRADRHOC)( int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
				       double*, double*, double*, double*, double*, double*, double*,
				       double*, double*, int* );
void F77_FUNC(addgradmula,ADDGRADMULA)( int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
				       double*, double*, double*, double*, double*, double*, double*,
					double*, double*, double*, int*, int*, int*, double* );
}

using namespace std;

#define SQR(x) ((x)*(x))

// constructor
EW::EW(const string& fileName, vector<Source*> & a_GlobalSources,
       vector<TimeSeries*> & a_GlobalTimeSeries, bool a_invproblem ): 
  m_topo_zmax(0.0),
  m_topoInputStyle(UNDEFINED), 
  mTopoImageFound(false),
  m_nx_base(0), m_ny_base(0), m_nz_base(0), m_h_base(0.0),
  mSourcesOK(false),
  mIsInitialized(false),
  mParsingSuccessful(false),
  mNumberOfGrids(0),
  mName(fileName),
  m_scenario(" "),
  mPath("./"),
  mObsPath("./"),
  mWriteGMTOutput(false),
  mPlotFrequency(80),
  mNumFiles(0),
  mVerbose(0),
  mQuiet(false),
  mDebugIO(false),
  mHomogeneous(false),
  m_iotiming(false),
  m_pfs(false),
  m_nwriters(8),
  mTimeIsSet(false),
  mTmax(0.0),
  mTstart(0.0),
  mDt(0.0),
  mNumberOfTimeSteps(-1),
  m_testing(false),
  m_moment_test(false),
  m_twilight_forcing(NULL),
  m_point_source_test(0),
  m_energy_test(0),
  m_lamb_test(0),
  m_rayleigh_wave_test(0),
  m_update_boundary_function(0),
  m_EFileResolution(-1.0),
  m_maxIter(10),
  m_topoFileName("NONE"),
  m_topoExtFileName("NONE"),
  m_QueryType("MAXRES"),
  //  mTestingEnergy(false),
  //  mTestSource(false),
  //  mTestLamb(false),
  mOrder(4),
  mCFL(1.15), // 1.15 is necessary for the rayleigh wave test when Cp/Cs=10
  //  mCFL(1.3),
  // m_d4coeff(0.0),
  // m_d4_cfl(0.2),
  // m_curlcoeff(0.0),

  // mRestartFilePrefix(""),
  // mRestartFromCycle(0),
  // mRestartDumpInterval(0),

  m_doubly_periodic(false),
  mbcsSet(false),

  m_analytical_topo(false),
  m_use_analytical_metric(false),
  m_GaussianAmp(0.05),
  m_GaussianLx(0.15),
  m_GaussianLy(0.15),
  m_GaussianXc(0.5),
  m_GaussianYc(0.5),

  m_use_supergrid(false),
  m_sg_gp_thickness(30),
//  m_sg_gp_transition(30), // always the same as the thickness
  m_supergrid_damping_coefficient(0.04),

  m_minJacobian(0.),
  m_maxJacobian(0.),

  m_energy_log(false),
  m_energy_print(false),
  m_energy_logfile("energy.dat"),
  m_saved_energy(0.0),

  m_do_timing(false),
  m_timing_print(0),
  m_output_detailed_timing(false),
  m_output_load(false),

  m_projection_cycle(1000),
  m_checkfornan(false),

  m_error_log_file("twilight_errors.dat"),
  m_error_log(false),
  m_error_print(true),
  m_inner_loop(9),
  m_topography_exists(false),
  m_useVelocityThresholds(false),
  m_vpMin(0.),
  m_vsMin(0.),
  m_grid_interpolation_order(0),
  m_zetaBreak(0.95),
  m_global_xmax(0.),
  m_global_ymax(0.),
  m_global_zmax(0.),
  m_global_zmin(0.),
  m_ghost_points(2), // for 4th order stencils
  m_ppadding(2),
  m_ext_ghost_points(0), // extra width for 6th order 
                         // discretization of metric at a source 
  //  m_ghost_points(3), // for 6th order stencils
  //  m_ppadding(3),

  mLonOrigin(-118.0), // NTS
  mLatOrigin(37.0), // NTS
  mGeoAz(0.0), // x=North, y=East
  //  mDefaultLocation(true),
  mMetersPerDegree(111319.5), // per degree in Latitude...
  mMetersPerLongitude(87721.0), // approximate for Lat=38 degrees
  mConstMetersPerLongitude(false),

// command limitfrequency
  m_limit_frequency(false), 
  m_ppw(15), 
  m_frequency_limit(1e38), // will hold min(Vs/h)/PPW

// command prefilter
  m_prefilter_sources(false), 
  m_filter_observations(false), 
  m_filter_ptr(0),
  m_filterobs_ptr(0),

  mPrintInterval(100),
  m_matrices_decomposed(false),
  m_citol(1e-3),
  m_cimaxiter(20),
  m_intp_conservative(true),
  mMaterialExtrapolate(0),

  m_use_attenuation(false),
  m_number_mechanisms(0),
  m_velo_omega(-1.0),
  m_min_omega(-1.0),
  m_max_omega(-1.0),
//  m_do_geodynbc(false),
  m_att_use_max_frequency(false),
  m_att_ppw(15),
  m_inverse_problem(a_invproblem),
  m_maxit(0),
  m_maxrestart(0),
  m_iniguess_pos(false),
  m_iniguess_t0fr(false),
  m_iniguess_mom(false),
  m_output_initial_seismograms(false),
  m_compute_scalefactors(false),
  m_cgstepselection(0),
  m_cgvarcase(0),
  m_cgfletcherreeves(true),
  m_do_linesearch(true),
  m_utc0set(false),
  m_utc0isrefevent(false),
  m_opttest(0),
  mEtreeFile(NULL),
  m_perturb(0),
  m_iperturb(1),
  m_jperturb(1),
  m_kperturb(1)
{
   MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
   MPI_Comm_size(MPI_COMM_WORLD, &m_nProcs);

// initialize the boundary condition array
   for (int i=0; i<6; i++)
   {
     mbcGlobalType[i] = bNone;
   }

   m_proc_array[0]=0;
   m_proc_array[1]=0;

// read the input file and setup the simulation object
   if (parseInputFile( a_GlobalSources, a_GlobalTimeSeries ))
     mParsingSuccessful = true;

// AP: need to figure out a better way of handling these error log files
   // char fname[100];
   // sprintf(fname,"sw4-error-log-p%i.txt", m_myRank);
   // msgStream.open(fname);
   
}

// Destructor
EW::
~EW()
{
//  msgStream.close();
}

//-----------------------------------
bool EW::isInitialized()
{
  return (mIsInitialized && mSourcesOK);
}
  
//-----------------------------------
bool EW::wasParsingSuccessful()
{
  return mParsingSuccessful;
}
  
//-----------------------------------------------------------------------
void EW::printTime( int cycle, double t, bool force ) const 
{
   if (!mQuiet && proc_zero() && (force || mPrintInterval == 1 ||
			(cycle % mPrintInterval) == 1 ||
			cycle == 1) )
// string big enough for >1 million time steps 
      printf("Time step %7i  t = %15.7e\n", cycle, t);
}
//-----------------------------------------------------------------------
void EW::printPreamble(vector<Source*> & a_Sources) const 
{
   stringstream msg;

   if (!mQuiet && proc_zero())
   {
      msg << "============================================================" << endl
          << " Running program on " << m_nProcs << " MPI tasks" << " using the following data: " << endl << endl
          << " Start Time = " << mTstart << " Goal Time = ";

      if (mTimeIsSet)
         msg << mTmax << endl;
      else
         msg << mNumberOfTimeSteps*mDt << endl;
      
      msg << " Number of time steps = " << mNumberOfTimeSteps << " dt: " << mDt << endl;
      
      if (mVerbose)
      {
	msg << endl;
	msg << "============================================================" << endl;
	msg << " Global boundary conditions " << endl;
	const char* side_names[6]={"x=0   ","x=xMax","y=0   ","y=yMax","z=topo","z=zMax"};
	for( int side = 0 ; side < 6 ; side++ )
	{
	  msg << "      ";
	  msg << side_names[side] << " " << bc_name(mbcGlobalType[side]) << "\n";
	}
	msg << endl;

	if (mHomogeneous)
	  msg << " Assuming Mu and Lambda to be constant within each grid  " << endl;
         
	 //         if (mForcing == 1 || mForcing == 2 || mForcing == 5)
	 //            msg << endl << " Second order Dirichlet boundary condition, gamma=" << mEBDirichletRegCoeff << endl;
	 //         else if (mForcing == 3 || mForcing == 4 || mForcing == 6)
	 //            msg << endl << " Second order Neumann boundary condition" << endl;
         

	if ( mVerbose >= 4 )
	  cout << " The following point sources are used: " << endl;
      }
      cout << msg.str();
   }
   MPI_Barrier(MPI_COMM_WORLD);

   cout.flush(); cerr.flush();
      
   // m0 values in each source command gets added up. This number is called the "Total seismic moment" 
   // and should be printed to stdout with the unit Nm (Newton-meter). If that number is >0, you should 
   // also print Mw = 2/3 *(log10(M0) - 9.1). That is the moment magnitude (dimensionless). 
      
   if( proc_zero() )
   {
     if (m_twilight_forcing)
     {
       cout << "-----------------------------------------------------------" << endl;
       cout << "Twilight zone testing (aka method of manufactured solution)" << endl;
       cout << "Parameters:" << endl;
       cout << "  omega = " << m_twilight_forcing->m_omega << endl;
       cout << "  c = " << m_twilight_forcing->m_c << endl;
       cout << "  phase = " << m_twilight_forcing->m_phase << endl;
       cout << "  mat-omega = " << m_twilight_forcing->m_momega << endl;
       cout << "  mat-phase = " << m_twilight_forcing->m_mphase << endl;
       cout << "  amprho = " << m_twilight_forcing->m_amprho << endl;
       cout << "  amplambda = " << m_twilight_forcing->m_amplambda << endl;
       cout << "  ampmu = " << m_twilight_forcing->m_ampmu << endl;
       cout << "-----------------------------------------------------------" << endl;
     }
     else if (m_lamb_test)
     {
       double fx, fy, fz, xs, ys, zs;
       a_Sources[0]->getForces( fx, fy, fz );
       xs = a_Sources[0]->getX0( );
       ys = a_Sources[0]->getY0( );
       zs = a_Sources[0]->getZ0( );
       string tfun;
       if( a_Sources[0]->getTfunc() == iVerySmoothBump )
	 tfun = "VerySmoothBump";
       else if( a_Sources[0]->getTfunc() == iC6SmoothBump )
	 tfun = "C6SmoothBump";

       cout << "-----------------------------------------------------------" << endl;
       cout << "Lamb's problem testing" << endl;
       cout << "Parameters:" << endl;
       cout << "  Cp = " << m_lamb_test->m_cp << endl;
       cout << "  Cs = " << m_lamb_test->m_cs << endl;
       cout << "  Rho = " << m_lamb_test->m_rho << endl;       
       cout << "  (xs, ys, zs) = " << xs << ", " << ys << ", " << zs << endl;       
       cout << "  (fx, fy, fz) = " << fx << ", " << fy << ", " << fz << endl;       
       cout << "  Source time fcn = " << tfun << endl;       
       cout << "-----------------------------------------------------------" << endl;
     }
     else
     {
       double myM0Sum = 0;
       int numsrc = 0; //, ignoredSources=0;
       for (unsigned int i=0; i < a_Sources.size(); ++i)
       {
         if (a_Sources[i]->isMomentSource())
	 {
// Note that proc 0 doen't know of all sources that need to be ignored
//	   if (!a_Sources[i]->ignore() ) 
	 {
	   numsrc++;
	   myM0Sum += a_Sources[i]->getAmplitude();
	 }
	 //	   else
	 //	     ignoredSources++;
	 }
	 
       }
       if (!mQuiet)
       {
	 stringstream msg2;
	 msg2 << endl
	      << "-----------------------------------------------------------------------" << endl
	      << "  Total seismic moment (M0): " << myM0Sum << " Nm " << endl;
	 if (myM0Sum > 0)
	   msg2 <<  "  Moment magnitude     (Mw): " << (2./3.)*(log10(myM0Sum) - 9.1)  << endl;
	 msg2 << "  Number of sources " << numsrc << endl;
	 msg2 << "-----------------------------------------------------------------------" << endl;
	 cout << msg2.str();
       }
     } // standard run
   } // end if proc_zero()
   
}

//-----------------------------------------------------------------------
void EW::switch_on_checkfornan()
{
   m_checkfornan = true;
}

//-----------------------------------------------------------------------
void EW::assign_local_bcs( )
{
// This routine assigns m_bcType[g][b], b=0,1,2,3, based on mbcGlobalType, taking parallel overlap boundaries into account

  int g, b, side;
  int top=mNumberOfGrids-1; // index of the top grid in the arrays m_iStart, m_iEnd, etc
  
// horizontal bc's are the same for all grids
  for( g= 0 ; g < mNumberOfGrids ; g++ )
  {
// start by copying the global bc's
    for (b=0; b<=3; b++)
      m_bcType[g][b] = mbcGlobalType[b];
  
    if (m_iStart[top]+m_ghost_points > 1)
    {
      m_bcType[g][0] = bProcessor;
    }
    if (m_iEnd[top]-m_ghost_points < m_global_nx[top])
    {
      m_bcType[g][1] = bProcessor;
    }

// for a periodic domain, we need to change all bc to bProcessor if more than 1 process is used in that direction
    if (m_doubly_periodic && m_proc_array[0] > 1)
    {
      m_bcType[g][0] = bProcessor;
      m_bcType[g][1] = bProcessor;
    }    

    if (m_jStart[top]+m_ghost_points > 1)
    {
      m_bcType[g][2] = bProcessor;
    }
    if (m_jEnd[top]-m_ghost_points < m_global_ny[top])
    {
      m_bcType[g][3] = bProcessor;
    }

// for a periodic domain, we need to change all bc to bProcessor if more than 1 process is used in that direction
    if (m_doubly_periodic && m_proc_array[1] > 1)
    {
      m_bcType[g][2] = bProcessor;
      m_bcType[g][3] = bProcessor;
    }

  }
  
// vertical bc's are interpolating except at the bottom and the top, where they equal the global conditions
  for( g= 0 ; g < mNumberOfGrids ; g++ )
  {
    m_bcType[g][4] = bInterpolate;
    m_bcType[g][5] = bInterpolate;
  }
  m_bcType[top][4] = mbcGlobalType[4];
  m_bcType[0][5] = mbcGlobalType[5];

// Find out which boundaries need one sided approximation in mixed derivatives
  for( g= 0 ; g < mNumberOfGrids ; g++ )
  {
    for(side=0 ; side < 6 ; side++ )
    {
// add energy conserving mesh coupling condition
      m_onesided[g][side] = (m_bcType[g][side] == bStressFree); 
    }
  }

// one-sided cross-terms at conservative interpolation
  if( m_intp_conservative ) 
  {
    for( g= 0 ; g < mNumberOfCartesianGrids ; g++ )
    {
      for(side=0 ; side < 6 ; side++ )
      {
	if(m_bcType[g][side] == bInterpolate ) // THIS MAKES THE STENCIL ONE-SIDED AT THE CURVILINEAR-CARTESIAN INTERFACE!
	  m_onesided[g][side] = 1;
      }
    }

// Must be careful not to use one-sided formulas at the curvilinear-cartesian interface!
    g = mNumberOfCartesianGrids - 1; // top Cartesian grid
    side = 4; // low-k
    if (topographyExists())
    {
      if(m_bcType[g][side] == bInterpolate ) // THIS MAKES THE STENCIL CENTERD AT THE CURVILINEAR-CARTESIAN INTERFACE!
      {
	m_onesided[g][side] = 0;
// tmp
// 	if (proc_zero() && mVerbose >= 1)
// 	  printf("******* Reverting the cross-term to CENTERED on the Cartesian-curvilinear interface!!!!!!!!!\n");
	
      }
      
    }
        
  } // end if
  
}

//-----------------------------------------------------------------------
// Note that the padding cell array is no longer needed.
// use m_iStartInt[g], m_iEndInt[g] to get the range of interior points
void EW::initializePaddingCells()
{
  int g = mNumberOfGrids-1;
  
   for (int aa = 0; aa < 4; aa++)
   {
     if (m_bcType[g][aa] == bProcessor)
     {
       m_paddingCells[aa] = m_ppadding;
     }
     else
     {
       m_paddingCells[aa] = m_ghost_points;
     }
   }
}

//-----------------------------------------------------------------------
bool EW::proc_zero() const
{
  return (m_myRank == 0);
}

//-----------------------------------------------------------------------
int EW::no_of_procs() const
{
  return m_nProcs;
}


//-----------------------------------------------------------------------
string EW::bc_name( const boundaryConditionType bc ) const
{
   string retval;
   if( bc == bStressFree )
      retval = "free surface";
   else if( bc == bDirichlet )
      retval = "dirichlet";
   else if( bc == bSuperGrid )
      retval = "supergrid";
   else if( bc == bPeriodic )
      retval = "periodic";
   else if( bc == bInterpolate )
      retval = "interpolation";
   else if( bc == bProcessor )
      retval = "processor";
   else if( bc == bNone )
      retval = "none";
   return retval;
}

//-----------------------------------------------------------------------
bool EW::getDepth( double x, double y, double z, double & depth)
{
// get the depth below the free surface
  bool success=false;
  
  if (!topographyExists())
  {
    depth = z;
    success = true;
  }
  else
  {
// topography 
     double zMinTilde;
     int gCurv = mNumberOfGrids - 1;
     double h = mGridSize[gCurv];
     double q = x/h + 1.0;
     double r = y/h + 1.0;

// define the depth for ghost points (in x or y) to equal the depth on the nearest boundary point
     double qMin = 1.0;
     double qMax = (double) m_global_nx[gCurv];
     double rMin = 1.0;
     double rMax = (double) m_global_ny[gCurv];

     if (q<qMin) q=qMin;
     if (q>qMax) q=qMax;
     if (r<rMin) r=rMin;
     if (r>rMax) r=rMax;

// // evaluate elevation of topography on the grid (smoothed topo)
    success=true;
    if (!interpolate_topography(q, r, zMinTilde, true))
    {
      cerr << "ERROR: getDepth: Unable to evaluate topography for x=" << x << " y= " << y << " on proc # " << getRank() << endl;
      // cerr << "q=" << q << " r=" << r << " qMin=" << qMin << " qMax=" << qMax << " rMin=" << rMin << " rMax=" << rMax << endl;
      // cerr << "Setting elevation of topography to ZERO" << endl;
      success = false;
//      zMinTilde = 0;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    depth = z-zMinTilde;
  }
  return success;
}

//-----------------------------------------------------------------------
void EW::computeCartesianCoord(double &x, double &y, double lon, double lat)
{
  // -----------------------------------------------------------------
  // Compute the cartesian coordinate given the geographic coordinate
  // -----------------------------------------------------------------
  if( m_geoproj == 0 )
   //  // compute x and y
  {
     double deg2rad = M_PI/180.0;
     double phi = mGeoAz * deg2rad;
     //     x = mMetersPerDegree*(cos(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*sin(phi));
     //     y = mMetersPerDegree*(-sin(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*cos(phi));
     if (mConstMetersPerLongitude)
     {
	x = mMetersPerDegree*cos(phi)*(lat-mLatOrigin)    + mMetersPerLongitude*(lon-mLonOrigin)*sin(phi);
	y = mMetersPerDegree*(-sin(phi))*(lat-mLatOrigin) + mMetersPerLongitude*(lon-mLonOrigin)*cos(phi);
     }
     else
     {
	x = mMetersPerDegree*(cos(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*sin(phi));
	y = mMetersPerDegree*(-sin(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*cos(phi));
     }
  }
  else
     m_geoproj->computeCartesianCoord(x,y,lon,lat);
   // Test with proj4
  //  projUV lonlat, xy;
  //  lonlat.u = lon*deg2rad;
  //  lonlat.v = lat*deg2rad;
  //  xy = pj_fwd( lonlat, m_projection );
  //  if( xy.u == HUGE_VAL )
  //     cout << "ERROR: computeCartesianCoord pj_fwd failed with message " << pj_strerrno(pj_errno) << endl;
  //  xy.u -= m_xoffset;
  //  xy.v -= m_yoffset;
  //  x =  xy.u*sin(phi) + cos(phi)*xy.v;
  //  y =  xy.u*cos(phi) - sin(phi)*xy.v;
}

//-----------------------------------------------------------------------
void EW::computeGeographicCoord(double x, double y, double & longitude, double & latitude)
{
  // conversion factor between degrees and radians
   if( m_geoproj == 0 )
   {
      double deg2rad = M_PI/180.0;
      double phi = mGeoAz * deg2rad;
      // Compute the latitude
      latitude = mLatOrigin + 
	 (x*cos(phi) - y*sin(phi))/mMetersPerDegree;
      // Compute the longitude
      if (mConstMetersPerLongitude)
      {
	 longitude = mLonOrigin + 
	    (x*sin(phi) + y*cos(phi))/(mMetersPerLongitude);
      }
      else
      {
	 longitude = mLonOrigin + 
	    (x*sin(phi) + y*cos(phi))/(mMetersPerDegree*cos(latitude*deg2rad));
      }
   }
   else
      m_geoproj->computeGeographicCoord( x, y, longitude, latitude );
  // Test with proj4
   // projUV lonlat, xy;
   // xy.u = x*sin(phi) + y*cos(phi) + m_xoffset;
   // xy.v = x*cos(phi) - y*sin(phi) + m_yoffset;
   // lonlat = pj_inv( xy, m_projection );
   //  if( lonlat.u == HUGE_VAL )
   //     cout << "ERROR: computeGeographicCoord, pj_inv failed with message " << pj_strerrno(pj_errno) << endl;
   // longitude = lonlat.u/deg2rad;
   // latitude  = lonlat.v/deg2rad;
}

//-------------------------------------------------------
void EW::computeNearestGridPoint(int & a_i, 
                                   int & a_j, 
                                   int & a_k,
                                   int & a_g, // grid on which indices are located
                                   double a_x, 
                                   double a_y, 
                                   double a_z)
{
  bool breakLoop = false;
  
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      if (a_z > m_zmin[g] || g == mNumberOfGrids-1) // We can not trust zmin for the curvilinear grid, since it doesn't mean anything
        {
          a_i = (int)floor(a_x/mGridSize[g])+1;
          if (a_x-((a_i-0.5)*mGridSize[g]) > 0.) (a_i)++;
          
          a_j = (int)floor(a_y/mGridSize[g])+1;
          if (a_y-((a_j-0.5)*mGridSize[g]) > 0.) (a_j)++;
          
          a_k    = (int)floor((a_z-m_zmin[g])/mGridSize[g])+1;  //Note: this component will be garbage for g=curvilinear grid
          if (a_z-(m_zmin[g]+(a_k-0.5)*mGridSize[g]) > 0.)   (a_k)++;
          
          a_g = g                                        ;
          
          breakLoop = true;
        }
      else if (a_z == m_zmin[g]) // testing for equality between doubles is kind of pointless...
        {
           // Point is located on top surface if g=finest grid, else the location is on
	   // a grid/grid interface, and point is flagged as located on the finer (upper) grid.
          if (g == mNumberOfGrids-1)
            {
              a_i = (int)floor(a_x/mGridSize[g])+1;
              if (a_x-((a_i-0.5)*mGridSize[g]) > 0.) (a_i)++;
              
              a_j = (int)floor(a_y/mGridSize[g])+1;
              if (a_y-((a_j-0.5)*mGridSize[g]) > 0.) (a_j)++;
              
              a_k = 1;
              
              a_g = g;
            }
          else
            {
              a_i = (int)floor(a_x/mGridSize[g+1])+1;
              if (a_x-((a_i-0.5)*mGridSize[g+1]) > 0.) (a_i)++;
              
              a_j = (int)floor(a_y/mGridSize[g+1])+1;
              if (a_y-((a_j-0.5)*mGridSize[g+1]) > 0.) (a_j)++;
              
              a_k = (int)floor((a_z-m_zmin[g+1])/mGridSize[g+1])+1; // Here, I know I am on a grid line
              
              a_g = g+1                                    ;
            }
          breakLoop = true;
        }
      
      if (breakLoop)
        {
              break;
        } 
    }
  
//  if (m_topography_exists && (a_g == mNumberOfGrids-1)) // The curvilinear grid will always be the one with the highest number. 
//    {
// tmp
//      printf("EW/computeNearestGridPt: You are in the curvilinear part of the grid, but we do compute the gridpt index using only the Cartesian grid\n");
//    }

  if (!m_topography_exists || (m_topography_exists && a_g < mNumberOfCartesianGrids))
    {
      VERIFY2(a_i >= 1-m_ghost_points && a_i <= m_global_nx[a_g]+m_ghost_points,
              "Grid Error: i (" << a_i << ") is out of bounds: ( " << 1 << "," 
              << m_global_nx[a_g] << ")" << " x,y,z = " << a_x << " " << a_y << " " << a_z);
      VERIFY2(a_j >= 1-m_ghost_points && a_j <= m_global_ny[a_g]+m_ghost_points,
              "Grid Error: j (" << a_j << ") is out of bounds: ( " << 1 << ","
              << m_global_ny[a_g] << ")" << " x,y,z = " << a_x << " " << a_y << " " << a_z);
      VERIFY2(a_k >= m_kStart[a_g] && a_k <= m_kEnd[a_g],
              "Grid Error: k (" << a_k << ") is out of bounds: ( " << 1 << "," 
              << m_kEnd[a_g]-m_ghost_points << ")" << " x,y,z = " << a_x << " " << a_y << " " << a_z);
    }
}

void EW::computeNearestLowGridPoint(int & a_i, 
                                      int & a_j, 
                                      int & a_k,
                                      int & a_g, // grid on which indices are located
                                      double a_x, 
                                      double a_y, 
                                      double a_z)
{
  bool breakLoop = false;
  
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      if (a_z > m_zmin[g] || g == mNumberOfGrids-1) // We can not trust zmin for the curvilinear grid, since it doesn't mean anything
        {
          a_i = (int)floor(a_x/mGridSize[g])+1;
	  //          VERIFY(a_x-((a_i-0.5)*mGridSize[g]) <= 0.);
          
          a_j = (int)floor(a_y/mGridSize[g])+1;
	  //          VERIFY(a_y-((a_j-0.5)*mGridSize[g]) <= 0.);
          
          a_k    = (int)floor((a_z-m_zmin[g])/mGridSize[g])+1;
	  //          VERIFY(a_z-(m_zmin[g]+(a_k-0.5)*mGridSize[g]) <= 0.);
          
          a_g = g                                        ;
          
          breakLoop = true;
        }
      else if (a_z == m_zmin[g])
        {
          if (g == mNumberOfGrids-1)
            {
              a_i = (int)floor(a_x/mGridSize[g])+1;
	      //              VERIFY(a_x-((a_i-0.5)*mGridSize[g]) <= 0.);
              
              a_j = (int)floor(a_y/mGridSize[g])+1;
	      //              VERIFY(a_y-((a_j-0.5)*mGridSize[g]) <= 0.);
              
              a_k = 1;
              
              a_g = g;
            }
          else
            {
              a_i = (int)floor(a_x/mGridSize[g+1])+1;
	      //              VERIFY(a_x-((a_i-0.5)*mGridSize[g+1]) <= 0.);
              
              a_j = (int)floor(a_y/mGridSize[g+1])+1;
	      //              VERIFY(a_y-((a_j-0.5)*mGridSize[g+1]) <= 0.);
              
              a_k = (int)floor((a_z-m_zmin[g+1])/mGridSize[g+1])+1; // Here, I know I am on a grid line
              
              a_g = g+1                                    ;
            }
          breakLoop = true;
        }
      
      if (breakLoop)
        {
              break;
        } 
    }
  
  VERIFY2(a_i >= 1 && a_i <= m_global_nx[a_g],
          "Grid Error: i (" << a_i << ") is out of bounds: ( " << 1 << "," 
          << m_global_nx[a_g] << ")");
  VERIFY2(a_j >= 1 && a_j <= m_global_ny[a_g],
          "Grid Error: j (" << a_j << ") is out of bounds: ( " << 1 << ","
          << m_global_ny[a_g] << ")");
  if( a_k > m_kEnd[a_g]-m_ghost_points )
     a_k = m_kEnd[a_g]-m_ghost_points;
  if( a_k < 1 )
     a_k = 1;
  //  VERIFY2(a_k >= 1 && a_k <= m_kEnd[a_g]-m_ghost_points,
  //          "Grid Error: k (" << a_k << ") is out of bounds: ( " << 1 << "," 
  //          << m_kEnd[a_g]-m_ghost_points << ")");
}


//-----------------------------------------------------------------------
bool EW::interior_point_in_proc(int a_i, int a_j, int a_g)
{
// NOT TAKING PARALLEL GHOST POINTS INTO ACCOUNT!
// Determine if grid point with index (a_i, a_j) on grid a_g is an interior grid point on this processor 

   bool retval = false;
   if (a_g >=0 && a_g < mNumberOfGrids){
     retval = (a_i >= m_iStartInt[a_g]) && (a_i <= m_iEndInt[a_g]) &&   
       (a_j >= m_jStartInt[a_g]) && (a_j <= m_jEndInt[a_g]);
   }
   return retval; 
}

//-----------------------------------------------------------------------
bool EW::point_in_proc(int a_i, int a_j, int a_g)
{
// TAKING PARALLEL GHOST POINTS INTO ACCOUNT!
// Determine if grid point with index (a_i, a_j) on grid a_g is a grid point on this processor 

   bool retval = false;
   if (a_g >=0 && a_g < mNumberOfGrids){
     retval = (a_i >= m_iStart[a_g] && a_i <= m_iEnd[a_g] &&   
               a_j >= m_jStart[a_g] && a_j <= m_jEnd[a_g] );
   }

   return retval; 
}

//-----------------------------------------------------------------------
bool EW::point_in_proc_ext(int a_i, int a_j, int a_g)
{
// TAKING PARALLEL GHOST POINTS+EXTRA GHOST POINTS INTO ACCOUNT!
// Determine if grid point with index (a_i, a_j) on grid a_g is a grid point on this processor 

   bool retval = false;
   if (a_g >=0 && a_g < mNumberOfGrids){
     retval = (a_i >= m_iStart[a_g]-m_ext_ghost_points && a_i <= m_iEnd[a_g]+m_ext_ghost_points &&   
               a_j >= m_jStart[a_g]-m_ext_ghost_points && a_j <= m_jEnd[a_g]+m_ext_ghost_points );
   }
   return retval; 
}

//-----------------------------------------------------------------------
void EW::getGlobalBoundingBox(double bbox[6])
{
  bbox[0] = 0.;
  bbox[1] = m_global_xmax;
  bbox[2] = 0.;
  bbox[3] = m_global_ymax;
  bbox[4] = m_global_zmin;
  bbox[5] = m_global_zmax;
}

//-----------------------------------------------------------------------
void EW::setGMTOutput(string filename, string wppfilename)
{
  mGMTFileName = filename;
  mWriteGMTOutput = true;

//  mWPPFileName = wppfilename;
}

//-----------------------------------------------------------------------
void EW::saveGMTFile( vector<Source*> & a_GlobalUniqueSources )
{
// this routine needs to be updated
   if (!mWriteGMTOutput) return;

   if (proc_zero())
   {
      stringstream contents;
      contents << "#!/bin/csh\n\n" 
               << "gmtset PLOT_DEGREE_FORMAT D\n"
               << "gmtset COLOR_MODEL HSV\n"
               << "gmtset PAPER_MEDIA letter\n"
               << "gmtset PAGE_ORIENTATION portrait\n"
               << "gmtset MEASURE_UNIT inch\n" 
               << endl;
      // grab these from grid
      double latNE,lonNE,latSW,lonSW,latSE,lonSE,latNW,lonNW;
      computeGeographicCoord(0.0,           0.0,           lonSW, latSW);
      computeGeographicCoord(m_global_xmax, 0.0,           lonSE, latSE);
      computeGeographicCoord(m_global_xmax, m_global_ymax, lonNE, latNE);
      computeGeographicCoord(0.0,           m_global_ymax, lonNW, latNW);
     
      // Round up/down
      double minx = min(lonSW, min(lonSE, min(lonNE, lonNW)));
      double maxx = max(lonSW, max(lonSE, max(lonNE, lonNW)));
      double miny = min(latSW, min(latSE, min(latNE, latNW)));
      double maxy = max(latSW, max(latSE, max(latNE, latNW))); 
      double margin = 0.1*fabs(maxy-miny);
      
// tmp
   printf("margin = %e\n", margin);

//      GeographicCoord eNW, eNE, eSW, eSE;
      
#ifdef ENABLE_ETREE
      if (mEtreeFile != NULL)
      {
	 //        mEtreeFile->getGeoBox()->getBounds(eNW, eNE, eSW, eSE);
// correct these as above (remove +/- 1        
//        minx = (min(eSW.getLongitude()-1,min(eSE.getLongitude()-1,min(eNE.getLongitude()-1,eNW.getLongitude()-1))));
//        maxx = (max(eSW.getLongitude()+1,max(eSE.getLongitude()+1,max(eNE.getLongitude()+1,eNW.getLongitude()+1))));
//        miny = (min(eSW.getLatitude()-1,min(eSE.getLatitude()-1,min(eNE.getLatitude()-1,eNW.getLatitude()-1))));
//        maxy = (max(eSW.getLatitude()+1,max(eSE.getLatitude()+1,max(eNE.getLatitude()+1,eNW.getLatitude()+1)))); 
         mEtreeFile->getbox( miny, maxy, minx, maxx );
	 minx -= 1;
	 maxx += 1;
	 miny -= 1;
	 maxy += 1;
      }
#endif
      
      contents << "# Region will need to be adjusted based on etree/grid values" << endl
               << "set REGION = " << minx-margin << "/" << maxx+margin << "/" << miny-margin << "/" << maxy+margin << endl
               << endl
               << "set SCALE = 6.0" << endl
               << endl
               << "# These commands are good if you have access to " << endl
               << "# a topography database file for the region modeled " << endl
               << "# Note:  if you uncomment these, adjust the -O -K, etc." << endl
               <<" #######################################################" << endl
               << "#grdraster 2 -R$REGION -I0.5m -Gwpp_topo.grd" << endl
               << "#grdgradient wpp_topo.grd -Gwpp_topo_shade.grd -A270 -Nt -M " << endl
               << "#grd2cpt wpp_topo.grd -Ctopo -Z >! wpptopo.cpt" << endl
               << "#grdimage wpp_topo.grd -R$REGION -JM$SCALE -Cwpptopo.cpt -Iwpp_topo_shade.grd -P -K >! plot.ps" << endl
               <<" #######################################################" << endl
               << "pscoast -R$REGION -JM$SCALE -Bf0.025a0.05 -Dfull -S100,200,255 -A2000 -W3 -N1t3 -N2t2a -K >! plot.ps" << endl << endl
               << "# computational grid region..." << endl;
      
      // Write out gridlines
      contents << "psxy -R$REGION -JM$SCALE -W10/255/255/0ta -O -K <<EOF>> plot.ps" << endl
               << lonSW << " " << latSW << endl
               << lonSE << " " << latSE << endl
               << lonNE << " " << latNE << endl
               << lonNW << " " << latNW << endl  
               << lonSW << " " << latSW << endl
               << "EOF" << endl << endl;
      
#ifdef ENABLE_ETREE
      if (mEtreeFile != NULL)
      {
// Consider Etree bounds also
//         GeographicCoord eNW, eNE, eSW, eSE;
//         mEtreeFile->getGeoBox()->getBounds(eNW, eNE, eSW, eSE);
         double elatSE, elonSE, elatSW, elonSW, elatNE, elonNE, elatNW, elonNW;
         mEtreeFile->getcorners( elatSE, elonSE, elatSW, elonSW, elatNE, elonNE, elatNW, elonNW );
         contents << "# Etree region: " << mEtreeFile->getFileName() << endl
                  << "psxy -R$REGION -JM$SCALE -W5/255/255/0ta -O -K <<EOF>> plot.ps" << endl
                  << elonNW << " " << elatNW << endl
                  << elonNE << " " << elatNE << endl
                  << elonSE << " " << elatSE << endl
                  << elonSW << " " << elatSW << endl
                  << elonNW << " " << elatNW << endl
                  << "EOF" << endl << endl;
	 //         contents << "# Etree region: " << mEtreeFile->getFileName() << endl
	 //                  << "psxy -R$REGION -JM$SCALE -W5/255/255/0ta -O -K <<EOF>> plot.ps" << endl
	 //                  << eNW.getLongitude() << " " << eNW.getLatitude() << endl
	 //                  << eNE.getLongitude() << " " << eNE.getLatitude() << endl
	 //                  << eSE.getLongitude() << " " << eSE.getLatitude() << endl
	 //                  << eSW.getLongitude() << " " << eSW.getLatitude() << endl
	 //	             << eNW.getLongitude() << " " << eNW.getLatitude() << endl
	 //                  << "EOF" << endl << endl;
      }
#endif
      
      if (a_GlobalUniqueSources.size() > 0)
      {
         contents << "# Sources... " << endl
	          << "cat << EOF >! event.d" << endl;
         
         for (int i=0; i < a_GlobalUniqueSources.size(); ++i)
         {
           double latSource,lonSource;

           computeGeographicCoord(a_GlobalUniqueSources[i]->getX0(), a_GlobalUniqueSources[i]->getY0(),
                                  lonSource ,latSource);
//  should name the event better
	   contents << lonSource << " " << latSource << " EVENT-NAME  CB" << endl;
         }
         contents << "EOF" << endl;
	 contents << "psxy -R -J -O -K -Sc0.1 -Gred -W0.25p event.d >> plot.ps" << endl;
         contents << "awk '{print $1, $2, 12, 1, 9, $4, $3}' event.d | pstext -R -J -O -D0.2/0.2v -Gred -N -K >> plot.ps" 
	   << endl << endl;
      }
      
      int numStations = 0;
      stringstream stationstr;
      stationstr << "# Stations... " << endl;  
      stationstr << "cat << EOF >! stations.d " << endl;
      // Write stations by rereading the WPP input file, since some might
      // live outside the grid...
      ifstream sw4InputFile(mName.c_str());
      if (!sw4InputFile.is_open())
         contents << "# Error re-opening input file, skipping stations" << endl;
      else
      {
         char buffer[256];
         while (!sw4InputFile.eof())
         { 
            sw4InputFile.getline(buffer, 256);
            if (startswith("receiver", buffer))
            {
               numStations += 1;
               bool cartCoordSet = false;
               bool gridPointSet = false;
               bool geoCoordSet = false;
               bool statSet = false;
               string name="null";
               double x=0.0, y=0.0, z=0.0;
               double lat=0.0, lon=0.0;
               int i=0,j=0,k=0;
               // Get location and write to file
               char* token = strtok(buffer, " \t");   
               token = strtok(NULL, " \t"); // skip sac
               while (token != NULL)
               {
                  // while there are tokens in the string still
                  // NOTE: we skip all verify stmts as these have
                  //       already been checked during initial parsing
                  if (startswith("#", token) || startswith(" ", buffer))
                     // Ignore commented lines and lines with just a space.
                     break;
                  if (startswith("x=", token))
                  {
                     token += 2; // skip x=
                     cartCoordSet = true;
                     x = atof(token);
                  }
                  else if (startswith("y=", token))
                  {
                     token += 2; // skip y=
                     cartCoordSet = true;
                     y = atof(token);
                  }
                  else if (startswith("z=", token))
                  {
                     token += 2; // skip z=
                     cartCoordSet = true;
                     z = atof(token);
                  }
                  else if (startswith("lat=", token))
                  {
                     token += 4; // skip lat=
                     lat = atof(token);
                     geoCoordSet = true;
                  }
                  else if (startswith("lon=", token))
                  {
                     token += 4; // skip lon=
                     lon = atof(token);
                     geoCoordSet = true;
                  }
                  else if (startswith("depth=", token))
                  {
                     token += 6; // skip depth=
                     z = atof(token);
                     geoCoordSet = true;
                  }
                  else if (startswith("sta=", token))
                  {
                     token += 4;
                     name = token;
                     statSet = true;
                  }
                  else if (startswith("file=", token) && !statSet)
                  {
                     token += 5;
                     name = token;
                  }
                  
                  token = strtok(NULL, " \t");
               }
               
               VERIFY(cartCoordSet || geoCoordSet);

               if (!geoCoordSet && cartCoordSet)
               {
                 computeGeographicCoord(x, y, lon, lat);
               }
               
               // Now have location
               stationstr << lon << " " << lat << " " << name << " CB" << endl; 
            } // token on sac line
         } // line in ew file
      }
      
      stationstr << "EOF" << endl << endl;
      
      stationstr << "# plot station names" << endl
                 << "psxy -R -J -O -K -St0.1 -Gblue -W0.25p stations.d >> plot.ps" << endl
                 << "awk '{print $1, $2, 12, 1, 9, $4, $3}' stations.d | pstext -R -J -O -Dj0.3/0.3v -Gblue -N >> plot.ps" << endl;
      
      // Only write station info if there are stations.
      if (numStations > 0) contents << stationstr.str() << endl;

      contents << "/bin/mv plot.ps " << mName << ".ps" << endl;

      stringstream filename;
      filename << mPath << mGMTFileName;
      ofstream gmtfile(filename.str().c_str());
      if (gmtfile.is_open())
      {
	cout << "GMT file is open, about to write" << endl;
	gmtfile << contents.str();
	cout << "Wrote GMT file: " << filename.str() << endl;
      }
      else
      {
	cout << "Unable to open GMT file: " << filename.str() << endl;
      }
      
   } // proc 0
}

//-----------------------------------------------------------------------
void EW::print_execution_time( double t1, double t2, string msg )
{
   if( !mQuiet && proc_zero() )
   {
      double s = t2 - t1;
      int h = static_cast<int>(s/3600.0);
      s = s - h*3600;
      int m = static_cast<int>(s/60.0);
      s = s - m*60;
      cout << "   Execution time, " << msg << " ";
      if( h > 1 )
	 cout << h << " hours ";
      else if( h > 0 )
	 cout << h << " hour  ";

      if( m > 1 )
	 cout << m << " minutes ";
      else if( m > 0 )
	 cout << m << " minute  ";

      if( s > 0 )
	 cout << s << " seconds " ;
      cout << endl;
   }
}


//-----------------------------------------------------------------------
void EW::print_execution_times( double times[7] )
{
   double* time_sums =new double[7*no_of_procs()];
   MPI_Gather( times, 7, MPI_DOUBLE, time_sums, 7, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   if( !mQuiet && proc_zero() )
   {
      cout << "\n----------------------------------------" << endl;
      cout << "          Execution time summary " << endl;
      cout << "Processor  Total      BC total   Step       Image&Sac  Comm.ref   Comm.bndry BC impose  "
	   <<endl;
      cout.setf(ios::left);
      cout.precision(3);
      for( int p= 0 ; p < no_of_procs() ; p++ )
      {
         cout.width(11);
         cout << p;
         cout.width(11);
	 cout << time_sums[7*p+3];
	 cout.width(11);
	 cout << time_sums[7*p+1];
	 cout.width(11);
	 cout << time_sums[7*p];
	 cout.width(11);
	 cout << time_sums[7*p+2];
	 cout.width(11);
	 cout << time_sums[7*p+4];
	 cout.width(11);
	 cout << time_sums[7*p+5];
	 cout.width(11);
	 cout << time_sums[7*p+6];
         cout << endl;
      }
      cout.setf(ios::right);
      cout.precision(6);
      //
      // << "|" << time_sums[p*7+3] << "|\t" << time_sums[p*7+1] << "|\t" << time_sums[p*7]
      //	      << "|\t " << time_sums[7*p+2] << "|\t" << time_sums[p*7+4] << "|\t" << time_sums[p*7+5]
      //	      << "|\t" << time_sums[7*p+6]<<endl;
      cout << "----------------------------------------\n" << endl;
   }
   delete[] time_sums;
}

//-----------------------------------------------------------------------
void EW::finalizeIO()
{
  //  //  if (proc_zero() && mDebugIO )
  //  //    mEnergyFile.close();
  // for (unsigned int i = 0; i < mSACOutputFiles.size(); ++i)
  //   mSACOutputFiles[i].writeFile();
}

//-----------------------------------------------------------------------
void EW::default_bcs( )
{
   for( int side=0 ; side < 6 ; side++ )
      mbcGlobalType[side] = bSuperGrid;
   mbcGlobalType[4] = bStressFree; // low-z is normally free surface
}

//---------------------------------------------------------------------------
void EW::normOfDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, 
                           double &diffL2, double &xInf, vector<Source*>& a_globalSources )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast;
  int imin, imax, jmin, jmax, kmin, kmax;
  
  double *uex_ptr, *u_ptr, h, linfLocal=0, l2Local=0, diffInfLocal=0, diffL2Local=0;
  double xInfLocal=0, xInfGrid=0;
  double radius =-1, x0=0, y0=0, z0=0;

//tmp  
//   if (proc_zero())
//     printf("Inside normOfDifference\n");
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    uex_ptr  = a_Uex[g].c_ptr();
    u_ptr    = a_U[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];  

    if (mbcGlobalType[0] == bSuperGrid)
      imin = max(m_iStartInt[g], m_sg_gp_thickness+1);
    else
      imin = m_iStartInt[g];
  
    if (mbcGlobalType[1] == bSuperGrid)
      imax = min(m_iEndInt[g], m_global_nx[g] - m_sg_gp_thickness);
    else
      imax = m_iEndInt[g];

    if (mbcGlobalType[2] == bSuperGrid)
      jmin = max(m_jStartInt[g], m_sg_gp_thickness+1);
    else
      jmin = m_jStartInt[g];

    if (mbcGlobalType[3] == bSuperGrid)
      jmax = min(m_jEndInt[g], m_global_ny[g] - m_sg_gp_thickness);
    else
      jmax = m_jEndInt[g];

    if (mbcGlobalType[4] == bSuperGrid)
      kmin = max(m_kStartInt[g], m_sg_gp_thickness+1);
    else
      kmin = m_kStartInt[g];

    if (mbcGlobalType[5] == bSuperGrid)
      kmax = min(m_kEndInt[g], m_global_nz[g] - m_sg_gp_thickness);
    else
      kmax = m_kEndInt[g];

// tmp
//     printf("proc=%i, iS= %i, iE=%i, jS=%i, jE=%i, kS=%i, kE=%i\n", m_myRank, 
// 	   m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g]);
//     printf("proc=%i, if= %i, il=%i, jf=%i, jl=%i, kf=%i, kl=%i\n", m_myRank, 
// 	   ifirst, ilast, jfirst, jlast, kfirst, klast);

    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?

    if( m_point_source_test )
    {
       radius = 4*h;
       x0 = a_globalSources[0]->getX0();
       y0 = a_globalSources[0]->getY0();
       z0 = a_globalSources[0]->getZ0();
    }

// need to exclude parallel overlap from L2 calculation
    if( topographyExists() && g == mNumberOfGrids-1 )
       F77_FUNC(solerr3c, SOLERR3C)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
				     uex_ptr, u_ptr, mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), mJ.c_ptr(),
				     &linfLocal, &l2Local, &xInfGrid, &x0, &y0, &z0, &radius,
				   &imin, &imax, &jmin, &jmax, &kmin, &kmax );
    else
       F77_FUNC(solerr3, SOLERR3)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &h,
				   uex_ptr, u_ptr, &linfLocal, &l2Local, &xInfGrid, &m_zmin[g], &x0,
				   &y0, &z0, &radius,
				   &imin, &imax, &jmin, &jmax, &kmin, &kmax );
    if (linfLocal > diffInfLocal) diffInfLocal = linfLocal;
    if (xInfGrid > xInfLocal) xInfLocal = xInfGrid;
    diffL2Local += l2Local;
  }
// communicate local results for global errors
  MPI_Allreduce( &diffInfLocal, &diffInf, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
  MPI_Allreduce( &xInfLocal,    &xInf,    1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
  MPI_Allreduce( &diffL2Local,  &diffL2,  1, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

  diffL2 = sqrt(diffL2);
}

//---------------------------------------------------------------------------
void EW::normOfDifferenceGhostPoints( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, double &diffL2 )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *uex_ptr, *u_ptr, h, linfLocal=0, l2Local=0, diffInfLocal=0, diffL2Local=0;

//tmp  
  if (proc_zero())
    printf("Inside normOfDifferenceGhostPoints\n");
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    uex_ptr  = a_Uex[g].c_ptr();
    u_ptr    = a_U[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];  

    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    
// need to exclude parallel overlap from L2 calculation
    F77_FUNC(solerrgp, SOLERRGP)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &h,
				uex_ptr, u_ptr, &linfLocal, &l2Local);
    if (linfLocal > diffInfLocal) diffInfLocal = linfLocal;
    diffL2Local += l2Local;
  }
// communicate local results for global errors
  MPI_Allreduce( &diffInfLocal, &diffInf, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
  MPI_Allreduce( &diffL2Local,  &diffL2,  1, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

//   diffL2 = diffL2Local;
//   diffInf = diffInfLocal;
    
  diffL2 = sqrt(diffL2);
}

//---------------------------------------------------------------------------
void EW::normOfSurfaceDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, 
				  double &diffL2, double &solInf, double &solL2, vector<Source*> & a_globalSources)
{
  int g;
  double absDiff, absSol;
  double *uex_ptr, *u_ptr, h, diffInfLocal=0, diffL2Local=0, solInfLocal=0, solL2Local=0;

  g = mNumberOfCartesianGrids-1;
  int k = 1;
  
  h = mGridSize[g];

// only evaluate error on the surface, not including ghost or parallel overlap points

// also exclude points in the super grid damping layer
  int imin, imax, jmin, jmax;
  
  if (mbcGlobalType[0] == bSuperGrid)
    imin = max(m_iStartInt[g], m_sg_gp_thickness+1);
  else
    imin = m_iStartInt[g];
  
  if (mbcGlobalType[1] == bSuperGrid)
    imax = min(m_iEndInt[g], m_global_nx[g] - m_sg_gp_thickness);
  else
    imax = m_iEndInt[g];

  if (mbcGlobalType[2] == bSuperGrid)
    jmin = max(m_jStartInt[g], m_sg_gp_thickness+1);
  else
    jmin = m_jStartInt[g];

  if (mbcGlobalType[3] == bSuperGrid)
    jmax = min(m_jEndInt[g], m_global_ny[g] - m_sg_gp_thickness);
  else
    jmax = m_jEndInt[g];
  
// also need to exclude grid points near the point source
  h = mGridSize[g];

  double radius2, x0, y0, dist2;
  
  if( m_lamb_test )
  {
    radius2 = SQR(4*h);
    x0 = a_globalSources[0]->getX0();
    y0 = a_globalSources[0]->getY0();
  }
  else
  {
    radius2 = -1;
    x0 = 0;
    y0 = 0;
  }
  

  for (int j=jmin; j<=jmax; j++)
    for (int i=imin; i<=imax; i++)
    {
      dist2 = SQR((i-1)*h-x0)+ SQR((j-1)*h-y0);
      
      if( dist2 > radius2 )
      {
	absDiff = fabs(a_Uex[g](3,i,j,k) - a_U[g](3,i,j,k));
	if (absDiff > diffInfLocal) diffInfLocal = absDiff;
	diffL2Local += h*h*absDiff*absDiff;
// exact sol norm
	absSol = fabs(a_Uex[g](3,i,j,k));
	if (absSol > solInfLocal) solInfLocal = absSol;
	solL2Local += h*h*absSol*absSol;
      }
    }
  
// communicate local results for global errors
  MPI_Allreduce( &diffInfLocal, &diffInf, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
  MPI_Allreduce( &diffL2Local,  &diffL2,  1, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

  MPI_Allreduce( &solInfLocal, &solInf, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
  MPI_Allreduce( &solL2Local,  &solL2,  1, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );

  diffL2 = sqrt(diffL2);
  solL2 = sqrt(solL2);
}

//---------------------------------------------------------------------------
void EW::bndryInteriorDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, 
				  double* lowZ, double* interiorZ, double* highZ )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nz;
  double *uex_ptr, *u_ptr, h, li, l2;
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    uex_ptr = a_Uex[g].c_ptr();
    u_ptr   = a_U[g].c_ptr();
    ifirst  = m_iStart[g];
    ilast   = m_iEnd[g];
    jfirst  = m_jStart[g];
    jlast   = m_jEnd[g];
    kfirst  = m_kStart[g];
    klast   = m_kEnd[g];
    h       = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    nz      = m_global_nz[g];
    
// need to do a gather over all processors
    F77_FUNC(rhserrfort, RHSERRFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz, &h,
				      uex_ptr, u_ptr, &lowZ[3*g], &interiorZ[3*g], &highZ[3*g] );
  }
}

//---------------------------------------------------------------------------
void EW::test_RhoUtt_Lu( vector<Sarray> & a_Uacc,  vector<Sarray> & a_Lu,   vector<Sarray> & a_F, 
			 double* lowZ, double* interiorZ, double* highZ )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nz;
  double *rho_ptr, *uacc_ptr, *lu_ptr, *f_ptr, h, li, l2;
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    rho_ptr = mRho[g].c_ptr();
    uacc_ptr= a_Uacc[g].c_ptr();
    lu_ptr  = a_Lu[g].c_ptr();
    f_ptr   = a_F[g].c_ptr();
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    nz = m_global_nz[g];
    
// evaluate rho*uacc - lu - f in fortran routine
     //  subroutine rhouttlumf(ifirst, ilast, jfirst, jlast, kfirst, klast,
     // +     nz, uacc, lu, fo, rho)
    F77_FUNC(rhouttlumf, RHOUTTLUMF)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
				      &nz, uacc_ptr, lu_ptr, f_ptr, rho_ptr,
				      &lowZ[3*g], &interiorZ[3*g], &highZ[3*g]);
  }
}


//---------------------------------------------------------------------------
void EW::initialData(double a_t, vector<Sarray> & a_U, vector<Sarray*> & a_AlphaVE)
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *u_ptr, om, ph, cv, h, zmin;
  
  if (m_twilight_forcing)
  {
     for(int g=0 ; g<mNumberOfCartesianGrids; g++ )
     {
	u_ptr    = a_U[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
	zmin = m_zmin[g];
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	F77_FUNC(twilightfort,TWILIGHTFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					     &klast, u_ptr, &a_t, &om, &cv, &ph, &h, &zmin );
     }
     if( topographyExists() )
     {
	int g = mNumberOfGrids-1;
	u_ptr    = a_U[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	F77_FUNC(twilightfortc,TWILIGHTFORTC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					       &klast, u_ptr, &a_t, &om, &cv, &ph,
					       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
     }
  }
  else if( m_rayleigh_wave_test )
  {
    double cr, lambda, mu, rho, alpha;
    for(int g=0 ; g<mNumberOfCartesianGrids; g++ ) // This case does not make sense with topography
    {
      u_ptr    = a_U[g].c_ptr();
      ifirst = m_iStart[g];
      ilast  = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast  = m_jEnd[g];
      kfirst = m_kStart[g];
      klast  = m_kEnd[g];
      h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
      zmin = m_zmin[g];
      om = m_rayleigh_wave_test->m_omega;
      cr = m_rayleigh_wave_test->m_cr;
      rho = m_rayleigh_wave_test->m_rho;
      lambda = m_rayleigh_wave_test->m_lambda;
      mu = m_rayleigh_wave_test->m_mu;
      alpha = m_rayleigh_wave_test->m_alpha;
      F77_FUNC(rayleighfort,RAYLEIGHFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
					   u_ptr, &a_t, &lambda, &mu, &rho, &cr, &om, &alpha, &h, &zmin );
    }
  }
  else if( m_energy_test )
  {
     //    for(int g=0 ; g<mNumberOfCartesianGrids; g++ ) // This case does not make sense with topography
    for(int g=0 ; g<mNumberOfGrids; g++ ) // This case does not make sense with topography
    {
       u_ptr    = a_U[g].c_ptr();
       for( size_t i=0 ; i < 3*static_cast<size_t>((m_iEnd[g]-m_iStart[g]+1))*(m_jEnd[g]-m_jStart[g]+1)*(m_kEnd[g]-m_kStart[g]+1); i++ )
	  u_ptr[i] = drand48();
    }
  }
  else
// homogeneous initial data is the default
    for(int g=0 ; g<mNumberOfGrids; g++ )
    {
      a_U[g].set_to_zero();
      for( int a=0 ; a < m_number_mechanisms ; a++ )
	a_AlphaVE[g][a].set_to_zero();
    }
}

//---------------------------------------------------------------------------
bool EW::exactSol(double a_t, vector<Sarray> & a_U, vector<Sarray*> & a_AlphaVE,
		  vector<Source*>& sources )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *u_ptr, om, ph, cv, h, zmin;
  bool retval;
  
  if (m_twilight_forcing)
  {
     for(int g=0 ; g<mNumberOfCartesianGrids; g++ ) // curvilinear case needs to be implemented
     {
	u_ptr    = a_U[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
	zmin = m_zmin[g];
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	F77_FUNC(twilightfort,TWILIGHTFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					     &klast, u_ptr, &a_t, &om, &cv, &ph, &h, &zmin );
     }
     if( topographyExists() )
     {
        int g = mNumberOfGrids-1;
	u_ptr    = a_U[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	F77_FUNC(twilightfortc,TWILIGHTFORTC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					       &klast, u_ptr, &a_t, &om, &cv, &ph, 
					       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
     }
     retval = true;
  }
  else if( m_point_source_test )
  {
    for(int g=0 ; g < mNumberOfGrids; g++ ) 
       get_exact_point_source( a_U[g].c_ptr(), a_t, g, *sources[0] );
     retval = true;
  }
  else if( m_lamb_test )
  {
    get_exact_lamb2( a_U, a_t, *sources[0] );
    retval = true;
  }
  else if( m_rayleigh_wave_test ) 
  {
    double cr, lambda, mu, rho, alpha;
    for(int g=0 ; g<mNumberOfCartesianGrids; g++ ) // This case does not make sense with topography
    {
      u_ptr    = a_U[g].c_ptr();
      ifirst = m_iStart[g];
      ilast  = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast  = m_jEnd[g];
      kfirst = m_kStart[g];
      klast  = m_kEnd[g];
      h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
      zmin = m_zmin[g];
      om = m_rayleigh_wave_test->m_omega;
      cr = m_rayleigh_wave_test->m_cr;
      rho = m_rayleigh_wave_test->m_rho;
      lambda = m_rayleigh_wave_test->m_lambda;
      mu = m_rayleigh_wave_test->m_mu;
      alpha = m_rayleigh_wave_test->m_alpha;
      F77_FUNC(rayleighfort,RAYLEIGHFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
					   u_ptr, &a_t, &lambda, &mu, &rho, &cr, &om, &alpha, &h, &zmin );
    }
    
    retval = true;
  }
  else // In general, the exact solution is unknown (m_energy_test falls into this category)
  {
     
     retval = false;
  }
  return retval;
}

//-----------------------------------------------------------------------
// smooth wave for time dependence to test point force term with 
double EW::SmoothWave(double t, double R, double c)
{
  double temp = R;
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1, (c0*pow(t-R/c,3)+c1*pow(t-R/c,4)+c2*pow(t-R/c,5)+c3*pow(t-R/c,6)+c4*pow(t-R/c,7)), 0);
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     temp = (c0*pow(t-R/c,3)+c1*pow(t-R/c,4)+c2*pow(t-R/c,5)+c3*pow(t-R/c,6)+c4*pow(t-R/c,7));
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// very smooth bump for time dependence for further testing of point force 
double EW::VerySmoothBump(double t, double R, double c)
{
  double temp = R;
  double c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120, c5 = -1024;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1, (c0*pow(t-R/c,5)+c1*pow(t-R/c,6)+c2*pow(t-R/c,7)+c3*pow(t-R/c,8)+c4*pow(t-R/c,9)+c5*pow(t-R/c,10)), 0);
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     temp = (c0*pow(t-R/c,5)+c1*pow(t-R/c,6)+c2*pow(t-R/c,7)+c3*pow(t-R/c,8)+c4*pow(t-R/c,9)+c5*pow(t-R/c,10));
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// C6 smooth bump for time dependence for further testing of point force 
double EW::C6SmoothBump(double t, double R, double c)
{
  double retval = 0;
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     retval = 51480.0*pow( (t-R/c)*(1-t+R/c), 7 );
  return retval;
}

//-----------------------------------------------------------------------
// derivative of smooth wave 
double EW::d_SmoothWave_dt(double t, double R, double c)
{
  double temp = R;
  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1, (3*c0*pow(t-R/c,2)+4*c1*pow(t-R/c,3)+5*c2*pow(t-R/c,4)+6*c3*pow(t-R/c,5)+7*c4*pow(t-R/c,6)), 0);
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     temp = (3*c0*pow(t-R/c,2)+4*c1*pow(t-R/c,3)+5*c2*pow(t-R/c,4)+6*c3*pow(t-R/c,5)+7*c4*pow(t-R/c,6));
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// very smooth bump for time dependence to further testing of point force 
double EW::d_VerySmoothBump_dt(double t, double R, double c)
{
  double temp = R;
  double c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120, c5 = -1024;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1, (5*c0*pow(t-R/c,4)+6*c1*pow(t-R/c,5)+7*c2*pow(t-R/c,6)+8*c3*pow(t-R/c,7)+9*c4*pow(t-R/c,8))+10*c5*pow(t-R/c,9), 0);
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     temp = (5*c0*pow(t-R/c,4)+6*c1*pow(t-R/c,5)+7*c2*pow(t-R/c,6)+8*c3*pow(t-R/c,7)+9*c4*pow(t-R/c,8))+10*c5*pow(t-R/c,9);
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// C6 smooth bump for time dependence to further testing of point force 
double EW::d_C6SmoothBump_dt(double t, double R, double c)
{
  double retval=0;
  if( (t-R/c) > 0 && (t-R/c) < 1 )
     retval = 51480.0*7*(1-2*(t-R/c))*pow((t-R/c)*(1-t+R/c),6);
  return retval;
}

//-----------------------------------------------------------------------
// Primitive function (for T) of SmoothWave(t-T)*T
double EW::SWTP(double Lim, double t)
{
  double temp = Lim;

  double c0 = 2187./8., c1 = -10935./8., c2 = 19683./8., c3 = -15309./8., c4 = 2187./4.;

  temp = (pow(t,3)*(c0 + c1*t + c2*pow(t,2) + c3*pow(t,3) + c4*pow(t,4))*pow(Lim,2))/2. - 
    (pow(t,2)*(3*c0 + 4*c1*t + 5*c2*pow(t,2) + 6*c3*pow(t,3) + 7*c4*pow(t,4))*pow(Lim,3))/3. + 
    (t*(3*c0 + 6*c1*t + 10*c2*pow(t,2) + 15*c3*pow(t,3) + 21*c4*pow(t,4))*pow(Lim,4))/4. + 
    ((-c0 - 4*c1*t - 10*c2*pow(t,2) - 20*c3*pow(t,3) - 35*c4*pow(t,4))*pow(Lim,5))/5. + 
    ((c1 + 5*c2*t + 15*c3*pow(t,2) + 35*c4*pow(t,3))*pow(Lim,6))/6. + 
    ((-c2 - 6*c3*t - 21*c4*pow(t,2))*pow(Lim,7))/7. + ((c3 + 7*c4*t)*pow(Lim,8))/8. - (c4*pow(Lim,9))/9.;

  return temp;
}

//-----------------------------------------------------------------------
// Primitive function (for T) of VerySmoothBump(t-T)*T
double EW::VSBTP(double Lim, double t)
{
  double temp = Lim;
  double f = 1024., g = -5120., h = 10240., i = -10240., j = 5120., k = -1024.;

  temp = (pow(Lim,11)*(-25200*k*t-2520*j)+2310*k*pow(Lim,12)+(124740*k*pow(t,2)
							  +24948*j*t+2772*i)*pow(Lim,10)+(-369600*k*pow(t,3)-110880*j*pow(t,2)-24640*i*t-3080*h)*pow(Lim,9)+(727650*k*pow(t,4)+291060*j*pow(t,3)+97020*i*pow(t,2)+24255*h*t+3465*g)*pow(Lim,8)+(-997920*k*pow(t,5)-498960*j*pow(t,4)-221760*i*pow(t,3)-83160*h*pow(t,2)-23760*g*t-3960*f)*pow(Lim,7)+(970200*k*pow(t,6)+582120*j*pow(t,5)+323400*i*pow(t,4)+161700*h*pow(t,3)+69300*g*pow(t,2)+23100*f*t)*pow(Lim,6)+(-665280*k*pow(t,7)-465696*j*pow(t,6)-310464*i*pow(t,5)-194040*h*pow(t,4)-110880*g*pow(t,3)-55440*f*pow(t,2))*pow(Lim,5)+
	  (311850*k*pow(t,8)+249480*j*pow(t,7)+194040*i*pow(t,6)+145530*h*pow(t,5)+103950*g*pow(t,4)+69300*f*pow(t,3))*pow(Lim,4)+(-92400*
																   k*pow(t,9)-83160*j*pow(t,8)-73920*i*pow(t,7)-64680*h*pow(t,6)-55440*g*pow(t,5)-46200*f*pow(t,4))*pow(Lim,3)+(13860*k*pow(t,10)+13860*j*pow(t,9)+13860*i*pow(t,8)+13860*h*pow(t,7)+13860*g*pow(t,6)+13860*f*pow(t,5))*pow(Lim,2))/27720.0;

  return temp;
}
//-----------------------------------------------------------------------
// Primitive function (for T) of C6SmoothBump(t-T)*T
double EW::C6SBTP(double Lim, double t)
{
  double x = t-Lim;
  return pow(x,8)*(-3217.5*pow(x,8)+3432.0*(7+t)*pow(x,7)-25740.0*(3+t)*pow(x,6)
		   +27720.0*(5+3*t)*pow(x,5)-150150.0*(t+1)*x*x*x*x +
		   32760.0*(3+5*t)*x*x*x-36036.0*(1+3*t)*x*x+5720.0*(1+7*t)*x-6435.0*t);
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*SmoothWave(t-T)*T from R/alpha to R/beta
double EW::SmoothWave_x_T_Integral(double t, double R, double alpha, double beta)
{
  double temp = R;

  double lowL, hiL;
  
  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t, R / beta, t);
  if( (R / alpha > t - 1 ) )
     lowL = R/alpha;
  else
     lowL = t-1;
  if( R / beta < t )
     hiL = R/beta;
  else
     hiL = t;
  
  //  temp = where (lowL < t && hiL > t - 1, SWTP(hiL, t) - SWTP(lowL, t), 0.0);
  if( lowL < t && hiL > t - 1 )
     temp = SWTP(hiL, t) - SWTP(lowL, t);
  else
     temp = 0;
  
  return temp;
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*VerySmoothBump(t-T)*T from R/alpha to R/beta
double EW::VerySmoothBump_x_T_Integral(double t, double R, double alpha, double beta)
{
  double temp = R;

  double lowL, hiL;
  
  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t, R / beta, t);
  if( R / alpha > t - 1 )
     lowL = R/alpha;
  else
     lowL = t-1;
  if( R / beta < t )
     hiL = R/beta;
  else
     hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t), 0.0);
  if( lowL < t && hiL > t - 1 )
     temp = VSBTP(hiL, t) - VSBTP(lowL, t);
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*C6SmoothBump(t-T)*T from R/alpha to R/beta
double EW::C6SmoothBump_x_T_Integral(double t, double R, double alpha, double beta)
{
  double temp = R;

  double lowL, hiL;
  
  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t, R / beta, t);
  if( R / alpha > t - 1 )
     lowL = R/alpha;
  else
     lowL = t-1;
  if( R / beta < t )
     hiL = R/beta;
  else
     hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t), 0.0);
  if( lowL < t && hiL > t - 1 )
     temp = C6SBTP(hiL, t) - C6SBTP(lowL, t);
  else
     temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
double EW::Gaussian(double t, double R, double c, double f )
{
  double temp = R;
  temp = 1 /(f* sqrt(2*M_PI))*exp(-pow(t-R/c,2) / (2*f*f));
  return temp;
}

//-----------------------------------------------------------------------
double EW::d_Gaussian_dt(double t, double R, double c, double f)
{
  double temp = R;
  temp = 1 /(f* sqrt(2*M_PI))*(-exp(-pow(t-R/c,2)/(2*f*f))*(t-R/c))/pow(f,2);
  return temp;
}

//-----------------------------------------------------------------------
double EW::Gaussian_x_T_Integral(double t, double R, double f, double alpha, double beta)
{
  double temp = R;
  temp = -0.5*t*(erf( (t-R/beta)/(sqrt(2.0)*f))     - erf( (t-R/alpha)/(sqrt(2.0)*f)) ) -
     f/sqrt(2*M_PI)*( exp(-pow(t-R/beta,2)/(2*f*f) ) - exp( -pow(t-R/alpha,2)/(2*f*f) )  ) ;
     //  temp = 1/(f*sqrt(2*M_PI))*( f*f*(-exp(-pow(t-R/beta,2)/(2*f*f))+exp(-pow(t-R/alpha,2)/(2*f*f)) ) +
     //	     t*0.5*sqrt(M_PI*2)*f*( erf((t-R/alpha)/(sqrt(2.0)*f)) - erf((t-R/beta)/(sqrt(2.0)*f)) ) );
  //  temp = 1 /(f*sqrt(2*M_PI))*(f*( (-exp(-pow(t-R / alpha,2)/pow(f,2)) + exp(-pow(t-R / beta,2)/pow(f,2)) )*f + sqrt(M_PI)*t*(-erf((t-R / alpha) / f) + erf(R / beta / f))))/2.;
  return temp;
}

//-----------------------------------------------------------------------
//void EW::get_exact_point_source( Sarray& u, double t, int g, Source& source )
void EW::get_exact_point_source( double* up, double t, int g, Source& source, int* wind )
{
   timeDep tD;
   if(!( source.getName() == "SmoothWave" || source.getName() == "VerySmoothBump" ||
	 source.getName() == "C6SmoothBump" || source.getName()== "Gaussian") )
   {
      cout << "EW::get_exact_point_source: Error, time dependency must be SmoothWave, VerySmoothBump, C6SmoothBump, or Gaussian, not "
	   << source.getName() << endl;
      return;
   }
   else if( source.getName() == "SmoothWave" )
      tD = iSmoothWave;
   else if( source.getName() == "VerySmoothBump" )
      tD = iVerySmoothBump;
   else if( source.getName() == "C6SmoothBump" )
      tD = iC6SmoothBump;
   else
      tD = iGaussian;

   //   u.set_to_zero();
   double alpha = m_point_source_test->m_cp;
   double beta  = m_point_source_test->m_cs;
   double rho   = m_point_source_test->m_rho;

   double x0    = source.getX0();
   double y0    = source.getY0();
   double z0    = source.getZ0();
   double fr=source.getFrequency();
   double time = (t-source.getOffset()) * source.getFrequency();
   if( tD == iGaussian )
   {
      fr = 1/fr;
      time = time*fr;
   }
   bool ismomentsource = source.isMomentSource();
   double fx, fy, fz;
   double mxx, myy, mzz, mxy, mxz, myz, m0;

   if( !ismomentsource )
   {
      source.getForces( fx, fy, fz );
   }
   else
   {
      source.getMoments( mxx, mxy, mxz, myy, myz, mzz );
      //      m0  = source.getAmplitude();
      m0 = 1;
   }
   bool curvilinear = topographyExists() && g == mNumberOfGrids-1;
   //   double* up = u.c_ptr();
   double h   = mGridSize[g];
   double eps = 1e-3*h;
   size_t ind = 0;
   int imax, imin, jmax, jmin, kmax, kmin;
   if( wind == 0 )
   {
      imin = m_iStart[g];
      imax = m_iEnd[g];
      jmin = m_jStart[g];
      jmax = m_jEnd[g];
      kmin = m_kStart[g];
      kmax = m_kEnd[g];
   }
   else
   {
      imin = wind[0];
      imax = wind[1];
      jmin = wind[2];
      jmax = wind[3];
      kmin = wind[4];
      kmax = wind[5];
   }
   // Note: Use of ind, assumes loop is over the domain over which u is defined.
   //   for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
   //      for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
   //	 for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
   for( int k=kmin ; k <= kmax ; k++ )
      for( int j=jmin ; j <= jmax ; j++ )
	 for( int i=imin ; i <= imax ; i++ )
	 {
            double x,y,z;
	    if( curvilinear )
	    {
               x = mX(i,j,k);
	       y = mY(i,j,k);
	       z = mZ(i,j,k);
	    }
	    else
	    {
	       x = (i-1)*h;
	       y = (j-1)*h;
	       z = (k-1)*h + m_zmin[g];
	    }
	    if( !ismomentsource )
	    {
	       double R = sqrt( (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0) );
	       if( R < eps )
		  up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
	       else
	       {
		  double A, B;
		  if (tD == iSmoothWave)
		  {
		     A = ( 1/pow(alpha,2) * SmoothWave(time, fr*R, alpha) - 1/pow(beta,2) * SmoothWave(time, fr*R, beta) +
			   3/pow(fr*R,2) * SmoothWave_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
	  
		     B = ( 1/pow(beta,2) * SmoothWave(time, fr*R, beta) -
			   1/pow(fr*R,2) * SmoothWave_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  else if (tD == iVerySmoothBump)
		  {
		     A = ( 1/pow(alpha,2) * VerySmoothBump(time, fr*R, alpha) - 1/pow(beta,2) * VerySmoothBump(time, fr*R, beta) +
			   3/pow(fr*R,2) * VerySmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * VerySmoothBump(time, fr*R, beta) -
			   1/pow(fr*R,2) * VerySmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  else if (tD == iC6SmoothBump)
		  {
		     A = ( 1/pow(alpha,2) * C6SmoothBump(time, fr*R, alpha) - 1/pow(beta,2) * C6SmoothBump(time, fr*R, beta) +
			   3/pow(fr*R,2) * C6SmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * C6SmoothBump(time, fr*R, beta) -
			   1/pow(fr*R,2) * C6SmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
                  else if( tD == iGaussian )
		  {
		     A = ( 1/pow(alpha,2) * Gaussian(time, R, alpha,fr) - 1/pow(beta,2) * Gaussian(time, R, beta,fr) +
			   3/pow(R,2) * Gaussian_x_T_Integral(time, R, fr, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
		     
		     B = ( 1/pow(beta,2) * Gaussian(time, R, beta,fr) -
			   1/pow(R,2) * Gaussian_x_T_Integral(time, R, fr, alpha, beta) ) / (4*M_PI*rho*R) ;
		  }
		  up[3*ind]   = ( (x - x0)*(x - x0)*fx + (x - x0)*(y - y0)*fy + (x - x0)*(z - z0)*fz )*A + fx*B;
		  up[3*ind+1] = ( (y - y0)*(x - x0)*fx + (y - y0)*(y - y0)*fy + (y - y0)*(z - z0)*fz )*A + fy*B;
		  up[3*ind+2] = ( (z - z0)*(x - x0)*fx + (z - z0)*(y - y0)*fy + (z - z0)*(z - z0)*fz )*A + fz*B;
	       }
	    }
	    else 
	    {
	       up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
	       // Here, ismomentsource == true
	       double R = sqrt( (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0) );
	       if( R < eps )
	       {
		  up[3*ind] = up[3*ind+1] = up[3*ind+2] = 0;
	       }
	       else
	       {
		  double A, B, C, D, E;
		  if (tD == iSmoothWave)
		  {
		     A = SmoothWave(time, R, alpha);
		     B = SmoothWave(time, R, beta);
		     C = SmoothWave_x_T_Integral(time, R, alpha, beta);
		     D = d_SmoothWave_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_SmoothWave_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iVerySmoothBump)
		  {
		     A = VerySmoothBump(time, R, alpha);
		     B = VerySmoothBump(time, R, beta);
		     C = VerySmoothBump_x_T_Integral(time, R, alpha, beta);
		     D = d_VerySmoothBump_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_VerySmoothBump_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iC6SmoothBump)
		  {
		     A = C6SmoothBump(time, R, alpha);
		     B = C6SmoothBump(time, R, beta);
		     C = C6SmoothBump_x_T_Integral(time, R, alpha, beta);
		     D = d_C6SmoothBump_dt(time, R, alpha) / pow(alpha,3) / R;
		     E = d_C6SmoothBump_dt(time, R, beta) / pow(beta,3) / R;
		  }
		  else if (tD == iGaussian)
		  {
		     A = Gaussian(time, R, alpha,fr);
		     B = Gaussian(time, R, beta,fr);
		     C = Gaussian_x_T_Integral(time, R, fr,alpha, beta);
		     D = d_Gaussian_dt(time, R, alpha,fr) / pow(alpha,3) / R;
		     E = d_Gaussian_dt(time, R, beta,fr) / pow(beta,3) / R;
		  }
		  up[3*ind] += 
	// m_xx*G_xx,x
		     + m0*mxx/(4*M_PI*rho)*
		     ( 
		      + 3*(x-x0)*(x-x0)*(x-x0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	 
		      + ( 15*(x-x0)*(x-x0)*(x-x0) / pow(R,7) - 6*(x-x0) / pow(R,5) ) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	 
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      - 3*(x-x0) / pow(R,5) * C

		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind] +=
		     // m_yy*G_xy,y
		     + m0*myy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_zz*G_xz,z
		     + m0*mzz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_xy*G_xy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_xy*G_xx,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind] +=
		     // m_xz*G_xz,x
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind] +=
		     // m_yz*G_xz,y
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind] +=
		     // m_xz*G_xx,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(x-x0)*(x-x0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) * C
	 
		      + (x-x0)*(x-x0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      - 3*(z-z0) / pow(R,5) * C

		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (z-z0)*E
		      );
		  up[3*ind] +=
		     // m_yz*G_yx,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  //------------------------------------------------------------
		  up[3*ind+1] += 
		     // m_xx*G_xy,x
		     m0*mxx/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_yy**G_yy,y
		     + m0*myy/(4*M_PI*rho)*
		     ( 
		      + 3*(y-y0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + ( 15*(y-y0)*(y-y0)*(y-y0) / pow(R,7) - 6*(y-y0) / pow(R,5) ) * C
	 
		      + (y-y0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind+1] += 
		     // m_zz*G_zy,z
		     + m0*mzz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(z-z0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (z-z0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(z-z0)*(z-z0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_xy*G_yy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) * C
	  
		      + (y-y0)*(y-y0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      - 3*(x-x0) / pow(R,5) * C
	  
		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B
	  
		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind+1] += 
		     // m_xz*G_zy,x
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (y-y0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+1] += 
		     // m_xy*G_xy,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (x-x0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	  
		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	  
		      + ( 15*(x-x0)*(y-y0)*(y-y0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_yz*G_zy,y
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (z-z0)*(y-y0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	  
		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	  
		      + ( 15*(z-z0)*(y-y0)*(y-y0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+1] += 
		     // m_xz*G_xy,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (x-x0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	  
		      + 3*(x-x0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+1] += 
		     // m_yz*G_yy,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(y-y0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(y-y0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      + 15*(z-z0)*(y-y0)*(y-y0) / pow(R,7) * C
	  
		      + (y-y0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	  
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	  
		      - 3*(z-z0) / pow(R,5) * C
	  
		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B
	  
		      + 1 / R * (z-z0)*E
		      );
		  //------------------------------------------------------------
		  up[3*ind+2] += 
		     // m_xx*G_zx,x
		     + m0*mxx/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(x-x0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      + ( 15*(x-x0)*(x-x0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_yy*G_zy,y
		     + m0*myy/(4*M_PI*rho)*
		     (
		      + 3*(y-y0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (y-y0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + ( 15*(y-y0)*(y-y0)*(z-z0) / pow(R,7) - 3*(z-z0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_zz**G_zz,z
		     + m0*mzz/(4*M_PI*rho)*
		     ( 
		      + 3*(z-z0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - 2*(z-z0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + ( 15*(z-z0)*(z-z0)*(z-z0) / pow(R,7) - 6*(z-z0) / pow(R,5) ) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      - 1 / pow(R,3) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      - 3*(z-z0) / pow(R,5) * C

		      + (z-z0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (z-z0)*E
		      );
		  up[3*ind+2] += 
		     // m_xy*G_zy,x
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	  
		      + (y-y0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	  
		      + 3*(y-y0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	  
		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+2] += 
		     // m_xz**G_zz,x
		     + m0*mxz/(4*M_PI*rho)*
		     ( 
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))
	 
		      + 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((x-x0)*D - (x-x0)*E)
	 
		      - 1 / pow(R,3) * ((x-x0)*A/pow(alpha,2) - (x-x0)*B/pow(beta,2))

		      - 3*(x-x0) / pow(R,5) * C

		      + (x-x0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (x-x0)*E
		      );
		  up[3*ind+2] += 
		     // m_xy*G_xz,y
		     + m0*mxy/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(y-y0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (x-x0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)

		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      + 15*(x-x0)*(y-y0)*(z-z0) / pow(R,7) * C
		      );
		  up[3*ind+2] += 
		     // m_yz*G_zz,y
		     + m0*myz/(4*M_PI*rho)*
		     ( 
		      + 3*(y-y0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + 3*(z-z0)*(z-z0) / pow(R,5) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))
	 
		      + 15*(y-y0)*(z-z0)*(z-z0) / pow(R,7) * C
	 
		      + (z-z0)*(z-z0) / pow(R,3)* ((y-y0)*D - (y-y0)*E)
	 
		      - 1 / pow(R,3) * ((y-y0)*A/pow(alpha,2) - (y-y0)*B/pow(beta,2))

		      - 3*(y-y0) / pow(R,5) * C

		      + (y-y0) / (pow(R,3)*pow(beta,2)) * B

		      + 1 / R * (y-y0)*E
		      );
		  up[3*ind+2] += 
		     // m_xz*G_xz,z
		     + m0*mxz/(4*M_PI*rho)*
		     (
		      + 3*(x-x0)*(z-z0)*(z-z0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      - (x-x0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))
	 
		      + (x-x0)*(z-z0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)
	 
		      + 3*(x-x0)*(z-z0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))
	 
		      + ( 15*(x-x0)*(z-z0)*(z-z0) / pow(R,7) - 3*(x-x0) / pow(R,5) ) * C
		      );
		  up[3*ind+2] += 
		     // m_yz*G_yz,z
		     + m0*myz/(4*M_PI*rho)*
		     (
		      + 3*(z-z0)*(z-z0)*(y-y0) / pow(R,5) * (A/pow(alpha,2) - B/pow(beta,2))

		      - (y-y0) / pow(R,3) * (A/pow(alpha,2) - B/pow(beta,2))

		      + (z-z0)*(y-y0) / pow(R,3)* ((z-z0)*D - (z-z0)*E)

		      + 3*(z-z0)*(y-y0) / pow(R,5) * ((z-z0)*A/pow(alpha,2) - (z-z0)*B/pow(beta,2))

		      + ( 15*(z-z0)*(z-z0)*(y-y0) / pow(R,7) - 3*(y-y0) / pow(R,5) ) * C
		      );
	       }
	    }
	    ind++;
	 }
}

#include <cmath>
#include <complex>

//-----------------------------------------------------------------------
complex<double> asin(complex<double> z)
{
  complex<double> I(0,1);
  return -I*log(I*z + sqrt(1. - pow(z,2)));
}
 
//-----------------------------------------------------------------------
complex<double> atan(complex<double> z)
{
  complex<double> I(0,1);
  return I/2.*log((I + z)/(I - z));
}
 
//-----------------------------------------------------------------------
complex<double> atan2(complex<double> z, complex<double> w)
{
  complex<double> I(0,1);
  complex<double> Zero(0,0);
  
  if (w == Zero)
    {
      if (z.real() > 0)
        return M_PI/2.;
      else
        return -M_PI/2.;
    }
  else
    {
      complex<double> retval = I/2.*log((I + z/w)/(I - z/w));
      if( retval.real() < 0 && z.real() > 0 )
         retval = retval + M_PI;
      if( retval.real() > 0 && z.real() < 0 )
         retval = retval - M_PI;
      return retval;
      //      return I/2.*log((I + z/w)/(I - z/w));
    }
}

//-----------------------------------------------------------------------
void EW::get_exact_lamb2( vector<Sarray> & a_U, double a_t, Source& a_source )
{
  int g;
// initialize
  for (g=0; g<mNumberOfGrids; g++)
    a_U[g].set_to_zero();
  double x0 = a_source.getX0();
  double y0 = a_source.getY0();
  double z0 = a_source.getZ0();
  
  double fx, fy, fz;
  a_source.getForces( fx, fy, fz );
  double cs  = m_lamb_test->m_cs;
  double mu  = m_lamb_test->m_mu;
  g = mNumberOfCartesianGrids - 1; // top Cartesian grid
  double h = mGridSize[g];
  int ifirst = m_iStart[g];
  int ilast  = m_iEnd[g];
  int jfirst = m_jStart[g];
  int jlast  = m_jEnd[g];
  int kfirst = m_kStart[g];
  int klast  = m_kEnd[g];
  int tfun = 0;
  if( a_source.getTfunc() == iVerySmoothBump )
     tfun = 1;
  else if( a_source.getTfunc() == iC6SmoothBump )
     tfun = 2;
  F77_FUNC(lambexact,LAMBEXACT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
				 a_U[g].c_ptr(), &a_t, &mu, &cs, &x0, &y0, &fz, &h, &tfun );
}

//-----------------------------------------------------------------------
void EW::get_exact_lamb( vector<Sarray> & a_U, double a_t, Source& a_source )
{
  int g;
  
// initialize
  for (g=0; g<mNumberOfGrids; g++)
    a_U[g].set_to_zero();
  
  double x, y, z, uz, h, t=a_t;
  
  double tau ,r;
  double gamma = sqrt(3. + sqrt(3.))/2.;

  double alpha = m_lamb_test->m_cp;
  double beta  = m_lamb_test->m_cs;
  double mu    = m_lamb_test->m_mu;

  double x0 = a_source.getX0();
  double y0 = a_source.getY0();
  double z0 = a_source.getZ0();
  
  double fx, fy, fz;
  a_source.getForces( fx, fy, fz );
     
// Only the z-component of solution on the flat surface (z=0) is known by this routine
  int k = 1; 
  double R;

  g = mNumberOfCartesianGrids - 1; // top Cartesian grid
  h = mGridSize[g];
  z = 0.0;

//loop over all points in the horizontal plane
  for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
    {
      x = (i-1)*h;
      y = (j-1)*h;

      R = sqrt( (x-x0)*(x-x0)+(y-y0)*(y-y0));
      if( R < h )
      {
	uz = 0;
      }
      else
      {
	uz = 0;
	if ( t <= R/alpha )
	{
	  uz = 0.0;
	}
	else
	{
	  tau = t*beta/R;
	  r = R;
	  if (tau > gamma)
	  {
	    uz += G4_Integral(min(max(0.0,tau - gamma),beta/r), tau, r, beta) - G4_Integral(0.0, tau, r, beta);
	  }
	  if (tau > 1 && tau < beta/r+gamma)
	  {
	    uz += G3_Integral(min(tau - 1,beta/r), tau, r, beta) - G3_Integral(max(0.0,tau - gamma), tau, r, beta);
	  }
	  if (tau > 1/sqrt(3.) && tau < beta/r+1)
	  {
	    uz += G2_Integral(min(tau - 1/sqrt(3.),beta/r), tau, r, beta) - G2_Integral(max(tau - 1,0.0), tau, r, beta);
	  }
	  uz *= -fz/(M_PI*M_PI*mu)*alpha*alpha/(beta*beta*beta);
	}
      } // end if R<h
// assign Sarray
      a_U[g](3,i,j,k) = uz;
    } // end for i,j
} // end get_exact_lamb()


//-----------------------------------------------------------------------
double EW::G4_Integral(double T, double t, double r, double beta)
{
  double c0 = 1024., c1 = -5120., c2 = 10240., c3 = -10240., c4 = 5120., c5 = -1024.;
 
  return -(M_PI*(  (c5*pow(r,9)*pow(T,10))/pow(beta,9) 
		 + (c4*pow(r,8)*pow(T,9))/pow(beta,8) 
		 + (c3*pow(r,7)*pow(T,8))/pow(beta,7) 
		 + (c2*pow(r,6)*pow(T,7))/pow(beta,6) 
		 + (c1*pow(r,5)*pow(T,6))/pow(beta,5) 
		 + (c0*pow(r,4)*pow(T,5))/pow(beta,4)
	     ) ) /8.;
}

//-----------------------------------------------------------------------
double EW::G3_Integral(double iT, double it, double ir, double ibeta)
{
  complex<double> T=iT, t=it, r=ir, beta=ibeta;
  complex<double> c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120, c5 = -1024;
  complex<double> gamma = sqrt(3. + sqrt(3.))/2.;
  complex<double> tmp;
 
  tmp = -(M_PI*((c5*pow(r,9)*pow(T,10))/pow(beta,9) + (c4*pow(r,8)*pow(T,9))/pow(beta,8) + (c3*pow(r,7)*pow(T,8))/pow(beta,7) +
        (c2*pow(r,6)*pow(T,7))/pow(beta,6) + (c1*pow(r,5)*pow(T,6))/pow(beta,5) + (c0*pow(r,4)*pow(T,5))/pow(beta,4)))
    /8.;
 
  tmp += (sqrt(5. + 3.*sqrt(3.))*M_PI*pow(r,4)*(-(sqrt(-pow(t,2) + 2.*t*T - pow(T,2) + pow(gamma,2))*
           (10.*c5*pow(r,5)*(114064.*pow(t,8) + 73744.*pow(t,7)*T + 8.*pow(t,6)*(6698.*pow(T,2) + 150373.*pow(gamma,2)) +
                8.*pow(t,5)*(5018.*pow(T,3) + 68871.*T*pow(gamma,2)) +
                2.*pow(t,4)*(15032.*pow(T,4) + 139272.*pow(T,2)*pow(gamma,2) + 961437.*pow(gamma,4)) +
                2.*pow(t,3)*(11000.*pow(T,5) + 68536.*pow(T,3)*pow(gamma,2) + 284361.*T*pow(gamma,4)) +
                pow(t,2)*(15280.*pow(T,6) + 61032.*pow(T,4)*pow(gamma,2) + 170190.*pow(T,2)*pow(gamma,4) +
                   572519.*pow(gamma,6)) + t*(9520.*pow(T,7) + 22200.*pow(T,5)*pow(gamma,2) + 41574.*pow(T,3)*pow(gamma,4) +
                   82841.*T*pow(gamma,6)) + 128.*(35.*pow(T,8) + 40.*pow(T,6)*pow(gamma,2) + 48.*pow(T,4)*pow(gamma,4) +
                   64.*pow(T,2)*pow(gamma,6) + 128.*pow(gamma,8))) +
             3.*beta*(9.*c4*pow(r,4)*(36528.*pow(t,7) + 23088.*pow(t,6)*T + 24.*pow(t,5)*(682.*pow(T,2) + 11989.*pow(gamma,2)) +
                   8.*pow(t,4)*(1486.*pow(T,3) + 15333.*T*pow(gamma,2)) +
                   pow(t,3)*(8528.*pow(T,4) + 56496.*pow(T,2)*pow(gamma,2) + 305934.*pow(gamma,4)) +
                   pow(t,2)*(5840.*pow(T,5) + 24272.*pow(T,3)*pow(gamma,2) + 75798.*T*pow(gamma,4)) +
                   t*(3600.*pow(T,6) + 8632.*pow(T,4)*pow(gamma,2) + 17226.*pow(T,2)*pow(gamma,4) + 45477.*pow(gamma,6)) +
                   35.*(48.*pow(T,7) + 56.*pow(T,5)*pow(gamma,2) + 70.*pow(T,3)*pow(gamma,4) + 105.*T*pow(gamma,6)))
+
                8.*beta*(8.*c3*pow(r,3)*(4356.*pow(t,6) + 2676.*pow(t,5)*T + 12.*pow(t,4)*(153.*pow(T,2) + 2033.*pow(gamma,2)) +
                      4.*pow(t,3)*(319.*pow(T,3) + 2358.*T*pow(gamma,2)) +
                      pow(t,2)*(856.*pow(T,4) + 3786.*pow(T,2)*pow(gamma,2) + 15525.*pow(gamma,4)) +
                      t*(520.*pow(T,5) + 1298.*pow(T,3)*pow(gamma,2) + 2907.*T*pow(gamma,4)) +
                      48.*(5.*pow(T,6) + 6.*pow(T,4)*pow(gamma,2) + 8.*pow(T,2)*pow(gamma,4) + 16.*pow(gamma,6))) +
                   7.*beta*(2.*beta*(25.*c0*beta*(50.*pow(t,3) + 26.*pow(t,2)*T + 14.*t*pow(T,2) + 6.*pow(T,3) + 55.*t*pow(gamma,2) +
                            9.*T*pow(gamma,2)) + 6.*c1*r*
                          (274.*pow(t,4) + 154.*pow(t,3)*T + 94.*pow(t,2)*pow(T,2) + 54.*t*pow(T,3) + 24.*pow(T,4) +
                            607.*pow(t,2)*pow(gamma,2) + 161.*t*T*pow(gamma,2) + 32.*pow(T,2)*pow(gamma,2) + 64.*pow(gamma,4))) +
                      7.*c2*pow(r,2)*(588.*pow(t,5) + 348.*pow(t,4)*T + 40.*pow(T,5) + 50.*pow(T,3)*pow(gamma,2) +
                         75.*T*pow(gamma,4) + 12.*pow(t,3)*(19.*pow(T,2) + 182.*pow(gamma,2)) +
                         4.*pow(t,2)*(37.*pow(T,3) + 183.*T*pow(gamma,2)) +
                         t*(88.*pow(T,4) + 234.*pow(T,2)*pow(gamma,2) + 693.*pow(gamma,4))))))))/40320. -
       ((10.*c5*pow(r,5)*t*(128.*pow(t,8) + 2304.*pow(t,6)*pow(gamma,2) + 6048.*pow(t,4)*pow(gamma,4) +
               3360.*pow(t,2)*pow(gamma,6) + 315.*pow(gamma,8)) +
            beta*(9.*c4*pow(r,4)*(128.*pow(t,8) + 1792.*pow(t,6)*pow(gamma,2) + 3360.*pow(t,4)*pow(gamma,4) +
                  1120.*pow(t,2)*pow(gamma,6) + 35.*pow(gamma,8)) +
               8.*beta*(8.*c3*pow(r,3)*t*(16.*pow(t,6) + 168.*pow(t,4)*pow(gamma,2) + 210.*pow(t,2)*pow(gamma,4) +
                     35.*pow(gamma,6)) + beta*(7.*c2*pow(r,2)*
                      (16.*pow(t,6) + 120.*pow(t,4)*pow(gamma,2) + 90.*pow(t,2)*pow(gamma,4) + 5.*pow(gamma,6)) +
                     2.*beta*(5.*c0*beta*(8.*pow(t,4) + 24.*pow(t,2)*pow(gamma,2) + 3.*pow(gamma,4)) +
                        6.*c1*r*t*(8.*pow(t,4) + 40.*pow(t,2)*pow(gamma,2) + 15.*pow(gamma,4)))))))*
          atan2((t - T),sqrt(-pow(t,2) + 2.*t*T - pow(T,2) + pow(gamma,2))))/128.))/(48.*pow(beta,9));
 
  //  cout << "ArcTan(Arg) = " << atan((t - T)/sqrt(-pow(t,2) + 2.*t*T - pow(T,2) + pow(gamma,2))) << ". Arg = " << (t - T) << "/" << sqrt(-pow(t,2) + 2.*t*T - pow(T,2) + pow(gamma,2)) << endl;
 
  return tmp.real();
}

//-----------------------------------------------------------------------
double EW::G2_Integral(double iT, double it, double ir, double ibeta)
{
  complex<double> T=iT, t=it, r=ir, beta=ibeta;
  complex<double> c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120, c5 = -1024;
  complex<double> gamma = sqrt(3. + sqrt(3.))/2.;
  complex<double> tmp;

  tmp = (-(M_PI*((c5*pow(r,9)*pow(T,10))/pow(beta,9) + (c4*pow(r,8)*pow(T,9))/pow(beta,8) + (c3*pow(r,7)*pow(T,8))/pow(beta,7) + 
        (c2*pow(r,6)*pow(T,7))/pow(beta,6) + (c1*pow(r,5)*pow(T,6))/pow(beta,5) + (c0*pow(r,4)*pow(T,5))/pow(beta,4)))
    /8.)/2.;

  tmp += ((sqrt(5. + 3.*sqrt(3.))*M_PI*pow(r,4)*(-(sqrt(-pow(t,2) + 2.*t*T - pow(T,2) + pow(gamma,2))*
           (10.*c5*pow(r,5)*(114064.*pow(t,8) + 73744.*pow(t,7)*T + 8.*pow(t,6)*(6698.*pow(T,2) + 150373.*pow(gamma,2)) + 
                8.*pow(t,5)*(5018.*pow(T,3) + 68871.*T*pow(gamma,2)) + 
                2.*pow(t,4)*(15032.*pow(T,4) + 139272.*pow(T,2)*pow(gamma,2) + 961437.*pow(gamma,4)) + 
                2.*pow(t,3)*(11000.*pow(T,5) + 68536.*pow(T,3)*pow(gamma,2) + 284361.*T*pow(gamma,4)) + 
                pow(t,2)*(15280.*pow(T,6) + 61032.*pow(T,4)*pow(gamma,2) + 170190.*pow(T,2)*pow(gamma,4) + 
                   572519.*pow(gamma,6)) + t*(9520.*pow(T,7) + 22200.*pow(T,5)*pow(gamma,2) + 41574.*pow(T,3)*pow(gamma,4) + 
                   82841.*T*pow(gamma,6)) + 128.*(35.*pow(T,8) + 40.*pow(T,6)*pow(gamma,2) + 48.*pow(T,4)*pow(gamma,4) + 
                   64.*pow(T,2)*pow(gamma,6) + 128.*pow(gamma,8))) + 
             3.*beta*(9.*c4*pow(r,4)*(36528.*pow(t,7) + 23088.*pow(t,6)*T + 24.*pow(t,5)*(682.*pow(T,2) + 11989.*pow(gamma,2)) + 
                   8.*pow(t,4)*(1486.*pow(T,3) + 15333.*T*pow(gamma,2)) + 
                   pow(t,3)*(8528.*pow(T,4) + 56496.*pow(T,2)*pow(gamma,2) + 305934.*pow(gamma,4)) + 
                   pow(t,2)*(5840.*pow(T,5) + 24272.*pow(T,3)*pow(gamma,2) + 75798.*T*pow(gamma,4)) + 
                   t*(3600.*pow(T,6) + 8632.*pow(T,4)*pow(gamma,2) + 17226.*pow(T,2)*pow(gamma,4) + 45477.*pow(gamma,6)) + 
                   35.*(48.*pow(T,7) + 56.*pow(T,5)*pow(gamma,2) + 70.*pow(T,3)*pow(gamma,4) + 105.*T*pow(gamma,6))) + 
                8.*beta*(8.*c3*pow(r,3)*(4356.*pow(t,6) + 2676.*pow(t,5)*T + 12.*pow(t,4)*(153.*pow(T,2) + 2033.*pow(gamma,2)) + 
                      4.*pow(t,3)*(319.*pow(T,3) + 2358.*T*pow(gamma,2)) + 
                      pow(t,2)*(856.*pow(T,4) + 3786.*pow(T,2)*pow(gamma,2) + 15525.*pow(gamma,4)) + 
                      t*(520.*pow(T,5) + 1298.*pow(T,3)*pow(gamma,2) + 2907.*T*pow(gamma,4)) + 
                      48.*(5.*pow(T,6) + 6.*pow(T,4)*pow(gamma,2) + 8.*pow(T,2)*pow(gamma,4) + 16.*pow(gamma,6))) + 
                   7.*beta*(2.*beta*(25.*c0*beta*(50.*pow(t,3) + 26.*pow(t,2)*T + 14.*t*pow(T,2) + 6.*pow(T,3) + 55.*t*pow(gamma,2) + 
                            9.*T*pow(gamma,2)) + 6.*c1*r*
                          (274.*pow(t,4) + 154.*pow(t,3)*T + 94.*pow(t,2)*pow(T,2) + 54.*t*pow(T,3) + 24.*pow(T,4) + 
                            607.*pow(t,2)*pow(gamma,2) + 161.*t*T*pow(gamma,2) + 32.*pow(T,2)*pow(gamma,2) + 64.*pow(gamma,4))) + 
                      7.*c2*pow(r,2)*(588.*pow(t,5) + 348.*pow(t,4)*T + 40.*pow(T,5) + 50.*pow(T,3)*pow(gamma,2) + 
                         75.*T*pow(gamma,4) + 12.*pow(t,3)*(19.*pow(T,2) + 182.*pow(gamma,2)) + 
                         4.*pow(t,2)*(37.*pow(T,3) + 183.*T*pow(gamma,2)) + 
                         t*(88.*pow(T,4) + 234.*pow(T,2)*pow(gamma,2) + 693.*pow(gamma,4))))))))/40320. - 
       ((10.*c5*pow(r,5)*t*(128.*pow(t,8) + 2304.*pow(t,6)*pow(gamma,2) + 6048.*pow(t,4)*pow(gamma,4) + 
               3360.*pow(t,2)*pow(gamma,6) + 315.*pow(gamma,8)) + 
            beta*(9.*c4*pow(r,4)*(128.*pow(t,8) + 1792.*pow(t,6)*pow(gamma,2) + 3360.*pow(t,4)*pow(gamma,4) + 
                  1120.*pow(t,2)*pow(gamma,6) + 35.*pow(gamma,8)) + 
               8.*beta*(8.*c3*pow(r,3)*t*(16.*pow(t,6) + 168.*pow(t,4)*pow(gamma,2) + 210.*pow(t,2)*pow(gamma,4) + 
                     35.*pow(gamma,6)) + beta*(7.*c2*pow(r,2)*
                      (16.*pow(t,6) + 120.*pow(t,4)*pow(gamma,2) + 90.*pow(t,2)*pow(gamma,4) + 5.*pow(gamma,6)) + 
                     2.*beta*(5.*c0*beta*(8.*pow(t,4) + 24.*pow(t,2)*pow(gamma,2) + 3.*pow(gamma,4)) + 
                        6.*c1*r*t*(8.*pow(t,4) + 40.*pow(t,2)*pow(gamma,2) + 15.*pow(gamma,4)))))))*
          atan2((t - T),sqrt(-pow(t,2) + 2.*t*T - pow(T,2) + pow(gamma,2))))/128.))/(48.*pow(beta,9)))/2.;

    
    tmp += -(sqrt(-5. + 3.*sqrt(3.))*M_PI*pow(r,4)*(sqrt(-0.75 + sqrt(3.)/4. + pow(t - T,2))*
         (640.*pow(-3. + sqrt(3.),4)*c5*pow(r,5) - 
           (pow(-3. + sqrt(3.),3)*pow(r,3)*(10.*c5*pow(r,2)*(572519.*pow(t,2) + 82841.*t*T + 8192.*pow(T,2)) + 
                9.*beta*(9.*c4*r*(15159.*t + 1225.*T) + 16384.*c3*beta)))/64. + 
           (3.*pow(-3. + sqrt(3.),2)*r*(10.*c5*pow(r,4)*
                 (320479.*pow(t,4) + 94787.*pow(t,3)*T + 28365.*pow(t,2)*pow(T,2) + 6929.*t*pow(T,3) + 1024.*pow(T,4))
                 + 3.*beta*(3.*c4*pow(r,3)*(152967.*pow(t,3) + 37899.*pow(t,2)*T + 8613.*t*pow(T,2) + 1225.*pow(T,3)) + 
                   4.*beta*(8.*c3*pow(r,2)*(5175.*pow(t,2) + 969.*t*T + 128.*pow(T,2)) + 7.*beta*(7.*c2*r*(231.*t + 25.*T) + 256.*c1*beta)))
                ))/8. + 2.*(3. - sqrt(3.))*(10.*c5*pow(r,5)*
               (150373.*pow(t,6) + 68871.*pow(t,5)*T + 34818.*pow(t,4)*pow(T,2) + 17134.*pow(t,3)*pow(T,3) + 
                 7629.*pow(t,2)*pow(T,4) + 2775.*t*pow(T,5) + 640.*pow(T,6)) + 
              3.*beta*(9.*c4*pow(r,4)*(35967.*pow(t,5) + 15333.*pow(t,4)*T + 7062.*pow(t,3)*pow(T,2) + 
                    3034.*pow(t,2)*pow(T,3) + 1079.*t*pow(T,4) + 245.*pow(T,5)) + 
                 2.*beta*(8.*c3*pow(r,3)*(12198.*pow(t,4) + 4716.*pow(t,3)*T + 1893.*pow(t,2)*pow(T,2) + 649.*t*pow(T,3) + 
                       144.*pow(T,4)) + 7.*beta*(7.*c2*pow(r,2)*
                        (1092.*pow(t,3) + 366.*pow(t,2)*T + 117.*t*pow(T,2) + 25.*pow(T,3)) + 
                       beta*(6.*c1*r*(607.*pow(t,2) + 161.*t*T + 32.*pow(T,2)) + 25.*c0*(55.*t + 9.*T)*beta))))) + 
           16.*(10.*c5*pow(r,5)*(7129.*pow(t,8) + 4609.*pow(t,7)*T + 3349.*pow(t,6)*pow(T,2) + 
                 2509.*pow(t,5)*pow(T,3) + 1879.*pow(t,4)*pow(T,4) + 1375.*pow(t,3)*pow(T,5) + 
                 955.*pow(t,2)*pow(T,6) + 595.*t*pow(T,7) + 280.*pow(T,8)) + 
              3.*beta*(9.*c4*pow(r,4)*(2283.*pow(t,7) + 1443.*pow(t,6)*T + 1023.*pow(t,5)*pow(T,2) + 
                    743.*pow(t,4)*pow(T,3) + 533.*pow(t,3)*pow(T,4) + 365.*pow(t,2)*pow(T,5) + 225.*t*pow(T,6) + 
                    105.*pow(T,7)) + 2.*beta*(8.*c3*pow(r,3)*
                     (1089.*pow(t,6) + 669.*pow(t,5)*T + 459.*pow(t,4)*pow(T,2) + 319.*pow(t,3)*pow(T,3) + 
                       214.*pow(t,2)*pow(T,4) + 130.*t*pow(T,5) + 60.*pow(T,6)) + 
                    7.*beta*(7.*c2*pow(r,2)*(147.*pow(t,5) + 87.*pow(t,4)*T + 57.*pow(t,3)*pow(T,2) + 
                          37.*pow(t,2)*pow(T,3) + 22.*t*pow(T,4) + 10.*pow(T,5)) + 
                       beta*(6.*c1*r*(137.*pow(t,4) + 77.*pow(t,3)*T + 47.*pow(t,2)*pow(T,2) + 27.*t*pow(T,3) + 
                             12.*pow(T,4)) + 25.*c0*(25.*pow(t,3) + 13.*pow(t,2)*T + 7.*t*pow(T,2) + 3.*pow(T,3))*beta))))))
         + 315.*((315.*pow(-3. + sqrt(3.),4)*pow(r,4)*(10.*c5*r*t + c4*beta))/256. - 
           (35.*pow(-3. + sqrt(3.),3)*pow(r,2)*(120.*c5*pow(r,3)*pow(t,3) + 
                beta*(36.*c4*pow(r,2)*pow(t,2) + beta*(8.*c3*r*t + c2*beta))))/8. - 
           90.*(-2. + sqrt(3.))*(252.*c5*pow(r,5)*pow(t,5) + 
              beta*(126.*c4*pow(r,4)*pow(t,4) + beta*(56.*c3*pow(r,3)*pow(t,3) + 21.*c2*pow(r,2)*pow(t,2)*beta + 
                    6.*c1*r*t*pow(beta,2) + c0*pow(beta,3)))) + 
           128.*pow(t,4)*(10.*c5*pow(r,5)*pow(t,5) + 
              beta*(9.*c4*pow(r,4)*pow(t,4) + beta*(8.*c3*pow(r,3)*pow(t,3) + 7.*c2*pow(r,2)*pow(t,2)*beta + 
                    6.*c1*r*t*pow(beta,2) + 5.*c0*pow(beta,3)))) - 
           48.*(-3. + sqrt(3.))*pow(t,2)*(120.*c5*pow(r,5)*pow(t,5) + 
              beta*(84.*c4*pow(r,4)*pow(t,4) + beta*(56.*c3*pow(r,3)*pow(t,3) + 35.*c2*pow(r,2)*pow(t,2)*beta + 
                    20.*c1*r*t*pow(beta,2) + 10.*c0*pow(beta,3)))))*log(-t + sqrt(-0.75 + sqrt(3.)/4. + pow(t - T,2)) + T)))/
   (3.87072e6*pow(beta,9));
    

    tmp += (M_PI*pow(r,4)*(4.*sqrt(-0.25 + pow(t - T,2))*(10.*c5*pow(r,5)*
           (7300096.*pow(t,8) + 4719616.*pow(t,7)*T + 128.*pow(t,5)*T*(68871. + 20072.*pow(T,2)) + 
             128.*pow(t,6)*(150373. + 26792.*pow(T,2)) + 8.*pow(t,3)*T*(284361. + 274144.*pow(T,2) + 176000.*pow(T,4)) + 
             8.*pow(t,4)*(961437. + 557088.*pow(T,2) + 240512.*pow(T,4)) + 
             t*T*(82841. + 166296.*pow(T,2) + 355200.*pow(T,4) + 609280.*pow(T,6)) + 
             pow(t,2)*(572519. + 680760.*pow(T,2) + 976512.*pow(T,4) + 977920.*pow(T,6)) + 
             4096.*(1. + 2.*pow(T,2) + 6.*pow(T,4) + 20.*pow(T,6) + 70.*pow(T,8))) + 
          3.*beta*(9.*c4*pow(r,4)*(2337792.*pow(t,7) + 1477632.*pow(t,6)*T + 384.*pow(t,5)*(11989. + 2728.*pow(T,2)) + 
                128.*pow(t,4)*T*(15333. + 5944.*pow(T,2)) + 8.*pow(t,2)*T*(37899. + 48544.*pow(T,2) + 46720.*pow(T,4)) + 
                8.*pow(t,3)*(152967. + 112992.*pow(T,2) + 68224.*pow(T,4)) + 
                35.*T*(105. + 280.*pow(T,2) + 896.*pow(T,4) + 3072.*pow(T,6)) + 
                t*(45477. + 68904.*pow(T,2) + 138112.*pow(T,4) + 230400.*pow(T,6))) + 
             32.*beta*(8.*c3*pow(r,3)*(69696.*pow(t,6) + 42816.*pow(t,5)*T + 48.*pow(t,4)*(2033. + 612.*pow(T,2)) + 
                   32.*pow(t,3)*T*(1179. + 638.*pow(T,2)) + t*T*(2907. + 5192.*pow(T,2) + 8320.*pow(T,4)) + 
                   pow(t,2)*(15525. + 15144.*pow(T,2) + 13696.*pow(T,4)) + 
                   192.*(1. + 2.*pow(T,2) + 6.*pow(T,4) + 20.*pow(T,6))) + 
                7.*beta*(7.*c2*pow(r,2)*(9408.*pow(t,5) + 5568.*pow(t,4)*T + 96.*pow(t,3)*(91. + 38.*pow(T,2)) + 
                      16.*pow(t,2)*T*(183. + 148.*pow(T,2)) + 5.*T*(15. + 40.*pow(T,2) + 128.*pow(T,4)) + 
                      t*(693. + 936.*pow(T,2) + 1408.*pow(T,4))) + 
                   8.*beta*(6.*c1*r*(1096.*pow(t,4) + 616.*pow(t,3)*T + t*T*(161. + 216.*pow(T,2)) + 
                         pow(t,2)*(607. + 376.*pow(T,2)) + 16.*(1. + 2.*pow(T,2) + 6.*pow(T,4))) + 
                      25.*c0*(55.*t + 200.*pow(t,3) + 9.*T + 104.*pow(t,2)*T + 56.*t*pow(T,2) + 24.*pow(T,3))*beta))))) + 
       315.*(10.*c5*pow(r,5)*t*(315. + 13440.*pow(t,2) + 96768.*pow(t,4) + 147456.*pow(t,6) + 32768.*pow(t,8)) + 
          beta*(9.*c4*pow(r,4)*(35. + 4480.*pow(t,2) + 53760.*pow(t,4) + 114688.*pow(t,6) + 32768.*pow(t,8)) + 
             32.*beta*(8.*c3*pow(r,3)*t*(35. + 840.*pow(t,2) + 2688.*pow(t,4) + 1024.*pow(t,6)) + 
                beta*(7.*c2*pow(r,2)*(5. + 360.*pow(t,2) + 1920.*pow(t,4) + 1024.*pow(t,6)) + 
                   8.*beta*(6.*c1*r*t*(15. + 160.*pow(t,2) + 128.*pow(t,4)) + 5.*c0*(3. + 96.*pow(t,2) + 128.*pow(t,4))*beta)))))*
       log(-t + sqrt(-0.25 + pow(t - T,2)) + T)))/(3.3030144e8*sqrt(3.)*pow(beta,9));

  return tmp.real();
}

//---------------------------------------------------------------------------
void EW::exactRhsTwilight(double a_t, vector<Sarray> & a_F)
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *f_ptr, om, ph, cv, h, zmin, omm, phm, amprho, ampmu, ampla;
  
  int g;
  
  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
  {
    f_ptr    = a_F[g].c_ptr();
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    zmin = m_zmin[g];
    if (m_twilight_forcing)
    {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
      omm = m_twilight_forcing->m_momega;
      phm = m_twilight_forcing->m_mphase;
      amprho = m_twilight_forcing->m_amprho;
      ampmu = m_twilight_forcing->m_ampmu;
      ampla = m_twilight_forcing->m_amplambda;
    }
    if( usingSupergrid() )
    {
       double omstrx = m_supergrid_taper_x.get_tw_omega();
       double omstry = m_supergrid_taper_y.get_tw_omega();
       double omstrz = m_supergrid_taper_z.get_tw_omega();
       F77_FUNC(exactrhsfortsg,EXACTRHSFORTSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
						&klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
						&amprho, &ampmu, &ampla, &h, &zmin,
						&omstrx, &omstry, &omstrz );
    }
    else
       F77_FUNC(exactrhsfort,EXACTRHSFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					    &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, 
					    &amprho, &ampmu, &ampla, &h, &zmin );
  }
  if( topographyExists() )
  {
     g = mNumberOfGrids-1;
     f_ptr  = a_F[g].c_ptr();
     ifirst = m_iStart[g];
     ilast  = m_iEnd[g];
     jfirst = m_jStart[g];
     jlast  = m_jEnd[g];
     kfirst = m_kStart[g];
     klast  = m_kEnd[g];
     if (m_twilight_forcing)
     {
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	omm = m_twilight_forcing->m_momega;
	phm = m_twilight_forcing->m_mphase;
	amprho = m_twilight_forcing->m_amprho;
	ampmu = m_twilight_forcing->m_ampmu;
	ampla = m_twilight_forcing->m_amplambda;
     }
     //  subroutine exactaccfort( ifirst, ilast, jfirst, jlast, kfirst, 
     // +     klast, utt, t, om, c, ph, h, zmin )
     if( usingSupergrid() )
     {
	double omstrx = m_supergrid_taper_x.get_tw_omega();
	double omstry = m_supergrid_taper_y.get_tw_omega();
	double omstrz = m_supergrid_taper_z.get_tw_omega();
	F77_FUNC(exactrhsfortsgc,EXACTRHSFORTSGC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
						   &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
						   &amprho, &ampmu, &ampla, 
						   mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(), 
						   &omstrx, &omstry, &omstrz );
     }
     else
	F77_FUNC(exactrhsfortc,EXACTRHSFORTC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					       &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
					       &amprho, &ampmu, &ampla,
					       mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
  }
}

//---------------------------------------------------------------------------
void EW::exactAccTwilight(double a_t, vector<Sarray> & a_Uacc)
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *uacc_ptr, om, ph, cv, h, zmin;
  
  int g;
  
  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
  {
    uacc_ptr    = a_Uacc[g].c_ptr();
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    zmin = m_zmin[g];
    if (m_twilight_forcing)
    {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
    }

     //  subroutine exactaccfort( ifirst, ilast, jfirst, jlast, kfirst, 
     // +     klast, utt, t, om, c, ph, h, zmin )
    F77_FUNC(exactaccfort,EXACTACCFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					 &klast, uacc_ptr, &a_t, &om, &cv, &ph,
					 &h, &zmin );
  }
  if( topographyExists() )
  {
     g = mNumberOfGrids-1;
     uacc_ptr    = a_Uacc[g].c_ptr();
     ifirst = m_iStart[g];
     ilast  = m_iEnd[g];
     jfirst = m_jStart[g];
     jlast  = m_jEnd[g];
     kfirst = m_kStart[g];
     klast  = m_kEnd[g];
     if (m_twilight_forcing)
     {
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
     }
     //  subroutine exactaccfort( ifirst, ilast, jfirst, jlast, kfirst, 
     // +     klast, utt, t, om, c, ph, h, zmin )
    F77_FUNC(exactaccfortc,EXACTACCFORTC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					 &klast, uacc_ptr, &a_t, &om, &cv, &ph,
					 mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
  }
}

//---------------------------------------------------------------------------
void EW::Force(double a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *f_ptr, om, ph, cv, h, zmin, omm, phm, amprho, ampmu, ampla;
  
  int g;
  
  if (m_twilight_forcing)
  {
     for(g=0 ; g<mNumberOfCartesianGrids; g++ )
     {
	f_ptr    = a_F[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
	zmin = m_zmin[g];
    
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	omm = m_twilight_forcing->m_momega;
	phm = m_twilight_forcing->m_mphase;
	amprho = m_twilight_forcing->m_amprho;
	ampmu = m_twilight_forcing->m_ampmu;
	ampla = m_twilight_forcing->m_amplambda;
        if( usingSupergrid() )
	{
	   double omstrx = m_supergrid_taper_x.get_tw_omega();
	   double omstry = m_supergrid_taper_y.get_tw_omega();
	   double omstrz = m_supergrid_taper_z.get_tw_omega();
	   F77_FUNC(forcingfortsg,FORCINGFORTSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					      &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
						  &h, &zmin, &omstrx, &omstry, &omstrz );
	}
        else
	   F77_FUNC(forcingfort,FORCINGFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					      &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
					      &h, &zmin );
     }
     if( topographyExists() )
     {
	g = mNumberOfGrids-1;
	f_ptr    = a_F[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	omm = m_twilight_forcing->m_momega;
	phm = m_twilight_forcing->m_mphase;
	amprho = m_twilight_forcing->m_amprho;
	ampmu = m_twilight_forcing->m_ampmu;
	ampla = m_twilight_forcing->m_amplambda;
        if( usingSupergrid() )
	{
	   double omstrx = m_supergrid_taper_x.get_tw_omega();
	   double omstry = m_supergrid_taper_y.get_tw_omega();
	   double omstrz = m_supergrid_taper_z.get_tw_omega();
	   F77_FUNC(forcingfortcsg,FORCINGFORTCSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					      &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
						    &amprho, &ampmu, &ampla,
						    mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(),
						    &omstrx, &omstry, &omstrz );
	}
        else
	   F77_FUNC(forcingfortc,FORCINGFORTC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
						&klast, f_ptr, &a_t, &om, &cv, &ph, &omm, 
						&phm, &amprho, &ampmu, &ampla,
						mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
     }
  }
  else if( m_rayleigh_wave_test )
  {
     for( int g =0 ; g < mNumberOfGrids ; g++ )
	a_F[g].set_to_zero();
  }
  else if( m_energy_test )
  {
     for( int g =0 ; g < mNumberOfGrids ; g++ )
	a_F[g].set_to_zero();
  }
  else 
  {
     // Default: m_point_source_test, m_lamb_test or full seismic case
     for( int g =0 ; g < mNumberOfGrids ; g++ )
	a_F[g].set_to_zero();

     for( int s = 0 ; s < point_sources.size() ; s++ )
     {
	int g = point_sources[s]->m_grid;
        double fxyz[3];
	point_sources[s]->getFxyz(a_t,fxyz);
	a_F[g](1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
	a_F[g](2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
	a_F[g](3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
     }
  }
}

//---------------------------------------------------------------------------
void EW::Force_tt(double a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *f_ptr, om, ph, cv, h, zmin, omm, phm, amprho, ampmu, ampla;
  
  int g;
  
  if (m_twilight_forcing)
  {
     for(g=0 ; g<mNumberOfCartesianGrids; g++ )
     {
	f_ptr    = a_F[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
	zmin = m_zmin[g];
    
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	omm = m_twilight_forcing->m_momega;
	phm = m_twilight_forcing->m_mphase;
	amprho = m_twilight_forcing->m_amprho;
	ampmu = m_twilight_forcing->m_ampmu;
	ampla = m_twilight_forcing->m_amplambda;
        if( usingSupergrid() )
	{
	   double omstrx = m_supergrid_taper_x.get_tw_omega();
	   double omstry = m_supergrid_taper_y.get_tw_omega();
	   double omstrz = m_supergrid_taper_z.get_tw_omega();
	   F77_FUNC(forcingttfortsg,FORCINGTTFORTSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				      &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
			  	      &h, &zmin, &omstrx, &omstry, &omstrz );
	}
	else
	   F77_FUNC(forcingttfort,FORCINGTTFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
						  &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
						  &amprho, &ampmu, &ampla, &h, &zmin );
     }
     if( topographyExists() )
     {
        g = mNumberOfGrids-1;
	f_ptr    = a_F[g].c_ptr();
	ifirst = m_iStart[g];
	ilast  = m_iEnd[g];
	jfirst = m_jStart[g];
	jlast  = m_jEnd[g];
	kfirst = m_kStart[g];
	klast  = m_kEnd[g];
	om = m_twilight_forcing->m_omega;
	ph = m_twilight_forcing->m_phase;
	cv = m_twilight_forcing->m_c;
	omm = m_twilight_forcing->m_momega;
	phm = m_twilight_forcing->m_mphase;
	amprho = m_twilight_forcing->m_amprho;
	ampmu = m_twilight_forcing->m_ampmu;
	ampla = m_twilight_forcing->m_amplambda;
        if( usingSupergrid() )
	{
	   double omstrx = m_supergrid_taper_x.get_tw_omega();
	   double omstry = m_supergrid_taper_y.get_tw_omega();
	   double omstrz = m_supergrid_taper_z.get_tw_omega();
	   F77_FUNC(forcingttfortcsg,FORCINGTTFORTCSG)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
							&klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
							&amprho, &ampmu, &ampla,
							mX.c_ptr(), mY.c_ptr(), mZ.c_ptr(),
							&omstrx, &omstry, &omstrz );
	}
	else
	   F77_FUNC(forcingttfortc,FORCINGTTFORTC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
						    &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
						    &amprho, &ampmu, &ampla,
						    mX.c_ptr(), mY.c_ptr(), mZ.c_ptr() );
     }
  }
  else if( m_rayleigh_wave_test )
  {
     for( int g =0 ; g < mNumberOfGrids ; g++ )
	a_F[g].set_to_zero();
  }
  else if( m_energy_test )
  {
     for( int g =0 ; g < mNumberOfGrids ; g++ )
	a_F[g].set_to_zero();
  }
  else
  {
     // Default: m_point_source_test, m_lamb_test or full seismic case
     for( int g =0 ; g < mNumberOfGrids ; g++ )
	a_F[g].set_to_zero();

     for( int s = 0 ; s < point_sources.size() ; s++ )
     {
	int g = point_sources[s]->m_grid;
        double fxyz[3];
	point_sources[s]->getFxyztt(a_t,fxyz);
	a_F[g](1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[0];
	a_F[g](2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[1];
	a_F[g](3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0) += fxyz[2];
     }
  }
}

//---------------------------------------------------------------------------
// perhaps a better name would be evalLu ??
void EW::evalRHS(vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		 vector<Sarray> & a_Uacc )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *uacc_ptr, *u_ptr, *mu_ptr, *la_ptr, h;
  
  int *onesided_ptr;
  
  int g, nz;
  
  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
  {
    uacc_ptr = a_Uacc[g].c_ptr();
    u_ptr   = a_U[g].c_ptr();
    mu_ptr  = a_Mu[g].c_ptr();
    la_ptr  = a_Lambda[g].c_ptr();
    //    rho_ptr = mRho[g].c_ptr();
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    nz = m_global_nz[g];
    onesided_ptr = m_onesided[g];
    
    if( usingSupergrid() )
       F77_FUNC(rhs4th3fortsgstr,RHS4TH3FORTSGSTR)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, 
  			          &klast, &nz, onesided_ptr, m_acof, m_bope, m_ghcof,
				   uacc_ptr, u_ptr, mu_ptr, la_ptr, &h,
				   m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g] );
    else
       F77_FUNC(rhs4th3fort,RHS4TH3FORT)(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
				      &klast, &nz, onesided_ptr, m_acof, m_bope, m_ghcof,
				      uacc_ptr, u_ptr, mu_ptr, la_ptr, &h );
  }
  if( topographyExists() )
  {
     g = mNumberOfGrids-1;
     uacc_ptr = a_Uacc[g].c_ptr();
     u_ptr    = a_U[g].c_ptr();
     mu_ptr   = a_Mu[g].c_ptr();
     la_ptr   = a_Lambda[g].c_ptr();
     double* met_ptr = mMetric.c_ptr();
     double* jac_ptr = mJ.c_ptr();
     ifirst   = m_iStart[g];
     ilast    = m_iEnd[g];
     jfirst   = m_jStart[g];
     jlast    = m_jEnd[g];
     kfirst   = m_kStart[g];
     klast    = m_kEnd[g];
     onesided_ptr = m_onesided[g];
     if( usingSupergrid() )
	F77_FUNC(curvilinear4sg,CURVILINEAR4SG)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
					    u_ptr, mu_ptr, la_ptr, met_ptr, jac_ptr,
					    uacc_ptr, onesided_ptr, m_acof, m_bope, m_ghcof,
                                            m_sg_str_x[g], m_sg_str_y[g] );
     else
	F77_FUNC(curvilinear4,CURVILINEAR4)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
					    u_ptr, mu_ptr, la_ptr, met_ptr, jac_ptr,
					    uacc_ptr, onesided_ptr, m_acof, m_bope, m_ghcof );

  }
}

//---------------------------------------------------------------------------
void EW::evalPredictor(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		       vector<Sarray>& a_Rho, vector<Sarray> & a_Lu, vector<Sarray> & a_F )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *up_ptr, *u_ptr, *um_ptr, *lu_ptr, *fo_ptr, *rho_ptr, dt2;
  
  int *onesided_ptr;
  
  int g, nz;
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    u_ptr   = a_U[g].c_ptr();
    um_ptr  = a_Um[g].c_ptr();
    lu_ptr  = a_Lu[g].c_ptr();
    fo_ptr  = a_F[g].c_ptr();
    rho_ptr = a_Rho[g].c_ptr();
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    dt2 = mDt*mDt;
    
     //  subroutine predfort(ifirst, ilast, jfirst, jlast, kfirst, klast,
     // +     up, u, um, lu, fo, rho, dt2 )

    F77_FUNC(predfort,PREDFORT)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
				up_ptr, u_ptr, um_ptr, lu_ptr, fo_ptr, rho_ptr, &dt2 );    

  }
}

//---------------------------------------------------------------------------
void EW::evalCorrector(vector<Sarray> & a_Up, vector<Sarray>& a_Rho,
		       vector<Sarray> & a_Lu, vector<Sarray> & a_F )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *up_ptr, *lu_ptr, *fo_ptr, *rho_ptr, dt4;
  
  int g;
  
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    lu_ptr  = a_Lu[g].c_ptr();
    fo_ptr  = a_F[g].c_ptr();
    rho_ptr = a_Rho[g].c_ptr();
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    dt4 = mDt*mDt*mDt*mDt;
    
     //  subroutine corrfort(ifirst, ilast, jfirst, jlast, kfirst, klast,
     // +     up, lu, fo, rho, dt4 )
    F77_FUNC(corrfort,CORRFORT)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
				up_ptr, lu_ptr, fo_ptr, rho_ptr, &dt4 );

  }
}

//---------------------------------------------------------------------------
void EW::evalDpDmInTime(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
			vector<Sarray> & a_Uacc )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *up_ptr, *u_ptr, *um_ptr, *uacc_ptr, dt2i;
  
  int g;
  
  //  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
  for(g=0 ; g<mNumberOfGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    u_ptr   = a_U[g].c_ptr();
    um_ptr  = a_Um[g].c_ptr();
    uacc_ptr  = a_Uacc[g].c_ptr();

    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    dt2i = 1./(mDt*mDt);
    
     //  subroutine dpdmtfort(ifirst, ilast, jfirst, jlast, kfirst, klast,
     // +     up, u, um, u2, dt2i)

    F77_FUNC(dpdmtfort,DPDMTFORT)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, 
				  up_ptr, u_ptr, um_ptr, uacc_ptr, &dt2i );    

  }
}

//-----------------------------------------------------------------------
// side_plane returns the index of the ghost points along side =0,1,2,3,4,5 (low-i, high-i, low-j, high-j, low-k, high-k)
void EW::side_plane( int g, int side, int wind[6], int nGhost )
{
   wind[0] = m_iStart[g];
   wind[1] = m_iEnd[g];
   wind[2] = m_jStart[g];
   wind[3] = m_jEnd[g];
   wind[4] = m_kStart[g];
   wind[5] = m_kEnd[g];
   // wind[0] = m_ib;
   // wind[1] = m_ie;
   // wind[2] = m_jb;
   // wind[3] = m_je;
   // wind[4] = m_kb;
   // wind[5] = m_ke;
   if( side == 0 )
     wind[1] = wind[0] + (nGhost-1);
   else if( side == 1 )
     wind[0] = wind[1] - (nGhost-1);
   else if( side == 2 )
     wind[3] = wind[2] + (nGhost-1);
   else if( side == 3 )
     wind[2] = wind[3] - (nGhost-1);
   else if( side == 4 )
     wind[5] = wind[4] + (nGhost-1);
   else
     wind[4] = wind[5] - (nGhost-1);
}

//-----------------------------------------------------------------------
void EW::update_images( int currentTimeStep, double time, vector<Sarray> & a_Up,
			vector<Sarray>& a_U, vector<Sarray>& a_Um,
			vector<Sarray>& a_Rho, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
			vector<Source*> & a_sources, int dminus )
{
   double maxerr;
   for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
   {
      Image* img = mImageFiles[fIndex];

      if( img->mMode == Image::HMAXDUDT)
      {
	 if( dminus )
	    img->update_maxes_hVelMax( a_Up, a_U, mDt );
	 else
	    img->update_maxes_hVelMax( a_Up, a_Um, 2*mDt );
      }
      if( img->mMode == Image::HMAX )
	    img->update_maxes_hMax( a_Up );

      if( img->mMode == Image::VMAXDUDT)
      {
	 if( dminus )
	    img->update_maxes_vVelMax( a_Up, a_U, mDt );
	 else
	    img->update_maxes_vVelMax( a_Up, a_Um, 2*mDt );
      }
      if( img->mMode == Image::VMAX )
	 img->update_maxes_vMax( a_Up );

      // Center time derivatives around t-dt, i.e., (up-um)/(2*dt), except when dminus
      // is set. Use (up-u)/dt assumed centered at t, when dminus is true.
      int td = 0;
      if( !dminus )
	 td = img->is_time_derivative();

      if (img->timeToWrite(time-td*mDt , currentTimeStep-td, mDt )) 
      {
	 if(img->mMode == Image::UX ) 
	    img->computeImageQuantity(a_Up, 1);
	 else if(img->mMode == Image::UY )
	    img->computeImageQuantity(a_Up, 2);
	 else if(img->mMode == Image::UZ )
	    img->computeImageQuantity(a_Up, 3);
	 //         else if(img->mMode == Image::FX ) 
	 //            img->computeImageQuantity(mF, 1);
	 //         else if(img->mMode == Image::FY )
	 //            img->computeImageQuantity(mF, 2);
	 //         else if(img->mMode == Image::FZ )
	 //            img->computeImageQuantity(mF, 3);
	 else if(img->mMode == Image::RHO )
	    img->computeImageQuantity(a_Rho, 1);
	 else if(img->mMode == Image::MU )
	    img->computeImageQuantity(a_Mu, 1);
	 else if(img->mMode == Image::LAMBDA )
	    img->computeImageQuantity(a_Lambda, 1);
	 else if(img->mMode == Image::P )
	    img->computeImagePvel(a_Mu, a_Lambda, a_Rho);
	 else if(img->mMode == Image::S )
	    img->computeImageSvel(a_Mu, a_Rho);
	 else if(img->mMode == Image::DIV || img->mMode == Image::DIVDT 
		 || img->mMode == Image::CURLMAG || img->mMode == Image::CURLMAGDT )
	    img->computeImageDivCurl( a_Up, a_U, a_Um, mDt, dminus );
	 else if(img->mMode == Image::LAT || img->mMode == Image::LON )
	    img->computeImageLatLon( mX, mY, mZ );
	 else if(img->mMode == Image::TOPO )
	 {
	    if (topographyExists())
	       img->copy2DArrayToImage(mTopo); // save the raw topography; the smoothed is saved by the mode=grid with z=0
	 }
	 else if( img->mMode == Image::UZEXACT || img->mMode == Image::UXEXACT ||
		  img->mMode == Image::UYEXACT || img->mMode == Image::UXERR   ||
		  img->mMode == Image::UYERR   || img->mMode == Image::UZERR   )
	 {
	    // Note: this is inefficient, the exact solution is computed everywhere, and once for each
	    //   EXACT or ERR image mode.
	    vector<Sarray> Uex(mNumberOfGrids);
	    vector<Sarray*> alpha; //dummy, the array is not used in routine exactSol.
	    for( int g=0 ; g < mNumberOfGrids ; g++ )
	       Uex[g].define(3,m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
	    exactSol( time, Uex, alpha, a_sources );
            if( img->mMode == Image::UXERR || img->mMode == Image::UYERR || img->mMode == Image::UZERR )
	    {
	       for( int g=0 ; g < mNumberOfGrids ; g++ )
	       {
		  size_t n=static_cast<size_t>(Uex[g].npts());
                  int nc  = Uex[g].ncomp();
                  double* uxp = Uex[g].c_ptr();
		  double* up  = a_Up[g].c_ptr();
		  for( size_t i=0 ; i < n*nc ; i++ )
		     uxp[i] = up[i]-uxp[i];
	       }
	    }
	    if( img->mMode == Image::UXEXACT || img->mMode == Image::UXERR )
	       img->computeImageQuantity(Uex,1);
	    else if( img->mMode == Image::UYEXACT || img->mMode == Image::UYERR )
	       img->computeImageQuantity(Uex,2);
	    else if( img->mMode == Image::UZEXACT || img->mMode == Image::UZERR )
	       img->computeImageQuantity(Uex,3);
	    Uex.clear();
	 }
         else if( img->mMode == Image::GRIDX || img->mMode == Image::GRIDY || img->mMode == Image::GRIDZ )
	    img->computeImageGrid(mX, mY, mZ );
         else if(img->mMode == Image::MAGDUDT )
	 {
            if( dminus )
	       img->computeImageMagdt( a_Up, a_U, mDt );
	    else
	       img->computeImageMagdt( a_Up, a_Um, 2*mDt );
	 }
	 else if(img->mMode == Image::HMAGDUDT )
	 {
            if( dminus )
	       img->computeImageHmagdt( a_Up, a_U, mDt );
	    else
	       img->computeImageHmagdt( a_Up, a_Um, 2*mDt );
	 }
         else if(img->mMode == Image::MAG )
	    img->computeImageMag( a_Up );
         else if(img->mMode == Image::HMAG )
	    img->computeImageHmag( a_Up );
	 //         else if(img->mMode == Image::QS )
	 //	 { 
	 //	    if (usingAttenuation())
	 //	       img->computeImageQuantity(mQs, 1);
	 //       }
	 //       else if(img->mMode == Image::QP )
	 //       { 
	 // 	if (usingAttenuation())
	 // 	  img->computeImageQuantity(mQp, 2);
	 //       }
// 	img->computeDivergence();
// 	if (m_forcing->knows_exact())
//         {
//           maxerr = img->computeImageErrorDebug(0);
//           if (proc_zero())
//             printf("maxErr DIV %f @ %fs\n",maxerr,time );
//         }
//       }
//       else if(img->mMode == Image::CURL )
//       {
//         img->computeCurlMagnitude();
// 	if (m_forcing->knows_exact())
//         {
//           maxerr = img->computeImageErrorDebug(2);
//           if (proc_zero())
//             printf("maxErr CURL %f @ %f s\n",maxerr,time );
//         }
//       }
//       else if(img->mMode == Image::VELDIV )
//       {
// 	img->computeVelDivergence();
// 	if (m_forcing->knows_exact())
//         {
//           maxerr = img->computeImageErrorDebug(1);
//           if (proc_zero())
//             printf("maxErr VELDIV %f @ %f s\n",maxerr,time );
//         }
//       }
//       else if(img->mMode == Image::VELCURL )
//       {
//         img->computeVelCurlMagnitude();
// 	if (m_forcing->knows_exact())
//         {
//           maxerr = img->computeImageErrorDebug(3);
//           if (proc_zero())
//             printf("maxErr VELCURL %f @ %f s\n",maxerr,time );
//         }
//       }
//       else if(img->mMode == Image::VELMAG )
//       {
//         img->computeVelocityMagnitude();
//       }
//       else if(img->mMode == Image::HVEL )
//       {
//         img->computeHorizontalVelocityMagnitude();
//       }
// // the following two modes converts mu,lambda,rho into vp and vs
//       else if(img->mMode == Image::S )
//       {
// 	img->computeImageS(mMu, mLambda, mRho);
//       }
//       else if(img->mMode == Image::TOPO )
//       {
//         if (topographyExists())
// 	  img->copy2DArrayToImage(mTopo); // save the raw topography; the smoothed is saved by the mode=grid with z=0
//       }
//       else if(img->mMode == Image::LAT )
//       {
// 	img->evaluateLatLonImage(mX, mY, mZ, 1);
//       }
//       else if(img->mMode == Image::LON )
//       {
// 	img->evaluateLatLonImage(mX, mY, mZ, 2);
//       }
//       else if(img->mMode == Image::GRID )
//       {
// // get the base name
// 	string filePrefix = img->mFilePrefix;
// 	if (img->getOrientation() == Image::X) // save y and z-coordinates
// 	{
// // copy the y-component
// 	  img->evaluateGridImage(mX, mY, mZ, 2); // save the y-component 
// // append a "Y" to the file name
// 	  img->mFilePrefix = filePrefix + "Y"; // will the "Y" accumulate?
// // save the file
//        img->writeImagePlane_2(currentTimeStep, mPath); // save the grid image
// 	}
// 	else if (img->getOrientation() == Image::Y) // save x and z-coordinates
// 	{
// // copy the x-component
// 	  img->evaluateGridImage(mX, mY, mZ, 1); // save the y-component 
// // append a "X" to the file name
// 	  img->mFilePrefix = filePrefix + "X";
// // save the file
//        img->writeImagePlane_2(currentTimeStep, mPath); // save the grid image
// 	}
// // copy the z-component
// 	img->evaluateGridImage(mX, mY, mZ, 3); // save the z-component 
// // append a "Z" to the file name
// 	img->mFilePrefix = filePrefix + "Z";
// // the file is saved below

//       }
//      else
      else if (!img->mMode == Image::HMAXDUDT || !img->mMode == Image::VMAXDUDT
	      || !img->mMode == Image::HMAX   || !img->mMode == Image::VMAX )
      {
	if (proc_zero())
	{
//	  printf("Can only write ux, uy, uz, mu, rho, lambda, uxerr, uyerr, uzerr- remove once completely implemented\n");
	  printf("Can only write ux, uy, uz, mu, rho, lambda: - remove once completely implemented\n");
	  printf("I can not print data of type %i\n", img->mMode );
	}
	MPI_Abort(MPI_COMM_WORLD,1);
      }

// write the image plane on file    
      double t[3];
      t[0] = t[1] = t[2] = MPI_Wtime();

      img->writeImagePlane_2( currentTimeStep-td, mPath, time-td*mDt );
      t[2] = MPI_Wtime();
      
// output timing info?
      if (m_iotiming)
      {
	t[0] = t[1]-t[0];
	t[1] = t[2]-t[1];

	double tmp[2];
	tmp[0] = t[0];
	tmp[1] = t[1];
	MPI_Reduce( tmp, t, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	if( proc_zero() )
	{
	  cout << "Maximum write time:";
	  cout << " (using Bjorn's I/O library) " << t[1] << " seconds. (<=" << m_nwriters << " procs writing)";
	  cout << endl;
	} // end if proc_zero
      } // end if iotiming      
	 
    } // end if time to write
  } // end for all images


  // Volume images
  // for (unsigned int fIndex = 0; fIndex < mImage3DFiles.size(); ++fIndex)
  // {
  //   Image3D* img = mImage3DFiles[fIndex];
  //   if(img->timeToWrite(time, currentTimeStep, mDt ) ) 
  //   {
  //     img->compute_image( );
  //     img->write_images( currentTimeStep, mPath );
  //   }
  // }
} // end EW::update_images()

//-----------------------------------------------------------------------
void EW::addImage(Image* i)
{
  mImageFiles.push_back(i);
}

//-----------------------------------------------------------------------
// void EW::addImage3D(Image3D* i)
// {
//   mImage3DFiles.push_back(i);
// }

//-----------------------------------------------------------------------
void EW::initialize_image_files( )
{
   Image::setSteps(mNumberOfTimeSteps);
//   Image3D::setSteps(mNumberOfTimeSteps);
   for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
   {
     mImageFiles[fIndex]->computeGridPtIndex();
     mImageFiles[fIndex]->allocatePlane();
   }
   // for (unsigned int fIndex = 0; fIndex < mImage3DFiles.size(); ++fIndex)
   // {
   //    mImage3DFiles[fIndex]->setup_images( );
   // }
   for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
   {
      mImageFiles[fIndex]->associate_gridfiles( mImageFiles );
   }
}

//-----------------------------------------------------------------------
void EW::set_sg_thickness(int gp_thickness)
{
  m_sg_gp_thickness = gp_thickness;
  if (m_myRank==0)
    cout << "Default Supergrid thickness has been tuned; thickness = " << m_sg_gp_thickness << " grid sizes" << endl;
}

//-----------------------------------------------------------------------
// void EW::set_sg_transition(int gp_trans)
// {
//   m_sg_gp_transition = gp_trans;
//   if (m_myRank==0)
//     cout << "Default Supergrid transition width has been tuned; width = " << m_sg_gp_transition << " grid sizes" << endl;
// }

//-----------------------------------------------------------------------
void EW::set_sg_damping(double damp_coeff)
{
  m_supergrid_damping_coefficient = damp_coeff;
  if (m_myRank==0)
    cout << "Default Supergrid damping coefficient has been tuned; damping coefficient = " << m_supergrid_damping_coefficient << endl;
}

//-----------------------------------------------------------------------
void EW::set_global_bcs(boundaryConditionType bct[6])
{
  for (int i=0; i<6; i++) 
    mbcGlobalType[i] = bct[i]; 
  mbcsSet = true; 

  //  cout << "mbcGlobalType = " << mbcGlobalType[0] << "," << mbcGlobalType[1] << "," << mbcGlobalType[2] << "," << mbcGlobalType[3] << "," << mbcGlobalType[4] << "," << mbcGlobalType[5] << endl;
}

//-----------------------------------------------------------------------
void EW::set_prefilter( FilterType passband, int order, int passes, double fc1, double fc2 )
{
  m_prefilter_sources = true;
// we could build the filter object right here...
  m_filter_ptr    = new Filter( passband, order, passes, fc1, fc2);
  m_filterobs_ptr = new Filter( passband, order, passes, fc1, fc2);
}

//-----------------------------------------------------------------------
void EW::set_threshold_velocities( double vpmin, double vsmin )
{
   m_useVelocityThresholds = true;
   m_vpMin = vpmin;
   m_vsMin = vsmin;
}

//-----------------------------------------------------------------------
void EW::average_speeds( double& cp, double& cs )
{
   cp = 0;
   cs = 0;
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      int istart = m_iStart[g];
      int iend   = m_iEnd[g];
      int jstart = m_jStart[g];
      int jend   = m_jEnd[g];
      int kstart = m_kStart[g];
      int kend   = m_kEnd[g];
      double cpgrid, csgrid, npts;
      if( m_bcType[g][0] == bSuperGrid )
	 istart = istart+m_sg_gp_thickness;
      if( m_bcType[g][1] == bSuperGrid )
	 iend = iend-m_sg_gp_thickness;
      if( m_bcType[g][2] == bSuperGrid )
	 jstart = jstart+m_sg_gp_thickness;
      if( m_bcType[g][3] == bSuperGrid )
	 jend = jend-m_sg_gp_thickness;
      if( m_bcType[g][4] == bSuperGrid )
	 kstart = kstart+m_sg_gp_thickness;
      if( m_bcType[g][5] == bSuperGrid )
	 kend = kend-m_sg_gp_thickness;
      double* mu_ptr     = mMu[g].c_ptr();
      double* lambda_ptr = mLambda[g].c_ptr();
      double* rho_ptr    = mRho[g].c_ptr();
      F77_FUNC(velsum,VELSUM)(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
			      &istart, &iend, &jstart, &jend, &kstart, &kend, 
			      mu_ptr, lambda_ptr, rho_ptr, &cpgrid, &csgrid, &npts );
      double cpgridtmp = cpgrid;
      double csgridtmp = csgrid;
      double nptstmp   = npts;
      MPI_Allreduce( &cpgridtmp, &cpgrid, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( &csgridtmp, &csgrid, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( &nptstmp, &npts, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      cp = cp + cpgrid/npts;
      cs = cs + csgrid/npts;
   }
   cp = cp/mNumberOfGrids;
   cs = cs/mNumberOfGrids;
}

//-----------------------------------------------------------------------
void EW::layered_speeds( vector<double>& cp, vector<double>& z )
{
   // Sample material on a uniform grid in the z direction with N points:
   int N=200;
   double* zv = new double[N+1];
   double* cpv = new double[N+1];
   double h = (m_global_zmax-m_global_zmin)/N;
   for( int i=0 ; i <= N ; i++ )
      zv[i] = m_global_zmin + i*h;
   double x0 = 0.5*((m_iStart[0]-1)*mGridSize[0] + (m_iEnd[0]-1)*mGridSize[0]);
   double y0 = 0.5*((m_jStart[0]-1)*mGridSize[0] + (m_jEnd[0]-1)*mGridSize[0]);
   int i0, j0, k0, g0;

   for( int i=0; i <= N ; i++ )
   {
      computeNearestGridPoint( i0, j0, k0, g0, x0, y0, zv[i] );
      cpv[i] = sqrt( (2*mMu[g0](i0,j0,k0)+mLambda[g0](i0,j0,k0))/mRho[g0](i0,j0,k0));
   }
   //   for( int b=0 ; b < m_mtrlblocks.size() ; b++ )
   //   {
   //      for( int i=0 ; i <= N ; i++ )
   //	 m_mtrlblocks[b]->get_material_pt( x0, y0, zv[i], rho, cs, cpv[i], qs, qp );
   //   }
   cp.push_back(cpv[0]);
   int j = 1;
   int l = 0;
   double tol = 1e-4;
   while( j < N )
   {
      while( fabs(cp[l]-cpv[j])<tol*cp[l] && (j<N) )
	 j++;
      if( j < N )
      {
         cp.push_back(cpv[j]);
	 z.push_back((zv[j]+zv[j-1])/2);
	 l++;
      }
   }
   delete[] cpv;   
   delete[] zv;   
}

//-----------------------------------------------------------------------
void EW::testsourcediff( vector<Source*> GlobalSources, double gradient[11],
			 double hessian[121] )
{
   for( int m=0 ; m < 11 ; m++ )
   {
      gradient[m] = 0;
      for( int j=0 ; j<11; j++ )
	 hessian[m+11*j] = 0;
   }
   vector<GridPointSource*> gpsources;
   GlobalSources[0]->set_grid_point_sources4( this, gpsources );

   vector<Sarray> kappa, eta;
   int ifirst, ilast, jfirst, jlast, kfirst, klast;
   kappa.resize(mNumberOfGrids);
   eta.resize(mNumberOfGrids);
   for( int g = 0; g <mNumberOfGrids; g++ )
   {
      ifirst = m_iStart[g];
      ilast  = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast  = m_jEnd[g];
      kfirst = m_kStart[g];
      klast  = m_kEnd[g];

      kappa[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      eta[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      kappa[g].set_value(1.0);
      eta[g].set_value(1.0);
   }
   //   for( int m=0 ; m < gpsources.size() ; m++ )
   //   cout << "size of sources " << gpsources.size() << endl;
   for( int m = 0 ; m < gpsources.size()-1 ; m++ )
   {
      gpsources[m]->add_to_gradient( kappa, eta, 0.63, mDt, gradient, mGridSize );
      gpsources[m]->add_to_hessian( kappa, eta, 0.63, mDt, hessian, mGridSize );
   }
}

//-----------------------------------------------------------------------
void EW::get_scalefactors( double sf[11] )
{
   for( int i=0 ; i < 11 ; i++ )
      sf[i] = m_scalefactors[i];
}

//-----------------------------------------------------------------------
bool EW::compute_sf(){return m_compute_scalefactors;}

//-----------------------------------------------------------------------
void EW::compute_guess(bool& guesspos, bool& guesst0fr, bool& guessmom,
		       bool& output_seismograms )
{
   guesspos = m_iniguess_pos;
   guesst0fr = m_iniguess_t0fr;
   guessmom = m_iniguess_mom;
   output_seismograms = m_output_initial_seismograms;
}

//-----------------------------------------------------------------------
void EW::get_cgparameters( int& maxit, int& maxrestart, double& tolerance,
			   bool& fletcherreeves, int& stepselection, bool& do_linesearch,
			   int& varcase, bool& testing )
{
   maxit = m_maxit;
   maxrestart = m_maxrestart;
   tolerance = m_tolerance;
   fletcherreeves = m_cgfletcherreeves;
   stepselection = m_cgstepselection;
   varcase = m_cgvarcase;
   do_linesearch = m_do_linesearch;
   testing = m_opt_testing;
}

//-----------------------------------------------------------------------
void EW::compute_energy( double dt, bool write_file, vector<Sarray>& Um,
			 vector<Sarray>& U, vector<Sarray>& Up, int step )
{
// Compute energy
   double energy    = 0;
   for( int g=0; g < mNumberOfGrids ; g++ )
   {
      int istart = m_iStartInt[g];
      int iend   = m_iEndInt[g];
      int jstart = m_jStartInt[g];
      int jend   = m_jEndInt[g];
      int kstart = m_kStartInt[g];
      int kend   = m_kEndInt[g];
      double* up_ptr  = Up[g].c_ptr(); 
      double* u_ptr   = U[g].c_ptr();
      double* um_ptr  = Um[g].c_ptr();
      double* rho_ptr = mRho[g].c_ptr();
      double locenergy;
      int* onesided_ptr = m_onesided[g];
      if( topographyExists() && g == mNumberOfGrids-1 )
	 F77_FUNC(energy4c,ENERGY4C)(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
			        &istart, &iend, &jstart, &jend, &kstart, &kend, onesided_ptr,
				     um_ptr, u_ptr, up_ptr, rho_ptr, mJ.c_ptr(), &locenergy );
      else
	 F77_FUNC(energy4,ENERGY4)(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g], &m_kEnd[g],
			        &istart, &iend, &jstart, &jend, &kstart, &kend, onesided_ptr,
			        um_ptr, u_ptr, up_ptr, rho_ptr, &mGridSize[g], &locenergy );
      energy += locenergy;
   }
   energy /= (dt*dt);
   double energytmp = energy;
   MPI_Allreduce( &energytmp, &energy, 1, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );
   m_energy_test->record_data( energy, step, write_file, m_myRank, mPath );
}

//-----------------------------------------------------------------------
void EW::set_utcref( TimeSeries& ts )
{
   if( m_utc0set )
      ts.set_station_utc( m_utc0 );
}

//-----------------------------------------------------------------------
void EW::print_utc()
{
   if( proc_zero() )
   {
      if( m_utc0set )
	 printf("EW reference UTC is  %02i/%02i/%i:%i:%i:%i.%i\n", m_utc0[1], m_utc0[2], m_utc0[0], m_utc0[3],
		m_utc0[4], m_utc0[5], m_utc0[6] );
      else
	 printf("EW reference UTC is not defined\n");
   }
}

//-----------------------------------------------------------------------
void EW::extractTopographyFromGridFile( string a_topoFileName )
{
   if (proc_zero())
      cout << "***inside extractTopographyFromGridFile***"<< endl;

//----------------------------------------------
// Check user specified file name. Abort if they are not there or not readable
//----------------------------------------------
   CHECK_INPUT(access(a_topoFileName.c_str(), R_OK) == 0,
	      "No read permission on topography grid file: " << a_topoFileName);

   int topLevel = mNumberOfGrids-1;
  
   double x, y;
   double lat, lon, elev;

// 1. read the grid file
   int Nlon, Nlat, i, j, dum;
   Sarray gridElev;
   double *latv, *lonv;
  
   FILE *gridfile = fopen(a_topoFileName.c_str(),"r");
  
   fscanf(gridfile, "%i %i", &Nlon, &Nlat);
   gridElev.define(1,1,Nlon,1,Nlat,1,1);
   latv = new double[Nlat+1];
   lonv = new double[Nlon+1];

   for (j=1; j<=Nlat; j++)
      for (i=1; i<=Nlon; i++)
	 fscanf(gridfile, "%le %le %le", &lonv[i], &latv[j], &gridElev(1,i,j,1));
   fclose(gridfile);
  
   if (proc_zero())
      printf("Nlon=%i Nlat=%i\n", Nlon, Nlat);

   double lonMax=-180.0, lonMin=180.0, latMax=-90.0, latMin=90.0, elevMax=-1e10, elevMin=1e10;
   for (i=1; i<=Nlon; i++)
   {
      if (lonv[i] < lonMin) lonMin=lonv[i];
      if (lonv[i] > lonMax) lonMax=lonv[i];
   }
   for (i=1; i<=Nlat; i++)
   {
      if (latv[i] < latMin) latMin=latv[i];
      if (latv[i] > latMax) latMax=latv[i];
   }
  
   for (i=1; i<=Nlon; i++)
      for (j=1; j<=Nlat; j++)
      {
	 if (gridElev(1,i,j,1) < elevMin) elevMin=gridElev(1,i,j,1);
	 if (gridElev(1,i,j,1) > elevMax) elevMax=gridElev(1,i,j,1);
      }
   if (proc_zero())
      printf("lonMin=%e, lonMax=%e\nlatMin=%e, latMax=%e\nelevMin=%e, evalMax=%e\n", 
	     lonMin, lonMax, latMin, latMax, elevMin, elevMax);
  
// If the lat vector is not in increasing order, we need to reorder it
   if (latv[1] > latv[Nlat])
   {
// tmp
      if (proc_zero()) printf("Reordering the latitude vector...\n");

      for (j=1; j<=Nlat/2; j++)
      {
	 lat=latv[Nlat+1-j];
	 latv[Nlat+1-j] = latv[j];
	 latv[j] = lat;
      
	 for (i=1; i<=Nlon; i++)
	 {
	    elev = gridElev(1,i,Nlat+1-j,1);
	    gridElev(1,i,Nlat+1-j,1) = gridElev(1,i,j,1);
	    gridElev(1,i,j,1) = elev;
	 }
      }// end for j    
   } // end if latv[1] > latv[Nlat]
  
// If the lon vector is not in increasing order, we need to reorder it
   if (lonv[1] > lonv[Nlon])
   {
// tmp
      if (proc_zero()) printf("Reordering the longitude vector...\n");

      for (i=1; i<=Nlon/2; i++)
      {
	 lon=lonv[Nlon+1-i];
	 lonv[Nlon+1-i] = lonv[i];
	 lonv[i] = lon;
      
	 for (j=1; j<=Nlat; j++)
	 {
	    elev = gridElev(1,Nlon+1-i,j,1);
	    gridElev(1,Nlon+1-i,j,1) = gridElev(1,i,j,1);
	    gridElev(1,i,j,1) = elev;
	 }
      }// end for i    
   } // end if lonv[1] > lonv[Nlon]
  

// 2. interpolate in the grid file to get elevations on the computational grid
   double deltaLat = (latMax-latMin)/Nlat;
   double deltaLon = (lonMax-lonMin)/Nlon;
   double eInterp, xi, eta;

   int i0, j0;
  
   for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
   {
      for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
      {
	 x = (i-1)*mGridSize[topLevel];
	 y = (j-1)*mGridSize[topLevel];
        
	 computeGeographicCoord( x, y, lon, lat ); // the grid command defines the parameters in this mapping
      
	 if (lat > latMax || lat < latMin || lon > lonMax || lon < lonMin)
	 {
	    printf("ERROR: x=%e, y=%e corresponds to lon=%e, lat=%e which are outside the topography grid\n", 
		   x, y, lon, lat);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	 }
	 i0 = 1+(int)((lon-lonMin)/deltaLon);
	 j0 = 1+(int)((lat-latMin)/deltaLat);

	 while ( lon < lonv[i0] || lonv[i0+1] < lon ) // should stop loop if i0 is out of bounds
	 {
	    if (lon<lonv[i0]) 
	       i0--;
	    else if (lon>lonv[i0+1])
	       i0++;
	 }

	 while (  lat < latv[j0] || latv[j0+1] < lat ) // should stop loop if j0 is out of bounds
	 {
	    if (lat<latv[j0]) 
	       j0--;
	    else if (lat>latv[j0+1])
	       j0++;
	 }
      
	 if (i0 > Nlon-1) i0 = Nlon-1;

	 if (j0 > Nlat-1) j0 = Nlat-1;
      
// test that we are inside the interval
	 if (!(lonv[i0] <= lon && lon < lonv[i0+1]))
	 {
	    printf("EW::extractTopographyFromGridFile: Fatal error: Unable to interpolate topography for lon=%e\n"
	       "because it is outside the cell (lonv[%i]=%e, lonv[%i]=%e)\n", lon, i0, lonv[i0], i0+1, lonv[i0+1]);
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
      
      
	 if (!(latv[j0] <= lat && lat < latv[j0+1]))
	 {
	    printf("EW::extractTopographyFromGridFile: Fatal error: Unable to interpolate topography for lat=%e\n"
	       "because it is outside the cell (latv[%i]=%e, latv[%i]=%e)\n", lat, j0, latv[j0], j0+1, latv[j0+1]);
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
      
// bi-cubic interpolation should make the surface smoother
// shift the stencil if it is too close to the boundaries
	 if (i0 < 2) i0 = 2;
	 if (i0 > Nlon-2) i0 = Nlon-2;

	 if (j0 < 2) j0 = 2;
	 if (j0 > Nlat-2) j0 = Nlat-2;

// local step sizes
	 double q = i0 + (lon - lonv[i0])/(lonv[i0+1]-lonv[i0]);
	 double r = j0 + (lat - latv[j0])/(latv[j0+1]-latv[j0]);
      
	 double Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1, tjp2;
	 Qim1 = (q-i0)*(q-i0-1)*(q-i0-2)/(-6.);
	 Qi   = (q-i0+1)*(q-i0-1)*(q-i0-2)/(2.);
	 Qip1 = (q-i0+1)*(q-i0)*(q-i0-2)/(-2.);
	 Qip2 = (q-i0+1)*(q-i0)*(q-i0-1)/(6.);

	 Rjm1 = (r-j0)*(r-j0-1)*(r-j0-2)/(-6.);
	 Rj   = (r-j0+1)*(r-j0-1)*(r-j0-2)/(2.);
	 Rjp1 = (r-j0+1)*(r-j0)*(r-j0-2)/(-2.);
	 Rjp2 = (r-j0+1)*(r-j0)*(r-j0-1)/(6.);

	 tjm1 = Qim1*gridElev(i0-1,j0-1,1) + Qi*gridElev(i0,j0-1,1) +  Qip1*gridElev(i0+1,j0-1,1) +  Qip2*gridElev(i0+2,j0-1,1);
	 tj   = Qim1*gridElev(i0-1,j0,1) + Qi*gridElev(i0,j0,1) +  Qip1*gridElev(i0+1,j0,1) +  Qip2*gridElev(i0+2,j0,1);
	 tjp1 = Qim1*gridElev(i0-1,j0+1,1) + Qi*gridElev(i0,j0+1,1) +  Qip1*gridElev(i0+1,j0+1,1) +  Qip2*gridElev(i0+2,j0+1,1);
	 tjp2 = Qim1*gridElev(i0-1,j0+2,1) + Qi*gridElev(i0,j0+2,1) +  Qip1*gridElev(i0+1,j0+2,1) +  Qip2*gridElev(i0+2,j0+2,1);

	 mTopo(i,j,1) = Rjm1*tjm1 + Rj*tj + Rjp1*tjp1 + Rjp2*tjp2;

// bi-linear interpolation
// // local step sizes
//       xi = (lon - lonv[i0])/(lonv[i0+1]-lonv[i0]);
//       eta = (lat - latv[j0])/(latv[j0+1]-latv[j0]);
//       mTopo(i,j,1) = (1.0-eta)*( (1.0-xi)*gridElev(1,i0,j0,1) + xi*gridElev(1,i0+1,j0,1) ) +
// 	eta*( (1.0-xi)*gridElev(1,i0,j0+1,1) + xi*gridElev(1,i0+1,j0+1,1) );
      
      }
   }
   delete[] latv;
   delete[] lonv;
}

//-----------------------------------------------------------------------
void EW::extractTopographyFromCartesianFile(string a_topoFileName)
{
   if (proc_zero())
      cout << "***inside extractTopographyFromCartesianFile***"<< endl;

//----------------------------------------------
// Check user specified file name. Abort if they are not there or not readable
//----------------------------------------------
   VERIFY2(access(a_topoFileName.c_str(), R_OK) == 0,
	       "No read permission on topography grid file: " << a_topoFileName);

// 1. read the grid file
   int Nx, Ny, i, j;
   Sarray gridElev;
   double *yv, *xv;
  
   FILE *gridfile = fopen(a_topoFileName.c_str(),"r");
  
   fscanf(gridfile, "%i %i", &Nx, &Ny);
   gridElev.define(1,1,Nx,1,Ny,1,1);
   yv = new double[Ny+1];
   xv = new double[Nx+1];

   for (j=1; j<=Ny; j++)
      for (i=1; i<=Nx; i++)
	 fscanf(gridfile, "%le %le %le", &xv[i], &yv[j], &gridElev(1,i,j,1));
   fclose(gridfile);
  
   if (proc_zero())
      printf("Nx=%i Ny=%i\n", Nx, Ny);

   double xMax=-999e10, xMin=999e10, yMax=-999e10, yMin=999e10, elevMax=-1e10, elevMin=1e10;
   for (i=1; i<=Nx; i++)
   {
      if (xv[i] < xMin) xMin=xv[i];
      if (xv[i] > xMax) xMax=xv[i];
   }
   for (i=1; i<=Ny; i++)
   {
      if (yv[i] < yMin) yMin=yv[i];
      if (yv[i] > yMax) yMax=yv[i];
   }
// make sure that the topography grid covers the whole computational domain
   if (xMin > 0 || yMin > 0 || xMax < m_global_xmax || yMax < m_global_ymax)
   {
      if (proc_zero()) printf("ERROR: Cartesian topography grid with %e<=x<=%e and %e<=y<=%e\n"
			      "does not cover the computational domain: 0<=x<=%e, 0<=y<=%e\n", 
			    xMin, xMax, yMin, yMax, m_global_xmax, m_global_ymax);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   for (i=1; i<=Nx; i++)
      for (j=1; j<=Ny; j++)
      {
	 if (gridElev(1,i,j,1) < elevMin) elevMin=gridElev(1,i,j,1);
	 if (gridElev(1,i,j,1) > elevMax) elevMax=gridElev(1,i,j,1);
      }
   if (proc_zero())
      printf("xMin=%e, xMax=%e\nyMin=%e, yMax=%e\nelevMin=%e, evalMax=%e\n", 
	     xMin, xMax, yMin, yMax, elevMin, elevMax);
  
   double xP, yP, elev;
// If the yv vector is not in increasing order, we need to reorder it
   if (yv[1] > yv[Ny])
   {
// tmp
      if (proc_zero()) printf("Reordering the yv vector...\n");
      for (j=1; j<=Ny/2; j++)
      {
	 yP=yv[Ny+1-j];
	 yv[Ny+1-j] = yv[j];
	 yv[j] = yP;
      
	 for (i=1; i<=Nx; i++)
	 {
	    elev = gridElev(1,i,Ny+1-j,1);
	    gridElev(1,i,Ny+1-j,1) = gridElev(1,i,j,1);
	    gridElev(1,i,j,1) = elev;
	 }
      }
   }
  
// If the xv vector is not in increasing order, we need to reorder it
   if (xv[1] > xv[Nx])
   {
// tmp
      if (proc_zero()) printf("Reordering the xv vector...\n");
      for (i=1; i<=Nx/2; i++)
      {
	 xP=xv[Nx+1-i];
	 xv[Nx+1-i] = xv[i];
	 xv[i] = xP;
      
	 for (j=1; j<=Ny; j++)
	 {
	    elev = gridElev(1,Nx+1-i,j,1);
	    gridElev(1,Nx+1-i,j,1) = gridElev(1,i,j,1);
	    gridElev(1,i,j,1) = elev;
	 }
      }
   }

// 2. interpolate in the grid file to get elevations on the computational grid
   double deltaY = (yMax-yMin)/Ny;
   double deltaX = (xMax-xMin)/Nx;

   int i0, j0;
   int topLevel = mNumberOfGrids-1;
   double hp = 1.01*mGridSize[topLevel];
   bool xGhost, yGhost;
  
   for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
   {
      for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
      {
	 xP = (i-1)*mGridSize[topLevel];
	 yP = (j-1)*mGridSize[topLevel];
	 xGhost=true;
	 yGhost=true;
	 if (yP > yMax+hp || yP < yMin-hp || xP > xMax+hp || xP < xMin-hp)
	 {
	    printf("ERROR: xP=%e, yP=%e is outside the topography grid by more than a grid step\n", 
	       xP, yP);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	 }
// Compute i0
	 if (xP < xMin)
	    i0=1;
	 else if (xP > xMax)
	    i0=Nx-1;
	 else
	 {
	    xGhost=false;
	    i0 = 1+(int)((xP-xMin)/deltaX);
	    if (i0 < 1)
	       i0 = 1;
	    if (i0 > Nx-1)
	       i0 = Nx-1;
// should stop loop if i0 is out of bounds
	    while ( i0>=1 && i0 <= Nx-1 && ( xP < xv[i0] || xv[i0+1] < xP ) ) 
	    {
	       if(xP<xv[i0]) 
		  i0--;
	       else if (xP>xv[i0+1])
		  i0++;
	    }
	 }
      
// Compute j0
	 if (yP < yMin)
	    j0=1;
	 else if (yP > yMax)
	    j0=Ny-1;
	 else
	 {
	    yGhost=false;
	    j0 = 1+(int)((yP-yMin)/deltaY);
	    if (j0 < 1)
	       j0 = 1;
	    if (j0 > Ny-1)
	       j0 = Ny-1;
 // should stop loop if j0 is out of bounds
	    while ( j0>=1 && j0 <= Ny-1 && ( yP < yv[j0] || yv[j0+1] < yP ) )
	    {
	       if (yP<yv[j0]) 
		  j0--;
	       else if (yP>yv[j0+1])
		  j0++;
	    }
	 }
      
// enforce bounds again
	 if (i0 < 1) i0 = 1;
	 if (i0 > Nx-1) i0 = Nx-1;
	 if (j0 < 1) j0 = 1;
	 if (j0 > Ny-1) j0 = Ny-1;
      
// test that we are inside the interval
	 if (!xGhost && !(xv[i0] <= xP && xP <= xv[i0+1]))
	 {
	    printf("EW::extractTopographyFromCartesianFile: Fatal error: Unable to interpolate topography for xP=%e\n"
	       "because it is outside the cell (xv[%i]=%e, xv[%i]=%e)\n", xP, i0, xv[i0], i0+1, xv[i0+1]);
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
      
      
	 if (!yGhost && !(yv[j0] <= yP && yP <= yv[j0+1]))
	 {
	    printf("EW::extractTopographyFromCartesianFile: Fatal error: Unable to interpolate topography for yP=%e\n"
	       "because it is outside the cell (yv[%i]=%e, yv[%i]=%e)\n", yP, j0, yv[j0], j0+1, yv[j0+1]);
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
      
// bi-cubic interpolation should make the surface smoother
// shift the stencil if it is too close to the boundaries
	 if (i0 < 2) i0 = 2;
	 if (i0 > Nx-2) i0 = Nx-2;

	 if (j0 < 2) j0 = 2;
	 if (j0 > Ny-2) j0 = Ny-2;

// local step sizes
	 double q = i0 + (xP - xv[i0])/(xv[i0+1]-xv[i0]);
	 double r = j0 + (yP - yv[j0])/(yv[j0+1]-yv[j0]);
      
	 double Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1, tjp2;
	 Qim1 = (q-i0)*(q-i0-1)*(q-i0-2)/(-6.);
	 Qi   = (q-i0+1)*(q-i0-1)*(q-i0-2)/(2.);
	 Qip1 = (q-i0+1)*(q-i0)*(q-i0-2)/(-2.);
	 Qip2 = (q-i0+1)*(q-i0)*(q-i0-1)/(6.);

	 Rjm1 = (r-j0)*(r-j0-1)*(r-j0-2)/(-6.);
	 Rj   = (r-j0+1)*(r-j0-1)*(r-j0-2)/(2.);
	 Rjp1 = (r-j0+1)*(r-j0)*(r-j0-2)/(-2.);
	 Rjp2 = (r-j0+1)*(r-j0)*(r-j0-1)/(6.);

	 tjm1 = Qim1*gridElev(i0-1,j0-1,1) +    Qi*gridElev(i0,  j0-1,1)
	     +  Qip1*gridElev(i0+1,j0-1,1) +  Qip2*gridElev(i0+2,j0-1,1);
	 tj   = Qim1*gridElev(i0-1,j0,  1) +    Qi*gridElev(i0,  j0,  1)
	     +  Qip1*gridElev(i0+1,j0,  1) +  Qip2*gridElev(i0+2,j0,  1);
	 tjp1 = Qim1*gridElev(i0-1,j0+1,1) +    Qi*gridElev(i0,  j0+1,1)
	     +  Qip1*gridElev(i0+1,j0+1,1) +  Qip2*gridElev(i0+2,j0+1,1);
	 tjp2 = Qim1*gridElev(i0-1,j0+2,1) +    Qi*gridElev(i0,  j0+2,1)
	     +  Qip1*gridElev(i0+1,j0+2,1) +  Qip2*gridElev(i0+2,j0+2,1);

	 mTopo(i,j,1) = Rjm1*tjm1 + Rj*tj + Rjp1*tjp1 + Rjp2*tjp2;
      }
   }
   delete[] yv;
   delete[] xv;  
}

//-----------------------------------------------------------------------
void EW::extractTopographyFromEfile(std::string a_topoFileName, std::string a_topoExtFileName,
				    std::string a_QueryType, double a_EFileResolution )
{
#ifdef ENABLE_ETREE
   if (proc_zero())
      cout << endl <<
	 "*** extracting TOPOGRAPHY from efile ***"<< endl << endl;
   cencalvm::query::VMQuery query;
   cencalvm::storage::ErrorHandler* pErrHandler = query.errorHandler();

// Check user specified file names. Abort if they are not there or not readable
   VERIFY2( access(a_topoFileName.c_str(), R_OK) == 0,
	    "No read permission on etree file: " << a_topoFileName);
   query.filename(a_topoFileName.c_str());

   if (a_topoExtFileName != "NONE")
   {
      // User specified, if it is not there, abort
      VERIFY2(access(a_topoExtFileName.c_str(), R_OK) == 0,
	      "No read permission on xefile: " << a_topoExtFileName);
      query.filenameExt(a_topoExtFileName.c_str());
   }
   int topLevel = mNumberOfGrids-1;
   if (a_QueryType == "MAXRES")
      query.queryType(cencalvm::query::VMQuery::MAXRES);
   else if (a_QueryType == "FIXEDRES")
   {
      query.queryType(cencalvm::query::VMQuery::FIXEDRES);
      if (a_EFileResolution < 0.)
	 a_EFileResolution = mGridSize[topLevel];
      if (proc_zero())
	 printf("Fixedres resolution = %e\n", a_EFileResolution);
      query.queryRes(a_EFileResolution);
   }

   const char* queryKeys[] = {"elevation", "Vp", "Vs"};
   int payloadSize = 3;
   query.queryVals(queryKeys, payloadSize);

   double x, y;
   double lat, lon, elev, elevDelta=25.;

   query.open();
   double *pVals=new double[payloadSize];

   for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
      for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
      {
	 x = (i-1)*mGridSize[topLevel];
	 y = (j-1)*mGridSize[topLevel];
	 computeGeographicCoord( x, y, lon, lat ); 
// initial query for elevation just below sealevel
	 elev = -25.0;
	 query.query(&pVals, payloadSize, lon, lat, elev);
// Make sure the query didn't generated a warning or error
	 if (pErrHandler->status() != cencalvm::storage::ErrorHandler::OK) 
	 {
// If query generated an error, then bail out, otherwise reset status
	    pErrHandler->resetStatus();
	    cout << "WARNING: Etree query failed for initial elevation of topography at grid point (i,j)= ("
		 << i << ", " << j << ") in curvilinear grid g = " << topLevel << endl
		 << " lat= " << lat << " lon= " << lon << " query elevation= " << elev << endl;
	    mTopo(i,j,1) = 0.;
	    continue;
	 } 
// save the actual topography which will be the starting point for computing smoother the grid topography
	 mTopo(i,j,1) = pVals[0];
      } 
   query.close();
#endif
}

//-----------------------------------------------------------------------
bool EW::is_onesided( int g, int side ) const
{
   return m_onesided[g][side] == 1;
}

//-----------------------------------------------------------------------
void EW::get_gridgen_info( int& order, double& zetaBreak ) const
{
   order = m_grid_interpolation_order;
   zetaBreak = m_zetaBreak;
}

//-----------------------------------------------------------------------
void EW::get_nr_of_material_parameters( int& nmvar )
{
   nmvar = 0;
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      nmvar += (m_iEndAct[g]-m_iStartAct[g]+1)*(m_jEndAct[g]-m_jStartAct[g]+1)*(m_kEndAct[g]-m_kStartAct[g]+1)*3;
   }
}

//-----------------------------------------------------------------------
void EW::parameters_to_material( int nmpar, double* xm, vector<Sarray>& rho,
				 vector<Sarray>& mu, vector<Sarray>& lambda )
{
   int gp, ind;
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      rho[g].copy( mRho[g] );
      mu[g].copy( mMu[g] );
      lambda[g].copy( mLambda[g] );
      if( g == 0 )
	 gp = 0;
      else
	 gp = gp + 3*ind;
      ind =0;
      for( int k=m_kStartAct[g]; k <= m_kEndAct[g]; k++ )
	 for( int j=m_jStartAct[g]; j <= m_jEndAct[g]; j++ )
	    for( int i=m_iStartAct[g]; i <= m_iEndAct[g]; i++ )
	    {
	       rho[g](i,j,k)    = xm[gp+ind*3];
	       mu[g](i,j,k)     = xm[gp+ind*3+1];
	       lambda[g](i,j,k) = xm[gp+ind*3+2];
	       ind++;
	    }
   }
}

//-----------------------------------------------------------------------
void EW::get_material_parameter( int nmpar, double* xm )
{
   int gp, ind;
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      if( g == 0 )
	 gp = 0;
      else
	 gp = gp + 3*ind;
      ind =0;
      for( int k=m_kStartAct[g]; k <= m_kEndAct[g]; k++ )
	 for( int j=m_jStartAct[g]; j <= m_jEndAct[g]; j++ )
	    for( int i=m_iStartAct[g]; i <= m_iEndAct[g]; i++ )
	    {
	       xm[gp+ind*3]   = mRho[g](i,j,k);
	       xm[gp+ind*3+1] = mMu[g](i,j,k);
	       xm[gp+ind*3+2] = mLambda[g](i,j,k);
	       ind++;
	    }
   }
}

//-----------------------------------------------------------------------
void EW::add_to_grad( vector<Sarray>& K, vector<Sarray>& Kacc, vector<Sarray>& Um, 
		      vector<Sarray>& U, vector<Sarray>& Up, vector<Sarray>& Uacc,
		      vector<Sarray>& gRho, vector<Sarray>& gMu, vector<Sarray>& gLambda )
{
   for( int g=0 ; g < mNumberOfGrids ; g++ )
   {
      int ifirst = m_iStart[g];
      int ilast  = m_iEnd[g];
      int jfirst = m_jStart[g];
      int jlast  = m_jEnd[g];
      int kfirst = m_kStart[g];
      int klast  = m_kEnd[g];
      int ifirstact = m_iStartAct[g];
      int ilastact  = m_iEndAct[g];
      int jfirstact = m_jStartAct[g];
      int jlastact  = m_jEndAct[g];
      int kfirstact = m_kStartAct[g];
      int klastact  = m_kEndAct[g];
      double* k_ptr = K[g].c_ptr();
      double* ka_ptr = Kacc[g].c_ptr();
      double* um_ptr = Um[g].c_ptr();
      double* u_ptr = U[g].c_ptr();
      double* up_ptr = Up[g].c_ptr();
      double* ua_ptr = Uacc[g].c_ptr();
      double* grho_ptr = gRho[g].c_ptr();
      double* gmu_ptr = gMu[g].c_ptr();
      double* glambda_ptr = gLambda[g].c_ptr();
      double h = mGridSize[g];
      int* onesided_ptr = m_onesided[g];
      if( topographyExists() && g == mNumberOfGrids-1 )
      {
	 F77_FUNC(addgradrhoc,ADDGRADRHOC)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
			    &ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact, &klastact,
					k_ptr, ka_ptr, um_ptr, u_ptr, up_ptr, ua_ptr, grho_ptr,
					     &mDt, mJ.c_ptr(), onesided_ptr );

      }
      else
      {
	 F77_FUNC(addgradrho,ADDGRADRHO)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                         &ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact, &klastact,
					k_ptr, ka_ptr, um_ptr, u_ptr, up_ptr, ua_ptr, grho_ptr,
					&mDt, &h, onesided_ptr );
         int nb = 4, wb=6;
	 F77_FUNC(addgradmula,ADDGRADMULA)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                         &ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact, &klastact,
					k_ptr, ka_ptr, um_ptr, u_ptr, up_ptr, ua_ptr, gmu_ptr,
				    glambda_ptr, &mDt, &h, onesided_ptr, &nb, &wb, m_bop );
      }
   }
}

//-----------------------------------------------------------------------
void EW::perturb_mtrl()
{
   int g=0;
   if( m_perturb != 0 && point_in_proc(m_iperturb,m_jperturb,g) )
   {
      cout << "per = " << m_perturb << " " << m_iperturb << " " << m_jperturb << " " << m_kperturb << endl;
      if( m_iperturb < m_iStartAct[g] || m_iperturb > m_iEndAct[g] )
	 cout << "warning i-index outside active domain " << endl;
      if( m_jperturb < m_jStartAct[g] || m_jperturb > m_jEndAct[g] )
	 cout << "warning j-index outside active domain " << endl;
      if( m_kperturb < m_kStartAct[g] || m_kperturb > m_kEndAct[g] )
	 cout << "warning k-index outside active domain " << endl;
      mMu[g](m_iperturb,m_jperturb,m_kperturb) += m_perturb;
      //      mRho[g](m_iperturb,m_jperturb,m_kperturb) += m_perturb;
   }
}

//-----------------------------------------------------------------------
void EW::get_optmethod( int& method, int& bfgs_m )
{
   method = m_opt_method;
   bfgs_m = m_lbfgs_m;
}
