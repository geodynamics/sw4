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
void F77_FUNC(forcingttfort,FORCINGTTFORT)(int*, int*, int*, int*, int*, 
				       int*, double*, double*, double*, double*, double*, double*, double*, 
				       double*, double*, double*, double*, double* );
void F77_FUNC(exactaccfort,EXACTACCFORT)(int*, int*, int*, int*, int*, int*, double*, double*, double*, 
					 double*, double*, double*, double* );
void F77_FUNC(rhserrfort, RHSERRFORT)(int*, int*, int*, int*, int*, int*, int*, double*,
				      double*, double*, double*, double*, double*);
void F77_FUNC(rhs4th3fort,RHS4TH3FORT)(int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*,
				       double*, double*, double*, double*, double*, double* );
void F77_FUNC(exactrhsfort,EXACTRHSFORT)( int*, int*, int*, int*, int*, int*, double*, double*, 
					  double*, double*, double*, double*, double*, double*, double*, double*,
					  double*, double* );
void F77_FUNC(solerr3, SOLERR3)(int*, int*, int*, int*, int*, int*, double*, double*, double*, double *li,
				double *, double* , double* , double* , double* , double* );
void F77_FUNC(solerrgp, SOLERRGP)(int*, int*, int*, int*, int*, int*, double*, double*, double*, double *li,
				double *l2 );
void F77_FUNC(twilightfort,TWILIGHTFORT)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, 
					  double*, double*, double* );
}

using namespace std;

#define SQR(x) ((x)*(x))

// constructor
EW::EW(const string& fileName, vector<Source*> & a_GlobalSources,
       vector<TimeSeries*> & a_GlobalTimeSeries): 
  m_topo_zmax(0.0),
  m_topoInputStyle(UNDEFINED), 
  mTopoImageFound(false),
  m_nx_base(0), m_ny_base(0), m_nz_base(0), m_h_base(0.0),
  mIsInitialized(false),
  mParsingSuccessful(false),
  mNumberOfGrids(0),
  mName(fileName),
  m_scenario(" "),
  mPath("./"),
  mWriteGMTOutput(false),
  mPlotFrequency(80),
  mNumFiles(0),
  mVerbose(0),
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
  m_twilight_forcing(false),
  m_point_source_test(0),
  m_energy_test(0),
  m_lamb_test(0),
  m_rayleigh_wave_test(0),
  m_update_boundary_function(0),

  //  mTestingEnergy(false),
  //  mTestSource(false),
  //  mTestLamb(false),
  mOrder(4),
  mCFL(1.3), // 0.8 for 2nd order
  // m_d4coeff(0.0),
  // m_d4_cfl(0.2),
  // m_curlcoeff(0.0),

  // mRestartFilePrefix(""),
  // mRestartFromCycle(0),
  // mRestartDumpInterval(0),

  mbcsSet(false),

  m_analytical_topo(false),
  m_GaussianAmp(0.),
  m_GaussianLx(1.0),
  m_GaussianLy(1.0),
  m_GaussianXc(0.0),
  m_GaussianYc(0.0),

  m_use_supergrid(false),
  m_sg_thickness_set(false),
  m_supergrid_thickness(-1.0),
  m_supergrid_damping_coefficient(0.15),

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
  m_number_material_surfaces(0),
  m_Nlon(0),
  m_Nlat(0),
  m_materialLon(NULL),
  m_materialLat(NULL),
  m_useVelocityThresholds(false),
  m_vpMin(0.),
  m_vsMin(0.),

  m_grid_interpolation_order(0),
  m_global_xmax(0.),
  m_global_ymax(0.),
  m_global_zmax(0.),
  m_global_zmin(0.),
  m_ghost_points(2), // for 4th order stencils
  m_ppadding(2),

  mLonOrigin(-118.0), // NTS
  mLatOrigin(37.0), // NTS
  mGeoAz(0.0), // x=North, y=East
  //  mDefaultLocation(true),
  mMetersPerDegree(111319.5),

// command limitfrequency
  m_limit_frequency(false), 
  m_ppw(15), 
  m_frequency_limit(1e38), // will hold min(Vs/h)/PPW

// command prefilter
  m_prefilter_sources(false), 
  m_fc(1.0),
  m_limit_source_freq(false),
  m_source_freq_max(1.0),

  mPrintInterval(100),
  m_t0Shift(0.0),
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
  m_att_ppw(15)
{
   MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
   MPI_Comm_size(MPI_COMM_WORLD, &m_nProcs);

// initialize the boundary condition array
   for (int i=0; i<6; i++)
   {
     mbcGlobalType[i] = bNone;
   }
   
// read the input file and setup the simulation object
   if (parseInputFile( a_GlobalSources, a_GlobalTimeSeries ))
     mParsingSuccessful = true;

}

// Destructor
EW::
~EW()
{
}

//-----------------------------------
bool EW::isInitialized()
{
  return mIsInitialized;
}
  
//-----------------------------------
bool EW::wasParsingSuccessful()
{
  return mParsingSuccessful;
}
  
//-----------------------------------------------------------------------
void EW::printTime( int cycle, double t, bool force ) const 
{
   if ( proc_zero() && (force || mPrintInterval == 1 ||
			(cycle % mPrintInterval) == 1 ||
			cycle == 1) )
// string big enough for >1 million time steps 
      printf("Time step %7i  t = %15.7e\n", cycle, t-m_t0Shift);
}
//-----------------------------------------------------------------------
void EW::printPreamble() const 
{
   stringstream msg;

   if ( proc_zero())
   {
      msg << "============================================================" << endl << endl 
          << " Running program on " << m_nProcs << " MPI tasks" << " using the following data: " << endl << endl
          << " Start Time = " << mTstart-m_t0Shift << " Goal Time = ";
      
      if (mTimeIsSet)
         msg << mTmax-m_t0Shift << endl;
      else
         msg << mNumberOfTimeSteps*mDt-m_t0Shift << endl;
      
      msg << " Number of time steps = " << mNumberOfTimeSteps << " dt: " << mDt << endl;
      
      if (mVerbose)
      {
//	msg << " Forcing = " << m_forcing->name() << endl;

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

   if( mVerbose >= 4 )
   {
// print out source info
// Note: Point sources are stored on the process that owns the grid point
     // for( int i=0 ; i < m_point_sources.size() ; i++ )
     // {
     //   cout << *m_point_sources[i];
     // }
   }

   cout.flush(); cerr.flush();
   MPI_Barrier(MPI_COMM_WORLD);
      
   // m0 values in each source command gets added up. This number is called the "Total seismic moment" 
   // and should be printed to stdout with the unit Nm (Newton-meter). If that number is M0, you should 
   // also print Mw = 2/3 *(log10(M0) - 9.1). That is the moment magnitude (dimensionless). 
      
   if( proc_zero() )
   {
      double myM0Sum = 0;
      int numsrc = 0, ignoredSources=0;
//       for (unsigned int i=0; i < mGlobalUniqueSources.size(); ++i)
//       {
//          if (mGlobalUniqueSources[i]->isMomentSource())
// 	 {
// // Note that proc 0 doen't know of all sources which need to be ignored
// 	   if (!mGlobalUniqueSources[i]->ignore() ) 
// 	   {
// 	     numsrc++;
//              myM0Sum += mGlobalUniqueSources[i]->getAmplitude();
// 	   }
// 	   else
// 	     ignoredSources++;
// 	 }
	 
//       }
      stringstream msg2;
      msg2 << endl
           << "-----------------------------------------------------------------------" << endl
           << "  Total seismic moment (M0): " << myM0Sum << " Nm " << endl;
      if (myM0Sum != 0)
         msg2 <<  "  Moment magnitude     (Mw): " << (2./3.)*(log10(myM0Sum) - 9.1)  << endl;
      msg2 << "  Number of sources " << numsrc << endl;
      msg2 << "-----------------------------------------------------------------------" << endl;
      cout << msg2.str();
  }
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

    if (m_jStart[top]+m_ghost_points > 1)
    {
      m_bcType[g][2] = bProcessor;
    }
    if (m_jEnd[top]-m_ghost_points < m_global_ny[top])
    {
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
// topography is not yet implemented
//
//     double zMinTilde;
//     int gCurv = mNumberOfGrids - 1;
//     double h = mGridSize[gCurv];
//     double q = x/h + 1.0;
//     double r = y/h + 1.0;

// // define the depth for ghost points to equal the depth on the nearest boundary point
//     double qMin = 1.0;
//     double qMax = (double) m_global_nx[gCurv];
//     double rMin = 1.0;
//     double rMax = (double) m_global_ny[gCurv];

//     if (q<qMin) q=qMin;
//     if (q>qMax) q=qMax;
//     if (r<rMin) r=rMin;
//     if (r>rMax) r=rMax;

// // evaluate elevation of topography on the grid (smoothed topo)
//     success=true;
//     if (!interpolate_topography(q, r, zMinTilde, true))
    {
      cerr << "ERROR: getDepth: Unable to evaluate topography for x=" << x << " y= " << y << " on proc # " << getRank() << endl;
      // cerr << "q=" << q << " r=" << r << " qMin=" << qMin << " qMax=" << qMax << " rMin=" << rMin << " rMax=" << rMax << endl;
      // cerr << "Setting elevation of topography to ZERO" << endl;
      success = false;
//      zMinTilde = 0;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
//    depth = z-zMinTilde;
  }
  return success;
}

//-----------------------------------------------------------------------
void EW::computeCartesianCoord(double &x, double &y, double lon, double lat)
{
  // -----------------------------------------------------------------
  // Compute the cartesian coordinate give the geographic coodinate
  // -----------------------------------------------------------------
  double deg2rad = M_PI/180.0;

  double phi = mGeoAz * deg2rad;

  // compute x and y
  x = mMetersPerDegree*(cos(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*sin(phi));
  y = mMetersPerDegree*(-sin(phi)*(lat-mLatOrigin) + cos(lat*deg2rad)*(lon-mLonOrigin)*cos(phi));

}


//-----------------------------------------------------------------------
void EW::computeGeographicCoord(double x, double y, double & longitude, double & latitude)
{
  // conversion factor between degrees and radians
  double deg2rad = M_PI/180.0;
  double phi = mGeoAz * deg2rad;

  // Compute the latitude
 latitude = mLatOrigin + 
    (x*cos(phi) - y*sin(phi))/mMetersPerDegree;

  // Compute the longitude
 longitude = mLonOrigin + 
     (x*sin(phi) + y*cos(phi))/(mMetersPerDegree*cos(latitude*deg2rad));

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
              << m_global_nx[a_g] << ")");
      VERIFY2(a_j >= 1-m_ghost_points && a_j <= m_global_ny[a_g]+m_ghost_points,
              "Grid Error: j (" << a_j << ") is out of bounds: ( " << 1 << ","
              << m_global_ny[a_g] << ")");
      VERIFY2(a_k >= m_kStart[a_g] && a_k <= m_kEnd[a_g],
              "Grid Error: k (" << a_k << ") is out of bounds: ( " << 1 << "," 
              << m_kEnd[a_g]-m_ghost_points << ")");
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
void EW::getGMTOutput( vector<Source*> & a_GlobalUniqueSources )
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
      int minx = static_cast<int>(min(lonSW-1,min(lonSE-1,min(lonNE-1,lonNW-1))));
      int maxx = static_cast<int>(max(lonSW+1,max(lonSE+1,max(lonNE+1,lonNW+1))));
      int miny = static_cast<int>(min(latSW-1,min(latSE-1,min(latNE-1,latNW-1))));
      int maxy = static_cast<int>(max(latSW+1,max(latSE+1,max(latNE+1,latNW+1)))); 

//      GeographicCoord eNW, eNE, eSW, eSE;
      
#ifdef ENABLE_ETREE
      if (mEtreeFile != NULL)
      {
//        mEtreeFile->getGeoBox()->getBounds(eNW, eNE, eSW, eSE);
        
        minx = static_cast<int>(min(eSW.getLongitude()-1,min(eSE.getLongitude()-1,min(eNE.getLongitude()-1,eNW.getLongitude()-1))));
        maxx = static_cast<int>(max(eSW.getLongitude()+1,max(eSE.getLongitude()+1,max(eNE.getLongitude()+1,eNW.getLongitude()+1))));
        miny = static_cast<int>(min(eSW.getLatitude()-1,min(eSE.getLatitude()-1,min(eNE.getLatitude()-1,eNW.getLatitude()-1))));
        maxy = static_cast<int>(max(eSW.getLatitude()+1,max(eSE.getLatitude()+1,max(eNE.getLatitude()+1,eNW.getLatitude()+1)))); 
      }
#endif
      
      contents << "# Region will need to be adjusted based on etree/grid values" << endl
               << "set REGION = " << minx << "/" << maxx << "/" << miny << "/" << maxy << endl
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
               << "pscoast -R$REGION -JM$SCALE -Ba1 -Dhigh -S200,200,255 -A2000 -W3 -N1t3 -N2t2a -K >! plot.ps" << endl << endl
               << "# WPP grid region..." << endl;
      
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
//          GeographicCoord eNW, eNE, eSW, eSE;
//          mEtreeFile->getGeoBox()->getBounds(eNW, eNE, eSW, eSE);
         contents << "# Etree region: " << mEtreeFile->getFileName() << endl
                  << "psxy -R$REGION -JM$SCALE -W5/255/255/0ta -O -K <<EOF>> plot.ps" << endl
                  << eNW.getLongitude() << " " << eNW.getLatitude() << endl
                  << eNE.getLongitude() << " " << eNE.getLatitude() << endl
                  << eSE.getLongitude() << " " << eSE.getLatitude() << endl
                  << eSW.getLongitude() << " " << eSW.getLatitude() << endl
                  << eNW.getLongitude() << " " << eNW.getLatitude() << endl
                  << "EOF" << endl << endl;
      }
#endif
      
      if (a_GlobalUniqueSources.size() > 0)
      {
         contents << "# Sources... " << endl
                  << "psxy -R$REGION -JM$SCALE -Sd0.1 -Gred -O -K <<EOF>> plot.ps" << endl;
         
         for (int i=0; i < a_GlobalUniqueSources.size(); ++i)
         {
           double latSource,lonSource;

           computeGeographicCoord(a_GlobalUniqueSources[i]->getX0(), a_GlobalUniqueSources[i]->getY0(),
                                  lonSource ,latSource);
// is this the correct syntax???
	   contents << latSource << " " << lonSource << endl;
         }
         contents << "EOF" << endl << endl;
      }
      
      int numStations = 0;
      stringstream stationstr;
      stationstr << "# Stations... " << endl;  
      stationstr << "cat << EOF >! stations.d " << endl;
      // Write stations by rereading the WPP input file, since some might
      // live outside the grid...
      ifstream wppfile(mWPPFileName.c_str());
      if (!wppfile.is_open())
         contents << "# Error opening wpp input file, skipping stations" << endl;
      else
      {
         char buffer[256];
         while (!wppfile.eof())
         { 
            wppfile.getline(buffer, 256);
            if (startswith("sac", buffer))
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
         } // line in wpp file
      }
      
      stationstr << "EOF" << endl << endl;
      
      stationstr << "# plot station names" << endl
                 << "psxy -R -J -O -K -St0.1 -Gblue -W0.25p stations.d >> plot.ps" << endl
                 << "awk '{print $1, $2, 12, 1, 9, $4, $3}' stations.d | pstext -R -J -O -Dj0.3/0.3v -Gblue -N >> plot.ps" << endl;
      
      // Only write station info if there are stations.
      if (numStations > 0) contents << stationstr.str() << endl;

      contents << "/bin/mv plot.ps " << mWPPFileName << ".ps" << endl;

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
   if( proc_zero() )
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
   if( proc_zero() )
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
void EW::default_bcs( boundaryConditionType bcs[6] )
{
   for( int side=0 ; side < 6 ; side++ )
      bcs[side] = bDirichlet;
   bcs[4] = bStressFree;
}

//---------------------------------------------------------------------------
void EW::normOfDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, 
                           double &diffL2, vector<Source*>& a_globalUniqueSources )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *uex_ptr, *u_ptr, h, linfLocal=0, l2Local=0, diffInfLocal=0, diffL2Local=0;
  double radius =-1, x0, y0, z0;

//tmp  
  // if (proc_zero())
  //   printf("Inside normOfDifference\n");
  
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

// tmp
//     printf("proc=%i, iS= %i, iE=%i, jS=%i, jE=%i, kS=%i, kE=%i\n", m_myRank, 
// 	   m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g]);
//     printf("proc=%i, if= %i, il=%i, jf=%i, jl=%i, kf=%i, kl=%i\n", m_myRank, 
// 	   ifirst, ilast, jfirst, jlast, kfirst, klast);

    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?

    if( m_point_source_test )
    {
       radius = 4*h;
       x0 = a_globalUniqueSources[0]->getX0();
       y0 = a_globalUniqueSources[0]->getY0();
       z0 = a_globalUniqueSources[0]->getZ0();
    }

// need to exclude parallel overlap from L2 calculation
    F77_FUNC(solerr3, SOLERR3)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &h,
				uex_ptr, u_ptr, &linfLocal, &l2Local, &m_zmin[g], &x0,
				&y0, &z0, &radius );
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
void EW::bndryInteriorDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, 
				  double lowZ[3], double interiorZ[3], double highZ[3] )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nz;
  double *uex_ptr, *u_ptr, h, li, l2;
  
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
    nz = m_global_nz[g];
    
// need to do a gather over all processors
    F77_FUNC(rhserrfort, RHSERRFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz, &h,
				      uex_ptr, u_ptr, lowZ, interiorZ, highZ);
  }
}

//---------------------------------------------------------------------------
void EW::test_RhoUtt_Lu( vector<Sarray> & a_Uacc,  vector<Sarray> & a_Lu,   vector<Sarray> & a_F, 
			 double lowZ[3], double interiorZ[3], double highZ[3] )
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
				      &nz, uacc_ptr, lu_ptr, f_ptr, rho_ptr, lowZ, interiorZ, highZ);
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
  }
  else if( m_rayleigh_wave_test )
  {
     // NYI
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
     retval = true;
  }
  else if( m_point_source_test )
  {
     for(int g=0 ; g < mNumberOfCartesianGrids; g++ )
        get_exact_point_source( a_U[g], a_t, g, *sources[0] );
     retval = true;
  }
  else if( m_lamb_test )
  {
     retval = false;
  }
  else if( m_rayleigh_wave_test )
  {
     retval = false;
  }
  else
  {
     // Exact solution unknown
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
void EW::get_exact_point_source( Sarray& u, double t, int g, Source& source )
{
   timeDep tD;
   if(!( source.getName() == "SmoothWave" || source.getName() == "VerySmoothBump" || source.getName()== "Gaussian") )
   {
      cout << "EW::get_exact_point_source: Error, time dependency must be SmoothWave, VerySmoothBump, or Gaussian, not "
	   << source.getName() << endl;
      return;
   }
   else if( source.getName() == "SmoothWave" )
      tD = iSmoothWave;
   else if( source.getName() == "VerySmoothBump" )
      tD = iVerySmoothBump;
   else
      tD = iGaussian;

   u.set_to_zero();
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
      source.getMoments( mxx, myy, mzz, mxy, mxz, myz );
      m0  = source.getAmplitude();
   }
   double* up = u.c_ptr();
   double h   = mGridSize[g];
   double eps = 1e-3*h;
   size_t ind = 0;
   // Note: Use of ind, assumes loop is over the domain over which u is defined.
   for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
      for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	 for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	 {
	    double x = (i-1)*h;
	    double y = (j-1)*h;
	    double z = (k-1)*h + m_zmin[g];
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

    F77_FUNC(exactrhsfort,EXACTRHSFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					 &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
					 &h, &zmin );
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
}

//---------------------------------------------------------------------------
void EW::exactForce(double a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources )
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
	F77_FUNC(forcingfort,FORCINGFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					   &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
					   &h, &zmin );
     }
  }
  else if( m_rayleigh_wave_test )
  {
     // NYI
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
void EW::exactForce_tt(double a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources )
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
     }

     //  subroutine forcingfort( ifirst, ilast, jfirst, jlast, kfirst, 
     // +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     // +     h, zmin)
     F77_FUNC(forcingttfort,FORCINGTTFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					    &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
					    &amprho, &ampmu, &ampla, &h, &zmin );
  }
  else if( m_rayleigh_wave_test )
  {
     // NYI
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
void EW::evalRHS(vector<Sarray> & a_U, vector<Sarray> & a_Uacc )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *uacc_ptr, *u_ptr, *mu_ptr, *la_ptr, *rho_ptr, h;
  
  int *onesided_ptr;
  
  int g, nz;
  
  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
  {
    uacc_ptr = a_Uacc[g].c_ptr();
    u_ptr   = a_U[g].c_ptr();
    mu_ptr  = mMu[g].c_ptr();
    la_ptr  = mLambda[g].c_ptr();
    rho_ptr = mRho[g].c_ptr();
    ifirst = m_iStart[g];
    ilast  = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast  = m_jEnd[g];
    kfirst = m_kStart[g];
    klast  = m_kEnd[g];
    h = mGridSize[g]; // how do we define the grid size for the curvilinear grid?
    nz = m_global_nz[g];
    onesided_ptr = m_onesided[g];
    
    F77_FUNC(rhs4th3fort,RHS4TH3FORT)(&ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				      &klast, &nz, onesided_ptr, m_acof, m_bope, m_ghcof,
				      uacc_ptr, u_ptr, mu_ptr, la_ptr, rho_ptr, &h );    

  }
}

//---------------------------------------------------------------------------
void EW::evalPredictor(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		       vector<Sarray> & a_Lu, vector<Sarray> & a_F )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *up_ptr, *u_ptr, *um_ptr, *lu_ptr, *fo_ptr, *rho_ptr, dt2;
  
  int *onesided_ptr;
  
  int g, nz;
  
  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    u_ptr   = a_U[g].c_ptr();
    um_ptr  = a_Um[g].c_ptr();
    lu_ptr  = a_Lu[g].c_ptr();
    fo_ptr  = a_F[g].c_ptr();
    rho_ptr = mRho[g].c_ptr();
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
void EW::evalCorrector(vector<Sarray> & a_Up, vector<Sarray> & a_Lu, vector<Sarray> & a_F )
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *up_ptr, *lu_ptr, *fo_ptr, *rho_ptr, dt4;
  
  int g;
  
  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
  {
    up_ptr  = a_Up[g].c_ptr();
    lu_ptr  = a_Lu[g].c_ptr();
    fo_ptr  = a_F[g].c_ptr();
    rho_ptr = mRho[g].c_ptr();
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
  
  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
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
void EW::update_images( int currentTimeStep, double time, vector<Sarray> & mUp )
{
  double maxerr;
  for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
  {
    Image* img = mImageFiles[fIndex];
    
    // if( img->mMode == Image::HVELMAX)
    //   img->update_maxes_hVelMax();
    // if( img->mMode == Image::VVELMAX)
    //   img->update_maxes_vVelMax();

    if (img->timeToWrite(time - m_t0Shift, currentTimeStep, mDt )) // subtract m_t0Shift from time to get actual time
    {
      if(img->mMode == Image::UX ) 
      {
	img->computeImageQuantity(mUp, 1);
      }
      else if(img->mMode == Image::UY )
      { 
	img->computeImageQuantity(mUp, 2);
      }
      else if(img->mMode == Image::UZ )
      {
	img->computeImageQuantity(mUp, 3);
      }
      // else if(img->mMode == Image::FX ) 
      // {
      // 	img->computeImageQuantity(mF, 1);
      // }
      // else if(img->mMode == Image::FY )
      // { 
      // 	img->computeImageQuantity(mF, 2);
      // }
      // else if(img->mMode == Image::FZ )
      // {
      // 	img->computeImageQuantity(mF, 3);
      // }
      else if(img->mMode == Image::RHO )
      {
	img->computeImageQuantity(mRho, 1);
      }
      else if(img->mMode == Image::MU )
      {
	img->computeImageQuantity(mMu, 1);
      }
      else if(img->mMode == Image::LAMBDA )
      {
	img->computeImageQuantity(mLambda, 1);
      }
//       else if(img->mMode == Image::QS )
//       { 
// 	if (usingAttenuation())
// 	  img->computeImageQuantity(mQs, 2);
//       }
//       else if(img->mMode == Image::QP )
//       { 
// 	if (usingAttenuation())
// 	  img->computeImageQuantity(mQp, 2);
//       }
// // the next 3 use Image::computeImageError to fill in values, which calls the forcing to get the exact solution
//       else if(img->mMode == Image::UXERR ) 
//       {
// 	img->computeImageError(mUp, 1);
//       }
//       else if(img->mMode == Image::UYERR )
//       { 
// 	img->computeImageError(mUp, 2);
//       }
//       else if(img->mMode == Image::UZERR )
//       {
// 	img->computeImageError(mUp, 3);
//       }
//       else if(img->mMode == Image::DIV )
//       {
// 	img->computeDivergence();
// 	if (m_forcing->knows_exact())
//         {
//           maxerr = img->computeImageErrorDebug(0);
//           if (proc_zero())
//             printf("maxErr DIV %f @ %fs\n",maxerr,time - m_t0Shift);
//         }
//       }
//       else if(img->mMode == Image::CURL )
//       {
//         img->computeCurlMagnitude();
// 	if (m_forcing->knows_exact())
//         {
//           maxerr = img->computeImageErrorDebug(2);
//           if (proc_zero())
//             printf("maxErr CURL %f @ %f s\n",maxerr,time - m_t0Shift);
//         }
//       }
//       else if(img->mMode == Image::VELDIV )
//       {
// 	img->computeVelDivergence();
// 	if (m_forcing->knows_exact())
//         {
//           maxerr = img->computeImageErrorDebug(1);
//           if (proc_zero())
//             printf("maxErr VELDIV %f @ %f s\n",maxerr,time - m_t0Shift);
//         }
//       }
//       else if(img->mMode == Image::VELCURL )
//       {
//         img->computeVelCurlMagnitude();
// 	if (m_forcing->knows_exact())
//         {
//           maxerr = img->computeImageErrorDebug(3);
//           if (proc_zero())
//             printf("maxErr VELCURL %f @ %f s\n",maxerr,time - m_t0Shift);
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
//       else if(img->mMode == Image::P )
//       {
// 	img->computeImageP(mMu, mLambda, mRho);
//       }
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
      else
//      else if (!img->mMode==Image::HVELMAX||!img->mMode==Image::VVELMAX)
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

      img->writeImagePlane_2(currentTimeStep, mPath );
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
  //   if(img->timeToWrite(time - m_t0Shift, currentTimeStep, mDt ) ) // subtract m_t0Shift from time to get actual time
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
}
