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
				double *l2 );
void F77_FUNC(solerrgp, SOLERRGP)(int*, int*, int*, int*, int*, int*, double*, double*, double*, double *li,
				double *l2 );
void F77_FUNC(twilightfort,TWILIGHTFORT)( int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, 
					  double*, double*, double* );
}

using namespace std;

#define SQR(x) ((x)*(x))

// constructor
EW::EW(const string& fileName, vector<Source*> & a_GlobalUniqueSources): 
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
  m_twilight_forcing(false),
  m_update_boundary_function(0),

  mTestingEnergy(false),
  mTestSource(false),
  mTestLamb(false),
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

//   mGeoCoord(0.0,0.0,0.0),
  mGeoCoord(37.0,-118.0,0.0), // NTS
  mGeoAz(135.0),
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
   if (parseInputFile( a_GlobalUniqueSources ))
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
void EW::computeCartesianCoord(double &x, double &y, double &z,const GeographicCoord& g)
{
   double lat = g.getLatitude();
   double lon = g.getLongitude();
   z = g.getDepth();

//    REQUIRE2(lat >= -90.0 && lat <= 90.0,
// 	    "latitude must be between -180 and 180 deg, not " << lat);
//    REQUIRE2(lon >= -180.0 && lon <= 180.0,
// 	    "longitude must be between -180 and 180 deg, not " << lon);
//    REQUIRE2(z >= 0.0,
// 	    "depth must be positive, not " << z);
  // -----------------------------------------------------------------
  // Compute the cartesian coordinate give the geographic coodinate
  // -----------------------------------------------------------------
  double deg2rad = M_PI/180.0;

  // First figure out the grid offset
  double gridLat = mGeoCoord.getLatitude();
  double gridLon = mGeoCoord.getLongitude();

  double phi = mGeoAz * deg2rad;

  // compute x and y
  x = mMetersPerDegree*(cos(phi)*(lat-gridLat) + cos(lat*deg2rad)*(lon-gridLon)*sin(phi));
  y = mMetersPerDegree*(-sin(phi)*(lat-gridLat) + cos(lat*deg2rad)*(lon-gridLon)*cos(phi));

//   REQUIRE2(x >= 0.0,
// 	  "x must be positive, not: " << x);
//   REQUIRE2(y >= 0.0,
// 	  "y must be positive, not: " << y);
//   REQUIRE2(z >= 0.0,
// 	  "z must be positive, not: " << z);

}


//-----------------------------------------------------------------------
void EW::computeGeographicCoord(double x, double y, double z, double & latitude, double & longitude)
{
  // Put the angle in radians
  double deg2rad = M_PI/180.0;
  double phi = mGeoAz * deg2rad;

  // Compute the latitude
 latitude = mGeoCoord.getLatitude() + 
    (x*cos(phi) - y*sin(phi))/mMetersPerDegree;

  // Compute the longitude
 longitude = mGeoCoord.getLongitude() + 
     (x*sin(phi) + y*cos(phi))/(mMetersPerDegree*cos(latitude*deg2rad));

// what about the depth?
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
     retval = (a_i >= (m_iStart[a_g]+m_paddingCells[0])) && (a_i <= (m_iEnd[a_g]-m_paddingCells[1])) &&   
              (a_j >= (m_jStart[a_g]+m_paddingCells[2])) && (a_j <= (m_jEnd[a_g]-m_paddingCells[3]));
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
      computeGeographicCoord(0.0,0.0,0.0,latSW,lonSW);
      computeGeographicCoord((m_global_nx[0]-1)*mGridSize[0],0.0,0.0,latSE,lonSE);
      computeGeographicCoord((m_global_nx[0]-1)*mGridSize[0],(m_global_ny[0]-1)*mGridSize[0],0.0,latNE,lonNE);
      computeGeographicCoord(0.0,(m_global_ny[0]-1)*mGridSize[0],0.0,latNW,lonNW);
     
      // Round up/down
      int minx = static_cast<int>(min(lonSW-1,min(lonSE-1,min(lonNE-1,lonNW-1))));
      int maxx = static_cast<int>(max(lonSW+1,max(lonSE+1,max(lonNE+1,lonNW+1))));
      int miny = static_cast<int>(min(latSW-1,min(latSE-1,min(latNE-1,latNW-1))));
      int maxy = static_cast<int>(max(latSW+1,max(latSE+1,max(latNE+1,latNW+1)))); 

      GeographicCoord eNW, eNE, eSW, eSE;
      
#ifdef ENABLE_ETREE
      if (mEtreeFile != NULL)
      {
        mEtreeFile->getGeoBox()->getBounds(eNW, eNE, eSW, eSE);
        
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

           computeGeographicCoord(a_GlobalUniqueSources[i]->getX0(),
                                  a_GlobalUniqueSources[i]->getY0(),
                                  a_GlobalUniqueSources[i]->getZ0(),
                                  latSource                       ,
                                  lonSource                       );
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
                 computeGeographicCoord(x, y, z,lat,lon);
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
void EW::normOfDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, double &diffL2 )
{
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *uex_ptr, *u_ptr, h, linfLocal=0, l2Local=0, diffInfLocal=0, diffL2Local=0;

//tmp  
  if (proc_zero())
    printf("Inside normOfDifference\n");
  
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
    
// need to exclude parallel overlap from L2 calculation
    F77_FUNC(solerr3, SOLERR3)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &h,
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
void EW::exactSolTwilight(double a_t, vector<Sarray> & a_U, vector<Sarray*> & a_AlphaVE)
{
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  double *u_ptr, om, ph, cv, h, zmin;
  
  int g;
  
  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
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
    
    if (m_twilight_forcing)
    {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
    }

    F77_FUNC(twilightfort,TWILIGHTFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					 &klast, u_ptr, &a_t, &om, &cv, &ph, &h, &zmin );

    // m_twilight_forcing->get_initial_data_Cartesian( m_zmin[g], mGridSize[g], a_t, a_U[g] );
    // for( int a=0 ; a < m_number_mechanisms ; a++ )
    // {
    //   for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
    // 	for (int j=m_jStart[g]; j<=m_jEnd[g]; j++)
    // 	  for (int i=m_iStart[g]; i<=m_iEnd[g]; i++)
    // 	  {
    // 	    xp = (i-1)*mGridSize[g];
    // 	    yp = (j-1)*mGridSize[g];
    // 	    zp = m_zmin[g] + (k-1)*mGridSize[g];
    // 	    m_twilight_forcing->get_exact_att( xp, yp, zp, a_t, al );
    // 	    a_AlphaVE[g][a](1,i,j,k) = al[0];
    // 	    a_AlphaVE[g][a](2,i,j,k) = al[1];
    // 	    a_AlphaVE[g][a](3,i,j,k) = al[2];
    // 	  }
    // }
  } // end for all Cartesian grids
      
// Exact solution for the curvilinear grid
  // if ( topographyExists() )
  // {
  //   g = mNumberOfGrids-1;
  //   m_twilight_forcing->get_initial_data_Curvilinear( mX, mY, mZ, a_t, a_U[g] );
  //   for( int a=0 ; a < m_number_mechanisms ; a++ )
  //   {
  //     double xp, yp, zp, al[3];
  //     for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
  // 	for (int j=m_jStart[g]; j<=m_jEnd[g]; j++)
  // 	  for (int i=m_iStart[g]; i<=m_iEnd[g]; i++)
  // 	  {
  // 	    xp = mX(i,j,k);
  // 	    yp = mY(i,j,k);
  // 	    zp = mZ(i,j,k);
  // 	    m_twilight_forcing->get_exact_att( xp, yp, zp, a_t, al );
  // 	    a_AlphaVE[g][a](1,i,j,k) = al[0];
  // 	    a_AlphaVE[g][a](2,i,j,k) = al[1];
  // 	    a_AlphaVE[g][a](3,i,j,k) = al[2];
  // 	  }
  //   }
  // }
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
void EW::exactForceTwilight(double a_t, vector<Sarray> & a_F)
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

     //  subroutine forcingfort( ifirst, ilast, jfirst, jlast, kfirst, 
     // +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     // +     h, zmin)
    F77_FUNC(forcingfort,FORCINGFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
				       &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
				       &h, &zmin );
  }
}

//---------------------------------------------------------------------------
void EW::exactForce_ttTwilight(double a_t, vector<Sarray> & a_F)
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

     //  subroutine forcingfort( ifirst, ilast, jfirst, jlast, kfirst, 
     // +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     // +     h, zmin)
    F77_FUNC(forcingttfort,FORCINGTTFORT)( &ifirst, &ilast, &jfirst, &jlast, &kfirst, 
					   &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
					    &amprho, &ampmu, &ampla, &h, &zmin );
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
void EW::addSAC(SAC s)
{
  mSACOutputFiles.push_back(s);
}


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
