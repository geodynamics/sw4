#include "mpi.h"

#include "EW.h"

#include "version.h"
#include "Require.h"
#include "nearlyEqual.h"
#include "boundaryConditionTypes.h"
#include "MaterialBlock.h"
#include "TimeSeries.h"

// #include "Image3D.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>

using namespace std;

#define SQR(x) ((x)*(x))

int gcd( int a, int b )
{
   // Euclidean algorithm
   while( b != 0 )
   {
      int t = b;
      b = a % b;
      a = t;
   }
   return a;
}

int computeEndGridPoint( double maxval, double dh )
{
  // We round up one, so that the end point
  // specified by the user is always included
  // in the domain.  i.e. if z was specified
  // as 15.0, and dh was computed to be 3.33,
  // the number of grid points would be 15.0/3.33 + 1
  // or 5, giving us a new max z = 16.64.
  int pts = 0;
  double x = 0.0;

  while (x < maxval && !dbc::nearlyEqual(x, maxval) )
    {
      x += dh;
      pts++;
    }
  
  // 1 based indexing
  pts++;

  return pts;
}

//-----------------------------------------------------------------------
bool endswith(string end, string& mystr)
{
   int lenEnd = end.length();
   int lenStr = mystr.length();

   if (lenEnd > lenStr) return false;

   cout << "mystr: " << mystr << " end: " << end << " " << mystr.substr(lenStr-lenEnd, lenEnd) << endl;

   if (mystr.substr(lenStr-lenEnd, lenEnd) == end)
      return true;
   else 
      return false;
}
//-----------------------------------------------------------------------
bool startswith(const char begin[], char *line)
{
  int lenb = strlen(begin);

  // We ignore any preceeding whitespace
  while (strncmp(line, " ", 1) == 0 || strncmp(line, "\t", 1) == 0)
    line++;
    
  if (strncmp(begin, line, lenb) == 0)
     return true;
  else 
     return false;
}

void EW::deprecatedOption(const string& command, 
		      const string& oldone, 
		      const string& newone)
{
  if (m_myRank == 0)
    cout << "DeprecationWarning: " 
	 << command << " option " << oldone << " is no longer supported.  Use "
	 << newone << " instead." << endl;
}

void unchecked(const char* cmd)
{
   cout << "*** Not yet checking command: " << cmd << endl;
}

void checked(const char* cmd)
{
   cout << "*** " << cmd << " command checked! " << endl;
}

//-----------------------------------
// 
// Note that parseInputFile() calls a lot of member functions of the EW class that
// should not be called after the initialization of the EW object is completed. 
// Make all these functions private!
//
bool EW::parseInputFile( vector<Source*> & a_GlobalUniqueSources,
			 vector<TimeSeries*> & a_GlobalTimeSeries )
{
  char buffer[256];
  ifstream inputFile;
  int blockCount=0;

  MPI_Barrier(MPI_COMM_WORLD);
  double time_start = MPI_Wtime();

  inputFile.open(mName.c_str());
  if (!inputFile.is_open())
  {
    if (m_myRank == 0)
      cerr << endl << "ERROR OPENING INPUT FILE: " << mName << endl << endl;
    return false;
  }
  
  bool foundGrid = false;

// First process Geodyn input for restrictions of allowable grid sizes.
  // while (!inputFile.eof())
  // {
  //    inputFile.getline(buffer, 256);
  //    if( startswith("geodynbc",buffer ) )
  // 	geodynFindFile(buffer);
  // }
  // inputFile.clear();
  // inputFile.seekg(0, ios::beg);
  
//---------------------------------------------------------------
// First process the grid, refinement and topography commands so
// we know how big the solution arrays need to be.
//
// Also, if we are enabling attenuation turn it on now so it
// can be read in as the material command blocks are initialized
//---------------------------------------------------------------

// these commands can enter data directly the object (this->)
  while (!inputFile.eof())
  {    
     inputFile.getline(buffer, 256);
     if( startswith("grid", buffer) )
     {
       foundGrid = true;
       processGrid(buffer);
     }
// read fileio here to enable verbose warnings in other commands
     else if (startswith("fileio", buffer))
     {
       processFileIO(buffer);
     }
     
     // else if (startswith("refinement", buffer))
     // {
     //   processRefinement(buffer);
     // }
     // else if (startswith("topography", buffer))
     // {
     //   processTopography(buffer);
     // }
     // else if (startswith("attenuation", buffer))
     // {
     //   processAttenuation(buffer);
     // }
     // else if (startswith("efile", buffer))
     // {
     //   getEfileInfo(buffer); // read efile name, etc for setting up topography from the efile
     // }
  }

// make sure there was a grid command
  if (!foundGrid)
  {
    if (m_myRank == 0)
    {
      cerr << "Error: No grid found in input file: " << mName << endl;
      return false; // unsuccessful
    }
    
  }

// sort and correct vector 'm_refinementBoundaries'. Initialize if not already available
  cleanUpRefinementLevels();
  
  inputFile.clear();
  inputFile.seekg(0, ios::beg); // reset file pointer to the beginning of the input file

// we only allocate solution arrays for the Cartesian grids and call a placeholder grid generator for the curvilinear grid
  allocateCartesianSolverArrays(m_global_zmax); 

  if (m_use_attenuation)
    setAttenuationParams(m_number_mechanisms, m_velo_omega, m_att_ppw, m_att_max_frequency );

// setup 2D communicators on the finest grid so that we can smooth the topography
  setup2D_MPICommunications();

// deal with topography
  if (m_topography_exists)
  {
    if (m_myRank == 0)
      cout << endl << 
	"*** ERROR: Topography not yet implemented ***" << endl << endl;
    return false;

// // 1. read topography from efile 
//       if (m_topoInputStyle == FileInput::Efile)
//       {
// 	extractTopographyFromEfile(m_topoFileName, m_topoExtFileName, m_QueryType,
//                                                 m_EFileResolution);
// // 2. smooth the topo
// 	smoothTopography(m_maxIter);
//       }
//       else if (m_topoInputStyle == FileInput::GridFile)
//       {
// // 1. read topography from grid file
// 	extractTopographyFromGridFile(m_topoFileName);
// // 2. smooth the topo
// 	smoothTopography(m_maxIter);
//       }
//       else if (m_topoInputStyle == FileInput::CartesianGrid)
//       {
// // 1. read topography from Cartesian grid file
// 	extractTopographyFromCartesianFile(m_topoFileName);
// // 2. smooth the topo
// 	smoothTopography(m_maxIter);
//       }
//       else if (m_topoInputStyle == FileInput::GaussianHill)
//       {
// 	buildGaussianHillTopography(m_GaussianAmp, m_GaussianLx, m_GaussianLy, m_GaussianXc, m_GaussianYc);
//       }
      
// // 3. Figure out the number of grid points in the vertical direction and allocate solution arrays on the curvilinear grid
//       allocateCurvilinearArrays(); // need to assign  m_global_nz[g] = klast - m_ghost_points; + allocate mUacc

  }
  else
  {
    if (m_myRank == 0)
      cout << endl << 
	"*** No topography command found in input file. Using z=0 as free surface boundary ***" << endl << endl;
  }
    
// setup communicators for 3D solutions on all grids
  setupMPICommunications();

// make the grid, allocate arrays for the curvilinear grid
  if (m_topography_exists)
  {
    generate_grid();
    setup_metric();
// note that the topo image can not be made until WPP2::InitializePaddingCells() has been called
  }

// output grid size info
  if (m_myRank == 0)
  {
    int nx, ny, nz;
    double nTot=0.;
    printf("\nGlobal grid sizes (without ghost points)\n");
//             1234  12345679  12345679  12345679  12345679
    printf("Grid         h        Nx        Ny        Nz       Points\n");
    for (int g = 0; g < mNumberOfGrids; g++)
    {
      nx = m_global_nx[g];
      ny = m_global_ny[g];
      nz = m_kEnd[g] - m_ghost_points;
      nTot += ((long long int)nx)*ny*nz;
      printf("%4i %9g %9i %9i %9i %12lld\n", g, mGridSize[g], nx, ny, nz, ((long long int)nx)*ny*nz);
    }
    printf("Total number of grid points (without ghost points): %g\n\n", nTot);
      
  }

  //----------------------------------------------------------
  // Now onto the rest of the input file...
  //----------------------------------------------------------
  while (!inputFile.eof())
  {
     inputFile.getline(buffer, 256);

     if (strlen(buffer) > 0) // empty lines produce this
     {
       if (startswith("#", buffer) || 
	   startswith("grid", buffer) ||
	   startswith("refinement", buffer) || 
	   startswith("topography", buffer) || 
	   startswith("attenuation", buffer) || 
	   startswith("fileio", buffer) ||
	   startswith("\n", buffer) || buffer == "\r" || buffer == "\0")
       {
	 // Ignore commented lines, newlines,
	 // grid, refinement, fileio, topography, and attenuation since we already processed those commands.
       }
       // else if (startswith("gmt", buffer))
       //   processGMT(buffer);
       else if (startswith("time", buffer))
	 processTime(buffer);
       // else if (startswith("globalmaterial", buffer))
       //   processGlobalMaterial(buffer);
       else if (startswith("receiver", buffer)) // was called "sac" in WPP
	 processReceiver(buffer, a_GlobalTimeSeries);
       // else if (startswith("energy", buffer))
       //   processEnergy(buffer);
       else if (startswith("twilight", buffer))
	 processTwilight(buffer);
       else if (startswith("testpointsource", buffer))
	 processTestPointSource(buffer);
       else if (startswith("testlamb", buffer))
         processTestLamb(buffer);
       else if (startswith("source", buffer))
	 processSource(buffer, a_GlobalUniqueSources);
       else if (startswith("block", buffer))
	 processMaterialBlock(buffer, blockCount);
//    else if (startswith("efile", buffer))
//    {
// #ifndef ENABLE_ETREE
//      if (m_myRank==0) 
//        cout << "Error: WPP was not built with Etree (efile) support" << endl;
//      return false;
// #endif
//      processMaterialEtree(buffer);
//    }
       // else if (startswith("pfile", buffer))
       //   processMaterialPfile(buffer);
       // else if (startswith("ifile", buffer))
       //   processMaterialIfile(buffer);
       // else if (startswith("material", buffer))
       //   processMaterial(buffer);
       else if (startswith("image", buffer))
         processImage(buffer);
       // else if (startswith("volimage", buffer))
       //   processImage3D(buffer);
//    else if (startswith("restart", buffer))
//       processRestart(buffer);
       else if (startswith("boundary_conditions", buffer))
         processBoundaryConditions(buffer);
       // else if (startswith("supergrid", buffer))
       //   processSupergrid(buffer);
       // else if (startswith("prefilter", buffer))
       //   processPrefilter(buffer);
//    else if( startswith("checkfornan", buffer ) )
//       processCheckfornan(buffer);
       // else if( startswith("developer", buffer ) )
       //   processDeveloper(buffer);
       // else if( startswith("geodynbc", buffer ) )
       //   processGeodynbc(buffer);
       else if (!inputFile.eof() && m_myRank == 0)
       {
	 // Maybe just reached eof, don't want to echo
	 // the ignoring command line for nothing
	 cout << "*** Ignoring command: '" << buffer << "'" << endl;
       }
     } // end if strlen(buffer) > 0
     
  } // end while !inputFile.eof() 
  
  if (m_myRank == 0)
     cout << endl;

  inputFile.close();

// wait until all processes have read the input file
  MPI_Barrier(MPI_COMM_WORLD);


  print_execution_time( time_start, MPI_Wtime(), "reading input file" );

  // ---------------------------------------------
  // cross command line checks
  // ---------------------------------------------
  if( mTopoImageFound && !m_topography_exists)
  {
    if (m_myRank == 0)
      cerr << "Error:  The input file is requesting a topo image but there is no topography command" << endl;
    return false;
  }
// if we made it this far, the object should be ready for time stepping
  return true;
}


void EW::processGrid(char* buffer)
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  int nx=0, ny=0, nz=0;
  double h = 0.0;

  double spherex=0, spherey=0, spherez=0, spherer=1.0;
  
  //-----------------------------------------------------------------
  // default geographical coordinates will be the 
  // nevada test site (see:  en.wikipedia.org/wiki/Nevada_Test_Site
  //-----------------------------------------------------------------
  double lat, lon;
  bool latSet = false, lonSet = false;
// default azimuth
  mGeoAz=0;
  
  //------------------------------------------------------------------
  // mesh refinement info
  //------------------------------------------------------------------
  int refineLevel = 0;
  // used to help with load balance issues
  int minBoxLength = 4;
  int maxBoxLength = 50;
  float psratio = 3.0;

  char* token = strtok(buffer, " \t");

  REQUIRE2(strcmp("grid", token) == 0, "ERROR: not a grid...: " << token);
  token = strtok(NULL, " \t");

  string err = "Grid Error: ";


  stringstream gridSetupErrStream;
  gridSetupErrStream << endl
		     << "----------------------------------------" << endl
		     << " Only five ways to setup grid: " << endl
		     << "  1. provide h and nx, ny, nz " << endl
		     << "  2. provide h and x, y, z " << endl
		     << "  3. provide x,y,z and nx " << endl
		     << "  4. provide x,y,z and ny " << endl
		     << "  5. provide x,y,z and nz " << endl
		     << "----------------------------------------" << endl
		     << endl;

  string gridSetupErr = gridSetupErrStream.str();

  if (m_myRank == 0)
    cout << endl << "* Processing the grid command..." << endl;

  while (token != NULL)
  {
     // while there are tokens in the string still
     if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
        break;
     if (startswith("ny=", token))
     {
	token += 3; // skip ny=

	CHECK_INPUT(atoi(token) > 0, 
		err << "ny is not a positive integer: " << token);
        ny = atoi(token);
     }
     else if (startswith("nx=", token))
     {
        token += 3; // skip nx=
           
	CHECK_INPUT(atoi(token) > 0, 
		err << "nx is not a positive integer: " << token);
        nx = atoi(token);
     }
     else if (startswith("nz=", token))
     {
        token += 3; // skip nz=
        
	CHECK_INPUT(atoi(token) >= 0, 
		err << "nz is not a positive integer: " << token);
        nz = atoi(token);
     }
     else if (startswith("x=", token))
     {
        token += 2; // skip x=
	CHECK_INPUT(atof(token) > 0.0, err << "x is not a positive float: " << token);
        x = atof(token);
     }
     else if (startswith("y=", token))
     {
      token += 2; // skip y=
      CHECK_INPUT(atof(token) >= 0.0, err << "y is negative: " << token);
      y = atof(token);
     }
     else if (startswith("z=", token))
     {
        token += 2; // skip z=
	CHECK_INPUT(atof(token) > 0.0, err << "z is not a positive float: " << token);
        z = atof(token);
     }
     else if (startswith("h=", token))
     {
       token += 2; // skip h=
       CHECK_INPUT(atof(token) > 0.0, 
 	       err << "h is not a positive float: " << token);
        h = atof(token);
     }
     else if (startswith("az=", token))
     {
        token += 3; // skip az=
        mGeoAz = atof(token);
        CHECK_INPUT(mGeoAz >= 0.0,
                err << "az must be greater than or equal to zero degrees, not: " << mGeoAz);
        CHECK_INPUT(mGeoAz <= 360.0,
                err << "az must be less than or equal to 360 degrees, not " << mGeoAz);
     }
     else if (startswith("lat=", token))
     {
        token += 4;
        lat = atof(token);
        CHECK_INPUT(lat >= -90.0,
                err << "lat must be greater than or equal to -90 degrees, not " 
                << lat);
        CHECK_INPUT(lat <= 90.0,
                err << "lat must be less than or equal to 90 degrees, not "
                << lat);
        latSet = true;
     }
     else if (startswith("lon=", token))
     {
        token += 4;
        lon = atof(token);
        CHECK_INPUT(lon >= -180.0,
                err << "lon must be greater or equal to -180 degrees, not " 
                << lon);
        CHECK_INPUT(lon <= 180.0,
                err << "lon must be less than or equal to 180 degrees, not "
                << lon);
        lonSet = true;
     }
//                        1234567890123456  
     else if (startswith("extrapolate=", token))
     {
        token += 12;
        int extrapolate = atoi(token);
        CHECK_INPUT(extrapolate >= 0 && extrapolate <= 2,
                err << "extrapolate must be an integer between 0 and 2, not " 
                << extrapolate);
	mMaterialExtrapolate = extrapolate;
     }
     else
     {
        badOption("grid", token);
     }
     token = strtok(NULL, " \t");
  }
  
  //--------------------------------------------------------------------
  // There are only three ways to specify a grid.
  //--------------------------------------------------------------------
  if (h != 0.0)
  {
    if (nx > 0 || nz > 0 || ny > 0)
    {
      //----------------------------------------------------------------
      // 1.  nx, [ny], nz and h
      //----------------------------------------------------------------
      CHECK_INPUT(nx && nz, gridSetupErr);
      CHECK_INPUT(x == 0.0 && y == 0.0 && z == 0.0, gridSetupErr);
    }
    else
    {
      //--------------------------------------------------------------
      // 2.  x, [y], z and h
      //--------------------------------------------------------------
      CHECK_INPUT(x > 0.0 && z > 0.0, gridSetupErr);
      CHECK_INPUT(nx == 0 && ny == 0 && nz == 0, gridSetupErr);
    }
  }
  else
  {
    //--------------------------------------------------------------------
    // 3.  x, [y], z and nx|ny|nz
    //--------------------------------------------------------------------
    CHECK_INPUT(x > 0.0 && z > 0.0, gridSetupErr);
    CHECK_INPUT((nx > 0) + (ny > 0) + (nz > 0) == 1, gridSetupErr);
  }
  
  int nxprime, nyprime, nzprime;
  // -------------------------------------------------------------
  // Make sure all the bounds are consistent.
  //
  // In order to divide up the space properly, we must take the 
  // coordinate dimension, say x, and divide by the number of grid
  // points requested minus one.  This is because we'd like to 
  // include the end points in the spatial data arrays.  For example,
  // if x = 1000. and nx = 10, you'd think h would be 100.  However,
  // if we'd like the bounds to go from -500 to 500, we actually 
  // need x divided up into 9 cells so that both x = -500 and x = 500
  // will be included in the data array.  In this case, h would be
  // 111.11, giving:
  //
  // x[1] = -500, x[2] = -388 x[3] = -277 x[4] = -166 x[5] = -55
  // x[6] = 55 x[7] = 166 x[8] = 277 x[9] = 388 x[10] = 500.
  // -------------------------------------------------------------

  if (latSet && lonSet)
  {
     mLatOrigin = lat;
     mLonOrigin = lon;
  }
  else
  {
// Default is NTS
     mLatOrigin = 37.0;
     mLonOrigin = 118.0;
  }

  double cubelen, zcubelen;
//   if( m_geodynbc_found )  {
// // Set WPP grid spacing based on Geodyn cube data

//      double origin[3]={0,0,0}, ibclat, ibclon, ibcaz;

//      bool found_latlon;
//      int adjust;
//      geodynbcGetSizes( m_geodynbc_filename, origin, cubelen, zcubelen, found_latlon, ibclat,
// 		       ibclon, ibcaz, adjust );
//      // Use approximate h
//      if( h == 0.0 )
//      {
// 	if( nx > 0 )
// 	   h = x/(nx-1);
// 	else if( nz > 0 )
// 	   h = z/(nz-1);
// 	else
// 	   h = y/(ny-1);
//      }

//      // rounding of cube position to two decimals (prec=100), three (prec=1000) etc..
//      double prec = 100;

//      if( found_latlon )
//      {
//         CHECK_INPUT( fabs(ibcaz - mGeoAz) < 1e-5, "Error: Az in Geodyn file, "
// 		 << ibcaz << " is different from Az in WPP, " << mGeoAz );
	   
// 	// lat-lon corner of cube given
// 	if( adjust == 1 || origin[2] == 0 )
// 	{
// 	   // h based on cube length only, adjust z-position of cube
// 	   int nc = static_cast<int>(round(cubelen)/h);
// 	   h = cubelen/nc;
// 	   origin[2] -= h*( origin[2]/h-round(origin[2]/h) );
// 	}
//         else
// 	{
// 	   // h based on cube length and z-position of cube
//            int a = static_cast<int>(round(origin[2]*prec));
// 	   int b = static_cast<int>(round((origin[2]+zcubelen)*prec));
// 	   // 
//            int d  = gcd(a,b);
// 	   int n1 = a/d;
// 	   int k  = static_cast<int>(round(origin[2]/(n1*h)));
//            h = origin[2]/(k*n1);
// 	}
// 	// Geographic origin adjustment:
//         double gridLat = mLatOrigin;
// 	   double gridLon = mLonOrigin;
//         double metersPerDegree = mMetersPerDegree;
//         double deg2rad = M_PI/180;
//         double phi = mGeoAz*deg2rad;
// 	double x = metersPerDegree*( cos(phi)*(ibclat-gridLat) + cos(ibclat*deg2rad)*(ibclon-gridLon)*sin(phi));
// 	double y = metersPerDegree*(-sin(phi)*(ibclat-gridLat) + cos(ibclat*deg2rad)*(ibclon-gridLon)*cos(phi));
//         x -= h*(x/h-round(x/h));
//         y -= h*(y/h-round(y/h));
//         gridLat = ibclat - (x*cos(phi) - y*sin(phi))/metersPerDegree;
// 	gridLon = ibclon - (x*sin(phi) + y*cos(phi))/(metersPerDegree*cos(ibclat*deg2rad));
     // mLatOrigin = gridLat;
     // mLonOrigin = gridLon;
//         origin[0] = x;
// 	origin[1] = y;
//      }     
//      else
//      {
// 	// lat-lon corner of cube not given, interpret origin realtive (0,0,0)
//         if( m_geodynbc_center )
// 	{
// 	   // Center cube in the middle of the domain (in x,y), discarding input origin.
// 	   double xlen = x;
// 	   double ylen = y;
// 	   if( xlen == 0 )
// 	      xlen = h*(nx-1);
// 	   if( ylen == 0 )
// 	      ylen = h*(ny-1);
// 	   origin[0] = 0.5*(xlen-cubelen);
// 	   origin[1] = 0.5*(ylen-cubelen);
// 	}
// 	if( adjust == 1 )
// 	{
// 	   // h based on cube length only, adjust cube position
// 	   int nc = static_cast<int>(round(cubelen/h));
// 	   h = cubelen/nc;
// 	   origin[0] -= h*( origin[0]/h-round(origin[0]/h) );
// 	   origin[1] -= h*( origin[1]/h-round(origin[1]/h) );
// 	   origin[2] -= h*( origin[2]/h-round(origin[2]/h) );
// 	}
// 	else
// 	{
// 	   // h based on cube length and cube position, might be very restrictive
// 	   CHECK_INPUT( false, "Error: cube position without lat/long position must be adjustable");
// 	}
//      }
     
     
//      if (nx == 0 && x != 0.0)
// 	nxprime = computeEndGridPoint(x, h);
//      else if (nx != 0)
// 	nxprime = nx;
//      else
// 	CHECK_INPUT(0, gridSetupErr);

//      if (nz == 0 && z != 0.0)
// 	nzprime = computeEndGridPoint(z, h);
//      else if (nz != 0)
// 	nzprime = nz;
//      else
// 	CHECK_INPUT(0, gridSetupErr);

//      if (ny == 0 && y != 0.0)
// 	nyprime = computeEndGridPoint(y, h);
//      else if (ny != 0)
// 	nyprime = ny;
//      else
// 	CHECK_INPUT(0, gridSetupErr);
//      m_ibc_origin[0] = origin[0];
//      m_ibc_origin[1] = origin[1];
//      m_ibc_origin[2] = origin[2];
//      //     cout << "Cube origin " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
//      //     cout << "Cube length " << cubelen << endl;
//      //     cout << "nx,ny,nz " << nxprime << " " << nyprime << " " << nzprime << endl;     
//   } // end if m_geodynbc_found
//   else
  {
     if (nx > 0 && h == 0.0)
     {
    // we set the number grid points in the x direction
    // so we'll compute the grid spacing from that.
	h = x / (nx-1);
	if (m_myRank == 0)
	   cout << "* Setting h to " << h << " from  x/(nx-1) (x=" << x << ", nx=" << nx << ")" << endl;
      
	nxprime = nx;
	nzprime = computeEndGridPoint(z, h);
	nyprime = computeEndGridPoint(y, h);
     }
     else if (nz > 0 && h == 0.0)
     {
    // set the number of grid points from z direction and nz
	h = z/(nz-1);
	if (m_myRank == 0)
	   cout << "* Setting h to " << h << " from  z/(nz-1) (z=" << z << ", nz=" << nz << ")" << endl;
	nzprime = nz;
	nxprime = computeEndGridPoint(x, h);
	nyprime = computeEndGridPoint(y, h);
     }
     else if (ny > 0 && h == 0.0)
     {
    // set hte number of grid points from y direction and ny
	h = y/(ny-1);
	if (m_myRank == 0)
	   cout << "* Setting h to " << h << " from  y/(ny-1) (y=" << y << ", ny=" << ny << ")" << endl;
	nyprime = ny;
	nxprime = computeEndGridPoint(x, h);
	nzprime = computeEndGridPoint(z, h);
     }
     else
     {
	//----------------------------------------------------
	// h was set by the user, so compute the appropriate
	// nx, ny, and nz or x, y, z.
	//----------------------------------------------------
	if (nx == 0 && x != 0.0)
	   nxprime = computeEndGridPoint(x, h);
	else if (nx != 0)
	   nxprime = nx;
	else
	   CHECK_INPUT(0, gridSetupErr);

	if (nz == 0 && z != 0.0)
	   nzprime = computeEndGridPoint(z, h);
	else if (nz != 0)
	   nzprime = nz;
	else
	   CHECK_INPUT(0, gridSetupErr);

	if (ny == 0 && y != 0.0)
	   nyprime = computeEndGridPoint(y, h);
	else if (ny != 0)
	   nyprime = ny;
	else
	   CHECK_INPUT(0, gridSetupErr);
     }
  }
  
  
  if (nxprime != nx && m_myRank == 0)
    cout << "* Setting nx to " << nxprime << " to be consistent with h=" << h << endl;
  if (nyprime != ny && m_myRank == 0)
    cout << "* Setting ny to " << nyprime << " to be consistent with h=" << h << endl;
  if (nzprime != nz && m_myRank == 0)
    cout << "* Setting nz to " << nzprime << " to be consistent with h=" << h << endl;

  // -------------------------------------------------------------
  // Now we adjust the geometry bounds based on the actual 
  // number of grid points used in each dimension.
  // -------------------------------------------------------------
  double xprime, yprime, zprime;
  xprime = (nxprime-1)*h;
  zprime = (nzprime-1)*h;
  yprime = (nyprime-1)*h;
  
  double eps = 1.e-9*sqrt(SQR(xprime)+SQR(yprime)+SQR(zprime));
  
  if (fabs(xprime-x) > eps && m_myRank == 0)
    cout << "* Changing x from " << x << " to " << xprime << " to be consistent with h=" << h << endl;
  if (fabs(zprime-z) > eps && m_myRank == 0)
    cout << "* Changing z from " << z << " to " << zprime << " to be consistent with h=" << h << endl;
  if (fabs(yprime-y) > eps && m_myRank == 0)
    cout << "* Changing y from " << y << " to " << yprime << " to be consistent with h=" << h << endl;

  // if( m_geodynbc_found )
  // {
  //    CHECK_INPUT( m_ibc_origin[0]>0 && m_ibc_origin[0]+cubelen<xprime , "Error: Cube x-dimension [" 
  // 	      << m_ibc_origin[0] << "," << m_ibc_origin[0]+cubelen << 
  // 	      "] not inside domain of length "<< xprime );
  //    CHECK_INPUT( m_ibc_origin[1]>0 && m_ibc_origin[1]+cubelen<yprime , "Error: Cube y-dimension ["
  // 	      << m_ibc_origin[1] << "," << m_ibc_origin[1]+cubelen << 
  // 	      "] not inside domain of length "<< yprime );
  //    CHECK_INPUT( m_ibc_origin[2]>=0 && m_ibc_origin[2]+zcubelen<zprime , "Error: Cube z-dimension ["
  // 	      << m_ibc_origin[2] << "," << m_ibc_origin[2]+zcubelen << 
  // 	      "] not inside domain of length "<< zprime );
  // }

  if (!(latSet && lonSet) && (latSet || lonSet))
  {
    stringstream msg;
    if (m_myRank == 0)
    {
      msg << " \n* Improper grid location specification, must specify both lat and lon variables " << endl
	  << " * Missing... ";
      if (!latSet)
	msg << " lat= keyword ";
      if (!lonSet)
	msg << " lon= keyword ";
    }
    CHECK_INPUT(0, msg.str());
  }

  m_nx_base = nxprime;
  m_ny_base = nyprime;
  m_nz_base = nzprime;
  m_h_base = h;
  m_global_xmax = xprime;
  m_global_ymax = yprime;
  m_global_zmax = zprime;

}

//----------------------------------------------------------
void EW::cleanUpRefinementLevels()
{
  CHECK_INPUT(m_topo_zmax < m_global_zmax-m_h_base,"The topography is extending too deep into the ground and there is no space for the Cartesian grid.");

// Add a top zMin level
// Here zMin = m_topo_zmax if m_topography_exists, otherwise zMin = 0;
  double zMin;
  
  if (m_topography_exists)
  {
    m_refinementBoundaries.push_back(m_topo_zmax);
    zMin = m_topo_zmax;
  }
  else
  {
    m_refinementBoundaries.push_back(0.0);
    zMin = 0.;
  }

// need to sort m_refinementBoundaries in decreasing order
  int nRef = m_refinementBoundaries.size();
  double *zValues = new double[nRef];
  int q;
// tmp
//   cout << "Original order"<<endl;
//   for (q=0; q<nRef; q++)
//     cout<< m_refinementBoundaries[q] << endl;

  for (q=0; q<nRef; q++)
    zValues[q] = m_refinementBoundaries[q];
  
  sort(zValues, zValues+nRef);
  
// reverse the ordering to get decreasing order
  for (q=0; q<nRef; q++)
    m_refinementBoundaries[q] = zValues[nRef-q-1];
  
// tmp
//   cout << "Sorted order"<<endl;
//   for (q=0; q<nRef; q++)
//     cout<< m_refinementBoundaries[q] << endl;

// cleanup
  delete [] zValues;

// remove any refinement levels outside zMin < z < m_zmax. 
// declare an iterator
  vector<double>::iterator it;
// tmp
//  cout << "Removing items outside the range zMin = " << zMin << " < z < " << " zMax=" << m_zmax << "..." << endl;

  for (it=m_refinementBoundaries.begin(); it!=m_refinementBoundaries.end(); it++)
  {
    if (*it < zMin || *it >= m_global_zmax)
    {
// remove this entry from the vector
//      cout << "Removing out-of-range refinement level="<< *it << endl;
      it = m_refinementBoundaries.erase(it); // returns next element
// need to back up one step to undo the it++ that always happens at the end of the for loop
      it--;
    }
  }
  
// need to remove any duplicate entries in the m_refinementBoundaries array
// tmp
//  cout << "Removing duplicate items..."<< endl;
  double z0 = m_refinementBoundaries[0]; // first item
  for (it=++m_refinementBoundaries.begin(); it!=m_refinementBoundaries.end(); it++)
  {
    if (*it == z0)
    {
// remove this entry from the vector
//      cout << "Removing duplicate refinement level="<< *it << endl;
      it = m_refinementBoundaries.erase(it); // returns next element
// need to back up one step to undo the it++ that always happens at the end of the for loop
      it--;
    }
    z0 = *it; // remember the current level
  }

// tmp
//   if (m_myRank==0)
//   {
//     cout << "Input refinement levels (z=):"<<endl;
//     for (it=m_refinementBoundaries.begin(); it!=m_refinementBoundaries.end(); it++)
//       cout<< *it << endl;
//   }
  

}


// //-----------------------------------------------------------------------
// void FileInput::processRefinement(char* buffer)
// {
//   char* token = strtok(buffer, " \t");
//   CHECK_INPUT(strcmp("refinement", token) == 0, 
// 	      "ERROR: not a refinement line...: " << token);
//   token = strtok(NULL, " \t");

//   string err = "Refinement error ";

//   while (token != NULL)
//   {
//     // while there are tokens in the string still
//     if (startswith("#", token) || startswith(" ", buffer))
//       // Ignore commented lines and lines with just a space.
//       break;
//     else if( startswith("zmax=", token) )
//     {
//       token += 5; // skip zmax=
//       double z1 = atof(token);
//       m_refinementBoundaries.push_back(z1);
// //       if (m_myRank==0)
// // 	cout <<"Adding refinement boundary at z=" << z1 << endl;
//     }
//     else
//     {
//       badOption("refinement", token);
//     }
//     token = strtok(NULL, " \t");
//   }
  
// //  mSimulation->add_refinement_block( refb );

// }

// //-----------------------------------------------------------------------
// void FileInput::processAttenuation(char* buffer)
// {
//   char* token = strtok(buffer, " \t");
//   CHECK_INPUT(strcmp("attenuation", token) == 0, 
// 	      "ERROR: not a attenuation line...: " << token);
//   token = strtok(NULL, " \t");

//   string err = "Attenuation error ";
//   int nmech=3;
//   double velofreq=1;
//   bool foundppw = false, foundfreq=false;
// // Default is max frequency 2 Hz, 
//   m_att_ppw = -1;
//   m_att_max_frequency = 2.0;
  
//   while (token != NULL)
//   {
//     // while there are tokens in the string still
//     if (startswith("#", token) || startswith(" ", buffer))
//       // Ignore commented lines and lines with just a space.
//       break;
// //                       123456
//     else if( startswith("nmech=", token) )
//     {
//       token += 6; // skip nmech=
//       nmech = atoi(token);
//       CHECK_INPUT(nmech > 0 && nmech <= 8, "ERROR: Number of attenuation mechanisms must be > 0 and <= 8, not " << nmech);
//     }
// //                       1234567890123
//     else if( startswith("phasefreq=", token) )
//     {
//       token += 10; // skip phasefreq=
//       velofreq = atof(token);
//       CHECK_INPUT(velofreq >= 0 && velofreq <= 1000, "ERROR: Velocity frequency must be >= 0 and <= 1000 [Hz], not " << velofreq);
//     }
//     else if( startswith("maxfreq=",token) )
//     {
//        token += 8;
//        m_att_ppw = -1;
//        m_att_max_frequency = atof( token );
//        foundfreq = true;
//        CHECK_INPUT(m_att_max_frequency >= 0,"ERROR: maximum frequency must be >= 0, not " << m_att_max_frequency);
//     }
//     else if( startswith("minppw=",token) )
//     {
//        token += 7;
//        m_att_ppw = atoi( token );
//        foundppw = true;
//        CHECK_INPUT(m_att_ppw >= 0, "ERROR: minimum ppw must be >= 0, not " << m_att_ppw);
//     }
//     else
//     {
//       badOption("attenuation", token);
//     }
//     token = strtok(NULL, " \t");
//   }
//   if( foundppw && foundfreq )
//   {
//      if (m_myRank == 0)
// 	cout << "ERROR: Can not give both minppw and maxfreq for attenuation " << endl;
//      MPI_Abort(MPI_COMM_WORLD, 1);
//   }

//   m_nmech = nmech;
//   m_velo_omega = velofreq*2*M_PI;
//   m_use_attenuation=true;
  
// // tmp
//   if (m_myRank==0)
//     printf("* Processing the attenuation command: m_nmech=%i, m_velo_omega=%e\n", m_nmech, m_velo_omega);
// }

// void FileInput::processTopography(char* buffer)
// {
//    char* token = strtok(buffer, " \t");
//    CHECK_INPUT(strcmp("topography", token) == 0, 
// 	    "ERROR: not a topography line...: " << token);
//    string topoFile="surf.tp", style, fileName;
//    bool needFileName=false, gotFileName=false;

//    token = strtok(NULL, " \t");

//    while (token != NULL)
//    {
//      // while there are still tokens in the string 
//      if (startswith("#", token) || startswith(" ", buffer))
//        // Ignore commented lines and lines with just a space.
//        break;
//      if (startswith("zmax=", token))
//      {
//        token += 5; // skip logfile=
//        m_topo_zmax = atof(token);
// //        if (m_myRank==0)
// // 	 cout << "Setting topo zmax="<<m_topo_zmax<<endl;
//      }
// //                       1234567890
//      else if (startswith("order=", token))
//      {
//        token += 6; // skip logfile=
//        m_grid_interpolation_order = atoi(token);
//        if (m_grid_interpolation_order < 2 || m_grid_interpolation_order > 4)
//        {
// 	 if (m_myRank == 0)
// 	   cout << "order needs to be 2,3, or 4, not: " << m_grid_interpolation_order << endl;
// 	 MPI_Abort(MPI_COMM_WORLD, 1);
//        }
       
// //        if (m_myRank==0)
// // 	 cout << "Setting interpolation order to=" << m_grid_interpolation_order << endl;
//      }
// //                        123456789
//      else if (startswith("smooth=", token))
//      {
//        token += 7; // skip smooth=
//        m_maxIter = atoi(token);
//        if (m_maxIter < 0 || m_maxIter > 1000)
//        {
// 	 if (m_myRank == 0)
// 	   cout << "Number of smoothing iterations needs to be >=0 and <=1000, not: "<< m_maxIter << endl;
// 	 MPI_Abort(MPI_COMM_WORLD, 1);
//        }
//      }
//      else if( startswith("input=", token ) )
//      {
//        token += 6;
//        style = token;
//        if (strcmp("grid", token) == 0)
//        {
// 	 m_topoInputStyle=GridFile;
// 	 m_topography_exists=true;
// 	 needFileName=true;
//        }
//        else if (strcmp("cartesian", token) == 0)
//        {
// 	 m_topoInputStyle=CartesianGrid;
// 	 m_topography_exists=true;
// 	 needFileName=true;
//        }
//        else if (strcmp("efile", token) == 0)
//        {
// 	 m_topoInputStyle=Efile;
// 	 m_topography_exists=true;
//        }
//        else if (strcmp("gaussian", token) == 0)
//        {
// 	 m_topoInputStyle=GaussianHill;
// 	 m_topography_exists=true;
//        }
//        else
//        {
// 	 badOption("topography> input", token);
//        }
//      }
//      else if( startswith("file=", token ) )
//      {
//        token += 5;
//        m_topoFileName = token;
//        gotFileName=true;
// //        if (m_myRank==0)
// // 	 cout << "read topo file name=" << m_topoFileName <<endl;
//      }
// //                        12345678901
//      else if( startswith("resolution=", token ) )
//      {
//        token += 11;
//        m_EFileResolution = atof(token);
//        CHECK_INPUT(m_EFileResolution>0.,"Resolution must be positive, not " << m_EFileResolution);
//      }
// //                        123456789012
//      else if( startswith("gaussianAmp=", token ) )
//      {
//        token += 12;
//        m_GaussianAmp = atof(token);
//      }
// //                        123456789012
//      else if( startswith("gaussianXc=", token ) )
//      {
//        token += 11;
//        m_GaussianXc = atof(token);
//      }
// //                        123456789012
//      else if( startswith("gaussianYc=", token ) )
//      {
//        token += 11;
//        m_GaussianYc = atof(token);
//      }
// //                        123456789012
//      else if( startswith("gaussianLx=", token ) )
//      {
//        token += 11;
//        m_GaussianLx = atof(token);
//      }
// //                        123456789012
//      else if( startswith("gaussianLy=", token ) )
//      {
//        token += 11;
//        m_GaussianLy = atof(token);
//      }

//      else
//      {
//        badOption("topography", token);
//      }
//      token = strtok(NULL, " \t");
//    }
//    if (needFileName)
//      CHECK_INPUT(gotFileName, 
// 	      "ERROR: no topography file name specified...: " << token);

// }

// void FileInput::processDamping(char* buffer)
// {
//   char* token = strtok(buffer, " \t");
//   CHECK_INPUT(strcmp("damping", token) == 0, 
// 	   "ERROR: not a damping line...: " << token);
//   token = strtok(NULL, " \t");

// // default: no damping

// //  Coefficients are in precentage of max allowed by CFL constraint
//   double d4cof= 0, curlcof=0, atacof=0;
//   bool d4set=false, curlset=false;
//   string err = "damping error ";

//   while (token != NULL)
//     {
//       // while there are tokens in the string still
//        if (startswith("#", token) || startswith(" ", buffer))
//         // Ignore commented lines and lines with just a space.
//         break;
// //        else if (startswith("ata=", token))
// //        {
// //           token += 4; // skip ata=
// //           atacof = atof(token);
// //        }
// //        else if (startswith("curlcurl=", token))
// //        {
// //           token += 9; // skip curlcurl=
// //           curlcof = atof(token);
// //           curlset = true;
// //        }
//        else if( startswith("d4=",token) )
//        {
//           token += 3;
// 	  d4cof = atof(token);
// 	  d4set = true;
//        }
//        else
//        {
//           badOption("damping", token);
//        }
//        token = strtok(NULL, " \t");
//     }
  
// //   if( curlset )
// //      mSimulation->turnOnCurlCurlDamping( curlcof );
  
// //   if( atacof != 0 )
// //      mSimulation->turnOnATADamping( atacof );

//   if( d4set )
//      mSimulation->setDampingCFL( d4cof );
// }

// //-----------------------------------------------------------------------
// void FileInput::processEnergy(char* buffer)
// {
//     char* token = strtok(buffer, " \t");
//     CHECK_INPUT(strcmp("energy", token) == 0, 
//  	    "ERROR: not a energy test line...: " << token);
//     int seed;
//     string logfile="energy.dat";
//     double cpcsratio=3;
//     bool print=false;
//     bool const_coeff = false;

//     token = strtok(NULL, " \t");

//     while (token != NULL)
//     {
//        // while there are tokens in the string still
//         if (startswith("#", token) || startswith(" ", buffer))
//            // Ignore commented lines and lines with just a space.
//            break;
//         if (startswith("logfile=", token))
//         {
//            token += 8; // skip logfile=
//            logfile = token;
//         }
//         else if( startswith("seed=", token ) )
//         {
//  	  token += 5;
//  	  seed = atoi(token);
//         }
//         else if( startswith("cpcsratio=", token ) )
//         {
//  	  token += 10;
//  	  cpcsratio = atof(token);
//         }
//         else if( startswith("constant_coeff=", token ) )
//         {
//  	  token += 15;
//  	  const_coeff = atoi(token) == 1;
//         }
//         else if( startswith("print=", token ) )
//         {
//  	  token += 6;
//  	  print = atoi(token) == 1;
//         }
//         else
//         {
//            badOption("energy", token);
//         }
//         token = strtok(NULL, " \t");
//     }
//     Forcing* force = new ForcingEnergy( seed, cpcsratio, const_coeff );
//     mSimulation->set_forcing( force );
//     mSimulation->set_energylog( logfile, print, true );
// }

//-----------------------------------------------------------------------
void EW::processTwilight(char* buffer)
{
   //  if (m_myRank == 0)
   //    cout << "Entering twilight mode..." << endl;

  double omega = 1.;
  double momega= 1.;
  double phase = 0;
  double mphase= 0.4;
  double c = 1.3;
  double amprho = 1;
  double ampmu  = 1;
  double amplambda = 1;

  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("twilight", token) == 0, "ERROR: not a twilight line...: " << token);
  token = strtok(NULL, " \t");

  string err = "Twilight Error: ";

  while (token != NULL)
    {
      // while there are tokens in the string still
       if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
       if( startswith("errorlog=",token) )
       {
          token += 9;
	  bool errorlog = atoi(token)==1;
	  if( errorlog )
	     switch_on_error_log();
       }
       else if( startswith("omega=",token) )
       {
          token += 6;
	  omega = atof(token);
       }
       else if( startswith("c=",token) )
       {
          token += 2;
	  c = atof(token);
       }
       else if( startswith("phase=",token) )
       {
          token += 6;
	  phase = atof(token);
       }
       else if( startswith("momega=",token) )
       {
          token += 7;
	  momega = atof(token);
       }
       else if( startswith("mphase=",token) )
       {
          token += 7;
	  mphase = atof(token);
       }
       else if( startswith("amprho=",token) )
       {
          token += 7;
	  amprho = atof(token);
       }
       else if( startswith("ampmu=",token) )
       {
          token += 6;
	  ampmu = atof(token);
       }
       else if( startswith("amplambda=",token) )
       {
          token += 10;
	  amplambda = atof(token);
       }
       else
       {
          badOption("twilight", token);
       }
       token = strtok(NULL, " \t");
    }

  ForcingTwilight* forcing;
  forcing = new ForcingTwilight( omega, c, phase, momega, mphase, amprho, ampmu, amplambda );

  set_twilight_forcing( forcing );
}

// void FileInput::processDeveloper(char* buffer)
// {
//    //   if (m_myRank == 0)
//    //      cout << "Entering developer mode..." << endl;

//   int ilno = 5;
//   int update_boundary_function = 1;
//   double cfl=-1;
//   bool cflset = false;
//   bool output_load = false;
//   bool output_timing = false;
//   bool use_alltoallv = true;
//   bool logenergy = false;
//   bool printenergy = false;
//   string energyfile = "energy.dat";
//   bool use_mpiio = false;
//   bool use_iotiming = false;

//   bool cons = true;
//   double ctol = 1e-3;
//   int cmaxit = 20;
//   char* token = strtok(buffer, " \t");
//   CHECK_INPUT(strcmp("developer", token) == 0, "ERROR: not a developer line...: " << token);
//   token = strtok(NULL, " \t");
//   while (token != NULL)
//   {
//     // while there are tokens in the string still
//     if (startswith("#", token) || startswith(" ", buffer))
//       // Ignore commented lines and lines with just a space.
//       break;
// //     if (startswith("inner_loop=", token))
// //     {
// //       token += 11; // skip name=
// //       ilno = atoi(token);
// //     }
//     if (startswith("cfl_number=", token))
//     {
//       token += 11; // skip name=
//       cfl = atof(token);
//       CHECK_INPUT( cfl > 0, "Error negative CFL number");
//       cflset = true;
//     }
// //     if (startswith("update_processor_boundary=", token))
// //     {
// //       token += 26;

// //       if( strcmp(token,"wpp_buffer") == 0 )
// // 	update_boundary_function = 2;
// //       else if( strcmp(token,"mpi_datastructure")==0 )
// // 	update_boundary_function = 1;
// //       else
// // 	CHECK_INPUT( false, "token x" << token << "x not valid for update_processor_boundary " <<
// // 		 "should be `wpp_buffer' or `mpi_datastructure'");
// //     }
//     if (startswith("output_load=", token))
//     {
//       token += 12;
//       output_load = (atoi(token) == 1);
//     }
//     else if (startswith("output_timing=", token))
//     {
//       token += 14;
//       output_timing = (atoi(token) == 1);
//     }
//     else if (startswith("interpolation=", token))
//     {
//       token += 14;
//       cons = strcmp(token,"conservative") == 0;
//     }
//     else if (startswith("ctol=", token))
//     {
//       token += 5;
//       ctol = atof(token);
//     }
//     else if (startswith("cmaxit=", token))
//     {
//       token += 7;
//       cmaxit = atoi(token);
//     }
//     else if (startswith("log_energy=", token))
//     {
//        logenergy = true;
//        token += 11;
//        energyfile=token;
//     }
//     else if (startswith("print_energy=", token))
//     {
//        token += 13;
//        printenergy = (atoi(token) == 1);
//     }
// //                       123456789
//        else if (startswith("mpiio=", token))
//        {
// 	 token += 6;
// 	 use_mpiio = (atoi(token) == 1);
//        }
// //                          123456789
//        else if (startswith("iotiming=", token))
//        {
// 	 token += 9;
// 	 use_iotiming = (atoi(token) == 1);
//        }
// //     if( startswith("use_alltoallv=", token ) )
// //     {
// //       token += 14;
// //       use_alltoallv = (atoi(token) == 1);
// //     }
//     else
//     {
//       badOption("developer", token);
//     }
//     token = strtok(NULL, " \t");
//   }
// //   if( ilno != 5 )
// //     mSimulation->set_inner_loop( ilno );
//   if( cflset )
//     mSimulation->set_cflnumber( cfl );
// //   mSimulation->set_update_boundary_function( update_boundary_function );
//    mSimulation->set_output_options( output_load, output_timing );
// //  mSimulation->set_alltoallv( use_alltoallv );
//    mSimulation->set_conservative_interpolation( cons, ctol, cmaxit );
//    if( logenergy || printenergy )
//       mSimulation->set_energylog( energyfile, printenergy, logenergy );
//   mSimulation->setIO_method(use_mpiio, use_iotiming);
// }

//-----------------------------------------------------------------------
 void EW::processTestPointSource(char* buffer)
 {
    char* token = strtok(buffer, " \t");
    CHECK_INPUT(strcmp("testpointsource", token) == 0, "ERROR: not a testpointsource line...: " << token);
    token = strtok(NULL, " \t");
    double cs = 1.0, rho=1.0, cp=sqrt(3.0);
    while (token != NULL)
    {
       if (startswith("#", token) || startswith(" ", buffer))
          break;

       if (startswith("cp=", token))
       {
          token += 3; 
          cp = atof(token);
       }
       else if (startswith("cs=", token))
       {
          token += 3; 
          cs = atof(token);
       }
       else if (startswith("rho=", token))
       {
          token += 4; 
          rho = atof(token);
       }
       else if (startswith("diractest=", token))
       {
          token += 10; 
          if( strcmp(token,"1")==0 || strcmp(token,"true")==0 )
	     m_moment_test = true;
       }
       else
       {
 	 badOption("testpointsource", token);
       }
       token = strtok(NULL, " \t");
    }
    m_point_source_test = new TestPointSource( rho, cs, cp );
 }

//-----------------------------------------------------------------------
void EW::processTestLamb(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("testlamb", token) == 0, "ERROR: not a testlamb line...: " << token);
   token = strtok(NULL, " \t");

   string err = "Testlamb Error: ";
   double x0=0.0, y0=0.0, z0=0.0;
   double cs = 1.0, rho=1.0, cp=sqrt(3.0), fz=1.0, freq=1.0, f0=1.0; // the exact solution assumes freq = 1

   while (token != NULL)
   {
      if (startswith("#", token) || startswith(" ", buffer))
         break;

      // if (startswith("x=", token))
      // {
      //    token += 2; // skip x=
      //    x0 = atof(token);
      // }
      // else if (startswith("y=", token))
      // {
      //    token += 2; // skip y=
      //    y0 = atof(token);
      // }
      if (startswith("cp=", token))
      {
         token += 3; 
         cp = atof(token);
      }
// exact solution assumes cs=cp/sqrt(3), so we will hard-wire this ratio below
//       else if (startswith("cs=", token))
//       {
//          token += 3; 
//          cs = atof(token);
//       }
      else if (startswith("rho=", token))
      {
         token += 4; 
         rho = atof(token);
      }
      // else if (startswith("fz=", token))
      // {
      //    token += 3; 
      //    fz = atof(token);
      // }
      else
      {
	 badOption("testlamb", token);
      }
      token = strtok(NULL, " \t");
   }
// point forcing is now set with the source command
   // double fx=0, fy=0, t0=0;
   // timeDep tdep = iVerySmoothBump;
   // Source* source = new Source( mSimulation, f0, freq, t0, x0, y0, z0, fx, fy, fz,
   // 				tdep, "lambsource", 0 );
   m_lamb_test = new TestLamb( rho, cp );
}

//-----------------------------------------------------------------------
void EW::processFileIO(char* buffer)
{
   int printcycle = 100;
   char* path = 0;
   char* scenario = 0;
   int nwriters=8;
   bool pfs=false;
   
   int verbose = 0;
   bool debug = false;

   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("fileio", token) == 0, "ERROR: not a fileio line...: " << token);
   token = strtok(NULL, " \t");

   string err = "FileIO Error: ";

  while (token != NULL)
    {
       if (startswith("#", token) || startswith(" ", buffer))
          break;
       if(startswith("path=", token)){
          token += 5; // skip path=
          path = token;
       }
//                          123456789
       else if (startswith("verbose=", token))
       {
          token += 8; // skip verbose=
          CHECK_INPUT(atoi(token) >= 0, err << "verbose must be non-negative, not: " << token);
          verbose = atoi(token);
       }
       else if (startswith("printcycle=", token))
       {
          token += 11; // skip printcycle=
          CHECK_INPUT(atoi(token) > -1,
	         err << "printcycle must be zero or greater, not: " << token);
          printcycle = atoi(token);
       }
       else if (startswith("pfs=", token))
       {
          token += 4; // skip pfs=
          pfs = (atoi(token) == 1);
       }
//                          1234567890
       else if (startswith("nwriters=", token))
       {
          token += 9; // skip nwriters=
          CHECK_INPUT(atoi(token) > 0,
	         err << "nwriters must be positive, not: " << token);
          nwriters = atoi(token);
       }
       else
       {
          badOption("fileio", token);
       }
       token = strtok(NULL, " \t");
    }

  if (path != 0) setOutputPath(path);
  setPrintCycle(printcycle);
  setVerbosity(verbose);
  setParallel_IO(pfs, nwriters);
}

// void
// FileInput::
// processGMT(char* buffer)
// {
//   string filename = "wpp.gmt.csh";
//   char* token = strtok(buffer, " \t");
//   CHECK_INPUT(strcmp("gmt", token) == 0, "ERROR: not a gmt line...: " << token);
//   token = strtok(NULL, " \t");

//   string err = "GMT Error: ";

//   while (token != NULL)
//     {
//       // while there are tokens in the string still
//       if (startswith("#", token) || startswith(" ", buffer))
// 	// Ignore commented lines and lines with just a space.
// 	break;
//       if (startswith("file=", token))
// 	{
//           token += 5; // skip file=
//           filename = token;
//        }
//       else
// 	{
//           badOption("gmt", token);
//        }
//       token = strtok(NULL, " \t");
//     }

//   mSimulation->setGMTOutput(filename, mFileName);

// }

void EW::processTime(char* buffer)
{
  double t=0.0;
  int steps = -1;
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("time", token) == 0, "ERROR: not a time line...: " << token);
  token = strtok(NULL, " \t");

  string err = "Time Error: ";

 while (token != NULL)
    {
      // while there are still tokens in the string
       if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
       if (startswith("t=", token))
       {
          // If t==0, just go one step (dt is going to be 0)
          token += 2; // skip t=
          CHECK_INPUT(atof(token) >= 0.0, err << "t is not a positive float: " << token);
          t = atof(token);
       }
       else if (startswith("steps=", token))
       {
          token += 6; // skip steps=
          CHECK_INPUT(atoi(token) >= 0, err << "steps is not a non-negative integer: " << token);
          steps = atoi(token);
       }
       else
       {
          badOption("time", token);
       }
       token = strtok(NULL, " \t");
    }
  CHECK_INPUT(!( (t > 0.0) && (steps >= 0) ),
          "Time Error: Cannot set both t and steps for time");
  
  if (t > 0.0)
    setGoalTime(t);
  else if (steps >= 0)
    setNumberSteps(steps);
}

void EW::processBoundaryConditions(char *buffer)
{
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("boundary_conditions", token) == 0, "ERROR: not a boundary condition line...: " << token);
  token = strtok(NULL, " \t");
  
  boundaryConditionType bct[6]={bDirichlet, bDirichlet, bDirichlet, bDirichlet, bStressFree, bDirichlet};
  
  int type, side;
  while (token != NULL)
  {
    if (startswith("#", token) || startswith(" ", buffer))
      // Ignore commented lines and lines with just a space.
      break;

    // low x
    if (startswith("lx=", token))
    {
      side = 0;
      token += 3;
      type = atoi(token);
    }
    else if (startswith("hx=", token))
    {
      side = 1;
      token += 3;
      type = atoi(token);
    }
    else if (startswith("ly=", token))
    {
      side = 2;
      token += 3;
      type = atoi(token);
    }
    else if (startswith("hy=", token))
    {
      side = 3;
      token += 3;
      type = atoi(token);
    }
    else if (startswith("lz=", token))
    {
      side = 4;
      token += 3;
      type = atoi(token);
    }
    else if (startswith("hz=", token))
    {
      side = 5;
      token += 3;
      type = atoi(token);
    }
    else
    {
      badOption("boundary_conditions", token);
    }
     
    switch (type) {
    case 0:
      bct[side] = bStressFree;
      break;
    case 1:
      bct[side] = bDirichlet;
      break;
    case 2:
      bct[side] = bSuperGrid;
      break;
    default:
      if (m_myRank==0)
      {
	printf("processBoundaryConditions:: Ignoring unknown boundary condition type = %i\n", type);
      }
    }
    token = strtok(NULL, " \t");
  } 
  set_global_bcs(bct);
}

void EW::processSupergrid(char *buffer)
{
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("supergrid", token) == 0, "ERROR: not a supergrid line...: " << token);
  token = strtok(NULL, " \t");
  int sg_thickness, sg_transition;
  double sg_coeff;
  bool thicknessSet=false, transitionSet=false, dampingCoeffSet=false;
  
  while (token != NULL)
  {
    if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
      break;

//                  1234567890
    if (startswith("thickness=", token)) // in number of grid sizes (different from WPP)
    {
      token += 10;
      sg_thickness = atoi(token);
      CHECK_INPUT(sg_thickness>0, "The number of grid points in the supergrid damping layer must be positive, not: "<< sg_thickness);
      thicknessSet = true;
    }
//                  12345678901
    if (startswith("transition=", token)) // in number of grid sizes (different from WPP)
    {
      token += 11;
      sg_transition = atoi(token);
      CHECK_INPUT(sg_transition>0, "The number of grid points in the supergrid transition layer must be positive, not: "<< sg_transition);
      transitionSet = true;
  }
//                       12345678901234567890
    else if (startswith("damping_coefficient=", token))
    {
      token += 20;
      sg_coeff = atof(token);
      CHECK_INPUT(sg_coeff>=0., "The supergrid damping coefficient must be non-negative, not: "<<sg_coeff);
      dampingCoeffSet=true;
    }
    else
    {
      badOption("supergrid", token);
    }
    token = strtok(NULL, " \t");
  } // end while token
  
  if (thicknessSet)
    set_sg_thickness(sg_thickness);

  if (transitionSet)
    set_sg_transition(sg_transition);

  if (dampingCoeffSet)
    set_sg_damping(sg_coeff);
}

// // void 
// // FileInput::
// // processRestart(char* buffer)
// // {
// //   int fromCycle=0, dumpInterval=0;
// //   string filePrefix = "wpp_restart";
  
// //   char* token = strtok(buffer, " \t");
// //   CHECK_INPUT(strcmp("restart", token) == 0, "ERROR: not a restart line...: " 
// // 	   << token);
// //   token = strtok(NULL, " \t");

// //   string err = "Restart Error: ";

// //   while (token != NULL)
// //     {
// //       // while there are tokens in the string still
// //      if (startswith("#", token) || startswith(" ", buffer))
// //         // Ignore commented lines and lines with just a space.
// //         break;
// //      if (startswith("fromCycle=", token) )
// //      {
// //         token += 10; // skip fromCycle=
// //         CHECK_INPUT(atoi(token) >= 0,
// //                 err << 
// //                 "fromCycle must be an integer >= to zero,  not: " << token);
// //         fromCycle = atoi(token);
// //      }
// //      else if (startswith("dumpInterval=", token) )
// //      {
// // 	     token += 13; // skip dumpInterval=
// //         CHECK_INPUT(atoi(token) >= 0,
// //                 err << "dumpInterval must be a positive int, not: " << token);
// //         dumpInterval = atoi(token);
// //      }
// //      else if (startswith("file=", token) )
// //      {
// //         token += 5; // skip file=
// //         filePrefix = token;
// //      }
// //      else
// //      {
// //         badOption("restart", token);
// //      }
// //      token = strtok(NULL, " \t");
// //     }
  
// //   mSimulation->setRestartInfo(fromCycle, dumpInterval, filePrefix);
// // }


//-----------------------------------------------------------------------
void EW::badOption(string name, char* option) const
{
   if (m_myRank == 0)
      cout << "\tWarning: ignoring " << name << " line option '" << option << "'" << endl;
}


// // //-----------------------------------------------------------------------
// // void FileInput::processCheckfornan( char* buffer )
// // {
// //    char* token = strtok(buffer, " \t");
// //    CHECK_INPUT(strcmp("checkfornan", token) == 0, "ERROR: not a check for nan line...: " << token);
// //    mSimulation->switch_on_checkfornan();
// // }


// //-----------------------------------------------------------------------

// void FileInput::processGlobalMaterial(char* buffer)
// {
//   char* token = strtok(buffer, " \t");
//   CHECK_INPUT(strcmp("globalmaterial", token) == 0, "ERROR: not an globalmaterial line...: " << token);
//   token = strtok(NULL, " \t");

//   string err = "globalmaterial error: ";
//   int modelnr = 0;
//   double frequency = 1;
//   bool useAttenuation = false;
//   double vpmin=0, vsmin=0;
  
//   while (token != NULL)
//   {
//      if ( startswith("vpmin=", token) )
//      {
// 	token += 6;
// 	vpmin = atof(token);
//      }
//      else if ( startswith("vsmin=", token) )
//      {
// 	token += 6;
// 	vsmin = atof(token);
//      }
//      else
// 	badOption("globalmaterial", token);
//      token = strtok(NULL, " \t");
//   }

//   mSimulation->set_threshold_velocities(vpmin, vsmin);
// }

// void FileInput::getEfileInfo(char* buffer)
// {
// #ifdef ENABLE_ETREE
//   // Used only for efiles
//   string accessmode = "parallel";

//   char* token = strtok(buffer, " \t");
//   CHECK_INPUT(strcmp("efile", token) == 0,
//            "ERROR: efile info can only be obtained from an efile line, not: " << token);

//   string commandName = token;

//   string err = token;
//   err += " Error: ";

//   token = strtok(NULL, " \t");

//   while (token != NULL)
//     {
//       // while there are tokens in the string still
//        if (startswith("#", token) || startswith(" ", buffer))
//           // Ignore commented lines and lines with just a space.
//           break;
// //       else if (startswith("model=", token))
// //       {
// //          token += 6; // skip model=
// //         model = token;
// //       }
//       else if (startswith("etree=", token))
//       {
//          token += 6; // skip etree=
//          m_topoFileName = token;
//       }
//       else if (startswith("xetree=", token))
//       {
//          token += 7; // skip xetree=
//          m_topoExtFileName = token;
//       }
//       else if (startswith("logfile=", token))
//       {
//          token += 8; // skip logfile=
//       }
//       else if (startswith("query=", token))
//       {
//          token += strlen("query=");
// 	 m_QueryType = token;
//          CHECK_INPUT(strcmp(token, "FIXEDRES") == 0 || 
//                  strcmp(token, "MAXRES") == 0,
//                  err << "query can only be set to FIXEDRES or MAXRES, not: " << m_QueryType);
//       }
//       else if (startswith("vsmin=", token))
//       {
//          token += strlen("vsmin=");
//       }
//       else if (startswith("vpmin=", token))
//       {
//          token += strlen("vpmin=");
//       }
// //      else if (startswith("water=", token))
// //      {
// //         token += strlen("water=");
// //          CHECK_INPUT(strcmp(token, "poisson") == 0 ||
// //                  strcmp(token, "noshear") == 0 ||
// //                  strcmp(token, "seafloor") == 0,
// //                  err << "water attribute can only be set to poisson,noshear or seafloor, not: " << token);
// //          if (strcmp(token, "noshear") == 0)
// //             waterMode = EtreeFile::NOSHEAR;
// //          else if (strcmp(token, "seafloor") == 0)
// //             waterMode = EtreeFile::SEAFLOOR;
// //          else if (strcmp(token, "poisson") == 0)
// //             waterMode = EtreeFile::POISSON;
// //      }
//       else if (startswith("access=", token))
//       {
//          token += strlen("access=");
//          CHECK_INPUT(strcmp(token, "parallel") == 0 ||
// 		     strcmp(token, "serial") == 0,
// 		     err << "access attribute can only be set to serial, or parallel, not: " << token);
//          accessmode = token;
//       }
//       else if( startswith("resolution=", token ) )
//       {
//         token += 11;
//         m_EFileResolution = atof(token);
//         CHECK_INPUT(m_EFileResolution>0.,"Resolution must be positive, not " << m_EFileResolution);
//       }
//       else
//       {
//          badOption(commandName, token);
//       }
//       token = strtok(NULL, " \t");
//     }
//   // End parsing...
  
// #else
//   CHECK_INPUT(0, "Error: Etree support not compiled into WPP (-DENABLE_ETREE)");
// #endif
// }

// //-----------------------------------------------------------------------
// void FileInput::processLimitfrequency(char* buffer)
// {
//    char* token = strtok(buffer, " \t");
//    CHECK_INPUT(strcmp("limitfrequency", token) == 0, "ERROR: not a limitfrequency line...: " << token);
//    token = strtok(NULL, " \t");

//    string err = "limitfrequency Error: ";
//    string commandName = token;
//    int ppw = 15;
//    while (token != NULL)
//    {
//       if (startswith("#", token) || startswith(" ", buffer))
//          break;

//       if (startswith("ppw=", token))
//       {
//         token += 4;
//         ppw = atoi(token);
//         CHECK_INPUT(ppw>0.,"points per wavelength must be positive, not " << ppw );
//       }
//       else
//       {
//          badOption(commandName, token);
//       }
//       token = strtok(NULL, " \t");
//    }
//    mSimulation->set_resolution( ppw );
// }

// //-----------------------------------------------------------------------
// void FileInput::processPrefilter(char* buffer)
// {
//    char* token = strtok(buffer, " \t");
//    CHECK_INPUT(strcmp("prefilter", token) == 0, "ERROR: not a prefilter line...: " << token);
//    token = strtok(NULL, " \t");

//    string err = "prefilter Error: ";
//    string commandName = token;
//    double corner_freq = 1.0, max_freq=1.0;
//    bool limit_source_freq=false, enable_prefilter=false;
//    while (token != NULL)
//    {
//       if (startswith("#", token) || startswith(" ", buffer))
//          break;

//       if (startswith("fc=", token))
//       {
//         token += 3;
//         corner_freq = atof(token);
//         CHECK_INPUT(corner_freq>0.,"corner frequency must be positive, not " << corner_freq );
// 	enable_prefilter=true;
//       }
// //                         1234567890
//       else if (startswith("maxfreq=", token))
//       {
//         token += 8;
//         max_freq = atof(token);
//         CHECK_INPUT(max_freq>0.,"max source freq parameter must be positive, not " << max_freq );
// 	limit_source_freq=true;	
//       }
//       else
//       {
//          badOption(commandName, token);
//       }
//       token = strtok(NULL, " \t");
//    }
//    mSimulation->set_prefilter( enable_prefilter, corner_freq, limit_source_freq, max_freq );
// }

// //-----------------------------------------------------------------------
// void FileInput::processGeodynbc(char* buf)
// {
//    // At this point, the geodyn file has already been read into m_geodyn_filename
//    char* token = strtok(buf, " \t");
//    CHECK_INPUT(strcmp("geodynbc", token) == 0, "ERROR: not a geodynbc line...: " << token);
//    CHECK_INPUT(m_geodynbc_found,"Error: geodynbc not obtained"<<token);

//    ifstream geodynfile(m_geodynbc_filename.c_str());
//    CHECK_INPUT(geodynfile.is_open(), "Error: opening geodyn file " << m_geodynbc_filename );


//    string err = "geodynbc Error: ";
//    string commandName = "geodynbc";

//    int faces=6, nx=0, ny=0, nz=0, nsteps=0, filter=0, adjust=1;
//    double x0, y0, z0, lat, lon, elev, az, timestep, rho=0, vs=0, vp=0, freq;
//    double srcx0, srcy0, srcz0, h, toff;

//    bool timestepset = false, nstepsset=false, toffset=false;
//    char buffer[256];
//    bool done = false;
//    while (!geodynfile.eof() && !done )
//    {
//       geodynfile.getline(buffer,256);
//       if (startswith("#", buffer) || startswith("\n", buffer) || buffer == "\0" )
//          break;
//       if( startswith("begindata",buffer) )
//       {
// 	 done = true;
//          break;
//       }

//       if( startswith("grid", buffer) )
//       {
// 	 char* token = strtok(buffer, " \t");
// 	 token = strtok(NULL, " \t");
// 	 while (token != NULL)
// 	 {
// 	    if (startswith("#", token) || startswith(" ", buffer))
// 	       break;
// 	    if (startswith("faces=", token))
// 	    {
// 	       token += 6;
// 	       faces = atoi(token);
// 	    }
// 	    else if( startswith("nx=",token))
// 	    {
// 	       token += 3;
// 	       nx = atoi(token);
// 	    }
// 	    else if( startswith("ny=",token))
// 	    {
// 	       token += 3;
// 	       ny = atoi(token);
// 	    }
// 	    else if( startswith("nz=",token))
// 	    {
// 	       token += 3;
// 	       nz = atoi(token);
// 	    }
// 	    else if( startswith("stepsize=",token))
// 	    {
// 	       token += 9;
// 	       h = atof(token);
// 	    }
// 	    else if( startswith("x0=",token))
// 	    {
// 	       token += 3;
// 	       x0 = atof(token);
// 	    }
// 	    else if( startswith("y0=",token))
// 	    {
// 	       token += 3;
// 	       y0 = atof(token);
// 	    }
// 	    else if( startswith("z0=",token))
// 	    {
// 	       token += 3;
// 	       z0 = atof(token);
// 	    }
// 	    else if( startswith("lat=",token))
// 	    {
// 	       token += 4;
// 	       lat = atof(token);
// 	    }
// 	    else if( startswith("lon=",token))
// 	    {
// 	       token += 4;
// 	       lon = atof(token);
// 	    }
// 	    else if( startswith("elev=",token))
// 	    {
// 	       token += 5;
// 	       elev = atof(token);
// 	    }
// 	    else if( startswith("az=",token))
// 	    {
// 	       token += 3;
// 	       az = atof(token);
// 	    }
// 	    else if( startswith("adjust=",token))
// 	    {
// 	       token += 7;
// 	       adjust = strcmp(token,"yes")==0;
// 	    }
// 	    else
// 	    {
// 	       badOption("geodyn-grid", token);
// 	    }
// 	    token = strtok(NULL, " \t");
// 	 }
//       }
//       else if( startswith("time", buffer) )
//       {
// 	 char* token = strtok(buffer, " \t");
// 	 token = strtok(NULL, " \t");
// 	 while (token != NULL)
// 	 {
// 	    if (startswith("#", token) || startswith(" ", buffer))
// 	       break;
// 	    if (startswith("timestep=", token))
// 	    {
// 	       token += 9;
// 	       timestep = atof(token);
// 	       timestepset=true;
// 	    }
// 	    else if( startswith("nsteps=",token))
// 	    {
// 	       token += 7;
// 	       nsteps = atoi(token);
// 	       nstepsset=true;
// 	    }
// 	    else if( startswith("toff=",token))
// 	    {
// 	       token += 5;
// 	       toff = atof(token);
// 	       toffset=true;
// 	    }
// 	    else
// 	    {
// 	       badOption("geodyn-time", token);
// 	    }
// 	    token = strtok(NULL, " \t");
// 	 }
//       }
//       else if( startswith("material",buffer) )
//       {
// 	 char* token = strtok(buffer, " \t");
// 	 token = strtok(NULL, " \t");
// 	 while (token != NULL)
// 	 {
// 	    if (startswith("#", token) || startswith(" ", buffer))
// 	       break;
// 	    if (startswith("rho=", token))
// 	    {
// 	       token += 4;
// 	       rho = atof(token);
// 	    }
// 	    else if( startswith("vs=",token))
// 	    {
// 	       token += 3;
// 	       vs = atof(token);
// 	    }
// 	    else if( startswith("vp=",token))
// 	    {
// 	       token += 3;
// 	       vp = atof(token);
// 	    }
// 	    else
// 	    {
// 	       badOption("geodyn-material", token);
// 	    }
// 	    token = strtok(NULL, " \t");
// 	 }
//       }
//       else if( startswith("source",buffer) )
//       {
// 	 char* token = strtok(buffer, " \t");
// 	 token = strtok(NULL, " \t");
// 	 while (token != NULL)
// 	 {
// 	    if (startswith("#", token) || startswith(" ", buffer))
// 	       break;
// 	    if (startswith("filter=", token))
// 	    {
// 	       token += 7;
// 	       if( strcmp(token,"butterworth")== 0 )
// 		  filter = 1;
// 	       else
// 		  filter = 0;
// 	    }
// 	    else if( startswith("frequency=",token))
// 	    {
// 	       token += 10;
// 	       freq = atof(token);
// 	    }
// 	    else if( startswith("x0=",token))
// 	    {
// 	       token += 3;
// 	       srcx0 = atof(token);
// 	    }
// 	    else if( startswith("y0=",token))
// 	    {
// 	       token += 3;
// 	       srcy0 = atof(token);
// 	    }
// 	    else if( startswith("z0=",token))
// 	    {
// 	       token += 3;
// 	       srcz0 = atof(token);
// 	    }
// 	    else
// 	    {
// 	       badOption("geodyn-source", token);
// 	    }
// 	    token = strtok(NULL, " \t");
// 	 }
//       }
//    }
//    geodynfile.close();
//    CHECK_INPUT( nx == ny, "Geodyn file error: x-y Cube dimensions must be equal, not "
//                   <<nx << " " << ny );
//    CHECK_INPUT( faces == 5 || faces == 6 , "Geodyn file error: Faces must be 5 or 6, not" << faces );
//    CHECK_INPUT( timestepset, "Geodyn file error: No time step given");
//    CHECK_INPUT( nstepsset, "Geodyn file error: Number of steps not given");
//    mSimulation->set_geodyn_data( m_geodynbc_filename, nx, nz, h, m_ibc_origin, timestep,
// 				 nsteps, faces );
// }

// //-----------------------------------------------------------------------
// void FileInput::geodynFindFile(char* buffer)
// {
//    char* token = strtok(buffer, " \t");
//    CHECK_INPUT(strcmp("geodynbc", token) == 0, "ERROR: not a geodynbc line...: " << token);
//    token = strtok(NULL, " \t");

//    string err = "geodynbc Error: ";
//    string commandName = token;

//    while (token != NULL)
//    {
//       if (startswith("#", token) || startswith(" ", buffer))
//          break;

//       if (startswith("file=", token))
//       {
//          token += 5;
//          m_geodynbc_filename = token;
//          m_geodynbc_found = true;
//       }
//       else if( startswith("center=",token))
//       {
//          token += 7;
//          if (atoi(token) == 1 || strcmp(token,"yes")==0 )
// 	    m_geodynbc_center = true;
//       }
//       else
//       {
//          badOption(commandName, token);
//       }
//       token = strtok(NULL, " \t");
//    }
// }

// //-----------------------------------------------------------------------
// void FileInput::geodynbcGetSizes( string filename, double origin[3], double &cubelen,
// 				  double& zcubelen, bool &found_latlon, double& lat, 
// 				  double& lon, double& az, int& adjust )
// {
//    ifstream geodynfile(m_geodynbc_filename.c_str());
//    CHECK_INPUT( geodynfile.is_open(), "Error: opening geodyn file " << m_geodynbc_filename );


//    string err = "geodynbc Error: ";
//    string commandName = "geodynbc";

//    int nx=0, ny=0, nz=0, faces=6;
//    double x0, y0, z0, elev, h;
//    adjust=1;

//    char buffer[256];
//    bool done = false;
//    bool nxfound=false, nyfound=false, nzfound=false, x0found=false, y0found=false, z0found=false;
//    bool latfound=false, lonfound=false, azfound=false, hfound=false, elevfound=false;
//    while (!geodynfile.eof() && !done )
//    {
//       geodynfile.getline(buffer,256);
//       if (startswith("#", buffer) || startswith("\n", buffer) || buffer == "\0" )
//          break;
//       if( startswith("begindata",buffer) )
//       {
// 	 done = true;
//          break;
//       }

//       if( startswith("grid", buffer) )
//       {
// 	 char* token = strtok(buffer, " \t");
// 	 token = strtok(NULL, " \t");
// 	 while (token != NULL)
// 	 {
// 	    if (startswith("#", token) || startswith(" ", buffer))
// 	       break;
// 	    if (startswith("faces=", token))
// 	    {
// 	       token += 6;
// 	       faces = atoi(token);
// 	    }
// 	    else if( startswith("nx=",token))
// 	    {
// 	       token += 3;
// 	       nx = atoi(token);
//                nxfound = true;
// 	    }
// 	    else if( startswith("ny=",token))
// 	    {
// 	       token += 3;
// 	       ny = atoi(token);
//                nyfound = true;
// 	    }
// 	    else if( startswith("nz=",token))
// 	    {
// 	       token += 3;
// 	       nz = atoi(token);
//                nzfound = true;
// 	    }
// 	    else if( startswith("stepsize=",token))
// 	    {
// 	       token += 9;
// 	       h = atof(token);
//                hfound = true;
// 	    }
// 	    else if( startswith("x0=",token))
// 	    {
// 	       token += 3;
// 	       x0 = atof(token);
//                x0found = true;
// 	    }
// 	    else if( startswith("y0=",token))
// 	    {
// 	       token += 3;
// 	       y0 = atof(token);
//                y0found = true;
// 	    }
// 	    else if( startswith("z0=",token))
// 	    {
// 	       token += 3;
// 	       z0 = atof(token);
//                z0found = true;
// 	    }
// 	    else if( startswith("lat=",token))
// 	    {
// 	       token += 4;
// 	       lat = atof(token);
//                latfound = true;
// 	    }
// 	    else if( startswith("lon=",token))
// 	    {
// 	       token += 4;
// 	       lon = atof(token);
//                lonfound = true;
// 	    }
// 	    else if( startswith("elev=",token))
// 	    {
// 	       token += 5;
// 	       elev = atof(token);
//                elevfound = true;
// 	    }
// 	    else if( startswith("az=",token))
// 	    {
// 	       token += 3;
// 	       az = atof(token);
//                azfound = true;
// 	    }
// 	    else if( startswith("adjust=",token))
// 	    {
// 	       token += 7;
// 	       //	       adjust = strcmp(token,"yes")==0;
//                adjust = atoi(token);
// 	    }
// 	    else
// 	    {
// 	       badOption("geodyn-grid", token);
// 	    }
// 	    token = strtok(NULL, " \t");
// 	 }
//       }
//    }
//    geodynfile.close();

//    if( nxfound && !nyfound )
//    {
//       ny = nx;
//       nyfound = true;
//    }
//    if( nxfound && !nzfound )
//    {
//       nz = nx;
//       nzfound = true;
//    }
//    CHECK_INPUT( (nxfound || nyfound || nzfound) && hfound, "Error in geodyn file: dimensions not specified");

//    if( nxfound )
//       cubelen = (nx-1)*h;
//    else if( nyfound )
//       cubelen = (ny-1)*h;
//    else
//       cubelen = (nz-1)*h;

//    found_latlon = latfound && lonfound && azfound;
//    if( elevfound )
//       origin[2] = -elev;
//    else
//       origin[2] = z0;

//    origin[0] = x0;
//    origin[1] = y0;

//    zcubelen = cubelen;
//    if( nzfound )
//       zcubelen = (nz-1)*h;
// }

// void FileInput::processMaterial( char* buffer )
// {
//   string name = "Material";

//   char* token = strtok(buffer, " \t");
//   CHECK_INPUT(strcmp("material", token) == 0,
// 	      "ERROR: material properties can be set by a material line, not: " << token);

//   string err = token;
//   err += " Error: ";

//   token = strtok(NULL, " \t");

//   int materialID=-1;
//   double vp0=-1, vs0=-1, rho0=-1, qp=-1, qs=-1;
//   double vp1=0, vs1=0, rho1=0;
//   double vp2=0, vs2=0, rho2=0;
//   double vp1o2=0, vs1o2=0, rho1o2=0;

//   bool gotID = false;
  
//   while (token != NULL)
//   {
//     // while there are tokens in the string still
//     if (startswith("#", token) || startswith(" ", buffer))
//       // Ignore commented lines and lines with just a space.
//       break;
// //                  1234567890
//     if (startswith("id=", token) )
//     {
//       token += 3; // skip id=
//       materialID = atoi(token);
//       gotID=true;
//     }
//     else if (startswith("vp=", token) )
//     {
//       token += 3; // skip vp=
//       vp0 = atof(token);
//     }
//     else if (startswith("vs=", token) )
//     {
//       token += 3; // skip vs=
//       vs0 = atof(token);
//     }
//     else if (startswith("rho=", token))
//     {
//       token += 4; // skip rho=
//       rho0 = atof(token);
//     }
// // linear variation
//     else if (startswith("rhograd=", token))
//     {
//       token += 8; // skip rhograd=
//       rho1 = atof(token);
//     }
//     else if (startswith("vpgrad=", token))
//     {
//       token += 7; // skip vpgrad=
//       vp1 = atof(token);
//     }
//     else if (startswith("vsgrad=", token))
//     {
//       token += 7; // skip vsgrad=
//       vs1 = atof(token);
//     }
// // quadratic variation
//     else if (startswith("rho2=", token))
//     {
//       token += 5; // skip rho2=
//       rho2 = atof(token);
//     }
//     else if (startswith("vp2=", token))
//     {
//       token += 4; // skip vp2=
//       vp2 = atof(token);
//     }
//     else if (startswith("vs2=", token))
//     {
//       token += 4; // skip vs2=
//       vs2 = atof(token);
//     }
// // sqrt variation
//     else if (startswith("rhosqrt=", token))
//     {
//       token += 8; // skip rhosqrt=
//       rho1o2 = atof(token);
//     }
//     else if (startswith("vpsqrt=", token))
//     {
//       token += 7; // skip vpsqrt=
//       vp1o2 = atof(token);
//     }
//     else if (startswith("vssqrt=", token))
//     {
//       token += 7; // skip vssqrt=
//       vs1o2 = atof(token);
//     }
// // attenuation variables
//     else if (startswith("Qp=", token) || startswith("qp=", token))
//     {
//       token += 3; // skip qp=
//       qp = atof(token);
//     }
//     else if (startswith("Qs=", token) || startswith("qs=", token))
//     {
//       token += 3; // skip qs=
//       qs = atof(token);
//     }
//     else
//     {
//       badOption("material", token);
//     }
//     token = strtok(NULL, " \t");
//   }
//   // End parsing...
  
//   CHECK_INPUT( gotID, "No id specified in material command");
  

//   CHECK_INPUT( (vs0 > 0 || vs1 != 0 || vs2 != 0) , 
// 	       "Error in material command: vs0, vs1, vs2 are " << vs0 << " " << vs1 << " " << vs2 );

//   CHECK_INPUT( (vp0 > 0 || vp1 != 0 || vp2 != 0) , 
// 	       "Error in material command: vp0, vp1, vp2 are " << vp0 << " " << vp1 << " " << vp2 );

//   CHECK_INPUT( (rho0 > 0 || rho1 != 0 || rho2 != 0) , 
// 	       "Error in material command: rho0, rho1, rho2 are " << rho0 << " " << rho1 << " " << rho2 );

//   if(mSimulation->getVerbosity() >=2 &&  m_myRank == 0 )
//   {
//     cout << "**** Material parameters: *****" << endl;
//     cout << "materialID=" << materialID << endl;
//     cout << "vp=" << vp0 <<  " vpgrad=" << vp1 << " vp2=" << vp2 << " vpsqrt=" << vp1o2 << endl;
//     cout << "vs=" << vs0 <<  " vsgrad=" << vs1 << " vs2=" << vs2 << " vssqrt=" << vs1o2 << endl;
//     cout << "rho=" << rho0 <<  " rhograd=" << rho1 << " rho2=" << rho2 << " rhosqrt=" << rho1o2 <<endl;
//     cout << "qp=" << qp <<  " qs=" << qs << endl;
//   }

// // add material to WPP2 object
//   MaterialProperty *mat=new MaterialProperty(materialID, vp0, vp1, vp2, vs0, vs1, vs2, rho0, rho1, rho2, qp, qs);
//   mat->setSqrtCoefficients( vp1o2, vs1o2, rho1o2 );
//   mSimulation->addMaterialProperty(mat);
// }

// //-----------------------------------------------------------------------
// void FileInput::processImage3D( char* buffer )
// {
//    int cycle=-1, cycleInterval=0, sample=1;
//    Image3D::Image3DMode mode=Image3D::RHO;
//    float time=0.0, timeInterval=0.0;
//    bool timingSet = false;
//    bool showsponge = false;
//    bool startTimeSet = false;
//    double tStart = 0;
//    string filePrefix="volimage";

//    bool boxset=false;
//    double x1=0, x2=m_xmax, y1=0, y2=m_ymax, z1=0, z2=m_zmax;
  

//   bool use_double = false;
  
//   char* token = strtok(buffer, " \t");
//   if ( strcmp("volimage", token) != 0 )
//   {
//     cerr << "Processing volimage command: " << "ERROR: not a volimage line...: " << token;
//     MPI_Abort( MPI_COMM_WORLD, 1 );
//   }
  
//   token = strtok(NULL, " \t");

//   string err = "volimage Error: ";

//   while (token != NULL)
//   {
//     // while there are tokens in the string still
//     if (startswith("#", token) || startswith(" ", buffer))
//       // Ignore commented lines and lines with just a space.
//       break;
//     if (startswith("time=", token) )
//     {
//       token += 5; // skip time=
//       if (atof(token) < 0.)	      
//       {
// 	cerr << "Processing volimage command: " << "time must be a non-negative number, not: " << token;
// 	MPI_Abort( MPI_COMM_WORLD, 1 );
//       }
//       time = atof(token);
//       timingSet = true;
//     }
// //                       1234567890123
//     else if (startswith("timeInterval=", token) )
//     {
//       token += 13; // skip timeInterval=
//       if (atof(token) <= 0.)	      
//       {
// 	cerr << "Processing volimage command: " << "timeInterval must be a positive number, not: " << token;
// 	MPI_Abort( MPI_COMM_WORLD, 1 );
//       }
//       timeInterval = atof(token);
//       timingSet = true;
//     }
// //                       1234567890
//     else if (startswith("startTime=", token) )
//     {
//       token += 10; // skip startTime=
//       tStart = atof(token);
//       startTimeSet = true;
//     }
//     else if (startswith("cycle=", token) )
//     {
//       token += 6; // skip cycle=
//       if (atoi(token) < 0)	      
//       {
// 	cerr << "Processing volimage command: " << "cycle must be a non-negative integer, not: " << token;
// 	MPI_Abort( MPI_COMM_WORLD, 1 );
//       }
//       cycle = atoi(token);
//       timingSet = true;
//     }
//     else if (startswith("cycleInterval=", token) )
//     {
//       token += 14; // skip cycleInterval=
//       if (atoi(token) <= 0)	      
//       {
// 	cerr << "Processing volimage command: " << "cycleInterval must be a positive integer, not: " << token;
// 	MPI_Abort( MPI_COMM_WORLD, 1 );
//       }
//       cycleInterval = atoi(token);
//       timingSet = true;
//     }
//     else if (startswith("file=", token))
//     {
//       token += 5; // skip file=
//       filePrefix = token;
//     }
//     else if (startswith("mode=", token))
//     {
//       token += 5; // skip mode=
//       if (strcmp(token, "ux") == 0)        mode = Image3D::UX;
//       else if (strcmp(token, "uy") == 0)   mode = Image3D::UY;
//       else if (strcmp(token, "uz") == 0)   mode = Image3D::UZ;
//       else if (strcmp(token, "rho") == 0)   mode = Image3D::RHO;
//       else if (strcmp(token, "p") == 0)   mode = Image3D::P;
//       else if (strcmp(token, "s") == 0)   mode = Image3D::S;
//       else if (strcmp(token, "div") == 0)   mode = Image3D::DIV;
//       else if (strcmp(token, "curl") == 0)   mode = Image3D::CURL;
//       else if (strcmp(token, "mag") == 0)   mode = Image3D::MAG;
//       else if (strcmp(token, "veldiv") == 0)   mode = Image3D::VELDIV;
//       else if (strcmp(token, "velcurl") == 0)   mode = Image3D::VELCURL;
//       else if (strcmp(token, "velmag") == 0)   mode = Image3D::VELMAG;
// //      else if (strcmp(token, "z") == 0)   mode = Image3D::ZCOORD; // AP: this mode can not be written explicitly
//       else
//       {
// 	mode = static_cast<Image3D::Image3DMode>(atoi(token));
// 	  cerr << "Processing image command: " << "mode must be one of the following: " << endl << endl
// 	       << "\t ux|uy|uz|rho|p|s|div|curl|mag|veldiv|velcurl|velmag" << endl
// 	       << "*not: " << token;
// 	  MPI_Abort( MPI_COMM_WORLD, 1 );
//       }
//     }
//     else if( startswith("precision=",token) )
//     {
//       token += 10;
// //      if ( !(strcmp(token,"double")==0 || strcmp(token,"float")==0) )
//       if ( !(startswith("double",token) || startswith("float",token)) )
//       {
// 	cerr << "Processing volimage command: " << " precision must be float or double, not '" << token << "'" << endl;
// 	MPI_Abort( MPI_COMM_WORLD, 1 );
//       }
//       use_double =  startswith("double",token)==0;
//     }
//     else if( startswith("sample=",token) )
//     {
//       token += 7;
//       sample = atoi(token);
//       if( sample < 1 )
//       {
// 	cerr << "Processing volimage command: Sample must be greater than one, not " << sample << endl;
// 	MPI_Abort( MPI_COMM_WORLD, 1 );
//       }
//     }
//     else if( startswith("savelayer=",token) )
//     {
//        token += 10;
//        if( strcmp(token,"yes")==0 || strcmp(token,"on")==0 )
// 	  showsponge = true;
//     }
//     else if (startswith("x1=", token))
//     {
//        token += 3; // skip x1=
//        x1 = atof(token);
//        boxset = true;
//     }
//     else if (startswith("x2=", token))
//     {
//        token += 3; // skip x2=
//        x2 = atof(token);
//        boxset = true;
//     }
//     else if (startswith("y1=", token))
//     {
//        token += 3; // skip y1=
//        y1 = atof(token);
//        boxset = true;
//     }
//     else if (startswith("y2=", token))
//     {
//        token += 3; // skip y2=
//        y2 = atof(token);
//        boxset = true;
//     }
//     else if (startswith("z1=", token))
//     {
//        token += 3; // skip z1=
//        z1 = atof(token);
//        boxset = true;
//     }
//     else if (startswith("z2=", token))
//     {
//        token += 3; // skip z2=
//        z2 = atof(token);
//        boxset = true;
//     }
    
//     else
//     {
//       badOption("volimage", token);
//     }
//     token = strtok(NULL, " \t");
//   }
  
//   if ( !timingSet )
//   {
//     cerr << "Processing volimage command: " << "at least one timing mechanism must be set: cycle, time, cycleInterval or timeInterval" << endl;
//     MPI_Abort( MPI_COMM_WORLD, 1 );
//   }
//   CHECK_INPUT(x1 >= 0.,
// 	      err << "x1 is less than the minimum x, " 
// 	      << x1 << " < " << 0.);
//   CHECK_INPUT(x1 <= m_xmax,
// 	      err << "x1 is greater than the maximum x, " 
// 	      << x1 << " > " << m_xmax);
//   CHECK_INPUT(x2 >= 0.,
// 	      err << "x2 is less than the minimum x, " 
// 	      << x2 << " < " << 0.);
//   CHECK_INPUT(x2 <= m_xmax,
// 	      err << "x2 is greater than the maximum x, " 
// 	      << x2 << " > " << m_xmax);
//   CHECK_INPUT(x2 > x1,
// 	      err << "x2 is smaller than x1 "
// 	      << x2 << " < " << x1 );
//   CHECK_INPUT(y1 >= 0.,
// 	     err << "y1 is less than the minimum y, " << y1 << " < " << 0.);
//   CHECK_INPUT(y1 <= m_ymax,
// 	      err << "y1 is greater than the maximum y, " << y1 << " > " << m_ymax);
//   CHECK_INPUT(y2 >= 0.,
// 	      err << "y2 is less than the minimum y, " << y2 << " < " << 0.);
//   CHECK_INPUT(y2 <= m_ymax,
// 	      err << "y2 is greater than the maximum y, " << y2 << " > " << m_ymax);
//   CHECK_INPUT(y2 > y1,
// 	      err << "y2 is smaller than y1 "
// 	      << y2 << " < " << y1 );
//   CHECK_INPUT(mSimulation->topographyExists() || z1 >= 0.,
// 	      err << "z1 is less than the minimum z, " << z1 << " < " << 0.);
//   CHECK_INPUT(z1 <= m_zmax, 
// 	      err << "z1 is greater than the maximum z, " << z1 << " > " << m_zmax);
//   CHECK_INPUT(z2 >= 0.,
// 	      err << "z2 is less than the minimum z, " << z2 << " < " << 0.);
//   CHECK_INPUT(z2 <= m_zmax,
// 	      err << "z2 is greater than the maximum z, " << z2 << " > " << m_zmax);

// // add gridfile name to constructor (only with topography)
//   Image3D* im3 = new Image3D( mSimulation, time, timeInterval, cycle, cycleInterval, 
// 			      filePrefix, sample, mode, use_double, mSimulation->topographyExists() );
//   if( boxset )
//     im3->set_boundingbox( x1, x2, y1, y2, z1, z2 );
//   if( showsponge )
//     im3->show_spongelayer();
//   if ( startTimeSet )
//     im3->set_start_time( tStart );
  
//   mSimulation->addImage3D( im3 );

//   if (mSimulation->topographyExists())
//   {
// // make a 3D grid file holding the z-coordinates. This file will only be saved once (before the timestepping starts)
//     cycle=0;
//     cycleInterval=0;
//     time=0.0;
//     timeInterval=0.0;
//     mode=Image3D::ZCOORD;
  
//     Image3D* grid_im3 = new Image3D( mSimulation, time, timeInterval, cycle, cycleInterval, 
// 				     filePrefix, sample, mode, use_double, false );
//     if( boxset )
//       grid_im3->set_boundingbox( x1, x2, y1, y2, z1, z2 );
//     if( showsponge )
//       grid_im3->show_spongelayer();
//     mSimulation->addImage3D( grid_im3 );
//   }

// }
  
//-----------------------------------------------------------------------
void EW::setOutputPath(const string& path) 
{ 
  stringstream s;
  s << path << "/";
  mPath = s.str();
}

//-----------------------------------------------------------------------
void EW::setIO_timing(bool iotiming) 
{ 
  m_iotiming = iotiming;
}

//-----------------------------------------------------------------------
void EW::setParallel_IO(bool pfs, int nwriters) 
{ 
  m_pfs = pfs;
  m_nwriters = nwriters;
}

//-----------------------------------------------------------------------
void EW::setGoalTime(double t) 
{ 
  mTmax = t; 
  mTstart = 0.0; 
  mTimeIsSet = true;
}

//-----------------------------------------------------------------------
void EW::setNumberSteps(int steps)
{
  mNumberOfTimeSteps = steps;
  mTimeIsSet = false;
}


//-----------------------------------------------------------------------
int EW::getNumberOfSteps() const
{
  return mNumberOfTimeSteps;
}

//-----------------------------------------------------------------------
void EW::switch_on_error_log()
{
   m_error_log = true;
}

//-----------------------------------------------------------------------
void EW::set_energylog( string logfile, bool print, bool elog )
{
   m_energy_log = elog;
   m_energy_logfile = logfile;
   m_energy_print = print;
}


//-----------------------------------------------------------------------
void EW::set_cflnumber( double cfl )
{
   mCFL = cfl;
}

//-----------------------------------------------------------------------
void EW::set_twilight_forcing( ForcingTwilight* a_forcing )
{
   m_twilight_forcing = a_forcing;
   set_testing_mode(true);
}

//-----------------------------------------------------------------------
void EW::setAttenuationParams(int numberOfMechanisms, double velocityOmega, 
			      int ppw, double maxfrequency )
{
  m_number_mechanisms = numberOfMechanisms;
  m_velo_omega = velocityOmega;
  m_att_use_max_frequency = (ppw <= 0);
  m_att_ppw = ppw;
  m_att_max_frequency = maxfrequency;
}


//-----------------------------------------------------------------------
void EW::allocateCartesianSolverArrays(double a_global_zmax)
{
//
// note that this routine might modify the value of m_global_zmax
//
   if (mVerbose && proc_zero())
     printf("allocateCartesianSolverArrays: #ghost points=%i, #parallel padding points=%i\n", m_ghost_points, m_ppadding);

   int nCartGrids = m_refinementBoundaries.size();
   int refFact = 1;
   for( int r = 0 ; r < nCartGrids-1 ; r++ )
   {
      refFact *= 2;
      //      cout << "refinement boundary " << r << " is " << m_refinementBoundaries[r] << endl;
   }

//    m_topography_exists tells if there was a topography command in the input file?   

// is there an attenuation command in the file?
   if (!m_use_attenuation)
      m_number_mechanisms = 0;

   int nx_finest_w_ghost = refFact*(m_nx_base-1)+1+2*m_ghost_points;
   int ny_finest_w_ghost = refFact*(m_ny_base-1)+1+2*m_ghost_points;
   int proc_max[2];
// this info is obtained by the contructor
//   MPI_Comm_size( MPI_COMM_WORLD, &nprocs  );
   proc_decompose_2d( nx_finest_w_ghost, ny_finest_w_ghost, m_nProcs, proc_max );
   int is_periodic[2]={0,0};
   MPI_Cart_create( MPI_COMM_WORLD, 2, proc_max, is_periodic, true, &m_cartesian_communicator );
   int my_proc_coords[2];
   MPI_Cart_get( m_cartesian_communicator, 2, proc_max, is_periodic, my_proc_coords );
   MPI_Cart_shift( m_cartesian_communicator, 0, 1, m_neighbor, m_neighbor+1 );
   MPI_Cart_shift( m_cartesian_communicator, 1, 1, m_neighbor+2, m_neighbor+3 );
   //   if( proc_zero() )
   //   {
   //      cout << " Grid distributed on " << m_nProcs << " processors " << endl;
   //      cout << " Finest grid size    " << nx_finest_w_ghost << " x " << ny_finest_w_ghost << endl;
   //      cout << " Processor array     " << proc_max[0] << " x " << proc_max[1] << endl;
   //   }
   int ifirst, ilast, jfirst, jlast;
   decomp1d( nx_finest_w_ghost, my_proc_coords[0], proc_max[0], ifirst, ilast );
   decomp1d( ny_finest_w_ghost, my_proc_coords[1], proc_max[1], jfirst, jlast );

   ifirst -= m_ghost_points;
   ilast  -= m_ghost_points;
   jfirst -= m_ghost_points;
   jlast  -= m_ghost_points;
   int nx = nx_finest_w_ghost-2*m_ghost_points;
   int ny = ny_finest_w_ghost-2*m_ghost_points;
   int kfirst, klast;

   //   cout << "nCartGrids = " << nCartGrids << endl;
   mNumberOfCartesianGrids = nCartGrids;
   mNumberOfGrids = mNumberOfCartesianGrids;
   if( m_topography_exists )
      mNumberOfGrids++;

   // mU.resize(mNumberOfGrids);
   // mUm.resize(mNumberOfGrids);
   // mUp.resize(mNumberOfGrids);
   // mUacc.resize(mNumberOfGrids);
//   mLu.resize(mNumberOfGrids);
//   mF.resize(mNumberOfGrids);
   mMu.resize(mNumberOfGrids);
   mLambda.resize(mNumberOfGrids);
   mRho.resize(mNumberOfGrids);
   m_zmin.resize(mNumberOfGrids);
   mGridSize.resize(mNumberOfGrids);
   mMinVsOverH.resize(mNumberOfGrids);

// always allocate the pointer arrays for the viscoelastic properties (allows. e.g., mQs[g].is_defined() to be called)
   mQs.resize(mNumberOfGrids);
   mQp.resize(mNumberOfGrids);

// viscoelastic material coefficients and memory variables are only allocated when attenuation is enabled
// Allocate pointers, even if attenuation not used, for avoid segfault in parameter list with mMuVE[g], etc...
   mMuVE.resize(mNumberOfGrids);
   mLambdaVE.resize(mNumberOfGrids);
// // Allocate pointers, even if attenuation not used, for avoid segfault in parameter list with mMuVE[g], etc...
//    mAlphaVE.resize(mNumberOfGrids);
//    mAlphaVEm.resize(mNumberOfGrids);
//    mAlphaVEp.resize(mNumberOfGrids);
   if (m_use_attenuation)
   {
     mOmegaVE.resize(m_number_mechanisms); // global relaxation frequencies (1 per mechanism)
     
// muVE and lambdaVE are vectors of vectors
     for (int g=0; g<mNumberOfGrids; g++)
     {
       mMuVE[g]     = new Sarray[m_number_mechanisms];
       mLambdaVE[g] = new Sarray[m_number_mechanisms];
       // mAlphaVE[g]  = new Sarray[m_number_mechanisms];
       // mAlphaVEp[g] = new Sarray[m_number_mechanisms];
       // mAlphaVEm[g] = new Sarray[m_number_mechanisms];
     }
   }
   
   m_iStart.resize(mNumberOfGrids);
   m_iEnd.resize(mNumberOfGrids);
   m_jStart.resize(mNumberOfGrids);
   m_jEnd.resize(mNumberOfGrids);
   m_kStart.resize(mNumberOfGrids);
   m_kEnd.resize(mNumberOfGrids);

   m_iStartInt.resize(mNumberOfGrids);
   m_iEndInt.resize(mNumberOfGrids);
   m_jStartInt.resize(mNumberOfGrids);
   m_jEndInt.resize(mNumberOfGrids);
   m_kStartInt.resize(mNumberOfGrids);
   m_kEndInt.resize(mNumberOfGrids);

   m_global_nx.resize(mNumberOfGrids);
   m_global_ny.resize(mNumberOfGrids);
   m_global_nz.resize(mNumberOfGrids);

   m_onesided.resize(mNumberOfGrids);
   for( int g= 0 ;g < mNumberOfGrids ; g++ )
      m_onesided[g] = new int[6];

   m_bcType.resize(mNumberOfGrids);
   for( int g= 0 ;g < mNumberOfGrids ; g++ )
      m_bcType[g] = new boundaryConditionType[6];

   // m_BCForcing.resize(mNumberOfGrids);
   m_NumberOfBCPoints.resize(mNumberOfGrids);
   m_BndryWindow.resize(mNumberOfGrids);

   int *wind;
   
   for( int g= 0 ;g < mNumberOfGrids ; g++ )
   {
     m_NumberOfBCPoints[g] = new int[6];
     // m_BCForcing[g] = new double *[6];
     m_BndryWindow[g] = new int [36]; // 6 by 6 array in Fortran
     for (int side=0; side < 6; side++)
     {
       m_NumberOfBCPoints[g][side]=0;
       // m_BCForcing[g][side]=NULL;
       for (int qq=0; qq<6; qq+=2) // 0, 2, 4
	 m_BndryWindow[g][qq + side*6]= 999;
       for (int qq=1; qq<6; qq+=2) // 1, 3, 5
	 m_BndryWindow[g][qq + side*6]= -999;
     }
   }
   
   double h = m_h_base;

// save the grid spacing for all Cartesian grids
   for (int g=0; g<nCartGrids; g++)
   {
     mGridSize[g] = h;
     h = h/2.;
   }
   
   m_global_nx[nCartGrids-1] = nx_finest_w_ghost-2*m_ghost_points;
   m_global_ny[nCartGrids-1] = ny_finest_w_ghost-2*m_ghost_points;
   
   for (int g=nCartGrids-2; g >=0; g--)
   {
     m_global_nx[g] = 1 + (m_global_nx[g+1]-1)/2;
     m_global_ny[g] = 1 + (m_global_ny[g+1]-1)/2;
   }

// the curvilinear grid has a variable grid size, but matches the finest Cartesian grid where they meet
   if( m_topography_exists )
   {
     mGridSize[mNumberOfGrids-1]   = mGridSize[mNumberOfGrids-2];
     m_global_nx[mNumberOfGrids-1] = m_global_nx[mNumberOfGrids-2];
     m_global_ny[mNumberOfGrids-1] = m_global_ny[mNumberOfGrids-2];
   }
   
// Define grid in z-direction, by formula z_k = (k-1)*h + zmin
   vector<int> nz;
   nz.resize(nCartGrids);

// don't change the zmin of the finest cartesian grid
   m_zmin[nCartGrids-1] = m_refinementBoundaries[nCartGrids-1];
   for( int g = nCartGrids-1; g >= 0; g-- )
   {
     double zmax = (g>0? m_refinementBoundaries[g-1]:a_global_zmax);
     nz[g]     = 1 + static_cast<int> ( (zmax - m_zmin[g])/mGridSize[g] + 0.5 );

     zmax = m_zmin[g] + (nz[g]-1)*mGridSize[g];
      
     m_global_nz[g] = nz[g]; // save the number of grid points in the z-direction

     if( g>0 )
       m_zmin[g-1] = zmax;
     else
       a_global_zmax = zmax;
   }

// extent of computational grid (without ghost points)
   m_global_xmax = mGridSize[nCartGrids-1]*(m_global_nx[nCartGrids-1] - 1);
   m_global_ymax = mGridSize[nCartGrids-1]*(m_global_ny[nCartGrids-1] - 1);
   m_global_zmax = m_zmin[0] + (nz[0]-1)*mGridSize[0];
   if (mVerbose >= 1 && proc_zero())
     cout << "Extent of the computational domain xmax=" << m_global_xmax << " ymax=" << m_global_ymax << " zmax=" << 
       m_global_zmax << endl;
   

// tmp
   if( mVerbose >= 1 && proc_zero() )
   {
     cout << "Corrected global_zmax = " << m_global_zmax << endl;
     cout << "Refinement levels after correction: " << endl;
     for( int g=0; g<nCartGrids; g++ )
     {
       cout << "grid=" << g << " min Z=" << m_zmin[g] << endl;
     }
   }
   
// Define grid arrays, loop from finest to coarsest
   for( int g = nCartGrids-1 ; g >= 0 ; g-- )
   {
// NOTE: same number of ghost points in all directions
      kfirst     = 1-m_ghost_points;
      klast      = nz[g] + m_ghost_points;
      //      cout << "defining " << ifirst << " " << ilast << " " << jfirst << " " << jlast << " " << kfirst << " " << klast << endl;
      //      cout << " zmin " << m_zmin[g] << endl;
      // mU[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      // mUm[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      // mUp[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      // mUacc[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
      // mLu[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);

      // mF[g].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);

// elastic material
      mMu[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
      mRho[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
      mLambda[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);

// initialize the material coefficients to -1
      mMu[g].set_to_minusOne();
      mRho[g].set_to_minusOne();
      mLambda[g].set_to_minusOne();

// viscoelastic material coefficients & memory variables
      if (m_use_attenuation)
      {
	mQs[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	mQp[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	for (int a=0; a<m_number_mechanisms; a++)
	{
	  mMuVE[g][a].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	  mLambdaVE[g][a].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	  // mAlphaVE[g][a].define( 3,ifirst,ilast,jfirst,jlast,kfirst,klast);
	  // mAlphaVEp[g][a].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
	  // mAlphaVEm[g][a].define(3,ifirst,ilast,jfirst,jlast,kfirst,klast);
	  
// initialize the viscoelastic material coefficients to -1
	  mMuVE[g][a].set_to_minusOne();
	  mLambdaVE[g][a].set_to_minusOne();
	}
// initialize Qp and Qs to -1
	mQs[g].set_to_minusOne();
	mQp[g].set_to_minusOne();
      }
	

// save the local index bounds
      m_iStart[g] = ifirst;
      m_iEnd[g]   = ilast;
      m_jStart[g] = jfirst;
      m_jEnd[g]   = jlast;
      m_kStart[g] = kfirst;
      m_kEnd[g]   = klast;

// local index bounds for interior points (= no ghost or parallel padding points)
      if (ifirst == 1-m_ghost_points)
	m_iStartInt[g] = 1;
      else
	m_iStartInt[g] = ifirst+m_ppadding;

      if (ilast == nx + m_ghost_points)
	m_iEndInt[g]   = nx;
      else
	m_iEndInt[g]   = ilast - m_ppadding;

      if (jfirst == 1-m_ghost_points)
	m_jStartInt[g] = 1;
      else
	m_jStartInt[g] = jfirst+m_ppadding;

      if (jlast == nx + m_ghost_points)
	m_jEndInt[g]   = nx;
      else
	m_jEndInt[g]   = jlast - m_ppadding;

      m_kStartInt[g] = 1;
      m_kEndInt[g]   = nz[g];

      // go to next coarser grid 
      coarsen1d( nx, ifirst, ilast );
      coarsen1d( ny, jfirst, jlast );
      //      cout << g << " " << my_proc_coords[0] << " I Split into " << ifirst << " , " << ilast << endl;
      //      cout << g << " " << my_proc_coords[1] << " J Split into " << jfirst << " , " << jlast << endl;
      //      cout << "grid " << g << " zmin = " << m_zmin[g] << " nz = " << nz[g] << " kinterval " << kfirst << " , " << klast << endl;
      
   }

   if( m_topography_exists )
   {
// get the size from the top Cartesian grid
      ifirst = m_iStart[nCartGrids-1];
      ilast  = m_iEnd[nCartGrids-1];
      jfirst = m_jStart[nCartGrids-1];
      jlast  = m_jEnd[nCartGrids-1];

// save the local index bounds
      m_iStart[mNumberOfGrids-1] = ifirst;
      m_iEnd[mNumberOfGrids-1]   = ilast;
      m_jStart[mNumberOfGrids-1] = jfirst;
      m_jEnd[mNumberOfGrids-1]   = jlast;

// interior points (what about k-direction???)
      m_iStartInt[mNumberOfGrids-1] = m_iStartInt[nCartGrids-1];
      m_iEndInt[mNumberOfGrids-1]   = m_iEndInt[nCartGrids-1];
      m_jStartInt[mNumberOfGrids-1] = m_jStartInt[nCartGrids-1];
      m_jEndInt[mNumberOfGrids-1]   = m_jEndInt[nCartGrids-1];

// 3 versions of the topography:
      mTopo.define(ifirst,ilast,jfirst,jlast,1,1); // true topography/bathymetry, read directly from etree
      mTopoGrid.define(ifirst,ilast,jfirst,jlast,1,1); // smoothed version of true topography
// At this point, we don't know the number of grid points in the k-direction of the curvi-linear grid.
// the arrays mX, mY, mZ must be allocated by the grid generator
   }
// the topoMat array is needed for "squishing" the material properties near the free surface when a
// Cartesian grid is used on top
   mTopoMat.define(m_iStart[mNumberOfGrids-1], m_iEnd[mNumberOfGrids-1], 
		   m_jStart[mNumberOfGrids-1], m_jEnd[mNumberOfGrids-1], 1, 1); 
// mTopoMat holds the highest elevation where the etree returns solid material properties

// tmp
//   int myRank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//   for (int q=0; q<mNumberOfCartesianGrids; q++)
//   {
//     printf("Proc #%i, m_iEnd[%i]=%i, m_global_nx[%i]=%i, m_jEnd[%i]=%i, m_global_ny[%i]=%i\n", 
// 	   myRank, q, m_iEnd[q], q, m_global_nx[q], q, m_jEnd[q], q, m_global_ny[q]);
   
//   }

}

void EW::deprecatedImageMode(int value, const char* name) const
{
  if (m_myRank == 0)
    cout << "***Warning image mode using integers is deprecated, mode="
	 << value << " should be mode=" << name << " instead." << endl;
}

void EW::processImage(char* buffer)
{
   int cycle=-1, cycleInterval=0, sample=0;
//   int pfs = 0, nwriters=1;
   Image::ImageMode mode=Image::RHO;
   float time=0.0, timeInterval=0.0;
   bool timingSet = false;
   string filePrefix="image";

   // -----------------------------------------------------
   // It is only valid to set one of the following:
   //   x, or y, or z
   // 
   // The selection of one coordinate specifies plane which will be output to disk.
   // It is an error to select more than one.
   // -----------------------------------------------------
  Image::ImageOrientation locationType=Image::UNDEFINED;
  float coordValue;
  int gridPointValue;
  bool coordWasSet = false;
  bool use_double = false;
  
  char* token = strtok(buffer, " \t");
  if ( strcmp("image", token) != 0 )
  {
    cerr << "Processing image command: " << "ERROR: not an image line...: " << token;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  
  token = strtok(NULL, " \t");

  string err = "Image Error: ";

  while (token != NULL)
  {
    // while there are tokens in the string still
    if (startswith("#", token) || startswith(" ", buffer))
      // Ignore commented lines and lines with just a space.
      break;
    if (startswith("time=", token) )
    {
      token += 5; // skip time=
      if (atof(token) < 0.)	      
      {
	cerr << "Processing image command: " << "time must be a non-negative number, not: " << token;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      time = atof(token);
      timingSet = true;
    }
    else if (startswith("timeInterval=", token) )
    {
      token += 13; // skip timeInterval=
      if (atof(token) <= 0.)	      
      {
	cerr << "Processing image command: " << "timeInterval must be a positive number, not: " << token;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      timeInterval = atof(token);
      timingSet = true;
    }
    else if (startswith("cycle=", token) )
    {
      token += 6; // skip cycle=
      if (atoi(token) < 0)	      
      {
	cerr << "Processing image command: " << "cycle must be a non-negative integer, not: " << token;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      cycle = atoi(token);
      timingSet = true;
    }
    else if (startswith("cycleInterval=", token) )
    {
      token += 14; // skip cycleInterval=
      if (atoi(token) <= 0)	      
      {
	cerr << "Processing image command: " << "cycleInterval must be a positive integer, not: " << token;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      cycleInterval = atoi(token);
      timingSet = true;
    }
    else if (startswith("file=", token))
    {
      token += 5; // skip file=
      filePrefix = token;
    }
    else if (startswith("mode=", token))
    {
      token += 5; // skip mode=
      if (strcmp(token, "ux") == 0)        mode = Image::UX;
      else if (strcmp(token, "uy") == 0)   mode = Image::UY;
      else if (strcmp(token, "uz") == 0)   mode = Image::UZ;
      else if (strcmp(token, "rho") == 0)   mode = Image::RHO;
      else if (strcmp(token, "lambda") == 0)   mode = Image::LAMBDA;
      else if (strcmp(token, "mu") == 0)   mode = Image::MU;
      // else if (strcmp(token, "p") == 0)   mode = Image::P;
      // else if (strcmp(token, "s") == 0)   mode = Image::S;
      // else if (strcmp(token, "div") == 0)   mode = Image::DIV;
      // else if (strcmp(token, "curl") == 0)   mode = Image::CURL;
      // else if (strcmp(token, "veldiv") == 0)   mode = Image::VELDIV;
      // else if (strcmp(token, "velcurl") == 0)   mode = Image::VELCURL;
      // else if (strcmp(token, "lat") == 0)   mode = Image::LAT;
      // else if (strcmp(token, "lon") == 0)   mode = Image::LON;
      // else if (strcmp(token, "hvelmax") == 0) mode = Image::HVELMAX;
      // else if (strcmp(token, "vvelmax") == 0) mode = Image::VVELMAX;
      // else if (strcmp(token, "topo") == 0) mode = Image::TOPO;
      // else if (strcmp(token, "grid") == 0) mode = Image::GRID;
      // else if (strcmp(token, "uxerr") == 0)   mode = Image::UXERR;
      // else if (strcmp(token, "uyerr") == 0)   mode = Image::UYERR;
      // else if (strcmp(token, "uzerr") == 0)   mode = Image::UZERR;
      // else if (strcmp(token, "fx") == 0)   mode = Image::FX;
      // else if (strcmp(token, "fy") == 0)   mode = Image::FY;
      // else if (strcmp(token, "fz") == 0)   mode = Image::FZ;
      // else if (strcmp(token, "velmag") == 0)   mode = Image::VELMAG;
      // else if (strcmp(token, "qs") == 0) mode = Image::QS;
      // else if (strcmp(token, "qp") == 0) mode = Image::QP;
      // else if (strcmp(token, "hvel") == 0) mode = Image::HVEL;
      else
      {
// deprecated, but still possible to alternatively give a numeric value for mode
	mode = static_cast<Image::ImageMode>(atoi(token));
	if (mode <= 0 || mode > 28)
	{
	  cerr << "Processing image command: " << "mode must be one of the following: " << endl
	       << "ux|uy|uz|rho|lambda|mu" << endl 
//"|p|s|div|curl|veldiv|velcurl " << endl
//	       << "lat|lon|hvelmax|vvelmax|topo|grid|uxerr|uyerr|uzerr " << endl
//	       << "fx|fy|fz|velmag|qs|qp|hvel " << endl << endl
	       << "*not: " << token << endl;
	  
	  MPI_Abort( MPI_COMM_WORLD, 1 );
	}
// what is the purpose of the following if statement?           
	if (mode == Image::UX) deprecatedImageMode(mode, "ux"); 
	else if (mode == Image::UY) deprecatedImageMode(mode, "uy"); 
	else if (mode == Image::UZ) deprecatedImageMode(mode, "uz");
	else if (mode == Image::RHO) deprecatedImageMode(mode, "rho");
	else if (mode == Image::LAMBDA) deprecatedImageMode(mode, "lambda");
	else if (mode == Image::MU) deprecatedImageMode(mode, "mu");
	// else if (mode == Image::P) deprecatedImageMode(mode, "p");
	// else if (mode == Image::S) deprecatedImageMode(mode, "s");
	// else if (mode == Image::DIV) deprecatedImageMode(mode, "div");
	// else if (mode == Image::CURL) deprecatedImageMode(mode, "curl");
	// else if (mode == Image::VELDIV) deprecatedImageMode(mode, "veldiv");
	// else if (mode == Image::VELCURL) deprecatedImageMode(mode, "velcurl");
	// else if (mode == Image::LAT)deprecatedImageMode(mode, "lat");
	// else if (mode == Image::LON)deprecatedImageMode(mode, "lon");
	// else if (mode == Image::HVELMAX)deprecatedImageMode(mode, "hvelmax");
	// else if (mode == Image::VVELMAX)deprecatedImageMode(mode, "vvelmax");
	// else if (mode == Image::TOPO)deprecatedImageMode(mode, "topo");
	// else if (mode == Image::GRID)deprecatedImageMode(mode, "grid");
	// else if (mode == Image::UXERR)deprecatedImageMode(mode, "uxerr");
	// else if (mode == Image::UYERR)deprecatedImageMode(mode, "uyerr");
	// else if (mode == Image::UZERR)deprecatedImageMode(mode, "uzerr");
	// else if (mode == Image::FX)deprecatedImageMode(mode, "fx");
	// else if (mode == Image::FY)deprecatedImageMode(mode, "fy");
	// else if (mode == Image::FZ)deprecatedImageMode(mode, "fz");
	// else if (mode == Image::VELMAG)deprecatedImageMode(mode, "velmag");
	// else if (mode == Image::QS)deprecatedImageMode(mode, "qs");
	// else if (mode == Image::QP)deprecatedImageMode(mode, "qp");
      }
    }
    else if( startswith("precision=",token) )
    {
      token += 10;
      if ( !(strcmp(token,"double")==0 || strcmp(token,"float")==0) )
      {
	cerr << "Processing image command: " << " precision must be float or double, not " << token << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      use_double =  strcmp(token,"double")==0;
    }
    else if (startswith("x=", token))
    {
      token += 2; // skip x=
      if ( coordWasSet )
      {
	cerr << "Processing image command: " << "cannot set a coordinate location twice, x and " << locationType << " were both set." << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      coordWasSet=true;
      locationType=Image::X;
      coordValue = atof(token);
      if ( coordValue < 0.0 || coordValue > m_global_xmax )
      {
	cerr << "Processing image command: " << "x value must be within the computational domain, not: " << coordValue << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
    }
     
    else if (startswith("y=", token))
    {
      token += 2; // skip y=
      if ( coordWasSet )
      {
	cerr << "Processing image command: " << "cannot set a coordinate location twice, y and " << locationType << " were both set." << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      coordWasSet=true;
      locationType=Image::Y;
      coordValue = atof(token);
      if ( coordValue < 0.0 || coordValue > m_global_ymax )
      {
	cerr << "Processing image command: " << "y value must be within the computational domain, not: " << coordValue << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
    }
    else if (startswith("z=", token))
    {
      token += 2; // skip z=
      if ( coordWasSet )
      {
	cerr << "Processing image command: " << "cannot set a coordinate location twice, z and " << locationType << " were both set." << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      coordWasSet=true;
      locationType=Image::Z;
      coordValue = atof(token);
      if ( coordValue < 0.0 || coordValue > m_global_zmax )
      {
	cerr << "Processing image command: " << "z value must be within the computational domain, not: " << coordValue << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
    }
    else
    {
      badOption("image", token);
    }
    token = strtok(NULL, " \t");
  }
  
  if ( !timingSet )
  {
    cerr << "Processing image command: " << "at least one timing mechanism must be set: cycle, time, cycleInterval or timeInterval" << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  
  // if topographic image, set flag
  // if (mode == Image::TOPO)
  // {
  //   mTopoImageFound = true;
  // }

  Image* i;

  // Set up the image object
  if (coordWasSet)
  {
    i = new Image(this, time, timeInterval, cycle, cycleInterval, 
                  filePrefix, sample, mode, locationType, coordValue, use_double);
    addImage(i);    
  }
  else 
  {
    cerr << "Processing image command: " << "one of the coordinate (x,y,z) option must be set to determine the image's 2D plane" << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
}

// int sgn(double arg)
// {
//   if (arg < 0)
//     return -1;
//   else 
//     return 1;
// }

//----------------------------------------------------------------------------
void EW::processSource(char* buffer, vector<Source*> & a_GlobalUniqueSources )
{

  Source* sourcePtr;
  
  double m0 = 1.0;
  double t0=0.0, f0=1.0, freq=1.0;
  // Should be center of the grid
  double x = 0.0, y = 0.0, z = 0.0;
  int i = 0, j = 0, k = 0;
  double mxx=0.0, mxy=0.0, mxz=0.0, myy=0.0, myz=0.0, mzz=0.0;
  double strike=0.0, dip=0.0, rake=0.0;
  double fx=0.0, fy=0.0, fz=0.0;
  int isMomentType = -1;
  
  double lat = 0.0, lon = 0.0, depth = 0.0;
  bool topodepth = false;
  
  bool cartCoordSet = false;
  bool geoCoordSet = false;
  bool strikeDipRake = false;

  int ncyc = 0;
  bool ncyc_set = false;

  timeDep tDep = iRickerInt;
  char formstring[100];
  strcpy(formstring, "Ricker");

  char* token = strtok(buffer, " \t");
  REQUIRE2(strcmp("source", token) == 0, "ERROR: not a source line...: " << token);
  token = strtok(NULL, " \t");

  string err = "Source Error: ";

  string cartAndGeoErr = "source command: Cannot set both a geographical (lat,lon,depth) and cartesian coordinate (x,y,z)";
  string pointAndMomentErr = "source command: Cannot set both a point source and moment tensor formulation";

  while (token != NULL)
    {
      // while there are tokens in the string still
       if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
       if (startswith("m0=", token) )
       {
          token += 3; // skip m0=
          CHECK_INPUT(atof(token) >= 0.0, 
                  err << "source command: scalar moment term must be positive, not: " << token);
          m0 = atof(token);
       }
       else if (startswith("x=", token))
       {
         CHECK_INPUT(!geoCoordSet, err << cartAndGeoErr);
         token += 2; // skip x=
         x = atof(token);
         cartCoordSet = true; 
      }
      else if (startswith("y=", token))
      {
         CHECK_INPUT(!geoCoordSet, err << cartAndGeoErr);
         token += 2; // skip y=
         y = atof(token);
         cartCoordSet = true;
      }
      else if (startswith("z=", token))
      {
         CHECK_INPUT(!geoCoordSet, err << cartAndGeoErr);
         token += 2; // skip z=
         CHECK_INPUT(atoi(token) >= 0.0, 
                 err << "source command: z coord must be zero or positive, not: " << token);
         z = atof(token);
         cartCoordSet = true;
      }
      else if (startswith("lat=", token))
      {
         CHECK_INPUT(!cartCoordSet, err << cartAndGeoErr);
         token += 4; // skip lat=
         lat = atof(token);
         CHECK_INPUT(lat >= -90.0,
                 "source command: lat must be greater than or equal to -90 degrees, not " 
                 << lat);
         CHECK_INPUT(lat <= 90.0,
                 "source command: lat must be less than or equal to 90 degrees, not "
                 << lat);
         geoCoordSet = true;
      }
      else if (startswith("lon=", token))
      {
         CHECK_INPUT(!cartCoordSet, err << cartAndGeoErr);
         token += 4; // skip lon=
         lon = atof(token);
         CHECK_INPUT(lon >= -180.0,
                 "source command: lon must be greater or equal to -180 degrees, not " 
                 << lon);
         CHECK_INPUT(lon <= 180.0,
                 "source command: lon must be less than or equal to 180 degrees, not "
                 << lon);
         geoCoordSet = true;
      }
      else if (startswith("depth=", token))
      {
         CHECK_INPUT(!cartCoordSet, err << cartAndGeoErr);
         token += 6; // skip depth=
         depth = atof(token);
	 topodepth = false;
// by depth we mean depth below mean sea level, i.e., z-coordinate. With topography, some sources can have negative z-coordinate
         geoCoordSet = true;
      }
//                         1234567890
      else if (startswith("topodepth=", token))
      {
         CHECK_INPUT(!cartCoordSet, err << cartAndGeoErr);
         token += 10; // skip depth=
         depth = atof(token);
	 topodepth = true;
         CHECK_INPUT(depth >= 0.0,
                 err << "source command: Depth below topography must be greater than or equal to zero");
// by depth we here mean depth below topography
         geoCoordSet = true;
      }
      else if (startswith("Mxx=", token) || startswith("mxx=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Mxx=
         mxx = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Mxy=", token) || startswith("mxy=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Mxy=
         mxy = atof(token);
	  isMomentType = 1;
      }
      else if (startswith("Mxz=", token) || startswith("mxz=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Mxz=
         mxz = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Myy=", token) || startswith("myy=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Myy=
         myy = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Myz=", token) || startswith("myz=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Myz=
         myz = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Mzz=", token) || startswith("mzz=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Mzz=
         mzz = atof(token);
         isMomentType = 1;
      }
      else if (startswith("Fz=", token) || startswith("fz=", token))
      {
         CHECK_INPUT(isMomentType != 1, err << pointAndMomentErr);
         token += 3; // skip Fz=
         fz = atof(token);
         isMomentType = 0;
      }
      else if (startswith("Fx=", token) || startswith("fx=", token))
      {
         CHECK_INPUT(isMomentType != 1, err << pointAndMomentErr);
         token += 3; // skip Fx=
         fx = atof(token);
         isMomentType = 0;
      }
      else if (startswith("Fy=", token) || startswith("fy=", token))
      {
         CHECK_INPUT(isMomentType != 1, err << pointAndMomentErr);
         token += 3; // skip Fy=
         fy = atof(token);
         isMomentType = 0;
      }
      else if (startswith("Rake=", token) || startswith("rake=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 5; // skip Rake=
         rake = atof(token);
	 strikeDipRake = true;
         isMomentType = 1;
      }
      else if (startswith("Strike=", token) || startswith("strike=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 7; // skip Strike=
         strike = atof(token);
	 strikeDipRake = true;
         isMomentType = 1;
      }
      else if (startswith("Dip=", token) || startswith("dip=", token))
      {
         CHECK_INPUT(isMomentType != 0, err << pointAndMomentErr);
         token += 4; // skip Dip=
         dip = atof(token);
	 strikeDipRake = true;
         isMomentType = 1;
      }
      else if (startswith("t0=", token))
      {
         token += 3; // skip t0=
         t0 = atof(token);
      }
      else if (startswith("freq=", token))
      {
         token += 5; // skip freq=
         freq = atof(token);
         CHECK_INPUT(freq > 0,
                 err << "source command: Frequency must be > 0");
      }
      else if (startswith("amp=", token) || startswith("f0=", token))
      {
         CHECK_INPUT(isMomentType != 1,
                 err << "source command: Cannot set force amplitude for moment tensor terms");
	 if (startswith("amp=", token))
	   {
	     deprecatedOption("source","amp","f0");
	     token += strlen("amp=");
	   }
	 else
	   token += strlen("f0=");
         f0 = atof(token);
      }
      else if (startswith("type=",token))
      {
         token += 5;
         strcpy(formstring, token);
         if (!strcmp("Ricker",formstring))
            tDep = iRicker;
         else if (!strcmp("Gaussian",formstring))
            tDep = iGaussian;
         else if (!strcmp("Ramp",formstring))
            tDep = iRamp;
         else if (!strcmp("Triangle",formstring))
            tDep = iTriangle;
         else if (!strcmp("Sawtooth",formstring))
            tDep = iSawtooth;
         else if (!strcmp("SmoothWave",formstring))
            tDep = iSmoothWave;
         else if (!strcmp("Erf",formstring) || !strcmp("GaussianInt",formstring) )
            tDep = iErf;
         else if (!strcmp("VerySmoothBump",formstring))
            tDep = iVerySmoothBump;
         else if (!strcmp("RickerInt",formstring) )
            tDep = iRickerInt;
         else if (!strcmp("Brune",formstring) )
	    tDep = iBrune;
         else if (!strcmp("BruneSmoothed",formstring) )
	    tDep = iBruneSmoothed;
         else if (!strcmp("DBrune",formstring) )
	    tDep = iDBrune;
         else if (!strcmp("GaussianWindow",formstring) )
	    tDep = iGaussianWindow;
         else if (!strcmp("Liu",formstring) )
	    tDep = iLiu;
	 else
            if (m_myRank == 0)
	      cout << "unknown time function: " << formstring << endl << " using default RickerInt function." << endl;
      }
      else if (startswith("ncyc=", token))
      {
         token += 5; // skip ncyc=
         ncyc = atoi(token);
         CHECK_INPUT(ncyc > 0,
                 err << "source command: Number of cycles must be > 0");
         ncyc_set = true;
      }
      else
      {
         badOption("source", token);
      }
      token = strtok(NULL, " \t");
    }

  // Set up source on wpp object.
  CHECK_INPUT(cartCoordSet || geoCoordSet,
	  err << "source command: cartesian or geographic coordinate must be specified");

  if( tDep == iGaussianWindow )
      CHECK_INPUT( ncyc_set, err << "source command: ncyc must be set for Gaussian Window function");

  // --------------------------------------------------------------------------- 
  // However the location for the source was specified, we are going to
  // find the grid points associated with the location. (i.e., assign
  // i, j, k to valid values)
  // --------------------------------------------------------------------------- 
  if (geoCoordSet)
  {
    computeCartesianCoord(x, y, lon, lat);
    cartCoordSet = true;
  }
  if (cartCoordSet)
  {
    double xmin = 0.;
    double ymin = 0.;
    double zmin;

// only check the z>zmin when we have topography. For a flat free surface, we will remove sources too 
// close or above the surface in the call to mGlobalUniqueSources[i]->correct_Z_level()
    if (topographyExists()) // topography command must be read before the source command
      zmin = m_global_zmin;
    else
      zmin = -m_global_zmax;

    if (x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax || z < zmin || z > m_global_zmax)
    {
      stringstream sourceposerr;
      sourceposerr << endl
		   << "***************************************************" << endl
		   << " FATAL ERROR:  Source positioned outside grid!  " << endl
		   << endl
		   << " Source Type: " << formstring << endl
		   << "              @ x=" << x 
		   << " y=" << y << " z=" << z << endl 
		   << endl;
	    
      if ( x < xmin )
	sourceposerr << " x is " << xmin - x << 
	  " meters away from min x (" << xmin << ")" << endl;
      else if ( x > m_global_xmax)
	sourceposerr << " x is " << x - m_global_xmax << 
	  " meters away from max x (" << m_global_xmax << ")" << endl;
      if ( y < ymin )
	sourceposerr << " y is " << ymin - y << 
	  " meters away from min y (" << ymin << ")" << endl;
      else if ( y > m_global_ymax)
	sourceposerr << " y is " << y - m_global_ymax << 
	  " meters away from max y (" << m_global_ymax << ")" << endl;
      if ( z < zmin )
	sourceposerr << " z is " << zmin - z << 
	  " meters away from min z (" << zmin << ")" << endl;
      else if ( z > m_global_zmax)
	sourceposerr << " z is " << z - m_global_zmax << 
	  " meters away from max z (" << m_global_zmax << ")" << endl;
      sourceposerr << "***************************************************" << endl;
      if (m_myRank == 0)
	cout << sourceposerr.str();
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  // if strike, dip and rake have been given we need to convert into M_{ij} form
  if ( strikeDipRake )
    {
      double radconv = M_PI / 180.;
      double S, D, R;
      strike -= mGeoAz; // subtract off the grid azimuth
      S = strike*radconv; D = dip*radconv; R = rake*radconv;
      
      mxx = -1.0 * ( sin(D) * cos(R) * sin (2*S) + sin(2*D) * sin(R) * sin(S)*sin(S) );
      myy =        ( sin(D) * cos(R) * sin (2*S) - sin(2*D) * sin(R) * cos(S)*cos(S) );
      mzz = -1.0 * ( mxx + myy );	
      mxy =        ( sin(D) * cos(R) * cos (2*S) + 0.5 * sin(2*D) * sin(R) * sin(2*S) );
      mxz = -1.0 * ( cos(D) * cos(R) * cos (S)   + cos(2*D) * sin(R) * sin(S) );
      myz = -1.0 * ( cos(D) * cos(R) * sin (S)   - cos(2*D) * sin(R) * cos(S) );
    }
  
  if (isMomentType)
    {
      // these have global location since they will be used by all processors
      sourcePtr = new Source(this, m0, freq, t0, x, y, z, mxx, mxy, mxz, myy, myz, mzz,
                             tDep, formstring, ncyc);
// relative depth?
      sourcePtr->set_z_is_relative_to_topography( topodepth );
      
//      addGlobalSorceTerm(sourcePtr);
      a_GlobalUniqueSources.push_back(sourcePtr);
    }
  else // point forcing
    {
      // global version (gets real coordinates)
      sourcePtr = new Source(this, f0, freq, t0, x, y, z, fx, fy, fz, tDep, formstring, ncyc);
// relative depth?
      sourcePtr->set_z_is_relative_to_topography( topodepth );
      //...and add it to the list of forcing terms
//      addGlobalSorceTerm(sourcePtr);
      a_GlobalUniqueSources.push_back(sourcePtr);
    }	  
}

//------------------------------------------------------------------------
void EW::processMaterialBlock( char* buffer, int & blockCount )
{
  double vpgrad=0.0, vsgrad=0.0, rhograd=0.0;
  bool x1set=false, x2set=false, y1set=false, y2set=false, 
    z1set=false, z2set=false;

  double x1=0.0, x2=0.0, y1=0.0, y2=0.0, z1=0.0, z2=0.0;
  int i1=-1, i2=-1, j1=-1, j2=-1, k1=-1, k2=-1;

  string name = "Block";

  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("block", token) == 0,
	      "ERROR: material block can be set by a block line, not: " << token);

  string err = token;
  err += " Error: ";

  token = strtok(NULL, " \t");

  double vp=-1, vs=-1, rho=-1, qp=-1, qs=-1, freq=1;
  bool absDepth=false;

  while (token != NULL)
    {
      // while there are tokens in the string still
       if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
       if (startswith("vp=", token) )
      {
         token += 3; // skip vp=
         vp = atof(token);
      }
      else if (startswith("vs=", token) )
      {
         token += 3; // skip vs=
         vs = atof(token);
      }
      else if (startswith("r=", token)) // superseded by rho=, but keep for backward compatibility
      {
         token += 2; // skip r=
         rho = atof(token);
      }
      else if (startswith("rho=", token))
      {
         token += 4; // skip rho=
         rho = atof(token);
      }
      else if (startswith("rhograd=", token))
      {
         token += 8; // skip rhograd=
         rhograd = atof(token);
      }
      else if (startswith("vpgrad=", token))
      {
         token += 7; // skip vpgrad=
         vpgrad = atof(token);
      }
      else if (startswith("vsgrad=", token))
      {
         token += 7; // skip vsgrad=
         vsgrad = atof(token);
      }
      else if (startswith("Qs=", token) || startswith("qs=",token) )
      {
         token += 3; // skip qs=
         qs = atof(token);
      }
      else if (startswith("Qp=", token) || startswith("qp=",token) )
      {
         token += 3; // skip qp=
         qp = atof(token);
      }
//                         1234567890
      else if (startswith("absdepth=", token) )
      {
	token += 9; // skip absdepth=
	absDepth = (bool) atoi(token);
      }
      else if (startswith("x1=", token))
      {
         token += 3; // skip x1=
         x1 = atof(token);
         x1set = true;
      }
      else if (startswith("x2=", token))
      {
         token += 3; // skip x2=
         x2 = atof(token);
         x2set = true;
      }
      else if (startswith("y1=", token))
      {
         token += 3; // skip y1=
         y1 = atof(token);
         y1set = true;
      }
      else if (startswith("y2=", token))
      {
         token += 3; // skip y2=
         y2 = atof(token);
         y2set = true;
      }
      else if (startswith("z1=", token))
      {
         token += 3; // skip z1=
         z1 = atof(token);
         z1set = true;
      }
      else if (startswith("z2=", token))
      {
         token += 3; // skip z2=
         z2 = atof(token);
         z2set = true;
      }
      else
      {
         badOption("block", token);
      }
      token = strtok(NULL, " \t");
    }
  // End parsing...
  
  blockCount++;
  stringstream blockname;
  blockname << name << " " << blockCount;
  name = blockname.str();

  // Set up a block on wpp object.

  if (x1set)
  {
     // CHECK_INPUT(x1 >= 0.,
     // 	     err << "x1 is less than the minimum x, " 
     // 	     << x1 << " < " << 0.);
     CHECK_INPUT(x1 <= m_global_xmax,
	     err << "x1 is greater than the maximum x, " 
	     << x1 << " > " << m_global_xmax);
  }
  else
    x1 = -m_global_xmax; //x1 = 0.;

  if (x2set)
  {
     CHECK_INPUT(x2 >= 0.,
             err << "x2 is less than the minimum x, " 
             << x2 << " < " << 0.);
     // CHECK_INPUT(x2 <= m_global_xmax,
     //         err << "x2 is greater than the maximum x, " 
     //         << x2 << " > " << m_global_xmax);
  }
  else
    x2 = 2.*m_global_xmax;//x2 = m_global_xmax;

  CHECK_INPUT( x2 >= x1, " (x1..x2), upper bound is smaller than lower bound");
  
  //--------------------------------------------------------
  // Set j bounds, goes with Y in WPP
  //--------------------------------------------------------
  if (y1set)
  {
     // CHECK_INPUT(y1 >= 0.,
     // 	     err << "y1 is less than the minimum y, " << y1 << " < " << 0.);

     CHECK_INPUT(y1 <= m_global_ymax,
		  err << "y1 is greater than the maximum y, " << y1 << " > " << m_global_ymax);
  }
  else
    y1 = -m_global_ymax;//y1 = 0.;
      
  if (y2set)
  {
     CHECK_INPUT(y2 >= 0.,
	     err << "y2 is less than the minimum y, " << y2 << " < " << 0.);
  }
  else
    y2 = 2.*m_global_ymax;//y2 = m_global_ymax;

  CHECK_INPUT( y2 >= y1, " (y1..y2), upper bound is smaller than lower bound");

  if (z1set)
  {
    // CHECK_INPUT(topographyExists() || z1 >= 0.,
    //         err << "z1 is less than the minimum z, " << z1 << " < " << 0.);
    CHECK_INPUT(z1 <= m_global_zmax, 
            err << "z1 is greater than the maximum z, " << z1 << " > " << m_global_zmax);
  }
  else
    z1 = -m_global_zmax;//z1 = 0.;

  if (z2set)
  {
    CHECK_INPUT(topographyExists() || z2 >= 0.,
            err << "z2 is less than the minimum z, " << z2 << " < " << 0.);
    // CHECK_INPUT(z2 <= m_global_zmax,
    // 		err << "z2 is greater than the maximum z, " << z2 << " > " << m_global_zmax);
  }
  else
    z2 = 2.*m_global_zmax;//z2 = m_global_zmax;

  CHECK_INPUT( z2 >= z1, " (z1..z2), upper bound is smaller than lower bound");

  if(getVerbosity() >=2 &&  m_myRank == 0 )
     cout << name << " has bounds " << x1 << " " << x2 << " " << y1 << " "
	  << y2 << " " << z1 << " " << z2 << endl;

  CHECK_INPUT( vs > 0 && vp > 0 && rho > 0 , "Error in block " << name << " vp vs rho are   "
	       << vs << " " << vp << " " << rho );

  MaterialBlock* bl = new MaterialBlock( this ,rho, vs, vp, x1, x2, y1, y2, z1, z2, qs, qp, freq );
  bl->set_gradients( rhograd, vsgrad, vpgrad );
  bl->set_absoluteDepth( absDepth );
  add_mtrl_block( bl );
}

void EW::processReceiver(char* buffer, vector<TimeSeries*> & a_GlobalTimeSeries)
{
  double x=0.0, y=0.0, z=0.0;
  double lat = 0.0, lon = 0.0, depth = 0.0;
  bool cartCoordSet = false, geoCoordSet = false;
  string name = "rec";
  int writeEvery = 1000;

  bool dateSet = false;
  bool timeSet = false;
  bool topodepth = false;

  string date = "";
  string time = "";

  bool usgsformat = 0, sacformat=0; // default is now to not write any file: you have to specify a format
  TimeSeries::receiverMode mode=TimeSeries::Displacement;

  char* token = strtok(buffer, " \t");
//  int nsew=0, vel=0;

// tmp
//  cerr << "******************** INSIDE process receiver *********************" << endl;

  CHECK_INPUT(strcmp("receiver", token) == 0, "ERROR: not a receiver line...: " << token);
  token = strtok(NULL, " \t");

  string err = "RECEIVER Error: ";

  //  cout << "start receiver " << m_myRank << " " << token <<"x"<<endl;
  while (token != NULL)
  {
     // while there are tokens in the string still
     //     cout << m_myRank << " token " << token <<"x"<<endl;

     if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
        break;
     if (startswith("x=", token))
     {
        CHECK_INPUT(!geoCoordSet,
                err << "receiver command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 2; // skip x=
        cartCoordSet = true;
        x = atof(token);
        CHECK_INPUT(x >= 0.0,
		    "receiver command: x must be greater than or equal to 0, not " << x);
        CHECK_INPUT(x <= m_global_xmax,
		    "receiver command: x must be less than or equal to xmax, not " << x);
     }
     else if (startswith("y=", token))
     {
        CHECK_INPUT(!geoCoordSet,
                err << "receiver command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 2; // skip y=
        cartCoordSet = true;
        y = atof(token);
        CHECK_INPUT(y >= 0.0,
                "receiver command: y must be greater than or equal to 0, not " << y);
        CHECK_INPUT(y <= m_global_ymax,
		    "receiver command: y must be less than or equal to ymax, not " << y);
     }
     else if (startswith("lat=", token))
     {
        CHECK_INPUT(!cartCoordSet,
                err << "receiver command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 4; // skip lat=
        lat = atof(token);
        CHECK_INPUT(lat >= -90.0,
                "receiver command: lat must be greater than or equal to -90 degrees, not " 
                << lat);
        CHECK_INPUT(lat <= 90.0,
                "receiver command: lat must be less than or equal to 90 degrees, not "
                << lat);
        geoCoordSet = true;
     }
     else if (startswith("lon=", token))
     {
        CHECK_INPUT(!cartCoordSet,
                err << "receiver command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 4; // skip lon=
        lon = atof(token);
        CHECK_INPUT(lon >= -180.0,
                "receiver command: lon must be greater or equal to -180 degrees, not " 
                << lon);
        CHECK_INPUT(lon <= 180.0,
                "receiver command: lon must be less than or equal to 180 degrees, not "
                << lon);
        geoCoordSet = true;
     }
     else if (startswith("z=", token))
     {
       token += 2; // skip z=
// depth is currently the same as z
       depth = z = atof(token);
       topodepth = false;
       CHECK_INPUT(z <= m_global_zmax,
		   "receiver command: z must be less than or equal to zmax, not " << z);
     }
     else if (startswith("depth=", token))
     {
        token += 6; // skip depth=
       z = depth = atof(token);
       topodepth = true;
       CHECK_INPUT(depth >= 0.0,
	       err << "receiver command: depth must be greater than or equal to zero");
       CHECK_INPUT(depth <= m_global_zmax,
		   "receiver command: depth must be less than or equal to zmax, not " << depth);
// by depth we here mean depth below topography
     }
     else if(startswith("file=", token))
     {
        token += 5; // skip file=
        name = token;
     }
     // else if( startswith("nsew=", token) )
     // {
     //    token += strlen("nsew=");
     //    nsew = atoi(token);
     // }
     // else if( startswith("velocity=", token) )
     // {
     //    token += strlen("velocity=");
     //    vel = atoi(token);
     // }
     else if (startswith("writeEvery=", token))
     {
       token += strlen("writeEvery=");
       writeEvery = atoi(token);
       CHECK_INPUT(writeEvery >= 0,
	       err << "sac command: writeEvery must be set to a non-negative integer, not: " << token);
     }
     else if( startswith("usgsformat=", token) )
     {
        token += strlen("usgsformat=");
        usgsformat = atoi(token);
     }
     else if( startswith("sacformat=", token) )
     {
        token += strlen("sacformat=");
        sacformat = atoi(token);
     }
     else if( startswith("variables=", token) )
     {
       token += strlen("variables=");
       if( strcmp("displacement",token)==0 )
       {
	 mode = TimeSeries::Displacement;
       }
       else if( strcmp("velocity",token)==0 )
       {
	 mode = TimeSeries::Velocity;
       }
       else if( strcmp("div",token)==0 )
       {
	 mode = TimeSeries::Div;
       }
       else if( strcmp("curl",token)==0 )
       {
	 mode = TimeSeries::Curl;
       }
       else if( strcmp("strains",token)==0 )
       {
	 mode = TimeSeries::Strains;
       }
       else
       {
	 if (proc_zero())
	   cerr << "receiver command: variables=" << token << " not understood" << endl
		<< "using default mode (displacement)" << endl << endl;
	 mode = TimeSeries::Displacement;
       }
       
     }
     else
     {
        badOption("receiver", token);
     }
     token = strtok(NULL, " \t");
  }  
  //  cout << "end receiver " << m_myRank << endl;

  if (geoCoordSet)
  {
    computeCartesianCoord(x, y, lon, lat);
// check if (x,y) is within the computational domain
  }

  bool inCurvilinear=false;
// we are in or above the curvilinear grid 
  if ( topographyExists() && z < m_zmin[mNumberOfCartesianGrids-1])
  {
    inCurvilinear = true;
  }
      
// check if (x,y,z) is not in the global bounding box
    if ( !( (inCurvilinear || z >= 0) && x>=0 && x<=m_global_xmax && y>=0 && y<=m_global_ymax))
    {
// The location of this station was outside the domain, so don't include it in the global list
    if (m_myRank == 0 && getVerbosity() > 0)
    {
      stringstream receivererr;
  
      receivererr << endl 
	     << "***************************************************" << endl
	     << " WARNING:  RECEIVER positioned outside grid!" << endl;
      receivererr << " No RECEIVER file will be generated for file = " << name << endl;
      if (geoCoordSet)
      {
	receivererr << " @ lon=" << lon << " lat=" << lat << " depth=" << depth << endl << endl;
      }
      else
      {
	receivererr << " @ x=" << x << " y=" << y << " z=" << z << endl << endl;
      }
      
      receivererr << "***************************************************" << endl;
      cerr << receivererr.str();
      cerr.flush();
    }
  }
  else
  {
    TimeSeries *ts_ptr = new TimeSeries(this, name, mode, sacformat, usgsformat, x, y, depth, 
					topodepth, writeEvery);
// include the receiver in the global list
    a_GlobalTimeSeries.push_back(ts_ptr);
    
  }
    
}

