//-*-c++-*-
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
#include "mpi.h"

#include "EW.h"

#include "version.h"
#include "Require.h"
#include "nearlyEqual.h"
#include "boundaryConditionTypes.h"
#include "AnisotropicMaterialBlock.h"
#include "MaterialBlock.h"
#include "MaterialPfile.h"
#include "MaterialIfile.h"
#include "MaterialVolimagefile.h"
#include "MaterialRfile.h"
#include "MaterialSfile.h"
#include "MaterialInvtest.h"
#include "EtreeFile.h"
#include "TimeSeries.h"
#include "Filter.h"
#include "Image3D.h"
#include "ESSI3D.h"
#include "SfileOutput.h"
#include "sacutils.h"
//#include "TestGrid.h"
#include "GridGeneratorGaussianHill.h"
#include "GridGeneratorGeneral.h"

#if USE_HDF5
#include "readhdf5.h"
#endif

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <time.h>

using namespace std;

#define SQR(x) ((x)*(x))

//int gcd( int a, int b )
//{
//   // Euclidean algorithm
//   while( b != 0 )
//   {
//      int t = b;
//      b = a % b;
//      a = t;
//   }
//   return a;
//}

void EW::revvector( int npts, float_sw4* v )
{
   for( int i=0 ; i < npts/2; i++ )
   {
      float_sw4 sl=v[i];
      v[i]=v[npts-1-i];
      v[npts-1-i]=sl;
   }
}

int computeEndGridPoint( float_sw4 maxval, float_sw4 dh )
{
  // We round up one, so that the end point
  // specified by the user is always included
  // in the domain.  i.e. if z was specified
  // as 15.0, and dh was computed to be 3.33,
  // the number of grid points would be 15.0/3.33 + 1
  // or 5, giving us a new max z = 16.64.
  int pts = 0;
  float_sw4 x = 0.0;

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

//-----------------------------------------------------------------------
//bool endswith(string end, string& mystr)
//{
//   int lenEnd = end.length();
//   int lenStr = mystr.length();

//   if (lenEnd > lenStr) return false;

//   cout << "mystr: " << mystr << " end: " << end << " " << mystr.substr(lenStr-lenEnd, lenEnd) << endl;

//   if (mystr.substr(lenStr-lenEnd, lenEnd) == end)
//      return true;
//   else 
//      return false;
//}

//-----------------------------------------------------------------------
bool EW::startswith(const char begin[], char *line)
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

//-----------------------------------------------------------------------
void EW::deprecatedOption(const string& command, 
		      const string& oldone, 
		      const string& newone)
{
  if (m_myRank == 0)
    cout << "DeprecationWarning: " 
	 << command << " option " << oldone << " is no longer supported.  Use "
	 << newone << " instead." << endl;
}

//void unchecked(const char* cmd)
//{
//   cout << "*** Not yet checking command: " << cmd << endl;
//}

//void checked(const char* cmd)
//{
//   cout << "*** " << cmd << " command checked! " << endl;
//}

//-----------------------------------
// 
// Note that parseInputFile() calls a lot of member functions of the EW class that
// should not be called after the initialization of the EW object is completed. 
// Make all these functions private!
//
bool EW::parseInputFile( vector<vector<Source*> > & a_GlobalUniqueSources,
			 vector< vector<TimeSeries*> > & a_GlobalTimeSeries )
{
  char buffer[256];
  ifstream inputFile;
  int blockCount=0;
  int ablockCount=0;

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

// tmp (the fileio command has not yet been parsed, so we don't know mVerbose
//  cout << "********Reading the input file, proc=" << m_myRank << endl;

// First process Geodyn input for restrictions of allowable grid sizes.
 while (!inputFile.eof())
 {
    inputFile.getline(buffer, 256);
    if( startswith("geodynbc",buffer ) )
       geodynFindFile(buffer);
 }
 inputFile.clear();
 inputFile.seekg(0, ios::beg);

// process the testrayleigh command to enable a periodic domain in the (x,y)-directions
// these commands can enter data directly the object (this->)
  while (!inputFile.eof())
  {    
     inputFile.getline(buffer, 256);
     if (startswith("testrayleigh", buffer) )
     {
       m_doubly_periodic = true;
     }
     else if( startswith("testenergy",buffer) )
     {
	m_doubly_periodic = checkTestEnergyPeriodic(buffer);
     }
     else if (startswith("refinement",buffer) )
     {
	// mesh refinements require 3 ghost points, must know 
	// before processing grid command.
	m_mesh_refinements = true;
     }
     else if( startswith("supergrid",buffer) )
     {
	// If supergrid damping is 6th order, 3 ghost points are needed, must know 
	// before processing grid command.
	processSupergrid(buffer);
     }
     else if( startswith("developer",buffer) )
	processDeveloper(buffer); // Need this early to determine array index order before any arrays are used.
  }

  inputFile.clear();
  inputFile.seekg(0, ios::beg); // reset file pointer to the beginning of the input file

//---------------------------------------------------------------
// Then process the grid, fileio, and topography commands so
// we know how big the solution arrays need to be.
//
// Also, if we are using attenuation, enable it on now so it
// can be read in with the other material properties
//---------------------------------------------------------------

// these commands can enter data directly into the object (this->)
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
     else if (startswith("refinement", buffer))
     {
	processRefinement(buffer);
     }
     else if (startswith("topography", buffer))
     {
        processTopography(buffer);
     }
     else if (startswith("attenuation", buffer))
     {
        processAttenuation(buffer);
     }
     else if (startswith("anisotropy", buffer))
     {
        m_anisotropic = true;
     }
     else if (startswith("time", buffer))
     {
        processTime(buffer); // process time command to set reference UTC before reading stations.
     }
     else if (startswith("prefilter", buffer))
     {
       // before reading any rupture command, we need to know 
       // if they need to be prefiltered  
       processPrefilter(buffer);
     }
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

  if( m_anisotropic && m_use_attenuation )
  {
    if (m_myRank == 0)
    {
      cerr << "Error: Attenuation not implemented with anisotropy " << endl;
      return false; // unsuccessful
    }
  }  

//  if( m_mesh_refinements && (m_anisotropic || (m_use_attenuation && m_number_mechanisms>0) ) )
  if( m_mesh_refinements && m_anisotropic )
  {
    if (m_myRank == 0)
    {
//      cerr << "Error: Grid refinements not implemented with attenuation or anisotropy " << endl;
      cerr << "Error: Grid refinements not implemented with anisotropy " << endl;
      return false; // unsuccessful
    }
  }

// sort and correct vector 'm_refinementBoundaries'. Initialize if not already available
  cleanUpRefinementLevels();
  
  inputFile.clear();
  inputFile.seekg(0, ios::beg); // reset file pointer to the beginning of the input file

// At this point we only allocate solution arrays for the Cartesian grids 
// Need to read the topography information before we can decide on sizes for the
// curvilinear grid.
  allocateCartesianSolverArrays(m_global_zmax); 

// setup 2D communicators on the finest grid so that we can smooth the topography
  setup2D_MPICommunications();

// deal with topography
  if (m_topography_exists)
  {
// 1. read topography from efile 
     if (m_topoInputStyle == EW::Efile)
     {
 	extractTopographyFromEfile(m_topoFileName, m_topoExtFileName, m_QueryType,
				   m_EFileResolution);
     }
     else if (m_topoInputStyle == EW::GridFile)
     {
 	extractTopographyFromGridFile(m_topoFileName);
     }
     else if (m_topoInputStyle == EW::CartesianGrid)
     {
 	extractTopographyFromCartesianFile(m_topoFileName);
     }
     else if (m_topoInputStyle == EW::TopoImage)
     {
 	extractTopographyFromImageFile(m_topoFileName);
     }
     else if (m_topoInputStyle == EW::GaussianHill) // assumed to populate all grid points
     {
        m_gridGenerator->fill_topo( mTopo, mGridSize[mNumberOfGrids-1] );
        m_gridGenerator->fill_topo( mTopoGridExt, mGridSize[mNumberOfGrids-1] );
        // 	buildGaussianHillTopography(m_GaussianAmp, m_GaussianLx, m_GaussianLy, m_GaussianXc, m_GaussianYc);
     }      
     else if( m_topoInputStyle == EW::Rfile )
	extractTopographyFromRfile( m_topoFileName );
     else if( m_topoInputStyle == EW::Sfile )
	extractTopographyFromSfile( m_topoFileName );

// preprocess the mTopo array
     if (m_topoInputStyle != EW::GaussianHill) // no smoothing or extrapolation for a gaussian hill
     {
// 1. fill in any undefined ghost point values by extrapolation
	extrapolateTopo(mTopo);
// 2. check that all values are defined...
	checkTopo(mTopo);
// 3. smooth the topo
 	smoothTopography(m_maxIter);

// Assign interface surfaces (needed when there is MR in the curvilinear portion of the grid)
        m_gridGenerator->assignInterfaceSurfaces( this, mTopoGridExt );
     }
     
// // 3. Figure out the number of grid points in the vertical direction and allocate solution arrays on the curvilinear grid
     allocateCurvilinearArrays(); // need to assign  m_global_nz[g] = klast - m_ghost_points; + allocate mUacc
  }
  else
  {
     if (proc_zero_evzero())
	cout << endl << 
	   "*** No topography command found in input file. Using z=0 as free surface boundary ***" << endl << endl;
  }

// setup communicators for 3D solutions on all grids
  setupMPICommunications();


// Make curvilinear grid and compute metric
  for( int g=mNumberOfCartesianGrids ; g < mNumberOfGrids ; g++ )
     m_gridGenerator->generate_grid_and_met( this, g, mX[g], mY[g], mZ[g], mJ[g], mMetric[g] );

  //  if (m_topography_exists)
  //  {
  //
  //     if( m_topoInputStyle == EW::GaussianHill && mNumberOfGrids-mNumberOfCartesianGrids > 1 )
  //     {
  //        TestGrid* gh = create_gaussianHill();
  //        for( int g=mNumberOfCartesianGrids ; g < mNumberOfGrids ; g++ )
  //           gh->generate_grid_and_met( this, g, mX[g], mY[g], mZ[g], mJ[g], mMetric[g] );
  //        delete gh;
  //     }
  //     else
  //     {
  //        generate_grid();
  //        setup_metric();
  //     }
  //  }

// output grid size info
  if (proc_zero_evzero())
  {
    int nx, ny, nz;
    double nTot=0.;
    printf("\nGlobal grid sizes (without ghost points)\n");
//             1234  12345679  12345679  12345679  12345679
    printf("Grid         h        Nx        Ny        Nz       Points      Type\n");
    for (int g = 0; g < mNumberOfGrids; g++)
    {
      nx = m_global_nx[g];
      ny = m_global_ny[g];
      nz = m_kEnd[g] - m_ghost_points;
      nTot += ((long long int)nx)*ny*nz;
      printf("%4i %9g %9i %9i %9i %12lld     %s\n", g, mGridSize[g], nx, ny, nz, ((long long int)nx)*ny*nz,
             (g < mNumberOfCartesianGrids) ? "Cartesian": "Curvilinear");
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
	   startswith("anisotropy", buffer) || 
	   startswith("fileio", buffer) ||
	   startswith("supergrid", buffer) ||
	   startswith("prefilter", buffer) ||
	   startswith("developer", buffer) ||
	   startswith("time", buffer) ||
// ignore material optimizer commands
 	   startswith("event", buffer) ||
 	   startswith("mparcart", buffer) ||
	   startswith("mpallpts", buffer) ||
 	   startswith("mrun", buffer) ||
 	   startswith("mscalefactors", buffer) ||
 	   startswith("lbfgs", buffer) ||
 	   startswith("nlcg", buffer) ||
 	   startswith("mfsurf", buffer) ||
 	   startswith("mimage", buffer) ||	   	   
 	   startswith("m3dimage", buffer) ||	   	   
 	   startswith("regularize", buffer) ||	   	   
 	   startswith("mtypx", buffer) ||
	   startswith("\n", buffer) || startswith("\r", buffer) )
// || startswith("\r", buffer) || startswith("\0", buffer))
       {
// Ignore commented lines, newlines,
// grid, fileio, and topography, since we have already processed those commands.
       }
       else if (startswith("gmt", buffer))
         processGMT(buffer);
       else if (startswith("checkpoint",buffer))
	  processCheckPoint(buffer);
       else if (startswith("globalmaterial", buffer))
         processGlobalMaterial(buffer);
       else if (!m_inverse_problem && (startswith("rechdf5", buffer) || startswith("sachdf5", buffer)) ) // was called "sac" in WPP
	 processReceiverHDF5(buffer, a_GlobalTimeSeries);
       else if (!m_inverse_problem && (startswith("rec", buffer) || startswith("sac", buffer)) ) // was called "sac" in WPP
	 processReceiver(buffer, a_GlobalTimeSeries);
       else if (m_inverse_problem && (startswith("obshdf5", buffer) || startswith("observationhdf5", buffer))) // 
	  processObservationHDF5(buffer, a_GlobalTimeSeries);
       else if (m_inverse_problem && startswith("obs", buffer)) // 
	  processObservation(buffer, a_GlobalTimeSeries);
       else if (m_inverse_problem && startswith("scalefactors", buffer)) // 
	  processScaleFactors(buffer);
       else if (m_inverse_problem && startswith("cg", buffer)) // 
	  processCG(buffer);
       // else if (startswith("energy", buffer))
       //   processEnergy(buffer);
       else if (startswith("twilight", buffer))
	 processTwilight(buffer);
       else if (startswith("testpointsource", buffer))
	 processTestPointSource(buffer);
       else if (startswith("testlamb", buffer))
         processTestLamb(buffer);
       else if (startswith("testrayleigh", buffer))
         processTestRayleigh(buffer);
       else if (startswith("testenergy", buffer))
         processTestEnergy(buffer);
       else if (startswith("source", buffer))
	 processSource(buffer, a_GlobalUniqueSources);
       else if (startswith("rupturehdf5", buffer))
	 processRuptureHDF5(buffer, a_GlobalUniqueSources);
       else if (startswith("rupture", buffer))
	 processRupture(buffer, a_GlobalUniqueSources);
       else if (startswith("block", buffer))
	 processMaterialBlock(buffer, blockCount);
       else if (startswith("ablock", buffer) && m_anisotropic )
	 processAnisotropicMaterialBlock(buffer, ablockCount);
       else if (startswith("pfile", buffer))
	 processMaterialPfile( buffer );
       else if (startswith("rfile", buffer))
	 processMaterialRfile( buffer );
       else if (startswith("sfileoutput", buffer))
          processSfileOutput(buffer);
       else if (startswith("sfile", buffer))
	 processMaterialSfile( buffer );
       else if (startswith("vimaterial", buffer))
	 processMaterialVimaterial( buffer );
       else if (startswith("invtestmaterial", buffer))
	  processMaterialInvtest(buffer);
       else if (startswith("efile", buffer))
       {
#ifndef ENABLE_ETREE
          if (m_myRank==0) 
             cout << "Error: SW4 was not built with Etree (efile) support" << endl;
          return false;
#endif
	  processMaterialEtree(buffer);
       }
       else if (startswith("ifile", buffer))
          processMaterialIfile(buffer);
       else if (startswith("material", buffer))
          processMaterial(buffer);
       else if (startswith("imagehdf5", buffer))
         processImage(buffer, true);
       else if (startswith("image", buffer))
         processImage(buffer, false);
       else if (startswith("volimage", buffer))
          processImage3D(buffer);
       else if (startswith("ssioutput", buffer))
          processESSI3D(buffer);
       else if (startswith("essioutput", buffer))
          processESSI3D(buffer);
       else if (startswith("boundary_conditions", buffer))
         processBoundaryConditions(buffer);
       //       else if (startswith("supergrid", buffer))
       //         processSupergrid(buffer);
       // else if (startswith("prefilter", buffer))
       // 	 processPrefilter(buffer);
       else if( startswith("developer", buffer ) )
          processDeveloper(buffer);
       else if( startswith("geodynbc", buffer ) )
          processGeodynbc(buffer);
       else if( startswith("randomize", buffer ) )
       {
	  //          processRandomize(buffer);
	  if( m_myRank == 0 )
	     cout << "randomize command is no longer supported. Use `randomblock' instead" <<endl;
       }
       else if( startswith("randomblock", buffer ) )
          processRandomBlock(buffer);
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

// tmp:
  // if (m_myRank == 0)
  // {
  //   cout << "INFO: m_mesh_refinements=" << m_mesh_refinements << " m_use_attenuation=" << m_use_attenuation << " mOrder=" << mOrder << endl;
  // }
  
  if (mVerbose >=3 && proc_zero())
    cout << "********Done reading the input file*********" << endl;

// wait until all processes have read the input file
  MPI_Barrier(MPI_COMM_WORLD);

  if (proc_zero())
    if (a_GlobalTimeSeries.size() > 0 && a_GlobalTimeSeries[0].size() > 0) 
      cout << "Read station input, took " << a_GlobalTimeSeries[0][0]->getReadTime() << "seconds." << endl;

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


//-----------------------------------------------------------------------
void EW::processGrid(char* buffer)
{
  float_sw4 x = 0.0;
  float_sw4 y = 0.0;
  float_sw4 z = 0.0;
  int nx=0, ny=0, nz=0;
  float_sw4 h = 0.0;

  //-----------------------------------------------------------------
  // default geographical coordinates will be the 
  // nevada test site (see:  en.wikipedia.org/wiki/Nevada_Test_Site
  //-----------------------------------------------------------------
  double lat, lon;
  bool latSet = false, lonSet = false, lon_p_set=false, lat_p_set=false, datum_set=false;
  bool ellps_set=false, proj_set=false;
  bool use_geoprojection=false;
  
  stringstream proj0;

// hard code units to be in meters
  proj0 << "+units=m";

// default azimuth
  mGeoAz=0;
  
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

  if (proc_zero_evzero() )
    cout << endl << "* Processing the grid command..." << endl;

  // Assume presence of mesh refinements has already been checked.
  // Assume supergrid command has already been processed.
  if( m_mesh_refinements || m_sg_damping_order == 6 )
  {
     m_ghost_points = 3;
     m_ppadding = 3;
  }

  // if (m_myRank == 0)
  //    cout << endl << "* number of ghost points = " << m_ghost_points << endl;

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
//                        1234567890
     else if (startswith("mlon=", token))
     {
        token += 5;
        mMetersPerLongitude = atof(token);
        CHECK_INPUT(mMetersPerLongitude > 0.0,
                err << "mMetersPerLongitude must be greater than 0, not " 
                << mMetersPerLongitude);
        mConstMetersPerLongitude = true;
     }
//                        1234567890
     else if (startswith("mlat=", token))
     {
        token += 5;
        mMetersPerDegree = atof(token);
        CHECK_INPUT(mMetersPerDegree > 0.0,
                err << "mMetersPerDegree must be greater than 0, not " 
                << mMetersPerDegree);
     }
//                        1234567890123456  
     else if (startswith("extrapolate=", token))
     {
        token += 12;
        int extrapolate = atoi(token);
        CHECK_INPUT(extrapolate >= 0 && extrapolate <= 5,
                err << "extrapolate must be an integer between 0 and 5, not " 
                << extrapolate);
	mMaterialExtrapolate = extrapolate;
     }
     //     else if( startswith("ghostpts=",token))
     //     {
     //	token += 9;
     //        int ghost = atoi(token);
     //	CHECK_INPUT( ghost == 2 || ghost == 3, err << "Number of ghost points must be 2 or 3, not " << ghost );

     //	if( m_mesh_refinements && ghost == 2 )
     //	   CHECK_INPUT( false, err << "Number of ghost points must be 3 when using mesh refinement  ");
     //	m_ghost_points = ghost;
     //        m_ppadding = ghost;
     //     }
//                        123456789
     else if( startswith("proj=",token))
     {
        token +=5;
// accumulate new style string
        proj0 << " +proj=" << token;
	use_geoprojection = true;
        proj_set=true;
     }
//                        123456789
     else if( startswith("ellps=",token))
     {
        token +=6;
// accumulate new style string
        proj0 << " +ellps=" << token;
	use_geoprojection = true;
        ellps_set=true;
     }
//                        123456789
     else if( startswith("datum=",token))
     {
        token +=6;
        proj0 << " +datum=" << token;
        datum_set=true;
	use_geoprojection = true;
     }
//                        123456789
     else if( startswith("lon_p=",token))
     {
        token +=6;
        proj0 << " +lon_0=" << atof(token);
	use_geoprojection = true;
        lon_p_set=true;
     }
//                        123456789
     else if( startswith("lat_p=",token))
     {
        token +=6;
        proj0 << " +lat_0=" << atof(token);
	use_geoprojection = true;
        lat_p_set=true;
     }
//                        123456789
     else if( startswith("scale=",token))
     {
        token +=6;
        proj0 << " +scale=" << atof(token);
	use_geoprojection = true;
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
  float_sw4 xprime, yprime, zprime;
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

// check syntax for lat & lon
  if (!(latSet && lonSet) && (latSet || lonSet))
  {
    stringstream msg;
    if (m_myRank == 0)
    {
      msg << " \n* Improper grid location specification, must specify both lat and lon variables " << endl
	  << " * Missing... ";
      if (!latSet)
	msg << " lat=value ";
      if (!lonSet)
	msg << " lon=value ";
    }
    CHECK_INPUT(0, msg.str());
  }

  if (latSet && lonSet)
  {
     mLatOrigin = lat;
     mLonOrigin = lon;
  }
  else
  {
// Default is NTS
     mLatOrigin = 37.0;
     mLonOrigin = -118.0;
  }

// default arguments for proj4 projection
  if (use_geoprojection)
  {
     if (!proj_set)
     {
// Default projection: Universal Transverse Mercator (UTM)
        proj0 << " +proj=utm";
     }

     if (!ellps_set && !datum_set)
     {
// default ellipse
        proj0 << " +ellps=WGS84";
     }
     
// if lon_p not given, use lon
     if (!lon_p_set)
     {
        proj0 << " +lon_0=" << mLonOrigin;
     }

// if lat_p not given, use lat
     if (!lat_p_set)
     {
        proj0 << " +lat_0=" << mLatOrigin;
     }
  }

  float_sw4 cubelen, zcubelen, hcube;
  if( m_geodynbc_found )
  {
// Set SW4 grid spacing based on Geodyn cube data

     float_sw4 origin[3]={0,0,0};
     double ibclat, ibclon, ibcaz;

      bool found_latlon;
      int adjust;
      geodynbcGetSizes( m_geodynbc_filename, origin, cubelen, zcubelen, hcube, found_latlon,
 		       ibclat, ibclon, ibcaz, adjust );
// Use approximate h
/*
      if( h == 0.0 )
      {
 	if( nx > 0 )
 	   h = x/(nx-1);
 	else if( nz > 0 )
 	   h = z/(nz-1);
 	else
 	   h = y/(ny-1);
      }
*/

      // rounding of cube position to two decimals (prec=100), three (prec=1000) etc..
      float_sw4 prec = 100;

      if( found_latlon )
      {
         CHECK_INPUT( fabs(ibcaz - mGeoAz) < 1e-5, "Error: Az in Geodyn file, "
 		 << ibcaz << " is different from Az in WPP, " << mGeoAz );
	   
 	// lat-lon corner of cube given
 	if( adjust == 1 || origin[2] == 0 )
 	{
 	   // h based on cube length only, adjust z-position of cube
 	   /*
 	   int nc = static_cast<int>(round(cubelen)/h);
 	   h = cubelen/nc;
           */
 	   origin[2] -= hcube*( origin[2]/hcube-round(origin[2]/hcube) );
 	}
        else
 	{
 	   // h based on cube length and z-position of cube
 	   /*
	   int a = static_cast<int>(round(origin[2]*prec));
 	   int b = static_cast<int>(round((origin[2]+zcubelen)*prec));
 	   // 
           int d  = gcd(a,b);
 	   int n1 = a/d;
 	   int k  = static_cast<int>(round(origin[2]/(n1*h)));
	   h = origin[2]/(k*n1);
           */
 	}
 	// Geographic origin adjustment:
	double gridLat = mLatOrigin;
	double gridLon = mLonOrigin;
	double metersPerDegree = mMetersPerDegree;
	double deg2rad = M_PI/180;
	double phi = mGeoAz*deg2rad;
 	float_sw4 x = metersPerDegree*( cos(phi)*(ibclat-gridLat) + cos(ibclat*deg2rad)*(ibclon-gridLon)*sin(phi));
 	float_sw4 y = metersPerDegree*(-sin(phi)*(ibclat-gridLat) + cos(ibclat*deg2rad)*(ibclon-gridLon)*cos(phi));
	x -= hcube*(x/hcube-round(x/hcube));
	y -= hcube*(y/hcube-round(y/hcube));
	gridLat = ibclat - (x*cos(phi) - y*sin(phi))/metersPerDegree;
 	gridLon = ibclon - (x*sin(phi) + y*cos(phi))/(metersPerDegree*cos(ibclat*deg2rad));
	mLatOrigin = gridLat;
	mLonOrigin = gridLon;
        origin[0] = x;
 	origin[1] = y;
      } // end if latlon given
      else
      {
 	// lat-lon corner of cube not given, interpret origin realtive (0,0,0)
         if( m_geodynbc_center )
	 {
 	   // Center cube in the middle of the domain (in x,y), discarding input origin.
	    float_sw4 xlen = x;
	    float_sw4 ylen = y;
	    if( xlen == 0 )
	       xlen = hcube*(nx-1);
	    if( ylen == 0 )
	       ylen = hcube*(ny-1);
	    origin[0] = 0.5*(xlen-cubelen);
	    origin[1] = 0.5*(ylen-cubelen);
	 }
	 if( adjust == 1 )
	 {
	   // h based on cube length only, adjust cube position
	   /*
 	   int nc = static_cast<int>(round(cubelen/h));
 	   h = cubelen/nc;
           */
	   //	   cout << "nc= " << nc << " cubelen= " << cubelen << " origin before " <<
	   //	      origin[0] << " " << origin[1] << " " << origin[2] << endl;
 	   origin[0] -= hcube*( origin[0]/hcube-round(origin[0]/hcube) );
 	   origin[1] -= hcube*( origin[1]/hcube-round(origin[1]/hcube) );
 	   origin[2] -= hcube*( origin[2]/hcube-round(origin[2]/hcube) );
	   //	   cout << " origin after " <<
	   //	      origin[0] << " " << origin[1] << " " << origin[2] << endl;

	 }
	 else
	 {
 	   // h based on cube length and cube position, might be very restrictive
	    CHECK_INPUT( false, "Error: cube position without lat/long position must be adjustable");
	 }
      } // end if latlon not given
      
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
      m_ibc_origin[0] = origin[0];
      m_ibc_origin[1] = origin[1];
      m_ibc_origin[2] = origin[2];
      //     cout << "Cube origin " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
      //     cout << "Cube length " << cubelen << endl;
      //     cout << "nx,ny,nz " << nxprime << " " << nyprime << " " << nzprime << endl;     
  } // end if m_geodynbc_found
  else
  {

     if (!m_doubly_periodic)
     {
	if (nx > 0 && h == 0.0)
	{
    // we set the number grid points in the x direction
    // so we'll compute the grid spacing from that.
	   h = x / (nx-1);
	   if (proc_zero_evzero())
	      cout << "* Setting h to " << h << " from  x/(nx-1) (x=" << x << ", nx=" << nx << ")" << endl;
      
	   nxprime = nx;
	   nzprime = computeEndGridPoint(z, h);
	   nyprime = computeEndGridPoint(y, h);
	}
	else if (ny > 0 && h == 0.0)
	{
    // set hte number of grid points from y direction and ny
	   h = y/(ny-1);
	   if (proc_zero_evzero())
	      cout << "* Setting h to " << h << " from  y/(ny-1) (y=" << y << ", ny=" << ny << ")" << endl;
	   nyprime = ny;
	   nxprime = computeEndGridPoint(x, h);
	   nzprime = computeEndGridPoint(z, h);
	}
	else if (nz > 0 && h == 0.0)
	{
    // set the number of grid points from z direction and nz
	   h = z/(nz-1);
	   if (proc_zero_evzero())
	      cout << "* Setting h to " << h << " from  z/(nz-1) (z=" << z << ", nz=" << nz << ")" << endl;
	   nzprime = nz;
	   nxprime = computeEndGridPoint(x, h);
	   nyprime = computeEndGridPoint(y, h);
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
  }
  if (!m_doubly_periodic)
  {
    if (proc_zero_evzero() && mVerbose >=3)
      printf("**** Setting up the grid for a non-periodic problem\n");
    
    if (nxprime != nx && proc_zero_evzero())
      cout << "* Setting nx to " << nxprime << " to be consistent with h=" << h << endl;
    if (nyprime != ny && proc_zero_evzero())
      cout << "* Setting ny to " << nyprime << " to be consistent with h=" << h << endl;
    if (nzprime != nz && proc_zero_evzero())
      cout << "* Setting nz to " << nzprime << " to be consistent with h=" << h << endl;

    // -------------------------------------------------------------
    // Now we adjust the geometry bounds based on the actual 
    // number of grid points used in each dimension.
    // -------------------------------------------------------------
    xprime = (nxprime-1)*h;
    zprime = (nzprime-1)*h;
    yprime = (nyprime-1)*h;
  
    float_sw4 eps = 1.e-9*sqrt(SQR(xprime)+SQR(yprime)+SQR(zprime));
    if( sizeof(float_sw4)==4)
       eps=eps*1e4;
  
    if (fabs(xprime-x) > eps && proc_zero_evzero())
      cout << "* Changing x from " << x << " to " << xprime << " to be consistent with h=" << h << endl;
    if (fabs(zprime-z) > eps && proc_zero_evzero())
      cout << "* Changing z from " << z << " to " << zprime << " to be consistent with h=" << h << endl;
    if (fabs(yprime-y) > eps && proc_zero_evzero())
      cout << "* Changing y from " << y << " to " << yprime << " to be consistent with h=" << h << endl;
  }
  else // special treatment of the doubly periodic case
  {
    if (proc_zero_evzero() && mVerbose >=3)
      printf("**** Setting up the grid for a PERIODIC problem\n");

// for the doubly periodic case, we only support the following style:
// grid x=... y=... z=... nx=...    
    CHECK_INPUT(nx > 0 && x>0. && y>0. && z>0., 
		"Period case: Must specify grid using x, y, z, nx");

    // we set the number grid points in the x direction
    // so we'll compute the grid spacing from that.
    h = x / nx;
    if (proc_zero_evzero())
      cout << "* Setting h to " << h << " from  x/nx (x=" << x << ", nx=" << nx << ")" << endl;
      
    nxprime = nx;
    nyprime = (int) (y/h + 0.5);
    nzprime = computeEndGridPoint(z, h); // non-periodic in z

    // -------------------------------------------------------------
    // Now we adjust the geometry bounds based on the actual 
    // number of grid points used in each dimension.
    // -------------------------------------------------------------
    xprime = nxprime*h;
    yprime = nyprime*h;
    zprime = (nzprime-1)*h; // non-periodic in z
  
    float_sw4 eps = 1.e-9*sqrt(SQR(xprime)+SQR(yprime)+SQR(zprime));
    if( sizeof(float_sw4)==4)
       eps=eps*1e4;
  
    if (fabs(xprime-x) > eps && proc_zero_evzero())
      cout << "* Changing x from " << x << " to " << xprime << " to be consistent with h=" << h << endl;
    if (fabs(yprime-y) > eps && proc_zero_evzero())
      cout << "* Changing y from " << y << " to " << yprime << " to be consistent with h=" << h << endl;
    if (fabs(zprime-z) > eps && proc_zero_evzero())
      cout << "* Changing z from " << z << " to " << zprime << " to be consistent with h=" << h << endl;
  }
  
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

  m_nx_base = nxprime;
  m_ny_base = nyprime;
  m_nz_base = nzprime;
  m_h_base = h;
  m_global_xmax = xprime;
  m_global_ymax = yprime;
  m_global_zmax = zprime;

#ifndef ENABLE_PROJ4
  CHECK_INPUT( !use_geoprojection, "ERROR: need to configure SW4 with proj=yes to use projections "
               "from the Proj4 library");
#endif
  if( use_geoprojection )
  {
// tmp
//     cout << "New proj4 string: '" << proj0.str() << "'" << endl;
     
     m_geoproj = new GeographicProjection( mLonOrigin, mLatOrigin, proj0.str(), mGeoAz );
  }
  else
     m_geoproj = static_cast<GeographicProjection*>(0);

}

//----------------------------------------------------------
void EW::cleanUpRefinementLevels()
{


// Add a top zMin level
// Here zMin = m_topo_zmax if m_topography_exists, otherwise zMin = 0;
   float_sw4 zMin, topo_zmax=0;
  
// NOW: allowing refinements in the curvilinear portion of the grid
   if (m_topography_exists)
   {
      topo_zmax = m_gridGenerator->get_topo_zmax();
      CHECK_INPUT(topo_zmax < m_global_zmax-m_h_base,"The topography is extending too deep into the ground and there is no space for the Cartesian grid.");

      m_curviRefLev.push_back(0.0); // for the curvilinear refinements
      m_refinementBoundaries.push_back(topo_zmax); // for the Cartesian refinements
      zMin = topo_zmax; 
   }
   else
   {
      m_refinementBoundaries.push_back(0.0); // flat free surface boundary
      zMin = 0.;
   }

// need to sort m_refinementBoundaries in decreasing order
   int nRef = m_refinementBoundaries.size();
   float_sw4 *zValues = new float_sw4[nRef];
   int q;

   for (q=0; q<nRef; q++)
      zValues[q] = m_refinementBoundaries[q];
   sort(zValues, zValues+nRef);
// reverse the ordering to get decreasing order
   for (q=0; q<nRef; q++)
      m_refinementBoundaries[q] = zValues[nRef-q-1];

// cleanup
  delete [] zValues;

  vector<float_sw4>::iterator it;
//  cout << "Removing items outside the range zMin = " << zMin << " < z < " << " zMax=" << m_zmax << "..." << endl;
  for (it=m_refinementBoundaries.begin(); it!=m_refinementBoundaries.end(); it++)
  {
    if (*it < zMin || *it >= m_global_zmax)
    {
// if 0 < zLev < zMin: add this as a curvilinear refinemtn level
       if (m_topography_exists)
       {
          float_sw4 zLev = *it;
          if (zLev > 0 && zLev < topo_zmax)
             m_curviRefLev.push_back(zLev);
       }
// remove this entry from the vector
//      cout << "Removing out-of-range refinement level="<< *it << endl;
      it = m_refinementBoundaries.erase(it); // returns next element
// need to back up one step to undo the it++ that always happens at the end of the for loop
      it--;
    }
  }

  // sort m_curviRefLev in decreasing order
  nRef = m_curviRefLev.size();
  zValues = new float_sw4[nRef];
  
  for (q=0; q<nRef; q++)
     zValues[q] = m_curviRefLev[q];
  sort(zValues, zValues+nRef);
// reverse the ordering to get decreasing order
   for (q=0; q<nRef; q++)
      m_curviRefLev[q] = zValues[nRef-q-1];

   delete[] zValues;  
// need to remove any duplicate entries in the m_refinementBoundaries array
// tmp
//  cout << "Removing duplicate items..."<< endl;
  float_sw4 z0 = m_refinementBoundaries[0]; // first item
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
  if (mVerbose >= 1 && proc_zero_evzero())
  {
     cout << "cleanupRefinementLevels: topo_zmax = " << topo_zmax << endl;

    cout << " Cartesian refinement levels (z=):" << endl;
    for (it=m_refinementBoundaries.begin(); it!=m_refinementBoundaries.end(); it++)
      cout<< *it << endl;

    if (m_topography_exists)
    {
       cout << " Curvilinear refinement levels (z=):" << endl;
       for (it=m_curviRefLev.begin(); it!=m_curviRefLev.end(); it++)
          cout<< *it << endl;
    }
    
  }
}


//-----------------------------------------------------------------------
void EW::processRefinement(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("refinement", token) == 0, 
	      "ERROR: not a refinement line...: " << token);
   token = strtok(NULL, " \t");
   string err = "Refinement error ";

   while (token != NULL)
   {
     // while there are tokens in the string still
     if (startswith("#", token) || startswith(" ", buffer))
       // Ignore commented lines and lines with just a space.
       break;
     else if( startswith("zmax=", token) )
     {
       token += 5; // skip zmax=
       float_sw4 z1 = atof(token);
       m_refinementBoundaries.push_back(z1);
 //       if (m_myRank==0)
 // 	cout <<"Adding refinement boundary at z=" << z1 << endl;
     }
     else
     {
       badOption("refinement", token);
     }
     token = strtok(NULL, " \t");
   }
}

//-----------------------------------------------------------------------
void EW::processAttenuation(char* buffer)
{
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("attenuation", token) == 0, 
	      "ERROR: not a attenuation line...: " << token);
  token = strtok(NULL, " \t");

   string err = "Attenuation error ";
   int nmech=3;
//   int nmech=-1;
   float_sw4 velofreq=1;
   bool foundppw = false, foundfreq=false;

// Default is max frequency 2 Hz, 
   m_att_ppw = -1;
   m_att_max_frequency = 2.0;
  
   while (token != NULL)
   {
    // while there are tokens in the string still
     if (startswith("#", token) || startswith(" ", buffer))
       // Ignore commented lines and lines with just a space.
       break;
 //                       123456
     else if( startswith("nmech=", token) )
     {
       token += 6; // skip nmech=
       nmech = atoi(token);
       CHECK_INPUT(nmech >= 0 && nmech <= 8, "ERROR: Number of attenuation mechanisms must be >= 0 and <= 8, not " << nmech);
     }
 //                       1234567890123
     else if( startswith("phasefreq=", token) )
     {
       token += 10; // skip phasefreq=
       velofreq = atof(token);
       CHECK_INPUT(velofreq >= 0 && velofreq <= 1000, "ERROR: Velocity frequency must be >= 0 and <= 1000 [Hz], not " << velofreq);
     }
     else if( startswith("maxfreq=",token) )
     {
        token += 8;
        m_att_ppw = -1;
        m_att_max_frequency = atof( token );
        CHECK_INPUT( !foundfreq, "ERROR: can not give both centerfreq and maxfreq");
        CHECK_INPUT(m_att_max_frequency >= 0,"ERROR: maximum frequency must be >= 0, not " << m_att_max_frequency);
        foundfreq = true;
     }
     // else if( startswith("centerfreq=",token) )
     // {
     //    token += 11;
     //    m_att_ppw = -1;
     //    m_att_max_frequency = atof( token );
     //    CHECK_INPUT( !foundfreq, "ERROR: can not give both centerfreq and maxfreq");
     //    CHECK_INPUT(m_att_max_frequency >= 0,"ERROR: maximum frequency must be >= 0, not " << m_att_max_frequency);
     //    foundfreq = true;
     // }
     else if( startswith("minppw=",token) )
     {
        token += 7;
        m_att_ppw = atof( token ); // AP: changed from atoi
        foundppw = true;
        CHECK_INPUT(m_att_ppw >= 0, "ERROR: minimum ppw must be >= 0, not " << m_att_ppw);
     }
     else if( startswith("qmultiplier=",token) )
     {
        token += 12;
        m_qmultiplier = atof( token ); //
        CHECK_INPUT(m_qmultiplier > 0, "ERROR: qmultiplier must be positive, not " << m_qmultiplier);	
     }
     else
     {
       badOption("attenuation", token);
     }
     token = strtok(NULL, " \t");
   }
   if( foundppw && foundfreq )
   {
      if (m_myRank == 0)
 	cout << "ERROR: Can not give both minppw and maxfreq for attenuation " << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   m_number_mechanisms = nmech;
   m_velo_omega = velofreq*2*M_PI;
   m_use_attenuation=true;
   m_att_use_max_frequency = (m_att_ppw <= 0);
  
 // tmp
   //   if (m_myRank==0)
   //     printf("* Processing the attenuation command: m_nmech=%i, m_velo_omega=%e\n", m_nmech, m_velo_omega);
}

//-----------------------------------------------------------------------
void EW::processTopography(char* buffer)
{
   //
   // Note, m_topoFileName, m_topoExtFileName, m_maxIter, m_EFileResolution, m_QueryTyp could
   // have been declared local variables in EW::parseInputFile, and transfered as
   // procedure parameters to smoothTopography and getEfileInfo
   //
    char* token = strtok(buffer, " \t");
    CHECK_INPUT(strcmp("topography", token) == 0, 
 	    "ERROR: not a topography line...: " << token);
    string topoFile="surf.tp", style, fileName;
    bool needFileName=false, gotFileName=false;
    
    float_sw4 zetaBreak=0.95, topo_zmax=0;
    float_sw4 GaussianAmp=0.05, GaussianLx=0.15, GaussianLy=0.15, GaussianXc=0.5, GaussianYc=0.5;
    int grid_interpolation_order = 3;
    bool use_analytical_metric = false, topo_zmax_given=false;
    bool always_new = false;


    token = strtok(NULL, " \t");

    while (token != NULL)
    {
      // while there are still tokens in the string 
       if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
	  break;
       if (startswith("zmax=", token))
       {
	  token += 5; // skip logfile=
	  topo_zmax = atof(token);
       }
       else if (startswith("order=", token))
       {
	  token += 6; // skip logfile=
	  grid_interpolation_order = atoi(token);
	  if (grid_interpolation_order < 2 || grid_interpolation_order > 7)
	  {
	     if (m_myRank == 0)
		cout << "order needs to be 2,3,4,5,6,or 7 not: " << grid_interpolation_order << endl;
	     MPI_Abort(MPI_COMM_WORLD, 1);
	  }
       }
       else if( startswith("zetabreak=", token) ) // developer option: not documented in user's guide
       {
	  token += 10;
	  zetaBreak = atof(token);
	  CHECK_INPUT( zetaBreak > 0 && zetaBreak <= 1, "Error: zetabreak must be in [0,1], not " << zetaBreak);
       }
       else if (startswith("smooth=", token))
       {
	  token += 7; // skip smooth=
	  m_maxIter = atoi(token);
	  if (m_maxIter < 0 || m_maxIter > 1000)
	  {
	     if (m_myRank == 0)
		cout << "Number of smoothing iterations needs to be >=0 and <=1000, not: "<< m_maxIter << endl;
	     MPI_Abort(MPI_COMM_WORLD, 1);
	  }
       }
       else if( startswith("input=", token ) )
       {
	  token += 6;
	  style = token;
// new keyword: geographic, but keeping grid for backwards compatibility with WPP
	  if (strcmp("grid", token) == 0 || strcmp("geographic", token) == 0)
	  {
	     m_topoInputStyle=GridFile;
	     m_topography_exists=true;
	     needFileName=true;
	  }
	  else if (strcmp("cartesian", token) == 0)
	  {
	     m_topoInputStyle=CartesianGrid;
	     m_topography_exists=true;
	     needFileName=true;
	  }
	  else if (strcmp("efile", token) == 0)
	  {
	     m_topoInputStyle=Efile;
	     m_topography_exists=true;
	     needFileName=true; // we require the file name to be given on the topography command line
	  }
	  else if (strcmp("rfile", token) == 0)
	  {
	     m_topoInputStyle=Rfile;
	     m_topography_exists=true;
	     needFileName=true; // we require the file name to be given on the topography command line
	  }
	  else if (strcmp("sfile", token) == 0)
	  {
	     m_topoInputStyle=Sfile;
	     m_topography_exists=true;
	     needFileName=true; // we require the file name to be given on the topography command line
	  }
	  else if (strcmp("image", token) == 0)
	  {
	     m_topoInputStyle=TopoImage;
	     m_topography_exists=true;
	     needFileName=true; // we require the file name to be given on the topography command line
	  }
	  else if (strcmp("gaussian", token) == 0)
	  {
	     m_topoInputStyle=GaussianHill;
	     m_topography_exists=true;
	  }
	  else
	  {
	     badOption("topography> input", token);
	  }
       }
       else if( startswith("file=", token ) )
       {
	  token += 5;
	  m_topoFileName = token;
	  gotFileName=true;
	//        if (m_myRank==0)
	// 	 cout << "read topo file name=" << m_topoFileName <<endl;
       }
#ifdef ENABLE_ETREE
       else if( startswith("etree=", token ) )
       {
	  token += 6;
	  m_topoFileName = token;
	  m_topoInputStyle=Efile;
	  m_topography_exists=true;
	  gotFileName=true;
       }
      else if (startswith("xetree=", token))
      {
	 token += 7; // skip xetree=
	 m_topoExtFileName = token;
      }
#endif
//                        12345678901
       else if( startswith("resolution=", token ) )
       {
	  token += 11;
	  m_EFileResolution = atof(token);
	  CHECK_INPUT(m_EFileResolution>0.,"Resolution must be positive, not " << m_EFileResolution);
       }
//                        123456789012
       else if( startswith("gaussianAmp=", token ) )
       {
	  token += 12;
	  GaussianAmp = atof(token);
       }
//                        123456789012
       else if( startswith("gaussianXc=", token ) )
       {
	  token += 11;
	  GaussianXc = atof(token);
       }
//                        123456789012
       else if( startswith("gaussianYc=", token ) )
       {
	  token += 11;
	  GaussianYc = atof(token);
       }
// //                        123456789012
       else if( startswith("gaussianLx=", token ) )
       {
	  token += 11;
	  GaussianLx = atof(token);
       }
       else if( startswith("gaussianLy=", token ) )
       {
	  token += 11;
	  GaussianLy = atof(token);
       }
       else if( startswith("analyticalMetric=", token ) )
       {
	  token += 17;
	  use_analytical_metric = strcmp(token,"1")==0 ||
	     strcmp(token,"true")==0 || strcmp(token,"yes")==0;
       }
       else if (startswith("gridgenerator=", token) )
       {
          token += 14;
	  always_new =  strcmp(token,"new")==0 || strcmp(token,"NEW")==0;
       }
       else
       {
	  badOption("topography", token);
       }
       token = strtok(NULL, " \t");
    }
    if (needFileName)
       CHECK_INPUT(gotFileName, 
		   "ERROR: no topography file name specified...: " << token);

    if( m_topoInputStyle != GaussianHill && use_analytical_metric )
    {
       use_analytical_metric = false;
       if( m_myRank == 0 )
	  cout << "Analytical metric only defined for Gaussian Hill topography" <<
	     " topography analyticalMetric option will be ignored " << endl;
    }

    if( m_topoInputStyle == GaussianHill )
       m_gridGenerator = new GridGeneratorGaussianHill( topo_zmax, always_new, use_analytical_metric,
                                                        grid_interpolation_order, zetaBreak, GaussianAmp,
                                                        GaussianXc, GaussianYc, GaussianLx, GaussianLy );
    else
       m_gridGenerator = new GridGeneratorGeneral( topo_zmax, always_new,
                                                   grid_interpolation_order, zetaBreak );
}

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

  float_sw4 omega = 1.;
  float_sw4 momega= 1.;
  float_sw4 phase = 0;
  float_sw4 mphase= 0.4;
  float_sw4 c = 1.3;
  float_sw4 amprho = 1;
  float_sw4 ampmu  = 1;
  float_sw4 amplambda = 1;
  float_sw4 omstrx=1.1, omstry=0.8, omstrz=0.9;

  int sgstretch = 0;
  int frsurfu=1, frsurfl=0;

  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("twilight", token) == 0, "ERROR: not a twilight line...: " << token);
  token = strtok(NULL, " \t");

  string err = "Twilight command syntax error: ";

  while (token != NULL)
    {
      // while there are tokens in the string still
       if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
//                     123456789
       if( startswith("errorlog=",token) )
       {
          token += 9;
	  bool errorlog = (atoi(token)==1);
	  if( errorlog )
	     switch_on_error_log();
       }
       else if( startswith("sgstretching=",token) )
       {
          token += 13;
          if( strcmp(token,"1") == 0 || strcmp(token,"yes") == 0 || strcmp(token,"true") == 0 )
	     sgstretch = 1;
	  else
	     sgstretch = 0;
       }
       else if( startswith("freeupper=",token) )
       {
          token += 10;
          if( strcmp(token,"1") == 0 || strcmp(token,"yes") == 0 || strcmp(token,"true") == 0 )
	     frsurfu = 1;
	  else
	     frsurfu = 0;
       }
       else if( startswith("freelower=",token) )
       {
          token += 10;
          if( strcmp(token,"1") == 0 || strcmp(token,"yes") == 0 || strcmp(token,"true") == 0 )
	     frsurfl = 1;
	  else
	     frsurfl = 0;
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
       else if( startswith("omstrx=",token) )
       {
          token += 7;
	  omstrx = atof(token);
       }
       else if( startswith("omstry=",token) )
       {
          token += 7;
	  omstry = atof(token);
       }
       else if( startswith("omstrz=",token) )
       {
          token += 7;
	  omstrz = atof(token);
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

  // Default to Dirichlet conditions on all sides
  boundaryConditionType bct[6]={bDirichlet, bDirichlet, bDirichlet, bDirichlet, bDirichlet, bDirichlet};

  // Free surface conditions upper side,
  if( frsurfu == 1 )
     bct[4] = bStressFree;

  // Free surface conditions lower side,
  if( frsurfl == 1 )
     bct[5] = bStressFree;

  // Use supergrid stretching
  if( sgstretch == 1 )
  {
     for( int side=0 ; side < 4 ; side++ )
	if( bct[side] == bDirichlet )
	   bct[side] = bSuperGrid;

     for( int side=4 ; side < 5 ; side++ )
	if( bct[side] == bDirichlet && !topographyExists() )
	   bct[side] = bSuperGrid;
     
     if( bct[0] == bSuperGrid || bct[1] == bSuperGrid )
     {
	for( int g=0 ; g < mNumberOfGrids ; g++ )
	   m_supergrid_taper_x[g].set_twilight(omstrx);
     }
     if( bct[2] == bSuperGrid || bct[3] == bSuperGrid )
     {
	for( int g=0 ; g < mNumberOfGrids ; g++ )
	   m_supergrid_taper_y[g].set_twilight(omstry);
     }
     CHECK_INPUT( (bct[4] == bSuperGrid && bct[5] == bSuperGrid) || (bct[4] == bStressFree && bct[5] == bStressFree) || (bct[4]==bDirichlet || bct[5] == bDirichlet),
	   "Error: Twilight testing with supergrid stretching must have the same b.c. (stress free or supergrid) on the z=low and z=high boundaries" );
     if( bct[4] == bSuperGrid && bct[5] == bSuperGrid )
	CHECK_INPUT( !topographyExists(), "Error: Twilight testing, supergrid stretching can not be used in the z-direction when topography is present");
	
     for( int g=0 ; g < mNumberOfGrids ; g++ )
	m_supergrid_taper_z[g].set_twilight(0.0);
     
     if( bct[4] == bSuperGrid && bct[5] == bSuperGrid )
     {
	m_supergrid_taper_z[mNumberOfGrids-1].set_twilight(omstrz);
	m_supergrid_taper_z[0].set_twilight(omstrz);
     }
// set the damping coefficient to zero
     set_sg_damping(0.0);
     set_sg_thickness(1); // just to keep the routine assign_supergrid_damping_arrays() happy
     
  } // end if sgstretch == 1
  set_global_bcs(bct);
}

//-----------------------------------------------------------------------
void EW::processDeveloper(char* buffer)
{
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
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("developer", token) == 0, "ERROR: not a developer line...: " << token);
   token = strtok(NULL, " \t");
   while (token != NULL)
   {
     // while there are tokens in the string still
     if (startswith("#", token) || startswith(" ", buffer))
       // Ignore commented lines and lines with just a space.
        break;
     if (startswith("opttest=", token))
     {
        token += 8; // skip opttest=
        if( strcmp(token,"source")==0 )
	   m_opttest = 1;
        else if( strcmp(token,"gradient")== 0 )
	   m_opttest = 2;
	else if( strcmp(token,"hessian") == 0 )
	   m_opttest = 3;
	else if( strcmp(token,"func1d") == 0 )
	   m_opttest = 4;
	else if( strcmp(token,"funcsurf") == 0 )
	   m_opttest = 5;
        else
	   CHECK_INPUT( false, "ERROR: opttest=" << token << " not understood");
        if( !m_inverse_problem )
	   CHECK_INPUT( false, "WARNING: developer opttest option does not apply to forward solver");
     }
     else if( startswith("cfl=",token) )
     {
	token += 4;
	float_sw4 cfl = atof(token);
	CHECK_INPUT( cfl > 0, "Error negative CFL number");
	set_cflnumber( cfl );
     }
     else if( startswith("time_order=",token) )
     {
	token += 11;
	int newOrder = atoi(token);
	CHECK_INPUT( newOrder == 2 || newOrder == 4, "Error unknown time-order");
	mOrder = newOrder;
     }
     else if( startswith("perturb=",token) )
     {
        token += 8;
	m_perturb = atof(token);
     }
     else if( startswith("peri=",token) )
     {
        token += 5;
	m_iperturb = atoi(token);
     }
     else if( startswith("perj=",token) )
     {
        token += 5;
	m_jperturb = atoi(token);
     }
     else if( startswith("perk=",token) )
     {
        token += 5;
	m_kperturb = atof(token);
     }
     else if( startswith("pervar=",token) )
     {
        token += 7;
        if( strcmp(token,"mu")==0 )
	   m_pervar = 1;
	else if( strcmp(token,"lambda")==0 )
	   m_pervar = 2;
	else if( strcmp(token,"rho")==0 )
	   m_pervar = 0;
	else
	   CHECK_INPUT(false," pervar must be , mu, lambda, or rho, not " << token << endl);
     }
     else if( startswith("checkfornan=",token) )
     {
	token += 12;
	m_checkfornan = strcmp(token,"1")==0 || strcmp(token,"on")==0 || strcmp(token,"yes")==0;
     }

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
     else if( startswith("reporttiming=",token) )
     {
	token += 13;
	m_output_detailed_timing = strcmp(token,"1")==0 || strcmp(token,"on")==0
	   || strcmp(token,"yes")==0;
     }
//     else if (startswith("interpolation=", token))
//     {
//       token += 14;
//       cons = strcmp(token,"conservative") == 0;
//     }
     else if (startswith("ctol=", token))
     {
       token += 5;
       m_citol = atof(token);
     }
     else if (startswith("cmaxit=", token))
     {
       token += 7;
       m_cimaxiter = atoi(token);
     }
     else if (startswith("crelax=", token))
     {
       token += 7;
       m_cirelfact = atof(token);
     }
     //     else if (startswith("ckernels=", token))
     //     {
     //       token += 9;
     //       m_croutines = atoi(token)==1;
     //       Sarray::m_corder = m_croutines;
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
     else
     {
       badOption("developer", token);
     }
     token = strtok(NULL, " \t");
   }
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
}

//-----------------------------------------------------------------------
void EW::processTestPointSource(char* buffer)
{
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("testpointsource", token) == 0, "ERROR: not a testpointsource line...: " << token);
  token = strtok(NULL, " \t");
  float_sw4 cs = 1.0, rho=1.0, cp=sqrt(3.0);
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

  boundaryConditionType bct[6]={bSuperGrid, bSuperGrid, bSuperGrid, bSuperGrid, bSuperGrid, bSuperGrid};
  set_global_bcs(bct);
}

//-----------------------------------------------------------------------
void EW::processTestRayleigh(char* buffer)
{
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("testrayleigh", token) == 0, "ERROR: not a testrayleigh line...: " << token);
  token = strtok(NULL, " \t");
  float_sw4 cs = 1.0, rho=1.0, cp=sqrt(3.0);
  int nwl = 1;
  
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
//                       1234567
    else if (startswith("nwl=", token))
    {
      token += 4; 
      nwl = atoi(token);
      CHECK_INPUT(nwl >= 1, 
		  "Parameter nwl must be >= 1, not: " << nwl);
    }
    else
    {
      badOption("testrayleigh", token);
    }
    token = strtok(NULL, " \t");
  }
// make a test object
  m_rayleigh_wave_test = new TestRayleighWave( rho, cs, cp, nwl, m_global_xmax);

  boundaryConditionType bct[6]={bPeriodic, bPeriodic, bPeriodic, bPeriodic, bStressFree, bDirichlet};
  set_global_bcs(bct);
  
  if (proc_zero())
  {
    float_sw4 Lwave = 2*M_PI/m_rayleigh_wave_test->m_omega;
    float_sw4 Period = Lwave/m_rayleigh_wave_test->m_cr;
    
    printf("TestRayleigh: rho=%e, cp=%e, cs=%e, cr=%e, Wave length=%e, Period=%e\n", 
	   m_rayleigh_wave_test->m_rho, m_rayleigh_wave_test->m_cp, m_rayleigh_wave_test->m_cs, 
	   m_rayleigh_wave_test->m_cr, Lwave, Period );
  }
  
}

//-----------------------------------------------------------------------
void EW::processTestLamb(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("testlamb", token) == 0, "ERROR: not a testlamb line...: " << token);
   token = strtok(NULL, " \t");

   string err = "Testlamb Error: ";
   float_sw4 x0=0.0, y0=0.0, z0=0.0;
   float_sw4 cs = 1.0, rho=1.0, cp=sqrt(3.0), fz=1.0, freq=1.0, f0=1.0; // the exact solution assumes freq = 1

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

   boundaryConditionType bct[6]={bSuperGrid, bSuperGrid, bSuperGrid, bSuperGrid, bStressFree, bSuperGrid};
   set_global_bcs(bct);
}

//-----------------------------------------------------------------------
void EW::processTestEnergy(char* buffer)
{
  char* token = strtok(buffer, " \t");
  string err = "Testenergy Error: ";
  CHECK_INPUT(strcmp("testenergy", token) == 0, "ERROR: not a testenergy line...: " << token);
  token = strtok(NULL, " \t");
  bool use_dirichlet = false;
  bool use_supergrid=false;
  float_sw4 stochastic_amp = 1;
  float_sw4 sg_eps = 1e-4;
  
  int seed=2934839, write_every=1000;
  string filename("energy.log");

  float_sw4 cpcsratio = sqrt(3.0);
  
  while (token != NULL)
  {
    if (startswith("#", token) || startswith(" ", buffer))
      break;

    if (startswith("cpcsratio=", token))
    {
      token += 10; 
      cpcsratio = atof(token);
    }
    else if (startswith("seed=", token))
    {
      token += 5; 
      seed = atoi(token);
    }
    else if (startswith("amplitude=", token))
    {
      token += 10; 
      stochastic_amp = atof(token);
    }
    else if (startswith("sg_eps=", token))
    {
      token += 7; 
      sg_eps = atof(token);
      CHECK_INPUT(sg_eps > 0,
                  err << "testenergy command: sg_eps must be positive, not: " << token);
    }
     else if (startswith("writeEvery=", token))
     {
       token += strlen("writeEvery=");
       write_every = atoi(token);
       CHECK_INPUT(write_every >= 0,
	       err << "testenergy command: writeEvery must be set to a non-negative integer, not: " << token);
     }
    else if (startswith("filename=", token))
    {
      token += 9; 
      filename = token;
    }
    else if( startswith("bchorizontal=",token))
    {
       token += 13;
       use_dirichlet =  strcmp(token,"dirichlet")==0 || strcmp(token,"Dirichlet")==0;
       use_supergrid =  strcmp(token,"supergrid")==0 || strcmp(token,"Supergrid")==0;
    }
    else
    {
       badOption("testenergy", token);
    }
    token = strtok(NULL, " \t");
  }
  m_energy_test = new TestEnergy( seed, cpcsratio, write_every, filename, stochastic_amp, sg_eps );
  // default bc is periodic in the horizontal directions
  boundaryConditionType bct[6]={bPeriodic, bPeriodic, bPeriodic, bPeriodic, bStressFree, bDirichlet};

  if (use_dirichlet)
  {
     for( int side=0 ; side < 4 ; side++ )
	bct[side] = bDirichlet;
  }
  else if (use_supergrid) // supergrid on all sides, except low-z, where we use a free surface bc 
  {
     for( int side=0 ; side < 6 ; side++ )
	bct[side] = bSuperGrid;

     bct[4] = bStressFree;
  }
  
  set_global_bcs(bct);
}

//-----------------------------------------------------------------------
bool EW::checkTestEnergyPeriodic(char* buffer)
{
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("testenergy", token) == 0, "ERROR: not a testenergy line...: " << token);
  token = strtok(NULL, " \t");
  bool use_periodic = true;
  while (token != NULL)
  {
    if (startswith("#", token) || startswith(" ", buffer))
      break;

    if( startswith("bchorizontal=",token))
    {
       token += 13;
       use_periodic =  strcmp(token,"periodic")==0 || strcmp(token,"Periodic")==0;
    }
    token = strtok(NULL, " \t");
  }
  return use_periodic;
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
       if(startswith("path=", token)) {
          token += 5; // skip path=
	  // If path already specified from event lines, skip this path specification.
	  if( m_nevents_specified == 0 )
	  {
             string path=token;
             path += '/';
             mPath.push_back(path);
             //	     mPath[0] = token;
             //	     mPath[0] += '/';
	  }
	  //          path = token;
       }
       else if (startswith("obspath=", token))
       {
          token += 8; // skip obspath=
	  // If obspath already specified from event lines, skip this path specification.
	  if( m_nevents_specified == 0 )
	  {
             string path=token;
             path += '/';
             mObsPath.push_back(path);
             //	     mObsPath[0] = token;
             //	     mObsPath[0] += '/';
	  }
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
       else if (startswith("temppath=", token))
       {
          token += 9; // skip temppath=
          mTempPath = token;
	  mTempPath += '/';
       }
       else
       {
          badOption("fileio", token);
       }
       token = strtok(NULL, " \t");
    }

//  if (path != 0) setOutputPath(path);
  setPrintCycle(printcycle);
  setVerbosity(verbose);
  setParallel_IO(pfs, nwriters);
}

//-----------------------------------------------------------------------
void EW::processGMT(char* buffer)
{
  string filename = "sw4.gmt.csh";
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("gmt", token) == 0, "ERROR: not a gmt line...: " << token);
  token = strtok(NULL, " \t");

  string err = "GMT Error: ";

  while (token != NULL)
    {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	break;
      if (startswith("file=", token))
	{
          token += 5; // skip file=
          filename = token;
       }
      else
	{
          badOption("gmt", token);
       }
      token = strtok(NULL, " \t");
    }
  setGMTOutput(filename, mName); // mName holds the name of the sw4 input file
}

//-----------------------------------------------------------------------
void EW::parsedate( char* datestr, int& year, int& month, int& day, int& hour, int& minute,
		    int& second, int& msecond, int& fail )
{
	  // Format: 01/04/2012:17:34:45.2343 (Month/Day/Year:Hour:Min:Sec.fraction)
   fail = 0;
   int n = strlen(datestr);
   //      cout << "x" << datestr << "x" << endl;
   //   cout << "strlen = " << n << endl;
   int i = 0;
   string buf="";
   while( i<n )
   {
      if( datestr[i]=='/' || datestr[i]==':' || datestr[i] == '.' )
	 buf += datestr[i];
      i++;
   }
   //   if( buf == "//:::." && isdigit(datestr[ifirst]) && isdigit(datestr[n-1]) )
   //   cout << "buf = x" << buf << "x" << endl;
   //   i = 0;
   //   while( !isdigit(datestr[i]) && i < n )
   //      i++;
   if( buf == "//:::." )
   {
      float fsec;
      //      cout << "x" << datestr << "x" << endl;
      sscanf(datestr,"%d/%d/%d:%d:%d:%f",&month,&day,&year,&hour,&minute,&fsec);
      //      cout << " mon " << month << " day " << day << " year " << year << endl;
      //      cout << " hour " << hour<< " minute " << minute << " fsec = " << fsec << endl;
      if( year < 1000 || year > 3000 )
	 fail = 2;
      if( month < 1 || month > 12 )
	 fail = 3;
      if( day < 1 || day > 31 )
	 fail = 4;
      if( hour < 0 || hour > 24 )
	 fail = 5;
      if( minute < 0 || minute > 60 )
	 fail = 6;
      if( fsec < 0 )
	 fail = 8;
      second = static_cast<int>(trunc(fsec));
      msecond = static_cast<int>( round((fsec-second)*1000));
      if( second < 0 || second > 60 )
	 fail = 7;
      //      cout << " second = " << second << " msecond = " << msecond <<endl;
   }
   else 
      fail = 1;
}

//-----------------------------------------------------------------------
void EW::processTime(char* buffer)
{
  float_sw4 t=0.0;
  int steps = -1;
  int year, month, day, hour, minute, second, msecond, fail;
  bool refdateset=false, refeventdateset=false;
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("time", token) == 0, "ERROR: not a time line...: " << token);
  token = strtok(NULL, " \t");
  int event=0;
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
     // Only care about 'event' if event lines are present in input file
       else if(startswith("event=",token))
       {
	  token += 6;
	// Ignore if no events given
	  if( m_nevents_specified > 0 )
	  {
	     map<string,int>::iterator it = m_event_names.find(token);
             //	     CHECK_INPUT( it != m_event_names.end(), 
             //		       err << "event with name "<< token << " not found" );
             if( it != m_event_names.end() )
                event = it->second;
             else if( proc_zero() )
                std::cout << "Time warning: event with name " << token << " not found" << std::endl;
	  }
       }
       else if( startswith("utcstart=",token) )
       {
          token += 9;
	  // Format: 01/04/2012:17:34:45.2343  (Month/Day/Year:Hour:Min:Sec.fraction)
          parsedate( token, year, month, day, hour, minute, second, msecond, fail );
          if( fail == 0 )
	     refdateset = true;
	  else
	     CHECK_INPUT(fail == 0 , "processTime: Error in utcstart format. Give as mm/dd/yyyy:hh:mm:ss.ms, not  " << token );
       }
       else
       {
          badOption("time", token);
       }
       token = strtok(NULL, " \t");
    }
  CHECK_INPUT(!( (t > 0.0) && (steps >= 0) ),
          "Time Error: Cannot set both t and steps for time");
  event = global_to_local_event(event);
  if( event >= 0 )
  {
 // event is handled by this processor group
     if (t > 0.0)
        setGoalTime(t,event);
     else if (steps >= 0)
        setNumberSteps(steps,event);
 
     if( refdateset )
     {
        m_utc0[event][0] = year;
        m_utc0[event][1] = month;
        m_utc0[event][2] = day;
        m_utc0[event][3] = hour;
        m_utc0[event][4] = minute;
        m_utc0[event][5] = second;
        m_utc0[event][6] = msecond;
     }
     else
     {
     // Set UTC as current date
        time_t tsec;
        time( &tsec );
        struct tm *utctime = gmtime( &tsec );
        m_utc0[event][0] = utctime->tm_year+1900;
        m_utc0[event][1] = utctime->tm_mon+1;
        m_utc0[event][2] = utctime->tm_mday;
        m_utc0[event][3] = utctime->tm_hour;
        m_utc0[event][4] = utctime->tm_min;
        m_utc0[event][5] = utctime->tm_sec;
        m_utc0[event][6] = 0; //milliseconds not given by 'time', not needed here.
     }
  }
}

//-----------------------------------------------------------------------
void EW::processBoundaryConditions(char *buffer)
{
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("boundary_conditions", token) == 0, "ERROR: not a boundary condition line...: " << token);
  token = strtok(NULL, " \t");
  
  boundaryConditionType bct[6]={bSuperGrid, bSuperGrid, bSuperGrid, bSuperGrid, bStressFree, bSuperGrid};
  
  int type;
  int side;
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
    case 3:
      bct[side] = bPeriodic;
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

//-----------------------------------------------------------------------
void EW::processSupergrid(char *buffer)
{
  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("supergrid", token) == 0, "ERROR: not a supergrid line...: " << token);
  token = strtok(NULL, " \t");
  int sg_n_gp; // sg_transition;
  float_sw4 sg_coeff, sg_width;
  bool gpSet=false, dampingCoeffSet=false, widthSet=false; // , transitionSet=false
  
  while (token != NULL)
  {
    if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
      break;

//                  1234567890
    if (startswith("gp=", token)) // in number of grid sizes (different from WPP)
    {
      token += 3;
      sg_n_gp = atoi(token);
      CHECK_INPUT(sg_n_gp>0, "The number of grid points in the supergrid damping layer must be positive, not: "<< sg_n_gp);
      gpSet = true;
    }
//                  12345678901
    // else if (startswith("transition=", token)) // in number of grid sizes (different from WPP)
    // {
    //   token += 11;
    //   sg_transition = atoi(token);
    //   CHECK_INPUT(sg_transition>0, "The number of grid points in the supergrid transition layer must be positive, not: "<< sg_transition);
    //   transitionSet = true;
    // }
    else if (startswith("width=", token))
    {
      token += 6;
      sg_width = atof(token);
      CHECK_INPUT(sg_width>0, "The width of the supergrid damping layer must be positive, not: "<< sg_width);
      widthSet = true;
    }
    else if (startswith("dc=", token))
    {
      token += 3;
      sg_coeff = atof(token);
      CHECK_INPUT(sg_coeff>=0., "The supergrid damping coefficient must be non-negative, not: "<<sg_coeff);
      dampingCoeffSet=true;
    }
    else if (startswith("order=", token))
    {
       token += 6;
       int damping = atoi(token);
       CHECK_INPUT( damping == 4 || damping == 6, "The supergrid dissipation order must be 4 or 6, not:" << damping);
       m_sg_damping_order = damping;
    }
    else
    {
      badOption("supergrid", token);
    }
    token = strtok(NULL, " \t");
  } // end while token
  
  CHECK_INPUT( !(gpSet && widthSet), "EW::Processsupergrid, ERROR, both gp and width of supergrid set\n");
     
  if (gpSet)// gp specified
     set_sg_thickness(sg_n_gp);

  if( widthSet )
     set_sg_width( sg_width );

  if (dampingCoeffSet)
    set_sg_damping(sg_coeff);
  else if( m_sg_damping_order == 4 )
     set_sg_damping(0.02);
  else if( m_sg_damping_order == 6 )
     set_sg_damping(0.005);
}

//-----------------------------------------------------------------------
void EW::badOption(string name, char* option) const
{
   if (m_myRank == 0)
      cout << "\tWarning: ignoring " << name << " line option '" << option << "'" << endl;
}

//-----------------------------------------------------------------------
void EW::processGlobalMaterial(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("globalmaterial", token) == 0, "ERROR: not an globalmaterial line...: " << token);
   token = strtok(NULL, " \t");

   string err = "globalmaterial error: ";
   int modelnr = 0;
   float_sw4 frequency = 1;
   float_sw4 vpmin=0, vsmin=0;
  
   while (token != NULL)
   {
      if ( startswith("vpmin=", token) )
      {
 	token += 6;
 	vpmin = atof(token);
      }
      else if ( startswith("vsmin=", token) )
      {
 	token += 6;
 	vsmin = atof(token);
      }
      else
 	badOption("globalmaterial", token);
      token = strtok(NULL, " \t");
   }

   set_threshold_velocities(vpmin, vsmin);
}

//-----------------------------------------------------------------------
// void EW::getEfileInfo(char* buffer)
// {
// #ifdef ENABLE_ETREE
//   // Used only for efiles
//    string accessmode = "parallel";

//    char* token = strtok(buffer, " \t");
//    CHECK_INPUT(strcmp("efile", token) == 0,
// 	       "ERROR: efile info can only be obtained from an efile line, not: " << token);

//    string commandName = token;

//    string err = token;
//    err += " Error: ";

//    token = strtok(NULL, " \t");

//    while (token != NULL)
//    {
//       // while there are tokens in the string still
//       if (startswith("#", token) || startswith(" ", buffer))
// 	 // Ignore commented lines and lines with just a space.
// 	 break;
// // //       else if (startswith("model=", token))
// // //       {
// // //          token += 6; // skip model=
// // //         model = token;
// // //       }
//       else if (startswith("etree=", token))
//       {
// 	 token += 6; // skip etree=
// 	 m_topoFileName = token;
//       }
//       else if (startswith("xetree=", token))
//       {
// 	 token += 7; // skip xetree=
// 	 m_topoExtFileName = token;
//       }
//       else if (startswith("logfile=", token))
//       {
// 	 token += 8; // skip logfile=
//       }
//       else if (startswith("query=", token))
//       {
// 	 token += strlen("query=");
//  	 m_QueryType = token;
// 	 CHECK_INPUT(strcmp(token, "FIXEDRES") == 0 || 
// 		     strcmp(token, "MAXRES") == 0,
// 		     err << "query can only be set to FIXEDRES or MAXRES, not: " << m_QueryType);
//       }
//       else if (startswith("vsmin=", token))
//       {
// 	 token += strlen("vsmin=");
//       }
//       else if (startswith("vpmin=", token))
//       {
// 	 token += strlen("vpmin=");
//       }
//       else if (startswith("access=", token))
//       {
// 	 token += strlen("access=");
// 	 CHECK_INPUT(strcmp(token, "parallel") == 0 ||
//  		     strcmp(token, "serial") == 0,
//  		     err << "access attribute can only be set to serial, or parallel, not: " << token);
// 	 accessmode = token;
//       }
//       else if( startswith("resolution=", token ) )
//       {
//          token += 11;
//          m_EFileResolution = atof(token);
//          CHECK_INPUT(m_EFileResolution>0.,"Resolution must be positive, not " << m_EFileResolution);
//       }
//       else
//       {
// 	 badOption(commandName, token);
//       }
//       token = strtok(NULL, " \t");
//    }
//    // End parsing...
 
// #else
//    CHECK_INPUT(0, "Error: Etree support not compiled into EW (-DENABLE_ETREE)");
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

//-----------------------------------------------------------------------
void EW::processPrefilter(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("prefilter", token) == 0, "ERROR: not a prefilter line...: " << token);
   token = strtok(NULL, " \t");

   string err = "prefilter Error: ";
   string commandName = token;
   float_sw4 fc1=0.1, fc2 = 1.0; // only fc2 is used for low-pass
   FilterType passband = bandPass; // 
   int passes=2; // forwards and backwards gives a zero-phase filter
   int order=2;
   
   while (token != NULL)
   {
      if (startswith("#", token) || startswith(" ", buffer))
         break;

//                    1234567890
      if (startswith("fc1=", token))
      {
        token += 4;
        fc1 = atof(token);
        CHECK_INPUT(fc1>0.,"corner frequency 1 must be positive, not " << fc1 );
      }
//                         1234567890
      else if (startswith("fc2=", token))
      {
        token += 4;
        fc2 = atof(token);
        CHECK_INPUT(fc2>0.,"corner frequency 2 must be positive, not " << fc2 );
      }
//                         1234567890
      else if (startswith("type=", token))
      {
        token += 5;
	if( strcmp(token,"lowpass") == 0 )
	  passband = lowPass;
	else if( strcmp(token,"bandpass") == 0 )
	  passband = bandPass;
	else
	  CHECK_INPUT( false, "processPrefilter: Error: type= " << token << 
		       " Only lowpass or bandpass are recognized\n" );
      }
//                         1234567890
      else if (startswith("passes=", token))
      {
        token += 7;
        passes = atoi(token);
	CHECK_INPUT( passes == 1 || passes == 2, "processPrefilter: Error: passes must be 1 or 2, not = " 
		     << token );
      }
//                         1234567890
      else if (startswith("order=", token))
      {
        token += 6;
        order = atoi(token);
	CHECK_INPUT( order > 0 && order <= 10, "processPrefilter: Error: order = " 
		     << token << " out of bounds\n" );
      }
      else
      {
         badOption(commandName, token);
      }
      token = strtok(NULL, " \t");
   }
   set_prefilter( passband, order, passes, fc1, fc2 );
}

//-----------------------------------------------------------------------
void EW::processGeodynbc(char* buf)
{
   // At this point, the geodyn file has already been read into m_geodyn_filename
   char* token = strtok(buf, " \t");
   CHECK_INPUT(strcmp("geodynbc", token) == 0, "ERROR: not a geodynbc line...: " << token);
   CHECK_INPUT(m_geodynbc_found,"Error: geodynbc not obtained"<<token);

   ifstream geodynfile(m_geodynbc_filename.c_str());
   CHECK_INPUT(geodynfile.is_open(), "Error: opening geodyn file " << m_geodynbc_filename );


   string err = "geodynbc Error: ";
   string commandName = "geodynbc";

   int faces=6, nx=0, ny=0, nz=0, nsteps=0, filter=0, adjust=1;
   float_sw4 x0, y0, z0, lat, lon, elev, az, timestep, rho=0, vs=0, vp=0, freq;
   float_sw4 srcx0, srcy0, srcz0, h, toff;

   bool timestepset = false, nstepsset=false, toffset=false;
   char buffer[256];
   bool done = false;
   while (!geodynfile.eof() && !done )
   {
      geodynfile.getline(buffer,256);
      if (startswith("#", buffer) || startswith("\n", buffer) || buffer == "\0" )
         break;
      if( startswith("begindata",buffer) )
      {
	 done = true;
         break;
      }

      if( startswith("grid", buffer) )
      {
	 char* token = strtok(buffer, " \t");
	 token = strtok(NULL, " \t");
	 while (token != NULL)
	 {
	    if (startswith("#", token) || startswith(" ", buffer))
	       break;
	    if (startswith("faces=", token))
	    {
	       token += 6;
	       faces = atoi(token);
	    }
	    else if( startswith("nx=",token))
	    {
	       token += 3;
	       nx = atoi(token);
	    }
	    else if( startswith("ny=",token))
	    {
	       token += 3;
	       ny = atoi(token);
	    }
	    else if( startswith("nz=",token))
	    {
	       token += 3;
	       nz = atoi(token);
	    }
	    else if( startswith("stepsize=",token))
	    {
	       token += 9;
	       h = atof(token);
	    }
	    else if( startswith("x0=",token))
	    {
	       token += 3;
	       x0 = atof(token);
	    }
	    else if( startswith("y0=",token))
	    {
	       token += 3;
	       y0 = atof(token);
	    }
	    else if( startswith("z0=",token))
	    {
	       token += 3;
	       z0 = atof(token);
	    }
	    else if( startswith("lat=",token))
	    {
	       token += 4;
	       lat = atof(token);
	    }
	    else if( startswith("lon=",token))
	    {
	       token += 4;
	       lon = atof(token);
	    }
	    else if( startswith("elev=",token))
	    {
	       token += 5;
	       elev = atof(token);
	    }
	    else if( startswith("az=",token))
	    {
	       token += 3;
	       az = atof(token);
	    }
	    else if( startswith("adjust=",token))
	    {
	       token += 7;
	       adjust = strcmp(token,"yes")==0;
	    }
	    else
	    {
	       badOption("geodyn-grid", token);
	    }
	    token = strtok(NULL, " \t");
	 }
      }
      else if( startswith("time", buffer) )
      {
	 char* token = strtok(buffer, " \t");
	 token = strtok(NULL, " \t");
	 while (token != NULL)
	 {
	    if (startswith("#", token) || startswith(" ", buffer))
	       break;
	    if (startswith("timestep=", token))
	    {
	       token += 9;
	       timestep = atof(token);
	       timestepset=true;
	    }
	    else if( startswith("nsteps=",token))
	    {
	       token += 7;
	       nsteps = atoi(token);
	       nstepsset=true;
	    }
	    else if( startswith("toff=",token))
	    {
	       token += 5;
	       toff = atof(token);
	       toffset=true;
	    }
	    else
	    {
	       badOption("geodyn-time", token);
	    }
	    token = strtok(NULL, " \t");
	 }
      }
      else if( startswith("material",buffer) )
      {
	 char* token = strtok(buffer, " \t");
	 token = strtok(NULL, " \t");
	 while (token != NULL)
	 {
	    if (startswith("#", token) || startswith(" ", buffer))
	       break;
	    if (startswith("rho=", token))
	    {
	       token += 4;
	       rho = atof(token);
	    }
	    else if( startswith("vs=",token))
	    {
	       token += 3;
	       vs = atof(token);
	    }
	    else if( startswith("vp=",token))
	    {
	       token += 3;
	       vp = atof(token);
	    }
	    else
	    {
	       badOption("geodyn-material", token);
	    }
	    token = strtok(NULL, " \t");
	 }
      }
      else if( startswith("source",buffer) )
      {
	 char* token = strtok(buffer, " \t");
	 token = strtok(NULL, " \t");
	 while (token != NULL)
	 {
	    if (startswith("#", token) || startswith(" ", buffer))
	       break;
	    if (startswith("filter=", token))
	    {
	       token += 7;
	       if( strcmp(token,"butterworth")== 0 )
		  filter = 1;
	       else
		  filter = 0;
	    }
	    else if( startswith("frequency=",token))
	    {
	       token += 10;
	       freq = atof(token);
	    }
	    else if( startswith("x0=",token))
	    {
	       token += 3;
	       srcx0 = atof(token);
	    }
	    else if( startswith("y0=",token))
	    {
	       token += 3;
	       srcy0 = atof(token);
	    }
	    else if( startswith("z0=",token))
	    {
	       token += 3;
	       srcz0 = atof(token);
	    }
	    else
	    {
	       badOption("geodyn-source", token);
	    }
	    token = strtok(NULL, " \t");
	 }
      }
   }
   geodynfile.close();
   CHECK_INPUT( nx == ny, "Geodyn file error: x-y Cube dimensions must be equal, not "
                  <<nx << " " << ny );
   CHECK_INPUT( faces == 5 || faces == 6 , "Geodyn file error: Faces must be 5 or 6, not" << faces );
   CHECK_INPUT( timestepset, "Geodyn file error: No time step given");
   CHECK_INPUT( nstepsset, "Geodyn file error: Number of steps not given");
   set_geodyn_data( m_geodynbc_filename, nx, nz, h, m_ibc_origin, timestep,
		    nsteps, faces );
}

//-----------------------------------------------------------------------
void EW::geodynFindFile(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("geodynbc", token) == 0, "ERROR: not a geodynbc line...: " << token);
   token = strtok(NULL, " \t");
   if( m_events_parallel )
   {
      if( proc_zero() )
      {
         std::cout << "WARNING: geodynbc command does not work together with parallel seismic events" << std::endl;
         std::cout << "geodynbc command will be ignored " << std::endl;
      }
      return;
   }

   string err = "geodynbc Error: ";
   string commandName = token;

   while (token != NULL)
   {
      if (startswith("#", token) || startswith(" ", buffer))
         break;

      if (startswith("file=", token))
      {
         token += 5;
         m_geodynbc_filename = token;
         m_geodynbc_found = true;
      }
      else if( startswith("center=",token))
      {
         token += 7;
         if (atoi(token) == 1 || strcmp(token,"yes")==0 )
	    m_geodynbc_center = true;
      }
      else
      {
         badOption(commandName, token);
      }
      token = strtok(NULL, " \t");
   }
}

//-----------------------------------------------------------------------
void EW::geodynbcGetSizes( string filename, float_sw4 origin[3], float_sw4 &cubelen,
			   float_sw4& zcubelen, float_sw4& hcube, bool& found_latlon,
			   double& lat, double& lon, double& az, int& adjust )
{
   ifstream geodynfile(m_geodynbc_filename.c_str());
   CHECK_INPUT( geodynfile.is_open(), "Error: opening geodyn file " << m_geodynbc_filename );

   string err = "geodynbc Error: ";
   string commandName = "geodynbc";

   int nx=0, ny=0, nz=0, faces=6;
   double x0, y0, z0, elev, h;
   adjust=1;

   char buffer[256];
   bool done = false;
   bool nxfound=false, nyfound=false, nzfound=false, x0found=false, y0found=false, z0found=false;
   bool latfound=false, lonfound=false, azfound=false, hfound=false, elevfound=false;
   while (!geodynfile.eof() && !done )
   {
      geodynfile.getline(buffer,256);
      if (startswith("#", buffer) || startswith("\n", buffer) || buffer == "\0" )
         break;
      if( startswith("begindata",buffer) )
      {
	 done = true;
         break;
      }

      if( startswith("grid", buffer) )
      {
	 char* token = strtok(buffer, " \t");
	 token = strtok(NULL, " \t");
	 while (token != NULL)
	 {
	    if (startswith("#", token) || startswith(" ", buffer))
	       break;
	    if (startswith("faces=", token))
	    {
	       token += 6;
	       faces = atoi(token);
	    }
	    else if( startswith("nx=",token))
	    {
	       token += 3;
	       nx = atoi(token);
               nxfound = true;
	    }
	    else if( startswith("ny=",token))
	    {
	       token += 3;
	       ny = atoi(token);
               nyfound = true;
	    }
	    else if( startswith("nz=",token))
	    {
	       token += 3;
	       nz = atoi(token);
               nzfound = true;
	    }
	    else if( startswith("stepsize=",token))
	    {
	       token += 9;
	       h = atof(token);
               hfound = true;
	    }
	    else if( startswith("x0=",token))
	    {
	       token += 3;
	       x0 = atof(token);
               x0found = true;
	    }
	    else if( startswith("y0=",token))
	    {
	       token += 3;
	       y0 = atof(token);
               y0found = true;
	    }
	    else if( startswith("z0=",token))
	    {
	       token += 3;
	       z0 = atof(token);
               z0found = true;
	    }
	    else if( startswith("lat=",token))
	    {
	       token += 4;
	       lat = atof(token);
               latfound = true;
	    }
	    else if( startswith("lon=",token))
	    {
	       token += 4;
	       lon = atof(token);
               lonfound = true;
	    }
	    else if( startswith("elev=",token))
	    {
	       token += 5;
	       elev = atof(token);
               elevfound = true;
	    }
	    else if( startswith("az=",token))
	    {
	       token += 3;
	       az = atof(token);
               azfound = true;
	    }
	    else if( startswith("adjust=",token))
	    {
	       token += 7;
	       //	       adjust = strcmp(token,"yes")==0;
               adjust = atoi(token);
	    }
	    else
	    {
	       badOption("geodyn-grid", token);
	    }
	    token = strtok(NULL, " \t");
	 }
      }
   }
   geodynfile.close();

   if( nxfound && !nyfound )
   {
      ny = nx;
      nyfound = true;
   }
   if( nxfound && !nzfound )
   {
      nz = nx;
      nzfound = true;
   }
   CHECK_INPUT( (nxfound || nyfound || nzfound) && hfound, "Error in geodyn file: dimensions not specified");

   if( nxfound )
      cubelen = (nx-1)*h;
   else if( nyfound )
      cubelen = (ny-1)*h;
   else
      cubelen = (nz-1)*h;

   found_latlon = latfound && lonfound && azfound;
   if( elevfound )
      origin[2] = -elev;
   else
      origin[2] = z0;

   origin[0] = x0;
   origin[1] = y0;

   zcubelen = cubelen;
   if( nzfound )
      zcubelen = (nz-1)*h;

   if( hfound )
      hcube = h;
}

//-----------------------------------------------------------------------
void EW::processMaterial( char* buffer )
{
   string name = "Material";

   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("material", token) == 0,
 	      "ERROR: material properties can be set by a material line, not: " << token);

   string err = token;
   err += " Error: ";

   token = strtok(NULL, " \t");

   int materialID=-1;
   float_sw4 vp0=-1, vs0=-1, rho0=-1, qp=-1, qs=-1;
   float_sw4 vp1=0, vs1=0, rho1=0;
   float_sw4 vp2=0, vs2=0, rho2=0;
   float_sw4 vp1o2=0, vs1o2=0, rho1o2=0;

   bool gotID = false;
  
   while (token != NULL)
   {
     // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
       // Ignore commented lines and lines with just a space.
	 break;
 //                  1234567890
      if (startswith("id=", token) )
      {
	 token += 3; // skip id=
	 materialID = atoi(token);
	 gotID=true;
      }
// linear variation
      else if (startswith("rhograd=", token))
      {
	 token += 8; // skip rhograd=
	 rho1 = atof(token);
      }
      else if (startswith("vpgrad=", token))
      {
	 token += 7; // skip vpgrad=
	 vp1 = atof(token);
      }
      else if (startswith("vsgrad=", token))
      {
	 token += 7; // skip vsgrad=
	 vs1 = atof(token);
      }
// quadratic variation
      else if (startswith("rho2=", token))
      {
	 token += 5; // skip rho2=
	 rho2 = atof(token);
      }
      else if (startswith("vp2=", token))
      {
	 token += 4; // skip vp2=
	 vp2 = atof(token);
      }
      else if (startswith("vs2=", token))
      {
	 token += 4; // skip vs2=
	 vs2 = atof(token);
      }
// sqrt variation
      else if (startswith("rhosqrt=", token))
      {
	 token += 8; // skip rhosqrt=
	 rho1o2 = atof(token);
      }
      else if (startswith("vpsqrt=", token))
      {
	 token += 7; // skip vpsqrt=
	 vp1o2 = atof(token);
      }
      else if (startswith("vssqrt=", token))
      {
	 token += 7; // skip vssqrt=
	 vs1o2 = atof(token);
      }
// plain vp, vs, rho come last because they start with the same letters as those above...
      else if (startswith("vp=", token) )
      {
	 token += 3; // skip vp=
         vp0 = atof(token);
      }
      else if (startswith("vs=", token) )
      {
         token += 3; // skip vs=
         vs0 = atof(token);
      }
      else if (startswith("rho=", token))
      {
         token += 4; // skip rho=
         rho0 = atof(token);
      }
// attenuation variables
      else if (startswith("Qp=", token) || startswith("qp=", token))
      {
	 token += 3; // skip qp=
	 qp = atof(token);
      }
      else if (startswith("Qs=", token) || startswith("qs=", token))
      {
	 token += 3; // skip qs=
	 qs = atof(token);
      }
      else
      {
	 badOption("material", token);
      }
      token = strtok(NULL, " \t");
   }
 // End parsing...
  
   CHECK_INPUT( gotID, "No id specified in material command");
  

   CHECK_INPUT( (vs0 > 0 || vs1 != 0 || vs2 != 0) , 
 	       "Error in material command: vs0, vs1, vs2 are " << vs0 << " " << vs1 << " " << vs2 );

   CHECK_INPUT( (vp0 > 0 || vp1 != 0 || vp2 != 0) , 
 	       "Error in material command: vp0, vp1, vp2 are " << vp0 << " " << vp1 << " " << vp2 );

   CHECK_INPUT( (rho0 > 0 || rho1 != 0 || rho2 != 0) , 
 	       "Error in material command: rho0, rho1, rho2 are " << rho0 << " " << rho1 << " " << rho2 );

   if( mVerbose >= 2 &&  m_myRank == 0 )
   {
     cout << "**** Material parameters: *****" << endl;
     cout << "materialID=" << materialID << endl;
     cout << "vp=" << vp0 <<  " vpgrad=" << vp1 << " vp2=" << vp2 << " vpsqrt=" << vp1o2 << endl;
     cout << "vs=" << vs0 <<  " vsgrad=" << vs1 << " vs2=" << vs2 << " vssqrt=" << vs1o2 << endl;
     cout << "rho=" << rho0 <<  " rhograd=" << rho1 << " rho2=" << rho2 << " rhosqrt=" << rho1o2 <<endl;
     cout << "qp=" << qp <<  " qs=" << qs << endl;
   }

 // add material to EW object
   MaterialProperty *mat=new MaterialProperty(materialID, vp0, vp1, vp2, vs0, vs1, vs2, rho0, rho1, rho2, qp, qs);
   mat->setSqrtCoefficients( vp1o2, vs1o2, rho1o2 );
   addMaterialProperty( mat );
}

//-----------------------------------------------------------------------
void EW::processSfileOutput( char* buffer )
{
   int cycle=-1, cycleInterval=0;
   int sampleFactorH = 1;
   int sampleFactorV = 1;
   float_sw4 time=0.0, timeInterval=0.0;
   bool timingSet = false;
   float_sw4 tStart = -999.99;
   string filePrefix="sfileoutput";
   bool use_double = false;
  
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("sfileoutput", token) == 0, "ERROR: Not a sfileoutput line...: " << token );

   token = strtok(NULL, " \t");
   string err = "sfileoutput Error: ";
   while (token != NULL)
   {
     // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	 // Ignore commented lines and lines with just a space.
	 break;
      /* if (startswith("time=", token) ) */
      /* { */
	 /* token += 5; // skip time= */
	 /* CHECK_INPUT( atof(token) >= 0.,"Processing sfileoutput command: time must be a non-negative number, not: " << token); */
	 /* time = atof(token); */
	 /* timingSet = true; */
      /* } */
      /* else if (startswith("timeInterval=", token) ) */
      /* { */
	 /* token += 13; // skip timeInterval= */
	 /* CHECK_INPUT( atof(token) >= 0.,"Processing sfileoutput command: timeInterval must be a non-negative number, not: " << token); */
	 /* timeInterval = atof(token); */
	 /* timingSet = true; */
      /* } */
      /* else if (startswith("startTime=", token) ) */
      /* { */
	 /* token += 10; // skip startTime= */
	 /* tStart = atof(token); */
      /* } */
      else if (startswith("sampleFactorH=", token) )
      {
	 token += 14; 
	 CHECK_INPUT( atoi(token) >= 1,"Processing sfileoutput command: sampleFactorH must be a positive integer, not: " << token);
	 sampleFactorH = atoi(token);
      }
      else if (startswith("sampleFactorV=", token) )
      {
	 token += 14; 
	 CHECK_INPUT( atoi(token) >= 1,"Processing sfileoutput command: sampleFactorV must be a positive integer, not: " << token);
	 sampleFactorV= atoi(token);
      }
      else if (startswith("sampleFactor=", token) )
      {
	 token += 13; 
	 CHECK_INPUT( atoi(token) >= 1,"Processing sfileoutput command: sampleFactor must be a positive integer, not: " << token);
	 sampleFactorH = sampleFactorV = atoi(token);
      }
      /* else if (startswith("cycle=", token) ) */
      /* { */
	 /* token += 6; // skip cycle= */
	 /* CHECK_INPUT( atoi(token) >= 0.,"Processing sfileoutput command: cycle must be a non-negative integer, not: " << token); */
	 /* cycle = atoi(token); */
	 /* timingSet = true; */
      /* } */
      /* else if (startswith("cycleInterval=", token) ) */
      /* { */
	 /* token += 14; // skip cycleInterval= */
	 /* CHECK_INPUT( atoi(token) >= 0.,"Processing sfileoutput command: cycleInterval must be a non-negative integer, not: " << token); */
	 /* cycleInterval = atoi(token); */
	 /* timingSet = true; */
      /* } */
      else if (startswith("file=", token))
      {
	 token += 5; // skip file=
	 filePrefix = token;
      }
      else if( startswith("precision=",token) )
      {
	 token += 10;
	 CHECK_INPUT( startswith("double",token) || startswith("float",token),
		      "Processing sfileoutput command: precision must be float or double, not '" << token );
	 use_double =  startswith("double",token);
      }
      else
      {
	 badOption("sfileoutput", token);
      }
      token = strtok(NULL, " \t");
   }

   if( !m_inverse_problem)
   {
      /* CHECK_INPUT( timingSet, "Processing sfileoutput command: " << */ 
		   /* "at least one timing mechanism must be set: cycle, time, cycleInterval or timeInterval"  << endl ); */
      SfileOutput* sfile = new SfileOutput( this, time, timeInterval, cycle, cycleInterval, 
 			       tStart, filePrefix, sampleFactorH, sampleFactorV, use_double);
      addSfileOutput( sfile );
   }
}


//-----------------------------------------------------------------------
void EW::processImage3D( char* buffer )
{
   int cycle=-1, cycleInterval=0;
   Image3D::Image3DMode mode=Image3D::RHO;
   float_sw4 time=0.0, timeInterval=0.0;
   bool timingSet = false;
   float_sw4 tStart = -999.99;
   string filePrefix="volimage";
   bool use_double = false;
  
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("volimage", token) == 0, "ERROR: Not a volimage line...: " << token );

   token = strtok(NULL, " \t");
   string err = "volimage Error: ";
   while (token != NULL)
   {
     // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	 // Ignore commented lines and lines with just a space.
	 break;
      if (startswith("time=", token) )
      {
	 token += 5; // skip time=
	 CHECK_INPUT( atof(token) >= 0.,"Processing volimage command: time must be a non-negative number, not: " << token);
	 time = atof(token);
	 timingSet = true;
      }
      else if (startswith("timeInterval=", token) )
      {
	 token += 13; // skip timeInterval=
	 CHECK_INPUT( atof(token) >= 0.,"Processing volimage command: timeInterval must be a non-negative number, not: " << token);
	 timeInterval = atof(token);
	 timingSet = true;
      }
      else if (startswith("startTime=", token) )
      {
	 token += 10; // skip startTime=
	 tStart = atof(token);
      }
      else if (startswith("cycle=", token) )
      {
	 token += 6; // skip cycle=
	 CHECK_INPUT( atoi(token) >= 0.,"Processing volimage command: cycle must be a non-negative integer, not: " << token);
	 cycle = atoi(token);
	 timingSet = true;
      }
      else if (startswith("cycleInterval=", token) )
      {
	 token += 14; // skip cycleInterval=
	 CHECK_INPUT( atoi(token) >= 0.,"Processing volimage command: cycleInterval must be a non-negative integer, not: " << token);
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
	 if (strcmp(token, "ux") == 0)        mode = Image3D::UX;
	 else if (strcmp(token, "uy") == 0)   mode = Image3D::UY;
	 else if (strcmp(token, "uz") == 0)   mode = Image3D::UZ;
	 else if (strcmp(token, "rho") == 0)   mode = Image3D::RHO;
	 else if (strcmp(token, "p") == 0)   mode = Image3D::P;
	 else if (strcmp(token, "s") == 0)   mode = Image3D::S;
	 else if (strcmp(token, "mu") == 0)   mode = Image3D::MU;
	 else if (strcmp(token, "lambda") == 0)   mode = Image3D::LAMBDA;
	 else if (strcmp(token, "gradrho") == 0)   mode = Image3D::GRADRHO;
	 else if (strcmp(token, "gradp") == 0)   mode = Image3D::GRADP;
	 else if (strcmp(token, "grads") == 0)   mode = Image3D::GRADS;
	 else if (strcmp(token, "gradmu") == 0)   mode = Image3D::GRADMU;
	 else if (strcmp(token, "gradlambda") == 0)   mode = Image3D::GRADLAMBDA;
	 else if (strcmp(token, "qs") == 0)   mode = Image3D::QS;
	 else if (strcmp(token, "qp") == 0)   mode = Image3D::QP;
	 else
	 {
	    //	    mode = static_cast<Image3D::Image3DMode>(atoi(token));
	    CHECK_INPUT(0,"Processing image3D command: " << "mode must be one of the following: " << endl <<
			"\t ux|uy|uz|rho|p|s|mu|lambda|gradrho|gradp|grads|gradmu|gradlambda|qs|qp *not: "<< token );
	 }
      }
      else if( startswith("precision=",token) )
      {
	 token += 10;
	 CHECK_INPUT( startswith("double",token) || startswith("float",token),
		      "Processing volimage command: precision must be float or double, not '" << token );
	 use_double =  startswith("double",token);
      }
      else
      {
	 badOption("volimage", token);
      }
      token = strtok(NULL, " \t");
   }
   bool forwardgrad = !m_inverse_problem && (mode == Image3D::GRADRHO ||mode == Image3D::GRADMU ||mode == Image3D::GRADLAMBDA ||
					     mode == Image3D::GRADP   ||mode == Image3D::GRADS);
   if( forwardgrad && proc_zero() )
   {
      cout << "WARNING: volume images of material gradients can not be computed by the forward solver" << endl;
      cout << "   volimage will not be created " << endl;
   }

   bool attenuation = (mode == Image3D::QS || mode == Image3D::QP );
  
   if (attenuation && !m_use_attenuation && proc_zero() )
   {
     cout << "ERROR: volume images of Qs or Qp can only be generated when attenuation is enabled" << endl;
     MPI_Abort( MPI_COMM_WORLD, 1 );
   }
   

   if( !forwardgrad )
   {
      CHECK_INPUT( timingSet, "Processing volimage command: " << 
		   "at least one timing mechanism must be set: cycle, time, cycleInterval or timeInterval"  << endl );
      Image3D* im3 = new Image3D( this, time, timeInterval, cycle, cycleInterval, 
 			       tStart, filePrefix, mode, use_double );
      addImage3D( im3 );
   }
}

//-----------------------------------------------------------------------
void EW::processESSI3D( char* buffer )
{
   int dumpInterval=-1, bufferInterval=1;
   string filePrefix="ssioutput";
   float_sw4 coordValue;
   float_sw4 coordBox[4];
   const float_sw4 zero=0.0;
   int precision = 8;
   int compressionMode = 0;
   double compressionPar;
   
   // Default is whole domain
   coordBox[0] = zero;
   coordBox[1] = m_global_xmax;
   coordBox[2] = zero;

   coordBox[3] = m_global_ymax;
   float_sw4 depth = -999.99; // default not specified

   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("ssioutput", token) == 0 || strcmp("essioutput", token) == 0, "ERROR: Not a essioutput/ssioutput line...: " << token );

   token = strtok(NULL, " \t");
   string err = "ssioutput Error: ";
   while (token != NULL)
   {
     // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
         // Ignore commented lines and lines with just a space.
         break;
      else if (startswith("file=", token))
      {
         token += 5; // skip file=
         filePrefix = token;
      }
      else if (startswith("dumpInterval=", token))
      {
          token += 13; // skip dumpInterval=
          dumpInterval = atoi(token);
      }
      else if (startswith("bufferInterval=", token))
      {
          token += 15; // skip bufferInterval=
          bufferInterval = atoi(token);
      }
      else if (startswith("xmin=", token))
      {
          token += 5; // skip xmin=
          coordValue = atof(token);
          coordBox[0] = min(m_global_xmax, coordValue);
          coordBox[0] = max(zero, coordBox[0]);
      }
      else if (startswith("xmax=", token))
      {
          token += 5; // skip xmax=
          coordValue = atof(token);
          coordBox[1] = min(m_global_xmax, coordValue);
          coordBox[1] = max(zero, coordBox[1]);
      }
      else if (startswith("ymin=", token))
      {
          token += 5; // skip ymin=
          coordValue = atof(token);
          coordBox[2] = min(m_global_ymax, coordValue);
          coordBox[2] = max(zero, coordBox[2]);
      }
      else if (startswith("ymax=", token))
      {
          token += 5; // skip ymax=
          coordValue = atof(token);
          coordBox[3] = min(m_global_ymax, coordValue);
          coordBox[3] = max(zero, coordBox[3]);
      }
      else if (startswith("depth=", token))
      {
          token += 6; // skip depth=
          coordValue = atof(token);
          depth = min(m_global_zmax, coordValue);
          depth = max(zero, depth);
      }
      else if (startswith("precision=", token))
      {
          token += 10; // skip precision=
          precision = atoi(token);
          if (precision != 4 && precision != 8) 
              badOption("ssioutput precision", token);
          if (proc_zero())
            cout << "SSI ouput will use " << precision*8 << "-bit floating point values." << endl;
      }
      else if (startswith("zfp-rate=", token))
      {
          token += 9;
          compressionMode = SW4_ZFP_MODE_RATE;
          compressionPar = atof(token);
          if (proc_zero())
            cout << "SSI ouput will use ZFP rate=" << compressionPar << endl;
      }
      else if (startswith("zfp-precision=", token))
      {
          token += 14;
          compressionMode = SW4_ZFP_MODE_PRECISION;
          compressionPar = atof(token);
          if (proc_zero())
            cout << "SSI ouput will use ZFP precision=" << compressionPar << endl;
      }
      else if (startswith("zfp-accuracy=", token))
      {
          token += 13;
          compressionMode = SW4_ZFP_MODE_ACCURACY;
          compressionPar = atof(token);
          if (proc_zero())
            cout << "SSI ouput will use ZFP accuracy=" << compressionPar << endl;
      }
      else if (startswith("zfp-reversible=", token))
      {
          token += 15;
          compressionMode = SW4_ZFP_MODE_REVERSIBLE;
          if (proc_zero())
            cout << "SSI ouput will use ZFP reversible mode" << endl;
      }
      else if (startswith("zlib=", token))
      {
          token += 5;
          compressionMode = SW4_ZLIB;
          compressionPar = atof(token);
          if (proc_zero())
            cout << "SSI ouput will use ZLIB level=" << (int)compressionPar << endl;
      }
      else if (startswith("szip=", token))
      {
          token += 5;
          compressionMode = SW4_SZIP;
          if (proc_zero())
            cout << "SSI ouput will use SZIP" << endl;
      }
      else if (startswith("sz=", token))
      {
          token += 5;
          compressionMode = SW4_SZ;
          if (proc_zero())
            cout << "SSI ouput will use SZ with configuration file [" << getenv("SZ_CONFIG_FILE") << "]" << endl;
      }
      else
      {
          badOption("ssioutput", token);
      }
      token = strtok(NULL, " \t");
   }

#ifndef USE_ZFP
   if (compressionMode != 0 && proc_zero()) 
      cout << "WARNING: SW4 is not compiled with ZFP but ZFP command is used " << endl;
#endif

#ifndef USE_SZ
   if (compressionMode != 0 && proc_zero()) 
      cout << "WARNING: SW4 is not compiled with SZ but SZ command is used " << endl;
#endif

   if (compressionMode != 0 && bufferInterval == 1) 
     bufferInterval = 100;

   // Check the specified min/max values make sense
   for (int d=0; d < 2*2; d+=2)
   {
      if (coordBox[d+1] < coordBox[d])
      {
         char coordName[2] = {'x','y'};
         if (proc_zero())
           cout << "ERROR: ssioutput subdomain " << coordName[d] <<
              " coordinate max value " << coordBox[d+1] <<
              " is less than min value " << coordBox[d] << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }

   // Use depth if zmin/zmax values are specified
   if ((depth < 0) && proc_zero())
   {
      cout << "WARNING: ssioutput depth not specified or less than zero, setting to zero" << endl;
      depth=0;
   }

   ESSI3D* essi3d = new ESSI3D( this, filePrefix, dumpInterval, bufferInterval, coordBox, depth, precision, compressionMode, compressionPar);
   addESSI3D( essi3d );
}

//-----------------------------------------------------------------------
void EW::processCheckPoint(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("checkpoint", token) == 0, "ERROR: not a checkpoint line...: " << token);
   token = strtok(NULL, " \t");
   string err = "CheckPoint Error: ";
   int cycle=-1, cycleInterval=0;
   float_sw4 time=0.0, timeInterval=0.0;
   bool timingSet=false;
   string filePrefix = "checkpoint";

   string restartFileName, restartPath;
   bool   restartFileGiven = false, restartPathGiven=false;

   size_t bufsize=10000000;

   while (token != NULL)
    {
       if (startswith("#", token) || startswith(" ", buffer))
          break;
       //      if (startswith("cycle=", token) )
       //      {
       //	 token += 6; // skip cycle=
       //	 CHECK_INPUT( atoi(token) >= 0., err << "cycle must be a non-negative integer, not: " << token);
       //	 cycle = atoi(token);
       //	 timingSet = true;
       //      }
      if (startswith("cycleInterval=", token) )
      {
	 token += 14; // skip cycleInterval=
	 CHECK_INPUT( atoi(token) >= 0., err << "cycleInterval must be a non-negative integer, not: " << token);
	 cycleInterval = atoi(token);
	 timingSet = true;
      }
      else if (startswith("file=", token))
      {
	 token += 5; // skip file=
	 filePrefix = token;
      }
      else if (startswith("restartfile=", token))
      {
	 token += 12; // skip file=
	 restartFileName = token;
	 restartFileGiven = true;
      }
      else if (startswith("restartpath=", token))
      {
	 token += 12;
	 restartPath = token;
	 restartPathGiven = true;
      }
      else if (startswith("bufsize=", token))
      {
	 token += 8; // skip bufsize=
	 bufsize = atoi(token);
      }
      else
      {
	 badOption("checkpoint", token);
      }
      token = strtok(NULL, " \t");
   }
   if( m_check_point == CheckPoint::nil )
      m_check_point = new CheckPoint(this);
   if( cycleInterval > 0 )
      m_check_point->set_checkpoint_file( filePrefix, cycle, cycleInterval, bufsize );
   if( restartFileGiven )
   {
      m_check_point->set_restart_file( restartFileName, bufsize );
   }
   if( restartPathGiven )
   {
     m_check_point->set_restart_path( restartPath );
   }
}

//-----------------------------------------------------------------------
void EW::setOutputPath(const string& path) 
{ 
  stringstream s;
  s << path << "/";
  mPath[0] = s.str();
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
void EW::setGoalTime(float_sw4 t, int event ) 
{ 
  mTmax[event] = t; 
  mTstart = 0.0; 
  mTimeIsSet[event] = true;
}

//-----------------------------------------------------------------------
void EW::setNumberSteps(int steps,int event)
{
  mNumberOfTimeSteps[event] = steps;
  mTimeIsSet[event] = false;
}

//-----------------------------------------------------------------------
int EW::getNumberOfSteps(int event) const
{
  return mNumberOfTimeSteps[event];
}

//-----------------------------------------------------------------------
int EW::getNumberOfEvents() const
{
  return m_nevent;
}

//-----------------------------------------------------------------------
int EW::getNumberOfLocalEvents() const
{
  return m_eEnd-m_eStart+1;
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
void EW::set_cflnumber( float_sw4 cfl )
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
void EW::allocateCartesianSolverArrays(float_sw4 a_global_zmax)
{
//
// note that this routine might modify the value of m_global_zmax
//
   if (mVerbose>=2 && proc_zero())
     printf("allocateCartesianSolverArrays: #ghost points=%d, #parallel padding points=%d, topoExists=%s\n", m_ghost_points,
            m_ppadding, m_topography_exists? "true":"false");

// z=0 is the last element in m_refinementBoundaries[]
   int nCurvilinearGrids = 0;
   int nCartGrids = 0;
   
//    m_topography_exists indicates if there was a topography command in the input file   
   if (m_topography_exists)
   {
      nCurvilinearGrids= m_curviRefLev.size();
   }

   nCartGrids = m_refinementBoundaries.size(); // There is always one ref boundary (at z=0)
   
   if (mVerbose>=2 && proc_zero())
      printf("refBndrSize= %lu, nCartGrids=%d, nCurviGrids=%d \n", m_refinementBoundaries.size(), nCartGrids, nCurvilinearGrids);
      
   int refFact = 1;
// Cartesian refinements
   for( int r = 0 ; r < nCartGrids-1 ; r++ )
   {
      refFact *= 2;
      //      cout << "refinement boundary " << r << " is " << m_refinementBoundaries[r] << endl;
   }

// Curvilinear refinements
   for( int r = 0 ; r < nCurvilinearGrids-1 ; r++ )
   {
      refFact *= 2;
      //      cout << "refinement boundary " << r << " is " << m_curviRefLev[r] << endl;
   }
   
// is there an attenuation command in the file?
   if (!m_use_attenuation)
      m_number_mechanisms = 0;

   int is_periodic[2]={0,0};

// some test cases, such as testrayleigh uses periodic boundary conditions in the x and y directions
   if (m_doubly_periodic)
   {
      is_periodic[0]=1;
      is_periodic[1]=1;
   }

// m_nx_base, m_ny_base: number of grid points in the coarsest grid: assigned by processGrid()
   int nx_finest_w_ghost = refFact*(m_nx_base-1)+1+2*m_ghost_points;
   int ny_finest_w_ghost = refFact*(m_ny_base-1)+1+2*m_ghost_points;
   if( is_periodic[0] )
      nx_finest_w_ghost = refFact*m_nx_base + 2*m_ghost_points;
   if( is_periodic[1] )
      ny_finest_w_ghost = refFact*m_ny_base + 2*m_ghost_points;

   int proc_max[2];
   bool old_decomp=true;
// this info is obtained by the contructor
//   MPI_Comm_size( MPI_COMM_WORLD, &nprocs  );
   if( old_decomp )
   {
      proc_decompose_2d( nx_finest_w_ghost, ny_finest_w_ghost, m_nProcs, proc_max );
      MPI_Cart_create( m_1d_communicator, 2, proc_max, is_periodic, true, &m_cartesian_communicator );
   }
   else
   {
      int mynewid=my_node_core_id( nx_finest_w_ghost, ny_finest_w_ghost, proc_max );  
      MPI_Comm renumbered_world;
      MPI_Comm_split( MPI_COMM_WORLD, 0, mynewid, &renumbered_world );   
      MPI_Cart_create( renumbered_world, 2, proc_max, is_periodic, true, &m_cartesian_communicator );
      std::cout << "old/new ranks " << getRank() << "/" << mynewid << std::endl;
   }


   int my_proc_coords[2];
   MPI_Cart_get( m_cartesian_communicator, 2, proc_max, is_periodic, my_proc_coords );
   MPI_Cart_shift( m_cartesian_communicator, 0, 1, m_neighbor, m_neighbor+1 );
   MPI_Cart_shift( m_cartesian_communicator, 1, 1, m_neighbor+2, m_neighbor+3 );

   if( proc_zero_evzero() && mVerbose >= 1 /*3*/) // tmp
   {
     cout << " Grid distributed on " << m_nProcs << " processors " << endl;
     cout << " Finest grid size    " << nx_finest_w_ghost << " x " << ny_finest_w_ghost << endl;
     cout << " Processor array     " << proc_max[0] << " x " << proc_max[1] << endl;
   }
// save the cartesian processor decomposition
   m_proc_array[0] = proc_max[0];
   m_proc_array[1] = proc_max[1];

// the domain decomposition is done for the finest grid
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
      mNumberOfGrids+=nCurvilinearGrids;
//      mNumberOfGrids++;

// tmp
   if (proc_zero_evzero())
   {
      cout << " Number of curvilinear grids = " << nCurvilinearGrids << endl;
      cout << " Number of Cartesian grids = " << mNumberOfCartesianGrids << endl;
      cout << " Total number of grids = " << mNumberOfGrids << endl;
   }   

   m_iscurvilinear.resize(mNumberOfGrids);
   for( int g=0 ; g < mNumberOfCartesianGrids ; g++ )
      m_iscurvilinear[g] = false;
//   if( m_topography_exists )
   for( int g=mNumberOfCartesianGrids ; g < mNumberOfGrids ; g++ )
      m_iscurvilinear[g] = true;

   m_supergrid_taper_x.resize(mNumberOfGrids);
   m_supergrid_taper_y.resize(mNumberOfGrids);
   m_supergrid_taper_z.resize(mNumberOfGrids);
   mMu.resize(mNumberOfGrids);
   mLambda.resize(mNumberOfGrids);
   mRho.resize(mNumberOfGrids);
   mC.resize(mNumberOfGrids);
   m_zmin.resize(mNumberOfGrids);
   mGridSize.resize(mNumberOfGrids);
   mMinVsOverH.resize(mNumberOfGrids);

// curvilinear arrays
   mX.resize(mNumberOfGrids);
   mY.resize(mNumberOfGrids);
   mZ.resize(mNumberOfGrids);
   mMetric.resize(mNumberOfGrids);
   mJ.resize(mNumberOfGrids);

// coefficients for Mesh refinement
   m_Morf.resize(mNumberOfGrids);
   m_Mlrf.resize(mNumberOfGrids);
   m_Mufs.resize(mNumberOfGrids);
   m_Mlfs.resize(mNumberOfGrids);

   m_Morc.resize(mNumberOfGrids);
   m_Mlrc.resize(mNumberOfGrids);
   m_Mucs.resize(mNumberOfGrids);
   m_Mlcs.resize(mNumberOfGrids);

// always allocate the pointer arrays for the viscoelastic properties (allows. e.g., mQs[g].is_defined() to be called)
   mQs.resize(mNumberOfGrids);
   mQp.resize(mNumberOfGrids);

// viscoelastic material coefficients and memory variables are only allocated when attenuation is enabled
// Allocate pointers, even if attenuation not used, to avoid segfault in parameter list with mMuVE[g], etc...
   mMuVE.resize(mNumberOfGrids);
   mLambdaVE.resize(mNumberOfGrids);
   if (m_use_attenuation && m_number_mechanisms > 0) // the simplest model only uses Q, not MuVe, LambdaVE, or OmegaVE
   {
     mOmegaVE.resize(m_number_mechanisms); // global relaxation frequencies (1 per mechanism)
     
// muVE and lambdaVE are vectors of vectors
     for (int g=0; g<mNumberOfGrids; g++)
     {
       mMuVE[g]     = new Sarray[m_number_mechanisms];
       mLambdaVE[g] = new Sarray[m_number_mechanisms];
     }
   }
   
   m_iStart.resize(mNumberOfGrids);
   m_iEnd.resize(mNumberOfGrids);
   m_jStart.resize(mNumberOfGrids);
   m_jEnd.resize(mNumberOfGrids);
   m_kStart.resize(mNumberOfGrids);
   m_kEnd.resize(mNumberOfGrids);

   m_iStartAct.resize(mNumberOfGrids);
   m_iEndAct.resize(mNumberOfGrids);
   m_jStartAct.resize(mNumberOfGrids);
   m_jEndAct.resize(mNumberOfGrids);
   m_kStartAct.resize(mNumberOfGrids);
   m_kEndAct.resize(mNumberOfGrids);

   m_iStartActGlobal.resize(mNumberOfGrids);
   m_iEndActGlobal.resize(mNumberOfGrids);
   m_jStartActGlobal.resize(mNumberOfGrids);
   m_jEndActGlobal.resize(mNumberOfGrids);
   m_kStartActGlobal.resize(mNumberOfGrids);
   m_kEndActGlobal.resize(mNumberOfGrids);

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

   m_NumberOfBCPoints.resize(mNumberOfGrids);
   m_BndryWindow.resize(mNumberOfGrids);

   int *wind;
   
   for( int g=0; g < mNumberOfGrids ; g++ )
   {
     m_NumberOfBCPoints[g] = new int[6];
     m_BndryWindow[g] = new int [36]; // 6 by 6 array in Fortran
     for (int side=0; side < 6; side++)
     {
       m_NumberOfBCPoints[g][side]=0;
       for (int qq=0; qq<6; qq+=2) // 0, 2, 4
	 m_BndryWindow[g][qq + side*6]= 999;
       for (int qq=1; qq<6; qq+=2) // 1, 3, 5
	 m_BndryWindow[g][qq + side*6]= -999;
     }
   }
   
   float_sw4 h = m_h_base;

// save the grid spacing for all Cartesian grids
   for (int g=0; g<nCartGrids; g++)
   {
     mGridSize[g] = h;
     h = h/2.;
   }

// NEW 3/13/2018
   // save the horizontal grid spacing for all curvilinear grids
   h = mGridSize[nCartGrids-1]; // bottom curvilinear grid has same grid size as the top Cartesian
   for (int g=nCartGrids; g<mNumberOfGrids; g++)
   {
     mGridSize[g] = h;
     h = h/2.;
   }
   
   m_global_nx[mNumberOfGrids-1] = nx_finest_w_ghost-2*m_ghost_points;
   m_global_ny[mNumberOfGrids-1] = ny_finest_w_ghost-2*m_ghost_points;
   
// Grid size in the curvilinear portion
   for (int g=mNumberOfGrids-2; g >=nCartGrids; g--) // the coarsest curvilinear grid has grid number 'nCartGrids'
   {
      if( is_periodic[0] )
	 m_global_nx[g] = m_global_nx[g+1]/2;
      else
	 m_global_nx[g] = 1 + (m_global_nx[g+1]-1)/2;
      if( is_periodic[1] )
	 m_global_ny[g] = m_global_ny[g+1]/2;
      else
	 m_global_ny[g] = 1 + (m_global_ny[g+1]-1)/2;
   }

   if (!m_topography_exists)
   { // finest grid on top
      m_global_nx[nCartGrids-1] = nx_finest_w_ghost-2*m_ghost_points;
      m_global_ny[nCartGrids-1] = ny_finest_w_ghost-2*m_ghost_points;
   }
   else
   { // top Cartesian grid has the same grid size as the bottom curvilinear grid
      m_global_nx[nCartGrids-1] = m_global_nx[nCartGrids];
      m_global_ny[nCartGrids-1] = m_global_ny[nCartGrids];
   }
   
// previous code for Cartesian MR
   for (int g=nCartGrids-2; g >=0; g--)
   {
      if( is_periodic[0] )
	 m_global_nx[g] = m_global_nx[g+1]/2;
      else
	 m_global_nx[g] = 1 + (m_global_nx[g+1]-1)/2;
      if( is_periodic[1] )
	 m_global_ny[g] = m_global_ny[g+1]/2;
      else
	 m_global_ny[g] = 1 + (m_global_ny[g+1]-1)/2;
   }

// the curvilinear grid has a variable grid size, but matches the finest Cartesian grid where they meet
   // if( m_topography_exists )
   // {
   //   mGridSize[mNumberOfGrids-1]   = mGridSize[mNumberOfGrids-2];
   //   m_global_nx[mNumberOfGrids-1] = m_global_nx[mNumberOfGrids-2];
   //   m_global_ny[mNumberOfGrids-1] = m_global_ny[mNumberOfGrids-2];
   // }
   
// Define grid in z-direction, by formula z_k = (k-1)*h + zmin
   vector<int> nz;
   nz.resize(nCartGrids);

// don't change the zmin of the finest cartesian grid
   m_zmin[nCartGrids-1] = m_refinementBoundaries[nCartGrids-1];
   for( int g = nCartGrids-1; g >= 0; g-- )
   {
     float_sw4 zmax = (g>0? m_refinementBoundaries[g-1]:a_global_zmax);
     nz[g]     = 1 + static_cast<int> ( (zmax - m_zmin[g])/mGridSize[g] + 0.5 );

     zmax = m_zmin[g] + (nz[g]-1)*mGridSize[g];
      
     m_global_nz[g] = nz[g]; // save the number of grid points in the z-direction

     if( g>0 )
        m_zmin[g-1] = zmax; // save m_zmin[g-1] for next iteration
     else
       a_global_zmax = zmax;
   }

// extent of computational grid (without ghost points)
   if (!m_doubly_periodic)
   {
     m_global_xmax = mGridSize[nCartGrids-1]*(m_global_nx[nCartGrids-1] - 1);
     m_global_ymax = mGridSize[nCartGrids-1]*(m_global_ny[nCartGrids-1] - 1);
   }
   else
   {
     m_global_xmax = mGridSize[nCartGrids-1]*(m_global_nx[nCartGrids-1]);
     m_global_ymax = mGridSize[nCartGrids-1]*(m_global_ny[nCartGrids-1]);
   }
   
   m_global_zmax = m_zmin[0] + (nz[0]-1)*mGridSize[0];
   if (mVerbose >= 1 && proc_zero_evzero())
     cout << "Extent of the computational domain xmax=" << m_global_xmax << " ymax=" << m_global_ymax << " zmax=" << 
       m_global_zmax << endl;

// detailed grid refinement info
   if( proc_zero_evzero() && mVerbose >= 1 /*2*/) // tmp
   {
      cout << "Cartesian refinement levels after correction: " << endl;

      for( int g=0; g<nCartGrids; g++ )
      {
         cout << "Grid=" << g << " z-min=" << m_zmin[g] << endl;
      }
      cout << "Corrected global_zmax = " << m_global_zmax << endl << endl;
   }
   
// Allocate the topography arrays and coarsen out the grid in the curvilinear portion of the grid
   if( m_topography_exists ) // UPDATED  for more than 1 curvilinear grid
   {
// Allocate elements in the m_curviInterface vector
      m_curviInterface.resize(mNumberOfGrids - mNumberOfCartesianGrids);

// NEW
      for (int g=mNumberOfGrids-1; g>=mNumberOfCartesianGrids; g--) // g=mNumberOfGrids-1 is the finest curvilinear grid
      {
// save the local index bounds
         m_iStart[g] = ifirst;
         m_iEnd[g]   = ilast; // finest grid size in x from above
         m_jStart[g] = jfirst;
         m_jEnd[g]   = jlast; // finest grid size in y from above
// k-bounds must be determined by the grid generator

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

         if (jlast == ny + m_ghost_points)
            m_jEndInt[g]   = ny;
         else
            m_jEndInt[g]   = jlast - m_ppadding;

         // m_kStartInt[g] = 1;
         // m_kEndInt[g]   = nz[g];

// check that there are more interior points than padding points
         if (m_iEndInt[g] - m_iStartInt[g] + 1 < m_ppadding)
         {
            printf("WARNING: less interior points than padding in proc=%d, grid=%d, m_iStartInt=%d, "
                   "m_iEndInt=%d, padding=%d\n", m_myRank, g, m_iStartInt[g], m_iEndInt[g], m_ppadding);
         }
         if (m_jEndInt[g] - m_jStartInt[g] + 1 < m_ppadding)
         {
            printf("WARNING: less interior points than padding in proc=%d, grid=%d, m_jStartInt=%d, "
                   "m_jEndInt=%d, padding=%d\n", m_myRank, g, m_jStartInt[g], m_jEndInt[g], m_ppadding);
         }
      
// output bounds
         if (proc_zero_evzero() && mVerbose >=1 /*3*/) // tmp
         {
            printf("Rank=%d, Grid #%d (curvilinear), iInterior=[%d,%d], jInterior=[%d,%d]\n", m_myRank, g, m_iStartInt[g], m_iEndInt[g],
                   m_jStartInt[g], m_jEndInt[g]);
         }

         // number of extra ghost points to allow highly accurate interpolation; needed for the source discretization
         m_ext_ghost_points = 8;

// Allocate interface the interface surface for this curvilinear grid
         m_curviInterface[g-mNumberOfCartesianGrids].define(m_iStart[g]-m_ext_ghost_points, m_iEnd[g]+m_ext_ghost_points,
                                                            m_jStart[g]-m_ext_ghost_points, m_jEnd[g]+m_ext_ghost_points,1,1);

// Allocate topo arrays for the top (finest) curvilinear grid
         if (g==mNumberOfGrids-1)
         {
// 2 versions of the topography:
            mTopo.define(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],1,1); // true topography/bathymetry, read directly from etree/rfile
//            mTopo.define(ifirst,ilast,jfirst,jlast,1,1); // true topography/bathymetry, read directly from etree
// smoothed version of true topography, with an extended number (4 instead of 2 ) of ghost points.
            mTopoGridExt.define(m_iStart[g]-m_ext_ghost_points, m_iEnd[g]+m_ext_ghost_points,
                                m_jStart[g]-m_ext_ghost_points, m_jEnd[g]+m_ext_ghost_points,1,1);
            // mTopoGridExt.define(ifirst-m_ext_ghost_points,ilast+m_ext_ghost_points,
            //                     jfirst-m_ext_ghost_points,jlast+m_ext_ghost_points,1,1);
         }

// At this point, we don't know the number of grid points in the k-direction of the curvi-linear grid.
// The arrays mX, mY, mZ, etc, must be allocated by the grid generator

// Any other arrays that we need to allocate here?

// Go to the next coarser grid, unless this is the coarsest curvilinear grid
         if (g > mNumberOfCartesianGrids)
         {
            coarsen1d( nx, ifirst, ilast, is_periodic[0] );
            coarsen1d( ny, jfirst, jlast, is_periodic[1] );
         }
      } // end for all curvilinear grids
       
   } // end if m_topography_exists

// Define Cartesian grid arrays, loop from finest to coarsest

// On entry to this loop: (nx, ifirst, ilast) and (ny, jfirst, jlast) are the local number of grid points, starting, and ending indices
// of the finest Cartesian grif 
   for( int g = nCartGrids-1 ; g >= 0 ; g-- )
   {
// NOTE: same number of ghost points in all directions
      kfirst     = 1-m_ghost_points;
      klast      = nz[g] + m_ghost_points;

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

      if (jlast == ny + m_ghost_points)
	m_jEndInt[g]   = ny;
      else
	m_jEndInt[g]   = jlast - m_ppadding;

      m_kStartInt[g] = 1;
      m_kEndInt[g]   = nz[g];

// check that there are more interior points than padding points
      if (m_iEndInt[g] - m_iStartInt[g] + 1 < m_ppadding)
      {
         printf("WARNING: less interior points than padding in proc=%d, grid=%d, m_iStartInt=%d, "
                "m_iEndInt=%d, padding=%d\n", m_myRank, g, m_iStartInt[g], m_iEndInt[g], m_ppadding);
      }
      if (m_jEndInt[g] - m_jStartInt[g] + 1 < m_ppadding)
      {
         printf("WARNING: less interior points than padding in proc=%d, grid=%d, m_jStartInt=%d, "
                "m_jEndInt=%d, padding=%d\n", m_myRank, g, m_jStartInt[g], m_jEndInt[g], m_ppadding);
      }
      
// output bounds
      if (proc_zero() && mVerbose >=1 /*3*/) // tmp
      {
         printf("Rank=%d, Grid #%d (Cartesian), iInterior=[%d,%d], jInterior=[%d,%d], kInterior=[%d,%d]\n", m_myRank, g, m_iStartInt[g], m_iEndInt[g],
                m_jStartInt[g], m_jEndInt[g], m_kStartInt[g], m_kEndInt[g]);
      }
      
//
// Allocate arrays as needed by the use case
//
      mRho[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
      mRho[g].set_to_minusOne();
      if( m_anisotropic )
      {
	 mC[g].define(21,ifirst,ilast,jfirst,jlast,kfirst,klast);
	 mC[g].set_to_minusOne();
      }
      else
      {
// elastic material
	 mMu[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	 mLambda[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
// initialize the material coefficients to -1
	 mMu[g].set_to_minusOne();
	 mLambda[g].set_to_minusOne();
// allocate space for material coefficient arrays needed by MR
	 m_Morc[g].define(ifirst,ilast,jfirst,jlast,1,1);
	 m_Mlrc[g].define(ifirst,ilast,jfirst,jlast,1,1);
	 m_Mucs[g].define(ifirst,ilast,jfirst,jlast,1,1);
	 m_Mlcs[g].define(ifirst,ilast,jfirst,jlast,1,1);

	 int nkf = m_global_nz[g];
	 m_Morf[g].define(ifirst,ilast,jfirst,jlast,nkf,nkf);
	 m_Mlrf[g].define(ifirst,ilast,jfirst,jlast,nkf,nkf);
	 m_Mufs[g].define(ifirst,ilast,jfirst,jlast,nkf,nkf);
	 m_Mlfs[g].define(ifirst,ilast,jfirst,jlast,nkf,nkf);
      }
// viscoelastic material coefficients & memory variables
      if (m_use_attenuation)
      {
	mQs[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	mQp[g].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	for (int a=0; a<m_number_mechanisms; a++) // the simplest attenuation model only uses Q, not MuVE or LambdaVE
	{
	  mMuVE[g][a].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
	  mLambdaVE[g][a].define(ifirst,ilast,jfirst,jlast,kfirst,klast);
// initialize the viscoelastic material coefficients to -1
	  mMuVE[g][a].set_to_minusOne();
	  mLambdaVE[g][a].set_to_minusOne();
	}
// initialize Qp and Qs to -1
	mQs[g].set_to_minusOne();
	mQp[g].set_to_minusOne();
      }
	
      // go to the next coarser grid 
      coarsen1d( nx, ifirst, ilast, is_periodic[0] );
      coarsen1d( ny, jfirst, jlast, is_periodic[1] );
      //      cout << g << " " << my_proc_coords[0] << " I Split into " << ifirst << " , " << ilast << endl;
      //      cout << g << " " << my_proc_coords[1] << " J Split into " << jfirst << " , " << jlast << endl;
      //      cout << "grid " << g << " zmin = " << m_zmin[g] << " nz = " << nz[g] << " kinterval " << kfirst << " , " << klast << endl;

   }// end for all Cartesian grids

// tmp
//   int myRank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//   for (int q=0; q<mNumberOfCartesianGrids; q++)
//   {
//     printf("Proc #%i, m_iEnd[%i]=%i, m_global_nx[%i]=%i, m_jEnd[%i]=%i, m_global_ny[%i]=%i\n", 
// 	   myRank, q, m_iEnd[q], q, m_global_nx[q], q, m_jEnd[q], q, m_global_ny[q]);
   
//   }

} // end allocateCartesianSolverArrays()


//-----------------------------------------------------------------------
void EW::allocateCurvilinearArrays()
{
   // This routine should define sizes and allocate arrays on the curvilinear grid
   if (!m_topography_exists ) return;

   if (mVerbose >= 1 && proc_zero())
      cout << "***inside allocateCurvilinearArrays***"<< endl;

// 1: get the min and max elevation from mTopoGridExt

//   int gTop = mNumberOfGrids-1;
//   int ifirst = m_iStart[gTop];
//   int ilast  = m_iEnd[gTop];
//   int jfirst = m_jStart[gTop];
//   int jlast  = m_jEnd[gTop];
   float_sw4 h = mGridSize[mNumberOfGrids-1]; // grid size must agree with top cartesian grid
//   float_sw4 zTopCart = m_zmin[g]; // bottom z-level for curvilinear grid
//   float_sw4 zTopCart = m_topo_zmax; // bottom z-level for curvilinear grid


// decide on the number of grid point in the k-direction (evaluate mTopoGrid...)
   float_sw4 zMinLocal, zMinGlobal, zMaxLocal, zMaxGlobal;
   //   int i=m_iStart[gTop], j=m_jEnd[gTop];
// note that the z-coordinate points downwards, so positive elevation (above sea level)
// has negative z-values
//   zMaxLocal = zMinLocal = -mTopoGridExt(i,j,1);
// tmp
//   int i_min_loc=i, i_max_loc=i;
//   int j_min_loc=j, j_max_loc=j;
// end tmp
// the mTopoGridExt array was allocated in allocateCartesianSolverArrays()   
   zMinLocal = -mTopoGridExt.maximum();
   zMaxLocal = -mTopoGridExt.minimum();
   MPI_Allreduce( &zMinLocal, &zMinGlobal, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
   MPI_Allreduce( &zMaxLocal, &zMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator);


   //   for (i= imin ; i<=imax ; i++)
   //      for (j=jmin; j<=jmax ; j++)
   //      {
   //	 if (-mTopoGridExt(i,j,1) > zMaxLocal)
   //	 {
   //	    zMaxLocal = -mTopoGridExt(i,j,1);
   //            i_max_loc = i;
   //            j_max_loc = j;
   //	 }
      
   //	 if (-mTopoGridExt(i,j,1) < zMinLocal)
   //	 {
   //	    zMinLocal = -mTopoGridExt(i,j,1);
   //            i_min_loc = i;
   //            j_min_loc = j;
   //	 }
   //      }
// tmp
//   printf("Proc #%i: zMaxLocal = %e at (%i %i), zMinLocal = %e at (%i %i)\n", m_myRank, zMaxLocal, i_max_loc, j_max_loc,
// 	 zMinLocal, i_min_loc, j_min_loc);
// end tmp

// Compute some un-divided differences of the topographic surface to evaluate its smoothness
   int imin = mTopoGridExt.m_ib;
   int imax = mTopoGridExt.m_ie;
   int jmin = mTopoGridExt.m_jb;
   int jmax = mTopoGridExt.m_je;
   float_sw4 maxd2zh=0, maxd2z2h=0, maxd3zh=1.e-20, maxd3z2h=1.e-20, d2h, d3h, h3=h*h*h;
// grid size h
   for (int i=imin+1; i<=imax-1; i++)
      for (int j=jmin+1; j<=jmax-1; j++)
      {
         d2h = sqrt( SQR((mTopoGridExt(i-1,j,1) - 2*mTopoGridExt(i,j,1) + mTopoGridExt(i+1,j,1))/1.) + 
                     SQR((mTopoGridExt(i,j-1,1) - 2*mTopoGridExt(i,j,1) + mTopoGridExt(i,j+1,1))/1.) + 
                     SQR((mTopoGridExt(i+1,j+1,1) - mTopoGridExt(i+1,j,1) - mTopoGridExt(i,j+1,1) + mTopoGridExt(i,j,1))/1.) );
         if (d2h > maxd2zh) maxd2zh = d2h;
      }
// 3rd differences
   for (int i=imin+1; i<=imax-2; i++)
      for (int j=jmin+1; j<=jmax-2; j++)
      {
         d3h = sqrt( SQR((mTopoGridExt(i-1,j,1) - 3*mTopoGridExt(i,j,1) + 3*mTopoGridExt(i+1,j,1) - mTopoGridExt(i+2,j,1))/1.) + 
                     SQR((mTopoGridExt(i,j-1,1) - 3*mTopoGridExt(i,j,1) + 3*mTopoGridExt(i,j+1,1) - mTopoGridExt(i,j+2,1))/1.) + 
                     SQR((mTopoGridExt(i+2,j+1,1) - mTopoGridExt(i+2,j,1) - 2*mTopoGridExt(i+1,j+1,1) + 2*mTopoGridExt(i+1,j,1) + 
                          mTopoGridExt(i,j+1,1) - mTopoGridExt(i,j,1))/1.) +
                     SQR((mTopoGridExt(i+1,j+2,1) - 2*mTopoGridExt(i+1,j+1,1) - mTopoGridExt(i,j+2,1) + 2*mTopoGridExt(i,j+1,1) + 
                          mTopoGridExt(i+1,j,1) - mTopoGridExt(i,j,1))/1.) );

         if (d3h > maxd3zh) maxd3zh = d3h;
      }
// grid size 2h
   for (int i=imin+2; i<=imax-2; i+=2)
      for (int j=jmin+2; j<=jmax-2; j+=2)
      {
         d2h = sqrt( SQR((mTopoGridExt(i-2,j,1) - 2*mTopoGridExt(i,j,1) + mTopoGridExt(i+2,j,1))/1.) + 
                     SQR((mTopoGridExt(i,j-2,1) - 2*mTopoGridExt(i,j,1) + mTopoGridExt(i,j+2,1))/1.) + 
                     SQR((mTopoGridExt(i+2,j+2,1) - mTopoGridExt(i+2,j,1) - mTopoGridExt(i,j+2,1) + mTopoGridExt(i,j,1))/1.) );
         if (d2h > maxd2z2h) maxd2z2h = d2h;
      }
// 3rd differences
   for (int i=imin+2; i<=imax-4; i+=2)
      for (int j=jmin+2; j<=jmax-4; j+=2)
      {
         d3h = sqrt( SQR((mTopoGridExt(i-2,j,1) - 3*mTopoGridExt(i,j,1) + 3*mTopoGridExt(i+2,j,1) - mTopoGridExt(i+4,j,1))/1.) + 
                     SQR((mTopoGridExt(i,j-2,1) - 3*mTopoGridExt(i,j,1) + 3*mTopoGridExt(i,j+2,1) - mTopoGridExt(i,j+4,1))/1.) + 
                     SQR((mTopoGridExt(i+4,j+2,1) - mTopoGridExt(i+4,j,1) - 2*mTopoGridExt(i+2,j+2,1) + 2*mTopoGridExt(i+2,j,1) + 
                          mTopoGridExt(i,j+2,1) - mTopoGridExt(i,j,1))/1.) +
                     SQR((mTopoGridExt(i+2,j+4,1) - 2*mTopoGridExt(i+2,j+2,1) - mTopoGridExt(i,j+4,1) + 2*mTopoGridExt(i,j+2,1) + 
                          mTopoGridExt(i+2,j,1) - mTopoGridExt(i,j,1))/1.) );
         if (d3h > maxd3z2h) maxd3z2h = d3h;
      }
   float_sw4 d2zh_global=0, d2z2h_global=0, d3zh_global=0, d3z2h_global=0;
   MPI_Allreduce( &maxd2zh,  &d2zh_global,  1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
   MPI_Allreduce( &maxd2z2h, &d2z2h_global, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
   MPI_Allreduce( &maxd3zh,  &d3zh_global,  1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
   MPI_Allreduce( &maxd3z2h, &d3z2h_global, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);

   float_sw4 topo_zmax = m_gridGenerator->get_topo_zmax();
   if(proc_zero() )
   {
      printf("\n");
      printf("***Topography grid: min z = %e, max z = %e, top Cartesian z = %e\n", zMinGlobal, zMaxGlobal, topo_zmax );
      if (mVerbose >= 3)
      {
         printf("***Un-divided differences of grid surface (ratio h*D2/D3 should be close to the same for h and 2h):\n"
                "h*max D2z(h)   = %e, max D3z(h)  = %e, Ratio: h*max D2z(h) / max D3z(h)    = %e\n"
                "2h*max D2z(2h) = %e, max D3z(2h) = %e, Ratio: 2h*max D2z(2h) / max D3z(2h) = %e\n", 
                h*d2zh_global, d3zh_global, h*d2zh_global/d3zh_global, 2*h*d2z2h_global, d3z2h_global,
                2*h*d2z2h_global/d3z2h_global);
         printf("\n");
      }
   }
// remember the global zmin
   m_global_zmin = zMinGlobal; // = -(highest elevation)
   CHECK_INPUT(topo_zmax>zMaxGlobal,"allocateCurvilinearArrays: Negative thickness of curvilinear grid.\n"
               "Increase topography zmax to exceed zMaxGlobal = "<< zMaxGlobal <<", preferrably by at least "
               << zMaxGlobal-zMinGlobal);

// 2: determine the number of grid points in each curvilinear grid

// set last element of m_curviRefLev to equal average
   float_sw4 avg_minZ = 0.5*(zMaxGlobal+zMinGlobal);

   // if (proc_zero())
   // {
   //    for (int q=0; q<m_refinementBoundaries.size(); q++)
   //       printf("m_refBndry[%d] = %e\n", q, m_refinementBoundaries[q]);
   // }

// Use the m_zmin array to help keep track of the top (min z) coordinate for each grid
// assigned for the Cartesian portion of the grid in allocateCartesianSolverArrays()   
   
   if (mVerbose >= 3 && proc_zero())
   {
      for (size_t q=0; q<m_curviRefLev.size(); q++)
         printf("m_curviRefLev[%lu]=%e\n", q, m_curviRefLev[q]);
   }
   
// scale the refinement levels to take the average topographic elevation into account
   for (int g=mNumberOfCartesianGrids; g<mNumberOfGrids; g++)
   {
      m_zmin[g] = avg_minZ + m_curviRefLev[g - mNumberOfCartesianGrids] * (topo_zmax - avg_minZ)/topo_zmax;
   }
   
   if (mVerbose >= 3 && proc_zero())
   {
      for (int g=0; g<mNumberOfGrids; g++)
         printf("m_zmin[%d] = %e\n", g, m_zmin[g]);
   }

//
// loop over all curvilinear grids and allocate space + estimate the number of grid points in z   
//
   for (int g = mNumberOfCartesianGrids; g <mNumberOfGrids; g++)
   {
// on average the same gridsize in z
      int Nz = 1+ (int) ((m_zmin[g-1] - m_zmin[g])/mGridSize[g]); 
      m_kStart[g] = 1 - m_ghost_points;
      m_kEnd[g]  = Nz + m_ghost_points;
      m_global_nz[g] = Nz;
      m_kStartInt[g] = 1;
      m_kEndInt[g]   = Nz;
      if(mVerbose >= 3 && proc_zero() )
         printf("allocateCurvilinearArrays: Number of grid points in curvilinear grid[%d] = %i, kStart = %i, kEnd = %i\n", 
                g, Nz, m_kStart[g], m_kEnd[g]);
//
// NOTE: mX, mY, mZ, etc are of type vector<Sarray>
// allocate mX, mY, and mZ arrays
      mX[g].define(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g]);
      mY[g].define(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g]);
      mZ[g].define(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g]);
// Allocate array for the metric
      mMetric[g].define(4,m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
// and the Jacobian of the transformation
      mJ[g].define(m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
      mX[g].set_to_zero();
      mY[g].set_to_zero();
      mZ[g].set_to_zero();      
      mMetric[g].set_to_zero();
      mJ[g].set_to_zero();

// and material properties, initialize to -1
      mRho[g].define(m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
      mRho[g].set_to_minusOne();

      if( m_anisotropic ) // NEED TO UPDATE FOR SEVERAL CURVILINEAR GRIDS!!!
      {
         mC[g].define(21,m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
         mCcurv.define(45,m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
         mC[g].set_to_minusOne();
         mCcurv.set_to_minusOne();
      }
      else
      {
         mMu[g].define(m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
         mMu[g].set_to_minusOne();
         mLambda[g].define(m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
         mLambda[g].set_to_minusOne();
      }
// viscoelastic material coefficients
      if (m_use_attenuation)
      {
// initialize the viscoelastic material coefficients to -1
         mQs[g].define(m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
         mQs[g].set_to_minusOne();
         mQp[g].define(m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
         mQp[g].set_to_minusOne();
         for (int a=0; a<m_number_mechanisms; a++) // the simplest attenuation model has m_number_mechanisms = 0
         {
            mMuVE[g][a].define(m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
            mMuVE[g][a].set_to_minusOne();
            mLambdaVE[g][a].define(m_iStart[g],m_iEnd[g],m_jStart[g],m_jEnd[g],m_kStart[g],m_kEnd[g]);
            mLambdaVE[g][a].set_to_minusOne();
         } // end for a...
      }// end if attenuation
   } // end for g...
}// end allocateCurvilinearArrays()

//-----------------------------------------------------------------------
void EW::deprecatedImageMode(int value, const char* name) const
{
  if (m_myRank == 0)
    cout << "***Warning specifying the mode using integers is deprecated, mode="
	 << value << " should be mode=" << name << " instead." << endl;
}

//-----------------------------------------------------------------------
void EW::processImage(char* buffer, bool use_hdf5)
{
   int cycle=-1, cycleInterval=0;
//   int pfs = 0, nwriters=1;
   Image::ImageMode mode=Image::RHO;
   float_sw4 time=0.0, timeInterval=0.0;
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
  float_sw4 coordValue;
  int gridPointValue;
  bool coordWasSet = false;
  bool use_double = false;
  bool mode_is_grid = false;
  
  char* token = strtok(buffer, " \t");
  if ( strcmp("image", token) != 0 && strcmp("imagehdf5", token) != 0 )
  {
    cerr << "Processing image command: " << "ERROR: not an image line...: " << token;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  
  token = strtok(NULL, " \t");

  string err = "Image Error: ";

  if (mVerbose >=4 && proc_zero())
    cout << "********Parsing image command*********" << endl;
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
      else if (strcmp(token, "uxexact") == 0)   mode = Image::UXEXACT;
      else if (strcmp(token, "uyexact") == 0)   mode = Image::UYEXACT;
      else if (strcmp(token, "uzexact") == 0)   mode = Image::UZEXACT;
      else if (strcmp(token, "p") == 0)   mode = Image::P;
      else if (strcmp(token, "s") == 0)   mode = Image::S;
      else if (strcmp(token, "div") == 0)   mode = Image::DIV;
      else if (strcmp(token, "curl") == 0)   mode = Image::CURLMAG;
//      else if (strcmp(token, "curlmag") == 0)   mode = Image::CURLMAG; // strange name, remove ??
      else if (strcmp(token, "veldiv") == 0)   mode = Image::DIVDT;
      else if (strcmp(token, "divdudt") == 0)   mode = Image::DIVDT;
      else if (strcmp(token, "velcurl") == 0)   mode = Image::CURLMAGDT;
      else if (strcmp(token, "curldudt") == 0)   mode = Image::CURLMAGDT;
      else if (strcmp(token, "lat") == 0)   mode = Image::LAT;
      else if (strcmp(token, "lon") == 0)   mode = Image::LON;
      else if (strcmp(token, "hvelmax") == 0) mode = Image::HMAXDUDT;
      else if (strcmp(token, "hmaxdudt") == 0) mode = Image::HMAXDUDT;
      else if (strcmp(token, "hmax") == 0) mode = Image::HMAX;
      else if (strcmp(token, "vvelmax") == 0) mode = Image::VMAXDUDT;
      else if (strcmp(token, "vmaxdudt") == 0) mode = Image::VMAXDUDT;
      else if (strcmp(token, "vmax") == 0) mode = Image::VMAX;
      else if (strcmp(token, "topo") == 0) mode = Image::TOPO;
      else if (strcmp(token, "grid") == 0) mode_is_grid=true;
      else if (strcmp(token, "gridx") == 0) mode = Image::GRIDX;
      else if (strcmp(token, "gridy") == 0) mode = Image::GRIDY;
      else if (strcmp(token, "gridz") == 0) mode = Image::GRIDZ;
      else if (strcmp(token, "uxerr") == 0)   mode = Image::UXERR;
      else if (strcmp(token, "uyerr") == 0)   mode = Image::UYERR;
      else if (strcmp(token, "uzerr") == 0)   mode = Image::UZERR;
      // else if (strcmp(token, "fx") == 0)   mode = Image::FX;
      // else if (strcmp(token, "fy") == 0)   mode = Image::FY;
      // else if (strcmp(token, "fz") == 0)   mode = Image::FZ;
      else if (strcmp(token, "velmag") == 0)   mode = Image::MAGDUDT;
      else if (strcmp(token, "magdudt") == 0)   mode = Image::MAGDUDT;
      else if (strcmp(token, "mag") == 0)   mode = Image::MAG;
      else if (strcmp(token, "hvelmag") == 0)   mode = Image::HMAGDUDT;
      else if (strcmp(token, "hmagdudt") == 0)   mode = Image::HMAGDUDT;
      else if (strcmp(token, "hmag") == 0)   mode = Image::HMAG;
      else if (strcmp(token, "gradrho") == 0)   mode = Image::GRADRHO;
      else if (strcmp(token, "gradmu") == 0)   mode = Image::GRADMU;
      else if (strcmp(token, "gradlambda") == 0)   mode = Image::GRADLAMBDA;
      else if (strcmp(token, "gradp") == 0)   mode = Image::GRADP;
      else if (strcmp(token, "grads") == 0)   mode = Image::GRADS;
      else if (strcmp(token, "qp") == 0) mode = Image::QP;
      else if (strcmp(token, "qs") == 0) mode = Image::QS;
      //      else if (strcmp(token, "hvel") == 0) mode = Image::HVEL;
      else
      {
	  cerr << "Processing image command: " << "mode must be one of the following: " << endl
	       << "ux|uy|uz|rho|lambda|mu" << endl 
               << "|p|s|div|curl|veldiv|divdudt|velcurl|curldudt " << endl
	       << "|lat|lon|hmaxdudt|hvelmax|hmax|vmaxdudt|vvelmax|vmax|topo|grid|gridx|gridy|gridz " << endl
	       << "|magdudt|velmag|mag|hvelmag|hmagdudt|hmag" << endl
	       << "|uxexact|uyexact|uzexact|uxerr|uyerr|uzerr|gradrho|gradmu|gradlambda|gradp|grads|qp|qs|" << endl
	       << "*not: " << token << endl;
	  MPI_Abort( MPI_COMM_WORLD, 1 );
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
	cerr << "Processing image command: " << "cannot set a coordinate location twice, x and " << 
	  locationType << " were both set." << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      coordWasSet=true;
      locationType=Image::X;
      coordValue = atof(token);
      if ( coordValue < 0.0 || coordValue > m_global_xmax )
      {
	cerr << "Processing image command: " << "x value must be within the computational domain 0<=x<=" 
	     << m_global_xmax << ", not x=: " << coordValue << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
    }
     
    else if (startswith("y=", token))
    {
      token += 2; // skip y=
      if ( coordWasSet )
      {
	cerr << "Processing image command: " << "cannot set a coordinate location twice, y and " 
	     << locationType << " were both set." << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      coordWasSet=true;
      locationType=Image::Y;
      coordValue = atof(token);
      if ( coordValue < 0.0 || coordValue > m_global_ymax )
      {
	cerr << "Processing image command: " << "y value must be within the computational domain 0<=y<=" 
	     << m_global_ymax << ", not y= " << coordValue << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
    }
    else if (startswith("z=", token))
    {
      token += 2; // skip z=
      if ( coordWasSet )
      {
	cerr << "Processing image command: " << "cannot set a coordinate location twice, z and " 
	     << locationType << " were both set." << endl;
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      coordWasSet=true;
      locationType=Image::Z;
      coordValue = atof(token);
      if ( coordValue < 0.0 || coordValue > m_global_zmax )
      {
	cerr << "Processing image command: " << "z value must be within the computational domain 0<=z<=" 
	     << m_global_zmax << ", not z= " << coordValue << endl;
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



   bool forwardgrad = !m_inverse_problem && (mode == Image::GRADRHO ||mode == Image::GRADMU ||mode == Image::GRADLAMBDA ||
					     mode == Image::GRADP   ||mode == Image::GRADS);
   if( forwardgrad && proc_zero() )
   {
      cout << "WARNING: images of material gradients can not be computed by the forward solver" << endl;
      cout << "   image will not be created " << endl;
   }
   if( !forwardgrad )
   {
  // Set up the image object
      Image* i;
      if (coordWasSet)
      {
	 if( mode_is_grid )
	 {
	    if( locationType == Image::X )
	    {
	       i = new Image(this, time, timeInterval, cycle, cycleInterval, 
			 filePrefix, Image::GRIDY, locationType, coordValue, use_double, use_hdf5);
	       addImage(i);
	       i = new Image(this, time, timeInterval, cycle, cycleInterval, 
			 filePrefix, Image::GRIDZ, locationType, coordValue, use_double, use_hdf5);
	       addImage(i);
	    }
	    else if( locationType == Image::Y )
	    {
	       i = new Image(this, time, timeInterval, cycle, cycleInterval, 
			 filePrefix, Image::GRIDX, locationType, coordValue, use_double, use_hdf5);
	       addImage(i);
	       i = new Image(this, time, timeInterval, cycle, cycleInterval, 
			 filePrefix, Image::GRIDZ, locationType, coordValue, use_double, use_hdf5);
	       addImage(i);
	    }
	    else if( locationType == Image::Z )
	    {
	       i = new Image(this, time, timeInterval, cycle, cycleInterval, 
			 filePrefix, Image::GRIDX, locationType, coordValue, use_double, use_hdf5);
	       addImage(i);
	       i = new Image(this, time, timeInterval, cycle, cycleInterval, 
			     filePrefix, Image::GRIDY, locationType, coordValue, use_double, use_hdf5);
	       addImage(i);
	    }
	 }
	 else
	 {
	    i = new Image(this, time, timeInterval, cycle, cycleInterval, 
			  filePrefix, mode, locationType, coordValue, use_double, use_hdf5);
	    addImage(i);
	 }
      }
      else 
      {
	 cerr << "Processing image command: " << "one of the coordinate (x,y,z) option must be set to determine the image's 2D plane" << endl;
	 MPI_Abort( MPI_COMM_WORLD, 1 );
      }
   }

  if (mVerbose >=4 && proc_zero())
    cout << "********Done parsing image command*********" << endl;
}

// int sgn(double arg)
// {
//   if (arg < 0)
//     return -1;
//   else 
//     return 1;
// }

//----------------------------------------------------------------------------
void EW::processSource(char* buffer, vector<vector<Source*> > & a_GlobalUniqueSources )
{

  Source* sourcePtr;
  
  float_sw4 m0 = 1.0;
  float_sw4 t0=0.0, f0=1.0, freq=1.0;
  // Should be center of the grid
  double x = 0.0, y = 0.0, z = 0.0;
  int i = 0, j = 0, k = 0;
  float_sw4 mxx=0.0, mxy=0.0, mxz=0.0, myy=0.0, myz=0.0, mzz=0.0;
  float_sw4 strike=0.0, dip=0.0, rake=0.0;
  float_sw4 fx=0.0, fy=0.0, fz=0.0;
  int isMomentType = -1;
  
  double lat = 0.0, lon = 0.0, depth = 0.0;
  bool topodepth = false, depthSet=false, zSet=false;
  
  bool cartCoordSet = false;
  bool geoCoordSet = false;
  bool strikeDipRake = false;
  bool dfileset=false;
  bool sacbaseset = false;

  int ncyc = 0;
  int event=0;
  bool ncyc_set = false;

  timeDep tDep = iRickerInt;
  char formstring[100];
  char dfile[1000];

  strcpy(formstring, "Ricker");

  char* token = strtok(buffer, " \t");
  REQUIRE2(strcmp("source", token) == 0, "ERROR: not a source line...: " << token);
  token = strtok(NULL, " \t");

  string err = "Source Error: ";

  string cartAndGeoErr = "source command: Cannot set both a geographical (lat,lon) and cartesian coordinate (x,y)";
  string pointAndMomentErr = "source command: Cannot set both a point source and moment tensor formulation";

  if (mVerbose >=4 && proc_zero())
    cout << "********Parsing source command*********" << endl;
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
         token += 2; // skip z=
// with topography, the z-coordinate can have both signs!
         z = atof(token);
	 topodepth=false; // this is absolute depth
	 
         zSet = true;
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
      else if (startswith("depth=", token)) // this is the same as topodepth: different from WPP
      {
         token += 6; // skip depth=
         depth = atof(token);
	 topodepth = true;
         CHECK_INPUT(depth >= 0.0,
		     err << "source command: Depth below topography must be greater than or equal to zero");
	 depthSet=true;
      }
//                         1234567890
      else if (startswith("topodepth=", token))
      {
         token += 10; // skip depth=
         depth = atof(token);
	 topodepth = true;
         CHECK_INPUT(depth >= 0.0,
                 err << "source command: Depth below topography must be greater than or equal to zero");
// by depth we here mean depth below topography
	 depthSet=true;
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
      else if(startswith("event=",token))
      {
	 token += 6;
	 //	 event = atoi(token);
	 //	 CHECK_INPUT( 0 <= event && event < m_nevent, err << "event no. "<< event << " out of range" );
	// Ignore if no events given
	 if( m_nevents_specified > 0 )
	 {
	    map<string,int>::iterator it = m_event_names.find(token);
            //	    CHECK_INPUT( it != m_event_names.end(), 
            //		     err << "event with name "<< token << " not found" );
             if( it != m_event_names.end() )
                event = it->second;
             else if( proc_zero() )
                std::cout << "Source warning: event with name " << token << " not found" << std::endl;

	 }
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
         strncpy(formstring, token,100);
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
         else if (!strcmp("Dirac",formstring) )
	    tDep = iDirac;
         else if (!strcmp("C6SmoothBump",formstring) )
	    tDep = iC6SmoothBump;
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
      else if (startswith("dfile=",token))
      {
         token += 6;
         strncpy(dfile, token,1000);
	 dfileset = true;
      }
      else if (startswith("sacbase=",token))
      {
         token += 8;
         strncpy(dfile, token,1000);
	 sacbaseset = true;
	 isMomentType = 1;
      }
      else if (startswith("sacbasedisp=",token))
      {
         token += 12;
         strncpy(dfile, token,1000);
	 sacbaseset = true;
	 isMomentType = 0;
      }
      else
      {
         badOption("source", token);
      }
      token = strtok(NULL, " \t");
    }
  if( event_is_in_proc(event) )
  {
    int elocal=global_to_local_event(event);
  // Set up source on wpp object.
  CHECK_INPUT(cartCoordSet || geoCoordSet,
	  err << "source command: cartesian or geographic coordinate must be specified");

  CHECK_INPUT(depthSet || zSet,
	  err << "source command: depth, topodepth or z-coordinate must be specified");

  if( tDep == iGaussianWindow )
      CHECK_INPUT( ncyc_set, err << "source command: ncyc must be set for Gaussian Window function");

  // Discrete source time function
  float_sw4* par=NULL;
  int* ipar=NULL;
  int npar=0, nipar=0;
  if( dfileset )
  {
     tDep = iDiscrete;
     //  g(t) defined by spline points on a uniform grid, read from file.
     //  Format: t0, dt, npts
     //          g_1
     //          g_2
     //         ....
     FILE* fd=fopen(dfile, "r" );
     CHECK_INPUT( fd !=NULL , err << "Source time function file " << dfile << " not found" );
     float_sw4 t0, dt;
     int npts;
     fscanf(fd," %lg %lg %i", &t0, &dt, &npts );
     par = new float_sw4[npts+1];
     par[0]  = t0;
     freq    = 1/dt;
     ipar    = new int[1];
     ipar[0] = npts;
     for( int i=0 ; i < npts ; i++ )
	fscanf(fd,"%lg", &par[i+1] );
     npar = npts+1;
     nipar = 1;
     //     cout << "Read disc source: t0=" << t0 << " dt="  << dt << " npts= " << npts << endl;
     fclose(fd);
  }
  if( sacbaseset )
  {
     // Read moment tensor components or forcing components from sac files.

     // Set constant tensor to identity. m0 can give some scaling.
     mxx = myy = mzz = 1;
     mxy = mxz = myz = 0;

     // Set first force components to identity. f0 can give some scaling.
     fx = 1;
     fy = fz = 0;
     
     bool timereverse = false; // Reverse the SAC data. Set true for testing purpose only, the users want to do this themselves outside SW4.
     bool useB = false; // Use sac header begin time parameter B.


     float_sw4 dt, t0, latsac, lonsac,cmpazsac, cmpincsac;
     int utcsac[7], npts;
     string basename = dfile;
     string fname;
     if( isMomentType )
     {
	tDep = iDiscrete6moments;
	fname = basename + ".xx";
	npar  = 6*(npts+1);
     }
     else
     {
	tDep = iDiscrete3forces;
	fname = basename + ".x";
	npar  = 3*(npts+1);
     }
     bool byteswap;
     readSACheader( fname.c_str(), dt, t0, latsac, lonsac, cmpazsac, cmpincsac, utcsac, npts, byteswap );
     if( !useB )
	t0 = 0;
     
     if( geoCoordSet )
     {
	double laterr = fabs((latsac-lat)/lat);
	double lonerr = fabs((lonsac-lon)/lon);
	if( laterr > 1e-6 || lonerr > 1e-6 )
	{
	   if( proc_zero() )
	      cout << "WARNING in processSource: reading sac files: (lat,lon) location on sac file different from (lat,lon) on command line" << endl;
	}
     }
     freq = 1/dt;
     nipar = 1;
     par = new float_sw4 [npar];
     ipar = new int[1];
     ipar[0] = npts;
     size_t offset = 0;
     par[offset] = t0;
     if( tDep == iDiscrete6moments )
     {
	fname = basename + ".xx";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );
	if( timereverse )
	   revvector( npts, &par[offset+1]);
	offset += npts+1;
	par[offset] = t0;
	fname = basename + ".xy";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );     
	if( timereverse )
	   revvector( npts, &par[offset+1]);
	offset += npts+1;
	par[offset] = t0;
	fname = basename + ".xz";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );     
	if( timereverse )
	   revvector( npts, &par[offset+1]);
	offset += npts+1;
	par[offset] = t0;
	fname = basename + ".yy";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );     
	if( timereverse )
	   revvector( npts, &par[offset+1]);
	offset += npts+1;
	par[offset] = t0;
	fname = basename + ".yz";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );     
	if( timereverse )
	   revvector( npts, &par[offset+1]);
	offset += npts+1;
	par[offset] = t0;
	fname = basename + ".zz";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );     
	if( timereverse )
	   revvector( npts, &par[offset+1]);
     }
     else
     {
	fname = basename + ".x";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );
	if( timereverse )
	   revvector( npts, &par[offset+1]);
	offset += npts+1;
	par[offset] = t0;
	fname = basename + ".y";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );     
	if( timereverse )
	   revvector( npts, &par[offset+1]);
	offset += npts+1;
	par[offset] = t0;
	fname = basename + ".z";
	readSACdata( fname.c_str(), npts, &par[offset+1], byteswap );     
	if( timereverse )
	   revvector( npts, &par[offset+1]);
     }
  }

  // --------------------------------------------------------------------------- 
  // Regardless of how the location for the source was specified, we are going to
  // find the grid points associated with the location. (i.e., assign
  // i, j, k to valid values)
  // --------------------------------------------------------------------------- 
  if (geoCoordSet)
  {
    computeCartesianCoord(x, y, lon, lat);
    cartCoordSet = true;
    if( mVerbose >= 1 && proc_zero() )
    {
      printf("Cartesian coordinates of source at (lon, lat)=(%e, %e) is (x,y)=(%g, %g)\n", 
	     lon, lat, x, y);
    }
  }

  if (depthSet)
  {
    z = depth;
  }
  
  if (cartCoordSet)
  {
    float_sw4 xmin = 0.;
    float_sw4 ymin = 0.;
    float_sw4 zmin;

// only check the z>zmin when we have topography. For a flat free surface, we will remove sources too 
// close or above the surface in the call to mGlobalUniqueSources[i]->correct_Z_level()
    if (topographyExists()) // topography command must be read before the source command
      zmin = m_global_zmin;
    else
      zmin = 0;

    if ( (topographyExists() && (x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax )) ||
	 (!topographyExists() && (x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax || 
				  z < zmin || z > m_global_zmax)) )
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
      float_sw4 radconv = M_PI / 180.;
      float_sw4 S, D, R;
      strike -= mGeoAz; // subtract off the grid azimuth
      S = strike*radconv; D = dip*radconv; R = rake*radconv;
      
      mxx = -1.0 * ( sin(D) * cos(R) * sin (2*S) + sin(2*D) * sin(R) * sin(S)*sin(S) );
      myy =        ( sin(D) * cos(R) * sin (2*S) - sin(2*D) * sin(R) * cos(S)*cos(S) );
      mzz = -1.0 * ( mxx + myy );	
      mxy =        ( sin(D) * cos(R) * cos (2*S) + 0.5 * sin(2*D) * sin(R) * sin(2*S) );
      mxz = -1.0 * ( cos(D) * cos(R) * cos (S)   + cos(2*D) * sin(R) * sin(S) );
      myz = -1.0 * ( cos(D) * cos(R) * sin (S)   - cos(2*D) * sin(R) * cos(S) );
      //      if( m_myRank == 0 )
      //      {
      //	 cout << "Mxx = " << mxx << endl;
      //	 cout << "Myy = " << myy << endl;
      //	 cout << "Mzz = " << mzz << endl;
      //	 cout << "Mxy = " << mxy << endl;
      //	 cout << "Mxz = " << mxz << endl;
      //	 cout << "Myz = " << myz << endl;
      //      }
    }
  
  if (isMomentType)
  {
    // Remove amplitude variable
    mxx *= m0;
    mxy *= m0;
    mxz *= m0;
    myy *= m0;
    myz *= m0;
    mzz *= m0;
    //       m0 = 1;
    //       if( m_myRank == 0 )
    //       {
    //	  cout << "m0 = " << m0 << endl;
    //	  cout << "x  = " << x << endl;
    //	  cout << "y  = " << y << endl;
    //	  cout << "z  = " << z << endl;
    //	  cout << "t0 = " << t0 << endl;
    //	  cout << "freq=" << freq << endl;
    //       	  cout << "Mxx = " << mxx << endl;
    //      	  cout << "Myy = " << myy << endl;
    //      	  cout << "Mzz = " << mzz << endl;
    //      	  cout << "Mxy = " << mxy << endl;
    //      	  cout << "Mxz = " << mxz << endl;
    //      	  cout << "Myz = " << myz << endl;
    //       }
    // these have global location since they will be used by all processors
    sourcePtr = new Source(this, freq, t0, x, y, z, mxx, mxy, mxz, myy, myz, mzz,
			   tDep, formstring, topodepth, ncyc, par, npar, ipar, nipar, false ); // false is correctStrengthForMu
    if (sourcePtr->ignore())
    {
      delete sourcePtr;
    }
    else
    {
      a_GlobalUniqueSources[elocal].push_back(sourcePtr);
    }
      
  }
  else // point forcing
  {
    // Remove amplitude variable
    fx *= f0;
    fy *= f0;
    fz *= f0;
    //       f0 = 1;
    // global version (gets real coordinates)
    sourcePtr = new Source(this, freq, t0, x, y, z, fx, fy, fz, tDep, formstring, topodepth, ncyc,
			   par, npar, ipar, nipar, false ); // false is correctStrengthForMu
    //...and add it to the list of forcing terms
    if (sourcePtr->ignore())
    {
      delete sourcePtr;
    }
    else
    {
      a_GlobalUniqueSources[elocal].push_back(sourcePtr);
    }

  }	  
  if( npar > 0 )
     delete[] par;
  if( nipar > 0 )
     delete[] ipar;
  if (mVerbose >=4 && proc_zero())
    cout << "********Done parsing source command*********" << endl;
  }
}


//----------------------------------------------------------------------------
void EW::processRuptureHDF5(char* buffer, vector<vector<Source*> > & a_GlobalUniqueSources )
{
#ifdef USE_HDF5
  int event = 0;
  bool rfileset=false;
  char rfile[100];
  double stime, etime;
  stime = MPI_Wtime();

// bounding box
// only check the z>zmin when we have topography. For a flat free surface, we will remove sources too 
// close or above the surface in the call to mGlobalUniqueSources[i]->correct_Z_level()
  float_sw4 xmin = 0.;
  float_sw4 ymin = 0.;
  float_sw4 zmin;
  if (topographyExists()) // topography command must be read before the source command
    zmin = m_global_zmin;
  else
    zmin = -m_global_zmax;

  string err = "Rupture Error: ";

  char* token = strtok(buffer, " \t");
  REQUIRE2(strcmp("rupturehdf5", token) == 0, "ERROR: not a rupturehdf5 line...: " << token);
  token = strtok(NULL, " \t");

  while (token != NULL)
    {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
      if (startswith("file=",token))
      {
	token += 5; // read past 'file='
         strncpy(rfile, token,100);
	 rfileset = true;
      }
      else if(startswith("event=",token))
      {
	 token += 6;
	 //	 event = atoi(token);
	 //	 CHECK_INPUT( 0 <= event && event < m_nevent, err << "event no. "<< event << " out of range" );
	// Ignore if no events given
	 if( m_nevents_specified > 0 )
	 {
	    map<string,int>::iterator it = m_event_names.find(token);
	    CHECK_INPUT( it != m_event_names.end(), 
		     err << "event with name "<< token << " not found" );
	    event = it->second;
	 }
      }
      else
      {
         badOption("rupturehdf5", token);
      }
      token = strtok(NULL, " \t");
    }


  if( event_is_in_proc(event) )
  {
     event = global_to_local_event(event);
  if( rfileset)
    readRuptureHDF5(rfile, a_GlobalUniqueSources, this, event, m_global_xmax, m_global_ymax, m_global_zmax, mGeoAz, xmin, ymin, zmin, mVerbose, m_nwriters);
  }

  etime = MPI_Wtime();
  
  if (proc_zero())
      cout << "Process rupture data, took " << etime-stime << "seconds." << endl;
#else
  if (proc_zero())
    cout << "Using HDF5 rupture input but sw4 is not compiled with HDF5!"<< endl;
#endif

} // end processRupture()


//----------------------------------------------------------------------------
void EW::processRupture(char* buffer, vector<vector<Source*> > & a_GlobalUniqueSources )
{
// the rupture command reads a file with
// point moment tensor sources in the strike, dip, rake format
// for each source, the slip velocity time function is defined by a discrete time function
  Source* sourcePtr;
  double stime, etime;
  stime = MPI_Wtime();
  
  float_sw4 m0 = 1.0;
  float_sw4 t0=0.0, f0=1.0, freq=1.0;
  // Should be center of the grid
  double x = 0.0, y = 0.0, z = 0.0;
  int i = 0, j = 0, k = 0;
  float_sw4 mxx=0.0, mxy=0.0, mxz=0.0, myy=0.0, myz=0.0, mzz=0.0;
  float_sw4 strike=0.0, dip=0.0, rake=0.0;
  float_sw4 fx=0.0, fy=0.0, fz=0.0;
  int event = 0;
  
  double lat = 0.0, lon = 0.0;
  bool topodepth = true;
  
  bool rfileset=false;

  timeDep tDep = iDiscrete;
  char formstring[100];
  strcpy(formstring, "Discrete");
  char rfile[100];

// bounding box
// only check the z>zmin when we have topography. For a flat free surface, we will remove sources too 
// close or above the surface in the call to mGlobalUniqueSources[i]->correct_Z_level()
  float_sw4 xmin = 0.;
  float_sw4 ymin = 0.;
  float_sw4 zmin;
  if (topographyExists()) // topography command must be read before the source command
    zmin = m_global_zmin;
  else
    zmin = -m_global_zmax;

  string err = "Rupture Error: ";

  char* token = strtok(buffer, " \t");
  REQUIRE2(strcmp("rupture", token) == 0, "ERROR: not a rupture line...: " << token);
  token = strtok(NULL, " \t");

  while (token != NULL)
    {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
      if (startswith("file=",token))
      {
	token += 5; // read past 'file='
         strncpy(rfile, token,100);
	 rfileset = true;
      }
      else if(startswith("event=",token))
      {
	 token += 6;
	 //	 event = atoi(token);
	 //	 CHECK_INPUT( 0 <= event && event < m_nevent, err << "event no. "<< event << " out of range" );
	// Ignore if no events given
	 if( m_nevents_specified > 0 )
	 {
	    map<string,int>::iterator it = m_event_names.find(token);
            //	    CHECK_INPUT( it != m_event_names.end(), 
            //		     err << "event with name "<< token << " not found" );
            //	    event = it->second;
             if( it != m_event_names.end() )
                event = it->second;
             else if( proc_zero() )
                std::cout << "Rupture warning: event with name " << token << " not found" << std::endl;
	 }
      }
      else
      {
         badOption("rupture", token);
      }
      token = strtok(NULL, " \t");
    }

  float_sw4 rVersion;

  const int bufsize=1024;
  char buf[bufsize];
  
// Discrete source time function
  float_sw4* par=NULL;
  int* ipar=NULL;
  int npar=0, nipar=0, ncyc=0;
  if( rfileset )
  {
     //  g(t) defined by spline points on a uniform grid, read from file.
     //  Format: t0, dt, npts
     //          g_1
     //          g_2
     //         ....

    FILE* fd=fopen(rfile, "r" );
    CHECK_INPUT( fd !=NULL , err << "Rupture file " << rfile << " not found" );
    if (proc_zero())
      printf("Opened rupture file '%s'\n", rfile);
// read 1st line
    fgets(buf,bufsize,fd);
    sscanf(buf," %lg", &rVersion );
    if (proc_zero())
      printf("Version = %.1f\n", rVersion);
// read 2nd line, starting header block
    fgets(buf,bufsize,fd);
    char* token = strtok(buf, " \t");
//    printf("token: '%s'\n", token);
    REQUIRE2(strcmp("PLANE", token) == 0, "ERROR: not a HEADER BLOCK line...: " << token);
// read the number of planes
    int nseg;
    token = strtok(NULL, " \t");
    nseg = atoi(token);
    if (proc_zero())
      printf("Number of segments in header block: %i\n", nseg);
// read each header block
    for (int seg=0; seg<nseg; seg++)
    {
      double elon, elat, len, wid, stk, dip, dtop, shyp, dhyp;
      int nstk, ndip;
      fgets(buf,bufsize,fd);
      sscanf(buf,"%lg %lg %i %i %lg %lg", &elon, &elat, &nstk, &ndip, &len, &wid);
      fgets(buf,bufsize,fd);
      sscanf(buf,"%lg %lg %lg %lg %lg", &stk, &dip, &dtop, &shyp, &dhyp);
      if (proc_zero())
      {
	printf("Seg #%i: elon=%g, elat=%g, nstk=%i, ndip=%i, len=%g, wid=%g\n", 
	       seg+1, elon, elat, nstk, ndip, len, wid);
	printf("        stk=%g, dip=%g, dtop=%g, shyp=%g, dhyp=%g\n", stk, dip, dtop, shyp, dhyp);
      }
      
    }
// read header for data block
    fgets(buf,bufsize,fd);
    token = strtok(buf, " \t");
//    printf("token: '%s'\n", token);
    REQUIRE2(strcmp("POINTS", token) == 0, "ERROR: not a DATA BLOCK line...: " << token);
// read the number of points
    int npts;
    token = strtok(NULL, " \t");
    npts = atoi(token);
    if (proc_zero())
      printf("Number of point sources in data block: %i\n", npts);

// read all point sources
    int nSources=0, nu1=0, nu2=0, nu3=0;
    for (int pts=0; pts<npts; pts++) 
    {
      double lon, lat, dep, stk, dip, area, tinit, dt, rake, slip1, slip2, slip3;
      int nt1=0, nt2=0, nt3=0;
      fgets(buf,bufsize,fd);
      sscanf(buf,"%lg %lg %lg %lg %lg %lg %lg %lg", &lon, &lat, &dep, &stk, &dip, &area, 
	     &tinit, &dt);
      fgets(buf,bufsize,fd);
      sscanf(buf,"%lg %lg %i %lg %i %lg %i", &rake, &slip1, &nt1, &slip2, &nt2, &slip3, &nt3);
// nothing to do if nt1=nt2=nt3=0
      if (nt1<=0 && nt2<=0 && nt3<=0) continue;
      if (proc_zero() && mVerbose >= 2)
      {
	printf("point #%i: lon=%g, lat=%g, dep=%g, stk=%g, dip=%g, area=%g, tinit=%g, dt=%g\n", 
	       pts+1, lon, lat, dep, stk, dip, area, tinit, dt);
	printf("          rake=%g, slip1=%g, nt1=%i, slip2=%g, nt2=%i, slip3=%g, nt3=%i\n", 
	       rake, slip1, nt1, slip2, nt2, slip3, nt3);
      }
      
// read discrete time series for u1
      if (nt1>0)
      {
	nu1++;
// note that the first data point is always zero, but the last is not
// for this reason we always pad the time zeries with a '0' 
// also note that we need at least 7 data points, i.e. nt1>=6
	int nt1dim = max(6,nt1);
	par = new float_sw4[nt1dim+2];
	par[0]  = tinit;
	t0      = tinit;
	freq    = 1/dt;
	ipar    = new int[1];
	ipar[0] = nt1dim+1; // add an extra point 
	fgets(buf,bufsize,fd);
	token = strtok(buf, " \t");
//	printf("buf='%s'\n", buf);
	for( int i=0 ; i < nt1 ; i++ )
	{
// read another line if there are no more tokens
	  if (token == NULL)
	  {
	    fgets(buf,bufsize,fd);
	    token = strtok(buf, " \t");
	  }
//	  printf("token='%s'\n", token);
	  sscanf(token,"%lg", &par[i+1] );
// read next token
	  token = strtok(NULL, " \t");
	}
// pad with 0
	if (nt1 < 6)
	{
	  for (int j=nt1; j<6; j++)
	    par[j+1]=0.;
	}
	
// last 0
	par[nt1dim+1]= 0.0;

// scale cm/s to m/s
	for (int i=1; i<=nt1dim+1; i++)
	{
	  par[i] *= 1e-2;
	}

// AP: Mar. 1, 2016: Additional scaling is needed to make the integral of the time function = 1
        float_sw4 slip_m=slip1*1e-2;
        float_sw4 slip_sum=0;
	for (int i=1; i<=nt1dim+1; i++)
	{
	  slip_sum += par[i];
	}
        slip_sum *=dt;

        if (proc_zero() && mVerbose >= 2)
        {
           printf("INFO: SRF file: dt*sum(slip_vel)=%e [m], total slip (from header)=%e [m]\n", slip_sum, slip_m);
        }
// scale time series to sum to integrate to one        
	for (int i=1; i<=nt1dim+1; i++)
	{
           par[i] /= slip_sum;
	}
        if (proc_zero() && mVerbose >= 2)
        {
           slip_sum=0;
           for (int i=1; i<=nt1dim+1; i++)
           {
              slip_sum += par[i];
           }
           slip_sum *=dt;
           printf("INFO: SRF file: After scaling time series: dt*sum(par)=%e [m]\n", slip_sum);
        }
//done scaling        
        
	npar = nt1dim+2;
	nipar = 1;

        // printf("Read discrete time series: tinit=%g, dt=%g, nt1=%i\n", tinit, dt, nt1);
	// for (int i=0; i<nt1+1; i++)
	//   printf("Sv1[%i]=%g\n", i+1, par[i+1]);

// convert lat, lon, depth to (x,y,z)
	computeCartesianCoord(x, y, lon, lat);
// convert depth in [km] to [m]
	z = dep * 1e3;

// convert strike, dip, rake to Mij
	float_sw4 radconv = M_PI / 180.;
	float_sw4 S, D, R;
	stk -= mGeoAz; // subtract off the grid azimuth
	S = stk*radconv; D = dip*radconv; R = rake*radconv;
      
	mxx = -1.0 * ( sin(D) * cos(R) * sin (2*S) + sin(2*D) * sin(R) * sin(S)*sin(S) );
	myy =        ( sin(D) * cos(R) * sin (2*S) - sin(2*D) * sin(R) * cos(S)*cos(S) );
	mzz = -1.0 * ( mxx + myy );	
	mxy =        ( sin(D) * cos(R) * cos (2*S) + 0.5 * sin(2*D) * sin(R) * sin(2*S) );
	mxz = -1.0 * ( cos(D) * cos(R) * cos (S)   + cos(2*D) * sin(R) * sin(S) );
	myz = -1.0 * ( cos(D) * cos(R) * sin (S)   - cos(2*D) * sin(R) * cos(S) );

// scale (note that the shear modulus is not yet available. Also note that we convert [cm] to [m])
	m0 = area*1e-4 * slip1*1e-2;
      
	mxx *= m0;
	mxy *= m0;
	mxz *= m0;
	myy *= m0;
	myz *= m0;
	mzz *= m0;

// before creating the source, make sure (x,y,z) is inside the computational domain

// only check the z>zmin when we have topography. For a flat free surface, we will remove sources too 
// close or above the surface in the call to mGlobalUniqueSources[i]->correct_Z_level()

	if (x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax || z < zmin || z > m_global_zmax)
	{
	  stringstream sourceposerr;
	  sourceposerr << endl
		       << "***************************************************" << endl
		       << " ERROR:  Source positioned outside grid!  " << endl
		       << endl
		       << " Source from rupture file @" << endl
		       << "  x=" << x << " y=" << y << " z=" << z << endl 
		       << "  lat=" << lat << " lon=" << lon << " dep=" << dep << endl 
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
	}
	else if( event_is_in_proc(event) )
	{
           event = global_to_local_event(event);
           sourcePtr = new Source(this, freq, t0, x, y, z, mxx, mxy, mxz, myy, myz, mzz,
				 tDep, formstring, topodepth, ncyc, par, npar, ipar, nipar, true ); // true is correctStrengthForMu

	  if (sourcePtr->ignore())
	  {
	    delete sourcePtr;
	  }
	  else
	  {
	    a_GlobalUniqueSources[event].push_back(sourcePtr);
	    nSources++;
	  }
	}
	

// deallocate temporary arrays...
	delete[] par;
	delete[] ipar;

      } // end if nt1 >0

// read past discrete time series for u2
      if (nt2>0)
      {
	nu2++;
	double dum;
	if (proc_zero())
	  printf("WARNING nt2=%i > 0 will be ignored\n", nt2);
	fgets(buf,bufsize,fd);
	token = strtok(buf, " \t");
//	printf("buf='%s'\n", buf);
	for( int i=0 ; i < nt2 ; i++ )
	{
// read another line if there are no more tokens
	  if (token == NULL)
	  {
	    fgets(buf,bufsize,fd);
	    token = strtok(buf, " \t");
	  }
//	  printf("token='%s'\n", token);
	  sscanf(token,"%lg", &dum );
// read next token
	  token = strtok(NULL, " \t");
	}
      } // end if nt2 > 0

// read past discrete time series for u3
      if (nt3>0)
      {
	nu3++;
	double dum;
	if (proc_zero())
	  printf("WARNING nt3=%i > 0 will be ignored\n", nt3);
	fgets(buf,bufsize,fd);
	token = strtok(buf, " \t");
//	printf("buf='%s'\n", buf);
	for( int i=0 ; i < nt3 ; i++ )
	{
// read another line if there are no more tokens
	  if (token == NULL)
	  {
	    fgets(buf,bufsize,fd);
	    token = strtok(buf, " \t");
	  }
//	  printf("token='%s'\n", token);
	  sscanf(token,"%lg", &dum );
// read next token
	  token = strtok(NULL, " \t");
	}
      } // end if nt3 > 0
      
    } // end for all sources
    if (proc_zero())
      printf("Read npts=%i, made %i point moment tensor sources, nu1=%i, nu2=%i, nu3=%i\n", 
	     npts, nSources, nu1, nu2, nu3);
    
    fclose(fd);
  }

  etime = MPI_Wtime();
  if (proc_zero())
      cout << "Process rupture data, took " << etime-stime << "seconds." << endl;
} // end processRupture()


//------------------------------------------------------------------------
void EW::processMaterialBlock( char* buffer, int & blockCount )
{
  float_sw4 vpgrad=0.0, vsgrad=0.0, rhograd=0.0;
  bool x1set=false, x2set=false, y1set=false, y2set=false, 
    z1set=false, z2set=false;

  float_sw4 x1=0.0, x2=0.0, y1=0.0, y2=0.0, z1=0.0, z2=0.0;
  int i1=-1, i2=-1, j1=-1, j2=-1, k1=-1, k2=-1;

  string name = "Block";

  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("block", token) == 0,
	      "ERROR: material block can be set by a block line, not: " << token);

  string err = token;
  err += " Error: ";

  token = strtok(NULL, " \t");

  float_sw4 vp=-1, vs=-1, rho=-1, qp=-1, qs=-1, freq=1;
  bool absDepth=false;

  while (token != NULL)
    {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
	break;
// the xygrad keywords must occur before the corresponding xy keywords
      if (startswith("rhograd=", token))
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
      else if (startswith("vp=", token) )
      {
         token += 3; // skip vp=
         vp = atof(token);
      }
      else if (startswith("vs=", token) )
      {
         token += 3; // skip vs=
         vs = atof(token);
      }
      else if (startswith("rho=", token))
      {
         token += 4; // skip rho=
         rho = atof(token);
      }
      else if (startswith("r=", token)) // superseded by rho=, but keep for backward compatibility
      {
         token += 2; // skip r=
         rho = atof(token);
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

  // Set up a block on the EW object.

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
     //    z1 = m_global_zmin - (m_global_zmax-m_global_zmin);
     z1 = m_global_zmin - 1e10;

  if (z2set)
  {
    CHECK_INPUT(topographyExists() || z2 >= 0.,
            err << "z2 is less than the minimum z, " << z2 << " < " << 0.);
    // CHECK_INPUT(z2 <= m_global_zmax,
    // 		err << "z2 is greater than the maximum z, " << z2 << " > " << m_global_zmax);
  }
  else
     //    z2 = m_global_zmax + (m_global_zmax-m_global_zmin);
     z2 = m_global_zmax + 1e10;

  CHECK_INPUT( z2 >= z1, " (z1..z2), upper bound is smaller than lower bound");

  if(getVerbosity() >=2 &&  m_myRank == 0 )
     cout << name << " has bounds " << x1 << " " << x2 << " " << y1 << " "
	  << y2 << " " << z1 << " " << z2 << endl;

  CHECK_INPUT( vs > 0 && vp > 0 && rho > 0 , "Error in block " << name << " vp vs rho are   "
	       << vp << " " << vs << " " << rho );

  MaterialBlock* bl = new MaterialBlock( this ,rho, vs, vp, x1, x2, y1, y2, z1, z2, qs, qp, freq );
  bl->set_gradients( rhograd, vsgrad, vpgrad );
  bl->set_absoluteDepth( absDepth );
  add_mtrl_block( bl );
}

//-----------------------------------------------------------------------
void EW::processAnisotropicMaterialBlock( char* buffer,  int & blockCount )
{
   float_sw4 rho=-1, rhograd=0.0;
   float_sw4 c[21], cgrad[21];
   for( int m=0 ; m < 21 ; m++ )
   {
      c[m] = -1;
      cgrad[m] = 0;
   }

   bool x1set=false, x2set=false, y1set=false, y2set=false, 
      z1set=false, z2set=false;

   float_sw4 x1=0.0, x2=0.0, y1=0.0, y2=0.0, z1=0.0, z2=0.0;
   int i1=-1, i2=-1, j1=-1, j2=-1, k1=-1, k2=-1;


   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("ablock", token) == 0,
	       "ERROR: material block can be set by a ablock line, not: " << token);

   string err = token;
   err += " Error: ";

   token = strtok(NULL, " \t");
   bool absDepth=false;

   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
	 break;
// the xygrad keywords must occur before the corresponding xy keywords
      if (startswith("rhograd=", token))
      {
         token += 8; // skip rhograd=
         rhograd = atof(token);
      }
      if (startswith("rho=", token))
      {
         token += 4; // skip rho=
         rho = atof(token);
      }
      else if (startswith("c11grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[0] = atof(token);
      }
      else if (startswith("c11=", token))
      {
         token += 4; // skip c1=
         c[0] = atof(token);
      }
      else if (startswith("c16grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[1] = atof(token);
      }
      else if (startswith("c16=", token))
      {
         token += 4; // skip c1=
         c[1] = atof(token);
      }
      else if (startswith("c15grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[2] = atof(token);
      }
      else if (startswith("c15=", token))
      {
         token += 4; // skip c1=
         c[2] = atof(token);
      }
      else if (startswith("c12grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[3] = atof(token);
      }
      else if (startswith("c12=", token))
      {
         token += 4; // skip c1=
         c[3] = atof(token);
      }
      else if (startswith("c14grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[4] = atof(token);
      }
      else if (startswith("c14=", token))
      {
         token += 4; // skip c1=
         c[4] = atof(token);
      }
      else if (startswith("c13grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[5] = atof(token);
      }
      else if (startswith("c13=", token))
      {
         token += 4; // skip c1=
         c[5] = atof(token);
      }
      else if (startswith("c66grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[6] = atof(token);
      }
      else if (startswith("c66=", token))
      {
         token += 4; // skip c1=
         c[6] = atof(token);
      }
      else if (startswith("c56grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[7] = atof(token);
      }
      else if (startswith("c56=", token))
      {
         token += 4; // skip c1=
         c[7] = atof(token);
      }
      else if (startswith("c26grad=", token))
      {
         token += 8; // skip c1grad=
         cgrad[8] = atof(token);
      }
      else if (startswith("c26=", token))
      {
         token += 4; // skip c1=
         c[8] = atof(token);
      }
      else if (startswith("c46grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[9] = atof(token);
      }
      else if (startswith("c46=", token))
      {
         token += 4; // skip cxx=
         c[9] = atof(token);
      }
      else if (startswith("c36grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[10] = atof(token);
      }
      else if (startswith("c36=", token))
      {
         token += 4; // skip cxx=
         c[10] = atof(token);
      }
      else if (startswith("c55grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[11] = atof(token);
      }
      else if (startswith("c55=", token))
      {
         token += 4; // skip cxx=
         c[11] = atof(token);
      }
      else if (startswith("c25grad=", token))
      {
         token += 8; // skip c10grad=
         cgrad[12] = atof(token);
      }
      else if (startswith("c25=", token))
      {
         token += 4; // skip c10=
         c[12] = atof(token);
      }
      else if (startswith("c45grad=", token))
      {
         token += 8; // skip c10grad=
         cgrad[13] = atof(token);
      }
      else if (startswith("c45=", token))
      {
         token += 4; // skip c10=
         c[13] = atof(token);
      }
      else if (startswith("c35grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[14] = atof(token);
      }
      else if (startswith("c35=", token))
      {
         token += 4; // skip cxx=
         c[14] = atof(token);
      }
      else if (startswith("c22grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[15] = atof(token);
      }
      else if (startswith("c22=", token))
      {
         token += 4; // skip cxx=
         c[15] = atof(token);
      }
      else if (startswith("c24grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[16] = atof(token);
      }
      else if (startswith("c24=", token))
      {
         token += 4; // skip cxx=
         c[16] = atof(token);
      }
      else if (startswith("c23grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[17] = atof(token);
      }
      else if (startswith("c23=", token))
      {
         token += 4; // skip cxx=
         c[17] = atof(token);
      }
      else if (startswith("c44grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[18] = atof(token);
      }
      else if (startswith("c44=", token))
      {
         token += 4; // skip cxx=
         c[18] = atof(token);
      }
      else if (startswith("c34grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[19] = atof(token);
      }
      else if (startswith("c34=", token))
      {
         token += 4; // skip cxx=
         c[19] = atof(token);
      }
      else if (startswith("c33grad=", token))
      {
         token += 8; // skip cxxgrad=
         cgrad[20] = atof(token);
      }
      else if (startswith("c33=", token))
      {
         token += 4; // skip cxx=
         c[20] = atof(token);
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
      else if (startswith("absdepth=", token) )
      {
	token += 9; // skip absdepth=
	absDepth = (bool) atoi(token);
      }
      else
      {
         badOption("ablock",token);
      }
      token = strtok(NULL, " \t");
   }
   string name = "ABlock";
   blockCount++;
   stringstream blockname;
   blockname << name << " " << blockCount;
   name = blockname.str();

  // Set up a block on the EW object.

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
      z1 = m_global_zmin - (m_global_zmax-m_global_zmin);

   if (z2set)
   {
      CHECK_INPUT(topographyExists() || z2 >= 0.,
            err << "z2 is less than the minimum z, " << z2 << " < " << 0.);
    // CHECK_INPUT(z2 <= m_global_zmax,
    // 		err << "z2 is greater than the maximum z, " << z2 << " > " << m_global_zmax);
   }
   else
      z2 = m_global_zmax + (m_global_zmax-m_global_zmin);

   CHECK_INPUT( z2 >= z1, " (z1..z2), upper bound is smaller than lower bound");

   if(getVerbosity() >=2 &&  m_myRank == 0 )
      cout << name << " has bounds " << x1 << " " << x2 << " " << y1 << " "
	   << y2 << " " << z1 << " " << z2 << endl;

   CHECK_INPUT( rho > 0 , "Error in ablock " << name << " rho is " << rho );

   AnisotropicMaterialBlock* bl = new AnisotropicMaterialBlock( this, rho, c, x1, x2, y1, y2, z1, z2 );

   bl->set_gradients( rhograd, cgrad );
   bl->set_absoluteDepth( absDepth );
   m_anisotropic_mtrlblocks.push_back(bl);
}   

//-----------------------------------------------------------------------
void EW::processReceiverHDF5(char* buffer, vector<vector<TimeSeries*> > & a_GlobalTimeSeries)
{
  string inFileName = "station";
  string fileName   = "station_out";
  string staName    = "station";
  int writeEvery    = 1000;
  int downSample    = 1;
  int event         = 0;
  TimeSeries::receiverMode mode=TimeSeries::Displacement;
  double stime, etime;
  stime = MPI_Wtime();

  char* token = strtok(buffer, " \t");

  CHECK_INPUT(strcmp("rechdf5", token) == 0 || strcmp("sachdf5", token) == 0, "ERROR: not a rechdf5 line...: " << token);
  token = strtok(NULL, " \t");

  string err = "RECEIVER Error: ";

  while (token != NULL)
  {
     if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
        break;

     if(startswith("infile=", token))
     {
        token += 7; // skip infile=
        inFileName= token;
     }
     else if(startswith("outfile=", token))
     {
        token += 8; // skip outfile=
        fileName = token;
     }
     else if (startswith("writeEvery=", token))
     {
       token += strlen("writeEvery=");
       writeEvery = atoi(token);
       CHECK_INPUT(writeEvery >= 0,
	       err << "rechdf5 command: writeEvery must be set to a non-negative integer, not: " << token);
     }
     else if (startswith("downSample=", token) || startswith("downsample=", token))
     {
       token += strlen("downsample=");
       downSample = atoi(token);
       CHECK_INPUT(downSample >= 1,
	       err << "rechdf5 command: downsample must be set to an integer greater or equal than 1, not: " << token);
     }
     else if(startswith("event=",token))
     {
	token += 6;
	// Ignore if no events given
	if( m_nevents_specified > 0 )
	{
	   map<string,int>::iterator it = m_event_names.find(token);
           //	   CHECK_INPUT( it != m_event_names.end(), 
           //		     err << "event with name "<< token << " not found" );
           //	   event = it->second;
             if( it != m_event_names.end() )
                event = it->second;
             else if( proc_zero() )
                std::cout << "Receiver warning: event with name " << token << " not found" << std::endl;
	}
     }
     else if( startswith("variables=", token) )
     {
       token += strlen("variables=");

       if( strcmp("displacement",token)==0 )
	 mode = TimeSeries::Displacement;
       else if( strcmp("velocity",token)==0 )
	 mode = TimeSeries::Velocity;
       else if( strcmp("div",token)==0 )
	 mode = TimeSeries::Div;
       else if( strcmp("curl",token)==0 )
	 mode = TimeSeries::Curl;
       else if( strcmp("strains",token)==0 )
	 mode = TimeSeries::Strains;
       else if( strcmp("displacementgradient",token)==0 )
	 mode = TimeSeries::DisplacementGradient;
       else
       {
	 if (proc_zero())
	   cout << "receiver command: variables=" << token << " not understood" << endl
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

#ifdef USE_HDF5
  bool is_obs = false;
  if (writeEvery % downSample != 0) {
    writeEvery = (int)writeEvery / downSample;
    writeEvery *= downSample;
    if (proc_zero())
      cout << "receiver command: writeEvery=" << writeEvery << " is not a multiple of downsample, "
          << downSample << "adjustding writeEvery to " << writeEvery << endl;
  }
  event = global_to_local_event( event );
  if( event >= 0 )
     readStationHDF5(this, inFileName, fileName, writeEvery, downSample, mode, event, &a_GlobalTimeSeries, m_global_xmax, m_global_ymax, is_obs, false, false, 0, 0, false, false, false, 0, false, 0);
#else
  if (proc_zero())
    cout << "Using HDF5 station input but sw4 is not compiled with HDF5!"<< endl;
#endif
  etime = MPI_Wtime();
  if (a_GlobalTimeSeries.size() > 0 && a_GlobalTimeSeries[0].size() > 0) 
    a_GlobalTimeSeries[0][0]->addReadTime(etime-stime);
}

//-----------------------------------------------------------------------
void EW::processReceiver(char* buffer, vector<vector<TimeSeries*> > & a_GlobalTimeSeries)
{
  double x=0.0, y=0.0, z=0.0;
  double lat = 0.0, lon = 0.0, depth = 0.0;
  bool cartCoordSet = false, geoCoordSet = false;
  string fileName = "station";
  string hdf5FileName = "station.hdf5";
  string staName = "station";
  bool staNameGiven=false;
  double stime, etime;
  stime = MPI_Wtime();
  
  int writeEvery = 1000;
  int downSample = 1;

  bool topodepth = false;

  bool usgsformat = 0, sacformat = 1, hdf5format = 0; // default is to write sac files
  TimeSeries::receiverMode mode=TimeSeries::Displacement;

  char* token = strtok(buffer, " \t");
  bool nsew=false; 
  int event=0;
  //int vel=0;

// tmp
//  cerr << "******************** INSIDE process receiver *********************" << endl;

  CHECK_INPUT(strcmp("rec", token) == 0 || strcmp("sac", token) == 0, "ERROR: not a rec line...: " << token);
  token = strtok(NULL, " \t");

  string err = "RECEIVER Error: ";

//* testing
  // if (proc_zero())
  //   cout << "start parsing of receiver command, token:" << token << "(end token)" << endl;

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
       depth = z = atof(token);
       topodepth = false; // absolute depth (below mean sea level)
       CHECK_INPUT(z <= m_global_zmax,
		   "receiver command: z must be less than or equal to zmax, not " << z);
     }
     else if (startswith("depth=", token))
     {
        token += 6; // skip depth=
       z = depth = atof(token);
       topodepth = true; // by depth we here mean depth below topography
       CHECK_INPUT(depth >= 0.0,
	       err << "receiver command: depth must be greater than or equal to zero");
       CHECK_INPUT(depth <= m_global_zmax,
		   "receiver command: depth must be less than or equal to zmax, not " << depth);
     }
//                        1234567890
     else if (startswith("topodepth=", token))
     {
        token += 10; // skip topodepth=
       z = depth = atof(token);
       topodepth = true; // by depth we here mean depth below topography
       CHECK_INPUT(depth >= 0.0,
	       err << "receiver command: depth must be greater than or equal to zero");
       CHECK_INPUT(depth <= m_global_zmax,
		   "receiver command: depth must be less than or equal to zmax, not " << depth);
     }
     else if(startswith("hdf5file=", token))
     {
        token += 9; // skip file=
        hdf5FileName = token;
     }
     else if(startswith("file=", token))
     {
        token += 5; // skip file=
        fileName = token;
     }
     else if (startswith("sta=", token))
     {
        token += strlen("sta=");
        staName = token;
	staNameGiven=true;
     }
     else if( startswith("nsew=", token) )
     {
        token += strlen("nsew=");
        nsew = atoi(token) == 1;
     }
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
     else if (startswith("downSample=", token) || startswith("downsample=", token))
     {
       token += strlen("downSample=");
       downSample = atoi(token);
       CHECK_INPUT(downSample >= 1,
	       err << "sac command: downSample must be set to an integer greater or equal than 1, not: " << token);
     }
     else if(startswith("event=",token))
     {
	token += 6;
	//	event = atoi(token);
	//	CHECK_INPUT( 0 <= event && event < m_nevent, err << "event no. "<< event << " out of range" );
	// Ignore if no events given
	if( m_nevents_specified > 0 )
	{
	   map<string,int>::iterator it = m_event_names.find(token);
           //	   CHECK_INPUT( it != m_event_names.end(), 
           //		     err << "event with name "<< token << " not found" );
           //	   event = it->second;
             if( it != m_event_names.end() )
                event = it->second;
             else if( proc_zero() )
                std::cout << "Receiver warning: event with name " << token << " not found" << std::endl;
	}
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
     else if( startswith("hdf5format=", token) )
     {
        token += strlen("hdf5format=");
        hdf5format = atoi(token);
     }
     else if( startswith("variables=", token) )
     {
//* testing
       // if (proc_zero())
       // 	 printf("Inside rec command, before parsing 'variables=', token:'%s'(end token)\n", token);
       
       token += strlen("variables=");

//* testing
       // if (proc_zero())
       // 	 printf("Inside rec command, after parsing 'variables=', token:'%s'(end token)\n", token);

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
       else if( strcmp("displacementgradient",token)==0 )
       {
	 mode = TimeSeries::DisplacementGradient;
       }
       else
       {
	 if (proc_zero())
	   cout << "receiver command: variables=" << token << " not understood" << endl
		<< "using default mode (displacement)" << endl << endl;
	 mode = TimeSeries::Displacement;
       }
       
     }
     else
     {
        badOption("receiver", token);
     }
     token = strtok(NULL, " \t");
//* testing
     // if (proc_zero())
     //   cout << "rec command: Bottom of while loop, token:" << token << "(end token)" << endl;
     
  }  
  //  cout << "end receiver " << m_myRank << endl;

  if (geoCoordSet)
  {
    computeCartesianCoord(x, y, lon, lat);
// check if (x,y) is within the computational domain
  }

  if (!staNameGiven)
    staName = fileName;

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
      receivererr << " No RECEIVER file will be generated for file = " << fileName << endl;
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
  else if( event_is_in_proc(event) )
  {
    if (writeEvery % downSample != 0) {
      writeEvery = (int)writeEvery / downSample;
      writeEvery *= downSample;
      if (proc_zero())
        cout << "receiver command: writeEvery=" << writeEvery << " is not a multiple of downsample, "
            << downSample << "adjustding writeEvery to " << writeEvery << endl;
    }

    event = global_to_local_event(event);
    TimeSeries *ts_ptr = new TimeSeries(this, fileName, staName, mode, sacformat, usgsformat, hdf5format, hdf5FileName, x, y, depth, 
					topodepth, writeEvery, downSample, !nsew, event );
#if USE_HDF5
    if (hdf5format) {
      if(a_GlobalTimeSeries[event].size() == 0) {
        ts_ptr->allocFid();
        ts_ptr->setTS0Ptr(ts_ptr);
      }
      else {
        ts_ptr->setFidPtr(a_GlobalTimeSeries[event][0]->getFidPtr());
        ts_ptr->setTS0Ptr(a_GlobalTimeSeries[event][0]);
      }
    }
#endif

// include the receiver in the global list
    a_GlobalTimeSeries[event].push_back(ts_ptr);
  }

  etime = MPI_Wtime();
  if (a_GlobalTimeSeries.size() > 0 && a_GlobalTimeSeries[0].size() > 0) 
    a_GlobalTimeSeries[0][0]->addReadTime(etime-stime);
  
}

//-----------------------------------------------------------------------
void EW::processObservationHDF5( char* buffer, vector<vector<TimeSeries*> > & a_GlobalTimeSeries)
{
  double x=0.0, y=0.0, z=0.0;
  double lat = 0.0, lon = 0.0, depth = 0.0;
  float_sw4 t0 = 0;
  float_sw4 scalefactor=1;
  bool cartCoordSet = false, geoCoordSet = false;
  string inhdf5file = "";
  string outhdf5file = "station";
  int writeEvery = 0;
  int downSample = 1;
  TimeSeries::receiverMode mode=TimeSeries::Displacement;
  float_sw4 winl, winr;
  bool winlset=false, winrset=false;
  char exclstr[4]={'\0','\0','\0','\0'};
  bool usex=true, usey=true, usez=true;
  bool scalefactor_set=false;
  int event=0;

  char* token = strtok(buffer, " \t");
  m_filter_observations = true;

  CHECK_INPUT(strcmp("observationhdf5", token) == 0 || strcmp("obshdf5", token) == 0, "ERROR: not an observation line...: " << token);
  token = strtok(NULL, " \t");

  string err = "OBSERVATION Error: ";

  while (token != NULL)
  {
     if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
        break;
     else if(startswith("hdf5file=", token))
     {
        token += 9; // skip hdf5file=
        inhdf5file = token;
     }
     else if(startswith("outhdf5file=", token))
     {
        token += 12; // skip outhdf5file=
        outhdf5file = token;
     }
     else if(startswith("event=",token))
     {
	token += 6;
	//	event = atoi(token);
	//	CHECK_INPUT( 0 <= event && event < m_nevent, err << "event no. "<< event << " out of range" );
	// Ignore if no events given
	if( m_nevents_specified > 0 )
	{
	   map<string,int>::iterator it = m_event_names.find(token);
           //	   CHECK_INPUT( it != m_event_names.end(), 
           //		     err << "event with name "<< token << " not found" );
           //	   event = it->second;
             if( it != m_event_names.end() )
                event = it->second;
             else if( proc_zero() )
                std::cout << "Observation warning: event with name " << token << " not found" << std::endl;
	}
     }
     // (small) shifts of the observation in time can be used to compensate for incorrect velocites
     // in the material model
     else if(startswith("shift=", token))
     {
        token += 6; // skip shift=
        t0 = atof(token);
     }
     //     else if( startswith("utc=",token) )
     //     {
     //	token += 4;
     //        if( strcmp("ignore",token)==0 || strcmp("off",token)==0 )
     //	   ignore_utc = true;
     //	else
     //	{
     //	   int year,month,day,hour,minute,second,msecond, fail;
     //	   // Format: 01/04/2012:17:34:45.2343  (Month/Day/Year:Hour:Min:Sec.fraction)
     //	   parsedate( token, year, month, day, hour, minute, second, msecond, fail );
     //	   if( fail == 0 )
     //	   {
     //              utcset = true;
     //	      utc[0] = year;
     //	      utc[1] = month;
     //	      utc[2] = day;
     //	      utc[3] = hour;
     //	      utc[4] = minute;
     //	      utc[5] = second;
     //	      utc[6] = msecond;
     //	   }
     //	   else
     //	      CHECK_INPUT(fail == 0 , "processObservation: Error in utc format. Give as mm/dd/yyyy:hh:mm:ss.ms "
     //			  << " or use utc=ignore" );
     //	}
     //     }
     else if( startswith("windowL=",token))
     {
        token += 8;
        winl = atof(token);
        winlset = true;
     }
     else if( startswith("windowR=",token))
     {
        token += 8;
        winr = atof(token);
        winrset = true;
     }
     else if( startswith("exclude=",token) )
     {
        token += 8;
	strncpy(exclstr,token,4);

	int c=0;
	while( c < 3 && exclstr[c] != '\0' )
	{
	   if( exclstr[c] == 'x' || exclstr[c] == 'e' )
	      usex=false;
	   if( exclstr[c] == 'y' || exclstr[c] == 'n' )
	      usey=false;
	   if( exclstr[c] == 'z' || exclstr[c] == 'u' )
	      usez=false;
	   c++;
	}
     }
     else if(startswith("filter=", token))
     {
        token += 7; // skip filter=
        if( strcmp(token,"0")==0 || strcmp(token,"no")==0 )
	   m_filter_observations = false;
     }
     else if( startswith("scalefactor=",token) )
     {
	token += 12;
	scalefactor = atof(token);
	scalefactor_set = true;
     }
     else
     {
        badOption("observation", token);
     }
     token = strtok(NULL, " \t");
  }  

  // Read from HDF5 file, and create time series data
#ifdef USE_HDF5
  bool is_obs = true;
  if( event_is_in_proc(event) )
  {
     event=global_to_local_event(event);
     readStationHDF5(this, inhdf5file, outhdf5file, writeEvery, downSample, mode, event, &a_GlobalTimeSeries, m_global_xmax, m_global_ymax, is_obs, winlset, winrset, winl, winr, usex, usey, usez, t0, scalefactor_set,  scalefactor);
  }
#else
  if (proc_zero())
    cout << "Using HDF5 station input but sw4 is not compiled with HDF5!"<< endl;
  return;
#endif

}

//-----------------------------------------------------------------------
void EW::processObservation( char* buffer, vector<vector<TimeSeries*> > & a_GlobalTimeSeries)
{
  double x=0.0, y=0.0, z=0.0;
  double lat = 0.0, lon = 0.0, depth = 0.0;
  float_sw4 t0 = 0;
  float_sw4 scalefactor=1;
  bool cartCoordSet = false, geoCoordSet = false;
  string fileName = "rec";
  string staName = "station";
  bool staNameGiven=false;
  
  int writeEvery = 0;

  bool dateSet = false;
  bool timeSet = false;
  bool topodepth = false;

  //  int utc[7];
  //  bool utcset = false;

  string date = "";
  string time = "";
  string sacfile1, sacfile2, sacfile3;
  string hdf5file = "";

  bool usgsformat = 1, sacformat=0, hdf5format = 0;
  TimeSeries::receiverMode mode=TimeSeries::Displacement;
  float_sw4 winl, winr;
  bool winlset=false, winrset=false;
  char exclstr[4]={'\0','\0','\0','\0'};
  bool usex=true, usey=true, usez=true;
  bool usgsfileset=false, sf1set=false, sf2set=false, sf3set=false;
  bool scalefactor_set=false;
  int event=0;

  char* token = strtok(buffer, " \t");
  m_filter_observations = true;

  CHECK_INPUT(strcmp("observation", token) == 0, "ERROR: not an observation line...: " << token);
  token = strtok(NULL, " \t");

  string err = "OBSERVATION Error: ";

  while (token != NULL)
  {
     if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
        break;
     if (startswith("x=", token))
     {
        CHECK_INPUT(!geoCoordSet,
                err << "observation command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 2; // skip x=
        cartCoordSet = true;
        x = atof(token);
        CHECK_INPUT(x >= 0.0,
		    "observation command: x must be greater than or equal to 0, not " << x);
        CHECK_INPUT(x <= m_global_xmax,
		    "observation command: x must be less than or equal to xmax, not " << x);
     }
     else if (startswith("y=", token))
     {
        CHECK_INPUT(!geoCoordSet,
                err << "observation command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 2; // skip y=
        cartCoordSet = true;
        y = atof(token);
        CHECK_INPUT(y >= 0.0,
                "observation command: y must be greater than or equal to 0, not " << y);
        CHECK_INPUT(y <= m_global_ymax,
		    "observation command: y must be less than or equal to ymax, not " << y);
     }
     else if (startswith("lat=", token))
     {
        CHECK_INPUT(!cartCoordSet,
                err << "observation command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 4; // skip lat=
        lat = atof(token);
        CHECK_INPUT(lat >= -90.0,
                "observation command: lat must be greater than or equal to -90 degrees, not " 
                << lat);
        CHECK_INPUT(lat <= 90.0,
                "observation command: lat must be less than or equal to 90 degrees, not "
                << lat);
        geoCoordSet = true;
     }
     else if (startswith("lon=", token))
     {
        CHECK_INPUT(!cartCoordSet,
                err << "observation command: Cannot set both a geographical (lat, lon) and a cartesian (x,y) coordinate");
        token += 4; // skip lon=
        lon = atof(token);
        CHECK_INPUT(lon >= -180.0,
                "observation command: lon must be greater or equal to -180 degrees, not " 
                << lon);
        CHECK_INPUT(lon <= 180.0,
                "observation command: lon must be less than or equal to 180 degrees, not "
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
		   "observation command: z must be less than or equal to zmax, not " << z);
     }
     else if (startswith("depth=", token))
     {
        token += 6; // skip depth=
       z = depth = atof(token);
       topodepth = true;
       CHECK_INPUT(depth >= 0.0,
	       err << "observation command: depth must be greater than or equal to zero");
       CHECK_INPUT(depth <= m_global_zmax,
		   "observation command: depth must be less than or equal to zmax, not " << depth);
// by depth we here mean depth below topography
     }
     else if(startswith("file=", token))
     {
        token += 5; // skip file=
        fileName = token;
        usgsfileset = true;
     }
     else if(startswith("event=",token))
     {
	token += 6;
	//	event = atoi(token);
	//	CHECK_INPUT( 0 <= event && event < m_nevent, err << "event no. "<< event << " out of range" );
	// Ignore if no events given
	if( m_nevents_specified > 0 )
	{
	   map<string,int>::iterator it = m_event_names.find(token);
           //	   CHECK_INPUT( it != m_event_names.end(), 
           //		     err << "event with name "<< token << " not found" );
           //	   event = it->second;
             if( it != m_event_names.end() )
                event = it->second;
             else if( proc_zero() )
                std::cout << "Observation warning: event with name " << token << " not found" << std::endl;
	}
     }
     else if (startswith("sta=", token))
     {
        token += strlen("sta=");
        staName = token;
	staNameGiven=true;
     }
// (small) shifts of the observation in time can be used to compensate for incorrect velocites
// in the material model
     else if(startswith("shift=", token))
     {
        token += 6; // skip shift=
        t0 = atof(token);
     }
     //     else if( startswith("utc=",token) )
     //     {
     //	token += 4;
     //        if( strcmp("ignore",token)==0 || strcmp("off",token)==0 )
     //	   ignore_utc = true;
     //	else
     //	{
     //	   int year,month,day,hour,minute,second,msecond, fail;
     //	   // Format: 01/04/2012:17:34:45.2343  (Month/Day/Year:Hour:Min:Sec.fraction)
     //	   parsedate( token, year, month, day, hour, minute, second, msecond, fail );
     //	   if( fail == 0 )
     //	   {
     //              utcset = true;
     //	      utc[0] = year;
     //	      utc[1] = month;
     //	      utc[2] = day;
     //	      utc[3] = hour;
     //	      utc[4] = minute;
     //	      utc[5] = second;
     //	      utc[6] = msecond;
     //	   }
     //	   else
     //	      CHECK_INPUT(fail == 0 , "processObservation: Error in utc format. Give as mm/dd/yyyy:hh:mm:ss.ms "
     //			  << " or use utc=ignore" );
     //	}
     //     }
     else if( startswith("windowL=",token))
     {
        token += 8;
        winl = atof(token);
        winlset = true;
     }
     else if( startswith("windowR=",token))
     {
        token += 8;
        winr = atof(token);
        winrset = true;
     }
     else if( startswith("exclude=",token) )
     {
        token += 8;
	strncpy(exclstr,token,4);

	int c=0;
	while( c < 3 && exclstr[c] != '\0' )
	{
	   if( exclstr[c] == 'x' || exclstr[c] == 'e' )
	      usex=false;
	   if( exclstr[c] == 'y' || exclstr[c] == 'n' )
	      usey=false;
	   if( exclstr[c] == 'z' || exclstr[c] == 'u' )
	      usez=false;
	   c++;
	}
     }
     else if(startswith("filter=", token))
     {
        token += 7; // skip filter=
        if( strcmp(token,"0")==0 || strcmp(token,"no")==0 )
	   m_filter_observations = false;
     }
     else if( startswith("sacfile1=",token) )
     {
        token += 9;
        sacfile1 += token;
	sf1set = true;
     }
     else if( startswith("sacfile2=",token) )
     {
        token += 9;
        sacfile2 += token;
	sf2set = true;
     }
     else if( startswith("sacfile3=",token) )
     {
        token += 9;
        sacfile3 += token;
	sf3set = true;
     }
     else if( startswith("scalefactor=",token) )
     {
	token += 12;
	scalefactor = atof(token);
	scalefactor_set = true;
     }
     else
     {
        badOption("observation", token);
     }
     token = strtok(NULL, " \t");
  }  

  // Make sure either one usgsfile or three sac files are input.
  if( event_is_in_proc(event) )
  {
     int eglobal = event;
     event = global_to_local_event(event);
     if( usgsfileset )
     {
        CHECK_INPUT( !sf1set && !sf2set && !sf3set, "processObservation, Error: can not give both usgs file and sacfiles" );
     }
     else 
     {
        CHECK_INPUT( sf1set && sf2set && sf3set, "processObservation, Error: must give at least three sac files" );
// Find a name for the SAC station
        int l = sacfile1.length();
        if( sacfile1.substr(l-4,4) == ".sac" )
           fileName = sacfile1.substr(0,l-4);
        else
           fileName = sacfile1;

// Read sac header to figure out the position
// Use only one of the files, more thorough checking later, in TimeSeries.readSACfile.
        float_sw4 latlon[2];
        if( m_myRank == 0 )
        {
           string fname= mObsPath[eglobal];
           fname += sacfile1;
           FILE* fd=fopen(fname.c_str(),"r");
           CHECK_INPUT( fd != NULL, "processObservation: ERROR: sac file " << sacfile1 << " could not be opened" );
           float float70[70];
           size_t nr = fread(float70, sizeof(float), 70, fd );
           CHECK_INPUT( nr == 70, "processObservation: ERROR, could not read float part of header of " << sacfile1 );
           latlon[0] = float70[31];
           latlon[1] = float70[32];
           CHECK_INPUT( latlon[0] != -12345 && latlon[1] != -12345, 
                       "processObservation: ERROR, sac file does not contain station coordinates " << sacfile1);
           fclose(fd);
        }
        MPI_Bcast( latlon, 2, MPI_DOUBLE, 0, m_1d_communicator );
        if( geoCoordSet && ( fabs(lat-latlon[0])<1e-10 && fabs(lon-latlon[1])<1e-10 ))
        {
           if( m_myRank == 0 )
              cout << "processObservation: WARNING station (lat,lon) on sac file do not match input (lat,lon)" << endl;
        }
        if( !cartCoordSet && !geoCoordSet )
        {
           geoCoordSet = true;
           lat = latlon[0];
           lon = latlon[1];
        }
     }
     if (geoCoordSet)
     {
        computeCartesianCoord(x, y, lon, lat);
// check if (x,y) is within the computational domain
     }

     if (!staNameGiven)
        staName = fileName;

     bool inCurvilinear=false;
//
// AP: This test is incorrect because we don't know the elevation of the observation
//
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
           stringstream observationerr;
  
           observationerr << endl 
                          << "***************************************************" << endl
                          << " WARNING:  OBSERVATION positioned outside grid!" << endl;
           observationerr << " No OBSERVATION file will be generated for file = " << fileName << endl;
           if (geoCoordSet)
           {
              observationerr << " @ lon=" << lon << " lat=" << lat << " depth=" << depth << endl << endl;
           }
           else
           {
              observationerr << " @ x=" << x << " y=" << y << " z=" << z << endl << endl;
           }
      
           observationerr << "***************************************************" << endl;
           cerr << observationerr.str();
           cerr.flush();
        }
     }
     else
     {
        TimeSeries *ts_ptr = new TimeSeries(this, fileName, staName, mode, sacformat, usgsformat, hdf5format, hdf5file, x, y, depth, 
					topodepth, writeEvery, 1, true, event );
    // Read in file. 
    // ignore_utc=true, ignores UTC read from file, instead uses the default utc = simulation utc as reference.
    //        This is useful for synthetic data.

        if( usgsfileset )
           ts_ptr->readFile( this, false );
        else
           ts_ptr->readSACfiles( this, sacfile1.c_str(), sacfile2.c_str(), sacfile3.c_str(), false );

// Set reference UTC to simulation UTC, for easier plotting.
        ts_ptr->set_utc_to_simulation_utc();

// Set window, in simulation time
        if( winlset || winrset )
        {
           if( winlset && !winrset )
              winr = 1e38;
           if( !winlset && winrset )
              winl = -1;
           ts_ptr->set_window( winl, winr );
        }

// Exclude some components
        if( !usex || !usey || !usez )
           ts_ptr->exclude_component( usex, usey, usez );

// Add extra shift from command line, use with care.
        if( t0 != 0 )
           ts_ptr->add_shift( t0 );

        // Set scale factor if given
        if( scalefactor_set )
           ts_ptr->set_scalefactor( scalefactor );

// include the observation in the global list
        a_GlobalTimeSeries[event].push_back(ts_ptr);
     }
  }
}

//-----------------------------------------------------------------------
void EW::processScaleFactors( char* buffer )
{
  char* token = strtok(buffer, " \t");

  CHECK_INPUT(strcmp("scalefactors", token) == 0, "ERROR: not a scalefactors line...: " << token);
  token = strtok(NULL, " \t");
  float_sw4 x0=1, y0=1, z0=1, mxx=1, mxy=1, mxz=1, myy=1, myz=1, mzz=1, t0=1, freq=1; 

  string err = "SCALEFACTORS Error: ";

  while (token != NULL)
  {
     if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
        break;
     if (startswith("x0=", token))
     {
        token += 3;
        x0 = atof(token);
        CHECK_INPUT(x0 > 0.0,
		    "scalefactors command: x0 must be greater than 0, not " << x0);
     }
     else if( startswith("y0=",token))
     {
        token += 3;
        y0 = atof(token);
        CHECK_INPUT(y0 > 0.0,
		    "scalefactors command: y0 must be greater than 0, not " << y0);
     }
     else if( startswith("z0=",token))
     {
        token += 3;
        z0 = atof(token);
        CHECK_INPUT(z0 > 0.0,
		    "scalefactors command: y0 must be greater than 0, not " << z0);
     }
     else if( startswith("Mxx=",token))
     {
        token += 4;
        mxx = atof(token);
        CHECK_INPUT(mxx > 0.0,
		    "scalefactors command: Mxx must be greater than 0, not " << mxx);
     }
     else if( startswith("Mxy=",token))
     {
        token += 4;
        mxy = atof(token);
        CHECK_INPUT(mxy > 0.0,
		    "scalefactors command: Mxy must be greater than 0, not " << mxy);
     }
     else if( startswith("Mxz=",token))
     {
        token += 4;
        mxz = atof(token);
        CHECK_INPUT(mxz > 0.0,
		    "scalefactors command: Mxz must be greater than 0, not " << mxz);
     }
     else if( startswith("Myy=",token))
     {
        token += 4;
        myy = atof(token);
        CHECK_INPUT(myy > 0.0,
		    "scalefactors command: Myy must be greater than 0, not " << myy);
     }
     else if( startswith("Myz=",token))
     {
        token += 4;
        myz = atof(token);
        CHECK_INPUT(myz > 0.0,
		    "scalefactors command: Myz must be greater than 0, not " << myz);
     }
     else if( startswith("Mzz=",token))
     {
        token += 4;
        mzz = atof(token);
        CHECK_INPUT(mzz > 0.0,
		    "scalefactors command: Mzz must be greater than 0, not " << mzz);
     }
     else if( startswith("t0=",token))
     {
        token += 3;
        t0 = atof(token);
        CHECK_INPUT(t0 > 0.0,
		    "scalefactors command: t0  must be greater than 0, not " << t0 );
     }
     else if( startswith("freq=",token))
     {
        token += 5;
        freq = atof(token);
        CHECK_INPUT(freq > 0.0,
		    "scalefactors command: freq  must be greater than 0, not " << freq );
     }
     else
     {
        badOption("scalefactors", token);
     }
     token = strtok(NULL, " \t");
  }  
  m_scalefactors[0] = x0;
  m_scalefactors[1] = y0;
  m_scalefactors[2] = z0;
  m_scalefactors[3] = mxx;
  m_scalefactors[4] = mxy;
  m_scalefactors[5] = mxz;
  m_scalefactors[6] = myy;
  m_scalefactors[7] = myz;
  m_scalefactors[8] = mzz;
  m_scalefactors[9] = t0;
  m_scalefactors[10]= freq;
}

//-----------------------------------------------------------------------
void EW::processCG( char* buffer )
{
  char* token = strtok(buffer, " \t");

  CHECK_INPUT(strcmp("cg", token) == 0, "ERROR: not a cg line...: " << token);
  token = strtok(NULL, " \t");
  m_maxit = 11;
  m_maxrestart = 0;
  m_tolerance = 1e-6;
  m_iniguess_pos = false;
  m_iniguess_t0fr = false;
  m_iniguess_mom = false;
  m_compute_scalefactors=false;
  m_cgstepselection = 0;
  m_cgvarcase = 0;
  m_cgfletcherreeves = true;
  m_do_linesearch = true;
  m_opt_method = 1;
  m_lbfgs_m = 4;
  m_opt_testing = false;

  string err = "CG Error: ";

  while (token != NULL)
  {
     if (startswith("#", token) || startswith(" ", buffer))
        // Ignore commented lines and lines with just a space.
        break;
     //     if (startswith("maxit=", token))
     //     {
     //        token += 6;
     //        m_maxit = atoi(token);
     //        CHECK_INPUT(m_maxit >= 0,
     //		    "cg command: maxit must be greater than or equal to 0, not " << m_maxit );
     //     }
     //     else if( startswith("maxouterit=",token) )
     if( startswith("maxouterit=",token) )
     {
        token += 11;
        m_maxrestart = atoi(token);
        CHECK_INPUT(m_maxrestart >= 0,
		    "cg command: maxouterit must be greater than or equal to 0, not " << m_maxrestart );
     }
     else if( startswith("tolerance=",token) )
     {
        token += 10;
        m_tolerance = atof(token);
        CHECK_INPUT(m_tolerance >= 0,
		    "cg command: tolerance must be greater than or equal to 0, not " << m_tolerance );
     }
     else if( startswith("initialguess=",token) )
     {
        token += 13;
        if( strcmp(token,"useSource")==0 || strcmp(token,"usesource")==0 )
   	   m_iniguess_pos = m_iniguess_t0fr = m_iniguess_mom = false;
	else if( strcmp(token,"estimate") == 0 )
	   m_iniguess_pos = m_iniguess_t0fr = m_iniguess_mom = true;
        else if( strcmp(token,"estimatePos")==0 || strcmp(token,"estimatepos")==0 )
	{
           m_iniguess_pos  = true;
	   m_iniguess_t0fr = false;
	   m_iniguess_mom  = false;
	}
        else if( strcmp(token,"estimateT0Pos")==0 )
	{
           m_iniguess_pos  = true;
	   m_iniguess_t0fr = true;
	   m_iniguess_mom  = false;
	}
        else if( strcmp(token,"estimateM")==0 )
	{
           m_iniguess_pos  = false;
	   m_iniguess_t0fr = false;
	   m_iniguess_mom  = true;
	}
        else
	   CHECK_INPUT( false,
		     "cg command: initialguess value " << token << " not understood");
     }
     else if( startswith("estimateshifts=",token) )
     {
        token += 15;
	m_iniguess_shifts = strcmp("yes",token)==0||strcmp("true",token)==0||strcmp("1",token)==0;
     }
     else if( startswith("write_initial_ts=",token) )
     {
        token += 17;
        int val=atoi(token);
	if( val == 1 )
	   m_output_initial_seismograms = true;
	else if( val == 0 )
	   m_output_initial_seismograms = false;
        CHECK_INPUT( (val == 1) || (val== 0) ,
		    "cg command: write_initial_ts must be equal to 0 or 1, not " << val );

     }
     else if( startswith("scalefactors=",token) )
     {
        token += 13;
        if( strcmp(token,"useinput")==0 )
	   m_compute_scalefactors = false;
	else if( strcmp(token,"estimate") == 0 )
	   m_compute_scalefactors = true;
        else
	   CHECK_INPUT( false,
		     "cg command: scalefactors value " << token << " not understood");
     }
     else if( startswith("linesearch=",token) )
     {
        token += 11;
        if( strcmp(token,"on")==0 )
	   m_do_linesearch = true;
	else if( strcmp(token,"off") == 0 )
	   m_do_linesearch = false;
        else
	   CHECK_INPUT( false,
		     "cg command: linesearch value " << token << " not understood");
     }
     else if( startswith("steptype=",token) )
     {
        token += 9;
        if( strcmp(token,"misfit")==0 )
	   m_cgstepselection = 0;
	else if( strcmp(token,"hessian") == 0 )
	   m_cgstepselection = 1;
        else
	   CHECK_INPUT( false,
		     "cg command: steptype value " << token << " not understood");
     }
     else if( startswith("optmethod=",token) )
     {
        token += 10;
        if( strcmp(token,"fletcher-reeves")==0 )
	{
           m_opt_method = 1;
	   m_cgfletcherreeves = true;
	}
	else if( strcmp(token,"polak-ribiere") == 0 )
	{
           m_opt_method = 1;
	   m_cgfletcherreeves = false;
	}
        else if( strcmp(token,"l-BFGS") == 0 )
	   m_opt_method = 2;
        else if( strcmp(token,"BFGS") == 0 )
	   m_opt_method = 3;
        else if( strcmp(token,"steepest-descent") == 0 )
	   m_opt_method = 4;
        else
	   CHECK_INPUT( false,
		     "cg command: optmethod value " << token << " not understood");
     }
     else if( startswith("lbfgsvectors=",token) )
     {
        token += 13;
	m_lbfgs_m = atoi(token);
	CHECK_INPUT( m_lbfgs_m > 0,
		     "Number of l-BFGS vectors must be positive. Input value = " << m_lbfgs_m );
     }
     else if( startswith("solvefor=" , token) )
     {
	token += 9;
        if( strcmp(token,"posMt0freq")==0 )
	   m_cgvarcase = 0;
	else if( strcmp(token,"posMt0") == 0 )
	   m_cgvarcase = 1;
	else if( strcmp(token,"posM") == 0 )
	   m_cgvarcase = 2;
	else if( strcmp(token,"posMobs") == 0 )
	   m_cgvarcase = 3;
        else
	   CHECK_INPUT( false,
		     "cg command: solvefor value " << token << " not understood");
     }
     else if( startswith("opttest=",token) )
     {
        token += 8;
	m_opt_testing =  strcmp(token,"yes")== 0 || strcmp(token,"1")==0;
     }
     else
     {
        badOption("cg", token);
     }
     token = strtok(NULL, " \t");
  }  
}

//-----------------------------------------------------------------------
void EW::processMaterialPfile(char* buffer)
{
  string name = "pfile";

  // Used for pfiles
  string filename = "NONE";
  string directory = "NONE";
  float_sw4 a_ppm=0.,vpmin_ppm=0.,vsmin_ppm=0,rhomin_ppm=0.;
  string cflatten = "NONE";
  bool flatten = false;
  bool coords_geographic = true;
  int nstenc = 5;

  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("pfile", token) == 0,
	      "ERROR: material data can only be set by an pfile line, not: " << token);

  string err = token;
  err += " Error: ";
  token = strtok(NULL, " \t");


  while (token != NULL)
    {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	break;
      //      else if (startswith("a=", token))
      //      {
      //         token += 2; // skip a=
      //         a_ppm = atof(token);
      //      }
      else if( startswith("smoothingsize=",token) )
      {
	token += 14;
	nstenc = atoi(token);
	VERIFY2( nstenc >= 1 ,
		 "processMaterialPfile Error: nstenc is " << nstenc << "but should be >= 1\n" );
      }
      else if (startswith("vpmin=", token))
      {
	token += 6; // skip vpmin=
	vpmin_ppm = atof(token);
      }
      else if (startswith("vsmin=", token))
      {
	token += 6; // skip vsmin=
	vsmin_ppm = atof(token);
      }
      else if (startswith("rhomin=", token))
      {
	token += 7; // skip rhomin=
	rhomin_ppm = atof(token);
      }
      else if (startswith("flatten=", token))
      {
	token += 8; // skip flatten=
	cflatten = token;
	VERIFY2( (int)cflatten.find('T')>=0 || (int)cflatten.find('t')>=0 ||
		 (int)cflatten.find('F')>=0 || (int)cflatten.find('f')>=0,
		 "processMaterialPfile Error: value of flatten unclear\n" );
	if ((int)cflatten.find('T')>=0||(int)cflatten.find('t')>=0)
	  flatten=true;
	else if ((int)cflatten.find('F')>=0||(int)cflatten.find('f')>=0)
	  flatten=false;
	else
	  flatten=false;
	 
      }
      else if (startswith("filename=", token))
      {
	token += 9; // skip filename=
	filename = token;
      }
      else if (startswith("directory=", token))
      {
	token += 10; // skip directory=
	directory = token;
      }
      else if (startswith("style=", token))
      {
	token += 6; // skip style=
	if( strcmp(token,"geographic") == 0 || strcmp(token,"Geographic")==0 )
	  coords_geographic = true;
	else if( strcmp(token,"cartesian") == 0 || strcmp(token,"Cartesian")==0 )
	  coords_geographic = false;
	else
	  CHECK_INPUT( false, "processMaterialPfile Error: style= " << token << " not recognized\n" );
      }
      else
      {
	cout << token << " is not a pfile option " << endl;
      }
      token = strtok(NULL, " \t");
    }
  // End parsing...

  //----------------------------------------------------------------
  // Check parameters
  //----------------------------------------------------------------
  if (strcmp(directory.c_str(),"NONE")==0)
  {
      directory = string("./");
  }

  if (m_myRank == 0)
  {
     cout << "*** Reading data from Pfile " << filename << " in directory " << directory << endl;
  }

  MaterialPfile* pf = new MaterialPfile( this,
					filename, directory,nstenc,vpmin_ppm,vsmin_ppm,rhomin_ppm,flatten,
					coords_geographic );

  add_mtrl_block( pf  );
     
}

//-----------------------------------------------------------------------
void EW::processMaterialEtree(char* buffer)
{
#ifdef ENABLE_ETREE
  string name = "Efile";

  // Used only for efiles
  string model = "SF";
  string etreefile = "NONE";
  string xetreefile = "NONE";
  string etreelogfile = "NONE";
  const char* queryStyle = "MAXRES";

  double eFileResolution = -1;
  //  if (!mCheckMode)
  //    resolution = mGridSize[0];

  //  float vpmin = 0.0;
  //  bool vpMinSet = false;
  //  float vsmin = 0.0;
  //  bool vsMinSet = false;
  string accessmode = "parallel";
  cencalvm::storage::Geometry* efileGeom = 0;

  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("efile", token) == 0,
           "ERROR: material data can only be set by an efile line, not: " << token);
  string materialCommandType = token;

  string err = token;
  err += " Error: ";

  token = strtok(NULL, " \t");

//  InitialConditionBlock* b = new InitialConditionBlock;

  while (token != NULL)
    {
      // while there are tokens in the string still
       if (startswith("#", token) || startswith(" ", buffer))
          // Ignore commented lines and lines with just a space.
          break;
//       else if (startswith("model=", token))
//       {
//          token += 6; // skip model=
//          model = token;
//       }
      else if (startswith("etree=", token))
      {
         token += 6; // skip etree=
         etreefile = token;
      }
      else if (startswith("xetree=", token))
      {
         token += 7; // skip xetree=
         xetreefile = token;
      }
      else if (startswith("logfile=", token))
      {
         token += 8; // skip logfile=
         etreelogfile = token;
      }
      else if (startswith("query=", token))
      {
         token += strlen("query=");
         queryStyle = token;
         CHECK_INPUT(strcmp(queryStyle, "FIXEDRES") == 0 || 
                 strcmp(queryStyle, "MAXRES") == 0,
                 err << "query can only be set to FIXEDRES or MAXRES, not: " << queryStyle);
      }
       //      else if (startswith("vsmin=", token))
       //      {
       //         token += strlen("vsmin=");
       //         vsmin = atof(token);
       //         vsMinSet = true;
       //      }
       //      else if (startswith("vpmin=", token))
       //      {
       //         token += strlen("vpmin=");
       //         vpmin = atof(token);
       //vpMinSet = true;
       //      }
     else if( startswith("resolution=", token ) )
     {
       token += 11;
       eFileResolution = atof(token);
       CHECK_INPUT(eFileResolution>0.,"Resolution must be positive, not " << eFileResolution);
     }
      else if (startswith("access=", token))
      {
         token += strlen("access=");
         CHECK_INPUT(strcmp(token, "parallel") == 0 ||
                 strcmp(token, "serial") == 0,
                 err << "access attribute can only be set to serial, or parallel, not: " << token);
         accessmode = token;
      }
      else
      {
         badOption(materialCommandType, token);
      }
      token = strtok(NULL, " \t");
    }
  // End parsing...
  
  CHECK_INPUT(model != "",
          err << "the model attribute must be set");


  //----------------------------------------------------------------
  // Check parameters
  //----------------------------------------------------------------
  if (m_myRank == 0)
  {
     cout << "\t Efile " << model << " Model specified..." << endl;
     cout << "\t Looking for etree model: " << etreefile << endl;
     if (xetreefile != "NONE")
     {
        cout << "\t Looking for xetree model: " << xetreefile << endl;
     }
  }
  //  if (mCheckMode) 
  //  {
  //     if (access(etreefile.c_str(), R_OK) != 0)
  //        cout << "***ERROR: No read permission on etree file: " << etreefile << endl;
  //     if (xetreefile != "NONE")
  //     {
  //        if (access(xetreefile.c_str(), R_OK) != 0)
  //           cout << "***ERROR: No read permission on xefile: " << xetreefile << endl;
  //     }
  //     return;
  //  }
  
  // checkmode does the copy from above, if needed, nothing else
  //  if (mCheckMode) return;

  name = model + " " + materialCommandType;

  std::vector<int> temp3D;
  std::vector<int> temp2D;

  if (eFileResolution < 0.)
     eFileResolution = mGridSize[mNumberOfGrids - 1];

  EtreeFile* ef = new EtreeFile( this,
                                accessmode, etreefile, xetreefile, model,
                                etreelogfile, queryStyle, eFileResolution);
  
  double latse, lonse, latsw, lonsw, latne, lonne, latnw, lonnw;
  ef->getcorners( latse, lonse, latsw, lonsw, latne, lonne, latnw, lonnw );
  //  ef->getGeoBox()->getBounds(NW, NE, SW, SE);
  double nwxyz[3];
  double nexyz[3];
  double swxyz[3];
  double sexyz[3];
  //  computeCartesianCoord(nwxyz[0],nwxyz[1],nwxyz[2],NW);
  //  computeCartesianCoord(nexyz[0],nexyz[1],nexyz[2],NE);
  //  computeCartesianCoord(swxyz[0],swxyz[1],swxyz[2],SW);
  //  computeCartesianCoord(sexyz[0],sexyz[1],sexyz[2],SE);

  computeCartesianCoord(sexyz[0],sexyz[1],lonse,latse);
  computeCartesianCoord(swxyz[0],swxyz[1],lonsw,latsw);
  computeCartesianCoord(nexyz[0],nexyz[1],lonne,latne);
  computeCartesianCoord(nwxyz[0],nwxyz[1],lonnw,latnw);

  // constrain on x domain
  double etreeXmin = min(min(min(nexyz[0], nwxyz[0]), sexyz[0]), swxyz[0]);
  double etreeXmax = max(max(max(nexyz[0], nwxyz[0]), sexyz[0]), swxyz[0]);
  double etreeYmin = min(min(min(nexyz[1], nwxyz[1]), sexyz[1]), swxyz[1]);
  double etreeYmax = max(max(max(nexyz[1], nwxyz[1]), sexyz[1]), swxyz[1]);
  //  double etreeZmin = min(min(min(nexyz[2], nwxyz[2]), sexyz[2]), swxyz[2]);
  //  double etreeZmax = max(max(max(nexyz[2], nwxyz[2]), sexyz[2]), swxyz[2]);
  
  if( getVerbosity() >=2 && m_myRank == 0 )
    printf("Horizontal extent of etree %s xmin=%e, xmax=%e, ymin=%e, ymax=%e\n", name.c_str(), etreeXmin, etreeXmax, etreeYmin, etreeYmax);
  
//   // needed to write topographic image slices
  mEtreeFile = ef;

  add_mtrl_block( ef  );  
     
#else
  CHECK_INPUT(0, "Error: Etree support not compiled into SW4 (use -DENABLE_ETREE)");
#endif
}

void EW::processMaterialIfile( char* buffer )
{
  bool x1set=false, x2set=false, y1set=false, y2set=false, 
    z1set=false, z2set=false;

  float_sw4 x1=0.0, x2=0.0, y1=0.0, y2=0.0, z1=0.0, z2=0.0;
  int i1=-1, i2=-1, j1=-1, j2=-1, k1=-1, k2=-1;

  string name = "Ifile";

  char* token = strtok(buffer, " \t");
  CHECK_INPUT(strcmp("ifile", token) == 0,
	      "ERROR: material ifile can be set by a ifile line, not: " << token);

  string err = token, filename="NONE";
  bool CartesianFormat=false;
  
  err += " Error: ";

  token = strtok(NULL, " \t");

  float_sw4 vp=-1, vs=-1, rho=-1, ps=-1, materialID=-1, freq=1;
  float_sw4 vpgrad=0, vsgrad=0, rhograd=0;
  float_sw4 vp2=0, vs2=0, rho2=0;
  
  bool gotFileName=false;
  
  while (token != NULL)
  {
    // while there are tokens in the string still
    if (startswith("#", token) || startswith(" ", buffer))
      // Ignore commented lines and lines with just a space.
      break;
//                  1234567890
    if (startswith("filename=", token) )
    {
      token += 9; // skip filename=
      filename = token;
      gotFileName = true;
    }
// Cartesian or (lon,lat) format?
//                       1234567890
    else if (startswith("input=", token))
    {
      token += 6; // skip input=
      if (strcmp("cartesian", token) == 0)
      {
	CartesianFormat=true;
      }
// change option to geographic, but keeping grid for backwards compatibility
      else if (strcmp("geographic", token) == 0 || strcmp("grid", token) == 0)
      {
	CartesianFormat=false;
      }
      else
      {
	badOption("ifile> input", token);
      }

    }
    else
    {
      badOption("ifile", token);
    }
    token = strtok(NULL, " \t");
  }
  // End parsing...
  
  CHECK_INPUT(gotFileName, "ERROR: no filename specified in ifile command. ");

  if(mVerbose >=2 &&  m_myRank == 0 )
  {
    cout << "**** Ifile parameters: *****" << endl;
    cout << "filename=" << filename << endl;
    cout << "CartesianFormat=" << CartesianFormat << endl;
  }

// add this material specificaiton to the global model
  MaterialIfile* bl = new MaterialIfile( this, filename, CartesianFormat );
  add_mtrl_block( bl );
}

//-----------------------------------------------------------------------
void EW::processMaterialVimaterial(char* buffer)
{
   string name = "vimtrl";
   string path=".";
   string rho, mu, lambda, qp, qs;
   bool cpset=false, csset=false, rhoset=false, muset=false, lambdaset=false;
   bool rhomula;
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("vimaterial", token) == 0,
	       "ERROR: not a vimaterial line: " << token);
   string err = token;
   err += " Error: ";
   token = strtok(NULL, " \t");
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      //      else if (startswith("a=", token))
      //      {
      //         token += 2; // skip a=
      //         a_ppm = atof(token);
      //      }
      else if( startswith("path=",token) )
      {
	 token += 5;
	 path = token;
      }
      else if( startswith("rho=",token) )
      {
	 token += 4;
	 rho = token;
	 rhoset = true;
      }
      else if( startswith("mu=",token) )
      {
	 token += 3;
	 mu = token;
	 muset = true;
      }
      else if( startswith("lambda=",token) )
      {
	 token += 7;
	 lambda = token;
	 lambdaset = true;
      }
      else if( startswith("vs=",token) )
      {
	 token += 3;
	 mu = token;
	 csset = true;
      }
      else if( startswith("vp=",token) )
      {
	 token += 3;
	 lambda = token;
	 cpset = true;
      }
      else if( startswith("qs=",token) )
      {
	 token += 3;
	 qs = token;
      }
      else if( startswith("qp=",token) )
      {
	 token += 3;
	 qp = token;
      }
      else
      {
	 badOption("vimaterial", token);
      }
      token = strtok(NULL, " \t");
   }
   if( rhoset && muset && lambdaset )
      rhomula = true;
   else if( rhoset && csset && cpset )
      rhomula = false;
   else
   {
      CHECK_INPUT( 0, "Error parsing vimaterial, must set (rho,mu,lambda) or (rho,vs,vp) ");
   }
   MaterialData *mdata = new MaterialVolimagefile( this, rhomula, path, rho, mu, lambda, qp, qs );
   add_mtrl_block(mdata);
}

//-----------------------------------------------------------------------
void EW::processMaterialRfile(char* buffer)
{
   string name = "rfile";

  // Used for pfiles
   string filename = "NONE";
   string directory = "NONE";
   float_sw4 a_ppm=0.,vpmin_ppm=0.,vsmin_ppm=0,rhomin_ppm=0.;
   string cflatten = "NONE";
   bool flatten = false;
   bool coords_geographic = true;
   int nstenc = 5;
   int bufsize = 200000;  // Parallel IO buffer, in number of grid points.

   char* token = strtok(buffer, " \t");
  //  CHECK_INPUT(strcmp("rfile", token) == 0,
  //	      "ERROR: material data can only be set by an rfile line, not: " << token);

   string err = token;
   err += " Error: ";
   token = strtok(NULL, " \t");

   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      //      else if (startswith("a=", token))
      //      {
      //         token += 2; // skip a=
      //         a_ppm = atof(token);
      //      }
      else if (startswith("filename=", token))
      {
	 token += 9; // skip filename=
	 filename = token;
      }
      else if (startswith("directory=", token))
      {
	 token += 10; // skip directory=
	 directory = token;
      }
      else if (startswith("bufsize=", token))
      {
	 token += 8;
	 bufsize = atoi(token);
	 CHECK_INPUT( bufsize > 10 && bufsize < 1e9, "ParseInputfile: rfile bufsize = " <<
		      bufsize << " out of allowed range" );
      }
      else
      {
	 cout << token << " is not a rfile option " << endl;
      }
      token = strtok(NULL, " \t");
   }
  // End parsing...

  //----------------------------------------------------------------
  // Check parameters
  //----------------------------------------------------------------
   if (strcmp(directory.c_str(),"NONE")==0)
      directory = string("./");

   if (m_myRank == 0)
      cout << "*** Reading data from Rfile " << filename << " in directory " << directory << endl;

   MaterialRfile* rf = new MaterialRfile( this, filename, directory, bufsize );
   add_mtrl_block( rf  );
}

//-----------------------------------------------------------------------
void EW::processMaterialSfile(char* buffer)
{
   string name = "sfile";
   string filename = "NONE";
   string directory = "NONE";
   float_sw4 a_ppm=0.,vpmin_ppm=0.,vsmin_ppm=0,rhomin_ppm=0.;
   string cflatten = "NONE";
   bool flatten = false;
   bool coords_geographic = true;
   int nstenc = 5;
   int bufsize = 200000;  // Parallel IO buffer, in number of grid points.

   char* token = strtok(buffer, " \t");
  //  CHECK_INPUT(strcmp("rfile", token) == 0,
  //	      "ERROR: material data can only be set by an rfile line, not: " << token);

   string err = token;
   err += " Error: ";
   token = strtok(NULL, " \t");

   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      //      else if (startswith("a=", token))
      //      {
      //         token += 2; // skip a=
      //         a_ppm = atof(token);
      //      }
      else if (startswith("filename=", token))
      {
	 token += 9; // skip filename=
	 filename = token;
      }
      else if (startswith("directory=", token))
      {
	 token += 10; // skip directory=
	 directory = token;
      }
      else
      {
	 cout << token << " is not a sfile option " << endl;
      }
      token = strtok(NULL, " \t");
   }
  // End parsing...

  //----------------------------------------------------------------
  // Check parameters
  //----------------------------------------------------------------
  if (strcmp(directory.c_str(),"NONE")==0)
     directory = string("./");

  if (m_myRank == 0)
     cout << "*** Using Sfile " << filename << " in directory " << directory << endl;

  MaterialSfile* sf = new MaterialSfile(this, filename, directory);
  add_mtrl_block( sf  );
}


//-----------------------------------------------------------------------
void EW::processMaterialInvtest(char* buffer)
{
   int nr=1;
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("invtestmaterial", token) == 0,
	       "ERROR: not an invtestmaterial line: " << token);
   token = strtok(NULL, " \t");
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("type=",token) )
      {
	 token += 5;
         if( strcmp(token,"sineperturbation")==0 )
	    nr = 1;
	 else if( strcmp(token,"box")==0 )
	    nr = 2;
	 else if( strcmp(token,"lohsine")==0 )
	    nr = 3;
	 else if( strcmp(token,"smoothlayer")==0 )
	    nr = 4;
         else
	    CHECK_INPUT( 0, "Error invtestmaterial, type = " << token << " not recognized" );
      }
      else
      {
	 badOption("invtestmaterial", token);
      }
      token = strtok(NULL, " \t");
   }
   MaterialData *mdata = new MaterialInvtest( this, nr );
   add_mtrl_block(mdata);
}

//-----------------------------------------------------------------------
//void EW::processRandomize(char* buffer)
//{
//    char* token = strtok(buffer, " \t");
//    CHECK_INPUT(strcmp("randomize", token) == 0,
// 	       "ERROR: not a randomize line: " << token);
//    token = strtok(NULL, " \t");
//    bool lengthscaleset=false, lengthscalezset=false;
//    m_random_dist = m_random_distz = 100;
//    m_random_sdlimit  = 3;
//    m_random_amp      = 0.1;
//    m_random_amp_grad = 0;
//    m_random_sdlimit  = 3;
//    m_random_seed[0]  = 1234;
//    m_random_seed[1]  = 5678;
//    m_random_seed[2]  = 9876;

//    m_randomize = true;
//    while (token != NULL)
//    {
//       // while there are tokens in the string still
//       if (startswith("#", token) || startswith(" ", buffer))
// 	// Ignore commented lines and lines with just a space.
// 	 break;
//       else if( startswith("amplitude=",token) )
//       {
// 	 token += 10;
// 	 m_random_amp = atof(token);
// 	 CHECK_INPUT( m_random_amp>0, "Error randomize, amplitude must be > 0, not " << token);
//       }
//       else if( startswith("gradient=",token) )
//       {
// 	 token += 9;
// 	 m_random_amp_grad = atof(token);
//       }
//       else if( startswith("lengthscale=",token) )
//       {
// 	 token += 12;
// 	 m_random_dist = atof(token);
// 	 lengthscaleset = true;
// 	 CHECK_INPUT( m_random_dist>0, "Error randomize, dist must be > 0, not " << token);
//       }
//       else if( startswith("lengthscalez=",token) )
//       {
// 	 token += 13;
// 	 m_random_distz = atof(token);
// 	 lengthscalezset = true;
// 	 CHECK_INPUT( m_random_distz>0, "Error randomize, distz must be > 0, not " << token);
//       }
//       else if( startswith("sdthreshold=",token) )
//       {
// 	 token += 12;
// 	 m_random_sdlimit = atof(token);
// 	 CHECK_INPUT( m_random_sdlimit>0, "Error sdthreshold > 0, not " << token);
//       }
//       else if( startswith("seed1=",token) )
//       {
// 	 token += 6;
//          m_random_seed[0] = atoi(token);
//       }
//       else if( startswith("seed2=",token) )
//       {
// 	 token += 6;
//          m_random_seed[1] = atoi(token);
//       }
//       else if( startswith("seed3=",token) )
//       {
// 	 token += 6;
//          m_random_seed[2] = atoi(token);
//       }
//       else
//       {
// 	 badOption("randomize", token);
//       }
//       token = strtok(NULL, " \t");
//    }
//    if( lengthscaleset && !lengthscalezset )
//       m_random_distz = m_random_dist;
//    //   if( !lengthscaleset && lengthscalezset )
//    //      m_random_dist = m_random_distz;
// }

//-----------------------------------------------------------------------
void EW::processRandomBlock(char* buffer)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("randomblock", token) == 0,
	       "ERROR: not a randomblock line: " << token);
   if( m_events_parallel )
   {
      if( proc_zero() )
      {
         std::cout << "WARNING: randomblock command does not work together with parallel seismic events" << std::endl;
         std::cout << "randomblock command will be ignored " << std::endl;
      }
      return;
   }
   token = strtok(NULL, " \t");
   bool lengthscaleset=false, lengthscalezset=false, vsmaxset=false;
   float_sw4 corrlen=1000, corrlenz=1000, sigma=0.1, hurst=0.3, zmin=-1e38, zmax=1e38, vsmax=1e38;
   unsigned int seed=0;

   m_randomize = true;
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("corrlen=",token) )
      {
	 token += 8;
	 corrlen = atof(token);
	 CHECK_INPUT(corrlen>0, "Error randomblock, corrlen must be > 0, not " << token);
	 lengthscaleset = true;
      }
      else if( startswith("corrlenz=",token) )
      {
	 token += 9;
	 corrlenz = atof(token);
	 CHECK_INPUT(corrlenz>0, "Error randomblock, corrlenz must be > 0, not " << token);
	 lengthscalezset = true;
      }
      else if( startswith("sigma=",token) )
      {
	 token += 6;
	 sigma = atof(token);
	 CHECK_INPUT(sigma>0, "Error randomblock, sigma must be > 0, not " << token);
      }
      else if( startswith("hurst=",token) )
      {
	 token += 6;
	 hurst = atof(token);
      }
      else if( startswith("seed=",token) )
      {
	 token += 5;
         seed = atoi(token);
      }
      else if( startswith("zmin=",token) )
      {
	 token += 5;
         zmin = atof(token);
      }
      else if( startswith("zmax=",token) )
      {
	 token += 5;
         zmax = atof(token);
      }
      else if( startswith("vsmax=",token) )
      {
	 token += 6;
         vsmax = atof(token);
	 vsmaxset = true;
      }
      else
      {
	 badOption("randomblock", token);
      }
      token = strtok(NULL, " \t");
   }
   if( lengthscaleset && !lengthscalezset )
      corrlenz = corrlen;
   RandomizedMaterial* mtrl = new RandomizedMaterial( this, zmin, zmax, corrlen, 
						      corrlenz, hurst, sigma, seed );
   if( vsmaxset )
      mtrl->set_vsmax(vsmax);
   
   m_random_blocks.push_back(mtrl);
   //   if( !lengthscaleset && lengthscalezset )
   //      m_random_dist = m_random_distz;
}

//-----------------------------------------------------------------------
void EW::processEvent( char* buffer, int enr )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("event", token) == 0,
	       "ERROR: not an event line: " << token);
   token = strtok(NULL, " \t");

   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("path=",token) )
      {
	 token += 5; // skip path=
	 string path = token;
	 path += '/';
	 mPath.push_back(path);
	 //	 mPath[enr] = token;
	 //	 mPath[enr] += '/';
      }
      else if (startswith("obspath=", token))
      {
	 token += 8; // skip obspath=
	 string path = token;
	 path += '/';
	 mObsPath.push_back(path);
	 //	 mObsPath[enr] = token;
	 //	 mObsPath[enr] += '/';
      }
      else if( startswith("name=",token) )
      {
	 token += 5;
	 map<string,int>::iterator it=m_event_names.find(token);
	 CHECK_INPUT(it == m_event_names.end(), "ERROR: processEvent, name = " << token << " multiply defined");
	 m_event_names[token]=enr;
      }
      else if( startswith("parallel=",token) )
      {
         token += 9;
         std::string p=token;
         m_events_parallel = (p =="1" || p == "yes" || p=="on") || m_events_parallel;
      }
      else
      {
	 badOption("randomblock", token);
      }
      token = strtok(NULL, " \t");
   }
}

//-----------------------------------------------------------------------
int EW::findNumberOfEvents()
{
   char buffer[256];
   ifstream inputFile;
   MPI_Barrier(MPI_COMM_WORLD);
   inputFile.open(mName.c_str());
   if (!inputFile.is_open())
   {
      if (m_myRank == 0)
	 cerr << endl << "ERROR OPENING INPUT FILE: " << mName << endl << endl;
      CHECK_INPUT(false,"ERROR opening input file : " << mName << endl << endl);
   }
   int events=0;
   while (!inputFile.eof())
   {
      inputFile.getline(buffer, 256);
      if( startswith("event",buffer ) )
      {
	 processEvent( buffer, events );
	 events++;
      }
   }
   inputFile.close();
   if( events == 0 )
   {
      // use path variables from fileio, and set a default name
      m_event_names["Default event"]=0;
   }
   else
   {
      CHECK_INPUT( events > 0, "ERROR: no events found for optimization, nevent= " << events );

   }
   MPI_Barrier(MPI_COMM_WORLD); // Make sure all processes close the file before continue.
   return events;
}
