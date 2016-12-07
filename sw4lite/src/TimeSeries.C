// -*-c++-*-
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
#include <mpi.h>

#include <iostream>
#include <sstream>
#include <cstdlib>

#include "TimeSeries.h"
//#include "mpi.h"
#include "sacsubc.h"
//#include "csstime.h"

#include "Require.h"
#include "Filter.h"

#include "EW.h"

using namespace std;

//void parsedate( char* datestr, int& year, int& month, int& day, int& hour, int& minute,
//		int& second, int& msecond, int& fail );

TimeSeries::TimeSeries( EW* a_ew, std::string fileName, std::string staName, receiverMode mode, bool sacFormat, bool usgsFormat, 
			float_sw4 x, float_sw4 y, float_sw4 depth, bool topoDepth, int writeEvery, bool xyzcomponent ):
  m_ew(a_ew),
  m_mode(mode),
  m_nComp(0),
  m_myPoint(false),
  m_fileName(fileName),
  m_staName(staName),
  m_path(a_ew->mPath),
  mX(x),
  mY(y),
  mZ(depth),
  mGPX(0.0),
  mGPY(0.0),
  mGPZ(0.0),
  m_zRelativeToTopography(topoDepth),
  m_zTopo(0.0),
  mWriteEvery(writeEvery),
  m_usgsFormat(usgsFormat),
  m_sacFormat(sacFormat),
  m_xyzcomponent(xyzcomponent),
  m_i0(-999),
  m_j0(-999),
  m_k0(-999),
  m_grid0(-999),
  m_t0(0.0),
  m_shift(0.0),
  m_dt(1.0),
  mAllocatedSize(-1),
  mLastTimeStep(-1),
  mRecordedSol(NULL),
  mRecordedFloats(NULL),
  //  mIgnore(true), // are we still using this flag???No.
  //  mEventYear(2012),
  //  mEventMonth(2),
  //  mEventDay(8),
  //  mEventHour(10),
  //  mEventMinute(28),
  //  mEventSecond(0.0),
  m_rec_lat(38.0),
  m_rec_lon(-122.5),
  m_epi_lat(38.0),
  m_epi_lon(-122.5),
  m_epi_depth(0.0),
  m_epi_time_offset(0.0),
  m_x_azimuth(0.0),
  mBinaryMode(true)
  //  m_utc_set(false),
  //  m_utc_offset_computed(false),
  //  m_use_win(false),
  //  m_use_x(true),
  //  m_use_y(true),
  //  m_use_z(true),
  //  mQuietMode(false),
  //  m_compute_scalefactor(true)
{

// preliminary determination of nearest grid point ( before topodepth correction to mZ )
// Defines m_i0, m_j0, m_grid0, and m_k0 when the grid is not curvilinear
   a_ew->computeNearestGridPoint(m_i0, m_j0, m_k0, m_grid0, mX, mY, mZ);

// does this processor write this station?
   m_myPoint = a_ew->interior_point_in_proc(m_i0, m_j0, m_grid0);

// The following is a safety check to make sure only one processor writes each time series.
   int iwrite = m_myPoint ? 1 : 0;
   int counter;
   MPI_Allreduce( &iwrite, &counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
   REQUIRE2(counter == 1,"Exactly one processor must be writing each SAC, but counter = " << counter <<
	    " for receiver station " << m_fileName );
   if (!m_myPoint)
      return;

// from here on this processor writes this sac station and knows about its topography
   // Check that station is not above free surface
   if ( !a_ew->topographyExists() && (mZ < -1.0e-9 ) )
   {
      printf("Ignoring SAC station %s mX=%g, mY=%g, mZ=%g, because it is above the topography z=%g\n", 
	     m_fileName.c_str(),  mX,  mY, mZ, 0.0);
      m_myPoint=false;
      return;
   }

   if( a_ew->topographyExists() && (m_grid0 == a_ew->mNumberOfGrids-1 || m_zRelativeToTopography) )
   {
// Topography exists and point is located in the curvilinear grid, or point needs be corrected for topodepth.
//             Need to find z-coordinate of topography to check that the station is below it.
//             The z-coordinate of topography is also needed for topodepth correction
// 1. Evaluate z-coordinate of topography
      find_topo_zcoord( mX, mY, m_zTopo );
// 2. If location was specified with topodepth, correct z-level  
      if (m_zRelativeToTopography)
      {
	 mZ += m_zTopo;
	 m_zRelativeToTopography = false; // set to false so the correction isn't repeated (e.g. by the copy function)
      }
// 3. Make sure the station is below the topography, allow for a small roundoff (z is positive downwards)
      if ( mZ < m_zTopo - 1.0e-9 )
      {
	 printf("Ignoring SAC station %s mX=%g, mY=%g, mZ=%g, because it is above the topography z=%g\n", 
		m_fileName.c_str(),  mX,  mY, mZ, m_zTopo);
	 m_myPoint=false;
	 return;
      }
// 4. Recompute m_grid0 and k0 with corrected mZ
      if( computeNearestGridPointZ( mX, mY, mZ, m_grid0, m_k0 ) )
      {
	 cerr << "Can't invert curvilinear grid mapping for recevier station" << m_fileName << " mX= " << mX << " mY= " 
	      << mY << " mZ= " << mZ << endl;
	 cerr << "Placing the station on the surface (depth=0)." << endl;
	 m_grid0 = m_ew->mNumberOfGrids-1;
	 m_k0 = 1;

      }
   }
   
// actual location of station (nearest grid point)
   float_sw4 xG, yG, zG;
   xG = (m_i0-1)*a_ew->mGridSize[m_grid0];
   yG = (m_j0-1)*a_ew->mGridSize[m_grid0];
   if (m_grid0 < a_ew->mNumberOfCartesianGrids)
   {
      zG = a_ew->m_zmin[m_grid0] + (m_k0-1)*a_ew->mGridSize[m_grid0];
   }
   else
   {
      zG = a_ew->mZ(m_i0, m_j0, m_k0);
   }
   
// remember corrected location
   mGPX = xG;
   mGPY = yG;
   mGPZ = zG;

   if ( a_ew->getVerbosity()>=2 )
   {
      cout << "Receiver INFO for station " << m_fileName << ":" << endl <<
      "     initial location (x,y,z) = " << mX << " " << mY << " " << mZ << " zTopo= " << m_zTopo << endl <<
      "     nearest grid point (x,y,z) = " << mGPX << " " << mGPY << " " << mGPZ << " h= " << a_ew->mGridSize[m_grid0] << 
      " with indices (i,j,k)= " << m_i0 << " " << m_j0 << " " << m_k0 << " in grid " << m_grid0 << endl;
   }

// get number of components from m_mode
   if (m_mode == Displacement || m_mode == Velocity)
      m_nComp=3;
   else if (m_mode == Div)
      m_nComp=1;
   else if (m_mode == Curl)
      m_nComp=3;
   else if (m_mode == Strains)
      m_nComp=6;
   else if (m_mode == DisplacementGradient )
      m_nComp=9;

// allocate handles to solution array pointers
   mRecordedSol = new float_sw4*[m_nComp];
   for (int q=0; q<m_nComp; q++)
      mRecordedSol[q] = static_cast<float_sw4*>(0);

// keep a copy for saving on a sac file
   if (m_sacFormat)
   {
      mRecordedFloats = new float*[m_nComp];
      for (int q=0; q<m_nComp; q++)
	 mRecordedFloats[q] = static_cast<float*>(0);
   }
   else
      mRecordedFloats = static_cast<float**>(0);
  
// do some misc pre computations
   m_x_azimuth = a_ew->getGridAzimuth(); // degrees
  
   a_ew->computeGeographicCoord(mX, mY, m_rec_lon, m_rec_lat);
   a_ew->computeGeographicCoord(mGPX, mGPY, m_rec_gp_lon, m_rec_gp_lat);

   m_calpha = cos(M_PI*m_x_azimuth/180.0);
   m_salpha = sin(M_PI*m_x_azimuth/180.0);

   float_sw4 cphi   = cos(M_PI*m_rec_lat/180.0);
   float_sw4 sphi   = sin(M_PI*m_rec_lat/180.0);

   float_sw4 metersperdegree = a_ew->getMetersPerDegree();

//
// NOTE: this calculation assumes a spheroidal mapping
//
   m_thxnrm = m_salpha + (mX*m_salpha+mY*m_calpha)/cphi/metersperdegree * (M_PI/180.0) * sphi * m_calpha;
   m_thynrm = m_calpha - (mX*m_salpha+mY*m_calpha)/cphi/metersperdegree * (M_PI/180.0) * sphi * m_salpha;
   float_sw4 nrm = sqrt( m_thxnrm*m_thxnrm + m_thynrm*m_thynrm );
   m_thxnrm /= nrm;
   m_thynrm /= nrm;

// Set station ref utc = simulation ref utc
// m_t0 = 0 is set by default above.
   a_ew->get_utc( m_utc );

} // end constructor

//--------------------------------------------------------------
TimeSeries::~TimeSeries()
{
// deallocate the recording arrays
   if (mRecordedSol)
   {
      for (int q=0; q<m_nComp; q++)
      {
	 if (mRecordedSol[q])
	    delete [] mRecordedSol[q];
      }
      delete [] mRecordedSol;
   }
   if (mRecordedFloats)
   {
      for (int q=0; q<m_nComp; q++)
      {
	 if (mRecordedFloats[q])
	    delete [] mRecordedFloats[q];
      }
      delete [] mRecordedFloats;
   }
}

//--------------------------------------------------------------
void TimeSeries::allocateRecordingArrays( int numberOfTimeSteps, float_sw4 startTime, float_sw4 timeStep )
{
   if (!m_myPoint) return; // only one processor saves each time series
   if (numberOfTimeSteps > 0)
   {
       mAllocatedSize = numberOfTimeSteps+1;
       mLastTimeStep = -1;
       for (int q=0; q<m_nComp; q++)
       {
	  if (mRecordedSol[q]) delete [] mRecordedSol[q];
	  mRecordedSol[q] = new float_sw4[mAllocatedSize];
       }
       // Sac uses float only
       if (m_sacFormat)
       {
	  for (int q=0; q<m_nComp; q++)
	  {
	     if (mRecordedFloats[q]) delete [] mRecordedFloats[q];
	     mRecordedFloats[q] = new float[mAllocatedSize];
	  }
       }
   }
  // Move this time series to always start at 'startTime'. Perhaps this should be done elsewhere ?
   m_shift = startTime-m_t0;
   m_dt = timeStep;
}

//-----------------------------------------------------------------------
bool TimeSeries::find_topo_zcoord( float_sw4 X, float_sw4 Y, float_sw4& Ztopo )
{
   float_sw4 q, r;
   bool success = true;
   if (m_ew->topographyExists())
   {
      float_sw4 h = m_ew->mGridSize[m_ew->mNumberOfGrids-1];
      q = X/h + 1.0;
      r = Y/h + 1.0;
// evaluate elevation of topography on the grid
      if (!m_ew->interpolate_topography(q, r, Ztopo, true))
      {
	 cerr << "Unable to evaluate topography for receiver station" << m_fileName << " mX= " << mX << " mY= " << mY << endl;
	 cerr << "Setting topography to ZERO" << endl;
	 Ztopo = 0;
	 success = false;
      }
   }
   else
   {
      Ztopo = 0; // no topography
   }
   return success;
}

//-----------------------------------------------------------------------
bool TimeSeries::computeNearestGridPointZ( float_sw4 X, float_sw4 Y, float_sw4 Z,
					   int& grid0, int& k0 )
{
   bool success=true;
   int i0, j0;
   m_ew->computeNearestGridPoint(i0, j0, k0, grid0, X, Y, Z);
   if( grid0 == m_ew->mNumberOfGrids-1 && m_ew->topographyExists() )
   {
      float_sw4 q, r, s;
   // This only works if (x,y,z) is in my processor
      success = m_ew->invert_grid_mapping( grid0, X, Y, Z, q, r, s );
      if( success )
      {
	 k0 = static_cast<int>(round(s));
	 // Limit k0 to be inside the domain:
	 if( k0 < m_ew->m_kStartInt[m_grid0] )
	    k0 = m_ew->m_kStartInt[m_grid0] ;
	 if( k0 > m_ew->m_kEndInt[m_grid0] )
	    k0 = m_ew->m_kEndInt[m_grid0];
      }
   }
   return success;
}

//--------------------------------------------------------------
void TimeSeries::recordData(vector<float_sw4> & u) 
{
   if (!m_myPoint) return;

// better pass the right amount of data!
   if (u.size() != m_nComp)
   {
     printf("Error: TimeSeries::recordData: passing a vector of size=%i but nComp=%i\n", (int) u.size(), m_nComp);
     return;
   }
   
// ---------------------------------------------------------------
// This routine only knows how to push the nComp doubles on the array stack.
// The calling routine need to figure out what needs to be saved
// and do any necessary pre-calculations
// ---------------------------------------------------------------

   mLastTimeStep++;
   if (mLastTimeStep < mAllocatedSize)
   {
      //      if( m_xyzcomponent || (m_nComp != 3) )
      //      {
	 for (int q=0; q<m_nComp; q++)
	    mRecordedSol[q][mLastTimeStep] = u[q];
	 if (m_sacFormat)
	 {
	    for (int q=0; q<m_nComp; q++)
	       mRecordedFloats[q][mLastTimeStep] = (float) u[q];
	 }
// AP: The transformation to east-north-up components is now done just before the file is written
//      }
	 //      else
	 //      {
	 //// Transform to North-South, East-West, and Up components
	 //	 double uns = m_thynrm*u[0]-m_thxnrm*u[1];
	 //	 double uew = m_salpha*u[0]+m_calpha*u[1];
	 //         mRecordedSol[0][mLastTimeStep] = uew;
	 //         mRecordedSol[1][mLastTimeStep] = uns;
	 //         mRecordedSol[2][mLastTimeStep] =-u[2];
	 //	 if( m_sacFormat )
	 //	 {
	 //	    mRecordedFloats[0][mLastTimeStep] = static_cast<float>(uew);
	 //	    mRecordedFloats[1][mLastTimeStep] = static_cast<float>(uns);
	 //	    mRecordedFloats[2][mLastTimeStep] =-static_cast<float>(u[2]);
	 //	 }
	 //      }
   }
   else
   {
     printf("Ran out of recording space for the receiver station at (i,j,k,grid) = (%i, %i, %i, %i)\n",
	    m_i0, m_j0, m_k0, m_grid0);
     return;
   }
   if (mWriteEvery > 0 && mLastTimeStep > 0 && mLastTimeStep % mWriteEvery == 0)
      writeFile();
}

   
//----------------------------------------------------------------------
void TimeSeries::writeFile( string suffix )
{
  if (!m_myPoint) return;

// ------------------------------------------------------------------
// We should add an argument to this function that describes how the
// header and filename should be constructed
// ------------------------------------------------------------------

  stringstream filePrefix;

//building the file name...
  if( m_path != "." )
    filePrefix << m_path;
  if( suffix == "" )
     filePrefix << m_fileName << "." ;
  else
     filePrefix << m_fileName << suffix.c_str() << "." ;
  
// get the epicenter from EW object (note that the epicenter is not always known when this object is created)
// Right now, have no epicenter
  m_epi_lat = m_epi_lon=m_epi_depth=m_epi_time_offset = 0;
//  m_ew->get_epicenter( m_epi_lat, m_epi_lon, m_epi_depth, m_epi_time_offset );

  stringstream ux, uy, uz, uxy, uxz, uyz, uyx, uzx, uzy;
  
// Write out displacement components (ux, uy, uz)

  if( m_sacFormat )
  {
    string mode = "ASCII";
    if (mBinaryMode)
      mode = "BINARY";
    inihdr();

    stringstream msg;
    msg << "Writing " << mode << " SAC files, "
	<< "of size " << mLastTimeStep+1 << ": "
	<< filePrefix.str();

    string xfield, yfield, zfield, xyfield, xzfield, yzfield, yxfield, zxfield, zyfield;
     float azimx, azimy, updownang;
     if( m_mode == Displacement )
     {
	if( m_xyzcomponent )
	{
	   xfield = "X";
	   yfield = "Y";
	   zfield = "Z";
	   ux << filePrefix.str() << "x";
	   uy << filePrefix.str() << "y";
	   uz << filePrefix.str() << "z";
	   azimx = m_x_azimuth;
	   azimy = m_x_azimuth+90.;
	   updownang = 180;
	   msg << "[x|y|z]" << endl;
	}
	else
	{
 	   xfield = "EW";
 	   yfield = "NS";
 	   zfield = "UP";
 	   ux << filePrefix.str() << "e";
 	   uy << filePrefix.str() << "n";
 	   uz << filePrefix.str() << "u";
 	   azimx = 90.;// UX is east if !m_xycomponent
 	   azimy = 0.; // UY is north if !m_xycomponent
 	   updownang = 0;
 	   msg << "[e|n|u]" << endl;

	}
     }
     else if( m_mode == Velocity )
     {
        if( m_xyzcomponent )
	{
	   xfield = "Vx";
	   yfield = "Vy";
	   zfield = "Vz";
	   ux << filePrefix.str() << "xv";
	   uy << filePrefix.str() << "yv";
	   uz << filePrefix.str() << "zv";
	   azimx = m_x_azimuth;
	   azimy = m_x_azimuth+90.;
	   updownang = 180;
	   msg << "[xv|yv|zv]" << endl;
	}
	else
	{
 	   xfield = "Vew";
 	   yfield = "Vns";
 	   zfield = "Vup";
 	   ux << filePrefix.str() << "ev";
 	   uy << filePrefix.str() << "nv";
 	   uz << filePrefix.str() << "uv";
 	   azimx = 90.;// UX is east if !m_xycomponent
 	   azimy = 0.; // UY is north if !m_xycomponent
 	   updownang = 0;
 	   msg << "[ev|nv|uv]" << endl;
	}
     }
     else if( m_mode == Div )
     {
     	xfield = "Div";
     	ux << filePrefix.str() << "div";
	azimx = m_x_azimuth;
	azimy = m_x_azimuth+90.;
	updownang = 180;
     	msg << "[div]" << endl;
     }
     else if( m_mode == Curl )
     {
     	xfield = "Curlx";
     	yfield = "Curly";
     	zfield = "Curlz";
     	ux << filePrefix.str() << "curlx";
     	uy << filePrefix.str() << "curly";
     	uz << filePrefix.str() << "curlz";
	azimx = m_x_azimuth;
	azimy = m_x_azimuth+90.;
	updownang = 180;
     	msg << "[curlx|curly|curlz]" << endl;
     }
     else if( m_mode == Strains )
     {
     	xfield = "Uxx";
     	yfield = "Uyy";
     	zfield = "Uzz";
     	xyfield = "Uxy";
     	xzfield = "Uxz";
     	yzfield = "Uyz";
     	ux << filePrefix.str() << "xx";
     	uy << filePrefix.str() << "yy";
     	uz << filePrefix.str() << "zz";
     	uxy << filePrefix.str() << "xy";
     	uxz << filePrefix.str() << "xz";
     	uyz << filePrefix.str() << "yz";
	azimx = m_x_azimuth;
	azimy = m_x_azimuth+90.;
     	updownang = 180;
     	msg << "[xx|yy|zz|xy|xz|yz]" << endl;
     }
     else if( m_mode == DisplacementGradient )
     {
     	xfield  = "DUXDX";
     	xyfield = "DUXDY";
     	xzfield = "DUXDZ";

     	yxfield = "DUYDX";
     	yfield  = "DUYDY";
     	yzfield = "DUYDZ";

     	zxfield = "DUZDX";
     	zyfield = "DUZDY";
     	zfield  = "DUZDZ";

     	ux  << filePrefix.str() << "duxdx";
     	uxy << filePrefix.str() << "duxdy";
     	uxz << filePrefix.str() << "duxdz";

     	uyx << filePrefix.str() << "duydx";
     	uy << filePrefix.str()  << "duydy";
     	uyz << filePrefix.str() << "duydz";

     	uzx << filePrefix.str() << "duzdx";
     	uzy << filePrefix.str() << "duzdy";
     	uz << filePrefix.str()  << "duzdz";

	azimx = m_x_azimuth;
	azimy = m_x_azimuth+90.;
     	updownang = 180;
     	msg << "[duxdx|duxdy|duxdz|duydx|duydy|duydz|duzdx|duzdy|duzdz]" << endl;
     }
     // 	else if( !m_xycomponent && !m_velocities )
     // 	{
     // 	   xfield = "EW";
     // 	   yfield = "NS";
     // 	   zfield = "UP";
     // 	   ux << filePrefix.str() << "e";
     // 	   uy << filePrefix.str() << "n";
     // 	   uz << filePrefix.str() << "u";
     // 	   azimx = 90.;// UX is east if !m_xycomponent
     // 	   azimy = 0.; // UY is north if !m_xycomponent
     // 	   updownang = 0;
     // 	   msg << "[e|n|u]" << endl;
     // 	}
     // 	else if( !m_xycomponent && m_velocities )
     // 	{
     // 	   xfield = "Vew";
     // 	   yfield = "Vns";
     // 	   zfield = "Vup";
     // 	   ux << filePrefix.str() << "ev";
     // 	   uy << filePrefix.str() << "nv";
     // 	   uz << filePrefix.str() << "uv";
     // 	   azimx = 90.;// UX is east if !m_xycomponent
     // 	   azimy = 0.; // UY is north if !m_xycomponent
     // 	   updownang = 0;
     // 	   msg << "[ev|nv|uv]" << endl;
     // 	}
     // }
     // else if( m_div && m_velocities )
     // {
     // 	xfield = "VelDiv";
     // 	ux << filePrefix.str() << "vdiv";
     // 	azimx = a_ew->mGeoAz;
     // 	azimy = a_ew->mGeoAz+90.;
     // 	updownang = 180;
     // 	msg << "[vdiv]" << endl;
     // }
     // else if( m_curl && m_velocities && m_xycomponent )
     // {
     // 	xfield = "VelCurlx";
     // 	yfield = "VelCurly";
     // 	zfield = "VelCurlz";
     // 	ux << filePrefix.str() << "vcurlx";
     // 	uy << filePrefix.str() << "vcurly";
     // 	uz << filePrefix.str() << "vcurlz";
     // 	azimx = a_ew->mGeoAz;
     // 	azimy = a_ew->mGeoAz+90.;
     // 	updownang = 180;
     // 	msg << "[vcurlx|vcurly|vcurlz]" << endl;
     // }
     // else if( m_curl && !m_velocities && !m_xycomponent )
     // {
     // 	xfield = "CurlEW";
     // 	yfield = "CurlNS";
     // 	zfield = "CurlUP";
     // 	ux << filePrefix.str() << "curle";
     // 	uy << filePrefix.str() << "curln";
     // 	uz << filePrefix.str() << "curlu";
     // 	azimx = a_ew->mGeoAz;
     // 	azimy = a_ew->mGeoAz+90.;
     // 	updownang = 180;
     // 	msg << "[curle|curln|curlu]" << endl;
     // }
     // else if( m_curl && m_velocities && !m_xycomponent )
     // {
     // 	xfield = "VelCurlEW";
     // 	yfield = "VelCurlNS";
     // 	zfield = "VelCurlUP";
     // 	ux << filePrefix.str() << "vcurle";
     // 	uy << filePrefix.str() << "vcurln";
     // 	uz << filePrefix.str() << "vcurlu";
     // 	azimx = a_ew->mGeoAz;
     // 	azimy = a_ew->mGeoAz+90.;
     // 	updownang = 180;
     // 	msg << "[vcurle|vcurln|vcurlu]" << endl;
     // }
     // else if( m_strains && m_velocities )
     // {
     // 	xfield = "Velxx";
     // 	yfield = "Velyy";
     // 	zfield = "Velzz";
     // 	xyfield = "Velxy";
     // 	xzfield = "Velxz";
     // 	yzfield = "Velyz";
     // 	ux << filePrefix.str() << "vxx";
     // 	uy << filePrefix.str() << "vyy";
     // 	uz << filePrefix.str() << "vzz";
     // 	uxy << filePrefix.str() << "vxy";
     // 	uxz << filePrefix.str() << "vxz";
     // 	uyz << filePrefix.str() << "vyz";
     // 	azimx = a_ew->mGeoAz;
     // 	azimy = a_ew->mGeoAz+90.;
     // 	updownang = 180;
     // 	msg << "[vxx|vyy|vzz|vxy|vxz|vyz]" << endl;
     // }

// time to write the SAC files
     cout << msg.str();
     if (m_mode == Displacement || m_mode == Velocity || m_mode == Curl) // 3 components
     {
	if( m_xyzcomponent )
	{
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(ux.str().c_str()), 
			mRecordedFloats[0], (float) m_shift, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uy.str().c_str()), 
			mRecordedFloats[1], (float) m_shift, (float) m_dt,
			const_cast<char*>(yfield.c_str()), 90.0, azimy); 
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uz.str().c_str()), 
			mRecordedFloats[2], (float) m_shift, (float) m_dt,
			const_cast<char*>(zfield.c_str()), updownang, 0.0);
	}
	else
	{
           float** geographic = new float*[3];
	   geographic[0] = new float[mLastTimeStep+1];
	   geographic[1] = new float[mLastTimeStep+1];
	   geographic[2] = new float[mLastTimeStep+1];
	   for( int i=0 ; i <= mLastTimeStep ; i++ )
	   {
	      geographic[1][i] = m_thynrm*mRecordedFloats[0][i]-m_thxnrm*mRecordedFloats[1][i]; //ns
	      geographic[0][i] = m_salpha*mRecordedFloats[0][i]+m_calpha*mRecordedFloats[1][i]; //ew
	      geographic[2][i] = -mRecordedFloats[2][i];

	   }
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(ux.str().c_str()), 
			geographic[0], (float) m_shift, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uy.str().c_str()), 
			geographic[1], (float) m_shift, (float) m_dt,
			const_cast<char*>(yfield.c_str()), 90.0, azimy); 
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uz.str().c_str()), 
			geographic[2], (float) m_shift, (float) m_dt,
			const_cast<char*>(zfield.c_str()), updownang, 0.0);
           delete[] geographic[0];
           delete[] geographic[1];
           delete[] geographic[2];
	   delete[] geographic;
	}
     }
     else if (m_mode == Div) // 1 component
     {
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(ux.str().c_str()), 
			mRecordedFloats[0], (float) m_shift, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
     }
     else if (m_mode == Strains) // 6 components
     {
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(ux.str().c_str()), 
			mRecordedFloats[0], (float) m_shift, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uy.str().c_str()), 
			mRecordedFloats[1], (float) m_shift, (float) m_dt,
			const_cast<char*>(yfield.c_str()), 90.0, azimy); 
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uz.str().c_str()), 
			mRecordedFloats[2], (float) m_shift, (float) m_dt,
			const_cast<char*>(zfield.c_str()), updownang, 0.0); 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uxy.str().c_str()), 
			mRecordedFloats[3], (float) m_shift, (float) m_dt,
			const_cast<char*>(xyfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uxz.str().c_str()), 
			mRecordedFloats[4], (float) m_shift, (float) m_dt,
			const_cast<char*>(xzfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uyz.str().c_str()), 
			mRecordedFloats[5], (float) m_shift, (float) m_dt,
			const_cast<char*>(yzfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
     }
     else if (m_mode == DisplacementGradient ) // 9 components
     {
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(ux.str().c_str()), 
			mRecordedFloats[0], (float) m_shift, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uxy.str().c_str()), 
			mRecordedFloats[1], (float) m_shift, (float) m_dt,
			const_cast<char*>(xyfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uxz.str().c_str()), 
			mRecordedFloats[2], (float) m_shift, (float) m_dt,
			const_cast<char*>(xzfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 

       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uyx.str().c_str()), 
			mRecordedFloats[3], (float) m_shift, (float) m_dt,
			const_cast<char*>(yxfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uy.str().c_str()), 
			mRecordedFloats[4], (float) m_shift, (float) m_dt,
			const_cast<char*>(yfield.c_str()), 90.0, azimy); 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uyz.str().c_str()), 
			mRecordedFloats[5], (float) m_shift, (float) m_dt,
			const_cast<char*>(yzfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 

       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uzx.str().c_str()), 
			mRecordedFloats[6], (float) m_shift, (float) m_dt,
			const_cast<char*>(zxfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uzy.str().c_str()), 
			mRecordedFloats[7], (float) m_shift, (float) m_dt,
			const_cast<char*>(zyfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uz.str().c_str()), 
			mRecordedFloats[8], (float) m_shift, (float) m_dt,
			const_cast<char*>(zfield.c_str()), updownang, 0.0); 
     }

  } // end if m_sacFormat
  
  if( m_usgsFormat )
  {
    filePrefix << "txt";
    //    if (!m_ew->getQuiet())
      cout << "Writing ASCII USGS file, "
	   << "of size " << mLastTimeStep+1 << ": "
	   << filePrefix.str() << endl;

    write_usgs_format( filePrefix.str() );
  }

}

//-----------------------------------------------------------------------
void TimeSeries::
write_sac_format(int npts, char *ofile, float *y, float btime, float dt, char *var,
		 float cmpinc, float cmpaz)
{
  /*
    PURPOSE: SAVE RECEIVER DATA ON A SAC FILE
    
    	ofile	Char	name of file
    	y	R	array of values
    	npts	I	number of points in data
    	btime	R	start time
    	dt	R	sample interval
    	maxpts	I	maximum number of points to read
    	nerr	I	error return
    -----
  */
  float e;
  float depmax, depmin, depmen;
  int* nerr = 0;
// assign all names in a string array
//               0         1         2         3          4         5           6          7          8          9
  const char *nm[]={"DEPMAX", "DEPMIN", "DEPMEN", "NPTS    ","DELTA   ","B       ", "E       ","LEVEN   ","LOVROK  ","LCALDA  ",
//              10          11          12          13          14          15           16          17          18
	       "NZYEAR  ", "NZJDAY  ", "NZHOUR  ", "NZMIN   ", "NZSEC   ", "NZMSEC   ", "KCMPNM  ", "STLA    ", "STLO    ",
//              19          20          21          22          23          24          25
	       "EVLA    ", "EVLO    ", "EVDP    ", "O       ", "CMPINC  ", "CMPAZ   ", "KSTNM   "
  };

  newhdr();
  scmxmn(y,npts,&depmax,&depmin,&depmen);
//  setfhv("DEPMAX", depmax, nerr);
  setfhv( nm[0], depmax, nerr);
  setfhv( nm[1], depmin, nerr);
  setfhv( nm[2], depmen, nerr);
  setnhv( nm[3], npts,nerr);
  setfhv( nm[4], dt  ,nerr);
  setfhv( nm[5], btime  ,nerr);
  e = btime + (npts -1 )*dt;
  setfhv( nm[6], e, nerr);
  setlhv( nm[7], 1, nerr);
  setlhv( nm[8], 1, nerr);
  setlhv( nm[9], 1, nerr);

  // setup time info
  //  if( m_utc_set )
  //  {
     int days = 0;
     for( int m=1 ; m<m_utc[1]; m++ )
	days += lastofmonth(m_utc[0],m);
     days += m_utc[2];

     setnhv( nm[10], m_utc[0], nerr);
     setnhv( nm[11], days, nerr);
     setnhv( nm[12], m_utc[3], nerr);
     setnhv( nm[13], m_utc[4], nerr);
     setnhv( nm[14], m_utc[5], nerr);
     setnhv( nm[15], m_utc[6], nerr);
     //  }
     //  else
     //  {
     //     setnhv( nm[10], mEventYear, nerr);
     //     setnhv( nm[11], mEventDay, nerr);
     //     setnhv( nm[12], mEventHour, nerr);
     //     setnhv( nm[13], mEventMinute, nerr);
     //     setnhv( nm[14], static_cast<int>(mEventSecond), nerr);
     //     setnhv( nm[15], 0, nerr);
     //  }

  // field we're writing
  setkhv( nm[16], var, nerr);

  // location of the receiver
//   double lat, lon;
//   a_ew->computeGeographicCoord(mX, mY, lon, lat); //(C.B: I think that this is the point we want)
  setfhv( nm[17], m_rec_lat, nerr);
  setfhv( nm[18], m_rec_lon, nerr);
  // location of epicenter
  setfhv( nm[19], m_epi_lat, nerr);
  setfhv( nm[20], m_epi_lon, nerr);
  setfhv( nm[21], m_epi_depth/1000.0, nerr); // in km, not meters
  // time offset for epicenter source
  setfhv( nm[22], m_epi_time_offset, nerr);

  // set inclination and azimuthal angle
  setfhv( nm[23], cmpinc, nerr);
  setfhv( nm[24], cmpaz, nerr);

  // set the station name
  setkhv( nm[25], const_cast<char*>(m_staName.c_str()), nerr);


  if (!mBinaryMode)
    awsac(npts, ofile, y);
  else
    bwsac(npts, ofile, y);
}

//-----------------------------------------------------------------------
void TimeSeries::write_usgs_format(string a_fileName)
{
   string mname[] = {"Zero","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
   FILE *fd=fopen(a_fileName.c_str(),"w");
   float_sw4 lat, lon;
   float_sw4 x, y, z;

   if( fd == NULL )
      cout << "ERROR: opening USGS file " << a_fileName << " for writing" <<  endl;
   else
   {

      //      cout << "IN write_usgs " << " mt0 = " << m_t0 << " mshift = " << m_shift << endl;
// frequency resolution
//    double freq_limit=-999;
//    if (a_ew->m_prefilter_sources)
//      freq_limit = a_ew->m_fc;
//    else if (a_ew->m_limit_frequency)
//      freq_limit = a_ew->m_frequency_limit;

// write the header
   fprintf(fd, "# Author: SW4\n");
   fprintf(fd, "# Scenario: %s\n", "test"/*a_ew->m_scenario.c_str()*/);
//   if( m_utc_set )
// AP: micro-second field is padded from left with 0, i.e., 1 micro sec gets written as 001, second is also padded by a zero, if needed
   fprintf(fd, "# Date: UTC  %02i/%02i/%i:%i:%i:%02i.%.3i\n", m_utc[1], m_utc[2], m_utc[0], m_utc[3],
	   m_utc[4], m_utc[5], m_utc[6] );
      //   else
      //      fprintf(fd, "# Date: %i-%s-%i\n", mEventDay, mname[mEventMonth].c_str(), mEventYear);

   fprintf(fd, "# Bandwith (Hz): %e\n", 1.234 /*freq_limit*/);
   fprintf(fd, "# Station: %s\n", m_fileName.c_str() /*mStationName.c_str()*/ );
   fprintf(fd, "# Target location (WGS84 longitude, latitude) (deg): %e %e\n", m_rec_lon, m_rec_lat);
   fprintf(fd, "# Actual location (WGS84 longitude, latitude) (deg): %e %e\n", m_rec_gp_lon, m_rec_gp_lat);
// distance in horizontal plane
   fprintf(fd, "# Distance from target to actual location (m): %e\n", sqrt( (mX-mGPX)*(mX-mGPX)+(mY-mGPY)*(mY-mGPY) ) );
   fprintf(fd, "# nColumns: %i\n", m_nComp+1);
   
   fprintf(fd, "# Column 1: Time (s)\n");
   if (m_mode == Displacement && m_xyzcomponent )
   {
     fprintf(fd, "# Column 2: X displacement (m)\n");
     fprintf(fd, "# Column 3: Y displacement (m)\n");
     fprintf(fd, "# Column 4: Z displacement (m)\n");
   }
   else if (m_mode == Displacement && !m_xyzcomponent )
   {
     fprintf(fd, "# Column 2: East-west displacement (m)\n");
     fprintf(fd, "# Column 3: North-sourth displacement (m)\n");
     fprintf(fd, "# Column 4: Up-down displacement (m)\n");
   }
   else if( m_mode == Velocity && m_xyzcomponent )
   {
     fprintf(fd, "# Column 2: X velocity (m/s)\n");
     fprintf(fd, "# Column 3: Y velocity (m/s)\n");
     fprintf(fd, "# Column 4: Z velocity (m/s)\n");
   }
   else if( m_mode == Velocity && !m_xyzcomponent )
   {
     fprintf(fd, "# Column 2: East-west velocity (m/s)\n");
     fprintf(fd, "# Column 3: Nort-south velocity (m/s)\n");
     fprintf(fd, "# Column 4: Up-down velocity (m/s)\n");
   }
   else if( m_mode == Div )
   {
     fprintf(fd, "# Column 2: divergence of displacement ()\n");
   }
   else if( m_mode == Curl )
   {
     fprintf(fd, "# Column 2: curl of displacement, component 1 ()\n");
     fprintf(fd, "# Column 3: curl of displacement, component 2 ()\n");
     fprintf(fd, "# Column 4: curl of displacement, component 3 ()\n");
//       }
   }
   else if( m_mode == Strains )
   {
     fprintf(fd, "# Column 2: xx strain component ()\n");
     fprintf(fd, "# Column 3: yy strain component ()\n");
     fprintf(fd, "# Column 4: zz strain component ()\n");
     fprintf(fd, "# Column 5: xy strain component ()\n");
     fprintf(fd, "# Column 6: xz strain component ()\n");
     fprintf(fd, "# Column 7: yz strain component ()\n");
   }
   else if( m_mode == DisplacementGradient )
   {
     fprintf(fd, "# Column 2 : dux/dx component ()\n");
     fprintf(fd, "# Column 3 : dux/dy component ()\n");
     fprintf(fd, "# Column 4 : dux/dz component ()\n");
     fprintf(fd, "# Column 5 : duy/dx component ()\n");
     fprintf(fd, "# Column 6 : duy/dy component ()\n");
     fprintf(fd, "# Column 7 : duy/dz component ()\n");
     fprintf(fd, "# Column 8 : duz/dx component ()\n");
     fprintf(fd, "# Column 9 : duz/dy component ()\n");
     fprintf(fd, "# Column 10: duz/dz component ()\n");
   }

// write the data

   if( m_xyzcomponent || (!m_xyzcomponent && m_mode == Div) )
   {
      for( int i = 0 ; i <= mLastTimeStep ; i++ )
      {
	 fprintf(fd, "%e", m_shift + i*m_dt);
	 for (int q=0; q<m_nComp; q++)
// AP (not always enough resolution)	    fprintf(fd, " %20.12g", mRecordedSol[q][i]);
	    fprintf(fd, " %24.17e", mRecordedSol[q][i]);
	 fprintf(fd, "\n");
      }
   }
   else if( m_mode == Displacement || m_mode == Velocity )
   {
      for( int i = 0 ; i <= mLastTimeStep ; i++ )
      {
	 fprintf(fd, "%e", m_shift + i*m_dt);
	 float_sw4 uns = m_thynrm*mRecordedSol[0][i]-m_thxnrm*mRecordedSol[1][i];
	 float_sw4 uew = m_salpha*mRecordedSol[0][i]+m_calpha*mRecordedSol[1][i];
// AP (not always enough resolution)
	 // fprintf(fd, " %20.12g", uew );
	 // fprintf(fd, " %20.12g", uns );
	 // fprintf(fd, " %20.12g", -mRecordedSol[2][i] );
	 fprintf(fd, " %24.17e", uew );
	 fprintf(fd, " %24.17e", uns );
	 fprintf(fd, " %24.17e", -mRecordedSol[2][i] );
	 fprintf(fd, "\n");
      }
   }
   else
   {
      printf("TimeSeries::write_usgs_format, Can not write ");
      if( m_mode == Strains )
	 printf("strains");
      else if( m_mode == Curl )
	 printf("curl");
      else if( m_mode == DisplacementGradient )
	 printf("displacement gradient");
      printf(" in geographic coordinates\n" );
   }
   fclose(fd);
   }
}

//-----------------------------------------------------------------------
int TimeSeries::lastofmonth( int year, int month )
{
   int days;
   //   int leapyear=0;
   if( month == 2 )
   {
      int leapyear = ( year % 400 == 0 ) || ( (year % 4 == 0) && !(year % 100 == 0) );
      days = 28 + leapyear;
   }
   else if( month==4 || month==6 || month==9 || month==11 )
      days = 30;
   else
      days = 31;
   return days;
}
