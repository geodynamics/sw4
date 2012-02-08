// #  WPP LICENSE
// # ----------------------------------------------------------------------
// # WPP - Wave propagation Program
// # ----------------------------------------------------------------------
// # Copyright (C) 2011, Lawrence Livermore National Security, LLC.  
// # Produced at the Lawrence Livermore National Laboratory
// # 
// # Written by:
// # 
// # Bjorn Sjogreen   (sjogreen2@llnl.gov)
// # Anders Petersson  (andersp@llnl.gov)
// # 
// # Alums:
// # Stefan Nilsson      
// # Daniel Appelo
// # Kathleen McCandless (mccandless2@llnl.gov)
// # Caroline Bono
// # 
// # CODE-227123 All rights reserved.
// # 
// # This file is part of WPP, v2.1
// # 
// # Please also read docs/GPLLICENSE.txt which contains 
// # "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License as published by
// # the Free Software Foundation; version 2, dated June 1991.
// # 
// # This program is distributed in the hope that it will be useful,
// # but WITHOUT ANY WARRANTY; without even the implied warranty of
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// # terms and conditions of the GNU General Public License for more details.
// # 
// # You should have received a copy of the GNU General Public License along with
// # this program; if not, write to the Free Software Foundation, Inc., 59 Temple
// # Place, Suite 330, Boston MA 02111-1307 USA.
// # ----------------------------------------------------------------------
#include "mpi.h"

#include "SAC.h"
#include "sacsubc.h"
#include "csstime.h"

#include "Require.h"
// #include "Grid.h"
// #include "Geometry.h"

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "WPP2.h"

using namespace std;

// initializing static member
struct tm* SAC::mTimePtr = NULL;

// int  SAC::m_ghost_points = 1;
// int  SAC::m_padding      = 1;

// std::vector<double> SAC::m_gridSize;
// std::vector<double> SAC::m_zmin    ;

SAC::
SAC(WPP2 *wpp,
    string fileName,
    string stationName,
    int frequency,
    int writeEvery,
    bool binaryMode,
    bool momentMode):
  mHeaderName(fileName),
  mStationName(stationName),
  mMomentMode(momentMode),
  mDateSet(false), mTimeSet(false),
  m_writeSAC(false),
  mEpicenterInitialized(false),
  mEpicenterTimeOffset(0.0),
  mBinaryMode(binaryMode),
  mFrequency(frequency),
  mWriteEvery(writeEvery),
  mWPP(wpp),
  m_previous(0),
  mRestarted(false),
  mEventYear(0),
  mEventMonth(0),
  mEventDay(0),
  mEventHour(0),
  mEventMinute(0),
  mEventSecond(0.0),
  m_xycomponent(true),
  m_velocities(false),
  m_sacformat(true),
  m_usgsformat(false),
  m_zRelativeToTopography(false),
  mIgnore(false),
  m_div(false),
  m_curl(false),
  m_strains(false)
{
   if (mStationName == "")
      // default to file name
      mStationName = mHeaderName;
}

//-----------------------------------------------------------------------
void SAC::setBinaryMode(bool onoff)
{
  mBinaryMode = onoff;
}
  
//-----------------------------------------------------------------------
void SAC::set_formats( int usgs, int sac )
{
   m_sacformat = sac==1;
   m_usgsformat = usgs==1;
}

//-----------------------------------------------------------------------
void SAC::set_coordinate( double x, double y, double z )
{
   //   REQUIRE2(!mIOInitialized, "Error: Cannot set_coordinate after initialization");
   mX = x;
   mY = y;
   mZ = z;
}

//-----------------------------------------------------------------------
void SAC::initialize()
{
// preliminary determination of nearest grid point ( before topodepth correction to mZ)
   mWPP->computeNearestGridPoint(m_iSAC,m_jSAC,m_kSAC,m_gridSAC,mX,mY,mZ);

// does this processor write this station?
   m_writeSAC = mWPP->interior_point_in_proc(m_iSAC,m_jSAC,m_gridSAC);

// The following is a safety check to make sure only one processor writes each SAC record.
// We could remove this check if we are sure interior_point_in_proc() never lies
  int iwrite = m_writeSAC ? 1 : 0;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  std::vector<int> whoIsOne(size);
  int counter = 0;
  MPI_Allgather(&iwrite, 1, MPI_INT, &whoIsOne[0], 1, MPI_INT,MPI_COMM_WORLD);
  for (unsigned int i = 0; i < whoIsOne.size(); ++i)
    if (whoIsOne[i] == 1)
      {
        counter++;
      }
  REQUIRE2(counter == 1,"Exactly one processor must be writing each SAC, but counter = " << counter <<
  "for SAC station" << mHeaderName );

  if (!m_writeSAC) return;

// from here on this processor writes this sac station and knows about its topography

  double zMinTilde, q, r, s;
  if (mWPP->topographyExists())
  {
    int gCurv = mWPP->mNumberOfGrids - 1;
    double h = mWPP->mGridSize[gCurv];
    q = mX/h + 1.0;
    r = mY/h + 1.0;
// evaluate elevation of topography on the grid
    if (!mWPP->interpolate_topography(q, r, zMinTilde, true))
    {
      cerr << "Unable to evaluate topography for SAC station" << mHeaderName << " mX= " << mX << " mY= " << mY << endl;
      cerr << "Setting topography to ZERO" << endl;
      zMinTilde = 0;
    }
  }
  else
  {
    zMinTilde = 0; // no topography
  }

// if location was specified with topodepth, correct z-level  
   if (m_zRelativeToTopography)
   {
     mZ += zMinTilde;
   }
   
   double zMin = zMinTilde - 1.e-9; // allow for a little roundoff
   
// make sure the station is below the topography (z is positive downwards)
   if ( mZ < zMin)
   {
     mIgnore = true;
     printf("Ignoring SAC station %s mX=%g, mY=%g, mZ=%g, because it is above the topography z=%g\n", 
     mHeaderName.c_str(),  mX,  mY, mZ, zMinTilde);
// don't write this station
     m_writeSAC=false;
     return;
   }
     
// now we can find the closest grid point  
   mWPP->computeNearestGridPoint(m_iSAC,m_jSAC,m_kSAC,m_gridSAC,mX,mY,mZ);

   if( m_gridSAC == mWPP->mNumberOfGrids-1 && mWPP->topographyExists() )
   {
// Curvilinear
     bool canBeInverted = mWPP->invert_curvilinear_grid_mapping( mX, mY, mZ, q, r, s );
     if (mWPP->invert_curvilinear_grid_mapping( mX, mY, mZ, q, r, s )) // the inversion was successful
     {
       m_kSAC = (int)floor(s);
       if (s-(m_kSAC+0.5) > 0.) m_kSAC++;
       m_kSAC = max(1, m_kSAC);
       int Nz = mWPP->m_kEnd[m_gridSAC] - mWPP->m_ghost_points;
       m_kSAC = min(Nz, m_kSAC);
     }
     else
     {
       cerr << "Can't invert curvilinear grid mapping for SAC station" << mHeaderName << " mX= " << mX << " mY= " 
	    << mY << " mZ= " << mZ << endl;
       cerr << "Placing the station on the surface (depth=0)." << endl;
       m_kSAC = 1;
     }
   }
   
// actual location of station (nearest grid point)
   double xG, yG, zG;
   xG = (m_iSAC-1)*mWPP->mGridSize[m_gridSAC];
   yG = (m_jSAC-1)*mWPP->mGridSize[m_gridSAC];
   if (m_gridSAC < mWPP->mNumberOfCartesianGrids)
   {
     zG = mWPP->m_zmin[m_gridSAC] + (m_kSAC-1)*mWPP->mGridSize[m_gridSAC];
   }
   else
   {
     zG = mWPP->mZ(m_iSAC, m_jSAC, m_kSAC);
   }
   
    if (mWPP->getVerbosity()>=1)
    {
      cout << "SAC info for station " << mHeaderName << ":" << 
	" initial location (x,y,z) = " << mX << " " << mY << " " << mZ << 
	" moved to nearest grid point (x,y,z) = " << xG << " " << yG << " " << zG << 
	" h= " << mWPP->mGridSize[m_gridSAC] <<
	" with indices (i,j,k)= " << m_iSAC << " " << m_jSAC << " " << m_kSAC << " in grid " << m_gridSAC << endl;
    }
    

// correct location
   mX = xG;
   mY = yG;
   mZ = zG;
   
// do some misc pre computations
   mWPP->computeGeographicCoord(mX, mY, m_lon, m_lat);

   m_calpha = cos(M_PI*mWPP->mGeoAz/180.0);
   m_salpha = sin(M_PI*mWPP->mGeoAz/180.0);

   double cphi   = cos(M_PI*m_lat/180.0);
   double sphi   = sin(M_PI*m_lat/180.0);
   double metersperdegree = mWPP->mMetersPerDegree;

   m_thxnrm = m_salpha + (mX*m_salpha+mY*m_calpha)/cphi/metersperdegree * (M_PI/180.0) * sphi * m_calpha;
   m_thynrm = m_calpha - (mX*m_salpha+mY*m_calpha)/cphi/metersperdegree * (M_PI/180.0) * sphi * m_salpha;
   double nrm = sqrt( m_thxnrm*m_thxnrm + m_thynrm*m_thynrm );
   m_thxnrm /= nrm;
   m_thynrm /= nrm;
   m_dthi = 1.0/(2*mWPP->mDt);

}

//-----------------------------------------------------------------------
void SAC::set_nsew( )
{
   m_xycomponent = false;
}

//-----------------------------------------------------------------------
void SAC::set_div()
{
   m_div = true;
   m_curl = false;
   m_strains = false;
}

//-----------------------------------------------------------------------
void SAC::set_curl()
{
   m_div = false;
   m_curl = true;
   m_strains = false;
}

//-----------------------------------------------------------------------
void SAC::set_velocities( )
{
   m_velocities = true;
}

//-----------------------------------------------------------------------
void SAC::set_strains()
{
   m_strains    = true;
   m_div        = false;
   m_curl       = false;
   m_sacformat  = false;
   m_usgsformat = true;
   m_xycomponent= true;
}

//-----------------------------------------------------------------------
void
SAC::
initializeEpicenter(const GeographicCoord& sourceLocation, double time)
{
  mEpicenter = sourceLocation;
  mEpicenterTimeOffset = time;
  mEpicenterInitialized = true;

//   // ---------------------------------------------------------
//   // Set up time information to date the SAC files in 
//   // the header.
//   //
//   // 1. Gather time information into class variables, using
//   //    values set by user (eventDate/Time) or the default
//   //    which is system time.
//   // 2. Convert this representation to seconds
//   // 3. Subtract off the epicenter time.
//   // 4. Store values (day/month/year) back into class variables
//   // ---------------------------------------------------------

//   if (!mRestarted)
//   {
//      // Fill in the datetime struct to convert to seconds...
//      struct date_time dateTime, dateTime2;
//      if (mDateSet)
//      {
//         dateTime.year = mEventYear;
//         dateTime.month = mEventMonth;
//         dateTime.day = mEventDay;
//      }
//      else
//      {
//         dateTime.year = 1900 + mTimePtr->tm_year;
//         dateTime.month = mTimePtr->tm_mon+1;
//         dateTime.day = mTimePtr->tm_mday;
//      }
//      if (mTimeSet)
//      {
//         dateTime.hour = mEventHour;
//         dateTime.minute = mEventMinute;
//         dateTime.second = mEventSecond;
//      }
//      else 
//      {
//         dateTime.hour = mTimePtr->tm_hour;
//         dateTime.minute = mTimePtr->tm_min;
//         dateTime.second = mTimePtr->tm_sec;
//      }
     
//      // Now convert to seconds...
//      mdtodate(&dateTime); /* fills in date field */
//      htoe(&dateTime); /* uses date field */
     
//      // Subtract off timeoffset
//      dateTime.epoch -= mEpicenterTimeOffset;
  
//      // Now back...
//      etoh(&dateTime);

//      mEventYear = dateTime.year;
//      mEventMonth = dateTime.month;
//      mEventDay = dateTime.doy;
//      mEventHour = dateTime.hour;
//      mEventMinute = dateTime.minute;
//      mEventSecond = dateTime.second;
//   }
}

void SAC::recordData(int cycle ) 
{
//   REQUIRE2(mPathSet && mEpicenterInitialized, 
// 	   "Error: Must call initialize before recordData");
  
   if (!m_writeSAC) return;

   bool curvilinear = mWPP->topographyExists() && m_gridSAC == mWPP->mNumberOfGrids-1;

    // Time to record data for sac file
   if( !m_div && !m_curl && !m_strains )
   {
      if( m_xycomponent )
      {
	// Cartesian components
	 mRecordedUX.push_back((float)(mWPP->mU)[m_gridSAC](1,m_iSAC,m_jSAC,m_kSAC));
	 mRecordedUY.push_back((float)(mWPP->mU)[m_gridSAC](2,m_iSAC,m_jSAC,m_kSAC));
	 mRecordedUZ.push_back((float)(mWPP->mU)[m_gridSAC](3,m_iSAC,m_jSAC,m_kSAC));
      }
      else
      {
	// North-South, East-West, and Up components
	 double ux  = (mWPP->mU)[m_gridSAC](1,m_iSAC,m_jSAC,m_kSAC);
	 double uy  = (mWPP->mU)[m_gridSAC](2,m_iSAC,m_jSAC,m_kSAC);

	 double uns = m_thynrm*ux-m_thxnrm*uy;
	 double uew = m_salpha*ux+m_calpha*uy;
	 
	 mRecordedUX.push_back( (float)(uew) ); //E-W is stored in UX
	 mRecordedUY.push_back( (float)(uns) ); //N-S is stored in UY
	 mRecordedUZ.push_back(-(float)(mWPP->mU)[m_gridSAC](3,m_iSAC,m_jSAC,m_kSAC));
      }
   }
   // For curl and div, don't worry about one sided formulas at boundaries, because of ghost points.
   else if( m_div && !curvilinear )
   {

      int i=m_iSAC, j=m_jSAC, k=m_kSAC, g=m_gridSAC;
      double factor = 1.0/(2*mWPP->mGridSize[g]);
      mRecordedUX.push_back((mWPP->mU[g](1,i+1, j, k) - mWPP->mU[g](1,i-1, j, k)+
			     mWPP->mU[g](2,i, j+1, k) - mWPP->mU[g](2,i, j-1, k)+
			     mWPP->mU[g](3,i, j, k+1) - mWPP->mU[g](3,i, j, k-1))*factor);
      // Record dummy variables ==> can keep the velocity computation below unchanged.
      mRecordedUY.push_back(0);
      mRecordedUZ.push_back(0);
   }
   else if( m_div && curvilinear )
   {
      int i=m_iSAC, j=m_jSAC, k=m_kSAC, g=m_gridSAC;
      double factor = 1.0/(2);
      mRecordedUX.push_back( ( mWPP->mQ(1,i,j,k)*(mWPP->mU[g](1,i+1,j,k) - mWPP->mU[g](1,i-1,j,k))+
			       mWPP->mQ(2,i,j,k)*(mWPP->mU[g](2,i+1,j,k) - mWPP->mU[g](2,i-1,j,k))+
			       mWPP->mQ(3,i,j,k)*(mWPP->mU[g](3,i+1,j,k) - mWPP->mU[g](3,i-1,j,k))+
			       mWPP->mR(1,i,j,k)*(mWPP->mU[g](1,i,j+1,k) - mWPP->mU[g](1,i,j-1,k))+
			       mWPP->mR(2,i,j,k)*(mWPP->mU[g](2,i,j+1,k) - mWPP->mU[g](2,i,j-1,k))+
			       mWPP->mR(3,i,j,k)*(mWPP->mU[g](3,i,j+1,k) - mWPP->mU[g](3,i,j-1,k))+
			       mWPP->mS(1,i,j,k)*(mWPP->mU[g](1,i,j,k+1) - mWPP->mU[g](1,i,j,k-1))+
			       mWPP->mS(2,i,j,k)*(mWPP->mU[g](2,i,j,k+1) - mWPP->mU[g](2,i,j,k-1))+
			       mWPP->mS(3,i,j,k)*(mWPP->mU[g](3,i,j,k+1) - mWPP->mU[g](3,i,j,k-1))  )*factor);
      // Record dummy variables ==> can keep the velocity computation below unchanged.
      mRecordedUY.push_back(0);
      mRecordedUZ.push_back(0);
   }
   else if( m_curl && !curvilinear )
   {
      int i=m_iSAC, j=m_jSAC, k=m_kSAC, g=m_gridSAC;
      double factor = 1.0/(2*mWPP->mGridSize[g]);
      double duydx = (mWPP->mU[g](2,i+1,j,k) - mWPP->mU[g](2,i-1,j,k))*factor;
      double duzdx = (mWPP->mU[g](3,i+1,j,k) - mWPP->mU[g](3,i-1,j,k))*factor;
      double duxdy = (mWPP->mU[g](1,i,j+1,k) - mWPP->mU[g](1,i,j-1,k))*factor;
      double duzdy = (mWPP->mU[g](3,i,j+1,k) - mWPP->mU[g](3,i,j-1,k))*factor;
      double duxdz = (mWPP->mU[g](1,i,j,k+1) - mWPP->mU[g](1,i,j,k-1))*factor;
      double duydz = (mWPP->mU[g](2,i,j,k+1) - mWPP->mU[g](2,i,j,k-1))*factor;
      if( m_xycomponent )
      {
	 mRecordedUX.push_back( duzdy-duydz );
	 mRecordedUY.push_back( duxdz-duzdx );
	 mRecordedUZ.push_back( duydx-duxdy );
      }
      else
      {
	 double uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
	 double uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
	 mRecordedUX.push_back( uew );
	 mRecordedUY.push_back( uns );
	 mRecordedUZ.push_back( -(duydx-duxdy) );
      }
   }
   else if( m_curl && curvilinear )
   {
      int i=m_iSAC, j=m_jSAC, k=m_kSAC, g=m_gridSAC;
      double factor = 1.0/(2);
      double duxdq = (mWPP->mU[g](1,i+1,j,k) - mWPP->mU[g](1,i-1,j,k));
      double duydq = (mWPP->mU[g](2,i+1,j,k) - mWPP->mU[g](2,i-1,j,k));
      double duzdq = (mWPP->mU[g](3,i+1,j,k) - mWPP->mU[g](3,i-1,j,k));
      double duxdr = (mWPP->mU[g](1,i,j+1,k) - mWPP->mU[g](1,i,j-1,k));
      double duydr = (mWPP->mU[g](2,i,j+1,k) - mWPP->mU[g](2,i,j-1,k));
      double duzdr = (mWPP->mU[g](3,i,j+1,k) - mWPP->mU[g](3,i,j-1,k));
      double duxds = (mWPP->mU[g](1,i,j,k+1) - mWPP->mU[g](1,i,j,k-1));
      double duyds = (mWPP->mU[g](2,i,j,k+1) - mWPP->mU[g](2,i,j,k-1));
      double duzds = (mWPP->mU[g](3,i,j,k+1) - mWPP->mU[g](3,i,j,k-1));
      double duzdy = mWPP->mQ(2,i,j,k)*duzdq+mWPP->mR(2,i,j,k)*duzdr+mWPP->mS(2,i,j,k)*duzds;
      double duydz = mWPP->mQ(3,i,j,k)*duydq+mWPP->mR(3,i,j,k)*duydr+mWPP->mS(3,i,j,k)*duyds;
      double duxdz = mWPP->mQ(3,i,j,k)*duxdq+mWPP->mR(3,i,j,k)*duxdr+mWPP->mS(3,i,j,k)*duxds;
      double duzdx = mWPP->mQ(1,i,j,k)*duzdq+mWPP->mR(1,i,j,k)*duzdr+mWPP->mS(1,i,j,k)*duzds;
      double duydx = mWPP->mQ(1,i,j,k)*duydq+mWPP->mR(1,i,j,k)*duydr+mWPP->mS(1,i,j,k)*duyds;
      double duxdy = mWPP->mQ(2,i,j,k)*duxdq+mWPP->mR(2,i,j,k)*duxdr+mWPP->mS(2,i,j,k)*duxds;
      if( m_xycomponent )
      {
	 mRecordedUX.push_back( (duzdy-duydz)*factor );
	 mRecordedUY.push_back( (duxdz-duzdx)*factor );
	 mRecordedUZ.push_back( (duydx-duxdy)*factor );
      }
      else
      {
	 double uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
	 double uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
	 mRecordedUX.push_back( uew*factor );
	 mRecordedUY.push_back( uns*factor );
	 mRecordedUZ.push_back( -(duydx-duxdy)*factor );
      }
   }
   else if( m_strains && !curvilinear )
   {
      int i=m_iSAC, j=m_jSAC, k=m_kSAC, g=m_gridSAC;
      double factor = 1.0/(2*mWPP->mGridSize[g]);
      double duydx = (mWPP->mU[g](2,i+1,j,k) - mWPP->mU[g](2,i-1,j,k))*factor;
      double duzdx = (mWPP->mU[g](3,i+1,j,k) - mWPP->mU[g](3,i-1,j,k))*factor;
      double duxdy = (mWPP->mU[g](1,i,j+1,k) - mWPP->mU[g](1,i,j-1,k))*factor;
      double duzdy = (mWPP->mU[g](3,i,j+1,k) - mWPP->mU[g](3,i,j-1,k))*factor;
      double duxdz = (mWPP->mU[g](1,i,j,k+1) - mWPP->mU[g](1,i,j,k-1))*factor;
      double duydz = (mWPP->mU[g](2,i,j,k+1) - mWPP->mU[g](2,i,j,k-1))*factor;
      double duxdx = (mWPP->mU[g](1,i+1,j,k) - mWPP->mU[g](1,i-1,j,k))*factor;
      double duydy = (mWPP->mU[g](2,i,j+1,k) - mWPP->mU[g](2,i,j-1,k))*factor;
      double duzdz = (mWPP->mU[g](3,i,j,k+1) - mWPP->mU[g](3,i,j,k-1))*factor;
      mRecordedUX.push_back( duxdx );
      mRecordedUY.push_back( duydy );
      mRecordedUZ.push_back( duzdz );
      mRecordedUXY.push_back( 0.5*(duydx+duxdy) );
      mRecordedUXZ.push_back( 0.5*(duzdx+duxdz) );
      mRecordedUYZ.push_back( 0.5*(duydz+duzdy) );
   }
   else if( m_strains && curvilinear )
   {
      int i=m_iSAC, j=m_jSAC, k=m_kSAC, g=m_gridSAC;
      double factor = 1.0/(2);
      double duxdq = (mWPP->mU[g](1,i+1,j,k) - mWPP->mU[g](1,i-1,j,k));
      double duydq = (mWPP->mU[g](2,i+1,j,k) - mWPP->mU[g](2,i-1,j,k));
      double duzdq = (mWPP->mU[g](3,i+1,j,k) - mWPP->mU[g](3,i-1,j,k));
      double duxdr = (mWPP->mU[g](1,i,j+1,k) - mWPP->mU[g](1,i,j-1,k));
      double duydr = (mWPP->mU[g](2,i,j+1,k) - mWPP->mU[g](2,i,j-1,k));
      double duzdr = (mWPP->mU[g](3,i,j+1,k) - mWPP->mU[g](3,i,j-1,k));
      double duxds = (mWPP->mU[g](1,i,j,k+1) - mWPP->mU[g](1,i,j,k-1));
      double duyds = (mWPP->mU[g](2,i,j,k+1) - mWPP->mU[g](2,i,j,k-1));
      double duzds = (mWPP->mU[g](3,i,j,k+1) - mWPP->mU[g](3,i,j,k-1));
      double duzdy = (mWPP->mQ(2,i,j,k)*duzdq+mWPP->mR(2,i,j,k)*duzdr+mWPP->mS(2,i,j,k)*duzds)*factor;
      double duydz = (mWPP->mQ(3,i,j,k)*duydq+mWPP->mR(3,i,j,k)*duydr+mWPP->mS(3,i,j,k)*duyds)*factor;
      double duxdz = (mWPP->mQ(3,i,j,k)*duxdq+mWPP->mR(3,i,j,k)*duxdr+mWPP->mS(3,i,j,k)*duxds)*factor;
      double duzdx = (mWPP->mQ(1,i,j,k)*duzdq+mWPP->mR(1,i,j,k)*duzdr+mWPP->mS(1,i,j,k)*duzds)*factor;
      double duydx = (mWPP->mQ(1,i,j,k)*duydq+mWPP->mR(1,i,j,k)*duydr+mWPP->mS(1,i,j,k)*duyds)*factor;
      double duxdy = (mWPP->mQ(2,i,j,k)*duxdq+mWPP->mR(2,i,j,k)*duxdr+mWPP->mS(2,i,j,k)*duxds)*factor;
      double duxdx = ( mWPP->mQ(1,i,j,k)*(mWPP->mU[g](1,i+1,j,k) - mWPP->mU[g](1,i-1,j,k))+
		       mWPP->mR(1,i,j,k)*(mWPP->mU[g](1,i,j+1,k) - mWPP->mU[g](1,i,j-1,k))+
		       mWPP->mS(1,i,j,k)*(mWPP->mU[g](1,i,j,k+1) - mWPP->mU[g](1,i,j,k-1)) )*factor;
      double duydy = ( mWPP->mQ(2,i,j,k)*(mWPP->mU[g](2,i+1,j,k) - mWPP->mU[g](2,i-1,j,k))+
		       mWPP->mR(2,i,j,k)*(mWPP->mU[g](2,i,j+1,k) - mWPP->mU[g](2,i,j-1,k))+
		       mWPP->mS(2,i,j,k)*(mWPP->mU[g](2,i,j,k+1) - mWPP->mU[g](2,i,j,k-1)) )*factor;
      double duzdz = ( mWPP->mQ(3,i,j,k)*(mWPP->mU[g](3,i+1,j,k) - mWPP->mU[g](3,i-1,j,k))+
		       mWPP->mR(3,i,j,k)*(mWPP->mU[g](3,i,j+1,k) - mWPP->mU[g](3,i,j-1,k))+
		       mWPP->mS(3,i,j,k)*(mWPP->mU[g](3,i,j,k+1) - mWPP->mU[g](3,i,j,k-1)) )*factor;
      mRecordedUX.push_back( duxdx );
      mRecordedUY.push_back( duydy );
      mRecordedUZ.push_back( duzdz );
      mRecordedUXY.push_back( 0.5*(duydx+duxdy) );
      mRecordedUXZ.push_back( 0.5*(duzdx+duxdz) );
      mRecordedUYZ.push_back( 0.5*(duydz+duzdy) );
   }

   if( m_velocities )
   {
     int n = mRecordedUX.size();
     double vx, vy, vz;
     if( n == 1 )
     {
	m_dmx = mRecordedUX[n-1];
	m_dmy = mRecordedUY[n-1];
	m_dmz = mRecordedUZ[n-1];
     }
     else if( n == 2 )
     {
	m_d0x = mRecordedUX[n-1];
	m_d0y = mRecordedUY[n-1];
	m_d0z = mRecordedUZ[n-1];
     }
     else if( n == 3 )
     {
// one sided formula for v(1)
        vx = (-3*m_dmx+4*m_d0x-mRecordedUX[n-1] )*m_dthi; 
        mRecordedUX[n-3] = vx;
        vy = (-3*m_dmy+4*m_d0y-mRecordedUY[n-1] )*m_dthi; 
        mRecordedUY[n-3] = vy;
        vz = (-3*m_dmz+4*m_d0z-mRecordedUZ[n-1] )*m_dthi;
        mRecordedUZ[n-3] = vz;
     }
     if( n >= 3 )
     {
// Standard case;
        vx = (mRecordedUX[n-1]-m_dmx)*m_dthi;
        vy = (mRecordedUY[n-1]-m_dmy)*m_dthi;
        vz = (mRecordedUZ[n-1]-m_dmz)*m_dthi;
        mRecordedUX[n-2] = vx;
        mRecordedUY[n-2] = vy;
        mRecordedUZ[n-2] = vz;

// Need one sided forward, in case we are at the last point.
// Otherwise, this is overwritten in the next time step.
        vx = (m_dmx-4*m_d0x+3*mRecordedUX[n-1])*m_dthi;
        m_dmx = m_d0x;
	m_d0x = mRecordedUX[n-1];
	mRecordedUX[n-1] = vx;
        vy = (m_dmy-4*m_d0y+3*mRecordedUY[n-1])*m_dthi;
        m_dmy = m_d0y;
	m_d0y = mRecordedUY[n-1];
	mRecordedUY[n-1] = vy;
        vz = (m_dmz-4*m_d0z+3*mRecordedUZ[n-1])*m_dthi;
        m_dmz = m_d0z;
	m_d0z = mRecordedUZ[n-1];
	mRecordedUZ[n-1] = vz;
     }
   }
   if( m_velocities && m_strains )
   {
     int n = mRecordedUXY.size();
     double vx, vy, vz;
     if( n == 1 )
     {
	m_dmxy = mRecordedUXY[n-1];
	m_dmxz = mRecordedUXZ[n-1];
	m_dmyz = mRecordedUYZ[n-1];
     }
     else if( n == 2 )
     {
	m_d0xy = mRecordedUXY[n-1];
	m_d0xz = mRecordedUXZ[n-1];
	m_d0yz = mRecordedUYZ[n-1];
     }
     else if( n == 3 )
     {
// one sided formula for v(1)
        vx = (-3*m_dmxy+4*m_d0xy-mRecordedUXY[n-1] )*m_dthi; 
        mRecordedUXY[n-3] = vx;
        vy = (-3*m_dmxz+4*m_d0xz-mRecordedUXZ[n-1] )*m_dthi; 
        mRecordedUXZ[n-3] = vy;
        vz = (-3*m_dmyz+4*m_d0yz-mRecordedUYZ[n-1] )*m_dthi;
        mRecordedUYZ[n-3] = vz;
     }
     if( n >= 3 )
     {
// Standard case;
        vx = (mRecordedUXY[n-1]-m_dmxy)*m_dthi;
        vy = (mRecordedUXZ[n-1]-m_dmxz)*m_dthi;
        vz = (mRecordedUYZ[n-1]-m_dmyz)*m_dthi;
        mRecordedUXY[n-2] = vx;
        mRecordedUXZ[n-2] = vy;
        mRecordedUYZ[n-2] = vz;

// Need one sided forward, in case we are at the last point.
// Otherwise, this is overwritten in the next time step.
        vx = (m_dmxy-4*m_d0xy+3*mRecordedUXY[n-1])*m_dthi;
        m_dmxy = m_d0xy;
	m_d0xy = mRecordedUXY[n-1];
	mRecordedUXY[n-1] = vx;
        vy = (m_dmxz-4*m_d0xz+3*mRecordedUXZ[n-1])*m_dthi;
        m_dmxz = m_d0xz;
	m_d0xz = mRecordedUXZ[n-1];
	mRecordedUYZ[n-1] = vy;
        vz = (m_dmyz-4*m_d0yz+3*mRecordedUYZ[n-1])*m_dthi;
        m_dmyz = m_d0yz;
	m_d0yz = mRecordedUYZ[n-1];
	mRecordedUYZ[n-1] = vz;
     }
   }
   if (mWriteEvery != 0 && cycle % mWriteEvery == 0 && 
      cycle >= mWriteEvery) writeFile();
}
   
void SAC::writeFile()
{
//   REQUIRE2(mPathSet && mEpicenterInitialized, 
// 	   "Error: Must call initialize before writeFile");
//   REQUIRE2(mTimingDataSet, "Error: Must set timing data before writeFile");
  
  if (!m_writeSAC) return;

  stringstream filePrefix;
  if (mWPP->mPath != ".") filePrefix << mWPP->mPath;

//building the file name...
  filePrefix << mHeaderName << "." ;//<< userN << ".";
  
  stringstream ux, uy, uz, uxy, uxz, uyz;
  
// Write out displacement components (ux, uy, uz)

  if( m_sacformat )
  {
     string mode = "ASCII";
     if (mBinaryMode)
	mode = "BINARY";
     inihdr();

     stringstream msg;
     msg << "Writing " << mode << " SAC Files, "
	 << "of size " << mRecordedUX.size() << ": "
	 << filePrefix.str();

     string xfield, yfield, zfield, xyfield, xzfield, yzfield;
     double azimx, azimy, updownang;
     if( !m_div && !m_curl && !m_strains )
     {
	if( m_xycomponent && !m_velocities )
	{
	   xfield = "X";
	   yfield = "Y";
	   zfield = "Z";
	   ux << filePrefix.str() << "x";
	   uy << filePrefix.str() << "y";
	   uz << filePrefix.str() << "z";
	   azimx = mWPP->mGeoAz;
	   azimy = mWPP->mGeoAz+90.;
	   updownang = 180;
	   msg << "[x|y|z]" << endl;
	}
	else if( m_xycomponent && m_velocities )
	{
	   xfield = "Vx";
	   yfield = "Vy";
	   zfield = "Vz";
	   ux << filePrefix.str() << "xv";
	   uy << filePrefix.str() << "yv";
	   uz << filePrefix.str() << "zv";
	   azimx = mWPP->mGeoAz;
	   azimy = mWPP->mGeoAz+90.;
	   updownang = 180;
	   msg << "[xv|yv|zv]" << endl;
	}
	else if( !m_xycomponent && !m_velocities )
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
	else if( !m_xycomponent && m_velocities )
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
     else if( m_div && !m_velocities )
     {
	xfield = "Div";
	ux << filePrefix.str() << "div";
	azimx = mWPP->mGeoAz;
	azimy = mWPP->mGeoAz+90.;
	updownang = 180;
	msg << "[div]" << endl;
     }
     else if( m_div && m_velocities )
     {
	xfield = "VelDiv";
	ux << filePrefix.str() << "vdiv";
	azimx = mWPP->mGeoAz;
	azimy = mWPP->mGeoAz+90.;
	updownang = 180;
	msg << "[vdiv]" << endl;
     }
     else if( m_curl && !m_velocities && m_xycomponent )
     {
	xfield = "Curlx";
	yfield = "Curly";
	zfield = "Curlz";
	ux << filePrefix.str() << "curlx";
	uy << filePrefix.str() << "curly";
	uz << filePrefix.str() << "curlz";
	azimx = mWPP->mGeoAz;
	azimy = mWPP->mGeoAz+90.;
	updownang = 180;
	msg << "[curlx|curly|curlz]" << endl;
     }
     else if( m_curl && m_velocities && m_xycomponent )
     {
	xfield = "VelCurlx";
	yfield = "VelCurly";
	zfield = "VelCurlz";
	ux << filePrefix.str() << "vcurlx";
	uy << filePrefix.str() << "vcurly";
	uz << filePrefix.str() << "vcurlz";
	azimx = mWPP->mGeoAz;
	azimy = mWPP->mGeoAz+90.;
	updownang = 180;
	msg << "[vcurlx|vcurly|vcurlz]" << endl;
     }
     else if( m_curl && !m_velocities && !m_xycomponent )
     {
	xfield = "CurlEW";
	yfield = "CurlNS";
	zfield = "CurlUP";
	ux << filePrefix.str() << "curle";
	uy << filePrefix.str() << "curln";
	uz << filePrefix.str() << "curlu";
	azimx = mWPP->mGeoAz;
	azimy = mWPP->mGeoAz+90.;
	updownang = 180;
	msg << "[curle|curln|curlu]" << endl;
     }
     else if( m_curl && m_velocities && !m_xycomponent )
     {
	xfield = "VelCurlEW";
	yfield = "VelCurlNS";
	zfield = "VelCurlUP";
	ux << filePrefix.str() << "vcurle";
	uy << filePrefix.str() << "vcurln";
	uz << filePrefix.str() << "vcurlu";
	azimx = mWPP->mGeoAz;
	azimy = mWPP->mGeoAz+90.;
	updownang = 180;
	msg << "[vcurle|vcurln|vcurlu]" << endl;
     }
     else if( m_strains && m_velocities )
     {
	xfield = "Velxx";
	yfield = "Velyy";
	zfield = "Velzz";
	xyfield = "Velxy";
	xzfield = "Velxz";
	yzfield = "Velyz";
	ux << filePrefix.str() << "vxx";
	uy << filePrefix.str() << "vyy";
	uz << filePrefix.str() << "vzz";
	uxy << filePrefix.str() << "vxy";
	uxz << filePrefix.str() << "vxz";
	uyz << filePrefix.str() << "vyz";
	azimx = mWPP->mGeoAz;
	azimy = mWPP->mGeoAz+90.;
	updownang = 180;
	msg << "[vxx|vyy|vzz|vxy|vxz|vyz]" << endl;
     }
     else if( m_strains )
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
	azimx = mWPP->mGeoAz;
	azimy = mWPP->mGeoAz+90.;
	updownang = 180;
	msg << "[xx|yy|zz|xy|xz|yz]" << endl;
     }

     cout << msg.str() ;
     writeSACFile(mRecordedUX.size(), 
               const_cast<char*>(ux.str().c_str()), 
               &mRecordedUX[0], (float)(mWPP->mTstart - mWPP->m_t0Shift), (float)(mWPP->mDt),
		  const_cast<char*>(xfield.c_str()), 90.0, azimx); 
     if( !m_div )
     {
	writeSACFile(mRecordedUZ.size(), 
               const_cast<char*>(uz.str().c_str()), 
               &mRecordedUZ[0], (float)(mWPP->mTstart - mWPP->m_t0Shift), (float)(mWPP->mDt),
               const_cast<char*>(zfield.c_str()), updownang, 0.0);      
	writeSACFile(mRecordedUY.size(), 
               const_cast<char*>(uy.str().c_str()), 
               &mRecordedUY[0], (float)(mWPP->mTstart - mWPP->m_t0Shift), (float)(mWPP->mDt),
               const_cast<char*>(yfield.c_str()), 90.0, azimy); 
     }
     if( m_strains )
     {
	writeSACFile(mRecordedUXY.size(), 
		     const_cast<char*>(uxy.str().c_str()), 
		     &mRecordedUXY[0], (float)(mWPP->mTstart - mWPP->m_t0Shift), (float)(mWPP->mDt),
		  const_cast<char*>(xyfield.c_str()), 90.0, azimx); 
	writeSACFile(mRecordedUXZ.size(), 
		     const_cast<char*>(uxz.str().c_str()), 
		     &mRecordedUXZ[0], (float)(mWPP->mTstart - mWPP->m_t0Shift), (float)(mWPP->mDt),
		  const_cast<char*>(xzfield.c_str()), 90.0, azimx); 
	writeSACFile(mRecordedUYZ.size(), 
		     const_cast<char*>(uyz.str().c_str()), 
		     &mRecordedUYZ[0], (float)(mWPP->mTstart - mWPP->m_t0Shift), (float)(mWPP->mDt),
		  const_cast<char*>(yzfield.c_str()), 90.0, azimx); 

     }

  }
  if( m_usgsformat )
  {
    if (mWPP->getVerbosity()>0)
    {
     cout << "Writing USGS receiever data file, "
	 << "of size " << mRecordedUX.size() << ": "
	  << filePrefix.str() << "txt" << endl;
    }
    filePrefix << "txt";
    write_usgs_format( filePrefix.str() );
  }

}

void 
SAC::
writeSACFile(int npts, char *ofile, float *y, float btime, float dt, char *var,
             float cmpinc, float cmpaz)
{
  /*
    PURPOSE: WRITE AN EVENLY SPACED SAC FILE
    
    	write a binary sac file with evenly sampled data
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
  string nm[]={"DEPMAX", "DEPMIN", "DEPMEN", "NPTS    ","DELTA   ","B       ", "E       ","LEVEN   ","LOVROK  ","LCALDA  ",
//              10          11          12          13          14          15           16          17          18
	       "NZYEAR  ", "NZJDAY  ", "NZHOUR  ", "NZMIN   ", "NZSEC   ", "NZMSEC   ", "KCMPNM  ", "STLA    ", "STLO    ",
//              19          20          21          22          23          24          25
	       "EVLA    ", "EVLO    ", "EVDP    ", "O       ", "CMPINC  ", "CMPAZ   ", "KSTNM   "
  };

  newhdr();
  scmxmn(y,npts,&depmax,&depmin,&depmen);
//  setfhv("DEPMAX", depmax, nerr);
  setfhv((char *) nm[0].c_str(), depmax, nerr);
  setfhv((char *) nm[1].c_str(), depmin, nerr);
  setfhv((char *) nm[2].c_str(), depmen, nerr);
  setnhv((char *) nm[3].c_str(), npts,nerr);
  setfhv((char *) nm[4].c_str(), dt  ,nerr);
  setfhv((char *) nm[5].c_str(), btime  ,nerr);
  e = btime + (npts -1 )*dt;
  setfhv((char *) nm[6].c_str(), e, nerr);
  setlhv((char *) nm[7].c_str(), 1, nerr);
  setlhv((char *) nm[8].c_str(), 1, nerr);
  setlhv((char *) nm[9].c_str(), 1, nerr);

  // setup time info
  setnhv((char *) nm[10].c_str(), mEventYear, nerr);
  setnhv((char *) nm[11].c_str(), mEventDay, nerr);
  setnhv((char *) nm[12].c_str(), mEventHour, nerr);
  setnhv((char *) nm[13].c_str(), mEventMinute, nerr);
  setnhv((char *) nm[14].c_str(), static_cast<int>(mEventSecond), nerr);
  setnhv((char *) nm[15].c_str(), 0, nerr);

  // field we're writing
  setkhv((char *) nm[16].c_str(), var, nerr);

  // location of the receiver
//   double lat, lon;
//   mWPP->computeGeographicCoord(mX, mY, lon, lat); //(C.B: I think that this is the point we want)
  setfhv((char *) nm[17].c_str(), m_lat, nerr);
  setfhv((char *) nm[18].c_str(), m_lon, nerr);
  // location of epicenter
  setfhv((char *) nm[19].c_str(), mEpicenter.getLatitude(), nerr);
  setfhv((char *) nm[20].c_str(), mEpicenter.getLongitude(), nerr);
  setfhv((char *) nm[21].c_str(), mEpicenter.getDepth()/1000.0, nerr);
  // time offset for epicenter source
  setfhv((char *) nm[22].c_str(), mEpicenterTimeOffset, nerr);

  // set inclination and azimuthal angle
  setfhv((char *) nm[23].c_str(), cmpinc, nerr);
  setfhv((char *) nm[24].c_str(), cmpaz, nerr);

  // set the station name
  setkhv((char *) nm[25].c_str(), const_cast<char*>(mStationName.c_str()), nerr);


  if (!mBinaryMode)
    awsac(npts, ofile, y);
  else
    bwsac(npts, ofile, y);
}

void 
SAC::
setEventDate(std::string date)
{
  string err = "SAC command error: ";
  string delimiter = "/";

  //---------------------------------------------------
  // parse the year
  //---------------------------------------------------
  string::size_type pos1 = date.find(delimiter, 0);
  VERIFY2(pos1 != string::npos,
	  err << "eventDate format must be YYYY/MM/DD, no '/' found,"
	  << date);
  VERIFY2(pos1 == 4,
	  err << "eventDate format must be YYYY/MM/DD,"
	  << date);
  string yearstring(date, 0, pos1);
  mEventYear = atoi(yearstring.c_str());
  VERIFY2(mEventYear >= 0,
	  err << "eventDate year must be positive, not: "
	  << mEventYear);
  
  //---------------------------------------------------
  // parse the month
  //---------------------------------------------------
  VERIFY2(date.length() > pos1,
	  err << "eventDate format must be YYYY/MM/DD, no 'MM/DD' found, "
	  << date);

  string::size_type pos2 = date.find(delimiter, pos1+1);
  VERIFY2(pos2 != string::npos,
	  err << "eventDate format must be YYYY/MM/DD, no 'MM/DD' found, "
	  << date);
  VERIFY2(pos2 == 7,
	  err << "eventDate format must be YYYY/MM/DD, "
	  << date);

  string monthstring(date, pos1+1, pos2-1);
  mEventMonth = atoi(monthstring.c_str());
  VERIFY2(mEventMonth > 0 && mEventMonth < 13,
	  err << "eventDate month must be between 1 and 12, not: "
	  << mEventMonth);
  
  //---------------------------------------------------
  // parse the day
  //---------------------------------------------------
  string daystring(date, pos2+1, date.length());
  mEventDay = atoi(daystring.c_str());
  VERIFY2(mEventDay > 0 && mEventDay < 32,
	  err << "eventDate day must be between 1 and 31, not: "
	  << mEventDay);

  mDateSet = true;
}
 
void 
SAC::
setEventTime(std::string time)
{
  string err = "SAC command error: ";
  string delimiter = ":";

  //---------------------------------------------------
  // parse the hour
  //---------------------------------------------------
  string::size_type pos1 = time.find(delimiter, 0);
  VERIFY2(pos1 != string::npos,
	  err << "eventTime format must be HH:MM:SS, no ':' found,"
	  << time);
  string hourstring(time, 0, pos1);
  mEventHour = atoi(hourstring.c_str());
  VERIFY2(mEventHour >= 0 && mEventHour < 24,
	  err << "eventTime hour must be between 0 and 24 hours, not: "
	  << mEventHour);

  //---------------------------------------------------
  // parse the minutes
  //---------------------------------------------------
  VERIFY2(time.length() > pos1,
	  err << "eventTime format must be HH:MM:SS, no 'MM:SS' found, "
	  << time);

  string::size_type pos2 = time.find(delimiter, pos1+1);
  VERIFY2(pos2 != string::npos,
	  err << "eventTime format must be HH:MM:SS, no 'MM:SS' found, "
	  << time);

  string minutestring(time, pos1+1, pos2-1);
  mEventMinute = atoi(minutestring.c_str());
  VERIFY2(mEventMinute >= 0 && mEventMinute < 60,
	  err << "eventTime minute must be between 0 and 59 minutes, not: "
	  << mEventMinute);
  
  //---------------------------------------------------
  // parse the seconds
  //---------------------------------------------------
  string secondstring(time, pos2+1, time.length());
  mEventSecond = atoi(secondstring.c_str());
  VERIFY2(mEventSecond >=0 && mEventSecond < 60,
	  err << "eventTime second must be between 0 and 59 seconds, not: "
	  << mEventSecond);

  mTimeSet = true;
}


// void SAC::rescale( std::vector<double> invs, double xc, double yc, double zc )
// {
//   if (!m_writeSAC) return;
//   double h = mGrid->getGridSpacing();
//   blitz::TinyVector<double,3> coord = mGrid->getCoordFromGlobalIndex(mGlobalN,mGlobalM,mGlobalL);
//   double x = coord[0];
//   double y = coord[1];
//   double z = coord[2];
//   int now  = mRecordedUX.size()-1;
//   double k = 1.0/(now-m_previous);
//   double m = -m_previous/(now-m_previous);
//   for( int t=m_previous+1 ; t <= now ; t++ )
//   {
//      double ipfact = k*t+m;
//      mRecordedUX[t] += ipfact*(-invs[0] -invs[3]*(y-yc)-invs[4]*(z-zc));
//      mRecordedUY[t] += ipfact*(-invs[1] +invs[3]*(x-xc)               -invs[5]*(z-zc));
//      mRecordedUZ[t] += ipfact*(-invs[2]                +invs[4]*(x-xc)+invs[5]*(y-yc));
//   }
// //  m_previous = now+1;
// // AP reordered recording and projection steps so that 'now' corresponds to 'up' which was recorded before
// // the invariants were removed. Next time around, m_previous should therefore point to 'now'
//   m_previous = now;
// }

void 
SAC::
initializeSystemTime(tm* tptr)
{
   mTimePtr = tptr;
}


void SAC::write_usgs_format( string fname )
{
//   char* mname[12] ={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
   string mname[12] ={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
//   ofstream outfile(fname.c_str(), ios_base::out);
   FILE *fd=fopen(fname.c_str(),"w");
   double lat, lon;
   double x, y, z;
   mWPP->coord_from_index( m_iSAC, m_jSAC, m_kSAC, m_gridSAC, x, y, z );
   mWPP->computeGeographicCoord( x, y, lon, lat );

// tmp
//    cout << "write USGS format, SAC station: " << mHeaderName << " mX= " << mX << " mY= " << mY << " mZ= " << mZ << endl
//  	<< " m_gridSAC= " << m_gridSAC << " m_iSAC= " << m_iSAC << " m_jSAC= " << m_jSAC << " m_kSAC= " << m_kSAC << endl
//  	<< " x= " << x << " y= " << y << " z= " << z << endl;

   double freq_limit=-999;
   if (mWPP->m_prefilter_sources)
     freq_limit = mWPP->m_fc;
   else if (mWPP->m_limit_frequency)
     freq_limit = mWPP->m_frequency_limit;

//   outfile << "# Author: WPP" << endl;
//   outfile << "# Scenario: " << mWPP->m_scenario << endl;
//   outfile << "# Date: " << mTimePtr->tm_mday << "-" << mname[mTimePtr->tm_mon] <<"-"<< 1900+mTimePtr->tm_year << endl;
//   outfile << "# Bandwith (Hz): " << freq_limit << endl;
//   outfile << "# Station: " << mStationName << endl;
//   outfile << "# Target location (WGS84 longitude, latitude) (deg): " << m_lon << " " << m_lat << endl;
//   outfile << "# Actual location (WGS84 longitude, latitude) (deg): " << lon   << " " << lat << endl;
//   outfile << "# Distance from target to actual location (m): " 
//	   << sqrt( (x-mX)*(x-mX)+(y-mY)*(y-mY) ) << endl; // distance in horizontal plane
//   outfile << "# Column 1: Time (s)" << endl;
   fprintf(fd, "# Author: WPP\n");
   fprintf(fd, "# Scenario: %s\n", mWPP->m_scenario.c_str());
   fprintf(fd, "# Date: %i-%s-%i\n", mTimePtr->tm_mday, mname[mTimePtr->tm_mon].c_str(), 1900+mTimePtr->tm_year);
   fprintf(fd, "# Bandwith (Hz): %e\n", freq_limit);
   fprintf(fd, "# Station: %s\n", mStationName.c_str());
   fprintf(fd, "# Target location (WGS84 longitude, latitude) (deg): %e %e\n", m_lon, m_lat);
   fprintf(fd, "# Actual location (WGS84 longitude, latitude) (deg): %e %e\n", lon, lat);
// distance in horizontal plane
   fprintf(fd, "# Distance from target to actual location (m): %e\n", sqrt( (x-mX)*(x-mX)+(y-mY)*(y-mY) ));
   if( m_div )
      fprintf(fd, "# nColumns: 2\n");
   else if( m_strains )
      fprintf(fd, "# nColumns: 7\n");
   else
      fprintf(fd, "# nColumns: 4\n");
   
   fprintf(fd, "# Column 1: Time (s)\n");
   if( !m_div && !m_curl && !m_strains )
   {
      if( !m_xycomponent && m_velocities )
      {
//       outfile << "# Column 2: East-west velocity (m/s)" << endl;
//       outfile << "# Column 3: North-south velocity (m/s)" << endl;
//       outfile << "# Column 4: Up-down velocity (m/s)" << endl;
	 fprintf(fd, "# Column 2: East-west velocity (m/s)\n");
	 fprintf(fd, "# Column 3: North-south velocity (m/s)\n");
	 fprintf(fd, "# Column 4: Up-down velocity (m/s)\n");
      }
      else if( !m_xycomponent && !m_velocities )
      {
//       outfile << "# Column 2: East-west displacement (m)" << endl;
//       outfile << "# Column 3: North-south displacement (m)" << endl;
//       outfile << "# Column 4: Up-down displacement (m)" << endl;
	 fprintf(fd, "# Column 2: East-west displacement (m)\n");
	 fprintf(fd, "# Column 3: North-south displacement (m)\n");
	 fprintf(fd, "# Column 4: Up-down displacement (m)\n");
      }
      else if( m_xycomponent && m_velocities )
      {
//       outfile << "# Column 2: X velocity (m/s)" << endl;
//       outfile << "# Column 3: Y velocity (m/s)" << endl;
//       outfile << "# Column 4: Z velocity (m/s)" << endl;
	 fprintf(fd, "# Column 2: X velocity (m/s)\n");
	 fprintf(fd, "# Column 3: Y velocity (m/s)\n");
	 fprintf(fd, "# Column 4: Z velocity (m/s)\n");
      }
      else
      {
//       outfile << "# Column 2: X displacement (m)" << endl;
//       outfile << "# Column 3: Y displacement (m)" << endl;
//       outfile << "# Column 4: Z displacement (m)" << endl;
	 fprintf(fd, "# Column 2: X displacement (m)\n");
	 fprintf(fd, "# Column 3: Y displacement (m)\n");
	 fprintf(fd, "# Column 4: Z displacement (m)\n");
      }
   }
   else if( m_div )
   {
      if( m_velocities )
	 fprintf(fd, "# Column 2: divergence of velocity (/s)\n");
      else
	 fprintf(fd, "# Column 2: divergence of displacement ()\n");
   }
   else if( m_curl )
   {
      if( m_velocities )
      {
	 fprintf(fd, "# Column 2: curl of velocity, component 1 (/s)\n");
	 fprintf(fd, "# Column 3: curl of velocity, component 2 (/s)\n");
	 fprintf(fd, "# Column 4: curl of velocity, component 3 (/s)\n");
      }
      else
      {
	 fprintf(fd, "# Column 2: curl of displacement, component 1 ()\n");
	 fprintf(fd, "# Column 3: curl of displacement, component 2 ()\n");
	 fprintf(fd, "# Column 4: curl of displacement, component 3 ()\n");
      }
   }
   else if( m_strains )
   {
      if( m_velocities )
      {
	 fprintf(fd, "# Column 2: xx velocity strain component (/s)\n");
	 fprintf(fd, "# Column 3: yy velocity strain component (/s)\n");
	 fprintf(fd, "# Column 4: zz velocity strain component (/s)\n");
	 fprintf(fd, "# Column 5: xy velocity strain component (/s)\n");
	 fprintf(fd, "# Column 6: xz velocity strain component (/s)\n");
	 fprintf(fd, "# Column 7: yz velocity strain component (/s)\n");
      }
      else
      {
	 fprintf(fd, "# Column 2: xx strain component ()\n");
	 fprintf(fd, "# Column 3: yy strain component ()\n");
	 fprintf(fd, "# Column 4: zz strain component ()\n");
	 fprintf(fd, "# Column 5: xy strain component ()\n");
	 fprintf(fd, "# Column 6: xz strain component ()\n");
	 fprintf(fd, "# Column 7: yz strain component ()\n");
      }
   }
// write the data
   if( m_strains )
      for( int i = 0 ; i < mRecordedUX.size() ; i++ )
      {
//      outfile << mWPP->mTstart - mWPP->m_t0Shift + i*mWPP->mDt << " " << mRecordedUX[i] << " " <<
//        mRecordedUY[i] << " " << mRecordedUZ[i] << endl;
	 fprintf(fd, "%e %e %e %e %e %e %e\n", mWPP->mTstart - mWPP->m_t0Shift + i*mWPP->mDt, mRecordedUX[i], 
		 mRecordedUY[i], mRecordedUZ[i], mRecordedUXY[i], mRecordedUXZ[i], mRecordedUYZ[i] );
      }
   else if( m_div )
      for( int i = 0 ; i < mRecordedUX.size() ; i++ )
      {
//      outfile << mWPP->mTstart - mWPP->m_t0Shift + i*mWPP->mDt << " " << mRecordedUX[i] << " " <<
//        mRecordedUY[i] << " " << mRecordedUZ[i] << endl;
	 fprintf(fd, "%e %e \n", mWPP->mTstart - mWPP->m_t0Shift + i*mWPP->mDt, mRecordedUX[i] ); 
      }
   else
      for( int i = 0 ; i < mRecordedUX.size() ; i++ )
      {
//      outfile << mWPP->mTstart - mWPP->m_t0Shift + i*mWPP->mDt << " " << mRecordedUX[i] << " " <<
//        mRecordedUY[i] << " " << mRecordedUZ[i] << endl;
	 fprintf(fd, "%e %e %e %e\n", mWPP->mTstart - mWPP->m_t0Shift + i*mWPP->mDt, mRecordedUX[i], 
		 mRecordedUY[i], mRecordedUZ[i]);
      }
   
//   outfile.close();
   fclose(fd);
   
}
