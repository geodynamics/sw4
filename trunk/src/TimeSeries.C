// -*-c++-*-
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "TimeSeries.h"
#include "mpi.h"
#include "sacsubc.h"
#include "csstime.h"

#include "Require.h"
#include "Filter.h"

#include "EW.h"

using namespace std;

void parsedate( char* datestr, int& year, int& month, int& day, int& hour, int& minute,
		int& second, int& msecond, int& fail );

TimeSeries::TimeSeries( EW* a_ew, std::string fileName, receiverMode mode, bool sacFormat, bool usgsFormat, 
			double x, double y, double depth, bool topoDepth, int writeEvery, bool xyzcomponent ):
  m_mode(mode),
  m_nComp(0),
  m_myPoint(false),
  m_fileName(fileName),
  m_path(a_ew->getOutputPath()),
  mX(x),
  mY(y),
  mZ(depth),
  mGPX(0.0),
  mGPY(0.0),
  mGPZ(0.0),
  m_zRelativeToTopography(topoDepth),
  mWriteEvery(writeEvery),
  m_usgsFormat(usgsFormat),
  m_sacFormat(sacFormat),
  m_xyzcomponent(xyzcomponent),
  m_i0(-999),
  m_j0(-999),
  m_k0(-999),
  m_grid0(-999),
  m_t0(0.0),
  m_dt(1.0),
  mAllocatedSize(-1),
  mLastTimeStep(-1),
  mRecordedSol(NULL),
  mRecordedFloats(NULL),
  mIgnore(true), // are we still using this flag???
  mEventYear(2012),
  mEventMonth(2),
  mEventDay(8),
  mEventHour(10),
  mEventMinute(28),
  mEventSecond(0.0),
  m_rec_lat(38.0),
  m_rec_lon(-122.5),
  m_epi_lat(38.0),
  m_epi_lon(-122.5),
  m_epi_depth(0.0),
  m_epi_time_offset(0.0),
  m_x_azimuth(0.0),
  mBinaryMode(true),
  m_utc_set(false),
  m_utc_offset_computed(false),
  m_use_win(false),
  m_use_x(true),
  m_use_y(true),
  m_use_z(true)
{
// preliminary determination of nearest grid point ( before topodepth correction to mZ)
   a_ew->computeNearestGridPoint(m_i0, m_j0, m_k0, m_grid0, mX, mY, mZ);

// does this processor write this station?
   m_myPoint = a_ew->interior_point_in_proc(m_i0, m_j0, m_grid0);

// tmp
//   printf("TimeSeries constructor, rank=%i, myPoint=%i\n", a_ew->getRank(), m_myPoint);

// The following is a safety check to make sure only one processor writes each time series.
// We could remove this check if we were certain that interior_point_in_proc() never lies
  int iwrite = m_myPoint ? 1 : 0;
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
  " for receiver station " << m_fileName );

  if (!m_myPoint) return;

// from here on this processor writes this sac station and knows about its topography

  double zMinTilde, q, r, s;
  if (a_ew->topographyExists())
  {
// NOT FULLY IMPLEMENTED
    int gCurv = a_ew->mNumberOfGrids - 1;
    double h = a_ew->mGridSize[gCurv];
    q = mX/h + 1.0;
    r = mY/h + 1.0;
// evaluate elevation of topography on the grid
    // if (!a_ew->interpolate_topography(q, r, zMinTilde, true))
    // {
    //   cerr << "Unable to evaluate topography for SAC station" << m_fileName << " mX= " << mX << " mY= " << mY << endl;
    //   cerr << "Setting topography to ZERO" << endl;
    //   zMinTilde = 0;
    // }
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
     m_fileName.c_str(),  mX,  mY, mZ, zMinTilde);
// don't write this station
     m_myPoint=false;
     return;
   }
     
// now we can find the closest grid point  
   a_ew->computeNearestGridPoint(m_i0, m_j0, m_k0, m_grid0, mX, mY, mZ);

   if( m_grid0 == a_ew->mNumberOfGrids-1 && a_ew->topographyExists() )
   {
// Curvilinear
     bool canBeInverted = a_ew->invert_curvilinear_grid_mapping( mX, mY, mZ, q, r, s );
     if (a_ew->invert_curvilinear_grid_mapping( mX, mY, mZ, q, r, s )) // the inversion was successful
     {
       m_k0 = (int)floor(s);
       if (s-(m_k0+0.5) > 0.) m_k0++;
       m_k0 = max(a_ew->m_kStartInt[m_grid0], m_k0);
       int Nz = a_ew->m_kEndInt[m_grid0];
       m_k0 = min(Nz, m_k0);
     }
     else
     {
       cerr << "Can't invert curvilinear grid mapping for recevier station" << m_fileName << " mX= " << mX << " mY= " 
	    << mY << " mZ= " << mZ << endl;
       cerr << "Placing the station on the surface (depth=0)." << endl;
       m_k0 = 1;
     }
   }
   
// actual location of station (nearest grid point)
   double xG, yG, zG;
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
   
   if (a_ew->getVerbosity()>=2 && fabs(mX-xG)+fabs(mY-yG)+fabs(mZ-zG) > 0.001*a_ew->mGridSize[m_grid0] )
   {
     cout << "recevier info for station " << m_fileName << ":" << 
       " initial location (x,y,z) = " << mX << " " << mY << " " << mZ << 
       " moved to nearest grid point (x,y,z) = " << xG << " " << yG << " " << zG << 
       " h= " << a_ew->mGridSize[m_grid0] <<
       " with indices (i,j,k)= " << m_i0 << " " << m_j0 << " " << m_k0 << " in grid " << m_grid0 << endl;
   }

// remember corrected location
   mGPX = xG;
   mGPY = yG;
   mGPZ = zG;

// remember file prefix
   if (a_ew->getPath() != ".")
   {
     m_filePrefix = a_ew->getPath();
   }
   else
   {
     m_filePrefix = "";
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

// allocate handles to solution array pointers
  mRecordedSol = new double* [m_nComp];
  
  for (int q=0; q<m_nComp; q++)
    mRecordedSol[q] = NULL;

// keep a copy for saving on a sac file
  if (m_sacFormat)
  {
    mRecordedFloats = new float* [m_nComp];
  
    for (int q=0; q<m_nComp; q++)
      mRecordedFloats[q] = NULL;
  }
  
  
// do some misc pre computations
  m_x_azimuth = a_ew->getGridAzimuth(); // degrees
  
  a_ew->computeGeographicCoord(mX, mY, m_rec_lon, m_rec_lat);
  a_ew->computeGeographicCoord(mGPX, mGPY, m_rec_gp_lon, m_rec_gp_lat);

  m_calpha = cos(M_PI*m_x_azimuth/180.0);
  m_salpha = sin(M_PI*m_x_azimuth/180.0);

  double cphi   = cos(M_PI*m_rec_lat/180.0);
  double sphi   = sin(M_PI*m_rec_lat/180.0);

  double metersperdegree = a_ew->getMetersPerDegree();

  m_thxnrm = m_salpha + (mX*m_salpha+mY*m_calpha)/cphi/metersperdegree * (M_PI/180.0) * sphi * m_calpha;
  m_thynrm = m_calpha - (mX*m_salpha+mY*m_calpha)/cphi/metersperdegree * (M_PI/180.0) * sphi * m_salpha;
  double nrm = sqrt( m_thxnrm*m_thxnrm + m_thynrm*m_thynrm );
  m_thxnrm /= nrm;
  m_thynrm /= nrm;

   // m_dthi = 1.0/(2*a_ew->mDt);

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
void TimeSeries::allocateRecordingArrays( int numberOfTimeSteps, double startTime, double timeStep )
{
  if (!m_myPoint) return; // only one processor saves each time series
  //  cout << "Time series allocaterecarr " <<  numberOfTimeSteps << " " << startTime << " " << timeStep << endl;
  if (numberOfTimeSteps > 0)
  {
    mAllocatedSize = numberOfTimeSteps+1;
    mLastTimeStep = -1;

    for (int q=0; q<m_nComp; q++)
    {
      if (mRecordedSol[q]) delete [] mRecordedSol[q];
      mRecordedSol[q] = new double[mAllocatedSize];
    }

    if (m_sacFormat)
    {
      for (int q=0; q<m_nComp; q++)
      {
	if (mRecordedFloats[q]) delete [] mRecordedFloats[q];
	mRecordedFloats[q] = new float[mAllocatedSize];
      }
    }
    
  }
  m_t0 = startTime;
  m_dt = timeStep;
}

//--------------------------------------------------------------
void TimeSeries::recordData(vector<double> & u) 
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
  
  stringstream ux, uy, uz, uxy, uxz, uyz;
  
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

     string xfield, yfield, zfield, xyfield, xzfield, yzfield;
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
			mRecordedFloats[0], (float) m_t0, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uy.str().c_str()), 
			mRecordedFloats[1], (float) m_t0, (float) m_dt,
			const_cast<char*>(yfield.c_str()), 90.0, azimy); 
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uz.str().c_str()), 
			mRecordedFloats[2], (float) m_t0, (float) m_dt,
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
			geographic[0], (float) m_t0, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uy.str().c_str()), 
			geographic[1], (float) m_t0, (float) m_dt,
			const_cast<char*>(yfield.c_str()), 90.0, azimy); 
	   write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uz.str().c_str()), 
			geographic[2], (float) m_t0, (float) m_dt,
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
			mRecordedFloats[0], (float) m_t0, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
     }
     else if (m_mode == Strains) // 6 components
     {
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(ux.str().c_str()), 
			mRecordedFloats[0], (float) m_t0, (float) m_dt,
			const_cast<char*>(xfield.c_str()), 90.0, azimx); 
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uy.str().c_str()), 
			mRecordedFloats[1], (float) m_t0, (float) m_dt,
			const_cast<char*>(yfield.c_str()), 90.0, azimy); 
       write_sac_format(mLastTimeStep+1, 
			const_cast<char*>(uz.str().c_str()), 
			mRecordedFloats[2], (float) m_t0, (float) m_dt,
			const_cast<char*>(zfield.c_str()), updownang, 0.0); 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uxy.str().c_str()), 
			mRecordedFloats[3], (float) m_t0, (float) m_dt,
			const_cast<char*>(xyfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uxz.str().c_str()), 
			mRecordedFloats[4], (float) m_t0, (float) m_dt,
			const_cast<char*>(xzfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
       write_sac_format(mLastTimeStep+1,
			const_cast<char*>(uyz.str().c_str()), 
			mRecordedFloats[5], (float) m_t0, (float) m_dt,
			const_cast<char*>(yzfield.c_str()), 90.0, azimx); // not sure what the updownang or azimuth should be here 
     }
  } // end if m_sacFormat
  
  if( m_usgsFormat )
  {
    filePrefix << "txt";
    cout << "Writing ASCII USGS file, "
    	 << "of size " << mLastTimeStep+1 << ": "
	 << filePrefix.str() << endl;

    write_usgs_format( filePrefix.str() );
  }

}

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
  if( m_utc_set )
  {
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
  }
  else
  {
     setnhv( nm[10], mEventYear, nerr);
     setnhv( nm[11], mEventDay, nerr);
     setnhv( nm[12], mEventHour, nerr);
     setnhv( nm[13], mEventMinute, nerr);
     setnhv( nm[14], static_cast<int>(mEventSecond), nerr);
     setnhv( nm[15], 0, nerr);
  }

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
  setkhv( nm[25], const_cast<char*>(m_fileName.c_str()), nerr);


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
   double lat, lon;
   double x, y, z;

   if( fd == NULL )
      cout << "ERROR: opening USGS file " << a_fileName << " for writing" <<  endl;
   else
   {

// frequency resolution
//    double freq_limit=-999;
//    if (a_ew->m_prefilter_sources)
//      freq_limit = a_ew->m_fc;
//    else if (a_ew->m_limit_frequency)
//      freq_limit = a_ew->m_frequency_limit;

// write the header
   fprintf(fd, "# Author: SBP4W\n");
   fprintf(fd, "# Scenario: %s\n", "test"/*a_ew->m_scenario.c_str()*/);
   if( m_utc_set )
      fprintf(fd, "# Date: UTC  %02i/%02i/%i:%i:%i:%i.%i\n", m_utc[1], m_utc[2], m_utc[0], m_utc[3],
	                                      m_utc[4], m_utc[5], m_utc[6] );
   else
      fprintf(fd, "# Date: %i-%s-%i\n", mEventDay, mname[mEventMonth].c_str(), mEventYear);

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

// write the data

   if( m_xyzcomponent || (!m_xyzcomponent && m_mode == Div) )
   {
      for( int i = 0 ; i <= mLastTimeStep ; i++ )
      {
	 fprintf(fd, "%e", m_t0 + i*m_dt);
	 for (int q=0; q<m_nComp; q++)
	    fprintf(fd, " %20.12g", mRecordedSol[q][i]);
	 fprintf(fd, "\n");
      }
   }
   else if( m_mode == Displacement || m_mode == Velocity )
   {
      for( int i = 0 ; i <= mLastTimeStep ; i++ )
      {
	 fprintf(fd, "%e", m_t0 + i*m_dt);
	 double uns = m_thynrm*mRecordedSol[0][i]-m_thxnrm*mRecordedSol[1][i];
	 double uew = m_salpha*mRecordedSol[0][i]+m_calpha*mRecordedSol[1][i];
	 fprintf(fd, " %20.12g", uew );
	 fprintf(fd, " %20.12g", uns );
	 fprintf(fd, " %20.12g", -mRecordedSol[2][i] );
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
      printf(" in geographic coordinates\n" );
   }
   fclose(fd);
   }
}

//-----------------------------------------------------------------------
void TimeSeries::readFile( EW *ew, double timeshift )
{
//building the file name...
// 
   stringstream filePrefix;
   if( ew->getObservationPath() != "./" )
      filePrefix << ew->getObservationPath();
   filePrefix << m_fileName << ".txt" ;

   if( m_myPoint && m_usgsFormat )
   {
      bool debug = false;
      FILE *fd=fopen(filePrefix.str().c_str(),"r");
      if( fd == NULL )
	 cout << "ERROR: observed data file " << filePrefix.str() << " not found " << endl;
      else
      {
	 int bufsize=1024;
	 char* buf = new char[bufsize];

      // Read header
	 for( int line=0 ; line < 13 ; line++ )
	 {
	    fgets(buf,bufsize,fd);
            if( line == 2 )
	    {
	       // set UTC time, if defined in file
	       string timestr(buf);
	       size_t utcind = timestr.find("UTC");
               if( utcind != string::npos  )
	       {
                  int fail;
                  char *utcstr=new char[1024];
		  // Skip characters 'UTC'
		  utcind += 3;

		  timestr.copy(utcstr, timestr.size()-utcind , utcind );
		  ew->parsedate( utcstr, m_utc[0], m_utc[1], m_utc[2], m_utc[3],
			     m_utc[4], m_utc[5], m_utc[6], fail );
                  if( fail == 0 )
		     m_utc_set = true;
                  else
		     cout << "ERROR reading observation " << m_fileName << " , UTC parse failure no. " << fail << endl;
                  delete[] utcstr;
                  if( debug )
		  {
		     cout << "found observation utc time " ;
                     for( int u=0 ; u < 7 ; u++ )
			cout << m_utc[u] << " ";
		     cout << endl;
		  }
	       }
	    }
	 }
         string bufstr(buf);
         bool foundd = (bufstr.find("displacement") != string::npos);
         bool foundv = (bufstr.find("velocity") != string::npos);
         if( foundd || foundv )
	 {
	    // The file contains velocities or displacements. 
	    // The last line read contains the z-component. Check whether it is a ENU or XYZ file.
            bool cartesian = (bufstr.find("Z") != string::npos);
	    
            m_xyzcomponent = cartesian;

            if( debug )
	    {
	       cout << "Found observed ";
	       if( foundd )
		  cout << "displacement ";
	       else
		  cout << "velocity ";
	       cout << "file with ";
	       if( cartesian )
		  cout << "Cartesian ";
	       else
		  cout << "geographic ";
	       cout << "components " << endl;
	    }
	    double tstart, dt, td, ux, uy, uz;
	    int nlines = 0;
	    if( fscanf(fd,"%le %le %le %le",&tstart,&ux,&uy,&uz) != EOF )
	       nlines++;
	    if( fscanf(fd,"%le %le %le %le",&dt,&ux,&uy,&uz) != EOF )
	       nlines++;
	    dt = dt-tstart;
	    while( fscanf(fd,"%le %le %le %le",&td,&ux,&uy,&uz) != EOF )
	       nlines++;

	    fclose(fd);
	    // Use offset in time column, and add a possible timeshift.
	    if( nlines > 2 )
		  allocateRecordingArrays( nlines, tstart+timeshift, dt );
	    else
	    {
	       cout << "ERROR: observed data is too short" << endl;
	       cout << "    File " << filePrefix.str() << " not read." << endl;
	    }
	    fd=fopen(filePrefix.str().c_str(),"r");	 
	    if( fd == NULL )
	       cout << "ERROR: observed data file " << filePrefix.str() << " could not be reopened" << endl;

	 // Read past header
	    for( int line=0 ; line < 13 ; line++ )
	       fgets(buf,bufsize,fd);
	    // Mapping to invert (e,n) to (x,y) components, Only needed in the non-cartesian case.
            double deti = 1.0/(m_thynrm*m_calpha+m_thxnrm*m_salpha);
            double a11 = m_calpha*deti;
	    double a12 = m_thxnrm*deti;
            double a21 =-m_salpha*deti;
	    double a22 = m_thynrm*deti;
	 // Read the data on file	 
            if( debug )
	       cout << "Found " << nlines << " lines of observation data " << endl;
            for( int line=0 ; line < nlines ; line++ )
	    {
	       int nr=fscanf(fd, "%lf %lf %lf %lf \n", &tstart,&ux,&uy,&uz);
       //               cout << "nr = " << nr << " tstart " << tstart << " ux " << ux << " uy " << uy << " uz " << uz << endl;
               if( nr != 4 )
	       {
		  cout << "ERROR: could not read observed data file " << endl;
	       }
               if( cartesian )
	       {
		  mRecordedSol[0][line]=ux;
		  mRecordedSol[1][line]=uy;
		  mRecordedSol[2][line]=uz;
	       }
	       else
	       {
		  // Geographic coordinates, read East, North, and Up velocities/displacements
                  double uns = uy;
                  double uew = ux;
		  mRecordedSol[2][line] = -uz;
		  // Transform to Cartesian
                  mRecordedSol[0][line] = a11*uns + a12*uew;
		  mRecordedSol[1][line] = a21*uns + a22*uew;
	       }

	    }
	    fclose(fd);
            mLastTimeStep = nlines-1;
	 }
	 else
	 {
	    cout << "ERROR: observed data must contain displacements" << endl;
	    cout << "    File " << filePrefix.str() << " not read " << endl;
	 }
	 delete[] buf;
      }
   }
   else if( !m_usgsFormat )
   {
      cout << "ERROR: observed data must be an ASCII USGS file" << endl;
      cout << "    File " << filePrefix.str() << " not read " << endl;
   }
}

//-----------------------------------------------------------------------
void TimeSeries::interpolate( TimeSeries& intpfrom )
{
// Interpolate data to the grid in this object
   int order = 4;

   double dtfr  = intpfrom.m_dt;
   double t0fr  = intpfrom.m_t0;
   int nfrsteps = intpfrom.mLastTimeStep+1;
   for( int i= 0 ; i <= mLastTimeStep ; i++ )
   {
      double t  = m_t0 + i*m_dt;
      double ir = (t-t0fr)/dtfr;
      int ie   = static_cast<int>(ir);
      int mmin = ie-order/2+1;
      int mmax = ie+order/2;
      if( m_usgsFormat )
	 mRecordedSol[0][i] = mRecordedSol[1][i] = mRecordedSol[2][i] = 0;
      else
	 mRecordedFloats[0][i] = mRecordedFloats[1][i] = mRecordedFloats[2][i] = 0;
      if( mmax-mmin+1 > nfrsteps )
      {
         cout << "Error in TimeSeries::interpolate : Can not interpolate, " <<
	    "because the grid is too coarse " << endl;
	 return;
      }

// If too far past the end of intpfrom, extrapolate constant value
      if( ie > nfrsteps + order/2 )
      {
	 if( m_usgsFormat )
	 {
	    mRecordedSol[0][i] = mRecordedSol[0][i-1];
	    mRecordedSol[1][i] = mRecordedSol[1][i-1];
	    mRecordedSol[2][i] = mRecordedSol[2][i-1];
	 }
	 else
	 {
	    mRecordedFloats[0][i] = mRecordedFloats[0][i-1];
            mRecordedFloats[1][i] = mRecordedFloats[1][i-1];
            mRecordedFloats[2][i] = mRecordedFloats[2][i-1];
	 }
      }
      else if( ie < 1 )
      {
// Before start of intpfrom, use first value
	 if( m_usgsFormat )
	 {
	    mRecordedSol[0][i] = intpfrom.mRecordedSol[0][0];
	    mRecordedSol[1][i] = intpfrom.mRecordedSol[1][0];
	    mRecordedSol[2][i] = intpfrom.mRecordedSol[2][0];
	 }
	 else
	 {
	    mRecordedFloats[0][i] = intpfrom.mRecordedFloats[0][0];
            mRecordedFloats[1][i] = intpfrom.mRecordedFloats[1][0];
            mRecordedFloats[2][i] = intpfrom.mRecordedFloats[2][0];
	 }

      }
      else
      {
         int off;
	 if( mmin < 1 )
	 {
	    off = 1-mmin;
	    mmin = mmin + off;
	    mmax = mmax + off;
	 }
	 if( mmax > nfrsteps )
	 {
	    off = mmax - nfrsteps;
	    mmin = mmin - off;
	    mmax = mmax - off;
	 }
	 for( int m = mmin ; m <= mmax ; m++ )
	 {
	    double cof=1;
	    for( int k = mmin ; k <= mmax ; k++ )
	       if( k != m )
		  cof *= (ir-k)/(m-k);

	    if( m_usgsFormat && intpfrom.m_usgsFormat )
	    {
	       mRecordedSol[0][i] += cof*intpfrom.mRecordedSol[0][m];
	       mRecordedSol[1][i] += cof*intpfrom.mRecordedSol[1][m];
	       mRecordedSol[2][i] += cof*intpfrom.mRecordedSol[2][m];
	    }
	    else if( m_usgsFormat && !intpfrom.m_usgsFormat )
	    {
	       mRecordedSol[0][i] += cof*intpfrom.mRecordedFloats[0][m];
	       mRecordedSol[1][i] += cof*intpfrom.mRecordedFloats[1][m];
	       mRecordedSol[2][i] += cof*intpfrom.mRecordedFloats[2][m];
	    }
	    else if( !m_usgsFormat && intpfrom.m_usgsFormat )
	    {
	       mRecordedFloats[0][i] += cof*intpfrom.mRecordedSol[0][m];
	       mRecordedFloats[1][i] += cof*intpfrom.mRecordedSol[1][m];
	       mRecordedFloats[2][i] += cof*intpfrom.mRecordedSol[2][m];
	    }
	    else if( !m_usgsFormat && !intpfrom.m_usgsFormat )
	    {
	       mRecordedFloats[0][i] += cof*intpfrom.mRecordedFloats[0][m];
	       mRecordedFloats[1][i] += cof*intpfrom.mRecordedFloats[1][m];
	       mRecordedFloats[2][i] += cof*intpfrom.mRecordedFloats[2][m];
	    }
	 }
      }
   }
}

//-----------------------------------------------------------------------
double TimeSeries::misfit( TimeSeries& observed, TimeSeries* diff )
{
   // Computes  diff := this - observed
   //       where 'diff' has the same grid points as 'this'. 'observed' is
   //       interpolated to this grid, and is set to zero outside its interval
   //       of definition.
   double misfit = 0;
   if( m_myPoint )
   {
// Interpolate data to this object
      int order = 4;
      double mf[3];
      double dtfr  = observed.m_dt;
      double t0fr  = observed.m_t0;
      int nfrsteps = observed.mLastTimeStep+1;

      bool compute_difference = (diff!=NULL);
      double** misfitsource;
      if( compute_difference )
      {
	 if( diff->mLastTimeStep < mLastTimeStep )
	    diff->allocateRecordingArrays(mLastTimeStep,m_t0,m_dt);
	 misfitsource = diff->getRecordingArray();
      }

      // Weight to ramp down the end of misfit.
      double wghv;
      int p =20 ; // Number of points in ramp;
      int istart = 1; // Starting index for downward ramp.
      if( mLastTimeStep-p+1 > 1 )
	 istart = mLastTimeStep-p+1;

      //      cout << "in misfit, " ;
      //      if( m_use_win )
      //	 cout << "using window= " << m_winL << " to " << m_winR << endl;
      //      else
      //	 cout << "not using window " << endl;

      //      if( m_use_x && m_use_y && m_use_z )
      //	 cout << "not excluding any component " << endl;
      //      if( !m_use_x )
      //	 cout << "excluding component x" << endl;
      //      if( !m_use_y )
      //	 cout << "excluding component y" << endl;
      //      if( !m_use_z )
      //	 cout << "excluding component z" << endl;

      for( int i= 0 ; i <= mLastTimeStep ; i++ )
      {
	 wghv = 1;
	 if( i >= istart )
	 {
	    double arg = (mLastTimeStep-i)/(p-1.0);
	    wghv = arg*arg*arg*arg*(35-84*arg+70*arg*arg-20*arg*arg*arg);
	 }


	 double t  = m_t0 + i*m_dt;
	 double ir = (t-t0fr)/dtfr;
	 int ie   = static_cast<int>(ir);
	 int mmin = ie-order/2+1;
	 int mmax = ie+order/2;
	 if( mmax-mmin+1 > nfrsteps )
	 {
	    cout << "Error in TimeSeries::misfit : Can not interpolate, " <<
	       "because the grid is too coarse " << endl;
            cout << "mmin = " << mmin << endl;
            cout << "mmax = " << mmax << endl;
	    cout << "nfrsteps = " << nfrsteps << endl;
	    return 0.0;
	 }

// Windowing and component selection
         double wghx, wghy, wghz;
	 wghx = wghy = wghz = wghv;
         if( m_use_win && (t < m_winL || t > m_winR) )
	    wghx = wghy = wghz = 0;
         if( !m_use_x )
	    wghx = 0;
         if( !m_use_y )
	    wghy = 0;
         if( !m_use_z )
	    wghz = 0;

	 mf[0] = mf[1] = mf[2] = 0;

	 // If too far past the end of observed, set to zero.
	 if( ie > nfrsteps + order/2 )
	 {
	    mf[0] = mf[1] = mf[2] = 0;
	 }
         else if( ie < 1 )
	 {
	    // Before the starting point of the observations.
	    mf[0] = mf[1] = mf[2] = 0;
	 }
	 else
	 {
	    int off;
	    if( mmin < 1 )
	    {
	       off = 1-mmin;
	       mmin = mmin + off;
	       mmax = mmax + off;
	    }
	    if( mmax > nfrsteps )
	    {
	       off = mmax - nfrsteps;
	       mmin = mmin - off;
	       mmax = mmax - off;
	    }
	    for( int m = mmin ; m <= mmax ; m++ )
	    {
	       double cof=1;
	       for( int k = mmin ; k <= mmax ; k++ )
		  if( k != m )
		     cof *= (ir-k)/(m-k);

	       if( observed.m_usgsFormat )
	       {
		  mf[0] += cof*observed.mRecordedSol[0][m];
		  mf[1] += cof*observed.mRecordedSol[1][m];
		  mf[2] += cof*observed.mRecordedSol[2][m];
	       }
	       else
	       {
		  mf[0] += cof*observed.mRecordedFloats[0][m];
		  mf[1] += cof*observed.mRecordedFloats[1][m];
		  mf[2] += cof*observed.mRecordedFloats[2][m];
	       }
	    }

	    if( m_usgsFormat )
	    {
	       misfit += ( (mf[0]-mRecordedSol[0][i])*(mf[0]-mRecordedSol[0][i])*wghx + 
			   (mf[1]-mRecordedSol[1][i])*(mf[1]-mRecordedSol[1][i])*wghy + 
			   (mf[2]-mRecordedSol[2][i])*(mf[2]-mRecordedSol[2][i])*wghz    );
	       if( compute_difference )
	       {
		  misfitsource[0][i] = wghx*(mRecordedSol[0][i]-mf[0]);
		  misfitsource[1][i] = wghy*(mRecordedSol[1][i]-mf[1]);
		  misfitsource[2][i] = wghz*(mRecordedSol[2][i]-mf[2]);
	       }
	    }
	    else
	    {
	       misfit += ( (mf[0]-mRecordedFloats[0][i])*(mf[0]-mRecordedFloats[0][i])*wghx + 
			   (mf[1]-mRecordedFloats[1][i])*(mf[1]-mRecordedFloats[1][i])*wghy + 
			   (mf[2]-mRecordedFloats[2][i])*(mf[2]-mRecordedFloats[2][i])*wghz );
	       if( compute_difference )
	       {
		  misfitsource[0][i] = (mRecordedFloats[0][i]-mf[0])*wghx;
		  misfitsource[1][i] = (mRecordedFloats[1][i]-mf[1])*wghy;
		  misfitsource[2][i] = (mRecordedFloats[2][i]-mf[2])*wghz;
	       }
	    }
	 }
      }
   }
   else
      misfit = 0;
   return 0.5*misfit;
}

//-----------------------------------------------------------------------
TimeSeries* TimeSeries::copy( EW* a_ew, string filename, bool addname )
{
   if( addname )
      filename = m_fileName + filename;

   TimeSeries* retval = new TimeSeries( a_ew, filename, m_mode, m_sacFormat, m_usgsFormat,
					mX, mY, mZ, m_zRelativeToTopography, mWriteEvery, m_xyzcomponent );
   retval->m_t0 = m_t0;
   retval->m_dt = m_dt;
   retval->mAllocatedSize = mAllocatedSize;
   retval->mLastTimeStep  = mLastTimeStep;
   //   if( m_myPoint )
   //      cout << "In copy, xyz = " << m_xyzcomponent << endl;

// Component rotation:
   retval->m_calpha = m_calpha;
   retval->m_salpha = m_salpha;
   retval->m_thxnrm = m_thxnrm;
   retval->m_thynrm = m_thynrm;
//   retval->m_xyzcomponent = m_xyzcomponent;

// UTC time reference point:
   retval->m_utc_set = m_utc_set;
   for( int c=0; c < 7; c++ )
      retval->m_utc[c] = m_utc[c];
	 
// windows and exclusions
   retval->m_use_win = m_use_win;
   retval->m_winL = m_winL;
   retval->m_winR = m_winR;
   retval->m_use_x = m_use_x;
   retval->m_use_y = m_use_y;
   retval->m_use_z = m_use_z;
   if( m_myPoint )
   {
      if( m_sacFormat )
      {
	 // Overwrite pointers, don't want to copy them.
	 retval->mRecordedFloats = new float*[m_nComp];
	 if( mAllocatedSize > 0 )
	 {
	    for( int q=0 ; q < m_nComp ; q++ )
	       retval->mRecordedFloats[q] = new float[mAllocatedSize];
	    for( int q=0 ; q < m_nComp ; q++ )
	       for( int i=0 ; i < mAllocatedSize ; i++ )
		  retval->mRecordedFloats[q][i] = mRecordedFloats[q][i];
	 }
	 else
	 {
	    for( int q=0 ; q < m_nComp ; q++ )
	       retval->mRecordedFloats[q] = NULL;
	 }
      }
      else
      {
	 retval->mRecordedSol = new double*[m_nComp];
	 if( mAllocatedSize > 0 )
	 {
	    for( int q=0 ; q < m_nComp ; q++ )
	       retval->mRecordedSol[q] = new double[mAllocatedSize];
	    for( int q=0 ; q < m_nComp ; q++ )
	       for( int i=0 ; i < mAllocatedSize ; i++ )
		  retval->mRecordedSol[q][i] = mRecordedSol[q][i];
	 }
	 else
	 {
	    for( int q=0 ; q < m_nComp ; q++ )
	       retval->mRecordedSol[q] = NULL;
	 }
      }
   }
   return retval;
}

//-----------------------------------------------------------------------
double TimeSeries::arrival_time( double lod )
{
  // Assume three components
   if( m_nComp != 3 )
   {
      cout << "ERROR: TimeSeries::arrival_time: Number of components must be three";
      cout << " not " << m_nComp << endl;
      return -1;
   }

   double* maxes = new double[m_nComp];
   for( int c=0 ; c < m_nComp ; c++ )
      maxes[c] = 0;
   
   int n;
   if( m_usgsFormat )
   {
      for( int i=0 ; i <= mLastTimeStep ; i++ )
	 for( int c=0 ; c < m_nComp ; c++ )
	    if( fabs(mRecordedSol[c][i]) > maxes[c] )
	       maxes[c] = fabs(mRecordedSol[c][i]);
      n=0;
  //            cout << "max = " << maxes[0] << " " << maxes[1] << " " << maxes[2] << endl;
      while( fabs(mRecordedSol[0][n])<maxes[0]*lod &&
	     fabs(mRecordedSol[1][n])<maxes[1]*lod &&
	     fabs(mRecordedSol[2][n])<maxes[2]*lod &&
	     n < mLastTimeStep )
	 n++;
   }
   else
   {
      for( int i=0 ; i <= mLastTimeStep ; i++ )
	 for( int c=0 ; c < m_nComp ; c++ )
	    if( fabs(mRecordedFloats[c][i]) > maxes[c] )
	       maxes[c] = fabs(mRecordedFloats[c][i]);
      n=0;
      while( fabs(mRecordedFloats[0][n])<maxes[0]*lod &&
	     fabs(mRecordedFloats[1][n])<maxes[1]*lod &&
	     fabs(mRecordedFloats[2][n])<maxes[2]*lod &&
	     n < mLastTimeStep )
	 n++;
   }

   delete[] maxes;
   return m_t0 + n*m_dt;
}

//-----------------------------------------------------------------------
void TimeSeries::use_as_forcing( int n, std::vector<Sarray>& f,
				 std::vector<double> & h, double dt )
{
   // Use at grid point, n, in the grid of this object.
   if( m_myPoint )
   {
      double normwgh[4]={17.0/48.0, 59.0/48.0, 43.0/48.0, 49.0/48.0 };
      double ih3 = 1.0/(h[m_grid0]*h[m_grid0]*h[m_grid0]);
      //      double ih3 = 1.0;
      double iwgh = 1.0;
      if( 1 <= m_k0 && m_k0 <= 4  )
	 iwgh = 1.0/normwgh[m_k0-1];
      ih3 *= iwgh;
      // Compensate for  dt^4/12 factor in forward corrector step.
      ih3 *= 12/(dt*dt);
      //      n = static_cast<int>(round( (t-m_t0)/m_dt ));
      if( n >= 0 && n <= mLastTimeStep )
      {
	 f[m_grid0](1,m_i0,m_j0,m_k0) -= mRecordedSol[0][n]*ih3;
	 f[m_grid0](2,m_i0,m_j0,m_k0) -= mRecordedSol[1][n]*ih3;
	 f[m_grid0](3,m_i0,m_j0,m_k0) -= mRecordedSol[2][n]*ih3;
      }
   }
}

//-----------------------------------------------------------------------
double TimeSeries::product( TimeSeries& ts ) const
{
   // No weighting, use if one of the time series already has
   // been multiplied by wgh, such as returned by the mistfit function
   double prod = 0;
   if( mLastTimeStep == ts.mLastTimeStep )
   {
      for( int i= 0 ; i <= mLastTimeStep ; i++ )
      {
	 prod += ts.mRecordedSol[0][i]*mRecordedSol[0][i] +
	    ts.mRecordedSol[1][i]*mRecordedSol[1][i] + 
	    ts.mRecordedSol[2][i]*mRecordedSol[2][i];
      }
   }
   else
      cout << "TimeSeries::product: Error time series have incompatible sizes" << endl;
   return prod;
}

//-----------------------------------------------------------------------
double TimeSeries::product_wgh( TimeSeries& ts ) const
{
   // Product which uses weighting, for computing Hessian
   double prod = 0;
   if( mLastTimeStep == ts.mLastTimeStep )
   {
   // Weight to ramp down the end of misfit.
      double wghv;
      int p =20 ; // Number of points in ramp;
      int istart = 1;
      if( mLastTimeStep-p+1 > 1 )
	 istart = mLastTimeStep-p+1;


      for( int i= 0 ; i <= mLastTimeStep ; i++ )
      {
	 wghv = 1;
	 if( i >= istart )
	 {
	    double arg = (mLastTimeStep-i)/(p-1.0);
	    wghv = arg*arg*arg*arg*(35-84*arg+70*arg*arg-20*arg*arg*arg);
	 }

// Windowing and component selection
         double wghx, wghy, wghz;
	 wghx = wghy = wghz = wghv;
         double t = m_t0 + i*m_dt;
         if( m_use_win && (t < m_winL || t > m_winR) )
	    wghx = wghy = wghz = 0;
         if( !m_use_x )
	    wghx = 0;
         if( !m_use_y )
	    wghy = 0;
         if( !m_use_z )
	    wghz = 0;

	 prod += (ts.mRecordedSol[0][i]*mRecordedSol[0][i]*wghx +
  	          ts.mRecordedSol[1][i]*mRecordedSol[1][i]*wghy + 
	          ts.mRecordedSol[2][i]*mRecordedSol[2][i]*wghz );
      }
   }
   else
      cout << "TimeSeries::product_wgh: Error time series have incompatible sizes" << endl;
   return prod;
}

//-----------------------------------------------------------------------
void TimeSeries::set_station_utc( int utc[7] )
{
   if( m_utc_set )
      m_t0 = utc_distance(utc,m_utc) + m_t0;
   for( int k=0; k < 7 ; k++ )
      m_utc[k] = utc[k];
   m_utc_set = true;
}

//-----------------------------------------------------------------------
//void TimeSeries::set_station_utc( int utc[7] )
//{
// // Reference UTC for the station
//   for( int k=0; k < 7 ; k++ )
//      m_utc[k] = utc[k];
//   m_utc_set = true;
//}

//-----------------------------------------------------------------------
//void TimeSeries::offset_ref_utc( int utcref[7] )
//{
 // Evaluate distance to the reference UTC for the computation.
 // (The reference UTC corresponds to t=0 simulation time)
//   if( m_utc_set && !m_utc_offset_computed )
//   {
//      m_t0 = utc_distance(utcref,m_utc) + m_t0;
//      m_utc_offset_computed = true;
//      for( int c=0 ; c < 7 ; c++ )
//	 m_utc[c] = utcref[c];
//   }
//}   

//-----------------------------------------------------------------------
double TimeSeries::utc_distance( int utc1[7], int utc2[7] )
{
   // Compute time in seconds between two [y,M,d,h,m,s,ms] times
   // returns utc2-utc1 in seconds.
   // Discovered afterwards: UTC occasionally adds a leap second, so this 
   // function is not always completely accurate.

   int start[7], finish[7];
   int c=0;
   while( utc1[c] == utc2[c] && c <= 6 )
      c++;
   if( c == 7 )
      // Identical times
      return 0;
   else
   {
      bool onesmallest;
      if( utc1[c] < utc2[c] )
         for( int k= 0 ; k <7 ;k++)
	 {
	    start[k] = utc1[k];
            finish[k] = utc2[k];
            onesmallest = true;
	 }
      else
         for( int k= 0 ; k <7 ;k++)
	 {
	    start[k] = utc2[k];
            finish[k] = utc1[k];
            onesmallest = false;
	 }
      double d = 0;
      if( c <= 1 )
      {
	 // different month or year, count days
         d = 0;
         while( !(start[0]==finish[0] && start[1]==finish[1]) )
	 {
	    dayinc( start );
	    d++; 
	 }
      }
      d = d + finish[2]-start[2];
      // Convert days,min,secs, and msecs to seconds
      double sg = 1;
      if( !onesmallest )
         sg = -1;
      int ls = leap_second_correction(start,finish);
      return sg*(86400.0*d + (finish[3]-start[3])*3600.0 + (finish[4]-start[4])*60.0 +
		 (finish[5]-start[5]) +(finish[6]-start[6])*1e-3 + ls);
   }
}

//-----------------------------------------------------------------------
void TimeSeries::dayinc( int date[7] )
{
   date[2]++;
   if( date[2] > lastofmonth( date[0], date[1] ) )
   {
      date[2] = 1;
      date[1]++;
   }
   if( date[1] > 12 )
   {
      date[2] = date[1] = 1;
      date[0]++;
   }
}

//-----------------------------------------------------------------------
int TimeSeries::lastofmonth( int year, int month )
{
   int days;
   int leapyear=0;
   leapyear = ( year % 400 == 0 ) || ( (year % 4 == 0) && !(year % 100 == 0) );
   if( month == 2 )
      days = 28 + leapyear;
   else if( month==4 || month==6 || month==9 || month==11 )
      days = 30;
   else
      days = 31;
   return days;
}

//-----------------------------------------------------------------------
int TimeSeries::utccompare( int utc1[7], int utc2[7] )
{
   int c = 0;
   int retval;
   while( utc1[c] == utc2[c] && c <= 6 )
      c++;
   if( c == 7 )
      retval=0;
   else
   {
      if( utc1[c] < utc2[c] )
	 retval=-1;
      else
	 retval=1;
   }
   return retval;
}   

//-----------------------------------------------------------------------
int TimeSeries::leap_second_correction( int utc1[7], int utc2[7] )
// Count the number of leap seconds between two utc times.
{
   int* leap_sec_y, *leap_sec_m;
   int nls = 25;
   leap_sec_y = new int[nls];
   leap_sec_m = new int[nls];
   // Table of UTC leap seconds, added at the end of june (6) or december (12) 
   leap_sec_y[0] = 1972;
   leap_sec_m[0] = 6;
   leap_sec_y[1] = 1972;
   leap_sec_m[1] = 12;
   leap_sec_y[2] = 1973;
   leap_sec_m[2] = 12;
   leap_sec_y[3] = 1974;
   leap_sec_m[3] = 12;
   leap_sec_y[4] = 1975;
   leap_sec_m[4] = 12;
   leap_sec_y[5] = 1976;
   leap_sec_m[5] = 12;
   leap_sec_y[6] = 1977;
   leap_sec_m[6] = 12;
   leap_sec_y[7] = 1978;
   leap_sec_m[7] = 12;
   leap_sec_y[8] = 1979;
   leap_sec_m[8] = 12;
   leap_sec_y[9] = 1981;
   leap_sec_m[9] = 6;
   leap_sec_y[10] = 1982;
   leap_sec_m[10] = 6;
   leap_sec_y[11] = 1983;
   leap_sec_m[11] = 6;
   leap_sec_y[12] = 1985;
   leap_sec_m[12] = 6;
   leap_sec_y[13] = 1987;
   leap_sec_m[13] = 12;
   leap_sec_y[14] = 1989;
   leap_sec_m[14] = 12;
   leap_sec_y[15] = 1990;
   leap_sec_m[15] = 12;
   leap_sec_y[16] = 1992;
   leap_sec_m[16] = 6;
   leap_sec_y[17] = 1993;
   leap_sec_m[17] = 6;
   leap_sec_y[18] = 1994;
   leap_sec_m[18] = 6;
   leap_sec_y[19] = 1995;
   leap_sec_m[19] = 12;
   leap_sec_y[20] = 1997;
   leap_sec_m[20] = 6;
   leap_sec_y[21] = 1998;
   leap_sec_m[21] = 12;
   leap_sec_y[22] = 2005;
   leap_sec_m[22] = 12;
   leap_sec_y[23] = 2008;
   leap_sec_m[23] = 12;
   leap_sec_y[24] = 2012;
   leap_sec_m[24] = 6;
   int leaps[7];
   leaps[2] = 30;
   leaps[3] = 23;
   leaps[4] = 59;
   leaps[5] = 60;
   leaps[6] = 00;
   int start[7], end[7];
   if( utccompare(utc1,utc2) <= 0 )
      for( int c=0 ; c < 7 ;c++ )
      {
	 start[c] = utc1[c];
         end[c]   = utc2[c];
      }
   else
      for( int c=0 ; c < 7 ;c++ )
      {
	 start[c] = utc2[c];
         end[c]   = utc1[c];
      }

   int l, lstart, lend;
   l = 0;
   leaps[0] = leap_sec_y[l];
   leaps[1] = leap_sec_m[l];
   leaps[2] = leap_sec_m[l] == 6 ? 30 : 31;
   while( l <= nls-1 && utccompare(start,leaps) >= 0 )
   {
      l++;
      if( l <= nls-1 )
      {
	 leaps[0] = leap_sec_y[l];
	 leaps[1] = leap_sec_m[l];
	 leaps[2] = leap_sec_m[l] == 6 ? 30 : 31;
      }
   }
   lstart = l;
   //   cout << "lstart = " << lstart << endl;
   l = nls-1;
   leaps[0] = leap_sec_y[l];
   leaps[1] = leap_sec_m[l];
   leaps[2] = leap_sec_m[l] == 6 ? 30 : 31;
   while( l >= 0 && utccompare(end,leaps) <= 0 )
   {
      l--;
      if( l >= 0 )
      {
	 leaps[0] = leap_sec_y[l];
	 leaps[1] = leap_sec_m[l];
	 leaps[2] = leap_sec_m[l] == 6 ? 30 : 31;
      }
   }
   lend = l;
   //   cout << "lend = " << lend << endl;
   int corr = lend-lstart+1;
   delete[] leap_sec_y;
   delete[] leap_sec_m;
   return corr;
}

//-----------------------------------------------------------------------
void TimeSeries::filter_data( Filter* filter_ptr )
{
   if( m_myPoint )
   {
      filter_ptr->computeSOS( m_dt );

      filter_ptr->evaluate( mLastTimeStep+1, &mRecordedSol[0][0], &mRecordedSol[0][0] );
      filter_ptr->evaluate( mLastTimeStep+1, &mRecordedSol[1][0], &mRecordedSol[1][0] );
      filter_ptr->evaluate( mLastTimeStep+1, &mRecordedSol[2][0], &mRecordedSol[2][0] );

// Give the source time function a smooth start if this is a 2-pass (forward + backward) bandpass filter
      if( filter_ptr->get_passes() == 2 && filter_ptr->get_type() == bandPass )
      {    
	 double wghv, xi;
	 int p0=3, p=20 ; // First non-zero time level, and number of points in ramp;

	 for( int i=1 ; i<=p0-1 ; i++ )
	 {
	    mRecordedSol[0][i-1] = 0;
	    mRecordedSol[1][i-1] = 0;
	    mRecordedSol[2][i-1] = 0;
	 }
	 for( int i=p0 ; i<=p0+p ; i++ )
	 {
	    wghv = 0;
	    xi = (i-p0)/((double) p);
	 // polynomial P(xi), P(0) = 0, P(1)=1
	    wghv = xi*xi*xi*xi*(35-84*xi+70*xi*xi-20*xi*xi*xi);
	    mRecordedSol[0][i-1] *= wghv;
	    mRecordedSol[1][i-1] *= wghv;
	    mRecordedSol[2][i-1] *= wghv;
	 }
      }
   }
}

//-----------------------------------------------------------------------
void TimeSeries::print_timeinfo() const
{
   if( m_myPoint )
   {
      cout << "Observation/TimeSeries at grid point " << m_i0 << " " << m_j0 << " " << m_k0 << endl;
      cout << "   t0 = " << m_t0 << " dt= " << m_dt << endl;
      cout << "   Observation interval  [ " << m_t0 << " , " << m_t0 + m_dt*mLastTimeStep << " ] simulation time " << endl;
      if( m_utc_set )
      {
        printf("   Observation reference UTC  %02i/%02i/%i:%i:%i:%i.%03i\n", m_utc[1], m_utc[2], m_utc[0], m_utc[3],
	       m_utc[4], m_utc[5], m_utc[6] );

     }
      else
	 cout << "   Observation UTC reference time not defined" << endl;
   }
}

//-----------------------------------------------------------------------
void TimeSeries::set_window( double winl, double winr )
{
   m_use_win = true;
   m_winL = winl;
   m_winR = winr;
}

//-----------------------------------------------------------------------
void TimeSeries::exclude_component( bool usex, bool usey, bool usez )
{
   m_use_x = usex;
   m_use_y = usey;
   m_use_z = usez;
}

//-----------------------------------------------------------------------
void TimeSeries::readSACfiles( EW *ew, double timeshift, const char* sac1,
			       const char* sac2, const char* sac3 )
{
   string file1, file2, file3;
   if( ew->getObservationPath() != "./" )
   {
      file1 += ew->getObservationPath();
      file2 += ew->getObservationPath();
      file3 += ew->getObservationPath();
   }
   file1 += sac1;
   file2 += sac2;
   file3 += sac3;

   if( m_myPoint )
   {
      bool debug = false;
      double dt1, dt2, dt3, t01, t02, t03, lat1, lat2, lat3, lon1, lon2, lon3;
      double cmpaz1, cmpaz2, cmpaz3, cmpinc1, cmpinc2, cmpinc3;
      int utc1[7], utc2[7], utc3[7], npts1, npts2, npts3;

// Read header information
      readSACheader( file1.c_str(), dt1, t01, lat1, lon1, cmpaz1, cmpinc1, utc1, npts1 );
      readSACheader( file2.c_str(), dt2, t02, lat2, lon2, cmpaz2, cmpinc2, utc2, npts2 );
      readSACheader( file3.c_str(), dt3, t03, lat3, lon3, cmpaz3, cmpinc3, utc3, npts3 );

// Check that files are consistent with each other
      int eflag = 0;
      if( dt1 != dt2 || dt1 != dt3 || dt2 != dt3 )
         eflag = 1;
      if( t01 != t02 || t01 != t03 || t02 != t03 )
         eflag = 1;
      if( lat1 != lat2 || lat1 != lat3 || lat2 != lat3 )
         eflag = 1;
      if( lon1 != lon2 || lon1 != lon3 || lon2 != lon3 )
         eflag = 1;
      if( npts1 != npts2 || npts1 != npts3 || npts2 != npts3 )
         eflag = 1;

      bool utcequal=true;
      for( int c=0 ; c < 7 ; c++ )
	 if( utc1[c] != utc2[c] )
	    utcequal = false;
      for( int c=0 ; c < 7 ; c++ )
	 if( utc1[c] != utc3[c] )
	    utcequal = false;
      for( int c=0 ; c < 7 ; c++ )
	 if( utc2[c] != utc3[c] )
	    utcequal = false;
      if( !utcequal )
	 eflag = 1;
      
      if( eflag == 0 )
// Headers are ok, get the data
      {
	 // Check that all data are available
         
         bool azfail = false, incfail = false;
         if( cmpaz1 == -12345 || cmpaz2 == -12345 || cmpaz3 == -12345 )
            azfail = true;
         if( cmpinc1 == -12345 || cmpinc2 == -12345 || cmpinc3 == -12345 )
            incfail = true;

         if( !azfail && !incfail )
	 {
	    double* u1 = new double[npts1];
	    double* u2 = new double[npts1];
	    double* u3 = new double[npts1];
	    readSACdata( file1.c_str(), npts1, u1 );
	    readSACdata( file2.c_str(), npts1, u2 );
	    readSACdata( file3.c_str(), npts1, u3 );

	    allocateRecordingArrays( npts1, t01+timeshift, dt1 );

	    m_utc_set = true;
	    for( int c=0 ; c < 7 ; c++ )
	       m_utc[c] = utc1[c];

	    if( debug )
	    {
	       cout << "Read sac files " << file1 << " " << file2 << " " << file3 << endl;
	       cout << "UTC = " << utc1[1] << "/" << utc1[2] << "/" << utc1[0] << ":" << utc1[3]
		    << ":" << utc1[4] << ":" << utc1[5] << "." << utc1[6] << endl;
	       cout << " lat = " << lat1 << " lon = " << lon1 << endl;
	       cout << " dt = " << dt1 << " t0= " << t01  << " npts = " << npts1 << endl;
	       cout << " az1 = " << cmpaz1 << " inc1 = " << cmpinc1 << endl;
	       cout << " az2 = " << cmpaz2 << " inc2 = " << cmpinc2 << endl;
	       cout << " az3 = " << cmpaz3 << " inc3 = " << cmpinc3 << endl;
	    }
// Assume that we are using geographic coordinates, transform to (east,north,up) components.
	    const double convfactor = M_PI/180.0;
	    cmpaz1  *= convfactor;
	    cmpaz2  *= convfactor;
	    cmpaz3  *= convfactor;
	    cmpinc1 *= convfactor;
	    cmpinc2 *= convfactor;
	    cmpinc3 *= convfactor;

// Convert from station azimut to (e,n,u) components
	    double tmat[9];
	    tmat[0] = sin(cmpinc1)*cos(cmpaz1);
	    tmat[1] = sin(cmpinc2)*cos(cmpaz2);
	    tmat[2] = sin(cmpinc3)*cos(cmpaz3);
	    tmat[3] = sin(cmpinc1)*sin(cmpaz1);
	    tmat[4] = sin(cmpinc2)*sin(cmpaz2);
	    tmat[5] = sin(cmpinc3)*sin(cmpaz3);
	    tmat[6] = cos(cmpinc1);
	    tmat[7] = cos(cmpinc2);
	    tmat[8] = cos(cmpinc3);

	    m_xyzcomponent = false; //note this is format on output file, 
 	                         //internally, we always use (x,y,z) during computation.

// Convert (e,n,u) to (x,y,z) components.
	    double deti = 1.0/(m_thynrm*m_calpha+m_thxnrm*m_salpha);
	    double a11 = m_calpha*deti;
	    double a12 = m_thxnrm*deti;
	    double a21 =-m_salpha*deti;
	    double a22 = m_thynrm*deti;
	    for( int i=0 ; i < npts1 ; i++ )
	    {
	       double ncomp = tmat[0]*u1[i] + tmat[1]*u2[i] + tmat[2]*u3[i];
	       double ecomp = tmat[3]*u1[i] + tmat[4]*u2[i] + tmat[5]*u3[i];
	       double ucomp = tmat[6]*u1[i] + tmat[7]*u2[i] + tmat[8]*u3[i];
	       mRecordedSol[0][i] = a11*ncomp + a12*ecomp;
	       mRecordedSol[1][i] = a21*ncomp + a22*ecomp;
	       mRecordedSol[2][i] = -ucomp;
	    }
	    mLastTimeStep = npts1-1;
	 }
	 else
	 {
	    cout << "readSACfile, ERROR: no information about ";
            if( azfail )
	       cout << "component azimut ";
	    if( incfail )
	       cout << "component inclination ";
	    cout << " found on sac file" << endl;
  	    cout << "  station not read " << endl;
	 }
      }
      else
      {
         cout << "readSACfile, ERROR: found inconsistent meta data for files " << file1 << ", "
	      << file2 << ", " << file3 << endl;
	 cout << "  station not read " << endl;
         cout << "dt = " << dt1 << " " << dt2 << " " << dt3 << endl;
         cout << "t0 = " << t01 << " " << t02 << " " << t03 << endl;
         cout << "lat= " << lat1 << " " << lat2 << " " << lat3 << endl;
         cout << "lon= " << lon1 << " " << lon2 << " " << lon3 << endl;
         cout << "npt= " << npts1 << " " << npts2 << " " << npts3 << endl;
	 cout << "utc1 = " ;
	 for( int c=0 ; c < 7 ; c++ )
	    cout << utc1[c] << " ";
	 cout << endl;
	 cout << "utc2 = " ;
	 for( int c=0 ; c < 7 ; c++ )
	    cout << utc2[c] << " ";
	 cout << endl;
	 cout << "utc3 = " ;
	 for( int c=0 ; c < 7 ; c++ )
	    cout << utc3[c] << " ";
	 cout << endl;
      }
   }
}

//-----------------------------------------------------------------------
void TimeSeries::readSACheader( const char* fname, double& dt, double& t0,
				double& lat, double& lon, double& cmpaz,
				double& cmpinc, int utc[7], int& npts )
{

   float float70[70];
   int int35[35], logical[5];
   char kvalues[192];

   if( !(sizeof(float)==4) || !(sizeof(int)==4) || !(sizeof(char)==1) )
   {
      cout << "readSACheader: ERROR, size of datatypes do not match the SAC specification. Can not read SAC file "
	   << fname << endl;
      return;
   }

// Open SAC file
   FILE* fd=fopen(fname,"r");
   if( fd == NULL )
   {
      cout << "readSACheader: ERROR, observed data file " << fname << " could not be opened" << endl;
      return;
   }

// Read header data blocks
   size_t nr = fread(float70, sizeof(float), 70, fd );
   if( nr != 70 )
   {
      cout << "readSACheader: ERROR, could not read float part of header of " << fname << endl;
      fclose(fd);
      return;
   }
   nr = fread(int35, sizeof(int), 35, fd );
   if( nr != 35 )
   {
      cout << "readSACheader: ERROR, could not read int part of header of " << fname << endl;
      fclose(fd);
      return;
   }
   nr = fread(logical, sizeof(int), 5, fd );   
   if( nr != 5 )
   {
      cout << "readSACheader: ERROR, could not read bool part of header of " << fname << endl;
      fclose(fd);
      return;
   }
   nr = fread(kvalues, sizeof(char), 192, fd );   
   if( nr != 192 )
   {
      cout << "readSACheader: ERROR, could not read character part of header of " << fname << endl;
      fclose(fd);
      return;
   }

// Take out wanted information
   dt     = float70[0];
   t0     = float70[5];
   lat    = float70[31];
   lon    = float70[32];
   cmpaz  = float70[57];
   cmpinc = float70[58];
   utc[0] = int35[0];
   int jday=int35[1];
   convertjday( jday, utc[0], utc[2], utc[1] );
   utc[3] = int35[2];
   utc[4] = int35[3];
   utc[5] = int35[4];
   utc[6] = int35[5];
   npts   = int35[9];
}

//-----------------------------------------------------------------------
void TimeSeries::readSACdata( const char* fname, int npts, double* u )
{
   if( !(sizeof(float)==4) || !(sizeof(int)==4) || !(sizeof(char)==1) )
   {
      cout << "readSACdata: ERROR, size of datatypes do not match the SAC specification. Can not read SAC file "
	   << fname << endl;
      return;
   }

// Open SAC file
   FILE* fd=fopen(fname,"r");
   if( fd == NULL )
   {
      cout << "readSACdata: ERROR, observed data file " << fname << " could not be opened" << endl;
      return;
   }

// Skip header
   int retcode = fseek(fd,158*4,SEEK_SET);
   if( retcode != 0 )
   {
      cout << "readSACdata: ERROR, could not skip header in file " << fname << endl;
      fclose(fd);
      return;
   }

// Read data
   float* uf = new float[npts];
   size_t nr = fread( uf, sizeof(float), npts, fd );
   if( nr != npts )
   {
      cout << "readSACdata: ERROR, could not read float array of " << fname << endl;
      delete[] uf;
      fclose(fd);
      return;
   }

// Return floats as doubles
   for( int i=0 ; i < npts ; i++ )
      u[i] = static_cast<double>(uf[i]);
   delete[] uf;
   fclose(fd);
}

//-----------------------------------------------------------------------
void TimeSeries::convertjday( int jday, int year, int& day, int& month )
{
   if( jday > 0 && jday < 367 )
   {
      day = 1;
      int jd  = 1;
      month = 1;
      while( jd < jday )
      {
	 jd++;
         day++;
	 if( day > lastofmonth(year,month) )
	 {
	    day = 1;
	    month++;
	 }
      }
   }
   else
      cout << "convertjday: ERROR, jday is outside range " << endl;
}
