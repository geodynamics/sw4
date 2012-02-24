// -*-c++-*-
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "TimeSeries.h"
#include "mpi.h"
#include "sacsubc.h"
#include "csstime.h"

#include "Require.h"


#include "EW.h"

using namespace std;

TimeSeries::TimeSeries( EW* a_ew, std::string fileName, receiverMode mode, bool sacFormat, bool usgsFormat, 
			double x, double y, double depth, bool topoDepth, int writeEvery ):
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
  mBinaryMode(true)
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
   
   if (a_ew->getVerbosity()>=1)
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

   // m_calpha = cos(M_PI*a_ew->mGeoAz/180.0);
   // m_salpha = sin(M_PI*a_ew->mGeoAz/180.0);

   // double cphi   = cos(M_PI*m_lat/180.0);
   // double sphi   = sin(M_PI*m_lat/180.0);
   // double metersperdegree = a_ew->mMetersPerDegree;

   // m_thxnrm = m_salpha + (mX*m_salpha+mY*m_calpha)/cphi/metersperdegree * (M_PI/180.0) * sphi * m_calpha;
   // m_thynrm = m_calpha - (mX*m_salpha+mY*m_calpha)/cphi/metersperdegree * (M_PI/180.0) * sphi * m_salpha;
   // double nrm = sqrt( m_thxnrm*m_thxnrm + m_thynrm*m_thynrm );
   // m_thxnrm /= nrm;
   // m_thynrm /= nrm;
   // m_dthi = 1.0/(2*a_ew->mDt);

} // end constructor

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
// This routine only knows how to push the 3 doubles on the array stack.
// The calling routine need to figure out what needs to be saved
// and do any necessary pre-calculations
// ---------------------------------------------------------------

   mLastTimeStep++;
   if (mLastTimeStep < mAllocatedSize)
   {
     for (int q=0; q<m_nComp; q++)
       mRecordedSol[q][mLastTimeStep] = u[q];
     if (m_sacFormat)
     {
       for (int q=0; q<m_nComp; q++)
	 mRecordedFloats[q][mLastTimeStep] = (float) u[q];
     }     
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
void TimeSeries::writeFile()
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
  filePrefix << m_fileName << "." ;
  
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
     else if( m_mode == Velocity )
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
  setnhv( nm[10], mEventYear, nerr);
  setnhv( nm[11], mEventDay, nerr);
  setnhv( nm[12], mEventHour, nerr);
  setnhv( nm[13], mEventMinute, nerr);
  setnhv( nm[14], static_cast<int>(mEventSecond), nerr);
  setnhv( nm[15], 0, nerr);

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


void TimeSeries::write_usgs_format(string a_fileName)
{
  string mname[] = {"Zero","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
   FILE *fd=fopen(a_fileName.c_str(),"w");
   double lat, lon;
   double x, y, z;

// frequency resolution
//    double freq_limit=-999;
//    if (a_ew->m_prefilter_sources)
//      freq_limit = a_ew->m_fc;
//    else if (a_ew->m_limit_frequency)
//      freq_limit = a_ew->m_frequency_limit;

// write the header
   fprintf(fd, "# Author: SBP4W\n");
   fprintf(fd, "# Scenario: %s\n", "test"/*a_ew->m_scenario.c_str()*/);
   fprintf(fd, "# Date: %i-%s-%i\n", mEventDay, mname[mEventMonth].c_str(), mEventYear);
   fprintf(fd, "# Bandwith (Hz): %e\n", 1.234 /*freq_limit*/);
   fprintf(fd, "# Station: %s\n", m_fileName.c_str() /*mStationName.c_str()*/ );
   fprintf(fd, "# Target location (WGS84 longitude, latitude) (deg): %e %e\n", m_rec_lon, m_rec_lat);
   fprintf(fd, "# Actual location (WGS84 longitude, latitude) (deg): %e %e\n", m_rec_gp_lon, m_rec_gp_lat);
// distance in horizontal plane
   fprintf(fd, "# Distance from target to actual location (m): %e\n", sqrt( (mX-mGPX)*(mX-mGPX)+(mY-mGPY)*(mY-mGPY) ) );
   fprintf(fd, "# nColumns: %i\n", m_nComp+1);
   
   fprintf(fd, "# Column 1: Time (s)\n");
   if (m_mode == Displacement)
   {
     fprintf(fd, "# Column 2: X displacement (m)\n");
     fprintf(fd, "# Column 3: Y displacement (m)\n");
     fprintf(fd, "# Column 4: Z displacement (m)\n");
   }
   else if( m_mode == Velocity )
   {
     fprintf(fd, "# Column 2: X velocity (m/s)\n");
     fprintf(fd, "# Column 3: Y velocity (m/s)\n");
     fprintf(fd, "# Column 4: Z velocity (m/s)\n");
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

   for( int i = 0 ; i <= mLastTimeStep ; i++ )
   {
     fprintf(fd, "%e", m_t0 + i*m_dt);
     for (int q=0; q<m_nComp; q++)
       fprintf(fd, " %e", mRecordedSol[q][i]);
     fprintf(fd, "\n");
   }
   
   fclose(fd);
   
}
