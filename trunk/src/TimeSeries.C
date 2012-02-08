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
  mIgnore(true) // are we still using this flag???
{
// preliminary determination of nearest grid point ( before topodepth correction to mZ)
   a_ew->computeNearestGridPoint(m_i0, m_j0, m_k0, m_grid0, mX, mY, mZ);

// does this processor write this station?
   m_myPoint = a_ew->interior_point_in_proc(m_i0, m_j0, m_grid0);

// tmp
   printf("TimeSeries constructor, rank=%i, myPoint=%i\n", a_ew->getRank(), m_myPoint);

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
  if (m_mode == Solution)
    m_nComp=3;
  else if (m_mode == Div)
    m_nComp=1;
  else if (m_mode == Curl)
    m_nComp=3;
  else if (m_mode == Strains)
    m_nComp=6;

// allocate pointers to solution arrya pointer
  mRecordedSol = new double* [m_nComp];
  
  for (int q=0; q<m_nComp; q++)
    mRecordedSol[q] = NULL;
  
// do some misc pre computations
   // a_ew->computeGeographicCoord(mX, mY, m_lon, m_lat);

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
// reallocate the recording arrays
  if (mRecordedSol)
  {
    for (int q=0; q<m_nComp; q++)
    {
      if (mRecordedSol[q])
	delete [] mRecordedSol[q];
    }
    delete [] mRecordedSol;
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
     printf("Error: TimeSeries::recordData: passing a vector of size=%i but nComp=%i\n", u.size(), m_nComp);
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
   }
   else
   {
     printf("Ran out of recording space for the receiver station at (i,j,k,grid) = (%i, %i, %i, %i)\n",
	    m_i0, m_j0, m_k0, m_grid0);
     return;
   }

   if (mLastTimeStep > 0 && mLastTimeStep % mWriteEvery == 0)
     writeFile();
}

   

//
// the following is from the old SAC class
//

//    bool curvilinear = a_ew->topographyExists() && m_grid0 == a_ew->mNumberOfGrids-1;

//     // Time to record data for sac file
//    if( !m_div && !m_curl && !m_strains )
//    {
//       if( m_xycomponent )
//       {
// 	// Cartesian components
// 	 mRecordedUX.push_back((float)(a_ew->mU)[m_grid0](1,m_i0,m_j0,m_k0));
// 	 mRecordedUY.push_back((float)(a_ew->mU)[m_grid0](2,m_i0,m_j0,m_k0));
// 	 mRecordedUZ.push_back((float)(a_ew->mU)[m_grid0](3,m_i0,m_j0,m_k0));
//       }
//       else
//       {
// 	// North-South, East-West, and Up components
// 	 double ux  = (a_ew->mU)[m_grid0](1,m_i0,m_j0,m_k0);
// 	 double uy  = (a_ew->mU)[m_grid0](2,m_i0,m_j0,m_k0);

// 	 double uns = m_thynrm*ux-m_thxnrm*uy;
// 	 double uew = m_salpha*ux+m_calpha*uy;
	 
// 	 mRecordedUX.push_back( (float)(uew) ); //E-W is stored in UX
// 	 mRecordedUY.push_back( (float)(uns) ); //N-S is stored in UY
// 	 mRecordedUZ.push_back(-(float)(a_ew->mU)[m_grid0](3,m_i0,m_j0,m_k0));
//       }
//    }
//    // For curl and div, don't worry about one sided formulas at boundaries, because of ghost points.
//    else if( m_div && !curvilinear )
//    {

//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
//       double factor = 1.0/(2*a_ew->mGridSize[g]);
//       mRecordedUX.push_back((a_ew->mU[g](1,i+1, j, k) - a_ew->mU[g](1,i-1, j, k)+
// 			     a_ew->mU[g](2,i, j+1, k) - a_ew->mU[g](2,i, j-1, k)+
// 			     a_ew->mU[g](3,i, j, k+1) - a_ew->mU[g](3,i, j, k-1))*factor);
//       // Record dummy variables ==> can keep the velocity computation below unchanged.
//       mRecordedUY.push_back(0);
//       mRecordedUZ.push_back(0);
//    }
//    else if( m_div && curvilinear )
//    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
//       double factor = 1.0/(2);
//       mRecordedUX.push_back( ( a_ew->mQ(1,i,j,k)*(a_ew->mU[g](1,i+1,j,k) - a_ew->mU[g](1,i-1,j,k))+
// 			       a_ew->mQ(2,i,j,k)*(a_ew->mU[g](2,i+1,j,k) - a_ew->mU[g](2,i-1,j,k))+
// 			       a_ew->mQ(3,i,j,k)*(a_ew->mU[g](3,i+1,j,k) - a_ew->mU[g](3,i-1,j,k))+
// 			       a_ew->mR(1,i,j,k)*(a_ew->mU[g](1,i,j+1,k) - a_ew->mU[g](1,i,j-1,k))+
// 			       a_ew->mR(2,i,j,k)*(a_ew->mU[g](2,i,j+1,k) - a_ew->mU[g](2,i,j-1,k))+
// 			       a_ew->mR(3,i,j,k)*(a_ew->mU[g](3,i,j+1,k) - a_ew->mU[g](3,i,j-1,k))+
// 			       a_ew->mS(1,i,j,k)*(a_ew->mU[g](1,i,j,k+1) - a_ew->mU[g](1,i,j,k-1))+
// 			       a_ew->mS(2,i,j,k)*(a_ew->mU[g](2,i,j,k+1) - a_ew->mU[g](2,i,j,k-1))+
// 			       a_ew->mS(3,i,j,k)*(a_ew->mU[g](3,i,j,k+1) - a_ew->mU[g](3,i,j,k-1))  )*factor);
//       // Record dummy variables ==> can keep the velocity computation below unchanged.
//       mRecordedUY.push_back(0);
//       mRecordedUZ.push_back(0);
//    }
//    else if( m_curl && !curvilinear )
//    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
//       double factor = 1.0/(2*a_ew->mGridSize[g]);
//       double duydx = (a_ew->mU[g](2,i+1,j,k) - a_ew->mU[g](2,i-1,j,k))*factor;
//       double duzdx = (a_ew->mU[g](3,i+1,j,k) - a_ew->mU[g](3,i-1,j,k))*factor;
//       double duxdy = (a_ew->mU[g](1,i,j+1,k) - a_ew->mU[g](1,i,j-1,k))*factor;
//       double duzdy = (a_ew->mU[g](3,i,j+1,k) - a_ew->mU[g](3,i,j-1,k))*factor;
//       double duxdz = (a_ew->mU[g](1,i,j,k+1) - a_ew->mU[g](1,i,j,k-1))*factor;
//       double duydz = (a_ew->mU[g](2,i,j,k+1) - a_ew->mU[g](2,i,j,k-1))*factor;
//       if( m_xycomponent )
//       {
// 	 mRecordedUX.push_back( duzdy-duydz );
// 	 mRecordedUY.push_back( duxdz-duzdx );
// 	 mRecordedUZ.push_back( duydx-duxdy );
//       }
//       else
//       {
// 	 double uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
// 	 double uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
// 	 mRecordedUX.push_back( uew );
// 	 mRecordedUY.push_back( uns );
// 	 mRecordedUZ.push_back( -(duydx-duxdy) );
//       }
//    }
//    else if( m_curl && curvilinear )
//    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
//       double factor = 1.0/(2);
//       double duxdq = (a_ew->mU[g](1,i+1,j,k) - a_ew->mU[g](1,i-1,j,k));
//       double duydq = (a_ew->mU[g](2,i+1,j,k) - a_ew->mU[g](2,i-1,j,k));
//       double duzdq = (a_ew->mU[g](3,i+1,j,k) - a_ew->mU[g](3,i-1,j,k));
//       double duxdr = (a_ew->mU[g](1,i,j+1,k) - a_ew->mU[g](1,i,j-1,k));
//       double duydr = (a_ew->mU[g](2,i,j+1,k) - a_ew->mU[g](2,i,j-1,k));
//       double duzdr = (a_ew->mU[g](3,i,j+1,k) - a_ew->mU[g](3,i,j-1,k));
//       double duxds = (a_ew->mU[g](1,i,j,k+1) - a_ew->mU[g](1,i,j,k-1));
//       double duyds = (a_ew->mU[g](2,i,j,k+1) - a_ew->mU[g](2,i,j,k-1));
//       double duzds = (a_ew->mU[g](3,i,j,k+1) - a_ew->mU[g](3,i,j,k-1));
//       double duzdy = a_ew->mQ(2,i,j,k)*duzdq+a_ew->mR(2,i,j,k)*duzdr+a_ew->mS(2,i,j,k)*duzds;
//       double duydz = a_ew->mQ(3,i,j,k)*duydq+a_ew->mR(3,i,j,k)*duydr+a_ew->mS(3,i,j,k)*duyds;
//       double duxdz = a_ew->mQ(3,i,j,k)*duxdq+a_ew->mR(3,i,j,k)*duxdr+a_ew->mS(3,i,j,k)*duxds;
//       double duzdx = a_ew->mQ(1,i,j,k)*duzdq+a_ew->mR(1,i,j,k)*duzdr+a_ew->mS(1,i,j,k)*duzds;
//       double duydx = a_ew->mQ(1,i,j,k)*duydq+a_ew->mR(1,i,j,k)*duydr+a_ew->mS(1,i,j,k)*duyds;
//       double duxdy = a_ew->mQ(2,i,j,k)*duxdq+a_ew->mR(2,i,j,k)*duxdr+a_ew->mS(2,i,j,k)*duxds;
//       if( m_xycomponent )
//       {
// 	 mRecordedUX.push_back( (duzdy-duydz)*factor );
// 	 mRecordedUY.push_back( (duxdz-duzdx)*factor );
// 	 mRecordedUZ.push_back( (duydx-duxdy)*factor );
//       }
//       else
//       {
// 	 double uns = m_thynrm*(duzdy-duydz)-m_thxnrm*(duxdz-duzdx);
// 	 double uew = m_salpha*(duzdy-duydz)+m_calpha*(duxdz-duzdx);
// 	 mRecordedUX.push_back( uew*factor );
// 	 mRecordedUY.push_back( uns*factor );
// 	 mRecordedUZ.push_back( -(duydx-duxdy)*factor );
//       }
//    }
//    else if( m_strains && !curvilinear )
//    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
//       double factor = 1.0/(2*a_ew->mGridSize[g]);
//       double duydx = (a_ew->mU[g](2,i+1,j,k) - a_ew->mU[g](2,i-1,j,k))*factor;
//       double duzdx = (a_ew->mU[g](3,i+1,j,k) - a_ew->mU[g](3,i-1,j,k))*factor;
//       double duxdy = (a_ew->mU[g](1,i,j+1,k) - a_ew->mU[g](1,i,j-1,k))*factor;
//       double duzdy = (a_ew->mU[g](3,i,j+1,k) - a_ew->mU[g](3,i,j-1,k))*factor;
//       double duxdz = (a_ew->mU[g](1,i,j,k+1) - a_ew->mU[g](1,i,j,k-1))*factor;
//       double duydz = (a_ew->mU[g](2,i,j,k+1) - a_ew->mU[g](2,i,j,k-1))*factor;
//       double duxdx = (a_ew->mU[g](1,i+1,j,k) - a_ew->mU[g](1,i-1,j,k))*factor;
//       double duydy = (a_ew->mU[g](2,i,j+1,k) - a_ew->mU[g](2,i,j-1,k))*factor;
//       double duzdz = (a_ew->mU[g](3,i,j,k+1) - a_ew->mU[g](3,i,j,k-1))*factor;
//       mRecordedUX.push_back( duxdx );
//       mRecordedUY.push_back( duydy );
//       mRecordedUZ.push_back( duzdz );
//       mRecordedUXY.push_back( 0.5*(duydx+duxdy) );
//       mRecordedUXZ.push_back( 0.5*(duzdx+duxdz) );
//       mRecordedUYZ.push_back( 0.5*(duydz+duzdy) );
//    }
//    else if( m_strains && curvilinear )
//    {
//       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
//       double factor = 1.0/(2);
//       double duxdq = (a_ew->mU[g](1,i+1,j,k) - a_ew->mU[g](1,i-1,j,k));
//       double duydq = (a_ew->mU[g](2,i+1,j,k) - a_ew->mU[g](2,i-1,j,k));
//       double duzdq = (a_ew->mU[g](3,i+1,j,k) - a_ew->mU[g](3,i-1,j,k));
//       double duxdr = (a_ew->mU[g](1,i,j+1,k) - a_ew->mU[g](1,i,j-1,k));
//       double duydr = (a_ew->mU[g](2,i,j+1,k) - a_ew->mU[g](2,i,j-1,k));
//       double duzdr = (a_ew->mU[g](3,i,j+1,k) - a_ew->mU[g](3,i,j-1,k));
//       double duxds = (a_ew->mU[g](1,i,j,k+1) - a_ew->mU[g](1,i,j,k-1));
//       double duyds = (a_ew->mU[g](2,i,j,k+1) - a_ew->mU[g](2,i,j,k-1));
//       double duzds = (a_ew->mU[g](3,i,j,k+1) - a_ew->mU[g](3,i,j,k-1));
//       double duzdy = (a_ew->mQ(2,i,j,k)*duzdq+a_ew->mR(2,i,j,k)*duzdr+a_ew->mS(2,i,j,k)*duzds)*factor;
//       double duydz = (a_ew->mQ(3,i,j,k)*duydq+a_ew->mR(3,i,j,k)*duydr+a_ew->mS(3,i,j,k)*duyds)*factor;
//       double duxdz = (a_ew->mQ(3,i,j,k)*duxdq+a_ew->mR(3,i,j,k)*duxdr+a_ew->mS(3,i,j,k)*duxds)*factor;
//       double duzdx = (a_ew->mQ(1,i,j,k)*duzdq+a_ew->mR(1,i,j,k)*duzdr+a_ew->mS(1,i,j,k)*duzds)*factor;
//       double duydx = (a_ew->mQ(1,i,j,k)*duydq+a_ew->mR(1,i,j,k)*duydr+a_ew->mS(1,i,j,k)*duyds)*factor;
//       double duxdy = (a_ew->mQ(2,i,j,k)*duxdq+a_ew->mR(2,i,j,k)*duxdr+a_ew->mS(2,i,j,k)*duxds)*factor;
//       double duxdx = ( a_ew->mQ(1,i,j,k)*(a_ew->mU[g](1,i+1,j,k) - a_ew->mU[g](1,i-1,j,k))+
// 		       a_ew->mR(1,i,j,k)*(a_ew->mU[g](1,i,j+1,k) - a_ew->mU[g](1,i,j-1,k))+
// 		       a_ew->mS(1,i,j,k)*(a_ew->mU[g](1,i,j,k+1) - a_ew->mU[g](1,i,j,k-1)) )*factor;
//       double duydy = ( a_ew->mQ(2,i,j,k)*(a_ew->mU[g](2,i+1,j,k) - a_ew->mU[g](2,i-1,j,k))+
// 		       a_ew->mR(2,i,j,k)*(a_ew->mU[g](2,i,j+1,k) - a_ew->mU[g](2,i,j-1,k))+
// 		       a_ew->mS(2,i,j,k)*(a_ew->mU[g](2,i,j,k+1) - a_ew->mU[g](2,i,j,k-1)) )*factor;
//       double duzdz = ( a_ew->mQ(3,i,j,k)*(a_ew->mU[g](3,i+1,j,k) - a_ew->mU[g](3,i-1,j,k))+
// 		       a_ew->mR(3,i,j,k)*(a_ew->mU[g](3,i,j+1,k) - a_ew->mU[g](3,i,j-1,k))+
// 		       a_ew->mS(3,i,j,k)*(a_ew->mU[g](3,i,j,k+1) - a_ew->mU[g](3,i,j,k-1)) )*factor;
//       mRecordedUX.push_back( duxdx );
//       mRecordedUY.push_back( duydy );
//       mRecordedUZ.push_back( duzdz );
//       mRecordedUXY.push_back( 0.5*(duydx+duxdy) );
//       mRecordedUXZ.push_back( 0.5*(duzdx+duxdz) );
//       mRecordedUYZ.push_back( 0.5*(duydz+duzdy) );
//    }

//    if( m_velocities )
//    {
//      int n = mRecordedUX.size();
//      double vx, vy, vz;
//      if( n == 1 )
//      {
// 	m_dmx = mRecordedUX[n-1];
// 	m_dmy = mRecordedUY[n-1];
// 	m_dmz = mRecordedUZ[n-1];
//      }
//      else if( n == 2 )
//      {
// 	m_d0x = mRecordedUX[n-1];
// 	m_d0y = mRecordedUY[n-1];
// 	m_d0z = mRecordedUZ[n-1];
//      }
//      else if( n == 3 )
//      {
// // one sided formula for v(1)
//         vx = (-3*m_dmx+4*m_d0x-mRecordedUX[n-1] )*m_dthi; 
//         mRecordedUX[n-3] = vx;
//         vy = (-3*m_dmy+4*m_d0y-mRecordedUY[n-1] )*m_dthi; 
//         mRecordedUY[n-3] = vy;
//         vz = (-3*m_dmz+4*m_d0z-mRecordedUZ[n-1] )*m_dthi;
//         mRecordedUZ[n-3] = vz;
//      }
//      if( n >= 3 )
//      {
// // Standard case;
//         vx = (mRecordedUX[n-1]-m_dmx)*m_dthi;
//         vy = (mRecordedUY[n-1]-m_dmy)*m_dthi;
//         vz = (mRecordedUZ[n-1]-m_dmz)*m_dthi;
//         mRecordedUX[n-2] = vx;
//         mRecordedUY[n-2] = vy;
//         mRecordedUZ[n-2] = vz;

// // Need one sided forward, in case we are at the last point.
// // Otherwise, this is overwritten in the next time step.
//         vx = (m_dmx-4*m_d0x+3*mRecordedUX[n-1])*m_dthi;
//         m_dmx = m_d0x;
// 	m_d0x = mRecordedUX[n-1];
// 	mRecordedUX[n-1] = vx;
//         vy = (m_dmy-4*m_d0y+3*mRecordedUY[n-1])*m_dthi;
//         m_dmy = m_d0y;
// 	m_d0y = mRecordedUY[n-1];
// 	mRecordedUY[n-1] = vy;
//         vz = (m_dmz-4*m_d0z+3*mRecordedUZ[n-1])*m_dthi;
//         m_dmz = m_d0z;
// 	m_d0z = mRecordedUZ[n-1];
// 	mRecordedUZ[n-1] = vz;
//      }
//    }
//    if( m_velocities && m_strains )
//    {
//      int n = mRecordedUXY.size();
//      double vx, vy, vz;
//      if( n == 1 )
//      {
// 	m_dmxy = mRecordedUXY[n-1];
// 	m_dmxz = mRecordedUXZ[n-1];
// 	m_dmyz = mRecordedUYZ[n-1];
//      }
//      else if( n == 2 )
//      {
// 	m_d0xy = mRecordedUXY[n-1];
// 	m_d0xz = mRecordedUXZ[n-1];
// 	m_d0yz = mRecordedUYZ[n-1];
//      }
//      else if( n == 3 )
//      {
// // one sided formula for v(1)
//         vx = (-3*m_dmxy+4*m_d0xy-mRecordedUXY[n-1] )*m_dthi; 
//         mRecordedUXY[n-3] = vx;
//         vy = (-3*m_dmxz+4*m_d0xz-mRecordedUXZ[n-1] )*m_dthi; 
//         mRecordedUXZ[n-3] = vy;
//         vz = (-3*m_dmyz+4*m_d0yz-mRecordedUYZ[n-1] )*m_dthi;
//         mRecordedUYZ[n-3] = vz;
//      }
//      if( n >= 3 )
//      {
// // Standard case;
//         vx = (mRecordedUXY[n-1]-m_dmxy)*m_dthi;
//         vy = (mRecordedUXZ[n-1]-m_dmxz)*m_dthi;
//         vz = (mRecordedUYZ[n-1]-m_dmyz)*m_dthi;
//         mRecordedUXY[n-2] = vx;
//         mRecordedUXZ[n-2] = vy;
//         mRecordedUYZ[n-2] = vz;

// // Need one sided forward, in case we are at the last point.
// // Otherwise, this is overwritten in the next time step.
//         vx = (m_dmxy-4*m_d0xy+3*mRecordedUXY[n-1])*m_dthi;
//         m_dmxy = m_d0xy;
// 	m_d0xy = mRecordedUXY[n-1];
// 	mRecordedUXY[n-1] = vx;
//         vy = (m_dmxz-4*m_d0xz+3*mRecordedUXZ[n-1])*m_dthi;
//         m_dmxz = m_d0xz;
// 	m_d0xz = mRecordedUXZ[n-1];
// 	mRecordedUYZ[n-1] = vy;
//         vz = (m_dmyz-4*m_d0yz+3*mRecordedUYZ[n-1])*m_dthi;
//         m_dmyz = m_d0yz;
// 	m_d0yz = mRecordedUYZ[n-1];
// 	mRecordedUYZ[n-1] = vz;
//      }
//    }

//   if (mWriteEvery != 0 && mLastTimeStep % mWriteEvery == 0 && 
//      mLastTimeStep >= mWriteEvery) writeFile();

   
void TimeSeries::writeFile()
{
  if (!m_myPoint) return;

// ------------------------------------------------------------------
// We should add an argument to this function that describes how the
// header and filename should be constructed
// ------------------------------------------------------------------

  stringstream filePrefix;

//building the file name...
  filePrefix << m_fileName << "." ;
  
  stringstream ux, uy, uz, uxy, uxz, uyz;
  
// Write out displacement components (ux, uy, uz)

  if( m_sacFormat )
  {
     // string mode = "ASCII";
     // if (mBinaryMode)
     // 	mode = "BINARY";
     // inihdr();

     // stringstream msg;
     // msg << "Writing " << mode << " SAC Files, "
     // 	 << "of size " << mRecordedUX.size() << ": "
     // 	 << filePrefix.str();

     // string xfield, yfield, zfield, xyfield, xzfield, yzfield;
     // double azimx, azimy, updownang;
     // if( !m_div && !m_curl && !m_strains )
     // {
     // 	if( m_xycomponent && !m_velocities )
     // 	{
     // 	   xfield = "X";
     // 	   yfield = "Y";
     // 	   zfield = "Z";
     // 	   ux << filePrefix.str() << "x";
     // 	   uy << filePrefix.str() << "y";
     // 	   uz << filePrefix.str() << "z";
     // 	   azimx = a_ew->mGeoAz;
     // 	   azimy = a_ew->mGeoAz+90.;
     // 	   updownang = 180;
     // 	   msg << "[x|y|z]" << endl;
     // 	}
     // 	else if( m_xycomponent && m_velocities )
     // 	{
     // 	   xfield = "Vx";
     // 	   yfield = "Vy";
     // 	   zfield = "Vz";
     // 	   ux << filePrefix.str() << "xv";
     // 	   uy << filePrefix.str() << "yv";
     // 	   uz << filePrefix.str() << "zv";
     // 	   azimx = a_ew->mGeoAz;
     // 	   azimy = a_ew->mGeoAz+90.;
     // 	   updownang = 180;
     // 	   msg << "[xv|yv|zv]" << endl;
     // 	}
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
     // else if( m_div && !m_velocities )
     // {
     // 	xfield = "Div";
     // 	ux << filePrefix.str() << "div";
     // 	azimx = a_ew->mGeoAz;
     // 	azimy = a_ew->mGeoAz+90.;
     // 	updownang = 180;
     // 	msg << "[div]" << endl;
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
     // else if( m_curl && !m_velocities && m_xycomponent )
     // {
     // 	xfield = "Curlx";
     // 	yfield = "Curly";
     // 	zfield = "Curlz";
     // 	ux << filePrefix.str() << "curlx";
     // 	uy << filePrefix.str() << "curly";
     // 	uz << filePrefix.str() << "curlz";
     // 	azimx = a_ew->mGeoAz;
     // 	azimy = a_ew->mGeoAz+90.;
     // 	updownang = 180;
     // 	msg << "[curlx|curly|curlz]" << endl;
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
     // else if( m_strains )
     // {
     // 	xfield = "Uxx";
     // 	yfield = "Uyy";
     // 	zfield = "Uzz";
     // 	xyfield = "Uxy";
     // 	xzfield = "Uxz";
     // 	yzfield = "Uyz";
     // 	ux << filePrefix.str() << "xx";
     // 	uy << filePrefix.str() << "yy";
     // 	uz << filePrefix.str() << "zz";
     // 	uxy << filePrefix.str() << "xy";
     // 	uxz << filePrefix.str() << "xz";
     // 	uyz << filePrefix.str() << "yz";
     // 	azimx = a_ew->mGeoAz;
     // 	azimy = a_ew->mGeoAz+90.;
     // 	updownang = 180;
     // 	msg << "[xx|yy|zz|xy|xz|yz]" << endl;
     // }

     // cout << msg.str() ;
     // writeSACFile(mRecordedUX.size(), 
     //           const_cast<char*>(ux.str().c_str()), 
     //           &mRecordedUX[0], (float)(a_ew->mTstart - a_ew->m_t0Shift), (float)(a_ew->mDt),
     // 		  const_cast<char*>(xfield.c_str()), 90.0, azimx); 
     // if( !m_div )
     // {
     // 	writeSACFile(mRecordedUZ.size(), 
     //           const_cast<char*>(uz.str().c_str()), 
     //           &mRecordedUZ[0], (float)(a_ew->mTstart - a_ew->m_t0Shift), (float)(a_ew->mDt),
     //           const_cast<char*>(zfield.c_str()), updownang, 0.0);      
     // 	writeSACFile(mRecordedUY.size(), 
     //           const_cast<char*>(uy.str().c_str()), 
     //           &mRecordedUY[0], (float)(a_ew->mTstart - a_ew->m_t0Shift), (float)(a_ew->mDt),
     //           const_cast<char*>(yfield.c_str()), 90.0, azimy); 
     // }
     // if( m_strains )
     // {
     // 	writeSACFile(mRecordedUXY.size(), 
     // 		     const_cast<char*>(uxy.str().c_str()), 
     // 		     &mRecordedUXY[0], (float)(a_ew->mTstart - a_ew->m_t0Shift), (float)(a_ew->mDt),
     // 		  const_cast<char*>(xyfield.c_str()), 90.0, azimx); 
     // 	writeSACFile(mRecordedUXZ.size(), 
     // 		     const_cast<char*>(uxz.str().c_str()), 
     // 		     &mRecordedUXZ[0], (float)(a_ew->mTstart - a_ew->m_t0Shift), (float)(a_ew->mDt),
     // 		  const_cast<char*>(xzfield.c_str()), 90.0, azimx); 
     // 	writeSACFile(mRecordedUYZ.size(), 
     // 		     const_cast<char*>(uyz.str().c_str()), 
     // 		     &mRecordedUYZ[0], (float)(a_ew->mTstart - a_ew->m_t0Shift), (float)(a_ew->mDt),
     // 		  const_cast<char*>(yzfield.c_str()), 90.0, azimx); 

     // }

  }
  if( m_usgsFormat )
  {
    filePrefix << "txt";
    cout << "Writing receiever data file using USGS ascii format, "
    	 << "of size " << mLastTimeStep+1 << ": "
	 << filePrefix.str() << endl;

    write_usgs_format( filePrefix.str() );
  }

}

// void TimeSeries::
// writeSACFile(int npts, char *ofile, float *y, float btime, float dt, char *var,
//              float cmpinc, float cmpaz)
// {
//   /*
//     PURPOSE: WRITE AN EVENLY SPACED SAC FILE
    
//     	write a binary sac file with evenly sampled data
//     	ofile	Char	name of file
//     	y	R	array of values
//     	npts	I	number of points in data
//     	btime	R	start time
//     	dt	R	sample interval
//     	maxpts	I	maximum number of points to read
//     	nerr	I	error return
//     -----
//   */
//   float e;
//   float depmax, depmin, depmen;
//   int* nerr = 0;
// // assign all names in a string array
// //               0         1         2         3          4         5           6          7          8          9
//   string nm[]={"DEPMAX", "DEPMIN", "DEPMEN", "NPTS    ","DELTA   ","B       ", "E       ","LEVEN   ","LOVROK  ","LCALDA  ",
// //              10          11          12          13          14          15           16          17          18
// 	       "NZYEAR  ", "NZJDAY  ", "NZHOUR  ", "NZMIN   ", "NZSEC   ", "NZMSEC   ", "KCMPNM  ", "STLA    ", "STLO    ",
// //              19          20          21          22          23          24          25
// 	       "EVLA    ", "EVLO    ", "EVDP    ", "O       ", "CMPINC  ", "CMPAZ   ", "KSTNM   "
//   };

//   newhdr();
//   scmxmn(y,npts,&depmax,&depmin,&depmen);
// //  setfhv("DEPMAX", depmax, nerr);
//   setfhv((char *) nm[0].c_str(), depmax, nerr);
//   setfhv((char *) nm[1].c_str(), depmin, nerr);
//   setfhv((char *) nm[2].c_str(), depmen, nerr);
//   setnhv((char *) nm[3].c_str(), npts,nerr);
//   setfhv((char *) nm[4].c_str(), dt  ,nerr);
//   setfhv((char *) nm[5].c_str(), btime  ,nerr);
//   e = btime + (npts -1 )*dt;
//   setfhv((char *) nm[6].c_str(), e, nerr);
//   setlhv((char *) nm[7].c_str(), 1, nerr);
//   setlhv((char *) nm[8].c_str(), 1, nerr);
//   setlhv((char *) nm[9].c_str(), 1, nerr);

//   // setup time info
//   setnhv((char *) nm[10].c_str(), mEventYear, nerr);
//   setnhv((char *) nm[11].c_str(), mEventDay, nerr);
//   setnhv((char *) nm[12].c_str(), mEventHour, nerr);
//   setnhv((char *) nm[13].c_str(), mEventMinute, nerr);
//   setnhv((char *) nm[14].c_str(), static_cast<int>(mEventSecond), nerr);
//   setnhv((char *) nm[15].c_str(), 0, nerr);

//   // field we're writing
//   setkhv((char *) nm[16].c_str(), var, nerr);

//   // location of the receiver
// //   double lat, lon;
// //   a_ew->computeGeographicCoord(mX, mY, lon, lat); //(C.B: I think that this is the point we want)
//   setfhv((char *) nm[17].c_str(), m_lat, nerr);
//   setfhv((char *) nm[18].c_str(), m_lon, nerr);
//   // location of epicenter
//   setfhv((char *) nm[19].c_str(), mEpicenter.getLatitude(), nerr);
//   setfhv((char *) nm[20].c_str(), mEpicenter.getLongitude(), nerr);
//   setfhv((char *) nm[21].c_str(), mEpicenter.getDepth()/1000.0, nerr);
//   // time offset for epicenter source
//   setfhv((char *) nm[22].c_str(), mEpicenterTimeOffset, nerr);

//   // set inclination and azimuthal angle
//   setfhv((char *) nm[23].c_str(), cmpinc, nerr);
//   setfhv((char *) nm[24].c_str(), cmpaz, nerr);

//   // set the station name
//   setkhv((char *) nm[25].c_str(), const_cast<char*>(mStationName.c_str()), nerr);


//   if (!mBinaryMode)
//     awsac(npts, ofile, y);
//   else
//     bwsac(npts, ofile, y);
// }


void TimeSeries::write_usgs_format(string a_fileName)
{
   string mname[12] ={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
   FILE *fd=fopen(a_fileName.c_str(),"w");
   double lat, lon;
   double x, y, z;
//    a_ew->coord_from_index( m_i0, m_j0, m_kSAC, m_grid0, x, y, z );
//    a_ew->computeGeographicCoord( x, y, lon, lat );

//    double freq_limit=-999;
//    if (a_ew->m_prefilter_sources)
//      freq_limit = a_ew->m_fc;
//    else if (a_ew->m_limit_frequency)
//      freq_limit = a_ew->m_frequency_limit;

   fprintf(fd, "# Author: WPP\n");
   fprintf(fd, "# Scenario: %s\n", "test"/*a_ew->m_scenario.c_str()*/);
   fprintf(fd, "# Date: %i-%s-%i\n", 7, "Feb", 2012); //mTimePtr->tm_mday, mname[mTimePtr->tm_mon].c_str(), 1900+mTimePtr->tm_year
   fprintf(fd, "# Bandwith (Hz): %e\n", 1.234 /*freq_limit*/);
   fprintf(fd, "# Station: %s\n", m_fileName.c_str() /*mStationName.c_str()*/ );
   fprintf(fd, "# Target location (WGS84 longitude, latitude) (deg): %e %e\n", 0.0, 0.0); // m_lon, m_lat);
   fprintf(fd, "# Actual location (WGS84 longitude, latitude) (deg): %e %e\n", 0.0, 0.0); // lon, lat);
// distance in horizontal plane
   fprintf(fd, "# Distance from target to actual location (m): %e\n", 0.0); // sqrt( (x-mX)*(x-mX)+(y-mY)*(y-mY) ) );
   fprintf(fd, "# nColumns: %i\n", m_nComp+1);
   
   fprintf(fd, "# Column 1: Time (s)\n");
//    if( !m_div && !m_curl && !m_strains )
//    {
//       if( !m_xycomponent && m_velocities )
//       {
// 	 fprintf(fd, "# Column 2: East-west velocity (m/s)\n");
// 	 fprintf(fd, "# Column 3: North-south velocity (m/s)\n");
// 	 fprintf(fd, "# Column 4: Up-down velocity (m/s)\n");
//       }
//       else if( !m_xycomponent && !m_velocities )
//       {
// 	 fprintf(fd, "# Column 2: East-west displacement (m)\n");
// 	 fprintf(fd, "# Column 3: North-south displacement (m)\n");
// 	 fprintf(fd, "# Column 4: Up-down displacement (m)\n");
//       }
//       else if( m_xycomponent && m_velocities )
//       {
// 	 fprintf(fd, "# Column 2: X velocity (m/s)\n");
// 	 fprintf(fd, "# Column 3: Y velocity (m/s)\n");
// 	 fprintf(fd, "# Column 4: Z velocity (m/s)\n");
//       }
//       else
//       {
   if (m_mode == Solution)
   {
     fprintf(fd, "# Column 2: X displacement (m)\n");
     fprintf(fd, "# Column 3: Y displacement (m)\n");
     fprintf(fd, "# Column 4: Z displacement (m)\n");
   }
   else if( m_mode == Div )
   {
//       if( m_velocities )
// 	 fprintf(fd, "# Column 2: divergence of velocity (/s)\n");
//       else
     fprintf(fd, "# Column 2: divergence of displacement ()\n");
   }
   else if( m_mode == Curl )
   {
//       if( m_velocities )
//       {
// 	 fprintf(fd, "# Column 2: curl of velocity, component 1 (/s)\n");
// 	 fprintf(fd, "# Column 3: curl of velocity, component 2 (/s)\n");
// 	 fprintf(fd, "# Column 4: curl of velocity, component 3 (/s)\n");
//       }
//       else
//       {
     fprintf(fd, "# Column 2: curl of displacement, component 1 ()\n");
     fprintf(fd, "# Column 3: curl of displacement, component 2 ()\n");
     fprintf(fd, "# Column 4: curl of displacement, component 3 ()\n");
//       }
   }
   else if( m_mode == Strains )
   {
//       if( m_velocities )
//       {
// 	 fprintf(fd, "# Column 2: xx velocity strain component (/s)\n");
// 	 fprintf(fd, "# Column 3: yy velocity strain component (/s)\n");
// 	 fprintf(fd, "# Column 4: zz velocity strain component (/s)\n");
// 	 fprintf(fd, "# Column 5: xy velocity strain component (/s)\n");
// 	 fprintf(fd, "# Column 6: xz velocity strain component (/s)\n");
// 	 fprintf(fd, "# Column 7: yz velocity strain component (/s)\n");
//       }
//       else
//       {
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
