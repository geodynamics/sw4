//-*-c++-*-
#include "mpi.h"

#include "EW.h"
#include "Require.h"
#include "Image.h"
#include <cstdio>
#include <cmath>
#include <fcntl.h>

// initializing static member
int Image::mPreceedZeros=0;

//int Image::MODES=29;
int Image::MODES=7;

using namespace std;

Image::Image(EW * a_ew,
             float time, 
             float timeInterval, 
             int cycle, 
             int cycleInterval,
             const std::string& filePrefix, 
             int sample,
             ImageMode mode,
	     ImageOrientation locationType, 
             float locationValue,
	     bool doubleMode ):
  mTime(time),
  mEW(a_ew),
  m_time_done(false),
  mTimeInterval(timeInterval),
  mWritingCycle(cycle),
  mCycleInterval(cycleInterval),
  mFilePrefix(filePrefix),
  mImageSamplingFactor(sample),
  mMode(mode),
  mFileName(""),
  m_gridPtValueInitialized(false),
  //   mIOReady(false),
  mCycle(-1),
  mWriting(false),
  mReadyToWrite(false),
  mLocationType(locationType),
  mCoordValue(locationValue),
  m_isDefined(false),
  m_isDefinedMPIWriters(false),
  m_double(doubleMode)
{
  mMode2Suffix.resize(MODES);
  mMode2Suffix[NONE] = "unknown";
  mMode2Suffix[UX] = "ux";
  mMode2Suffix[UY] = "uy";
  mMode2Suffix[UZ] = "uz";
  mMode2Suffix[RHO] = "rho";
  mMode2Suffix[LAMBDA] = "lambda";
  mMode2Suffix[MU] = "mu";
  // mMode2Suffix[P] = "p";
  // mMode2Suffix[S] = "s";
  // mMode2Suffix[DIV] = "div";
  // mMode2Suffix[CURL] = "curl";
  // mMode2Suffix[VELDIV] = "veldiv";
  // mMode2Suffix[VELCURL] = "velcurl";
  // mMode2Suffix[LAT] = "lat";
  // mMode2Suffix[LON] = "lon";
  // mMode2Suffix[HVELMAX] = "hvelmax";
  // mMode2Suffix[VVELMAX] = "vvelmax";
  // mMode2Suffix[TOPO] = "topo";
  // mMode2Suffix[GRID] = "grid";
  // mMode2Suffix[UXERR] = "uxerr";
  // mMode2Suffix[UYERR] = "uyerr";
  // mMode2Suffix[UZERR] = "uzerr";
  // mMode2Suffix[FX] = "fx";
  // mMode2Suffix[FY] = "fy";
  // mMode2Suffix[FZ] = "fz";
  // mMode2Suffix[VELMAG] = "velmag";
  // mMode2Suffix[QS] = "qs";
  // mMode2Suffix[QP] = "qp";
  // mMode2Suffix[HVEL] = "hvel";

  mOrientationString.resize(4);
  mOrientationString[UNDEFINED] = "undefined";
  mOrientationString[X] = "x";
  mOrientationString[Y] = "y";
  mOrientationString[Z] = "z";

  initializeTime();

  mWindow.resize(mEW->mNumberOfGrids);
  for (int g=0; g<mEW->mNumberOfGrids; g++)
  {
    mWindow[g] = new int[6];
    mWindow[g][0] = 1;
    mWindow[g][1] = 0;
    mWindow[g][2] = 1;
    mWindow[g][3] = 0;
    mWindow[g][4] = 1;
    mWindow[g][5] = 0;
  }
  
  if (m_double)
  {
    m_doubleField.resize(mEW->mNumberOfGrids);
  }
  else
  {
    m_floatField.resize(mEW->mNumberOfGrids);
  }
  
}


//-----------------------------------------------------------------------
bool Image::proc_write()
{
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  return (myRank == m_rankWriter);
}

//-------------------------------------
void Image::setSteps(int a_steps)
{
  char buffer[50];
  mPreceedZeros = snprintf( buffer, 50, "%d", a_steps );
}

//-------------------------------------

/* Here, we compute the index --in the local grid-- of the coordinate value at which we have to plot. 
   For x, y, the local grid is the same as the global grid, but for z, k resets at each refinement boundary. */

// the purpose of this routine is to assign the vector<int> m_gridPtIndex
void Image::computeGridPtIndex()
{
//   ASSERT(m_isInitializedGridSize);
//   ASSERT(m_isInitializedZMin)    ;

// Seems to work as follows: 
//           For X and Y images,  m_gridPtIndex[g] is the coordinate index of the 
//              part of the image plane in grid g.
//           For Z images, m_gridPtIndex[0] is the coordinate index, and m_gridPtIndex[1] is
//              the grid no. of the image plane.
//
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//   if (myRank == 0)
//     printf("======== Initializing Image ==========\n");

  int n = mEW->mNumberOfCartesianGrids; 
  int nTotal = mEW->mNumberOfGrids;

  if (mLocationType == Image::X || mLocationType == Image::Y)
  { 
     m_gridPtIndex.resize(nTotal);
      
    /* I store the indices for i on the local grid of all levels: index(iCoarse,iFiner,iFinest) */

    for (int g = 0; g < n; g++)
    {
      if (g == 0) // We find the closest ***COARSE*** grid line
      {
	m_gridPtIndex[g] = (int)floor(mCoordValue/mEW->mGridSize[g])+1;
              
	if (mCoordValue-((m_gridPtIndex[g]-0.5)*mEW->mGridSize[g]) > 0.) (m_gridPtIndex[g])++;
      }
      else
      {
	m_gridPtIndex[g] = 2*m_gridPtIndex[g-1]-1;
      }
//       if (myRank == 0)
// 	printf("The closest grid line is located at %s = %.2f; index = %i on grid %i\n",
// 	       mOrientationString[mLocationType].c_str(), (m_gridPtIndex[g]-1)*mEW->mGridSize[g],m_gridPtIndex[g],g);
    }
// curvilinear grid on top: copy location from top Cartesian grid
    if (nTotal > n)
    {
      m_gridPtIndex[nTotal-1] = m_gridPtIndex[n-1];
    }
  } // end if X or Y
  else if (mLocationType == Image::Z)
  {
// I store the index as follows index(k,g): k in local reference frame and g as the grid it belongs to*/

    m_gridPtIndex.resize(2)  ;
    bool breakLoop = false;

// are we in the curvilinear grid?
    if (mCoordValue < mEW->m_zmin[n-1])
    {
// C.B> 09-10-08. Changed to plot the top surface
      m_gridPtIndex[0] = 1       ;
      m_gridPtIndex[1] = nTotal-1;
    }
    else
    {
      for (int g = 0; g < n; g++) // currently only looking through the Cartesian grids
      {
	if (mCoordValue > mEW->m_zmin[g])
	{
	  m_gridPtIndex[0] = (int)floor((mCoordValue-mEW->m_zmin[g])/mEW->mGridSize[g])+1;
	  m_gridPtIndex[1] = g                                   ;

	  if (mCoordValue-(mEW->m_zmin[g]+(m_gridPtIndex[0]-0.5)*mEW->mGridSize[g]) > 0.)  m_gridPtIndex[0]++;

	  breakLoop = true;
	}
	else if (mCoordValue == mEW->m_zmin[g])
	{
	  if (g == n-1)
	  {
	    m_gridPtIndex[0] = 1;
	    m_gridPtIndex[1] = g;
	  }
	  else
	  {
	    m_gridPtIndex[0] = (int)floor((mCoordValue-mEW->m_zmin[g+1])/mEW->mGridSize[g+1])+1; // Here, I know I am on a grid line
	    m_gridPtIndex[1] = g+1;
	  }
	  breakLoop = true;
	}

	if (breakLoop)
	{
	  break;
	} 
      }
    } // end if in the Cartesian
  } // end if Z
  
  int iwrite = plane_in_proc(m_gridPtIndex[0]) ? 1 : 0;
  
  MPI_Group origGroup, newGroup; 
  MPI_Comm_group(MPI_COMM_WORLD, &origGroup); 

//   int myRank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//   cout<<"myRank "<<myRank<<" and I write "<<iwrite<<endl;
//   MPI_Barrier;
  
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size); 
  std::vector<int> writers(size);
  MPI_Allgather(&iwrite, 1, MPI_INT, &writers[0], 1, MPI_INT,MPI_COMM_WORLD);
  std::vector<int> fileWriterIDs;
  for (unsigned int i = 0; i < writers.size(); ++i)
    if (writers[i] == 1)
      {
        fileWriterIDs.push_back(i);
      }

  ASSERT2(fileWriterIDs.size() > 0,
          "ImageSlice write error, no processors in communicator")          ;

  m_rankWriter = fileWriterIDs[0];

  MPI_Group_incl(origGroup,fileWriterIDs.size(),&fileWriterIDs[0],&newGroup);
  MPI_Comm_create(MPI_COMM_WORLD,newGroup,&m_mpiComm_writers)               ;
  
//   int newRank;
//   MPI_Group_rank(newGroup,&newRank);
//   MPI_Group_size(newGroup,&size);

  MPI_Group_free(&origGroup);
  MPI_Group_free(&newGroup) ;
  
//   getchar();
  
//  ASSERT(m_mpiComm_writers != MPI_COMM_NULL);
  
  m_isDefinedMPIWriters = true;
  m_gridPtValueInitialized = true;
}


//-----------------------------------------------------------------------
bool Image::plane_in_proc( int a_gridIndexCoarsest)
{
// Find intersection of image with local processor grid block, all computed
// in global indices.
  bool retval = false;
  int a_iStart = mEW->m_iStart[0];
  int a_iEnd = mEW->m_iEnd[0];
  int a_jStart =  mEW->m_jStart[0];
  int a_jEnd = mEW->m_jEnd[0];
  
  if (mLocationType == Image::X)
  {
    retval = (a_gridIndexCoarsest >= (a_iStart+mEW->m_paddingCells[0])) && (a_gridIndexCoarsest <= (a_iEnd-mEW->m_paddingCells[1]));
  }
  if (mLocationType == Image::Y)
  {
    retval = (a_gridIndexCoarsest >= (a_jStart+mEW->m_paddingCells[2])) && (a_gridIndexCoarsest <= (a_jEnd-mEW->m_paddingCells[3]));
  }
  if (mLocationType == Image::Z)
  {
    retval = true;
  }
  return retval;
}


//-------------------------------------
void Image::initializeTime()
{
  mNextTime = 0.; 
  m_time_done = false;
// with the option timeInterval=..., first time is always t=0
// is this what we want when m_t0Shift >> 0?
}

const std::string Image::fieldSuffix(ImageMode mode) const
{
  REQUIRE2(mode >= 0 && mode < MODES, "mode=" << mode << " out of bounds");
  return mMode2Suffix[mode];
}

//-----------------------------------------------------------------------
bool Image::timeToWrite(double time, int cycle, double dt )
{
  // -----------------------------------------------
  // Check based on cycle
  // -----------------------------------------------
  bool do_it=false;
  if( cycle == mWritingCycle )
    do_it = true;
  if( mCycleInterval !=  0 && cycle%mCycleInterval == 0 ) 
    do_it = true;
  // ---------------------------------------------------
  // Check based on time
  // ---------------------------------------------------
  if(mTime > 0.0 && (  mTime <= time + dt*0.5 ) && !m_time_done )
  {
     m_time_done = true;
     do_it = true;
  }
  if( mTimeInterval != 0.0 && mNextTime <= time + dt*0.5 )
  {
      // we're going to write it, so increase the time interval

    //     while( mNextTime < time )
	mNextTime += mTimeInterval;
     do_it =  true;
  }
  return do_it;
}

//-----------------------------------------------------------------------
void Image::define_pio( )
{
   int glow = 0, ghigh = mEW->mNumberOfGrids;
   if( mLocationType == Image::Z )
   {
      glow = m_gridPtIndex[1];
      ghigh = glow+1;
   }
   m_pio = new Parallel_IO*[ghigh-glow+1];
   for( int g=glow ; g < ghigh ; g++ )
   {
      int global[3] = {mEW->m_global_nx[g],
		       mEW->m_global_ny[g],
		       mEW->m_global_nz[g]} ;
//		       mEW->m_kEnd[g] - mEW->m_kStart[g] - 2*mEW->m_ghost_points+1} ;
      int local[3];
      local[0] = mWindow[g][1] - mWindow[g][0] + 1;
      local[1] = mWindow[g][3] - mWindow[g][2] + 1;
      local[2] = mWindow[g][5] - mWindow[g][4] + 1;
	  
// subtracting off 1 because C-arrays are base 0
      int start[3];
      start[0]   = mWindow[g][0] - 1;
      start[1]   = mWindow[g][2] - 1;
      start[2]   = mWindow[g][4] - 1;
          
      if (mLocationType == Image::X)
      {
	 global[0]  = 1;
	 start[0]   = 0;
      } 
      if (mLocationType == Image::Y)
      {
	 global[1]  = 1;
	 start[1]   = 0;
      }
      if (mLocationType == Image::Z)
      {
	 global[2]  = 1;
	 start[2]   = 0;
      }
      if( !plane_in_proc(m_gridPtIndex[0]) )
	 local[0]=local[1]=local[2]=0;

      int iwrite = 0;
      if( m_mpiComm_writers != MPI_COMM_NULL )
      {
	 int nproc=0, myid=0;
	 // Select group of writing processors as 
	 // subset of the processors that own the plane.
	 MPI_Comm_size(m_mpiComm_writers, &nproc);
	 MPI_Comm_rank(m_mpiComm_writers, &myid);
	 int nrwriters= mEW->getNumberOfWritersPFS();
	 if( nrwriters > nproc )
	    nrwriters = nproc;
         int q, r;
	 if( nproc == 1 || nrwriters == 1 )
	 {
	    q = 0;
            r = 0;
	 }
	 else
	 {
	    q = (nproc-1)/(nrwriters-1);
	    r = (nproc-1) % (nrwriters-1);
	 }
         for( int w=0 ; w < nrwriters ; w++ )
	    if( q*w+r == myid )
	       iwrite = 1;
      }
      m_pio[g-glow] = new Parallel_IO( iwrite, mEW->usingParallelFS(), global, local, start );
   }
}

//-----------------------------------------------------------------------
void Image::allocatePlane()
{
  ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
  bool iwrite   = plane_in_proc(m_gridPtIndex[0]);

  if (iwrite)
    {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      bool breakLoop      = (mLocationType == Image::Z);
      
// write the header
//       if (proc_write())
//       {
// 	printf("Allocating plane for image '%s'...(msg from proc #%i)\n", mFilePrefix.c_str(), m_rankWriter); 
//       }

      int nPatches = mEW->mNumberOfGrids; // both Cartesian and curvilinear 
      
      if (breakLoop) nPatches = 1;
      
// figure out the sizes...
      for (int g = 0; g < mEW->mNumberOfGrids; g++)
        {
          if (breakLoop) 
            g = m_gridPtIndex[1];
          
	  mWindow[g][0] = mEW->m_iStart[g] + mEW->m_paddingCells[0];
	  mWindow[g][1] = mEW->m_iEnd[g]   - mEW->m_paddingCells[1];
	  mWindow[g][2] = mEW->m_jStart[g] + mEW->m_paddingCells[2];
	  mWindow[g][3] = mEW->m_jEnd[g]   - mEW->m_paddingCells[3];
	  // mWindow[g][4] = mEW->m_kStart[g] + mEW->m_ghost_points;
	  // mWindow[g][5] = mEW->m_kEnd[g]   - mEW->m_ghost_points;
	  mWindow[g][4] = 1;	  
	  mWindow[g][5] = mEW->m_global_nz[g];

          if (mLocationType == Image::X)
	  {
	    mWindow[g][0] = m_gridPtIndex[g];
	    mWindow[g][1] = m_gridPtIndex[g];
	  } 
          if (mLocationType == Image::Y)
	  {
	    mWindow[g][2] = m_gridPtIndex[g];
	    mWindow[g][3] = m_gridPtIndex[g];
	  }
          if (mLocationType == Image::Z)
	  {
	    mWindow[g][4] = m_gridPtIndex[0]; // note [0] and not [g] (only one grid for z=const images)
	    mWindow[g][5] = m_gridPtIndex[0];
	  }
          
          int npts     = (mWindow[g][1] - mWindow[g][0] + 1)*(mWindow[g][3] - mWindow[g][2] + 1)*(mWindow[g][5] - mWindow[g][4] + 1);
          
          if( m_double )
	  {
	    m_doubleField[g] = new double[npts];
	  }
          else // float
	  {
	    m_floatField[g] = new float[npts];
	  }
          
          if (breakLoop) break;
        } // for g...
      
    } // end if iwrite
   define_pio();
} // end Image::allocatePlane()

//-----------------------------------------------------------------------
void Image::computeImageQuantity(std::vector<Sarray> &a_mu, int a_nComp)
{
  ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
  bool iwrite   = plane_in_proc(m_gridPtIndex[0]);

  if (iwrite)
  {
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    bool breakLoop      = (mLocationType == Image::Z);

// tmp
//     if (proc_write())
//     {
//       printf("computing image quantity for image '%s' of mode '%s'...(msg from proc #%i)\n", mFilePrefix.c_str(), 
// 	     const_cast<char*>(mMode2Suffix[mMode].c_str()), m_rankWriter); 
//     }

// copy the data for all Cartesian patches
    int g;
    for (g = 0; g < mEW->mNumberOfCartesianGrids; g++)
    {
      if (breakLoop) 
	g = m_gridPtIndex[1];
          
      int iField = 0;
              
      for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
      {
	for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	{
	  for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	  {                   
	    if( m_double )
	    {
	      m_doubleField[g][iField] = a_mu[g](a_nComp,ii,jj,kk);
	    }
	    else
	    {
	      m_floatField[g][iField] = (float) a_mu[g](a_nComp,ii,jj,kk);
	    }
		  
	    iField++;      
	  }
	}
      }
	    
          
      if (breakLoop) break;
    } // for g...
      
// curvilinear grid
    if (mEW->topographyExists())
    {
      g = mEW->mNumberOfGrids - 1;

      int iField = 0;
              
      for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
      {
	for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	{
	  for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	  {                   
	    if( m_double )
	    {
	      m_doubleField[g][iField] = a_mu[g](a_nComp,ii,jj,kk);
	    }
	    else
	    {
	      m_floatField[g][iField] = (float) a_mu[g](a_nComp,ii,jj,kk);
	    }
		  
	    iField++;      
	  }
	}
      }
    } // end if curvilinear grid

  } // end if iwrite
} // end Image::computeImageQuantity()

//-----------------------------------------------------------------------
void Image::evaluateGridImage(Sarray &a_X, Sarray &a_Y, Sarray &a_Z, int a_component)
{
  ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
  bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
  double val;
  int g;
  
  if (iwrite)
  {
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    bool breakLoop      = (mLocationType == Image::Z); // not sure how the z=const case is meaningful for grid images

    int topCartesian = mEW->mNumberOfCartesianGrids - 1;

// write the data...
    for (g = 0; g < mEW->mNumberOfGrids; g++)
    {
      if (breakLoop) 
	g = m_gridPtIndex[1];
          
      int iField = 0;
              
      for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
      {
	for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	{
	  for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	  {                   
	    if (g > topCartesian)
	    {
	      if (a_component==1)
	      {
		val = a_X(ii,jj,kk);
	      }
	      else if (a_component == 2)
	      {
		val = a_Y(ii,jj,kk);
	      }
	      else if (a_component == 3)
	      {
		val = a_Z(ii,jj,kk);
	      }
	    }
	    else
	    {
	      if (a_component==1)
	      {
		val = (ii-1)*mEW->mGridSize[g];
	      }
	      else if (a_component == 2)
	      {
		val = (jj-1)*mEW->mGridSize[g];
	      }
	      else if (a_component == 3)
	      {
		val = (kk-1)*mEW->mGridSize[g] + mEW->m_zmin[g];
	      }
	    }
	    
	    if( m_double )
	    {
	      m_doubleField[g][iField] = val;
	    }
	    else
	    {
	      m_floatField[g][iField] = (float) val;
	    }
		  
	    iField++;      
	  }
	}
      }
	    
      if (breakLoop) break;
    } // for g...

  } // end if iwrite
} // end Image::computeImageQuantity()

//-----------------------------------------------------------------------
void Image::evaluateLatLonImage(Sarray &a_X, Sarray &a_Y, Sarray &a_Z, int a_component) 
{
// a_component=1: Latitude
// a_component=2: Longitude
  ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
  bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
  double latP, lonP, xP, yP, zP, val;
  int g;
  
  if (iwrite)
  {
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

// lat, lon images are only meaningful on Image::Z planes, but we will write lat and lon on any plane

    bool breakLoop = (mLocationType == Image::Z);

    int topCartesian = mEW->mNumberOfCartesianGrids - 1;
// write the data...
    for (g = 0; g < mEW->mNumberOfGrids; g++)
    {
      if (breakLoop) 
	g = m_gridPtIndex[1];
          
      int iField = 0;
              
      for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
      {
	for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	{
	  for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	  {                   
	    if (g > topCartesian)// curvilinear grid
    	    {
	      xP = a_X(ii,jj,kk);
	      yP = a_Y(ii,jj,kk);
	      zP = a_Z(ii,jj,kk);
	    }
	    else
	    {
	      xP = (ii-1)*mEW->mGridSize[g];
	      yP = (jj-1)*mEW->mGridSize[g];
	      zP = (kk-1)*mEW->mGridSize[g] + mEW->m_zmin[g];
	    }
// evaluate the mapping to (lon,lat)	    
	    mEW->computeGeographicCoord(xP, yP, lonP, latP);
	    if (a_component==1)
	    {
	      val = latP;
	    }
	    else if (a_component == 2)
	    {
	      val = lonP;
	    }
	    
	    if( m_double )
	    {
	      m_doubleField[g][iField] = val;
	    }
	    else
	    {
	      m_floatField[g][iField] = (float) val;
	    }
		  
	    iField++;      
	  } // end for ii	  
	} // end for jj	
      } // end for kk
	    
      if (breakLoop) break;
    } // for g...

  } // end if iwrite
} // end Image::evaluateLatLonImage()

//-----------------------------------------------------------------------
// void Image::computeImageError(std::vector<Sarray> &mu, int nComp)
// {
//   ASSERT(m_isDefinedMPIWriters);
//   if (!mEW->m_forcing->knows_exact())
//   {
//     if (proc_write())
//       cout << "Image::computeImageError can only be called when the exact solution is known" << endl;
//     MPI_Abort(MPI_COMM_WORLD, 1);
//   }


// // plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
//   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);

//   if (iwrite)
//     {
//       int myRank;
//       MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

//       bool breakLoop      = (mLocationType == Image::Z);

//       double err, x, y, z, t=mEW->getCurrentTime();
//       double uExact[3];
      
// // tmp
// //       if (proc_write())
// //       {
// // 	printf("computing image error for image '%s' of mode '%s' time=%e...(msg from proc #%i)\n", mFilePrefix.c_str(), 
// // 	       const_cast<char*>(mMode2Suffix[mMode].c_str()), t, m_rankWriter); 
// //       }

// // tmp
// //      double maxerr=0.;
      
// // write the data...
//       int g;
      
//       for (g = 0; g < mEW->mNumberOfCartesianGrids; g++)
//       {
// 	if (breakLoop) 
// 	  g = m_gridPtIndex[1];
          
// 	int iField = 0;
              
// 	for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
// 	{
// 	  for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
// 	  {
// 	    for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
// 	    {                   
// 	      x = (ii - 1)*mEW->mGridSize[g];
// 	      y = (jj - 1)*mEW->mGridSize[g];
// 	      z = (kk - 1)*mEW->mGridSize[g] + mEW->m_zmin[g];
// 	      mEW->m_forcing->get_exact(x, y, z, t, uExact, mEW->mGridSize[g]);
		  
// 	      err = ( mu[g](nComp,ii,jj,kk) - uExact[nComp-1] );
// //	      maxerr = (fabs(err) > maxerr ? fabs(err) : maxerr);

// 	      if( m_double )
// 	      {
// 		m_doubleField[g][iField] = err;
// 	      }
// 	      else
// 	      {
// 		m_floatField[g][iField] = (float) err;
// 	      }
		  
// 	      iField++;      
// 	    }
// 	  }
// 	}
	    
          
// 	if (breakLoop) break;
//       } // for g...

// // curvilinear grid
//       if (mEW->topographyExists())
//       {
// 	g = mEW->mNumberOfGrids - 1;
	
// 	int iField = 0;
              
// 	for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
// 	{
// 	  for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
// 	  {
// 	    for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
// 	    {                   
// 	      x = mEW->mX(ii,jj,kk);
// 	      y = mEW->mY(ii,jj,kk);
// 	      z = mEW->mZ(ii,jj,kk);

// 	      mEW->m_forcing->get_exact(x, y, z, t, uExact, mEW->mGridSize[g]);
		  
// 	      err = ( mu[g](nComp,ii,jj,kk) - uExact[nComp-1] );
// //	      maxerr = (fabs(err) > maxerr ? fabs(err) : maxerr);

// 	      if( m_double )
// 	      {
// 		m_doubleField[g][iField] = err;
// 	      }
// 	      else
// 	      {
// 		m_floatField[g][iField] = (float) err;
// 	      }
		  
// 	      iField++;      
// 	    }
// 	  }
// 	}
//       } // end if topographyExists
      
	
// // tmp
// //      printf("Image::maxerr=%e\n", maxerr);
      
//     } // end if iwrite
// } // end Image::computeImageError()

//-----------------------------------------------------------------------
void Image::computeImageP(std::vector<Sarray> &mu, std::vector<Sarray> &lambda, std::vector<Sarray> &rho)
{
  ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
  bool iwrite   = plane_in_proc(m_gridPtIndex[0]);

  if (iwrite)
    {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      bool breakLoop      = (mLocationType == Image::Z);

      double vp;
      
// tmp
//       if (proc_write())
//       {
// 	printf("computing P-velocity for image '%s' of mode '%s'...(msg from proc #%i)\n", mFilePrefix.c_str(), 
// 	       const_cast<char*>(mMode2Suffix[mMode].c_str()), m_rankWriter); 
//       }

// write the data...
      for (int g = 0; g < mEW->mNumberOfGrids; g++)
      {
	if (breakLoop) 
	  g = m_gridPtIndex[1];
          
	int iField = 0;
              
	for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
	{
	  for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	  {
	    for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	    {                   
	      vp = sqrt( (2*mu[g](ii,jj,kk) + lambda[g](ii,jj,kk))/rho[g](ii,jj,kk) );

	      if( m_double )
	      {
		m_doubleField[g][iField] = vp;
	      }
	      else
	      {
		m_floatField[g][iField] = (float) vp;
	      }
		  
	      iField++;      
	    }
	  }
	}
	    
          
	if (breakLoop) break;
      } // for g...

// tmp
//      printf("Image::maxerr=%e\n", maxerr);
      
    } // end if iwrite
} // end Image::computeImageP()

//-----------------------------------------------------------------------
void Image::computeImageS(std::vector<Sarray> &mu, std::vector<Sarray> &lambda, std::vector<Sarray> &rho)
{
  ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
  bool iwrite   = plane_in_proc(m_gridPtIndex[0]);

  if (iwrite)
    {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      bool breakLoop      = (mLocationType == Image::Z);

      double vs;
      
// tmp
//       if (proc_write())
//       {
// 	printf("computing S-velocity for image '%s' of mode '%s'...(msg from proc #%i)\n", mFilePrefix.c_str(), 
// 	       const_cast<char*>(mMode2Suffix[mMode].c_str()), m_rankWriter); 
//       }

// write the data...
      for (int g = 0; g < mEW->mNumberOfGrids; g++)
      {
	if (breakLoop) 
	  g = m_gridPtIndex[1];
          
	int iField = 0;
              
	for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
	{
	  for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	  {
	    for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	    {                   
	      vs = sqrt( mu[g](ii,jj,kk)/rho[g](ii,jj,kk) );

	      if( m_double )
	      {
		m_doubleField[g][iField] = vs;
	      }
	      else
	      {
		m_floatField[g][iField] = (float) vs;
	      }
		  
	      iField++;      
	    }
	  }
	}
	    
          
	if (breakLoop) break;
      } // for g...

    } // end if iwrite
} // end Image::computeImageS()

//-----------------------------------------------------------------------
void Image::copy2DArrayToImage(Sarray &u2)
{
  ASSERT(m_isDefinedMPIWriters);

   REQUIRE2( u2.m_nc == 1, "Image::copy2DArrayToImage, only implemented for one-component arrays" );

   int ie = u2.m_ie, ib=u2.m_ib, je=u2.m_je, jb=u2.m_jb;
   REQUIRE2( ib == mEW->m_iStart[mEW->mNumberOfGrids-1] && ie == mEW->m_iEnd[mEW->mNumberOfGrids-1] &&
             jb == mEW->m_jStart[mEW->mNumberOfGrids-1] && je == mEW->m_jEnd[mEW->mNumberOfGrids-1] , 
             "Communicate array 2d: Can only use it on the finest grid, grid sizes don't match");

// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
  bool iwrite   = plane_in_proc(m_gridPtIndex[0]);

  if (iwrite)
    {
      ASSERT2(mLocationType == Image::Z, "Image::copy2DArrayToImage only works for z=const");

      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      bool breakLoop      = (mLocationType == Image::Z);

      double val;
      
// tmp
//       if (proc_write())
//       {
// 	printf("copying 2D array to image '%s' of mode '%s'...(msg from proc #%i)\n", mFilePrefix.c_str(), 
// 	       const_cast<char*>(mMode2Suffix[mMode].c_str()), m_rankWriter); 
//       }

// write the data...
      int g = m_gridPtIndex[1];
          
      int iField = 0;
              
      for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
      {
	for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	{                   
	  val = u2(ii,jj,1);

	  if( m_double )
	  {
	    m_doubleField[g][iField] = val;
	  }
	  else
	  {
	    m_floatField[g][iField] = (float) val;
	  }
		  
	  iField++;      
	}
      }

    } // end if iwrite
} // end Image::computeImageS()


//-----------------------------------------------------------------------
void Image::writeImagePlane_2(int cycle, std::string &path )
{
  ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
  bool ihavearray = plane_in_proc(m_gridPtIndex[0]);

  int glow = 0, ghigh = mEW->mNumberOfGrids;
  if( mLocationType == Image::Z )
  {
     glow = m_gridPtIndex[1];
     ghigh = glow+1;
  }

  bool iwrite = false;
  for( int g=glow ; g < ghigh ; g++ )
     iwrite = iwrite || m_pio[g-glow]->i_write();
  
  // Header: [precision(int) npatches(int) h_1 sizes_1 h_2 sizes_2 ... h_ng sizes_ng]
  // offset initialized to header size:
  off_t offset = 2*sizeof(int) + (ghigh-glow)*(sizeof(double)+4*sizeof(int));
  int fid=-1;
  stringstream s, fileSuffix;
  int prec, nPatches, globalPlaneSize[4];

  if( iwrite )
  {
     mCycle = cycle; // mCycle must be set before: compute_file_suffix
     compute_file_suffix( fileSuffix );
     if( path != "." )
        s << path;
     s << fileSuffix.str();
//     s << fileSuffix.str() << "new";
  }

  //  if( m_pio[0]->i_write() )
  //  {
     // Need to synchronize writing of header if file system is not parallel:
     //     m_pio[0]->begin_sequential( );
     //     if (proc_write())
     if( m_pio[0]->proc_zero() )
     {
	fid = open( const_cast<char*>(s.str().c_str()), O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
	if (fid == -1 )
	{
	   //	char buffer[MPI_MAX_ERROR_STRING];
	   //	int resultlen;
	   //	MPI_Error_string( ierr, buffer, &resultlen);
	   VERIFY2(0, "write_images:: Error opening: " << s );
	}  
	cout << "writing image plane on file " << s.str() << " (msg from proc # " << m_rankWriter << ")" << endl;

        prec = m_double ? 8 : 4;
        int ret=write(fid,&prec,sizeof(int));

        nPatches = mLocationType == Image::Z ? 1 : mEW->mNumberOfGrids;
        ret = write(fid,&nPatches,sizeof(int));
	for(int g = glow; g < ghigh ;g++)
        {
	   //	   h = mEW->mGridSize[g];
           ret = write(fid,&mEW->mGridSize[g],sizeof(double));

	   //	   int globalPlaneSize[4];
	   // should hold the global number of interior points
	   if(mLocationType == Image::X)
	   {
	      globalPlaneSize[0] = 1;
	      globalPlaneSize[1] = mEW->m_global_ny[g];
	      globalPlaneSize[2] = 1;
//	      globalPlaneSize[3] = mEW->m_kEnd[g] - mEW->m_kStart[g] - 2*mEW->m_ghost_points + 1;
	      globalPlaneSize[3] = mEW->m_global_nz[g];
	   }
	   if(mLocationType == Image::Y)
	   {
	      globalPlaneSize[0] = 1;
	      globalPlaneSize[1] = mEW->m_global_nx[g];
	      globalPlaneSize[2] = 1;
//	      globalPlaneSize[3] = mEW->m_kEnd[g] - mEW->m_kStart[g] - 2*mEW->m_ghost_points + 1;
	      globalPlaneSize[3] = mEW->m_global_nz[g];
	   }
	   if(mLocationType == Image::Z)
	   {
	      globalPlaneSize[0] = 1;
	      globalPlaneSize[1] = mEW->m_global_nx[g];
	      globalPlaneSize[2] = 1;
	      globalPlaneSize[3] = mEW->m_global_ny[g];
	   }
// tmp
	   // printf("Image::writeImagePlane_2: globalPlaneSize=(%i %i %i %i)\n", globalPlaneSize[0], 
	   // 	  globalPlaneSize[1], globalPlaneSize[2], globalPlaneSize[3]);
	   
           ret = write(fid,globalPlaneSize,4*sizeof(int));
	}
	//	fsync(fid);
     }
     //     m_pio[0]->end_sequential();
     //  }

      // write the data...
  if( ihavearray )
  {
     for( int g = glow ; g < ghigh; g++)
     {
	int globalSizes[3] = {mEW->m_global_nx[g],
			      mEW->m_global_ny[g],
			      mEW->m_global_nz[g]} ;
//			      mEW->m_kEnd[g] - mEW->m_kStart[g] - 2*mEW->m_ghost_points+1} ;

	if(mLocationType == Image::X)
	   globalSizes[0]    = 1;
	if (mLocationType == Image::Y)
	   globalSizes[1]    = 1;
	if (mLocationType == Image::Z)
	   globalSizes[2]    = 1;

        if( !mEW->usingParallelFS() || g == glow )
	   MPI_Barrier( m_mpiComm_writers );

	if( g==glow && iwrite && !m_pio[0]->proc_zero() )
	{
	   fid = open( const_cast<char*>(s.str().c_str()), O_WRONLY );
	   if (fid == -1 )
	   {
	      VERIFY2(0, "write_images:: Error opening: " << s );
	   }  
	}

	if( m_double )
	{
	  char dblStr[]="double";	  
	  m_pio[g-glow]->write_array( &fid, 1, m_doubleField[g], offset, dblStr );
	  offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]*sizeof(double));
	}
	else
	{
	  char fltStr[]="float";
	  m_pio[g-glow]->write_array( &fid, 1, m_floatField[g], offset, fltStr );
	  offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]*sizeof(float));
	}
     }
  }
  if( iwrite )
     close(fid);
}


//-----------------------------------------------------------------------
void Image::compute_file_suffix( stringstream& fileSuffix )
{
   fileSuffix << mFilePrefix << ".cycle=";
   int temp = static_cast<int>(pow(10.0, mPreceedZeros - 1));
   int testcycle = mCycle;
   if (mCycle == 0)
      testcycle=1;
   while (testcycle < temp)
   {
      fileSuffix << "0";
      temp /= 10;
   }
   fileSuffix << mCycle << "." << mOrientationString[mLocationType] << "=";
   fileSuffix << mCoordValue;
   fileSuffix << "." << fieldSuffix(mMode);
}

//C.B> with N-S and E-W max: max(|u_{N-S}|,|u_{E-W}|))

// void Image::update_maxes_hVelMax()
// {
//   static bool firstHVM = true;
//   bool breakLoop = (mLocationType == Image::Z)    ;
//   double velNS, velEW;
  
//   for( int g = 0; g < mEW->mNumberOfGrids; g++ )
//   {
//     double di  = 1/(2*mEW->mDt)  ;
//     double di2 = 1/(mEW->mDt*mEW->mDt);
      
//     if (breakLoop) 
//       g = m_gridPtIndex[1];

//     int iField = 0;
//     double Ux, Uy, nrm, xP, yP, latitude, latOrigin=mEW->getLatOrigin();
//     double deg2rad = M_PI/180.0, cphi, sphi, calpha, salpha, thxnrm, thynrm;
//     double az = deg2rad*mEW->getGridAzimuth(); // Note that mGeoAz is in degrees
//     double metersPerDeg = mEW->getMetersPerDegree();

//     calpha = cos(az);
//     salpha = sin(az);

//     for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//       for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
// 	for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
// 	{
// // first get velocities in the (x,y) directions
// 	  Ux = (mEW->mUp[g](1,i,j,k) - mEW->mUm[g](1,i,j,k))*di;
// 	  Uy = (mEW->mUp[g](2,i,j,k) - mEW->mUm[g](2,i,j,k))*di;
	  
// // note: the angle between (x,y) and East is constant throughout the mesh, but the angle to North 
// // varies throughout the grid. Here we use the chain rule on the mapping between (x,y) and (lon,lat)
// // see EW:computeGeographicCoord for the mapping
// 	  xP = (i-1)*mEW->mGridSize[g];
// 	  yP = (j-1)*mEW->mGridSize[g];

// 	  latitude = latOrigin + (xP*calpha - yP*salpha)/metersPerDeg;
	  
// 	  cphi = cos(latitude*deg2rad);
// 	  sphi = sin(latitude*deg2rad);

// 	  thxnrm = salpha + (xP*salpha+yP*calpha)/cphi/metersPerDeg * (M_PI/180.0) * sphi * calpha;
// 	  thynrm = calpha - (xP*salpha+yP*calpha)/cphi/metersPerDeg * (M_PI/180.0) * sphi * salpha;

// 	  nrm = sqrt( thxnrm*thxnrm + thynrm*thynrm );
// 	  thxnrm /= nrm;
// 	  thynrm /= nrm;

// 	  velNS = fabs( thynrm*Ux  - thxnrm*Uy );
// 	  velEW = fabs( salpha*Ux  + calpha*Uy );

// //	  velNS = fabs( Ux*cos(az) - Uy*sin(az) );
// //	  velEW = fabs( Ux*sin(az) + Uy*cos(az) );

// 	  if (m_double)
// 	  {
// 	    if( firstHVM ||  m_doubleField[g][iField] < velNS)
// 	      m_doubleField[g][iField] = velNS;
                  
// 	    if( m_doubleField[g][iField] < velEW )
// 	      m_doubleField[g][iField] = velEW;
// 	  }
// 	  else
// 	  {
// 	    if( firstHVM ||  m_floatField[g][iField] < (float) velNS)
// 	      m_floatField[g][iField] = (float) velNS;
// 	    if(m_floatField[g][iField] < (float) velEW)
// 	      m_floatField[g][iField] = (float) velEW;
// 	  }
// 	  iField++;  
// 	}

//     if (breakLoop) break;
//   }
//   if( firstHVM )
//     firstHVM = false;
// }

// void Image::update_maxes_vVelMax()
// {
//   static bool firstVVM = true;
//   bool breakLoop = (mLocationType == Image::Z)    ;
  
//   for( int g = 0; g < mEW->mNumberOfGrids; g++ )
//   {
//     int iField = 0;

//     if (breakLoop) 
//       g = m_gridPtIndex[1];
      
//     double di  = 1/(2*mEW->mDt)  ;
//     double di2 = 1/(mEW->mDt*mEW->mDt);
      
//     for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//       for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
// 	for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
// 	{
// 	  double vel = fabs( mEW->mUp[g](3,i,j,k)-mEW->mUm[g](3,i,j,k) )*di;
// 	  if (m_double)
// 	  {
// 	    if( firstVVM ||  m_doubleField[g][iField] < vel )
// 	      m_doubleField[g][iField] = vel;
// 	  }
// 	  else
// 	  {
// 	    if( firstVVM ||  m_floatField[g][iField] < (float) vel )
// 	      m_floatField[g][iField] = (float) vel;

// 	  }
// 	  iField++; 
// 	}
//     if (breakLoop) break;
//   }
     
//   if( firstVVM )
//     firstVVM = false;
// }

// //-----------------------------------------------------------------------
// void Image::computeDivergence()
// {
//   ASSERT(m_isDefinedMPIWriters);

// // plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
//   bool iwrite    = plane_in_proc(m_gridPtIndex[0]);
//   bool breakLoop = (mLocationType == Image::Z)    ;
//   bool done      = false                          ;

//   if (iwrite)
//     {
//       for (int g = 0; g < mEW->mNumberOfCartesianGrids; g++)
//         {
//           if (breakLoop) 
//             {
//               g = m_gridPtIndex[1];
              
//               if (g >= mEW->mNumberOfCartesianGrids)
//                 break;
//             }

//           double factor = 1.0/(2*mEW->mGridSize[g]);
//           int iField = 0;
          
//           for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//             for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//               for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//                 {
//                   double val =0.;
                  
//                   if (i == 1)
//                     {
//                       val  = (-mEW->mUp[g](1,i+2, j, k)+4.*mEW->mUp[g](1,i+1, j, k)-3.*mEW->mUp[g](1,i, j, k));
//                     }
//                   else if (i == mEW->m_global_nx[g])
//                     {
//                       val  = (3.*mEW->mUp[g](1,i, j, k)-4.*mEW->mUp[g](1,i-1, j, k)+mEW->mUp[g](1,i-2, j, k));
//                     }
//                   else
//                     {
//                       val  =  (mEW->mUp[g](1,i+1, j, k) - mEW->mUp[g](1,i-1, j, k));
//                     }
                  
//                   if (j == 1)
//                     {
//                       val  += (-mEW->mUp[g](2,i, j+2, k)+4.*mEW->mUp[g](2,i, j+1, k)-3.*mEW->mUp[g](2,i,j,k));
//                     }
//                   else if (j == mEW->m_global_ny[g])
//                     {
//                       val  += (3.*mEW->mUp[g](2,i, j, k)-4.*mEW->mUp[g](2,i,j-1,k)+mEW->mUp[g](2,i,j-2,k));
//                     }
//                   else
//                     {
//                       val  += (mEW->mUp[g](2,i, j+1, k) - mEW->mUp[g](2,i,j-1,k));
//                     }
                  
//                   if (k == 1 && g==mEW->mNumberOfGrids-1)
//                     {
//                       val  += (-mEW->mUp[g](3,i,j,k+2)+4.*mEW->mUp[g](3,i,j,k+1)-3.*mEW->mUp[g](3,i,j,k));
//                     }
// //                  else if (k == mEW->m_kEnd[g]-mEW->m_ghost_points && g == 0)
//                   else if (k == mEW->m_global_nz[g] && g == 0)
//                     {
//                       val  += (3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i,j,k-1)+mEW->mUp[g](3,i,j,k-2));
//                     }
//                   else
//                     {
//                       val  += (mEW->mUp[g](3,i,j,k+1) - mEW->mUp[g](3,i,j,k-1));
//                     }
                  
//                   val *= factor;
                  
//                   if (m_double)
//                     {
//                       m_doubleField[g][iField] =val;
//                     }
//                   else
//                     {
//                       m_floatField[g][iField] = (float)val;
//                     }
//                   iField++; 
//                 }
//           if (breakLoop) 
//             {
//               done = true;
//               break;
//             }
//         }
      
//       // curvilinear grid
//       if (mEW->topographyExists() && !done)
//         {
//           int g = mEW->mNumberOfGrids - 1;
//           double factor = 1.0/(2.); //C.B: metric terms are undivided differences. So, we suppress the dx in factor.
          
//           int iField = 0;
          
//           for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//             for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//               for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//                 {
//                   double val = 0.;
                  
//                   if (i == 1)
//                     {
//                       val  = (mEW->mQ(1,i,j,k)*(-mEW->mUp[g](1,i+2,j,k)+4.*mEW->mUp[g](1,i+1,j,k)-3.*mEW->mUp[g](1,i,j,k))
//                              +mEW->mQ(2,i,j,k)*(-mEW->mUp[g](2,i+2,j,k)+4.*mEW->mUp[g](2,i+1,j,k)-3.*mEW->mUp[g](2,i,j,k))
//                              +mEW->mQ(3,i,j,k)*(-mEW->mUp[g](3,i+2,j,k)+4.*mEW->mUp[g](3,i+1,j,k)-3.*mEW->mUp[g](3,i,j,k)));
//                     }
//                   else if (i == mEW->m_global_nx[g])
                    
//                     {
//                       val  = (mEW->mQ(1,i,j,k)*(3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i-1,j,k)+mEW->mUp[g](1,i-2,j,k))
//                              +mEW->mQ(2,i,j,k)*(3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i-1,j,k)+mEW->mUp[g](2,i-2,j,k))
//                              +mEW->mQ(3,i,j,k)*(3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i-1,j,k)+mEW->mUp[g](3,i-2,j,k)));
//                     }
//                   else
//                     {
//                       val  = (mEW->mQ(1,i,j,k)*(mEW->mUp[g](1,i+1,j,k) - mEW->mUp[g](1,i-1,j,k))
//                              +mEW->mQ(2,i,j,k)*(mEW->mUp[g](2,i+1,j,k) - mEW->mUp[g](2,i-1,j,k))
//                              +mEW->mQ(3,i,j,k)*(mEW->mUp[g](3,i+1,j,k) - mEW->mUp[g](3,i-1,j,k)));
//                     }
                  
//                   if (j == 1)
//                     {
//                       val  += (mEW->mR(1,i,j,k)*(-mEW->mUp[g](1,i,j+2,k)+4.*mEW->mUp[g](1,i,j+1,k)-3.*mEW->mUp[g](1,i,j,k))
//                               +mEW->mR(2,i,j,k)*(-mEW->mUp[g](2,i,j+2,k)+4.*mEW->mUp[g](2,i,j+1,k)-3.*mEW->mUp[g](2,i,j,k))
//                               +mEW->mR(3,i,j,k)*(-mEW->mUp[g](3,i,j+2,k)+4.*mEW->mUp[g](3,i,j+1,k)-3.*mEW->mUp[g](3,i,j,k)));
//                     }
//                   else if (j == mEW->m_global_ny[g])
//                     {
//                       val  += (mEW->mR(1,i,j,k)*(3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i,j-1,k)+mEW->mUp[g](1,i,j-2,k))
//                               +mEW->mR(2,i,j,k)*(3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i,j-1,k)+mEW->mUp[g](2,i,j-2,k))
//                               +mEW->mR(3,i,j,k)*(3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i,j-1,k)+mEW->mUp[g](3,i,j-2,k)));
//                     }
//                   else
//                     {
//                       val  += (mEW->mR(1,i,j,k)*(mEW->mUp[g](1,i,j+1,k) - mEW->mUp[g](1,i,j-1,k))
//                               +mEW->mR(2,i,j,k)*(mEW->mUp[g](2,i,j+1,k) - mEW->mUp[g](2,i,j-1,k))
//                               +mEW->mR(3,i,j,k)*(mEW->mUp[g](3,i,j+1,k) - mEW->mUp[g](3,i,j-1,k)));
//                     }
                  
//                   if (k == 1)
//                     {
//                       val  += (mEW->mS(1,i,j,k)*(-mEW->mUp[g](1,i,j,k+2)+4.*mEW->mUp[g](1,i,j,k+1)-3.*mEW->mUp[g](1,i,j,k))
//                               +mEW->mS(2,i,j,k)*(-mEW->mUp[g](2,i,j,k+2)+4.*mEW->mUp[g](2,i,j,k+1)-3.*mEW->mUp[g](2,i,j,k))
//                               +mEW->mS(3,i,j,k)*(-mEW->mUp[g](3,i,j,k+2)+4.*mEW->mUp[g](3,i,j,k+1)-3.*mEW->mUp[g](3,i,j,k)));
//                     }
//                   else // no k=N because we are in the curvilinear grid
//                     {
//                       val  += (mEW->mS(1,i,j,k)*(mEW->mUp[g](1,i,j,k+1) - mEW->mUp[g](1,i,j,k-1))
//                               +mEW->mS(2,i,j,k)*(mEW->mUp[g](2,i,j,k+1) - mEW->mUp[g](2,i,j,k-1))
//                               +mEW->mS(3,i,j,k)*(mEW->mUp[g](3,i,j,k+1) - mEW->mUp[g](3,i,j,k-1)));
//                     }
                  
//                   val *= factor;  

//                   if( m_double )
//                     {
//                       m_doubleField[g][iField] = val;
//                     }
//                   else
//                     {
//                       m_floatField[g][iField] = (float)val;  
//                     }
//                   iField++;  
//                 }    
//         } // end if curvilinear grid
//     }
// }

// //-----------------------------------------------------------------------
// void Image::computeCurlMagnitude()
// {
//   ASSERT(m_isDefinedMPIWriters);

// // plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
//   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
//   bool breakLoop = (mLocationType == Image::Z)    ;
//   bool done      = false                          ;

//   if (iwrite)
//     {
//       for (int g = 0; g < mEW->mNumberOfCartesianGrids; g++)
//         {
//           if (breakLoop) 
//             {
//               g = m_gridPtIndex[1];
              
//               if (g >= mEW->mNumberOfCartesianGrids)
//                 break;
//             }

//           double factor = 1.0/(2*mEW->mGridSize[g]);
//           int iField = 0;
          
//           for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//             for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//               for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//                 {
//                   double duydx = 0.;
//                   double duzdx = 0.;
//                   double duxdy = 0.;
//                   double duzdy = 0.;
//                   double duxdz = 0.;
//                   double duydz = 0.;
                  
//                   if (i == 1)
//                     {
//                       duydx = (-mEW->mUp[g](2,i+2,j,k)+4.*mEW->mUp[g](2,i+1,j,k)-3.*mEW->mUp[g](2,i,j,k));
//                       duzdx = (-mEW->mUp[g](3,i+2,j,k)+4.*mEW->mUp[g](3,i+1,j,k)-3.*mEW->mUp[g](3,i,j,k));
//                     }
//                   else if (i == mEW->m_global_nx[g])
//                     {
//                       duydx = (3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i-1,j,k)+mEW->mUp[g](2,i-2,j,k));
//                       duzdx = (3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i-1,j,k)+mEW->mUp[g](3,i-2,j,k));
//                     }
//                   else
//                     {
//                       duydx = (mEW->mUp[g](2,i+1,j,k) - mEW->mUp[g](2,i-1,j,k));
//                       duzdx = (mEW->mUp[g](3,i+1,j,k) - mEW->mUp[g](3,i-1,j,k));
//                     }
                  
//                   if (j == 1)
//                     {   
//                       duxdy = (-mEW->mUp[g](1,i,j+2,k)+4.*mEW->mUp[g](1,i,j+1,k)-3.*mEW->mUp[g](1,i,j,k));
//                       duzdy = (-mEW->mUp[g](3,i,j+2,k)+4.*mEW->mUp[g](3,i,j+1,k)-3.*mEW->mUp[g](3,i,j,k));
//                     }
//                   else if (j == mEW->m_global_ny[g])
//                     { 
//                       duxdy = (3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i,j-1,k)+mEW->mUp[g](1,i,j-2,k));
//                       duzdy = (3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i,j-1,k)+mEW->mUp[g](3,i,j-2,k));
//                     }
//                   else
//                     {
//                       duxdy = (mEW->mUp[g](1,i,j+1,k) - mEW->mUp[g](1,i,j-1,k));
//                       duzdy = (mEW->mUp[g](3,i,j+1,k) - mEW->mUp[g](3,i,j-1,k));
//                     }
                  
//                   if (k == 1 && g==mEW->mNumberOfGrids-1)
//                     {
//                       duxdz = (-mEW->mUp[g](1,i,j,k+2)+4.*mEW->mUp[g](1,i,j,k+1)-3.*mEW->mUp[g](1,i,j,k));
//                       duydz = (-mEW->mUp[g](2,i,j,k+2)+4.*mEW->mUp[g](2,i,j,k+1)-3.*mEW->mUp[g](2,i,j,k));
//                     }
// //                  else if (k == mEW->m_kEnd[g]-mEW->m_ghost_points && g == 0)
//                   else if (k == mEW->m_global_nz[g] && g == 0)
//                     {
//                       duxdz = (3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i,j,k-1)+mEW->mUp[g](1,i,j,k-2));
//                       duydz = (3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i,j,k-1)+mEW->mUp[g](2,i,j,k-2));
//                     }
//                   else // no k=N because we are in the curvilinear grid
//                     {
//                       duxdz = (mEW->mUp[g](1,i,j,k+1) - mEW->mUp[g](1,i,j,k-1));
//                       duydz = (mEW->mUp[g](2,i,j,k+1) - mEW->mUp[g](2,i,j,k-1));
//                     }
                                    
//                   if (m_double)
//                     {
//                       m_doubleField[g][iField] = (duzdy-duydz)*(duzdy-duydz)+(duxdz-duzdx)*(duxdz-duzdx)+(duydx-duxdy)*(duydx-duxdy);
                      
//                       m_doubleField[g][iField] = factor*sqrt(m_doubleField[g][iField]);
//                     }
//                   else
//                     {
//                       m_floatField[g][iField] = (duzdy-duydz)*(duzdy-duydz)+(duxdz-duzdx)*(duxdz-duzdx)+(duydx-duxdy)*(duydx-duxdy);
                      
//                       m_floatField[g][iField] = factor*sqrt(m_floatField[g][iField]);
//                     }
//                   iField++; 
//                 }
                  
//           if (breakLoop) 
//             {
//               done = true;
//               break;
//             }
//         }
      
//       if (mEW->topographyExists() && !done)
//         {
//           int g = mEW->mNumberOfGrids - 1;
//           int iField = 0;
//           double factor = 1.0/(2);
          
//           for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//             for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//               for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//                 {
//                   double duxdq = 0.;
//                   double duydq = 0.;
//                   double duzdq = 0.;
//                   double duxdr = 0.;
//                   double duydr = 0.;
//                   double duzdr = 0.;
//                   double duxds = 0.;
//                   double duyds = 0.;
//                   double duzds = 0.;
                  
//                   if (i == 1)
//                     {      
//                       duxdq = (-mEW->mUp[g](1,i+2,j,k)+4.*mEW->mUp[g](1,i+1,j,k)-3.*mEW->mUp[g](1,i,j,k));
//                       duydq = (-mEW->mUp[g](2,i+2,j,k)+4.*mEW->mUp[g](2,i+1,j,k)-3.*mEW->mUp[g](2,i,j,k));
//                       duzdq = (-mEW->mUp[g](3,i+2,j,k)+4.*mEW->mUp[g](3,i+1,j,k)-3.*mEW->mUp[g](3,i,j,k));
//                     }
//                   else if (i == mEW->m_global_nx[g])
//                     { 
//                       duxdq = (3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i-1,j,k)+mEW->mUp[g](1,i-2,j,k));
//                       duydq = (3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i-1,j,k)+mEW->mUp[g](2,i-2,j,k));
//                       duzdq = (3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i-1,j,k)+mEW->mUp[g](3,i-2,j,k));
//                     }
//                   else
//                     { 
//                       duxdq = (mEW->mUp[g](1,i+1,j,k) - mEW->mUp[g](1,i-1,j,k));
//                       duydq = (mEW->mUp[g](2,i+1,j,k) - mEW->mUp[g](2,i-1,j,k));
//                       duzdq = (mEW->mUp[g](3,i+1,j,k) - mEW->mUp[g](3,i-1,j,k));
//                     }
                  
//                   if (j == 1)
//                     {   
//                       duxdr = (-mEW->mUp[g](1,i,j+2,k)+4.*mEW->mUp[g](1,i,j+1,k)-3.*mEW->mUp[g](1,i,j,k));
//                       duydr = (-mEW->mUp[g](2,i,j+2,k)+4.*mEW->mUp[g](2,i,j+1,k)-3.*mEW->mUp[g](2,i,j,k));
//                       duzdr = (-mEW->mUp[g](3,i,j+2,k)+4.*mEW->mUp[g](3,i,j+1,k)-3.*mEW->mUp[g](3,i,j,k));
//                     }
//                   else if (j == mEW->m_global_ny[g])
//                     { 
//                       duxdr = (3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i,j-1,k)+mEW->mUp[g](1,i,j-2,k));
//                       duydr = (3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i,j-1,k)+mEW->mUp[g](2,i,j-2,k));
//                       duzdr = (3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i,j-1,k)+mEW->mUp[g](3,i,j-2,k));
//                     }
//                   else
//                     {
//                       duxdr = (mEW->mUp[g](1,i,j+1,k) - mEW->mUp[g](1,i,j-1,k));
//                       duydr = (mEW->mUp[g](2,i,j+1,k) - mEW->mUp[g](2,i,j-1,k));
//                       duzdr = (mEW->mUp[g](3,i,j+1,k) - mEW->mUp[g](3,i,j-1,k));
//                     }
                  
//                   /*S.D*/
// //                   if ((j==51&&k==1&&mEW->mGridSize[g]==0.02) || (j==101&&k==1&&mEW->mGridSize[g]==0.01))
// //                     {
// //                       printf("duzdr %f %f %f %f \n",3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i,j-1,k)+mEW->mUp[g](3,i,j-2,k),mEW->mUp[g](3,i,j,k),mEW->mUp[g](3,i,j-1,k),mEW->mUp[g](3,i,j-2,k));
// //                     }         
                  
//                   /*E.D*/

//                   if (k == 1)
//                     {
//                       duxds = (-mEW->mUp[g](1,i,j,k+2)+4.*mEW->mUp[g](1,i,j,k+1)-3.*mEW->mUp[g](1,i,j,k));
//                       duyds = (-mEW->mUp[g](2,i,j,k+2)+4.*mEW->mUp[g](2,i,j,k+1)-3.*mEW->mUp[g](2,i,j,k));
//                       duzds = (-mEW->mUp[g](3,i,j,k+2)+4.*mEW->mUp[g](3,i,j,k+1)-3.*mEW->mUp[g](3,i,j,k));
//                     }
//                   else // no k=N because we are in the curvilinear grid
//                     {
//                       duxds = (mEW->mUp[g](1,i,j,k+1) - mEW->mUp[g](1,i,j,k-1));
//                       duyds = (mEW->mUp[g](2,i,j,k+1) - mEW->mUp[g](2,i,j,k-1));
//                       duzds = (mEW->mUp[g](3,i,j,k+1) - mEW->mUp[g](3,i,j,k-1));
//                     }                  

//                   /*S.D*/
// //                   if ((j==51&&k==1&&mEW->mGridSize[g]==0.02)|| (j==101&&k==1&&mEW->mGridSize[g]==0.01))
// //                     {
// //                       printf("duzdq  duzdr duzds %f %f %f \n",duzdq,duzdr,duzds);
// //                       printf("duydq  duydr duyds %f %f %f \n",duydq,duydr,duyds);
// //                     }          
//                   /*E.D*/

//                   double duzdy = mEW->mQ(2,i,j,k)*duzdq+mEW->mR(2,i,j,k)*duzdr+mEW->mS(2,i,j,k)*duzds;
                  
//                   double duydz = mEW->mQ(3,i,j,k)*duydq+mEW->mR(3,i,j,k)*duydr+mEW->mS(3,i,j,k)*duyds;
                  
//                   double duxdz = mEW->mQ(3,i,j,k)*duxdq+mEW->mR(3,i,j,k)*duxdr+mEW->mS(3,i,j,k)*duxds;
                  
//                   double duzdx = mEW->mQ(1,i,j,k)*duzdq+mEW->mR(1,i,j,k)*duzdr+mEW->mS(1,i,j,k)*duzds;

//                   double duydx = mEW->mQ(1,i,j,k)*duydq+mEW->mR(1,i,j,k)*duydr+mEW->mS(1,i,j,k)*duyds;
                  
//                   double duxdy = mEW->mQ(2,i,j,k)*duxdq+mEW->mR(2,i,j,k)*duxdr+mEW->mS(2,i,j,k)*duxds;
                      
//                   /*S.D*/
// //                   if (j==51&&k==1|| (j==101&&k==1&&mEW->mGridSize[g]==0.01))
// //                     {
// //                       printf("duydz  %f %f %f %f %f %f %f\n",mEW->mQ(3,i,j,k)*duydq+mEW->mR(3,i,j,k)*duydr+mEW->mS(3,i,j,k)*duyds,mEW->mQ(3,i,j,k),duydq,mEW->mR(3,i,j,k),duydr,mEW->mS(3,i,j,k),duyds);
// //                     }       
                  
//                   /*E.D*/

//                    if (m_double)
//                      {
//                        // ----------------------------------------------------
//                        // x component: duz/dy - duy/dz
//                        // ----------------------------------------------------
                       
//                        m_doubleField[g][iField] = pow(duzdy-duydz,2);
                       
//                        // ------------------------------------------
//                        // y component: dux/dz - duz/dx
//                        // ----------------------------------------------
                       
//                        m_doubleField[g][iField] += pow(duxdz-duzdx,2);
                       
//                        // ------------------------------------------
//                        // z component: duy/dx - dux/dy
//                        // ------------------------------------------
                       
//                        m_doubleField[g][iField] += pow(duydx-duxdy,2);
                       
//                        m_doubleField[g][iField] = factor*sqrt(m_doubleField[g][iField]);
//                      }
//                    else
//                      {
//                        // ----------------------------------------------------
//                        // x component: duz/dy - duy/dz
//                        // ----------------------------------------------------
                       
//                        m_floatField[g][iField] = (float)pow(duzdy-duydz,2);
                       
// //                        if ((j==51&&k==1&&mEW->mGridSize[g]==0.02) || (j==101&&k==1&&mEW->mGridSize[g]==0.01))
// //                          {
// //                            printf("duzdy duydz duzdy-duydz %f %f %f %i %i %i %f %f %f %f %f %f %f %f\n",mEW->mX(i,j,k),mEW->mY(i,j,k),mEW->mZ(i,j,k),i,j,k,duzdy,duydz,duzdy-duydz,duzdq,duzdr,duzds,duydq,duydr,duyds);
// //                            printf("Q R S %f %f %f %f %f %f %f %f %f \n",mEW->mQ(1,i,j,k),mEW->mQ(2,i,j,k),mEW->mQ(3,i,j,k),mEW->mR(1,i,j,k),mEW->mR(2,i,j,k),mEW->mR(3,i,j,k),mEW->mS(1,i,j,k),mEW->mS(2,i,j,k),mEW->mS(3,i,j,k));
// //                          }         
                       
                       

//                        // ------------------------------------------
//                        // y component: dux/dz - duz/dx
//                        // ----------------------------------------------
                       
//                        m_floatField[g][iField] += (float)pow(duxdz-duzdx,2);
                       
// //                        if ((j==51&&k==1&&mEW->mGridSize[g]==0.02) || (j==101&&k==1&&mEW->mGridSize[g]==0.01))
// //                          printf("duxdz duzdx duxdz-duzdx %f %f %f %i %i %i %f %f %f %f %f %f %f %f %f\n",mEW->mX(i,j,k),mEW->mY(i,j,k),mEW->mZ(i,j,k),i,j,k,duxdz,duzdx,duxdz-duzdx,duxdq,duxdr,duxds,duzdq,duzdr,duzds);

//                        // ------------------------------------------
//                        // z component: duy/dx - dux/dy
//                        // ------------------------------------------
                       
//                        m_floatField[g][iField] += (float)pow(duydx-duxdy,2);
                       
// //                        if ((j==51&&k==1&&mEW->mGridSize[g]==0.02) || (j==101&&k==1&&mEW->mGridSize[g]==0.01))
// //                          printf("duydx duxdy duydx-duxdy %f %f %f %i %i %i %f %f %f %f %f %f %f %f %f\n",mEW->mX(i,j,k),mEW->mY(i,j,k),mEW->mZ(i,j,k),i,j,k,duydx,duxdy,duydx-duxdy,duydq,duydr,duyds,duxdq,duxdr,duxds);

//                        m_floatField[g][iField] = factor*sqrt(m_floatField[g][iField]);

// //                        if ((j==51&&k==1&&mEW->mGridSize[g]==0.02) || (j==101&&k==1&&mEW->mGridSize[g]==0.01))
// //                          printf("floatField %i %i %i %f \n",i,j,k,m_floatField[g][iField]);
//                      }
//                    iField++; 
//                 }
//         }
//     }
// }


// //-----------------------------------------------------------------------
// void Image::computeVelocityMagnitude()
// {
//   ASSERT(m_isDefinedMPIWriters);

// // plane_in_proc returns true for z=const planes, because all processors have a part in these planes
//   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
//   bool breakLoop = (mLocationType == Image::Z)    ;
//   bool done      = false                          ;

//   if (iwrite)
//   {
//       double factor = 1.0/(2.*mEW->mDt);
//       double v1, v2, v3, vmag;
//       for (int g = 0; g < mEW->mNumberOfGrids; g++)
//       {
//         if (breakLoop) 
//         {
//           g = m_gridPtIndex[1];
//         }              

//         int iField = 0;
          
//         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//           for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//             for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//             {
//               v1 = factor*(mEW->mUp[g](1,i,j,k) - mEW->mUm[g](1,i,j,k));
//               v2 = factor*(mEW->mUp[g](2,i,j,k) - mEW->mUm[g](2,i,j,k));
//               v3 = factor*(mEW->mUp[g](3,i,j,k) - mEW->mUm[g](3,i,j,k));

//               vmag = sqrt(v1*v1 + v2*v2 + v3*v3);

//               if (m_double)
//               {
//                 m_doubleField[g][iField] = vmag;
//               }
//               else
//               {
//                 m_floatField[g][iField] = (float) vmag;
//               }
//               iField++; 
//             } // end for i,j,k
//         if (breakLoop) 
//         {
//           break;
//         }
//       } // end for g...
//     } // end if iwrite
// } // end computeVelocityMagnitude()

// //-----------------------------------------------------------------------
// void Image::computeHorizontalVelocityMagnitude()
// {
//   ASSERT(m_isDefinedMPIWriters);

// // plane_in_proc returns true for z=const planes, because all processors have a part in these planes
//   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
//   bool breakLoop = (mLocationType == Image::Z)    ;
//   bool done      = false                          ;

//   if (iwrite)
//   {
//       double factor = 1.0/(2.*mEW->mDt);
//       double v1, v2, v3, vmag;
//       for (int g = 0; g < mEW->mNumberOfGrids; g++)
//       {
//         if (breakLoop) 
//         {
//           g = m_gridPtIndex[1];
//         }              

//         int iField = 0;
          
//         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//           for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//             for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//             {
//               v1 = factor*(mEW->mUp[g](1,i,j,k) - mEW->mUm[g](1,i,j,k));
//               v2 = factor*(mEW->mUp[g](2,i,j,k) - mEW->mUm[g](2,i,j,k));

//               vmag = sqrt(v1*v1 + v2*v2);

//               if (m_double)
//               {
//                 m_doubleField[g][iField] = vmag;
//               }
//               else
//               {
//                 m_floatField[g][iField] = (float) vmag;
//               }
//               iField++; 
//             } // end for i,j,k
//         if (breakLoop) 
//         {
//           break;
//         }
//       } // end for g...
//     } // end if iwrite
// } // end computeHorizontalVelocityMagnitude()

// //-----------------------------------------------------------------------
// void Image::computeVelCurlMagnitude()
// {
//   ASSERT(m_isDefinedMPIWriters);

// // plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
//   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
//   bool breakLoop = (mLocationType == Image::Z)    ;
//   bool done      = false                          ;

//   if (iwrite)
//     {
//       for (int g = 0; g < mEW->mNumberOfCartesianGrids; g++)
//         {
//           if (breakLoop) 
//             {
//               g = m_gridPtIndex[1];
              
//               if (g >= mEW->mNumberOfCartesianGrids) // m_gridPtIndex[1] could equal the curvilinear grid (messy logic)
//                 break;
//             }              

//           int iField = 0;
//           double factor = 1.0/(4*mEW->mGridSize[g]*mEW->mDt);
          
//           for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//             for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//               for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//                 {
//                   double duydx = 0.;
//                   double duzdx = 0.;
//                   double duxdy = 0.;
//                   double duzdy = 0.;
//                   double duxdz = 0.;
//                   double duydz = 0.;
                  
//                   if (i == 1)
//                     {
//                       duydx = ( (-mEW->mUp[g](2,i+2,j,k)+4.*mEW->mUp[g](2,i+1,j,k)-3.*mEW->mUp[g](2,i,j,k))
// 			       -(-mEW->mUm[g](2,i+2,j,k)+4.*mEW->mUm[g](2,i+1,j,k)-3.*mEW->mUm[g](2,i,j,k)));
//                       duzdx = ( (-mEW->mUp[g](3,i+2,j,k)+4.*mEW->mUp[g](3,i+1,j,k)-3.*mEW->mUp[g](3,i,j,k))
// 			       -(-mEW->mUm[g](3,i+2,j,k)+4.*mEW->mUm[g](3,i+1,j,k)-3.*mEW->mUm[g](3,i,j,k)));
//                     }
//                   else if (i == mEW->m_global_nx[g])
//                     {
//                       duydx = ( (3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i-1,j,k)+mEW->mUp[g](2,i-2,j,k))
// 			       -(3.*mEW->mUm[g](2,i,j,k)-4.*mEW->mUm[g](2,i-1,j,k)+mEW->mUm[g](2,i-2,j,k)));
//                       duzdx = ( (3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i-1,j,k)+mEW->mUp[g](3,i-2,j,k))
// 			       -(3.*mEW->mUm[g](3,i,j,k)-4.*mEW->mUm[g](3,i-1,j,k)+mEW->mUm[g](3,i-2,j,k)));
//                     }
//                   else
//                     {
//                       duydx = (mEW->mUp[g](2,i+1,j,k) - mEW->mUp[g](2,i-1,j,k))-(mEW->mUm[g](2,i+1,j,k) - mEW->mUm[g](2,i-1,j,k));
//                       duzdx = (mEW->mUp[g](3,i+1,j,k) - mEW->mUp[g](3,i-1,j,k))-(mEW->mUm[g](3,i+1,j,k) - mEW->mUm[g](3,i-1,j,k));
//                     }
                  
//                   if (j == 1)
//                     {   
//                       duxdy = ( (-mEW->mUp[g](1,i,j+2,k)+4.*mEW->mUp[g](1,i,j+1,k)-3.*mEW->mUp[g](1,i,j,k))
// 			       -(-mEW->mUm[g](1,i,j+2,k)+4.*mEW->mUm[g](1,i,j+1,k)-3.*mEW->mUm[g](1,i,j,k)));
//                       duzdy = ( (-mEW->mUp[g](3,i,j+2,k)+4.*mEW->mUp[g](3,i,j+1,k)-3.*mEW->mUp[g](3,i,j,k))
// 			       -(-mEW->mUm[g](3,i,j+2,k)+4.*mEW->mUm[g](3,i,j+1,k)-3.*mEW->mUm[g](3,i,j,k)));
//                     }
//                   else if ( j == mEW->m_global_ny[g])
//                     { 
//                       duxdy = ( (3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i,j-1,k)+mEW->mUp[g](1,i,j-2,k))
// 			       -(3.*mEW->mUm[g](1,i,j,k)-4.*mEW->mUm[g](1,i,j-1,k)+mEW->mUm[g](1,i,j-2,k)));
//                       duzdy = ( (3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i,j-1,k)+mEW->mUp[g](3,i,j-2,k))
// 			       -(3.*mEW->mUm[g](3,i,j,k)-4.*mEW->mUm[g](3,i,j-1,k)+mEW->mUm[g](3,i,j-2,k)));
//                     }
//                   else
//                     {
//                       duxdy = (mEW->mUp[g](1,i,j+1,k) - mEW->mUp[g](1,i,j-1,k))-(mEW->mUm[g](1,i,j+1,k) - mEW->mUm[g](1,i,j-1,k));
//                       duzdy = (mEW->mUp[g](3,i,j+1,k) - mEW->mUp[g](3,i,j-1,k))-(mEW->mUm[g](3,i,j+1,k) - mEW->mUm[g](3,i,j-1,k));
//                     }
                  
//                   if (k == 1 && g==mEW->mNumberOfGrids-1)
//                     {
//                       duxdz = ( (-mEW->mUp[g](1,i,j,k+2)+4.*mEW->mUp[g](1,i,j,k+1)-3.*mEW->mUp[g](1,i,j,k))
// 			       -(-mEW->mUm[g](1,i,j,k+2)+4.*mEW->mUm[g](1,i,j,k+1)-3.*mEW->mUm[g](1,i,j,k)));
//                       duydz = ( (-mEW->mUp[g](2,i,j,k+2)+4.*mEW->mUp[g](2,i,j,k+1)-3.*mEW->mUp[g](2,i,j,k))
// 			       -(-mEW->mUm[g](2,i,j,k+2)+4.*mEW->mUm[g](2,i,j,k+1)-3.*mEW->mUm[g](2,i,j,k)));
//                     }
// //                  else if (k == mEW->m_kEnd[g]-mEW->m_ghost_points && g == 0)
//                   else if (k == mEW->m_global_nz[g] && g == 0)
//                     {
//                       duxdz = ( (mEW->mUp[g](1,i,j,k-2)-4.*mEW->mUp[g](1,i,j,k-1)+3.*mEW->mUp[g](1,i,j,k))
// 			       -(mEW->mUm[g](1,i,j,k-2)-4.*mEW->mUm[g](1,i,j,k-1)+3.*mEW->mUm[g](1,i,j,k)));
//                       duydz = ( (mEW->mUp[g](2,i,j,k-2)-4.*mEW->mUp[g](2,i,j,k-1)+3.*mEW->mUp[g](2,i,j,k))
// 			       -(mEW->mUm[g](2,i,j,k-2)-4.*mEW->mUm[g](2,i,j,k-1)+3.*mEW->mUm[g](2,i,j,k)));
//                     }
//                   else // no k=N because we are in the curvilinear grid
//                     {
//                       duxdz = (mEW->mUp[g](1,i,j,k+1) - mEW->mUp[g](1,i,j,k-1))-(mEW->mUm[g](1,i,j,k+1) - mEW->mUm[g](1,i,j,k-1));
//                       duydz = (mEW->mUp[g](2,i,j,k+1) - mEW->mUp[g](2,i,j,k-1))-(mEW->mUm[g](2,i,j,k+1) - mEW->mUm[g](2,i,j,k-1));
//                     }
             
//                   if (m_double)
//                     {
//                       m_doubleField[g][iField] = (duzdy-duydz)*(duzdy-duydz)+(duxdz-duzdx)*(duxdz-duzdx)+(duydx-duxdy)*(duydx-duxdy);
                      
//                       m_doubleField[g][iField] = factor*sqrt(m_doubleField[g][iField]);
//                     }
//                   else
//                     {
//                       m_floatField[g][iField] = (duzdy-duydz)*(duzdy-duydz)+(duxdz-duzdx)*(duxdz-duzdx)+(duydx-duxdy)*(duydx-duxdy);
                      
//                       m_floatField[g][iField] = factor*sqrt(m_floatField[g][iField]);

// //                       printf("curlSq curl1 curl2 curl3 %f %f %f %f %f %f %f \n",mEW->mX(i,j,k),mEW->mY(i,j,k),mEW->mZ(i,j,k),m_floatField[g][iField],factor*(duzdy-duydz),factor*(duxdz-duzdx),factor*(duydx-duxdy));
//                     }
//                   iField++; 
//                 }
//           if (breakLoop) 
//             {
//               done = true;
//               break;
//             }
//         }
//       if (mEW->topographyExists() && !done)
//         {   
//           int g = mEW->mNumberOfGrids - 1;
//           int iField = 0;
//           double factor = 1.0/(4*mEW->mDt);
          
//           for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//             for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//               for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//                 {
//                   double duxdq = 0.;
//                   double duydq = 0.;
//                   double duzdq = 0.;
//                   double duxdr = 0.;
//                   double duydr = 0.;
//                   double duzdr = 0.;
//                   double duxds = 0.;
//                   double duyds = 0.;
//                   double duzds = 0.;

//                   if (i == 1)
//                     {      
//                       duxdq = ((-mEW->mUp[g](1,i+2,j,k)+4.*mEW->mUp[g](1,i+1,j,k)-3.*mEW->mUp[g](1,i,j,k))
// 			       -(-mEW->mUm[g](1,i+2,j,k)+4.*mEW->mUm[g](1,i+1,j,k)-3.*mEW->mUm[g](1,i,j,k)));
//                       duydq = ((-mEW->mUp[g](2,i+2,j,k)+4.*mEW->mUp[g](2,i+1,j,k)-3.*mEW->mUp[g](2,i,j,k))
// 			       -(-mEW->mUm[g](2,i+2,j,k)+4.*mEW->mUm[g](2,i+1,j,k)-3.*mEW->mUm[g](2,i,j,k)));
//                       duzdq = ((-mEW->mUp[g](3,i+2,j,k)+4.*mEW->mUp[g](3,i+1,j,k)-3.*mEW->mUp[g](3,i,j,k))
// 			       -(-mEW->mUm[g](3,i+2,j,k)+4.*mEW->mUm[g](3,i+1,j,k)-3.*mEW->mUm[g](3,i,j,k)));
//                     }
//                   else if (i == mEW->m_global_nx[g])
//                     { 
//                       duxdq = ((3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i-1,j,k)+mEW->mUp[g](1,i-2,j,k))
// 			       -(3.*mEW->mUm[g](1,i,j,k)-4.*mEW->mUm[g](1,i-1,j,k)+mEW->mUm[g](1,i-2,j,k)));
//                       duydq = ((3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i-1,j,k)+mEW->mUp[g](2,i-2,j,k))
// 			       -(3.*mEW->mUm[g](2,i,j,k)-4.*mEW->mUm[g](2,i-1,j,k)+mEW->mUm[g](2,i-2,j,k)));
//                       duzdq = ((3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i-1,j,k)+mEW->mUp[g](3,i-2,j,k))
// 			       -(3.*mEW->mUm[g](3,i,j,k)-4.*mEW->mUm[g](3,i-1,j,k)+mEW->mUm[g](3,i-2,j,k)));
//                     }
//                   else
//                     { 
//                       duxdq = (mEW->mUp[g](1,i+1,j,k) - mEW->mUp[g](1,i-1,j,k))-(mEW->mUm[g](1,i+1,j,k) - mEW->mUm[g](1,i-1,j,k));
//                       duydq = (mEW->mUp[g](2,i+1,j,k) - mEW->mUp[g](2,i-1,j,k))-(mEW->mUm[g](2,i+1,j,k) - mEW->mUm[g](2,i-1,j,k));
//                       duzdq = (mEW->mUp[g](3,i+1,j,k) - mEW->mUp[g](3,i-1,j,k))-(mEW->mUm[g](3,i+1,j,k) - mEW->mUm[g](3,i-1,j,k));
//                     }
                  
//                   if (j == 1)
//                     {   
//                       duxdr = ((-mEW->mUp[g](1,i,j+2,k)+4.*mEW->mUp[g](1,i,j+1,k)-3.*mEW->mUp[g](1,i,j,k))
// 			       -(-mEW->mUm[g](1,i,j+2,k)+4.*mEW->mUm[g](1,i,j+1,k)-3.*mEW->mUm[g](1,i,j,k)));
//                       duydr = ((-mEW->mUp[g](2,i,j+2,k)+4.*mEW->mUp[g](2,i,j+1,k)-3.*mEW->mUp[g](2,i,j,k))
// 			       -(-mEW->mUm[g](2,i,j+2,k)+4.*mEW->mUm[g](2,i,j+1,k)-3.*mEW->mUm[g](2,i,j,k)));
//                       duzdr = ((-mEW->mUp[g](3,i,j+2,k)+4.*mEW->mUp[g](3,i,j+1,k)-3.*mEW->mUp[g](3,i,j,k))
// 			       -(-mEW->mUm[g](3,i,j+2,k)+4.*mEW->mUm[g](3,i,j+1,k)-3.*mEW->mUm[g](3,i,j,k)));
//                     }
//                   else if (j == mEW->m_global_ny[g])
//                     { 
//                       duxdr = ((3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i,j-1,k)+mEW->mUp[g](1,i,j-2,k))
// 			       -(3.*mEW->mUm[g](1,i,j,k)-4.*mEW->mUm[g](1,i,j-1,k)+mEW->mUm[g](1,i,j-2,k)));
//                       duydr = ((3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i,j-1,k)+mEW->mUp[g](2,i,j-2,k))
// 			       -(3.*mEW->mUm[g](2,i,j,k)-4.*mEW->mUm[g](2,i,j-1,k)+mEW->mUm[g](2,i,j-2,k)));
//                       duzdr = ((3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i,j-1,k)+mEW->mUp[g](3,i,j-2,k))
// 			       -(3.*mEW->mUm[g](3,i,j,k)-4.*mEW->mUm[g](3,i,j-1,k)+mEW->mUm[g](3,i,j-2,k)));
//                     }
//                   else
//                     {
//                       duxdr = (mEW->mUp[g](1,i,j+1,k) - mEW->mUp[g](1,i,j-1,k))-(mEW->mUm[g](1,i,j+1,k) - mEW->mUm[g](1,i,j-1,k));
//                       duydr = (mEW->mUp[g](2,i,j+1,k) - mEW->mUp[g](2,i,j-1,k))-(mEW->mUm[g](2,i,j+1,k) - mEW->mUm[g](2,i,j-1,k));
//                       duzdr = (mEW->mUp[g](3,i,j+1,k) - mEW->mUp[g](3,i,j-1,k))-(mEW->mUm[g](3,i,j+1,k) - mEW->mUm[g](3,i,j-1,k));
//                     }
                  
//                   if (k == 1)
//                     {
//                       duxds = ((-mEW->mUp[g](1,i,j,k+2)+4.*mEW->mUp[g](1,i,j,k+1)-3.*mEW->mUp[g](1,i,j,k))
// 			       -(-mEW->mUm[g](1,i,j,k+2)+4.*mEW->mUm[g](1,i,j,k+1)-3.*mEW->mUm[g](1,i,j,k)));
//                       duyds = ((-mEW->mUp[g](2,i,j,k+2)+4.*mEW->mUp[g](2,i,j,k+1)-3.*mEW->mUp[g](2,i,j,k))
// 			       -(-mEW->mUm[g](2,i,j,k+2)+4.*mEW->mUm[g](2,i,j,k+1)-3.*mEW->mUm[g](2,i,j,k)));
//                       duzds = ((-mEW->mUp[g](3,i,j,k+2)+4.*mEW->mUp[g](3,i,j,k+1)-3.*mEW->mUp[g](3,i,j,k))
// 			       -(-mEW->mUm[g](3,i,j,k+2)+4.*mEW->mUm[g](3,i,j,k+1)-3.*mEW->mUm[g](3,i,j,k)));
//                     }
//                   else // no k=N because we are in the curvilinear grid
//                     {
//                       duxds = (mEW->mUp[g](1,i,j,k+1) - mEW->mUp[g](1,i,j,k-1))-(mEW->mUm[g](1,i,j,k+1) - mEW->mUm[g](1,i,j,k-1));
//                       duyds = (mEW->mUp[g](2,i,j,k+1) - mEW->mUp[g](2,i,j,k-1))-(mEW->mUm[g](2,i,j,k+1) - mEW->mUm[g](2,i,j,k-1));
//                       duzds = (mEW->mUp[g](3,i,j,k+1) - mEW->mUp[g](3,i,j,k-1))-(mEW->mUm[g](3,i,j,k+1) - mEW->mUm[g](3,i,j,k-1));
//                     }                  

//                   double duzdy = mEW->mQ(2,i,j,k)*duzdq+mEW->mR(2,i,j,k)*duzdr+mEW->mS(2,i,j,k)*duzds;
                  
//                   double duydz = mEW->mQ(3,i,j,k)*duydq+mEW->mR(3,i,j,k)*duydr+mEW->mS(3,i,j,k)*duyds;
                  
//                   double duxdz = mEW->mQ(3,i,j,k)*duxdq+mEW->mR(3,i,j,k)*duxdr+mEW->mS(3,i,j,k)*duxds;
                  
//                   double duzdx = mEW->mQ(1,i,j,k)*duzdq+mEW->mR(1,i,j,k)*duzdr+mEW->mS(1,i,j,k)*duzds;

//                   double duydx = mEW->mQ(1,i,j,k)*duydq+mEW->mR(1,i,j,k)*duydr+mEW->mS(1,i,j,k)*duyds;
                  
//                   double duxdy = mEW->mQ(2,i,j,k)*duxdq+mEW->mR(2,i,j,k)*duxdr+mEW->mS(2,i,j,k)*duxds;
                      
//                    if (m_double)
//                      {
//                        // ----------------------------------------------------
//                        // x component: duz/dy - duy/dz
//                        // ----------------------------------------------------
                       
//                        m_doubleField[g][iField] = pow(duzdy-duydz,2);
                       
//                        // ------------------------------------------
//                        // y component: dux/dz - duz/dx
//                        // ----------------------------------------------
                       
//                        m_doubleField[g][iField] += pow(duxdz-duzdx,2);
                       
//                        // ------------------------------------------
//                        // z component: duy/dx - dux/dy
//                        // ------------------------------------------
                       
//                        m_doubleField[g][iField] += pow(duydx-duxdy,2);
                       
//                        m_doubleField[g][iField] = factor*sqrt(m_doubleField[g][iField]);
//                      }
//                    else
//                      {
//                        // ----------------------------------------------------
//                        // x component: duz/dy - duy/dz
//                        // ----------------------------------------------------
                       
//                        m_floatField[g][iField] = (float)pow(duzdy-duydz,2);
                       
//                        // ------------------------------------------
//                        // y component: dux/dz - duz/dx
//                        // ----------------------------------------------
                       
//                        m_floatField[g][iField] += (float)pow(duxdz-duzdx,2);
                       
//                        // ------------------------------------------
//                        // z component: duy/dx - dux/dy
//                        // ------------------------------------------
                       
//                        m_floatField[g][iField] += (float)pow(duydx-duxdy,2);
                       
//                        m_floatField[g][iField] = factor*sqrt(m_floatField[g][iField]); 
                       
// //                        if (i==31&&j==1&&k==1)
// //                          {
// //                            printf("sanity check %i %i %i %f %f %f %f \n",i,j,k,pow(duzdy-duydz,2),pow(duxdz-duzdx,2),pow(duydx-duxdy,2),m_floatField[g][iField]);
// //                            printf("sanity check %i %i %i %f %f %f %f %f %f %f %f \n",i,j,k,duzdy,duydz,duzdq,duzdr,duzds,duydq,duydr,duyds);
// //                          }
//                      }
//                    iField++;
//                 }
//         }
//     }
// }

// //-----------------------------------------------------------------------
// void Image::computeVelDivergence()
// {
//   ASSERT(m_isDefinedMPIWriters);

// // plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
//   bool iwrite    = plane_in_proc(m_gridPtIndex[0]);
//   bool breakLoop = (mLocationType == Image::Z)    ;
//   bool done      = false                          ;

//   if (iwrite)
//     {
//       for (int g = 0; g < mEW->mNumberOfCartesianGrids; g++)
//         {
//           if (breakLoop) 
//             {
//               g = m_gridPtIndex[1];
              
//               if (g >= mEW->mNumberOfCartesianGrids)
//                 break;
//             }

//           double factor = 1.0/(4*mEW->mDt*mEW->mGridSize[g]);
//           int iField = 0;
          
//           for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//             for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//               for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//                 {
//                   double val = 0.;

//                   if (i == 1)
//                     {
//                       val  = ((-mEW->mUp[g](1,i+2,j,k)+4.*mEW->mUp[g](1,i+1,j,k)-3.*mEW->mUp[g](1,i,j,k))
//                              -(-mEW->mUm[g](1,i+2,j,k)+4.*mEW->mUm[g](1,i+1,j,k)-3.*mEW->mUm[g](1,i,j,k)));
//                     }
//                   else if (i == mEW->m_global_nx[g])
//                     {
//                       val  = ((3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i-1,j,k)+mEW->mUp[g](1,i-2,j,k))
//                              -(3.*mEW->mUm[g](1,i,j,k)-4.*mEW->mUm[g](1,i-1,j,k)+mEW->mUm[g](1,i-2,j,k)));
//                     }
//                   else
//                     {
//                       val  = ((mEW->mUp[g](1,i+1, j, k) - mEW->mUp[g](1,i-1, j, k))
//                              -(mEW->mUm[g](1,i+1, j, k) - mEW->mUm[g](1,i-1, j, k)));
//                     }
                  
//                   if (j == 1)
//                     {
//                       val  += ((-mEW->mUp[g](2,i,j+2,k)+4.*mEW->mUp[g](2,i,j+1,k)-3.*mEW->mUp[g](2,i,j,k))
//                               -(-mEW->mUm[g](2,i,j+2,k)+4.*mEW->mUm[g](2,i,j+1,k)-3.*mEW->mUm[g](2,i,j,k)));
//                     }
//                   else if (j == mEW->m_global_ny[g])
//                     {
//                       val  += ((3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i,j-1,k)+mEW->mUp[g](2,i,j-2,k))
//                               -(3.*mEW->mUm[g](2,i,j,k)-4.*mEW->mUm[g](2,i,j-1,k)+mEW->mUm[g](2,i,j-2,k)));
//                     }
//                   else
//                     {
//                       val  += ((mEW->mUp[g](2,i, j+1, k) - mEW->mUp[g](2,i,j-1,k))
//                               -(mEW->mUm[g](2,i, j+1, k) - mEW->mUm[g](2,i,j-1,k)));
//                     }
                  
//                   if (k == 1&& g==mEW->mNumberOfGrids-1)
//                     {
//                       val  += ((-mEW->mUp[g](3,i,j,k+2)+4.*mEW->mUp[g](3,i,j,k+1)-3.*mEW->mUp[g](3,i,j,k))
//                               -(-mEW->mUm[g](3,i,j,k+2)+4.*mEW->mUm[g](3,i,j,k+1)-3.*mEW->mUm[g](3,i,j,k)));
//                     }
// //                  else if (k == mEW->m_kEnd[g]-mEW->m_ghost_points && g == 0)
//                   else if (k == mEW->m_global_nz[g] && g == 0)
//                     {
//                       val  += ((mEW->mUp[g](3,i,j,k-2)-4.*mEW->mUp[g](3,i,j,k-1)+3.*mEW->mUp[g](3,i,j,k))
//                               -(mEW->mUm[g](3,i,j,k-2)-4.*mEW->mUm[g](3,i,j,k-1)+3.*mEW->mUm[g](3,i,j,k)));
//                     }
//                   else // no k=N because we are in the curvilinear grid
//                     {
//                       val  += ((mEW->mUp[g](3,i,j,k+1) - mEW->mUp[g](3,i,j,k-1))
//                               -(mEW->mUm[g](3,i,j,k+1) - mEW->mUm[g](3,i,j,k-1)));
//                     }
                      
//                   val *= factor;
                  

//                   if (m_double)
//                     {
//                       m_doubleField[g][iField] = val;
//                     }
//                   else
//                     {
//                       m_floatField[g][iField] = val;
//                     }
//                   iField++; 
//                 }
//           if (breakLoop) 
//             {
//               done = true;
//               break;
//             }
//         }
      
//       // curvilinear grid
//       if (mEW->topographyExists() && !done)
//         {
//           int g = mEW->mNumberOfGrids - 1;
//           double factor = 1.0/(4*mEW->mDt);
          
//           int iField = 0;
          
//           for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
//             for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
//               for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
//                 {
//                   double val = 0.;

//                   if (i == 1)
//                     {
//                       val  = ((mEW->mQ(1,i,j,k)*((-mEW->mUp[g](1,i+2,j,k)+4.*mEW->mUp[g](1,i+1,j,k)-3.*mEW->mUp[g](1,i,j,k))
//                                                  -(-mEW->mUm[g](1,i+2,j,k)+4.*mEW->mUm[g](1,i+1,j,k)-3.*mEW->mUm[g](1,i,j,k)))
//                               +mEW->mQ(2,i,j,k)*((-mEW->mUp[g](2,i+2,j,k)+4.*mEW->mUp[g](2,i+1,j,k)-3.*mEW->mUp[g](2,i,j,k))
//                                                  -(-mEW->mUm[g](2,i+2,j,k)+4.*mEW->mUm[g](2,i+1,j,k)-3.*mEW->mUm[g](2,i,j,k)))
//                               +mEW->mQ(3,i,j,k)*((-mEW->mUp[g](3,i+2,j,k)+4.*mEW->mUp[g](3,i+1,j,k)-3.*mEW->mUp[g](3,i,j,k))
//                                                  -(-mEW->mUm[g](3,i+2,j,k)+4.*mEW->mUm[g](3,i+1,j,k)-3.*mEW->mUm[g](3,i,j,k)))));
//                     }
//                   else if (i == mEW->m_global_nx[g])
//                     {
//                       val  = ((mEW->mQ(1,i,j,k)*((3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i-1,j,k)+mEW->mUp[g](1,i-2,j,k))
//                                                  -(3.*mEW->mUm[g](1,i,j,k)-4.*mEW->mUm[g](1,i-1,j,k)+mEW->mUm[g](1,i-2,j,k)))
//                               +mEW->mQ(2,i,j,k)*((3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i-1,j,k)+mEW->mUp[g](2,i-2,j,k))
//                                                  -(3.*mEW->mUm[g](2,i,j,k)-4.*mEW->mUm[g](2,i-1,j,k)+mEW->mUm[g](2,i-2,j,k)))
//                               +mEW->mQ(3,i,j,k)*((3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i-1,j,k)+mEW->mUp[g](3,i-2,j,k))
//                                                  -(3.*mEW->mUm[g](3,i,j,k)-4.*mEW->mUm[g](3,i-1,j,k)+mEW->mUm[g](3,i-2,j,k)))));
//                     }
//                   else
//                     {
//                       val  = (mEW->mQ(1,i,j,k)*((mEW->mUp[g](1,i+1,j,k) - mEW->mUp[g](1,i-1,j,k))
//                                                 -(mEW->mUm[g](1,i+1,j,k) - mEW->mUm[g](1,i-1,j,k)))
//                              +mEW->mQ(2,i,j,k)*((mEW->mUp[g](2,i+1,j,k) - mEW->mUp[g](2,i-1,j,k))
//                                                 -(mEW->mUm[g](2,i+1,j,k) - mEW->mUm[g](2,i-1,j,k)))
//                              +mEW->mQ(3,i,j,k)*((mEW->mUp[g](3,i+1,j,k) - mEW->mUp[g](3,i-1,j,k))
//                                                 -(mEW->mUm[g](3,i+1,j,k) - mEW->mUm[g](3,i-1,j,k))));
//                     }
                  
//                   if (j == 1)
//                     {
//                       val  += ((mEW->mR(1,i,j,k)*((-mEW->mUp[g](1,i,j+2,k)+4.*mEW->mUp[g](1,i,j+1,k)-3.*mEW->mUp[g](1,i,j,k))
//                                                   -(-mEW->mUm[g](1,i,j+2,k)+4.*mEW->mUm[g](1,i,j+1,k)-3.*mEW->mUm[g](1,i,j,k)))
//                                +mEW->mR(2,i,j,k)*((-mEW->mUp[g](2,i,j+2,k)+4.*mEW->mUp[g](2,i,j+1,k)-3.*mEW->mUp[g](2,i,j,k))
//                                                   -(-mEW->mUm[g](2,i,j+2,k)+4.*mEW->mUm[g](2,i,j+1,k)-3.*mEW->mUm[g](2,i,j,k)))
//                                +mEW->mR(3,i,j,k)*((-mEW->mUp[g](3,i,j+2,k)+4.*mEW->mUp[g](3,i,j+1,k)-3.*mEW->mUp[g](3,i,j,k))
//                                                   -(-mEW->mUm[g](3,i,j+2,k)+4.*mEW->mUm[g](3,i,j+1,k)-3.*mEW->mUm[g](3,i,j,k)))));
//                     }
//                   else if (j == mEW->m_global_ny[g])
//                     {
//                       val  += ((mEW->mR(1,i,j,k)*((3.*mEW->mUp[g](1,i,j,k)-4.*mEW->mUp[g](1,i,j-1,k)+mEW->mUp[g](1,i,j-2,k))
//                                                   -(3.*mEW->mUm[g](1,i,j,k)-4.*mEW->mUm[g](1,i,j-1,k)+mEW->mUm[g](1,i,j-2,k)))
//                                +mEW->mR(2,i,j,k)*((3.*mEW->mUp[g](2,i,j,k)-4.*mEW->mUp[g](2,i,j-1,k)+mEW->mUp[g](2,i,j-2,k))
//                                                   -(3.*mEW->mUm[g](2,i,j,k)-4.*mEW->mUm[g](2,i,j-1,k)+mEW->mUm[g](2,i,j-2,k)))
//                                +mEW->mR(3,i,j,k)*((3.*mEW->mUp[g](3,i,j,k)-4.*mEW->mUp[g](3,i,j-1,k)+mEW->mUp[g](3,i,j-2,k))
//                                                   -(3.*mEW->mUm[g](3,i,j,k)-4.*mEW->mUm[g](3,i,j-1,k)+mEW->mUm[g](3,i,j-2,k)))));
//                     }
//                   else
//                     {
//                       val  += (mEW->mR(1,i,j,k)*((mEW->mUp[g](1,i,j+1,k) - mEW->mUp[g](1,i,j-1,k))
//                                                  -(mEW->mUm[g](1,i,j+1,k) - mEW->mUm[g](1,i,j-1,k)))
//                               +mEW->mR(2,i,j,k)*((mEW->mUp[g](2,i,j+1,k) - mEW->mUp[g](2,i,j-1,k))
//                                                  -(mEW->mUm[g](2,i,j+1,k) - mEW->mUm[g](2,i,j-1,k)))
//                               +mEW->mR(3,i,j,k)*((mEW->mUp[g](3,i,j+1,k) - mEW->mUp[g](3,i,j-1,k))
//                                                  -(mEW->mUm[g](3,i,j+1,k) - mEW->mUm[g](3,i,j-1,k))));
//                     }
                  
//                   if (k == 1)
//                     {
//                       val  += ((mEW->mS(1,i,j,k)*((-mEW->mUp[g](1,i,j,k+2)+4.*mEW->mUp[g](1,i,j,k+1)-3.*mEW->mUp[g](1,i,j,k))
//                                                   -(-mEW->mUm[g](1,i,j,k+2)+4.*mEW->mUm[g](1,i,j,k+1)-3.*mEW->mUm[g](1,i,j,k)))
//                                +mEW->mS(2,i,j,k)*((-mEW->mUp[g](2,i,j,k+2)+4.*mEW->mUp[g](2,i,j,k+1)-3.*mEW->mUp[g](2,i,j,k))
//                                                   -(-mEW->mUm[g](2,i,j,k+2)+4.*mEW->mUm[g](2,i,j,k+1)-3.*mEW->mUm[g](2,i,j,k)))
//                                +mEW->mS(3,i,j,k)*((-mEW->mUp[g](3,i,j,k+2)+4.*mEW->mUp[g](3,i,j,k+1)-3.*mEW->mUp[g](3,i,j,k))
//                                                   -(-mEW->mUm[g](3,i,j,k+2)+4.*mEW->mUm[g](3,i,j,k+1)-3.*mEW->mUm[g](3,i,j,k)))));
//                         }
//                   else // no k=N because we are in the curvilinear grid
//                     {
//                       val  += (mEW->mS(1,i,j,k)*((mEW->mUp[g](1,i,j,k+1) - mEW->mUp[g](1,i,j,k-1))
//                                                  -(mEW->mUm[g](1,i,j,k+1) - mEW->mUm[g](1,i,j,k-1)))
//                               +mEW->mS(2,i,j,k)*((mEW->mUp[g](2,i,j,k+1) - mEW->mUp[g](2,i,j,k-1))
//                                                  -(mEW->mUm[g](2,i,j,k+1) - mEW->mUm[g](2,i,j,k-1)))
//                               +mEW->mS(3,i,j,k)*((mEW->mUp[g](3,i,j,k+1) - mEW->mUp[g](3,i,j,k-1))
//                                                  -(mEW->mUm[g](3,i,j,k+1) - mEW->mUm[g](3,i,j,k-1))));
//                     }
                      
//                   val *= factor;  

//                   if( m_double )
//                     {
//                       m_doubleField[g][iField] = val;
//                     }
//                   else
//                     {
//                       m_floatField[g][iField] = (float)val;  
//                     }
//                   iField++;
//                 }
//         } // end if curvilinear grid
//     }
// }


//C.B> For debugging purposes.

//-----------------------------------------------------------------------
// double Image::computeImageErrorDebug(int whichCase)
// {
//   ASSERT(m_isDefinedMPIWriters);
//   if (!mEW->m_forcing->knows_exact())
//   {
//     if (proc_write())
//       cout << "Image::computeImageError can only be called when the exact solution is known" << endl;
//     MPI_Abort(MPI_COMM_WORLD, 1);
//   }
  
//   double maxerrlocal = 0.;
//   double maxerrlocalreg = 0.;
//   double maxerrlocalcurv = 0.;
//   double maxdivlocal = 0.;
//   double maxdivClocal = 0.;
//   double maxdivCCurvlocal = 0.;
//   double maxdivReglocal = 0.;
//   double errl2 = 0.;
  
//   // plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
//   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);

//   if (iwrite)
//     {
//       int myRank;
//       MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

//       bool breakLoop      = (mLocationType == Image::Z);

//       double err, x, y, z, t=mEW->getCurrentTime();
//       double uExact;
//       if (whichCase == 1)
//         t -= mEW->mDt;
//       else if (whichCase == 3)
//         t -= mEW->mDt;
// // tmp
//       if (proc_write())
//       {
//  	cout << endl << "Computing image error for image " << mFilePrefix << " of mode " << mMode2Suffix[mMode] << 
// 	  " at location " << mOrientationString[mLocationType] << "=" << mCoordValue << " and time=" << t << "..." << endl;
// //	  "...msg from proc # " << m_rankWriter << endl; 
//       }

//       // write the data...
//       int g;
      
//       for (g = 0; g < mEW->mNumberOfCartesianGrids; g++)
//       {
//         if (breakLoop) break;
        
// 	int iField = 0;
              
// 	for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
//           for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
//             for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
// 	    {                   
// 	      x = (ii - 1)*mEW->mGridSize[g];
// 	      y = (jj - 1)*mEW->mGridSize[g];
// 	      z = (kk - 1)*mEW->mGridSize[g] + mEW->m_zmin[g];
              
//               if (whichCase == 0)
//                 mEW->m_forcing->get_exactDiv(x, y, z, t, uExact, mEW->mGridSize[g])    ;   
//               else if (whichCase == 1)
//                   mEW->m_forcing->get_exactVelDiv(x, y, z, t, uExact, mEW->mGridSize[g]) ;
//               else if (whichCase == 2)
//                 mEW->m_forcing->get_exactCurl(x, y, z, t, uExact, mEW->mGridSize[g])   ;
//               else if (whichCase == 3)
//                 mEW->m_forcing->get_exactVelCurl(x, y, z, t, uExact, mEW->mGridSize[g]);
//               else
//                 {
//                   printf("Case not accounted for.");
//                   MPI_Abort(MPI_COMM_WORLD,1);
//                 }

//               if (m_double)
//                 {
//                   err = (m_doubleField[g][iField] - uExact );
//                   maxdivClocal = (fabs(m_doubleField[g][iField]) > maxdivClocal ? fabs(m_doubleField[g][iField]) : maxdivClocal);
//                   maxdivReglocal = (fabs(m_doubleField[g][iField]) > maxdivReglocal ? fabs(m_doubleField[g][iField]) : maxdivReglocal);
//                   maxdivlocal = (fabs(uExact) > maxdivlocal ? fabs(uExact) : maxdivlocal);
//                 }              
//               else
//                 {
//                   err = m_floatField[g][iField]-uExact;
//                   maxdivClocal = (fabs(m_floatField[g][iField]) > maxdivClocal ? fabs(m_floatField[g][iField]) : maxdivClocal);
//                   maxdivReglocal = (fabs(m_floatField[g][iField]) > maxdivReglocal ? fabs(m_floatField[g][iField]) : maxdivReglocal);
//                   maxdivlocal = (fabs(uExact) > maxdivlocal ? fabs(uExact) : maxdivlocal);
//                 }
//               maxerrlocal = (fabs(err) > maxerrlocal ? fabs(err) : maxerrlocal);
//               maxerrlocalreg = (fabs(err) > maxerrlocalreg ? fabs(err) : maxerrlocalreg);
//               errl2 += err*err*pow(mEW->mGridSize[g],3);
             
// //               printf("reg err  %f %f %f %i %i %i %i\n",err,m_floatField[g][iField],uExact,ii,jj,kk,g);

//               iField++;
// 	    }
//       } // for g...

// // curvilinear grid
//       if (mEW->topographyExists())
//       {
// 	g = mEW->mNumberOfGrids - 1;
	
// 	int iField = 0;
              
// 	for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
// 	  for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
// 	    for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
// 	    {                   
// 	      x = mEW->mX(ii,jj,kk);
// 	      y = mEW->mY(ii,jj,kk);
// 	      z = mEW->mZ(ii,jj,kk);

//               if (whichCase == 0)
//                 mEW->m_forcing->get_exactDiv(x, y, z, t, uExact, mEW->mGridSize[g])    ;   
//               else if (whichCase == 1)
//                 mEW->m_forcing->get_exactVelDiv(x, y, z, t, uExact, mEW->mGridSize[g]) ;
//               else if (whichCase == 2)
//                 mEW->m_forcing->get_exactCurl(x, y, z, t, uExact, mEW->mGridSize[g])   ;
//               else if (whichCase == 3)
//                 mEW->m_forcing->get_exactVelCurl(x, y, z, t, uExact, mEW->mGridSize[g]);
//               else
//                 {
//                   printf("Case not accounted for.");
//                   MPI_Abort(MPI_COMM_WORLD,1);
//                 }
		  
//               maxdivlocal = (fabs(uExact) > maxdivlocal ? fabs(uExact) : maxdivlocal);
              
//               if (m_double)
//                 {
//                   err = (m_doubleField[g][iField] - uExact );
//                   maxdivClocal = (fabs(m_doubleField[g][iField]) > maxdivClocal ? fabs(m_doubleField[g][iField]) : maxdivClocal);
//                   maxdivCCurvlocal = (fabs(m_doubleField[g][iField]) > maxdivCCurvlocal ? fabs(m_doubleField[g][iField]) : maxdivCCurvlocal);
//                 }
//               else
//                 {
//                   err = m_floatField[g][iField]-uExact;
//                   maxdivClocal = (fabs(m_floatField[g][iField]) > maxdivClocal ? fabs(m_floatField[g][iField]) : maxdivClocal);
//                   maxdivCCurvlocal = (fabs(m_floatField[g][iField]) > maxdivCCurvlocal ? fabs(m_floatField[g][iField]) : maxdivCCurvlocal);
//                 }                  
              
//               maxerrlocal = (fabs(err) > maxerrlocal ? fabs(err) : maxerrlocal);
//               maxerrlocalcurv = (fabs(err) > maxerrlocalcurv ? fabs(err) : maxerrlocalcurv);
//               errl2 += err*err*pow(mEW->mGridSize[g],3);
              
//               iField++;
//             }
//       } // end if topographyExists
      
      
//       // tmp
// //      printf("Image::maxerr=%e\n", maxerr);
      
//     } // end if iwrite

//   double maxerr     = 0.;
//   double maxerrreg  = 0.;
//   double maxerrcurv = 0.;
//   double maxdiv = 0.;
//   double maxdivC = 0.;
//   double maxdivCCurv = 0.;
//   double maxdivReg = 0.;
//   MPI_Allreduce( &maxerrlocal, &maxerr, 1, MPI_DOUBLE, MPI_MAX, mEW->m_cartesian_communicator );  
//   MPI_Allreduce( &maxerrlocalreg, &maxerrreg, 1, MPI_DOUBLE, MPI_MAX, mEW->m_cartesian_communicator );  
//   MPI_Allreduce( &maxerrlocalcurv, &maxerrcurv, 1, MPI_DOUBLE, MPI_MAX, mEW->m_cartesian_communicator );  
//   MPI_Allreduce( &maxdivlocal, &maxdiv, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)                 ;
//   MPI_Allreduce( &maxdivClocal, &maxdivC, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)                 ;
//   MPI_Allreduce( &maxdivCCurvlocal, &maxdivCCurv, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)                 ;
//   MPI_Allreduce( &maxdivReglocal, &maxdivReg, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)                 ;
  
//   if (proc_write())
//     {
//       printf("max Error Cartesian:   %e \n",maxerrreg);
//       printf("max Error Curvilinear: %e \n",maxerrcurv);
// //       printf("maxDivExact          %f \n",maxdiv);
// //       printf("maxDivComputed       %f \n",maxdivC);
// //       printf("maxDivCCurvComputed  %f \n",maxdivCCurv);
// //       printf("maxDivRegComputed    %f \n",maxdivReg);
//       printf("L2-Error:              %e \n",sqrt(errl2));
//     }
  
//   return maxerr;
// } // end Image::computeImageError()

