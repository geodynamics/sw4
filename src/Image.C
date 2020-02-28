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
#include "Require.h"
#include "Image.h"
#include <cstdio>
#include <cmath>
#include <fcntl.h>
#include <ctime>
#include <cstring>
#include <unistd.h>

#ifdef USE_HDF5
#include "sachdf5.h"
#endif

// initializing static member
int Image::mPreceedZeros=0;

int Image::MODES=40;

Image* Image::nil=static_cast<Image*>(0);

using namespace std;

Image::Image(EW * a_ew,
             float_sw4 time, 
             float_sw4 timeInterval, 
             int cycle, 
             int cycleInterval,
             const std::string& filePrefix, 
             ImageMode mode,
	     ImageOrientation locationType, 
             float_sw4 locationValue,
	     bool doubleMode,
	     bool usehdf5,
	     bool userCreated ):
  mTime(time),
  mEW(a_ew),
  m_time_done(false),
  mTimeInterval(timeInterval),
  mWritingCycle(cycle),
  mCycleInterval(cycleInterval),
  mFilePrefix(filePrefix),
  mMode(mode),
  //  mFileName(""),
  //  m_gridPtValueInitialized(false),
  //  mWriting(false),
  //  mReadyToWrite(false),
  mLocationType(locationType),
  mCoordValue(locationValue),
  //  m_isDefined(false),
  m_mpiComm_writers(MPI_COMM_NULL),
  m_isDefinedMPIWriters(false),
  m_double(doubleMode),
  m_usehdf5(usehdf5),
  mGridinfo(-1),
  mStoreGrid(true),
  m_gridimage(Image::nil),
  m_write_time(0.0),
  m_user_created(userCreated)
{
  mMode2Suffix.resize(MODES);
  mMode2Suffix[NONE] = "unknown";
  mMode2Suffix[UX] = "ux";
  mMode2Suffix[UY] = "uy";
  mMode2Suffix[UZ] = "uz";
  mMode2Suffix[RHO] = "rho";
  mMode2Suffix[LAMBDA] = "lambda";
  mMode2Suffix[MU] = "mu";
  mMode2Suffix[UXEXACT] = "uxexact";
  mMode2Suffix[UYEXACT] = "uyexact";
  mMode2Suffix[UZEXACT] = "uzexact";
  mMode2Suffix[P] = "p";
  mMode2Suffix[S] = "s";
  mMode2Suffix[DIV] = "div";
  mMode2Suffix[CURLMAG] = "curlmag";
  mMode2Suffix[DIVDT] = "divdt";
  mMode2Suffix[CURLMAGDT] = "curlmagdt";
  mMode2Suffix[LAT] = "lat";
  mMode2Suffix[LON] = "lon";
  mMode2Suffix[HMAXDUDT] = "hmaxdudt";
  mMode2Suffix[VMAXDUDT] = "vmaxdudt";
  mMode2Suffix[TOPO] = "topo";
  mMode2Suffix[GRIDX] = "gridx";
  mMode2Suffix[GRIDY] = "gridy";
  mMode2Suffix[GRIDZ] = "gridz";
  mMode2Suffix[UXERR] = "uxerr";
  mMode2Suffix[UYERR] = "uyerr";
  mMode2Suffix[UZERR] = "uzerr";
  mMode2Suffix[MAGDUDT] = "magdudt";
  mMode2Suffix[HMAGDUDT] = "hmagdudt";
  mMode2Suffix[MAG] = "mag";
  mMode2Suffix[HMAG] = "hmag";
  mMode2Suffix[HMAX] = "hmax";
  mMode2Suffix[VMAX] = "vmax";
  mMode2Suffix[GRADRHO] = "gradrho";
  mMode2Suffix[GRADMU] = "gradmu";
  mMode2Suffix[GRADLAMBDA] = "gradlambda";
  mMode2Suffix[GRADP] = "gradp";
  mMode2Suffix[GRADS] = "grads";
  mMode2Suffix[QP] = "qp";
  mMode2Suffix[QS] = "qs";

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
void Image::associate_gridfiles( vector<Image*>& imgs )
{
   // Only if a curvilinear grid is needed
   if( mEW->topographyExists() && (mLocationType == Image::X || mLocationType == Image::Y) &&
       !(mMode == Image::GRIDX || mMode == Image::GRIDY || mMode == Image::GRIDZ) )
   {
   // if my location type is X find Gridz with equal x-coord.
   //                        Y find Gridz with equal y-coord.
   // Require same precision for now, simplifies file format.
      for( int i = 0 ; i < imgs.size() ; i++ )
      {
         if( mLocationType == Image::X && imgs[i]->mMode == Image::GRIDZ && imgs[i]->mLocationType == Image::X &&
				       fabs(mCoordValue-imgs[i]->mCoordValue) < 1e-3 && m_double == imgs[i]->m_double )
	    m_gridimage = imgs[i];
	 if( mLocationType == Image::Y && imgs[i]->mMode == Image::GRIDZ && imgs[i]->mLocationType == Image::Y &&
	     fabs(mCoordValue-imgs[i]->mCoordValue) < 1e-3 && m_double == imgs[i]->m_double )
	    m_gridimage = imgs[i];
      }
      if( m_gridimage == Image::nil )
      {
         m_gridimage = new Image( mEW, 0.0, 0.0, 0, 0, mFilePrefix, Image::GRIDZ, mLocationType, mCoordValue,
				  m_double, m_usehdf5, false );
	 m_gridimage->computeGridPtIndex();
	 m_gridimage->allocatePlane();
         imgs.push_back(m_gridimage);
      }
      if( m_gridimage == Image::nil )
      {
	 mGridinfo = -1; // Failed to find grid
         cout << "WARNING: Image::associate_gridfiles did not find compatible grid images" << endl;
      }
      else
	 mGridinfo = mStoreGrid ? 1 : 2 ; // Found grid_images
   }
   else
      mGridinfo = 0; // Grid is Cartesian, or image is a grid.
}

//-----------------------------------------------------------------------
//bool Image::proc_write()
//{
//  int myRank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//  return (myRank == m_rankWriter);
//}

//-------------------------------------
void Image::setSteps(int a_steps)
{
  char buffer[50];
  mPreceedZeros = snprintf( buffer, 50, "%d", a_steps );
}

//-----------------------------------------------------------------------
void Image::computeGridPtIndex()
{
// the purpose of this routine is to assign the vector<int> m_gridPtIndex
/* Here, we compute the index --in the local grid-- of the coordinate value at which we have to plot. 
   For x, y, the local grid is the same as the global grid, but for z, k resets at each refinement boundary. */

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

  //  m_rankWriter = fileWriterIDs[0];

  if (m_mpiComm_writers != MPI_COMM_NULL) 
      MPI_Comm_free(&m_mpiComm_writers);

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
  //  m_gridPtValueInitialized = true;
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
void Image::initializeTime(double t)
{
  mNextTime = t; 
  m_time_done = false;
// with the option timeInterval=..., first time is always t=0
}

//-----------------------------------------------------------------------
const std::string Image::fieldSuffix(ImageMode mode) const
{
  REQUIRE2(mode >= 0 && mode < MODES, "mode=" << mode << " out of bounds");
  return mMode2Suffix[mode];
}

//-----------------------------------------------------------------------
bool Image::is_time_derivative() const
{
   return mMode == Image::DIVDT    || mMode == Image::CURLMAGDT || 
          mMode == Image::MAGDUDT  || mMode == Image::HMAGDUDT  || 
          mMode == Image::HMAXDUDT || mMode == Image::VMAXDUDT;
}

//-----------------------------------------------------------------------
bool Image::timeToWrite(float_sw4 time, int cycle, float_sw4 dt )
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
bool Image::timeToWrite( int cycle )
{
  // -----------------------------------------------
  // Check based on cycle only
  // -----------------------------------------------
  bool do_it=false;
  if( cycle == mWritingCycle )
    do_it = true;
  if( mCycleInterval !=  0 && cycle%mCycleInterval == 0 ) 
    do_it = true;
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
      // mpiComm_writers consists of all processors that
      // own some part of the image.
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
void Image::computeImageDivCurl( vector<Sarray> &a_Up, vector<Sarray>& a_U,
				 vector<Sarray> &a_Um, float_sw4 dt, int dminus )
{
   ASSERT(m_isDefinedMPIWriters);
   ASSERT( mMode == Image::DIV || mMode == Image::CURLMAG || mMode == Image::DIVDT
	   || mMode == Image::CURLMAGDT );
// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
	 gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      if( mMode == Image::DIV )
      {
	 vector<float_sw4*> div(mEW->mNumberOfGrids);
	 for( int g=gmin ; g <= gmax ; g++ )
	 {
	    int npts=(mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1)*(mWindow[g][5]-mWindow[g][4]+1);
	    div[g] = new float_sw4[npts];
	 }
	 computeDivergence( a_Up, div );
	 for( int g=gmin ; g <= gmax ; g++ )
	 {
	    int npts=(mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1)*(mWindow[g][5]-mWindow[g][4]+1);
	    if( m_double )
	       for( size_t i=0 ; i < npts ; i++ )
		  m_doubleField[g][i] = (double)div[g][i];
	    else
	       for( size_t i=0 ; i < npts ; i++ )
		  m_floatField[g][i] = (float) div[g][i];
	    delete[] div[g];
	 }
      }
      else if( mMode == Image::CURLMAG )
      {
	 vector<float_sw4*> curl(mEW->mNumberOfGrids);
	 for( int g=gmin ; g <= gmax ; g++ )
	 {
	    int npts=(mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1)*(mWindow[g][5]-mWindow[g][4]+1);
	    curl[g] = new float_sw4[3*npts];
	 }
	 computeCurl( a_Up, curl );
	 for( int g=gmin ; g <= gmax ; g++ )
	 {
	    int npts=(mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1)*(mWindow[g][5]-mWindow[g][4]+1);
	    if( m_double )
	       for( size_t i=0 ; i < npts ; i++ )
		  m_doubleField[g][i] = (double)sqrt(curl[g][3*i]*curl[g][3*i]+curl[g][3*i+1]*curl[g][3*i+1]+
						curl[g][3*i+2]*curl[g][3*i+2]);
	    else
	       for( size_t i=0 ; i < npts ; i++ )
		  m_floatField[g][i] = (float) sqrt(curl[g][3*i]*curl[g][3*i]+curl[g][3*i+1]*curl[g][3*i+1]+
						       curl[g][3*i+2]*curl[g][3*i+2]);
	    delete[] curl[g];
	 }
      }
      else if( mMode == Image::DIVDT )
      {
	 vector<float_sw4*> div(mEW->mNumberOfGrids), divm(mEW->mNumberOfGrids);
	 for( int g=gmin ; g <= gmax ; g++ )
	 {
	    int npts=(mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1)*(mWindow[g][5]-mWindow[g][4]+1);
	    div[g]  = new float_sw4[npts];
	    divm[g] = new float_sw4[npts];
	 }
	 computeDivergence( a_Up, div );
	 float_sw4 idt;
	 if( dminus )
	 {
	    computeDivergence(a_U,divm);
	    idt = 1/dt;
	 }
	 else
	 {
	    computeDivergence(a_Um,divm);
	    idt = 1/(2*dt);
	 }
	 for( int g=gmin ; g <= gmax ; g++ )
	 {
	    int npts=(mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1)*(mWindow[g][5]-mWindow[g][4]+1);
	    if( m_double )
	       for( size_t i=0 ; i < npts ; i++ )
		  m_doubleField[g][i] = (double) idt*(div[g][i]-divm[g][i]);
	    else
	       for( size_t i=0 ; i < npts ; i++ )
		  m_floatField[g][i] = (float) idt*(div[g][i]-divm[g][i]);
	    delete[] div[g];
	    delete[] divm[g];
	 } 
      }
      else if( mMode == Image::CURLMAGDT )
      {
	 vector<float_sw4*> curl(mEW->mNumberOfGrids), curlm(mEW->mNumberOfGrids);
	 for( int g=gmin ; g <= gmax ; g++ )
	 {
	    int npts=(mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1)*(mWindow[g][5]-mWindow[g][4]+1);
	    curl[g]  = new float_sw4[3*npts];
	    curlm[g] = new float_sw4[3*npts];
	 }
	 computeCurl( a_Up, curl );
	 float_sw4 idt;
	 if( dminus )
	 {
	    computeCurl( a_U, curlm );
	    idt = 1/dt;
	 }
	 else
	 {
	    computeCurl( a_Um, curlm );
	    idt = 1/(2*dt);
	 }
	 for( int g=gmin ; g <= gmax ; g++ )
	 {
	    int npts=(mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1)*(mWindow[g][5]-mWindow[g][4]+1);
	    if( m_double )
	       for( size_t i=0 ; i < npts ; i++ )
		  m_doubleField[g][i] = (double) idt*sqrt( (curl[g][3*i]-curlm[g][3*i])*(curl[g][3*i]-curlm[g][3*i]) +
						  (curl[g][3*i+1]-curlm[g][3*i+1])*(curl[g][3*i+1]-curlm[g][3*i+1]) +
						  (curl[g][3*i+2]-curlm[g][3*i+2])*(curl[g][3*i+2]-curlm[g][3*i+2]));
	    else
	       for( size_t i=0 ; i < npts ; i++ )
		  m_floatField[g][i] = (float) idt*sqrt( (curl[g][3*i]-curlm[g][3*i])*(curl[g][3*i]-curlm[g][3*i]) +
							 (curl[g][3*i+1]-curlm[g][3*i+1])*(curl[g][3*i+1]-curlm[g][3*i+1]) +
							 (curl[g][3*i+2]-curlm[g][3*i+2])*(curl[g][3*i+2]-curlm[g][3*i+2]));
	    delete[] curl[g];
	    delete[] curlm[g];
	 }
      }
   }
}

//-----------------------------------------------------------------------
void Image::computeImageQuantity(std::vector<Sarray> &a_mu, int a_nComp )
{
  ASSERT(m_isDefinedMPIWriters);
// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
  bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
  if (iwrite)
  {
    int gmin, gmax;
    if( mLocationType == Image::Z )
       gmin = gmax = m_gridPtIndex[1];
    else
    {
       gmin = 0;
       gmax = mEW->mNumberOfGrids-1;
    }    
    for( int g=gmin ; g <= gmax ; g++ )
    {
       //       size_t iField=0;
       size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
       size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
       for( int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
	  for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	     for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	     {
		size_t iField = (ii-mWindow[g][0])+ni*(jj-mWindow[g][2])+nij*(kk-mWindow[g][4]);
		if( m_double )
		   m_doubleField[g][iField] = (double) a_mu[g](a_nComp,ii,jj,kk);
		else
		   m_floatField[g][iField] = (float) a_mu[g](a_nComp,ii,jj,kk);
		//		iField++;      
	     }
    }
  }
} // end Image::computeImageQuantity()

//-----------------------------------------------------------------------
void Image::computeImageGrid( Sarray &a_X, Sarray &a_Y, Sarray &a_Z )
{
  ASSERT(m_isDefinedMPIWriters);
  ASSERT( mMode == Image::GRIDX ||mMode == Image::GRIDY ||mMode == Image::GRIDZ );
  if (plane_in_proc(m_gridPtIndex[0]))
  {
     int topCartesian = mEW->mNumberOfCartesianGrids - 1;
     int component = 1;
     if( mMode == Image::GRIDY )
	component = 2;
     else if( mMode == Image::GRIDZ )
	component = 3;

     int gmin, gmax;
     if( mLocationType == Image::Z )
	gmin = gmax = m_gridPtIndex[1];
     else
     {
	gmin = 0;
	gmax = mEW->mNumberOfGrids-1;
     }
     for (int g = gmin; g <= gmax ; g++ )
     {
//	size_t iField = 0;
       size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
       size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
	for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
	   for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	      for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	      {                   
		 size_t iField = (ii-mWindow[g][0])+ni*(jj-mWindow[g][2])+nij*(kk-mWindow[g][4]);
		 float_sw4 val;
		 if (g > topCartesian)
		 {
		    if(component==1 )
		       val = a_X(ii,jj,kk);
		    else if (component == 2)
		       val = a_Y(ii,jj,kk);
		    else if (component == 3)
		       val = a_Z(ii,jj,kk);
		 }
		 else
		 {
		    if (component==1)
		       val = (ii-1)*mEW->mGridSize[g];
		    else if (component == 2)
		       val = (jj-1)*mEW->mGridSize[g];
		    else if (component == 3)
		       val = (kk-1)*mEW->mGridSize[g] + mEW->m_zmin[g];
		 }
		 if( m_double )
		    m_doubleField[g][iField] = (double) val;
		 else
		    m_floatField[g][iField] = (float) val;
		 //		 iField++;      
	      }
     }
  }
}

//-----------------------------------------------------------------------
void Image::computeImageLatLon(Sarray &a_X, Sarray &a_Y, Sarray &a_Z ) 
{
   ASSERT(m_isDefinedMPIWriters);
   ASSERT( mMode == Image::LAT || mMode == Image::LON );
// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
   if (plane_in_proc(m_gridPtIndex[0]))
   {

      int g;
// lat, lon images are only meaningful on Image::Z planes, but we will write lat and lon on any plane
      int topCartesian = mEW->mNumberOfCartesianGrids - 1;
// write the data...
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
	 gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      for (g = gmin; g <= gmax ; g++ )
      {
//	 size_t iField = 0;
       size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
       size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
	 for (int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
#pragma omp parallel for
	    for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	       for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	       {
		 size_t iField = (ii-mWindow[g][0])+ni*(jj-mWindow[g][2])+nij*(kk-mWindow[g][4]);
		 double latP, lonP, xP, yP, zP, val;
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
		  mEW->computeGeographicCoord(xP, yP, lonP, latP);
		  if( mMode == Image::LAT )
		     val = latP;
		  else 
		     val = lonP;
		  if( m_double )
		     m_doubleField[g][iField] = val;
		  else
		     m_floatField[g][iField] = (float) val;
		  //		  iField++;      
	       } // end for ii	  
      }
   }
}


//-----------------------------------------------------------------------
void Image::computeImagePvel(std::vector<Sarray> &mu, std::vector<Sarray> &lambda,
			     std::vector<Sarray> &rho )
{
   ASSERT(m_isDefinedMPIWriters);
   ASSERT( mMode == Image::P );
   if (plane_in_proc(m_gridPtIndex[0]))
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
	 gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }    
      for( int g=gmin ; g <= gmax ; g++ )
      {
	 //	 size_t iField=0;
       size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
       size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
	 for( int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
	    for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	       for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	       {
		  size_t iField = (ii-mWindow[g][0])+ni*(jj-mWindow[g][2])+nij*(kk-mWindow[g][4]);
		  if( m_double )
		     m_doubleField[g][iField] = (double) sqrt( (2*mu[g](ii,jj,kk)+lambda[g](ii,jj,kk))/rho[g](ii,jj,kk));
		  else
		     m_floatField[g][iField] = (float) sqrt( (2*mu[g](ii,jj,kk)+lambda[g](ii,jj,kk))/rho[g](ii,jj,kk));
		  //		  iField++;      
	       }
      }
   }
}

//-----------------------------------------------------------------------
void Image::computeImageSvel(std::vector<Sarray> &mu, std::vector<Sarray> &rho )
{
   ASSERT(m_isDefinedMPIWriters);
   ASSERT( mMode == Image::S );
   if (plane_in_proc(m_gridPtIndex[0]))
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
	 gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }    
      for( int g=gmin ; g <= gmax ; g++ )
      {
	 //	 size_t iField=0;
       size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
       size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
	 for( int kk = mWindow[g][4]; kk <= mWindow[g][5]; kk++)
	    for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	       for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	       {
		  size_t iField = (ii-mWindow[g][0])+ni*(jj-mWindow[g][2])+nij*(kk-mWindow[g][4]);
		  if( m_double )
		     m_doubleField[g][iField] = (double) sqrt( mu[g](ii,jj,kk)/rho[g](ii,jj,kk));
		  else
		     m_floatField[g][iField] = (float) sqrt( mu[g](ii,jj,kk)/rho[g](ii,jj,kk));
		  //		  iField++;
	       }
      }
   }
}

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
// write the data...
      int g = m_gridPtIndex[1];
      //      size_t iField = 0;
       size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
#pragma omp parallel for
      for (int jj = mWindow[g][2]; jj <= mWindow[g][3]; jj++)
	for (int ii = mWindow[g][0]; ii <= mWindow[g][1]; ii++)
	{                   
	   size_t iField = (ii-mWindow[g][0])+ni*(jj-mWindow[g][2]);
	  if( m_double )
	     m_doubleField[g][iField] = (double) u2(ii,jj,1);
	  else
	    m_floatField[g][iField] = (float) u2(ii,jj,1);
 //	  iField++;      
	}
   }
}

//-----------------------------------------------------------------------
void Image::writeImagePlane_2(int cycle, std::string &path, float_sw4 t )
{
   if( !m_user_created )
      return;
   
   ASSERT(m_isDefinedMPIWriters);

   double stime, etime;
   stime = MPI_Wtime();
#ifdef USE_HDF5
   hid_t h5_fid, grd, dset, attr, dtype, attr_space1, dset_space, fapl;
#endif
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
  
  // Header: [precision(int) npatches(int) time, plane, coordinate, imagetype, gridinfo, timeofday,
  //                            h_1 zmin_1 sizes_1 h_2 zmin_2 sizes_2 ... h_ng zmin_ng sizes_ng ]
  // plane - 0 x=const
   //        1 y=const
   //        2 z=const
   // coordinate: value of constant variable, e.g., x=10.4 if plane=0.
   // imagetype - mode nr. (ux=1, uy=2, etc...)
   // gridinfo - -1 No grid info given.
   //             0 All patches are on Cartesian grids, no curvilinear grid needed.
   //             1 Curvilinear grid is stored after the data patches.
   //             2 Grid file names are stored after the data patches.
   // timeofday - String describing the creation date of the file, max 25 bytes.

  // offset initialized to header size:
   off_t offset = 2*sizeof(int) + 2*sizeof(double) + 3*sizeof(int) + 25*sizeof(char) +
      (ghigh-glow)*(2*sizeof(double)+4*sizeof(int));
   int fid=-1;
   stringstream s, fileSuffix;
   int prec, nPatches, globalPlaneSize[4];
   if( iwrite )
   {
      if( mMode == Image::GRIDX || mMode == Image::GRIDY || mMode == Image::GRIDZ )
	 compute_file_suffix( fileSuffix, 0 );
      else
	 compute_file_suffix( fileSuffix, cycle );

      if( path != "." )
	 s << path;
      s << fileSuffix.str();
   }

   if( m_usehdf5 == false && m_pio[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
      if (fid == -1 )
      {
	 VERIFY2(0, "ERROR: Image::writeImagePlane_2, error opening file " << s.str() << " for writing header");
      }  

      cout << "writing image plane on file " << s.str() << endl;// " (msg from proc # " << m_rankWriter << ")" << endl;
      prec = m_double ? 8 : 4;
      size_t ret=write(fid,&prec,sizeof(int));
      if( ret != sizeof(int) )
	 cout << "ERROR: Image::writeImagePlane_2 could not write precision" << endl;

      nPatches = mLocationType == Image::Z ? 1 : mEW->mNumberOfGrids;
      ret = write(fid,&nPatches,sizeof(int));
      if( ret != sizeof(int) )
	 cout << "ERROR: Image::writeImagePlane_2 could not write number of patches" << endl;

      double dblevar=static_cast<double>(t);
      ret = write(fid,&dblevar,sizeof(double));
      if( ret != sizeof(double) )
	 cout << "ERROR: Image::writeImagePlane_2 could not write time" << endl;

      int ltype;
      if( mLocationType == Image::X )
         ltype = 0;
      else if( mLocationType == Image::Y )
	 ltype = 1;
      else
	 ltype = 2;

      ret = write(fid,&ltype,sizeof(int));
      if( ret != sizeof(int) )
	 cout << "ERROR: Image::writeImagePlane_2 could not write plane type" << endl;

      float_sw4 coordvalue = (m_gridPtIndex[0]-1)*mEW->mGridSize[glow];
      if( mLocationType == Image::Z )
         coordvalue += mEW->m_zmin[glow];

      dblevar=static_cast<double>(coordvalue);
      ret = write(fid,&dblevar,sizeof(double));
      if( ret != sizeof(double) )
	 cout << "ERROR: Image::writeImagePlane_2 could not write coordinate value" << endl;

      int imode = static_cast<int>(mMode);
      ret = write(fid,&imode,sizeof(int));
      if( ret != sizeof(int) )
	 cout << "ERROR: Image::writeImagePlane_2 could not write imode" << endl;

      ret = write(fid,&mGridinfo,sizeof(int));
      if( ret != sizeof(int) )
	 cout << "ERROR: Image::writeImagePlane_2 could not write gridinfo" << endl;

      time_t realtime;
      time(&realtime);
      string strtime;
      strtime += asctime(localtime(&realtime));
      char strtimec[25];
      strncpy(strtimec,strtime.c_str(),25);
      strtimec[24] ='\0';
      ret = write(fid,strtimec,25*sizeof(char));
      if( ret != 25*sizeof(char) )
	 cout << "ERROR: Image::writeImagePlane_2 could not write strtimec" << endl;
      
      for(int g = glow; g < ghigh ;g++)
      {
	 dblevar = static_cast<double>(mEW->mGridSize[g]);
	 ret = write(fid,&dblevar,sizeof(double));
         if( ret != sizeof(double) )
	    cout << "ERROR: Image::writeImagePlane_2 could not write h for grid " << g << endl;

	 dblevar = static_cast<double>(mEW->m_zmin[g]);
	 ret = write(fid,&dblevar,sizeof(double));
         if( ret != sizeof(double) )
	    cout << "ERROR: Image::writeImagePlane_2 could not write zmin for grid " << g << endl;

	   // should hold the global number of interior points
	 if(mLocationType == Image::X)
	 {
	    globalPlaneSize[0] = 1;
	    globalPlaneSize[1] = mEW->m_global_ny[g];
	    globalPlaneSize[2] = 1;
	    globalPlaneSize[3] = mEW->m_global_nz[g];
	 }
	 if(mLocationType == Image::Y)
	 {
	    globalPlaneSize[0] = 1;
	    globalPlaneSize[1] = mEW->m_global_nx[g];
	    globalPlaneSize[2] = 1;
	    globalPlaneSize[3] = mEW->m_global_nz[g];
	 }
	 if(mLocationType == Image::Z)
	 {
	    globalPlaneSize[0] = 1;
	    globalPlaneSize[1] = mEW->m_global_nx[g];
	    globalPlaneSize[2] = 1;
	    globalPlaneSize[3] = mEW->m_global_ny[g];
	 }
	 ret = write(fid,globalPlaneSize,4*sizeof(int));
         if( ret != 4*sizeof(int) )
	    cout << "ERROR: Image::writeImagePlane_2 could not write dimensions of grid " << g << endl;
      }
      fsync(fid);
   }
   else if( m_usehdf5 == true && m_pio[0]->proc_zero() ) 
   {
#ifdef USE_HDF5
      int ret, ltype;
      hsize_t dims, dims1 = 1, total_elem = 0;
      setenv("HDF5_USE_FILE_LOCKING", "FALSE", 1);
      h5_fid = H5Fcreate((const char*)(s.str().c_str()), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (h5_fid < 0) 
	VERIFY2(0, "ERROR: Image::writeImagePlane_2, error opening HDF5 file " << s.str() << " for writing header");

      /* cout << "Rank " << mEW->getRank() << " created new file [" << s.str() << "]" << endl; */
      cout << "writing image plane on file " << s.str() << endl;// " (msg from proc # " << m_rankWriter << ")" << endl;

      attr_space1 = H5Screate_simple(1, &dims1, NULL);

      nPatches = mLocationType == Image::Z ? 1 : mEW->mNumberOfGrids;
      ret = createWriteAttr(h5_fid, "npatch", H5T_NATIVE_INT, attr_space1, &nPatches);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write number of patches" << endl;

      double dblevar=static_cast<double>(t);
      ret = createWriteAttr(h5_fid, "time", H5T_NATIVE_DOUBLE, attr_space1, &dblevar);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 time" << endl;

      if( mLocationType == Image::X )
         ltype = 0;
      else if( mLocationType == Image::Y )
	 ltype = 1;
      else
	 ltype = 2;

      ret = createWriteAttr(h5_fid, "plane", H5T_NATIVE_INT, attr_space1, &ltype);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 plane type" << endl;

      float_sw4 coordvalue = (m_gridPtIndex[0]-1)*mEW->mGridSize[glow];
      if( mLocationType == Image::Z )
         coordvalue += mEW->m_zmin[glow];

      dblevar=static_cast<double>(coordvalue);
      ret = createWriteAttr(h5_fid, "coordinate", H5T_NATIVE_DOUBLE, attr_space1, &dblevar);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 coordinate value" << endl;

      int imode = static_cast<int>(mMode);
      ret = createWriteAttr(h5_fid, "mode", H5T_NATIVE_INT, attr_space1, &imode);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 imode" << endl;

      ret = createWriteAttr(h5_fid, "gridinfo", H5T_NATIVE_INT, attr_space1, &mGridinfo);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 gridinfo" << endl;

      H5Sclose(attr_space1);

      time_t realtime;
      time(&realtime);
      string strtime;
      strtime += asctime(localtime(&realtime));
      char strtimec[25];
      strncpy(strtimec,strtime.c_str(),25);
      strtimec[24] ='\0';

      createWriteAttrStr(h5_fid, "creationtime", (const char*)strtimec);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 strtimec" << endl;

      double *grid_size = new double[nPatches];
      double *zmin      = new double[nPatches];
      int    *ni        = new int[nPatches];
      int    *nj        = new int[nPatches];
      /* int    *ib        = new int[nPatches]; */
      /* int    *jb        = new int[nPatches]; */
      
      for(int g = glow; g < ghigh ;g++)
      {
          grid_size[g-glow] = static_cast<double>(mEW->mGridSize[g]);
          zmin[g-glow]      = static_cast<double>(mEW->m_zmin[g]); 

         // should hold the global number of interior points
	 if(mLocationType == Image::X)
	 {
	    ni[g-glow] = mEW->m_global_ny[g];
	    nj[g-glow] = mEW->m_global_nz[g];
	 }
	 if(mLocationType == Image::Y)
	 {
	    ni[g-glow] = mEW->m_global_nx[g];
	    nj[g-glow] = mEW->m_global_nz[g];
	 }
	 if(mLocationType == Image::Z)
	 {
	    ni[g-glow] = mEW->m_global_nx[g];
	    nj[g-glow] = mEW->m_global_ny[g];
	 }
         total_elem += (ni[g-glow]*nj[g-glow]);
      }

      dims = nPatches;
      dset_space = H5Screate_simple(1, &dims, NULL);

      dset = H5Dcreate(h5_fid, "grid_size", H5T_NATIVE_DOUBLE, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ret  = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid_size);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 grid_size dataset" << endl;
      H5Dclose(dset);

      dset = H5Dcreate(h5_fid, "zmin", H5T_NATIVE_DOUBLE, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ret  = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, zmin);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 zmin dataset" << endl;
      H5Dclose(dset);


      dset = H5Dcreate(h5_fid, "ni", H5T_NATIVE_INT, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ret  = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ni);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 ni dataset" << endl;
      H5Dclose(dset);


      dset = H5Dcreate(h5_fid, "nj", H5T_NATIVE_INT, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ret  = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nj);
      if( ret < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not write HDF5 nj dataset" << endl;
      H5Dclose(dset);

      H5Sclose(dset_space);

      if (m_double) 
        dtype = H5T_NATIVE_DOUBLE;
      else
        dtype = H5T_NATIVE_FLOAT;

      dims = total_elem;
      dset_space = H5Screate_simple(1, &dims, NULL);
      /* cout << "Rank " << mEW->getRank() << " creating patches array with length " << dims << endl; */
      dset = H5Dcreate(h5_fid, "patches", dtype, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if( dset < 0 )
	 cout << "ERROR: Image::writeImagePlane_2 could not create HDF5 patches dataset" << endl;
      H5Sclose(dset_space);
      H5Dclose(dset);

      if( mGridinfo == 1 )
      {
        int g=mEW->mNumberOfGrids-1;
        int globalSizes[3] = {mEW->m_global_nx[g], mEW->m_global_ny[g], mEW->m_global_nz[g]} ;
        if(mLocationType == Image::X)
  	 globalSizes[0]    = 1;
        if (mLocationType == Image::Y)
  	 globalSizes[1]    = 1;
        if (mLocationType == Image::Z)
  	 globalSizes[2]    = 1;

         dims = globalSizes[0]*globalSizes[1]*globalSizes[2];
         dset_space = H5Screate_simple(1, &dims, NULL);
         /* cout << "Rank " << mEW->getRank() << " creating grid array with length " << dims << endl; */
         dset = H5Dcreate(h5_fid, "grid", dtype, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         if( dset < 0 )
            cout << "ERROR: Image::writeImagePlane_2 could not create HDF5 grid dataset" << endl;
         H5Sclose(dset_space);
         H5Dclose(dset);
      }

      delete [] grid_size;
      delete [] zmin;
      delete [] ni;
      delete [] nj;

      /* H5Fflush(h5_fid, H5F_SCOPE_LOCAL); */
      H5Fclose(h5_fid);
#else
     cout << "ERROR: cannot write image in HDF5 format without sw4 compiled with HDF5 library!" << endl;
#endif
   } // End proc 0 write metadata

// write the data...
   if(m_usehdf5 == false) 
   {
     if( ihavearray )
     {
        for( int g = glow ; g < ghigh; g++)
        {
           int globalSizes[3] = {mEW->m_global_nx[g],
           		       mEW->m_global_ny[g],
           		       mEW->m_global_nz[g]} ;
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
                 VERIFY2(0, "ERROR: Image::writeImagePlane2, Error opening file: " << s.str() << " for writing data patches" );
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
               //	    bool dbg=false;
               //	    if( mMode == S && mEW->getRank() == 10 )
               //	    {
               //	       cout << "Before write, value is  " << m_floatField[g][618] << " g= " << g << "off= " << offset << " glow= " << glow << " ghigh= "<< ghigh << endl;
               //	       dbg = true;
               //	    }
               //	    if( mMode == S && mEW->getRank() == 63 && mLocationType == Z )
               //	       dbg = true;
              m_pio[g-glow]->write_array( &fid, 1, m_floatField[g], offset, fltStr );
              offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]*sizeof(float));
           }
        }
     }
     if( iwrite )
        close(fid);

     if( mGridinfo == 1 )
        add_grid_to_file( s.str().c_str(), iwrite, offset );
     if( mGridinfo == 2 )
        add_grid_filenames_to_file( s.str().c_str() );
   }
   else
   {
#ifdef USE_HDF5
      hsize_t dims, dims1 = 1;
      if( ihavearray )
      {
         offset = 0;
         for( int g = glow ; g < ghigh; g++)
         {
            int globalSizes[3] = {mEW->m_global_nx[g], mEW->m_global_ny[g], mEW->m_global_nz[g]} ;

            if(mLocationType == Image::X)
               globalSizes[0]    = 1;
            if (mLocationType == Image::Y)
               globalSizes[1]    = 1;
            if (mLocationType == Image::Z)
               globalSizes[2]    = 1;

            /* if( !mEW->usingParallelFS() || g == glow ) */
            if( g == glow )
            {
               MPI_Barrier( m_mpiComm_writers );
               /* cout << "Rank " << mEW->getRank() <<" after barrier" << endl; */
            }
            if( m_double )
            {
               char dblStr[]="double";	  
               m_pio[g-glow]->write_array_hdf5( (const char*)(s.str().c_str()), "patches", 1, m_doubleField[g], offset, dblStr );
               /* m_pio[g-glow]->write_array_hdf5( dset, 1, m_doubleField[g], offset, dblStr ); */
               offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]);
            }
            else
            {
               char fltStr[]="float";
               m_pio[g-glow]->write_array_hdf5( (const char*)(s.str().c_str()), "patches", 1, m_floatField[g], offset, fltStr );
               offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]);
            }
         }
      } // End if ihavearray

      if( mGridinfo == 1 )
         add_grid_to_file_hdf5( s.str().c_str(), iwrite, 0);
      /* if( mGridinfo == 2 ) */
      /*    add_grid_filenames_to_file( s.str().c_str() ); */
#endif
   } // End use_hdf5 == true
   etime = MPI_Wtime();
   m_write_time += (etime - stime);
}

//-----------------------------------------------------------------------
void Image::add_grid_to_file_hdf5( const char* fname, bool iwrite, size_t offset )
{
   bool ihavearray = plane_in_proc(m_gridPtIndex[0]);
   if( ihavearray )
   {
#ifdef USE_HDF5
      int g=mEW->mNumberOfGrids-1;
      int globalSizes[3] = {mEW->m_global_nx[g], mEW->m_global_ny[g], mEW->m_global_nz[g]} ;
      if(mLocationType == Image::X)
	 globalSizes[0]    = 1;
      if (mLocationType == Image::Y)
	 globalSizes[1]    = 1;
      if (mLocationType == Image::Z)
	 globalSizes[2]    = 1;

      if( mEW->usingParallelFS() )
	 MPI_Barrier( m_mpiComm_writers );

      if( m_double )
      {
	 char dblStr[]="double";	  
	 m_pio[g]->write_array_hdf5(fname, "grid", 1, m_gridimage->m_doubleField[g], offset, dblStr );
	 /* m_pio[g]->write_array_hdf5(dset, 1, m_gridimage->m_doubleField[g], offset, dblStr ); */
	 offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]*sizeof(double));
      }
      else
      {
	 char fltStr[]="float";
	 m_pio[g]->write_array_hdf5(fname, "grid", 1, m_gridimage->m_floatField[g], offset, fltStr );
	 /* m_pio[g]->write_array_hdf5(dset, 1, m_gridimage->m_floatField[g], offset, fltStr ); */
	 offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]*sizeof(float));
      }
#endif
   }
}

//-----------------------------------------------------------------------
void Image::add_grid_to_file( const char* fname, bool iwrite, size_t offset )
{
   bool ihavearray = plane_in_proc(m_gridPtIndex[0]);
   if( ihavearray )
   {
      int fid;
      if( iwrite )
      {
	 fid = open( fname, O_WRONLY, 0660 ); 
	 if (fid == -1 )
	    VERIFY2(0, "ERROR: Image::add_grid_to_file, error opening file " << fname );
      }

      int g=mEW->mNumberOfGrids-1;
      int globalSizes[3] = {mEW->m_global_nx[g],
			    mEW->m_global_ny[g],
			    mEW->m_global_nz[g]} ;
      if(mLocationType == Image::X)
	 globalSizes[0]    = 1;
      if (mLocationType == Image::Y)
	 globalSizes[1]    = 1;
      if (mLocationType == Image::Z)
	 globalSizes[2]    = 1;

      if( !mEW->usingParallelFS() )
	 MPI_Barrier( m_mpiComm_writers );

      if( m_double )
      {
	 char dblStr[]="double";	  
	 m_pio[g]->write_array( &fid, 1, m_gridimage->m_doubleField[g], offset, dblStr );
	 offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]*sizeof(double));
      }
      else
      {
	 char fltStr[]="float";
	 m_pio[g]->write_array( &fid, 1, m_gridimage->m_floatField[g], offset, fltStr );
	 offset += (globalSizes[0]*globalSizes[1]*globalSizes[2]*sizeof(float));
      }
      if( iwrite )
	 close(fid);
   }
}

//-----------------------------------------------------------------------
void Image::add_grid_filenames_to_file( const char* fname )
{
   if( m_pio[0]->proc_zero() )
   {
      stringstream str1;
      m_gridimage->compute_file_suffix( str1, 0 );
      string img1str = str1.str();
      int fid = open( fname, O_WRONLY, 0660 ); 
      if (fid == -1 )
	 VERIFY2(0, "ERROR: Image::add_grid_filenames_to_file, error opening file " << fname );
      size_t nr=lseek(fid,0,SEEK_END);
      int n = img1str.length();
      nr = write(fid,&n,sizeof(int) );
      if( nr != sizeof(int) )
	 cout << "ERROR: Image::add_grid_filenames_to file, could not write n1 " << endl;

      nr = write(fid,str1.str().c_str(),sizeof(char)*n);
      if( nr != n*sizeof(char) )
	 cout << "ERROR: Image::add_grid_filenames_to file, could not write name1 " << endl;
      close(fid);
   }
}

//-----------------------------------------------------------------------
void Image::compute_file_suffix( stringstream& fileSuffix, int cycle )
{
   fileSuffix << mFilePrefix << ".cycle=";
   int temp = static_cast<int>(pow(10.0, mPreceedZeros - 1));
   int testcycle = cycle;
   if (cycle == 0)
      testcycle=1;
   
   while (testcycle < temp)
   {
      fileSuffix << "0";
      temp /= 10;
   }
   fileSuffix << cycle << "." << mOrientationString[mLocationType] << "=";
   fileSuffix << mCoordValue;
   fileSuffix << "." << fieldSuffix(mMode) << ".sw4img";
   // Add .h5 suffix
   if (m_usehdf5) 
      fileSuffix << ".h5";
}

//-----------------------------------------------------------------------
void Image::update_maxes_hVelMax( vector<Sarray> &a_Up, vector<Sarray> &a_Um,
				  float_sw4 dt )
{
   static bool firstHVM = true;
   bool geographic = false;
   // geographic = true  --> Max of N- and E-components (max(|u_{N-S}|,|u_{E-W}|))
   // geograpic  = false --> L2 norm of x- and y-components (sqrt(ux*ux+uy*uy)).
   int gmin, gmax;
   if( mLocationType == Image::Z )
      gmin = gmax = m_gridPtIndex[1];
   else
   {
      gmin = 0;
      gmax = mEW->mNumberOfGrids-1;
   }
   float_sw4 di  = 1/dt  ;
   if( geographic )
   {
//C.B> with N-S and E-W max: max(|u_{N-S}|,|u_{E-W}|))
      for( int g = gmin; g <= gmax; g++ )
      {
	 //	 size_t iField = 0;
	 //	 double velNS, velEW;
	 //	 double Ux, Uy, nrm;
	 //	 double xP, yP, latitude;
	 double latOrigin=mEW->getLatOrigin();
	 double deg2rad = M_PI/180.0;//, cphi, sphi, calpha, salpha, thxnrm, thynrm;
	 double az = deg2rad*mEW->getGridAzimuth(); // Note that mGeoAz is in degrees
	 double metersPerDeg = mEW->getMetersPerDegree();
	 double calpha = cos(az);
	 double salpha = sin(az);
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
	 for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
#pragma omp parallel for
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
 // first get velocities in the (x,y) directions
		  float_sw4 Ux = (a_Up[g](1,i,j,k) - a_Um[g](1,i,j,k))*di;
		  float_sw4 Uy = (a_Up[g](2,i,j,k) - a_Um[g](2,i,j,k))*di;

 // note: the angle between (x,y) and East is constant throughout the mesh, but the angle to North 
 // varies throughout the grid. Here we use the chain rule on the mapping between (x,y) and (lon,lat)
 // see EW:computeGeographicCoord for the mapping
               
		  double xP = (i-1)*mEW->mGridSize[g];
		  double yP = (j-1)*mEW->mGridSize[g];

		  double latitude = latOrigin + (xP*calpha - yP*salpha)/metersPerDeg;
	  
		  double cphi = cos(latitude*deg2rad);
		  double sphi = sin(latitude*deg2rad);

		  double thxnrm = salpha + (xP*salpha+yP*calpha)/cphi/metersPerDeg * (M_PI/180.0) * sphi * calpha;
		  double thynrm = calpha - (xP*salpha+yP*calpha)/cphi/metersPerDeg * (M_PI/180.0) * sphi * salpha;

		  double nrm = sqrt( thxnrm*thxnrm + thynrm*thynrm );
		  thxnrm /= nrm;
		  thynrm /= nrm;

		  float_sw4 velNS = fabs( thynrm*Ux  - thxnrm*Uy );
		  float_sw4 velEW = fabs( salpha*Ux  + calpha*Uy );
		  if (m_double)
		  {
		     if( firstHVM ||  m_doubleField[g][iField] < velNS)
			m_doubleField[g][iField] = (double) velNS;
		     if( m_doubleField[g][iField] < velEW )
			m_doubleField[g][iField] = (double) velEW;
		  }
		  else
		  {
		     if( firstHVM ||  m_floatField[g][iField] < (float) velNS)
			m_floatField[g][iField] = (float) velNS;
		     if(m_floatField[g][iField] < (float) velEW)
			m_floatField[g][iField] = (float) velEW;
		  }
		  //		  iField++;  
	       }
      }
   }
   else
   {
      for( int g = gmin; g <= gmax; g++ )
      {
	 //	 size_t iField = 0;
	 //	 double Ux, Uy, nrm, hvmag;
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
	 for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
#pragma omp parallel for
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		  float_sw4 Ux = (a_Up[g](1,i,j,k) - a_Um[g](1,i,j,k))*di;
		  float_sw4 Uy = (a_Up[g](2,i,j,k) - a_Um[g](2,i,j,k))*di;
                  float_sw4 hvmag = sqrt(Ux*Ux+Uy*Uy);
		  if (m_double)
		  {
		     if( firstHVM ||  m_doubleField[g][iField] < hvmag )
			m_doubleField[g][iField] = (double) hvmag;
		  }
		  else
		  {
		     if( firstHVM ||  m_floatField[g][iField] < (float) hvmag)
			m_floatField[g][iField] = (float) hvmag;
		  }
//		  iField++;  
	       }
      }
   }
   if( firstHVM )
      firstHVM = false;
}

//-----------------------------------------------------------------------
void Image::update_maxes_vVelMax(std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, float_sw4 dt )
{
   static bool firstVVM = true;
   int gmax, gmin;
   if( mLocationType == Image::Z )
      gmax = gmin = m_gridPtIndex[1];
   else
   {
      gmin = 0;
      gmax = mEW->mNumberOfGrids-1;
   }
   float_sw4 di  = 1/dt;
   for( int g = gmin; g <= gmax; g++ )
   {
	 //      size_t iField = 0;
      size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
      size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
      for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
#pragma omp parallel for
	 for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	    for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	    {
	       float_sw4 vel = fabs( a_Up[g](3,i,j,k)-a_Um[g](3,i,j,k) )*di;
	       size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
	       if (m_double)
	       {
		  if( firstVVM ||  m_doubleField[g][iField] < vel )
		     m_doubleField[g][iField] = (double) vel;
	       }
	       else
	       {
		  if( firstVVM ||  m_floatField[g][iField] < (float) vel )
		     m_floatField[g][iField] = (float) vel;
	       }
   //	       iField++; 
	    }
   }
   if( firstVVM )
      firstVVM = false;
}

//-----------------------------------------------------------------------
void Image::update_maxes_hMax( vector<Sarray> &a_U )
			
{
   static bool firstHM = true;
   int gmin, gmax;
   if( mLocationType == Image::Z )
      gmin = gmax = m_gridPtIndex[1];
   else
   {
      gmin = 0;
      gmax = mEW->mNumberOfGrids-1;
   }
   for( int g = gmin; g <= gmax; g++ )
   {
      //      size_t iField = 0;
      //      float_sw4 hmag;
      size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
      size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
      for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
#pragma omp parallel for
	 for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	    for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	    {
	       size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
	       float_sw4 hmag = sqrt(a_U[g](1,i,j,k)*a_U[g](1,i,j,k) + a_U[g](2,i,j,k)*a_U[g](2,i,j,k));
	       if (m_double)
	       {
		  if( firstHM ||  m_doubleField[g][iField] < hmag )
		     m_doubleField[g][iField] = (double) hmag;
	       }
	       else
	       {
		  if( firstHM ||  m_floatField[g][iField] < (float) hmag)
		     m_floatField[g][iField] = (float) hmag;
	       }
	       //	       iField++;  
	    }
   }
   if( firstHM )
      firstHM = false;
}

//-----------------------------------------------------------------------
void Image::update_maxes_vMax(std::vector<Sarray> &a_U )
{
   static bool firstVM = true;
   int gmax, gmin;
   if( mLocationType == Image::Z )
      gmax = gmin = m_gridPtIndex[1];
   else
   {
      gmin = 0;
      gmax = mEW->mNumberOfGrids-1;
   }
   for( int g = gmin; g <= gmax; g++ )
   {
      //      size_t iField = 0;
      size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
      size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
      for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
#pragma omp parallel for
	 for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	    for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	    {
	       float_sw4 uz = fabs(a_U[g](3,i,j,k));
	       size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
	       if (m_double)
	       {
		  if( firstVM ||  m_doubleField[g][iField] < uz )
		     m_doubleField[g][iField] = (double) uz;
	       }
	       else
	       {
		  if( firstVM ||  m_floatField[g][iField] < (float) uz )
		     m_floatField[g][iField] = (float) uz;
	       }
	       //	       iField++; 
	    }
   }
   if( firstVM )
      firstVM = false;
}

//-----------------------------------------------------------------------
void Image::computeDivergence( std::vector<Sarray> &a_U, std::vector<float_sw4*>& a_div )
{
   ASSERT(m_isDefinedMPIWriters);
   ASSERT( a_U.size() == mEW->mNumberOfGrids )
   ASSERT( a_div.size() == a_U.size() )

// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
   bool iwrite    = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      const float_sw4 c1 = 2.0/3;
      const float_sw4 c2 =-1.0/12;
      const float_sw4 fs = 5.0/6;
      const float_sw4 ot = 1.0/12;
      const float_sw4 ft = 4.0/3;
      const float_sw4 os = 1.0/6;
      const float_sw4 d3 = 14.0/3;

      int gmin,gmax;
      bool dotop;
      if( mLocationType == Image::Z )
      {
         if( m_gridPtIndex[1] < mEW->mNumberOfCartesianGrids )
	 {	    
// Z-plane in Cartesian mesh
	    gmin = gmax = m_gridPtIndex[1];
            dotop = false;
	 }
	 else
	 {
// Z-plane in curvilinear mesh, set gmin,gmax so that 
// the Cartesian loop is not executed
	    gmin = 0;
            gmax =-1;
	    dotop = true;
	 }
      }
      else
      {
// X- or Y- plane. All grids contribute.
         gmin = 0;
	 gmax = mEW->mNumberOfCartesianGrids - 1;
	 dotop = true;
      }

      int gh = mEW->getNumberOfGhostPoints();
// Do the Cartesian grids.
      for (int g = gmin; g <= gmax ; g++ )
      {
          float_sw4 factor = 1.0/(mEW->mGridSize[g]);
	  //          size_t iField = 0;
	  size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	  size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
          for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
            for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
              for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
              {
		 float_sw4 val =0.;
		 size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		 // I-direction
		 if( i == 1 )
		 {
                    val =-2.25*a_U[g](1,i,j,k)+(4+fs)*a_U[g](1,i+1,j,k)-d3*a_U[g](1,i+2,j,k)+
		       3*a_U[g](1,i+3,j,k)-(1+ot)*a_U[g](1,i+4,j,k) +os*a_U[g](1,i+5,j,k);
		 }
                 else if( i == 2 )
		 {
                    val  = -os*a_U[g](1,i-1,j,k) -1.25*a_U[g](1,i,j,k)+(1+ft)*a_U[g](1,i+1,j,k)
		       - ft*a_U[g](1,i+2,j,k) + 0.5*a_U[g](1,i+3,j,k) - ot*a_U[g](1,i+4,j,k);
		 }
                 else if (i == mEW->m_global_nx[g] )
		 {
                    val = 2.25*a_U[g](1,i,j,k)-(4+fs)*a_U[g](1,i-1,j,k)+d3*a_U[g](1,i-2,j,k)-
		       3*a_U[g](1,i-3,j,k)+(1+ot)*a_U[g](1,i-4,j,k) -os*a_U[g](1,i-5,j,k);
		 }
                 else if (i == mEW->m_global_nx[g]-1 )
		 {
                    val  = os*a_U[g](1,i+1,j,k) +1.25*a_U[g](1,i,j,k)-(1+ft)*a_U[g](1,i-1,j,k) +
		        ft*a_U[g](1,i-2,j,k) - 0.5*a_U[g](1,i-3,j,k) + ot*a_U[g](1,i-4,j,k);
		 }
		 else
		 {
		    val  =  c1*(a_U[g](1,i+1,j,k)-a_U[g](1,i-1,j,k)) + c2*(a_U[g](1,i+2,j,k)-a_U[g](1,i-2,j,k));
		 }

		 // J-direction
		 if( j == 1 )
		 {
                    val +=-2.25*a_U[g](2,i,j,k)+(4+fs)*a_U[g](2,i,j+1,k)-d3*a_U[g](2,i,j+2,k)+
		       3*a_U[g](2,i,j+3,k)-(1+ot)*a_U[g](2,i,j+4,k) +os*a_U[g](2,i,j+5,k);
		 }
                 else if( j == 2 )
		 {
                    val  += -os*a_U[g](2,i,j-1,k) -1.25*a_U[g](2,i,j,k)+(1+ft)*a_U[g](2,i,j+1,k)
		       - ft*a_U[g](2,i,j+2,k) + 0.5*a_U[g](2,i,j+3,k) - ot*a_U[g](2,i,j+4,k);
		 }
                 else if (j == mEW->m_global_ny[g] )
		 {
                    val += 2.25*a_U[g](2,i,j,k)-(4+fs)*a_U[g](2,i,j-1,k)+d3*a_U[g](2,i,j-2,k)-
		       3*a_U[g](2,i,j-3,k)+(1+ot)*a_U[g](2,i,j-4,k) -os*a_U[g](2,i,j-5,k);
		 }
                 else if (j == mEW->m_global_ny[g]-1 )
		 {
                    val  += os*a_U[g](2,i,j+1,k) +1.25*a_U[g](2,i,j,k)-(1+ft)*a_U[g](2,i,j-1,k) +
		        ft*a_U[g](2,i,j-2,k) - 0.5*a_U[g](2,i,j-3,k) + ot*a_U[g](2,i,j-4,k);
		 }
		 else
		 {
		    val  +=  c1*(a_U[g](2,i,j+1,k)-a_U[g](2,i,j-1,k)) + c2*(a_U[g](2,i,j+2,k)-a_U[g](2,i,j-2,k));
		 }

		 // K-direction
                 if( k == 1 && g==mEW->mNumberOfGrids-1 )
                 {
                    val +=-2.25*a_U[g](3,i,j,k)+(4+fs)*a_U[g](3,i,j,k+1)-d3*a_U[g](3,i,j,k+2)+
		       3*a_U[g](3,i,j,k+3)-(1+ot)*a_U[g](3,i,j,k+4) +os*a_U[g](3,i,j,k+5);
		 }
		 else if( k == 2 && g == mEW->mNumberOfGrids-1 )
		 {
                    val  += -os*a_U[g](3,i,j,k-1) -1.25*a_U[g](3,i,j,k)+(1+ft)*a_U[g](3,i,j,k+1)
		       - ft*a_U[g](3,i,j,k+2) + 0.5*a_U[g](3,i,j,k+3) - ot*a_U[g](3,i,j,k+4);
		 }
                 else if( k == mEW->m_kEnd[g]-gh && g == 0 )
		 {
                    val += 2.25*a_U[g](3,i,j,k)-(4+fs)*a_U[g](3,i,j,k-1)+d3*a_U[g](3,i,j,k-2)-
		       3*a_U[g](3,i,j,k-3)+(1+ot)*a_U[g](3,i,j,k-4) -os*a_U[g](3,i,j,k-5);
		 }
                 else if( k == mEW->m_kEnd[g]-gh-1 && g == 0 )
		 {
                    val  += os*a_U[g](3,i,j,k+1) +1.25*a_U[g](3,i,j,k)-(1+ft)*a_U[g](3,i,j,k-1) +
		        ft*a_U[g](3,i,j,k-2) - 0.5*a_U[g](3,i,j,k-3) + ot*a_U[g](3,i,j,k-4);
		 }
                 else
		 {
		    val  +=  c1*(a_U[g](3,i,j,k+1)-a_U[g](3,i,j,k-1)) + c2*(a_U[g](3,i,j,k+2)-a_U[g](3,i,j,k-2));
		 }

                 val *= factor;
                 a_div[g][iField] = val;
		 //		 iField++; 
	      }
      }
       // curvilinear grid
      if ( mEW->topographyExists() && dotop )
      {
          int g = mEW->mNumberOfGrids - 1;
          float_sw4 factor = 1.0/2.0; 
	  //          size_t iField = 0;
	  size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	  size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
          for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
#pragma omp parallel for
            for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
              for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
                {
                  float_sw4 val = 0.;
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
                  if (i == 1)
		  {
		     val  = mEW->mMetric(1,i,j,k)*(
                               -2.25*a_U[g](1,i,j,k)+(4+fs)*a_U[g](1,i+1,j,k)-d3*a_U[g](1,i+2,j,k)+
			       3*a_U[g](1,i+3,j,k)-(1+ot)*a_U[g](1,i+4,j,k) +os*a_U[g](1,i+5,j,k) );
		  }
                  else if( i== 2 )
		  {
		     val  = mEW->mMetric(1,i,j,k)*(
				 -os*a_U[g](1,i-1,j,k) -1.25*a_U[g](1,i,j,k)+(1+ft)*a_U[g](1,i+1,j,k)
				- ft*a_U[g](1,i+2,j,k) + 0.5*a_U[g](1,i+3,j,k) - ot*a_U[g](1,i+4,j,k) );
		  }
                  else if (i == mEW->m_global_nx[g]-1)
		  {
		     val  = mEW->mMetric(1,i,j,k)*(os*a_U[g](1,i+1,j,k) +1.25*a_U[g](1,i,j,k)-(1+ft)*a_U[g](1,i-1,j,k) +
						   ft*a_U[g](1,i-2,j,k) - 0.5*a_U[g](1,i-3,j,k) + ot*a_U[g](1,i-4,j,k));
		  }
                  else if (i == mEW->m_global_nx[g])
		  {
		     val  = mEW->mMetric(1,i,j,k)*(2.25*a_U[g](1,i,j,k)-(4+fs)*a_U[g](1,i-1,j,k)+d3*a_U[g](1,i-2,j,k)-
						   3*a_U[g](1,i-3,j,k)+(1+ot)*a_U[g](1,i-4,j,k) -os*a_U[g](1,i-5,j,k) );
		  }
                  else
		  {
		       val  = mEW->mMetric(1,i,j,k)*(c1*(a_U[g](1,i+1,j,k) - a_U[g](1,i-1,j,k)) +
						     c2*(a_U[g](1,i+2,j,k) - a_U[g](1,i-2,j,k) )  );
		  }
                  if (j == 1)
                    {
		       val  += mEW->mMetric(1,i,j,k)*(-2.25*a_U[g](2,i,j,k)+(4+fs)*a_U[g](2,i,j+1,k)-d3*a_U[g](2,i,j+2,k)+
						      3*a_U[g](2,i,j+3,k)-(1+ot)*a_U[g](2,i,j+4,k) +os*a_U[g](2,i,j+5,k));
                    }
		  else if( j== 2 )
                    {
		       val  += mEW->mMetric(1,i,j,k)*(-os*a_U[g](2,i,j-1,k) -1.25*a_U[g](2,i,j,k)+(1+ft)*a_U[g](2,i,j+1,k)
						      - ft*a_U[g](2,i,j+2,k) + 0.5*a_U[g](2,i,j+3,k) - ot*a_U[g](2,i,j+4,k));
                    }
                  else if (j == mEW->m_global_ny[g]-1)
                    {
		       val  += mEW->mMetric(1,i,j,k)*( os*a_U[g](2,i,j+1,k) +1.25*a_U[g](2,i,j,k)-(1+ft)*a_U[g](2,i,j-1,k) +
		        ft*a_U[g](2,i,j-2,k) - 0.5*a_U[g](2,i,j-3,k) + ot*a_U[g](2,i,j-4,k));
                    }
                  else if (j == mEW->m_global_ny[g])
                    {
		       val  += mEW->mMetric(1,i,j,k)*(2.25*a_U[g](2,i,j,k)-(4+fs)*a_U[g](2,i,j-1,k)+d3*a_U[g](2,i,j-2,k)-
						      3*a_U[g](2,i,j-3,k)+(1+ot)*a_U[g](2,i,j-4,k) -os*a_U[g](2,i,j-5,k) );
                    }
                  else
                    {
		       val  += mEW->mMetric(1,i,j,k)*(c1*(a_U[g](2,i,j+1,k) - a_U[g](2,i,j-1,k)) +
						      c2*(a_U[g](2,i,j+2,k) - a_U[g](2,i,j-2,k) ));
                    }
                  
                  if (k == 1)
                    {
		       val  += mEW->mMetric(2,i,j,k)*(-2.25*a_U[g](1,i,j,k)+(4+fs)*a_U[g](1,i,j,k+1)-d3*a_U[g](1,i,j,k+2)+
		       3*a_U[g](1,i,j,k+3)-(1+ot)*a_U[g](1,i,j,k+4) +os*a_U[g](1,i,j,k+5))
			      +mEW->mMetric(3,i,j,k)*(-2.25*a_U[g](2,i,j,k)+(4+fs)*a_U[g](2,i,j,k+1)-d3*a_U[g](2,i,j,k+2)+
		       3*a_U[g](2,i,j,k+3)-(1+ot)*a_U[g](2,i,j,k+4) +os*a_U[g](2,i,j,k+5))
			      +mEW->mMetric(4,i,j,k)*(-2.25*a_U[g](3,i,j,k)+(4+fs)*a_U[g](3,i,j,k+1)-d3*a_U[g](3,i,j,k+2)+
							3*a_U[g](3,i,j,k+3)-(1+ot)*a_U[g](3,i,j,k+4) +os*a_U[g](3,i,j,k+5));
                    }
                  if (k == 2)
                    {
                      val  += mEW->mMetric(2,i,j,k)*(-os*a_U[g](1,i,j,k-1) -1.25*a_U[g](1,i,j,k)+(1+ft)*a_U[g](1,i,j,k+1)
		       - ft*a_U[g](1,i,j,k+2) + 0.5*a_U[g](1,i,j,k+3) - ot*a_U[g](1,i,j,k+4))
                              +mEW->mMetric(3,i,j,k)*(-os*a_U[g](2,i,j,k-1) -1.25*a_U[g](2,i,j,k)+(1+ft)*a_U[g](2,i,j,k+1)
		       - ft*a_U[g](2,i,j,k+2) + 0.5*a_U[g](2,i,j,k+3) - ot*a_U[g](2,i,j,k+4))
                              +mEW->mMetric(4,i,j,k)*(-os*a_U[g](3,i,j,k-1) -1.25*a_U[g](3,i,j,k)+(1+ft)*a_U[g](3,i,j,k+1)
		       - ft*a_U[g](3,i,j,k+2) + 0.5*a_U[g](3,i,j,k+3) - ot*a_U[g](3,i,j,k+4));
                    }
                  else // no k=N because we are in the curvilinear grid
                    {
		       val  += mEW->mMetric(2,i,j,k)*(c1*(a_U[g](1,i,j,k+1) - a_U[g](1,i,j,k-1)) +
						      c2*(a_U[g](1,i,j,k+2) - a_U[g](1,i,j,k-2)) )
			      +mEW->mMetric(3,i,j,k)*(c1*(a_U[g](2,i,j,k+1) - a_U[g](2,i,j,k-1)) +
						      c2*(a_U[g](2,i,j,k+2) - a_U[g](2,i,j,k-2)) )
			      +mEW->mMetric(4,i,j,k)*(c1*(a_U[g](3,i,j,k+1) - a_U[g](3,i,j,k-1)) + 
						      c2*(a_U[g](3,i,j,k+2) - a_U[g](3,i,j,k-2))   );
                    }
                  val *= factor/sqrt(mEW->mJ(i,j,k));
		  a_div[g][iField] = val;
		  //                  if( m_double )
		  //		     m_doubleField[g][iField] = (double) val;
		  //                  else
		  //                      m_floatField[g][iField] = (float)val;  
		  //                  iField++;  
                }   
      }
   }
}

//-----------------------------------------------------------------------
void Image::computeCurl( std::vector<Sarray>& a_U, std::vector<float_sw4*>& a_curl )
{
   ASSERT(m_isDefinedMPIWriters);
   ASSERT( a_U.size()    == mEW->mNumberOfGrids )
   ASSERT( a_curl.size() == a_U.size() )

// plane_in_proc returns true for z=const lpanes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      const float_sw4 c1 = 2.0/3;
      const float_sw4 c2 =-1.0/12;
      const float_sw4 fs = 5.0/6;
      const float_sw4 ot = 1.0/12;
      const float_sw4 ft = 4.0/3;
      const float_sw4 os = 1.0/6;
      const float_sw4 d3 = 14.0/3;

      int gmin,gmax;
      bool dotop;
      if( mLocationType == Image::Z )
      {
         if( m_gridPtIndex[1] < mEW->mNumberOfCartesianGrids )
	 {	    
// Z-plane in Cartesian mesh
	    gmin = gmax = m_gridPtIndex[1];
            dotop = false;
	 }
	 else
	 {
// Z-plane in curvilinear mesh, set gmin,gmax so that 
// the Cartesian loop is not executed
	    gmin = 0;
            gmax =-1;
	    dotop = true;
	 }
      }
      else
      {
// X- or Y- plane. All grids contribute.
         gmin = 0;
	 gmax = mEW->mNumberOfCartesianGrids - 1;
	 dotop = true;
      }
      int gh = mEW->getNumberOfGhostPoints();
// Do the Cartesian grids.
      for (int g = gmin; g <= gmax ; g++ )
      {
          float_sw4 factor = 1.0/(mEW->mGridSize[g]);
  //          size_t iField = 0;
	  size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	  size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
          for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
            for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
              for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
              {
		 float_sw4 duydx, duzdx, duxdy, duzdy, duxdz, duydz;
		 size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		 // I-direction
		 if( i == 1 )
		 {
                    duydx =-2.25*a_U[g](2,i,j,k)+(4+fs)*a_U[g](2,i+1,j,k)-d3*a_U[g](2,i+2,j,k)+
		       3*a_U[g](2,i+3,j,k)-(1+ot)*a_U[g](2,i+4,j,k) +os*a_U[g](2,i+5,j,k);
                    duzdx =-2.25*a_U[g](3,i,j,k)+(4+fs)*a_U[g](3,i+1,j,k)-d3*a_U[g](3,i+2,j,k)+
		       3*a_U[g](3,i+3,j,k)-(1+ot)*a_U[g](3,i+4,j,k) +os*a_U[g](3,i+5,j,k);
		 }
                 else if( i == 2 )
		 {
                    duydx  = -os*a_U[g](2,i-1,j,k) -1.25*a_U[g](2,i,j,k)+(1+ft)*a_U[g](2,i+1,j,k)
		       - ft*a_U[g](2,i+2,j,k) + 0.5*a_U[g](2,i+3,j,k) - ot*a_U[g](2,i+4,j,k);
                    duzdx  = -os*a_U[g](3,i-1,j,k) -1.25*a_U[g](3,i,j,k)+(1+ft)*a_U[g](3,i+1,j,k)
		       - ft*a_U[g](3,i+2,j,k) + 0.5*a_U[g](3,i+3,j,k) - ot*a_U[g](3,i+4,j,k);
		 }
                 else if (i == mEW->m_global_nx[g] )
		 {
                    duydx = 2.25*a_U[g](2,i,j,k)-(4+fs)*a_U[g](2,i-1,j,k)+d3*a_U[g](2,i-2,j,k)-
		       3*a_U[g](2,i-3,j,k)+(1+ot)*a_U[g](2,i-4,j,k) -os*a_U[g](2,i-5,j,k);
                    duzdx = 2.25*a_U[g](3,i,j,k)-(4+fs)*a_U[g](3,i-1,j,k)+d3*a_U[g](3,i-2,j,k)-
		       3*a_U[g](3,i-3,j,k)+(1+ot)*a_U[g](3,i-4,j,k) -os*a_U[g](3,i-5,j,k);
		 }
                 else if (i == mEW->m_global_nx[g]-1 )
		 {
                    duydx  = os*a_U[g](2,i+1,j,k) +1.25*a_U[g](2,i,j,k)-(1+ft)*a_U[g](2,i-1,j,k) +
		        ft*a_U[g](2,i-2,j,k) - 0.5*a_U[g](2,i-3,j,k) + ot*a_U[g](2,i-4,j,k);
                    duzdx  = os*a_U[g](3,i+1,j,k) +1.25*a_U[g](3,i,j,k)-(1+ft)*a_U[g](3,i-1,j,k) +
		        ft*a_U[g](3,i-2,j,k) - 0.5*a_U[g](3,i-3,j,k) + ot*a_U[g](3,i-4,j,k);
		 }
		 else
		 {
		    duydx  =  c1*(a_U[g](2,i+1,j,k)-a_U[g](2,i-1,j,k)) + c2*(a_U[g](2,i+2,j,k)-a_U[g](2,i-2,j,k));
		    duzdx  =  c1*(a_U[g](3,i+1,j,k)-a_U[g](3,i-1,j,k)) + c2*(a_U[g](3,i+2,j,k)-a_U[g](3,i-2,j,k));
		 }

		 // J-direction
		 if( j == 1 )
		 {
                    duxdy =-2.25*a_U[g](1,i,j,k)+(4+fs)*a_U[g](1,i,j+1,k)-d3*a_U[g](1,i,j+2,k)+
		       3*a_U[g](1,i,j+3,k)-(1+ot)*a_U[g](1,i,j+4,k) +os*a_U[g](1,i,j+5,k);
                    duzdy =-2.25*a_U[g](3,i,j,k)+(4+fs)*a_U[g](3,i,j+1,k)-d3*a_U[g](3,i,j+2,k)+
		       3*a_U[g](3,i,j+3,k)-(1+ot)*a_U[g](3,i,j+4,k) +os*a_U[g](3,i,j+5,k);
		 }
                 else if( j == 2 )
		 {
                    duxdy = -os*a_U[g](1,i,j-1,k) -1.25*a_U[g](1,i,j,k)+(1+ft)*a_U[g](1,i,j+1,k)
		       - ft*a_U[g](1,i,j+2,k) + 0.5*a_U[g](1,i,j+3,k) - ot*a_U[g](1,i,j+4,k);
                    duzdy = -os*a_U[g](3,i,j-1,k) -1.25*a_U[g](3,i,j,k)+(1+ft)*a_U[g](3,i,j+1,k)
		       - ft*a_U[g](3,i,j+2,k) + 0.5*a_U[g](3,i,j+3,k) - ot*a_U[g](3,i,j+4,k);
		 }
                 else if (j == mEW->m_global_ny[g] )
		 {
                    duxdy = 2.25*a_U[g](1,i,j,k)-(4+fs)*a_U[g](1,i,j-1,k)+d3*a_U[g](1,i,j-2,k)-
		       3*a_U[g](1,i,j-3,k)+(1+ot)*a_U[g](1,i,j-4,k) -os*a_U[g](1,i,j-5,k);
                    duzdy = 2.25*a_U[g](3,i,j,k)-(4+fs)*a_U[g](3,i,j-1,k)+d3*a_U[g](3,i,j-2,k)-
		       3*a_U[g](3,i,j-3,k)+(1+ot)*a_U[g](3,i,j-4,k) -os*a_U[g](3,i,j-5,k);
		 }
                 else if (j == mEW->m_global_ny[g]-1 )
		 {
                    duxdy = os*a_U[g](1,i,j+1,k) +1.25*a_U[g](1,i,j,k)-(1+ft)*a_U[g](1,i,j-1,k) +
		        ft*a_U[g](1,i,j-2,k) - 0.5*a_U[g](1,i,j-3,k) + ot*a_U[g](1,i,j-4,k);
                    duzdy = os*a_U[g](3,i,j+1,k) +1.25*a_U[g](3,i,j,k)-(1+ft)*a_U[g](3,i,j-1,k) +
		        ft*a_U[g](3,i,j-2,k) - 0.5*a_U[g](3,i,j-3,k) + ot*a_U[g](3,i,j-4,k);
		 }
		 else
		 {
		    duxdy = c1*(a_U[g](1,i,j+1,k)-a_U[g](1,i,j-1,k)) + c2*(a_U[g](1,i,j+2,k)-a_U[g](1,i,j-2,k));
		    duzdy = c1*(a_U[g](3,i,j+1,k)-a_U[g](3,i,j-1,k)) + c2*(a_U[g](3,i,j+2,k)-a_U[g](3,i,j-2,k));
		 }

		 // K-direction
                 if( k == 1 && g==mEW->mNumberOfGrids-1 )
                 {
                    duxdz =-2.25*a_U[g](1,i,j,k)+(4+fs)*a_U[g](1,i,j,k+1)-d3*a_U[g](1,i,j,k+2)+
		       3*a_U[g](1,i,j,k+3)-(1+ot)*a_U[g](1,i,j,k+4) +os*a_U[g](1,i,j,k+5);
                    duydz =-2.25*a_U[g](2,i,j,k)+(4+fs)*a_U[g](2,i,j,k+1)-d3*a_U[g](2,i,j,k+2)+
		       3*a_U[g](2,i,j,k+3)-(1+ot)*a_U[g](2,i,j,k+4) +os*a_U[g](2,i,j,k+5);
		 }
		 else if( k == 2 && g == mEW->mNumberOfGrids-1 )
		 {
                    duxdz = -os*a_U[g](1,i,j,k-1) -1.25*a_U[g](1,i,j,k)+(1+ft)*a_U[g](1,i,j,k+1)
		       - ft*a_U[g](1,i,j,k+2) + 0.5*a_U[g](1,i,j,k+3) - ot*a_U[g](1,i,j,k+4);
                    duydz = -os*a_U[g](2,i,j,k-1) -1.25*a_U[g](2,i,j,k)+(1+ft)*a_U[g](2,i,j,k+1)
		       - ft*a_U[g](2,i,j,k+2) + 0.5*a_U[g](2,i,j,k+3) - ot*a_U[g](2,i,j,k+4);
		 }
                 else if( k == mEW->m_kEnd[g]-gh && g == 0 )
		 {
                    duxdz = 2.25*a_U[g](1,i,j,k)-(4+fs)*a_U[g](1,i,j,k-1)+d3*a_U[g](1,i,j,k-2)-
		       3*a_U[g](1,i,j,k-3)+(1+ot)*a_U[g](1,i,j,k-4) -os*a_U[g](1,i,j,k-5);
                    duydz = 2.25*a_U[g](2,i,j,k)-(4+fs)*a_U[g](2,i,j,k-1)+d3*a_U[g](2,i,j,k-2)-
		       3*a_U[g](2,i,j,k-3)+(1+ot)*a_U[g](2,i,j,k-4) -os*a_U[g](2,i,j,k-5);
		 }
                 else if( k == mEW->m_kEnd[g]-gh-1 && g == 0 )
		 {
                    duxdz = os*a_U[g](1,i,j,k+1) +1.25*a_U[g](1,i,j,k)-(1+ft)*a_U[g](1,i,j,k-1) +
		        ft*a_U[g](1,i,j,k-2) - 0.5*a_U[g](1,i,j,k-3) + ot*a_U[g](1,i,j,k-4);
                    duydz = os*a_U[g](2,i,j,k+1) +1.25*a_U[g](2,i,j,k)-(1+ft)*a_U[g](2,i,j,k-1) +
		        ft*a_U[g](2,i,j,k-2) - 0.5*a_U[g](2,i,j,k-3) + ot*a_U[g](2,i,j,k-4);
		 }
                 else
		 {
		    duxdz = c1*(a_U[g](1,i,j,k+1)-a_U[g](1,i,j,k-1)) + c2*(a_U[g](1,i,j,k+2)-a_U[g](1,i,j,k-2));
		    duydz = c1*(a_U[g](2,i,j,k+1)-a_U[g](2,i,j,k-1)) + c2*(a_U[g](2,i,j,k+2)-a_U[g](2,i,j,k-2));
		 }
                 a_curl[g][3*iField]   = factor*(duzdy-duydz); 
                 a_curl[g][3*iField+1] = factor*(duxdz-duzdx);
                 a_curl[g][3*iField+2] = factor*(duydx-duxdy);
		 //		 iField++; 
	      }
      }
       // curvilinear grid
      if ( mEW->topographyExists() && dotop )
      {
          int g = mEW->mNumberOfGrids - 1;
	  //          int iField = 0;
	  //          float_sw4 factor = 1.0/(2);
	  size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	  size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
          for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
#pragma omp parallel for
	     for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
		for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
                {
	           size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		   float_sw4 duxdq = 0., duydq=0.0, duzdq=0.0;
		   float_sw4 duxdr = 0., duydr=0.0, duzdr=0.0;
		   float_sw4 duxds = 0., duyds=0.0, duzds=0.0;
		   if (i == 1)
		   {      
                      duxdq = -2.25*a_U[g](1,i,j,k)+(4+fs)*a_U[g](1,i+1,j,k)-d3*a_U[g](1,i+2,j,k)+
		       3*a_U[g](1,i+3,j,k)-(1+ot)*a_U[g](1,i+4,j,k) +os*a_U[g](1,i+5,j,k);
                      duydq = -2.25*a_U[g](2,i,j,k)+(4+fs)*a_U[g](2,i+1,j,k)-d3*a_U[g](2,i+2,j,k)+
		       3*a_U[g](2,i+3,j,k)-(1+ot)*a_U[g](2,i+4,j,k) +os*a_U[g](2,i+5,j,k);
                      duzdq = -2.25*a_U[g](3,i,j,k)+(4+fs)*a_U[g](3,i+1,j,k)-d3*a_U[g](3,i+2,j,k)+
		       3*a_U[g](3,i+3,j,k)-(1+ot)*a_U[g](3,i+4,j,k) +os*a_U[g](3,i+5,j,k);
		   }
		   if (i == 2)
		   {      
		      duxdq  = -os*a_U[g](1,i-1,j,k) -1.25*a_U[g](1,i,j,k)+(1+ft)*a_U[g](1,i+1,j,k)
		       - ft*a_U[g](1,i+2,j,k) + 0.5*a_U[g](1,i+3,j,k) - ot*a_U[g](1,i+4,j,k);
		      duydq  = -os*a_U[g](2,i-1,j,k) -1.25*a_U[g](2,i,j,k)+(1+ft)*a_U[g](2,i+1,j,k)
		       - ft*a_U[g](2,i+2,j,k) + 0.5*a_U[g](2,i+3,j,k) - ot*a_U[g](2,i+4,j,k);
		      duzdq  = -os*a_U[g](3,i-1,j,k) -1.25*a_U[g](3,i,j,k)+(1+ft)*a_U[g](3,i+1,j,k)
		       - ft*a_U[g](3,i+2,j,k) + 0.5*a_U[g](3,i+3,j,k) - ot*a_U[g](3,i+4,j,k);
		   }
		   else if (i == mEW->m_global_nx[g]-1)
		   { 
		      duxdq  = os*a_U[g](1,i+1,j,k) +1.25*a_U[g](1,i,j,k)-(1+ft)*a_U[g](1,i-1,j,k) +
		        ft*a_U[g](1,i-2,j,k) - 0.5*a_U[g](1,i-3,j,k) + ot*a_U[g](1,i-4,j,k);
		      duydq  = os*a_U[g](2,i+1,j,k) +1.25*a_U[g](2,i,j,k)-(1+ft)*a_U[g](2,i-1,j,k) +
		        ft*a_U[g](2,i-2,j,k) - 0.5*a_U[g](2,i-3,j,k) + ot*a_U[g](2,i-4,j,k);
		      duzdq  = os*a_U[g](3,i+1,j,k) +1.25*a_U[g](3,i,j,k)-(1+ft)*a_U[g](3,i-1,j,k) +
		        ft*a_U[g](3,i-2,j,k) - 0.5*a_U[g](3,i-3,j,k) + ot*a_U[g](3,i-4,j,k);
		   }
		   else if (i == mEW->m_global_nx[g])
		   { 
		      duxdq = 2.25*a_U[g](1,i,j,k)-(4+fs)*a_U[g](1,i-1,j,k)+d3*a_U[g](1,i-2,j,k)-
		       3*a_U[g](1,i-3,j,k)+(1+ot)*a_U[g](1,i-4,j,k) -os*a_U[g](1,i-5,j,k);
		      duydq = 2.25*a_U[g](2,i,j,k)-(4+fs)*a_U[g](2,i-1,j,k)+d3*a_U[g](2,i-2,j,k)-
		       3*a_U[g](2,i-3,j,k)+(1+ot)*a_U[g](2,i-4,j,k) -os*a_U[g](2,i-5,j,k);
		      duzdq = 2.25*a_U[g](3,i,j,k)-(4+fs)*a_U[g](3,i-1,j,k)+d3*a_U[g](3,i-2,j,k)-
		       3*a_U[g](3,i-3,j,k)+(1+ot)*a_U[g](3,i-4,j,k) -os*a_U[g](3,i-5,j,k);
		   }
		   else
		   { 
                      duxdq = c1*(a_U[g](1,i+1,j,k) - a_U[g](1,i-1,j,k))+c2*(a_U[g](1,i+2,j,k) - a_U[g](1,i-2,j,k));
                      duydq = c1*(a_U[g](2,i+1,j,k) - a_U[g](2,i-1,j,k))+c2*(a_U[g](2,i+2,j,k) - a_U[g](2,i-2,j,k));
                      duzdq = c1*(a_U[g](3,i+1,j,k) - a_U[g](3,i-1,j,k))+c2*(a_U[g](3,i+2,j,k) - a_U[g](3,i-2,j,k));
		   }

		   if (j == 1)
		   {   
		      duxdr =-2.25*a_U[g](1,i,j,k)+(4+fs)*a_U[g](1,i,j+1,k)-d3*a_U[g](1,i,j+2,k)+
		       3*a_U[g](1,i,j+3,k)-(1+ot)*a_U[g](1,i,j+4,k) +os*a_U[g](1,i,j+5,k);
		      duydr =-2.25*a_U[g](2,i,j,k)+(4+fs)*a_U[g](2,i,j+1,k)-d3*a_U[g](2,i,j+2,k)+
		       3*a_U[g](2,i,j+3,k)-(1+ot)*a_U[g](2,i,j+4,k) +os*a_U[g](2,i,j+5,k);
		      duzdr =-2.25*a_U[g](3,i,j,k)+(4+fs)*a_U[g](3,i,j+1,k)-d3*a_U[g](3,i,j+2,k)+
		       3*a_U[g](3,i,j+3,k)-(1+ot)*a_U[g](3,i,j+4,k) +os*a_U[g](3,i,j+5,k);
		   }
		   if (j == 2)
		   {   
		      duxdr = -os*a_U[g](1,i,j-1,k) -1.25*a_U[g](1,i,j,k)+(1+ft)*a_U[g](1,i,j+1,k)
		       - ft*a_U[g](1,i,j+2,k) + 0.5*a_U[g](1,i,j+3,k) - ot*a_U[g](1,i,j+4,k);
		      duydr = -os*a_U[g](2,i,j-1,k) -1.25*a_U[g](2,i,j,k)+(1+ft)*a_U[g](2,i,j+1,k)
		       - ft*a_U[g](2,i,j+2,k) + 0.5*a_U[g](2,i,j+3,k) - ot*a_U[g](2,i,j+4,k);
		      duzdr = -os*a_U[g](3,i,j-1,k) -1.25*a_U[g](3,i,j,k)+(1+ft)*a_U[g](3,i,j+1,k)
		       - ft*a_U[g](3,i,j+2,k) + 0.5*a_U[g](3,i,j+3,k) - ot*a_U[g](3,i,j+4,k);
		   }
		   else if (j == mEW->m_global_ny[g]-1)
		   { 
		      duxdr = os*a_U[g](1,i,j+1,k) +1.25*a_U[g](1,i,j,k)-(1+ft)*a_U[g](1,i,j-1,k) +
		        ft*a_U[g](1,i,j-2,k) - 0.5*a_U[g](1,i,j-3,k) + ot*a_U[g](1,i,j-4,k);
		      duydr = os*a_U[g](2,i,j+1,k) +1.25*a_U[g](2,i,j,k)-(1+ft)*a_U[g](2,i,j-1,k) +
		        ft*a_U[g](2,i,j-2,k) - 0.5*a_U[g](2,i,j-3,k) + ot*a_U[g](2,i,j-4,k);
		      duzdr = os*a_U[g](3,i,j+1,k) +1.25*a_U[g](3,i,j,k)-(1+ft)*a_U[g](3,i,j-1,k) +
		        ft*a_U[g](3,i,j-2,k) - 0.5*a_U[g](3,i,j-3,k) + ot*a_U[g](3,i,j-4,k);
		   }
		   else if (j == mEW->m_global_ny[g])
		   { 
		      duxdr = 2.25*a_U[g](1,i,j,k)-(4+fs)*a_U[g](1,i,j-1,k)+d3*a_U[g](1,i,j-2,k)-
		       3*a_U[g](1,i,j-3,k)+(1+ot)*a_U[g](1,i,j-4,k) -os*a_U[g](1,i,j-5,k);
		      duydr = 2.25*a_U[g](2,i,j,k)-(4+fs)*a_U[g](2,i,j-1,k)+d3*a_U[g](2,i,j-2,k)-
		       3*a_U[g](2,i,j-3,k)+(1+ot)*a_U[g](2,i,j-4,k) -os*a_U[g](2,i,j-5,k);
		      duzdr = 2.25*a_U[g](3,i,j,k)-(4+fs)*a_U[g](3,i,j-1,k)+d3*a_U[g](3,i,j-2,k)-
		       3*a_U[g](3,i,j-3,k)+(1+ot)*a_U[g](3,i,j-4,k) -os*a_U[g](3,i,j-5,k);
		   }
		   else
		   {
                      duxdr = c1*(a_U[g](1,i,j+1,k) - a_U[g](1,i,j-1,k))+c2*(a_U[g](1,i,j+2,k) - a_U[g](1,i,j-2,k));
                      duydr = c1*(a_U[g](2,i,j+1,k) - a_U[g](2,i,j-1,k))+c2*(a_U[g](2,i,j+2,k) - a_U[g](2,i,j-2,k));
                      duzdr = c1*(a_U[g](3,i,j+1,k) - a_U[g](3,i,j-1,k))+c2*(a_U[g](3,i,j+2,k) - a_U[g](3,i,j-2,k));
		   }

		   if (k == 1)
		   {
		      duxds =-2.25*a_U[g](1,i,j,k)+(4+fs)*a_U[g](1,i,j,k+1)-d3*a_U[g](1,i,j,k+2)+
		       3*a_U[g](1,i,j,k+3)-(1+ot)*a_U[g](1,i,j,k+4) +os*a_U[g](1,i,j,k+5);
		      duyds =-2.25*a_U[g](2,i,j,k)+(4+fs)*a_U[g](2,i,j,k+1)-d3*a_U[g](2,i,j,k+2)+
		       3*a_U[g](2,i,j,k+3)-(1+ot)*a_U[g](2,i,j,k+4) +os*a_U[g](2,i,j,k+5);
		      duzds =-2.25*a_U[g](3,i,j,k)+(4+fs)*a_U[g](3,i,j,k+1)-d3*a_U[g](3,i,j,k+2)+
		       3*a_U[g](3,i,j,k+3)-(1+ot)*a_U[g](3,i,j,k+4) +os*a_U[g](3,i,j,k+5);
		   }
		   if (k == 2)
		   {
		      duxds = -os*a_U[g](1,i,j,k-1) -1.25*a_U[g](1,i,j,k)+(1+ft)*a_U[g](1,i,j,k+1)
		       - ft*a_U[g](1,i,j,k+2) + 0.5*a_U[g](1,i,j,k+3) - ot*a_U[g](1,i,j,k+4);
		      duyds = -os*a_U[g](2,i,j,k-1) -1.25*a_U[g](2,i,j,k)+(1+ft)*a_U[g](2,i,j,k+1)
		       - ft*a_U[g](2,i,j,k+2) + 0.5*a_U[g](2,i,j,k+3) - ot*a_U[g](2,i,j,k+4);
		      duzds = -os*a_U[g](3,i,j,k-1) -1.25*a_U[g](3,i,j,k)+(1+ft)*a_U[g](3,i,j,k+1)
		       - ft*a_U[g](3,i,j,k+2) + 0.5*a_U[g](3,i,j,k+3) - ot*a_U[g](3,i,j,k+4);
		   }
		   else // no k=N because we are in the curvilinear grid
		   {
                      duxds = c1*(a_U[g](1,i,j,k+1) - a_U[g](1,i,j,k-1))+c2*(a_U[g](1,i,j,k+2) - a_U[g](1,i,j,k-2));
                      duyds = c1*(a_U[g](2,i,j,k+1) - a_U[g](2,i,j,k-1))+c2*(a_U[g](2,i,j,k+2) - a_U[g](2,i,j,k-2));
                      duzds = c1*(a_U[g](3,i,j,k+1) - a_U[g](3,i,j,k-1))+c2*(a_U[g](3,i,j,k+2) - a_U[g](3,i,j,k-2));
		   }
		   float_sw4 duzdy = mEW->mMetric(1,i,j,k)*duzdr+mEW->mMetric(3,i,j,k)*duzds;
		   float_sw4 duydz = mEW->mMetric(4,i,j,k)*duyds;
		   float_sw4 duxdz = mEW->mMetric(4,i,j,k)*duxds;
		   float_sw4 duzdx = mEW->mMetric(1,i,j,k)*duzdq+mEW->mMetric(2,i,j,k)*duzds;
		   float_sw4 duydx = mEW->mMetric(1,i,j,k)*duydq+mEW->mMetric(2,i,j,k)*duyds;
		   float_sw4 duxdy = mEW->mMetric(1,i,j,k)*duxdr+mEW->mMetric(3,i,j,k)*duxds;
		   float_sw4 factor = 1.0/sqrt(mEW->mJ(i,j,k));
		   a_curl[g][3*iField]   = factor*(duzdy-duydz); 
		   a_curl[g][3*iField+1] = factor*(duxdz-duzdx);
		   a_curl[g][3*iField+2] = factor*(duydx-duxdy);
		   //                   iField++; 
                }
      }
   }
}

//-----------------------------------------------------------------------
void Image::computeImageMagdt( vector<Sarray> &a_Up, vector<Sarray> &a_Um,
			       float_sw4 dt )
{
   // dt is distance between Up and Um.
   ASSERT(m_isDefinedMPIWriters);
// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
         gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      float_sw4 factor = 1.0/dt;
      for (int g = gmin; g <= gmax; g++)
      {
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		  float_sw4 vmag = factor*sqrt( 
                     (a_Up[g](1,i,j,k) - a_Um[g](1,i,j,k))*(a_Up[g](1,i,j,k) - a_Um[g](1,i,j,k)) +
		     (a_Up[g](2,i,j,k) - a_Um[g](2,i,j,k))*(a_Up[g](2,i,j,k) - a_Um[g](2,i,j,k)) +
		     (a_Up[g](3,i,j,k) - a_Um[g](3,i,j,k))*(a_Up[g](3,i,j,k) - a_Um[g](3,i,j,k)) );
		  if (m_double)
		     m_doubleField[g][iField] = (double) vmag;
		  else
		     m_floatField[g][iField] = (float) vmag;
		  //		  iField++; 
	       }
      }
   }
} 

//-----------------------------------------------------------------------
void Image::computeImageMag( vector<Sarray> &a_U )
{
   // dt is distance between Up and Um.
   ASSERT(m_isDefinedMPIWriters);
// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
         gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      for (int g = gmin; g <= gmax; g++)
      {
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		  float_sw4 mag = sqrt(  a_U[g](1,i,j,k)*a_U[g](1,i,j,k) +
			       a_U[g](2,i,j,k)*a_U[g](2,i,j,k) +
			       a_U[g](3,i,j,k)*a_U[g](3,i,j,k) );
		  if (m_double)
		     m_doubleField[g][iField] = (double) mag;
		  else
		     m_floatField[g][iField] = (float) mag;
		  //		  iField++; 
	       }
      }
   }
} 

//-----------------------------------------------------------------------
void Image::computeImageHmagdt(vector<Sarray> &a_Up, vector<Sarray> &a_Um,
			       float_sw4 dt )
{
   // dt is distance between Up and Um.
   ASSERT(m_isDefinedMPIWriters);
// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
         gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      float_sw4 factor = 1.0/dt;
      for (int g = gmin; g <= gmax; g++)
      {
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		  float_sw4 vmag = factor*sqrt( 
                     (a_Up[g](1,i,j,k) - a_Um[g](1,i,j,k))*(a_Up[g](1,i,j,k) - a_Um[g](1,i,j,k)) +
		     (a_Up[g](2,i,j,k) - a_Um[g](2,i,j,k))*(a_Up[g](2,i,j,k) - a_Um[g](2,i,j,k))); 
		  if (m_double)
		     m_doubleField[g][iField] = (double) vmag;
		  else
		     m_floatField[g][iField] = (float) vmag;
		  //		  iField++; 
	       }
      }
   }
} 

//-----------------------------------------------------------------------
void Image::computeImageHmag( vector<Sarray> &a_U )
{
   // dt is distance between Up and Um.
   ASSERT(m_isDefinedMPIWriters);
// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
         gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      for (int g = gmin; g <= gmax; g++)
      {
	 //         size_t iField = 0;
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		  float_sw4 mag = sqrt( a_U[g](1,i,j,k)*a_U[g](1,i,j,k)+a_U[g](2,i,j,k)*a_U[g](2,i,j,k));
		  if (m_double)
		     m_doubleField[g][iField] = (double) mag;
		  else
		     m_floatField[g][iField] = (float) mag;
		  //		  iField++; 
	       }
      }
   }
} 

//-----------------------------------------------------------------------
void Image::compute_image_gradp( vector<Sarray>& a_gLambda, vector<Sarray>& a_Mu,
				 vector<Sarray>& a_Lambda, vector<Sarray>& a_Rho )
{
   ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
         gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      for (int g = gmin; g <= gmax; g++)
      {
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		  float_sw4 gradp = a_gLambda[g](i,j,k)*2*sqrt( (2*a_Mu[g](i,j,k)+a_Lambda[g](i,j,k))*a_Rho[g](i,j,k) );
		  if (m_double)
		     m_doubleField[g][iField] = gradp;
		  else
		     m_floatField[g][iField] = (float) gradp;
		  //		  iField++; 
	       }
      }
   }
}

//-----------------------------------------------------------------------
void Image::compute_image_grads( vector<Sarray>& a_gMu, vector<Sarray>& a_gLambda, 
				 vector<Sarray>& a_Mu, vector<Sarray>& a_Rho )
{
   ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
         gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      for (int g = gmin; g <= gmax; g++)
      {
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		  float_sw4 grads = (2*a_gMu[g](i,j,k)-4*a_gLambda[g](i,j,k))*2*sqrt( a_Mu[g](i,j,k)*a_Rho[g](i,j,k));
		  if (m_double)
		     m_doubleField[g][iField] = (double) grads;
		  else
		     m_floatField[g][iField] = (float) grads;
		  //		  iField++; 
	       }
      }
   }
}

//-----------------------------------------------------------------------
void Image::update_image( int a_cycle, float_sw4 a_time, float_sw4 a_dt,
			  vector<Sarray>& a_Up,  vector<Sarray>& a_U, vector<Sarray>& a_Um,
			  vector<Sarray>& a_Rho, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
			  vector<Sarray>& a_gRho, vector<Sarray>& a_gMu, vector<Sarray>& a_gLambda,
                          vector<Source*>& a_sources, int a_dminus )
{
   if( mMode == HMAXDUDT)
   {
      if( a_dminus )
	 update_maxes_hVelMax( a_Up, a_U, a_dt );
      else
	 update_maxes_hVelMax( a_Up, a_Um, 2*a_dt );
   }
   if( mMode == HMAX )
      update_maxes_hMax( a_Up );

   if( mMode == VMAXDUDT)
   {
      if( a_dminus )
	 update_maxes_vVelMax( a_Up, a_U, a_dt );
      else
	 update_maxes_vVelMax( a_Up, a_Um, 2*a_dt );
   }
   if( mMode == VMAX )
      update_maxes_vMax( a_Up );

      // Center time derivatives around t-dt, i.e., (up-um)/(2*dt), except when dminus
      // is set. Use (up-u)/dt assumed centered at t, when dminus is true.
   int td = 0;
   if( !a_dminus )
      td = is_time_derivative();

   if( timeToWrite( a_time-td*a_dt , a_cycle-td, a_dt ) )
      output_image( a_cycle, a_time, a_dt, a_Up, a_U, a_Um, a_Rho, a_Mu, a_Lambda,
		    a_gRho, a_gMu, a_gLambda, a_sources, a_dminus );
}

//-----------------------------------------------------------------------
void Image::output_image( int a_cycle, float_sw4 a_time, float_sw4 a_dt,
			  vector<Sarray>& a_Up,  vector<Sarray>& a_U, vector<Sarray>& a_Um,
			  vector<Sarray>& a_Rho, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
			  vector<Sarray>& a_gRho, vector<Sarray>& a_gMu, vector<Sarray>& a_gLambda,
                          vector<Source*>& a_sources, int a_dminus )
{


   int td = 0;
   if( !a_dminus )
      td = is_time_derivative();

   if(mMode == UX ) 
      computeImageQuantity(a_Up, 1);
   else if(mMode == UY )
      computeImageQuantity(a_Up, 2);
   else if(mMode == UZ )
      computeImageQuantity(a_Up, 3);
   else if(mMode == RHO )
      computeImageQuantity(a_Rho, 1);
   else if(mMode == MU )
      computeImageQuantity(a_Mu, 1);
   else if(mMode == LAMBDA )
      computeImageQuantity(a_Lambda, 1);
   else if(mMode == P )
      computeImagePvel(a_Mu, a_Lambda, a_Rho);
   else if(mMode == S )
      computeImageSvel(a_Mu, a_Rho);
   else if(mMode == DIV || mMode == DIVDT 
	   || mMode == CURLMAG || mMode == CURLMAGDT )
      computeImageDivCurl( a_Up, a_U, a_Um, a_dt, a_dminus );
   else if(mMode == LAT || mMode == LON )
      computeImageLatLon( mEW->mX, mEW->mY, mEW->mZ );
   else if(mMode == TOPO )
   {
      if (mEW->topographyExists())
	 copy2DArrayToImage(mEW->mTopo); // save the raw topography; the smoothed is saved by the mode=grid with z=0
   }
   if( mMode == UXEXACT || mMode == UYEXACT || mMode == UZEXACT ||
       mMode == UXERR || mMode == UYERR || mMode == UZERR )
   {
      vector<Sarray> Uex(mEW->mNumberOfGrids);
      vector<Sarray*> alpha; //dummy, the array is not used in routine exactSol.
      for( int g=0 ; g < mEW->mNumberOfGrids ; g++ )
	 Uex[g].define(3,mEW->m_iStart[g],mEW->m_iEnd[g],mEW->m_jStart[g],mEW->m_jEnd[g],mEW->m_kStart[g],mEW->m_kEnd[g]);
      mEW->exactSol( a_time, Uex, alpha, a_sources );
      if( mMode == UXERR )
	 computeImageQuantityDiff(a_Up,Uex,1);
      else if( mMode == UYERR )
	 computeImageQuantityDiff(a_Up,Uex,2);
      else if( mMode == UZERR ) 
	 computeImageQuantityDiff(a_Up,Uex,3);
      else if( mMode == UXEXACT )
	 computeImageQuantity(Uex,1);
      else if( mMode == UYEXACT )
	 computeImageQuantity(Uex,2);
      else if( mMode == UZEXACT )
	 computeImageQuantity(Uex,3);
      Uex.clear();
   }
   else if( mMode == GRIDX || mMode == GRIDY || mMode == GRIDZ )
      computeImageGrid(mEW->mX, mEW->mY, mEW->mZ );
   else if(mMode == MAGDUDT )
   {
      if( a_dminus )
	 computeImageMagdt( a_Up, a_U, a_dt );
      else
	 computeImageMagdt( a_Up, a_Um, 2*a_dt );
   }
   else if(mMode == HMAGDUDT )
   {
      if( a_dminus )
	 computeImageHmagdt( a_Up, a_U, a_dt );
      else
	 computeImageHmagdt( a_Up, a_Um, 2*a_dt );
   }
   else if(mMode == MAG )
      computeImageMag( a_Up );
   else if(mMode == HMAG )
      computeImageHmag( a_Up );
   else if(mMode == GRADRHO )
      computeImageQuantity( a_gRho, 1 );
   else if(mMode == GRADMU )
      computeImageQuantity( a_gMu, 1 );
   else if(mMode == GRADLAMBDA )
      computeImageQuantity( a_gLambda, 1 );
   else if(mMode == GRADP )
      compute_image_gradp( a_gLambda, a_Mu, a_Lambda, a_Rho );
   else if(mMode == GRADS )
      compute_image_grads( a_gMu, a_gLambda, a_Mu, a_Rho );
   else if (mMode != HMAXDUDT || mMode != VMAXDUDT
	    || mMode != HMAX   || mMode != VMAX )
   {
      if (mEW->proc_zero())
      {
	    //	  printf("Can only write ux, uy, uz, mu, rho, lambda, uxerr, uyerr, uzerr- remove once completely implemented\n");
	 printf("Can only write ux, uy, uz, mu, rho, lambda: - remove once completely implemented\n");
	 printf("I can not print data of type %i\n", mMode );
      }
      MPI_Abort(MPI_COMM_WORLD,1);
   }

// write the image plane on file    
   double t[2];
   t[0] = t[1] = MPI_Wtime();
   string path = mEW->getPath();
   writeImagePlane_2( a_cycle-td, path, a_time-td*a_dt );
   t[1] = MPI_Wtime();
      
// output timing info?
   bool iotiming=false;
   if( iotiming )
   {
      t[0] = t[1]-t[0];
      double tmp[2];
      tmp[0] = t[0];
      tmp[1] = t[1];
      MPI_Reduce( tmp, t, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
      if( mEW->proc_zero() )
      {
	 cout << "Maximum write time:";
	 cout << " (using Bjorn's I/O library) " << t[1] << " seconds. (<=" << m_pio[0]->n_writers() << " procs writing)";
	 cout << endl;
      } // end if proc_zero
   } // end if iotiming      
}


//-----------------------------------------------------------------------
void Image::computeImageQuantityDiff( vector<Sarray>& a_U, vector<Sarray>& a_Uex,
				      int comp )
{
   ASSERT(m_isDefinedMPIWriters);

// plane_in_proc returns true for z=const planes, because all processors have a part in these planes
   bool iwrite   = plane_in_proc(m_gridPtIndex[0]);
   if (iwrite)
   {
      int gmin, gmax;
      if( mLocationType == Image::Z )
	 gmin = gmax = m_gridPtIndex[1];
      else
      {
         gmin = 0;
	 gmax = mEW->mNumberOfGrids-1;
      }
      for (int g = gmin; g <= gmax; g++)
      {
	 //         size_t iField = 0;
	 size_t ni = (mWindow[g][1]-mWindow[g][0]+1);
	 size_t nij= (mWindow[g][1]-mWindow[g][0]+1)*(mWindow[g][3]-mWindow[g][2]+1);
#pragma omp parallel for
         for (int k = mWindow[g][4]; k <= mWindow[g][5]; k++)
	    for (int j = mWindow[g][2]; j <= mWindow[g][3]; j++)
	       for (int i = mWindow[g][0]; i <= mWindow[g][1]; i++)
	       {
		  size_t iField = (i-mWindow[g][0])+ni*(j-mWindow[g][2])+nij*(k-mWindow[g][4]);
		  if (m_double)
		     m_doubleField[g][iField] = (double) a_U[g](comp,i,j,k)-a_Uex[g](comp,i,j,k);
		  else
		     m_floatField[g][iField] = (float) a_U[g](comp,i,j,k)-a_Uex[g](comp,i,j,k);
		  //		  iField++; 
	       }
      }
   }
}

//-----------------------------------------------------------------------
bool Image::needs_mgrad() const
{
   return mMode==GRADRHO || mMode==GRADMU || mMode==GRADLAMBDA || 
      mMode==GRADP || mMode==GRADS;
}
