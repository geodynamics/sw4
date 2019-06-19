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
#ifndef EW_IMAGE_H
#define EW_IMAGE_H

#include <string>
#include <sstream>
#include <vector>

#include "boundaryConditionTypes.h"
#include "Sarray.h"
#include "Parallel_IO.h"

class EW;

class Image
{
public:

   enum ImageMode { NONE, UX, UY, UZ, RHO, LAMBDA, MU, P, S, UXEXACT, UYEXACT, UZEXACT, DIV, CURLMAG,
		    DIVDT, CURLMAGDT, LAT, LON, TOPO, GRIDX, GRIDY, GRIDZ, UXERR, UYERR, UZERR, 
		    MAGDUDT, HMAGDUDT, HMAXDUDT, VMAXDUDT, MAG, HMAG, HMAX, VMAX,
		    GRADRHO, GRADMU, GRADLAMBDA, GRADP, GRADS, QP, QS }; 

static int MODES;
static Image* nil;

enum ImageOrientation {UNDEFINED, X, Y, Z};

Image(EW * a_ew,
      float_sw4 time, 
      float_sw4 timeInterval, 
      int cycle, 
      int cycleInterval,
      const std::string& filePrefix, 
      ImageMode mode,
      ImageOrientation locationType, 
      float_sw4 locationValue,
      bool doubleMode, bool userCreated=true );

static  void setSteps(int a_steps);

// static void setTiming(float startTime, float dt);
// static void setGridAttributes(std::vector<double> a_gridSize      ,
// 			      std::vector<double> a_zmin          ,
// 			      const int           a_n_ghost_points,
// 			      const int           a_padding       );

void set_double( bool val=true ); 
  /* Here, we compute the index --in the local grid-- of the coordinate value at which we have to plot. 
     For x, y, the local grid is the same as the global grid, but for z, k resets at each refinement boundary. */
void computeGridPtIndex();
void allocatePlane();

void computeImageQuantity(std::vector<Sarray> &a_mu, int a_nComp);
void computeImagePvel(std::vector<Sarray> &mu, std::vector<Sarray> &lambda,
		      std::vector<Sarray> &rho );
void computeImageSvel(std::vector<Sarray> &mu, std::vector<Sarray> &rho );
   void computeImageGrid( std::vector<Sarray> &a_X, std::vector<Sarray> &a_Y, std::vector<Sarray> &a_Z );
   void computeImageLatLon( std::vector<Sarray> &a_X, std::vector<Sarray> &a_Y, std::vector<Sarray> &a_Z );
void computeImageDivCurl( std::vector<Sarray> &a_Up, std::vector<Sarray>& a_U,
			  std::vector<Sarray> &a_Um, float_sw4 dt, int dminus );
void computeImageMagdt( std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, float_sw4 dt );
void computeImageMag( std::vector<Sarray> &a_U );
void computeImageHmagdt( std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, float_sw4 dt );
void computeImageHmag( std::vector<Sarray> &a_U );
void compute_image_gradp( vector<Sarray>& a_gLambda, vector<Sarray>& a_Mu,
			  vector<Sarray>& a_Lambda, vector<Sarray>& a_Rho );
void compute_image_grads( vector<Sarray>& a_gMu, vector<Sarray>& a_gLambda, 
			  vector<Sarray>& a_Mu, vector<Sarray>& a_Rho );

void computeImageQuantityDiff( vector<Sarray>& a_U, vector<Sarray>& a_Uex,
			       int comp );

void output_image( int a_cycle, float_sw4 a_time, float_sw4 a_dt,
		   vector<Sarray>& a_Up,  vector<Sarray>& a_U, vector<Sarray>& a_Um,
		   vector<Sarray>& a_Rho, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		   vector<Sarray>& a_gRho, vector<Sarray>& a_gMu, vector<Sarray>& a_gLambda,
		   vector<Source*>& a_sources, int a_dminus );

void update_image( int a_cycle, float_sw4 a_time, float_sw4 a_dt,
		   vector<Sarray>& a_Up,  vector<Sarray>& a_U, vector<Sarray>& a_Um,
		   vector<Sarray>& a_Rho, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		   vector<Sarray>& a_gRho, vector<Sarray>& a_gMu, vector<Sarray>& a_gLambda,
		   vector<Source*>& a_sources, int a_dminus );

//void computeImageError(std::vector<Sarray> &a_mu, int a_nComp);

void copy2DArrayToImage(Sarray &twoDimensionalArray);

bool is_time_derivative() const;

void associate_gridfiles( vector<Image*>& imgs );
   
void writeImagePlane_2(int cycle, std::string &a_path, float_sw4 time );
void add_grid_filenames_to_file( const char* fname );
void add_grid_to_file( const char* fname, bool iwrite, size_t offset );

// several curvilinear grids (MR)
void add_grids_to_file( const char* fname, bool iwrite, size_t offset );

bool plane_in_proc(int a_gridIndexCoarsest);
void initializeIO();

void computeNearestGridPoint(std::vector<int> & a_arrayIndex, float_sw4 x) const;
  
void update_maxes_vVelMax( std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, float_sw4 dt );
void update_maxes_hVelMax( std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, float_sw4 dt );
void update_maxes_vMax( std::vector<Sarray> &a_U );
void update_maxes_hMax( std::vector<Sarray> &a_U );

const std::string fieldSuffix(ImageMode mode) const;

bool timeToWrite(float_sw4 time, int cycle, float_sw4 dt );
bool timeToWrite( int cycle );
void compute_file_suffix( std::stringstream & fileSuffix, int cycle );

ImageOrientation getOrientation() const {return mLocationType;};

ImageMode mMode;
std::string mFilePrefix;
void initializeTime(double t=0.0);
bool needs_mgrad() const;

protected:

void define_pio();  
   //bool proc_write();

float_sw4 mTime;
bool m_time_done;
float_sw4 mTimeInterval;
float_sw4 mNextTime;
int mWritingCycle;
int mCycleInterval;
//std::string mFileName;
  
std::vector<std::string> mMode2Suffix;
std::vector<std::string> mOrientationString;
  
static int mPreceedZeros; // number of digits for unique time step in file names

   //bool m_gridPtValueInitialized;
  
private:

Image(); // make it impossible to call default constructor
Image(const Image &im); // hide copy constructor 

void computeDivergence( std::vector<Sarray> &a_U, std::vector<float_sw4*>& a_div );
void computeCurl( std::vector<Sarray> &a_U, std::vector<float_sw4*>& a_curl );

   //bool mWriting;
   //bool mReadyToWrite;

ImageOrientation mLocationType;
float_sw4 mCoordValue;
std::vector<int> m_gridPtIndex;

MPI_Comm m_mpiComm_writers;
bool m_isDefinedMPIWriters;

   //int mNonempty;

int mGridinfo;   // -1 = undefined, 0=Cartesian patches, 1=Curvilinear grid appended, 2=grid file names appended.
bool mStoreGrid; // true=append curvilinear grid to image, false=append grid file names to image.
Image* m_gridimage; // Curvilinear grid z-coordinate 
bool m_user_created; // true --> This image was created from the input file

   //int m_rankWriter;
bool m_isDefined;
bool m_double;
EW* mEW;
Parallel_IO** m_pio;

// moved to class EW
//int m_pfs, m_nwriters;

std::vector<int*> mWindow; // start + end indices for (i,j,k) for each grid level

// pointers to in-plane data (only one will get allocated by Image::allocatePlane()
std::vector<double*> m_doubleField;
std::vector<float*> m_floatField;

};

#endif
