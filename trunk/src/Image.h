//-*-c++-*-
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

// WPP modes:
// enum ImageMode { NONE, UX, UY, UZ, RHO, LAMBDA, MU, P, S,
// 		 DIV, CURL, VELDIV, VELCURL, LAT, LON,
// 		 HVELMAX, VVELMAX, TOPO, GRID, UXERR, UYERR, UZERR, 
// 		 FX, FY, FZ, VELMAG, QS, QP, HVEL }; // NONE + 28 modes = 29 are currently defined

enum ImageMode { NONE, UX, UY, UZ, RHO, LAMBDA, MU, P, S, UXEXACT, UYEXACT, UZEXACT, DIV, CURLMAG,
                 DIVDT, CURLMAGDT, LAT, LON, TOPO, GRIDX, GRIDY, GRIDZ, UXERR, UYERR, UZERR, 
                 MAGDUDT, HMAGDUDT, HMAXDUDT, VMAXDUDT, MAG, HMAG, HMAX, VMAX }; 

static int MODES;
static Image* nil;

enum ImageOrientation {UNDEFINED, X, Y, Z};

Image(EW * a_ew,
      double time, 
      double timeInterval, 
      int cycle, 
      int cycleInterval,
      const std::string& filePrefix, 
      ImageMode mode,
      ImageOrientation locationType, 
      double locationValue,
      bool doubleMode );

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
void computeImageGrid( Sarray &a_X, Sarray &a_Y, Sarray &a_Z );
void computeImageLatLon( Sarray &a_X, Sarray &a_Y, Sarray &a_Z );
void computeImageDivCurl( std::vector<Sarray> &a_Up, std::vector<Sarray>& a_U,
			  std::vector<Sarray> &a_Um, double dt, int dminus );
void computeImageMagdt( std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, double dt );
void computeImageMag( std::vector<Sarray> &a_U );
void computeImageHmagdt( std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, double dt );
void computeImageHmag( std::vector<Sarray> &a_U );
//void computeImageError(std::vector<Sarray> &a_mu, int a_nComp);

void copy2DArrayToImage(Sarray &twoDimensionalArray);

bool is_time_derivative() const;

void associate_gridfiles( vector<Image*> imgs );
   
void writeImagePlane_2(int cycle, std::string &a_path, double time );
void add_grid_filenames_to_file( const char* fname );
void add_grid_to_file( const char* fname, bool iwrite, size_t offset );

bool plane_in_proc(int a_gridIndexCoarsest);
void initializeIO();

void computeNearestGridPoint(std::vector<int> & a_arrayIndex, double x) const;
  
void update_maxes_vVelMax( std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, double dt );
void update_maxes_hVelMax( std::vector<Sarray> &a_Up, std::vector<Sarray> &a_Um, double dt );
void update_maxes_vMax( std::vector<Sarray> &a_U );
void update_maxes_hMax( std::vector<Sarray> &a_U );

const std::string fieldSuffix(ImageMode mode) const;

bool timeToWrite(double time, int cycle, double dt );
void compute_file_suffix( std::stringstream & fileSuffix, int cycle );

ImageOrientation getOrientation() const {return mLocationType;};

ImageMode mMode;
std::string mFilePrefix;
void initializeTime();

protected:
void define_pio();  

bool proc_write();

double mTime;
bool m_time_done;
double mTimeInterval;
double mNextTime;
int mWritingCycle;
int mCycleInterval;
std::string mFileName;
  
std::vector<std::string> mMode2Suffix;
std::vector<std::string> mOrientationString;
  
static int mPreceedZeros; // number of digits for unique time step in file names

bool m_gridPtValueInitialized;
  
private:
Image(); // make it impossible to call default constructor
Image(const Image &im); // hide copy constructor 

void computeDivergence( std::vector<Sarray> &a_U, std::vector<double*>& a_div );
void computeCurl( std::vector<Sarray> &a_U, std::vector<double*>& a_curl );


bool mWriting;
bool mReadyToWrite;

ImageOrientation mLocationType;
double mCoordValue;
std::vector<int> m_gridPtIndex;

MPI_Comm m_mpiComm_writers;
bool m_isDefinedMPIWriters;

int mNonempty;

int mGridinfo;   // -1 = undefined, 0=Cartesian patches, 1=Curvilinear grid appended, 2=grid file names appended.
bool mStoreGrid; // true=append curvilinear grid to image, false=append grid file names to image.
Image* m_gridimage1; // Curvilinear grid coordinate 1
Image* m_gridimage2;  // Curvilinear grid coordinate 2

int m_rankWriter;
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
