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

// enum ImageMode { NONE, UX, UY, UZ, RHO, LAMBDA, MU, P, S,
// 		 DIV, CURL, VELDIV, VELCURL, LAT, LON,
// 		 HVELMAX, VVELMAX, TOPO, GRID, UXERR, UYERR, UZERR, 
// 		 FX, FY, FZ, VELMAG, QS, QP, HVEL }; // NONE + 28 modes = 29 are currently defined

// start with the most basic modes
enum ImageMode { NONE, UX, UY, UZ, RHO, LAMBDA, MU}; // NONE + 6 modes = 7 are currently defined
static int MODES;
   

enum ImageOrientation {UNDEFINED, X, Y, Z};

Image(EW * a_ew,
      float time, 
      float timeInterval, 
      int cycle, 
      int cycleInterval,
      const std::string& filePrefix, 
      int sample,
      ImageMode mode,
      ImageOrientation locationType, 
      float locationValue,
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
//void computeImageError(std::vector<Sarray> &a_mu, int a_nComp);
void computeImageP(std::vector<Sarray> &mu, std::vector<Sarray> &lambda, std::vector<Sarray> &rho);
void computeImageS(std::vector<Sarray> &mu, std::vector<Sarray> &lambda, std::vector<Sarray> &rho);
void copy2DArrayToImage(Sarray &twoDimensionalArray);
void evaluateGridImage(Sarray &a_X, Sarray &a_Y, Sarray &a_Z, int a_component);
void evaluateLatLonImage(Sarray &a_X, Sarray &a_Y, Sarray &a_Z, int a_component);

void writeImagePlane_2(int cycle, std::string &a_path);

bool plane_in_proc(int a_gridIndexCoarsest);

void initializeIO();

void computeNearestGridPoint(std::vector<int> & a_arrayIndex, double x) const;
  
// we need to add vector<Sarray> arguments to these functions
// void update_maxes_hVelMax();
// void update_maxes_vVelMax();

// void computeDivergence();
// void computeCurlMagnitude();
// void computeVelocityMagnitude();
// void computeVelCurlMagnitude();
// void computeVelDivergence();
// void computeHorizontalVelocityMagnitude();

//double  computeImageErrorDebug(int whichCase);

const std::string fieldSuffix(ImageMode mode) const;

bool timeToWrite(double time, int cycle, double dt );

void compute_file_suffix( std::stringstream & fileSuffix );
ImageOrientation getOrientation(){return mLocationType;};

ImageMode mMode;
std::string mFilePrefix;

protected:
void define_pio();  
void initializeTime();
bool proc_write();

double mTime;
bool m_time_done;
double mTimeInterval;
double mNextTime;
int mWritingCycle;
int mCycleInterval;
int mImageSamplingFactor;
std::string mFileName;
  
std::vector<std::string> mMode2Suffix;
std::vector<std::string> mOrientationString;
  
static int mPreceedZeros; // number of digits for unique time step in file names

bool m_gridPtValueInitialized;
int mCycle;
  
private:
Image(); // make it impossible to call default constructor
Image(const Image &im); // hide copy constructor 

bool mWriting;
bool mReadyToWrite;

ImageOrientation mLocationType;
float mCoordValue;
std::vector<int> m_gridPtIndex;

MPI_Comm m_mpiComm_writers;
bool m_isDefinedMPIWriters;

int mNonempty;

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
