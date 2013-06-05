#ifndef SW4_IMAGE3D_H
#define SW4_IMAGE3D_H

#include <string>
#include <sstream>
#include <vector>

#include "boundaryConditionTypes.h"
#include "Sarray.h"
#include "Parallel_IO.h"

class EW;

class Image3D
{
public:
  enum Image3DMode { NONE, UX, UY, UZ, RHO, LAMBDA, MU, P, S, GRADRHO, GRADMU, GRADLAMBDA, GRADP, GRADS, QP, QS };
// 15 modes are currently defined

   static Image3D* nil;

   Image3D( EW * a_ew,
	    double time, 
	    double timeInterval, 
	    int cycle, 
	    int cycleInterval,
	    double tstart,
	    const std::string& filePrefix, 
	    Image3DMode mode,
	    bool doubleMode );
   ~Image3D();

   void setup_images( );

   static void setSteps(int a_steps);

   //   void set_double( bool val=true ); 

   void update_image( int a_cycle, double a_time, double a_dt, std::vector<Sarray>& a_U, 
		      std::vector<Sarray>& a_Rho, std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda,
		      std::vector<Sarray>& a_gRho, std::vector<Sarray>& a_gMu, std::vector<Sarray>& a_gLambda,
		      std::vector<Sarray>& a_Qp, std::vector<Sarray>& a_Qs,
		      std::string a_path, Sarray& a_Z );

   void force_write_image( double a_time, int a_cycle, vector<Sarray>& a_U, 
			   vector<Sarray>& a_Rho, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
			   vector<Sarray>& a_gRho, vector<Sarray>& a_gMu, vector<Sarray>& a_gLambda,
			   vector<Sarray>& a_Qp, vector<Sarray>& a_Qs,
			   std::string a_path, Sarray& a_Z );

   //   void set_start_time(double tStart);

protected:
   bool timeToWrite( double time, int cycle, double dt );

   void compute_image( std::vector<Sarray>& a_U, std::vector<Sarray>& a_Rho,
		       std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda,
		       std::vector<Sarray>& a_gRho, std::vector<Sarray>& a_gMu,
		       std::vector<Sarray>& a_gLambda,
     		       std::vector<Sarray>& a_Qp, std::vector<Sarray>& a_Qs );

   void write_image( int cycle, std::string &path, double t, Sarray& a_Z );

   void define_pio( );

   void compute_file_suffix( int cycle, std::stringstream & fileSuffix );

   Image3DMode mMode;
   std::string mFilePrefix;
   double mTime;
   bool m_time_done;
   double mTimeInterval;
   double mNextTime;
   double mStartTime;

   int mWritingCycle;
   int mCycleInterval;
   int mImageSamplingFactor;

   std::string mFileName;
   std::string m_modestring;

   bool m_isDefinedMPIWriters;
   bool m_winallocated;
   bool m_memallocated;

   static int mPreceedZeros; // number of digits for unique time step in file names
  
private:
   Image3D(); // make it impossible to call default constructor
   Image3D(const Image &im); // hide copy constructor 

   bool m_double;
   EW* mEW;
   Parallel_IO ** m_parallel_io;

   std::vector<int*> mWindow; // Local in processor start + end indices for (i,j,k) for each grid level
   std::vector<int*> mGlobalDims; // Global start + end indices for (i,j,k) for each grid level

   std::vector<double*> m_doubleField;
   std::vector<float*> m_floatField;
   std::vector<int> m_extraz;
   std::vector<bool> m_ihavearray;
};

#endif
