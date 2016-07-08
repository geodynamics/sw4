#ifndef SW4_CHECKPOINT_H
#define SW4_CHECKPOINT_H

#include <string>
#include <sstream>
#include <vector>

#include "sw4.h"
//#include "boundaryConditionTypes.h"
#include "Sarray.h"
#include "Parallel_IO.h"

class EW;

class CheckPoint
{
public:
   static CheckPoint* nil; // nil pointer
   CheckPoint( EW * a_ew,
	       float_sw4 time, 
	       float_sw4 timeInterval, 
	       int cycle, 
	       int cycleInterval,
	       string fname,
	       size_t bufsize=10000000 );
   CheckPoint( EW * a_ew, string fname, size_t bufsize=10000000 );
   ~CheckPoint();
   void write_checkpoint( float_sw4 a_time, int a_cycle, std::vector<Sarray>& a_Um,
			  std::vector<Sarray>& a_U );
   void read_checkpoint( float_sw4& a_time, int& a_cycle, std::vector<Sarray>& a_Um,
			 std::vector<Sarray>& a_U );
   void setup_sizes();
   bool timeToWrite( float_sw4 time, int cycle, float_sw4 dt );

protected:
   void define_pio( );
   void setSteps( int a_steps );
   void compute_file_suffix( int cycle, std::stringstream & fileSuffix );

   std::string mFilePrefix;
   std::string mRestartFile;
   float_sw4 mTime;
   float_sw4 mTimeInterval;
   bool m_time_done;
   float_sw4 mNextTime;
   float_sw4 mStartTime;

   int mWritingCycle;
   int mCycleInterval;
   bool m_winallocated;
   size_t m_bufsize;

private:
   CheckPoint(); // make it impossible to call default constructor
   CheckPoint(const CheckPoint &cp ); // hide copy constructor 
   int mPreceedZeros; // number of digits for unique time step in file names
   bool m_double;
   EW* mEW;
   Parallel_IO** m_parallel_io;
   std::vector<int*> mWindow; // Local in processor start + end indices for (i,j,k) for each grid level
   std::vector<int*> mGlobalDims; // Global start + end indices for (i,j,k) for each grid level
   std::vector<bool> m_ihavearray;
};

#endif
