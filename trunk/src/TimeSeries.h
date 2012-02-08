// -*-c++-*-
#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <vector>

class EW;

using namespace std;

class TimeSeries{

public:

// support for the time derivtives are not yet implemented (mainly header info is missing)
enum receiverMode{Solution, Div, Curl, Strains /*, SolutionTderiv, DivTderiv, CurlTderiv, StrainsTderiv */ };

TimeSeries( EW* a_ew, std::string name, receiverMode mode, bool sacFormat, bool usgsFormat, 
	    double x, double y, double z, bool topoDepth, int writeEvery );
~TimeSeries();

void allocateRecordingArrays( int numberOfTimeSteps, double startTime, double timeStep );
  
void recordData(vector<double> & u);

void writeFile( );

double **getRecordingArray(){ return mRecordedSol; }

bool myPoint(){ return m_myPoint; }

// for simplicity, make the grid point location public
int m_i0;
int m_j0;
int m_k0;
int m_grid0;

private:   
TimeSeries();
void write_usgs_format( string a_fileName);

receiverMode m_mode;
int m_nComp;

bool m_myPoint; // set to true if this processor writes to the arrays

std::string m_fileName, m_filePrefix;

double mX, mY, mZ, mGPX, mGPY, mGPZ; // original and actual location
bool m_zRelativeToTopography; // location is given relative to topography

int mWriteEvery;
bool m_usgsFormat, m_sacFormat;


// start time and time step 
double m_t0, m_dt;

// size of recording arrays
int mAllocatedSize;

// index of last recorded element
int mLastTimeStep;

// recording arrays
double** mRecordedSol;

// ignore this station if it is above the topography or outside the computational domain
bool mIgnore;

// what are all these numbers used for? SAC header?
// bool m_xycomponent, m_velocities;
// double m_calpha, m_salpha, m_thxnrm, m_thynrm, m_lat, m_lon, m_dthi;
// double m_dmx, m_dmy, m_dmz, m_d0x, m_d0y, m_d0z;
// double m_dmxy, m_dmxz, m_dmyz, m_d0xy, m_d0xz, m_d0yz;
};


#endif
