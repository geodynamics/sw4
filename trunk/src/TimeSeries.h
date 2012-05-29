// -*-c++-*-
#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <vector>
#include <string>

class EW;
class Sarray;

using namespace std;

class TimeSeries{

public:

// support for derived quantities of the time derivative are not yet implemented
enum receiverMode{Displacement, Div, Curl, Strains, Velocity /*, DivVelo, CurlVelo, StrainsVelo */ };

TimeSeries( EW* a_ew, std::string name, receiverMode mode, bool sacFormat, bool usgsFormat, 
	    double x, double y, double z, bool topoDepth, int writeEvery, bool xyzcomponent=true );
~TimeSeries();

void allocateRecordingArrays( int numberOfTimeSteps, double startTime, double timeStep );
  
void recordData(vector<double> & u);

void writeFile( );

void readFile( );

double **getRecordingArray(){ return mRecordedSol; }

int getNsteps() const {return mLastTimeStep+1;}

bool myPoint(){ return m_myPoint; }

receiverMode getMode(){ return m_mode; }

double getX() const {return mX;}
double getY() const {return mY;}
double getZ() const {return mZ;}

double arrival_time( double lod );

TimeSeries* copy( EW* a_ew, string filename, bool addname=false );

double misfit( TimeSeries& observed, TimeSeries* diff );

void interpolate( TimeSeries& intpfrom );

void use_as_forcing( int n, std::vector<Sarray>& f, std::vector<double> & h, double dt );

double product( TimeSeries& ts );
double product_wgh( TimeSeries& ts );

void set_station_utc( int utc[7] );
void offset_ref_utc( int utc[7] );
double utc_distance( int utc1[7], int utc2[7] );
void dayinc( int date[7] );
int lastofmonth( int year, int month );

// for simplicity, make the grid point location public
int m_i0;
int m_j0;
int m_k0;
int m_grid0;

private:   
TimeSeries();
void write_usgs_format( string a_fileName);
void write_sac_format( int npts, char *ofile, float *y, float btime, float dt, char *var,
		       float cmpinc, float cmpaz);

receiverMode m_mode;
int m_nComp;

bool m_myPoint; // set to true if this processor writes to the arrays

std::string m_fileName, m_filePrefix;

double mX, mY, mZ, mGPX, mGPY, mGPZ; // original and actual location
bool m_zRelativeToTopography; // location is given relative to topography

int mWriteEvery;
bool m_usgsFormat, m_sacFormat;
string m_path;

// start time and time step 
double m_t0, m_dt;

// size of recording arrays
int mAllocatedSize;

// index of last recorded element
int mLastTimeStep;

// recording arrays
double** mRecordedSol;
float** mRecordedFloats;

// ignore this station if it is above the topography or outside the computational domain
bool mIgnore;

// sac header data
int mEventYear, mEventMonth, mEventDay, mEventHour, mEventMinute;
double mEventSecond, m_rec_lat, m_rec_lon, m_rec_gp_lat, m_rec_gp_lon;
double m_epi_lat, m_epi_lon, m_epi_depth, m_epi_time_offset, m_x_azimuth;

// sac ascii or binary?
bool mBinaryMode;

// UTC time for start of seismogram, 
//     m_t0 is start of seismogram in simulation time =  m_utc - utc reference time,
//           where utc reference time corresponds to simulation time zero.
bool m_utc_set;
int m_utc[7];

// Variables for rotating the output displacement or velocity components when Nort-East-UP is 
// selected (m_xyzcomponent=false) instead of Cartesian components (m_xyzcomponent=true).
   bool m_xyzcomponent;
   double m_calpha, m_salpha, m_thxnrm, m_thynrm;
//  double m_dthi, m_velocities;
// double m_dmx, m_dmy, m_dmz, m_d0x, m_d0y, m_d0z;
// double m_dmxy, m_dmxz, m_dmyz, m_d0xy, m_d0xz, m_d0yz;
};


#endif
