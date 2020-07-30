//  solving trave time given a velocity model
//  Wei Liu
//  May 2020
//

#include "EW.h"
#include "impose_cartesian_bc.h"
#include "cf_interface.h"
#include "f_interface.h"

#ifdef USE_HDF5
#include "sachdf5.h"
#endif

#define SQR(x) ((x)*(x))

//--------------------------------------------------------------------
void EW::solveTT( vector<Source*> & a_Sources, vector<TimeSeries*> & a_TimeSeries,
		double* xs, int nmpars, int event, int myRank)
{

// Set the number of time steps, allocate the recording arrays, and set reference time in all time series objects  
//#pragma omp parallel for



  for (int ts=0; ts<a_TimeSeries.size(); ts++)
  {
     a_TimeSeries[ts]->allocateRecordingArrays( mNumberOfTimeSteps[event]+1, mTstart, mDt); // AP: added one to mNumber...
     // In forward solve, the output receivers will use the same UTC as the
     // global reference utc0, therefore, set station utc equal reference utc.
     //     if( m_utc0set )
     //	a_TimeSeries[ts]->set_station_utc( m_utc0 );
  }
  // the Source objects get discretized into GridPointSource objects
  vector<GridPointSource*> point_sources;

// Transfer source terms to each individual grid as point sources at grid points.
  for( unsigned int i=0 ; i < a_Sources.size() ; i++ )
      a_Sources[i]->set_grid_point_sources4( this, point_sources );


std::cout << "source x0=" << a_Sources[0]->getX0() << " y0=" << a_Sources[0]->getY0() << " z0=" << a_Sources[0]->getZ0() << std::endl;

float_sw4* cp =(float_sw4*)malloc(nmpars/3*sizeof(float_sw4));
float_sw4* cs =(float_sw4*)malloc(nmpars/3*sizeof(float_sw4));


float_sw4 cpmin, cpmax;
float_sw4 csmin, csmax;

cpmin=1e20; cpmax=-1e20;
csmin=1e20; csmax=-1e20;

//#pragma omp parallel for
	 for( size_t i = 0 ; i < nmpars/3 ; i++ )
	 {
	       cs[i] = xs[3*i+1];
          cp[i] = xs[3*i+2];

          if(cp[i]<cpmin) cpmin=cp[i];
          if(cp[i]>cpmax) cpmax=cp[i];
          if(cs[i]<csmin) csmin=cs[i];
          if(cs[i]>csmax) csmax=cs[i];
	 }

std::cout << "solveTT cpmin=" << cpmin << " cpmax=" << cpmax << std::endl;
std::cout << "solveTT csmin=" << csmin << " csmax=" << csmax << std::endl;

int ntr = a_TimeSeries.size();
std::cout << "number of traces=" << ntr << std::endl;

// each rank has the complete trace info
if(myRank==0) {
for(int ig=0; ig<ntr; ig++) {
   std::cout << "trace ig=" << ig << " x=" << a_TimeSeries[ig]->getX() << " y=" << a_TimeSeries[ig]->getY() <<  std::endl;
}}

free(cp);
free(cs);


   MPI_Abort(MPI_COMM_WORLD,0);
}