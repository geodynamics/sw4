//  solving trave time given a velocity model
//  Wei Liu
//  July 2020
//

#include "EW.h"
#include "impose_cartesian_bc.h"
#include "cf_interface.h"
#include "f_interface.h"

#ifdef USE_HDF5
#include "sachdf5.h"
#endif

#define SQR(x) ((x)*(x))

// traveltime estimation for time windowing
void fastmarch_init (int n3,int n2,int n1);

void fastmarch (float* time                /* time */,
        float* v                   /* slowness squared */,
        int* in                    /* in/front/out flag */,
        bool* plane                /* if plane source */,
        int   n3,  int n2,  int n1 /* dimensions */,
        float o3,float o2,float o1 /* origin */,
        float d3,float d2,float d1 /* sampling */,
        float s3,float s2,float s1 /* source */,
        int   b3,  int b2,  int b1 /* box around the source */,
        int order                  /* accuracy order (1,2,3) */);

void fastmarch_close (void);
void checkMinMax(float* data, int n);

//--------------------------------------------------------------------
void EW::solveTT( vector<Source*> & a_Sources, vector<TimeSeries*> & a_TimeSeries,
		double* xs, int nmpars, MaterialParameterization* mp, int event)
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

float *cp =(float*)malloc(nmpars/3*sizeof(float));
float *cs =(float*)malloc(nmpars/3*sizeof(float));

float *sp =(float*)malloc(nmpars/3*sizeof(float));
float *ss =(float*)malloc(nmpars/3*sizeof(float));
float* timep =(float*)malloc(nmpars/3*sizeof(float));
float* times =(float*)malloc(nmpars/3*sizeof(float));
int *inflag =(int*)malloc(nmpars/3*sizeof(int));

float_sw4 cpmin, cpmax;
float_sw4 csmin, csmax;

cpmin=1e20; cpmax=-1e20;
csmin=1e20; csmax=-1e20;

//#pragma omp parallel for
	 for( size_t i = 0 ; i < nmpars/3 ; i++ )
	 {
	       cs[i] = xs[3*i+1]; 
          cp[i] = xs[3*i+2];

          ss[i] = 1./xs[3*i+1];
          ss[i] = ss[i]*ss[i];
          sp[i] = 1./xs[3*i+2];
          sp[i] = sp[i]*sp[i];
          
          if(cp[i]<cpmin) cpmin=cp[i];
          if(cp[i]>cpmax) cpmax=cp[i];
          if(cs[i]<csmin) csmin=cs[i];
          if(cs[i]>csmax) csmax=cs[i];
	 }

std::cout << "solveTT cpmin=" << cpmin << " cpmax=" << cpmax << std::endl;
std::cout << "solveTT csmin=" << csmin << " csmax=" << csmax << std::endl;

int ntr = a_TimeSeries.size();
std::cout << "nmpars=" << nmpars << " number of traces=" << ntr << std::endl;
   int ix, iy, iz, nx, ny, nz;
   nx = mp->getNX();
   ny = mp->getNY();
   nz = mp->getNZ();

std::cout << "mp xmin=" << mp->getXmin() << " hx=" << mp->getDx() << " nx=" << mp->getNX() << std::endl;
    bool plane[3]={false,false,false};

       fastmarch_init(mp->getNZ(),mp->getNY(),mp->getNX());     // nx fast 
       fastmarch(timep,               // time 
                sp,                   // slowness squared
                inflag,               // in/front/out flag
                plane,                // if plane source 
                nz, ny, nx,           // dimensions
                mp->getZmin(),mp->getYmin(),mp->getXmin(),     // origin
                mp->getDz(),mp->getDy(),mp->getDx(),           // sampling
                a_Sources[0]->getZ0(),a_Sources[0]->getY0(),a_Sources[0]->getX0(),     // source
                1,1,1,                // box around the source
                2);                   // accuracy order (1,2,3) 
         
         fastmarch(times,               // time 
                ss,                   // slowness squared
                inflag,               // in/front/out flag
                plane,                // if plane source 
                nz, ny, nx,           // dimensions
                mp->getZmin(),mp->getYmin(),mp->getXmin(),     // origin
                mp->getDz(),mp->getDy(),mp->getDx(),           // sampling
                a_Sources[0]->getZ0(),a_Sources[0]->getY0(),a_Sources[0]->getX0(),     // source
                1,1,1,                // box around the source
                2);                   // accuracy order (1,2,3) 


       fastmarch_close();

   //checkMinMax(timep, nmpars/3);

   for(int ig=0; ig<ntr; ig++) {
   //
    ix = (a_TimeSeries[ig]->getX() - mp->getXmin())/mp->getDx()+0.5;
    iy = (a_TimeSeries[ig]->getY() - mp->getYmin())/mp->getDy()+0.5;
    iz = (a_TimeSeries[ig]->getZ() - mp->getZmin())/mp->getDz()+0.5;
    std::cout << "trace ig=" << ig << " x=" << a_TimeSeries[ig]->getX() << " y=" << a_TimeSeries[ig]->getY() <<  " tp=" 
       << timep[iz*nx*ny+iy*nx+ix] << " ts=" << times[iz*nx*ny+iy*nx+ix] << std::endl;
   } 

   if(timep) free(timep);
   if(times) free(times);
   if(sp) free(sp);
   if(ss) free(ss);
   if(inflag) free(inflag);

std::cout << "solveTT completed" << std::endl;
   
   MPI_Abort(MPI_COMM_WORLD,0);
}

void checkMinMax(float* data, int n)
{
   float min=1e20;
   float max=-1e20;
   for(int i=0; i<n; i++) {
      if(data[i]<min) min=data[i];
      if(data[i]>max) max=data[i];
   }
   std::cout << "min=" << min << " max=" << max << std::endl;
}