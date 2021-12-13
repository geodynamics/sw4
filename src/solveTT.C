//  solving trave time given a velocity model
//  Wei Liu
//  July 2020
//

#include "EW.h"
#include "impose_cartesian_bc.h"
#include "cf_interface.h"
#include "f_interface.h"
#include "MaterialParameterization.h"

#define STRINGSIZE 128

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


//--------------------------------------------------------------------
void EW::solveTT( Source* a_Source, vector<TimeSeries*> & a_TimeSeries,
		double* xs, int nmpars, MaterialParameterization* mp, int wave_mode, float twinshift, float twinscale, float freq, int event, int myrank)
{
   int verbose=0;
   int ix, iy, iz, nx, ny, nz, n;
   nx = mp->getNX();
   ny = mp->getNY();
   nz = mp->getNZ();

   if( verbose >= 1 && myrank == 0 )
      std::cout << "source x0=" << a_Source->getX0() << " y0=" << a_Source->getY0() << " z0=" << a_Source->getZ0() << std::endl;
   n= nmpars/nx/ny/nz;

   if( verbose >= 1 && myrank == 0 )
   {
      std::cout << "mp->xmin=" << mp->getXmin() << " ymin=" << mp->getYmin() << " zmin=" << mp->getZmin() << std::endl;
   std::cout << "nx=" << nx << " ny=" << ny << " nz=" << nz << " n=" << n << std::endl;
   }

float *cp =(float*)malloc(nmpars/n*sizeof(float));
float *cs =(float*)malloc(nmpars/n*sizeof(float));

float *sp =(float*)malloc(nmpars/n*sizeof(float));
float *ss =(float*)malloc(nmpars/n*sizeof(float));
float* timep =(float*)malloc(nmpars/n*sizeof(float));
float* times =(float*)malloc(nmpars/n*sizeof(float));
int *inflag =(int*)malloc(nmpars/n*sizeof(int));


float_sw4 cpmin, cpmax;
float_sw4 csmin, csmax;

cpmin=1e20; cpmax=-1e20;
csmin=1e20; csmax=-1e20;


//#pragma omp parallel for
	 for( size_t i = 0 ; i < nmpars/n ; i++ )
	 {
	      cs[i] = xs[n*i+n-2]; 
          cp[i] = xs[n*i+n-1];
          
          ss[i] = 1./xs[n*i+n-2];
          ss[i] = ss[i]*ss[i];
          sp[i] = 1./xs[n*i+n-1];
          sp[i] = sp[i]*sp[i];
          
          if(cp[i]<cpmin) cpmin=cp[i];
          if(cp[i]>cpmax) cpmax=cp[i];
          if(cs[i]<csmin) csmin=cs[i];
          if(cs[i]>csmax) csmax=cs[i];
	 }

         if( verbose >= 1 && myrank == 0 )
         {
            std::cout << "solveTT cpmin=" << cpmin << " cpmax=" << cpmax << std::endl;
            std::cout << "solveTT csmin=" << csmin << " csmax=" << csmax << std::endl;
         }

int ntr = a_TimeSeries.size();
//std::cout << "nmpars=" << nmpars << " number of traces=" << ntr << std::endl;


//std::cout << "mp xmin=" << mp->getXmin() << " hx=" << mp->getDx() << " nx=" << mp->getNX() << std::endl;
    bool plane[3]={false,false,false};

       fastmarch_init(mp->getNZ(),mp->getNY(),mp->getNX());     // nx fast 
       fastmarch(timep,               // time 
                sp,                   // slowness squared
                inflag,               // in/front/out flag
                plane,                // if plane source 
                nz, ny, nx,           // dimensions
                mp->getZmin(),mp->getYmin(),mp->getXmin(),     // origin
                mp->getDz(),mp->getDy(),mp->getDx(),           // sampling
                a_Source->getZ0(),a_Source->getY0(),a_Source->getX0(),     // source
                1,1,1,                // box around the source
                2);                   // accuracy order (1,2,3) 
         
         fastmarch(times,               // time 
                ss,                   // slowness squared
                inflag,               // in/front/out flag
                plane,                // if plane source 
                nz, ny, nx,           // dimensions
                mp->getZmin(),mp->getYmin(),mp->getXmin(),     // origin
                mp->getDz(),mp->getDy(),mp->getDx(),           // sampling
                a_Source->getZ0(),a_Source->getY0(),a_Source->getX0(),     // source
                1,1,1,                // box around the source
                2);                   // accuracy order (1,2,3) 

       fastmarch_close();

   //checkMinMax(nmpars/3, timep, "timep");mopt->m_misfit == Mopt::CROSSCORR
   // output to a file for QC
   char file[STRINGSIZE];
   FILE *fd;

   if(myrank==0) {     
      sprintf(file, "%s/time_event_%d.txt", a_TimeSeries[0]->getPath().c_str(),event);
   fd = fopen(file, "w");
   }
   

    float_sw4 win;

    win = 1*(a_Source->getFrequency()<freq? 1.0/a_Source->getFrequency() : 1.0/freq);  // one cycle of peak frequency
    //cout << "original win=" << win << endl;
    win *= twinscale;  // expand or shrink 
    //cout << "dialated win=" << win << endl;

   for(int ig=0; ig<ntr; ig++) {
   //
    ix = (a_TimeSeries[ig]->getX() - mp->getXmin())/mp->getDx()+0.5;
    iy = (a_TimeSeries[ig]->getY() - mp->getYmin())/mp->getDy()+0.5;
    iz = (a_TimeSeries[ig]->getZ() - mp->getZmin())/mp->getDz()+0.5;
    //std::cout << "trace ig=" << ig+1 << " x=" << a_TimeSeries[ig]->getX() << " y=" << a_TimeSeries[ig]->getY() <<  " tp=" 
    //   << timep[iz*nx*ny+iy*nx+ix] << " ts=" << times[iz*nx*ny+iy*nx+ix] << std::endl;

   // set window 
   switch(wave_mode) {
   case 0:  // P-wave only
        a_TimeSeries[ig]->set_window(timep[iz*nx*ny+iy*nx+ix]+win*twinshift, timep[iz*nx*ny+iy*nx+ix]+win, 
       timep[iz*nx*ny+iy*nx+ix]+win*twinshift, timep[iz*nx*ny+iy*nx+ix]+win);

       if(myrank==0) fprintf(fd, "%d   %d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
               event, ig+1, a_TimeSeries[ig]->getStationName().c_str(), a_TimeSeries[ig]->getX(),a_TimeSeries[ig]->getY(),a_TimeSeries[ig]->getZ(),
               timep[iz*nx*ny+iy*nx+ix]+win*twinshift, timep[iz*nx*ny+iy*nx+ix]+win, 
               timep[iz*nx*ny+iy*nx+ix]+win*twinshift, timep[iz*nx*ny+iy*nx+ix]+win);
         break;
   case 1:  // S-wave only
       a_TimeSeries[ig]->set_window(times[iz*nx*ny+iy*nx+ix]+win*twinshift, times[iz*nx*ny+iy*nx+ix]+win, 
       times[iz*nx*ny+iy*nx+ix]+win*twinshift, times[iz*nx*ny+iy*nx+ix]+win);

       if(myrank==0) fprintf(fd, "%d   %d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
               event, ig+1, a_TimeSeries[ig]->getStationName().c_str(), a_TimeSeries[ig]->getX(),a_TimeSeries[ig]->getY(),a_TimeSeries[ig]->getZ(),
               times[iz*nx*ny+iy*nx+ix]+win*twinshift, times[iz*nx*ny+iy*nx+ix]+win, 
               times[iz*nx*ny+iy*nx+ix]+win*twinshift, times[iz*nx*ny+iy*nx+ix]+win);

         break;
   default:
       a_TimeSeries[ig]->set_window(timep[iz*nx*ny+iy*nx+ix]+win*twinshift, timep[iz*nx*ny+iy*nx+ix]+win, 
       times[iz*nx*ny+iy*nx+ix]+win*twinshift, times[iz*nx*ny+iy*nx+ix]+win);

       if(myrank==0) fprintf(fd, "%d   %d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
               event, ig+1, a_TimeSeries[ig]->getStationName().c_str(), a_TimeSeries[ig]->getX(),a_TimeSeries[ig]->getY(),a_TimeSeries[ig]->getZ(),
               timep[iz*nx*ny+iy*nx+ix]+win*twinshift, timep[iz*nx*ny+iy*nx+ix]+win, 
               times[iz*nx*ny+iy*nx+ix]+win*twinshift, times[iz*nx*ny+iy*nx+ix]+win);
     }
   }  // loop over traces

   if(myrank==0) fclose(fd);


   if(timep) free(timep);
   if(times) free(times);
   if(cp) free(cp);
   if(cs) free(cs);
   if(sp) free(sp);
   if(ss) free(ss);
   if(inflag) free(inflag);

   if( verbose >= 1 && myrank == 0 )
      std::cout << "solveTT completed" << std::endl;
 
}
