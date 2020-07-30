/*-------------------------------------------------------------------------------------------------
 * Forward and Adjoint Simulation via Elastic Staggered grid Finite Difference
 *
 * Frechet Kernels are computed via Adjoint State Method
 *
 * Updating flow:
 * Forward modeling:                 compute v (t) -> compute sigma (t) -> add source on sigma (t)
 * Source wavefield reconstruction:  substract source from sigma (t + 1) -> compute sigma (t) -> compute v (t)  (Exact time reversal)
 *
 * Backward propagation:             compute v (t) -> Add source on v (t) -> compute sigma (t)
 * Born modeling:                    compute v (t) -> Add source on v (t) -> compute sigma (t)
 --------------------------------------------------------------------------------------------------*/

#include "efwi.h"

void traveltime_simulation(int VERBOSE, int FREE_SURFACE, int SOURCE_TYPE,
                          float **Vp, float **Vs,
                          float dx, float dz, double ox, double oz, float dt, int NX, int NZ, int NXP, int NZP,
                          int nx, int nz, int npad, int npml, int nzshift, int nt,
                          float *w, int isx, int isz, int **ireceivers, int nreceivers, float sfreq,
                          int ishot, int USE_WINDOW, char *WINDOW_TYPE, char *tag, int myid)
{

    int i, ig, j, b1, b2, b3, irx, irz;
    int shiftx, shiftz;
    int *in;
    float s1, s2, t0, t1, t2, t3;
    float smin, smax;
    FILE *fp, *ft;


    float *timep, *times, *sp, *ss;
    bool plane[3];

    plane[0]=false;
    plane[1]=false;
    plane[2]=false;
    b1=1;
    b2=1;
    b3=1;

    timep=alloc1float(nz*nx);  // nz fastest
    times=alloc1float(nz*nx);
    sp=alloc1float(nz*nx);
    ss=alloc1float(nz*nx);
    in=alloc1int(nz*nx);


    shiftx = npad+npml;
    shiftz = nzshift;

    s1 = (isz-shiftz)*dz;
    s2 = (isx-shiftx)*dx;


    fprintf(stderr, "fast marching sx=%g sz=%g nx=%d nz=%d NX=%d NZ=%d WINDOW_TYPE=%c\n", s2, s1, nx, nz, NX, NZ, *WINDOW_TYPE);
        smin = 1e20;
        smax = -1e20;

        for(i=0; i<nx; i++) {
            for(j=0; j<nz; j++) {
                sp[i*nz+j] = 1./(Vp[i+shiftx][j+shiftz]*Vp[i+shiftx][j+shiftz]);  // slowness^2
                if(Vp[i+shiftx][j+shiftz] < smin) smin = Vp[i+shiftx][j+shiftz];
                if(Vp[i+shiftx][j+shiftz] > smax) smax = Vp[j+shiftx][i+shiftz];
            }
        }


         for(i=0; i<nx; i++) {
             for(j=0; j<nz; j++) {
                 ss[i*nz+j] = 1./(Vs[i+shiftx][j+shiftz]*Vs[i+shiftx][j+shiftz]);  // slowness^2
                 if(Vs[i+shiftx][j+shiftz] < smin) smin = Vs[i+shiftx][j+shiftz];
                 if(Vs[i+shiftx][j+shiftz] > smax) smax = Vs[j+shiftx][i+shiftz];
             }
         }

       //fprintf(stderr, "vmin=%g vmax=%g\n", smin, smax);


       if(*WINDOW_TYPE=='P' || *WINDOW_TYPE=='B') {
       fastmarch_init(1,nx,nz);
       fastmarch(timep              /* time */,
                sp                  /* slowness squared */,
                in                  /* in/front/out flag */,
                plane               /* if plane source */,
                1,nx,nz             /* dimensions */,
                0.0,ox,oz           /* origin */,
                1.0,dx,dz           /* sampling */,
                0.0,s2,s1           /* source */,
                b3,b2,b1            /* box around the source */,
                2                   /* accuracy order (1,2,3) */);
       fastmarch_close();
       }

       if(*WINDOW_TYPE=='S' || *WINDOW_TYPE=='B') {
       fastmarch_init(1,nx,nz);

       fastmarch(times              /* time */,
                ss                  /* slowness squared */,
                in                  /* in/front/out flag */,
                plane               /* if plane source */,
                1,nx,nz             /* dimensions */,
                0.0,ox,oz           /* origin */,
                1.0,dx,dz           /* sampling */,
                0.0,s2,s1           /* source */,
                b3,b2,b1            /* box around the source */,
                2                   /* accuracy order (1,2,3) */);
       fastmarch_close();

       }



       //if(myid==0) {
       //fp=fopen("tt.bin", "wb");
       // fwrite(time, sizeof(float), nx*nz, fp);
       //fclose(fp);
       //}

       char filename[STRINGSIZE];
       if(USE_WINDOW) {
       if(strcmp(tag, "obs")==0) sprintf(filename, "./input/window_%s_%d.txt", tag, ishot);
       else sprintf(filename, "window_%s_%d.txt", tag, ishot);

       ft=fopen(filename, "w");
       if(ft==NULL) {
           fprintf(stderr, "unable to open file=%s\n", filename);
           exit(-1);
       }

       for (ig = 0; ig < nreceivers; ig++)
       {

           irx = ireceivers[ig][0];
           irz = ireceivers[ig][1];

           if(*WINDOW_TYPE=='P') {
           t0 = timep[(irx-shiftx)*nz+irz-shiftz];
           t1 = t0 + 1./(sfreq*0.50);

           fprintf(ft, "%g\t%g\n", t0, t1);
           }
           else if(*WINDOW_TYPE=='S') {
               t0 = times[(irx-shiftx)*nz+irz-shiftz];
               t1 = t0 + 1./(sfreq*0.50);

               fprintf(ft, "%g\t%g\n", t0, t1);
           }
           else if(*WINDOW_TYPE=='B') {
               t0 = timep[(irx-shiftx)*nz+irz-shiftz];
               t1 = t0 + 1./(sfreq*0.50);
               t2 = times[(irx-shiftx)*nz+irz-shiftz];
               t3 = t2 + 1./(sfreq*0.50);

               fprintf(ft, "%g\t%g\t\t%g\t%g\n", t0, t1, t2, t3);
           }
        } // ig
       fclose(ft);
       } // WINDOW


   if(timep) free(timep);
   if(times) free(times);
   if(sp) free(sp);
   if(ss) free(ss);
   if(in) free(in);


}//End function
