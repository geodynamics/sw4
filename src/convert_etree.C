#include <stdio.h>
#include <string>
#include <proj_api.h>

#include "cencalvm/storage/Payload.h"
#include "cencalvm/query/VMQuery.h"
#include "cencalvm/storage/ErrorHandler.h"
#include "cencalvm/storage/Geometry.h"

#define SQR(x) ((x)*(x))
int
main(int argc, char **argv) {
   
   int nprojpars = 6;
   double lon_orig = -123.0;
   double lat_orig = 35.0;
   double az_deg = 53.638, az;
// corner 0 and 3 have been corrected
   double dm_lon[]={-120.6440513713, -121.922036, -123.858493, -122.5623650935};
   double dm_lat[]={37.05006160038, 36.320331, 38.424179, 39.17450494695};
// azimuth of our x-axis has been corrected
   double alpha_deg=143.6380001671;
   
   double dm_x[4], dm_y[4], len[4], dx, dy, alpha, beta;

   projPJ pj_merc, pj_latlong;
   double x, y, xlon, ylat, x0, y0, xr, yr, xc, yc, xabs, yabs, xmap, ymap;
   int status, merc_len;
   char merc_def[256];
   
// hard coded for the USGS Bay Area model
   sprintf(merc_def, "+proj=tmerc +datum=NAD83 +units=m +lon_0=-123.0 +lat_0=35.0 +scale=0.9996");
   merc_len = strlen(merc_def);
   printf("Cartographic projection string: '%s' of length=%i\n", merc_def, merc_len);
   
   if (!(pj_merc = pj_init_plus(merc_def)) )
   {
      printf("Init of cartographic projection failed\n");
      exit(1);
   }
   
//   if (!(pj_latlong = pj_init_plus("+proj=latlong +ellps=WGS84")) )
   if (!(pj_latlong = pj_init_plus("+proj=latlong")) )
   {
      printf("Init of latlong projection failed\n");      
      exit(1);
   }
   
// compute offset in x and y for lon_0, lat_0
   x0 = lon_orig*DEG_TO_RAD;
   y0 = lat_orig*DEG_TO_RAD;
   status = pj_transform(pj_latlong, pj_merc, 1, 1, &x0, &y0, NULL );
   printf("Origin mapped from (lon,lat)=(%e, %e) to (x0,y0)=(%e, %e)\n", lon_orig, lat_orig, x0, y0);

// for all practical purposes, we can take (x0, y0)=(0,0)

   az = DEG_TO_RAD*az_deg;
   for (int q=0; q<4; q++)
   {
      xlon = dm_lon[q]*DEG_TO_RAD;
      ylat = dm_lat[q]*DEG_TO_RAD;
      status = pj_transform(pj_latlong, pj_merc, 1, 1, &xlon, &ylat, NULL );
      dm_x[q] = xlon-x0;
      dm_y[q] = ylat-y0;
   }
//  lenghts of sides
   for (int q=0; q<3; q++)
   {
      len[q]=sqrt(SQR(dm_x[q+1]-dm_x[q]) + SQR(dm_y[q+1]-dm_y[q]));
   }
   len[3]=sqrt(SQR(dm_x[0]-dm_x[3]) + SQR(dm_y[0]-dm_y[3]));
   
   for (int q=0; q<4; q++)
   {
      printf("Rel coord, corner %i: (x,y)=(%e, %e)\n", q, dm_x[q], dm_y[q]);
   }
   for (int q=0; q<4; q++)
   {
      printf("Edge %i length = %e\n", q, len[q]);
   }
// calculate azimuths
   dx = dm_x[0] - dm_x[1];
   dy = dm_y[0] - dm_y[1];
   beta = atan2(dy,dx);

   alpha = DEG_TO_RAD*alpha_deg;
   
   printf("Azimuth angle (x-axis) Alpha = %e\n", RAD_TO_DEG*alpha);
   
// tmp
//   exit(0);
   

// calculate (lon,lat) coordinates along cell centers in detailed model
// use SW corner as origin (dm_*[1])

   const char* mQueryKeys[]={"Density", "Vp", "Vs", "elevation", "Qp", "Qs","FaultBlock"};
   int mPayloadSize=7;
   double *mPayload, elev;
   cencalvm::query::VMQuery mQuery;
   cencalvm::storage::Geometry* mQueryGeom;
   char filename[]={"/Users/petersson1/USGSBayAreaVM-08.3.0.etree"}; // Hard-coded filename!
   
   mPayload = new double[7];
   
// open up the data base file
   mQuery.filename(filename);
   mQuery.queryType(cencalvm::query::VMQuery::MAXRES);
// Set values to be returned in queries
   mQuery.queryVals(mQueryKeys, mPayloadSize);
   mQuery.open();

// which subregion to query?
   int imax[]={2896, 1448, 724, 362}; // model defined for i=2897 on finest grid, but needs to be even for coarser models
   int jmax[]={1400,  700, 350, 175};

// cell size in horizontal directions
   double cl[] = {100.0, 200.0, 400.0, 800.0};

// which level to store? higher is coarser
   int lev=3;
   
// extent of model in etree database
   double xmin=0.0, xmax = 289715.88;
   double ymin=0.0, ymax = 140056.01;

// origin
   xc = 0.0;
   yc = 0.0;

// magic number
   int magic=1;
// float=4, double=8
   int precision=8; 

// how many blocks are stored?
   int nblocks=1;
// how many points in the vertical direction?
   int kmax[]={1, 1, 1, 1};
   
   FILE *fp=fopen("topography.dat","wb"); // Another hard-coded file name

// write header
   fwrite(&magic, sizeof(int), 1, fp);
   fwrite(&precision, sizeof(int), 1, fp);
// azimuth
   fwrite(&alpha_deg, sizeof(double), 1, fp);
// number of characters in projection string
   fwrite(&merc_len, sizeof(int), 1, fp);
// projection string
   fwrite(&merc_def, sizeof(char), merc_len, fp);

// number of blocks
   fwrite(&nblocks, sizeof(int), 1, fp);
// block sizes (base zero)
   fwrite(&(imax[lev]), sizeof(int), 1, fp);
   fwrite(&(jmax[lev]), sizeof(int), 1, fp);
   fwrite(&(kmax[lev]), sizeof(int), 1, fp);

// check material model on a coarse mesh along the surface
   int nQuery = 0;
   printf("Querying etree for topography at %i by %i points\n", imax[lev], jmax[lev]);
   

   for (int j=0; j<=jmax[lev]; j++)
//   int j = jmax; 
   {
      printf("j=%i out of %i\n", j, jmax[lev]);
      for (int i=0; i<=imax[lev]; i++)
//      int i=imax;
      
      {
         xr = xc + i*cl[lev];
         yr = yc + j*cl[lev];
      
// coordinate system centered at NE corner (#3)
         xabs = xmap = dm_x[3] + xr*sin(alpha) + yr*cos(alpha);
         yabs = ymap = dm_y[3] + xr*cos(alpha) - yr*sin(alpha);

// inverse projection to get (lon,lat)   
         status = pj_transform(pj_merc, pj_latlong, 1, 1, &xmap, &ymap, NULL );

         xlon = xmap*RAD_TO_DEG;
         ylat = ymap*RAD_TO_DEG;
   
//      printf("cell center (x,y)=(%e, %e), (lon, lat)=(%e, %e)\n", xabs, yabs, xlon, ylat);

// query the detailed material model at depth 12.5 m to avoid the air/water interface
         elev = -12.5;
         mQuery.query(&mPayload, mPayloadSize, xlon, ylat, elev);
         nQuery++;
         
// Make sure the query didn't generated a warning or error
         if (mQuery.errorHandler()->status() != cencalvm::storage::ErrorHandler::OK) 
         {
            printf("Something went wrong for rel. coords (xr,yr)=(%e, %e)\n", xr, yr);
// reset status for next query
	    mQuery.errorHandler()->resetStatus();
// phony value
            mPayload[3] = -999.999;
         }
         fwrite(&(mPayload[3]), sizeof(double), 1, fp);
      }
   }
   
// done
   fclose(fp);
   mQuery.close();

   printf("Called query() %i times\n", nQuery);
   exit(0);
}
