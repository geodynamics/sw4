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
   double dm_lon[]={-120.644051, -121.922036, -123.858493, -122.562365};
   double dm_lat[]={37.050062, 36.320331, 38.424179, 39.174505};
   double dm_x[4], dm_y[4], len[4], dx, dy, beta;
   
   projPJ pj_merc, pj_latlong;
   double x, y, xlon, ylat, x0, y0, cl, xr, yr, xc, yc, xabs, yabs, xmap, ymap;
   int status;
   char merc_def[256];
   
   sprintf(merc_def, "+proj=tmerc +datum=NAD83 +units=m +lon_0=-123.0 +lat_0=35.0 +scale=0.9996");
// sprintf(merc_def, "+proj=utm +ellps=WGS84 +units=m +lon_0=-123.0 +lat_0=35.0");
//   sprintf(merc_def, "+proj=utm +ellps=WGS84 +units=m +lon_0=0.0 +lat_0=0.0");
   printf("Mercator string: '%s'\n", merc_def);
   
   if (!(pj_merc = pj_init_plus(merc_def)) )
   {
      printf("Init of mercator projection failed\n");
      exit(1);
   }
   
   if (!(pj_latlong = pj_init_plus("+proj=latlong +ellps=WGS84")) )
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
   
   printf("Beta = %e, Azimuth (90-beta) = %e\n", RAD_TO_DEG*beta, 90-RAD_TO_DEG*beta);
   
// calculate (lon,lat) coordinates along cell centers in detailed model
// use SW corner as origin (dm_*[1])

   const char* mQueryKeys[]={"Density", "Vp", "Vs", "elevation", "Qp", "Qs","FaultBlock"};
   int mPayloadSize=7;
   double *mPayload, elev;
   cencalvm::query::VMQuery mQuery;
   cencalvm::storage::Geometry* mQueryGeom;
   char filename[]={"/Users/petersson1/USGSBayAreaVM-08.3.0.etree"}; // Rather hard-coded
   
   mPayload = new double[7];
   
// open up the data base file
   mQuery.filename(filename);
   mQuery.queryType(cencalvm::query::VMQuery::MAXRES);
// Set values to be returned in queries
   mQuery.queryVals(mQueryKeys, mPayloadSize);
   mQuery.open();

// cell size
   cl = 100.0;
   
   FILE *fp=fopen("elev.dat","w");

// cell center near SW corner
   xr = 0.5*cl;
   yr = 0.5*cl;

// check material model on a coarse mesh along the surface
   for (int j=0; j<2897; j+=4)
   {
      for (int i=0; i<1400; i+=4)
      {
         xr = 0.5*cl + i*cl;
         yr = 0.5*cl + j*cl;
      
         xc = xr*cos(beta) - yr*sin(beta);
         yc = xr*sin(beta) + yr*cos(beta);
   
         xabs = xmap = dm_x[1] + xc;
         yabs = ymap = dm_y[1] + yc;

// inverse projection to get (lon,lat)   
         status = pj_transform(pj_merc, pj_latlong, 1, 1, &xmap, &ymap, NULL );

         xlon = xmap*RAD_TO_DEG;
         ylat = ymap*RAD_TO_DEG;
   
//      printf("cell center (x,y)=(%e, %e), (lon, lat)=(%e, %e)\n", xabs, yabs, xlon, ylat);

// query the detailed material model at depth 12.5 m to avoid the air/water interface
         elev = -12.5;
         mQuery.query(&mPayload, mPayloadSize, xlon, ylat, elev);
                  
// Make sure the query didn't generated a warning or error
         if (mQuery.errorHandler()->status() != cencalvm::storage::ErrorHandler::OK) 
         {
            printf("Something went wrong for rel. coords (xr,yr)=(%e, %e)\n", xr, yr);
         }
         else
         {
            fprintf(fp, "%e %e %e\n", xlon, ylat, mPayload[3]);
            // fprintf(fp, "%e %e ", xr, yr);
            // for (int q=0; q<mPayloadSize; q++)
            //    fprintf(fp, "%e ", mPayload[q]);
            // fprintf(fp, "\n");
         }
      }
   }
   
// done
   fclose(fp);
   
   mQuery.close();
   exit(0);
}
