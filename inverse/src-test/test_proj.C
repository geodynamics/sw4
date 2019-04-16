//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
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
   double x, y, xlon, ylat, x0, y0, cl, xr, yr, xc, yc, xabs, yabs, xmap, ymap;
   int status, q;
   char merc_def[256];
   
// hard coded for the USGS Bay Area model
   sprintf(merc_def, "+proj=tmerc +datum=NAD83 +units=m +lon_0=-123.0 +lat_0=35.0 +scale=0.9996");
   printf("Mercator string: '%s'\n", merc_def);
   
   if (!(pj_merc = pj_init_plus(merc_def)) )
   {
      printf("Init of mercator projection failed\n");
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
   for (q=0; q<4; q++)
   {
      xlon = dm_lon[q]*DEG_TO_RAD;
      ylat = dm_lat[q]*DEG_TO_RAD;
      status = pj_transform(pj_latlong, pj_merc, 1, 1, &xlon, &ylat, NULL );
      dm_x[q] = xlon-x0;
      dm_y[q] = ylat-y0;
   }
//  lenghts of sides
   for (q=0; q<3; q++)
   {
      len[q]=sqrt(SQR(dm_x[q+1]-dm_x[q]) + SQR(dm_y[q+1]-dm_y[q]));
   }
   len[3]=sqrt(SQR(dm_x[0]-dm_x[3]) + SQR(dm_y[0]-dm_y[3]));
   
   for (q=0; q<4; q++)
   {
      printf("Rel coord, corner %i: (x,y)=(%e, %e)\n", q, dm_x[q], dm_y[q]);
   }
   for (q=0; q<4; q++)
   {
      printf("Edge %i length = %e\n", q, len[q]);
   }
// calculate azimuths
   dx = dm_x[0] - dm_x[1];
   dy = dm_y[0] - dm_y[1];
   beta = atan2(dy,dx);

   alpha = DEG_TO_RAD*alpha_deg;
   
   printf("Beta(c1-c0) = %e, Alpha (c3-c0) = %e\n", RAD_TO_DEG*beta, RAD_TO_DEG*alpha);
   
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

// cell size in horizontal directions
   cl = 100.0;
   
// test for location of SE corner of etree model   
//    xr = 289715.88;
//    yr = 0.0;

// // coordinate system centered at NE corner (#3)
//    xabs = xmap = dm_x[3] + xr*sin(alpha) + yr*cos(alpha);
//    yabs = ymap = dm_y[3] + xr*cos(alpha) - yr*sin(alpha);

// // inverse projection to get (lon,lat)   
//    status = pj_transform(pj_merc, pj_latlong, 1, 1, &xmap, &ymap, NULL );

//    xlon = xmap*RAD_TO_DEG;
//    ylat = ymap*RAD_TO_DEG;

//    printf("Origin lon = %20.12e, SE-corner lon = %20.12e\n", xlon, dm_lon[0]);
//    printf("Origin lat = %20.12e, SE-corner lat = %20.12e\n", ylat, dm_lat[0]);

// // correct the (lat,lon) coordinates of the SE corner
//    dm_lon[0] = xlon;
//    dm_lat[0] = ylat;
// // re-compute the projection of the SE corner
//    q = 0;
//    xlon = dm_lon[q]*DEG_TO_RAD;
//    ylat = dm_lat[q]*DEG_TO_RAD;
//    status = pj_transform(pj_latlong, pj_merc, 1, 1, &xlon, &ylat, NULL );
//    dm_x[q] = xlon-x0;
//    dm_y[q] = ylat-y0;

// // Correct the azimuth of the x-axis...
//    double az2 = 0.5*M_PI + atan2(dm_y[3]-dm_y[0], dm_x[0]-dm_x[3]);
   
//    printf("Alpha = %20.12e, Alpha-corrected = %20.12e\n", RAD_TO_DEG*alpha, RAD_TO_DEG*az2);
   
   int i,j, k;
   
   int nQuery = 0;
// test for location of Sw corner of etree model   
//    xc = 0.0;
//    yc = 140056.01;
// // step size
//    cl = 1e-2;
//    for (j=-10; j<=10; j++)
//       for (i=-10; i<=10; i++)
//       {
//          xr = xc + cl*i;
//          yr = yc + cl*j;
// // coordinate system centered at NE corner (#3)
//          xabs = xmap = dm_x[3] + xr*sin(alpha) + yr*cos(alpha);
//          yabs = ymap = dm_y[3] + xr*cos(alpha) - yr*sin(alpha);

// // inverse projection to get (lon,lat)   
//          status = pj_transform(pj_merc, pj_latlong, 1, 1, &xmap, &ymap, NULL );

//          xlon = xmap*RAD_TO_DEG;
//          ylat = ymap*RAD_TO_DEG;
   
// //      printf("cell center (x,y)=(%e, %e), (lon, lat)=(%e, %e)\n", xabs, yabs, xlon, ylat);

// // query the detailed material model at depth 12.5 m to avoid the air/water interface
//          elev = -12.5;
//          mQuery.query(&mPayload, mPayloadSize, xlon, ylat, elev);
//          nQuery++;
         
// // Make sure the query didn't generated a warning or error
//          if (mQuery.errorHandler()->status() != cencalvm::storage::ErrorHandler::OK) 
//          {
//             printf("Outside point with rel. index (i,j)=(%i, %i)\n", i, j);
// // reset status for next query
// 	    mQuery.errorHandler()->resetStatus();
//          }
//          else
//          {
//             printf("Inside point with rel. index (i,j)=(%i, %i)\n", i, j);
//          }
//       }
//    printf("Step size cl=%e\n", cl);
   
//    exit(0);

// check material model on a coarse mesh along the surface
   nQuery = 0;

   cl = 800.0; // cell size in horizontal directions
//   cl = 100.0;
   
//   double clv = 800.0; // cell size in vertical dir
   double clv = 200.0; // cell size in vertical dir

   double xmax = 289715.88;
   double ymax = 140056.01;
   double zmax = 45000.0;

// which subregion to query?
   int imax=2897;
   int jmax=1400;

// bottom block
//   int kmax=56;
// top block
   int kmax=11;

//   FILE *fp=fopen("elev-x1.dat","w");

// origin
   xc = 0;
   yc = 0;

   double z;
// bottom block
//   double z0=412.5;
// top block
   double z0=-1587.5;
   
//   for (j=0; j<=jmax; j+=1)
//   j = jmax; 
   j = 175; // should be over the ocean
//   j=0;
   {
//      for (i=0; i<=imax; i+=1)
//      i=imax;
//      i=0;    
      i=30;
      {
//         xr = 50.0 + i*cl; // cell center ?
//         yr = 50.0 + j*cl;
         xr = xc + i*cl; // full extent of model (0.1 is necessary to add!)
         yr = yc + j*cl;
      
// coordinate system centered at SW corner (#1)
// xabs = xmap = dm_x[1] + xr*cos(beta) - yr*sin(beta);
// yabs = ymap = dm_y[1] + xr*sin(beta) + yr*cos(beta);

// coordinate system centered at NE corner (#3)
         xabs = xmap = dm_x[3] + xr*sin(alpha) + yr*cos(alpha);
         yabs = ymap = dm_y[3] + xr*cos(alpha) - yr*sin(alpha);

// inverse projection to get (lon,lat)   
         status = pj_transform(pj_merc, pj_latlong, 1, 1, &xmap, &ymap, NULL );

         xlon = xmap*RAD_TO_DEG;
         ylat = ymap*RAD_TO_DEG;
   
//      printf("cell center (x,y)=(%e, %e), (lon, lat)=(%e, %e)\n", xabs, yabs, xlon, ylat);

// query the detailed material model at depth 12.5 m to avoid the air/water interface
         for (k=0; k<kmax; k++)
         {
//            elev = -12.5;
            z = z0 + k*clv;
            elev = -z;

            mQuery.query(&mPayload, mPayloadSize, xlon, ylat, elev);
            nQuery++;
         
// Make sure the query didn't generated a warning or error
            if (mQuery.errorHandler()->status() == cencalvm::storage::ErrorHandler::ERROR) 
            {
               printf("Fatal query error for point (ix, jy, kz, elev)=(%i, %i, %i, %e)\n", i, j, k, elev);
// reset status for next query
               mQuery.errorHandler()->resetStatus();
            }
            else
            {
               printf("index=(%i %i %i), elev=%e, dens=%e, topo=%e, status=%i\n", i, j, k, elev, 
                      mPayload[0], mPayload[3], mQuery.errorHandler()->status()); // payload[0] is density
// reset status for next query
               mQuery.errorHandler()->resetStatus();
            // fprintf(fp, "%e %e ", xr, yr);
            // for (int q=0; q<mPayloadSize; q++)
            //    fprintf(fp, "%e ", mPayload[q]);
            // fprintf(fp, "\n");
            }
         } // end for k
         
      }
   }
   
// done
//   fclose(fp);
   mQuery.close();

   printf("Called query() %i times\n", nQuery);
   exit(0);
}
