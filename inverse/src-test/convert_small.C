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
#include <cstring>
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
   
   double dm_x[4], dm_y[4], len[4], dx, dy, alpha, beta, origin_x, origin_y, origin_lon, origin_lat;

   projPJ pj_merc, pj_latlong;
   double x, y, z, xlon, ylat, x0, y0, xr, yr, xc, yc, xabs, yabs, xmap, ymap;
   int status, merc_len;
   char merc_def[256];
   
/* should parse the file name from the arguments... */
   std::string filename="/Users/petersson1/USGSBayAreaVM-08.3.0.efile"; // Hard-coded filename!

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
   if (!(pj_latlong = pj_init_plus("+proj=latlong +datum=NAD83")) )
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
// calculate azimuths (not used anymore)
   // dx = dm_x[0] - dm_x[1];
   // dy = dm_y[0] - dm_y[1];
   // beta = atan2(dy,dx);

   alpha = DEG_TO_RAD*alpha_deg;
   
   printf("Azimuth angle (x-axis) Alpha = %e\n", RAD_TO_DEG*alpha);
   
// tmp
//   exit(0);
   
//
   origin_lon =-122.25;
   origin_lat = 37.93;
   
   // origin_lon =-122.56;
   // origin_lat = 39.012585;
   
   xlon = origin_lon*DEG_TO_RAD;
   ylat = origin_lat*DEG_TO_RAD;
   status = pj_transform(pj_latlong, pj_merc, 1, 1, &xlon, &ylat, NULL );
// these are the (x,y) coordinates of our origin
   origin_x = xlon-x0;
   origin_y = ylat-y0;

// calculate (lon,lat) coordinates along cell centers in detailed model

   const char* mQueryKeys[]={"Density", "Vp", "Vs", "elevation", "Qp", "Qs","FaultBlock"};
   int mPayloadSize=7;
   double *mPayload, elev;
   cencalvm::query::VMQuery mQuery;
   cencalvm::storage::Geometry* mQueryGeom;
   
   mPayload = new double[7];
   
// open up the data base file
   mQuery.filename(filename.c_str());
   mQuery.queryType(cencalvm::query::VMQuery::MAXRES);
// Set values to be returned in queries
   mQuery.queryVals(mQueryKeys, mPayloadSize);
   mQuery.open();

// how many blocks are stored?
   int nblocks=4;
   int attenuation=1; // saving attenuation (Qp & Qs)?

// float=4, double=8
   int prec=4; 
   printf("Saving data with %i bytes of precision\n", prec);
   
// number of points for different grid sizes
   int nimax[]={121, 121, 61, 31}; 
   int njmax[]={121, 121, 61, 31};
// how many points in the vertical direction?
//   int nkmax[]={1, 74, 57, 33, 194};
   int nkmax[]={1, 39, 57, 33, 194};

// how many components in each block?
   int nc[]={1, 5, 5, 5, 5};
   
// block header info
// cell size in horizontal directions
   double clh[] = {100.0, 100.0, 200.0, 400.0, 800.0};
// cell size in vertical direction
   double clv[] = {25.0, 25.0, 50.0, 100.0, 200.0};
// topo does not use z0, but needed by format
//   double z0[] = {0.0, -1437.5, 387.5, 3187.5, 6387.5};
   double z0[] = {0.0, -562.5, 387.5, 3187.5, 6387.5};
   
// level = block
   int lev;
   
// tmp storage for mat prop
   double mat[5];
   float felev, fmat[5];

// origin offset
   xc = 0.0;
   yc = 0.0;

// magic number
   int magic=1;

   FILE *fp=fopen("rfile.dat","wb"); // Another hard-coded file name
   FILE *eh=fopen("warnings.txt","w"); // Another hard-coded file name

// write header
   fwrite(&magic, sizeof(int), 1, fp);
   fwrite(&prec, sizeof(int), 1, fp);
   fwrite(&attenuation, sizeof(int), 1, fp);
// azimuth
   fwrite(&alpha_deg, sizeof(double), 1, fp);
// origin longitude (NE corner)
//   fwrite(&(dm_lon[3]), sizeof(double), 1, fp);
   fwrite(&(origin_lon), sizeof(double), 1, fp);
// origin latitude (NE corner)
//   fwrite(&(dm_lat[3]), sizeof(double), 1, fp);
   fwrite(&(origin_lat), sizeof(double), 1, fp);
// number of characters in projection string
   fwrite(&merc_len, sizeof(int), 1, fp);
// projection string
   fwrite(&merc_def, sizeof(char), merc_len, fp);

// number of blocks
   fwrite(&nblocks, sizeof(int), 1, fp);

// block # 1: topography on a 2-D grid
   for (lev=0; lev<nblocks; lev++)
   {
      fwrite(&(clh[lev]), sizeof(double), 1, fp); // horizontal grid size
      fwrite(&(clv[lev]), sizeof(double), 1, fp); // vertical grid size
      fwrite(&z0[lev], sizeof(double), 1, fp); // starting z-level

// block sizes (number of grid points)
      fwrite(&(nc[lev]), sizeof(int), 1, fp); // one component per grid point
      fwrite(&(nimax[lev]), sizeof(int), 1, fp);
      fwrite(&(njmax[lev]), sizeof(int), 1, fp);
      fwrite(&(nkmax[lev]), sizeof(int), 1, fp);
   }
   
   int nQuery = 0;
   
// save topography
   lev = 0;
   printf("Querying etree for topography at %i by %i points\n", nimax[lev], njmax[lev]);

   double mintopo=1e6, maxtopo=-1e6;
   double minlat=1e6, maxlat=-1e6;
   double minlon=1e6, maxlon=-1e6;
   
// store data in "C" order
   for (int i=0; i<nimax[lev]; i++)
   {
      printf("..%i", i);
      for (int j=0; j<njmax[lev]; j++)
      {
         xr = xc + i*clh[lev];
         yr = yc + j*clh[lev];
      
// coordinate system with origin at NE corner (dm_x[3], dm_y[3])
         // xabs = xmap = dm_x[3] + xr*sin(alpha) + yr*cos(alpha);
         // yabs = ymap = dm_y[3] + xr*cos(alpha) - yr*sin(alpha);

// coordinate system with origin at (origin_x, origin_y)
         xabs = xmap = origin_x + xr*sin(alpha) + yr*cos(alpha);
         yabs = ymap = origin_y + xr*cos(alpha) - yr*sin(alpha);

// inverse projection to get (lon,lat)   
         status = pj_transform(pj_merc, pj_latlong, 1, 1, &xmap, &ymap, NULL );

         xlon = xmap*RAD_TO_DEG;
         ylat = ymap*RAD_TO_DEG;
   
//      printf("cell center (x,y)=(%e, %e), (lon, lat)=(%e, %e)\n", xabs, yabs, xlon, ylat);

// min/max
      if (ylat < minlat) minlat=ylat;
      if (ylat > maxlat) maxlat=ylat;

      if (xlon < minlon) minlon=xlon;
      if (xlon > maxlon) maxlon=xlon;

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
            mPayload[3] = -999.0;
         }

// min/max
         if (mPayload[3] < mintopo) mintopo=mPayload[3];
         if (mPayload[3] > maxtopo) maxtopo=mPayload[3];

	 if (prec==8)
	 {
	   fwrite(&(mPayload[3]), sizeof(double), 1, fp);
	 }
	 else
	 {
	   felev = (float) mPayload[3];
	   fwrite(&felev, sizeof(float), 1, fp);
	 }
      }
   }
   printf("\n");

   printf("min/max topo= %e, %e\n", mintopo, maxtopo);
   printf("min/max lat= %e, %e\n", minlat, maxlat);
   printf("min/max lon= %e, %e\n", minlon, maxlon);
   
   double minmat[3], maxmat[3];

// save material properties for block 1...,nblocks
   for (lev = 1; lev<nblocks; lev++)
   {
      printf("Querying the material properties for block=%i\n", lev);

      for (int q=0; q<3; q++) minmat[q] = 1e5;
      for (int q=0; q<3; q++) maxmat[q] = 0;
      
// store data in "C" order
      for (int i=0; i<nimax[lev]; i++)
      {
//         printf("Querying the material properties for block=%i, i=%i\n", lev, i);
   
         for (int j=0; j<njmax[lev]; j++)
         {
	   for (int k=0; k<nkmax[lev]; k++)
	   {
	     z = z0[lev] + k*clv[lev];
	     elev = -z;

	     xr = xc + i*clh[lev];
	     yr = yc + j*clh[lev];
      
// coordinate system centered at NE corner (#3)
	     // xabs = xmap = dm_x[3] + xr*sin(alpha) + yr*cos(alpha);
	     // yabs = ymap = dm_y[3] + xr*cos(alpha) - yr*sin(alpha);
// coordinate system with origin at (origin_x, origin_y)
         xabs = xmap = origin_x + xr*sin(alpha) + yr*cos(alpha);
         yabs = ymap = origin_y + xr*cos(alpha) - yr*sin(alpha);

// inverse projection to get (lon,lat)   
	     status = pj_transform(pj_merc, pj_latlong, 1, 1, &xmap, &ymap, NULL );

	     xlon = xmap*RAD_TO_DEG;
	     ylat = ymap*RAD_TO_DEG;
   
//      printf("cell center (x,y)=(%e, %e), (lon, lat)=(%e, %e)\n", xabs, yabs, xlon, ylat);

// query the material model
	     mQuery.query(&mPayload, mPayloadSize, xlon, ylat, elev);
	     nQuery++;
         
// Make sure the query didn't generated a warning or error
	     if (mQuery.errorHandler()->status() == cencalvm::storage::ErrorHandler::ERROR) 
	     {
	       printf("Query error for indices (ix, jy, kz)=(%i, %i, %i), elev=%e\n", i, j, k, elev);
// phony values
	       for (int q=0; q<6; q++)
		 mPayload[q] = -1.0;
	     }
// SHOULD ALSO check for warnings
	     if (mQuery.errorHandler()->status() != cencalvm::storage::ErrorHandler::OK) 
	     {
	       fprintf(eh,"lev=%i, coord=(%e %e) elev=%e, status=%i\n", lev, xr, yr, 
		       elev, mQuery.errorHandler()->status() );
	     }

// save data
	     mat[0]=mPayload[0]; // dens
	     mat[1]=mPayload[1]; // Vp
	     mat[2]=mPayload[2]; // Vs

	     for (int q=0; q<3; q++)
	     {
	       if (mat[q] < minmat[q]) minmat[q] = mat[q];
	       if (mat[q] > maxmat[q]) maxmat[q] = mat[q];
	       fmat[q] = (float) mat[q];
	     }
               
               
	     if (attenuation)
	     {
	       mat[3]=mPayload[4]; // Qp
	       mat[4]=mPayload[5]; // Qs
	       fmat[3] = (float) mat[3];
	       fmat[4] = (float) mat[4];
	       if (prec == 8)
	       {
		 fwrite(mat, sizeof(double), 5, fp);
	       }
	       else
	       {
		 fwrite(fmat, sizeof(float), 5, fp);
	       }
	     }
	     else
	     {
	       if (prec == 8)
	       {
		 fwrite(mat, sizeof(double), 3, fp);
	       }
	       else
	       {
		 fwrite(fmat, sizeof(float), 3, fp);
	       }               
	     }
            
// reset status for next query
	     mQuery.errorHandler()->resetStatus();
	   } // end for i...
         } // end for j...
      } // end for k...
      printf("Block %i: Min dens, Vp, Vs= %e %e %e\n", lev, minmat[0], minmat[1], minmat[2]);
      printf("Block %i: Max dens, Vp, Vs= %e %e %e\n", lev, maxmat[0], maxmat[1], maxmat[2]);
      z = z0[lev] + (nkmax[lev]-1)*clv[lev];
      printf("Block %i: Max z= %e\n", lev, z);
   } // end for block
   
   
// done
   fclose(fp);
   fclose(eh);
   mQuery.close();

   printf("Called query() %i times\n", nQuery);
   exit(0);
}
