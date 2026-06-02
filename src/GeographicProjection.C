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
#include <mpi.h>

#include <cstring>
#include <cstdio>
#include "GeographicProjection.h"
#include "Require.h"

using namespace std;

static const char* sw4_geographic_crs()
{
   return "+proj=latlong +datum=NAD83";
}

#ifdef ENABLE_PROJ
static PJ* create_lonlat_xy_transform(const char* crs_from,
                                      const char* crs_to)
{
  PJ* transform = proj_create_crs_to_crs(PJ_DEFAULT_CTX, crs_from, crs_to, NULL);
  if (transform == NULL)
    return NULL;
  PJ* normalized = proj_normalize_for_visualization(PJ_DEFAULT_CTX, transform);
  if (normalized != NULL) {
    proj_destroy(transform);
    transform = normalized;
  }
  return transform;
}
#endif

//-----------------------------------------------------------------------
GeographicProjection::GeographicProjection( double lon_origin, double lat_origin,
					    string projection, double az )
{
   m_geographic_crs = sw4_geographic_crs();
   m_gmg_crs_from_cache.clear();
   m_gmg_crs_to_cache.clear();
/* #ifdef ENABLE_PROJ4 */
/*    m_projection = pj_init_plus(projection.c_str()); */
/*    CHECK_INPUT( m_projection != 0, "ERRROR: Init of cartographic projection failed with message: " << pj_strerrno(pj_errno) ); */

/*    m_latlong = pj_init_plus("+proj=latlong +datum=NAD83"); */
/*    CHECK_INPUT( m_latlong != 0, "ERRROR: Init of latlong projection failed with message: " << pj_strerrno(pj_errno) ); */

/*    m_deg2rad = M_PI/180; */

/*    double x0 = lon_origin*DEG_TO_RAD; */
/*    double y0 = lat_origin*DEG_TO_RAD; */
/*    int status = pj_transform(m_latlong, m_projection, 1, 1, &x0, &y0, NULL ); */
/* // tmp */
/* //   printf("Origin mapped from (lon,lat)=(%e, %e) to (x0,y0)=(%e, %e)\n", lon_origin, lat_origin, x0, y0); */
/*    m_xoffset = x0; */
/*    m_yoffset = y0; */

/*    m_az = az*m_deg2rad; */
/* #endif */

#ifdef ENABLE_PROJ
   PJ_COORD c, c_out;
   const char *crs_from = m_geographic_crs.c_str();
   const char *crs_to = projection.c_str();

   m_P = create_lonlat_xy_transform(crs_from, crs_to);
   /* printf("projection: [%s]\n", crs_to); */
   ASSERT(m_P);

   c = proj_coord(lon_origin, lat_origin, 0.0, 0.0);
   c_out = proj_trans(m_P, PJ_FWD, c);

   m_xoffset = c_out.xyzt.x;
   m_yoffset = c_out.xyzt.y;

   m_deg2rad = M_PI/180;
   m_az = az*m_deg2rad;

   m_Pgmg = NULL;
#endif
   /* printf("GeographicProjection origin: %f %f -> %f %f\n", lon_origin, lat_origin, m_xoffset, m_yoffset); */
}

//-----------------------------------------------------------------------
GeographicProjection::~GeographicProjection(void)
{
#ifdef ENABLE_PROJ
   if (m_P)
       proj_destroy(m_P);
   if (m_Pgmg)
       proj_destroy(m_Pgmg);
   m_P = NULL;
   m_Pgmg = NULL;
#endif
}

//-----------------------------------------------------------------------
void GeographicProjection::computeGeographicCoord(double x, double y,
						  double & longitude, double & latitude)
{
/* #ifdef ENABLE_PROJ4 */
/*    double xmap, ymap; */
/*    int status; */

/*    xmap = x*sin(m_az) + y*cos(m_az) + m_xoffset; */
/*    ymap = x*cos(m_az) - y*sin(m_az) + m_yoffset; */
/*    status = pj_transform(m_projection, m_latlong, 1, 1, &xmap, &ymap, NULL ); */
/*    longitude = xmap*RAD_TO_DEG; */
/*    latitude  = ymap*RAD_TO_DEG; */

/* #endif */

#ifdef ENABLE_PROJ
   PJ_COORD c, c_out;
   ASSERT(m_P);

   double xmap, ymap;
   xmap = x*sin(m_az) + y*cos(m_az) + m_xoffset;
   ymap = x*cos(m_az) - y*sin(m_az) + m_yoffset;

   /* No need to convert to radian */
   c = proj_coord(xmap, ymap, 0.0, 0.0);
   c_out = proj_trans(m_P, PJ_INV, c);

   longitude = c_out.xyzt.x;
   latitude  = c_out.xyzt.y;

#endif

   /* printf("computeGeographicCoord: %f %f -> %f\t%f\n", x, y, longitude, latitude); */
}

//-----------------------------------------------------------------------
void GeographicProjection::computeCartesianCoord(double &x, double &y,
						 double lon, double lat )
{
/* #ifdef ENABLE_PROJ4 */
/*   double xlon, ylat; */
/*   int status; */
  
/*   xlon = lon*DEG_TO_RAD; */
/*   ylat = lat*DEG_TO_RAD; */
/*   status = pj_transform(m_latlong, m_projection, 1, 1, &xlon, &ylat, NULL ); */
/*   xlon -= m_xoffset; */
/*   ylat -= m_yoffset; */
/*   x =  xlon*sin(m_az) + ylat*cos(m_az); */
/*   y =  xlon*cos(m_az) - ylat*sin(m_az); */
/* #endif */

#ifdef ENABLE_PROJ
  PJ_COORD c, c_out;
  double xlon, ylat;

  ASSERT(m_P);

  c = proj_coord(lon, lat, 0.0, 0.0);
  c_out = proj_trans(m_P, PJ_FWD, c);

  xlon = c_out.xyzt.x;
  ylat = c_out.xyzt.y;

  xlon -= m_xoffset;
  ylat -= m_yoffset;
  x =  xlon*sin(m_az) + ylat*cos(m_az);
  y =  xlon*cos(m_az) - ylat*sin(m_az);
#endif

   /* printf("computeCartesianCoord: %f %f -> %f\t%f\n", lon, lat, x, y); */
}

//-----------------------------------------------------------------------
void GeographicProjection::computeCartesianCoordGMG(double &x, double &y,
						 double lon, double lat, char* crs_to )
{
  x = 0.0, y = 0.0;
#ifdef ENABLE_PROJ
  PJ_COORD c, c_out;
  const char *crs_from = m_geographic_crs.c_str();
  const string gmg_crs_to = crs_to ? string(crs_to) : string();

  CHECK_INPUT(!gmg_crs_to.empty(), "ERROR: empty GMG target coordinate reference system");

  /* printf("computeCartesianCoordGMG: crs_to %s\n", crs_to); */

  if (m_Pgmg == NULL || m_gmg_crs_from_cache != m_geographic_crs ||
      m_gmg_crs_to_cache != gmg_crs_to) {
    if (m_Pgmg) proj_destroy(m_Pgmg);
    m_Pgmg = create_lonlat_xy_transform(crs_from, gmg_crs_to.c_str());
    m_gmg_crs_from_cache = m_geographic_crs;
    m_gmg_crs_to_cache = gmg_crs_to;
  }

  ASSERT(m_Pgmg);

  c = proj_coord(lon, lat, 0.0, 0.0);
  c_out = proj_trans(m_Pgmg, PJ_FWD, c);

  x = c_out.xyzt.x;
  y = c_out.xyzt.y;

#else
  printf("GMG format only works with proj 6+, abort!\n");
  ASSERT(0);
#endif

  /* printf("computeCartesianCoordGMG: %f %f -> %f\t%f\n", lon, lat, x, y); */
}
