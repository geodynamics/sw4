#include <mpi.h>

#include <cstring>
#include "GeographicProjection.h"
#include "Require.h"

using namespace std;

//-----------------------------------------------------------------------
GeographicProjection::GeographicProjection( double lon_origin, double lat_origin,
					    string projection, double az )
{
#ifdef ENABLE_PROJ4

   m_projection = pj_init_plus(projection.c_str());
   CHECK_INPUT( m_projection != 0, "ERRROR: Init of cartographic projection failed with message: " << pj_strerrno(pj_errno) );

   m_latlong = pj_init_plus("+proj=latlong");
   CHECK_INPUT( m_latlong != 0, "ERRROR: Init of latlong projection failed with message: " << pj_strerrno(pj_errno) );

   m_deg2rad = M_PI/180;

   double x0 = lon_origin*DEG_TO_RAD;
   double y0 = lat_origin*DEG_TO_RAD;
   int status = pj_transform(m_latlong, m_projection, 1, 1, &x0, &y0, NULL );
   printf("Origin mapped from (lon,lat)=(%e, %e) to (x0,y0)=(%e, %e)\n", lon_origin, lat_origin, x0, y0);
   m_xoffset = x0;
   m_yoffset = y0;

   m_az = az*m_deg2rad;
#endif
}

//-----------------------------------------------------------------------
void GeographicProjection::computeGeographicCoord(double x, double y,
						  double & longitude, double & latitude)
{
#ifdef ENABLE_PROJ4
   double xmap, ymap;
   int status;

   xmap = x*sin(m_az) + y*cos(m_az) + m_xoffset;
   ymap = x*cos(m_az) - y*sin(m_az) + m_yoffset;
   status = pj_transform(m_projection, m_latlong, 1, 1, &xmap, &ymap, NULL );
   longitude = xmap*RAD_TO_DEG;
   latitude  = ymap*RAD_TO_DEG;

#endif
}

//-----------------------------------------------------------------------
void GeographicProjection::computeCartesianCoord(double &x, double &y,
						 double lon, double lat )
{
#ifdef ENABLE_PROJ4
  double xlon, ylat;
  int status;
  
  xlon = lon*DEG_TO_RAD;
  ylat = lat*DEG_TO_RAD;
  status = pj_transform(m_latlong, m_projection, 1, 1, &xlon, &ylat, NULL );
  xlon -= m_xoffset;
  ylat -= m_yoffset;
  x =  xlon*sin(m_az) + ylat*cos(m_az);
  y =  xlon*cos(m_az) - ylat*sin(m_az);
#endif
}

