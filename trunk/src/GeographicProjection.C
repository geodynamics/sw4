#include <cstring>
#include "GeographicProjection.h"
#include "Require.h"

using namespace std;

//-----------------------------------------------------------------------
GeographicProjection::GeographicProjection( double lon_origin, double lat_origin,
					    string projection, string ellipse, double az )
{
#ifdef ENABLE_ETREE
   int nprojpars = 6;
   char** projpars = new char*[nprojpars];
   for( int n=0 ; n < nprojpars ; n++ )
      projpars[n] = new char[35];
   strncpy(projpars[0],projection.c_str(),35);
   strncpy(projpars[1],ellipse.c_str(),35);
   snprintf(projpars[2],35,"lon_0=%25.17g",lon_origin);
   snprintf(projpars[3],35,"lat_0=%25.17g",lat_origin);
   strncpy(projpars[4],"units=m",35);
   strncpy(projpars[5],"no_defs",35);

   m_projection = pj_init( nprojpars, projpars );
   CHECK_INPUT( m_projection != 0, "ERRROR: first call to pj_init failed with message: " << pj_strerrno(pj_errno) );
     
   m_deg2rad = M_PI/180;
   projUV xy, lonlat;
   lonlat.u = lon_origin*m_deg2rad;
   lonlat.v = lat_origin*m_deg2rad;
   xy = pj_fwd( lonlat, m_projection );
   CHECK_INPUT( xy.u != HUGE_VAL, "ERROR: first call to pj_fwd failed with message " << pj_strerrno(pj_errno) );
   m_xoffset = xy.u;
   m_yoffset = xy.v;

   for( int n=0 ; n < nprojpars ; n++ )
      delete[] projpars[n];
   delete[] projpars;
   m_az = az*m_deg2rad;
#endif
}

//-----------------------------------------------------------------------
void GeographicProjection::computeGeographicCoord(double x, double y,
						  double & longitude, double & latitude)
{
#ifdef ENABLE_ETREE
 projUV lonlat, xy;
 xy.u = x*sin(m_az) + y*cos(m_az) + m_xoffset;
 xy.v = x*cos(m_az) - y*sin(m_az) + m_yoffset;
 lonlat = pj_inv( xy, m_projection );
  if( lonlat.u == HUGE_VAL )
     cout << "ERROR: computeGeographicCoord, pj_inv failed with message " << pj_strerrno(pj_errno) << endl;
 longitude = lonlat.u/m_deg2rad;
 latitude  = lonlat.v/m_deg2rad;
#endif
}

//-----------------------------------------------------------------------
void GeographicProjection::computeCartesianCoord(double &x, double &y,
						 double lon, double lat )
{
#ifdef ENABLE_ETREE
  projUV lonlat, xy;
  lonlat.u = lon*m_deg2rad;
  lonlat.v = lat*m_deg2rad;
  xy = pj_fwd( lonlat, m_projection );
  if( xy.u == HUGE_VAL )
     cout << "ERROR: computeCartesianCoord pj_fwd failed with message " << pj_strerrno(pj_errno) << endl;
  xy.u -= m_xoffset;
  xy.v -= m_yoffset;
  x =  xy.u*sin(m_az) + cos(m_az)*xy.v;
  y =  xy.u*cos(m_az) - sin(m_az)*xy.v;
#endif
}

