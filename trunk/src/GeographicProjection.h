#ifndef GEOGRAPHIC_PROJECTION_H
#define GEOGRAPHIC_PROJECTION_H

#include <string>

#ifdef ENABLE_ETREE
#include "projects.h"
#endif

class GeographicProjection
{
  public:
     GeographicProjection( double lon_origin, double lat_origin, std::string projection, std::string ellipse, double az );
     void computeGeographicCoord(double x, double y, double & longitude, double & latitude );
     void computeCartesianCoord( double &x, double &y, double longitude, double latitude);
  private:
#ifdef ENABLE_ETREE
     PJ* m_projection;
#endif
     double m_xoffset, m_yoffset, m_az, m_deg2rad;
};

#endif
