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
#ifndef GEOGRAPHIC_PROJECTION_H
#define GEOGRAPHIC_PROJECTION_H

#include <string>

#ifdef ENABLE_PROJ4
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H 
#include "proj_api.h"
#endif

#ifdef ENABLE_PROJ_6
#include "proj.h"
#endif

class GeographicProjection
{
  public:
     GeographicProjection( double lon_origin, double lat_origin, std::string projection, double az );
     ~GeographicProjection();
     void computeGeographicCoord(double x, double y, double & longitude, double & latitude );
     void computeCartesianCoord( double & x, double & y, double longitude, double latitude );
     void computeCartesianCoordGMG( double & x, double & y, double longitude, double latitude, char* crs_to );
  private:
#ifdef ENABLE_PROJ4
     projPJ m_projection, m_latlong;
#endif
#ifdef ENABLE_PROJ_6
   PJ *m_P;
   PJ *m_Pgmg;
#endif
   double m_xoffset, m_yoffset, m_az, m_deg2rad;
};

#endif
