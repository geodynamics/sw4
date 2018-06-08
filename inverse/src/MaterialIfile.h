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
#ifndef MATERIALIFILE_H
#define MATERIALIFILE_H

#include <string>

#include "MaterialData.h"
#include "MaterialProperty.h"

class EW;

class MaterialIfile : public MaterialData
{
public:
   MaterialIfile( EW * a_ew, std::string fileName, bool CartesianFormat );

   void set_material_properties( std::vector<Sarray> &rho, std::vector<Sarray> &cs,
				 std::vector<Sarray> &cp,
				 std::vector<Sarray>& xis, std::vector<Sarray>& xip);

protected:
   void extractSurfaceFromGridFile(std::string a_surfaceFileName);
   void extractSurfaceFromCartesianFile(std::string a_surfaceFileName);
   int getMaterialID(double lat, double lon, double depth );
   int getCartesianMaterialID(double xP, double yP, double depth );

   inline bool inside_material_surfaces( double lat, double lon )
   {
      return (lat <= m_materialLatMax && lat >= m_materialLatMin && 
	      lon <= m_materialLonMax && lon >= m_materialLonMin);
   }
   inline bool inside_cartesian_material_surfaces( double xP, double yP )
   {
      return (yP <= m_mat_Ymax && yP >= m_mat_Ymin && 
	      xP <= m_mat_Xmax && xP >= m_mat_Xmin);
   }
   inline double lookup_Rho( MaterialProperty* prop, double depth )
   {
      return prop->m_rho0 + prop->m_rho1*depth + prop->m_rho2*depth*depth +
	 prop->m_rho1o2*sqrt(fabs(depth));
   }
   inline double lookup_Vs( MaterialProperty* prop, double depth )
   {
      return prop->m_vs0 + prop->m_vs1*depth + prop->m_vs2*depth*depth +
	 prop->m_vs1o2*sqrt(fabs(depth));
   }
   inline double lookup_Vp( MaterialProperty* prop, double depth )
   {
      return prop->m_vp0 + prop->m_vp1*depth + prop->m_vp2*depth*depth +
	 prop->m_vp1o2*sqrt(fabs(depth));
   }
   inline double lookup_Qs( MaterialProperty* prop, double depth )
   {
      return prop->m_qs;
   }
   inline double lookup_Qp( MaterialProperty* prop, double depth )
   {
      return prop->m_qp;
   }
   // General variables
   bool m_mat_Cartesian;
   int m_number_material_surfaces;
   Sarray m_materialDepth;
   EW *mEw;

   // Variables for geographic coords
   double m_Nlon, m_Nlat;
   double m_materialLonMax, m_materialLonMin, m_materialLatMax, m_materialLatMin;
   double *m_materialLon, *m_materialLat;

   // Variables for Cartesian coords
   int m_mat_Nx, m_mat_Ny;
   double m_mat_Xmax, m_mat_Xmin, m_mat_Ymax, m_mat_Ymin;
   double *m_mat_Xvec, *m_mat_Yvec;

};

#endif
