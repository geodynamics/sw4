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
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
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
#include "sw4.h"

class EW;

class MaterialIfile : public MaterialData {
 public:
  MaterialIfile(EW* a_ew, std::string fileName, bool CartesianFormat);

  void set_material_properties(std::vector<Sarray>& rho,
                               std::vector<Sarray>& cs, std::vector<Sarray>& cp,
                               std::vector<Sarray>& xis,
                               std::vector<Sarray>& xip);

 protected:
  void extractSurfaceFromGridFile(std::string a_surfaceFileName);
  void extractSurfaceFromCartesianFile(std::string a_surfaceFileName);
  sw4_type getMaterialID(double lat, double lon, float_sw4 depth);
  sw4_type getCartesianMaterialID(float_sw4 xP, float_sw4 yP, float_sw4 depth);

  inline bool inside_material_surfaces(double lat, double lon) {
    return (lat <= m_materialLatMax && lat >= m_materialLatMin &&
            lon <= m_materialLonMax && lon >= m_materialLonMin);
  }
  inline bool inside_cartesian_material_surfaces(float_sw4 xP, float_sw4 yP) {
    return (yP <= m_mat_Ymax && yP >= m_mat_Ymin && xP <= m_mat_Xmax &&
            xP >= m_mat_Xmin);
  }
  inline float_sw4 lookup_Rho(MaterialProperty* prop, float_sw4 depth) {
    return prop->m_rho0 + prop->m_rho1 * depth + prop->m_rho2 * depth * depth +
           prop->m_rho1o2 * sqrt(fabs(depth));
  }
  inline float_sw4 lookup_Vs(MaterialProperty* prop, float_sw4 depth) {
    return prop->m_vs0 + prop->m_vs1 * depth + prop->m_vs2 * depth * depth +
           prop->m_vs1o2 * sqrt(fabs(depth));
  }
  inline float_sw4 lookup_Vp(MaterialProperty* prop, float_sw4 depth) {
    return prop->m_vp0 + prop->m_vp1 * depth + prop->m_vp2 * depth * depth +
           prop->m_vp1o2 * sqrt(fabs(depth));
  }
  inline float_sw4 lookup_Qs(MaterialProperty* prop, float_sw4 depth) {
    return prop->m_qs;
  }
  inline float_sw4 lookup_Qp(MaterialProperty* prop, float_sw4 depth) {
    return prop->m_qp;
  }
  // General variables
  bool m_mat_Cartesian;
  sw4_type m_number_material_surfaces;
  Sarray m_materialDepth;
  EW* mEw;

  // Variables for geographic coords
  double m_Nlon, m_Nlat;
  double m_materialLonMax, m_materialLonMin, m_materialLatMax, m_materialLatMin;
  double *m_materialLon, *m_materialLat;

  // Variables for Cartesian coords
  sw4_type m_mat_Nx, m_mat_Ny;
  float_sw4 m_mat_Xmax, m_mat_Xmin, m_mat_Ymax, m_mat_Ymin;
  float_sw4 *m_mat_Xvec, *m_mat_Yvec;
};

#endif
