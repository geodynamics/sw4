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
#ifndef MATERIALSFILE_FILE
#define MATERIALSFILE_FILE

#include <string>

#include "MaterialData.h"
#include "sw4.h"

class EW;

using namespace std;

class MaterialSfile : public MaterialData
{
  friend class SfileHDF5;
 
  public:
    
    MaterialSfile( EW * a_ew, const std::string file, const std::string directory);
 
   ~MaterialSfile();
 
   void set_material_properties(std::vector<Sarray> & rho, 
 			       std::vector<Sarray> & cs,
 			       std::vector<Sarray> & cp, 
 			       std::vector<Sarray> & xis, 
 			       std::vector<Sarray> & xip);
 
   bool use_attenuation() {return m_use_attenuation;}
   
  protected:
    inline bool inside( float_sw4 x, float_sw4 y, float_sw4 z ) {
      return m_xminloc <= x && x <= m_xmaxloc && m_yminloc <= y && y <= m_ymaxloc 
        && m_zminloc <= z && z <= m_zmaxloc;
    }
 
    void read_sfile( );
    void fill_in_fluids( );
    void material_check( bool water );
    
    EW* mEW;
 
    std::string m_model_file, m_model_dir;
 
    bool m_use_attenuation;
    int m_npatches;

    vector<int> m_ifirst, m_ilast, m_jfirst, m_jlast, m_kfirst, m_klast, m_ni, m_nj, m_nk;
    vector<float_sw4> m_hh;
    double m_x0, m_y0, m_lon0, m_lat0, m_azim;
 
    // xminloc, xmaxloc, etc. is the bounding box for the set of data patches in this processor.
    float_sw4 m_xminloc, m_xmaxloc, m_yminloc, m_ymaxloc, m_zminloc, m_zmaxloc;
    bool m_outside;

    vector<Sarray> mMaterial_rho;
    vector<Sarray> mMaterial_cp;
    vector<Sarray> mMaterial_cs;
    vector<Sarray> mMaterial_qp;
    vector<Sarray> mMaterial_qs;

    vector<Sarray> mInterface;
    vector<bool> m_isempty;
};
#endif
