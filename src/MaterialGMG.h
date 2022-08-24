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
#ifndef MATERIALGMG_FILE
#define MATERIALGMG_FILE

#include <string>

#include "MaterialData.h"
#include "sw4.h"

class EW;

using namespace std;

class MaterialGMG : public MaterialData
{
  friend class SfileHDF5;
 
  public:
    
    MaterialGMG( EW * a_ew, const std::string file, const std::string directory);
 
   ~MaterialGMG();
 
   void set_material_properties(std::vector<Sarray> & rho, 
 			       std::vector<Sarray> & cs,
 			       std::vector<Sarray> & cp, 
 			       std::vector<Sarray> & xis, 
 			       std::vector<Sarray> & xip);
 
   bool use_attenuation() {return m_use_attenuation;}
   
  protected:
 
    void read_gmg( );
    void fill_in_fluids( );
    void material_check( bool water );

    float mat(int gr, int c, int i, int j, int k) {
#ifdef BZ_DEBUG
      ASSERT(m_npatches > 0);
      ASSERT(gr >= 0 && gr < m_npatches);
      ASSERT( i >= 0 &&  i < m_ni[gr]);
      ASSERT( j >= 0 &&  j < m_nj[gr]);
      ASSERT( k >= 0 &&  k < m_nk[gr]);
      ASSERT( c >= 0 &&  c < m_nc[gr]);
#endif
      return m_Material[gr][i * m_nj[gr] * m_nk[gr] * m_nc[gr] +
                            j * m_nk[gr] * m_nc[gr] +
                            k * m_nc[gr] + 
                            c];
    };

    void mat_assign(int gr, int c, int i, int j, int k, float v) {
#ifdef BZ_DEBUG
      ASSERT(m_npatches > 0);
      ASSERT(gr >= 0 && gr < m_npatches);
      ASSERT( i >= 0 &&  i < m_ni[gr]);
      ASSERT( j >= 0 &&  j < m_nj[gr]);
      ASSERT( k >= 0 &&  k < m_nk[gr]);
      ASSERT( c >= 0 &&  c < m_nc[gr]);
#endif

      m_Material[gr][i * m_nj[gr] * m_nk[gr] * m_nc[gr] +
                     j * m_nk[gr] * m_nc[gr] +
                     k * m_nc[gr] + 
                     c] = v;
      return;
    };
    
    EW* mEW;
 
    std::string m_model_file, m_model_dir;
 
    bool m_use_attenuation;
    int m_npatches;
    double m_Origin_x, m_Origin_y, m_Yaz, m_Zmax, m_Zmin;
    char *m_CRS;
#ifdef USE_HDF5
    hsize_t m_Top_dims[2];
#else
    size_t m_Top_dims[2];
#endif
    float* m_Top_surface;
    vector<double> m_hv, m_hh, m_ztop;
    vector<int> m_ni, m_nj, m_nk, m_nc;
    vector<float*> m_Material;

};
#endif
