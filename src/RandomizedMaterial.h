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
#ifndef RANDOMIZEDMATERIAL_FILE
#define RANDOMIZEDMATERIAL_FILE

#include <complex>
#include <string>
#include <vector>

#include "sw4.h"
//#include "AllDims.h"

class AllDims;
class EW;

using namespace std;

class RandomizedMaterial {
 public:
  RandomizedMaterial(EW* a_ew, float_sw4 zmin, float_sw4 zmax,
                     float_sw4 corrlen, float_sw4 corrlenz, float_sw4 hurst,
                     float_sw4 sigma, float_sw4 rhoamplitude = 0.0,
                     bool randomrho = false, unsigned int seed = 0);
  ~RandomizedMaterial();
  void perturb_velocities(int g, Sarray& cs, Sarray& cp, double h, double zmin,
                          double zmax);

  void perturb_velocities(std::vector<Sarray>& cs, std::vector<Sarray>& cp);
  void assign_perturbation(int g, Sarray& pert, Sarray& cs, double h,
                           double zmin, double zmax, bool rho);

  void set_vsmax(float_sw4 vsmax);
  double get_vsmax();
  void set_vsmin(float_sw4 vsmin);
  double get_vsmin();
  bool randomize_rho() { return m_random_rho; }
  float_sw4 rhoamplitude() { return m_rhoamplitude; }

 private:
  void gen_random_mtrl_fft3d_fftw(int n1g, int n2g, int n3g, float_sw4 Lx,
                                  float_sw4 Ly, float_sw4 Lz, float_sw4 hurst);

  void get_fourier_modes(std::complex<float_sw4>* uhat, int n1, int ib1,
                         int n1g, int n2, int n3, float_sw4 l1, float_sw4 l2,
                         float_sw4 l3, float_sw4 hurst, unsigned int seed);
  void rescale_perturbation();
  template <class T>
  void redistribute_array(AllDims& src, AllDims& dest, T* src_array,
                          T* dest_array);

  void repad_sarray(Sarray& sar, int old_padding, int new_padding);
  void comm_sarray(Sarray& sar, int neigh[4], int padding);

  inline bool inside(float_sw4 x, float_sw4 y, float_sw4 z) {
    return m_xminloc <= x && x <= m_xmaxloc && m_yminloc <= y &&
           y <= m_ymaxloc && m_zminloc <= z && z <= m_zmaxloc;
  }

  EW* mEW;

  // Index range of patch in file coordinate system
  int m_ifirst, m_ilast, m_jfirst, m_jlast, m_kfirst, m_klast, m_ni, m_nj, m_nk;
  int m_nig, m_njg, m_nkg;

  // file coordinate system is x=(i-1)*m_hx[gr] + m_xmin[gr], in SW4
  // coordinates.
  float_sw4 m_hh, m_hv, m_zmin, m_zmax, m_corrlen, m_corrlenz, m_hurst, m_sigma,
      m_vsmax;
  float_sw4 m_vsmin;
  //   float_sw4 m_x0, m_y0;
  unsigned int m_seed;
  int m_nproc2d[2];

  // xminloc, xmaxloc, etc. is the bounding box for the set of data patches in
  // this processor.
  float_sw4 m_xminloc, m_xmaxloc, m_yminloc, m_ymaxloc, m_zminloc, m_zmaxloc;
  float_sw4 m_rhoamplitude;
  bool m_outside, m_random_rho;

  // 3-dimensional Sarrays
  Sarray mRndMaterial;
  bool m_isempty;
};
#endif
