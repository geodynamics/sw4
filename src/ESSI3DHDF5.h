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
#ifndef SW4_ESSI3DHDF5_H
#define SW4_ESSI3DHDF5_H

#include <string>
#include <sstream>
#include <vector>

#ifdef USE_HDF5
#include "hdf5.h"
#endif

class ESSI3DHDF5
{
public:
  static ESSI3DHDF5* nil;

  ESSI3DHDF5(const std::string& filename, int (&global)[3],
    int (&window)[6], bool ihavearray, int precision);
  ~ESSI3DHDF5();

  void create_file(bool is_open);
  void close_file();
  void write_header(double h, double (&lonlat_origin)[2], double az,
    double (&origin)[3], int cycle, double t, double dt);
  void write_topo(void* window_array);

  void write_vel(void* window_array, int comp, int cycle, int nstep);

  void init_write_vel(bool m_isRestart, int ntimestep, int ZFPmode, double ZFPpar, int dumpInterval);

  const std::string& filename() {return m_filename;};
  void set_ihavearray(bool ihavearray) {m_ihavearray=ihavearray;};
  void finalize_hdf5();

protected:

private:
  ESSI3DHDF5(); // make it impossible to call default constructor
  ESSI3DHDF5(const ESSI3DHDF5 &in); // hide copy constructor

  std::string m_filename;
  bool m_ihavearray;
  int m_end_cycle;
  int m_window[6];
  int m_global[3];
  int m_precision;

#ifdef USE_HDF5
  hsize_t m_window_dims[4]; // for just this proc, this cycle
  hsize_t m_global_dims[4]; // unlimited
  hsize_t m_cycle_dims[4]; // for all cycles until now
  hsize_t m_slice_dims[4]; // for just a column

  hid_t m_file_id;
  hid_t m_es_id;
#endif // def USE_HDF5
};

#endif
