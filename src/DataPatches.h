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
#ifndef DATAPATCHES_H
#define DATAPATCHES_H

#include <string>
#include <vector>

#include "Sarray.h"

class DataPatches {
  // Number of subcubes in domain to save
  sw4_type m_npatches;

  // dimensions of cube p, stored as:
  //  imin  = m_dims[6*p],   imax = m_dims[6*p+1]
  //  jmin  = m_dims[6*p+2], jmax = m_dims[6*p+3]
  //  kmin  = m_dims[6*p+4], kmax = m_dims[6*p+5]
  std::vector<size_t> m_dims;

  // pointer to data array. Start index of subcube p is m_dataptr[p]
  size_t* m_dataptr;

  // Number of components per grid point to store, typically 3 for the 3D
  // elastic wave eq
  sw4_type m_ncomp;

  // First and last step stored on the file.
  sw4_type m_nmin, m_nmax;

  // Number of time steps to hold in memory
  sw4_type m_nsteps;

  // Number of time steps currently in memory
  sw4_type m_ncurrent;

  size_t m_ni, m_nj;

  // Name of storage file
  std::string m_filename;

  // Stored data m_data[n] is the n:th time step in memory
  std::vector<double*> m_data;

  // Translate step number in memory to step number in solver
  // i.e., m_steps[n] is the step number of the n:th time step in memory.
  sw4_type* m_steps;

  // Check if it the first time we save to file, need to write header then.
  bool m_startedsave;

  // Error flag.
  bool m_error;
  bool m_isnonempty;

  void save_to_file();
  void read_from_file(sw4_type n);
  void add_patch(sw4_type wind[6]);
  void prsw4_type_openerr(sw4_type ecode) const;

 public:
  DataPatches(std::string fname, Sarray& u, sw4_type imin, sw4_type imax, sw4_type jmin,
              sw4_type jmax, sw4_type kmax, sw4_type nlayers, sw4_type ntsteps, double dt);
  ~DataPatches();
  void push(Sarray& u, sw4_type n);
  void pop(Sarray& u, sw4_type n);
  size_t get_noofpoints() const;
};

#endif
