//-*-c++-*-
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
#ifndef EW_WPPPIO_H
#define EW_WPPPIO_H

#include "Byteswapper.h"
#include "mpi.h"
#include "sw4.h"
#ifdef USE_HDF5
#include "hdf5.h"
#endif

class Comminfo {
 public:
  Comminfo();
  ~Comminfo();
  void prsw4_type(sw4_type recv);
  bool m_has_values;
  sw4_type m_steps;
  sw4_type* m_ncomm;
  sw4_type** m_comm_id;
  sw4_type** m_comm_index[6];
  sw4_type m_maxbuf;
  sw4_type m_maxiobuf;
  sw4_type* m_ilow;
  sw4_type* m_jlow;
  sw4_type* m_klow;
  sw4_type* m_niblock;
  sw4_type* m_njblock;
  sw4_type* m_nkblock;

  // Communication substeps
  sw4_type* m_nsubcomm;
  sw4_type** m_subcomm;
  sw4_type** m_subcommlabel;
};

class Parallel_IO {
 public:
  Parallel_IO(sw4_type iwrite, sw4_type pfs, sw4_type globalsizes[3], sw4_type localsizes[3],
              sw4_type starts[3], sw4_type nptsbuf = 8000000, sw4_type padding = 0);
  void write_array(sw4_type* fid, sw4_type nc, void* array, off_t pos0, char* type);
#ifdef USE_HDF5
  void write_array_hdf5(const char* fname, const char* gname, const char* dname,
                        sw4_type nc, void* array, hsize_t pos0, char* type);
#endif
  void read_array(sw4_type* fid, sw4_type nc, float_sw4* array, off_t pos0,
                  const char* typ, bool swap_bytes = false);

  void prsw4_type();
  void begin_sequential(MPI_Comm comm);
  void end_sequential(MPI_Comm comm);
  sw4_type proc_zero();
  sw4_type i_write() const { return m_iwrite == 1; }
  void writer_barrier();
  sw4_type n_writers() const { return m_nwriters; }
  sw4_type proc_zero_rank_in_comm_world();

 private:
  void init_pio(sw4_type iwrite, sw4_type pfs, sw4_type ihave_array = -1);
  void init_array(sw4_type globalsizes[3], sw4_type localsizes[3], sw4_type starts[3],
                  sw4_type nptsbuf, sw4_type padding = 0);
  void setup_substeps();
  //   size_t read_dble_wlim( sw4_type* fid, double* rbuf, size_t nelem, size_t limit
  //   ); size_t write_dble_wlim( sw4_type* fid, double* rbuf, size_t nelem, size_t
  //   limit );
  template <class T>
  size_t read_with_limit(sw4_type* fid, T* rbuf, size_t nelem, size_t limit);
  template <class T>
  size_t write_with_limit(sw4_type* fid, T* rbuf, size_t nelem, size_t limit);

  int m_iwrite;
  sw4_type m_nwriters, m_parallel_file_system;
  sw4_type m_csteps;
  sw4_type* m_writer_ids;
  sw4_type ni, nj, nk, nig, njg, nkg, oi, oj, ok;
  sw4_type m_zerorank_in_commworld;
  Byteswapper m_bswap;

  MPI_Comm m_write_comm;
  MPI_Comm m_data_comm;

  Comminfo m_isend;
  Comminfo m_irecv;
};

#endif
