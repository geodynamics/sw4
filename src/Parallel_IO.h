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

#include "mpi.h"
#include "sw4.h"

#include "Byteswapper.h"
#ifdef USE_HDF5
#include "hdf5.h"
#endif

class Comminfo {
 public:
  Comminfo();
  ~Comminfo();
  void print(int recv);
  bool m_has_values;
  int m_steps;
  int* m_ncomm;
  int** m_comm_id;
  int** m_comm_index[6];
  int m_maxbuf;
  int m_maxiobuf;
  int* m_ilow;
  int* m_jlow;
  int* m_klow;
  int* m_niblock;
  int* m_njblock;
  int* m_nkblock;

  // Communication substeps
  int* m_nsubcomm;
  int** m_subcomm;
  int** m_subcommlabel;
};

class Parallel_IO {
 public:
  Parallel_IO(int iwrite, int pfs, int globalsizes[3], int localsizes[3],
              int starts[3], MPI_Comm ewcomm, int nptsbuf = 8000000,
              int padding = 0);
  void write_array(int* fid, int nc, void* array, off_t pos0, char* type);
#ifdef USE_HDF5
  void write_array_hdf5(const char* fname, const char* gname, const char* dname,
                        int nc, void* array, hsize_t pos0, char* type);
#endif
  void read_array(int* fid, int nc, float_sw4* array, off_t pos0,
                  const char* typ, bool swap_bytes = false);

  void print();
  void begin_sequential(MPI_Comm comm);
  void end_sequential(MPI_Comm comm);
  int proc_zero();
  int i_write() const { return m_iwrite == 1; }
  void writer_barrier();
  int n_writers() const { return m_nwriters; }
  int proc_zero_rank_in_comm_world();

 private:
  void init_pio(int iwrite, int pfs, int ihave_array = -1);
  void init_array(int globalsizes[3], int localsizes[3], int starts[3],
                  int nptsbuf, int padding = 0);
  void setup_substeps();
  //   size_t read_dble_wlim( int* fid, double* rbuf, size_t nelem, size_t limit
  //   ); size_t write_dble_wlim( int* fid, double* rbuf, size_t nelem, size_t
  //   limit );
  template <class T>
  size_t read_with_limit(int* fid, T* rbuf, size_t nelem, size_t limit);
  template <class T>
  size_t write_with_limit(int* fid, T* rbuf, size_t nelem, size_t limit);

  int m_iwrite, m_nwriters, m_parallel_file_system;
  int m_csteps;
  int* m_writer_ids;
  int ni, nj, nk, nig, njg, nkg, oi, oj, ok;
  int m_zerorank_in_commworld, m_gproc;

  Byteswapper m_bswap;

  MPI_Comm m_write_comm;
  MPI_Comm m_data_comm;
  MPI_Comm m_ewcomm;
  Comminfo m_isend;
  Comminfo m_irecv;
};

#endif
