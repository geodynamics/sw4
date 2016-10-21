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
#ifndef EW_WPPPIO_H
#define EW_WPPPIO_H

#include "sw4.h"
#include "Byteswapper.h"

class Comminfo
{
public:
   Comminfo();
   ~Comminfo();
   void print( int recv );
   bool m_has_values;
   int m_steps;
   int*  m_ncomm; 
   int** m_comm_id;
   int** m_comm_index[6];
   size_t  m_maxbuf;
   size_t  m_maxiobuf;
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

class Parallel_IO
{
public:
   Parallel_IO( int iwrite, int pfs, int globalsizes[3], int localsizes[3],
	    int starts[3], int nptsbuf=1000000, int padding=0 );
   void write_array( int* fid, int nc, void* array, off_t pos0, char* type );
   void read_array( int* fid, int nc, float_sw4* array, off_t pos0, const char* typ, bool swap_bytes=false );
			      
   void print( );
   void begin_sequential( MPI_Comm comm );
   void end_sequential( MPI_Comm comm );
   int proc_zero();
   int i_write() const {return m_iwrite==1;}
   void writer_barrier( );
   int n_writers() const {return m_nwriters;}
private:
   void init_pio( int iwrite, int pfs, int ihave_array=-1 );
   void init_array( int globalsizes[3], int localsizes[3], 
		    int starts[3], int nptsbuf, int padding=0 );
   void setup_substeps( );
   int m_iwrite, m_nwriters, m_parallel_file_system;
   int m_csteps;
   int* m_writer_ids;
   int ni, nj, nk, nig, njg, nkg, oi, oj, ok;
   Byteswapper m_bswap;

   MPI_Comm m_write_comm; 
   MPI_Comm m_data_comm;

   Comminfo m_isend;
   Comminfo m_irecv;
};

#endif
