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
#ifndef EW_SARRAY_H
#define EW_SARRAY_H

#include <sys/types.h>
//#include <iostream>
//#include <mpi.h>
#include "Require.h"
#include <vector>
#include <string>
#include "sw4.h"

using std::string;

class EWCuda;

class Sarray
{
public:
//   Sarray( CartesianProcessGrid* cartcomm, int nc=1 );
   Sarray( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend );
   Sarray( int ibeg, int iend, int jbeg, int jend, int kbeg, int kend );
   Sarray( int nc, int iend, int jend, int kend );
   Sarray( int iend, int jend, int kend );
   Sarray( const Sarray& u );
   Sarray( Sarray& u, int nc=-1 );
   Sarray();
   ~Sarray() {if( m_data != 0 ) delete[] m_data; }
//   void define( CartesianProcessGrid* cartcomm, int nc );
   void define( int iend, int jend, int kend );
   void define( int nc, int iend, int jend, int kend );
   void define( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg,
		int kend );
   void define( int ibeg, int iend, int jbeg, int jend, int kbeg,
		int kend );
   void define( const Sarray& u );
   inline float_sw4* c_ptr() {return m_data;}
#ifdef SW4_CUDA
   __host__ __device__ inline float_sw4* dev_ptr() {return dev_data;}
#else
   inline float_sw4* dev_ptr() {return dev_data;}
#endif
   void reference( float_sw4* new_data ){m_data = new_data; }
   void reference_dev( float_sw4* new_data ){dev_data = new_data; }

   //   inline float_sw4& operator()( int c, int i, int j, int k )
   //   {return m_data[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}
   inline bool in_range( int c, int i, int j, int k )
     {return 1 <= c && c <= m_nc && m_ib <= i && i <= m_ie && m_jb <= j 
	&& j <= m_je && m_kb <= k && k <= m_ke ;}
   inline float_sw4& operator()( int c, int i, int j, int k )
   {
#ifdef BZ_DEBUG
      VERIFY2( in_range(c,i,j,k), "Error Index (c,i,j,k) = (" << c << "," << i << "," << j << "," << k
       << ") not in range 1<= c <= " << m_nc << " " << m_ib << " <= i <= " << m_ie << " " << m_jb 
	       << " <= j <= " << m_je << " " << m_kb << " <=  k <= " << m_ke );
#endif
//      return m_data[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}
      return m_data[m_base+m_offc*c+m_offi*i+m_offj*j+m_offk*k];}
   inline float_sw4& operator()( int i, int j, int k )
      {
#ifdef BZ_DEBUG
         if (!in_range(1,i,j,k))
            VERIFY2(0, 
                    "Error Index (c,i,j,k) = (" << 1 << "," << i << "," << j << "," << k
                    << ") not in range 1<= c <= "
                    << m_nc << " "
                    << m_ib << " <= i <= " << m_ie << " "
                    << m_jb << " <= j <= " << m_je << " "
                    << m_kb << " <= k <= " << m_ke );
            
#endif
//      return m_data[m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}
      return m_data[m_base+m_offi*i+m_offj*j+m_offk*k+m_offc];}
   inline bool is_defined()
      {return m_data != NULL;}
   int m_ib, m_ie, m_jb, m_je, m_kb, m_ke;
   static bool m_corder;
   ssize_t m_base;
   size_t m_offi, m_offj, m_offk, m_offc, m_npts;
//   int index( int i, int j, int k ) {return (i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb);}
   size_t index( int i, int j, int k ) {return m_base+m_offc+m_offi*i+m_offj*j+m_offk*k;}
#ifdef SW4_CUDA
   __host__ __device__ size_t index( int c, int i, int j, int k ) {return m_base+m_offc*c+m_offi*i+m_offj*j+m_offk*k;}
#else
   size_t index( int c, int i, int j, int k ) {return m_base+m_offc*c+m_offi*i+m_offj*j+m_offk*k;}
#endif
   void intersection( int ib, int ie, int jb, int je, int kb, int ke, int wind[6] );
   void side_plane( int side, int wind[6], int nGhost=1 );
   void side_plane_fortran( int side, int wind[6], int nGhost=1 );
   bool in_domain( int i, int j, int k );
   void set_to_zero();
   void checknan(char* tag);
   void set_to_minusOne();
   void set_value( float_sw4 scalar );
   void set_to_random( float_sw4 llim =0.0, float_sw4 ulim = 1.0 );
   void read_from_disk( const char* fname );
   void save_to_disk( const char* fname );
   int ncomp() const {return m_nc;}
   int npts() const  {return m_ni*m_nj*m_nk;}
   void copy( const Sarray& u );
   float_sw4 maximum( int c=1 );
   float_sw4 minimum( int c=1 );
   float_sw4 sum( int c=1 );
   size_t count_nans();
   size_t count_nans( int& cfirst, int& ifirst, int& jfirst, int& kfirst );
   size_t check_match_cpu_gpu( EWCuda* cu, string name );
   size_t check_match_cpu_gpu( EWCuda* cu,  int& cfirst, int& ifirst, int& jfirst, int& kfirst, string name );
   void extract_subarray( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* ar );
   void extract_subarrayIK( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* ar );

   void insert_subarray( int ib, int ie, int jb, int je, int kb, int ke, double* ar );
   void insert_subarray( int ib, int ie, int jb, int je, int kb, int ke, float* ar );
   void insert_subarrayIK( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* ar );
   void copy_kplane( Sarray& u, int k );
   void assign( const float* ar, int corder );
   void assign(  const double* ar, int corder );
   void extract( double* ar, int corder );
   void assign( const float* ar );
   void assign( const double* ar );
   void transposeik();
   void copy_to_device( EWCuda* cu, bool async=false, int st=0 );
   void copy_from_device( EWCuda* cu, bool async=false, int st=0 );
   void allocate_on_device( EWCuda* cu );
   void page_lock( EWCuda* cu );
   void page_unlock( EWCuda* cu );
   Sarray* create_copy_on_device( EWCuda* cu );
   void define_offsets();
   void gaussian_smooth_v1(int width, float decay);
   void gaussian_smooth(int width, float decay);

//   void write( char* filename, CartesianProcessGrid* cartcomm, std::vector<float_sw4> pars );
   int m_nc, m_ni, m_nj, m_nk;
private:
   float_sw4* m_data;
   float_sw4* dev_data;
   inline int min(int i1,int i2){if( i1<i2 ) return i1;else return i2;}
   inline int max(int i1,int i2){if( i1>i2 ) return i1;else return i2;}
//   void init_mpi_datatype( CartesianProcessGrid* cartcomm );
//    bool m_mpi_datatype_initialized;
//    MPI_Datatype m_local_block_type;
};

#endif
