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
#include "Sarray.h"

#include <iostream>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>

using namespace std;

#ifdef SW4_CUDA
#include <cuda_runtime.h>
#endif
#include "EWCuda.h"

// Default value 
bool Sarray::m_corder = false;

//-----------------------------------------------------------------------
Sarray::Sarray( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend )
{
   m_nc = nc;
   m_ib = ibeg;
   m_ie = iend;
   m_jb = jbeg;
   m_je = jend;
   m_kb = kbeg;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
//   m_mpi_datatype_initialized = false;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
Sarray::Sarray( int ibeg, int iend, int jbeg, int jend, int kbeg, int kend )
{
   m_nc = 1;
   m_ib = ibeg;
   m_ie = iend;
   m_jb = jbeg;
   m_je = jend;
   m_kb = kbeg;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
//   m_mpi_datatype_initialized = false;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
Sarray::Sarray( int nc, int iend, int jend, int kend )
{
   m_nc = nc;
   m_ib = 1;
   m_ie = iend;
   m_jb = 1;
   m_je = jend;
   m_kb = 1;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
//   m_mpi_datatype_initialized = false;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
Sarray::Sarray( int iend, int jend, int kend )
{
   m_nc = 1;
   m_ib = 1;
   m_ie = iend;
   m_jb = 1;
   m_je = jend;
   m_kb = 1;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
//   m_mpi_datatype_initialized = false;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
Sarray::Sarray()
{
//   m_mpi_datatype_initialized = false;
   m_nc = m_ib = m_ie = m_jb = m_je = m_kb = m_ke = 0;
   m_data = NULL;
   dev_data = NULL;
}

//-----------------------------------------------------------------------
Sarray::Sarray( const Sarray& u )
{
   m_nc = u.m_nc;
   m_ib = u.m_ib;
   m_ie = u.m_ie;
   m_jb = u.m_jb;
   m_je = u.m_je;
   m_kb = u.m_kb;
   m_ke = u.m_ke;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
Sarray::Sarray( Sarray& u, int nc )
{
   if( nc == -1 )
      m_nc = u.m_nc;
   else
      m_nc = nc;
   m_ib = u.m_ib;
   m_ie = u.m_ie;
   m_jb = u.m_jb;
   m_je = u.m_je;
   m_kb = u.m_kb;
   m_ke = u.m_ke;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
//   m_mpi_datatype_initialized = false;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
// void Sarray::define( CartesianProcessGrid* cartcomm, int nc )
// {
//   if( m_data != NULL )
//      delete[] m_data;
//    m_nc = nc;

//    // global index ranges:
//    //   m_ib = cartcomm->gNstart;
//    //   m_ie = cartcomm->getGlobalN()+m_ib-1;
//    //   m_jb = cartcomm->gMstart;
//    //   m_je = cartcomm->getGlobalM()+m_jb-1;
//    //   m_kb = cartcomm->gLstart;
//    //   m_ke = cartcomm->getGlobalL()+m_kb-1;

//    // local index ranges in global coordinates:
//    //   m_ib = cartcomm->Nstart - cartcomm->NlowParallelPadding;
//    //   m_ie = cartcomm->Nstart + cartcomm->getLocalN() - cartcomm->NlowParallelPadding - 1;
//    //   m_jb = cartcomm->Mstart - cartcomm->MlowParallelPadding;
//    //   m_je = cartcomm->Mstart + cartcomm->getLocalM() - cartcomm->MlowParallelPadding - 1;
//    //   m_kb = cartcomm->Lstart - cartcomm->LlowParallelPadding;
//    //   m_ke = cartcomm->Lstart + cartcomm->getLocalL() - cartcomm->LlowParallelPadding - 1;

//    // local index ranges in normalized coordinates
//    m_ib = 1;
//    m_ie = cartcomm->getLocalN();
//    m_jb = 1;
//    m_je = cartcomm->getLocalM();
//    m_kb = 1;
//    m_ke = cartcomm->getLocalL();

//    m_ni = m_ie-m_ib+1;
//    m_nj = m_je-m_jb+1;
//    m_nk = m_ke-m_kb+1;
//    //   std::cout << "Sarray dims " << m_nc << " " << m_ni << " " 
//    //	     << m_nj << " " << m_nk << std::endl;
//    m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
// }

//-----------------------------------------------------------------------
void Sarray::define( int nc, int iend, int jend, int kend )
{
   if( m_data != NULL )
      delete[] m_data;

   m_nc = nc;
   m_ib = 1;
   m_ie = iend;
   m_jb = 1;
   m_je = jend;
   m_kb = 1;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
//   m_mpi_datatype_initialized = false;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
void Sarray::define( int iend, int jend, int kend )
{
   if( m_data != NULL )
      delete[] m_data;

   m_nc = 1;
   m_ib = 1;
   m_ie = iend;
   m_jb = 1;
   m_je = jend;
   m_kb = 1;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
//   m_mpi_datatype_initialized = false;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
void Sarray::define( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg,
		     int kend )
{
   if( m_data != NULL )
      delete[] m_data;
   m_nc = nc;
   m_ib = ibeg;
   m_ie = iend;
   m_jb = jbeg;
   m_je = jend;
   m_kb = kbeg;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
void Sarray::define( int ibeg, int iend, int jbeg, int jend, int kbeg,
		     int kend )
{
   if( m_data != NULL )
      delete[] m_data;
   m_nc = 1;
   m_ib = ibeg;
   m_ie = iend;
   m_jb = jbeg;
   m_je = jend;
   m_kb = kbeg;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
void Sarray::define( const Sarray& u ) 
{
   if( m_data != NULL )
      delete[] m_data;
   m_nc = u.m_nc;
   m_ib = u.m_ib;
   m_ie = u.m_ie;
   m_jb = u.m_jb;
   m_je = u.m_je;
   m_kb = u.m_kb;
   m_ke = u.m_ke;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
   dev_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
void Sarray::intersection( int ib, int ie, int jb, int je, int kb, int ke, int wind[6] )
{
   wind[0] = max(ib,m_ib);
   wind[1] = min(ie,m_ie);
   wind[2] = max(jb,m_jb);
   wind[3] = min(je,m_je);
   wind[4] = max(kb,m_kb);
   wind[5] = min(ke,m_ke);
}

//-----------------------------------------------------------------------
// side_plane returns the index of the ghost points along side =0,1,2,3,4,5 (low-i, high-i, low-j, high-j, low-k, high-k)
void Sarray::side_plane( int side, int wind[6], int nGhost )
{
   wind[0] = m_ib;
   wind[1] = m_ie;
   wind[2] = m_jb;
   wind[3] = m_je;
   wind[4] = m_kb;
   wind[5] = m_ke;
   if( side == 0 )
     wind[1] = wind[0] + (nGhost-1);
   else if( side == 1 )
     wind[0] = wind[1] - (nGhost-1);
   else if( side == 2 )
     wind[3] = wind[2] + (nGhost-1);
   else if( side == 3 )
     wind[2] = wind[3] - (nGhost-1);
   else if( side == 4 )
     wind[5] = wind[4] + (nGhost-1);
   else
     wind[4] = wind[5] - (nGhost-1);
   // if( side == 0 )
   //    wind[1] = wind[0];
   // else if( side == 1 )
   //    wind[0] = wind[1];
   // else if( side == 2 )
   //    wind[3] = wind[2];
   // else if( side == 3 )
   //    wind[2] = wind[3];
   // else if( side == 4 )
   //    wind[5] = wind[4];
   // else
   //    wind[4] = wind[5];
}

//-----------------------------------------------------------------------
void Sarray::side_plane_fortran( int side, int wind[6], int nGhost )
{
// Fortran arrays are base 1
   wind[0] = 1;
   wind[1] = m_ni;
   wind[2] = 1;
   wind[3] = m_nj;
   wind[4] = 1;
   wind[5] = m_nk;
   if( side == 0 )
     wind[1] = wind[0] + (nGhost-1);
   else if( side == 1 )
     wind[0] = wind[1] - (nGhost-1);
   else if( side == 2 )
     wind[3] = wind[2] + (nGhost-1);
   else if( side == 3 )
     wind[2] = wind[3] - (nGhost-1);
   else if( side == 4 )
     wind[5] = wind[4] + (nGhost-1);
   else
     wind[4] = wind[5] - (nGhost-1);
}

//-----------------------------------------------------------------------
void Sarray::set_to_zero()
{
#pragma omp parallel for
   for( size_t i=0 ; i < m_npts ; i++ )
      m_data[i] = 0;
}

//-----------------------------------------------------------------------
void Sarray::set_to_minusOne()
{
#pragma omp parallel for
   for( size_t i=0 ; i < m_npts ; i++ )
      m_data[i] = -1.;
}

//-----------------------------------------------------------------------
void Sarray::set_value( float_sw4 scalar )
{
#pragma omp parallel for
   for( size_t i=0 ; i < m_npts ; i++ )
      m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void Sarray::set_to_random( float_sw4 llim, float_sw4 ulim )
{
   // drand48 is not thread-safe; you will probably not get what you expect
#pragma omp parallel for
   for( size_t i=0 ; i<m_npts ; i++ )
      m_data[i] = llim + (ulim-llim)*drand48();
}

//-----------------------------------------------------------------------
bool Sarray::in_domain( int i, int j, int k )
{
   return m_ib <= i && i <= m_ie && m_jb <= j && j <= m_je
      && m_kb <= k && k <= m_ke;
}

//-----------------------------------------------------------------------
float_sw4 Sarray::maximum( int c )
{
   ///   int cm = c-1;
   //   float_sw4 mx = m_data[cm];
   //   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
   //      mx = mx > m_data[cm+i*m_nc] ? mx : m_data[cm+i*m_nc];
   //   size_t first = m_base+m_offc*c+m_offi*m_ib+m_offj*m_jb+m_offk*m_kb;
   size_t npts = static_cast<size_t>(m_ni)*m_nj*m_nk;
   float_sw4 mx;
   if( m_corder )
   {
      size_t first = (c-1)*npts;
      mx = m_data[first];
#pragma omp parallel for reduction(max:mx)
      for( int i=0 ; i<npts ; i++ )
	 mx = mx > m_data[first+i] ? mx : m_data[first+i];
   }
   else
   {
      size_t first = (c-1);
      mx = m_data[first];
#pragma omp parallel for reduction(max:mx)
      for( int i=0 ; i<npts ; i++ )
	 mx = mx > m_data[first+i*m_nc] ? mx : m_data[first+i*m_nc];
   }
   return mx;
}

//-----------------------------------------------------------------------
float_sw4 Sarray::minimum( int c )
{
   //   int cm = c-1;
   //   float_sw4 mn = m_data[cm];
   //   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
   //      mn = mn < m_data[cm+i*m_nc] ? mn : m_data[cm+i*m_nc];
   size_t npts = static_cast<size_t>(m_ni)*m_nj*m_nk;
   float_sw4 mn;
   if( m_corder )
   {
      size_t first = (c-1)*npts;
      mn = m_data[first];
#pragma omp parallel for reduction(min:mn)
      for( int i=0 ; i<npts ; i++ )
	 mn = mn < m_data[first+i] ? mn : m_data[first+i];
   }
   else
   {
      size_t first = (c-1);
      mn = m_data[first];
#pragma omp parallel for reduction(min:mn)
      for( int i=0 ; i<npts ; i++ )
	 mn = mn < m_data[first+i*m_nc] ? mn : m_data[first+i*m_nc];
   }
   return mn;
}

//-----------------------------------------------------------------------
float_sw4 Sarray::sum( int c )
{
   //   int cm = c-1;
   //   float_sw4 s = 0;
   //   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
   //      s += m_data[cm+i*m_nc];
   //   size_t first = m_base+m_offc*c+m_offi*m_ib+m_offj*m_jb+m_offk*m_kb;
   size_t npts = static_cast<size_t>(m_ni)*m_nj*m_nk;
   float_sw4 s = 0;
   if( m_corder )
   {
      size_t first = (c-1)*npts;
#pragma omp parallel for reduction(+:s)
      for( int i=0 ; i<npts ; i++ )
	 s += m_data[first+i];
   }
   else
   {
      size_t first = (c-1);
#pragma omp parallel for reduction(+:s)
      for( int i=0 ; i<npts ; i++ )
	 s += m_data[first+i*m_nc];
   }
   return s;
}

//-----------------------------------------------------------------------
size_t Sarray::count_nans()
{
   size_t retval = 0;
   size_t npts = m_nc*m_ni*static_cast<size_t>(m_nj)*m_nk;
#pragma omp parallel for reduction(+:retval)
   for( size_t ind = 0; ind < npts ; ind++)
      if( std::isnan(m_data[ind]) )
	 retval++;
   return retval;
}

//-----------------------------------------------------------------------
size_t Sarray::count_nans( int& cfirst, int& ifirst, int& jfirst, int& kfirst )
{
   cfirst = ifirst = jfirst = kfirst = 0;
   size_t retval = 0, ind=0;
   // Note: you're going to get various threads racing to set the "first" values. This won't work.
#pragma omp parallel for reduction(+:retval)
   for( int k=m_kb ; k<=m_ke ; k++ )
      for( int j=m_jb ; j<=m_je ; j++ )
	 for( int i=m_ib ; i <= m_ie ; i++ )
	    for( int c=1 ; c <= m_nc ; c++ )
	    {
	       if( std::isnan(m_data[ind]) )
	       {
		  if( retval == 0 )
		  {
		     ifirst = i;
		     jfirst = j;
		     kfirst = k;
		     cfirst = c;
		  }
		  retval++;
	       }
	       ind++;
	    }
   return retval;
}

//-----------------------------------------------------------------------
void Sarray::copy( const Sarray& u )
{
   if( m_data != NULL )
      delete[] m_data;

   m_nc = u.m_nc;
   m_ib = u.m_ib;
   m_ie = u.m_ie;
   m_jb = u.m_jb;
   m_je = u.m_je;
   m_kb = u.m_kb;
   m_ke = u.m_ke;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
   {
      m_data = new float_sw4[m_nc*m_ni*m_nj*m_nk];
#pragma omp parallel for 
      for( int i=0 ; i < m_nc*m_ni*m_nj*m_nk ; i++ )
	 m_data[i] = u.m_data[i];
   }
   else
      m_data = NULL;
   define_offsets();
}

//-----------------------------------------------------------------------
void Sarray::extract_subarray( int ib, int ie, int jb, int je, int kb,
			       int ke, float_sw4* ar )
{
   // Assuming nc is the same for m_data and subarray ar.
   int nis = ie-ib+1;
   int njs = je-jb+1;
   //   int nks = ke-kb+1;
   size_t sind=0, ind=0;
   for( int k=kb ; k<=ke ; k++ )
      for( int j=jb ; j<=je ; j++ )
	 for( int i=ib ; i <= ie ; i++ )
	    {
               sind = (i-ib)  +  nis*(j-jb)   +  nis*njs*(k-kb);
               ind = (i-m_ib) + m_ni*(j-m_jb) + m_ni*m_nj*(k-m_kb);
	       for( int c=1 ; c <= m_nc ; c++ )
		  ar[sind*m_nc+c-1] = m_data[ind*m_nc+c-1];
	    }
}

//-----------------------------------------------------------------------
void Sarray::insert_subarray( int ib, int ie, int jb, int je, int kb,
			      int ke, double* ar )
{
   // Assuming nc is the same for m_data and subarray ar.
   int nis = ie-ib+1;
   int njs = je-jb+1;
   //   int nks = ke-kb+1;
   size_t sind=0, ind=0;
   for( int k=kb ; k<=ke ; k++ )
      for( int j=jb ; j<=je ; j++ )
	 for( int i=ib ; i <= ie ; i++ )
	    {
               sind = (i-ib)  +  nis*(j-jb)   +  nis*njs*(k-kb);
               ind = (i-m_ib) + m_ni*(j-m_jb) + m_ni*m_nj*(k-m_kb);
	       for( int c=1 ; c <= m_nc ; c++ )
		  m_data[ind*m_nc+c-1] = ar[sind*m_nc+c-1];
	    }
}

//-----------------------------------------------------------------------
void Sarray::insert_subarray( int ib, int ie, int jb, int je, int kb,
			      int ke, float* ar )
{
   // Assuming nc is the same for m_data and subarray ar.
   int nis = ie-ib+1;
   int njs = je-jb+1;
   //   int nks = ke-kb+1;
   size_t sind=0, ind=0;
   for( int k=kb ; k<=ke ; k++ )
      for( int j=jb ; j<=je ; j++ )
	 for( int i=ib ; i <= ie ; i++ )
	    {
               sind = (i-ib) + nis*(j-jb) + nis*njs*(k-kb);
               ind = (i-m_ib) + m_ni*(j-m_jb) + m_ni*m_nj*(k-m_kb);
	       for( int c=1 ; c <= m_nc ; c++ )
		  m_data[ind*m_nc+c-1] = (float_sw4)ar[sind*m_nc+c-1];
	    }
}

//-----------------------------------------------------------------------
void Sarray::save_to_disk( const char* fname )
{
   int fd = open(fname, O_CREAT | O_TRUNC | O_WRONLY, 0660 );
   if( fd == -1 )
      std::cout << "ERROR opening file" << fname << " for writing " << std::endl;
   size_t nr = write(fd,&m_nc,sizeof(int));
   if( nr != sizeof(int) )
      std::cout << "Error saving nc to " << fname << std::endl;
   nr = write(fd,&m_ni,sizeof(int));
   if( nr != sizeof(int) )
      std::cout << "Error saving ni to " << fname << std::endl;
   nr = write(fd,&m_nj,sizeof(int));
   if( nr != sizeof(int) )
      std::cout << "Error saving nj to " << fname << std::endl;
   nr = write(fd,&m_nk,sizeof(int));
   if( nr != sizeof(int) )
      std::cout << "Error saving nk to " << fname << std::endl;
   size_t npts = m_nc*( (size_t)m_ni)*m_nj*( (size_t)m_nk);
   if( m_corder )
   {
      float_sw4* ar = new float_sw4[npts];
      for( int k = 0 ; k < m_nk ; k++ )
	 for( int j = 0 ; j < m_nj ; j++ )
	    for( int i = 0 ; i < m_ni ; i++ )
	       for( int c=0 ; c < m_nc ; c++ )
		  ar[c+m_nc*i+m_nc*m_ni*j+m_nc*m_ni*m_nj*k] = m_data[i+m_ni*j+m_ni*m_nj*k+m_ni*m_nj*m_nk*c];
      nr = write(fd,ar,sizeof(float_sw4)*npts);
      delete[] ar;
   }
   else
      nr = write(fd,m_data,sizeof(float_sw4)*npts);
   if( nr != sizeof(float_sw4)*npts )
      std::cout << "Error saving data array to " << fname << std::endl;
   close(fd);
}

//-----------------------------------------------------------------------
void Sarray::assign( const float_sw4* ar, int corder )
{
   if( corder == m_corder || corder == -1 )
   {
      // Both arrays in the same order
#pragma omp parallel for
      for( size_t i=0 ; i < m_ni*((size_t) m_nj)*m_nk*m_nc ; i++ )
	 m_data[i] = ar[i];
   }
   else if( m_corder )
   {
      // Class array in corder, input array in fortran order, 
#pragma omp parallel for
      for( int i=0 ; i <m_ni ; i++ )
	 for( int j=0 ; j <m_nj ; j++ )
	    for( int k=0 ; k <m_nk ; k++ )
	       for( int c=0 ; c < m_nc ; c++ )
		  m_data[i+m_ni*j+m_ni*m_nj*k+m_ni*m_nj*m_nk*c] = ar[c+m_nc*i+m_nc*m_ni*j+m_nc*m_ni*m_nj*k];
   }
   else
   {
  // Class array in fortran order, input array in corder, 
#pragma omp parallel for
      for( int i=0 ; i <m_ni ; i++ )
	 for( int j=0 ; j <m_nj ; j++ )
	    for( int k=0 ; k <m_nk ; k++ )
	       for( int c=0 ; c < m_nc ; c++ )
		  m_data[c+m_nc*i+m_nc*m_ni*j+m_nc*m_ni*m_nj*k] = ar[i+m_ni*j+m_ni*m_nj*k+m_ni*m_nj*m_nk*c];
   }
}

//-----------------------------------------------------------------------
void Sarray::assign( const float* ar )
{
#pragma omp parallel for
   for( size_t i=0 ; i < m_ni*((size_t) m_nj)*m_nk*m_nc ; i++ )
     m_data[i] = (float_sw4) ar[i];
}

//-----------------------------------------------------------------------
void Sarray::define_offsets()
{
   m_npts = static_cast<size_t>(m_ni)*m_nj*m_nk*m_nc;
   if( m_corder )
   {
      // (i,j,k,c)=i-ib+ni*(j-jb)+ni*nj*(k-kb)+ni*nj*nk*(c-1)
      m_base = -m_ib-m_ni*m_jb-m_ni*m_nj*m_kb-m_ni*m_nj*m_nk;
      m_offc = m_ni*m_nj*m_nk;
      m_offi = 1;
      m_offj = m_ni;
      m_offk = m_ni*m_nj;
      // Can use zero based array internally in class, i.e.,
      // (i,j,k,c) = i + ni*j+ni*nj*k+ni*nj*nk*c
   }
   else
   {
      // (c,i,j,k)=c-1+nc*(i-ib)+nc*ni*(j-jb)+nc*ni*nj*(k-kb)
      m_base = -1-m_nc*m_ib-m_nc*m_ni*m_jb-m_nc*m_ni*m_nj*m_kb;
      m_offc = 1;
      m_offi = m_nc;
      m_offj = m_nc*m_ni;
      m_offk = m_nc*m_ni*m_nj;
      // Can use zero based array internally in class, i.e.,
      // (i,j,k,c) = c + nc*i + nc*ni*j+nc*ni*nj*k
   }
}
//-----------------------------------------------------------------------
void Sarray::transposeik( )
{
   // Transpose a_{i,j,k} := a_{k,j,i}
   float_sw4* tmpar = new float_sw4[m_nc*m_ni*m_nj*m_nk];
   if( m_corder )
   {
      size_t npts = static_cast<size_t>(m_ni)*m_nj*m_nk;
#pragma omp parallel for   
      for( int i=0 ; i <m_ni ; i++ )
	 for( int j=0 ; j <m_nj ; j++ )
	    for( int k=0 ; k <m_nk ; k++ )
	    {
	       size_t ind  = i + m_ni*j + m_ni*m_nj*k;
	       size_t indr = k + m_nk*j + m_nk*m_nj*i;
	       for( int c=0 ; c < m_nc ; c++ )
		  tmpar[npts*c+ind] = m_data[npts*c+indr];
	    }
   }
   else
   {
#pragma omp parallel for   
      for( int i=0 ; i <m_ni ; i++ )
	 for( int j=0 ; j <m_nj ; j++ )
	    for( int k=0 ; k <m_nk ; k++ )
	    {
	       size_t ind  = i + m_ni*j + m_ni*m_nj*k;
	       size_t indr = k + m_nk*j + m_nk*m_nj*i;
	       for( int c=0 ; c < m_nc ; c++ )
		  tmpar[c+m_nc*ind] = m_data[c+m_nc*indr];
	    }
   }
#pragma omp parallel for   
   for( size_t i=0 ; i < m_ni*((size_t) m_nj)*m_nk*m_nc ; i++ )
      m_data[i] = tmpar[i];
   delete[] tmpar;
}

//-----------------------------------------------------------------------
void Sarray::copy_to_device( EWCuda* cu, bool async, int st )
{
#ifdef SW4_CUDA
   if( cu->has_gpu() )
   {
      cudaError_t retcode;
      if( dev_data == NULL )
      {
	 retcode = cudaMalloc( (void**)&dev_data, m_ni*((size_t) m_nj)*m_nk*m_nc*sizeof(float_sw4) );
	 if( retcode != cudaSuccess )
	    cout << "Error Sarray::copy_to_device, cudaMalloc returned " <<
	       cudaGetErrorString(retcode) << endl;
      }
      if( !async )
      {
	 retcode = cudaMemcpy( dev_data, m_data, m_ni*((size_t) m_nj)*m_nk*m_nc*sizeof(float_sw4),
			    cudaMemcpyHostToDevice );
	 if( retcode != cudaSuccess )
	 {
	    cout << "Error Sarray::copy_to_device, cudaMemcpy returned " <<
	       cudaGetErrorString(retcode) << endl;
	    exit(2);
	 }
      }
      else
      {
	 if( st < cu->m_nstream )
	 {
	    retcode = cudaMemcpyAsync( dev_data, m_data, m_ni*((size_t) m_nj)*m_nk*m_nc*sizeof(float_sw4),
				    cudaMemcpyHostToDevice, cu->m_stream[st] );
	    if( retcode != cudaSuccess )
	       cout << "Error Sarray::copy_to_device, cudaMemcpyAsync returned " <<
		  cudaGetErrorString(retcode) << endl;
	 }
	 else
	    cout << "Error Sarray::copy_to_device, stream number " << st << " does not exist " << endl;
      }
   }
#endif
}

//-----------------------------------------------------------------------
void Sarray::copy_from_device( EWCuda* cu, bool async, int st )
{
#ifdef SW4_CUDA
   if( cu->has_gpu() )
   {
      cudaError_t retcode;
      if( dev_data == NULL )
	 cout << "Error Sarray::copy_from_device: Device memory is not allocated " << endl;
      if( !async )
      {
	 retcode = cudaMemcpy( m_data, dev_data, m_ni*((size_t) m_nj)*m_nk*m_nc*sizeof(float_sw4),
			    cudaMemcpyDeviceToHost );
	 if( retcode != cudaSuccess )
	 {
	    cout << "Error Sarray::copy_from_device, cudaMemcpy returned " <<
	       cudaGetErrorString(retcode) << endl;
	    exit(2);
	 }
      }
      else
      {
	 if( st < cu->m_nstream )
	 {
	    retcode =cudaMemcpyAsync( m_data, dev_data, m_ni*((size_t) m_nj)*m_nk*m_nc*sizeof(float_sw4),
				   cudaMemcpyDeviceToHost, cu->m_stream[st] );
	    if( retcode != cudaSuccess )
	       cout << "Error Sarray::copy_from_device, cudaMemcpyAsync returned " <<
		  cudaGetErrorString(retcode) << endl;
	 }
	 else
	    cout << "Error Sarray::copy_from_device, stream number " << st << " does not exist " << endl;
      }
   }
#endif
}

//-----------------------------------------------------------------------
void Sarray::allocate_on_device( EWCuda* cu )
{
#ifdef SW4_CUDA
   if( cu->has_gpu() )
   {
      cudaError_t retcode;
      if( dev_data != NULL )
      {
	 retcode = cudaFree(dev_data);
	 if( retcode != cudaSuccess )
	    cout << "Error Sarray::allocate_on_device, cudaFree returned " <<
	       cudaGetErrorString(retcode) << endl;
      }
      retcode = cudaMalloc( (void**)&dev_data, m_ni*((size_t) m_nj)*m_nk*m_nc*sizeof(float_sw4) );
      if( retcode != cudaSuccess )
	 cout << "Error Sarray::allocate_on_device, cudaMalloc returned " <<
	    cudaGetErrorString(retcode) << endl;
   }
#endif
}

//-----------------------------------------------------------------------
void Sarray::page_lock( EWCuda* cu )
{
#ifdef SW4_CUDA
   if( cu->has_gpu() )
   {
      cudaError_t retcode;
      retcode = cudaHostRegister( m_data, m_ni*((size_t) m_nj)*m_nk*m_nc*sizeof(float_sw4),
				  cudaHostRegisterDefault );
      if( retcode != cudaSuccess )
	 cout << "Error Sarray::page_lock, cudaHostRegister returned " <<
	    cudaGetErrorString(retcode) << endl;
   }
#endif
}

//-----------------------------------------------------------------------
void Sarray::page_unlock( EWCuda* cu )
{
#ifdef SW4_CUDA
   if( cu->has_gpu() )
   {
      cudaError_t retcode;
      retcode = cudaHostUnregister( m_data );
      if( retcode != cudaSuccess )
	 cout << "Error Sarray::page_unlock, cudaHostUnregister returned " <<
	    cudaGetErrorString(retcode) << endl;
   }
#endif
}

//-----------------------------------------------------------------------
size_t Sarray::check_match_cpu_gpu( EWCuda* cu, string name )
{

#ifdef SW4_CUDA

   size_t retval = 0;
   size_t npts = m_nc*m_ni*static_cast<size_t>(m_nj)*m_nk;
   size_t nsize_bytes = npts*sizeof(double);
   double* m_data_test;
   cudaError_t retcode;

   if( npts > 0 )
      m_data_test = new double [npts];
   else
      m_data_test = NULL;

   if(m_data_test == NULL)
   {
      cout  << name << " error: no allocaiton to m_data_test "  << endl;
      exit(-1);
   }

   if(m_data == NULL)
   {
      cout << name << " error: no allocaiton to m_data "  << endl;
      exit(-1);
   }

   if(dev_data == NULL)
   {
       cout << name << " error: no allocaiton to dev_data "  << endl;
       exit(-1);
   }

   retcode = cudaMemcpy( m_data_test, dev_data, nsize_bytes, cudaMemcpyDeviceToHost );
   if( cudaSuccess != retcode )
   {
      if( m_data_test != NULL)
          delete[] m_data_test;
      cout << name << " error Sarray:: check_match_cpu_gpu (*), cudaMemcpy returned " << cudaGetErrorString(retcode) << endl;
      exit(-1);
   }

   for( size_t ind = 0; ind < npts; ind++ )
   {
      if( fabs(m_data[ind] - m_data_test[ind]) >= (1.0e-4) )
          retval++;
   }

   if( m_data_test != NULL)
       delete[] m_data_test;

   return retval;

#endif

}

//-----------------------------------------------------------------------
size_t Sarray::check_match_cpu_gpu( EWCuda* cu, int& cfirst, int& ifirst, int& jfirst, int& kfirst , string name )
{

#ifdef SW4_CUDA

   size_t retval;

   size_t npts = m_nc*m_ni*static_cast<size_t>(m_nj)*m_nk;
   size_t nsize_bytes = npts*sizeof(double);
   double* m_data_test;
   cudaError_t retcode;

   if( npts > 0 )
      m_data_test = new double[npts];
   else
      m_data_test = NULL;

   if(m_data_test == NULL)
   {
       cout << name << " error: no allocaiton to m_data_test "  << endl;
       exit(-1);
   }

   if(m_data == NULL)
   {
       cout << name << " error: no allocaiton to m_data "  << endl;
       exit(-1);
   }

   if(dev_data == NULL)
   {
       cout << name << " error: no allocaiton to dev_data "  << endl;
       exit(-1);
   }

   retcode = cudaMemcpy( m_data_test, dev_data, nsize_bytes, cudaMemcpyDeviceToHost );
   if( retcode != cudaSuccess )
   {
      cout << name << " error Sarray:: check_match_cpu_gpu(*,*,*,*)  :  cudaMemcpy returned" << cudaGetErrorString(retcode) << endl;
      exit(-1);
   }

   cfirst = ifirst = jfirst = kfirst = 0;
   size_t  ind=0;
   retval = 0;
   for( int k=m_kb ; k<=m_ke ; k++ )
          for( int j=m_jb ; j<=m_je ; j++ )
             for( int i=m_ib ; i <= m_ie ; i++ )
                for( int c=1 ; c <= m_nc ; c++ )
                {
                   if( abs(m_data[ind] - m_data_test[ind]) >= (1.0e-4) )
                   {
                      if( retval == 0 )
                      {
                         ifirst = i;
                         jfirst = j;
                         kfirst = k;
                         cfirst = c;
                     }
                     retval++;
                  }
                  ind++;
               }

    delete[] m_data_test;

    return retval;

#endif

}

//-----------------------------------------------------------------------
Sarray* Sarray::create_copy_on_device( EWCuda* cu )
{
#ifdef SW4_CUDA
   copy_to_device(cu);
   Sarray* dev_array;
   cudaError_t retcode = cudaMalloc( (void**)&dev_array,sizeof(Sarray));
   if( retcode != cudaSuccess )
      cout << "Error creating Sarray on device. retval = " << cudaGetErrorString(retcode) << endl;
   retcode = cudaMemcpy(dev_array,this,sizeof(Sarray),cudaMemcpyHostToDevice);
   if( retcode != cudaSuccess )
      cout << "Error create_copy Sarray to device. retval = " << cudaGetErrorString(retcode) << endl;
   return dev_array;
#endif
}
   
