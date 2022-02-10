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
// # WITHOUT ANY WARRNTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details.
// #
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#include <fcntl.h>
#include <unistd.h>

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "Mspace.h"
#include "Sarray.h"
#include "caliper.h"
#include "foralls.h"
#include "policies.h"
using namespace std;
std::unordered_map<std::string, float_sw4*> Sarray::static_map = {
    {std::string("Initializer"), (float_sw4*)nullptr}};
// Default value
bool Sarray::m_corder = false;
// This allocator keeps the m_data allocation around and re-uses it for all
// subsequent calls. The allocations are deleted in the EW dtor. It reduces
// runtime at the cost of additional memory usage. A memory pool would be a
// better way to do this.
Sarray::Sarray(int nc, int ibeg, int iend, int jbeg, int jend, int kbeg,
               int kend, const char* file, int line) {
  // static size_t total = 0;
  // SW4_MARK_FUNCTION;
  m_nc = nc;
  m_ib = ibeg;
  m_ie = iend;
  m_jb = jbeg;
  m_je = jend;
  m_kb = kbeg;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  if (m_nc * m_ni * m_nj * m_nk > 0) {
#ifndef SW4_USE_UMPIRE
    ostringstream ss;
    ss << file << line << "." << m_nc * m_ni * m_nj * m_nk;
    auto found = static_map.find(ss.str());
    if (found != static_map.end())
      m_data = found->second;
    else {
      // total += m_nc * m_ni * m_nj * m_nk*8;
      // std::cout<<"NEW MEM "<<ss.str()<<" total = "<<total<<"\n"<<std::flush;
      space = Space::Managed;
      m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
      static_map[ss.str()] = m_data;
    }
#else
    space = Space::Managed_temps;
    m_data =
        SW4_NEW(Space::Managed_temps,
                float_sw4[m_nc * m_ni * m_nj *
                          m_nk]);  // THIS NEEDS TO MATCH THE SPACE IN THE DTOR
    // umpire::ResourceManager &rma = umpire::ResourceManager::getInstance();
    // auto allocator = rma.getAllocator("UM_pool_small");
    // m_data=
    // static_cast<float_sw4*>(allocator.allocate(sizeof(float_sw4)*m_nc*m_ni*m_nj*m_nk));
#endif
  } else
    m_data = NULL;
  //   m_mpi_datatype_initialized = false;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
#ifndef SW4_USE_UMPIRE
  static_alloc = true;
#else
  static_alloc = true;
#endif
}

//-----------------------------------------------------------------------
Sarray::Sarray(int nc, int ibeg, int iend, int jbeg, int jend, int kbeg,
               int kend)
    : static_alloc(false) {
  m_nc = nc;
  m_ib = ibeg;
  m_ie = iend;
  m_jb = jbeg;
  m_je = jend;
  m_kb = kbeg;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  //   m_mpi_datatype_initialized = false;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray(int ibeg, int iend, int jbeg, int jend, int kbeg, int kend)
    : static_alloc(false) {
  m_nc = 1;
  m_ib = ibeg;
  m_ie = iend;
  m_jb = jbeg;
  m_je = jend;
  m_kb = kbeg;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  //   m_mpi_datatype_initialized = false;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray(int nc, int iend, int jend, int kend) : static_alloc(false) {
  m_nc = nc;
  m_ib = 1;
  m_ie = iend;
  m_jb = 1;
  m_je = jend;
  m_kb = 1;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  //   m_mpi_datatype_initialized = false;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray(int iend, int jend, int kend) : static_alloc(false) {
  m_nc = 1;
  m_ib = 1;
  m_ie = iend;
  m_jb = 1;
  m_je = jend;
  m_kb = 1;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  //   m_mpi_datatype_initialized = false;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray() : static_alloc(false) {
  //   m_mpi_datatype_initialized = false;
  m_nc = m_ib = m_ie = m_jb = m_je = m_kb = m_ke = 0;
  m_data = NULL;
  dev_data = NULL;
  prefetched = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray(const Sarray& u) : static_alloc(false) {
  m_nc = u.m_nc;
  m_ib = u.m_ib;
  m_ie = u.m_ie;
  m_jb = u.m_jb;
  m_je = u.m_je;
  m_kb = u.m_kb;
  m_ke = u.m_ke;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
  // for HIP use __HIP_DEVICE_COMPILE__
  // #ifdef __CUDA_ARCH__
  //   printf("Creating a device object \n");
  //   bool device_object=true;
  // #endif
}
//-----------------------------------------------------------------------
Sarray::Sarray(const Sarray& u, Space space_in) {
  space = space_in;
  m_nc = u.m_nc;
  m_ib = u.m_ib;
  m_ie = u.m_ie;
  m_jb = u.m_jb;
  m_je = u.m_je;
  m_kb = u.m_kb;
  m_ke = u.m_ke;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(space, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
  if (space == Space::Managed)
    static_alloc = false;
  else if (space == Space::Managed_temps)
    static_alloc = false;  // Used to be true
  else if (space == Space::Host)
    static_alloc = false;
  else
    std::cout
        << "Warning :: Sarrays outside Space::Managed not fully supported \n";
}

//-----------------------------------------------------------------------
Sarray::Sarray(Sarray& u, int nc) : static_alloc(false) {
  if (nc == -1)
    m_nc = u.m_nc;
  else
    m_nc = nc;
  m_ib = u.m_ib;
  m_ie = u.m_ie;
  m_jb = u.m_jb;
  m_je = u.m_je;
  m_kb = u.m_kb;
  m_ke = u.m_ke;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  //   m_mpi_datatype_initialized = false;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
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
//    //   m_ie = cartcomm->Nstart + cartcomm->getLocalN() -
//    cartcomm->NlowParallelPadding - 1;
//    //   m_jb = cartcomm->Mstart - cartcomm->MlowParallelPadding;
//    //   m_je = cartcomm->Mstart + cartcomm->getLocalM() -
//    cartcomm->MlowParallelPadding - 1;
//    //   m_kb = cartcomm->Lstart - cartcomm->LlowParallelPadding;
//    //   m_ke = cartcomm->Lstart + cartcomm->getLocalL() -
//    cartcomm->LlowParallelPadding - 1;

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
void Sarray::define(int nc, int iend, int jend, int kend) {
  if (m_data != NULL) ::operator delete[](m_data, Space::Managed);

  m_nc = nc;
  m_ib = 1;
  m_ie = iend;
  m_jb = 1;
  m_je = jend;
  m_kb = 1;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  //   m_mpi_datatype_initialized = false;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
}

//-----------------------------------------------------------------------
void Sarray::define(int iend, int jend, int kend) {
  if (m_data != NULL) ::operator delete[](m_data, space);

  m_nc = 1;
  m_ib = 1;
  m_ie = iend;
  m_jb = 1;
  m_je = jend;
  m_kb = 1;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  //   m_mpi_datatype_initialized = false;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
}

//-----------------------------------------------------------------------
void Sarray::define(int nc, int ibeg, int iend, int jbeg, int jend, int kbeg,
                    int kend, Space space_in) {
  if (m_data != NULL) ::operator delete[](m_data, space);
  m_nc = nc;
  m_ib = ibeg;
  m_ie = iend;
  m_jb = jbeg;
  m_je = jend;
  m_kb = kbeg;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = space_in;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(space, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
  if (space == Space::Managed)
    static_alloc = false;
  else if (space == Space::Managed_temps)
    static_alloc = false;  // Used to be true
  else if (space == Space::Host)
    static_alloc = false;
  else
    std::cout
        << "ERROR :: Sarrays outside Space::Managed not fully supported \n";
}

//-----------------------------------------------------------------------
void Sarray::define(int ibeg, int iend, int jbeg, int jend, int kbeg, int kend,
                    Space space_in) {
  if (m_data != NULL)
    ::operator delete[](
        m_data, space);  // This is potential bug since we dont know the space
                         // that m_data was originally allocated in. PBUGS
  m_nc = 1;
  m_ib = ibeg;
  m_ie = iend;
  m_jb = jbeg;
  m_je = jend;
  m_kb = kbeg;
  m_ke = kend;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = space_in;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(space, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
  if (space == Space::Managed)
    static_alloc = false;
  else if (space == Space::Managed_temps)
    static_alloc = true;
  else if (space == Space::Host)
    static_alloc = false;
  else
    std::cout
        << "ERROR :: Sarrays outside Space::Managed not fully supported \n";
}

//-----------------------------------------------------------------------
void Sarray::define(const Sarray& u) {
  if (m_data != NULL) ::operator delete[](m_data, space);
  m_nc = u.m_nc;
  m_ib = u.m_ib;
  m_ie = u.m_ie;
  m_jb = u.m_jb;
  m_je = u.m_je;
  m_kb = u.m_kb;
  m_ke = u.m_ke;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0)
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  else
    m_data = NULL;
  dev_data = NULL;
  define_offsets();
  prefetched = false;
}

//-----------------------------------------------------------------------
void Sarray::intersection(int ib, int ie, int jb, int je, int kb, int ke,
                          int wind[6]) {
  SW4_MARK_FUNCTION;
  wind[0] = max(ib, m_ib);
  wind[1] = min(ie, m_ie);
  wind[2] = max(jb, m_jb);
  wind[3] = min(je, m_je);
  wind[4] = max(kb, m_kb);
  wind[5] = min(ke, m_ke);
}

//-----------------------------------------------------------------------
// side_plane returns the index of the ghost points along side =0,1,2,3,4,5
// (low-i, high-i, low-j, high-j, low-k, high-k)
void Sarray::side_plane(int side, int wind[6], int nGhost) {
  wind[0] = m_ib;
  wind[1] = m_ie;
  wind[2] = m_jb;
  wind[3] = m_je;
  wind[4] = m_kb;
  wind[5] = m_ke;
  if (side == 0)
    wind[1] = wind[0] + (nGhost - 1);
  else if (side == 1)
    wind[0] = wind[1] - (nGhost - 1);
  else if (side == 2)
    wind[3] = wind[2] + (nGhost - 1);
  else if (side == 3)
    wind[2] = wind[3] - (nGhost - 1);
  else if (side == 4)
    wind[5] = wind[4] + (nGhost - 1);
  else
    wind[4] = wind[5] - (nGhost - 1);
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
void Sarray::side_plane_fortran(int side, int wind[6], int nGhost) {
  // Fortran arrays are base 1
  wind[0] = 1;
  wind[1] = m_ni;
  wind[2] = 1;
  wind[3] = m_nj;
  wind[4] = 1;
  wind[5] = m_nk;
  if (side == 0)
    wind[1] = wind[0] + (nGhost - 1);
  else if (side == 1)
    wind[0] = wind[1] - (nGhost - 1);
  else if (side == 2)
    wind[3] = wind[2] + (nGhost - 1);
  else if (side == 3)
    wind[2] = wind[3] - (nGhost - 1);
  else if (side == 4)
    wind[5] = wind[4] + (nGhost - 1);
  else
    wind[4] = wind[5] - (nGhost - 1);
}

//-----------------------------------------------------------------------
void Sarray::set_to_zero_async() {
  SW4_MARK_FUNCTION;
  prefetch();
  float_sw4* lm_data = m_data;
  RAJA::forall<DEFAULT_LOOP1_ASYNC>(
      RAJA::RangeSegment(0, m_npts),
      [=] RAJA_DEVICE(size_t i) { lm_data[i] = 0; });
}
//-----------------------------------------------------------------------
void Sarray::set_to_zero() {
  SW4_MARK_FUNCTION;
  prefetch();
  float_sw4* lm_data = m_data;
#ifdef RAJA_ONLY
  RAJA::forall<DEFAULT_LOOP1>(RAJA::RangeSegment(0, m_npts),
                              [=] RAJA_DEVICE(size_t i) { lm_data[i] = 0; });
#else
  forall(0, m_npts, [=] RAJA_DEVICE(size_t i) { lm_data[i] = 0.0; });
#endif
}

//-----------------------------------------------------------------------
void Sarray::set_to_minusOne() {
  SW4_MARK_FUNCTION;

  prefetch();
  float_sw4* lm_data = m_data;

#ifdef PEEKS_GALORE
  SYNC_STREAM;
  SW4_PEEK;
  SYNC_STREAM;
#endif

#ifdef RAJA_ONLY

  RAJA::forall<DEFAULT_LOOP1>(RAJA::RangeSegment(0, m_npts),
                              [=] RAJA_DEVICE(size_t i) { lm_data[i] = -1.; });

#else
  forall(0, m_npts, [=] RAJA_DEVICE(size_t i) {
    lm_data[i] = -1.;
  });  // cuda-memcheck fails here with -M -gpu. 2/26/2021. Works on X86
       // machines
#endif

#ifdef PEEKS_GALORE
  SYNC_STREAM;
  SW4_PEEK;
  SYNC_STREAM;
#endif
}
//-----------------------------------------------------------------------
void Sarray::set_to_minusOneHost() {
  SW4_MARK_FUNCTION;
#pragma omp parallel for
  for (size_t i = 0; i < m_npts; i++) m_data[i] = -1.0;
}

//-----------------------------------------------------------------------
void Sarray::set_valueHost(float_sw4 scalar) {
  SW4_MARK_FUNCTION;
#pragma omp parallel for
  for (size_t i = 0; i < m_npts; i++) m_data[i] = scalar;
}
//-----------------------------------------------------------------------
void Sarray::set_value(float_sw4 scalar) {
  SW4_MARK_FUNCTION;
  float_sw4* lm_data = m_data;
  RAJA::forall<DEFAULT_LOOP1>(
      RAJA::RangeSegment(0, m_npts),
      [=] RAJA_DEVICE(size_t i) { lm_data[i] = scalar; });
}
//-----------------------------------------------------------------------
void Sarray::set_value_async(float_sw4 scalar) {
  SW4_MARK_FUNCTION;
  float_sw4* lm_data = m_data;
  RAJA::forall<DEFAULT_LOOP1_ASYNC>(
      RAJA::RangeSegment(0, m_npts),
      [=] RAJA_DEVICE(size_t i) { lm_data[i] = scalar; });
}

//-----------------------------------------------------------------------
void Sarray::set_to_random(float_sw4 llim, float_sw4 ulim) {
  SW4_MARK_FUNCTION;
  // drand48 is not thread-safe; you will probably not get what you expect
#pragma omp parallel for
  for (size_t i = 0; i < m_npts; i++)
    m_data[i] = llim + (ulim - llim) * drand48();
}

//-----------------------------------------------------------------------
bool Sarray::in_domain(int i, int j, int k) {
  return m_ib <= i && i <= m_ie && m_jb <= j && j <= m_je && m_kb <= k &&
         k <= m_ke;
}

//-----------------------------------------------------------------------
float_sw4 Sarray::maximum(int c) {
  ///   int cm = c-1;
  //   float_sw4 mx = m_data[cm];
  //   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
  //      mx = mx > m_data[cm+i*m_nc] ? mx : m_data[cm+i*m_nc];
  //   size_t first = m_base+m_offc*c+m_offi*m_ib+m_offj*m_jb+m_offk*m_kb;
  size_t npts = static_cast<size_t>(m_ni) * m_nj * m_nk;
  float_sw4 mx;
  if (m_corder) {
    size_t first = (c - 1) * npts;
    mx = m_data[first];
#pragma omp parallel for reduction(max : mx)
    for (int i = 0; i < npts; i++)
      mx = mx > m_data[first + i] ? mx : m_data[first + i];
  } else {
    size_t first = (c - 1);
    mx = m_data[first];
#pragma omp parallel for reduction(max : mx)
    for (int i = 0; i < npts; i++)
      mx = mx > m_data[first + i * m_nc] ? mx : m_data[first + i * m_nc];
  }
  return mx;
}

//-----------------------------------------------------------------------
float_sw4 Sarray::minimum(int c) {
  //   int cm = c-1;
  //   float_sw4 mn = m_data[cm];
  //   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
  //      mn = mn < m_data[cm+i*m_nc] ? mn : m_data[cm+i*m_nc];
  size_t npts = static_cast<size_t>(m_ni) * m_nj * m_nk;
  float_sw4 mn;
  if (m_corder) {
    size_t first = (c - 1) * npts;
    mn = m_data[first];
#pragma omp parallel for reduction(min : mn)
    for (int i = 0; i < npts; i++)
      mn = mn < m_data[first + i] ? mn : m_data[first + i];
  } else {
    size_t first = (c - 1);
    mn = m_data[first];
#pragma omp parallel for reduction(min : mn)
    for (int i = 0; i < npts; i++)
      mn = mn < m_data[first + i * m_nc] ? mn : m_data[first + i * m_nc];
  }
  return mn;
}

//-----------------------------------------------------------------------
float_sw4 Sarray::sum(int c) {
  //   int cm = c-1;
  //   float_sw4 s = 0;
  //   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
  //      s += m_data[cm+i*m_nc];
  //   size_t first = m_base+m_offc*c+m_offi*m_ib+m_offj*m_jb+m_offk*m_kb;
  size_t npts = static_cast<size_t>(m_ni) * m_nj * m_nk;
  float_sw4 s = 0;
  if (m_corder) {
    size_t first = (c - 1) * npts;
#pragma omp parallel for reduction(+ : s)
    for (int i = 0; i < npts; i++) s += m_data[first + i];
  } else {
    size_t first = (c - 1);
#pragma omp parallel for reduction(+ : s)
    for (int i = 0; i < npts; i++) s += m_data[first + i * m_nc];
  }
  return s;
}

//-----------------------------------------------------------------------
size_t Sarray::count_nans() {
  size_t retval = 0;
  size_t npts = m_nc * m_ni * static_cast<size_t>(m_nj) * m_nk;
#pragma omp parallel for reduction(+ : retval)
  for (size_t ind = 0; ind < npts; ind++)
    if (std::isnan(m_data[ind])) retval++;
  return retval;
}

//-----------------------------------------------------------------------
size_t Sarray::count_nans(int& cfirst, int& ifirst, int& jfirst, int& kfirst) {
  cfirst = ifirst = jfirst = kfirst = 0;
  size_t retval = 0, ind = 0;
  // Note: you're going to get various threads racing to set the "first" values.
  // This won't work.
  //#pragma omp parallel for reduction(+:retval)
  if (m_corder) {
    for (int c = 1; c <= m_nc; c++)
      for (int k = m_kb; k <= m_ke; k++)
        for (int j = m_jb; j <= m_je; j++)
          for (int i = m_ib; i <= m_ie; i++) {
            if (std::isnan(m_data[ind])) {
              if (retval == 0) {
                ifirst = i;
                jfirst = j;
                kfirst = k;
                cfirst = c;
              }
              retval++;
            }
            ind++;
          }
  } else {
    for (int k = m_kb; k <= m_ke; k++)
      for (int j = m_jb; j <= m_je; j++)
        for (int i = m_ib; i <= m_ie; i++)
          for (int c = 1; c <= m_nc; c++) {
            if (std::isnan(m_data[ind])) {
              if (retval == 0) {
                ifirst = i;
                jfirst = j;
                kfirst = k;
                cfirst = c;
              }
              retval++;
            }
            ind++;
          }
  }
  return retval;
}

//-----------------------------------------------------------------------
void Sarray::copy(const Sarray& u) {
  if (m_data != NULL) ::operator delete[](m_data, space);

  m_nc = u.m_nc;
  m_ib = u.m_ib;
  m_ie = u.m_ie;
  m_jb = u.m_jb;
  m_je = u.m_je;
  m_kb = u.m_kb;
  m_ke = u.m_ke;
  m_ni = m_ie - m_ib + 1;
  m_nj = m_je - m_jb + 1;
  m_nk = m_ke - m_kb + 1;
  space = Space::Managed;
  if (m_nc * m_ni * m_nj * m_nk > 0) {
    m_data = SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
#pragma omp parallel for
    for (int i = 0; i < m_nc * m_ni * m_nj * m_nk; i++) m_data[i] = u.m_data[i];
  } else
    m_data = NULL;
  define_offsets();
}

//-----------------------------------------------------------------------
void Sarray::extract_subarray(int ib, int ie, int jb, int je, int kb, int ke,
                              float_sw4* ar) {
  // Assuming nc is the same for m_data and subarray ar.
  int nis = ie - ib + 1;
  int njs = je - jb + 1;
  size_t sind = 0, ind = 0;
  if (m_corder) {
    size_t totpts = static_cast<size_t>(m_ni) * m_nj * m_nk;
    size_t totptss = static_cast<size_t>(nis) * njs * (ke - kb + 1);
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (i - ib) + nis * (j - jb) + nis * njs * (k - kb);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            ar[sind + totptss * (c - 1)] = m_data[ind + totpts * (c - 1)];
        }
  } else {
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (i - ib) + nis * (j - jb) + nis * njs * (k - kb);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            ar[sind * m_nc + c - 1] = m_data[ind * m_nc + c - 1];
        }
  }
}

//-----------------------------------------------------------------------
void Sarray::insert_subarray(int ib, int ie, int jb, int je, int kb, int ke,
                             double* ar) {
  // Assuming nc is the same for m_data and subarray ar.
  int nis = ie - ib + 1;
  int njs = je - jb + 1;
  //   int nks = ke-kb+1;
  size_t sind = 0, ind = 0;
  if (m_corder) {
    size_t totpts = static_cast<size_t>(m_ni) * m_nj * m_nk;
    size_t totptss = static_cast<size_t>(nis) * njs * (ke - kb + 1);
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (i - ib) + nis * (j - jb) + nis * njs * (k - kb);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            m_data[ind + totpts * (c - 1)] = ar[sind + totptss * (c - 1)];
        }
  } else {
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (i - ib) + nis * (j - jb) + nis * njs * (k - kb);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            m_data[ind * m_nc + c - 1] = ar[sind * m_nc + c - 1];
        }
  }
}

//-----------------------------------------------------------------------
void Sarray::insert_subarray(int ib, int ie, int jb, int je, int kb, int ke,
                             float* ar) {
  // Assuming nc is the same for m_data and subarray ar.
  int nis = ie - ib + 1;
  int njs = je - jb + 1;
  //   int nks = ke-kb+1;
  size_t sind = 0, ind = 0;
  if (m_corder) {
    size_t totpts = static_cast<size_t>(m_ni) * m_nj * m_nk;
    size_t totptss = static_cast<size_t>(nis) * njs * (ke - kb + 1);
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (i - ib) + nis * (j - jb) + nis * njs * (k - kb);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            m_data[ind + totpts * (c - 1)] =
                (float_sw4)ar[sind + totptss * (c - 1)];
        }
  } else {
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (i - ib) + nis * (j - jb) + nis * njs * (k - kb);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            m_data[ind * m_nc + c - 1] = (float_sw4)ar[sind * m_nc + c - 1];
        }
  }
}

//-----------------------------------------------------------------------
void Sarray::extract_subarrayIK(int ib, int ie, int jb, int je, int kb, int ke,
                                float_sw4* ar) {
  // Assuming nc is the same for m_data and subarray ar.
  // Return `ar' in order suitable for storing array on file.

  //   int nis = ie-ib+1;
  int njs = je - jb + 1;
  int nks = ke - kb + 1;
  size_t sind = 0, ind = 0;
  if (m_corder) {
    size_t totpts = static_cast<size_t>(m_ni) * m_nj * m_nk;
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (k - kb) + nks * (j - jb) + nks * njs * (i - ib);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            ar[sind * m_nc + c - 1] = m_data[ind + totpts * (c - 1)];
        }
  } else {
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (k - kb) + nks * (j - jb) + nks * njs * (i - ib);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            ar[sind * m_nc + c - 1] = m_data[ind * m_nc + c - 1];
        }
  }
}

//-----------------------------------------------------------------------
void Sarray::insert_subarrayIK(int ib, int ie, int jb, int je, int kb, int ke,
                               float_sw4* ar) {
  // Assuming nc is the same for m_data and subarray ar.
  // Insert array `ar', where `ar' is in order suitable for storing array on
  // file.

  //   int nis = ie-ib+1;
  int njs = je - jb + 1;
  int nks = ke - kb + 1;
  //   int nks = ke-kb+1;
  size_t sind = 0, ind = 0;
  if (m_corder) {
    size_t totpts = static_cast<size_t>(m_ni) * m_nj * m_nk;
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (k - kb) + nks * (j - jb) + nks * njs * (i - ib);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            m_data[ind + totpts * (c - 1)] = ar[sind * m_nc + c - 1];
        }
  } else {
    for (int k = kb; k <= ke; k++)
      for (int j = jb; j <= je; j++)
        for (int i = ib; i <= ie; i++) {
          sind = (k - kb) + nks * (j - jb) + nks * njs * (i - ib);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            m_data[ind * m_nc + c - 1] = ar[sind * m_nc + c - 1];
        }
  }
}

//-----------------------------------------------------------------------
void Sarray::copy_kplane(Sarray& u, int k) {
  SW4_MARK_FUNCTION;
  if (!(u.m_ib == m_ib && u.m_ie == m_ie && u.m_jb == m_jb && u.m_je == m_je)) {
    cout << "Sarray::copy_kplane, ERROR arrays must have same (i,j) dimensions"
         << endl;
    return;
  }
  if (m_corder) {
    size_t nijk = m_ni * m_nj * m_nk;
    size_t unijk = u.m_ni * u.m_nj * u.m_nk;

    float_sw4* um_data = u.m_data;
    float_sw4* lm_data = m_data;
    int um_kb = u.m_kb;
    ASSERT_MANAGED(m_data);
    ASSERT_MANAGED(um_data);
    //    int mib = m_ib;
    // int mjb = m_jb;
    int mkb = m_kb;
    int mni = m_ni;
    int mnj = m_nj;
    // int mnc = m_nc;
    // SW4_MARK_BEGIN("CK_PREF");
    // prefetch();
    // u.prefetch();
    // SW4_MARK_END("CK_PREF");
    size_t ind_start = mni * mnj * (k - mkb);
    size_t uind_start = mni * mnj * (k - um_kb);

#if defined(RAJA_ONLY)
    RAJA::RangeSegment c_range(0, m_nc);
    RAJA::RangeSegment j_range(0, m_je - m_jb + 1);
    RAJA::RangeSegment i_range(0, m_ie - m_ib + 1);

    RAJA::kernel<COPY_KPLANE_EXEC_POL>(
        RAJA::make_tuple(c_range, j_range, i_range),
        [=] RAJA_DEVICE(int c, int j, int i) {

#else
    Range<16> I(0, m_ie + 1 - m_ib);
    Range<4> J(0, m_je + 1 - m_jb);
    Range<4> C(0, m_nc);
    forall3(I, J, C, [=] RAJA_DEVICE(int i, int j, int c) {
    // forall2(I, J, [=] RAJA_DEVICE(int i, int j ) {

#endif
          // for( int c=0 ; c < m_nc ; c++ )
          // 	 for( int j=m_jb ; j<=m_je ; j++ )
          // 	    for( int i=m_ib ; i <= m_ie ; i++ )
          // 	    {
          size_t ind = (i) + mni * (j) + ind_start;    // mni*mnj*(k-mkb);
          size_t uind = (i) + mni * (j) + uind_start;  // mni*mnj*(k-um_kb);

          // size_t ind = i+ mni*j + ind_start; // mni*mnj*(k-mkb);
          // size_t uind = i + mni*j + uind_start; // mni*mnj*(k-um_kb);
          // for (int c=0;c<3;c++)
          lm_data[ind + c * nijk] = um_data[uind + c * unijk];
        });  // SYNC_STREAM;
  } else {
    SW4_MARK_BEGIN("RUNNING ON HOST");
    for (int j = m_jb; j <= m_je; j++)
      for (int i = m_ib; i <= m_ie; i++) {
        size_t ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
        size_t uind =
            (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - u.m_kb);
        for (int c = 0; c < m_nc; c++)
          m_data[c + m_nc * ind] = u.m_data[c + m_nc * uind];
      }
    SW4_MARK_END("RUNNING ON HOST");
  }
}

//-----------------------------------------------------------------------
void Sarray::save_to_disk(const char* fname) {
  int fd = open(fname, O_CREAT | O_TRUNC | O_WRONLY, 0660);
  if (fd == -1)
    std::cout << "ERROR opening file" << fname << " for writing " << std::endl;
  size_t nr = write(fd, &m_nc, sizeof(int));
  if (nr != sizeof(int))
    std::cout << "Error saving nc to " << fname << std::endl;
  nr = write(fd, &m_ni, sizeof(int));
  if (nr != sizeof(int))
    std::cout << "Error saving ni to " << fname << std::endl;
  nr = write(fd, &m_nj, sizeof(int));
  if (nr != sizeof(int))
    std::cout << "Error saving nj to " << fname << std::endl;
  nr = write(fd, &m_nk, sizeof(int));
  if (nr != sizeof(int))
    std::cout << "Error saving nk to " << fname << std::endl;
  size_t npts = m_nc * ((size_t)m_ni) * m_nj * ((size_t)m_nk);
  if (m_corder) {
    float_sw4* ar = SW4_NEW(Space::Host, float_sw4[npts]);
    for (int k = 0; k < m_nk; k++)
      for (int j = 0; j < m_nj; j++)
        for (int i = 0; i < m_ni; i++)
          for (int c = 0; c < m_nc; c++)
            ar[c + m_nc * i + m_nc * m_ni * j + m_nc * m_ni * m_nj * k] =
                m_data[i + m_ni * j + m_ni * m_nj * k + m_ni * m_nj * m_nk * c];
    nr = write(fd, ar, sizeof(float_sw4) * npts);
    delete[] ar;
  } else
    nr = write(fd, m_data, sizeof(float_sw4) * npts);
  if (nr != sizeof(float_sw4) * npts)
    std::cout << "Error saving data array to " << fname << std::endl;
  close(fd);
}

//-----------------------------------------------------------------------
void Sarray::assign(const float* ar, int corder) {
  std::cout << "WARNING:: Float version of Sarray::assign not offloaded \n";
  if (corder == m_corder || corder == -1) {
    // Both arrays in the same order
#pragma omp parallel for
    for (size_t i = 0; i < m_ni * ((size_t)m_nj) * m_nk * m_nc; i++)
      m_data[i] = ar[i];
  } else if (m_corder) {
    // Class array in corder, input array in fortran order,
#pragma omp parallel for
    for (int i = 0; i < m_ni; i++)
      for (int j = 0; j < m_nj; j++)
        for (int k = 0; k < m_nk; k++)
          for (int c = 0; c < m_nc; c++)
            m_data[i + m_ni * j + m_ni * m_nj * k + m_ni * m_nj * m_nk * c] =
                ar[c + m_nc * i + m_nc * m_ni * j + m_nc * m_ni * m_nj * k];
  } else {
    // Class array in fortran order, input array in corder,
#pragma omp parallel for
    for (int i = 0; i < m_ni; i++)
      for (int j = 0; j < m_nj; j++)
        for (int k = 0; k < m_nk; k++)
          for (int c = 0; c < m_nc; c++)
            m_data[c + m_nc * i + m_nc * m_ni * j + m_nc * m_ni * m_nj * k] =
                ar[i + m_ni * j + m_ni * m_nj * k + m_ni * m_nj * m_nk * c];
  }
}

//-----------------------------------------------------------------------
void Sarray::assign(const double* ar, int corder) {
  SW4_MARK_FUNCTION;
  if (corder == m_corder || corder == -1) {
    // Both arrays in the same order
    // #pragma omp parallel for
    //       for( size_t i=0 ; i < m_ni*((size_t) m_nj)*m_nk*m_nc ; i++ )
    const size_t end = m_ni * ((size_t)m_nj) * m_nk * m_nc;
    auto mdata = m_data;
    RAJA::forall<SARRAY_LOOP_POL1>(
        RAJA::RangeSegment(0, end),
        [=] RAJA_DEVICE(size_t i) { mdata[i] = ar[i]; });
  } else if (m_corder) {
    // Class array in corder, input array in fortran order,
    // #pragma omp parallel for
    //       for( int i=0 ; i <m_ni ; i++ )
    // 	 for( int j=0 ; j <m_nj ; j++ )
    // 	    for( int k=0 ; k <m_nk ; k++ )
    // 	       for( int c=0 ; c < m_nc ; c++ )

    // RAJA::RangeSegment i_range(0,m_ni);
    // RAJA::RangeSegment j_range(0,m_nj);
    // RAJA::RangeSegment k_range(0,m_nk);
    // RAJA::nested::forall(SARRAY_LOOP_POL2{},
    // 				      RAJA::make_tuple(i_range,j_range,k_range),
    // 				      [=]RAJA_DEVICE (int i,int j, int k) {
    // 					for(int c=0;c<m_nc;c++)
    // 					m_data[i+m_ni*j+m_ni*m_nj*k+m_ni*m_nj*m_nk*c]
    // = ar[c+m_nc*i+m_nc*m_ni*j+m_nc*m_ni*m_nj*k];});
    int mni = m_ni;
    int mnj = m_nj;
    int mnk = m_nk;
    int mnc = m_nc;
    ASSERT_MANAGED(m_data);
    ASSERT_MANAGED((void*)ar);
    prefetch();
    PREFETCH(ar);
    float_sw4* mdata = m_data;

    RAJA::RangeSegment i_range(0, mni);
    RAJA::RangeSegment j_range(0, mnj);
    RAJA::RangeSegment k_range(0, mnk);
    RAJA::RangeSegment c_range(0, mnc);
    RAJA::kernel<SAA_POL>(
        RAJA::make_tuple(c_range, k_range, j_range, i_range),
        [=] RAJA_DEVICE(int c, int k, int j, int i) {
          mdata[i + mni * j + mni * mnj * k + mni * mnj * mnk * c] =
              ar[c + mnc * i + mnc * mni * j + mnc * mni * mnj * k];
        });  // SYNC_STREAM;
  } else {
    // Class array in fortran order, input array in corder,
    // #pragma omp parallel for
    //       for( int i=0 ; i <m_ni ; i++ )
    // 	 for( int j=0 ; j <m_nj ; j++ )
    // 	    for( int k=0 ; k <m_nk ; k++ )
    // 	       for( int c=0 ; c < m_nc ; c++ )
    float_sw4* mdata = m_data;
    int mni = m_ni;
    int mnj = m_nj;
    int mnk = m_nk;
    int mnc = m_nc;
    RAJA::RangeSegment i_range(0, m_ni);
    RAJA::RangeSegment j_range(0, m_nj);
    RAJA::RangeSegment k_range(0, m_nk);
    RAJA::RangeSegment c_range(0, m_nc);
    RAJA::kernel<SARRAY_LOOP_POL2>(
        RAJA::make_tuple(i_range, j_range, k_range, c_range),
        [=] RAJA_DEVICE(int i, int j, int k, int c) {
          mdata[c + mnc * i + mnc * mni * j + mnc * mni * mnj * k] =
              ar[i + mni * j + mni * mnj * k + mni * mnj * mnk * c];
        });  // SYNC_STREAM;
  }
}

//-----------------------------------------------------------------------
void Sarray::extract(double* ar, int corder) {
  if (corder == m_corder || corder == -1) {
    // Both arrays in the same order
#pragma omp parallel for
    for (size_t i = 0; i < m_ni * ((size_t)m_nj) * m_nk * m_nc; i++)
      ar[i] = m_data[i];
  } else if (m_corder) {
    // Class array in corder, ar array in fortran order,
#pragma omp parallel for
    for (int i = 0; i < m_ni; i++)
      for (int j = 0; j < m_nj; j++)
        for (int k = 0; k < m_nk; k++)
          for (int c = 0; c < m_nc; c++)
            ar[c + m_nc * i + m_nc * m_ni * j + m_nc * m_ni * m_nj * k] =
                m_data[i + m_ni * j + m_ni * m_nj * k + m_ni * m_nj * m_nk * c];
  } else {
    // Class array in fortran order, ar array in corder,
#pragma omp parallel for
    for (int i = 0; i < m_ni; i++)
      for (int j = 0; j < m_nj; j++)
        for (int k = 0; k < m_nk; k++)
          for (int c = 0; c < m_nc; c++)
            ar[i + m_ni * j + m_ni * m_nj * k + m_ni * m_nj * m_nk * c] =
                m_data[c + m_nc * i + m_nc * m_ni * j + m_nc * m_ni * m_nj * k];
  }
}

//-----------------------------------------------------------------------
void Sarray::assign(const float* ar) {
#pragma omp parallel for
  for (size_t i = 0; i < m_ni * ((size_t)m_nj) * m_nk * m_nc; i++)
    m_data[i] = (float_sw4)ar[i];
}

//-----------------------------------------------------------------------
void Sarray::assign(const double* ar) {
#pragma omp parallel for
  for (size_t i = 0; i < m_ni * ((size_t)m_nj) * m_nk * m_nc; i++)
    m_data[i] = (float_sw4)ar[i];
}

//-----------------------------------------------------------------------
void Sarray::define_offsets() {
  m_npts = static_cast<size_t>(m_ni) * m_nj * m_nk * m_nc;
  if (m_corder) {
    // (i,j,k,c)=i-ib+ni*(j-jb)+ni*nj*(k-kb)+ni*nj*nk*(c-1)
    m_base = -m_ib - m_ni * m_jb - m_ni * m_nj * m_kb - m_ni * m_nj * m_nk;
    m_offc = m_ni * m_nj * m_nk;
    m_offi = 1;
    m_offj = m_ni;
    m_offk = m_ni * m_nj;
    // Can use zero based array internally in class, i.e.,
    // (i,j,k,c) = i + ni*j+ni*nj*k+ni*nj*nk*c
  } else {
    // (c,i,j,k)=c-1+nc*(i-ib)+nc*ni*(j-jb)+nc*ni*nj*(k-kb)
    m_base = -1 - m_nc * m_ib - m_nc * m_ni * m_jb - m_nc * m_ni * m_nj * m_kb;
    m_offc = 1;
    m_offi = m_nc;
    m_offj = m_nc * m_ni;
    m_offk = m_nc * m_ni * m_nj;
    // Can use zero based array internally in class, i.e.,
    // (i,j,k,c) = c + nc*i + nc*ni*j+nc*ni*nj*k
  }
  view.set(*this);
}
//-----------------------------------------------------------------------
void Sarray::transposeik() {
  // Transpose a_{i,j,k} := a_{k,j,i}
  float_sw4* tmpar =
      SW4_NEW(Space::Managed, float_sw4[m_nc * m_ni * m_nj * m_nk]);
  if (m_corder) {
    size_t npts = static_cast<size_t>(m_ni) * m_nj * m_nk;
#pragma omp parallel for
    for (int i = 0; i < m_ni; i++)
      for (int j = 0; j < m_nj; j++)
        for (int k = 0; k < m_nk; k++) {
          size_t ind = i + m_ni * j + m_ni * m_nj * k;
          size_t indr = k + m_nk * j + m_nk * m_nj * i;
          for (int c = 0; c < m_nc; c++)
            tmpar[npts * c + ind] = m_data[npts * c + indr];
        }
  } else {
#pragma omp parallel for
    for (int i = 0; i < m_ni; i++)
      for (int j = 0; j < m_nj; j++)
        for (int k = 0; k < m_nk; k++) {
          size_t ind = i + m_ni * j + m_ni * m_nj * k;
          size_t indr = k + m_nk * j + m_nk * m_nj * i;
          for (int c = 0; c < m_nc; c++)
            tmpar[c + m_nc * ind] = m_data[c + m_nc * indr];
        }
  }
#pragma omp parallel for
  for (size_t i = 0; i < m_ni * ((size_t)m_nj) * m_nk * m_nc; i++)
    m_data[i] = tmpar[i];
  //::operator delete[] tmpar;
  ::operator delete[](tmpar, Space::Managed);
}

void Sarray::prefetch(int device) {
#if defined(DISABLE_PREFETCH)
  return;
#else
#if defined(ENABLE_CUDA)
  if (!prefetched) {
    SW4_MARK_BEGIN("PREFETCH");
    SW4_CheckDeviceError(cudaMemPrefetchAsync(
        m_data, m_nc * m_ni * m_nj * m_nk * sizeof(float_sw4), device, 0));
    SW4_MARK_END("PREFETCH");
    prefetched = true;
  }
#endif

#endif  // #if defined(DISABLE_PREFETCH)
}
void Sarray::forceprefetch(int device) {
#if defined(DISABLE_PREFETCH2)
  return;
#else
#if defined(ENABLE_CUDA)
  SW4_MARK_BEGIN("FORCE_PREFETCH");
  SW4_CheckDeviceError(cudaMemPrefetchAsync(
      m_data, m_nc * m_ni * m_nj * m_nk * sizeof(float_sw4), device, 0));
  SW4_MARK_END("FORCE_PREFETCH");
  prefetched = true;

#endif

#endif  // #if defined(DISABLE_PREFETCH)
}
SView::SView() : data{NULL} {
  // std:cerr<<"SVIEW default ctor, should never be called \n";
}
SView::SView(float_sw4* data, ssize_t base, size_t offc, size_t offi,
             size_t offj, size_t offk)
    : data{data}, base{base}, offc{offc}, offi{offi}, offj{offj}, offk{offk} {}
SView::SView(Sarray& x) {
  data = x.c_ptr();
  base = x.m_base;
  offc = x.m_offc;
  offi = x.m_offi;
  offj = x.m_offj;
  offk = x.m_offk;
  // std::cout<<"Sview created with "<<data<<" base = "<<base<<" offc =
  // "<<offc<<" offi = "<<offi<<" offj = "<<offj<<" offk = "<<offk<<"\n";
}
void SView::set(Sarray& x) {
  data = x.c_ptr();
  base = x.m_base;
  offc = x.m_offc;
  offi = x.m_offi;
  offj = x.m_offj;
  offk = x.m_offk;
  // std::cout<<"Sview created with "<<data<<" base = "<<base<<" offc =
  // "<<offc<<" offi = "<<offi<<" offj = "<<offj<<" offk = "<<offk<<"\n";
}
void SarrayVectorPrefetch(vector<Sarray>& v) {
  SW4_MARK_FUNCTION;
  for (int i = 0; i < v.size(); i++) v[i].prefetch();
}
void SarrayVectorPrefetch(vector<Sarray*>& v) {
  SW4_MARK_FUNCTION;
  for (int i = 0; i < v.size(); i++) v[i]->prefetch();
}
void SarrayVectorPrefetch(vector<Sarray*>& v, int n) {
  SW4_MARK_FUNCTION;
  for (int i = 0; i < v.size(); i++)
    for (int j = 0; j < n; j++) v[i][j].prefetch();
}

typedef std::tuple<uintptr_t, int, int, int, int> mkey_t;

struct key_hash : public std::unary_function<mkey_t, std::size_t> {
  std::size_t operator()(const mkey_t& k) const {
    return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k) ^ std::get<3>(k) ^
           std::get<4>(k);
  }
};
float_sw4* memoize(Sarray& u, int c, int i, int j, int k) {
  static std::unordered_map<std::tuple<uintptr_t, int, int, int, int>,
                            float_sw4*, key_hash>
      map;
  static int fc = 0;
  static int nfc = 0;
  auto index = std::make_tuple((uintptr_t)u.c_ptr(), c, i, j, k);
  auto found = map.find(index);
  // std::cout<<" Stats "<<fc<<" "<<nfc<<"\n";
  if (found != map.end()) {
    fc++;
    return found->second;
  } else {
    auto retval = &u(c, i, j, k);
    map[index] = retval;
    nfc++;
    return retval;
  }
}
//-----------------------------------------------------------------------
void Sarray::insert_intersection(Sarray& a_U) {
  SW4_MARK_FUNCTION;
  // Assuming nc is the same for m_data and a_U.m_data.
  int wind[6];
  int ib = a_U.m_ib, ie = a_U.m_ie, jb = a_U.m_jb, je = a_U.m_je, kb = a_U.m_kb,
      ke = a_U.m_ke;
  intersection(ib, ie, jb, je, kb, ke, wind);
  int nis = ie - ib + 1;
  int njs = je - jb + 1;
  int nks = ke - kb + 1;

  if (m_corder) {
    const size_t totpts = static_cast<size_t>(m_ni) * m_nj * m_nk;
    const size_t totptss = static_cast<size_t>(nis) * njs * (nks);
    float_sw4* dst_m_data = m_data;
    float_sw4* src_m_data = a_U.m_data;
    // for (int k = wind[4]; k <= wind[5]; k++)
    //   for (int j = wind[2]; j <= wind[3]; j++)
    //     for (int i = wind[0]; i <= wind[1]; i++) {
    const int lm_ib = m_ib;
    const int lm_jb = m_jb;
    const int lm_kb = m_kb;
    const int lm_ni = m_ni;
    const int lm_nj = m_nj;
    // const int lm_nk = m_nk;
    const int lm_nc = m_nc;
    // std::cout<<"Calling interest \n"<<std::flush;

#if !defined(RAJA_ONLY)

    Range<16> I(wind[0], wind[1] + 1);
    Range<16> J(wind[2], wind[3] + 1);
    Range<1> K(wind[4], wind[5] + 1);

#pragma forceinline
    forall3async(I, J, K, [=] RAJA_DEVICE(int i, int j, int k) {
#else

    RAJA::RangeSegment k_range(wind[4], wind[5] + 1);
    RAJA::RangeSegment j_range(wind[2], wind[3] + 1);
    RAJA::RangeSegment i_range(wind[0], wind[1] + 1);
    RAJA::kernel<SII_POL>(
        RAJA::make_tuple(k_range, j_range, i_range),
        [=] RAJA_DEVICE(int k, int j, int i) {
#endif
      size_t sind = (i - ib) + nis * (j - jb) + nis * njs * (k - kb);
      size_t ind =
          (i - lm_ib) + lm_ni * (j - lm_jb) + lm_ni * lm_nj * (k - lm_kb);

      for (int c = 1; c <= lm_nc; c++)
        dst_m_data[ind + totpts * (c - 1)] =
            src_m_data[sind + totptss * (c - 1)];
    });
    // SYNC_STREAM;
  } else {
    size_t sind = 0, ind = 0;
    for (int k = wind[4]; k <= wind[5]; k++)
      for (int j = wind[2]; j <= wind[3]; j++)
        for (int i = wind[0]; i <= wind[1]; i++) {
          sind = (i - ib) + nis * (j - jb) + nis * njs * (k - kb);
          ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
          for (int c = 1; c <= m_nc; c++)
            m_data[ind * m_nc + c - 1] = a_U.m_data[sind * m_nc + c - 1];
        }
  }
}

void Sarray::extrapolij(int npts) {
  // Extrapolate to layer npts thick at outermost points in i- and j-directions.
  if (m_corder) {
    size_t nijk = static_cast<size_t>(m_ni) * m_nj * m_nk;
    for (int c = 0; c < m_nc; c++) {
      for (int k = 0; k <= m_nk - 1; k++)
        for (int j = 0; j <= m_nj - 1; j++)
          for (int i = 0; i <= npts - 1; i++) {
            size_t ind = i + m_ni * j + m_ni * m_nj * k;
            size_t indf = npts + m_ni * j + m_ni * m_nj * k;
            m_data[ind + c * nijk] = m_data[indf + c * nijk];
          }
      // side i-high
      for (int k = 0; k <= m_nk - 1; k++)
        for (int j = 0; j <= m_nj - 1; j++)
          for (int i = m_ni - npts; i <= m_ni - 1; i++) {
            size_t ind = i + m_ni * j + m_ni * m_nj * k;
            size_t indf = m_ni - npts - 1 + m_ni * j + m_ni * m_nj * k;
            m_data[ind + c * nijk] = m_data[indf + c * nijk];
          }
      // side j-low
      for (int k = 0; k <= m_nk - 1; k++)
        for (int j = 0; j <= npts - 1; j++)
          for (int i = 0; i <= m_ni - 1; i++) {
            size_t ind = i + m_ni * j + m_ni * m_nj * k;
            size_t indf = i + m_ni * npts + m_ni * m_nj * k;
            m_data[ind + c * nijk] = m_data[indf + c * nijk];
          }
      // side j-high
      for (int k = 0; k <= m_nk - 1; k++)
        for (int j = m_nj - npts; j <= m_nj - 1; j++)
          for (int i = 0; i <= m_ni - 1; i++) {
            size_t ind = i + m_ni * j + m_ni * m_nj * k;
            size_t indf = i + m_ni * (m_nj - npts - 1) + m_ni * m_nj * k;
            m_data[ind + c * nijk] = m_data[indf + c * nijk];
          }
    }
  } else {
    // side i-low
    for (int k = 0; k <= m_nk - 1; k++)
      for (int j = 0; j <= m_nj - 1; j++)
        for (int i = 0; i <= npts - 1; i++)
          for (int c = 0; c < m_nc; c++) {
            size_t ind = i + m_ni * j + m_ni * m_nj * k;
            size_t indf = npts + m_ni * j + m_ni * m_nj * k;
            m_data[c + m_nc * ind] = m_data[c + m_nc * indf];
          }
    // side i-high
    for (int k = 0; k <= m_nk - 1; k++)
      for (int j = 0; j <= m_nj - 1; j++)
        for (int i = m_ni - npts; i <= m_ni - 1; i++)
          for (int c = 0; c < m_nc; c++) {
            size_t ind = i + m_ni * j + m_ni * m_nj * k;
            size_t indf = m_ni - 1 - npts + m_ni * j + m_ni * m_nj * k;
            m_data[c + m_nc * ind] = m_data[c + m_nc * indf];
          }
    // side j-low
    for (int k = 0; k <= m_nk - 1; k++)
      for (int j = 0; j <= npts - 1; j++)
        for (int i = 0; i <= m_ni - 1; i++)
          for (int c = 0; c < m_nc; c++) {
            size_t ind = i + m_ni * j + m_ni * m_nj * k;
            size_t indf = i + m_ni * npts + m_ni * m_nj * k;
            m_data[c + m_nc * ind] = m_data[c + m_nc * indf];
          }
    // side j-high
    for (int k = 0; k <= m_nk - 1; k++)
      for (int j = m_nj - npts; j <= m_nj - 1; j++)
        for (int i = 0; i <= m_ni - 1; i++)
          for (int c = 0; c < m_nc; c++) {
            size_t ind = i + m_ni * j + m_ni * m_nj * k;
            size_t indf = i + m_ni * (m_nj - 1 - npts) + m_ni * m_nj * k;
            m_data[c + m_nc * ind] = m_data[c + m_nc * indf];
          }
  }
}
//-----------------------------------------------------------------------
void Sarray::copy_kplane2(Sarray& u, int k) {
  // Only check k-dimension, other dims do not have to match, only copy the
  // intersecting part.
  SW4_MARK_FUNCTION;
  if (!(u.m_kb <= k && k <= u.m_ke && m_kb <= k && k <= m_ke)) {
    cout << "Sarray::copy_kplane, ERROR k index " << k << " not in range "
         << endl;
    return;
  }
  int wind[6];
  intersection(u.m_ib, u.m_ie, u.m_jb, u.m_je, u.m_kb, u.m_ke, wind);
  if (m_corder) {
    int lm_nc = m_nc;
    int lm_ib = m_ib;
    int lm_ni = m_ni;
    int lm_jb = m_jb;
    int lm_nj = m_nj;
    int lm_kb = m_kb;

    int ulm_ib = u.m_ib;
    int ulm_ni = u.m_ni;
    int ulm_jb = u.m_jb;
    int ulm_nj = u.m_nj;
    int ulm_kb = u.m_kb;

    float_sw4* lm_data = m_data;
    float_sw4* ulm_data = u.m_data;
    size_t nijk = m_ni * m_nj * m_nk;
    size_t unijk = u.m_ni * u.m_nj * u.m_nk;
    // for (int c = 0; c < m_nc; c++)
    //   for (int j = wind[2]; j <= wind[3]; j++)
    //     for (int i = wind[0]; i <= wind[1]; i++) {
    RAJA::RangeSegment j_range(wind[2], wind[3] + 1);
    RAJA::RangeSegment i_range(wind[0], wind[1] + 1);
    RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
        RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
          size_t ind =
              (i - lm_ib) + lm_ni * (j - lm_jb) + lm_ni * lm_nj * (k - lm_kb);

          size_t uind = (i - ulm_ib) + ulm_ni * (j - ulm_jb) +
                        ulm_ni * ulm_nj * (k - ulm_kb);

          for (int c = 0; c < lm_nc; c++)
            lm_data[ind + c * nijk] = ulm_data[uind + c * unijk];
        });
    // SYNC_STREAM;
  } else {
    std::cout << "WARNING Sarray::copy_kplane2 running on CPU !!\n"
              << std::flush;
    for (int j = wind[2]; j <= wind[3]; j++)
      for (int i = wind[0]; i <= wind[1]; i++) {
        size_t ind = (i - m_ib) + m_ni * (j - m_jb) + m_ni * m_nj * (k - m_kb);
        size_t uind = (i - u.m_ib) + u.m_ni * (j - u.m_jb) +
                      u.m_ni * u.m_nj * (k - u.m_kb);
        for (int c = 0; c < m_nc; c++)
          m_data[c + m_nc * ind] = u.m_data[c + m_nc * uind];
      }
  }
}
//----------------------------------------------------------------------
void Sarray::swrite(std::string filename) {
  if (!of.is_open()) {
    int myRank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    std::stringstream ss;
    ss << filename << "_" << myRank << ".dat";
    std::cout << "Opening " << ss.str() << "\n";
    of.open(ss.str());
  } else {
    std::cout << "File is open\n";
  }

  of << " \n\n++++ " << filename << " \n\n";
  for (int i = 0; i < m_nc * m_ni * m_nj * m_nk; i++) {
    of << i << " " << m_data[i] << "\n";
  }
}
//----------------------------------------------------------------------
void Sarray::GetAtt(char* file, int line) {
  short int data = -999;
#ifdef ENABLE_CUDA
  if (cuMemRangeGetAttribute(
          &data, 4, CU_MEM_RANGE_ATTRIBUTE_PREFERRED_LOCATION,
          (CUdeviceptr)m_data,
          m_nc * m_ni * m_nj * m_nk * sizeof(float_sw4)) != CUDA_SUCCESS) {
    std::cerr << " cuMemRangeGetAttributes failed\n";
  } else {
    if (data != 0) {
      std::cout << m_data << " PREF LOCATION IS DEVICE " << data << " in "
                << file << " line " << line << " CPU = " << CU_DEVICE_CPU
                << " INV = " << CU_DEVICE_INVALID << "\n"
                << std::flush;
      // abort();
    }
  }
#endif
}
//----------------------------------------------------------------------
void Sarray::switch_space(Space new_space) {
#ifdef ENABLE_GPU

  // std::cout<<"Switching from "<<as_int(space)<<" to
  // "<<as_int(new_space)<<"\n"<<std::flush;
  if (space == new_space) return;

  float_sw4* n_data = SW4_NEW(new_space, float_sw4[m_nc * m_ni * m_nj * m_nk]);

  spacecopy(n_data, m_data, new_space, space, m_nc * m_ni * m_nj * m_nk);

  ::operator delete[](m_data, space);
  reference(n_data);
  space = new_space;
  // std::cout<<"Switching from "<<as_int(space)<<" to "<<as_int(new_space)<<"
  // DONE\n"<<std::flush;
#endif
}
//----------------------------------------------------------------------
void mset_to_zero_async(Sarray& S0, Sarray& S1, Sarray& S2, Sarray& S3) {
#if defined(RAJA_ONLY)
  S0.set_to_zero_async();
  S1.set_to_zero_async();
  S2.set_to_zero_async();
  S3.set_to_zero_async();
#else
  float_sw4* m0 = S0.m_data;
  float_sw4* m1 = S1.m_data;
  float_sw4* m2 = S2.m_data;
  float_sw4* m3 = S3.m_data;

  size_t zero = 0;
  multiforall<512>(
      zero, S0.m_npts, [=] RAJA_DEVICE(size_t i) { m0[i] = 0; }, zero,
      S1.m_npts, [=] RAJA_DEVICE(size_t i) { m1[i] = 0; }, zero, S2.m_npts,
      [=] RAJA_DEVICE(size_t i) { m2[i] = 0; }, zero, S3.m_npts,
      [=] RAJA_DEVICE(size_t i) { m3[i] = 0; });
  // gmforall<512>(
  // 		  zero, S0.m_npts,[=] RAJA_DEVICE(size_t i) {
  // 		    m0[i] = 0;
  // 		  },
  // 		  zero, S1.m_npts,[=] RAJA_DEVICE(size_t i) {
  // 		    m1[i] = 0;
  // 		   },
  // 		  zero, S2.m_npts,[=] RAJA_DEVICE(size_t i) { m2[i] = 0;
  // 		   },
  // 		  zero, S3.m_npts,[=] RAJA_DEVICE(size_t i) { m3[i] = 0;
  // 		   });
#endif
}
int aligned(double* p) {
  for (int i = 16; i <= 2048; i *= 2)
    if ((long int)p % i != 0) return i / 2;
  return 2048;
}
void vset_to_zero_async(std::vector<Sarray>& v, int N) {
  // for(int i=0;i<N;i++)
  // std::cout<<"SIZES "<<v[i].m_npts<<"\n";

#if defined(RAJA_ONLY)
  for (int g = 0; g < N; g++) v[g].set_to_zero_async();

#else
  float_sw4 *m0, *m1, *m2, *m3, *m4, *m5, *m6, *m7;
  size_t zero = 0;

  switch (N) {
    case 1:
      m0 = v[0].m_data;

      gmforall<512>(zero, v[0].m_npts,
                    [=] RAJA_DEVICE(size_t i) { m0[i] = 0; });

      break;

    case 2:
      m0 = v[0].m_data;
      m1 = v[1].m_data;

      gmforall<512>(
          zero, v[0].m_npts, [=] RAJA_DEVICE(size_t i) { m0[i] = 0; }, zero,
          v[1].m_npts, [=] RAJA_DEVICE(size_t i) { m1[i] = 0; });
      break;

    case 3:
      m0 = v[0].m_data;
      m1 = v[1].m_data;
      m2 = v[2].m_data;
      gmforall<512>(
          zero, v[0].m_npts, [=] RAJA_DEVICE(size_t i) { m0[i] = 0; }, zero,
          v[1].m_npts, [=] RAJA_DEVICE(size_t i) { m1[i] = 0; }, zero,
          v[2].m_npts, [=] RAJA_DEVICE(size_t i) { m2[i] = 0; });

      break;

    case 4:

      m0 = v[0].m_data;
      m1 = v[1].m_data;
      m2 = v[2].m_data;
      m3 = v[3].m_data;
      gmforall<512>(
          zero, v[0].m_npts, [=] RAJA_DEVICE(size_t i) { m0[i] = 0; }, zero,
          v[1].m_npts, [=] RAJA_DEVICE(size_t i) { m1[i] = 0; }, zero,
          v[2].m_npts, [=] RAJA_DEVICE(size_t i) { m2[i] = 0; }, zero,
          v[3].m_npts, [=] RAJA_DEVICE(size_t i) { m3[i] = 0; });

      break;
    case 5:

      m0 = v[0].m_data;
      m1 = v[1].m_data;
      m2 = v[2].m_data;
      m3 = v[3].m_data;
      m4 = v[4].m_data;
      gmforall<512>(
          zero, v[0].m_npts, [=] RAJA_DEVICE(size_t i) { m0[i] = 0; }, zero,
          v[1].m_npts, [=] RAJA_DEVICE(size_t i) { m1[i] = 0; }, zero,
          v[2].m_npts, [=] RAJA_DEVICE(size_t i) { m2[i] = 0; }, zero,
          v[3].m_npts, [=] RAJA_DEVICE(size_t i) { m3[i] = 0; }, zero,
          v[4].m_npts, [=] RAJA_DEVICE(size_t i) { m4[i] = 0; });

      break;
    case 6:

      m0 = v[0].m_data;
      m1 = v[1].m_data;
      m2 = v[2].m_data;
      m3 = v[3].m_data;
      m4 = v[4].m_data;
      m5 = v[5].m_data;
      gmforall<512>(
          zero, v[0].m_npts, [=] RAJA_DEVICE(size_t i) { m0[i] = 0; }, zero,
          v[1].m_npts, [=] RAJA_DEVICE(size_t i) { m1[i] = 0; }, zero,
          v[2].m_npts, [=] RAJA_DEVICE(size_t i) { m2[i] = 0; }, zero,
          v[3].m_npts, [=] RAJA_DEVICE(size_t i) { m3[i] = 0; }, zero,
          v[4].m_npts, [=] RAJA_DEVICE(size_t i) { m4[i] = 0; }, zero,
          v[5].m_npts, [=] RAJA_DEVICE(size_t i) { m5[i] = 0; });
      break;
    case 7:

      m0 = v[0].m_data;
      m1 = v[1].m_data;
      m2 = v[2].m_data;
      m3 = v[3].m_data;
      m4 = v[4].m_data;
      m5 = v[5].m_data;
      m6 = v[6].m_data;
      gmforall<512>(
          zero, v[0].m_npts, [=] RAJA_DEVICE(size_t i) { m0[i] = 0; }, zero,
          v[1].m_npts, [=] RAJA_DEVICE(size_t i) { m1[i] = 0; }, zero,
          v[2].m_npts, [=] RAJA_DEVICE(size_t i) { m2[i] = 0; }, zero,
          v[3].m_npts, [=] RAJA_DEVICE(size_t i) { m3[i] = 0; }, zero,
          v[4].m_npts, [=] RAJA_DEVICE(size_t i) { m4[i] = 0; }, zero,
          v[5].m_npts, [=] RAJA_DEVICE(size_t i) { m5[i] = 0; }, zero,
          v[6].m_npts, [=] RAJA_DEVICE(size_t i) { m6[i] = 0; });
      break;

    case 8:

      m0 = v[0].m_data;
      m1 = v[1].m_data;
      m2 = v[2].m_data;
      m3 = v[3].m_data;
      m4 = v[4].m_data;
      m5 = v[5].m_data;
      m6 = v[6].m_data;
      m7 = v[7].m_data;
      gmforall<512>(
          zero, v[0].m_npts, [=] RAJA_DEVICE(size_t i) { m0[i] = 0; }, zero,
          v[1].m_npts, [=] RAJA_DEVICE(size_t i) { m1[i] = 0; }, zero,
          v[2].m_npts, [=] RAJA_DEVICE(size_t i) { m2[i] = 0; }, zero,
          v[3].m_npts, [=] RAJA_DEVICE(size_t i) { m3[i] = 0; }, zero,
          v[4].m_npts, [=] RAJA_DEVICE(size_t i) { m4[i] = 0; }, zero,
          v[5].m_npts, [=] RAJA_DEVICE(size_t i) { m5[i] = 0; }, zero,
          v[6].m_npts, [=] RAJA_DEVICE(size_t i) { m6[i] = 0; }, zero,
          v[7].m_npts, [=] RAJA_DEVICE(size_t i) { m7[i] = 0; });
      break;

    default:
      std::cerr << "ERROR:: vset_to_zero_async not implemented for " << N
                << " grids\n";
      abort();
  }

#endif
}
// Old norm based on atomics left as reference.
// Does not work well for comparisons and is slower 
// float_sw4 Sarray::norm() {
//   float_sw4* sum;
//   sum = SW4_NEW(Space::Managed_temps, float_sw4[1]);
//   sum[0] = 0.0;
//   float_sw4* lm_data = m_data;
//   RAJA::forall<DEFAULT_LOOP1>(
//       RAJA::RangeSegment(0, m_npts), [=] RAJA_DEVICE(size_t i) {
//         RAJA::atomicAdd<RAJA::auto_atomic>(sum, lm_data[i] * lm_data[i]);
//       });
//   float_sw4 retval = *sum;
//   ::operator delete[](sum, Space::Managed_temps);
//   return retval;
// }
float_sw4 Sarray::norm() {

  RAJA::ReduceSum<REDUCTION_POLICY,float_sw4> rsum(0);
  float_sw4* lm_data = m_data;
  RAJA::forall<DEFAULT_LOOP1>(
      RAJA::RangeSegment(0, m_npts), [=] RAJA_DEVICE(size_t i) {
	rsum+=lm_data[i] * lm_data[i];
      });
  float_sw4 retval = static_cast<float_sw4>(rsum.get());
  return retval;
}
 size_t Sarray::fwrite(FILE *file){
   size_t size = m_nc*m_ni*m_nj*m_nk;
   if (std::fwrite(m_data,sizeof(float_sw4),m_nc*m_ni*m_nj*m_nk,file)!=size){
     std::cerr<<"Write failed in Sarray::fwrite\n";
     abort();
   }
   return size;
 }
 size_t Sarray::fread(FILE *file){
   size_t size = m_nc*m_ni*m_nj*m_nk;
   if (std::fread(m_data,sizeof(float_sw4),m_nc*m_ni*m_nj*m_nk,file)!=size){
     std::cerr<<"Read failed in Sarray::fread\n";
     abort();
   }
   return size;
 }
