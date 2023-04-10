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
#ifndef EW_SARRAY_H
#define EW_SARRAY_H

#include <sys/types.h>
//#include <iostream>
//#include <mpi.h>
#include <cstdlib>
#include <string>
#include <vector>

#include "Mspace.h"
#include "Require.h"
#include "policies.h"
#include "sw4.h"

using std::string;

class EWCuda;
class Sarray;

class SView {
 public:
  float_sw4* data;
  ssize_t base;
  size_t offc, offi, offj, offk;
  SView(float_sw4* data, ssize_t base, size_t offc, size_t offi, size_t offj,
        size_t offk);
  // SView(SView &in){}
  SView(Sarray& x);
  SView();
  bool compare(SView& in) {
    if ((in.data != data) || (in.base != base) || (in.offi != offi) ||
        (in.offj != offj) || (in.offk != offk)) {
      return true;
    } else
      return false;
  }
  void set(Sarray& x);
  RAJA_HOST_DEVICE inline float_sw4& operator()(sw4_type c, sw4_type i, sw4_type j,
                                                sw4_type k) const {
    return data[base + c * offc + i * offi + j * offj +
                k * offk];  // When will this overflow if at all PBUGS
  }

  RAJA_HOST_DEVICE inline float_sw4& operator()(sw4_type i, sw4_type j, sw4_type k) const {
    return data[base + offc + i * offi + j * offj + k * offk];
  }
  RAJA_HOST_DEVICE void print(bool cond) const {
    if (cond)
      printf("SView pointer = %p base = %d offi = %d %d %d\n", data, int(base),
             int(offi), int(offj), int(offk));
  }
};

class Sarray {
 public:
  Space space;
  //   Sarray( CartesianProcessGrid* cartcomm, sw4_type nc=1 );
  Sarray(sw4_type nc, sw4_type ibeg, sw4_type iend, sw4_type jbeg, sw4_type jend, sw4_type kbeg, sw4_type kend,
         const char* file, sw4_type line);
  Sarray(sw4_type nc, sw4_type ibeg, sw4_type iend, sw4_type jbeg, sw4_type jend, sw4_type kbeg, sw4_type kend);
  Sarray(sw4_type ibeg, sw4_type iend, sw4_type jbeg, sw4_type jend, sw4_type kbeg, sw4_type kend);
  Sarray(sw4_type nc, sw4_type iend, sw4_type jend, sw4_type kend);
  Sarray(sw4_type iend, sw4_type jend, sw4_type kend);
  Sarray(const Sarray& u);
  Sarray(const Sarray& u, Space space);
  Sarray(Sarray& u, sw4_type nc = -1);
  Sarray();
  ~Sarray() {
#ifndef SW4_USE_UMPIRE
    if ((m_data != 0) && (!static_alloc)) ::operator delete[](m_data, space);
      // else std::cout<<"Skipped delete\n"<<std::flush;
#else
    if (m_data != 0) {
      if (static_alloc) {
        ::operator delete[](
            m_data,
            Space::Managed_temps);  // THIS NEEDS TO MATCH THE SPACE IN THE CTOR
                                    // umpire::ResourceManager &rma =
        // umpire::ResourceManager::getInstance(); auto allocator =
        // rma.getAllocator("UM_pool_small"); allocator.deallocate(m_data);
      } else
        ::operator delete[](m_data, space);
    }
#endif
  }
  //   void define( CartesianProcessGrid* cartcomm, sw4_type nc );
  inline Sarray& operator=(Sarray& rhs) {
    // std::cout<<"opertaot=\n"<<std::flush;
    if (this == &rhs) return *this;
    if ((this->space == Space::Device) || (rhs.space == Space::Device)) {
      std::cerr << "Sarray::operator= node implemented from device memory\n";
      abort();
    }
    for (sw4_type i = 0; i < m_npts; i++)
      this->m_data[i] = rhs.m_data[i];  // Assumes Space is host or Managed,
                                        // never device PBUGS
    return *this;
  }
  void define(sw4_type iend, sw4_type jend, sw4_type kend);
  void define(sw4_type nc, sw4_type iend, sw4_type jend, sw4_type kend);
  void define(sw4_type nc, sw4_type ibeg, sw4_type iend, sw4_type jbeg, sw4_type jend, sw4_type kbeg,
              sw4_type kend, Space space = Space::Managed);
  void define(sw4_type ibeg, sw4_type iend, sw4_type jbeg, sw4_type jend, sw4_type kbeg, sw4_type kend,
              Space space = Space::Managed);
  void define(const Sarray& u);
  inline float_sw4* c_ptr() { return m_data; }
#ifdef SW4_CUDA
  __host__ __device__ inline float_sw4* dev_ptr() { return dev_data; }
#else
  inline float_sw4* dev_ptr() { return dev_data; }
#endif
  void reference(float_sw4* new_data) {
    m_data = new_data;
    view.set(*this);
  }
  void reference_dev(float_sw4* new_data) { dev_data = new_data; }

  //   inline float_sw4& operator()( sw4_type c, sw4_type i, sw4_type j, sw4_type k )
  //   {return
  //   m_data[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}
  inline bool in_range(sw4_type c, sw4_type i, sw4_type j, sw4_type k) const {
    return 1 <= c && c <= m_nc && m_ib <= i && i <= m_ie && m_jb <= j &&
           j <= m_je && m_kb <= k && k <= m_ke;
  }
  inline float_sw4& operator()(sw4_type c, sw4_type i, sw4_type j, sw4_type k) const {
#ifdef BZ_DEBUG
    VERIFY2(in_range(c, i, j, k), "Error Index (c,i,j,k) = ("
                                      << c << "," << i << "," << j << "," << k
                                      << ") not in range 1<= c <= " << m_nc
                                      << " " << m_ib << " <= i <= " << m_ie
                                      << " " << m_jb << " <= j <= " << m_je
                                      << " " << m_kb << " <=  k <= " << m_ke);
#endif
    //      return
    //      m_data[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}
    return m_data[m_base + m_offc * c + m_offi * i + m_offj * j + m_offk * k];
  }
  inline float_sw4& operator()(sw4_type i, sw4_type j, sw4_type k) const {
#ifdef BZ_DEBUG
    if (!in_range(1, i, j, k))
      VERIFY2(0, "Error Index (c,i,j,k) = ("
                     << 1 << "," << i << "," << j << "," << k
                     << ") not in range 1<= c <= " << m_nc << " " << m_ib
                     << " <= i <= " << m_ie << " " << m_jb << " <= j <= "
                     << m_je << " " << m_kb << " <= k <= " << m_ke);

#endif
    //      return
    //      m_data[m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}
    return m_data[m_base + m_offi * i + m_offj * j + m_offk * k + m_offc];
  }
  inline bool is_defined() { return m_data != NULL; }
  sw4_type m_ib, m_ie, m_jb, m_je, m_kb, m_ke;
  static bool m_corder;
  ssize_t m_base;
  size_t m_offi, m_offj, m_offk, m_offc, m_npts;
  //   sw4_type index( sw4_type i, sw4_type j, sw4_type k ) {return
  //   (i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb);}
  size_t index(sw4_type i, sw4_type j, sw4_type k) {
    return m_base + m_offc + m_offi * i + m_offj * j + m_offk * k;
  }
#ifdef SW4_CUDA
  __host__ __device__ size_t index(sw4_type c, sw4_type i, sw4_type j, sw4_type k) {
    return m_base + m_offc * c + m_offi * i + m_offj * j + m_offk * k;
  }
#else
  size_t index(sw4_type c, sw4_type i, sw4_type j, sw4_type k) {
    return m_base + m_offc * c + m_offi * i + m_offj * j + m_offk * k;
  }
#endif
  void getnonzero() {
    for (sw4_type i = 0; i < m_nc * m_ni * m_nj * m_nk; i++)
      if (m_data[i] > 0.0)
        std::cout << "FORCE " << i << " " << m_data[i] << "\n";
  }
  void intersection(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                    sw4_type wind[6]);
  void side_plane(sw4_type side, sw4_type wind[6], sw4_type nGhost = 1);
  void side_plane_fortran(sw4_type side, sw4_type wind[6], sw4_type nGhost = 1);
  bool in_domain(sw4_type i, sw4_type j, sw4_type k);
  void set_to_zero();
  void set_to_zero_async();
  void set_to_minusOne();
  void set_to_minusOneHost();
  void set_value(float_sw4 scalar);
  void set_valueHost(float_sw4 scalar);
  void set_value_async(float_sw4 scalar);
  void set_to_random(float_sw4 llim = 0.0, float_sw4 ulim = 1.0);
  void save_to_disk(const char* fname);
  sw4_type ncomp() const { return m_nc; }
  sw4_type npts() const { return m_ni * m_nj * m_nk; }
  void copy(const Sarray& u);
  float_sw4 maximum(sw4_type c = 1);
  float_sw4 minimum(sw4_type c = 1);
  float_sw4 sum(sw4_type c = 1);
  size_t count_nans();
  size_t count_nans(sw4_type& cfirst, sw4_type& ifirst, sw4_type& jfirst, sw4_type& kfirst);
  size_t check_match_cpu_gpu(EWCuda* cu, string name);
  size_t check_match_cpu_gpu(EWCuda* cu, sw4_type& cfirst, sw4_type& ifirst, sw4_type& jfirst,
                             sw4_type& kfirst, string name);
  void extract_subarray(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                        float_sw4* ar);
  void extract_subarrayIK(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                          float_sw4* ar);

  void insert_subarray(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                       double* ar);
  void insert_subarray(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                       float* ar);
  void insert_intersection(Sarray& a_U);
  void insert_subarrayIK(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                         float_sw4* ar);
  void copy_kplane(Sarray& u, sw4_type k);
  void copy_kplane2(Sarray& u, sw4_type k);
  void assign(const float* ar, sw4_type corder);
  void assign(const double* ar, sw4_type corder);
  void extract(double* ar, sw4_type corder);
  void assign(const float* ar);
  void assign(const double* ar);
  void transposeik();
  void extrapolij(sw4_type npts);
  void copy_to_device(EWCuda* cu, bool async = false, sw4_type st = 0);
  void copy_from_device(EWCuda* cu, bool async = false, sw4_type st = 0);
  void allocate_on_device(EWCuda* cu);
  void page_lock(EWCuda* cu);
  void page_unlock(EWCuda* cu);
  void swrite(std::string filename);
  size_t fwrite(FILE *file);
  size_t fread(FILE *file);
  Sarray* create_copy_on_device(EWCuda* cu);
  void define_offsets();
  void GetAtt(char* file, sw4_type line);
  //   void write( char* filename, CartesianProcessGrid* cartcomm,
  //   std::vector<float_sw4> pars );
  sw4_type m_nc, m_ni, m_nj, m_nk;
  void prefetch(sw4_type device = 0);
  void forceprefetch(sw4_type device = 0);
  void switch_space(Space space_in);
  float_sw4 norm();
  inline SView& getview() {
    // prefetch();
    return view;
  }
  SView view;
  static std::unordered_map<std::string, float_sw4*> static_map;
  bool static_alloc;

 private:
  float_sw4* m_data;
  bool prefetched;
  friend void mset_to_zero_async(Sarray& A, Sarray& B, Sarray& C, Sarray& D);
  friend void vset_to_zero_async(std::vector<Sarray>& v, sw4_type N);
  std::ofstream of;
  float_sw4* dev_data;
  inline sw4_type min(sw4_type i1, sw4_type i2) {
    if (i1 < i2)
      return i1;
    else
      return i2;
  }
  inline sw4_type max(sw4_type i1, sw4_type i2) {
    if (i1 > i2)
      return i1;
    else
      return i2;
  }
  //   void init_mpi_datatype( CartesianProcessGrid* cartcomm );
  //    bool m_mpi_datatype_initialized;
  //    MPI_Datatype m_local_block_type;
};
void SarrayVectorPrefetch(std::vector<Sarray>& v);
void SarrayVectorPrefetch(std::vector<Sarray*>& v);
void SarrayVectorPrefetch(std::vector<Sarray*>& v, sw4_type n);

float_sw4* memoize(Sarray& u, sw4_type c, sw4_type i, sw4_type j, sw4_type k);

void mset_to_zero_async(Sarray& S0, Sarray& S1, Sarray& S2, Sarray& S3);
#endif
