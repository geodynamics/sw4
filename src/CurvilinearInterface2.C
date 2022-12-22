#include <iostream>

#include "CurvilinearInterface2.h"
#include "EW.h"
#include "F77_FUNC.h"
#include "GridGenerator.h"
#include "Sarray.h"
#include "TestEcons.h"
#include "TestTwilight.h"
#include "caliper.h"
#include "foralls.h"
#include "policies.h"
extern "C" {
void F77_FUNC(dgetrf, DGETRF)(int*, int*, double*, int*, int*, int*);
void F77_FUNC(dgetrs, DGETRS)(char*, int*, int*, double*, int*, int*, double*,
                              int*, int*);
}

void bndryOpNoGhostc(double* acof_no_gp, double* ghcof_no_gp,
                     double* sbop_no_gp);

void curvilinear4sgwind(int, int, int, int, int, int, int, int, float_sw4*,
                        float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                        float_sw4*, int*, float_sw4*, float_sw4*, float_sw4*,
                        float_sw4*, float_sw4*, float_sw4*, float_sw4*, int,
                        char);

//-----------------------------------------------------------------------
CurvilinearInterface2::CurvilinearInterface2(int a_gc, EW* a_ew) {
  SW4_MARK_FUNCTION;
#if defined(ENABLE_CUDA)
#if defined(SW4_MANAGED_MPI_BUFFERS)
  if (!mpi_supports_device_buffers()) {
    std::cerr << "**** ERROR::   SW4_MANAGED_MPI_BUFFERS selected at compile "
                 "time. Cases must be "
                 "run with the -M -gpu flag****************\n";
    abort();
  }
#endif
#endif

  m_reltol = 1e-6;
  m_abstol = 1e-6;
  m_maxit = 30;
  m_gc = a_gc;
  m_gf = a_gc + 1;
  m_ew = a_ew;
  m_etest = a_ew->create_energytest();
  m_tw = a_ew->create_twilight();
  m_psource = a_ew->get_point_source_test();
  m_nghost = 5;

#if defined(ENABLE_GPU)
  float_sw4* tmpa =
      SW4_NEW(Space::Managed, float_sw4[6 + 384 + 24 + 48 + 6 + 384 + 6 + 6]);
  m_sbop = tmpa;  // PTR_PUSH(Space::Managed,m_sbop);
  m_acof = m_sbop + 6;
  PTR_PUSH(Space::Managed, m_acof, 384 * sizeof(float_sw4));
  m_bop = m_acof + 384;
  PTR_PUSH(Space::Managed, m_bop, 24 * sizeof(float_sw4));
  m_bope = m_bop + 24;
  PTR_PUSH(Space::Managed, m_bope, 48 * sizeof(float_sw4));
  m_ghcof = m_bope + 48;
  PTR_PUSH(Space::Managed, m_ghcof, 6 * sizeof(float_sw4));

  m_acof_no_gp = m_ghcof + 6;
  PTR_PUSH(Space::Managed, m_acof_no_gp, 384 * sizeof(float_sw4));
  m_ghcof_no_gp = m_acof_no_gp + 384;
  PTR_PUSH(Space::Managed, m_ghcof_no_gp, 6 * sizeof(float_sw4));
  m_sbop_no_gp = m_ghcof_no_gp + 6;
  PTR_PUSH(Space::Managed, m_sbop_no_gp, 6 * sizeof(float_sw4));
#endif

  a_ew->GetStencilCoefficients(m_acof, m_ghcof, m_bop, m_bope, m_sbop);
  bndryOpNoGhostc(m_acof_no_gp, m_ghcof_no_gp, m_sbop_no_gp);
  for (int s = 0; s < 4; s++) m_isbndry[s] = true;
  m_use_attenuation = a_ew->usingAttenuation();
  m_number_mechanisms = a_ew->getNumberOfMechanisms();
  // Allocate space for MPI buffers
}

CurvilinearInterface2::~CurvilinearInterface2() {
  std::cout << "~CurvilinearInterface2()...\n" << std::flush;
#if defined(ENABLE_CUDA)
  ::operator delete[](m_sbop, Space::Managed);
#endif
  ::operator delete[](m_strx_c, Space::Managed);
  ::operator delete[](m_stry_c, Space::Managed);
  ::operator delete[](m_strx_f, Space::Managed);
  ::operator delete[](m_stry_f, Space::Managed);
  ::operator delete[](m_mpi_buffer_space, Space::Pinned);
#ifdef USE_MAGMA
  ::operator delete[](m_mass_block, Space::Managed);
  ::operator delete[](dA_array, Space::Managed);
  ::operator delete[](m_ipiv_block, Space::Managed);
  ::operator delete[](piv_array, Space::Managed);
  //::operator delete[](info,Space::Managed);
  ::operator delete[](dB_array, Space::Managed);
  ::operator delete[](x, Space::Managed);
#endif
#ifdef USE_DIRECT_INVERSE
  ::operator delete[](m_mass_block, Space::Managed);
  ::operator delete[](x, Space::Managed);
#endif
  std::cout << "~CurvilinearInterface2().. Done\n" << std::flush;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::bnd_zero_host(Sarray& u, int npts) {
  // Homogeneous Dirichet at boundaries on sides. Do not apply at upper and
  // lower boundaries.
  for (int s = 0; s < 4; s++)
    if (m_isbndry[s]) {
      int kb = u.m_kb, ke = u.m_ke, jb = u.m_jb, je = u.m_je, ib = u.m_ib,
          ie = u.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;
      for (int c = 1; c <= u.m_nc; c++)
        for (int k = kb; k <= ke; k++)
          for (int j = jb; j <= je; j++)
            for (int i = ib; i <= ie; i++) u(c, i, j, k) = 0;
    }
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::bnd_zero(Sarray& u, int npts) {
  SW4_MARK_FUNCTION;
  // Homogeneous Dirichet at boundaries on sides. Do not apply at upper and
  // lower boundaries.

  int nc = u.m_nc;
  SView& uV = u.getview();

#if defined(RAJA_ONLY) || defined(__cray__) || not defined(ENABLE_GPU)
  for (int s = 0; s < 4; s++)
    if (m_isbndry[s]) {
      int kb = u.m_kb, ke = u.m_ke, jb = u.m_jb, je = u.m_je, ib = u.m_ib,
          ie = u.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;
      // for (int c = 1; c <= u.m_nc; c++)
      // for (int k = kb; k <= ke; k++)
      //  for (int j = jb; j <= je; j++)
      //     for (int i = ib; i <= ie; i++) u(c,i,j,k)=0;
      RAJA::RangeSegment k_range(kb, ke + 1);
      RAJA::RangeSegment j_range(jb, je + 1);
      RAJA::RangeSegment i_range(ib, ie + 1);
      RAJA::kernel<BZ_POL_ASYNC>(RAJA::make_tuple(k_range, j_range, i_range),
                                 [=] RAJA_DEVICE(int k, int j, int i) {
                                   for (int c = 1; c <= nc; c++)
                                     uV(c, i, j, k) = 0;
                                 });
    }
    // SYNC_STREAM;
#else
  class MRange start[4], end[4];
  std::vector<MRange> st, en;
  for (int s = 0; s < 4; s++)
    if (m_isbndry[s]) {
      int kb = u.m_kb, ke = u.m_ke, jb = u.m_jb, je = u.m_je, ib = u.m_ib,
          ie = u.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;

      start[s] = {ib, jb, kb};
      st.push_back(start[s]);
      end[s] = {ie + 1, je + 1, ke + 1};
      en.push_back(end[s]);
    }
  if (st.size() == 0) {
  } else if (st.size() == 1) {
    gmforall3async<16, 16, 1>(st[0], en[0],
                              [=] RAJA_DEVICE(int i, int j, int k) {
                                for (int c = 1; c <= nc; c++)
                                  uV(c, i, j, k) = 0;
                              });
  } else if (st.size() == 2) {
    gmforall3async<16, 16, 1>(
        st[0], en[0],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        },
        st[1], en[1],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        });
  } else if (st.size() == 3) {  // NOT TESTED PBUGS
    gmforall3async<16, 16, 1>(
        st[0], en[0],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        },
        st[1], en[1],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        },
        st[2], en[2],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        });
  } else if (st.size() == 4) {  // NOT TESTED PBUGS
    gmforall3async<16, 16, 1>(
        st[0], en[0],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        },
        st[1], en[1],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        },
        st[2], en[2],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        },
        st[3], en[3],
        [=] RAJA_DEVICE(int i, int j, int k) {
          for (int c = 1; c <= nc; c++) uV(c, i, j, k) = 0;
        });

  } else {
    std::cerr << " ERROR SIZE in CurvilinearInterface2::bnd_zero 0 "
              << st.size() << "\n";
    abort();
  }

#endif
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::copy_str(float_sw4* dest, float_sw4* src,
                                     int offset, int n, int nsw) {
  SW4_MARK_FUNCTION;
  //
  // Copy supergrid stretching function array into an array with different
  // number of ghost points. The new values are filled in by constant
  // extrapolation. Input:  src    - Old array
  //         offset - Number of added points points in new array at the
  //                  lower end (<0 means fewer points)
  //         n      - Size of new array
  //         nsw    - Size of old array
  //
  // Output: dest   - New array
  // Note: Giving n as input implicitly determines the number of new points at
  // the upper end.
  //
  if (offset >= 0) {
    for (int i = 0; i < nsw; i++) dest[i + offset] = src[i];
    for (int i = 0; i < offset; i++) dest[i] = dest[offset];
    for (int i = nsw; i < n; i++) dest[i] = dest[nsw - 1];
  } else {
    for (int i = 0; i < min(n, nsw - offset); i++) dest[i] = src[i + offset];
  }
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::init_arrays(vector<float_sw4*>& a_strx,
                                        vector<float_sw4*>& a_stry) {
  SW4_MARK_FUNCTION;
  // std::cout << "void CurvilinearInterface2::init_arrays \n";
  for (int s = 0; s < 4; s++)
    m_isbndry[s] = m_ew->getLocalBcType(m_gc, s) != bProcessor;

  m_ib = m_ew->m_iStartInt[m_gc] - m_nghost;
  m_ie = m_ew->m_iEndInt[m_gc] + m_nghost;
  m_jb = m_ew->m_jStartInt[m_gc] - m_nghost;
  m_je = m_ew->m_jEndInt[m_gc] + m_nghost;
  m_ibf = m_ew->m_iStartInt[m_gf] - m_nghost;
  m_ief = m_ew->m_iEndInt[m_gf] + m_nghost;
  m_jbf = m_ew->m_jStartInt[m_gf] - m_nghost;
  m_jef = m_ew->m_jEndInt[m_gf] + m_nghost;
  m_nkf = m_ew->m_global_nz[m_gf];

  m_kb = 0;
  m_ke = 8;
  m_kbf = m_nkf - 7;
  m_kef = m_nkf + 1;

  m_strx_c = SW4_NEW(Space::Managed, float_sw4[m_ie - m_ib + 1]);
  m_stry_c = SW4_NEW(Space::Managed, float_sw4[m_je - m_jb + 1]);
  m_strx_f = SW4_NEW(Space::Managed, float_sw4[m_ief - m_ibf + 1]);
  m_stry_f = SW4_NEW(Space::Managed, float_sw4[m_jef - m_jbf + 1]);
#ifdef ENABLE_GPU
  float_sw4* lm_strx_c = new float_sw4[m_ie - m_ib + 1];
  float_sw4* lm_stry_c = new float_sw4[m_je - m_jb + 1];
  float_sw4* lm_strx_f = new float_sw4[m_ief - m_ibf + 1];
  float_sw4* lm_stry_f = new float_sw4[m_jef - m_jbf + 1];
#else

  float_sw4* lm_strx_c = m_strx_c;
  float_sw4* lm_stry_c = m_stry_c;
  float_sw4* lm_strx_f = m_strx_f;
  float_sw4* lm_stry_f = m_stry_f;
  
#endif

  allocate_mpi_buffers();

  int ndif = m_nghost - (m_ew->m_iStartInt[m_gc] - m_ew->m_iStart[m_gc]);
  int nsw = m_ew->m_iEnd[m_gc] - m_ew->m_iStart[m_gc] + 1;
  copy_str(lm_strx_c, a_strx[m_gc], ndif, m_ie - m_ib + 1, nsw);
  communicate_array1d(lm_strx_c, m_ie - m_ib + 1, 0, m_nghost);
#ifdef ENABLE_CUDA
  cudaMemcpyAsync(m_strx_c, lm_strx_c, (m_ie - m_ib + 1) * sizeof(double),
                  cudaMemcpyHostToDevice, 0);
#endif
#ifdef ENABLE_HIP
  hipMemcpyAsync(m_strx_c, lm_strx_c, (m_ie - m_ib + 1) * sizeof(double),
                 hipMemcpyHostToDevice, 0);
#endif

  ndif = m_nghost - (m_ew->m_iStartInt[m_gf] - m_ew->m_iStart[m_gf]);
  nsw = m_ew->m_iEnd[m_gf] - m_ew->m_iStart[m_gf] + 1;
  copy_str(lm_strx_f, a_strx[m_gf], ndif, m_ief - m_ibf + 1, nsw);
  communicate_array1d(lm_strx_f, m_ief - m_ibf + 1, 0, m_nghost);
#ifdef ENABLE_CUDA
  cudaMemcpyAsync(m_strx_f, lm_strx_f, (m_ief - m_ibf + 1) * sizeof(double),
                  cudaMemcpyHostToDevice, 0);
#endif
#ifdef ENABLE_HIP
  hipMemcpyAsync(m_strx_f, lm_strx_f, (m_ief - m_ibf + 1) * sizeof(double),
                 hipMemcpyHostToDevice, 0);
#endif

  ndif = m_nghost - (m_ew->m_jStartInt[m_gc] - m_ew->m_jStart[m_gc]);
  nsw = m_ew->m_jEnd[m_gc] - m_ew->m_jStart[m_gc] + 1;
  copy_str(lm_stry_c, a_stry[m_gc], ndif, m_je - m_jb + 1, nsw);
  communicate_array1d(lm_stry_c, m_je - m_jb + 1, 1, m_nghost);
#ifdef ENABLE_CUDA
  cudaMemcpyAsync(m_stry_c, lm_stry_c, (m_je - m_jb + 1) * sizeof(double),
                  cudaMemcpyHostToDevice, 0);
#endif
#ifdef ENABLE_HIP
  hipMemcpyAsync(m_stry_c, lm_stry_c, (m_je - m_jb + 1) * sizeof(double),
                 hipMemcpyHostToDevice, 0);
#endif

  ndif = m_nghost - (m_ew->m_jStartInt[m_gf] - m_ew->m_jStart[m_gf]);
  nsw = m_ew->m_jEnd[m_gf] - m_ew->m_jStart[m_gf] + 1;
  copy_str(lm_stry_f, a_stry[m_gf], ndif, m_jef - m_jbf + 1, nsw);
  communicate_array1d(lm_stry_f, m_jef - m_jbf + 1, 1, m_nghost);
#ifdef ENABLE_CUDA
  cudaMemcpyAsync(m_stry_f, lm_stry_f, (m_jef - m_jbf + 1) * sizeof(double),
                  cudaMemcpyHostToDevice, 0);
#endif
#ifdef ENABLE_HIP
  hipMemcpyAsync(m_stry_f, lm_stry_f, (m_jef - m_jbf + 1) * sizeof(double),
                 hipMemcpyHostToDevice, 0);
#endif

  SYNC_STREAM;

  m_rho_c.define(m_ib, m_ie, m_jb, m_je, 1, 1);
  m_rho_f.define(m_ibf, m_ief, m_jbf, m_jef, m_nkf, m_nkf);

  m_mu_c.define(m_ib, m_ie, m_jb, m_je, m_kb, m_ke);
  m_lambda_c.define(m_ib, m_ie, m_jb, m_je, m_kb, m_ke);
  m_jac_c.define(m_ib, m_ie, m_jb, m_je, m_kb, m_ke);

  m_mu_f.define(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);
  m_lambda_f.define(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);
  m_jac_f.define(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);

  m_x_c.define(m_ib, m_ie, m_jb, m_je, m_kb, m_ke);
  m_y_c.define(m_ib, m_ie, m_jb, m_je, m_kb, m_ke);
  m_z_c.define(m_ib, m_ie, m_jb, m_je, m_kb, m_ke);
  m_met_c.define(4, m_ib, m_ie, m_jb, m_je, m_kb, m_ke);
  m_ew->m_gridGenerator->generate_grid_and_met(m_ew, m_gc, m_x_c, m_y_c, m_z_c,
                                               m_jac_c, m_met_c, false);
  // std::cout<<"HERE 1\n";
  m_met_c.insert_intersection(m_ew->mMetric[m_gc]);
  // std::cout<<"HERE 1.1\n"<<std::flush;
  m_jac_c.insert_intersection(m_ew->mJ[m_gc]);
  // std::cout<<"HERE 1.2\n"<<std::flush;
  communicate_array(m_met_c, true);
  communicate_array(m_jac_c, true);

  m_x_f.define(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);
  m_y_f.define(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);
  m_z_f.define(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);
  m_met_f.define(4, m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);
  m_ew->m_gridGenerator->generate_grid_and_met(m_ew, m_gf, m_x_f, m_y_f, m_z_f,
                                               m_jac_f, m_met_f, false);
  m_met_f.insert_intersection(m_ew->mMetric[m_gf]);
  m_jac_f.insert_intersection(m_ew->mJ[m_gf]);

  communicate_array(m_met_f, true);
  communicate_array(m_jac_f, true);

  if (m_tw != 0) {
    m_tw->get_rho(m_rho_c, m_x_c, m_y_c, m_z_c);
    m_tw->get_rho(m_rho_f, m_x_f, m_y_f, m_z_f);
    m_tw->get_mula(m_mu_c, m_lambda_c, m_x_c, m_y_c, m_z_c);
    m_tw->get_mula(m_mu_f, m_lambda_f, m_x_f, m_y_f, m_z_f);
  } else {
    m_rho_c.insert_intersection(m_ew->mRho[m_gc]);
    m_rho_f.insert_intersection(m_ew->mRho[m_gf]);
    m_mu_c.insert_intersection(m_ew->mMu[m_gc]);
    m_mu_f.insert_intersection(m_ew->mMu[m_gf]);
    m_lambda_c.insert_intersection(m_ew->mLambda[m_gc]);
    m_lambda_f.insert_intersection(m_ew->mLambda[m_gf]);
    SYNC_STREAM;
    int extra_ghost = m_nghost - m_ew->getNumberOfGhostPoints();
    if (extra_ghost > 0) {
      if (m_etest != 0) {
        int sides[6] = {1, 1, 1, 1, 0, 0};
        m_etest->get_rhobnd(m_rho_c, extra_ghost, sides);
        m_etest->get_rhobnd(m_rho_f, extra_ghost, sides);
        m_etest->get_mulabnd(m_mu_c, m_lambda_c, extra_ghost, sides);
        m_etest->get_mulabnd(m_mu_f, m_lambda_f, extra_ghost, sides);
      } else {
        m_rho_c.extrapolij(extra_ghost);
        m_rho_f.extrapolij(extra_ghost);
        m_mu_c.extrapolij(extra_ghost);
        m_mu_f.extrapolij(extra_ghost);
        m_lambda_c.extrapolij(extra_ghost);
        m_lambda_f.extrapolij(extra_ghost);
      }
    }
  }
  communicate_array(m_rho_c, true);
  communicate_array(m_mu_c, true);
  communicate_array(m_lambda_c, true);

  communicate_array(m_rho_f, true);
  communicate_array(m_mu_f, true);
  communicate_array(m_lambda_f, true);

  if (m_use_attenuation) init_arrays_att();

  // Matrix only defined at interior points
  m_Mass_block.define(9, m_ib + m_nghost, m_ie - m_nghost, m_jb + m_nghost,
                      m_je - m_nghost, 1, 1);
  interface_block(m_Mass_block);

  // Repackage Mass_block into array of fortran order.
  int nimb = (m_Mass_block.m_ie - m_Mass_block.m_ib + 1);
  size_t msize = nimb * (m_Mass_block.m_je - m_Mass_block.m_jb + 1);
  // std::cout << "Batch size in setup is " << msize << "\n";
  //#define USE_DIRECT_INVERSE 1
#ifdef USE_DIRECT_INVERSE
  // std::cout << " USING DIRECT INVERSE \n";
  m_mass_block = SW4_NEW(Space::Managed, float_sw4[9 * msize]);
  x = SW4_NEW(Space::Managed, float_sw4[3 * msize]);
  for (int j = m_jb + m_nghost; j <= m_je - m_nghost; j++)
    for (int i = m_ib + m_nghost; i <= m_ie - m_nghost; i++) {
      size_t ind = (i - (m_ib + m_nghost)) + nimb * (j - (m_jb + m_nghost));
      for (int c = 1; c <= 9; c++)
        m_mass_block[c - 1 + 9 * ind] = m_Mass_block(c, i, j, 1);
    }
  invert(m_mass_block, msize);

#endif

#ifdef USE_MAGMA
  std::cout << " USING MAGMA FOR DGETRF WITH BATCH SIZE " << msize << "\n";
  m_mass_block = SW4_NEW(Space::Managed, float_sw4[9 * msize]);
  dA_array = SW4_NEW(Space::Managed, float_sw4* [msize]);
  m_ipiv_block = SW4_NEW(Space::Managed, int[3 * msize]);
  piv_array = SW4_NEW(Space::Managed, int* [msize]);
  magma_int_t* info = SW4_NEW(Space::Managed, int[msize]);
  // For use in solve
  dB_array = SW4_NEW(Space::Managed, float_sw4* [msize]);
  x = SW4_NEW(Space::Managed, float_sw4[3 * msize]);
  // End for use in solver
  for (int j = m_jb + m_nghost; j <= m_je - m_nghost; j++)
    for (int i = m_ib + m_nghost; i <= m_ie - m_nghost; i++) {
      size_t ind = (i - (m_ib + m_nghost)) + nimb * (j - (m_jb + m_nghost));
      for (int c = 1; c <= 9; c++)
        m_mass_block[c - 1 + 9 * ind] = m_Mass_block(c, i, j, 1);
    }
  for (int i = 0; i < msize; i++) {
    dA_array[i] = &m_mass_block[i * 9];
    piv_array[i] = &m_ipiv_block[3 * i];
    info[i] = 0;
    dB_array[i] = &x[3 * i];
  }

  magma_queue_create(0, &queue);

  magma_int_t retval;
  magma_int_t m = 3;
  SW4_MARK_BEGIN("MAGMA_DGETRF");
  retval =
      magma_dgetrf_batched(m, m, dA_array, m, piv_array, info, msize, queue);
  if (retval != MAGMA_SUCCESS) {
    std::cerr << "ERROR MAGMA DGETRF BATCHED call failed" << std::flush;
    abort();
  }
  SW4_MARK_END("MAGMA_DGETRF");
  int nfails = 0;
  for (int i = 0; i < msize; i++)
    if (info[i] != 0) {
      std::cerr << "ERROR :: MAGMA DGETRF FAILED for i =" << i << " " << info[i]
                << "\n";
      nfails++;
    }
  if (nfails != 0) {
    std::cerr << "ERRROR Stopping due to failed LU decomposition\n";
    abort();
  }
  ::operator delete[](info, Space::Managed);
  // Calculate the sub-batch sizes
  const int maxbatchsize = 65535;
  for (int i = 0; i < msize; i += maxbatchsize) {
    subbatchsize.push_back(maxbatchsize);
    subbatchoffset.push_back(i);
  }
  // The tail
  subbatchsize.back() = msize % (subbatchoffset.back());
  int sum = 0;
  for (int i = 0; i < subbatchsize.size(); i++) sum += subbatchsize[i];
  if (sum != msize) {
    std::cerr << "Magma sub-batchsize calc is incorrect " << sum
              << "!=" << msize << "\n";
    abort();
  }
  for (int i = 0; i < subbatchsize.size(); i++) {
    std::cout << "MAGMA Sub-batching " << i << " " << subbatchoffset[i] << " "
              << subbatchsize[i] << "\n";
  }
#endif
#ifdef USE_LAPACK_ON_CPU
  m_mass_block = new float_sw4[9 * msize];
  for (int j = m_jb + m_nghost; j <= m_je - m_nghost; j++)
    for (int i = m_ib + m_nghost; i <= m_ie - m_nghost; i++) {
      size_t ind = (i - (m_ib + m_nghost)) + nimb * (j - (m_jb + m_nghost));
      for (int c = 1; c <= 9; c++)
        m_mass_block[c - 1 + 9 * ind] = m_Mass_block(c, i, j, 1);
    }
  int three = 3;
  int info = 0;
  m_ipiv_block = new int[3 * msize];
  SW4_MARK_BEGIN("DGETRF");
  for (size_t ind = 0; ind < msize; ind++) {
    F77_FUNC(dgetrf, DGETRF)
    (&three, &three, &m_mass_block[9 * ind], &three, &m_ipiv_block[3 * ind],
     &info);
    if (info != 0) {
      int j = ind / m_Mass_block.m_ni + m_Mass_block.m_jb;
      int i =
          ind + m_Mass_block.m_ib - m_Mass_block.m_ni * (j - m_Mass_block.m_jb);
      std::cerr << "LU Fails at (i,j) equals" << i << "," << j
                << " info = " << info << " "
                << m_Mass_block(info + 3 * (info - 1), i, j, 1) << "\n";
      for (int l = 1; l <= 3; l++)
        for (int m = 1; m <= 3; m++)
          std::cerr << m_Mass_block(m + 3 * (l - 1), i, j, 1) << ",";
      std::cerr << "\n";
    }
  }
  SW4_MARK_END("DGETRF");
#endif
#ifdef ENABLE_GPU
  delete [] lm_strx_c;
  delete [] lm_stry_c;
  delete [] lm_strx_f;
  delete [] lm_stry_f;
#endif
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::init_arrays_att() {
  SW4_MARK_FUNCTION;
  // Attenuation material setup
  if (m_use_attenuation) {
    m_muve_c.resize(m_number_mechanisms);
    m_lambdave_c.resize(m_number_mechanisms);
    for (int a = 0; a < m_number_mechanisms; a++) {
      m_muve_c[a].define(m_ib, m_ie, m_jb, m_je, m_kb, m_ke);
      m_lambdave_c[a].define(m_ib, m_ie, m_jb, m_je, m_kb, m_ke);
    }
    m_muve_f.resize(m_number_mechanisms);
    m_lambdave_f.resize(m_number_mechanisms);
    for (int a = 0; a < m_number_mechanisms; a++) {
      m_muve_f[a].define(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);
      m_lambdave_f[a].define(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef);
    }
    if (m_tw != 0) {
      // Note, twilight uses only one attenuation mechanism
      m_tw->get_mula_att(m_muve_c[0], m_lambdave_c[0], m_x_c, m_y_c, m_z_c);
      m_tw->get_mula_att(m_muve_f[0], m_lambdave_f[0], m_x_f, m_y_f, m_z_f);
    } else {
      for (int a = 0; a < m_number_mechanisms; a++) {
        m_muve_c[a].insert_intersection(m_ew->mMuVE[m_gc][a]);
        m_lambdave_c[a].insert_intersection(m_ew->mLambdaVE[m_gc][a]);
        m_muve_f[a].insert_intersection(m_ew->mMuVE[m_gf][a]);
        m_lambdave_f[a].insert_intersection(m_ew->mLambdaVE[m_gf][a]);
      }
      SYNC_STREAM;
      int extra_ghost = m_nghost - m_ew->getNumberOfGhostPoints();
      if (extra_ghost > 0) {
        if (m_etest != 0) {
          int sides[6] = {1, 1, 1, 1, 0, 0};
          for (int a = 0; a < m_number_mechanisms; a++) {
            m_etest->get_mulabnd(m_muve_c[a], m_lambdave_c[a], extra_ghost,
                                 sides);
            m_etest->get_mulabnd(m_muve_f[a], m_lambdave_f[a], extra_ghost,
                                 sides);
          }
        } else {
          for (int a = 0; a < m_number_mechanisms; a++) {
            m_muve_c[a].extrapolij(extra_ghost);
            m_muve_f[a].extrapolij(extra_ghost);
            m_lambdave_c[a].extrapolij(extra_ghost);
            m_lambdave_f[a].extrapolij(extra_ghost);
          }
        }
      }
    }
    for (int a = 0; a < m_number_mechanisms; a++) {
      communicate_array(m_muve_c[a], true);
      communicate_array(m_lambdave_c[a], true);
      communicate_array(m_muve_f[a], true);
      communicate_array(m_lambdave_f[a], true);
    }
  }
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::impose_ic(std::vector<Sarray>& a_U, float_sw4 t,
                                      std::vector<Sarray*>& a_AlphaVE) {
  SW4_MARK_FUNCTION;
#ifdef SW4_NORM_TRACE
  static int myRank = -1;
  static bool first = true;
  if (myRank==-1)
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  static std::ofstream ofile;
  if ((!myRank) &&(first)){
    ofile.open("CIC.dat");
    first=false;
  }
#endif
  // SYNC_STREAM;                   // CURVI_CPU
  bool force_dirichlet = false;  //, check_stress_cont=false;
  //   int fg=0;
  //   if( force_dirichlet )
  //      fg = 1;
  SW4_MARK_BEGIN("IMPOSE_IC_1");

  Sarray U_f(3, m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef, __FILE__, __LINE__);
  Sarray U_c(3, m_ib, m_ie, m_jb, m_je, m_kb, m_ke, __FILE__, __LINE__);
  vector<Sarray> Alpha_c, Alpha_f;

  //  1. copy   a_U into U_f and U_c
  // std::cout<<"HERE 3\n";
  U_f.insert_intersection(a_U[m_gf]);
  U_c.insert_intersection(a_U[m_gc]);
  // std::cout<<"Initial UC "<<U_c.norm()<<"\n"<<std::flush;
  if (m_use_attenuation) {
    Alpha_c.resize(m_number_mechanisms);
    Alpha_f.resize(m_number_mechanisms);
    //      Alpha_c = new Sarray[m_number_mechanisms];
    //      Alpha_f = new Sarray[m_number_mechanisms];
    // std::cout<<"BRACKET OPEN\n"<<std::flush;
    for (int a = 0; a < m_number_mechanisms; a++) {
      Alpha_f[a].define(3, m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef,
                        Space::Managed_temps);
      Alpha_c[a].define(3, m_ib, m_ie, m_jb, m_je, m_kb, m_ke,
                        Space::Managed_temps);
      Alpha_f[a].insert_intersection(a_AlphaVE[m_gf][a]);
      Alpha_c[a].insert_intersection(a_AlphaVE[m_gc][a]);
    }
    // std::cout<<"BRACKET CLOSE\n"<<std::flush;
  }
  SW4_MARK_END("IMPOSE_IC_1");
  SW4_MARK_BEGIN("IMPOSE_IC_2");
  // 2a. Impose dirichlet conditions at ghost points
  int sides[6] = {m_isbndry[0], m_isbndry[1], m_isbndry[2], m_isbndry[3], 0, 0};
  if (m_tw != 0) {
    // std::cout<<"THIS \n"<<std::flush;
    m_tw->get_ubnd(U_f, m_x_f, m_y_f, m_z_f, t, m_nghost, sides);
    m_tw->get_ubnd(U_c, m_x_c, m_y_c, m_z_c, t, m_nghost, sides);
    if (m_use_attenuation) {
      // std::cout<<"AND THIS \n"<<std::flush;
      m_tw->get_bnd_att(Alpha_f[0], m_x_f, m_y_f, m_z_f, t, m_nghost, sides);
      m_tw->get_bnd_att(Alpha_c[0], m_x_c, m_y_c, m_z_c, t, m_nghost, sides);
    }
    if (force_dirichlet) {
      // Debug
      sides[0] = sides[1] = sides[2] = sides[3] = sides[4] = 0;
      sides[5] = 1;
      m_tw->get_ubnd(U_f, m_x_f, m_y_f, m_z_f, t, m_nghost + 1, sides);
      sides[0] = sides[1] = sides[2] = sides[3] = sides[5] = 0;
      sides[4] = 1;
      m_tw->get_ubnd(U_c, m_x_c, m_y_c, m_z_c, t, m_nghost + 1, sides);
      a_U[m_gc].copy_kplane2(U_c, 0);      // have computed U_c:s ghost points
      a_U[m_gc].copy_kplane2(U_c, 1);      // have computed U_c:s ghost points
      a_U[m_gf].copy_kplane2(U_f, m_nkf);  // .. and U_f:s interface points
      a_U[m_gf].copy_kplane2(U_f, m_nkf + 1);  // .. and U_f:s ghost point
      return;
      // End debug
    }
  } else if (m_etest != 0) {
    m_etest->get_ubnd(U_f, m_nghost, sides);
    m_etest->get_ubnd(U_c, m_nghost, sides);
  } else if (m_psource != 0) {
    m_psource->ubnd(U_f, m_x_f, m_y_f, m_z_f, t, m_ew->mGridSize[m_gf],
                    m_nghost, sides);
    m_psource->ubnd(U_c, m_x_c, m_y_c, m_z_c, t, m_ew->mGridSize[m_gc],
                    m_nghost, sides);
  } else {
    bnd_zero(U_c, m_nghost);
    bnd_zero(U_f, m_nghost);
    if (m_use_attenuation)
      for (int a = 0; a < m_number_mechanisms; a++) {
        bnd_zero(Alpha_c[a], m_nghost);
        bnd_zero(Alpha_f[a], m_nghost);
      }
  }
  SW4_MARK_END("IMPOSE_IC_2");
  SW4_MARK_BEGIN("IMPOSE_IC_3");
  // 3. Inject U_f := U_c on interface
  communicate_array(U_c, true);
  injection(U_f, U_c);
  communicate_array(U_f, true);

  if (m_use_attenuation)
    for (int a = 0; a < m_number_mechanisms; a++) {
      communicate_array(Alpha_c[a], true);
      injection(Alpha_f[a], Alpha_c[a]);
      communicate_array(Alpha_f[a], true);
    }

  // 4. Solve equation for stress continuity, formulated as lhs*x+rhs=0, where x
  // are uc's ghost points at k=0.

  // 4.a Form right hand side of equation
  Sarray rhs(3, m_ib, m_ie, m_jb, m_je, 1, 1, __FILE__, __LINE__);
  //  int asize = 3*(m_ie-m_ib+1)*(m_je-m_jb+1);
  // std::cout<<"RHS SIZE "<<m_ib<<" "<<m_ie<<" "<<m_jb<<" "<<m_je<<" size =
  // "<<asize<<"\n";
  interface_rhs(rhs, U_c, U_f, Alpha_c, Alpha_f);
  SW4_MARK_END("IMPOSE_IC_3");
  SW4_MARK_BEGIN("IMPOSE_IC_4");
  // 4.b Left hand side, lhs*x
  Sarray lhs(rhs, Space::Managed_temps), residual(rhs, Space::Managed_temps);
  interface_lhs(lhs, U_c);

  // Initial residual
  float_sw4 maxresloc = 0;
  // for (int c = 1; c <= 3; c++)
  //   for (int j = lhs.m_jb + 5; j <= lhs.m_je - 5; j++)
  //     for (int i = lhs.m_ib + 5; i <= lhs.m_ie - 5; i++) {
  //       residual(c, i, j, 1) = lhs(c, i, j, 1) + rhs(c, i, j, 1);
  //       if (abs(residual(c, i, j, 1)) > maxresloc)
  //         maxresloc = abs(residual(c, i, j, 1));
  //     }

  RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> rmax(0);
  SView& residualV = residual.getview();
  SView& lhsV = lhs.getview();
  SView& rhsV = rhs.getview();
  RAJA::RangeSegment j_range(lhs.m_jb + 5, lhs.m_je - 4);
  RAJA::RangeSegment i_range(lhs.m_ib + 5, lhs.m_ie - 4);
  RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
        for (int c = 1; c <= 3; c++) {
          residualV(c, i, j, 1) = lhsV(c, i, j, 1) + rhsV(c, i, j, 1);
          rmax.max(fabs(residualV(c, i, j, 1)));
        }
      });
// float_sw4 nmax=-1.234;
//for (int i=lhs.m_ib + 5;i<lhs.m_ie - 4;i++) for(int j=lhs.m_jb + 5;j<lhs.m_je - 4;j++) for(int c=1;c<=3;c++) nmax = std::max(nmax,fabs(residual(c, i, j, 1)));
  // std::cout<<"LHS RHS "<<lhs.norm()<<" "<<rhs.norm()<<"\n"<<std::flush;
  maxresloc = static_cast<float_sw4>(rmax.get());
  // std::cout<<"MAX RES LOC "<<maxresloc<<"\n"<<std::flush;
  float_sw4 maxres = maxresloc;
  MPI_Allreduce(&maxresloc, &maxres, 1, m_ew->m_mpifloat, MPI_MAX,
                m_ew->m_cartesian_communicator);
  SW4_MARK_END("IMPOSE_IC_4");

  SW4_MARK_BEGIN("IMPOSE_IC_JACOBI");
  // 4.c Jacobi iteration
  float_sw4 scalef =
      (m_ew->m_global_nx[m_gc] - 1) *
      (m_ew->m_global_ny[m_gc] - 1);  // scale residual to be size O(1).
  //   int maxit = 30;
  //   float_sw4 reltol=1e-6, abstol=1e-6;
  int iter = 0;
  // int three = 3, one = 1;

  int nimb = m_Mass_block.m_ie - m_Mass_block.m_ib + 1;
  // Block Jacobi, lhs*x+rhs=0 and lhs=M+N --> M*xp+N*x+rhs=0 -->
  // M*(xp-x)+lhs*x+rhs=0
  //        --> xp-x=-inv(M)*(lhs*x+rhs) --> xp = x - inv(M)*(lhs*x+rhs)
  float_sw4 maxres0 = maxres;
  float_sw4 relax = 1.0;
  //   if(  convhist.size() == 0 )
  //   {
  //      convhist.push_back(sfact);
  //      convhist.push_back(reltol);
  //      convhist.push_back(abstol);
  //   }
#ifdef SW4_NORM_TRACE
  if (!myRank) ofile<<"WHILE "<<maxres<<" "<<maxresloc<<" "<<lhs.norm()<<" "<<rhs.norm()<<"\n";
#endif
  while (maxres > m_reltol * maxres0 && scalef * maxres > m_abstol &&
         iter <= m_maxit) {
    iter++;
    //std::cout << "Iteration " << iter << " " << scalef*maxres << "\n";

#ifdef USE_DIRECT_INVERSE
    int l_ib = m_Mass_block.m_ib;
    int l_ie = m_Mass_block.m_ie;
    int l_jb = m_Mass_block.m_jb;
    int l_je = m_Mass_block.m_je;
    SView& residualV = residual.getview();
    RAJA::RangeSegment j_range(l_jb, l_je + 1);
    RAJA::RangeSegment i_range(l_ib, l_ie + 1);

    float_sw4* lx = x;
    float_sw4* lm_mass_block = m_mass_block;
    RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
        RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
          // for (int j = m_Mass_block.m_jb; j <= m_Mass_block.m_je; j++)
          //   for (int i = m_Mass_block.m_ib; i <= m_Mass_block.m_ie; i++) {
          size_t ind = (i - l_ib) + nimb * (j - l_jb);
          for (int l = 1; l < 4; l++)
            lx[l - 1 + 3 * ind] = residualV(l, i, j, 1);
        });
    size_t batchsize = (l_ie - l_ib) + nimb * (l_je - l_jb) + 1;
    // std::cout<<"Batch size in solve is "<<batchsize<<"\n";
    SW4_MARK_BEGIN("MATVEC");
    RAJA::forall<DEFAULT_LOOP1_ASYNC>(
        RAJA::RangeSegment(0, batchsize), [=] RAJA_DEVICE(int l) {
          int base = l * 9;
          float_sw4 sum = 0.0;
          float_sw4 res[3];
          for (int k = 0; k < 3; k++) {
            sum = 0.0;
            for (int i = 0; i < 3; i++)
              sum += lm_mass_block[base + k * 3 + i] * lx[l * 3 + i];
            res[k] = sum;
          }
          for (int k = 0; k < 3; k++) lx[l * 3 + k] = res[k];
        });
    SW4_MARK_END("MATVEC");
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
#endif
#ifdef USE_MAGMA
    //    int lc=0;
    int l_ib = m_Mass_block.m_ib;
    int l_ie = m_Mass_block.m_ie;
    int l_jb = m_Mass_block.m_jb;
    int l_je = m_Mass_block.m_je;
    SView& residualV = residual.getview();
    
    RAJA::RangeSegment j_range(l_jb, l_je + 1);
    RAJA::RangeSegment i_range(l_ib, l_ie + 1);

    float_sw4* lx = x;
    RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
        RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
          // for (int j = m_Mass_block.m_jb; j <= m_Mass_block.m_je; j++)
          //   for (int i = m_Mass_block.m_ib; i <= m_Mass_block.m_ie; i++) {
          size_t ind = (i - l_ib) + nimb * (j - l_jb);
          for (int l = 1; l < 4; l++)
            lx[l - 1 + 3 * ind] = residualV(l, i, j, 1);
        });
    //#define PEEKS_GALORE
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
    // SYNC_STREAM;

    // magma_int_t m_three = 3;
    // magma_int_t m_one = 1;
    //  magma_int_t batchsize=(l_ie-l_ib)+nimb*(l_je-l_jb)+1;
    // std::cout<<"BATCH SIZE IN SOLVE IS "<<batchsize<<"\n"<<std::flush;
    int three = 3, one = 1;
    magma_trans_t m_trans = MagmaNoTrans;
    SW4_MARK_BEGIN("MAGMA_DGETRS");
    // magma_int_t retval =
    // magma_dgetrs_batched(m_trans,three,one,dA_array,three,piv_array,dB_array,three,batchsize,queue);
    magma_int_t retval;
    for (int i = 0; i < subbatchsize.size(); i++)
      retval = magma_dgetrs_batched(
          m_trans, three, one, dA_array + subbatchoffset[i], three,
          piv_array + subbatchoffset[i], dB_array + subbatchoffset[i], three,
          subbatchsize[i], queue);
    if (retval != MAGMA_SUCCESS) {
      std::cerr << "MAGMA DGETRS BATCHED call failed" << std::flush;
      abort();
    }
    SW4_MARK_END("MAGMA_DGETRS");
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
#endif
#ifndef USE_LAPACK_ON_CPU
    //   lc=0;
    // for (int j = m_Mass_block.m_jb; j <= m_Mass_block.m_je; j++)
    //       for (int i = m_Mass_block.m_ib; i <= m_Mass_block.m_ie; i++) {
    SView& U_cV = U_c.getview();

    RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
        RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
          // for (int j = m_Mass_block.m_jb; j <= m_Mass_block.m_je; j++)
          //   for (int i = m_Mass_block.m_ib; i <= m_Mass_block.m_ie; i++) {
          size_t ind = (i - l_ib) + nimb * (j - l_jb);
          residualV(1, i, j, 1) = lx[3 * ind + 0];
          residualV(2, i, j, 1) = lx[3 * ind + 1];
          residualV(3, i, j, 1) = lx[3 * ind + 2];
          U_cV(1, i, j, 0) -= relax * residualV(1, i, j, 1);
          U_cV(2, i, j, 0) -= relax * residualV(2, i, j, 1);
          U_cV(3, i, j, 0) -= relax * residualV(3, i, j, 1);
        });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
#endif
#ifdef USE_LAPACK_ON_CPU
    int one = 1;
    int three = 3;

    // WARNING THIS IS RUNNING ON THE HOST
    int info = 0;
    char trans = 'N';
    for (int j = m_Mass_block.m_jb; j <= m_Mass_block.m_je; j++)
      for (int i = m_Mass_block.m_ib; i <= m_Mass_block.m_ie; i++) {
        size_t ind = (i - m_Mass_block.m_ib) + nimb * (j - m_Mass_block.m_jb);
        float_sw4 x[3] = {residual(1, i, j, 1), residual(2, i, j, 1),
                          residual(3, i, j, 1)};
        SW4_MARK_BEGIN("DGETRS");
        F77_FUNC(dgetrs, DGETRS)
        (&trans, &three, &one, &m_mass_block[9 * ind], &three,
         &m_ipiv_block[3 * ind], x, &three, &info);
        SW4_MARK_END("DGETRS");
        if (info != 0) {
          std::cerr << "SOLVE Fails at (i,j) equals" << i << "," << j
                    << " INFO = " << info << " "
                    << m_Mass_block(info + 3 * (info - 1), i, j, 1) << "\n";
          abort();
        }
        residual(1, i, j, 1) = x[0];
        residual(2, i, j, 1) = x[1];
        residual(3, i, j, 1) = x[2];
        U_c(1, i, j, 0) -= relax * residual(1, i, j, 1);
        U_c(2, i, j, 0) -= relax * residual(2, i, j, 1);
        U_c(3, i, j, 0) -= relax * residual(3, i, j, 1);
      }
#endif

#ifdef SW4_NORM_TRACE
    if (!myRank) ofile<<"PRE-INTERFACE U_C"<<U_c.norm()<<"\n";
#endif

    // 4.d Communicate U_c here (only k=0 plane)
    communicate_array(U_c, false, 0);
    interface_lhs(lhs, U_c);

    // 4.e. Compute residual and its norm
    maxresloc = 0;
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> raja_maxresloc(0);

    // for (int j = lhs.m_jb + 5; j <= lhs.m_je - 5; j++)
    //   for (int i = lhs.m_ib + 5; i <= lhs.m_ie - 5; i++) {
    // RAJA::RangeSegment j_range(lhs.m_jb + 5,lhs.m_je - 4);
    // RAJA::RangeSegment i_range(lhs.m_ib + 5,lhs.m_ie - 4);
    RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
        RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
          for (int c = 1; c <= 3; c++) {
            residualV(c, i, j, 1) = lhsV(c, i, j, 1) + rhsV(c, i, j, 1);
            raja_maxresloc.max(fabs(residualV(c, i, j, 1)));
            // if (abs(residual(c, i, j, 1)) > maxresloc)
            //   maxresloc = abs(residual(c, i, j, 1));
          }
        });
    maxresloc = static_cast<float_sw4>(raja_maxresloc.get());
    MPI_Allreduce(&maxresloc, &maxres, 1, m_ew->m_mpifloat, MPI_MAX,
                  m_ew->m_cartesian_communicator);
  }
  SW4_MARK_END("IMPOSE_IC_JACOBI");
  //   convhist.push_back(maxres0);
  //   convhist.push_back(maxres);
  //   convhist.push_back(it);
  if (maxres > m_reltol * maxres0 && scalef * maxres > m_abstol) {
    std::cerr << "WARNING, no convergence in curvilinear interface, res = "
              << maxres << " reltol= " << m_reltol
              << " initial res = " << maxres0 << std::endl;
    std::cerr << "     scaled res = " << scalef * maxres
              << " abstol= " << m_abstol << std::endl;
  }
#ifdef SW4_NORM_TRACE
  //int myRank = 0;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  //std::ofstream ofile; 
  if (!myRank){
    //ofile.open("CIC.dat", std::ios_base::app);
    ofile<< "CIC_END U_c="<<U_c.norm()<<" U_f "<<U_f.norm()<<"\n\n"<<std::flush;
    //ofile.close();
  }
#endif
  
  SW4_MARK_BEGIN("IMPOSE_IC_5");
  // 5. Copy U_c and U_f back to a_U, only k=0 for U_c and k=n3f for U_f.
  // std::cout<<"IMPOSIC 1 "<<a_U[m_gc].norm()<<" "<<a_U[m_gf].norm()<<"\n";
  a_U[m_gc].copy_kplane2(U_c, 0);      // have computed U_c:s ghost points
  a_U[m_gf].copy_kplane2(U_f, m_nkf);  // .. and U_f:s interface points
  // std::cout<<"IMPOSIC 2 "<<a_U[m_gc].norm()<<" "<<a_U[m_gf].norm()<<"\n";
  if (m_use_attenuation) {
    for (int a = 0; a < m_number_mechanisms; a++)
      a_AlphaVE[m_gf][a].copy_kplane2(Alpha_f[a], m_nkf);
  }
  SW4_MARK_END("IMPOSE_IC_5");
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::injection(Sarray& u_fA, Sarray& u_cA) {
  SW4_MARK_FUNCTION;
  // Injection at the interface

  const float_sw4 a = 9.0 / 16;
  const float_sw4 b = -1.0 / 16;
  //  const int ngh = m_nghost;
  int i1 = u_cA.m_ib + m_nghost - 1, i2 = u_cA.m_ie - m_nghost + 1;
  int j1 = u_cA.m_jb + m_nghost - 1, j2 = u_cA.m_je - m_nghost + 1;
  if (m_isbndry[0]) i1++;
  if (m_isbndry[1]) i2 -= 2;
  if (m_isbndry[2]) j1++;
  if (m_isbndry[3]) j2 -= 2;
  SView& u_f = u_fA.getview();
  SView& u_c = u_cA.getview();

  // for (int l = 1; l <= u_cA.m_nc; l++)
  //   //    for (int j = u_c.m_jb+ngh-1; j <= u_c.m_je-ngh+1; j++)
  //   //      for (int i = u_c.m_ib+ngh-1; i <= u_c.m_ie-ngh+1; i++)
  //   for (int j = j1; j <= j2; j++)
  //     for (int i = i1; i <= i2; i++) {

  auto& lm_nkf = m_nkf;
  RAJA::RangeSegment l_range(1, u_cA.m_nc + 1);
  RAJA::RangeSegment j_range(j1, j2 + 1);
  RAJA::RangeSegment i_range(i1, i2 + 1);
  RAJA::kernel<INJ_POL_ASYNC>(
      RAJA::make_tuple(l_range, j_range, i_range),
      [=] RAJA_DEVICE(int l, int j, int i) {
        u_f(l, 2 * i - 1, 2 * j - 1, lm_nkf) = u_c(l, i, j, 1);
        u_f(l, 2 * i, 2 * j - 1, lm_nkf) =
            b * u_c(l, i - 1, j, 1) + a * u_c(l, i, j, 1) +
            a * u_c(l, i + 1, j, 1) + b * u_c(l, i + 2, j, 1);
        u_f(l, 2 * i - 1, 2 * j, lm_nkf) =
            b * u_c(l, i, j - 1, 1) + a * u_c(l, i, j, 1) +
            a * u_c(l, i, j + 1, 1) + b * u_c(l, i, j + 2, 1);
        u_f(l, 2 * i, 2 * j, lm_nkf) =
            b * (b * u_c(l, i - 1, j - 1, 1) + a * u_c(l, i, j - 1, 1) +
                 a * u_c(l, i + 1, j - 1, 1) + b * u_c(l, i + 2, j - 1, 1)) +
            a * (b * u_c(l, i - 1, j, 1) + a * u_c(l, i, j, 1) +
                 a * u_c(l, i + 1, j, 1) + b * u_c(l, i + 2, j, 1)) +
            a * (b * u_c(l, i - 1, j + 1, 1) + a * u_c(l, i, j + 1, 1) +
                 a * u_c(l, i + 1, j + 1, 1) + b * u_c(l, i + 2, j + 1, 1)) +
            b * (b * u_c(l, i - 1, j + 2, 1) + a * u_c(l, i, j + 2, 1) +
                 a * u_c(l, i + 1, j + 2, 1) + b * u_c(l, i + 2, j + 2, 1));
      });

  if (m_isbndry[1]) {
    int i = i2 + 1;
    // for (int l = 1; l <= u_c.m_nc; l++)
    //   for (int j = j1; j <= j2; j++) {

    RAJA::kernel<INJ_POL2_ASYNC>(
        RAJA::make_tuple(l_range, j_range), [=] RAJA_DEVICE(int l, int j) {
          u_f(l, 2 * i - 1, 2 * j - 1, lm_nkf) = u_c(l, i, j, 1);
          u_f(l, 2 * i - 1, 2 * j, lm_nkf) =
              b * u_c(l, i, j - 1, 1) + a * u_c(l, i, j, 1) +
              a * u_c(l, i, j + 1, 1) + b * u_c(l, i, j + 2, 1);
        });
  }
  if (m_isbndry[3]) {
    int j = j2 + 1;
    // for (int l = 1; l <= u_cA.m_nc; l++)
    //   for (int i = i1; i <= i2; i++) {
    RAJA::kernel<INJ_POL2_ASYNC>(
        RAJA::make_tuple(l_range, i_range), [=] RAJA_DEVICE(int l, int i) {
          u_f(l, 2 * i - 1, 2 * j - 1, lm_nkf) = u_c(l, i, j, 1);
          u_f(l, 2 * i, 2 * j - 1, lm_nkf) =
              b * u_c(l, i - 1, j, 1) + a * u_c(l, i, j, 1) +
              a * u_c(l, i + 1, j, 1) + b * u_c(l, i + 2, j, 1);
        });
  }
  if (m_isbndry[3] && m_isbndry[1]) {
    int i = i2 + 1;
    int j = j2 + 1;
    // for (int l = 1; l <= u_c.m_nc; l++)
    RAJA::forall<DEFAULT_LOOP1_ASYNC>(
        RAJA::RangeSegment(1, u_cA.m_nc + 1), [=] RAJA_DEVICE(int l) {
          u_f(l, 2 * i - 1, 2 * j - 1, lm_nkf) = u_c(l, i, j, 1);
        });
  }
  // SYNC_STREAM;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::interface_block(Sarray& matrix) {
  SW4_MARK_FUNCTION;
  const float_sw4 w1 = 17.0 / 48;
  matrix_Lu(m_ib, m_jb, matrix, m_met_c, m_jac_c, m_mu_c, m_lambda_c, m_strx_c,
            m_stry_c, m_ghcof[0]);

  for (int c = 1; c <= 9; c++)
    for (int j = matrix.m_jb; j <= matrix.m_je; j++)
      for (int i = matrix.m_ib; i <= matrix.m_ie; i++)
        matrix(c, i, j, 1) /= m_rho_c(i, j, 1);

  Sarray alpha(m_ibf, m_ief, m_jbf, m_jef, m_nkf, m_nkf);
  for (int j = alpha.m_jb; j <= alpha.m_je; j++)
    for (int i = alpha.m_ib; i <= alpha.m_ie; i++)
      alpha(i, j, m_nkf) = w1 * m_jac_f(i, j, m_nkf) * m_rho_f(i, j, m_nkf) /
                           (m_strx_f[i - m_ibf] * m_stry_f[j - m_jbf]);
  if (!m_tw && !m_psource)
    bnd_zero_host(alpha, m_nghost);  // This should be done on host

  restprol2D(matrix, alpha, 1, m_nkf);

  // Add -B(uc) contribution to block matrix
  mat_icstresses_curv(m_ib, m_jb, matrix, 1, m_met_c, m_mu_c, m_lambda_c,
                      m_strx_c, m_stry_c, m_sbop);
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::interface_lhs(Sarray& lhs, Sarray& uc) {
  SW4_MARK_FUNCTION;

  const float_sw4 w1 = 17.0 / 48;
  lhs_Lu(uc, lhs, m_met_c, m_jac_c, m_mu_c, m_lambda_c, m_strx_c, m_stry_c,
         m_ghcof[0]);

  // for (int c = 1; c <= 3; c++)
  //   for (int j = lhs.m_jb; j <= lhs.m_je; j++)
  //     for (int i = lhs.m_ib; i <= lhs.m_ie; i++)
  SView& lhsV = lhs.getview();
  SView& m_rho_cV = m_rho_c.getview();
  RAJA::RangeSegment j_range(lhs.m_jb, lhs.m_je + 1);
  RAJA::RangeSegment i_range(lhs.m_ib, lhs.m_ie + 1);
  RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
#pragma unroll
        for (int c = 1; c <= 3; c++) lhsV(c, i, j, 1) /= m_rho_cV(i, j, 1);
      });
  // std::cout<<"interface_lhs 1 "<<lhs.norm()<<"
  // "<<m_rho_c.norm()<<"\n"<<std::flush;
  if (!m_tw && !m_psource) bnd_zero(lhs, m_nghost);

  Sarray prollhs(3, m_ibf, m_ief, m_jbf, m_jef, m_nkf, m_nkf, __FILE__,
                 __LINE__);
  prolongate2D(lhs, prollhs, 1, m_nkf);
  // std::cout<<"interface_lhs 2 "<<lhs.norm()<<"
  // "<<prollhs.norm()<<"\n"<<std::flush;
  // for (int c = 1; c <= 3; c++)
  //   for (int j = prollhs.m_jb; j <= prollhs.m_je; j++)
  //     for (int i = prollhs.m_ib; i <= prollhs.m_ie; i++)
  auto& prollhsV = prollhs.getview();
  auto& m_jac_fV = m_jac_f.getview();
  auto& m_rho_fV = m_rho_f.getview();
  float_sw4* lm_strx_f = m_strx_f;
  float_sw4* lm_stry_f = m_stry_f;
  int lm_ibf = m_ibf;
  int lm_jbf = m_jbf;
  int lm_nkf = m_nkf;
 // std::cout<<"PRE PROLL"<<prollhs.norm()<<"\n";
  RAJA::RangeSegment j_range2(prollhs.m_jb, prollhs.m_je + 1);
  RAJA::RangeSegment i_range2(prollhs.m_ib, prollhs.m_ie + 1);
  RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
      RAJA::make_tuple(j_range2, i_range2), [=] RAJA_DEVICE(int j, int i) {
#pragma unroll
        for (int c = 1; c <= 3; c++)
          prollhsV(c, i, j, lm_nkf) =
              w1 * m_jac_fV(i, j, lm_nkf) * m_rho_fV(i, j, lm_nkf) *
              prollhsV(c, i, j, lm_nkf) /
              (lm_strx_f[i - lm_ibf] * lm_stry_f[j - lm_jbf]);
	//	printf("LOOP %d %d %lf %lf\n",i,j,lm_strx_f[i - lm_ibf] ,lm_stry_f[j - lm_jbf]);
      });
 // std::cout<<"JAC = "<<m_jac_f.norm()<<" RHOF = "<<m_rho_f.norm()<<"\n";
  //std::cout<<"interface_lhs 3 "<<lhs.norm()<<"   "<<prollhs.norm()<<"\n"<<std::flush;
  if (!m_tw && !m_psource) bnd_zero(prollhs, m_nghost);
  restrict2D(lhs, prollhs, 1, m_nkf);
  // std::cout<<"interface_lhs 4 "<<lhs.norm()<<"
  // "<<prollhs.norm()<<"\n"<<std::flush;
  Sarray Bc(lhs, Space::Managed_temps);
  // Bc.set_to_zero(); // For debugging
  lhs_icstresses_curv(uc, Bc, 1, m_met_c, m_mu_c, m_lambda_c, m_strx_c,
                      m_stry_c, m_sbop);
  // for (int c = 1; c <= 3; c++)
  //   for (int j = lhs.m_jb; j <= lhs.m_je; j++)
  //     for (int i = lhs.m_ib; i <= lhs.m_ie; i++)
  auto& BcV = Bc.getview();
  RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
#pragma unroll
        for (int c = 1; c <= 3; c++) lhsV(c, i, j, 1) -= BcV(c, i, j, 1);
      });
  // std::cout<<"interface_lhs 5 "<<lhs.norm()<<"
  // "<<Bc.norm()<<"\n"<<std::flush;
  // SYNC_STREAM;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::interface_rhs(Sarray& rhs, Sarray& uc, Sarray& uf,
                                          vector<Sarray>& Alpha_c,
                                          vector<Sarray>& Alpha_f) {
  SW4_MARK_FUNCTION;

  Sarray utmp(3, uc.m_ib, uc.m_ie, uc.m_jb, uc.m_je, 0, 0, __FILE__, __LINE__);

  const float_sw4 w1 = 17.0 / 48;
  //  1. Set ghost points to zero, and save the old value to restore after, so
  //     that the routine does not change uc.
  // for (int c = 1; c <= 3; c++)
  //   for (int j = uc.m_jb; j <= uc.m_je; j++)
  //     for (int i = uc.m_ib; i <= uc.m_ie; i++) {
  SW4_MARK_BEGIN("SAVE_GHOSTS");
  SView& ucV = uc.getview();
  SView& utmpV = utmp.getview();
  RAJA::RangeSegment j_range(uc.m_jb, uc.m_je + 1);
  RAJA::RangeSegment i_range(uc.m_ib, uc.m_ie + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(RAJA::make_tuple(j_range, i_range),
                                     [=] RAJA_DEVICE(int j, int i) {
                                       for (int c = 1; c <= 3; c++) {
                                         utmpV(c, i, j, 0) = ucV(c, i, j, 0);
                                         ucV(c, i, j, 0) = 0;
                                       }
                                     });
  SW4_MARK_END("SAVE_GHOSTS");
  // std::cout<<"CALL TO CURV4SGWIND A\n";
  int onesided[6] = {0, 0, 0, 0, 1, 1};
  // 2. Compute L(uc)/rhoc
  curvilinear4sgwind(m_ib, m_ie, m_jb, m_je, m_kb, m_ke, 1, 1, uc.c_ptr(),
                     m_mu_c.c_ptr(), m_lambda_c.c_ptr(), m_met_c.c_ptr(),
                     m_jac_c.c_ptr(), rhs.c_ptr(), onesided, m_acof, m_bope,
                     m_ghcof, m_acof_no_gp, m_ghcof_no_gp, m_strx_c, m_stry_c,
                     8, '=');
  if (m_use_attenuation)
    for (int a = 0; a < m_number_mechanisms; a++)
      curvilinear4sgwind(m_ib, m_ie, m_jb, m_je, m_kb, m_ke, 1, 1,
                         Alpha_c[a].c_ptr(), m_muve_c[a].c_ptr(),
                         m_lambdave_c[a].c_ptr(), m_met_c.c_ptr(),
                         m_jac_c.c_ptr(), rhs.c_ptr(), onesided, m_acof_no_gp,
                         m_bope, m_ghcof_no_gp, m_acof_no_gp, m_ghcof_no_gp,
                         m_strx_c, m_stry_c, 8, '-');

  // for (int c = 1; c <= 3; c++)
  //   for (int j = rhs.m_jb; j <= rhs.m_je; j++)
  //     for (int i = rhs.m_ib; i <= rhs.m_ie; i++)
  auto& rhsV = rhs.getview();
  auto& m_rho_cV = m_rho_c.getview();

  RAJA::RangeSegment j_range2(rhs.m_jb, rhs.m_je + 1);
  RAJA::RangeSegment i_range2(rhs.m_ib, rhs.m_ie + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(RAJA::make_tuple(j_range2, i_range2),
                                     [=] RAJA_DEVICE(int j, int i) {
                                       for (int c = 1; c <= 3; c++)
                                         rhsV(c, i, j, 1) /= m_rho_cV(i, j, 1);
                                     });
  // SYNC_STREAM;
  if (!m_tw && !m_psource) bnd_zero(rhs, m_nghost);

  // 3. Compute prolrhs := p(L(uc)/rhoc)
  Sarray prolrhs(3, m_ibf, m_ief, m_jbf, m_jef, m_nkf, m_nkf, __FILE__,
                 __LINE__);
  prolongate2D(rhs, prolrhs, 1, m_nkf);

  // 4. Compute L(uf)
  Sarray Luf(prolrhs, Space::Managed_temps);
  curvilinear4sgwind(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef, m_nkf, m_nkf,
                     uf.c_ptr(), m_mu_f.c_ptr(), m_lambda_f.c_ptr(),
                     m_met_f.c_ptr(), m_jac_f.c_ptr(), Luf.c_ptr(), onesided,
                     m_acof, m_bope, m_ghcof, m_acof_no_gp, m_ghcof_no_gp,
                     m_strx_f, m_stry_f, m_nkf, '=');
  if (m_use_attenuation)
    for (int a = 0; a < m_number_mechanisms; a++)
      curvilinear4sgwind(m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef, m_nkf, m_nkf,
                         Alpha_f[a].c_ptr(), m_muve_f[a].c_ptr(),
                         m_lambdave_f[a].c_ptr(), m_met_f.c_ptr(),
                         m_jac_f.c_ptr(), Luf.c_ptr(), onesided, m_acof_no_gp,
                         m_bope, m_ghcof_no_gp, m_acof_no_gp, m_ghcof_no_gp,
                         m_strx_f, m_stry_f, m_nkf, '-');

  // 5. Compute B(uf)
  Sarray Bf(prolrhs, Space::Managed_temps);
  compute_icstresses_curv(uf, Bf, m_nkf, m_met_f, m_mu_f, m_lambda_f, m_strx_f,
                          m_stry_f, m_sbop_no_gp, '=');
  if (m_use_attenuation)
    for (int a = 0; a < m_number_mechanisms; a++)
      compute_icstresses_curv(Alpha_f[a], Bf, m_nkf, m_met_f, m_muve_f[a],
                              m_lambdave_f[a], m_strx_f, m_stry_f, m_sbop_no_gp,
                              '-');

  // 6. Form term prolrhs := r(w1*J[gf]*(rhof*p(L(uc)/rhoc-L(uf))+B(uf))
  // for (int c = 1; c <= 3; c++)
  //   for (int j = prolrhs.m_jb; j <= prolrhs.m_je; j++)
  //     for (int i = prolrhs.m_ib; i <= prolrhs.m_ie; i++)
  auto& prolrhsV = prolrhs.getview();
  auto& m_jac_fV = m_jac_f.getview();
  auto& m_rho_fV = m_rho_f.getview();
  auto& LufV = Luf.getview();
  auto& BfV = Bf.getview();

  float_sw4* lm_strx_f = m_strx_f;
  float_sw4* lm_stry_f = m_stry_f;
  auto& lm_nkf = m_nkf;
  auto& lm_ibf = m_ibf;
  auto& lm_jbf = m_jbf;
  RAJA::RangeSegment j_range3(prolrhs.m_jb, prolrhs.m_je + 1);
  RAJA::RangeSegment i_range3(prolrhs.m_ib, prolrhs.m_ie + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(
      RAJA::make_tuple(j_range3, i_range3), [=] RAJA_DEVICE(int j, int i) {
        for (int c = 1; c <= 3; c++)
          prolrhsV(c, i, j, lm_nkf) =
              w1 * m_jac_fV(i, j, lm_nkf) *
                  (m_rho_fV(i, j, lm_nkf) * prolrhsV(c, i, j, lm_nkf) -
                   LufV(c, i, j, lm_nkf)) /
                  (lm_strx_f[i - lm_ibf] * lm_stry_f[j - lm_jbf]) +
              BfV(c, i, j, lm_nkf);
      });
  // SYNC_STREAM;
  if (!m_tw && !m_psource) bnd_zero(prolrhs, m_nghost);
  restrict2D(rhs, prolrhs, 1, m_nkf);

  // 7. Compute B(uc), and form rhs := rhs - B(uc) =
  // r(w1*J[gf]*(rhof*p(L(uc)/rhoc-L(uf))+B(uf))-B(uc)
  Sarray Bc(rhs, Space::Managed_temps);
  compute_icstresses_curv(uc, Bc, 1, m_met_c, m_mu_c, m_lambda_c, m_strx_c,
                          m_stry_c, m_sbop, '=');
  if (m_use_attenuation)
    for (int a = 0; a < m_number_mechanisms; a++)
      compute_icstresses_curv(Alpha_c[a], Bc, 1, m_met_c, m_muve_c[a],
                              m_lambdave_c[a], m_strx_c, m_stry_c, m_sbop_no_gp,
                              '-');

  // for (int c = 1; c <= 3; c++)
  //   for (int j = rhs.m_jb; j <= rhs.m_je; j++)
  //     for (int i = rhs.m_ib; i <= rhs.m_ie; i++)
  auto& BcV = Bc.getview();
  RAJA::RangeSegment j_range4(rhs.m_jb, rhs.m_je + 1);
  RAJA::RangeSegment i_range4(rhs.m_ib, rhs.m_ie + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(RAJA::make_tuple(j_range4, i_range4),
                                     [=] RAJA_DEVICE(int j, int i) {
#pragma unroll
                                       for (int c = 1; c <= 3; c++)
                                         rhsV(c, i, j, 1) -= BcV(c, i, j, 1);
                                     });

  // 8. Restore ghost point values to U.
  // for (int c = 1; c <= 3; c++)
  //   for (int j = uc.m_jb; j <= uc.m_je; j++)
  //     for (int i = uc.m_ib; i <= uc.m_ie; i++)
  RAJA::RangeSegment j_range5(uc.m_jb, uc.m_je + 1);
  RAJA::RangeSegment i_range5(uc.m_ib, uc.m_ie + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(RAJA::make_tuple(j_range5, i_range5),
                                     [=] RAJA_DEVICE(int j, int i) {
#pragma unroll
                                       for (int c = 1; c <= 3; c++)
                                         ucV(c, i, j, 0) = utmpV(c, i, j, 0);
                                     });

  SYNC_STREAM;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::compute_icstresses_curv_host(
    Sarray& a_Up, Sarray& B, int kic, Sarray& a_metric, Sarray& a_mu,
    Sarray& a_lambda, float_sw4* a_str_x, float_sw4* a_str_y, float_sw4* sbop,
    char op) {
  SW4_MARK_FUNCTION;
  const float_sw4 a1 = 2.0 / 3, a2 = -1.0 / 12;
  const bool upper = (kic == 1);
  const int k = kic;
  const int kl = upper ? 1 : -1;
  const int ifirst = a_Up.m_ib;
  const int jfirst = a_Up.m_jb;
#define str_x(i) a_str_x[(i - ifirst)]
#define str_y(j) a_str_y[(j - jfirst)]
  float_sw4 sgn = 1;
  if (op == '=') {
    B.set_value(0.0);
    SYNC_STREAM;
    sgn = 1;
  }
  if (op == '-') {
    sgn = -1;
  }

#pragma omp parallel for
  for (int j = B.m_jb + 2; j <= B.m_je - 2; j++)
#pragma omp simd
    for (int i = B.m_ib + 2; i <= B.m_ie - 2; i++) {
      float_sw4 uz, vz, wz;
      uz = vz = wz = 0;
      for (int m = 0; m <= 5; m++) {
        uz += sbop[m] * a_Up(1, i, j, k + kl * (m - 1));
        vz += sbop[m] * a_Up(2, i, j, k + kl * (m - 1));
        wz += sbop[m] * a_Up(3, i, j, k + kl * (m - 1));
      }
      uz *= kl;
      vz *= kl;
      wz *= kl;

      // Normal terms
      float_sw4 m2 = str_x(i) * a_metric(2, i, j, k);
      float_sw4 m3 = str_y(j) * a_metric(3, i, j, k);
      float_sw4 m4 = a_metric(4, i, j, k);
      float_sw4 un = m2 * uz + m3 * vz + m4 * wz;
      float_sw4 mnrm = m2 * m2 + m3 * m3 + m4 * m4;
      float_sw4 B1, B2, B3;

      B1 = a_mu(i, j, k) * mnrm * uz +
           (a_mu(i, j, k) + a_lambda(i, j, k)) * m2 * un;
      B2 = a_mu(i, j, k) * mnrm * vz +
           (a_mu(i, j, k) + a_lambda(i, j, k)) * m3 * un;
      B3 = a_mu(i, j, k) * mnrm * wz +
           (a_mu(i, j, k) + a_lambda(i, j, k)) * m4 * un;

      // Tangential terms
      // p-derivatives
      float_sw4 up1 =
          str_x(i) * (a2 * (a_Up(1, i + 2, j, k) - a_Up(1, i - 2, j, k)) +
                      a1 * (a_Up(1, i + 1, j, k) - a_Up(1, i - 1, j, k)));
      float_sw4 up2 =
          str_x(i) * (a2 * (a_Up(2, i + 2, j, k) - a_Up(2, i - 2, j, k)) +
                      a1 * (a_Up(2, i + 1, j, k) - a_Up(2, i - 1, j, k)));
      float_sw4 up3 =
          str_x(i) * (a2 * (a_Up(3, i + 2, j, k) - a_Up(3, i - 2, j, k)) +
                      a1 * (a_Up(3, i + 1, j, k) - a_Up(3, i - 1, j, k)));
      B1 += a_metric(1, i, j, k) *
            ((2 * a_mu(i, j, k) + a_lambda(i, j, k)) * m2 * up1 +
             a_mu(i, j, k) * (m3 * up2 + m4 * up3));
      B2 += a_metric(1, i, j, k) *
            (a_lambda(i, j, k) * m3 * up1 + a_mu(i, j, k) * m2 * up2);
      B3 += a_metric(1, i, j, k) *
            (a_lambda(i, j, k) * m4 * up1 + a_mu(i, j, k) * m2 * up3);

      // q-derivatives
      float_sw4 uq1 =
          str_y(j) * (a2 * (a_Up(1, i, j + 2, k) - a_Up(1, i, j - 2, k)) +
                      a1 * (a_Up(1, i, j + 1, k) - a_Up(1, i, j - 1, k)));
      float_sw4 uq2 =
          str_y(j) * (a2 * (a_Up(2, i, j + 2, k) - a_Up(2, i, j - 2, k)) +
                      a1 * (a_Up(2, i, j + 1, k) - a_Up(2, i, j - 1, k)));
      float_sw4 uq3 =
          str_y(j) * (a2 * (a_Up(3, i, j + 2, k) - a_Up(3, i, j - 2, k)) +
                      a1 * (a_Up(3, i, j + 1, k) - a_Up(3, i, j - 1, k)));
      B1 += a_metric(1, i, j, k) *
            (a_lambda(i, j, k) * m2 * uq2 + a_mu(i, j, k) * m3 * uq1);
      B2 += a_metric(1, i, j, k) *
            ((2 * a_mu(i, j, k) + a_lambda(i, j, k)) * m3 * uq2 +
             a_mu(i, j, k) * (m2 * uq1 + m4 * uq3));
      B3 += a_metric(1, i, j, k) *
            (a_lambda(i, j, k) * m4 * uq2 + a_mu(i, j, k) * m3 * uq3);

      float_sw4 isgxy = 1.0 / (str_x(i) * str_y(j));
      B1 *= isgxy;
      B2 *= isgxy;
      B3 *= isgxy;

      B(1, i, j, k) += sgn * B1;
      B(2, i, j, k) += sgn * B2;
      B(3, i, j, k) += sgn * B3;
    }
#undef str_x
#undef str_y
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::lhs_icstresses_curv(
    Sarray& a_Up, Sarray& a_lhs, int kic, Sarray& a_metric, Sarray& a_mu,
    Sarray& a_lambda, float_sw4* a_str_x, float_sw4* a_str_y, float_sw4* sbop) {
  SW4_MARK_FUNCTION;
  // As compute_icstresses_curv, but evaluates the ghost point part only
  //   const float_sw4 a1=2.0/3, a2=-1.0/12;
  const bool upper = (kic == 1);
  const int k = kic;
  // const int kl = upper ? 1 :-1;
  const int ifirst = a_Up.m_ib;
  const int jfirst = a_Up.m_jb;
#define str_x(i) a_str_x[(i - ifirst)]
#define str_y(j) a_str_y[(j - jfirst)]
  // std::cout<<" lhs_icstresses_curv 1 "<<a_lhs.norm()<<"
  // "<<sbop[0]<<"\n"<<std::flush; for(int i=0;i<5;i++) std::cout<<"
  // str["<<i<<"] = "<<str_x(i+ifirst)<<" "<<str_y(i+jfirst)<<"\n"<<std::flush;
  // #pragma omp parallel for
  //   for (int j = a_lhs.m_jb; j <= a_lhs.m_je; j++)
  // #pragma omp simd
  //     for (int i = a_lhs.m_ib; i <= a_lhs.m_ie; i++) {
  auto& a_UpV = a_Up.getview();
  auto& a_metricV = a_metric.getview();
  auto& a_muV = a_mu.getview();
  auto& a_lambdaV = a_lambda.getview();
  auto& a_lhsV = a_lhs.getview();
  RAJA::RangeSegment j_range(a_lhs.m_jb, a_lhs.m_je + 1);
  RAJA::RangeSegment i_range(a_lhs.m_ib, a_lhs.m_ie + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
        float_sw4 uz, vz, wz;
        uz = vz = wz = 0;
        if (upper) {
          uz = sbop[0] * a_UpV(1, i, j, k - 1);
          vz = sbop[0] * a_UpV(2, i, j, k - 1);
          wz = sbop[0] * a_UpV(3, i, j, k - 1);
        } else {
          uz = -sbop[0] * a_UpV(1, i, j, k + 1);
          vz = -sbop[0] * a_UpV(2, i, j, k + 1);
          wz = -sbop[0] * a_UpV(3, i, j, k + 1);
        }

        // Normal terms
        float_sw4 m2 = str_x(i) * a_metricV(2, i, j, k);
        float_sw4 m3 = str_y(j) * a_metricV(3, i, j, k);
        float_sw4 m4 = a_metricV(4, i, j, k);
        float_sw4 un = m2 * uz + m3 * vz + m4 * wz;
        float_sw4 mnrm = m2 * m2 + m3 * m3 + m4 * m4;

        a_lhsV(1, i, j, k) = a_muV(i, j, k) * mnrm * uz +
                             (a_muV(i, j, k) + a_lambdaV(i, j, k)) * m2 * un;
        a_lhsV(2, i, j, k) = a_muV(i, j, k) * mnrm * vz +
                             (a_muV(i, j, k) + a_lambdaV(i, j, k)) * m3 * un;
        a_lhsV(3, i, j, k) = a_muV(i, j, k) * mnrm * wz +
                             (a_muV(i, j, k) + a_lambdaV(i, j, k)) * m4 * un;

        float_sw4 isgxy = 1.0 / (str_x(i) * str_y(j));
        a_lhsV(1, i, j, k) *= isgxy;
        a_lhsV(2, i, j, k) *= isgxy;
        a_lhsV(3, i, j, k) *= isgxy;
      });
  // std::cout<<" lhs_icstresses_curv 2 "<<a_lhs.norm()<<"\n"<<std::flush;
  // SYNC_STREAM;
#undef str_x
#undef str_y
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::lhs_Lu(Sarray& a_Ua, Sarray& a_lhsa, Sarray& meta,
                                   Sarray& jaca, Sarray& mua, Sarray& laa,
                                   float_sw4* a_str_x, float_sw4* a_str_y,
                                   float_sw4 ghcof) {
  SW4_MARK_FUNCTION;
  const int ifirst = a_Ua.m_ib;
  const int jfirst = a_Ua.m_jb;
#define strx(i) a_str_x[(i - ifirst)]
#define stry(j) a_str_y[(j - jfirst)]

  SView& a_lhs = a_lhsa.getview();
  SView& jac = jaca.getview();
  SView& mu = mua.getview();
  SView& la = laa.getview();
  SView& met = meta.getview();
  SView& a_U = a_Ua.getview();
  RAJA::RangeSegment j_range(a_lhsa.m_jb, a_lhsa.m_je + 1);
  RAJA::RangeSegment i_range(a_lhsa.m_ib, a_lhsa.m_ie + 1);
  RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
        // for (int j = a_lhs.m_jb; j <= a_lhs.m_je; j++)
        //   for (int i = a_lhs.m_ib; i <= a_lhs.m_ie; i++) {
        float_sw4 ijac = ghcof / jac(i, j, 1);
        float_sw4 mucofu2 = ((2 * mu(i, j, 1) + la(i, j, 1)) * met(2, i, j, 1) *
                                 strx(i) * met(2, i, j, 1) * strx(i) +
                             mu(i, j, 1) * (met(3, i, j, 1) * stry(j) *
                                                met(3, i, j, 1) * stry(j) +
                                            met(4, i, j, 1) * met(4, i, j, 1)));
        float_sw4 mucofv2 = ((2 * mu(i, j, 1) + la(i, j, 1)) * met(3, i, j, 1) *
                                 stry(j) * met(3, i, j, 1) * stry(j) +
                             mu(i, j, 1) * (met(2, i, j, 1) * strx(i) *
                                                met(2, i, j, 1) * strx(i) +
                                            met(4, i, j, 1) * met(4, i, j, 1)));
        float_sw4 mucofw2 =
            ((2 * mu(i, j, 1) + la(i, j, 1)) * met(4, i, j, 1) *
                 met(4, i, j, 1) +
             mu(i, j, 1) *
                 (met(2, i, j, 1) * strx(i) * met(2, i, j, 1) * strx(i) +
                  met(3, i, j, 1) * stry(j) * met(3, i, j, 1) * stry(j)));
        float_sw4 mucofuv = (mu(i, j, 1) + la(i, j, 1)) * met(2, i, j, 1) *
                            met(3, i, j, 1) * strx(i) * stry(j);
        float_sw4 mucofuw = (mu(i, j, 1) + la(i, j, 1)) * met(2, i, j, 1) *
                            met(4, i, j, 1) * strx(i);
        float_sw4 mucofvw = (mu(i, j, 1) + la(i, j, 1)) * met(3, i, j, 1) *
                            met(4, i, j, 1) * stry(j);
        a_lhs(1, i, j, 1) =
            (mucofu2 * a_U(1, i, j, 0) + mucofuv * a_U(2, i, j, 0) +
             mucofuw * a_U(3, i, j, 0)) *
            ijac;
        a_lhs(2, i, j, 1) =
            (mucofuv * a_U(1, i, j, 0) + mucofv2 * a_U(2, i, j, 0) +
             mucofvw * a_U(3, i, j, 0)) *
            ijac;
        a_lhs(3, i, j, 1) =
            (mucofuw * a_U(1, i, j, 0) + mucofvw * a_U(2, i, j, 0) +
             mucofw2 * a_U(3, i, j, 0)) *
            ijac;
        //	       r1 += istrxy*mucofu2*u(1,i,j,0) + mucofuv*u(2,i,j,0) +
        // istry*mucofuw*u(3,i,j,0); 	       r2 += mucofuv*u(1,i,j,0) +
        // istrxy*mucofv2*u(2,i,j,0) + istrx*mucofvw*u(3,i,j,0); r3
        // += istry*mucofuw*u(1,i,j,0) + istrx*mucofvw*u(2,i,j,0) +
        // istrxy*mucofw2*u(3,i,j,0);
      });
#undef strx
#undef stry
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::mat_icstresses_curv(
    int ib, int jb, Sarray& a_mat, int kic, Sarray& a_metric, Sarray& a_mu,
    Sarray& a_lambda, float_sw4* a_str_x, float_sw4* a_str_y, float_sw4* sbop) {
  SW4_MARK_FUNCTION;
  // As compute_icstresses_curv, but evaluates the matrix multiplying the ghost
  // point part
  //   const float_sw4 a1=2.0/3, a2=-1.0/12;
  const bool upper = (kic == 1);
  const int k = kic;
  // const int kl = upper ? 1 :-1;
#define str_x(i) a_str_x[(i - ib)]
#define str_y(j) a_str_y[(j - jb)]

  float_sw4 sb = sbop[0];
  if (!upper) sb = -sb;

#pragma omp parallel for
  for (int j = a_mat.m_jb; j <= a_mat.m_je; j++)
#pragma omp simd
    for (int i = a_mat.m_ib; i <= a_mat.m_ie; i++) {
      float cof = sb / (str_x(i) * str_y(j));
      // Normal terms
      float_sw4 m2 = str_x(i) * a_metric(2, i, j, k);
      float_sw4 m3 = str_y(j) * a_metric(3, i, j, k);
      float_sw4 m4 = a_metric(4, i, j, k);
      //         float_sw4 un   = m2*uz + m3*vz + m4*wz;
      float_sw4 mnrm = m2 * m2 + m3 * m3 + m4 * m4;

      //         a_lhs(1,i,j,k) = a_mu(i,j,k)*mnrm*uz +
      //         (a_mu(i,j,k)+a_lambda(i,j,k))*m2*un;
      a_mat(1, i, j, k) -=
          cof * (a_mu(i, j, k) * mnrm +
                 (a_mu(i, j, k) + a_lambda(i, j, k)) * m2 * m2);  // dB1/du1
      a_mat(4, i, j, k) -=
          cof * (a_mu(i, j, k) + a_lambda(i, j, k)) * m2 * m3;  // dB1/du2
      a_mat(7, i, j, k) -=
          cof * (a_mu(i, j, k) + a_lambda(i, j, k)) * m2 * m4;  // dB1/du3

      //         a_lhs(2,i,j,k) = a_mu(i,j,k)*mnrm*vz +
      //         (a_mu(i,j,k)+a_lambda(i,j,k))*m3*un;
      a_mat(2, i, j, k) -=
          cof * (a_mu(i, j, k) + a_lambda(i, j, k)) * m3 * m2;  // dB2/du1
      a_mat(5, i, j, k) -=
          cof * (a_mu(i, j, k) * mnrm +
                 (a_mu(i, j, k) + a_lambda(i, j, k)) * m3 * m3);  // dB2/du2
      a_mat(8, i, j, k) -=
          cof * (a_mu(i, j, k) + a_lambda(i, j, k)) * m3 * m4;  // dB2/du3

      //         a_lhs(3,i,j,k) = a_mu(i,j,k)*mnrm*wz +
      //         (a_mu(i,j,k)+a_lambda(i,j,k))*m4*un;
      a_mat(3, i, j, k) -=
          cof * (a_mu(i, j, k) + a_lambda(i, j, k)) * m4 * m2;  // dB3/du1
      a_mat(6, i, j, k) -=
          cof * (a_mu(i, j, k) + a_lambda(i, j, k)) * m4 * m3;  // dB3/du2
      a_mat(9, i, j, k) -=
          cof * (a_mu(i, j, k) * mnrm +
                 (a_mu(i, j, k) + a_lambda(i, j, k)) * m4 * m4);  // dB3/du3

      //         float_sw4 isgxy = 1.0/(str_x(i)*str_y(j));
      //         a_lhs(1,i,j,k) *= isgxy;
      //         a_lhs(2,i,j,k) *= isgxy;
      //         a_lhs(3,i,j,k) *= isgxy;
    }
#undef str_x
#undef str_y
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::matrix_Lu(int strib, int strjb, Sarray& a_mat,
                                      Sarray& met, Sarray& jac, Sarray& mu,
                                      Sarray& la, float_sw4* a_str_x,
                                      float_sw4* a_str_y, float_sw4 ghcof) {
  SW4_MARK_FUNCTION;
#define strx(i) a_str_x[(i - strib)]
#define stry(j) a_str_y[(j - strjb)]
#pragma omp parallel for
  for (int j = a_mat.m_jb; j <= a_mat.m_je; j++)
#pragma omp simd
    for (int i = a_mat.m_ib; i <= a_mat.m_ie; i++) {
      float_sw4 ijac = ghcof / jac(i, j, 1);
      float_sw4 mucofu2 = ((2 * mu(i, j, 1) + la(i, j, 1)) * met(2, i, j, 1) *
                               strx(i) * met(2, i, j, 1) * strx(i) +
                           mu(i, j, 1) * (met(3, i, j, 1) * stry(j) *
                                              met(3, i, j, 1) * stry(j) +
                                          met(4, i, j, 1) * met(4, i, j, 1)));
      float_sw4 mucofv2 = ((2 * mu(i, j, 1) + la(i, j, 1)) * met(3, i, j, 1) *
                               stry(j) * met(3, i, j, 1) * stry(j) +
                           mu(i, j, 1) * (met(2, i, j, 1) * strx(i) *
                                              met(2, i, j, 1) * strx(i) +
                                          met(4, i, j, 1) * met(4, i, j, 1)));
      float_sw4 mucofw2 =
          ((2 * mu(i, j, 1) + la(i, j, 1)) * met(4, i, j, 1) * met(4, i, j, 1) +
           mu(i, j, 1) *
               (met(2, i, j, 1) * strx(i) * met(2, i, j, 1) * strx(i) +
                met(3, i, j, 1) * stry(j) * met(3, i, j, 1) * stry(j)));
      float_sw4 mucofuv = (mu(i, j, 1) + la(i, j, 1)) * met(2, i, j, 1) *
                          met(3, i, j, 1) * strx(i) * stry(j);
      float_sw4 mucofuw = (mu(i, j, 1) + la(i, j, 1)) * met(2, i, j, 1) *
                          met(4, i, j, 1) * strx(i);
      float_sw4 mucofvw = (mu(i, j, 1) + la(i, j, 1)) * met(3, i, j, 1) *
                          met(4, i, j, 1) * stry(j);
      a_mat(1, i, j, 1) = mucofu2 * ijac;
      a_mat(4, i, j, 1) = mucofuv * ijac;
      a_mat(7, i, j, 1) = mucofuw * ijac;

      a_mat(2, i, j, 1) = mucofuv * ijac;
      a_mat(5, i, j, 1) = mucofv2 * ijac;
      a_mat(8, i, j, 1) = mucofvw * ijac;

      a_mat(3, i, j, 1) = mucofuw * ijac;
      a_mat(6, i, j, 1) = mucofvw * ijac;
      a_mat(9, i, j, 1) = mucofw2 * ijac;
    }
#undef strx
#undef stry
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::restprol2D(Sarray& Uc, Sarray& alpha, int kc,
                                       int kf) {
  SW4_MARK_FUNCTION;
  SYNC_STREAM;
  //
  // Multiplies the diagonal element of the operator P^T*diag(alpha)*P
  //
  const float_sw4 a = (9.0 / 16);
  const float_sw4 b = -1.0 / 16;
  const float_sw4 a2 = a * a;
  const float_sw4 b2 = b * b;
  const float_sw4 a4 = a2 * a2;
  const float_sw4 b4 = b2 * b2;
  const float_sw4 a2b2 = a2 * b2;
  for (int c = 1; c <= Uc.m_nc; c++)
  //   for( int j=Uc.m_jb+5; j <= Uc.m_je-5 ;j++ )
  //      for( int i=Uc.m_ib+5; i <= Uc.m_ie-5 ;i++ )
  // Note, this routine is called with Uc = matrix block, which is
  // declared without ghost points, so need to loop over the
  // full i- and j- ranges.
#pragma omp parallel for
    for (int j = Uc.m_jb; j <= Uc.m_je; j++)
#pragma omp simd
      for (int i = Uc.m_ib; i <= Uc.m_ie; i++) {
        Uc(c, i, j, kc) =
            (alpha(2 * i - 1, 2 * j - 1, kf) +
             a2 * (alpha(2 * i - 1, 2 * j, kf) +
                   alpha(2 * i - 1, 2 * j - 2, kf) +
                   alpha(2 * i, 2 * j - 1, kf) +
                   alpha(2 * i - 2, 2 * j - 1, kf)) +
             b2 * (alpha(2 * i + 2, 2 * j - 1, kf) +
                   alpha(2 * i - 4, 2 * j - 1, kf) +
                   alpha(2 * i - 1, 2 * j + 2, kf) +
                   alpha(2 * i - 1, 2 * j - 4, kf)) +
             a4 * (alpha(2 * i, 2 * j, kf) + alpha(2 * i, 2 * j - 2, kf) +
                   alpha(2 * i - 2, 2 * j, kf) +
                   alpha(2 * i - 2, 2 * j - 2, kf)) +
             b4 * (alpha(2 * i + 2, 2 * j + 2, kf) +
                   alpha(2 * i + 2, 2 * j - 4, kf) +
                   alpha(2 * i - 4, 2 * j + 2, kf) +
                   alpha(2 * i - 4, 2 * j - 4, kf)) +
             a2b2 *
                 (alpha(2 * i, 2 * j + 2, kf) + alpha(2 * i + 2, 2 * j, kf) +
                  alpha(2 * i + 2, 2 * j - 2, kf) +
                  alpha(2 * i - 2, 2 * j + 2, kf) +
                  alpha(2 * i - 2, 2 * j - 4, kf) +
                  alpha(2 * i - 4, 2 * j - 2, kf) +
                  alpha(2 * i, 2 * j - 4, kf) + alpha(2 * i - 4, 2 * j, kf))) *
            Uc(c, i, j, kc);
      }
}
//-----------------------------------------------------------------------
void CurvilinearInterface2::prolongate2D(Sarray& Uc, Sarray& Uf, int kc,
                                         int kf) {
  SW4_MARK_FUNCTION;
  const float_sw4 i16 = 1.0 / 16;
  const float_sw4 i256 = 1.0 / 256;
  int ib1, ie1, ib2, ie2;
  if (Uf.m_ib % 2 == 0)
    ib1 = Uf.m_ib / 2 + 1;
  else
    ib1 = (Uf.m_ib + 1) / 2;
  ib1 = max(Uc.m_ib, ib1);
  if (Uf.m_ie % 2 == 0)
    ie1 = Uf.m_ie / 2;
  else
    ie1 = (Uf.m_ie + 1) / 2;
  ie1 = min(Uc.m_ie, ie1);

  if (Uf.m_ib % 2 == 0)
    ib2 = Uf.m_ib / 2;
  else
    ib2 = (Uf.m_ib + 1) / 2;
  ib2 = max(Uc.m_ib + 1, ib2);
  if (Uf.m_ie % 2 == 0)
    ie2 = Uf.m_ie / 2;
  else
    ie2 = (Uf.m_ie - 1) / 2;
  ie2 = min(Uc.m_ie - 2, ie2);

  int jb1, je1, jb2, je2;
  if (Uf.m_jb % 2 == 0)
    jb1 = Uf.m_jb / 2 + 1;
  else
    jb1 = (Uf.m_jb + 1) / 2;
  jb1 = max(Uc.m_jb, jb1);
  if (Uf.m_je % 2 == 0)
    je1 = Uf.m_je / 2;
  else
    je1 = (Uf.m_je + 1) / 2;
  je1 = min(Uc.m_je, je1);

  if (Uf.m_jb % 2 == 0)
    jb2 = Uf.m_jb / 2;
  else
    jb2 = (Uf.m_jb + 1) / 2;
  jb2 = max(Uc.m_jb + 1, jb2);
  if (Uf.m_je % 2 == 0)
    je2 = Uf.m_je / 2;
  else
    je2 = (Uf.m_je - 1) / 2;
  je2 = min(Uc.m_je - 2, je2);
  auto& lm_nc = Uf.m_nc;
  auto& UcV = Uc.getview();
  auto& UfV = Uf.getview();
  //#pragma omp parallel
  //  {
  //     for (int c = 1; c <= Uf.m_nc; c++)
  // #pragma omp for
  //       for (int j = jb1; j <= je1; j++)
  // #pragma omp simd
  //         for (int i = ib1; i <= ie1; i++)
  RAJA::RangeSegment j_range1(jb1, je1 + 1);
  RAJA::RangeSegment i_range1(ib1, ie1 + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(
      RAJA::make_tuple(j_range1, i_range1), [=] RAJA_DEVICE(int j, int i) {
        for (int c = 1; c <= lm_nc; c++)
          UfV(c, 2 * i - 1, 2 * j - 1, kf) = UcV(c, i, j, kc);
      });

  //     for (int c = 1; c <= Uf.m_nc; c++)
  // #pragma omp for
  //       for (int j = jb2; j <= je2; j++)
  // #pragma omp simd
  //         for (int i = ib1; i <= ie1; i++)
  RAJA::RangeSegment j_range2(jb2, je2 + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(
      RAJA::make_tuple(j_range2, i_range1), [=] RAJA_DEVICE(int j, int i) {
        for (int c = 1; c <= lm_nc; c++)
          UfV(c, 2 * i - 1, 2 * j, kf) =
              i16 * (-UcV(c, i, j - 1, kc) +
                     9 * (UcV(c, i, j, kc) + UcV(c, i, j + 1, kc)) -
                     UcV(c, i, j + 2, kc));
      });

  //     for (int c = 1; c <= Uf.m_nc; c++)
  // #pragma omp for
  //       for (int j = jb1; j <= je1; j++)
  // #pragma omp simd
  //         for (int i = ib2; i <= ie2; i++)
  RAJA::RangeSegment i_range2(ib2, ie2 + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(
      RAJA::make_tuple(j_range1, i_range2), [=] RAJA_DEVICE(int j, int i) {
        for (int c = 1; c <= lm_nc; c++)
          UfV(c, 2 * i, 2 * j - 1, kf) =
              i16 * (-UcV(c, i - 1, j, kc) +
                     9 * (UcV(c, i, j, kc) + UcV(c, i + 1, j, kc)) -
                     UcV(c, i + 2, j, kc));
      });

  //     for (int c = 1; c <= Uf.m_nc; c++)
  // #pragma omp for
  //       for (int j = jb2; j <= je2; j++)
  // #pragma omp simd
  //         for (int i = ib2; i <= ie2; i++)
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(
      RAJA::make_tuple(j_range2, i_range2), [=] RAJA_DEVICE(int j, int i) {
        for (int c = 1; c <= lm_nc; c++)
          UfV(c, 2 * i, 2 * j, kf) =
              i256 *
              (UcV(c, i - 1, j - 1, kc) -
               9 * (UcV(c, i, j - 1, kc) + UcV(c, i + 1, j - 1, kc)) +
               UcV(c, i + 2, j - 1, kc) +
               9 * (-UcV(c, i - 1, j, kc) +
                    9 * (UcV(c, i, j, kc) + UcV(c, i + 1, j, kc)) -
                    UcV(c, i + 2, j, kc) - UcV(c, i - 1, j + 1, kc) +
                    9 * (UcV(c, i, j + 1, kc) + UcV(c, i + 1, j + 1, kc)) -
                    UcV(c, i + 2, j + 1, kc)) +
               UcV(c, i - 1, j + 2, kc) -
               9 * (UcV(c, i, j + 2, kc) + UcV(c, i + 1, j + 2, kc)) +
               UcV(c, i + 2, j + 2, kc));
      });
  // SYNC_STREAM;
  //}
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::restrict2D(Sarray& Uc, Sarray& Uf, int kc, int kf) {
  SW4_MARK_FUNCTION;
  int icb, ice, jcb, jce;
  if (Uf.m_ib % 2 == 0)
    icb = Uf.m_ib / 2 + 2;
  else
    icb = (Uf.m_ib + 1) / 2 + 2;
  icb = max(Uc.m_ib, icb);
  if (Uf.m_ie % 2 == 0)
    ice = Uf.m_ie / 2 - 1;
  else
    ice = (Uf.m_ie - 1) / 2 - 1;
  ice = min(Uc.m_ie, ice);

  if (Uf.m_jb % 2 == 0)
    jcb = Uf.m_jb / 2 + 2;
  else
    jcb = (Uf.m_jb + 1) / 2 + 2;
  jcb = max(Uc.m_jb, jcb);
  if (Uf.m_je % 2 == 0)
    jce = Uf.m_je / 2 - 1;
  else
    jce = (Uf.m_je - 1) / 2 - 1;
  jce = min(Uc.m_je, jce);

  const float_sw4 i1024 = 4.0 / 1024;  // Multiply r:=4*r
  int lm_nc = Uf.m_nc;
  auto& UfV = Uf.getview();
  auto& UcV = Uc.getview();
  // #pragma omp parallel
  //   for (int c = 1; c <= Uf.m_nc; c++)
  // #pragma omp for
  //     for (int jc = jcb; jc <= jce; jc++)
  // #pragma omp simd
  //       for (int ic = icb; ic <= ice; ic++) {
  RAJA::RangeSegment j_range(jcb, jce + 1);
  RAJA::RangeSegment i_range(icb, ice + 1);
  RAJA::kernel<DEFAULT_LOOP2X_ASYNC>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int jc, int ic) {
        int i = 2 * ic - 1, j = 2 * jc - 1;
        for (int c = 1; c <= lm_nc; c++)
          UcV(c, ic, jc, kc) =
              i1024 *
              (UfV(c, i - 3, j - 3, kf) - 9 * UfV(c, i - 3, j - 1, kf) -
               16 * UfV(c, i - 3, j, kf) - 9 * UfV(c, i - 3, j + 1, kf) +
               UfV(c, i - 3, j + 3, kf) +
               9 * (-UfV(c, i - 1, j - 3, kf) + 9 * UfV(c, i - 1, j - 1, kf) +
                    16 * UfV(c, i - 1, j, kf) + 9 * UfV(c, i - 1, j + 1, kf) -
                    UfV(c, i - 1, j + 3, kf)) +
               16 * (-UfV(c, i, j - 3, kf) + 9 * UfV(c, i, j - 1, kf) +
                     16 * UfV(c, i, j, kf) + 9 * UfV(c, i, j + 1, kf) -
                     UfV(c, i, j + 3, kf)) +
               9 * (-UfV(c, i + 1, j - 3, kf) + 9 * UfV(c, i + 1, j - 1, kf) +
                    16 * UfV(c, i + 1, j, kf) + 9 * UfV(c, i + 1, j + 1, kf) -
                    UfV(c, i + 1, j + 3, kf)) +
               UfV(c, i + 3, j - 3, kf) - 9 * UfV(c, i + 3, j - 1, kf) -
               16 * UfV(c, i + 3, j, kf) - 9 * UfV(c, i + 3, j + 1, kf) +
               UfV(c, i + 3, j + 3, kf));
      });
  // SYNC_STREAM;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::communicate_array_org(Sarray& u, bool allkplanes,
                                                  int kplane) {
  SW4_MARK_FUNCTION;
  //
  // General ghost point exchange at processor boundaries.
  //
  // Excplicit copy to buffers, not using fancy MPI-datatypes or sendrecv.
  //
  int kb = u.m_kb;
  int ke = u.m_ke;
  if (!allkplanes) ke = kb = kplane;
  const int ng = m_nghost;
  const int ni = (u.m_ie - u.m_ib + 1);
  const int nj = (u.m_je - u.m_jb + 1);
  const int nk = ke - kb + 1;
  float_sw4 *sbuf1, *sbuf2, *rbuf1, *rbuf2;

  MPI_Request req1, req2, req3, req4;
  MPI_Status status;
  int tag1 = 203, tag2 = 204;

  size_t npts1 = ng * nj * nk;
  size_t npts2 = ni * ng * nk;
  size_t nptsmax = max(npts1, npts2);
  // float_sw4* tmp = new float_sw4[4 * nptsmax * u.m_nc];
  float_sw4* tmp = m_mpi_buffer_space;
  if ((4 * nptsmax * u.m_nc) > m_mpi_buffer_size) {
    std::cerr << "MPI buffer size exceeded\n Aborting\n";
    abort();
  }
  sbuf1 = &tmp[0];
  rbuf1 = &tmp[nptsmax * u.m_nc];
  sbuf2 = &tmp[2 * nptsmax * u.m_nc];
  rbuf2 = &tmp[3 * nptsmax * u.m_nc];

  // i-direction communication
  MPI_Irecv(rbuf1, npts1 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag1,
            m_ew->m_cartesian_communicator, &req1);
  MPI_Irecv(rbuf2, npts1 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag2,
            m_ew->m_cartesian_communicator, &req2);
  if (m_ew->m_neighbor[0] != MPI_PROC_NULL)
    for (int c = 1; c <= u.m_nc; c++)
      for (int k = kb; k <= ke; k++)
        for (int j = u.m_jb; j <= u.m_je; j++)
          for (int i = u.m_ib + ng; i <= u.m_ib + 2 * ng - 1; i++) {
            size_t ind =
                i - (u.m_ib + ng) + ng * (j - u.m_jb) + ng * nj * (k - kb);
            sbuf1[ind + npts1 * (c - 1)] = u(c, i, j, k);
          }
  MPI_Isend(sbuf1, npts1 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag1,
            m_ew->m_cartesian_communicator, &req3);
  if (m_ew->m_neighbor[1] != MPI_PROC_NULL)
    for (int c = 1; c <= u.m_nc; c++)
      for (int k = kb; k <= ke; k++)
        for (int j = u.m_jb; j <= u.m_je; j++)
          for (int i = u.m_ie - 2 * ng + 1; i <= u.m_ie - ng; i++) {
            size_t ind = i - (u.m_ie - 2 * ng + 1) + ng * (j - u.m_jb) +
                         ng * nj * (k - kb);
            sbuf2[ind + npts1 * (c - 1)] = u(c, i, j, k);
          }
  MPI_Isend(sbuf2, npts1 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag2,
            m_ew->m_cartesian_communicator, &req4);
  MPI_Wait(&req1, &status);
  if (m_ew->m_neighbor[1] != MPI_PROC_NULL)
    for (int c = 1; c <= u.m_nc; c++)
      for (int k = kb; k <= ke; k++)
        for (int j = u.m_jb; j <= u.m_je; j++)
          for (int i = u.m_ie - ng + 1; i <= u.m_ie; i++) {
            size_t ind =
                i - (u.m_ie - ng + 1) + ng * (j - u.m_jb) + ng * nj * (k - kb);
            u(c, i, j, k) = rbuf1[ind + npts1 * (c - 1)];
          }
  MPI_Wait(&req2, &status);
  if (m_ew->m_neighbor[0] != MPI_PROC_NULL)
    for (int c = 1; c <= u.m_nc; c++)
      for (int k = kb; k <= ke; k++)
        for (int j = u.m_jb; j <= u.m_je; j++)
          for (int i = u.m_ib; i <= u.m_ib + ng - 1; i++) {
            size_t ind = i - u.m_ib + ng * (j - u.m_jb) + ng * nj * (k - kb);
            u(c, i, j, k) = rbuf2[ind + npts1 * (c - 1)];
          }

  MPI_Wait(&req3, &status);
  MPI_Wait(&req4, &status);

  // j-direction communication
  MPI_Irecv(rbuf1, npts2 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag1,
            m_ew->m_cartesian_communicator, &req1);
  MPI_Irecv(rbuf2, npts2 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag2,
            m_ew->m_cartesian_communicator, &req2);
  if (m_ew->m_neighbor[2] != MPI_PROC_NULL)
    for (int c = 1; c <= u.m_nc; c++)
      for (int k = kb; k <= ke; k++)
        for (int j = u.m_jb + ng; j <= u.m_jb + 2 * ng - 1; j++)
          for (int i = u.m_ib; i <= u.m_ie; i++) {
            size_t ind =
                i - u.m_ib + ni * (j - (u.m_jb + ng)) + ng * ni * (k - kb);
            sbuf1[ind + npts2 * (c - 1)] = u(c, i, j, k);
          }
  MPI_Isend(sbuf1, npts2 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag1,
            m_ew->m_cartesian_communicator, &req3);
  if (m_ew->m_neighbor[3] != MPI_PROC_NULL)
    for (int c = 1; c <= u.m_nc; c++)
      for (int k = kb; k <= ke; k++)
        for (int j = u.m_je - 2 * ng + 1; j <= u.m_je - ng; j++)
          for (int i = u.m_ib; i <= u.m_ie; i++) {
            size_t ind = i - u.m_ib + ni * (j - (u.m_je - 2 * ng + 1)) +
                         ng * ni * (k - kb);
            sbuf2[ind + npts2 * (c - 1)] = u(c, i, j, k);
          }
  MPI_Isend(sbuf2, npts2 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag2,
            m_ew->m_cartesian_communicator, &req4);
  MPI_Wait(&req1, &status);
  if (m_ew->m_neighbor[3] != MPI_PROC_NULL)
    for (int c = 1; c <= u.m_nc; c++)
      for (int k = kb; k <= ke; k++)
        for (int j = u.m_je - ng + 1; j <= u.m_je; j++)
          for (int i = u.m_ib; i <= u.m_ie; i++) {
            size_t ind =
                i - u.m_ib + ni * (j - (u.m_je - ng + 1)) + ng * ni * (k - kb);
            u(c, i, j, k) = rbuf1[ind + npts2 * (c - 1)];
          }
  MPI_Wait(&req2, &status);
  if (m_ew->m_neighbor[2] != MPI_PROC_NULL)
    for (int c = 1; c <= u.m_nc; c++)
      for (int k = kb; k <= ke; k++)
        for (int j = u.m_jb; j <= u.m_jb + ng - 1; j++)
          for (int i = u.m_ib; i <= u.m_ie; i++) {
            size_t ind = i - u.m_ib + ni * (j - u.m_jb) + ng * ni * (k - kb);
            u(c, i, j, k) = rbuf2[ind + npts2 * (c - 1)];
          }

  MPI_Wait(&req3, &status);
  MPI_Wait(&req4, &status);
  // delete[] tmp;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::communicate_array1d(float_sw4* u, int n, int dir,
                                                int ngh) {
  // WARNING :: This is using managed memory. Needs -M -gpu flag for safety
  SW4_MARK_FUNCTION;
  //
  // Communicate one dimensional array in i- or j-direction
  // Input:  u - The array
  //         n - Number of elements in array
  //         dir- Direction, 0 is i, 1 is j
  //         ngh- Number of points to communicate (ghost points)
  //
  // The first and last ngh points of u are updated.
  //
  MPI_Request req1, req2, req3, req4;
  MPI_Status status;
  int tag1 = 302, tag2 = 303;
  int no = 2 * dir;

  MPI_Irecv(&u[n - ngh], ngh, m_ew->m_mpifloat, m_ew->m_neighbor[1 + no], tag1,
            m_ew->m_cartesian_communicator, &req1);

  MPI_Irecv(&u[0], ngh, m_ew->m_mpifloat, m_ew->m_neighbor[no], tag2,
            m_ew->m_cartesian_communicator, &req2);

  MPI_Isend(&u[ngh], ngh, m_ew->m_mpifloat, m_ew->m_neighbor[no], tag1,
            m_ew->m_cartesian_communicator, &req3);

  MPI_Isend(&u[n - 2 * ngh], ngh, m_ew->m_mpifloat, m_ew->m_neighbor[1 + no],
            tag2, m_ew->m_cartesian_communicator, &req4);

  MPI_Wait(&req1, &status);
  MPI_Wait(&req2, &status);
}
//-----------------------------------------------------------------------
void CurvilinearInterface2::compute_icstresses_curv(
    Sarray& a_Up, Sarray& B, int kic, Sarray& a_metric, Sarray& a_mu,
    Sarray& a_lambda, float_sw4* a_str_x, float_sw4* a_str_y, float_sw4* sbop,
    char op) {
  const float_sw4 a1 = 2.0 / 3, a2 = -1.0 / 12;
  const bool upper = (kic == 1);
  const int k = kic;
  const int kl = upper ? 1 : -1;
  const int ifirst = a_Up.m_ib;
  const int jfirst = a_Up.m_jb;
#define str_x(i) a_str_x[(i - ifirst)]
#define str_y(j) a_str_y[(j - jfirst)]
  float_sw4 sgn = 1;
  if (op == '=') {
    B.set_value_async(0.0);
    sgn = 1;
  }
  if (op == '-') {
    sgn = -1;
  }

  SView& a_UpV = a_Up.getview();
  SView& BV = B.getview();
  SView& a_metricV = a_metric.getview();
  SView& a_muV = a_mu.getview();
  SView& a_lambdaV = a_lambda.getview();

  RAJA::RangeSegment j_range(B.m_jb + 2, B.m_je - 2 + 1);
  RAJA::RangeSegment i_range(B.m_ib + 2, B.m_ie - 2 + 1);
  RAJA::kernel<ODDIODDJ_EXEC_POL1_ASYNC>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
        // #pragma omp parallel for
        //   for (int j = B.m_jb + 2; j <= B.m_je - 2; j++)
        // #pragma omp simd
        //     for (int i = B.m_ib + 2; i <= B.m_ie - 2; i++) {
        float_sw4 uz, vz, wz;
        uz = vz = wz = 0;
        for (int m = 0; m <= 5; m++) {
          uz += sbop[m] * a_UpV(1, i, j, k + kl * (m - 1));
          vz += sbop[m] * a_UpV(2, i, j, k + kl * (m - 1));
          wz += sbop[m] * a_UpV(3, i, j, k + kl * (m - 1));
#ifdef CURVI_DEBUG
          std::cout << "UZ" << m << " " << uz << " " << sbop[m] << " "
                    << a_Up(1, i, j, k + kl * (m - 1)) << "\n";
#endif
        }
        uz *= kl;
        vz *= kl;
        wz *= kl;

        // Normal terms
        float_sw4 m2 = str_x(i) * a_metricV(2, i, j, k);
        float_sw4 m3 = str_y(j) * a_metricV(3, i, j, k);
        float_sw4 m4 = a_metricV(4, i, j, k);
        float_sw4 un = m2 * uz + m3 * vz + m4 * wz;
        float_sw4 mnrm = m2 * m2 + m3 * m3 + m4 * m4;
        float_sw4 B1, B2, B3;

        B1 = a_muV(i, j, k) * mnrm * uz +
             (a_muV(i, j, k) + a_lambdaV(i, j, k)) * m2 * un;
        B2 = a_muV(i, j, k) * mnrm * vz +
             (a_muV(i, j, k) + a_lambdaV(i, j, k)) * m3 * un;
        B3 = a_muV(i, j, k) * mnrm * wz +
             (a_muV(i, j, k) + a_lambdaV(i, j, k)) * m4 * un;

        // Tangential terms
        // p-derivatives
        float_sw4 up1 =
            str_x(i) * (a2 * (a_UpV(1, i + 2, j, k) - a_UpV(1, i - 2, j, k)) +
                        a1 * (a_UpV(1, i + 1, j, k) - a_UpV(1, i - 1, j, k)));
        float_sw4 up2 =
            str_x(i) * (a2 * (a_UpV(2, i + 2, j, k) - a_UpV(2, i - 2, j, k)) +
                        a1 * (a_UpV(2, i + 1, j, k) - a_UpV(2, i - 1, j, k)));
        float_sw4 up3 =
            str_x(i) * (a2 * (a_UpV(3, i + 2, j, k) - a_UpV(3, i - 2, j, k)) +
                        a1 * (a_UpV(3, i + 1, j, k) - a_UpV(3, i - 1, j, k)));
        B1 += a_metricV(1, i, j, k) *
              ((2 * a_muV(i, j, k) + a_lambdaV(i, j, k)) * m2 * up1 +
               a_muV(i, j, k) * (m3 * up2 + m4 * up3));
        B2 += a_metricV(1, i, j, k) *
              (a_lambdaV(i, j, k) * m3 * up1 + a_muV(i, j, k) * m2 * up2);
        B3 += a_metricV(1, i, j, k) *
              (a_lambdaV(i, j, k) * m4 * up1 + a_muV(i, j, k) * m2 * up3);

        // q-derivatives
        float_sw4 uq1 =
            str_y(j) * (a2 * (a_UpV(1, i, j + 2, k) - a_UpV(1, i, j - 2, k)) +
                        a1 * (a_UpV(1, i, j + 1, k) - a_UpV(1, i, j - 1, k)));
        float_sw4 uq2 =
            str_y(j) * (a2 * (a_UpV(2, i, j + 2, k) - a_UpV(2, i, j - 2, k)) +
                        a1 * (a_UpV(2, i, j + 1, k) - a_UpV(2, i, j - 1, k)));
        float_sw4 uq3 =
            str_y(j) * (a2 * (a_UpV(3, i, j + 2, k) - a_UpV(3, i, j - 2, k)) +
                        a1 * (a_UpV(3, i, j + 1, k) - a_UpV(3, i, j - 1, k)));
        B1 += a_metricV(1, i, j, k) *
              (a_lambdaV(i, j, k) * m2 * uq2 + a_muV(i, j, k) * m3 * uq1);
        B2 += a_metricV(1, i, j, k) *
              ((2 * a_muV(i, j, k) + a_lambdaV(i, j, k)) * m3 * uq2 +
               a_muV(i, j, k) * (m2 * uq1 + m4 * uq3));
        B3 += a_metricV(1, i, j, k) *
              (a_lambdaV(i, j, k) * m4 * uq2 + a_muV(i, j, k) * m3 * uq3);

        float_sw4 isgxy = 1.0 / (str_x(i) * str_y(j));
        B1 *= isgxy;
        B2 *= isgxy;
        B3 *= isgxy;
        //
        BV(1, i, j, k) += sgn * B1;
        BV(2, i, j, k) += sgn * B2;
        BV(3, i, j, k) += sgn * B3;
        // std::cout<<"STTRESS C"<<i<<j<<k<<" "<<B(1,i,j,k)<<"
        // "<<a_muV(i,j,k)<<"
        // "<<mnrm<<" "<<uz<<"\n";
      });
#ifdef PEEKS_GALORE
  SW4_PEEK;
  SYNC_DEVICE;
#endif
  SYNC_STREAM;
  // std::cout<<"END OF STRESS C CALL\n";
#undef str_x
#undef str_y
}
void CurvilinearInterface2::allocate_mpi_buffers() {
  int ni = m_ief - m_ibf + 1;
  int nj = m_jef - m_jbf + 1;
  int nk = m_kef - m_kbf + 1;
  int ng = m_nghost;

  size_t npts1 = ng * nj * nk;
  size_t npts2 = ni * ng * nk;

  m_mpi_buffer_size =
      max(npts1, npts2) * 4 * 4;  // Asumming u.m_nc is never more than 4

  m_mpi_buffer_space = SW4_NEW(Space::Pinned, float_sw4[m_mpi_buffer_size]);
}
//-----------------------------------------------------------------------
void CurvilinearInterface2::communicate_array(Sarray& u, bool allkplanes,
                                              int kplane) {
  SW4_MARK_FUNCTION;
  //
  // General ghost point exchange at processor boundaries.
  //
  // Excplicit copy to buffers, not using fancy MPI-datatypes or sendrecv.
  //
  int kb = u.m_kb;
  int ke = u.m_ke;
  if (!allkplanes) ke = kb = kplane;
  const int ng = m_nghost;
  const int ni = (u.m_ie - u.m_ib + 1);
  const int nj = (u.m_je - u.m_jb + 1);
  const int nk = ke - kb + 1;
  float_sw4 *sbuf1, *sbuf2, *rbuf1, *rbuf2;

  MPI_Request req1, req2, req3, req4;
  MPI_Status status;
  int tag1 = 203, tag2 = 204;

  size_t npts1 = ng * nj * nk;
  size_t npts2 = ni * ng * nk;
  size_t nptsmax = max(npts1, npts2);
  // float_sw4* tmp = new float_sw4[4 * nptsmax * u.m_nc];
  float_sw4* tmp = m_mpi_buffer_space;
  if ((4 * nptsmax * u.m_nc) > m_mpi_buffer_size) {
    std::cerr << "MPI buffer size exceeded\n Aborting\n";
    abort();
  }
  sbuf1 = &tmp[0];
  rbuf1 = &tmp[nptsmax * u.m_nc];
  sbuf2 = &tmp[2 * nptsmax * u.m_nc];
  rbuf2 = &tmp[3 * nptsmax * u.m_nc];

  int ib = u.m_ib;
  int ie = u.m_ie;
  int jb = u.m_jb;
  int lm_nc = u.m_nc;
  auto& uV = u.getview();
  RAJA::RangeSegment k_range(kb, ke + 1);

#ifdef PEEKS_GALORE
  SW4_PEEK;
  SYNC_DEVICE;
#endif

  // i-direction communication
  MPI_Irecv(rbuf1, npts1 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag1,
            m_ew->m_cartesian_communicator, &req1);
  MPI_Irecv(rbuf2, npts1 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag2,
            m_ew->m_cartesian_communicator, &req2);
  if (m_ew->m_neighbor[0] != MPI_PROC_NULL) {
    // for (int c = 1; c <= u.m_nc; c++)
    //   for (int k = kb; k <= ke; k++)
    //     for (int j = u.m_jb; j <= u.m_je; j++)
    //       for (int i = u.m_ib + ng; i <= u.m_ib + 2 * ng - 1; i++) {
    RAJA::RangeSegment j_range1(u.m_jb, u.m_je + 1);
    RAJA::RangeSegment i_range1(u.m_ib + ng, u.m_ib + 2 * ng - 1 + 1);
    RAJA::kernel<CA_POL>(RAJA::make_tuple(k_range, j_range1, i_range1),
                         [=] RAJA_DEVICE(int k, int j, int i) {
                           size_t ind = i - (ib + ng) + ng * (j - jb) +
                                        ng * nj * (k - kb);
                           for (int c = 1; c <= lm_nc; c++) {
                             sbuf1[ind + npts1 * (c - 1)] = uV(c, i, j, k);
                           }
                         });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }
  MPI_Isend(sbuf1, npts1 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag1,
            m_ew->m_cartesian_communicator, &req3);
  if (m_ew->m_neighbor[1] != MPI_PROC_NULL) {
    // for (int c = 1; c <= u.m_nc; c++)
    //   for (int k = kb; k <= ke; k++)
    //     for (int j = u.m_jb; j <= u.m_je; j++)
    //       for (int i = u.m_ie - 2 * ng + 1; i <= u.m_ie - ng; i++) {
    RAJA::RangeSegment j_range1(u.m_jb, u.m_je + 1);
    RAJA::RangeSegment i_range1(u.m_ie - 2 * ng + 1, u.m_ie - ng + 1);
    RAJA::kernel<CA_POL>(RAJA::make_tuple(k_range, j_range1, i_range1),
                         [=] RAJA_DEVICE(int k, int j, int i) {
                           size_t ind = i - (ie - 2 * ng + 1) + ng * (j - jb) +
                                        ng * nj * (k - kb);
                           for (int c = 1; c <= lm_nc; c++) {
                             sbuf2[ind + npts1 * (c - 1)] = uV(c, i, j, k);
                           }
                         });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }
  MPI_Isend(sbuf2, npts1 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag2,
            m_ew->m_cartesian_communicator, &req4);
  MPI_Wait(&req1, &status);
  if (m_ew->m_neighbor[1] != MPI_PROC_NULL) {
    // for (int c = 1; c <= u.m_nc; c++)
    //   for (int k = kb; k <= ke; k++)
    //     for (int j = u.m_jb; j <= u.m_je; j++)
    //       for (int i = u.m_ie - ng + 1; i <= u.m_ie; i++) {
    RAJA::RangeSegment j_range1(u.m_jb, u.m_je + 1);
    RAJA::RangeSegment i_range1(u.m_ie - ng + 1, u.m_ie + 1);
    RAJA::kernel<CA_POL>(RAJA::make_tuple(k_range, j_range1, i_range1),
                         [=] RAJA_DEVICE(int k, int j, int i) {
                           size_t ind = i - (ie - ng + 1) + ng * (j - jb) +
                                        ng * nj * (k - kb);
                           for (int c = 1; c <= lm_nc; c++) {
                             uV(c, i, j, k) = rbuf1[ind + npts1 * (c - 1)];
                           }
                         });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }
  MPI_Wait(&req2, &status);
  if (m_ew->m_neighbor[0] != MPI_PROC_NULL) {
    // for (int c = 1; c <= u.m_nc; c++)
    //   for (int k = kb; k <= ke; k++)
    //     for (int j = u.m_jb; j <= u.m_je; j++)
    //       for (int i = u.m_ib; i <= u.m_ib + ng - 1; i++) {
    RAJA::RangeSegment j_range1(u.m_jb, u.m_je + 1);
    RAJA::RangeSegment i_range1(u.m_ib, u.m_ib + ng - 1 + 1);
    RAJA::kernel<CA_POL>(RAJA::make_tuple(k_range, j_range1, i_range1),
                         [=] RAJA_DEVICE(int k, int j, int i) {
                           size_t ind =
                               i - ib + ng * (j - jb) + ng * nj * (k - kb);
                           for (int c = 1; c <= lm_nc; c++) {
                             uV(c, i, j, k) = rbuf2[ind + npts1 * (c - 1)];
                           }
                         });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }

  MPI_Wait(&req3, &status);
  MPI_Wait(&req4, &status);

  // j-direction communication
  MPI_Irecv(rbuf1, npts2 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag1,
            m_ew->m_cartesian_communicator, &req1);
  MPI_Irecv(rbuf2, npts2 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag2,
            m_ew->m_cartesian_communicator, &req2);

  if (m_ew->m_neighbor[2] != MPI_PROC_NULL) {
    // for (int c = 1; c <= u.m_nc; c++)
    //   for (int k = kb; k <= ke; k++)
    //     for (int j = u.m_jb + ng; j <= u.m_jb + 2 * ng - 1; j++)
    //       for (int i = u.m_ib; i <= u.m_ie; i++) {

    RAJA::RangeSegment j_range1(u.m_jb + ng, u.m_jb + 2 * ng - 1 + 1);
    RAJA::RangeSegment i_range1(u.m_ib, u.m_ie + 1);
    RAJA::kernel<CA_POL>(RAJA::make_tuple(k_range, j_range1, i_range1),
                         [=] RAJA_DEVICE(int k, int j, int i) {
                           size_t ind = i - ib + ni * (j - (jb + ng)) +
                                        ng * ni * (k - kb);
                           for (int c = 1; c <= lm_nc; c++) {
                             sbuf1[ind + npts2 * (c - 1)] = uV(c, i, j, k);
                           }
                         });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }

  MPI_Isend(sbuf1, npts2 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag1,
            m_ew->m_cartesian_communicator, &req3);

  int je = u.m_je;
  if (m_ew->m_neighbor[3] != MPI_PROC_NULL) {
    // for (int c = 1; c <= u.m_nc; c++)
    //   for (int k = kb; k <= ke; k++)
    //     for (int j = u.m_je - 2 * ng + 1; j <= u.m_je - ng; j++)
    //       for (int i = u.m_ib; i <= u.m_ie; i++) {
    RAJA::RangeSegment j_range2(u.m_je - 2 * ng + 1, u.m_je - ng + 1);
    RAJA::RangeSegment i_range2(u.m_ib, u.m_ie + 1);
    RAJA::kernel<CA_POL>(RAJA::make_tuple(k_range, j_range2, i_range2),
                         [=] RAJA_DEVICE(int k, int j, int i) {
                           size_t ind = i - ib + ni * (j - (je - 2 * ng + 1)) +
                                        ng * ni * (k - kb);
                           for (int c = 1; c <= lm_nc; c++) {
                             sbuf2[ind + npts2 * (c - 1)] = uV(c, i, j, k);
                           }
                         });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }
  MPI_Isend(sbuf2, npts2 * u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag2,
            m_ew->m_cartesian_communicator, &req4);
  MPI_Wait(&req1, &status);
  if (m_ew->m_neighbor[3] != MPI_PROC_NULL) {
    // for (int c = 1; c <= u.m_nc; c++){
    //   for (int k = kb; k <= ke; k++)
    //     for (int j = u.m_je - ng + 1; j <= u.m_je; j++)
    //       for (int i = u.m_ib; i <= u.m_ie; i++) {
    RAJA::RangeSegment j_range3(u.m_je - ng + 1, u.m_je + 1);
    RAJA::RangeSegment i_range3(u.m_ib, u.m_ie + 1);
    RAJA::kernel<CA_POL>(RAJA::make_tuple(k_range, j_range3, i_range3),
                         [=] RAJA_DEVICE(int k, int j, int i) {
                           size_t ind = i - ib + ni * (j - (je - ng + 1)) +
                                        ng * ni * (k - kb);
                           for (int c = 1; c <= lm_nc; c++) {
                             uV(c, i, j, k) = rbuf1[ind + npts2 * (c - 1)];
                           }
                         });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }
  MPI_Wait(&req2, &status);
  if (m_ew->m_neighbor[2] != MPI_PROC_NULL) {
    // for (int c = 1; c <= u.m_nc; c++)
    //   for (int k = kb; k <= ke; k++)
    //     for (int j = u.m_jb; j <= u.m_jb + ng - 1; j++)
    //       for (int i = u.m_ib; i <= u.m_ie; i++) {
    RAJA::RangeSegment j_range4(u.m_jb, u.m_jb + ng - 1 + 1);
    RAJA::RangeSegment i_range4(u.m_ib, u.m_ie + 1);
    RAJA::kernel<CA_POL>(RAJA::make_tuple(k_range, j_range4, i_range4),
                         [=] RAJA_DEVICE(int k, int j, int i) {
                           size_t ind =
                               i - ib + ni * (j - jb) + ng * ni * (k - kb);
                           for (int c = 1; c <= lm_nc; c++) {
                             uV(c, i, j, k) = rbuf2[ind + npts2 * (c - 1)];
                           }
                         });
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }

  MPI_Wait(&req3, &status);
  MPI_Wait(&req4, &status);
  // delete[] tmp;
}
