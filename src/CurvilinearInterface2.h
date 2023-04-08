#ifndef SW4_CURVILINEARSW4_TYPEERFACE2
#define SW4_CURVILINEARSW4_TYPEERFACE2

#include <vector>

#include "TestEcons.h"
#include "TestPointSource.h"
#include "TestTwilight.h"
#ifdef USE_MAGMA
#include "magma_dbatched.h"
#include "magma_v2.h"
#endif

class Sarray;
class EW;

class CurvilinearSw4_Typeerface2 {
  EW* m_ew;
  TestTwilight* m_tw;
  TestEcons* m_etest;
  TestPointSource* m_psource;

  sw4_type m_nghost, m_ib, m_ie, m_jb, m_je, m_ibf, m_ief, m_jbf, m_jef, m_nkf;
  sw4_type m_kb, m_ke, m_kbf, m_kef;
  sw4_type m_gc, m_gf;
  bool m_isbndry[4];  // side is physical boundary

  float_sw4 m_reltol, m_abstol;
  sw4_type m_maxit;

  float_sw4 *m_strx_c, *m_strx_f, *m_stry_c, *m_stry_f;
  Sarray m_Mass_block, m_rho_c, m_rho_f, m_mu_c, m_mu_f, m_lambda_c, m_lambda_f;
  Sarray m_met_c, m_met_f, m_jac_c, m_jac_f, m_x_c, m_x_f, m_y_c, m_y_f, m_z_c,
      m_z_f;
  bool m_use_attenuation;
  sw4_type m_number_mechanisms;
  std::vector<Sarray> m_muve_f, m_lambdave_f, m_muve_c, m_lambdave_c;

  float_sw4* m_mass_block;
  sw4_type* m_ipiv_block;
  float_sw4* m_mpi_buffer_space;
  size_t m_mpi_buffer_size;

#ifdef USE_MAGMA
  float_sw4* m_mass_block_gpu;
  sw4_type* m_ipiv_block_gpu;

  float_sw4** dA_array;
  magma_sw4_type_t** piv_array;
  magma_queue_t queue;
  float_sw4** dB_array;
  float_sw4* x;
  std::vector<sw4_type> subbatchsize, subbatchoffset;
#endif

#ifdef USE_DIRECT_INVERSE
  float_sw4* m_mass_block_gpu;
  float_sw4* x;
#endif

#if defined(ENABLE_GPU)
  float_sw4 *m_sbop, *m_acof, *m_bop, *m_bope, *m_ghcof;
  float_sw4 *m_acof_no_gp, *m_ghcof_no_gp, *m_sbop_no_gp;
#else
  float_sw4 m_acof[384], m_bope[48], m_ghcof[6], m_acof_no_gp[384],
      m_ghcof_no_gp[6];
  float_sw4 m_sbop[6], m_sbop_no_gp[6], m_bop[24];
#endif

  void sw4_typeerface_block(Sarray& matrix);

  void compute_icstresses_curv_host(Sarray& a_Up, Sarray& B, sw4_type kic,
                                    Sarray& a_metric, Sarray& a_mu,
                                    Sarray& a_lambda, float_sw4* a_str_x,
                                    float_sw4* a_str_y, float_sw4* sbop,
                                    char op);

  void mat_icstresses_curv(sw4_type ib, sw4_type jb, Sarray& a_mat, sw4_type kic,
                           Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
                           float_sw4* a_str_x, float_sw4* a_str_y,
                           float_sw4* sbop);

  void matrix_Lu(sw4_type ib, sw4_type jb, Sarray& a_mat, Sarray& met, Sarray& jac,
                 Sarray& mu, Sarray& la, float_sw4* a_str_x, float_sw4* a_str_y,
                 float_sw4 ghcof);

  void restprol2D(Sarray& Uc, Sarray& alpha, sw4_type kc, sw4_type kf);

  void copy_str(float_sw4* dest, float_sw4* src, sw4_type offset, sw4_type n, sw4_type nsw);
  void communicate_array1d(float_sw4* u, sw4_type n, sw4_type dir, sw4_type ngh);

  void init_arrays_att();
  void allocate_mpi_buffers();

 public:
  CurvilinearSw4_Typeerface2(sw4_type a_gc, EW* a_ew);
  CurvilinearSw4_Typeerface2() {}
  ~CurvilinearSw4_Typeerface2();
  void init_arrays(std::vector<float_sw4*>& a_strx,
                   std::vector<float_sw4*>& a_stry);
  //   void test1( EW* a_ew, sw4_type gc, std::vector<Sarray>& a_U );
  //   void test2( EW* a_ew, sw4_type gc, std::vector<Sarray>& a_U );

  void impose_ic(std::vector<Sarray>& a_U, float_sw4 t,
                 std::vector<Sarray*>& a_AlphaVE);
  void sw4_typeerface_lhs(Sarray& lhs, Sarray& uc);
  void compute_icstresses_curv(Sarray& a_Up, Sarray& B, sw4_type kic,
                               Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
                               float_sw4* a_str_x, float_sw4* a_str_y,
                               float_sw4* sbop, char op);
  void lhs_Lu(Sarray& a_U, Sarray& a_lhs, Sarray& metric, Sarray& jac,
              Sarray& mu, Sarray& lambda, float_sw4* a_str_x,
              float_sw4* a_str_y, float_sw4 ghcof);
  void bnd_zero(Sarray& u, sw4_type npts);
  void bnd_zero_host(Sarray& u, sw4_type npts);
  void injection(Sarray& u_f, Sarray& u_c);
  void sw4_typeerface_rhs(Sarray& rhs, Sarray& uc, Sarray& uf,
                     std::vector<Sarray>& Alpha_c,
                     std::vector<Sarray>& Alpha_f);
  void prolongate2D(Sarray& Uc, Sarray& Uf, sw4_type kc, sw4_type kf);
  void restrict2D(Sarray& Uc, Sarray& Uf, sw4_type kc, sw4_type kf);
  void lhs_icstresses_curv(Sarray& a_Up, Sarray& a_lhs, sw4_type kic,
                           Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
                           float_sw4* a_str_x, float_sw4* a_str_y,
                           float_sw4* sbop);
  void communicate_array(Sarray& u, bool allkplanes = true, sw4_type kplane = 0);
  void communicate_array_org(Sarray& u, bool allkplanes = true, sw4_type kplane = 0);
};

#endif
