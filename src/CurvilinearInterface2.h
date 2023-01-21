#ifndef SW4_CURVILINEARINTERFACE2
#define SW4_CURVILINEARINTERFACE2

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

class CurvilinearInterface2 {
  EW* m_ew;
  TestTwilight* m_tw;
  TestEcons* m_etest;
  TestPointSource* m_psource;

  int m_nghost, m_ib, m_ie, m_jb, m_je, m_ibf, m_ief, m_jbf, m_jef, m_nkf;
  int m_kb, m_ke, m_kbf, m_kef;
  int m_gc, m_gf;
  bool m_isbndry[4];  // side is physical boundary

  float_sw4 m_reltol, m_abstol;
  int m_maxit;

  float_sw4 *m_strx_c, *m_strx_f, *m_stry_c, *m_stry_f;
  Sarray m_Mass_block, m_rho_c, m_rho_f, m_mu_c, m_mu_f, m_lambda_c, m_lambda_f;
  Sarray m_met_c, m_met_f, m_jac_c, m_jac_f, m_x_c, m_x_f, m_y_c, m_y_f, m_z_c,
      m_z_f;
  bool m_use_attenuation;
  int m_number_mechanisms;
  std::vector<Sarray> m_muve_f, m_lambdave_f, m_muve_c, m_lambdave_c;

  float_sw4* m_mass_block;
  int* m_ipiv_block;
  bool m_memory_is_allocated;
  
  float_sw4* m_mpi_buffer_space;
  size_t m_mpi_buffer_size;

#ifdef USE_MAGMA
  float_sw4* m_mass_block_gpu;
  int* m_ipiv_block_gpu;

  float_sw4** dA_array;
  magma_int_t** piv_array;
  magma_queue_t queue;
  float_sw4** dB_array;
  float_sw4* x;
  std::vector<int> subbatchsize, subbatchoffset;
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

  void interface_block(Sarray& matrix);

  void compute_icstresses_curv_host(Sarray& a_Up, Sarray& B, int kic,
                                    Sarray& a_metric, Sarray& a_mu,
                                    Sarray& a_lambda, float_sw4* a_str_x,
                                    float_sw4* a_str_y, float_sw4* sbop,
                                    char op);

  void mat_icstresses_curv(int ib, int jb, Sarray& a_mat, int kic,
                           Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
                           float_sw4* a_str_x, float_sw4* a_str_y,
                           float_sw4* sbop);

  void matrix_Lu(int ib, int jb, Sarray& a_mat, Sarray& met, Sarray& jac,
                 Sarray& mu, Sarray& la, float_sw4* a_str_x, float_sw4* a_str_y,
                 float_sw4 ghcof);

  void restprol2D(Sarray& Uc, Sarray& alpha, int kc, int kf);

  void copy_str(float_sw4* dest, float_sw4* src, int offset, int n, int nsw);
  void communicate_array1d(float_sw4* u, int n, int dir, int ngh);

  void init_arrays_att();
  void allocate_mpi_buffers();

 public:
  CurvilinearInterface2(int a_gc, EW* a_ew);
  CurvilinearInterface2() {}
  ~CurvilinearInterface2();
  void init_arrays(std::vector<float_sw4*>& a_strx,
                   std::vector<float_sw4*>& a_stry, std::vector<Sarray>& a_rho,
                   std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda);
  //   void test1( EW* a_ew, int gc, std::vector<Sarray>& a_U );
  //   void test2( EW* a_ew, int gc, std::vector<Sarray>& a_U );

  void impose_ic(std::vector<Sarray>& a_U, float_sw4 t,
                 std::vector<Sarray*>& a_AlphaVE);
  void interface_lhs(Sarray& lhs, Sarray& uc);
  void compute_icstresses_curv(Sarray& a_Up, Sarray& B, int kic,
                               Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
                               float_sw4* a_str_x, float_sw4* a_str_y,
                               float_sw4* sbop, char op);
  void lhs_Lu(Sarray& a_U, Sarray& a_lhs, Sarray& metric, Sarray& jac,
              Sarray& mu, Sarray& lambda, float_sw4* a_str_x,
              float_sw4* a_str_y, float_sw4 ghcof);
  void bnd_zero(Sarray& u, int npts);
  void bnd_zero_host(Sarray& u, int npts);
  void injection(Sarray& u_f, Sarray& u_c);
  void interface_rhs(Sarray& rhs, Sarray& uc, Sarray& uf,
                     std::vector<Sarray>& Alpha_c,
                     std::vector<Sarray>& Alpha_f);
  void prolongate2D(Sarray& Uc, Sarray& Uf, int kc, int kf);
  void restrict2D(Sarray& Uc, Sarray& Uf, int kc, int kf);
  void lhs_icstresses_curv(Sarray& a_Up, Sarray& a_lhs, int kic,
                           Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
                           float_sw4* a_str_x, float_sw4* a_str_y,
                           float_sw4* sbop);
  void communicate_array(Sarray& u, bool allkplanes = true, int kplane = 0);
  void communicate_array_org(Sarray& u, bool allkplanes = true, int kplane = 0);
};

#endif
