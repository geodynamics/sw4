#ifndef SW4_CURVILINEARINTERFACE
#define SW4_CURVILINEARINTERFACE

#include "EW.h"
#include "TestGrid.h"
#include "TestTwilight.h"
// Coefficients: Rop, Sb, ghcof, P
#include "Farray.h"
//class Farray;

class CurvilinearInterface
{
 private:
   //   Farray Rop(-4,2), P(-1,2), ghcof(1,6), Sb(0,4), acof(1,6,1,8,1,8);
   //   Farray bof(1,4,1,6), acof_no_gp(1,6,1,8,1,8), sbop_no_gp(0,5);
   //   Farray ux_cof(-2,2);
   Farray Rop, P, ghcof, Sb, acof;
   Farray bof, acof_no_gp, sbop_no_gp;
   Farray ux_cof;
   Farray Mass_block;

   Sarray rho_c, rho_f, Jacobian_c, Jacobian_f;
   Sarray mu_c, lambda_c, mu_f, lambda_f;
   Sarray XI13_c, XI23_c, XI33_c;
   Sarray XI13_f, XI23_f, XI33_f;
   Sarray x_c, y_c, z_c;
   Sarray x_f, y_f, z_f;
   int* IPIV_block;
   int gc, gf, nghost;
   float_sw4 hc, hf;
   //   std::vector<float_sw4> m_scale_factors;

   PackArgs a;
   TestGrid* m_test_grid;
   TestTwilight* m_tw;
   EW* m_ew;

   void define_coeffs();
   void varcoeffs4( Farray& acof, Farray& ghcof, Farray& Sb );
   void dx_46( Farray &bop );
   void varcoeff_noghost();

   void convert_metric( Sarray& jac, Sarray& met, Sarray& XI13, Sarray& XI23, 
                        Sarray& XI33, Sarray& Jacobian, float_sw4 h, int n1, int n2, int n3 );

   void interface_block(Farray &Rop, Farray &ghcof, Farray &Sb, Sarray &rho_c,
                        Sarray &lambda_c, Sarray &rho_f, Sarray &Jacobian_c,
                        Sarray &mu_c, Farray &P, Sarray &XI13_c, Sarray &XI23_c,
                        Sarray &XI33_c, Farray &Mass_block, PackArgs &a);

   void interface_rhs(Farray &Vass, Farray &lh_c, Farray &lh_f, Sarray &Jacobian_c,
                      Sarray &Jacobian_f, Sarray &mu_c, Sarray &mu_f,
                      Sarray &lambda_c, Sarray &lambda_f, Sarray &rho_c,
                      Sarray &rho_f, Sarray &XI13_c, Sarray &XI23_c,
                      Sarray &XI33_c, Sarray &XI13_f, Sarray &XI23_f,
                      Sarray &XI33_f, Farray &P, Farray &Sb, Farray &Rop,
                      Farray &sbop_no_gp, Farray &acof_no_gp, Sarray &u_c,
                      Sarray &u_f, Farray &Mass_f1, Farray &ux_cof, Farray &ghcof,
                      Farray &acof, Farray &bof, PackArgs &a );

   void interface_lhs(Farray &LHS, Sarray &Jacobian_c,
                      Sarray &mu_c, Sarray &lambda_c, Sarray &rho_c,
                      Sarray &rho_f, Sarray &XI13_c, Sarray &XI23_c,
                      Sarray &XI33_c, Farray &P, Farray &Sb, Farray &Rop,
                      Sarray &u_c, Farray &ghcof, PackArgs &a );

   void injection(Sarray &u_f, Sarray &u_c, Farray &P, PackArgs &a );
   void scale_rho( Sarray& rho, Sarray& Jacobian );
public:
   CurvilinearInterface( int a_gc, EW* a_ew  );
   void init_arrays( std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda, 
                     std::vector<Sarray>& a_rho, std::vector<Sarray>& a_metric, 
                     std::vector<Sarray>& a_jac );
   void impose_ic( std::vector<Sarray>& a_U, float_sw4 t );
};

#endif
