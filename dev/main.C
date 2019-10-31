#include <math.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "Farray.h"
#include "Sarray.h"
#include "cfunctions.h"
extern "C" {
// void interpolation_restriction(double *P, double *Rop,double *RPop);
// void central_difference(double *ux_cof,double *uxx_cof);
// void varcoeff_noghost(double *acof_no_gp, double *ghcof_no_gp, double
// *sbop_no_gp ); void varcoeffs4( double *acof, double *ghcof, double *Sb );
// void dx_46(double *);
// void equation_cof(double *, double *, double *,double*,double*,
// 		    double*,double*,double*,double*,double*,
// 		    double*,double*,double*,double*,double*,
// 		    double*,double*,double*,double*,double*);
// void kappa(double*,double*,double*,double*,double*,double*,double*,double*);
// void exact_solution(double *, double *, double *, double *, double *, double
// *, double *, int *); void interface_block(double *, double *, double*, double
// *, double *, double *, double *, 		       double *, double *,
// double
// *, double *, double *, double *);
void dgetrf_wrap(int *, int *, double *, int *, int *, int *);
void dgetrs_wrap(int, int, double *, int, int *, double *, int, int *);
// void injection(double*, double*, double*,int);
// void interface_rhs(double*,double*,double*,double*,double*,double*,double*,
// 		     double*,double*,double*,double*,double*,double*,double*,
// 		     double*,double*,double*,double*,double*,double*,double*,
// 		     double*,double*,double*,double*,double*,double*,double*,double*,int);
// void interface_lhs(double*,double*,double*,double*,double*,double*,double*,
// 		     double*,double*,double*,double*,double*,double*,double*,
// 		     double*,double*,double*,double*,double*,double*,double*,
// 		     double*,double*,double*,double*,double*,double*,double*,double*,int);
// void update_interior(double*,double*,double*,double*,double*,
// 		       double*,double*,double*,double*,double*,
// 		       double*,double*,double*,double*,double*,
// 		       double*,double*,double*,double*,double*);
// void forcing(double*,double*,double*,double*,double*,double*,double*);
// void forcing_tt(double*,double*,double*,double,double*,double*,double*);
// void
// update_gp(double*,double*,double*,double*,double*,double*,double*,double*,double*,int);
// void update_traction(double*,double*,double*,double*,double*,
// 		       double*,double*,double*,double*,double*,
// 		       double*,double*,double*,double*,int);
// void update_dirichlet_bc(double*,double*,double*,double*,double*,int);
void print_uf(double *);
void print_f(double *);
}
int Farray::count = 0;
int main(int argc, char *argv[]) {
  const int dim = 3;
  const float_sw4 tn = 1.0;
  const float_sw4 pi = M_PI;
  const float_sw4 l1 = 2.0 * pi, l2 = 2.0 * pi,
                  l3 = 2.0 * pi;  // space interval
  const float_sw4 int_pos =
      0.50;  // This is the position in r3 where the interface is located.
  const int nrg = 5;
  const int n1_c = 25;
  const int n2_c = 25;
  const float_sw4 h1phy_c = l1 / (n1_c - 1),
                  h1phy_f = h1phy_c * 0.50;  // mesh size in physical space, x
  const float_sw4 h2phy_c = l2 / (n2_c - 1),
                  h2phy_f = h2phy_c * 0.50;  // mesh size in physical space, y
  const int n3_c = std::ceil(int_pos * l3 / h1phy_c) +
                   1;  //  number of grid points in direction-3
  const int n1_f = n1_c * 2 - 1, n2_f = n2_c * 2 - 1,
            n3_f = std::ceil((1.0 - int_pos) * l3 / (h1phy_f)) + 1;
  const float_sw4 h1_c = 1.0 / (n1_c - 1), h2_c = 1.0 / (n2_c - 1),
                  h3_c = 1.0 / (n3_c - 1);
  const float_sw4 h1_f = 1.0 / (n1_f - 1), h2_f = 1.0 / (n2_f - 1),
                  h3_f = 1.0 / (n3_f - 1);
  const float_sw4 amp = 0.20, peak = 0.040;

  PackArgs a;
  a.dim = dim;
  a.tn = tn;
  a.pi = pi;

  a.l1 = l1;
  a.l2 = l2;
  a.l3 = l3;

  a.int_pos = int_pos;

  a.nrg = nrg;

  a.n1_c = n1_c;
  a.n2_c = n2_c;
  a.n3_c = n3_c;

  a.h1phy_c = h1phy_c;
  a.h2phy_c = h2phy_c;

  a.h1phy_f = h1phy_f;
  a.h2phy_f = h2phy_f;

  a.n1_f = n1_f;
  a.n2_f = n2_f;
  a.n3_f = n3_f;

  a.h1_c = h1_c;
  a.h2_c = h2_c;
  a.h3_c = h3_c;

  a.h1_f = h1_f;
  a.h2_f = h2_f;
  a.h3_f = h3_f;

  a.amp = amp;
  a.peak = peak;

  std::chrono::high_resolution_clock::time_point t1, t2;

  Farray RPop(-3, 3);
  Farray P(-1, 2);
  Farray Rop(-4, 2);
  Farray ux_cof(-2, 2);
  Farray uxx_cof(-2, 2);
  Farray acof_no_gp(1, 6, 1, 8, 1, 8);
  Farray ghcof_no_gp(1, 6);
  Farray sbop_no_gp(0, 5);
  Farray acof(1, 6, 1, 8, 1, 8);
  Farray ghcof(1, 6);
  Farray Sb(0, 4);
  Farray bof(1, 4, 1, 6);

  // allocate memory for Physical domain X
  Farray Xgrid_c_1(1 - nrg, n1_c + nrg);
  Farray Xgrid_c_2(1 - nrg, n2_c + nrg);
  Farray Xgrid_f_1(1 - nrg, n1_f + nrg);
  Farray Xgrid_f_2(1 - nrg, n2_f + nrg);
  Farray Xgrid_c_3(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg,
                   n3_c + nrg);
  Farray Xgrid_f_3(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg,
                   n3_f + nrg);

  // allocate memory for forcing function and spatial discretization

  Farray force_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg,
                 1, dim);
  Farray force_tt_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg,
                    n3_c + nrg, 1, dim);
  Farray lh_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg, 1,
              dim);
  Farray force_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg,
                 1, dim);
  Farray force_tt_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg,
                    n3_f + nrg, 1, dim);
  Farray lh_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg, 1,
              dim);

  // allocate memory for metric derivatives
  Sarray XI13_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg);
  Sarray XI23_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg);
  Sarray XI33_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg);
  Sarray Jacobian_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg,
                    n3_c + nrg);
  Sarray XI13_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg);
  Sarray XI23_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg);
  Sarray XI33_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg);
  Sarray Jacobian_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg,
                    n3_f + nrg);

  // allocate memory for material parameters
  Sarray rho_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg);
  Sarray mu_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg);
  Sarray lambda_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg,
                  n3_c + nrg);
  Sarray rho_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg);
  Sarray mu_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg);
  Sarray lambda_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg,
                  n3_f + nrg);
  Garray<2> g = {1 - nrg, n1_c + nrg, 3, 4};

  // allocate memory for solutions
  //Farray u_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg, 1,
  //          dim, 1, 4);
  //Farray u_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg, 1,
  //          dim, 1, 4);
  // Split arrays as in SW4 
  Farray Up_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg, 1,dim);
  Farray Um_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg, 1,dim);
  Farray U_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg, 1,dim);
  Farray Ux_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg, 1,dim);

  Farray Up_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg, 1,
             dim);
  Farray U_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg, 1,
             dim);
  Farray Um_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg, 1,
             dim);
  Farray Ux_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg, 1,
             dim);
  // Temps not in original code
  Farray u_ct(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1 - nrg, n3_c + nrg, 1,
              dim);
  Farray u_ft(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, 1 - nrg, n3_f + nrg, 1,
              dim);

  // allocate memory for error
  Sarray err_f(1, n1_f, 1, n2_f, 1, n3_f);
  Sarray err_c(1, n1_c, 1, n2_c, 1, n3_c);

  // allocate memory for interface linear system
  Farray Vass(1, n1_c * n2_c * 3);
  Farray LHS(1, n1_c * n2_c * 3);
  Farray residual(1, n1_c * n2_c * 3);
  Farray Mass_f1(-2, n1_f + 3, -2, n2_f + 3);
  Farray Mass_block(1, 3, 1, 3, 1, n1_c, 1, n2_c);

  int *IPIV_block = new int[3 * n1_c * n2_c];

  // traction B.C. on the top
  float_sw4 traction_rhs[3] = {0, 0, 0};
  Sarray traction_data(1, n1_f, 1, n2_f, 1, dim);

  // Difference stencils and interpolations
  interpolation_restriction(P, Rop, RPop);
  central_difference(ux_cof, uxx_cof);
  varcoeff_noghost(acof_no_gp, ghcof_no_gp, sbop_no_gp);
  varcoeffs4(acof, ghcof, Sb);
  dx_46(bof);

  RPop.print();
  P.print();
  Rop.print();
  bof.print();

  // metric derivatives and coefficients of the equation
  // equation_cof(Xgrid_c_1.get(),Xgrid_c_2.get(),Xgrid_c_3.get(),Xgrid_f_1.get(),Xgrid_f_2.get(),Xgrid_f_3.get(),
  // 	       XI13_c.get(),XI23_c.get(),XI33_c.get(),Jacobian_c.get(),rho_c.get(),rho_f.get(),
  // 	       XI13_f.get(),XI23_f.get(),XI33_f.get(),Jacobian_f.get(),
  // 	       mu_c.get(),mu_f.get(),lambda_c.get(),lambda_f.get());
  equation_cof(Xgrid_c_1, Xgrid_c_2, Xgrid_c_3, Xgrid_f_1, Xgrid_f_2, Xgrid_f_3,
               XI13_c, XI23_c, XI33_c, Jacobian_c, rho_c, rho_f, XI13_f, XI23_f,
               XI33_f, Jacobian_f, mu_c, mu_f, lambda_c, lambda_f, a);

  // estimate time step
  float_sw4 tv;
  // kappa(mu_f.get(),lambda_f.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),Jacobian_f.get(),rho_f.get(),&tv);
  kappa(mu_f, lambda_f, XI13_f, XI23_f, XI33_f, Jacobian_f, rho_f, tv, a);
  float_sw4 cfl = sqrt(3.0 / 1.40 / tv);
  float_sw4 dt0 = std::min(h1_f, std::min(h2_f, h3_f)) * cfl * 0.9;
  int nt = std::ceil(tn / dt0) + 1;

  if (argc >= 2) nt = 2;  // Shorter run for debuggin
  float_sw4 dt = tn / nt;
  std::cout << " TV IS " << tv << " NT is " << nt << " should be 70 \n";
  // exact solution at t = -dt, 0, dt
  int flag = 0;
  float_sw4 zero = 0.0;
  float_sw4 mdt = -dt;
  for (int k = 1 - nrg; k <= n3_c + nrg; k++) {
    for (int i = 1 - nrg; i <= n2_c + nrg; i++) {
      for (int j = 1 - nrg; j <= n1_c + nrg; j++) {
        exact_solution(Xgrid_c_1(j), Xgrid_c_2(i), Xgrid_c_3(j, i, k), mdt,
                       Um_c(j, i, k, 1), Um_c(j, i, k, 2),
                       Um_c(j, i, k, 3), flag);
        exact_solution(Xgrid_c_1(j), Xgrid_c_2(i), Xgrid_c_3(j, i, k), zero,
                       U_c(j, i, k, 1), U_c(j, i, k, 2),
                       U_c(j, i, k, 3), flag);
        exact_solution(Xgrid_c_1(j), Xgrid_c_2(i), Xgrid_c_3(j, i, k), dt,
                       Up_c(j, i, k, 1), Up_c(j, i, k, 2),
                       Up_c(j, i, k, 3), flag);
      }
    }
  }
  flag = 1;
  for (int k = 1 - nrg; k <= n3_f + nrg; k++) {
    for (int i = 1 - nrg; i <= n2_f + nrg; i++) {
      for (int j = 1 - nrg; j <= n1_f + nrg; j++) {
        exact_solution(Xgrid_f_1(j), Xgrid_f_2(i), Xgrid_f_3(j, i, k), mdt,
                       Um_f(j, i, k, 1), Um_f(j, i, k, 2),
                       Um_f(j, i, k, 3), flag);
        exact_solution(Xgrid_f_1(j), Xgrid_f_2(i), Xgrid_f_3(j, i, k), zero,
                       U_f(j, i, k, 1), U_f(j, i, k, 2),
                       U_f(j, i, k, 3), flag);
        exact_solution(Xgrid_f_1(j), Xgrid_f_2(i), Xgrid_f_3(j, i, k), dt,
                       Up_f(j, i, k, 1), Up_f(j, i, k, 2),
                       Up_f(j, i, k, 3), flag);
      }
    }
  }

  //  ! Construct the system matrix for computing ghost points values on the
  //  interface
  // ! We have 3*n1_c*n2_c equations in 3D
  // ! There are three sets of equations, one for the first component of u, one
  // for the ! second component and another for the third component of u ! Each
  // set consists of n1_c*n2_c equations

  // ! construct the diagnoal block jacobian matrix (block is 3x3)

  // interface_block(Rop.get(),ghcof.get(),Sb.get(),rho_c.get(),lambda_c.get(),rho_f.get(),Jacobian_c.get(),mu_c.get(),
  // 		  P.get(),XI13_c.get(),XI23_c.get(),XI33_c.get(),Mass_block.get());
  interface_block(Rop, ghcof, Sb, rho_c, lambda_c, rho_f, Jacobian_c, mu_c, P,
                  XI13_c, XI23_c, XI33_c, Mass_block, a);

  // std::cout<<"POST INTERFACE BLOCK  ";
  // for(int ii=0;ii<9;ii++) std::cout<<Mass_block[ii]<<" ";
  // std::cout<<"\n";
  int three = 3;
  int INFO;
  for (int j = 1; j <= n2_c; j++) {
    for (int i = 1; i <= n1_c; i++) {
      dgetrf_wrap(&three, &three, &Mass_block(1, 1, i, j), &three,
                  &IPIV_block[0 + (i - 1) * 3 + (j - 1) * 3 * n1_c], &INFO);
      if (INFO != 0) {
        std::cerr << "LU Fails at (i,j) equals" << i << "," << j
                  << " INFO = " << INFO << " " << Mass_block(INFO, INFO, i, j)
                  << "\n";
        for (int l = 1; l <= 3; l++) {
          for (int m = 1; m <= 3; m++)
            std::cerr << Mass_block(l, m, i, j) << ",";
          std::cerr << "\n";
        }
      }
    }
  }

  // Before the time loop, we make the initial conditions compatible with
  // interface conditions
  injection(Um_f, Um_c, P, a, 1);
  injection(U_f, U_c, P, a, 2);
  // u_f matches here

  // update ghost points value for the interface with block jacobian iterative
  // method
  // interface_rhs(Vass.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
  // mu_c.get(),
  // 		mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
  // XI13_c.get(),
  // 		XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
  // 		Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
  // 		ux_cof.get(),ghcof.get(),acof.get(),bof.get(),1);

  interface_rhs(Vass, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f, lambda_c,
                lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c, XI13_f, XI23_f,
                XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp, Um_c, Um_f, Mass_f1,
                ux_cof, ghcof, acof, bof, a, 1);

  // interface_lhs(LHS.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
  // mu_c.get(),
  // 		 mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
  // XI13_c.get(),
  // 		 XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
  // 		 Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
  // 		 ux_cof.get(),ghcof.get(),acof.get(),bof.get(),1);

  interface_lhs(LHS, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f, lambda_c,
                lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c, XI13_f, XI23_f,
                XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp, Um_c, Um_f, Mass_f1,
                ux_cof, ghcof, acof, bof, a, 1);

  for (int i = 1; i <= n1_c * n2_c * 3; i++) residual(i) = Vass(i) - LHS(i);

  float_sw4 tol = 1e-7;
  int iter = 0;
  float_sw4 res;
  while ((res = residual.maxabs()) > tol) {
    // std::cout<<"Iteration "<<iter++<<" "<<std::setprecision(17)<<res<<"\n";
    for (int j = 1; j <= n2_c; j++) {
      for (int i = 1; i <= n1_c; i++) {
        // if ((i==1) && (j==1)) {
        //   std::cout<<"PRE RHS"<<residual((j-1)*3*n1_c+3*(i-1)+1)<<"
        //   "<<residual((j-1)*3*n1_c+3*(i-1)+2)<<"
        //   "<<residual((j-1)*3*n1_c+3*(i-1)+3)<<"\n"; std::cout<<"M "; for(int
        //   jj=1;jj<=3;jj++) for(int ii=1;ii<=3;ii++)
        //   std::cout<<Mass_block(ii,jj,i,j)<<" "; std::cout<<"\n"; for(int
        //   ii=0;ii<9;ii++) std::cout<<Mass_block[ii]<<" "; std::cout<<"\n";
        // }
        dgetrs_wrap(3, 1, &Mass_block(1, 1, i, j), 3,
                    &IPIV_block[0 + (i - 1) * 3 + (j - 1) * 3 * n1_c],
                    &residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 1), 3, &INFO);
        // if ((i==1) && (j==1))
        // std::cout<<"RHS"<<residual((j-1)*3*n1_c+3*(i-1)+1)<<"
        // "<<residual((j-1)*3*n1_c+3*(i-1)+2)<<"
        // "<<residual((j-1)*3*n1_c+3*(i-1)+3)<<"\n";
        if (INFO != 0) {
          std::cerr << "SOLVE Fails at (i,j) equals" << i << "," << j
                    << " INFO = " << INFO << " " << Mass_block(INFO, INFO, i, j)
                    << "\n";
        }
        Um_c(i, j, n3_c + 1, 1) =
            Um_c(i, j, n3_c + 1, 1) +
            residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 1);
        Um_c(i, j, n3_c + 1, 2) =
            Um_c(i, j, n3_c + 1, 2) +
            residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 2);
        Um_c(i, j, n3_c + 1, 3) =
	  Um_c(i, j, n3_c + 1, 3) +
            residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 3);
      }
    }
    // interface_lhs(LHS.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
    // mu_c.get(),
    // 		   mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
    // XI13_c.get(),
    // 		   XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
    // 		   Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
    // 		   ux_cof.get(),ghcof.get(),acof.get(),bof.get(),1);

    interface_lhs(LHS, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f, lambda_c,
                  lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c, XI13_f,
                  XI23_f, XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp, Um_c, Um_f,
                  Mass_f1, ux_cof, ghcof, acof, bof, a, 1);

    for (int i = 1; i <= n1_c * n2_c * 3; i++) residual(i) = Vass(i) - LHS(i);
  }
  std::cout
      << "LAST ITERATION SHOULD BE 8 WITH RESIDUAL 6.6760102157559231E-007\n";

  tv = 0.0;
  // time stepping
  // nt = 2; // FOR DEBUGGIN
  // std::cerr<<" ****** WARNING ***********NT HAS BEEN RESTET FOR TESTING\n";

  //  std::shared_ptr<Farray> u_c_t = u_c.subset(2); // PBUGS USE U_C DIRECTLY
  //std::shared_ptr<Farray> u_f_t = u_f.subset(2);

  t1 = std::chrono::high_resolution_clock::now();
  for (int time_index = 1; time_index <= nt; time_index++) {
    std::cout << "Time step is " << time_index << "\n";
    // Predictor step

    // update ghost points value for the interface with block jacobian iterative
    // method
    // interface_rhs(Vass.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
    // mu_c.get(),
    // 		   mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
    // XI13_c.get(),
    // 		   XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
    // 		   Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
    // 		   ux_cof.get(),ghcof.get(),acof.get(),bof.get(),2);

    interface_rhs(Vass, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f,
                  lambda_c, lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c,
                  XI13_f, XI23_f, XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp,
                  U_c, U_f, Mass_f1, ux_cof, ghcof, acof, bof, a, 2);

    // interface_lhs(LHS.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
    // mu_c.get(),
    // 		   mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
    // XI13_c.get(),
    // 		   XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
    // 		   Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
    // 		   ux_cof.get(),ghcof.get(),acof.get(),bof.get(),2);

    interface_lhs(LHS, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f, lambda_c,
                  lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c, XI13_f,
                  XI23_f, XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp, U_c, U_f,
                  Mass_f1, ux_cof, ghcof, acof, bof, a, 2);

    for (int i = 1; i <= n1_c * n2_c * 3; i++) residual(i) = Vass(i) - LHS(i);

    iter = 0;
    while ((res = residual.maxabs()) > tol) {
      std::cout << "Iteration predictor " << iter++ << " "
                << std::setprecision(17) << res << "\n";
      for (int j = 1; j <= n2_c; j++) {
        for (int i = 1; i <= n1_c; i++) {
          dgetrs_wrap(3, 1, &Mass_block(1, 1, i, j), 3,
                      &IPIV_block[0 + (i - 1) * 3 + (j - 1) * 3 * n1_c],
                      &residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 1), 3,
                      &INFO);

          if (INFO != 0) {
            std::cerr << "SOLVE 2 Fails at (i,j) equals" << i << "," << j
                      << " INFO = " << INFO << " "
                      << Mass_block(INFO, INFO, i, j) << "\n";
          }

          U_c(i, j, n3_c + 1, 1) =
              U_c(i, j, n3_c + 1, 1) +
              residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 1);
          U_c(i, j, n3_c + 1, 2) =
              U_c(i, j, n3_c + 1, 2) +
              residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 2);
          U_c(i, j, n3_c + 1, 3) =
              U_c(i, j, n3_c + 1, 3) +
              residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 3);
        }
      }
      // interface_lhs(LHS.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
      // mu_c.get(),
      // 		     mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
      // XI13_c.get(),
      // 		     XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
      // 		     Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
      // 		     ux_cof.get(),ghcof.get(),acof.get(),bof.get(),2);

      interface_lhs(LHS, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f,
                    lambda_c, lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c,
                    XI13_f, XI23_f, XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp,
                    U_c, U_f, Mass_f1, ux_cof, ghcof, acof, bof, a, 2);

      for (int i = 1; i <= n1_c * n2_c * 3; i++) residual(i) = Vass(i) - LHS(i);
    }

    // Evaluate the difference operators
    int s = 1 - nrg;
    // update_interior(&u_c(s,s,s,1,2),&u_f(s,s,s,1,2),bof.get(),ghcof.get(),acof.get(),acof_no_gp.get(),
    // 		     lh_f.get(),Jacobian_f.get(),mu_f.get(),lambda_f.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),
    // 		     lh_c.get(),Jacobian_c.get(),mu_c.get(),lambda_c.get(),XI13_c.get(),XI23_c.get(),XI33_c.get());

    // std::cout<<"USE COUNT IS "<<u_c_t.use_count()<<"\n";
    update_interior(U_c, U_f, bof, ghcof, acof, acof_no_gp, lh_f,
                    Jacobian_f, mu_f, lambda_f, XI13_f, XI23_f, XI33_f, lh_c,
                    Jacobian_c, mu_c, lambda_c, XI13_c, XI23_c, XI33_c, a);

    // Compute forcing functions
    for (int k = 1 - nrg; k <= n3_c + nrg; k++) {
      for (int i = 1 - nrg; i <= n2_c + nrg; i++) {
        for (int j = 1 - nrg; j <= n1_c + nrg; j++) {
          // forcing(&Xgrid_c_1(j),&Xgrid_c_2(i),&Xgrid_c_3(j,i,k),&tv,
          // 	   &force_c(j,i,k,1),&force_c(j,i,k,2),&force_c(j,i,k,3));
          forcing(Xgrid_c_1(j), Xgrid_c_2(i), Xgrid_c_3(j, i, k), tv,
                  force_c(j, i, k, 1), force_c(j, i, k, 2),
                  force_c(j, i, k, 3));
        }
      }
    }

    for (int k = 1 - nrg; k <= n3_f + nrg; k++) {
      for (int i = 1 - nrg; i <= n2_f + nrg; i++) {
        for (int j = 1 - nrg; j <= n1_f + nrg; j++) {
          // forcing(&Xgrid_f_1(j),&Xgrid_f_2(i),&Xgrid_f_3(j,i,k),&tv,
          //  	   &force_f(j,i,k,1),&force_f(j,i,k,2),&force_f(j,i,k,3));
          forcing(Xgrid_f_1(j), Xgrid_f_2(i), Xgrid_f_3(j, i, k), tv,
                  force_f(j, i, k, 1), force_f(j, i, k, 2),
                  force_f(j, i, k, 3));
        }
      }
    }
    // Update the solution in the predictor step

    for (int k = 1; k <= n3_f; k++) {
      for (int j = 1; j <= n2_f; j++) {
        for (int i = 1; i <= n1_f; i++) {
          Up_f(i, j, k, 1) = 2.0 * U_f(i, j, k, 1) - Um_f(i, j, k, 1) +
                               dt * dt *
                                   (lh_f(i, j, k, 1) +
                                    Jacobian_f(i, j, k) * force_f(i, j, k, 1)) /
                                   rho_f(i, j, k);
          Up_f(i, j, k, 2) = 2.0 * U_f(i, j, k, 2) - Um_f(i, j, k, 2) +
                               dt * dt *
                                   (lh_f(i, j, k, 2) +
                                    Jacobian_f(i, j, k) * force_f(i, j, k, 2)) /
                                   rho_f(i, j, k);
          Up_f(i, j, k, 3) = 2.0 * U_f(i, j, k, 3) - Um_f(i, j, k, 3) +
                               dt * dt *
                                   (lh_f(i, j, k, 3) +
                                    Jacobian_f(i, j, k) * force_f(i, j, k, 3)) /
                                   rho_f(i, j, k);
        }
      }
    }
    for (int k = 1; k <= n3_c; k++) {
      for (int j = 1; j <= n2_c; j++) {
        for (int i = 1; i <= n1_c; i++) {
          Up_c(i, j, k, 1) = 2.0 * U_c(i, j, k, 1) - Um_c(i, j, k, 1) +
                               dt * dt *
                                   (lh_c(i, j, k, 1) +
                                    Jacobian_c(i, j, k) * force_c(i, j, k, 1)) /
                                   rho_c(i, j, k);
          Up_c(i, j, k, 2) = 2.0 * U_c(i, j, k, 2) - Um_c(i, j, k, 2) +
                               dt * dt *
                                   (lh_c(i, j, k, 2) +
                                    Jacobian_c(i, j, k) * force_c(i, j, k, 2)) /
                                   rho_c(i, j, k);
          Up_c(i, j, k, 3) = 2.0 * U_c(i, j, k, 3) - Um_c(i, j, k, 3) +
                               dt * dt *
                                   (lh_c(i, j, k, 3) +
                                    Jacobian_c(i, j, k) * force_c(i, j, k, 3)) /
                                   rho_c(i, j, k);
        }
      }
    }
    // print_uf(u_f.get());
    // print_f(force_f.get());
    // Update time
    tv = tv + dt;

    // Update ghost points outside the left and right boundaries
    // The argument '3' means time level star
    //'1', '2' and '4' mean time level n-1, n and n+1, respectively.
    // update_gp(Xgrid_c_1.get(),Xgrid_c_2.get(),Xgrid_c_3.get(),u_c.get(),Xgrid_f_1.get(),Xgrid_f_2.get(),Xgrid_f_3.get(),u_f.get(),&tv,3);
    update_gp(Xgrid_c_1, Xgrid_c_2, Xgrid_c_3, Up_c, Xgrid_f_1, Xgrid_f_2,
              Xgrid_f_3, Up_f, tv, a, 3);
    // Update ghost point values for the traction boundary
    // update_traction(traction_rhs,traction_data.get(),Xgrid_f_1.get(),Xgrid_f_2.get(),Xgrid_f_3.get(),
    // 		     u_f.get(),mu_f.get(),lambda_f.get(),Jacobian_f.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),Sb.get(),&tv,3);

    // float_sw4 tmps[3];
    // Farray u_f_tmp(1-nrg,n1_f+nrg,1-nrg,n2_f+nrg,1-nrg,n3_f+nrg,1,dim,1,4);
    // u_f_tmp=u_f;
    update_traction(traction_data, Xgrid_f_1, Xgrid_f_2, Xgrid_f_3, Up_f, mu_f,
                    lambda_f, Jacobian_f, XI13_f, XI23_f, XI33_f, Sb, tv, a, 3);

    // u_f_tmp.compare(u_f);
    // Injection at the interface
    injection(Up_f, Up_c, P, a, 3);

    // Update Dirichlet boundary condition
    // update_dirichlet_bc(Xgrid_c_1.get(),Xgrid_c_2.get(),Xgrid_c_3.get(),u_c.get(),&tv,3);
    update_dirichlet_bc(Xgrid_c_1, Xgrid_c_2, Xgrid_c_3, Up_c, tv, a, 3);

    // Corrector step

    // update ghost points value for the interface with block jacobian iterative
    // method
    // interface_rhs(Vass.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
    // mu_c.get(),
    // 		   mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
    // XI13_c.get(),
    // 		   XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
    // 		   Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
    // 		   ux_cof.get(),ghcof.get(),acof.get(),bof.get(),3);

    interface_rhs(Vass, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f,
                  lambda_c, lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c,
                  XI13_f, XI23_f, XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp,
                  Up_c, Up_f, Mass_f1, ux_cof, ghcof, acof, bof, a, 3);

    // interface_lhs(LHS.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
    // mu_c.get(),
    // 		   mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
    // XI13_c.get(),
    // 		   XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
    // 		   Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
    // 		   ux_cof.get(),ghcof.get(),acof.get(),bof.get(),3);
    interface_lhs(LHS, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f, lambda_c,
                  lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c, XI13_f,
                  XI23_f, XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp, Up_c, Up_f,
                  Mass_f1, ux_cof, ghcof, acof, bof, a, 3);

    for (int i = 1; i <= n1_c * n2_c * 3; i++) residual(i) = Vass(i) - LHS(i);

    iter = 0;
    while ((res = residual.maxabs()) > tol) {
      std::cout << "Iteration corrector  " << iter++ << " "
                << std::setprecision(17) << res << "\n";
      for (int j = 1; j <= n2_c; j++) {
        for (int i = 1; i <= n1_c; i++) {
          dgetrs_wrap(3, 1, &Mass_block(1, 1, i, j), 3,
                      &IPIV_block[0 + (i - 1) * 3 + (j - 1) * 3 * n1_c],
                      &residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 1), 3,
                      &INFO);

          if (INFO != 0) {
            std::cerr << "SOLVE 3 Fails at (i,j) equals" << i << "," << j
                      << " INFO = " << INFO << " "
                      << Mass_block(INFO, INFO, i, j) << "\n";
          }

          Up_c(i, j, n3_c + 1, 1) =
              Up_c(i, j, n3_c + 1, 1) +
              residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 1);
          Up_c(i, j, n3_c + 1, 2) =
              Up_c(i, j, n3_c + 1, 2) +
              residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 2);
          Up_c(i, j, n3_c + 1, 3) =
              Up_c(i, j, n3_c + 1, 3) +
              residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 3);
        }
      }
      // interface_lhs(LHS.get(),lh_c.get(),lh_f.get(),Jacobian_c.get(),Jacobian_f.get(),
      // mu_c.get(),
      // 		     mu_f.get(),lambda_c.get(),lambda_f.get(),rho_c.get(),rho_f.get(),
      // XI13_c.get(),
      // 		     XI23_c.get(),XI33_c.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),P.get(),Sb.get(),
      // 		     Rop.get(),sbop_no_gp.get(),acof_no_gp.get(),u_c.get(),u_f.get(),Mass_f1.get(),
      // 		     ux_cof.get(),ghcof.get(),acof.get(),bof.get(),3);

      interface_lhs(LHS, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f,
                    lambda_c, lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c,
                    XI13_f, XI23_f, XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp,
                    Up_c, Up_f, Mass_f1, ux_cof, ghcof, acof, bof, a, 3);

      for (int i = 1; i <= n1_c * n2_c * 3; i++) residual(i) = Vass(i) - LHS(i);
    }

    // Evaluate the difference operators

    for (int l = 1; l <= dim; l++)
      for (int k = 1 - nrg; k <= n3_f + nrg; k++)
        for (int j = 1 - nrg; j <= n2_f + nrg; j++)
          for (int i = 1 - nrg; i <= n1_f + nrg; i++)
            u_ft(i, j, k, l) = (Up_f(i, j, k, l) - 2 * U_f(i, j, k, l) +
                                Um_f(i, j, k, l)) /
                               (dt * dt);

    for (int l = 1; l <= dim; l++)
      for (int k = 1 - nrg; k <= n3_c + nrg; k++)
        for (int j = 1 - nrg; j <= n2_c + nrg; j++)
          for (int i = 1 - nrg; i <= n1_c + nrg; i++)
            u_ct(i, j, k, l) = (Up_c(i, j, k, l) - 2 * U_c(i, j, k, l) +
                                Um_c(i, j, k, l)) /
                               (dt * dt);

    s = 1 - nrg;
    // update_interior(&u_ct(s,s,s,1),&u_ft(s,s,s,1),bof.get(),ghcof.get(),acof.get(),acof_no_gp.get(),
    // 		     lh_f.get(),Jacobian_f.get(),mu_f.get(),lambda_f.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),
    // 		     lh_c.get(),Jacobian_c.get(),mu_c.get(),lambda_c.get(),XI13_c.get(),XI23_c.get(),XI33_c.get());

    update_interior(u_ct, u_ft, bof, ghcof, acof, acof_no_gp, lh_f, Jacobian_f,
                    mu_f, lambda_f, XI13_f, XI23_f, XI33_f, lh_c, Jacobian_c,
                    mu_c, lambda_c, XI13_c, XI23_c, XI33_c, a);

    // Compute the second time derivative of the forcing fucntions
    for (int k = 1 - nrg; k <= n3_c + nrg; k++) {
      for (int i = 1 - nrg; i <= n2_c + nrg; i++) {
        for (int j = 1 - nrg; j <= n1_c + nrg; j++) {
          // forcing_tt(&Xgrid_c_1(j),&Xgrid_c_2(i),&Xgrid_c_3(j,i,k),tv-dt,
          // 	   &force_tt_c(j,i,k,1),&force_tt_c(j,i,k,2),&force_tt_c(j,i,k,3));
          forcing_tt(Xgrid_c_1(j), Xgrid_c_2(i), Xgrid_c_3(j, i, k), tv - dt,
                     force_tt_c(j, i, k, 1), force_tt_c(j, i, k, 2),
                     force_tt_c(j, i, k, 3));
        }
      }
    }

    for (int k = 1 - nrg; k <= n3_f + nrg; k++) {
      for (int i = 1 - nrg; i <= n2_f + nrg; i++) {
        for (int j = 1 - nrg; j <= n1_f + nrg; j++) {
          // forcing_tt(&Xgrid_f_1(j),&Xgrid_f_2(i),&Xgrid_f_3(j,i,k),tv-dt,
          //  	   &force_tt_f(j,i,k,1),&force_tt_f(j,i,k,2),&force_tt_f(j,i,k,3));
          forcing_tt(Xgrid_f_1(j), Xgrid_f_2(i), Xgrid_f_3(j, i, k), tv - dt,
                     force_tt_f(j, i, k, 1), force_tt_f(j, i, k, 2),
                     force_tt_f(j, i, k, 3));
        }
      }
    }

    // Evaluate the difference operators
    float_sw4 dt4 = dt * dt * dt * dt;
    for (int k = 1; k <= n3_f; k++) {
      for (int j = 1; j <= n2_f; j++) {
        for (int i = 1; i <= n1_f; i++) {
          Ux_f(i, j, k, 1) =
              Up_f(i, j, k, 1) +
              dt4 / 12.0 *
                  (lh_f(i, j, k, 1) +
                   Jacobian_f(i, j, k) * force_tt_f(i, j, k, 1)) /
                  rho_f(i, j, k);
          Ux_f(i, j, k, 2) =
              Up_f(i, j, k, 2) +
              dt4 / 12.0 *
                  (lh_f(i, j, k, 2) +
                   Jacobian_f(i, j, k) * force_tt_f(i, j, k, 2)) /
                  rho_f(i, j, k);
          Ux_f(i, j, k, 3) =
              Up_f(i, j, k, 3) +
              dt4 / 12.0 *
                  (lh_f(i, j, k, 3) +
                   Jacobian_f(i, j, k) * force_tt_f(i, j, k, 3)) /
                  rho_f(i, j, k);
        }
      }
    }
    for (int k = 1; k <= n3_c; k++) {
      for (int j = 1; j <= n2_c; j++) {
        for (int i = 1; i <= n1_c; i++) {
          Ux_c(i, j, k, 1) =
              Up_c(i, j, k, 1) +
              dt4 / 12.0 *
                  (lh_c(i, j, k, 1) +
                   Jacobian_c(i, j, k) * force_tt_c(i, j, k, 1)) /
                  rho_c(i, j, k);
          Ux_c(i, j, k, 2) =
              Up_c(i, j, k, 2) +
              dt4 / 12.0 *
                  (lh_c(i, j, k, 2) +
                   Jacobian_c(i, j, k) * force_tt_c(i, j, k, 2)) /
                  rho_c(i, j, k);
          Ux_c(i, j, k, 3) =
              Up_c(i, j, k, 3) +
              dt4 / 12.0 *
                  (lh_c(i, j, k, 3) +
                   Jacobian_c(i, j, k) * force_tt_c(i, j, k, 3)) /
                  rho_c(i, j, k);
        }
      }
    }

    // Update ghost point values outside the left and right boundaries
    // update_gp(Xgrid_c_1.get(),Xgrid_c_2.get(),Xgrid_c_3.get(),u_c.get(),Xgrid_f_1.get(),Xgrid_f_2.get(),Xgrid_f_3.get(),u_f.get(),&tv,4);

    update_gp(Xgrid_c_1, Xgrid_c_2, Xgrid_c_3, Ux_c, Xgrid_f_1, Xgrid_f_2,
              Xgrid_f_3, Ux_f, tv, a, 4);
    // Update ghost point values for the traction boundary
    // update_traction(traction_rhs,traction_data.get(),Xgrid_f_1.get(),Xgrid_f_2.get(),Xgrid_f_3.get(),
    // 		     u_f.get(),mu_f.get(),lambda_f.get(),Jacobian_f.get(),XI13_f.get(),XI23_f.get(),XI33_f.get(),Sb.get(),&tv,4);

    update_traction(traction_data, Xgrid_f_1, Xgrid_f_2, Xgrid_f_3, Ux_f, mu_f,
                    lambda_f, Jacobian_f, XI13_f, XI23_f, XI33_f, Sb, tv, a, 4);
    // injection at the interface
    injection(Ux_f, Ux_c, P, a, 4);

    // Update Dirichlet boundary condition
    // update_dirichlet_bc(Xgrid_c_1.get(),Xgrid_c_2.get(),Xgrid_c_3.get(),u_c.get(),&tv,4);
    update_dirichlet_bc(Xgrid_c_1, Xgrid_c_2, Xgrid_c_3, Ux_c, tv, a, 4);

    // Update solutions
    for (int l = 1; l <= dim; l++)
      for (int k = 1 - nrg; k <= n3_f + nrg; k++)
        for (int j = 1 - nrg; j <= n2_f + nrg; j++)
          for (int i = 1 - nrg; i <= n1_f + nrg; i++) {
            Um_f(i, j, k, l) = U_f(i, j, k, l);
            U_f(i, j, k, l) = Ux_f(i, j, k, l);
          }

    for (int l = 1; l <= dim; l++)
      for (int k = 1 - nrg; k <= n3_c + nrg; k++)
        for (int j = 1 - nrg; j <= n2_c + nrg; j++)
          for (int i = 1 - nrg; i <= n1_c + nrg; i++) {
            Um_c(i, j, k, l) = U_c(i, j, k, l);
            U_c(i, j, k, l) = Ux_c(i, j, k, l);
          }

  }  // TIME STEP LOOP

  t2 = std::chrono::high_resolution_clock::now();

  float_sw4 l2_err = 0.0;
  flag = 1;
  for (int k = 1; k <= n3_f; k++) {
    for (int i = 1; i <= n2_f; i++) {
      for (int j = 1; j <= n1_f; j++) {
        exact_solution(Xgrid_f_1(j), Xgrid_f_2(i), Xgrid_f_3(j, i, k), tv,
                       Up_f(j, i, k, 1), Up_f(j, i, k, 2),
                       Up_f(j, i, k, 3), flag);
        err_f(j, i, k) = std::max(
            std::fabs(Ux_f(j, i, k, 1) - Up_f(j, i, k, 1)),
            std::max(std::fabs(Ux_f(j, i, k, 2) - Up_f(j, i, k, 2)),
                     std::fabs(Ux_f(j, i, k, 3) - Up_f(j, i, k, 3))));
        l2_err =
            l2_err + h1_f * h2_f * h3_f *
                         (pow((Ux_f(j, i, k, 1) - Up_f(j, i, k, 1)), 2) +
                          pow((Ux_f(j, i, k, 2) - Up_f(j, i, k, 2)), 2) +
                          pow((Ux_f(j, i, k, 3) - Up_f(j, i, k, 3)), 2));
      }
    }
  }
  flag = 0;
  for (int k = 1; k <= n3_c; k++) {
    for (int i = 1; i <= n2_c; i++) {
      for (int j = 1; j <= n1_c; j++) {
        exact_solution(Xgrid_c_1(j), Xgrid_c_2(i), Xgrid_c_3(j, i, k), tv,
                       Up_c(j, i, k, 1), Up_c(j, i, k, 2),
                       Up_c(j, i, k, 3), flag);
        err_c(j, i, k) = std::max(
            std::fabs(Ux_c(j, i, k, 1) - Up_c(j, i, k, 1)),
            std::max(std::fabs(Ux_c(j, i, k, 2) - Up_c(j, i, k, 2)),
                     std::fabs(Ux_c(j, i, k, 3) - Up_c(j, i, k, 3))));
        l2_err =
            l2_err + h1_c * h2_c * h3_c *
                         (pow((Ux_c(j, i, k, 1) - Up_c(j, i, k, 1)), 2) +
                          pow((Ux_c(j, i, k, 2) - Up_c(j, i, k, 2)), 2) +
                          pow((Ux_c(j, i, k, 3) - Up_c(j, i, k, 3)), 2));
      }
    }
  }

  auto rtime =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  l2_err = sqrt(l2_err);
  std::cout << "No of grid points " << n1_c << " " << n2_c << " " << n3_c << " "
            << n1_f << " " << n2_f << " " << n3_f << "\n";
  std::cout << "error " << err_c.maximum() << " " << err_f.maximum() << " "
            << l2_err << "\n";
  std::cout << "error/16 " << err_c.maximum() / 16 << " "
            << err_f.maximum() / 16 << " " << l2_err / 16 << "\n";
  std::cout << "computational time " << rtime << " ms \n";
  std::cout << "cfl " << cfl << "\n";
  float_sw4 erra[3];
  erra[0] = err_c.maximum();
  erra[1] = err_f.maximum();
  erra[2] = l2_err;
  bool baseline = false;
  if (baseline) {
    std::fstream myfile;
    myfile.open("errors.bin", std::ios::out | std::ios::binary);
    myfile.write(reinterpret_cast<char *>(&erra[0]), sizeof(float_sw4) * 3);
    myfile.close();
  } else {
    std::fstream myfile;
    float_sw4 err_base[3];
    std::cout << "ERROR ANALYSES COMPARED TO BASELINE\n";
    myfile.open("baseline.bin", std::ios::in | std::ios::binary);
    myfile.read(reinterpret_cast<char *>(&err_base[0]), sizeof(float_sw4) * 3);
    myfile.close();
    for (int i = 0; i < 3; i++) {
      float_sw4 diff = (erra[i] - err_base[i]) / err_base[i] * 100.0;
      std::cout << i << " " << erra[i] << " " << err_base[i] << " " << diff
                << " %\n";
    }
  }

  // std::cerr<<" ****** WARNING ***********NT HAS BEEN RESTET FOR TESTING\n";
}
