void interpolation_restriction(Farray &P, Farray &Rop, Farray &RPop);
void central_difference(Farray &ux_cof, Farray &uxx_cof);
void varcoeffs4(Farray &acof, Farray &ghcof, Farray &Sb);
void varcoeff_noghost(Farray &acof_no_gp, Farray &ghcof_no_gp,
                      Farray &sbop_no_gp);
void dx_46(Farray &bop);
void equation_cof(Farray &Xgrid_c_1, Farray &Xgrid_c_2, Farray &Xgrid_c_3,
                  Farray &Xgrid_f_1, Farray &Xgrid_f_2, Farray &Xgrid_f_3,
                  Sarray &XI13_c, Sarray &XI23_c, Sarray &XI33_c,
                  Sarray &Jacobian_c, Sarray &rho_c, Sarray &rho_f,
                  Sarray &XI13_f, Sarray &XI23_f, Sarray &XI33_f,
                  Sarray &Jacobian_f, Sarray &mu_c, Sarray &mu_f,
                  Sarray &lambda_c, Sarray &lambda_f, PackArgs &a);
/* void equation_cof(Farray &Xgrid_c_1, Farray &Xgrid_c_2, Farray &Xgrid_c_3,
 * Farray &Xgrid_f_1, Farray &Xgrid_f_2, Farray &Xgrid_f_3,  */
/* 		  Farray &XI13_c, Farray &XI23_c, Farray &XI33_c, Farray
 * &Jacobian_c, Farray &rho_c, Farray &rho_f, */
/* 		  Farray &XI13_f, Farray &XI23_f, Farray &XI33_f,Farray
 * &Jacobian_f, Farray &mu_c, Farray &mu_f, Farray &lambda_c, */
/* 		  Farray &lambda_f,PackArgs &a); */
void generate_grid(Farray &X1_c, Farray &X2_c, Farray &X3_c, Farray &X1_f,
                   Farray &X2_f, Farray &X3_f, PackArgs &a);
void kappa(Sarray &mu_f, Sarray &lambda_f, Sarray &XI13_f, Sarray &XI23_f,
           Sarray &XI33_f, Sarray &Jacobian_f, Sarray &rho_f, float_sw4 &kappa2,
           PackArgs &a);
void exact_solution(float_sw4 x1, float_sw4 x2, float_sw4 x3, float_sw4 t,
                    float_sw4 &u1, float_sw4 &u2, float_sw4 &u3, int flag);
void interface_block(Farray &Rop, Farray &ghcof, Farray &Sb, Sarray &rho_c,
                     Sarray &lambda_c, Sarray &rho_f, Sarray &Jacobian_c,
                     Sarray &mu_c, Farray &P, Sarray &XI13_c, Sarray &XI23_c,
                     Sarray &XI33_c, Farray &Mass_block, PackArgs &a);
void interface_lhs(Farray &LHS, Farray &lh_c, Farray &lh_f, Sarray &Jacobian_c,
                   Sarray &Jacobian_f, Sarray &mu_c, Sarray &mu_f,
                   Sarray &lambda_c, Sarray &lambda_f, Sarray &rho_c,
                   Sarray &rho_f, Sarray &XI13_c, Sarray &XI23_c,
                   Sarray &XI33_c, Sarray &XI13_f, Sarray &XI23_f,
                   Sarray &XI33_f, Farray &P, Farray &Sb, Farray &Rop,
                   Farray &sbop_no_gp, Farray &acof_no_gp, Sarray &u_c,
                   Sarray &u_f, Farray &Mass_f1, Farray &ux_cof, Farray &ghcof,
                   Farray &acof, Farray &bof, PackArgs &a, int index);
void injection(Sarray &u_f, Sarray &u_c, Farray &P, PackArgs &a, int index);

void interface_rhs(Farray &Vass, Farray &lh_c, Farray &lh_f, Sarray &Jacobian_c,
                   Sarray &Jacobian_f, Sarray &mu_c, Sarray &mu_f,
                   Sarray &lambda_c, Sarray &lambda_f, Sarray &rho_c,
                   Sarray &rho_f, Sarray &XI13_c, Sarray &XI23_c,
                   Sarray &XI33_c, Sarray &XI13_f, Sarray &XI23_f,
                   Sarray &XI33_f, Farray &P, Farray &Sb, Farray &Rop,
                   Farray &sbop_no_gp, Farray &acof_no_gp, Sarray &u_c,
                   Sarray &u_f, Farray &Mass_f1, Farray &ux_cof, Farray &ghcof,
                   Farray &acof, Farray &bof, PackArgs &a, int index);
void update_interior(Sarray &u_c_t, Sarray &u_f_t, Farray &bof, Farray &ghcof,
                     Farray &acof, Farray &acof_no_gp, Farray &lh_f,
                     Sarray &Jacobian_f, Sarray &mu_f, Sarray &lambda_f,
                     Sarray &XI13_f, Sarray &XI23_f, Sarray &XI33_f,
                     Farray &lh_c, Sarray &Jacobian_c, Sarray &mu_c,
                     Sarray &lambda_c, Sarray &XI13_c, Sarray &XI23_c,
                     Sarray &XI33_c, PackArgs &a);

void forcing(float_sw4 x1, float_sw4 x2, float_sw4 x3, float_sw4 t,
             float_sw4 &f1, float_sw4 &f2, float_sw4 &f3);
void forcing_tt(float_sw4 x1, float_sw4 x2, float_sw4 x3, float_sw4 t,
                float_sw4 &f1tt, float_sw4 &f2tt, float_sw4 &f3tt);
void update_gp(Farray &Xgrid_c_1, Farray &Xgrid_c_2, Farray &Xgrid_c_3,
               Sarray &u_c, Farray &Xgrid_f_1, Farray &Xgrid_f_2,
               Farray &Xgrid_f_3, Sarray &u_f, float_sw4 tv, PackArgs &a,
               int index);

void update_traction(Sarray &traction_data, Farray &Xgrid_f_1,
                     Farray &Xgrid_f_2, Farray &Xgrid_f_3, Sarray &u_f,
                     Sarray &mu_f, Sarray &lambda_f, Sarray &Jacobian_f,
                     Sarray &XI13_f, Sarray &XI23_f, Sarray &XI33_f, Farray &Sb,
                     float_sw4 tv, PackArgs &a, int index);
void top_normal_data(float_sw4 x1, float_sw4 x2, float_sw4 x3, float_sw4 l1,
                     float_sw4 l2, float_sw4 t, int i, int j, Sarray &g,
                     PackArgs &a);
void update_dirichlet_bc(Farray &Xgrid_c_1, Farray &Xgrid_c_2,
                         Farray &Xgrid_c_3, Sarray &u_c, float_sw4 tv,
                         PackArgs &a, int index);
