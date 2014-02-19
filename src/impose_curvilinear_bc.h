#ifndef impose_curvilinear_bc_h
#define impose_curvilinear_bc_h
void
impose_curvilinear_bc(Sarray & a_u, double ** bcForcing, Sarray & a_x, Sarray & a_y, Sarray & a_z,
		      Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J,
		      double t, boundaryConditionType bcType[] );
void eval_curvilinear_bc_stress(Sarray & a_u, double ** bcForcing, Sarray & a_x, Sarray & a_y, Sarray & a_z,
				Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J);
/* void curvilinear_bc_forcing(double ** bcForcing, int * numberOfBCPoints, Sarray & a_x, Sarray & a_y, Sarray & a_z, */
/* 			    Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J, */
/* 			    double t, boundaryConditionType bcType[], Forcing * f, */
/* 			    int noof_mechanisms, Sarray* alphaVE, Sarray* muVE, Sarray* lambdaVE ); */

#endif
