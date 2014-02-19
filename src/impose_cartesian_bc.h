#ifndef impose_cartesian_bc_h
#define impose_cartesian_bc_h
void
impose_cartesian_bc(Sarray & u,
		    Sarray & mu,
		    Sarray & lam,
		    double ** bcForcing,
		    boundaryConditionType bcType[],
		    int onesided[], double curlcoeff, double h );


/* void */
/* cartesian_bc_forcing(Sarray & mu, */
/* 		     double ** bcForcing, */
/* 		     int * numberOfBCPoints, */
/* 		     boundaryConditionType bcType[], */
/* 		     Forcing* forcing, double zmin, double t, double h, int onesided[], double curlcoeff, */
/* 		     int numberof_mechanisms, Sarray* AlphaVE, Sarray* MuVE, Sarray* LambdaVE ); */
#endif
