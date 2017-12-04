#ifndef SW4_RHS4TH3C_H
#include "sw4.h"

void rhs4th3fortwind_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			 int nk, int* __restrict__ onesided, float_sw4* __restrict__ a_acof, 
			 float_sw4 *__restrict__ a_bope, float_sw4* __restrict__ a_ghcof, 
			 float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
			 float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, 
			 float_sw4 h, float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
			 float_sw4* __restrict__ a_strz, char op, int kfirstu, int klastu, int kfirstw,
			 int klastw );

void rhs4th3fortsgstr_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			  int nk, int* __restrict__ onesided, float_sw4* __restrict__ a_acof,
			  float_sw4 *__restrict__ a_bope, float_sw4* __restrict__ a_ghcof,
			  float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
			  float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, 
			  float_sw4 h, float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
			  float_sw4* __restrict__ a_strz, char op );

void rhs4th3fort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		     int nk, int* __restrict__ onesided, float_sw4* __restrict__ a_acof, 
		     float_sw4 *__restrict__ a_bope, float_sw4* __restrict__ a_ghcof, 
		     float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
		     float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, 
		     float_sw4 h, char op );

void predfort_ci( int ib, int ie, int jb, int je, int kb, int ke,
		  float_sw4* __restrict__ up, float_sw4* __restrict__ u,
		  float_sw4* __restrict__ um, float_sw4* __restrict__ lu,
		  float_sw4* __restrict__ fo, float_sw4* __restrict__ rho,
		  float_sw4 dt2 );

void corrfort_ci( int ib, int ie, int jb, int je, int kb, int ke,
		  float_sw4* __restrict__ up, float_sw4* __restrict__ lu,
		  float_sw4* __restrict__ fo, float_sw4* __restrict__ rho,
		  float_sw4 dt4 );

void dpdmtfort_ci( int ib, int ie, int jb, int je, int kb, int ke,
		   float_sw4* __restrict__ up, float_sw4* __restrict__ u,
		   float_sw4* __restrict__ um, float_sw4* __restrict__ u2,
		   float_sw4 dt2i );

void rhouttlumf_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		    int nz, float_sw4* __restrict__ a_uacc, float_sw4* __restrict__ a_lu,
		    float_sw4* __restrict__ a_fo, float_sw4* __restrict__ a_rho,
		    float_sw4 lowZ[3], float_sw4 interZ[3], float_sw4 highZ[3] );

void rhserrfort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		    int nz, float_sw4 h,
		    float_sw4* __restrict__ a_fo, float_sw4* __restrict__ a_u2,
		    float_sw4 lowZ[3], float_sw4 interZ[3], float_sw4 highZ[3] );

void dpdmtfortatt_ci( int ib, int ie, int jb, int je, int kb, int ke,
		      float_sw4* __restrict__ up, float_sw4* __restrict__ u,
		      float_sw4* __restrict__ um, float_sw4 dt2i );

void satt_ci( float_sw4* __restrict__ up, float_sw4* __restrict__ qs,
	      float_sw4 dt, float_sw4 cfreq, int ifirst, int ilast,
	      int jfirst, int jlast, int kfirst, int klast );

void solveattfreeac_ci( int ifirst, int ilast, int jfirst, int jlast,
			int kfirst, int klast,
			float_sw4* __restrict__ a_alpha, float_sw4 cof,
			float_sw4* __restrict__ a_up );

void solveattfreec_ci( int ifirst, int ilast, int jfirst, int jlast,
		       int kfirst, int klast, float_sw4* __restrict__ a_u,
		       float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_la,
		       float_sw4* __restrict__ a_muve, float_sw4* __restrict__ a_lave,
		       float_sw4* __restrict__ a_bforcerhs, float_sw4* __restrict__ a_met,
		       float_sw4 s[5], int usesg, float_sw4* __restrict__ a_strx,
		       float_sw4* __restrict__ a_stry );

void addbstresswresc_ci( int ifirst, int ilast, int jfirst, int jlast,
			 int kfirst, int klast, int nz, float_sw4* __restrict__ a_alphap,
			 float_sw4* __restrict__ a_alpham, float_sw4* __restrict__ a_muve,
			 float_sw4* __restrict__ a_lave, float_sw4* __restrict__ a_bforcerhs,
			 float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_um,
			 float_sw4* __restrict__ a_met, int side, float_sw4 dt, float_sw4 omegave,
			 float_sw4* __restrict__ a_memforce, float_sw4* __restrict__ a_muvebnd, 
			 float_sw4* __restrict__ a_lambdavebnd, float_sw4 s[5], float_sw4& cof,
			 int usesg, float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry );

#endif
