#ifndef SW4_CF_INTERFACE_H
#include "sw4.h"
#include "Sarray.h"

// the extern "c" is only needed for linking the Fortran version of these routines
#ifdef SW4_NOC
extern "C" {
#endif

void evenIevenJinterpJacobi(float_sw4 rmax[6], Sarray &Uf, Sarray &UfNew, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
			    Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
			    Sarray &Morc, Sarray &Mlrc,
			    Sarray &Unextf, Sarray &Bf, Sarray &UnextcInterp, Sarray &Bc,
			    int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
			    int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
			    float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
			    float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void evenIevenJinterpOpt(float_sw4 rmax[6], float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_muf, 
			 float_sw4* __restrict__ a_lambdaf, float_sw4* __restrict__ a_rhof, 
			 float_sw4* __restrict__ a_uc, float_sw4* __restrict__ a_muc, 
			 float_sw4* __restrict__ a_lambdac, float_sw4* __restrict__ a_rhoc,
			 float_sw4* __restrict__ a_morc, float_sw4* __restrict__ a_mlrc,
			 float_sw4* __restrict__ a_unextf, float_sw4* __restrict__ a_bf, 
			 float_sw4* __restrict__ a_unextc, float_sw4* __restrict__ a_bc,
			 int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[], 
			 int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
			 int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
			 float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
			 float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void evenIoddJinterpOpt(float_sw4 rmax[6], float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_muf, 
			float_sw4* __restrict__ a_lambdaf, float_sw4* __restrict__ a_rhof, 
			float_sw4* __restrict__ a_uc, float_sw4* __restrict__ a_muc, 
			float_sw4* __restrict__ a_lambdac, float_sw4* __restrict__ a_rhoc,
			float_sw4* __restrict__ a_morc, float_sw4* __restrict__ a_mlrc,
			float_sw4* __restrict__ a_unextf, float_sw4* __restrict__ a_bf, 
			float_sw4* __restrict__ a_unextc, float_sw4* __restrict__ a_bc,
			int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[], 
			int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
			int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
			float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
			float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIevenJinterpOpt(float_sw4 rmax[6], float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_muf, 
			float_sw4* __restrict__ a_lambdaf, float_sw4* __restrict__ a_rhof, 
			float_sw4* __restrict__ a_uc, float_sw4* __restrict__ a_muc, 
			float_sw4* __restrict__ a_lambdac, float_sw4* __restrict__ a_rhoc,
			float_sw4* __restrict__ a_morc, float_sw4* __restrict__ a_mlrc,
			float_sw4* __restrict__ a_unextf, float_sw4* __restrict__ a_bf, 
			float_sw4* __restrict__ a_unextc, float_sw4* __restrict__ a_bc,
			int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[], 
			int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
			int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
			float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
			float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterpOpt(float_sw4 rmax[3], float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_muf, 
		       float_sw4* __restrict__ a_lambdaf, float_sw4* __restrict__ a_rhof, 
		       float_sw4* __restrict__ a_uc, float_sw4* __restrict__ amuc, 
		       float_sw4* __restrict__ a_lambdac, float_sw4* __restrict__ a_rhoc,
		       float_sw4* __restrict__ a_mufs, float_sw4* __restrict__ a_mlfs,
		       float_sw4* __restrict__ a_unextf, float_sw4* __restrict__ a_bf, 
		       float_sw4* __restrict__ a_unextc, float_sw4* __restrict__ a_bc,
		       int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[], 
		       int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		       int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		       float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		       float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterpJacobi(float_sw4 rmax[3], Sarray &Uf, Sarray &UfJacobi, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
                      Sarray &Uc, Sarray &UcNew, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
                      Sarray &Mufs, Sarray &Mlfs,
                      Sarray &Unextf, Sarray &BfRestrict, Sarray &Unextc, Sarray &Bc,
                      int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
                      int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
                      float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
                      float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterpRO(float_sw4 rmax[3], Sarray &Uf, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
                      Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
                      Sarray &Mufs, Sarray &Mlfs,
                      Sarray &Unextf, Sarray &BfRestrict, Sarray &Unextc, Sarray &Bc,
                      int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
                      int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
                      float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
                      float_sw4 a_sbop[], float_sw4 a_ghcof[]);
   
void evenIevenJinterp(float_sw4 rmax[6], Sarray &Uf, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		     Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		     Sarray &Morc, Sarray &Mlrc,
		     Sarray &Unextf, Sarray &Bf, Sarray &Unextc, Sarray &Bc,
		     int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		     float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		      float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIevenJinterpJacobi(float_sw4 rmax[6], Sarray &Uf, Sarray &UfNew, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		     Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		     Sarray &Morc, Sarray &Mlrc,
		     Sarray &Unextf, Sarray &Bf, Sarray &Unextc, Sarray &Bc,
		     int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		     float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		     float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIevenJinterp(float_sw4 rmax[6], Sarray &Uf, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		     Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		     Sarray &Morc, Sarray &Mlrc,
		     Sarray &Unextf, Sarray &Bf, Sarray &Unextc, Sarray &Bc,
		     int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		     float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		     float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void evenIoddJinterp(float_sw4 rmax[6], Sarray &Uf, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		     Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		     Sarray &Morc, Sarray &Mlrc,
		     Sarray &Unextf, Sarray &Bf, Sarray &Unextc, Sarray &Bc,
		     int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		     float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		     float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void evenIoddJinterpJacobi(float_sw4 rmax[6], Sarray &Uf, Sarray &UfNext, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		     Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		     Sarray &Morc, Sarray &Mlrc,
		     Sarray &Unextf, Sarray &Bf, Sarray &Unextc, Sarray &Bc,
		     int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		     float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		     float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterp(float_sw4 rmax[3], Sarray &Uf, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		    Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		    Sarray &Mufs, Sarray &Mlfs,
		    Sarray &Unextf, Sarray &Bf, Sarray &Unextc, Sarray &Bc,
		    int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		    int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		    float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		    float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void update_unext( int ib, int ie, int jb, int je, int kb, int ke,
		   float_sw4* __restrict__ a_unext, float_sw4* __restrict__ a_up,
		   float_sw4* __restrict__ a_lutt, float_sw4* __restrict__ a_force,
		   float_sw4* __restrict__ a_rho, float_sw4 cof, int kic );

void dpdmt_wind( int ib, int ie, int jb, int je, int kb_tt, int ke_tt, int kb_u, int ke_u,
		 float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		 float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_utt,
		 float_sw4 dt2i );

void rhs4th3wind( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
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
			  float_sw4 h, float_sw4* __restrict__ a_strx, 
			  float_sw4* __restrict__ a_stry,  float_sw4* __restrict__ a_strz, 
			  char op );

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

   void addsg4wind_ci( float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* ,
		       float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* ,
		       float_sw4*, int, int, int, int, int, int, float_sw4, int, int, int, int );
   void ve_bndry_stress_curvi_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nz,
                              float_sw4 *alphap, float_sw4 *muve, float_sw4 *lave, float_sw4 *bforcerhs, 
			      float_sw4 *met, int side, float_sw4 *sbop, int usesg, float_sw4 *sgstrx,
			      float_sw4 *sgstry );
   void att_free_curvi_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                        int klast, float_sw4 *u, float_sw4 *mu, float_sw4 *la, float_sw4 *bforcerhs, 
			float_sw4 *met, float_sw4 *sbop, int usesg, float_sw4 *sgstrx, float_sw4 *sgstry );

#ifdef SW4_NOC
}
#endif

#endif
