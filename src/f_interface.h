#include "sw4.h"
extern "C" { // Fortran prototypes (to be removed once arguments are made equivalent)
   void tw_aniso_free_surf_z(int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                             int kz, float_sw4 t, float_sw4 om, float_sw4 cv, float_sw4 ph, float_sw4 omm,
			     float_sw4* phc, float_sw4* bforce, float_sw4 h, float_sw4 zmin );

   void twfrsurfz_wind( int *ifirst, int *ilast, int *jfirst, int *jlast, int *kfirst, int *klast,
                        float_sw4 *h, int* kz, float_sw4 *t, float_sw4 *omega, float_sw4 *c, float_sw4 *phase, 
			float_sw4 *bforce,
                        float_sw4 *mu, float_sw4 *lambda, float_sw4 *zmin,
                        int *i1, int *i2, int *j1, int *j2 );

   void twfrsurfz_att_wind( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                          float_sw4 h, int kz, float_sw4 t, float_sw4 omega, float_sw4 c, float_sw4 phase,
                          float_sw4 *bforce, float_sw4 *mu, float_sw4 *lambda, float_sw4 zmin,
                          int i1, int i2, int j1, int j2 );

   void twfrsurfzsg_wind( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                          float_sw4 h, int kz, float_sw4 t, float_sw4 omega, float_sw4 c, float_sw4 phase, 
			  float_sw4 omstrx, float_sw4 omstry,
                          float_sw4 *bforce, float_sw4 *mu, float_sw4 *lambda, float_sw4 zmin,
                          int i1, int i2, int j1, int j2 );

   void twfrsurfzsg_att_wind( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			      float_sw4 h, int kz, float_sw4 t, float_sw4 omega, float_sw4 c, float_sw4 phase, 
			      float_sw4 omstrx, float_sw4 omstry,
			      float_sw4 *bforce, float_sw4 *mu, float_sw4 *lambda, float_sw4 zmin,
			      int i1, int i2, int j1, int j2 );

   void satt(float_sw4 *up, float_sw4 *qs, float_sw4 *dt, float_sw4 *cfreq, int *ifirst, int *ilast, 
	     int *jfirst, int *jlast, int *kfirst, int *klast);


   void bcfort( int*, int*, int*, int*, int*, int*, 
		int *, int*, int*, int*,
		float_sw4*, float_sw4*, boundaryConditionType*, float_sw4 *, float_sw4*, float_sw4*, float_sw4*,
		float_sw4* bf0_p, float_sw4* bf1_p, 
		float_sw4* bf2_p, float_sw4*bf3_p, 
		float_sw4*bf4_p, float_sw4*bf5_p, 
		float_sw4*, float_sw4*, float_sw4*, int* );
   void freesurfcurvi(int*, int*, int*, int*, int*, int*, int*, int*,
		      float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );
   void freesurfcurvisg(int*, int*, int*, int*, int*, int*, int*, int*,
			float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
			float_sw4*, float_sw4*, float_sw4* );
   void bcfortsg( int*, int*, int*, int*, int*, int*, 
		  int *, int*, int*, int*,
		  float_sw4*, float_sw4*, boundaryConditionType*, float_sw4 *, float_sw4*, float_sw4*, float_sw4*,
		  float_sw4* bf0_p, float_sw4* bf1_p, 
		  float_sw4* bf2_p, float_sw4*bf3_p, 
		  float_sw4*bf4_p, float_sw4*bf5_p, 
		  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );

   void bcfortanisg( int*, int*, int*, int*, int*, int*,  int*, int*, int*, int*,
		     float_sw4*, float_sw4*, boundaryConditionType*, float_sw4*, float_sw4*,
		     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
		     float_sw4*, float_sw4* );

   void twfrsurfz(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p, int * kfirst_p, 
		  int * klast_p, float_sw4* h_p, int * k_p,
		  float_sw4* t_p, float_sw4* om_p, float_sw4* cv_p, float_sw4* ph_p,
		  float_sw4* bforce_side5_ptr, float_sw4* mu_ptr, float_sw4* la_ptr, float_sw4* zmin );
   void twfrsurfzatt(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p, int * kfirst_p, 
		     int * klast_p, float_sw4* h_p, int * k_p,
		     float_sw4* t_p, float_sw4* om_p, float_sw4* cv_p, float_sw4* ph_p,
		     float_sw4* bforce_side5_ptr, float_sw4* mu_ptr, float_sw4* la_ptr, float_sw4* zmin );
   void twfrsurfzsgstr(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p, int * kfirst_p, 
		       int * klast_p, float_sw4* h_p, int * k_p,
		       float_sw4* t_p, float_sw4* om_p, float_sw4* cv_p, float_sw4* ph_p,
		       float_sw4* omstrx_p, float_sw4* omstry_p,
		       float_sw4* bforce_side5_ptr, float_sw4* mu_ptr, float_sw4* la_ptr, float_sw4* zmin );
   void twfrsurfzsgstratt(int * ifirst_p, int * ilast_p, int * jfirst_p, int * jlast_p,
			  int * kfirst_p, int * klast_p, float_sw4* h_p, int * k_p,
			  float_sw4* t_p, float_sw4* om_p, float_sw4* cv_p, float_sw4* ph_p,
			  float_sw4* omstrx_p, float_sw4* omstry_p,
			  float_sw4* bforce_side5_ptr, float_sw4* mu_ptr, float_sw4* la_ptr, float_sw4* zmin );
   void memvarforcesurf( int*, int*, int*, int*, int*, float_sw4*, float_sw4*, float_sw4*, 
			 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );
   void memvarforcesurfc( int*, int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*, float_sw4*, 
			  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );
   void twdirbdry( int *wind_ptr, float_sw4 *h_p, float_sw4 *t_p, float_sw4 *om_p, float_sw4 * cv_p, 
		   float_sw4 *ph_p,  float_sw4 * bforce_side_ptr, float_sw4* zmin );
   void twdirbdryc( int* ifirst, int* ilast, int* jfirst, int* jlast, int* kfirst, int* klast,
		    int *wind_ptr, float_sw4 *t_p, float_sw4 *om_p, float_sw4 * cv_p, 
		    float_sw4 *ph_p,  float_sw4 * bforce_side_ptr, float_sw4* x, float_sw4* y, float_sw4* z );
   void testsrc( float_sw4* f_ptr, int* ifirst, int* ilast, int* jfirst, int* jlast, int* kfirst,
		 int* klast, int* nz, int* wind, float_sw4* m_zmin, float_sw4* h, int* kx, int* ky, int* kz,
		 float_sw4* momgrid );

   void addsg4wind( float_sw4*, float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* ,
		    float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* , float_sw4* ,
		    float_sw4*, int, int, int, int, int, int, float_sw4, int, int, int, int );

   void addsgd4(float_sw4* dt, float_sw4 *h, float_sw4 *a_Up, float_sw4*a_U, float_sw4*a_Um, float_sw4* Rho,
		float_sw4 *sg_dc_x, float_sw4* sg_dc_y, float_sw4* sg_dc_z, float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* sg_str_z,
		float_sw4* sg_corner_x, float_sw4* sg_corner_y, float_sw4* sg_corner_z,
		int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, float_sw4* damping_coefficient );
   void addsgd6(float_sw4* dt, float_sw4 *h, float_sw4 *a_Up, float_sw4*a_U, float_sw4*a_Um, float_sw4* Rho,
		float_sw4 *sg_dc_x, float_sw4* sg_dc_y, float_sw4* sg_dc_z, float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* sg_str_z,
		float_sw4* sg_corner_x, float_sw4* sg_corner_y, float_sw4* sg_corner_z,
		int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, float_sw4* damping_coefficient );
   void addsgd4c(float_sw4* dt, float_sw4 *a_Up, float_sw4*a_U, float_sw4*a_Um, float_sw4* Rho,
		 float_sw4 *sg_dc_x, float_sw4* sg_dc_y, float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* jac,
		 float_sw4* sg_corner_x, float_sw4* sg_corner_y, 
		 int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, float_sw4* damping_coefficient );
   void addsgd6c(float_sw4* dt, float_sw4 *a_Up, float_sw4*a_U, float_sw4*a_Um, float_sw4* Rho,
		 float_sw4 *sg_dc_x, float_sw4* sg_dc_y, float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* jac,
		 float_sw4* sg_corner_x, float_sw4* sg_corner_y, 
		 int *ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, float_sw4* damping_coefficient );
   //  subroutine RAYDIRBDRY( bforce, wind, t, lambda, mu, rho, cr, 
   // +     omega, alpha, h, zmin )
   void raydirbdry( float_sw4 *bforce_side_ptr, int *wind_ptr, float_sw4 *t, float_sw4 *lambda,
		    float_sw4 *mu, float_sw4 *rho,
		    float_sw4 *cr, float_sw4 *omega, float_sw4 *alpha, float_sw4 *h, float_sw4 *zmin );

   void twstensor( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* kz,
		   float_sw4* t, float_sw4* om, float_sw4* c, float_sw4* ph, float_sw4* xx, float_sw4* yy, float_sw4* zz,
		   float_sw4* tau, float_sw4* mu, float_sw4* lambda );
   void twstensoratt( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* kz,
		      float_sw4* t, float_sw4* om, float_sw4* c, float_sw4* ph, float_sw4* xx, float_sw4* yy, float_sw4* zz,
		      float_sw4* tau, float_sw4* mu, float_sw4* lambda );
   void twstensorsg( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* kz,
		     float_sw4* t, float_sw4* om, float_sw4* c, float_sw4* ph, float_sw4* xx, float_sw4* yy, float_sw4* zz,
		     float_sw4* tau, float_sw4* mu, float_sw4* lambda, float_sw4* omstrx, float_sw4* omstry );
   void twstensorsgatt( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, 
			int* kz, float_sw4* t, float_sw4* om, float_sw4* c, float_sw4* ph, float_sw4* xx, float_sw4* yy, 
			float_sw4* zz, float_sw4* tau, float_sw4* mu, float_sw4* lambda, float_sw4* omstrx,
			float_sw4* omstry );
   void getsurfforcing( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast,
			int* k, float_sw4* met, float_sw4* jac, float_sw4* tau, float_sw4* forcing );
   void subsurfforcing( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast,
			int* k, float_sw4* met, float_sw4* jac, float_sw4* tau, float_sw4* forcing );
   //   void getsurfforcinggh( int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast,
   //			  int* k, float_sw4* h, float_sw4* tau, float_sw4* forcing, float_sw4* amp, float_sw4* xc,
   //			  float_sw4* yc, float_sw4* xl, float_sw4* yl );
   void getsurfforcingsg( 
			 int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* k,
			 float_sw4* met, float_sw4* jac, float_sw4* tau, float_sw4* strx, float_sw4* stry,
			 float_sw4* forcing );
   void subsurfforcingsg( 
			 int*ifirst, int *ilast, int *jfirst, int* jlast, int* kfirst, int* klast, int* k,
			 float_sw4* met, float_sw4* jac, float_sw4* tau, float_sw4* strx, float_sw4* stry,
			 float_sw4* forcing );
   void addbstressc( int*, int*, int*, int*, int*, int*, int*, float_sw4*,  
		     float_sw4*,  float_sw4*,  float_sw4*, float_sw4*, int*,  float_sw4*, char*,
		     int*, int*,  float_sw4*,  float_sw4* );
   void addbstresswresc( int*, int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,    
			 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, int*, 
			 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
			 int*, float_sw4*, float_sw4* );
   void solveattfreec( int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,    
		       float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
		       float_sw4*, int*, float_sw4*, float_sw4* );
   void solveattfreeac( int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*, float_sw4*);
   void bcfreesurfcurvani( int*, int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*, int*,
			   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );
   void twilightfortwind( int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
			  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, int*, int*,
			  int*, int*, int*, int* );

   void forcingfort( int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, 
		     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );
   void forcingfortsg(int*, int*, int*, int*, int*, 
                      int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, 
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,float_sw4*,float_sw4*,float_sw4* );
   void forcingfortsgatt(int*, int*, int*, int*, int*, 
                      int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, 
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,float_sw4*,float_sw4*,float_sw4* );
   void forcingttfort( int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, 
		       float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4* );
   void att_free_curvi( int *ifirst, int *ilast, int *jfirst, int *jlast, int *kfirst,
                        int *klast, float_sw4 *u, float_sw4 *mu, float_sw4 *la, float_sw4 *bforce_rhs, 
			float_sw4 *met, float_sw4 *sbop, int *usesg, float_sw4 *sgstrx, float_sw4 *sgstry );

   void ve_bndry_stress_curvi(int *ifirst, int *ilast, int *jfirst, int *jlast, int *kfirst, int *klast, int *nz,
                              float_sw4 *alphap, float_sw4 *muve, float_sw4 *lave, float_sw4 *bforcerhs, 
			      float_sw4 *met, int *side, float_sw4 *sbop, int *usesg, float_sw4 *sgstrx,
			      float_sw4 *sgstry );
}
