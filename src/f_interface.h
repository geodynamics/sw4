#include "sw4.h"
extern "C" {  // Fortran prototypes (to be removed once arguments are made
              // equivalent)
void tw_aniso_free_surf_z(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                          sw4_type kfirst, sw4_type klast, sw4_type kz, float_sw4 t,
                          float_sw4 om, float_sw4 cv, float_sw4 ph,
                          float_sw4 omm, float_sw4* phc, float_sw4* bforce,
                          float_sw4 h, float_sw4 zmin);

void twfrsurfz_wind(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
                    sw4_type* kfirst, sw4_type* klast, float_sw4* h, sw4_type* kz,
                    float_sw4* t, float_sw4* omega, float_sw4* c,
                    float_sw4* phase, float_sw4* bforce, float_sw4* mu,
                    float_sw4* lambda, float_sw4* zmin, sw4_type* i1, sw4_type* i2,
                    sw4_type* j1, sw4_type* j2);

void twfrsurfz_att_wind(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                        sw4_type kfirst, sw4_type klast, float_sw4 h, sw4_type kz, float_sw4 t,
                        float_sw4 omega, float_sw4 c, float_sw4 phase,
                        float_sw4* bforce, float_sw4* mu, float_sw4* lambda,
                        float_sw4 zmin, sw4_type i1, sw4_type i2, sw4_type j1, sw4_type j2);

void twfrsurfzsg_wind(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst,
                      sw4_type klast, float_sw4 h, sw4_type kz, float_sw4 t,
                      float_sw4 omega, float_sw4 c, float_sw4 phase,
                      float_sw4 omstrx, float_sw4 omstry, float_sw4* bforce,
                      float_sw4* mu, float_sw4* lambda, float_sw4 zmin, sw4_type i1,
                      sw4_type i2, sw4_type j1, sw4_type j2);

void twfrsurfzsg_att_wind(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                          sw4_type kfirst, sw4_type klast, float_sw4 h, sw4_type kz,
                          float_sw4 t, float_sw4 omega, float_sw4 c,
                          float_sw4 phase, float_sw4 omstrx, float_sw4 omstry,
                          float_sw4* bforce, float_sw4* mu, float_sw4* lambda,
                          float_sw4 zmin, sw4_type i1, sw4_type i2, sw4_type j1, sw4_type j2);

void satt(float_sw4* up, float_sw4* qs, float_sw4* dt, float_sw4* cfreq,
          sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst,
          sw4_type* klast);

void bcfort(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*,
            float_sw4*, float_sw4*, boundaryConditionType*, float_sw4*,
            float_sw4*, float_sw4*, float_sw4*, float_sw4* bf0_p,
            float_sw4* bf1_p, float_sw4* bf2_p, float_sw4* bf3_p,
            float_sw4* bf4_p, float_sw4* bf5_p, float_sw4*, float_sw4*,
            float_sw4*, sw4_type*);
void freesurfcurvi(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void freesurfcurvisg(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*);
void bcfortsg(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*,
              float_sw4*, float_sw4*, boundaryConditionType*, float_sw4*,
              float_sw4*, float_sw4*, float_sw4*, float_sw4* bf0_p,
              float_sw4* bf1_p, float_sw4* bf2_p, float_sw4* bf3_p,
              float_sw4* bf4_p, float_sw4* bf5_p, float_sw4*, float_sw4*,
              float_sw4*, float_sw4*, float_sw4*);

void bcfortanisg(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*,
                 float_sw4*, float_sw4*, boundaryConditionType*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*);

void twfrsurfz(sw4_type* ifirst_p, sw4_type* ilast_p, sw4_type* jfirst_p, sw4_type* jlast_p,
               sw4_type* kfirst_p, sw4_type* klast_p, float_sw4* h_p, sw4_type* k_p,
               float_sw4* t_p, float_sw4* om_p, float_sw4* cv_p,
               float_sw4* ph_p, float_sw4* bforce_side5_ptr, float_sw4* mu_ptr,
               float_sw4* la_ptr, float_sw4* zmin);
void twfrsurfzatt(sw4_type* ifirst_p, sw4_type* ilast_p, sw4_type* jfirst_p, sw4_type* jlast_p,
                  sw4_type* kfirst_p, sw4_type* klast_p, float_sw4* h_p, sw4_type* k_p,
                  float_sw4* t_p, float_sw4* om_p, float_sw4* cv_p,
                  float_sw4* ph_p, float_sw4* bforce_side5_ptr,
                  float_sw4* mu_ptr, float_sw4* la_ptr, float_sw4* zmin);
void twfrsurfzsgstr(sw4_type* ifirst_p, sw4_type* ilast_p, sw4_type* jfirst_p, sw4_type* jlast_p,
                    sw4_type* kfirst_p, sw4_type* klast_p, float_sw4* h_p, sw4_type* k_p,
                    float_sw4* t_p, float_sw4* om_p, float_sw4* cv_p,
                    float_sw4* ph_p, float_sw4* omstrx_p, float_sw4* omstry_p,
                    float_sw4* bforce_side5_ptr, float_sw4* mu_ptr,
                    float_sw4* la_ptr, float_sw4* zmin);
void twfrsurfzsgstratt(sw4_type* ifirst_p, sw4_type* ilast_p, sw4_type* jfirst_p, sw4_type* jlast_p,
                       sw4_type* kfirst_p, sw4_type* klast_p, float_sw4* h_p, sw4_type* k_p,
                       float_sw4* t_p, float_sw4* om_p, float_sw4* cv_p,
                       float_sw4* ph_p, float_sw4* omstrx_p,
                       float_sw4* omstry_p, float_sw4* bforce_side5_ptr,
                       float_sw4* mu_ptr, float_sw4* la_ptr, float_sw4* zmin);
void memvarforcesurf(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*);
void memvarforcesurfc(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*);
void twdirbdry(sw4_type* wind_ptr, float_sw4* h_p, float_sw4* t_p, float_sw4* om_p,
               float_sw4* cv_p, float_sw4* ph_p, float_sw4* bforce_side_ptr,
               float_sw4* zmin);
void twdirbdryc(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst,
                sw4_type* klast, sw4_type* wind_ptr, float_sw4* t_p, float_sw4* om_p,
                float_sw4* cv_p, float_sw4* ph_p, float_sw4* bforce_side_ptr,
                float_sw4* x, float_sw4* y, float_sw4* z);
void testsrc(float_sw4* f_ptr, sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
             sw4_type* kfirst, sw4_type* klast, sw4_type* nz, sw4_type* wind, float_sw4* m_zmin,
             float_sw4* h, sw4_type* kx, sw4_type* ky, sw4_type* kz, float_sw4* momgrid);

void addsg4wind(float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, sw4_type,
                sw4_type, sw4_type, sw4_type, sw4_type, sw4_type, float_sw4, sw4_type, sw4_type, sw4_type, sw4_type);

void addsgd4(float_sw4* dt, float_sw4* h, float_sw4* a_Up, float_sw4* a_U,
             float_sw4* a_Um, float_sw4* Rho, float_sw4* sg_dc_x,
             float_sw4* sg_dc_y, float_sw4* sg_dc_z, float_sw4* sg_str_x,
             float_sw4* sg_str_y, float_sw4* sg_str_z, float_sw4* sg_corner_x,
             float_sw4* sg_corner_y, float_sw4* sg_corner_z, sw4_type* ifirst,
             sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst, sw4_type* klast,
             float_sw4* damping_coefficient);
void addsgd6(float_sw4* dt, float_sw4* h, float_sw4* a_Up, float_sw4* a_U,
             float_sw4* a_Um, float_sw4* Rho, float_sw4* sg_dc_x,
             float_sw4* sg_dc_y, float_sw4* sg_dc_z, float_sw4* sg_str_x,
             float_sw4* sg_str_y, float_sw4* sg_str_z, float_sw4* sg_corner_x,
             float_sw4* sg_corner_y, float_sw4* sg_corner_z, sw4_type* ifirst,
             sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst, sw4_type* klast,
             float_sw4* damping_coefficient);
void addsgd4c(float_sw4* dt, float_sw4* a_Up, float_sw4* a_U, float_sw4* a_Um,
              float_sw4* Rho, float_sw4* sg_dc_x, float_sw4* sg_dc_y,
              float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* jac,
              float_sw4* sg_corner_x, float_sw4* sg_corner_y, sw4_type* ifirst,
              sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst, sw4_type* klast,
              float_sw4* damping_coefficient);
void addsgd6c(float_sw4* dt, float_sw4* a_Up, float_sw4* a_U, float_sw4* a_Um,
              float_sw4* Rho, float_sw4* sg_dc_x, float_sw4* sg_dc_y,
              float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* jac,
              float_sw4* sg_corner_x, float_sw4* sg_corner_y, sw4_type* ifirst,
              sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst, sw4_type* klast,
              float_sw4* damping_coefficient);
//  subroutine RAYDIRBDRY( bforce, wind, t, lambda, mu, rho, cr,
// +     omega, alpha, h, zmin )
void raydirbdry(float_sw4* bforce_side_ptr, sw4_type* wind_ptr, float_sw4* t,
                float_sw4* lambda, float_sw4* mu, float_sw4* rho, float_sw4* cr,
                float_sw4* omega, float_sw4* alpha, float_sw4* h,
                float_sw4* zmin);

void twstensor(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst,
               sw4_type* klast, sw4_type* kz, float_sw4* t, float_sw4* om, float_sw4* c,
               float_sw4* ph, float_sw4* xx, float_sw4* yy, float_sw4* zz,
               float_sw4* tau, float_sw4* mu, float_sw4* lambda);
void twstensoratt(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst,
                  sw4_type* klast, sw4_type* kz, float_sw4* t, float_sw4* om,
                  float_sw4* c, float_sw4* ph, float_sw4* xx, float_sw4* yy,
                  float_sw4* zz, float_sw4* tau, float_sw4* mu,
                  float_sw4* lambda);
void twstensorsg(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast, sw4_type* kfirst,
                 sw4_type* klast, sw4_type* kz, float_sw4* t, float_sw4* om, float_sw4* c,
                 float_sw4* ph, float_sw4* xx, float_sw4* yy, float_sw4* zz,
                 float_sw4* tau, float_sw4* mu, float_sw4* lambda,
                 float_sw4* omstrx, float_sw4* omstry);
void twstensorsgatt(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
                    sw4_type* kfirst, sw4_type* klast, sw4_type* kz, float_sw4* t,
                    float_sw4* om, float_sw4* c, float_sw4* ph, float_sw4* xx,
                    float_sw4* yy, float_sw4* zz, float_sw4* tau, float_sw4* mu,
                    float_sw4* lambda, float_sw4* omstrx, float_sw4* omstry);
void getsurfforcing(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
                    sw4_type* kfirst, sw4_type* klast, sw4_type* k, float_sw4* met,
                    float_sw4* jac, float_sw4* tau, float_sw4* forcing);
void subsurfforcing(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
                    sw4_type* kfirst, sw4_type* klast, sw4_type* k, float_sw4* met,
                    float_sw4* jac, float_sw4* tau, float_sw4* forcing);
//   void getsurfforcinggh( sw4_type*ifirst, sw4_type *ilast, sw4_type *jfirst, sw4_type* jlast,
//   sw4_type* kfirst, sw4_type* klast,
//			  sw4_type* k, float_sw4* h, float_sw4* tau, float_sw4*
// forcing, float_sw4* amp, float_sw4* xc, 			  float_sw4* yc,
// float_sw4* xl, float_sw4* yl );
void getsurfforcingsg(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
                      sw4_type* kfirst, sw4_type* klast, sw4_type* k, float_sw4* met,
                      float_sw4* jac, float_sw4* tau, float_sw4* strx,
                      float_sw4* stry, float_sw4* forcing);
void subsurfforcingsg(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
                      sw4_type* kfirst, sw4_type* klast, sw4_type* k, float_sw4* met,
                      float_sw4* jac, float_sw4* tau, float_sw4* strx,
                      float_sw4* stry, float_sw4* forcing);
void addbstressc(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, sw4_type*,
                 float_sw4*, char*, sw4_type*, sw4_type*, float_sw4*, float_sw4*);
void addbstresswresc(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, sw4_type*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     sw4_type*, float_sw4*, float_sw4*);
void solveattfreec(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, sw4_type*, float_sw4*, float_sw4*);
void solveattfreeac(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*, float_sw4*,
                    float_sw4*);
void bcfreesurfcurvani(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*,
                       float_sw4*, sw4_type*, float_sw4*, float_sw4*, float_sw4*,
                       float_sw4*, float_sw4*);
void twilightfortwind(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*,
                      sw4_type*);

void forcingfort(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void forcingfortsg(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*);
void forcingfortsgatt(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*);
void forcingttfort(sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, sw4_type*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void att_free_curvi(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
                    sw4_type* kfirst, sw4_type* klast, float_sw4* u, float_sw4* mu,
                    float_sw4* la, float_sw4* bforce_rhs, float_sw4* met,
                    float_sw4* sbop, sw4_type* usesg, float_sw4* sgstrx,
                    float_sw4* sgstry);

void ve_bndry_stress_curvi(sw4_type* ifirst, sw4_type* ilast, sw4_type* jfirst, sw4_type* jlast,
                           sw4_type* kfirst, sw4_type* klast, sw4_type* nz, float_sw4* alphap,
                           float_sw4* muve, float_sw4* lave,
                           float_sw4* bforcerhs, float_sw4* met, sw4_type* side,
                           float_sw4* sbop, sw4_type* usesg, float_sw4* sgstrx,
                           float_sw4* sgstry);
}
