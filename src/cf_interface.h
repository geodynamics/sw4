#ifndef SW4_CF_INTERFACE_H
#include "Sarray.h"
#include "sw4.h"

// the extern "c" is only needed for linking the Fortran version of these
// routines
#ifdef SW4_NOC
extern "C" {
#endif

void update_unext(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                  float_sw4* __restrict__ a_unext, float_sw4* __restrict__ a_up,
                  float_sw4* __restrict__ a_lutt,
                  float_sw4* __restrict__ a_force,
                  float_sw4* __restrict__ a_rho, float_sw4 cof, sw4_type kic);

void rhs4th3wind(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst,
                 sw4_type klast, sw4_type nk, sw4_type* __restrict__ onesided,
                 float_sw4* __restrict__ a_acof, float_sw4* __restrict__ a_bope,
                 float_sw4* __restrict__ a_ghcof, float_sw4* __restrict__ a_lu,
                 float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_mu,
                 float_sw4* __restrict__ a_lambda, float_sw4 h,
                 float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
                 float_sw4* __restrict__ a_strz, char op, sw4_type kfirstu,
                 sw4_type klastu, sw4_type kfirstw, sw4_type klastw);
void rhs4th3wind_host(
    sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst, sw4_type klast, sw4_type nk,
    sw4_type* __restrict__ onesided, float_sw4* __restrict__ a_acof,
    float_sw4* __restrict__ a_bope, float_sw4* __restrict__ a_ghcof,
    float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, float_sw4 h,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
    float_sw4* __restrict__ a_strz, char op, sw4_type kfirstu, sw4_type klastu,
    sw4_type kfirstw, sw4_type klastw);

#ifdef SW4_NOC
}
#endif

void dpdmt_wind(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb_tt, sw4_type ke_tt, sw4_type kb_u,
                sw4_type ke_u, float_sw4* __restrict__ a_up,
                float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_um,
                float_sw4* __restrict__ a_utt, float_sw4 dt2i);

void evenIevenJinterpJacobi(float_sw4 rmax[6], Sarray& Uf, Sarray& UfNew,
                            Sarray& Uc, Sarray& Morc, Sarray& Mlrc,
                            Sarray& Morf, Sarray& Mlrf, Sarray& Unextf,
                            Sarray& UnextcInterp, sw4_type a_iStart[], sw4_type a_iEnd[],
                            sw4_type a_jStart[], sw4_type a_jEnd[], sw4_type a_kStart[],
                            sw4_type a_kEnd[], sw4_type a_iStartInt[], sw4_type a_iEndInt[],
                            sw4_type a_jStartInt[], sw4_type a_jEndInt[], sw4_type gf, sw4_type gc,
                            sw4_type nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
                            float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[],
                            float_sw4 a_ghcof[]);

void evenIevenJinterpJacobiOpt(
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax1,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax2,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax3,
    float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_ufnew,
    float_sw4* __restrict__ a_uc, float_sw4* __restrict__ a_morc,
    float_sw4* __restrict__ a_mlrc, float_sw4* __restrict__ a_morf,
    float_sw4* __restrict__ a_mlrf, float_sw4* __restrict__ a_unextf,
    float_sw4* __restrict__ a_uncint, sw4_type a_iStart[], sw4_type a_iEnd[],
    sw4_type a_jStart[], sw4_type a_jEnd[], sw4_type a_kStart[], sw4_type a_kEnd[],
    int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
    sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
    float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterpJacobi(
    float_sw4 rmax[3], Sarray& Uf, Sarray& UfNew, Sarray& Uc, Sarray& UcNew,
    Sarray& Mufs, Sarray& Mlfs, Sarray& Morc, Sarray& Mlrc, Sarray& Mucs,
    Sarray& Mlcs, Sarray& Morf, Sarray& Mlrf, Sarray& Unextf,
    Sarray& BfRestrict, Sarray& Unextc, Sarray& Bc, sw4_type a_iStart[],
    sw4_type a_iEnd[], sw4_type a_jStart[], sw4_type a_jEnd[], sw4_type a_kStart[], sw4_type a_kEnd[],
    sw4_type a_iStartInt[], sw4_type a_iEndInt[], sw4_type a_jStartInt[], sw4_type a_jEndInt[],
    sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
    float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterpJacobiOpt(
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax1,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax2,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax3,
    float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_ufnew,
    float_sw4* __restrict__ a_uc, float_sw4* __restrict__ a_ucnew,
    float_sw4* __restrict__ a_mufs, float_sw4* __restrict__ a_mlfs,
    float_sw4* __restrict__ a_morc, float_sw4* __restrict__ a_mlrc,
    float_sw4* __restrict__ a_mucs, float_sw4* __restrict__ a_mlcs,
    float_sw4* __restrict__ a_morf, float_sw4* __restrict__ a_mlrf,
    float_sw4* __restrict__ a_unextf, float_sw4* __restrict__ a_bfr,
    float_sw4* __restrict__ a_unextc, float_sw4* __restrict__ a_bc,
    sw4_type a_iStart[], sw4_type a_iEnd[], sw4_type a_jStart[], sw4_type a_jEnd[], sw4_type a_kStart[],
    sw4_type a_kEnd[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[],
    int a_jEndInt[], sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf,
    float_sw4 hc, float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[],
    float_sw4 a_ghcof[]);

void evenIevenJinterp(float_sw4 rmax[6], Sarray& Uf, Sarray& Muf,
                      Sarray& Lambdaf, Sarray& Rhof, Sarray& Uc, Sarray& Muc,
                      Sarray& Lambdac, Sarray& Rhoc, Sarray& Morc, Sarray& Mlrc,
                      Sarray& Unextf, Sarray& Bf, Sarray& Unextc, Sarray& Bc,
                      sw4_type a_iStart[], sw4_type a_jStart[], int a_iStartInt[],
                      int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
                      sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf,
                      float_sw4 hc, float_sw4 cof, float_sw4 relax,
                      float_sw4* a_strf_x, float_sw4* a_strf_y,
                      float_sw4* a_strc_x, float_sw4* a_strc_y,
                      float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIevenJinterpJacobi(float_sw4 rmax[6], Sarray& Uf, Sarray& UfNew,
                           Sarray& Uc, Sarray& Morc, Sarray& Mlrc, Sarray& Morf,
                           Sarray& Mlrf, Sarray& Unextf, Sarray& Unextc,
                           sw4_type a_iStart[], sw4_type a_iEnd[], sw4_type a_jStart[],
                           sw4_type a_jEnd[], sw4_type a_kStart[], sw4_type a_kEnd[],
                           int a_iStartInt[], int a_iEndInt[],
                           int a_jStartInt[], int a_jEndInt[], sw4_type gf, sw4_type gc,
                           sw4_type nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
                           float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[],
                           float_sw4 a_ghcof[]);

void oddIevenJinterpJacobiOpt(
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax1,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax2,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax3,
    float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_ufnew,
    float_sw4* __restrict__ a_uc, float_sw4* __restrict__ a_morc,
    float_sw4* __restrict__ a_mlrc, float_sw4* __restrict__ a_morf,
    float_sw4* __restrict__ a_mlrf, float_sw4* __restrict__ a_unextf,
    float_sw4* __restrict__ a_uncint,
    //			      Sarray &UnextcInterp,
    sw4_type a_iStart[], sw4_type a_iEnd[], sw4_type a_jStart[], sw4_type a_jEnd[], sw4_type a_kStart[],
    sw4_type a_kEnd[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[],
    int a_jEndInt[], sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf,
    float_sw4 hc, float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[],
    float_sw4 a_ghcof[]);

void oddIevenJinterp(float_sw4 rmax[6], Sarray& Uf, Sarray& Muf,
                     Sarray& Lambdaf, Sarray& Rhof, Sarray& Uc, Sarray& Muc,
                     Sarray& Lambdac, Sarray& Rhoc, Sarray& Morc, Sarray& Mlrc,
                     Sarray& Unextf, Sarray& Bf, Sarray& Unextc, Sarray& Bc,
                     sw4_type a_iStart[], sw4_type a_jStart[], sw4_type a_iStartInt[],
                     sw4_type a_iEndInt[], sw4_type a_jStartInt[], sw4_type a_jEndInt[],
                     sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf,
                     float_sw4 hc, float_sw4 cof, float_sw4 relax,
                     float_sw4* a_strf_x, float_sw4* a_strf_y,
                     float_sw4* a_strc_x, float_sw4* a_strc_y,
                     float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void evenIoddJinterp(float_sw4 rmax[6], Sarray& Uf, Sarray& Muf,
                     Sarray& Lambdaf, Sarray& Rhof, Sarray& Uc, Sarray& Muc,
                     Sarray& Lambdac, Sarray& Rhoc, Sarray& Morc, Sarray& Mlrc,
                     Sarray& Unextf, Sarray& Bf, Sarray& Unextc, Sarray& Bc,
                     sw4_type a_iStart[], sw4_type a_jStart[], sw4_type a_iStartInt[],
                     sw4_type a_iEndInt[], sw4_type a_jStartInt[], sw4_type a_jEndInt[],
                     sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf,
                     float_sw4 hc, float_sw4 cof, float_sw4 relax,
                     float_sw4* a_strf_x, float_sw4* a_strf_y,
                     float_sw4* a_strc_x, float_sw4* a_strc_y,
                     float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void evenIoddJinterpJacobi(float_sw4 rmax[6], Sarray& Uf, Sarray& UfNext,
                           Sarray& Uc, Sarray& Morc, Sarray& Mlrc, Sarray& Morf,
                           Sarray& Mlrf, Sarray& Unextf, Sarray& Unextc,
                           sw4_type a_iStart[], sw4_type a_iEnd[], sw4_type a_jStart[],
                           sw4_type a_jEnd[], sw4_type a_kStart[], sw4_type a_kEnd[],
                           sw4_type a_iStartInt[], sw4_type a_iEndInt[],
                           sw4_type a_jStartInt[], sw4_type a_jEndInt[], sw4_type gf, sw4_type gc,
                           sw4_type nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
                           float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[],
                           float_sw4 a_ghcof[]);

void evenIoddJinterpJacobiOpt(
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax1,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax2,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4>& rmax3,
    float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_ufnew,
    float_sw4* __restrict__ a_uc, float_sw4* __restrict__ a_morc,
    float_sw4* __restrict__ a_mlrc, float_sw4* __restrict__ a_morf,
    float_sw4* __restrict__ a_mlrf, float_sw4* __restrict__ a_unextf,
    float_sw4* __restrict__ a_uncint, sw4_type a_iStart[], sw4_type a_iEnd[],
    sw4_type a_jStart[], sw4_type a_jEnd[], sw4_type a_kStart[], sw4_type a_kEnd[],
    int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
    sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
    float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterp(float_sw4 rmax[3], Sarray& Uf, Sarray& Muf, Sarray& Lambdaf,
                    Sarray& Rhof, Sarray& Uc, Sarray& Muc, Sarray& Lambdac,
                    Sarray& Rhoc, Sarray& Mufs, Sarray& Mlfs, Sarray& Unextf,
                    Sarray& Bf, Sarray& Unextc, Sarray& Bc, sw4_type a_iStart[],
                    sw4_type a_jStart[], sw4_type a_iStartInt[], sw4_type a_iEndInt[],
                    sw4_type a_jStartInt[], sw4_type a_jEndInt[], sw4_type gf, sw4_type gc, sw4_type nkf,
                    float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof,
                    float_sw4 relax, float_sw4* a_strf_x, float_sw4* a_strf_y,
                    float_sw4* a_strc_x, float_sw4* a_strc_y,
                    float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void rhs4th3fortsgstr_ci(
    sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst, sw4_type klast, sw4_type nk,
    sw4_type* __restrict__ onesided, float_sw4* __restrict__ a_acof,
    float_sw4* __restrict__ a_bope, float_sw4* __restrict__ a_ghcof,
    float_sw4* __restrict__ a_lu, 
    float_sw4* __restrict__ a_u1,float_sw4* __restrict__ a_u2,float_sw4* __restrict__ a_u3,
    float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, float_sw4 h,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
    float_sw4* __restrict__ a_strz, char op);

void rhs4th3fort_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst,
                    sw4_type klast, sw4_type nk, sw4_type* __restrict__ onesided,
                    float_sw4* __restrict__ a_acof,
                    float_sw4* __restrict__ a_bope,
                    float_sw4* __restrict__ a_ghcof,
                    float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
                    float_sw4* __restrict__ a_mu,
                    float_sw4* __restrict__ a_lambda, float_sw4 h, char op);

void predfort_ci(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                 float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                 float_sw4* __restrict__ um, float_sw4* __restrict__ lu,
                 float_sw4* __restrict__ fo, float_sw4* __restrict__ rho,
                 float_sw4 dt2);

void corrfort_ci(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                 float_sw4* __restrict__ up, float_sw4* __restrict__ lu,
                 float_sw4* __restrict__ fo, float_sw4* __restrict__ rho,
                 float_sw4 dt4);

void dpdmtfort_ci(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                  float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                  float_sw4* __restrict__ um, float_sw4* __restrict__ u2,
                  float_sw4 dt2i, sw4_type rank);

void rhouttlumf_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst,
                   sw4_type klast, sw4_type nz, float_sw4* __restrict__ a_uacc,
                   float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_fo,
                   float_sw4* __restrict__ a_rho, float_sw4 lowZ[3],
                   float_sw4 sw4_typeerZ[3], float_sw4 highZ[3]);

void rhserrfort_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst,
                   sw4_type klast, sw4_type nz, float_sw4 h, float_sw4* __restrict__ a_fo,
                   float_sw4* __restrict__ a_u2, float_sw4 lowZ[3],
                   float_sw4 interZ[3], float_sw4 highZ[3]);

void dpdmtfortatt_ci(sw4_type ib, sw4_type ie, sw4_type jb, sw4_type je, sw4_type kb, sw4_type ke,
                     float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                     float_sw4* __restrict__ um, float_sw4 dt2i);

void satt_ci(float_sw4* __restrict__ up, float_sw4* __restrict__ qs,
             float_sw4 dt, float_sw4 cfreq, sw4_type ifirst, sw4_type ilast, sw4_type jfirst,
             sw4_type jlast, sw4_type kfirst, sw4_type klast);

void solveattfreeac_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst,
                       sw4_type klast, float_sw4* __restrict__ a_alpha,
                       float_sw4 cof, float_sw4* __restrict__ a_up);

void solveattfreec_ci(
    sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst, sw4_type klast,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_mu,
    float_sw4* __restrict__ a_la, float_sw4* __restrict__ a_muve,
    float_sw4* __restrict__ a_lave, float_sw4* __restrict__ a_bforcerhs,
    float_sw4* __restrict__ a_met, float_sw4 s[5], sw4_type usesg,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry);

void addbstresswresc_ci(
    sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst, sw4_type klast, sw4_type nz,
    float_sw4* __restrict__ a_alphap, float_sw4* __restrict__ a_alpham,
    float_sw4* __restrict__ a_muve, float_sw4* __restrict__ a_lave,
    float_sw4* __restrict__ a_bforcerhs, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_met, sw4_type side,
    float_sw4 dt, float_sw4 omegave, float_sw4* __restrict__ a_memforce,
    float_sw4* __restrict__ a_muvebnd, float_sw4* __restrict__ a_lambdavebnd,
    float_sw4 s[5], float_sw4& cof, sw4_type usesg, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry);

void addsg4wind_ci(float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, sw4_type, sw4_type, sw4_type, sw4_type, sw4_type,
                   sw4_type, float_sw4, sw4_type, sw4_type, sw4_type, sw4_type);
void ve_bndry_stress_curvi_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast,
                              sw4_type kfirst, sw4_type klast, sw4_type nz, float_sw4* alphap,
                              float_sw4* muve, float_sw4* lave,
                              float_sw4* bforcerhs, float_sw4* met, sw4_type side,
                              float_sw4* sbop, sw4_type usesg, float_sw4* sgstrx,
                              float_sw4* sgstry);
void att_free_curvi_ci(sw4_type ifirst, sw4_type ilast, sw4_type jfirst, sw4_type jlast, sw4_type kfirst,
                       sw4_type klast, float_sw4* u, float_sw4* mu, float_sw4* la,
                       float_sw4* bforcerhs, float_sw4* met, float_sw4* sbop,
                       sw4_type usesg, float_sw4* sgstrx, float_sw4* sgstry);

#endif
