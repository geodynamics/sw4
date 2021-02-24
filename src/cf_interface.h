#ifndef SW4_CF_INTERFACE_H
#include "Sarray.h"
#include "sw4.h"

// the extern "c" is only needed for linking the Fortran version of these
// routines
#ifdef SW4_NOC
extern "C" {
#endif

void update_unext(int ib, int ie, int jb, int je, int kb, int ke,
                  float_sw4* __restrict__ a_unext, float_sw4* __restrict__ a_up,
                  float_sw4* __restrict__ a_lutt,
                  float_sw4* __restrict__ a_force,
                  float_sw4* __restrict__ a_rho, float_sw4 cof, int kic);

void rhs4th3wind(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                 int klast, int nk, int* __restrict__ onesided,
                 float_sw4* __restrict__ a_acof, float_sw4* __restrict__ a_bope,
                 float_sw4* __restrict__ a_ghcof, float_sw4* __restrict__ a_lu,
                 float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_mu,
                 float_sw4* __restrict__ a_lambda, float_sw4 h,
                 float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
                 float_sw4* __restrict__ a_strz, char op, int kfirstu,
                 int klastu, int kfirstw, int klastw);
void rhs4th3wind_host(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nk,
    int* __restrict__ onesided, float_sw4* __restrict__ a_acof,
    float_sw4* __restrict__ a_bope, float_sw4* __restrict__ a_ghcof,
    float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, float_sw4 h,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
    float_sw4* __restrict__ a_strz, char op, int kfirstu, int klastu,
    int kfirstw, int klastw);

#ifdef SW4_NOC
}
#endif

void dpdmt_wind(int ib, int ie, int jb, int je, int kb_tt, int ke_tt, int kb_u,
                int ke_u, float_sw4* __restrict__ a_up,
                float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_um,
                float_sw4* __restrict__ a_utt, float_sw4 dt2i);

void evenIevenJinterpJacobi(float_sw4 rmax[6], Sarray& Uf, Sarray& UfNew,
                            Sarray& Uc, Sarray& Morc, Sarray& Mlrc,
                            Sarray& Morf, Sarray& Mlrf, Sarray& Unextf,
                            Sarray& UnextcInterp, int a_iStart[], int a_iEnd[],
                            int a_jStart[], int a_jEnd[], int a_kStart[],
                            int a_kEnd[], int a_iStartInt[], int a_iEndInt[],
                            int a_jStartInt[], int a_jEndInt[], int gf, int gc,
                            int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
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
    float_sw4* __restrict__ a_uncint, int a_iStart[], int a_iEnd[],
    int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[],
    int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
    int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
    float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterpJacobi(
    float_sw4 rmax[3], Sarray& Uf, Sarray& UfNew, Sarray& Uc, Sarray& UcNew,
    Sarray& Mufs, Sarray& Mlfs, Sarray& Morc, Sarray& Mlrc, Sarray& Mucs,
    Sarray& Mlcs, Sarray& Morf, Sarray& Mlrf, Sarray& Unextf,
    Sarray& BfRestrict, Sarray& Unextc, Sarray& Bc, int a_iStart[],
    int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[],
    int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
    int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
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
    int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[],
    int a_kEnd[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[],
    int a_jEndInt[], int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf,
    float_sw4 hc, float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[],
    float_sw4 a_ghcof[]);

void evenIevenJinterp(float_sw4 rmax[6], Sarray& Uf, Sarray& Muf,
                      Sarray& Lambdaf, Sarray& Rhof, Sarray& Uc, Sarray& Muc,
                      Sarray& Lambdac, Sarray& Rhoc, Sarray& Morc, Sarray& Mlrc,
                      Sarray& Unextf, Sarray& Bf, Sarray& Unextc, Sarray& Bc,
                      int a_iStart[], int a_jStart[], int a_iStartInt[],
                      int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
                      int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf,
                      float_sw4 hc, float_sw4 cof, float_sw4 relax,
                      float_sw4* a_strf_x, float_sw4* a_strf_y,
                      float_sw4* a_strc_x, float_sw4* a_strc_y,
                      float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIevenJinterpJacobi(float_sw4 rmax[6], Sarray& Uf, Sarray& UfNew,
                           Sarray& Uc, Sarray& Morc, Sarray& Mlrc, Sarray& Morf,
                           Sarray& Mlrf, Sarray& Unextf, Sarray& Unextc,
                           int a_iStart[], int a_iEnd[], int a_jStart[],
                           int a_jEnd[], int a_kStart[], int a_kEnd[],
                           int a_iStartInt[], int a_iEndInt[],
                           int a_jStartInt[], int a_jEndInt[], int gf, int gc,
                           int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
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
    int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[],
    int a_kEnd[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[],
    int a_jEndInt[], int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf,
    float_sw4 hc, float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[],
    float_sw4 a_ghcof[]);

void oddIevenJinterp(float_sw4 rmax[6], Sarray& Uf, Sarray& Muf,
                     Sarray& Lambdaf, Sarray& Rhof, Sarray& Uc, Sarray& Muc,
                     Sarray& Lambdac, Sarray& Rhoc, Sarray& Morc, Sarray& Mlrc,
                     Sarray& Unextf, Sarray& Bf, Sarray& Unextc, Sarray& Bc,
                     int a_iStart[], int a_jStart[], int a_iStartInt[],
                     int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
                     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf,
                     float_sw4 hc, float_sw4 cof, float_sw4 relax,
                     float_sw4* a_strf_x, float_sw4* a_strf_y,
                     float_sw4* a_strc_x, float_sw4* a_strc_y,
                     float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void evenIoddJinterp(float_sw4 rmax[6], Sarray& Uf, Sarray& Muf,
                     Sarray& Lambdaf, Sarray& Rhof, Sarray& Uc, Sarray& Muc,
                     Sarray& Lambdac, Sarray& Rhoc, Sarray& Morc, Sarray& Mlrc,
                     Sarray& Unextf, Sarray& Bf, Sarray& Unextc, Sarray& Bc,
                     int a_iStart[], int a_jStart[], int a_iStartInt[],
                     int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
                     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf,
                     float_sw4 hc, float_sw4 cof, float_sw4 relax,
                     float_sw4* a_strf_x, float_sw4* a_strf_y,
                     float_sw4* a_strc_x, float_sw4* a_strc_y,
                     float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void evenIoddJinterpJacobi(float_sw4 rmax[6], Sarray& Uf, Sarray& UfNext,
                           Sarray& Uc, Sarray& Morc, Sarray& Mlrc, Sarray& Morf,
                           Sarray& Mlrf, Sarray& Unextf, Sarray& Unextc,
                           int a_iStart[], int a_iEnd[], int a_jStart[],
                           int a_jEnd[], int a_kStart[], int a_kEnd[],
                           int a_iStartInt[], int a_iEndInt[],
                           int a_jStartInt[], int a_jEndInt[], int gf, int gc,
                           int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
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
    float_sw4* __restrict__ a_uncint, int a_iStart[], int a_iEnd[],
    int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[],
    int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
    int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
    float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void oddIoddJinterp(float_sw4 rmax[3], Sarray& Uf, Sarray& Muf, Sarray& Lambdaf,
                    Sarray& Rhof, Sarray& Uc, Sarray& Muc, Sarray& Lambdac,
                    Sarray& Rhoc, Sarray& Mufs, Sarray& Mlfs, Sarray& Unextf,
                    Sarray& Bf, Sarray& Unextc, Sarray& Bc, int a_iStart[],
                    int a_jStart[], int a_iStartInt[], int a_iEndInt[],
                    int a_jStartInt[], int a_jEndInt[], int gf, int gc, int nkf,
                    float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof,
                    float_sw4 relax, float_sw4* a_strf_x, float_sw4* a_strf_y,
                    float_sw4* a_strc_x, float_sw4* a_strc_y,
                    float_sw4 a_sbop[], float_sw4 a_ghcof[]);

void rhs4th3fortsgstr_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nk,
    int* __restrict__ onesided, float_sw4* __restrict__ a_acof,
    float_sw4* __restrict__ a_bope, float_sw4* __restrict__ a_ghcof,
    float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_mu, float_sw4* __restrict__ a_lambda, float_sw4 h,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
    float_sw4* __restrict__ a_strz, char op);

void rhs4th3fort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                    int klast, int nk, int* __restrict__ onesided,
                    float_sw4* __restrict__ a_acof,
                    float_sw4* __restrict__ a_bope,
                    float_sw4* __restrict__ a_ghcof,
                    float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_u,
                    float_sw4* __restrict__ a_mu,
                    float_sw4* __restrict__ a_lambda, float_sw4 h, char op);

void predfort_ci(int ib, int ie, int jb, int je, int kb, int ke,
                 float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                 float_sw4* __restrict__ um, float_sw4* __restrict__ lu,
                 float_sw4* __restrict__ fo, float_sw4* __restrict__ rho,
                 float_sw4 dt2);

void corrfort_ci(int ib, int ie, int jb, int je, int kb, int ke,
                 float_sw4* __restrict__ up, float_sw4* __restrict__ lu,
                 float_sw4* __restrict__ fo, float_sw4* __restrict__ rho,
                 float_sw4 dt4);

void dpdmtfort_ci(int ib, int ie, int jb, int je, int kb, int ke,
                  float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                  float_sw4* __restrict__ um, float_sw4* __restrict__ u2,
                  float_sw4 dt2i, int rank);

void rhouttlumf_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, int nz, float_sw4* __restrict__ a_uacc,
                   float_sw4* __restrict__ a_lu, float_sw4* __restrict__ a_fo,
                   float_sw4* __restrict__ a_rho, float_sw4 lowZ[3],
                   float_sw4 interZ[3], float_sw4 highZ[3]);

void rhserrfort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, int nz, float_sw4 h, float_sw4* __restrict__ a_fo,
                   float_sw4* __restrict__ a_u2, float_sw4 lowZ[3],
                   float_sw4 interZ[3], float_sw4 highZ[3]);

void dpdmtfortatt_ci(int ib, int ie, int jb, int je, int kb, int ke,
                     float_sw4* __restrict__ up, float_sw4* __restrict__ u,
                     float_sw4* __restrict__ um, float_sw4 dt2i);

void satt_ci(float_sw4* __restrict__ up, float_sw4* __restrict__ qs,
             float_sw4 dt, float_sw4 cfreq, int ifirst, int ilast, int jfirst,
             int jlast, int kfirst, int klast);

void solveattfreeac_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* __restrict__ a_alpha,
                       float_sw4 cof, float_sw4* __restrict__ a_up);

void solveattfreec_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_mu,
    float_sw4* __restrict__ a_la, float_sw4* __restrict__ a_muve,
    float_sw4* __restrict__ a_lave, float_sw4* __restrict__ a_bforcerhs,
    float_sw4* __restrict__ a_met, float_sw4 s[5], int usesg,
    float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry);

void addbstresswresc_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nz,
    float_sw4* __restrict__ a_alphap, float_sw4* __restrict__ a_alpham,
    float_sw4* __restrict__ a_muve, float_sw4* __restrict__ a_lave,
    float_sw4* __restrict__ a_bforcerhs, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_met, int side,
    float_sw4 dt, float_sw4 omegave, float_sw4* __restrict__ a_memforce,
    float_sw4* __restrict__ a_muvebnd, float_sw4* __restrict__ a_lambdavebnd,
    float_sw4 s[5], float_sw4& cof, int usesg, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry);

void addsg4wind_ci(float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, int, int, int, int, int,
                   int, float_sw4, int, int, int, int);
void ve_bndry_stress_curvi_ci(int ifirst, int ilast, int jfirst, int jlast,
                              int kfirst, int klast, int nz, float_sw4* alphap,
                              float_sw4* muve, float_sw4* lave,
                              float_sw4* bforcerhs, float_sw4* met, int side,
                              float_sw4* sbop, int usesg, float_sw4* sgstrx,
                              float_sw4* sgstry);
void att_free_curvi_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* u, float_sw4* mu, float_sw4* la,
                       float_sw4* bforcerhs, float_sw4* met, float_sw4* sbop,
                       int usesg, float_sw4* sgstrx, float_sw4* sgstry);

#endif
