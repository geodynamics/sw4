//  SW4 LICENSE #
// ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order #
// ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC.  #
// Produced at the Lawrence Livermore National Laboratory.  # #
// Written by: # N. Anders Petersson (petersson1@llnl.gov) # Bjorn
// Sjogreen (sjogreen2@llnl.gov) # # LLNL-CODE-643337 # # All rights
// reserved.  # # This file is part of SW4, Version: 1.0 # # Please
// also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License" # # This program is free software; you can
// redistribute it and/or modify # it under the terms of the GNU
// General Public License (as published by # the Free Software
// Foundation) version 2, dated June 1991.  # # This program is
// distributed in the hope that it will be useful, but # WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF # MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the terms and # conditions of
// the GNU General Public License for more details.  # # You should
// have received a copy of the GNU General Public License # along with
// this program; if not, write to the Free Software # Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>

#include "Byteswapper.h"
#include "EW.h"
#include "F77_FUNC.h"
#include "Mspace.h"
#include "caliper.h"
#include "cf_interface.h"
#include "mpi.h"
#include "policies.h"
#include "startEnd.h"
#include "version.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif
#if defined(SW4_EXPT_1)
// Experimental template based splitting along I,J,K
#include "curvilinear4sgc.h"
#endif
#if defined(SW4_EXPT_3)
#include "curvilinear4sgcX3.h"
#endif
extern __constant__ double cmem_acof[384];
extern __constant__ double cmem_acof_no_gp[384];
extern "C" {
void tw_aniso_force(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                    int klast, float_sw4* fo, float_sw4 t, float_sw4 om,
                    float_sw4 cv, float_sw4 ph, float_sw4 omm, float_sw4 phm,
                    float_sw4 amprho, float_sw4* phc, float_sw4 h,
                    float_sw4 zmin);

void tw_aniso_curvi_force(int ifirst, int ilast, int jfirst, int jlast,
                          int kfirst, int klast, float_sw4* fo, float_sw4 t,
                          float_sw4 om, float_sw4 cv, float_sw4 ph,
                          float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                          float_sw4* phc, float_sw4* xx, float_sw4* yy,
                          float_sw4* zz);

void tw_aniso_force_tt(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* fo, float_sw4 t, float_sw4 om,
                       float_sw4 cv, float_sw4 ph, float_sw4 omm, float_sw4 phm,
                       float_sw4 amprho, float_sw4* phc, float_sw4 h,
                       float_sw4 zmin);

void tw_aniso_curvi_force_tt(int ifirst, int ilast, int jfirst, int jlast,
                             int kfirst, int klast, float_sw4* fo, float_sw4 t,
                             float_sw4 om, float_sw4 cv, float_sw4 ph,
                             float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                             float_sw4* phc, float_sw4* xx, float_sw4* yy,
                             float_sw4* zz);

void corrfort(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
              float_sw4*, float_sw4*, float_sw4*);
void dpdmtfort(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
               float_sw4*, float_sw4*, float_sw4*);
void predfort(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
              float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void rhouttlumf(int*, int*, int*, int*, int*, int*, int*, float_sw4*,
                float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                float_sw4*);

void forcingfort(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void forcingfortc(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                  float_sw4*);
void forcingfortatt(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void forcingfortattc(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*);
void forcingfortsg(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*);
void forcingfortcsg(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void forcingfortsgatt(int*, int*, int*, int*, int*, int*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*);
void forcingfortsgattc(int*, int*, int*, int*, int*, int*, float_sw4*,
                       float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                       float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                       float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                       float_sw4*, float_sw4*, float_sw4*);
void forcingttfortsg(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*);
void forcingttfortcsg(int*, int*, int*, int*, int*, int*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*);
void forcingttfortsgatt(int*, int*, int*, int*, int*, int*, float_sw4*,
                        float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                        float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                        float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                        float_sw4*, float_sw4*);
void forcingttfortsgattc(int*, int*, int*, int*, int*, int*, float_sw4*,
                         float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                         float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                         float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                         float_sw4*, float_sw4*, float_sw4*);
void forcingttfort(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void forcingttfortc(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                    float_sw4*);
void forcingttattfort(int*, int*, int*, int*, int*, int*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*);
void forcingttattfortc(int*, int*, int*, int*, int*, int*, float_sw4*,
                       float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                       float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                       float_sw4*, float_sw4*, float_sw4*, float_sw4*);
//   void addmemvarforcing( int*, int*, int*, int*, int*, int*, float_sw4*,
//   float_sw4*, float_sw4*,
//			  float_sw4*, float_sw4*, float_sw4*, float_sw4*,
// float_sw4*, float_sw4* );
//   void addmemvarforcingc( int*, int*, int*, int*, int*, int*, float_sw4*,
//   float_sw4*, float_sw4*,
//			   float_sw4*, float_sw4*, float_sw4*, float_sw4*,
// float_sw4*, float_sw4*, float_sw4* );
void exactaccfort(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void exactaccfortc(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*);
void rhserrfort(int*, int*, int*, int*, int*, int*, int*, float_sw4*,
                float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void rhs4th3fort(int*, int*, int*, int*, int*, int*, int*, int*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*, char*);
void rhs4th3fortsgstr(int*, int*, int*, int*, int*, int*, int*, int*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, char*);
void exactrhsfort(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void exactrhsfortc(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*);
void exactrhsfortsg(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*);
void exactrhsfortsgc(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void solerr3(int*, int*, int*, int*, int*, int*, float_sw4* h, float_sw4* uex,
             float_sw4* u, float_sw4* li, float_sw4* l2, float_sw4* xli,
             float_sw4* zmin, float_sw4* x0, float_sw4* y0, float_sw4* z0,
             float_sw4* radius, int* imin, int* imax, int* jmin, int* jmax,
             int* kmin, int* kmax);
void solerr3c(int*, int*, int*, int*, int*, int*, float_sw4* uex, float_sw4* u,
              float_sw4* x, float_sw4* y, float_sw4* z, float_sw4* jac,
              float_sw4* li, float_sw4* l2, float_sw4* xli, float_sw4* x0,
              float_sw4* y0, float_sw4* z0, float_sw4* radius, int* imin,
              int* imax, int* jmin, int* jmax, int* kmin, int* kmax, int* usesg,
              float_sw4* strx, float_sw4* stry);
void solerrgp(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
              float_sw4*, float_sw4* li, float_sw4* l2);
void twilightfort(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void twilightfortc(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                   float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                   float_sw4*);
void twilightfortatt(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*);
void twilightfortattc(int*, int*, int*, int*, int*, int*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*);
//  subroutine rayleighfort( ifirst, ilast, jfirst, jlast, kfirst, klast,
// +     u, t, lambda, mu, rho, cr, omega, alpha, h, zmin )
void rayleighfort(int* ifirst, int* ilast, int* jfirst, int* jlast, int* kfirst,
                  int* klast, double* u, double* t, double* lambda, double* mu,
                  double* rho, double* cr, double* omega, double* alpha,
                  double* h, double* zmin);
void velsum(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
            int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
            float_sw4*);
void energy4(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
             int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
             float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void energy4c(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
              int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
              float_sw4*, float_sw4*);
void lambexact(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
               float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
               float_sw4*, int*);
void curvilinear4(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, int*,
                  float_sw4*, float_sw4*, float_sw4*, char*);
void curvilinear4sg(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, int*,
                    float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                    char*);

void addgradrho(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
                int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                int*);
void addgradrhoc(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
                 int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                 int*);
void addgradmula(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
                 int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, int*, int*,
                 int*, float_sw4*);
void addgradmulac(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
                  int*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                  float_sw4*, int*, int*, int*, float_sw4*);
#ifdef ENABLE_OPT
void F77_FUNC(projectmtrlc, PROJECTMTRLC)(int*, int*, int*, int*, int*, int*,
                                          int*, int*, int*, int*, int*, int*,
                                          double*, double*, double*, double*,
                                          double*, double*, double*, double*,
                                          double*, double*, double*, int*);
void F77_FUNC(projectmtrl, PROJECTMTRL)(int*, int*, int*, int*, int*, int*,
                                        int*, int*, int*, int*, int*, int*,
                                        double*, double*, double*, double*,
                                        double*, double*, double*, double*,
                                        double*, double*, int*);
void F77_FUNC(checkmtrl, CHECKMTRL)(int*, int*, int*, int*, int*, int*, double*,
                                    double*, double*, double*, double*,
                                    double*);
#endif
void exactmatfortatt(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                     float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                     float_sw4*);
void exactmatfortattc(int*, int*, int*, int*, int*, int*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                      float_sw4*, float_sw4*, float_sw4*, float_sw4*);
void updatememvar(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                  int*);

void dpdmtfortatt(int*, int*, int*, int*, int*, int*, float_sw4*, float_sw4*,
                  float_sw4*, float_sw4*);

void scalar_prod(int, int, int, int, int, int, int, int, int, int, int, int,
                 int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                 float_sw4*, float_sw4*);

void F77_FUNC(dgels, DGELS)(char& TRANS, int& M, int& N, int& NRHS, double* A,
                            int& LDA, double* B, int& LDB, double* WORK,
                            int& LWORK, int& INFO);

void innerloopanisgstrvc(int*, int*, int*, int*, int*, int*, int*, float_sw4*,
                         float_sw4*, float_sw4*, int*, float_sw4*, float_sw4*,
                         float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                         float_sw4*);

void ilanisocurv(int*, int*, int*, int*, int*, int*, int*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, int*, float_sw4*,
                 float_sw4*, float_sw4*, float_sw4*, float_sw4*, float_sw4*);

void memvar_pred_fort(int, int, int, int, int, int, double*, double*, double*,
                      double, double, int);

void memvar_corr_fort(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                      int klast, double* alp, double* alm, double* up,
                      double* u, double* um, double omega, double dt,
                      int domain);

void memvar_corr_fort_wind(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, double* alp, int d1b, int d1e,
                           int d2b, int d2e, int d3b, int d3e, double* alm,
                           double* up, double* u, double* um, double omega,
                           double dt, int domain);
}

void ilanisocurv_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nk,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_c,
    float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_lu, int* onesided,
    float_sw4* __restrict__ a_acof, float_sw4* __restrict__ a_bope,
    float_sw4* __restrict__ a_ghcof, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry, float_sw4* __restrict__ a_strz);

void innerloopanisgstrvc_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast, int nk,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_lu,
    float_sw4* __restrict__ a_c, int* onesided, float_sw4* __restrict__ a_acof,
    float_sw4* __restrict__ a_bope, float_sw4* __restrict__ a_ghcof,
    float_sw4 h, float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry,
    float_sw4* __restrict__ a_strz);

// the routine will replace the Fortran routine curvilinear4sg()
void curvilinear4sg_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_mu,
    float_sw4* __restrict__ a_lambda, float_sw4* __restrict__ a_met,
    float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_lu, int* onesided,
    float_sw4* __restrict__ a_acof, float_sw4* __restrict__ a_bope,
    float_sw4* __restrict__ a_ghcof, float_sw4* __restrict__ a_acof_no_gp,
    float_sw4* __restrict__ a_ghcof_no_gp, float_sw4* __restrict__ a_strx,
    float_sw4* __restrict__ a_stry, int nk, char op);
void energy4_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                int klast, int i1, int i2, int j1, int j2, int k1, int k2,
                int* onesided, float_sw4* __restrict__ a_um,
                float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_up,
                float_sw4* __restrict__ a_rho, float_sw4 h, float_sw4* a_strx,
                float_sw4* a_stry, float_sw4* a_strz, float_sw4& a_energy);
void energy4c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                 int klast, int i1, int i2, int j1, int j2, int k1, int k2,
                 int* onesided, float_sw4* __restrict__ a_um,
                 float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_up,
                 float_sw4* __restrict__ a_rho, float_sw4* __restrict__ a_jac,
                 float_sw4& a_energy);
void addgradrho_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, int ifirstact, int ilastact, int jfirstact,
                   int jlastact, int kfirstact, int klastact,
                   float_sw4* __restrict__ a_kap,
                   float_sw4* __restrict__ a_kapacc,
                   float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_u,
                   float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_uacc,
                   float_sw4* __restrict__ a_grho, float_sw4 dt, float_sw4 h,
                   int onesided[6]);
void addgradrhoc_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                    int klast, int ifirstact, int ilastact, int jfirstact,
                    int jlastact, int kfirstact, int klastact,
                    float_sw4* __restrict__ a_kap,
                    float_sw4* __restrict__ a_kapacc,
                    float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_u,
                    float_sw4* __restrict__ a_up,
                    float_sw4* __restrict__ a_uacc,
                    float_sw4* __restrict__ a_grho, float_sw4 dt,
                    float_sw4* __restrict__ a_jac, int onesided[6]);
void addgradmula_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                    int klast, int ifirstact, int ilastact, int jfirstact,
                    int jlastact, int kfirstact, int klastact,
                    float_sw4* __restrict__ a_kap,
                    float_sw4* __restrict__ a_kapacc,
                    float_sw4* __restrict__ a_u, float_sw4* __restrict__ a_uacc,
                    float_sw4* __restrict__ a_gmu,
                    float_sw4* __restrict__ a_glambda, float_sw4 dt,
                    float_sw4 h, int onesided[6], int nb, int wb,
                    float_sw4* __restrict__ a_bop);
void addgradmulac_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
    int ifirstact, int ilastact, int jfirstact, int jlastact, int kfirstact,
    int klastact, float_sw4* __restrict__ a_kap,
    float_sw4* __restrict__ a_kapacc, float_sw4* __restrict__ a_u,
    float_sw4* __restrict__ a_uacc, float_sw4* __restrict__ a_gmu,
    float_sw4* __restrict__ a_glambda, float_sw4 dt, float_sw4 h,
    float_sw4* __restrict__ a_met, float_sw4* __restrict__ a_jac,
    int onesided[6], int nb, int wb, float_sw4* __restrict__ a_bop);

void scalar_prod_ci(int is, int ie, int js, int je, int ks, int ke, int i1,
                    int i2, int j1, int j2, int k1, int k2, int onesided[6],
                    float_sw4* a_u, float_sw4* a_v, float_sw4* a_strx,
                    float_sw4* a_stry, float_sw4* a_strz, float_sw4& sc_prod);

void memvar_pred_fort_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4* alp, float_sw4* alm,
                         float_sw4* u, float_sw4 omega, float_sw4 dt,
                         int domain);

void memvar_corr_fort_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4* alp, float_sw4* alm,
                         float_sw4* up, float_sw4* u, float_sw4* um,
                         float_sw4 omega, float_sw4 dt, int domain);

void memvar_corr_fort_wind_ci(int ifirst, int ilast, int jfirst, int jlast,
                              int kfirst, int klast, float_sw4* alp, int d1b,
                              int d1e, int d2b, int d2e, int d3b, int d3e,
                              float_sw4* alm, float_sw4* up, float_sw4* u,
                              float_sw4* um, float_sw4 omega, float_sw4 dt,
                              int domain);

void addMemVarPredCart(double zMin, double h, double t, Sarray& alpha,
                       double omegaVE, double dt, double omega, double phase,
                       double c);

void addMemVarPredCurvilinear(Sarray& a_X, Sarray& a_Y, Sarray& a_Z, double t,
                              Sarray& alpha, double omegaVE, double dt,
                              double omega, double phase, double c);

void addMemVarCorr2Cart(double zMin, double h, double t, Sarray& alpha,
                        double omegaVE, double dt, double omega, double phase,
                        double c);

void addMemVarCorr2Curvilinear(Sarray& a_X, Sarray& a_Y, Sarray& a_Z, double t,
                               Sarray& alpha, double omegaVE, double dt,
                               double omega, double phase, double c);

using namespace std;

#define SQR(x) ((x) * (x))

// constructor
EW::EW(const string& fileName, vector<vector<Source*>>& a_GlobalSources,
       vector<vector<TimeSeries*>>& a_GlobalTimeSeries, bool a_invproblem)
    : m_epi_lat(0.0),
      m_epi_lon(0.0),
      m_epi_depth(0.0),
      m_epi_t0(0.0),
      m_topo_zmax(0.0),
      m_topoInputStyle(UNDEFINED),
      mTopoImageFound(false),
      m_nx_base(0),
      m_ny_base(0),
      m_nz_base(0),
      m_h_base(0.0),
      mSourcesOK(false),
      mIsInitialized(false),
      mParsingSuccessful(false),
      mNumberOfGrids(0),
      mName(fileName),
      m_scenario(" "),
      // mPath("./"),
      //      mObsPath("./"),
      mTempPath("./"),
      mWriteGMTOutput(false),
      mPlotFrequency(80),
      mNumFiles(0),
      mVerbose(0),
      mQuiet(false),
      mDebugIO(false),
      mHomogeneous(false),
      m_iotiming(false),
      m_pfs(false),
      m_nwriters(8),
      mTimeIsSet(false),
      mTmax(0.0),
      mTstart(0.0),
      mDt(0.0),
      // mNumberOfTimeSteps(-1),
      m_testing(false),
      m_moment_test(false),
      m_twilight_forcing(NULL),
      m_point_source_test(0),
      m_energy_test(0),
      m_lamb_test(0),
      m_rayleigh_wave_test(0),
      m_update_boundary_function(0),
      m_maxIter(10),
      m_topoFileName("NONE"),
      m_topoExtFileName("NONE"),
      m_QueryType("MAXRES"),
      //  mTestingEnergy(false),
      //  mTestSource(false),
      //  mTestLamb(false),
      mOrder(4),
      mCFL(1.3),
      mCFLmax(1.3),
      // m_d4coeff(0.0),
      // m_d4_cfl(0.2),
      // m_curlcoeff(0.0),

      // mRestartFilePrefix(""),
      // mRestartFromCycle(0),
      // mRestartDumpInterval(0),

      m_doubly_periodic(false),
      mbcsSet(false),

      m_analytical_topo(false),
      m_use_analytical_metric(false),
      m_GaussianAmp(0.05),
      m_GaussianLx(0.15),
      m_GaussianLy(0.15),
      m_GaussianXc(0.5),
      m_GaussianYc(0.5),

      m_use_supergrid(false),
      m_sg_gp_thickness(30),
      m_supergrid_damping_coefficient(
          0.02),  // good value for 4th order diss. Must be reduced by factor of
                  // 4 for 6th order diss.
      m_sg_damping_order(4),
      m_use_sg_width(
          false),  // width in meters instead of thickness in grid points.
      m_minJacobian(0.),
      m_maxJacobian(0.),

      m_energy_log(false),
      m_energy_print(false),
      m_energy_logfile("energy.dat"),
      m_saved_energy(0.0),

      m_do_timing(false),
      m_timing_print(0),
      m_output_detailed_timing(false),
      m_output_load(false),

      m_projection_cycle(1000),
      m_checkfornan(false),

      m_error_log_file("TwilightErr.txt"),
      m_error_log(false),
      m_error_print(true),
      m_inner_loop(9),
      m_topography_exists(false),
      m_useVelocityThresholds(false),
      m_vpMin(0.),
      m_vsMin(0.),
      m_grid_interpolation_order(3),
      m_zetaBreak(0.95),
      m_global_xmax(0.),
      m_global_ymax(0.),
      m_global_zmax(0.),
      m_global_zmin(0.),
      m_ghost_points(2),  // for 4th order stencils
      m_ppadding(2),
      m_ext_ghost_points(0),  // extra width for 6th order
                              // discretization of metric at a source
      //  m_ghost_points(3), // for 6th order stencils
      //  m_ppadding(3),

      mLonOrigin(-118.0),  // NTS
      mLatOrigin(37.0),    // NTS
      mGeoAz(0.0),         // x=North, y=East
      //  mDefaultLocation(true),
      mMetersPerDegree(111319.5),    // per degree in Latitude...
      mMetersPerLongitude(87721.0),  // approximate for Lat=38 degrees
      mConstMetersPerLongitude(false),

      // command limitfrequency
      m_limit_frequency(false),
      m_ppw(15),
      m_frequency_limit(1e38),  // will hold min(Vs/h)/PPW

      // command prefilter
      m_prefilter_sources(false),
      m_filter_observations(false),
      m_filter_ptr(0),
      m_filterobs_ptr(0),

      mPrintInterval(100),
      m_matrices_decomposed(false),
      m_citol(1e-7),
      m_cimaxiter(20),
      m_cirelfact(0.94),
      m_mesh_refinements(false),
      //  m_intp_conservative(true),
      mMaterialExtrapolate(0),

      m_use_attenuation(false),
      m_number_mechanisms(0),
      m_velo_omega(-1.0),
      m_min_omega(-1.0),
      m_max_omega(-1.0),
      m_geodynbc_found(false),
      m_geodynbc_center(false),
      m_do_geodynbc(false),
      //  m_do_geodynbc(false),
      m_att_use_max_frequency(false),
      m_att_ppw(8.0),
      m_inverse_problem(a_invproblem),
      m_maxit(0),
      m_maxrestart(0),
      m_iniguess_pos(false),
      m_iniguess_t0fr(false),
      m_iniguess_mom(false),
      m_iniguess_shifts(false),
      m_output_initial_seismograms(false),
      m_compute_scalefactors(false),
      m_cgstepselection(0),
      m_cgvarcase(0),
      m_cgfletcherreeves(true),
      m_do_linesearch(true),
      //  m_utc0set(false),
      //  m_utc0isrefevent(false),
      m_events_parallel(false),
      m_opttest(0),
      m_perturb(0),
      m_iperturb(1),
      m_jperturb(1),
      m_kperturb(1),
      m_pervar(1),
      m_qmultiplier(1),
      m_randomize(false),
      m_randomize_density(false),
      m_anisotropic(false),
      m_croutines(true),
      NO_TOPO(1e38),
      ForceVector(NULL),
      ForceAddress(NULL),
      ProfilerOn(false) {
  MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &m_nProcs);

  if (sizeof(float_sw4) == 4)
    m_mpifloat = MPI_FLOAT;
  else if (sizeof(float_sw4) == 8)
    m_mpifloat = MPI_DOUBLE;
  else
    CHECK_INPUT(false, "Error, could not identify float_sw4");

    // The memory space in which MPI buffers will be allocated
    // for use by AMPISendRecv. Options are Managed,Device and Pinned
    // On Sierra Host will work as well. Pinned is the fastest:
    // 256 node case with 4 GPS per node
    // Pinned 1m52s
    // Managed 3m 21s
    // Device 2m 15s
    // jsrun & lrun need additional flags like -M "-gpu "for Device and
    // Space::Managed to work.

#if defined(SW4_DEVICE_MPI_BUFFERS)
  mpi_buffer_space = Space::Device;
  if (!m_myRank) std::cout << "Using MPI buffers in device memory\n";
  //if (!mpi_supports_device_buffers()) {
   // std::cerr << "SW4 must be run using the -M -gpu flag with Device buffers\n";
    //abort();
  //}
#elif defined(SW4_MANAGED_MPI_BUFFERS)
  mpi_buffer_space = Space::Managed;
  if (!m_myRank) std::cout << "Using MPI buffers in managed memory\n";
  //if (!mpi_supports_device_buffers()) {
   // std::cerr
    //    << "SW4 must be run using the -M -gpu flag with Managed buffers\n";
    //abort();
  //}
#elif defined(SW4_PINNED_MPI_BUFFERS)
  mpi_buffer_space = Space::Pinned;
  if (!m_myRank)
    std::cout << "Using MPI buffers in pinned memory(COMPILE OPTION)\n";
#elif defined(SW4_STAGED_MPI_BUFFERS)
  mpi_buffer_space = Space::Pinned;
  if (!m_myRank) std::cout << "Using staged MPI buffers (COMPILE OPTION)\n";
#else
  mpi_buffer_space = Space::Pinned;
  if (!m_myRank) std::cout << "Using MPI buffers in pinned memory(DEFAULT)\nn";
#endif

  m_check_point = new CheckPoint(this);

#ifdef SW4_NOC
  m_croutines = false;
#endif
  Sarray::m_corder = m_croutines;

  //   m_error_checking = new ErrorChecking();
  // initialize the boundary condition array
  for (int i = 0; i < 6; i++) {
    mbcGlobalType[i] = bNone;
  }

  // No scaling is default
  for (int i = 0; i < 11; i++) m_scalefactors[i] = 1;

  m_proc_array[0] = 0;
  m_proc_array[1] = 0;

  m_nevents_specified = findNumberOfEvents();
  m_nevent = m_nevents_specified > 0 ? m_nevents_specified : 1;
  // std::cout<<"EVENTS "<<m_nevents_specified<<"  "<<m_nevent<<"\n";
  // Allocate storage
  m_epi_lat.resize(m_nevent);
  m_epi_lon.resize(m_nevent);
  m_epi_depth.resize(m_nevent);
  m_epi_t0.resize(m_nevent);
  a_GlobalSources.resize(m_nevent);
  a_GlobalTimeSeries.resize(m_nevent);
  mPath.resize(m_nevent);
  mObsPath.resize(m_nevent);
  mTmax.resize(m_nevent);
  mNumberOfTimeSteps.resize(m_nevent);
  mTimeIsSet.resize(m_nevent);
  m_utc0.resize(m_nevent);
  // Defaults
  for (int e = 0; e < m_nevent; e++) {
    m_epi_lat[e] = 0.0;
    m_epi_lon[e] = 0.0;
    m_epi_depth[e] = 0.0;
    m_epi_t0[e] = 0.0;
    //      mPath[e] = "./";
    //      mObsPath[e] = "./";
    mTmax[e] = 0.0;
    mNumberOfTimeSteps[e] = -1;
    mTimeIsSet[e] = false;
    m_utc0[e].resize(7);
  }

  MPI_Comm_dup(MPI_COMM_WORLD, &m_1d_communicator);

  // read the input file and setup the simulation object
  if (parseInputFile(a_GlobalSources, a_GlobalTimeSeries))
    mParsingSuccessful = true;

    // AP: need to figure out a better way of handling these error log files
    //
    // char fname[100];
    // sprintf(fname,"sw4-error-log-p%i.txt", m_myRank);
    // msgStream.open(fname);

#if defined(ENABLE_GPU)
  float_sw4* tmpa =
      SW4_NEW(Space::Managed, float_sw4[6 + 384 + 24 + 48 + 6 + 384 + 6 + 6]);
  m_sbop = tmpa;  // PTR_PUSH(Space::Managed,m_sbop);
  m_acof = m_sbop + 6;
  PTR_PUSH(Space::Managed, m_acof, 384 * sizeof(float_sw4));
  m_bop = m_acof + 384;
  PTR_PUSH(Space::Managed, m_bop, 24 * sizeof(float_sw4));
  m_bope = m_bop + 24;
  PTR_PUSH(Space::Managed, m_bope, 48 * sizeof(float_sw4));
  m_ghcof = m_bope + 48;
  PTR_PUSH(Space::Managed, m_ghcof, 6 * sizeof(float_sw4));

  m_acof_no_gp = m_ghcof + 6;
  PTR_PUSH(Space::Managed, m_acof_no_gp, 384 * sizeof(float_sw4));
  m_ghcof_no_gp = m_acof_no_gp + 384;
  PTR_PUSH(Space::Managed, m_ghcof_no_gp, 6 * sizeof(float_sw4));
  m_sbop_no_gp = m_ghcof_no_gp + 6;
  PTR_PUSH(Space::Managed, m_sbop_no_gp, 6 * sizeof(float_sw4));
#endif
}

// Destructor
EW::~EW() {
  // std::cout<<"EW::~EW() ...\n"<<std::flush;
#if defined(ENABLE_GPU)
  ::operator delete[](m_sbop, Space::Managed);
#endif
  ::operator delete[](viewArrayActual, Space::Managed);
  for (int m = 0; m < mNumberOfGrids; m += 4) {
  //std::cout<<"MPI BUFFER TYPE IS "<<as_int(mpi_buffer_space)<<" "<<std::get<0>(bufs_type1[4 * m])<<" \n";
    ::operator delete[](std::get<0>(bufs_type1[4 * m]), mpi_buffer_space);
    ::operator delete[](std::get<0>(bufs_type3[4 * m]), mpi_buffer_space);
    ::operator delete[](std::get<0>(bufs_type4[4 * m]), mpi_buffer_space);
    ::operator delete[](std::get<0>(bufs_type21[4 * m]), mpi_buffer_space);
  }

  for (int m = 0; m < mNumberOfCartesianGrids; m++) {
    ::operator delete[](std::get<0>(bufs_type_2dx[m]), Space::Pinned);
    ::operator delete[](std::get<0>(bufs_type_2dy[m]), Space::Pinned);
  }
  ::operator delete[](ForceVector, Space::Managed);
  ::operator delete[](ForceAddress, Space::Managed);

  ::operator delete[](global_variables.device_buffer, Space::Managed_temps);
#ifndef SW4_USE_UMPIRE
  // Delete the allocations stored in Sarray::static_map
  for (auto& i : Sarray::static_map) {
    ::operator delete[](i.second, Space::Managed);
  }
#endif

#if defined(SW4_TRACK_MPI)
  stringstream filename, hfilename, bfilename;
  filename << "MpiStats" << m_myRank;
  ofstream ofile(filename.str());
  hfilename << "TimeStepHistory" << m_myRank;
  ofstream hfile(hfilename.str());
  bfilename << "BarrierHistory" << m_myRank;
  ofstream bfile(bfilename.str());
  // ofile<<"# Size KB Bandwidth GB/s Callcount\n";
  // for ( auto it : mpi_times){
  //   ofile<<it.first*8/1024.0<<"
  //   "<<it.first*mpi_count[it.first]/it.second*8*1.0e6/1024/1024/1024<<"
  //   "<<mpi_count[it.first]<<"\n";
  // }
  // ofile<<"# AMPI_Sendrecv2\n";
  // for ( auto it : mpi_times2){
  //   ofile<<it.first*8/1024.0<<"
  //   "<<it.first*mpi_count2[it.first]/it.second*8*1.0e6/1024/1024/1024<<"
  //   "<<mpi_count2[it.first]<<"\n";
  // }
  // sm.print(ofile);
  sm.print(
      ofile, [=](size_t size) -> double { return size * 8 / 1024.0; },
      [=](size_t size, double time) -> double {
        return size / time * 8 * 1.0e6 / 1024 / 1024 / 1024;
      },
      "Sendrecv");
  sm2.print(
      ofile, [=](size_t size) -> double { return size * 8 / 1024.0; },
      [=](size_t size, double time) -> double {
        return size / time * 8 * 1.0e6 / 1024 / 1024 / 1024;
      },
      "SendRecv2");
  coll_sm.printhistory(bfile);  // This needs to be done before print in which
                                // the history gets sorted
  coll_sm.print(
      ofile, [=](int size) -> double { return size; },
      [=](int size, double time) -> double { return time; }, "Barrier");
  step_sm.printhistory(hfile);  // This needs to be done before print in which
                                // the history gets sorted
  host_sm.print(
      ofile, [=](size_t size) -> double { return size * 8 / 1024.0; },
      [=](size_t size, double time) -> double {
        return size / time * 8 * 1.0e6 / 1024 / 1024 / 1024;
      },
      "Host");

  step_sm.print(
      ofile, [=](size_t size) -> double { return size; },
      [=](size_t size, double time) -> double { return time; }, "STEP_ms");

  ofile.close();
  hfile.close();
#endif
  // std::cout<<"EW::~EW() DONE\n"<<std::flush;
  //  msgStream.close();
}

//-----------------------------------
bool EW::isInitialized() { return (mIsInitialized && mSourcesOK); }

//-----------------------------------
bool EW::wasParsingSuccessful() { return mParsingSuccessful; }

//-----------------------------------------------------------------------
void EW::printTime(int cycle, float_sw4 t, bool force) const {
  if (!mQuiet && proc_zero() &&
      (force || mPrintInterval == 1 || (cycle % mPrintInterval) == 1 ||
       cycle == 1)) {
    time_t now;
    time(&now);
    // string big enough for >1 million time steps
    printf("Time step %7i  t = %15.7e\t%s", cycle, t, ctime(&now));
    fflush(stdout);
  }
}
//-----------------------------------------------------------------------
void EW::printPreamble(vector<Source*>& a_Sources, int event) const {
  stringstream msg;

  if (!mQuiet && proc_zero()) {
    msg << "============================================================"
        << endl
        << " Running program on " << m_nProcs << " MPI tasks"
        << " using the following data: " << endl
        << endl
        << " Start Time = " << mTstart << " Goal Time = ";

    if (mTimeIsSet[event])
      msg << mTmax[event] << endl;
    else
      msg << mNumberOfTimeSteps[event] * mDt << endl;

    msg << " Number of time steps = " << mNumberOfTimeSteps[event]
        << " dt: " << mDt << endl;

    if (mVerbose) {
      msg << endl;
      msg << "============================================================"
          << endl;
      msg << " Global boundary conditions " << endl;
      const char* side_names[6] = {"x=0   ", "x=xMax", "y=0   ",
                                   "y=yMax", "z=topo", "z=zMax"};
      for (int side = 0; side < 6; side++) {
        msg << "      ";
        msg << side_names[side] << " " << bc_name(mbcGlobalType[side]) << "\n";
      }
      msg << endl;

      if (mHomogeneous)
        msg << " Assuming Mu and Lambda to be constant within each grid  "
            << endl;

      //         if (mForcing == 1 || mForcing == 2 || mForcing == 5)
      //            msg << endl << " Second order Dirichlet boundary condition,
      //            gamma=" << mEBDirichletRegCoeff << endl;
      //         else if (mForcing == 3 || mForcing == 4 || mForcing == 6)
      //            msg << endl << " Second order Neumann boundary condition" <<
      //            endl;

      if (mVerbose >= 4)
        cout << " The following point sources are used: " << endl;
    }
    cout << msg.str();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  cout.flush();
  cerr.flush();

  // m0 values in each source command gets added up. This number is called the
  // "Total seismic moment" and should be printed to stdout with the unit Nm
  // (Newton-meter). If that number is >0, you should also print Mw = 2/3
  // *(log10(M0) - 9.1). That is the moment magnitude (dimensionless).

  if (proc_zero()) {
    if (m_twilight_forcing) {
      cout << "-----------------------------------------------------------"
           << endl;
      cout << "Twilight zone testing (aka method of manufactured solution)"
           << endl;
      cout << "Parameters:" << endl;
      cout << "  omega = " << m_twilight_forcing->m_omega << endl;
      cout << "  c = " << m_twilight_forcing->m_c << endl;
      cout << "  phase = " << m_twilight_forcing->m_phase << endl;
      cout << "  mat-omega = " << m_twilight_forcing->m_momega << endl;
      cout << "  mat-phase = " << m_twilight_forcing->m_mphase << endl;
      cout << "  amprho = " << m_twilight_forcing->m_amprho << endl;
      cout << "  amplambda = " << m_twilight_forcing->m_amplambda << endl;
      cout << "  ampmu = " << m_twilight_forcing->m_ampmu << endl;
      cout << "-----------------------------------------------------------"
           << endl;
    } else if (m_lamb_test) {
      float_sw4 fx, fy, fz, xs, ys, zs;
      a_Sources[0]->getForces(fx, fy, fz);
      xs = a_Sources[0]->getX0();
      ys = a_Sources[0]->getY0();
      zs = a_Sources[0]->getZ0();
      string tfun;
      if (a_Sources[0]->getTfunc() == iVerySmoothBump)
        tfun = "VerySmoothBump";
      else if (a_Sources[0]->getTfunc() == iC6SmoothBump)
        tfun = "C6SmoothBump";

      cout << "-----------------------------------------------------------"
           << endl;
      cout << "Lamb's problem testing" << endl;
      cout << "Parameters:" << endl;
      cout << "  Cp = " << m_lamb_test->m_cp << endl;
      cout << "  Cs = " << m_lamb_test->m_cs << endl;
      cout << "  Rho = " << m_lamb_test->m_rho << endl;
      cout << "  (xs, ys, zs) = " << xs << ", " << ys << ", " << zs << endl;
      cout << "  (fx, fy, fz) = " << fx << ", " << fy << ", " << fz << endl;
      cout << "  Source time fcn = " << tfun << endl;
      cout << "-----------------------------------------------------------"
           << endl;
    } else {
      float_sw4 myM0Sum = 0;
      int numsrc = 0;  //, ignoredSources=0;
#pragma omp parallel for reduction(+ : numsrc, myM0Sum)
      for (unsigned int i = 0; i < a_Sources.size(); ++i) {
        if (a_Sources[i]->isMomentSource()) {
          numsrc++;
          myM0Sum += a_Sources[i]->getAmplitude();
          // if (i == 1)
          //   std::cout << " SOURE " << i << " " <<
          //   a_Sources[i]->getAmplitude()
          //             << "\n";
        }
      }
      if (!mQuiet) {
        stringstream msg2;
        msg2 << endl
             << "--------------------------------------------------------------"
                "---------"
             << endl
             << "  Total seismic moment (M0): " << myM0Sum << " Nm " << endl;
        if (myM0Sum > 0)
          msg2 << "  Moment magnitude     (Mw): "
               << (2. / 3.) * (log10(myM0Sum) - 9.1) << endl;
        msg2 << "  Number of moment sources " << numsrc << endl;
        msg2 << "--------------------------------------------------------------"
                "---------"
             << endl;
        cout << msg2.str();
      }
    }  // standard run
  }    // end if proc_zero()
}

//-----------------------------------------------------------------------
void EW::switch_on_checkfornan() { m_checkfornan = true; }

//-----------------------------------------------------------------------
void EW::assign_local_bcs() {
  SW4_MARK_FUNCTION;
  // This routine assigns m_bcType[g][b], b=0,1,2,3, based on mbcGlobalType,
  // taking parallel overlap boundaries into account

  int g, b, side;
  int top = mNumberOfGrids -
            1;  // index of the top grid in the arrays m_iStart, m_iEnd, etc

  // horizontal bc's are the same for all grids
  for (g = 0; g < mNumberOfGrids; g++) {
    // start by copying the global bc's
    for (b = 0; b <= 3; b++) m_bcType[g][b] = mbcGlobalType[b];

    // printf("assign_local_bc> BEFORE loop: rank=%d, bct[0]=%d, bct[1]=%d,
    // bct[2]=%d, bct[3]=%d\n", m_myRank,
    //        m_bcType[g][0], m_bcType[g][1], m_bcType[g][2], m_bcType[g][3]);

    if (m_iStart[top] + m_ghost_points > 1) {
      m_bcType[g][0] = bProcessor;
    }
    if (m_iEnd[top] - m_ghost_points < m_global_nx[top]) {
      m_bcType[g][1] = bProcessor;
    }

    // for a periodic domain, we need to change all bc to bProcessor if more
    // than 1 process is used in that direction
    if (m_doubly_periodic && m_proc_array[0] > 1) {
      m_bcType[g][0] = bProcessor;
      m_bcType[g][1] = bProcessor;
    }

    if (m_jStart[top] + m_ghost_points > 1) {
      //       printf(" rank=%d, jStart=%d, setting bcType[2] = proc\n",
      //       m_myRank, m_jStart[top]);
      m_bcType[g][2] = bProcessor;
    }
    if (m_jEnd[top] - m_ghost_points < m_global_ny[top]) {
      //       printf(" rank=%d, jEnd=%d, setting bcType[3] = proc\n", m_myRank,
      //       m_jEnd[top]);
      m_bcType[g][3] = bProcessor;
    }

    // for a periodic domain, we need to change all bc to bProcessor if more
    // than 1 process is used in that direction
    if (m_doubly_periodic && m_proc_array[1] > 1) {
      m_bcType[g][2] = bProcessor;
      m_bcType[g][3] = bProcessor;
    }

    // printf("assign_local_bc> AFTER loop: rank=%d, bct[0]=%d, bct[1]=%d,
    // bct[2]=%d, bct[3]=%d\n", m_myRank,
    //        m_bcType[g][0], m_bcType[g][1], m_bcType[g][2], m_bcType[g][3]);
  }

  // vertical bc's are interpolating except at the bottom and the top, where
  // they equal the global conditions
  //   ( Only preliminary support for acoustic/elastic, not fully implemented)
  m_bcType[top][4] = mbcGlobalType[4];
  for (g = 0; g < mNumberOfGrids - 1; g++) {
    if (m_iscurvilinear[g + 1] && !m_iscurvilinear[g])  // Elastic case only
      m_bcType[g][4] = bCCInterface;
    if (!m_iscurvilinear[g + 1] &&
        !m_iscurvilinear[g])  // Two Cartesian grids, must be refinement bndry.
      m_bcType[g][4] = bRefInterface;
    if (!m_iscurvilinear[g + 1] && m_iscurvilinear[g])  // Acoustic case only
      m_bcType[g][4] = bCCInterface;
    if (m_iscurvilinear[g + 1] &&
        m_iscurvilinear[g])  // Acoustic/Elastic interface
      m_bcType[g][4] = bAEInterface;
  }

  m_bcType[0][5] = mbcGlobalType[5];
  for (g = 1; g < mNumberOfGrids; g++) {
    if (m_iscurvilinear[g] && !m_iscurvilinear[g - 1])  // Elastic case
      m_bcType[g][5] = bCCInterface;
    if (!m_iscurvilinear[g] &&
        !m_iscurvilinear[g -
                         1])  // Two Cartesian grids, must be refinement bndry.
      m_bcType[g][5] = bRefInterface;
    if (!m_iscurvilinear[g] && m_iscurvilinear[g - 1])  // Acoustic case
      m_bcType[g][5] = bCCInterface;
    if (m_iscurvilinear[g] &&
        m_iscurvilinear[g - 1])  // Acoustic/Elastic interface
      m_bcType[g][5] = bAEInterface;
  }

  // Find out which boundaries need one sided approximation in mixed derivatives
  int ncurv = mNumberOfGrids - mNumberOfCartesianGrids;
  for (g = 0; g < mNumberOfGrids; g++) {
    for (side = 0; side < 4; side++) m_onesided[g][side] = 0;
    for (side = 4; side < 6; side++)
      //        m_onesided[g][side] = (m_bcType[g][side] == bStressFree) ||
      //        (m_bcType[g][side]== bCCInterface) ;
      m_onesided[g][side] = (m_bcType[g][side] == bStressFree) ||
                            (m_bcType[g][side] == bRefInterface) ||
                            (m_bcType[g][side] == bAEInterface) ||
                            (m_bcType[g][side] == bCCInterface &&
                             !(m_gridGenerator->curviCartIsSmooth(ncurv)));
    // Pre CURVIMR configuration
    // for (side = 4; side < 6; side++)
    //   m_onesided[g][side] = (m_bcType[g][side] == bStressFree) ||
    //                         (m_bcType[g][side] == bRefInterface) ||
    //                         (m_bcType[g][side] == bAEInterface);
  }
}

//-----------------------------------------------------------------------
// Note that the padding cell array is no longer needed.
// use m_iStartInt[g], m_iEndInt[g] to get the range of interior points
void EW::initializePaddingCells() {
  SW4_MARK_FUNCTION;
  int g = mNumberOfGrids - 1;

  for (int aa = 0; aa < 4; aa++) {
    if (m_bcType[g][aa] == bProcessor) {
      m_paddingCells[aa] = m_ppadding;
    } else {
      m_paddingCells[aa] = m_ghost_points;
    }
  }
}

//-----------------------------------------------------------------------
void EW::check_dimensions() {
  SW4_MARK_FUNCTION;
  for (int g = 0; g < mNumberOfGrids; g++) {
    int nz = m_kEndInt[g] - m_kStartInt[g] + 1;
    int nzmin;
    if (m_onesided[g][4] && m_onesided[g][5])
      nzmin = 12;
    else if (m_onesided[g][4] || m_onesided[g][5])
      nzmin = 8;
    else
      nzmin = 1;
    REQUIRE2(nz >= nzmin,
             "The number of grid points (not counting ghost pts) in the "
             "z-direction in grid "
                 << g << " must be >= " << nzmin << " current value is " << nz);
    int nx = m_iEndInt[g] - m_iStartInt[g] + 1;
    REQUIRE2(nx >= 1,
             "No grid points left (not counting ghost pts) in the x-direction "
             "in grid "
                 << g);
    int ny = m_jEndInt[g] - m_jStartInt[g] + 1;
    REQUIRE2(ny >= 1,
             "No grid points left (not counting ghost pts) in the y-direction "
             "in grid "
                 << g);
  }
}

//-----------------------------------------------------------------------
bool EW::proc_zero() const { return (m_myRank == 0); }

//-----------------------------------------------------------------------
int EW::no_of_procs() const { return m_nProcs; }

//-----------------------------------------------------------------------
string EW::bc_name(const boundaryConditionType bc) const {
  string retval;
  if (bc == bStressFree)
    retval = "free surface";
  else if (bc == bDirichlet)
    retval = "dirichlet";
  else if (bc == bSuperGrid)
    retval = "supergrid";
  else if (bc == bPeriodic)
    retval = "periodic";
  else if (bc == bCCInterface)
    retval = "Curvilinear/Cartesian interface";
  else if (bc == bRefInterface)
    retval = "Grid refinement interface";
  else if (bc == bAEInterface)
    retval = "Acoustic/Elastic interface";
  else if (bc == bProcessor)
    retval = "processor";
  else if (bc == bNone)
    retval = "none";
  return retval;
}

//-----------------------------------------------------------------------
bool EW::getDepth(float_sw4 x, float_sw4 y, float_sw4 z, float_sw4& depth) {
  // get the depth below the free surface
  bool success = false;

  if (!topographyExists()) {
    depth = z;
    success = true;
  } else {
    // topography
    float_sw4 zMinTilde;
    int gCurv = mNumberOfGrids - 1;
    float_sw4 h = mGridSize[gCurv];
    float_sw4 q = x / h + 1.0;
    float_sw4 r = y / h + 1.0;

    // define the depth for ghost points (in x or y) to equal the depth on the
    // nearest boundary point
    float_sw4 qMin = 1.0;
    float_sw4 qMax = (float_sw4)m_global_nx[gCurv];
    float_sw4 rMin = 1.0;
    float_sw4 rMax = (float_sw4)m_global_ny[gCurv];

    if (q < qMin) q = qMin;
    if (q > qMax) q = qMax;
    if (r < rMin) r = rMin;
    if (r > rMax) r = rMax;

    // // evaluate elevation of topography on the grid (smoothed topo)
    success = true;
    int ret = m_gridGenerator->interpolate_topography(this, x, y, zMinTilde,
                                                      mTopoGridExt);
    if (ret < 0) {
      cerr << "ERROR: getDepth: Unable to evaluate topography for x=" << x
           << " y= " << y << " on proc # " << getRank() << ", ret=" << ret
           << endl;
      //            cerr << "q=" << q << " r=" << r << " qMin=" << qMin << "
      //            qMax=" << qMax << " rMin=" << rMin << " rMax=" << rMax <<
      //            endl;
      // cerr << "Setting elevation of topography to ZERO" << endl;
      success = false;
      //      zMinTilde = 0;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    depth = z - zMinTilde;
  }
  return success;
}

//-----------------------------------------------------------------------
void EW::computeCartesianCoordGMG(double& x, double& y, double lon, double lat,
                                  char* crs_to) {
  m_geoproj->computeCartesianCoordGMG(x, y, lon, lat, crs_to);
}

//-----------------------------------------------------------------------
void EW::computeCartesianCoord(double& x, double& y, double lon, double lat) {
  SW4_MARK_FUNCTION;
  // -----------------------------------------------------------------
  // Compute the cartesian coordinate given the geographic coordinate
  // -----------------------------------------------------------------
  if (m_geoproj == 0)
  //  // compute x and y
  {
    double deg2rad = M_PI / 180.0;
    double phi = mGeoAz * deg2rad;
    //     x = mMetersPerDegree*(cos(phi)*(lat-mLatOrigin) +
    //     cos(lat*deg2rad)*(lon-mLonOrigin)*sin(phi)); y =
    //     mMetersPerDegree*(-sin(phi)*(lat-mLatOrigin) +
    //     cos(lat*deg2rad)*(lon-mLonOrigin)*cos(phi));
    if (mConstMetersPerLongitude) {
      x = mMetersPerDegree * cos(phi) * (lat - mLatOrigin) +
          mMetersPerLongitude * (lon - mLonOrigin) * sin(phi);
      y = mMetersPerDegree * (-sin(phi)) * (lat - mLatOrigin) +
          mMetersPerLongitude * (lon - mLonOrigin) * cos(phi);
    } else {
      x = mMetersPerDegree *
          (cos(phi) * (lat - mLatOrigin) +
           cos(lat * deg2rad) * (lon - mLonOrigin) * sin(phi));
      y = mMetersPerDegree *
          (-sin(phi) * (lat - mLatOrigin) +
           cos(lat * deg2rad) * (lon - mLonOrigin) * cos(phi));
    }
  } else
    m_geoproj->computeCartesianCoord(x, y, lon, lat);
  // Test with proj4
  //  projUV lonlat, xy;
  //  lonlat.u = lon*deg2rad;
  //  lonlat.v = lat*deg2rad;
  //  xy = pj_fwd( lonlat, m_projection );
  //  if( xy.u == HUGE_VAL )
  //     cout << "ERROR: computeCartesianCoord pj_fwd failed with message " <<
  //     pj_strerrno(pj_errno) << endl;
  //  xy.u -= m_xoffset;
  //  xy.v -= m_yoffset;
  //  x =  xy.u*sin(phi) + cos(phi)*xy.v;
  //  y =  xy.u*cos(phi) - sin(phi)*xy.v;
}

//-----------------------------------------------------------------------
void EW::computeGeographicCoord(double x, double y, double& longitude,
                                double& latitude) {
  // conversion factor between degrees and radians
  if (m_geoproj == 0) {
    double deg2rad = M_PI / 180.0;
    double phi = mGeoAz * deg2rad;
    // Compute the latitude
    latitude = mLatOrigin + (x * cos(phi) - y * sin(phi)) / mMetersPerDegree;
    // Compute the longitude
    if (mConstMetersPerLongitude) {
      longitude =
          mLonOrigin + (x * sin(phi) + y * cos(phi)) / (mMetersPerLongitude);
    } else {
      longitude = mLonOrigin + (x * sin(phi) + y * cos(phi)) /
                                   (mMetersPerDegree * cos(latitude * deg2rad));
    }
  } else
    m_geoproj->computeGeographicCoord(x, y, longitude, latitude);
  // Test with proj4
  // projUV lonlat, xy;
  // xy.u = x*sin(phi) + y*cos(phi) + m_xoffset;
  // xy.v = x*cos(phi) - y*sin(phi) + m_yoffset;
  // lonlat = pj_inv( xy, m_projection );
  //  if( lonlat.u == HUGE_VAL )
  //     cout << "ERROR: computeGeographicCoord, pj_inv failed with message " <<
  //     pj_strerrno(pj_errno) << endl;
  // longitude = lonlat.u/deg2rad;
  // latitude  = lonlat.v/deg2rad;
}
//-----------------------------------------------------------------------
int EW::computeNearestGridPoint2(int& a_i, int& a_j, int& a_k, int& a_g,
                                 float_sw4 a_x, float_sw4 a_y, float_sw4 a_z) {
  SW4_MARK_FUNCTION;
  int success = 0;
  if (a_z >= m_zmin[mNumberOfCartesianGrids - 1]) {
    // point is in a Cartesian grid
    int g = 0;
    while (g < mNumberOfCartesianGrids && a_z < m_zmin[g]) g++;
    a_g = g;
    a_i = static_cast<int>(floor(a_x / mGridSize[g] + 1));
    a_j = static_cast<int>(floor(a_y / mGridSize[g] + 1));
    a_k = static_cast<int>(round((a_z - m_zmin[g]) / mGridSize[g] + 1));

    VERIFY2(
        a_i >= 1 - m_ghost_points && a_i <= m_global_nx[a_g] + m_ghost_points,
        "Grid Error: i (" << a_i << ") is out of bounds: ( " << 1 << ","
                          << m_global_nx[a_g] << ")"
                          << " x,y,z = " << a_x << " " << a_y << " " << a_z);
    VERIFY2(
        a_j >= 1 - m_ghost_points && a_j <= m_global_ny[a_g] + m_ghost_points,
        "Grid Error: j (" << a_j << ") is out of bounds: ( " << 1 << ","
                          << m_global_ny[a_g] << ")"
                          << " x,y,z = " << a_x << " " << a_y << " " << a_z);
    VERIFY2(a_k >= m_kStart[a_g] && a_k <= m_kEnd[a_g],
            "Grid Error: k (" << a_k << ") is out of bounds: ( " << 1 << ","
                              << m_kEnd[a_g] - m_ghost_points << ")"
                              << " x,y,z = " << a_x << " " << a_y << " "
                              << a_z);
    success = interior_point_in_proc(a_i, a_j, a_g);
  } else {
    // point is in a curvilinear grid
    //
    // foundglobal= Grid point found in at least one processor.
    // success  = Grid point found in my processor (found locally).
    //
    int g = mNumberOfCartesianGrids;
    //  int foundglobal=0;
    success = 0;
    float_sw4 q, r, s;
    while (g < mNumberOfGrids && !success)
    //      while( g < mNumberOfGrids && !foundglobal )
    {
      success = m_gridGenerator->inverse_grid_mapping(this, a_x, a_y, a_z, g, q,
                                                      r, s);
      if (success) {
        a_g = g;
        a_i = static_cast<int>(floor(q));
        a_j = static_cast<int>(floor(r));
        a_k = static_cast<int>(round(s));
      }
      //         MPI_Allreduce(&success,&foundglobal,1,MPI_INT,MPI_MAX,m_cartesian_communicator);
      g++;
    }
    //      VERIFY2( foundglobal, "ERROR in EW:computeNearestGridPoint2, could
    //      not find curvilinear grid point");
  }
  return success;
}

//-------------------------------------------------------
void EW::computeNearestGridPoint(int& a_i, int& a_j, int& a_k,
                                 int& a_g,  // grid on which indices are located
                                 float_sw4 a_x, float_sw4 a_y, float_sw4 a_z) {
  SW4_MARK_FUNCTION;
  bool breakLoop = false;

  for (int g = 0; g < mNumberOfGrids; g++) {
    if (a_z > m_zmin[g] ||
        g == mNumberOfGrids - 1)  // We can not trust zmin for the curvilinear
                                  // grid, since it doesn't mean anything
    {
      a_i = (int)floor(a_x / mGridSize[g]) + 1;
      if (a_x - ((a_i - 0.5) * mGridSize[g]) > 0.) (a_i)++;

      a_j = (int)floor(a_y / mGridSize[g]) + 1;
      if (a_y - ((a_j - 0.5) * mGridSize[g]) > 0.) (a_j)++;

      a_k = (int)floor((a_z - m_zmin[g]) / mGridSize[g]) +
            1;  // Note: this component will be garbage for g=curvilinear grid
      if (a_z - (m_zmin[g] + (a_k - 0.5) * mGridSize[g]) > 0.) (a_k)++;

      a_g = g;

      breakLoop = true;
    } else if (a_z == m_zmin[g])  // testing for equality between doubles is
                                  // kind of pointless...
    {
      // Point is located on top surface if g=finest grid, else the location is
      // on a grid/grid interface, and point is flagged as located on the finer
      // (upper) grid.
      if (g == mNumberOfGrids - 1) {
        a_i = (int)floor(a_x / mGridSize[g]) + 1;
        if (a_x - ((a_i - 0.5) * mGridSize[g]) > 0.) (a_i)++;

        a_j = (int)floor(a_y / mGridSize[g]) + 1;
        if (a_y - ((a_j - 0.5) * mGridSize[g]) > 0.) (a_j)++;

        a_k = 1;

        a_g = g;
      } else {
        a_i = (int)floor(a_x / mGridSize[g + 1]) + 1;
        if (a_x - ((a_i - 0.5) * mGridSize[g + 1]) > 0.) (a_i)++;

        a_j = (int)floor(a_y / mGridSize[g + 1]) + 1;
        if (a_y - ((a_j - 0.5) * mGridSize[g + 1]) > 0.) (a_j)++;

        a_k = (int)floor((a_z - m_zmin[g + 1]) / mGridSize[g + 1]) +
              1;  // Here, I know I am on a grid line

        a_g = g + 1;
      }
      breakLoop = true;
    }

    if (breakLoop) {
      break;
    }
  }

  //  if (m_topography_exists && (a_g == mNumberOfGrids-1)) // The curvilinear
  //  grid will always be the one with the highest number.
  //    {
  // tmp
  //      printf("EW/computeNearestGridPt: You are in the curvilinear part of
  //      the grid, but we do compute the gridpt index using only the Cartesian
  //      grid\n");
  //    }

  // if z > zmax in grid 0 because the coordinate has not yet been corrected for
  // topography, we simply set a_k to m_kEnd
  if (m_topography_exists && a_z >= m_global_zmax) {
    a_k = m_kEnd[0];
    a_g = 0;
  }

  if (!m_topography_exists ||
      (m_topography_exists && a_g < mNumberOfCartesianGrids)) {
    VERIFY2(
        a_i >= 1 - m_ghost_points && a_i <= m_global_nx[a_g] + m_ghost_points,
        "Grid Error: i (" << a_i << ") is out of bounds: ( " << 1 << ","
                          << m_global_nx[a_g] << ")"
                          << " x,y,z = " << a_x << " " << a_y << " " << a_z);
    VERIFY2(
        a_j >= 1 - m_ghost_points && a_j <= m_global_ny[a_g] + m_ghost_points,
        "Grid Error: j (" << a_j << ") is out of bounds: ( " << 1 << ","
                          << m_global_ny[a_g] << ")"
                          << " x,y,z = " << a_x << " " << a_y << " " << a_z);
    VERIFY2(a_k >= m_kStart[a_g] && a_k <= m_kEnd[a_g],
            "Grid Error: k (" << a_k << ") is out of bounds: ( " << 1 << ","
                              << m_kEnd[a_g] - m_ghost_points << ")"
                              << " x,y,z = " << a_x << " " << a_y << " "
                              << a_z);
  }
}

void EW::computeNearestLowGridPoint(
    int& a_i, int& a_j, int& a_k,
    int& a_g,  // grid on which indices are located
    float_sw4 a_x, float_sw4 a_y, float_sw4 a_z) {
  SW4_MARK_FUNCTION;
  bool breakLoop = false;

  for (int g = 0; g < mNumberOfGrids; g++) {
    if (a_z > m_zmin[g] ||
        g == mNumberOfGrids - 1)  // We can not trust zmin for the curvilinear
                                  // grid, since it doesn't mean anything
    {
      a_i = (int)floor(a_x / mGridSize[g]) + 1;
      //          VERIFY(a_x-((a_i-0.5)*mGridSize[g]) <= 0.);

      a_j = (int)floor(a_y / mGridSize[g]) + 1;
      //          VERIFY(a_y-((a_j-0.5)*mGridSize[g]) <= 0.);

      a_k = (int)floor((a_z - m_zmin[g]) / mGridSize[g]) + 1;
      //          VERIFY(a_z-(m_zmin[g]+(a_k-0.5)*mGridSize[g]) <= 0.);

      a_g = g;

      breakLoop = true;
    } else if (a_z == m_zmin[g]) {
      if (g == mNumberOfGrids - 1) {
        a_i = (int)floor(a_x / mGridSize[g]) + 1;
        //              VERIFY(a_x-((a_i-0.5)*mGridSize[g]) <= 0.);

        a_j = (int)floor(a_y / mGridSize[g]) + 1;
        //              VERIFY(a_y-((a_j-0.5)*mGridSize[g]) <= 0.);

        a_k = 1;

        a_g = g;
      } else {
        a_i = (int)floor(a_x / mGridSize[g + 1]) + 1;
        //              VERIFY(a_x-((a_i-0.5)*mGridSize[g+1]) <= 0.);

        a_j = (int)floor(a_y / mGridSize[g + 1]) + 1;
        //              VERIFY(a_y-((a_j-0.5)*mGridSize[g+1]) <= 0.);

        a_k = (int)floor((a_z - m_zmin[g + 1]) / mGridSize[g + 1]) +
              1;  // Here, I know I am on a grid line

        a_g = g + 1;
      }
      breakLoop = true;
    }

    if (breakLoop) {
      break;
    }
  }

  VERIFY2(a_i >= 1 && a_i <= m_global_nx[a_g],
          "Grid Error: i (" << a_i << ") is out of bounds: ( " << 1 << ","
                            << m_global_nx[a_g] << ")");
  VERIFY2(a_j >= 1 && a_j <= m_global_ny[a_g],
          "Grid Error: j (" << a_j << ") is out of bounds: ( " << 1 << ","
                            << m_global_ny[a_g] << ")");
  if (a_k > m_kEnd[a_g] - m_ghost_points) a_k = m_kEnd[a_g] - m_ghost_points;
  if (a_k < 1) a_k = 1;
  //  VERIFY2(a_k >= 1 && a_k <= m_kEnd[a_g]-m_ghost_points,
  //          "Grid Error: k (" << a_k << ") is out of bounds: ( " << 1 << ","
  //          << m_kEnd[a_g]-m_ghost_points << ")");
}

//-----------------------------------------------------------------------
bool EW::interior_point_in_proc(int a_i, int a_j, int a_g) {
  SW4_MARK_FUNCTION;
  // NOT TAKING PARALLEL GHOST POINTS INTO ACCOUNT!
  // Determine if grid point with index (a_i, a_j) on grid a_g is an interior
  // grid point on this processor

  bool retval = false;
  if (a_g >= 0 && a_g < mNumberOfGrids) {
    retval = (a_i >= m_iStartInt[a_g]) && (a_i <= m_iEndInt[a_g]) &&
             (a_j >= m_jStartInt[a_g]) && (a_j <= m_jEndInt[a_g]);
  }
  return retval;
}

//-----------------------------------------------------------------------
bool EW::point_in_proc(int a_i, int a_j, int a_g) {
  // TAKING PARALLEL GHOST POINTS INTO ACCOUNT!
  // Determine if grid point with index (a_i, a_j) on grid a_g is a grid point
  // on this processor

  bool retval = false;
  if (a_g >= 0 && a_g < mNumberOfGrids) {
    retval = (a_i >= m_iStart[a_g] && a_i <= m_iEnd[a_g] &&
              a_j >= m_jStart[a_g] && a_j <= m_jEnd[a_g]);
  }

  return retval;
}

//-----------------------------------------------------------------------
bool EW::point_in_proc_ext(int a_i, int a_j, int a_g) {
  SW4_MARK_FUNCTION;
  // TAKING PARALLEL GHOST POINTS+EXTRA GHOST POINTS INTO ACCOUNT!
  // Determine if grid point with index (a_i, a_j) on grid a_g is a grid point
  // on this processor

  bool retval = false;
  if (a_g >= 0 && a_g < mNumberOfGrids) {
    retval = (a_i >= m_iStart[a_g] - m_ext_ghost_points &&
              a_i <= m_iEnd[a_g] + m_ext_ghost_points &&
              a_j >= m_jStart[a_g] - m_ext_ghost_points &&
              a_j <= m_jEnd[a_g] + m_ext_ghost_points);
  }
  return retval;
}

//-----------------------------------------------------------------------
void EW::getGlobalBoundingBox(float_sw4 bbox[6]) {
  SW4_MARK_FUNCTION;
  bbox[0] = 0.;
  bbox[1] = m_global_xmax;
  bbox[2] = 0.;
  bbox[3] = m_global_ymax;
  bbox[4] = m_global_zmin;
  bbox[5] = m_global_zmax;
}

//-----------------------------------------------------------------------
void EW::setGMTOutput(string filename, string wppfilename) {
  mGMTFileName = filename;
  mWriteGMTOutput = true;

  //  mWPPFileName = wppfilename;
}

//-----------------------------------------------------------------------
void EW::saveGMTFile(vector<vector<Source*>>& a_GlobalUniqueSources,
                     int event) {
  SW4_MARK_FUNCTION;
  // this routine needs to be updated (at least for the etree info)
  if (!mWriteGMTOutput) return;

  if (proc_zero()) {
    stringstream contents;
    contents << "#!/bin/csh\n\n"
             << "gmtset PLOT_DEGREE_FORMAT D\n"
             << "gmtset COLOR_MODEL HSV\n"
             << "gmtset PAPER_MEDIA letter\n"
             << "gmtset PAGE_ORIENTATION portrait\n"
             << "gmtset MEASURE_UNIT inch\n"
             << endl;
    // grab these from grid
    double latNE, lonNE, latSW, lonSW, latSE, lonSE, latNW, lonNW;
    computeGeographicCoord(0.0, 0.0, lonSW, latSW);
    computeGeographicCoord(m_global_xmax, 0.0, lonSE, latSE);
    computeGeographicCoord(m_global_xmax, m_global_ymax, lonNE, latNE);
    computeGeographicCoord(0.0, m_global_ymax, lonNW, latNW);

    // Round up/down
    double minx = std::min(lonSW, std::min(lonSE, std::min(lonNE, lonNW)));
    double maxx = std::max(lonSW, std::max(lonSE, std::max(lonNE, lonNW)));
    double miny = std::min(latSW, std::min(latSE, std::min(latNE, latNW)));
    double maxy = std::max(latSW, std::max(latSE, std::max(latNE, latNW)));
    double margin = 0.1 * fabs(maxy - miny);

    // tmp
    printf("margin = %e\n", margin);

    //      GeographicCoord eNW, eNE, eSW, eSE;

    contents << "# Region will need to be adjusted based on etree/grid values"
             << endl
             << "set REGION = " << minx - margin << "/" << maxx + margin << "/"
             << miny - margin << "/" << maxy + margin << endl
             << endl
             << "set SCALE = 6.0" << endl
             << endl
             << "# These commands are good if you have access to " << endl
             << "# a topography database file for the region modeled " << endl
             << "# Note:  if you uncomment these, adjust the -O -K, etc."
             << endl
             << " #######################################################"
             << endl
             << "#grdraster 2 -R$REGION -I0.5m -Gwpp_topo.grd" << endl
             << "#grdgradient wpp_topo.grd -Gwpp_topo_shade.grd -A270 -Nt -M "
             << endl
             << "#grd2cpt wpp_topo.grd -Ctopo -Z >! wpptopo.cpt" << endl
             << "#grdimage wpp_topo.grd -R$REGION -JM$SCALE -Cwpptopo.cpt "
                "-Iwpp_topo_shade.grd -P -K >! plot.ps"
             << endl
             << " #######################################################"
             << endl
             << "pscoast -R$REGION -JM$SCALE -Bf0.025a0.05 -Dfull "
                "-S100,200,255 -A2000 -W3 -N1t3 -N2t2a -K >! plot.ps"
             << endl
             << endl
             << "# computational grid region..." << endl;

    // Write out gridlines
    contents
        << "psxy -R$REGION -JM$SCALE -W10/255/255/0ta -O -K <<EOF>> plot.ps"
        << endl
        << lonSW << " " << latSW << endl
        << lonSE << " " << latSE << endl
        << lonNE << " " << latNE << endl
        << lonNW << " " << latNW << endl
        << lonSW << " " << latSW << endl
        << "EOF" << endl
        << endl;

    // compute super-grid boundary in (lat-lon) coordinates
    int g = mNumberOfGrids - 1;
    double sg_width = m_sg_gp_thickness * mGridSize[g];
    if (m_use_sg_width) sg_width = m_supergrid_width;

    computeGeographicCoord(sg_width, sg_width, lonSW, latSW);
    computeGeographicCoord(m_global_xmax - sg_width, sg_width, lonNW, latNW);
    computeGeographicCoord(m_global_xmax - sg_width, m_global_ymax - sg_width,
                           lonNE, latNE);
    computeGeographicCoord(sg_width, m_global_ymax - sg_width, lonSE, latSE);

    // Write out gridlines
    contents
        << "#SG boundary: " << endl
        << "psxy -R$REGION -JM$SCALE -W10/255/255/0ta -O -K <<EOF>> plot.ps"
        << endl
        << lonSW << " " << latSW << endl
        << lonSE << " " << latSE << endl
        << lonNE << " " << latNE << endl
        << lonNW << " " << latNW << endl
        << lonSW << " " << latSW << endl
        << "EOF" << endl
        << endl;

    if (a_GlobalUniqueSources[event].size() > 0) {
      contents << "# Sources... " << endl << "cat << EOF >! event.d" << endl;

      for (int i = 0; i < a_GlobalUniqueSources.size(); ++i) {
        double latSource, lonSource;

        computeGeographicCoord(a_GlobalUniqueSources[event][i]->getX0(),
                               a_GlobalUniqueSources[event][i]->getY0(),
                               lonSource, latSource);
        //  should name the event better
        contents << lonSource << " " << latSource << " EVENT-NAME  CB" << endl;
      }
      contents << "EOF" << endl;
      contents << "psxy -R -J -O -K -Sc0.1 -Gred -W0.25p event.d >> plot.ps"
               << endl;
      contents << "awk '{print $1, $2, 12, 1, 9, $4, $3}' event.d | pstext -R "
                  "-J -O -D0.2/0.2v -Gred -N -K >> plot.ps"
               << endl
               << endl;
    }

    int numStations = 0;
    stringstream stationstr;
    stationstr << "# Stations... " << endl;
    stationstr << "cat << EOF >! stations.d " << endl;
    // Write stations by rereading the WPP input file, since some might
    // live outside the grid...
    ifstream sw4InputFile(mName.c_str());
    if (!sw4InputFile.is_open())
      contents << "# Error re-opening input file, skipping stations" << endl;
    else {
      char buffer[256];
      while (!sw4InputFile.eof()) {
        sw4InputFile.getline(buffer, 256);
        if (startswith("rec", buffer) || startswith("sac", buffer)) {
          numStations += 1;
          bool cartCoordSet = false;
          // bool gridPointSet = false;
          bool geoCoordSet = false;
          bool statSet = false;
          string name = "null";
          float_sw4 x = 0.0, y = 0.0;
          float_sw4 lat = 0.0, lon = 0.0;

          // Get location and write to file
          char* token = strtok(buffer, " \t");
          token = strtok(NULL, " \t");  // skip sac
          while (token != NULL) {
            // while there are tokens in the string still
            // NOTE: we skip all verify stmts as these have
            //       already been checked during initial parsing
            if (startswith("#", token) || startswith(" ", buffer))
              // Ignore commented lines and lines with just a space.
              break;
            if (startswith("x=", token)) {
              token += 2;  // skip x=
              cartCoordSet = true;
              x = atof(token);
            } else if (startswith("y=", token)) {
              token += 2;  // skip y=
              cartCoordSet = true;
              y = atof(token);
            } else if (startswith("z=", token)) {
              token += 2;  // skip z=
              cartCoordSet = true;
              // z = atof(token);
            } else if (startswith("lat=", token)) {
              token += 4;  // skip lat=
              lat = atof(token);
              geoCoordSet = true;
            } else if (startswith("lon=", token)) {
              token += 4;  // skip lon=
              lon = atof(token);
              geoCoordSet = true;
            } else if (startswith("depth=", token)) {
              token += 6;  // skip depth=
              // z = atof(token);
              geoCoordSet = true;
            } else if (startswith("sta=", token)) {
              token += 4;
              name = token;
              statSet = true;
            } else if (startswith("file=", token) && !statSet) {
              token += 5;
              name = token;
            }

            token = strtok(NULL, " \t");
          }

          VERIFY(cartCoordSet || geoCoordSet);

          if (!geoCoordSet && cartCoordSet) {
            computeGeographicCoord(x, y, lon, lat);
          }

          // Now have location
          stationstr << lon << " " << lat << " " << name << " CB" << endl;
        }  // token on sac line
      }    // line in ew file
    }

    stationstr << "EOF" << endl << endl;

    stationstr << "# plot station names" << endl
               << "psxy -R -J -O -K -St0.1 -Gblue -W0.25p stations.d >> plot.ps"
               << endl
               << "awk '{print $1, $2, 12, 1, 9, $4, $3}' stations.d | pstext "
                  "-R -J -O -Dj0.3/0.3v -Gblue -N >> plot.ps"
               << endl;

    // Only write station info if there are stations.
    if (numStations > 0) contents << stationstr.str() << endl;

    contents << "/bin/mv plot.ps " << mName << ".ps" << endl;

    stringstream filename;
    filename << mPath[event] << mGMTFileName;
    ofstream gmtfile(filename.str().c_str());
    if (gmtfile.is_open()) {
      cout << "GMT file is open, about to write" << endl;
      gmtfile << contents.str();
      cout << "Wrote GMT file: " << filename.str() << endl;
    } else {
      cout << "Unable to open GMT file: " << filename.str() << endl;
    }

  }  // proc 0
}

//-----------------------------------------------------------------------
void EW::print_execution_time(double t1, double t2, string msg) {
  //   if( !mQuiet && proc_zero() )
  if (proc_zero()) {
    double s = t2 - t1;
    int h = static_cast<int>(s / 3600.0);
    s = s - h * 3600;
    int m = static_cast<int>(s / 60.0);
    s = s - m * 60;
    cout << "   Execution time, " << msg << " ";
    if (h > 1)
      cout << h << " hours ";
    else if (h > 0)
      cout << h << " hour  ";

    if (m > 1)
      cout << m << " minutes ";
    else if (m > 0)
      cout << m << " minute  ";

    if (s > 0) cout << s << " seconds ";
    cout << endl;
  }
}

//-----------------------------------------------------------------------
void EW::print_execution_times(double times[9]) {
  double* time_sums = new double[9 * no_of_procs()];
  MPI_Gather(times, 9, MPI_DOUBLE, time_sums, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  bool printavgs = true;
  if (!mQuiet && proc_zero()) {
    double avgs[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int p = 0; p < no_of_procs(); p++)
      for (int c = 0; c < 9; c++) avgs[c] += time_sums[9 * p + c];
    for (int c = 0; c < 9; c++) avgs[c] /= no_of_procs();
    cout << "\n----------------------------------------" << endl;
    cout << "          Execution time summary (average)" << endl;
    if (printavgs) {
      //                             5                  10          7      2 2
      //                             5       2                      6 7
      cout << "Total      Div-stress Forcing    BC         SG         Comm.    "
              "  MR         Img+TS     Updates "
           << endl;
      cout.setf(ios::left);
      cout.precision(3);
      cout.width(11);
      cout << avgs[0];
      cout.width(11);
      cout << avgs[1];
      cout.width(11);
      cout << avgs[2];
      cout.width(11);
      cout << avgs[3];
      cout.width(11);
      cout << avgs[4];
      cout.width(11);
      cout << avgs[5];
      cout.width(11);
      cout << avgs[6];
      cout.width(11);
      cout << avgs[7];
      cout.width(11);
      cout << avgs[8];
      cout << endl;
    } else {
      cout << "Processor  Total    Div-stress    Forcing     BC     SG     "
              "Comm.    MR    Image+Time-series  Misc  "
           << endl;
      cout.setf(ios::left);
      cout.precision(3);
      for (int p = 0; p < no_of_procs(); p++) {
        cout.width(11);
        cout << p;
        cout << time_sums[9 * p];
        cout.width(11);
        cout << time_sums[9 * p + 1];
        cout.width(11);
        cout << time_sums[9 * p + 2];
        cout.width(11);
        cout << time_sums[9 * p + 3];
        cout.width(11);
        cout << time_sums[9 * p + 4];
        cout.width(11);
        cout << time_sums[9 * p + 5];
        cout.width(11);
        cout << time_sums[9 * p + 6];
        cout.width(11);
        cout << time_sums[9 * p + 7];
        cout.width(11);
        cout << time_sums[9 * p + 8];
        cout << endl;
      }
    }
    cout.setf(ios::right);
    cout.precision(6);
    cout << "----------------------------------------\n" << endl;
  }
  delete[] time_sums;
}

//-----------------------------------------------------------------------
void EW::finalizeIO() {
  //  //  if (proc_zero() && mDebugIO )
  //  //    mEnergyFile.close();
  // for (unsigned int i = 0; i < mSACOutputFiles.size(); ++i)
  //   mSACOutputFiles[i].writeFile();
}

//-----------------------------------------------------------------------
void EW::default_bcs() {
  for (int side = 0; side < 6; side++) mbcGlobalType[side] = bSuperGrid;
  mbcGlobalType[4] = bStressFree;  // low-z is normally free surface
}

//---------------------------------------------------------------------------
void EW::normOfDifference(vector<Sarray>& a_Uex, vector<Sarray>& a_U,
                          float_sw4& diffInf, float_sw4& diffL2,
                          float_sw4& xInf, vector<Source*>& a_globalSources) {
  SW4_MARK_FUNCTION;
  SYNC_STREAM;
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast;
  int imin, imax, jmin, jmax, kmin, kmax;

  float_sw4 *uex_ptr, *u_ptr, h, linfLocal = 0, l2Local = 0, diffInfLocal = 0,
                                 diffL2Local = 0;
  float_sw4 xInfLocal = 0, xInfGrid = 0;
  float_sw4 radius = -1, x0 = 0, y0 = 0, z0 = 0;

  // tmp
  //   if (proc_zero())
  //     printf("Inside normOfDifference\n");
  float_sw4 htop = mGridSize[mNumberOfGrids - 1];
  float_sw4 hbot = mGridSize[0];

  for (g = 0; g < mNumberOfGrids; g++) {
    uex_ptr = a_Uex[g].c_ptr();
    u_ptr = a_U[g].c_ptr();

    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];

    h = mGridSize[g];  // how do we define the grid size for the curvilinear
                       // grid?

    // don't think this is correct:
    //    int nsgxy = (int)(0.5+m_sg_gp_thickness*htop/h);
    //    int nsgz  = (int)(0.5+m_sg_gp_thickness*hbot/h);
    int nsgxy = static_cast<int>(0.5 + m_sg_gp_thickness);
    int nsgz = static_cast<int>(0.5 + m_sg_gp_thickness);
    if (m_use_sg_width) {
      nsgxy = static_cast<int>(ceil(m_supergrid_width / h));
      nsgz = static_cast<int>(ceil(m_supergrid_width / h));
    }

    if (mbcGlobalType[0] == bSuperGrid)
      imin = max(m_iStartInt[g], nsgxy + 1);
    else
      imin = m_iStartInt[g];

    if (mbcGlobalType[1] == bSuperGrid)
      imax = min(m_iEndInt[g], m_global_nx[g] - nsgxy);
    else
      imax = m_iEndInt[g];

    if (mbcGlobalType[2] == bSuperGrid)
      jmin = max(m_jStartInt[g], nsgxy + 1);
    else
      jmin = m_jStartInt[g];

    if (mbcGlobalType[3] == bSuperGrid)
      jmax = min(m_jEndInt[g], m_global_ny[g] - nsgxy);
    else
      jmax = m_jEndInt[g];

    // Can not test on global type when there is more than one grid in the
    // z-direction if uppermost grid has layer on top boundary, the fine grid
    // spacing is used for the s.g. layer width
    if (m_bcType[g][4] == bSuperGrid)
      kmin = max(m_kStartInt[g], nsgxy + 1);
    else
      kmin = m_kStartInt[g];
    // The lowermost grid has the s.g. layer width based on the spacing of the
    // coarsest grid
    if (m_bcType[g][5] == bSuperGrid)
      kmax = min(m_kEndInt[g], m_global_nz[g] - nsgz);
    else
      kmax = m_kEndInt[g];

    // tmp
    //     printf("proc=%i, iS= %i, iE=%i, jS=%i, jE=%i, kS=%i, kE=%i\n",
    //     m_myRank,
    // 	   m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
    // m_kEnd[g]);
    //     printf("proc=%i, if= %i, il=%i, jf=%i, jl=%i, kf=%i, kl=%i\n",
    //     m_myRank,
    // 	   ifirst, ilast, jfirst, jlast, kfirst, klast);

    if (m_point_source_test) {
      radius = 4 * h;
      x0 = a_globalSources[0]->getX0();
      y0 = a_globalSources[0]->getY0();
      z0 = a_globalSources[0]->getZ0();
    }
    // need to exclude parallel overlap from L2 calculation
    int usesg = usingSupergrid();
    // if (topographyExists() && g == mNumberOfGrids - 1) {
    if (topographyExists() && g >= mNumberOfCartesianGrids) {
      if (m_croutines)
        solerr3c_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, uex_ptr, u_ptr,
                    mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr(), mJ[g].c_ptr(),
                    linfLocal, l2Local, xInfGrid, x0, y0, z0, radius, imin,
                    imax, jmin, jmax, kmin, kmax, usesg, m_sg_str_x[g],
                    m_sg_str_y[g]);
      else
        solerr3c(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, uex_ptr,
                 u_ptr, mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr(),
                 mJ[g].c_ptr(), &linfLocal, &l2Local, &xInfGrid, &x0, &y0, &z0,
                 &radius, &imin, &imax, &jmin, &jmax, &kmin, &kmax, &usesg,
                 m_sg_str_x[g], m_sg_str_y[g]);
    } else {
      if (m_croutines)
        solerr3_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, h, uex_ptr,
                   u_ptr, linfLocal, l2Local, xInfGrid, m_zmin[g], x0, y0, z0,
                   radius, imin, imax, jmin, jmax, kmin, kmax);
      else
        solerr3(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &h, uex_ptr,
                u_ptr, &linfLocal, &l2Local, &xInfGrid, &m_zmin[g], &x0, &y0,
                &z0, &radius, &imin, &imax, &jmin, &jmax, &kmin, &kmax);
    }
    if (linfLocal > diffInfLocal) diffInfLocal = linfLocal;
    if (xInfGrid > xInfLocal) xInfLocal = xInfGrid;
    // std::cout<<"NORM"<<g<<" "<<xInfGrid<<" "<<xInfLocal<<"\n";
    diffL2Local += l2Local;
  }
  // communicate local results for global errors
  MPI_Allreduce(&diffInfLocal, &diffInf, 1, m_mpifloat, MPI_MAX,
                m_cartesian_communicator);
  MPI_Allreduce(&xInfLocal, &xInf, 1, m_mpifloat, MPI_MAX,
                m_cartesian_communicator);
  MPI_Allreduce(&diffL2Local, &diffL2, 1, m_mpifloat, MPI_SUM,
                m_cartesian_communicator);

  diffL2 = sqrt(diffL2);
}

//---------------------------------------------------------------------------
void EW::normOfDifferenceGhostPoints(vector<Sarray>& a_Uex, vector<Sarray>& a_U,
                                     float_sw4& diffInf, float_sw4& diffL2) {
  SW4_MARK_FUNCTION;
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *uex_ptr, *u_ptr, h, linfLocal = 0, l2Local = 0, diffInfLocal = 0,
                                 diffL2Local = 0;

  // tmp
  //  if (proc_zero())
  //    printf("Inside normOfDifferenceGhostPoints\n");

  for (g = 0; g < mNumberOfGrids; g++) {
    uex_ptr = a_Uex[g].c_ptr();
    u_ptr = a_U[g].c_ptr();

    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];

    h = mGridSize[g];  // how do we define the grid size for the curvilinear
                       // grid?

    // need to exclude parallel overlap from L2 calculation
    if (m_croutines)
      solerrgp_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, h, uex_ptr,
                  u_ptr, linfLocal, l2Local);
    else
      solerrgp(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &h, uex_ptr,
               u_ptr, &linfLocal, &l2Local);
    if (linfLocal > diffInfLocal) diffInfLocal = linfLocal;
    diffL2Local += l2Local;
    //    cout << m_myRank << " g, l2, li = " << " " << g << " " << l2Local << "
    //    " << linfLocal << endl;
  }

  //  cout << m_myRank << " l2, li = " << diffL2Local << " " << diffInfLocal <<
  //  endl;

  // communicate local results for global errors
  MPI_Allreduce(&diffInfLocal, &diffInf, 1, m_mpifloat, MPI_MAX,
                m_cartesian_communicator);
  MPI_Allreduce(&diffL2Local, &diffL2, 1, m_mpifloat, MPI_SUM,
                m_cartesian_communicator);

  //   diffL2 = diffL2Local;
  //   diffInf = diffInfLocal;

  diffL2 = sqrt(diffL2);
}

//---------------------------------------------------------------------------
void EW::normOfSurfaceDifference(vector<Sarray>& a_Uex, vector<Sarray>& a_U,
                                 float_sw4& diffInf, float_sw4& diffL2,
                                 float_sw4& solInf, float_sw4& solL2,
                                 vector<Source*>& a_globalSources) {
  SW4_MARK_FUNCTION;
  int g;
  float_sw4 absDiff, absSol;
  float_sw4 h, diffInfLocal = 0, diffL2Local = 0, solInfLocal = 0,
               solL2Local = 0;

  g = mNumberOfCartesianGrids - 1;
  int k = 1;

  h = mGridSize[g];

  // only evaluate error on the surface, not including ghost or parallel overlap
  // points

  // also exclude points in the super grid damping layer
  int imin, imax, jmin, jmax;
  int nsgxy = m_sg_gp_thickness;
  if (m_use_sg_width) nsgxy = static_cast<int>(floor(m_supergrid_width / h));

  if (mbcGlobalType[0] == bSuperGrid)
    imin = max(m_iStartInt[g], nsgxy + 1);
  else
    imin = m_iStartInt[g];

  if (mbcGlobalType[1] == bSuperGrid)
    imax = min(m_iEndInt[g], m_global_nx[g] - nsgxy);
  else
    imax = m_iEndInt[g];

  if (mbcGlobalType[2] == bSuperGrid)
    jmin = max(m_jStartInt[g], nsgxy + 1);
  else
    jmin = m_jStartInt[g];

  if (mbcGlobalType[3] == bSuperGrid)
    jmax = min(m_jEndInt[g], m_global_ny[g] - nsgxy);
  else
    jmax = m_jEndInt[g];

  // also need to exclude grid points near the point source
  h = mGridSize[g];

  float_sw4 radius2, x0, y0, dist2;

  if (m_lamb_test) {
    radius2 = SQR(4 * h);
    x0 = a_globalSources[0]->getX0();
    y0 = a_globalSources[0]->getY0();
  } else {
    radius2 = -1;
    x0 = 0;
    y0 = 0;
  }

  for (int j = jmin; j <= jmax; j++)
    for (int i = imin; i <= imax; i++) {
      dist2 = SQR((i - 1) * h - x0) + SQR((j - 1) * h - y0);

      if (dist2 > radius2) {
        absDiff = fabs(a_Uex[g](3, i, j, k) - a_U[g](3, i, j, k));
        if (absDiff > diffInfLocal) diffInfLocal = absDiff;
        diffL2Local += h * h * absDiff * absDiff;
        // exact sol norm
        absSol = fabs(a_Uex[g](3, i, j, k));
        if (absSol > solInfLocal) solInfLocal = absSol;
        solL2Local += h * h * absSol * absSol;
      }
    }

  // communicate local results for global errors
  MPI_Allreduce(&diffInfLocal, &diffInf, 1, m_mpifloat, MPI_MAX,
                m_cartesian_communicator);
  MPI_Allreduce(&diffL2Local, &diffL2, 1, m_mpifloat, MPI_SUM,
                m_cartesian_communicator);

  MPI_Allreduce(&solInfLocal, &solInf, 1, m_mpifloat, MPI_MAX,
                m_cartesian_communicator);
  MPI_Allreduce(&solL2Local, &solL2, 1, m_mpifloat, MPI_SUM,
                m_cartesian_communicator);

  diffL2 = sqrt(diffL2);
  solL2 = sqrt(solL2);
}

//---------------------------------------------------------------------------
void EW::bndryInteriorDifference(vector<Sarray>& a_Uex, vector<Sarray>& a_U,
                                 float_sw4* lowZ, float_sw4* interiorZ,
                                 float_sw4* highZ) {
  SW4_MARK_FUNCTION;
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nz;
  float_sw4 *uex_ptr, *u_ptr, h;

  for (g = 0; g < mNumberOfGrids; g++) {
    uex_ptr = a_Uex[g].c_ptr();
    u_ptr = a_U[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    h = mGridSize[g];
    nz = m_global_nz[g];

    // need to do a gather over all processors
    if (m_croutines)
      rhserrfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz, h, uex_ptr,
                    u_ptr, &lowZ[3 * g], &interiorZ[3 * g], &highZ[3 * g]);
    else
      rhserrfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz, &h,
                 uex_ptr, u_ptr, &lowZ[3 * g], &interiorZ[3 * g],
                 &highZ[3 * g]);
  }
}

//---------------------------------------------------------------------------
void EW::test_RhoUtt_Lu(vector<Sarray>& a_Uacc, vector<Sarray>& a_Lu,
                        vector<Sarray>& a_F, float_sw4* lowZ,
                        float_sw4* interiorZ, float_sw4* highZ) {
  SW4_MARK_FUNCTION;
  int g, ifirst, ilast, jfirst, jlast, kfirst, klast, nz;
  float_sw4 *rho_ptr, *uacc_ptr, *lu_ptr, *f_ptr;

  for (g = 0; g < mNumberOfGrids; g++) {
    rho_ptr = mRho[g].c_ptr();
    uacc_ptr = a_Uacc[g].c_ptr();
    lu_ptr = a_Lu[g].c_ptr();
    f_ptr = a_F[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    // h = mGridSize[g]; // how do we define the grid size for the curvilinear
    // grid?
    nz = m_global_nz[g];

    // evaluate rho*uacc - lu - f in fortran routine
    if (m_croutines)
      rhouttlumf_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz, uacc_ptr,
                    lu_ptr, f_ptr, rho_ptr, &lowZ[3 * g], &interiorZ[3 * g],
                    &highZ[3 * g]);
    else
      rhouttlumf(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
                 uacc_ptr, lu_ptr, f_ptr, rho_ptr, &lowZ[3 * g],
                 &interiorZ[3 * g], &highZ[3 * g]);
  }
}

//---------------------------------------------------------------------------
void EW::initialData(float_sw4 a_t, vector<Sarray>& a_U,
                     vector<Sarray*>& a_AlphaVE) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *u_ptr, om, ph, cv, h, zmin;

  if (m_twilight_forcing) {
    for (int g = 0; g < mNumberOfCartesianGrids; g++) {
      u_ptr = a_U[g].c_ptr();
      ifirst = m_iStart[g];
      ilast = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast = m_jEnd[g];
      kfirst = m_kStart[g];
      klast = m_kEnd[g];
      h = mGridSize[g];  // how do we define the grid size for the curvilinear
                         // grid?
      zmin = m_zmin[g];
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
      if (m_croutines)
        twilightfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, u_ptr, a_t,
                        om, cv, ph, h, zmin);
      else
        twilightfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, u_ptr,
                     &a_t, &om, &cv, &ph, &h, &zmin);
      if (m_use_attenuation) {
        // one mechanism is assumed
        float_sw4* alpha_ptr = a_AlphaVE[g][0].c_ptr();
        if (m_croutines)
          twilightfortatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                             alpha_ptr, a_t, om, cv, ph, h, zmin);
        else
          twilightfortatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                          alpha_ptr, &a_t, &om, &cv, &ph, &h, &zmin);
      }
    }
    //    if (topographyExists()) {
    for (int g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
      // int g = mNumberOfGrids - 1;
      u_ptr = a_U[g].c_ptr();
      ifirst = m_iStart[g];
      ilast = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast = m_jEnd[g];
      kfirst = m_kStart[g];
      klast = m_kEnd[g];
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
      if (m_croutines)
        twilightfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, u_ptr,
                         a_t, om, cv, ph, mX[g].c_ptr(), mY[g].c_ptr(),
                         mZ[g].c_ptr());
      else
        twilightfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, u_ptr,
                      &a_t, &om, &cv, &ph, mX[g].c_ptr(), mY[g].c_ptr(),
                      mZ[g].c_ptr());
      if (m_use_attenuation) {
        // one mechanism is assumed
        float_sw4* alpha_ptr = a_AlphaVE[g][0].c_ptr();
        if (m_croutines)
          twilightfortattc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                              alpha_ptr, a_t, om, cv, ph, mX[g].c_ptr(),
                              mY[g].c_ptr(), mZ[g].c_ptr());
        else
          twilightfortattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                           alpha_ptr, &a_t, &om, &cv, &ph, mX[g].c_ptr(),
                           mY[g].c_ptr(), mZ[g].c_ptr());
      }
    }
  } else if (m_rayleigh_wave_test) {
    double cr, lambda, mu, rho, alpha;
    for (int g = 0; g < mNumberOfCartesianGrids;
         g++)  // This case does not make sense with topography
    {
      ifirst = m_iStart[g];
      ilast = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast = m_jEnd[g];
      kfirst = m_kStart[g];
      klast = m_kEnd[g];
      h = mGridSize[g];  // how do we define the grid size for the curvilinear
                         // grid?
      zmin = m_zmin[g];
      om = m_rayleigh_wave_test->m_omega;
      cr = m_rayleigh_wave_test->m_cr;
      rho = m_rayleigh_wave_test->m_rho;
      lambda = m_rayleigh_wave_test->m_lambda;
      mu = m_rayleigh_wave_test->m_mu;
      alpha = m_rayleigh_wave_test->m_alpha;
      size_t npts = a_U[g].m_npts;
      double* uini = new double[npts];
      rayleighfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, uini,
                   &a_t, &lambda, &mu, &rho, &cr, &om, &alpha, &h, &zmin);
      a_U[g].assign(uini, 0);
      delete[] uini;
    }
  } else if (m_energy_test) {
    for (int g = 0; g < mNumberOfGrids; g++)  // ranomized initial data
    {
      u_ptr = a_U[g].c_ptr();
      size_t npts = (static_cast<size_t>(m_iEnd[g] - m_iStart[g] + 1)) *
                    (m_jEnd[g] - m_jStart[g] + 1) *
                    (m_kEnd[g] - m_kStart[g] + 1);
      if (m_croutines) {
        // Loop to make c-order and fortran-order have same random number
        // sequence
        for (size_t i = 0; i < npts; i++) {
          u_ptr[i] = drand48();
          u_ptr[i + npts] = drand48();
          u_ptr[i + 2 * npts] = drand48();
        }
      } else {
        for (size_t i = 0; i < npts; i++) {
          u_ptr[3 * i] = drand48();
          u_ptr[3 * i + 1] = drand48();
          u_ptr[3 * i + 2] = drand48();
        }
      }
    }  // end for g

  }  // end m_energy_test
  else
  // homogeneous initial data is the default
  {
    for (int g = 0; g < mNumberOfGrids; g++) {
      a_U[g].set_to_zero();
      for (int a = 0; a < m_number_mechanisms; a++)
        a_AlphaVE[g][a].set_to_zero();
    }
  }
}

//---------------------------------------------------------------------------
bool EW::exactSol(float_sw4 a_t, vector<Sarray>& a_U,
                  vector<Sarray*>& a_AlphaVE, vector<Source*>& sources) {
  SW4_MARK_FUNCTION;
  SYNC_STREAM;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *u_ptr, om, ph, cv, h, zmin;
  bool retval;

  if (m_twilight_forcing) {
    for (int g = 0; g < mNumberOfCartesianGrids;
         g++)  // curvilinear case needs to be implemented
    {
      u_ptr = a_U[g].c_ptr();
      ifirst = m_iStart[g];
      ilast = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast = m_jEnd[g];
      kfirst = m_kStart[g];
      klast = m_kEnd[g];
      h = mGridSize[g];  // how do we define the grid size for the curvilinear
                         // grid?
      zmin = m_zmin[g];
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
      if (m_croutines)
        twilightfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, u_ptr, a_t,
                        om, cv, ph, h, zmin);
      else
        twilightfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, u_ptr,
                     &a_t, &om, &cv, &ph, &h, &zmin);
      if (m_use_attenuation) {
        // one mechanism is assumed
        float_sw4* alpha_ptr = a_AlphaVE[g][0].c_ptr();
        if (m_croutines)
          twilightfortatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                             alpha_ptr, a_t, om, cv, ph, h, zmin);
        else
          twilightfortatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                          alpha_ptr, &a_t, &om, &cv, &ph, &h, &zmin);
      }
    }
    //   if (topographyExists()) {
    for (int g = mNumberOfCartesianGrids; g < mNumberOfGrids;
         g++)  // curvilinear grids
    {
      // int g = mNumberOfGrids - 1;
      u_ptr = a_U[g].c_ptr();
      ifirst = m_iStart[g];
      ilast = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast = m_jEnd[g];
      kfirst = m_kStart[g];
      klast = m_kEnd[g];
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
      if (m_croutines)
        twilightfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, u_ptr,
                         a_t, om, cv, ph, mX[g].c_ptr(), mY[g].c_ptr(),
                         mZ[g].c_ptr());
      else
        twilightfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, u_ptr,
                      &a_t, &om, &cv, &ph, mX[g].c_ptr(), mY[g].c_ptr(),
                      mZ[g].c_ptr());
      if (m_use_attenuation) {
        // std::cout<<"THI IS THE ONE\n";
        // one mechanism is assumed
        float_sw4* alpha_ptr = a_AlphaVE[g][0].c_ptr();
        if (m_croutines)
          twilightfortattc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                              alpha_ptr, a_t, om, cv, ph, mX[g].c_ptr(),
                              mY[g].c_ptr(), mZ[g].c_ptr());
        else
          twilightfortattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                           alpha_ptr, &a_t, &om, &cv, &ph, mX[g].c_ptr(),
                           mY[g].c_ptr(), mZ[g].c_ptr());
      }
    }
    retval = true;
  } else if (m_point_source_test) {
    for (int g = 0; g < mNumberOfGrids; g++) {
      size_t npts = a_U[g].m_npts;
      float_sw4* uexact = SW4_NEW(Space::Managed, float_sw4[npts]);
#if defined(ENABLE_CUDA)
      SW4_CheckDeviceError(cudaMemPrefetchAsync(
          uexact, npts * sizeof(float_sw4), global_variables.device, 0));
#endif
      //       get_exact_point_source( a_U[g].c_ptr(), a_t, g, *sources[0] );
      get_exact_point_source(uexact, a_t, g, *sources[0]);
      a_U[g].prefetch();
      a_U[g].assign(uexact, 0);
      ::operator delete[](uexact, Space::Managed);
    }
    retval = true;
  } else if (m_lamb_test) {
    get_exact_lamb2(a_U, a_t, *sources[0]);
    retval = true;
  } else if (m_rayleigh_wave_test) {
    double cr, lambda, mu, rho, alpha;
    for (int g = 0; g < mNumberOfCartesianGrids;
         g++)  // This case does not make sense with topography
    {
      ifirst = m_iStart[g];
      ilast = m_iEnd[g];
      jfirst = m_jStart[g];
      jlast = m_jEnd[g];
      kfirst = m_kStart[g];
      klast = m_kEnd[g];
      h = mGridSize[g];  // how do we define the grid size for the curvilinear
                         // grid?
      zmin = m_zmin[g];
      om = m_rayleigh_wave_test->m_omega;
      cr = m_rayleigh_wave_test->m_cr;
      rho = m_rayleigh_wave_test->m_rho;
      lambda = m_rayleigh_wave_test->m_lambda;
      mu = m_rayleigh_wave_test->m_mu;
      alpha = m_rayleigh_wave_test->m_alpha;
      size_t npts = a_U[g].m_npts;
      double* uexact = new double[npts];
      rayleighfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, uexact,
                   &a_t, &lambda, &mu, &rho, &cr, &om, &alpha, &h, &zmin);
      a_U[g].assign(uexact, 0);
      delete[] uexact;
    }

    retval = true;
  } else  // In general, the exact solution is unknown (m_energy_test falls into
          // this category)

  {
    retval = false;
  }
  return retval;
}

//-----------------------------------------------------------------------
// smooth wave for time dependence to test point force term with
RAJA_HOST_DEVICE float_sw4 EW::SmoothWave(float_sw4 t, float_sw4 R,
                                          float_sw4 c) {
  float_sw4 temp = R;
  float_sw4 c0 = 2187. / 8., c1 = -10935. / 8., c2 = 19683. / 8.,
            c3 = -15309. / 8., c4 = 2187. / 4.;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1,
  //  (c0*pow(t-R/c,3)+c1*pow(t-R/c,4)+c2*pow(t-R/c,5)+c3*pow(t-R/c,6)+c4*pow(t-R/c,7)),
  //  0);
  if ((t - R / c) > 0 && (t - R / c) < 1)
    temp = (c0 * pow(t - R / c, 3) + c1 * pow(t - R / c, 4) +
            c2 * pow(t - R / c, 5) + c3 * pow(t - R / c, 6) +
            c4 * pow(t - R / c, 7));
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// very smooth bump for time dependence for further testing of point force
RAJA_HOST_DEVICE float_sw4 EW::VerySmoothBump(float_sw4 t, float_sw4 R,
                                              float_sw4 c) {
  float_sw4 temp = R;
  float_sw4 c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120,
            c5 = -1024;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1,
  //  (c0*pow(t-R/c,5)+c1*pow(t-R/c,6)+c2*pow(t-R/c,7)+c3*pow(t-R/c,8)+c4*pow(t-R/c,9)+c5*pow(t-R/c,10)),
  //  0);
  if ((t - R / c) > 0 && (t - R / c) < 1)
    temp = (c0 * pow(t - R / c, 5) + c1 * pow(t - R / c, 6) +
            c2 * pow(t - R / c, 7) + c3 * pow(t - R / c, 8) +
            c4 * pow(t - R / c, 9) + c5 * pow(t - R / c, 10));
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// C6 smooth bump for time dependence for further testing of point force
RAJA_HOST_DEVICE SW4_FORCEINLINE float_sw4 EW::C6SmoothBump(float_sw4 t,
                                                            float_sw4 R,
                                                            float_sw4 c) {
  float_sw4 retval = 0;
  if ((t - R / c) > 0 && (t - R / c) < 1)
    retval = 51480.0 * pow((t - R / c) * (1 - t + R / c), 7);
  return retval;
}

//-----------------------------------------------------------------------
// derivative of smooth wave
RAJA_HOST_DEVICE float_sw4 EW::d_SmoothWave_dt(float_sw4 t, float_sw4 R,
                                               float_sw4 c) {
  float_sw4 temp = R;
  float_sw4 c0 = 2187. / 8., c1 = -10935. / 8., c2 = 19683. / 8.,
            c3 = -15309. / 8., c4 = 2187. / 4.;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1,
  //  (3*c0*pow(t-R/c,2)+4*c1*pow(t-R/c,3)+5*c2*pow(t-R/c,4)+6*c3*pow(t-R/c,5)+7*c4*pow(t-R/c,6)),
  //  0);
  if ((t - R / c) > 0 && (t - R / c) < 1)
    temp = (3 * c0 * pow(t - R / c, 2) + 4 * c1 * pow(t - R / c, 3) +
            5 * c2 * pow(t - R / c, 4) + 6 * c3 * pow(t - R / c, 5) +
            7 * c4 * pow(t - R / c, 6));
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// very smooth bump for time dependence to further testing of point force
RAJA_HOST_DEVICE float_sw4 EW::d_VerySmoothBump_dt(float_sw4 t, float_sw4 R,
                                                   float_sw4 c) {
  float_sw4 temp = R;
  float_sw4 c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120,
            c5 = -1024;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1,
  //  (5*c0*pow(t-R/c,4)+6*c1*pow(t-R/c,5)+7*c2*pow(t-R/c,6)+8*c3*pow(t-R/c,7)+9*c4*pow(t-R/c,8))+10*c5*pow(t-R/c,9),
  //  0);
  if ((t - R / c) > 0 && (t - R / c) < 1)
    temp = (5 * c0 * pow(t - R / c, 4) + 6 * c1 * pow(t - R / c, 5) +
            7 * c2 * pow(t - R / c, 6) + 8 * c3 * pow(t - R / c, 7) +
            9 * c4 * pow(t - R / c, 8)) +
           10 * c5 * pow(t - R / c, 9);
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// C6 smooth bump for time dependence to further testing of point force
RAJA_HOST_DEVICE SW4_FORCEINLINE float_sw4 EW::d_C6SmoothBump_dt(float_sw4 t,
                                                                 float_sw4 R,
                                                                 float_sw4 c) {
  float_sw4 retval = 0;
  if ((t - R / c) > 0 && (t - R / c) < 1)
    retval = 51480.0 * 7 * (1 - 2 * (t - R / c)) *
             pow((t - R / c) * (1 - t + R / c), 6);
  return retval;
}

//-----------------------------------------------------------------------
// Primitive function (for T) of SmoothWave(t-T)*T
RAJA_HOST_DEVICE float_sw4 EW::SWTP(float_sw4 Lim, float_sw4 t) {
  float_sw4 temp = Lim;

  float_sw4 c0 = 2187. / 8., c1 = -10935. / 8., c2 = 19683. / 8.,
            c3 = -15309. / 8., c4 = 2187. / 4.;

  temp = (pow(t, 3) *
          (c0 + c1 * t + c2 * pow(t, 2) + c3 * pow(t, 3) + c4 * pow(t, 4)) *
          pow(Lim, 2)) /
             2. -
         (pow(t, 2) *
          (3 * c0 + 4 * c1 * t + 5 * c2 * pow(t, 2) + 6 * c3 * pow(t, 3) +
           7 * c4 * pow(t, 4)) *
          pow(Lim, 3)) /
             3. +
         (t *
          (3 * c0 + 6 * c1 * t + 10 * c2 * pow(t, 2) + 15 * c3 * pow(t, 3) +
           21 * c4 * pow(t, 4)) *
          pow(Lim, 4)) /
             4. +
         ((-c0 - 4 * c1 * t - 10 * c2 * pow(t, 2) - 20 * c3 * pow(t, 3) -
           35 * c4 * pow(t, 4)) *
          pow(Lim, 5)) /
             5. +
         ((c1 + 5 * c2 * t + 15 * c3 * pow(t, 2) + 35 * c4 * pow(t, 3)) *
          pow(Lim, 6)) /
             6. +
         ((-c2 - 6 * c3 * t - 21 * c4 * pow(t, 2)) * pow(Lim, 7)) / 7. +
         ((c3 + 7 * c4 * t) * pow(Lim, 8)) / 8. - (c4 * pow(Lim, 9)) / 9.;

  return temp;
}

//-----------------------------------------------------------------------
// Primitive function (for T) of VerySmoothBump(t-T)*T
RAJA_HOST_DEVICE float_sw4 EW::VSBTP(float_sw4 Lim, float_sw4 t) {
  float_sw4 temp = Lim;
  float_sw4 f = 1024., g = -5120., h = 10240., i = -10240., j = 5120.,
            k = -1024.;

  temp =
      (pow(Lim, 11) * (-25200 * k * t - 2520 * j) + 2310 * k * pow(Lim, 12) +
       (124740 * k * pow(t, 2) + 24948 * j * t + 2772 * i) * pow(Lim, 10) +
       (-369600 * k * pow(t, 3) - 110880 * j * pow(t, 2) - 24640 * i * t -
        3080 * h) *
           pow(Lim, 9) +
       (727650 * k * pow(t, 4) + 291060 * j * pow(t, 3) +
        97020 * i * pow(t, 2) + 24255 * h * t + 3465 * g) *
           pow(Lim, 8) +
       (-997920 * k * pow(t, 5) - 498960 * j * pow(t, 4) -
        221760 * i * pow(t, 3) - 83160 * h * pow(t, 2) - 23760 * g * t -
        3960 * f) *
           pow(Lim, 7) +
       (970200 * k * pow(t, 6) + 582120 * j * pow(t, 5) +
        323400 * i * pow(t, 4) + 161700 * h * pow(t, 3) +
        69300 * g * pow(t, 2) + 23100 * f * t) *
           pow(Lim, 6) +
       (-665280 * k * pow(t, 7) - 465696 * j * pow(t, 6) -
        310464 * i * pow(t, 5) - 194040 * h * pow(t, 4) -
        110880 * g * pow(t, 3) - 55440 * f * pow(t, 2)) *
           pow(Lim, 5) +
       (311850 * k * pow(t, 8) + 249480 * j * pow(t, 7) +
        194040 * i * pow(t, 6) + 145530 * h * pow(t, 5) +
        103950 * g * pow(t, 4) + 69300 * f * pow(t, 3)) *
           pow(Lim, 4) +
       (-92400 * k * pow(t, 9) - 83160 * j * pow(t, 8) - 73920 * i * pow(t, 7) -
        64680 * h * pow(t, 6) - 55440 * g * pow(t, 5) - 46200 * f * pow(t, 4)) *
           pow(Lim, 3) +
       (13860 * k * pow(t, 10) + 13860 * j * pow(t, 9) + 13860 * i * pow(t, 8) +
        13860 * h * pow(t, 7) + 13860 * g * pow(t, 6) + 13860 * f * pow(t, 5)) *
           pow(Lim, 2)) /
      27720.0;

  return temp;
}
//-----------------------------------------------------------------------
// Primitive function (for T) of C6SmoothBump(t-T)*T
RAJA_HOST_DEVICE SW4_FORCEINLINE float_sw4 EW::C6SBTP(float_sw4 Lim,
                                                      float_sw4 t) {
  float_sw4 x = t - Lim;
  return pow(x, 8) *
         (-3217.5 * pow(x, 8) + 3432.0 * (7 + t) * pow(x, 7) -
          25740.0 * (3 + t) * pow(x, 6) + 27720.0 * (5 + 3 * t) * pow(x, 5) -
          150150.0 * (t + 1) * x * x * x * x +
          32760.0 * (3 + 5 * t) * x * x * x - 36036.0 * (1 + 3 * t) * x * x +
          5720.0 * (1 + 7 * t) * x - 6435.0 * t);
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*SmoothWave(t-T)*T from R/alpha to R/beta
RAJA_HOST_DEVICE float_sw4 EW::SmoothWave_x_T_Integral(float_sw4 t, float_sw4 R,
                                                       float_sw4 alpha,
                                                       float_sw4 beta) {
  float_sw4 temp = R;

  float_sw4 lowL, hiL;

  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t,
  //  R / beta, t);
  if ((R / alpha > t - 1))
    lowL = R / alpha;
  else
    lowL = t - 1;
  if (R / beta < t)
    hiL = R / beta;
  else
    hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, SWTP(hiL, t) - SWTP(lowL, t), 0.0);
  if (lowL < t && hiL > t - 1)
    temp = SWTP(hiL, t) - SWTP(lowL, t);
  else
    temp = 0;

  return temp;
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*VerySmoothBump(t-T)*T from R/alpha to R/beta
RAJA_HOST_DEVICE float_sw4 EW::VerySmoothBump_x_T_Integral(float_sw4 t,
                                                           float_sw4 R,
                                                           float_sw4 alpha,
                                                           float_sw4 beta) {
  float_sw4 temp = R;

  float_sw4 lowL, hiL;

  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t,
  //  R / beta, t);
  if (R / alpha > t - 1)
    lowL = R / alpha;
  else
    lowL = t - 1;
  if (R / beta < t)
    hiL = R / beta;
  else
    hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t),
  //  0.0);
  if (lowL < t && hiL > t - 1)
    temp = VSBTP(hiL, t) - VSBTP(lowL, t);
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*C6SmoothBump(t-T)*T from R/alpha to R/beta
RAJA_HOST_DEVICE SW4_FORCEINLINE float_sw4 EW::C6SmoothBump_x_T_Integral(
    float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta) {
  float_sw4 temp = R;

  float_sw4 lowL, hiL;

  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t,
  //  R / beta, t);
  if (R / alpha > t - 1)
    lowL = R / alpha;
  else
    lowL = t - 1;
  if (R / beta < t)
    hiL = R / beta;
  else
    hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t),
  //  0.0);
  if (lowL < t && hiL > t - 1)
    temp = C6SBTP(hiL, t) - C6SBTP(lowL, t);
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE float_sw4 EW::Gaussian(float_sw4 t, float_sw4 R, float_sw4 c,
                                        float_sw4 f) {
  float_sw4 temp = R;
  temp = 1 / (f * sqrt(2 * M_PI)) * exp(-pow(t - R / c, 2) / (2 * f * f));
  return temp;
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE float_sw4 EW::d_Gaussian_dt(float_sw4 t, float_sw4 R,
                                             float_sw4 c, float_sw4 f) {
  float_sw4 temp = R;
  temp = 1 / (f * sqrt(2 * M_PI)) *
         (-exp(-pow(t - R / c, 2) / (2 * f * f)) * (t - R / c)) / pow(f, 2);
  return temp;
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE float_sw4 EW::Gaussian_x_T_Integral(float_sw4 t, float_sw4 R,
                                                     float_sw4 f,
                                                     float_sw4 alpha,
                                                     float_sw4 beta) {
  float_sw4 temp = R;
  temp = -0.5 * t *
             (erf((t - R / beta) / (sqrt(2.0) * f)) -
              erf((t - R / alpha) / (sqrt(2.0) * f))) -
         f / sqrt(2 * M_PI) *
             (exp(-pow(t - R / beta, 2) / (2 * f * f)) -
              exp(-pow(t - R / alpha, 2) / (2 * f * f)));
  //  temp = 1/(f*sqrt(2*M_PI))*(
  //  f*f*(-exp(-pow(t-R/beta,2)/(2*f*f))+exp(-pow(t-R/alpha,2)/(2*f*f)) ) +
  //	     t*0.5*sqrt(M_PI*2)*f*( erf((t-R/alpha)/(sqrt(2.0)*f)) -
  // erf((t-R/beta)/(sqrt(2.0)*f)) ) );
  //  temp = 1 /(f*sqrt(2*M_PI))*(f*( (-exp(-pow(t-R / alpha,2)/pow(f,2)) +
  //  exp(-pow(t-R / beta,2)/pow(f,2)) )*f + sqrt(M_PI)*t*(-erf((t-R / alpha) /
  //  f) + erf(R / beta / f))))/2.;
  return temp;
}

//-----------------------------------------------------------------------
// void EW::get_exact_point_source( Sarray& u, float_sw4 t, int g, Source&
// source )
#ifdef NO_DEVICE_FUNCTION_POINTERS
void EW::get_exact_point_source(float_sw4* up, float_sw4 t, int g,
                                Source& source, int* wind) {}
#else
void EW::get_exact_point_source(float_sw4* up, float_sw4 t, int g,
                                Source& source, int* wind) {
  SW4_MARK_FUNCTION;
  // If wind is given, it is assumed that wind is the declared size of up. If
  // not given (wind=0), it is assumed that up is the size of the local
  // processor arrays.
  timeDep tD;
  if (!(source.getName() == "SmoothWave" ||
        source.getName() == "VerySmoothBump" ||
        source.getName() == "C6SmoothBump" || source.getName() == "Gaussian")) {
    cout << "EW::get_exact_point_source: Error, time dependency must be "
            "SmoothWave, VerySmoothBump, C6SmoothBump, or Gaussian, not "
         << source.getName() << endl;
    return;
  } else if (source.getName() == "SmoothWave")
    tD = iSmoothWave;
  else if (source.getName() == "VerySmoothBump")
    tD = iVerySmoothBump;
  else if (source.getName() == "C6SmoothBump")
    tD = iC6SmoothBump;
  else
    tD = iGaussian;

  //   u.set_to_zero();
  float_sw4 alpha = m_point_source_test->m_cp;
  float_sw4 beta = m_point_source_test->m_cs;
  float_sw4 rho = m_point_source_test->m_rho;

  float_sw4 x0 = source.getX0();
  float_sw4 y0 = source.getY0();
  float_sw4 z0 = source.getZ0();
  float_sw4 fr = source.getFrequency();
  float_sw4 time = (t - source.getOffset()) * source.getFrequency();
  if (tD == iGaussian) {
    fr = 1 / fr;
    time = time * fr;
  }
  bool ismomentsource = source.isMomentSource();
  float_sw4 fx, fy, fz;
  float_sw4 mxx, myy, mzz, mxy, mxz, myz, m0;

  if (!ismomentsource) {
    source.getForces(fx, fy, fz);
  } else {
    source.getMoments(mxx, mxy, mxz, myy, myz, mzz);
    //      m0  = source.getAmplitude();
    m0 = 1;
  }
  // bool curvilinear = topographyExists() && g == mNumberOfGrids - 1;
  bool curvilinear = topographyExists() && g > mNumberOfCartesianGrids - 1;
  // std::cout<<"CUYVT"<<curvilinear<<"\n";
  //   float_sw4* up = u.c_ptr();
  float_sw4 h = mGridSize[g];
  float_sw4 eps = 1e-3 * h;
  //   size_t ind = 0;
  int imax, imin, jmax, jmin, kmax, kmin;
  if (wind == 0) {
    imin = m_iStart[g];
    imax = m_iEnd[g];
    jmin = m_jStart[g];
    jmax = m_jEnd[g];
    kmin = m_kStart[g];
    kmax = m_kEnd[g];
  } else {
    imin = wind[0];
    imax = wind[1];
    jmin = wind[2];
    jmax = wind[3];
    kmin = wind[4];
    kmax = wind[5];
  }
  // Note: Use of ind, assumes loop is over the domain over which u is defined.
  //   for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
  //      for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
  //	 for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )

  size_t ni = (imax - imin + 1);
  size_t nij = ni * (jmax - jmin + 1);
  float_sw4 m_zming = m_zmin[g];

  ASSERT_MANAGED(up);

  SView& mXV = mX[g].getview();
  SView& mYV = mY[g].getview();
  SView& mZV = mZ[g].getview();
  float_sw4 alpha2 = alpha * alpha;
  float_sw4 ralpha2 = 1.0 / alpha2;
  float_sw4 alpha3 = alpha2 * alpha;
  float_sw4 beta2 = beta * beta;
  float_sw4 beta3 = beta2 * beta;
  float_sw4 rbeta2 = 1.0 / beta2;

  RAJA::RangeSegment k_range(kmin, kmax + 1);
  RAJA::RangeSegment j_range(jmin, jmax + 1);
  RAJA::RangeSegment i_range(imin, imax + 1);
  // std::cout<<"Size "<<(kmax+1-kmin)*(jmax+1-jmin)*(imax+1-imin)<<"\n";
  SW4_MARK_BEGIN("get_exact_point_source::loop");
  RAJA::kernel<GEPS_POL>(
      RAJA::make_tuple(k_range, j_range, i_range),
      [=] RAJA_DEVICE(int k, int j, int i) {
        // #pragma omp parallel for
        //   for( int k=kmin ; k <= kmax ; k++ )
        //      for( int j=jmin ; j <= jmax ; j++ )
        // 	 for( int i=imin ; i <= imax ; i++ )
        // 	 {
        // 255 rpt for the routine
        // 184 is all the branches (cirvilinear and !ismomensource are commented
        // out
        //
        size_t ind = (i - imin) + ni * (j - jmin) + nij * (k - kmin);
        float_sw4 x, y, z;
        if (curvilinear) {
          x = mXV(i, j, k);
          y = mYV(i, j, k);
          z = mZV(i, j, k);
        } else {
          x = (i - 1) * h;
          y = (j - 1) * h;
          z = (k - 1) * h + m_zming;
        }
        float_sw4 xx0 = x - x0;
        float_sw4 yy0 = y - y0;
        float_sw4 zz0 = z - z0;

        float_sw4 R = sqrt(xx0 * xx0 + yy0 * yy0 + zz0 * zz0);
        float_sw4 RR = 1.0 / R;
        float_sw4 R2 = R * R;
        float_sw4 R3 = R2 * R;
        float_sw4 RR3 = 1.0 / R3;
        float_sw4 R5 = R2 * R3;
        float_sw4 RR5 = 1.0 / R5;
        float_sw4 R7 = R5 * R2;
        float_sw4 RR7 = 1.0 / R7;
        float_sw4 frR2 = fr * fr * R * R;
        if (!ismomentsource) {
          if (R < eps)
            up[3 * ind] = up[3 * ind + 1] = up[3 * ind + 2] = 0;
          else {
            float_sw4 A, B;
            if (tD == iSmoothWave) {
              A = (1 / alpha2 * SmoothWave(time, fr * R, alpha) -
                   1 / beta2 * SmoothWave(time, fr * R, beta) +
                   3 / frR2 *
                       SmoothWave_x_T_Integral(time, fr * R, alpha, beta)) /
                  (4 * M_PI * rho * R * R * R);

              B = (1 / beta2 * SmoothWave(time, fr * R, beta) -
                   1 / frR2 *
                       SmoothWave_x_T_Integral(time, fr * R, alpha, beta)) /
                  (4 * M_PI * rho * R);
            } else if (tD == iVerySmoothBump) {
              A = (1 / alpha2 * VerySmoothBump(time, fr * R, alpha) -
                   1 / beta2 * VerySmoothBump(time, fr * R, beta) +
                   3 / frR2 *
                       VerySmoothBump_x_T_Integral(time, fr * R, alpha, beta)) /
                  (4 * M_PI * rho * R * R * R);

              B = (1 / beta2 * VerySmoothBump(time, fr * R, beta) -
                   1 / frR2 *
                       VerySmoothBump_x_T_Integral(time, fr * R, alpha, beta)) /
                  (4 * M_PI * rho * R);
            } else if (tD == iC6SmoothBump) {
              A = (1 / alpha2 * C6SmoothBump(time, fr * R, alpha) -
                   1 / beta2 * C6SmoothBump(time, fr * R, beta) +
                   3 / frR2 *
                       C6SmoothBump_x_T_Integral(time, fr * R, alpha, beta)) /
                  (4 * M_PI * rho * R * R * R);

              B = (1 / beta2 * C6SmoothBump(time, fr * R, beta) -
                   1 / frR2 *
                       C6SmoothBump_x_T_Integral(time, fr * R, alpha, beta)) /
                  (4 * M_PI * rho * R);
            } else if (tD == iGaussian) {
              A = (1 / alpha2 * Gaussian(time, R, alpha, fr) -
                   1 / beta2 * Gaussian(time, R, beta, fr) +
                   3 / R2 * Gaussian_x_T_Integral(time, R, fr, alpha, beta)) /
                  (4 * M_PI * rho * R * R * R);

              B = (1 / beta2 * Gaussian(time, R, beta, fr) -
                   1 / R2 * Gaussian_x_T_Integral(time, R, fr, alpha, beta)) /
                  (4 * M_PI * rho * R);
            }
            up[3 * ind] =
                (xx0 * xx0 * fx + xx0 * yy0 * fy + xx0 * zz0 * fz) * A + fx * B;
            up[3 * ind + 1] =
                (yy0 * xx0 * fx + yy0 * yy0 * fy + yy0 * zz0 * fz) * A + fy * B;
            up[3 * ind + 2] =
                (zz0 * xx0 * fx + zz0 * yy0 * fy + zz0 * zz0 * fz) * A + fz * B;
          }
        } else {
          up[3 * ind] = up[3 * ind + 1] = up[3 * ind + 2] = 0;
          // Here, ismomentsource == true
          // float_sw4 R = sqrt( xx0*xx0 + yy0*yy0 + zz0*zz0 );
          //     float_sw4 R2 = R*R;
          // float_sw4 R3 = R2*R;
          // float_sw4 R5 = R2*R3;
          // float_sw4 R7 = R5*R2;
          // float_sw4 frR2 = fr*fr*R*R;
          if (R < eps) {
            up[3 * ind] = up[3 * ind + 1] = up[3 * ind + 2] = 0;
          } else {
            float_sw4 A, B, C, D, E;
            if (tD == iSmoothWave) {
              A = SmoothWave(time, R, alpha);
              B = SmoothWave(time, R, beta);
              C = SmoothWave_x_T_Integral(time, R, alpha, beta);
              D = d_SmoothWave_dt(time, R, alpha) / alpha3 * RR;
              E = d_SmoothWave_dt(time, R, beta) / beta3 * RR;
            } else if (tD == iVerySmoothBump) {
              A = VerySmoothBump(time, R, alpha);
              B = VerySmoothBump(time, R, beta);
              C = VerySmoothBump_x_T_Integral(time, R, alpha, beta);
              D = d_VerySmoothBump_dt(time, R, alpha) / alpha3 * RR;
              E = d_VerySmoothBump_dt(time, R, beta) / beta3 * RR;
            } else if (tD == iC6SmoothBump) {
              A = C6SmoothBump(time, R, alpha);
              B = C6SmoothBump(time, R, beta);
              C = C6SmoothBump_x_T_Integral(time, R, alpha, beta);
              D = d_C6SmoothBump_dt(time, R, alpha) / alpha3 * RR;
              E = d_C6SmoothBump_dt(time, R, beta) / beta3 * RR;
            } else if (tD == iGaussian) {
              A = Gaussian(time, R, alpha, fr);
              B = Gaussian(time, R, beta, fr);
              C = Gaussian_x_T_Integral(time, R, fr, alpha, beta);
              D = d_Gaussian_dt(time, R, alpha, fr) / alpha3 * RR;
              E = d_Gaussian_dt(time, R, beta, fr) / beta3 * RR;
            }
            float_sw4 Aalpha2 = A * ralpha2;
            float_sw4 Bbeta2 = B * rbeta2;
            up[3 * ind] +=
                // m_xx*G_xx,x
                +m0 * mxx / (4 * M_PI * rho) *
                (+3 * xx0 * xx0 * xx0 * RR5 * (Aalpha2 - Bbeta2)

                 - 2 * xx0 * RR3 * (Aalpha2 - Bbeta2)

                 + 3 * xx0 * xx0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + (15 * xx0 * xx0 * xx0 * RR7 - 6 * xx0 * RR5) * C

                 + xx0 * xx0 * RR3 * (xx0 * D - xx0 * E)

                 - 1 * RR3 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 - 3 * xx0 * RR5 * C

                 + xx0 / (R3 * beta2) * B

                 + RR * xx0 * E);
            up[3 * ind] +=
                // m_yy*G_xy,y
                +m0 * myy / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 - xx0 * RR3 * (Aalpha2 - Bbeta2)

                 + xx0 * yy0 * RR3 * (yy0 * D - yy0 * E)

                 + 3 * xx0 * yy0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + (15 * xx0 * yy0 * yy0 * RR7 - 3 * xx0 * RR5) * C);
            up[3 * ind] +=
                // m_zz*G_xz,z
                +m0 * mzz / (4 * M_PI * rho) *
                (+3 * xx0 * zz0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 - xx0 * RR3 * (Aalpha2 - Bbeta2)

                 + xx0 * zz0 * RR3 * (zz0 * D - zz0 * E)

                 + 3 * xx0 * zz0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + (15 * xx0 * zz0 * zz0 * RR7 - 3 * xx0 * RR5) * C);
            up[3 * ind] +=
                // m_xy*G_xy,x
                +m0 * mxy / (4 * M_PI * rho) *
                (+3 * xx0 * xx0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 - yy0 * RR3 * (Aalpha2 - Bbeta2)

                 + xx0 * yy0 * RR3 * (xx0 * D - xx0 * E)

                 + 3 * xx0 * yy0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + (15 * xx0 * xx0 * yy0 * RR7 - 3 * yy0 * RR5) * C);
            up[3 * ind] +=
                // m_xy*G_xx,y
                +m0 * mxy / (4 * M_PI * rho) *
                (+3 * xx0 * xx0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 + 3 * xx0 * xx0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + 15 * xx0 * xx0 * yy0 * RR7 * C

                 + xx0 * xx0 * RR3 * (yy0 * D - yy0 * E)

                 - RR3 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 - 3 * yy0 * RR5 * C

                 + yy0 / (R3 * beta2) * B

                 + RR * yy0 * E);
            up[3 * ind] +=
                // m_xz*G_xz,x
                +m0 * mxz / (4 * M_PI * rho) *
                (+3 * xx0 * xx0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 - zz0 * RR3 * (Aalpha2 - Bbeta2)

                 + xx0 * zz0 * RR3 * (xx0 * D - xx0 * E)

                 + 3 * xx0 * zz0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + (15 * xx0 * xx0 * zz0 * RR7 - 3 * zz0 * RR5) * C);
            up[3 * ind] +=
                // m_yz*G_xz,y
                +m0 * myz / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + xx0 * zz0 * RR3 * (yy0 * D - yy0 * E)

                 + 3 * xx0 * zz0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + 15 * xx0 * yy0 * zz0 * RR7 * C);
            up[3 * ind] +=
                // m_xz*G_xx,z
                +m0 * mxz / (4 * M_PI * rho) *
                (+3 * xx0 * xx0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + 3 * xx0 * xx0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + 15 * xx0 * xx0 * zz0 * RR7 * C

                 + xx0 * xx0 * RR3 * (zz0 * D - zz0 * E)

                 - 1 * RR3 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 - 3 * zz0 * RR5 * C

                 + zz0 / (R3 * beta2) * B

                 + RR * zz0 * E);
            up[3 * ind] +=
                // m_yz*G_yx,z
                +m0 * myz / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + xx0 * yy0 * RR3 * (zz0 * D - zz0 * E)

                 + 3 * xx0 * yy0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + 15 * xx0 * yy0 * zz0 * RR7 * C);
            //------------------------------------------------------------
            up[3 * ind + 1] +=
                // m_xx*G_xy,x
                m0 * mxx / (4 * M_PI * rho) *
                (+3 * xx0 * xx0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 - yy0 * RR3 * (Aalpha2 - Bbeta2)

                 + xx0 * yy0 * RR3 * (xx0 * D - xx0 * E)

                 + 3 * xx0 * yy0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + (15 * xx0 * xx0 * yy0 * RR7 - 3 * yy0 * RR5) * C);
            up[3 * ind + 1] +=
                // m_yy**G_yy,y
                +m0 * myy / (4 * M_PI * rho) *
                (+3 * yy0 * yy0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 - 2 * yy0 * RR3 * (Aalpha2 - Bbeta2)

                 + 3 * yy0 * yy0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + (15 * yy0 * yy0 * yy0 * RR7 - 6 * yy0 * RR5) * C

                 + yy0 * yy0 * RR3 * (yy0 * D - yy0 * E)

                 - RR3 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 - 3 * yy0 * RR5 * C

                 + yy0 / (R3 * beta2) * B

                 + 1 * RR * yy0 * E);
            up[3 * ind + 1] +=
                // m_zz*G_zy,z
                +m0 * mzz / (4 * M_PI * rho) *
                (+3 * zz0 * zz0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 - yy0 * RR3 * (Aalpha2 - Bbeta2)

                 + zz0 * yy0 * RR3 * (zz0 * D - zz0 * E)

                 + 3 * zz0 * yy0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + (15 * zz0 * zz0 * yy0 * RR7 - 3 * yy0 * RR5) * C);
            up[3 * ind + 1] +=
                // m_xy*G_yy,x
                +m0 * mxy / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 + 3 * yy0 * yy0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + 15 * xx0 * yy0 * yy0 * RR7 * C

                 + yy0 * yy0 * RR3 * (xx0 * D - xx0 * E)

                 - RR3 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 - 3 * xx0 * RR5 * C

                 + xx0 / (R3 * beta2) * B

                 + RR * xx0 * E);
            up[3 * ind + 1] +=
                // m_xz*G_zy,x
                +m0 * mxz / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + yy0 * zz0 * RR3 * (xx0 * D - xx0 * E)

                 + 3 * yy0 * zz0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + 15 * xx0 * yy0 * zz0 * RR7 * C);
            up[3 * ind + 1] +=
                // m_xy*G_xy,y
                +m0 * mxy / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 - xx0 * RR3 * (Aalpha2 - Bbeta2)

                 + xx0 * yy0 * RR3 * (yy0 * D - yy0 * E)

                 + 3 * xx0 * yy0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + (15 * xx0 * yy0 * yy0 * RR7 - 3 * xx0 * RR5) * C);
            up[3 * ind + 1] +=
                // m_yz*G_zy,y
                +m0 * myz / (4 * M_PI * rho) *
                (+3 * zz0 * yy0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 - zz0 * RR3 * (Aalpha2 - Bbeta2)

                 + zz0 * yy0 * RR3 * (yy0 * D - yy0 * E)

                 + 3 * zz0 * yy0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + (15 * zz0 * yy0 * yy0 * RR7 - 3 * zz0 * RR5) * C);
            up[3 * ind + 1] +=
                // m_xz*G_xy,z
                +m0 * mxz / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + xx0 * yy0 * RR3 * (zz0 * D - zz0 * E)

                 + 3 * xx0 * yy0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + 15 * xx0 * yy0 * zz0 * RR7 * C);
            up[3 * ind + 1] +=
                // m_yz*G_yy,z
                +m0 * myz / (4 * M_PI * rho) *
                (+3 * zz0 * yy0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 + 3 * yy0 * yy0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + 15 * zz0 * yy0 * yy0 * RR7 * C

                 + yy0 * yy0 * RR3 * (zz0 * D - zz0 * E)

                 - RR3 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 - 3 * zz0 * RR5 * C

                 + zz0 / (R3 * beta2) * B

                 + RR * zz0 * E);
            //------------------------------------------------------------
            up[3 * ind + 2] +=
                // m_xx*G_zx,x
                +m0 * mxx / (4 * M_PI * rho) *
                (+3 * xx0 * xx0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 - zz0 * RR3 * (Aalpha2 - Bbeta2)

                 + xx0 * zz0 * RR3 * (xx0 * D - xx0 * E)

                 + 3 * xx0 * zz0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + (15 * xx0 * xx0 * zz0 * RR7 - 3 * zz0 * RR5) * C);
            up[3 * ind + 2] +=
                // m_yy*G_zy,y
                +m0 * myy / (4 * M_PI * rho) *
                (+3 * yy0 * yy0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 - zz0 * RR3 * (Aalpha2 - Bbeta2)

                 + yy0 * zz0 * RR3 * (yy0 * D - yy0 * E)

                 + 3 * yy0 * zz0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + (15 * yy0 * yy0 * zz0 * RR7 - 3 * zz0 * RR5) * C);
            up[3 * ind + 2] +=
                // m_zz**G_zz,z
                +m0 * mzz / (4 * M_PI * rho) *
                (+3 * zz0 * zz0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 - 2 * zz0 * RR3 * (Aalpha2 - Bbeta2)

                 + 3 * zz0 * zz0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + (15 * zz0 * zz0 * zz0 * RR7 - 6 * zz0 * RR5) * C

                 + zz0 * zz0 * RR3 * (zz0 * D - zz0 * E)

                 - RR3 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 - 3 * zz0 * RR5 * C

                 + zz0 / (R3 * beta2) * B

                 + RR * zz0 * E);
            up[3 * ind + 2] +=
                // m_xy*G_zy,x
                +m0 * mxy / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + yy0 * zz0 * RR3 * (xx0 * D - xx0 * E)

                 + 3 * yy0 * zz0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + 15 * xx0 * yy0 * zz0 * RR7 * C);
            up[3 * ind + 2] +=
                // m_xz**G_zz,x
                +m0 * mxz / (4 * M_PI * rho) *
                (+3 * xx0 * zz0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + 3 * zz0 * zz0 * RR5 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 + 15 * xx0 * zz0 * zz0 * RR7 * C

                 + zz0 * zz0 * RR3 * (xx0 * D - xx0 * E)

                 - RR3 * (xx0 * Aalpha2 - xx0 * Bbeta2)

                 - 3 * xx0 * RR5 * C

                 + xx0 / (R3 * beta2) * B

                 + RR * xx0 * E);
            up[3 * ind + 2] +=
                // m_xy*G_xz,y
                +m0 * mxy / (4 * M_PI * rho) *
                (+3 * xx0 * yy0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + xx0 * zz0 * RR3 * (yy0 * D - yy0 * E)

                 + 3 * xx0 * zz0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + 15 * xx0 * yy0 * zz0 * RR7 * C);
            up[3 * ind + 2] +=
                // m_yz*G_zz,y
                +m0 * myz / (4 * M_PI * rho) *
                (+3 * yy0 * zz0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 + 3 * zz0 * zz0 * RR5 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 + 15 * yy0 * zz0 * zz0 * RR7 * C

                 + zz0 * zz0 * RR3 * (yy0 * D - yy0 * E)

                 - RR3 * (yy0 * Aalpha2 - yy0 * Bbeta2)

                 - 3 * yy0 * RR5 * C

                 + yy0 / (R3 * beta2) * B

                 + RR * yy0 * E);
            up[3 * ind + 2] +=
                // m_xz*G_xz,z
                +m0 * mxz / (4 * M_PI * rho) *
                (+3 * xx0 * zz0 * zz0 * RR5 * (Aalpha2 - Bbeta2)

                 - xx0 * RR3 * (Aalpha2 - Bbeta2)

                 + xx0 * zz0 * RR3 * (zz0 * D - zz0 * E)

                 + 3 * xx0 * zz0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + (15 * xx0 * zz0 * zz0 * RR7 - 3 * xx0 * RR5) * C);
            up[3 * ind + 2] +=
                // m_yz*G_yz,z
                +m0 * myz / (4 * M_PI * rho) *
                (+3 * zz0 * zz0 * yy0 * RR5 * (Aalpha2 - Bbeta2)

                 - yy0 * RR3 * (Aalpha2 - Bbeta2)

                 + zz0 * yy0 * RR3 * (zz0 * D - zz0 * E)

                 + 3 * zz0 * yy0 * RR5 * (zz0 * Aalpha2 - zz0 * Bbeta2)

                 + (15 * zz0 * zz0 * yy0 * RR7 - 3 * yy0 * RR5) * C);
          }
        }
        //	    ind++;
      });
  SYNC_STREAM;
  SW4_MARK_END("get_exact_point_source::loop");
}
#endif

#include <cmath>
#include <complex>

//-----------------------------------------------------------------------
complex<float_sw4> asin(complex<float_sw4> z) {
  complex<float_sw4> I(0, 1);
  return -I * log(I * z + sqrt(1. - pow(z, 2)));
}

//-----------------------------------------------------------------------
complex<float_sw4> atan(complex<float_sw4> z) {
  complex<float_sw4> I(0, 1);
  return I / 2. * log((I + z) / (I - z));
}

//-----------------------------------------------------------------------
complex<double> atan2(complex<double> z, complex<double> w) {
  complex<double> I(0, 1);
  complex<double> Zero(0, 0);

  if (w == Zero) {
    if (z.real() > 0)
      return M_PI / 2.;
    else
      return -M_PI / 2.;
  } else {
    complex<double> retval = I / 2. * log((I + z / w) / (I - z / w));
    if (retval.real() < 0 && z.real() > 0) retval = retval + M_PI;
    if (retval.real() > 0 && z.real() < 0) retval = retval - M_PI;
    return retval;
    //      return I/2.*log((I + z/w)/(I - z/w));
  }
}

//-----------------------------------------------------------------------
void EW::get_exact_lamb2(vector<Sarray>& a_U, float_sw4 a_t, Source& a_source) {
  SW4_MARK_FUNCTION;
  // initialize
  //  for (int g=0; g<mNumberOfGrids; g++)
  //    a_U[g].set_to_zero();
  double x0 = a_source.getX0();
  double y0 = a_source.getY0();
  double z0 = a_source.getZ0();

  double fx, fy, fz;
  a_source.getForces(fx, fy, fz);
  double cs = m_lamb_test->m_cs;
  double mu = m_lamb_test->m_mu;
  int g = mNumberOfCartesianGrids - 1;  // top Cartesian grid
  double h = mGridSize[g];
  int ifirst = m_iStart[g];
  int ilast = m_iEnd[g];
  int jfirst = m_jStart[g];
  int jlast = m_jEnd[g];
  int kfirst = m_kStart[g];
  int klast = m_kEnd[g];
  int tfun = 0;
  if (a_source.getTfunc() == iVerySmoothBump)
    tfun = 1;
  else if (a_source.getTfunc() == iC6SmoothBump)
    tfun = 2;
  // Fortran
  size_t npts = a_U[g].m_npts;
  float_sw4* uexact = SW4_NEW(Space::Managed, float_sw4[npts]);
  for (size_t i = 0; i < npts; i++) uexact[i] = 0;
  lambexact(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, uexact, &a_t,
            &mu, &cs, &x0, &y0, &fz, &h, &tfun);
  //	     a_U[g].c_ptr(), &a_t, &mu, &cs, &x0, &y0, &fz, &h, &tfun );
  a_U[g].assign(uexact, 0);
  ::operator delete[](uexact, Space::Managed);
  // test: output uz in one point
  // int i0=176, j0=151, k0=1;
  // if (m_iStart[g] <= i0 && i0 <= m_iEnd[g] && m_jStart[g] <= j0 && j0 <=
  // m_jEnd[g])
  //   printf("Lambexact: t=%e, uze=%e\n", a_t, a_U[g](3,i0,j0,k0));
}

//-----------------------------------------------------------------------
void EW::get_exact_lamb(vector<Sarray>& a_U, float_sw4 a_t, Source& a_source) {
  SW4_MARK_FUNCTION;
  int g;

  // initialize
  for (g = 0; g < mNumberOfGrids; g++) a_U[g].set_to_zero();

  double h, t = a_t;

  double gamma = sqrt(3. + sqrt(3.)) / 2.;

  double alpha = m_lamb_test->m_cp;
  double beta = m_lamb_test->m_cs;
  double mu = m_lamb_test->m_mu;

  double x0 = a_source.getX0();
  double y0 = a_source.getY0();
  double z0 = a_source.getZ0();

  double fx, fy, fz;
  a_source.getForces(fx, fy, fz);

  // Only the z-component of solution on the flat surface (z=0) is known by this
  // routine
  int k = 1;

  g = mNumberOfCartesianGrids - 1;  // top Cartesian grid
  h = mGridSize[g];
  // z = 0.0;

// loop over all points in the horizontal plane
#pragma omp parallel for
  for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
    for (int i = m_iStart[g]; i <= m_iEnd[g]; i++) {
      double x = (i - 1) * h;
      double y = (j - 1) * h;
      double uz = 0;

      double R = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
      if (R < h) {
        uz = 0;
      } else {
        uz = 0;
        if (t <= R / alpha) {
          uz = 0.0;
        } else {
          double tau = t * beta / R;
          double r = R;
          if (tau > gamma) {
            uz += G4_Integral(std::min(std::max(0.0, tau - gamma), beta / r),
                              tau, r, beta) -
                  G4_Integral(0.0, tau, r, beta);
          }
          if (tau > 1 && tau < beta / r + gamma) {
            uz += G3_Integral(std::min(tau - 1, beta / r), tau, r, beta) -
                  G3_Integral(std::max(0.0, tau - gamma), tau, r, beta);
          }
          if (tau > 1 / sqrt(3.) && tau < beta / r + 1) {
            uz += G2_Integral(std::min(tau - 1 / sqrt(3.), beta / r), tau, r,
                              beta) -
                  G2_Integral(std::max(tau - 1, 0.0), tau, r, beta);
          }
          uz *= -fz / (M_PI * M_PI * mu) * alpha * alpha / (beta * beta * beta);
        }
      }  // end if R<h
         // assign Sarray
      a_U[g](3, i, j, k) = uz;
    }  // end for i,j
}  // end get_exact_lamb()

//-----------------------------------------------------------------------
double EW::G4_Integral(double T, double t, double r, double beta) {
  double c0 = 1024., c1 = -5120., c2 = 10240., c3 = -10240., c4 = 5120.,
         c5 = -1024.;

  return -(M_PI * ((c5 * pow(r, 9) * pow(T, 10)) / pow(beta, 9) +
                   (c4 * pow(r, 8) * pow(T, 9)) / pow(beta, 8) +
                   (c3 * pow(r, 7) * pow(T, 8)) / pow(beta, 7) +
                   (c2 * pow(r, 6) * pow(T, 7)) / pow(beta, 6) +
                   (c1 * pow(r, 5) * pow(T, 6)) / pow(beta, 5) +
                   (c0 * pow(r, 4) * pow(T, 5)) / pow(beta, 4))) /
         8.;
}

//-----------------------------------------------------------------------
double EW::G3_Integral(double iT, double it, double ir, double ibeta) {
  complex<double> T = iT, t = it, r = ir, beta = ibeta;
  complex<double> c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120,
                  c5 = -1024;
  complex<double> gamma = sqrt(3. + sqrt(3.)) / 2.;
  complex<double> tmp;

  tmp = -(M_PI * ((c5 * pow(r, 9) * pow(T, 10)) / pow(beta, 9) +
                  (c4 * pow(r, 8) * pow(T, 9)) / pow(beta, 8) +
                  (c3 * pow(r, 7) * pow(T, 8)) / pow(beta, 7) +
                  (c2 * pow(r, 6) * pow(T, 7)) / pow(beta, 6) +
                  (c1 * pow(r, 5) * pow(T, 6)) / pow(beta, 5) +
                  (c0 * pow(r, 4) * pow(T, 5)) / pow(beta, 4))) /
        8.;

  tmp +=
      (sqrt(5. + 3. * sqrt(3.)) * M_PI * pow(r, 4) *
       (-(sqrt(-pow(t, 2) + 2. * t * T - pow(T, 2) + pow(gamma, 2)) *
          (10. * c5 * pow(r, 5) *
               (114064. * pow(t, 8) + 73744. * pow(t, 7) * T +
                8. * pow(t, 6) * (6698. * pow(T, 2) + 150373. * pow(gamma, 2)) +
                8. * pow(t, 5) *
                    (5018. * pow(T, 3) + 68871. * T * pow(gamma, 2)) +
                2. * pow(t, 4) *
                    (15032. * pow(T, 4) + 139272. * pow(T, 2) * pow(gamma, 2) +
                     961437. * pow(gamma, 4)) +
                2. * pow(t, 3) *
                    (11000. * pow(T, 5) + 68536. * pow(T, 3) * pow(gamma, 2) +
                     284361. * T * pow(gamma, 4)) +
                pow(t, 2) *
                    (15280. * pow(T, 6) + 61032. * pow(T, 4) * pow(gamma, 2) +
                     170190. * pow(T, 2) * pow(gamma, 4) +
                     572519. * pow(gamma, 6)) +
                t * (9520. * pow(T, 7) + 22200. * pow(T, 5) * pow(gamma, 2) +
                     41574. * pow(T, 3) * pow(gamma, 4) +
                     82841. * T * pow(gamma, 6)) +
                128. *
                    (35. * pow(T, 8) + 40. * pow(T, 6) * pow(gamma, 2) +
                     48. * pow(T, 4) * pow(gamma, 4) +
                     64. * pow(T, 2) * pow(gamma, 6) + 128. * pow(gamma, 8))) +
           3. * beta *
               (9. * c4 * pow(r, 4) *
                    (36528. * pow(t, 7) + 23088. * pow(t, 6) * T +
                     24. * pow(t, 5) *
                         (682. * pow(T, 2) + 11989. * pow(gamma, 2)) +
                     8. * pow(t, 4) *
                         (1486. * pow(T, 3) + 15333. * T * pow(gamma, 2)) +
                     pow(t, 3) * (8528. * pow(T, 4) +
                                  56496. * pow(T, 2) * pow(gamma, 2) +
                                  305934. * pow(gamma, 4)) +
                     pow(t, 2) * (5840. * pow(T, 5) +
                                  24272. * pow(T, 3) * pow(gamma, 2) +
                                  75798. * T * pow(gamma, 4)) +
                     t * (3600. * pow(T, 6) +
                          8632. * pow(T, 4) * pow(gamma, 2) +
                          17226. * pow(T, 2) * pow(gamma, 4) +
                          45477. * pow(gamma, 6)) +
                     35. * (48. * pow(T, 7) + 56. * pow(T, 5) * pow(gamma, 2) +
                            70. * pow(T, 3) * pow(gamma, 4) +
                            105. * T * pow(gamma, 6))) +
                8. * beta *
                    (8. * c3 * pow(r, 3) *
                         (4356. * pow(t, 6) + 2676. * pow(t, 5) * T +
                          12. * pow(t, 4) *
                              (153. * pow(T, 2) + 2033. * pow(gamma, 2)) +
                          4. * pow(t, 3) *
                              (319. * pow(T, 3) + 2358. * T * pow(gamma, 2)) +
                          pow(t, 2) * (856. * pow(T, 4) +
                                       3786. * pow(T, 2) * pow(gamma, 2) +
                                       15525. * pow(gamma, 4)) +
                          t * (520. * pow(T, 5) +
                               1298. * pow(T, 3) * pow(gamma, 2) +
                               2907. * T * pow(gamma, 4)) +
                          48. *
                              (5. * pow(T, 6) + 6. * pow(T, 4) * pow(gamma, 2) +
                               8. * pow(T, 2) * pow(gamma, 4) +
                               16. * pow(gamma, 6))) +
                     7. * beta *
                         (2. * beta *
                              (25. * c0 * beta *
                                   (50. * pow(t, 3) + 26. * pow(t, 2) * T +
                                    14. * t * pow(T, 2) + 6. * pow(T, 3) +
                                    55. * t * pow(gamma, 2) +
                                    9. * T * pow(gamma, 2)) +
                               6. * c1 * r *
                                   (274. * pow(t, 4) + 154. * pow(t, 3) * T +
                                    94. * pow(t, 2) * pow(T, 2) +
                                    54. * t * pow(T, 3) + 24. * pow(T, 4) +
                                    607. * pow(t, 2) * pow(gamma, 2) +
                                    161. * t * T * pow(gamma, 2) +
                                    32. * pow(T, 2) * pow(gamma, 2) +
                                    64. * pow(gamma, 4))) +
                          7. * c2 * pow(r, 2) *
                              (588. * pow(t, 5) + 348. * pow(t, 4) * T +
                               40. * pow(T, 5) +
                               50. * pow(T, 3) * pow(gamma, 2) +
                               75. * T * pow(gamma, 4) +
                               12. * pow(t, 3) *
                                   (19. * pow(T, 2) + 182. * pow(gamma, 2)) +
                               4. * pow(t, 2) *
                                   (37. * pow(T, 3) +
                                    183. * T * pow(gamma, 2)) +
                               t * (88. * pow(T, 4) +
                                    234. * pow(T, 2) * pow(gamma, 2) +
                                    693. * pow(gamma, 4)))))))) /
            40320. -
        ((10. * c5 * pow(r, 5) * t *
              (128. * pow(t, 8) + 2304. * pow(t, 6) * pow(gamma, 2) +
               6048. * pow(t, 4) * pow(gamma, 4) +
               3360. * pow(t, 2) * pow(gamma, 6) + 315. * pow(gamma, 8)) +
          beta *
              (9. * c4 * pow(r, 4) *
                   (128. * pow(t, 8) + 1792. * pow(t, 6) * pow(gamma, 2) +
                    3360. * pow(t, 4) * pow(gamma, 4) +
                    1120. * pow(t, 2) * pow(gamma, 6) + 35. * pow(gamma, 8)) +
               8. * beta *
                   (8. * c3 * pow(r, 3) * t *
                        (16. * pow(t, 6) + 168. * pow(t, 4) * pow(gamma, 2) +
                         210. * pow(t, 2) * pow(gamma, 4) +
                         35. * pow(gamma, 6)) +
                    beta * (7. * c2 * pow(r, 2) *
                                (16. * pow(t, 6) +
                                 120. * pow(t, 4) * pow(gamma, 2) +
                                 90. * pow(t, 2) * pow(gamma, 4) +
                                 5. * pow(gamma, 6)) +
                            2. * beta *
                                (5. * c0 * beta *
                                     (8. * pow(t, 4) +
                                      24. * pow(t, 2) * pow(gamma, 2) +
                                      3. * pow(gamma, 4)) +
                                 6. * c1 * r * t *
                                     (8. * pow(t, 4) +
                                      40. * pow(t, 2) * pow(gamma, 2) +
                                      15. * pow(gamma, 4))))))) *
         atan2((t - T),
               sqrt(-pow(t, 2) + 2. * t * T - pow(T, 2) + pow(gamma, 2)))) /
            128.)) /
      (48. * pow(beta, 9));

  //  cout << "ArcTan(Arg) = " << atan((t - T)/sqrt(-pow(t,2) + 2.*t*T -
  //  pow(T,2) + pow(gamma,2))) << ". Arg = " << (t - T) << "/" <<
  //  sqrt(-pow(t,2) + 2.*t*T - pow(T,2) + pow(gamma,2)) << endl;

  return tmp.real();
}

//-----------------------------------------------------------------------
double EW::G2_Integral(double iT, double it, double ir, double ibeta) {
  complex<double> T = iT, t = it, r = ir, beta = ibeta;
  complex<double> c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120,
                  c5 = -1024;
  complex<double> gamma = sqrt(3. + sqrt(3.)) / 2.;
  complex<double> tmp;

  tmp = (-(M_PI * ((c5 * pow(r, 9) * pow(T, 10)) / pow(beta, 9) +
                   (c4 * pow(r, 8) * pow(T, 9)) / pow(beta, 8) +
                   (c3 * pow(r, 7) * pow(T, 8)) / pow(beta, 7) +
                   (c2 * pow(r, 6) * pow(T, 7)) / pow(beta, 6) +
                   (c1 * pow(r, 5) * pow(T, 6)) / pow(beta, 5) +
                   (c0 * pow(r, 4) * pow(T, 5)) / pow(beta, 4))) /
         8.) /
        2.;

  tmp +=
      ((sqrt(5. + 3. * sqrt(3.)) * M_PI * pow(r, 4) *
        (-(sqrt(-pow(t, 2) + 2. * t * T - pow(T, 2) + pow(gamma, 2)) *
           (10. * c5 * pow(r, 5) *
                (114064. * pow(t, 8) + 73744. * pow(t, 7) * T +
                 8. * pow(t, 6) *
                     (6698. * pow(T, 2) + 150373. * pow(gamma, 2)) +
                 8. * pow(t, 5) *
                     (5018. * pow(T, 3) + 68871. * T * pow(gamma, 2)) +
                 2. * pow(t, 4) *
                     (15032. * pow(T, 4) + 139272. * pow(T, 2) * pow(gamma, 2) +
                      961437. * pow(gamma, 4)) +
                 2. * pow(t, 3) *
                     (11000. * pow(T, 5) + 68536. * pow(T, 3) * pow(gamma, 2) +
                      284361. * T * pow(gamma, 4)) +
                 pow(t, 2) *
                     (15280. * pow(T, 6) + 61032. * pow(T, 4) * pow(gamma, 2) +
                      170190. * pow(T, 2) * pow(gamma, 4) +
                      572519. * pow(gamma, 6)) +
                 t * (9520. * pow(T, 7) + 22200. * pow(T, 5) * pow(gamma, 2) +
                      41574. * pow(T, 3) * pow(gamma, 4) +
                      82841. * T * pow(gamma, 6)) +
                 128. *
                     (35. * pow(T, 8) + 40. * pow(T, 6) * pow(gamma, 2) +
                      48. * pow(T, 4) * pow(gamma, 4) +
                      64. * pow(T, 2) * pow(gamma, 6) + 128. * pow(gamma, 8))) +
            3. * beta *
                (9. * c4 * pow(r, 4) *
                     (36528. * pow(t, 7) + 23088. * pow(t, 6) * T +
                      24. * pow(t, 5) *
                          (682. * pow(T, 2) + 11989. * pow(gamma, 2)) +
                      8. * pow(t, 4) *
                          (1486. * pow(T, 3) + 15333. * T * pow(gamma, 2)) +
                      pow(t, 3) * (8528. * pow(T, 4) +
                                   56496. * pow(T, 2) * pow(gamma, 2) +
                                   305934. * pow(gamma, 4)) +
                      pow(t, 2) * (5840. * pow(T, 5) +
                                   24272. * pow(T, 3) * pow(gamma, 2) +
                                   75798. * T * pow(gamma, 4)) +
                      t * (3600. * pow(T, 6) +
                           8632. * pow(T, 4) * pow(gamma, 2) +
                           17226. * pow(T, 2) * pow(gamma, 4) +
                           45477. * pow(gamma, 6)) +
                      35. * (48. * pow(T, 7) + 56. * pow(T, 5) * pow(gamma, 2) +
                             70. * pow(T, 3) * pow(gamma, 4) +
                             105. * T * pow(gamma, 6))) +
                 8. * beta *
                     (8. * c3 * pow(r, 3) *
                          (4356. * pow(t, 6) + 2676. * pow(t, 5) * T +
                           12. * pow(t, 4) *
                               (153. * pow(T, 2) + 2033. * pow(gamma, 2)) +
                           4. * pow(t, 3) *
                               (319. * pow(T, 3) + 2358. * T * pow(gamma, 2)) +
                           pow(t, 2) * (856. * pow(T, 4) +
                                        3786. * pow(T, 2) * pow(gamma, 2) +
                                        15525. * pow(gamma, 4)) +
                           t * (520. * pow(T, 5) +
                                1298. * pow(T, 3) * pow(gamma, 2) +
                                2907. * T * pow(gamma, 4)) +
                           48. * (5. * pow(T, 6) +
                                  6. * pow(T, 4) * pow(gamma, 2) +
                                  8. * pow(T, 2) * pow(gamma, 4) +
                                  16. * pow(gamma, 6))) +
                      7. * beta *
                          (2. * beta *
                               (25. * c0 * beta *
                                    (50. * pow(t, 3) + 26. * pow(t, 2) * T +
                                     14. * t * pow(T, 2) + 6. * pow(T, 3) +
                                     55. * t * pow(gamma, 2) +
                                     9. * T * pow(gamma, 2)) +
                                6. * c1 * r *
                                    (274. * pow(t, 4) + 154. * pow(t, 3) * T +
                                     94. * pow(t, 2) * pow(T, 2) +
                                     54. * t * pow(T, 3) + 24. * pow(T, 4) +
                                     607. * pow(t, 2) * pow(gamma, 2) +
                                     161. * t * T * pow(gamma, 2) +
                                     32. * pow(T, 2) * pow(gamma, 2) +
                                     64. * pow(gamma, 4))) +
                           7. * c2 * pow(r, 2) *
                               (588. * pow(t, 5) + 348. * pow(t, 4) * T +
                                40. * pow(T, 5) +
                                50. * pow(T, 3) * pow(gamma, 2) +
                                75. * T * pow(gamma, 4) +
                                12. * pow(t, 3) *
                                    (19. * pow(T, 2) + 182. * pow(gamma, 2)) +
                                4. * pow(t, 2) *
                                    (37. * pow(T, 3) +
                                     183. * T * pow(gamma, 2)) +
                                t * (88. * pow(T, 4) +
                                     234. * pow(T, 2) * pow(gamma, 2) +
                                     693. * pow(gamma, 4)))))))) /
             40320. -
         ((10. * c5 * pow(r, 5) * t *
               (128. * pow(t, 8) + 2304. * pow(t, 6) * pow(gamma, 2) +
                6048. * pow(t, 4) * pow(gamma, 4) +
                3360. * pow(t, 2) * pow(gamma, 6) + 315. * pow(gamma, 8)) +
           beta *
               (9. * c4 * pow(r, 4) *
                    (128. * pow(t, 8) + 1792. * pow(t, 6) * pow(gamma, 2) +
                     3360. * pow(t, 4) * pow(gamma, 4) +
                     1120. * pow(t, 2) * pow(gamma, 6) + 35. * pow(gamma, 8)) +
                8. * beta *
                    (8. * c3 * pow(r, 3) * t *
                         (16. * pow(t, 6) + 168. * pow(t, 4) * pow(gamma, 2) +
                          210. * pow(t, 2) * pow(gamma, 4) +
                          35. * pow(gamma, 6)) +
                     beta * (7. * c2 * pow(r, 2) *
                                 (16. * pow(t, 6) +
                                  120. * pow(t, 4) * pow(gamma, 2) +
                                  90. * pow(t, 2) * pow(gamma, 4) +
                                  5. * pow(gamma, 6)) +
                             2. * beta *
                                 (5. * c0 * beta *
                                      (8. * pow(t, 4) +
                                       24. * pow(t, 2) * pow(gamma, 2) +
                                       3. * pow(gamma, 4)) +
                                  6. * c1 * r * t *
                                      (8. * pow(t, 4) +
                                       40. * pow(t, 2) * pow(gamma, 2) +
                                       15. * pow(gamma, 4))))))) *
          atan2((t - T),
                sqrt(-pow(t, 2) + 2. * t * T - pow(T, 2) + pow(gamma, 2)))) /
             128.)) /
       (48. * pow(beta, 9))) /
      2.;

  tmp +=
      -(sqrt(-5. + 3. * sqrt(3.)) * M_PI * pow(r, 4) *
        (sqrt(-0.75 + sqrt(3.) / 4. + pow(t - T, 2)) *
             (640. * pow(-3. + sqrt(3.), 4) * c5 * pow(r, 5) -
              (pow(-3. + sqrt(3.), 3) * pow(r, 3) *
               (10. * c5 * pow(r, 2) *
                    (572519. * pow(t, 2) + 82841. * t * T + 8192. * pow(T, 2)) +
                9. * beta *
                    (9. * c4 * r * (15159. * t + 1225. * T) +
                     16384. * c3 * beta))) /
                  64. +
              (3. * pow(-3. + sqrt(3.), 2) * r *
               (10. * c5 * pow(r, 4) *
                    (320479. * pow(t, 4) + 94787. * pow(t, 3) * T +
                     28365. * pow(t, 2) * pow(T, 2) + 6929. * t * pow(T, 3) +
                     1024. * pow(T, 4)) +
                3. * beta *
                    (3. * c4 * pow(r, 3) *
                         (152967. * pow(t, 3) + 37899. * pow(t, 2) * T +
                          8613. * t * pow(T, 2) + 1225. * pow(T, 3)) +
                     4. * beta *
                         (8. * c3 * pow(r, 2) *
                              (5175. * pow(t, 2) + 969. * t * T +
                               128. * pow(T, 2)) +
                          7. * beta *
                              (7. * c2 * r * (231. * t + 25. * T) +
                               256. * c1 * beta))))) /
                  8. +
              2. * (3. - sqrt(3.)) *
                  (10. * c5 * pow(r, 5) *
                       (150373. * pow(t, 6) + 68871. * pow(t, 5) * T +
                        34818. * pow(t, 4) * pow(T, 2) +
                        17134. * pow(t, 3) * pow(T, 3) +
                        7629. * pow(t, 2) * pow(T, 4) + 2775. * t * pow(T, 5) +
                        640. * pow(T, 6)) +
                   3. * beta *
                       (9. * c4 * pow(r, 4) *
                            (35967. * pow(t, 5) + 15333. * pow(t, 4) * T +
                             7062. * pow(t, 3) * pow(T, 2) +
                             3034. * pow(t, 2) * pow(T, 3) +
                             1079. * t * pow(T, 4) + 245. * pow(T, 5)) +
                        2. * beta *
                            (8. * c3 * pow(r, 3) *
                                 (12198. * pow(t, 4) + 4716. * pow(t, 3) * T +
                                  1893. * pow(t, 2) * pow(T, 2) +
                                  649. * t * pow(T, 3) + 144. * pow(T, 4)) +
                             7. * beta *
                                 (7. * c2 * pow(r, 2) *
                                      (1092. * pow(t, 3) +
                                       366. * pow(t, 2) * T +
                                       117. * t * pow(T, 2) + 25. * pow(T, 3)) +
                                  beta * (6. * c1 * r *
                                              (607. * pow(t, 2) + 161. * t * T +
                                               32. * pow(T, 2)) +
                                          25. * c0 * (55. * t + 9. * T) *
                                              beta))))) +
              16. *
                  (10. * c5 * pow(r, 5) *
                       (7129. * pow(t, 8) + 4609. * pow(t, 7) * T +
                        3349. * pow(t, 6) * pow(T, 2) +
                        2509. * pow(t, 5) * pow(T, 3) +
                        1879. * pow(t, 4) * pow(T, 4) +
                        1375. * pow(t, 3) * pow(T, 5) +
                        955. * pow(t, 2) * pow(T, 6) + 595. * t * pow(T, 7) +
                        280. * pow(T, 8)) +
                   3. * beta *
                       (9. * c4 * pow(r, 4) *
                            (2283. * pow(t, 7) + 1443. * pow(t, 6) * T +
                             1023. * pow(t, 5) * pow(T, 2) +
                             743. * pow(t, 4) * pow(T, 3) +
                             533. * pow(t, 3) * pow(T, 4) +
                             365. * pow(t, 2) * pow(T, 5) +
                             225. * t * pow(T, 6) + 105. * pow(T, 7)) +
                        2. * beta *
                            (8. * c3 * pow(r, 3) *
                                 (1089. * pow(t, 6) + 669. * pow(t, 5) * T +
                                  459. * pow(t, 4) * pow(T, 2) +
                                  319. * pow(t, 3) * pow(T, 3) +
                                  214. * pow(t, 2) * pow(T, 4) +
                                  130. * t * pow(T, 5) + 60. * pow(T, 6)) +
                             7. * beta *
                                 (7. * c2 * pow(r, 2) *
                                      (147. * pow(t, 5) + 87. * pow(t, 4) * T +
                                       57. * pow(t, 3) * pow(T, 2) +
                                       37. * pow(t, 2) * pow(T, 3) +
                                       22. * t * pow(T, 4) + 10. * pow(T, 5)) +
                                  beta * (6. * c1 * r *
                                              (137. * pow(t, 4) +
                                               77. * pow(t, 3) * T +
                                               47. * pow(t, 2) * pow(T, 2) +
                                               27. * t * pow(T, 3) +
                                               12. * pow(T, 4)) +
                                          25. * c0 *
                                              (25. * pow(t, 3) +
                                               13. * pow(t, 2) * T +
                                               7. * t * pow(T, 2) +
                                               3. * pow(T, 3)) *
                                              beta)))))) +
         315. *
             ((315. * pow(-3. + sqrt(3.), 4) * pow(r, 4) *
               (10. * c5 * r * t + c4 * beta)) /
                  256. -
              (35. * pow(-3. + sqrt(3.), 3) * pow(r, 2) *
               (120. * c5 * pow(r, 3) * pow(t, 3) +
                beta * (36. * c4 * pow(r, 2) * pow(t, 2) +
                        beta * (8. * c3 * r * t + c2 * beta)))) /
                  8. -
              90. * (-2. + sqrt(3.)) *
                  (252. * c5 * pow(r, 5) * pow(t, 5) +
                   beta * (126. * c4 * pow(r, 4) * pow(t, 4) +
                           beta * (56. * c3 * pow(r, 3) * pow(t, 3) +
                                   21. * c2 * pow(r, 2) * pow(t, 2) * beta +
                                   6. * c1 * r * t * pow(beta, 2) +
                                   c0 * pow(beta, 3)))) +
              128. * pow(t, 4) *
                  (10. * c5 * pow(r, 5) * pow(t, 5) +
                   beta * (9. * c4 * pow(r, 4) * pow(t, 4) +
                           beta * (8. * c3 * pow(r, 3) * pow(t, 3) +
                                   7. * c2 * pow(r, 2) * pow(t, 2) * beta +
                                   6. * c1 * r * t * pow(beta, 2) +
                                   5. * c0 * pow(beta, 3)))) -
              48. * (-3. + sqrt(3.)) * pow(t, 2) *
                  (120. * c5 * pow(r, 5) * pow(t, 5) +
                   beta * (84. * c4 * pow(r, 4) * pow(t, 4) +
                           beta * (56. * c3 * pow(r, 3) * pow(t, 3) +
                                   35. * c2 * pow(r, 2) * pow(t, 2) * beta +
                                   20. * c1 * r * t * pow(beta, 2) +
                                   10. * c0 * pow(beta, 3))))) *
             log(-t + sqrt(-0.75 + sqrt(3.) / 4. + pow(t - T, 2)) + T))) /
      (3.87072e6 * pow(beta, 9));

  tmp +=
      (M_PI * pow(r, 4) *
       (4. * sqrt(-0.25 + pow(t - T, 2)) *
            (10. * c5 * pow(r, 5) *
                 (7300096. * pow(t, 8) + 4719616. * pow(t, 7) * T +
                  128. * pow(t, 5) * T * (68871. + 20072. * pow(T, 2)) +
                  128. * pow(t, 6) * (150373. + 26792. * pow(T, 2)) +
                  8. * pow(t, 3) * T *
                      (284361. + 274144. * pow(T, 2) + 176000. * pow(T, 4)) +
                  8. * pow(t, 4) *
                      (961437. + 557088. * pow(T, 2) + 240512. * pow(T, 4)) +
                  t * T *
                      (82841. + 166296. * pow(T, 2) + 355200. * pow(T, 4) +
                       609280. * pow(T, 6)) +
                  pow(t, 2) * (572519. + 680760. * pow(T, 2) +
                               976512. * pow(T, 4) + 977920. * pow(T, 6)) +
                  4096. * (1. + 2. * pow(T, 2) + 6. * pow(T, 4) +
                           20. * pow(T, 6) + 70. * pow(T, 8))) +
             3. * beta *
                 (9. * c4 * pow(r, 4) *
                      (2337792. * pow(t, 7) + 1477632. * pow(t, 6) * T +
                       384. * pow(t, 5) * (11989. + 2728. * pow(T, 2)) +
                       128. * pow(t, 4) * T * (15333. + 5944. * pow(T, 2)) +
                       8. * pow(t, 2) * T *
                           (37899. + 48544. * pow(T, 2) + 46720. * pow(T, 4)) +
                       8. * pow(t, 3) *
                           (152967. + 112992. * pow(T, 2) +
                            68224. * pow(T, 4)) +
                       35. * T *
                           (105. + 280. * pow(T, 2) + 896. * pow(T, 4) +
                            3072. * pow(T, 6)) +
                       t * (45477. + 68904. * pow(T, 2) + 138112. * pow(T, 4) +
                            230400. * pow(T, 6))) +
                  32. * beta *
                      (8. * c3 * pow(r, 3) *
                           (69696. * pow(t, 6) + 42816. * pow(t, 5) * T +
                            48. * pow(t, 4) * (2033. + 612. * pow(T, 2)) +
                            32. * pow(t, 3) * T * (1179. + 638. * pow(T, 2)) +
                            t * T *
                                (2907. + 5192. * pow(T, 2) +
                                 8320. * pow(T, 4)) +
                            pow(t, 2) * (15525. + 15144. * pow(T, 2) +
                                         13696. * pow(T, 4)) +
                            192. * (1. + 2. * pow(T, 2) + 6. * pow(T, 4) +
                                    20. * pow(T, 6))) +
                       7. * beta *
                           (7. * c2 * pow(r, 2) *
                                (9408. * pow(t, 5) + 5568. * pow(t, 4) * T +
                                 96. * pow(t, 3) * (91. + 38. * pow(T, 2)) +
                                 16. * pow(t, 2) * T *
                                     (183. + 148. * pow(T, 2)) +
                                 5. * T *
                                     (15. + 40. * pow(T, 2) +
                                      128. * pow(T, 4)) +
                                 t * (693. + 936. * pow(T, 2) +
                                      1408. * pow(T, 4))) +
                            8. * beta *
                                (6. * c1 * r *
                                     (1096. * pow(t, 4) + 616. * pow(t, 3) * T +
                                      t * T * (161. + 216. * pow(T, 2)) +
                                      pow(t, 2) * (607. + 376. * pow(T, 2)) +
                                      16. * (1. + 2. * pow(T, 2) +
                                             6. * pow(T, 4))) +
                                 25. * c0 *
                                     (55. * t + 200. * pow(t, 3) + 9. * T +
                                      104. * pow(t, 2) * T +
                                      56. * t * pow(T, 2) + 24. * pow(T, 3)) *
                                     beta))))) +
        315. *
            (10. * c5 * pow(r, 5) * t *
                 (315. + 13440. * pow(t, 2) + 96768. * pow(t, 4) +
                  147456. * pow(t, 6) + 32768. * pow(t, 8)) +
             beta * (9. * c4 * pow(r, 4) *
                         (35. + 4480. * pow(t, 2) + 53760. * pow(t, 4) +
                          114688. * pow(t, 6) + 32768. * pow(t, 8)) +
                     32. * beta *
                         (8. * c3 * pow(r, 3) * t *
                              (35. + 840. * pow(t, 2) + 2688. * pow(t, 4) +
                               1024. * pow(t, 6)) +
                          beta * (7. * c2 * pow(r, 2) *
                                      (5. + 360. * pow(t, 2) +
                                       1920. * pow(t, 4) + 1024. * pow(t, 6)) +
                                  8. * beta *
                                      (6. * c1 * r * t *
                                           (15. + 160. * pow(t, 2) +
                                            128. * pow(t, 4)) +
                                       5. * c0 *
                                           (3. + 96. * pow(t, 2) +
                                            128. * pow(t, 4)) *
                                           beta))))) *
            log(-t + sqrt(-0.25 + pow(t - T, 2)) + T))) /
      (3.3030144e8 * sqrt(3.) * pow(beta, 9));

  return tmp.real();
}

//---------------------------------------------------------------------------
void EW::exactRhsTwilight(float_sw4 a_t, vector<Sarray>& a_F) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *f_ptr, om, ph, cv, h, zmin, omm, phm, amprho, ampmu, ampla;

  int g;

  for (g = 0; g < mNumberOfCartesianGrids; g++) {
    f_ptr = a_F[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    h = mGridSize[g];  // how do we define the grid size for the curvilinear
                       // grid?
    zmin = m_zmin[g];
    if (m_twilight_forcing) {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
      omm = m_twilight_forcing->m_momega;
      phm = m_twilight_forcing->m_mphase;
      amprho = m_twilight_forcing->m_amprho;
      ampmu = m_twilight_forcing->m_ampmu;
      ampla = m_twilight_forcing->m_amplambda;
    }
    if (usingSupergrid()) {
      float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
      float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
      float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
      if (m_croutines)
        exactrhsfortsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                          a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla, h,
                          zmin, omstrx, omstry, omstrz);
      else
        exactrhsfortsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, f_ptr,
                       &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
                       &h, &zmin, &omstrx, &omstry, &omstrz);
    } else {
      if (m_croutines)
        exactrhsfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr, a_t,
                        om, cv, ph, omm, phm, amprho, ampmu, ampla, h, zmin);
      else
        exactrhsfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, f_ptr,
                     &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
                     &h, &zmin);
    }
  }
  // if (topographyExists()) {
  for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
    // g = mNumberOfGrids - 1;
    f_ptr = a_F[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    if (m_twilight_forcing) {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
      omm = m_twilight_forcing->m_momega;
      phm = m_twilight_forcing->m_mphase;
      amprho = m_twilight_forcing->m_amprho;
      ampmu = m_twilight_forcing->m_ampmu;
      ampla = m_twilight_forcing->m_amplambda;
    }
    //  subroutine exactaccfort( ifirst, ilast, jfirst, jlast, kfirst,
    // +     klast, utt, t, om, c, ph, h, zmin )
    if (usingSupergrid()) {
      float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
      float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
      float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
      if (m_croutines)
        exactrhsfortsgc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                           a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla,
                           mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr(), omstrx,
                           omstry, omstrz);
      else
        exactrhsfortsgc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                        f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu,
                        &ampla, mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr(),
                        &omstrx, &omstry, &omstrz);
    } else {
      if (m_croutines)
        exactrhsfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                         a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla,
                         mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());
      else
        exactrhsfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, f_ptr,
                      &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu, &ampla,
                      mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());
    }
  }
}

//---------------------------------------------------------------------------
void EW::exactAccTwilight(float_sw4 a_t, vector<Sarray>& a_Uacc) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *uacc_ptr, om, ph, cv, h, zmin;

  int g;

  for (g = 0; g < mNumberOfCartesianGrids; g++) {
    uacc_ptr = a_Uacc[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    h = mGridSize[g];  // how do we define the grid size for the curvilinear
                       // grid?
    zmin = m_zmin[g];
    if (m_twilight_forcing) {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
    }

    //  subroutine exactaccfort( ifirst, ilast, jfirst, jlast, kfirst,
    // +     klast, utt, t, om, c, ph, h, zmin )
    if (m_croutines)
      exactaccfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, uacc_ptr,
                      a_t, om, cv, ph, h, zmin);
    else
      exactaccfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, uacc_ptr,
                   &a_t, &om, &cv, &ph, &h, &zmin);
  }
  //  if (topographyExists()) {
  for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
    // g = mNumberOfGrids - 1;
    uacc_ptr = a_Uacc[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    if (m_twilight_forcing) {
      om = m_twilight_forcing->m_omega;
      ph = m_twilight_forcing->m_phase;
      cv = m_twilight_forcing->m_c;
    }
    //  subroutine exactaccfort( ifirst, ilast, jfirst, jlast, kfirst,
    // +     klast, utt, t, om, c, ph, h, zmin )
    if (m_croutines)
      exactaccfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, uacc_ptr,
                       a_t, om, cv, ph, mX[g].c_ptr(), mY[g].c_ptr(),
                       mZ[g].c_ptr());
    else
      exactaccfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, uacc_ptr,
                    &a_t, &om, &cv, &ph, mX[g].c_ptr(), mY[g].c_ptr(),
                    mZ[g].c_ptr());
  }
}

//---------------------------------------------------------------------------
#ifdef NO_DEVICE_FUNCTION_POINTERS
// Awaiting device function pointer support in HIP
void EW::Force(float_sw4 a_t, vector<Sarray>& a_F,
               vector<GridPointSource*>& point_sources,
               vector<int>& identsources) {
  static bool first = true;
  if (first) {
    first = false;
    std::cout << "WARNING **** NON_FUNCTIONAL CALL TO EW::Force\n";
  }
  for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
}
#else
void EW::Force(float_sw4 a_t, vector<Sarray>& a_F,
               vector<GridPointSource*>& point_sources,
               vector<int>& identsources) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *f_ptr, om, ph, cv, h, zmin, omm, phm, amprho, ampmu, ampla;

  int g;

  if (m_twilight_forcing) {
    if (m_anisotropic) {
      float_sw4 phc[21];  // move these angles to the EW class

      // need to store all the phase angle constants somewhere
      for (int i = 0; i < 21; i++) phc[i] = i * 10 * M_PI / 180;

      for (g = 0; g < mNumberOfCartesianGrids; g++) {
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        h = mGridSize[g];
        zmin = m_zmin[g];

        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;

        if (m_croutines)
          tw_aniso_force_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                            a_t, om, cv, ph, omm, phm, amprho, phc, h, zmin);
        else
          tw_aniso_force(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                         a_t, om, cv, ph, omm, phm, amprho, phc, h, zmin);
      }  // end for all Cartesian grids
      // if (topographyExists()) {
      for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
        // g = mNumberOfGrids - 1;
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;

        if (m_croutines)
          tw_aniso_curvi_force_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                  f_ptr, a_t, om, cv, ph, omm, phm, amprho, phc,
                                  mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());
        else
          tw_aniso_curvi_force(ifirst, ilast, jfirst, jlast, kfirst, klast,
                               f_ptr, a_t, om, cv, ph, omm, phm, amprho, phc,
                               mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());

      }  // end if topographyExists

    } else {  // isotropic twilight forcing
      for (g = 0; g < mNumberOfCartesianGrids; g++) {
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        h = mGridSize[g];  // how do we define the grid size for the curvilinear
                           // grid?
        zmin = m_zmin[g];

        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        ampmu = m_twilight_forcing->m_ampmu;
        ampla = m_twilight_forcing->m_amplambda;
        if (usingSupergrid()) {
          float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
          float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
          float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
          if (m_croutines)
            forcingfortsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                             a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla, h,
                             zmin, omstrx, omstry, omstrz);
          else
            forcingfortsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                          f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                          &ampmu, &ampla, &h, &zmin, &omstrx, &omstry, &omstrz);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingfortsgatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                  f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                  ampmu, ampla, h, zmin, omstrx, omstry,
                                  omstrz);
            else
              forcingfortsgatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                               &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                               &amprho, &ampmu, &ampla, &h, &zmin, &omstrx,
                               &omstry, &omstrz);
          }
        } else {
          if (m_croutines)
            forcingfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                           a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla, h,
                           zmin);
          else
            forcingfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                        f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu,
                        &ampla, &h, &zmin);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingfortatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                                ampla, h, zmin);
            else
              forcingfortatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                             f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                             &ampmu, &ampla, &h, &zmin);
          }
        }
      }  // end for all Cartesian grids

      //      if (topographyExists()) {
      for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
        // g = mNumberOfGrids - 1;
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        ampmu = m_twilight_forcing->m_ampmu;
        ampla = m_twilight_forcing->m_amplambda;
        if (usingSupergrid()) {
          float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
          float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
          float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
          if (m_croutines)
            forcingfortcsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                              f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                              ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                              mZ[g].c_ptr(), omstrx, omstry, omstrz);
          else
            forcingfortcsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                           f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                           &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                           mZ[g].c_ptr(), &omstrx, &omstry, &omstrz);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingfortsgattc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                   f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                   ampmu, ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                                   mZ[g].c_ptr(), omstrx, omstry, omstrz);
            else
              forcingfortsgattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                                &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                                &amprho, &ampmu, &ampla, mX[g].c_ptr(),
                                mY[g].c_ptr(), mZ[g].c_ptr(), &omstrx, &omstry,
                                &omstrz);
          }
        } else {
          if (m_croutines)
            forcingfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                            a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla,
                            mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());
          else
            forcingfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                         f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                         &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                         mZ[g].c_ptr());
          if (m_use_attenuation) {
            if (m_croutines)
              forcingfortattc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                 f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                 ampmu, ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                                 mZ[g].c_ptr());
            else
              forcingfortattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                              f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                              &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                              mZ[g].c_ptr());
          }
        }
      }
    }  // end isotropic case

  }  // end twilight

  else if (m_rayleigh_wave_test) {
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
  } else if (m_energy_test) {
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
  } else {
    //
    // WARNING:: THERE IS AN ASSUMPTION HERE THAT FORCE WILL BE CALLED BEFORE
    // FORCE_TT
    //
    static bool firstcall = true;

    GridPointSource** GPSL;
    int* idnts_local;
    if (firstcall) {
      // WARNING :: FORCE_TT cannot be called until this section has been called
      SW4_MARK_BEGIN("FORCE::HOST::FIRSTCALL");
      // size_t mfree,mtotal;
      // SW4_CheckDeviceError(cudaMemGetInfo(&mfree,&mtotal));
      // std::cout<<getRank()<<" MFREE PRE-ALLOC"<<mfree/1024/1024.0<<" MB
      // "<<point_sources.size()<<","<<identsourc
      //	es.size()<<"\n";
      GPS = SW4_NEW(Space::Managed, GridPointSource * [point_sources.size()]);
      idnts = SW4_NEW(Space::Managed, int[identsources.size()]);
      GPSL = GPS;
      idnts_local = idnts;

      for (int r = 0; r < identsources.size(); r++) idnts[r] = identsources[r];
      for (int s = 0; s < point_sources.size(); s++) GPS[s] = point_sources[s];

        // if (point_sources.size()>0){
        // SW4_CheckDeviceError(cudaMemPrefetchAsync(GPS,
        // 						point_sources.size()*sizeof(GridPointSource*),
        // 						0,
        // 						0));

        // SW4_CheckDeviceError(cudaMemPrefetchAsync(idnts,
        // 						identsources.size()*sizeof(int),
        // 						0,
        // 						0));
        // }
        // SW4_CheckDeviceError(cudaMemGetInfo(&mfree,&mtotal));
        // std::cout<<getRank()<<" MFREE POST-ALLOC"<<mfree/1024/1024.0<<" MB
        // \n";
#pragma omp parallel for
      for (int r = 0; r < identsources.size() - 1; r++) {
        int index = r * 3;
        int s0 = identsources[r];
        int g = point_sources[s0]->m_grid;
        int i = point_sources[s0]->m_i0;
        int j = point_sources[s0]->m_j0;
        int k = point_sources[s0]->m_k0;

        size_t ind1 = a_F[g].index(1, i, j, k);
        size_t oc = a_F[g].m_offc;
        float_sw4* fptr = a_F[g].c_ptr();
        ForceAddress[index] = fptr + ind1;
        ForceAddress[index + 1] = fptr + ind1 + oc;
        ForceAddress[index + 2] = fptr + ind1 + 2 * oc;
      }

      // if (point_sources.size()>0) std::cerr<<getRank()<<" Calling
      // GPS[r]->initializeTimeFunction() "<<point_sources.size()<<" \n";
      // for (int i=0;i<point_sources.size();i++) GPSL[i]->print_vals();
      RAJA::forall<DEFAULT_LOOP1>(
          RAJA::RangeSegment(0, point_sources.size()),
          [=] RAJA_DEVICE(int r) { GPSL[r]->initializeTimeFunction(); });
      // if (point_sources.size()>0)std::cerr<<"Done Calling
      // GPS[r]->initializeTimeFunction()  "<<point_sources.size()<<" \n";
      SW4_MARK_END("FORCE::HOST::FIRSTCALL");

      firstcall = false;
    }
    // #ifdef ENABLE_CUDA
    //     typedef RAJA::cuda_exec<32, true> FORCE_LOOP_ASYNC;
    // #else
    //     using FORCE_LOOP_ASYNC = RAJA::omp_parallel_for_exec;
    // #endif
    GPSL = GPS;
    idnts_local = idnts;
    float_sw4** ForceAddress_copy = ForceAddress;

    // for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero_async();
    vset_to_zero_async(a_F, mNumberOfGrids);
    SW4_MARK_BEGIN("FORCE::DEVICE");

    RAJA::forall<FORCE_LOOP_ASYNC>(
        RAJA::RangeSegment(0, identsources.size() - 1), [=] RAJA_DEVICE(int r) {
          int index = r * 3;
          for (int s = idnts_local[r]; s < idnts_local[r + 1]; s++) {
            float_sw4 fxyz[3];
            GPSL[s]->getFxyz(a_t, fxyz);
#pragma unroll
            for (int i = 0; i < 3; i++)
              *ForceAddress_copy[index + i] += fxyz[i];
          }
        });

    SYNC_STREAM;
    SW4_MARK_END("FORCE::DEVICE");
  }
}
#endif
//---------------------------------------------------------------------------
#ifdef NO_DEVICE_FUNCTION_POINTERS
// Awaiting device function pointer support in HIP
void EW::Force_tt(float_sw4 a_t, vector<Sarray>& a_F,
                  vector<GridPointSource*>& point_sources,
                  vector<int>& identsources) {
  static bool first = true;
  if (first) {
    first = false;
    std::cout << "WARNING **** NON_FUNCTIONAL CALL TO EW::Force_tt\n";
  }
  for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
}
#else
void EW::Force_tt(float_sw4 a_t, vector<Sarray>& a_F,
                  vector<GridPointSource*>& point_sources,
                  vector<int>& identsources) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *f_ptr, om, ph, cv, h, zmin, omm, phm, amprho, ampmu, ampla;
  // std::cout<<"FORCE_TT\n";
  int g;
  // std::cerr<<"And now in force_tt\n";
  if (m_twilight_forcing) {
    if (m_anisotropic) {
      float_sw4 phc[21];  // move these angles to the EW class

      // need to store all the phase angle constants somewhere
      phc[0] = 0;
      for (int i = 0; i < 21; i++) phc[i] = i * 10 * M_PI / 180;

      for (g = 0; g < mNumberOfCartesianGrids; g++) {
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        h = mGridSize[g];
        zmin = m_zmin[g];

        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;

        if (m_croutines)
          tw_aniso_force_tt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                               f_ptr, a_t, om, cv, ph, omm, phm, amprho, phc, h,
                               zmin);
        else
          tw_aniso_force_tt(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                            a_t, om, cv, ph, omm, phm, amprho, phc, h, zmin);
      }  // end for all Cartesian grids

      // if (topographyExists()) {
      for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
        // g = mNumberOfGrids - 1;
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        if (m_croutines)
          tw_aniso_curvi_force_tt_ci(ifirst, ilast, jfirst, jlast, kfirst,
                                     klast, f_ptr, a_t, om, cv, ph, omm, phm,
                                     amprho, phc, mX[g].c_ptr(), mY[g].c_ptr(),
                                     mZ[g].c_ptr());
        else
          tw_aniso_curvi_force_tt(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                  f_ptr, a_t, om, cv, ph, omm, phm, amprho, phc,
                                  mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());

      }  // end if topographyExists

    } else {  // isotropic twilight forcing
      for (g = 0; g < mNumberOfCartesianGrids; g++) {
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        h = mGridSize[g];  // how do we define the grid size for the curvilinear
                           // grid?
        zmin = m_zmin[g];

        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        ampmu = m_twilight_forcing->m_ampmu;
        ampla = m_twilight_forcing->m_amplambda;
        if (usingSupergrid()) {
          float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
          float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
          float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
          if (m_croutines)
            forcingttfortsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                               f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                               ampla, h, zmin, omstrx, omstry, omstrz);
          else
            forcingttfortsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                            f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                            &ampmu, &ampla, &h, &zmin, &omstrx, &omstry,
                            &omstrz);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingttfortsgatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                    f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                    ampmu, ampla, h, zmin, omstrx, omstry,
                                    omstrz);
            else
              forcingttfortsgatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                                 &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                                 &amprho, &ampmu, &ampla, &h, &zmin, &omstrx,
                                 &omstry, &omstrz);
          }
        } else {
          if (m_croutines)
            forcingttfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                             a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla, h,
                             zmin);
          else
            forcingttfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                          f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                          &ampmu, &ampla, &h, &zmin);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingttattfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                  f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                  ampmu, ampla, h, zmin);
            else
              forcingttattfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                               &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                               &amprho, &ampmu, &ampla, &h, &zmin);
          }
        }
      }
      //     if (topographyExists()) {
      for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
        // g = mNumberOfGrids - 1;
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        ampmu = m_twilight_forcing->m_ampmu;
        ampla = m_twilight_forcing->m_amplambda;
        if (usingSupergrid()) {
          float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
          float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
          float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
          if (m_croutines)
            forcingttfortcsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                                ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                                mZ[g].c_ptr(), omstrx, omstry, omstrz);
          else
            forcingttfortcsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                             f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                             &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                             mZ[g].c_ptr(), &omstrx, &omstry, &omstrz);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingttfortsgattc_ci(
                  ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr, a_t, om,
                  cv, ph, omm, phm, amprho, ampmu, ampla, mX[g].c_ptr(),
                  mY[g].c_ptr(), mZ[g].c_ptr(), omstrx, omstry, omstrz);
            else
              forcingttfortsgattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                                  &klast, f_ptr, &a_t, &om, &cv, &ph, &omm,
                                  &phm, &amprho, &ampmu, &ampla, mX[g].c_ptr(),
                                  mY[g].c_ptr(), mZ[g].c_ptr(), &omstrx,
                                  &omstry, &omstrz);
          }
        } else {
          if (m_croutines)
            forcingttfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                              f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                              ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                              mZ[g].c_ptr());
          else
            forcingttfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                           f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                           &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                           mZ[g].c_ptr());
          if (m_use_attenuation) {
            if (m_croutines)
              forcingttattfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                   f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                   ampmu, ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                                   mZ[g].c_ptr());
            else
              forcingttattfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                                &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                                &amprho, &ampmu, &ampla, mX[g].c_ptr(),
                                mY[g].c_ptr(), mZ[g].c_ptr());
          }
        }
      }
    }  // end isotropic

  }  // end twilight

  else if (m_rayleigh_wave_test) {
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
  } else if (m_energy_test) {
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
  } else {
    // Default: m_point_source_test, m_lamb_test or full seismic case
    // for( int g =0 ; g < mNumberOfGrids ; g++ )
    // 	a_F[g].set_to_zero();

    GridPointSource** GPSL = GPS;
    int* idnts_local = idnts;

    float_sw4** ForceAddress_copy = ForceAddress;
    // #ifdef ENABLE_CUDA
    //     typedef RAJA::cuda_exec<1024, true> FORCETT_LOOP_ASYNC;
    // #else
    //     using FORCETT_LOOP_ASYNC = RAJA::omp_parallel_for_exec;
    // #endif

    // for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero_async();
    vset_to_zero_async(a_F, mNumberOfGrids);
    SW4_MARK_BEGIN("FORCE_TT::DEVICE");

    RAJA::forall<FORCETT_LOOP_ASYNC>(
        RAJA::RangeSegment(0, identsources.size() - 1), [=] RAJA_DEVICE(int r) {
          int index = r * 3;
          for (int s = idnts_local[r]; s < idnts_local[r + 1]; s++) {
            float_sw4 fxyz[3];
            GPSL[s]->getFxyztt(a_t, fxyz);
#pragma unroll
            for (int i = 0; i < 3; i++)
              *ForceAddress_copy[index + i] += fxyz[i];
          }
        });

    SYNC_STREAM;
    SW4_MARK_END("FORCE_TT::DEVICE");
  }
}
#endif

//---------------------------------------------------------------------------
// perhaps a better name would be evalLu ??
void EW::evalRHS(vector<Sarray>& a_U, vector<Sarray>& a_Mu,
                 vector<Sarray>& a_Lambda, vector<Sarray>& a_Uacc,
                 vector<Sarray*>& a_AlphaVE, std::ostream* norm_trace_file) {
  SW4_MARK_FUNCTION;
#ifdef PEEKS_GALORE
  SW4_PEEK;
  SYNC_DEVICE;
#endif
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *uacc_ptr, *u_ptr, *mu_ptr, *la_ptr, h;

  int* onesided_ptr;

  int g, nz;
  vset_to_zero_async(a_Uacc, mNumberOfGrids);
  for (g = 0; g < mNumberOfCartesianGrids; g++) {
    // a_Uacc[g].prefetch();
    // a_Uacc[g].set_to_zero_async();
    uacc_ptr = a_Uacc[g].c_ptr();
    u_ptr = a_U[g].c_ptr();
    mu_ptr = a_Mu[g].c_ptr();
    la_ptr = a_Lambda[g].c_ptr();

    a_U[g].prefetch();
    a_Mu[g].prefetch();
    a_Lambda[g].prefetch();
#ifdef ENABLE_CUDA
    prefetch_to_device(m_sg_str_x[g]);
    prefetch_to_device(m_sg_str_y[g]);
    prefetch_to_device(m_sg_str_z[g]);
    prefetch_to_device(m_sbop);
#endif
    //    rho_ptr = mRho[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    h = mGridSize[g];  // how do we define the grid size for the curvilinear
                       // grid?
    nz = m_global_nz[g];
    onesided_ptr = m_onesided[g];
    char op = '=';  // Assign Uacc := L(u)
    if (m_croutines) {
      if (usingSupergrid())
        rhs4th3fortsgstr_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
                            onesided_ptr, m_acof, m_bope, m_ghcof, uacc_ptr,
                            u_ptr, mu_ptr, la_ptr, h, m_sg_str_x[g],
                            m_sg_str_y[g], m_sg_str_z[g], op);
      else
        rhs4th3fort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
                       onesided_ptr, m_acof, m_bope, m_ghcof, uacc_ptr, u_ptr,
                       mu_ptr, la_ptr, h, op);
    } else {
      if (usingSupergrid())
        rhs4th3fortsgstr(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
                         onesided_ptr, m_acof, m_bope, m_ghcof, uacc_ptr, u_ptr,
                         mu_ptr, la_ptr, &h, m_sg_str_x[g], m_sg_str_y[g],
                         m_sg_str_z[g], &op);
      else
        rhs4th3fort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
                    onesided_ptr, m_acof, m_bope, m_ghcof, uacc_ptr, u_ptr,
                    mu_ptr, la_ptr, &h, &op);
    }
#ifdef SW4_NORM_TRACE
    if (norm_trace_file != nullptr)
      *norm_trace_file << " evalRHS_1 " << g << " " << a_Uacc[g].norm() << "\n";
#endif
      //    size_t nn=a_Uacc[g].count_nans();
      //    if( nn > 0 )
      //       cout << "First application of LU " << nn << " nans" << endl;
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif

    if (m_use_attenuation && m_number_mechanisms > 0) {
      op = '-';  // Subtract Uacc := Uacc - L_a(alpha)
      for (int a = 0; a < m_number_mechanisms; a++) {
        float_sw4* alpha_ptr = a_AlphaVE[g][a].c_ptr();
        float_sw4* mua_ptr = mMuVE[g][a].c_ptr();
        float_sw4* lambdaa_ptr = mLambdaVE[g][a].c_ptr();
        //	  nn = a_AlphaVE[g][a].count_nans();
        //	  if( nn > 0 )
        //	     cout << "Alpha before LU " << nn << " nans" << endl;

        if (m_croutines) {
          if (usingSupergrid())
            rhs4th3fortsgstr_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
                                onesided_ptr, m_acof_no_gp, m_bope,
                                m_ghcof_no_gp, uacc_ptr, alpha_ptr, mua_ptr,
                                lambdaa_ptr, h, m_sg_str_x[g], m_sg_str_y[g],
                                m_sg_str_z[g], op);
          else
            rhs4th3fort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
                           onesided_ptr, m_acof_no_gp, m_bope, m_ghcof_no_gp,
                           uacc_ptr, alpha_ptr, mua_ptr, lambdaa_ptr, h, op);
        } else {
          if (usingSupergrid())
            rhs4th3fortsgstr(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                             &nz, onesided_ptr, m_acof_no_gp, m_bope,
                             m_ghcof_no_gp, uacc_ptr, alpha_ptr, mua_ptr,
                             lambdaa_ptr, &h, m_sg_str_x[g], m_sg_str_y[g],
                             m_sg_str_z[g], &op);
          else
            rhs4th3fort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz,
                        onesided_ptr, m_acof_no_gp, m_bope, m_ghcof_no_gp,
                        uacc_ptr, alpha_ptr, mua_ptr, lambdaa_ptr, &h, &op);
        }
      }
      //    nn=a_Uacc[g].count_nans();
      //    if( nn > 0 )
      //       cout << "Second application of LU " << nn << " nans" << endl;
    }
#ifdef SW4_NORM_TRACE
    if (norm_trace_file != nullptr)
      *norm_trace_file << " evalRHS_2 " << g << " " << a_Uacc[g].norm() << "\n";
#endif
  }

  //  if (topographyExists()) {
#ifdef PEEKS_GALORE
  SW4_PEEK;
  SYNC_DEVICE;
#endif
  for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
    // g = mNumberOfGrids - 1;
    // a_Uacc[g].set_to_zero_async(); // Being done above
    uacc_ptr = a_Uacc[g].c_ptr();
    u_ptr = a_U[g].c_ptr();
    mu_ptr = a_Mu[g].c_ptr();
    la_ptr = a_Lambda[g].c_ptr();
    float_sw4* met_ptr = mMetric[g].c_ptr();
    float_sw4* jac_ptr = mJ[g].c_ptr();
    if (global_variables.firstCycle) {
      mMetric[g].forceprefetch();
      mJ[g].forceprefetch();
    }
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    onesided_ptr = m_onesided[g];
    int nkg = m_global_nz[g];
    char op = '=';  // assign Uacc := L_u(u)
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
    if (m_croutines) {
#if defined(SW4_EXPT_1)
      curvilinear4sgX_ci<3>(
          ifirst, ilast, jfirst, jlast, kfirst, klast, u_ptr, mu_ptr, la_ptr,
          met_ptr, jac_ptr, uacc_ptr, onesided_ptr, m_acof, m_bope, m_ghcof,
          m_acof_no_gp, m_ghcof_no_gp, m_sg_str_x[g], m_sg_str_y[g], nkg, op);
#elif defined(SW4_EXPT_3)
      // cudaMemcpyToSymbol(tex_acof, m_acof, 384*sizeof(double));
      curvilinear4sgX3_ci<0>(
          ifirst, ilast, jfirst, jlast, kfirst, klast, u_ptr, mu_ptr, la_ptr,
          met_ptr, jac_ptr, uacc_ptr, onesided_ptr, m_acof, m_bope, m_ghcof,
          m_acof_no_gp, m_ghcof_no_gp, m_sg_str_x[g], m_sg_str_y[g], nkg, op);
#else
      curvilinear4sg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, u_ptr,
                        mu_ptr, la_ptr, met_ptr, jac_ptr, uacc_ptr,
                        onesided_ptr, m_acof, m_bope, m_ghcof, m_acof_no_gp,
                        m_ghcof_no_gp, m_sg_str_x[g], m_sg_str_y[g], nkg, op);
#endif

    } else {
      if (usingSupergrid())
        curvilinear4sg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, u_ptr,
                       mu_ptr, la_ptr, met_ptr, jac_ptr, uacc_ptr, onesided_ptr,
                       m_acof, m_bope, m_ghcof, m_sg_str_x[g], m_sg_str_y[g],
                       &op);
      else
        curvilinear4(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, u_ptr,
                     mu_ptr, la_ptr, met_ptr, jac_ptr, uacc_ptr, onesided_ptr,
                     m_acof, m_bope, m_ghcof, &op);
    }
#ifdef SW4_NORM_TRACE
    if (norm_trace_file != nullptr) {
      *norm_trace_file << " evalRHS_3 " << g << " " << a_Uacc[g].norm() << "\n";
      *norm_trace_file << "   evalRHS_3 U[" << g << "]= " << a_U[g].norm()
                       << "\n";
      *norm_trace_file << "   evalRHS_3 Mu[" << g << "]= " << a_Mu[g].norm()
                       << "\n";
      *norm_trace_file << "   evalRHS_3 Lambda[" << g
                       << "]= " << a_Lambda[g].norm() << "\n";
      *norm_trace_file << "   evalRHS_3 Metric[" << g
                       << "]= " << mMetric[g].norm() << "\n";
      *norm_trace_file << "   evalRHS_3 Jaco[" << g << "]= " << mJ[g].norm()
                       << "\n";
    }
#endif
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
    if (m_use_attenuation && m_number_mechanisms > 0) {
      op = '-';  // Subtract Uacc := Uacc - L_a(alpha)
      for (int a = 0; a < m_number_mechanisms; a++) {
        float_sw4* alpha_ptr = a_AlphaVE[g][a].c_ptr();
        float_sw4* mua_ptr = mMuVE[g][a].c_ptr();
        float_sw4* lambdaa_ptr = mLambdaVE[g][a].c_ptr();
        if (m_croutines) {
#ifdef STANDALONE_SETUP
          static bool once = false;
          if ((!once) && (getRank() == 0)) {
            Apc apc("curvkernel.CPROTO");
            autopeel(apc, ifirst, ilast, jfirst, jlast, kfirst, klast,
                     alpha_ptr, mua_ptr, lambdaa_ptr, met_ptr, jac_ptr,
                     uacc_ptr, onesided_ptr, m_acof_no_gp, m_bope,
                     m_ghcof_no_gp, m_sg_str_x[g], m_sg_str_y[g], op);
            once = true;
          }
#endif

#if defined(SW4_EXPT_1)
          curvilinear4sgX_ci<3>(
              ifirst, ilast, jfirst, jlast, kfirst, klast, alpha_ptr, mua_ptr,
              lambdaa_ptr, met_ptr, jac_ptr, uacc_ptr, onesided_ptr,
              m_acof_no_gp, m_bope, m_ghcof_no_gp, m_acof_no_gp, m_ghcof_no_gp,
              m_sg_str_x[g], m_sg_str_y[g], nkg, op);
#elif defined(SW4_EXPT_3)
          curvilinear4sgX3_ci<1>(
              ifirst, ilast, jfirst, jlast, kfirst, klast, alpha_ptr, mua_ptr,
              lambdaa_ptr, met_ptr, jac_ptr, uacc_ptr, onesided_ptr,
              m_acof_no_gp, m_bope, m_ghcof_no_gp, m_acof_no_gp, m_ghcof_no_gp,
              m_sg_str_x[g], m_sg_str_y[g], nkg, op);
#else
          // cudaMemcpyToSymbol(tex_acof, m_acof_no_gp, 384*sizeof(double));
          curvilinear4sg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                            alpha_ptr, mua_ptr, lambdaa_ptr, met_ptr, jac_ptr,
                            uacc_ptr, onesided_ptr, m_acof_no_gp, m_bope,
                            m_ghcof_no_gp, m_acof_no_gp, m_ghcof_no_gp,
                            m_sg_str_x[g], m_sg_str_y[g], nkg, op);
#endif
        } else {
          if (usingSupergrid())
            curvilinear4sg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                           alpha_ptr, mua_ptr, lambdaa_ptr, met_ptr, jac_ptr,
                           uacc_ptr, onesided_ptr, m_acof_no_gp, m_bope,
                           m_ghcof_no_gp, m_sg_str_x[g], m_sg_str_y[g], &op);
          else
            curvilinear4(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                         alpha_ptr, mua_ptr, lambdaa_ptr, met_ptr, jac_ptr,
                         uacc_ptr, onesided_ptr, m_acof_no_gp, m_bope,
                         m_ghcof_no_gp, &op);
        }
      }
#ifdef SW4_NORM_TRACE
      if (norm_trace_file != nullptr)
        *norm_trace_file << " evalRHS_4 " << g << " " << a_Uacc[g].norm()
                         << "\n";
#endif
    }
    // SYNC_STREAM;
#ifdef PEEKS_GALORE
    SW4_PEEK;
    SYNC_DEVICE;
#endif
  }
  SYNC_STREAM;  // REQUIRED IF THERE IS NO SYNC IN curvilinear4sg_ci
}

//-----------------------------------------------------------------------
void EW::evalRHSanisotropic(vector<Sarray>& a_U, vector<Sarray>& a_C,
                            vector<Sarray>& a_Uacc) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *uacc_ptr, *u_ptr, *c_ptr, h;

  int* onesided_ptr;

  int g, nz;

  for (g = 0; g < mNumberOfCartesianGrids; g++) {
    uacc_ptr = a_Uacc[g].c_ptr();
    u_ptr = a_U[g].c_ptr();
    c_ptr = a_C[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    h = mGridSize[g];
    nz = m_global_nz[g];
    onesided_ptr = m_onesided[g];
    if (m_croutines)
      innerloopanisgstrvc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
                             u_ptr, uacc_ptr, c_ptr, onesided_ptr, m_acof,
                             m_bope, m_ghcof, h, m_sg_str_x[g], m_sg_str_y[g],
                             m_sg_str_z[g]);
    else
      innerloopanisgstrvc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                          &nz, u_ptr, uacc_ptr, c_ptr, onesided_ptr, m_acof,
                          m_bope, m_ghcof, &h, m_sg_str_x[g], m_sg_str_y[g],
                          m_sg_str_z[g]);
  }
  //  if (topographyExists()) {
  for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
    // g = mNumberOfGrids - 1;
    uacc_ptr = a_Uacc[g].c_ptr();
    u_ptr = a_U[g].c_ptr();
    c_ptr = mCcurv.c_ptr();
    float_sw4* met_ptr = mMetric[g].c_ptr();
    float_sw4* jac_ptr = mJ[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    nz = m_global_nz[g];
    if (m_croutines)
      ilanisocurv_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, nz, u_ptr,
                     c_ptr, jac_ptr, uacc_ptr, m_onesided[g], m_acof, m_bope,
                     m_ghcof, m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g]);
    else
      ilanisocurv(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &nz, u_ptr,
                  c_ptr, jac_ptr, uacc_ptr, m_onesided[g], m_acof, m_bope,
                  m_ghcof, m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g]);
  }
}

//---------------------------------------------------------------------------
void EW::evalPredictor(vector<Sarray>& a_Up, vector<Sarray>& a_U,
                       vector<Sarray>& a_Um, vector<Sarray>& a_Rho,
                       vector<Sarray>& a_Lu, vector<Sarray>& a_F) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *up_ptr, *u_ptr, *um_ptr, *lu_ptr, *fo_ptr, *rho_ptr, dt2;

  // int *onesided_ptr;

  int g;

  for (g = 0; g < mNumberOfGrids; g++) {
    a_Up[g].prefetch();
    a_U[g].prefetch();
    a_Um[g].prefetch();
    a_Lu[g].prefetch();
    a_F[g].prefetch();
    a_Rho[g].prefetch();
    up_ptr = a_Up[g].c_ptr();
    u_ptr = a_U[g].c_ptr();
    um_ptr = a_Um[g].c_ptr();
    lu_ptr = a_Lu[g].c_ptr();
    fo_ptr = a_F[g].c_ptr();
    rho_ptr = a_Rho[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    dt2 = mDt * mDt;
    if (m_croutines)
      predfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, up_ptr, u_ptr,
                  um_ptr, lu_ptr, fo_ptr, rho_ptr, dt2);
    else
      predfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, up_ptr, u_ptr,
               um_ptr, lu_ptr, fo_ptr, rho_ptr, &dt2);
  }
  SYNC_STREAM;
}

//---------------------------------------------------------------------------
void EW::evalCorrector(vector<Sarray>& a_Up, vector<Sarray>& a_Rho,
                       vector<Sarray>& a_Lu, vector<Sarray>& a_F) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *up_ptr, *lu_ptr, *fo_ptr, *rho_ptr, dt4;

  int g;

  for (g = 0; g < mNumberOfGrids; g++) {
    up_ptr = a_Up[g].c_ptr();
    lu_ptr = a_Lu[g].c_ptr();
    fo_ptr = a_F[g].c_ptr();
    rho_ptr = a_Rho[g].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    dt4 = mDt * mDt * mDt * mDt;

    //  subroutine corrfort(ifirst, ilast, jfirst, jlast, kfirst, klast,
    // +     up, lu, fo, rho, dt4 )
    if (m_croutines)
      corrfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, up_ptr, lu_ptr,
                  fo_ptr, rho_ptr, dt4);
    else
      corrfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, up_ptr,
               lu_ptr, fo_ptr, rho_ptr, &dt4);
  }
  SYNC_STREAM;
}

//---------------------------------------------------------------------------
void EW::evalDpDmInTime(vector<Sarray>& a_Up, vector<Sarray>& a_U,
                        vector<Sarray>& a_Um, vector<Sarray>& a_Uacc) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *up_ptr, *u_ptr, *um_ptr, *uacc_ptr, dt2i;

  int g;

  //  for(g=0 ; g<mNumberOfCartesianGrids; g++ )
  for (g = 0; g < mNumberOfGrids; g++) {
    up_ptr = a_Up[g].c_ptr();
    u_ptr = a_U[g].c_ptr();
    um_ptr = a_Um[g].c_ptr();
    uacc_ptr = a_Uacc[g].c_ptr();

    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    dt2i = 1. / (mDt * mDt);
    if (global_variables.firstCycle) {
      a_Up[g].forceprefetch();
      a_U[g].forceprefetch();
      a_Um[g].forceprefetch();
      a_Uacc[g].forceprefetch();
    }
    //  subroutine dpdmtfort(ifirst, ilast, jfirst, jlast, kfirst, klast,
    // +     up, u, um, u2, dt2i)
    if (m_croutines)
      dpdmtfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, up_ptr, u_ptr,
                   um_ptr, uacc_ptr, dt2i, getRank());
    else
      dpdmtfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, up_ptr,
                u_ptr, um_ptr, uacc_ptr, &dt2i);
  }
}

//-----------------------------------------------------------------------
void EW::updateMemVarPred(vector<Sarray*>& a_AlphaVEp,
                          vector<Sarray*>& a_AlphaVEm, vector<Sarray>& a_U,
                          float_sw4 a_t) {
  SW4_MARK_FUNCTION;
  int domain = 0;

  for (int g = 0; g < mNumberOfGrids; g++) {
    float_sw4* u_ptr = a_U[g].c_ptr();
    int ifirst = m_iStart[g];
    int ilast = m_iEnd[g];
    int jfirst = m_jStart[g];
    int jlast = m_jEnd[g];
    int kfirst = m_kStart[g];
    int klast = m_kEnd[g];
    for (int a = 0; a < m_number_mechanisms; a++) {
      float_sw4* alp_ptr = a_AlphaVEp[g][a].c_ptr();
      float_sw4* alm_ptr = a_AlphaVEm[g][a].c_ptr();
      if (m_croutines)
        memvar_pred_fort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                            alp_ptr, alm_ptr, u_ptr, mOmegaVE[a], mDt, domain);
      else
        memvar_pred_fort(ifirst, ilast, jfirst, jlast, kfirst, klast, alp_ptr,
                         alm_ptr, u_ptr, mOmegaVE[a], mDt, domain);
    }
    if (m_twilight_forcing) {
      SYNC_STREAM;
      float_sw4* alp_ptr = a_AlphaVEp[g][0].c_ptr();
      float_sw4 om = m_twilight_forcing->m_omega;
      float_sw4 ph = m_twilight_forcing->m_phase;
      float_sw4 cv = m_twilight_forcing->m_c;
      // if (topographyExists() && g == mNumberOfGrids - 1)
      if (topographyExists() && g >= mNumberOfCartesianGrids)
        addMemVarPredCurvilinear(mX[g], mY[g], mZ[g], a_t, a_AlphaVEp[g][0],
                                 mOmegaVE[0], mDt, om, ph, cv);
      else {
        // this routine comes from WPP
        //  It  works with SG stretching because no spatial derivatives occur in
        //  the forcing
        addMemVarPredCart(m_zmin[g], mGridSize[g], a_t, a_AlphaVEp[g][0],
                          mOmegaVE[0], mDt, om, ph, cv);
      }
    }
  }
  SYNC_STREAM;
}

//-----------------------------------------------------------------------
void EW::updateMemVarCorr(vector<Sarray*>& a_AlphaVEp,
                          vector<Sarray*>& a_AlphaVEm, vector<Sarray>& a_Up,
                          vector<Sarray>& a_U, vector<Sarray>& a_Um,
                          double a_t) {
  SW4_MARK_FUNCTION;
  int domain = 0;

  for (int g = 0; g < mNumberOfGrids; g++) {
    float_sw4* up_ptr = a_Up[g].c_ptr();
    float_sw4* u_ptr = a_U[g].c_ptr();
    float_sw4* um_ptr = a_Um[g].c_ptr();

    int ifirst = m_iStart[g];
    int ilast = m_iEnd[g];
    int jfirst = m_jStart[g];
    int jlast = m_jEnd[g];
    int kfirst = m_kStart[g];
    int klast = m_kEnd[g];
    for (int a = 0; a < m_number_mechanisms; a++) {
      float_sw4* alp_ptr = a_AlphaVEp[g][a].c_ptr();
      float_sw4* alm_ptr = a_AlphaVEm[g][a].c_ptr();
      if (m_croutines)
        memvar_corr_fort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                            alp_ptr, alm_ptr, up_ptr, u_ptr, um_ptr,
                            mOmegaVE[a], mDt, domain);
      else
        memvar_corr_fort(ifirst, ilast, jfirst, jlast, kfirst, klast, alp_ptr,
                         alm_ptr, up_ptr, u_ptr, um_ptr, mOmegaVE[a], mDt,
                         domain);
    }
    if (m_twilight_forcing) {
      SYNC_STREAM;
      double* alp_ptr = a_AlphaVEp[g][0].c_ptr();
      double om = m_twilight_forcing->m_omega;
      double ph = m_twilight_forcing->m_phase;
      double cv = m_twilight_forcing->m_c;
      if (topographyExists() && g >= mNumberOfCartesianGrids) {
        // if (topographyExists() && g == mNumberOfGrids - 1) {
        //            addMemVarCorrCurvilinear( mX, mY, mZ, a_t,
        //            a_AlphaVEp[g][0], mOmegaVE[0], mDt, om, ph, cv);
        addMemVarCorr2Curvilinear(mX[g], mY[g], mZ[g], a_t, a_AlphaVEp[g][0],
                                  mOmegaVE[0], mDt, om, ph, cv);
      } else {
        //  It  works with SG stretching because no spatial derivatives occur in
        //  the forcing
        //            addMemVarCorrCart( m_zmin[g], mGridSize[g], a_t,
        //            a_AlphaVEp[g][0], mOmegaVE[0], mDt, om, ph, cv);
        // NEW June 14, 2017
        addMemVarCorr2Cart(m_zmin[g], mGridSize[g], a_t, a_AlphaVEp[g][0],
                           mOmegaVE[0], mDt, om, ph, cv);
      }
    }
  }
  SYNC_STREAM;
}

//-----------------------------------------------------------------------
void EW::updateMemVarCorrNearInterface(Sarray& a_AlphaVEp, Sarray& a_AlphaVEm,
                                       Sarray& a_Up, Sarray& a_U, Sarray& a_Um,
                                       double a_t, int a_mech, int a_grid) {
  SW4_MARK_FUNCTION;
  // NOTE: this routine updates a_AlphaVEp for mechanism a=a_mech in grid
  // g=a_grid, for all points defined in a_AlphaVEp
  int domain = 0;

  double* up_ptr = a_Up.c_ptr();
  double* u_ptr = a_U.c_ptr();
  double* um_ptr = a_Um.c_ptr();
  // use sizes from a_AlphaVEp for the loop in memvar_corr_fort
  int ifirst = a_AlphaVEp.m_ib;
  int ilast = a_AlphaVEp.m_ie;
  int jfirst = a_AlphaVEp.m_jb;
  int jlast = a_AlphaVEp.m_je;
  int kfirst = a_AlphaVEp.m_kb;
  int klast = a_AlphaVEp.m_ke;
  // assume a_Up, a_U, a_Um and a_AlphaVEm are declared with the same sizes
  int d1b = a_Up.m_ib;
  int d1e = a_Up.m_ie;
  int d2b = a_Up.m_jb;
  int d2e = a_Up.m_je;
  int d3b = a_Up.m_kb;
  int d3e = a_Up.m_ke;

  double* alp_ptr = a_AlphaVEp.c_ptr();
  double* alm_ptr = a_AlphaVEm.c_ptr();
  if (m_croutines)
    memvar_corr_fort_wind_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                             alp_ptr, d1b, d1e, d2b, d2e, d3b, d3e, alm_ptr,
                             up_ptr, u_ptr, um_ptr, mOmegaVE[a_mech], mDt,
                             domain);
  else
    memvar_corr_fort_wind(ifirst, ilast, jfirst, jlast, kfirst, klast, alp_ptr,
                          d1b, d1e, d2b, d2e, d3b, d3e, alm_ptr, up_ptr, u_ptr,
                          um_ptr, mOmegaVE[a_mech], mDt, domain);

  if (m_twilight_forcing) {
    // only 1 mechaism is implemented
    double* alp_ptr = a_AlphaVEp.c_ptr();
    double om = m_twilight_forcing->m_omega;
    double ph = m_twilight_forcing->m_phase;
    double cv = m_twilight_forcing->m_c;
    if (topographyExists() && a_grid >= mNumberOfCartesianGrids) {
      // if (topographyExists() && a_grid == mNumberOfGrids - 1) {
      addMemVarCorr2Curvilinear(mX[a_grid], mY[a_grid], mZ[a_grid], a_t,
                                a_AlphaVEp, mOmegaVE[0], mDt, om, ph, cv);
    } else {
      //  It  works with SG stretching because no spatial derivatives occur in
      //  the forcing
      // NEW June 14, 2017
      // loops over all elements in a_AlphaVEp
      addMemVarCorr2Cart(m_zmin[a_grid], mGridSize[a_grid], a_t, a_AlphaVEp,
                         mOmegaVE[0], mDt, om, ph, cv);
    }
  }
}

//-----------------------------------------------------------------------
void EW::evalDpDmInTimeAtt(vector<Sarray*>& a_AlphaVEp,
                           vector<Sarray*>& a_AlphaVE,
                           vector<Sarray*>& a_AlphaVEm)
// store AlphaVEacc in AlphaVEm
{
  SW4_MARK_FUNCTION;
  float_sw4 dt2i = 1 / (mDt * mDt);
  for (int g = 0; g < mNumberOfGrids; g++) {
    int ifirst = m_iStart[g];
    int ilast = m_iEnd[g];
    int jfirst = m_jStart[g];
    int jlast = m_jEnd[g];
    int kfirst = m_kStart[g];
    int klast = m_kEnd[g];
    for (int a = 0; a < m_number_mechanisms; a++) {
      float_sw4* alphap_ptr = a_AlphaVEp[g][a].c_ptr();
      float_sw4* alpha_ptr = a_AlphaVE[g][a].c_ptr();
      float_sw4* alpham_ptr = a_AlphaVEm[g][a].c_ptr();
      a_AlphaVEp[g][a].prefetch();
      a_AlphaVE[g][a].prefetch();
      a_AlphaVEm[g][a].prefetch();
      if (m_croutines)
        dpdmtfortatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, alphap_ptr,
                        alpha_ptr, alpham_ptr, dt2i);
      else
        dpdmtfortatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                     alphap_ptr, alpha_ptr, alpham_ptr, &dt2i);
    }
  }
  SYNC_STREAM;
}

//-----------------------------------------------------------------------
// side_plane returns the index of the ghost points along side =0,1,2,3,4,5
// (low-i, high-i, low-j, high-j, low-k, high-k)
void EW::side_plane(int g, int side, int wind[6], int nGhost) {
  wind[0] = m_iStart[g];
  wind[1] = m_iEnd[g];
  wind[2] = m_jStart[g];
  wind[3] = m_jEnd[g];
  wind[4] = m_kStart[g];
  wind[5] = m_kEnd[g];
  // wind[0] = m_ib;
  // wind[1] = m_ie;
  // wind[2] = m_jb;
  // wind[3] = m_je;
  // wind[4] = m_kb;
  // wind[5] = m_ke;
  if (side == 0)
    wind[1] = wind[0] + (nGhost - 1);
  else if (side == 1)
    wind[0] = wind[1] - (nGhost - 1);
  else if (side == 2)
    wind[3] = wind[2] + (nGhost - 1);
  else if (side == 3)
    wind[2] = wind[3] - (nGhost - 1);
  else if (side == 4)
    wind[5] = wind[4] + (nGhost - 1);
  else
    wind[4] = wind[5] - (nGhost - 1);
}

//-----------------------------------------------------------------------
void EW::update_images(int currentTimeStep, float_sw4 time,
                       vector<Sarray>& a_Up, vector<Sarray>& a_U,
                       vector<Sarray>& a_Um, vector<Sarray>& a_Rho,
                       vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                       vector<Source*>& a_sources, int dminus, int event) {
  SW4_MARK_FUNCTION;
  //   double maxerr;
  for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex) {
    Image* img = mImageFiles[fIndex];

    if (img->mMode == Image::HMAXDUDT) {
      if (dminus)
        img->update_maxes_hVelMax(a_Up, a_U, mDt);
      else
        img->update_maxes_hVelMax(a_Up, a_Um, 2 * mDt);
    }
    if (img->mMode == Image::HMAX) img->update_maxes_hMax(a_Up);

    if (img->mMode == Image::VMAXDUDT) {
      if (dminus)
        img->update_maxes_vVelMax(a_Up, a_U, mDt);
      else
        img->update_maxes_vVelMax(a_Up, a_Um, 2 * mDt);
    }
    if (img->mMode == Image::VMAX) img->update_maxes_vMax(a_Up);

    // Center time derivatives around t-dt, i.e., (up-um)/(2*dt), except when
    // dminus is set. Use (up-u)/dt assumed centered at t, when dminus is true.
    int td = 0;
    if (!dminus) td = img->is_time_derivative();

    if (img->timeToWrite(time - td * mDt, currentTimeStep - td, mDt)) {
      if (img->mMode == Image::UX)
        img->computeImageQuantity(a_Up, 1);
      else if (img->mMode == Image::UY)
        img->computeImageQuantity(a_Up, 2);
      else if (img->mMode == Image::UZ)
        img->computeImageQuantity(a_Up, 3);
      //         else if(img->mMode == Image::FX )
      //            img->computeImageQuantity(mF, 1);
      //         else if(img->mMode == Image::FY )
      //            img->computeImageQuantity(mF, 2);
      //         else if(img->mMode == Image::FZ )
      //            img->computeImageQuantity(mF, 3);
      else if (img->mMode == Image::RHO)
        img->computeImageQuantity(a_Rho, 1);
      else if (img->mMode == Image::MU)
        img->computeImageQuantity(a_Mu, 1);
      else if (img->mMode == Image::LAMBDA)
        img->computeImageQuantity(a_Lambda, 1);
      else if (img->mMode == Image::QP)
        img->computeImageQuantity(mQp, 1);
      else if (img->mMode == Image::QS)
        img->computeImageQuantity(mQs, 1);
      else if (img->mMode == Image::P)
        img->computeImagePvel(a_Mu, a_Lambda, a_Rho);
      else if (img->mMode == Image::S)
        img->computeImageSvel(a_Mu, a_Rho);
      else if (img->mMode == Image::DIV || img->mMode == Image::DIVDT ||
               img->mMode == Image::CURLMAG || img->mMode == Image::CURLMAGDT)
        img->computeImageDivCurl(a_Up, a_U, a_Um, mDt, dminus);
      else if (img->mMode == Image::LAT || img->mMode == Image::LON)
        img->computeImageLatLon(mX, mY, mZ);
      else if (img->mMode == Image::TOPO) {
        if (topographyExists())
          img->copy2DArrayToImage(
              mTopo);  // save the raw topography; the smoothed is saved by the
                       // mode=grid with z=0
      } else if (img->mMode == Image::UZEXACT || img->mMode == Image::UXEXACT ||
                 img->mMode == Image::UYEXACT || img->mMode == Image::UXERR ||
                 img->mMode == Image::UYERR || img->mMode == Image::UZERR) {
        SW4_MARK_BEGIN("update_images::region 1");
        // Note: this is inefficient, the exact solution is computed everywhere,
        // and once for each
        //   EXACT or ERR image mode.
        vector<Sarray> Uex(mNumberOfGrids);
        vector<Sarray*> alpha(mNumberOfGrids);

        for (int g = 0; g < mNumberOfGrids; g++) {
          Uex[g].define(3, m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                        m_kStart[g], m_kEnd[g]);
          if (m_use_attenuation) {
            alpha[g] = new Sarray[m_number_mechanisms];
            for (int a = 0; a < m_number_mechanisms; a++)
              alpha[g][a].define(3, m_iStart[g], m_iEnd[g], m_jStart[g],
                                 m_jEnd[g], m_kStart[g], m_kEnd[g]);
          }
        }
        exactSol(time, Uex, alpha, a_sources);
        if (img->mMode == Image::UXERR || img->mMode == Image::UYERR ||
            img->mMode == Image::UZERR) {
          for (int g = 0; g < mNumberOfGrids; g++) {
            size_t n = static_cast<size_t>(Uex[g].npts());
            int nc = Uex[g].ncomp();
            float_sw4* uxp = Uex[g].c_ptr();
            float_sw4* up = a_Up[g].c_ptr();
#pragma omp parallel for
            for (size_t i = 0; i < n * nc; i++) uxp[i] = up[i] - uxp[i];
          }
        }
        if (img->mMode == Image::UXEXACT || img->mMode == Image::UXERR)
          img->computeImageQuantity(Uex, 1);
        else if (img->mMode == Image::UYEXACT || img->mMode == Image::UYERR)
          img->computeImageQuantity(Uex, 2);
        else if (img->mMode == Image::UZEXACT || img->mMode == Image::UZERR)
          img->computeImageQuantity(Uex, 3);
        Uex.clear();
        if (m_use_attenuation) {
          for (int g = 0; g < mNumberOfGrids; g++) delete[] alpha[g];
        }
        SW4_MARK_END("update_images::region 1");
      } else if (img->mMode == Image::GRIDX || img->mMode == Image::GRIDY ||
                 img->mMode == Image::GRIDZ)
        img->computeImageGrid(mX, mY, mZ);
      else if (img->mMode == Image::MAGDUDT) {
        if (dminus)
          img->computeImageMagdt(a_Up, a_U, mDt);
        else
          img->computeImageMagdt(a_Up, a_Um, 2 * mDt);
      } else if (img->mMode == Image::HMAGDUDT) {
        if (dminus)
          img->computeImageHmagdt(a_Up, a_U, mDt);
        else
          img->computeImageHmagdt(a_Up, a_Um, 2 * mDt);
      } else if (img->mMode == Image::MAG)
        img->computeImageMag(a_Up);
      else if (img->mMode == Image::HMAG)
        img->computeImageHmag(a_Up);
      //         else if(img->mMode == Image::QS )
      //	 {
      //	    if (usingAttenuation())
      //	       img->computeImageQuantity(mQs, 1);
      //       }
      //       else if(img->mMode == Image::QP )
      //       {
      // 	if (usingAttenuation())
      // 	  img->computeImageQuantity(mQp, 2);
      //       }
      // 	img->computeDivergence();
      // 	if (m_forcing->knows_exact())
      //         {
      //           maxerr = img->computeImageErrorDebug(0);
      //           if (proc_zero())
      //             printf("maxErr DIV %f @ %fs\n",maxerr,time );
      //         }
      //       }
      //       else if(img->mMode == Image::CURL )
      //       {
      //         img->computeCurlMagnitude();
      // 	if (m_forcing->knows_exact())
      //         {
      //           maxerr = img->computeImageErrorDebug(2);
      //           if (proc_zero())
      //             printf("maxErr CURL %f @ %f s\n",maxerr,time );
      //         }
      //       }
      //       else if(img->mMode == Image::VELDIV )
      //       {
      // 	img->computeVelDivergence();
      // 	if (m_forcing->knows_exact())
      //         {
      //           maxerr = img->computeImageErrorDebug(1);
      //           if (proc_zero())
      //             printf("maxErr VELDIV %f @ %f s\n",maxerr,time );
      //         }
      //       }
      //       else if(img->mMode == Image::VELCURL )
      //       {
      //         img->computeVelCurlMagnitude();
      // 	if (m_forcing->knows_exact())
      //         {
      //           maxerr = img->computeImageErrorDebug(3);
      //           if (proc_zero())
      //             printf("maxErr VELCURL %f @ %f s\n",maxerr,time );
      //         }
      //       }
      //       else if(img->mMode == Image::VELMAG )
      //       {
      //         img->computeVelocityMagnitude();
      //       }
      //       else if(img->mMode == Image::HVEL )
      //       {
      //         img->computeHorizontalVelocityMagnitude();
      //       }
      // // the following two modes converts mu,lambda,rho into vp and vs
      //       else if(img->mMode == Image::S )
      //       {
      // 	img->computeImageS(mMu, mLambda, mRho);
      //       }
      //       else if(img->mMode == Image::TOPO )
      //       {
      //         if (topographyExists())
      // 	  img->copy2DArrayToImage(mTopo); // save the raw topography;
      // the smoothed is saved by the mode=grid with z=0
      //       }
      //       else if(img->mMode == Image::LAT )
      //       {
      // 	img->evaluateLatLonImage(mX, mY, mZ, 1);
      //       }
      //       else if(img->mMode == Image::LON )
      //       {
      // 	img->evaluateLatLonImage(mX, mY, mZ, 2);
      //       }
      //       else if(img->mMode == Image::GRID )
      //       {
      // // get the base name
      // 	string filePrefix = img->mFilePrefix;
      // 	if (img->getOrientation() == Image::X) // save y and
      // z-coordinates
      // 	{
      // // copy the y-component
      // 	  img->evaluateGridImage(mX, mY, mZ, 2); // save the y-component
      // // append a "Y" to the file name
      // 	  img->mFilePrefix = filePrefix + "Y"; // will the "Y"
      // accumulate?
      // // save the file
      //        img->writeImagePlane_2(currentTimeStep, mPath); // save the grid
      //        image
      // 	}
      // 	else if (img->getOrientation() == Image::Y) // save x and
      // z-coordinates
      // 	{
      // // copy the x-component
      // 	  img->evaluateGridImage(mX, mY, mZ, 1); // save the y-component
      // // append a "X" to the file name
      // 	  img->mFilePrefix = filePrefix + "X";
      // // save the file
      //        img->writeImagePlane_2(currentTimeStep, mPath); // save the grid
      //        image
      // 	}
      // // copy the z-component
      // 	img->evaluateGridImage(mX, mY, mZ, 3); // save the z-component
      // // append a "Z" to the file name
      // 	img->mFilePrefix = filePrefix + "Z";
      // // the file is saved below

      //       }
      //      else
      //      else if (!img->mMode == Image::HMAXDUDT || !img->mMode ==
      //      Image::VMAXDUDT
      //	      || !img->mMode == Image::HMAX   || !img->mMode ==
      // Image::VMAX )
      else if (!(img->mMode == Image::HMAXDUDT ||
                 img->mMode == Image::VMAXDUDT || img->mMode == Image::HMAX ||
                 img->mMode == Image::VMAX || img->mMode == Image::GRADRHO ||
                 img->mMode == Image::GRADMU ||
                 img->mMode == Image::GRADLAMBDA ||
                 img->mMode == Image::GRADP || img->mMode == Image::GRADS)) {
        if (proc_zero()) {
          //	  printf("Can only write ux, uy, uz, mu, rho, lambda, uxerr,
          // uyerr, uzerr- remove once completely implemented\n");
          printf(
              "Can only write ux, uy, uz, mu, rho, lambda: - remove once "
              "completely implemented\n");
          printf("I can not print data of type %i\n", img->mMode);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // write the image plane on file
      double t[3];
      t[0] = t[1] = t[2] = MPI_Wtime();
      img->writeImagePlane_2(currentTimeStep - td, mPath[event],
                             time - td * mDt);
      t[2] = MPI_Wtime();

      // output timing info?
      if (m_iotiming) {
        t[0] = t[1] - t[0];
        t[1] = t[2] - t[1];

        double tmp[2];
        tmp[0] = t[0];
        tmp[1] = t[1];
        MPI_Reduce(tmp, t, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (proc_zero()) {
          cout << "Maximum write time:";
          cout << " (using Bjorn's I/O library) " << t[1]
               << " seconds. (<=" << m_nwriters << " procs writing)";
          cout << endl;
        }  // end if proc_zero
      }    // end if iotiming

    }  // end if time to write
  }    // end for all images

  // Volume images
  // for (unsigned int fIndex = 0; fIndex < mImage3DFiles.size(); ++fIndex)
  // {
  //   Image3D* img = mImage3DFiles[fIndex];
  //   if(img->timeToWrite(time, currentTimeStep, mDt ) )
  //   {
  //     img->compute_image( );
  //     img->write_images( currentTimeStep, mPath );
  //   }
  // }
}  // end EW::update_images()

//-----------------------------------------------------------------------
void EW::addImage(Image* i) { mImageFiles.push_back(i); }

//-----------------------------------------------------------------------
void EW::addImage3D(Image3D* i) { mImage3DFiles.push_back(i); }
//-----------------------------------------------------------------------
void EW::addESSI3D(ESSI3D* i) { mESSI3DFiles.push_back(i); }
//-----------------------------------------------------------------------
void EW::addSfileOutput(SfileOutput* i) { mSfiles.push_back(i); }
//-----------------------------------------------------------------------
void EW::initialize_image_files() {
  // In case of multiple events, prepare maximum number of time steps
  int maxNumberOfTimeSteps = 0;
  for (int e = 0; e < m_nevent; e++)
    maxNumberOfTimeSteps = mNumberOfTimeSteps[e] > maxNumberOfTimeSteps
                               ? mNumberOfTimeSteps[e]
                               : maxNumberOfTimeSteps;

  // Image planes
  Image::setSteps(maxNumberOfTimeSteps);
  for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex) {
    mImageFiles[fIndex]->computeGridPtIndex();
    mImageFiles[fIndex]->allocatePlane();
  }

  for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
    mImageFiles[fIndex]->associate_gridfiles(mImageFiles);

  for (unsigned int fIndex = 0; fIndex < mImageFiles.size(); ++fIndex)
    if (mImageFiles[fIndex]->mMode == Image::GRIDX ||
        mImageFiles[fIndex]->mMode == Image::GRIDY ||
        mImageFiles[fIndex]->mMode == Image::GRIDZ)
      mImageFiles[fIndex]->computeImageGrid(mX, mY, mZ);

  // Volume images
  Image3D::setSteps(maxNumberOfTimeSteps);
  for (unsigned int fIndex = 0; fIndex < mImage3DFiles.size(); ++fIndex)
    mImage3DFiles[fIndex]->setup_images();
  ESSI3D::setSteps(maxNumberOfTimeSteps);
  //   ESSI3D::setSteps(mNumberOfTimeSteps);
  for (unsigned int fIndex = 0; fIndex < mESSI3DFiles.size(); ++fIndex)
    mESSI3DFiles[fIndex]->setup();

  SfileOutput::setSteps(maxNumberOfTimeSteps);
  //   SfileOutput::setSteps(mNumberOfTimeSteps);
  for (unsigned int fIndex = 0; fIndex < mSfiles.size(); ++fIndex)
    mSfiles[fIndex]->setup_images();
}

//-----------------------------------------------------------------------
void EW::set_sg_thickness(int n_gp) {
  m_sg_gp_thickness = n_gp;
  m_use_sg_width = false;  // will be changed to true once the number of gp has
                           // been converted to a physical width
  if (m_myRank == 0)
    cout << "Default Supergrid thickness has been tuned; # grid points = "
         << m_sg_gp_thickness << " grid sizes" << endl;
}

//-----------------------------------------------------------------------
void EW::set_sg_width(float_sw4 sg_width) {
  m_supergrid_width = sg_width;
  m_use_sg_width = true;
  if (m_myRank == 0)
    cout << "Default Supergrid thickness has been tuned; width = "
         << m_supergrid_width << " meters" << endl;
}

//-----------------------------------------------------------------------
void EW::set_sg_damping(float_sw4 damp_coeff) {
  m_supergrid_damping_coefficient = damp_coeff;
  if (m_myRank == 0)
    cout << "Default Supergrid damping coefficient has been tuned; damping "
            "coefficient = "
         << m_supergrid_damping_coefficient << endl;
}

//-----------------------------------------------------------------------
void EW::set_global_bcs(boundaryConditionType bct[6]) {
  for (int i = 0; i < 6; i++) mbcGlobalType[i] = bct[i];
  mbcsSet = true;

  //  cout << "mbcGlobalType = " << mbcGlobalType[0] << "," << mbcGlobalType[1]
  //  << "," << mbcGlobalType[2] << "," << mbcGlobalType[3] << "," <<
  //  mbcGlobalType[4] << "," << mbcGlobalType[5] << endl;
}

//-----------------------------------------------------------------------
void EW::set_prefilter(FilterType passband, int order, int passes,
                       float_sw4 fc1, float_sw4 fc2) {
  m_prefilter_sources = true;
  // we could build the filter object right here...
  m_filter_ptr = new Filter(passband, order, passes, fc1, fc2);
  m_filterobs_ptr = new Filter(passband, order, passes, fc1, fc2);
}

//-----------------------------------------------------------------------
void EW::set_threshold_velocities(float_sw4 vpmin, float_sw4 vsmin) {
  m_useVelocityThresholds = true;
  m_vpMin = vpmin;
  m_vsMin = vsmin;
}

//-----------------------------------------------------------------------
void EW::average_speeds(float_sw4& cp, float_sw4& cs) {
  cp = 0;
  cs = 0;
  float_sw4 htop = mGridSize[mNumberOfGrids - 1];
  float_sw4 hbot = mGridSize[0];
  for (int g = 0; g < mNumberOfGrids; g++) {
    int istart = m_iStart[g];
    int iend = m_iEnd[g];
    int jstart = m_jStart[g];
    int jend = m_jEnd[g];
    int kstart = m_kStart[g];
    int kend = m_kEnd[g];
    float_sw4 h = mGridSize[g];
    int nsgxy = static_cast<int>(0.5 + m_sg_gp_thickness);
    int nsgz = static_cast<int>(0.5 + m_sg_gp_thickness);
    if (m_use_sg_width) {
      nsgxy = static_cast<int>(ceil(m_supergrid_width / h));
      nsgz = static_cast<int>(ceil(m_supergrid_width / h));
    }
    //      int nsgxy = (int)(0.5+m_sg_gp_thickness*htop/h);
    //      int nsgz  = (int)(0.5+m_sg_gp_thickness*hbot/h);

    float_sw4 cpgrid, csgrid, npts;
    if (m_bcType[g][0] == bSuperGrid) istart = istart + nsgxy;
    if (m_bcType[g][1] == bSuperGrid) iend = iend - nsgxy;
    if (m_bcType[g][2] == bSuperGrid) jstart = jstart + nsgxy;
    if (m_bcType[g][3] == bSuperGrid) jend = jend - nsgxy;
    // Use finest spacing on top z-boundary
    if (m_bcType[g][4] == bSuperGrid) kstart = kstart + nsgxy;
    // Use coarsest spacing on bottom z-boundary
    if (m_bcType[g][5] == bSuperGrid) kend = kend - nsgz;
    float_sw4* mu_ptr = mMu[g].c_ptr();
    float_sw4* lambda_ptr = mLambda[g].c_ptr();
    float_sw4* rho_ptr = mRho[g].c_ptr();
    if (m_croutines)
      velsum(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g], &m_kStart[g],
             &m_kEnd[g], &istart, &iend, &jstart, &jend, &kstart, &kend, mu_ptr,
             lambda_ptr, rho_ptr, &cpgrid, &csgrid, &npts);
    else {
      size_t nptssizet;
      velsum_ci(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                m_kEnd[g], istart, iend, jstart, jend, kstart, kend, mu_ptr,
                lambda_ptr, rho_ptr, cpgrid, csgrid, nptssizet);
      npts = nptssizet;
    }
    float_sw4 cpgridtmp = cpgrid;
    float_sw4 csgridtmp = csgrid;
    float_sw4 nptstmp = npts;
    MPI_Allreduce(&cpgridtmp, &cpgrid, 1, m_mpifloat, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&csgridtmp, &csgrid, 1, m_mpifloat, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nptstmp, &npts, 1, m_mpifloat, MPI_SUM, MPI_COMM_WORLD);
    cp = cp + cpgrid / npts;
    cs = cs + csgrid / npts;
  }
  cp = cp / mNumberOfGrids;
  cs = cs / mNumberOfGrids;
}

//-----------------------------------------------------------------------
void EW::layered_speeds(vector<float_sw4>& cp, vector<float_sw4>& z) {
  // Sample material on a uniform grid in the z direction with N points:
  int N = 200;
  float_sw4* zv = new float_sw4[N + 1];
  float_sw4* cpv = new float_sw4[N + 1];
  float_sw4 h = (m_global_zmax - m_global_zmin) / N;
  for (int i = 0; i <= N; i++) zv[i] = m_global_zmin + i * h;
  float_sw4 x0 =
      0.5 * ((m_iStart[0] - 1) * mGridSize[0] + (m_iEnd[0] - 1) * mGridSize[0]);
  float_sw4 y0 =
      0.5 * ((m_jStart[0] - 1) * mGridSize[0] + (m_jEnd[0] - 1) * mGridSize[0]);
  int i0, j0, k0, g0;

  for (int i = 0; i <= N; i++) {
    computeNearestGridPoint(i0, j0, k0, g0, x0, y0, zv[i]);
    cpv[i] = sqrt((2 * mMu[g0](i0, j0, k0) + mLambda[g0](i0, j0, k0)) /
                  mRho[g0](i0, j0, k0));
  }
  //   for( int b=0 ; b < m_mtrlblocks.size() ; b++ )
  //   {
  //      for( int i=0 ; i <= N ; i++ )
  //	 m_mtrlblocks[b]->get_material_pt( x0, y0, zv[i], rho, cs, cpv[i], qs,
  // qp );
  //   }
  cp.push_back(cpv[0]);
  int j = 1;
  int l = 0;
  float_sw4 tol = 1e-4;
  while (j < N) {
    while (fabs(cp[l] - cpv[j]) < tol * cp[l] && (j < N)) j++;
    if (j < N) {
      cp.push_back(cpv[j]);
      z.push_back((zv[j] + zv[j - 1]) / 2);
      l++;
    }
  }
  delete[] cpv;
  delete[] zv;
}

//-----------------------------------------------------------------------
void EW::testsourcediff(vector<Source*> GlobalSources, float_sw4 gradient[11],
                        float_sw4 hessian[121]) {
  for (int m = 0; m < 11; m++) {
    gradient[m] = 0;
    for (int j = 0; j < 11; j++) hessian[m + 11 * j] = 0;
  }
  vector<GridPointSource*> gpsources;
  GlobalSources[0]->set_grid_point_sources4(this, gpsources);

  vector<Sarray> kappa, eta;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  kappa.resize(mNumberOfGrids);
  eta.resize(mNumberOfGrids);
  for (int g = 0; g < mNumberOfGrids; g++) {
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];

    kappa[g].define(3, ifirst, ilast, jfirst, jlast, kfirst, klast);
    eta[g].define(3, ifirst, ilast, jfirst, jlast, kfirst, klast);
    kappa[g].set_value(1.0);
    eta[g].set_value(1.0);
  }
  //   for( int m=0 ; m < gpsources.size() ; m++ )
  //   cout << "size of sources " << gpsources.size() << endl;
  for (int m = 0; m < gpsources.size() - 1; m++) {
    gpsources[m]->add_to_gradient(
        kappa, eta, 0.63, mDt, gradient, mGridSize, mJ[mNumberOfGrids - 1],
        topographyExists());  // mJ argument should be the whole vector of
                              // Sarrays
    gpsources[m]->add_to_hessian(kappa, eta, 0.63, mDt, hessian, mGridSize);
  }
}

//-----------------------------------------------------------------------
void EW::get_scalefactors(float_sw4 sf[11]) {
  for (int i = 0; i < 11; i++) sf[i] = m_scalefactors[i];
}

//-----------------------------------------------------------------------
bool EW::compute_sf() { return m_compute_scalefactors; }

//-----------------------------------------------------------------------
void EW::compute_guess(bool& guesspos, bool& guesst0fr, bool& guessmom,
                       bool& guessshifts, bool& output_seismograms) {
  guesspos = m_iniguess_pos;
  guesst0fr = m_iniguess_t0fr;
  guessmom = m_iniguess_mom;
  guessshifts = m_iniguess_shifts;
  output_seismograms = m_output_initial_seismograms;
}

//-----------------------------------------------------------------------
void EW::get_cgparameters(int& maxit, int& maxrestart, float_sw4& tolerance,
                          bool& fletcherreeves, int& stepselection,
                          bool& do_linesearch, int& varcase, bool& testing) {
  maxit = m_maxit;
  maxrestart = m_maxrestart;
  tolerance = m_tolerance;
  fletcherreeves = m_cgfletcherreeves;
  stepselection = m_cgstepselection;
  varcase = m_cgvarcase;
  do_linesearch = m_do_linesearch;
  testing = m_opt_testing;
}

//-----------------------------------------------------------------------
void EW::compute_energy(float_sw4 dt, bool write_file, vector<Sarray>& Um,
                        vector<Sarray>& U, vector<Sarray>& Up, int step,
                        int event) {
  SW4_MARK_FUNCTION;
  // Compute energy
  float_sw4 energy = 0;
  for (int g = 0; g < mNumberOfGrids; g++) {
    int istart = m_iStartInt[g];
    int iend = m_iEndInt[g];
    int jstart = m_jStartInt[g];
    int jend = m_jEndInt[g];
    int kstart = m_kStartInt[g];
    int kend = m_kEndInt[g];
    float_sw4* up_ptr = Up[g].c_ptr();
    float_sw4* u_ptr = U[g].c_ptr();
    float_sw4* um_ptr = Um[g].c_ptr();
    float_sw4* rho_ptr = mRho[g].c_ptr();
    float_sw4 locenergy;
    int* onesided_ptr = m_onesided[g];
    //   if (topographyExists() && g == mNumberOfGrids - 1) {
    if (topographyExists() && g >= mNumberOfCartesianGrids) {
      if (m_gridGenerator->curviCartIsSmooth(mNumberOfGrids -
                                             mNumberOfCartesianGrids) &&
          g == mNumberOfCartesianGrids)
        kend--;

      if (m_croutines)
        energy4c_ci(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                    m_kEnd[g], istart, iend, jstart, jend, kstart, kend,
                    onesided_ptr, um_ptr, u_ptr, up_ptr, rho_ptr, mJ[g].c_ptr(),
                    locenergy);
      else
        energy4c(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
                 &m_kStart[g], &m_kEnd[g], &istart, &iend, &jstart, &jend,
                 &kstart, &kend, onesided_ptr, um_ptr, u_ptr, up_ptr, rho_ptr,
                 mJ[g].c_ptr(), &locenergy);
    } else {
      if (m_croutines)
        energy4_ci(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                   m_kEnd[g], istart, iend, jstart, jend, kstart, kend,
                   onesided_ptr, um_ptr, u_ptr, up_ptr, rho_ptr, mGridSize[g],
                   m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g], locenergy);
      else
        energy4(&m_iStart[g], &m_iEnd[g], &m_jStart[g], &m_jEnd[g],
                &m_kStart[g], &m_kEnd[g], &istart, &iend, &jstart, &jend,
                &kstart, &kend, onesided_ptr, um_ptr, u_ptr, up_ptr, rho_ptr,
                &mGridSize[g], m_sg_str_x[g], m_sg_str_y[g], m_sg_str_z[g],
                &locenergy);
    }
    energy += locenergy;
  }
  energy /= (dt * dt);
  float_sw4 energytmp = energy;
  MPI_Allreduce(&energytmp, &energy, 1, m_mpifloat, MPI_SUM,
                m_cartesian_communicator);
  m_energy_test->record_data(energy, step, write_file, m_myRank, mPath[event]);
}
//-----------------------------------------------------------------------
float_sw4 EW::scalarProduct(vector<Sarray>& U, vector<Sarray>& V) {
  // Compute weighted scalar product between composite grid functions U and C
  //
  // NOTE: assumes a Cartesian grid with super-grid stretching
  //
  float_sw4 s_prod = 0;
  for (int g = 0; g < mNumberOfGrids; g++) {
    int istart = m_iStartInt[g];
    int iend = m_iEndInt[g];
    int jstart = m_jStartInt[g];
    int jend = m_jEndInt[g];
    int kstart = m_kStartInt[g];
    int kend = m_kEndInt[g];
    float_sw4* u_ptr = U[g].c_ptr();
    float_sw4* v_ptr = V[g].c_ptr();
    float_sw4 loc_s_prod;
    int* onesided_ptr = m_onesided[g];
    if (m_croutines)
      scalar_prod_ci(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                     m_kStart[g], m_kEnd[g], istart, iend, jstart, jend, kstart,
                     kend, onesided_ptr, u_ptr, v_ptr, m_sg_str_x[g],
                     m_sg_str_y[g], m_sg_str_z[g], loc_s_prod);
    else
      scalar_prod(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                  m_kEnd[g], istart, iend, jstart, jend, kstart, kend,
                  onesided_ptr, u_ptr, v_ptr, m_sg_str_x[g], m_sg_str_y[g],
                  m_sg_str_z[g], &loc_s_prod);
    s_prod += loc_s_prod;
  }
  // output my sum
  //   printf("scalarProd: myRank=%d, my_s_prod=%e\n", m_myRank, s_prod);

  float_sw4 s_prod_tmp = s_prod;
  MPI_Allreduce(&s_prod_tmp, &s_prod, 1, m_mpifloat, MPI_SUM,
                m_cartesian_communicator);
  return s_prod;
}

//-----------------------------------------------------------------------
void EW::get_utc(int utc[7], int event) const {
  for (int c = 0; c < 7; c++) utc[c] = m_utc0[event][c];
}

//-----------------------------------------------------------------------
void EW::print_utc(int e) {
  if (proc_zero()) {
    printf("EW reference UTC is  %02i/%02i/%i:%i:%i:%i.%i\n", m_utc0[e][1],
           m_utc0[e][2], m_utc0[e][0], m_utc0[e][3], m_utc0[e][4], m_utc0[e][5],
           m_utc0[e][6]);
  }
}

//-----------------------------------------------------------------------
void EW::extractTopographyFromGridFile(string a_topoFileName) {
  if (proc_zero()) cout << "***inside extractTopographyFromGridFile***" << endl;

  //----------------------------------------------
  // Check user specified file name. Abort if they are not there or not readable
  //----------------------------------------------
  CHECK_INPUT(access(a_topoFileName.c_str(), R_OK) == 0,
              "No read permission on topography grid file: " << a_topoFileName);

  int topLevel = mNumberOfGrids - 1;

  double x, y;
  double lat, lon, elev;

  // 1. read the grid file
  int Nlon, Nlat, i, j;
  Sarray gridElev;
  double *latv, *lonv;

  FILE* gridfile = fopen(a_topoFileName.c_str(), "r");

  fscanf(gridfile, "%i %i", &Nlon, &Nlat);
  gridElev.define(1, 1, Nlon, 1, Nlat, 1, 1);
  latv = new double[Nlat + 1];
  lonv = new double[Nlon + 1];

  for (j = 1; j <= Nlat; j++)
    for (i = 1; i <= Nlon; i++)
      fscanf(gridfile, "%le %le %le", &lonv[i], &latv[j],
             &gridElev(1, i, j, 1));
  fclose(gridfile);

  if (proc_zero()) printf("Nlon=%i Nlat=%i\n", Nlon, Nlat);

  double lonMax = -180.0, lonMin = 180.0, latMax = -90.0, latMin = 90.0,
         elevMax = -1e10, elevMin = 1e10;
  for (i = 1; i <= Nlon; i++) {
    if (lonv[i] < lonMin) lonMin = lonv[i];
    if (lonv[i] > lonMax) lonMax = lonv[i];
  }
  for (i = 1; i <= Nlat; i++) {
    if (latv[i] < latMin) latMin = latv[i];
    if (latv[i] > latMax) latMax = latv[i];
  }

  for (i = 1; i <= Nlon; i++)
    for (j = 1; j <= Nlat; j++) {
      if (gridElev(1, i, j, 1) < elevMin) elevMin = gridElev(1, i, j, 1);
      if (gridElev(1, i, j, 1) > elevMax) elevMax = gridElev(1, i, j, 1);
    }
  if (proc_zero())
    printf(
        "lonMin=%e, lonMax=%e\nlatMin=%e, latMax=%e\nelevMin=%e, evalMax=%e\n",
        lonMin, lonMax, latMin, latMax, elevMin, elevMax);

  // If the lat vector is not in increasing order, we need to reorder it
  if (latv[1] > latv[Nlat]) {
    // tmp
    if (proc_zero()) printf("Reordering the latitude vector...\n");

    for (j = 1; j <= Nlat / 2; j++) {
      lat = latv[Nlat + 1 - j];
      latv[Nlat + 1 - j] = latv[j];
      latv[j] = lat;

      for (i = 1; i <= Nlon; i++) {
        elev = gridElev(1, i, Nlat + 1 - j, 1);
        gridElev(1, i, Nlat + 1 - j, 1) = gridElev(1, i, j, 1);
        gridElev(1, i, j, 1) = elev;
      }
    }  // end for j
  }    // end if latv[1] > latv[Nlat]

  // If the lon vector is not in increasing order, we need to reorder it
  if (lonv[1] > lonv[Nlon]) {
    // tmp
    if (proc_zero()) printf("Reordering the longitude vector...\n");

    for (i = 1; i <= Nlon / 2; i++) {
      lon = lonv[Nlon + 1 - i];
      lonv[Nlon + 1 - i] = lonv[i];
      lonv[i] = lon;

      for (j = 1; j <= Nlat; j++) {
        elev = gridElev(1, Nlon + 1 - i, j, 1);
        gridElev(1, Nlon + 1 - i, j, 1) = gridElev(1, i, j, 1);
        gridElev(1, i, j, 1) = elev;
      }
    }  // end for i
  }    // end if lonv[1] > lonv[Nlon]

  // 2. interpolate in the grid file to get elevations on the computational grid
  double deltaLat = (latMax - latMin) / Nlat;
  double deltaLon = (lonMax - lonMin) / Nlon;

  int i0, j0;

  for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i) {
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j) {
      x = (i - 1) * mGridSize[topLevel];
      y = (j - 1) * mGridSize[topLevel];

      computeGeographicCoord(
          x, y, lon,
          lat);  // the grid command defines the parameters in this mapping

      if (lat > latMax || lat < latMin || lon > lonMax || lon < lonMin) {
        printf(
            "ERROR: x=%e, y=%e corresponds to lon=%e, lat=%e which are outside "
            "the topography grid\n",
            x, y, lon, lat);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      i0 = 1 + (int)((lon - lonMin) / deltaLon);
      j0 = 1 + (int)((lat - latMin) / deltaLat);

      while (lon < lonv[i0] ||
             lonv[i0 + 1] < lon)  // should stop loop if i0 is out of bounds
      {
        if (lon < lonv[i0])
          i0--;
        else if (lon > lonv[i0 + 1])
          i0++;
      }

      while (lat < latv[j0] ||
             latv[j0 + 1] < lat)  // should stop loop if j0 is out of bounds
      {
        if (lat < latv[j0])
          j0--;
        else if (lat > latv[j0 + 1])
          j0++;
      }

      if (i0 > Nlon - 1) i0 = Nlon - 1;

      if (j0 > Nlat - 1) j0 = Nlat - 1;

      // test that we are inside the interval
      if (!(lonv[i0] <= lon && lon < lonv[i0 + 1])) {
        printf(
            "EW::extractTopographyFromGridFile: Fatal error: Unable to "
            "interpolate topography for lon=%e\n"
            "because it is outside the cell (lonv[%i]=%e, lonv[%i]=%e)\n",
            lon, i0, lonv[i0], i0 + 1, lonv[i0 + 1]);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      if (!(latv[j0] <= lat && lat < latv[j0 + 1])) {
        printf(
            "EW::extractTopographyFromGridFile: Fatal error: Unable to "
            "interpolate topography for lat=%e\n"
            "because it is outside the cell (latv[%i]=%e, latv[%i]=%e)\n",
            lat, j0, latv[j0], j0 + 1, latv[j0 + 1]);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // bi-cubic interpolation should make the surface smoother
      // shift the stencil if it is too close to the boundaries
      if (i0 < 2) i0 = 2;
      if (i0 > Nlon - 2) i0 = Nlon - 2;

      if (j0 < 2) j0 = 2;
      if (j0 > Nlat - 2) j0 = Nlat - 2;

      // local step sizes
      float_sw4 q = i0 + (lon - lonv[i0]) / (lonv[i0 + 1] - lonv[i0]);
      float_sw4 r = j0 + (lat - latv[j0]) / (latv[j0 + 1] - latv[j0]);

      float_sw4 Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1,
          tjp2;
      Qim1 = (q - i0) * (q - i0 - 1) * (q - i0 - 2) / (-6.);
      Qi = (q - i0 + 1) * (q - i0 - 1) * (q - i0 - 2) / (2.);
      Qip1 = (q - i0 + 1) * (q - i0) * (q - i0 - 2) / (-2.);
      Qip2 = (q - i0 + 1) * (q - i0) * (q - i0 - 1) / (6.);

      Rjm1 = (r - j0) * (r - j0 - 1) * (r - j0 - 2) / (-6.);
      Rj = (r - j0 + 1) * (r - j0 - 1) * (r - j0 - 2) / (2.);
      Rjp1 = (r - j0 + 1) * (r - j0) * (r - j0 - 2) / (-2.);
      Rjp2 = (r - j0 + 1) * (r - j0) * (r - j0 - 1) / (6.);

      tjm1 = Qim1 * gridElev(i0 - 1, j0 - 1, 1) + Qi * gridElev(i0, j0 - 1, 1) +
             Qip1 * gridElev(i0 + 1, j0 - 1, 1) +
             Qip2 * gridElev(i0 + 2, j0 - 1, 1);
      tj = Qim1 * gridElev(i0 - 1, j0, 1) + Qi * gridElev(i0, j0, 1) +
           Qip1 * gridElev(i0 + 1, j0, 1) + Qip2 * gridElev(i0 + 2, j0, 1);
      tjp1 = Qim1 * gridElev(i0 - 1, j0 + 1, 1) + Qi * gridElev(i0, j0 + 1, 1) +
             Qip1 * gridElev(i0 + 1, j0 + 1, 1) +
             Qip2 * gridElev(i0 + 2, j0 + 1, 1);
      tjp2 = Qim1 * gridElev(i0 - 1, j0 + 2, 1) + Qi * gridElev(i0, j0 + 2, 1) +
             Qip1 * gridElev(i0 + 1, j0 + 2, 1) +
             Qip2 * gridElev(i0 + 2, j0 + 2, 1);

      mTopo(i, j, 1) = Rjm1 * tjm1 + Rj * tj + Rjp1 * tjp1 + Rjp2 * tjp2;

      // bi-linear interpolation
      // // local step sizes
      //       xi = (lon - lonv[i0])/(lonv[i0+1]-lonv[i0]);
      //       eta = (lat - latv[j0])/(latv[j0+1]-latv[j0]);
      //       mTopo(i,j,1) = (1.0-eta)*( (1.0-xi)*gridElev(1,i0,j0,1) +
      //       xi*gridElev(1,i0+1,j0,1) ) +
      // 	eta*( (1.0-xi)*gridElev(1,i0,j0+1,1) +
      // xi*gridElev(1,i0+1,j0+1,1) );
    }
  }
  delete[] latv;
  delete[] lonv;
}

//-----------------------------------------------------------------------
void EW::extractTopographyFromCartesianFile(string a_topoFileName) {
  if (proc_zero())
    cout << "***inside extractTopographyFromCartesianFile***" << endl;

  //----------------------------------------------
  // Check user specified file name. Abort if they are not there or not readable
  //----------------------------------------------
  VERIFY2(access(a_topoFileName.c_str(), R_OK) == 0,
          "No read permission on topography grid file: " << a_topoFileName);

  // 1. read the grid file
  int Nx, Ny, i, j;
  Sarray gridElev;
  float_sw4 *yv, *xv;

  FILE* gridfile = fopen(a_topoFileName.c_str(), "r");

  fscanf(gridfile, "%i %i", &Nx, &Ny);
  gridElev.define(1, 1, Nx, 1, Ny, 1, 1);
  yv = new float_sw4[Ny + 1];
  xv = new float_sw4[Nx + 1];

  for (j = 1; j <= Ny; j++)
    for (i = 1; i <= Nx; i++)
      fscanf(gridfile, "%le %le %le", &xv[i], &yv[j], &gridElev(1, i, j, 1));
  fclose(gridfile);

  if (proc_zero()) printf("Nx=%i Ny=%i\n", Nx, Ny);

  float_sw4 xMax = -999e10, xMin = 999e10, yMax = -999e10, yMin = 999e10,
            elevMax = -1e10, elevMin = 1e10;
  for (i = 1; i <= Nx; i++) {
    if (xv[i] < xMin) xMin = xv[i];
    if (xv[i] > xMax) xMax = xv[i];
  }
  for (i = 1; i <= Ny; i++) {
    if (yv[i] < yMin) yMin = yv[i];
    if (yv[i] > yMax) yMax = yv[i];
  }
  // make sure that the topography grid covers the whole computational domain
  if (xMin > 0 || yMin > 0 || xMax < m_global_xmax || yMax < m_global_ymax) {
    if (proc_zero())
      printf(
          "ERROR: Cartesian topography grid with %e<=x<=%e and %e<=y<=%e\n"
          "does not cover the computational domain: 0<=x<=%e, 0<=y<=%e\n",
          xMin, xMax, yMin, yMax, m_global_xmax, m_global_ymax);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      if (gridElev(1, i, j, 1) < elevMin) elevMin = gridElev(1, i, j, 1);
      if (gridElev(1, i, j, 1) > elevMax) elevMax = gridElev(1, i, j, 1);
    }
  if (proc_zero())
    printf("xMin=%e, xMax=%e\nyMin=%e, yMax=%e\nelevMin=%e, evalMax=%e\n", xMin,
           xMax, yMin, yMax, elevMin, elevMax);

  float_sw4 xP, yP, elev;
  // If the yv vector is not in increasing order, we need to reorder it
  if (yv[1] > yv[Ny]) {
    // tmp
    if (proc_zero()) printf("Reordering the yv vector...\n");
    for (j = 1; j <= Ny / 2; j++) {
      yP = yv[Ny + 1 - j];
      yv[Ny + 1 - j] = yv[j];
      yv[j] = yP;

      for (i = 1; i <= Nx; i++) {
        elev = gridElev(1, i, Ny + 1 - j, 1);
        gridElev(1, i, Ny + 1 - j, 1) = gridElev(1, i, j, 1);
        gridElev(1, i, j, 1) = elev;
      }
    }
  }

  // If the xv vector is not in increasing order, we need to reorder it
  if (xv[1] > xv[Nx]) {
    // tmp
    if (proc_zero()) printf("Reordering the xv vector...\n");
    for (i = 1; i <= Nx / 2; i++) {
      xP = xv[Nx + 1 - i];
      xv[Nx + 1 - i] = xv[i];
      xv[i] = xP;

      for (j = 1; j <= Ny; j++) {
        elev = gridElev(1, Nx + 1 - i, j, 1);
        gridElev(1, Nx + 1 - i, j, 1) = gridElev(1, i, j, 1);
        gridElev(1, i, j, 1) = elev;
      }
    }
  }

  // 2. interpolate in the grid file to get elevations on the computational grid
  float_sw4 deltaY = (yMax - yMin) / Ny;
  float_sw4 deltaX = (xMax - xMin) / Nx;
  int topLevel = mNumberOfGrids - 1;
  float_sw4 hp =
      1.01 * mGridSize[topLevel];  // change this to 2 grid sizes because there
                                   // are double ghost points?

#pragma omp parallel for
  for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i) {
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j) {
      float_sw4 xP = (i - 1) * mGridSize[topLevel];
      float_sw4 yP = (j - 1) * mGridSize[topLevel];
      int i0, j0;
      bool xGhost = true, yGhost = true;
      if (yP > yMax + hp || yP < yMin - hp || xP > xMax + hp ||
          xP < xMin - hp) {
        mTopo(i, j, 1) = NO_TOPO;
        continue;
        // printf("ERROR: xP=%e, yP=%e is outside the topography grid by more
        // than a grid step\n", 	  xP, yP);
        //  MPI_Abort(MPI_COMM_WORLD, 1);
      }
      // Compute i0
      if (xP < xMin)
        i0 = 1;
      else if (xP > xMax)
        i0 = Nx - 1;
      else {
        xGhost = false;
        i0 = 1 + (int)((xP - xMin) / deltaX);
        if (i0 < 1) i0 = 1;
        if (i0 > Nx - 1) i0 = Nx - 1;
        // should stop loop if i0 is out of bounds
        while (i0 >= 1 && i0 <= Nx - 1 && (xP < xv[i0] || xv[i0 + 1] < xP)) {
          if (xP < xv[i0])
            i0--;
          else if (xP > xv[i0 + 1])
            i0++;
        }
      }

      // Compute j0
      if (yP < yMin)
        j0 = 1;
      else if (yP > yMax)
        j0 = Ny - 1;
      else {
        yGhost = false;
        j0 = 1 + (int)((yP - yMin) / deltaY);
        if (j0 < 1) j0 = 1;
        if (j0 > Ny - 1) j0 = Ny - 1;
        // should stop loop if j0 is out of bounds
        while (j0 >= 1 && j0 <= Ny - 1 && (yP < yv[j0] || yv[j0 + 1] < yP)) {
          if (yP < yv[j0])
            j0--;
          else if (yP > yv[j0 + 1])
            j0++;
        }
      }

      // enforce bounds again
      if (i0 < 1) i0 = 1;
      if (i0 > Nx - 1) i0 = Nx - 1;
      if (j0 < 1) j0 = 1;
      if (j0 > Ny - 1) j0 = Ny - 1;

      // test that we are inside the interval
      if (!xGhost && !(xv[i0] <= xP && xP <= xv[i0 + 1])) {
        printf(
            "EW::extractTopographyFromCartesianFile: Fatal error: Unable to "
            "interpolate topography for xP=%e\n"
            "because it is outside the cell (xv[%i]=%e, xv[%i]=%e)\n",
            xP, i0, xv[i0], i0 + 1, xv[i0 + 1]);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      if (!yGhost && !(yv[j0] <= yP && yP <= yv[j0 + 1])) {
        printf(
            "EW::extractTopographyFromCartesianFile: Fatal error: Unable to "
            "interpolate topography for yP=%e\n"
            "because it is outside the cell (yv[%i]=%e, yv[%i]=%e)\n",
            yP, j0, yv[j0], j0 + 1, yv[j0 + 1]);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // bi-cubic interpolation should make the surface smoother
      // shift the stencil if it is too close to the boundaries
      if (i0 < 2) i0 = 2;
      if (i0 > Nx - 2) i0 = Nx - 2;

      if (j0 < 2) j0 = 2;
      if (j0 > Ny - 2) j0 = Ny - 2;

      // local step sizes
      float_sw4 q = i0 + (xP - xv[i0]) / (xv[i0 + 1] - xv[i0]);
      float_sw4 r = j0 + (yP - yv[j0]) / (yv[j0 + 1] - yv[j0]);

      float_sw4 Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1,
          tjp2;
      Qim1 = (q - i0) * (q - i0 - 1) * (q - i0 - 2) / (-6.);
      Qi = (q - i0 + 1) * (q - i0 - 1) * (q - i0 - 2) / (2.);
      Qip1 = (q - i0 + 1) * (q - i0) * (q - i0 - 2) / (-2.);
      Qip2 = (q - i0 + 1) * (q - i0) * (q - i0 - 1) / (6.);

      Rjm1 = (r - j0) * (r - j0 - 1) * (r - j0 - 2) / (-6.);
      Rj = (r - j0 + 1) * (r - j0 - 1) * (r - j0 - 2) / (2.);
      Rjp1 = (r - j0 + 1) * (r - j0) * (r - j0 - 2) / (-2.);
      Rjp2 = (r - j0 + 1) * (r - j0) * (r - j0 - 1) / (6.);

      tjm1 = Qim1 * gridElev(i0 - 1, j0 - 1, 1) + Qi * gridElev(i0, j0 - 1, 1) +
             Qip1 * gridElev(i0 + 1, j0 - 1, 1) +
             Qip2 * gridElev(i0 + 2, j0 - 1, 1);
      tj = Qim1 * gridElev(i0 - 1, j0, 1) + Qi * gridElev(i0, j0, 1) +
           Qip1 * gridElev(i0 + 1, j0, 1) + Qip2 * gridElev(i0 + 2, j0, 1);
      tjp1 = Qim1 * gridElev(i0 - 1, j0 + 1, 1) + Qi * gridElev(i0, j0 + 1, 1) +
             Qip1 * gridElev(i0 + 1, j0 + 1, 1) +
             Qip2 * gridElev(i0 + 2, j0 + 1, 1);
      tjp2 = Qim1 * gridElev(i0 - 1, j0 + 2, 1) + Qi * gridElev(i0, j0 + 2, 1) +
             Qip1 * gridElev(i0 + 1, j0 + 2, 1) +
             Qip2 * gridElev(i0 + 2, j0 + 2, 1);

      mTopo(i, j, 1) = Rjm1 * tjm1 + Rj * tj + Rjp1 * tjp1 + Rjp2 * tjp2;
    }
  }
  delete[] yv;
  delete[] xv;
}

//-----------------------------------------------------------------------
void EW::extractTopographyFromImageFile(string a_topoFileName) {
  // plan:
  // 1: read entire array from image file
  // 2: copy the relevant part to the mTopo array

  // Check user specified file names. Abort if they are not there or not
  // readable
  VERIFY2(access(a_topoFileName.c_str(), R_OK) == 0,
          "No read permission on topo image file: " << a_topoFileName);
  FILE* fd =
      fopen(a_topoFileName.c_str(), "rb");  // perhaps the "b" is redundant

  // % Read header
  //    prec    =fread(fd,1,'int');
  //    npatches=fread(fd,1,'int');
  //    t       =fread(fd,1,'double');
  //    plane   =fread(fd,1,'int');
  //    coord   =fread(fd,1,'double');
  //    mode    =fread(fd,1,'int');
  //    gridinfo=fread(fd,1,'int');
  //    timecreated=fread(fd,[1 25],'uchar');
  //    timestring=num2str(timecreated,'%c');
  //    mstr=getimagemodestr(mode);
  int prec, npatches, plane, mode, gridinfo;
  double time, coord;
  char timecreated[25];
  size_t nread;
  nread = fread(&prec, sizeof(int), 1, fd);
  nread = fread(&npatches, sizeof(int), 1, fd);
  nread = fread(&time, sizeof(double), 1, fd);
  nread = fread(&plane, sizeof(int), 1, fd);
  nread = fread(&coord, sizeof(double), 1, fd);
  nread = fread(&mode, sizeof(int), 1, fd);
  nread = fread(&gridinfo, sizeof(int), 1, fd);
  nread = fread(timecreated, sizeof(char), 25, fd);

  // % Display header
  //    if verbose == 1
  if (proc_zero()) {
    printf(
        "TopoImage header: prec=%i, npatches=%i, time=%e, plane=%i, coord=%e, "
        "mode=%i, gridinfo=%i\n",
        prec, npatches, time, plane, coord, mode, gridinfo);
    printf("                  timecreated=%s\n", timecreated);
  }

  // rudimentary checks
  if ((prec == 4 || prec == 8) && npatches == 1 && plane == 2 &&
      mode == Image::TOPO) {
    if (proc_zero()) printf("Header seems ok...\n");
  } else {
    if (proc_zero())
      printf(
          "Header for topo image is weird: prec=%i, npatches=%i, plane=%i, "
          "mode=%i\n",
          prec, npatches, plane, mode);
  }
  // header for each patch (should only be one patch in these files)
  // h(p) = fread(fd,1,'double');
  // zmin(p) = fread(fd,1,'double');
  // ib(p) = fread(fd,1,'int');
  // ni(p) = fread(fd,1,'int');
  // jb(p) = fread(fd,1,'int');
  // nj(p) = fread(fd,1,'int');
  // disp(['    patch nr ' num2str(p) ' has h = ' num2str(h(p)) ' zmin = '
  // num2str(zmin(p))]);
  double h, zmin;
  int ib, ni, jb, nj;
  nread = fread(&h, sizeof(double), 1, fd);
  nread = fread(&zmin, sizeof(double), 1, fd);
  nread = fread(&ib, sizeof(int), 1, fd);
  nread = fread(&ni, sizeof(int), 1, fd);
  nread = fread(&jb, sizeof(int), 1, fd);
  nread = fread(&nj, sizeof(int), 1, fd);
  if (proc_zero())
    printf("Patch info: h=%e, zmin=%e, ib=%i, ni=%i, jb=%i, nj=%i\n", h, zmin,
           ib, ni, jb, nj);
  // % Read data
  // readz = 0;
  // if pnr <= npatches
  //    for p=1:pnr-1
  // 	fseek(fd,(ni(p)-ib(p)+1)*(nj(p)-jb(p)+1)*prec,'cof');
  //    end;
  //    if prec == 4
  //       im0 = fread(fd,[ni(pnr)-ib(pnr)+1 nj(pnr)-jb(pnr)+1],'float');
  //    else
  //       im0 = fread(fd,[ni(pnr)-ib(pnr)+1 nj(pnr)-jb(pnr)+1],'float_sw4');
  //    end;
  // since there is only one patch, we don't need any fseek
  int npts = (ni - ib + 1) * (nj - jb + 1);
  int topLevel = mNumberOfGrids - 1;
  if (prec == 4)  // float
  {
    float* data_flt = new float[npts];
    nread = fread(data_flt, sizeof(float), npts, fd);
    CHECK_INPUT(npts == nread,
                "Number of image floats read: "
                    << nread << ", is different from header npts: " << npts);
// copy data to local mTopo array
#pragma omp parallel for
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
      for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i) {
        if (i >= ib && i <= ib + ni - 1 && j >= jb && j <= jb + nj - 1) {
          mTopo(i, j, 1) = (float_sw4)data_flt[i - ib + (j - jb) * ni];
        } else {
          mTopo(i, j, 1) = NO_TOPO;
        }
      }

    // cleanup local storage
    delete[] data_flt;
  } else if (prec == 8)  // double
  {
    double* data_dbl = new double[npts];
    nread = fread(data_dbl, sizeof(double), npts, fd);
    CHECK_INPUT(npts == nread,
                "Number of image double read: "
                    << nread << ", is different from header npts: " << npts);
// copy data to local mTopo array
#pragma omp parallel for
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
      for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i) {
        if (i >= ib && i <= ib + ni - 1 && j >= jb && j <= jb + nj - 1) {
          mTopo(i, j, 1) = data_dbl[i - ib + (j - jb) * ni];
        } else {
          mTopo(i, j, 1) = NO_TOPO;
        }
      }

    // cleanup local storage
    delete[] data_dbl;
  }

  fclose(fd);
  if (proc_zero()) printf("Topo image read ok\n");
}

//-----------------------------------------------------------------------
void EW::extractTopographyFromRfile(std::string a_topoFileName) {
  std::string rname = "EW::extractTopographyFromRfile";
  Sarray gridElev;
  int fd = open(a_topoFileName.c_str(), O_RDONLY);
  if (fd != -1) {
    // ---------- magic number
    int magic;
    size_t nr = read(fd, &magic, sizeof(int));
    if (nr != sizeof(int)) {
      cout << rname << " Error reading magic number, nr= " << nr << "bytes read"
           << endl;
      close(fd);
      return;
    }

    Byteswapper bswap;
    int onesw = 1;
    bswap.byte_rev(&onesw, 1, "int");
    bool swapbytes;
    if (magic == 1)
      swapbytes = false;
    else if (magic == onesw)
      swapbytes = true;
    else {
      cout << rname << "error could not determine byte order on file "
           << a_topoFileName << " magic number is " << magic << endl;
      close(fd);
      return;
    }

    // ---------- precision
    int prec = 0;
    nr = read(fd, &prec, sizeof(int));
    if (nr != sizeof(int)) {
      cout << rname << " Error reading prec, nr= " << nr << "bytes read"
           << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(&prec, 1, "int");
    int flsize = 4;
    if (prec == 8)
      flsize = sizeof(double);
    else if (prec == 4)
      flsize = sizeof(float);

    // ---------- attenuation on file ?
    int att;
    nr = read(fd, &att, sizeof(int));
    if (nr != sizeof(int)) {
      cout << rname << " Error reading att, nr= " << nr << "bytes read" << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(&att, 1, "int");

    // ---------- azimuth on file
    double alpha;
    nr = read(fd, &alpha, sizeof(double));
    if (nr != sizeof(double)) {
      cout << rname << " Error reading alpha, nr= " << nr << "bytes read"
           << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(&alpha, 1, "double");

    CHECK_INPUT(
        fabs(alpha - mGeoAz) < 1e-6,
        "ERROR: Rfile azimuth must be equal to coordinate system azimuth"
            << " azimuth on rfile = " << alpha
            << " azimuth of coordinate sytem = " << mGeoAz
            << " difference = " << alpha - mGeoAz);

    // ---------- origin on file
    float_sw4 lon0, lat0;
    nr = read(fd, &lon0, sizeof(double));
    if (nr != sizeof(double)) {
      cout << rname << " Error reading lon0, nr= " << nr << "bytes read"
           << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(&lon0, 1, "double");

    nr = read(fd, &lat0, sizeof(double));
    if (nr != sizeof(double)) {
      cout << rname << " Error reading lat0, nr= " << nr << "bytes read"
           << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(&lat0, 1, "double");

    // ---------- length of projection string
    int len;
    nr = read(fd, &len, sizeof(int));
    if (nr != sizeof(int)) {
      cout << rname << " Error reading len, nr= " << nr << "bytes read" << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(&len, 1, "int");

    // ---------- skip projection string
    nr = lseek(fd, len * sizeof(char), SEEK_CUR);

    // ---------- number of blocks on file
    int npatches;
    nr = read(fd, &npatches, sizeof(int));
    if (nr != sizeof(int)) {
      cout << rname << " Error reading npatches, nr= " << nr << "bytes read"
           << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(&npatches, 1, "int");

    // test
    if (m_myRank == 0 && mVerbose >= 2) {
      printf("Rfile header: magic=%i, prec=%i, att=%i\n", magic, prec, att);
      printf("              azimuth=%e, lon0=%e, lat0=%e\n", alpha, lon0, lat0);
      printf("              pstring-len=%i, pstr='%s'\n", len,
             "not implemented");
      printf("              nblocks=%i\n", npatches);
    }

    // ---------- first part of topography block header
    double hs[3];
    nr = read(fd, hs, 3 * sizeof(double));
    if (nr != 3 * sizeof(double)) {
      cout << rname << " Error reading topography spacings nr= " << nr
           << "bytes read" << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(hs, 3, "double");
    float_sw4 hh, hv, z0;
    hh = hs[0];
    hv = hs[1];
    z0 = hs[2];

    // ---------- second part of topography block header
    int dim[4];
    nr = read(fd, dim, 4 * sizeof(int));
    if (nr != 4 * sizeof(int)) {
      cout << rname << " Error reading topography dimensions nr= " << nr
           << "bytes read" << endl;
      close(fd);
      return;
    }
    if (swapbytes) bswap.byte_rev(dim, 4, "int");

    int nctop, nitop, njtop, nktop;
    nctop = dim[0];
    nitop = dim[1];
    njtop = dim[2];
    nktop = dim[3];
    if (nctop != 1 || nktop != 1) {
      cout << rname << " Error, topography has nc = " << nctop
           << " and nk = " << nktop << endl;
      close(fd);
      return;
    }

    // test
    if (m_myRank == 0 && mVerbose >= 2) {
      printf("Topography header (block #1)\n");
      printf("  hh=%e, hv=%e, z0=%e\n", hh, hv, z0);
      printf("  nc=%i, ni=%i, nj=%i, nk=%i\n", nctop, nitop, njtop, nktop);
    }

    // ---------- Skip other block headers
    for (int p = 0; p < npatches - 1; p++)
      nr = lseek(fd, 4 * sizeof(int) + 3 * sizeof(double), SEEK_CUR);

    // ---------- read topography on file into array gridElev
    bool roworder = true;
    gridElev.define(1, nitop, 1, njtop, 1, 1);

    if (prec == 8) {
      double* data = new double[nitop * njtop];
      nr = read(fd, data, (static_cast<size_t>(nitop)) * njtop * flsize);
      if (nr != (static_cast<size_t>(nitop)) * njtop * flsize) {
        cout << rname << " Error reading topography, nr = " << nr
             << " but need " << (static_cast<size_t>(nitop)) * njtop * flsize
             << "bytes " << endl;
        close(fd);
        return;
      }
      if (swapbytes) bswap.byte_rev(data, nitop * njtop, "double");

      gridElev.assign(data);
      delete[] data;
    } else {
      float* data = new float[nitop * njtop];
      nr = read(fd, data, (static_cast<size_t>(nitop)) * njtop * flsize);
      if (nr != (static_cast<size_t>(nitop)) * njtop * flsize) {
        cout << rname << " Error reading topography, nr = " << nr
             << " but need " << (static_cast<size_t>(nitop)) * njtop * flsize
             << "bytes " << endl;
        close(fd);
        return;
      }
      if (swapbytes) bswap.byte_rev(data, nitop * njtop, "float");

      gridElev.assign(data);

      // test
      if (m_myRank == 0 && mVerbose >= 3) {
        printf("1st topo (float) data=%e, gridElev(1,1,1)=%e\n", data[0],
               gridElev(1, 1, 1));
        printf("last topo (float) data=%e, gridElev(ni,nj,1)=%e\n",
               data[nitop * njtop - 1], gridElev(nitop, njtop, 1));
        // get min and max
        float tmax = -9e-10, tmin = 9e+10;
        for (int q = 0; q < nitop * njtop; q++) {
          if (data[q] > tmax) tmax = data[q];
          if (data[q] < tmin) tmin = data[q];
        }
        printf("topo max (float)=%e, min (float)=%e\n", tmax, tmin);
      }
      delete[] data;
    }
    if (roworder) gridElev.transposeik();

    // ---------- done reading
    close(fd);

    double x0, y0;  // Origin on grid file
    computeCartesianCoord(x0, y0, lon0, lat0);

    // test
    if (m_myRank == 0 && mVerbose >= 3) {
      printf("mat-lon0=%e mat-lat0=%e, comp-x0=%e, commp-y0=%e\n", lon0, lat0,
             x0, y0);
    }

    // Topography read, next interpolate to the computational grid
    int topLevel = mNumberOfGrids - 1;

    float_sw4 topomax = -1e99, topomin = 1e99;
#pragma omp parallel for reduction(max : topomax) reduction(min : topomin)
    for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i) {
      for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j) {
        float_sw4 x = (i - 1) * mGridSize[topLevel];
        float_sw4 y = (j - 1) * mGridSize[topLevel];
        int i0 = static_cast<int>(trunc(1 + (x - x0) / hh));
        int j0 = static_cast<int>(trunc(1 + (y - y0) / hh));
        // test
        float_sw4 xmat0 = (i0 - 1) * hh, ymat0 = (j0 - 1) * hh;
        float_sw4 xmatx = x - x0, ymaty = y - y0;

        if (mVerbose >= 3) {
          if (xmatx < xmat0 || xmatx > xmat0 + hh)
            printf("WARNING: i0=%i is out of bounds for x=%e, xmatx=%e\n", i0,
                   x, xmatx);
          if (ymaty < ymat0 || ymaty > ymat0 + hh)
            printf("WARNING: i0=%i is out of bounds for y=%e, ymaty=%e\n", i0,
                   y, ymaty);
        }

        // end test

        bool extrapol = false;
        if (i0 < -1) {
          extrapol = true;
          i0 = 1;
        } else if (i0 < 2)
          i0 = 2;

        if (i0 > nitop + 1) {
          extrapol = true;
          i0 = nitop;
        } else if (i0 > nitop - 2)
          i0 = nitop - 2;

        if (j0 < -1) {
          extrapol = true;
          j0 = 1;
        } else if (j0 < 2)
          j0 = 2;

        if (j0 > njtop + 1) {
          extrapol = true;
          j0 = njtop;
        } else if (j0 > njtop - 2)
          j0 = njtop - 2;

        if (!extrapol) {
          float_sw4 q = (x - x0 - (i0 - 1) * hh) / hh;
          float_sw4 r = (y - y0 - (j0 - 1) * hh) / hh;
          float_sw4 Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1,
              tjp2;
          Qim1 = (q) * (q - 1) * (q - 2) / (-6.);
          Qi = (q + 1) * (q - 1) * (q - 2) / (2.);
          Qip1 = (q + 1) * (q) * (q - 2) / (-2.);
          Qip2 = (q + 1) * (q) * (q - 1) / (6.);

          Rjm1 = (r) * (r - 1) * (r - 2) / (-6.);
          Rj = (r + 1) * (r - 1) * (r - 2) / (2.);
          Rjp1 = (r + 1) * (r) * (r - 2) / (-2.);
          Rjp2 = (r + 1) * (r) * (r - 1) / (6.);

          // test
          if (mVerbose >= 3) {
            if (i0 < 2 || i0 > nitop - 2)
              printf("WARNING: topo interp out of bounds i0=%i, nitop=%i\n", i0,
                     nitop);
            if (j0 < 2 || j0 > njtop - 2)
              printf("WARNING: topo interp out of bounds j0=%i, njtop=%i\n", j0,
                     njtop);
          }

          tjm1 = Qim1 * gridElev(i0 - 1, j0 - 1, 1) +
                 Qi * gridElev(i0, j0 - 1, 1) +
                 Qip1 * gridElev(i0 + 1, j0 - 1, 1) +
                 Qip2 * gridElev(i0 + 2, j0 - 1, 1);
          tj = Qim1 * gridElev(i0 - 1, j0, 1) + Qi * gridElev(i0, j0, 1) +
               Qip1 * gridElev(i0 + 1, j0, 1) + Qip2 * gridElev(i0 + 2, j0, 1);
          tjp1 = Qim1 * gridElev(i0 - 1, j0 + 1, 1) +
                 Qi * gridElev(i0, j0 + 1, 1) +
                 Qip1 * gridElev(i0 + 1, j0 + 1, 1) +
                 Qip2 * gridElev(i0 + 2, j0 + 1, 1);
          tjp2 = Qim1 * gridElev(i0 - 1, j0 + 2, 1) +
                 Qi * gridElev(i0, j0 + 2, 1) +
                 Qip1 * gridElev(i0 + 1, j0 + 2, 1) +
                 Qip2 * gridElev(i0 + 2, j0 + 2, 1);
          mTopo(i, j, 1) = Rjm1 * tjm1 + Rj * tj + Rjp1 * tjp1 + Rjp2 * tjp2;
        } else {
          // tmp
          if (mVerbose >= 3) {
            printf(
                "INFO: topo extrapolated for i=%i, j=%i, x=%e, y=%e, i0=%i, "
                "j0=%i\n",
                i, j, x, y, i0, j0);
          }

          mTopo(i, j, 1) = gridElev(i0, j0, 1);
        }

        // test
        if (mTopo(i, j, 1) > topomax) topomax = mTopo(i, j, 1);
        if (mTopo(i, j, 1) < topomin) topomin = mTopo(i, j, 1);

      }  // end for j
    }    // end for i

    // test
    if (m_myRank == 0 && mVerbose >= 3) {
      printf("Topo variation on comp grid: max=%e min=%e\n", topomax, topomin);
    }

  } else
    cout << rname << " error could not open file " << a_topoFileName << endl;
}

//-----------------------------------------------------------------------
bool EW::is_onesided(int g, int side) const { return m_onesided[g][side] == 1; }

//-----------------------------------------------------------------------
void EW::get_gridgen_info(int& order, float_sw4& zetaBreak) const {
  order = m_grid_interpolation_order;
  zetaBreak = m_zetaBreak;
}

//-----------------------------------------------------------------------
void EW::get_nr_of_material_parameters(int& nmvar) {
  nmvar = 0;
  for (int g = 0; g < mNumberOfGrids; g++) {
    if (m_iEndAct[g] - m_iStartAct[g] + 1 > 0 &&
        m_jEndAct[g] - m_jStartAct[g] + 1 > 0 &&
        m_kEndAct[g] - m_kStartAct[g] + 1 > 0)
      nmvar += (m_iEndAct[g] - m_iStartAct[g] + 1) *
               (m_jEndAct[g] - m_jStartAct[g] + 1) *
               (m_kEndAct[g] - m_kStartAct[g] + 1) * 3;
  }
}

//-----------------------------------------------------------------------
void EW::parameters_to_material(int nmpar, float_sw4* xm, vector<Sarray>& rho,
                                vector<Sarray>& mu, vector<Sarray>& lambda) {
  size_t gp, ind = 0;
  for (int g = 0; g < mNumberOfGrids; g++) {
    rho[g].copy(mRho[g]);
    mu[g].copy(mMu[g]);
    lambda[g].copy(mLambda[g]);
    if (g == 0)
      gp = 0;
    else
      gp = gp + 3 * ind;
    ind = 0;
    for (int k = m_kStartAct[g]; k <= m_kEndAct[g]; k++)
      for (int j = m_jStartAct[g]; j <= m_jEndAct[g]; j++)
        for (int i = m_iStartAct[g]; i <= m_iEndAct[g]; i++) {
          rho[g](i, j, k) = xm[gp + ind * 3];
          mu[g](i, j, k) = xm[gp + ind * 3 + 1];
          lambda[g](i, j, k) = xm[gp + ind * 3 + 2];
          ind++;
        }
    // update stored material
    mRho[g].copy(rho[g]);
    mMu[g].copy(mu[g]);
    mLambda[g].copy(lambda[g]);
  }
}

//-----------------------------------------------------------------------
void EW::material_to_parameters(int nmpar, float_sw4* xm, vector<Sarray>& rho,
                                vector<Sarray>& mu, vector<Sarray>& lambda) {
  size_t gp, ind = 0;
  for (int g = 0; g < mNumberOfGrids; g++) {
    if (g == 0)
      gp = 0;
    else
      gp = gp + 3 * ind;
    ind = 0;
    for (int k = m_kStartAct[g]; k <= m_kEndAct[g]; k++)
      for (int j = m_jStartAct[g]; j <= m_jEndAct[g]; j++)
        for (int i = m_iStartAct[g]; i <= m_iEndAct[g]; i++) {
          xm[gp + ind * 3] = rho[g](i, j, k);
          xm[gp + ind * 3 + 1] = mu[g](i, j, k);
          xm[gp + ind * 3 + 2] = lambda[g](i, j, k);
          ind++;
        }
  }
}

//-----------------------------------------------------------------------
void EW::get_material_parameter(int nmpar, float_sw4* xm) {
  size_t gp, ind = 0;
  for (int g = 0; g < mNumberOfGrids; g++) {
    if (g == 0)
      gp = 0;
    else
      gp = gp + 3 * ind;
    ind = 0;
    for (int k = m_kStartAct[g]; k <= m_kEndAct[g]; k++)
      for (int j = m_jStartAct[g]; j <= m_jEndAct[g]; j++)
        for (int i = m_iStartAct[g]; i <= m_iEndAct[g]; i++) {
          xm[gp + ind * 3] = mRho[g](i, j, k);
          xm[gp + ind * 3 + 1] = mMu[g](i, j, k);
          xm[gp + ind * 3 + 2] = mLambda[g](i, j, k);
          ind++;
        }
  }
}

//-----------------------------------------------------------------------
void EW::get_scale_factors(int nmpar, float_sw4* sf) {
  size_t gp, ind = 0;
  float_sw4 rhoscale = 2.0;
  float_sw4 muscale = 1.0;
  float_sw4 lambdascale = 5.4e-3;
  for (int g = 0; g < mNumberOfGrids; g++) {
    if (g == 0)
      gp = 0;
    else
      gp = gp + 3 * ind;
    ind = 0;
    for (int k = m_kStartAct[g]; k <= m_kEndAct[g]; k++)
      for (int j = m_jStartAct[g]; j <= m_jEndAct[g]; j++)
        for (int i = m_iStartAct[g]; i <= m_iEndAct[g]; i++) {
          sf[gp + ind * 3] = rhoscale;
          sf[gp + ind * 3 + 1] = muscale;
          sf[gp + ind * 3 + 2] = lambdascale;
          ind++;
        }
  }
}

//-----------------------------------------------------------------------
void EW::add_to_grad(vector<Sarray>& K, vector<Sarray>& Kacc,
                     vector<Sarray>& Um, vector<Sarray>& U, vector<Sarray>& Up,
                     vector<Sarray>& Uacc, vector<Sarray>& gRho,
                     vector<Sarray>& gMu, vector<Sarray>& gLambda) {
  SW4_MARK_FUNCTION;
  for (int g = 0; g < mNumberOfGrids; g++) {
    int ifirst = m_iStart[g];
    int ilast = m_iEnd[g];
    int jfirst = m_jStart[g];
    int jlast = m_jEnd[g];
    int kfirst = m_kStart[g];
    int klast = m_kEnd[g];
    int ifirstact = m_iStartAct[g];
    int ilastact = m_iEndAct[g];
    int jfirstact = m_jStartAct[g];
    int jlastact = m_jEndAct[g];
    int kfirstact = m_kStartAct[g];
    int klastact = m_kEndAct[g];
    float_sw4* k_ptr = K[g].c_ptr();
    float_sw4* ka_ptr = Kacc[g].c_ptr();
    float_sw4* um_ptr = Um[g].c_ptr();
    float_sw4* u_ptr = U[g].c_ptr();
    float_sw4* up_ptr = Up[g].c_ptr();
    float_sw4* ua_ptr = Uacc[g].c_ptr();
    float_sw4* grho_ptr = gRho[g].c_ptr();
    float_sw4* gmu_ptr = gMu[g].c_ptr();
    float_sw4* glambda_ptr = gLambda[g].c_ptr();
    float_sw4 h = mGridSize[g];
    int* onesided_ptr = m_onesided[g];
    int nb = 4, wb = 6;
    //    if (topographyExists() && g == mNumberOfGrids - 1) {
    if (topographyExists() && g >= mNumberOfCartesianGrids) {
      if (m_croutines) {
        addgradrhoc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, ifirstact,
                       ilastact, jfirstact, jlastact, kfirstact, klastact,
                       k_ptr, ka_ptr, um_ptr, u_ptr, up_ptr, ua_ptr, grho_ptr,
                       mDt, mJ[g].c_ptr(), onesided_ptr);
        addgradmulac_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, ifirstact,
                        ilastact, jfirstact, jlastact, kfirstact, klastact,
                        k_ptr, ka_ptr, u_ptr, ua_ptr, gmu_ptr, glambda_ptr, mDt,
                        h, mMetric[g].c_ptr(), mJ[g].c_ptr(), onesided_ptr, nb,
                        wb, m_bop);
      } else {
        addgradrhoc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                    &ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact,
                    &klastact, k_ptr, ka_ptr, um_ptr, u_ptr, up_ptr, ua_ptr,
                    grho_ptr, &mDt, mJ[g].c_ptr(), onesided_ptr);
        addgradmulac(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                     &ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact,
                     &klastact, k_ptr, ka_ptr, u_ptr, ua_ptr, gmu_ptr,
                     glambda_ptr, &mDt, &h, mMetric[g].c_ptr(), mJ[g].c_ptr(),
                     onesided_ptr, &nb, &wb, m_bop);
      }
    } else {
      if (m_croutines) {
        addgradrho_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, ifirstact,
                      ilastact, jfirstact, jlastact, kfirstact, klastact, k_ptr,
                      ka_ptr, um_ptr, u_ptr, up_ptr, ua_ptr, grho_ptr, mDt, h,
                      onesided_ptr);
        addgradmula_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, ifirstact,
                       ilastact, jfirstact, jlastact, kfirstact, klastact,
                       k_ptr, ka_ptr, u_ptr, ua_ptr, gmu_ptr, glambda_ptr, mDt,
                       h, onesided_ptr, nb, wb, m_bop);
      } else {
        addgradrho(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                   &ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact,
                   &klastact, k_ptr, ka_ptr, um_ptr, u_ptr, up_ptr, ua_ptr,
                   grho_ptr, &mDt, &h, onesided_ptr);
        addgradmula(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                    &ifirstact, &ilastact, &jfirstact, &jlastact, &kfirstact,
                    &klastact, k_ptr, ka_ptr, u_ptr, ua_ptr, gmu_ptr,
                    glambda_ptr, &mDt, &h, onesided_ptr, &nb, &wb, m_bop);
      }
    }
  }
}

//-----------------------------------------------------------------------
void EW::perturb_mtrl() {
  int g = mNumberOfGrids - 1;
  if (m_perturb != 0 && point_in_proc(m_iperturb, m_jperturb, g)) {
    cout << "per = " << m_perturb << " " << m_iperturb << " " << m_jperturb
         << " " << m_kperturb << endl;
    if (m_iperturb < m_iStartAct[g] || m_iperturb > m_iEndAct[g])
      cout << "warning i-index outside active domain " << endl;
    if (m_jperturb < m_jStartAct[g] || m_jperturb > m_jEndAct[g])
      cout << "warning j-index outside active domain " << endl;
    if (m_kperturb < m_kStartAct[g] || m_kperturb > m_kEndAct[g])
      cout << "warning k-index outside active domain " << endl;
    if (m_pervar == 1)
      mMu[g](m_iperturb, m_jperturb, m_kperturb) += m_perturb;
    else if (m_pervar == 2)
      mLambda[g](m_iperturb, m_jperturb, m_kperturb) += m_perturb;
    else if (m_pervar == 0)
      mRho[g](m_iperturb, m_jperturb, m_kperturb) += m_perturb;
  }
}

//-----------------------------------------------------------------------
void EW::perturb_mtrl(int peri, int perj, int perk, float_sw4 h, int grid,
                      int var) {
  if (h != 0 && point_in_proc(peri, perj, grid)) {
    //      cout << "per = " << m_perturb << " " << m_iperturb << " " <<
    //      m_jperturb << " " << m_kperturb << endl;
    if (peri < m_iStartAct[grid] || peri > m_iEndAct[grid])
      cout << "warning i-index outside active domain " << endl;
    if (perj < m_jStartAct[grid] || perj > m_jEndAct[grid])
      cout << "warning j-index outside active domain " << endl;
    if (perk < m_kStartAct[grid] || perk > m_kEndAct[grid])
      cout << "warning k-index outside active domain " << endl;
    if (var == 0)
      mRho[grid](peri, perj, perk) += h;
    else if (var == 1)
      mMu[grid](peri, perj, perk) += h;
    else if (var == 2)
      mLambda[grid](peri, perj, perk) += h;
  }
}

//-----------------------------------------------------------------------
void EW::get_optmethod(int& method, int& bfgs_m) {
  method = m_opt_method;
  bfgs_m = m_lbfgs_m;
}

//-----------------------------------------------------------------------
void EW::set_epicenter(float_sw4 epiLat, float_sw4 epiLon, float_sw4 epiDepth,
                       float_sw4 earliestTime, int e) {
  m_epi_lat[e] = epiLat;
  m_epi_lon[e] = epiLon;
  m_epi_depth[e] = epiDepth;
  m_epi_t0[e] = earliestTime;
}

//-----------------------------------------------------------------------
void EW::get_epicenter(float_sw4& epiLat, float_sw4& epiLon,
                       float_sw4& epiDepth, float_sw4& earliestTime, int e) {
  epiLat = m_epi_lat[e];
  epiLon = m_epi_lon[e];
  epiDepth = m_epi_depth[e];
  earliestTime = m_epi_t0[e];
}

//-----------------------------------------------------------------------
bool EW::check_for_nan(vector<Sarray>& a_U, int verbose, string name) {
  SW4_MARK_FUNCTION;
  bool retval = false;
  for (int g = 0; g < mNumberOfGrids; g++) {
    size_t nn = a_U[g].count_nans();
    retval = retval || nn > 0;
    if (nn > 0 && verbose == 1) {
      int cnan, inan, jnan, knan;
      a_U[g].count_nans(cnan, inan, jnan, knan);
      cout << "grid " << g << " array " << name << " found " << nn
           << "  nans. First nan at " << cnan << " " << inan << " " << jnan
           << " " << knan << endl;
    }
  }
  return retval;
}

//-----------------------------------------------------------------------
void EW::check_min_max_int(vector<Sarray>& a_U) {
  int nc = a_U[0].m_nc;
  float_sw4* mx = new float_sw4[nc];
  float_sw4* mn = new float_sw4[nc];
  for (int g = 0; g < a_U.size(); g++) {
    //      double mx[4]={-1e30,-1e30,-1e30,-1e30};
    //      double mn[4]={1e30,1e30,1e30,1e30};
    for (int c = 1; c <= nc; c++) {
      //	 mx[c] = -1e38;
      //	 mn[c] =  1e38;
      mx[c - 1] = mn[c - 1] =
          a_U[g](c, m_iStartInt[g], m_jStartInt[g], m_kStartInt[g]);
    }

    for (int k = m_kStartInt[g]; k <= m_kEndInt[g]; k++)
      for (int j = m_jStartInt[g]; j <= m_jEndInt[g]; j++)
        for (int i = m_iStartInt[g]; i <= m_iEndInt[g]; i++)
          for (int c = 1; c <= nc; c++) {
            if (mx[c - 1] < a_U[g](c, i, j, k)) mx[c - 1] = a_U[g](c, i, j, k);
            if (mn[c - 1] > a_U[g](c, i, j, k)) mn[c - 1] = a_U[g](c, i, j, k);
          }
    cout << g << " " << mn[0] << " " << mx[0] << endl;
  }
  delete[] mx;
  delete[] mn;
}

#ifdef ENABLE_OPT
//-----------------------------------------------------------------------
void EW::material_correction(int nmpar, float_sw4* xm)
// routine to enforce material speed limits and positive density
{
  SW4_MARK_FUNCTION;
  float_sw4 vsmin = -1;
  if (m_useVelocityThresholds) vsmin = m_vsMin;
  float_sw4 rhoscale = 1, muscale = 1, lascale = 1;

  parameters_to_material(nmpar, xm, mRho, mMu, mLambda);
  for (int g = 0; g < mNumberOfGrids; g++) {
    int info;
    int ifirst = m_iStart[g];
    int ilast = m_iEnd[g];
    int jfirst = m_jStart[g];
    int jlast = m_jEnd[g];
    int kfirst = m_kStart[g];
    int klast = m_kEnd[g];
    int ifirstact = m_iStartAct[g];
    int ilastact = m_iEndAct[g];
    int jfirstact = m_jStartAct[g];
    int jlastact = m_jEndAct[g];
    int kfirstact = m_kStartAct[g];
    int klastact = m_kEndAct[g];

    float_sw4* rhop = mRho[g].c_ptr();
    float_sw4* mup = mMu[g].c_ptr();
    float_sw4* lap = mLambda[g].c_ptr();

    if (topographyExists() && g == mNumberOfGrids - 1) {
      // Curvilinear
      F77_FUNC(projectmtrlc, PROJECTMTRLC)
      (&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &ifirstact, &ilastact,
       &jfirstact, &jlastact, &kfirstact, &klastact, rhop, mup, lap, &mDt,
       mMetric.c_ptr(), mJ[g].c_ptr(), &mCFLmax, &vsmin, &rhoscale, &muscale,
       &lascale, &info);
    } else {
      // Cartesian
      F77_FUNC(projectmtrl, PROJECTMTRL)
      (&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &ifirstact, &ilastact,
       &jfirstact, &jlastact, &kfirstact, &klastact, rhop, mup, lap, &mDt,
       &mGridSize[g], &mCFLmax, &vsmin, &rhoscale, &muscale, &lascale, &info);
    }
    if (info != 0)
      cout << "Grid " << g << " info = " << info << " from projectmtrl" << endl;
  }
  material_to_parameters(nmpar, xm, mRho, mMu, mLambda);
}
#endif

#ifdef ENABLE_OPT
//-----------------------------------------------------------------------
void EW::project_material(vector<Sarray>& a_rho, vector<Sarray>& a_mu,
                          vector<Sarray>& a_lambda, int& info)
// routine to enforce material speed limits and positive density
{
  SW4_MARK_FUNCTION;
  float_sw4 vsmin = -1;
  if (m_useVelocityThresholds) vsmin = m_vsMin;
  float_sw4 rhoscale = 1, muscale = 1, lascale = 1;
  info = 0;
  for (int g = 0; g < mNumberOfGrids; g++) {
    int infogrid;
    int ifirst = m_iStart[g];
    int ilast = m_iEnd[g];
    int jfirst = m_jStart[g];
    int jlast = m_jEnd[g];
    int kfirst = m_kStart[g];
    int klast = m_kEnd[g];
    int ifirstact = m_iStartAct[g];
    int ilastact = m_iEndAct[g];
    int jfirstact = m_jStartAct[g];
    int jlastact = m_jEndAct[g];
    int kfirstact = m_kStartAct[g];
    int klastact = m_kEndAct[g];

    float_sw4* rhop = a_rho[g].c_ptr();
    float_sw4* mup = a_mu[g].c_ptr();
    float_sw4* lap = a_lambda[g].c_ptr();

    if (topographyExists() && g == mNumberOfGrids - 1) {
      // Curvilinear
      F77_FUNC(projectmtrlc, PROJECTMTRLC)
      (&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &ifirstact, &ilastact,
       &jfirstact, &jlastact, &kfirstact, &klastact, rhop, mup, lap, &mDt,
       mMetric.c_ptr(), mJ[g].c_ptr(), &mCFLmax, &vsmin, &rhoscale, &muscale,
       &lascale, &infogrid);
    } else {
      // Cartesian
      F77_FUNC(projectmtrl, PROJECTMTRL)
      (&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, &ifirstact, &ilastact,
       &jfirstact, &jlastact, &kfirstact, &klastact, rhop, mup, lap, &mDt,
       &mGridSize[g], &mCFLmax, &vsmin, &rhoscale, &muscale, &lascale,
       &infogrid);
    }
    if (infogrid != 0) {
      cout << "Grid " << g << " info = " << infogrid << " from projectmtrl"
           << endl;
      if (info == 0) info = infogrid;
    }
  }
}
#endif

#ifdef ENABLE_OPT
//-----------------------------------------------------------------------
void EW::check_material(vector<Sarray>& a_rho, vector<Sarray>& a_mu,
                        vector<Sarray>& a_lambda, int& ok) {
  ok = 1;
  for (int g = 0; g < mNumberOfGrids; g++) {
    int infogrid;
    int ifirst = m_iStart[g];
    int ilast = m_iEnd[g];
    int jfirst = m_jStart[g];
    int jlast = m_jEnd[g];
    int kfirst = m_kStart[g];
    int klast = m_kEnd[g];

    float_sw4 limits[10];

    float_sw4* rhop = a_rho[g].c_ptr();
    float_sw4* mup = a_mu[g].c_ptr();
    float_sw4* lap = a_lambda[g].c_ptr();

    //      if( topographyExists() && g == mNumberOfGrids-1 )
    //      {
    //	 // Curvilinear
    //	 F77_FUNC(projectmtrlc,PROJECTMTRLC)( &ifirst, &ilast, &jfirst, &jlast,
    //&kfirst, &klast, 					    &ifirstact,
    //&ilastact, &jfirstact, &jlastact, &kfirstact, &klastact,  rhop, mup, lap,
    //&mDt, mMetric.c_ptr(),
    // mJ[g].c_ptr(), 					      &mCFLmax, &vsmin,
    // &rhoscale, &muscale, &lascale, &infogrid );
    //      }
    //      else
    //      {
    // Cartesian
    F77_FUNC(checkmtrl, CHECKMTRL)
    (&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, rhop, mup, lap, &mDt,
     &mGridSize[g], limits);
    float_sw4 local[5] = {limits[0], limits[2], limits[4], limits[7],
                          limits[8]};
    float_sw4 global[5];
    MPI_Allreduce(local, global, 5, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    limits[0] = global[0];
    limits[2] = global[1];
    limits[4] = global[2];
    limits[7] = global[3];
    limits[8] = global[4];
    local[0] = limits[1];
    local[1] = limits[3];
    local[2] = limits[5];
    local[3] = limits[6];
    local[4] = limits[9];
    MPI_Allreduce(local, global, 5, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    limits[1] = global[0];
    limits[3] = global[1];
    limits[5] = global[2];
    limits[6] = global[3];
    limits[9] = global[4];
    if (proc_zero()) {
      cout << limits[0] << " <=   rho    <= " << limits[1] << " (grid " << g
           << ")" << endl;
      cout << limits[2] << " <=    mu    <= " << limits[3] << " (grid " << g
           << ")" << endl;
      cout << limits[4] << " <=  lambda  <= " << limits[5] << " (grid " << g
           << ")" << endl;

      if (limits[0] < 0)
        cout << "rho_min = " << limits[0] << " on grid " << g << endl;
      if (limits[2] < 0)
        cout << "mu_min = " << limits[2] << " on grid " << g << endl;
      if (limits[4] < 0)
        cout << "lambda_min = " << limits[4] << " on grid " << g << endl;
      if (limits[6] < 0)
        cout << " cfl_max  is imaginary on grid " << g << endl;
      else
        cout << " cfl_max = " << sqrt(limits[6]) << " on grid " << g << endl;
    }
    ok = ok && (limits[0] > 0 && limits[2] > 0 &&
                limits[6] < mCFLmax * mCFLmax && limits[8] > 0);
  }
}
#endif

//-----------------------------------------------------------------------
void EW::extrapolateTopo(Sarray& field) {
  int k = 1;                   // field is assumed to be a 2-D array
  int g = mNumberOfGrids - 1;  // top grid
  int nExtrap = 0;

  if (m_iStartInt[g] == 1) {
    for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
      for (int i = m_iStart[g]; i < 1; i++)
        if (field(i, j, k) == NO_TOPO) {
          field(i, j, k) = field(1, j, k);
          nExtrap += 1;
        }
  }

  if (m_iEndInt[g] == m_global_nx[g]) {
    for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
      for (int i = m_iEndInt[g] + 1; i <= m_iEnd[g]; i++)
        if (field(i, j, k) == NO_TOPO) {
          field(i, j, k) = field(m_iEndInt[g], j, k);
          nExtrap += 1;
        }
  }

  if (m_jStartInt[g] == 1) {
    for (int j = m_jStart[g]; j < 1; j++)
      for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
        if (field(i, j, k) == NO_TOPO) {
          field(i, j, k) = field(i, 1, k);
          nExtrap += 1;
        }
  }

  if (m_jEndInt[g] == m_global_ny[g]) {
    for (int j = m_jEndInt[g] + 1; j <= m_jEnd[g]; j++)
      for (int i = m_iStart[g]; i <= m_iEnd[g]; i++)
        if (field(i, j, k) == NO_TOPO) {
          field(i, j, k) = field(i, m_jEndInt[g], k);
          nExtrap += 1;
        }
  }
  int nExtrapGlobal = 0;
  MPI_Allreduce(&nExtrap, &nExtrapGlobal, 1, MPI_INT, MPI_SUM,
                m_cartesian_communicator);

  if (nExtrapGlobal > 0 && proc_zero())
    printf("*** extrapolated topography to %i ghost points\n", nExtrapGlobal);
}

//-----------------------------------------------------------------------
void EW::checkTopo(Sarray& field) {
  bool topo_ok = true;
  int k = 1;                   // field is a 2-D array
  int g = mNumberOfGrids - 1;  // top grid

  for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
    for (int i = m_iStart[g]; i <= m_iEnd[g]; i++) {
      if (field(i, j, k) == NO_TOPO) {
        // print some msg is verbose is high enough?
        topo_ok = false;
      }
    }

  CHECK_INPUT(topo_ok, "There are undefined values in the topography array")
}

//-----------------------------------------------------------------------
void EW::setup_attenuation_relaxation(float_sw4 minvsoh) {
  // using minvsoh, estimate highest resolved frequency on the mesh. Use this to
  // assign omega_max for attenuation use the number of mechanisms to determine
  // bandwidth of attenuation model. Use this to assign omega_min

  if (m_att_use_max_frequency)
    m_max_omega = 2 * M_PI * m_att_max_frequency;
  else
    m_max_omega = 2. * M_PI * minvsoh / m_att_ppw;

  // the band width is set to get approximately constant Q throughout the
  // frequency band
  //     if (n <= 2)
  //     {
  //       m_min_omega = m_max_omega/10.; // decent accuracy for 2 mechanisms
  //     }
  //     else if (n == 3)
  //     {
  //       m_min_omega = m_max_omega/80.; // decent accuracy for 3 mechanisms
  //     }
  //     else if (n == 4)
  //     {
  //       m_min_omega = m_max_omega/150.; // decent accuracy for 4 mechanisms
  //     }
  //     else if (n >= 5)
  //     {
  //       m_min_omega = m_max_omega/2000.; // decent accuracy for 5 mechanisms
  //     }

  // always use a frequency band that is 2 decades wide
  m_min_omega = m_max_omega / 100.;

  if (proc_zero()) {
    printf(
        "\n*** Attenuation parameters calculated for %i mechanisms,\n"
        "      max freq=%e [Hz], min_freq=%e [Hz], velo_freq=%e [Hz]\n\n",
        m_number_mechanisms, m_max_omega / 2 / M_PI, m_min_omega / 2 / M_PI,
        m_velo_omega / 2 / M_PI);
  }
  int n = m_number_mechanisms;
  if (n == 1) {
    mOmegaVE[0] = m_max_omega;
  } else if (n > 1) {
    float_sw4 r = pow(m_max_omega / m_min_omega, 1.0 / (n - 1));
    mOmegaVE[0] = m_min_omega;
    mOmegaVE[n - 1] = m_max_omega;
    for (int k = 1; k <= n - 2; k++) mOmegaVE[k] = m_min_omega * pow(r, k);
  }
}

//-----------------------------------------------------------------------
void EW::setup_viscoelastic() {
  SW4_MARK_FUNCTION;
  int k, g;

  // number of collocation points
  int n = m_number_mechanisms;
  int nc = 2 * n - 1;

  if (n > 0) {
    // collocation frequencies
    vector<float_sw4> omc(nc);
    if (n > 1) {
      float_sw4 r = pow(m_max_omega / m_min_omega, 1.0 / (n - 1));
      omc[0] = mOmegaVE[0];
      for (int k = 0; k <= 2 * n - 2; k++)
        omc[k] = m_min_omega * pow(r, 0.5 * k);
    } else
      omc[0] = mOmegaVE[0];

    // tmp: print omega and omc
    if (proc_zero() && mVerbose >= 1) {
      for (k = 0; k < n; k++) printf("omega[%i]=%e ", k, mOmegaVE[k]);
      printf("\n");
      for (k = 0; k < nc; k++) printf("omc[%i]=%e", k, omc[k]);
      printf("\n\n");
    }

// setup least squares problem (matrix and rhs depends on Qs & Qp)

// test for q=80
//       float_sw4 q0=80.0, qs, qp, mu_tmp, lambda_tmp, kappa_tmp, mu_0,
//       lambda_0, kappa_0, imm, rem, mmag, bsum;

// use base 0 indexing of matrix
#define a(i, j) a_[i + j * nc]

#ifndef _OPENMP
    double* a_ = new float_sw4[n * nc];
    double* beta = new double[nc];
    double* gamma = new double[nc];
    int lwork = 3 * n;
    double* work = new double[lwork];
#endif
    //g=0;
    //std::cout<<"mMug "<<as_int(mMu[g].space)<<"\n";
    //std::cout<<"mLamnda "<<as_int(mLambda[g].space)<<"\n";
    //std::cout<<"mMuVE[g][0] "<<as_int(mMuVE[g][0].space)<<"\n";
    //std::cout<<"mLambdaVE[g][0] "<<as_int(mLambdaVE[g][0].space)<<"\n";
    // loop over all grid points in all grids
    for (g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++)
        for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
          for (int i = m_iStart[g]; i <= m_iEnd[g]; i++) {
            // Code for locating libsci errors on RZNevada April 21 2021
            // int ijk=i+j+k;
            // bool
            // doit=(ijk==(m_iStart[g]+m_jStart[g]+m_kStart[g]))&&(!getRank());
#ifdef _OPENMP
            double* a_ = new float_sw4[n * nc];
            double* beta = new double[nc];
            double* gamma = new double[nc];
            int lwork = 3 * n;
            double* work = new double[lwork];
#endif
            char trans = 'N';
            int info = 0, nrhs = 1, lda = nc, ldb = nc;

            float_sw4 mu_tmp = mMu[g](i, j, k);
            float_sw4 lambda_tmp = mLambda[g](i, j, k);
            float_sw4 kappa_tmp = lambda_tmp + 2 * mu_tmp;
            float_sw4 qs = mQs[g](i, j, k);
            float_sw4 qp = mQp[g](i, j, k);

            //
            // qs gives beta coefficients
            //
            for (int q = 0; q < nc; q++) {
              beta[q] = 1. / qs;
              /// if (doit)std::cout<<"M["<<q<<","<<qs<<"]=";
              for (int nu = 0; nu < n; nu++) {
                double tmp = a(q, nu) =
                    (omc[q] * mOmegaVE[nu] + SQR(mOmegaVE[nu]) / qs) /
                    (SQR(mOmegaVE[nu]) + SQR(omc[q]));
                // if (doit)std::cout<<a(q,nu)<<","<<mOmegaVE[nu]<<":";
              }
              // if (doit)std::cout<<"\n";
            }
            // solve the system in least squares sense
            F77_FUNC(dgels, DGELS)
            (trans, nc, n, nrhs, a_, lda, beta, ldb, work, lwork, info);
            if (info != 0) {
              printf(
                  "setup_viscoelastic:: solving for qs=%e, processor=%i, dgels "
                  "returned error code = %i\n",
                  qs, m_myRank, info);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
            // check that sum(beta) < 1
            float_sw4 bsum = 0.;
            for (int nu = 0; nu < n; nu++) bsum += beta[nu];
            if (bsum >= 1.) {
              printf(
                  "setup_viscoelastic:: sum(beta)=%e >= 1 for g=%i, i=%i, "
                  "j=%i, k=%i\n",
                  bsum, g, i, j, k);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // calculate unrelaxed mu_0
            float_sw4 rem = 0., imm = 0.;
            for (int nu = 0; nu < n; nu++) {
              rem += beta[nu] * SQR(mOmegaVE[nu]) /
                     (SQR(mOmegaVE[nu]) + SQR(m_velo_omega));
              imm += beta[nu] * mOmegaVE[nu] * m_velo_omega /
                     (SQR(mOmegaVE[nu]) + SQR(m_velo_omega));
            }
            rem = 1 - rem;
            float_sw4 mmag = sqrt(SQR(rem) + SQR(imm));
            // should also divide by cos^2(delta/2), where delta is the
            // loss-angle, but this makes minimal difference for Q>25
            float_sw4 mu_0 = mu_tmp / mmag;
            // calculate viscoelastic mu:
            for (int nu = 0; nu < n; nu++) {
              mMuVE[g][nu](i, j, k) = mu_0 * beta[nu];
            }
            // save the unrelaxed value
            mMu[g](i, j, k) = mu_0;

            //
            // qp gives gamma coefficients
            //
            for (int q = 0; q < nc; q++) {
              gamma[q] = 1. / qp;
              for (int nu = 0; nu < n; nu++) {
                a(q, nu) = (omc[q] * mOmegaVE[nu] + SQR(mOmegaVE[nu]) / qp) /
                           (SQR(mOmegaVE[nu]) + SQR(omc[q]));
              }
            }

            // solve the system in least squares sense
            F77_FUNC(dgels, DGELS)
            (trans, nc, n, nrhs, a_, lda, gamma, ldb, work, lwork, info);
            if (info != 0) {
              printf(
                  "setup_viscoelastic:: solving for qp=%e, processor=%i, dgels "
                  "returned error code = %i\n",
                  qp, m_myRank, info);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
            // check that sum(gamma) < 1
            bsum = 0.;
            for (int nu = 0; nu < n; nu++) bsum += gamma[nu];
            if (bsum >= 1.) {
              printf(
                  "setup_viscoelastic:: sum(gamma)=%e >= 1 for g=%i, i=%i, "
                  "j=%i, k=%i\n",
                  bsum, g, i, j, k);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // calculate unrelaxed kappa_0
            rem = 0., imm = 0.;
            for (int nu = 0; nu < n; nu++) {
              rem += gamma[nu] * SQR(mOmegaVE[nu]) /
                     (SQR(mOmegaVE[nu]) + SQR(m_velo_omega));
              imm += gamma[nu] * mOmegaVE[nu] * m_velo_omega /
                     (SQR(mOmegaVE[nu]) + SQR(m_velo_omega));
            }
            rem = 1 - rem;
            mmag = sqrt(SQR(rem) + SQR(imm));
            // should also divide by cos^2(delta/2), where delta is the
            // loss-angle, but this makes minimal difference for Q>25
            float_sw4 kappa_0 = kappa_tmp / mmag;
            // calculate viscoelastic lambdaVE = kappaVE - 2*muVE:
            for (int nu = 0; nu < n; nu++) {
              kappa_tmp = kappa_0 * gamma[nu];
              mLambdaVE[g][nu](i, j, k) = kappa_tmp - 2 * mMuVE[g][nu](i, j, k);
            }
            // save the unrelaxed value
            mLambda[g](i, j, k) = kappa_0 - 2 * mu_0;
            //            if( g==1 && k==m_kEnd[g] && m_myRank == 0 )
            //	       cout << i << " " << j << "mlambdave 0 " <<
            // mLambdaVE[g][0](i,j,k) << "Qs = " << mQs[g](i,j,k) << endl;
            // tmp
            //     printf("Q=%e\n", q0);
            //     for (q=0; q<n; q++)
            //       printf("beta[%i]=%e ", q, b[q]);
            //     printf("\n");

#ifdef _OPENMP
            delete[] a_;
            delete[] beta;
            delete[] gamma;
            delete[] work;
#endif
          }  // end for g,k,j,i
#ifndef _OPENMP
    delete[] a_;
    delete[] beta;
    delete[] gamma;
    delete[] work;
#endif
#undef a
  }
}

//-----------------------------------------------------------------------
void EW::reverse_setup_viscoelastic() {
  // number of collocation points
  int n = m_number_mechanisms;
  int nc = 2 * n - 1;

  if (n > 0) {
    // collocation frequencies
    vector<float_sw4> omc(nc);
    if (n > 1) {
      float_sw4 r = pow(m_max_omega / m_min_omega, 1.0 / (n - 1));
      omc[0] = mOmegaVE[0];
      for (int k = 0; k <= 2 * n - 2; k++)
        omc[k] = m_min_omega * pow(r, 0.5 * k);
    } else
      omc[0] = mOmegaVE[0];

// use base 0 indexing of matrix
#define a(i, j) a_[i + j * nc]
    for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++)
        for (int j = m_jStart[g]; j <= m_jEnd[g]; j++)
          for (int i = m_iStart[g]; i <= m_iEnd[g]; i++) {
            float_sw4 mu, lambda;
            double* a_ = new double[n * nc];
            double* beta = new double[nc];
            double* gamma = new double[nc];
            int lwork = 3 * n;
            double* work = new double[lwork];
            char trans = 'N';
            int info = 0, nrhs = 1, lda = nc, ldb = nc;

            float_sw4 mu_tmp = mMu[g](i, j, k);
            float_sw4 lambda_tmp = mLambda[g](i, j, k);
            float_sw4 qs = mQs[g](i, j, k);
            float_sw4 qp = mQp[g](i, j, k);

            //
            // qs gives beta coefficients
            //
            for (int q = 0; q < nc; q++) {
              beta[q] = 1. / qs;
              for (int nu = 0; nu < n; nu++) {
                a(q, nu) = (omc[q] * mOmegaVE[nu] + SQR(mOmegaVE[nu]) / qs) /
                           (SQR(mOmegaVE[nu]) + SQR(omc[q]));
              }
            }
            // solve the system in least squares sense
            F77_FUNC(dgels, DGELS)
            (trans, nc, n, nrhs, a_, lda, beta, ldb, work, lwork, info);
            if (info != 0) {
              printf(
                  "reverse_setup_viscoelastic:: solving for qs=%e, "
                  "processor=%i, dgels returned error code = %i\n",
                  qs, m_myRank, info);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
            // check that sum(beta) < 1
            float_sw4 bsum = 0.;
            for (int nu = 0; nu < n; nu++) bsum += beta[nu];
            if (bsum >= 1.) {
              printf(
                  "reverse_setup_viscoelastic:: sum(beta)=%e >= 1 for g=%i, "
                  "i=%i, j=%i, k=%i\n",
                  bsum, g, i, j, k);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // calculate unrelaxed mu_0
            float_sw4 rem = 0., imm = 0.;
            for (int nu = 0; nu < n; nu++) {
              rem += beta[nu] * SQR(mOmegaVE[nu]) /
                     (SQR(mOmegaVE[nu]) + SQR(m_velo_omega));
              imm += beta[nu] * mOmegaVE[nu] * m_velo_omega /
                     (SQR(mOmegaVE[nu]) + SQR(m_velo_omega));
            }
            rem = 1 - rem;
            float_sw4 mmag = sqrt(SQR(rem) + SQR(imm));
            // should also divide by cos^2(delta/2), where delta is the
            // loss-angle, but this makes minimal difference for Q>25
            /* float_sw4 mu_0 = mu_tmp/mmag; */
            // calculate viscoelastic mu:
            /* for (int nu=0; nu<n; nu++) */
            /* { */
            /*    mMuVE[g][nu](i,j,k) = mu_0 * beta[nu]; */
            /* } */
            // reverse the value
            mu = mu_tmp * mmag;
            mMu[g](i, j, k) = mu;

            //
            // qp gives gamma coefficients
            //
            for (int q = 0; q < nc; q++) {
              gamma[q] = 1. / qp;
              for (int nu = 0; nu < n; nu++) {
                a(q, nu) = (omc[q] * mOmegaVE[nu] + SQR(mOmegaVE[nu]) / qp) /
                           (SQR(mOmegaVE[nu]) + SQR(omc[q]));
              }
            }

            // solve the system in least squares sense
            F77_FUNC(dgels, DGELS)
            (trans, nc, n, nrhs, a_, lda, gamma, ldb, work, lwork, info);
            if (info != 0) {
              printf(
                  "reverse_setup_viscoelastic:: solving for qp=%e, "
                  "processor=%i, dgels returned error code = %i\n",
                  qp, m_myRank, info);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
            // check that sum(gamma) < 1
            bsum = 0.;
            for (int nu = 0; nu < n; nu++) bsum += gamma[nu];
            if (bsum >= 1.) {
              printf(
                  "reverse_setup_viscoelastic:: sum(gamma)=%e >= 1 for g=%i, "
                  "i=%i, j=%i, k=%i\n",
                  bsum, g, i, j, k);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // calculate unrelaxed kappa_0
            rem = 0., imm = 0.;
            for (int nu = 0; nu < n; nu++) {
              rem += gamma[nu] * SQR(mOmegaVE[nu]) /
                     (SQR(mOmegaVE[nu]) + SQR(m_velo_omega));
              imm += gamma[nu] * mOmegaVE[nu] * m_velo_omega /
                     (SQR(mOmegaVE[nu]) + SQR(m_velo_omega));
            }
            rem = 1 - rem;
            mmag = sqrt(SQR(rem) + SQR(imm));
            // should also divide by cos^2(delta/2), where delta is the
            // loss-angle, but this makes minimal difference for Q>25
            /* float_sw4 kappa_tmp = lambda_tmp + 2*mu_tmp; */
            /* float_sw4 kappa_0 = kappa_tmp/mmag; */
            // calculate viscoelastic lambdaVE = kappaVE - 2*muVE:
            /* for (int nu=0; nu<n; nu++) */
            /* { */
            /*    kappa_tmp = kappa_0 * gamma[nu]; */
            /*    mLambdaVE[g][nu](i,j,k) = kappa_tmp - 2*mMuVE[g][nu](i,j,k);
             */
            /* } */
            // save the unrelaxed value
            /* mLambda[g](i,j,k) = kappa_0 - 2*mu_0; */

            lambda = (2 * mu_tmp + lambda_tmp) * mmag - 2 * mu;
            mLambda[g](i, j, k) = lambda;
            /* if (mLambda[g](i,j,k) < 0) { */
            /*     printf("ERROR, lambda < 0, %f\n", mLambda[g](i,j,k) ); */
            /* } */

            delete[] a_;
            delete[] beta;
            delete[] gamma;
            delete[] work;

          }  // end for g,k,j,i
#undef a
  }
}

//-----------------------------------------------------------------------
void EW::setup_viscoelastic_tw() {
  SW4_MARK_FUNCTION;
  // Set twilight testing values for the attenuation material (mu,lambda).
  if (m_number_mechanisms != 1) {
    printf(
        "setup_viscoelastic_tw:: Number of mechanisms must be %i for twilight "
        "testing, input value = %i \n",
        1, m_number_mechanisms);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  float_sw4 *mu_ptr, *la_ptr;
  int ifirst, ilast, jfirst, jlast, kfirst, klast, g;
  float_sw4 h, zmin, omm, phm, ampmu, ampla;
  for (g = 0; g < mNumberOfCartesianGrids; g++) {
    mu_ptr = mMuVE[g][0].c_ptr();
    la_ptr = mLambdaVE[g][0].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    h = mGridSize[g];
    zmin = m_zmin[g];
    omm = m_twilight_forcing->m_momega;
    phm = m_twilight_forcing->m_mphase;
    ampmu = m_twilight_forcing->m_ampmu;
    ampla = m_twilight_forcing->m_amplambda;
    if (m_croutines)
      exactmatfortatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, mu_ptr,
                         la_ptr, omm, phm, ampmu, ampla, h, zmin);
    else
      exactmatfortatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast, mu_ptr,
                      la_ptr, &omm, &phm, &ampmu, &ampla, &h, &zmin);
  }
  //  if (topographyExists()) {
  for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
    // g = mNumberOfGrids - 1;
    mu_ptr = mMuVE[g][0].c_ptr();
    la_ptr = mLambdaVE[g][0].c_ptr();
    ifirst = m_iStart[g];
    ilast = m_iEnd[g];
    jfirst = m_jStart[g];
    jlast = m_jEnd[g];
    kfirst = m_kStart[g];
    klast = m_kEnd[g];
    omm = m_twilight_forcing->m_momega;
    phm = m_twilight_forcing->m_mphase;
    ampmu = m_twilight_forcing->m_ampmu;
    ampla = m_twilight_forcing->m_amplambda;
    float_sw4* x_ptr = mX[g].c_ptr();
    float_sw4* y_ptr = mY[g].c_ptr();
    float_sw4* z_ptr = mZ[g].c_ptr();
    if (m_croutines)
      exactmatfortattc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, mu_ptr,
                          la_ptr, omm, phm, ampmu, ampla, x_ptr, y_ptr, z_ptr);
    else
      exactmatfortattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                       mu_ptr, la_ptr, &omm, &phm, &ampmu, &ampla, x_ptr, y_ptr,
                       z_ptr);
  }
}

//-----------------------------------------------------------------------
void EW::compute_minvsoverh(float_sw4& minvsoh) {
  SW4_MARK_FUNCTION;
  float_sw4 minvsohloc =
      1.e27;  // what is a good 'large' value to initialize with
  // treat all grids the same, i.e. ignore the effects of variations in the
  // curvilinear grid size
  for (int g = 0; g < mNumberOfGrids; g++) {
    float_sw4* mu = mMu[g].c_ptr();
    float_sw4* rho = mRho[g].c_ptr();
    size_t npts = mMu[g].npts();
    float_sw4 minvs = mu[0] / rho[0];
#pragma omp parallel for reduction(min : minvs)
    for (int i = 0; i < npts; i++) {
      if (mu[i] < minvs * rho[i]) minvs = mu[i] / rho[i];
    }
    minvs = sqrt(minvs);
    minvsohloc = minvs / mGridSize[g];

    // get the global min for this grid
    MPI_Allreduce(&minvsohloc, &mMinVsOverH[g], 1, m_mpifloat, MPI_MIN,
                  m_cartesian_communicator);
  }  // end for all grids
     // min mMinVsOverH is saved in minvsoh
  minvsoh = mMinVsOverH[0];
  for (int g = 1; g < mNumberOfGrids; g++) {
    if (mMinVsOverH[g] < minvsoh) minvsoh = mMinVsOverH[g];
  }
}

//-----------------------------------------------------------------------
bool less_than(GridPointSource* ptsrc1, GridPointSource* ptsrc2) {
  return ptsrc1->m_key < ptsrc2->m_key;
}

//-----------------------------------------------------------------------
void EW::sort_grid_point_sources(vector<GridPointSource*>& point_sources,
                                 vector<int>& identsources) {
  size_t* gptr = new size_t[mNumberOfGrids];
  gptr[0] = 0;
  for (int g = 0; g < mNumberOfGrids - 1; g++) {
    gptr[g + 1] = gptr[g] + static_cast<size_t>((m_iEnd[g] - m_iStart[g] + 1)) *
                                (m_jEnd[g] - m_jStart[g] + 1) *
                                (m_kEnd[g] - m_kStart[g] + 1);
  }
  size_t* ni = new size_t[mNumberOfGrids];
  size_t* nij = new size_t[mNumberOfGrids];
  for (int g = 0; g < mNumberOfGrids; g++) {
    ni[g] = (m_iEnd[g] - m_iStart[g] + 1);
    nij[g] = ni[g] * (m_jEnd[g] - m_jStart[g] + 1);
  }
  for (int s = 0; s < point_sources.size(); s++) {
    int g = point_sources[s]->m_grid;
    size_t key = gptr[g] + (point_sources[s]->m_i0 - m_iStart[g]) +
                 ni[g] * (point_sources[s]->m_j0 - m_jStart[g]) +
                 nij[g] * (point_sources[s]->m_k0 - m_kStart[g]);
    point_sources[s]->set_sort_key(key);
  }
  delete[] gptr;
  delete[] ni;
  delete[] nij;

  std::sort(point_sources.begin(), point_sources.end(), less_than);
  // set up array detecting sources belonging to idential points
  identsources.resize(1);
  identsources[0] = 0;
  int k = 0;
  while (identsources[k] < point_sources.size()) {
    int m = identsources[k];
    size_t key = point_sources[m]->m_key;
    while (m + 1 < point_sources.size() && point_sources[m + 1]->m_key == key)
      m++;
    identsources.push_back(m + 1);
    k++;
  }

  // Test
  int nrsrc = point_sources.size();
  int nrunique = identsources.size() - 1;
  int nrsrctot, nruniquetot;
  MPI_Reduce(&nrsrc, &nrsrctot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nrunique, &nruniquetot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (m_myRank == 0) {
    cout << "number of grid point  sources = " << nrsrctot << endl;
    cout << "number of unique g.p. sources = " << nruniquetot << endl;
  }

  ForceVector = SW4_NEW(Space::Managed, float_sw4[nrunique * 3]);
  ForceAddress = SW4_NEW(Space::Managed, float_sw4 * [nrunique * 3]);

  //   for( int s=0 ; s<m_identsources.size()-1 ; s++ )
  //      for( int i=m_identsources[s]; i< m_identsources[s+1] ; i++ )
  //	 std::cout << "src= " << i << " key=" << m_point_sources[i]->m_key <<
  //	    "grid= " << m_point_sources[i]->m_grid << " (i,j,k) = " <<
  //	    m_point_sources[i]->m_i0 << " " <<
  //	    m_point_sources[i]->m_j0 << " " <<
  //	    m_point_sources[i]->m_k0 << std::endl;
}
#ifdef FORCE_OMP
void EW::Force(float_sw4 a_t, vector<Sarray>& a_F,
               vector<GridPointSource*> point_sources,
               vector<int> identsources) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *f_ptr, om, ph, cv, h, zmin, omm, phm, amprho, ampmu, ampla;

  int g;

  if (m_twilight_forcing) {
    if (m_anisotropic) {
      float_sw4 phc[21];  // move these angles to the EW class

      // need to store all the phase angle constants somewhere
      for (int i = 0; i < 21; i++) phc[i] = i * 10 * M_PI / 180;

      for (g = 0; g < mNumberOfCartesianGrids; g++) {
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        h = mGridSize[g];
        zmin = m_zmin[g];

        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;

        if (m_croutines)
          tw_aniso_force_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                            a_t, om, cv, ph, omm, phm, amprho, phc, h, zmin);
        else
          tw_aniso_force(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                         a_t, om, cv, ph, omm, phm, amprho, phc, h, zmin);
      }  // end for all Cartesian grids
         // if (topographyExists()) {
      for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
        // g = mNumberOfGrids - 1;
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;

        if (m_croutines)
          tw_aniso_curvi_force_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                  f_ptr, a_t, om, cv, ph, omm, phm, amprho, phc,
                                  mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());
        else
          tw_aniso_curvi_force(ifirst, ilast, jfirst, jlast, kfirst, klast,
                               f_ptr, a_t, om, cv, ph, omm, phm, amprho, phc,
                               mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());

      }  // end if topographyExists

    } else {  // isotropic twilight forcing
      for (g = 0; g < mNumberOfCartesianGrids; g++) {
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        h = mGridSize[g];  // how do we define the grid size for the curvilinear
                           // grid?
        zmin = m_zmin[g];

        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        ampmu = m_twilight_forcing->m_ampmu;
        ampla = m_twilight_forcing->m_amplambda;
        if (usingSupergrid()) {
          float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
          float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
          float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
          if (m_croutines)
            forcingfortsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                             a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla, h,
                             zmin, omstrx, omstry, omstrz);
          else
            forcingfortsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                          f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                          &ampmu, &ampla, &h, &zmin, &omstrx, &omstry, &omstrz);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingfortsgatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                  f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                  ampmu, ampla, h, zmin, omstrx, omstry,
                                  omstrz);
            else
              forcingfortsgatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                               &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                               &amprho, &ampmu, &ampla, &h, &zmin, &omstrx,
                               &omstry, &omstrz);
          }
        } else {
          if (m_croutines)
            forcingfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                           a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla, h,
                           zmin);
          else
            forcingfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                        f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho, &ampmu,
                        &ampla, &h, &zmin);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingfortatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                                ampla, h, zmin);
            else
              forcingfortatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                             f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                             &ampmu, &ampla, &h, &zmin);
          }
        }
      }  // end for all Cartesian grids

      // if (topographyExists()) {
      for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
        // g = mNumberOfGrids - 1;
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        ampmu = m_twilight_forcing->m_ampmu;
        ampla = m_twilight_forcing->m_amplambda;
        if (usingSupergrid()) {
          float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
          float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
          float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
          if (m_croutines)
            forcingfortcsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                              f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                              ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                              mZ[g].c_ptr(), omstrx, omstry, omstrz);
          else
            forcingfortcsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                           f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                           &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                           mZ[g].c_ptr(), &omstrx, &omstry, &omstrz);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingfortsgattc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                   f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                   ampmu, ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                                   mZ[g].c_ptr(), omstrx, omstry, omstrz);
            else
              forcingfortsgattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                                &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                                &amprho, &ampmu, &ampla, mX[g].c_ptr(),
                                mY[g].c_ptr(), mZ[g].c_ptr(), &omstrx, &omstry,
                                &omstrz);
          }
        } else {
          if (m_croutines)
            forcingfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                            a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla,
                            mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());
          else
            forcingfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                         f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                         &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                         mZ[g].c_ptr());
          if (m_use_attenuation) {
            if (m_croutines)
              forcingfortattc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                 f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                 ampmu, ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                                 mZ[g].c_ptr());
            else
              forcingfortattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                              f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                              &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                              mZ[g].c_ptr());
          }
        }
      }
    }  // end isotropic case

  }  // end twilight

  else if (m_rayleigh_wave_test) {
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
  } else if (m_energy_test) {
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
  } else {
    // Default: m_point_source_test, m_lamb_test or full seismic case
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();

#pragma omp parallel for
    for (int r = 0; r < identsources.size() - 1; r++) {
      int s0 = identsources[r];
      int g = point_sources[s0]->m_grid;
      int i = point_sources[s0]->m_i0;
      int j = point_sources[s0]->m_j0;
      int k = point_sources[s0]->m_k0;
      float_sw4 f1 = 0, f2 = 0, f3 = 0;
      for (int s = identsources[r]; s < identsources[r + 1]; s++) {
        float_sw4 fxyz[3];
        point_sources[s]->getFxyz(a_t, fxyz);
        f1 += fxyz[0];
        f2 += fxyz[1];
        f3 += fxyz[2];
      }
      a_F[g](1, i, j, k) += f1;
      a_F[g](2, i, j, k) += f2;
      a_F[g](3, i, j, k) += f3;
    }
  }
}

//---------------------------------------------------------------------------
void EW::Force_tt(float_sw4 a_t, vector<Sarray>& a_F,
                  vector<GridPointSource*> point_sources,
                  vector<int> identsources) {
  SW4_MARK_FUNCTION;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  float_sw4 *f_ptr, om, ph, cv, h, zmin, omm, phm, amprho, ampmu, ampla;
  // std::cout<<"FORCE_TT CALLED\n";
  int g;

  if (m_twilight_forcing) {
    if (m_anisotropic) {
      float_sw4 phc[21];  // move these angles to the EW class

      // need to store all the phase angle constants somewhere
      phc[0] = 0;
      for (int i = 0; i < 21; i++) phc[i] = i * 10 * M_PI / 180;

      for (g = 0; g < mNumberOfCartesianGrids; g++) {
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        h = mGridSize[g];
        zmin = m_zmin[g];

        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;

        if (m_croutines)
          tw_aniso_force_tt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                               f_ptr, a_t, om, cv, ph, omm, phm, amprho, phc, h,
                               zmin);
        else
          tw_aniso_force_tt(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                            a_t, om, cv, ph, omm, phm, amprho, phc, h, zmin);
      }  // end for all Cartesian grids

      //      if (topographyExists()) {
      for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
        // g = mNumberOfGrids - 1;
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        if (m_croutines)
          tw_aniso_curvi_force_tt_ci(ifirst, ilast, jfirst, jlast, kfirst,
                                     klast, f_ptr, a_t, om, cv, ph, omm, phm,
                                     amprho, phc, mX[g].c_ptr(), mY[g].c_ptr(),
                                     mZ[g].c_ptr());
        else
          tw_aniso_curvi_force_tt(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                  f_ptr, a_t, om, cv, ph, omm, phm, amprho, phc,
                                  mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr());

      }  // end if topographyExists

    } else {  // isotropic twilight forcing
      for (g = 0; g < mNumberOfCartesianGrids; g++) {
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        h = mGridSize[g];  // how do we define the grid size for the curvilinear
                           // grid?
        zmin = m_zmin[g];

        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        ampmu = m_twilight_forcing->m_ampmu;
        ampla = m_twilight_forcing->m_amplambda;
        if (usingSupergrid()) {
          float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
          float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
          float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
          if (m_croutines)
            forcingttfortsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                               f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                               ampla, h, zmin, omstrx, omstry, omstrz);
          else
            forcingttfortsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                            f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                            &ampmu, &ampla, &h, &zmin, &omstrx, &omstry,
                            &omstrz);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingttfortsgatt_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                    f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                    ampmu, ampla, h, zmin, omstrx, omstry,
                                    omstrz);
            else
              forcingttfortsgatt(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                                 &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                                 &amprho, &ampmu, &ampla, &h, &zmin, &omstrx,
                                 &omstry, &omstrz);
          }
        } else {
          if (m_croutines)
            forcingttfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr,
                             a_t, om, cv, ph, omm, phm, amprho, ampmu, ampla, h,
                             zmin);
          else
            forcingttfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                          f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                          &ampmu, &ampla, &h, &zmin);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingttattfort_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                  f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                  ampmu, ampla, h, zmin);
            else
              forcingttattfort(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                               &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                               &amprho, &ampmu, &ampla, &h, &zmin);
          }
        }
      }
      // if (topographyExists()) {
      for (g = mNumberOfCartesianGrids; g < mNumberOfGrids; g++) {
        // g = mNumberOfGrids - 1;
        f_ptr = a_F[g].c_ptr();
        ifirst = m_iStart[g];
        ilast = m_iEnd[g];
        jfirst = m_jStart[g];
        jlast = m_jEnd[g];
        kfirst = m_kStart[g];
        klast = m_kEnd[g];
        om = m_twilight_forcing->m_omega;
        ph = m_twilight_forcing->m_phase;
        cv = m_twilight_forcing->m_c;
        omm = m_twilight_forcing->m_momega;
        phm = m_twilight_forcing->m_mphase;
        amprho = m_twilight_forcing->m_amprho;
        ampmu = m_twilight_forcing->m_ampmu;
        ampla = m_twilight_forcing->m_amplambda;
        if (usingSupergrid()) {
          float_sw4 omstrx = m_supergrid_taper_x[g].get_tw_omega();
          float_sw4 omstry = m_supergrid_taper_y[g].get_tw_omega();
          float_sw4 omstrz = m_supergrid_taper_z[g].get_tw_omega();
          if (m_croutines)
            forcingttfortcsg_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                                ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                                mZ[g].c_ptr(), omstrx, omstry, omstrz);
          else
            forcingttfortcsg(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                             f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                             &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                             mZ[g].c_ptr(), &omstrx, &omstry, &omstrz);
          if (m_use_attenuation) {
            if (m_croutines)
              forcingttfortsgattc_ci(
                  ifirst, ilast, jfirst, jlast, kfirst, klast, f_ptr, a_t, om,
                  cv, ph, omm, phm, amprho, ampmu, ampla, mX[g].c_ptr(),
                  mY[g].c_ptr(), mZ[g].c_ptr(), omstrx, omstry, omstrz);
            else
              forcingttfortsgattc(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                                  &klast, f_ptr, &a_t, &om, &cv, &ph, &omm,
                                  &phm, &amprho, &ampmu, &ampla, mX[g].c_ptr(),
                                  mY[g].c_ptr(), mZ[g].c_ptr(), &omstrx,
                                  &omstry, &omstrz);
          }
        } else {
          if (m_croutines)
            forcingttfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                              f_ptr, a_t, om, cv, ph, omm, phm, amprho, ampmu,
                              ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                              mZ[g].c_ptr());
          else
            forcingttfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst, &klast,
                           f_ptr, &a_t, &om, &cv, &ph, &omm, &phm, &amprho,
                           &ampmu, &ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                           mZ[g].c_ptr());
          if (m_use_attenuation) {
            if (m_croutines)
              forcingttattfortc_ci(ifirst, ilast, jfirst, jlast, kfirst, klast,
                                   f_ptr, a_t, om, cv, ph, omm, phm, amprho,
                                   ampmu, ampla, mX[g].c_ptr(), mY[g].c_ptr(),
                                   mZ[g].c_ptr());
            else
              forcingttattfortc(&ifirst, &ilast, &jfirst, &jlast, &kfirst,
                                &klast, f_ptr, &a_t, &om, &cv, &ph, &omm, &phm,
                                &amprho, &ampmu, &ampla, mX[g].c_ptr(),
                                mY[g].c_ptr(), mZ[g].c_ptr());
          }
        }
      }
    }  // end isotropic

  }  // end twilight

  else if (m_rayleigh_wave_test) {
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
  } else if (m_energy_test) {
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();
  } else {
    // Default: m_point_source_test, m_lamb_test or full seismic case
    for (int g = 0; g < mNumberOfGrids; g++) a_F[g].set_to_zero();

#pragma omp parallel for
    for (int r = 0; r < identsources.size() - 1; r++) {
      int s0 = identsources[r];
      int g = point_sources[s0]->m_grid;
      int i = point_sources[s0]->m_i0;
      int j = point_sources[s0]->m_j0;
      int k = point_sources[s0]->m_k0;
      float_sw4 f1 = 0, f2 = 0, f3 = 0;
      for (int s = identsources[r]; s < identsources[r + 1]; s++) {
        float_sw4 fxyz[3];
        point_sources[s]->getFxyztt(a_t, fxyz);
        f1 += fxyz[0];
        f2 += fxyz[1];
        f3 += fxyz[2];
      }
      a_F[g](1, i, j, k) += f1;
      a_F[g](2, i, j, k) += f2;
      a_F[g](3, i, j, k) += f3;
    }

    //     for( int s = 0 ; s < point_sources.size() ; s++ )
    //     {
    //	int g = point_sources[s]->m_grid;
    //        float_sw4 fxyz[3];
    //	point_sources[s]->getFxyztt(a_t,fxyz);
    //	a_F[g](1,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0)
    //+= fxyz[0];
    //	a_F[g](2,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0)
    //+= fxyz[1];
    //	a_F[g](3,point_sources[s]->m_i0,point_sources[s]->m_j0,point_sources[s]->m_k0)
    //+= fxyz[2];
    //     }
  }
}
#endif
#include "AllDims.h"
AllDims* EW::get_fine_alldimobject() {
  int g = mNumberOfCartesianGrids - 1;
  AllDims* fine = new AllDims(m_proc_array[0], m_proc_array[1], 1, 1,
                              m_global_nx[g], 1, m_global_ny[g], 1,
                              m_global_nz[g], m_ghost_points, m_ppadding);
  return fine;
}
//-----------------------------------------------------------------------
void EW::extractTopographyFromSfile(std::string a_topoFileName) {
  double start_time, end_time;
  start_time = MPI_Wtime();
#ifdef USE_HDF5
  /* int verbose = mVerbose; */
  std::string rname = "EW::extractTopographyFromSfile";
  Sarray gridElev;
  herr_t ierr;
  hid_t file_id, dataset_id, datatype_id, h5_dtype, group_id, dataspace_id,
      attr_id;  //, plist_id;
  int prec;
  double lonlataz[3];

  /* plist_id = H5Pcreate(H5P_FILE_ACCESS); */
  /* H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL); */
  /* file_id = H5Fopen(a_topoFileName.c_str(), H5F_ACC_RDONLY, plist_id); */
  if (m_myRank == 0) {
    file_id = H5Fopen(a_topoFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "Could not open hdf5 file: " << a_topoFileName.c_str() << endl;
      MPI_Abort(MPI_COMM_WORLD, file_id);
    }
  }

  // Origin longitude, latitude, azimuth
  if (m_myRank == 0) {
    attr_id =
        H5Aopen(file_id, "Origin longitude, latitude, azimuth", H5P_DEFAULT);
    ASSERT(attr_id >= 0);
    ierr = H5Aread(attr_id, H5T_NATIVE_DOUBLE, lonlataz);
    ASSERT(ierr >= 0);
    H5Aclose(attr_id);
  }
  MPI_Bcast(&lonlataz, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double alpha = lonlataz[2], lon0 = lonlataz[0], lat0 = lonlataz[1];

  CHECK_INPUT(fabs(alpha - mGeoAz) < 1e-6,
              "ERROR: Sfile azimuth must be equal to coordinate system azimuth"
                  << " azimuth on sfile = " << alpha
                  << " azimuth of coordinate sytem = " << mGeoAz
                  << " difference = " << alpha - mGeoAz);

  // Ngrids - int, number of 3D grids in the file
  int npatches;
  if (m_myRank == 0) {
    attr_id = H5Aopen(file_id, "ngrids", H5P_DEFAULT);
    ierr = H5Aread(attr_id, H5T_NATIVE_INT, &npatches);
    ASSERT(ierr >= 0);
    H5Aclose(attr_id);
  }
  MPI_Bcast(&npatches, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (m_myRank == 0 && mVerbose >= 2) {
    printf("Sfile header: \n");
    printf("              azimuth=%e, lon0=%e, lat0=%e\n", alpha, lon0, lat0);
    printf("              nblocks=%i\n", npatches);
  }

  double hh;
  if (m_myRank == 0) {
    attr_id = H5Aopen(file_id, "Coarsest horizontal grid spacing", H5P_DEFAULT);
    ASSERT(attr_id >= 0);
    ierr = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &hh);
    ASSERT(ierr >= 0);
    H5Aclose(attr_id);
  }
  MPI_Bcast(&hh, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ---------- read topography on file into array gridElev
  hsize_t dims[2];
  char intf_name[32];

  if (m_myRank == 0) {
    group_id = H5Gopen(file_id, "Z_interfaces", H5P_DEFAULT);
    ASSERT(group_id >= 0);

    sprintf(intf_name, "z_values_%d", 0);
    dataset_id = H5Dopen(group_id, intf_name, H5P_DEFAULT);
    ASSERT(dataset_id >= 0);

    dataspace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    H5Sclose(dataspace_id);

    datatype_id = H5Dget_type(dataset_id);
    prec = (int)H5Tget_size(datatype_id);
    H5Tclose(datatype_id);
  }
  MPI_Bcast(dims, 2, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&prec, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int nitop = (int)dims[0], njtop = (int)dims[1];

  if (m_myRank == 0 && mVerbose >= 2) {
    printf("Topography header\n");
    printf("  hh=%e\n", hh);
    printf("  ni=%i, nj=%i\n", nitop, njtop);
  }

  bool roworder = true;

  // Depending on the precision of the sfile and sw4, need to convert between
  // float and doulbe
  void* in_data;
  float* f_data = new float[nitop * njtop];
  double* d_data = new double[nitop * njtop];
  if (prec == 4) {
    h5_dtype = H5T_NATIVE_FLOAT;
    in_data = (void*)f_data;
  } else if (prec == 8) {
    h5_dtype = H5T_NATIVE_DOUBLE;
    in_data = (void*)d_data;
  }

  if (m_myRank == 0) {
    ierr =
        H5Dread(dataset_id, h5_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, in_data);
    ASSERT(ierr >= 0);
    H5Dclose(dataset_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
  }
  MPI_Bcast(in_data, nitop * njtop * prec, MPI_CHAR, 0, MPI_COMM_WORLD);

  float_sw4* data = new float_sw4[nitop * njtop];
  for (int i = 0; i < nitop * njtop; i++) {
    if (prec == 4)
      data[i] = -(float_sw4)f_data[i];
    else if (prec == 8)
      data[i] = -(float_sw4)d_data[i];
  }

  delete[] f_data;
  delete[] d_data;

  gridElev.define(1, nitop, 1, njtop, 1, 1);
  gridElev.assign(data);

  if (m_myRank == 0 && mVerbose >= 2) {
    printf("1st topo (float) data=%e, gridElev(1,1,1)=%e\n", data[0],
           gridElev(1, 1, 1));
    printf("last topo (float) data=%e, gridElev(ni,nj,1)=%e\n",
           data[nitop * njtop - 1], gridElev(nitop, njtop, 1));
    // get min and max
    float tmax = -9e-10, tmin = 9e+10;
    for (int q = 0; q < nitop * njtop; q++) {
      if (data[q] > tmax) tmax = data[q];
      if (data[q] < tmin) tmin = data[q];
    }
    printf("topo max (float)=%e, min (float)=%e\n", tmax, tmin);
  }
  delete[] data;

  if (roworder) gridElev.transposeik();

  // ---------- done reading

  // ---------- origin on file
  double x0, y0;  // Origin on grid file
  computeCartesianCoord(x0, y0, lon0, lat0);
  if (m_myRank == 0 && mVerbose >= 2) {
    printf("mat-lon0=%e mat-lat0=%e, comp-x0=%e, commp-y0=%e\n", lon0, lat0, x0,
           y0);
  }

  // Topography read, next interpolate to the computational grid
  int topLevel = mNumberOfGrids - 1;

  float_sw4 topomax = -1e30, topomin = 1e30;
#pragma omp parallel for reduction(max : topomax) reduction(min : topomin)
  for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i) {
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j) {
      float_sw4 x = (i - 1) * mGridSize[topLevel];
      float_sw4 y = (j - 1) * mGridSize[topLevel];
      int i0 = static_cast<int>(trunc(1 + (x - x0) / hh));
      int j0 = static_cast<int>(trunc(1 + (y - y0) / hh));
      // test
      float_sw4 xmat0 = (i0 - 1) * hh, ymat0 = (j0 - 1) * hh;
      float_sw4 xmatx = x - x0, ymaty = y - y0;

      if (mVerbose >= 3) {
        if (xmatx < xmat0 || xmatx > xmat0 + hh)
          printf("WARNING: i0=%i is out of bounds for x=%e, xmatx=%e\n", i0, x,
                 xmatx);
        if (ymaty < ymat0 || ymaty > ymat0 + hh)
          printf("WARNING: i0=%i is out of bounds for y=%e, ymaty=%e\n", i0, y,
                 ymaty);
      }

      bool extrapol = false;
      if (i0 < -1) {
        extrapol = true;
        i0 = 1;
      } else if (i0 < 2)
        i0 = 2;

      if (i0 > nitop + 1) {
        extrapol = true;
        i0 = nitop;
      } else if (i0 > nitop - 2)
        i0 = nitop - 2;

      if (j0 < -1) {
        extrapol = true;
        j0 = 1;
      } else if (j0 < 2)
        j0 = 2;

      if (j0 > njtop + 1) {
        extrapol = true;
        j0 = njtop;
      } else if (j0 > njtop - 2)
        j0 = njtop - 2;

      if (!extrapol) {
        float_sw4 q = (x - x0 - (i0 - 1) * hh) / hh;
        float_sw4 r = (y - y0 - (j0 - 1) * hh) / hh;
        float_sw4 Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1,
            tjp2;
        Qim1 = (q) * (q - 1) * (q - 2) / (-6.);
        Qi = (q + 1) * (q - 1) * (q - 2) / (2.);
        Qip1 = (q + 1) * (q) * (q - 2) / (-2.);
        Qip2 = (q + 1) * (q) * (q - 1) / (6.);

        Rjm1 = (r) * (r - 1) * (r - 2) / (-6.);
        Rj = (r + 1) * (r - 1) * (r - 2) / (2.);
        Rjp1 = (r + 1) * (r) * (r - 2) / (-2.);
        Rjp2 = (r + 1) * (r) * (r - 1) / (6.);

        if (mVerbose >= 3) {
          if (i0 < 2 || i0 > nitop - 2)
            printf("WARNING: topo interp out of bounds i0=%i, nitop=%i\n", i0,
                   nitop);
          if (j0 < 2 || j0 > njtop - 2)
            printf("WARNING: topo interp out of bounds j0=%i, njtop=%i\n", j0,
                   njtop);
        }

        tjm1 = Qim1 * gridElev(i0 - 1, j0 - 1, 1) +
               Qi * gridElev(i0, j0 - 1, 1) +
               Qip1 * gridElev(i0 + 1, j0 - 1, 1) +
               Qip2 * gridElev(i0 + 2, j0 - 1, 1);
        tj = Qim1 * gridElev(i0 - 1, j0, 1) + Qi * gridElev(i0, j0, 1) +
             Qip1 * gridElev(i0 + 1, j0, 1) + Qip2 * gridElev(i0 + 2, j0, 1);
        tjp1 = Qim1 * gridElev(i0 - 1, j0 + 1, 1) +
               Qi * gridElev(i0, j0 + 1, 1) +
               Qip1 * gridElev(i0 + 1, j0 + 1, 1) +
               Qip2 * gridElev(i0 + 2, j0 + 1, 1);
        tjp2 = Qim1 * gridElev(i0 - 1, j0 + 2, 1) +
               Qi * gridElev(i0, j0 + 2, 1) +
               Qip1 * gridElev(i0 + 1, j0 + 2, 1) +
               Qip2 * gridElev(i0 + 2, j0 + 2, 1);
        mTopo(i, j, 1) = Rjm1 * tjm1 + Rj * tj + Rjp1 * tjp1 + Rjp2 * tjp2;
      } else {
        if (mVerbose >= 3)
          printf(
              "INFO: topo extrapolated for i=%i, j=%i, x=%e, y=%e, i0=%i, "
              "j0=%i\n",
              i, j, x, y, i0, j0);

        mTopo(i, j, 1) = gridElev(i0, j0, 1);
      }

      // test
      if (mTopo(i, j, 1) > topomax) topomax = mTopo(i, j, 1);
      if (mTopo(i, j, 1) < topomin) topomin = mTopo(i, j, 1);

    }  // end for j
  }    // end for i

  MPI_Barrier(MPI_COMM_WORLD);
  end_time = MPI_Wtime();
  if (m_myRank == 0)
    printf("Read topography from sfile time=%e seconds\n",
           end_time - start_time);

  // test
  if (m_myRank == 0 && mVerbose >= 2) {
    printf("Topo variation on comp grid: max=%e min=%e\n", topomax, topomin);
  }
#else
  if (m_myRank == 0)
    printf(
        "WARNING: sw4 not compiled with hdf5=yes, ignoring read sfile, "
        "abort!\n");
  MPI_Abort(MPI_COMM_WORLD, -1);
#endif  // #ifdef USE_HDF5
}

#ifdef USE_HDF5
static void read_hdf5_attr(hid_t loc, hid_t dtype, const char* name,
                           void* data) {
  hid_t attr_id;
  int ierr;
  attr_id = H5Aopen(loc, name, H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  ierr = H5Aread(attr_id, dtype, data);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);
}

static char* read_hdf5_attr_str(hid_t loc, const char* name) {
  hid_t attr_id, dtype;
  int ierr;
  char* data = NULL;

  attr_id = H5Aopen(loc, name, H5P_DEFAULT);
  ASSERT(attr_id >= 0);

  dtype = H5Aget_type(attr_id);

  ierr = H5Aread(attr_id, dtype, &data);
  ASSERT(ierr >= 0);

  H5Tclose(dtype);
  H5Aclose(attr_id);

  /* fprintf(stderr, "Read data: [%s]\n", data); */
  return data;
}
#endif

//-----------------------------------------------------------------------
void EW::extractTopographyFromGMG(std::string a_topoFileName) {
  double start_time, end_time;
  start_time = MPI_Wtime();
#ifdef USE_HDF5
  Sarray gridElev;
  herr_t ierr;
  hid_t file_id, dataset_id, datatype_id, group_id, dataspace_id;
  int prec, str_len;
  double az = 0, origin_x = 0, origin_y = 0, hh = 0, alpha = 0;
  hsize_t dims[2];
  char* crs_to = NULL;

  if (m_myRank == 0) {
    file_id = H5Fopen(a_topoFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "Could not open hdf5 file: " << a_topoFileName.c_str() << endl;
      MPI_Abort(MPI_COMM_WORLD, file_id);
    }

    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "origin_x", &origin_x);
    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "origin_y", &origin_y);
    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "y_azimuth", &az);

    group_id = H5Gopen(file_id, "surfaces", H5P_DEFAULT);
    ASSERT(group_id >= 0);

    dataset_id = H5Dopen(group_id, "topography_bathymetry", H5P_DEFAULT);
    ASSERT(dataset_id >= 0);

    dataspace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    H5Sclose(dataspace_id);

    datatype_id = H5Dget_type(dataset_id);
    prec = (int)H5Tget_size(datatype_id);
    H5Tclose(datatype_id);

    read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "resolution_horiz", &hh);

    crs_to = read_hdf5_attr_str(file_id, "crs");
    str_len = (int)(strlen(crs_to) + 1);
  }

  MPI_Bcast(&origin_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&origin_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&az, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hh, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(dims, 2, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&prec, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&str_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (m_myRank != 0) crs_to = (char*)malloc(str_len * sizeof(char));

  MPI_Bcast(crs_to, str_len, MPI_CHAR, 0, MPI_COMM_WORLD);

  // For some reason origin_x is not correctly read sometimes
  if (origin_x < 1.0) {
    origin_x = 99286.2;
    if (m_myRank == 0)
      printf("GMG origin_x read zero value, correct to 99286.2 \n");
  }

  ASSERT(origin_x > 0);
  ASSERT(origin_y > 0);
  ASSERT(az > 0);

  // Convert GMG az to SW4 az
  alpha = az - 180.0;

  CHECK_INPUT(fabs(alpha - mGeoAz) < 1e-6,
              "ERROR: GMG azimuth must be equal to coordinate system azimuth"
                  << " azimuth on GMG = " << alpha
                  << " azimuth of coordinate sytem = " << mGeoAz
                  << " difference = " << alpha - mGeoAz);

  if (m_myRank == 0 && mVerbose >= 2) {
    printf("GMG header: azimuth=%e, origin_x=%f, origin_y=%f\n", az, origin_x,
           origin_y);
    printf("            hh=%e, ni=%llu, nj=%llu\n", hh, dims[0], dims[1]);
  }

  float* f_data = new float[dims[0] * dims[1]];

  if (m_myRank == 0) {
    ierr = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   f_data);
    ASSERT(ierr >= 0);
    H5Dclose(dataset_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
  }
  MPI_Bcast(f_data, dims[0] * dims[1], MPI_FLOAT, 0, MPI_COMM_WORLD);

  // Topography read, next interpolate to the computational grid
  int topLevel = mNumberOfGrids - 1;

  float_sw4 topomax = -1e30, topomin = 1e30;

  /* printf("x0=%f, y0=%f\n", x0, y0); */
  /* printf("topoGMG: m_iStart %d, m_iEnd %d, m_jStart %d, m_jEnd %d\n", */
  /*         m_iStart[topLevel],  m_iEnd[topLevel],  m_jStart[topLevel],
   * m_jEnd[topLevel]); */

  const double yazimuthRad = az * M_PI / 180.0;
  const double cosAz = cos(yazimuthRad);
  const double sinAz = sin(yazimuthRad);

  /* fprintf(stderr, "origin xy: %f %f\n", origin_x, origin_y); */

  // proj is not thread safe, so no omp pragma here
  for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i) {
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j) {
      float_sw4 x = (i - 1) * mGridSize[topLevel];
      float_sw4 y = (j - 1) * mGridSize[topLevel];

      double sw4_lon, sw4_lat, gmg_x, gmg_y, gmg_x0, gmg_y0;
      computeGeographicCoord(x, y, sw4_lon, sw4_lat);
      /* printf("\ncomputeGeographicCoord: %f %f %f %f\n", x, y, sw4_lon,
       * sw4_lat); */

      // GMG x/y, lat/lon is switched from sw4 CRS
      computeCartesianCoordGMG(gmg_y0, gmg_x0, sw4_lon, sw4_lat, crs_to);
      /* printf("computeCartesianCoordGMG : %f %f %f %f\n", gmg_x0, gmg_y0,
       * sw4_lon, sw4_lat); */

      const double xRel = gmg_x0 - origin_x;
      const double yRel = gmg_y0 - origin_y;
      gmg_x = xRel * cosAz - yRel * sinAz;
      gmg_y = xRel * sinAz + yRel * cosAz;
      /* printf("converted gmg xy: %f, %f, origin: %f %f\n", gmg_x, gmg_y,
       * origin_x, origin_y); */

      int i0 = static_cast<int>(floor(gmg_x / hh));
      int j0 = static_cast<int>(floor(gmg_y / hh));

      double fac0 = (gmg_y - j0 * hh) / hh;
      double fac1 = (gmg_x - i0 * hh) / hh;

      /* printf("x=%f, y=%f, i0=%d, j0=%d\n", x, y, i0, j0); */
      /* printf("interp points: %f %f %f %f\n", f_data[i0*dims[1]+j0],
       * f_data[(i0+1)*dims[1]+j0], f_data[i0*dims[1]+j0+1],
       * f_data[(i0+1)*dims[1]+j0+1]); */

      // Linear interpolation with 4 surrounding points
      float_sw4 mytopo =
          f_data[i0 * dims[1] + j0] +
          (f_data[i0 * dims[1] + j0 + 1] - f_data[i0 * dims[1] + j0]) * fac0 +
          (f_data[(i0 + 1) * dims[1] + j0] +
           (f_data[(i0 + 1) * dims[1] + j0 + 1] -
            f_data[(i0 + 1) * dims[1] + j0]) *
               fac0 -
           (f_data[i0 * dims[1] + j0] +
            (f_data[i0 * dims[1] + j0 + 1] - f_data[i0 * dims[1] + j0]) *
                fac0)) *
              fac1;
      /* printf("Calculated topo: %f, fac %f %f\n", mytopo, fac0, fac1); */
      if (mytopo > topomax) topomax = mytopo;
      if (mytopo < topomin) topomin = mytopo;

      mTopo(i, j, 1) = mytopo;
    }  // end for j
  }    // end for i

  if (crs_to) free(crs_to);

  delete[] f_data;

  MPI_Barrier(MPI_COMM_WORLD);
  end_time = MPI_Wtime();

  if (m_myRank == 0) {
    printf("Read topography from GeoModelGrids time=%e seconds\n",
           end_time - start_time);
    if (mVerbose >= 2) {
      printf("Topo corners %f, %f, %f, %f\n",
             mTopo(m_iStart[topLevel], m_jStart[topLevel], 1),
             mTopo(m_iEnd[topLevel], m_jStart[topLevel], 1),
             mTopo(m_iEnd[topLevel], m_jStart[topLevel], 1),
             mTopo(m_iEnd[topLevel], m_jEnd[topLevel], 1));
      printf("Topo variation on comp grid: max=%e min=%e\n", topomax, topomin);
    }
  }
#else
  if (m_myRank == 0)
    printf(
        "WARNING: sw4 not compiled with hdf5=yes, ignoring read GMG, abort!\n");
  MPI_Abort(MPI_COMM_WORLD, -1);
#endif
}

// CURVI_MR_CODE ADDED HERE
#include "TestTwilight.h"

TestTwilight* EW::create_twilight() {
  if (m_twilight_forcing != 0)
    return new TestTwilight(
        m_twilight_forcing->m_omega, m_twilight_forcing->m_c,
        m_twilight_forcing->m_phase, m_twilight_forcing->m_momega,
        m_twilight_forcing->m_mphase, m_twilight_forcing->m_amprho,
        m_twilight_forcing->m_ampmu, m_twilight_forcing->m_amplambda);
  else
    return 0;
}

#include "TestEcons.h"

TestEcons* EW::create_energytest() {
  if (m_energy_test != 0)
    return new TestEcons(m_energy_test->m_stochastic_amp,
                         m_energy_test->m_cpcsratio);
  else
    return 0;
}
TestPointSource* EW::get_point_source_test() { return m_point_source_test; }
void EW::load_balance() {
  for (int g = 0; g < mNumberOfGrids; g++) {
    size_t local_size = (m_kEndInt[g] - m_kStartInt[g] + 1) *
                        (m_jEndInt[g] - m_jStartInt[g] + 1) *
                        (m_iEndInt[g] - m_iStartInt[g] + 1);
    size_t global_size = m_global_nx[g] * m_global_ny[g] * m_global_nz[g];

    float perc = static_cast<float>(local_size) / global_size * 100.0;
    std::cout << getRank() << " Grid # " << g << " " << perc << "%\n";
  }
}
