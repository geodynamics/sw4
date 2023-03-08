//  SW4 LICENSE
// ----------------------------------------------------------------------
// SW4 - Seismic Waves, 4th order
// ----------------------------------------------------------------------
// Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// Written by:
// N. Anders Petersson (petersson1@llnl.gov)
// Bjorn Sjogreen      (sjogreen2@llnl.gov)
//
// LLNL-CODE-643337
//
// All rights reserved.
//
// This file is part of SW4, Version: 1.0
//
// Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License (as published by
// the Free Software Foundation) version 2, dated June 1991.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// conditions of the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#include "EW.h"
#include "caliper.h"

#define SQR(x) ((x) * (x))

//-----------------------------------------------------------------------
void projectmtrl(int ib, int ie, int jb, int je, int kb, int ke, int iba,
                 int iea, int jba, int jea, int kba, int kea,
                 float_sw4* __restrict__ rho, float_sw4* __restrict__ mu,
                 float_sw4* __restrict__ lambda, float_sw4 dt, float_sw4 h,
                 float_sw4 cfl, float_sw4 vsmin, float_sw4 rhoscale,
                 float_sw4 muscale, float_sw4 lascale, int& info) {
  SW4_MARK_FUNCTION;
  //*********************************************************************
  //*
  //* Projection of the material arrays to satisfy the following constraints
  //*
  //*    1. CFL condition max(vp) < h*cflmax/dt
  //*    2. Smallest resolvable speed, min(vs) > vsmin
  //*    3. Positive density rho > eps
  //*    4. Positive mu,  mu > eps
  //*    5. Non-negative lambda, lambda >= 0
  //*
  //* Input: ib, ie, jb, je, kb, ke - Declared dimensions of arrays
  //*        iba, iea, jba, jea, kba, kea - Apply projection to this subarray
  //*        rho - Density
  //*        mu  - Lam\'e parameter
  //*        lambda - Lam\'e parameter
  //*        h - Grid spacing
  //*        cfl - Maximum stable cfl number
  //*        vsmin - Smallest resolvable wave speed, condition is
  //*                not enforced if vsmin < 0
  //*        rhoscale, muscale, lascale - Scaling parameters for rho, mu, lambda
  //*
  //* Output: rho, mu, lambda - Material is projected to satisfy the
  // constraints.
  //*         info - 0   --> No correction was needed
  //*               -1   --> Projection was unsuccessful
  //*                1   --> Negative density or mu was corrected.
  //*                10  --> CFL limit was enforced
  //*                100 --> smallest resolvable wave speed was enforced.
  //*       ( can be a sum of 1, 10, 100 if more than one condition was enforced
  //)
  //*
  //*********************************************************************

  float_sw4 acof = (cfl * h / dt) * (cfl * h / dt);
  float_sw4 bcof = -1;
  if (vsmin > 0) bcof = vsmin * vsmin;
  float_sw4 zcorr = 0;
  float_sw4 cflcorr = 0;
  float_sw4 lowvscorr = 0;
  // index = i-ib + ni*(j-jb) + ni*nj*(k-kb)
  const size_t ni = (ie - ib + 1);
  const size_t nij = ni * (je - jb + 1);
  const ssize_t base = -ib - ni * jb - nij * kb;

  const float_sw4 muscale2 = muscale * muscale;
  const float_sw4 lascale2 = lascale * lascale;
  const float_sw4 rhoscale2 = rhoscale * rhoscale;
  for (int k = kba; k <= kea; k++)
    for (int j = jba; j <= jea; j++)
      for (int i = iba; i <= iea; i++) {
        size_t ind = base + i + j * ni + k * nij;
        if (mu[ind] <= 0) {
          mu[ind] = 1e-6 * muscale;
          zcorr = 1;
        }
        if (rho[ind] <= 0) {
          rho[ind] = 1e-6 * rhoscale;
          zcorr = 1;
        }
        if (lambda[ind] < 0) {
          lambda[ind] = 0;
          zcorr = 1;
        }
        float_sw4 k1, k2;
        float_sw4 c1 = acof * rho[ind] - 4 * mu[ind] - lambda[ind];
        float_sw4 c2 = mu[ind] - bcof * rho[ind];
        if (c1 < 0 && c2 >= 0) {
          //* correct first condition
          k1 = c1 / (SQR(acof * rhoscale) * +16 * muscale2 + lascale2);
          mu[ind] = mu[ind] + muscale2 * 4 * k1;
          lambda[ind] = lambda[ind] + lascale2 * k1;
          rho[ind] = rho[ind] - rhoscale2 * acof * k1;
          c1 = acof * rho[ind] - 4 * mu[ind] - lambda[ind];
          c2 = mu[ind] - bcof * rho[ind];
          cflcorr = 1;
        }
        if (c1 >= 0 && c2 < 0) {
          //* correct second condition
          k2 = c2 / (muscale2 + SQR(bcof * rhoscale));
          mu[ind] = mu[ind] - muscale2 * k2;
          rho[ind] = rho[ind] + rhoscale2 * bcof * k2;
          c1 = acof * rho[ind] - 4 * mu[ind] - lambda[ind];
          c2 = mu[ind] - bcof * rho[ind];
          lowvscorr = 1;
        }
        if (c1 < 0 || c2 < 0) {
          //* correct both conditions
          float_sw4 a11 = -(SQR(acof * rhoscale) + 16 * muscale2 + lascale2);
          float_sw4 a12 = acof * bcof * rhoscale2 + 4 * muscale2;
          float_sw4 a22 = -(muscale2 + SQR(bcof * rhoscale));
          float_sw4 idet = 1 / (a11 * a22 - a12 * a12);
          k1 = idet * (-a22 * c1 + a12 * c2);
          k2 = idet * (a12 * c1 - a11 * c2);
          cflcorr = 1;
          lowvscorr = 1;
        }
        //* Was it successful ?
        if (mu[ind] <= 0 || lambda[ind] < 0 || rho[ind] <= 0) {
          info = -1;
          return;
        }
      }
  info = 0;
  if (zcorr == 1) info = 1;
  if (cflcorr == 1) info = info + 10;
  if (lowvscorr == 1) info = info + 100;
}

//-----------------------------------------------------------------------
void projectmtrlc(int ib, int ie, int jb, int je, int kb, int ke, int iba,
                  int iea, int jba, int jea, int kba, int kea,
                  float_sw4* __restrict__ rho, float_sw4* __restrict__ mu,
                  float_sw4* __restrict__ lambda, float_sw4 dt,
                  float_sw4* __restrict__ met, float_sw4* __restrict__ jac,
                  float_sw4 cfl, float_sw4 vsmin, float_sw4 rhoscale,
                  float_sw4 muscale, float_sw4 lascale, int& info) {
  SW4_MARK_FUNCTION;
  //*********************************************************************
  //*
  //* Same as PROJECTMTRL, but for the curvilinear grid.
  //*
  //*********************************************************************

  //      implicit none
  //      integer ib, ie, jb, je, kb, ke, iba, iea, jba, jea, kba, kea
  //      integer i, j, k, info, zcorr, cflcorr, lowvscorr
  //      real*8 rho(ib:ie,jb:je,kb:ke), mu(ib:ie,jb:je,kb:ke)
  //      real*8 lambda(ib:ie,jb:je,kb:ke), jac(ib:ie,jb:je,kb:ke)
  //      real*8 met(4,ib:ie,jb:je,kb:ke)
  //      real*8 acof, bcof, c1, c2, acofbase, acof1, acof2, as2, d
  //      real*8 dt, cfl, rhoscale, muscale, lascale
  //      real*8 vsmin, k1, k2, idet, a11, a12, a22

  float_sw4 acofbase = cfl * cfl / (dt * dt);
  float_sw4 bcof = -1;
  if (vsmin > 0) bcof = vsmin * vsmin;
  float_sw4 zcorr = 0;
  float_sw4 cflcorr = 0;
  float_sw4 lowvscorr = 0;
  // index = i-ib + ni*(j-jb) + ni*nj*(k-kb)
  const size_t ni = (ie - ib + 1);
  const size_t nij = ni * (je - jb + 1);
  const size_t nijk = nij * (ke - kb + 1);
  const ssize_t base = -ib - ni * jb - nij * kb;

  const float_sw4 muscale2 = muscale * muscale;
  const float_sw4 lascale2 = lascale * lascale;
  const float_sw4 rhoscale2 = rhoscale * rhoscale;
  for (int k = kba; k <= kea; k++)
    for (int j = jba; j <= jea; j++)
      for (int i = iba; i <= iea; i++) {
        size_t ind = base + i + j * ni + k * nij;
        if (mu[ind] <= 0) {
          mu[ind] = 1e-6 * muscale;
          zcorr = 1;
        }
        if (rho[ind] <= 0) {
          rho[ind] = 1e-6 * rhoscale;
          zcorr = 1;
        }
        if (lambda[ind] < 0) {
          lambda[ind] = 0;
          zcorr = 1;
        }

        float_sw4 as2 = SQR(met[ind + nijk]) + SQR(met[ind + 2 * nijk]) +
                        SQR(met[ind + 3 * nijk]);
        float_sw4 d = sqrt(SQR(met[ind] * met[ind] + as2) -
                           4 * SQR(met[ind] * met[ind + 3 * nijk]));
        float_sw4 acof1, acof2, acof, c1, c2, k1, k2;
        if (as2 + d > SQR(met[ind])) {
          acof1 = 0.5 * (5 * SQR(met[ind]) + 3 * as2 + d);
          acof2 = 0.5 * (SQR(met[ind]) + as2 + d);
        } else {
          acof1 = 3 * SQR(met[ind]) + as2;
          acof2 = SQR(met[ind]);
        }
        acof = acofbase * jac[ind];
        c1 = acof * rho[ind] - acof1 * mu[ind] - acof2 * lambda[ind];
        c2 = mu[ind] - bcof * rho[ind];
        if (c1 < 0 && c2 >= 0) {
          //* correct first condition
          k1 = c1 / (SQR(acof * rhoscale) + SQR(acof1 * muscale) +
                     SQR(acof2 * lascale));
          mu[ind] = mu[ind] + muscale2 * acof1 * k1;
          lambda[ind] = lambda[ind] + lascale2 * acof2 * k1;
          rho[ind] = rho[ind] - rhoscale2 * acof * k1;
          c1 = acof * rho[ind] - acof1 * mu[ind] - acof2 * lambda[ind];
          c2 = mu[ind] - bcof * rho[ind];
          cflcorr = 1;
        }
        if (c1 >= 0 && c2 < 0) {
          //* correct second condition
          k2 = c2 / (muscale2 + SQR(bcof * rhoscale));
          mu[ind] = mu[ind] - muscale2 * k2;
          rho[ind] = rho[ind] + rhoscale2 * bcof * k2;
          c1 = acof * rho[ind] - acof1 * mu[ind] - acof2 * lambda[ind];
          c2 = mu[ind] - bcof * rho[ind];
          lowvscorr = 1;
        }
        if (c1 < 0 || c2 < 0) {
          //* correct both conditions
          float_sw4 a11 = -(SQR(acof * rhoscale) + SQR(acof1 * muscale) +
                            SQR(acof2 * lascale));
          float_sw4 a12 = acof * bcof * rhoscale2 + acof1 * muscale2;
          float_sw4 a22 = -(muscale2 + SQR(bcof * rhoscale));
          float_sw4 idet = 1 / (a11 * a22 - a12 * a12);
          k1 = idet * (-a22 * c1 + a12 * c2);
          k2 = idet * (a12 * c1 - a11 * c2);
          cflcorr = 1;
          lowvscorr = 1;
        }
        //* Was it successful ?
        if (mu[ind] <= 0 || lambda[ind] < 0 || rho[ind] <= 0) {
          info = -1;
          return;
        }
      }
  info = 0;
  if (zcorr == 1) info = 1;
  if (cflcorr == 1) info = info + 10;
  if (lowvscorr == 1) info = info + 100;
}

//-----------------------------------------------------------------------
void checkmtrl(int ib, int ie, int jb, int je, int kb, int ke,
               float_sw4* __restrict__ rho, float_sw4* __restrict__ mu,
               float_sw4* __restrict__ lambda, float_sw4 dt, float_sw4 h,
               float_sw4 limits[10]) {
  SW4_MARK_FUNCTION;
  //-----------------------------------------------------------------------
  //      subroutine CHECKMTRL( ib,  ie,  jb,  je,  kb,  ke,
  //     *         rho, mu, lambda, dt, h, limits )

  //***********************************************************************
  //*
  //* Check if the material satisfies the following constraints
  //*
  //*    1. CFL condition max(vp) < h*cflmax/dt
  //*    2. Smallest resolvable speed, min(vs) > vsmin
  //*    3. Positive density rho > eps
  //*    4. Positive mu,  mu > eps
  //*    5. Non-negative lambda, lambda >= 0
  //*
  //* Input: ib, ie, jb, je, kb, ke - Declared dimensions of arrays
  //*        rho - Density
  //*        mu  - Lam\'e parameter
  //*        lambda - Lam\'e parameter
  //*        h - Grid spacing
  //*        dt - time step
  //*
  //* Output: limits - vector of 8 elements:
  //*            (rhomin,rhomax,mumin,mumax,lamin,lamax,cfl2max,vs2min)
  //*
  //*********************************************************************

  //      implicit none
  //      integer ib, ie, jb, je, kb, ke
  //      integer i, j, k
  //      real*8 rho(ib:ie,jb:je,kb:ke), mu(ib:ie,jb:je,kb:ke)
  //      real*8 lambda(ib:ie,jb:je,kb:ke), acof, c1, c2
  //      real*8 dt, h
  //      real*8 rhmin, rhmax, mumin, mumax, lamin, lamax
  //      real*8 cfl2max, vs2min, limits(10), tmthlmin, tmthlmax

  // index = i-ib + ni*(j-jb) + ni*nj*(k-kb)
  const size_t ni = (ie - ib + 1);
  const size_t nij = ni * (je - jb + 1);
  const size_t nijk = nij * (ke - kb + 1);
  //   const ssize_t base = -ib-ni*jb-nij*kb;

  float_sw4 rhmin = 1e38;
  float_sw4 rhmax = -1e38;
  float_sw4 mumin = 1e38;
  float_sw4 mumax = -1e38;
  float_sw4 lamin = 1e38;
  float_sw4 lamax = -1e38;
  float_sw4 cfl2max = -1e38;
  float_sw4 vs2min = 1e38;
  float_sw4 bulkmin = 1e38;
  float_sw4 bulkmax = -1e38;

  float_sw4 acof = dt * dt / (h * h);
  size_t nans = 0;
  for (size_t ind = 0; ind < nijk; ind++) {
    if (std::isnan(rho[ind]) || std::isnan(mu[ind]) || std::isnan(lambda[ind]))
      nans++;

    if (rho[ind] < rhmin) rhmin = rho[ind];
    if (rho[ind] > rhmax) rhmax = rho[ind];
    if (mu[ind] < mumin) mumin = mu[ind];
    if (mu[ind] > mumax) mumax = mu[ind];
    if (2 * mu[ind] + 3 * lambda[ind] < bulkmin)
      bulkmin = 2 * mu[ind] + 3 * lambda[ind];
    if (2 * mu[ind] + 3 * lambda[ind] > bulkmax)
      bulkmax = 2 * mu[ind] + 3 * lambda[ind];
    if (lambda[ind] < lamin) lamin = lambda[ind];
    if (lambda[ind] > lamax) lamax = lambda[ind];
    float_sw4 c1 = acof * (4 * mu[ind] + lambda[ind]) / rho[ind];
    if (c1 > cfl2max) cfl2max = c1;
    float_sw4 c2 = mu[ind] / rho[ind];
    if (c2 < vs2min) vs2min = c2;
    //               c1 = acof*rho(i,j,k)-4*mu(i,j,k)-lambda(i,j,k)
    //               if( c1.lt.c1min )then
    //                  c1min   = c1
    //                  cfl2max = dt*dt*(c1+4*mu(i,j,k)+lambda(i,j,k))/
    //     *                                    (h*h*rho(i,j,k))
    //               endif
    //               c2 = mu(i,j,k)-bcof*rho(i,j,k)
    //               if( c2.lt.c2min )then
    //                  c2min = c2
    //                  vs2min = mu(i,j,k)/rho(i,j,k)
    //               endif
  }
  if (nans > 0)
    cout << "ERROR in EW::checkmtrl. There are " << nans
         << " NaN values in material" << endl;
  limits[0] = rhmin;
  limits[1] = rhmax;
  limits[2] = mumin;
  limits[3] = mumax;
  limits[4] = lamin;
  limits[5] = lamax;
  limits[6] = cfl2max;
  limits[7] = vs2min;
  limits[8] = bulkmin;
  limits[9] = bulkmax;
}

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

    if (topographyExists() && g >= mNumberOfCartesianGrids) {
      // Curvilinear
      projectmtrlc(ifirst, ilast, jfirst, jlast, kfirst, klast, ifirstact,
                   ilastact, jfirstact, jlastact, kfirstact, klastact, rhop,
                   mup, lap, mDt, mMetric[g].c_ptr(), mJ[g].c_ptr(), mCFLmax,
                   vsmin, rhoscale, muscale, lascale, info);
    } else {
      // Cartesian
      projectmtrl(ifirst, ilast, jfirst, jlast, kfirst, klast, ifirstact,
                  ilastact, jfirstact, jlastact, kfirstact, klastact, rhop, mup,
                  lap, mDt, mGridSize[g], mCFLmax, vsmin, rhoscale, muscale,
                  lascale, info);
    }
    if (info != 0)
      cout << "Grid " << g << " info = " << info << " from projectmtrl" << endl;
  }
  material_to_parameters(nmpar, xm, mRho, mMu, mLambda);
}

//-----------------------------------------------------------------------
void correct_material(vector<Sarray>& a_rho, vector<Sarray>& a_mu,
                      vector<Sarray>& a_lambda, float_sw4 vsmin,
                      float_sw4 vpvsminratio) {
  SW4_MARK_FUNCTION;
  // Find (mu,lambda) minimizing (mu-mu0)^2 + (lambda-lamda0)^2
  // with inequality constraints vs >= vsmin, vp/vs >= vpvsminratio
  // Solved as a quadratic programming problem at each grid point.
  float_sw4 alpha = vpvsminratio * vpvsminratio - 2;
  float_sw4 vsmin2 = vsmin * vsmin;
  float_sw4 alfactor = 1 / (1 + alpha * alpha);
  for (int g = 0; g < a_rho.size(); g++) {
    //  int ifirst = m_iStart[g];
    //      int ilast  = m_iEnd[g];
    //      int jfirst = m_jStart[g];
    //      int jlast  = m_jEnd[g];
    //      int kfirst = m_kStart[g];
    //      int klast  = m_kEnd[g];
    //      int ifirstact = m_iStartAct[g];
    //      int ilastact  = m_iEndAct[g];
    //      int jfirstact = m_jStartAct[g];
    //      int jlastact  = m_jEndAct[g];
    //      int kfirstact = m_kStartAct[g];
    //      int klastact  = m_kEndAct[g];

    float_sw4* rhop = a_rho[g].c_ptr();
    float_sw4* mup = a_mu[g].c_ptr();
    float_sw4* lap = a_lambda[g].c_ptr();
    for (int ind = 0; ind <= a_mu[g].npts(); ind++) {
      float_sw4 beta = rhop[ind] * vsmin2;
      float_sw4 kap1 = -lap[ind] + alpha * beta;
      float_sw4 kap2 = alpha * kap1 + beta - mup[ind];
      if (kap1 >= 0 && kap2 >= 0) {
        lap[ind] = alpha * beta;
        mup[ind] = beta;
      } else {
        kap1 = -lap[ind] + alpha * mup[ind];
        kap2 = -mup[ind] + beta;
        if (kap2 >= 0)
          mup[ind] = beta;
        else if (kap1 >= 0) {
          kap1 = (-lap[ind] + alpha * mup[ind]) * alfactor;
          mup[ind] = mup[ind] - alpha * kap1;
          lap[ind] = lap[ind] + kap1;
        }
      }
    }
  }
}

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

    if (topographyExists() && g >= mNumberOfCartesianGrids) {
      // Curvilinear
      projectmtrlc(ifirst, ilast, jfirst, jlast, kfirst, klast, ifirstact,
                   ilastact, jfirstact, jlastact, kfirstact, klastact, rhop,
                   mup, lap, mDt, mMetric[g].c_ptr(), mJ[g].c_ptr(), mCFLmax,
                   vsmin, rhoscale, muscale, lascale, infogrid);
    } else {
      // Cartesian
      projectmtrl(ifirst, ilast, jfirst, jlast, kfirst, klast, ifirstact,
                  ilastact, jfirstact, jlastact, kfirstact, klastact, rhop, mup,
                  lap, mDt, mGridSize[g], mCFLmax, vsmin, rhoscale, muscale,
                  lascale, infogrid);
    }
    if (infogrid != 0) {
      cout << "Grid " << g << " info = " << infogrid << " from projectmtrl"
           << endl;
      if (info == 0) info = infogrid;
    }
  }
}

//-----------------------------------------------------------------------
int EW::check_material(vector<Sarray>& a_rho, vector<Sarray>& a_mu,
                       vector<Sarray>& a_lambda, int& ok, int verbose) {
  SW4_MARK_FUNCTION;
  int err_code = 0;
  ok = 1;
  for (int g = 0; g < mNumberOfGrids; g++) {
    //      int infogrid;
    int numnan = a_rho[g].count_nans();
    bool nansfound = numnan > 0;
    if (numnan > 0)
      std::cout << "check_material found " << numnan << " NaNs in rho "
                << std::endl;
    numnan = a_mu[g].count_nans();
    nansfound = nansfound || numnan > 0;
    if (numnan > 0)
      std::cout << "check_material found " << numnan << " NaNs in mu "
                << std::endl;
    numnan = a_lambda[g].count_nans();
    nansfound = nansfound || numnan > 0;
    if (numnan > 0)
      std::cout << "check_material found " << numnan << " NaNs in lambda "
                << std::endl;
    VERIFY2(!nansfound, "ERROR: check_material found NaNs");

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
    //&kfirst, &klast, 					    &ifirstact, &ilastact, &jfirstact,
    //&jlastact, &kfirstact, 					    &klastact,  rhop, mup, lap, &mDt, mMetric.c_ptr(),
    // mJ.c_ptr(), 					      &mCFLmax, &vsmin, &rhoscale, &muscale, &lascale,
    // &infogrid );
    //      }
    //      else
    //      {
    // Cartesian grid
    checkmtrl(ifirst, ilast, jfirst, jlast, kfirst, klast, rhop, mup, lap, mDt,
              mGridSize[g], limits);
    // Output: limits - vector of 10 elements:
    // (0=rhomin, 1=rhomax, 2=mumin, 3=mumax, 4=lamin, 5=lamax, 6=cfl2max,
    // 7=vs2min, 8=bulkmin, 9=bulkmax)
    //
    // bulk-modulus = 2*mu+3*lambda,  bulk-modulus >=0 gives cp/cs >= sqrt(4/3)

    float_sw4 local[5] = {limits[0], limits[2], limits[4], limits[7],
                          limits[8]};
    float_sw4 global[5];
    MPI_Allreduce(local, global, 5, m_mpifloat, MPI_MIN, m_1d_communicator);
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
    MPI_Allreduce(local, global, 5, m_mpifloat, MPI_MAX, m_1d_communicator);
    limits[1] = global[0];
    limits[3] = global[1];
    limits[5] = global[2];
    limits[6] = global[3];
    limits[9] = global[4];
    if (proc_zero() && verbose > 1) {
      cout << limits[0] << " <=   rho    <= " << limits[1] << " (grid " << g
           << ")" << endl;
      cout << limits[2] << " <=    mu    <= " << limits[3] << " (grid " << g
           << ")" << endl;
      cout << limits[4] << " <=  lambda  <= " << limits[5] << " (grid " << g
           << ")" << endl;

      if (limits[0] <= 0)
        cout << "rho_min = " << limits[0] << " on grid " << g << endl;
      if (limits[2] <= 0)
        cout << "mu_min = " << limits[2] << " on grid " << g << endl;
      if (limits[4] < 0)
        cout << "lambda_min = " << limits[4] << " on grid " << g << endl;
      if (limits[6] <= 0)
        cout << " cfl_max  is imaginary on grid " << g << endl;
      else
        cout << " cfl_max = " << sqrt(limits[6]) << " on grid " << g << endl;
    }  // end if proc_zero && verbose > 1
    ok = ok && (limits[0] > 0 && limits[2] > 0 &&
                limits[6] <= mCFLmax * mCFLmax && limits[8] > 0);
    if (!ok) {
      if (limits[0] <= 0) err_code += 1;
      if (limits[2] <= 0) err_code += 2;
      if (limits[6] > mCFLmax * mCFLmax) err_code += 4;
      if (limits[8] <= 0) err_code += 8;
    }  // end if !ok

  }  // for all grids
  return err_code;
}
