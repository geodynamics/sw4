//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// # Produced at the Lawrence Livermore National Laboratory.
// #
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// #
// # LLNL-CODE-643337
// #
// # All rights reserved.
// #
// # This file is part of SW4, Version: 1.0
// #
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
// #
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991.
// #
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details.
// #
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#include <fenv.h>

#include <cmath>

#include "EW.h"
#include "GridPointSource.h"
#include "Require.h"
#include "mpi.h"

using namespace std;

#include "time_functions.h"

//-----------------------------------------------------------------------
GridPointSource::GridPointSource(float_sw4 frequency, float_sw4 t0, int N,
                                 int M, int L, int G, float_sw4 Fx,
                                 float_sw4 Fy, float_sw4 Fz, timeDep tDep,
                                 int ncyc, float_sw4* pars, int npar,
                                 int* ipars, int nipar, float_sw4* jacobian,
                                 float_sw4* dddp, float_sw4* hess1,
                                 float_sw4* hess2, float_sw4* hess3)
    : mFreq(frequency),
      mT0(t0),
      m_i0(N),
      m_j0(M),
      m_k0(L),
      m_grid(G),
      mTimeDependence(tDep)
#if defined(SOURCE_INVERSION)
      ,
      mNcyc(ncyc) {
#else
{
#endif

  // Copy only pointers, not memory.
  // mPar, mIpar will no longer be changed by this class, they should
  // be set correctly in Source.

  mPar = pars;
  mNpar = npar;
  mIpar = ipars;
  mNipar = nipar;

  mForces[0] = Fx;
  mForces[1] = Fy;
  mForces[2] = Fz;

  initializeTimeFunction();

#if defined(SOURCE_INVERSION)
  m_derivative = -1;

  m_jacobian_known = jacobian != NULL;
  if (jacobian != NULL)
    for (int m = 0; m < 27; m++) m_jacobian[m] = jacobian[m];

  m_hessian_known = false;
  if (dddp != NULL && hess1 != NULL && hess2 != NULL && hess3 != NULL) {
    for (int m = 0; m < 9; m++) {
      m_dddp[m] = dddp[m];
      m_hesspos1[m] = hess1[m];
      m_hesspos2[m] = hess2[m];
      m_hesspos3[m] = hess3[m];
    }
    m_hessian_known = true;
  }
#endif
}

//-----------------------------------------------------------------------
GridPointSource::~GridPointSource() {
  //   if( mPar != NULL )
  //  delete[] mPar;
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE
void GridPointSource::initializeTimeFunction() {
  //   if( mTimeDependence != iDiscrete )
  //      mPar[0] = m_min_exponent;
  // printf("mTimeDepe %d \n",mTimeDependence);
  switch (mTimeDependence) {
    case iRicker:
      mTimeFunc = RickerWavelet;
      mTimeFunc_t = RickerWavelet_t;
      mTimeFunc_tt = RickerWavelet_tt;
      mTimeFunc_ttt = RickerWavelet_ttt;
      mTimeFunc_om = RickerWavelet_om;
      mTimeFunc_omtt = RickerWavelet_omtt;
      break;
    case iGaussian:
      mTimeFunc = Gaussian;
      mTimeFunc_t = Gaussian_t;
      mTimeFunc_tt = Gaussian_tt;
      mTimeFunc_ttt = Gaussian_ttt;
      mTimeFunc_om = Gaussian_om;
      mTimeFunc_omtt = Gaussian_omtt;
      break;
    case iRamp:
      mTimeFunc = Ramp;
      mTimeFunc_t = Ramp_t;
      mTimeFunc_tt = Ramp_tt;
      mTimeFunc_ttt = Ramp_ttt;
      mTimeFunc_om = Ramp_om;
      mTimeFunc_omtt = Ramp_omtt;
      break;
    case iTriangle:
      mTimeFunc = Triangle;
      mTimeFunc_t = Triangle_t;
      mTimeFunc_tt = Triangle_tt;
      mTimeFunc_ttt = Triangle_ttt;
      mTimeFunc_om = Triangle_om;
      mTimeFunc_omtt = Triangle_omtt;
      break;
    case iSawtooth:
      mTimeFunc = Sawtooth;
      mTimeFunc_t = Sawtooth_t;
      mTimeFunc_tt = Sawtooth_tt;
      mTimeFunc_ttt = Sawtooth_ttt;
      mTimeFunc_om = Sawtooth_om;
      mTimeFunc_omtt = Sawtooth_omtt;
      break;
    case iSmoothWave:
      mTimeFunc = SmoothWave;
      mTimeFunc_t = SmoothWave_t;
      mTimeFunc_tt = SmoothWave_tt;
      mTimeFunc_ttt = SmoothWave_ttt;
      mTimeFunc_om = SmoothWave_om;
      mTimeFunc_omtt = SmoothWave_omtt;
      break;
    case iErf:
      mTimeFunc = Erf;
      mTimeFunc_t = Erf_t;
      mTimeFunc_tt = Erf_tt;
      mTimeFunc_ttt = Erf_ttt;
      mTimeFunc_om = Erf_om;
      mTimeFunc_omtt = Erf_omtt;
      break;
    case iVerySmoothBump:
      mTimeFunc = VerySmoothBump;
      mTimeFunc_t = VerySmoothBump_t;
      mTimeFunc_tt = VerySmoothBump_tt;
      mTimeFunc_ttt = VerySmoothBump_ttt;
      mTimeFunc_om = VerySmoothBump_om;
      mTimeFunc_omtt = VerySmoothBump_omtt;
      break;
    case iRickerInt:
      mTimeFunc = RickerInt;
      mTimeFunc_t = RickerInt_t;
      mTimeFunc_tt = RickerInt_tt;
      mTimeFunc_ttt = RickerInt_ttt;
      mTimeFunc_om = RickerInt_om;
      mTimeFunc_omtt = RickerInt_omtt;
      break;
    case iBrune:
      mTimeFunc = Brune;
      mTimeFunc_t = Brune_t;
      mTimeFunc_tt = Brune_tt;
      mTimeFunc_ttt = Brune_ttt;
      mTimeFunc_om = Brune_om;
      mTimeFunc_omtt = Brune_omtt;
      break;
    case iBruneSmoothed:
      mTimeFunc = BruneSmoothed;
      mTimeFunc_t = BruneSmoothed_t;
      mTimeFunc_tt = BruneSmoothed_tt;
      mTimeFunc_ttt = BruneSmoothed_ttt;
      mTimeFunc_om = BruneSmoothed_om;
      mTimeFunc_omtt = BruneSmoothed_omtt;
      break;
    case iDBrune:
      mTimeFunc = DBrune;
      mTimeFunc_t = DBrune_t;
      mTimeFunc_tt = DBrune_tt;
      mTimeFunc_ttt = DBrune_ttt;
      mTimeFunc_om = DBrune_om;
      mTimeFunc_omtt = DBrune_omtt;
      break;
    case iGaussianWindow:
      //      mPar[1] = mNcyc;
      mTimeFunc = GaussianWindow;
      mTimeFunc_t = GaussianWindow_t;
      mTimeFunc_tt = GaussianWindow_tt;
      mTimeFunc_ttt = GaussianWindow_ttt;
      mTimeFunc_om = GaussianWindow_om;
      mTimeFunc_omtt = GaussianWindow_omtt;
      break;
    case iLiu:
      mTimeFunc = Liu;
      mTimeFunc_t = Liu_t;
      mTimeFunc_tt = Liu_tt;
      mTimeFunc_ttt = Liu_ttt;
      mTimeFunc_om = Liu_om;
      mTimeFunc_omtt = Liu_omtt;
      break;
    case iDirac:
      mTimeFunc = Dirac;
      mTimeFunc_t = Dirac_t;
      mTimeFunc_tt = Dirac_tt;
      mTimeFunc_ttt = Dirac_ttt;
      mTimeFunc_om = Dirac_om;
      mTimeFunc_omtt = Dirac_omtt;
      break;
    case iDiscrete:
    case iDiscrete6moments:
    case iDiscrete3forces:
      mTimeFunc = Discrete;
      mTimeFunc_t = Discrete_t;
      mTimeFunc_tt = Discrete_tt;
      mTimeFunc_ttt = Discrete_ttt;
      mTimeFunc_om = Discrete_om;
      mTimeFunc_omtt = Discrete_omtt;
      break;
    case iC6SmoothBump:
      mTimeFunc = C6SmoothBump;
      mTimeFunc_t = C6SmoothBump_t;
      mTimeFunc_tt = C6SmoothBump_tt;
      mTimeFunc_ttt = C6SmoothBump_ttt;
      mTimeFunc_om = C6SmoothBump_om;
      mTimeFunc_omtt = C6SmoothBump_omtt;
      break;
    default:
      // std::cout << "incorrect argument to GridPointSource constructor :
      // default RickerWavelet used " << std::endl;
      mTimeFunc = RickerWavelet;
      mTimeFunc_t = RickerWavelet_t;
      mTimeFunc_tt = RickerWavelet_tt;
      mTimeFunc_ttt = RickerWavelet_ttt;
      mTimeFunc_om = RickerWavelet_om;
      mTimeFunc_omtt = RickerWavelet_omtt;
  }
  // Treat fourth derivatives in special 'switch', because not (yet?)
  // implemented for all time functions
  switch (mTimeDependence) {
    case iVerySmoothBump:
      mTimeFunc_tttt = VerySmoothBump_tttt;
      mTimeFunc_tttom = VerySmoothBump_tttom;
      mTimeFunc_ttomom = VerySmoothBump_ttomom;
      mTimeFunc_tom = VerySmoothBump_tom;
      mTimeFunc_omom = VerySmoothBump_omom;
      break;
    case iGaussian:
      mTimeFunc_tttt = Gaussian_tttt;
      mTimeFunc_tttom = Gaussian_tttom;
      mTimeFunc_ttomom = Gaussian_ttomom;
      mTimeFunc_tom = Gaussian_tom;
      mTimeFunc_omom = Gaussian_omom;
      break;
    case iDirac:
      mTimeFunc_tttt = Dirac_tttt;
      mTimeFunc_tttom = Dirac_tttom;
      mTimeFunc_ttomom = Dirac_ttomom;
      mTimeFunc_tom = Dirac_tom;
      mTimeFunc_omom = Dirac_omom;
      break;
    case iDiscrete:
    case iDiscrete6moments:
    case iDiscrete3forces:
      mTimeFunc_tttt = Discrete_tttt;
      mTimeFunc_tttom = Discrete_tttom;
      mTimeFunc_ttomom = Discrete_ttomom;
      mTimeFunc_tom = Discrete_tom;
      mTimeFunc_omom = Discrete_omom;
      break;
    default:
      // tmp
      // std::cout << "High derivatives not implemented for time fuction:" <<
      // mTimeDependence <<
      //   " default Gaussian used for tttt, ttt-omega derivatives, etc " <<
      //   std::endl;
      mTimeFunc_tttt = Gaussian_tttt;
      mTimeFunc_tttom = Gaussian_tttom;
      mTimeFunc_ttomom = Gaussian_ttomom;
      mTimeFunc_tom = Gaussian_tom;
      mTimeFunc_omom = Gaussian_omom;
  }
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE
void GridPointSource::getFxyz(float_sw4 t, float_sw4* fxyz) {
  float_sw4 afun, afunv[6];
  //printf("mTimeDependence =%d %d\n",mTimeDependence,iRicker);

  if (mNpar==0) {
    std::cout<<"FXYZ "<<mNpar<<" "<<mNipar<<"\n"<<std::flush;
    afun = 0.0;
    return;
  }
  if (mTimeDependence != iDiscrete6moments &&
      mTimeDependence != iDiscrete3forces)
    afun = mTimeFunc(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);

  else if (mTimeDependence == iDiscrete6moments) {
    int npts = mIpar[0];
    int size = 6 * (npts - 1) + 1;
    size_t pos = 0;
    afunv[0] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[1] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[2] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[3] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[4] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[5] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
  } else if (mTimeDependence == iDiscrete3forces) {
    int npts = mIpar[0];
    int size = 6 * (npts - 1) + 1;
    size_t pos = 0;
    afunv[0] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[1] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[2] = mTimeFunc(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
  }
  //printf("AFUN]=%f %f %f %f\n",afun,mForces[0],mForces[1],mForces[2]);
  int lm_derivative;
#if defined(SOURCE_INVERSION)
  lm_derivative = m_derivative;
#else
  lm_derivative = -1;
#endif
  if (lm_derivative == -1) {
    if (mTimeDependence != iDiscrete6moments &&
        mTimeDependence != iDiscrete3forces) {
      fxyz[0] = mForces[0] * afun;
      fxyz[1] = mForces[1] * afun;
      fxyz[2] = mForces[2] * afun;
    } else if (mTimeDependence == iDiscrete6moments) {
      fxyz[0] =
          mForces[0] * afunv[0] + mForces[1] * afunv[1] + mForces[2] * afunv[2];
      fxyz[1] =
          mForces[0] * afunv[1] + mForces[1] * afunv[3] + mForces[2] * afunv[4];
      fxyz[2] =
          mForces[0] * afunv[2] + mForces[1] * afunv[4] + mForces[2] * afunv[5];
    } else if (mTimeDependence == iDiscrete3forces) {
      fxyz[0] = mForces[0] * afunv[0];
      fxyz[1] = mForces[0] * afunv[1];
      fxyz[2] = mForces[0] * afunv[2];
    }
  } else if (lm_derivative >= 0 && lm_derivative <= 8) {
#if defined(SOURCE_INVERSION)
    fxyz[0] = m_jacobian[m_derivative * 3] * afun;
    fxyz[1] = m_jacobian[m_derivative * 3 + 1] * afun;
    fxyz[2] = m_jacobian[m_derivative * 3 + 2] * afun;
#endif
  } else if (lm_derivative == 9) {
    afun = -mTimeFunc_t(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    fxyz[0] = mForces[0] * afun;
    fxyz[1] = mForces[1] * afun;
    fxyz[2] = mForces[2] * afun;
  } else if (lm_derivative == 10) {
    afun = mTimeFunc_om(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    fxyz[0] = mForces[0] * afun;
    fxyz[1] = mForces[1] * afun;
    fxyz[2] = mForces[2] * afun;
  } else if (lm_derivative == 11) {
#if defined(SOURCE_INVERSION)
    fxyz[0] = fxyz[1] = fxyz[2] = 0;
    int i;
    for (i = 0; i < 9; i++) {
      fxyz[0] += afun * m_jacobian[i * 3] * m_dir[i];
      fxyz[1] += afun * m_jacobian[i * 3 + 1] * m_dir[i];
      fxyz[2] += afun * m_jacobian[i * 3 + 2] * m_dir[i];
    }
    i = 9;
    afun = -mTimeFunc_t(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    fxyz[0] += afun * mForces[0] * m_dir[i];
    fxyz[1] += afun * mForces[1] * m_dir[i];
    fxyz[2] += afun * mForces[2] * m_dir[i];
    i = 10;
    afun = mTimeFunc_om(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    fxyz[0] += afun * mForces[0] * m_dir[i];
    fxyz[1] += afun * mForces[1] * m_dir[i];
    fxyz[2] += afun * mForces[2] * m_dir[i];
#endif
  }
}

//-----------------------------------------------------------------------
void GridPointSource::getFxyz_notime(float_sw4* fxyz) const {
  // For source spatial discretization testing
  fxyz[0] = mForces[0];
  fxyz[1] = mForces[1];
  fxyz[2] = mForces[2];
}

//-----------------------------------------------------------------------
RAJA_HOST_DEVICE
void GridPointSource::getFxyztt(float_sw4 t, float_sw4* fxyz) const {
  float_sw4 afun, afunv[6];
  if (mTimeDependence != iDiscrete6moments &&
      mTimeDependence != iDiscrete3forces)
    afun = mTimeFunc_tt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
  else if (mTimeDependence == iDiscrete6moments) {
    int npts = mIpar[0];
    int size = 6 * (npts - 1) + 1;
    size_t pos = 0;
    afunv[0] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[1] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[2] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[3] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[4] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[5] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
  } else if (mTimeDependence == iDiscrete3forces) {
    int npts = mIpar[0];
    int size = 6 * (npts - 1) + 1;
    size_t pos = 0;
    afunv[0] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[1] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
    pos += size;
    afunv[2] = mTimeFunc_tt(mFreq, t - mT0, mPar + pos, mNpar, mIpar, mNipar);
  }

  //  float_sw4 afun = mTimeFunc_tt(mFreq,t-mT0,mPar, mNpar, mIpar, mNipar);
#if defined(SOURCE_INVERSION)
  if (m_derivative == -1) {
    if (mTimeDependence != iDiscrete6moments &&
        mTimeDependence != iDiscrete3forces) {
      fxyz[0] = mForces[0] * afun;
      fxyz[1] = mForces[1] * afun;
      fxyz[2] = mForces[2] * afun;
    } else if (mTimeDependence == iDiscrete6moments) {
      fxyz[0] =
          mForces[0] * afunv[0] + mForces[1] * afunv[1] + mForces[2] * afunv[2];
      fxyz[1] =
          mForces[0] * afunv[1] + mForces[1] * afunv[3] + mForces[2] * afunv[4];
      fxyz[2] =
          mForces[0] * afunv[2] + mForces[1] * afunv[4] + mForces[2] * afunv[5];
    } else if (mTimeDependence == iDiscrete3forces) {
      fxyz[0] = mForces[0] * afunv[0];
      fxyz[1] = mForces[0] * afunv[1];
      fxyz[2] = mForces[0] * afunv[2];
    }
  } else if (m_derivative >= 0 && m_derivative <= 8) {
    fxyz[0] = m_jacobian[m_derivative * 3] * afun;
    fxyz[1] = m_jacobian[m_derivative * 3 + 1] * afun;
    fxyz[2] = m_jacobian[m_derivative * 3 + 2] * afun;
  } else if (m_derivative == 9) {
    afun = -mTimeFunc_ttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    fxyz[0] = mForces[0] * afun;
    fxyz[1] = mForces[1] * afun;
    fxyz[2] = mForces[2] * afun;
  } else if (m_derivative == 10) {
    afun = mTimeFunc_omtt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    fxyz[0] = mForces[0] * afun;
    fxyz[1] = mForces[1] * afun;
    fxyz[2] = mForces[2] * afun;
  } else if (m_derivative == 11) {
    fxyz[0] = fxyz[1] = fxyz[2] = 0;
    int i;
    for (i = 0; i < 9; i++) {
      fxyz[0] += afun * m_jacobian[i * 3] * m_dir[i];
      fxyz[1] += afun * m_jacobian[i * 3 + 1] * m_dir[i];
      fxyz[2] += afun * m_jacobian[i * 3 + 2] * m_dir[i];
    }
    i = 9;
    afun = -mTimeFunc_ttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    fxyz[0] += afun * mForces[0] * m_dir[i];
    fxyz[1] += afun * mForces[1] * m_dir[i];
    fxyz[2] += afun * mForces[2] * m_dir[i];
    i = 10;
    afun = mTimeFunc_omtt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    fxyz[0] += afun * mForces[0] * m_dir[i];
    fxyz[1] += afun * mForces[1] * m_dir[i];
    fxyz[2] += afun * mForces[2] * m_dir[i];
  }
#else
  if (mTimeDependence != iDiscrete6moments) {
    fxyz[0] = mForces[0] * afun;
    fxyz[1] = mForces[1] * afun;
    fxyz[2] = mForces[2] * afun;
  } else {
    fxyz[0] =
        mForces[0] * afunv[0] + mForces[1] * afunv[1] + mForces[2] * afunv[2];
    fxyz[1] =
        mForces[0] * afunv[1] + mForces[1] * afunv[3] + mForces[2] * afunv[4];
    fxyz[2] =
        mForces[0] * afunv[2] + mForces[1] * afunv[4] + mForces[2] * afunv[5];
  }
#endif
}

//-----------------------------------------------------------------------
void GridPointSource::set_derivative(int der, const float_sw4 dir[11]) {
#if defined(SOURCE_INVERSION)
  if (der >= 0 && der <= 11) m_derivative = der;
  for (int i = 0; i < 11; i++) m_dir[i] = dir[i];
#endif
}

//-----------------------------------------------------------------------
void GridPointSource::set_noderivative() {
#if defined(SOURCE_INVERSION)
  m_derivative = -1;
#endif
}

//-----------------------------------------------------------------------
void GridPointSource::limitFrequency(float_sw4 max_freq) {
  if (mFreq > max_freq) mFreq = max_freq;
}

//-----------------------------------------------------------------------
float_sw4 GridPointSource::getTimeFunc(float_sw4 t) const {
  return mTimeFunc(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
float_sw4 GridPointSource::evalTimeFunc_t(float_sw4 t) const {
  return mTimeFunc_t(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
float_sw4 GridPointSource::evalTimeFunc_tt(float_sw4 t) const {
  return mTimeFunc_tt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
float_sw4 GridPointSource::evalTimeFunc_ttt(float_sw4 t) const {
  return mTimeFunc_ttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
float_sw4 GridPointSource::evalTimeFunc_tttt(float_sw4 t) const {
  return mTimeFunc_tttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
}

//-----------------------------------------------------------------------
ostream& operator<<(ostream& output, const GridPointSource& s) {
  output << "GridPointSource at (i,j,k) = " << s.m_i0 << "," << s.m_j0 << ","
         << s.m_k0 << " in grid no " << s.m_grid << endl;
  //   output << "   Strength " << s.mAmp;
  output << " Fx Fy Fz = " << s.mForces[0] << " " << s.mForces[1] << " "
         << s.mForces[2] << endl;
  output << " freq = " << s.mFreq << " t0 = " << s.mT0 << endl;
  output << " npar = " << s.mNpar << " nipar= " << s.mNipar << endl;
  if (s.mNpar > 0) output << " mpar[0] = " << s.mPar[0];
  if (s.mNipar > 0) output << " mipar[0] = " << s.mIpar[0];
  if (s.mNpar > 0 || s.mNipar > 0) cout << endl;

  return output;
}

//-----------------------------------------------------------------------
void GridPointSource::add_to_gradient(std::vector<Sarray>& kappa,
                                      std::vector<Sarray>& eta, float_sw4 t,
                                      float_sw4 dt, float_sw4 gradient[11],
                                      std::vector<float_sw4>& h,
                                      vector<Sarray>& Jac,
                                      bool topography_exists) {
#if defined(SOURCE_INVERSION)
  if (m_jacobian_known) {
    float_sw4 normwgh[4] = {17.0 / 48.0, 59.0 / 48.0, 43.0 / 48.0, 49.0 / 48.0};
    float_sw4 dt2o12 = dt * dt / 12.0;
    float_sw4 g0 = mTimeFunc(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 g =
        g0 + dt2o12 * mTimeFunc_tt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);

    // save some work by accessing array elements only once:
    float_sw4 kap1 = kappa[m_grid](1, m_i0, m_j0, m_k0);
    float_sw4 kap2 = kappa[m_grid](2, m_i0, m_j0, m_k0);
    float_sw4 kap3 = kappa[m_grid](3, m_i0, m_j0, m_k0);
    float_sw4 eta1 = eta[m_grid](1, m_i0, m_j0, m_k0);
    float_sw4 eta2 = eta[m_grid](2, m_i0, m_j0, m_k0);
    float_sw4 eta3 = eta[m_grid](3, m_i0, m_j0, m_k0);
    float_sw4 h3 = h[m_grid] * h[m_grid] * h[m_grid];
    if (topography_exists && m_grid == h.size() - 1)
      h3 = Jac[m_grid](m_i0, m_j0, m_k0);
    //      float_sw4 h3 = 1.0;
    if (1 <= m_k0 && m_k0 <= 4) h3 *= normwgh[m_k0 - 1];

    // derivative wrt. position (m=0,1,2) and moment tensor components
    // (m=3,..,8)
    for (int m = 0; m < 9; m++) {
      gradient[m] -= g *
                     (kap1 * m_jacobian[3 * m] + kap2 * m_jacobian[3 * m + 1] +
                      kap3 * m_jacobian[3 * m + 2]) *
                     h3;
      gradient[m] -= dt2o12 * g0 *
                     (eta1 * m_jacobian[3 * m] + eta2 * m_jacobian[3 * m + 1] +
                      eta3 * m_jacobian[3 * m + 2]) *
                     h3;
    }

    // derivative wrt. (t0, freq)
    float_sw4 dgt0 = -mTimeFunc_t(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 dgom = mTimeFunc_om(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    gradient[9] -= dt2o12 * dgt0 *
                   (eta1 * mForces[0] + eta2 * mForces[1] + eta3 * mForces[2]) *
                   h3;
    gradient[10] -=
        dt2o12 * dgom *
        (eta1 * mForces[0] + eta2 * mForces[1] + eta3 * mForces[2]) * h3;

    dgt0 = dgt0 -
           dt2o12 * mTimeFunc_ttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    dgom = dgom +
           dt2o12 * mTimeFunc_omtt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    gradient[9] -=
        dgt0 * (kap1 * mForces[0] + kap2 * mForces[1] + kap3 * mForces[2]) * h3;
    gradient[10] -=
        dgom * (kap1 * mForces[0] + kap2 * mForces[1] + kap3 * mForces[2]) * h3;
  }
#else
  std::cerr << "This routine is disabled. Compile with -DSOURCE_INVERSION=1\n"
            << std::flush;
  abort();
#endif
}

//-----------------------------------------------------------------------
void GridPointSource::add_to_hessian(std::vector<Sarray>& kappa,
                                     std::vector<Sarray>& eta, float_sw4 t,
                                     float_sw4 dt, float_sw4 hessian[121],
                                     std::vector<float_sw4>& h)
// Add upper part of symmetric matrix
{
#if defined(SOURCE_INVERSION)
  if (m_hessian_known && m_jacobian_known) {
    float_sw4 normwgh[4] = {17.0 / 48.0, 59.0 / 48.0, 43.0 / 48.0, 49.0 / 48.0};
    float_sw4 dt2o12 = dt * dt / 12.0;
    float_sw4 g0 = mTimeFunc(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 g =
        g0 + dt2o12 * mTimeFunc_tt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);

    // save some work by accessing array elements only once:
    float_sw4 kap1 = kappa[m_grid](1, m_i0, m_j0, m_k0);
    float_sw4 kap2 = kappa[m_grid](2, m_i0, m_j0, m_k0);
    float_sw4 kap3 = kappa[m_grid](3, m_i0, m_j0, m_k0);
    float_sw4 eta1 = eta[m_grid](1, m_i0, m_j0, m_k0);
    float_sw4 eta2 = eta[m_grid](2, m_i0, m_j0, m_k0);
    float_sw4 eta3 = eta[m_grid](3, m_i0, m_j0, m_k0);
    float_sw4 h3 = h[m_grid] * h[m_grid] * h[m_grid];

    if (1 <= m_k0 && m_k0 <= 4) h3 *= normwgh[m_k0 - 1];

    float_sw4 c1 = g * h3;
    float_sw4 c2 = g0 * dt2o12 * h3;

    // (pos,pos)
    for (int m = 0; m < 3; m++)
      for (int j = m; j < 3; j++) {
        hessian[m + 11 * j] -= (kap1 * c1 + eta1 * c2) * m_hesspos1[m + 3 * j] +
                               (kap2 * c1 + eta2 * c2) * m_hesspos2[m + 3 * j] +
                               (kap3 * c1 + eta3 * c2) * m_hesspos3[m + 3 * j];
      }
    // (pos,mij)
    for (int m = 0; m < 3; m++) {
      int j = 3;
      hessian[m + 11 * j] -= (c1 * kap1 + eta1 * c2) * m_dddp[m];
      j = 4;
      hessian[m + 11 * j] -= (c1 * kap1 + eta1 * c2) * m_dddp[m + 3] +
                             (c1 * kap2 + eta2 * c2) * m_dddp[m];
      j = 5;
      hessian[m + 11 * j] -= (c1 * kap1 + eta1 * c2) * m_dddp[m + 6] +
                             (c1 * kap3 + eta3 * c2) * m_dddp[m];
      j = 6;
      hessian[m + 11 * j] -= (c1 * kap2 + eta2 * c2) * m_dddp[m + 3];
      j = 7;
      hessian[m + 11 * j] -= (c1 * kap2 + eta2 * c2) * m_dddp[m + 6] +
                             (c1 * kap3 + eta3 * c2) * m_dddp[m + 3];
      j = 8;
      hessian[m + 11 * j] -= (c1 * kap3 + eta3 * c2) * m_dddp[m + 6];
    }
    // (pos,t0)
    float_sw4 dgt0 = -mTimeFunc_t(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 dgom = mTimeFunc_om(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);

    float_sw4 c2t0 = dgt0 * dt2o12;
    float_sw4 c2om0 = dgom * dt2o12;

    dgt0 = dgt0 -
           dt2o12 * mTimeFunc_ttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    dgom = dgom +
           dt2o12 * mTimeFunc_omtt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 c1t0 = dgt0;
    float_sw4 c1om0 = dgom;

    // define c1t0,c2t0, c1om0, c2om0
    for (int m = 0; m < 9; m++) {
      int j = 9;
      hessian[m + 11 * j] -=
          ((c1t0 * kap1 + c2t0 * eta1) * m_jacobian[3 * m] +
           (c1t0 * kap2 + c2t0 * eta2) * m_jacobian[3 * m + 1] +
           (c1t0 * kap3 + c2t0 * eta3) * m_jacobian[3 * m + 2]) *
          h3;
      j = 10;
      hessian[m + 11 * j] -=
          ((c1om0 * kap1 + c2om0 * eta1) * m_jacobian[3 * m] +
           (c1om0 * kap2 + c2om0 * eta2) * m_jacobian[3 * m + 1] +
           (c1om0 * kap3 + c2om0 * eta3) * m_jacobian[3 * m + 2]) *
          h3;
    }
    float_sw4 cmfact0 =
        (kap1 * mForces[0] + kap2 * mForces[1] + kap3 * mForces[2]);
    float_sw4 cmfact =
        ((kap1 * mForces[0] + kap2 * mForces[1] + kap3 * mForces[2]) +
         dt2o12 * (eta1 * mForces[0] + eta2 * mForces[1] + eta3 * mForces[2]));
    // Second derivatives of time function
    float_sw4 d2gdt02 =
        mTimeFunc_tt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 d2gdt0dom =
        -mTimeFunc_tom(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 d2gdomdom =
        mTimeFunc_omom(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);

    float_sw4 dgdttt0t0 =
        dt2o12 * mTimeFunc_tttt(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 dgdttt0om =
        -dt2o12 * mTimeFunc_tttom(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);
    float_sw4 dgdttomom =
        dt2o12 * mTimeFunc_ttomom(mFreq, t - mT0, mPar, mNpar, mIpar, mNipar);

    int m = 9;
    int j = 9;
    hessian[m + 11 * j] -= (cmfact * d2gdt02 + cmfact0 * dgdttt0t0) * h3;
    j = 10;
    hessian[m + 11 * j] -= (cmfact * d2gdt0dom + cmfact0 * dgdttt0om) * h3;
    m = 10;
    j = 10;
    hessian[m + 11 * j] -= (cmfact * d2gdomdom + cmfact0 * dgdttomom) * h3;
  }
#endif
}

//-----------------------------------------------------------------------
void GridPointSource::print_info() const {
  cout << "--------------------------------------------------------------------"
          "---"
       << endl;
  cout << " position " << m_i0 << " " << m_j0 << " " << m_k0 << endl;
  cout << "Forces = " << mForces[0] << " " << mForces[1] << " " << mForces[2]
       << endl;
  cout << " jac = \n";
#if defined(SOURCE_INVERSION)
  for (int c = 0; c < 3; c++) {
    for (int m = 0; m < 9; m++) cout << m_jacobian[c + 3 * m] << " ";
    cout << endl;
  }
#endif
  cout << "Time dep " << mTimeDependence << endl;
  cout << "--------------------------------------------------------------------"
          "---"
       << endl;
}

//-----------------------------------------------------------------------
void GridPointSource::set_sort_key(size_t key) { m_key = key; }
