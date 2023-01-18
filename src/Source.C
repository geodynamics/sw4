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
#include "F77_FUNC.h"
#include "Filter.h"
#include "GridPointSource.h"
#include "Mspace.h"
#include "Qspline.h"
#include "Require.h"
#include "Source.h"
#include "caliper.h"
#include "mpi.h"
#include "time_functions.h"
//#include <fcntl.h>
//#include <unistd.h>

using namespace std;

#define SQR(x) ((x) * (x))

//-----------------------------------------------------------------------
// Constructor,
//
//    ncyc is only used in the 'GaussianWindow' time function
//
//    pars is only used in the 'Discrete' time function
//        when pars[1],..pars[npts] should contain the discrete function on a
//        uniform grid with spacing dt=1/freq, and pars[0] is the first time,
//        thus the grid is
//            t_k = pars[0] + dt*k, k=0,1,..,npts-1
//    ipar should have size 1, with ipar[0] containing npts.
//
//    When the source time function is not 'Discrete', the input pars and ipars
//    will not be used.
//
Source::Source(EW* a_ew, float_sw4 frequency, float_sw4 t0, float_sw4 x0,
               float_sw4 y0, float_sw4 z0, float_sw4 Mxx, float_sw4 Mxy,
               float_sw4 Mxz, float_sw4 Myy, float_sw4 Myz, float_sw4 Mzz,
               timeDep tDep, const char* name, bool topodepth, int ncyc,
               float_sw4* pars, int npar, int* ipars, int nipar,
               bool correctForMu)
    : mIsMomentSource(true),
      mFreq(frequency),
      mT0(t0),
      mX0(x0),
      mY0(y0),
      mZ0(z0),
      m_zTopo(-1e38),
      mIgnore(false),
      //  mGridPointSet(false),
      mTimeDependence(tDep),
      mNcyc(ncyc),
      m_zRelativeToTopography(topodepth),
      mShearModulusFactor(correctForMu),
      m_derivative(-1),
      m_is_filtered(false),
      m_myPoint(false),
      m_timeFuncIsReady(false) {
 
  mForces.resize(6);
  mForces[0] = Mxx;
  mForces[1] = Mxy;
  mForces[2] = Mxz;
  mForces[3] = Myy;
  mForces[4] = Myz;
  mForces[5] = Mzz;
  mName = name;

  mNpar = npar;
  if (mNpar > 0) {
    mPar = SW4_NEW(Space::Managed, float_sw4[mNpar]);
    for (int i = 0; i < mNpar; i++) mPar[i] = pars[i];
  } else {
    mNpar = 2;
    mPar = SW4_NEW(Space::Managed, float_sw4[2]);
  }
  mNipar = nipar;
  if (mNipar > 0) {
    mIpar = SW4_NEW(Space::Managed, int[mNipar]);
    for (int i = 0; i < mNipar; i++) mIpar[i] = ipars[i];
  } else {
    mNipar = 1;
    mIpar = SW4_NEW(Space::Managed, int[1]);
  }

  // if( mTimeDependence == iDiscrete || mTimeDependence == iDiscrete6moments )
  //    spline_interpolation();
  // else
  if (mTimeDependence != iDiscrete && mTimeDependence != iDiscrete6moments &&
      mTimeDependence != iDiscrete3forces)  // not sure about iDiscrete6moments
  {
    mPar[0] = find_min_exponent();
    mPar[1] = mNcyc;
  }

  //   a_ew->computeNearestGridPoint2(m_i0,m_j0,m_k0,m_grid,mX0,mY0,mZ0);

  // Correct source location for discrepancy between raw and smoothed topography
  correct_Z_level(a_ew);
  compute_grid_point(a_ew);
  if (a_ew->getVerbosity() >= 3 && a_ew->proc_zero()) {
    printf(
        "Moment source at x=%e, y=%e, z=%e is centered at grid point i=%d, "
        "j=%d, k=%d, in grid=%d\n",
        mX0, mY0, mZ0, m_i0, m_j0, m_k0, m_grid);
  }
}

//-----------------------------------------------------------------------
Source::Source(EW* a_ew, float_sw4 frequency, float_sw4 t0, float_sw4 x0,
               float_sw4 y0, float_sw4 z0, float_sw4 Fx, float_sw4 Fy,
               float_sw4 Fz, timeDep tDep, const char* name, bool topodepth,
               int ncyc, float_sw4* pars, int npar, int* ipars, int nipar,
               bool correctForMu)
    : mIsMomentSource(false),
      mFreq(frequency),
      mT0(t0),
      mX0(x0),
      mY0(y0),
      mZ0(z0),
      m_zTopo(-1e38),
      mIgnore(false),
      //  mGridPointSet(false),
      mTimeDependence(tDep),
      mNcyc(ncyc),
      m_zRelativeToTopography(topodepth),
      m_derivative(-1),
      m_is_filtered(false),
      mShearModulusFactor(correctForMu),
      m_myPoint(false),
      m_timeFuncIsReady(false) {
  //printf("MXX MZY MXZ %f %f %f \n",Fx,Fy,Fz);
  mForces.resize(3);
  mForces[0] = Fx;
  mForces[1] = Fy;
  mForces[2] = Fz;
  mName = name;

  mNpar = npar;
  if (mNpar > 0) {
    mPar = SW4_NEW(Space::Managed, float_sw4[mNpar]);
    for (int i = 0; i < mNpar; i++) mPar[i] = pars[i];
  } else {
    mNpar = 2;
    mPar = SW4_NEW(Space::Managed, float_sw4[2]);
  }

  mNipar = nipar;
  if (mNipar > 0) {
    mIpar = new int[mNipar];
    for (int i = 0; i < mNipar; i++) mIpar[i] = ipars[i];
  } else {
    mNipar = 1;
    mIpar = new int[1];
  }

  if (mTimeDependence == iDiscrete || mTimeDependence == iDiscrete6moments ||
      mTimeDependence == iDiscrete3forces)
    spline_interpolation();
  else {
    mPar[0] = find_min_exponent();
    mPar[1] = mNcyc;
  }

  //  a_ew->computeNearestGridPoint2(m_i0,m_j0,m_k0,m_grid,mX0,mY0,mZ0);

  // Correct source location for discrepancy between raw and smoothed topography
  correct_Z_level(a_ew);  // also sets the ignore flag for sources that are
                          // above the topography
  compute_grid_point(a_ew);
}

//-----------------------------------------------------------------------
Source::Source() {
  mNpar = 0;
  mNipar = 0;
}

//-----------------------------------------------------------------------
Source::~Source() {
  if (mNpar > 0) ::operator delete[](mPar, Space::Managed);
  if (mNipar > 0) ::operator delete[](mIpar, Space::Managed);
}

//-----------------------------------------------------------------------
float_sw4 Source::getX0() const { return mX0; }

//-----------------------------------------------------------------------
float_sw4 Source::getY0() const { return mY0; }

//-----------------------------------------------------------------------
float_sw4 Source::getZ0() const { return mZ0; }

//-----------------------------------------------------------------------
float_sw4 Source::getDepth() const { return mZ0 - m_zTopo; }

//-----------------------------------------------------------------------
float_sw4 Source::getOffset() const { return mT0; }

//-----------------------------------------------------------------------
float_sw4 Source::getFrequency() const { return mFreq; }

//-----------------------------------------------------------------------
void Source::setMaxFrequency(float_sw4 max_freq) {
  if (mFreq > max_freq) mFreq = max_freq;
}

//-----------------------------------------------------------------------
bool Source::isMomentSource() const { return mIsMomentSource; }

//-----------------------------------------------------------------------
void Source::getForces(float_sw4& fx, float_sw4& fy, float_sw4& fz) const {
  SW4_MARK_FUNCTION;
  if (!mIsMomentSource) {
    fx = mForces[0];
    fy = mForces[1];
    fz = mForces[2];
  } else
    fx = fy = fz = 0;
}

//-----------------------------------------------------------------------
void Source::getMoments(float_sw4& mxx, float_sw4& mxy, float_sw4& mxz,
                        float_sw4& myy, float_sw4& myz, float_sw4& mzz) const {
  SW4_MARK_FUNCTION;
  if (mIsMomentSource) {
    mxx = mForces[0];
    mxy = mForces[1];
    mxz = mForces[2];
    myy = mForces[3];
    myz = mForces[4];
    mzz = mForces[5];
  } else
    mxx = mxy = mxz = myy = myz = mzz = 0;
}

//-----------------------------------------------------------------------
void Source::setMoments(float_sw4 mxx, float_sw4 mxy, float_sw4 mxz,
                        float_sw4 myy, float_sw4 myz, float_sw4 mzz) {
  if (mIsMomentSource) {
    mForces[0] = mxx;
    mForces[1] = mxy;
    mForces[2] = mxz;
    mForces[3] = myy;
    mForces[4] = myz;
    mForces[5] = mzz;
  } else {
    mForces[0] = mxx;
    mForces[1] = myy;
    mForces[2] = mzz;
  }
}

//-----------------------------------------------------------------------
float_sw4 Source::getAmplitude() const {
  float_sw4 amplitude = 0;
  if (mIsMomentSource) {
    float_sw4 msqr = 0;
    for (int q = 0; q < 6; q++) msqr += SQR(mForces[q]);
    //    amplitude = mAmp*sqrt(msqr/2.);
    msqr += SQR(mForces[1]) + SQR(mForces[2]) + SQR(mForces[4]);
    amplitude = sqrt(0.5 * msqr);
  } else {
    float_sw4 fsqr = 0;
    for (int q = 0; q < 3; q++) fsqr += SQR(mForces[q]);
    //    amplitude = mAmp*sqrt(fsqr);
    amplitude = sqrt(fsqr);
  }
  return amplitude;
}

//-----------------------------------------------------------------------
ostream& operator<<(ostream& output, const Source& s) {
  output << s.mName << (s.isMomentSource() ? " moment" : " force")
         << " source term" << endl;
  output << "   Location (X,Y,Z) = " << s.mX0 << "," << s.mY0 << "," << s.mZ0
         << " at grid point " << s.m_i0 << " " << s.m_j0 << " " << s.m_k0
         << " in grid no " << s.m_grid << endl;
  output << "   Strength " << s.getAmplitude();
  output << "   t0 = " << s.mT0 << " freq = " << s.mFreq << endl;
  if (s.mIsMomentSource) {
    output << " Mxx Mxy Myy Mxz Myz Mzz = " << s.mForces[0] << " "
           << s.mForces[1] << " " << s.mForces[3] << " " << s.mForces[2] << " "
           << s.mForces[4] << " " << s.mForces[5] << endl;
  } else {
    output << " Fx Fy Fz = " << s.mForces[0] << " " << s.mForces[1] << " "
           << s.mForces[2] << endl;
  }
  return output;
}

//-----------------------------------------------------------------------
void Source::limit_frequency(int ppw, float_sw4 minvsoh) {
  float_sw4 freqlim = minvsoh / (ppw);

  if (mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed ||
      mTimeDependence == iDBrune || mTimeDependence == iGaussian ||
      mTimeDependence == iErf || mTimeDependence == iVerySmoothBump ||
      mTimeDependence == iSmoothWave || mTimeDependence == iLiu ||
      mTimeDependence == iC6SmoothBump) {
    if (mFreq > 2 * M_PI * freqlim) mFreq = 2 * M_PI * freqlim;
  } else {
    if (mFreq > freqlim) mFreq = freqlim;
  }
}

//-----------------------------------------------------------------------
float_sw4 Source::compute_t0_increase(float_sw4 t0_min) const {
  // Gaussian, GaussianInt=Erf, Ricker and RickerInt are all centered around mT0
  if (mTimeDependence == iGaussian || mTimeDependence == iErf)
    return t0_min + 6.0 / mFreq -
           mT0;  // translating these by at least 6*sigma = 6/freq
  else if (mTimeDependence == iRicker || mTimeDependence == iRickerInt)
    return t0_min + 1.9 / mFreq - mT0;  // 1.9 ?
  else
    return t0_min - mT0;  // the rest of the time functions are zero for t<mT0
}

//-----------------------------------------------------------------------
void Source::adjust_t0(float_sw4 dt0) {
  if (dt0 > 0 && !m_is_filtered) mT0 += dt0;
}

//-----------------------------------------------------------------------
float_sw4 Source::dt_to_resolve(int ppw) const {
  float_sw4 dt_resolved = 0;
  if (mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed ||
      mTimeDependence == iDBrune) {
    const float_sw4 t95 = 4.744 / mFreq;
    dt_resolved = t95 / ppw;
  } else {
  }
  return dt_resolved;
}

//-----------------------------------------------------------------------
int Source::ppw_to_resolve(float_sw4 dt) const {
  int ppw = 1;
  if (mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed ||
      mTimeDependence == iDBrune) {
    const float_sw4 t95 = 4.744 / mFreq;
    ppw = static_cast<int>(t95 / dt);
  } else {
  }
  return ppw;
}

//-----------------------------------------------------------------------
void Source::set_derivative(int der) {
  if (der >= 0 && der <= 10) m_derivative = der;
}

//-----------------------------------------------------------------------
void Source::set_noderivative() { m_derivative = -1; }

//-----------------------------------------------------------------------
void Source::set_dirderivative(float_sw4 dir[11]) {
  for (int i = 0; i < 11; i++) m_dir[i] = dir[i];
  m_derivative = 11;
}

//-----------------------------------------------------------------------
void Source::set_parameters(float_sw4 x[11]) {
  if (mIsMomentSource) {
    mX0 = x[0];
    mY0 = x[1];
    mZ0 = x[2];
    mForces[0] = x[3];
    mForces[1] = x[4];
    mForces[2] = x[5];
    mForces[3] = x[6];
    mForces[4] = x[7];
    mForces[5] = x[8];
    mT0 = x[9];
    mFreq = x[10];
  } else
    cout << "Error in Source::set_parameters(), "
         << "function only implemented for moment sources" << endl;
}

//-----------------------------------------------------------------------
void Source::get_parameters(float_sw4 x[11]) const {
  if (mIsMomentSource) {
    x[0] = mX0;
    x[1] = mY0;
    x[2] = mZ0;
    x[3] = mForces[0];
    x[4] = mForces[1];
    x[5] = mForces[2];
    x[6] = mForces[3];
    x[7] = mForces[4];
    x[8] = mForces[5];
    x[9] = mT0;
    x[10] = mFreq;
  } else
    cout << "Error in Source::get_parameters(), "
         << "function only implemented for moment sources" << endl;
}

//-----------------------------------------------------------------------
void Source::setFrequency(float_sw4 freq) { mFreq = freq; }

//-----------------------------------------------------------------------
void Source::correct_Z_level(EW* a_ew) {
  // this routine
  // 1. calculates the z-coordinate of the topography right above the source and
  // saves it in m_zTopo
  // 2. if m_relativeToTopography == true, it adds m_zTopo to mZ0
  // 3. checks if the source is inside the computational domain. If not, set
  // mIgnore=true

  // tmp
  //   printf("Entering correct_Z_level()\n");

  int i, j, k, g;
  int success = a_ew->computeNearestGridPoint2(i, j, k, g, mX0, mY0, mZ0);
  m_myPoint = success && a_ew->interior_point_in_proc(i, j, g);

  // The following is a safety check to make sure only one processor considers
  // this (i,j) to be interior We could remove this check if we were certain
  // that interior_point_in_proc() never lies
  int iwrite = m_myPoint ? 1 : 0;
  int size;
  MPI_Comm_size(a_ew->m_1d_communicator, &size);
  std::vector<int> whoIsOne(size);
  int counter = 0;
  MPI_Allgather(&iwrite, 1, MPI_INT, &whoIsOne[0], 1, MPI_INT,
                a_ew->m_1d_communicator);
  for (unsigned int p = 0; p < whoIsOne.size(); ++p)
    if (whoIsOne[p] == 1) {
      counter++;
    }
  REQUIRE2(counter == 1,
           "Source error: the nearest grid point should only be interior to "
           "one proc, but counter = "
               << counter << " for source station at (x,y,depth)=" << mX0
               << ", " << mY0 << ", " << mZ0);

  if (!a_ew->topographyExists()) {
    // This is the easy case w/o topography
    m_zTopo = 0.0;
  } else {
    // With topography, compute z-coordinate at topography directly above the
    // source
    float_sw4 zTopoLocal;
    if (!a_ew->m_gridGenerator->interpolate_topography(
            a_ew, mX0, mY0, zTopoLocal, a_ew->mTopoGridExt))
      zTopoLocal = -1e38;
    MPI_Allreduce(&zTopoLocal, &m_zTopo, 1, a_ew->m_mpifloat, MPI_MAX,
                  a_ew->m_1d_communicator);
    if (m_zRelativeToTopography) {
      // If location was specified with topodepth, correct z-level
      mZ0 += m_zTopo;
      m_zRelativeToTopography = false;  // set to false so the correction isn't
                                        // repeated (for whatever reason)
    }
  }

  // Make sure the station is below or on the topography (z is positive
  // downwards)
  if (mZ0 < m_zTopo - 1.0e-9)  // allow for a little roundoff
  {
    mIgnore = true;
    printf(
        "Ignoring Source at X=%g, Y=%g, Z=%g, because it is above the "
        "topography z=%g\n",
        mX0, mY0, mZ0, m_zTopo);
  }
  // tmp
  //   printf("Exiting correct_Z_level()\n");
}

//-----------------------------------------------------------------------
void Source::compute_grid_point(EW* a_ew) {
  // Sets values (m_i0,m_j0,m_k0) and m_grid.
  // Should be called after the topographic correction of mZ0
  int i, j, k, g;
  int success = a_ew->computeNearestGridPoint2(i, j, k, g, mX0, mY0, mZ0);
  m_myPoint = success && a_ew->interior_point_in_proc(i, j, g);
  int inds[4] = {-9999, -9999, -9999, -9999};
  if (m_myPoint) {
    inds[0] = i;
    inds[1] = j;
    inds[2] = k;
    inds[3] = g;
  }
  int indsg[4] = {0, 0, 0, 0};
  MPI_Allreduce(inds, indsg, 4, MPI_INT, MPI_MAX, a_ew->m_1d_communicator);
  m_i0 = indsg[0];
  m_j0 = indsg[1];
  m_k0 = indsg[2];
  m_grid = indsg[3];
  //   std::cout << a_ew->getRank() << " mypoint = " << m_myPoint
  //             << " (i,j,k)= " << m_i0 <<" " << m_j0 << " " << m_k0 << " grid=
  //             "
  //             << m_grid << std::endl;
}

//-----------------------------------------------------------------------
void Source::getsourcewgh(float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6],
                          float_sw4 ddwghda[6]) const {
  // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position
  float_sw4 p5 = ai * ai * ai * ai * ai *
                 (5.0 / 3 - 7.0 / 24 * ai - 17 / 12.0 * ai * ai +
                  1.125 * ai * ai * ai - 0.25 * ai * ai * ai * ai);
  wgh[0] = 1.0 / 24 *
               (2 * ai - ai * ai - 2 * ai * ai * ai - 19 * ai * ai * ai * ai) +
           p5;
  wgh[1] = 1.0 / 6 * (-4 * ai + 4 * ai * ai + ai * ai * ai) +
           4 * ai * ai * ai * ai - 5 * p5;
  wgh[2] = 1 - 1.25 * ai * ai - 97.0 / 12 * ai * ai * ai * ai + 10 * p5;
  wgh[3] =
      1.0 / 6 * (4 * ai + 4 * ai * ai - ai * ai * ai + 49 * ai * ai * ai * ai) -
      10 * p5;
  wgh[4] = 1.0 / 24 * (-2 * ai - ai * ai + 2 * ai * ai * ai) -
           4.125 * ai * ai * ai * ai + 5 * p5;
  wgh[5] = 5.0 / 6 * ai * ai * ai * ai - p5;

  // Derivatives of wgh wrt. ai:
  p5 = 5 * ai * ai * ai * ai *
           (5.0 / 3 - 7.0 / 24 * ai - 17 / 12.0 * ai * ai +
            1.125 * ai * ai * ai - 0.25 * ai * ai * ai * ai) +
       ai * ai * ai * ai * ai *
           (-7.0 / 24 - 17 / 6.0 * ai + 3 * 1.125 * ai * ai - ai * ai * ai);
  dwghda[0] =
      1.0 / 24 * (2 - 2 * ai - 6 * ai * ai - 19 * 4 * ai * ai * ai) + p5;
  dwghda[1] =
      1.0 / 6 * (-4 + 8 * ai + 3 * ai * ai) + 16 * ai * ai * ai - 5 * p5;
  dwghda[2] = -2.5 * ai - 97.0 / 3 * ai * ai * ai + 10 * p5;
  dwghda[3] =
      1.0 / 6 * (4 + 8 * ai - 3 * ai * ai + 49 * 4 * ai * ai * ai) - 10 * p5;
  dwghda[4] = 1.0 / 24 * (-2 - 2 * ai + 6 * ai * ai) -
              4.125 * 4 * ai * ai * ai + 5 * p5;
  dwghda[5] = 20.0 / 6 * ai * ai * ai - p5;

  // Second derivatives of wgh wrt. ai:
  p5 = ai * ai * ai *
       (100.0 / 3 - 8.75 * ai - 59.5 * ai * ai + 63 * ai * ai * ai -
        18 * ai * ai * ai * ai);

  ddwghda[0] = -1.0 / 12 - 0.5 * ai - 9.5 * ai * ai + p5;
  ddwghda[1] = 4.0 / 3 + ai + 48 * ai * ai - 5 * p5;
  ddwghda[2] = -2.5 - 97 * ai * ai + 10 * p5;
  ddwghda[3] = 4.0 / 3 - ai + 98 * ai * ai - 10 * p5;
  ddwghda[4] = -1.0 / 12 + 0.5 * ai - 49.5 * ai * ai + 5 * p5;
  ddwghda[5] = 10 * ai * ai - p5;
}

//-----------------------------------------------------------------------
void Source::getsourcedwgh(float_sw4 ai, float_sw4 wgh[6], float_sw4 dwghda[6],
                           float_sw4 ddwghda[6]) const {
  // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position
  float_sw4 p5 = ai * ai * ai * ai *
                 (-25.0 / 12 - 0.75 * ai + 59.0 / 12 * ai * ai -
                  4 * ai * ai * ai + ai * ai * ai * ai);
  wgh[0] = 1.0 / 12 * (-1 + ai + 3 * ai * ai + 8 * ai * ai * ai) + p5;
  wgh[1] = 2.0 / 3 * (1 - 2 * ai) - 0.5 * ai * ai - 3.5 * ai * ai * ai - 5 * p5;
  wgh[2] = 2.5 * ai + 22.0 / 3 * ai * ai * ai + 10 * p5;
  wgh[3] = 2.0 / 3 * (-1 - 2 * ai) + 0.5 * ai * ai - 23.0 / 3 * ai * ai * ai -
           10 * p5;
  wgh[4] = (1 + ai) / 12 - 0.25 * ai * ai + 4 * ai * ai * ai + 5 * p5;
  wgh[5] = -5.0 / 6 * ai * ai * ai - p5;

  // Derivatives of wgh wrt. ai:
  p5 = 4 * ai * ai * ai *
           (-25.0 / 12 - 0.75 * ai + 59.0 / 12 * ai * ai - 4 * ai * ai * ai +
            ai * ai * ai * ai) +
       ai * ai * ai * ai *
           (-0.75 + 59.0 / 6 * ai - 12 * ai * ai + 4 * ai * ai * ai);
  dwghda[0] = 1.0 / 12 * (1 + 6 * ai + 24 * ai * ai) + p5;
  dwghda[1] = 2.0 / 3 * (-2) - ai - 3 * 3.5 * ai * ai - 5 * p5;
  dwghda[2] = 2.5 + 22.0 * ai * ai + 10 * p5;
  dwghda[3] = 2.0 / 3 * (-2) + ai - 23.0 * ai * ai - 10 * p5;
  dwghda[4] = 1.0 / 12 - 0.5 * ai + 12 * ai * ai + 5 * p5;
  dwghda[5] = -5.0 / 2 * ai * ai - p5;

  // Second derivatives of wgh wrt. ai:
  p5 = ai * ai *
       (-25 - 15 * ai + 147.5 * ai * ai - 168 * ai * ai * ai +
        56 * ai * ai * ai * ai);

  ddwghda[0] = 0.5 + 4 * ai + p5;
  ddwghda[1] = -1 - 21 * ai - 5 * p5;
  ddwghda[2] = 44 * ai + 10 * p5;
  ddwghda[3] = 1 - 46 * ai - 10 * p5;
  ddwghda[4] = -0.5 + 24 * ai + 5 * p5;
  ddwghda[5] = -5 * ai - p5;
}

//-----------------------------------------------------------------------
void Source::getsourcewghlow(float_sw4 ai, float_sw4 wgh[6],
                             float_sw4 dwghda[6], float_sw4 ddwghda[6]) const {
  // Lower component stencil, to use at lower boundaries
  // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position

  wgh[0] = (2 * ai - ai * ai - 2 * ai * ai * ai + ai * ai * ai * ai) / 24;
  wgh[1] = (-4 * ai + 4 * ai * ai + ai * ai * ai - ai * ai * ai * ai) / 6;
  wgh[2] = 1 - 1.25 * ai * ai + 0.25 * ai * ai * ai * ai;
  wgh[3] = (4 * ai + 4 * ai * ai - ai * ai * ai - ai * ai * ai * ai) / 6;
  wgh[4] = (-2 * ai - ai * ai + 2 * ai * ai * ai + ai * ai * ai * ai) / 24;
  wgh[5] = 0;

  // Derivatives of wgh wrt. ai:
  dwghda[0] = (1 - ai - 3 * ai * ai + 2 * ai * ai * ai) / 12;
  dwghda[1] = (-2 + 4 * ai + 1.5 * ai * ai - 2 * ai * ai * ai) / 3;
  dwghda[2] = -2.5 * ai + ai * ai * ai;
  dwghda[3] = (2 + 4 * ai - 1.5 * ai * ai - 2 * ai * ai * ai) / 3;
  dwghda[4] = (-1 - ai + 3 * ai * ai + 2 * ai * ai * ai) / 12;
  dwghda[5] = 0;

  // Second derivatives of wgh wrt. ai:
  ddwghda[0] = -1.0 / 12 - 0.5 * ai + 0.5 * ai * ai;
  ddwghda[1] = 4.0 / 3 + ai - 2 * ai * ai;
  ddwghda[2] = -2.5 + 3 * ai * ai;
  ddwghda[3] = 4.0 / 3 - ai - 2 * ai * ai;
  ddwghda[4] = -1.0 / 12 + 0.5 * ai + 0.5 * ai * ai;
  ddwghda[5] = 0;
}

//-----------------------------------------------------------------------
void Source::getsourcedwghlow(float_sw4 ai, float_sw4 wgh[6],
                              float_sw4 dwghda[6], float_sw4 ddwghda[6]) const {
  // Lower component stencil, to use at lower boundaries, dirac derivative
  // weights. Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position

  // same as derivatives of dirac weights.
  wgh[0] = (-1 + ai + 3 * ai * ai - 2 * ai * ai * ai) / 12;
  wgh[1] = (2 - 4 * ai - 1.5 * ai * ai + 2 * ai * ai * ai) / 3;
  wgh[2] = 2.5 * ai - ai * ai * ai;
  wgh[3] = (-2 - 4 * ai + 1.5 * ai * ai + 2 * ai * ai * ai) / 3;
  wgh[4] = (1 + ai - 3 * ai * ai - 2 * ai * ai * ai) / 12;
  wgh[5] = 0;

  // Derivatives of wgh wrt. ai:
  dwghda[0] = 1.0 / 12 + 0.5 * ai - 0.5 * ai * ai;
  dwghda[1] = -4.0 / 3 - ai + 2 * ai * ai;
  dwghda[2] = 2.5 - 3 * ai * ai;
  dwghda[3] = -4.0 / 3 + ai + 2 * ai * ai;
  dwghda[4] = 1.0 / 12 - 0.5 * ai - 0.5 * ai * ai;
  dwghda[5] = 0;

  // Second derivatives of wgh wrt. ai:
  ddwghda[0] = 0.5 - ai;
  ddwghda[1] = -1 + 4 * ai;
  ddwghda[2] = -6 * ai;
  ddwghda[3] = 1 + 4 * ai;
  ddwghda[4] = -0.5 - ai;
  ddwghda[5] = 0;
}

//-----------------------------------------------------------------------
void Source::getmetwgh(float_sw4 ai, float_sw4 wgh[8], float_sw4 dwgh[8],
                       float_sw4 ddwgh[8], float_sw4 dddwgh[8]) const {
  float_sw4 pol = ai * ai * ai * ai * ai * ai * ai *
                  (-251 + 135 * ai + 25 * ai * ai - 33 * ai * ai * ai +
                   6 * ai * ai * ai * ai) /
                  720;

  wgh[0] = -1.0 / 60 * ai + 1.0 / 180 * ai * ai + 1.0 / 48 * ai * ai * ai +
           23.0 / 144 * ai * ai * ai * ai -
           (17.0 * ai + 223.0) * ai * ai * ai * ai * ai / 720 - pol;
  wgh[1] = 3.0 / 20 * ai - 3.0 / 40 * ai * ai - 1.0 / 6 * ai * ai * ai -
           13.0 / 12 * ai * ai * ai * ai + 97.0 / 45 * ai * ai * ai * ai * ai +
           1.0 / 6 * ai * ai * ai * ai * ai * ai + 7 * pol;
  wgh[2] = -0.75 * ai + 0.75 * ai * ai + (13.0 + 155 * ai) * ai * ai * ai / 48 -
           103.0 / 16 * ai * ai * ai * ai * ai -
           121.0 / 240 * ai * ai * ai * ai * ai * ai - 21 * pol;
  wgh[3] = 1 - 49.0 / 36 * ai * ai - 49.0 / 9 * ai * ai * ai * ai +
           385.0 / 36 * ai * ai * ai * ai * ai +
           61.0 / 72 * ai * ai * ai * ai * ai * ai + 35 * pol;
  wgh[4] = 0.75 * ai + 0.75 * ai * ai - 13.0 / 48 * ai * ai * ai +
           89.0 / 16 * ai * ai * ai * ai -
           1537.0 / 144 * ai * ai * ai * ai * ai -
           41.0 / 48 * ai * ai * ai * ai * ai * ai - 35 * pol;
  wgh[5] = -3.0 / 20 * ai - 3.0 / 40 * ai * ai + 1.0 / 6 * ai * ai * ai -
           41.0 / 12 * ai * ai * ai * ai + 6.4 * ai * ai * ai * ai * ai +
           31.0 / 60 * ai * ai * ai * ai * ai * ai + 21 * pol;
  wgh[6] = 1.0 / 60 * ai + 1.0 / 180 * ai * ai - 1.0 / 48 * ai * ai * ai +
           167.0 / 144 * ai * ai * ai * ai -
           1537.0 / 720 * ai * ai * ai * ai * ai -
           25.0 / 144 * ai * ai * ai * ai * ai * ai - 7 * pol;
  wgh[7] = -1.0 / 6 * ai * ai * ai * ai + 11.0 / 36 * ai * ai * ai * ai * ai +
           1.0 / 40 * ai * ai * ai * ai * ai * ai + pol;

  // Derivative wrt. ai
  pol = ai * ai * ai * ai * ai * ai *
        (-1757.0 / 720 + 1.5 * ai + 0.31250 * ai * ai -
         (1.375 * ai * ai * ai - 0.275 * ai * ai * ai * ai) / 3);
  dwgh[0] = -1.0 / 60 + 1.0 / 90 * ai + ai * ai / 16 +
            23.0 / 36 * ai * ai * ai - 223.0 / 144 * ai * ai * ai * ai -
            17.0 / 120 * ai * ai * ai * ai * ai - pol;
  dwgh[1] = 3.0 / 20 - 3.0 / 20 * ai - 0.5 * ai * ai - 13.0 / 3 * ai * ai * ai +
            97.0 / 9 * ai * ai * ai * ai + ai * ai * ai * ai * ai + 7 * pol;
  dwgh[2] = -0.75 + 1.5 * ai + 13.0 / 16 * ai * ai + 155.0 * ai * ai * ai / 12 -
            103.0 * 5.0 / 16 * ai * ai * ai * ai -
            121.0 / 40 * ai * ai * ai * ai * ai - 21 * pol;
  dwgh[3] = -49.0 / 18 * ai - 4 * 49.0 / 9.0 * ai * ai * ai +
            385.0 * 5.0 / 36 * ai * ai * ai * ai +
            61.0 / 12 * ai * ai * ai * ai * ai + 35 * pol;
  dwgh[4] = 0.75 + 1.5 * ai - 13.0 / 16 * ai * ai + 89.0 / 4 * ai * ai * ai -
            1537.0 * 5 / 144.0 * ai * ai * ai * ai -
            41.0 / 8 * ai * ai * ai * ai * ai - 35 * pol;
  dwgh[5] = -3.0 / 20 - 3.0 / 20 * ai + 0.5 * ai * ai -
            41.0 / 3 * ai * ai * ai + 32 * ai * ai * ai * ai +
            3.1 * ai * ai * ai * ai * ai + 21 * pol;
  dwgh[6] = 1.0 / 60 + 1.0 / 90 * ai - 1.0 / 16 * ai * ai +
            167.0 / 36 * ai * ai * ai - 1537.0 / 144 * ai * ai * ai * ai -
            25.0 / 24 * ai * ai * ai * ai * ai - 7 * pol;
  dwgh[7] = -2.0 / 3 * ai * ai * ai + 55.0 / 36 * ai * ai * ai * ai +
            3.0 / 20 * ai * ai * ai * ai * ai + pol;

  // Second derivative wrt. ai
  pol = ai * ai * ai * ai * ai *
        (-1757.0 / 120 + 10.5 * ai + 2.5 * ai * ai - 4.125 * ai * ai * ai +
         11.0 / 12 * ai * ai * ai * ai);
  ddwgh[0] = 1.0 / 90 + 0.125 * ai + 23.0 / 12 * ai * ai -
             223.0 / 36 * ai * ai * ai - 17.0 / 24 * ai * ai * ai * ai - pol;
  ddwgh[1] = -3.0 / 20 - ai - 13.0 * ai * ai + 4 * 97.0 / 9.0 * ai * ai * ai +
             5 * ai * ai * ai * ai + 7 * pol;
  ddwgh[2] = 1.5 + 13.0 / 8 * ai + 155.0 / 4 * ai * ai -
             103.0 * 5.0 / 4 * ai * ai * ai - 121.0 / 8 * ai * ai * ai * ai -
             21 * pol;
  ddwgh[3] = -49.0 / 18 - 4 * 49.0 / 3.0 * ai * ai +
             385.0 * 5.0 / 9.0 * ai * ai * ai +
             5 * 61.0 / 12 * ai * ai * ai * ai + 35 * pol;
  ddwgh[4] = 1.5 - 13.0 / 8 * ai + 89.0 * 3.0 / 4 * ai * ai -
             1537.0 * 5.0 / 36 * ai * ai * ai - 205.0 / 8 * ai * ai * ai * ai -
             35 * pol;
  ddwgh[5] = -3.0 / 20 + ai - 41.0 * ai * ai + 128 * ai * ai * ai +
             15.5 * ai * ai * ai * ai + 21 * pol;
  ddwgh[6] = 1.0 / 90 - 0.125 * ai + 167.0 / 12 * ai * ai -
             1537.0 / 36 * ai * ai * ai - 125.0 / 24 * ai * ai * ai * ai -
             7 * pol;
  ddwgh[7] =
      -2 * ai * ai + 220.0 / 36 * ai * ai * ai + 0.75 * ai * ai * ai * ai + pol;

  // Third derivative wrt. ai
  pol = ai * ai * ai * ai *
        (-1757.0 / 24 + 63 * ai + 17.5 * ai * ai - 33 * ai * ai * ai +
         8.25 * ai * ai * ai * ai);
  dddwgh[0] = 0.125 + 23.0 / 6 * ai - 223.0 / 12 * ai * ai -
              17.0 / 6 * ai * ai * ai - pol;
  dddwgh[1] =
      -1 - 26.0 * ai + 4 * 97.0 / 3 * ai * ai + 20 * ai * ai * ai + 7 * pol;
  dddwgh[2] =
      1.625 + 77.5 * ai - 386.25 * ai * ai - 60.5 * ai * ai * ai - 21 * pol;
  dddwgh[3] = -392.0 / 3 * ai + 1925.0 / 3 * ai * ai +
              305.0 / 3 * ai * ai * ai + 35 * pol;
  dddwgh[4] = -1.625 + 133.5 * ai - 7685.0 / 12 * ai * ai -
              102.5 * ai * ai * ai - 35 * pol;
  dddwgh[5] = 1 - 82.0 * ai + 384.0 * ai * ai + 62.0 * ai * ai * ai + 21 * pol;
  dddwgh[6] = -0.125 + 167.0 / 6 * ai - 1537.0 / 12 * ai * ai -
              125.0 / 6 * ai * ai * ai - 7 * pol;
  dddwgh[7] = -4 * ai + 220.0 / 12 * ai * ai + 3 * ai * ai * ai + pol;
}

//-----------------------------------------------------------------------
void Source::getmetdwgh(float_sw4 ai, float_sw4 wgh[8]) const {
  float_sw4 pol = ai * ai * ai * ai * ai * ai *
                  (-827 + 420 * ai + 165 * ai * ai - 180 * ai * ai * ai +
                   36 * ai * ai * ai * ai) /
                  720;

  wgh[0] = -1.0 / 60 + 1.0 / 90 * ai + 1.0 / 16 * ai * ai +
           5.0 / 36 * ai * ai * ai - 55.0 / 144 * ai * ai * ai * ai -
           7.0 / 20 * ai * ai * ai * ai * ai - pol;
  wgh[1] = 3.0 / 20 * (1 - ai) - 0.5 * ai * ai - 5.0 / 6 * ai * ai * ai +
           47.0 / 18 * ai * ai * ai * ai + 59.0 / 24 * ai * ai * ai * ai * ai +
           7 * pol;
  wgh[2] = -0.75 + 1.5 * ai + 13.0 / 16 * ai * ai + 29.0 / 12 * ai * ai * ai -
           123.0 / 16 * ai * ai * ai * ai - 7.4 * ai * ai * ai * ai * ai -
           21 * pol;
  wgh[3] = (-49.0 * ai - 77.0 * ai * ai * ai) / 18 +
           455.0 / 36 * ai * ai * ai * ai + 99.0 / 8 * ai * ai * ai * ai * ai +
           35 * pol;
  wgh[4] = 0.75 + 1.5 * ai - 13.0 / 16 * ai * ai + 4.75 * ai * ai * ai -
           1805.0 / 144 * ai * ai * ai * ai -
           149.0 / 12 * ai * ai * ai * ai * ai - 35 * pol;
  wgh[5] = -3.0 / 20 * (1 + ai) + 0.5 * ai * ai - 19.0 / 6 * ai * ai * ai +
           7.5 * ai * ai * ai * ai + 299.0 / 40 * ai * ai * ai * ai * ai +
           21 * pol;
  wgh[6] = 1.0 / 60 + 1.0 / 90 * ai - 1.0 / 16 * ai * ai +
           41.0 / 36 * ai * ai * ai - 361.0 / 144 * ai * ai * ai * ai -
           2.5 * ai * ai * ai * ai * ai - 7 * pol;
  wgh[7] = -1.0 / 6 * ai * ai * ai + 13.0 / 36 * ai * ai * ai * ai +
           43.0 / 120 * ai * ai * ai * ai * ai + pol;
}
//-----------------------------------------------------------------------
void Source::getmetwgh7(float_sw4 ai, float_sw4 wgh[7]) const {
  wgh[0] = ai * ai / 180.0 - ai * ai * ai * ai / 144.0 +
           ai * ai * ai * ai * ai * ai / 720.0 - ai / 60.0 +
           ai * ai * ai / 48.0 - ai * ai * ai * ai * ai / 240.0;
  wgh[1] = -3.0 / 40.0 * ai * ai + ai * ai * ai * ai / 12.0 -
           ai * ai * ai * ai * ai * ai / 120.0 + 3.0 / 20.0 * ai -
           ai * ai * ai / 6.0 + ai * ai * ai * ai * ai / 60.0;
  wgh[2] = 3.0 / 4.0 * ai * ai - 13.0 / 48.0 * ai * ai * ai * ai +
           ai * ai * ai * ai * ai * ai / 48.0 - 3.0 / 4.0 * ai +
           13.0 / 48.0 * ai * ai * ai - ai * ai * ai * ai * ai / 48.0;
  wgh[3] = 1.0 - 49.0 / 36.0 * ai * ai + 7.0 / 18.0 * ai * ai * ai * ai -
           ai * ai * ai * ai * ai * ai / 36.0;
  wgh[4] = 3.0 / 4.0 * ai * ai - 13.0 / 48.0 * ai * ai * ai * ai +
           3.0 / 4.0 * ai - 13.0 / 48.0 * ai * ai * ai +
           ai * ai * ai * ai * ai / 48.0 + ai * ai * ai * ai * ai * ai / 48.0;
  wgh[5] = -3.0 / 40.0 * ai * ai + ai * ai * ai * ai / 12.0 - 3.0 / 20.0 * ai +
           ai * ai * ai / 6.0 - ai * ai * ai * ai * ai / 60.0 -
           ai * ai * ai * ai * ai * ai / 120.0;
  wgh[6] = ai * ai / 180.0 - ai * ai * ai * ai / 144.0 +
           ai * ai * ai * ai * ai * ai / 720.0 + ai / 60.0 +
           ai * ai * ai * ai * ai / 240.0 - ai * ai * ai / 48.0;
}

//-----------------------------------------------------------------------
void Source::getmetdwgh7(float_sw4 ai, float_sw4 wgh[7]) const {
  wgh[0] = -1.0 / 60.0 + ai * ai / 16.0 - ai * ai * ai * ai / 48.0 + ai / 90.0 -
           ai * ai * ai / 36.0 + ai * ai * ai * ai * ai / 120.0;
  wgh[1] = 3.0 / 20.0 - ai * ai / 2.0 + ai * ai * ai * ai / 12.0 -
           3.0 / 20.0 * ai + ai * ai * ai / 3.0 - ai * ai * ai * ai * ai / 20.0;
  wgh[2] = -3.0 / 4.0 + 13.0 / 16.0 * ai * ai - 5.0 / 48.0 * ai * ai * ai * ai +
           3.0 / 2.0 * ai - 13.0 / 12.0 * ai * ai * ai +
           ai * ai * ai * ai * ai / 8.0;
  wgh[3] = -49.0 / 18.0 * ai + 14.0 / 9.0 * ai * ai * ai -
           ai * ai * ai * ai * ai / 6.0;
  wgh[4] = 3.0 / 4.0 - 13.0 / 16.0 * ai * ai + 3.0 / 2.0 * ai -
           13.0 / 12.0 * ai * ai * ai + 5.0 / 48.0 * ai * ai * ai * ai +
           ai * ai * ai * ai * ai / 8.0;
  wgh[5] = -3.0 / 20.0 + ai * ai / 2.0 - ai * ai * ai * ai / 12.0 -
           3.0 / 20.0 * ai + ai * ai * ai / 3.0 - ai * ai * ai * ai * ai / 20.0;
  wgh[6] = 1.0 / 60.0 - ai * ai / 16.0 + ai * ai * ai * ai / 48.0 + ai / 90.0 -
           ai * ai * ai / 36.0 + ai * ai * ai * ai * ai / 120.0;
}

//-----------------------------------------------------------------------
// Sources at grid refinement boundary
void Source::getsourcewghNM2sm6(float_sw4 alph, float_sw4 wghk[6]) const {
  wghk[0] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
      (alph * alph * alph * alph / 24.0 + alph / 12.0 - alph * alph / 24.0 -
       alph * alph * alph / 12.0);
  wghk[1] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (-alph * alph * alph * alph / 6.0 - 2.0 / 3.0 * alph +
           2.0 / 3.0 * alph * alph + alph * alph * alph / 6.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (-pow(alph - 1.0, 2.0) / 30.0 + alph / 10.0 - 1.0 / 10.0 +
           pow(alph - 1.0, 4.0) / 30.0 - pow(alph - 1.0, 3.0) / 10.0);
  wghk[2] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (1.0 + alph * alph * alph * alph / 4.0 - 5.0 / 4.0 * alph * alph) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (5.0 / 8.0 * pow(alph - 1.0, 2.0) - 3.0 / 4.0 * alph + 3.0 / 4.0 -
           pow(alph - 1.0, 4.0) / 8.0 + pow(alph - 1.0, 3.0) / 4.0);
  wghk[3] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (-alph * alph * alph / 6.0 - alph * alph * alph * alph / 6.0 +
           2.0 / 3.0 * alph + 2.0 / 3.0 * alph * alph) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (5.0 / 6.0 - 7.0 / 6.0 * pow(alph - 1.0, 2.0) + alph / 6.0 +
           pow(alph - 1.0, 4.0) / 6.0 - pow(alph - 1.0, 3.0) / 6.0);
  wghk[4] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (alph * alph * alph * alph / 24.0 - alph / 12.0 - alph * alph / 24.0 +
           alph * alph * alph / 12.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (7.0 / 12.0 * pow(alph - 1.0, 2.0) + alph / 2.0 - 1.0 / 2.0 -
           pow(alph - 1.0, 4.0) / 12.0);
  wghk[5] = (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
             6.0 * alph * alph * alph * alph * alph) *
            (-pow(alph - 1.0, 2.0) / 120.0 - alph / 60.0 + 1.0 / 60.0 +
             pow(alph - 1.0, 4.0) / 120.0 + pow(alph - 1.0, 3.0) / 60.0);
}

//-----------------------------------------------------------------------
void Source::getsourcedwghNM2sm6(float_sw4 alph, float_sw4 dwghk[6]) const {
  dwghk[0] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
      (-1.0 / 12.0 - alph * alph * alph / 6.0 + alph * alph / 4.0 +
       alph / 12.0);
  dwghk[1] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (2.0 / 3.0 + 2.0 / 3.0 * alph * alph * alph - alph * alph / 2.0 -
           4.0 / 3.0 * alph) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (-1.0 / 6.0 - 2.0 / 15.0 * pow(alph - 1.0, 3.0) +
           3.0 / 10.0 * pow(alph - 1.0, 2.0) + alph / 15.0);
  dwghk[2] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (-alph * alph * alph + 5.0 / 2.0 * alph) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (2.0 + pow(alph - 1.0, 3.0) / 2.0 - 3.0 / 4.0 * pow(alph - 1.0, 2.0) -
           5.0 / 4.0 * alph);
  dwghk[3] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (-2.0 / 3.0 + 2.0 / 3.0 * alph * alph * alph + alph * alph / 2.0 -
           4.0 / 3.0 * alph) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (-5.0 / 2.0 - 2.0 / 3.0 * pow(alph - 1.0, 3.0) + 7.0 / 3.0 * alph +
           pow(alph - 1.0, 2.0) / 2.0);
  dwghk[4] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (1.0 / 12.0 - alph * alph * alph / 6.0 - alph * alph / 4.0 +
           alph / 12.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (2.0 / 3.0 + pow(alph - 1.0, 3.0) / 3.0 - 7.0 / 6.0 * alph);
  dwghk[5] = (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
              6.0 * alph * alph * alph * alph * alph) *
             (-pow(alph - 1.0, 3.0) / 30.0 - pow(alph - 1.0, 2.0) / 20.0 +
              alph / 60.0);
}

//-----------------------------------------------------------------------
void Source::getsourcewghNM1sm6(float_sw4 alph, float_sw4 wghk[6]) const {
  wghk[0] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
      (-alph * alph / 30.0 + alph / 10.0 + alph * alph * alph * alph / 30.0 -
       alph * alph * alph / 10.0);
  wghk[1] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (5.0 / 8.0 * alph * alph - 3.0 / 4.0 * alph -
           alph * alph * alph * alph / 8.0 + alph * alph * alph / 4.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (-5.0 / 48.0 * pow(alph - 1.0, 3.0) + pow(alph - 1.0, 4.0) / 48.0 +
           alph / 6.0 - 1.0 / 6.0 + pow(alph - 1.0, 2.0) / 24.0);
  wghk[2] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (1.0 - 7.0 / 6.0 * alph * alph + alph / 6.0 +
           alph * alph * alph * alph / 6.0 - alph * alph * alph / 6.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (4.0 / 15.0 * pow(alph - 1.0, 3.0) - pow(alph - 1.0, 4.0) / 15.0 -
           16.0 / 15.0 * alph + 16.0 / 15.0 +
           4.0 / 15.0 * pow(alph - 1.0, 2.0));
  wghk[3] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (7.0 / 12.0 * alph * alph + alph / 2.0 -
           alph * alph * alph * alph / 12.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (1.0 / 4.0 + 3.0 / 4.0 * alph - 3.0 / 16.0 * pow(alph - 1.0, 3.0) -
           pow(alph - 1.0, 2.0) / 2.0 + pow(alph - 1.0, 4.0) / 16.0);
  wghk[4] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (-alph * alph / 120.0 - alph / 60.0 +
           alph * alph * alph * alph / 120.0 + alph * alph * alph / 60.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (pow(alph - 1.0, 3.0) / 48.0 - pow(alph - 1.0, 4.0) / 48.0 +
           alph / 6.0 - 1.0 / 6.0 + 5.0 / 24.0 * pow(alph - 1.0, 2.0));
  wghk[5] = (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
             6.0 * alph * alph * alph * alph * alph) *
            (pow(alph - 1.0, 3.0) / 240.0 + pow(alph - 1.0, 4.0) / 240.0 -
             alph / 60.0 + 1.0 / 60.0 - pow(alph - 1.0, 2.0) / 60.0);
}

//-----------------------------------------------------------------------
void Source::getsourcedwghNM1sm6(float_sw4 alph, float_sw4 dwghk[6]) const {
  dwghk[0] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
      (-1.0 / 10.0 - 2.0 / 15.0 * alph * alph * alph +
       3.0 / 10.0 * alph * alph + alph / 15.0);
  dwghk[1] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (3.0 / 4.0 + alph * alph * alph / 2.0 - 3.0 / 4.0 * alph * alph -
           5.0 / 4.0 * alph) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (-alph / 12.0 - 1.0 / 12.0 + 5.0 / 16.0 * pow(alph - 1.0, 2.0) -
           pow(alph - 1.0, 3.0) / 12.0);
  dwghk[2] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (-1.0 / 6.0 - 2.0 / 3.0 * alph * alph * alph + 7.0 / 3.0 * alph +
           alph * alph / 2.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (-8.0 / 15.0 * alph + 8.0 / 5.0 - 4.0 / 5.0 * pow(alph - 1.0, 2.0) +
           4.0 / 15.0 * pow(alph - 1.0, 3.0));
  dwghk[3] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (-1.0 / 2.0 + alph * alph * alph / 3.0 - 7.0 / 6.0 * alph) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (-7.0 / 4.0 + 9.0 / 16.0 * pow(alph - 1.0, 2.0) + alph -
           pow(alph - 1.0, 3.0) / 4.0);
  dwghk[4] =
      (1.0 - 10.0 * alph * alph * alph + 15.0 * alph * alph * alph * alph -
       6.0 * alph * alph * alph * alph * alph) *
          (1.0 / 60.0 - alph * alph * alph / 30.0 - alph * alph / 20.0 +
           alph / 60.0) +
      (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
       6.0 * alph * alph * alph * alph * alph) *
          (1.0 / 4.0 - pow(alph - 1.0, 2.0) / 16.0 - 5.0 / 12.0 * alph +
           pow(alph - 1.0, 3.0) / 12.0);
  dwghk[5] = (10.0 * alph * alph * alph - 15.0 * alph * alph * alph * alph +
              6.0 * alph * alph * alph * alph * alph) *
             (-1.0 / 60.0 - pow(alph - 1.0, 2.0) / 80.0 + alph / 30.0 -
              pow(alph - 1.0, 3.0) / 60.0);
}

//-----------------------------------------------------------------------
void Source::getsourcewghNsm6(float_sw4 alph, float_sw4 wghk[6]) const {
  wghk[0] =
      (1.0 - 5.0 / 4.0 * alph * alph * alph +
       15.0 / 16.0 * alph * alph * alph * alph -
       3.0 / 16.0 * alph * alph * alph * alph * alph) *
      (-5.0 / 48.0 * alph * alph * alph + alph * alph * alph * alph / 48.0 +
       alph / 6.0 + alph * alph / 24.0);
  wghk[1] =
      (1.0 - 5.0 / 4.0 * alph * alph * alph +
       15.0 / 16.0 * alph * alph * alph * alph -
       3.0 / 16.0 * alph * alph * alph * alph * alph) *
          (4.0 / 15.0 * alph * alph * alph - alph * alph * alph * alph / 15.0 -
           16.0 / 15.0 * alph + 4.0 / 15.0 * alph * alph) +
      (5.0 / 4.0 * alph * alph * alph -
       15.0 / 16.0 * alph * alph * alph * alph +
       3.0 / 16.0 * alph * alph * alph * alph * alph) *
          (-4.0 / 105.0 * pow(alph - 2.0, 2.0) + pow(alph - 2.0, 4.0) / 105.0 +
           16.0 / 105.0 * alph - 32.0 / 105.0 -
           4.0 / 105.0 * pow(alph - 2.0, 3.0));
  wghk[2] =
      (1.0 - 5.0 / 4.0 * alph * alph * alph +
       15.0 / 16.0 * alph * alph * alph * alph -
       3.0 / 16.0 * alph * alph * alph * alph * alph) *
          (1.0 + 3.0 / 4.0 * alph - 3.0 / 16.0 * alph * alph * alph -
           alph * alph / 2.0 + alph * alph * alph * alph / 16.0) +
      (5.0 / 4.0 * alph * alph * alph -
       15.0 / 16.0 * alph * alph * alph * alph +
       3.0 / 16.0 * alph * alph * alph * alph * alph) *
          (5.0 / 24.0 * pow(alph - 2.0, 2.0) - pow(alph - 2.0, 4.0) / 48.0 -
           alph / 2.0 + 1.0 + pow(alph - 2.0, 3.0) / 16.0);
  wghk[3] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
             15.0 / 16.0 * alph * alph * alph * alph -
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (alph * alph * alph / 48.0 - alph * alph * alph * alph / 48.0 +
                 alph / 6.0 + 5.0 / 24.0 * alph * alph) +
            (5.0 / 4.0 * alph * alph * alph -
             15.0 / 16.0 * alph * alph * alph * alph +
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (5.0 / 6.0 + pow(alph - 2.0, 4.0) / 48.0 + alph / 12.0 -
                 pow(alph - 2.0, 3.0) / 48.0 - pow(alph - 2.0, 2.0) / 3.0);
  wghk[4] =
      (1.0 - 5.0 / 4.0 * alph * alph * alph +
       15.0 / 16.0 * alph * alph * alph * alph -
       3.0 / 16.0 * alph * alph * alph * alph * alph) *
          (alph * alph * alph / 240.0 + alph * alph * alph * alph / 240.0 -
           alph / 60.0 - alph * alph / 60.0) +
      (5.0 / 4.0 * alph * alph * alph -
       15.0 / 16.0 * alph * alph * alph * alph +
       3.0 / 16.0 * alph * alph * alph * alph * alph) *
          (-pow(alph - 2.0, 3.0) / 80.0 + 7.0 / 40.0 * pow(alph - 2.0, 2.0) -
           pow(alph - 2.0, 4.0) / 80.0 + 3.0 / 10.0 * alph - 3.0 / 5.0);
  wghk[5] = (5.0 / 4.0 * alph * alph * alph -
             15.0 / 16.0 * alph * alph * alph * alph +
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
            (-pow(alph - 2.0, 2.0) / 84.0 + pow(alph - 2.0, 4.0) / 336.0 -
             alph / 28.0 + 1.0 / 14.0 + pow(alph - 2.0, 3.0) / 112.0);
}

//-----------------------------------------------------------------------
void Source::getsourcedwghNsm6(float_sw4 alph, float_sw4 dwghk[6]) const {
  //      dwghk[0] =
  //      (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
  // 16.0*alph*alph*alph*alph*alph)*(-alph/12.0-1.0/6.0+5.0/16.0*alph*alph-alph*alph
  //*alph/12.0);
  //      dwghk[1] =
  //      (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
  // 16.0*alph*alph*alph*alph*alph)*(-8.0/15.0*alph+16.0/15.0-4.0/5.0*alph*alph+4.0/
  // 15.0*alph*alph*alph)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/
  // 16.0*alph*alph*alph*alph*alph)*(-alph/12.0+1.0/12.0+5.0/16.0*pow(alph-3.0,2.0)-
  // pow(alph-3.0,3.0)/12.0);
  //      dwghk[2] =
  //      (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
  // 16.0*alph*alph*alph*alph*alph)*(-3.0/4.0+9.0/16.0*alph*alph+alph-alph*alph*alph
  /// 4.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph*alph*
  // alph*alph*alph)*(-8.0/15.0*alph+8.0/3.0-4.0/5.0*pow(alph-3.0,2.0)+4.0/15.0*pow(
  // alph-3.0,3.0));
  //      dwghk[3] =
  //      (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
  // 16.0*alph*alph*alph*alph*alph)*(-1.0/6.0-alph*alph/16.0-5.0/12.0*alph+alph*alph
  //*alph/12.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph
  //*alph*alph*alph*alph)*(-15.0/4.0+9.0/16.0*pow(alph-3.0,2.0)+alph-pow(alph-3.0,
  // 3.0)/4.0);
  //      dwghk[4] =
  //      (1.0-5.0/4.0*alph*alph*alph+15.0/16.0*alph*alph*alph*alph-3.0/
  // 16.0*alph*alph*alph*alph*alph)*(1.0/60.0-alph*alph/80.0+alph/30.0-alph*alph*
  // alph/60.0)+(5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*alph*
  // alph*alph*alph*alph)*(13.0/12.0-pow(alph-3.0,2.0)/16.0-5.0/12.0*alph+pow(alph
  //-3.0,3.0)/12.0);
  //      dwghk[5] =
  //      (5.0/4.0*alph*alph*alph-15.0/16.0*alph*alph*alph*alph+3.0/16.0*
  // alph*alph*alph*alph*alph)*(-1.0/12.0-pow(alph-3.0,2.0)/80.0+alph/30.0-pow(alph
  //-3.0,3.0)/60.0);

  dwghk[0] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
             (-alph / 12.0 - 1.0 / 6.0 + 5.0 / 16.0 * alph * alph -
              alph * alph * alph / 12.0);
  dwghk[1] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (-8.0 / 15.0 * alph + 16.0 / 15.0 - 4.0 / 5.0 * alph * alph +
                  4.0 / 15.0 * alph * alph * alph) +
             (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (-32.0 / 105.0 - 4.0 / 105.0 * pow(alph - 2.0, 3.0) +
                  4.0 / 35.0 * pow(alph - 2.0, 2.0) + 8.0 / 105.0 * alph);
  dwghk[2] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (-3.0 / 4.0 + 9.0 / 16.0 * alph * alph + alph -
                  alph * alph * alph / 4.0) +
             (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (4.0 / 3.0 + pow(alph - 2.0, 3.0) / 12.0 -
                  3.0 / 16.0 * pow(alph - 2.0, 2.0) - 5.0 / 12.0 * alph);
  dwghk[3] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (-1.0 / 6.0 - alph * alph / 16.0 - 5.0 / 12.0 * alph +
                  alph * alph * alph / 12.0) +
             (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (-17.0 / 12.0 - pow(alph - 2.0, 3.0) / 12.0 +
                  pow(alph - 2.0, 2.0) / 16.0 + 2.0 / 3.0 * alph);
  dwghk[4] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (1.0 / 60.0 - alph * alph / 80.0 + alph / 30.0 -
                  alph * alph * alph / 60.0) +
             (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (2.0 / 5.0 + pow(alph - 2.0, 3.0) / 20.0 +
                  3.0 / 80.0 * pow(alph - 2.0, 2.0) - 7.0 / 20.0 * alph);
  dwghk[5] = (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
             (-1.0 / 84.0 - pow(alph - 2.0, 3.0) / 84.0 -
              3.0 / 112.0 * pow(alph - 2.0, 2.0) + alph / 42.0);
}

//-----------------------------------------------------------------------
void Source::getsourcewghP1sm6(float_sw4 alph, float_sw4 wghk[6]) const {
  wghk[0] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
             15.0 / 16.0 * alph * alph * alph * alph -
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
            (-4.0 / 105.0 * alph * alph + alph * alph * alph * alph / 105.0 +
             16.0 / 105.0 * alph - 4.0 / 105.0 * alph * alph * alph);
  wghk[1] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
             15.0 / 16.0 * alph * alph * alph * alph -
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (5.0 / 24.0 * alph * alph - alph * alph * alph * alph / 48.0 -
                 alph / 2.0 + alph * alph * alph / 16.0) +
            (5.0 / 4.0 * alph * alph * alph -
             15.0 / 16.0 * alph * alph * alph * alph +
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (-pow(alph - 2.0, 2.0) / 96.0 + pow(alph - 2.0, 4.0) / 384.0 +
                 alph / 24.0 - 1.0 / 12.0 - pow(alph - 2.0, 3.0) / 96.0);
  wghk[2] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
             15.0 / 16.0 * alph * alph * alph * alph -
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (1.0 + alph * alph * alph * alph / 48.0 + alph / 12.0 -
                 alph * alph * alph / 48.0 - alph * alph / 3.0) +
            (5.0 / 4.0 * alph * alph * alph -
             15.0 / 16.0 * alph * alph * alph * alph +
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (pow(alph - 2.0, 2.0) / 6.0 - pow(alph - 2.0, 4.0) / 96.0 -
                 alph / 3.0 + 2.0 / 3.0 + pow(alph - 2.0, 3.0) / 48.0);
  wghk[3] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
             15.0 / 16.0 * alph * alph * alph * alph -
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (-alph * alph * alph / 80.0 + 7.0 / 40.0 * alph * alph -
                 alph * alph * alph * alph / 80.0 + 3.0 / 10.0 * alph) +
            (5.0 / 4.0 * alph * alph * alph -
             15.0 / 16.0 * alph * alph * alph * alph +
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (1.0 - 5.0 / 16.0 * pow(alph - 2.0, 2.0) +
                 pow(alph - 2.0, 4.0) / 64.0);
  wghk[4] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
             15.0 / 16.0 * alph * alph * alph * alph -
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (-alph * alph / 84.0 + alph * alph * alph * alph / 336.0 -
                 alph / 28.0 + alph * alph * alph / 112.0) +
            (5.0 / 4.0 * alph * alph * alph -
             15.0 / 16.0 * alph * alph * alph * alph +
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
                (-pow(alph - 2.0, 3.0) / 48.0 + pow(alph - 2.0, 2.0) / 6.0 -
                 pow(alph - 2.0, 4.0) / 96.0 + alph / 3.0 - 2.0 / 3.0);
  wghk[5] = (5.0 / 4.0 * alph * alph * alph -
             15.0 / 16.0 * alph * alph * alph * alph +
             3.0 / 16.0 * alph * alph * alph * alph * alph) *
            (-pow(alph - 2.0, 2.0) / 96.0 + pow(alph - 2.0, 4.0) / 384.0 -
             alph / 24.0 + 1.0 / 12.0 + pow(alph - 2.0, 3.0) / 96.0);
}

//-----------------------------------------------------------------------
void Source::getsourcedwghP1sm6(float_sw4 alph, float_sw4 dwghk[6]) const {
  dwghk[0] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
             (-16.0 / 105.0 - 4.0 / 105.0 * alph * alph * alph +
              4.0 / 35.0 * alph * alph + 8.0 / 105.0 * alph);
  dwghk[1] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (1.0 / 2.0 + alph * alph * alph / 12.0 -
                  3.0 / 16.0 * alph * alph - 5.0 / 12.0 * alph) +
             (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (-1.0 / 12.0 - pow(alph - 2.0, 3.0) / 96.0 +
                  pow(alph - 2.0, 2.0) / 32.0 + alph / 48.0);
  dwghk[2] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (-1.0 / 12.0 - alph * alph * alph / 12.0 + alph * alph / 16.0 +
                  2.0 / 3.0 * alph) +
             (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (1.0 + pow(alph - 2.0, 3.0) / 24.0 -
                  pow(alph - 2.0, 2.0) / 16.0 - alph / 3.0);
  dwghk[3] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (-3.0 / 10.0 + alph * alph * alph / 20.0 +
                  3.0 / 80.0 * alph * alph - 7.0 / 20.0 * alph) +
             (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (5.0 / 8.0 * alph - 5.0 / 4.0 - pow(alph - 2.0, 3.0) / 16.0);
  dwghk[4] = (1.0 - 5.0 / 4.0 * alph * alph * alph +
              15.0 / 16.0 * alph * alph * alph * alph -
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (1.0 / 28.0 - alph * alph * alph / 84.0 -
                  3.0 / 112.0 * alph * alph + alph / 42.0) +
             (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
                 (1.0 / 3.0 + pow(alph - 2.0, 3.0) / 24.0 +
                  pow(alph - 2.0, 2.0) / 16.0 - alph / 3.0);
  dwghk[5] = (5.0 / 4.0 * alph * alph * alph -
              15.0 / 16.0 * alph * alph * alph * alph +
              3.0 / 16.0 * alph * alph * alph * alph * alph) *
             (-pow(alph - 2.0, 3.0) / 96.0 - pow(alph - 2.0, 2.0) / 32.0 +
              alph / 48.0);
}

//------ filter and fix up any discrete time functions ----------------
void Source::prepareTimeFunc(bool doFilter, float_sw4 sw4TimeStep,
                             int sw4TimeSamples, Filter* sw4_filter) {
  if (mTimeDependence == iDiscrete || mTimeDependence == iDiscrete6moments ||
      mTimeDependence == iDiscrete3forces) {
    // new approach for Discrete time functions
    if (!doFilter) {
      spline_interpolation();
      m_timeFuncIsReady = true;
    } else {
      // 0. build filter for the discrete time series
      // 1. compute pre and post durations of filter tail
      // 2. extend time grid by appropriate number of time samples, pad with
      // zeros
      // 3. Perform filtering with the dt of the time series
      // 4. Interpolate by spline
      Filter my_filter(sw4_filter->get_type(), sw4_filter->get_order(),
                       sw4_filter->get_passes(), sw4_filter->get_corner_freq1(),
                       sw4_filter->get_corner_freq2());
      // time step in the time series
      float_sw4 dt = 1 / mFreq;
      // setup the filter for time step dt
      my_filter.computeSOS(dt);
      float_sw4 preCursorTime = my_filter.estimatePrecursor();
      int nPadding = static_cast<int>(ceil(preCursorTime / dt));
      // current number of points
      int npts = mIpar[0];
      float_sw4 tstart = mPar[0];
      int ext_npts = npts + 2 * nPadding;
      float_sw4* ext_par = SW4_NEW(Space::Managed, float_sw4[ext_npts + 1]);
      float_sw4 ext_tstart = tstart - nPadding * dt;
      // setup ext_par
      ext_par[0] = ext_tstart;
      int i;
      // initial zeros
      for (i = 0; i < nPadding; i++) ext_par[i + 1] = 0.0;
// copy values from origianla time series
#pragma omp parallel for
      for (int i = 0; i < npts; i++) ext_par[i + nPadding + 1] = mPar[i + 1];
      // trailing zeros
      for (i = nPadding + npts; i < ext_npts; i++) ext_par[i + 1] = 0.0;

      // Filter the extended time series (ext_par[0] = ext_tstart and should not
      // be filtered)
      my_filter.evaluate(ext_npts, &ext_par[1], &ext_par[1]);

      // Give the source time function a smooth start if this is a 2-pass
      // (forward + backward) bandpass filter (copied from below)
      if (my_filter.get_passes() == 2 && my_filter.get_type() == bandPass) {
        float_sw4 wghv, xi;
        int p0 = 3,
            p = 20;  // First non-zero time level, and number of points in ramp;
        if (p0 + p <=
            ext_npts)  // only do this if the time series is long enough
        {
          for (int i = 1; i <= p0 - 1; i++) {
            ext_par[i] = 0;
          }
          for (int i = p0; i <= p0 + p; i++) {
            wghv = 0;
            xi = (i - p0) / ((float_sw4)p);
            // polynomial P(xi), P(0) = 0, P(1)=1
            wghv = xi * xi * xi * xi *
                   (35 - 84 * xi + 70 * xi * xi - 20 * xi * xi * xi);
            ext_par[i] *= wghv;
          }
        }
      }

      // update parameters for the discrete function
      mNipar = 1;
      //      mIpar = new int[mNipar];
      mIpar[0] = ext_npts;
      //      mFreq = 1./dt;
      ::operator delete[](
          mPar, Space::Managed);  // return memory for the previous time series
      mNpar = ext_npts + 1;
      mPar = ext_par;
      mPar[0] = ext_tstart;  // regular (like Gaussian) time functions are
                             // defined from t=tstart=0
      mT0 = ext_tstart;
      // for( int i=0 ; i < nsteps; i++ )
      //    mPar[i+1] = discfunc[i];
      // delete[] discfunc;

      // Build the spline representation
      spline_interpolation();
      m_is_filtered = true;
      // only do the filtering once
      m_timeFuncIsReady = true;
    }  // end doFilter

  }     // end if iDiscrete or iDiscrete6Moments
  else  // all other time functions
  {
    // Modify the time functions if prefiltering is enabled
    if (doFilter) {
      // the sw4 filter is assumed to already been initialized for the sw4 time
      // step

      // 1. Make sure the smallest time offset is at least t0_min + (timeFcn
      // dependent offset for centered fcn's)
      float_sw4 t0_inc = 0;
      float_sw4 t0_min;
      t0_min = sw4_filter->estimatePrecursor();

      // here we only have one time function
      // // old estimate for 2-pole low-pass Butterworth
      // //	t0_min = 4./sw4_filter->get_corner_freq2();

      t0_inc = compute_t0_increase(t0_min);

      // If t0_inc is positive, the t0 field in the source command should be
      // incremented by at least this amount. Otherwise, there might be
      // significant artifacts from a sudden start of some source.
      if (t0_inc > 0.) {
        // Don't mess with t0.
        // Instead, warn the user of potential transients due to unsmooth start
        printf(
            "\n*** WARNING: the 2 pass prefilter has an estimated precursor of "
            "length %e s\n"
            "*** To avoid artifacts due to sudden startup, increase t0 in the "
            "source named '%s' by at least %e\n\n",
            t0_min, getName().c_str(), t0_inc);
      }

      // Do the actual filtering
      // this function no longer handles discrete time functions
      filter_timefunc(sw4_filter, 0.0, sw4TimeStep, sw4TimeSamples);
    }  // end if doFilter

    // set the flag to indicate that the filtering is complete
    m_timeFuncIsReady = true;
  }  // end if timeDep != iDiscrete
}

//-----------------------------------------------------------------------
void Source::set_grid_point_sources4(EW* a_EW,
                                     vector<GridPointSource*>& point_sources) {
  // This routine is called from all processors, for each input source.
  int i, j, k, g;
  int success = a_EW->computeNearestGridPoint2(i, j, k, g, mX0, mY0, mZ0);
  int gg = -1;
  if (success) gg = g;
  MPI_Allreduce(&gg, &g, 1, MPI_INT, MPI_MAX, a_EW->m_1d_communicator);
  float_sw4 q, r, s;
  float_sw4 h = a_EW->mGridSize[g];
  bool canBeInverted, curvilinear;
  float_sw4 normwgh[4] = {17.0 / 48.0, 59.0 / 48.0, 43.0 / 48.0, 49.0 / 48.0};

  //   cout << "source location " << i << " " << j << " g= " << g  << " " << mX0
  //   << " "
  //        << mY0 << " " << mZ0 << endl;

  if (g > a_EW->mNumberOfCartesianGrids - 1 && a_EW->topographyExists()) {
    // Curvilinear
    // Problem when the curvilinear mapping is not analytic:
    // This routine can only compute the 's' coordinate if (mX0, mY0) is owned
    // by this processor
    canBeInverted = a_EW->m_gridGenerator->inverse_grid_mapping(
        a_EW, mX0, mY0, mZ0, g, q, r, s);
    // canBeInverted =
    //     a_EW->invert_curvilinear_grid_mapping(mX0, mY0, mZ0, q, r, s);

    // Broadcast the computed s to all processors.
    // First find out the ID of a processor that defines s ...
    int s_owner = -1;
    if (canBeInverted) MPI_Comm_rank(a_EW->m_1d_communicator, &s_owner);
    int s_owner_tmp = s_owner;
    MPI_Allreduce(&s_owner_tmp, &s_owner, 1, MPI_INT, MPI_MAX,
                  a_EW->m_1d_communicator);
    // ...then broadcast s
    if (s_owner > -1)
      MPI_Bcast(
          &s, 1, a_EW->m_mpifloat, s_owner,
          a_EW->m_1d_communicator);  // s_owner is sender, all others receive
    else {
      printf(
          "ERROR in Source::set_grid_point_sources4, no processor could invert "
          "the grid mapping \n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // if s < 0, the source is located above the grid and the call to
    // find_curvilinear_derivatives_at_point will fail
    if (s < 0.) {
      float_sw4 xTop, yTop, zTop;
      a_EW->m_gridGenerator->grid_mapping(a_EW, q, r, 0., g, xTop, yTop, zTop);
      double lat, lon;
      a_EW->computeGeographicCoord(mX0, mY0, lon, lat);
      printf(
          "Found a source above the curvilinear grid! Lat=%e, Lon=%e, source "
          "Z-level = %e, grid boundary Z = %e\n",
          lat, lon, mZ0, zTop);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    curvilinear = true;
    canBeInverted = true;
  } else {
    // Cartesian case
    q = mX0 / h + 1;
    r = mY0 / h + 1;
    s = (mZ0 - a_EW->m_zmin[g]) / h + 1;
    canBeInverted = true;
    curvilinear = false;
  }
  bool curvilineargp1 = g + 1 >= a_EW->mNumberOfCartesianGrids;
  bool curvilineargm1 = g - 1 >= a_EW->mNumberOfCartesianGrids;

  int Ni = a_EW->m_global_nx[g];
  int Nj = a_EW->m_global_ny[g];
  int Nz = a_EW->m_global_nz[g];

  int ic, jc, kc;
  bool upperbndry, lowerbndry, ccbndry, gridrefbndry;
  float_sw4 ai, bi, ci;

  // Delta distribution
  float_sw4 wghi[6], wghj[6], wghk[6], wghix[6], wghjy[6], wghkz[6];
  float_sw4 wghixx[6], wghjyy[6], wghkzz[6];
  // Delta' distribution
  float_sw4 dwghi[6], dwghj[6], dwghk[6], dwghix[6], dwghjy[6], dwghkz[6];
  float_sw4 dwghixx[6], dwghjyy[6], dwghkzz[6];

  // k-weights across mesh refinement boundary
  float_sw4 wghkref[6], dwghkref[6], wghrefkz[6], wghrefkzz[6];

  // at this point canBeInverted == true
  //   if( canBeInverted )
  //   {
  // Compute source location and weights in source discretization
  ic = static_cast<int>(floor(q));
  jc = static_cast<int>(floor(r));
  kc = static_cast<int>(floor(s));

  // Bias stencil away from boundary, no source at ghost/padding points
  if (ic <= 2) ic = 3;
  if (ic >= Ni - 2) ic = Ni - 3;
  if (jc <= 2) jc = 3;
  if (jc >= Nj - 2) jc = Nj - 3;

  // Six point stencil, with points, kc-2,..kc+3, Interior in domain
  // if kc-2>=1, kc+3 <= Nz --> kc >= 3 and kc <= Nz-3
  // Can evaluate with two ghost points if kc-2>=-1 and kc+3 <= Nz+2
  //  --->  kc >=1 and kc <= Nz-1
  //
  if (kc >= Nz) kc = Nz - 1;
  if (kc < 1) kc = 1;

  // upper(surface) and lower boundaries , when the six point stencil
  // kc-2,..kc+3 make use of the first (k=1) or the last (k=Nz) interior point.
  upperbndry = (kc == 1 || kc == 2 || kc == 3);
  lowerbndry = (kc == Nz - 1 || kc == Nz - 2 || kc == Nz - 3);

  // ccbndry=true if at the interface between the curvilinear grid and the
  // cartesian grid. Defined as the six point stencil uses values from both
  // grids.

  //   ccbndry = a_EW->topographyExists() &&  ( (upperbndry && g ==
  //   a_EW->mNumberOfGrids-2) ||
  //                                            (lowerbndry && g ==
  //                                            a_EW->mNumberOfGrids-1)    );

  ccbndry = a_EW->topographyExists() &&
            ((upperbndry && g == a_EW->mNumberOfCartesianGrids - 1) ||
             (lowerbndry && g == a_EW->mNumberOfCartesianGrids));

  // gridrefbndry=true if at the interface between two grids of different
  // refinements.
  gridrefbndry = (upperbndry && g < a_EW->mNumberOfGrids - 1 && !ccbndry) ||
                 (lowerbndry && g > 0 && !ccbndry);

  bool curvilinear_refbndry =
      gridrefbndry && g >= a_EW->mNumberOfCartesianGrids;
  int ncurv = a_EW->mNumberOfGrids - a_EW->mNumberOfCartesianGrids;
  bool cc_ic_bndry =
      ccbndry && !a_EW->m_gridGenerator->curviCartIsSmooth(ncurv);
  //
  // ********* do the filtering of the time function here as needed based on
  // (ic, jc, kc) and gridrefbndry ******
  //

  // AP: Jun 29, 2017: Since interior_point_in_proc only takes the (i,j)
  // indices into account, there is no point looping over
  // k. Furthermore, since this source belongs to this MPI task, it
  // needs to get initialized (filtered), unless it has already been done
  //
  for (int j = jc - 2; j <= jc + 3; j++)
    for (int i = ic - 2; i <= ic + 3; i++) {
      // check if (i,j) belongs to this processor
      if (a_EW->interior_point_in_proc(i, j, g) && !m_timeFuncIsReady) {
        // (optionally) filter and spline interpolate the time function in
        // this Source object, unless already done before
        prepareTimeFunc(a_EW->m_prefilter_sources, a_EW->getTimeStep(),
                        a_EW->getNumberOfTimeSteps(), a_EW->m_filter_ptr);
      }
    }

  // AP: Jun 29, 2017: The purpose of the following code is to make sure
  // the time function is properly initialized. This is important for
  // iDiscrete time functions, which have to be (optionally) filtered
  // and interpolated by a spline before they can be evaluated.
  //
  // Consider a source that is near a processor boundary and a grid
  // refinement boundary, but belongs to a neighboring processor. It is
  // possible that the source function also is needed on this processes,
  // because the 6-point source stencil extends further on a coarser
  // grid. In that case, we must also initialize the time function on
  // this process.
  if (gridrefbndry)  // near MR interface
  {
    // compute (icref, jcref) (duplicated from below)
    int icref, jcref;
    int Niref, Njref;
    float_sw4 qref, rref;
    if (kc - 1 < Nz - kc) {
      // kc closer to upper boundary. Source spread to finer grid above.
      qref = mX0 / (0.5 * h) + 1;
      rref = mY0 / (0.5 * h) + 1;
      Niref = a_EW->m_global_nx[g + 1];
      Njref = a_EW->m_global_ny[g + 1];
    } else {
      // kc closer to lower boundary. Source spread to coarser grid below.
      qref = mX0 / (2 * h) + 1;
      rref = mY0 / (2 * h) + 1;
      Niref = a_EW->m_global_nx[g - 1];
      Njref = a_EW->m_global_ny[g - 1];
    }
    icref = static_cast<int>(qref);
    jcref = static_cast<int>(rref);
    if (icref <= 2) icref = 3;
    if (icref >= Niref - 2) icref = Niref - 3;
    if (jcref <= 2) jcref = 3;
    if (jcref >= Njref - 2) jcref = Njref - 3;
    // done computing (icref, jcref)
    for (int k = kc - 2; k <= kc + 3; k++) {
      if (k <= 1) {
        // Finer grid above
        for (int j = jcref - 2; j <= jcref + 3; j++)
          for (int i = icref - 2; i <= icref + 3; i++) {
            if (a_EW->interior_point_in_proc(i, j, g + 1) &&
                !m_timeFuncIsReady)  // checks if (i,j) belongs to this
                                     // processor
            {
              // filter the time function in this Source object, unless already
              // done
              prepareTimeFunc(a_EW->m_prefilter_sources, a_EW->getTimeStep(),
                              a_EW->getNumberOfTimeSteps(), a_EW->m_filter_ptr);
            }
          }
      }  // end if k<=1
      if (k >= Nz) {
        // Coarser grid below
        for (int j = jcref - 2; j <= jcref + 3; j++)
          for (int i = icref - 2; i <= icref + 3; i++) {
            if (a_EW->interior_point_in_proc(i, j, g - 1) &&
                !m_timeFuncIsReady)  // checks if (i,j) belongs to this
                                     // processor
            {
              // filter the time function in this Source object, unless already
              // done
              prepareTimeFunc(a_EW->m_prefilter_sources, a_EW->getTimeStep(),
                              a_EW->getNumberOfTimeSteps(), a_EW->m_filter_ptr);
            }
          }
      }  // end if k >= Nz

    }  // end for kc
  }    // end if near MR interface

  // If not at the interface between different grids, bias stencil away
  // from the boundary.
  if (!ccbndry && !gridrefbndry) {
    if (kc <= 2) kc = 3;
    if (kc >= Nz - 2) kc = Nz - 3;
  }
  ai = q - ic, bi = r - jc, ci = s - kc;

  // Delta distribution
  getsourcewgh(ai, wghi, wghix, wghixx);
  getsourcewgh(bi, wghj, wghjy, wghjyy);
  getsourcewgh(ci, wghk, wghkz, wghkzz);

  // Delta' distribution
  getsourcedwgh(ai, dwghi, dwghix, dwghixx);
  getsourcedwgh(bi, dwghj, dwghjy, dwghjyy);
  getsourcedwgh(ci, dwghk, dwghkz, dwghkzz);

  // Special boundary stencil at free surface
  if (!ccbndry && !gridrefbndry && (kc == 3 && ci <= 0)) {
    getsourcewghlow(ci, wghk, wghkz, wghkzz);
    getsourcedwghlow(ci, dwghk, dwghkz, dwghkzz);
  }

  // Special source discretization across grid refinement boundary
  if (gridrefbndry) {
    float_sw4 sw = 1.0 / 3;
    if (lowerbndry) {
      if (kc == Nz - 1) {
        getsourcedwghNM1sm6(ci, dwghk);
        getsourcewghNM1sm6(ci, wghk);
        //            wghkref[3]  = wghk[3]*0.5;
        //            wghk[3]     = wghk[3]*0.5;
        wghkref[3] = wghk[3] * (1 - sw);
        wghk[3] = wghk[3] * sw;

        wghkref[4] = wghk[4];
        wghkref[5] = wghk[5];

        //            dwghkref[3] = dwghk[3]*0.5;
        //            dwghk[3]    = dwghk[3]*0.5;
        dwghkref[3] = dwghk[3] * (1 - sw);
        dwghk[3] = dwghk[3] * sw;

        dwghkref[4] = dwghk[4];
        dwghkref[5] = dwghk[5];

        wghkref[3] /= normwgh[0];
        wghkref[4] /= normwgh[1];
        wghkref[5] /= normwgh[2];
        wghk[3] /= normwgh[0];
        wghk[2] /= normwgh[1];
        wghk[1] /= normwgh[2];
        wghk[0] /= normwgh[3];
        dwghkref[3] /= normwgh[0];
        dwghkref[4] /= normwgh[1];
        dwghkref[5] /= normwgh[2];
        dwghk[3] /= normwgh[0];
        dwghk[2] /= normwgh[1];
        dwghk[1] /= normwgh[2];
        dwghk[0] /= normwgh[3];
      } else if (kc == Nz - 2) {
        getsourcedwghNM2sm6(ci, dwghk);
        getsourcewghNM2sm6(ci, wghk);
        //            wghkref[4]  = wghk[4]*0.5;
        //            wghk[4]     = wghk[4]*0.5;
        wghkref[4] = wghk[4] * (1 - sw);
        wghk[4] = wghk[4] * sw;

        wghkref[5] = wghk[5];

        //            dwghkref[4] = dwghk[4]*0.5;
        //            dwghk[4]    = dwghk[4]*0.5;
        dwghkref[4] = dwghk[4] * (1 - sw);
        dwghk[4] = dwghk[4] * sw;
        dwghkref[5] = dwghk[5];

        wghkref[4] /= normwgh[0];
        wghkref[5] /= normwgh[1];
        wghk[4] /= normwgh[0];
        wghk[3] /= normwgh[1];
        wghk[2] /= normwgh[2];
        wghk[1] /= normwgh[3];
        dwghkref[4] /= normwgh[0];
        dwghkref[5] /= normwgh[1];
        dwghk[4] /= normwgh[0];
        dwghk[3] /= normwgh[1];
        dwghk[2] /= normwgh[2];
        dwghk[1] /= normwgh[3];
      } else if (kc == Nz - 3) {
        getsourcedwgh(ci, dwghk, wghrefkz, wghrefkzz);
        getsourcewgh(ci, wghk, wghrefkz, wghrefkzz);
        //            wghkref[5]  = wghk[5]*0.5;
        //            wghk[5]     = wghk[5]*0.5;
        //            dwghkref[5] = dwghk[5]*0.5;
        //            dwghk[5]    = dwghk[5]*0.5;
        wghkref[5] = wghk[5] * (1 - sw);
        wghk[5] = wghk[5] * sw;
        dwghkref[5] = dwghk[5] * (1 - sw);
        dwghk[5] = dwghk[5] * sw;
        //	       cout << " sumwgh = " <<
        // dwghk[0]+dwghk[1]+dwghk[2]+dwghk[3]+dwghk[4]+dwghk[5]+dwghkref[5] <<
        // endl;

        wghkref[5] /= normwgh[0];
        wghk[5] /= normwgh[0];
        wghk[4] /= normwgh[1];
        wghk[3] /= normwgh[2];
        wghk[2] /= normwgh[3];
        dwghkref[5] /= normwgh[0];
        dwghk[5] /= normwgh[0];
        dwghk[4] /= normwgh[1];
        dwghk[3] /= normwgh[2];
        dwghk[2] /= normwgh[3];
      }
    } else {
      if (kc == 1) {
        getsourcedwghNsm6(2 * ci, dwghk);
        getsourcewghNsm6(2 * ci, wghk);
        wghkref[0] = wghk[0];
        wghkref[1] = wghk[1];
        wghkref[2] = wghk[2] * sw;
        wghk[2] = wghk[2] * (1 - sw);

        dwghkref[0] = dwghk[0];
        dwghkref[1] = dwghk[1];
        dwghkref[2] = dwghk[2] * sw;
        dwghk[2] = dwghk[2] * (1 - sw);

        //	       cout << "kc = 1  ref:  " << wghkref[0] << " " <<
        // wghkref[1] << " " << wghkref[2] << endl; 	       cout << " this: "
        // << wghk[2] << " " << wghk[3] << " " << wghk[4] << " " << wghk[5] <<
        // endl; 	       cout << "  middle sum: " << wghk[2]+wghkref[2] <<
        // endl; 	       cout <<
        //" 2*ci = " << 2*ci << endl;

        wghkref[0] /= normwgh[2];
        wghkref[1] /= normwgh[1];
        wghkref[2] /= normwgh[0];
        wghk[2] /= normwgh[0];
        wghk[3] /= normwgh[1];
        wghk[4] /= normwgh[2];
        wghk[5] /= normwgh[3];

        dwghkref[0] /= normwgh[2];
        dwghkref[1] /= normwgh[1];
        dwghkref[2] /= normwgh[0];
        dwghk[2] /= normwgh[0];
        dwghk[3] /= normwgh[1];
        dwghk[4] /= normwgh[2];
        dwghk[5] /= normwgh[3];

        dwghk[2] *= 2;
        dwghk[3] *= 2;
        dwghk[4] *= 2;
        dwghk[5] *= 2;
      } else if (kc == 2) {
        getsourcedwghP1sm6(2 * ci, dwghk);
        getsourcewghP1sm6(2 * ci, wghk);

        wghkref[0] = wghk[0];
        wghkref[1] = wghk[1] * sw;
        wghk[1] = wghk[1] * (1 - sw);

        dwghkref[0] = dwghk[0];
        dwghkref[1] = dwghk[1] * sw;
        dwghk[1] = dwghk[1] * (1 - sw);

        wghkref[0] /= normwgh[1];
        wghkref[1] /= normwgh[0];
        wghk[1] /= normwgh[0];
        wghk[2] /= normwgh[1];
        wghk[3] /= normwgh[2];
        wghk[4] /= normwgh[3];

        dwghkref[0] /= normwgh[1];
        dwghkref[1] /= normwgh[0];
        dwghk[1] /= normwgh[0];
        dwghk[2] /= normwgh[1];
        dwghk[3] /= normwgh[2];
        dwghk[4] /= normwgh[3];

        dwghk[1] *= 2;
        dwghk[2] *= 2;
        dwghk[3] *= 2;
        dwghk[4] *= 2;
        dwghk[5] *= 2;

      } else if (kc == 3) {
        getsourcedwgh(ci, dwghk, wghrefkz, wghrefkzz);
        for (int k = 0; k <= 5; k++) dwghk[k] *= 0.5;
        getsourcewgh(ci, wghk, wghrefkz, wghrefkzz);
        wghkref[0] = wghk[0] * sw;
        wghk[0] = wghk[0] * (1 - sw);
        dwghkref[0] = dwghk[0] * sw;
        dwghk[0] = dwghk[0] * (1 - sw);

        wghkref[0] /= normwgh[0];
        wghk[0] /= normwgh[0];
        wghk[1] /= normwgh[1];
        wghk[2] /= normwgh[2];
        wghk[3] /= normwgh[3];
        dwghkref[0] /= normwgh[0];
        dwghk[0] /= normwgh[0];
        dwghk[1] /= normwgh[1];
        dwghk[2] /= normwgh[2];
        dwghk[3] /= normwgh[3];
        dwghk[0] *= 2;
        dwghk[1] *= 2;
        dwghk[2] *= 2;
        dwghk[3] *= 2;
        dwghk[4] *= 2;
        dwghk[5] *= 2;
      }
    }
  }

  // Boundary correction, at upper boundary, but only if SBP operators are used
  // there
  //      if( !gridrefbndry && (g == a_EW->mNumberOfGrids-1) &&
  //      a_EW->is_onesided(g,4)  )
  if (!gridrefbndry && a_EW->is_onesided(g, 4)) {
    for (int k = 0; k <= 5; k++) {
      if ((1 <= k + kc - 2) && (k + kc - 2 <= 4)) {
        wghk[k] /= normwgh[k + kc - 3];
        dwghk[k] /= normwgh[k + kc - 3];
        wghkz[k] /= normwgh[k + kc - 3];
        dwghkz[k] /= normwgh[k + kc - 3];
        wghkzz[k] /= normwgh[k + kc - 3];
        dwghkzz[k] /= normwgh[k + kc - 3];
      }
    }
  }
  if (!gridrefbndry && a_EW->is_onesided(g, 5)) {
    for (int k = 0; k <= 5; k++) {
      if ((Nz - 3 <= k + kc - 2) && (k + kc - 2 <= Nz)) {
        wghk[k] /= normwgh[Nz - k - kc + 2];
        dwghk[k] /= normwgh[Nz - k - kc + 2];
        wghkz[k] /= normwgh[Nz - k - kc + 2];
        dwghkz[k] /= normwgh[Nz - k - kc + 2];
        wghkzz[k] /= normwgh[Nz - k - kc + 2];
        dwghkzz[k] /= normwgh[Nz - k - kc + 2];
      }
    }
  }

  int myid;
  MPI_Comm_rank(a_EW->m_1d_communicator, &myid);
  //   cout << myid << " SOURCE at " << ic << " " << jc << " "  << kc ;
  //   if( canBeInverted )
  //      cout << " can be inverted";
  //   else
  //      cout << " can not be inverted";

  // If source at grid refinement interface, set up variables for
  // discretization on grid on the other side of the interface
  int icref, jcref;
  float_sw4 airef, biref, wghiref[6], wghirefx[6], wghirefxx[6];
  float_sw4 wghjref[6], wghjrefy[6], wghjrefyy[6];
  float_sw4 dwghiref[6], dwghjref[6];
  if (gridrefbndry) {
    int Niref, Njref;
    float_sw4 qref, rref;
    if (kc - 1 < Nz - kc) {
      // kc closer to upper boundary. Source spread to finer grid above.
      qref = mX0 / (0.5 * h) + 1;
      rref = mY0 / (0.5 * h) + 1;
      Niref = a_EW->m_global_nx[g + 1];
      Njref = a_EW->m_global_ny[g + 1];
    } else {
      // kc closer to lower boundary. Source spread to coarser grid below.
      qref = mX0 / (2 * h) + 1;
      rref = mY0 / (2 * h) + 1;
      Niref = a_EW->m_global_nx[g - 1];
      Njref = a_EW->m_global_ny[g - 1];
    }
    icref = static_cast<int>(qref);
    jcref = static_cast<int>(rref);
    if (icref <= 2) icref = 3;
    if (icref >= Niref - 2) icref = Niref - 3;
    if (jcref <= 2) jcref = 3;
    if (jcref >= Njref - 2) jcref = Njref - 3;
    airef = qref - icref;
    biref = rref - jcref;
    getsourcewgh(airef, wghiref, wghirefx, wghirefxx);
    getsourcewgh(biref, wghjref, wghjrefy, wghjrefyy);
    // reuse wghirefx,wghirefxx, these are assumed not to be used with grid
    // refinement.
    getsourcedwgh(airef, dwghiref, wghirefx, wghirefxx);
    getsourcedwgh(biref, dwghjref, wghjrefy, wghjrefyy);
  }

  // Point source. NOTE: Derivatives needed for source inversion not implemented
  // for this case.
  if (!mIsMomentSource /*&& canBeInverted*/) {
    if (curvilinear_refbndry) {
      //         get_mr_psources(  a_EW, g, ic, jc, kc, icref, jcref, normwgh,
      //         point_sources );
      get_mr_psources(a_EW, g, q, r, s, false, normwgh, point_sources);
    } else if (cc_ic_bndry) {
      get_cc_psources(a_EW, g, q, r, s, false, normwgh, point_sources);
    } else {
      for (int k = kc - 2; k <= kc + 3; k++)
        for (int j = jc - 2; j <= jc + 3; j++)
          for (int i = ic - 2; i <= ic + 3; i++) {
            float_sw4 wF =
                wghi[i - ic + 2] * wghj[j - jc + 2] * wghk[k - kc + 2];
            if ((wF != 0) &&
                (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) &&
                a_EW->interior_point_in_proc(
                    i, j, g))  // checks if (i,j) belongs to this processor
            {
              if (curvilinear)
                wF /= a_EW->mJ[g](i, j, k);
              else
                wF /= h * h * h;

              if (1 <= k && k <= Nz) {
                GridPointSource* sourcePtr = new GridPointSource(
                    mFreq, mT0, i, j, k, g, wF * mForces[0], wF * mForces[1],
                    wF * mForces[2], mTimeDependence, mNcyc, mPar, mNpar, mIpar,
                    mNipar);
                point_sources.push_back(sourcePtr);
              }
              if (k <= 1 && ccbndry && upperbndry) {
                int Nzp = a_EW->m_global_nz[g + 1];
                int kk = Nzp - 1 + k;
                float_sw4 wF =
                    wghi[i - ic + 2] * wghj[j - jc + 2] * wghk[k - kc + 2];
                if (curvilineargp1)
                  wF /= a_EW->mJ[g + 1](i, j, kk);
                else
                  wF /= 0.125 * h * h * h;

                GridPointSource* sourcePtr = new GridPointSource(
                    mFreq, mT0, i, j, kk, g + 1, wF * mForces[0],
                    wF * mForces[1], wF * mForces[2], mTimeDependence, mNcyc,
                    mPar, mNpar, mIpar, mNipar);
                point_sources.push_back(sourcePtr);
              }
              if (k >= Nz && ccbndry && lowerbndry) {
                int kk = k - Nz + 1;
                float_sw4 wF =
                    wghi[i - ic + 2] * wghj[j - jc + 2] * wghk[k - kc + 2];
                if (curvilineargm1)
                  wF /= a_EW->mJ[g - 1](i, j, kk);
                else
                  wF /= 8 * h * h * h;
                GridPointSource* sourcePtr = new GridPointSource(
                    mFreq, mT0, i, j, kk, g - 1, wF * mForces[0],
                    wF * mForces[1], wF * mForces[2], mTimeDependence, mNcyc,
                    mPar, mNpar, mIpar, mNipar);
                point_sources.push_back(sourcePtr);
              }
            }
          }
      if (gridrefbndry) {
        for (int k = kc - 2; k <= kc + 3; k++) {
          if (k <= 1) {
            // Finer grid above
            for (int j = jcref - 2; j <= jcref + 3; j++)
              for (int i = icref - 2; i <= icref + 3; i++) {
                float_sw4 wF = wghiref[i - icref + 2] * wghjref[j - jcref + 2] *
                               wghkref[k - kc + 2];
                if ((wF != 0) &&
                    (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) &&
                    a_EW->interior_point_in_proc(
                        i, j,
                        g + 1))  // checks if (i,j) belongs to this processor
                {
                  int Nzp = a_EW->m_global_nz[g + 1];
                  int kk = Nzp - 1 + k;
                  if (curvilineargp1)
                    wF /= a_EW->mJ[g + 1](i, j, kk);
                  else
                    wF /= 0.125 * h * h * h;
                  GridPointSource* sourcePtr = new GridPointSource(
                      mFreq, mT0, i, j, kk, g + 1, wF * mForces[0],
                      wF * mForces[1], wF * mForces[2], mTimeDependence, mNcyc,
                      mPar, mNpar, mIpar, mNipar);
                  point_sources.push_back(sourcePtr);
                }
              }
          }  // end if k<=1
          if (k >= Nz) {
            // Coarser grid below
            for (int j = jcref - 2; j <= jcref + 3; j++)
              for (int i = icref - 2; i <= icref + 3; i++) {
                float_sw4 wF = wghiref[i - icref + 2] * wghjref[j - jcref + 2] *
                               wghkref[k - kc + 2];
                if ((wF != 0) &&
                    (mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0) &&
                    a_EW->interior_point_in_proc(
                        i, j,
                        g - 1))  // checks if (i,j) belongs to this processor
                {
                  int kk = k - Nz + 1;
                  if (curvilineargm1)
                    wF /= a_EW->mJ[g - 1](i, j, kk);
                  else
                    wF /= 8 * h * h * h;
                  GridPointSource* sourcePtr = new GridPointSource(
                      mFreq, mT0, i, j, kk, g - 1, wF * mForces[0],
                      wF * mForces[1], wF * mForces[2], mTimeDependence, mNcyc,
                      mPar, mNpar, mIpar, mNipar);
                  point_sources.push_back(sourcePtr);
                }
              }
          }
        }
      }  // Grid refinement boundary
    }
  }  // if !mIsMomentSource (i.e. pointForce)

  else if (mIsMomentSource)  // Moment source.
  {
    if (curvilinear_refbndry) {
      //         get_mr_psources(  a_EW, g, ic, jc, kc, icref, jcref, normwgh,
      //         point_sources );
      get_mr_psources(a_EW, g, q, r, s, true, normwgh, point_sources);
    } else if (cc_ic_bndry) {
      get_cc_psources(a_EW, g, q, r, s, true, normwgh, point_sources);
    } else {
      float_sw4 qX0[3], rX0[3], sX0[3];

      // Gradients of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s)
      float_sw4 dsX0[3], dsY0[3], dsZ0[3];
      // Hessians of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s), in order
      // qq,qr,qs,rr,rs,ss
      float_sw4 d2sX0[6], d2sY0[6], d2sZ0[6];
      if (!curvilinear) {
        // Cartesian case, constant metric
        qX0[0] = 1 / h;
        qX0[1] = 0;
        qX0[2] = 0;
        rX0[0] = 0;
        rX0[1] = 1 / h;
        rX0[2] = 0;
        sX0[0] = 0;
        sX0[1] = 0;
        sX0[2] = 1 / h;
        dsX0[0] = dsX0[1] = dsX0[2] = 0;
        dsY0[0] = dsY0[1] = dsY0[2] = 0;
        dsZ0[0] = dsZ0[1] = dsZ0[2] = 0;
        for (int i = 0; i < 6; i++) d2sX0[i] = d2sY0[i] = d2sZ0[i] = 0;
      } else {
        // Compute the curvilinear metric in the processor that owns the source.
        //   (ic, jc are undefined if canBeInverted is false.)
        float_sw4 zdertmp[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, zq, zr, zs;
        float_sw4 zqq, zqr, zqs, zrr, zrs, zss;
        if (a_EW->interior_point_in_proc(ic, jc, g) /*&& canBeInverted*/) {
          // TESTING grid generator functions:
          //            float_sw4 xx,yy,zz;
          //            a_EW->m_gridGenerator->grid_mapping(a_EW,q,r,s,g,xx,yy,zz);

          //            a_EW->m_gridGenerator->grid_mapping_diff( a_EW, q, r, s,
          //            g, ic, jc, kc,
          //                                                      zq, zr, zs,
          //                                                      zqq, zqr, zqs,
          //                                                       zrr, zrs, zss
          //                                                       );
          // Second derivatives
          //            float_sw4 h=1e-3, zqh, zrh, zsh, zqqh, zqrh, zqsh, zrrh,
          //            zrsh, zssh; a_EW->m_gridGenerator->grid_mapping_diff(
          //            a_EW, q+h, r, s, g, ic, jc, kc,
          //                                                      zqh, zrh, zsh,
          //                                                      zqqh, zqrh,
          //                                                      zqsh, zrrh,
          //                                                      zrsh, zssh );
          //            cout << " zqqa= " << zqq << " zqqn= " << (zqh-zq)/h <<
          //            endl; cout << " zqra= " << zqr << " zqrn= " <<
          //            (zrh-zr)/h << endl; cout << " zqsa= " << zqs << " zqsn=
          //            " << (zsh-zs)/h << endl;
          //            a_EW->m_gridGenerator->grid_mapping_diff( a_EW, q, r+h,
          //            s, g, ic, jc, kc,
          //                                                      zqh, zrh, zsh,
          //                                                      zqqh, zqrh,
          //                                                      zqsh, zrrh,
          //                                                      zrsh, zssh );
          //            cout << " zrra= " << zrr << " zrrn= " << (zrh-zr)/h <<
          //            endl; cout << " zrsa= " << zrs << " zrsn= " <<
          //            (zsh-zs)/h << endl;
          //            a_EW->m_gridGenerator->grid_mapping_diff( a_EW, q, r,
          //            s+h, g, ic, jc, kc,
          //                                        zqh, zrh, zsh, zqqh, zqrh,
          //                                        zqsh, zrrh, zrsh, zssh );
          //            cout << " zssa= " << zss << " zssn= " << (zsh-zs)/h <<
          //            endl;                    // First derivatives
          //            a_EW->m_gridGenerator->grid_mapping(a_EW,q+h,r,s,g,xxh,yyh,zzh);
          //            cout << " zqa= " << zq << " zqn= " << (zzh-zz)/h <<
          //            endl;
          //            a_EW->m_gridGenerator->grid_mapping(a_EW,q,r+h,s,g,xxh,yyh,zzh);
          //            cout << " zra= " << zr << " zrn= " << (zzh-zz)/h <<
          //            endl;
          //            a_EW->m_gridGenerator->grid_mapping(a_EW,q,r,s+h,g,xxh,yyh,zzh);
          //            cout << " zsa= " << zs << " zsn= " << (zzh-zz)/h <<
          //            endl;
          //
          //            compute_metric_at_source( a_EW, q, r, s, ic, jc, kc, g,
          //            zq, zr, zs,
          //                                      zqq, zqr, zqs, zrr, zrs, zss
          //                                      );
          //            cout << "Old zq,zr,zs " << zq << " " << zr << " " << zs
          //            << endl; cout << endl;
          // END testing
          a_EW->m_gridGenerator->grid_mapping_diff(a_EW, q, r, s, g, ic, jc, kc,
                                                   zq, zr, zs, zqq, zqr, zqs,
                                                   zrr, zrs, zss);
          zdertmp[0] = zq;
          zdertmp[1] = zr;
          zdertmp[2] = zs;
          zdertmp[3] = zqq;
          zdertmp[4] = zqr;
          zdertmp[5] = zqs;
          zdertmp[6] = zrr;
          zdertmp[7] = zrs;
          zdertmp[8] = zss;
        }
        //	 // Broadcast the computed metric to all processors.
        //	 // First find out the ID of the processor that computed the
        //metric... 	 int owner = -1; 	 if( a_EW->interior_point_in_proc( ic, jc, g
        //) && canBeInverted )
        //            MPI_Comm_rank(MPI_COMM_WORLD, &owner );
        //	 int owntmp = owner;
        //	 MPI_Allreduce( &owntmp, &owner, 1, MPI_INT, MPI_MAX,
        //MPI_COMM_WORLD );
        //	 // ...then broadcast the derivatives
        //	 MPI_Bcast( zder, 3, MPI_DOUBLE, owner, MPI_COMM_WORLD );

        // Simpler solution
        float_sw4 zder[9];
        MPI_Allreduce(zdertmp, zder, 9, a_EW->m_mpifloat, MPI_SUM,
                      a_EW->m_1d_communicator);
        zq = zder[0];
        zr = zder[1];
        zs = zder[2];
        zqq = zder[3];
        zqr = zder[4];
        zqs = zder[5];
        zrr = zder[6];
        zrs = zder[7];
        zss = zder[8];

        qX0[0] = 1 / h;
        qX0[1] = 0;
        qX0[2] = 0;
        rX0[0] = 0;
        rX0[1] = 1 / h;
        rX0[2] = 0;
        sX0[0] = -zq / (h * zs);
        sX0[1] = -zr / (h * zs);
        sX0[2] = 1 / zs;

        float_sw4 deni = 1 / (h * zs * zs);
        dsX0[0] = -(zqq * zs - zqs * zq) * deni;
        dsX0[1] = -(zqr * zs - zrs * zq) * deni;
        dsX0[2] = -(zqs * zs - zss * zq) * deni;

        dsY0[0] = -(zqr * zs - zqs * zr) * deni;
        dsY0[1] = -(zrr * zs - zrs * zr) * deni;
        dsY0[2] = -(zrs * zs - zss * zr) * deni;

        deni *= h;
        dsZ0[0] = -zqs * deni;
        dsZ0[1] = -zrs * deni;
        dsZ0[2] = -zss * deni;
      }  // end curvilinear

      // Gradients of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s)
      // NYI
      //            double dsX0[3], dsY0[3], dsZ0[3], d2sX0[6], d2sY0[6],
      //            d2sZ0[6];
      //	    dsX0[0] = 0;
      //	    dsX0[1] = 0;
      //	    dsX0[2] = 0;

      //	    dsY0[0] = 0;
      //	    dsY0[1] = 0;
      //	    dsY0[2] = 0;

      //	    dsZ0[0] = 0;
      //	    dsZ0[1] = 0;
      //	    dsZ0[2] = 0;
      // Hessians of sX0[0]=sX, sX0[1]=sY, and sX0[2]=sZ wrt. (q,r,s), in order
      // qq,qr,qs,rr,rs,ss
      //	    d2sX0[0] = 0;
      //	    d2sX0[1] =0;
      //	    d2sX0[2] =0;
      //	    d2sX0[3] =0;
      //	    d2sX0[4] =0;
      //	    d2sX0[5] =0;

      //	    d2sY0[0] =0;
      //	    d2sY0[1] =0;
      //	    d2sY0[2] =0;
      //	    d2sY0[3] =0;
      //	    d2sY0[4] =0;
      //	    d2sY0[5] =0;

      //	    d2sZ0[0] =0;
      //	    d2sZ0[1] =0;
      //	    d2sZ0[2] =0;
      //	    d2sZ0[3] =0;
      //	    d2sZ0[4] =0;
      //	    d2sZ0[5] =0;

      //      if( canBeInverted )
      //   {
      for (int k = kc - 2; k <= kc + 3; k++)
        for (int j = jc - 2; j <= jc + 3; j++)
          for (int i = ic - 2; i <= ic + 3; i++) {
            float_sw4 wFx = 0, wFy = 0, wFz = 0, dsdp[27];
            if (a_EW->interior_point_in_proc(i, j, g)) {
              //                     cout << " src at " << i << " " << j << " "
              //                     << k << endl;
              wFx += qX0[0] * dwghi[i - ic + 2] * wghj[j - jc + 2] *
                     wghk[k - kc + 2];
              //		  wFy += qX0[1]*dwghi[i-ic+2]* wghj[j-jc+2]*
              //wghk[k-kc+2]; 		  wFz += qX0[2]*dwghi[i-ic+2]* wghj[j-jc+2]*
              //wghk[k-kc+2];

              //		  wFx +=  wghi[i-ic+2]*rX0[0]*dwghj[j-jc+2]*
              //wghk[k-kc+2];
              wFy += wghi[i - ic + 2] * rX0[1] * dwghj[j - jc + 2] *
                     wghk[k - kc + 2];
              //		  wFz +=  wghi[i-ic+2]*rX0[2]*dwghj[j-jc+2]*
              //wghk[k-kc+2];

              wFx += wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[0] *
                     dwghk[k - kc + 2];
              //		     wFx +=  sX0[0];
              wFy += wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[1] *
                     dwghk[k - kc + 2];
              wFz += wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[2] *
                     dwghk[k - kc + 2];

              float_sw4 hi = 1.0 / h;
              float_sw4 hi2 = hi * hi;
              float_sw4 wFxdx0 = dwghix[i - ic + 2] * wghj[j - jc + 2] *
                                 wghk[k - kc + 2] * hi * qX0[0];
              float_sw4 wFxdy0 = dwghi[i - ic + 2] * wghjy[j - jc + 2] *
                                 wghk[k - kc + 2] * hi * rX0[1];
              //                     float_sw4 wFxdy0 = 0;
              float_sw4 wFxdz0 = dwghi[i - ic + 2] * wghj[j - jc + 2] *
                                 wghkz[k - kc + 2] * hi * sX0[2];

              float_sw4 wFydx0 = wghix[i - ic + 2] * dwghj[j - jc + 2] *
                                 wghk[k - kc + 2] * hi * qX0[0];
              float_sw4 wFydy0 = wghi[i - ic + 2] * dwghjy[j - jc + 2] *
                                 wghk[k - kc + 2] * hi * rX0[1];
              float_sw4 wFydz0 = wghi[i - ic + 2] * dwghj[j - jc + 2] *
                                 wghkz[k - kc + 2] * hi * sX0[2];

              float_sw4 wFzdx0 = wghix[i - ic + 2] * wghj[j - jc + 2] *
                                 dwghk[k - kc + 2] * sX0[2] * qX0[0];
              float_sw4 wFzdy0 = wghi[i - ic + 2] * wghjy[j - jc + 2] *
                                 dwghk[k - kc + 2] * sX0[2] * rX0[1];
              float_sw4 wFzdz0 = wghi[i - ic + 2] * wghj[j - jc + 2] *
                                 dwghkz[k - kc + 2] * sX0[2] * sX0[2];

              if (curvilinear && kc <= Nz - 3) {
                wFxdx0 += dwghi[i - ic + 2] * wghj[j - jc + 2] *
                              wghkz[k - kc + 2] * hi * sX0[0] +
                          wghix[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] * hi * sX0[0] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghkz[k - kc + 2] * sX0[0] * sX0[0] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] *
                              (hi * dsX0[0] + sX0[0] * dsX0[2]);

                wFxdy0 += dwghi[i - ic + 2] * wghj[j - jc + 2] *
                              wghkz[k - kc + 2] * hi * sX0[1] +
                          wghi[i - ic + 2] * wghjy[j - jc + 2] *
                              dwghk[k - kc + 2] * hi * sX0[0] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghkz[k - kc + 2] * sX0[0] * sX0[1] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] *
                              (hi * dsX0[1] + sX0[1] * dsX0[2]);
                //			wFxdy0 += (hi*dsX0[1]+sX0[1]*dsX0[2]);

                wFxdz0 += wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghkz[k - kc + 2] * sX0[0] * sX0[2] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] * sX0[2] * dsX0[2];

                wFydx0 += wghi[i - ic + 2] * dwghj[j - jc + 2] *
                              wghkz[k - kc + 2] * hi * sX0[0] +
                          wghix[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] * hi * sX0[1] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghkz[k - kc + 2] * sX0[0] * sX0[1] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] *
                              (hi * dsY0[0] + sX0[0] * dsY0[2]);

                wFydy0 += wghi[i - ic + 2] * dwghj[j - jc + 2] *
                              wghkz[k - kc + 2] * hi * sX0[1] +
                          wghi[i - ic + 2] * wghjy[j - jc + 2] *
                              dwghk[k - kc + 2] * hi * sX0[1] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghkz[k - kc + 2] * sX0[1] * sX0[1] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] *
                              (hi * dsY0[1] + sX0[1] * dsY0[2]);

                wFydz0 += wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghkz[k - kc + 2] * sX0[1] * sX0[2] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] * sX0[2] * dsY0[2];

                wFzdx0 += wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghkz[k - kc + 2] * sX0[0] * sX0[2] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] *
                              (hi * dsZ0[0] + sX0[0] * dsZ0[2]);

                wFzdy0 += wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghkz[k - kc + 2] * sX0[1] * sX0[2] +
                          wghi[i - ic + 2] * wghj[j - jc + 2] *
                              dwghk[k - kc + 2] *
                              (hi * dsZ0[1] + sX0[1] * dsZ0[2]);

                wFzdz0 += wghi[i - ic + 2] * wghj[j - jc + 2] *
                          dwghk[k - kc + 2] * (sX0[2] * dsZ0[2]);
              }
              // NOTE:  Source second derivatives wrt. (x0,y0,z0) currently not
              // yet implemented for curvilinear grids.

              // Second derivatives

              float_sw4 wFxdx0dx0 = dwghixx[i - ic + 2] * wghj[j - jc + 2] *
                                    wghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFxdx0dy0 = dwghix[i - ic + 2] * wghjy[j - jc + 2] *
                                    wghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFxdx0dz0 = dwghix[i - ic + 2] * wghj[j - jc + 2] *
                                    wghkz[k - kc + 2] * hi2 * hi;
              float_sw4 wFxdy0dy0 = dwghi[i - ic + 2] * wghjyy[j - jc + 2] *
                                    wghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFxdy0dz0 = dwghi[i - ic + 2] * wghjy[j - jc + 2] *
                                    wghkz[k - kc + 2] * hi2 * hi;
              float_sw4 wFxdz0dz0 = dwghi[i - ic + 2] * wghj[j - jc + 2] *
                                    wghkzz[k - kc + 2] * hi2 * hi;

              float_sw4 wFydx0dx0 = wghixx[i - ic + 2] * dwghj[j - jc + 2] *
                                    wghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFydx0dy0 = wghix[i - ic + 2] * dwghjy[j - jc + 2] *
                                    wghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFydx0dz0 = wghix[i - ic + 2] * dwghj[j - jc + 2] *
                                    wghkz[k - kc + 2] * hi2 * hi;
              float_sw4 wFydy0dy0 = wghi[i - ic + 2] * dwghjyy[j - jc + 2] *
                                    wghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFydy0dz0 = wghi[i - ic + 2] * dwghjy[j - jc + 2] *
                                    wghkz[k - kc + 2] * hi2 * hi;
              float_sw4 wFydz0dz0 = wghi[i - ic + 2] * dwghj[j - jc + 2] *
                                    wghkzz[k - kc + 2] * hi2 * hi;

              float_sw4 wFzdx0dx0 = wghixx[i - ic + 2] * wghj[j - jc + 2] *
                                    dwghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFzdx0dy0 = wghix[i - ic + 2] * wghjy[j - jc + 2] *
                                    dwghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFzdx0dz0 = wghix[i - ic + 2] * wghj[j - jc + 2] *
                                    dwghkz[k - kc + 2] * hi2 * hi;
              float_sw4 wFzdy0dy0 = wghi[i - ic + 2] * wghjyy[j - jc + 2] *
                                    dwghk[k - kc + 2] * hi2 * hi;
              float_sw4 wFzdy0dz0 = wghi[i - ic + 2] * wghjy[j - jc + 2] *
                                    dwghkz[k - kc + 2] * hi2 * hi;
              float_sw4 wFzdz0dz0 = wghi[i - ic + 2] * wghj[j - jc + 2] *
                                    dwghkzz[k - kc + 2] * hi2 * hi;
              //                  if( curvilinear && kc <= Nz-3 )
              //		  {

              //		  }

              float_sw4 jaci;
              if (curvilinear)
                jaci = 1 / a_EW->mJ[g](i, j, k);
              else
                jaci = 1.0 / (h * h * h);

              float_sw4 fx =
                  -(mForces[0] * wFx + mForces[1] * wFy + mForces[2] * wFz) *
                  jaci;
              float_sw4 fy =
                  -(mForces[1] * wFx + mForces[3] * wFy + mForces[4] * wFz) *
                  jaci;
              float_sw4 fz =
                  -(mForces[2] * wFx + mForces[4] * wFy + mForces[5] * wFz) *
                  jaci;

              // Derivatives with respect to (x0,y0,z0,mxx,mxy,mxz,myy,myz,mzz)
              dsdp[0] = -(mForces[0] * wFxdx0 + mForces[1] * wFydx0 +
                          mForces[2] * wFzdx0) *
                        jaci;
              dsdp[1] = -(mForces[1] * wFxdx0 + mForces[3] * wFydx0 +
                          mForces[4] * wFzdx0) *
                        jaci;
              dsdp[2] = -(mForces[2] * wFxdx0 + mForces[4] * wFydx0 +
                          mForces[5] * wFzdx0) *
                        jaci;
              dsdp[3] = -(mForces[0] * wFxdy0 + mForces[1] * wFydy0 +
                          mForces[2] * wFzdy0) *
                        jaci;
              dsdp[4] = -(mForces[1] * wFxdy0 + mForces[3] * wFydy0 +
                          mForces[4] * wFzdy0) *
                        jaci;
              dsdp[5] = -(mForces[2] * wFxdy0 + mForces[4] * wFydy0 +
                          mForces[5] * wFzdy0) *
                        jaci;
              dsdp[6] = -(mForces[0] * wFxdz0 + mForces[1] * wFydz0 +
                          mForces[2] * wFzdz0) *
                        jaci;
              dsdp[7] = -(mForces[1] * wFxdz0 + mForces[3] * wFydz0 +
                          mForces[4] * wFzdz0) *
                        jaci;
              dsdp[8] = -(mForces[2] * wFxdz0 + mForces[4] * wFydz0 +
                          mForces[5] * wFzdz0) *
                        jaci;
              dsdp[9] = -wFx * jaci;
              dsdp[10] = 0;
              dsdp[11] = 0;
              dsdp[12] = -wFy * jaci;
              dsdp[13] = -wFx * jaci;
              dsdp[14] = 0;
              dsdp[15] = -wFz * jaci;
              dsdp[16] = 0;
              dsdp[17] = -wFx * jaci;
              dsdp[18] = 0;
              dsdp[19] = -wFy * jaci;
              dsdp[20] = 0;
              dsdp[21] = 0;
              dsdp[22] = -wFz * jaci;
              dsdp[23] = -wFy * jaci;
              dsdp[24] = 0;
              dsdp[25] = 0;
              dsdp[26] = -wFz * jaci;

              // Matrices needed for computing the Hessian wrt
              // (x0,y0,z0,mxx,mxy,mxz,myy,myz,mzz)
              float_sw4 dddp[9], dh1[9], dh2[9], dh3[9];
              dddp[0] = -wFxdx0 * jaci;
              dddp[1] = -wFxdy0 * jaci;
              dddp[2] = -wFxdz0 * jaci;
              dddp[3] = -wFydx0 * jaci;
              dddp[4] = -wFydy0 * jaci;
              dddp[5] = -wFydz0 * jaci;
              dddp[6] = -wFzdx0 * jaci;
              dddp[7] = -wFzdy0 * jaci;
              dddp[8] = -wFzdz0 * jaci;

              //                     if( i == ic && j == jc && k == kc )
              //		     {
              //                        cout.precision(16);
              //                        cout << "forcing " << endl;
              //			cout << fx << " " << fy << " " << fz <<
              // endl; 			cout << "gradient dsdp- = " << endl;
              // for( int dd=0 ; dd <
              // 9 ; dd++ ) 			   cout << dsdp[dd] << endl;
              //                        cout << "wFz and dwFzdx0 at " << ic << "
              //                        "
              //                        << jc << " " << kc << endl;
              //			cout << wFz << endl;
              //			cout << wFzdx0 << endl;
              //			cout << wFzdy0 << endl;
              //			cout << wFzdz0 << endl;
              //		     }
              // derivative of (dsdp[0],dsdp[3],dsdp[6]) (first component)
              dh1[0] = -(mForces[0] * wFxdx0dx0 + mForces[1] * wFydx0dx0 +
                         mForces[2] * wFzdx0dx0) *
                       jaci;
              dh1[1] = -(mForces[0] * wFxdx0dy0 + mForces[1] * wFydx0dy0 +
                         mForces[2] * wFzdx0dy0) *
                       jaci;
              dh1[2] = -(mForces[0] * wFxdx0dz0 + mForces[1] * wFydx0dz0 +
                         mForces[2] * wFzdx0dz0) *
                       jaci;

              dh1[3] = dh1[1];
              dh1[4] = -(mForces[0] * wFxdy0dy0 + mForces[1] * wFydy0dy0 +
                         mForces[2] * wFzdy0dy0) *
                       jaci;
              dh1[5] = -(mForces[0] * wFxdy0dz0 + mForces[1] * wFydy0dz0 +
                         mForces[2] * wFzdy0dz0) *
                       jaci;

              dh1[6] = dh1[2];
              dh1[7] = dh1[5];
              dh1[8] = -(mForces[0] * wFxdz0dz0 + mForces[1] * wFydz0dz0 +
                         mForces[2] * wFzdz0dz0) *
                       jaci;

              // derivative of (dsdp[1],dsdp[4],dsdp[7]) (second component)
              dh2[0] = -(mForces[1] * wFxdx0dx0 + mForces[3] * wFydx0dx0 +
                         mForces[4] * wFzdx0dx0) *
                       jaci;
              dh2[1] = -(mForces[1] * wFxdx0dy0 + mForces[3] * wFydx0dy0 +
                         mForces[4] * wFzdx0dy0) *
                       jaci;
              dh2[2] = -(mForces[1] * wFxdx0dz0 + mForces[3] * wFydx0dz0 +
                         mForces[4] * wFzdx0dz0) *
                       jaci;

              dh2[3] = dh2[1];
              dh2[4] = -(mForces[1] * wFxdy0dy0 + mForces[3] * wFydy0dy0 +
                         mForces[4] * wFzdy0dy0) *
                       jaci;
              dh2[5] = -(mForces[1] * wFxdy0dz0 + mForces[3] * wFydy0dz0 +
                         mForces[4] * wFzdy0dz0) *
                       jaci;

              dh2[6] = dh2[2];
              dh2[7] = dh2[5];
              dh2[8] = -(mForces[1] * wFxdz0dz0 + mForces[3] * wFydz0dz0 +
                         mForces[4] * wFzdz0dz0) *
                       jaci;

              // derivative of (dsdp[2],dsdp[5],dsdp[8]) (third component)
              dh3[0] = -(mForces[2] * wFxdx0dx0 + mForces[4] * wFydx0dx0 +
                         mForces[5] * wFzdx0dx0) *
                       jaci;
              dh3[1] = -(mForces[2] * wFxdx0dy0 + mForces[4] * wFydx0dy0 +
                         mForces[5] * wFzdx0dy0) *
                       jaci;
              dh3[2] = -(mForces[2] * wFxdx0dz0 + mForces[4] * wFydx0dz0 +
                         mForces[5] * wFzdx0dz0) *
                       jaci;

              dh3[3] = dh3[1];
              dh3[4] = -(mForces[2] * wFxdy0dy0 + mForces[4] * wFydy0dy0 +
                         mForces[5] * wFzdy0dy0) *
                       jaci;
              dh3[5] = -(mForces[2] * wFxdy0dz0 + mForces[4] * wFydy0dz0 +
                         mForces[5] * wFzdy0dz0) *
                       jaci;

              dh3[6] = dh3[2];
              dh3[7] = dh3[5];
              dh3[8] = -(mForces[2] * wFxdz0dz0 + mForces[4] * wFydz0dz0 +
                         mForces[5] * wFzdz0dz0) *
                       jaci;

              //                  if( i==42 && j==55 && k==39 )
              //		  {
              //		     cout.precision(16);
              //                  cout <<
              //                  "-----------------------------------------------------------------------\n";
              //		  cout << "     " << i <<  " " << j << " " << k
              //<< endl;
              //                  cout << "dsp = " << dsdp[2] << " " << dsdp[5]
              //                  << "  " << dsdp[8] << endl; cout << "dh  = "
              //                  << dh3[0] << " " << dh3[1] << "  " << dh3[2]
              //                  << endl; cout << "      " << dh3[3] << " " <<
              //                  dh3[4] << "  " << dh3[5] << endl; cout << " "
              //                  << dh3[6] << " " << dh3[7] << "  " << dh3[8]
              //                  << endl; cout <<
              //                  "-----------------------------------------------------------------------"
              //                  << endl;
              //		  }

              //		  if( mAmp != 0 && (fx != 0 || fy != 0 || fz !=
              //0) )
              if (1 <= k && k <= Nz) {
                GridPointSource* sourcePtr = new GridPointSource(
                    mFreq, mT0, i, j, k, g, fx, fy, fz, mTimeDependence, mNcyc,
                    mPar, mNpar, mIpar, mNipar, dsdp, dddp, dh1, dh2, dh3);
                if (m_derivative >= 0)
                  sourcePtr->set_derivative(m_derivative, m_dir);
                point_sources.push_back(sourcePtr);
              }
              if (k <= 1 && ccbndry && upperbndry) {
                int Nzp = a_EW->m_global_nz[g + 1];
                int kk = Nzp - 1 + k;

                wFx = qX0[0] * dwghi[i - ic + 2] * wghj[j - jc + 2] *
                      wghk[k - kc + 2];
                wFy = wghi[i - ic + 2] * rX0[1] * dwghj[j - jc + 2] *
                      wghk[k - kc + 2];
                wFx += wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[0] *
                       dwghk[k - kc + 2];
                wFy += wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[1] *
                       dwghk[k - kc + 2];
                wFz = wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[2] *
                      dwghk[k - kc + 2];

                float_sw4 jaci;
                if (curvilineargp1)
                  jaci = 1 / a_EW->mJ[g + 1](i, j, kk);
                else
                  jaci = 1.0 / (0.125 * h * h * h);

                float_sw4 fx =
                    -(mForces[0] * wFx + mForces[1] * wFy + mForces[2] * wFz) *
                    jaci;
                float_sw4 fy =
                    -(mForces[1] * wFx + mForces[3] * wFy + mForces[4] * wFz) *
                    jaci;
                float_sw4 fz =
                    -(mForces[2] * wFx + mForces[4] * wFy + mForces[5] * wFz) *
                    jaci;

                GridPointSource* sourcePtr = new GridPointSource(
                    mFreq, mT0, i, j, kk, g + 1, fx, fy, fz, mTimeDependence,
                    mNcyc, mPar, mNpar, mIpar, mNipar, dsdp, dddp, dh1, dh2,
                    dh3);
                if (m_derivative >= 0)
                  sourcePtr->set_derivative(m_derivative, m_dir);
                point_sources.push_back(sourcePtr);
              }
              if (k >= Nz && ccbndry && lowerbndry) {
                int kk = k - Nz + 1;
                wFx = qX0[0] * dwghi[i - ic + 2] * wghj[j - jc + 2] *
                      wghk[k - kc + 2];
                wFy = wghi[i - ic + 2] * rX0[1] * dwghj[j - jc + 2] *
                      wghk[k - kc + 2];
                wFx += wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[0] *
                       dwghk[k - kc + 2];
                wFy += wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[1] *
                       dwghk[k - kc + 2];
                wFz = wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[2] *
                      dwghk[k - kc + 2];

                float_sw4 jaci;
                if (curvilineargm1)
                  jaci = 1 / a_EW->mJ[g - 1](i, j, kk);
                else
                  jaci = 1.0 / (8 * h * h * h);

                float_sw4 fx =
                    -(mForces[0] * wFx + mForces[1] * wFy + mForces[2] * wFz) *
                    jaci;
                float_sw4 fy =
                    -(mForces[1] * wFx + mForces[3] * wFy + mForces[4] * wFz) *
                    jaci;
                float_sw4 fz =
                    -(mForces[2] * wFx + mForces[4] * wFy + mForces[5] * wFz) *
                    jaci;

                GridPointSource* sourcePtr = new GridPointSource(
                    mFreq, mT0, i, j, kk, g - 1, fx, fy, fz, mTimeDependence,
                    mNcyc, mPar, mNpar, mIpar, mNipar, dsdp, dddp, dh1, dh2,
                    dh3);
                if (m_derivative >= 0)
                  sourcePtr->set_derivative(m_derivative, m_dir);
                point_sources.push_back(sourcePtr);
              }
            }
          }  // end for k,j,i

      if (gridrefbndry) {
        // These arrays are currently undefined across the mesh refinement
        // boundary.
        // --> Source inversion can not be done if the source is located at the
        // interface.
        float_sw4 dddp[9], dh1[9], dh2[9], dh3[9], dsdp[27];
        for (int k = kc - 2; k <= kc + 3; k++) {
          if (k <= 1) {
            //               float_sw4 hi = 1.0/(0.5*h);
            //               float_sw4 jaci = 1.0/(0.125*h*h*h);
            for (int j = jcref - 2; j <= jcref + 3; j++)
              for (int i = icref - 2; i <= icref + 3; i++) {
                float_sw4 wFx = 0, wFy = 0, wFz = 0;
                if (a_EW->interior_point_in_proc(i, j, g + 1)) {
                  wFx += qX0[0] * dwghiref[i - icref + 2] *
                         wghjref[j - jcref + 2] * wghkref[k - kc + 2] * 2;
                  wFy += wghiref[i - icref + 2] * rX0[1] *
                         dwghjref[j - jcref + 2] * wghkref[k - kc + 2] * 2;
                  wFx += wghiref[i - icref + 2] * wghjref[j - jcref + 2] *
                         sX0[0] * dwghkref[k - kc + 2] * 2;
                  wFy += wghiref[i - icref + 2] * wghjref[j - jcref + 2] *
                         sX0[1] * dwghkref[k - kc + 2] * 2;
                  wFz += wghiref[i - icref + 2] * wghjref[j - jcref + 2] *
                         sX0[2] * dwghkref[k - kc + 2] * 2;

                  //                        wFx = dwghiref[i-icref+2]*
                  //                        wghjref[j-jcref+2]*
                  //                        wghkref[k-kc+2]*hi; wFy =
                  //                        wghiref[i-icref+2]*dwghjref[j-jcref+2]*
                  //                        wghkref[k-kc+2]*hi; wFz =
                  //                        wghiref[i-icref+2]*
                  //                        wghjref[j-jcref+2]*dwghkref[k-kc+2]*hi;

                  int Nzp = a_EW->m_global_nz[g + 1];
                  int kk = Nzp - 1 + k;
                  float_sw4 jaci;
                  if (curvilineargp1)
                    jaci = 1 / a_EW->mJ[g + 1](i, j, kk);
                  else
                    jaci = 1 / (0.125 * h * h * h);

                  float_sw4 fx = -(mForces[0] * wFx + mForces[1] * wFy +
                                   mForces[2] * wFz) *
                                 jaci;
                  float_sw4 fy = -(mForces[1] * wFx + mForces[3] * wFy +
                                   mForces[4] * wFz) *
                                 jaci;
                  float_sw4 fz = -(mForces[2] * wFx + mForces[4] * wFy +
                                   mForces[5] * wFz) *
                                 jaci;
                  GridPointSource* sourcePtr = new GridPointSource(
                      mFreq, mT0, i, j, kk, g + 1, fx, fy, fz, mTimeDependence,
                      mNcyc, mPar, mNpar, mIpar, mNipar, dsdp, dddp, dh1, dh2,
                      dh3);
                  point_sources.push_back(sourcePtr);
                }
              }
          }  // end if k <= 1

          if (k >= Nz) {
            //               float_sw4 jaci = 1.0/(8*h*h*h);
            //               float_sw4 hi = 1.0/(2*h);
            for (int j = jcref - 2; j <= jcref + 3; j++)
              for (int i = icref - 2; i <= icref + 3; i++) {
                float_sw4 wFx = 0, wFy = 0, wFz = 0;
                if (a_EW->interior_point_in_proc(i, j, g - 1)) {
                  wFx += qX0[0] * dwghiref[i - icref + 2] *
                         wghjref[j - jcref + 2] * wghkref[k - kc + 2] * 0.5;
                  wFy += wghiref[i - icref + 2] * rX0[1] *
                         dwghjref[j - jcref + 2] * wghkref[k - kc + 2] * 0.5;
                  wFx += wghiref[i - icref + 2] * wghjref[j - jcref + 2] *
                         sX0[0] * dwghkref[k - kc + 2];
                  wFy += wghiref[i - icref + 2] * wghjref[j - jcref + 2] *
                         sX0[1] * dwghkref[k - kc + 2];
                  wFz += wghiref[i - icref + 2] * wghjref[j - jcref + 2] *
                         sX0[2] * dwghkref[k - kc + 2];
                  //                        wFx = dwghiref[i-icref+2]*
                  //                        wghjref[j-jcref+2]*
                  //                        wghkref[k-kc+2]*hi; wFy =
                  //                        wghiref[i-icref+2]*dwghjref[j-jcref+2]*
                  //                        wghkref[k-kc+2]*hi; wFz =
                  //                        wghiref[i-icref+2]*
                  //                        wghjref[j-jcref+2]*dwghkref[k-kc+2]*2*hi;

                  int kk = k - Nz + 1;
                  float_sw4 jaci;
                  if (curvilineargm1)
                    jaci = 1.0 / a_EW->mJ[g - 1](i, j, kk);
                  else
                    jaci = 1.0 / (8 * h * h * h);
                  float_sw4 fx = -(mForces[0] * wFx + mForces[1] * wFy +
                                   mForces[2] * wFz) *
                                 jaci;
                  float_sw4 fy = -(mForces[1] * wFx + mForces[3] * wFy +
                                   mForces[4] * wFz) *
                                 jaci;
                  float_sw4 fz = -(mForces[2] * wFx + mForces[4] * wFy +
                                   mForces[5] * wFz) *
                                 jaci;

                  GridPointSource* sourcePtr = new GridPointSource(
                      mFreq, mT0, i, j, kk, g - 1, fx, fy, fz, mTimeDependence,
                      mNcyc, mPar, mNpar, mIpar, mNipar, dsdp, dddp, dh1, dh2,
                      dh3);
                  point_sources.push_back(sourcePtr);
                }
              }  // end for i

          }  // end if kz >= Nz

        }  // end for kc

      }  // end if gridrefbndry
    }
  }  // end if momentTensorSource

}  // end set_grid_point_sources4

//-----------------------------------------------------------------------
void Source::exact_testmoments(int kx[3], int ky[3], int kz[3],
                               float_sw4 momex[3]) {
  // Integrals over the domain of a polynomial of degree (kx,ky,kz) times the
  // source
  if (!mIsMomentSource) {
    float_sw4 x1, y1, z1;
    for (int c = 0; c < 3; c++) {
      if (kx[c] == 0)
        x1 = 1;
      else
        x1 = pow(mX0, kx[c]);
      if (ky[c] == 0)
        y1 = 1;
      else
        y1 = pow(mY0, ky[c]);
      if (kz[c] == 0)
        z1 = 1;
      else
        z1 = pow(mZ0, kz[c]);
      momex[c] = mForces[c] * x1 * y1 * z1;
    }
  } else {
    float_sw4 x1, y1, z1, xp1, yp1, zp1;
    for (int c = 0; c < 3; c++) {
      if (kx[c] == 0)
        x1 = 1;
      else
        x1 = pow(mX0, kx[c]);
      if (kx[c] == 0)
        xp1 = 0;
      else if (kx[c] == 1)
        xp1 = -1;
      else
        xp1 = -kx[c] * pow(mX0, (kx[c] - 1));

      if (ky[c] == 0)
        y1 = 1;
      else
        y1 = pow(mY0, ky[c]);
      if (ky[c] == 0)
        yp1 = 0;
      else if (ky[c] == 1)
        yp1 = -1;
      else
        yp1 = -ky[c] * pow(mY0, (ky[c] - 1));

      if (kz[c] == 0)
        z1 = 1;
      else
        z1 = pow(mZ0, kz[c]);
      if (kz[c] == 0)
        zp1 = 0;
      else if (kz[c] == 1)
        zp1 = -1;
      else
        zp1 = -kz[c] * pow(mZ0, (kz[c] - 1));
      if (c == 0)
        momex[c] = -(mForces[0] * xp1 * y1 * z1 + mForces[1] * x1 * yp1 * z1 +
                     mForces[2] * x1 * y1 * zp1);
      else if (c == 1)
        momex[c] = -(mForces[1] * xp1 * y1 * z1 + mForces[3] * x1 * yp1 * z1 +
                     mForces[4] * x1 * y1 * zp1);
      else
        momex[c] = -(mForces[2] * xp1 * y1 * z1 + mForces[4] * x1 * yp1 * z1 +
                     mForces[5] * x1 * y1 * zp1);
    }
  }
}

//-----------------------------------------------------------------------
void Source::perturb(float_sw4 h, int comp) {
  if (comp == 0)
    mX0 += h;
  else if (comp == 1)
    mY0 += h;
  else if (comp == 2)
    mZ0 += h;
  else if (comp >= 3 && comp <= 8)
    mForces[comp - 3] += h;
  else if (comp == 9)
    mT0 += h;
  else
    mFreq += h;
}

//-----------------------------------------------------------------------
void Source::filter_timefunc(Filter* filter_ptr, float_sw4 tstart, float_sw4 dt,
                             int nsteps) {
  if (!m_is_filtered) {
    float_sw4 (*timeFunc)(float_sw4 f, float_sw4 t, float_sw4 * par, int npar,
                          int* ipar, int nipar);
    switch (mTimeDependence) {
      case iRicker:
        timeFunc = RickerWavelet;
        break;
      case iGaussian:
        timeFunc = Gaussian;
        break;
      case iRamp:
        timeFunc = Ramp;
        break;
      case iTriangle:
        timeFunc = Triangle;
        break;
      case iSawtooth:
        timeFunc = Sawtooth;
        break;
      case iSmoothWave:
        timeFunc = SmoothWave;
        break;
      case iErf:
        timeFunc = Erf;
        break;
      case iVerySmoothBump:
        timeFunc = VerySmoothBump;
        break;
      case iC6SmoothBump:
        timeFunc = C6SmoothBump;
        break;
      case iRickerInt:
        timeFunc = RickerInt;
        break;
      case iBrune:
        timeFunc = Brune;
        break;
      case iBruneSmoothed:
        timeFunc = BruneSmoothed;
        break;
      case iDBrune:
        timeFunc = DBrune;
        break;
      case iGaussianWindow:
        timeFunc = GaussianWindow;
        break;
      case iLiu:
        timeFunc = Liu;
        break;
      case iDirac:
        timeFunc = Dirac;
        break;
        // discrete time functions are now handled differently
        // case iDiscrete :
        //    timeFunc = Discrete;
        //    break;
      default:
        cout << "ERROR in Source::filter_timefunc, source type not recoginzed"
             << endl;
    }

    // Convert to discrete representation
    mTimeDependence = iDiscrete;

    float_sw4* discfunc = new float_sw4[nsteps];
    for (int k = 0; k < nsteps; k++)
      discfunc[k] =
          timeFunc(mFreq, tstart + k * dt - mT0, mPar, mNpar, mIpar, mNipar);

    // Filter the discretized function
    filter_ptr->evaluate(nsteps, &discfunc[0], &discfunc[0]);

    // Give the source time function a smooth start if this is a 2-pass (forward
    // + backward) bandpass filter
    if (filter_ptr->get_passes() == 2 && filter_ptr->get_type() == bandPass) {
      float_sw4 wghv, xi;
      int p0 = 3,
          p = 20;  // First non-zero time level, and number of points in ramp;
      if (p0 + p <= nsteps) {
        for (int i = 1; i <= p0 - 1; i++) {
          discfunc[i - 1] = 0;
        }
        for (int i = p0; i <= p0 + p; i++) {
          wghv = 0;
          xi = (i - p0) / ((float_sw4)p);
          // polynomial P(xi), P(0) = 0, P(1)=1
          wghv = xi * xi * xi * xi *
                 (35 - 84 * xi + 70 * xi * xi - 20 * xi * xi * xi);
          discfunc[i - 1] *= wghv;
        }
      }
    }

    // Save discrete function
    mNipar = 1;
    mIpar = SW4_NEW(Space::Managed, int[mNipar]);
    mIpar[0] = nsteps;

    mFreq = 1. / dt;
    ::operator delete[](mPar, Space::Managed);
    mNpar = nsteps + 1;
    mPar = SW4_NEW(Space::Managed, float_sw4[mNpar]);
    mPar[0] = tstart;  // regular (like Gaussian) time functions are defined
                       // from t=tstart=0
    mT0 = tstart;
    //      mPar[0] = tstart-mT0;
    for (int i = 0; i < nsteps; i++) mPar[i + 1] = discfunc[i];
    delete[] discfunc;

    // Build the spline representation
    spline_interpolation();
    m_is_filtered = true;
  }
}

//-----------------------------------------------------------------------
int Source::spline_interpolation() {
  // Assume mPar[1], to mPar[npts] contain the function
  // Assume mIpar[0] contains npts
  // Assume mFreq contains 1/dt, and mPar[0] is tstart.
  // Compute the six spline coefficients for each interval and return in
  // mPar[1],to mPar[6*(npts-1)]
  if (mTimeDependence == iDiscrete) {
    int npts = mIpar[0];

    // tmp
    // int myRank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // if (myRank == 0)
    // {
    // 	cout << "before spline interp" << endl;
    // 	cout << "npts = " << npts << " t0 = " << mPar[0] << " dt= " << 1/mFreq
    // << endl; 	for( int i=0 ; i < npts ; i++ ) 	  cout << "fun["
    // <<
    // i
    // <<
    // "]
    // =
    // "<< mPar[i+1] << endl;
    // }

    Qspline quinticspline(npts, &mPar[1], mPar[0], 1 / mFreq);
    float_sw4 tstart = mPar[0];
    ::operator delete[](mPar, Space::Managed);
    mPar = SW4_NEW(Space::Managed, float_sw4[6 * (npts - 1) + 1]);
    mNpar = 6 * (npts - 1) + 1;
    mPar[0] = tstart;
    float_sw4* qsppt = quinticspline.get_polycof_ptr();
    for (int i = 0; i < 6 * (npts - 1); i++) mPar[i + 1] = qsppt[i];
    //      cout << "after spline interp" << endl;
    //      for( int i=0 ; i < npts ; i++ )
    //	 cout << "fun[" << i << "] = "<< mPar[6*i+1] << endl;

    return 1;
  } else if (mTimeDependence == iDiscrete3forces) {
    int npts = mIpar[0];
    // Three different time function, one for each displacement component.
    // Store sequentially in mPar, i.e.,
    // mPar = [tstart, first time func, tstart, second time func, ...]
    // tstart before each function ---> Can reuse time function iDiscrete.

    float_sw4* parin = new float_sw4[(npts + 1) * 3];
    for (int i = 0; i < (npts + 1) * 3; i++) parin[i] = mPar[i];
    float_sw4 tstart = mPar[0];
    delete[] mPar;
    mNpar = 3 * (6 * (npts - 1) + 1);
    mPar = new float_sw4[mNpar];

    size_t pos_in = 0, pos_out = 0;
    for (int tf = 0; tf < 3; tf++) {
      Qspline quinticspline(npts, &parin[pos_in + 1], tstart, 1 / mFreq);
      pos_in += npts + 1;
      mPar[pos_out] = tstart;
      float_sw4* qsppt = quinticspline.get_polycof_ptr();
      for (int i = 0; i < 6 * (npts - 1); i++) mPar[pos_out + i + 1] = qsppt[i];
      pos_out += 6 * (npts - 1) + 1;
    }
    delete[] parin;
    return 1;
  } else if (mTimeDependence == iDiscrete6moments) {
    int npts = mIpar[0];
    // Six different time function, one for each momentum component.
    // Store sequentially in mPar, i.e.,
    // mPar = [tstart, first time func, tstart, second time func, ...]
    // tstart before each function ---> Can reuse time function iDiscrete.

    float_sw4* parin = new float_sw4[(npts + 1) * 6];
    for (int i = 0; i < (npts + 1) * 6; i++) parin[i] = mPar[i];
    float_sw4 tstart = mPar[0];
    ::operator delete[](mPar, Space::Managed);
    mNpar = 6 * (6 * (npts - 1) + 1);
    mPar = SW4_NEW(Space::Managed, float_sw4[mNpar]);

    size_t pos_in = 0, pos_out = 0;
    for (int tf = 0; tf < 6; tf++) {
      Qspline quinticspline(npts, &parin[pos_in + 1], tstart, 1 / mFreq);
      pos_in += npts + 1;
      mPar[pos_out] = tstart;
      float_sw4* qsppt = quinticspline.get_polycof_ptr();
      for (int i = 0; i < 6 * (npts - 1); i++) mPar[pos_out + i + 1] = qsppt[i];
      pos_out += 6 * (npts - 1) + 1;
    }
    delete[] parin;
    return 1;
  } else
    return 0;
}

//-----------------------------------------------------------------------
Source* Source::copy(std::string a_name) {
  if (a_name == " ") a_name = mName;

  Source* retval = new Source();
  //   retval->m_i0 = m_i0;
  //   retval->m_j0 = m_j0;
  //   retval->m_k0 = m_k0;
  //   retval->m_grid = m_grid;
  retval->mName = a_name;
  retval->mIsMomentSource = mIsMomentSource;
  retval->mForces.push_back(mForces[0]);
  retval->mForces.push_back(mForces[1]);
  retval->mForces.push_back(mForces[2]);
  if (mIsMomentSource) {
    retval->mForces.push_back(mForces[3]);
    retval->mForces.push_back(mForces[4]);
    retval->mForces.push_back(mForces[5]);
  }
  retval->mFreq = mFreq;
  retval->mT0 = mT0;
  //   retval->mGridPointSet = mGridPointSet;
  retval->m_myPoint = m_myPoint;
  retval->m_zRelativeToTopography = m_zRelativeToTopography;
  retval->mX0 = mX0;
  retval->mY0 = mY0;
  retval->mZ0 = mZ0;

  retval->mNpar = mNpar;
  retval->mPar = SW4_NEW(Space::Managed, float_sw4[mNpar]);
  for (int i = 0; i < mNpar; i++) retval->mPar[i] = mPar[i];

  retval->mNipar = mNipar;
  retval->mIpar = new int[mNipar];
  for (int i = 0; i < mNipar; i++) retval->mIpar[i] = mIpar[i];

  retval->mNcyc = mNcyc;
  retval->m_derivative = m_derivative;
  retval->mTimeDependence = mTimeDependence;
  for (int i = 0; i < 11; i++) retval->m_dir[i] = m_dir[i];
  retval->m_is_filtered = m_is_filtered;

  retval->m_zTopo = m_zTopo;
  retval->mIgnore = mIgnore;
  retval->mShearModulusFactor = mShearModulusFactor;

  return retval;
}

//-----------------------------------------------------------------------
float_sw4 Source::find_min_exponent() const {
  // smallest number x, such that exp(x) does not cause underflow
  return -700.0;
}

//-----------------------------------------------------------------------
void Source::compute_metric_at_source(EW* a_EW, float_sw4 q, float_sw4 r,
                                      float_sw4 s, int ic, int jc, int kc,
                                      int g, float_sw4& zq, float_sw4& zr,
                                      float_sw4& zs, float_sw4& zqq,
                                      float_sw4& zqr, float_sw4& zqs,
                                      float_sw4& zrr, float_sw4& zrs,
                                      float_sw4& zss) const {
  int Nz = a_EW->m_global_nz[g];
  float_sw4 h = a_EW->mGridSize[g];
  if (kc > Nz - 3) {
    // Treat downmost grid lines as cartesian
    zq = zr = 0;
    zs = h;
    zqq = zqr = zqs = zrr = zrs = zss = 0;
  } else {
    bool analytic_derivative = true;
    // Derivative of metric wrt. source position. Not yet fully implemented.
    //      double zqdx0, zqdy0, zqsz0, zrdx0, zrdy0, zrdz0;

    // 3. Recompute metric to sixth order accuracy. Increased accuracy needed
    // because
    //    of the multiplication with a singular (Dirac) function.
    //    compute only in the processor where the point is interior

    //      bool eightptstencil = true;
    float_sw4 ai = q - ic;
    float_sw4 bi = r - jc;
    float_sw4 ci = s - kc;

    //      double d6cofi[8], d6cofj[8], d6cofk[8];
    //      double dd6cofi[8], dd6cofj[8], dd6cofk[8];
    //      if( eightptstencil )
    //      {
    //	 // Eight point stencil, smooth wrt. source position,
    //	 // needed for source optimization
    //	 getmetdwgh( ai, d6cofi );
    //	 getmetdwgh( bi, d6cofj );
    //	 if( kc <= 3 && ci < 0 )
    //	 {
    //	    getmetdwgh7( ci, d6cofk );
    ///	    d6cofk[7] = 0;
    //	 }
    //	 else
    //	    getmetdwgh( ci, d6cofk );
    //      }
    //      else
    //      {
    //	 // Seven point stencil, ok for forward solver.
    //	 getmetdwgh7( ai, d6cofi );
    //	 d6cofi[7] = 0;
    //	 getmetdwgh7( bi, d6cofj );
    //	 d6cofj[7] = 0;
    //	 getmetdwgh7( ci, d6cofk );
    //	 d6cofk[7] = 0;
    //      }

    float_sw4 a6cofi[8], a6cofj[8], a6cofk[8];
    float_sw4 d6cofi[8], d6cofj[8], d6cofk[8];
    float_sw4 dd6cofi[8], dd6cofj[8], dd6cofk[8];
    float_sw4 ddd6cofi[8], ddd6cofj[8], ddd6cofk[8];
    //      if( eightptstencil )
    //      {
    getmetwgh(ai, a6cofi, d6cofi, dd6cofi, ddd6cofi);
    getmetwgh(bi, a6cofj, d6cofj, dd6cofj, ddd6cofj);
    getmetwgh(ci, a6cofk, d6cofk, dd6cofk, ddd6cofk);
    //	 if( kc <= 3 && ci < 0 )
    //	 {
    //	    getmetwgh7( ci, a6cofk );
    //	    a6cofk[7] = 0;
    //	 }
    //	 else
    //	    getmetwgh( ci, a6cofk );
    //      }
    //      else
    //      {
    //	 getmetwgh7( ai, a6cofi );
    //	 a6cofi[7] = 0;
    //	 getmetwgh7( bi, a6cofj );
    //	 a6cofj[7] = 0;
    //	 getmetwgh7( ci, a6cofk );
    //	 a6cofk[7] = 0;
    //      }

    // Assume grid uniform in x and y, compute metric with z=z(q,r,s)
    zq = zr = zs = 0;
    int order;
    float_sw4 zetaBreak;
    //    a_EW->get_gridgen_info(order, zetaBreak);
    a_EW->m_gridGenerator->get_gridgen_info(order, zetaBreak);
    float_sw4 zpar = (s - 1) / (zetaBreak * (Nz - 1));
    float_sw4 kBreak = 1 + zetaBreak * (Nz - 1);

    if (zpar >= 1) {
      zq = 0;
      zr = 0;
      zs = h;
      zqq = zqr = zqs = zrr = zrs = zss = 0;
    } else {
      float_sw4 pp = pow(1 - zpar, order - 1);
      float_sw4 powo = (1 - zpar) * pp;
      float_sw4 dpowo = -order * pp / zetaBreak;
      float_sw4 tauavg = 0;
      float_sw4 tauq = 0, taur = 0;
      float_sw4 tauqq = 0, tauqr = 0, taurr = 0;
      for (int j = jc - 3; j <= jc + 4; j++)
        for (int i = ic - 3; i <= ic + 4; i++) {
          tauavg += a6cofi[i - (ic - 3)] * a6cofj[j - (jc - 3)] *
                    a_EW->mTopoGridExt(i, j, 1);
          tauq += d6cofi[i - (ic - 3)] * a6cofj[j - (jc - 3)] *
                  a_EW->mTopoGridExt(i, j, 1);
          taur += a6cofi[i - (ic - 3)] * d6cofj[j - (jc - 3)] *
                  a_EW->mTopoGridExt(i, j, 1);
          tauqq += dd6cofi[i - (ic - 3)] * a6cofj[j - (jc - 3)] *
                   a_EW->mTopoGridExt(i, j, 1);
          tauqr += d6cofi[i - (ic - 3)] * d6cofj[j - (jc - 3)] *
                   a_EW->mTopoGridExt(i, j, 1);
          taurr += a6cofi[i - (ic - 3)] * dd6cofj[j - (jc - 3)] *
                   a_EW->mTopoGridExt(i, j, 1);
        }
      //	 double powo = pow(1-zpar,order);
      //      double dpowo = -order*pow(1-zpar,order-1);
      zq = (-tauq) * powo;
      zr = (-taur) * powo;
      zqq = (-tauqq) * powo;
      zqr = (-tauqr) * powo;
      zrr = (-taurr) * powo;
      zqs = (-tauq) * dpowo / (Nz - 1);
      zrs = (-taur) * dpowo / (Nz - 1);

      // Compute dz/ds directly from the grid mapping, the explicit expression
      // here should be the same as in EW::curvilinear_grid_mapping
      float_sw4 zMax =
          a_EW->m_zmin[a_EW->mNumberOfCartesianGrids - 1] - (Nz - kBreak) * h;
      float_sw4 c1 = zMax + tauavg - h * (kBreak - 1);

      // Divide by Nz-1 to make consistent with undivided differences
      if (analytic_derivative) {
        zs = h + c1 * (-dpowo) / (Nz - 1);
        zss = -c1 * order * (order - 1) * pow(1 - zpar, order - 2) /
              (zetaBreak * zetaBreak * (Nz - 1) * (Nz - 1));
        //         cout << "AN: zs = " << zs << " zss= " << zss << endl;
      } else {
        zs = 0;
        zss = 0;
        float_sw4 z1d = 0;
        for (int k = kc - 3; k <= kc + 4; k++) {
          zpar = (k - 1) / (zetaBreak * (Nz - 1));
          if (zpar >= 1)
            z1d = zMax + (k - kBreak) * h;
          else {
            z1d = (1 - zpar) * (-tauavg) + zpar * (zMax + c1 * (1 - zpar));
            for (int o = 2; o < order; o++) z1d += zpar * c1 * pow(1 - zpar, o);
          }
          zs += d6cofk[k - (kc - 3)] * z1d;
          zss += dd6cofk[k - (kc - 3)] * z1d;
        }
        //         cout << "NU: zs = " << zs << " zss= " << zss << endl;
      }
    }
  }
}

extern "C" {
void F77_FUNC(dgels, DGELS)(char&, int&, int&, int&, double*, int&, double*,
                            int&, double*, int&, int&);
}

//-----------------------------------------------------------------------
void Source::get_mr_psources(EW* a_EW, int g, float_sw4 q, float_sw4 r,
                             float_sw4 s, bool gradient, float_sw4 normwgh[4],
                             vector<GridPointSource*>& point_sources) {
#define CUB(x) (x) * (x) * (x)
#define BISQR(x) (x) * (x) * (x) * (x)
  //   ic, jc on current grid
  //   icref, jcref on adjacent grid
  // kc
  // Stencil is always -2 <= m <= 2, -2 <= n <= 2, -2 <= o <= 2
  //   std::cout <<a_EW->getRank() << " in get_mr_psources " << std::endl;
  int ncond = 35;  // Number of moment conditions
#define a(i, j) a_[(i - 1) + ncond * (j - 1)]
#define b(i, j) b_[(i - 1) + 125 * (j - 1)]
#define x(i, j) b_[(i - 1) + 125 * (j - 1)]
  // setup stencil centers etc.
  int ic = static_cast<int>(round(q));
  int jc = static_cast<int>(round(r));
  int kc = static_cast<int>(round(s));
  int Nz = a_EW->m_global_nz[g];
  int icref, jcref, gref;
  if (kc > Nz / 2) {
    // (q,r,s) in fine grid, stencil extends to coarse grid below
    icref = static_cast<int>(round(0.5 * (q + 1)));
    jcref = static_cast<int>(round(0.5 * (r + 1)));
    gref = g - 1;
  } else {
    // (q,r,s) in coarse grid, stencil extends to fine grid above
    icref = static_cast<int>(round(2 * q - 1));
    jcref = static_cast<int>(round(2 * r - 1));
    gref = g + 1;
  }
  float_sw4 h = a_EW->mGridSize[g], href = a_EW->mGridSize[gref];

  // Processor owns part of source stencil if (ic,jc) is inside processor.
  // If (ic,jc) is in the processor overlap, then all the stencil widths
  // ic-2,..,ic+2 and jc-2,..,jc+2 are not guaranteed to be owned by this
  // processor. Need to supply grid z and Jacobian for all stencil points,
  // extended arrays are being used.
  if (a_EW->interior_point_in_proc(ic - 2, jc - 2, g) ||
      a_EW->interior_point_in_proc(ic + 2, jc - 2, g) ||
      a_EW->interior_point_in_proc(ic - 2, jc + 2, g) ||
      a_EW->interior_point_in_proc(ic + 2, jc + 2, g) ||
      a_EW->interior_point_in_proc(icref - 2, jcref - 2, gref) ||
      a_EW->interior_point_in_proc(icref + 2, jcref - 2, gref) ||
      a_EW->interior_point_in_proc(icref - 2, jcref + 2, gref) ||
      a_EW->interior_point_in_proc(icref + 2, jcref + 2, gref)) {
    Sarray zg(ic - 2, ic + 2, jc - 2, jc + 2, kc - 2, kc + 2);
    Sarray Jg(zg);
    if (!a_EW->interior_point_in_proc(ic, jc, g))
      a_EW->m_gridGenerator->generate_z_and_j(a_EW, g, zg, Jg);
    else {
      zg.insert_intersection(a_EW->mZ[g]);
      Jg.insert_intersection(a_EW->mJ[g]);
    }
    int kll, kul;
    if (gref == g + 1) {
      int Nzp = a_EW->m_global_nz[g + 1];
      kll = Nzp - 4;
      kul = Nzp;
    } else {
      kll = 1;
      kul = 5;
    }
    Sarray zgref(icref - 2, icref + 2, jcref - 2, jcref + 2, kll, kul);
    Sarray Jref(zgref);
    if (!a_EW->interior_point_in_proc(icref, jcref, gref))
      a_EW->m_gridGenerator->generate_z_and_j(a_EW, gref, zgref, Jref);
    else {
      zgref.insert_intersection(a_EW->mZ[gref]);
      Jref.insert_intersection(a_EW->mJ[gref]);
    }

    int nrhs = 1;
    if (gradient) nrhs = 4;

    float_sw4* a_ = new float_sw4[ncond * 125];
    int ldb = ncond > 125 ? ncond : 125;
    float_sw4* b_ = new float_sw4[ldb * nrhs];
    //   float_sw4* x_ = new float_sw4[125];

    //   std::cout << "SOURCE at interface, g= " << g << " (ic,jc,kc)= " << ic
    //   <<
    //      ", " << jc << ", " << kc << std::endl;
    float_sw4 xc = h * (ic - 1);
    float_sw4 yc = h * (jc - 1);
    float_sw4 zc = zg(ic, jc, kc);

    for (int k = kc - 2; k <= kc + 2; k++) {
      if (1 <= k && k <= Nz - 1) {
        // Discretization on this grid
        for (int j = jc - 2; j <= jc + 2; j++)
          for (int i = ic - 2; i <= ic + 2; i++) {
            int ind = i - ic + 2 + 5 * (j - jc + 2) + 25 * (k - kc + 2) + 1;
            float_sw4 x = (i - 1) * h, y = (j - 1) * h, z = zg(i, j, k);
            a(1, ind) = 1;
            a(2, ind) = x - xc;
            a(3, ind) = y - yc;
            a(4, ind) = z - zc;
            a(5, ind) = SQR(x - xc);
            a(6, ind) = SQR(y - yc);
            a(7, ind) = SQR(z - zc);
            a(8, ind) = (x - xc) * (y - yc);
            a(9, ind) = (x - xc) * (z - zc);
            a(10, ind) = (y - yc) * (z - zc);
            a(11, ind) = CUB(x - xc);
            a(12, ind) = CUB(y - yc);
            a(13, ind) = CUB(z - zc);
            a(14, ind) = SQR(x - xc) * (y - yc);
            a(15, ind) = SQR(y - yc) * (x - xc);
            a(16, ind) = SQR(z - zc) * (x - xc);
            a(17, ind) = SQR(x - xc) * (z - zc);
            a(18, ind) = SQR(y - yc) * (z - zc);
            a(19, ind) = SQR(z - zc) * (y - yc);
            a(20, ind) = (x - xc) * (y - yc) * (z - zc);
            if (ncond == 35) {
              a(21, ind) = BISQR(x - xc);
              a(22, ind) = BISQR(y - yc);
              a(23, ind) = BISQR(z - zc);
              a(24, ind) = CUB(x - xc) * (y - yc);
              a(25, ind) = CUB(y - yc) * (x - xc);
              a(26, ind) = CUB(z - zc) * (x - xc);
              a(27, ind) = CUB(x - xc) * (z - zc);
              a(28, ind) = CUB(z - zc) * (y - yc);
              a(29, ind) = CUB(y - yc) * (z - zc);
              a(30, ind) = SQR(x - xc) * SQR(y - yc);
              a(31, ind) = SQR(y - yc) * SQR(z - zc);
              a(32, ind) = SQR(z - zc) * SQR(x - xc);
              a(33, ind) = SQR(x - xc) * (y - yc) * (z - zc);
              a(34, ind) = SQR(y - yc) * (x - xc) * (z - zc);
              a(35, ind) = SQR(z - zc) * (x - xc) * (y - yc);
            }
          }
      }
      if (k <= 0 || k >= Nz) {
        int kk;
        if (k <= 0) {
          // Into finer grid above, continue into fine grid from k=0 and
          // upward --> from kk=Nz-1
          int Nzp = a_EW->m_global_nz[g + 1];
          kk = k + Nzp - 1;
        }
        if (k >= Nz) {
          // Into coarser grid below, let k=Nz be the coarse grid
          // discretization, and continue on coarse grid from k=Nz+1 --> kk=2
          kk = k - Nz + 1;
        }
        for (int j = jcref - 2; j <= jcref + 2; j++)
          for (int i = icref - 2; i <= icref + 2; i++) {
            float_sw4 x = (i - 1) * href, y = (j - 1) * href,
                      z = zgref(i, j, kk);
            int ind =
                i - icref + 2 + 5 * (j - jcref + 2) + 25 * (k - kc + 2) + 1;
            a(1, ind) = 1;
            a(2, ind) = x - xc;
            a(3, ind) = y - yc;
            a(4, ind) = z - zc;
            a(5, ind) = SQR(x - xc);
            a(6, ind) = SQR(y - yc);
            a(7, ind) = SQR(z - zc);
            a(8, ind) = (x - xc) * (y - yc);
            a(9, ind) = (x - xc) * (z - zc);
            a(10, ind) = (y - yc) * (z - zc);
            a(11, ind) = CUB(x - xc);
            a(12, ind) = CUB(y - yc);
            a(13, ind) = CUB(z - zc);
            a(14, ind) = SQR(x - xc) * (y - yc);
            a(15, ind) = SQR(y - yc) * (x - xc);
            a(16, ind) = SQR(z - zc) * (x - xc);
            a(17, ind) = SQR(x - xc) * (z - zc);
            a(18, ind) = SQR(y - yc) * (z - zc);
            a(19, ind) = SQR(z - zc) * (y - yc);
            a(20, ind) = (x - xc) * (y - yc) * (z - zc);
            if (ncond == 35) {
              a(21, ind) = BISQR(x - xc);
              a(22, ind) = BISQR(y - yc);
              a(23, ind) = BISQR(z - zc);
              a(24, ind) = CUB(x - xc) * (y - yc);
              a(25, ind) = CUB(y - yc) * (x - xc);
              a(26, ind) = CUB(z - zc) * (x - xc);
              a(27, ind) = CUB(x - xc) * (z - zc);
              a(28, ind) = CUB(z - zc) * (y - yc);
              a(29, ind) = CUB(y - yc) * (z - zc);
              a(30, ind) = SQR(x - xc) * SQR(y - yc);
              a(31, ind) = SQR(y - yc) * SQR(z - zc);
              a(32, ind) = SQR(z - zc) * SQR(x - xc);
              a(33, ind) = SQR(x - xc) * (y - yc) * (z - zc);
              a(34, ind) = SQR(y - yc) * (x - xc) * (z - zc);
              a(35, ind) = SQR(z - zc) * (x - xc) * (y - yc);
            }
          }
      }
      if (k == 1 || k == Nz) {
        // On interface, impose special condition on the fine side
        int gf, icf, jcf, icc, jcc, Nzf;
        const float_sw4 i256 = 1.0 / 256;
        //         std::cout << "SOURCE:  k= " << k << " ic,jc= " << ic << " "
        //         << jc <<
        //            " icref,jcref= " << icref << " " << jcref << std::endl;
        if (k == 1) {
          // If k==1, icref,jcref are in fine grid
          gf = g + 1;
          //            icf=icref;
          //            jcf=jcref;
          icc = ic;
          jcc = jc;
          icf = 2 * ic - 1;
          jcf = 2 * jc - 1;
          Nzf = a_EW->m_global_nz[g + 1];
        } else {
          // If k==Nz, icref,jcref are in the coarse grid
          gf = g;
          //            icf=ic;
          //            jcf=jc;
          icc = icref;
          jcc = jcref;
          icf = 2 * icc - 1;
          jcf = 2 * jcc - 1;
          Nzf = Nz;
        }
        float_sw4 hf = a_EW->mGridSize[gf];
        float_sw4 h = a_EW->mGridSize[g];
        Sarray zfa(icf - 7, icf + 7, jcf - 7, jcf + 7, Nzf - 5, Nzf);
        Sarray Jfa(icf - 7, icf + 7, jcf - 7, jcf + 7, Nzf - 5, Nzf);
        a_EW->m_gridGenerator->generate_z_and_j(a_EW, gf, zfa, Jfa);
        float_sw4* _mom_f = new double[15 * 15 * ncond];
#define mom_f(c, i, j) \
  _mom_f[(c - 1) + ncond * (i - icf + 7) + 15 * ncond * (j - jcf + 7)]
        for (int j = jcf - 7; j <= jcf + 7; j++)
          for (int i = icf - 7; i <= icf + 7; i++) {
            float_sw4 dx = hf * (i - 1) - xc;
            float_sw4 dy = hf * (j - 1) - yc;
            float_sw4 dz = zfa(i, j, Nzf) - zc;
            float_sw4 Jf = Jfa(i, j, Nzf);
            mom_f(1, i, j) = 1 * Jf;
            mom_f(2, i, j) = dx * Jf;
            mom_f(3, i, j) = dy * Jf;
            mom_f(4, i, j) = dz * Jf;
            mom_f(5, i, j) = dx * dx * Jf;
            mom_f(6, i, j) = dy * dy * Jf;
            mom_f(7, i, j) = dz * dz * Jf;
            mom_f(8, i, j) = dx * dy * Jf;
            mom_f(9, i, j) = dx * dz * Jf;
            mom_f(10, i, j) = dy * dz * Jf;
            mom_f(11, i, j) = dx * dx * dx * Jf;
            mom_f(12, i, j) = dy * dy * dy * Jf;
            mom_f(13, i, j) = dz * dz * dz * Jf;
            mom_f(14, i, j) = dx * dx * dy * Jf;
            mom_f(15, i, j) = dx * dy * dy * Jf;
            mom_f(16, i, j) = dx * dz * dz * Jf;
            mom_f(17, i, j) = dx * dx * dz * Jf;
            mom_f(18, i, j) = dy * dy * dz * Jf;
            mom_f(19, i, j) = dy * dz * dz * Jf;
            mom_f(20, i, j) = dx * dy * dz * Jf;
            if (ncond == 35) {
              mom_f(21, i, j) = BISQR(dx) * Jf;
              mom_f(22, i, j) = BISQR(dy) * Jf;
              mom_f(23, i, j) = BISQR(dz) * Jf;
              mom_f(24, i, j) = CUB(dx) * (dy)*Jf;
              mom_f(25, i, j) = CUB(dy) * (dx)*Jf;
              mom_f(26, i, j) = CUB(dz) * (dx)*Jf;
              mom_f(27, i, j) = CUB(dx) * (dz)*Jf;
              mom_f(28, i, j) = CUB(dz) * (dy)*Jf;
              mom_f(29, i, j) = CUB(dy) * (dz)*Jf;
              mom_f(30, i, j) = SQR(dx) * SQR(dy) * Jf;
              mom_f(31, i, j) = SQR(dy) * SQR(dz) * Jf;
              mom_f(32, i, j) = SQR(dz) * SQR(dx) * Jf;
              mom_f(33, i, j) = SQR(dx) * (dy) * (dz)*Jf;
              mom_f(34, i, j) = SQR(dy) * (dx) * (dz)*Jf;
              mom_f(35, i, j) = SQR(dz) * (dx) * (dy)*Jf;
            }
          }

        for (int jo = jcc - 2; jo <= jcc + 2; jo++)
          for (int io = icc - 2; io <= icc + 2; io++) {
            int i = 2 * io - 1, j = 2 * jo - 1;
            for (int c = 1; c <= ncond; c++) {
              float_sw4 mom_r =
                  i256 *
                  (mom_f(c, i - 3, j - 3) - 9 * mom_f(c, i - 3, j - 1) -
                   16 * mom_f(c, i - 3, j) - 9 * mom_f(c, i - 3, j + 1) +
                   mom_f(c, i - 3, j + 3) +
                   9 * (-mom_f(c, i - 1, j - 3) + 9 * mom_f(c, i - 1, j - 1) +
                        16 * mom_f(c, i - 1, j) + 9 * mom_f(c, i - 1, j + 1) -
                        mom_f(c, i - 1, j + 3)) +
                   16 * (-mom_f(c, i, j - 3) + 9 * mom_f(c, i, j - 1) +
                         16 * mom_f(c, i, j) + 9 * mom_f(c, i, j + 1) -
                         mom_f(c, i, j + 3)) +
                   9 * (-mom_f(c, i + 1, j - 3) + 9 * mom_f(c, i + 1, j - 1) +
                        16 * mom_f(c, i + 1, j) + 9 * mom_f(c, i + 1, j + 1) -
                        mom_f(c, i + 1, j + 3)) +
                   mom_f(c, i + 3, j - 3) - 9 * mom_f(c, i + 3, j - 1) -
                   16 * mom_f(c, i + 3, j) - 9 * mom_f(c, i + 3, j + 1) +
                   mom_f(c, i + 3, j + 3));
              int ind =
                  io - icc + 2 + 5 * (jo - jcc + 2) + 25 * (k - kc + 2) + 1;
              if (g == gf)
                a(c, ind) += mom_r / Jref(io, jo, 1);
              else
                a(c, ind) += mom_r / Jg(io, jo, 1);
              //                  a(c,ind) += mom_r/a_EW->mJ[gf-1](io,jo,1);
              //                  a(c,ind) += mom_r/Jref(io,jo,1);
            }
          }

        delete[] _mom_f;
      }
    }

#undef mom_f
    // Right hand sides for moment conditions

    b(1, 1) = 1;
    b(2, 1) = mX0 - xc;
    b(3, 1) = mY0 - yc;
    b(4, 1) = mZ0 - zc;
    b(5, 1) = SQR(mX0 - xc);
    b(6, 1) = SQR(mY0 - yc);
    b(7, 1) = SQR(mZ0 - zc);
    b(8, 1) = (mX0 - xc) * (mY0 - yc);
    b(9, 1) = (mX0 - xc) * (mZ0 - zc);
    b(10, 1) = (mY0 - yc) * (mZ0 - zc);
    b(11, 1) = CUB(mX0 - xc);
    b(12, 1) = CUB(mY0 - yc);
    b(13, 1) = CUB(mZ0 - zc);
    b(14, 1) = SQR(mX0 - xc) * (mY0 - yc);
    b(15, 1) = SQR(mY0 - yc) * (mX0 - xc);
    b(16, 1) = SQR(mZ0 - zc) * (mX0 - xc);
    b(17, 1) = SQR(mX0 - xc) * (mZ0 - zc);
    b(18, 1) = SQR(mY0 - yc) * (mZ0 - zc);
    b(19, 1) = SQR(mZ0 - zc) * (mY0 - yc);
    b(20, 1) = (mX0 - xc) * (mY0 - yc) * (mZ0 - zc);
    if (ncond == 35) {
      b(21, 1) = BISQR(mX0 - xc);
      b(22, 1) = BISQR(mY0 - yc);
      b(23, 1) = BISQR(mZ0 - zc);
      b(24, 1) = CUB(mX0 - xc) * (mY0 - yc);
      b(25, 1) = CUB(mY0 - yc) * (mX0 - xc);
      b(26, 1) = CUB(mZ0 - zc) * (mX0 - xc);
      b(27, 1) = CUB(mX0 - xc) * (mZ0 - zc);
      b(28, 1) = CUB(mZ0 - zc) * (mY0 - yc);
      b(29, 1) = CUB(mY0 - yc) * (mZ0 - zc);
      b(30, 1) = SQR(mX0 - xc) * SQR(mY0 - yc);
      b(31, 1) = SQR(mY0 - yc) * SQR(mZ0 - zc);
      b(32, 1) = SQR(mZ0 - zc) * SQR(mX0 - xc);
      b(33, 1) = SQR(mX0 - xc) * (mY0 - yc) * (mZ0 - zc);
      b(34, 1) = SQR(mY0 - yc) * (mX0 - xc) * (mZ0 - zc);
      b(35, 1) = SQR(mZ0 - zc) * (mX0 - xc) * (mY0 - yc);
    }
    if (gradient) {
      // x-derivatives
      b(1, 2) = 0;

      b(2, 2) = -1;
      b(3, 2) = 0;
      b(4, 2) = 0;

      b(5, 2) = -2 * (mX0 - xc);
      b(6, 2) = 0;
      b(7, 2) = 0;
      b(8, 2) = -(mY0 - yc);
      b(9, 2) = -(mZ0 - zc);
      b(10, 2) = 0;

      b(11, 2) = -3 * SQR(mX0 - xc);
      b(12, 2) = 0;
      b(13, 2) = 0;
      b(14, 2) = -2 * (mX0 - xc) * (mY0 - yc);
      b(15, 2) = -SQR(mY0 - yc);
      b(16, 2) = -SQR(mZ0 - zc);
      b(17, 2) = -2 * (mX0 - xc) * (mZ0 - zc);
      b(18, 2) = 0;
      b(19, 2) = 0;
      b(20, 2) = -(mY0 - yc) * (mZ0 - zc);
      if (ncond == 35) {
        b(21, 2) = -4 * CUB(mX0 - xc);
        b(22, 2) = 0;
        b(23, 2) = 0;
        b(24, 2) = -3 * SQR(mX0 - xc) * (mY0 - yc);
        b(25, 2) = -CUB(mY0 - yc);
        b(26, 2) = -CUB(mZ0 - zc);
        b(27, 2) = -3 * SQR(mX0 - xc) * (mZ0 - zc);
        b(28, 2) = 0;
        b(29, 2) = 0;
        b(30, 2) = -2 * (mX0 - xc) * SQR(mY0 - yc);
        b(31, 2) = 0;
        b(32, 2) = -2 * SQR(mZ0 - zc) * (mX0 - xc);
        b(33, 2) = -2 * (mX0 - xc) * (mY0 - yc) * (mZ0 - zc);
        b(34, 2) = -SQR(mY0 - yc) * (mZ0 - zc);
        b(35, 2) = -SQR(mZ0 - zc) * (mY0 - yc);
      }

      // y-derivatives
      b(1, 3) = 0;

      b(2, 3) = 0;
      b(3, 3) = -1;
      b(4, 3) = 0;

      b(5, 3) = 0;
      b(6, 3) = -2 * (mY0 - yc);
      b(7, 3) = 0;
      b(8, 3) = -(mX0 - xc);
      b(9, 3) = 0;
      b(10, 3) = -(mZ0 - zc);

      b(11, 3) = 0;
      b(12, 3) = -3 * SQR(mY0 - yc);
      b(13, 3) = 0;
      b(14, 3) = -SQR(mX0 - xc);
      b(15, 3) = -2 * (mY0 - yc) * (mX0 - xc);
      b(16, 3) = 0;
      b(17, 3) = 0;
      b(18, 3) = -2 * (mY0 - yc) * (mZ0 - zc);
      b(19, 3) = -SQR(mZ0 - zc);
      b(20, 3) = -(mX0 - xc) * (mZ0 - zc);
      if (ncond == 35) {
        b(21, 3) = 0;
        b(22, 3) = -4 * CUB(mY0 - yc);
        b(23, 3) = 0;
        b(24, 3) = -CUB(mX0 - xc);
        b(25, 3) = -3 * SQR(mY0 - yc) * (mX0 - xc);
        b(26, 3) = 0;
        b(27, 3) = 0;
        b(28, 3) = -CUB(mZ0 - zc);
        b(29, 3) = -3 * SQR(mY0 - yc) * (mZ0 - zc);
        b(30, 3) = -2 * SQR(mX0 - xc) * (mY0 - yc);
        b(31, 3) = -2 * (mY0 - yc) * SQR(mZ0 - zc);
        b(32, 3) = 0;
        b(33, 3) = -SQR(mX0 - xc) * (mZ0 - zc);
        b(34, 3) = -2 * (mY0 - yc) * (mX0 - xc) * (mZ0 - zc);
        b(35, 3) = -SQR(mZ0 - zc) * (mX0 - xc);
      }

      // z-derivatives
      b(1, 4) = 0;

      b(2, 4) = 0;
      b(3, 4) = 0;
      b(4, 4) = -1;

      b(5, 4) = 0;
      b(6, 4) = 0;
      b(7, 4) = -2 * (mZ0 - zc);
      b(8, 4) = 0;
      b(9, 4) = -(mX0 - xc);
      b(10, 4) = -(mY0 - yc);

      b(11, 4) = 0;
      b(12, 4) = 0;
      b(13, 4) = -3 * SQR(mZ0 - zc);
      b(14, 4) = 0;
      b(15, 4) = 0;
      b(16, 4) = -2 * (mZ0 - zc) * (mX0 - xc);
      b(17, 4) = -SQR(mX0 - xc);
      b(18, 4) = -SQR(mY0 - yc);
      b(19, 4) = -2 * (mZ0 - zc) * (mY0 - yc);
      b(20, 4) = -(mX0 - xc) * (mY0 - yc);
      if (ncond == 35) {
        b(21, 4) = 0;
        b(22, 4) = 0;
        b(23, 4) = -4 * CUB(mZ0 - zc);
        b(24, 4) = 0;
        b(25, 4) = 0;
        b(26, 4) = -3 * SQR(mZ0 - zc) * (mX0 - xc);
        b(27, 4) = -CUB(mX0 - xc);
        b(28, 4) = -3 * SQR(mZ0 - zc) * (mY0 - yc);
        b(29, 4) = -CUB(mY0 - yc);
        b(30, 4) = 0;
        b(31, 4) = -2 * SQR(mY0 - yc) * (mZ0 - zc);
        b(32, 4) = -2 * (mZ0 - zc) * SQR(mX0 - xc);
        b(33, 4) = -SQR(mX0 - xc) * (mY0 - yc);
        b(34, 4) = -SQR(mY0 - yc) * (mX0 - xc);
        b(35, 4) = -2 * (mZ0 - zc) * (mX0 - xc) * (mY0 - yc);
      }
    }
    // Solve A*x=b for x

    char tr = 'N';
    int ssize = 125, one = 1, info = 0, nb = 20;
    int lwork = ncond + ncond * nb;
    float_sw4* work = new float_sw4[lwork];
    F77_FUNC(dgels, DGELS)
    (tr, ncond, ssize, nrhs, a_, ncond, b_, ldb, work, lwork, info);
    delete[] work;
    REQUIRE2(info == 0, "ERROR, info = " << info << " returned from DGELS ")

    // Define the point sources
    for (int k = kc - 2; k <= kc + 2; k++) {
      if (1 <= k && k <= Nz - 1) {
        float_sw4 nwgh = 1.0;
        if (k <= 4)
          nwgh = normwgh[k - 1];
        else if (k >= Nz - 3)
          nwgh = normwgh[Nz - k];

        // Discretization on this grid
        for (int j = jc - 2; j <= jc + 2; j++)
          for (int i = ic - 2; i <= ic + 2; i++)
            if (a_EW->interior_point_in_proc(i, j, g)) {
              int ind = i - ic + 2 + 5 * (j - jc + 2) + 25 * (k - kc + 2) + 1;
              float_sw4 fx, fy, fz;
              if (gradient) {
                fx = -(mForces[0] * x(ind, 2) + mForces[1] * x(ind, 3) +
                       mForces[2] * x(ind, 4));
                fy = -(mForces[1] * x(ind, 2) + mForces[3] * x(ind, 3) +
                       mForces[4] * x(ind, 4));
                fz = -(mForces[2] * x(ind, 2) + mForces[4] * x(ind, 3) +
                       mForces[5] * x(ind, 4));
              } else {
                fx = mForces[0] * x(ind, 1);
                fy = mForces[1] * x(ind, 1);
                fz = mForces[2] * x(ind, 1);
              }
              if (fx != 0 || fy != 0 || fz != 0) {
                //                     float_sw4
                //                     ijac=1.0/(nwgh*(a_EW->mJ[g](i,j,k)));
                float_sw4 ijac = 1.0 / (nwgh * (Jg(i, j, k)));
                GridPointSource* sourcePtr = new GridPointSource(
                    mFreq, mT0, i, j, k, g, fx * ijac, fy * ijac, fz * ijac,
                    mTimeDependence, mNcyc, mPar, mNpar, mIpar, mNipar);
                point_sources.push_back(sourcePtr);
              }
            }
      }
      if (k <= 0) {
        // Into finer grid above, keep k=1 on this (coarse) grid
        int Nzp = a_EW->m_global_nz[g + 1];
        int kk = k + Nzp - 1;
        float_sw4 nwgh = 1.0;
        if (kk >= Nzp - 3) nwgh = normwgh[Nzp - kk];
        for (int j = jcref - 2; j <= jcref + 2; j++)
          for (int i = icref - 2; i <= icref + 2; i++)
            if (a_EW->interior_point_in_proc(i, j, g + 1)) {
              int ind =
                  i - icref + 2 + 5 * (j - jcref + 2) + 25 * (k - kc + 2) + 1;
              float_sw4 fx, fy, fz;
              if (gradient) {
                fx = -(mForces[0] * x(ind, 2) + mForces[1] * x(ind, 3) +
                       mForces[2] * x(ind, 4));
                fy = -(mForces[1] * x(ind, 2) + mForces[3] * x(ind, 3) +
                       mForces[4] * x(ind, 4));
                fz = -(mForces[2] * x(ind, 2) + mForces[4] * x(ind, 3) +
                       mForces[5] * x(ind, 4));
              } else {
                fx = mForces[0] * x(ind, 1);
                fy = mForces[1] * x(ind, 1);
                fz = mForces[2] * x(ind, 1);
              }
              if (fx != 0 || fy != 0 || fz != 0) {
                //                  float_sw4 wF = x(ind,1);
                //                  if( wF != 0 )
                //                  {
                //                     wF = wF/(nwgh*(a_EW->mJ[g+1](i,j,kk)));
                float_sw4 ijac = 1.0 / (nwgh * Jref(i, j, kk));
                //                     float_sw4 ijac
                //                     = 1.0/(nwgh*(a_EW->mJ[g+1](i,j,kk)));
                //                     std::cout << g+1 << " (i,j,k)" << i << "
                //                     " << j << " " << kk << " " << wF <<
                //                     std::endl;
                GridPointSource* sourcePtr =
                    new GridPointSource(mFreq, mT0, i, j, kk, g + 1, fx * ijac,
                                        fy * ijac, fz * ijac, mTimeDependence,
                                        mNcyc, mPar, mNpar, mIpar, mNipar);
                point_sources.push_back(sourcePtr);
              }
            }
      }
      if (k >= Nz) {
        // Into coarser grid below, move this k=Nz to k=1 on coarse grid.
        int kk = k - Nz + 1;
        float_sw4 nwgh = 1.0;
        if (kk <= 4) nwgh = normwgh[kk - 1];
        for (int j = jcref - 2; j <= jcref + 2; j++)
          for (int i = icref - 2; i <= icref + 2; i++)
            if (a_EW->interior_point_in_proc(i, j, g - 1)) {
              int ind =
                  i - icref + 2 + 5 * (j - jcref + 2) + 25 * (k - kc + 2) + 1;
              float_sw4 fx, fy, fz;
              if (gradient) {
                fx = -(mForces[0] * x(ind, 2) + mForces[1] * x(ind, 3) +
                       mForces[2] * x(ind, 4));
                fy = -(mForces[1] * x(ind, 2) + mForces[3] * x(ind, 3) +
                       mForces[4] * x(ind, 4));
                fz = -(mForces[2] * x(ind, 2) + mForces[4] * x(ind, 3) +
                       mForces[5] * x(ind, 4));
              } else {
                fx = mForces[0] * x(ind, 1);
                fy = mForces[1] * x(ind, 1);
                fz = mForces[2] * x(ind, 1);
              }
              if (fx != 0 || fy != 0 || fz != 0) {
                //                  float_sw4 wF = x(ind,1);
                //                  if( wF != 0 )
                //                  {
                //                     wF = wF/(nwgh*(a_EW->mJ[g-1](i,j,kk)));
                float_sw4 ijac = 1.0 / (nwgh * Jref(i, j, kk));
                //                     float_sw4 ijac
                //                     = 1.0/(nwgh*(a_EW->mJ[g-1](i,j,kk)));
                //                     std::cout << g-1 << " (i,j,k)" << i << "
                //                     " << j << " " << kk << " " << wF <<
                //                     std::endl;
                GridPointSource* sourcePtr =
                    new GridPointSource(mFreq, mT0, i, j, kk, g - 1, fx * ijac,
                                        fy * ijac, fz * ijac, mTimeDependence,
                                        mNcyc, mPar, mNpar, mIpar, mNipar);
                point_sources.push_back(sourcePtr);
              }
            }
      }
    }
    delete[] a_;
    delete[] b_;
  }
  //   delete[] x_;
#undef a
#undef b
#undef x
#undef CUB
#undef BISQR
}

//-----------------------------------------------------------------------
void Source::get_cc_psources(EW* a_EW, int g, float_sw4 q, float_sw4 r,
                             float_sw4 s, bool gradient, float_sw4 normwgh[4],
                             vector<GridPointSource*>& point_sources) {
  // Curvilinear/Cartesian interface with interface conditions imposed.
#define CUB(x) (x) * (x) * (x)
#define BISQR(x) (x) * (x) * (x) * (x)
  //   ic, jc on current grid
  //   icref, jcref on adjacent grid
  // kc
  // Stencil is always -2 <= m <= 2, -2 <= n <= 2, -2 <= o <= 2
  //   std::cout <<"in get_cc_psources " << std::endl;
  int ncond = 35;  // Number of moment conditions
#define a(i, j) a_[(i - 1) + ncond * (j - 1)]
#define b(i, j) b_[(i - 1) + 125 * (j - 1)]
#define x(i, j) b_[(i - 1) + 125 * (j - 1)]
  // setup stencil centers etc.
  int ic = static_cast<int>(round(q));
  int jc = static_cast<int>(round(r));
  int kc = static_cast<int>(round(s));
  int Nz = a_EW->m_global_nz[g];
  int icref = ic, jcref = jc, gref;
  if (g == a_EW->mNumberOfCartesianGrids) {
    // (q,r,s) in curvilinear grid, stencil extends to Cartesian grid below
    gref = g - 1;
  } else {
    // (q,r,s) in Cartesian grid, stencil extends to curvilinear grid above
    gref = g + 1;
  }
  float_sw4 h = a_EW->mGridSize[g];  // h is the same on both grids.

  // Processor owns part of source stencil if (ic,jc) is inside processor.
  // If (ic,jc) is in the processor overlap, then all the stencil widths
  // ic-2,..,ic+2 and jc-2,..,jc+2 are not guaranteed to be owned by this
  // processor. Need to supply grid z and Jacobian for all stencil points,
  // extended arrays are being used.
  if (a_EW->interior_point_in_proc(ic - 2, jc - 2, g) ||
      a_EW->interior_point_in_proc(ic + 2, jc - 2, g) ||
      a_EW->interior_point_in_proc(ic - 2, jc + 2, g) ||
      a_EW->interior_point_in_proc(ic + 2, jc + 2, g) ||
      a_EW->interior_point_in_proc(icref - 2, jcref - 2, gref) ||
      a_EW->interior_point_in_proc(icref + 2, jcref - 2, gref) ||
      a_EW->interior_point_in_proc(icref - 2, jcref + 2, gref) ||
      a_EW->interior_point_in_proc(icref + 2, jcref + 2, gref)) {
    Sarray zg, Jg, zgref, Jgref;
    int kll, kul;
    if (g == a_EW->mNumberOfCartesianGrids) {
      // g is curvilinear, gref Cartesian
      zg.define(ic - 2, ic + 2, jc - 2, jc + 2, kc - 2, kc + 2);
      Jg.define(ic - 2, ic + 2, jc - 2, jc + 2, kc - 2, kc + 2);
      if (!a_EW->interior_point_in_proc(ic, jc, g))
        a_EW->m_gridGenerator->generate_z_and_j(a_EW, g, zg, Jg);
      else {
        zg.insert_intersection(a_EW->mZ[g]);
        Jg.insert_intersection(a_EW->mJ[g]);
      }
      int kll = 1;
      int kul = 5;
      zgref.define(icref - 2, icref + 2, jcref - 2, jcref + 2, kll, kul);
      Jgref.define(icref - 2, icref + 2, jcref - 2, jcref + 2, kll, kul);
      float_sw4 zmin = a_EW->m_zmin[gref];
      for (int k = kll; k <= kul; k++)
        for (int j = jcref - 2; j <= jcref + 2; j++)
          for (int i = icref - 2; i <= icref + 2; i++) {
            zgref(i, j, k) = zmin + (k - 1) * h;
            Jgref(i, j, k) = h * h * h;
          }

    } else {
      // g is Cartesian, gref curvilinear
      zg.define(ic - 2, ic + 2, jc - 2, jc + 2, kc - 2, kc + 2);
      Jg.define(ic - 2, ic + 2, jc - 2, jc + 2, kc - 2, kc + 2);
      float_sw4 zmin = a_EW->m_zmin[g];
      for (int k = kc - 2; k <= kc + 2; k++)
        for (int j = jc - 2; j <= jc + 2; j++)
          for (int i = ic - 2; i <= ic + 2; i++) {
            zg(i, j, k) = zmin + (k - 1) * h;
            Jg(i, j, k) = h * h * h;
          }

      int Nzp = a_EW->m_global_nz[gref];
      kll = Nzp - 4;
      kul = Nzp;
      zgref.define(icref - 2, icref + 2, jcref - 2, jcref + 2, kll, kul);
      Jgref.define(icref - 2, icref + 2, jcref - 2, jcref + 2, kll, kul);
      if (!a_EW->interior_point_in_proc(icref, jcref, gref))
        a_EW->m_gridGenerator->generate_z_and_j(a_EW, gref, zgref, Jgref);
      else {
        zgref.insert_intersection(a_EW->mZ[gref]);
        Jgref.insert_intersection(a_EW->mJ[gref]);
      }
    }

    int nrhs = 1;
    if (gradient) nrhs = 4;

    float_sw4* a_ = new float_sw4[ncond * 125];
    int ldb = ncond > 125 ? ncond : 125;
    float_sw4* b_ = new float_sw4[ldb * nrhs];
    //   float_sw4* x_ = new float_sw4[125];

    //   std::cout << "SOURCE at interface, g= " << g << " (ic,jc,kc)= " << ic
    //   <<
    //      ", " << jc << ", " << kc << std::endl;
    float_sw4 xc = h * (ic - 1);
    float_sw4 yc = h * (jc - 1);
    float_sw4 zc = zg(ic, jc, kc);

    for (int k = kc - 2; k <= kc + 2; k++) {
      if (1 <= k && k <= Nz - 1) {
        // Discretization on this grid
        for (int j = jc - 2; j <= jc + 2; j++)
          for (int i = ic - 2; i <= ic + 2; i++) {
            int ind = i - ic + 2 + 5 * (j - jc + 2) + 25 * (k - kc + 2) + 1;
            float_sw4 x = (i - 1) * h, y = (j - 1) * h, z = zg(i, j, k);
            a(1, ind) = 1;
            a(2, ind) = x - xc;
            a(3, ind) = y - yc;
            a(4, ind) = z - zc;
            a(5, ind) = SQR(x - xc);
            a(6, ind) = SQR(y - yc);
            a(7, ind) = SQR(z - zc);
            a(8, ind) = (x - xc) * (y - yc);
            a(9, ind) = (x - xc) * (z - zc);
            a(10, ind) = (y - yc) * (z - zc);
            a(11, ind) = CUB(x - xc);
            a(12, ind) = CUB(y - yc);
            a(13, ind) = CUB(z - zc);
            a(14, ind) = SQR(x - xc) * (y - yc);
            a(15, ind) = SQR(y - yc) * (x - xc);
            a(16, ind) = SQR(z - zc) * (x - xc);
            a(17, ind) = SQR(x - xc) * (z - zc);
            a(18, ind) = SQR(y - yc) * (z - zc);
            a(19, ind) = SQR(z - zc) * (y - yc);
            a(20, ind) = (x - xc) * (y - yc) * (z - zc);
            if (ncond == 35) {
              a(21, ind) = BISQR(x - xc);
              a(22, ind) = BISQR(y - yc);
              a(23, ind) = BISQR(z - zc);
              a(24, ind) = CUB(x - xc) * (y - yc);
              a(25, ind) = CUB(y - yc) * (x - xc);
              a(26, ind) = CUB(z - zc) * (x - xc);
              a(27, ind) = CUB(x - xc) * (z - zc);
              a(28, ind) = CUB(z - zc) * (y - yc);
              a(29, ind) = CUB(y - yc) * (z - zc);
              a(30, ind) = SQR(x - xc) * SQR(y - yc);
              a(31, ind) = SQR(y - yc) * SQR(z - zc);
              a(32, ind) = SQR(z - zc) * SQR(x - xc);
              a(33, ind) = SQR(x - xc) * (y - yc) * (z - zc);
              a(34, ind) = SQR(y - yc) * (x - xc) * (z - zc);
              a(35, ind) = SQR(z - zc) * (x - xc) * (y - yc);
            }
          }
      }
      if (k <= 0 || k >= Nz) {
        int kk;
        if (k <= 0) {
          // Into finer grid above, continue into fine grid from k=0 and
          // upward --> from kk=Nz-1
          int Nzp = a_EW->m_global_nz[g + 1];
          kk = k + Nzp - 1;
        }
        if (k >= Nz) {
          // Into coarser grid below, let k=Nz be the coarse grid
          // discretization, and continue on coarse grid from k=Nz+1 --> kk=2
          kk = k - Nz + 1;
        }
        for (int j = jcref - 2; j <= jcref + 2; j++)
          for (int i = icref - 2; i <= icref + 2; i++) {
            float_sw4 x = (i - 1) * h, y = (j - 1) * h, z = zgref(i, j, kk);
            int ind =
                i - icref + 2 + 5 * (j - jcref + 2) + 25 * (k - kc + 2) + 1;
            a(1, ind) = 1;
            a(2, ind) = x - xc;
            a(3, ind) = y - yc;
            a(4, ind) = z - zc;
            a(5, ind) = SQR(x - xc);
            a(6, ind) = SQR(y - yc);
            a(7, ind) = SQR(z - zc);
            a(8, ind) = (x - xc) * (y - yc);
            a(9, ind) = (x - xc) * (z - zc);
            a(10, ind) = (y - yc) * (z - zc);
            a(11, ind) = CUB(x - xc);
            a(12, ind) = CUB(y - yc);
            a(13, ind) = CUB(z - zc);
            a(14, ind) = SQR(x - xc) * (y - yc);
            a(15, ind) = SQR(y - yc) * (x - xc);
            a(16, ind) = SQR(z - zc) * (x - xc);
            a(17, ind) = SQR(x - xc) * (z - zc);
            a(18, ind) = SQR(y - yc) * (z - zc);
            a(19, ind) = SQR(z - zc) * (y - yc);
            a(20, ind) = (x - xc) * (y - yc) * (z - zc);
            if (ncond == 35) {
              a(21, ind) = BISQR(x - xc);
              a(22, ind) = BISQR(y - yc);
              a(23, ind) = BISQR(z - zc);
              a(24, ind) = CUB(x - xc) * (y - yc);
              a(25, ind) = CUB(y - yc) * (x - xc);
              a(26, ind) = CUB(z - zc) * (x - xc);
              a(27, ind) = CUB(x - xc) * (z - zc);
              a(28, ind) = CUB(z - zc) * (y - yc);
              a(29, ind) = CUB(y - yc) * (z - zc);
              a(30, ind) = SQR(x - xc) * SQR(y - yc);
              a(31, ind) = SQR(y - yc) * SQR(z - zc);
              a(32, ind) = SQR(z - zc) * SQR(x - xc);
              a(33, ind) = SQR(x - xc) * (y - yc) * (z - zc);
              a(34, ind) = SQR(y - yc) * (x - xc) * (z - zc);
              a(35, ind) = SQR(z - zc) * (x - xc) * (y - yc);
            }
          }
      }
    }

#undef mom_f
    // Right hand sides for moment conditions

    b(1, 1) = 1;
    b(2, 1) = mX0 - xc;
    b(3, 1) = mY0 - yc;
    b(4, 1) = mZ0 - zc;
    b(5, 1) = SQR(mX0 - xc);
    b(6, 1) = SQR(mY0 - yc);
    b(7, 1) = SQR(mZ0 - zc);
    b(8, 1) = (mX0 - xc) * (mY0 - yc);
    b(9, 1) = (mX0 - xc) * (mZ0 - zc);
    b(10, 1) = (mY0 - yc) * (mZ0 - zc);
    b(11, 1) = CUB(mX0 - xc);
    b(12, 1) = CUB(mY0 - yc);
    b(13, 1) = CUB(mZ0 - zc);
    b(14, 1) = SQR(mX0 - xc) * (mY0 - yc);
    b(15, 1) = SQR(mY0 - yc) * (mX0 - xc);
    b(16, 1) = SQR(mZ0 - zc) * (mX0 - xc);
    b(17, 1) = SQR(mX0 - xc) * (mZ0 - zc);
    b(18, 1) = SQR(mY0 - yc) * (mZ0 - zc);
    b(19, 1) = SQR(mZ0 - zc) * (mY0 - yc);
    b(20, 1) = (mX0 - xc) * (mY0 - yc) * (mZ0 - zc);
    if (ncond == 35) {
      b(21, 1) = BISQR(mX0 - xc);
      b(22, 1) = BISQR(mY0 - yc);
      b(23, 1) = BISQR(mZ0 - zc);
      b(24, 1) = CUB(mX0 - xc) * (mY0 - yc);
      b(25, 1) = CUB(mY0 - yc) * (mX0 - xc);
      b(26, 1) = CUB(mZ0 - zc) * (mX0 - xc);
      b(27, 1) = CUB(mX0 - xc) * (mZ0 - zc);
      b(28, 1) = CUB(mZ0 - zc) * (mY0 - yc);
      b(29, 1) = CUB(mY0 - yc) * (mZ0 - zc);
      b(30, 1) = SQR(mX0 - xc) * SQR(mY0 - yc);
      b(31, 1) = SQR(mY0 - yc) * SQR(mZ0 - zc);
      b(32, 1) = SQR(mZ0 - zc) * SQR(mX0 - xc);
      b(33, 1) = SQR(mX0 - xc) * (mY0 - yc) * (mZ0 - zc);
      b(34, 1) = SQR(mY0 - yc) * (mX0 - xc) * (mZ0 - zc);
      b(35, 1) = SQR(mZ0 - zc) * (mX0 - xc) * (mY0 - yc);
    }
    if (gradient) {
      // x-derivatives
      b(1, 2) = 0;

      b(2, 2) = -1;
      b(3, 2) = 0;
      b(4, 2) = 0;

      b(5, 2) = -2 * (mX0 - xc);
      b(6, 2) = 0;
      b(7, 2) = 0;
      b(8, 2) = -(mY0 - yc);
      b(9, 2) = -(mZ0 - zc);
      b(10, 2) = 0;

      b(11, 2) = -3 * SQR(mX0 - xc);
      b(12, 2) = 0;
      b(13, 2) = 0;
      b(14, 2) = -2 * (mX0 - xc) * (mY0 - yc);
      b(15, 2) = -SQR(mY0 - yc);
      b(16, 2) = -SQR(mZ0 - zc);
      b(17, 2) = -2 * (mX0 - xc) * (mZ0 - zc);
      b(18, 2) = 0;
      b(19, 2) = 0;
      b(20, 2) = -(mY0 - yc) * (mZ0 - zc);
      if (ncond == 35) {
        b(21, 2) = -4 * CUB(mX0 - xc);
        b(22, 2) = 0;
        b(23, 2) = 0;
        b(24, 2) = -3 * SQR(mX0 - xc) * (mY0 - yc);
        b(25, 2) = -CUB(mY0 - yc);
        b(26, 2) = -CUB(mZ0 - zc);
        b(27, 2) = -3 * SQR(mX0 - xc) * (mZ0 - zc);
        b(28, 2) = 0;
        b(29, 2) = 0;
        b(30, 2) = -2 * (mX0 - xc) * SQR(mY0 - yc);
        b(31, 2) = 0;
        b(32, 2) = -2 * SQR(mZ0 - zc) * (mX0 - xc);
        b(33, 2) = -2 * (mX0 - xc) * (mY0 - yc) * (mZ0 - zc);
        b(34, 2) = -SQR(mY0 - yc) * (mZ0 - zc);
        b(35, 2) = -SQR(mZ0 - zc) * (mY0 - yc);
      }

      // y-derivatives
      b(1, 3) = 0;

      b(2, 3) = 0;
      b(3, 3) = -1;
      b(4, 3) = 0;

      b(5, 3) = 0;
      b(6, 3) = -2 * (mY0 - yc);
      b(7, 3) = 0;
      b(8, 3) = -(mX0 - xc);
      b(9, 3) = 0;
      b(10, 3) = -(mZ0 - zc);

      b(11, 3) = 0;
      b(12, 3) = -3 * SQR(mY0 - yc);
      b(13, 3) = 0;
      b(14, 3) = -SQR(mX0 - xc);
      b(15, 3) = -2 * (mY0 - yc) * (mX0 - xc);
      b(16, 3) = 0;
      b(17, 3) = 0;
      b(18, 3) = -2 * (mY0 - yc) * (mZ0 - zc);
      b(19, 3) = -SQR(mZ0 - zc);
      b(20, 3) = -(mX0 - xc) * (mZ0 - zc);
      if (ncond == 35) {
        b(21, 3) = 0;
        b(22, 3) = -4 * CUB(mY0 - yc);
        b(23, 3) = 0;
        b(24, 3) = -CUB(mX0 - xc);
        b(25, 3) = -3 * SQR(mY0 - yc) * (mX0 - xc);
        b(26, 3) = 0;
        b(27, 3) = 0;
        b(28, 3) = -CUB(mZ0 - zc);
        b(29, 3) = -3 * SQR(mY0 - yc) * (mZ0 - zc);
        b(30, 3) = -2 * SQR(mX0 - xc) * (mY0 - yc);
        b(31, 3) = -2 * (mY0 - yc) * SQR(mZ0 - zc);
        b(32, 3) = 0;
        b(33, 3) = -SQR(mX0 - xc) * (mZ0 - zc);
        b(34, 3) = -2 * (mY0 - yc) * (mX0 - xc) * (mZ0 - zc);
        b(35, 3) = -SQR(mZ0 - zc) * (mX0 - xc);
      }

      // z-derivatives
      b(1, 4) = 0;

      b(2, 4) = 0;
      b(3, 4) = 0;
      b(4, 4) = -1;

      b(5, 4) = 0;
      b(6, 4) = 0;
      b(7, 4) = -2 * (mZ0 - zc);
      b(8, 4) = 0;
      b(9, 4) = -(mX0 - xc);
      b(10, 4) = -(mY0 - yc);

      b(11, 4) = 0;
      b(12, 4) = 0;
      b(13, 4) = -3 * SQR(mZ0 - zc);
      b(14, 4) = 0;
      b(15, 4) = 0;
      b(16, 4) = -2 * (mZ0 - zc) * (mX0 - xc);
      b(17, 4) = -SQR(mX0 - xc);
      b(18, 4) = -SQR(mY0 - yc);
      b(19, 4) = -2 * (mZ0 - zc) * (mY0 - yc);
      b(20, 4) = -(mX0 - xc) * (mY0 - yc);
      if (ncond == 35) {
        b(21, 4) = 0;
        b(22, 4) = 0;
        b(23, 4) = -4 * CUB(mZ0 - zc);
        b(24, 4) = 0;
        b(25, 4) = 0;
        b(26, 4) = -3 * SQR(mZ0 - zc) * (mX0 - xc);
        b(27, 4) = -CUB(mX0 - xc);
        b(28, 4) = -3 * SQR(mZ0 - zc) * (mY0 - yc);
        b(29, 4) = -CUB(mY0 - yc);
        b(30, 4) = 0;
        b(31, 4) = -2 * SQR(mY0 - yc) * (mZ0 - zc);
        b(32, 4) = -2 * (mZ0 - zc) * SQR(mX0 - xc);
        b(33, 4) = -SQR(mX0 - xc) * (mY0 - yc);
        b(34, 4) = -SQR(mY0 - yc) * (mX0 - xc);
        b(35, 4) = -2 * (mZ0 - zc) * (mX0 - xc) * (mY0 - yc);
      }
    }
    // Solve A*x=b for x

    char tr = 'N';
    int ssize = 125, one = 1, info = 0, nb = 20;
    int lwork = ncond + ncond * nb;
    float_sw4* work = new float_sw4[lwork];
    F77_FUNC(dgels, DGELS)
    (tr, ncond, ssize, nrhs, a_, ncond, b_, ldb, work, lwork, info);
    delete[] work;
    REQUIRE2(info == 0, "ERROR, info = " << info << " returned from DGELS ")

    // Define the point sources
    for (int k = kc - 2; k <= kc + 2; k++) {
      if (1 <= k && k <= Nz - 1) {
        float_sw4 nwgh = 1.0;
        if (k <= 4)
          nwgh = normwgh[k - 1];
        else if (k >= Nz - 3)
          nwgh = normwgh[Nz - k];

        // Discretization on this grid
        for (int j = jc - 2; j <= jc + 2; j++)
          for (int i = ic - 2; i <= ic + 2; i++)
            if (a_EW->interior_point_in_proc(i, j, g)) {
              int ind = i - ic + 2 + 5 * (j - jc + 2) + 25 * (k - kc + 2) + 1;
              float_sw4 fx, fy, fz;
              if (gradient) {
                fx = -(mForces[0] * x(ind, 2) + mForces[1] * x(ind, 3) +
                       mForces[2] * x(ind, 4));
                fy = -(mForces[1] * x(ind, 2) + mForces[3] * x(ind, 3) +
                       mForces[4] * x(ind, 4));
                fz = -(mForces[2] * x(ind, 2) + mForces[4] * x(ind, 3) +
                       mForces[5] * x(ind, 4));
              } else {
                fx = mForces[0] * x(ind, 1);
                fy = mForces[1] * x(ind, 1);
                fz = mForces[2] * x(ind, 1);
              }
              if (fx != 0 || fy != 0 || fz != 0) {
                float_sw4 ijac = 1.0 / (nwgh * (Jg(i, j, k)));
                if (k == 1)  // On interface, special
                  ijac = 1.0 / (nwgh *
                                (Jg(i, j, k) +
                                 Jgref(i, j, a_EW->m_global_nz[gref] + k - 1)));
                GridPointSource* sourcePtr = new GridPointSource(
                    mFreq, mT0, i, j, k, g, fx * ijac, fy * ijac, fz * ijac,
                    mTimeDependence, mNcyc, mPar, mNpar, mIpar, mNipar);
                point_sources.push_back(sourcePtr);
              }
            }
      }
      if (k <= 0) {
        // Into finer grid above, keep k=1 on this (coarse) grid
        int Nzp = a_EW->m_global_nz[g + 1];
        int kk = k + Nzp - 1;
        float_sw4 nwgh = 1.0;
        if (kk >= Nzp - 3) nwgh = normwgh[Nzp - kk];
        for (int j = jcref - 2; j <= jcref + 2; j++)
          for (int i = icref - 2; i <= icref + 2; i++)
            if (a_EW->interior_point_in_proc(i, j, g + 1)) {
              int ind =
                  i - icref + 2 + 5 * (j - jcref + 2) + 25 * (k - kc + 2) + 1;
              float_sw4 fx, fy, fz;
              if (gradient) {
                fx = -(mForces[0] * x(ind, 2) + mForces[1] * x(ind, 3) +
                       mForces[2] * x(ind, 4));
                fy = -(mForces[1] * x(ind, 2) + mForces[3] * x(ind, 3) +
                       mForces[4] * x(ind, 4));
                fz = -(mForces[2] * x(ind, 2) + mForces[4] * x(ind, 3) +
                       mForces[5] * x(ind, 4));
              } else {
                fx = mForces[0] * x(ind, 1);
                fy = mForces[1] * x(ind, 1);
                fz = mForces[2] * x(ind, 1);
              }
              if (fx != 0 || fy != 0 || fz != 0) {
                //                  float_sw4 wF = x(ind,1);
                //                  if( wF != 0 )
                //                  {
                //                     wF = wF/(nwgh*(a_EW->mJ[g+1](i,j,kk)));
                float_sw4 ijac = 1.0 / (nwgh * (a_EW->mJ[g + 1](i, j, kk)));
                //                     std::cout << g+1 << " (i,j,k)" << i << "
                //                     " << j << " " << kk << " " << wF <<
                //                     std::endl;
                GridPointSource* sourcePtr =
                    new GridPointSource(mFreq, mT0, i, j, kk, g + 1, fx * ijac,
                                        fy * ijac, fz * ijac, mTimeDependence,
                                        mNcyc, mPar, mNpar, mIpar, mNipar);
                point_sources.push_back(sourcePtr);
              }
            }
      }
      if (k >= Nz) {
        // Into coarser grid below, move this k=Nz to k=1 on coarse grid.
        int kk = k - Nz + 1;
        float_sw4 nwgh = 1.0;
        if (kk <= 4) nwgh = normwgh[kk - 1];
        for (int j = jcref - 2; j <= jcref + 2; j++)
          for (int i = icref - 2; i <= icref + 2; i++)
            if (a_EW->interior_point_in_proc(i, j, g - 1)) {
              int ind =
                  i - icref + 2 + 5 * (j - jcref + 2) + 25 * (k - kc + 2) + 1;
              float_sw4 fx, fy, fz;
              if (gradient) {
                fx = -(mForces[0] * x(ind, 2) + mForces[1] * x(ind, 3) +
                       mForces[2] * x(ind, 4));
                fy = -(mForces[1] * x(ind, 2) + mForces[3] * x(ind, 3) +
                       mForces[4] * x(ind, 4));
                fz = -(mForces[2] * x(ind, 2) + mForces[4] * x(ind, 3) +
                       mForces[5] * x(ind, 4));
              } else {
                fx = mForces[0] * x(ind, 1);
                fy = mForces[1] * x(ind, 1);
                fz = mForces[2] * x(ind, 1);
              }
              if (fx != 0 || fy != 0 || fz != 0) {
                //                  float_sw4 wF = x(ind,1);
                //                  if( wF != 0 )
                //                  {
                //                     wF = wF/(nwgh*(a_EW->mJ[g-1](i,j,kk)));
                float_sw4 ijac = 1.0 / (nwgh * (Jgref(i, j, kk)));
                //                     std::cout << g-1 << " (i,j,k)" << i << "
                //                     " << j << " " << kk << " " << wF <<
                //                     std::endl;
                if (k == Nz)  // on interface, special
                  ijac = 1.0 / (nwgh * (Jg(i, j, k) + Jgref(i, j, kk)));
                GridPointSource* sourcePtr =
                    new GridPointSource(mFreq, mT0, i, j, kk, g - 1, fx * ijac,
                                        fy * ijac, fz * ijac, mTimeDependence,
                                        mNcyc, mPar, mNpar, mIpar, mNipar);
                point_sources.push_back(sourcePtr);
              }
            }
      }
    }
    delete[] a_;
    delete[] b_;
  }
  //   delete[] x_;
#undef a
#undef b
#undef x
#undef CUB
#undef BISQR
}
