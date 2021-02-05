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
#include <fcntl.h>
#include <math.h>
#include <unistd.h>

#include <cstring>
#include <ctime>

#include "EW.h"
#include "Require.h"
#include "SfileOutput.h"
#include "mpi.h"

// static variable definition (in class only declaration):
int SfileOutput::mPreceedZeros = 0;

SfileOutput* SfileOutput::nil = static_cast<SfileOutput*>(0);

//-----------------------------------------------------------------------
SfileOutput::SfileOutput(EW* a_ew, float_sw4 time, float_sw4 timeInterval,
                         int cycle, int cycleInterval, float_sw4 tstart,
                         const std::string& filePrefix, int sampleFactorH,
                         int sampleFactorV, bool doubleMode)
    : mTime(time),
      mEW(a_ew),
      m_time_done(false),
      mTimeInterval(timeInterval),
      mWritingCycle(cycle),
      mCycleInterval(cycleInterval),
      mFilePrefix(filePrefix),
      mSampleH(sampleFactorH),
      mSampleV(sampleFactorV),
      mFileName(""),
      mNextTime(0.0),
      mStartTime(tstart),
      m_isDefinedMPIWriters(false),
      m_double(doubleMode),
      m_winallocated(false),
      m_isCreated(false),
      m_memallocated(false) {
  m_att = mEW->usingAttenuation();
}

//-----------------------------------------------------------------------
SfileOutput::~SfileOutput() {
  if (m_memallocated) {
    if (m_double)
      for (unsigned int g = 0; g < m_doubleField.size(); g++)
        delete[] m_doubleField[g];
    else
      for (unsigned int g = 0; g < m_floatField.size(); g++)
        delete[] m_floatField[g];
  }
  if (m_winallocated)
    for (int g = 0; g < mEW->mNumberOfGrids; g++) {
      delete[] mWindow[g];
      delete[] mGlobalDims[g];
    }
}

//-----------------------------------------------------------------------
// void SfileOutput::set_start_time(double tStart)
//{
//   mStartTime = tStart;
//}

//-----------------------------------------------------------------------
void SfileOutput::setup_images() {
  if (!m_winallocated) {
    mWindow.resize(mEW->mNumberOfGrids);
    mGlobalDims.resize(mEW->mNumberOfGrids);
    for (int g = 0; g < mEW->mNumberOfGrids; g++) {
      mWindow[g] = new int[6];
      mGlobalDims[g] = new int[6];
    }
    m_winallocated = true;
  }
  for (int g = 0; g < mEW->mNumberOfGrids; g++) {
    mWindow[g][0] = mEW->m_iStartInt[g];
    mWindow[g][1] = mEW->m_iEndInt[g];
    mWindow[g][2] = mEW->m_jStartInt[g];
    mWindow[g][3] = mEW->m_jEndInt[g];
    mWindow[g][4] = mEW->m_kStartInt[g];
    mWindow[g][5] = mEW->m_kEndInt[g];

    mGlobalDims[g][0] = 1;
    mGlobalDims[g][1] = mEW->m_global_nx[g];
    mGlobalDims[g][2] = 1;
    mGlobalDims[g][3] = mEW->m_global_ny[g];
    mGlobalDims[g][4] = mWindow[g][4];
    mGlobalDims[g][5] = mWindow[g][5];

    // For now, the 3D image is assumed to span the entire computational domain
    m_ihavearray.resize(mEW->mNumberOfGrids);
    m_ihavearray[g] = true;

    // Coarsen grid a number of levels
    if (mSampleH > 1 && m_ihavearray[g]) {
      for (int dim = 0; dim < 2; dim++) {
        int remainder =
            (mWindow[g][2 * dim] - mGlobalDims[g][2 * dim]) % mSampleH;
        if (remainder != 0) mWindow[g][2 * dim] += mSampleH - remainder;
        remainder =
            (mWindow[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) % mSampleH;
        if (remainder != 0) mWindow[g][2 * dim + 1] -= remainder;

        remainder =
            (mGlobalDims[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) % mSampleH;
        if (remainder != 0) mGlobalDims[g][2 * dim + 1] -= remainder;
      }
      for (int dim = 2; dim < 3; dim++) {
        int remainder =
            (mWindow[g][2 * dim] - mGlobalDims[g][2 * dim]) % mSampleV;
        if (remainder != 0) mWindow[g][2 * dim] += mSampleV - remainder;
        remainder =
            (mWindow[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) % mSampleV;
        if (remainder != 0) mWindow[g][2 * dim + 1] -= remainder;

        remainder =
            (mGlobalDims[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) % mSampleV;
        if (remainder != 0) mGlobalDims[g][2 * dim + 1] -= remainder;
      }
    }
  }

  // If there is a glitch, avoid it by adding an extra point to the z-direction
  m_extraz.resize(mEW->mNumberOfGrids);
  m_extraz[0] = 0;
  // Feature disabled
  //   for( int g=1 ; g < mEW->mNumberOfGrids ; g++ )
  //   {
  //      m_extraz[g] = 0;
  //      if( mEW->m_zmin[g] + (mWindow[g][5]-1)*mEW->mGridSize[g] <
  //          mEW->m_zmin[g-1] + (mWindow[g-1][4]-1)*mEW->mGridSize[g-1] )
  //	 m_extraz[g] = 1;
  //   }

  if (m_double) {
    m_doubleField.resize(mEW->mNumberOfGrids);
    for (int g = 0; g < mEW->mNumberOfGrids; g++) {
      if (m_ihavearray[g]) {
        size_t npts =
            ((size_t)(mWindow[g][1] - mWindow[g][0]) / mSampleH + 1) *
            ((mWindow[g][3] - mWindow[g][2]) / mSampleH + 1) *
            ((mWindow[g][5] - mWindow[g][4]) / mSampleV + 1 + m_extraz[g]);
        m_doubleField[g] = new double[npts];
      } else
        m_doubleField[g] = new double[1];
    }
  } else {
    m_floatField.resize(mEW->mNumberOfGrids);
    for (int g = 0; g < mEW->mNumberOfGrids; g++) {
      if (m_ihavearray[g]) {
        size_t npts =
            ((size_t)(mWindow[g][1] - mWindow[g][0]) / mSampleH + 1) *
            ((mWindow[g][3] - mWindow[g][2]) / mSampleH + 1) *
            ((mWindow[g][5] - mWindow[g][4]) / mSampleV + 1 + m_extraz[g]);
        m_floatField[g] = new float[npts];
      } else
        m_floatField[g] = new float[1];
    }
  }
  m_memallocated = true;
  define_pio();
}

//-----------------------------------------------------------------------
void SfileOutput::define_pio() {
  int glow = 0, ghigh = mEW->mNumberOfGrids;
  m_parallel_io = new Parallel_IO*[ghigh - glow + 1];
  for (int g = glow; g < ghigh; g++) {
    int global[3], local[3], start[3];
    for (int dim = 0; dim < 2; dim++) {
      global[dim] =
          (mGlobalDims[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) / mSampleH +
          1;
      local[dim] =
          (mWindow[g][2 * dim + 1] - mWindow[g][2 * dim]) / mSampleH + 1;
      start[dim] = (mWindow[g][2 * dim] - mGlobalDims[g][2 * dim]) / mSampleH;
    }
    for (int dim = 2; dim < 3; dim++) {
      global[dim] =
          (mGlobalDims[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) / mSampleV +
          1;
      local[dim] =
          (mWindow[g][2 * dim + 1] - mWindow[g][2 * dim]) / mSampleV + 1;
      start[dim] = (mWindow[g][2 * dim] - mGlobalDims[g][2 * dim]) / mSampleV;
    }
    global[2] += m_extraz[g];
    local[2] += m_extraz[g];

    int iwrite = 0;
    int nrwriters = mEW->getNumberOfWritersPFS();
    int nproc = 0, myid = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // new hack
    int* owners = new int[nproc];
    int i = 0;
    for (int p = 0; p < nproc; p++)
      if (m_ihavearray[g]) owners[i++] = p;
    if (nrwriters > i) nrwriters = i;

    if (nrwriters > nproc) nrwriters = nproc;
    int q, r;
    if (nproc == 1 || nrwriters == 1) {
      q = 0;
      r = 0;
    } else {
      q = (nproc - 1) / (nrwriters - 1);
      r = (nproc - 1) % (nrwriters - 1);
    }
    for (int w = 0; w < nrwriters; w++)
      if (q * w + r == myid) iwrite = 1;
    //      std::cout << "Define PIO: grid " << g << " myid = " << myid << "
    //      iwrite= " << iwrite << " start= "
    //		<< start[0] << " " << start[1] << " " << start[2] << std::endl;
    m_parallel_io[g - glow] =
        new Parallel_IO(iwrite, mEW->usingParallelFS(), global, local, start);
    delete[] owners;
  }
  m_isDefinedMPIWriters = true;
}

//-----------------------------------------------------------------------
void SfileOutput::setSteps(int a_steps) {
  char buffer[50];
  mPreceedZeros = snprintf(buffer, 50, "%d", a_steps);
}

//-----------------------------------------------------------------------
// void SfileOutput::set_double(bool val)
//{
//   m_double = val;
//}

//-----------------------------------------------------------------------
bool SfileOutput::timeToWrite(float_sw4 time, int cycle, float_sw4 dt) {
  // -----------------------------------------------
  // Check based on cycle
  // -----------------------------------------------
  //   cout << "in time to write " << mWritingCycle << " " << mCycleInterval <<
  //   " " << " " << mTime << " " <<  mTimeInterval << " " << endl;
  bool do_it = false;
  if (cycle == mWritingCycle) do_it = true;
  if (mCycleInterval != 0 && cycle % mCycleInterval == 0 && time >= mStartTime)
    do_it = true;

  // ---------------------------------------------------
  // Check based on time
  // ---------------------------------------------------
  if (mTime > 0.0 && (mTime <= time + dt * 0.5) && !m_time_done) {
    m_time_done = true;
    do_it = true;
  }
  if (mTimeInterval != 0.0 && mNextTime <= time + dt * 0.5) {
    mNextTime += mTimeInterval;
    if (time >= mStartTime) do_it = true;
  }
  return do_it;
}

//-----------------------------------------------------------------------
void SfileOutput::update_image(
    int a_cycle, float_sw4 a_time, float_sw4 a_dt, std::vector<Sarray>& a_U,
    std::vector<Sarray>& a_Rho, std::vector<Sarray>& a_Mu,
    std::vector<Sarray>& a_Lambda, std::vector<Sarray>& a_gRho,
    std::vector<Sarray>& a_gMu, std::vector<Sarray>& a_gLambda,
    std::vector<Sarray>& a_Qp, std::vector<Sarray>& a_Qs, std::string a_path,
    std::vector<Sarray>& a_Z) {
  if (timeToWrite(a_time, a_cycle, a_dt))
    force_write_image(a_time, a_cycle, a_U, a_Rho, a_Mu, a_Lambda, a_gRho,
                      a_gMu, a_gLambda, a_Qp, a_Qs, a_path, a_Z);
}

//-----------------------------------------------------------------------
void SfileOutput::force_write_image(
    float_sw4 a_time, int a_cycle, vector<Sarray>& a_U, vector<Sarray>& a_Rho,
    vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda, vector<Sarray>& a_gRho,
    vector<Sarray>& a_gMu, vector<Sarray>& a_gLambda, vector<Sarray>& a_Qp,
    vector<Sarray>& a_Qs, std::string a_path, std::vector<Sarray>& a_Z) {
  double stime, etime, ttime = 0.0;
  enum SfileOutputMode vars[5] = {RHO, P, S, QP, QS};
  int nvar = m_att ? 5 : 3;
  std::string fname;
  gen_fname(a_path, a_cycle, fname);

  for (int varid = 0; varid < nvar; varid++) {
    mMode = vars[varid];
    compute_image(a_U, a_Rho, a_Mu, a_Lambda, a_gRho, a_gMu, a_gLambda, a_Qp,
                  a_Qs, a_Z);
    stime = MPI_Wtime();
    write_image(fname.c_str(), a_Z);
    etime = MPI_Wtime();
    ttime += etime - stime;
  }

  if (m_parallel_io[0]->proc_zero())
    cout << "  Done writing sfile, took " << etime - stime << " seconds"
         << std::endl;
  m_isCreated = false;
}

//-----------------------------------------------------------------------
void SfileOutput::compute_image(vector<Sarray>& a_U, vector<Sarray>& a_Rho,
                                vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                                vector<Sarray>& a_gRho, vector<Sarray>& a_gMu,
                                vector<Sarray>& a_gLambda, vector<Sarray>& a_Qp,
                                vector<Sarray>& a_Qs,
                                std::vector<Sarray>& a_Z) {
  int stH = mSampleH;
  int stV = mSampleV;
  double my_z, up_z, down_z, up_v, down_v;

  /* if (mMode== RHO && m_parallel_io[0]->proc_zero()) { */
  /*     for( int g=mEW->mNumberOfCartesianGrids ; g < mEW->mNumberOfGrids ; g++
   * ) { */
  /*       cout << "g="<< g << ", a_Z.m_kb="<< a_Z[g].m_kb << ", a_Z.m_ke=" <<
   * a_Z[g].m_ke << ", gz=" << mWindow[g][5] << endl; */
  /*       cout << "kend, kstart: " << mEW->m_kEndInt[g] << ", " <<
   * mEW->m_kStartInt[g]; */
  /*       printf("Z values\n"); */
  /*       for (int j = 0; j <= a_Z[g].m_ke; j++) { */
  /*         printf("(1,1,%d): %f\n", j, a_Z[g](1,1,j)); */
  /*       } */
  /*     } */
  /* } */

  for (int g = 0; g < mEW->mNumberOfGrids; g++) {
    int nkw = (mWindow[g][5] - mWindow[g][4]) / stV + 1;
    int njkw = nkw * ((mWindow[g][3] - mWindow[g][2]) / stH + 1);
    int gz = g;
    /* int ku = mWindow[gz][5]; */
    int kl = mEW->m_kStartInt[gz];
    int ku = mEW->m_kEndInt[gz];

    Sarray *data1 = NULL, *data2 = NULL, *data3 = NULL;
    if (mMode == RHO)
      data1 = &a_Rho[g];
    else if (mMode == QP)
      data1 = &a_Qp[g];
    else if (mMode == QS)
      data1 = &a_Qs[g];
    else if (mMode == P || mMode == S) {
      data1 = &a_Rho[g];
      data2 = &a_Mu[g];
      data3 = &a_Lambda[g];
    }

    /* printf("g=%d, mWindow (%d, %d, %d, %d, %d, %d)\n", g, mWindow[g][0],
     * mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4],
     * mWindow[g][5]); */
    /* fflush(stdout); */

    if (mMode == RHO || mMode == QP ||
        mMode ==
            QS) {  // these modes just copy the values straight from the array
      if (m_double) {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += stV)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += stH)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += stH) {
              size_t ind = (k - mWindow[g][4]) / stV +
                           nkw * (j - mWindow[g][2]) / stH +
                           njkw * (i - mWindow[g][0]) / stH;
              if (g < mEW->mNumberOfCartesianGrids)
                m_doubleField[g][ind] = (double)((*data1)(1, i, j, k));
              else {
                double z_kl = -mEW->mTopo(i, j, 1);
                double z_ku = a_Z[gz](i, j, ku);
                my_z = z_kl + (z_ku - z_kl) * (k - 1) / (double)(ku - 1);
                int t = 1;
                while (my_z >= a_Z[gz](i, j, t) && t < mWindow[g][5]) {
                  t++;
                }
                up_z = (double)a_Z[gz](i, j, t - 1);
                down_z = (double)a_Z[gz](i, j, t);
                up_v = (double)((*data1)(1, i, j, t - 1));
                down_v = (double)((*data1)(1, i, j, t));
                // Linear interp
                m_doubleField[g][ind] =
                    up_v + (down_v - up_v) * (my_z - up_z) / (down_z - up_z);
              }
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += stV)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += stH)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += stH) {
              size_t ind = (k - mWindow[g][4]) / stV +
                           nkw * (j - mWindow[g][2]) / stH +
                           njkw * (i - mWindow[g][0]) / stH;
              if (g < mEW->mNumberOfCartesianGrids)
                m_floatField[g][ind] = (float)((*data1)(1, i, j, k));
              else {
                // debug
                double z_kl = -mEW->mTopo(i, j, 1);
                double z_ku = a_Z[gz](i, j, ku);
                my_z = z_kl + (z_ku - z_kl) * (k - 1) / (double)(ku - 1);
                int t = 1;
                while (my_z >= a_Z[gz](i, j, t) && t < mWindow[g][5]) {
                  t++;
                }
                up_z = (double)a_Z[gz](i, j, t - 1);
                down_z = (double)a_Z[gz](i, j, t);
                up_v = (double)((*data1)(1, i, j, t - 1));
                down_v = (double)((*data1)(1, i, j, t));
                // Linear interp
                m_floatField[g][ind] =
                    up_v + (down_v - up_v) * (my_z - up_z) / (down_z - up_z);
              }
            }
      }
    } else if (mMode == P) {
      if (m_double) {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += stV)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += stH)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += stH) {
              size_t ind = (k - mWindow[g][4]) / stV +
                           nkw * (j - mWindow[g][2]) / stH +
                           njkw * (i - mWindow[g][0]) / stH;
              if (g < mEW->mNumberOfCartesianGrids)
                m_doubleField[g][ind] = (double)sqrt(
                    (2 * ((*data2)(1, i, j, k)) + ((*data3)(1, i, j, k))) /
                    ((*data1)(1, i, j, k)));
              else {
                double z_kl = -mEW->mTopo(i, j, 1);
                double z_ku = a_Z[gz](i, j, ku);
                my_z = z_kl + (z_ku - z_kl) * (k - 1) / (double)(ku - 1);
                int t = 1;
                while (my_z >= a_Z[gz](i, j, t) && t < mWindow[g][5]) {
                  t++;
                }
                up_z = (double)a_Z[gz](i, j, t - 1);
                down_z = (double)a_Z[gz](i, j, t);

                double down_mu = (*data2)(1, i, j, t),
                       down_lambda = (*data3)(1, i, j, t),
                       up_mu = (*data2)(1, i, j, t - 1),
                       up_lambda = (*data3)(1, i, j, t - 1);
                /* if (mEW->usingAttenuation()) { */
                /*  if( NULL == mEW->use_twilight_forcing() ) { */
                /*     mEW->reverse_setup_viscoelastic(gz, i, j, t, down_mu,
                 * down_lambda); */
                /*     mEW->reverse_setup_viscoelastic(gz, i, j, t-1, up_mu,
                 * up_lambda); */
                /*  } */
                /* } */
                up_v = sqrt((2.0 * up_mu + up_lambda) /
                            ((*data1)(1, i, j, t - 1)));
                down_v = sqrt((2.0 * down_mu + down_lambda) /
                              ((*data1)(1, i, j, t)));

                // Linear interp
                m_doubleField[g][ind] =
                    up_v + (down_v - up_v) * (my_z - up_z) / (down_z - up_z);
              }
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += stV)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += stH)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += stH) {
              size_t ind = (k - mWindow[g][4]) / stV +
                           nkw * (j - mWindow[g][2]) / stH +
                           njkw * (i - mWindow[g][0]) / stH;
              if (g < mEW->mNumberOfCartesianGrids)
                m_floatField[g][ind] = (float)sqrt(
                    (2.0 * ((*data2)(1, i, j, k)) + ((*data3)(1, i, j, k))) /
                    ((*data1)(1, i, j, k)));
              else {
                double z_kl = -mEW->mTopo(i, j, 1);
                double z_ku = a_Z[gz](i, j, ku);
                my_z = z_kl + (z_ku - z_kl) * (k - 1) / (double)(ku - 1);
                int t = 1;
                while (my_z >= a_Z[gz](i, j, t) && t < mWindow[g][5]) {
                  t++;
                }
                up_z = (double)a_Z[gz](i, j, t - 1);
                down_z = (double)a_Z[gz](i, j, t);

                double down_mu = (*data2)(1, i, j, t),
                       down_lambda = (*data3)(1, i, j, t),
                       up_mu = (*data2)(1, i, j, t - 1),
                       up_lambda = (*data3)(1, i, j, t - 1);
                /* if (mEW->usingAttenuation()) { */
                /*  if( NULL == mEW->use_twilight_forcing() ) { */
                /*     mEW->reverse_setup_viscoelastic(gz, i, j, t, down_mu,
                 * down_lambda); */
                /*     mEW->reverse_setup_viscoelastic(gz, i, j, t-1, up_mu,
                 * up_lambda); */
                /*     /1* if (i==1&&j==1&&gz==1) *1/ */
                /*     /1*   printf("after reverse (1,1,%d) mu=%f, lambda=%f\n",
                 * t, down_mu, down_lambda); *1/ */
                /*  } */
                /* } */
                up_v = sqrt((2.0 * up_mu + up_lambda) /
                            ((*data1)(1, i, j, t - 1)));
                down_v = sqrt((2.0 * down_mu + down_lambda) /
                              ((*data1)(1, i, j, t)));
                // Linear interp
                m_floatField[g][ind] =
                    up_v + (down_v - up_v) * (my_z - up_z) / (down_z - up_z);
                /* // TODO: debug */
                /* if (i==1&&j==1&&gz==1) { */
                /*     printf("k=%d, my_z=%f, up_z=%f, down_z=%f, up_v[%d]=%f,
                 * down_v[%d]=%f, my_v=%f\n", */
                /*             k, my_z, up_z, down_z, t-1, up_v, t, down_v,
                 * m_floatField[g][ind]); */
                /* } */
              }
            }
      }
    } else if (mMode == S) {
      if (m_double) {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += stV)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += stH)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += stH) {
              size_t ind = (k - mWindow[g][4]) / stV +
                           nkw * (j - mWindow[g][2]) / stH +
                           njkw * (i - mWindow[g][0]) / stH;
              if (g < mEW->mNumberOfCartesianGrids)
                m_doubleField[g][ind] = (double)sqrt(((*data2)(1, i, j, k)) /
                                                     ((*data1)(1, i, j, k)));
              else {
                double z_kl = -mEW->mTopo(i, j, 1);
                double z_ku = a_Z[gz](i, j, ku);
                my_z = z_kl + (z_ku - z_kl) * (k - 1) / (double)(ku - 1);
                int t = 1;
                while (my_z >= a_Z[gz](i, j, t) && t < mWindow[g][5]) {
                  t++;
                }
                up_z = (double)a_Z[gz](i, j, t - 1);
                down_z = (double)a_Z[gz](i, j, t);

                double down_mu = (*data2)(1, i, j, t),
                       down_lambda = (*data3)(1, i, j, t),
                       up_mu = (*data2)(1, i, j, t - 1),
                       up_lambda = (*data3)(1, i, j, t - 1);
                /* if (mEW->usingAttenuation()) { */
                /*   if( NULL == mEW->use_twilight_forcing() ) { */
                /*      mEW->reverse_setup_viscoelastic(gz, i, j, t, down_mu,
                 * down_lambda); */
                /*      mEW->reverse_setup_viscoelastic(gz, i, j, t-1, up_mu,
                 * up_lambda); */
                /*   } */
                /* } */
                up_v = sqrt(up_mu / ((*data1)(1, i, j, t - 1)));
                down_v = sqrt(down_mu / ((*data1)(1, i, j, t)));
                // Linear interp
                m_doubleField[g][ind] =
                    up_v + (down_v - up_v) * (my_z - up_z) / (down_z - up_z);
              }
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += stV)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += stH)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += stH) {
              size_t ind = (k - mWindow[g][4]) / stV +
                           nkw * (j - mWindow[g][2]) / stH +
                           njkw * (i - mWindow[g][0]) / stH;
              if (g < mEW->mNumberOfCartesianGrids)
                m_floatField[g][ind] = (float)sqrt(((*data2)(1, i, j, k)) /
                                                   ((*data1)(1, i, j, k)));
              else {
                double z_kl = -mEW->mTopo(i, j, 1);
                double z_ku = a_Z[gz](i, j, ku);
                my_z = z_kl + (z_ku - z_kl) * (k - 1) / (double)(ku - 1);
                int t = 1;
                while (my_z >= a_Z[gz](i, j, t) && t < mWindow[g][5]) {
                  t++;
                }
                up_z = (double)a_Z[gz](i, j, t - 1);
                down_z = (double)a_Z[gz](i, j, t);

                double down_mu = (*data2)(1, i, j, t),
                       down_lambda = (*data3)(1, i, j, t),
                       up_mu = (*data2)(1, i, j, t - 1),
                       up_lambda = (*data3)(1, i, j, t - 1);
                /* if (mEW->usingAttenuation()) { */
                /*   if( NULL == mEW->use_twilight_forcing() ) { */
                /*      mEW->reverse_setup_viscoelastic(gz, i, j, t, down_mu,
                 * down_lambda); */
                /*      mEW->reverse_setup_viscoelastic(gz, i, j, t-1, up_mu,
                 * up_lambda); */
                /*   } */
                /* } */
                up_v = sqrt(up_mu / ((*data1)(1, i, j, t - 1)));
                down_v = sqrt(down_mu / ((*data1)(1, i, j, t)));

                // Linear interp
                m_floatField[g][ind] =
                    up_v + (down_v - up_v) * (my_z - up_z) / (down_z - up_z);
              }
            }
      }
    }

    // Extrapolate extra point in z
    if (m_extraz[g] == 1) {
      /* cout << "g=" << g << ", need extrapolate extra z" << endl; */
      int k = mWindow[g][5] + stV;
      size_t ind, ind_1, ind_2;  // npts;
#ifdef BZ_DEBUG
      size_t npts = ((size_t)(mWindow[g][1] - mWindow[g][0]) / stH + 1) *
                    ((mWindow[g][3] - mWindow[g][2]) / stH + 1) *
                    ((mWindow[g][5] - mWindow[g][4]) / stV + 1 + m_extraz[g]);
#endif
      // Linear extrapolation assumes that
      // mWindow[g][5]-st>=mWindow[g][4], i.e. mWindow[g][5] -mWindow[g][4]>=st.
      int ok = 2;
      if (mWindow[g][5] - mWindow[g][4] < stV) ok = 1;

      if (m_double) {
        for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += stH)
          for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += stH) {
            ind = (k - mWindow[g][4]) / stV + nkw * (j - mWindow[g][2]) / stH +
                  njkw * (i - mWindow[g][0]) / stH;
            ind_1 = (k - 1 - mWindow[g][4]) / stV +
                    nkw * (j - mWindow[g][2]) / stH +
                    njkw * (i - mWindow[g][0]) / stH;
            ind_2 = (k - ok - mWindow[g][4]) / stV +
                    nkw * (j - mWindow[g][2]) / stH +
                    njkw * (i - mWindow[g][0]) / stH;
#ifdef BZ_DEBUG
            ASSERT(ind < npts);
#endif
            m_doubleField[g][ind] =
                2 * m_doubleField[g][ind_1] - m_doubleField[g][ind_2];
          }

      } else {
        for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += stH)
          for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += stH) {
            ind = (k - mWindow[g][4]) / stV + nkw * (j - mWindow[g][2]) / stH +
                  njkw * (i - mWindow[g][0]) / stH;
            ind_1 = (k - 1 - mWindow[g][4]) / stV +
                    nkw * (j - mWindow[g][2]) / stH +
                    njkw * (i - mWindow[g][0]) / stH;
            ind_2 = (k - ok - mWindow[g][4]) / stV +
                    nkw * (j - mWindow[g][2]) / stH +
                    njkw * (i - mWindow[g][0]) / stH;
#ifdef BZ_DEBUG
            ASSERT(ind < npts);
#endif
            m_floatField[g][ind] =
                2 * m_floatField[g][ind_1] - m_floatField[g][ind_2];
          }
      }
    }
  }
}

//-----------------------------------------------------------------------
void SfileOutput::gen_fname(std::string& path, int cycle, std::string& fname) {
  fname = path;
  fname += mFilePrefix;
  /* fname += ".cycle="; */
  /* int temp = static_cast<int>(pow(10.0, mPreceedZeros - 1)); */
  /* int testcycle = cycle; */
  /* if (cycle == 0) */
  /*   testcycle=1; */
  /* while (testcycle < temp) { */
  /*   fname += "0"; */
  /*   temp /= 10; */
  /* } */
  /* fname += std::to_string(cycle); */
  fname += ".sfile";
}

//-----------------------------------------------------------------------
void SfileOutput::write_image(const char* fname, std::vector<Sarray>& a_Z) {
  std::string m_modestring;
  if (mMode == RHO)
    m_modestring = "Rho";
  else if (mMode == P)
    m_modestring = "Cp";
  else if (mMode == S)
    m_modestring = "Cs";
  else if (mMode == QP)
    m_modestring = "Qp";
  else if (mMode == QS)
    m_modestring = "Qs";

#ifdef USE_HDF5
  hid_t h5_fid, grp, grp2, dset, attr, dspace, attr_space1, attr_space2,
      attr_space3, fapl, dxpl, filespace, memspace;
  int ret;
  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  hsize_t offsets[3], counts[3];
  char dname[128], gname[128];

  int ng = mEW->mNumberOfGrids;

  ASSERT(m_isDefinedMPIWriters);
  int gridinfo = 0;
  if (mEW->topographyExists()) gridinfo = 1;

  bool iwrite = false;
  for (int g = 0; g < ng; g++) iwrite = iwrite || m_parallel_io[g]->i_write();

  int stV = mSampleV;
  int stH = mSampleH;

  int alignment = 65536;
  setenv("HDF5_USE_FILE_LOCKING", "FALSE", 1);

  // Open file from processor zero and write header.
  if (m_parallel_io[0]->proc_zero() && !m_isCreated) {
    hsize_t dims2 = 2, dims3 = 3;

    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_alignment(fapl, alignment, alignment);

    h5_fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    if (h5_fid < 0)
      VERIFY2(0, "ERROR: SfileOutput::write_image, error creating HDF5 file "
                     << fname << " for writing header");

    std::cout << "writing Sfile to " << fname
              << ", with horizontal sample factor=" << stH
              << ", vertical sample factor=" << stV << std::endl;

    attr_space1 = H5Screate(H5S_SCALAR);
    attr_space2 = H5Screate_simple(1, &dims2, NULL);
    attr_space3 = H5Screate_simple(1, &dims3, NULL);

    const char* aname;
    aname = "Origin longitude, latitude, azimuth";
    double lonlataz[3];
    lonlataz[0] = mEW->getLonOrigin();
    lonlataz[1] = mEW->getLatOrigin();
    lonlataz[2] = mEW->getGridAzimuth();
    attr = H5Acreate(h5_fid, aname, H5T_NATIVE_DOUBLE, attr_space3, H5P_DEFAULT,
                     H5P_DEFAULT);
    if (attr < 0)
      VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " << aname);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, lonlataz);
    H5Aclose(attr);

    aname = "Coarsest horizontal grid spacing";
    double spacing = mEW->mGridSize[ng - 1] * stH;
    attr = H5Acreate(h5_fid, aname, H5T_NATIVE_DOUBLE, attr_space1, H5P_DEFAULT,
                     H5P_DEFAULT);
    if (attr < 0)
      VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " << aname);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &spacing);
    H5Aclose(attr);

    aname = "Attenuation";
    int att = mEW->usingAttenuation();
    attr = H5Acreate(h5_fid, aname, H5T_NATIVE_INT, attr_space1, H5P_DEFAULT,
                     H5P_DEFAULT);
    if (attr < 0)
      VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " << aname);
    H5Awrite(attr, H5T_NATIVE_INT, &att);
    H5Aclose(attr);

    aname = "ngrids";
    attr = H5Acreate(h5_fid, aname, H5T_NATIVE_INT, attr_space1, H5P_DEFAULT,
                     H5P_DEFAULT);
    if (attr < 0)
      VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " << aname);
    H5Awrite(attr, H5T_NATIVE_INT, &ng);
    H5Aclose(attr);

    aname = "Min, max depth";
    double minmaxdp[2];
    minmaxdp[0] = mEW->getGlobalZmin();
    minmaxdp[1] = mEW->getGlobalZmax();
    attr = H5Acreate(h5_fid, aname, H5T_NATIVE_DOUBLE, attr_space2, H5P_DEFAULT,
                     H5P_DEFAULT);
    if (attr < 0)
      VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " << aname);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, minmaxdp);
    H5Aclose(attr);

    grp = H5Gcreate(h5_fid, "Z_interfaces", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
    if (grp < 0)
      VERIFY2(0,
              "ERROR: SfileOutput::write_image, error creating Z_interfaces");

    // Top interface (topo)
    hsize_t globalSize[3];
    globalSize[0] =
        (hsize_t)(mGlobalDims[ng - 1][1] - mGlobalDims[ng - 1][0]) / stH + 1;
    globalSize[1] =
        (hsize_t)(mGlobalDims[ng - 1][3] - mGlobalDims[ng - 1][2]) / stH + 1;
    sprintf(dname, "z_values_%d", 0);

    dspace = H5Screate_simple(2, globalSize, NULL);
    hid_t dcpl;
    float intf = 0.0;
    // No topo fill all 0s to top interface
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (gridinfo == 0) H5Pset_fill_value(dcpl, H5T_NATIVE_FLOAT, &intf);
    dset = H5Dcreate(grp, dname, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, dcpl,
                     H5P_DEFAULT);

    H5Pclose(dcpl);
    H5Sclose(dspace);
    H5Dclose(dset);

    // Create interfaces other than top one
    for (int g = ng - 1; g >= 0; g--) {
      sprintf(dname, "z_values_%d", ng - g);
      globalSize[0] =
          (hsize_t)(mGlobalDims[g][1] - mGlobalDims[g][0]) / stH + 1;
      globalSize[1] =
          (hsize_t)(mGlobalDims[g][3] - mGlobalDims[g][2]) / stH + 1;

      dspace = H5Screate_simple(2, globalSize, NULL);
      dcpl = H5Pcreate(H5P_DATASET_CREATE);
      // Cartisian grid, fill with a const value
      if (g <= mEW->mNumberOfCartesianGrids) {
        if (g == 0)
          intf = mEW->getGlobalZmax();
        else
          intf = mEW->m_zmin[g - 1];
        /* std::cout << "Setting const z value to intf " << ng-g << " with " <<
         * intf << std::endl; */
        H5Pset_fill_value(dcpl, H5T_NATIVE_FLOAT, &intf);
      }

      dset = H5Dcreate(grp, dname, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, dcpl,
                       H5P_DEFAULT);

      /* double h = (double) mEW->mGridSize[g]*stH; */
      /* aname = "Horizontal grid size"; */
      /* attr = H5Acreate(dset, aname, H5T_NATIVE_DOUBLE, attr_space1,
       * H5P_DEFAULT, H5P_DEFAULT); */
      /* if( attr < 0 ) */
      /*   VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " <<
       * aname); */
      /* H5Awrite(attr, H5T_NATIVE_DOUBLE, &h); */
      /* H5Aclose(attr); */

      H5Pclose(dcpl);
      H5Sclose(dspace);
      H5Dclose(dset);
    }
    H5Gclose(grp);

    // Create material dsets
    grp = H5Gcreate(h5_fid, "Material_model", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
    if (grp < 0)
      VERIFY2(0,
              "ERROR: SfileOutput::write_image, error creating Material model "
              "group");

    for (int g = 0; g < ng; g++) {
      sprintf(gname, "grid_%d", ng - g - 1);
      grp2 = H5Gcreate(grp, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (grp2 < 0)
        VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " << gname);

      double h = (double)mEW->mGridSize[g] * stH;
      aname = "Horizontal grid size";
      attr = H5Acreate(grp2, aname, H5T_NATIVE_DOUBLE, attr_space1, H5P_DEFAULT,
                       H5P_DEFAULT);
      if (attr < 0)
        VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " << aname);
      H5Awrite(attr, H5T_NATIVE_DOUBLE, &h);
      H5Aclose(attr);

      aname = "Number of components";
      int ncomp = m_att ? 5 : 3;
      attr = H5Acreate(grp2, aname, H5T_NATIVE_INT, attr_space1, H5P_DEFAULT,
                       H5P_DEFAULT);
      if (attr < 0)
        VERIFY2(0, "ERROR: SfileOutput::write_image, error creating " << aname);
      H5Awrite(attr, H5T_NATIVE_INT, &ncomp);
      H5Aclose(attr);

      globalSize[0] =
          (hsize_t)(mGlobalDims[g][1] - mGlobalDims[g][0]) / stH + 1;
      globalSize[1] =
          (hsize_t)(mGlobalDims[g][3] - mGlobalDims[g][2]) / stH + 1;
      globalSize[2] = (hsize_t)(mGlobalDims[g][5] - mGlobalDims[g][4]) / stV +
                      1 + m_extraz[g];
      dspace = H5Screate_simple(3, globalSize, NULL);

      dset = H5Dcreate(grp2, "Cp", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(dset);
      dset = H5Dcreate(grp2, "Cs", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(dset);
      dset = H5Dcreate(grp2, "Rho", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(dset);
      if (m_att) {
        dset = H5Dcreate(grp2, "Qp", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset);
        dset = H5Dcreate(grp2, "Qs", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset);
      }
      H5Sclose(dspace);
      H5Gclose(grp2);
    }

    H5Gclose(grp);

    // Top interface (topo)
    H5Sclose(attr_space1);
    H5Sclose(attr_space2);
    H5Sclose(attr_space3);
    H5Pclose(fapl);
    H5Fclose(h5_fid);
    m_isCreated = true;
  }

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_alignment(fapl, alignment, alignment);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  h5_fid = H5Fopen(fname, H5F_ACC_RDWR, fapl);
  if (h5_fid < 0) {
    cout << "Rank " << myid << " error opening file [" << fname << "]" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Write interfaces
  if (gridinfo == 1 && m_modestring == "Rho") {
    for (int g = mEW->mNumberOfCartesianGrids; g < ng; g++) {
      int real_g = g + 1;
      bool is_topo = (g == ng - 1);
      if (is_topo) real_g = g;

      size_t npts =
          ((size_t)(mGlobalDims[real_g][1] - mGlobalDims[real_g][0]) / stH +
           1) *
          ((mGlobalDims[real_g][3] - mGlobalDims[real_g][2]) / stH + 1);

      int nj = (int)(mWindow[real_g][3] - mWindow[real_g][2]) / stH + 1;
      int nk = mWindow[real_g][5];
      float* zfp = new float[npts];
#pragma omp parallel for
      for (int j = mWindow[real_g][2]; j <= mWindow[real_g][3]; j += stH)
        for (int i = mWindow[real_g][0]; i <= mWindow[real_g][1]; i += stH) {
          size_t ind = (size_t)(j - mWindow[real_g][2]) / stH +
                       nj * (i - mWindow[real_g][0]) / stH;
          /* ASSERT(ind < npts); */
          if (is_topo) zfp[ind] = (float)-mEW->mTopo(i, j, 1);
          /* zfp[ind] = (float) a_Z[real_g](i,j,1); */
          else
            zfp[ind] = (float)a_Z[real_g](i, j, nk);
        }

      sprintf(gname, "Z_interfaces");
      grp = H5Gopen(h5_fid, gname, H5P_DEFAULT);
      if (grp < 0) {
        cout << "Rank " << myid << " error opening [" << gname
             << "] group from file [" << fname << "]" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      sprintf(dname, "z_values_%d", ng - g - 1);
      dset = H5Dopen(grp, dname, H5P_DEFAULT);
      if (dset < 0) {
        cout << "Rank " << myid << " error opening [" << dname
             << "] dset from file [" << fname << "]" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      offsets[0] = mWindow[real_g][0] / stH == 0
                       ? 0
                       : ceil((double)mWindow[real_g][0] / stH) - 1;
      offsets[1] = mWindow[real_g][2] / stH == 0
                       ? 0
                       : ceil((double)mWindow[real_g][2] / stH) - 1;
      counts[0] = (hsize_t)(mWindow[real_g][1] - mWindow[real_g][0]) / stH + 1;
      counts[1] = (hsize_t)(mWindow[real_g][3] - mWindow[real_g][2]) / stH + 1;

      filespace = H5Dget_space(dset);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, counts,
                          NULL);

      memspace = H5Screate_simple(2, counts, NULL);

      ret = H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl, zfp);
      if (ret < 0) {
        cout << "Sfileoutput error writing interface!" << endl;
        cout << "Rank " << myid << ": offsets " << offsets[0] << ", "
             << offsets[1] << ", " << offsets[2] << endl;
        cout << "Rank " << myid << ": counts  " << counts[0] << ", "
             << counts[1] << ", " << counts[2] << endl;
        cout << "Rank " << myid << ": mWindow " << mWindow[real_g][0] << ", "
             << mWindow[real_g][1] << ", " << mWindow[real_g][2] << ", "
             << mWindow[real_g][3] << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      H5Sclose(memspace);
      H5Sclose(filespace);
      H5Dclose(dset);
      H5Gclose(grp);

      delete[] zfp;
    }  // end for g (curvilinear)

  }  // end if grid info

  for (int g = 0; g < ng; g++) {
    // Sfile grid order is reverse of sw4 grid order
    sprintf(gname, "/Material_model/grid_%d", ng - g - 1);

    grp = H5Gopen(h5_fid, gname, H5P_DEFAULT);
    if (grp < 0) {
      cout << "Rank " << myid << " error opening [" << gname
           << "] group from file [" << fname << "]" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    dset = H5Dopen(grp, m_modestring.c_str(), H5P_DEFAULT);
    if (dset < 0) {
      cout << "Rank " << myid << " error opening [" << m_modestring
           << "] dset from file [" << fname << "]" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    offsets[0] =
        mWindow[g][0] / stH == 0 ? 0 : ceil((double)mWindow[g][0] / stH) - 1;
    offsets[1] =
        mWindow[g][2] / stH == 0 ? 0 : ceil((double)mWindow[g][2] / stH) - 1;
    offsets[2] =
        mWindow[g][4] / stV == 0 ? 0 : ceil((double)mWindow[g][4] / stV) - 1;

    counts[0] = (hsize_t)(mWindow[g][1] - mWindow[g][0]) / stH + 1;
    counts[1] = (hsize_t)(mWindow[g][3] - mWindow[g][2]) / stH + 1;
    counts[2] = (hsize_t)(mWindow[g][5] - mWindow[g][4]) / stV + 1;

    filespace = H5Dget_space(dset);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, counts, NULL);

    memspace = H5Screate_simple(3, counts, NULL);

    /* if (g == 2) { */
    /*   cout << "Rank " << myid << ": offsets " << offsets[0]<< ", "
     * <<offsets[1]<< ", " <<offsets[2] << endl; */
    /*   cout << "Rank " << myid << ": counts  " << counts[0]<< ", "
     * <<counts[1]<< ", " <<counts[2] << endl; */
    /*   printf("Rank %d: grid %d, mWindow (%d, %d, %d, %d, %d, %d)\n", myid, g,
     * mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
     * mWindow[g][4], mWindow[g][5]); */
    /* } */

    if (m_double)
      ret = H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl,
                     m_doubleField[g]);
    else
      ret = H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl,
                     m_floatField[g]);

    if (ret < 0) {
      cout << "Sfileoutput error writing " << m_modestring << " dataset!"
           << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset);
    H5Gclose(grp);
  }

  H5Pclose(fapl);
  H5Pclose(dxpl);
  H5Fclose(h5_fid);
#else
  cout << "ERROR: cannot write sfile without sw4 compiled with HDF5 library!"
       << endl;
#endif
}
