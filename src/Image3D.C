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
#include <unistd.h>

#include <cstring>
#include <ctime>

#include "EW.h"
#include "Image3D.h"
#include "Require.h"
#include "caliper.h"
#include "mpi.h"
// static variable definition (in class only declaration):
int Image3D::mPreceedZeros = 0;

Image3D* Image3D::nil = static_cast<Image3D*>(0);

//-----------------------------------------------------------------------
Image3D::Image3D(EW* a_ew, float_sw4 time, float_sw4 timeInterval, int cycle,
                 int cycleInterval, float_sw4 tstart,
                 const std::string& filePrefix, Image3DMode mode,
                 bool doubleMode)
    : mTime(time),
      mEW(a_ew),
      m_time_done(false),
      mTimeInterval(timeInterval),
      mWritingCycle(cycle),
      mCycleInterval(cycleInterval),
      mFilePrefix(filePrefix),
      mImageSamplingFactor(1),
      mMode(mode),
      mFileName(""),
      //   mCycle(-1),
      mNextTime(0.0),
      mStartTime(tstart),
      m_isDefinedMPIWriters(false),
      m_double(doubleMode),
      m_winallocated(false),
      m_memallocated(false) {
  // note that the require2 macro doesn't generate any code unless compiled with
  // debugging options
  REQUIRE2(mode == UX || mode == UY || mode == UZ || mode == RHO || mode == P ||
               mode == S || mode == MU || mode == LAMBDA || mode == GRADRHO ||
               mode == GRADMU || mode == GRADLAMBDA || mode == GRADP ||
               mode == GRADS || mode == QP || mode == QS,
           "Image3D: mode=" << mode << " not supported");
  // the zcoord mode can only be save by an explicit call to the Image3D
  // constructor, and NOT from a volimage line in the command file
  if (mode == UX)
    m_modestring = "ux";
  else if (mode == UY)
    m_modestring = "uy";
  else if (mode == UZ)
    m_modestring = "uz";
  else if (mode == RHO)
    m_modestring = "rho";
  else if (mode == P)
    m_modestring = "p";
  else if (mode == S)
    m_modestring = "s";
  else if (mode == MU)
    m_modestring = "mu";
  else if (mode == LAMBDA)
    m_modestring = "lambda";
  else if (mode == GRADRHO)
    m_modestring = "gradrho";
  else if (mode == GRADMU)
    m_modestring = "gradmu";
  else if (mode == GRADLAMBDA)
    m_modestring = "gradlambda";
  else if (mode == GRADP)
    m_modestring = "gradp";
  else if (mode == GRADS)
    m_modestring = "grads";
  else if (mode == QP)
    m_modestring = "qp";
  else if (mode == QS)
    m_modestring = "qs";
}

//-----------------------------------------------------------------------
Image3D::~Image3D() {
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
// void Image3D::set_start_time(double tStart)
//{
//   mStartTime = tStart;
//}

//-----------------------------------------------------------------------
void Image3D::setup_images() {
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
    if (mImageSamplingFactor > 1 && m_ihavearray[g]) {
      for (int dim = 0; dim < 3; dim++) {
        int remainder = (mWindow[g][2 * dim] - mGlobalDims[g][2 * dim]) %
                        mImageSamplingFactor;
        if (remainder != 0)
          mWindow[g][2 * dim] += mImageSamplingFactor - remainder;
        remainder = (mWindow[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) %
                    mImageSamplingFactor;
        if (remainder != 0) mWindow[g][2 * dim + 1] -= remainder;

        remainder = (mGlobalDims[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) %
                    mImageSamplingFactor;
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
            ((size_t)(mWindow[g][1] - mWindow[g][0]) / mImageSamplingFactor +
             1) *
            ((mWindow[g][3] - mWindow[g][2]) / mImageSamplingFactor + 1) *
            ((mWindow[g][5] - mWindow[g][4]) / mImageSamplingFactor + 1 +
             m_extraz[g]);
        m_doubleField[g] = new double[npts];
      } else
        m_doubleField[g] = new double[1];
    }
  } else {
    m_floatField.resize(mEW->mNumberOfGrids);
    for (int g = 0; g < mEW->mNumberOfGrids; g++) {
      if (m_ihavearray[g]) {
        size_t npts =
            ((size_t)(mWindow[g][1] - mWindow[g][0]) / mImageSamplingFactor +
             1) *
            ((mWindow[g][3] - mWindow[g][2]) / mImageSamplingFactor + 1) *
            ((mWindow[g][5] - mWindow[g][4]) / mImageSamplingFactor + 1 +
             m_extraz[g]);
        m_floatField[g] = new float[npts];
      } else
        m_floatField[g] = new float[1];
    }
  }
  m_memallocated = true;
  define_pio();
}

//-----------------------------------------------------------------------
void Image3D::define_pio() {
  int glow = 0, ghigh = mEW->mNumberOfGrids;
  m_parallel_io = new Parallel_IO*[ghigh - glow + 1];
  for (int g = glow; g < ghigh; g++) {
    int global[3], local[3], start[3];
    for (int dim = 0; dim < 3; dim++) {
      global[dim] = (mGlobalDims[g][2 * dim + 1] - mGlobalDims[g][2 * dim]) /
                        mImageSamplingFactor +
                    1;
      local[dim] = (mWindow[g][2 * dim + 1] - mWindow[g][2 * dim]) /
                       mImageSamplingFactor +
                   1;
      start[dim] = (mWindow[g][2 * dim] - mGlobalDims[g][2 * dim]) /
                   mImageSamplingFactor;
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
void Image3D::setSteps(int a_steps) {
  char buffer[50];
  mPreceedZeros = snprintf(buffer, 50, "%d", a_steps);
}

//-----------------------------------------------------------------------
// void Image3D::set_double(bool val)
//{
//   m_double = val;
//}

//-----------------------------------------------------------------------
bool Image3D::timeToWrite(float_sw4 time, int cycle, float_sw4 dt) {
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
void Image3D::update_image(int a_cycle, float_sw4 a_time, float_sw4 a_dt,
                           vector<Sarray>& a_U, vector<Sarray>& a_Rho,
                           vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                           vector<Sarray>& a_gRho, vector<Sarray>& a_gMu,
                           vector<Sarray>& a_gLambda, vector<Sarray>& a_Qp,
                           vector<Sarray>& a_Qs, std::string a_path,
                           std::vector<Sarray>& a_Z) {
  SW4_MARK_FUNCTION;
  if (timeToWrite(a_time, a_cycle, a_dt)) {
    compute_image(a_U, a_Rho, a_Mu, a_Lambda, a_gRho, a_gMu, a_gLambda, a_Qp,
                  a_Qs);
    write_image(a_cycle, a_path, a_time, a_Z);
  }
}

//-----------------------------------------------------------------------
void Image3D::force_write_image(float_sw4 a_time, int a_cycle,
                                vector<Sarray>& a_U, vector<Sarray>& a_Rho,
                                vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                                vector<Sarray>& a_gRho, vector<Sarray>& a_gMu,
                                vector<Sarray>& a_gLambda, vector<Sarray>& a_Qp,
                                vector<Sarray>& a_Qs, std::string a_path,
                                std::vector<Sarray>& a_Z) {
  compute_image(a_U, a_Rho, a_Mu, a_Lambda, a_gRho, a_gMu, a_gLambda, a_Qp,
                a_Qs);
  write_image(a_cycle, a_path, a_time, a_Z);
}

//-----------------------------------------------------------------------
void Image3D::compute_image(vector<Sarray>& a_U, vector<Sarray>& a_Rho,
                            vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                            vector<Sarray>& a_gRho, vector<Sarray>& a_gMu,
                            vector<Sarray>& a_gLambda, vector<Sarray>& a_Qp,
                            vector<Sarray>& a_Qs) {
  SW4_MARK_FUNCTION;
  // Introduce 'st' to simplify the variable name
  int st = mImageSamplingFactor;
  for (int g = 0; g < mEW->mNumberOfGrids; g++) {
    int il = mEW->m_iStart[g];
    int iu = mEW->m_iEnd[g];
    int jl = mEW->m_jStart[g];
    int ju = mEW->m_jEnd[g];
    int kl = mEW->m_kStart[g];
    int ku = mEW->m_kEnd[g];
    int ni = (iu - il + 1);
    int nj = (ju - jl + 1);

    float_sw4* up = a_U[g].c_ptr();
    if (mMode == RHO)
      up = a_Rho[g].c_ptr();
    else if (mMode == MU)
      up = a_Mu[g].c_ptr();
    else if (mMode == LAMBDA)
      up = a_Lambda[g].c_ptr();
    else if (mMode == GRADRHO)
      up = a_gRho[g].c_ptr();
    else if (mMode == GRADMU)
      up = a_gMu[g].c_ptr();
    else if (mMode == GRADLAMBDA)
      up = a_gLambda[g].c_ptr();
    else if (mMode == QP)
      up = a_Qp[g].c_ptr();
    else if (mMode == QS)
      up = a_Qs[g].c_ptr();

    int niw = (mWindow[g][1] - mWindow[g][0]) / st + 1;
    int nijw = ni * ((mWindow[g][3] - mWindow[g][2]) / st + 1);
    if (mMode == UX || mMode == UY || mMode == UZ) {
      int c = 0;
      if (mMode == UY) c = 1;
      if (mMode == UZ) c = 2;
      if (m_double) {
        //	    for( int ks=0 ; ks <= (mWindow[g][5]-mWindow[g][4])/st ;
        // ks++ ) 	       for( int js=0 ; js <=
        // (mWindow[g][3]-mWindow[g][2])/st ; js++ ) 		  for( int is=0
        // ; is <=
        //(mWindow[g][1]-mWindow[g][0])/st ; is++ )
        //		  {
        //		     int k = mWindow[g][4]+ks*st;
        //		     int j = mWindow[g][2]+js*st;
        //		     int i = mWindow[g][0]+is*st;
        //		     size_t ind = is+ni*js+nij*ks;
        //		     m_doubleField[g][ind]= a_U(c,i,j,k);
        //		  }
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_doubleField[g][ind] = (double)a_U[g](c + 1, i, j, k);
              //		     m_doubleField[g][ind] = up[c+
              // 3*((i-il)+ni*(j-jl)+((size_t)ni)*nj*(k-kl))];
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_floatField[g][ind] = (float)a_U[g](c + 1, i, j, k);
              //		     m_floatField[g][ind++] = up[c+
              // 3*((i-il)+ni*(j-jl)+((size_t)ni)*nj*(k-kl))];
            }
      }
    } else if (mMode == RHO || mMode == MU || mMode == LAMBDA ||
               mMode == GRADRHO || mMode == GRADMU || mMode == GRADLAMBDA ||
               mMode == QP || mMode == QS)  // these modes just copy the values
                                            // straight from the up array
    {
      if (m_double) {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_doubleField[g][ind] =
                  up[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)];
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_floatField[g][ind] =
                  up[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)];
            }
      }
    } else if (mMode == P) {
      float_sw4* rho = a_Rho[g].c_ptr();
      float_sw4* mu = a_Mu[g].c_ptr();
      float_sw4* la = a_Lambda[g].c_ptr();
      if (m_double) {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_doubleField[g][ind] = sqrt(
                  (2 * mu[(i - il) + ni * (j - jl) +
                          ((size_t)ni) * nj * (k - kl)] +
                   la[(i - il) + ni * (j - jl) +
                      ((size_t)ni) * nj * (k - kl)]) /
                  rho[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)]);
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_floatField[g][ind] = sqrt(
                  (2 * mu[(i - il) + ni * (j - jl) +
                          ((size_t)ni) * nj * (k - kl)] +
                   la[(i - il) + ni * (j - jl) +
                      ((size_t)ni) * nj * (k - kl)]) /
                  rho[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)]);
            }
      }
    } else if (mMode == S) {
      float_sw4* rho = a_Rho[g].c_ptr();
      float_sw4* mu = a_Mu[g].c_ptr();
      if (m_double) {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_doubleField[g][ind] = sqrt(
                  mu[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)] /
                  rho[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)]);
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_floatField[g][ind] = sqrt(
                  mu[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)] /
                  rho[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)]);
            }
      }
    } else if (mMode == GRADP) {
      float_sw4* rho = a_Rho[g].c_ptr();
      float_sw4* mu = a_Mu[g].c_ptr();
      float_sw4* la = a_Lambda[g].c_ptr();
      float_sw4* gla = a_gLambda[g].c_ptr();

      if (m_double) {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_doubleField[g][ind] =
                  2 *
                  sqrt((2 * mu[(i - il) + ni * (j - jl) +
                               ((size_t)ni) * nj * (k - kl)] +
                        la[(i - il) + ni * (j - jl) +
                           ((size_t)ni) * nj * (k - kl)]) *
                       rho[(i - il) + ni * (j - jl) +
                           ((size_t)ni) * nj * (k - kl)]) *
                  gla[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)];
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_floatField[g][ind] =
                  2 *
                  sqrt((2 * mu[(i - il) + ni * (j - jl) +
                               ((size_t)ni) * nj * (k - kl)] +
                        la[(i - il) + ni * (j - jl) +
                           ((size_t)ni) * nj * (k - kl)]) *
                       rho[(i - il) + ni * (j - jl) +
                           ((size_t)ni) * nj * (k - kl)]) *
                  gla[(i - il) + ni * (j - jl) + ((size_t)ni) * nj * (k - kl)];
            }
      }
    } else if (mMode == GRADS) {
      float_sw4* rho = a_Rho[g].c_ptr();
      float_sw4* mu = a_Mu[g].c_ptr();
      float_sw4* gla = a_gLambda[g].c_ptr();
      float_sw4* gmu = a_gMu[g].c_ptr();

      if (m_double) {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_doubleField[g][ind] = sqrt(mu[(i - il) + ni * (j - jl) +
                                              ((size_t)ni) * nj * (k - kl)] *
                                           rho[(i - il) + ni * (j - jl) +
                                               ((size_t)ni) * nj * (k - kl)]) *
                                      (2 * gmu[(i - il) + ni * (j - jl) +
                                               ((size_t)ni) * nj * (k - kl)] -
                                       4 * gla[(i - il) + ni * (j - jl) +
                                               ((size_t)ni) * nj * (k - kl)]);
            }
      } else {
#pragma omp parallel for
        for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
          for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
            for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
              size_t ind = (i - mWindow[g][0]) / st +
                           niw * (j - mWindow[g][2]) / st +
                           nijw * (k - mWindow[g][4]) / st;
              m_floatField[g][ind] = sqrt(mu[(i - il) + ni * (j - jl) +
                                             ((size_t)ni) * nj * (k - kl)] *
                                          rho[(i - il) + ni * (j - jl) +
                                              ((size_t)ni) * nj * (k - kl)]) *
                                     (2 * gmu[(i - il) + ni * (j - jl) +
                                              ((size_t)ni) * nj * (k - kl)] -
                                      4 * gla[(i - il) + ni * (j - jl) +
                                              ((size_t)ni) * nj * (k - kl)]);
            }
      }
    }
    // Extrapolate extra point in z
    if (m_extraz[g] == 1) {
      int k = mWindow[g][5] + st;
      size_t koff = ((size_t)(mWindow[g][1] - mWindow[g][0]) / st + 1) *
                    ((mWindow[g][3] - mWindow[g][2]) / st + 1);
      size_t ind = (k - mWindow[g][4]) / st * koff;

      // Linear extrapolation assumes that
      // mWindow[g][5]-st>=mWindow[g][4], i.e. mWindow[g][5] -mWindow[g][4]>=st.
      int ok = 2;
      if (mWindow[g][5] - mWindow[g][4] < st) ok = 1;

      if (m_double) {
        for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
          for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
            m_doubleField[g][ind] = 2 * m_doubleField[g][ind - koff] -
                                    m_doubleField[g][ind - ok * koff];
            ind++;
          }

      } else {
        for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
          for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
            m_floatField[g][ind] = 2 * m_floatField[g][ind - koff] -
                                   m_floatField[g][ind - ok * koff];
            ind++;
          }
      }
    }
  }
}

//-----------------------------------------------------------------------
void Image3D::compute_file_suffix(int cycle, std::stringstream& fileSuffix) {
  fileSuffix << mFilePrefix << ".cycle=";
  int temp = static_cast<int>(pow(10.0, mPreceedZeros - 1));
  int testcycle = cycle;
  if (cycle == 0) testcycle = 1;
  while (testcycle < temp) {
    fileSuffix << "0";
    temp /= 10;
  }
  // AP changed the suffix to 3Dimg
  fileSuffix << cycle;
  fileSuffix << "." << m_modestring << ".3Dimg";
}

//-----------------------------------------------------------------------
// void Image3D::write_image(int cycle, std::string& path, float_sw4 t,
//                           Sarray& a_Z) {
//   SW4_MARK_FUNCTION;
//   // File format:
//   //
//   // [precision(int), npatches(int), time(double), plane(int=-1),
//   // coordinate(double=-1),
//   //  imagetype(int), gridinfo(int), timeofday(string),
//   //     h_1(double)  zmin_1(double)  sizes_1(6 int)
//   //     h_2(double)  zmin_2(double)  sizes_2(6 int)
//   //                    ....
//   //     h_ng(double) zmin_ng(double) sizes_ng(6 int)
//   //     data_1(float/double array) ... data_ng(float/double array)
//   //     [z(float/double array)] ]
//   //  plane and coordinate are always -1, exist here for compatibility with the
//   //  2D images imagetype = mode nr. (ux=1,uy=2, etc..) gridinfo  = -1 no grid
//   //  info given
//   //               0 all blocks are on Cartesian grids, no z-coordinate needed
//   //               1 Curvilinear grid is stored after the data blocks
//   // timeofday - String describing the creation date of the file, max 25 bytes.
//   //

//   int ng = mEW->mNumberOfGrids;
//   // offset initialized to header size:
//   off_t offset = 2 * sizeof(int) + 2 * sizeof(double) + 3 * sizeof(int) +
//                  25 * sizeof(char) +
//                  ng * (2 * sizeof(double) + 6 * sizeof(int));

//   ASSERT(m_isDefinedMPIWriters);
//   int gridinfo = 0;
//   if (mEW->topographyExists()) gridinfo = 1;

//   bool iwrite = false;
//   for (int g = 0; g < ng; g++) iwrite = iwrite || m_parallel_io[g]->i_write();

//   int fid = -1;
//   std::stringstream s, fileSuffix;

//   if (iwrite) {
//     compute_file_suffix(cycle, fileSuffix);
//     if (path != ".") s << path;
//     s << fileSuffix.str();  // string 's' is the file name including path
//   }

//   int st = mImageSamplingFactor;

//   // Open file from processor zero and write header.
//   if (m_parallel_io[0]->proc_zero()) {
//     fid = open(const_cast<char*>(s.str().c_str()), O_CREAT | O_TRUNC | O_WRONLY,
//                0660);
//     CHECK_INPUT(fid != -1, "Image3D::write_image: Error opening: " << s.str());
//     int myid;

//     MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//     std::cout << "writing volume image on file " << s.str();
//     std::cout << " (msg from proc # " << myid << ")" << std::endl;

//     int prec = m_double ? 8 : 4;
//     size_t ret = write(fid, &prec, sizeof(int));
//     CHECK_INPUT(ret == sizeof(int),
//                 "Image3D::write_image: Error writing precision");
//     ret = write(fid, &ng, sizeof(int));
//     CHECK_INPUT(ret == sizeof(int), "Image3D::write_image: Error writing ng");
//     double dblevar = (double)t;
//     ret = write(fid, &dblevar, sizeof(double));
//     CHECK_INPUT(ret == sizeof(double),
//                 "Image3D::write_image: Error writing time");
//     int mone = -1;
//     ret = write(fid, &mone, sizeof(int));
//     CHECK_INPUT(ret == sizeof(int), "Image3D::write_image: Error writing -1");
//     double dum = -1;
//     ret = write(fid, &dum, sizeof(double));
//     CHECK_INPUT(ret == sizeof(double),
//                 "Image3D::write_image: Error writing dummy-coord");

//     int imode = static_cast<int>(mMode);
//     ret = write(fid, &imode, sizeof(int));
//     CHECK_INPUT(ret == sizeof(int),
//                 "Image3D::write_image could not write imode");

//     ret = write(fid, &gridinfo, sizeof(int));
//     CHECK_INPUT(ret == sizeof(int),
//                 "Image3D::write_image could not write gridinfo");

//     time_t realtime;
//     time(&realtime);
//     string strtime;
//     strtime += asctime(localtime(&realtime));
//     char strtimec[25];

//     strncpy(strtimec, strtime.c_str(), 25);
//     ret = write(fid, strtimec, 25 * sizeof(char));
//     CHECK_INPUT(ret == 25 * sizeof(char),
//                 "Image3D::write_image could not write strtimec");

//     for (int g = 0; g < ng; g++) {
//       double h = (double)mEW->mGridSize[g] * mImageSamplingFactor;
//       ret = write(fid, &h, sizeof(double));
//       CHECK_INPUT(ret == sizeof(double),
//                   "Image3D::write_image: Error writing h");
//       dblevar = (double)mEW->m_zmin[g];
//       ret = write(fid, &dblevar, sizeof(double));
//       CHECK_INPUT(ret == sizeof(double),
//                   "Image3D::write_image: Error writing zmin");
//       int globalSize[6];
//       globalSize[0] = 1;
//       globalSize[1] = (mGlobalDims[g][1] - mGlobalDims[g][0]) / st + 1;
//       globalSize[2] = 1;
//       globalSize[3] = (mGlobalDims[g][3] - mGlobalDims[g][2]) / st + 1;
//       globalSize[4] = 1;
//       globalSize[5] =
//           (mGlobalDims[g][5] - mGlobalDims[g][4]) / st + 1 + m_extraz[g];
//       ret = write(fid, globalSize, 6 * sizeof(int));
//       CHECK_INPUT(ret == 6 * sizeof(int),
//                   "Image3D::write_image: Error writing global sizes");
//     }
//     fsync(fid);
//   }
//   m_parallel_io[0]->writer_barrier();

//   // Open file from all writers
//   if (iwrite && !m_parallel_io[0]->proc_zero()) {
//     fid = open(const_cast<char*>(s.str().c_str()), O_WRONLY);
//     CHECK_INPUT(fid != -1,
//                 "Image3D::write_images:: Error opening: " << s.str());
//   }

//   // Write data blocks
//   for (int g = 0; g < ng; g++) {
//     size_t npts =
//         ((size_t)(mGlobalDims[g][1] - mGlobalDims[g][0]) / st + 1) *
//         ((mGlobalDims[g][3] - mGlobalDims[g][2]) / st + 1) *
//         ((mGlobalDims[g][5] - mGlobalDims[g][4]) / st + 1 + m_extraz[g]);

//     if (!mEW->usingParallelFS() || g == 0) m_parallel_io[g]->writer_barrier();

//     if (m_double) {
//       char cprec[] = "double";
//       m_parallel_io[g]->write_array(&fid, 1, m_doubleField[g], offset, cprec);
//       offset += npts * sizeof(double);
//     } else {
//       char cprec[] = "float";
//       m_parallel_io[g]->write_array(&fid, 1, m_floatField[g], offset, cprec);
//       offset += npts * sizeof(float);
//     }
//   }

//   // Add curvilinear grid, if needed
//   if (gridinfo == 1) {
//     int g = ng - 1;
//     float_sw4* zp = a_Z.c_ptr();
//     size_t npts =
//         ((size_t)(mGlobalDims[g][1] - mGlobalDims[g][0]) / st + 1) *
//         ((mGlobalDims[g][3] - mGlobalDims[g][2]) / st + 1) *
//         ((mGlobalDims[g][5] - mGlobalDims[g][4]) / st + 1 + m_extraz[g]);

//     if (!mEW->usingParallelFS() || g == 0) m_parallel_io[g]->writer_barrier();

//     size_t nptsloc =
//         ((size_t)(mWindow[g][1] - mWindow[g][0]) / mImageSamplingFactor + 1) *
//         ((mWindow[g][3] - mWindow[g][2]) / mImageSamplingFactor + 1) *
//         ((mWindow[g][5] - mWindow[g][4]) / mImageSamplingFactor + 1 +
//          m_extraz[g]);

//     int ni = (mWindow[g][1] - mWindow[g][0]) / st + 1;
//     int nij = ni * ((mWindow[g][3] - mWindow[g][2]) / st + 1);
//     if (m_double) {
//       double* zfp = new double[nptsloc];
// #pragma omp parallel for
//       for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
//         for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
//           for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
//             size_t ind = (i - mWindow[g][0]) / st +
//                          ni * (j - mWindow[g][2]) / st +
//                          nij * (k - mWindow[g][4]) / st;
//             zfp[ind] = (double)a_Z(i, j, k);
//           }
//       //	 for( size_t i = 0; i < nptsloc ; i++ )
//       //	    zfp[i] = zp[i];
//       char cprec[] = "double";
//       m_parallel_io[g]->write_array(&fid, 1, zfp, offset, cprec);
//       offset += npts * sizeof(double);
//       delete[] zfp;
//     } else {
//       float* zfp = new float[nptsloc];
// #pragma omp parallel for
//       for (int k = mWindow[g][4]; k <= mWindow[g][5]; k += st)
//         for (int j = mWindow[g][2]; j <= mWindow[g][3]; j += st)
//           for (int i = mWindow[g][0]; i <= mWindow[g][1]; i += st) {
//             size_t ind = (i - mWindow[g][0]) / st +
//                          ni * (j - mWindow[g][2]) / st +
//                          nij * (k - mWindow[g][4]) / st;
//             zfp[ind] = (float)a_Z(i, j, k);
//           }
//       //	 for( size_t i = 0; i < nptsloc ; i++ )
//       //	    zfp[i] = zp[i];
//       char cprec[] = "float";
//       m_parallel_io[g]->write_array(&fid, 1, zfp, offset, cprec);
//       offset += npts * sizeof(float);
//       delete[] zfp;
//     }
//   }
//   if (iwrite) close(fid);
// }

//-----------------------------------------------------------------------
void EW::read_volimage(std::string& path, std::string& fname,
                       vector<Sarray>& data) {
  // File format:
  //
  // [precision(int), npatches(int), time(double), plane(int=-1),
  // coordinate(double=-1),
  //  imagetype(int), gridinfo(int), timeofday(string),
  //     h_1(double)  zmin_1(double)  sizes_1(6 int)
  //     h_2(double)  zmin_2(double)  sizes_2(6 int)
  //                    ....
  //     h_ng(double) zmin_ng(double) sizes_ng(6 int)
  //     data_1(float/double array) ... data_ng(float/double array)
  //     [z(float/double array)] ]
  //  plane and coordinate are always -1, exist here for compatibility with the
  //  2D images imagetype = mode nr. (ux=1,uy=2, etc..) gridinfo  = -1 no grid
  //  info given
  //               0 all blocks are on Cartesian grids, no z-coordinate needed
  //               1 Curvilinear grid is stored after the data blocks
  // timeofday - String describing the creation date of the file, max 25 bytes.
  //

  vector<Parallel_IO*> parallel_io;
  define_parallel_io(parallel_io);

  int ng = mNumberOfGrids;
  // offset initialized to header size:
  off_t offset = 2 * sizeof(int) + 2 * sizeof(double) + 3 * sizeof(int) +
                 25 * sizeof(char) +
                 ng * (2 * sizeof(double) + 6 * sizeof(int));

  bool iread = false;
  for (int g = 0; g < ng; g++) iread = iread || parallel_io[g]->i_write();

  int fid = -1;
  std::stringstream s;

  if (iread) {
    if (path != ".")
      s << path << "/" << fname;
    else
      s << fname;
  }
  // Open file from processor zero and read header.

  int prec = 0;
  if (parallel_io[0]->proc_zero()) {
    fid = open(const_cast<char*>(s.str().c_str()), O_RDONLY);
    CHECK_INPUT(fid != -1, "EW::read_image: Error opening: " << s.str());
    int myid;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "reading volume image on file " << s.str();
    std::cout << " (msg from proc # " << myid << ")" << std::endl;

    size_t ret = read(fid, &prec, sizeof(int));
    CHECK_INPUT(ret == sizeof(int), "EW::read_image: Error reading precision");
    int ngfile;
    ret = read(fid, &ngfile, sizeof(int));
    CHECK_INPUT(ret == sizeof(int), "EW::read_image: Error reading ng");
    CHECK_INPUT(ng == ngfile,
                "EW::read_image: Error number of grids on file "
                    << "does not match number of grids in computation");
    double t;
    ret = read(fid, &t, sizeof(double));
    CHECK_INPUT(ret == sizeof(double), "EW::read_image: Error reading time");
    int mone;
    ret = read(fid, &mone, sizeof(int));
    CHECK_INPUT(ret == sizeof(int), "EW::read_image: Error reading -1");
    CHECK_INPUT(mone == -1, "EW::read_image: Error reading -1 != " << mone);

    double dum;
    ret = read(fid, &dum, sizeof(double));
    CHECK_INPUT(ret == sizeof(double),
                "EW::read_image: Error reading dummy-coord");

    int imode;
    ret = read(fid, &imode, sizeof(int));
    CHECK_INPUT(ret == sizeof(int), "EW::read_image could not read imode");

    int gridinfo;
    ret = read(fid, &gridinfo, sizeof(int));
    CHECK_INPUT(ret == sizeof(int), "EW::read_image could not read gridinfo");

    char strtimec[25];
    ret = read(fid, strtimec, 25 * sizeof(char));
    CHECK_INPUT(ret == 25 * sizeof(char),
                "EW::read_image could not read strtimec");

    for (int g = 0; g < ng; g++) {
      double tol = 1e-12;
      if (sizeof(float_sw4) < 8) tol = 1e-6;

      double h;
      // = mEW->mGridSize[g]*mImageSamplingFactor;
      ret = read(fid, &h, sizeof(double));
      CHECK_INPUT(ret == sizeof(double), "EW::read_image: Error reading h");
      CHECK_INPUT(abs(h - mGridSize[g]) < tol,
                  "EW::read_image: Error h on file does not"
                      << " match h in computation"
                      << " hfile = " << h << " hcomp= " << mGridSize[g]);
      double zmin;
      ret = read(fid, &zmin, sizeof(double));
      CHECK_INPUT(ret == sizeof(double), "EW::read_image: Error reading zmin");
      CHECK_INPUT(abs(zmin - m_zmin[g]) < tol,
                  "EW::read_image: Error zmin on file does not"
                      << " match zmin in computation"
                      << " zmin file = " << zmin
                      << " zmin comp= " << m_zmin[g]);
      int globalSize[6];
      ret = read(fid, globalSize, 6 * sizeof(int));
      CHECK_INPUT(ret == 6 * sizeof(int),
                  "EW::read_image: Error reading global sizes");
      CHECK_INPUT(globalSize[0] == 1, "EW::read_image: Error in global sizes, "
                                          << "low i-index is "
                                          << globalSize[0]);
      CHECK_INPUT(globalSize[1] == m_global_nx[g],
                  "EW::read_image: Error in global sizes, "
                      << "upper i-index is " << globalSize[1]);
      CHECK_INPUT(globalSize[2] == 1, "EW::read_image: Error in global sizes, "
                                          << "low j-index is "
                                          << globalSize[2]);
      CHECK_INPUT(globalSize[3] == m_global_ny[g],
                  "EW::read_image: Error in global sizes, "
                      << "upper j-index is " << globalSize[3]);
      CHECK_INPUT(globalSize[4] == 1, "EW::read_image: Error in global sizes, "
                                          << "low k-index is "
                                          << globalSize[4]);
      CHECK_INPUT(globalSize[5] == m_global_nz[g],
                  "EW::read_image: Error in global sizes, "
                      << "upper k-index is " << globalSize[5]);
    }
  }

  parallel_io[0]->writer_barrier();
  int tmpprec = prec;
  MPI_Allreduce(&tmpprec, &prec, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // Open file from all readers
  if (iread && !parallel_io[0]->proc_zero()) {
    fid = open(const_cast<char*>(s.str().c_str()), O_RDONLY);
    CHECK_INPUT(fid != -1, "EW::read_image:: Error opening file: " << s.str());
  }

  // Read data blocks
  for (int g = 0; g < ng; g++) {
    size_t npts = (size_t)(m_global_nx[g]) * (size_t)(m_global_ny[g]) *
                  (size_t)(m_global_nz[g]);
    size_t nptsloc = (size_t)((m_iEndInt[g] - m_iStartInt[g] + 1)) *
                     (m_jEndInt[g] - m_jStartInt[g] + 1) *
                     (m_kEndInt[g] - m_kStartInt[g] + 1);

    if (!usingParallelFS() || g == 0) parallel_io[g]->writer_barrier();

    double* doubleField = new double[nptsloc];
    if (prec == 8) {
      parallel_io[g]->read_array(&fid, 1, doubleField, offset, "double");
      offset += npts * sizeof(double);
    } else {
      parallel_io[g]->read_array(&fid, 1, doubleField, offset, "float");
      offset += npts * sizeof(float);
    }
    data[g].insert_subarray(m_iStartInt[g], m_iEndInt[g], m_jStartInt[g],
                            m_jEndInt[g], m_kStartInt[g], m_kEndInt[g],
                            doubleField);
    delete[] doubleField;
  }
  if (iread) close(fid);
  for (int g = 0; g < mNumberOfGrids; g++) delete parallel_io[g];
}

//-----------------------------------------------------------------------
void EW::define_parallel_io(vector<Parallel_IO*>& parallel_io) {
  parallel_io.resize(mNumberOfGrids);
  for (int g = 0; g < mNumberOfGrids; g++) {
    int global[3], local[3], start[3];
    global[0] = m_global_nx[g];
    global[1] = m_global_ny[g];
    global[2] = m_global_nz[g];
    local[0] = m_iEndInt[g] - m_iStartInt[g] + 1;
    local[1] = m_jEndInt[g] - m_jStartInt[g] + 1;
    local[2] = m_kEndInt[g] - m_kStartInt[g] + 1;
    // Offset assume global array starts at (i,j,k)=(1,1,1)
    // and that the k-dimension is always in only one processor.
    start[0] = m_iStartInt[g] - 1;
    start[1] = m_jStartInt[g] - 1;
    start[2] = 0;

    int iwrite = 0;
    int nrwriters = getNumberOfWritersPFS();
    int nproc = 0, myid = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // Find out the I/O processors
    // Assume all processors have some part of the array.

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
    parallel_io[g] =
        new Parallel_IO(iwrite, usingParallelFS(), global, local, start);
  }
}
//-----------------------------------------------------------------------
void Image3D::write_image( int cycle, std::string &path, float_sw4 t,
			   std::vector<Sarray>& a_Z )
{
  //File format: 
  //
  // [precision(int), npatches(int), time(double), plane(int=-1), coordinate(double=-1), 
  //  imagetype(int), gridinfo(int), timeofday(string),
  //     h_1(double)  zmin_1(double)  sizes_1(6 int) 
  //     h_2(double)  zmin_2(double)  sizes_2(6 int) 
  //                    ....
  //     h_ng(double) zmin_ng(double) sizes_ng(6 int) 
  //     data_1(float/double array) ... data_ng(float/double array) [z(float/double array)] ]
  //  plane and coordinate are always -1, exist here for compatibility with the 2D images
  //  imagetype = mode nr. (ux=1,uy=2, etc..)
  //  gridinfo  = -1 no grid info given
  //               0 all blocks are on Cartesian grids, no z-coordinate needed
  //               1 Curvilinear grid is stored after the data blocks
  // timeofday - String describing the creation date of the file, max 25 bytes.
  //

   int ng = mEW->mNumberOfGrids;
  // offset initialized to header size:
   off_t offset = 2*sizeof(int) + 2*sizeof(double) + 3*sizeof(int) + 25*sizeof(char) +
      ng*(2*sizeof(double)+6*sizeof(int));

   ASSERT(m_isDefinedMPIWriters);
   int gridinfo = 0;
   if( mEW->topographyExists() )
      gridinfo = 1;

   bool iwrite = false;
   for( int g=0 ; g < ng ; g++ )
      iwrite = iwrite || m_parallel_io[g]->i_write();
  
   int fid=-1;
   std::stringstream s, fileSuffix;

   if( iwrite )
   {
      compute_file_suffix( cycle, fileSuffix );
      if( path != "." )
	 s << path;
      s << fileSuffix.str(); // string 's' is the file name including path
   }

   int st = mImageSamplingFactor;

   // Open file from processor zero and write header.
   if( m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
      CHECK_INPUT(fid != -1, "Image3D::write_image: Error opening: " << s.str() );
      int myid;

      MPI_Comm_rank( MPI_COMM_WORLD, &myid );
      std::cout << "writing volume image on file " << s.str();
      std::cout << " (msg from proc # " << myid << ")" << std::endl;

      int prec = m_double ? 8 : 4;
      size_t ret = write(fid,&prec,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"Image3D::write_image: Error writing precision" );
      ret = write(fid,&ng,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"Image3D::write_image: Error writing ng" );
      double dblevar = (double) t;
      ret = write(fid,&dblevar,sizeof(double));
      CHECK_INPUT( ret == sizeof(double),"Image3D::write_image: Error writing time" );
      int mone = -1;
      ret = write(fid,&mone,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"Image3D::write_image: Error writing -1" );
      double dum=-1;
      ret = write(fid,&dum,sizeof(double));
      CHECK_INPUT( ret == sizeof(double),"Image3D::write_image: Error writing dummy-coord" );

      int imode = static_cast<int>(mMode);
      ret = write(fid,&imode,sizeof(int));
      CHECK_INPUT( ret == sizeof(int), "Image3D::write_image could not write imode" );

      ret = write(fid,&gridinfo,sizeof(int));
      CHECK_INPUT( ret == sizeof(int), "Image3D::write_image could not write gridinfo");

      time_t realtime;
      time(&realtime);
      string strtime;
      strtime += asctime(localtime(&realtime));
      char strtimec[25];

      strncpy(strtimec,strtime.c_str(),25);
      ret = write(fid,strtimec,25*sizeof(char));
      CHECK_INPUT(ret == 25*sizeof(char),"Image3D::write_image could not write strtimec");

      for(int g = 0; g < ng ;g++ )
      {
         double h = (double) mEW->mGridSize[g]*mImageSamplingFactor;
	 ret = write( fid, &h, sizeof(double) );
	 CHECK_INPUT( ret == sizeof(double),"Image3D::write_image: Error writing h" );
	 dblevar = (double) mEW->m_zmin[g];
         ret = write( fid, &dblevar, sizeof(double) );
	 CHECK_INPUT( ret == sizeof(double),"Image3D::write_image: Error writing zmin" );
	 int globalSize[6];
         globalSize[0] = 1;
         globalSize[1] = (mGlobalDims[g][1]-mGlobalDims[g][0])/st+1;
         globalSize[2] = 1;
         globalSize[3] = (mGlobalDims[g][3]-mGlobalDims[g][2])/st+1;
         globalSize[4] = 1;
         globalSize[5] = (mGlobalDims[g][5]-mGlobalDims[g][4])/st+1+m_extraz[g];
	 ret = write( fid, globalSize, 6*sizeof(int) );
	 CHECK_INPUT( ret == 6*sizeof(int),"Image3D::write_image: Error writing global sizes" );
      }
      fsync(fid);
   }
   m_parallel_io[0]->writer_barrier();

   // Open file from all writers
   if( iwrite && !m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_WRONLY );
      CHECK_INPUT(fid != -1, "Image3D::write_images:: Error opening: " << s.str() );
   }

   // Write data blocks
   for( int g = 0 ; g < ng ; g++ )
   {
      size_t npts = ((size_t)(mGlobalDims[g][1]-mGlobalDims[g][0])/st+1)*
	 ((mGlobalDims[g][3]-mGlobalDims[g][2])/st+1)*
	 ((mGlobalDims[g][5]-mGlobalDims[g][4])/st+1+m_extraz[g]);

      if( !mEW->usingParallelFS() || g == 0 )
	 m_parallel_io[g]->writer_barrier();
      
      if( m_double )
      {
	char cprec[]="double";
	m_parallel_io[g]->write_array( &fid, 1, m_doubleField[g], offset, cprec );
	offset += npts*sizeof(double);
      }
      else
      {
	char cprec[]="float";
	m_parallel_io[g]->write_array( &fid, 1, m_floatField[g], offset, cprec );
	offset += npts*sizeof(float);
      }
   }

   // Add curvilinear grid, if needed
   if( gridinfo == 1 )
   {
      for (int g = mEW->mNumberOfCartesianGrids; g < mEW->mNumberOfGrids; g++)
      {
//      int g = ng-1;
         float_sw4* zp = a_Z[g].c_ptr();
         size_t npts = ((size_t)(mGlobalDims[g][1]-mGlobalDims[g][0])/st+1)*
            ((mGlobalDims[g][3]-mGlobalDims[g][2])/st+1)*
            ((mGlobalDims[g][5]-mGlobalDims[g][4])/st+1+m_extraz[g]);

         if( !mEW->usingParallelFS() || g == 0 )
            m_parallel_io[g]->writer_barrier();
      
         size_t nptsloc  = ((size_t)(mWindow[g][1] - mWindow[g][0])/mImageSamplingFactor + 1)*
            ( (mWindow[g][3] - mWindow[g][2])/mImageSamplingFactor + 1)*
            ( (mWindow[g][5] - mWindow[g][4])/mImageSamplingFactor + 1 + m_extraz[g]);

         int ni = (mWindow[g][1]-mWindow[g][0])/st+1;
         int nij=ni*((mWindow[g][3]-mWindow[g][2])/st+1);
         if( m_double )
         {
            double* zfp = new double[nptsloc];
#pragma omp parallel for
	    for( int k=mWindow[g][4] ; k <= mWindow[g][5] ; k+=st )
	       for( int j=mWindow[g][2] ; j <= mWindow[g][3] ; j+=st )
		  for( int i=mWindow[g][0] ; i <= mWindow[g][1] ; i+=st )
		  {
		     size_t ind = (i-mWindow[g][0])/st+ni*(j-mWindow[g][2])/st+nij*(k-mWindow[g][4])/st;
		     zfp[ind] = (double) a_Z[g](i,j,k);
		  }
	    //	 for( size_t i = 0; i < nptsloc ; i++ )
	    //	    zfp[i] = zp[i];
            char cprec[]="double";
            m_parallel_io[g]->write_array( &fid, 1, zfp, offset, cprec );
            offset += npts*sizeof(double);
            delete[] zfp;
         }
         else
         {
            float* zfp = new float[nptsloc];
#pragma omp parallel for
	    for( int k=mWindow[g][4] ; k <= mWindow[g][5] ; k+=st )
	       for( int j=mWindow[g][2] ; j <= mWindow[g][3] ; j+=st )
		  for( int i=mWindow[g][0] ; i <= mWindow[g][1] ; i+=st )
		  {
		     size_t ind = (i-mWindow[g][0])/st+ni*(j-mWindow[g][2])/st+nij*(k-mWindow[g][4])/st;
		     zfp[ind] = (float) a_Z[g](i,j,k);
		  }
	    //	 for( size_t i = 0; i < nptsloc ; i++ )
	    //	    zfp[i] = zp[i];
            char cprec[]="float";
            m_parallel_io[g]->write_array( &fid, 1, zfp, offset, cprec );
            offset += npts*sizeof(float);
            delete[] zfp;
         }
      } // end for g (curvilinear)
      
   } // end if grid info
   
   if( iwrite )
      close(fid);
}
