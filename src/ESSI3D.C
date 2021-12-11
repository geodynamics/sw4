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

#include "ESSI3D.h"
#include "ESSI3DHDF5.h"
#include "EW.h"
#include "Require.h"
#include "mpi.h"

int ESSI3D::mPreceedZeros = 0;
int ESSI3D::mNumberOfTimeSteps = -1;

ESSI3D* ESSI3D::nil = static_cast<ESSI3D*>(0);

//-----------------------------------------------------------------------
ESSI3D::ESSI3D(EW* a_ew, const std::string& filePrefix, int dumpInterval,
               int bufferInterval, float_sw4 coordBox[6], float_sw4 depth,
               int precision, int compressionMode, double compressionPar)
    : mEW(a_ew),
      mFilePrefix(filePrefix),
      mFileName(""),
      m_isDefinedMPIWriters(false),
      m_memallocated(false),
      m_fileOpen(false),
      m_ihavearray(false),
      m_hdf5_time(0),
      m_cycle(-1),
      m_dumpInterval(-1),
      m_bufferInterval(1),
      m_nbufstep(0),
      mDepth(depth),
      m_precision(precision),
      m_compressionMode(compressionMode),
      m_compressionPar(compressionPar),
      m_hdf5helper(NULL) {
  // volimage subdomain x,y corner coordinates
  for (int d = 0; d < 2 * 2; d++) mCoordBox[d] = coordBox[d];

  set_dump_interval(dumpInterval);
  set_buffer_interval(bufferInterval);

  MPI_Comm_rank(a_ew->m_cartesian_communicator, &m_rank);
}

//-----------------------------------------------------------------------
ESSI3D::~ESSI3D() {
  if (m_memallocated) {
    if (m_precision == 4) {
      for (int i = 0; i < 3; i++) delete[] m_floatField[i];
      delete[] m_floatField;
    } else if (m_precision == 8) {
      for (int i = 0; i < 3; i++) delete[] m_doubleField[i];
      delete[] m_doubleField;
    }
  }

  if (m_hdf5helper != NULL) delete m_hdf5helper;
}

//-----------------------------------------------------------------------
void ESSI3D::set_ntimestep(int ntimestep) {
  m_ntimestep = ntimestep;
  return;
}

//-----------------------------------------------------------------------
void ESSI3D::set_restart(bool is_restart) {
  m_isRestart = is_restart;
  return;
}

//-----------------------------------------------------------------------
void ESSI3D::set_dump_interval(int a_dumpInterval) {
  if (a_dumpInterval > 0) m_dumpInterval = a_dumpInterval;
  return;
}

//-----------------------------------------------------------------------
void ESSI3D::set_buffer_interval(int a_bufferInterval) {
  if (a_bufferInterval > 0) m_bufferInterval = a_bufferInterval;
  return;
}

//-----------------------------------------------------------------------
void ESSI3D::setup() {
  const bool debug = false;
  int g = mEW->mNumberOfGrids - 1;   // top curvilinear grid only
  m_ihavearray = true;               // gets negated if we don't, below
  for (int dim = 0; dim < 2; dim++)  // i,j dims only
  {
    int boxStart = (int)floor(mCoordBox[2 * dim] / mEW->mGridSize[g]) + 1;
    int boxEnd = (int)ceil(mCoordBox[2 * dim + 1] / mEW->mGridSize[g]) + 1;
    // Bound it into valid indices
    int indexMax = (dim == 0) ? mEW->m_global_nx[g] : mEW->m_global_ny[g];
    mGlobalDims[2 * dim] = max(1, min(boxStart, indexMax));
    mGlobalDims[2 * dim + 1] = max(1, min(boxEnd, indexMax));
    // cout << "ESSI3D index ranges: x[" << dim << "] = ("
    //   << boxStart << " , " << boxEnd << ")" << endl;
    int ixStart = (dim == 0) ? mEW->m_iStartInt[g] : mEW->m_jStartInt[g];
    int ixEnd = (dim == 0) ? mEW->m_iEndInt[g] : mEW->m_jEndInt[g];
    // Check if the output box indices intersect with this proc's
    if (m_ihavearray && (boxStart <= ixEnd) && (boxEnd >= ixStart)) {
      mWindow[2 * dim] = max(ixStart, min(boxStart, ixEnd));
      mWindow[2 * dim + 1] = max(ixStart, min(boxEnd, ixEnd));
    } else {
      m_ihavearray = false;
    }
  }

  // In k direction, guestimate the index for the requested depth
  int kmax = (int)ceil(mDepth / mEW->mGridSize[g]) + 1;
  if (mDepth == 1) kmax = 1;
  mGlobalDims[4] = mEW->m_kStartInt[g];
  mGlobalDims[5] = min(kmax, mEW->m_kEndInt[g]);
  if (m_ihavearray) {
    mWindow[4] = mGlobalDims[4];
    mWindow[5] = mGlobalDims[5];
  } else {
    for (int dim = 0; dim < 3; dim++) {
      mWindow[2 * dim] = 0;
      mWindow[2 * dim + 1] = -1;
    }
  }

  if (debug && (m_rank == 0)) {
    cout << "Global index range:"
         << " i=(" << mGlobalDims[0] << " , " << mGlobalDims[1] << "), j=("
         << mGlobalDims[2] << " , " << mGlobalDims[3] << "), k=("
         << mGlobalDims[4] << " , " << mGlobalDims[5] << ")" << endl;
  }

  if (debug)
    cout << "ESSI3D rank: " << m_rank << ", local index range: i=("
         << mWindow[0] << " , " << mWindow[1] << "), j=(" << mWindow[2] << " , "
         << mWindow[3] << "), k=(" << mWindow[4] << " , " << mWindow[5] << ")"
         << endl;

  if (m_ihavearray) {
    size_t npts = ((size_t)(mWindow[1] - mWindow[0]) + 1) *
                  ((mWindow[3] - mWindow[2]) + 1) *
                  ((mWindow[5] - mWindow[4]) + 1);
    npts *= m_bufferInterval;

    if (m_precision == 4) {
      m_floatField = new float*[3];
      for (int i = 0; i < 3; i++) m_floatField[i] = new float[npts];
    } else if (m_precision == 8) {
      m_doubleField = new double*[3];
      for (int i = 0; i < 3; i++) m_doubleField[i] = new double[npts];
    }
  } else {
    if (m_precision == 4) {
      m_floatField = new float*[3];
      for (int i = 0; i < 3; i++) m_floatField[i] = new float[1];
    } else if (m_precision == 8) {
      m_doubleField = new double*[3];
      for (int i = 0; i < 3; i++) m_doubleField[i] = new double[1];
    }
  }
  m_memallocated = true;
  return;
}

//-----------------------------------------------------------------------
void ESSI3D::define_pio() { return; }

//-----------------------------------------------------------------------
void ESSI3D::setSteps(int a_steps) {
  mNumberOfTimeSteps = a_steps;
  char buffer[50];
  mPreceedZeros = snprintf(buffer, 50, "%d", a_steps);
  return;
}

//-----------------------------------------------------------------------
double ESSI3D::getHDF5Timings() {
  return m_hdf5_time;  // Note just for this rank & file instance
}

//-----------------------------------------------------------------------
void ESSI3D::update_image(int a_cycle, float_sw4 a_time, float_sw4 a_dt,
                          vector<Sarray>& a_U, std::string& a_path,
                          Sarray& a_Z) {
#ifdef USE_HDF5
  int o_cycle = a_cycle;
  double hdf5_time = MPI_Wtime();
  if (!m_fileOpen)  // must be first call, open file and write
    open_vel_file(a_cycle, a_path, a_time, a_Z);

  if (m_dumpInterval != -1) {
    if (a_cycle != 0 && a_cycle % m_dumpInterval != 0 &&
        a_cycle != mNumberOfTimeSteps)
      return;
    a_cycle /= m_dumpInterval;
  }

  write_image_hdf5(a_cycle, a_path, a_time, a_U);

  if (o_cycle == mNumberOfTimeSteps)  // last time step
    close_vel_file();

  m_hdf5_time += (MPI_Wtime() - hdf5_time);
#else
  if (m_rank == 0)
    std::cout << "Using SSI output but SW4 is not compiled with HDF5"
              << std::endl;
#endif
  return;
}

//-----------------------------------------------------------------------
void ESSI3D::force_write_image(float_sw4 a_time, int a_cycle,
                               vector<Sarray>& a_U, std::string& a_path,
                               Sarray& a_Z) {
#ifdef USE_HDF5
  double hdf5_time = MPI_Wtime();
  open_vel_file(a_cycle, a_path, a_time, a_Z);
  write_image_hdf5(a_cycle, a_path, a_time, a_U);
  close_vel_file();
  m_hdf5_time += (MPI_Wtime() - hdf5_time);
#else
  if (m_rank == 0)
    std::cout << "Using SSI output but SW4 is not compiled with HDF5"
              << std::endl;
#endif
  return;
}

//-----------------------------------------------------------------------
void ESSI3D::compute_image(Sarray& a_A, int a_comp, int cycle) {
  int g = mEW->mNumberOfGrids - 1;

  int il = mEW->m_iStart[g];
  int iu = mEW->m_iEnd[g];
  int jl = mEW->m_jStart[g];
  int ju = mEW->m_jEnd[g];
  int kl = mEW->m_kStart[g];
  int ku = mEW->m_kEnd[g];

  // int niw = (mWindow[1]-mWindow[0])+1;
  // int nijw=niw*((mWindow[3]-mWindow[2])+1);
  size_t nkw = (mWindow[5] - mWindow[4]) + 1;
  size_t njkw = nkw * ((mWindow[3] - mWindow[2]) + 1);
  size_t nijkw = njkw * ((mWindow[1] - mWindow[0]) + 1);
  const int c = a_comp + 1;  // a_A array is 1-based
#pragma omp parallel for
  for (int k = mWindow[4]; k <= mWindow[5]; k++)
    for (int j = mWindow[2]; j <= mWindow[3]; j++)
      for (int i = mWindow[0]; i <= mWindow[1]; i++) {
        // HDF5 expects k major, then j, i
        // size_t ind = (i-mWindow[0])+niw*(j-mWindow[2])+nijw*(k-mWindow[4]);
        size_t ind =
            njkw * (i - mWindow[0]) + nkw * (j - mWindow[2]) + (k - mWindow[4]);
        ind += m_nbufstep * nijkw;
        if (m_precision == 4)
          m_floatField[a_comp][ind] = (float)a_A(c, i, j, k);
        else if (m_precision == 8)
          m_doubleField[a_comp][ind] = (double)a_A(c, i, j, k);
      }

  return;
}

//-----------------------------------------------------------------------
void ESSI3D::compute_file_suffix(int cycle, std::stringstream& fileSuffix) {
  fileSuffix << mFilePrefix << "."
             << "ssi";
  return;
}

//-----------------------------------------------------------------------
void ESSI3D::write_image(int cycle, std::string& path, float_sw4 t,
                         Sarray& a_Z) {
  return;
}

#ifdef USE_HDF5
void ESSI3D::open_vel_file(int a_cycle, std::string& a_path, float_sw4 a_time,
                           Sarray& a_Z) {
  bool debug = false;
  /* debug = true; */

  int g = mEW->mNumberOfGrids - 1;
  int window[6], global[3];
  for (int d = 0; d < 3; d++) {
    window[2 * d] = (m_ihavearray) ? (mWindow[2 * d] - mGlobalDims[2 * d]) : 0;
    window[2 * d + 1] =
        (m_ihavearray) ? (mWindow[2 * d + 1] - mGlobalDims[2 * d]) : -1;
    global[d] = mGlobalDims[2 * d + 1] - mGlobalDims[2 * d] + 1;
  }

  std::stringstream s, fileSuffix;
  compute_file_suffix(a_cycle, fileSuffix);
  if (a_path != ".") s << a_path;
  s << fileSuffix.str();  // string 's' is the file name including path
  m_hdf5helper =
      new ESSI3DHDF5(s.str(), global, window, m_ihavearray, m_precision);
  m_hdf5helper->set_ihavearray(m_ihavearray);

  if (debug && (m_rank == 0))
    cout << "Creating hdf5 file: " << m_hdf5helper->filename() << endl;

  double hdf5_time = MPI_Wtime();

  m_hdf5helper->create_file(m_isRestart);

  // Write header metadata
  double h = mEW->mGridSize[g];
  double lonlat_origin[2] = {mEW->getLonOrigin(), mEW->getLatOrigin()};
  double origin[3];
  for (int d = 0; d < 3; d++)
    origin[d] = (mGlobalDims[2 * d] - 1) * h;  // low end of each index range
  double az = mEW->getGridAzimuth();
  double dt = mEW->getTimeStep();

  if (!m_isRestart) {
    m_hdf5helper->write_header(h, lonlat_origin, az, origin, a_cycle, a_time,
                               dt);

    // Write z coodinates if necesito
    if (mEW->topographyExists()) {
      compute_image(a_Z, 0, 0);
      if (m_precision == 4)
        m_hdf5helper->write_topo(m_floatField[0]);
      else if (m_precision == 8)
        m_hdf5helper->write_topo(m_doubleField[0]);
    }
  }

  double hdf5_topo_time = MPI_Wtime();
  if (m_rank == 0)
    cout << "Create and write topo time " << MPI_Wtime() - hdf5_topo_time
         << endl;
  ;

  if (m_dumpInterval > 0) {
    int nstep = (int)ceil(m_ntimestep / m_dumpInterval);
    m_ntimestep = nstep;
    if (m_compressionMode > 0)
      m_hdf5helper->init_write_vel(m_isRestart, nstep, m_compressionMode,
                                   m_compressionPar, m_bufferInterval);
    else
      m_hdf5helper->init_write_vel(m_isRestart, nstep, 0, 0.0,
                                   m_bufferInterval);
  } else {
    if (m_compressionMode > 0)
      m_hdf5helper->init_write_vel(m_isRestart, m_ntimestep, m_compressionMode,
                                   m_compressionPar, m_bufferInterval);
    else
      m_hdf5helper->init_write_vel(m_isRestart, m_ntimestep, 0, 0.0,
                                   m_bufferInterval);
  }

  double hdf5_vel_time = MPI_Wtime();
  if (m_rank == 0)
    cout << "Create vel time " << hdf5_vel_time - hdf5_topo_time << endl;
  ;

  m_hdf5_time += (hdf5_vel_time - hdf5_time);

  m_nbufstep = 0;
  m_fileOpen = true;
  return;
}

void ESSI3D::close_vel_file() {
  if (m_fileOpen) {
    /*
    if (debug && (m_rank == 0))
      cout << "Closing hdf5 file: " << s.str() << endl;
    */

    m_hdf5helper->close_file();
    m_fileOpen = false;
  }
  return;
}

void ESSI3D::write_image_hdf5(int cycle, std::string& path, float_sw4 t,
                              vector<Sarray>& a_U) {
  // Top grid only
  int g = mEW->mNumberOfGrids - 1;
  int doWrite = 0;
  bool debug = false;
  /* debug = true; */

  for (int i = 0; i < 3; i++) {
    compute_image(a_U[g], i, cycle);
    if (cycle > 0 &&
        (m_nbufstep == m_bufferInterval - 1 || cycle == m_ntimestep)) {
      if (m_precision == 4)
        m_hdf5helper->write_vel((void*)m_floatField[i], i, cycle,
                                m_nbufstep + 1);
      else if (m_precision == 8)
        m_hdf5helper->write_vel((void*)m_doubleField[i], i, cycle,
                                m_nbufstep + 1);

      doWrite++;
    }
  }

  m_nbufstep++;
  m_nbufstep %= m_bufferInterval;

  if (!(doWrite == 3 || doWrite == 0))
    fprintf(stderr,
            "Error with essioutput write_image_hdf5, not all variables are "
            "written correctly!\n");

  if (doWrite == 3) {
    if (debug)
      fprintf(stderr, "Rank %d: write_image_hdf5 cycle=%d/%d, m_nbufstep=%d\n",
              m_rank, cycle, m_ntimestep, m_nbufstep);
    m_nbufstep = 0;
  }
  return;
}
#endif  // ifdef USE_HDF5
