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
               float_sw4 coordBox[6], float_sw4 depth)
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
      mDepth(depth),
      m_hdf5helper(NULL) {
  // volimage subdomain x,y corner coordinates
  for (int d = 0; d < 2 * 2; d++) mCoordBox[d] = coordBox[d];

  set_dump_interval(dumpInterval);
}

//-----------------------------------------------------------------------
ESSI3D::~ESSI3D() {
  if (m_memallocated) delete[] m_doubleField;

  if (m_hdf5helper != NULL) delete m_hdf5helper;
}

//-----------------------------------------------------------------------
void ESSI3D::set_ntimestep(int ntimestep) { m_ntimestep = ntimestep; }

//-----------------------------------------------------------------------
void ESSI3D::set_dump_interval(int a_dumpInterval) {
  if (a_dumpInterval > 0) m_dumpInterval = a_dumpInterval;
}

//-----------------------------------------------------------------------
void ESSI3D::setup() {
  const bool debug = false;
  //  MPI_Comm comm = MPI_COMM_WORLD;
  // MPI_Info info = MPI_INFO_NULL;
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

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

  if (debug && (myRank == 0)) {
    cout << "Global index range:"
         << " i=(" << mGlobalDims[0] << " , " << mGlobalDims[1] << "), j=("
         << mGlobalDims[2] << " , " << mGlobalDims[3] << "), k=("
         << mGlobalDims[4] << " , " << mGlobalDims[5] << ")" << endl;
  }

  // TODO - pseudocode:
  // - figure out depth vs. zmin/zmax indices
  // Figure out x,y index ranges that contain requested coord box
  // cout << "ESSI3D ptr: " << this << endl;
  if (debug)
    cout << "ESSI3D rank: " << myRank << ", local index range: i=("
         << mWindow[0] << " , " << mWindow[1] << "), j=(" << mWindow[2] << " , "
         << mWindow[3] << "), k=(" << mWindow[4] << " , " << mWindow[5] << ")"
         << endl;

  if (m_ihavearray) {
    size_t npts = ((size_t)(mWindow[1] - mWindow[0]) + 1) *
                  ((mWindow[3] - mWindow[2]) + 1) *
                  ((mWindow[5] - mWindow[4]) + 1);
    m_doubleField = new double[npts];
  } else {
    m_doubleField = new double[1];
  }
  m_memallocated = true;
#ifdef USE_HDF5
  define_pio_hdf5();
#else
  define_pio();
#endif
}

//-----------------------------------------------------------------------
void ESSI3D::define_pio() {
  // TODO - error out
}

//-----------------------------------------------------------------------
void ESSI3D::setSteps(int a_steps) {
  mNumberOfTimeSteps = a_steps;
  char buffer[50];
  mPreceedZeros = snprintf(buffer, 50, "%d", a_steps);
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
  double hdf5_time = MPI_Wtime();
  if (!m_fileOpen)  // must be first call, open file and write
    open_vel_file(a_cycle, a_path, a_time, a_Z);

  if ((a_cycle > 0) &&                  // don't close and open on first cycle
      (m_dumpInterval != -1) &&         // dump interval was set
      (a_cycle % m_dumpInterval == 0))  // close, open new
  {
    close_vel_file();
    open_vel_file(a_cycle, a_path, a_time, a_Z);
  }
  write_image_hdf5(a_cycle, a_path, a_time, a_U);
  if (a_cycle == mNumberOfTimeSteps)  // last time step
    close_vel_file();

  m_hdf5_time += (MPI_Wtime() - hdf5_time);
#endif
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
#endif
}

//-----------------------------------------------------------------------
void ESSI3D::compute_image(Sarray& a_A, int a_comp) {
  int g = mEW->mNumberOfGrids - 1;

  int il = mEW->m_iStart[g];
  int iu = mEW->m_iEnd[g];
  int jl = mEW->m_jStart[g];
  int ju = mEW->m_jEnd[g];
  int kl = mEW->m_kStart[g];
  int ku = mEW->m_kEnd[g];
  //  int ni = (iu - il + 1);
  // int nj = (ju - jl + 1);

  // int niw = (mWindow[1]-mWindow[0])+1;
  // int nijw=niw*((mWindow[3]-mWindow[2])+1);
  int nkw = (mWindow[5] - mWindow[4]) + 1;
  int njkw = nkw * ((mWindow[3] - mWindow[2]) + 1);
  const int c = a_comp + 1;  // why comp+1?
#pragma omp parallel for
  for (int k = mWindow[4]; k <= mWindow[5]; k++)
    for (int j = mWindow[2]; j <= mWindow[3]; j++)
      for (int i = mWindow[0]; i <= mWindow[1]; i++) {
        // HDF5 expects k major, then j, i
        // size_t ind = (i-mWindow[0])+niw*(j-mWindow[2])+nijw*(k-mWindow[4]);
        size_t ind =
            njkw * (i - mWindow[0]) + nkw * (j - mWindow[2]) + (k - mWindow[4]);
        m_doubleField[ind] = (double)a_A(c, i, j, k);
      }
}

//-----------------------------------------------------------------------
void ESSI3D::compute_file_suffix(int cycle, std::stringstream& fileSuffix) {
  fileSuffix << mFilePrefix << ".cycle=";
  int temp = static_cast<int>(pow(10.0, mPreceedZeros - 1));
  int testcycle = cycle;
  if (cycle == 0) testcycle = 1;
  while (testcycle < temp) {
    fileSuffix << "0";
    temp /= 10;
  }
  fileSuffix << cycle;
  fileSuffix << "."
             << "essi";
}

//-----------------------------------------------------------------------
void ESSI3D::write_image(int cycle, std::string& path, float_sw4 t,
                         Sarray& a_Z) {
  // TODO - error out
}

#ifdef USE_HDF5
//-----------------------------------------------------------------------
void ESSI3D::define_pio_hdf5() {
  // TODO - error out
}

void ESSI3D::open_vel_file(int a_cycle, std::string& a_path, float_sw4 a_time,
                           Sarray& a_Z) {
  //  herr_t ierr;
  // hid_t dataspace_id;
  // hid_t dataset_id;
  // hid_t prop_id;

  // Parameters for extendible dataset
  const int num_dims = 3 + 1;  // 3 space + 1 time that will be extendible
  hsize_t dims[num_dims] = {0, 0, 0, H5S_UNLIMITED};
  // hsize_t slice_dims[num_dims];
  //  hsize_t block_dims[num_dims];
  // hsize_t global_dims[num_dims];

  int g = mEW->mNumberOfGrids - 1;
  int window[6], global[3];
  for (int d = 0; d < 3; d++) {
    dims[d] = mWindow[2 * d + 1] - mWindow[2 * d] + 1;
    // slice_dims[d] = dims[d];
    // block_dims[d] = dims[d];
    // global_dims[d] = mGlobalDims[2 * d + 1] - mGlobalDims[2 * d] + 1;

    window[2 * d] = (m_ihavearray) ? (mWindow[2 * d] - mGlobalDims[2 * d]) : 0;
    window[2 * d + 1] =
        (m_ihavearray) ? (mWindow[2 * d + 1] - mGlobalDims[2 * d]) : -1;
    global[d] = mGlobalDims[2 * d + 1] - mGlobalDims[2 * d] + 1;
  }
  // slice_dims[0] = 1;              // Make sure we're chunking on slices
  // slice_dims[num_dims - 1] = 1;   // 1 time step per write
  // block_dims[num_dims - 1] = 1;   // 1 time step per write
  // global_dims[num_dims - 1] = 1;  // 1 time step per write

  // Just to see what's being copied correctly
  for (int i = 0; i < dims[0]; i++)
    for (int j = 0; j < dims[1]; j++)
      for (int k = 0; k < dims[2]; k++) {
        size_t ind = i + dims[0] * (j + dims[1] * k);
        m_doubleField[ind] = (double)ind;
      }

  /*
  if (debug && (myRank == 0))
  {
    cout << "Setting dims to: " << endl;
    for (int d=0; d < num_dims; d++)
      cout << "dims[" << d << "] = " << dims[d]
        << ", slice_dims[" << d << "] = " << slice_dims[d] << endl;
  }
  */

  std::stringstream s, fileSuffix;
  compute_file_suffix(a_cycle, fileSuffix);
  if (a_path != ".") s << a_path;
  s << fileSuffix.str();  // string 's' is the file name including path
  m_hdf5helper = new ESSI3DHDF5(s.str(), global, window, m_ihavearray);
  m_hdf5helper->set_ihavearray(m_ihavearray);
  /*
  if (debug && (myRank == 0))
     cout << "Creating hdf5 file: " << m_hdf5helper->filename() << endl;
  */
  double hdf5_time = MPI_Wtime();
  m_hdf5helper->create_file();

  // Write header metadata
  double h = mEW->mGridSize[g];
  double lonlat_origin[2] = {mEW->getLonOrigin(), mEW->getLatOrigin()};
  double origin[3];
  for (int d = 0; d < 3; d++)
    origin[d] = (mGlobalDims[2 * d] - 1) * h;  // low end of each index range
  double az = mEW->getGridAzimuth();
  double dt = mEW->getTimeStep();
  m_hdf5helper->write_header(h, lonlat_origin, az, origin, a_cycle, a_time, dt);

  // Write z coodinates if necesito
  if (mEW->topographyExists()) {
    compute_image(a_Z, 0);
    m_hdf5helper->write_topo(m_doubleField);
  }

  /*
  if (debug && (myRank == 0))
     cout << "Creating extendible hdf5 velocity fields..." << endl;
  */

  if (m_dumpInterval > 0)
    m_hdf5helper->init_write_vel(m_dumpInterval);
  else
    m_hdf5helper->init_write_vel(m_ntimestep);
  m_hdf5_time += (MPI_Wtime() - hdf5_time);

  m_fileOpen = true;
}

void ESSI3D::close_vel_file() {
  if (m_fileOpen) {
    /*
    if (debug && (myRank == 0))
      cout << "Closing hdf5 file: " << s.str() << endl;
    */

    m_hdf5helper->close_file();
    m_fileOpen = false;
  }
}

void ESSI3D::write_image_hdf5(int cycle, std::string& path, float_sw4 t,
                              vector<Sarray>& a_U) {
  // Top grid only
  int g = mEW->mNumberOfGrids - 1;
  compute_image(a_U[g], 0);
  m_hdf5helper->write_vel(m_doubleField, 0, cycle);
  compute_image(a_U[g], 1);
  m_hdf5helper->write_vel(m_doubleField, 1, cycle);
  compute_image(a_U[g], 2);
  m_hdf5helper->write_vel(m_doubleField, 2, cycle);

  /*
  MPI_Barrier(MPI_COMM_WORLD);
  if (m_fileOpen)
  {
    if (debug && (myRank == 0))
      cout << "Closing hdf5 file: " << s.str() << endl;

    m_hdf5helper->close_file();
    m_fileOpen = false;
  }
  */
}
#endif  // ifdef USE_HDF5
