#include <fcntl.h>
#include <unistd.h>

#include <cstdio>
#include <cstring>
#include <ctime>

#include "CheckPoint.h"
#include "EW.h"
#include "Require.h"
#include "mpi.h"

#ifdef USE_ZFP
#include "H5Zzfp_lib.h"
#include "H5Zzfp_props.h"
#endif

CheckPoint* CheckPoint::nil = static_cast<CheckPoint*>(0);

//-----------------------------------------------------------------------
// Null constructor, makes it possible to always query 'DoCheckPointing' and
// 'DoRestart'
CheckPoint::CheckPoint(EW* a_ew)
    : mEW(a_ew),
      mWritingCycle(-1),
      mCycleInterval(0),
      mStartTime(0.0),
      mCheckPointFile(" "),
      mPreceedZeros(0),
      m_double(true),
      mRestartFile(" "),
      m_winallocated(false),
      m_bufsize(0),
      m_fileno(0),
      mDoCheckPointing(false),
      mRestartPathSet(false),
      mDoRestart(false),
      m_es_id(0),
      m_kji_order(true) {}

//-----------------------------------------------------------------------
// Save check point files, but no restart
CheckPoint::CheckPoint(EW* a_ew, int cycle, int cycleInterval, string fname,
                       size_t bufsize)
    : mEW(a_ew),
      mWritingCycle(cycle),
      mCycleInterval(cycleInterval),
      mStartTime(0.0),
      mCheckPointFile(fname),
      mPreceedZeros(0),
      m_double(true),
      mRestartFile("restart"),
      m_winallocated(false),
      m_bufsize(bufsize),
      m_fileno(0),
      mDoCheckPointing(true),
      mRestartPathSet(false),
      mDoRestart(false),
      m_es_id(0),
      m_kji_order(true) {
  m_double = sizeof(float_sw4) == sizeof(double);
}

//-----------------------------------------------------------------------
// Restart, but not initialized for saving check points
CheckPoint::CheckPoint(EW* a_ew, string fname, size_t bufsize)
    : mEW(a_ew),
      mWritingCycle(-1),
      mCycleInterval(0),
      mStartTime(0.0),
      mCheckPointFile("chkpt"),
      mPreceedZeros(0),
      m_double(true),
      mRestartFile(fname),
      m_winallocated(false),
      m_bufsize(bufsize),
      m_fileno(0),
      mDoCheckPointing(false),
      mRestartPathSet(false),
      m_es_id(0),
      mDoRestart(true) {
  m_double = sizeof(float_sw4) == sizeof(double);
}

//-----------------------------------------------------------------------
CheckPoint::~CheckPoint() {
  if (m_winallocated)
    for (int g = 0; g < mEW->mNumberOfGrids; g++) {
      delete[] mWindow[g];
      delete[] mGlobalDims[g];
    }
}

//-----------------------------------------------------------------------
bool CheckPoint::do_checkpointing() { return mDoCheckPointing; }

//-----------------------------------------------------------------------
int CheckPoint::get_checkpoint_cycle_interval() { return mCycleInterval; }

//-----------------------------------------------------------------------
bool CheckPoint::do_restart() { return mDoRestart; }

//-----------------------------------------------------------------------
void CheckPoint::setup_sizes() {
  if (mDoRestart || mDoCheckPointing) {
    if (!m_winallocated) {
      mWindow.resize(mEW->mNumberOfGrids);
      mGlobalDims.resize(mEW->mNumberOfGrids);
      for (int g = 0; g < mEW->mNumberOfGrids; g++) {
        mWindow[g] = new int[6];
        mGlobalDims[g] = new int[6];
      }
      m_winallocated = true;
    }

    m_ihavearray.resize(mEW->mNumberOfGrids);

    int ghost_points = mEW->getNumberOfGhostPoints();
    // tmp
    //      printf("CheckPoint: Number of ghost points = %d\n", ghost_points);

    for (int g = 0; g < mEW->mNumberOfGrids; g++) {
      // With attenuation, the memory variables must be saved at all ghost
      // points. This is because the memory variables satisfy ODEs at all
      // points, and do not satify any boundary conditions

      if (mEW->getLocalBcType(g, 0) == bProcessor)
        mWindow[g][0] = mEW->m_iStartInt[g];
      else
        mWindow[g][0] = mEW->m_iStartInt[g] - ghost_points;
      ;

      if (mEW->getLocalBcType(g, 1) == bProcessor)
        mWindow[g][1] = mEW->m_iEndInt[g];
      else
        mWindow[g][1] = mEW->m_iEndInt[g] + ghost_points;

      if (mEW->getLocalBcType(g, 2) == bProcessor)
        mWindow[g][2] = mEW->m_jStartInt[g];
      else
        mWindow[g][2] = mEW->m_jStartInt[g] - ghost_points;
      ;

      if (mEW->getLocalBcType(g, 3) == bProcessor)
        mWindow[g][3] = mEW->m_jEndInt[g];
      else
        mWindow[g][3] = mEW->m_jEndInt[g] + ghost_points;

      // all points in k-dir are local to each proc
      mWindow[g][4] = mEW->m_kStartInt[g] -
                      ghost_points;  // need 1 ghost point at MR interface
      mWindow[g][5] =
          mEW->m_kEndInt[g] + ghost_points;  // and all ghost points at bottom

      // Need to store ghost point values outside the physical boundaries for
      // the attenuation memory variables, because they satisfy ODEs and not BC

      mGlobalDims[g][0] = 1 - ghost_points;
      mGlobalDims[g][1] = mEW->m_global_nx[g] + ghost_points;
      mGlobalDims[g][2] = 1 - ghost_points;
      mGlobalDims[g][3] = mEW->m_global_ny[g] + ghost_points;
      // k-dir is local
      mGlobalDims[g][4] = mWindow[g][4];
      mGlobalDims[g][5] = mWindow[g][5];

      // The 3D array is assumed to span the entire computational domain
      m_ihavearray[g] = true;
    }
    setSteps(mEW->getNumberOfTimeSteps());
    // tmp
    if (mEW->proc_zero())
      cout << "Checkpoint::setup_sizes: Calling define_pio()..." << endl;

    define_pio();
  }  // end if doRestart || doCheckpointing

  //   cout << "mwind = " << mWindow[0][4] << " " << mWindow[0][5] << endl;
  //   cout << "globaldims = " << mGlobalDims[0][4] << " " << mGlobalDims[0][5]
  //   << endl;
}

//-----------------------------------------------------------------------
void CheckPoint::define_pio() {
  int glow = 0, ghigh = mEW->mNumberOfGrids;

  double time_start = MPI_Wtime();
  double time_measure[12];
  time_measure[0] = time_start;

  // Create the restart directory if it doesn't exist
  //
  // AP: On the burst buffer at Cori it takes *forever* to build a directory
  // For now, assume the directory is already there
  //
  if (mRestartPathSet) mEW->create_directory(mRestartPath);

  m_parallel_io = new Parallel_IO*[ghigh - glow + 1];
  for (int g = glow; g < ghigh; g++) {
    // tmp
    if (mEW->proc_zero()) cout << "setup_sizes for grid g = " << g << endl;

    int global[3], local[3], start[3];
    for (int dim = 0; dim < 3; dim++) {
      global[dim] = mGlobalDims[g][2 * dim + 1] - mGlobalDims[g][2 * dim] + 1;
      local[dim] = mWindow[g][2 * dim + 1] - mWindow[g][2 * dim] + 1;
      start[dim] = mWindow[g][2 * dim] - mGlobalDims[g][2 * dim];
    }

    int iwrite = 0;
    int nrwriters = mEW->getNumberOfWritersPFS();
    int nproc = 0, myid = 0;
    MPI_Comm_size(mEW->m_cartesian_communicator, &nproc);
    MPI_Comm_rank(mEW->m_cartesian_communicator, &myid);

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
    if (m_kji_order) {
      // Swap i and k on file
      int tmp = global[0];
      global[0] = global[2];
      global[2] = tmp;
      tmp = local[0];
      local[0] = local[2];
      local[2] = tmp;
      tmp = start[0];
      start[0] = start[2];
      start[2] = tmp;
    }
    if (mEW->proc_zero())
      cout << "Creating a Parallel_IO object for grid g = " << g << endl;
    m_parallel_io[g - glow] =
        new Parallel_IO(iwrite, mEW->usingParallelFS(), global, local, start,
                        m_bufsize);
    // tmp
    if (mEW->proc_zero())
      cout << "Done creating the Parallel_IO object" << endl;
    delete[] owners;
  }
}

//-----------------------------------------------------------------------
void CheckPoint::setSteps(int a_steps) {
  char buffer[50];
  mPreceedZeros = snprintf(buffer, 50, "%d", a_steps);
}

//-----------------------------------------------------------------------
bool CheckPoint::timeToWrite(float_sw4 time, int cycle, float_sw4 dt) {
  if (!mDoCheckPointing) return false;

  // Will we write at this time step (cycle) ?
  bool do_it = false;
  if (cycle == mWritingCycle) do_it = true;
  if (mCycleInterval != 0 && cycle % mCycleInterval == 0 && time >= mStartTime)
    do_it = true;
  return do_it;
}

//-----------------------------------------------------------------------
void CheckPoint::compute_file_suffix(int cycle, std::stringstream& fileSuffix) {
  fileSuffix << mCheckPointFile << ".cycle=";
  int temp = static_cast<int>(pow(10.0, mPreceedZeros - 1));
  int testcycle = cycle;
  if (cycle == 0) testcycle = 1;
  while (testcycle < temp) {
    fileSuffix << "0";
    temp /= 10;
  }
  fileSuffix << cycle;
  fileSuffix << ".sw4checkpoint";
}

//-----------------------------------------------------------------------
void CheckPoint::write_checkpoint(float_sw4 a_time, int a_cycle,
                                  vector<Sarray>& a_Um, vector<Sarray>& a_U,
                                  vector<Sarray*>& a_AlphaVEm,
                                  vector<Sarray*>& a_AlphaVE) {
  //
  // File format:
  //
  //    header (see routine write_header)
  //    for g=1,number of grids
  //       Um(g) (3 component float/double array)
  //       U(g)  (3 component float/double array)
  //       for m=1,number of mechanisms
  //           AlphaVEm(g,m) (3 component float/double array)
  //           AlphaVE(g,m)  (3 component float/double array)
  //       endfor
  //    endfor
  //

  //
  // Would it be possible to save the entire input file to the restart file ?
  //
  int ng = mEW->mNumberOfGrids;
  //   off_t offset = (4+6*ng)*sizeof(int) + 2*sizeof(float_sw4);

  bool iwrite = false;
  for (int g = 0; g < ng; g++) iwrite = iwrite || m_parallel_io[g]->i_write();

  std::stringstream s;
  if (iwrite) {
    std::stringstream fileSuffix;
    compute_file_suffix(a_cycle, fileSuffix);
    if (mRestartPathSet)
      s << mRestartPath << "/";
    else if (mEW->getPath() != "./")
      s << mEW->getPath() << "/";
    s << fileSuffix.str();
  }

  // Keep track of the number of files, save previous file name, and delete the
  // second last.
  cycle_checkpoints(s.str());

  // Open file from processor zero and write header.
  int hsize;
  int fid = -1;
  if (m_parallel_io[0]->proc_zero()) {
    fid = open(const_cast<char*>(s.str().c_str()), O_CREAT | O_TRUNC | O_WRONLY,
               0660);
    CHECK_INPUT(fid != -1,
                "CheckPoint::write_file: Error opening: " << s.str());
    int myid;

    MPI_Comm_rank(mEW->m_cartesian_communicator, &myid);
    std::cout << "writing check point on file " << s.str() << " using "
              << m_parallel_io[0]->n_writers() << " writers" << std::endl;
    write_header(fid, a_time, a_cycle, hsize);
    fsync(fid);
  }
  //   m_parallel_io[0]->writer_barrier();
  int bcast_root = m_parallel_io[0]->proc_zero_rank_in_comm_world();
  MPI_Bcast(&hsize, 1, MPI_INT, bcast_root, mEW->m_cartesian_communicator);
  off_t offset = hsize;

  // Open file from all writers
  if (iwrite && !m_parallel_io[0]->proc_zero()) {
    fid = open(const_cast<char*>(s.str().c_str()), O_WRONLY);
    CHECK_INPUT(fid != -1,
                "CheckPoint::write_checkpoint:: Error opening: " << s.str());
  }

  // Write data blocks
  char cprec[] = "double";
  if (!m_double) strcpy(cprec, "float");
  for (int g = 0; g < ng; g++) {
    size_t npts = ((size_t)(mGlobalDims[g][1] - mGlobalDims[g][0] + 1)) *
                  ((size_t)(mGlobalDims[g][3] - mGlobalDims[g][2] + 1)) *
                  ((size_t)(mGlobalDims[g][5] - mGlobalDims[g][4] + 1));

    if (!mEW->usingParallelFS() || g == 0) m_parallel_io[g]->writer_barrier();

    size_t nptsloc = (size_t)(mWindow[g][1] - mWindow[g][0] + 1) *
                     (mWindow[g][3] - mWindow[g][2] + 1) *
                     (mWindow[g][5] - mWindow[g][4] + 1);

    // allocate local buffer array
    float_sw4* doubleField = new float_sw4[3 * nptsloc];
    if (m_kji_order) {
      a_Um[g].extract_subarrayIK(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                                 mWindow[g][3], mWindow[g][4], mWindow[g][5],
                                 doubleField);
      m_parallel_io[g]->write_array(&fid, 3, doubleField, offset, cprec);
      offset += 3 * npts * sizeof(float_sw4);

      a_U[g].extract_subarrayIK(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                                mWindow[g][3], mWindow[g][4], mWindow[g][5],
                                doubleField);
      m_parallel_io[g]->write_array(&fid, 3, doubleField, offset, cprec);
      offset += 3 * npts * sizeof(float_sw4);
      for (int m = 0; m < mEW->getNumberOfMechanisms(); m++) {
        a_AlphaVEm[g][m].extract_subarrayIK(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
        m_parallel_io[g]->write_array(&fid, 3, doubleField, offset, cprec);
        offset += 3 * npts * sizeof(float_sw4);

        a_AlphaVE[g][m].extract_subarrayIK(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
        m_parallel_io[g]->write_array(&fid, 3, doubleField, offset, cprec);
        offset += 3 * npts * sizeof(float_sw4);
      }
    } else {
      a_Um[g].extract_subarray(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                               mWindow[g][3], mWindow[g][4], mWindow[g][5],
                               doubleField);
      m_parallel_io[g]->write_array(&fid, 3, doubleField, offset, cprec);
      offset += 3 * npts * sizeof(float_sw4);

      a_U[g].extract_subarray(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                              mWindow[g][3], mWindow[g][4], mWindow[g][5],
                              doubleField);
      m_parallel_io[g]->write_array(&fid, 3, doubleField, offset, cprec);
      offset += 3 * npts * sizeof(float_sw4);

      for (int m = 0; m < mEW->getNumberOfMechanisms(); m++) {
        a_AlphaVEm[g][m].extract_subarray(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
        m_parallel_io[g]->write_array(&fid, 3, doubleField, offset, cprec);
        offset += 3 * npts * sizeof(float_sw4);

        a_AlphaVE[g][m].extract_subarray(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
        m_parallel_io[g]->write_array(&fid, 3, doubleField, offset, cprec);
        offset += 3 * npts * sizeof(float_sw4);
      }
    }
    delete[] doubleField;
  }
  if (iwrite) close(fid);
}  // end write_checkpoint()

//-----------------------------------------------------------------------
void CheckPoint::read_checkpoint(float_sw4& a_time, int& a_cycle,
                                 vector<Sarray>& a_Um, vector<Sarray>& a_U,
                                 vector<Sarray*>& a_AlphaVEm,
                                 vector<Sarray*>& a_AlphaVE) {
  //
  // It is assumed that the arrays are already declared with the right
  // dimensions. This routine will check that sizes match, but will not
  // allocate or resize the arrays Um and U, AlphaVEm, AlphaVE
  //
  int ng = mEW->mNumberOfGrids;
  //   off_t offset = (4+6*ng)*sizeof(int) + 2*sizeof(float_sw4);
  bool iread = false;
  for (int g = 0; g < ng; g++) iread = iread || m_parallel_io[g]->i_write();

  std::stringstream s;
  if (iread) {
    if (mRestartPathSet)
      s << mRestartPath << "/";
    else if (mEW->getPath() != "./")
      s << mEW->getPath();
    s << mRestartFile;
  }

  // Open file from processor zero and read header.
  int fid = -1;
  int hsize;
  if (m_parallel_io[0]->proc_zero()) {
    fid = open(const_cast<char*>(s.str().c_str()), O_RDONLY);
    CHECK_INPUT(fid != -1,
                "CheckPoint::read_checkpoint: Error opening: " << s.str());
    int myid;

    MPI_Comm_rank(mEW->m_cartesian_communicator, &myid);
    std::cout << "reading check point on file " << s.str() << endl;
    read_header(fid, a_time, a_cycle, hsize);
  }
  //   m_parallel_io[0]->writer_barrier();

  // Broadcast read information to all other processors.
  int bcast_root = m_parallel_io[0]->proc_zero_rank_in_comm_world();
  MPI_Bcast(&a_cycle, 1, MPI_INT, bcast_root, mEW->m_cartesian_communicator);
  MPI_Bcast(&a_time, 1, mEW->m_mpifloat, bcast_root,
            mEW->m_cartesian_communicator);
  MPI_Bcast(&hsize, 1, MPI_INT, bcast_root, mEW->m_cartesian_communicator);
  off_t offset = hsize;

  // Open file from all readers
  if (iread && !m_parallel_io[0]->proc_zero()) {
    fid = open(const_cast<char*>(s.str().c_str()), O_RDONLY);
    CHECK_INPUT(fid != -1, "CheckPoint::read_checkpoint:: Error opening file: "
                               << s.str());
  }

  // Read data blocks.
  char cprec[] = "double";
  if (!m_double) strcpy(cprec, "float");

  for (int g = 0; g < ng; g++) {
    size_t npts = ((size_t)(mGlobalDims[g][1] - mGlobalDims[g][0] + 1)) *
                  ((size_t)(mGlobalDims[g][3] - mGlobalDims[g][2] + 1)) *
                  ((size_t)(mGlobalDims[g][5] - mGlobalDims[g][4] + 1));

    // size_t nptsloc = ((size_t)(mEW->m_iEnd[g]-mEW->m_iStart[g]+1))*
    //                  ((size_t)(mEW->m_jEnd[g]-mEW->m_jStart[g]+1))*
    //                  ((size_t)(mEW->m_kEnd[g]-mEW->m_kStart[g]+1));

    size_t nptsloc = (size_t)(mWindow[g][1] - mWindow[g][0] + 1) *
                     (mWindow[g][3] - mWindow[g][2] + 1) *
                     (mWindow[g][5] - mWindow[g][4] + 1);

    if (!mEW->usingParallelFS() || g == 0) m_parallel_io[g]->writer_barrier();

    // array with ghost points in k-ONLY read into doubleField,
    float_sw4* doubleField = new float_sw4[3 * nptsloc];
    if (m_kji_order) {
      m_parallel_io[g]->read_array(&fid, 3, doubleField, offset, cprec);
      offset += 3 * npts * sizeof(float_sw4);
      a_Um[g].insert_subarrayIK(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                                mWindow[g][3], mWindow[g][4], mWindow[g][5],
                                doubleField);

      m_parallel_io[g]->read_array(&fid, 3, doubleField, offset, cprec);
      offset += 3 * npts * sizeof(float_sw4);
      a_U[g].insert_subarrayIK(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                               mWindow[g][3], mWindow[g][4], mWindow[g][5],
                               doubleField);

      for (int m = 0; m < mEW->getNumberOfMechanisms(); m++) {
        m_parallel_io[g]->read_array(&fid, 3, doubleField, offset, cprec);
        offset += 3 * npts * sizeof(float_sw4);
        a_AlphaVEm[g][m].insert_subarrayIK(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);

        m_parallel_io[g]->read_array(&fid, 3, doubleField, offset, cprec);
        offset += 3 * npts * sizeof(float_sw4);
        a_AlphaVE[g][m].insert_subarrayIK(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
      }
    } else {
      m_parallel_io[g]->read_array(&fid, 3, doubleField, offset, cprec);
      offset += 3 * npts * sizeof(float_sw4);
      a_Um[g].insert_subarray(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                              mWindow[g][3], mWindow[g][4], mWindow[g][5],
                              doubleField);

      m_parallel_io[g]->read_array(&fid, 3, doubleField, offset, cprec);
      offset += 3 * npts * sizeof(float_sw4);
      a_U[g].insert_subarray(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                             mWindow[g][3], mWindow[g][4], mWindow[g][5],
                             doubleField);

      for (int m = 0; m < mEW->getNumberOfMechanisms(); m++) {
        m_parallel_io[g]->read_array(&fid, 3, doubleField, offset, cprec);
        offset += 3 * npts * sizeof(float_sw4);
        a_AlphaVEm[g][m].insert_subarray(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);

        m_parallel_io[g]->read_array(&fid, 3, doubleField, offset, cprec);
        offset += 3 * npts * sizeof(float_sw4);
        a_AlphaVE[g][m].insert_subarray(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
      }
    }
    delete[] doubleField;
  }
  if (iread) close(fid);
}

//-----------------------------------------------------------------------
float_sw4 CheckPoint::getDt() {
  float_sw4 dt;
  if (mEW->getRank() == 0) {
    std::stringstream s;
    if (mRestartPathSet)
      s << mRestartPath << "/";
    else if (mEW->getPath() != "./")
      s << mEW->getPath();
    s << mRestartFile;  // string 's' is the file name including path

    if (mUseHDF5) {
#ifdef USE_HDF5
      hid_t fid = H5Fopen(const_cast<char*>(s.str().c_str()), H5F_ACC_RDONLY,
                          H5P_DEFAULT);
      CHECK_INPUT(fid > 0, "CheckPoint::read_checkpoint_hdf5: Error opening: "
                               << s.str());

      hid_t attr = H5Aopen(fid, "dt", H5P_DEFAULT);
      hid_t dtype = (m_double ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT);
      int ret = H5Aread(attr, dtype, &dt);
      H5Aclose(attr);
      CHECK_INPUT(ret >= 0, "CheckPoint::read_header_hdf5: Error reading dt");
      H5Fclose(fid);
#else
      if (m_parallel_io[0]->proc_zero())
        cout << "Configured to restart with HDF5 but SW4 is not compiled with "
                "HDF5!"
             << endl;
#endif
    } else {
      int fid = open(const_cast<char*>(s.str().c_str()), O_RDONLY);
      CHECK_INPUT(fid != -1, "CheckPoint::getDt: Error opening: " << s.str());
      lseek(fid, 3 * sizeof(int) + sizeof(float_sw4), SEEK_SET);
      size_t nr = read(fid, &dt, sizeof(float_sw4));
      CHECK_INPUT(
          nr == sizeof(float_sw4),
          "CheckPoint::getDt, error reading time step from restart file\n");
      close(fid);
    }
  }
  MPI_Bcast(&dt, 1, mEW->m_mpifloat, 0, mEW->m_cartesian_communicator);
  return dt;
}

//-----------------------------------------------------------------------
void CheckPoint::write_header(int& fid, float_sw4 a_time, int a_cycle,
                              int& hsize) {
  //
  // Header format: prec - precision 8--> double, 4--> single (int)
  //                ng   - Number of grids (int)
  //                time - Time (float)
  //                cycle- Time step number corresponding to Time (int)
  //                dt   - Time step size (float)
  //                nmech- Number of mechanisms in attenuation model
  //                dims(g,1:6) - Size of array on grid g,
  //                       dim(1) <= i <= dim(2), dim(3) <= j <= dim(4)
  //                       dim(5) <= k <= dim(6)
  //
  int prec = m_double ? 8 : 4;
  size_t ret = write(fid, &prec, sizeof(int));
  CHECK_INPUT(ret == sizeof(int),
              "CheckPoint::write_header: Error writing precision");

  int ng = mEW->mNumberOfGrids;
  ret = write(fid, &ng, sizeof(int));
  CHECK_INPUT(ret == sizeof(int), "CheckPoint::write_header: Error writing ng");

  ret = write(fid, &a_time, sizeof(float_sw4));
  CHECK_INPUT(ret == sizeof(float_sw4),
              "CheckPoint::write_header: Error writing time");

  ret = write(fid, &a_cycle, sizeof(int));
  CHECK_INPUT(ret == sizeof(int),
              "CheckPoint::write_header: Error writing cycle");

  float_sw4 dt = mEW->getTimeStep();
  ret = write(fid, &dt, sizeof(float_sw4));
  CHECK_INPUT(ret == sizeof(float_sw4),
              "CheckPoint::write_header: Error writing dt");

  int nmech = mEW->getNumberOfMechanisms();
  ret = write(fid, &nmech, sizeof(int));
  CHECK_INPUT(ret == sizeof(int),
              "CheckPoint::write_header: Error writing nmech");
  for (int g = 0; g < ng; g++) {
    int globalSize[6];
    globalSize[0] = 1;
    globalSize[1] = mGlobalDims[g][1] - mGlobalDims[g][0] + 1;
    globalSize[2] = 1;
    globalSize[3] = mGlobalDims[g][3] - mGlobalDims[g][2] + 1;
    globalSize[4] = 1;
    globalSize[5] = mGlobalDims[g][5] - mGlobalDims[g][4] + 1;
    ret = write(fid, globalSize, 6 * sizeof(int));
    CHECK_INPUT(ret == 6 * sizeof(int),
                "CheckPoint::write_header: Error writing global sizes");
    //      cout << "wrote global size " << globalSize[0] << " " <<
    //      globalSize[1] << " " << globalSize[2] << " "
    //	   << globalSize[3] << " " << globalSize[4] << " " << globalSize[5] <<
    //endl;
  }
  hsize = (4 + 6 * ng) * sizeof(int) + 2 * sizeof(float_sw4);
}

//-----------------------------------------------------------------------
void CheckPoint::read_header(int& fid, float_sw4& a_time, int& a_cycle,
                             int& hsize) {
  //
  // Header format: prec - precision 8--> double, 4--> single (int)
  //                ng   - Number of grids (int)
  //                time - Time (float)
  //                cycle- Time step number corresponding to Time (int)
  //                dt   - Time step size (float)
  //                nmech- Number of mechanisms in attenuation model
  //                dims(g,1:6) - Size of array on grid g,
  //                       dim(1) <= i <= dim(2), dim(3) <= j <= dim(4)
  //                       dim(5) <= k <= dim(6)
  //
  int prec;
  size_t ret = read(fid, &prec, sizeof(int));
  CHECK_INPUT(ret == sizeof(int),
              "CheckPoint::read_header: Error reading precision");
  CHECK_INPUT(
      (m_double && prec == 8) || (!m_double && prec == 4),
      "CheckPoint::read_header, floating point precision on restart file"
          << " does not match precision in solver");
  int ng;
  ret = read(fid, &ng, sizeof(int));
  CHECK_INPUT(ret == sizeof(int), "CheckPoint::read_header: Error reading ng");
  CHECK_INPUT(ng == mEW->mNumberOfGrids,
              "CheckPoint::read_header: Error number of grids on restart file"
                  << " does not match number of grids in solver");

  ret = read(fid, &a_time, sizeof(float_sw4));
  CHECK_INPUT(ret == sizeof(float_sw4),
              "CheckPoint::read_header: Error reading time");

  ret = read(fid, &a_cycle, sizeof(int));
  CHECK_INPUT(ret == sizeof(int),
              "CheckPoint::read_header: Error reading cycle");

  float_sw4 dt;
  ret = read(fid, &dt, sizeof(float_sw4));
  CHECK_INPUT(ret == sizeof(float_sw4),
              "CheckPoint::read_header: Error reading dt");

  int nmech;
  ret = read(fid, &nmech, sizeof(int));
  CHECK_INPUT(ret == sizeof(int),
              "CheckPoint::read_header: Error reading nmech");
  CHECK_INPUT(
      nmech == mEW->getNumberOfMechanisms(),
      "CheckPoint::read_header: Error number "
          << "of attenuation mechanisms on restart file"
          << " does not match number of attenuation mechanisms in solver");

  for (int g = 0; g < ng; g++) {
    int globalSize[6];
    ret = read(fid, globalSize, 6 * sizeof(int));
    CHECK_INPUT(ret == 6 * sizeof(int),
                "CheckPoint::read_header: Error reading global sizes");
    CHECK_INPUT(globalSize[0] == 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "low i-index is " << globalSize[0]);
    CHECK_INPUT(globalSize[1] == mGlobalDims[g][1] - mGlobalDims[g][0] + 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "upper i-index is " << globalSize[1]);
    CHECK_INPUT(globalSize[2] == 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "low j-index is " << globalSize[2]);
    CHECK_INPUT(globalSize[3] == mGlobalDims[g][3] - mGlobalDims[g][2] + 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "upper j-index is " << globalSize[3]);
    CHECK_INPUT(globalSize[4] == 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "low k-index is " << globalSize[4]);
    //      CHECK_INPUT( globalSize[5] == mGlobalDims[g][5],
    //      "CheckPoint::read_checkpoint: Error in global sizes, "
    CHECK_INPUT(globalSize[5] == mGlobalDims[g][5] - mGlobalDims[g][4] + 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "upper k-index is " << globalSize[5]);
  }
  hsize = (4 + 6 * ng) * sizeof(int) + 2 * sizeof(float_sw4);
}

//-----------------------------------------------------------------------
void CheckPoint::cycle_checkpoints(string CheckPointFile) {
  // Keep previous check point, remove the second previous.
  // m_fileno is initialized to zero in constructor.
  m_fileno++;
  if (m_fileno >= 3) {
    // Delete  mCheckPointFileMM
    if (m_parallel_io[0]->proc_zero()) {
      int er = unlink(const_cast<char*>(mCheckPointFileMM.c_str()));
      CHECK_INPUT(
          er == 0,
          "CheckPoint::cycle_checkpoints, error deleting old check point file "
              << mCheckPointFileMM);
    }
    mCheckPointFileMM = mCheckPointFileM;
    mCheckPointFileM = CheckPointFile;
  } else if (m_fileno == 1)
    mCheckPointFileMM = CheckPointFile;
  else if (m_fileno == 2)
    mCheckPointFileM = CheckPointFile;
}

//-----------------------------------------------------------------------
void CheckPoint::set_restart_file(string fname, size_t bufsize) {
  mRestartFile = fname;
  m_bufsize = bufsize;
  mDoRestart = true;
}

//-----------------------------------------------------------------------
void CheckPoint::set_restart_path(string restartPath) {
  mRestartPath = restartPath;
  mRestartPathSet = true;
}

//-----------------------------------------------------------------------
std::string CheckPoint::get_restart_path() {
  std::string retval;
  if (mRestartPathSet) {
    retval = mRestartPath;
    return retval;
  }
}

//-----------------------------------------------------------------------
void CheckPoint::set_checkpoint_file(string fname, int cycle, int cycleInterval,
                                     size_t bufsize, bool useHDF5,
                                     int compressionMode,
                                     double compressionPar) {
  mCheckPointFile = fname;
  mWritingCycle = cycle;
  mCycleInterval = cycleInterval;
  m_bufsize = bufsize;
  mDoCheckPointing = true;
  mUseHDF5 = useHDF5;
  mCompMode = compressionMode;
  mCompPar = compressionPar;
}

#ifdef USE_HDF5
//-----------------------------------------------------------------------
void CheckPoint::write_header_hdf5(hid_t fid, float_sw4 a_time, int a_cycle) {
  // Must be called by all participating process
  //
  // Header format: prec - precision 8--> double, 4--> single (int)
  //                ng   - Number of grids (int)
  //                time - Time (float)
  //                cycle- Time step number corresponding to Time (int)
  //                dt   - Time step size (float)
  //                nmech- Number of mechanisms in attenuation model
  //                dims(g,1:6) - Size of array on grid g,
  //                       dim(1) <= i <= dim(2), dim(3) <= j <= dim(4)
  //                       dim(5) <= k <= dim(6)
  //
  int ret;
  hid_t attr, attr_space, attr_space6;
  hsize_t dims = 1, dims6 = 6;

  attr_space = H5Screate_simple(1, &dims, NULL);
  attr_space6 = H5Screate_simple(1, &dims6, NULL);

  hid_t dtype = (m_double ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT);

  int prec = m_double ? 8 : 4;
  attr = H5Acreate(fid, "prec", H5T_NATIVE_INT, attr_space, H5P_DEFAULT,
                   H5P_DEFAULT);
  ret = H5Awrite(attr, H5T_NATIVE_INT, &prec);
  CHECK_INPUT(ret >= 0,
              "CheckPoint::write_header_hdf5: Error writing precision");
  H5Aclose(attr);

  int ng = mEW->mNumberOfGrids;
  attr = H5Acreate(fid, "ngrid", H5T_NATIVE_INT, attr_space, H5P_DEFAULT,
                   H5P_DEFAULT);
  ret = H5Awrite(attr, H5T_NATIVE_INT, &ng);
  CHECK_INPUT(ret >= 0, "CheckPoint::write_header_hdf5: Error writing ng");
  H5Aclose(attr);

  attr = H5Acreate(fid, "time", dtype, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  ret = H5Awrite(attr, dtype, &a_time);
  CHECK_INPUT(ret >= 0, "CheckPoint::write_header_hdf5: Error writing time");
  H5Aclose(attr);

  attr = H5Acreate(fid, "cycle", H5T_NATIVE_INT, attr_space, H5P_DEFAULT,
                   H5P_DEFAULT);
  ret = H5Awrite(attr, H5T_NATIVE_INT, &a_cycle);
  CHECK_INPUT(ret >= 0, "CheckPoint::write_header_hdf5: Error writing cycle");
  H5Aclose(attr);

  float_sw4 dt = mEW->getTimeStep();
  attr = H5Acreate(fid, "dt", dtype, attr_space, H5P_DEFAULT, H5P_DEFAULT);
  ret = H5Awrite(attr, dtype, &dt);
  CHECK_INPUT(ret >= 0, "CheckPoint::write_header_hdf5: Error writing dt");
  H5Aclose(attr);

  int nmech = mEW->getNumberOfMechanisms();
  attr = H5Acreate(fid, "nmesh", H5T_NATIVE_INT, attr_space, H5P_DEFAULT,
                   H5P_DEFAULT);
  ret = H5Awrite(attr, H5T_NATIVE_INT, &nmech);
  CHECK_INPUT(ret >= 0, "CheckPoint::write_header_hdf5: Error writing nmech");
  H5Aclose(attr);

  int nrank;
  MPI_Comm_size(mEW->m_cartesian_communicator, &nrank);
  attr = H5Acreate(fid, "nrank", H5T_NATIVE_INT, attr_space, H5P_DEFAULT,
                   H5P_DEFAULT);
  ret = H5Awrite(attr, H5T_NATIVE_INT, &nrank);
  CHECK_INPUT(ret >= 0, "CheckPoint::write_header_hdf5: Error writing nrank");
  H5Aclose(attr);

  int globalSize[6];
  char name[16];
  for (int g = 0; g < ng; g++) {
    sprintf(name, "mesh%d", g);
    globalSize[0] = 1;
    globalSize[1] = mGlobalDims[g][1] - mGlobalDims[g][0] + 1;
    globalSize[2] = 1;
    globalSize[3] = mGlobalDims[g][3] - mGlobalDims[g][2] + 1;
    globalSize[4] = 1;
    globalSize[5] = mGlobalDims[g][5] - mGlobalDims[g][4] + 1;
    attr = H5Acreate(fid, name, H5T_NATIVE_INT, attr_space6, H5P_DEFAULT,
                     H5P_DEFAULT);
    ret = H5Awrite(attr, H5T_NATIVE_INT, globalSize);
    CHECK_INPUT(ret >= 0,
                "CheckPoint::write_header_hdf5: Error writing global sizes");
    H5Aclose(attr);
  }
  H5Sclose(attr_space);
  H5Sclose(attr_space6);
}

//-----------------------------------------------------------------------
void CheckPoint::read_header_hdf5(hid_t fid, float_sw4& a_time, int& a_cycle) {
  //
  // Header format: prec - precision 8--> double, 4--> single (int)
  //                ng   - Number of grids (int)
  //                time - Time (float)
  //                cycle- Time step number corresponding to Time (int)
  //                dt   - Time step size (float)
  //                nmech- Number of mechanisms in attenuation model
  //                dims(g,1:6) - Size of array on grid g,
  //                       dim(1) <= i <= dim(2), dim(3) <= j <= dim(4)
  //                       dim(5) <= k <= dim(6)
  //
  hid_t attr;
  hid_t dtype = (m_double ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT);
  int ret;

  int prec;
  attr = H5Aopen(fid, "prec", H5P_DEFAULT);
  ret = H5Aread(attr, H5T_NATIVE_INT, &prec);
  H5Aclose(attr);
  CHECK_INPUT(ret >= 0,
              "CheckPoint::read_header_hdf5: Error reading precision");
  CHECK_INPUT(
      (m_double && prec == 8) || (!m_double && prec == 4),
      "CheckPoint::read_header_hdf5, floating point precision on restart file"
          << " does not match precision in solver");
  int ng;
  attr = H5Aopen(fid, "ngrid", H5P_DEFAULT);
  ret = H5Aread(attr, H5T_NATIVE_INT, &ng);
  H5Aclose(attr);
  CHECK_INPUT(ret >= 0, "CheckPoint::read_header_hdf5: Error reading ng");
  CHECK_INPUT(
      ng == mEW->mNumberOfGrids,
      "CheckPoint::read_header_hdf5: Error number of grids on restart file"
          << " does not match number of grids in solver");

  attr = H5Aopen(fid, "time", H5P_DEFAULT);
  ret = H5Aread(attr, dtype, &a_time);
  H5Aclose(attr);
  CHECK_INPUT(ret >= 0, "CheckPoint::read_header_hdf5: Error reading time");

  attr = H5Aopen(fid, "cycle", H5P_DEFAULT);
  ret = H5Aread(attr, H5T_NATIVE_INT, &a_cycle);
  H5Aclose(attr);
  CHECK_INPUT(ret >= 0, "CheckPoint::read_header_hdf5: Error reading cycle");

  /* float_sw4 dt; */
  /* attr = H5Aopen(fid, "dt", H5P_DEFAULT); */
  /* ret = H5Aread(attr, dtype, &dt); */
  /* H5Aclose(attr); */
  /* CHECK_INPUT( ret >= 0,"CheckPoint::read_header_hdf5: Error reading dt" );
   */

  int nmech;
  attr = H5Aopen(fid, "nmesh", H5P_DEFAULT);
  ret = H5Aread(attr, H5T_NATIVE_INT, &nmech);
  H5Aclose(attr);
  CHECK_INPUT(ret >= 0, "CheckPoint::read_header_hdf5: Error reading nmech");
  CHECK_INPUT(
      nmech == mEW->getNumberOfMechanisms(),
      "CheckPoint::read_header_hdf5: Error number "
          << "of attenuation mechanisms on restart file"
          << " does not match number of attenuation mechanisms in solver");

  int nrank, nrank_restart;
  MPI_Comm_size(mEW->m_cartesian_communicator, &nrank);

  attr = H5Aopen(fid, "nrank", H5P_DEFAULT);
  ret = H5Aread(attr, H5T_NATIVE_INT, &nrank_restart);
  H5Aclose(attr);
  CHECK_INPUT(ret >= 0, "CheckPoint::read_header_hdf5: Error reading nrank");
  CHECK_INPUT(nrank == nrank_restart,
              "CheckPoint::read_header_hdf5: Error number "
                  << "of ranks on restart file does not match number of ranks "
                     "in currrent run");

  int globalSize[6];
  char name[16];
  for (int g = 0; g < ng; g++) {
    sprintf(name, "mesh%d", g);
    attr = H5Aopen(fid, name, H5P_DEFAULT);
    ret = H5Aread(attr, H5T_NATIVE_INT, globalSize);
    H5Aclose(attr);
    CHECK_INPUT(ret >= 0,
                "CheckPoint::read_header_hdf5: Error reading global sizes");
    CHECK_INPUT(globalSize[0] == 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "low i-index is " << globalSize[0]);
    CHECK_INPUT(globalSize[1] == mGlobalDims[g][1] - mGlobalDims[g][0] + 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "upper i-index is " << globalSize[1]);
    CHECK_INPUT(globalSize[2] == 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "low j-index is " << globalSize[2]);
    CHECK_INPUT(globalSize[3] == mGlobalDims[g][3] - mGlobalDims[g][2] + 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "upper j-index is " << globalSize[3]);
    CHECK_INPUT(globalSize[4] == 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "low k-index is " << globalSize[4]);
    CHECK_INPUT(globalSize[5] == mGlobalDims[g][5] - mGlobalDims[g][4] + 1,
                "CheckPoint::read_checkpoint: Error in global sizes, "
                    << "upper k-index is " << globalSize[5]);
  }
}

void CheckPoint::finalize_hdf5() {
  size_t num_in_progress;
  hbool_t op_failed;
  int ret;
#ifdef USE_HDF5_ASYNC
  if (m_es_id > 0) {
    ret = H5ESwait(m_es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (ret < 0) fprintf(stderr, "Error with H5ESwait!\n");
    H5ESclose(m_es_id);
    m_es_id = 0;
  }
#endif
  return;
}

void CheckPoint::create_hdf5_dset(hid_t fid, char* dset_name, hid_t dtype,
                                  hid_t dspace, hid_t dcpl) {
  hid_t dset;
#ifdef USE_HDF5_ASYNC
  dset = H5Dcreate_async(fid, dset_name, dtype, dspace, H5P_DEFAULT, dcpl,
                         H5P_DEFAULT, m_es_id);
  H5Dclose_async(dset, m_es_id);
#else
  dset =
      H5Dcreate(fid, dset_name, dtype, dspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  H5Dclose(dset);
#endif
}

void CheckPoint::write_hdf5_dset(hid_t fid, char* dset_name, hid_t dtype,
                                 hid_t mspace, hid_t mydspace, hid_t dxpl,
                                 void* buf) {
  hid_t dset;
#ifdef USE_HDF5_ASYNC
  dset = H5Dopen_async(fid, dset_name, H5P_DEFAULT, m_es_id);
  H5Dwrite_async(dset, dtype, mspace, mydspace, dxpl, buf, m_es_id);
  H5Dclose_async(dset, m_es_id);
#else
  dset = H5Dopen(fid, dset_name, H5P_DEFAULT);
  H5Dwrite(dset, dtype, mspace, mydspace, dxpl, buf);
  H5Dclose(dset);
#endif
}

//-----------------------------------------------------------------------
void CheckPoint::write_checkpoint_hdf5(float_sw4 a_time, int a_cycle,
                                       vector<Sarray>& a_Um,
                                       vector<Sarray>& a_U,
                                       vector<Sarray*>& a_AlphaVEm,
                                       vector<Sarray*>& a_AlphaVE) {
  std::stringstream s;
  std::stringstream fileSuffix;
  compute_file_suffix(a_cycle, fileSuffix);
  if (mRestartPathSet)
    s << mRestartPath << "/";
  else if (mEW->getPath() != "./")
    s << mEW->getPath() << "/";
  s << fileSuffix.str();

  // Keep track of the number of files, save previous file name, and delete the
  // second last.
  cycle_checkpoints(s.str());

  hid_t fid, fapl, dxpl, dspace, mspace, mydspace, dtype, dcpl;
  int myrank, nrank;
  double stime, etime;
  hsize_t my_chunk[1];

  stime = MPI_Wtime();

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  if (myrank == 0)
    std::cout << "Writing checkpoint to file " << s.str() << std::endl;

  fapl = H5Pcreate(H5P_FILE_ACCESS);
#ifdef USE_HDF5_ASYNC
  if (m_es_id > 0) finalize_hdf5();
#endif

  // Each rank writes its data to its own dataset
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

  dtype = (m_double ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT);

  if (mCompMode > 0) {
    my_chunk[0] = (m_double ? 65536 : 4194304);
    char* env_char = NULL;
    env_char = getenv("CHECKPOINT_CHUNK_X");
    if (env_char != NULL) my_chunk[0] = atoi(env_char);

    H5Pset_chunk(dcpl, 1, my_chunk);
    if (myrank == 0) fprintf(stderr, "set chunk size to %llu\n", my_chunk[0]);
  }

  if (mCompMode == SW4_SZIP) {
    H5Pset_szip(dcpl, H5_SZIP_NN_OPTION_MASK, 32);
  } else if (mCompMode == SW4_ZLIB) {
    H5Pset_deflate(dcpl, (int)mCompPar);
  }
#ifdef USE_ZFP
  else if (mCompMode == SW4_ZFP_MODE_RATE) {
    H5Pset_zfp_rate(dcpl, mCompPar);
  } else if (mCompMode == SW4_ZFP_MODE_PRECISION) {
    H5Pset_zfp_precision(dcpl, (unsigned int)mCompPar);
  } else if (mCompMode == SW4_ZFP_MODE_ACCURACY) {
    H5Pset_zfp_accuracy(dcpl, mCompPar);
  } else if (mCompMode == SW4_ZFP_MODE_REVERSIBLE) {
    H5Pset_zfp_reversible(dcpl);
  }
#endif
#ifdef USE_SZ
  else if (mCompMode == SW4_SZ) {
    size_t cd_nelmts;
    unsigned int* cd_values = NULL;
    int dataType = SZ_DOUBLE;
    if (m_precision == 4) dataType = SZ_FLOAT;
    SZ_metaDataToCdArray(&cd_nelmts, &cd_values, dataType, 0, m_cycle_dims[3],
                         m_cycle_dims[2], m_cycle_dims[1], m_cycle_dims[0]);
    H5Pset_filter(dcpl, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts,
                  cd_values);
  }
#endif

  char dset_name[128];
  hsize_t* npts = new hsize_t[mEW->mNumberOfGrids];
  hsize_t* myoff = new hsize_t[mEW->mNumberOfGrids];
  hsize_t* total = new hsize_t[mEW->mNumberOfGrids];

  // Rank 0 create everything and write metadata first
  if (myrank == 0) {
    fid = H5Fcreate(const_cast<char*>(s.str().c_str()), H5F_ACC_TRUNC,
                    H5P_DEFAULT, H5P_DEFAULT);
    CHECK_INPUT(fid > 0, "CheckPoint::write_file: Error opening: " << s.str());

    write_header_hdf5(fid, a_time, a_cycle);
  }

  // Workaround for summit system with spectrum MPI and HDF5 1.10.4,
  // have to create all dsets before write to avoid runtime error
  for (int g = 0; g < mEW->mNumberOfGrids; g++) {
    npts[g] = (hsize_t)(mWindow[g][1] - mWindow[g][0] + 1) *
              (mWindow[g][3] - mWindow[g][2] + 1) *
              (mWindow[g][5] - mWindow[g][4] + 1);
    npts[g] *= 3;
    myoff[g] = 0;
    MPI_Exscan(&npts[g], &myoff[g], 1, MPI_LONG_LONG_INT, MPI_SUM,
               MPI_COMM_WORLD);
    total[g] = npts[g] + myoff[g];
    MPI_Bcast(&total[g], 1, MPI_LONG_LONG_INT, nrank - 1, MPI_COMM_WORLD);
    /* fprintf(stderr, "Create, g %d, rank %d: off %llu, size %llu, total
     * %llu\n", g, myrank, myoff[g], */
    /*         npts[g], total[g]); */
  }

  // Only rank 0 creates the dataset
  if (myrank == 0) {
    for (int g = 0; g < mEW->mNumberOfGrids; g++) {
      dspace = H5Screate_simple(1, &total[g], NULL);

      sprintf(dset_name, "Um%d", g);
      create_hdf5_dset(fid, dset_name, dtype, dspace, dcpl);

      sprintf(dset_name, "U%d", g);
      create_hdf5_dset(fid, dset_name, dtype, dspace, dcpl);

      for (int m = 0; m < mEW->getNumberOfMechanisms(); m++) {
        sprintf(dset_name, "VEM%d_m%d", g, m);
        create_hdf5_dset(fid, dset_name, dtype, dspace, dcpl);

        sprintf(dset_name, "VE%d_m%d", g, m);
        create_hdf5_dset(fid, dset_name, dtype, dspace, dcpl);
      }
      H5Sclose(dspace);
    }
    H5Fflush(fid, H5F_SCOPE_GLOBAL);
    H5Fclose(fid);
  }  // end if myrank=0

  MPI_Barrier(MPI_COMM_WORLD);

  fapl = H5Pcreate(H5P_FILE_ACCESS);

  int alignment = 16777216;
  char* env = getenv("HDF5_ALIGNMENT_SIZE");
  if (env != NULL) alignment = atoi(env);
  if (alignment < 65536) alignment = 65536;

  H5Pset_alignment(fapl, 65536, alignment);
  H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  H5Pset_all_coll_metadata_ops(fapl, 1);
  H5Pset_coll_metadata_write(fapl, 1);

  dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

#ifdef USE_HDF5_ASYNC
  if (m_es_id == 0) m_es_id = H5EScreate();
  fid = H5Fopen_async(const_cast<char*>(s.str().c_str()), H5F_ACC_RDWR, fapl,
                      m_es_id);
#else
  fid = H5Fopen(const_cast<char*>(s.str().c_str()), H5F_ACC_RDWR, fapl);
#endif

  // Start write
  for (int g = 0; g < mEW->mNumberOfGrids; g++) {
    mspace = H5Screate_simple(1, &npts[g], NULL);
    mydspace = H5Screate_simple(1, &total[g], NULL);
    H5Sselect_hyperslab(mydspace, H5S_SELECT_SET, &myoff[g], NULL, &npts[g],
                        NULL);

    // allocate local buffer array
    float_sw4* doubleField = new float_sw4[npts[g]];

    if (m_kji_order)
      a_Um[g].extract_subarrayIK(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                                 mWindow[g][3], mWindow[g][4], mWindow[g][5],
                                 doubleField);
    else
      a_Um[g].extract_subarray(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                               mWindow[g][3], mWindow[g][4], mWindow[g][5],
                               doubleField);

    sprintf(dset_name, "Um%d", g);
    write_hdf5_dset(fid, dset_name, dtype, mspace, mydspace, dxpl, doubleField);

    if (m_kji_order)
      a_U[g].extract_subarrayIK(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                                mWindow[g][3], mWindow[g][4], mWindow[g][5],
                                doubleField);
    else
      a_U[g].extract_subarray(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                              mWindow[g][3], mWindow[g][4], mWindow[g][5],
                              doubleField);
    sprintf(dset_name, "U%d", g);
    write_hdf5_dset(fid, dset_name, dtype, mspace, mydspace, dxpl, doubleField);

    for (int m = 0; m < mEW->getNumberOfMechanisms(); m++) {
      if (m_kji_order)
        a_AlphaVEm[g][m].extract_subarrayIK(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
      else
        a_AlphaVEm[g][m].extract_subarray(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
      sprintf(dset_name, "VEM%d_m%d", g, m);
      write_hdf5_dset(fid, dset_name, dtype, mspace, mydspace, dxpl,
                      doubleField);

      if (m_kji_order)
        a_AlphaVE[g][m].extract_subarrayIK(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
      else
        a_AlphaVE[g][m].extract_subarray(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
      sprintf(dset_name, "VE%d_m%d", g, m);
      write_hdf5_dset(fid, dset_name, dtype, mspace, mydspace, dxpl,
                      doubleField);
    }
    H5Sclose(mspace);
    H5Sclose(mydspace);

    delete[] doubleField;
  }

  delete[] npts;
  delete[] myoff;
  delete[] total;

  H5Pclose(dcpl);
  H5Pclose(fapl);
  H5Pclose(dxpl);

#ifdef USE_HDF5_ASYNC
  H5Fclose_async(fid, m_es_id);
#else
  H5Fclose(fid);
#endif
  etime = MPI_Wtime();

  if (myrank == 0)
    std::cout << "Written checkpoint, " << etime - stime << " seconds"
              << std::endl;

}  // end write_checkpoint_hdf5()

//-----------------------------------------------------------------------
void CheckPoint::read_checkpoint_hdf5(float_sw4& a_time, int& a_cycle,
                                      vector<Sarray>& a_Um, vector<Sarray>& a_U,
                                      vector<Sarray*>& a_AlphaVEm,
                                      vector<Sarray*>& a_AlphaVE) {
  int myrank, nrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nrank);

  std::stringstream s;
  if (mRestartPathSet)
    s << mRestartPath << "/";
  else if (mEW->getPath() != "./")
    s << mEW->getPath();
  s << mRestartFile;

  if (myrank == 0)
    std::cout << "Reading checkpoint from file " << s.str() << std::endl;

  // Open file from processor zero and read header.
  hid_t fid, fapl, dset, dxpl, dspace, mspace, mydspace;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  H5Pset_all_coll_metadata_ops(fapl, 1);
  H5Pset_coll_metadata_write(fapl, 1);
  H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

  dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  fid = H5Fopen(const_cast<char*>(s.str().c_str()), H5F_ACC_RDONLY, fapl);
  CHECK_INPUT(fid > 0,
              "CheckPoint::read_checkpoint_hdf5: Error opening: " << s.str());

  read_header_hdf5(fid, a_time, a_cycle);

  char dset_name[128];
  hsize_t npts, myoff, total;
  hid_t dtype = (m_double ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT);

  for (int g = 0; g < mEW->mNumberOfGrids; g++) {
    npts = (hsize_t)(mWindow[g][1] - mWindow[g][0] + 1) *
           (mWindow[g][3] - mWindow[g][2] + 1) *
           (mWindow[g][5] - mWindow[g][4] + 1);
    npts *= 3;
    myoff = 0;
    MPI_Exscan(&npts, &myoff, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    total = npts + myoff;
    MPI_Bcast(&total, 1, MPI_LONG_LONG_INT, nrank - 1, MPI_COMM_WORLD);
    /* fprintf(stderr, "rank %d: off %llu, size %llu, total %llu\n", myrank,
     * myoff, npts, total); */

    dspace = H5Screate_simple(1, &total, NULL);
    mspace = H5Screate_simple(1, &npts, NULL);
    mydspace = H5Screate_simple(1, &total, NULL);
    H5Sselect_hyperslab(mydspace, H5S_SELECT_SET, &myoff, NULL, &npts, NULL);

    // allocate local buffer array
    float_sw4* doubleField = new float_sw4[3 * npts];

    sprintf(dset_name, "Um%d", g);
    dset = H5Dopen(fid, dset_name, H5P_DEFAULT);
    H5Dread(dset, dtype, mspace, mydspace, dxpl, doubleField);
    H5Dclose(dset);
    if (m_kji_order)
      a_Um[g].insert_subarrayIK(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                                mWindow[g][3], mWindow[g][4], mWindow[g][5],
                                doubleField);
    else
      a_Um[g].insert_subarray(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                              mWindow[g][3], mWindow[g][4], mWindow[g][5],
                              doubleField);

    sprintf(dset_name, "U%d", g);
    dset = H5Dopen(fid, dset_name, H5P_DEFAULT);
    H5Dread(dset, dtype, mspace, mydspace, dxpl, doubleField);
    H5Dclose(dset);
    if (m_kji_order)
      a_U[g].insert_subarrayIK(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                               mWindow[g][3], mWindow[g][4], mWindow[g][5],
                               doubleField);
    else
      a_U[g].insert_subarray(mWindow[g][0], mWindow[g][1], mWindow[g][2],
                             mWindow[g][3], mWindow[g][4], mWindow[g][5],
                             doubleField);

    for (int m = 0; m < mEW->getNumberOfMechanisms(); m++) {
      sprintf(dset_name, "VEM%d_m%d", g, m);
      dset = H5Dopen(fid, dset_name, H5P_DEFAULT);
      H5Dread(dset, dtype, mspace, mydspace, dxpl, doubleField);
      H5Dclose(dset);
      if (m_kji_order)
        a_AlphaVEm[g][m].insert_subarrayIK(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
      else
        a_AlphaVEm[g][m].insert_subarray(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);

      sprintf(dset_name, "VE%d_m%d", g, m);
      dset = H5Dopen(fid, dset_name, H5P_DEFAULT);
      H5Dread(dset, dtype, mspace, mydspace, dxpl, doubleField);
      H5Dclose(dset);
      if (m_kji_order)
        a_AlphaVE[g][m].insert_subarrayIK(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
      else
        a_AlphaVE[g][m].insert_subarray(
            mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3],
            mWindow[g][4], mWindow[g][5], doubleField);
    }
    H5Sclose(dspace);
    H5Sclose(mspace);
    H5Sclose(mydspace);

    delete[] doubleField;
  }

  H5Pclose(fapl);
  H5Pclose(dxpl);
  H5Fclose(fid);
}
#endif  // End USE_HDF5
