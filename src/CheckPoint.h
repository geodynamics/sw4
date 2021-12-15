#ifndef SW4_CHECKPOINT_H
#define SW4_CHECKPOINT_H

#include <sstream>
#include <string>
#include <vector>

#include "sw4.h"
//#include "boundaryConditionTypes.h"
#include "Parallel_IO.h"
#include "Sarray.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

class EW;

class CheckPoint {
 public:
  static CheckPoint* nil;  // nil pointer
  CheckPoint(EW* a_ew);
  CheckPoint(EW* a_ew, int cycle, int cycleInterval, string fname,
             size_t bufsize = 10000000);
  CheckPoint(EW* a_ew, string fname, size_t bufsize = 10000000);
  ~CheckPoint();
  void set_restart_file(string fname, size_t bufsize);
  void set_checkpoint_file(string fname, int cycle, int cycleInterval,
                           size_t bufsize, bool useHDF5, int compressionMode,
                           double compressionPar);

  void write_checkpoint(float_sw4 a_time, int a_cycle,
                        std::vector<Sarray>& a_Um, std::vector<Sarray>& a_U,
                        std::vector<Sarray*>& a_AlphaVEm,
                        std::vector<Sarray*>& a_AlphaVE);
  void read_checkpoint(float_sw4& a_time, int& a_cycle,
                       std::vector<Sarray>& a_Um, std::vector<Sarray>& a_U,
                       std::vector<Sarray*>& a_AlphaVEm,
                       std::vector<Sarray*>& a_AlphaVE);
#ifdef USE_HDF5
  void write_checkpoint_hdf5(float_sw4 a_time, int a_cycle,
                             std::vector<Sarray>& a_Um,
                             std::vector<Sarray>& a_U,
                             std::vector<Sarray*>& a_AlphaVEm,
                             std::vector<Sarray*>& a_AlphaVE);
  void read_checkpoint_hdf5(float_sw4& a_time, int& a_cycle,
                            std::vector<Sarray>& a_Um, std::vector<Sarray>& a_U,
                            std::vector<Sarray*>& a_AlphaVEm,
                            std::vector<Sarray*>& a_AlphaVE);
  void create_hdf5_dset(hid_t fid, char* dset_name, hid_t dtype, hid_t dspace,
                        hid_t dcpl);
  void write_hdf5_dset(hid_t fid, char* dset_name, hid_t dtype, hid_t mspace,
                       hid_t mydspace, hid_t dxpl, void* buf);
  void finalize_hdf5();
#endif
  void setup_sizes();
  bool timeToWrite(float_sw4 time, int cycle, float_sw4 dt);
  float_sw4 getDt();
  bool do_checkpointing();
  int get_checkpoint_cycle_interval();
  bool do_restart();
  void set_restart_path(string restartPath);
  std::string get_restart_path();
  bool useHDF5() { return mUseHDF5; };

 protected:
  void define_pio();
  void setSteps(int a_steps);
  void compute_file_suffix(int cycle, std::stringstream& fileSuffix);
  void cycle_checkpoints(string fname);
  void write_header(int& fid, float_sw4 a_time, int a_cycle, int& hsize);
  void read_header(int& fid, float_sw4& a_time, int& a_cycle, int& hsize);
#ifdef USE_HDF5
  void write_header_hdf5(hid_t fid, float_sw4 a_time, int a_cycle);
  void read_header_hdf5(hid_t fid, float_sw4& a_time, int& a_cycle);
  hid_t m_es_id;
#endif

  std::string mCheckPointFile;
  std::string mRestartFile;

  // keep track of previous check point files
  std::string mCheckPointFileM;
  std::string mCheckPointFileMM;
  int m_fileno;

  float_sw4 mStartTime;
  int mWritingCycle;
  int mCycleInterval;
  bool m_winallocated;
  size_t m_bufsize;

  bool mDoCheckPointing;
  bool mDoRestart;
  bool mRestartPathSet;
  string mRestartPath;
  bool mUseHDF5;
  int mCompMode;
  double mCompPar;

 private:
  CheckPoint();  // make it impossible to call default constructor
  CheckPoint(const CheckPoint& cp);  // hide copy constructor
  int mPreceedZeros;  // number of digits for unique time step in file names
  bool m_double;
  bool m_kji_order;
  EW* mEW;
  Parallel_IO** m_parallel_io;
  std::vector<int*> mWindow;      // Local in processor start + end indices for
                                  // (i,j,k) for each grid level
  std::vector<int*> mGlobalDims;  // Global start + end indices for (i,j,k) for
                                  // each grid level
  std::vector<bool> m_ihavearray;
};

#endif
