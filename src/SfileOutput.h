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
#ifndef SW4_SFILEOUTPUT_H
#define SW4_SFILEOUTPUT_H

#include <sstream>
#include <string>
#include <vector>

#include "Parallel_IO.h"
#include "Sarray.h"
#include "boundaryConditionTypes.h"
#include "sw4.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

class EW;

class SfileOutput {
 public:
  enum SfileOutputMode {
    NONE,
    UX,
    UY,
    UZ,
    RHO,
    LAMBDA,
    MU,
    P,
    S,
    GRADRHO,
    GRADMU,
    GRADLAMBDA,
    GRADP,
    GRADS,
    QP,
    QS
  };
  // 15 modes are currently defined

  static SfileOutput* nil;

  SfileOutput(EW* a_ew, float_sw4 time, float_sw4 timeSw4_Typeerval, sw4_type cycle,
              sw4_type cycleSw4_Typeerval, float_sw4 tstart,
              const std::string& filePrefix, sw4_type sampleFactor_h,
              sw4_type sampleFactor_v, bool doubleMode);
  ~SfileOutput();

  void setup_images();

  static void setSteps(sw4_type a_steps);

  void update_image(sw4_type a_cycle, float_sw4 a_time, float_sw4 a_dt,
                    std::vector<Sarray>& a_U, std::vector<Sarray>& a_Rho,
                    std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda,
                    std::vector<Sarray>& a_gRho, std::vector<Sarray>& a_gMu,
                    std::vector<Sarray>& a_gLambda, std::vector<Sarray>& a_Qp,
                    std::vector<Sarray>& a_Qs, std::string a_path,
                    std::vector<Sarray>& a_Z);

  void force_write_image(float_sw4 a_time, sw4_type a_cycle,
                         std::vector<Sarray>& a_U, std::vector<Sarray>& a_Rho,
                         std::vector<Sarray>& a_Mu,
                         std::vector<Sarray>& a_Lambda,
                         std::vector<Sarray>& a_gRho,
                         std::vector<Sarray>& a_gMu,
                         std::vector<Sarray>& a_gLambda,
                         std::vector<Sarray>& a_Qp, std::vector<Sarray>& a_Qs,
                         std::string a_path, std::vector<Sarray>& a_Z);

 protected:
  bool timeToWrite(float_sw4 time, sw4_type cycle, float_sw4 dt);

  void compute_image(std::vector<Sarray>& a_U, std::vector<Sarray>& a_Rho,
                     std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda,
                     std::vector<Sarray>& a_gRho, std::vector<Sarray>& a_gMu,
                     std::vector<Sarray>& a_gLambda, std::vector<Sarray>& a_Qp,
                     std::vector<Sarray>& a_Qs, std::vector<Sarray>& a_Z);

  void write_image(const char* fname, std::vector<Sarray>& a_Z);

  void define_pio();

  void compute_file_suffix(sw4_type cycle, std::stringstream& fileSuffix);

  void gen_fname(std::string& path, sw4_type cycle, std::string& fname);

  SfileOutputMode mMode;
  std::string mFilePrefix;
  float_sw4 mTime;
  bool m_time_done;
  float_sw4 mTimeSw4_Typeerval;
  float_sw4 mNextTime;
  float_sw4 mStartTime;

  sw4_type mWritingCycle;
  sw4_type mCycleSw4_Typeerval;
  /* sw4_type mImageSamplingFactor; */
  sw4_type mSampleH;
  sw4_type mSampleV;

  std::string mFileName;
  /* std::string m_modestring; */

  bool m_isDefinedMPIWriters;
  bool m_winallocated;
  bool m_memallocated;

  static sw4_type
      mPreceedZeros;  // number of digits for unique time step in file names

 private:
  SfileOutput();  // make it impossible to call default constructor
  SfileOutput(const Image& im);  // hide copy constructor

  bool m_double;
  bool m_att;
  EW* mEW;
  Parallel_IO** m_parallel_io;

  std::vector<sw4_type*> mWindow;      // Local in processor start + end indices for
                                  // (i,j,k) for each grid level
  std::vector<sw4_type*> mGlobalDims;  // Global start + end indices for (i,j,k) for
                                  // each grid level

  std::vector<double*> m_doubleField;
  std::vector<float*> m_floatField;
  std::vector<sw4_type> m_extraz;
  std::vector<bool> m_ihavearray;
  bool m_isCreated;
};

#endif
