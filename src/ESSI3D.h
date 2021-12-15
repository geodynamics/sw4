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
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
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
#ifndef SW4_ESSI3D_H
#define SW4_ESSI3D_H

#include <string>
#include <sstream>
#include <vector>

#include "sw4.h"
#include "boundaryConditionTypes.h"
#include "Sarray.h"
#include "Parallel_IO.h"

#define SW4_ZFP_MODE_RATE       1
#define SW4_ZFP_MODE_PRECISION  2
#define SW4_ZFP_MODE_ACCURACY   3
#define SW4_ZFP_MODE_REVERSIBLE 4
#define SW4_SZIP                5
#define SW4_ZLIB                6
#define SW4_SZ                  7

class EW;
class ESSI3DHDF5;

class ESSI3D
{
public:
   static ESSI3D* nil;

   ESSI3D( EW * a_ew,
      const std::string& filePrefix,
      int dumpInterval,
      int bufferInterval,
      float_sw4 coordBox[4],
      float_sw4 depth,
      int precision,
      int compressionMode,
      double compressionPar
      );
   ~ESSI3D();

   void set_dump_interval( int a_dumpInterval );
   void set_buffer_interval( int a_bufferInterval );
   void setup( );

   double getHDF5Timings();
   void set_ntimestep(int ntimestep);
   void set_restart(bool is_restart);

   static void setSteps(int a_steps);

   void update_image( int a_cycle, float_sw4 a_time, float_sw4 a_dt,
       std::vector<Sarray>& a_U, std::string& a_path, Sarray& a_Z );

   void force_write_image( float_sw4 a_time, int a_cycle, std::vector<Sarray>& a_U,
       std::string& a_path, Sarray& a_Z );

   void finalize_hdf5();

protected:
   void compute_image( Sarray& a_U, int a_comp, int cycle );

   void write_image( int cycle, std::string& path, float_sw4 t, Sarray& a_Z );

   void define_pio( );

#ifdef USE_HDF5
   void open_vel_file(int a_cycle, std::string& a_path, float_sw4 a_time,
       Sarray& a_Z);
   void close_vel_file();
   void write_image_hdf5( int cycle, std::string& path, float_sw4 t,
       std::vector<Sarray>& a_U );
#endif

   void compute_file_suffix( int cycle, std::stringstream & fileSuffix );

   std::string mFilePrefix;
   float_sw4 mTime;
   float_sw4 mCoordBox[4];
   float_sw4 mDepth;
   int m_precision;

   std::string mFileName;

   bool m_isDefinedMPIWriters;
   bool m_memallocated;

   double m_hdf5_time;

   bool m_fileOpen;

   int m_compressionMode;
   double m_compressionPar;

   static int mPreceedZeros; // number of digits for unique time step in file names
   static int mNumberOfTimeSteps; // number of time steps for the whole sim

private:
   ESSI3D(); // make it impossible to call default constructor
   ESSI3D(const ESSI3D &in); // hide copy constructor

   EW* mEW;

   ESSI3DHDF5* m_hdf5helper;

   int m_cycle;
   int m_dumpInterval; // Note: this cycle interval to start a new file
   int m_bufferInterval; // Note: number of steps to buffer data before writting out
   int m_nbufstep;
   int mWindow[6]; // Local in processor start + end indices for (i,j,k) for last curvilinear grid
   int mGlobalDims[6]; // Global start + end indices for (i,j,k) for last curvilinear grid
   double** m_doubleField;
   float** m_floatField;
   bool m_ihavearray;
   int m_ntimestep;
   bool m_isRestart;
   int m_rank;
};

#endif
