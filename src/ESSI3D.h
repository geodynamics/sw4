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

class EW;
class ESSI3DHDF5;

class ESSI3D
{
public:
   static ESSI3D* nil;

   ESSI3D( EW * a_ew,
	    const std::string& filePrefix,
	    bool doubleMode,
      int cycleInterval,
 	    float_sw4 coordBox[4],
	    float_sw4 depth );
   ~ESSI3D();

   void setup( );

   static void setSteps(int a_steps);

   void update_image( int a_cycle, float_sw4 a_time, float_sw4 a_dt,
       std::vector<Sarray>& a_U, std::string a_path, Sarray& a_Z );

   void force_write_image( float_sw4 a_time, int a_cycle, vector<Sarray>& a_U,
			   std::string a_path, Sarray& a_Z );

   //   void set_start_time(double tStart);

protected:
   void compute_image( Sarray& a_U, int a_comp );

   void write_image( int cycle, std::string &path, float_sw4 t, Sarray& a_Z );

   void define_pio( );

#ifdef USE_HDF5
   void write_image_hdf5( int cycle, std::string &path, float_sw4 t,
      Sarray& a_Z, vector<Sarray>& a_U );
   void define_pio_hdf5( );
#endif

   void compute_file_suffix( int cycle, std::stringstream & fileSuffix );

   std::string mFilePrefix;
   float_sw4 mTime;
   float_sw4 mCoordBox[4];
   float_sw4 mDepth;

   std::string mFileName;

   bool m_isDefinedMPIWriters;
   bool m_memallocated;
   
   bool m_fileOpen;

   static int mPreceedZeros; // number of digits for unique time step in file names

private:
   ESSI3D(); // make it impossible to call default constructor
   ESSI3D(const ESSI3D &in); // hide copy constructor

   bool m_double;
   EW* mEW;

   ESSI3DHDF5* m_hdf5helper;

   int m_cycle;
   int m_cycleInterval; // Note: this cycle interval to start a new file
   int mWindow[6]; // Local in processor start + end indices for (i,j,k) for last curvilinear grid
   int mGlobalDims[6]; // Global start + end indices for (i,j,k) for last curvilinear grid
   double* m_doubleField;
   float* m_floatField;
   bool m_ihavearray;
};

#endif
