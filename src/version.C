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
#include <sstream>

#include "version.h"

using namespace std;
std::string compiler_options();

namespace ewversion {

const char* madeby = EW_MADEBY;
const char* when = EW_WHEN;
const char* hostname = EW_HOSTNAME;
const char* optimization = EW_OPT_LEVEL;
const char* compiler = EW_COMPILER;
//   const char* basedir = EW_BASEDIR;
const char* libdir = EW_LIBDIR;
const char* incdir = EW_INCDIR;
const char* version = "2.1-alpha";

std::string getVersionInfo() {
  std::stringstream versioninfo;

  versioninfo
      << "----------------------------------------------------------------"
      << endl
      << "            sw4 version " << version << endl
      << endl
      << " This program comes with ABSOLUTELY NO WARRANTY; released under GPL."
      << endl
      << " This is free software, and you are welcome to redistribute     "
      << endl
      << " it under certain conditions, see LICENSE.txt for more details  "
      << endl
      << "----------------------------------------------------------------"
      << endl
      << "  Compiled on: " << when << std::endl
      << "  By user:     " << madeby << std::endl
      << "  Machine:     " << hostname << std::endl
      << "  Compiler:    " << compiler << std::endl
      << "  3rd party include dir: " << incdir
      << ", and library dir: " << libdir << std::endl
      << "  Options:     " << compiler_options()<<std::endl
      << "----------------------------------------------------------------"
      << std::endl;
  return versioninfo.str();
}
};  // namespace ewversion




// Print library versions and compiler options
// ! Indicates options that slow down the code
//
#include "RAJA/RAJA.hpp"
#if defined(SW4_USE_UMPIRE)
#include "umpire/config.hpp"
#endif
std::string compiler_options(){

  
  std::stringstream opts;

  opts<<"RAJA("<<RAJA_VERSION_MAJOR<<"."<<RAJA_VERSION_MINOR<<"."<<RAJA_VERSION_PATCHLEVEL<<"):  ";

#if defined(SW4_USE_UMPIRE)
  opts<<"UMPIRE("<<UMPIRE_VERSION_MAJOR<<"."<<UMPIRE_VERSION_MINOR<<"."<<UMPIRE_VERSION_PATCH<<"):\n\t\t";
#else
  opts<<"NO_UMPIRE(!):";
#endif


#if defined(ENABLE_PROJ4)
  opts<<"PROJ4: ";
#endif

#if defined(RAJA_ONLY)
  opts<<"RAJA_ONLY(!): ";
#endif

#if defined(_OPENMP)
  opts<<"OPENMP: ";
#endif

#if defined(SW4_STAGED_MPI_BUFFERS)
  opts<<"STAGED_MPI_BUFFERS: ";
#endif

#if defined(SW4_MANAGED_MPI_BUFFERS)
  opts<<"MANAGED_MPI_BUFFERS: ";
#endif

#if defined(SW4_PINNED_MPI_BUFFERS)
  opts<<"PINNED_MPI_BUFFERS: ";
#endif

#if defined(SW4_DEVICE_MPI_BUFFERS)
  opts<<"DEVICE_MPI_BUFFERS: ";
#endif

#if defined(USE_HDF5)
  opts<<"HDF5: ";
#endif

#if defined(ENABLE_CALIPER)
  opts<<"CALIPER(!): ";
#endif

#if defined(ENABLE_FFTW)
  opts<<"FFTW: ";
#endif
  return opts.str().c_str();
}
