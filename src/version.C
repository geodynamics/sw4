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
#include "version.h"
#include <sstream>

using namespace std;

namespace ewversion 
{

   const char* madeby = EW_MADEBY;
   const char* when = EW_WHEN;
   const char* hostname = EW_HOSTNAME;
   const char* optimization = EW_OPT_LEVEL;
   const char* compiler = EW_COMPILER;
//   const char* basedir = EW_BASEDIR;
   const char* libdir = EW_LIBDIR;
   const char* incdir = EW_INCDIR;
   const char* version = "1.18";
   
   std::string getVersionInfo()
   {
      std::stringstream versioninfo;
      
      versioninfo << "----------------------------------------------------------------" << endl
                  << "            sw4 version " << version    << endl
		  << endl
		  << " This program comes with ABSOLUTELY NO WARRANTY; released under GPL." << endl
		  << " This is free software, and you are welcome to redistribute     " << endl
		  << " it under certain conditions, see LICENSE.txt for more details  " << endl
                  << "----------------------------------------------------------------" << endl
                  << "  Compiled on: " << when << std::endl
                  << "  By user:     " << madeby << std::endl
                  << "  Machine:     " << hostname << std::endl
                  << "  Compiler:    " << compiler << std::endl
                  << "  3rd party include dir: " << incdir << ", and library dir: " << libdir << std::endl
                  << "----------------------------------------------------------------" << std::endl;
      return versioninfo.str();
   }
};
