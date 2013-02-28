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
   const char* version = "1.0";
   
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
