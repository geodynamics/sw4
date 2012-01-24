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
   const char* basedir = EW_BASEDIR;
   const char* version = "0.0.1";
   
   std::string getVersionInfo()
   {
      std::stringstream versioninfo;
      
      versioninfo << "----------------------------------------------------------------" << endl
                  << "            sbp4f version " << version    << endl
		  << endl
		  << " This program comes with ABSOLUTELY NO WARRANTY; released under GPL." << endl
		  << " This is free software, and you are welcome to redistribute     " << endl
		  << " it under certain conditions, see LICENSE.txt for more details  " << endl
                  << "----------------------------------------------------------------" << endl
                  << "  Compiled: " << when << std::endl
                  << "  By:       " << madeby << std::endl
                  << "  Machine:  " << hostname << std::endl
                  << "  Compiler: " << compiler << std::endl
                  << "  3rd party software base directory: " << basedir << std::endl
                  << "----------------------------------------------------------------" << std::endl;
      return versioninfo.str();
   }
};
