#include <fstream>
using namespace std;
//
// Generate sac commands for WPP
// addoption = 1 : Input is WPP input file, the program adds new
//               sac keywords to each sac-line
// genfile = 1 : Input is list of three columns, (lat,lon,station)
//                the program generates a sac commands for each line
//                in the input file.
int main( int argc, char** argv )
{
   int addoption=0,genfile=1;
   ifstream infil(argv[1]);
   ofstream utfil("sacstations.in");
   char* buf = new char[500];
   if( addoption == 1 )
   {
      while( !infil.eof() )
      {
	 infil.getline(buf,500);
	 string line = buf;
	 if( buf[0] == 's' && buf[1] =='a' && buf[2] == 'c' )
	 {
	    line += " nsew=1 velocity=1 usgsformat=1 sacformat=0";
	 }
	 utfil << line.c_str() << endl;
      }
   }
   if( genfile == 1 )
   {
      while( !infil.eof() )
      {
         double lat, lon;
	 string station;
         infil >> lon >> lat >> station;
	 //	 infil.getline(buf,500);
	 //	 string line = buf;
         if( !infil.eof() )
         utfil << "sac lon="<<lon<< " lat="<<lat <<" depth=0.0 sta=" <<station
	       << " file="<<station << " nsew=1 velocity=1 usgsformat=1 sacformat=0" << endl;
      }
   }
}
