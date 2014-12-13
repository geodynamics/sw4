#include <iostream>
#include <cstdio>

#include "Require.h"
#include "Byteswapper.h"

using namespace std;
//-----------------------------------------------------------------------
int lastofmonth( int year, int month )
{
   int days;
   int leapyear=0;
   leapyear = ( year % 400 == 0 ) || ( (year % 4 == 0) && !(year % 100 == 0) );
   if( month == 2 )
      days = 28 + leapyear;
   else if( month==4 || month==6 || month==9 || month==11 )
      days = 30;
   else
      days = 31;
   return days;
}

//-----------------------------------------------------------------------
void convertjday( int jday, int year, int& day, int& month )
{
   if( jday > 0 && jday < 367 )
   {
      day = 1;
      int jd  = 1;
      month = 1;
      while( jd < jday )
      {
	 jd++;
         day++;
	 if( day > lastofmonth(year,month) )
	 {
	    day = 1;
	    month++;
	 }
      }
   }
   else
      cout << "convertjday: ERROR, jday is outside range " << endl;
}

//-----------------------------------------------------------------------
void readSACheader( const char* fname, double& dt, double& t0,
		    double& lat, double& lon, double& cmpaz,
		    double& cmpinc, int utc[7], int& npts,
		    bool& need_byte_reversal )
{
   float float70[70];
   int int35[35], logical[5];
   char kvalues[192];
   const int maxversion = 6;

   need_byte_reversal = false;
   if( !(sizeof(float)==4) || !(sizeof(int)==4) || !(sizeof(char)==1) )
   {
      cout << "readSACheader: ERROR, size of datatypes do not match the SAC specification. Can not read SAC file "
	   << fname << endl;
      return;
   }

// Open SAC file
   FILE* fd=fopen(fname,"r");
   if( fd == NULL )
   {
      cout << "readSACheader: ERROR, observed data file " << fname << " could not be opened" << endl;
      return;
   }

// Read header data blocks
   size_t nr = fread(float70, sizeof(float), 70, fd );
   if( nr != 70 )
   {
      cout << "readSACheader: ERROR, could not read float part of header of " << fname << endl;
      fclose(fd);
      return;
   }
   nr = fread(int35, sizeof(int), 35, fd );
   if( nr != 35 )
   {
      cout << "readSACheader: ERROR, could not read int part of header of " << fname << endl;
      fclose(fd);
      return;
   }
   nr = fread(logical, sizeof(int), 5, fd );   
   if( nr != 5 )
   {
      cout << "readSACheader: ERROR, could not read bool part of header of " << fname << endl;
      fclose(fd);
      return;
   }
   nr = fread(kvalues, sizeof(char), 192, fd );   
   if( nr != 192 )
   {
      cout << "readSACheader: ERROR, could not read character part of header of " << fname << endl;
      fclose(fd);
      return;
   }
   int nvhdr = int35[6];
   if( nvhdr <0 || nvhdr > maxversion )
   {
      // Something is wrong, try byte reversal
      int nvhdrold = nvhdr;
      Byteswapper bswap;
      bswap.byte_rev( &nvhdr, 1, "int");
      if( nvhdr >= 0 && nvhdr <= maxversion )
      {
	 need_byte_reversal = true;
	 bswap.byte_rev( float70, 70, "float");
	 bswap.byte_rev( int35,   35, "int"  );
	 bswap.byte_rev( logical, 5,  "int"  );
      }
      else
	 CHECK_INPUT(false,"readSACheader: ERROR reading sac header of file " << fname << " header version = "
		 << nvhdrold << " byte swapped header version = " << nvhdr );
   }

// Take out wanted information
   dt     = float70[0]; 
   t0     = float70[5];
   lat    = float70[31];
   lon    = float70[32];
   cmpaz  = float70[57];
   cmpinc = float70[58];
   utc[0] = int35[0];
   int jday=int35[1];
   if( utc[0] != -12345 )
      convertjday( jday, utc[0], utc[2], utc[1] );
   utc[3] = int35[2];
   utc[4] = int35[3];
   utc[5] = int35[4];
   utc[6] = int35[5];
   npts   = int35[9];
}

//-----------------------------------------------------------------------
void readSACdata( const char* fname, int npts, double* u, bool need_byte_reversal )
{
   if( !(sizeof(float)==4) || !(sizeof(int)==4) || !(sizeof(char)==1) )
   {
      cout << "readSACdata: ERROR, size of datatypes do not match the SAC specification. Can not read SAC file "
	   << fname << endl;
      return;
   }

// Open SAC file
   FILE* fd=fopen(fname,"r");
   if( fd == NULL )
   {
      cout << "readSACdata: ERROR, observed data file " << fname << " could not be opened" << endl;
      return;
   }

// Skip header
   int retcode = fseek(fd,158*4,SEEK_SET);
   if( retcode != 0 )
   {
      cout << "readSACdata: ERROR, could not skip header in file " << fname << endl;
      fclose(fd);
      return;
   }

// Read data
   float* uf = new float[npts];
   size_t nr = fread( uf, sizeof(float), npts, fd );
   if( nr != npts )
   {
      cout << "readSACdata: ERROR, could not read float array of " << fname << endl;
      delete[] uf;
      fclose(fd);
      return;
   }
   if( need_byte_reversal )
   {
      Byteswapper bswap;
      bswap.byte_rev( uf, npts, "float" );
   }

// Return floats as doubles
   for( int i=0 ; i < npts ; i++ )
      u[i] = static_cast<double>(uf[i]);
   delete[] uf;
   fclose(fd);
}



