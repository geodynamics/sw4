#ifndef SW4_SACUTILS
#define SW4_SACUTILS

int lastofmonth( int year, int month );

void convertjday( int jday, int year, int& day, int& month );

void readSACheader( const char* fname, double& dt, double& t0,
		    double& lat, double& lon, double& cmpaz,
		    double& cmpinc, int utc[7], int& npts,
		    bool& need_byte_reversal );

void readSACdata( const char* fname, int npts, double* u, bool need_byte_reversal=false );

#endif
