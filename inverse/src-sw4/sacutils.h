#ifndef SW4_SACUTILS
#define SW4_SACUTILS

int lastofmonth( int year, int month );

void convertjday( int jday, int year, int& day, int& month );

void readSACheader( const char* fname, float_sw4& dt, float_sw4& t0,
		    float_sw4& lat, float_sw4& lon, float_sw4& cmpaz,
		    float_sw4& cmpinc, int utc[7], int& npts,
		    bool& need_byte_reversal );

void readSACdata( const char* fname, int npts, float_sw4* u, bool need_byte_reversal=false );

#endif
