#ifndef SW4_SACUTILS
#define SW4_SACUTILS

sw4_type lastofmonth(sw4_type year, sw4_type month);

void convertjday(sw4_type jday, sw4_type year, sw4_type& day, sw4_type& month);

void readSACheader(const char* fname, double& dt, double& t0, double& lat,
                   double& lon, double& cmpaz, double& cmpinc, int utc[7],
                   int& npts, bool& need_byte_reversal);

void readSACdata(const char* fname, sw4_type npts, double* u,
                 bool need_byte_reversal = false);

#endif
