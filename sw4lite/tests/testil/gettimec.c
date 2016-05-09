#include <sys/time.h>

double gettimec_()
{
   struct timeval tv;
   struct timezone tz;
   gettimeofday( &tv, &tz );
   return tv.tv_sec + tv.tv_usec*1e-6;
}
