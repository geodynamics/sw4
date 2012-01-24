#include <sys/time.h>

int zeropt;
int first=1;

double gettimec()
{
   struct timeval tv;
   struct timezone tz;
   gettimeofday( &tv, &tz );
   if( first )
   {
     zeropt = tv.tv_sec;
     first = 0;
   }
   return tv.tv_sec-zeropt + tv.tv_usec*1e-6;
}
