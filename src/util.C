#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>

#include "util.h"

void checkMinMax(int n, double* a, char* name)
{
    float min=1e20;
    float max=-1e20;

   for(int i=0; i<n; i++) {
       if(a[i]<min) min=a[i];
       if(a[i]>max) max=a[i];
   }
   cout << ">>>>>>>>>> " << name << ":" << " min=" << min << " max=" << max << endl;
}

void save_array_to_disk(int n, double* a, char* fname)
{
  size_t nr;

  int fd = open(fname, O_CREAT | O_TRUNC | O_WRONLY, 0660 );
   if( fd == -1 )
    cerr << "ERROR opening file" << fname << " for writing " << endl;
   
      nr = write(fd,a,sizeof(double)*n);
      if( nr != sizeof(double)*n )
      cout << "Error saving data to " << fname << endl;
   close(fd);
}