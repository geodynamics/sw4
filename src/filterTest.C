#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Filter.h"
#include "SecondOrderSection.h"

using namespace std;

int
main(int argc, char **argv)
{
  unsigned int numberOfPasses = 1;
  double f1=0.5;
  double f2=2.0;
  double dt=1.0/100.0;
  
  Filter filter(bandPass, 3, numberOfPasses, f1, f2);
//  Filter filter(lowPass, 5, numberOfPasses, f1, f2);

 filter.computeSOS(dt);
 cout << filter;

  double u[1000], fu[1000];
  int N=1000, i;
// make a Heaviside function
  for (i=0; i<N; i++)
    u[i] = 0;
  for (i=N/2; i<N; i++)
    u[i] = 1;
  
  filter.evaluate(N, u, fu);

  FILE *fd=fopen("h1.dat","w");
  for (i=0; i<N; i++)
  {
    fprintf(fd,"%e %e %e\n", dt*i, u[i], fu[i]);
  }
  fclose(fd);
  
  
} // end of main
