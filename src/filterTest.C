//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
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
