//-*-c++-*-
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
#ifndef TEST_ENERGY_H
#define TEST_ENERGY_H

#include <stdlib.h>

class TestEnergy
{
public:

   TestEnergy( int seed, double cpcsratio, int write_every, std::string filename, double amp, double sg_eps ) :
      m_seed(seed), m_write_every(write_every), m_cpcsratio(cpcsratio), m_filename(filename)
   {
      m_stochastic_amp = amp;
      m_sg_epsL = sg_eps;
      
      srand48( m_seed );
   }

   void record_data( double energy, int step, bool write_file, int myrank, string path )
   {
      m_energyvector.push_back(energy);
      if( myrank == 0 )
      {
	 if( (m_write_every>0 && step % m_write_every == 0) || write_file )
	 {
 
            stringstream filewpath;
            if( path != "." )
	       filewpath << path;
	    filewpath << m_filename;
	    //	    FILE *fd=fopen(m_filename.c_str(),"w");
	    FILE *fd=fopen(filewpath.str().c_str(),"w");
	    //	 	 cout << "energy = " << energy << endl;
	    for( unsigned int i=0 ; i < m_energyvector.size() ; i++ )
	       fprintf(fd, "%.20e\n", m_energyvector[i] );
	    fclose(fd);
	 }
      }
   }

   int m_seed, m_write_every;
   double m_cpcsratio, m_stochastic_amp, m_sg_epsL;
   std::vector<double> m_energyvector;
   std::string m_filename;
   
private:
   TestEnergy(const TestEnergy&);
   TestEnergy& operator=(const TestEnergy&);

};

#endif
