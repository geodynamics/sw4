//-*-c++-*-
#ifndef TEST_ENERGY_H
#define TEST_ENERGY_H

#include <stdlib.h>

class TestEnergy
{
public:

   TestEnergy( int seed, double cpcsratio, int write_every, std::string filename ) :
      m_seed(seed), m_cpcsratio(cpcsratio), m_write_every(write_every), m_filename(filename)
   {
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
	    for( int i=0 ; i < m_energyvector.size() ; i++ )
	       fprintf(fd, "%20.12g\n", m_energyvector[i] );
	    fclose(fd);
	 }
      }
   }

int m_seed, m_write_every;
double m_cpcsratio;
std::vector<double> m_energyvector;
std::string m_filename;
   
private:
TestEnergy(const TestEnergy&);
TestEnergy& operator=(const TestEnergy&);

};

#endif
