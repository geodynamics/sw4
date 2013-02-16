#ifndef DATAPATCHES_H
#define DATAPATCHES_H

#include "Sarray.h"
#include <vector>
#include <string>

class DataPatches
{
   // Number of subcubes in domain to save
   int  m_npatches;

   // dimensions of cube p, stored as:
   //  imin  = m_dims[6*p],   imax = m_dims[6*p+1]
   //  jmin  = m_dims[6*p+2], jmax = m_dims[6*p+3]
   //  kmin  = m_dims[6*p+4], kmax = m_dims[6*p+5]
   std::vector<size_t> m_dims;

   // pointer to data array. Start index of subcube p is m_dataptr[p]
   size_t* m_dataptr;

   // Number of components per grid point to store, typically 3 for the 3D elastic wave eq
   int m_ncomp;

   // First and last step stored on the file.
   int m_nmin, m_nmax;

   // Number of time steps to hold in memory
   int m_nsteps;

   // Number of time steps currently in memory
   int m_ncurrent;

   size_t m_ni, m_nj;

   // Name of storage file
   std::string m_filename;

   // Stored data m_data[n] is the n:th time step in memory
   std::vector<double*> m_data;

   // Translate step number in memory to step number in solver
   // i.e., m_steps[n] is the step number of the n:th time step in memory.
   int* m_steps;

   // Check if it the first time we save to file, need to write header then.
   bool m_startedsave;

   // Error flag.
   bool m_error;
   bool m_isnonempty;

   void save_to_file( );
   void read_from_file( int n );
   void add_patch( int wind[6] );
public:
   DataPatches( std::string fname, Sarray& u, int imin, int imax, int jmin, int jmax, int kmax,
		int nlayers, int ntsteps, double dt );
   ~DataPatches();
   void push( Sarray& u, int n );
   void pop( Sarray& u, int n );
};

#endif
