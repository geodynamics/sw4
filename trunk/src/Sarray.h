//-*-c++-*-
#ifndef EW_SARRAY_H
#define EW_SARRAY_H

//#include <iostream>
//#include <mpi.h>
#include "Require.h"
#include <vector>

class Sarray
{
public:
//   Sarray( CartesianProcessGrid* cartcomm, int nc=1 );
   Sarray( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend );
   Sarray( int ibeg, int iend, int jbeg, int jend, int kbeg, int kend );
   Sarray( int nc, int iend, int jend, int kend );
   Sarray( int iend, int jend, int kend );
   Sarray( const Sarray& u );
   Sarray( Sarray& u, int nc=-1 );
   Sarray();
   ~Sarray() {if( m_data != 0 ) delete[] m_data;}
//   void define( CartesianProcessGrid* cartcomm, int nc );
   void define( int iend, int jend, int kend );
   void define( int nc, int iend, int jend, int kend );
   void define( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg,
		int kend );
   void define( int ibeg, int iend, int jbeg, int jend, int kbeg,
		int kend );
   inline double* c_ptr() {return m_data;}
   void reference( double* new_data ){m_data = new_data; }

   //   inline double& operator()( int c, int i, int j, int k )
   //   {return m_data[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}
   inline bool in_range( int c, int i, int j, int k )
     {return 1 <= c && c <= m_nc && m_ib <= i && i <= m_ie && m_jb <= j 
	&& j <= m_je && m_kb <= k && k <= m_ke ;}
   inline double& operator()( int c, int i, int j, int k )
   {
#ifdef BZ_DEBUG
      VERIFY2( in_range(c,i,j,k), "Error Index (c,i,j,k) = (" << c << "," << i << "," << j << "," << k
       << ") not in range 1<= c <= " << m_nc << " " << m_ib << " <= i <= " << m_ie << " " << m_jb 
	       << " <= j <= " << m_je << " " << m_kb << " <=  k <= " << m_ke );
#endif
      return m_data[c-1+m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}

   inline double& operator()( int i, int j, int k )
      {
#ifdef BZ_DEBUG
         if (!in_range(1,i,j,k))
            VERIFY2(0, 
                    "Error Index (c,i,j,k) = (" << 1 << "," << i << "," << j << "," << k
                    << ") not in range 1<= c <= "
                    << m_nc << " "
                    << m_ib << " <= i <= " << m_ie << " "
                    << m_jb << " <= j <= " << m_je << " "
                    << m_kb << " <= k <= " << m_ke );
            
#endif
      return m_data[m_nc*(i-m_ib)+m_nc*m_ni*(j-m_jb)+m_nc*m_ni*m_nj*(k-m_kb)];}
   inline bool is_defined()
      {return m_data != NULL;}
   int m_ib, m_ie, m_jb, m_je, m_kb, m_ke;
   int index( int i, int j, int k ) {return (i-m_ib)+m_ni*(j-m_jb)+m_ni*m_nj*(k-m_kb);}
   void intersection( int ib, int ie, int jb, int je, int kb, int ke, int wind[6] );
   void side_plane( int side, int wind[6], int nGhost=1 );
   void side_plane_fortran( int side, int wind[6], int nGhost=1 );
   bool in_domain( int i, int j, int k );
   void set_to_zero();
   void set_to_minusOne();
   void set_to_random( double llim =0.0, double ulim = 1.0 );
   int ncomp() const {return m_nc;}
   int npts() const  {return m_ni*m_nj*m_nk;}
   double maximum( int c=1 );
   double minimum( int c=1 );
//   void write( char* filename, CartesianProcessGrid* cartcomm, std::vector<double> pars );
   int m_nc, m_ni, m_nj, m_nk;
private:

   double* m_data;
   inline int min(int i1,int i2){if( i1<i2 ) return i1;else return i2;}
   inline int max(int i1,int i2){if( i1>i2 ) return i1;else return i2;}
//   void init_mpi_datatype( CartesianProcessGrid* cartcomm );
//    bool m_mpi_datatype_initialized;
//    MPI_Datatype m_local_block_type;
};

#endif
