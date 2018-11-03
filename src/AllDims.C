#include <mpi.h>
#include <vector>
#include <iostream>
#ifdef ENABLE_FFTW
#include <fftw3-mpi.h>
#endif

#include "AllDims.h"

AllDims::AllDims( int nproci, int nprocj, int nprock, int ibg, int ieg, 
		  int jbg, int jeg, int kbg, int keg, int nghost, int npad )
{
   // SW4 array distribution
   m_nproci=nproci;
   m_nprocj=nprocj;
   m_nprock=nprock;
   MPI_Comm_rank( MPI_COMM_WORLD, &m_myid1d );
   compute_myid3d();
   m_ibg = ibg;
   m_ieg = ieg;
   m_jbg = jbg;
   m_jeg = jeg;
   m_kbg = kbg;
   m_keg = keg;
   m_ib.resize(m_nproci);
   m_ie.resize(m_nproci);
   m_jb.resize(m_nprocj);
   m_je.resize(m_nprocj);
   m_kb.resize(m_nprock);
   m_ke.resize(m_nprock);
   for( int p1=0 ; p1 < m_nproci ; p1++ )
   {
      decomp1d( ieg-ibg+1, p1, m_nproci, m_ib[p1], m_ie[p1], nghost, npad );
      m_ib[p1] += ibg-1;
      m_ie[p1] += ibg-1;
   }
   for( int p2=0 ; p2 < m_nprocj ; p2++ )
   {
      decomp1d( jeg-jbg+1, p2, m_nprocj, m_jb[p2], m_je[p2], nghost, npad );
      m_jb[p2] += jbg-1;
      m_je[p2] += jbg-1;
   }
   for( int p3=0 ; p3 < m_nprock ; p3++ )
   {
      decomp1d( keg-kbg+1, p3, m_nprock, m_kb[p3], m_ke[p3], nghost, npad );
      m_kb[p3] += kbg-1;
      m_ke[p3] += kbg-1;
   }
   m_npad = npad;
   m_nghost = nghost;
   m_indrev = false;
}


AllDims::AllDims( int nprocs, int ibg, int ieg, int jbg, int jeg,
		  int kbg, int keg, int nghost )
{
   MPI_Comm_rank( MPI_COMM_WORLD, &m_myid1d );
   m_nproci = nprocs;
   m_nprocj = 1;
   m_nprock = 1;
   m_ibg = ibg;
   m_ieg = ieg;
   m_jbg = jbg;
   m_jeg = jeg;
   m_kbg = kbg;
   m_keg = keg;

   // FFTW array distribution
   int nig=ieg-ibg+1;
   int njg=jeg-jbg+1;
   int nkg=keg-kbg+1;
   ptrdiff_t ni, ib;
#ifdef ENABLE_FFTW
   ptrdiff_t fftw_alloc_local = fftw_mpi_local_size_3d( nig, njg, nkg, MPI_COMM_WORLD, &ni, &ib );
   m_fftw_alloc_local = static_cast<size_t>(fftw_alloc_local);
#endif
   
   std::vector<int> niloc(m_nproci), ibloc(m_nproci);
   niloc[m_myid1d] = ni;
   ibloc[m_myid1d] = ib;

   MPI_Allgather( &ni, 1, MPI_INT, &niloc[0], 1, MPI_INT, MPI_COMM_WORLD );
   MPI_Allgather( &ib, 1, MPI_INT, &ibloc[0], 1, MPI_INT, MPI_COMM_WORLD );

   m_ib.resize(m_nproci);
   m_ie.resize(m_nproci);
   for( int p1 = 0 ; p1 < m_nproci ; p1++ )
   {
      m_ib[p1] = ibloc[p1]+ibg;
      m_ie[p1] = niloc[p1]+m_ib[p1]-1;
   }
   m_myid3di = m_myid1d;
   m_myid3dj = 0;
   m_myid3dk = 0;
   m_jb.resize(1);
   m_je.resize(1);
   m_kb.resize(1);
   m_ke.resize(1);
   m_jb[0] = jbg;
   m_je[0] = jeg;
   m_kb[0] = kbg;
   m_ke[0] = keg;
   m_npad = 0;
   m_nghost=nghost;
   m_indrev = true;
   //      if( m_myid1d == 7 )
   //      {
   //	 cout << m_myid1d << "dims " << m_ib[m_myid3di] << " " << m_ie[m_myid3di] << " " 
   //	      << m_jb[m_myid3dj] << " " << m_je[m_myid3dj] << " " 
   //	      << m_kb[m_myid3dk] << " " << m_ke[m_myid3dk] << endl; 
   //      }

}

void AllDims::getdims_nopad( int dims[6], int p1, int p2, int p3 )
{
   if( p1 == -1 )
   {
      p1 = m_myid3di;
      p2 = m_myid3dj;
      p3 = m_myid3dk;
   }
   dims[0] = m_ib[p1];
   if( p1 != 0 )
      dims[0] += m_npad;
   dims[1] = m_ie[p1];
   if( p1 != m_nproci-1 )
      dims[1] -= m_npad;
   dims[2] = m_jb[p2];
   if( p2 != 0 )
      dims[2] += m_npad;
   dims[3] = m_je[p2];
   if( p2 != m_nprocj-1 )
      dims[3] -= m_npad;
   dims[4] = m_kb[p3];
   if( p3 != 0 )
      dims[4] += m_npad;
   dims[5] = m_ke[p3];
   if( p3 != m_nprock-1 )
      dims[5] -= m_npad;
}
      
bool AllDims::intersect( int p1, int p2, int p3, AllDims& other, int dims[6] )
{
   // Intersection of patch(procno) with other.patch(myid)
   // First remove padding points
   int dims1[6];
   getdims_nopad( dims1, p1, p2, p3 );
   int dims2[6];
   other.getdims_nopad( dims2 );
      
   if( dims1[0] <= dims2[1] && dims1[1] >= dims2[0] &&
       dims1[2] <= dims2[3] && dims1[3] >= dims2[2] &&
       dims1[4] <= dims2[5] && dims1[5] >= dims2[4] )
   {
      dims[0] = dims1[0] >= dims2[0] ? dims1[0] : dims2[0];
      dims[1] = dims1[1] <= dims2[1] ? dims1[1] : dims2[1];
      dims[2] = dims1[2] >= dims2[2] ? dims1[2] : dims2[2];
      dims[3] = dims1[3] <= dims2[3] ? dims1[3] : dims2[3];
      dims[4] = dims1[4] >= dims2[4] ? dims1[4] : dims2[4];
      dims[5] = dims1[5] <= dims2[5] ? dims1[5] : dims2[5];
      return true;
   }
   else
      return false;
}
   
void AllDims::compute_myid3d( )
{
   // Note: SW4 uses built-in mpi cartesian communicator,
   // which orders the processes as pk+npk*pj+npj*npk*pi
   //  
   //   m_myid3di = m_myid1d % m_nproci;
   //   int rem=(m_myid1d-m_myid3di)/m_nproci;
   //   m_myid3dj = rem % m_nprocj;
   //   m_myid3dk = (rem-m_myid3dj)/m_nprocj;
   m_myid3dk = m_myid1d % m_nprock;
   int rem=(m_myid1d-m_myid3dk)/m_nprock;
   m_myid3dj = rem % m_nprocj;
   m_myid3di = (rem-m_myid3dj)/m_nprocj;   
}
int AllDims::proc1d( int p1, int p2, int p3 )
{
   //   return p1 + m_nproci*(p2+m_nprocj*p3);
   return p3 + m_nprock*(p2+m_nprocj*p1);
}
bool AllDims::owner( int p1, int p2, int p3 )
{
   return p1==m_myid3di && p2==m_myid3dj && p3==m_myid3dk;
}

int AllDims::owner_i( int i )
{
   // which procssor along the i-dimension owns index i ?
   double rat=(i-m_ibg)/(m_ieg-m_ibg);
   int proc = rat*m_nproci;
   if( proc >= m_nproci )
      proc = m_nproci-1;
   //   std::cout << "init proc " << proc << " mm_nproci = " << m_nproci << std::endl;
   if( m_ib[proc] <= i && i <= m_ie[proc] )
      return proc;
   else
   {
      proc = 0;
      while( proc < m_nproci-1 && !(m_ib[proc]<= i && i <= m_ie[proc] ) )
	 proc++;
      if( proc > m_nproci-1 )
      {
	 std::cout << "ERROR in AllDims::owner_i, proc= " << proc << std::endl;
	 return -1;
      }
      else if( m_ib[proc] <= i && i <= m_ie[proc] )
	 return proc;
      else
      {
	 std::cout << "inexplacable ERROR in AllDims::owner_i " << proc <<
	    " " << m_ib[proc] << " " << m_ie[proc] << std::endl;
	 return -1;
      }
   }
}

//   }
//  else if( m_ie[proc] < i )
//   {
//      // search to the right
//      while( proc < m_nproci-1 && m_ie[proc] < i )
//	 proc++;
//      //   std::cout << "out+ proc " << proc << std::endl;
//      if( m_ib[proc] <= i && i <= m_ie[proc] )
//	 return proc;
//      else
//      {
//	 std::cout << "ERROR in AllDims::owner_i+" << proc << " " << m_ib[proc] << " " <<  m_ie[proc]  << std::endl;
//	 return -1;
//      }
//   }
//   else if( m_ib[proc] > i )
//   {
//      // search to the left	 
//      while( proc > 0 && m_ib[proc]>i )
//	 proc--;
      //   std::cout << "out- proc " << proc << std::endl;
//      if( m_ib[proc] <= i && i <= m_ie[proc] )
//	 return proc;
//      else
//      {
//	 std::cout << "ERROR in AllDims::owner_i- " << proc << " " << m_ib[proc] << " " <<  m_ie[proc]  << std::endl;
//	 return -1;
//      }
//   }
//   else
//   {
//      std::cout << "ERROR in AllDims::owner_i0 " << proc << " " << m_ib[proc] << " " <<  m_ie[proc]  << std::endl;
//      return -1;
//   }
//}

void AllDims::decomp1d( int nglobal, int myid, int nproc, int& s, int& e, int nghost, int npad )
//
// Decompose index space 1 <= i <= nglobal into nproc blocks
// returns start and end indices for block nr. myid, 
//          where 0 <= myid <= nproc-1
// It is assumed that nglobal includes nghost ghost points at each end.
// npad is number of padding points at processor interfaces.
//
{
   int olap    = 2*npad;
   int nlocal  = (nglobal + (nproc-1)*olap ) / nproc;
   int deficit = (nglobal + (nproc-1)*olap ) % nproc;

   if( myid < deficit )
      s = myid*(nlocal-olap) + myid+1;
   else
      s = myid*(nlocal-olap) + deficit+1;

   if (myid < deficit)
      nlocal = nlocal + 1;

   e = s + nlocal - 1;
   e -= nghost;
   s -= nghost;
}

