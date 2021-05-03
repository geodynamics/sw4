#include <mpi.h>
#include <vector>
#include <iostream>
#ifdef ENABLE_FFTW
#include <fftw3-mpi.h>
#endif

#include "AllDims.h"

//-----------------------------------------------------------------------
AllDims::AllDims( int nproci, int nprocj, int nprock, int ibg, int ieg, 
		  int jbg, int jeg, int kbg, int keg, int nghost, int npad,
                  MPI_Comm ewcomm )
{
//-----------------------------------------------------------------------
// General 3D array distribution.
//
// Input: nproci, nprocj, nprock - Number of processors along each dimension (i,j,k).
//        ibg, ieg - First and last i-index over total domain, excluding ghost points.
//        jbg, jeg - First and last j-index over total domain, excluding ghost points.
//        kbg, keg - First and last k-index over total domain, excluding ghost points.
//        nghost   - Number of ghost points at boundaries
//        npad     - Number of pad points between processors.
//
//  Computes the index cube in each processor, local index stored in
//        m_ib,m_ie,m_jb,m_je,m_kb,m_ke, 
//  these include padding points and ghost points.
//
//-----------------------------------------------------------------------   
   MPI_Comm_dup( ewcomm, &m_communicator );
   m_nproci=nproci;
   m_nprocj=nprocj;
   m_nprock=nprock;
   MPI_Comm_rank( m_communicator, &m_myid1d );
   compute_myid3d();

   m_ibg = ibg-nghost;
   m_ieg = ieg+nghost;
   m_jbg = jbg-nghost;
   m_jeg = jeg+nghost;
   m_kbg = kbg-nghost;
   m_keg = keg-nghost;
   m_ib.resize(m_nproci);
   m_ie.resize(m_nproci);
   m_jb.resize(m_nprocj);
   m_je.resize(m_nprocj);
   m_kb.resize(m_nprock);
   m_ke.resize(m_nprock);

   int Ntot = ieg-ibg+1+2*nghost;
   for( int p1=0 ; p1 < m_nproci ; p1++ )
   {
      decomp1d( Ntot, p1, m_nproci, m_ib[p1], m_ie[p1], nghost, npad );
      // Shift if array not 1-based.
      m_ib[p1] += ibg-1;
      m_ie[p1] += ibg-1;
   }
   Ntot = jeg-jbg+1+2*nghost;
   for( int p2=0 ; p2 < m_nprocj ; p2++ )
   {
      decomp1d( Ntot, p2, m_nprocj, m_jb[p2], m_je[p2], nghost, npad );
      m_jb[p2] += jbg-1;
      m_je[p2] += jbg-1;
   }
   Ntot = keg-kbg+1+2*nghost;
   for( int p3=0 ; p3 < m_nprock ; p3++ )
   {
      decomp1d( Ntot, p3, m_nprock, m_kb[p3], m_ke[p3], nghost, npad );
      m_kb[p3] += kbg-1;
      m_ke[p3] += kbg-1;
   }
   m_npad = npad;
   m_nghost = nghost;
   m_indrev = false;
}

//-----------------------------------------------------------------------
AllDims::AllDims( int nprocs, int ibg, int ieg, int jbg, int jeg,
		  int kbg, int keg, int nghost, MPI_Comm ewcomm )
{
// Use FFTW array distribution (i-direction split onto the processors without overlap).

   MPI_Comm_dup( ewcomm, &m_communicator );

   MPI_Comm_rank( m_communicator, &m_myid1d );
   m_nproci = nprocs;
   m_nprocj = 1;
   m_nprock = 1;

   m_ibg = ibg-nghost;
   m_ieg = ieg+nghost;
   m_jbg = jbg-nghost;
   m_jeg = jeg+nghost;
   m_kbg = kbg-nghost;
   m_keg = keg+nghost;

   // FFTW array distribution
   int nig=ieg-ibg+1+2*nghost;
   int njg=jeg-jbg+1+2*nghost;
   int nkg=keg-kbg+1+2*nghost;

#ifdef ENABLE_FFTW
   ptrdiff_t ni, ib;
   ptrdiff_t fftw_alloc_local = fftw_mpi_local_size_3d( nig, njg, nkg, m_communicator, &ni, &ib );
   m_fftw_alloc_local = static_cast<size_t>(fftw_alloc_local);
#else
   int ni, ib;
#endif

   std::vector<int> niloc(m_nproci), ibloc(m_nproci);
   niloc[m_myid1d] = ni;
   ibloc[m_myid1d] = ib;

   MPI_Allgather( &ni, 1, MPI_INT, &niloc[0], 1, MPI_INT, m_communicator );
   MPI_Allgather( &ib, 1, MPI_INT, &ibloc[0], 1, MPI_INT, m_communicator );

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
}

//-----------------------------------------------------------------------
AllDims::AllDims( AllDims* fine, int ibg, int ieg, int jbg, int jeg, 
		  int kbg, int keg, int nghost, int npad )
{
   // Construct coarser grid from already distributed fine grid
   // Number of processors and my id's are same as fine grid.
   m_nproci=fine->m_nproci;
   m_nprocj=fine->m_nprocj;
   m_nprock=fine->m_nprock;
   m_myid1d=fine->m_myid1d;
   m_myid3di = fine->m_myid3di;
   m_myid3dj = fine->m_myid3dj;
   m_myid3dk = fine->m_myid3dk;
   MPI_Comm_dup( fine->m_communicator, &m_communicator );

   m_ibg = ibg-nghost;
   m_ieg = ieg+nghost;
   m_jbg = jbg-nghost;
   m_jeg = jeg+nghost;
   m_kbg = kbg-nghost;
   m_keg = keg+nghost;

   m_ib.resize(m_nproci);
   m_ie.resize(m_nproci);
   m_jb.resize(m_nprocj);
   m_je.resize(m_nprocj);
   m_kb.resize(m_nprock);
   m_ke.resize(m_nprock);

   int nghostf = fine->m_nghost;
   int npadf   = fine->m_npad;
   int Nf = (fine->m_ieg-fine->m_ibg+1)-2*nghostf;
   for( int p1=0 ; p1 < m_nproci ; p1++ )
   {
      decomp1d_frf( ieg-ibg+1, p1, m_nproci, m_ib[p1], m_ie[p1], nghost, npad, 
		    Nf, fine->m_ib[p1], fine->m_ie[p1], nghostf, npadf );
      m_ib[p1] += ibg-1;
      m_ie[p1] += ibg-1;
   }
   Nf = (fine->m_jeg-fine->m_jbg+1)-2*nghostf;
   for( int p2=0 ; p2 < m_nprocj ; p2++ )
   {
      decomp1d_frf( jeg-jbg+1, p2, m_nprocj, m_jb[p2], m_je[p2], nghost, npad,
		    Nf, fine->m_jb[p2], fine->m_je[p2], nghostf, npadf );
      m_jb[p2] += jbg-1;
      m_je[p2] += jbg-1;
   }
   Nf = (fine->m_keg-fine->m_kbg+1)-2*nghostf;
   for( int p3=0 ; p3 < m_nprock ; p3++ )
   {
      decomp1d_frf( keg-kbg+1, p3, m_nprock, m_kb[p3], m_ke[p3], nghost, npad,
		    Nf, fine->m_kb[p3], fine->m_ke[p3], nghostf, npadf );
      m_kb[p3] += kbg-1;
      m_ke[p3] += kbg-1;
   }
   m_npad   = npad;
   m_nghost = nghost;
   m_indrev = false;
}

//-----------------------------------------------------------------------
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
      
//-----------------------------------------------------------------------
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
   
//-----------------------------------------------------------------------
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

//-----------------------------------------------------------------------
int AllDims::proc1d( int p1, int p2, int p3 )
{
   //   return p1 + m_nproci*(p2+m_nprocj*p3);
   return p3 + m_nprock*(p2+m_nprocj*p1);
}

//-----------------------------------------------------------------------
bool AllDims::owner( int p1, int p2, int p3 )
{
   return p1==m_myid3di && p2==m_myid3dj && p3==m_myid3dk;
}

//-----------------------------------------------------------------------
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

//-----------------------------------------------------------------------
void AllDims::decomp1d( int nglobal, int myid, int nproc, int& s, int& e, int nghost, int npad )
//
// Decompose index space 1-nghost <= i <= N+nghost into nproc blocks
//
// Input: nglobal - Total number of points = N+2*nghost.
//        myid    - Processor ID of current task,  0 <= myid <= nproc-1.
//        nproc   - Total number of processors (tasks).
//        nghost  - Number of ghost points at domain boundaries.
//        npad    - Number of overlap (padding) points at processor boundaries.
//   
// Output: s - Low index in this processor.
//         e - High index in this processor, ie, current task holds s <= i <= e
//
// The nglobal points are distributed as evenly as possible on the tasks.
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

//-----------------------------------------------------------------------
void AllDims::decomp1d_2( int N, int myid, int nproc, int& s, int& e, int nghost, int npad )
//
// Decompose index space 1-nghost <= i <= N+nghost into nproc blocks
//
// Input: N      - Number of points in domain.
//        myid   - Processor ID of current task,  0 <= myid <= nproc-1.
//        nproc  - Total number of processors (tasks).
//        nghost - Number of ghost points at domain boundaries.
//        npad   - Number of overlap (padding) points at processor boundaries.
//   
// Output: s - Low index in this processor.
//         e - High index in this processor, ie, current task holds s <= i <= e
//
// The N points are distributed as evenly as possible on the tasks. Ghost points and
// padding points are added after distribution.
//
{
   int nglobal = N+2*nghost;
   int olap    = 2*npad;

   int nlocal  = N / nproc;
   int deficit = N % nproc;

   if( myid < deficit )
      s = myid*nlocal + myid+1;
   else
      s = myid*nlocal + deficit+1;

   if (myid < deficit)
      nlocal = nlocal + 1;

   e = s + nlocal - 1;
   if( myid == nproc-1 )
      e += nghost;
   else
      e += npad;
   if( myid == 0 )
      s -= nghost;
   else
      s -= npad;
}

//-----------------------------------------------------------------------
void AllDims::decomp1d_frf( int N, int myid, int nproc, int& s, int& e, int nghost, 
			    int npad, int Nf, int sf, int ef, int nghostf, int npadf )
{
// Decompose index space 1-nghost <= i <= N+nghost into nproc blocks
// The decomposition is done subordinate to an already distributed finer grid
//
// Input: N      - Number of points in domain (without ghost points).
//        myid   - Processor ID of current task,  0 <= myid <= nproc-1.
//        nproc  - Total number of processors (tasks).
//        nghost - Number of ghost points at domain boundaries.
//        npad   - Number of overlap (padding) points at processor boundaries.
//        Nf     - Number of points in domain (without ghost points) of finer grid.
//        sf     - First index of fine grid in this processor (including ghost points and pad points)
//        ef     - Last index of fine grid in this processor (including ghost points and pad points)
//        nghostf- Number of ghost points used by fine grid.
//        npadf  - Number of pad points used by fine grid.
//
// Output: s - Low index in this processor.
//         e - High index in this processor, ie, current task holds s <= i <= e
//

   double hrat = static_cast<double>(N-1)/(Nf-1);   

   if(  myid == 0 )
      sf += nghostf;
   else
      sf += npadf;
   if( myid ==  nproc-1 )
      ef -= nghostf;
   else
      ef -= npadf;

   if( myid == 0 )
      s = 1;
   else
   {
      // First point in this processor (not counting pad and ghost points)
      // is first point to the right of interface pt xf_{sf-1/2} on fine grid.
      s = static_cast<int>( (sf-1-0.5)*hrat )+2;
   }
   
   if( myid ==  nproc-1 )
      e = N;
   else
   {
      // Last point in this processor (not counting pad and ghost points)
      // is last point to the left of interface pt xf_{ef+1/2} on fine grid.
      e = static_cast<int>((ef-1+0.5)*hrat)+1;
   }
   if( myid == 0 )
      s -= nghost;
   else
      s -= npad;
   if( myid == nproc-1 )
      e += nghost;
   else
      e += npad;
}
