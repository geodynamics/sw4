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
#include <mpi.h>
#include <iostream>
using namespace std;

#include <cstring>
#include <errno.h>
#include <unistd.h>
#include "Parallel_IO.h"
#include <fcntl.h>
//-----------------------------------------------------------------------
Comminfo::Comminfo()
{
   m_has_values = false;
   m_steps = 0;
   // Initialize pointers to nil
   m_ncomm = NULL;
   m_comm_id = NULL;
   for( int p=0 ; p < 6 ; p++ )
      m_comm_index[p] = NULL;
   m_ilow = NULL;
   m_jlow = NULL;
   m_klow = NULL;
   m_niblock = NULL;
   m_njblock = NULL;
   m_nkblock = NULL;

   m_nsubcomm = NULL;
   m_subcomm = NULL;
   m_subcommlabel = NULL;

   m_maxbuf = 0;
   m_maxiobuf=0;
}

//-----------------------------------------------------------------------
Comminfo::~Comminfo()
{
   if( m_ncomm != NULL )
      delete[] m_ncomm;
   if( m_ilow != NULL )
      delete[] m_ilow;
   if( m_jlow != NULL )
      delete[] m_jlow;
   if( m_klow != NULL )
      delete[] m_klow;
   if( m_niblock != NULL )
      delete[] m_niblock;
   if( m_njblock != NULL )
      delete[] m_njblock;
   if( m_nkblock != NULL )
      delete[] m_nkblock;
   if( m_comm_id != NULL )
   {
      for( int p=0 ; p < m_steps ; p++ )
	 if( m_comm_id[p] != NULL )
	    delete[] m_comm_id[p];
      delete[] m_comm_id;
   }   
   for( int p= 0 ; p < 6 ; p++ )
   {
      if( m_comm_index[p] != NULL )
      {
	 for( int s=0 ; s < m_steps ; s++ )
	    if( m_comm_index[p][s] != NULL )
	       delete[] m_comm_index[p][s];
	 delete[] m_comm_index[p];
      }
   }

   if( m_nsubcomm != NULL )
   {
      for( int p=0 ; p < m_steps ; p++ )
	 if( m_subcomm[p] != NULL )
	    delete[] m_subcomm[p];
      delete[] m_subcomm;
      delete[] m_nsubcomm;
   }
   if( m_subcommlabel != NULL )
   {
      for( int p=0 ; p < m_steps ; p++ )
	 if( m_subcommlabel[p] != NULL )
	    delete[] m_subcommlabel[p];
      delete[] m_subcommlabel;
   }
}

//-----------------------------------------------------------------------
void Comminfo::print( int recv )
{
   if( m_has_values )
   {
      bool printsends=true;
      cout << "maxbuf= " << m_maxbuf << endl;
      cout << "steps= " << m_steps << endl;
      if( printsends )
      {
	 for( int s=0 ; s < m_steps ; s++ )
	 {
	    cout << "step no " << s << endl;
	    cout << "      ncomm= " << m_ncomm[s] << endl;
	    for( int n=0 ; n < m_ncomm[s] ; n++ )
	    {
	       cout << "      ncomm_id= " << m_comm_id[s][n] << endl;
	       cout << "      comm inds ";
	       for( int side=0 ; side < 6 ; side++ )
		  cout << m_comm_index[side][s][n] << " ";
	       cout << endl;
	       if( !recv )
		  cout << " subcommlabel = " << m_subcommlabel[s][n] << endl;
	    }
	 }
      }
      if( recv == 1 )
      {
	 for( int s=0 ; s < m_steps ; s++ )
	 {
            cout << " step " << s << endl;
	    cout << "(i,j,k) wr block " << m_ilow[s] << " " << m_jlow[s] << " " << m_klow[s] <<endl;
	    cout << "size wr block    " << m_niblock[s] << " " << m_njblock[s] << " " << m_nkblock[s] <<endl;
	    if( m_ncomm[s] > 0 )
	    {
	       cout << "number of substeps = " << m_nsubcomm[s] << " intervals = " << endl;
	       for( int ss= 0 ; ss <= m_nsubcomm[s] ; ss++ )
		  cout << " " << m_subcomm[s][ss] << " ";
	       cout << endl;
	    }  
	 }
      }
   }
}

//-----------------------------------------------------------------------
Parallel_IO::Parallel_IO( int iwrite, int pfs, int globalsizes[3], int localsizes[3],
		  int starts[3], int nptsbuf, int padding )
{
   int ihave_array=1;
   if( localsizes[0] < 1 || localsizes[1] < 1 || localsizes[2] < 1 )
      ihave_array = 0;
   init_pio( iwrite, pfs, ihave_array );
   int myid;
   MPI_Comm_rank( MPI_COMM_WORLD, &myid );
   if( myid == 1 )
   {
      cout << myid << " Parallel_IO:: Done init_pio " << endl;
   cout << "gsizes " << globalsizes[0] <<  " " << globalsizes[1] << " " << globalsizes[2] << endl;
   cout << "lsizes " << localsizes[0] <<  " " << localsizes[1] << " " << localsizes[2] << endl;
   cout << "ssizes " << starts[0] <<  " " << starts[1] << " " << starts[2] << endl;
   }
   init_array( globalsizes, localsizes, starts, nptsbuf, padding );
   MPI_Barrier(MPI_COMM_WORLD);
   if( myid == 1 )
      cout << "Parallel_IO:: Done init_array " << endl;

   if( m_data_comm != MPI_COMM_NULL )
   {
      //      int mywid;
      //      MPI_Comm_rank( m_data_comm, &mywid );
      //      if( m_writer_ids[0]==mywid )
      //      int myid ;
      //      MPI_Comm_rank( MPI_COMM_WORLD, &myid );
      //      if( myid == 0 )
      //      {
      //         cout << "ISEND = " << endl;
      //         m_isend.print(0);
      //	 if( m_iwrite )
      //	 {
      //	    cout << "IRECV = " << endl;
      //	    m_irecv.print(1);
      //	 }
      //      }
   }
   //            if( myid == 0 )
   //            {
   //               cout << "IRECV = " << endl;
   //               m_irecv.print(1);
   //            }
   //            if( myid == 0 )
   //            {
   //               cout << "ISEND = " << endl;
   //               m_isend.print(0);
   //            }
}

//-----------------------------------------------------------------------
void Parallel_IO::init_pio( int iwrite, int pfs, int ihave_array )
//-----------------------------------------------------------------------
// Initialize for parallel I/O, 
// Input: iwrite - 0 this processor will not participate in I/O op.
//                 1 this processor will participate in I/O op.
//        pfs    - 0 I/O will be done to non-parallel file system.
//                 1 I/O will be done to parallel file system.
//        ihave_array - 0 this processor holds no part of the array.
//                      1 this processor holds some part of the array.
//                     -1 (default) the array is present in all processors.
//
// Constructs communicators from input iwrite and ihave_array:
//
//    m_data_comm  - All processors that holds part of the array or 
//                   writes to disc, or both, i.e., iwrite .or. ihave_array
//    m_write_comm - Subset of m_data_comm that writes to disk, i.e., iwrite.
//
// m_writer_ids[i] is id in m_data_comm of proc i in m_write_comm.
//-----------------------------------------------------------------------
{
   int* tmp;
   int nprocs, p, i, retcode;
   MPI_Group world_group, writer_group, array_group;

   // Global processor no for error messages
   int gproc;
   MPI_Comm_rank( MPI_COMM_WORLD, &gproc );

// 0. Default communicators are NULL.
   m_data_comm = m_write_comm = MPI_COMM_NULL;

// 1. Create communicator of all procs that either hold part of the array or will perform I/O.
//    Save as m_data_comm. This communicator will be used for most operations.
   if( ihave_array == -1 )
      m_data_comm = MPI_COMM_WORLD;
   else
   {
      MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
      try
      {
	 tmp      = new int[nprocs];
      }
      catch( bad_alloc& ba )
      {
	 cout << "Parallel_IO::init_pio, processor " << gproc << ". Allocation of tmp failed. "
	      << " nprocs = " << nprocs << endl;
	 MPI_Abort(MPI_COMM_WORLD,0);
      }

      int participator = ihave_array || iwrite;
      retcode = MPI_Allgather( &participator, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD );
      if( retcode != MPI_SUCCESS )
      {
	 cout << "Parallel_IO::init_pio, error from first call to MPI_Allgather, "
	      << "return code = " << retcode << " from processor " << gproc << endl;
      }
      int npartprocs = 0;
      for( p = 0 ; p < nprocs ; p++ )
	 if( tmp[p] == 1 )
	    npartprocs++;
      if( npartprocs < 1 )
      {
	 // error 
         cout << "ERROR 1, in Parallel_IO::init_pio, no processors selected, npartprocs =  " << npartprocs << endl;
         delete[] tmp;
         return;
      }
      int* array_holders;
      try
      {
	 array_holders = new int[npartprocs];
      }
      catch( bad_alloc &ba )
      {
	 cout << "Parallel_IO::init_pio, processor " << gproc << ". Allocation of array_holders failed. "
	      << " npartprocs = " << npartprocs << endl;
	 MPI_Abort(MPI_COMM_WORLD,0);
      }
      i = 0;
      for( p=0 ; p < nprocs ; p++ )
      {
	 if( tmp[p] == 1 )
	    array_holders[i++] = p;
      }
      delete[] tmp;
      retcode = MPI_Comm_group( MPI_COMM_WORLD, &world_group );
      if( retcode != MPI_SUCCESS )
      {
	 cout << "Parallel_IO::init_pio, error from first call to MPI_Comm_group, "
	      << "return code = " << retcode << " from processor " << gproc << endl;
      }
      retcode = MPI_Group_incl( world_group, npartprocs, array_holders, &array_group );
      if( retcode != MPI_SUCCESS )
      {
	 cout << "Parallel_IO::init_pio, error from first call to MPI_Group_incl, "
	      << "return code = " << retcode << " from processor " << gproc << endl;
      }
      retcode = MPI_Comm_create( MPI_COMM_WORLD, array_group, &m_data_comm );
      if( retcode != MPI_SUCCESS )
      {
	 cout << "Parallel_IO::init_pio, error from first call to MPI_Comm_create, "
	      << "return code = " << retcode << " from processor " << gproc << endl;
      }
      MPI_Group_free( &world_group );
      MPI_Group_free( &array_group );
      delete[] array_holders;
   }

   if( m_data_comm != MPI_COMM_NULL )
   {
      MPI_Comm_size( m_data_comm, &nprocs );
// 2. Create communicator of all I/O processors. Save as m_write_comm. This communicator is 
//    used in MPI_Barrier to synchronize read and writes on non-parallel file systems
      m_iwrite = iwrite;
      try
      {
	 tmp      = new int[nprocs];
      }
      catch( bad_alloc &ba )
      {
	 cout << "Parallel_IO::init_pio, processor " << gproc << ". Allocation of second tmp failed. "
	      << " nprocs = " << nprocs << endl;
	 MPI_Abort(MPI_COMM_WORLD,0);
      }
      retcode = MPI_Allgather( &m_iwrite, 1, MPI_INT, tmp, 1, MPI_INT, m_data_comm );
      if( retcode != MPI_SUCCESS )
      {
	 cout << "Parallel_IO::init_pio, error from second call to MPI_Allgather, "
	      << "return code = " << retcode << " from processor " << gproc << endl;
      }
      m_nwriters = 0;
      for( p = 0 ; p < nprocs ; p++ )
	 if( tmp[p] == 1 )
	    m_nwriters++;
      if( m_nwriters < 1 )
      {
         cout << "ERROR 2 in Parallel_IO::init_pio, there are no writing processors left, nwriters = " << m_nwriters << endl;
         delete[] tmp;
         return;
      }
 // writer_ids are the processor ids of the writer processors in the enumeration by the m_data_comm group.
      try
      {
	 m_writer_ids = new int[m_nwriters];
      }
      catch( bad_alloc &ba )
      {
	 cout << "Parallel_IO::init_pio, processor " << gproc << ". Allocation of m_writer_ids failed. "
	      << " nwriters = " << m_nwriters << endl;
	 MPI_Abort(MPI_COMM_WORLD,0);
      }
      i = 0;
      for( p=0 ; p < nprocs ; p++ )
      {
	 if( tmp[p] == 1 )
	    m_writer_ids[i++] = p;
      }
      delete[] tmp;

      retcode = MPI_Comm_group( m_data_comm, &world_group );
      if( retcode != MPI_SUCCESS )
      {
	 cout << "Parallel_IO::init_pio, error from second call to MPI_Comm_group, "
	      << "return code = " << retcode << " from processor " << gproc << endl;
      }
      retcode = MPI_Group_incl( world_group, m_nwriters, m_writer_ids, &writer_group );
      if( retcode != MPI_SUCCESS )
      {
	 cout << "Parallel_IO::init_pio, error from second call to MPI_Group_incl, "
	      << "return code = " << retcode << " from processor " << gproc << endl;
      }
      retcode = MPI_Comm_create( m_data_comm, writer_group, &m_write_comm );
      if( retcode != MPI_SUCCESS )
      {
	 cout << "Parallel_IO::init_pio, error from second call to MPI_Comm_create, "
	      << "return code = " << retcode << " from processor " << gproc << endl;
      }
      MPI_Group_free( &world_group );
      MPI_Group_free( &writer_group );
   }
   // 3. Save parallel file system info
   m_parallel_file_system = pfs;

   // 4. Initialize some variables to shut up the memory checker 
   if( m_data_comm == MPI_COMM_NULL )
      m_iwrite = false;
}

//-----------------------------------------------------------------------
void Parallel_IO::init_array( int globalsizes[3], int localsizes[3], 
			      int starts[3], int nptsbuf, int padding )
{
// Set up data structures for communication before I/O
// Input: globalsizes - global size of array ( [1..nig]x[1..njg]x[1..nkg] )
//        localsizes  - Size of array subblock in this processor
//        starts      - Location of subblock in global array
//                      This processor holds the subrange
//                        1+starts[0] <= i <= localsizes[0]+starts[0]
//                      of the global range 1 <= i <= globalsizes[0]
//                       and similarly for j and k.
//        nptsbuf     - Number of grid points in temporary buffer
//        padding     - If there is overlap between processors, setting
//                      padding avoids writing these twice.
   int blsize, s, blocks_in_writer, r, p, b, blnr, kb, ke, l;
   int ibl, iel, jbl, jel, kbl, kel, nsend;
   int found, i, j, q, lims[6], v[6], vr[6], nprocs, tag2, myid;
   int retcode, gproc;
   int* nrecvs;
   size_t nblocks, npts, maxpts;

   MPI_Status status;
   //   int nrw, flim = 600;
   double t0 = MPI_Wtime();

   if( m_data_comm != MPI_COMM_NULL )
   {
      MPI_Comm_rank( MPI_COMM_WORLD, &gproc );
      
   ni  = localsizes[0];
   nj  = localsizes[1];
   nk  = localsizes[2];

   nig = globalsizes[0];
   njg = globalsizes[1];
   nkg = globalsizes[2];

   oi  = starts[0];
   oj  = starts[1];
   ok  = starts[2];

   MPI_Comm_size( m_data_comm, &nprocs );
   MPI_Comm_rank( m_data_comm, &myid );

   if( m_iwrite == 1 )
   {
      try{
	 nrecvs = new int[nprocs];
      }
      catch( bad_alloc& ba )
      {
	 cout << "Parallel_IO::init_array, Processor " << gproc <<  ". Allocation of nrecvs failed "
	      << " nprocs = " << nprocs << " Exception= " << ba.what() << endl;
	 MPI_Abort(MPI_COMM_WORLD,0);
      }
   }

   // 3. Split the domain into strips
   // Takes care of 3D and all 2D cases
   // nkg > 1 --> split k
   // nkg = 1 --> split j

   nblocks = static_cast<size_t>((1.0*nig*((size_t)njg)*nkg)/nptsbuf);
   if( (((off_t)nig)*njg*nkg % ((off_t)nptsbuf) ) != 0 )
      nblocks++;

// blsize is 1D size of each block along the last dimension
   if( nkg == 1 )
   {
      blsize = njg/nblocks;
      s = njg % nblocks;
   }
   else
   {
      blsize = nkg/nblocks;
      s = nkg % nblocks;
   }
   blocks_in_writer = nblocks/m_nwriters;
   r = nblocks % m_nwriters;

   //      cout << "myid  " << myid << " nblocks = "<<nblocks << "ni,nj,nk " << ni << " " << nj << " " << nk << endl;
   //      cout << "blsize = " << blsize << " s " << s << endl;
      /** b= number of blocks written by each writer if s=0 */
   /* if s <> 0, then some writers write b+1 blocks. */

   /* m_csteps is maximum number of writes */
   m_csteps = blocks_in_writer;
   if( r > 0 )
      m_csteps++;
   
/* 4. Set up communication data structures */
/*   m_isend[fd].m_maxbuf = bufsize;*/
//   if( nkg == 1 )
//      m_isend[fd].m_maxbuf = blsize*typsize*nig;
//   else
//      m_isend[fd].m_maxbuf = blsize*typsize*nig*njg;

   try{
      m_isend.m_steps = m_csteps;
      m_isend.m_has_values = true;
      m_isend.m_ncomm = new int[m_csteps];
      m_isend.m_comm_id = new int*[m_csteps];
      for( p = 0 ; p < 6 ; p++ )
	 m_isend.m_comm_index[p] = new int*[m_csteps];
      // Substep communication
      m_isend.m_subcommlabel = new int*[m_csteps];
      //   unused for send
      m_isend.m_nsubcomm = NULL;
      m_isend.m_subcomm = NULL;
      for( p = 0 ; p < m_csteps ; p++ )
      {
	 m_isend.m_subcommlabel[p] = NULL;
	 m_isend.m_comm_id[p] = NULL;
	 for( int p2=0 ; p2 < 6 ; p2++ )
	    m_isend.m_comm_index[p2][p] = NULL;
      }
   }
   catch( bad_alloc &ba )
   {
      cout << "Parallel_IO::init_array, processor " << gproc <<  ". Initial allocation of m_isend failed "
	      << " csteps = " << m_csteps << " Exception= " << ba.what() << endl;
      MPI_Abort(MPI_COMM_WORLD,0);
   }
      // Initialize pointers to nil
   for( int p1=0 ; p1 < m_csteps ; p1++ )
   {
      m_isend.m_comm_id[p1] = NULL;
      for( int p2= 0 ; p2 < 6 ; p2++ )
	 m_isend.m_comm_index[p2][p1] = NULL;
   }

   int nglast=nkg, nlast=nk, olast=ok; 
   if( nkg == 1 )
   {
      nglast = njg;
      nlast = nj;
      olast = oj;
   }

   //   if( myid == 0 )
   //   {
   //      cout << " init_array " << nglast << " " << nlast << " " << olast << endl;
   //      cout << "            " << nblocks << " " << blsize << " " << m_csteps << endl;
   //   }

// Count the number of sends
   for( b = 1 ; b <= m_csteps ; b++ )
   {
      nsend = 0;
      for( p = 1 ; p <= m_nwriters ; p++ )
      {
	 if( p <= r )
            blnr = (p-1)*(blocks_in_writer+1)+b;
	 else
	    blnr = r*(blocks_in_writer+1)+(p-r-1)*blocks_in_writer+b;
         if( blnr <= nblocks )
	 {
	    if(  blnr <= s )
	    {
	       kb = (blnr-1)*(blsize+1)+1;
	       ke = blnr*(blsize+1);
	    }
	    else
	    {
	       kb = s*(blsize+1) + (blnr-s-1)*(blsize)+1;
	       ke = s*(blsize+1) + (blnr-s)*(blsize);	       
	    }
	    if( kb <= 1 )
	       kb = 1;
	    if( ke >= nglast )
	       ke = nglast;
	 // intersect my array patch [1+oi,ni+oi]x[1+oj,nj+oj]x..
	 // with the block in writer p, [1..nig]x[1..njg]x[kb,ke]
	    kbl = 1+olast;
	    kel = nlast+olast;
            if( kbl > 1 )
	       kbl += padding;
	    if( kel < nglast )
	       kel -= padding;
	    if( !(kel<kb || kbl>ke) )
	    {
	       nsend++;
	    }
	 }
      }
      m_isend.m_ncomm[b-1] = nsend;
      if( nsend > 0 )
      {
         try
	 {
	    m_isend.m_comm_id[b-1]  = new int[nsend];
	    for( p = 0 ; p < 6 ; p++ )
	       m_isend.m_comm_index[p][b-1] = new int[nsend];
	 }
	 catch( bad_alloc& ba )
	 {
	    cout << "Parallel_IO::init_array, processor " << gproc <<  
	       ". Allocation of m_isend.m_comm_id or m_comm_index failed " 
		 << " nsend = " << nsend << " b= " << b << " Exception= " << ba.what() << endl;
	    MPI_Abort(MPI_COMM_WORLD,0);
	 }
      }
   }

/* Setup send information */
   maxpts = 0;
   for( b = 1 ; b <= m_csteps ; b++ )
   { 
      nsend = 0;
      for( p = 1 ; p <= m_nwriters ; p++ )
      {
	 if( p <= r )
            blnr = (p-1)*(blocks_in_writer+1)+b;
	 else
	    blnr = r*(blocks_in_writer+1)+(p-r-1)*blocks_in_writer+b;
         if( blnr <= nblocks )
	 {
	    if(  blnr <= s )
	    {
	       kb = (blnr-1)*(blsize+1)+1;
	       ke = blnr*(blsize+1);
	    }
	    else
	    {
	       kb = s*(blsize+1) + (blnr-s-1)*(blsize)+1;
	       ke = s*(blsize+1) + (blnr-s)*(blsize);	       
	    }
	    if( kb <= 1 )
	       kb = 1;
	    if( ke >= nglast )
	       ke = nglast;

	 // intersect my array patch [1+lgmap(1),ni+lgmap(1)]x[1+lgmap(2),nj+lgmap(2)]x..
	 // with the block in writer p, [1..nig]x[1..njg]x[kb,ke]

	    ibl = 1+oi;
	    iel = ni+oi;
	    jbl = 1+oj;
	    jel = nj+oj;
	    kbl = 1+ok;
	    kel = nk+ok;
            if( ibl > 1 )
	       ibl += padding;
	    if( iel < nig )
	       iel -= padding;
            if( jbl > 1 )
	       jbl += padding;
	    if( jel < njg )
	       jel -= padding;
            if( kbl > 1 )
	       kbl += padding;
	    if( kel < nkg )
	       kel -= padding;

	    if( nkg > 1 && !(kel<kb || kbl>ke) )
	    {
	       m_isend.m_comm_index[0][b-1][nsend] = ibl;
	       m_isend.m_comm_index[1][b-1][nsend] = iel;
	       m_isend.m_comm_index[2][b-1][nsend] = jbl;
	       m_isend.m_comm_index[3][b-1][nsend] = jel;
	       if( kbl < kb )
		  kbl = kb;
	       if( kel > ke )
		  kel = ke;
	       m_isend.m_comm_index[4][b-1][nsend] = kbl;
	       m_isend.m_comm_index[5][b-1][nsend] = kel;
	       m_isend.m_comm_id[b-1][nsend] = m_writer_ids[p-1];
	       nsend++;
               npts = (iel-ibl+1)*static_cast<size_t>((jel-jbl+1))*(kel-kbl+1);
               maxpts = maxpts > npts ? maxpts : npts;
	    }
            else if( nkg == 1 && !(jel<kb || jbl>ke) )
	    {
	       m_isend.m_comm_index[0][b-1][nsend] = ibl;
	       m_isend.m_comm_index[1][b-1][nsend] = iel;
	       if( jbl < kb )
		  jbl = kb;
	       if( jel > ke )
		  jel = ke;
	       m_isend.m_comm_index[2][b-1][nsend] = jbl;
	       m_isend.m_comm_index[3][b-1][nsend] = jel;
	       m_isend.m_comm_index[4][b-1][nsend] = kbl;
	       m_isend.m_comm_index[5][b-1][nsend] = kel;
	       m_isend.m_comm_id[b-1][nsend] = m_writer_ids[p-1];
	       nsend++;
               npts = (iel-ibl+1)*static_cast<size_t>((jel-jbl+1))*(kel-kbl+1);
               maxpts = maxpts > npts ? maxpts : npts;
	    }
	 }
      }
   }
   m_isend.m_maxbuf = maxpts;

   tag2 = 336;
   /* Senders pass info to receievers */
   if( m_iwrite == 1 )
   {
      try
      {
	 m_irecv.m_steps = m_csteps;
	 m_irecv.m_has_values = true;
	 m_irecv.m_ncomm   = new int[m_csteps];
	 m_irecv.m_comm_id = new int*[m_csteps];
	 for( p= 0 ; p < 6 ; p++ )
	    m_irecv.m_comm_index[p] = new int*[m_csteps];
	 m_irecv.m_ilow    = new int[m_csteps];
	 m_irecv.m_jlow    = new int[m_csteps];
	 m_irecv.m_klow    = new int[m_csteps];
	 m_irecv.m_niblock = new int[m_csteps];
	 m_irecv.m_njblock = new int[m_csteps];
	 m_irecv.m_nkblock = new int[m_csteps];
	 // Substep communication
         m_irecv.m_nsubcomm = new int[m_csteps];
         m_irecv.m_subcomm = new int*[m_csteps];
	 // Unused for receive
         m_irecv.m_subcommlabel = NULL;
	 for( int p1=0 ; p1 < m_csteps ; p1++ )
	 {
	    m_irecv.m_comm_id[p1] = NULL;
	    m_irecv.m_subcomm[p1] = NULL;
	    for( int p2= 0 ; p2 < 6 ; p2++ )
	       m_irecv.m_comm_index[p2][p1] = NULL;
	 }
      }
      catch( bad_alloc& ba )
      {
	 cout << "Parallel_IO::init_array, processor " << gproc <<  ". Initial allocation of m_irecv failed "
	      << " csteps = " << m_csteps << " Exception= " << ba.what() << endl;
	 MPI_Abort(MPI_COMM_WORLD,0);
      }

      // Initialize pointers to nil, and initialize
      for( int p1=0 ; p1 < m_csteps ; p1++ )
      {
	 m_irecv.m_comm_id[p1] = NULL;
         for( int p2= 0 ; p2 < 6 ; p2++ )
	    m_irecv.m_comm_index[p2][p1] = NULL;
	 m_irecv.m_ilow[p1] = 0;
	 m_irecv.m_jlow[p1] = 0;
	 m_irecv.m_klow[p1] = 0;
	 m_irecv.m_niblock[p1] = -1;
	 m_irecv.m_njblock[p1] = -1;
	 m_irecv.m_nkblock[p1] = -1;
	 m_irecv.m_maxbuf = 0;
      }
   }
   maxpts = 0;

   //   MPI_Barrier(m_data_comm);
   //   int fd=open("dbginfo.bin", O_CREAT | O_TRUNC | O_WRONLY, 0660 );
   //   MPI_Barrier( m_data_comm );

   //   if( myid == 0 )
   //      cout << "init_array point 4, csteps= " << m_csteps << endl;

   for( b = 1; b <= m_csteps ; b++ )
   {
      for( p = 1 ; p <= m_nwriters ; p++ )
      {
	 tag2 = 336 + b-1 + m_csteps*(p-1);
	 MPI_Request req;
         found = -1;
         if( m_isend.m_ncomm[b-1] > 0 )
	 {
            for( i = 0 ; i < m_isend.m_ncomm[b-1] ; i++ )
	       if( m_isend.m_comm_id[b-1][i] == m_writer_ids[p-1] )
	       {
		  found = i;
		  break;
	       }
	 }
	 // Note: nrecvs only touched by m_writer_ids[p-1], hence do not need to allocate
	 //       it in non-writer procs.
         retcode = MPI_Gather( &found, 1, MPI_INT, nrecvs, 1, MPI_INT, m_writer_ids[p-1], m_data_comm );
	 if( retcode != MPI_SUCCESS )
	 {
	    cout << "Parallel_IO::init_array, error from call to MPI_Gather. "
	      << "Return code = " << retcode << " from processor " << gproc << endl;
	 }
         if( found != -1  )
	 {
            v[0] = m_isend.m_comm_index[0][b-1][found];
            v[1] = m_isend.m_comm_index[1][b-1][found];
            v[2] = m_isend.m_comm_index[2][b-1][found];
            v[3] = m_isend.m_comm_index[3][b-1][found];
            v[4] = m_isend.m_comm_index[4][b-1][found];
            v[5] = m_isend.m_comm_index[5][b-1][found];
	    if( myid != m_writer_ids[p-1] )
	    {
	       //	       retcode = MPI_Send( v, 6, MPI_INT, m_writer_ids[p-1], tag2, m_data_comm );
	       retcode = MPI_Isend( v, 6, MPI_INT, m_writer_ids[p-1], tag2, m_data_comm, &req );
	       if( retcode != MPI_SUCCESS )
	       {
		  cout << "Parallel_IO::init_array, error from call to MPI_Send. "
		       << "Return code = " << retcode << " from processor " << gproc << endl;
	       }

	    }
	 }

	 //	 MPI_Barrier(m_data_comm);
	 //	 if( myid == 0 )
	 //	    cout << "init_array point 4c, b= " << b << " p= " << p << endl;

	 if( m_writer_ids[p-1] == myid )
	 {
	    j = 0;
            for( i=0 ; i < nprocs ; i++ )
	       if( nrecvs[i] > -1 )
	       {
		  j++;
	       }
	    m_irecv.m_ncomm[b-1] = j;
	    if( j > 0 )
	    {
               try
	       {
		  m_irecv.m_comm_id[b-1] = new int[j];
		  l = 0;
		  for( i=0 ; i < nprocs ; i++ )
		  {
		     if( nrecvs[i]>-1)
			m_irecv.m_comm_id[b-1][l++] = i;
		  }
	       // l should be j here 
		  for( q = 0 ; q < 6 ; q++ )
		     m_irecv.m_comm_index[q][b-1] = new int[j];
	       }
	       catch( bad_alloc& ba )
	       {
		  cout << "Parallel_IO::init_array, processor " << gproc <<  
		     ". Allocation of m_irecv.m_comm_id or m_comm_index failed " 
		       << " j = " << j << " b= " << b << " Exception= " << ba.what() << endl;
		  MPI_Abort(MPI_COMM_WORLD,0);
	       }
	       lims[0] = nig+1;
	       lims[1] = -1;
	       lims[2] = njg+1;
	       lims[3] = -1;
	       lims[4] = nkg+1;
	       lims[5] = -1;
	       //	       if( b == 10 && p == 77 )
	       //	       {
	       //		  cout << " preparing to receive from " << j << " processors :" << endl;
	       //		  for( i=0 ; i < j ; i++ )
	       //		     cout << i << " from proc no " << m_irecv.m_comm_id[b-1][i] << endl;
	       //	       }
	       for( i=0 ; i<j ; i++ )
	       {
                  if( myid != m_irecv.m_comm_id[b-1][i] )
		  {
		     retcode = MPI_Recv( vr, 6, MPI_INT, m_irecv.m_comm_id[b-1][i], tag2, m_data_comm, &status );
		     if( retcode != MPI_SUCCESS )
		     {
			cout << "Parallel_IO::init_array, error from call to MPI_Recv. "
			     << "Return code = " << retcode << " from processor " << gproc << endl;
		     }
		  }
		  else
		  {
                     for( l=0 ; l < 6 ; l++ )
			vr[l] = v[l];
		  }
		  m_irecv.m_comm_index[0][b-1][i] = vr[0];
		  m_irecv.m_comm_index[1][b-1][i] = vr[1];
		  m_irecv.m_comm_index[2][b-1][i] = vr[2];
		  m_irecv.m_comm_index[3][b-1][i] = vr[3];
		  m_irecv.m_comm_index[4][b-1][i] = vr[4];
		  m_irecv.m_comm_index[5][b-1][i] = vr[5];
		  if( vr[0] < lims[0] )
		     lims[0] = vr[0];
		  if( vr[2] < lims[2] )
		     lims[2] = vr[2];
		  if( vr[4] < lims[4] )
		     lims[4] = vr[4];
		  if( vr[1] > lims[1] )
		     lims[1] = vr[1];
		  if( vr[3] > lims[3] )
		     lims[3] = vr[3];
		  if( vr[5] > lims[5] )
		     lims[5] = vr[5];
	       }
	       if( nkg > 1 )
	       {
		  //	       m_irecv.m_ilow[b-1] = lims[0];
		  //	       m_irecv.m_jlow[b-1] = lims[2];
		  m_irecv.m_ilow[b-1] = 1;
		  m_irecv.m_jlow[b-1] = 1;
		  m_irecv.m_klow[b-1] = lims[4];
		  //	       m_irecv.m_niblock[b-1] = lims[1]-lims[0]+1;
		  //	       m_irecv.m_njblock[b-1] = lims[3]-lims[2]+1;
		  m_irecv.m_niblock[b-1] = nig;
		  m_irecv.m_njblock[b-1] = njg;
		  m_irecv.m_nkblock[b-1] = lims[5]-lims[4]+1;
	       //               npts = (lims[1]-lims[0]+1)*(lims[3]-lims[2]+1)*(lims[5]-lims[4]+1);
		  npts = nig*static_cast<size_t>(njg)*(lims[5]-lims[4]+1);
	       }
	       else if( nkg == 1 )
	       {
		  m_irecv.m_ilow[b-1] = 1;
		  m_irecv.m_jlow[b-1] = lims[2];
		  m_irecv.m_klow[b-1] = lims[4];
		  m_irecv.m_niblock[b-1] = nig;
		  m_irecv.m_njblock[b-1] = lims[3]-lims[2]+1;
		  m_irecv.m_nkblock[b-1] = lims[5]-lims[4]+1;
		  npts = nig*static_cast<size_t>(lims[3]-lims[2]+1)*(lims[5]-lims[4]+1);
	       }
               maxpts = maxpts > npts ? maxpts : npts;
	    }
	 }
	 //         lseek(fd,myid*sizeof(int),SEEK_SET);
	 //         nrw=write(fd,&tag2,sizeof(int));
	 //	 if( MPI_Wtime()-t0 > flim )
	 //	 {
	 //	    fsync(fd);
	 //	    flim = flim + 500;
	 //	 }
	 //	 //	 MPI_Barrier(m_data_comm);
	 //	 //	 if( myid == 0 )
	 //	 //	    cout << "init_array point 4d, b= " << b << " p= " << p << endl;
      }
   }
   //   close(fd);

   if( m_iwrite == 1 )
      m_irecv.m_maxiobuf = maxpts;
   else
      m_irecv.m_maxiobuf = 0;

   //   MPI_Barrier(m_data_comm);
   //   if( myid == 0 )
   //      cout << "init_array point 5" << endl;

   //   if( m_iwrite == 1 )
   //      delete[] nrecvs;

   setup_substeps();

   //   if( m_iwrite == 1 )
   //   {
   //      size_t maxptsbz=0;
   //      for( b = 1; b <= m_csteps ; b++ )
   //      {
   //	 size_t bsize = 0;
   //	 for( i = 0  ; i < m_irecv.m_ncomm[b-1] ; i++ )
   //	 {
   //	    vr[0] = m_irecv.m_comm_index[0][b-1][i];
   //	    vr[1] = m_irecv.m_comm_index[1][b-1][i];
   //	    vr[2] = m_irecv.m_comm_index[2][b-1][i];
   //	    vr[3] = m_irecv.m_comm_index[3][b-1][i];
   //	    vr[4] = m_irecv.m_comm_index[4][b-1][i];
   //	    vr[5] = m_irecv.m_comm_index[5][b-1][i];
   //	    bsize += (vr[1]-vr[0]+1)*(vr[3]-vr[2]+1)*(vr[5]-vr[4]+1);
   //	    //	    if( (vr[1]-vr[0]+1)*(vr[3]-vr[2]+1)*(vr[5]-vr[4]+1) > bsizemax )
   //	    //	    {
   //	    //	       maxcomm = i;
   //	    //	       bsizemax = (vr[1]-vr[0]+1)*(vr[3]-vr[2]+1)*(vr[5]-vr[4]+1);
   //	    //	    }
   //	 }
   //	 maxptsbz = maxptsbz > bsize ? maxptsbz:bsize;
   //      }
      //      if( gproc == 263 )
      //      {
      //	 cout << "maxpts   = " << maxpts   << endl;
      //	 cout << "maxptsbz = " << maxptsbz << endl;
      //	 for( b=1 ; b<= m_csteps ; b++ )
      //	    cout << "b= " << b << " ncomm = " << m_irecv.m_ncomm[b-1] << endl;
	 //	 cout << "max bsize at i= "  << i << endl;
	 //	 cout << "   box= \n" ;
	 //	    for( b= 0 ; b < 6; b++ )
	 //	       cout << m_irecv.m_comm_index[b][0][maxcomm]  << " ";
	 //	    cout << endl;
      //      }
   //   m_irecv.m_maxbuf = maxpts;
      //      m_irecv.m_maxbuf = maxptsbz;
      //      if( maxptsbz > maxpts )
      //	 m_irecv.m_maxbuf = maxptsbz;
      //      else
      //	 m_irecv.m_maxbuf = maxpts;
      //      cout << "Maxpts buffer old " << maxpts << " new " << maxptsbz << endl;
      //      m_irecv.m_maxbuf = maxptsbz;

   if( m_iwrite == 1 )
   {
      int bufsize1, bufsize2;
      MPI_Allreduce(&(m_irecv.m_maxbuf),  &bufsize1,1,MPI_INT,MPI_MAX,m_write_comm);
      MPI_Allreduce(&(m_irecv.m_maxiobuf),&bufsize2,1,MPI_INT,MPI_MAX,m_write_comm);
// should check verbose value before printing this...
      // if( myid == m_writer_ids[0] )
      // 	 cout << "Parallel_IO::init_array maxiobuf = " << bufsize2 << " maxbuf = " << bufsize1 << endl;
      delete[] nrecvs;
   }
   }
}

//-----------------------------------------------------------------------
void Parallel_IO::setup_substeps( )
{
   size_t buflim = static_cast<size_t>(m_irecv.m_maxiobuf*1.3);
//   size_t buflim = 1512;

//   int fd=open("dbginfo.bin", O_CREAT | O_TRUNC | O_WRONLY, 0660 );
//   MPI_Barrier( m_data_comm );

//   int myid;
//   MPI_Comm_rank( m_data_comm, &myid );

   // 1. Make sure that bufsize is large enough to fit the largest 
   //    single message.
   if( m_iwrite )
   {
      size_t mxblock=0;
      for( int b=0; b < m_csteps ; b++ )
      {
	 if( m_irecv.m_ncomm[b]>0 )
	 {
	    for( int i=0 ; i < m_irecv.m_ncomm[b]; i++ )
	    {
	       int i1 = m_irecv.m_comm_index[0][b][i];
	       int i2 = m_irecv.m_comm_index[1][b][i];
	       int j1 = m_irecv.m_comm_index[2][b][i];
	       int j2 = m_irecv.m_comm_index[3][b][i];
	       int k1 = m_irecv.m_comm_index[4][b][i];
	       int k2 = m_irecv.m_comm_index[5][b][i];
	       int nri = i2-i1+1;
	       int nrj = j2-j1+1;
	       int nrk = k2-k1+1;
	       size_t blsize = nri*((size_t)nrj)*nrk ;
	       mxblock = blsize > mxblock ? blsize : mxblock;
	    }
	 }
      }
      if( buflim < mxblock )
      {
	 cout << "Parallel_IO::setup_substeps, Error: buflim must be at least " << mxblock
	      << " current value = " << buflim << endl;
	 cout << "....adjusting buflim to " << mxblock << "...done" << endl;
	 buflim=mxblock;
      }
   }

   //   if( myid == 0 )
   //      cout << "setup_substeps point 1 "  << endl;

   for( int b=0; b < m_csteps ; b++ )
   {
      int tag=665+b;
      MPI_Status status;
      MPI_Request* req;
      int* rbuf;
      // Count the number of substeps
      if( m_iwrite && m_irecv.m_ncomm[b]>0 )
      {
	 int nsub=1;
	 //	 size_t mxblock=0;
	 size_t size=0;

	 for( int i=0 ; i < m_irecv.m_ncomm[b]; i++ )
	 {
	    int i1 = m_irecv.m_comm_index[0][b][i];
	    int i2 = m_irecv.m_comm_index[1][b][i];
	    int j1 = m_irecv.m_comm_index[2][b][i];
	    int j2 = m_irecv.m_comm_index[3][b][i];
	    int k1 = m_irecv.m_comm_index[4][b][i];
	    int k2 = m_irecv.m_comm_index[5][b][i];
	    int nri = i2-i1+1;
	    int nrj = j2-j1+1;
	    int nrk = k2-k1+1;
	    size_t blsize = nri*((size_t)nrj)*nrk ;
	    //	    mxblock = blsize > mxblock ? blsize : mxblock;
	    if( size + blsize > buflim )
	    {
	       nsub++;
	       size = blsize;
	    }
	    else
	       size += blsize;
	 }
	 m_irecv.m_nsubcomm[b] = nsub;
	 // Set up intervalls for the substeps
	 m_irecv.m_subcomm[b] = new int[nsub+1];
	 size=0;
	 m_irecv.m_subcomm[b][0] = 0;
	 int s=1;
	 for( int i=0 ; i < m_irecv.m_ncomm[b]; i++ )
	 {
	    int i1 = m_irecv.m_comm_index[0][b][i];
	    int i2 = m_irecv.m_comm_index[1][b][i];
	    int j1 = m_irecv.m_comm_index[2][b][i];
	    int j2 = m_irecv.m_comm_index[3][b][i];
	    int k1 = m_irecv.m_comm_index[4][b][i];
	    int k2 = m_irecv.m_comm_index[5][b][i];
	    int nri = i2-i1+1;
	    int nrj = j2-j1+1;
	    int nrk = k2-k1+1;
	    if( size + nri*((size_t)nrj)*nrk > buflim )
	    {
	       m_irecv.m_subcomm[b][s] = i;
	       size = nri*((size_t)nrj)*nrk;
	       s++;
	    }
	    else
	       size += nri*((size_t)nrj)*nrk;
	 }
	 m_irecv.m_subcomm[b][s] = m_irecv.m_ncomm[b];
      
      // Communicate the step id to non-io processors
	 rbuf = new int[m_irecv.m_ncomm[b]];
         req  = new MPI_Request[m_irecv.m_ncomm[b]];
	 for( int ss= 0 ; ss < nsub ; ss++ )
	 {
	    for( int i= m_irecv.m_subcomm[b][ss] ; i < m_irecv.m_subcomm[b][ss+1] ; i++ )
	    {
	       rbuf[i] = ss;
	       MPI_Isend( &rbuf[i], 1, MPI_INT, m_irecv.m_comm_id[b][i], tag, m_data_comm, &req[i] );
	    }
	 }
      }
      else if( m_iwrite && m_irecv.m_ncomm[b]==0 )
      {
	 m_irecv.m_nsubcomm[b] = 0;
      }
      //      if( myid == 0 )
      //	 cout << "setup_substeps point 2, b= "  << b << " out of " << m_csteps << " steps " << endl;
  
      // Receive the step id
      if( m_isend.m_ncomm[b] > 0 )
	 m_isend.m_subcommlabel[b] = new int[m_isend.m_ncomm[b]];

      //      lseek(fd,myid*sizeof(int),SEEK_SET);
      //      int one=1;
      //      size_t nr=write(fd,&one,sizeof(int));

      for( int i = 0 ; i < m_isend.m_ncomm[b] ; i++ )
      {
	 int subid;
	 MPI_Recv( &subid, 1, MPI_INT, m_isend.m_comm_id[b][i], tag, m_data_comm, &status );
	 m_isend.m_subcommlabel[b][i] = subid;
      }
      //      lseek(fd,myid*sizeof(int),SEEK_SET);
      //      int two=2;
      //      nr=write(fd,&two,sizeof(int));

      if( m_iwrite == 1 )
	 for( int i = 0 ; i < m_irecv.m_ncomm[b]; i++ )
	    MPI_Wait( &req[i], &status );

      //      lseek(fd,myid*sizeof(int),SEEK_SET);
      //      int three=3;
      //      nr=write(fd,&three,sizeof(int));

      if( m_iwrite == 1 && m_irecv.m_ncomm[b]>0 )
      {
	 delete[] rbuf;
	 delete[] req;
      }
   }
   //   close(fd);
   // Compute exact bufsize needed, should now be < buflim.
   //   if( myid == 0 )
   //      cout << "setup_substeps point 3 "  << endl;

   if( m_iwrite )
   {
      size_t maxbuf=0;
      for( int b = 0 ; b < m_csteps ; b++ )
      {
	 for( int ss= 0 ; ss < m_irecv.m_nsubcomm[b]; ss++ )
	 {
	    size_t subbuf=0;
	    for( int i= m_irecv.m_subcomm[b][ss] ; i < m_irecv.m_subcomm[b][ss+1] ; i++ )
	    {
	       int i1 = m_irecv.m_comm_index[0][b][i];
	       int i2 = m_irecv.m_comm_index[1][b][i];
	       int j1 = m_irecv.m_comm_index[2][b][i];
	       int j2 = m_irecv.m_comm_index[3][b][i];
	       int k1 = m_irecv.m_comm_index[4][b][i];
	       int k2 = m_irecv.m_comm_index[5][b][i];
	       int nri = i2-i1+1;
	       int nrj = j2-j1+1;
	       int nrk = k2-k1+1;
	       size_t blsize = nri*((size_t)nrj)*nrk ;
	       subbuf += blsize;
	    }
	    maxbuf = maxbuf > subbuf ? maxbuf : subbuf;
	 }
      }
      m_irecv.m_maxbuf = maxbuf;
   }
}

//-----------------------------------------------------------------------
void Parallel_IO::write_array( int* fid, int nc, void* array, off_t pos0,
			       char* typ )
{
//
//  Write array previously set up by constructing object.
// Input: fid - File descriptor, obtained by calling open.
//        nc  - Number of components per grid point of array.
//        array - The data array, local in the processor
//        pos0  - Start writing the array at this byte position in file.
//        typ   - Declared type of 'array', possible values are "float" or "double".
//                Note the array on disk will be written with this type.
//
   int i1, i2, j1, j2, k1, k2, nsi, nsj, nsk, nri, nrj, nrk;
   int b, i, mxsize, ii, jj, kk, c, niblock, njblock, nkblock;
   int il, jl, kl, tag, myid, retcode, gproc;
   off_t ind, ptr, sizew;
   MPI_Status status;
   MPI_Request* req;
   double* rbuf, *ribuf;
   float* rfbuf, *ribuff;
   bool debug =false;
   if( m_data_comm != MPI_COMM_NULL )
   {
      MPI_Comm_rank( MPI_COMM_WORLD, &gproc );
      float* arf;
      double* ar;
      double* sbuf;
      float*  sbuff;
      int flt, typsize;
      try
      {
	 if( strcmp(typ,"float")==0)
	 {
	    arf = static_cast<float*>(array);
	    sbuff = new float[m_isend.m_maxbuf*nc];
	    flt = 1;
	    typsize = sizeof(float);
	 }
	 else if( strcmp(typ,"double")==0 )
	 {
	    ar = static_cast<double*>(array);
	    sbuf  = new double[m_isend.m_maxbuf*nc];
	    typsize = sizeof(double);
	    flt = 0;
	 }
	 else
	 {
	 // error return
	 }
      }
      catch( bad_alloc& ba )
      {
	 int gproc;
	 MPI_Comm_rank( MPI_COMM_WORLD, &gproc );
	 cout << "Parallel_IO::write_array, processor " << gproc <<  
	    ". Allocation of sbuf or sbuff failed. Tried to allocate " << m_isend.m_maxbuf*nc;
	 if( flt == 0 )
	    cout << "doubles";
	 else
	    cout << "floats";
	 cout << " Exception= " << ba.what() << endl;
	 MPI_Abort(MPI_COMM_WORLD,0);
      }

      if( debug )
      {
	 cout << "DEBUGPARIO iwrite= " << m_iwrite << endl;
	 m_isend.print(0);
	 if( m_iwrite==1 )
	    m_irecv.print(1);
      }
      bool really_writing;
      if( m_iwrite == 1 )
      {
         really_writing = false;
         for( b=0 ; b < m_csteps ; b++ )
	    if( m_irecv.m_ncomm[b] > 0 )
	       really_writing = true;
      }

      MPI_Comm_rank( m_data_comm, &myid );

      if( m_iwrite == 1 && really_writing )
      {
	 try
	 {
	    if( flt == 1 )
	    {
	       rfbuf  = new float[m_irecv.m_maxiobuf*nc];
	       ribuff = new float[m_irecv.m_maxbuf*nc];
	    }
	    else
	    {
	       rbuf  = new double[m_irecv.m_maxiobuf*nc];
	       ribuf = new double[m_irecv.m_maxbuf*nc];
	    }
	 }
	 catch( bad_alloc& ba )
	 {
	    int gproc;
	    MPI_Comm_rank( MPI_COMM_WORLD, &gproc );
	    cout << "Parallel_IO::write_array, processor " << gproc <<  
	    ". Allocation of rbuf and ribuff failed. Tried to allocate " << m_irecv.m_maxbuf*nc
		 << " + "  << m_irecv.m_maxiobuf*nc;
	    if( flt ==  0 )
	       cout  << " doubles";
	    else
	       cout  << " floats";
	    cout  <<  ". Exception= " << ba.what() << endl;
	    MPI_Abort(MPI_COMM_WORLD,0);
	 }
	 mxsize = 0;
	 for( b= 0; b < m_csteps ; b++ )
	    if( mxsize < m_irecv.m_ncomm[b] )
	       mxsize = m_irecv.m_ncomm[b];
	 try
	 {
	    req = new MPI_Request[mxsize];
	 }
	 catch( bad_alloc& ba )
	 {
	    int gproc;
	    MPI_Comm_rank( MPI_COMM_WORLD, &gproc );
	    cout << "Parallel_IO::write_array, processor " << gproc << ". Allocating req failed. " 
		 << "Tried to allocate " << mxsize << " MPI_Requests. " << "Exception = " << ba.what() << endl;
	    MPI_Abort(MPI_COMM_WORLD,0);
	 }

	 il = m_irecv.m_ilow[0];
	 jl = m_irecv.m_jlow[0];
	 kl = m_irecv.m_klow[0];
	 ind = il-1+nig*(jl-1)+((off_t)nig)*njg*(kl-1);
	 sizew = lseek( *fid, pos0+nc*ind*typsize, SEEK_SET );
	 if( sizew == -1 )
	 {
	    int eno = errno;
	    cout << "Error in write_array: could not go to write start position" << endl;
	    if( eno == EBADF )
	       cout << "errno = EBADF" << endl;
	    if( eno == EINVAL )
	       cout << "errno = EINVAL" << endl;
	    if( eno == EOVERFLOW )
	       cout << "errno = EOVERFLOW" << endl;
	    if( eno == ESPIPE )
	       cout << "errno = ESPIPE" << endl;
	    cout << "errno = " << eno << endl;
            cout << "Requested offset = " << pos0+nc*ind*typsize << endl;
            cout << "pos0 = " << pos0 << endl;
	    cout << "nc = " << nc << endl;
	    cout << "ind = " << ind << endl;
	    cout << "typsize = " << typsize << endl;
            cout << "m_csteps = " << m_csteps << endl;
	    cout << "nglobal = " << nig << " " << njg << " " << nkg << endl;
	    cout << "m_irecv.m_ilow " << m_irecv.m_ilow[0] << endl;
	    cout << "m_irecv.m_jlow " << m_irecv.m_jlow[0] << endl;
	    cout << "m_irecv.m_klow " << m_irecv.m_klow[0] << endl;
            cout << "m_irecv.m_ncomm[0] = " << m_irecv.m_ncomm[0] << endl;
	    //	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
      }

      tag = 334;
      for( b = 0; b < m_csteps ; b++ )
      {
      // Post receive
	 if( m_iwrite == 1 )
	 {
	    ptr = 0;
	    for( i = 0  ; i < m_irecv.m_ncomm[b] ; i++ )
	    {
	       i1 = m_irecv.m_comm_index[0][b][i];
	       i2 = m_irecv.m_comm_index[1][b][i];
	       j1 = m_irecv.m_comm_index[2][b][i];
	       j2 = m_irecv.m_comm_index[3][b][i];
	       k1 = m_irecv.m_comm_index[4][b][i];
	       k2 = m_irecv.m_comm_index[5][b][i];
	       nri = i2-i1+1;
	       nrj = j2-j1+1;
	       nrk = k2-k1+1;
               if( flt == 0 )
		  retcode = MPI_Irecv( ribuf+ptr, nri*nrj*nrk*nc, MPI_DOUBLE, m_irecv.m_comm_id[b][i],
			     tag, m_data_comm, &req[i] );
	       else
		  retcode = MPI_Irecv( ribuff+ptr, nri*nrj*nrk*nc, MPI_FLOAT, m_irecv.m_comm_id[b][i],
			     tag, m_data_comm, &req[i] );
               if( retcode != MPI_SUCCESS )
	       {
		  cout << "Parallel_IO::write_array, error from call to MPI_Irecv. "
		       << "Return code = " << retcode << " from processor " << gproc << endl;
	       }
	       ptr += ((off_t)nri)*nrj*nrk*nc;
	    }
	 }
      // Send 
	 for( i = 0 ; i < m_isend.m_ncomm[b] ; i++ )
	 {
	    i1 = m_isend.m_comm_index[0][b][i];
	    i2 = m_isend.m_comm_index[1][b][i];
	    j1 = m_isend.m_comm_index[2][b][i];
	    j2 = m_isend.m_comm_index[3][b][i];
	    k1 = m_isend.m_comm_index[4][b][i];
	    k2 = m_isend.m_comm_index[5][b][i];
	    nsi = i2-i1+1;
	    nsj = j2-j1+1;
	    nsk = k2-k1+1;
            if( flt == 0 )
	    {
	       for( kk=k1 ; kk <= k2 ; kk++ )
		  for( jj=j1 ; jj <= j2 ; jj++ )
		     for( ii=i1 ; ii <= i2 ; ii++ )
			for( c=0 ; c < nc ; c++ )
			{
			   sbuf[c+nc*(ii-i1)+nc*nsi*(jj-j1)+nc*nsi*nsj*(kk-k1)]
			      = ar[c+nc*(ii-1-oi)+ni*nc*(jj-1-oj)+((off_t)ni)*nj*nc*(kk-1-ok)];
			}
	       retcode = MPI_Send( sbuf, nsi*nsj*nsk*nc, MPI_DOUBLE, m_isend.m_comm_id[b][i], tag, m_data_comm );
	    }
	    else
	    {
	       for( kk=k1 ; kk <= k2 ; kk++ )
		  for( jj=j1 ; jj <= j2 ; jj++ )
		     for( ii=i1 ; ii <= i2 ; ii++ )
			for( c=0 ; c < nc ; c++ )
			{
			   sbuff[c+nc*(ii-i1)+nc*nsi*(jj-j1)+nc*nsi*nsj*(kk-k1)]
			      = arf[c+nc*(ii-1-oi)+ni*nc*(jj-1-oj)+((off_t)ni)*nj*nc*(kk-1-ok)];
			}
	       retcode = MPI_Send( sbuff, nsi*nsj*nsk*nc, MPI_FLOAT, m_isend.m_comm_id[b][i], tag, m_data_comm );
	    }
	    if( retcode != MPI_SUCCESS )
	    {
	       cout << "Parallel_IO::write_array, error from call to MPI_Send. "
		       << "Return code = " << retcode << " from processor " << gproc << endl;
	    }

	 }

      // Do actual receive
	 if( m_iwrite == 1 && m_irecv.m_ncomm[b] > 0 )
	 {
	    ptr = 0;
	    il = m_irecv.m_ilow[b];
	    jl = m_irecv.m_jlow[b];
	    kl = m_irecv.m_klow[b];
	    niblock = m_irecv.m_niblock[b];
	    njblock = m_irecv.m_njblock[b];
	    nkblock = m_irecv.m_nkblock[b];
	    for( i = 0  ; i < m_irecv.m_ncomm[b] ; i++ )
	    {
	       retcode = MPI_Wait( &req[i], &status );
               if( retcode != MPI_SUCCESS )
	       {
		  cout << "Parallel_IO::write_array, error from call to MPI_Wait. "
		       << "Return code = " << retcode << " from processor " << gproc << endl;
	       }
	       i1 = m_irecv.m_comm_index[0][b][i];
	       i2 = m_irecv.m_comm_index[1][b][i];
	       j1 = m_irecv.m_comm_index[2][b][i];
	       j2 = m_irecv.m_comm_index[3][b][i];
	       k1 = m_irecv.m_comm_index[4][b][i];
	       k2 = m_irecv.m_comm_index[5][b][i];

	       nri = i2-i1+1;
	       nrj = j2-j1+1;
	       nrk = k2-k1+1;

	       if( flt == 0 )
	       {
		  double* recbuf = ribuf+ptr;
		  for( kk=k1 ; kk <= k2 ; kk++ )
		     for( jj=j1 ; jj <= j2 ; jj++ )
			for( ii=i1 ; ii <= i2 ; ii++ )
			   for( c=0 ; c < nc ; c++ )
			   {
			      rbuf[c+nc*(ii-il)+nc*niblock*(jj-jl)+nc*((off_t)niblock)*njblock*(kk-kl)]
				 = recbuf[c+nc*(ii-i1)+nri*nc*(jj-j1)+nri*nrj*nc*(kk-k1)];
			   }
	       }
	       else
	       {
		  float* recbuf = ribuff+ptr;
		  for( kk=k1 ; kk <= k2 ; kk++ )
		     for( jj=j1 ; jj <= j2 ; jj++ )
			for( ii=i1 ; ii <= i2 ; ii++ )
			   for( c=0 ; c < nc ; c++ )
			   {
			      rfbuf[c+nc*(ii-il)+nc*niblock*(jj-jl)+nc*((off_t)niblock)*njblock*(kk-kl)]
				 = recbuf[c+nc*(ii-i1)+nri*nc*(jj-j1)+nri*nrj*nc*(kk-k1)];
			   }
	       }
	       ptr += ((off_t)nri)*nrj*nrk*nc;
	    }

// Write to disk
	    begin_sequential( m_write_comm );
	    if( flt == 0 )
	    {
	       sizew = write( *fid, rbuf, sizeof(double)*((off_t)nc)*niblock*njblock*nkblock );
               if( sizew != sizeof(double)*((off_t)nc)*niblock*njblock*nkblock )
	       {
                  cout << "Error in write_array: could not write requested array size";
		  cout << "  requested "<< sizeof(double)*((off_t)nc)*niblock*njblock*nkblock << " bytes\n";
		  cout << "  written "<< sizew << " bytes\n";
	          MPI_Abort(MPI_COMM_WORLD,1);
	       }
	    }
	    else
	    {
	       sizew = write( *fid, rfbuf, sizeof(float)*((off_t)nc)*niblock*njblock*nkblock );
	       if( sizew != sizeof(float)*((off_t)nc)*niblock*njblock*nkblock )
	       {
                  int eno = errno;
		  if( eno == EAGAIN )
		     cout << "errno = EAGAIN" << endl;
		  if( eno == EBADF )
		     cout << "errno = EBADF" << endl;
                  if( eno == EFBIG )
		     cout << "errno = EFBIG" << endl;
                  if( eno == EINTR )
		     cout << "errno = EINTR" << endl;
                  if( eno == EINVAL )
		     cout << "errno = EINVAL" << endl;
                  if( eno == EIO )
		     cout << "errno = EIO" << endl;
                  if( eno == ENOSPC )
		     cout << "errno = ENOSPC" << endl;
                  if( eno == EPIPE )
		     cout << "errno = EPIPE" << endl;
                  cout << "errno = " << eno << endl;
                  cout << "Error in write_array: could not write requested array size";
		  cout << "  requested "<< sizeof(float)*((off_t)nc)*niblock*njblock*nkblock << " bytes\n";
		  cout << "  written "<< sizew << " bytes\n";
	          MPI_Abort(MPI_COMM_WORLD,1);
	       }
	    }
	    end_sequential( m_write_comm );
// Is this really needed ?
// Need to sync before throwing away rbuf/rfbuf in next communication step
//	    fsync(*fid);

	 }
      }
      if( flt == 0 )
	 delete[] sbuf;
      else
	 delete[] sbuff;

      if( m_iwrite == 1 && really_writing )
      {
	 if( flt == 0 )
	 {
	    delete[] rbuf;
	    delete[] ribuf;
	 }
	 else
	 {
	    delete[] rfbuf;
	    delete[] ribuff;
	 }
	 delete[] req;
      }
   }
}

//-----------------------------------------------------------------------
void Parallel_IO::read_array( int* fid, int nc, float_sw4* array, off_t pos0,
			      const char* typ, bool swap_bytes )
{
//  Read array previously set up by constructing object.
// Input: fid   - File descriptor, obtained by calling open (before calling this function).
//        nc    - Number of components per grid point of array.
//        array - The data array, local in the processor
//        pos0  - Start reading the array at this byte position in file.
//        typ   - Type of array stored on file, possible values are "float" or "double".
//       Note: internal buffer uses double irrespective of the type of float_sw4,
//             however, this would have been a problem only if data on disk were stored with
//             higher precision than double.
   int i1, i2, j1, j2, k1, k2, nsi, nsj, nsk, nri, nrj, nrk;
   int b, i, mxsize, ii, jj, kk, c, niblock, njblock, nkblock;
   int il, jl, kl, tag, myid, retcode, gproc, s;
   size_t ind, ptr, sizew;
   off_t sizelseek;
   MPI_Status status;
   MPI_Request* req;

   if( m_data_comm != MPI_COMM_NULL )
   {
      MPI_Comm_rank( MPI_COMM_WORLD, &gproc );
      double* rbuf, *ribuf;
      float* rfbuf;
      double* sbuf;
      int flt, typsize;
      if( strcmp(typ,"float") == 0 )
      {
	 flt = 1;
         typsize = sizeof(float);
      }
      else if( strcmp(typ,"double")==0 )
      {
         flt = 0;
         typsize = sizeof(double);
      }
      else
      {
	 // error return
      }
      try{
	 sbuf  = new double[m_isend.m_maxbuf*nc];
      }
      catch( bad_alloc &ba )
      {
	 cout << "Parallel_IO::read_array, processor " << gproc << ". Allocating memory for sbuf failed. "
	      << "Tried to allocate " << m_isend.m_maxbuf*nc << " doubles. "  
	      <<  "Exception= " << ba.what() << endl;
	 MPI_Abort(MPI_COMM_WORLD,0);
      }
      bool really_reading=false;
      if( m_iwrite == 1 )
      {
         really_reading = false;
         for( b=0 ; b < m_csteps ; b++ )
	    if( m_irecv.m_ncomm[b] > 0 )
	       really_reading = true;
      }
      MPI_Comm_rank( m_data_comm, &myid );

      //      if( gproc == 263 )
      //      {
      //	 m_irecv.print(1);
      //      }

      if( m_iwrite == 1 && really_reading )
      {
	 try {
	    if( flt == 1 )
	       rfbuf  = new float[m_irecv.m_maxiobuf*nc];
	    else
	       rbuf  = new double[m_irecv.m_maxiobuf*nc];
	    ribuf = new double[m_irecv.m_maxbuf*nc];
	 }
	 catch( bad_alloc &ba )
	 {
	    cout << "Parallel_IO::read_array, processor " << gproc << ". Allocating memory for rbuf and ribuf failed. "
		 << "Tried to allocate " << m_irecv.m_maxiobuf*nc;
	    if( flt ==  0 )
	       cout  << " doubles ";
	    else
	       cout  << " floats ";
	    cout << "and " << m_irecv.m_maxbuf*nc << " doubles.";
	    cout  <<  " Exception= " << ba.what() << endl;
	    MPI_Abort(MPI_COMM_WORLD,0);
	 }
	 mxsize = 0;
	 for( b= 0; b < m_csteps ; b++ )
	    if( mxsize < m_irecv.m_ncomm[b] )
	       mxsize = m_irecv.m_ncomm[b];
	 try{
	    req = new MPI_Request[mxsize];
	 }
	 catch( bad_alloc &ba )
	 {
	    cout << "Parallel_IO::read_array, processor " << gproc << ". Allocating req failed. " 
		 << "Tried to allocate " << mxsize << " MPI_Requests. " << "Exception = " << ba.what() << endl;
	    MPI_Abort(MPI_COMM_WORLD,0);
	 }

         int bb=0;
         while( m_irecv.m_ncomm[bb] == 0 )
	    bb++;
	 
	 il = m_irecv.m_ilow[bb];
	 jl = m_irecv.m_jlow[bb];
	 kl = m_irecv.m_klow[bb];
	 
	 ind = il-1+((off_t)nig)*(jl-1)+((off_t)nig)*njg*(kl-1);
	 sizelseek = lseek( *fid, pos0+nc*ind*typsize, SEEK_SET );
	 if( sizelseek == -1 )
	 {
	    int eno = errno;
	    cout << "Error in read_array: could not go to read start position" << endl;
	    if( eno == EBADF )
	       cout << "errno = EBADF" << endl;
	    if( eno == EINVAL )
	       cout << "errno = EINVAL" << endl;
	    if( eno == EOVERFLOW )
	       cout << "errno = EOVERFLOW" << endl;
	    if( eno == ESPIPE )
	       cout << "errno = ESPIPE" << endl;
	    cout << "errno = " << eno << endl;
            cout << "Requested offset = " << pos0+nc*ind*typsize << endl;
            cout << "pos0 = " << pos0 << endl;
	    cout << "nc = " << nc << endl;
	    cout << "ind = " << ind << endl;
	    cout << "typsize = " << typsize << endl;
            cout << "m_csteps = " << m_csteps << endl;
	    cout << "nglobal = " << nig << " " << njg << " " << nkg << endl;
	    cout << "m_irecv.m_ilow " << m_irecv.m_ilow[0] << endl;
	    cout << "m_irecv.m_jlow " << m_irecv.m_jlow[0] << endl;
	    cout << "m_irecv.m_klow " << m_irecv.m_klow[0] << endl;
            cout << "m_irecv.m_ncomm[0] = " << m_irecv.m_ncomm[0] << endl;
	    //	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
      }

      tag = 334;
      for( b = 0; b < m_csteps ; b++ )
      {
	 if( m_iwrite == 1 )
	 {
	    ptr = 0;
	    il = m_irecv.m_ilow[b];
	    jl = m_irecv.m_jlow[b];
	    kl = m_irecv.m_klow[b];
	    niblock = m_irecv.m_niblock[b];
	    njblock = m_irecv.m_njblock[b];
      	    nkblock = m_irecv.m_nkblock[b];

// Read from disk
	    begin_sequential( m_write_comm );
	    if( m_irecv.m_ncomm[b] > 0 )
	    {
	       if( flt == 0 )
	       {
		  sizew = read( *fid, rbuf, sizeof(double)*nc*((size_t)niblock)*njblock*nkblock );
		  if( sizew != sizeof(double)*nc*((size_t)niblock)*njblock*nkblock )
		  {
		     cout << "Error in read_array: could not read requested array size";
		     cout << "  requested "<< sizeof(double)*((off_t)nc)*niblock*njblock*nkblock << " bytes\n";
		     cout << "  read "<< sizew << " bytes\n";
		  }
		  if( swap_bytes )
		     m_bswap.byte_rev( rbuf, nc*((size_t)niblock)*njblock*nkblock, "double");
	       }
	       else
	       {
		  sizew = read( *fid, rfbuf, sizeof(float)*nc*((size_t)niblock)*njblock*nkblock );
		  if( sizew != sizeof(float)*nc*((size_t)niblock)*njblock*nkblock )
		  {
		     cout << "Error in read_array: could not read requested array size";
		     cout << "  requested "<< sizeof(float)*((off_t)nc)*niblock*njblock*nkblock << " bytes\n";
		     cout << "  read "<< sizew << " bytes\n";
		  }
		  if( swap_bytes )
		     m_bswap.byte_rev( rfbuf, nc*((size_t)niblock)*njblock*nkblock, "float");
	       }
	    }
	    end_sequential( m_write_comm );
	 }
// Hand out to other processors
//	    if( m_iwrite == 1 && m_irecv.m_ncomm[b] > 0 )
//	    {
	 int mxstep=0;
	 if( m_iwrite )
	    mxstep = m_irecv.m_nsubcomm[b];
	 for( i= 0 ; i < m_isend.m_ncomm[b] ; i++ )
	    mxstep = mxstep > (m_isend.m_subcommlabel[b][i]+1) ? mxstep : (m_isend.m_subcommlabel[b][i]+1);
	 
//	 for( s = 0 ; s < m_irecv.m_nsubcomm[b] ; s++ )
	 for( s = 0 ; s < mxstep ; s++ )
	 {
	    if( m_iwrite == 1 && m_irecv.m_ncomm[b] > 0 && s < m_irecv.m_nsubcomm[b] )
	    {
	       ptr = 0;
	       for( i = m_irecv.m_subcomm[b][s]  ; i < m_irecv.m_subcomm[b][s+1] ; i++ )
	       {
		  i1 = m_irecv.m_comm_index[0][b][i];
		  i2 = m_irecv.m_comm_index[1][b][i];
		  j1 = m_irecv.m_comm_index[2][b][i];
		  j2 = m_irecv.m_comm_index[3][b][i];
		  k1 = m_irecv.m_comm_index[4][b][i];
		  k2 = m_irecv.m_comm_index[5][b][i];
		  nri = i2-i1+1;
		  nrj = j2-j1+1;
		  nrk = k2-k1+1;
		  double* recbuf = ribuf+ptr;
		  if( flt == 0 )
		  {
		     for( kk=k1 ; kk <= k2 ; kk++ )
			for( jj=j1 ; jj <= j2 ; jj++ )
			   for( ii=i1 ; ii <= i2 ; ii++ )
			      for( c=0 ; c < nc ; c++ )
			      {
				 recbuf[c+nc*(ii-i1)+nri*nc*(jj-j1)+nri*nrj*nc*(kk-k1)]=
				    rbuf[c+nc*(ii-il)+nc*niblock*(jj-jl)+nc*niblock*njblock*(kk-kl)];
			      }
		  }
		  else
		  {
		     for( kk=k1 ; kk <= k2 ; kk++ )
			for( jj=j1 ; jj <= j2 ; jj++ )
			   for( ii=i1 ; ii <= i2 ; ii++ )
			      for( c=0 ; c < nc ; c++ )
			      {
				 recbuf[c+nc*(ii-i1)+nri*nc*(jj-j1)+nri*nrj*nc*(kk-k1)] =
				    rfbuf[c+nc*(ii-il)+nc*niblock*(jj-jl)+nc*niblock*njblock*(kk-kl)];
			      }
		  }
	    //            printf("%d sending %d step %d to %d, size %d %d %d\n",myid,i,b,m_irecv[fd].m_comm_id[b][i],nri,nrj,nrk);
		  retcode = MPI_Isend( recbuf, nri*nrj*nrk*nc, MPI_DOUBLE, m_irecv.m_comm_id[b][i],
			  tag, m_data_comm, &req[i] );
		  if( retcode != MPI_SUCCESS )
		  {
		     cout << "Parallel_IO::read_array, error calling MPI_Isend, "
			  << "return code = " << retcode << " from processor " << gproc << endl;
		  }
		  ptr += nri*((size_t)nrj)*nrk*nc;
	       }
	    }
      // Do actual receive
	    for( i = 0 ; i < m_isend.m_ncomm[b] ; i++ )
	    {
	       if( m_isend.m_subcommlabel[b][i] == s )
	       {
		  i1 = m_isend.m_comm_index[0][b][i];
		  i2 = m_isend.m_comm_index[1][b][i];
		  j1 = m_isend.m_comm_index[2][b][i];
		  j2 = m_isend.m_comm_index[3][b][i];
		  k1 = m_isend.m_comm_index[4][b][i];
		  k2 = m_isend.m_comm_index[5][b][i];
		  nsi = i2-i1+1;
		  nsj = j2-j1+1;
		  nsk = k2-k1+1;
		  retcode = MPI_Recv( sbuf, nsi*nsj*nsk*nc, MPI_DOUBLE,
				      m_isend.m_comm_id[b][i], tag, m_data_comm, &status );
		  if( retcode != MPI_SUCCESS )
		  {
		     cout << "Parallel_IO::read_array, error calling MPI_Recv, "
			  << "return code = " << retcode << " from processor " << gproc << endl;
		  }
		  for( kk=k1 ; kk <= k2 ; kk++ )
		     for( jj=j1 ; jj <= j2 ; jj++ )
			for( ii=i1 ; ii <= i2 ; ii++ )
			   for( c=0 ; c < nc ; c++ )
			   {
			      array[c+nc*(ii-1-oi)+ni*nc*(jj-1-oj)+ni*nj*nc*(kk-1-ok)] =
				 sbuf[c+nc*(ii-i1)+nc*nsi*(jj-j1)+nc*nsi*nsj*(kk-k1)];
			   }
	       }
	    }
	    if( m_iwrite == 1 && s < m_irecv.m_nsubcomm[b])
	       for( i = m_irecv.m_subcomm[b][s]  ; i < m_irecv.m_subcomm[b][s+1] ; i++ )
	       {
		  retcode = MPI_Wait( &req[i], &status );
		  if( retcode != MPI_SUCCESS )
		  {
		     cout << "Parallel_IO::read_array, error calling MPI_Wait, "
			  << "return code = " << retcode << " from processor " << gproc << endl;
		  }
	       }
	 }
      }
      delete[] sbuf;

      if( m_iwrite == 1 && really_reading )
      {
	 delete[] req;
	 delete[] ribuf;
	 if( flt == 1 )
	    delete[] rfbuf;
	 else
	    delete[] rbuf;
      }
   }
}

//-----------------------------------------------------------------------
void Parallel_IO::begin_sequential( MPI_Comm comm )
{
   if( m_parallel_file_system == 0 )
   {
      int mtag, slask, myid;
      MPI_Status status;
      mtag = 10;
      MPI_Comm_rank( comm, &myid );
      if( myid > 0 )
	 MPI_Recv( &slask, 1, MPI_INT, myid-1, mtag, comm, &status );
   }
}

//-----------------------------------------------------------------------
void Parallel_IO::end_sequential( MPI_Comm comm )
{
   if( m_parallel_file_system == 0 )
   {
      int mtag, slask=0, myid, nproc;
      mtag = 10;
      MPI_Comm_rank( comm, &myid );
      MPI_Comm_size( comm, &nproc );
      if( myid < nproc-1 )
	 MPI_Send( &slask, 1, MPI_INT, myid+1, mtag, comm );
   }
}

//-----------------------------------------------------------------------
void Parallel_IO::writer_barrier( )
{
   if( m_iwrite == 1 )
      MPI_Barrier( m_write_comm );
}

//-----------------------------------------------------------------------
void Parallel_IO::print( )
{
   int myid, mydid, mywid;
   MPI_Comm_rank( MPI_COMM_WORLD, &myid );
   cout << myid << " printing " << endl;
   if( m_data_comm != MPI_COMM_NULL )
   {
      MPI_Comm_rank( m_data_comm, &mydid );
      cout << "past first " << endl;
      if( m_iwrite )
	 MPI_Comm_rank( m_write_comm, &mywid );
      else
	 mywid = -1;
      cout << "ID in world " << myid << endl;
      cout << "ID in data comm " << mydid << endl;
      cout << "ID in writer " << mywid << endl;
      cout << "iwrite = " << m_iwrite << " local sizes " << ni << " " << nj << " " << nk << endl;
      cout << " in global space [" << 1+oi << "," << ni+oi << "]x[" << 1+oj << "," << nj+oj << "]x["
	   << 1+ok << "," << nk+ok << "]" << endl;
      cout << "send info: " << endl;
      m_isend.print(0);
      if( m_iwrite )
      {
	 cout << "Recv info: " << endl;
	 m_irecv.print(1);
      }
   }
}

//-----------------------------------------------------------------------
int Parallel_IO::proc_zero()
// One unique processor in m_write_comm, to write file headers etc.
{
   int retval=0;
   if( m_write_comm != MPI_COMM_NULL )
   {
      int myid;
      MPI_Comm_rank( m_write_comm, &myid );
      //      cout << "PIO::Proc_zero myid= " << myid << " m_writer_ids[0]= " << m_writer_ids[0] << endl;
      //      if( myid == m_writer_ids[0] )
      if( myid == 0 )
	 retval = 1;
   }
   return retval;
}

