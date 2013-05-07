#include <mpi.h>
#include <iostream>
using namespace std;

#include <cstring>
#include <errno.h>

#include "Parallel_IO.h"

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
}

//-----------------------------------------------------------------------
void Comminfo::print( int recv )
{
   if( m_has_values )
   {
      cout << "maxbuf= " << m_maxbuf << endl;
      cout << "steps= " << m_steps << endl;
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
	 }
      }	 
      if( recv == 1 )
      {
	 for( int s=0 ; s < m_steps ; s++ )
	 {
            cout << " step " << s << endl;
	    cout << "(i,j,k) wr block " << m_ilow[s] << " " << m_jlow[s] << " " << m_klow[s] <<endl;
	    cout << "size wr block    " << m_niblock[s] << " " << m_njblock[s] << " " << m_nkblock[s] <<endl;
	 }
      }
   }
}

//-----------------------------------------------------------------------
Parallel_IO::Parallel_IO( int iwrite, int pfs, int globalsizes[3], int localsizes[3],
		  int starts[3], int nptsbuf, int padding )
{
   int ihave_array=1;
   m_data_comm = m_write_comm = MPI_COMM_NULL;
   if( localsizes[0] < 1 || localsizes[1] < 1 || localsizes[2] < 1 )
      ihave_array = 0;
   init_pio( iwrite, pfs, ihave_array );
   //   cout << "gsizes " << globalsizes[0] <<  " " << globalsizes[1] << " " << globalsizes[2] << endl;
   //   cout << "lsizes " << localsizes[0] <<  " " << localsizes[1] << " " << localsizes[2] << endl;
   //   cout << "ssizes " << starts[0] <<  " " << starts[1] << " " << starts[2] << endl;
   init_array( globalsizes, localsizes, starts, nptsbuf, padding );
   //   m_irecv.print(1);
}

//-----------------------------------------------------------------------
void Parallel_IO::init_pio( int iwrite, int pfs, int ihave_array )
// Initialize for parallel I/O, 
// Input: iwrite - 0 this processor will not participate in I/O op.
//                 1 this processor will participate in I/O op.
//        pfs    - 0 I/O will be done to non-parallel file system.
//                 1 I/O will be done to parallel file system.
//        ihave_array - 0 this processor holds no part of the array.
//                      1 this processor holds some part of the array.
//                     -1 (default) the array is present in all processors.
// Note: I/O processors are selected from the set of processors holding the array.
{
   int* tmp;
   int nprocs, p, i;
   MPI_Group world_group, writer_group, array_group;

   // 1. Create communicator of all procs. that hold part of the array
   // save as m_data_comm.
   if( ihave_array == -1 )
      m_data_comm = MPI_COMM_WORLD;
   else
   {
      MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
      tmp      = new int[nprocs];
      MPI_Allgather( &ihave_array, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD );
      int narray = 0;
      for( p = 0 ; p < nprocs ; p++ )
	 if( tmp[p] == 1 )
	    narray++;
      if( narray < 1 )
      {
	 // error 
      }
      int* array_holders = new int[narray];
      i = 0;
      for( p=0 ; p < nprocs ; p++ )
      {
	 if( tmp[p] == 1 )
	    array_holders[i++] = p;
      }
      delete[] tmp;
      MPI_Comm_group( MPI_COMM_WORLD, &world_group );
      MPI_Group_incl( world_group, narray, array_holders, &array_group );
      MPI_Comm_create( MPI_COMM_WORLD, array_group, &m_data_comm );
      MPI_Group_free( &world_group );
      MPI_Group_free( &array_group );
      delete[] array_holders;
   }

   if( m_data_comm != MPI_COMM_NULL )
   {
      MPI_Comm_size( m_data_comm, &nprocs );

   // 2. Create communicator of all I/O processors, save as m_write_comm
      m_iwrite = iwrite;
      tmp      = new int[nprocs];
      MPI_Allgather( &m_iwrite, 1, MPI_INT, tmp, 1, MPI_INT, m_data_comm );
      m_nwriters = 0;
      for( p = 0 ; p < nprocs ; p++ )
	 if( tmp[p] == 1 )
	    m_nwriters++;
      if( m_nwriters < 1 )
      {
         delete[] tmp;
      // error return
      }
      m_writer_ids = new int[m_nwriters];
      i = 0;
      for( p=0 ; p < nprocs ; p++ )
      {
	 if( tmp[p] == 1 )
	    m_writer_ids[i++] = p;
      }
      delete[] tmp;

      MPI_Comm_group( m_data_comm, &world_group );
      MPI_Group_incl( world_group, m_nwriters, m_writer_ids, &writer_group );
      MPI_Comm_create( m_data_comm, writer_group, &m_write_comm );
      MPI_Group_free( &world_group );
      MPI_Group_free( &writer_group );
   }

   // 3. Save parallel file system info
   m_parallel_file_system = pfs;
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
   int found, i, j, q, lims[6], v[6], vr[6], nprocs, tag, tag2, myid;
   int maxpts, npts;
   int* nrecvs;
   size_t nblocks;

   MPI_Status status;

   if( m_data_comm != MPI_COMM_NULL )
   {
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
      nrecvs = new int[nprocs];

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

   m_isend.m_steps = m_csteps;
   m_isend.m_has_values = true;
   m_isend.m_ncomm = new int[m_csteps];
   m_isend.m_comm_id = new int*[m_csteps];
   for( p = 0 ; p < 6 ; p++ )
      m_isend.m_comm_index[p] = new int*[m_csteps];
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
	 m_isend.m_comm_id[b-1]  = new int[nsend];
	 for( p = 0 ; p < 6 ; p++ )
	    m_isend.m_comm_index[p][b-1] = new int[nsend];
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
               npts = (iel-ibl+1)*(jel-jbl+1)*(kel-kbl+1);
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
               npts = (iel-ibl+1)*(jel-jbl+1)*(kel-kbl+1);
               maxpts = maxpts > npts ? maxpts : npts;
	    }
	 }
      }
   }
   m_isend.m_maxbuf = maxpts;

   tag = 335;
   tag2 = 336;
   /* Senders pass info to receievers */
   if( m_iwrite == 1 )
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
      // Initialize pointers to nil
      for( int p1=0 ; p1 < m_csteps ; p1++ )
      {
	 m_irecv.m_comm_id[p1] = NULL;
         for( int p2= 0 ; p2 < 6 ; p2++ )
	    m_irecv.m_comm_index[p2][p1] = NULL;
      }
   }
   maxpts = 0;

   for( b = 1; b <= m_csteps ; b++ )
   {
      for( p = 1 ; p <= m_nwriters ; p++ )
      {
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
         MPI_Gather( &found, 1, MPI_INT, nrecvs, 1, MPI_INT, m_writer_ids[p-1], m_data_comm );

         if( found != -1  )
	 {
            v[0] = m_isend.m_comm_index[0][b-1][found];
            v[1] = m_isend.m_comm_index[1][b-1][found];
            v[2] = m_isend.m_comm_index[2][b-1][found];
            v[3] = m_isend.m_comm_index[3][b-1][found];
            v[4] = m_isend.m_comm_index[4][b-1][found];
            v[5] = m_isend.m_comm_index[5][b-1][found];
	    if( myid != m_writer_ids[p-1] )
	       MPI_Send( v, 6, MPI_INT, m_writer_ids[p-1], tag2, m_data_comm );
	 }
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
	       lims[0] = nig+1;
	       lims[1] = -1;
	       lims[2] = njg+1;
	       lims[3] = -1;
	       lims[4] = nkg+1;
	       lims[5] = -1;
	       for( i=0 ; i<j ; i++ )
	       {
                  if( myid != m_irecv.m_comm_id[b-1][i] )
		     MPI_Recv( vr, 6, MPI_INT, m_irecv.m_comm_id[b-1][i], tag2, m_data_comm, &status );
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
	       m_irecv.m_ilow[b-1] = lims[0];
	       m_irecv.m_jlow[b-1] = lims[2];
	       m_irecv.m_klow[b-1] = lims[4];
	       m_irecv.m_niblock[b-1] = lims[1]-lims[0]+1;
	       m_irecv.m_njblock[b-1] = lims[3]-lims[2]+1;
	       m_irecv.m_nkblock[b-1] = lims[5]-lims[4]+1;
               npts = (lims[1]-lims[0]+1)*(lims[3]-lims[2]+1)*(lims[5]-lims[4]+1);
               maxpts = maxpts > npts ? maxpts : npts;
	    }
	 }
      }
   }
   m_irecv.m_maxbuf = maxpts;
   if( m_iwrite == 1 )
      delete[] nrecvs;
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
//
   int i1, i2, j1, j2, k1, k2, nsi, nsj, nsk, nri, nrj, nrk;
   int b, i, mxsize, ii, jj, kk, c, niblock, njblock, nkblock;
   int il, jl, kl, tag, myid;
   off_t ind, ptr, sizew;
   MPI_Status status;
   MPI_Request* req;
   double* rbuf, *ribuf;
   float* rfbuf, *ribuff;

   if( m_data_comm != MPI_COMM_NULL )
   {
      float* arf;
      double* ar;
      double* sbuf;
      float*  sbuff;
      int flt, typsize;
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
	 if( flt == 1 )
	 {
	    rfbuf  = new float[m_irecv.m_maxbuf*nc];
	    ribuff = new float[m_irecv.m_maxbuf*nc];
	 }
	 else
	 {
	    rbuf  = new double[m_irecv.m_maxbuf*nc];
	    ribuf = new double[m_irecv.m_maxbuf*nc];
	 }
	 mxsize = 0;
	 for( b= 0; b < m_csteps ; b++ )
	    if( mxsize < m_irecv.m_ncomm[b] )
	       mxsize = m_irecv.m_ncomm[b];
	 req = new MPI_Request[mxsize];

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
		  MPI_Irecv( ribuf+ptr, nri*nrj*nrk*nc, MPI_DOUBLE, m_irecv.m_comm_id[b][i],
			     tag, m_data_comm, &req[i] );
	       else
		  MPI_Irecv( ribuff+ptr, nri*nrj*nrk*nc, MPI_FLOAT, m_irecv.m_comm_id[b][i],
			     tag, m_data_comm, &req[i] );
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
	       MPI_Send( sbuf, nsi*nsj*nsk*nc, MPI_DOUBLE, m_isend.m_comm_id[b][i], tag, m_data_comm );
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
	       MPI_Send( sbuff, nsi*nsj*nsk*nc, MPI_FLOAT, m_isend.m_comm_id[b][i], tag, m_data_comm );
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
	       MPI_Wait( &req[i], &status );
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
void Parallel_IO::read_array( int* fid, int nc, double* array, off_t pos0,
			      char* typ )
{
//  Read array previously set up by constructing object.
// Input: fid - File descriptor, obtained by calling open.
//        nc  - Number of components per grid point of array.
//        array - The data array, local in the processor
//        pos0  - Start reading the array at this byte position in file.
//        typ   - Declared type of 'array', possible values are "float" or "double".
//
   int i1, i2, j1, j2, k1, k2, nsi, nsj, nsk, nri, nrj, nrk;
   int b, i, mxsize, ii, jj, kk, c, niblock, njblock, nkblock;
   int il, jl, kl, tag, myid;
   size_t ind, ptr, sizew;
   MPI_Status status;
   MPI_Request* req;

   if( m_data_comm != MPI_COMM_NULL )
   {
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
      sbuf  = new double[m_isend.m_maxbuf*nc];
      bool really_reading=false;
      if( m_iwrite == 1 )
      {
         really_reading = false;
         for( b=0 ; b < m_csteps ; b++ )
	    if( m_irecv.m_ncomm[b] > 0 )
	       really_reading = true;
      }
      MPI_Comm_rank( m_data_comm, &myid );

      if( m_iwrite == 1 && really_reading )
      {
	 if( flt == 1 )
	    rfbuf  = new float[m_irecv.m_maxbuf*nc];
	 else
	    rbuf  = new double[m_irecv.m_maxbuf*nc];
	 ribuf = new double[m_irecv.m_maxbuf*nc];
	 mxsize = 0;
	 for( b= 0; b < m_csteps ; b++ )
	    if( mxsize < m_irecv.m_ncomm[b] )
	       mxsize = m_irecv.m_ncomm[b];
	 req = new MPI_Request[mxsize];

	 il = m_irecv.m_ilow[0];
	 jl = m_irecv.m_jlow[0];
	 kl = m_irecv.m_klow[0];
	 ind = il-1+((off_t)nig)*(jl-1)+((off_t)nig)*njg*(kl-1);
	 sizew = lseek( *fid, pos0+nc*ind*typsize, SEEK_SET );
	 if( sizew == -1 )
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
	 if( m_iwrite == 1 && m_irecv.m_ncomm[b] > 0 )
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
	    if( flt == 0 )
	    {
	       sizew = read( *fid, rbuf, sizeof(double)*nc*((size_t)niblock)*njblock*nkblock );
	       if( sizew != sizeof(double)*nc*((size_t)niblock)*njblock*nkblock )
	       {
                  cout << "Error in read_array: could not read requested array size";
		  cout << "  requested "<< sizeof(double)*((off_t)nc)*niblock*njblock*nkblock << " bytes\n";
		  cout << "  read "<< sizew << " bytes\n";
	       }
	    }
	    else
	    {
	       sizew = read( *fid, rfbuf, sizeof(float)*nc*((size_t)niblock)*njblock*nkblock );
	       if( sizew != sizeof(float)*nc*((size_t)niblock)*njblock*nkblock )
	       {
                  cout << "Error in read_array: could not read requested array size";
		  cout << "  requested "<< sizeof(double)*((off_t)nc)*niblock*njblock*nkblock << " bytes\n";
		  cout << "  read "<< sizew << " bytes\n";
	       }
	    }
	    end_sequential( m_write_comm );

// Hand out to other processors
//	    if( m_iwrite == 1 && m_irecv.m_ncomm[b] > 0 )
//	    {
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
	       MPI_Isend( recbuf, nri*nrj*nrk*nc, MPI_DOUBLE, m_irecv.m_comm_id[b][i],
			  tag, m_data_comm, &req[i] );
	       ptr += nri*((size_t)nrj)*nrk*nc;
	    }
	 }
      // Do actual receive
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
	    //	    ni = m_dims[fd].ni;
	    //	    nj = m_dims[fd].nj;
	    MPI_Recv( sbuf, nsi*nsj*nsk*nc, MPI_DOUBLE,
		      m_isend.m_comm_id[b][i], tag, m_data_comm, &status );
	    for( kk=k1 ; kk <= k2 ; kk++ )
	       for( jj=j1 ; jj <= j2 ; jj++ )
		  for( ii=i1 ; ii <= i2 ; ii++ )
		     for( c=0 ; c < nc ; c++ )
		     {
			array[c+nc*(ii-1-oi)+ni*nc*(jj-1-oj)+ni*nj*nc*(kk-1-ok)] =
	                       sbuf[c+nc*(ii-i1)+nc*nsi*(jj-j1)+nc*nsi*nsj*(kk-k1)];
		     }
	 }
	 if( m_iwrite == 1 )
	    for( i = 0 ; i < m_irecv.m_ncomm[b] ; i++ )
	       MPI_Wait( &req[i], &status );
      }
      delete[] sbuf;

      if( m_iwrite == 1 && really_reading )
      {
	 delete[] req;
	 delete[] ribuf;
	 if( flt == 0 )
	    delete[] rbuf;
	 else
	    delete[] rfbuf;
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
      int mtag, slask, myid, nproc;
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
{
   int retval=0;
   if( m_write_comm != MPI_COMM_NULL )
   {
      int myid;
      MPI_Comm_rank( m_write_comm, &myid );
      if( myid == m_writer_ids[0] )
	 retval = 1;
   }
   return retval;
}

