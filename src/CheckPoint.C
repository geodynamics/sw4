#include "mpi.h"

#include "EW.h"
#include "Require.h"
#include "CheckPoint.h"
#include <fcntl.h>
#include <ctime>
#include <cstring>
#include <unistd.h>
#include <cstdio>

CheckPoint* CheckPoint::nil=static_cast<CheckPoint*>(0);

//-----------------------------------------------------------------------
// Null constructor, makes it possible to always query 'DoCheckPointing' and 'DoRestart'
CheckPoint::CheckPoint( EW* a_ew ) :
   mEW(a_ew),
   mWritingCycle(-1),
   mCycleInterval(0),
   mStartTime(0.0),
   mCheckPointFile(" "),
   mPreceedZeros(0),
   m_double(true),
   mRestartFile(" "),
   m_winallocated(false),
   m_bufsize(0),
   m_fileno(0),
   mDoCheckPointing(false),
   mRestartPathSet(false),
   mDoRestart(false),
   m_kji_order(true)
{

}

//-----------------------------------------------------------------------
// Save check point files, but no restart
CheckPoint::CheckPoint( EW* a_ew,
			int cycle, 
			int cycleInterval,
			string fname,
			size_t bufsize ) :
   mEW(a_ew),
   mWritingCycle(cycle),
   mCycleInterval(cycleInterval),
   mStartTime(0.0),
   mCheckPointFile(fname),
   mPreceedZeros(0),
   m_double(true),
   mRestartFile("restart"),
   m_winallocated(false),
   m_bufsize(bufsize),
   m_fileno(0),
   mDoCheckPointing(true),
   mRestartPathSet(false),
   mDoRestart(false),
   m_kji_order(true)
{
   m_double = sizeof(float_sw4)==sizeof(double);
}

//-----------------------------------------------------------------------
// Restart, but not initialized for saving check points
CheckPoint::CheckPoint( EW* a_ew, string fname, size_t bufsize ) :
   mEW(a_ew),
   mWritingCycle(-1),
   mCycleInterval(0),
   mStartTime(0.0),
   mCheckPointFile("chkpt"),
   mPreceedZeros(0),
   m_double(true),
   mRestartFile(fname),
   m_winallocated(false),
   m_bufsize(bufsize),
   m_fileno(0),
   mDoCheckPointing(false),
   mRestartPathSet(false),
   mDoRestart(true)
{
   m_double = sizeof(float_sw4)==sizeof(double);
}

//-----------------------------------------------------------------------
CheckPoint::~CheckPoint()
{
   if( m_winallocated )
      for( int g=0 ; g < mEW->mNumberOfGrids ; g++ )
      {
	 delete[] mWindow[g];
         delete[] mGlobalDims[g];
      }
}

//-----------------------------------------------------------------------
bool CheckPoint::do_checkpointing()
{
   return mDoCheckPointing;
}

//-----------------------------------------------------------------------
bool CheckPoint::do_restart()
{
   return mDoRestart;
}

//-----------------------------------------------------------------------
void CheckPoint::setup_sizes( )
{
   if( mDoRestart || mDoCheckPointing )
   {
      if( !m_winallocated )
      {
	 mWindow.resize(mEW->mNumberOfGrids);
	 mGlobalDims.resize(mEW->mNumberOfGrids);
	 for( int g=0; g < mEW->mNumberOfGrids ; g++ )
	 {
	    mWindow[g] = new int[6];
	    mGlobalDims[g] = new int[6];
	 }
	 m_winallocated = true;
      }

      m_ihavearray.resize( mEW->mNumberOfGrids );
         
      int ghost_points = mEW->getNumberOfGhostPoints();
// tmp
//      printf("CheckPoint: Number of ghost points = %d\n", ghost_points);
      
      for( int g=0 ; g < mEW->mNumberOfGrids ; g++ )
      {
// With attenuation, the memory variables must be saved at all ghost points. This is because the memory variables satisfy
// ODEs at all points, and do not satify any boundary conditions 

         if (mEW->getLocalBcType(g, 0) == bProcessor)
            mWindow[g][0] = mEW->m_iStartInt[g];
         else
            mWindow[g][0] = mEW->m_iStartInt[g] - ghost_points;;

         if (mEW->getLocalBcType(g, 1) == bProcessor)
            mWindow[g][1] = mEW->m_iEndInt[g];
         else
            mWindow[g][1] = mEW->m_iEndInt[g] + ghost_points;
         
         if (mEW->getLocalBcType(g, 2) == bProcessor)
            mWindow[g][2] = mEW->m_jStartInt[g];
         else
            mWindow[g][2] = mEW->m_jStartInt[g] - ghost_points;;
            
         if (mEW->getLocalBcType(g, 3) == bProcessor)
            mWindow[g][3] = mEW->m_jEndInt[g];
         else
            mWindow[g][3] = mEW->m_jEndInt[g] + ghost_points;
            
// all points in k-dir are local to each proc
	 mWindow[g][4] = mEW->m_kStartInt[g]-ghost_points; // need 1 ghost point at MR interface
	 mWindow[g][5] = mEW->m_kEndInt[g]+ghost_points; // and all ghost points at bottom


// Need to store ghost point values outside the physical boundaries for the
// attenuation memory variables, because they satisfy ODEs and not BC

         mGlobalDims[g][0] = 1 - ghost_points;
         mGlobalDims[g][1] = mEW->m_global_nx[g] + ghost_points;
         mGlobalDims[g][2] = 1 - ghost_points;
         mGlobalDims[g][3] = mEW->m_global_ny[g] + ghost_points;
         // k-dir is local
	 mGlobalDims[g][4] = mWindow[g][4];
	 mGlobalDims[g][5] = mWindow[g][5];

 // The 3D array is assumed to span the entire computational domain
	 m_ihavearray[g] = true;
      }
      setSteps( mEW->getNumberOfTimeSteps() );
// tmp
      if (mEW->proc_zero())
         cout << "Checkpoint::setup_sizes: Calling define_pio()..." << endl;
      
      define_pio();
   } // end if doRestart || doCheckpointing
   
   //   cout << "mwind = " << mWindow[0][4] << " " << mWindow[0][5] << endl;
   //   cout << "globaldims = " << mGlobalDims[0][4] << " " << mGlobalDims[0][5] << endl;
}

//-----------------------------------------------------------------------
void CheckPoint::define_pio( )
{
   int glow = 0, ghigh = mEW->mNumberOfGrids;

  double time_start = MPI_Wtime();
  double time_measure[12];
  time_measure[0] = time_start;

   // Create the restart directory if it doesn't exist
   if( mRestartPathSet )
     mEW->create_directory(mRestartPath);

   m_parallel_io = new Parallel_IO*[ghigh-glow+1];
   for( int g=glow ; g < ghigh ; g++ )
   {
// tmp
      if (mEW->proc_zero())
         cout << "setup_sizes for grid g = " << g << endl;

      int global[3], local[3], start[3];
      for( int dim=0 ; dim < 3 ; dim++ )
      {
         global[dim] = mGlobalDims[g][2*dim+1]-mGlobalDims[g][2*dim]+1;
	 local[dim]  = mWindow[g][2*dim+1]-mWindow[g][2*dim]+1;
	 start[dim]  = mWindow[g][2*dim]-mGlobalDims[g][2*dim];
      }

      int iwrite = 0;
      int nrwriters = mEW->getNumberOfWritersPFS();
      int nproc=0, myid=0;
      MPI_Comm_size( MPI_COMM_WORLD, &nproc );
      MPI_Comm_rank( MPI_COMM_WORLD, &myid);

      // new hack 
      int* owners = new int[nproc];
      int i=0;
      for( int p=0 ; p<nproc ; p++ )
	 if( m_ihavearray[g] )
	    owners[i++] = p;
      if( nrwriters > i )
	 nrwriters = i;

      if( nrwriters > nproc )
	 nrwriters = nproc;
      int q, r;
      if( nproc == 1 || nrwriters == 1 )
      {
	 q = 0;
	 r = 0;
      }
      else
      {
	 q = (nproc-1)/(nrwriters-1);
	 r = (nproc-1) % (nrwriters-1);
      }
      for( int w=0 ; w < nrwriters ; w++ )
	 if( q*w+r == myid )
	    iwrite = 1;
//      std::cout << "Define PIO: grid " << g << " myid = " << myid << " iwrite= " << iwrite << " start= "
      //		<< start[0] << " " << start[1] << " " << start[2] << std::endl;
      if( m_kji_order )
      {
// Swap i and k on file
	 int tmp=global[0];
	 global[0]=global[2];
	 global[2]=tmp;
	 tmp=local[0];
	 local[0]=local[2];
	 local[2]=tmp;
	 tmp=start[0];
	 start[0]=start[2];
	 start[2]=tmp;
      }      
      if (mEW->proc_zero())
         cout << "Creating a Parallel_IO object for grid g = " << g << endl;
      m_parallel_io[g-glow] = new Parallel_IO( iwrite, mEW->usingParallelFS(), global, local, start, m_bufsize );
// tmp
      if (mEW->proc_zero())
         cout << "Done creating the Parallel_IO object" << endl;
      delete[] owners;
   }
}

//-----------------------------------------------------------------------
void CheckPoint::setSteps(int a_steps)
{
  char buffer[50];
  mPreceedZeros = snprintf(buffer, 50, "%d", a_steps );
}

//-----------------------------------------------------------------------
bool CheckPoint::timeToWrite( float_sw4 time, int cycle, float_sw4 dt )
{
   if( !mDoCheckPointing )
      return false;

   // Will we write at this time step (cycle) ?
   bool do_it=false;
   if( cycle == mWritingCycle )
      do_it = true;
   if( mCycleInterval !=  0 && cycle%mCycleInterval == 0 && time >= mStartTime )
      do_it = true;
   return do_it;
}

//-----------------------------------------------------------------------
void CheckPoint::compute_file_suffix( int cycle, std::stringstream& fileSuffix )
{
   fileSuffix << mCheckPointFile << ".cycle=";
   int temp = static_cast<int>(pow(10.0, mPreceedZeros - 1));
   int testcycle = cycle;
   if (cycle == 0)
      testcycle=1;
   while (testcycle < temp)
   {
      fileSuffix << "0";
      temp /= 10;
   }
   fileSuffix << cycle ;
   fileSuffix << ".sw4checkpoint";
}

//-----------------------------------------------------------------------
void CheckPoint::write_checkpoint( float_sw4 a_time, int a_cycle, vector<Sarray>& a_Um,
				   vector<Sarray>& a_U, vector<Sarray*>& a_AlphaVEm,
				   vector<Sarray*>& a_AlphaVE )
{
   //
   //File format: 
   //
   //    header (see routine write_header)
   //    for g=1,number of grids
   //       Um(g) (3 component float/double array)
   //       U(g)  (3 component float/double array)
   //       for m=1,number of mechanisms
   //           AlphaVEm(g,m) (3 component float/double array)
   //           AlphaVE(g,m)  (3 component float/double array)
   //       endfor
   //    endfor
   //

   //
   // Would it be possible to save the entire input file to the restart file ? 
   //
   int ng       = mEW->mNumberOfGrids;
//   off_t offset = (4+6*ng)*sizeof(int) + 2*sizeof(float_sw4);

   bool iwrite = false;
   for( int g=0 ; g < ng ; g++ )
      iwrite = iwrite || m_parallel_io[g]->i_write();

   std::stringstream s;
   if( iwrite )
   {
      std::stringstream fileSuffix;
      compute_file_suffix( a_cycle, fileSuffix );
      if ( mRestartPathSet )
	      s << mRestartPath << "/";
      else
      if( mEW->getPath() != "./" )
	      s << mEW->getPath() << "/";
      s << fileSuffix.str();
   }

 // Keep track of the number of files, save previous file name, and delete the second last.
   cycle_checkpoints( s.str() );

   // Open file from processor zero and write header.
   int hsize;
   int fid=-1;
   if( m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
      CHECK_INPUT(fid != -1, "CheckPoint::write_file: Error opening: " << s.str() );
      int myid;

      MPI_Comm_rank( MPI_COMM_WORLD, &myid );
      std::cout << "writing check point on file " << s.str() << " using " <<
	 m_parallel_io[0]->n_writers() << " writers" << std::endl;
      write_header( fid, a_time, a_cycle, hsize );
      fsync(fid); 
   }
   //   m_parallel_io[0]->writer_barrier();
   int bcast_root = m_parallel_io[0]->proc_zero_rank_in_comm_world();
   MPI_Bcast( &hsize, 1, MPI_INT, bcast_root, MPI_COMM_WORLD );
   off_t offset = hsize;

   // Open file from all writers
   if( iwrite && !m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_WRONLY );
      CHECK_INPUT(fid != -1, "CheckPoint::write_checkpoint:: Error opening: " << s.str() );
   }

   // Write data blocks
   char cprec[]="double";
   if( !m_double )
      strcpy(cprec,"float");
   for( int g = 0 ; g < ng ; g++ )
   {
      size_t npts = ((size_t)(mGlobalDims[g][1]-mGlobalDims[g][0]+1))*
 	            ((size_t)(mGlobalDims[g][3]-mGlobalDims[g][2]+1))*
	            ((size_t)(mGlobalDims[g][5]-mGlobalDims[g][4]+1));

      if( !mEW->usingParallelFS() || g == 0 )
	 m_parallel_io[g]->writer_barrier();
      
      size_t nptsloc = (size_t) (mWindow[g][1] - mWindow[g][0]+1)*
         (mWindow[g][3] - mWindow[g][2]+1)*
         (mWindow[g][5] - mWindow[g][4]+1);
      
      // allocate local buffer array
      float_sw4* doubleField = new float_sw4[3*nptsloc];

      if( m_kji_order )
      {
	 a_Um[g].extract_subarrayIK( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
				     doubleField );
	 m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	 offset += 3*npts*sizeof(float_sw4);

	 a_U[g].extract_subarrayIK( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
				    doubleField );
	 m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	 offset += 3*npts*sizeof(float_sw4);
	 for( int m=0 ; m < mEW->getNumberOfMechanisms() ; m++ )
	 {
	    a_AlphaVEm[g][m].extract_subarrayIK( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
					      doubleField );
	    m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	    offset += 3*npts*sizeof(float_sw4);

	    a_AlphaVE[g][m].extract_subarrayIK(  mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
                                              doubleField );
	    m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	    offset += 3*npts*sizeof(float_sw4);
	 }
      }
      else
      {
	 a_Um[g].extract_subarray( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
				  doubleField );
	 m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	 offset += 3*npts*sizeof(float_sw4);

	 a_U[g].extract_subarray( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
				 doubleField );
	 m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	 offset += 3*npts*sizeof(float_sw4);

	 for( int m=0 ; m < mEW->getNumberOfMechanisms() ; m++ )
	 {
	    a_AlphaVEm[g][m].extract_subarray( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
					      doubleField );
	    m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	    offset += 3*npts*sizeof(float_sw4);

	    a_AlphaVE[g][m].extract_subarray(  mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
                                              doubleField );
	    m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	    offset += 3*npts*sizeof(float_sw4);
	 }
      }
      delete[] doubleField;
   }
   if( iwrite )
      close(fid);
} // end write_checkpoint()

//-----------------------------------------------------------------------
void CheckPoint::read_checkpoint( float_sw4& a_time, int& a_cycle,
				  vector<Sarray>& a_Um, vector<Sarray>& a_U,
				  vector<Sarray*>& a_AlphaVEm, vector<Sarray*>& a_AlphaVE )
{
   //
   // It is assumed that the arrays are already declared with the right 
   // dimensions. This routine will check that sizes match, but will not
   // allocate or resize the arrays Um and U, AlphaVEm, AlphaVE
   //
   int ng       = mEW->mNumberOfGrids;
   //   off_t offset = (4+6*ng)*sizeof(int) + 2*sizeof(float_sw4);
   bool iread = false;
   for( int g=0 ; g < ng ; g++ )
      iread = iread || m_parallel_io[g]->i_write();

   std::stringstream s;
   if( iread )
   {
      if( mRestartPathSet )
	 s << mRestartPath << "/";
      else if( mEW->getPath() != "./" )
	 s << mEW->getPath();
      s << mRestartFile;
   }

   // Open file from processor zero and read header.
   int fid=-1;
   int hsize;
   if( m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_RDONLY ); 
      CHECK_INPUT(fid != -1, "CheckPoint::read_checkpoint: Error opening: " << s.str() );
      int myid;

      MPI_Comm_rank( MPI_COMM_WORLD, &myid );
      std::cout << "reading check point on file " << s.str() << endl;
      read_header( fid, a_time, a_cycle, hsize );
   }
   //   m_parallel_io[0]->writer_barrier();

   // Broadcast read information to all other processors.
   int bcast_root = m_parallel_io[0]->proc_zero_rank_in_comm_world();
   MPI_Bcast( &a_cycle, 1, MPI_INT,         bcast_root, MPI_COMM_WORLD );
   MPI_Bcast( &a_time,  1, mEW->m_mpifloat, bcast_root, MPI_COMM_WORLD );
   MPI_Bcast( &hsize, 1, MPI_INT, bcast_root, MPI_COMM_WORLD );
   off_t offset = hsize;

   // Open file from all readers
   if( iread && !m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_RDONLY );
      CHECK_INPUT(fid != -1, "CheckPoint::read_checkpoint:: Error opening file: " << s.str() );
   }

   // Read data blocks.
   char cprec[]="double";
   if( !m_double )
      strcpy(cprec,"float");

   for( int g = 0 ; g < ng ; g++ )
   {
      size_t npts = ((size_t)(mGlobalDims[g][1]-mGlobalDims[g][0]+1))*
	            ((size_t)(mGlobalDims[g][3]-mGlobalDims[g][2]+1))*
	            ((size_t)(mGlobalDims[g][5]-mGlobalDims[g][4]+1));

      // size_t nptsloc = ((size_t)(mEW->m_iEnd[g]-mEW->m_iStart[g]+1))*
      //                  ((size_t)(mEW->m_jEnd[g]-mEW->m_jStart[g]+1))*
      //                  ((size_t)(mEW->m_kEnd[g]-mEW->m_kStart[g]+1));

      size_t nptsloc = (size_t) (mWindow[g][1] - mWindow[g][0]+1)*
         (mWindow[g][3] - mWindow[g][2]+1)*
         (mWindow[g][5] - mWindow[g][4]+1);

      if( !mEW->usingParallelFS() || g == 0 )
	 m_parallel_io[g]->writer_barrier();      

      // array with ghost points in k-ONLY read into doubleField, 
      float_sw4* doubleField = new float_sw4[3*nptsloc];
      if( m_kji_order )
      {
	 m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, cprec );
	 offset += 3*npts*sizeof(float_sw4);
	 a_Um[g].insert_subarrayIK( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
				  doubleField );

	 m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, cprec );
	 offset += 3*npts*sizeof(float_sw4);
	 a_U[g].insert_subarrayIK( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
				 doubleField );

	 for( int m=0 ; m < mEW->getNumberOfMechanisms() ; m++ )
	 {
	    m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, cprec );
	    offset += 3*npts*sizeof(float_sw4);
	    a_AlphaVEm[g][m].insert_subarrayIK( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
					      doubleField );

	    m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, cprec );
	    offset += 3*npts*sizeof(float_sw4);
	    a_AlphaVE[g][m].insert_subarrayIK( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
					     doubleField );
	 }
      }
      else
      {
	 m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, cprec );
	 offset += 3*npts*sizeof(float_sw4);
	 a_Um[g].insert_subarray( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
				  doubleField );

	 m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, cprec );
	 offset += 3*npts*sizeof(float_sw4);
	 a_U[g].insert_subarray( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
				 doubleField );

	 for( int m=0 ; m < mEW->getNumberOfMechanisms() ; m++ )
	 {
	    m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, cprec );
	    offset += 3*npts*sizeof(float_sw4);
	    a_AlphaVEm[g][m].insert_subarray( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
					      doubleField );

	    m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, cprec );
	    offset += 3*npts*sizeof(float_sw4);
	    a_AlphaVE[g][m].insert_subarray( mWindow[g][0], mWindow[g][1], mWindow[g][2], mWindow[g][3], mWindow[g][4], mWindow[g][5],
					     doubleField );
	 }
      }
      delete[] doubleField;
   }
   if( iread )
      close(fid);
}

//-----------------------------------------------------------------------
float_sw4 CheckPoint::getDt()
{
   float_sw4 dt;
   if( mEW->getRank() == 0 )
   {
      std::stringstream s;
      if( mRestartPathSet )
	 s << mRestartPath << "/";
      else if( mEW->getPath() != "./" )
	 s << mEW->getPath();
      s << mRestartFile; // string 's' is the file name including path
      int fid = open( const_cast<char*>(s.str().c_str()), O_RDONLY ); 
      CHECK_INPUT(fid != -1, "CheckPoint::getDt: Error opening: " << s.str() );
      lseek(fid,3*sizeof(int)+sizeof(float_sw4),SEEK_SET);
      size_t nr=read(fid,&dt,sizeof(float_sw4));
      CHECK_INPUT( nr == sizeof(float_sw4) ,
	 "CheckPoint::getDt, error reading time step from restart file\n");
      close(fid);
   }
   MPI_Bcast( &dt, 1, mEW->m_mpifloat, 0, MPI_COMM_WORLD );      
   return dt;
}

//-----------------------------------------------------------------------
void CheckPoint::write_header( int& fid, float_sw4 a_time, int a_cycle,
			       int& hsize )
{
   //
   // Header format: prec - precision 8--> double, 4--> single (int)
   //                ng   - Number of grids (int)
   //                time - Time (float)
   //                cycle- Time step number corresponding to Time (int)
   //                dt   - Time step size (float)
   //                nmech- Number of mechanisms in attenuation model
   //                dims(g,1:6) - Size of array on grid g,
   //                       dim(1) <= i <= dim(2), dim(3) <= j <= dim(4)
   //                       dim(5) <= k <= dim(6)
   //   
   int prec = m_double ? 8 : 4;
   size_t ret = write(fid,&prec,sizeof(int));
   CHECK_INPUT( ret == sizeof(int),"CheckPoint::write_header: Error writing precision" );

   int ng = mEW->mNumberOfGrids;
   ret = write(fid,&ng,sizeof(int));
   CHECK_INPUT( ret == sizeof(int),"CheckPoint::write_header: Error writing ng" );

   ret = write(fid,&a_time,sizeof(float_sw4));
   CHECK_INPUT( ret == sizeof(float_sw4),"CheckPoint::write_header: Error writing time" );

   ret = write(fid,&a_cycle,sizeof(int));
   CHECK_INPUT( ret == sizeof(int),"CheckPoint::write_header: Error writing cycle" );

   float_sw4 dt = mEW->getTimeStep();
   ret = write(fid,&dt,sizeof(float_sw4));
   CHECK_INPUT( ret == sizeof(float_sw4),"CheckPoint::write_header: Error writing dt" );

   int nmech = mEW->getNumberOfMechanisms();
   ret = write(fid,&nmech,sizeof(int));
   CHECK_INPUT( ret == sizeof(int),"CheckPoint::write_header: Error writing nmech" );
   for(int g = 0; g < ng ;g++ )
   {
      int globalSize[6];
      globalSize[0] = 1;
      globalSize[1] = mGlobalDims[g][1]-mGlobalDims[g][0]+1;
      globalSize[2] = 1;
      globalSize[3] = mGlobalDims[g][3]-mGlobalDims[g][2]+1;
      globalSize[4] = 1;
      globalSize[5] = mGlobalDims[g][5]-mGlobalDims[g][4]+1;
      ret = write( fid, globalSize, 6*sizeof(int) );
      CHECK_INPUT( ret == 6*sizeof(int),"CheckPoint::write_header: Error writing global sizes" );
      //      cout << "wrote global size " << globalSize[0] << " " << globalSize[1] << " " << globalSize[2] << " " 
      //	   << globalSize[3] << " " << globalSize[4] << " " << globalSize[5] << endl; 
   }
   hsize = (4+6*ng)*sizeof(int) + 2*sizeof(float_sw4);
}

//-----------------------------------------------------------------------
void CheckPoint::read_header( int& fid, float_sw4& a_time, int& a_cycle,
			      int& hsize )
{
   //
   // Header format: prec - precision 8--> double, 4--> single (int)
   //                ng   - Number of grids (int)
   //                time - Time (float)
   //                cycle- Time step number corresponding to Time (int)
   //                dt   - Time step size (float)
   //                nmech- Number of mechanisms in attenuation model
   //                dims(g,1:6) - Size of array on grid g,
   //                       dim(1) <= i <= dim(2), dim(3) <= j <= dim(4)
   //                       dim(5) <= k <= dim(6)
   //   
   int prec;
   size_t ret = read(fid,&prec,sizeof(int));
   CHECK_INPUT( ret == sizeof(int),"CheckPoint::read_header: Error reading precision" );
   CHECK_INPUT( (m_double && prec==8) || (!m_double && prec==4), 
		"CheckPoint::read_header, floating point precision on restart file" <<
		" does not match precision in solver");
   int ng;
   ret = read(fid,&ng,sizeof(int));
   CHECK_INPUT( ret == sizeof(int),"CheckPoint::read_header: Error reading ng" );
   CHECK_INPUT( ng == mEW->mNumberOfGrids, "CheckPoint::read_header: Error number of grids on restart file"
		<< " does not match number of grids in solver" );

   ret = read(fid,&a_time,sizeof(float_sw4));
   CHECK_INPUT( ret == sizeof(float_sw4),"CheckPoint::read_header: Error reading time" );

   ret = read(fid,&a_cycle,sizeof(int));
   CHECK_INPUT( ret == sizeof(int),"CheckPoint::read_header: Error reading cycle" );

   float_sw4 dt;
   ret = read(fid,&dt,sizeof(float_sw4));
   CHECK_INPUT( ret == sizeof(float_sw4),"CheckPoint::read_header: Error reading dt" );

   int nmech;
   ret = read(fid,&nmech,sizeof(int));
   CHECK_INPUT( ret == sizeof(int),"CheckPoint::read_header: Error reading nmech" );
   CHECK_INPUT( nmech == mEW->getNumberOfMechanisms(), "CheckPoint::read_header: Error number "
		<< "of attenuation mechanisms on restart file"
		<< " does not match number of attenuation mechanisms in solver" );

   for(int g = 0; g < ng ;g++ )
   {
      int globalSize[6];
      ret = read( fid, globalSize, 6*sizeof(int) );
      CHECK_INPUT( ret == 6*sizeof(int),"CheckPoint::read_header: Error reading global sizes" );
      CHECK_INPUT( globalSize[0] == 1, "CheckPoint::read_checkpoint: Error in global sizes, " <<
		   "low i-index is "<< globalSize[0]);
      CHECK_INPUT( globalSize[1] == mGlobalDims[g][1]-mGlobalDims[g][0]+1, "CheckPoint::read_checkpoint: Error in global sizes, "
		   << "upper i-index is " << globalSize[1]);
      CHECK_INPUT( globalSize[2] == 1, "CheckPoint::read_checkpoint: Error in global sizes, " <<
		   "low j-index is "<< globalSize[2]);
      CHECK_INPUT( globalSize[3] == mGlobalDims[g][3]-mGlobalDims[g][2]+1, "CheckPoint::read_checkpoint: Error in global sizes, "
		   << "upper j-index is " << globalSize[3]);
      CHECK_INPUT( globalSize[4] == 1, "CheckPoint::read_checkpoint: Error in global sizes, " <<
		   "low k-index is "<< globalSize[4]);
//      CHECK_INPUT( globalSize[5] == mGlobalDims[g][5], "CheckPoint::read_checkpoint: Error in global sizes, "
      CHECK_INPUT( globalSize[5] == mGlobalDims[g][5]-mGlobalDims[g][4]+1,
                   "CheckPoint::read_checkpoint: Error in global sizes, " << "upper k-index is " << globalSize[5]);
   }
   hsize = (4+6*ng)*sizeof(int) + 2*sizeof(float_sw4);
}

//-----------------------------------------------------------------------
void CheckPoint::cycle_checkpoints( string CheckPointFile )
{
   // Keep previous check point, remove the second previous. 
   // m_fileno is initialized to zero in constructor.
   m_fileno++;
   if( m_fileno >= 3 )
   {
      // Delete  mCheckPointFileMM
      if( m_parallel_io[0]->proc_zero() )
      {
	 int er= unlink(const_cast<char*>(mCheckPointFileMM.c_str()));
	 CHECK_INPUT(er==0, "CheckPoint::cycle_checkpoints, error deleting old check point file "
		     << mCheckPointFileMM );
      }
      mCheckPointFileMM = mCheckPointFileM;
      mCheckPointFileM  = CheckPointFile;
   }
   else if( m_fileno == 1 )
      mCheckPointFileMM = CheckPointFile;
   else if( m_fileno == 2 )
      mCheckPointFileM  = CheckPointFile;
}

//-----------------------------------------------------------------------
void CheckPoint::set_restart_file( string fname, size_t bufsize )
{
   mRestartFile = fname;
   m_bufsize    = bufsize;
   mDoRestart = true;

}

//-----------------------------------------------------------------------
void CheckPoint::set_restart_path( string restartPath )
{
   mRestartPath = restartPath;
   mRestartPathSet = true;
}

//-----------------------------------------------------------------------
std::string CheckPoint::get_restart_path()
{
  std::string retval;
  if (mRestartPathSet)
  {
    retval = mRestartPath;
    return retval;
  }
}

//-----------------------------------------------------------------------
void CheckPoint::set_checkpoint_file( string fname, int cycle, int cycleInterval,
				      size_t bufsize )
{
   mCheckPointFile  = fname;
   mWritingCycle    = cycle;
   mCycleInterval   = cycleInterval;
   m_bufsize        = bufsize;
   mDoCheckPointing = true;
}
