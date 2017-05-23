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
CheckPoint::CheckPoint( EW* a_ew,
			float_sw4 time, 
			float_sw4 timeInterval, 
			int cycle, 
			int cycleInterval,
			string fname,
			size_t bufsize ) :
   mEW(a_ew),
   mTime(time),
   m_time_done(false),
   mTimeInterval(timeInterval),
   mWritingCycle(cycle),
   mCycleInterval(cycleInterval),
   mStartTime(0.0),
   mNextTime(0.0),
   mFilePrefix(fname),
   mPreceedZeros(0),
   m_double(true),
   mRestartFile("restart"),
   m_winallocated(false),
   m_bufsize(bufsize)
{

}

//-----------------------------------------------------------------------
CheckPoint::CheckPoint( EW* a_ew, string fname, size_t bufsize ) :
   mEW(a_ew),
   mTime(0.0),
   m_time_done(false),
   mTimeInterval(0.0),
   mWritingCycle(10),
   mCycleInterval(0),
   mStartTime(0.0),
   mNextTime(0.0),
   mFilePrefix("chkpt"),
   mPreceedZeros(0),
   m_double(true),
   mRestartFile(fname),
   m_winallocated(false),
   m_bufsize(bufsize)
{

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
void CheckPoint::setup_sizes( )
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
   for( int g=0 ; g < mEW->mNumberOfGrids ; g++ )
   {
      mWindow[g][0] = mEW->m_iStartInt[g];
      mWindow[g][1] = mEW->m_iEndInt[g];
      mWindow[g][2] = mEW->m_jStartInt[g];
      mWindow[g][3] = mEW->m_jEndInt[g];
      mWindow[g][4] = mEW->m_kStartInt[g];
      mWindow[g][5] = mEW->m_kEndInt[g];


      mGlobalDims[g][0] = 1;
      mGlobalDims[g][1] = mEW->m_global_nx[g];
      mGlobalDims[g][2] = 1;
      mGlobalDims[g][3] = mEW->m_global_ny[g];
      mGlobalDims[g][4] = mWindow[g][4];
      mGlobalDims[g][5] = mWindow[g][5];

 // The 3D array is assumed to span the entire computational domain
      m_ihavearray.resize( mEW->mNumberOfGrids );
      m_ihavearray[g] = true;
   }
   setSteps( mEW->mNumberOfTimeSteps );
   define_pio();
   //   cout << "mwind = " << mWindow[0][4] << " " << mWindow[0][5] << endl;
   //   cout << "globaldims = " << mGlobalDims[0][4] << " " << mGlobalDims[0][5] << endl;
}

//-----------------------------------------------------------------------
void CheckPoint::define_pio( )
{
   int glow = 0, ghigh = mEW->mNumberOfGrids;
   m_parallel_io = new Parallel_IO*[ghigh-glow+1];
   for( int g=glow ; g < ghigh ; g++ )
   {
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
      m_parallel_io[g-glow] = new Parallel_IO( iwrite, mEW->usingParallelFS(), global, local, start, m_bufsize );
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
   // -----------------------------------------------
   // Check based on cycle
   // -----------------------------------------------
   //   cout << "in time to write " << mWritingCycle << " " << mCycleInterval << " " << " " << mTime << " " <<  mTimeInterval << " " << endl;
   bool do_it=false;
   if( cycle == mWritingCycle )
      do_it = true;
   if( mCycleInterval !=  0 && cycle%mCycleInterval == 0 && time >= mStartTime) 
      do_it = true;

   // ---------------------------------------------------
   // Check based on time
   // ---------------------------------------------------
   if(mTime > 0.0 && (  mTime <= time + dt*0.5 ) && !m_time_done )
   {
      m_time_done = true;
      do_it = true;
   }
   if( mTimeInterval != 0.0 && mNextTime <= time + dt*0.5 )
   {
      mNextTime += mTimeInterval;
      if (time >= mStartTime) do_it =  true;
   }
   return do_it;
}

//-----------------------------------------------------------------------
void CheckPoint::compute_file_suffix( int cycle, std::stringstream& fileSuffix )
{
   fileSuffix << mFilePrefix << ".cycle=";
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
				   vector<Sarray>& a_U )
{
  //File format: 
  //
   //   [precision(int), npatches(int), time(double), sizes_1(6 int), .. sizes_ng(6 int),
   //    Um(3 component float/double array), U(3 component float/double array)]
  // Would it be possible to save the entire input file to the restart file ? 
  //
   string path = mEW->mPath;
   int ng = mEW->mNumberOfGrids;
  // offset initialized to header size:
   off_t offset = 3*sizeof(int) + 1*sizeof(float_sw4) + ng*(6*sizeof(int));

   bool iwrite = false;
   for( int g=0 ; g < ng ; g++ )
      iwrite = iwrite || m_parallel_io[g]->i_write();
  
   int fid=-1;
   std::stringstream s, fileSuffix;

   if( iwrite )
   {
      compute_file_suffix( a_cycle, fileSuffix );
      if( path != "." )
	 s << path;
      s << fileSuffix.str(); // string 's' is the file name including path
   }
   // Open file from processor zero and write header.
   if( m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
      CHECK_INPUT(fid != -1, "CheckPoint::write_file: Error opening: " << s );
      int myid;

      MPI_Comm_rank( MPI_COMM_WORLD, &myid );
      std::cout << "writing check point on file " << s.str() << " using " <<
	 m_parallel_io[0]->n_writers() << " writers" << std::endl;
      //      std::cout << " (msg from proc # " << myid << ")" << std::endl;

      int prec = m_double ? 8 : 4;
      size_t ret = write(fid,&prec,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"CheckPoint::write_checkpoint: Error writing precision" );
      ret = write(fid,&ng,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"CheckPoint::write_checkpoint: Error writing ng" );
      ret = write(fid,&a_time,sizeof(float_sw4));
      CHECK_INPUT( ret == sizeof(float_sw4),"CheckPoint::write_checkpoint: Error writing time" );
      ret = write(fid,&a_cycle,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"CheckPoint::write_checkpoint: Error writing cycle" );
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
	 CHECK_INPUT( ret == 6*sizeof(int),"CheckPoint::write_checkpoint: Error writing global sizes" );
	 cout << "wrote global size " << globalSize[0] << " " << globalSize[1] << " " << globalSize[2] << " " 
	      << globalSize[3] << " " << globalSize[4] << " " << globalSize[5] << endl; 
      }
      fsync(fid);
   }
   m_parallel_io[0]->writer_barrier();

   // Open file from all writers
   if( iwrite && !m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(s.str().c_str()), O_WRONLY );
      CHECK_INPUT(fid != -1, "CheckPoint::write_checkpoint:: Error opening: " << s );
   }

   // Write data blocks
   for( int g = 0 ; g < ng ; g++ )
   {
      size_t npts = ((size_t)(mGlobalDims[g][1]-mGlobalDims[g][0]+1))*
 	            ((size_t)(mGlobalDims[g][3]-mGlobalDims[g][2]+1))*
	            ((size_t)(mGlobalDims[g][5]-mGlobalDims[g][4]+1));

      if( !mEW->usingParallelFS() || g == 0 )
	 m_parallel_io[g]->writer_barrier();
      
      size_t nptsloc = ((size_t)(mEW->m_iEndInt[g]-mEW->m_iStartInt[g]+1))*
 	               ((size_t)(mEW->m_jEndInt[g]-mEW->m_jStartInt[g]+1))*
	               ((size_t)(mEW->m_kEndInt[g]-mEW->m_kStartInt[g]+1));
      
      // Write without ghost points. Would probably work with ghost points too.
      float_sw4* doubleField = new float_sw4[3*nptsloc];
      if( m_double )
      {
	char cprec[]="double";
	a_Um[g].extract_subarray( mEW->m_iStartInt[g], mEW->m_iEndInt[g], mEW->m_jStartInt[g],
				  mEW->m_jEndInt[g], mEW->m_kStartInt[g], mEW->m_kEndInt[g],
				  doubleField );
	m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	offset += 3*npts*sizeof(float_sw4);

	a_U[g].extract_subarray( mEW->m_iStartInt[g], mEW->m_iEndInt[g], mEW->m_jStartInt[g],
				 mEW->m_jEndInt[g], mEW->m_kStartInt[g], mEW->m_kEndInt[g],
				 doubleField );
	m_parallel_io[g]->write_array( &fid, 3, doubleField, offset, cprec );
	offset += 3*npts*sizeof(float_sw4);
      }
      else
	 cout <<"CheckPoint::write_checkpoint, currently only implemented double precision" << endl;
      //      else
      //      {
      //	char cprec[]="float";
      //	m_parallel_io[g]->write_array( &fid, 3, a_Um[g].c_ptr(), offset, cprec );
      //	offset += 3*npts*sizeof(float);
      //	m_parallel_io[g]->write_array( &fid, 3, m_U[g].c_ptr(), offset, cprec );
      //	offset += 3*npts*sizeof(float);
      //      }
      delete[] doubleField;
   }
   if( iwrite )
      close(fid);
}

//-----------------------------------------------------------------------
void CheckPoint::read_checkpoint( float_sw4& a_time, int& a_cycle,
				  vector<Sarray>& a_Um, vector<Sarray>& a_U )
{
   // It is assumed that the arrays are already declared with the right 
   // dimensions. This routine will check that sizes match, but will not
   // allocate or resize the arrays Um and U.


   string path=mEW->mPath;

   //   vector<Parallel_IO*> parallel_io;
   //   define_parallel_io( parallel_io );
   int ng = mEW->mNumberOfGrids;

  // offset initialized to header size:
   off_t offset = 3*sizeof(int) + 1*sizeof(float_sw4) + ng*(6*sizeof(int));

   bool iread = false;
   for( int g=0 ; g < ng ; g++ )
      iread = iread || m_parallel_io[g]->i_write();
  
   int fid=-1;
   std::stringstream s;

   //   if( iread )
   //   {
   //      if( path != "." )
   //	 s << path << "/" << a_filename;
   //      else
   //	 s << a_filename; 
   //   }
   // Open file from processor zero and read header.

   int prec=0;
   a_cycle = 0;
   a_time  = 0;
   if( m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(mRestartFile.c_str()), O_RDONLY ); 
      CHECK_INPUT(fid != -1, "CheckPoint::read_checkpoint: Error opening: " << s.str() );
      int myid;

      MPI_Comm_rank( MPI_COMM_WORLD, &myid );
      std::cout << "reading check point on file " << s.str();
      std::cout << " (msg from proc # " << myid << ")" << std::endl;

      size_t ret = read(fid,&prec,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"CheckPoint::read_checkpoint: Error reading precision" );
      int ngfile;
      ret = read(fid,&ngfile,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"CheckPoint::read_checkpoint: Error reading ng" );
      CHECK_INPUT( ng == ngfile , "CheckPoint::read_checkpoint: Error number of grids on file " <<
		   "does not match number of grids in computation" );

      ret = read(fid,&a_time,sizeof(float_sw4));
      CHECK_INPUT( ret == sizeof(float_sw4),"CheckPoint::read_checkpoint: Error reading time" );
      ret = read(fid,&a_cycle,sizeof(int));
      CHECK_INPUT( ret == sizeof(int),"CheckPoint::read_checkpoint: Error reading cycle" );

      for(int g = 0; g < ng ;g++ )
      {
	 int globalSize[6];
	 ret = read( fid, globalSize, 6*sizeof(int) );
	 CHECK_INPUT( ret == 6*sizeof(int),"CheckPoint::read_checkpoint: Error reading global sizes" );
	 CHECK_INPUT( globalSize[0] == 1, "CheckPoint::read_checkpoint: Error in global sizes, " <<
		      "low i-index is "<< globalSize[0]);
	 CHECK_INPUT( globalSize[1] == mGlobalDims[g][1], "CheckPoint::read_checkpoint: Error in global sizes, "
		      << "upper i-index is " << globalSize[1]);
	 CHECK_INPUT( globalSize[2] == 1, "CheckPoint::read_checkpoint: Error in global sizes, " <<
		      "low j-index is "<< globalSize[2]);
	 CHECK_INPUT( globalSize[3] == mGlobalDims[g][3], "CheckPoint::read_checkpoint: Error in global sizes, "
		      << "upper j-index is " << globalSize[3]);
	 CHECK_INPUT( globalSize[4] == 1, "CheckPoint::read_checkpoint: Error in global sizes, " <<
		      "low k-index is "<< globalSize[4]);
	 CHECK_INPUT( globalSize[5] == mGlobalDims[g][5], "CheckPoint::read_checkpoint: Error in global sizes, "
		      << "upper k-index is " << globalSize[5]);
      }
   }

   m_parallel_io[0]->writer_barrier();
   int tmpprec=prec;
   MPI_Allreduce( &tmpprec, &prec, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
   int tmpcyc=a_cycle;
   MPI_Allreduce( &tmpcyc, &a_cycle, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
   float_sw4 tmpt  = a_time;
   MPI_Allreduce( &tmpt, &a_time, 1, mEW->m_mpifloat, MPI_MAX, MPI_COMM_WORLD );

   // Open file from all readers
   if( iread && !m_parallel_io[0]->proc_zero() )
   {
      fid = open( const_cast<char*>(mRestartFile.c_str()), O_RDONLY );
      CHECK_INPUT(fid != -1, "CheckPoint::read_checkpoint:: Error opening file: " << s.str() );
   }

   // Read data blocks
   for( int g = 0 ; g < ng ; g++ )
   {
      size_t npts = ((size_t)(mGlobalDims[g][1]-mGlobalDims[g][0]+1))*
	            ((size_t)(mGlobalDims[g][3]-mGlobalDims[g][2]+1))*
	            ((size_t)(mGlobalDims[g][5]-mGlobalDims[g][4]+1));
      //      size_t npts = ((size_t)(mEW->m_global_nx[g]))*
      //                    ((size_t)(mEW->m_global_ny[g]))*
      // 	            ((size_t)(mEW->m_global_nz[g]));
      size_t nptsloc = ((size_t)(mEW->m_iEndInt[g]-mEW->m_iStartInt[g]+1))*
 	               ((size_t)(mEW->m_jEndInt[g]-mEW->m_jStartInt[g]+1))*
	               ((size_t)(mEW->m_kEndInt[g]-mEW->m_kStartInt[g]+1));
      
      if( !mEW->usingParallelFS() || g == 0 )
	 m_parallel_io[g]->writer_barrier();
      
      // array without ghost points read into doubleField
      float_sw4* doubleField = new float_sw4[3*nptsloc];
      if( prec == 8 )
      {
	 m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, "double" );
	 offset += 3*npts*sizeof(float_sw4);
	 a_Um[g].insert_subarray( mEW->m_iStartInt[g], mEW->m_iEndInt[g], mEW->m_jStartInt[g], mEW->m_jEndInt[g],
			       mEW->m_kStartInt[g], mEW->m_kEndInt[g], doubleField );
	 m_parallel_io[g]->read_array( &fid, 3, doubleField, offset, "double" );
	 a_U[g].insert_subarray( mEW->m_iStartInt[g], mEW->m_iEndInt[g], mEW->m_jStartInt[g], mEW->m_jEndInt[g],
			       mEW->m_kStartInt[g], mEW->m_kEndInt[g], doubleField );
	 offset += 3*npts*sizeof(float_sw4);
      }
      //      else
      //      {
      //	 parallel_io[g]->read_array( &fid, 1, doubleField, offset, "float" );
      //	 offset += npts*sizeof(float);
      //      }
      delete[] doubleField;
   }
   if( iread )
      close(fid);
   //   for( int g=0 ; g < mNumberOfGrids ; g++ )
   //      delete parallel_io[g];
}
