#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>
#include <string>

#include "Parallel_IO-sw4.h"

using namespace std;

#ifdef IOTEST_HDF5
#include "hdf5.h"
//#include "H5Cpp.h"
//using namespace H5;
#endif

void getseed( int& ss )
{
   int fd = open("/dev/urandom",O_RDONLY);
   int mask = 0x0001FFFF;
   int nr;
   nr = read(fd,&ss,sizeof(int));
   //   *ss &= mask;
   close(fd);
}

void initrandom( int seed )
{
  time_t secs = time(0);
  srand48(seed);
  //  srand(secs);
}

double randomnr()
{
  return drand48();
 //  return rand();
}

//-----------------------------------------------------------------------
void decomp1d( int nglobal, int myid, int nproc, int padding, int& s, int& e )
//
// Decompose index space 1 <= i <= nglobal into nproc blocks
// returns start and end indices for block nr. myid, 
//          where 0 <= myid <= nproc-1
//
{
   int olap    = 2*padding;
   int nlocal  = (nglobal + (nproc-1)*olap ) / nproc;
   int deficit = (nglobal + (nproc-1)*olap ) % nproc;
   if( myid < deficit )
      s = myid*(nlocal-olap) + myid+1;
   else
      s = myid*(nlocal-olap) + deficit+1;

   if (myid < deficit)
      nlocal = nlocal + 1;

   e = s + nlocal - 1;
}

//-----------------------------------------------------------------------
bool proc_decompose_2d( int ni, int nj, int nproc, int padding, int proc_max[2] )
{
   // This routine determines a decomposition of nproc processors into
   // a 2D processor array  proc_max[0] x proc_max[1], which gives minimal 
   // communication boundary for a grid with ni x nj points.

   double fmin = ni+nj;
   bool first  = true;
   int p1max = ni, p2max=nj;
   if( padding > 1 )
   {
      p1max   = ni/padding;
      p2max   = nj/padding;
      p1max = p1max < 1 ? 1 : p1max;
      p2max = p2max < 1 ? 1 : p2max;
   }
   for( int p1 = 1 ; p1 <= nproc; p1++)
      if( nproc%p1 == 0 )
      {
        int p2 = nproc/p1;
        if( p1 <= p1max && p2 <= p2max )
        {
           // int w1 = p1==1?0:1;
           // int w2 = p2==1?0:1;
           // double f = w2*(double)(ni)/p1 + w1*(double)(nj)/p2;
// try to make each subdomain as square as possible
          double f = fabs((double)(ni)/p1 - (double)(nj)/p2);
           if( f < fmin || first )
           {
              fmin = f;
              proc_max[0]   = p1;
              proc_max[1]   = p2;
              first= false;
           }
        }
      }
   return !first;
}

//-----------------------------------------------------------------------
bool proc_decompose_3d( int ni, int nj, int nk, int nproc, int padding,
			int proc_max[3] )
{
   // This routine determines a decomposition of nproc processors into
   // a 3D processor array  proc_max[0] x proc_max[1] x proc_max[2], which gives minimal 
   // communication boundary for a grid with ni x nj x nk points.

   double fmin = ni+nj+nk;
   bool first  = true;
   int p1max = ni, p2max=nj, p3max=nk;
   if( padding > 1 )
   {
      p1max = ni/padding;
      p2max = nj/padding;
      p3max = nk/padding;
      p1max = p1max < 1 ? 1 : p1max;
      p2max = p2max < 1 ? 1 : p2max;
      p3max = p3max < 1 ? 1 : p3max;
   }
   for( int p3 = 1 ; p3 <= nproc; p3++)
      for( int p1 = 1 ; p1 <= nproc; p1++)
	 if( nproc%(p1*p3) == 0 )
	 {
	    int p2 = nproc/(p1*p3);
	    if( p1 <= p1max && p2 <= p2max && p3 <= p3max )
	    {
           // int w1 = p1==1?0:1;
           // int w2 = p2==1?0:1;
           // double f = w2*(double)(ni)/p1 + w1*(double)(nj)/p2;
// try to make each subdomain as square as possible
	       double f = fabs((double)(ni)/p1 - (double)(nj)/p2) + 
		  fabs((double)(nj)/p2 - (double)(nk)/p3);
	       if( f < fmin || first )
	       {
		  fmin = f;
		  proc_max[0]   = p1;
		  proc_max[1]   = p2;
		  proc_max[2]   = p3;
		  first= false;
           }
        }
      }
   return !first;
}

//-----------------------------------------------------------------------
void myproc_3d( int myname, int pmax[3], int proc[3] )
{
   proc[0] = myname % pmax[0] + 1;
   int remain  = ( myname - proc[0] + 1 )/pmax[0];
   proc[1] = remain % pmax[1] + 1;
   proc[2] = ( remain - proc[1] + 1 )/pmax[1] + 1;
}

//-----------------------------------------------------------------------
size_t read_dble_wlim( int* fid, double* rbuf, size_t nelem, size_t limit )
{
// Read a vector of `nelem' elements of type double into rbuf, with maximum of 
// elements per read is limited to `limit'.
//
// Input: fid   - File descriptor previously opened with `open'.
//        rbuf  - Pointer to vector of doubles.
//        nelem - Number of elements in rbuf.
//        limit - Maximum number of elements to read at each call to `read'.
// Output: Returns the number of bytes read.
//
   size_t ind=0;
   size_t nreads = (nelem % limit) == 0 ? nelem/limit:nelem/limit+1;
   for( int i= 0 ; i < nreads ; i++ )
   {
      size_t nrtoread = nelem-ind > limit ? limit:nelem-ind;
      size_t nr = read( *fid, &rbuf[ind], nrtoread*sizeof(double) );
      if( nr != nrtoread*sizeof(double) )
      {
	 cout << "ERROR in read_dble_wlim, read " << nr <<
	    " bytes, requested " << nrtoread*sizeof(double) << " bytes " << endl;
	 return ind*sizeof(double)+nr;
      }
      ind += nrtoread;
   }
   return ind*sizeof(double);
}

//-----------------------------------------------------------------------
template<class T> size_t read_with_limit( int* fid, T* rbuf, size_t nelem, size_t limit )
{
// Read a vector of `nelem' elements of type T into rbuf, with maximum of 
// elements per read is limited to `limit'.
//
// Input: fid   - File descriptor previously opened with `open'.
//        rbuf  - Pointer to vector of doubles.
//        nelem - Number of elements in rbuf.
//        limit - Maximum number of elements to read at each call to `read'.
// Output: Returns the number of bytes read.
//
   size_t ind=0;
   size_t nreads = (nelem % limit) == 0 ? nelem/limit:nelem/limit+1;
   for( int i= 0 ; i < nreads ; i++ )
   {
      size_t nrtoread = nelem-ind > limit ? limit:nelem-ind;
      size_t nr = read( *fid, &rbuf[ind], nrtoread*sizeof(T) );
      if( nr != nrtoread*sizeof(T) )
      {
	 cout << "ERROR in read_with_limit, read " << nr <<
	    " bytes, requested " << nrtoread*sizeof(T) << " bytes " << endl;
	 return ind*sizeof(T)+nr;
      }
      ind += nrtoread;
   }
   return ind*sizeof(T);
}


//-----------------------------------------------------------------------
void get_size( const char* fname, int& ni, int& nj, int& nk )
{
   int fd = open( fname, O_RDONLY );
   if( fd != -1 )
   {
      size_t nr = read(fd,&ni,sizeof(int));
      if( nr != sizeof(int) )
	 cout << "Error reading ni, nr=" << nr << endl;
      nr = read(fd,&nj,sizeof(int));
      if( nr != sizeof(int) )
	 cout << "Error reading nj, nr=" << nr << endl;
      nr = read(fd,&nk,sizeof(int));
      if( nr != sizeof(int) )
	 cout << "Error reading nk, nr=" << nr << endl;
      close(fd);
   }
   else
   {
      cout << "ERROR opening " << fname << " for reading"  << endl;
      MPI_Abort(MPI_COMM_WORLD,0);
   }
}

//-----------------------------------------------------------------------
bool check_result( int locdim[3], double* w, int start[3], double hx,
		   double hy, double hz, double& err, bool corder )
{
   err = 0;
   double pi= 4*atan(1.0);
   if( corder )
   {
      double* cpizp = new double[locdim[2]];
      for( int k=start[2] ; k <= locdim[2]+start[2]-1 ; k++ )
	 cpizp[k-start[2]] = cos(pi*(k-1)*hz);
      double* eyp = new double[locdim[1]];
      for( int j=start[1] ; j <= locdim[1]+start[1]-1 ; j++ )
	 eyp[j-start[1]]  = exp(-(j-1)*hy);

      for( int i=start[0] ; i <= locdim[0]+start[0]-1 ; i++ )
	 for( int j=start[1] ; j <= locdim[1]+start[1]-1 ; j++ )
	    for( int k=start[2] ; k <= locdim[2]+start[2]-1 ; k++ )
	    {
	       size_t ind = k-start[2]+locdim[2]*( j-start[1]+locdim[1]*(i-start[0]));
	       double locerr = fabs(w[ind]-(i-1)*hx*eyp[j-start[1]]*cpizp[k-start[2]]);
	       if( locerr > err )
		  err = locerr;
	    }
      delete[] eyp;
      delete[] cpizp;
   }
   else
   {
      for( int k=start[2] ; k <= locdim[2]+start[2]-1 ; k++ )
      {
	 double z = (k-1)*hz;
	 double cpiz = cos(pi*z);
	 for( int j=start[1] ; j <= locdim[1]+start[1]-1 ; j++ )
	 {
	    double y = (j-1)*hy;
	    double ey = exp(-y);
	    for( int i=start[0] ; i <= locdim[0]+start[0]-1 ; i++ )
	    {
	       double x = (i-1)*hx;
	       size_t ind = i-(start[0]) + locdim[0]*((j-(start[1])) + locdim[1]*(k-(start[2])));
	       double locerr = fabs(w[ind]-x*ey*cpiz);
	       if( locerr > err )
		  err = locerr;
	    }
	 }
      }
   }
   return err < 1e-10;
}

//-----------------------------------------------------------------------
int size_chk( size_t wanted, size_t actual, int ref )
{
   int retval = 0;
   if( wanted != actual )
   {
      cout << "read_wp_patch_f: " << ref << " , error read " << actual <<
	 " bytes, but wanted " << wanted << " bytes " << endl;
      MPI_Abort(MPI_COMM_WORLD,0);
      retval = -1;
      
   }
   return retval;
}

//-----------------------------------------------------------------------
#ifdef IOTEST_HDF5
void get_size_hdf5( std::string fname, int& ni, int& nj, int& nk, int myid )
{
   hid_t fapl_id = H5Pcreate( H5P_FILE_ACCESS );
   herr_t status = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL );
   hid_t file_id = H5Fopen( fname.c_str(), H5F_ACC_RDWR, fapl_id );
   if( file_id < 0 )
   {
      cout << "ERROR opening " << fname << " for reading"  << endl;
      MPI_Abort(MPI_COMM_WORLD,0);
   }
   hid_t dataset_id = H5Dopen( file_id, "testdata", H5P_DEFAULT );
   hid_t dataspace_id = H5Dget_space( dataset_id );
   hsize_t dims[3], maxdims[3];
   H5Sget_simple_extent_dims( dataspace_id, dims, maxdims );
   ni = dims[0];
   nj = dims[1];
   nk = dims[2];
   //   if( myid == 0 )
   //   {
   //      cout << "dims    = " << dims[0] << " " << dims[1] << " " << dims[2] << endl;
   //      cout << "maxdims = " << maxdims[0] << " " << maxdims[1] << " " << maxdims[2] << endl;
   //   }
   H5Sclose(dataspace_id);
   H5Dclose(dataset_id);
   H5Fclose(file_id);
   H5Pclose(fapl_id);
}

//-----------------------------------------------------------------------
int read_hdf5_patch_f( std::string fname, int nc, int ni, int nj, int nk, int nig,
                     int njg, int nkg, int ds, double* w, int* lgmap,
                     size_t pos0 )
{
   int myid;
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   hid_t fapl_id = H5Pcreate( H5P_FILE_ACCESS );
   herr_t status = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL );
   if( status == -1 )
      cout << "Error from H5Pset_fapl_mpio " << endl;
   hid_t file_id      = H5Fopen( fname.c_str(), H5F_ACC_RDWR, fapl_id );
   hid_t dataset_id   = H5Dopen( file_id, "testdata", H5P_DEFAULT );
   hid_t dataspace_id = H5Dget_space( dataset_id );
   hid_t plist_id     = H5Pcreate( H5P_DATASET_XFER );
   status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
   if( status == -1 )
      cout << "Error from H5Pset_dxpl_mpio " << endl;

   //   hid_t plist_id = H5Dget_access_plist( dataset_id );
   if( ni==nig && nj==njg && nk==nkg )
   {
      // Read full array into this processor
      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace_id,
		       plist_id, w );
      if( status == -1 )
	 cout << "Error from H5Dread, full array read " << endl;
   }
   else
   {
      hsize_t dims[3]={ni,nj,nk};
      hid_t subcube_id = H5Screate_simple( 3, dims, NULL );
      hsize_t stride[3]={1,1,1};
      hsize_t block[3] ={1,1,1};
      hsize_t start[3]={lgmap[0],lgmap[1],lgmap[2]};
      //      hsize_t dimsg[3]={nig,njg,nkg};
      //      hsize_t dataspace2_id = H5Screate_simple( 3, dimsg, NULL );
      //      status = H5Sselect_hyperslab(dataspace2_id,H5S_SELECT_SET, start, stride, dims, block); 
      //      cout << myid << " start " << start[0] << " " << start[1] << " " << start[2] << endl;
      //      cout << myid << " dims "  << dims[0] << " " << dims[1] << " " << dims[2] << endl;
      status = H5Sselect_hyperslab(dataspace_id,H5S_SELECT_SET, start, stride, dims, block);
      if( status == -1 )
	 cout << "Error from H5Sselect_hyperslab " << endl;

      if( ds == 8 )
      {
	 status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, subcube_id, dataspace_id,
				 plist_id, w );
	 if( status == -1 )
	    cout << "Error from H5Dread, part of array access " << endl;
      }
      else
      {
	 //	 float* ws = new float[ni*nj*nk];
	 //	 dataset.read(ws,PredType::NATIVE_FLOAT,subcube,dataspace);
	 //	 for( size_t i=0 ; ((size_t)ni)*nj*nk ;i++ )
	 //	    w[i] = ws[i];
	 //	 delete[] ws;
	 cout << "Float data type NYI" << endl;
      }
      H5Sclose(subcube_id);
      //      H5Sclose(dataspace2_id);
   }
   H5Pclose(plist_id);
   H5Sclose(dataspace_id);
   H5Dclose(dataset_id);
   H5Fclose(file_id);
   H5Pclose(fapl_id);
}
#else
// Dummy routines, will never be called
void get_size_hdf5( std::string fname, int& ni, int& nj, int& nk, int myid ){}
int read_hdf5_patch_f( std::string fname, int nc, int ni, int nj, int nk, int nig,
                     int njg, int nkg, int ds, double* w, int* lgmap,
		       size_t pos0 ){}
#endif
//-----------------------------------------------------------------------
int read_wp_patch_f( int* fd, int nc, int ni, int nj, int nk, int nig,
                     int njg, int nkg, int ds, double* w, int* lgmap,
                     size_t pos0 )
{
   // Old IO routine, use for comparison
   int datasize, ig, jg, kg, i, j, k;
   size_t pos, nr;
   float* ws;

   if( ds == 8 )
      datasize = sizeof(double);
   else if( ds == 4 )
      datasize = sizeof(float);
   else
      return -1;
   
   size_t nijg = static_cast<size_t>(nig)*njg;
   size_t nij  = static_cast<size_t>(ni)*nj;

   size_t linux_lim = 268434944; // Maximum number of doubles in one read

/* Full 2D slice of array goes into one processor, read whole block with one read */
   if( nig == ni && njg == nj )
   {
      size_t size=static_cast<size_t>(ni)*nj*nk*nc;
      if( ds == 4 )
	 ws = new float[size];

      kg  = 1+lgmap[2];
      jg  = 1+lgmap[1];
      ig  = 1+lgmap[0];

      pos = pos0+datasize*(nc*(ig-1) + nc*nig*(jg-1) + nc*nijg*(kg-1));
      nr  = lseek( *fd, pos, SEEK_SET );
      if( ds == 4 )
      {
	 if( size > 2*linux_lim )
	    nr = read_with_limit( fd, ws, size, 2*linux_lim );
	 else
	    nr  = read( *fd, ws, size*datasize );
	 size_chk(size*datasize,nr,1);
          /*          memcpy( w, ws, nc*(ni)*(nj)*(nk)*datasize ); */
	 for( size_t i = 1 ; i <= size ; i++ )
	    *(w+i-1) = ws[i-1]; 
      }
      else
      {
	 //	 linux_lim = 1000000;
	 //	 cout << "Resetting linux limit " << endl;
	 if( size > linux_lim )
	    //	    nr = read_dble_wlim( fd, w, size, linux_lim );
	    nr = read_with_limit( fd, w, size, linux_lim );
	 else
	    nr  = read( *fd, w, size*datasize );
	 size_chk(size*datasize,nr,2);
      }
   }
/* Last dimension goes into one processor, read 2D slices */
   else if( nig == ni )
   {
      if( ds == 4 )
	 ws = new float[nc*nij];
      for( k = 1 ; k <= nk ; k++ )
      {
	 kg  = k+lgmap[2];
	 jg  = 1+lgmap[1];
	 ig  = 1+lgmap[0];
	 pos = pos0+datasize*(nc*(ig-1)+static_cast<size_t>(nc)*nig*(jg-1)+nc*nijg*(kg-1));
	 nr  = lseek( *fd, pos, SEEK_SET );
	 if( ds == 4 )
	 {
	    if( nc*nij > 2*linux_lim )
	       nr = read_with_limit( fd, ws+nc*nij*(k-1), nc*nij, 2*linux_lim );
	    else
	       nr  = read( *fd, ws, nc*nij*sizeof(float) );
	    size_chk(nc*nij*sizeof(float),nr,3);
	    for( i=1 ; i <= nc*nij ; i++ )
               *(w+i-1+nc*nij*(k-1) ) = ws[i-1];
	 }
	 else
	 {
	    if( nc*nij > linux_lim )
	       nr = read_dble_wlim( fd, w+nc*nij*(k-1), nc*nij, linux_lim );
	    else
	       nr  = read( *fd, w+nc*nij*(k-1),nc*nij*datasize );
	    size_chk(nc*nij*datasize,nr,4);
	 }
      }
   }
/* General block distribution, read line by line */
   else    
   {
      if( ds == 4 )
	 ws = new float[nc*ni];
       
      for( k= 1 ; k <= nk ; k++ )
	 for( j=1 ; j<= nj ; j++ )
	 {
	    kg  = k+lgmap[2];
	    jg  = j+lgmap[1];
	    ig  = 1+lgmap[0];
	    pos = pos0 + datasize*( nc*(ig-1) + static_cast<size_t>(nc)*nig*(jg-1) +
                                     nc*nijg*(kg-1) );
	    nr  = lseek( *fd, pos, SEEK_SET );
	    if( ds == 4 )
	    {
	       nr  = read( *fd, ws, nc*(ni)*sizeof(float) );
	       size_chk(nc*ni*sizeof(float),nr,5);
	       for( i=1 ; i<=(ni)*nc ; i++ )
		  *(w+nc*nij*(k-1)+nc*(ni)*(j-1)+i-1) = ws[i-1];
	    }
	    else
	    {
	       nr = read( *fd, w+nc*ni*(j-1)+nc*nij*(k-1),
			   nc*(ni)*datasize );
	       size_chk(nc*ni*datasize,nr,6);
	    }
	 }
   }
   if( ds == 4 )
      delete[] ws;

   return 0;
}

void print_help()
{
   cout << "Usage: " << endl;
   cout << "iotest -file filename -wproc wfreq -old -bufsize npts -procdist 3d -procdim 2 2 4 -printdim -save -newfile" << endl;
   cout << endl;
   cout << "       -file filename --> Name of file to read, use extension *.hdf5 for hdf5 format (default testfile.bin)"  << endl;
   cout << "       -old           --> Use old IO routines, not MPI-parallel" << endl;
   cout << "       -wproc wfreq   --> Write from every wfreq:th processor only (SW4 IO only)" << endl;
   cout << "       -bufsize npts  --> Size of internal buffer in number of elements (SW4 IO only)"  << endl;
   cout << "       -procdist [3d|2d]   --> 3d means distribute 3d array on 3d processor grid. \n" <<
           "                               2d will use 2d processor grid to distribute the first two dimensions of the array"  << endl;
   cout << "       -procdim 2 2 4 --> Use processor grid of size indicated. The three numbers must have product\n" <<
           "                          equal to number of processors used. -procdim overrides -procdist" << endl;
   cout << "       -save          --> Save timing numbers on file timings.dat" << endl;
   cout << "       -newfile       --> The -save option will append to timings.dat. By giving -newfile iotest will truncate\n" <<
           "                          the file, starting new recording at the first line" << endl;
   cout << "       -printdim      --> Print out information about array dimensions from each processor, for debug purpose" << endl;
}


//-----------------------------------------------------------------------
int main( int argc, char** argv )
{

   string fname = "testfile.bin";


// Read input parameters
   int seed = 2343234;
   bool seedread = false;

   int ionr = 0;
   int wfreq = 1;
   int arraydist = 1;
   int procdist = 2;
   int p1, p2, p3;
   int olap = 0;
   int nptsbuf = 8000000;
   double efraction = 0.7;
   double relsize = 0.2;
   bool empty = false;
   bool printdim = false;
   bool savedata = false;
   bool appendtofile = true;

   int a=1;
   while( a < argc )
   {
      if( strcmp(argv[a],"-file") == 0 )
      {
	 fname = argv[a+1];
	 a += 2;
      }
      else if( strcmp(argv[a],"-old")== 0 )
      {
         ionr = 1;
	 a++;
      }
      else if( strcmp(argv[a],"-wproc")== 0 )
      {
         wfreq = atoi(argv[a+1]);
	 a += 2;
      }
      else if( strcmp(argv[a],"-olap")== 0 )
      {
         olap = atoi(argv[a+1]);
	 a += 2;
      }
      else if( strcmp(argv[a],"-bufsize")== 0 )
      {
         nptsbuf = atoi(argv[a+1]);
	 a += 2;
      }
      else if( strcmp(argv[a],"-seed")== 0 )
      {
         seed = atoi(argv[a+1]);
	 a += 2;
         seedread = true;
      }
      else if( strcmp(argv[a],"-procdist")== 0 )
      {
	 if( procdist != 3 )
	 {
	    if( strcmp(argv[a+1],"2d")==0 )
	       procdist = 1;
	    else if( strcmp(argv[a+1],"3d")==0 )
	       procdist = 2;
	    else
	       cout << "procdist value " << argv[a+1] << " not recognized" << endl;
	 }
	 else
	 {
	    cout << "procdim given, ignoring procdist option" << endl;
	 }
	 a += 2;	    
      }
      else if( strcmp(argv[a],"-procdim")== 0 )
      {
	 procdist = 3;
	 p1 = atoi(argv[a+1]);
	 p2 = atoi(argv[a+2]);
	 p3 = atoi(argv[a+3]);
	 a += 4;
      }
      else if( strcmp(argv[a],"-arraydist")== 0 )
      {
         if( strcmp(argv[a+1],"regular")==0 )
	    arraydist = 1;
	 else if( strcmp(argv[a+1],"randompert")==0 )
	    arraydist = 2;
	 else if( strcmp(argv[a+1],"random")==0 )
	    arraydist = 3;
	 else
	    cout << "arraydist value " << argv[a+1] << " not recognized" << endl;
	 a += 2;
      }
      else if( strcmp(argv[a],"-printdim")== 0 )
      {
	 printdim = true;
	 a++;
      }
      else if( strcmp(argv[a],"-save")== 0 )
      {
	 savedata = true;
	 a++;
      }
      else if( strcmp(argv[a],"-newfile")== 0 )
      {
	 savedata = true;
	 appendtofile = false;
	 a++;
      }
      else if( strcmp(argv[a],"-h")== 0 )
      {
	 print_help();
	 exit(0);
      }
      else
      {
	 cout << "unrecognized option " << argv[a] << endl;
	 a++;
      }
   }
   // Initialize pseudo-random number generator
   if( !seedread )
      getseed( seed );
   initrandom( seed );


   MPI_Init(&argc, &argv);
   double tstart = MPI_Wtime();
   int myid, nprocs;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   if( procdist == 3 )
   {
      if( nprocs != p1*p2*p3 )
      {
	 cout << "Error, the given procdist " << p1 << " " << p2 << " " << p3 <<
	    " is not equal to the total number of processors " << nprocs << endl;
	 exit(0);
      }
   }
   if( myid == 0 )
   {
      int version,subversion;
      MPI_Get_version(&version,&subversion);
      cout << "Running MPI version " << version <<"."<<subversion << endl;
   }
   bool hdf5 = fname.compare(fname.size()-5,5,".hdf5")==0;
#ifndef IOTEST_HDF5
   if( hdf5 )
   {
      cout << "Can not read file " << fname << " testio not compiled with HDF5 support" << endl;
      exit(0);
   }
#endif
   int ni, nj, nk;
   if( myid == 0 && !hdf5 )
      get_size( fname.c_str(), ni, nj, nk );
   if( hdf5 )
   {
      get_size_hdf5( fname, ni, nj, nk, myid );
      ionr = 2;
   }
   if( myid == 0 )
   {
      if( hdf5 )
	 cout << "Read *.hdf5 file with dimensions " << ni << " x " << nj << " x " << nk << endl;
      else
	 cout << "Read *.bin file with dimensions " << ni << " x " << nj << " x " << nk << endl;
      cout << " -->  size = " << 1e-9*ni*nj*nk*sizeof(double) << " G bytes " << endl;
   }
   int dim[3] ={ni,nj,nk}, ierr;
   MPI_Bcast( dim, 3, MPI_INT, 0, MPI_COMM_WORLD );
   ni = dim[0];
   nj = dim[1];
   nk = dim[2];

// Every wfreq:th processor will read
   int iread = 0;
   if( myid % wfreq == 0 )
      iread = 1;

// Distribute the dimensions
   int start[3], locdim[3], proc_max[3];
   if( procdist == 1 )
   {
	 // 2D WPP/SW4 type processor mesh
	 proc_decompose_2d( ni, nj, nprocs, olap, proc_max );
	 proc_max[2] = 1;
   }
   else if( procdist == 2 )
   {
	 // 3D general processor mesh
	 proc_decompose_3d( ni, nj, nk, nprocs, olap, proc_max );
   }
   else if( procdist == 3 )
   {
      proc_max[0] = p1;
      proc_max[1] = p2;
      proc_max[2] = p3;
   }
   if( myid == 0 )
      cout << "processor mesh has " << proc_max[0] << " x "<< proc_max[1] << " x " << proc_max[2] << " processors" << endl;
      
   empty = false;
   if( arraydist == 1 )
   {
    // full array w. or wo. overlap
      int myproc[3];
      myproc_3d( myid, proc_max, myproc );
      int s, e;
      decomp1d( ni, myproc[0]-1, proc_max[0], olap, s, e );
      start[0]  = s;
      locdim[0] = e-s+1;
      decomp1d( nj, myproc[1]-1, proc_max[1], olap, s, e );
      start[1]  = s;
      locdim[1] = e-s+1;
      decomp1d( nk, myproc[2]-1, proc_max[2], olap, s, e );
      start[2]  = s;
      locdim[2] = e-s+1;
   }
   else if( arraydist == 2 )
   {
      // Random 1. Perturbation of uniform array distribution
      int myproc[3];
      myproc_3d( myid, proc_max, myproc );
      //      cout << " myproc3d " << myproc[0] << " " << myproc[1] << " " << myproc[2] << endl;
      int s, e;
      decomp1d( ni, myproc[0]-1, proc_max[0], olap, s, e );
      start[0]  = s;
      locdim[0] = e-s+1;
      decomp1d( nj, myproc[1]-1, proc_max[1], olap, s, e );
      start[1]  = s;
      locdim[1] = e-s+1;
      decomp1d( nk, myproc[2]-1, proc_max[2], olap, s, e );
      start[2]  = s;
      locdim[2] = e-s+1;

      double th=1+relsize*(2*randomnr()-1);
      start[0] = static_cast<int>(round(th*start[0]));
      th=1+relsize*(2*randomnr()-1);
      start[1] = static_cast<int>(round(th*start[1]));
      th=1+relsize*(2*randomnr()-1);
      start[2] = static_cast<int>(round(th*start[2]));
      th=1+relsize*(2*randomnr()-1);
      locdim[0] = static_cast<int>(round(th*locdim[0]));
      th=1+relsize*(2*randomnr()-1);
      locdim[1] = static_cast<int>(round(th*locdim[1]));
      th=1+relsize*(2*randomnr()-1);
      locdim[2] = static_cast<int>(round(th*locdim[2]));

      if( locdim[0] < 1 )
	 locdim[0] = round(randomnr()*10)+1;
      if( locdim[1] < 1 )
	 locdim[1] = round(randomnr()*10)+1;
      if( locdim[2] < 1 )
	 locdim[2] = round(randomnr()*10)+1;

      if( start[0] < 1 )
	 start[0] = 1;
      if( start[0] > ni )
	 start[0] = ni;
      if( start[1] < 1 )
	 start[1] = 1;
      if( start[1] > nj )
	 start[1] = nj;
      if( start[2] < 1 )
	 start[2] = 1;
      if( start[2] > nk )
	 start[2] = nk;

      if( start[0]+locdim[0]-1 > ni )
	 locdim[0]= ni-start[0]+1;
      if( start[1]+locdim[1]-1 > nj )
	 locdim[1]= nj-start[1]+1;
      if( start[2]+locdim[2]-1 > nk )
	 locdim[2]= nk-start[2]+1;
   }
   else if( arraydist == 3 )
   {
      // Random 2. Completely random with some empty processors
      if( randomnr() > efraction )
      {
	 // Non-empty, do regular decomposition to get locdim
	 int myproc[3];
	 myproc_3d( myid, proc_max, myproc );
	 int s, e;
	 decomp1d( ni, myproc[0]-1, proc_max[0], olap, s, e );
	 start[0]  = s;
	 locdim[0] = e-s+1;
	 decomp1d( nj, myproc[1]-1, proc_max[1], olap, s, e );
	 start[1]  = s;
	 locdim[1] = e-s+1;
	 decomp1d( nk, myproc[2]-1, proc_max[2], olap, s, e );
	 start[2]  = s;
	 locdim[2] = e-s+1;
	 start[0] =static_cast<int>(round(1 + randomnr()*dim[0]));
	 start[1] =static_cast<int>(round(1 + randomnr()*dim[1]));
	 start[2] =static_cast<int>(round(1 + randomnr()*dim[2]));
	 double th=1+relsize*(2*randomnr()-1);
	 locdim[0] = static_cast<int>(round(th*locdim[0]));
	 th=1+relsize*(2*randomnr()-1);
	 locdim[1] = static_cast<int>(round(th*locdim[1]));
	 th=1+relsize*(2*randomnr()-1);
	 locdim[2] = static_cast<int>(round(th*locdim[2]));
         if( start[0] < 1 )
	    start[0] = 1;
	 if( start[0] > ni )
	    start[0] = ni;
         if( start[1] < 1 )
	    start[1] = 1;
	 if( start[1] > nj )
	    start[1] = nj;
         if( start[2] < 1 )
	    start[2] = 1;
	 if( start[2] > nk )
	    start[2] = nk;

         if( start[0]+locdim[0]-1 > ni )
	    locdim[0]= ni-start[0]+1;
         if( start[1]+locdim[1]-1 > nj )
	    locdim[1]= nj-start[1]+1;
         if( start[2]+locdim[2]-1 > nk )
	    locdim[2]= nk-start[2]+1;
      }
      else
      {
         start[0] = start[1] = start[2] = 0;
	 locdim[0] = locdim[1] = locdim[2] = -1;
	 empty = true;
      }
   }
   size_t npts  = (static_cast<size_t>(locdim[0]))*locdim[1]*locdim[2];
   if( printdim )
   {
      cout << "proc " << myid << "start  i,j,k: " << start[0] << " " << start[1] << " " << start[2] << endl;
      cout << "proc " << myid << "locdim i,j,k: " << locdim[0] << " " << locdim[1] << " " << locdim[2] << endl;
      cout << "proc " << myid << " number of points " << npts << endl;
   }

   double* w;
   try {
      if( !empty )
	 w = new double[npts];
      else
	 w = new double[1];
   }
   catch( bad_alloc& ba )
   {
      cout << "Processor " << myid << " memory allocation failed. npts = " << npts <<
	 " Exception= " << ba.what() << endl;
      MPI_Abort(MPI_COMM_WORLD,0);
   }

// Read the file
   MPI_Barrier(MPI_COMM_WORLD);
   double t[4];
   t[0] = tstart;
   int fd;
   if( !hdf5 )
      fd = open( fname.c_str(), O_RDONLY );
   size_t pos0=12;
   int lgmap[3]={start[0]-1,start[1]-1,start[2]-1};
   t[1] = MPI_Wtime();
   if( ionr == 0 )
   {
      Parallel_IO* pio = new Parallel_IO( iread, 1, dim, locdim, lgmap, nptsbuf );
      t[2] = MPI_Wtime();
      pio->read_array( &fd, 1, w, pos0, "double" );
      t[3] = MPI_Wtime();
   }
   else if( ionr == 1 )
   {
      t[2] = MPI_Wtime();
      read_wp_patch_f( &fd, 1, locdim[0], locdim[1], locdim[2], ni, nj, nk, 8, w, lgmap,
		       pos0 );
      t[3] = MPI_Wtime();
   }
   else if( ionr == 2 )
   {
      t[2] = MPI_Wtime();
      read_hdf5_patch_f( fname, 1, locdim[0], locdim[1], locdim[2], ni, nj, nk, 8, w, lgmap,
		       pos0 );
      t[3] = MPI_Wtime();
   }
   if( !hdf5 )
      close(fd);

   // Check the result
   double err;
   double hx=1.0/(ni-1), hy=1.0/(nj-1), hz=1.0/(nk-1);
   check_result( locdim, w, start, hx, hy, hz, err, hdf5 );
   double errtmp=err;

   MPI_Allreduce( &errtmp, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   if( myid == 0 )
   {
      char esc = 27;
      if( err < 1e-10 )
	 cout << esc << "[32m" << "Result is correct" << esc << "[0m" << ", max error " << err << endl;
      else
	 cout << esc << "[31m" << "Result is not correct" << esc << "[0m" << ", max error " << err << endl;
   }
   delete[] w;
   double* timing=0;
   if( myid == 0 )
      timing = new double[4*nprocs];
   MPI_Gather( t, 4, MPI_DOUBLE, timing, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   if( myid == 0 )
   {
      double t1min=1e38,t3max=0, t2max=0;
      for( int p = 0 ; p < nprocs ; p++ )
      {
	 if( timing[1+4*p] <t1min )
	    t1min = timing[1+4*p];
	 if( timing[3+4*p] > t3max )
	    t3max = timing[3+4*p];
      }
      cout << "Read time is " << t3max-t1min << " seconds " << endl;
      if( savedata )
      {
	 ofstream timefile;
	 if( appendtofile )
	    timefile.open("timings.dat",ofstream::app);
	 else
	    timefile.open("timings.dat");
	 timefile << nprocs << " " << t3max-t1min << " " << ni << " " << nj << " " << nk << " " << proc_max[0] << " " << proc_max[1] << " " << proc_max[2] << endl;
	 timefile.close();
	 
      }
   }
   MPI_Finalize();
}
