#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#ifdef IOTEST_HDF5
#include "hdf5.h"
#endif

using namespace std;

//extern "C" {
//   void getmat_( int*,int*,int*,int*,int*,int*,double*, double*, double*, double* );
//}

void getmatc( int ib, int ie, int jb, int je, int kb, int ke,
	      double hx, double hy, double hz, double* w, bool corder )
{
   double pi= 4*atan(1.0);
   if( corder )
   {
      int nk=ke-kb+1;
      int nj=je-jb+1;
      double* cpizp = new double[nk];
      for( int k=kb ; k <= ke ; k++ )
	 cpizp[k-kb] = cos(pi*(k-1)*hz);
      double* eyp = new double[nj];
      for( int j=jb ; j <= je ; j++ )
	 eyp[j-jb]  = exp(-(j-1)*hy);
      for( int i=ib ; i <= ie ; i++ )
	 for( int j=jb ; j <= je ; j++ )
	    for( int k=kb ; k <= ke ; k++ )
	    {
	       size_t ind = k-kb+nk*( j-jb+nj*(i-ib));
	       w[ind] = (i-1)*hx*eyp[j-jb]*cpizp[k-kb];
	    }
      delete[] eyp;
      delete[] cpizp;
   }
   else
   {
      int ni=ie-ib+1;
      int nj=je-jb+1;
      for( int k=kb ; k <= ke ; k++ )
      {
	 double z = (k-1)*hz;
	 double cpiz = cos(pi*z);
	 for( int j=jb ; j <= je ; j++ )
	 {
	    double y = (j-1)*hy;
	    double ey = exp(-y);
	    for( int i=ib ; i <= ie ; i++ )
	    {
	       double x = (i-1)*hx;
	       size_t ind = i-ib + ni*(j-jb + nj*(k-kb) );
	       w[ind] = (i-1)*hx*ey*cpiz;
	    }
	 }
      }
   }
}

void print_help()
{
   cout << "Usage: " << endl;
   cout << "gentestfile -nowrite -size mbytes -file filename -ndisp nfreq " << endl;
   cout << "         -nowrite  --> Only estimate dimensions, do not write a file" << endl;
   cout << "         -size mbytes  --> mbytes = size of generated file in megabytes (default 1000 = 1G byte)" << endl;
   cout << "         -file filename --> Name of file to generate, use extension *.hdf5 for hdf5 format (default testfile.bin)"  << endl;
   cout << "         -ndisp nfreq --> Print progress when generating file every nfreq:th grid plane (default 100)" << endl;
}

int main( int argc, char** argv )
{
   // Desired size of file in Mbyte:
   size_t target=1000;
   int ndisp = 100;
   bool do_write=true;
   string fname = "testfile.bin";
   
   int a=1;
   while( a < argc )
   {
      if( strcmp(argv[a],"-nowrite") == 0 )
      {
	 do_write = false;
         a++;
      }
      else if( strcmp(argv[a],"-size") == 0 )
      {
         target = atof(argv[a+1]);
	 a+=2;
      }
      else if( strcmp(argv[a],"-ndisp") == 0 )
      {
         ndisp = atoi(argv[a+1]);
	 a+=2;
      }
      else if( strcmp(argv[a],"-file") == 0 )
      {
         fname = argv[a+1];
	 a+=2;
      }
      else if( strcmp(argv[a],"-h") == 0 )
      {
	 print_help();
	 exit(0);
	 a++;
      }
      else
      {
	 cout << "unrecoginzed option " << argv[a] << endl;
	 a++;
      }
   }
   bool hdf5 = fname.compare(fname.size()-5,5,".hdf5")==0;
#ifndef IOTEST_HDF5
   if( hdf5 )
   {
      cout << "Can not create file " << fname << " gentestfile not compiled with HDF5 support" << endl;
      exit(0);
   }
#endif
      
   // Total number of grid points
   double nbase = pow(0.125*target,1.0/3.0)*100;

   // Split into three arrays with different size in each direction
   size_t nj = static_cast<int>(round(nbase));
   size_t ni = static_cast<int>(round(1.414*nbase));
   size_t nk = static_cast<int>(round(0.707*nbase));

   if( !do_write )
   {
      cout << "  Estimated file size is " << 1e-9*ni*nj*nk*sizeof(double) << " G bytes " << endl;
      cout << "  Array dimensions  " << ni << " x " << nj << " x " << nk << endl;
   }
   else
   {
      cout << "Generating array of dimensions " << ni << " x " << nj << " x " << nk << endl;
      double hx = 1.0/(ni-1);
      double hy = 1.0/(nj-1);
      double hz = 1.0/(nk-1);
      int ni32 = ni;
      int nj32 = nj;
      int nk32 = nk;
      struct timeval tv;
      struct timezone tz;

      if( !hdf5 )
      {
	 // Write plain binary format
	 int fd=open( fname.c_str(), O_CREAT | O_TRUNC | O_WRONLY , 0660 );
	 if( fd == -1 )
	 {
	    cout << "ERROR: could not open file " << fname << " for writing" << endl;
	    exit(-1);
	 }
	 else
	 {
	    cout << "Writing file of size " << 1e-9*ni*nj*nk*sizeof(double) << " G bytes " << endl;

	    double* w = new double[ni*nj];

	    gettimeofday( &tv, &tz );
	    double t1 = tv.tv_sec+tv.tv_usec*1e-6;
	 // Write header
	    size_t nr = write( fd, &ni32, sizeof(int) );
	    nr = write( fd, &nj32, sizeof(int) );
	    nr = write( fd, &nk32, sizeof(int) );
	    for( int k=1 ; k <= nk ; k++ )
	    {
	       //	       int o=1;
	       //	       getmat_( &o, &ni32, &o, &nj32, &k, &k, &hx, &hy, &hz, w );
	       getmatc(1,ni32,1,nj32,k,k,hx,hy,hz,w,false);
	       size_t nr=write(fd,w,sizeof(double)*ni*nj);
	       if( nr != sizeof(double)*ni*nj )
	       {
		  cout << "ERROR: could not write plane " << k << endl;
		  cout << " Number of bytes written= " << nr << " Number of bytes in plane= "
		       << sizeof(double)*ni*nj << endl;
		  break;
	       }
	       if( k % ndisp == 0 )
		  cout << " written to k= " << k << endl;
	    }
	    close(fd);
	    gettimeofday( &tv, &tz );
	    cout << " time writing is " << tv.tv_sec+tv.tv_usec*1e-6 - t1 << " seconds " << endl;
	 }
      }
      else
      {
#ifdef IOTEST_HDF5
	 // Write hdf5 format
	 //Create HDF5 file
	 hid_t file_id = H5Fcreate( fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
	 if( file_id < 0 )
	 {
	    cout << "H5Fcreate: Error opening file " << fname << " for writing " << endl;
	    exit(-1);
	 }
	 // Define HDF5 data space
	 hsize_t dims[3]={ni32,nj32,nk32};
	 hid_t dataspace = H5Screate_simple(3, dims, NULL );
	 // Create HDF5 data set      
	 hid_t dataset = H5Dcreate(file_id, "testdata", H5T_NATIVE_DOUBLE, dataspace,
			     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	 hsize_t slicedims[3]={1,nj32,nk32};
	 hid_t slice_id = H5Screate_simple( 3, slicedims, NULL );
	 hsize_t stride[3]={1,1,1};
	 hsize_t block[3] ={1,1,1};

	 double* w = new double[nj*nk];
	 gettimeofday( &tv, &tz );
	 double t1 = tv.tv_sec+tv.tv_usec*1e-6;
	 for( int i=1 ; i <= ni ; i++ )
	 {
	    getmatc(i,i,1,nj32,1,nk32,hx,hy,hz,w,true);
	    hsize_t start[3]={i-1,0,0};
	    herr_t status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, start, stride,
						slicedims, block);
	    if( status == -1 )
	    {
	       cout << "Error from H5Sselect_hyperslab " << endl;
	       exit(-1);
	    }
	    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, slice_id, dataspace,
			      H5P_DEFAULT, w );
	    if( status == -1 )
	    {
	       cout << "Error from H5Dwrite " << endl;
	       exit(-1);
	    }
	    if( i % ndisp == 0 )
	       cout << " written to i= " << i << endl;

	 }
	 H5Sclose(slice_id);
	 H5Dclose(dataset);
	 H5Sclose(dataspace);
	 H5Fclose(file_id);
	 gettimeofday( &tv, &tz );
	 cout << " time writing is " << tv.tv_sec+tv.tv_usec*1e-6 - t1 << " seconds " << endl;
#endif
      }
   }
}
