
#include <fcntl.h>
#include <iostream>
#include <fstream>


using namespace std;

int convert( string infname, string outfname, int verbose )
{
   int prec, npatches, plane,  mode, gridinfo;
   double t, coord;
   char timestr[25];


   int fd = open(infname.c_str(),O_RDONLY);
   if( fd == -1 )
   {
      cout << "Error: Could not open file " << infname << endl;
      return -1;
   }
   size_t nr = read(fd,&prec,sizeof(int));
   if( nr != sizeof(int) )
   {
      cout << "Error reading precision" << endl;
      close(fd);
      return -1;
   }
   nr = read(fd,&npatches,sizeof(int));
   if( nr != sizeof(int) )
   {
      cout << "Error reading npatches" << endl;
      close(fd);
      return -1;
   }
   nr = read(fd,&t,sizeof(double));
   if( nr != sizeof(double) )
   {
      cout << "Error reading time" << endl;
      close(fd);
      return -1;
   }
   nr = read(fd,&plane,sizeof(int));
   if( nr != sizeof(int) )
   {
      cout << "Error reading plane" << endl;
      close(fd);
      return -1;
   }
   nr = read(fd,&coord,sizeof(double));
   if( nr != sizeof(coord) )
   {
      cout << "Error reading coord" << endl;
      close(fd);
      return -1;
   }
   nr = read(fd,&mode,sizeof(int));   
   if( nr != sizeof(int) )
   {
      cout << "Error reading mode" << endl;
      close(fd);
      return -1;
   }
   nr = read(fd,&gridinfo,sizeof(int));   
   if( nr != sizeof(int) )
   {
      cout << "Error reading gridinfo" << endl;
      close(fd);
      return -1;
   }
   nr = read(fd,&timestr[0],25*sizeof(char));
   timestr[24] = '\0';
   if( nr != 25*sizeof(char) )
   {
      cout << "Error reading timestr" << endl;
      close(fd);
      return -1;
   }

   if( verbose == 1 )
   {
      cout << "Read: prec = " << prec << " t= " << t << " plane= " << plane << endl;
      //      cout << "     coord = " coord <<  " mode= " << mode <<  " npatches= " << npatches<< endl;
      cout << "      coord = " << coord <<  " npatches= " << npatches<< " gridinfo= " << gridinfo << endl;
      cout << "      file created " << timestr << endl;
   }

// Check values for reasonableness
   if( prec != 4 && prec != 8 )
   {
      cout << "Error precision must be 4 or 8, not " << prec << endl;
      close(fd);
      return -1;
   }
   if( npatches <= 0 )
   {
      cout << "Error npatches must be positive, not " << npatches << endl;
      close(fd);
      return -1;
   }
   if( !(plane==0 || plane == 1 || plane ==2 ) )
   {
      cout << "Error plane must be 0,1, or 2, not " << plane << endl;
      close(fd);
      return -1;
   }

   double *h = new double[npatches];
   double *zmin=new double[npatches];
   int *ib = new int[npatches];
   int *ni = new int[npatches];
   int *jb = new int[npatches];
   int *nj = new int[npatches];
   for( int p= 0 ; p < npatches ; p++ )
   {
      nr = read( fd, h+p, sizeof(double));
      if( nr != sizeof(double) )
      {
	 cout << "Error reading h on patch " << p+1 << endl;
	 close(fd);
	 return -1;
      }

      nr = read( fd, zmin+p, sizeof(double));
      if( nr != sizeof(double) )
      {
	 cout << "Error reading zmin on patch " << p+1 << endl;
	 close(fd);
	 return -1;
      }

      nr = read( fd, ib+p, sizeof(int));
      if( nr != sizeof(int) )
      {
	 cout << "Error reading ib on patch " << p+1 << endl;
	 close(fd);
	 return -1;
      }
      nr = read( fd, ni+p, sizeof(int));
      if( nr != sizeof(int) )
      {
	 cout << "Error reading ni on patch " << p+1 << endl;
	 close(fd);
	 return -1;
      }
      nr = read( fd, jb+p, sizeof(int));
      if( nr != sizeof(int) )
      {
	 cout << "Error reading jb on patch " << p+1 << endl;
	 close(fd);
	 return -1;
      }
      nr = read( fd, nj+p, sizeof(int));
      if( nr != sizeof(int) )
      {
	 cout << "Error reading nj on patch " << p+1 << endl;
	 close(fd);
	 return -1;
      }
      if( verbose == 1 )
         cout << "    patch nr " << p+1 << " has h= " << h[p] << " zmin = " << zmin[p] << endl;
   }
   ofstream ofile(outfname.c_str());
   ofile << npatches << endl;
   for( int p= 0 ; p < npatches ; p++ )
      ofile << ib[p] << " " << ni[p] << " " << jb[p] << " " << nj[p] << " " << h[p] << " " << zmin[p] << endl;

   if( prec == 4 )
   {
      float* im;
      float* zgrid;
      for( int p= 0 ; p < npatches ; p++ )
      {
         size_t npts = (ni[p]-ib[p]+1)*(nj[p]-jb[p]+1);
	 im = new float[npts];
	 nr = read(fd,im,npts*sizeof(float));
         if( nr != npts*sizeof(float) )
	 {
            cout << "Error reading image patch nr " << p+1 << endl;
	    close(fd);
	    return -1;
	 }
         bool gridknown = false;
	 if( p == npatches-1 && gridinfo == 1 )
	 {
            zgrid = new float[npts];
	    nr = read(fd,zgrid,npts*sizeof(float));
	    if( nr != npts*sizeof(float) )
	    {
	       cout << "Error reading grid patch nr " << p+1 << endl;
	       close(fd);
	       return -1;
	    }
            gridknown = true;
	 }
         size_t npti   = ni[p]-ib[p]+1;
	 size_t offset = -ib[p] - jb[p]*npti;
	 for( int j=jb[p] ; j <= nj[p]  ; j++ )
	    for( int i=ib[p] ; i <= ni[p]  ; i++ )
	    {
	       double x = (i-1)*h[p];
	       double y = (j-1)*h[p];
               if( plane == 0 || plane == 1 )
	       {
		  if( gridknown )
		     y = zgrid[i+npti*j+offset];
		  else
		     y = y + zmin[p];
	       }
	       ofile << x << " " << y << " " << im[i+npti*j+offset] << endl;
	    }
 	 delete[] im;
         if( gridknown )
	    delete[] zgrid;
      }
   }
   else if( prec == 8 )
   {
      double* im;
      double* zgrid;
      for( int p= 0 ; p < npatches ; p++ )
      {
         size_t npts = (ni[p]-ib[p]+1)*(nj[p]-jb[p]+1);
	 im = new double[npts];
	 nr = read(fd,im,npts*sizeof(double));
         if( nr != npts*sizeof(double) )
	 {
            cout << "Error reading image patch nr " << p+1 << endl;
	    close(fd);
	    return -1;
	 }
         bool gridknown = false;
	 if( p == npatches-1 && gridinfo == 1 )
	 {
            zgrid = new double[npts];
	    nr = read(fd,zgrid,npts*sizeof(double));
	    if( nr != npts*sizeof(double) )
	    {
	       cout << "Error reading grid patch nr " << p+1 << endl;
	       close(fd);
	       return -1;
	    }
            gridknown = true;
	 }
         size_t npti   = ni[p]-ib[p]+1;
	 size_t offset = -ib[p] - jb[p]*npti;
	 for( int j=jb[p] ; j <= nj[p]  ; j++ )
	    for( int i=ib[p] ; i <= ni[p]  ; i++ )
	    {
	       double x = (i-1)*h[p];
	       double y = (j-1)*h[p];
               if( plane == 0 || plane == 1 )
	       {
		  if( gridknown )
		     y = zgrid[i+npti*j+offset];
		  else
		     y = y + zmin[p];
	       }
	       ofile << x << " " << y << " " << im[i+npti*j+offset] << endl;
	    }
	 delete[] im;
         if( gridknown )
	    delete[] zgrid;
      }
   }
   ofile.close();  
   close(fd);
   delete[] h;   
   delete[] zmin;
   delete[] ib;
   delete[] ni;
   delete[] jb;
   delete[] nj;
   return 0;
}

//-----------------------------------------------------------------------
int main( int argc, char** argv )
{
   string infname, outfname;
   int verbose=1;

   if( argc > 1 )
   {
      infname = argv[1];
      int n = infname.length();
      if( n > 6 && infname.substr(n-6).compare("sw4img")==0 )
      {	 
         outfname = infname.substr(0,n-6) + "sw4ascimg";
         if( verbose == 1 )
	 {
	    cout << "converting sw4 image file: " << infname << endl;
            cout << "  to ascii sw4 image file: " << outfname << endl;
	 }
	 convert( infname, outfname, verbose );
      }
      else
         cout << "Error, the input file, " << infname << " , is not a sw4 image file" << endl;
   }
   else
   {
      cout << "Usage: readimage imagefile " << endl;
      cout << " The output file will have the same name as the input file, " <<
	 "but with extension *.sw4ascimg instead of the input *.sw4img" << endl;
   }
}
