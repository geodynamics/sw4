#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <math.h>

#include "Sarray.h"
using namespace std;

#include "DataPatches.h"
#include <errno.h>

//-----------------------------------------------------------------------
DataPatches::DataPatches( string fname, Sarray& u, int imin, int imax, int jmin, int jmax,
			  int kmax, int layers, int ntsteps, double dt )
{
   m_filename = fname;
   // Find the patches residing in this processor

 // My block
   int ib=u.m_ib;
   int ie=u.m_ie;
   int jb=u.m_jb;
   int je=u.m_je;
   int kb=u.m_kb;
   int ke=u.m_ke;

   int off = layers-1;

// Bottom block
   m_npatches = 0;
   int wind[6];
   u.intersection(imin-off,imax+off,jmin-off,jmax+off,kmax,kmax+off,wind);
   add_patch(wind);

// Side 3
   u.intersection(imin-off,imax+off,jmin-off,jmin,1,kmax,wind);
   add_patch(wind);

// Side 4
   u.intersection(imin-off,imax+off,jmax,jmax+off,1,kmax,wind);
   add_patch(wind);

// Side 1
   u.intersection(imin-off,imin,jmin,jmax,1,kmax,wind);
   add_patch(wind);

// Side 2
   u.intersection(imax,imax+off,jmin,jmax,1,kmax,wind);
   add_patch(wind);

   m_isnonempty = m_npatches > 0;
   if( m_isnonempty )
   {
      m_error = false;
   // Setup data structures
      m_ncomp = u.m_nc;
      m_dataptr = new size_t[m_npatches+1];
      m_dataptr[0] = 0;
      for( int p=0 ; p < m_npatches ; p++ )
      {
	 size_t npts = m_ncomp*(m_dims[6*p+1]-m_dims[6*p]+1)*(m_dims[6*p+3]-m_dims[6*p+2]+1)
	 *(m_dims[6*p+5]-m_dims[6*p+4]+1);
	 m_dataptr[p+1] = m_dataptr[p] + npts;
      }
      m_nsteps = ntsteps;
      m_data.resize(m_nsteps);
      size_t totdim = m_dataptr[m_npatches];
      for( int n=0 ; n < m_nsteps ; n++ )
	 m_data[n] = new double[totdim];
      m_ncurrent = 0;
      m_steps = new int[m_nsteps];

   // Create file and write header
      int fd = open( m_filename.c_str(), O_CREAT|O_TRUNC|O_WRONLY, 0660 );
      if( fd == -1 )
      {
	 printf("ERROR DataPatches: file %s could not be created \n",m_filename.c_str());
         print_openerr( errno );
	 m_error = true;
      }
      else
      {
	 m_nmin=0;
	 m_nmax=-1;
	 size_t nr = write(fd,&m_nmin,sizeof(int));
	 if( nr != sizeof(int) )
	 {
	    printf("ERROR DataPatches::Constructor : could not write nmin \n");
	    m_error = true;
	 }
	 
	 nr = write(fd,&m_nmax,sizeof(int));
	 if( nr != sizeof(int) )
	 {
	    printf("ERROR DataPatches::Constructor : could not write nmax \n");
	    m_error = true;
	 }
      
	 nr = write(fd,&m_npatches,sizeof(int));
	 if( nr != sizeof(int) )
	 {
	    printf("ERROR DataPatches::Constructor : could not write npatches \n");
	    m_error = true;
	 }

	 nr = write(fd,&m_ncomp,sizeof(int));
	 if( nr != sizeof(int) )
	 {
	    printf("ERROR DataPatches::Constructor : could not write ncomp \n");
	    m_error = true;
	 }

	 nr = write(fd,&m_dims[0], m_npatches*6*sizeof(size_t));
	 if( nr != sizeof(size_t)*6*m_npatches )
	 {
	    printf("ERROR DataPatches::Constructor : could not write dims \n");
	    m_error = true;
	 }

	 nr = write(fd,&dt,sizeof(double));
	 if( nr != sizeof(double) )
	 {
	    printf("ERROR DataPatches::Constructor : could not write dt \n");
	    m_error = true;
	 }
	 close(fd);
      }
   }
   m_startedsave = false;
}

//-----------------------------------------------------------------------
DataPatches::~DataPatches()
{
   if( m_isnonempty )
   {
      delete[] m_steps;
      for( int i=0 ; i < m_data.size() ; i++ )
	 delete[] m_data[i];
      delete[] m_dataptr;
      unlink(m_filename.c_str());
   // unlink deletes the file
   }
}
//-----------------------------------------------------------------------
size_t DataPatches::get_noofpoints() const
{
   if( m_isnonempty )
      return m_dataptr[m_npatches];
   else
      return 0;
}

//-----------------------------------------------------------------------
void DataPatches::add_patch( int wind[6] )
{
   if( wind[1]>=wind[0] && wind[3]>=wind[2] && wind[5]>=wind[4] )
   {
      m_dims.push_back(wind[0]);
      m_dims.push_back(wind[1]);
      m_dims.push_back(wind[2]);
      m_dims.push_back(wind[3]);
      m_dims.push_back(wind[4]);
      m_dims.push_back(wind[5]);
      m_npatches++;
   }
}

//-----------------------------------------------------------------------
void DataPatches::push( Sarray& u, int n )
{
   if( m_isnonempty && !m_error )
   {
      if( m_ncurrent == m_nsteps )
      {
	 save_to_file( );
	 m_ncurrent = 0;
      }
      int ib=u.m_ib;
      int jb=u.m_jb;
      int kb=u.m_kb;
      size_t ni=static_cast<size_t>(u.m_ni);
      size_t nj=static_cast<size_t>(u.m_nj);
      double* uptr=u.c_ptr();
      for( int p=0 ; p < m_npatches ; p++ )
      {
 	 size_t ptr = m_dataptr[p];
	 size_t ind=0, indu;
	 for( int k=m_dims[6*p+4] ; k <= m_dims[6*p+5] ; k++ )
	    for( int j=m_dims[6*p+2] ; j <= m_dims[6*p+3] ; j++ )
	       for( int i=m_dims[6*p] ; i <= m_dims[6*p+1] ; i++ )
	       {
		  indu = (i-ib)+ni*(j-jb)+ni*nj*(k-kb);
		  for( int c= 0 ; c<m_ncomp ; c++ )
		     m_data[m_ncurrent][ptr+m_ncomp*ind+c] = uptr[m_ncomp*indu+c];
		  ind++;
	       }
      }
      m_steps[m_ncurrent] = n;
      m_ncurrent++;
   }
}

//-----------------------------------------------------------------------
void DataPatches::save_to_file( )
{
   if( m_ncurrent > 0 && !m_error && m_isnonempty )
   {
      size_t nr;
      int fd = open( m_filename.c_str(), O_RDWR );
      if( fd == -1 )
      {
	 printf("ERROR DataPatches::save_to_file : file %s could not be opened \n",m_filename.c_str());
         print_openerr( errno );
         return;
      }
      
   // File format: Header,Data_1,Data_2,..,Data_N
   // Where:
   // Header is: nmin,nmax,npatches,ncomp,dims,dt
   // Data_i is the field data associated with time step i.

      if( !m_startedsave )
      {
	 // Set up nmin, nmax,
	 m_nmin = m_steps[0];
	 m_nmax = m_steps[m_ncurrent-1];
         m_startedsave = true;
	 nr = write(fd,&m_nmin,sizeof(int));
	 if( nr != sizeof(int) )
	 {
	    printf("ERROR DataPatches::save_to_file : Failed to write nmin\n");
	    return;
	 }
	 nr = write(fd,&m_nmax,sizeof(int));
	 if( nr != sizeof(int) )
	 {
	    printf("ERROR DataPatches::save_to_file : Failed to write nmax\n");
	    return;
	 }
      }
      else
      {
   // Check if add to sequence of steps is ok ? Update nmax.
	 if( m_steps[0] != m_nmax+1 )
	 {
	    printf("ERROR DataPatches::save_to_file : Trying to store non-consequtive time steps\n");
            return;
	 }
	 else
	    m_nmax = m_steps[m_ncurrent-1];

         nr = lseek(fd,sizeof(int),SEEK_SET);
	 if( nr == -1 )
	 {
	    printf("ERROR DataPatches::save_to_file : Could not find beginning of file\n");
	    return;
	 }
	 nr = write(fd,&m_nmax,sizeof(int));
	 if( nr != sizeof(int) )
	 {
	    printf("ERROR DataPatches::save_to_file : Failed to write nmax\n");
	    return;
	 }
      }

      // Skip to end of file...
      nr = lseek(fd,0,SEEK_END);
      if( nr == -1 )
      {
	 printf("ERROR DataPatches::save_to_file : Could not find end of file\n");
	 return;
      }
	 
      // ...and add new data set(s)
      size_t linuxlimit = static_cast<size_t>(round(pow(2.0,31)/sizeof(double)));
      for( int n=0; n < m_ncurrent ; n++ )
      {
         size_t npts      = m_dataptr[m_npatches];
         size_t nrwritten = 0;
         while( nrwritten < npts )
	 {
            size_t nwri = npts-nrwritten;
	    if( nwri > linuxlimit )
	       nwri = linuxlimit;
	    nr = write(fd,m_data[n]+nrwritten,nwri*sizeof(double));
	    if( nr != nwri*sizeof(double) )
	    {
	       printf("ERROR DataPatches::save_to_file : Could not write data\n");
	       printf("    Attempt to write %ld bytes, only wrote %ld \n",nwri*sizeof(double),nr );
	       return;
	    }
	    nrwritten += nwri;
	 }
      }
      close(fd);
   }
}

//-----------------------------------------------------------------------
void DataPatches::pop( Sarray& u, int n )
{
   if( m_isnonempty && !m_error )
   {
      int nloc;
      if( !(m_steps[0] <= n && n <= m_steps[m_ncurrent-1]) )
      {
      // Step not in memory, read from file
      // Assume solving backwards, read the entire interval:  n-nsteps+1 .., n into 0,..,nsteps-1
	 read_from_file( n );
	 //	 nloc = m_ncurrent-1;
      }
	 // Find index in memory
      nloc = m_ncurrent-1;
      while( nloc >= 0 && m_steps[nloc] != n )
	nloc--;
      //      printf("nloc = %i n= %i step[0] = %i \n",nloc,n,m_steps[0] );
      int ib=u.m_ib;
      int jb=u.m_jb;
      int kb=u.m_kb;
      size_t ni=static_cast<size_t>(u.m_ni);
      size_t nj=static_cast<size_t>(u.m_nj);
      for( int p=0 ; p < m_npatches ; p++ )
      {
	 size_t ptr = m_dataptr[p];
	 size_t ind = 0, indu;
         double* uptr = u.c_ptr();
	 for( int k=m_dims[6*p+4] ; k <= m_dims[6*p+5] ; k++ )
	    for( int j=m_dims[6*p+2] ; j <= m_dims[6*p+3] ; j++ )
	       for( int i=m_dims[6*p] ; i <= m_dims[6*p+1] ; i++ )
	       {
		  indu = (i-ib)+ni*(j-jb)+ni*nj*(k-kb);
		  for( int c= 0 ; c<m_ncomp ; c++ )
		     uptr[m_ncomp*indu+c] = m_data[nloc][ptr+m_ncomp*ind+c]; 
		  ind++;
	       }
      }
   }
}

//-----------------------------------------------------------------------
void DataPatches::read_from_file( int n )
{
   if( !m_error && m_isnonempty )
   {
      size_t nr;
      int fd = open( m_filename.c_str(), O_RDONLY );
      if( fd == -1 )
      {
	 printf("ERROR DataPatches::read_from_file : file %s could not be opened \n",m_filename.c_str());
         print_openerr( errno );
         return;
      }
      int nstart = n-m_nsteps+1;
      if( nstart < m_nmin )
	nstart = m_nmin;
      int nmax = nstart + m_nsteps-1;
      
      // header
      size_t off = 4*sizeof(int) + 6*m_npatches*sizeof(size_t) + sizeof(double);
      // skip data before nstart
      off += (nstart - m_nmin )*m_dataptr[m_npatches]*sizeof(double);
      nr = lseek( fd, off, SEEK_SET );

      size_t linuxlimit = static_cast<size_t>(round(pow(2.0,31)/sizeof(double)));

      m_ncurrent = m_nsteps;
      for( int s=0 ; s < m_ncurrent ; s++ )
      {
 	 m_steps[s] = nstart + s;
         size_t npts      = m_dataptr[m_npatches];
         size_t nrread = 0;
         while( nrread < npts )
	 {
            size_t nread = npts-nrread;
	    if( nread > linuxlimit )
	       nread = linuxlimit;
	    nr = read(fd,m_data[s]+nrread,nread*sizeof(double));
	    if( nr != nread*sizeof(double) )
	    {
	       printf("ERROR DataPatches::read_from_file : Could not read data\n");
	       printf("    Attempt to read %ld bytes, only read %ld \n",nread*sizeof(double),nr );
               printf("   file name = %s \n",m_filename.c_str() );
	       return;
	    }
	    nrread += nread;
	 }
      }
      close(fd);
   }
}

//-----------------------------------------------------------------------
void DataPatches::print_openerr( int ecode ) const
{
   string emsg;
   switch( ecode )
   {
   case EACCES: emsg = "EACCESS";break;
   case EEXIST: emsg = "EEXIST";break;
   case EFAULT: emsg = "EFAULT";break;
   case EISDIR: emsg = "EISDIR";break;
   case ELOOP: emsg = "ELOOP";break;
   case EMFILE: emsg = "EMFILE";break;
   case ENAMETOOLONG: emsg = "ENAMETOOLONG";break;
   case ENFILE: emsg = "ENFILE";break;
   case ENODEV: emsg = "ENODEV";break;
   case ENOENT: emsg = "ENOENT";break;
   case ENOMEM: emsg = "ENOMEM";break;
   case ENOSPC: emsg = "ENOSPC";break;
   case ENOTDIR: emsg = "ENOTDIR";break;
   case ENXIO: emsg = "ENOXIO";break;
   case EOVERFLOW: emsg = "EOVERFLOW";break;
   case EPERM: emsg = "EPERM";break;
   case EROFS: emsg = "EROFS";break;
   case ETXTBSY: emsg = "ETXTBSY";break;
   case EWOULDBLOCK:emsg = "EWOULDBLOCK";break;
   default : emsg = "unknown";
   }
   cout << "open failed with error code " << emsg << endl;
}
