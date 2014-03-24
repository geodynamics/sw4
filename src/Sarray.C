#include "Sarray.h"

#include <iostream>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>

//-----------------------------------------------------------------------
Sarray::Sarray( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend )
{
   m_nc = nc;
   m_ib = ibeg;
   m_ie = iend;
   m_jb = jbeg;
   m_je = jend;
   m_kb = kbeg;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
//   m_mpi_datatype_initialized = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray( int ibeg, int iend, int jbeg, int jend, int kbeg, int kend )
{
   m_nc = 1;
   m_ib = ibeg;
   m_ie = iend;
   m_jb = jbeg;
   m_je = jend;
   m_kb = kbeg;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
//   m_mpi_datatype_initialized = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray( int nc, int iend, int jend, int kend )
{
   m_nc = nc;
   m_ib = 1;
   m_ie = iend;
   m_jb = 1;
   m_je = jend;
   m_kb = 1;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
//   m_mpi_datatype_initialized = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray( int iend, int jend, int kend )
{
   m_nc = 1;
   m_ib = 1;
   m_ie = iend;
   m_jb = 1;
   m_je = jend;
   m_kb = 1;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
//   m_mpi_datatype_initialized = false;
}

//-----------------------------------------------------------------------
Sarray::Sarray()
{
//   m_mpi_datatype_initialized = false;
   m_nc = m_ib = m_ie = m_jb = m_je = m_kb = m_ke = 0;
   m_data = NULL;
}

//-----------------------------------------------------------------------
Sarray::Sarray( const Sarray& u )
{
   m_nc = u.m_nc;
   m_ib = u.m_ib;
   m_ie = u.m_ie;
   m_jb = u.m_jb;
   m_je = u.m_je;
   m_kb = u.m_kb;
   m_ke = u.m_ke;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new double[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
}

//-----------------------------------------------------------------------
Sarray::Sarray( Sarray& u, int nc )
{
   if( nc == -1 )
      m_nc = u.m_nc;
   else
      m_nc = nc;
   m_ib = u.m_ib;
   m_ie = u.m_ie;
   m_jb = u.m_jb;
   m_je = u.m_je;
   m_kb = u.m_kb;
   m_ke = u.m_ke;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
      m_data = new double[m_nc*m_ni*m_nj*m_nk];
   else
      m_data = NULL;
//   m_mpi_datatype_initialized = false;
}

//-----------------------------------------------------------------------
// void Sarray::define( CartesianProcessGrid* cartcomm, int nc )
// {
//   if( m_data != NULL )
//      delete[] m_data;
//    m_nc = nc;

//    // global index ranges:
//    //   m_ib = cartcomm->gNstart;
//    //   m_ie = cartcomm->getGlobalN()+m_ib-1;
//    //   m_jb = cartcomm->gMstart;
//    //   m_je = cartcomm->getGlobalM()+m_jb-1;
//    //   m_kb = cartcomm->gLstart;
//    //   m_ke = cartcomm->getGlobalL()+m_kb-1;

//    // local index ranges in global coordinates:
//    //   m_ib = cartcomm->Nstart - cartcomm->NlowParallelPadding;
//    //   m_ie = cartcomm->Nstart + cartcomm->getLocalN() - cartcomm->NlowParallelPadding - 1;
//    //   m_jb = cartcomm->Mstart - cartcomm->MlowParallelPadding;
//    //   m_je = cartcomm->Mstart + cartcomm->getLocalM() - cartcomm->MlowParallelPadding - 1;
//    //   m_kb = cartcomm->Lstart - cartcomm->LlowParallelPadding;
//    //   m_ke = cartcomm->Lstart + cartcomm->getLocalL() - cartcomm->LlowParallelPadding - 1;

//    // local index ranges in normalized coordinates
//    m_ib = 1;
//    m_ie = cartcomm->getLocalN();
//    m_jb = 1;
//    m_je = cartcomm->getLocalM();
//    m_kb = 1;
//    m_ke = cartcomm->getLocalL();

//    m_ni = m_ie-m_ib+1;
//    m_nj = m_je-m_jb+1;
//    m_nk = m_ke-m_kb+1;
//    //   std::cout << "Sarray dims " << m_nc << " " << m_ni << " " 
//    //	     << m_nj << " " << m_nk << std::endl;
//    m_data = new double[m_nc*m_ni*m_nj*m_nk];
// }

//-----------------------------------------------------------------------
void Sarray::define( int nc, int iend, int jend, int kend )
{
   if( m_data != NULL )
      delete[] m_data;

   m_nc = nc;
   m_ib = 1;
   m_ie = iend;
   m_jb = 1;
   m_je = jend;
   m_kb = 1;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
//   m_mpi_datatype_initialized = false;
}

//-----------------------------------------------------------------------
void Sarray::define( int iend, int jend, int kend )
{
   if( m_data != NULL )
      delete[] m_data;

   m_nc = 1;
   m_ib = 1;
   m_ie = iend;
   m_jb = 1;
   m_je = jend;
   m_kb = 1;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
//   m_mpi_datatype_initialized = false;
}

//-----------------------------------------------------------------------
void Sarray::define( int nc, int ibeg, int iend, int jbeg, int jend, int kbeg,
		     int kend )
{
   if( m_data != NULL )
      delete[] m_data;
   m_nc = nc;
   m_ib = ibeg;
   m_ie = iend;
   m_jb = jbeg;
   m_je = jend;
   m_kb = kbeg;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
}

//-----------------------------------------------------------------------
void Sarray::define( int ibeg, int iend, int jbeg, int jend, int kbeg,
		     int kend )
{
   if( m_data != NULL )
      delete[] m_data;
   m_nc = 1;
   m_ib = ibeg;
   m_ie = iend;
   m_jb = jbeg;
   m_je = jend;
   m_kb = kbeg;
   m_ke = kend;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
}

//-----------------------------------------------------------------------
void Sarray::define( const Sarray& u ) 
{
   if( m_data != NULL )
      delete[] m_data;
   m_nc = u.m_nc;
   m_ib = u.m_ib;
   m_ie = u.m_ie;
   m_jb = u.m_jb;
   m_je = u.m_je;
   m_kb = u.m_kb;
   m_ke = u.m_ke;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   m_data = new double[m_nc*m_ni*m_nj*m_nk];
}

//-----------------------------------------------------------------------
void Sarray::intersection( int ib, int ie, int jb, int je, int kb, int ke, int wind[6] )
{
   wind[0] = max(ib,m_ib);
   wind[1] = min(ie,m_ie);
   wind[2] = max(jb,m_jb);
   wind[3] = min(je,m_je);
   wind[4] = max(kb,m_kb);
   wind[5] = min(ke,m_ke);
}

//-----------------------------------------------------------------------
// side_plane returns the index of the ghost points along side =0,1,2,3,4,5 (low-i, high-i, low-j, high-j, low-k, high-k)
void Sarray::side_plane( int side, int wind[6], int nGhost )
{
   wind[0] = m_ib;
   wind[1] = m_ie;
   wind[2] = m_jb;
   wind[3] = m_je;
   wind[4] = m_kb;
   wind[5] = m_ke;
   if( side == 0 )
     wind[1] = wind[0] + (nGhost-1);
   else if( side == 1 )
     wind[0] = wind[1] - (nGhost-1);
   else if( side == 2 )
     wind[3] = wind[2] + (nGhost-1);
   else if( side == 3 )
     wind[2] = wind[3] - (nGhost-1);
   else if( side == 4 )
     wind[5] = wind[4] + (nGhost-1);
   else
     wind[4] = wind[5] - (nGhost-1);
   // if( side == 0 )
   //    wind[1] = wind[0];
   // else if( side == 1 )
   //    wind[0] = wind[1];
   // else if( side == 2 )
   //    wind[3] = wind[2];
   // else if( side == 3 )
   //    wind[2] = wind[3];
   // else if( side == 4 )
   //    wind[5] = wind[4];
   // else
   //    wind[4] = wind[5];
}

//-----------------------------------------------------------------------
void Sarray::side_plane_fortran( int side, int wind[6], int nGhost )
{
// Fortran arrays are base 1
   wind[0] = 1;
   wind[1] = m_ni;
   wind[2] = 1;
   wind[3] = m_nj;
   wind[4] = 1;
   wind[5] = m_nk;
   if( side == 0 )
     wind[1] = wind[0] + (nGhost-1);
   else if( side == 1 )
     wind[0] = wind[1] - (nGhost-1);
   else if( side == 2 )
     wind[3] = wind[2] + (nGhost-1);
   else if( side == 3 )
     wind[2] = wind[3] - (nGhost-1);
   else if( side == 4 )
     wind[5] = wind[4] + (nGhost-1);
   else
     wind[4] = wind[5] - (nGhost-1);
}

//-----------------------------------------------------------------------
void Sarray::set_to_zero()
{
   for( int i=0 ; i < m_nc*m_ni*m_nj*m_nk ; i++ )
      m_data[i] = 0;
}

//-----------------------------------------------------------------------
void Sarray::set_to_minusOne()
{
   for( int i=0 ; i < m_nc*m_ni*m_nj*m_nk ; i++ )
      m_data[i] = -1.;
}

//-----------------------------------------------------------------------
void Sarray::set_value( double scalar )
{
   for( int i=0 ; i < m_nc*m_ni*m_nj*m_nk ; i++ )
      m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void Sarray::set_to_random( double llim, double ulim )
{
   for( int i=0 ; i<m_nc*m_ni*m_nj*m_nk ; i++ )
      m_data[i] = llim + (ulim-llim)*drand48();
}

//-----------------------------------------------------------------------
bool Sarray::in_domain( int i, int j, int k )
{
   return m_ib <= i && i <= m_ie && m_jb <= j && j <= m_je
      && m_kb <= k && k <= m_ke;
}

//-----------------------------------------------------------------------
double Sarray::maximum( int c )
{
   int cm = c-1;
   double mx = m_data[cm];
   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
      mx = mx > m_data[cm+i*m_nc] ? mx : m_data[cm+i*m_nc];
   return mx;
}

//-----------------------------------------------------------------------
double Sarray::minimum( int c )
{
   int cm = c-1;
   double mn = m_data[cm];
   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
      mn = mn < m_data[cm+i*m_nc] ? mn : m_data[cm+i*m_nc];
   return mn;
}

//-----------------------------------------------------------------------
double Sarray::sum( int c )
{
   int cm = c-1;
   double s = 0;
   for( int i=0 ; i<m_ni*m_nj*m_nk ; i++ )
      s += m_data[cm+i*m_nc];
   return s;
}

//-----------------------------------------------------------------------
size_t Sarray::count_nans()
{
   size_t retval = 0;
   size_t npts = m_nc*m_ni*static_cast<size_t>(m_nj)*m_nk;
   for( size_t ind = 0; ind < npts ; ind++)
      if( std::isnan(m_data[ind]) )
	 retval++;
   return retval;
}

//-----------------------------------------------------------------------
size_t Sarray::count_nans( int& cfirst, int& ifirst, int& jfirst, int& kfirst )
{
   cfirst = ifirst = jfirst = kfirst = 0;
   size_t retval = 0, ind=0;
   for( int k=m_kb ; k<=m_ke ; k++ )
      for( int j=m_jb ; j<=m_je ; j++ )
	 for( int i=m_ib ; i <= m_ie ; i++ )
	    for( int c=1 ; c <= m_nc ; c++ )
	    {
	       if( std::isnan(m_data[ind]) )
	       {
		  if( retval == 0 )
		  {
		     ifirst = i;
		     jfirst = j;
		     kfirst = k;
		     cfirst = c;
		  }
		  retval++;
	       }
	       ind++;
	    }
   return retval;
}

//-----------------------------------------------------------------------
void Sarray::copy( const Sarray& u )
{
   if( m_data != NULL )
      delete[] m_data;

   m_nc = u.m_nc;
   m_ib = u.m_ib;
   m_ie = u.m_ie;
   m_jb = u.m_jb;
   m_je = u.m_je;
   m_kb = u.m_kb;
   m_ke = u.m_ke;
   m_ni = m_ie-m_ib+1;
   m_nj = m_je-m_jb+1;
   m_nk = m_ke-m_kb+1;
   if( m_nc*m_ni*m_nj*m_nk > 0 )
   {
      m_data = new double[m_nc*m_ni*m_nj*m_nk];
      for( int i=0 ; i < m_nc*m_ni*m_nj*m_nk ; i++ )
	 m_data[i] = u.m_data[i];
   }
   else
      m_data = NULL;

}

//-----------------------------------------------------------------------
void Sarray::insert_subarray( int ib, int ie, int jb, int je, int kb,
			      int ke, double* ar )
{
   // Assuming nc is the same for m_data and subarray ar.
   int nis = ie-ib+1;
   int njs = je-jb+1;
   int nks = ke-kb+1;
   size_t sind=0, ind=0;
   for( int k=kb ; k<=ke ; k++ )
      for( int j=jb ; j<=je ; j++ )
	 for( int i=ib ; i <= ie ; i++ )
	    {
               sind = (i-ib)  +  nis*(j-jb)   +  nis*njs*(k-kb);
               ind = (i-m_ib) + m_ni*(j-m_jb) + m_ni*m_nj*(k-m_kb);
	       for( int c=1 ; c <= m_nc ; c++ )
		  m_data[ind*m_nc+c-1] = ar[sind*m_nc+c-1];
	    }
}

//-----------------------------------------------------------------------
void Sarray::insert_subarray( int ib, int ie, int jb, int je, int kb,
			      int ke, float* ar )
{
   // Assuming nc is the same for m_data and subarray ar.
   int nis = ie-ib+1;
   int njs = je-jb+1;
   int nks = ke-kb+1;
   size_t sind=0, ind=0;
   for( int k=kb ; k<=ke ; k++ )
      for( int j=jb ; j<=je ; j++ )
	 for( int i=ib ; i <= ie ; i++ )
	    {
               sind = (i-ib) + nis*(j-jb) + nis*njs*(k-kb);
               ind = (i-m_ib) + m_ni*(j-m_jb) + m_ni*m_nj*(k-m_kb);
	       for( int c=1 ; c <= m_nc ; c++ )
		  m_data[ind*m_nc+c-1] = (double)ar[sind*m_nc+c-1];
	    }
}

//-----------------------------------------------------------------------
void Sarray::save_to_disk( const char* fname )
{
   int fd = open(fname, O_CREAT | O_TRUNC | O_WRONLY, 0660 );
   if( fd == -1 )
      std::cout << "ERROR opening file" << fname << " for writing " << std::endl;
   size_t nr = write(fd,&m_nc,sizeof(int));
   if( nr != sizeof(int) )
      std::cout << "Error saving nc to " << fname << std::endl;
   nr = write(fd,&m_ni,sizeof(int));
   if( nr != sizeof(int) )
      std::cout << "Error saving ni to " << fname << std::endl;
   nr = write(fd,&m_nj,sizeof(int));
   if( nr != sizeof(int) )
      std::cout << "Error saving nj to " << fname << std::endl;
   nr = write(fd,&m_nk,sizeof(int));
   if( nr != sizeof(int) )
      std::cout << "Error saving nk to " << fname << std::endl;
   size_t npts = m_nc*( (size_t)m_ni)*m_nj*( (size_t)m_nk);
   nr = write(fd,m_data,sizeof(double)*npts);
   if( nr != sizeof(double)*npts )
      std::cout << "Error saving data array to " << fname << std::endl;
   close(fd);
}

//-----------------------------------------------------------------------
void Sarray::assign( const double* ar )
{
   for( size_t i=0 ; i < m_ni*((size_t) m_nj)*m_nk*m_nc ; i++ )
      m_data[i] = ar[i];
}

//-----------------------------------------------------------------------
void Sarray::assign( const float* ar )
{
   for( size_t i=0 ; i < m_ni*((size_t) m_nj)*m_nk*m_nc ; i++ )
     m_data[i] = (double) ar[i];
}

//-----------------------------------------------------------------------
// void Sarray::write( char* filename, CartesianProcessGrid* cartcomm,
// 		    std::vector<double> pars )
// {

//  // Open file
//    MPI_File fh=0;
//    int ierr = MPI_File_open( MPI::COMM_WORLD, filename,
// 			     MPI_MODE_CREATE|MPI_MODE_WRONLY,
// 			     MPI_INFO_NULL, &fh);
//    if( ierr != MPI_SUCCESS )
//    {
//      char buffer[MPI_MAX_ERROR_STRING];
//      int resultlen;
//      MPI_Error_string( ierr, buffer, &resultlen);
//      VERIFY2(0, "Sarray::write MPI Open Error: " << buffer);
//    }
//    MPI_Status status;
//    MPI_Offset oh = 0;

//    // Header, precision,ni,nj,hr,ht,rmin,time
//    int iheader[4];
//    iheader[0] = cartcomm->getGlobalN();
//    iheader[1] = cartcomm->getGlobalM();
//    iheader[2] = 8;
//    iheader[3] = pars.size();

//    MPI_File_write_at( fh, oh, iheader, 4, MPI::INT, &status );
//    oh += 4*sizeof(int);

//    MPI_File_write_at( fh, oh, &pars[0], pars.size(), MPI::DOUBLE, &status );
//    oh += pars.size()*sizeof(double);

//    if( !m_mpi_datatype_initialized )
//       init_mpi_datatype( cartcomm );

//    MPI_File_set_view( fh, oh, MPI::DOUBLE, m_local_block_type, "native",
// 		      MPI_INFO_NULL );

//    MPI_File_write_all(fh, m_data, m_nc*npts(), MPI::DOUBLE, &status);
//    MPI_File_close(&fh);
// }

// //-----------------------------------------------------------------------
// void Sarray::init_mpi_datatype( CartesianProcessGrid* cartcomm )
// {
//    // local sizes
//    int subsizes[4] = { m_nc, m_ie-m_ib+1, m_je-m_jb+1, m_ke-m_kb+1 };

//    // global sizes
//    int sizes[4]     = {m_nc, cartcomm->getGlobalN(), cartcomm->getGlobalM(), 1 };

//    // staring point of local array in global index space
//    int starts[4]   = { 0, m_ib-1, m_jb-1, m_kb-1};

//    MPI_Type_create_subarray(4, sizes, subsizes, starts, 
// 			    MPI_ORDER_FORTRAN, MPI_DOUBLE,
// 			    &m_local_block_type);

//    MPI_Type_commit(&m_local_block_type);

//    m_mpi_datatype_initialized = true;
// }
