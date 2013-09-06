#include <mpi.h>
#include <fcntl.h>
#include <cstring>

#include "MaterialParameterization.h"

#include "EW.h"
#include "Require.h"

using namespace std;

//-----------------------------------------------------------------------
MaterialParameterization::MaterialParameterization( EW* a_ew, char* fname )
{
   MPI_Comm_rank( MPI_COMM_WORLD, &m_myrank );
   m_ew = a_ew;
   int n = strlen(fname)+1;
   m_filename = new char[n];
   strcpy(m_filename,fname);
}

//-----------------------------------------------------------------------
void MaterialParameterization::get_nr_of_parameters( int& nms, int& nmd,
						     int& nmd_global ) const
{
   nmd = m_nmd;
   nms = m_nms;
   nmd_global = m_nmd_global;
}

//-----------------------------------------------------------------------
void MaterialParameterization::store_material( vector<Sarray>& a_rho,
					       vector<Sarray>& a_mu,
					       vector<Sarray>& a_lambda )
{
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      m_ew->mRho[g].copy( a_rho[g] );
      m_ew->mMu[g].copy( a_mu[g] );
      m_ew->mLambda[g].copy( a_lambda[g] );
   }
}

//-----------------------------------------------------------------------
void MaterialParameterization::constrain_material( int nmd, double* xmd,
						   int nms, double* xms )
{
   VERIFY2( nmd== m_nmd && nms == m_nms, "MP::constrain_material: Error in sizes");
   int ng = m_ew->mNumberOfGrids;
   vector<Sarray> rho, mu, lambda;
   rho.resize(ng);
   mu.resize(ng);
   lambda.resize(ng);
   for( int g=0; g < ng ; g++ )
   {
      rho[g].define(m_ew->m_iStart[g],m_ew->m_iEnd[g],m_ew->m_jStart[g],
		    m_ew->m_jEnd[g], m_ew->m_kStart[g], m_ew->m_kEnd[g] );
      mu[g].define(m_ew->m_iStart[g],m_ew->m_iEnd[g],m_ew->m_jStart[g],
		    m_ew->m_jEnd[g], m_ew->m_kStart[g], m_ew->m_kEnd[g] );
      lambda[g].define(m_ew->m_iStart[g],m_ew->m_iEnd[g],m_ew->m_jStart[g],
		    m_ew->m_jEnd[g], m_ew->m_kStart[g], m_ew->m_kEnd[g] );
   }
   get_material( nmd, xmd, nms, xms, rho, mu, lambda );
   int info;
   m_ew->project_material( rho, mu, lambda, info );
   if( info != 0 )
      get_parameters( nmd, xmd, nms, xms, rho, mu, lambda );
}

//-----------------------------------------------------------------------
void MaterialParameterization::write_parameters( const char* filename,
						int nms, double* xms )
{
   VERIFY2( nms == m_nms, "MP::write_parameters: Error in sizes ");
   // Format: nms,xms[0],xms[1],..xms[nms-1]
   if( m_myrank == 0 )
   {
      int fd=open(filename,O_CREAT|O_TRUNC|O_WRONLY,0660);
      size_t nr=write(fd,&nms,sizeof(int));
      if( nr != sizeof(int) )
	 cout << "Error in MaterialParameterization::write_parameters "
	      << " could not write size " << endl;
      nr=write(fd,xms,nms*sizeof(double));
      if( nr != nms*sizeof(double) )
	 cout << "Error in MaterialParameterization::write_parameters "
	      << " could not write parameters " << endl;
      close(fd);
   }
}

//-----------------------------------------------------------------------
void MaterialParameterization::write_parameters( int nms, double* xms )
{
   VERIFY2( nms == m_nms, "MP::write_parameters: Error in sizes ");
   // Format: nms,xms[0],xms[1],..xms[nms-1]
   if( m_myrank == 0 )
   {
      int fd=open(m_filename,O_CREAT|O_TRUNC|O_WRONLY,0660);
      size_t nr=write(fd,&nms,sizeof(int));
      if( nr != sizeof(int) )
	 cout << "Error in MaterialParameterization::write_parameters "
	      << " could not write size " << endl;
      nr=write(fd,xms,nms*sizeof(double));
      if( nr != nms*sizeof(double) )
	 cout << "Error in MaterialParameterization::write_parameters "
	      << " could not write parameters " << endl;
      close(fd);
   }
}

//-----------------------------------------------------------------------
void MaterialParameterization::read_parameters( const char* filename,
						int nms, double* xms )
{
   int errflag = 0;
   if( m_myrank == 0 )
   {
      int fd=open(filename,O_RDONLY );
      int nms_read;
      size_t nr = read(fd,&nms_read,sizeof(int));
      if( nms_read == nms && nr == sizeof(int) && nms_read == m_nms )
      {
	 nr = read(fd,xms,nms*sizeof(double));
         if( nr != nms*sizeof(double) )
	    errflag = 2;
      }
      else if( nms_read != m_nms )
	 errflag = 3;
      else
	 errflag = 1;
      close(fd);
   }
   MPI_Bcast( xms, nms, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   VERIFY2( errflag == 0, "Error no " << errflag << " in MaterialParameterization::read_parameters");
}

//-----------------------------------------------------------------------
void MaterialParameterization::read_parameters( int nms, double* xms )
{
   int errflag = 0;
   if( m_myrank == 0 )
   {
      int fd=open(m_filename,O_RDONLY );
      int nms_read;
      size_t nr = read(fd,&nms_read,sizeof(int));
      if( nms_read == nms && nr == sizeof(int) && nms_read == m_nms )
      {
	 nr = read(fd,xms,nms*sizeof(double));
         if( nr != nms*sizeof(double) )
	    errflag = 2;
      }
      else if( nms_read != m_nms )
	 errflag = 3;
      else
	 errflag = 1;
      close(fd);
   }
   MPI_Bcast( xms, nms, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   VERIFY2( errflag == 0, "Error no " << errflag << " in MaterialParameterization::read_parameters");
}

//-----------------------------------------------------------------------
void MaterialParameterization::parameters_from_basematerial( int nmd, double* xmd,
							 int nms, double* xms )
{
   get_parameters( nmd, xmd, nms, xms, m_ew->mRho, m_ew->mMu, m_ew->mLambda );
}
