#include <mpi.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstring>
#include <fstream>

#include "MaterialParameterization.h"

#include "EW.h"
#include "Require.h"

using namespace std;

//-----------------------------------------------------------------------
MaterialParameterization::MaterialParameterization( EW* a_ew, char* fname )
{
   m_ew = a_ew;
   MPI_Comm_rank( m_ew->m_1d_communicator, &m_myrank );
   int n = strlen(fname)+1;
   m_filename = new char[n];
   strcpy(m_filename,fname);
   m_path = "./";
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
//void MaterialParameterization::store_material( vector<Sarray>& a_rho,
//					       vector<Sarray>& a_mu,
//					       vector<Sarray>& a_lambda )
//{
//   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
//   {
//      m_ew->mRho[g].copy( a_rho[g] );
//      m_ew->mMu[g].copy( a_mu[g] );
//      m_ew->mLambda[g].copy( a_lambda[g] );
//   }
//}

//-----------------------------------------------------------------------
//void MaterialParameterization::constrain_material( int nmd, double* xmd,
//						   int nms, double* xms )
//{
//   VERIFY2( nmd== m_nmd && nms == m_nms, "MP::constrain_material: Error in sizes");
//   int ng = m_ew->mNumberOfGrids;
//   vector<Sarray> rho, mu, lambda;
//   rho.resize(ng);
//   mu.resize(ng);
//   lambda.resize(ng);
//   for( int g=0; g < ng ; g++ )
//   {
//      rho[g].define(m_ew->m_iStart[g],m_ew->m_iEnd[g],m_ew->m_jStart[g],
//		    m_ew->m_jEnd[g], m_ew->m_kStart[g], m_ew->m_kEnd[g] );
//      mu[g].define(m_ew->m_iStart[g],m_ew->m_iEnd[g],m_ew->m_jStart[g],
//		    m_ew->m_jEnd[g], m_ew->m_kStart[g], m_ew->m_kEnd[g] );
//      lambda[g].define(m_ew->m_iStart[g],m_ew->m_iEnd[g],m_ew->m_jStart[g],
//		    m_ew->m_jEnd[g], m_ew->m_kStart[g], m_ew->m_kEnd[g] );
//   }
//   get_material( nmd, xmd, nms, xms, rho, mu, lambda );
//   int info;
//   m_ew->project_material( rho, mu, lambda, info );
//   if( info != 0 )
//      get_parameters( nmd, xmd, nms, xms, rho, mu, lambda );
//}

//-----------------------------------------------------------------------
void MaterialParameterization::write_parameters( const char* filename,
						int nms, double* xms )
{
  VERIFY2( nms == m_nms, "MP::write_parameters: Error in sizes nms = " << nms << " m_nms = " << m_nms );
   // Format: nms,xms[0],xms[1],..xms[nms-1]
   if( m_myrank == 0 && nms > 0 )
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
      bool write_ascii=true;
      if( write_ascii )
      {
 	 fstream fut("parameters.txt");
	 for( int i=0 ; i < nms ; i++ )
	   fut << i+1 << " " << xms[i] << endl;
	 fut.close();
      }
   }
}

//-----------------------------------------------------------------------
void MaterialParameterization::write_parameters( int nms, double* xms )
{
   VERIFY2( nms == m_nms, "MP::write_parameters: Error in sizes ");
   // Format: nms,xms[0],xms[1],..xms[nms-1]
   if( m_myrank == 0 && nms > 0 )
   {
      string fname = m_path + m_filename;
      int fd=open(fname.c_str(),O_CREAT|O_TRUNC|O_WRONLY,0660);
      //      int fd=open(m_filename,O_CREAT|O_TRUNC|O_WRONLY,0660);
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
void MaterialParameterization::write_parameters_dist( const char* outfile,
                                                      int nmd, double* xmd )
{
   if( nmd > 0 )
   {
      std::stringstream fid;
      fid << m_myrank << "\0";
      string fname = m_path + outfile + fid.str(); 
      int fd=open(fname.c_str(),O_CREAT|O_TRUNC|O_WRONLY,0660);
      //      int fd=open(m_filename,O_CREAT|O_TRUNC|O_WRONLY,0660);
      size_t nr=write(fd,&nmd,sizeof(int));
      if( nr != sizeof(int) )
	 cout << "Error in MaterialParameterization::write_parameters_dist "
	      << " could not write size " << endl;
      nr=write(fd,xmd,nmd*sizeof(double));
      if( nr != nmd*sizeof(double) )
	 cout << "Error in MaterialParameterization::write_parameters_dist "
	      << " could not write parameters " << endl;
      close(fd);
   }

}

//-----------------------------------------------------------------------
void MaterialParameterization::read_parameters( const char* filename,
						int npars, double* xptr )
{
   // Assumes global distribution of xptr. 
   // Read from one processor, and broadcast to all.
   int errflag = 0;
   if( m_myrank == 0 )
   {
      int fd=open(filename,O_RDONLY );
      int npars_read;
      size_t nr = read(fd,&npars_read,sizeof(int));
      if( npars_read == npars && nr == sizeof(int) )
      {
	 nr = read(fd,xptr,npars*sizeof(double));
         if( nr != npars*sizeof(double) )
	    errflag = 2;
      }
      else
	 errflag = 1;
      close(fd);
   }
   MPI_Bcast( xptr, npars, MPI_DOUBLE, 0, m_ew->m_1d_communicator );
   VERIFY2( errflag == 0, "Error no " << errflag << " in MaterialParameterization::read_parameters");
}

//-----------------------------------------------------------------------
void MaterialParameterization::read_parameters( int npars, double* xptr )
{
   read_parameters( m_filename, npars, xptr );
}

////-----------------------------------------------------------------------
//void MaterialParameterization::read_parameters( int nms, double* xms )
//{
//   int errflag = 0;
//   if( m_myrank == 0 )
//   {
//      //      string fname = m_path + m_filename;
//      string fname = m_filename;
//      int fd=open(fname.c_str(),O_RDONLY );
//      //      int fd=open(m_filename,O_RDONLY );
//      //      VERIFY2( fd != -1, "Error opening file " << m_filename << " in MaterialParameterization::read_parameters"<<endl);
//      VERIFY2( fd != -1, "Error opening file " << fname << " in MaterialParameterization::read_parameters"<<endl);
//      int nms_read;
//      size_t nr = read(fd,&nms_read,sizeof(int));
//      if( nms_read == nms && nr == sizeof(int) && nms_read == m_nms )
//      {
//	 nr = read(fd,xms,nms*sizeof(double));
//         if( nr != nms*sizeof(double) )
//	    errflag = 2;
//      }
//      else if( nms_read != m_nms )
//      {
//	 errflag = 3;
//      std::cout << "nms_read = " << nms_read << " m_nms= " << m_nms << " nms= " << nms << std::endl;
//      }
//      else
//	 errflag = 1;
//      close(fd);
//   }
//   MPI_Bcast( xms, nms, MPI_DOUBLE, 0, MPI_COMM_WORLD );
//   VERIFY2( errflag == 0, "Error no " << errflag << " in MaterialParameterization::read_parameters"<<endl);
//}

//-----------------------------------------------------------------------
//void MaterialParameterization::parameters_from_basematerial( int nmd, double* xmd,
//							 int nms, double* xms )
//{
//   get_parameters( nmd, xmd, nms, xms, m_ew->mRho, m_ew->mMu, m_ew->mLambda );
//}

//-----------------------------------------------------------------------
void MaterialParameterization::set_scalefactors( int nmpars, double* sfs, 
                  double rho_ref, double mu_ref, double lambda_ref, double vs_ref, double vp_ref )
{
   for( int i=0 ; i < nmpars ; i += 3 )
   {
      sfs[i]   = rho_ref;
      sfs[i+1] = mu_ref;
      sfs[i+2] = lambda_ref;
   }
}

//-----------------------------------------------------------------------
void  MaterialParameterization::get_regularizer( int nmd, double* xmd, int nms, double* xms, 
		      double* xmd0, double* xms0, double regcoeff,
		      std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
		      std::vector<Sarray>& a_lambda, double& mf_reg,
		      double* sfd, double* sfs, bool compute_derivative, 
		      double* dmfd_reg, double* dmfs_reg )
{
   mf_reg = 0;
   if( fabs(regcoeff) < 1e-12 )
      return;

// Default, Tikhonov regularizing term:
#define SQR(x) ((x)*(x))
   mf_reg=0;
   double tcoff = (1.0/(m_nms+m_nmd_global))*regcoeff;
   if( compute_derivative )
   {
 // Shared parameters
      double tikhonov=0;
      //      cout << "In get_regularizer " << m_nms << " " << m_nmd << endl;
      for (int q=0; q<m_nms; q++)
      {
	 tikhonov +=  SQR( (xms[q] - xms0[q])/sfs[q]);
	 dmfs_reg[q]   = 2*tcoff*(xms[q] - xms0[q])/SQR(sfs[q]);
         //	 tikhonov +=  SQR( (xms[q] - xms0[q]));
         //	 dmfs_reg[q]   = 2*tcoff*(xms[q] - xms0[q]);
	 //         cout << " q = " << q << " sfs = " << sfs[q] << " " << xms0[q] << " " << xms[q] << endl;
      }
      mf_reg += tcoff*tikhonov;

 // Distributed parameters
      double tikhonovd = 0;
      for (int q=0; q<m_nmd; q++)
      {
	 tikhonovd += SQR( (xmd[q] - xmd0[q])/sfd[q]);
	 dmfd_reg[q] = 2*tcoff*(xmd[q] - xmd0[q])/SQR(sfd[q]);
         //	 tikhonovd += SQR( (xmd[q] - xmd0[q]));
         //	 dmfd_reg[q] = 2*tcoff*(xmd[q] - xmd0[q]);
      }
      MPI_Allreduce( &tikhonovd, &tikhonov, 1, MPI_DOUBLE, MPI_SUM, m_ew->m_1d_communicator );
      mf_reg += tcoff*tikhonov;
   }
   else
   {
      double tikhonov=0;
      for (int q=0; q<m_nms; q++)
      {
         //	 tikhonov +=  SQR( (xms[q] - xms0[q])/sfs[q]);
	 tikhonov +=  SQR( (xms[q] - xms0[q]));
      }
      mf_reg += tcoff*tikhonov;

      double tikhonovd = 0;
      for (int q=0; q<m_nmd; q++)
	 tikhonovd += SQR( (xmd[q] - xmd0[q])/sfd[q]);
      //	 tikhonovd += SQR( (xmd[q] - xmd0[q]));
      MPI_Allreduce( &tikhonovd, &tikhonov, 1, MPI_DOUBLE, MPI_SUM, m_ew->m_1d_communicator );
      mf_reg += tcoff*tikhonov;
   }
#undef SQR
}
//-----------------------------------------------------------------------
void MaterialParameterization::interpolate_to_cartesian( int nmd, double* xmd, 
                                                         int nms, double* xms,
                                                         std::vector<Sarray>& a_rho,
                                                         std::vector<Sarray>& a_mu,
                                                         std::vector<Sarray>& a_lambda )
{
   get_parameters( nmd, xmd, nms, xms, a_rho, a_mu, a_lambda, 5 );
}

