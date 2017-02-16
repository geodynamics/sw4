#include "EW.h"

#include "EWCuda.h"

#ifdef SW4_CUDA

// IBM Comment: Uncomment following macro to check kernel error
#define DEBUG_CUDA
  #ifdef DEBUG_CUDA
    #define CHECK_ERROR(str) \
        cudaError err = cudaGetLastError(); \
        if ( cudaSuccess != err ) \
          cerr << "Error in " << str << " : " << cudaGetErrorString( err ) << endl;
  #else
    #define CHECK_ERROR(str)
  #endif

#include <cuda_runtime.h>



void copy_stencilcoefficients( float_sw4* acof, float_sw4* ghcof, float_sw4* bope );

__global__ void pred_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			  float_sw4* up, float_sw4* u, float_sw4* um, float_sw4* lu, float_sw4* fo,
			  float_sw4* rho, float_sw4 dt2, int ghost_points );

__global__ void pred_dev_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			  float_sw4* up, float_sw4* u, float_sw4* um, float_sw4* lu, float_sw4* fo,
			  float_sw4* rho, float_sw4 dt2, int ghost_points );

__global__ void corr_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			  float_sw4* up, float_sw4* lu, float_sw4* fo,
			  float_sw4* rho, float_sw4 dt4, int ghost_points );

__global__ void corr_dev_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			  float_sw4* up, float_sw4* lu, float_sw4* fo,
			  float_sw4* rho, float_sw4 dt4, int ghost_points );

__global__ void dpdmt_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			   float_sw4* up, float_sw4* u, float_sw4* um,
			   float_sw4* u2, float_sw4 dt2i, int ghost_points );

__global__ void addsgd4_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			     float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
			     float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
			     float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
			     float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
			     float_sw4 beta, int ghost_points );

__global__ void addsgd4_dev_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			     float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
			     float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
			     float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
			     float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
			     float_sw4 beta, int ghost_points );

__global__ void addsgd6_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			     float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
			     float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
			     float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
			     float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
			     float_sw4 beta, int ghost_points );

__global__ void addsgd6_dev_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			     float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
			     float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
			     float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
			     float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
			     float_sw4 beta, int ghost_points );

__global__ void rhs4center_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
				float_sw4* a_lu, float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, 
				float_sw4 h, float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
				int ghost_points );

__global__ void rhs4center_dev_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
				float_sw4* a_lu, float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, 
				float_sw4 h, float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
				int ghost_points );

__global__ void rhs4upper_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			       float_sw4* a_lu, float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, 
			       float_sw4 h, float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
			       int ghost_points );

__global__ void rhs4upper_dev_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			       float_sw4* a_lu, float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, 
			       float_sw4 h, float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
			       int ghost_points );

__global__ void rhs4lower_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			       int nk, float_sw4* a_lu, float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, 
			       float_sw4 h, float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
			       int ghost_points );

__global__ void rhs4lower_dev_rev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
				   int nk, float_sw4* a_lu, float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda, 
			       float_sw4 h, float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
			       int ghost_points );

__global__ void check_nan_dev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			       float_sw4* u, int* retval_dev, int* idx_dev );

__global__ void forcing_dev( float_sw4 t, Sarray* dev_F, int NumberOfGrids, GridPointSource** dev_point_sources,
			     int nptsrc, int* dev_identsources, int nident, bool tt );
__global__ void init_forcing_dev( GridPointSource** point_sources, int nsrc );
#endif

//-----------------------------------------------------------------------
void EW::evalRHSCU(vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		   vector<Sarray> & a_Uacc, int st )
{
#ifdef SW4_CUDA
   dim3 gridsize, blocksize;
   gridsize.x  = m_gpu_gridsize[0];
   gridsize.y  = m_gpu_gridsize[1];
   gridsize.z  = m_gpu_gridsize[2];
   blocksize.x = m_gpu_blocksize[0];
   blocksize.y = m_gpu_blocksize[1];
   blocksize.z = m_gpu_blocksize[2];
   for(int g=0 ; g<mNumberOfCartesianGrids; g++ )
   {
      if( m_corder )
	 rhs4center_dev_rev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>
	 ( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
	   m_kEnd[g], a_Uacc[g].dev_ptr(), a_U[g].dev_ptr(), a_Mu[g].dev_ptr(),
	   a_Lambda[g].dev_ptr(), mGridSize[g],
	   dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g], m_ghost_points );
      else
	 rhs4center_dev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>
	 ( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
	   m_kEnd[g], a_Uacc[g].dev_ptr(), a_U[g].dev_ptr(), a_Mu[g].dev_ptr(),
	   a_Lambda[g].dev_ptr(), mGridSize[g],
	   dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g], m_ghost_points );
      CHECK_ERROR("evalRHSCU")
   }
// Boundary operator at upper boundary
   blocksize.z = 1;
   gridsize.z  = 6;

   for(int g=0 ; g<mNumberOfCartesianGrids; g++ )
   {
      if( m_onesided[g][4] )
      {
	 if( m_corder )
	    rhs4upper_dev_rev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>
	    ( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],  m_kEnd[g],
	      a_Uacc[g].dev_ptr(), a_U[g].dev_ptr(), a_Mu[g].dev_ptr(),
	      a_Lambda[g].dev_ptr(), mGridSize[g],
	      dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g], m_ghost_points );
	 else
	    rhs4upper_dev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>
	    ( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],  m_kEnd[g],
	      a_Uacc[g].dev_ptr(), a_U[g].dev_ptr(), a_Mu[g].dev_ptr(),
	      a_Lambda[g].dev_ptr(), mGridSize[g],
	      dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g], m_ghost_points );
      }
      CHECK_ERROR("evalRHSCU_upper")
   }
#endif
}

//-----------------------------------------------------------------------
void EW::evalPredictorCU( vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
			  vector<Sarray>& a_Rho, vector<Sarray> & a_Lu, vector<Sarray> & a_F, int st )
//			  vector<Sarray>& a_Rho, vector<Sarray> & a_Lu, Sarray* a_F, int st )
{
#ifdef SW4_CUDA
   float_sw4 dt2 = mDt*mDt;
   dim3 gridsize, blocksize;
   gridsize.x  = m_gpu_gridsize[0] * m_gpu_gridsize[1] * m_gpu_gridsize[2];
   gridsize.y  = 1;
   gridsize.z  = 1;
   blocksize.x = m_gpu_blocksize[0] * m_gpu_blocksize[1] * m_gpu_blocksize[2];
   blocksize.y = 1;
   blocksize.z = 1;
   for( int g=0 ; g<mNumberOfGrids; g++ )
   {
      if( m_corder )
	 pred_dev_rev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								   m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), a_U[g].dev_ptr(),
								   a_Um[g].dev_ptr(), a_Lu[g].dev_ptr(), a_F[g].dev_ptr(),
								   a_Rho[g].dev_ptr(), dt2, m_ghost_points );
      else
	 pred_dev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								   m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), a_U[g].dev_ptr(),
								   a_Um[g].dev_ptr(), a_Lu[g].dev_ptr(), a_F[g].dev_ptr(),
								   a_Rho[g].dev_ptr(), dt2, m_ghost_points );
      CHECK_ERROR("evalPredictorCU")
   }
#endif
}

//---------------------------------------------------------------------------
void EW::evalCorrectorCU( vector<Sarray> & a_Up, vector<Sarray>& a_Rho,
			  vector<Sarray> & a_Lu, vector<Sarray> & a_F, int st )
//			  vector<Sarray> & a_Lu, Sarray* a_F, int st )
{
#ifdef SW4_CUDA
   float_sw4 dt4 = mDt*mDt*mDt*mDt;  
   dim3 gridsize, blocksize;
   gridsize.x  = m_gpu_gridsize[0] * m_gpu_gridsize[1] * m_gpu_gridsize[2];
   gridsize.y  = 1;
   gridsize.z  = 1;
   blocksize.x = m_gpu_blocksize[0] * m_gpu_blocksize[1] * m_gpu_blocksize[2];
   blocksize.y = 1;
   blocksize.z = 1;
   for( int g=0 ; g<mNumberOfGrids; g++ )
   {
      if( m_corder )
	 corr_dev_rev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), a_Lu[g].dev_ptr(),
								a_F[g].dev_ptr(), a_Rho[g].dev_ptr(), dt4,
								m_ghost_points );
      else
	 corr_dev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), a_Lu[g].dev_ptr(),
								a_F[g].dev_ptr(), a_Rho[g].dev_ptr(), dt4,
								m_ghost_points );
      CHECK_ERROR("evalCorrectorCU")
   }
#endif
}

//---------------------------------------------------------------------------
void EW::evalDpDmInTimeCU(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
			  vector<Sarray> & a_Uacc, int st )
{
#ifdef SW4_CUDA
   float_sw4 dt2i = 1./(mDt*mDt);
   dim3 gridsize, blocksize;
   gridsize.x  = m_gpu_gridsize[0] * m_gpu_gridsize[1] * m_gpu_gridsize[2];
   gridsize.y  = 1;
   gridsize.z  = 1;
   blocksize.x = m_gpu_blocksize[0] * m_gpu_blocksize[1] * m_gpu_blocksize[2];
   blocksize.y = 1;
   blocksize.z = 1;
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
      dpdmt_dev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								 m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), a_U[g].dev_ptr(),
								 a_Um[g].dev_ptr(), a_Uacc[g].dev_ptr(), dt2i,
								 m_ghost_points );
      CHECK_ERROR("evalDpDmInTimeCU")
   }
#endif
}

//-----------------------------------------------------------------------
void EW::addSuperGridDampingCU(vector<Sarray> & a_Up, vector<Sarray> & a_U,
			     vector<Sarray> & a_Um, vector<Sarray> & a_Rho, int st )
{
#ifdef SW4_CUDA
   dim3 gridsize, blocksize;
   gridsize.x  = m_gpu_gridsize[0];
   gridsize.y  = m_gpu_gridsize[1];
   gridsize.z  = m_gpu_gridsize[2];
   blocksize.x = m_gpu_blocksize[0];
   blocksize.y = m_gpu_blocksize[1];
   blocksize.z = m_gpu_blocksize[2];
   for(int g=0 ; g<mNumberOfGrids; g++ )
   {
      if( m_sg_damping_order == 4 )
      {
	 if( m_corder )
	    addsgd4_dev_rev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								      m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), 
                                                                      a_U[g].dev_ptr(), a_Um[g].dev_ptr(), a_Rho[g].dev_ptr(),
								      dev_sg_dc_x[g], dev_sg_dc_y[g], dev_sg_dc_z[g],
								      dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g],
								      dev_sg_corner_x[g], dev_sg_corner_y[g], dev_sg_corner_z[g],
								      m_supergrid_damping_coefficient, m_ghost_points );
	 else
	    addsgd4_dev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								      m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), 
                                                                      a_U[g].dev_ptr(), a_Um[g].dev_ptr(), a_Rho[g].dev_ptr(),
								      dev_sg_dc_x[g], dev_sg_dc_y[g], dev_sg_dc_z[g],
								      dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g],
								      dev_sg_corner_x[g], dev_sg_corner_y[g], dev_sg_corner_z[g],
								      m_supergrid_damping_coefficient, m_ghost_points );
	 CHECK_ERROR("addsgd4_dev")
      }
      else if(  m_sg_damping_order == 6 )
      {
	 if( m_corder )
	    addsgd6_dev_rev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								      m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), 
                                                                      a_U[g].dev_ptr(), a_Um[g].dev_ptr(), a_Rho[g].dev_ptr(),
								      dev_sg_dc_x[g], dev_sg_dc_y[g], dev_sg_dc_z[g],
								      dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g],
								      dev_sg_corner_x[g], dev_sg_corner_y[g], dev_sg_corner_z[g],
								      m_supergrid_damping_coefficient, m_ghost_points );
	 else
	    addsgd6_dev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
								      m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(), 
                                                                      a_U[g].dev_ptr(), a_Um[g].dev_ptr(), a_Rho[g].dev_ptr(),
								      dev_sg_dc_x[g], dev_sg_dc_y[g], dev_sg_dc_z[g],
								      dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g],
								      dev_sg_corner_x[g], dev_sg_corner_y[g], dev_sg_corner_z[g],
								      m_supergrid_damping_coefficient, m_ghost_points );
	 CHECK_ERROR("addsgd6_dev")
      }
   }
#endif
}

//-----------------------------------------------------------------------
void EW::setupSBPCoeff()
{
   //   float_sw4 gh2; // this coefficient is also stored in m_ghcof[0]
   if (mVerbose >=1 && m_myrank == 0)
      cout << "Setting up SBP boundary stencils" << endl;
// get coefficients for difference approximation of 2nd derivative with variable coefficients
//      call VARCOEFFS4( acof, ghcof )
//   F77_FUNC(varcoeffs4,VARCOEFFS4)(m_acof, m_ghcof);
// get coefficients for difference approximation of 1st derivative
//      call WAVEPROPBOP_4( iop, iop2, bop, bop2, gh2, hnorm, sbop )
//   F77_FUNC(wavepropbop_4,WAVEPROPBOP_4)(m_iop, m_iop2, m_bop, m_bop2, &gh2, m_hnorm, m_sbop);
// extend the definition of the 1st derivative tothe first 6 points
//      call BOPEXT4TH( bop, bope )
//   F77_FUNC(bopext4th,BOPEXT4TH)(m_bop, m_bope);
   GetStencilCoefficients( m_acof, m_ghcof, m_bope, m_sbop );
#ifdef SW4_CUDA
   copy_stencilcoefficients( m_acof, m_ghcof, m_bope );
#endif   
}

//-----------------------------------------------------------------------
void EW::copy_supergrid_arrays_to_device()
{
#ifdef SW4_CUDA
  dev_sg_str_x.resize(mNumberOfGrids);
  dev_sg_str_y.resize(mNumberOfGrids);
  dev_sg_str_z.resize(mNumberOfGrids);
  dev_sg_dc_x.resize(mNumberOfGrids);
  dev_sg_dc_y.resize(mNumberOfGrids);
  dev_sg_dc_z.resize(mNumberOfGrids);
  dev_sg_corner_x.resize(mNumberOfGrids);
  dev_sg_corner_y.resize(mNumberOfGrids);
  dev_sg_corner_z.resize(mNumberOfGrids);
  if( m_ndevice > 0 )
  {
     cudaError_t retcode;
     for( int g=0 ; g<mNumberOfGrids; g++) 
     {
	// sg_str
	retcode = cudaMalloc( (void**)&dev_sg_str_x[g], sizeof(float_sw4)*(m_iEnd[g]-m_iStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc x returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMalloc( (void**)&dev_sg_str_y[g], sizeof(float_sw4)*(m_jEnd[g]-m_jStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc y returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMalloc( (void**)&dev_sg_str_z[g], sizeof(float_sw4)*(m_kEnd[g]-m_kStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc z returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_str_x[g], m_sg_str_x[g], sizeof(float_sw4)*(m_iEnd[g]-m_iStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy x returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_str_y[g], m_sg_str_y[g], sizeof(float_sw4)*(m_jEnd[g]-m_jStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy y returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_str_z[g], m_sg_str_z[g], sizeof(float_sw4)*(m_kEnd[g]-m_kStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy z returned "
		<< cudaGetErrorString(retcode) << endl;

	// sg_dc
	retcode = cudaMalloc( (void**)&dev_sg_dc_x[g], sizeof(float_sw4)*(m_iEnd[g]-m_iStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc dc_x returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMalloc( (void**)&dev_sg_dc_y[g], sizeof(float_sw4)*(m_jEnd[g]-m_jStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc dc_y returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMalloc( (void**)&dev_sg_dc_z[g], sizeof(float_sw4)*(m_kEnd[g]-m_kStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc dc_z returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_dc_x[g], m_sg_dc_x[g], sizeof(float_sw4)*(m_iEnd[g]-m_iStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy dc_x returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_dc_y[g], m_sg_dc_y[g], sizeof(float_sw4)*(m_jEnd[g]-m_jStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy dc_y returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_dc_z[g], m_sg_dc_z[g], sizeof(float_sw4)*(m_kEnd[g]-m_kStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy dc_z returned "
		<< cudaGetErrorString(retcode) << endl;
	// sg_corner
	retcode = cudaMalloc( (void**)&dev_sg_corner_x[g], sizeof(float_sw4)*(m_iEnd[g]-m_iStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc corner_x returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMalloc( (void**)&dev_sg_corner_y[g], sizeof(float_sw4)*(m_jEnd[g]-m_jStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc corner_y returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMalloc( (void**)&dev_sg_corner_z[g], sizeof(float_sw4)*(m_kEnd[g]-m_kStart[g]+1));
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMalloc corner_z returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_corner_x[g], m_sg_corner_x[g], sizeof(float_sw4)*(m_iEnd[g]-m_iStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy corner_x returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_corner_y[g], m_sg_corner_y[g], sizeof(float_sw4)*(m_jEnd[g]-m_jStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy corner_y returned "
		<< cudaGetErrorString(retcode) << endl;
	retcode = cudaMemcpy( dev_sg_corner_z[g], m_sg_corner_z[g], sizeof(float_sw4)*(m_kEnd[g]-m_kStart[g]+1),
			      cudaMemcpyHostToDevice );
	if( retcode != cudaSuccess )
	   cout << "Error, EW::copy_supergrid_arrays_to_device cudaMemcpy corner_z returned "
		<< cudaGetErrorString(retcode) << endl;
     }
  }
#endif   
}

//-----------------------------------------------------------------------
void EW::copy_material_to_device()
{
#ifdef SW4_CUDA
   for( int g=0 ; g<mNumberOfGrids; g++)
   {
      mMu[g].copy_to_device( m_cuobj );
      mLambda[g].copy_to_device( m_cuobj );
      mRho[g].copy_to_device( m_cuobj );
   }
#endif   
}

//-----------------------------------------------------------------------
void EW::find_cuda_device()
{
#ifdef SW4_CUDA
   cudaError_t retcode;
   cudaDeviceProp prop;
   retcode = cudaGetDeviceCount(&m_ndevice);
   if( retcode != cudaSuccess )
   {
      cout << "Error from cudaGetDeviceCount: Error string = " <<
	 cudaGetErrorString(retcode) << endl;
   }
   // Note: This will not work if some nodes have GPU and others do not
   // It is assumed that all nodes are identical wrt GPU
   if( m_ndevice > 0 && m_myrank == 0 )
   {
      cout << m_ndevice << " Cuda devices found:" << endl;
      for( int i=0 ;  i < m_ndevice ; i++ )
      {
	 retcode = cudaGetDeviceProperties( &prop, i );
	 cout << "      Device " << i << ": name = " << prop.name <<
	    ",  Compute capability:" << prop.major << "." << prop.minor << 
	    ",  Memory (GB) " << (prop.totalGlobalMem  >> 30) << endl;
      }
   }
   //Added following line for all ranks 
   retcode = cudaGetDeviceProperties(&prop, 0 );

   // Check block size
   CHECK_INPUT( m_gpu_blocksize[0] <= prop.maxThreadsDim[0],
		"Error: EW::find_cuda_device max block x " << m_gpu_blocksize[0] << " too large\n");
   CHECK_INPUT( m_gpu_blocksize[1] <= prop.maxThreadsDim[1], 
		"Error: EW::find_cuda_device max block y " << m_gpu_blocksize[1] << " too large\n");
   CHECK_INPUT( m_gpu_blocksize[2] <= prop.maxThreadsDim[2],
		"Error: EW::find_cuda_device max block z " << m_gpu_blocksize[2] << " too large\n");
   CHECK_INPUT( m_gpu_blocksize[0]*m_gpu_blocksize[1]*m_gpu_blocksize[2] <= prop.maxThreadsPerBlock, 
   "Error: EW::find_cuda_device max number of threads per block " << prop.maxThreadsPerBlock <<
		" is exceeded \n");
   // Determine grid dimensions.
   int ghost = m_ghost_points;
   int ni=m_iEnd[0]-m_iStart[0]+1-2*ghost;
   int nj=m_jEnd[0]-m_jStart[0]+1-2*ghost;
   int nk=m_kEnd[0]-m_kStart[0]+1-2*ghost;
   bool m_gpu_overlap = false;
   if( m_gpu_overlap )
   {
      REQUIRE2( m_gpu_blocksize[0] > 2*ghost && m_gpu_blocksize[1] > 2*ghost 
	       && m_gpu_blocksize[2] > 2*ghost , "Error, need block size at least "
	       << 2*ghost+1 << " in each direction\n" );
      m_gpu_gridsize[0]=ni/(m_gpu_blocksize[0]-2*ghost);
      if( ni % (m_gpu_blocksize[0]-2*ghost)  != 0 )
	 m_gpu_gridsize[0]++;
      m_gpu_gridsize[1]=nj/(m_gpu_blocksize[1]-2*ghost);
      if( nj % (m_gpu_blocksize[1]-2*ghost)  != 0 )
	 m_gpu_gridsize[1]++;
      m_gpu_gridsize[2]=nk/(m_gpu_blocksize[2]-2*ghost);
      if( nk % (m_gpu_blocksize[2]-2*ghost)  != 0 )
	 m_gpu_gridsize[2]++;
   }
   else
   {
      m_gpu_gridsize[0]=ni/m_gpu_blocksize[0];
      if( ni % m_gpu_blocksize[0]  != 0 )
	 m_gpu_gridsize[0]++;
      m_gpu_gridsize[1]=nj/m_gpu_blocksize[1];
      if( nj % m_gpu_blocksize[1]  != 0 )
	 m_gpu_gridsize[1]++;
      m_gpu_gridsize[2]=nk/m_gpu_blocksize[2];
      if( nk % m_gpu_blocksize[2]  != 0 )
	 m_gpu_gridsize[2]++;
   }
   if( m_myrank == 0 )
   {
   cout << " size of domain " << ni << " x " << nj << " x " << nk << " grid points (excluding ghost points)\n";
   cout << " GPU block size " <<  m_gpu_blocksize[0] << " x " <<  m_gpu_blocksize[1] << " x " <<  m_gpu_blocksize[2] << "\n";
   cout << " GPU grid size  " <<  m_gpu_gridsize[0]  << " x " <<  m_gpu_gridsize[1]  << " x " <<  m_gpu_gridsize[2] << endl;
   cout << " GPU cores " << m_gpu_blocksize[0]*m_gpu_gridsize[0] << " x " 
	   << m_gpu_blocksize[1]*m_gpu_gridsize[1] << " x "
	   << m_gpu_blocksize[2]*m_gpu_gridsize[2] << endl;
   }
   if( m_ndevice > 0 )
   {
      // create two streams
      m_cuobj = new EWCuda( m_ndevice, 2 );
   }
#endif
   if( m_ndevice == 0 && m_myrank == 0 )
      cout << " no GPUs found" << endl;
   if( m_ndevice == 0 )
      m_cuobj = new EWCuda(0,0);
}

//-----------------------------------------------------------------------
bool EW::check_for_nan_GPU( vector<Sarray>& a_U, int verbose, string name )
{
#ifdef SW4_CUDA
   dim3 gridsize, blocksize;
   gridsize.x  = m_gpu_gridsize[0] * m_gpu_gridsize[1] * m_gpu_gridsize[2];
   gridsize.y  = 1;
   gridsize.z  = 1;
   blocksize.x = m_gpu_blocksize[0] * m_gpu_blocksize[1] * m_gpu_blocksize[2];
   blocksize.y = 1;
   blocksize.z = 1;

   int  *retval_dev, retval_host = 0;
   int  *idx_dev, idx_host = 0;
   int cnan, inan, jnan, knan;
   int  nijk, nij, ni;

   cudaMalloc( (void **)&retval_dev, sizeof(int) );
   cudaMemcpy( retval_dev, &retval_host, sizeof(int), cudaMemcpyHostToDevice );
   cudaMalloc( (void **)&idx_dev, sizeof (int) );
   cudaMemcpy( idx_dev, &idx_host, sizeof(int), cudaMemcpyHostToDevice );

   for( int g=0 ; g<mNumberOfGrids; g++ )
   {
      check_nan_dev<<<gridsize,blocksize>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
					     m_kStart[g], m_kEnd[g], a_U[g].dev_ptr(), retval_dev, idx_dev );

      CHECK_ERROR("check_for_nan_GPU")

      cudaMemcpy( &retval_host, retval_dev, sizeof(int), cudaMemcpyDeviceToHost );

      if ( retval_host != 0) 
      {
         cudaMemcpy(&idx_host, idx_dev, sizeof(int), cudaMemcpyDeviceToHost);

         nijk = (m_kEnd[g]-m_kStart[g]+1)*(m_jEnd[g]-m_jStart[g]+1)*(m_iEnd[g]-m_iStart[g]+1);
         nij  = (m_jEnd[g]-m_jStart[g]+1)*(m_iEnd[g]-m_iStart[g]+1);
         ni   = m_iEnd[g]-m_iStart[g]+1;

         cnan = idx_host/nijk;
         idx_host = idx_host - cnan*nijk; 
         knan = idx_host / nij; 
         idx_host = idx_host - knan*nij;
         jnan = idx_host/ni;
         inan = idx_host - jnan*ni; 

         cout << "grid " << g << " array " << name << " found " << retval_host << " nans. First nan at " <<
                    cnan << " " << inan << " " << jnan << " " << knan << endl;

         return false; 
      }
   }

  return true;
#endif
}

//-----------------------------------------------------------------------
void EW::ForceCU( float_sw4 t, Sarray* dev_F, bool tt, int st )
{
#ifdef SW4_CUDA
   dim3 blocksize, gridsize;
   gridsize.x  = m_gpu_gridsize[0] * m_gpu_gridsize[1] * m_gpu_gridsize[2];
   gridsize.y  = 1;
   gridsize.z  = 1;
   blocksize.x = m_gpu_blocksize[0] * m_gpu_blocksize[1] * m_gpu_blocksize[2];
   blocksize.y = 1;
   blocksize.z = 1;
   forcing_dev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>( t, dev_F, mNumberOfGrids, dev_point_sources, m_point_sources.size(),
					dev_identsources, m_identsources.size(), tt );
   CHECK_ERROR("ForceCU")
#endif
}

//-----------------------------------------------------------------------
void EW::init_point_sourcesCU( )
{
#ifdef SW4_CUDA
   dim3 blocksize, gridsize;
   gridsize.x  = m_gpu_gridsize[0] * m_gpu_gridsize[1] * m_gpu_gridsize[2];
   gridsize.y  = 1;
   gridsize.z  = 1;
   blocksize.x = m_gpu_blocksize[0] * m_gpu_blocksize[1] * m_gpu_blocksize[2];
   blocksize.y = 1;
   blocksize.z = 1;
   init_forcing_dev<<<gridsize,blocksize>>>( dev_point_sources, m_point_sources.size() );
   CHECK_ERROR("InitPointSourcesCU")
#endif
}
