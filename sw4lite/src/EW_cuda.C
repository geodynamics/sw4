#include "EW.h"
   
#include "EWCuda.h"

#ifdef SW4_CUDA

// IBM Comment: Uncomment following macro to check kernel error
#define DEBUG_CUDA
  #ifdef DEBUG_CUDA
    #define CHECK_ERROR(str) \
        cudaError err = cudaGetLastError(); \
        if ( cudaSuccess != err ) {                                     \
                                   cerr << "Error in " << str << " : " << cudaGetErrorString( err ) << __FILE__ << __LINE__ << endl; \
          exit(0); }
  #else
    #define CHECK_ERROR(str)
  #endif

#include <cuda_runtime.h>



//void copy_stencilcoefficients( float_sw4* acof, float_sw4* ghcof, float_sw4* bope );
void copy_stencilcoefficients1( float_sw4* acof, float_sw4* ghcof, float_sw4* bope, float_sw4*  );



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

__global__ void addsgd4_dev_v2( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
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

__global__ void addsgd4_dev_rev_v2( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
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

__global__ void rhs4center_dev_v2( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                                       float_sw4* a_lu, float_sw4* a_u, float_sw4* a_mu, float_sw4* a_lambda,
                                       float_sw4 h, float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
                                       int ghost_points );

__global__ void rhs4center_dev_rev_v2( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
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

__global__ void BufferToHaloKernel_dev(float_sw4* block_left, float_sw4* block_right, float_sw4* block_up, float_sw4* block_down,
                        float_sw4 * leftSideEdge, float_sw4 * rightSideEdge, float_sw4 * upSideEdge, float_sw4 * downSideEdge,
                        int ni, int nj, int nk, int m_padding, const int m_neighbor0 ,const int  m_neighbor1, const int m_neighbor2,
                        const int m_neighbor3, const int mpi_process_null_cuda);

__global__ void BufferToHaloKernel_dev_rev(float_sw4* block_left, float_sw4* block_right, float_sw4* block_up, float_sw4* block_down,
                        float_sw4 * leftSideEdge, float_sw4 * rightSideEdge, float_sw4 * upSideEdge, float_sw4 * downSideEdge,
                        int ni, int nj, int nk, int m_padding, const int m_neighbor0 ,const int  m_neighbor1, const int m_neighbor2,
                        const int m_neighbor3, const int mpi_process_null_cuda);

__global__ void HaloToBufferKernel_dev(float_sw4* block_left, float_sw4* block_right, float_sw4* block_up, float_sw4* block_down,
                        float_sw4 * leftSideEdge, float_sw4 * rightSideEdge, float_sw4 * upSideEdge, float_sw4 * downSideEdge,
                        int ni, int nj, int nk, int m_padding, const int m_neighbor0 ,const int  m_neighbor1, const int m_neighbor2,
                        const int m_neighbor3, const int mpi_process_null_cuda);

__global__ void HaloToBufferKernel_dev_rev(float_sw4* block_left, float_sw4* block_right, float_sw4* block_up, float_sw4* block_down,
                        float_sw4 * leftSideEdge, float_sw4 * rightSideEdge, float_sw4 * upSideEdge, float_sw4 * downSideEdge,
                        int ni, int nj, int nk, int m_padding, const int m_neighbor0 ,const int  m_neighbor1, const int m_neighbor2,
                        const int m_neighbor3, const int mpi_process_null_cuda);

__global__ void bcfortsg_dev( int ib, int ie, int jb, int je, int kb, int ke, int* wind,
                              int nx, int ny, int nz, float_sw4* u, float_sw4 h, boundaryConditionType *bccnd,
                              float_sw4* mu, float_sw4* la, float_sw4 t,
                              float_sw4* bforce1, float_sw4* bforce2, float_sw4* bforce3,
                              float_sw4* bforce4, float_sw4* bforce5, float_sw4* bforce6,
                              float_sw4 om, float_sw4 ph, float_sw4 cv,
                              float_sw4* strx, float_sw4* stry );


__global__ void bcfortsg_dev_indrev( int ib, int ie, int jb, int je, int kb, int ke, int* wind,
                                     int nx, int ny, int nz, float_sw4* u, float_sw4 h, boundaryConditionType *bccnd,
                                     float_sw4* mu, float_sw4* la, float_sw4 t,
                                     float_sw4* bforce1, float_sw4* bforce2, float_sw4* bforce3,
                                     float_sw4* bforce4, float_sw4* bforce5, float_sw4* bforce6,
                                     float_sw4 om, float_sw4 ph, float_sw4 cv,
                                     float_sw4* strx, float_sw4* stry );

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
      {
        int ni = m_iEnd[g] - m_iStart[g] + 1 - 2*m_ghost_points;
        int nj = m_jEnd[g] - m_jStart[g] + 1 - 2*m_ghost_points;
        int nk = m_kEnd[g] - m_kStart[g] + 1 - 2*m_ghost_points;
        dim3 block(RHS4_BLOCKX,RHS4_BLOCKY);
        dim3 grid;
        grid.x = (ni + block.x - 1) / block.x;
        grid.y = (nj + block.y - 1) / block.y;
        grid.z = 1;
        rhs4center_dev_rev_v2<<<grid,block,0,m_cuobj->m_stream[st]>>>
          ( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
            m_kEnd[g], a_Uacc[g].dev_ptr(), a_U[g].dev_ptr(), a_Mu[g].dev_ptr(),
            a_Lambda[g].dev_ptr(), mGridSize[g],
            dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g], m_ghost_points );
	//	rhs4center_dev_rev<<<gridsize,blocksize,0,m_cuobj->m_stream[st]>>>
	//	        ( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
	//	        m_kEnd[g], a_Uacc[g].dev_ptr(), a_U[g].dev_ptr(), a_Mu[g].dev_ptr(),
	//	          a_Lambda[g].dev_ptr(), mGridSize[g],
	//	          dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g], m_ghost_points );
      }
      else
      {
        int ni = m_iEnd[g] - m_iStart[g] + 1 - 2*m_ghost_points;
        int nj = m_jEnd[g] - m_jStart[g] + 1 - 2*m_ghost_points;
        int nk = m_kEnd[g] - m_kStart[g] + 1 - 2*m_ghost_points;
        dim3 block(RHS4_BLOCKX,RHS4_BLOCKY);
        dim3 grid;
        grid.x = (ni + block.x - 1) / block.x;
        grid.y = (nj + block.y - 1) / block.y;
        grid.z = 1;
        rhs4center_dev_v2<<<grid,block,0,m_cuobj->m_stream[st]>>>
	 ( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
	   m_kEnd[g], a_Uacc[g].dev_ptr(), a_U[g].dev_ptr(), a_Mu[g].dev_ptr(),
	   a_Lambda[g].dev_ptr(), mGridSize[g],
	   dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g], m_ghost_points );
      }
      CHECK_ERROR("rhs4center_dev")
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
         {
           int ni = m_iEnd[g] - m_iStart[g] + 1 - 2*m_ghost_points;
           int nj = m_jEnd[g] - m_jStart[g] + 1 - 2*m_ghost_points;
           int nk = m_kEnd[g] - m_kStart[g] + 1 - 2*m_ghost_points;
           const dim3 block(ADDSGD4_BLOCKX,ADDSGD4_BLOCKY,1);
           dim3 grid;
           grid.x = (ni + block.x - 1) / block.x;
           grid.y = (nj + block.y - 1) / block.y;
           grid.z = 1; 
           addsgd4_dev_rev_v2<<<grid,block,0,m_cuobj->m_stream[st]>>>(
             m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
             m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(),
             a_U[g].dev_ptr(), a_Um[g].dev_ptr(), a_Rho[g].dev_ptr(),
             dev_sg_dc_x[g], dev_sg_dc_y[g], dev_sg_dc_z[g],
             dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g],
             dev_sg_corner_x[g], dev_sg_corner_y[g], dev_sg_corner_z[g],
             m_supergrid_damping_coefficient, m_ghost_points );
           CHECK_ERROR("addsgd4_dev_rev");
         }
	 else
         {
           int ni = m_iEnd[g] - m_iStart[g] + 1 - 2*m_ghost_points;
           int nj = m_jEnd[g] - m_jStart[g] + 1 - 2*m_ghost_points;
           int nk = m_kEnd[g] - m_kStart[g] + 1 - 2*m_ghost_points;
           const dim3 block(ADDSGD4_BLOCKX,ADDSGD4_BLOCKY,1);
           dim3 grid;
           grid.x = (ni + block.x - 1) / block.x;
           grid.y = (nj + block.y - 1) / block.y;
           grid.z = 1; 
           addsgd4_dev_v2<<<grid,block,0,m_cuobj->m_stream[st]>>>(
             m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
             m_kStart[g], m_kEnd[g], a_Up[g].dev_ptr(),
             a_U[g].dev_ptr(), a_Um[g].dev_ptr(), a_Rho[g].dev_ptr(),
             dev_sg_dc_x[g], dev_sg_dc_y[g], dev_sg_dc_z[g],
             dev_sg_str_x[g], dev_sg_str_y[g], dev_sg_str_z[g],
             dev_sg_corner_x[g], dev_sg_corner_y[g], dev_sg_corner_z[g],
             m_supergrid_damping_coefficient, m_ghost_points );
           CHECK_ERROR("addsgd4_dev_rev");
          }
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
   //copy_stencilcoefficients( m_acof, m_ghcof, m_bope );
   copy_stencilcoefficients1( m_acof, m_ghcof, m_bope, m_sbop );
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

//-----------------------------------------------------------------------

void EW::setup_device_communication_array()
{
#ifdef SW4_CUDA
  dev_SideEdge_Send.resize(mNumberOfGrids);
  dev_SideEdge_Recv.resize(mNumberOfGrids);
#ifndef SW4_CUDA_AWARE_MPI
  m_SideEdge_Send.resize(mNumberOfGrids);
  m_SideEdge_Recv.resize(mNumberOfGrids);
#endif

  if( m_ndevice > 0 )
  {
     cudaError_t retcode;

     for( int g=0 ; g<mNumberOfGrids; g++)
     {

        int ni = m_iEnd[g] - m_iStart[g] + 1;
        int nj = m_jEnd[g] - m_jStart[g] + 1;
        int nk = m_kEnd[g] - m_kStart[g] + 1;
        int n_m_ppadding1 = 3*nj*nk*m_ppadding;
        int n_m_ppadding2 = 3*ni*nk*m_ppadding;

        retcode = cudaMalloc( (void**)&dev_SideEdge_Send[g], sizeof(float_sw4)*2*(n_m_ppadding1+n_m_ppadding2) );
        if( retcode != cudaSuccess )
           cout << "Error, EW::setup_device_communication_arra cudaMalloc returned "
                << cudaGetErrorString(retcode) << endl;

        retcode = cudaMalloc( (void**)&dev_SideEdge_Recv[g], sizeof(float_sw4)*2*(n_m_ppadding1+n_m_ppadding2) );
        if( retcode != cudaSuccess )
           cout << "Error, EW::setup_device_communication_arra cudaMalloc returned "
                << cudaGetErrorString(retcode) << endl;

        //retcode = cudaMemset(dev_SideEdge_Send[g], 0.0, sizeof(float_sw4)*2*(n_m_ppadding1+n_m_ppadding2) );
        //if( retcode != cudaSuccess )
        //   cout << "Error, EW::setup_device_communication_arra cudaMalloc returned "
        //        << cudaGetErrorString(retcode) << endl;

        //retcode = cudaMemset(dev_SideEdge_Recv[g], 0.0, sizeof(float_sw4)*2*(n_m_ppadding1+n_m_ppadding2) );
        //if( retcode != cudaSuccess )
        //   cout << "Error, EW::setup_device_communication_arra cudaMalloc returned "
        //        << cudaGetErrorString(retcode) << endl;


#ifndef SW4_CUDA_AWARE_MPI

        retcode = cudaMallocHost( (void**)&m_SideEdge_Send[g], sizeof(float_sw4)*2*(n_m_ppadding1+n_m_ppadding2) );
        if( retcode != cudaSuccess )
           cout << "Error, EW::setup_device_communication_arra cudaMallocHost returned "
                << cudaGetErrorString(retcode) << endl;

        retcode = cudaMallocHost( (void**)&m_SideEdge_Recv[g], sizeof(float_sw4)*2*(n_m_ppadding1+n_m_ppadding2) );
        if( retcode != cudaSuccess )
           cout << "Error, EW::setup_device_communication_arra cudaMallocHost returned "
                << cudaGetErrorString(retcode) << endl;

        //memset( m_SideEdge_Send[g], 0.0, sizeof(float_sw4)*2*(n_m_ppadding1+n_m_ppadding2) );

        //memset( m_SideEdge_Recv[g], 0.0, sizeof(float_sw4)*2*(n_m_ppadding1+n_m_ppadding2) );

#endif

     }
  }
#endif
}

//-----------------------------------------------------------------------

void EW::communicate_arrayCU( Sarray& u, int g , int st)
{
#ifdef SW4_CUDA
   REQUIRE2( u.m_nc == 3 || u.m_nc == 1, "Communicate array, only implemented for one- and three-component arrays"
             << " nc = " << u.m_nc );
   int ie = u.m_ie, ib=u.m_ib, je=u.m_je, jb=u.m_jb, kb=u.m_kb;//,ke=u.m_ke;
   MPI_Status status;
   cudaError_t retcode;
   dim3 gridsize, blocksize;
   gridsize.x  = m_gpu_gridsize[0] * m_gpu_gridsize[1] * m_gpu_gridsize[2];
   gridsize.y  = 1;
   gridsize.z  = 1;
   blocksize.x = m_gpu_blocksize[0] * m_gpu_blocksize[1] * m_gpu_blocksize[2];
   blocksize.y = 1;
   blocksize.z = 1;

   int ni = m_iEnd[g] - m_iStart[g] + 1;
   int nj = m_jEnd[g] - m_jStart[g] + 1;
   int nk = m_kEnd[g] - m_kStart[g] + 1;
   int n_m_ppadding1 = 3*nj*nk*m_ppadding;
   int n_m_ppadding2 = 3*ni*nk*m_ppadding;
   int idx_left = 0;
   int idx_right = n_m_ppadding2;
   int idx_up = 2*n_m_ppadding2;
   int idx_down = 2*n_m_ppadding2 + n_m_ppadding1;
   int n_m_ppadding_total = 2*(n_m_ppadding1+n_m_ppadding2);

   if( u.m_nc == 1 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      int grid = g;
      // X-direction communication
      MPI_Sendrecv( &u(ie-(2*m_ppadding-1),jb,kb,true), 1, m_send_type1[2*grid], m_neighbor[1], xtag1,
                    &u(ib,jb,kb,true), 1, m_send_type1[2*grid], m_neighbor[0], xtag1,
                    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(ib+m_ppadding,jb,kb,true), 1, m_send_type1[2*grid], m_neighbor[0], xtag2,
                    &u(ie-(m_ppadding-1),jb,kb,true), 1, m_send_type1[2*grid], m_neighbor[1], xtag2,
                    m_cartesian_communicator, &status );
      //Y-direction communication
      MPI_Sendrecv( &u(ib,je-(2*m_ppadding-1),kb,true), 1, m_send_type1[2*grid+1], m_neighbor[3], ytag1,
                    &u(ib,jb,kb,true), 1, m_send_type1[2*grid+1], m_neighbor[2], ytag1,
                    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(ib,jb+m_ppadding,kb,true), 1, m_send_type1[2*grid+1], m_neighbor[2], ytag2,
                    &u(ib,je-(m_ppadding-1),kb,true), 1, m_send_type1[2*grid+1], m_neighbor[3], ytag2,
                    m_cartesian_communicator, &status );
   }
   else if( u.m_nc == 3 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;


      if(m_corder)
         BufferToHaloKernel_dev_rev<<<gridsize, blocksize>>>( &u(1,ib,jb+m_ppadding,kb,true), &u(1,ib,je-(2*m_ppadding-1),kb,true),
                                                   &u(1,ie-(2*m_ppadding-1),jb,kb,true), &u(1,ib+m_ppadding,jb,kb,true),
                                                   &dev_SideEdge_Send[g][idx_left], &dev_SideEdge_Send[g][idx_right],
                                                   &dev_SideEdge_Send[g][idx_up], &dev_SideEdge_Send[g][idx_down],
                    ni, nj, nk, m_ppadding, m_neighbor[0],  m_neighbor[1], m_neighbor[2], m_neighbor[3], MPI_PROC_NULL );

      else
         BufferToHaloKernel_dev<<<gridsize, blocksize>>>( &u(1,ib,jb+m_ppadding,kb,true), &u(1,ib,je-(2*m_ppadding-1),kb,true),
                                                   &u(1,ie-(2*m_ppadding-1),jb,kb,true), &u(1,ib+m_ppadding,jb,kb,true),
                                                   &dev_SideEdge_Send[g][idx_left], &dev_SideEdge_Send[g][idx_right],
                                                   &dev_SideEdge_Send[g][idx_up], &dev_SideEdge_Send[g][idx_down],
                    ni, nj, nk, m_ppadding, m_neighbor[0],  m_neighbor[1], m_neighbor[2], m_neighbor[3], MPI_PROC_NULL );

      CheckCudaCall(cudaGetLastError(), "BufferToHaloKernel<<<,>>>(...)", __FILE__, __LINE__);
      SafeCudaCall(cudaStreamSynchronize(NULL));

#ifdef SW4_CUDA_AWARE_MPI

     MPI_Sendrecv(&dev_SideEdge_Send[g][idx_up], n_m_ppadding1, m_mpifloat, m_neighbor[1], xtag1, &dev_SideEdge_Recv[g][idx_down],
                        n_m_ppadding1, m_mpifloat, m_neighbor[0], xtag1, m_cartesian_communicator, &status);

     MPI_Sendrecv(&dev_SideEdge_Send[g][idx_down], n_m_ppadding1, m_mpifloat, m_neighbor[0], xtag2, &dev_SideEdge_Recv[g][idx_up],
                        n_m_ppadding1, m_mpifloat, m_neighbor[1], xtag2, m_cartesian_communicator, &status);

     MPI_Sendrecv(&dev_SideEdge_Send[g][idx_right], n_m_ppadding2, m_mpifloat, m_neighbor[3], ytag2, &dev_SideEdge_Recv[g][idx_left],
                        n_m_ppadding2, m_mpifloat, m_neighbor[2], ytag2, m_cartesian_communicator, &status);

     MPI_Sendrecv(&dev_SideEdge_Send[g][idx_left], n_m_ppadding2, m_mpifloat, m_neighbor[2], ytag1, &dev_SideEdge_Recv[g][idx_right],
                        n_m_ppadding2, m_mpifloat, m_neighbor[3], ytag1, m_cartesian_communicator, &status);


#else

     retcode = cudaMemcpy(m_SideEdge_Send[g], dev_SideEdge_Send[g], n_m_ppadding_total*sizeof(float_sw4), cudaMemcpyDeviceToHost);
     if( retcode != cudaSuccess )
     {
        cout << "Error cmmunicate_array cudaMemcpy returned (DeviceToHost) " << cudaGetErrorString(retcode) << endl;
        exit(1);
     }

     MPI_Sendrecv(&m_SideEdge_Send[g][idx_left], n_m_ppadding2, m_mpifloat, m_neighbor[2], ytag1, &m_SideEdge_Recv[g][idx_right],
                         n_m_ppadding2, m_mpifloat, m_neighbor[3], ytag1, m_cartesian_communicator, &status);

     MPI_Sendrecv(&m_SideEdge_Send[g][idx_right], n_m_ppadding2, m_mpifloat, m_neighbor[3], ytag2, &m_SideEdge_Recv[g][idx_left],
                        n_m_ppadding2, m_mpifloat, m_neighbor[2], ytag2, m_cartesian_communicator, &status);

     MPI_Sendrecv(&m_SideEdge_Send[g][idx_up], n_m_ppadding1, m_mpifloat, m_neighbor[1], xtag1, &m_SideEdge_Recv[g][idx_down],
                        n_m_ppadding1, m_mpifloat, m_neighbor[0], xtag1, m_cartesian_communicator, &status);

     MPI_Sendrecv(&m_SideEdge_Send[g][idx_down], n_m_ppadding1, m_mpifloat, m_neighbor[0], xtag2, &m_SideEdge_Recv[g][idx_up],
                        n_m_ppadding1, m_mpifloat, m_neighbor[1], xtag2, m_cartesian_communicator, &status);


     retcode = cudaMemcpy(dev_SideEdge_Recv[g], m_SideEdge_Recv[g], n_m_ppadding_total*sizeof(float_sw4), cudaMemcpyHostToDevice);
     if( retcode != cudaSuccess )
     {
        cout << "Error cmmunicate_array cudaMemcpy returned (Host2Device) " << cudaGetErrorString(retcode) << endl;
        exit(1);
     }

#endif

     if(m_corder)
        HaloToBufferKernel_dev_rev<<<gridsize, blocksize>>>( &u(1,ib,jb,kb,true), &u(1,ib,je-(m_ppadding-1),kb,true),
                                                  &u(1,ie-(m_ppadding-1),jb,kb,true), &u(1,ib,jb,kb,true),
                                                  &dev_SideEdge_Recv[g][idx_left], &dev_SideEdge_Recv[g][idx_right],
                                                  &dev_SideEdge_Recv[g][idx_up], &dev_SideEdge_Recv[g][idx_down], ni, nj, nk, m_ppadding,
                                     m_neighbor[0],  m_neighbor[1], m_neighbor[2], m_neighbor[3], MPI_PROC_NULL );

     else
        HaloToBufferKernel_dev<<<gridsize, blocksize>>>( &u(1,ib,jb,kb,true), &u(1,ib,je-(m_ppadding-1),kb,true),
                                                  &u(1,ie-(m_ppadding-1),jb,kb,true), &u(1,ib,jb,kb,true),
                                                  &dev_SideEdge_Recv[g][idx_left], &dev_SideEdge_Recv[g][idx_right],
                                                  &dev_SideEdge_Recv[g][idx_up], &dev_SideEdge_Recv[g][idx_down], ni, nj, nk, m_ppadding,
                                     m_neighbor[0],  m_neighbor[1], m_neighbor[2], m_neighbor[3], MPI_PROC_NULL );


     CheckCudaCall(cudaGetLastError(), "HaloToBufferKernel<<<,>>>(...)", __FILE__, __LINE__);
     SafeCudaCall(cudaStreamSynchronize(NULL));


   }
#endif
}


//-----------------------------------------------------------------------
void EW::enforceBCCU( vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                      float_sw4 t, vector<float_sw4**> & a_BCForcing, int st )
{
#ifdef SW4_CUDA
  dim3 gridsize, blocksize;
  gridsize.x  = m_gpu_gridsize[0];
  gridsize.y  = m_gpu_gridsize[1];
  gridsize.z  = m_gpu_gridsize[2];
  blocksize.x = m_gpu_blocksize[0];
  blocksize.y = m_gpu_blocksize[1];
  blocksize.z = m_gpu_blocksize[2];

  float_sw4 om=0, ph=0, cv=0;
  for(int g=0 ; g<mNumberOfGrids; g++ )
  {
    if( m_corder )
      bcfortsg_dev_indrev<<<gridsize, blocksize, 0, m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                                                              dev_BndryWindow[g], m_global_nx[g], m_global_ny[g], m_global_nz[g], a_U[g].dev_ptr(),
                                                                              mGridSize[g], dev_bcType[g], a_Mu[g].dev_ptr(), a_Lambda[g].dev_ptr(),
                                                                              t, dev_BCForcing[g][0], dev_BCForcing[g][1], dev_BCForcing[g][2],
                                                                              dev_BCForcing[g][3], dev_BCForcing[g][4], dev_BCForcing[g][5],
                                                                              om, ph, cv, dev_sg_str_x[g], dev_sg_str_y[g] );
    else
      bcfortsg_dev<<<gridsize, blocksize, 0, m_cuobj->m_stream[st]>>>( m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g], m_kEnd[g],
                                                                       dev_BndryWindow[g], m_global_nx[g], m_global_ny[g], m_global_nz[g], a_U[g].dev_ptr(),
                                                                       mGridSize[g], dev_bcType[g], a_Mu[g].dev_ptr(), a_Lambda[g].dev_ptr(),
                                                                       t, dev_BCForcing[g][0], dev_BCForcing[g][1], dev_BCForcing[g][2],
                                                                       dev_BCForcing[g][3], dev_BCForcing[g][4], dev_BCForcing[g][5],
                                                                       om, ph, cv, dev_sg_str_x[g], dev_sg_str_y[g] );
  }
#endif
}

//-----------------------------------------------------------------------

void EW::cartesian_bc_forcingCU( float_sw4 t, vector<float_sw4**> & a_BCForcing,
                                 vector<Source*>& a_sources , int st)
// assign the boundary forcing arrays dev_BCForcing[g][side]
{
#ifdef SW4_CUDA
  cudaError_t retcode;
  for(int g=0 ; g<mNumberOfGrids; g++ )
  {
    if( m_point_source_test )
    {
      for( int side=0 ; side < 6 ; side++ )
      {
        size_t nBytes = sizeof(float_sw4)*3*m_NumberOfBCPoints[g][side];
        if( m_bcType[g][side] == bDirichlet )
        {
          get_exact_point_source( a_BCForcing[g][side], t, g, *a_sources[0], &m_BndryWindow[g][6*side] );
          retcode = cudaMemcpyAsync( dev_BCForcing[g][side], a_BCForcing[g][side], nBytes, cudaMemcpyHostToDevice, m_cuobj->m_stream[st]);
          if( retcode != cudaSuccess )
            cout << "Error, EW::cartesian_bc_forcing_CU cudaMemcpy x returned " << cudaGetErrorString(retcode) << endl;
        }
        else
        {
          cudaMemsetAsync( dev_BCForcing[g][side], 0, nBytes , m_cuobj->m_stream[st]);
        }
      }
    }
    else
    {
      for( int side=0 ; side < 6 ; side++ )
      {
        size_t nBytes = sizeof(float_sw4)*3*m_NumberOfBCPoints[g][side];
        cudaMemsetAsync( dev_BCForcing[g][side], 0, nBytes , m_cuobj->m_stream[st]);
      }
    }
  }
#endif
}

//-----------------------------------------------------------------------

void EW::copy_bcforcing_arrays_to_device()
{

#ifdef SW4_CUDA
   //Set up boundary data array on the deivec
  if(m_ndevice > 0 )
  {
    cudaError_t retcode;
    dev_BCForcing.resize(mNumberOfGrids);
    for( int g = 0; g <mNumberOfGrids; g++ )
    {
      dev_BCForcing[g] = new float_sw4*[6];
      for (int side=0; side < 6; side++)
      {
        dev_BCForcing[g][side] = NULL;
        if (m_bcType[g][side] == bStressFree || m_bcType[g][side] == bDirichlet || m_bcType[g][side] == bSuperGrid)
        {
          size_t nBytes = sizeof(float_sw4)*3*m_NumberOfBCPoints[g][side];
          retcode  = cudaMalloc((void**) &dev_BCForcing[g][side], nBytes );
          if( retcode != cudaSuccess )
          {
             cout << "Error, EW::copy_bcforcing_arrays_to_device cudaMalloc x returned " << cudaGetErrorString(retcode) << endl;
             exit(-1);
          }
          retcode = cudaMemcpy( dev_BCForcing[g][side], BCForcing[g][side], nBytes, cudaMemcpyHostToDevice );
          if( retcode != cudaSuccess )
          {
            cout << "Error, EW::copy_bcforcing_arrays_to_device cudaMemcpy x returned " << cudaGetErrorString(retcode) << endl;
            exit(-1);
          }
        }
      }
    }
  }
#endif
}

//-----------------------------------------------------------------------

void EW::copy_bctype_arrays_to_device()
{
#ifdef SW4_CUDA
  // Set up boundary type array on the deivec
  if(m_ndevice > 0 )
  {
    cudaError_t retcode;
    dev_bcType.resize(mNumberOfGrids);
    for( int g = 0; g <mNumberOfGrids; g++ )
    {
      size_t nBytes = sizeof(boundaryConditionType)*6;
      retcode = cudaMalloc( (void**) &dev_bcType[g], nBytes );
      if( retcode != cudaSuccess )
        cout << "Error, EW::copy_bctype_arrays_to_device cudaMalloc x returned " << cudaGetErrorString(retcode) << endl;
      retcode = cudaMemcpy( dev_bcType[g], m_bcType[g], nBytes, cudaMemcpyHostToDevice );
      if( retcode != cudaSuccess )
        cout << "Error, EW::copy_bctype_arrays_to_device cudaMemcpy x returned " << cudaGetErrorString(retcode) << endl;
    }
  }
#endif
}

//-----------------------------------------------------------------------

void EW::copy_bndrywindow_arrays_to_device()
{
#ifdef SW4_CUDA
  //Set up boundary window array on the deivec
  if(m_ndevice > 0 )
  {
    cudaError_t retcode;
    dev_BndryWindow.resize(mNumberOfGrids);
    for( int g = 0; g <mNumberOfGrids; g++ )
    {
      size_t nBytes = sizeof(int)*36;
      retcode = cudaMalloc( (void**) &dev_BndryWindow[g], nBytes );
      if( retcode != cudaSuccess )
        cout << "Error, EW::copy_bndrywindow_arrays_to_device cudaMalloc x returned " << cudaGetErrorString(retcode) << endl;
      retcode = cudaMemcpy( dev_BndryWindow[g], m_BndryWindow[g], nBytes, cudaMemcpyHostToDevice );
      if( retcode != cudaSuccess )
        cout << "Error, EW::copy_bndrywindow_arrays_to_device cudaMemcpy x returned " << cudaGetErrorString(retcode) << endl;
    }
  }
#endif
}

