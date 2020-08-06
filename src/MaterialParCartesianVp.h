#ifndef MATERIALPARCARTESIANVP_H
#define MATERIALPARCARTESIANVP_H

#include "MaterialParameterization.h"

class MaterialParCartesianVp : public MaterialParameterization
{
private:
   double m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz;  // Grid given as x_i = xmin + (i-1)*hx,  y_j=ymin+(j-1)*hy, z_k=zmin+(k-1)*hz
   double m_xmax, m_ymax, m_zmax ; // Not really needed.
   int m_nx, m_ny, m_nz;                      // i=1,..,nx, j=1,..,ny, k=1,..,nz
   int m_init;
   Sarray m_rho, m_mu, m_lambda;
   Sarray m_cs, m_cp; // ?needed now ? 
   double m_ratio, m_gamma;
   bool m_fixrho;
   void interpolate_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
				std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
public:
   MaterialParCartesianVp( EW* a_ew, int nx, int ny, int nz, int init, char* fname, double ratio, double gamma, bool fixrho );
   void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
		      std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   void get_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
			std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   void get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfs, double* dfm,
		      std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
		      std::vector<Sarray>& a_lambda,
		      std::vector<Sarray>& a_gradrho, std::vector<Sarray>& a_gradmu,
		      std::vector<Sarray>& a_gradlambda);
   void smooth_gradient(double* dfs);
   void smooth_gradient( std::vector<Sarray>& a_grad ) {};
   void interpolate_pseudohessian(int nmpars, double* phs, int nmpard, double* phm, 
                                  std::vector<Sarray>& phgrid);

   //   void gradient_transformation( std::vector<Sarray>& a_rho,
   //				 std::vector<Sarray>& a_mu,
   //				 std::vector<Sarray>& a_lambda,
   //				 std::vector<Sarray>& a_gradrho,
   //				 std::vector<Sarray>& a_gradmu,
   //				 std::vector<Sarray>& a_gradlambda );
   void subtract_base_mtrl( int nms, double* xms );
   ssize_t parameter_index( int ip, int jp, int kp, int grid, int var );
   ssize_t local_index( size_t ind_global );
   void set_scalefactors( int nmpars, double* sfs, double rho_ref, double mu_ref, double lambda_ref, 
			  double vs_ref, double vp_ref );
   int get_varcase(){return 4;}
};

#endif
