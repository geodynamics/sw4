#ifndef MATERIALPARCARTESIANVSVP_H
#define MATERIALPARCARTESIANVSVP_H

#include "MaterialParameterization.h"

class MaterialParCartesianVsVp : public MaterialParameterization
{
private:
   double m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz;  // Grid given as x_i = xmin + (i-1)*hx,  y_j=ymin+(j-1)*hy, z_k=zmin+(k-1)*hz
   double m_xmax, m_ymax, m_zmax ; // Not really needed.
   int m_nx, m_ny, m_nz;                      // i=1,..,nx, j=1,..,ny, k=1,..,nz
   int m_init;
   Sarray m_rho, m_mu, m_lambda;
   Sarray m_cs, m_cp; // ?needed now ? 
   void interpolate_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
				std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   //void interpolate_base_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
	//			std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
public:
   MaterialParCartesianVsVp( EW* a_ew, int nx, int ny, int nz, int init, char* fname );
   void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
		      std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda,
            float_sw4 vp_min, float_sw4 vp_max, float_sw4 vs_min, float_sw4 vs_max, int wave_mode);
   void get_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
			std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda, int nr );
   void get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfs, double* dfm,
		      std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
		      std::vector<Sarray>& a_lambda,
		      std::vector<Sarray>& a_gradrho, std::vector<Sarray>& a_gradmu,
		      std::vector<Sarray>& a_gradlambda);
   void smooth_gradient(double* dfs, std::vector<Sarray>& a_Rho, std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda, float_sw4 freq, float_sw4 sz);


   //   void gradient_transformation( std::vector<Sarray>& a_rho,
   //				 std::vector<Sarray>& a_mu,
   //				 std::vector<Sarray>& a_lambda,
   //				 std::vector<Sarray>& a_gradrho,
   //				 std::vector<Sarray>& a_gradmu,
   //				 std::vector<Sarray>& a_gradlambda );
   void interpolate_pseudohessian(int nmpars, double* phs, int nmpard, double* phm, 
                                  std::vector<Sarray>& phgrid);

   ssize_t parameter_index( int ip, int jp, int kp, int grid, int var );
   ssize_t local_index( size_t ind_global );
   void set_scalefactors( int nmpars, double* sfs, double rho_ref, double mu_ref, double lambda_ref, 
			  double vs_ref, double vp_ref );
   void subtract_base_mtrl( int nms, double* xms );

   double getXmin() const { return m_xmin; }
   double getDx() const { return m_hx; }
   int getNX() const { return m_nx; }
   double getYmin() const { return m_ymin; }
   double getDy() const { return m_hy; }
   int getNY() const { return m_ny; }
   double getZmin() const { return m_zmin; }
   double getDz() const { return m_hz; }
   int getNZ() const { return m_nz; }  
   int get_varcase(){return 3;};
   void write_dfm_hdf5(double* dfm, std::string fname,  MPI_Comm comm) {printf("%s not supported!\n", __func__);}
};

#endif
