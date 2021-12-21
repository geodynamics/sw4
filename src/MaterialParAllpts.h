#ifndef MATERIALPARALLPTS_H
#define MATERIALPARALLPTS_H

#include "MaterialParameterization.h"

class MaterialParAllpts : public MaterialParameterization
{
   std::vector<size_t> m_npts_per_grid, m_npts_per_grid_local;
   int m_nc; // Number of components per grid point
       // RML --> (rho,mu,lambda), RVSVP--> (rho,cs,cp)
       // CSCP --> (cs,cp)   CP --> (cp)
   enum Variabletype {RML,RCSCP,CSCP,CP};
   Variabletype  m_variables;
   double m_ratio; // cp/cs ratio when m_variables == cp
public:
   MaterialParAllpts(EW* a_ew, char* fname, int variables );
   //void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
	//				 std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda,
   //              float_sw4 vp_min, float_sw4 vp_max, float_sw4 vs_min, float_sw4 vs_max,int wave_mode);
   //void get_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
	//				   std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   void get_base_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
				std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda ) {};
   void smooth_gradient(double* dfs, std::vector<Sarray>& a_Rho, std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda, float_sw4 freq, float_sw4 sz) {};

   virtual void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
					 std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   void get_material( int nmd, double* xmd, int nms,
					     double* xms, std::vector<Sarray>& a_rho,
					     std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda,
                    float_sw4 vp_min, float_sw4 vp_max, float_sw4 vs_min, float_sw4 vs_max, int wave_mode);

   virtual void get_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
                        std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda, int nr );
   //virtual void get_base_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
	//			std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   virtual void get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfms, double* dfmd,
		      std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
		      std::vector<Sarray>& a_lambda, 
		      std::vector<Sarray>& a_gradrho, std::vector<Sarray>& a_gradmu,
		      std::vector<Sarray>& a_gradlambda );

   virtual void set_scalefactors( int nmpars, double* sfs, 
			  double rho_ref, double mu_ref, double lambda_ref, double vs_ref, double vp_ref );
   //   void perturb_material( int ip, int jp, int kp, int grid, int var, double h, double* xs, double* xm );
   virtual void interpolate_pseudohessian(int nmpars, double* phs, int nmpard, double* phm, 
                                          std::vector<Sarray>& phgrid);
   virtual ssize_t parameter_index( int ip, int jp, int kp, int grid, int var );
   virtual ssize_t local_index( size_t ind_global );

   virtual void get_regularizer( int nmd, double* xmd, int nms, double* xms, 
			 double* xmd0, double* xms0, double regcoeff,
			 std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
			 std::vector<Sarray>& a_lambda, double& mf_reg,
			 double* sfd, double* sfs, bool compute_derivative, 
			 double* dmfd_reg, double* dmfs_reg );
   virtual int get_varcase();
   virtual void write_dfm_hdf5(double* dfm, std::string fname,  MPI_Comm comm) {printf("%s not supported!\n", __func__);}

   virtual void limit_x( int nmd, double* xmd, int nms, double* xms,
                         float_sw4 vsmin, float_sw4 vsmax, 
                         float_sw4 vpmin, float_sw4 vpmax ){};
   virtual void limit_df( int nmd, double* dfd, int nms, double* dfs ){};

   double getXmin() const {  }
   double getDx() const {  }
   int getNX() const {  }
   double getYmin() const {  }
   double getDy() const {  }
   int getNY() const {  }
   double getZmin() const {  }
   double getDz() const {  }
   int getNZ() const {  }
};

#endif
