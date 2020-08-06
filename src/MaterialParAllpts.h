#ifndef MATERIALPARALLPTS_H
#define MATERIALPARALLPTS_H

#include "MaterialParameterization.h"

class MaterialParAllpts : public MaterialParameterization
{
   std::vector<size_t> m_npts_per_grid;
   int m_nc; // Number of components per grid point
       // RML --> (rho,mu,lambda), RVSVP--> (rho,cs,cp)
       // CSCP --> (cs,cp)   CP --> (cp)
   enum Variabletype {RML,RCSCP,CSCP,CP};
   Variabletype  m_variables;
   double m_ratio; // cp/cs ratio when m_variables == cp
public:
   MaterialParAllpts(EW* a_ew, char* fname, int variables );
   void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
					 std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   void get_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
					   std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   void get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfms, double* dfmd,
		      std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
		      std::vector<Sarray>& a_lambda, 
		      std::vector<Sarray>& a_gradrho, std::vector<Sarray>& a_gradmu,
		      std::vector<Sarray>& a_gradlambda, int rank);
   void smooth_gradient( std::vector<Sarray>& a_grad);

   void set_scalefactors( int nmpars, double* sfs, 
			  double rho_ref, double mu_ref, double lambda_ref, double vs_ref, double vp_ref );
   //   void perturb_material( int ip, int jp, int kp, int grid, int var, double h, double* xs, double* xm );
   void interpolate_pseudohessian(int nmpars, double* phs, int nmpard, double* phm, 
                                  std::vector<Sarray>& phgrid);
   ssize_t parameter_index( int ip, int jp, int kp, int grid, int var );
   ssize_t local_index( size_t ind_global );

   void get_regularizer( int nmd, double* xmd, int nms, double* xms, 
			 double* xmd0, double* xms0, double regcoeff,
			 std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
			 std::vector<Sarray>& a_lambda, double& mf_reg,
			 double* sfd, double* sfs, bool compute_derivative, 
			 double* dmfd_reg, double* dmfs_reg );
   int get_varcase();
};

#endif
