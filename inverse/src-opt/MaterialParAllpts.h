#ifndef MATERIALPARALLPTS_H
#define MATERIALPARALLPTS_H

#include "MaterialParameterization.h"

class MaterialParAllpts : public MaterialParameterization
{
   std::vector<size_t> m_npts_per_grid;
public:
   MaterialParAllpts(EW* a_ew, char* fname );
   void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
					 std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   void get_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
					   std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   void get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfms, double* dfmd,
					 std::vector<Sarray>& a_gradrho, std::vector<Sarray>& a_gradmu,
					 std::vector<Sarray>& a_gradlambda );
   //   void perturb_material( int ip, int jp, int kp, int grid, int var, double h, double* xs, double* xm );
   ssize_t parameter_index( int ip, int jp, int kp, int grid, int var );
   ssize_t local_index( size_t ind_global );
};

#endif
