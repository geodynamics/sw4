#ifndef MATERIALPARCART_H
#define MATERIALPARCART_H

#include "MaterialParameterization.h"

// Represents material parameters on coarsened grid, given as
//    x_i= xmin + (i-1)*hx, i=1,..,nx
//    y_j= ymin + (j-1)*hy, j=1,..,ny
//    z_k= zmin + (k-1)*hz, k=1,..,nz
//
//
class MaterialParCart : public MaterialParameterization
{
protected:
   float_sw4 m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz;   
   float_sw4 m_xmax, m_ymax, m_zmax, m_ratio;
   float_sw4 m_amplitude, m_omega;

   bool m_global;
   int m_nx, m_ny, m_nz;  // Global dimensions of parameter grid.
   int m_ib, m_ie, m_jb, m_je, m_kb, m_ke; // Index limits of par. grid in this processor.
   int m_ibint, m_ieint, m_jbint, m_jeint, m_kbint, m_keint; // interior grid in this processor.
   int m_ibpp, m_iepm, m_jbpp, m_jepm; // Index limits at proc p+1 (pp) and p-1 (pm).
   int m_init;
   int m_variables;
   std::vector<bool> m_limited;

   //   Sarray m_rho, m_mu, m_lambda;
   //   Sarray m_cs, m_cp; 
   //   EW* m_ew;
   void interpolate_parameters( int nmd, double* xmd, int nms, double* xms, 
                                std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu, 
                                std::vector<Sarray>& a_lambda, bool update );
   void find_lims( int ib, int ie, int iepm, int ibpp, int& ibint, int& ieint );

   bool compute_overlap( bool force_shared );
   void getwgh( float_sw4 ai, float_sw4 wgh[2], int& sl, int& su );
   void communicate( Sarray& u );
   void communicate_add( Sarray& u );
   void interpolate_gradient( int g, Sarray& a_gradrho, Sarray& a_gradmu, 
                              Sarray& a_gradlambda, Sarray& gradc );
   void get_gradient_shared( int nms, double* xms, double* dfs, 
                             std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
                             std::vector<Sarray>& a_lambda, std::vector<Sarray>& a_gradrho, 
                             std::vector<Sarray>& a_gradmu, std::vector<Sarray>& a_gradlambda );
   void interpolate_to_coarse( std::vector<Sarray>& rhogrid, std::vector<Sarray>& mugrid,
                               std::vector<Sarray>& lambdagrid, Sarray& rho, Sarray& mu, 
                               Sarray& lambda, bool update );
   void interpolate_to_coarse_vel(std::vector<Sarray>& rhogrid, std::vector<Sarray>& mugrid,
                                  std::vector<Sarray>& lambdagrid,
                                  Sarray& rho, Sarray& cs, Sarray& cp, bool update );

public:
   MaterialParCart( EW* a_ew, int nx, int ny, int nz, int init, int varcase, char* fname,
                    float_sw4 amp, float_sw4 omega, bool force_shared );

   virtual void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
		      std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );

   virtual void interpolate( Sarray& matcart, int g, Sarray& rho, Sarray& mu, Sarray& lambda );

   virtual void get_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
                                std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda, int nr );

   virtual void get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfs, double* dfm,
		      std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
		      std::vector<Sarray>& a_lambda,
		      std::vector<Sarray>& a_gradrho, std::vector<Sarray>& a_gradmu,
		      std::vector<Sarray>& a_gradlambda );

   virtual void interpolate_pseudohessian(int nmpars, double* phs, int nmpard, double* phm, 
                                          std::vector<Sarray>& phgrid);

   virtual void set_scalefactors( int nmpars, double* sfs, double rho_ref, 
                                  double mu_ref, double lambda_ref, 
                                  double vs_ref, double vp_ref );

   ssize_t parameter_index( int ip, int jp, int kp, int grid, int var );
   ssize_t local_index( size_t ind_global );

   //   void subtract_base_mtrl( int nms, double* xms );
   int get_varcase(){return m_variables;};
   //   int get_varcase(){return 2;}
   void write_dfm_hdf5(double* dfm, std::string fname, MPI_Comm comm);
   void get_local_grid_begin_end(int* begin, int* end) {begin[0] = m_ibint; begin[1] = m_jbint; begin[2] = m_kbint; 
                                                        end[0] = m_ieint; end[1] = m_jeint; end[2] = m_keint;};
   void get_global_grid_size(int *grid_size) {printf("%s\n", __func__); grid_size[0] = m_nx; grid_size[1] = m_ny; grid_size[2] = m_nz;}


   virtual void limit_x( int nmd, double* xmd, int nms, double* xms,
                 float_sw4 vsmin, float_sw4 vsmax, 
                 float_sw4 vpmin, float_sw4 vpmax );
   virtual void limit_df( int nmd, double* dfd, int nms, double* dfs );

   double getXmin() const { return m_xmin; }
   double getDx() const { return m_hx; }
   int getNX() const { return m_nx; }
   double getYmin() const { return m_ymin; }
   double getDy() const { return m_hy; }
   int getNY() const { return m_ny; }
   double getZmin() const { return m_zmin; }
   double getDz() const { return m_hz; }
   int getNZ() const { return m_nz; }   
};

#endif
