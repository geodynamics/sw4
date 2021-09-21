#ifndef MATERIALPARAMETERIZATION_H
#define MATERIALPARAMETERIZATION_H

#include <vector>
#include <string>
#include <Sarray.h>
#include <Source.h>

class EW;

class MaterialParameterization
{
protected:
   int m_myrank, m_nmd, m_nms, m_nmd_global;
   EW* m_ew;
   char* m_filename;
   std::string m_path;
public:
   MaterialParameterization( EW* a_ew, char* fname );
   virtual void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
			      std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )=0;
   virtual void get_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
   				std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda, int nr )=0;

   //virtual void get_material( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho,
	//		      std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda, 
   //            float_sw4 vp_min=-100, float_sw4 vp_max=-100, float_sw4 vs_min=-100, float_sw4 vs_max=-100,int wave_mode=2)=0;

   virtual void get_base_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
				std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )=0;
            
   //virtual void get_base_parameters( int nmd, double* xmd, int nms, double* xms, std::vector<Sarray>& a_rho, 
	//			std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )=0;

   virtual void get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfs, double* dfm,
			      std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
			      std::vector<Sarray>& a_lambda,
			      std::vector<Sarray>& a_gradrho, std::vector<Sarray>& a_gradmu,
			      std::vector<Sarray>& a_gradlambda)=0;

   virtual void smooth_gradient(double* dfs, std::vector<Sarray>& a_Rho, std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda, float_sw4 freq, float_sw4 sz)=0;
   //   virtual void gradient_transformation( std::vector<Sarray>& a_rho,    std::vector<Sarray>& a_mu,
   //					 std::vector<Sarray>& a_lambda, std::vector<Sarray>& a_gradrho, 
   //					 std::vector<Sarray>& a_gradmu, std::vector<Sarray>& a_gradlambda ){};

   //   virtual void perturb_material( int ip, int jp, int kp, int grid, int var, double h, double* xs, double* xm ) = 0;
   virtual void interpolate_pseudohessian(int nmpars, double* phs, int nmpard, double* phm, 
                                          std::vector<Sarray>& phgrid)=0;

   virtual ssize_t parameter_index( int ip, int jp, int kp, int grid, int var )=0;
   virtual ssize_t local_index( size_t ind_global )=0;
   virtual void set_scalefactors( int nmpars, double* sfs, double rho_ref, double mu_ref, double lambda_ref, 
				  double vs_ref, double vp_ref );
   virtual void get_regularizer( int nmd, double* xmd, int nms, double* xms, 
				 double* xmd0, double* xms0, double regcoeff,
				 std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
				 std::vector<Sarray>& a_lambda, double& mf_reg,
				 double* sfd, double* sfs, bool compute_derivative, 
				 double* dmfd_reg, double* dmfs_reg );
   virtual int get_varcase()=0; 
   virtual void write_dfm_hdf5(double* dfm, std::string fname,  MPI_Comm comm)=0;
   void get_nr_of_parameters( int& nms, int& nmd, int& nmd_global ) const;
   //   void parameters_from_basematerial( int nmd, double* xmd, int nms, double* xms );
   //   void store_material( std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda );
   //   void constrain_material( int nmd, double* xmd, int nms, double* xms );
   void write_parameters( const char* filename, int nms, double* xms ); // Only shared parameters for now
   void read_parameters( const char* filename, int nms, double* xms ); // Only shared parameters for now
   void write_parameters( int nms, double* xms ); 
   void write_parameters_dist( const char* outfile, int nmd, double* xmd );
   void read_parameters( int nms, double* xms ); // Only shared parameters for now
   void set_path( std::string path ){m_path = path;}

   virtual void limit_x( int nmd, double* xmd, int nms, double* xms,
                 float_sw4 vsmin, float_sw4 vsmax, 
                 float_sw4 vpmin, float_sw4 vpmax )=0;
   virtual void limit_df( int nmd, double* dfd, int nms, double* dfs )=0;

   virtual double getXmin() const =0;
   virtual double getDx() const=0;
   virtual int    getNX() const=0;
   virtual double getYmin() const=0;
   virtual double getDy() const =0;
   virtual int getNY() const =0;
   virtual double getZmin() const =0;
   virtual double getDz() const =0;
   virtual int getNZ() const =0;
};

#endif
