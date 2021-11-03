#ifndef MOPT_H
#define MOPT_H

#include <vector>
#include "MaterialParameterization.h"
#include "Image.h"
#include "Image3D.h"
#include "SfileOutput.h"

class EW;
class MaterialParCartesian;

class Mopt
{
 private:
   EW* m_ew;
   int m_myrank;
   double m_rhoscale, m_muscale, m_lambdascale, m_misfitscale;
   double m_vsscale, m_vpscale;
   double m_typrho, m_typmu, m_typlambda;
   bool m_typrhoset, m_typmuset, m_typlambdaset;
   double m_rhosffactor, m_musffactor, m_lambdasffactor;
   bool m_use_pseudohessian;

   void badOption(string name, char* option) const;
   bool startswith(const char begin[], char *line);
   void processMaterialParCart( char* buffer );
   void processMaterialAllpts( char* buffer );
   void processMrun( char* buffer );
   void processMscalefactors( char* buffer );
   void processLBFGS( char* buffer );
   void processNLCG( char* buffer );
   void processMfsurf( char* buffer );
   void processMimage( char* buffer, bool use_hdf5 );
   void processM3Dimage( char* buffer );
   void processSfileoutput( char* buffer );
   void processMtypx( char* buffer );
   void processMfileio( char* buffer );
   void processMregularize( char* buffer );
 public:
   Mopt( EW* a_ew );
   bool parseInputFileOpt( std::string filename );
   void get_scalefactors( double& rhoscale, double& muscale,
			  double& lambdascale );

   void set_sscalefactors( );
   void set_sourcescalefactors( int nspar, double* sfs );
   void set_dscalefactors( );
   void set_typx( int nmpar, double* sf, double* typx );
   void init_pseudohessian( vector<Sarray>& ph );
   int get_pseudo_hessian_case();
   float_sw4 get_vp_min() const { return m_vp_min; }
   float_sw4 get_vp_max() const { return m_vp_max; }
   float_sw4 get_vs_min() const { return m_vs_min; }
   float_sw4 get_vs_max() const { return m_vs_max; }
   int get_wave_mode() const { return m_wave_mode; }
   float get_twin_shift() const { return m_twin_shift; }
   float get_twin_scale() const { return m_twin_scale; }

   const string& getPath() const {return m_path;}
   EW* get_EWptr() const {return m_ew;}
   void set_baseMat(double* xs, double* xm );

   enum Misfittype {L2,CROSSCORR};
   Misfittype m_misfit;

   int m_opttest, m_nspar;
   int m_maxit, m_maxsubit, m_nbfgs_vectors, m_optmethod, m_ihess_guess;
   bool m_dolinesearch, m_fletcher_reeves, m_wolfe, m_mcheck, m_output_ts;
   bool m_misfit1d_images;
   bool m_test_regularizer;
   bool m_do_profiling;
   double m_tolerance;
   // Test cases
   int m_var, m_var2, m_itest, m_jtest, m_ktest, m_itest2, m_jtest2, m_ktest2;
   int m_nsurfpts, m_nsurfpts2;
   double m_pmin, m_pmax, m_pmin2, m_pmax2;
   // FWI workflow options
   float_sw4 m_vp_min, m_vp_max, m_vs_min, m_vs_max;  // global velocity constraints
   int m_wave_mode; // 0: P  1: S  2: both
   float m_twin_shift, m_twin_scale;
   int m_win_mode;

   MaterialParameterization *m_mp;   
   //   MaterialParCartesian *m_mpcart0;   
   std::vector<Image*> m_image_files;
   std::vector<Image3D*> m_3dimage_files;
   std::vector<SfileOutput*> m_sfiles;
   std::string m_scales_fname, m_scalem_fname;
   bool m_scales_file_given;
   std::string m_path;
   double m_reg_coeff;
   int m_nsteps_in_memory; // Number of steps when saving solution on domain boundaries for backward solve.
   int m_nstot; // dimension of m_sfs and m_xs0
   double *m_sfs; // scale factors, shared
   double *m_sfm; // scale factors, distributed
   double *m_xs0; // initial material perturbation, shared
   double *m_xm0; // initial material perturbation, distributed
   bool m_write_dfm;
};

#endif
