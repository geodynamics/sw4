#ifndef SW4_EW
#define SW4_EW

#include <string>
#include <vector>
#include <iostream>

#include "sw4.h"
#include "Sarray.h"
#include "SuperGrid.h"

using namespace std;

class Source;
class GridPointSource;
class EWCuda;
class CheckPoint;

class EW
{
 public:
   // Methods ----------
   void timesteploop( vector<Sarray>& U, vector<Sarray>& Um);
   void timeStepLoopdGalerkin();
   EW( const string& filename );

   int computeEndGridPoint( float_sw4 maxval, float_sw4 h );
   bool startswith(const char begin[], char *line);
   void badOption(string name, char* option) const;
   void processGrid( char* buffer );
   void processTime(char* buffer);
   void processTestPointSource(char* buffer);
   void processSource( char* buffer );
   void processSuperGrid( char* buffer );
   void processDeveloper(char* buffer);
   void processFileIO( char* buffer );
   void processCheckPoint( char* buffer );
   void processRestart( char* buffer );
   void processdGalerkin( char* buffer );
   void defineDimensionsGXY( );
   void defineDimensionsZ();
   void allocateArrays();
   void printGridSizes() const;
   bool parseInputFile( const string& filename );
   void setupRun();
   bool proc_decompose_2d( int ni, int nj, int nproc, int proc_max[2] );
   void decomp1d( int nglobal, int myid, int nproc, int& s, int& e );
   void setupMPICommunications();
   bool check_for_nan( vector<Sarray>& a_U, int verbose, string name );
   bool check_for_nan_GPU( vector<Sarray>& a_U, int verbose, string name );
   bool check_for_match_on_cpu_gpu( vector<Sarray>& a_U, int verbose, string name );
   void cycleSolutionArrays(vector<Sarray> & a_Um, vector<Sarray> & a_U,
			    vector<Sarray> & a_Up ) ;
   void Force(float_sw4 a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources, bool tt );
   void evalRHS( vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		 vector<Sarray> & a_Uacc );
   void evalRHSCU( vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		   vector<Sarray> & a_Uacc, int st );
   void evalPredictor(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		      vector<Sarray>& a_Rho, vector<Sarray> & a_Lu, vector<Sarray> & a_F );
   void evalPredictorCU(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
			vector<Sarray>& a_Rho, vector<Sarray> & a_Lu, vector<Sarray> & a_F, int st );
   void evalCorrector(vector<Sarray> & a_Up, vector<Sarray>& a_Rho,
		      vector<Sarray> & a_Lu, vector<Sarray> & a_F );
   void evalCorrectorCU(vector<Sarray> & a_Up, vector<Sarray>& a_Rho,
			vector<Sarray> & a_Lu, vector<Sarray> & a_F, int st );
   void evalDpDmInTime(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		       vector<Sarray> & a_Uacc );
   void evalDpDmInTimeCU(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
			 vector<Sarray> & a_Uacc, int st );
   void communicate_array( Sarray& U, int g );
 
   void cartesian_bc_forcing( float_sw4 t, vector<float_sw4**> & a_BCForcing,
			      vector<Source*>& a_sources );
   void setup_boundary_arrays();
   void side_plane( int g, int side, int wind[6], int nGhost );   
   void enforceBC( vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		   float_sw4 t, vector<float_sw4**> & a_BCForcing );
   void addSuperGridDamping(vector<Sarray> & a_Up, vector<Sarray> & a_U,
			    vector<Sarray> & a_Um, vector<Sarray> & a_Rho );
   void addSuperGridDampingCU(vector<Sarray> & a_Up, vector<Sarray> & a_U,
			      vector<Sarray> & a_Um, vector<Sarray> & a_Rho, int st );
   void printTime( int cycle, float_sw4 t, bool force ) const;
   bool exactSol(float_sw4 a_t, vector<Sarray> & a_U, vector<Source*>& sources );

   float_sw4 SmoothWave(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 VerySmoothBump(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 C6SmoothBump(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 Gaussian(float_sw4 t, float_sw4 R, float_sw4 c, float_sw4 f );
   float_sw4 d_SmoothWave_dt(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 d_VerySmoothBump_dt(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 d_C6SmoothBump_dt(float_sw4 t, float_sw4 R, float_sw4 c);
   float_sw4 d_Gaussian_dt(float_sw4 t, float_sw4 R, float_sw4 c, float_sw4 f);
   float_sw4 SWTP(float_sw4 Lim, float_sw4 t);   
   float_sw4 VSBTP(float_sw4 Lim, float_sw4 t);
   float_sw4 C6SBTP(float_sw4 Lim, float_sw4 t);
   float_sw4 SmoothWave_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta);
   float_sw4 VerySmoothBump_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta);
   float_sw4 C6SmoothBump_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 alpha, float_sw4 beta);
   float_sw4 Gaussian_x_T_Integral(float_sw4 t, float_sw4 R, float_sw4 f, float_sw4 alpha, float_sw4 beta);
   void get_exact_point_source( float_sw4* up, float_sw4 t, int g, Source& source,
				int* wind=NULL );
   void normOfDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, float_sw4 &diffInf, 
			  float_sw4 &diffL2, float_sw4 &xInf, vector<Source*>& a_globalSources );
   void setupSBPCoeff();
   void check_dimensions();
   void setup_supergrid( );
   void assign_supergrid_damping_arrays();
   void assign_local_bcs( );
   void create_output_directory( );
   int mkdirs(const string& path);
   void computeDT();
   void computeNearestGridPoint(int & a_i, int & a_j, int & a_k, int & a_g, float_sw4 a_x, 
				float_sw4 a_y, float_sw4 a_z);
   bool interior_point_in_proc(int a_i, int a_j, int a_g);   
   bool is_onesided( int g, int side ) const;
   void print_execution_time( float_sw4 t1, float_sw4 t2, string msg );
   void print_execution_times( float_sw4 time[7] );
   void copy_supergrid_arrays_to_device();
   void copy_material_to_device();
   void find_cuda_device();
   //   void reset_gpu();
   void corrfort( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* up,
		  float_sw4* lu, float_sw4* fo, float_sw4* rho, float_sw4 dt4 );
   void predfort( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* up,
		   float_sw4* u, float_sw4* um, float_sw4* lu, float_sw4* fo,
		  float_sw4* rho, float_sw4 dt2 );
   void dpdmtfort( int ib, int ie, int jb, int je, int kb, int ke, float_sw4* up,
		   float_sw4* u, float_sw4* um, float_sw4* u2, float_sw4 dt2i );
   void solerr3fort( int ib, int ie, int jb, int je, int kb, int ke,
		      float_sw4 h, float_sw4* uex, float_sw4* u, float_sw4& li,
		      float_sw4& l2, float_sw4& xli, float_sw4 zmin, float_sw4 x0,
		      float_sw4 y0, float_sw4 z0, float_sw4 radius,
		     int imin, int imax, int jmin, int jmax, int kmin, int kmax );
   void bcfortsg( int ib, int ie, int jb, int je, int kb, int ke, int wind[36], 
		   int nx, int ny, int nz, float_sw4* u, float_sw4 h, boundaryConditionType bccnd[6],
		   float_sw4 sbop[5], float_sw4* mu, float_sw4* la, float_sw4 t,
		   float_sw4* bforce1, float_sw4* bforce2, float_sw4* bforce3, 
		   float_sw4* bforce4, float_sw4* bforce5, float_sw4* bforce6,
		   float_sw4 om, float_sw4 ph, float_sw4 cv,
		  float_sw4* strx, float_sw4* stry );
   void bcfortsg_indrev( int ib, int ie, int jb, int je, int kb, int ke, int wind[36], 
		   int nx, int ny, int nz, float_sw4* u, float_sw4 h, boundaryConditionType bccnd[6],
		   float_sw4 sbop[5], float_sw4* mu, float_sw4* la, float_sw4 t,
		   float_sw4* bforce1, float_sw4* bforce2, float_sw4* bforce3, 
		   float_sw4* bforce4, float_sw4* bforce5, float_sw4* bforce6,
		   float_sw4 om, float_sw4 ph, float_sw4 cv,
		  float_sw4* strx, float_sw4* stry );
   void addsgd4fort( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		      float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		      float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		     float_sw4 beta );
   void addsgd4fort_indrev( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		      float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		      float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		     float_sw4 beta );
   void addsgd6fort( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		      float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		      float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		      float_sw4 beta );
   void addsgd6fort_indrev( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx,  float_sw4* a_dcy,  float_sw4* a_dcz,
		      float_sw4* a_strx, float_sw4* a_stry, float_sw4* a_strz,
		      float_sw4* a_cox,  float_sw4* a_coy,  float_sw4* a_coz,
		      float_sw4 beta );
   void GetStencilCoefficients( float_sw4* _acof, float_sw4* _ghcof,
				float_sw4* _bope, float_sw4* _sbop );

   bool usingParallelFS(){ return m_pfs;}
   int getNumberOfWritersPFS(){ return m_nwriters;}
   // DG stuff
   void set_dg_orders( int qu, int qv);
   int m_qu;
   int m_qv; 

   // Variables ----------

   // Grid and domain
   int mNumberOfGrids, mNumberOfCartesianGrids;
   int m_ghost_points, m_ppadding;
   float_sw4 m_global_xmax, m_global_ymax, m_global_zmax, m_global_zmin;
   vector<float_sw4> m_zmin;
   vector<float_sw4> mGridSize;
   vector<int> m_global_nx, m_global_ny, m_global_nz; 
   vector<int> m_iStart, m_iEnd, m_jStart, m_jEnd, m_kStart, m_kEnd;
   vector<int> m_iStartInt, m_iEndInt, m_jStartInt, m_jEndInt, m_kStartInt, m_kEndInt;
   int m_nx_base, m_ny_base, m_nz_base;
   float_sw4 m_h_base;
   bool m_topography_exists;
   float_sw4 m_topo_zmax;
   vector<bool> m_is_curvilinear;
   
   // MPI information and data structures
   int m_myrank, m_nprocs, m_nprocs_2d[2], m_myrank_2d[2], m_neighbor[4];
   MPI_Comm  m_cartesian_communicator;
   vector<MPI_Datatype> m_send_type1;
   vector<MPI_Datatype> m_send_type3;
   //   vector<MPI_Datatype> m_send_type4; // metric
   //   vector<MPI_Datatype> m_send_type21; // anisotropic

   // Vectors of Sarrays hold material properties on all grids. 
   vector<Sarray> mMu;
   vector<Sarray> mLambda;
   vector<Sarray> mRho;

   // Vectors of solution at time t_n and t_{n-1}
   vector<Sarray> mU;
   vector<Sarray> mUm;

   // SBP boundary operator coefficients and info
   float_sw4 m_iop[5], m_iop2[5], m_bop2[24], m_sbop[5], m_acof[384], m_bop[24];
   float_sw4 m_bope[48], m_ghcof[6], m_hnorm[4];
   vector<int*> m_onesided; 

   // Time stepping variables
   float_sw4 mCFL, mTstart, mTmax, mDt;
   int mNumberOfTimeSteps;
   bool mTimeIsSet;

   // Storage for supergrid damping coefficients (1-D)
   vector<float_sw4*> m_sg_dc_x, m_sg_dc_y, m_sg_dc_z;
   vector<float_sw4*> dev_sg_dc_x, dev_sg_dc_y, dev_sg_dc_z;
   vector<float_sw4*> m_sg_str_x, m_sg_str_y, m_sg_str_z;
   vector<float_sw4*> dev_sg_str_x, dev_sg_str_y, dev_sg_str_z;
   vector<float_sw4*> m_sg_corner_x, m_sg_corner_y, m_sg_corner_z;
   vector<float_sw4*> dev_sg_corner_x, dev_sg_corner_y, dev_sg_corner_z;
   SuperGrid m_supergrid_taper_x, m_supergrid_taper_y;
   vector<SuperGrid> m_supergrid_taper_z;

   // Boundary conditions
   boundaryConditionType mbcGlobalType[6]; 
   vector<boundaryConditionType*> m_bcType;
   vector<int *> m_NumberOfBCPoints;
   vector<int *> m_BndryWindow;

   // Test modes
   bool m_point_source_test, m_moment_test;

   // diagnostic output, error checking
   int mPrintInterval, mVerbose;
   bool mQuiet;
   bool m_checkfornan, m_output_detailed_timing;
   string mPath;

   // File io
   bool m_pfs;
   int m_nwriters;

   // Sources
   vector<Source*> m_globalUniqueSources;
   vector<GridPointSource*> m_point_sources;

   // Supergrid boundary conditions
   float_sw4 m_supergrid_damping_coefficient;
   int m_sg_damping_order, m_sg_gp_thickness;
   bool m_use_supergrid;

   // GPU computing
   int m_ndevice;
   int m_gpu_gridsize[3], m_gpu_blocksize[3];
   EWCuda* m_cuobj;

   bool m_corder; // (i,j,k,c) order 

   // Output: Images, stations, checkpoints
   vector<CheckPoint*> m_check_points;
   CheckPoint* m_restart_check_point;

   // Discontinuous Galerkin stuff
   bool m_use_dg;
   
};

#endif
