#ifndef SW4_EW
#define SW4_EW

#include <string>
#include <vector>
#include <iostream>

#include "sw4.h"
#include "Sarray.h"
#include "SuperGrid.h"
#include "MaterialData.h"
#include "TimeSeries.h"
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
   void processMaterialBlock( char* buffer );
   void processdGalerkin( char* buffer );
   void processReceiver( char* buffer );
   void processTopography( char* buffer );
   void defineDimensionsGXY( );
   void defineDimensionsZ();
   void allocateTopoArrays();
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
   void enforceCartTopo( vector<Sarray>& a_U );
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
   void default_bcs();
   void assign_local_bcs( );
   void create_output_directory( );
   int mkdirs(const string& path);
   bool topographyExists() {return m_topography_exists;}
   bool interpolate_topography( float_sw4 q, float_sw4 r, float_sw4& Z0, bool smoothed );
   void buildGaussianHillTopography(double amp, double Lx, double Ly, double x0, double y0);
   void compute_minmax_topography( double& topo_zmin, double& topo_zmax );
   void gettopowgh( float_sw4 ai, float_sw4 wgh[8] ) const;
   bool find_topo_zcoord_owner( float_sw4 X, float_sw4 Y, float_sw4& Ztopo );
   bool find_topo_zcoord_all( float_sw4 X, float_sw4 Y, float_sw4& Ztopo );
   void generate_grid();
   void setup_metric();
   void grid_mapping( float_sw4 q, float_sw4 r, float_sw4 s, float_sw4& x,
		      float_sw4& y, float_sw4& z );
   bool invert_grid_mapping( int g, float_sw4 x, float_sw4 y, float_sw4 z, 
			     float_sw4& q, float_sw4& r, float_sw4& s );
   int metric(  int ib, int ie, int jb, int je, int kb, int ke, float_sw4* a_x,
		 float_sw4* a_y, float_sw4* a_z, float_sw4* a_met, float_sw4* a_jac );
   int metric_rev(  int ib, int ie, int jb, int je, int kb, int ke, float_sw4* a_x,
		 float_sw4* a_y, float_sw4* a_z, float_sw4* a_met, float_sw4* a_jac );
   void metricexgh( int ib, int ie, int jb, int je, int kb, int ke,
		    int nz, float_sw4* a_x, float_sw4* a_y, float_sw4* a_z, 
		    float_sw4* a_met, float_sw4* a_jac, int order,
		    float_sw4 sb, float_sw4 zmax, float_sw4 amp, float_sw4 xc,
		    float_sw4 yc, float_sw4 xl, float_sw4 yl );
   void metricexgh_rev( int ib, int ie, int jb, int je, int kb, int ke,
			int nz, float_sw4* a_x, float_sw4* a_y, float_sw4* a_z, 
			float_sw4* a_met, float_sw4* a_jac, int order,
			float_sw4 sb, float_sw4 zmax, float_sw4 amp, float_sw4 xc,
			float_sw4 yc, float_sw4 xl, float_sw4 yl );
   void freesurfcurvisg( int ib, int ie, int jb, int je, int kb, int ke,
			 int nz, int side, float_sw4* a_u, float_sw4* a_mu,
			 float_sw4* a_la, float_sw4* a_met, float_sw4* s,
			 float_sw4* a_forcing, float_sw4* a_strx, float_sw4* a_stry );
   void freesurfcurvisg_rev( int ib, int ie, int jb, int je, int kb, int ke,
			 int nz, int side, float_sw4* a_u, float_sw4* a_mu,
			 float_sw4* a_la, float_sw4* a_met, float_sw4* s,
			 float_sw4* a_forcing, float_sw4* a_strx, float_sw4* a_stry );
   void gridinfo( int ib, int ie, int jb, int je, int kb, int ke,
		  float_sw4* a_met, float_sw4* a_jac, float_sw4&  minj,
		  float_sw4& maxj );

   void computeDT();
   void computeNearestGridPoint(int & a_i, int & a_j, int & a_k, int & a_g, float_sw4 a_x, 
				float_sw4 a_y, float_sw4 a_z);
   bool interior_point_in_proc(int a_i, int a_j, int a_g);   
   bool point_in_proc(int a_i, int a_j, int a_g);
   bool point_in_proc_ext(int a_i, int a_j, int a_g);
   void getGlobalBoundingBox(float_sw4 bbox[6]);
   bool getDepth( float_sw4 x, float_sw4 y, float_sw4 z, float_sw4 & depth );
   void computeCartesianCoord(float_sw4 &x, float_sw4 &y, float_sw4 lon, float_sw4 lat);
   void computeGeographicCoord(float_sw4 x, float_sw4 y, float_sw4 & longitude, float_sw4 & latitude);
   float_sw4 getGridAzimuth() const {return mGeoAz;}
   float_sw4 getMetersPerDegree() const {return mMetersPerDegree;}
   bool is_onesided( int g, int side ) const;
   void print_execution_time( double t1, double t2, string msg );
   void print_execution_times( double time[7] );
   void copy_supergrid_arrays_to_device();
   void copy_material_to_device();
   void setup_materials();
   void convert_material_to_mulambda();
   void extrapolateInZ( int g, Sarray& field , bool lowk, bool highk);
   void extrapolateInXY( vector<Sarray>& fields );
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
   void addsgd4cfort( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx, float_sw4* a_dcy, float_sw4* a_strx, float_sw4* a_stry, 
		      float_sw4* a_jac, float_sw4* a_cox,  float_sw4* a_coy, float_sw4 beta );
   void addsgd4cfort_indrev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			     float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
			     float_sw4* a_dcx, float_sw4* a_dcy, float_sw4* a_strx, float_sw4* a_stry, 
			     float_sw4* a_jac, float_sw4* a_cox,  float_sw4* a_coy, float_sw4 beta );

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
   void addsgd6cfort( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
		      float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
		      float_sw4* a_dcx, float_sw4* a_dcy, float_sw4* a_strx, float_sw4* a_stry, 
		      float_sw4* a_jac, float_sw4* a_cox,  float_sw4* a_coy, float_sw4 beta );
   void addsgd6cfort_indrev( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			     float_sw4* a_up, float_sw4* a_u, float_sw4* a_um, float_sw4* a_rho,
			     float_sw4* a_dcx, float_sw4* a_dcy, float_sw4* a_strx, float_sw4* a_stry, 
			     float_sw4* a_jac, float_sw4* a_cox,  float_sw4* a_coy, float_sw4 beta );
   void GetStencilCoefficients( float_sw4* _acof, float_sw4* _ghcof,
				float_sw4* _bope, float_sw4* _sbop );
   int getVerbosity() const {return 0;}
   bool usingParallelFS() const { return m_pfs;}
   int getNumberOfWritersPFS() const { return m_nwriters;}
   void get_utc( int utc[7] ) const;
   void extractRecordData(TimeSeries::receiverMode mode, int i0, int j0, int k0, int g0, 
			  vector<float_sw4> &uRec, vector<Sarray> &Um2, vector<Sarray> &U);

   // DG stuff
   int m_qu;
   int m_qv; 
   int m_nint;
   int m_single_mode_problem;
   int m_dg_bc_type;
       // Methods
   void initialData(double* udg, double* vdg);
   void assemble(double* MU,double* MV,double* SU,double* SV,double* LU,double* LV);
   void set_dg_orders( int qu, int qv);   
   void build_w_and_v(double* udg, double* vdg, double* w_in_all_faces, double* v_in_all_faces);
   void numerical_fluxes(double* w_star_all_faces, double* v_star_all_faces,
                         double* w_in_all_faces, double* v_in_all_faces,
                         double* w_out_all_faces, double* v_out_all_faces);
   void computeError(double* udg, double* vdg, double t);
   void get_exact_point_source_dG(double* u, double t, double x, double y, double z);


   enum InputMode { UNDEFINED, Efile, GaussianHill, GridFile, CartesianGrid, TopoImage, Rfile};


   // Variables ----------

   // Grid and domain
   int mNumberOfGrids, mNumberOfCartesianGrids;
   int m_ghost_points, m_ppadding, m_ext_ghost_points;
   float_sw4 m_global_xmax, m_global_ymax, m_global_zmax, m_global_zmin;
   vector<float_sw4> m_zmin;
   vector<float_sw4> mGridSize;
   vector<int> m_global_nx, m_global_ny, m_global_nz; 
   vector<int> m_iStart, m_iEnd, m_jStart, m_jEnd, m_kStart, m_kEnd;
   vector<int> m_iStartInt, m_iEndInt, m_jStartInt, m_jEndInt, m_kStartInt, m_kEndInt;
   int m_nx_base, m_ny_base, m_nz_base;
   float_sw4 m_h_base;

 // Curvilinear grid and topography
   bool m_topography_exists;
   float_sw4 m_topo_zmax;
   vector<bool> m_is_curvilinear;
   Sarray mX, mY, mZ;
   Sarray mMetric, mJ;

   Sarray mTopo, mTopoGridExt;
   InputMode m_topoInputStyle;
   string m_topoFileName;
   int m_grid_interpolation_order;
   float_sw4 m_zetaBreak;
// For some simple topographies (e.g. Gaussian hill) there is an analytical expression for the top elevation
   bool m_analytical_topo, m_use_analytical_metric;
   float_sw4 m_GaussianAmp, m_GaussianLx, m_GaussianLy, m_GaussianXc, m_GaussianYc;

// Geographic coordinate system
   float_sw4 mGeoAz, mLatOrigin, mLonOrigin, mMetersPerLongitude, mMetersPerDegree;
   bool mConstMetersPerLongitude;
   
   // MPI information and data structures
   int m_myrank, m_nprocs, m_nprocs_2d[2], m_myrank_2d[2], m_neighbor[4];
   MPI_Comm  m_cartesian_communicator;
   vector<MPI_Datatype> m_send_type1;
   vector<MPI_Datatype> m_send_type3;
   vector<MPI_Datatype> m_send_type4; // metric
   MPI_Datatype m_mpifloat;

   //   vector<MPI_Datatype> m_send_type21; // anisotropic

   // Vectors of Sarrays hold material properties on all grids. 
   vector<Sarray> mMu;
   vector<Sarray> mLambda;
   vector<Sarray> mRho;
   // Material data is the material model used to populate mMu,mLambda,mRho
   vector<MaterialData*> m_mtrlblocks;
   
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
// UTC time corresponding to simulation time 0.
   int m_utc0[7];

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
   vector<TimeSeries*> m_GlobalTimeSeries;

   // Discontinuous Galerkin stuff
   bool m_use_dg;
   
};

#endif
