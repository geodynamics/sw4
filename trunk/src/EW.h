//-*-c++-*-
#ifndef EW_H
#define EW_H

#include <mpi.h>

#include <string>
#include <vector>
#include <fstream>
#include <list>

#include "Sarray.h"
#include "TimeSeries.h"
#include "Source.h"
#include "GridPointSource.h"
#include "Filter.h"

#include "Image.h"
// #include "Image3D.h"

#include "boundaryConditionTypes.h"
#include "ForcingTwilight.h"
#include "TestPointSource.h"
#include "TestEnergy.h"
#include "TestLamb.h"
#include "TestRayleighWave.h"

#include "MaterialData.h"
#include "EtreeFile.h"

#include "SuperGrid.h"
#include "MaterialProperty.h"

using namespace std;

class EW 
{
public:
EW(const string& name, vector<Source*> & a_GlobalUniqueSources, 
   vector<TimeSeries*> & a_GlobalTimeSeries, bool invproblem=false );
~EW();
bool wasParsingSuccessful();
bool isInitialized();

void set_output_options( bool output_load, bool output_detailed_timing );
void setGMTOutput(string filename, string wppfilename);
void saveGMTFile( vector<Source*> & a_GlobalUniqueSources );
void allocateCartesianSolverArrays(double a_global_zmax);
void setGoalTime(double t);
//double getCurrentTime(){return mTime;}

void setAttenuationParams(int numberOfMechanisms, double velocityOmega, int ppw, double maxfrequency );

void setNumberSteps(int steps); // remove???
int getNumberOfSteps() const;

void setupRun( );
void preprocessSources( vector<Source*> & a_GlobalSources );

void solve( vector<Source*> & a_GlobalSources, vector<TimeSeries*> & a_GlobalTimeSeries );
   void solve_backward( vector<Source*> & a_Sources, vector<TimeSeries*> & a_TimeSeries, double gradient[11], double hessian[121] );

bool parseInputFile( vector<Source*> & a_GlobalSources, vector<TimeSeries*> & a_GlobalTimeSeries );
void parsedate( char* datestr, int& year, int& month, int& day, int& hour, int& minute,
		int& second, int& msecond, int& fail );

void extractRecordData(TimeSeries::receiverMode mode, int i0, int j0, int k0, int grid0, 
		       vector<double> &uRec, vector<Sarray> &Um2, vector<Sarray> &U);

// some (all?) of these functions are called from parseInputFile() and should be made private
void badOption(string name, char* option) const;
void processGrid(char* buffer);
void deprecatedOption(const string& command, 
		      const string& oldone, 
		      const string& newone);
void processTime(char* buffer);
void processTwilight(char* buffer);
void processFileIO(char* buffer);
void processImage(char* buffer);
void deprecatedImageMode(int value, const char* name) const;
void processTestPointSource(char* buffer);
void processTestRayleigh(char* buffer);
void processTestLamb(char* buffer);
void processTestEnergy(char* buffer);
void processSource(char* buffer, vector<Source*> & a_GlobalUniqueSources);
void processMaterialBlock( char* buffer, int & blockCount );
void processMaterialPfile(char* buffer);
void processMaterialEtree(char* buffer);
void processReceiver(char* buffer, vector<TimeSeries*> & a_GlobalTimeSeries);
void processObservation(char* buffer, vector<TimeSeries*> & a_GlobalTimeSeries);
void processBoundaryConditions(char *buffer);
void processPrefilter(char* buffer);
void processGMT(char* buffer);
void processDeveloper(char* buffer);
void processGlobalMaterial(char* buffer);

void side_plane( int g, int side, int wind[6], int nGhost );
void setPrintCycle(int cycle) { mPrintInterval = cycle; }
void setVerbosity(int level) { mVerbose = level; };
int  getVerbosity() {return mVerbose; };
int  getRank() {return m_myRank; };
void setDebugIO(bool onoff) { mDebugIO = onoff; }
  
//void setDampingCFL(double d4_cfl) { m_d4_cfl = d4_cfl; }

void printTime(int cycle, double t, bool force=false ) const;
void printPreamble(vector<Source*> & a_Sources) const;
void switch_on_checkfornan();
void switch_on_error_log();
void set_energylog( string logfile, bool print, bool elog );
void set_inner_loop( int loopnr );
void set_cflnumber( double cfl );
void set_testing_mode(bool a_testing){m_testing = a_testing;}
bool get_testing_mode(){return m_testing;}

void default_bcs( );

void set_twilight_forcing( ForcingTwilight* a_forcing );
// perhaps these functions should be in the ForcingTwilight class? 
// but how will they get access to the material properties and grid sizes?
void initialData(double a_t, vector<Sarray> & a_U, vector<Sarray*> & a_AlphaVE);
bool exactSol(double a_t, vector<Sarray> & a_U, vector<Sarray*> & a_AlphaVE, vector<Source*>& source );
void exactRhsTwilight(double a_t, vector<Sarray> & a_F);
void exactAccTwilight(double a_t, vector<Sarray> & a_Uacc);
void Force(double a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources );
void Force_tt(double a_t, vector<Sarray> & a_F, vector<GridPointSource*> point_sources );

void normOfDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, double &diffL2, double &xInf,
		       vector<Source*>& a_globalSources );
void normOfDifferenceGhostPoints( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, double &diffL2 );
void normOfSurfaceDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, 
			      double &diffL2, double &solInf, double &solL2, vector<Source*> & a_globalSources);

void test_sources( vector<GridPointSource*>& a_point_sources, vector<Source*>& a_global_unique_sources,
		   vector<Sarray>& F );
void testSourceDiscretization( int kx[3], int ky[3], int kz[3],
			       double moments[3], vector<GridPointSource*>& point_sources, vector<Sarray>& F );

void setupSBPCoeff( );

// time stepping routines
void enforceBC( vector<Sarray> & a_U, double t, vector<double **> & a_BCForcing );
void cartesian_bc_forcing( double t, vector<double **> & a_BCForcing );
void evalRHS(vector<Sarray> & a_U, vector<Sarray> & a_Lu );
void evalPredictor(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		   vector<Sarray> & a_Lu, vector<Sarray> & a_F );
void evalDpDmInTime(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		    vector<Sarray> & a_Uacc );
void evalCorrector(vector<Sarray> & a_Up, vector<Sarray> & a_Lu, vector<Sarray> & a_F );
void addSuperGridDamping(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um );
void cycleSolutionArrays(vector<Sarray> & a_Um, vector<Sarray> & a_U, vector<Sarray> & a_Up, 
			 vector<Sarray*> & a_AlphaVEm, vector<Sarray*> & a_AlphaVE, vector<Sarray*> & a_AlphaVEp);
void cycleSolutionArrays(vector<Sarray> & a_Um, vector<Sarray> & a_U, vector<Sarray> & a_Up ); 

void bndryInteriorDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, 
			      double lowZ[3], double interiorZ[3], double highZ[3] );
void test_RhoUtt_Lu( vector<Sarray> & a_Uacc, vector<Sarray> & a_Lu, vector<Sarray> & a_F, 
		     double lowZ[3], double interiorZ[3], double highZ[3] );

void setRestartInfo(int fromCycle, int dumpInterval, const string& filePrefix);
void computeDT();
   //bool inTestSourceMode() { return mTestSource; }
   //bool inTestLambMode() { return mTestLamb; }
bool proc_zero() const;
int no_of_procs() const;
void create_output_directory();
void initialize_image_files();
void update_images( int Nsteps, double time, vector<Sarray> & a_Up, vector<Sarray>& a_U, 
		    vector<Sarray>& a_Um, vector<Source*> & a_sources, int dminus );

void initialize_SAC_files(); // going away
void update_SACs( int Nsteps ); // going away

void print_execution_times( double times[7] );
void print_execution_time( double t1, double t2, string msg );
void finalizeIO();
string bc_name( const boundaryConditionType bc ) const;
int mkdirs(const string& path);
void setOutputPath(const string& path);
const string& getOutputPath() { return mPath; };
const string& getObservationPath() { return mObsPath; };
const string& getName() { return mName; };
void set_global_bcs(boundaryConditionType bct[6]); // assigns the global boundary conditions

void add_mtrl_block( MaterialData* md ){ m_mtrlblocks.push_back( md ); };


void set_threshold_velocities(double vpmin, double vsmin);

// material properties by id
inline bool inside_material_surfaces( double lat, double lon )
    {
      return (lat <= m_materialLatMax && lat >= m_materialLatMin && 
	      lon <= m_materialLonMax && lon >= m_materialLonMin);
    }

void addMaterialProperty(MaterialProperty* mat){m_materials.push_back(mat);}


void getMaterialID(double lat, double lon, double depth, int &materialID);
bool knownMaterial(int materialID);
double lookup_Rho(int materialID, double depth);
double lookup_Vs(int materialID, double depth);
double lookup_Vp(int materialID, double depth);

// attenuation model
double lookup_Qp(int materialID, double depth);
double lookup_Qs(int materialID, double depth);

// super-grid functions
void processSupergrid(char *buffer);
void set_sg_damping(double coeff);
void set_sg_thickness(int gp_thickness);
void set_sg_transition(int gp_trans);
bool usingSupergrid(){return m_use_supergrid;};
void setup_supergrid( );
   //void supergrid_taper_material();
void assign_supergrid_damping_arrays();

void assign_local_bcs( );
bool timeSteppingSet();
bool proc_decompose_2d( int ni, int nj, int nproc, int proc_max[2] );
void decomp1d( int nglobal, int myid, int nproc, int& s, int& e );
void coarsen1d( int& n, int& ifirst, int& ilast );
void allocateCurvilinearArrays();
void generate_grid();
void setup_metric();

void setupMPICommunications();
void setup2D_MPICommunications();
void communicate_array( Sarray& u, int grid );
void communicate_arrays( vector<Sarray>& u );
void communicate_array_2dfinest( Sarray& u );
void communicate_array_2d( Sarray& u, int g, int k );
void communicate_array_2d_asym( Sarray& u, int g, int k );

void set_materials();
void setup_viscoelastic(double minvsoh );
void extrapolateInZ(Sarray& field, bool useThreshold, double thresHoldValue, bool linear);

void addImage(Image* i);
//void addImage3D(Image3D* i);
void setIO_timing(bool iotiming);
void setParallel_IO(bool pfs, int nwriters);

void extractTopographyFromGridFile(string a_topoFileName);
void extractTopographyFromCartesianFile(string a_topoFileName);

void setEtreeFile(EtreeFile* efile); 
void extractTopographyFromEfile(string a_topoFileName, string a_topoExtFileName, string a_QueryType,
                                double a_EFileResolution);
void smoothTopography(int maxIter);

void buildGaussianHillTopography(double amp, double Lx, double Ly, double x0, double y0);

void extractSurfaceFromGridFile(string a_surfaceFileName);

void computeCartesianCoord(double &x, double &y, double lon, double lat);
void computeGeographicCoord(double x, double y, double & longitude, double & latitude);

void initializeSystemTime();
void set_epicenter_in_SAC(); 
   
// void update_all_boundaries(vector<Sarray> &U, vector<Sarray> &UM, double t,
// 			   vector<Sarray*> &AlphaVE );

// void impose_physical_bc(vector<Sarray> &U, vector<Sarray> &UM, double t,
// 			vector<Sarray*> &AlphaVE );

/* void bc_dirichlet( Sarray& u, int g, double t, int side, */
/*  		   Forcing* forcing, double h ); */

/* void bc_free_surface( Sarray& u, int g, double t, int side, */
/* 		      Forcing* forcing, double h, int onesided[6] ); */

void computeNearestGridPoint(int & a_i, 
			     int & a_j, 
			     int & a_k,
			     int & a_g, // grid on which indices are located
			     double a_x, 
			     double a_y, 
			     double a_z);
  

void computeNearestSurfaceGridPoint(int & a_i, 
                                    int & a_j, 
                                    double a_x, 
                                    double a_y, 
                                    double a_z);
  
void coord_from_index( int i, int j, int k, int g, double& x, double& y, double& z );

double distance(double a_x1, double a_y1, double a_z1,
                double a_x0, double a_y0, double a_z0);

void computeNearestLowGridPoint(int & a_i, 
                                int & a_j, 
                                int & a_k,
                                int & a_g, // grid on which indices are located
                                double a_x, 
                                double a_y, 
                                double a_z);
  
bool interior_point_in_proc(int a_i, int a_j, int a_g); // only takes interior points into account
bool point_in_proc(int a_i, int a_j, int a_g);          // both interior and parallel ghost points
  
void initializePaddingCells();

void convert_material_to_mulambda();

void check_materials();

void computeSolutionError(vector<Sarray> &U, double t, vector<Sarray*> &Alpha );

double localMin(vector<Sarray> & a_field);
double localMax(vector<Sarray> & a_field);
double localMinVp();
double localMaxVp(); 
double localMinVs(); 
double localMaxVs(); 
double localMinVpOverVs();
double localMaxVpOverVs(); 

bool topographyExists(){return m_topography_exists;};
bool usingAttenuation(){return m_use_attenuation;};

void interpolate_between_grids( vector<Sarray>& u, vector<Sarray>& um, double t, 
  			        vector<Sarray*> &AlphaVE );

bool interpolate_topography( double q, double r, double & Z0, bool smoothed);

bool getDepth( double x, double y, double z, double & depth);

bool curvilinear_grid_mapping( double q, double r, double s, double & X0, double & Y0, double & Z0 );

bool invert_curvilinear_grid_mapping( double X0, double Y0, double Z0, double& q, double& r, double& s );

bool find_curvilinear_derivatives_at_point( double q, double r, double s,
					    double qX[], double rX[], double sX[]);
 
void save_errors( double max_error[3], double l2_error[3] );

void compute_minvsoverh( double& minvsoh );

void set_resolution( int ppw );

void set_prefilter( FilterType passband, int order, int passes, double fc1, double fc2 );

void set_scenario(const string& scenario );

void set_conservative_interpolation( bool onoff, double ctol, int cmaxit );

void set_geodyn_data( string filename, int nx, int nz, double h, double origin[3],
		      double dt, int nsteps, int faces );

void impose_geodyn_ibcdata( vector<Sarray> &u, vector<Sarray> &um, double t );

void get_geodyn_timelevel( vector<Sarray>& geodyndata );

void copy_geodyn_timelevel( vector<Sarray>& geodyndata1,
			    vector<Sarray>& geodyndata2 );

void consintp( Sarray& u_a, Sarray& um_a, Sarray& f_a, Sarray& mu_a, Sarray& la_a, Sarray& rho_a,
	       Sarray& uf_a, Sarray& ufm_a, Sarray& ff_a, Sarray& muf_a, Sarray& laf_a,
	       Sarray& rhof_a, Sarray* AlphaVE, Sarray* AlphaVEf, double hc, double hf, double dt, int g,
	       double* a1, int* ipiv1, double* a2, int* ipiv2, int bctype[4], double tp1 );
void check_consintp( Sarray& uc_a, Sarray& uf_a, Sarray* alphac_a, Sarray* alphaf_a );

void integrate_source( );

   void compute_energy( double dt, bool write_file, vector<Sarray>& Um,
			vector<Sarray>& U, vector<Sarray>& Up, int step );

//  void update_maxes_hVelMax();
//  void update_maxes_vVelMax();

//   void computeDivergence()   ;
//   void computeCurlMagnitude();

//   void setTestPointSourceMode(){ mTestSource = true; };

// functions from the old FileInput class
void cleanUpRefinementLevels();

enum InputMode { UNDEFINED, Efile, GaussianHill, GridFile, CartesianGrid};

// access functions needed by the Image (and perhaps other) classes
int getNumberOfCartesianGrids(){return mNumberOfCartesianGrids;};
int getNumberOfGrids(){return mNumberOfGrids;};
int getNumberOfGhostPoints(){return m_ghost_points;};
int getNumberOfParallelPaddingPoints(){return m_ppadding;};
double getLatOrigin(){ return mLatOrigin;};
double getGridAzimuth(){ return mGeoAz;};
double getMetersPerDegree(){ return mMetersPerDegree;};
bool usingParallelFS(){ return m_pfs;};
int getNumberOfWritersPFS(){ return m_nwriters;};

 // test point source
void get_exact_point_source( Sarray& u, double t, int g, Source& source );
double VerySmoothBump_x_T_Integral(double t, double R, double alpha, double beta);
double C6SmoothBump_x_T_Integral(double t, double R, double alpha, double beta);
double SmoothWave_x_T_Integral(double t, double R, double alpha, double beta);
double Gaussian_x_T_Integral(double t, double R, double f, double alpha, double beta);
double VSBTP(double Lim, double t);
double C6SBTP(double Lim, double t);
double SWTP(double Lim, double t);
double d_VerySmoothBump_dt(double t, double R, double c);
double d_C6SmoothBump_dt(double t, double R, double c);
double d_SmoothWave_dt(double t, double R, double c);
double d_Gaussian_dt(double t, double R, double c, double f);
double VerySmoothBump(double t, double R, double c);
double C6SmoothBump(double t, double R, double c);
double SmoothWave(double t, double R, double c);
double Gaussian(double t, double R, double c,double f);

// Lamb's problem
void get_exact_lamb( vector<Sarray> & a_U, double a_t, Source& a_source );
void get_exact_lamb2( vector<Sarray> & a_U, double a_t, Source& a_source );
double G4_Integral(double T, double t, double r, double beta);
double G3_Integral(double iT, double it, double ir, double ibeta);
double G2_Integral(double iT, double it, double ir, double ibeta);


void getGlobalBoundingBox(double bbox[6]);

string getPath(){ return mPath; }
void set_utcref( TimeSeries& ts );
void print_utc();

   // For inverse problem
void processCG(char* buffer );
void processScaleFactors(char* buffer );
void average_speeds( double& cp, double& cs );
void layered_speeds( vector<double>& cp, vector<double>& z );
void testsourcediff( vector<Source*> GlobalSources, double gradient[11], double hessian[121] );
void get_scalefactors( double sf[11] ); 
bool compute_sf();
void compute_guess( bool& guesspos, bool& guesst0fr, bool& guessmom, bool& output_seismograms );
void get_cgparameters( int& maxit, int& maxrestart, double& tolerance, bool& fletcherreeves,
		       int& stepselection, bool& do_linesearch, int& varcase );

//
// VARIABLES BEYOND THIS POINT
//

// ------------------------------------------
// Grid 
// ------------------------------------------

int mNumberOfGrids, mNumberOfCartesianGrids;

// grid sizes are needed by the Source and Image classes, so should be kept public
vector<double> mGridSize;

// part of global array on each processor, including ghost points = all points
vector<int> m_iStart, m_iEnd, m_jStart, m_jEnd, m_kStart, m_kEnd; 

// global number of grid points on each refinement level, without ghost points
vector<int> m_global_nx, m_global_ny, m_global_nz; 

// part of global array on each processor, excluding ghost points and parallel overlap points = interior points
vector<int> m_iStartInt, m_iEndInt, m_jStartInt, m_jEndInt, m_kStartInt, m_kEndInt; 

// Note that the m_paddingCells array is no longer needed to get the range of internal grid points 
// Instead use m_iStartInt[g], m_iEndInt[g], etc, 
int m_paddingCells[4]; // indexing is [0] = low-i, [1] = high-i, [2] = low-j, [3] = high-j

// For the Cartesian grid, we only need to offset in z
vector<double> m_zmin; // needed by the Source and Image classes

// for the curvilinear grid, we also store the cartesian coordinates of the grid points
Sarray mX, mY, mZ; // needed by the Source class, so must be public
Sarray mJ; // Jacobian also needed by the Source class

// command prefilter
bool m_prefilter_sources, m_filter_observations;
// filter setup
// Filter for time function
Filter *m_filter_ptr;
// Filter for observations
Filter *m_filterobs_ptr;
  // Test cases for optimizer, validate gradient, hessian, output function surface, etc...
int m_opttest;

// 2-D arrays with elevation-values (=-z) as function of horizontal indices
// mTopo holds the raw topography (according to the "elevation" field in the etree)
// mTopoMat holds the highest elevation where the etree returns solid material properties
// mTopoGrid holds the smoothed topography which follows the top surface of the curvilinear grid
Sarray mTopo, mTopoMat, mTopoGrid;

private:

ForcingTwilight* m_twilight_forcing;
TestPointSource* m_point_source_test;
bool m_moment_test;
TestEnergy* m_energy_test;
TestLamb* m_lamb_test;
TestRayleighWave* m_rayleigh_wave_test;

vector<MaterialData*> m_mtrlblocks;
// index convention: [0]: low-x, [1]: high-x, [2]: low-y, [3]: high-y; [4]: low-z, [5]: high-z  
boundaryConditionType mbcGlobalType[6]; // these are the boundary conditions for the global problem
vector<boundaryConditionType*> m_bcType;  // these are the boundary conditions for each grid on the local processor, with bProcessor conditions
double mTstart;
double mDt;
EtreeFile * mEtreeFile;

bool m_doubly_periodic;
MPI_Comm m_cartesian_communicator;
int m_proc_array[2];

bool mbcsSet;


// for some simple topographies (e.g. Gaussian hill) there is an analytical expression for the top elevation
bool m_analytical_topo;
double m_GaussianAmp, m_GaussianLx, m_GaussianLy, m_GaussianXc, m_GaussianYc;

// interface surfaces in the material model
int m_number_material_surfaces, m_Nlon, m_Nlat;
double m_materialLonMax, m_materialLonMin, m_materialLatMax, m_materialLatMin;
Sarray m_materialDepth;
double *m_materialLon, *m_materialLat;

// global material thresholds
bool m_useVelocityThresholds;
double m_vpMin, m_vsMin;

// material description used with material surfaces and the ifile command
vector<MaterialProperty*> m_materials;

// order of polynomial mapping in algebraic grid genenerator
int m_grid_interpolation_order;

// metric of the curvilinear grid
double m_minJacobian, m_maxJacobian;

string m_scenario;

// command limitfrequency
double m_frequency_limit;
bool m_limit_frequency;
int m_ppw;


// parallel io stuff
bool m_pfs;
int m_nwriters;

// supergrid
bool m_use_supergrid;
int m_sg_gp_thickness, m_sg_gp_transition;
double m_supergrid_damping_coefficient;
SuperGrid m_supergrid_taper_x, m_supergrid_taper_y, m_supergrid_taper_z;

string mPath, mObsPath;

// number of boundary points on each side
vector<int *> m_NumberOfBCPoints;

// ghost point index window for each side of the boundary on each grid
vector<int *> m_BndryWindow;

// vectors of Sarrays hold material properties on all grids. 
vector<Sarray> mMu;
vector<Sarray> mLambda;
vector<Sarray> mRho;

// attenuation variables (only allocated if attenuation is enabled)
bool m_use_attenuation, m_att_use_max_frequency;
int m_number_mechanisms, m_att_ppw;
double m_velo_omega, m_min_omega, m_max_omega, m_att_max_frequency;

vector<Sarray> mQp, mQs;
vector<Sarray*> mMuVE, mLambdaVE;
// relaxation frequencies
vector<double> mOmegaVE;

// and the metric derivatives as well as the jacobian
Sarray mQ, mR, mS;

// Vectors of pointers to hold boundary forcing arrays in each grid
// this is innner cube data for coupling with other codes
// bool m_do_geodynbc;
// vector<int*> m_geodyn_dims;
// vector<Sarray> m_geodyn_data1;
// vector<Sarray> m_geodyn_data2;
// double m_geodyn_origin[3], m_geodyn_h, m_geodyn_dt;
// int m_geodyn_step, m_geodyn_maxsteps, m_geodyn_blocksize;
//    int m_geodyn_ni, m_geodyn_nj, m_geodyn_nk, m_geodyn_faces;
// string m_geodyn_filename;
// ifstream m_geodynfile;
// bool m_geodyn_iwillread;   

// with topo, zmin might be different from 0
double m_global_xmax, m_global_ymax, m_global_zmin, m_global_zmax; 

// number of grid points near mesh refinement boundary, for  extrapolating material properties
int mMaterialExtrapolate; 

// variables from the old FileInput class
int m_nx_base, m_ny_base, m_nz_base;
double m_h_base;
vector<double> m_refinementBoundaries;
InputMode m_topoInputStyle;
string m_topoFileName;
bool mTopoImageFound;
double m_topo_zmax;

//-------------------------------------------
// IO data
//-------------------------------------------
int m_myRank, m_nProcs;

string mName;
//string mWPPFileName;
string mGMTFileName;


bool mWriteGMTOutput;
int mPlotFrequency;
int mNumFiles;
int mVerbose;
bool mDebugIO;
bool mHomogeneous;

// SAC file info (is this used anymore?)
bool mCompareSACFiles;
float mSACFileErrorTolerance;

// Image file info
vector<Image*> mImageFiles;
//vector<Image3D*> mImage3DFiles;
bool m_iotiming;

// time data
bool mTimeIsSet;
double mTmax;

int mNumberOfTimeSteps;

// Test modes
int m_update_boundary_function;

   //bool mTestSource;
   //bool mTestLamb;
   //bool mTestingEnergy;
int mOrder;
double mCFL;

vector<int*> m_onesided; 
double m_curlcoeff, m_d4coeff, m_d4_cfl; // these should go away

// storage for the 1-D damping coefficients
vector<double*> m_sg_dc_x, m_sg_dc_y, m_sg_dc_z;
vector<double*> m_sg_str_x, m_sg_str_y, m_sg_str_z;

//-------------------------------------------
// restart data
//-------------------------------------------
// string mRestartFilePrefix;
// int mRestartFromCycle;
// int mRestartDumpInterval;

//----------------------------------------
// Energy test data
//----------------------------------------
bool m_energy_log, m_energy_print;
double m_saved_energy;
string m_energy_logfile;
vector<double> m_energy; // *

//-------------------------------------------
// Measure wall clock time variables
//-------------------------------------------
bool m_do_timing;
int  m_timing_print;
bool m_output_detailed_timing;
bool m_output_load;

int m_projection_cycle;

bool m_checkfornan;
  
// testing
double m_max_error[3], m_l2_error[3];

string m_error_log_file;
bool m_error_log, m_error_print;
int m_inner_loop;

//  Conservative interface
bool m_intp_conservative;
bool m_matrices_decomposed;
double m_citol;
int m_cimaxiter;

vector<double*> m_cimat1;
vector<double*> m_cimat2;
vector<int*> m_ciipiv1;
vector<int*> m_ciipiv2;

EW(const EW&);
EW& operator=(const EW&);

int mPrintInterval;
// (lon, lat) origin of Grid as well as
double mGeoAz;
double mLonOrigin, mLatOrigin;

//GeographicCoord mGeoCoord;
double mMetersPerDegree;

// is this object ready for time-stepping?
bool mParsingSuccessful, mIsInitialized, mSourcesOK;
bool m_testing;

   // Parameters related to the inverse problem   
bool m_inverse_problem; // Will we solve the inverse problem?
bool m_iniguess_pos, m_iniguess_t0fr, m_iniguess_mom;// Estimate initial guess ?
bool m_output_initial_seismograms;
bool m_compute_scalefactors;
int m_maxit,m_maxrestart;
double m_tolerance;
double m_scalefactors[11];   
int m_cgstepselection, m_cgvarcase;
bool m_cgfletcherreeves, m_do_linesearch;

// Number of grid points per wave length, P = min Vs/(f*h) 
vector<double> mMinVsOverH;
int m_ghost_points;
int m_ppadding;

// coefficients for boundary modified 4th order SBP operators
double m_iop[5], m_iop2[5], m_bop2[24], m_sbop[5], m_acof[384], m_bop[24];
double m_bope[48], m_ghcof[6], m_hnorm[4];

int m_neighbor[4];
vector<MPI_Datatype> m_send_type1;
vector<MPI_Datatype> m_send_type3;
MPI_Datatype m_send_type_2dfinest[2];
vector<MPI_Datatype> m_send_type_2dx;
vector<MPI_Datatype> m_send_type_2dy;
vector<MPI_Datatype> m_send_type_2dx3p;
vector<MPI_Datatype> m_send_type_2dy3p;
vector<MPI_Datatype> m_send_type_2dx1p;
vector<MPI_Datatype> m_send_type_2dy1p;
bool m_topography_exists;

// UTC time corresponding to simulation time 0.
bool m_utc0set, m_utc0isrefevent;
int m_utc0[7];
};

#endif
