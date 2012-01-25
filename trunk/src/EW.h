//-*-c++-*-
#ifndef EW_H
#define EW_H

#include <mpi.h>

#include <string>
#include <vector>
#include <fstream>
#include <list>

#include "Sarray.h"
#include "SAC.h"
#include "Source.h"

#include "Image.h"
// #include "Image3D.h"

#include "boundaryConditionTypes.h"
#include "ForcingTwilight.h"

#include "MaterialData.h"
// #include "EtreeFile.h"

// #include "SuperGrid.h"
#include "MaterialProperty.h"

using namespace std;

class EW 
{
public:
EW(const string& name, vector<Source*> & a_GlobalUniqueSources);
~EW();
bool wasParsingSuccessful();
bool isInitialized();

void set_output_options( bool output_load, bool output_detailed_timing );
void setGMTOutput(string filename, string wppfilename);
void getGMTOutput( vector<Source*> & a_GlobalUniqueSources );
void allocateCartesianSolverArrays(double a_global_zmax);
void setGoalTime(double t);
//double getCurrentTime(){return mTime;}

void setAttenuationParams(int numberOfMechanisms, double velocityOmega, int ppw, double maxfrequency );

void setNumberSteps(int steps);
int getNumberOfSteps() const;
void setupRun( vector<Source*> & a_GlobalUniqueSources );
void solve( vector<Source*> & a_GlobalUniqueSources );
bool parseInputFile( vector<Source*> & a_GlobalUniqueSources );
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
void processSource(char* buffer, vector<Source*> & a_GlobalUniqueSources);

void side_plane( int g, int side, int wind[6], int nGhost );
void setPrintCycle(int cycle) { mPrintInterval = cycle; }
void setVerbosity(int level) { mVerbose = level; };
int  getVerbosity() {return mVerbose; };
int  getRank() {return m_myRank; };
void setDebugIO(bool onoff) { mDebugIO = onoff; }
  
//void setDampingCFL(double d4_cfl) { m_d4_cfl = d4_cfl; }

void printTime(int cycle, double t, bool force=false ) const;
void printPreamble() const;
void switch_on_checkfornan();
void switch_on_error_log();
void set_energylog( string logfile, bool print, bool elog );
void set_inner_loop( int loopnr );
void set_cflnumber( double cfl );
void set_testing_mode(bool a_testing){m_testing = a_testing;}
bool get_testing_mode(){return m_testing;}

void default_bcs( boundaryConditionType bcs[6] );

void set_twilight_forcing( ForcingTwilight* a_forcing );
// perhaps these functions should be in the ForcingTwilight class? 
// but how will they get access to the material properties and grid sizes?
void exactSolTwilight(double a_t, vector<Sarray> & a_U, vector<Sarray*> & a_AlphaVE);
void exactRhsTwilight(double a_t, vector<Sarray> & a_F);
void exactAccTwilight(double a_t, vector<Sarray> & a_Uacc);
void exactForceTwilight(double a_t, vector<Sarray> & a_F);
void exactForce_ttTwilight(double a_t, vector<Sarray> & a_F);

void normOfDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, double &diffL2 );
void normOfDifferenceGhostPoints( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, double &diffInf, double &diffL2 );
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
void cycleSolutionArrays(vector<Sarray> & a_Um, vector<Sarray> & a_U, vector<Sarray> & a_Up, 
			 vector<Sarray*> & a_AlphaVEm, vector<Sarray*> & a_AlphaVE, vector<Sarray*> & a_AlphaVEp);

void bndryInteriorDifference( vector<Sarray> & a_Uex,  vector<Sarray> & a_U, 
			      double lowZ[3], double interiorZ[3], double highZ[3] );
void test_RhoUtt_Lu( vector<Sarray> & a_Uacc, vector<Sarray> & a_Lu, vector<Sarray> & a_F, 
		     double lowZ[3], double interiorZ[3], double highZ[3] );

void setRestartInfo(int fromCycle, int dumpInterval, const string& filePrefix);
void computeDT();
bool inTestSourceMode() { return mTestSource; }
bool inTestLambMode() { return mTestLamb; }
bool proc_zero() const;
int no_of_procs() const;
void create_output_directory();
void initialize_image_files();
void initialize_SAC_files();
void addSAC(SAC s);
//void addPointSource( Source* source );
void update_SACs( int Nsteps );
void update_images( int Nsteps, double time, vector<Sarray> & a_U );
void print_execution_times( double times[7] );
void print_execution_time( double t1, double t2, string msg );
void finalizeIO();
string bc_name( const boundaryConditionType bc ) const;
int mkdirs(const string& path);
void setOutputPath(const string& path);
const string& getOutputPath() { return mPath; };
const string& getName() { return mName; };
void set_global_bcs(boundaryConditionType bct[6]); // assigns the global boundary conditions

void add_mtrl_block( MaterialData* md );

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

void tune_supergrid_damping(double thickness);
void tune_supergrid_thickness(double thickness);
bool usingSupergrid(){return m_use_supergrid;};
void setup_supergrid( boundaryConditionType a_bc[6] );
void supergrid_taper_material();
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

//void setEtreeFile(EtreeFile* efile); 
void extractTopographyFromEfile(string a_topoFileName, string a_topoExtFileName, string a_QueryType,
                                double a_EFileResolution);
void smoothTopography(int maxIter);

void buildGaussianHillTopography(double amp, double Lx, double Ly, double x0, double y0);

void extractSurfaceFromGridFile(string a_surfaceFileName);

void computeCartesianCoord(double &x, double &y, double &z, const GeographicCoord& g);
void computeGeographicCoord(double x, double y, double z, double & latitude, double & longitude);
   //GeographicCoord& getGeographicCoord();


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

void set_prefilter( bool enable_prefilter, double fc, bool limit_source_freq, double max_freq );

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

void compute_energy( double dt, bool write_file );

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
double getLatOrigin(){ return mGeoCoord.getLatitude();};
double getGridAzimuth(){ return mGeoAz;};
double getMetersPerDegree(){ return mMetersPerDegree;};
bool usingParallelFS(){ return m_pfs;};
int getNumberOfWritersPFS(){ return m_nwriters;};
    
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

private:

ForcingTwilight* m_twilight_forcing;
vector<MaterialData*> m_mtrlblocks;
// index convention: [0]: low-x, [1]: high-x, [2]: low-y, [3]: high-y; [4]: low-z, [5]: high-z  
boundaryConditionType mbcGlobalType[6]; // these are the boundary conditions for the global problem
vector<boundaryConditionType*> m_bcType;  // these are the boundary conditions for each grid on the local processor, with bProcessor conditions
double mTstart;
double mDt;
//EtreeFile * mEtreeFile;

MPI_Comm m_cartesian_communicator;

bool mbcsSet;

// 2-D arrays with elevation-values (=-z) as function of horizontal indices
// mTopo holds the raw topography (according to the "elevation" field in the etree)
// mTopoMat holds the highest elevation where the etree returns solid material properties
// mTopoGrid holds the smoothed topography which follows the top surface of the curvilinear grid
Sarray mTopo, mTopoMat, mTopoGrid;

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

// command prefilter
bool m_prefilter_sources, m_limit_source_freq;
double m_fc, m_source_freq_max;

double m_t0Shift;

// parallel io stuff
bool m_pfs;
int m_nwriters;

// supergrid
bool m_use_supergrid, m_sg_thickness_set;
double m_supergrid_thickness, m_supergrid_damping_coefficient;
//SuperGrid m_supergrid_taper_x, m_supergrid_taper_y, m_supergrid_taper_z;

string mPath;

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
string mWPPFileName;
string mGMTFileName;


bool mWriteGMTOutput;
int mPlotFrequency;
int mNumFiles;
int mVerbose;
bool mDebugIO;
bool mHomogeneous;

// SAC file info
vector<SAC> mSACOutputFiles;
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

bool mTestSource;
bool mTestLamb;
bool mTestingEnergy;
int mOrder;
double mCFL;

vector<int*> m_onesided; 
double m_curlcoeff, m_d4coeff, m_d4_cfl; // these should go away

// storage for the 1-D damping coefficients
vector<double*> m_sg_dc_x, m_sg_dc_y, m_sg_dc_z;

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
GeographicCoord mGeoCoord;
double mMetersPerDegree;

// is this object ready for time-stepping?
bool mParsingSuccessful, mIsInitialized;
bool m_testing;

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
};

#endif
