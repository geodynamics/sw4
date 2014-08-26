//-*-c++-*-
//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
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
#include "Image3D.h"

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
#include "GeographicProjection.h"
#include "DataPatches.h"


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

void setNumberSteps(int steps); // remove???
int getNumberOfSteps() const;

void setupRun( vector<Source*> & a_GlobalUniqueSources );

void solve( vector<Source*> & a_GlobalSources, vector<TimeSeries*> & a_GlobalTimeSeries );
void solve_backward( vector<Source*> & a_Sources, vector<TimeSeries*> & a_TimeSeries, double gradient[11], double hessian[121] );
void solve_allpars( vector<Source*> & a_GlobalSources, vector<Sarray>& a_Rho, vector<Sarray>& a_Mu,
		    vector<Sarray>& a_Lambda, vector<TimeSeries*> & a_GlobalTimeSeries,
		    vector<Sarray>& a_U, vector<Sarray>& a_Um, vector<DataPatches*>& Upred_saved_sides,
		    vector<DataPatches*>& Ucorr_saved_sides, bool save_sides );

void solve_backward_allpars( vector<Source*> & a_GlobalSources, vector<Sarray>& a_Rho, vector<Sarray>& a_Mu,
		    vector<Sarray>& a_Lambda, vector<TimeSeries*> & a_GlobalTimeSeries,
		    vector<Sarray>& a_U, vector<Sarray>& a_Um, vector<DataPatches*>& Upred_saved_sides,
			     vector<DataPatches*>& Ucorr_saved_sides, double gradients[11], 
			     vector<Sarray>& gRho, vector<Sarray>& gMu, vector<Sarray>& gLambda );
   //int nmpar, double* gradientm );

bool parseInputFile( vector<Source*> & a_GlobalSources, vector<TimeSeries*> & a_GlobalTimeSeries );
void parsedate( char* datestr, int& year, int& month, int& day, int& hour, int& minute,
		int& second, int& msecond, int& fail );

void extractRecordData(TimeSeries::receiverMode mode, int i0, int j0, int k0, int grid0, 
		       vector<double> &uRec, vector<Sarray> &Um2, vector<Sarray> &U);

// some (all?) of these functions are called from parseInputFile() and should be made private
void badOption(string name, char* option) const;
bool startswith(const char begin[], char *line);
void processGrid(char* buffer);
void deprecatedOption(const string& command, 
		      const string& oldone, 
		      const string& newone);
void processTime(char* buffer);
void processTwilight(char* buffer);
void processFileIO(char* buffer);
void processImage(char* buffer);
void processImage3D(char* buffer);
void deprecatedImageMode(int value, const char* name) const;
void processTestPointSource(char* buffer);
void processTestRayleigh(char* buffer);
void processTestLamb(char* buffer);
void processTestEnergy(char* buffer);
void processSource(char* buffer, vector<Source*> & a_GlobalUniqueSources);
void processRupture(char* buffer, vector<Source*> & a_GlobalUniqueSources);
void processMaterial( char* buffer );
void processMaterialIfile( char* buffer );
void processMaterialBlock( char* buffer, int & blockCount );
void processMaterialPfile(char* buffer);
void processMaterialEtree(char* buffer);
void processMaterialVimaterial(char* buffer);
void processMaterialInvtest(char* buffer);
void processMaterialRfile(char* buffer);
void processReceiver(char* buffer, vector<TimeSeries*> & a_GlobalTimeSeries);
void processObservation(char* buffer, vector<TimeSeries*> & a_GlobalTimeSeries);
void processBoundaryConditions(char *buffer);
void processPrefilter(char* buffer);
void processGMT(char* buffer);
void processDeveloper(char* buffer);
void processGlobalMaterial(char* buffer);
void processTopography(char* buffer);
void processAttenuation(char* buffer);
void processRandomize(char* buffer);

//void getEfileInfo(char* buffer);

void side_plane( int g, int side, int wind[6], int nGhost );
void setPrintCycle(int cycle) { mPrintInterval = cycle; }
void setVerbosity(int level) { mVerbose = level; };
void setQuiet(bool stealth) { mQuiet = stealth; };
bool getQuiet() {return mQuiet; };
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
void update_curvilinear_cartesian_interface( vector<Sarray>& a_U );

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
void simpleAttenuation( vector<Sarray> & a_Up );
void enforceBC( vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		double t, vector<double **> & a_BCForcing );
void enforceBCfreeAtt( vector<Sarray>& a_Up, vector<Sarray>& a_U, vector<Sarray>& a_Um, 
			   vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
			   vector<Sarray*>& a_AlphaVEp, vector<Sarray*>& a_AlphaVEm,
		       vector<double **>& a_BCForcing, double bop[5], double a_t );
void addAttToFreeBcForcing( vector<Sarray*>& AlphaVEp, vector<double**>& BCForcing, double bop[5] );

void cartesian_bc_forcing( double t, vector<double **> & a_BCForcing, vector<Source*>& a_Source );
void evalRHS(vector<Sarray> & a_U, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda, vector<Sarray> & a_Lu,
	     vector<Sarray*>& a_Alpha );

void evalPredictor(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		   vector<Sarray>& a_Rho, vector<Sarray> & a_Lu, vector<Sarray> & a_F );
void evalDpDmInTime(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
		    vector<Sarray> & a_Uacc );
void evalCorrector(vector<Sarray> & a_Up, vector<Sarray>& a_Rho, vector<Sarray> & a_Lu, vector<Sarray> & a_F );
void updateMemoryVariables( vector<Sarray*>& a_AlphaVEp,
			    vector<Sarray*>& a_AlphaVEm,
			    vector<Sarray>& a_Up, vector<Sarray>& a_U, vector<Sarray>& a_Um, double a_t );
void updateMemoryVariablesBndry( vector<Sarray*>& a_AlphaVEp,
			    vector<Sarray*>& a_AlphaVEm,
			    vector<Sarray>& a_Up, vector<Sarray>& a_U, vector<Sarray>& a_Um );
void evalDpDmInTimeAtt( vector<Sarray*>& a_AlphaVEp, vector<Sarray*>& a_AlphaVE,
                            vector<Sarray*>& a_AlphaVEm );

void addSuperGridDamping(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um, vector<Sarray>& a_Rho );

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
void update_images( int Nsteps, double time, vector<Sarray> & a_Up, vector<Sarray>& a_U, vector<Sarray>& a_Um,
		    vector<Sarray>& a_Rho, vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
		    vector<Source*> & a_sources, int dminus );

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
//inline bool inside_material_surfaces( double lat, double lon )
//    {
//      return (lat <= m_materialLatMax && lat >= m_materialLatMin && 
//	      lon <= m_materialLonMax && lon >= m_materialLonMin);
//    }

void addMaterialProperty(MaterialProperty* mat){m_materials.push_back(mat);}

   //void getMaterialID(double lat, double lon, double depth, int &materialID);
   //bool knownMaterial(int materialID);
   //double lookup_Rho(int materialID, double depth);
   //double lookup_Vs(int materialID, double depth);
   //double lookup_Vp(int materialID, double depth);

//// attenuation model
//double lookup_Qp(int materialID, double depth);
//double lookup_Qs(int materialID, double depth);

// super-grid functions
void processSupergrid(char *buffer);
void set_sg_damping(double coeff);
void set_sg_thickness(int gp_thickness);
//void set_sg_transition(int gp_trans);
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
void check_min_max_int( vector<Sarray>& a_U );

void setupMPICommunications();
void setup2D_MPICommunications();
void communicate_array( Sarray& u, int grid );
void communicate_arrays( vector<Sarray>& u );
void communicate_array_2dfinest( Sarray& u );
void communicate_array_2d( Sarray& u, int g, int k );
void communicate_array_2d_asym( Sarray& u, int g, int k );
void communicate_array_2d_ext( Sarray& u );

void set_materials();
void setup_attenuation_relaxation(double minvsoh );
void setup_viscoelastic();
void setup_viscoelastic_tw();

void extrapolateInZ(int g, Sarray& field, bool lowk, bool highk );
void extrapolateInXY( vector<Sarray>& field );
void extrapolateTopo(Sarray& field);
void checkTopo(Sarray& field);

void addImage(Image* i);
void addImage3D(Image3D* i);
void setIO_timing(bool iotiming);
void setParallel_IO(bool pfs, int nwriters);

void extractTopographyFromGridFile(string a_topoFileName);
void extractTopographyFromImageFile(string a_topoFileName);
void extractTopographyFromCartesianFile(string a_topoFileName);

void setEtreeFile(EtreeFile* efile); 
void extractTopographyFromEfile(string a_topoFileName, string a_topoExtFileName, string a_QueryType,
                                double a_EFileResolution);
void extractTopographyFromRfile( std::string a_topoFileName );

void smoothTopography(int maxIter);

void buildGaussianHillTopography(double amp, double Lx, double Ly, double x0, double y0);

void extractSurfaceFromGridFile(string a_surfaceFileName);
void extractSurfaceFromCartesianFile(string a_surfaceFileName);

void computeCartesianCoord(double &x, double &y, double lon, double lat);
void computeGeographicCoord(double x, double y, double & longitude, double & latitude);

void initializeSystemTime();
void compute_epicenter( vector<Source*> & a_GlobalUniqueSources );
void set_epicenter(double epiLat, double epiLon, double epiDepth, double earliestTime); 
void get_epicenter(double &epiLat, double &epiLon, double &epiDepth, double &earliestTime); 
   
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
bool point_in_proc_ext(int a_i, int a_j, int a_g);      // both interior and parallel ghost points+extra ghost points
  
void initializePaddingCells();

void convert_material_to_mulambda();

void check_materials(); // verify that the density is positive on the grid

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

bool is_onesided( int g, int side ) const;

void interpolate_between_grids( vector<Sarray>& u, vector<Sarray>& um, double t, 
  			        vector<Sarray*> &AlphaVE );

bool interpolate_topography( double q, double r, double & Z0, bool smoothed);

void copy_topo_to_topogridext();

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

void get_gridgen_info( int& order, double& zetaBreak ) const;

//  void update_maxes_hVelMax();
//  void update_maxes_vVelMax();

//   void computeDivergence()   ;
//   void computeCurlMagnitude();

//   void setTestPointSourceMode(){ mTestSource = true; };

// functions from the old FileInput class
void cleanUpRefinementLevels();

enum InputMode { UNDEFINED, Efile, GaussianHill, GridFile, CartesianGrid, TopoImage, Rfile};

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
double getTimeStep() const {return mDt;};

 // test point source
void get_exact_point_source( double* u, double t, int g, Source& source, int* wind=0 );
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
void compute_guess( bool& guesspos, bool& guesst0fr, bool& guessmom, bool& guessshifts, bool& output_seismograms );
void get_cgparameters( int& maxit, int& maxrestart, double& tolerance, bool& fletcherreeves,
		       int& stepselection, bool& do_linesearch, int& varcase, bool& testing );
void parameters_to_material( int nmpar, double* xm, vector<Sarray>& rho,
			     vector<Sarray>& mu, vector<Sarray>& lambda );
void material_to_parameters( int nmpar, double* xm, vector<Sarray>& rho,
			     vector<Sarray>& mu, vector<Sarray>& lambda );
void get_material_parameter( int nmpar, double* xm );
void get_scale_factors( int nmpar, double* xm );

#ifdef ENABLE_OPT
void material_correction( int nmpar, double* xm );

void project_material( vector<Sarray>& a_rho, vector<Sarray>& a_mu,
		       vector<Sarray>& a_lambda, int& info );

void check_material( vector<Sarray>& a_rho, vector<Sarray>& a_mu,
		     vector<Sarray>& a_lambda, int& ok );
#endif

void get_nr_of_material_parameters( int& nmvar );
void add_to_grad( vector<Sarray>& K, vector<Sarray>& Kacc, vector<Sarray>& Um, 
		  vector<Sarray>& U, vector<Sarray>& Up, vector<Sarray>& Uacc,
		  vector<Sarray>& gRho, vector<Sarray>& gMu, vector<Sarray>& gLambda );

void get_optmethod( int& method, int& bfgs_m );
void get_utc( int utc[7] ) const;

void perturb_mtrl();
void perturb_mtrl( int peri, int perj, int perk, double h, int grid, int var );

void perturb_velocities( vector<Sarray>& a_vs, vector<Sarray>& a_vp );

void metric_derivatives_test();

void material_ic( vector<Sarray>& a_mtrl );

void gettopowgh( double ai, double wgh[8] ) const;

void smooth_grid( int maxIter );

void enforceDirichlet5( vector<Sarray> & a_U );

bool check_for_nan( vector<Sarray>& a_U, int verbose, string name );

void define_parallel_io( vector<Parallel_IO*>& parallel_io );

void read_volimage( std::string &path, std::string &fname, vector<Sarray>& data );

void interpolate( int nx, int ny, int nz, double xmin, double ymin, double zmin, double hx,
		  double hy, double hz, Sarray& rho, Sarray& mu, Sarray& lambda,
		  int grid, Sarray& rhogrid, Sarray& mugrid, Sarray& lambdagrid );

void interpolate_to_coarse( int nx, int ny, int nz, double xmin, double ymin,
			    double zmin, double hx, double hy, double hz,
			    Sarray& rho, Sarray& mu, Sarray& lambda, vector<Sarray>& rhogrid, 
			    vector<Sarray>& mugrid, vector<Sarray>& lambdagrid );

void interpolation_gradient( int nx, int ny, int nz, double xmin, double ymin, double zmin, double hx,
			     double hy, double hz, Sarray& gradrho, Sarray& gradmu, Sarray& gradlambda,
			     int grid, Sarray& gradrhogrid, Sarray& gradmugrid, Sarray& gradlambdagrid );

//
// VARIABLES BEYOND THIS POINT
//
const double NO_TOPO;

// ------------------------------------------
// Grid 
// ------------------------------------------

int mNumberOfGrids, mNumberOfCartesianGrids;

// grid sizes are needed by the Source and Image classes, so should be kept public
vector<double> mGridSize;

// part of global array on each processor, including ghost points = all points
vector<int> m_iStart, m_iEnd, m_jStart, m_jEnd, m_kStart, m_kEnd; 

   // Active subcube is the part of the domain where the material is
   // variable in material inversion.
   vector<int> m_iStartAct, m_iEndAct, m_jStartAct, m_jEndAct, m_kStartAct, m_kEndAct; 
   vector<int> m_iStartActGlobal, m_iEndActGlobal, m_jStartActGlobal, m_jEndActGlobal;
   vector<int>  m_kStartActGlobal, m_kEndActGlobal; 

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
// and the metric derivatives as well as the jacobian
Sarray mMetric;

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
// topoMat holds the highest elevation where the etree returns solid material properties (now local to EtreeFile::readEFile() )
// mTopoGridExt holds the smoothed topography which follows the top surface of the curvilinear grid
   Sarray mTopo, mTopoGridExt;

// material description used with material surfaces and the ifile command
vector<MaterialProperty*> m_materials;
MPI_Comm m_cartesian_communicator;

ofstream msgStream;

// vectors of Sarrays hold material properties on all grids. 
vector<Sarray> mMu;
vector<Sarray> mLambda;
vector<Sarray> mRho;

private:
void preprocessSources( vector<Source*> & a_GlobalSources );

// epicenter
double m_epi_lat, m_epi_lon, m_epi_depth, m_epi_t0;

   //PJ *m_projection;
   //double m_xoffset, m_yoffset;
GeographicProjection* m_geoproj;

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

int m_proc_array[2];

bool mbcsSet;


// for some simple topographies (e.g. Gaussian hill) there is an analytical expression for the top elevation
bool m_analytical_topo, m_use_analytical_metric;
double m_GaussianAmp, m_GaussianLx, m_GaussianLy, m_GaussianXc, m_GaussianYc;

// interface surfaces in the material model
//int m_number_material_surfaces, m_Nlon, m_Nlat;
//double m_materialLonMax, m_materialLonMin, m_materialLatMax, m_materialLatMin;
//Sarray m_materialDepth;
//double *m_materialLon, *m_materialLat;

// global material thresholds
bool m_useVelocityThresholds;
double m_vpMin, m_vsMin;


// order of polynomial mapping in algebraic grid genenerator
int m_grid_interpolation_order;
double m_zetaBreak;

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
int m_sg_gp_thickness; //, m_sg_gp_transition;
double m_supergrid_damping_coefficient;
SuperGrid m_supergrid_taper_x, m_supergrid_taper_y, m_supergrid_taper_z;

string mPath, mObsPath, mTempPath;

// number of boundary points on each side
vector<int *> m_NumberOfBCPoints;

// ghost point index window for each side of the boundary on each grid
vector<int *> m_BndryWindow;


// attenuation variables (only allocated if attenuation is enabled)
bool m_use_attenuation, m_att_use_max_frequency;
int m_number_mechanisms;
double m_velo_omega, m_min_omega, m_max_omega, m_att_max_frequency, m_att_ppw;
double m_qmultiplier;

vector<Sarray> mQp, mQs;
vector<Sarray*> mMuVE, mLambdaVE;
// relaxation frequencies
vector<double> mOmegaVE;

// Randomization of the material
bool m_randomize;
int m_random_seed[3];
double m_random_dist, m_random_distz, m_random_amp, m_random_amp_grad;

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
string m_topoFileName, m_topoExtFileName, m_QueryType;
bool mTopoImageFound;
double m_topo_zmax;
int m_maxIter;
double m_EFileResolution;

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
bool mQuiet;

bool mDebugIO;
bool mHomogeneous;

// SAC file info (is this used anymore?)
bool mCompareSACFiles;
float mSACFileErrorTolerance;

// Image file info
vector<Image*> mImageFiles;
vector<Image3D*> mImage3DFiles;
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
   // mCFL actual cfl. Used to determine time step in forward solver.
   // mCFLmax, maximum possible cfl. Used for limiting 
   //          wave speeds during material inversion
double mCFL, mCFLmax;

// info on SBP boundary operators, or not.
vector<int*> m_onesided; 
double m_curlcoeff, m_d4coeff, m_d4_cfl; // these should go away

// storage for the 1-D damping coefficients
vector<double*> m_sg_dc_x, m_sg_dc_y, m_sg_dc_z;
vector<double*> m_sg_str_x, m_sg_str_y, m_sg_str_z;
vector<double*> m_sg_corner_x, m_sg_corner_y, m_sg_corner_z;

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
double mMetersPerDegree, mMetersPerLongitude;
bool mConstMetersPerLongitude;

// is this object ready for time-stepping?
bool mParsingSuccessful, mIsInitialized, mSourcesOK;
bool m_testing;

   // Parameters related to the inverse problem   
bool m_inverse_problem; // Will we solve the inverse problem?
bool m_iniguess_pos, m_iniguess_t0fr, m_iniguess_mom, m_iniguess_shifts;// Estimate initial guess ?
bool m_output_initial_seismograms;
bool m_compute_scalefactors;
int m_maxit,m_maxrestart;
double m_tolerance;
double m_scalefactors[11];   
int m_cgstepselection, m_cgvarcase;
bool m_cgfletcherreeves, m_do_linesearch;
bool m_opt_testing;
int m_opt_method, m_lbfgs_m;
   // perturbations for testing
double m_perturb;
int m_iperturb, m_jperturb, m_kperturb, m_pervar;

// Number of grid points per wave length, P = min Vs/(f*h) 
vector<double> mMinVsOverH;

int m_ext_ghost_points;
int m_ghost_points;
int m_ppadding;

// coefficients for boundary modified 4th order SBP operators
double m_iop[5], m_iop2[5], m_bop2[24], m_sbop[5], m_acof[384], m_bop[24];
double m_bope[48], m_ghcof[6], m_hnorm[4];

int m_neighbor[4];
vector<MPI_Datatype> m_send_type1;
vector<MPI_Datatype> m_send_type3;
vector<MPI_Datatype> m_send_type4; // metric
MPI_Datatype m_send_type_2dfinest[2];
MPI_Datatype m_send_type_2dfinest_ext[2];
vector<MPI_Datatype> m_send_type_2dx;
vector<MPI_Datatype> m_send_type_2dy;
vector<MPI_Datatype> m_send_type_2dx3p;
vector<MPI_Datatype> m_send_type_2dy3p;
vector<MPI_Datatype> m_send_type_2dx1p;
vector<MPI_Datatype> m_send_type_2dy1p;
bool m_topography_exists;

// UTC time corresponding to simulation time 0.
//bool m_utc0set, m_utc0isrefevent;
int m_utc0[7];

// Error handling facility
//ErrorChecking* m_error_checking;
};

#endif
