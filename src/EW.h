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
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
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

#include <fstream>
#include <list>
#include <string>
#include <vector>

#include "AnisotropicMaterial.h"
#include "CheckPoint.h"
#include "DataPatches.h"
#include "ESSI3D.h"
#include "Filter.h"
#include "ForcingTwilight.h"
#include "GeographicProjection.h"
#include "GridPointSource.h"
#include "Image.h"
#include "Image3D.h"
#include "MaterialData.h"
#include "MaterialProperty.h"
#include "Mspace.h"
#include "RandomizedMaterial.h"
#include "Sarray.h"
#include "SfileOutput.h"
#include "Source.h"
#include "SuperGrid.h"
#include "TestEnergy.h"
#include "TestLamb.h"
#include "TestPointSource.h"
#include "TestRayleighWave.h"
#include "TimeSeries.h"
#include "boundaryConditionTypes.h"
#include "policies.h"
#include "sw4.h"
#ifdef SW4_TRACK_MPI
#include "StatMachine.h"
#endif
#include "CurvilinearInterface2.h"
#include "GridGenerator.h"

using namespace std;

class EW {
 public:
  EW(const string& name, vector<vector<Source*>>& a_GlobalUniqueSources,
     vector<vector<TimeSeries*>>& a_GlobalTimeSeries, bool invproblem = false);
  ~EW();
  bool wasParsingSuccessful();
  bool isInitialized();

  void set_output_options(bool output_load, bool output_detailed_timing);
  void setGMTOutput(string filename, string wppfilename);
  void saveGMTFile(vector<vector<Source*>>& a_GlobalUniqueSources, int event);
  void allocateCartesianSolverArrays(float_sw4 a_global_zmax);
  // void setGoalTime(float_sw4 t);
  void setGoalTime(float_sw4 t, int event = 0);
  // double getCurrentTime(){return mTime;}

  // void setNumberSteps(int steps);  // remove???
  void setNumberSteps(int steps, int event = 0);  // remove???
  // int getNumberOfSteps() const;
  int getNumberOfSteps(int event = 0) const;
  int getNumberOfLocalEvents() const;
  int getNumberOfEvents() const;
  float_sw4 getGlobalZmin() { return m_global_zmin; }
  float_sw4 getGlobalZmax() { return m_global_zmax; }
  int findNumberOfEvents();
  bool event_is_in_proc( int e ) const;
  int global_to_local_event( int e ) const;
  int local_to_global_event( int e ) const;

  void setupRun(vector<vector<Source*>>& a_GlobalUniqueSources);

  void solve( vector<Source*> & a_GlobalSources, vector<TimeSeries*> & a_GlobalTimeSeries,
	    vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda, vector<Sarray>& a_Rho,
	    vector<Sarray>& U, vector<Sarray>& Um,
	    vector<DataPatches*>& Upred_saved_sides,
	    vector<DataPatches*>& Ucorr_saved_sides, bool save_sides, int event, int save_steps,
            int varcase, vector<Sarray>& pseudoHessian );
  void solve_backward(vector<Source*>& a_Sources,
                      vector<TimeSeries*>& a_TimeSeries, float_sw4 gradient[11],
                      float_sw4 hessian[121]);
  void solve_allpars(vector<Source*>& a_GlobalSources, vector<Sarray>& a_Rho,
                     vector<Sarray>& a_Mu, vector<Sarray>& a_Lambda,
                     vector<TimeSeries*>& a_GlobalTimeSeries,
                     vector<Sarray>& a_U, vector<Sarray>& a_Um,
                     vector<DataPatches*>& Upred_saved_sides,
                     vector<DataPatches*>& Ucorr_saved_sides, bool save_sides);

  void solve_backward_allpars(vector<Source*>& a_GlobalSources,
                              vector<Sarray>& a_Rho, vector<Sarray>& a_Mu,
                              vector<Sarray>& a_Lambda,
                              vector<TimeSeries*>& a_GlobalTimeSeries,
                              vector<Sarray>& a_U, vector<Sarray>& a_Um,
                              vector<DataPatches*>& Upred_saved_sides,
                              vector<DataPatches*>& Ucorr_saved_sides,
                              float_sw4 gradients[11], vector<Sarray>& gRho,
                              vector<Sarray>& gMu, vector<Sarray>& gLambda);
  // int nmpar, float_sw4* gradientm );

  bool parseInputFile(vector<vector<Source*>>& a_GlobalSources,
                      vector<vector<TimeSeries*>>& a_GlobalTimeSeries);
  void parsedate(char* datestr, int& year, int& month, int& day, int& hour,
                 int& minute, int& second, int& msecond, int& fail);

  void extractRecordData(TimeSeries::receiverMode mode, int i0, int j0, int k0,
                         int grid0, vector<float_sw4>& uRec,
                         vector<Sarray>& Um2, vector<Sarray>& U);

  // some (all?) of these functions are called from parseInputFile() and should
  // be made private
  void badOption(string name, char* option) const;
  bool startswith(const char begin[], char* line);
  void processGrid(char* buffer);
  void processRefinement(char* buffer);
  void deprecatedOption(const string& command, const string& oldone,
                        const string& newone);
  void processTime(char* buffer);
  void processTwilight(char* buffer);
  void processFileIO(char* buffer);
  // void processImage(char* buffer);
  void processImage(char* buffer, bool usehdf5);
  void processImage3D(char* buffer);
  void processESSI3D(char* buffer);
  void processSfileOutput(char* buffer);
  void deprecatedImageMode(int value, const char* name) const;
  void processTestPointSource(char* buffer);
  void processTestRayleigh(char* buffer);
  void processTestLamb(char* buffer);
  void processTestEnergy(char* buffer);
  bool checkTestEnergyPeriodic(char* buffer);
  // void processSource(char* buffer, vector<Source*>& a_GlobalUniqueSources);
  void processSource(char* buffer,
                     vector<vector<Source*>>& a_GlobalUniqueSources);
  // void processRupture(char* buffer, vector<Source*>& a_GlobalUniqueSources);
  void processRupture(char* buffer,
                      vector<vector<Source*>>& a_GlobalUniqueSources);
  void processRuptureHDF5(char* buffer,
                          vector<vector<Source*>>& a_GlobalUniqueSources);
  void processMaterial(char* buffer);
  void processMaterialIfile(char* buffer);
  void processMaterialBlock(char* buffer, int& blockCount);
  void processMaterialPfile(char* buffer);
  void processMaterialVimaterial(char* buffer);
  void processMaterialInvtest(char* buffer);
  void processMaterialRfile(char* buffer);
  void processMaterialSfile(char* buffer);
  void processMaterialGMG(char* buffer);
  void processAnisotropicMaterialBlock(char* buffer, int& ablockCount);
  // void processReceiver(char* buffer, vector<TimeSeries*>&
  // a_GlobalTimeSeries);
  void processReceiver(char* buffer,
                       vector<vector<TimeSeries*>>& a_GlobalTimeSeries);
  void processReceiverHDF5(char* buffer,
                           vector<vector<TimeSeries*>>& a_GlobalTimeSeries);
  // void processObservation(char* buffer,
  //                       vector<TimeSeries*>& a_GlobalTimeSeries);
  void processObservation(char* buffer,
                          vector<vector<TimeSeries*>>& a_GlobalTimeSeries);
  void processObservationHDF5(char* buffer,
                              vector<vector<TimeSeries*>>& a_GlobalTimeSeries);
  void processBoundaryConditions(char* buffer);
  void processPrefilter(char* buffer);
  void processGMT(char* buffer);
  void processDeveloper(char* buffer);
  void processGlobalMaterial(char* buffer);
  void processTopography(char* buffer);
  void processAttenuation(char* buffer);
  void processRandomize(char* buffer);
  void processRandomBlock(char* buffer);
  void processCheckPoint(char* buffer);
  void processGeodynbc(char* buffer);

  void processEvent(char* buffer, int enr);
  // void getEfileInfo(char* buffer);

  void side_plane(int g, int side, int wind[6], int nGhost);
  void setPrintCycle(int cycle) { mPrintInterval = cycle; }
  void setVerbosity(int level) { mVerbose = level; };
  void setQuiet(bool stealth) { mQuiet = stealth; };
  bool getQuiet() { return mQuiet; };
  int getVerbosity() { return mVerbose; };
  int getRank() { return m_myRank; };
  void setDebugIO(bool onoff) { mDebugIO = onoff; }

  // void setDampingCFL(float_sw4 d4_cfl) { m_d4_cfl = d4_cfl; }

  void printTime(int cycle, float_sw4 t, bool force = false) const;
  void printPreamble(vector<Source*>& a_Sources, int event) const;
  void switch_on_checkfornan();
  void switch_on_error_log();
  void set_energylog(string logfile, bool print, bool elog);
  void set_inner_loop(int loopnr);
  void set_cflnumber(float_sw4 cfl);
  void set_testing_mode(bool a_testing) { m_testing = a_testing; }
  bool get_testing_mode() { return m_testing; }

  void default_bcs();
  void update_curvilinear_cartesian_interface(vector<Sarray>& a_U);
  void update_curvilinear_cartesian_interface_org(vector<Sarray>& a_U);

  void set_twilight_forcing(ForcingTwilight* a_forcing);
  // perhaps these functions should be in the ForcingTwilight class?
  // but how will they get access to the material properties and grid sizes?
  void initialData(float_sw4 a_t, vector<Sarray>& a_U,
                   vector<Sarray*>& a_AlphaVE);
  bool exactSol(float_sw4 a_t, vector<Sarray>& a_U, vector<Sarray*>& a_AlphaVE,
                vector<Source*>& source);
  void exactRhsTwilight(float_sw4 a_t, vector<Sarray>& a_F);
  void exactAccTwilight(float_sw4 a_t, vector<Sarray>& a_Uacc);
  void Force(float_sw4 a_t, vector<Sarray>& a_F,
             vector<GridPointSource*>& point_sources,
             vector<int>& identsources);
  void Force_tt(float_sw4 a_t, vector<Sarray>& a_F,
                vector<GridPointSource*>& point_sources,
                vector<int>& identsources);
  // void ForceX(float_sw4 a_t, vector<Sarray> & a_F, vector<GridPointSource*>
  // point_sources, vector<int> identsources ); void ForceX_tt(float_sw4 a_t,
  // vector<Sarray> & a_F, vector<GridPointSource*> point_sources, vector<int>
  // identsources );
  void sort_grid_point_sources(vector<GridPointSource*>& point_sources,
                               vector<int>& identsources);

  void normOfDifference(vector<Sarray>& a_Uex, vector<Sarray>& a_U,
                        float_sw4& diffInf, float_sw4& diffL2, float_sw4& xInf,
                        vector<Source*>& a_globalSources);
  void normOfDifferenceGhostPoints(vector<Sarray>& a_Uex, vector<Sarray>& a_U,
                                   float_sw4& diffInf, float_sw4& diffL2);
  void normOfSurfaceDifference(vector<Sarray>& a_Uex, vector<Sarray>& a_U,
                               float_sw4& diffInf, float_sw4& diffL2,
                               float_sw4& solInf, float_sw4& solL2,
                               vector<Source*>& a_globalSources);

  void test_sources(vector<GridPointSource*>& a_point_sources,
                    vector<Source*>& a_global_unique_sources, vector<Sarray>& F,
                    vector<int>& identsources);
  void testSourceDiscretization(int kx[3], int ky[3], int kz[3],
                                float_sw4 moments[3],
                                vector<GridPointSource*>& point_sources,
                                vector<Sarray>& F, vector<int>& identsources);

  void setupSBPCoeff();

  // time stepping routines
  void simpleAttenuation(vector<Sarray>& a_Up);
  void enforceBC(vector<Sarray>& a_U, vector<Sarray>& a_Rho,
		 vector<Sarray>& a_Mu,
                 vector<Sarray>& a_Lambda, vector<Sarray*>& a_AlphaVE,
                 float_sw4 t, vector<float_sw4**>& a_BCForcing);

  void enforceBCfreeAtt(vector<Sarray>& a_Up, vector<Sarray>& a_U,
                        vector<Sarray>& a_Um, vector<Sarray>& a_Mu,
                        vector<Sarray>& a_Lambda, vector<Sarray*>& a_AlphaVEp,
                        vector<Sarray*>& a_AlphaVEm,
                        vector<float_sw4**>& a_BCForcing, float_sw4 bop[5],
                        float_sw4 a_t);

  void enforceBCfreeAtt2(vector<Sarray>& a_Up, vector<Sarray>& a_Mu,
                         vector<Sarray>& a_Lambda, vector<Sarray*>& a_AlphaVEp,
                         vector<double**>& a_BCForcing);

  void enforceBCanisotropic(vector<Sarray>& a_U, vector<Sarray>& a_C,
                            float_sw4 t, vector<float_sw4**>& a_BCForcing);

  void addAttToFreeBcForcing(vector<Sarray*>& AlphaVEp,
                             vector<float_sw4**>& BCForcing, float_sw4 bop[5]);

  void cartesian_bc_forcing(float_sw4 t, vector<float_sw4**>& a_BCForcing,
                            vector<Source*>& a_Source);

  void cartesian_bc_forcing_olde(float_sw4 t, vector<float_sw4**>& a_BCForcing,
                                 vector<Source*>& a_Source);

  void evalRHS(vector<Sarray>& a_U, vector<Sarray>& a_Mu,
               vector<Sarray>& a_Lambda, vector<Sarray>& a_Lu,
               vector<Sarray*>& a_Alpha,
               std::ostream* norm_trace_file = nullptr);

  void evalRHSanisotropic(vector<Sarray>& a_U, vector<Sarray>& a_C,
                          vector<Sarray>& a_Uacc);

  void evalPredictor(vector<Sarray>& a_Up, vector<Sarray>& a_U,
                     vector<Sarray>& a_Um, vector<Sarray>& a_Rho,
                     vector<Sarray>& a_Lu, vector<Sarray>& a_F);

  void evalDpDmInTime(vector<Sarray>& a_Up, vector<Sarray>& a_U,
                      vector<Sarray>& a_Um, vector<Sarray>& a_Uacc);

  void evalCorrector(vector<Sarray>& a_Up, vector<Sarray>& a_Rho,
                     vector<Sarray>& a_Lu, vector<Sarray>& a_F);

  void updateMemVarPred(vector<Sarray*>& a_AlphaVEp,
                        vector<Sarray*>& a_AlphaVEm, vector<Sarray>& a_U,
                        double a_t);

  void updateMemVarCorr(vector<Sarray*>& a_AlphaVEp,
                        vector<Sarray*>& a_AlphaVEm, vector<Sarray>& a_Up,
                        vector<Sarray>& a_U, vector<Sarray>& a_Um, double a_t);

  void updateMemVarCorrNearInterface(Sarray& a_AlphaVEp, Sarray& a_AlphaVEm,
                                     Sarray& a_Up, Sarray& a_U, Sarray& a_Um,
                                     double a_t, int a_mech, int a_grid);

  // void updateMemoryVariables( vector<Sarray*>& a_AlphaVEp,
  // 			    vector<Sarray*>& a_AlphaVEm,
  // 			    vector<Sarray>& a_Up, vector<Sarray>& a_U,
  // vector<Sarray>& a_Um, double a_t ); void updateMemoryVariablesBndry(
  // vector<Sarray*>&
  // a_AlphaVEp, 			    vector<Sarray*>& a_AlphaVEm,
  // vector<Sarray>& a_Up, vector<Sarray>& a_U, vector<Sarray>& a_Um );
  void evalDpDmInTimeAtt(vector<Sarray*>& a_AlphaVEp,
                         vector<Sarray*>& a_AlphaVE,
                         vector<Sarray*>& a_AlphaVEm);

  void addSuperGridDamping(vector<Sarray>& a_Up, vector<Sarray>& a_U,
                           vector<Sarray>& a_Um, vector<Sarray>& a_Rho);

  void cycleSolutionArrays(vector<Sarray>& a_Um, vector<Sarray>& a_U,
                           vector<Sarray>& a_Up, vector<Sarray*>& a_AlphaVEm,
                           vector<Sarray*>& a_AlphaVE,
                           vector<Sarray*>& a_AlphaVEp);
  void cycleSolutionArrays(vector<Sarray>& a_Um, vector<Sarray>& a_U,
                           vector<Sarray>& a_Up);

  void bndryInteriorDifference(vector<Sarray>& a_Uex, vector<Sarray>& a_U,
                               float_sw4 lowZ[3], float_sw4 interiorZ[3],
                               float_sw4 highZ[3]);
  void test_RhoUtt_Lu(vector<Sarray>& a_Uacc, vector<Sarray>& a_Lu,
                      vector<Sarray>& a_F, float_sw4 lowZ[3],
                      float_sw4 interiorZ[3], float_sw4 highZ[3]);

  void setRestartInfo(int fromCycle, int dumpInterval,
                      const string& filePrefix);
  void computeDT();
  void computeDTanisotropic();
  // bool inTestSourceMode() { return mTestSource; }
  // bool inTestLambMode() { return mTestLamb; }
  bool proc_zero() const;
  bool proc_zero_evzero() const;
  int no_of_procs() const;
  void create_output_directory();
  void create_directory(const string& path);
  void initialize_image_files();
  void update_images(int Nsteps, float_sw4 time, vector<Sarray>& a_Up,
                     vector<Sarray>& a_U, vector<Sarray>& a_Um,
                     vector<Sarray>& a_Rho, vector<Sarray>& a_Mu,
                     vector<Sarray>& a_Lambda, vector<Source*>& a_sources,
                     int dminus, int event = 0);

  void initialize_SAC_files();   // going away
  void update_SACs(int Nsteps);  // going away

  void print_execution_times(float_sw4 times[7]);
  void print_execution_time(float_sw4 t1, float_sw4 t2, string msg);
  void finalizeIO();
  string bc_name(const boundaryConditionType bc) const;
  int mkdirs(const string& path);
  void setOutputPath(const string& path);
  // const string& getOutputPath() {
  //   return mPath;
  // };  // Consider getPath instead! This function has caused grief in the past
  const string& getOutputPath(int event = 0) {
    return mPath[event];
  };  // Consider getPath instead! This function has caused grief in the past
  // const string& getObservationPath() { return mObsPath; };
  const string& getObservationPath(int event) { return mObsPath[event]; };
  const string& getName() { return mName; };
  void set_global_bcs(
      boundaryConditionType bct[6]);  // assigns the global boundary conditions
  boundaryConditionType getLocalBcType(int g, int side) {
    return m_bcType[g][side];
  };

  void add_mtrl_block(MaterialData* md) { m_mtrlblocks.push_back(md); };

  void set_threshold_velocities(float_sw4 vpmin, float_sw4 vsmin);

  // material properties by id
  // inline bool inside_material_surfaces( float_sw4 lat, float_sw4 lon )
  //    {
  //      return (lat <= m_materialLatMax && lat >= m_materialLatMin &&
  //	      lon <= m_materialLonMax && lon >= m_materialLonMin);
  //    }

  void addMaterialProperty(MaterialProperty* mat) {
    m_materials.push_back(mat);
  }

  // void getMaterialID(float_sw4 lat, float_sw4 lon, float_sw4 depth, int
  // &materialID); bool knownMaterial(int materialID); float_sw4 lookup_Rho(int
  // materialID, float_sw4 depth); float_sw4 lookup_Vs(int materialID, float_sw4
  // depth); float_sw4 lookup_Vp(int materialID, float_sw4 depth);

  //// attenuation model
  // float_sw4 lookup_Qp(int materialID, float_sw4 depth);
  // float_sw4 lookup_Qs(int materialID, float_sw4 depth);

  // super-grid functions
  void processSupergrid(char* buffer);
  void set_sg_damping(float_sw4 coeff);
  void set_sg_thickness(int gp_thickness);
  void set_sg_width(float_sw4 sg_width);
  // void set_sg_transition(int gp_trans);
  bool usingSupergrid() { return m_use_supergrid; };
  void setup_supergrid();
  // void supergrid_taper_material();
  void assign_supergrid_damping_arrays();

  // MR coefficients
  void setup_MR_coefficients( vector<Sarray>& Rho, vector<Sarray>& Mu, 
                            vector<Sarray>& Lambda );

  void assign_local_bcs();
  bool timeSteppingSet();
  bool proc_decompose_2d(int ni, int nj, int nproc, int proc_max[2]);
  void decomp1d(int nglobal, int myid, int nproc, int& s, int& e);
  void decomp1d_2( int N, int myid, int nproc, int& s, int& e, int nghost, int npad );
  void coarsen1d(int& n, int& ifirst, int& ilast, int periodic);

  bool node_core_decomp( int ni, int nj, int& Cx, int& Cy, int& Nx, int &Ny );
  void my_node_core_rank( int Cx, int Cy, int Nx, int Ny,
			  int& cx, int& cy, int& nx, int &ny );
  int  my_node_core_id( int ni,int nj, int proc_max[2] );
  
  void allocateCurvilinearArrays();
  void generate_grid();
  void setup_metric();
  void check_min_max_int(vector<Sarray>& a_U);

  void setupMPICommunications();
  void setup2D_MPICommunications();
  void communicate_array(Sarray& u, int grid);
  void communicate_array_host(Sarray& u, int grid);
  void communicate_arrays(vector<Sarray>& u);
  void communicate_host_arrays(vector<Sarray>& u);
  void communicate_array_2dfinest(Sarray& u);
  void communicate_array_2d(Sarray& u, int g, int k);
  void communicate_array_2d_asym(Sarray& u, int g, int k);
  void communicate_array_2d_ext(Sarray& u);
  void communicate_array_2d_ext_async(Sarray& u);
  void communicate_array_2d_isurf(Sarray& u, int iSurf);

  void set_materials();
  void set_anisotropic_materials();
  void setup_attenuation_relaxation(float_sw4 minvsoh);
  void setup_viscoelastic();
  void setup_viscoelastic_tw();
  void reverse_setup_viscoelastic();
  void* use_twilight_forcing() { return m_twilight_forcing; };

  void extrapolateInZ(int g, Sarray& field, bool lowk, bool highk);
  void extrapolateInXY(vector<Sarray>& field);
  void extrapolateInZvector(int g, Sarray& field, bool lowk, bool highk);
  void extrapolateInXYvector(vector<Sarray>& field);
  void extrapolateTopo(Sarray& field);
  void checkTopo(Sarray& field);

  void addImage(Image* i);
  void addImage3D(Image3D* i);
  void addESSI3D(ESSI3D* i);
  void addSfileOutput(SfileOutput* i);
  void setIO_timing(bool iotiming);
  void setParallel_IO(bool pfs, int nwriters);

  void extractTopographyFromGridFile(string a_topoFileName);
  void extractTopographyFromImageFile(string a_topoFileName);
  void extractTopographyFromCartesianFile(string a_topoFileName);

  void extractTopographyFromRfile(std::string a_topoFileName);
  void extractTopographyFromSfile(std::string a_topoFileName);
  void extractTopographyFromGMG(std::string a_topoFileName);

  void smoothTopography(int maxIter);

  void buildGaussianHillTopography(float_sw4 amp, float_sw4 Lx, float_sw4 Ly,
                                   float_sw4 x0, float_sw4 y0);

  void extractSurfaceFromGridFile(string a_surfaceFileName);
  void extractSurfaceFromCartesianFile(string a_surfaceFileName);

  void computeCartesianCoord(float_sw4& x, float_sw4& y, float_sw4 lon,
                             float_sw4 lat);
  void computeCartesianCoordGMG(double& x, double& y, double lon, double lat,
                                char* crs_to);
  void computeGeographicCoord(float_sw4 x, float_sw4 y, float_sw4& longitude,
                              float_sw4& latitude);

  void initializeSystemTime();
  void compute_epicenter(vector<Source*>& a_GlobalUniqueSources, int event = 0);
  void set_epicenter(float_sw4 epiLat, float_sw4 epiLon, float_sw4 epiDepth,
                     float_sw4 earliestTime, int e = 0);
  void get_epicenter(float_sw4& epiLat, float_sw4& epiLon, float_sw4& epiDepth,
                     float_sw4& earliestTime, int e = 0);

  // void update_all_boundaries(vector<Sarray> &U, vector<Sarray> &UM, float_sw4
  // t, 			   vector<Sarray*> &AlphaVE );

  // void impose_physical_bc(vector<Sarray> &U, vector<Sarray> &UM, float_sw4 t,
  // 			vector<Sarray*> &AlphaVE );

  /* void bc_dirichlet( Sarray& u, int g, float_sw4 t, int side, */
  /*  		   Forcing* forcing, float_sw4 h ); */

  /* void bc_free_surface( Sarray& u, int g, float_sw4 t, int side, */
  /* 		      Forcing* forcing, float_sw4 h, int onesided[6] ); */

  void computeNearestTopoGridPoint(int & iNear, 
                                 int & jNear, 
                                 float_sw4 a_x, 
                                 float_sw4 a_y);

void computeLowTopoGridPoint(int & iLow, 
                             int & jLow, 
                             float_sw4 a_x, 
                             float_sw4 a_y);
  
  void computeNearestGridPoint(int& a_i, int& a_j, int& a_k,
                               int& a_g,  // grid on which indices are located
                               float_sw4 a_x, float_sw4 a_y, float_sw4 a_z);
  int computeNearestGridPoint2(int& a_i, int& a_j, int& a_k,
                               int& a_g,  // grid on which indices are located
                               float_sw4 a_x, float_sw4 a_y, float_sw4 a_z);

  int computeInvGridMap( float_sw4& a_i, float_sw4& a_j, float_sw4& a_k, int& a_g,
                       float_sw4 a_x, float_sw4 a_y, float_sw4 a_z );

  void computeNearestSurfaceGridPoint(int& a_i, int& a_j, float_sw4 a_x,
                                      float_sw4 a_y, float_sw4 a_z);

  void coord_from_index(int i, int j, int k, int g, float_sw4& x, float_sw4& y,
                        float_sw4& z);

  float_sw4 distance(float_sw4 a_x1, float_sw4 a_y1, float_sw4 a_z1,
                     float_sw4 a_x0, float_sw4 a_y0, float_sw4 a_z0);

  void computeNearestLowGridPoint(
      int& a_i, int& a_j, int& a_k,
      int& a_g,  // grid on which indices are located
      float_sw4 a_x, float_sw4 a_y, float_sw4 a_z);

  bool interior_point_in_proc(
      int a_i, int a_j, int a_g);  // only takes interior points into account
  bool point_in_proc(int a_i, int a_j,
                     int a_g);  // both interior and parallel ghost points
  bool point_in_proc_ext(
      int a_i, int a_j,
      int a_g);  // both interior and parallel ghost points+extra ghost points

  void initializePaddingCells();
  void check_dimensions();

  void convert_material_to_mulambda();

  void check_materials();  // verify that the density is positive on the grid

  void computeSolutionError(vector<Sarray>& U, float_sw4 t,
                            vector<Sarray*>& Alpha);

  float_sw4 localMin(vector<Sarray>& a_field);
  float_sw4 localMax(vector<Sarray>& a_field);
  float_sw4 localMinVp();
  float_sw4 localMaxVp();
  float_sw4 localMinVs();
  float_sw4 localMaxVs();
  float_sw4 localMinVpOverVs();
  float_sw4 localMaxVpOverVs();

  bool topographyExists() { return m_topography_exists; };
  bool usingAttenuation() { return m_use_attenuation; };

  bool is_onesided(int g, int side) const;

  void interpolate_between_grids(vector<Sarray>& u, vector<Sarray>& um,
                                 float_sw4 t, vector<Sarray*>& AlphaVE);

  int interpolate_topography(float_sw4 q, float_sw4 r, float_sw4& Z0,
                             bool smoothed);

  void copy_topo_to_topogridext();

  bool getDepth(float_sw4 x, float_sw4 y, float_sw4 z, float_sw4& depth);

  bool curvilinear_grid_mapping(float_sw4 q, float_sw4 r, float_sw4 s,
                                float_sw4& X0, float_sw4& Y0, float_sw4& Z0);

  bool invert_curvilinear_grid_mapping(float_sw4 X0, float_sw4 Y0, float_sw4 Z0,
                                       float_sw4& q, float_sw4& r,
                                       float_sw4& s);

  bool find_curvilinear_derivatives_at_point(float_sw4 q, float_sw4 r,
                                             float_sw4 s, float_sw4 qX[],
                                             float_sw4 rX[], float_sw4 sX[]);

  void save_errors(float_sw4 max_error[3], float_sw4 l2_error[3]);

  void compute_minvsoverh(float_sw4& minvsoh);

  void set_resolution(int ppw);

  void set_prefilter(FilterType passband, int order, int passes, float_sw4 fc1,
                     float_sw4 fc2);

  void set_scenario(const string& scenario);

  void set_conservative_interpolation(bool onoff, float_sw4 ctol, int cmaxit);

  void set_geodyn_data(string filename, int nx, int nz, float_sw4 h,
                       float_sw4 origin[3], float_sw4 dt, int nsteps,
                       int faces);

  void impose_geodyn_ibcdata(vector<Sarray>& u, vector<Sarray>& um, float_sw4 t,
                             vector<float_sw4**>& bforcing);

  void advance_geodyn_time(float_sw4 t);

  void get_geodyn_timelevel(vector<Sarray>& geodyndata);

  void copy_geodyn_timelevel(vector<Sarray>& geodyndata1,
                             vector<Sarray>& geodyndata2);

  // void geodyn_second_ghost_point(vector<Sarray>& geodyndata1,
  //                                vector<Sarray>& geodyndata2,
  //                                vector<Sarray>& rho, vector<Sarray>& mu,
  //                                vector<Sarray>& lambda,
  //                                vector<Sarray>& forcing, double t,
  //                                vector<Sarray>& U, vector<Sarray>& Um,
  //                                int crf);

  void geodyn_second_ghost_point(vector<Sarray>& rho, vector<Sarray>& mu,
                                 vector<Sarray>& lambda,
                                 vector<Sarray>& forcing, float_sw4 t,
                                 vector<Sarray>& U, vector<Sarray>& Um,
                                 int crf);

  void geodyn_second_ghost_point_curvilinear(vector<Sarray>& rho,
                                             vector<Sarray>& mu,
                                             vector<Sarray>& lambda,
                                             vector<Sarray>& forcing,
                                             float_sw4 t, vector<Sarray>& U,
                                             vector<Sarray>& Um, int crf);

  void geodyn_up_from_uacc(vector<Sarray>& Up, vector<Sarray>& Uacc,
                           vector<Sarray>& U, vector<Sarray>& Um, float_sw4 dt);

  void save_geoghost(vector<Sarray>& U);
  void restore_geoghost(vector<Sarray>& U);

  void geodynbcGetSizes(string filename, float_sw4 origin[3],
                        float_sw4& cubelen, float_sw4& zcubelen,
                        float_sw4& hcube, bool& found_latlon, double& lat,
                        double& lon, double& az, int& adjust);
  void geodynFindFile(char* buffer);
  void bcsurf_curvilinear_2nd_order(int side, int i0, int i1, int j0, int j1,
                                    int k0, int g, Sarray& u,
                                    float_sw4* bforcing);
  void integrate_source();

  void compute_energy(float_sw4 dt, bool write_file, vector<Sarray>& Um,
                      vector<Sarray>& U, vector<Sarray>& Up, int step,
                      int event);

  float_sw4 scalarProduct(vector<Sarray>& U, vector<Sarray>& V);
  void get_gridgen_info(int& order, float_sw4& zetaBreak) const;

  //  void update_maxes_hVelMax();
  //  void update_maxes_vVelMax();

  //   void computeDivergence()   ;
  //   void computeCurlMagnitude();

  //   void setTestPointSourceMode(){ mTestSource = true; };

  // functions from the old FileInput class
  void cleanUpRefinementLevels();

  enum InputMode {
    UNDEFINED,
    GaussianHill,
    GridFile,
    CartesianGrid,
    TopoImage,
    Rfile,
    Sfile,
    GMG
  };

  // access functions needed by the Image (and perhaps other) classes
  int getNumberOfCartesianGrids() { return mNumberOfCartesianGrids; };
  int getNumberOfGrids() { return mNumberOfGrids; };
  int getNumberOfGhostPoints() { return m_ghost_points; };
  int getNumberOfParallelPaddingPoints() { return m_ppadding; };
  float_sw4 getLonOrigin() { return mLonOrigin; };
  float_sw4 getLatOrigin() { return mLatOrigin; };
  float_sw4 getGridAzimuth() { return mGeoAz; };
  float_sw4 getMetersPerDegree() { return mMetersPerDegree; };
  bool usingParallelFS() { return m_pfs; };
  int getNumberOfWritersPFS() { return m_nwriters; };
  float_sw4 getTimeStep() const { return mDt; };
  // int getNumberOfTimeSteps() const { return mNumberOfTimeSteps; };
  int getNumberOfTimeSteps(int event = 0) const {
    return mNumberOfTimeSteps[event];
  };
  int getNumberOfMechanisms() const { return m_number_mechanisms; };

  // test point source
  void get_exact_point_source(float_sw4* u, float_sw4 t, int g, Source& source,
                              int* wind = 0);
  RAJA_HOST_DEVICE static float_sw4 VerySmoothBump_x_T_Integral(float_sw4 t,
                                                                float_sw4 R,
                                                                float_sw4 alpha,
                                                                float_sw4 beta);
  RAJA_HOST_DEVICE static float_sw4 C6SmoothBump_x_T_Integral(float_sw4 t,
                                                              float_sw4 R,
                                                              float_sw4 alpha,
                                                              float_sw4 beta);
  RAJA_HOST_DEVICE static float_sw4 SmoothWave_x_T_Integral(float_sw4 t,
                                                            float_sw4 R,
                                                            float_sw4 alpha,
                                                            float_sw4 beta);
  RAJA_HOST_DEVICE static float_sw4 Gaussian_x_T_Integral(
      float_sw4 t, float_sw4 R, float_sw4 f, float_sw4 alpha, float_sw4 beta);
  RAJA_HOST_DEVICE static float_sw4 VSBTP(float_sw4 Lim, float_sw4 t);
  RAJA_HOST_DEVICE static float_sw4 C6SBTP(float_sw4 Lim, float_sw4 t);
  RAJA_HOST_DEVICE static float_sw4 SWTP(float_sw4 Lim, float_sw4 t);
  RAJA_HOST_DEVICE static float_sw4 d_VerySmoothBump_dt(float_sw4 t,
                                                        float_sw4 R,
                                                        float_sw4 c);
  RAJA_HOST_DEVICE static float_sw4 d_C6SmoothBump_dt(float_sw4 t, float_sw4 R,
                                                      float_sw4 c);
  RAJA_HOST_DEVICE static float_sw4 d_SmoothWave_dt(float_sw4 t, float_sw4 R,
                                                    float_sw4 c);
  RAJA_HOST_DEVICE static float_sw4 d_Gaussian_dt(float_sw4 t, float_sw4 R,
                                                  float_sw4 c, float_sw4 f);
  RAJA_HOST_DEVICE static float_sw4 VerySmoothBump(float_sw4 t, float_sw4 R,
                                                   float_sw4 c);
  RAJA_HOST_DEVICE static float_sw4 C6SmoothBump(float_sw4 t, float_sw4 R,
                                                 float_sw4 c);
  RAJA_HOST_DEVICE static float_sw4 SmoothWave(float_sw4 t, float_sw4 R,
                                               float_sw4 c);
  RAJA_HOST_DEVICE static float_sw4 Gaussian(float_sw4 t, float_sw4 R,
                                             float_sw4 c, float_sw4 f);

  // Lamb's problem
  void get_exact_lamb(vector<Sarray>& a_U, float_sw4 a_t, Source& a_source);
  void get_exact_lamb2(vector<Sarray>& a_U, float_sw4 a_t, Source& a_source);
  float_sw4 G4_Integral(float_sw4 T, float_sw4 t, float_sw4 r, float_sw4 beta);
  float_sw4 G3_Integral(float_sw4 iT, float_sw4 it, float_sw4 ir,
                        float_sw4 ibeta);
  float_sw4 G2_Integral(float_sw4 iT, float_sw4 it, float_sw4 ir,
                        float_sw4 ibeta);

  void getGlobalBoundingBox(float_sw4 bbox[6]);

  // string getPath() { return mPath; }
  string getPath(int event = 0) { return mPath[event]; }
  void set_utcref(TimeSeries& ts);
  void print_utc(int event = 0);

  // For inverse problem
  void processCG(char* buffer);
  void processScaleFactors(char* buffer);
  void average_speeds(float_sw4& cp, float_sw4& cs);
  void layered_speeds(vector<float_sw4>& cp, vector<float_sw4>& z);
  void testsourcediff(vector<Source*> GlobalSources, float_sw4 gradient[11],
                      float_sw4 hessian[121]);
  void get_scalefactors(float_sw4 sf[11]);
  bool compute_sf();
  void compute_guess(bool& guesspos, bool& guesst0fr, bool& guessmom,
                     bool& guessshifts, bool& output_seismograms);
  void get_cgparameters(int& maxit, int& maxrestart, float_sw4& tolerance,
                        bool& fletcherreeves, int& stepselection,
                        bool& do_linesearch, int& varcase, bool& testing);
  void parameters_to_material(int nmpar, float_sw4* xm, vector<Sarray>& rho,
                              vector<Sarray>& mu, vector<Sarray>& lambda);
  void material_to_parameters(int nmpar, float_sw4* xm, vector<Sarray>& rho,
                              vector<Sarray>& mu, vector<Sarray>& lambda);
  void get_material_parameter(int nmpar, float_sw4* xm);
  void get_scale_factors(int nmpar, float_sw4* xm);

#ifdef ENABLE_OPT
  void material_correction(int nmpar, float_sw4* xm);

  void project_material(vector<Sarray>& a_rho, vector<Sarray>& a_mu,
                        vector<Sarray>& a_lambda, int& info);

  void check_material(vector<Sarray>& a_rho, vector<Sarray>& a_mu,
                      vector<Sarray>& a_lambda, int& ok);
#endif

  void check_anisotropic_material(vector<Sarray>& rho, vector<Sarray>& c);

  void get_nr_of_material_parameters(int& nmvar);
  void add_to_grad(vector<Sarray>& K, vector<Sarray>& Kacc, vector<Sarray>& Um,
                   vector<Sarray>& U, vector<Sarray>& Up, vector<Sarray>& Uacc,
                   vector<Sarray>& gRho, vector<Sarray>& gMu,
                   vector<Sarray>& gLambda);

  void get_optmethod(int& method, int& bfgs_m);
  void get_utc(int utc[7], int event = 0) const;

  void perturb_mtrl();
  void perturb_mtrl(int peri, int perj, int perk, float_sw4 h, int grid,
                    int var);

  void perturb_velocities(vector<Sarray>& a_vs, vector<Sarray>& a_vp);

  void metric_derivatives_test();

  void material_ic(vector<Sarray>& a_mtrl);

  void gettopowgh(float_sw4 ai, float_sw4 wgh[8]) const;

  void smooth_grid(int maxIter);

  void enforceDirichlet5(vector<Sarray>& a_U);

  bool check_for_nan(vector<Sarray>& a_U, int verbose, string name);

  bool check_for_nan( vector<Sarray*>& a_U, int nmech, int verbose, string name );

  void define_parallel_io(vector<Parallel_IO*>& parallel_io);

  void read_volimage(std::string& path, std::string& fname,
                     vector<Sarray>& data);

  void interpolate(int nx, int ny, int nz, float_sw4 xmin, float_sw4 ymin,
                   float_sw4 zmin, float_sw4 hx, float_sw4 hy, float_sw4 hz,
                   Sarray& rho, Sarray& mu, Sarray& lambda, int grid,
                   Sarray& rhogrid, Sarray& mugrid, Sarray& lambdagrid);

  void interpolate_to_coarse(int nx, int ny, int nz, float_sw4 xmin,
                             float_sw4 ymin, float_sw4 zmin, float_sw4 hx,
                             float_sw4 hy, float_sw4 hz, Sarray& rho,
                             Sarray& mu, Sarray& lambda,
                             vector<Sarray>& rhogrid, vector<Sarray>& mugrid,
                             vector<Sarray>& lambdagrid);

  void interpolation_gradient(int nx, int ny, int nz, float_sw4 xmin,
                              float_sw4 ymin, float_sw4 zmin, float_sw4 hx,
                              float_sw4 hy, float_sw4 hz, Sarray& gradrho,
                              Sarray& gradmu, Sarray& gradlambda, int grid,
                              Sarray& gradrhogrid, Sarray& gradmugrid,
                              Sarray& gradlambdagrid);

  // Functions to impose conditions at grid refinement interface:
  // void enforceIC( std::vector<Sarray> & a_Up, std::vector<Sarray> & a_U,
  // std::vector<Sarray> & a_Um,
  //                 vector<Sarray*>& a_AlphaVEp,
  //      	   double t, bool predictor, std::vector<GridPointSource*>
  //      point_sources );
  // NEW June 14, 2017

  // void enforceIC(std::vector<Sarray>& a_Up, std::vector<Sarray>& a_U,
  //                std::vector<Sarray>& a_Um, float_sw4 t, bool predictor,
  //                std::vector<GridPointSource*>& point_sources);
  // void dirichlet_hom_ic( Sarray& U, int g, int k, bool inner );
  // void dirichlet_LRic( Sarray& U, int g, int kic, float_sw4 t, int adj );
  // void gridref_initial_guess( Sarray& u, int g, bool upper );
  // void compute_preliminary_corrector( Sarray& a_Up, Sarray& a_U, Sarray&
  // a_Um, Sarray& Unext, 				    int g, int kic,
  // float_sw4 t, std::vector<GridPointSource*> point_sources ); void
  // compute_preliminary_predictor( Sarray& a_Up, Sarray& a_U, Sarray& Unext,
  // int
  // g, int kic, float_sw4 t, std::vector<GridPointSource*> point_sources );
  // void
  // compute_icstresses( Sarray& a_Up, Sarray& B, int g, int kic, float_sw4*
  // a_str_x, float_sw4* a_str_y); void consintp( Sarray& Uf, Sarray& Unextf,
  // Sarray& Bf, Sarray& Muf, Sarray& Lambdaf, Sarray& Rhof, float_sw4 hf,
  //	       Sarray& Uc, Sarray& Unextc, Sarray& Bc, Sarray& Muc, Sarray&
  // Lambdac, Sarray& Rhoc, float_sw4 hc, 	       float_sw4 cof, int gc,
  // int gp, int is_periodic[2] ); void check_corrector( Sarray& Uf, Sarray& Uc,
  // Sarray&
  // Unextf, Sarray& Unextc, int kf, int kc );

  // Previous version fortran routines, now in C
  void addsgd4_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                  int klast, float_sw4* a_Up, float_sw4* a_U, float_sw4* a_Um,
                  float_sw4* Rho, float_sw4* sg_dc_x, float_sw4* sg_dc_y,
                  float_sw4* sg_dc_z, float_sw4* sg_str_x, float_sw4* sg_str_y,
                  float_sw4* sg_str_z, float_sw4* sg_corner_x,
                  float_sw4* sg_corner_y, float_sw4* sg_corner_z,
                  float_sw4 damping_coefficient);
  void addsgd6_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                  int klast, float_sw4* a_Up, float_sw4* a_U, float_sw4* a_Um,
                  float_sw4* Rho, float_sw4* sg_dc_x, float_sw4* sg_dc_y,
                  float_sw4* sg_dc_z, float_sw4* sg_str_x, float_sw4* sg_str_y,
                  float_sw4* sg_str_z, float_sw4* sg_corner_x,
                  float_sw4* sg_corner_y, float_sw4* sg_corner_z,
                  float_sw4 damping_coefficient);
  void addsgd4c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, float_sw4* a_Up, float_sw4* a_U, float_sw4* a_Um,
                   float_sw4* Rho, float_sw4* sg_dc_x, float_sw4* sg_dc_y,
                   float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* jac,
                   float_sw4* sg_corner_x, float_sw4* sg_corner_y,
                   float_sw4 damping_coefficient);
  void addsgd6c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, float_sw4* a_Up, float_sw4* a_U, float_sw4* a_Um,
                   float_sw4* Rho, float_sw4* sg_dc_x, float_sw4* sg_dc_y,
                   float_sw4* sg_str_x, float_sw4* sg_str_y, float_sw4* jac,
                   float_sw4* sg_corner_x, float_sw4* sg_corner_y,
                   float_sw4 damping_coefficient);

  void bcfort_ci(int ib, int ie, int jb, int je, int kb, int ke, int wind[36],
                 int nx, int ny, int nz, float_sw4* u, float_sw4 h,
                 boundaryConditionType bccnd[6], float_sw4 sbop[6],
                 float_sw4* mu, float_sw4* la, float_sw4 t, float_sw4* bforce1,
                 float_sw4* bforce2, float_sw4* bforce3, float_sw4* bforce4,
                 float_sw4* bforce5, float_sw4* bforce6, float_sw4 om,
                 float_sw4 ph, float_sw4 cv, int curvilinear);
  void bcfortsg_ci(int ib, int ie, int jb, int je, int kb, int ke, int wind[36],
                   int nx, int ny, int nz, float_sw4* u, float_sw4 h,
                   boundaryConditionType bccnd[6], float_sw4 sbop[6],
                   float_sw4* mu, float_sw4* la, float_sw4 t,
                   float_sw4* bforce1, float_sw4* bforce2, float_sw4* bforce3,
                   float_sw4* bforce4, float_sw4* bforce5, float_sw4* bforce6,
                   float_sw4 om, float_sw4 ph, float_sw4 cv, float_sw4* strx,
                   float_sw4* stry);
  void twdirbdry_ci(int wind[6], float_sw4 h, float_sw4 t, float_sw4 om,
                    float_sw4 cv, float_sw4 ph, float_sw4* bforce,
                    float_sw4 zmin);
  void twdirbdryc_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                     int klast, int wind[6], float_sw4 t, float_sw4 om,
                     float_sw4 cv, float_sw4 ph, float_sw4* bforce,
                     float_sw4* x, float_sw4* y, float_sw4* z);
  void twfrsurfz_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                    int klast, float_sw4 h, int kz, float_sw4 t,
                    float_sw4 omega, float_sw4 c, float_sw4 phase,
                    float_sw4* bforce, float_sw4* mu, float_sw4* lambda,
                    float_sw4 zmin);
  void twfrsurfzatt_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4 h, int kz, float_sw4 t,
                       float_sw4 omega, float_sw4 c, float_sw4 phase,
                       float_sw4* bforce, float_sw4* mua, float_sw4* lambdaa,
                       float_sw4 zmin);
  void twfrsurfzsgstr_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4 h, int kz,
                         float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                         float_sw4 omstrx, float_sw4 omstry, float_sw4* bforce,
                         float_sw4* mu, float_sw4* lambda, float_sw4 zmin);
  void twfrsurfzsgstratt_ci(int ifirst, int ilast, int jfirst, int jlast,
                            int kfirst, int klast, float_sw4 h, int kz,
                            float_sw4 t, float_sw4 omega, float_sw4 c,
                            float_sw4 phase, float_sw4 omstrx, float_sw4 omstry,
                            float_sw4* bforce, float_sw4* mua,
                            float_sw4* lambdaa, float_sw4 zmin);
  void twstensor_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                    int klast, int kz, float_sw4 t, float_sw4 om, float_sw4 c,
                    float_sw4 ph, float_sw4* xx, float_sw4* yy, float_sw4* zz,
                    float_sw4* tau, float_sw4* mu, float_sw4* lambda);
  void twstensorsg_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                      int klast, int kz, float_sw4 t, float_sw4 om, float_sw4 c,
                      float_sw4 ph, float_sw4* xx, float_sw4* yy, float_sw4* zz,
                      float_sw4* tau, float_sw4* mu, float_sw4* lambda,
                      float_sw4 omstrx, float_sw4 omstry);
  void twstensoratt_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, int kz, float_sw4 t, float_sw4 omega,
                       float_sw4 c, float_sw4 phase, float_sw4* xx,
                       float_sw4* yy, float_sw4* zz, float_sw4* tau,
                       float_sw4* mu, float_sw4* lambda);
  void twstensorsgatt_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, int kz, float_sw4 t,
                         float_sw4 omega, float_sw4 c, float_sw4 phase,
                         float_sw4* xx, float_sw4* yy, float_sw4* zz,
                         float_sw4* tau, float_sw4* mu, float_sw4* lambda,
                         float_sw4 omstrx, float_sw4 omstry);
  void bcfortanisg_ci(int ib, int ie, int jb, int je, int kb, int ke,
                      int wind[36], int nx, int ny, int nz, float_sw4* u,
                      float_sw4 h, boundaryConditionType bccnd[6],
                      float_sw4 sbop[6], float_sw4* c, float_sw4* bforce1,
                      float_sw4* bforce2, float_sw4* bforce3,
                      float_sw4* bforce4, float_sw4* bforce5,
                      float_sw4* bforce6, float_sw4* strx, float_sw4* stry);
  void bcfreesurfcurvani_ci(int ifirst, int ilast, int jfirst, int jlast,
                            int kfirst, int klast, int nz, float_sw4* u,
                            float_sw4* c, int side, float_sw4 sbop[6],
                            float_sw4* bforce5, float_sw4* bforce6,
                            float_sw4* strx, float_sw4* stry);
  void GetStencilCoefficients(float_sw4* _acof, float_sw4* _ghcof,
                              float_sw4* _bop, float_sw4* _bope,
                              float_sw4* _sbop);
  void checkanisomtrl_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4* rho, float_sw4* c,
                         float_sw4& rhomin, float_sw4& rhomax,
                         float_sw4& eigmin, float_sw4& eigmax);
  void maxwave(float_sw4 c[21], float_sw4 rho, float_sw4& eigestimate);
  void maxwavecurv(float_sw4 c[21], float_sw4 rho, float_sw4 jac,
                   float_sw4& eigestimate);
  void computedtaniso2_ci(int ifirst, int ilast, int jfirst, int jlast,
                          int kfirst, int klast, float_sw4* rho, float_sw4* c,
                          float_sw4 cfl, float_sw4 dx, float_sw4& a_dtloc);
  void computedtaniso2curv_ci(int ifirst, int ilast, int jfirst, int jlast,
                              int kfirst, int klast, float_sw4* rho,
                              float_sw4* c, float_sw4* jac, float_sw4 cfl,
                              float_sw4& a_dtloc);
  void randomfield3d_ci(int ifirst, int ilast, int jfirst, int jlast,
                        int kfirst, int klast, int nig, int njg, int nkg,
                        int gh, float_sw4* __restrict__ a_w,
                        float_sw4* __restrict__ a_wgh, float_sw4 dist,
                        float_sw4 distz, float_sw4 h, int* randw,
                        float_sw4* __restrict__ a_saverands, int p, int pz);
  void randomfield3dc_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, int nig, int njg, int nkg,
                         int gh, float_sw4* __restrict__ a_w,
                         float_sw4* __restrict__ a_wgh, float_sw4 dist,
                         float_sw4 distz, float_sw4 h,
                         float_sw4* __restrict__ a_z, int* randw,
                         float_sw4* __restrict__ a_saverands, int p, int pz);
  void perturbvelocity_ci(int ifirst, int ilast, int jfirst, int jlast,
                          int kfirst, int klast, float_sw4* __restrict__ a_vs,
                          float_sw4* __restrict__ a_vp,
                          float_sw4* __restrict__ a_per, float_sw4 amp,
                          float_sw4 grad, float_sw4 zmin, float_sw4 h,
                          float_sw4 plimit);
  void perturbvelocityc_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* a_vs,
                           float_sw4* a_vp, float_sw4* a_per, float_sw4 amp,
                           float_sw4 grad, float_sw4* a_z, float_sw4 plimit);

  void gridinfo_ci(int ib, int ie, int jb, int je, int kb, int ke,
                   float_sw4* met, float_sw4* jac, float_sw4& minj,
                   float_sw4& maxj);

  int metric_ci(int ib, int ie, int jb, int je, int kb, int ke, float_sw4* a_x,
                float_sw4* a_y, float_sw4* a_z, float_sw4* a_met,
                float_sw4* a_jac);

  void metricexgh_ci(int ib, int ie, int jb, int je, int kb, int ke, int nz,
                     float_sw4* a_x, float_sw4* a_y, float_sw4* a_z,
                     float_sw4* a_met, float_sw4* a_jac, int order,
                     float_sw4 sb, float_sw4 zmax, float_sw4 amp, float_sw4 xc,
                     float_sw4 yc, float_sw4 xl, float_sw4 yl);

  void freesurfcurvi_ci(int ib, int ie, int jb, int je, int kb, int ke, int nz,
                        int side, float_sw4* a_u, float_sw4* a_mu,
                        float_sw4* a_la, float_sw4* a_met, float_sw4* s,
                        float_sw4* a_forcing);

  void freesurfcurvisg_ci(int ib, int ie, int jb, int je, int kb, int ke,
                          int nz, int side, float_sw4* a_u, float_sw4* a_mu,
                          float_sw4* a_la, float_sw4* a_met, float_sw4* s,
                          float_sw4* a_forcing, float_sw4* a_strx,
                          float_sw4* a_stry);

  void getsurfforcing_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, int k, float_sw4* a_met,
                         float_sw4* a_jac, float_sw4* a_tau,
                         float_sw4* a_forcing);

  void getsurfforcingsg_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, int k, float_sw4* a_met,
                           float_sw4* a_jac, float_sw4* a_tau,
                           float_sw4* a_strx, float_sw4* a_stry,
                           float_sw4* a_forcing);

  void getsurfforcinggh_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, int k, float_sw4 h,
                           float_sw4* a_tau, float_sw4* a_forcing,
                           float_sw4 amp, float_sw4 xc, float_sw4 yc,
                           float_sw4 xl, float_sw4 yl);

  void subsurfforcing_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, int k, float_sw4* a_met,
                         float_sw4* a_jac, float_sw4* a_tau,
                         float_sw4* a_forcing);

  void subsurfforcingsg_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, int k, float_sw4* a_met,
                           float_sw4* a_jac, float_sw4* a_tau,
                           float_sw4* a_strx, float_sw4* a_stry,
                           float_sw4* a_forcing);

  void addbstressc_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                      int klast, int nz, float_sw4* a_u, float_sw4* a_mu,
                      float_sw4* a_la, float_sw4* a_bs, float_sw4* a_met,
                      int side, float_sw4* s, char op, int ghterm, int usesg,
                      float_sw4* a_sgstrx, float_sw4* a_sgstry);

  // void updatememvar_ci( int ifirst, int ilast, int jfirst, int jlast, int
  // kfirst, int klast, 		      float_sw4* __restrict__ a_alp,
  // float_sw4* __restrict__
  // a_alm, float_sw4* __restrict__ a_up, 		      float_sw4*
  // __restrict__ a_u, float_sw4*
  //__restrict__ a_um, 		      float_sw4 omega, float_sw4 dt, int domain
  //);

  void solerr3_ci( int ib, int ie, int jb, int je, int kb, int ke,
		 float_sw4 h, float_sw4* __restrict__ uex,
		 float_sw4* __restrict__ u, float_sw4& li,
		 float_sw4& l2, float_sw4& xli, float_sw4 zmin, float_sw4 x0,
		 float_sw4 y0, float_sw4 z0, float_sw4 radius,
		 int imin, int imax, int jmin, int jmax, int kmin, int kmax, int geocube,
		 int i0, int i1, int j0, int j1, int k0, int k1 );
  void solerrgp_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, float_sw4 h, float_sw4* __restrict__ uex,
                   float_sw4* __restrict__ u, float_sw4& li, float_sw4& l2);
  void solerr3c_ci(int ib, int ie, int jb, int je, int kb, int ke,
                   float_sw4* __restrict__ uex, float_sw4* __restrict__ u,
                   float_sw4* __restrict__ x, float_sw4* __restrict__ y,
                   float_sw4* __restrict__ z, float_sw4* __restrict__ jac,
                   float_sw4& li, float_sw4& l2, float_sw4& xli, float_sw4 x0,
                   float_sw4 y0, float_sw4 z0, float_sw4 radius, int imin,
                   int imax, int jmin, int jmax, int kmin, int kmax, int usesg,
                   float_sw4* __restrict__ strx, float_sw4* __restrict__ stry);
  void meterr4c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                   int klast, float_sw4* __restrict__ met,
                   float_sw4* __restrict__ metex, float_sw4* __restrict__ jac,
                   float_sw4* __restrict__ jacex, float_sw4 li[5],
                   float_sw4 l2[5], int imin, int imax, int jmin, int jmax,
                   int kmin, int kmax, float_sw4 h);
  void testsrc_ci(float_sw4* __restrict__ f, int ib, int ie, int jb, int je,
                  int kb, int ke, int nk, int wind[6], float_sw4 zmin,
                  float_sw4 h, int kx[3], int ky[3], int kz[3],
                  float_sw4 mom[3]);
  void testsrcc_ci(float_sw4* __restrict__ f, int ib, int ie, int jb, int je,
                   int kb, int ke, int nk, int g, int wind[6], int kx[3],
                   int ky[3], int kz[3], float_sw4 mom[3]);
  void tw_aniso_force_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4* fo, float_sw4 t,
                         float_sw4 om, float_sw4 cv, float_sw4 ph,
                         float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                         float_sw4 phc[21], float_sw4 h, float_sw4 zmin);
  void tw_aniso_curvi_force_ci(int ifirst, int ilast, int jfirst, int jlast,
                               int kfirst, int klast,
                               float_sw4* __restrict__ fo, float_sw4 t,
                               float_sw4 om, float_sw4 cv, float_sw4 ph,
                               float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                               float_sw4 phc[21], float_sw4* __restrict__ xx,
                               float_sw4* __restrict__ yy,
                               float_sw4* __restrict__ zz);
  void tw_aniso_free_surf_z_ci(int ifirst, int ilast, int jfirst, int jlast,
                               int kfirst, int klast, int kz, float_sw4 t,
                               float_sw4 om, float_sw4 cv, float_sw4 ph,
                               float_sw4 omm, float_sw4 phc[21],
                               float_sw4* __restrict__ bforce, float_sw4 h,
                               float_sw4 zmin);
  void tw_aniso_force_tt_ci(int ifirst, int ilast, int jfirst, int jlast,
                            int kfirst, int klast, float_sw4* __restrict__ fo,
                            float_sw4 t, float_sw4 om, float_sw4 cv,
                            float_sw4 ph, float_sw4 omm, float_sw4 phm,
                            float_sw4 amprho, float_sw4 phc[21], float_sw4 h,
                            float_sw4 zmin);
  void tw_aniso_curvi_force_tt_ci(
      int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
      float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om, float_sw4 cv,
      float_sw4 ph, float_sw4 omm, float_sw4 phm, float_sw4 amprho,
      float_sw4 phc[21], float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
      float_sw4* __restrict__ zz);

  void twilightfort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* __restrict__ u, float_sw4 t,
                       float_sw4 om, float_sw4 cv, float_sw4 ph, float_sw4 h,
                       float_sw4 zmin);
  void twilightfortwind_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* __restrict__ u,
                           float_sw4 t, float_sw4 om, float_sw4 cv,
                           float_sw4 ph, float_sw4 h, float_sw4 zmin, int i1,
                           int i2, int j1, int j2, int k1, int k2);
  void twilightfortc_ci(int ifirst, int ilast, int jfirst, int jlast,
                        int kfirst, int klast, float_sw4* __restrict__ u,
                        float_sw4 t, float_sw4 om, float_sw4 cv, float_sw4 ph,
                        float_sw4* __restrict__ x, float_sw4* __restrict__ y,
                        float_sw4* __restrict__ z);
  void twilightfortatt_ci(int ifirst, int ilast, int jfirst, int jlast,
                          int kfirst, int klast, float_sw4* __restrict__ alpha,
                          float_sw4 t, float_sw4 om, float_sw4 cv, float_sw4 ph,
                          float_sw4 h, float_sw4 zmin);
  void twilightfortattc_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* __restrict__ alpha,
                           float_sw4 t, float_sw4 om, float_sw4 cv,
                           float_sw4 ph, float_sw4* __restrict__ x,
                           float_sw4* __restrict__ y,
                           float_sw4* __restrict__ z);
  void exactrhsfort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* __restrict__ fo, float_sw4 t,
                       float_sw4 om, float_sw4 c, float_sw4 ph, float_sw4 omm,
                       float_sw4 phm, float_sw4 amprho, float_sw4 ampmu,
                       float_sw4 amplambda, float_sw4 h, float_sw4 zmin);
  void exactrhsfortc_ci(int ifirst, int ilast, int jfirst, int jlast,
                        int kfirst, int klast, float_sw4* __restrict__ fo,
                        float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                        float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                        float_sw4 ampmu, float_sw4 amplambda,
                        float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
                        float_sw4* __restrict__ zz);
  void exactaccfort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* __restrict__ utt, float_sw4 t,
                       float_sw4 om, float_sw4 c, float_sw4 ph, float_sw4 h,
                       float_sw4 zmin);
  void exactaccfortc_ci(int ifirst, int ilast, int jfirst, int jlast,
                        int kfirst, int klast, float_sw4* __restrict__ utt,
                        float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                        float_sw4* __restrict__ x, float_sw4* __restrict__ y,
                        float_sw4* __restrict__ z);
  void forcingfort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                      int klast, float_sw4* __restrict__ fo, float_sw4 t,
                      float_sw4 om, float_sw4 c, float_sw4 ph, float_sw4 omm,
                      float_sw4 phm, float_sw4 amprho, float_sw4 ampmu,
                      float_sw4 amplambda, float_sw4 h, float_sw4 zmin);
  void forcingttfort_ci(int ifirst, int ilast, int jfirst, int jlast,
                        int kfirst, int klast, float_sw4* __restrict__ fo,
                        float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                        float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                        float_sw4 ampmu, float_sw4 amplambda, float_sw4 h,
                        float_sw4 zmin);
  void forcingfortc_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* __restrict__ fo, float_sw4 t,
                       float_sw4 om, float_sw4 c, float_sw4 ph, float_sw4 omm,
                       float_sw4 phm, float_sw4 amprho, float_sw4 ampmu,
                       float_sw4 amplambda, float_sw4* __restrict__ xx,
                       float_sw4* __restrict__ yy, float_sw4* __restrict__ zz);
  void forcingttfortc_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4* __restrict__ fo,
                         float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                         float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                         float_sw4 ampmu, float_sw4 amplambda,
                         float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
                         float_sw4* __restrict__ zz);
  void exactmatfort_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4* __restrict__ rho,
                       float_sw4* __restrict__ mu, float_sw4* __restrict__ la,
                       float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                       float_sw4 ampmu, float_sw4 amplambda, float_sw4 h,
                       float_sw4 zmin);
  void exactmatfortc_ci(int ifirst, int ilast, int jfirst, int jlast,
                        int kfirst, int klast, float_sw4* __restrict__ rho,
                        float_sw4* __restrict__ mu, float_sw4* __restrict__ la,
                        float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                        float_sw4 ampmu, float_sw4 amplambda,
                        float_sw4* __restrict__ x, float_sw4* __restrict__ y,
                        float_sw4* __restrict__ z);
  void exactrhsfortsg_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4* __restrict__ fo,
                         float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                         float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                         float_sw4 ampmu, float_sw4 amplambda, float_sw4 h,
                         float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
                         float_sw4 omstrz);
  void exactrhsfortsgc_ci(int ifirst, int ilast, int jfirst, int jlast,
                          int kfirst, int klast, float_sw4* __restrict__ fo,
                          float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                          float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                          float_sw4 ampmu, float_sw4 amplambda,
                          float_sw4* __restrict__ xx,
                          float_sw4* __restrict__ yy,
                          float_sw4* __restrict__ zz, float_sw4 omstrx,
                          float_sw4 omstry, float_sw4 omstrz);
  void exactmatfortatt_ci(int ifirst, int ilast, int jfirst, int jlast,
                          int kfirst, int klast, float_sw4* __restrict__ mu,
                          float_sw4* __restrict__ la, float_sw4 momega,
                          float_sw4 mphase, float_sw4 ampmu,
                          float_sw4 amplambda, float_sw4 h, float_sw4 zmin);
  void exactmatfortattc_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* __restrict__ mu,
                           float_sw4* __restrict__ la, float_sw4 momega,
                           float_sw4 mphase, float_sw4 ampmu,
                           float_sw4 amplambda, float_sw4* __restrict__ xx,
                           float_sw4* __restrict__ yy,
                           float_sw4* __restrict__ zz);
  void forcingfortatt_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4* __restrict__ fo,
                         float_sw4 t, float_sw4 omega, float_sw4 c,
                         float_sw4 phase, float_sw4 momega, float_sw4 mphase,
                         float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
                         float_sw4 h, float_sw4 zmin);
  void forcingttattfort_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* __restrict__ fo,
                           float_sw4 t, float_sw4 omega, float_sw4 c,
                           float_sw4 phase, float_sw4 momega, float_sw4 mphase,
                           float_sw4 amprho, float_sw4 ampmu,
                           float_sw4 amplambda, float_sw4 h, float_sw4 zmin);
  void addmemvarforcing_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* __restrict__ alpha,
                           float_sw4 t, float_sw4 omega, float_sw4 c,
                           float_sw4 phase, float_sw4 omegaVE, float_sw4 dt,
                           float_sw4 h, float_sw4 zmin);
  void memvarforcesurf_ci(int ifirst, int ilast, int jfirst, int jlast, int k,
                          float_sw4* __restrict__ fo, float_sw4 t,
                          float_sw4 omega, float_sw4 c, float_sw4 phase,
                          float_sw4 omegaVE, float_sw4 dt, float_sw4 h,
                          float_sw4 zmin);
  void forcingfortattc_ci(int ifirst, int ilast, int jfirst, int jlast,
                          int kfirst, int klast, float_sw4* __restrict__ fo,
                          float_sw4 t, float_sw4 omega, float_sw4 c,
                          float_sw4 phase, float_sw4 momega, float_sw4 mphase,
                          float_sw4 amprho, float_sw4 ampmu,
                          float_sw4 amplambda, float_sw4* __restrict__ xx,
                          float_sw4* __restrict__ yy,
                          float_sw4* __restrict__ zz);
  void forcingttattfortc_ci(int ifirst, int ilast, int jfirst, int jlast,
                            int kfirst, int klast, float_sw4* __restrict__ fo,
                            float_sw4 t, float_sw4 omega, float_sw4 c,
                            float_sw4 phase, float_sw4 momega, float_sw4 mphase,
                            float_sw4 amprho, float_sw4 ampmu,
                            float_sw4 amplambda, float_sw4* __restrict__ xx,
                            float_sw4* __restrict__ yy,
                            float_sw4* __restrict__ zz);
  void addmemvarforcingc_ci(int ifirst, int ilast, int jfirst, int jlast,
                            int kfirst, int klast,
                            float_sw4* __restrict__ alpha, float_sw4 t,
                            float_sw4 omega, float_sw4 c, float_sw4 phase,
                            float_sw4 omegaVE, float_sw4 dt,
                            float_sw4* __restrict__ xx,
                            float_sw4* __restrict__ yy,
                            float_sw4* __restrict__ zz);
  void memvarforcesurfc_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, int k,
                           float_sw4* __restrict__ fo, float_sw4 t,
                           float_sw4 omega, float_sw4 c, float_sw4 phase,
                           float_sw4 omegaVE, float_sw4 dt,
                           float_sw4* __restrict__ xx,
                           float_sw4* __restrict__ yy,
                           float_sw4* __restrict__ zz);

  void forcingfortsg_ci(int ifirst, int ilast, int jfirst, int jlast,
                        int kfirst, int klast, float_sw4* __restrict__ fo,
                        float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                        float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                        float_sw4 ampmu, float_sw4 amplambda, float_sw4 h,
                        float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
                        float_sw4 omstrz);
  void forcingttfortsg_ci(int ifirst, int ilast, int jfirst, int jlast,
                          int kfirst, int klast, float_sw4* __restrict__ fo,
                          float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                          float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                          float_sw4 ampmu, float_sw4 amplambda, float_sw4 h,
                          float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
                          float_sw4 omstrz);
  void forcingfortcsg_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4* __restrict__ fo,
                         float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                         float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                         float_sw4 ampmu, float_sw4 amplambda,
                         float_sw4* __restrict__ xx, float_sw4* __restrict yy,
                         float_sw4* __restrict__ zz, float_sw4 omstrx,
                         float_sw4 omstry, float_sw4 omstrz);
  void forcingttfortcsg_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* __restrict__ fo,
                           float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                           float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                           float_sw4 ampmu, float_sw4 amplambda,
                           float_sw4* __restrict__ xx, float_sw4* __restrict yy,
                           float_sw4* __restrict__ zz, float_sw4 omstrx,
                           float_sw4 omstry, float_sw4 omstrz);
  void forcingfortsgatt_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* __restrict__ fo,
                           float_sw4 t, float_sw4 omega, float_sw4 c,
                           float_sw4 phase, float_sw4 momega, float_sw4 mphase,
                           float_sw4 amprho, float_sw4 ampmu,
                           float_sw4 amplambda, float_sw4 h, float_sw4 zmin,
                           float_sw4 omstrx, float_sw4 omstry,
                           float_sw4 omstrz);
  void forcingttfortsgatt_ci(int ifirst, int ilast, int jfirst, int jlast,
                             int kfirst, int klast, float_sw4* __restrict__ fo,
                             float_sw4 t, float_sw4 omega, float_sw4 c,
                             float_sw4 phase, float_sw4 momega,
                             float_sw4 mphase, float_sw4 amprho,
                             float_sw4 ampmu, float_sw4 amplambda, float_sw4 h,
                             float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
                             float_sw4 omstrz);
  void forcingfortsgattc_ci(int ifirst, int ilast, int jfirst, int jlast,
                            int kfirst, int klast, float_sw4* __restrict__ fo,
                            float_sw4 t, float_sw4 omega, float_sw4 c,
                            float_sw4 phase, float_sw4 momega, float_sw4 mphase,
                            float_sw4 amprho, float_sw4 ampmu,
                            float_sw4 amplambda, float_sw4* __restrict__ xx,
                            float_sw4* __restrict__ yy,
                            float_sw4* __restrict__ zz, float_sw4 omstrx,
                            float_sw4 omstry, float_sw4 omstrz);
  void forcingttfortsgattc_ci(
      int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
      float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega, float_sw4 c,
      float_sw4 phase, float_sw4 momega, float_sw4 mphase, float_sw4 amprho,
      float_sw4 ampmu, float_sw4 amplambda, float_sw4* __restrict__ xx,
      float_sw4* __restrict__ yy, float_sw4* __restrict__ zz, float_sw4 omstrx,
      float_sw4 omstry, float_sw4 omstrz);

  void twfrsurfz_wind_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4 h, int kz,
                         float_sw4 t, float_sw4 omega, float_sw4 c,
                         float_sw4 phase, float_sw4* __restrict__ bforce,
                         float_sw4* __restrict__ mu,
                         float_sw4* __restrict__ lambda, float_sw4 zmin, int i1,
                         int i2, int j1, int j2);
  void twfrsurfzsg_wind_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4 h, int kz,
                           float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                           float_sw4 omstrx, float_sw4 omstry,
                           float_sw4* __restrict__ bforce,
                           float_sw4* __restrict__ mu,
                           float_sw4* __restrict__ lambda, float_sw4 zmin,
                           int i1, int i2, int j1, int j2);
  void twfrsurfz_att_wind_ci(int ifirst, int ilast, int jfirst, int jlast,
                             int kfirst, int klast, float_sw4 h, int kz,
                             float_sw4 t, float_sw4 omega, float_sw4 c,
                             float_sw4 phase, float_sw4* __restrict__ bforce,
                             float_sw4* __restrict__ mua,
                             float_sw4* __restrict__ lambdaa, float_sw4 zmin,
                             int i1, int i2, int j1, int j2);
  void twfrsurfzsg_att_wind_ci(int ifirst, int ilast, int jfirst, int jlast,
                               int kfirst, int klast, float_sw4 h, int kz,
                               float_sw4 t, float_sw4 omega, float_sw4 c,
                               float_sw4 phase, float_sw4 omstrx,
                               float_sw4 omstry, float_sw4* __restrict__ bforce,
                               float_sw4* __restrict__ mua,
                               float_sw4* __restrict__ lambdaa, float_sw4 zmin,
                               int i1, int i2, int j1, int j2);

  void tw_ani_stiff_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, float_sw4 h, float_sw4 zmin, float_sw4 omm,
                       float_sw4 phm, float_sw4 amprho,
                       float_sw4* __restrict__ a_rho, float_sw4 a_phc[21],
                       float_sw4* __restrict__ a_cm);

  void tw_ani_curvi_stiff_ci(int ifirst, int ilast, int jfirst, int jlast,
                             int kfirst, int klast, float_sw4* __restrict__ xx,
                             float_sw4* __restrict__ yy,
                             float_sw4* __restrict__ zz, float_sw4 omm,
                             float_sw4 phm, float_sw4 amprho,
                             float_sw4* __restrict__ a_rho, float_sw4 a_phc[21],
                             float_sw4* __restrict__ a_cm);

  void anisomtrltocurvilinear_ci(int ifirst, int ilast, int jfirst, int jlast,
                                 int kfirst, int klast,
                                 float_sw4* __restrict__ a_met,
                                 float_sw4* __restrict__ a_c,
                                 float_sw4* __restrict__ a_cnew);

  void velsum_ci(int is, int ie, int js, int je, int ks, int ke, int i1, int i2,
                 int j1, int j2, int k1, int k2, float_sw4* __restrict__ mu,
                 float_sw4* __restrict__ lambda, float_sw4* __restrict__ rho,
                 float_sw4& cp, float_sw4& cs, size_t& npts);

  // void enforceIC(std::vector<Sarray>& a_Up, std::vector<Sarray>& a_U,
  //                std::vector<Sarray>& a_Um, vector<Sarray*>& a_AlphaVEp,
  //                vector<Sarray*>& a_AlphaVE, vector<Sarray*>& a_AlphaVEm,
  //                float_sw4 t, bool predictor, vector<Sarray>& F,
  //                std::vector<GridPointSource*>& point_sources);
  void enforceIC(std::vector<Sarray>& a_Up, std::vector<Sarray>& a_U,
                 std::vector<Sarray>& a_Um, vector<Sarray*>& a_AlphaVEp,
                 vector<Sarray*>& a_AlphaVE, vector<Sarray*>& a_AlphaVEm,
                 float_sw4 t, bool predictor, vector<Sarray>& F,
                 std::vector<GridPointSource*>&point_sources,
                 vector<Sarray>& a_Rho, vector<Sarray>& a_Mu,
		 vector<Sarray>& a_Lambda, bool backward = false);
  void enforceIC2(std::vector<Sarray>& a_Up, std::vector<Sarray>& a_U,
                  std::vector<Sarray>& a_Um, vector<Sarray*>& a_AlphaVEp,
                  float_sw4 t, vector<Sarray>& F,
                  std::vector<GridPointSource*>& point_sources,
		  vector<Sarray>& a_Rho, vector<Sarray>& a_Mu,
		  vector<Sarray>& a_Lambda);
  void CurviCartIC(int gcart, vector<Sarray>& a_U, vector<Sarray>& a_Rho,
		   vector<Sarray>& a_Mu,
                   vector<Sarray>& a_Lambda, vector<Sarray*>& a_AlphaVE,
                   float_sw4 t);
  void dirichlet_hom_ic(Sarray& U, int g, int k, bool inner);
  void dirichlet_twilight_ic(Sarray& U, int g, int kic, float_sw4 t);

  void dirichlet_LRic(Sarray& U, int g, int kic, float_sw4 t, int adj);
  void dirichlet_LRstress(Sarray& B, int g, int kic, float_sw4 t, int adj);

  void gridref_initial_guess(Sarray& u, int g, bool upper);
  void compute_preliminary_corrector(
      Sarray& a_Up, Sarray& a_U, Sarray& a_Um, Sarray* a_AlphaVEp,
      Sarray* a_AlphaVE, Sarray* a_AlphaVEm, Sarray& Utt, Sarray& Unext, int g,
      int kic, float_sw4 t, Sarray& Ftt,
      std::vector<GridPointSource*> &point_sources, Sarray& a_Rho, Sarray& a_Mu,
      Sarray& a_Lambda);
  // void compute_preliminary_corrector( Sarray& a_Up, Sarray& a_U, Sarray&
  // a_Um,
  //                                     Sarray& Utt, Sarray& Unext,
  //                                     int g, int kic, double t,
  //                                     std::vector<GridPointSource*>
  //                                     point_sources );

   void compute_preliminary_predictor(Sarray& a_Up, Sarray& a_U,
                                     Sarray* a_AlphaVEp, Sarray& Unext, int g,
                                     int kic, float_sw4 t, Sarray& F,
                                     vector<GridPointSource*>& point_sources,
                                     Sarray& a_Rho, Sarray& a_Mu,
                                     Sarray& a_Lambda);

  void compute_icstresses(Sarray& a_Up, Sarray& B, int g, int kic,
                          float_sw4* a_str_x, float_sw4* a_str_y);
  void compute_icstresses_cpu(Sarray& a_Up, Sarray& B, int g, int kic,
                              float_sw4* a_str_x, float_sw4* a_str_y,
                              float_sw4* sbop, char op);
  void compute_icstresses2(Sarray& a_Up, Sarray& B, int kic, float_sw4 h,
                           Sarray& a_mu, Sarray& a_lambda, float_sw4* a_str_x,
                           float_sw4* a_str_y, float_sw4* sbop, char op);

  void compute_icstresses_curv(Sarray& a_Up, Sarray& B, int kic,
                               Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
                               float_sw4* a_str_x, float_sw4* a_str_y,
                               float_sw4* sbop, char op);
  void add_ve_stresses(Sarray& a_Up, Sarray& B, int g, int kic, int a_a,
                       float_sw4* a_str_x, float_sw4* a_str_y);

  void consintp(Sarray& Uf, Sarray& Unextf, Sarray& Bf, Sarray& Muf,
                Sarray& Lambdaf, Sarray& Rhof, float_sw4 hf, Sarray& Uc,
                Sarray& Unextc, Sarray& Bc, Sarray& Muc, Sarray& Lambdac,
                Sarray& Rhoc, float_sw4 hc, float_sw4 cof, int gc, int gp,
                int is_periodic[2]);

  void checkintp(Sarray& Uf, Sarray& Unextf, Sarray& Bf, Sarray& Muf,
                 Sarray& Lambdaf, Sarray& Rhof, float_sw4 hf, Sarray& Uc,
                 Sarray& Unextc, Sarray& Bc, Sarray& Muc, Sarray& Lambdac,
                 Sarray& Rhoc, float_sw4 hc, float_sw4 cof, int gc, int gf,
                 int is_periodic[2], float_sw4 t);

  void check_displacement_continuity(Sarray& Uf, Sarray& Uc, int gf, int gc);
  void check_corrector(Sarray& Uf, Sarray& Uc, Sarray& Unextf, Sarray& Unextc,
                       int kf, int kc);
  void getDtFromRestartFile();
  void initial_tw_test(vector<Sarray>& U, vector<Sarray>& Up, vector<Sarray>& F,
                       vector<Sarray>& Mu, vector<Sarray>& Lambda,
                       vector<Sarray>& Lu, vector<Sarray>& Uacc,
                       vector<Sarray*> AlphaVE,
                       vector<GridPointSource*> point_sources,
                       vector<int> identsources, float_sw4 t);
  void checkpoint_twilight_test(vector<Sarray>& Um, vector<Sarray>& U,
                                vector<Sarray>& Up, vector<Sarray*> AlphaVEm,
                                vector<Sarray*> AlphaVE,
                                vector<Sarray*> AlphaVEp,
                                vector<Source*> a_Sources, float_sw4 t);
  void make_type(vector<std::tuple<int, int, int>>& send_type,
                 vector<std::tuple<float_sw4*, float_sw4*>>& bufs_type, int i1,
                 int j1, int k1, int i2, int j2, int k2, int g);
  void make_type_2d(vector<std::tuple<int, int, int>>& send_type,
                    vector<std::tuple<float_sw4*, float_sw4*>>& bufs_type,
                    int i1, int j1, int k1, int g);
  void communicate_array_async(Sarray& u, int grid);
  void communicate_array_2d_async(Sarray& u, int g, int k);
  void communicate_array_2d_async_memo(Sarray& u, int g, int k);
  void AMPI_Sendrecv(float_sw4* a, int scount, std::tuple<int, int, int>& sendt,
                     int sentto, int stag, float_sw4* b, int rcount,
                     std::tuple<int, int, int>& recvt, int recvfrom, int rtag,
                     std::tuple<float_sw4*, float_sw4*>& buf, MPI_Comm comm,
                     MPI_Status* status);
  void AMPI_Sendrecv2(float_sw4* a, int scount,
                      std::tuple<int, int, int>& sendt, int sentto, int stag,
                      float_sw4* b, int rcount,
                      std::tuple<int, int, int>& recvt, int recvfrom, int rtag,
                      std::tuple<float_sw4*, float_sw4*>& buf, MPI_Comm comm,
                      MPI_Status* status);
  void getbuffer_device(float_sw4* data, float_sw4* buf,
                        std::tuple<int, int, int>& mtype, bool async = false);
  void putbuffer_device(float_sw4* data, float_sw4* buf,
                        std::tuple<int, int, int>& mtype, bool async = false);
  void getbuffer_host(float_sw4* data, float_sw4* buf,
                      std::tuple<int, int, int>& mtype);
  void putbuffer_host(float_sw4* data, float_sw4* buf,
                      std::tuple<int, int, int>& mtype);
  void perturb_vels(Sarray& cs, Sarray& cp, Sarray& rndpert);
  void perturb_rho(Sarray& rho, Sarray& rndpert);

  TestTwilight* create_twilight();
  TestEcons* create_energytest();
  AllDims* get_fine_alldimobject();
  TestPointSource* get_point_source_test();
  //
  // VARIABLES BEYOND THIS POINT
  //
  const float_sw4 NO_TOPO;

  // ------------------------------------------
  // Grid
  // ------------------------------------------

  int mNumberOfGrids, mNumberOfCartesianGrids;

  // grid sizes are needed by the Source and Image classes, so should be kept
  // public
  vector<float_sw4> mGridSize;

  // part of global array on each processor, including ghost points = all points
  vector<int> m_iStart, m_iEnd, m_jStart, m_jEnd, m_kStart, m_kEnd;

  // Active subcube is the part of the domain where the material is
  // variable in material inversion.
  vector<int> m_iStartAct, m_iEndAct, m_jStartAct, m_jEndAct, m_kStartAct,
      m_kEndAct;
  vector<int> m_iStartActGlobal, m_iEndActGlobal, m_jStartActGlobal,
      m_jEndActGlobal;
  vector<int> m_kStartActGlobal, m_kEndActGlobal;

  // global number of grid points on each refinement level, without ghost points
  vector<int> m_global_nx, m_global_ny, m_global_nz;

  // part of global array on each processor, excluding ghost points and parallel
  // overlap points = interior points
  vector<int> m_iStartInt, m_iEndInt, m_jStartInt, m_jEndInt, m_kStartInt,
      m_kEndInt;

  // Note that the m_paddingCells array is no longer needed to get the range of
  // internal grid points Instead use m_iStartInt[g], m_iEndInt[g], etc,
  int m_paddingCells[4];  // indexing is [0] = low-i, [1] = high-i, [2] = low-j,
                          // [3] = high-j

  // For the Cartesian grid, we only need to offset in z
  vector<float_sw4> m_zmin;  // needed by the Source and Image classes

  // for the curvilinear grid, we also store the cartesian coordinates of the
  // grid points
  vector<Sarray> mX, mY, mZ;  // needed by the Source class, so must be public
  vector<Sarray> mJ;          // Jacobian also needed by the Source class
  // and the metric derivatives as well as the jacobian
  vector<Sarray> mMetric;

  GridGenerator* m_gridGenerator;

  // command prefilter
  bool m_prefilter_sources, m_filter_observations;
  // filter setup
  // Filter for time function
  Filter* m_filter_ptr;
  // Filter for observations
  Filter* m_filterobs_ptr;
  // Test cases for optimizer, validate gradient, hessian, output function
  // surface, etc...
  int m_opttest;

  // 2-D arrays with elevation-values (=-z) as function of horizontal indices
  // mTopo holds the raw topography (according to the "elevation" field in the
  // etree) topoMat holds the highest elevation where the etree returns solid
  // material properties (now local to EtreeFile::readEFile() ) mTopoGridExt
  // holds the smoothed topography which follows th`e top surface of the
  // curvilinear grid
  Sarray mTopo, mTopoGridExt;

  // 2-D arrays with interface surfaces (z-coordinates) for mesh refinement in
  // the curvilinear grid
  vector<Sarray> m_curviInterface;

  // material description used with material surfaces and the ifile command
  vector<MaterialProperty*> m_materials;
  MPI_Comm m_cartesian_communicator, m_1d_communicator,m_cross_communicator;

  ofstream msgStream;

  // vectors of Sarrays hold material properties on all grids.
  vector<Sarray> mMu;
  vector<Sarray> mLambda;
  vector<Sarray> mRho;
  vector<Sarray*> mMuVE, mLambdaVE;  // Attenuation material
  vector<Sarray> mC;                 // Anisotropic material parameters
  Sarray mCcurv;  // Anisotropic material with metric (on curvilinear grid).

  // Store coefficeints needed for Mesh refinement
  vector<Sarray> m_Morf, m_Mlrf, m_Mufs, m_Mlfs, m_Morc, m_Mlrc, m_Mucs, m_Mlcs;

  vector<float_sw4> m_curviRefLev;

 private:
  // void preprocessSources(vector<Source*>& a_GlobalSources);
  void preprocessSources(vector<vector<Source*>>& a_GlobalSources);
  void revvector(int npts, float_sw4* v);

  int m_nevent;  // Number of events, needed for multiple event material
                 // optimization.
  int m_nevents_specified;  // Number of event lines in input file
  bool m_events_parallel;   // Process events in parallel
  int m_eStart, m_eEnd;
  map<string, int> m_event_names;

  // epicenter
  vector<float_sw4> m_epi_lat, m_epi_lon, m_epi_depth, m_epi_t0;

  // PJ *m_projection;
  // float_sw4 m_xoffset, m_yoffset;
  GeographicProjection* m_geoproj;

  ForcingTwilight* m_twilight_forcing;
  TestPointSource* m_point_source_test;
  bool m_moment_test;
  TestEnergy* m_energy_test;
  TestLamb* m_lamb_test;
  TestRayleighWave* m_rayleigh_wave_test;

  vector<MaterialData*> m_mtrlblocks;
  vector<AnisotropicMaterial*> m_anisotropic_mtrlblocks;

  // index convention: [0]: low-x, [1]: high-x, [2]: low-y, [3]: high-y; [4]:
  // low-z, [5]: high-z
  boundaryConditionType mbcGlobalType[6];  // these are the boundary conditions
                                           // for the global problem
  vector<boundaryConditionType*>
      m_bcType;  // these are the boundary conditions for each grid on the local
                 // processor, with bProcessor conditions
  float_sw4 mTstart;
  float_sw4 mDt;

  bool m_doubly_periodic;

  int m_proc_array[2];

  bool mbcsSet;

  // for some simple topographies (e.g. Gaussian hill) there is an analytical
  // expression for the top elevation
  bool m_analytical_topo, m_use_analytical_metric;
  float_sw4 m_GaussianAmp, m_GaussianLx, m_GaussianLy, m_GaussianXc,
      m_GaussianYc;

  // interface surfaces in the material model
  // int m_number_material_surfaces, m_Nlon, m_Nlat;
  // float_sw4 m_materialLonMax, m_materialLonMin, m_materialLatMax,
  // m_materialLatMin; Sarray m_materialDepth; float_sw4 *m_materialLon,
  // *m_materialLat;

  // global material thresholds
  bool m_useVelocityThresholds;
  float_sw4 m_vpMin, m_vsMin;

  // order of polynomial mapping in algebraic grid genenerator
  int m_grid_interpolation_order;
  float_sw4 m_zetaBreak;

  // metric of the curvilinear grid
  float_sw4 m_minJacobian, m_maxJacobian;

  string m_scenario;

  // command limitfrequency
  float_sw4 m_frequency_limit;
  bool m_limit_frequency;
  int m_ppw;

  // parallel io stuff
  bool m_pfs;
  int m_nwriters;

  // supergrid
  bool m_use_supergrid;
  int m_sg_gp_thickness;   //, m_sg_gp_transition;
  int m_sg_damping_order;  // 4 or 6 order dissipation operator
  float_sw4 m_supergrid_damping_coefficient;
  float_sw4 m_supergrid_width;  // width in physical units
  bool m_use_sg_width;          // use width instead of gp
  vector<SuperGrid> m_supergrid_taper_x, m_supergrid_taper_y;
  vector<SuperGrid> m_supergrid_taper_z;

  // string mPath, mObsPath, mTempPath;
  vector<string> mPath, mObsPath;  // Nevent?
  string mTempPath;

  // number of boundary points on each side
  vector<int*> m_NumberOfBCPoints;

  // ghost point index window for each side of the boundary on each grid
  vector<int*> m_BndryWindow;

  // attenuation variables (only allocated if attenuation is enabled)
  bool m_use_attenuation, m_att_use_max_frequency;
  int m_number_mechanisms;
  float_sw4 m_velo_omega, m_min_omega, m_max_omega, m_att_max_frequency,
      m_att_ppw;
  float_sw4 m_qmultiplier;

  vector<Sarray> mQp, mQs;
  // vector<Sarray*> mMuVE, mLambdaVE;
  // relaxation frequencies
  vector<float_sw4> mOmegaVE;

  // Anisotropic material
  bool m_anisotropic;

  // Randomization of the material
  bool m_randomize, m_randomize_density;
  int m_random_seed[3];
  float_sw4 m_random_dist, m_random_distz, m_random_amp, m_random_amp_grad,
      m_random_sdlimit;
  vector<RandomizedMaterial*> m_random_blocks;

  // Vectors of pointers to hold boundary forcing arrays in each grid
  // this is innner cube data for coupling with other codes
  // bool m_do_geodynbc;
  // vector<int*> m_geodyn_dims;
  // vector<Sarray> m_geodyn_data1;
  // vector<Sarray> m_geodyn_data2;
  // float_sw4 m_geodyn_origin[3], m_geodyn_h, m_geodyn_dt;
  // int m_geodyn_step, m_geodyn_maxsteps, m_geodyn_blocksize;
  //    int m_geodyn_ni, m_geodyn_nj, m_geodyn_nk, m_geodyn_faces;
  // string m_geodyn_filename;
  // ifstream m_geodynfile;
  // bool m_geodyn_iwillread;

  // with topo, zmin might be different from 0
  float_sw4 m_global_xmax, m_global_ymax, m_global_zmin, m_global_zmax;

  // number of grid points near mesh refinement boundary, for  extrapolating
  // material properties
  int mMaterialExtrapolate;

  // variables from the old FileInput class
  int m_nx_base, m_ny_base, m_nz_base;
  float_sw4 m_h_base;
  vector<bool> m_iscurvilinear;
  vector<float_sw4> m_refinementBoundaries;
  InputMode m_topoInputStyle;
  string m_topoFileName, m_topoExtFileName, m_QueryType;
  bool mTopoImageFound;
  float_sw4 m_topo_zmax;
  int m_maxIter;
  float_sw4 m_EFileResolution;

  //-------------------------------------------
  // IO data
  //-------------------------------------------
  int m_myRank, m_nProcs;

  string mName;
  // string mWPPFileName;
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
  vector<ESSI3D*> mESSI3DFiles;
  vector<SfileOutput*> mSfiles;
  bool m_iotiming;

  // time data
  vector<bool> mTimeIsSet;
  vector<float_sw4> mTmax;

  vector<int> mNumberOfTimeSteps;

  // Test modes
  int m_update_boundary_function;

  // bool mTestSource;
  // bool mTestLamb;
  // bool mTestingEnergy;
  int mOrder;
  // mCFL actual cfl. Used to determine time step in forward solver.
  // mCFLmax, maximum possible cfl. Used for limiting
  //          wave speeds during material inversion
  float_sw4 mCFL, mCFLmax;

  // info on SBP boundary operators, or not.
  vector<int*> m_onesided;
  float_sw4 m_curlcoeff, m_d4coeff, m_d4_cfl;  // these should go away

  // storage for the 1-D damping coefficients
  vector<float_sw4*> m_sg_dc_x, m_sg_dc_y, m_sg_dc_z;
  vector<float_sw4*> m_sg_str_x, m_sg_str_y, m_sg_str_z;
  vector<float_sw4*> m_sg_corner_x, m_sg_corner_y, m_sg_corner_z;

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
  float_sw4 m_saved_energy;
  string m_energy_logfile;
  vector<float_sw4> m_energy;  // *

  //-------------------------------------------
  // Measure wall clock time variables
  //-------------------------------------------
  bool m_do_timing;
  int m_timing_print;
  bool m_output_detailed_timing;
  bool m_output_load;

  int m_projection_cycle;

  bool m_checkfornan;

  // testing
  float_sw4 m_max_error[3], m_l2_error[3];

  string m_error_log_file;
  bool m_error_log, m_error_print;
  int m_inner_loop;

  //  Conservative interface
  // bool m_intp_conservative;
  bool m_mesh_refinements;
  bool m_matrices_decomposed;
  float_sw4 m_citol, m_cirelfact;
  int m_cimaxiter;

  vector<float_sw4*> m_cimat1;
  vector<float_sw4*> m_cimat2;
  vector<int*> m_ciipiv1;
  vector<int*> m_ciipiv2;

  EW(const EW&);
  EW& operator=(const EW&);

  // Geodyn coupling
  bool m_do_geodynbc;
  std::vector<float_sw4*> m_geo_usgh;  // Save ghost point
  std::vector<int*> m_geodyn_dims;
  std::vector<Sarray> m_geodyn_data1;
  std::vector<Sarray> m_geodyn_data2;
  double m_geodyn_origin[3], m_geodyn_h, m_geodyn_dt;
  int m_geodyn_step, m_geodyn_maxsteps, m_geodyn_blocksize;
  int m_geodyn_ni, m_geodyn_nj, m_geodyn_nk, m_geodyn_faces;
  std::string m_geodyn_filename;
  std::ifstream m_geodynfile;
  bool m_geodyn_iwillread, m_geodyn_past_end;

  // From wpp FileInput class
  bool m_geodynbc_found, m_geodynbc_center;
  std::string m_geodynbc_filename;
  float_sw4 m_ibc_origin[3];

  int mPrintInterval;
  // (lon, lat) origin of Grid as well as
  double mGeoAz;
  double mLonOrigin, mLatOrigin;

  // GeographicCoord mGeoCoord;
  float_sw4 mMetersPerDegree, mMetersPerLongitude;
  bool mConstMetersPerLongitude;

  // is this object ready for time-stepping?
  bool mParsingSuccessful, mIsInitialized, mSourcesOK;
  bool m_testing;

  // Parameters related to the inverse problem
  bool m_inverse_problem;  // Will we solve the inverse problem?
  bool m_iniguess_pos, m_iniguess_t0fr, m_iniguess_mom,
      m_iniguess_shifts;  // Estimate initial guess ?
  bool m_output_initial_seismograms;
  bool m_compute_scalefactors;
  int m_maxit, m_maxrestart;
  float_sw4 m_tolerance;
  float_sw4 m_scalefactors[11];
  int m_cgstepselection, m_cgvarcase;
  bool m_cgfletcherreeves, m_do_linesearch;
  bool m_opt_testing;
  int m_opt_method, m_lbfgs_m;
  bool m_zerograd_at_src, m_filter_gradient, m_zerograd_at_rec;
  int m_zerograd_pad, m_gradfilter_it, m_zerogradrec_pad;
  float_sw4 m_gradfilter_ep;
  
  // perturbations for testing
  float_sw4 m_perturb;
  int m_iperturb, m_jperturb, m_kperturb, m_pervar;

  // Number of grid points per wave length, P = min Vs/(f*h)
  vector<float_sw4> mMinVsOverH;

  int m_ext_ghost_points;
  int m_ghost_points;
  int m_ppadding;

// coefficients for boundary modified 4th order SBP operators
#if defined(ENABLE_GPU)
  float_sw4 *m_sbop, *m_acof, *m_bop, *m_bope, *m_ghcof;
  float_sw4 *m_acof_no_gp, *m_ghcof_no_gp, *m_sbop_no_gp;
#else
  float_sw4 m_sbop[6], m_acof[384], m_bop[24], m_bope[48], m_ghcof[6];
  float_sw4 m_acof_no_gp[384], m_ghcof_no_gp[6], m_sbop_no_gp[6];
#endif
  // float_sw4 m_iop[5], m_iop2[5], m_bop2[24], m_sbop[6], m_acof[384],
  // m_bop[24]; float_sw4 m_hnorm[4], m_iop[5], m_iop2[5], m_bop2[24]; // unused

  // Array of sviews used in EW::enforceBCfreeAtt2 in solve.C
  SView* viewArrayActual;

  vector<CurvilinearInterface2*> m_cli2;

  // int m_neighbor[4];
  vector<MPI_Datatype> m_send_type1;
  vector<MPI_Datatype> m_send_type3;
  vector<MPI_Datatype> m_send_type4;   // metric
  vector<MPI_Datatype> m_send_type21;  // anisotropic
  MPI_Datatype m_send_type_2dfinest[2];
  MPI_Datatype m_send_type_2dfinest_ext[2];

  // for communicating interface surfaces
  vector<MPI_Datatype> m_send_type_isurfx;
  vector<MPI_Datatype> m_send_type_isurfy;

  vector<MPI_Datatype> m_send_type_2dx;
  vector<MPI_Datatype> m_send_type_2dy;
  vector<MPI_Datatype> m_send_type_2dx3p;
  vector<MPI_Datatype> m_send_type_2dy3p;
  vector<MPI_Datatype> m_send_type_2dx1p;
  vector<MPI_Datatype> m_send_type_2dy1p;

  // Data used by async send_recv with RAJA
  vector<std::tuple<int, int, int>> send_type1;
  vector<std::tuple<int, int, int>> send_type3;
  vector<std::tuple<int, int, int>> send_type4;
  vector<std::tuple<int, int, int>> send_type21;

  vector<std::tuple<int, int, int>> send_type_2dx;
  vector<std::tuple<int, int, int>> send_type_2dy;

  vector<std::tuple<int, int, int>> send_type_2dfinest_ext;

  vector<std::tuple<float_sw4*, float_sw4*>> bufs_type1;
  vector<std::tuple<float_sw4*, float_sw4*>> bufs_type3;
  vector<std::tuple<float_sw4*, float_sw4*>> bufs_type4;
  vector<std::tuple<float_sw4*, float_sw4*>> bufs_type21;

  vector<std::tuple<float_sw4*, float_sw4*>> bufs_type_2dx;
  vector<std::tuple<float_sw4*, float_sw4*>> bufs_type_2dy;

  vector<std::tuple<float_sw4*, float_sw4*>> bufs_type_2dfinest_ext;

  Space mpi_buffer_space;

  // Arrays used for offloading parts of EW::Force
  float_sw4* ForceVector;
  float_sw4** ForceAddress;
  int* idnts;
  GridPointSource** GPS;

  // std::unordered_map<size_t, double > mpi_times,mpi_times2;
  // std::unordered_map<size_t, unsigned long long > mpi_count,mpi_count2;
#ifdef SW4_TRACK_MPI
  StatMachine<size_t, double> sm;
  StatMachine<size_t, double> sm2;
  StatMachine<int, double> coll_sm;
  StatMachine<int, double> step_sm;
  StatMachine<int, double> host_sm;
#endif
 public:
  int m_neighbor[4];
  MPI_Datatype m_mpifloat;

  bool m_topography_exists;

  // UTC time corresponding to simulation time 0.
  // bool m_utc0set, m_utc0isrefevent;
  // int m_utc0[7];
  vector<vector<int>> m_utc0;  // Nevent?

  // Error handling facility
  // ErrorChecking* m_error_checking;

  // Use C-version of computational kernels
  bool m_croutines;

  // Checkpointing and restart
  //   CheckPoint* m_restart_check_point;
  CheckPoint* m_check_point;
  bool ProfilerOn;
  void load_balance();
};

#endif
