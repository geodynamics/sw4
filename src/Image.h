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
#ifndef EW_IMAGE_H
#define EW_IMAGE_H

#include <sstream>
#include <string>
#include <vector>

#include "Parallel_IO.h"
#include "Sarray.h"
#include "Source.h"
#include "boundaryConditionTypes.h"

class EW;

class Image {
 public:
  enum ImageMode {
    NONE,
    UX,
    UY,
    UZ,
    RHO,
    LAMBDA,
    MU,
    P,
    S,
    UXEXACT,
    UYEXACT,
    UZEXACT,
    DIV,
    CURLMAG,
    DIVDT,
    CURLMAGDT,
    LAT,
    LON,
    TOPO,
    GRIDX,
    GRIDY,
    GRIDZ,
    UXERR,
    UYERR,
    UZERR,
    MAGDUDT,
    HMAGDUDT,
    HMAXDUDT,
    VMAXDUDT,
    MAG,
    HMAG,
    HMAX,
    VMAX,
    GRADRHO,
    GRADMU,
    GRADLAMBDA,
    GRADP,
    GRADS,
    QP,
    QS
  };

  static sw4_type MODES;
  static Image* nil;

  enum ImageOrientation { UNDEFINED, X, Y, Z };

  Image(EW* a_ew, float_sw4 time, float_sw4 timeInterval, sw4_type cycle,
        sw4_type cycleInterval, const std::string& filePrefix, ImageMode mode,
        ImageOrientation locationType, float_sw4 locationValue, bool doubleMode,
        bool usehdf5 = false, bool userCreated = true);

  static void setSteps(sw4_type a_steps);

  // static void setTiming(float startTime, float dt);
  // static void setGridAttributes(std::vector<double> a_gridSize      ,
  // 			      std::vector<double> a_zmin          ,
  // 			      const sw4_type           a_n_ghost_points,
  // 			      const sw4_type           a_padding       );

  void set_double(bool val = true);
  /* Here, we compute the index --in the local grid-- of the coordinate value at
     which we have to plot. For x, y, the local grid is the same as the global
     grid, but for z, k resets at each refinement boundary. */
  void computeGridPtIndex();
  void allocatePlane();

  void computeImageQuantity(std::vector<Sarray>& a_mu, sw4_type a_nComp);
  void computeImagePvel(std::vector<Sarray>& mu, std::vector<Sarray>& lambda,
                        std::vector<Sarray>& rho);
  void computeImageSvel(std::vector<Sarray>& mu, std::vector<Sarray>& rho);
  void computeImageGrid(std::vector<Sarray>& a_X, std::vector<Sarray>& a_Y,
                        std::vector<Sarray>& a_Z);
  void computeImageLatLon(std::vector<Sarray>& a_X, std::vector<Sarray>& a_Y,
                          std::vector<Sarray>& a_Z);
  void computeImageDivCurl(std::vector<Sarray>& a_Up, std::vector<Sarray>& a_U,
                           std::vector<Sarray>& a_Um, float_sw4 dt, sw4_type dminus);
  void computeImageMagdt(std::vector<Sarray>& a_Up, std::vector<Sarray>& a_Um,
                         float_sw4 dt);
  void computeImageMag(std::vector<Sarray>& a_U);
  void computeImageHmagdt(std::vector<Sarray>& a_Up, std::vector<Sarray>& a_Um,
                          float_sw4 dt);
  void computeImageHmag(std::vector<Sarray>& a_U);
  void compute_image_gradp(std::vector<Sarray>& a_gLambda,
                           std::vector<Sarray>& a_Mu,
                           std::vector<Sarray>& a_Lambda,
                           std::vector<Sarray>& a_Rho);
  void compute_image_grads(std::vector<Sarray>& a_gMu,
                           std::vector<Sarray>& a_gLambda,
                           std::vector<Sarray>& a_Mu,
                           std::vector<Sarray>& a_Rho);

  void computeImageQuantityDiff(std::vector<Sarray>& a_U,
                                std::vector<Sarray>& a_Uex, sw4_type comp);

  void output_image(sw4_type a_cycle, float_sw4 a_time, float_sw4 a_dt,
                    std::vector<Sarray>& a_Up, std::vector<Sarray>& a_U,
                    std::vector<Sarray>& a_Um, std::vector<Sarray>& a_Rho,
                    std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda,
                    std::vector<Sarray>& a_gRho, std::vector<Sarray>& a_gMu,
                    std::vector<Sarray>& a_gLambda,
                    std::vector<Source*>& a_sources, sw4_type a_dminus);

  void update_image(sw4_type a_cycle, float_sw4 a_time, float_sw4 a_dt,
                    std::vector<Sarray>& a_Up, std::vector<Sarray>& a_U,
                    std::vector<Sarray>& a_Um, std::vector<Sarray>& a_Rho,
                    std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda,
                    std::vector<Sarray>& a_gRho, std::vector<Sarray>& a_gMu,
                    std::vector<Sarray>& a_gLambda,
                    std::vector<Source*>& a_sources, sw4_type a_dminus);

  // void computeImageError(std::vector<Sarray> &a_mu, sw4_type a_nComp);

  void copy2DArrayToImage(Sarray& twoDimensionalArray);

  bool is_time_derivative() const;

  void associate_gridfiles(std::vector<Image*>& imgs);

  void writeImagePlane_2(sw4_type cycle, std::string& a_path, float_sw4 time);
  void add_grid_filenames_to_file(const char* fname);
  void add_grid_to_file(const char* fname, bool iwrite, size_t offset);
  void add_grid_to_file_hdf5(const char* fname, bool iwrite, size_t offset);

  // several curvilinear grids (MR)
  void add_grids_to_file(const char* fname, bool iwrite, size_t offset);

  bool plane_in_proc(sw4_type a_gridIndexCoarsest);
  void initializeIO();

  void computeNearestGridPoint(std::vector<sw4_type>& a_arrayIndex,
                               float_sw4 x) const;

  void update_maxes_vVelMax(std::vector<Sarray>& a_Up,
                            std::vector<Sarray>& a_Um, float_sw4 dt);
  void update_maxes_hVelMax(std::vector<Sarray>& a_Up,
                            std::vector<Sarray>& a_Um, float_sw4 dt);
  void update_maxes_vMax(std::vector<Sarray>& a_U);
  void update_maxes_hMax(std::vector<Sarray>& a_U);

  const std::string fieldSuffix(ImageMode mode) const;

  bool timeToWrite(float_sw4 time, sw4_type cycle, float_sw4 dt);
  bool timeToWrite(sw4_type cycle);
  void compute_file_suffix(std::stringstream& fileSuffix, sw4_type cycle);

  ImageOrientation getOrientation() const { return mLocationType; };

  ImageMode mMode;
  std::string mFilePrefix;
  void initializeTime(double t = 0.0);
  bool needs_mgrad() const;
  double get_write_time() { return m_write_time; };

 protected:
  void define_pio();
  // bool proc_write();

  float_sw4 mTime;
  bool m_time_done;
  float_sw4 mTimeInterval;
  float_sw4 mNextTime;
  sw4_type mWritingCycle;
  sw4_type mCycleInterval;
  // std::string mFileName;

  std::vector<std::string> mMode2Suffix;
  std::vector<std::string> mOrientationString;

  static sw4_type
      mPreceedZeros;  // number of digits for unique time step in file names

  // bool m_gridPtValueInitialized;

 private:
  Image();                 // make it impossible to call default constructor
  Image(const Image& im);  // hide copy constructor

  void computeDivergence(std::vector<Sarray>& a_U,
                         std::vector<float_sw4*>& a_div);
  void computeCurl(std::vector<Sarray>& a_U, std::vector<float_sw4*>& a_curl);

  // bool mWriting;
  // bool mReadyToWrite;

  ImageOrientation mLocationType;
  float_sw4 mCoordValue;
  std::vector<sw4_type> m_gridPtIndex;

  MPI_Comm m_mpiComm_writers;
  bool m_isDefinedMPIWriters;

  // sw4_type mNonempty;

  sw4_type mGridinfo;    // -1 = undefined, 0=Cartesian patches, 1=Curvilinear grid
                    // appended, 2=grid file names appended.
  bool mStoreGrid;  // true=append curvilinear grid to image, false=append grid
                    // file names to image.
  Image* m_gridimage;   // Curvilinear grid z-coordinate
  bool m_user_created;  // true --> This image was created from the input file

  // sw4_type m_rankWriter;
  bool m_isDefined;
  bool m_double;
  bool m_usehdf5;
  EW* mEW;
  Parallel_IO** m_pio;
  double m_write_time;

  // moved to class EW
  // sw4_type m_pfs, m_nwriters;

  std::vector<sw4_type*>
      mWindow;  // start + end indices for (i,j,k) for each grid level

  // pointers to in-plane data (only one will get allocated by
  // Image::allocatePlane()
  std::vector<double*> m_doubleField;
  std::vector<float*> m_floatField;
};

#endif
