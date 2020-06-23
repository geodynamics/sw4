// -*-c++-*-
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

#include <math.h>
#include <stdlib.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "EW.h"
#include "MaterialPfile.h"
#include "Require.h"

using namespace std;

//-----------------------------------------------------------------------
MaterialPfile::MaterialPfile(EW *a_ew, const std::string file,
                             const std::string directory, const int nstenc,
                             const float_sw4 vpmin, const float_sw4 vsmin,
                             const float_sw4 rhomin, const bool flatten,
                             const bool coords_geographic)
    : mEW(a_ew),
      m_model_file(file),
      m_model_dir(directory),
      m_nstenc(nstenc),
      m_vpmin(vpmin),
      m_vsmin(vsmin),
      m_rhomin(rhomin),
      m_coords_geographic(coords_geographic) {
  read_pfile();

  mCoversAllPoints = false;
  if (m_coords_geographic) {
    double x1, x2, x3, x4, y1, y2, y3, y4;
    mEW->computeCartesianCoord(x1, y1, m_lonmin, m_latmin);
    mEW->computeCartesianCoord(x2, y2, m_lonmin, m_latmax);
    mEW->computeCartesianCoord(x3, y3, m_lonmax, m_latmin);
    mEW->computeCartesianCoord(x4, y4, m_lonmax, m_latmax);
    m_xmin = min(x1, min(x2, min(x3, x4)));
    m_xmax = max(x1, max(x2, max(x3, x4)));
    m_ymin = min(y1, min(y2, min(y3, y4)));
    m_ymax = max(y1, max(y2, max(y3, y4)));
  }
  float_sw4 bbox[6];
  mEW->getGlobalBoundingBox(bbox);

  if (m_xmin > bbox[0] || m_ymin > bbox[2] || m_depthmin > bbox[4] ||
      m_xmax < bbox[1] || m_ymax < bbox[3] || m_depthmax < bbox[5]) {
    mCoversAllPoints = false;
    //      cout << "This block does NOT cover all grid points" << endl;
  } else {
    mCoversAllPoints = true;
    // tmp
    //     cout << "This block COVERS all grid points" << endl;
  }
}

//-----------------------------------------------------------------------
void MaterialPfile::set_material_properties(std::vector<Sarray> &rho,
                                            std::vector<Sarray> &cs,
                                            std::vector<Sarray> &cp,
                                            std::vector<Sarray> &qs,
                                            std::vector<Sarray> &qp) {
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  // the qs and qp arrays are always allocated to allow qs[g].is_defined() to be
  // called
  //  bool use_attenuation = !(qs.size()==0);
  //  bool use_attenuation = m_ew->usingAttenuation();
  //  bool use_attenuation = false;

  int outside = 0;
  int material = 0;

  if (myRank == 0)
    cout << "Assigning material properties from pfile data..." << endl;

  // tmp
  //  printf("MPF: set_mat_prop: use_attenuation=%i\n",
  //  mWPP->usingAttenuation());

  //  float_sw4 x, y, z, depth;
  //  double lon, lat;

  //  bool foundcrust;

  int g, kLow, topLevel = mEW->mNumberOfGrids - 1;
  // first deal with the Cartesian grids
  for (g = 0; g < mEW->mNumberOfCartesianGrids; g++) {
    kLow = mEW->m_kStart[g];
    if (g == topLevel)  // no topography, so k=1 is at the top surface
    {
      kLow = 1;
    }
    if (m_coords_geographic) {
//       for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k)
#pragma omp parallel for reduction(+ : material, outside)
      for (int k = kLow; k <= mEW->m_kEnd[g]; ++k)
        for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j)
          for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) {
            float_sw4 x = (i - 1) * mEW->mGridSize[g];
            float_sw4 y = (j - 1) * mEW->mGridSize[g];
            float_sw4 z = mEW->m_zmin[g] + (k - 1) * mEW->mGridSize[g];
            double lon, lat;
            mEW->computeGeographicCoord(x, y, lon, lat);
            float_sw4 depth;
            mEW->getDepth(x, y, z, depth);

            if (inside(lat, lon, depth)) {
              //---------------------------------------------------------
              // Query the location...
              //---------------------------------------------------------
              // tmp
              //		  bool debug = (i==15 && j==15 && k==1);

              float_sw4 vp, vs, density, qup, qus;
              sample_latlon(lat, lon, depth, vp, vs, density, qup, qus, false);
              rho[g](i, j, k) = density;
              cp[g](i, j, k) = vp;
              cs[g](i, j, k) = vs;
              if (m_qf) {
                if (qp[g].is_defined()) qp[g](i, j, k) = qup;
                if (qs[g].is_defined()) qs[g](i, j, k) = qus;
              }
              material++;
            } else {
              if (mEW->getVerbosity() > 2) {
                printf(
                    "Point (i,j,k)=(%i, %i, %i) in grid g=%i\n"
                    "with (x,y,z)=(%e,%e,%e) and lat=%e, lon=%e, depth=%e\n"
                    "is outside the pfile domain: %e<= lat <= %e, %e <= lon <= "
                    "%e, %e <= depth <= %e\n",
                    i, j, k, g, x, y, z, lat, lon, depth, m_latmin, m_latmax,
                    m_lonmin, m_lonmax, m_depthmin, m_depthmax);
              }
              outside++;
            }

          }  // end for i, j, k
    } else {
      // Cartesian p-file
#pragma omp parallel for reduction(+ : material, outside)
      for (int k = kLow; k <= mEW->m_kEnd[g]; ++k)
        for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j)
          for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) {
            float_sw4 x = (i - 1) * mEW->mGridSize[g];
            float_sw4 y = (j - 1) * mEW->mGridSize[g];
            float_sw4 z = mEW->m_zmin[g] + (k - 1) * mEW->mGridSize[g];
            float_sw4 depth;
            mEW->getDepth(x, y, z, depth);

            if (inside_cart(x, y, depth)) {
              //---------------------------------------------------------
              // Query the location...
              //---------------------------------------------------------
              // tmp
              //		  bool debug = (i==15 && j==15 && k==1);
              float_sw4 vp, vs, density, qup, qus;
              sample_cart(x, y, depth, vp, vs, density, qup, qus, false);

              rho[g](i, j, k) = density;
              cp[g](i, j, k) = vp;
              cs[g](i, j, k) = vs;

              //		   cout << "x= " << x << " y= " << y << " depth=
              //" << depth << " vp = " << vp << " vs = " << vs << " rho = " <<
              // density << endl;

              if (m_qf) {
                if (qp[g].is_defined()) qp[g](i, j, k) = qup;
                if (qs[g].is_defined()) qs[g](i, j, k) = qus;
              }
              material++;
            } else
              outside++;
          }  // end for i, j, k

    }  // end cartesian pfile case

    // communicate material properties to ghost points (necessary on refined
    // meshes because ghost points don't have a well defined depth/topography)
    mEW->communicate_array(rho[g], g);
    mEW->communicate_array(cs[g], g);
    mEW->communicate_array(cp[g], g);
    if (m_qf) {
      if (qp[g].is_defined()) mEW->communicate_array(qp[g], g);
      if (qs[g].is_defined()) mEW->communicate_array(qs[g], g);
    }

  }  // end for all Cartesian grids

  // Now, the curvilinear grid
  if (mEW->topographyExists()) {
    int gTop = mEW->mNumberOfGrids - 1;
    for (int g = mEW->mNumberOfCartesianGrids; g < mEW->mNumberOfGrids; g++) {
      //    g = mEW->mNumberOfGrids-1;

      // the curvilinear grid is always on top
      if (g == gTop)
        kLow = 1;
      else
        kLow = mEW->m_kStart[g];

      if (m_coords_geographic) {
//       for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k)
#pragma omp parallel for reduction(+ : material, outside)
        for (int k = kLow; k <= mEW->m_kEnd[g];
             ++k)  // don't attempt querying the pfile above the topography
                   // (start at k=1 for top grid)
          for (int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; ++j)
            for (int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g]; ++i) {
              float_sw4 x = mEW->mX[g](i, j, k);
              float_sw4 y = mEW->mY[g](i, j, k);
              float_sw4 z = mEW->mZ[g](i, j, k);
              double lon, lat;
              mEW->computeGeographicCoord(x, y, lon, lat);
              float_sw4 depth = z - mEW->mZ[gTop](i, j, 1);

              if (inside(lat, lon, depth)) {
                //---------------------------------------------------------
                // Query the location...
                //---------------------------------------------------------
                float_sw4 vp, vs, density, qup, qus;
                sample_latlon(lat, lon, depth, vp, vs, density, qup, qus,
                              false);
                rho[g](i, j, k) = density;
                cp[g](i, j, k) = vp;
                cs[g](i, j, k) = vs;
                if (m_qf) {
                  if (qp[g].is_defined()) qp[g](i, j, k) = qup;
                  if (qs[g].is_defined()) qs[g](i, j, k) = qus;
                }
                material++;
              } else {
                if (mEW->getVerbosity() > 2) {
                  printf(
                      "Point (i,j,k)=(%i, %i, %i) in grid g=%i\n"
                      "with (x,y,z)=(%e,%e,%e) and lat=%e, lon=%e, depth=%e\n"
                      "is outside the pfile domain: %e<= lat <= %e, %e <= lon "
                      "<= %e, %e <= depth <= %e\n",
                      i, j, k, g, x, y, z, lat, lon, depth, m_latmin, m_latmax,
                      m_lonmin, m_lonmax, m_depthmin, m_depthmax);
                }
                outside++;
              }
            }  // end for i,j,k
      }        // end coords_geographic
      else {
//       for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) // don't
//       attempt querying the pfile above the topography (start at k=1)
#pragma omp parallel for reduction(+ : material, outside)
        for (int k = kLow; k <= mEW->m_kEnd[g];
             ++k)  // don't attempt querying the pfile above the topography
                   // (start at k=1)
          for (int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; ++j)
            for (int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g]; ++i) {
              float_sw4 x = mEW->mX[g](i, j, k);
              float_sw4 y = mEW->mY[g](i, j, k);
              float_sw4 z = mEW->mZ[g](i, j, k);
              float_sw4 depth = z - mEW->mZ[gTop](i, j, 1);
              if (inside_cart(x, y, depth)) {
                //---------------------------------------------------------
                // Query the location...
                //---------------------------------------------------------
                float_sw4 vp, vs, density, qup, qus;
                sample_cart(x, y, depth, vp, vs, density, qup, qus, false);
                rho[g](i, j, k) = density;
                cp[g](i, j, k) = vp;
                cs[g](i, j, k) = vs;
                if (m_qf) {
                  if (qp[g].is_defined()) qp[g](i, j, k) = qup;
                  if (qs[g].is_defined()) qs[g](i, j, k) = qus;
                }
                material++;
              } else
                outside++;
            }  // end for i,j,k
      }
    }  // end for g (curvilinear)

  }  // end if topographyExists()

  //  extrapolation is now done in WPP2:set_materials()

  int outsideSum, materialSum;
  MPI_Reduce(&outside, &outsideSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (mEW->proc_zero())
    cout << "outside = " << outsideSum << ", "
         << "material = " << materialSum << endl;
}

//-----------------------------------------------------------------------
void MaterialPfile::read_pfile() {
  //   int kk, m;

  //   m_qf = false;

  //   m_model_dir  = ppdir;
  //   m_model_file = ppfile;
  string ppmfile = m_model_dir + "/" + m_model_file;

  //   m_vpmin   = vpmin;
  //   m_vsmin   = vsmin;
  //   m_rhomin  = rhomin;

  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // Open file
  FILE *fd = fopen(ppmfile.c_str(), "r");
  if (fd == NULL) {
    if (myRank == 0)
      cerr << "Unable to open the pfile input file: '" << ppmfile << "'"
           << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int bufsize = 1024;
  int nread;
  char *buf = new char[bufsize];

  // Read pfile header

  // Line 1
  CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
              "Error line 1 in pfile header not found\n");
  string tok0 = strtok(buf, " \t");
  // strip off any white space
  size_t nWhite = tok0.find_first_of(" \t\n");
  m_model_name = tok0.substr(0, nWhite);

  // Line 2
  CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
              "Error line 2 in pfile header not found\n");
  nread = sscanf(buf, "%le", &m_h);
  CHECK_INPUT(nread == 1, "Error reading 2nd line of header, nread= "
                              << nread << " but expected 1\n");

  // Line 3
  CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
              "Error line 3 in pfile header not found\n");
  nread = sscanf(buf, "%i %le %le", &m_nlat, &m_latmin, &m_latmax);
  CHECK_INPUT(nread == 3, "Error reading 3rd line of header, nread= "
                              << nread << " but expected 3\n");

  // Line 4
  CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
              "Error line 4 in pfile header not found\n");
  nread = sscanf(buf, "%i %le %le", &m_nlon, &m_lonmin, &m_lonmax);
  CHECK_INPUT(nread == 3, "Error reading 4th line of header, nread= "
                              << nread << " but expected 3\n");

  // Line 5
  CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
              "Error line 5 in pfile header not found\n");
  nread = sscanf(buf, "%i %le %le", &m_nmaxdepth, &m_depthmin, &m_depthmax);
  CHECK_INPUT(nread == 3, "Error reading 5th line of header, nread= "
                              << nread << " but expected 3\n");

  // Line 6
  CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
              "Error line 6 in pfile header not found\n");
  nread = sscanf(buf, "%i %i %i %i", &m_ksed, &m_kmoho, &m_k410, &m_k660);
  CHECK_INPUT(nread == 4, "Error reading 6th line of header, nread= "
                              << nread << " but expected 4\n");

  // Line 7
  CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
              "Error line 7 in pfile header not found\n");
  char *tok = strtok(buf, " \t");
  CHECK_INPUT(tok != NULL,
              "Error on line 7 in pfile header, no Q-available flag\n");
  string cqf0 = tok;
  // strip off any white space
  nWhite = cqf0.find_first_of(" \t\n");
  string cqf = cqf0.substr(0, nWhite);

  // test
  //   printf("Q-flag string '%s'\n", cqf.c_str());

  m_qf = (cqf == "T") || (cqf == "t") || (cqf == ".TRUE.") || (cqf == ".true.");

  // done reading the header

  // make sure the stencil width does not exceed the number of grid points
  if (m_nstenc > m_nlat || m_nstenc > m_nlon) {
    m_nstenc = min(m_nlat, m_nlon);
    if (mEW->proc_zero())
      printf("Warning: pfile: stencil width reduced to %i\n", m_nstenc);
  }

  float_sw4 km = 1000;
  if (m_coords_geographic) {
    m_depthmin = m_depthmin * km;
    m_depthmax = m_depthmax * km;
    m_dlat = (m_latmax - m_latmin) / (m_nlat - 1);
    m_dlon = (m_lonmax - m_lonmin) / (m_nlon - 1);
  } else {
    m_nx = m_nlat;
    m_ny = m_nlon;
    m_xmin = m_latmin;
    m_xmax = m_latmax;
    m_ymin = m_lonmin;
    m_ymax = m_lonmax;
  }

  if (myRank == 0) {
    cout << "Pfile model name (string): '" << m_model_name << "'" << endl;

    if (!m_coords_geographic) {
      cout << "Step size in x and y: " << m_h << endl;
      cout << "Number of x-direction points: " << m_nx << endl;
      cout << "Min x: " << m_xmin << " Max x: " << m_xmax << endl;
      cout << "Number of y-direction points: " << m_ny << endl;
      cout << "Min y: " << m_ymin << " Max y: " << m_ymax << endl;
    } else {
      cout << "Grid length scale (Delta): " << m_h << endl;
      cout << "Number of latitude points: " << m_nlat << endl;
      cout << "Min Lat: " << m_latmin << " Max Lat: " << m_latmax << endl;
      cout << "Number of longitude points: " << m_nlon << endl;
      cout << "Min Lon: " << m_lonmin << " Max Lon: " << m_lonmax << endl;
    }
    cout << "Number of depth points: " << m_nmaxdepth << endl;
    cout << "Min depth: " << m_depthmin << " Max depth: " << m_depthmax << endl;
    cout << "Optional indices: Sediment: " << m_ksed << " MoHo: " << m_kmoho
         << " 410: " << m_k410 << " 660: " << m_k660 << endl;
    cout << "Attenuation Q-factors available: " << (m_qf ? "yes" : "no")
         << endl;
  }

  // Allocate arrays
  if (m_coords_geographic) {
    m_lat = new double[m_nlat];
    m_lon = new double[m_nlon];
    // new 3-D Sarrays
    mZ.define(m_nlon, m_nlat, m_nmaxdepth);
    mVp.define(m_nlon, m_nlat, m_nmaxdepth);
    mVs.define(m_nlon, m_nlat, m_nmaxdepth);
    mRho.define(m_nlon, m_nlat, m_nmaxdepth);
    if (m_qf) {
      // new 3-D Sarrays
      mQp.define(m_nlon, m_nlat, m_nmaxdepth);
      mQs.define(m_nlon, m_nlat, m_nmaxdepth);
    }
  } else {
    m_x = new float_sw4[m_nx];
    m_y = new float_sw4[m_ny];
    // new 3-D Sarrays
    mZ.define(m_nx, m_ny, m_nmaxdepth);
    mVp.define(m_nx, m_ny, m_nmaxdepth);
    mVs.define(m_nx, m_ny, m_nmaxdepth);
    mRho.define(m_nx, m_ny, m_nmaxdepth);
    if (m_qf) {
      // new 3-D Sarrays
      mQp.define(m_nx, m_ny, m_nmaxdepth);
      mQs.define(m_nx, m_ny, m_nmaxdepth);
    }
  }
  if (!m_qf) {
    if (myRank == 0) printf("ppmod: NOT allocating arrays for Qp and Qs\n");
  }

  int kk, ndepth, line = 7;

  float_sw4 zc, vp, vs, rho, qp, qs;

  if (!m_coords_geographic)  // cartesian
  {
    // note: nx = nlat, ny = nlon
    for (int jy = 0; jy < m_ny; jy++)
      for (int ix = 0; ix < m_nx; ix++) {
        CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
                    "Error in pfile profile header at coordinate "
                        << ix << " " << jy << "\n");
        nread = sscanf(buf, "%le %le %i", &m_x[ix], &m_y[jy], &ndepth);
        CHECK_INPUT(nread == 3, "Error reading 1st line of profile at "
                                    << ix << " " << jy << " nread= " << nread
                                    << " but expected 3\n");
        line++;

        // fundamental sanity checks
        if (!(m_y[jy] <= m_ymax && m_y[jy] >= m_ymin && m_x[ix] <= m_xmax &&
              m_x[ix] >= m_xmin)) {
          printf(
              "Error reading pfile: x profile #%i, y profile #%i: x=%e or y=%e "
              "out of bounds!"
              " min(x)=%e, max(x)=%e, min(y)=%e, max(y)=%e\n",
              ix + 1, jy + 1, m_x[ix], m_y[jy], m_xmin, m_xmax, m_ymin, m_ymax);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // sanity check 2
        if (ndepth != m_nmaxdepth) {
          if (myRank == 0) {
            cerr << "pfile reader error, ppmod file line=" << line << endl;
            cerr << "read ndepth=" << ndepth
                 << " which is different from header nmaxdepth=" << m_nmaxdepth
                 << endl;
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // check the range
        float_sw4 y = m_ymin + jy * m_h;
        float_sw4 x = m_xmin + ix * m_h;
        if (fabs(y - m_y[jy]) + fabs(x - m_x[ix]) > 0.1 * m_h) {
          if (myRank == 0) {
            cerr << "pfile reader error, ppmod file line=" << line << endl;
            cerr << "read x[" << ix << "]=" << m_x[ix]
                 << " but expected x=" << x << endl;
            cerr << "read y[" << jy << "]=" << m_y[jy]
                 << " but expected y=" << y << endl;
            cerr << "CHECK THE PPMOD FILE." << endl;
            cerr << "DEPTH PROFILES SHOULD BE ORDERED SUCH THAT X VARIES THE "
                    "FASTEST"
                 << endl;
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Read depth profile
        for (int k = 0; k < m_nmaxdepth; k++) {
          CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
                      "Error in pfile profile at coordinate "
                          << ix << " " << jy << " " << k << "\n");

          if (m_qf) {
            nread = sscanf(buf, "%i %le %le %le %le %le %le", &kk, &zc, &vp,
                           &vs, &rho, &qp, &qs);
            CHECK_INPUT(nread == 7, "Error reading pfile at "
                                        << ix << " " << jy << " " << k
                                        << " nread= " << nread
                                        << " but expected 7\n");
          } else {
            nread = sscanf(buf, "%i %le %le %le %le", &kk, &zc, &vp, &vs, &rho);
            CHECK_INPUT(nread == 5, "Error reading pfile at "
                                        << ix << " " << jy << " " << k
                                        << " nread= " << nread
                                        << " but expected 5\n");
          }

          mZ(ix + 1, jy + 1, k + 1) = zc;
          mVp(ix + 1, jy + 1, k + 1) = vp;
          mVs(ix + 1, jy + 1, k + 1) = vs;
          mRho(ix + 1, jy + 1, k + 1) = rho;

          if (m_qf) {
            mQp(ix + 1, jy + 1, k + 1) = qp;
            mQs(ix + 1, jy + 1, k + 1) = qs;
          }

          line++;

          // do we need to cap the values here?
          // m_vp[m] = max(m_vp[m], m_vpmin );
          // m_vs[m] = max(m_vs[m], m_vsmin );
          // m_rho[m]= max(m_rho[m],m_rhomin);
        }  // end for k
      }
  } else  // geographic coordinates
  {
    for (int j = 0; j < m_nlat; j++)
      for (int i = 0; i < m_nlon; i++) {
        CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
                    "Error in pfile profile header at coordinate "
                        << i << " " << j << "\n");

        nread = sscanf(buf, "%le %le %i", &m_lat[j], &m_lon[i], &ndepth);
        CHECK_INPUT(nread == 3, "Error reading 1st line of profile at "
                                    << i << " " << j << " nread= " << nread
                                    << " but expected 3\n");

        line++;

        // fundamental sanity checks
        if (!(m_lat[j] <= m_latmax && m_lat[j] >= m_latmin &&
              m_lon[i] <= m_lonmax && m_lon[i] >= m_lonmin)) {
          printf(
              "Error reading pfile: lat profile #%i, lon profile #%i: lat=%e "
              "or lon=%e out of bounds!"
              " min(lat)=%e, max(lat)=%e, min(lon)=%e, max(lon)=%e\n",
              j + 1, i + 1, m_lat[j], m_lon[i], m_latmin, m_latmax, m_lonmin,
              m_lonmax);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // sanity check 2
        if (ndepth != m_nmaxdepth) {
          if (myRank == 0) {
            cerr << "pfile reader error, ppmod file line=" << line << endl;
            cerr << "read ndepth=" << ndepth
                 << " which is different from header nmaxdepth=" << m_nmaxdepth
                 << endl;
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // check the coordinate (now allowing for different step sizes in lat
        // and lon)
        double lat_j = m_latmin + j * m_dlat;
        double lon_i = m_lonmin + i * m_dlon;
        if (fabs(lat_j - m_lat[j]) > 0.05 * m_dlat ||
            fabs(lon_i - m_lon[i]) > 0.05 * m_dlon) {
          if (myRank == 0) {
            cerr << "pfile reader error, ppmod file line=" << line << endl;
            cerr << "read lon_ppm[" << i << "]=" << m_lon[i]
                 << " but expected lon=" << lon_i << endl;
            cerr << "read lat_ppm[" << j << "]=" << m_lat[j]
                 << " but expected lat=" << lat_j << endl;
            cerr << "CHECK THE PPMOD FILE." << endl;
            cerr << "DEPTH PROFILES SHOULD BE ORDERED SUCH THAT LONGITUDE "
                    "VARIES THE FASTEST"
                 << endl;
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Read depth profile
        for (int k = 0; k < m_nmaxdepth; k++) {
          CHECK_INPUT(fgets(buf, bufsize, fd) != NULL,
                      "Error reading pfile buffer at " << i << " " << j << " "
                                                       << k << "\n");

          if (m_qf) {
            nread = sscanf(buf, "%i %le %le %le %le %le %le", &kk, &zc, &vp,
                           &vs, &rho, &qp, &qs);
            CHECK_INPUT(nread == 7, "Error reading pfile at "
                                        << i << " " << j << " " << k
                                        << " nread= " << nread
                                        << " but expected 7\n");
          } else {
            nread = sscanf(buf, "%i %le %le %le %le", &kk, &zc, &vp, &vs, &rho);
            CHECK_INPUT(nread == 5, "Error reading pfile at "
                                        << i << " " << j << " " << k
                                        << " nread= " << nread
                                        << " but expected 5\n");
          }

          mZ(i + 1, j + 1, k + 1) = km * zc;
          mVp(i + 1, j + 1, k + 1) = km * vp;
          mVs(i + 1, j + 1, k + 1) = km * vs;
          mRho(i + 1, j + 1, k + 1) = km * rho;

          if (m_qf) {
            mQp(i + 1, j + 1, k + 1) = qp;
            mQs(i + 1, j + 1, k + 1) = qs;
          }

          line++;
          // are these needed anymore???
          // m_vp[m] = max(m_vp[m], m_vpmin );
          // m_vs[m] = max(m_vs[m], m_vsmin );
          // m_rho[m]= max(m_rho[m],m_rhomin);
        }  // end for k
      }
  }
  // closing the pfile
  fclose(fd);

  if (myRank == 0) {
    cout << "******* Done reading Pfile **********" << endl << endl;
  }
  // tmp
  //   cout << "******* Done reading Pfile, proc=" << myRank << endl;

  delete[] buf;
}

//---------------------------------------------------------------------------------------
// int MaterialPfile::get_material_pt( double x, double y, double z, double&
// rho, double& cs, double& cp, 				    double& qs,
// double& qp
// )
//  {
//   int retval = 0;
//   double zsed, zmoho;
//   bool foundcrust;

//   if( m_coords_geographic )
//   {
//     double lon, lat, depth;

//     mEW->computeGeographicCoord( x, y, lon, lat );
//     mEW->getDepth(x,y,z,depth);

//     if( inside( lat, lon, depth )  )
//     {
//---------------------------------------------------------
// Query the location...
//---------------------------------------------------------
//       sample_latlon( lat, lon, depth, cp, cs, rho, qp, qs, false );
//     }
//     else
//       retval = -1;
//   }
//   else
//   {
//     if( inside_cart( x, y, z )  ) // elev = -depth
//     {
//---------------------------------------------------------
// Query the location...
//---------------------------------------------------------
//       sample_cart( x, y, z, cp, cs, rho, qp, qs );
//     }
//     else
//       retval = -1;
//   }
//   return retval;
//  }

//-----------------------------------------------------------------------
void MaterialPfile::sample_latlon(double lats, double lons, float_sw4 zs,
                                  float_sw4 &vp, float_sw4 &vs, float_sw4 &rho,
                                  float_sw4 &qp, float_sw4 &qs, bool debug)
//--------------------------------------------------------------------------
// return material properties (vp, vs, rho) at point (lats, lons, zs)
//--------------------------------------------------------------------------
{
  // tmp
  // if (debug) printf("DEBUG::sample_latlon: lats=%e, lons=%e, zs=%e, dlon=%e,
  // dlat=%e, m_h=%e\n", lats, lons, zs, 		    m_dlon, m_dlat,
  // m_h);

  //  Check if lats and lons are out of range
  if (lats < m_latmin) {
    cerr << "MaterialPfile::sample lats out of range (min): " << lats << ", "
         << m_latmin << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (lats > m_latmax) {
    cerr << "MaterialPfile::sample lats out of range (max): " << lats << ", "
         << m_latmax << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (lons < m_lonmin) {
    cerr << "MaterialPfile::sample lons out of range (min): " << lons << ", "
         << m_lonmin << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (lons > m_lonmax) {
    cerr << "MaterialPfile::sample lons out of range (max): " << lons << ", "
         << m_lonmax << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int ii, jj;
  if (m_nstenc % 2 == 1) {
    // odd number of points in stencil
    int s = (m_nstenc - 1) / 2;  // half stencil width
    ii = static_cast<int>(floor((lons - m_lonmin) / m_dlon + 0.5)) -
         s;  // updated to use different step lengths in lat and lon
    jj = static_cast<int>(floor((lats - m_latmin) / m_dlat + 0.5)) - s;
  } else {
    // even number of points in stencil
    int s = m_nstenc / 2;
    ii = static_cast<int>(floor((lons - m_lonmin) / m_dlon)) - s +
         1;  // updated to use different step lengths in lat and lon
    jj = static_cast<int>(floor((lats - m_latmin) / m_dlat)) - s + 1;
  }

  // make sure we stay within array boundaries
  // min indices
  if (ii < 0) ii = 0;
  if (jj < 0) jj = 0;

  // max indices
  int ii2 = ii + m_nstenc - 1;
  int jj2 = jj + m_nstenc - 1;

  if (ii2 >= m_nlon) {
    ii2 = m_nlon - 1;
    ii = ii2 - (m_nstenc - 1);
  }

  if (jj2 >= m_nlat) {
    jj2 = m_nlat - 1;
    jj = jj2 - (m_nstenc - 1);
  }

  // if (debug)
  // {
  //   printf("pfile: lon array:\n");
  //   for (int q=ii; q<=ii2; q++)
  //     printf("lon[%i]=%e\n", q, m_lon[q]);

  //   printf("pfile: lat array:\n");
  //   for (int q=jj; q<=jj2; q++)
  //     printf("lat[%i]=%e\n", q, m_lat[q]);
  // }

  float_sw4 w = 0;
  vp = vs = rho = qp = qs = 0;
  float_sw4 appm = 0.5 * m_nstenc * m_h / sqrt(-log(1e-6));
  float_sw4 appmi2 = 1.0 / (appm * appm);

  for (int j1 = jj; j1 <= jj2; j1++)
    for (int i1 = ii; i1 <= ii2; i1++) {
      float_sw4 wgh = exp(-((lons - m_lon[i1]) * (lons - m_lon[i1]) +
                            (lats - m_lat[j1]) * (lats - m_lat[j1])) *
                          appmi2);
      w += wgh;

      // depth index
      int kk;
      for (kk = 1; kk < m_nmaxdepth; kk++)  // AP changed from kk <= m_nmaxdepth
      {
        if (mZ(i1 + 1, j1 + 1, kk) > zs) break;
      }
      // at this point we should have mZ(kk-1) <= zs < mZ(kk), kk < m_nmaxdepth

      int k1 = kk - 1;
      // now we should have mZ(k1) <= zs < mZ(k1+1)
      if (k1 <= 0) k1 = 1;

      // linear interpolation factor ( what happens if two mZ values are
      // identical? )
      //	  double factor =
      //(zs-mZ(i1+1,j1+1,k1))/(mZ(i1+1,j1+1,k1+1)-mZ(i1+1,j1+1,k1)); if( factor
      //< 0 ) 	     factor = 0;
      float_sw4 dz = mZ(i1 + 1, j1 + 1, k1 + 1) - mZ(i1 + 1, j1 + 1, k1);
      float_sw4 factor = 0;
      if (dz != 0) {
        factor = (zs - mZ(i1 + 1, j1 + 1, k1)) / dz;
        if (factor < 0) factor = 0;
        if (factor > 1) factor = 1;
      }

      vp += (mVp(i1 + 1, j1 + 1, k1) +
             factor * (mVp(i1 + 1, j1 + 1, k1 + 1) - mVp(i1 + 1, j1 + 1, k1))) *
            wgh;
      vs += (mVs(i1 + 1, j1 + 1, k1) +
             factor * (mVs(i1 + 1, j1 + 1, k1 + 1) - mVs(i1 + 1, j1 + 1, k1))) *
            wgh;
      rho +=
          (mRho(i1 + 1, j1 + 1, k1) +
           factor * (mRho(i1 + 1, j1 + 1, k1 + 1) - mRho(i1 + 1, j1 + 1, k1))) *
          wgh;
      // tmp
      // if (debug) printf("DEBUG: i1+1=%i, j1+1=%i, k1=%i, vp=%e, wgh=%e\n",
      // i1+1, j1+1, k1, 		    mVp(i1+1,j1+1,k1) +
      // factor*(mVp(i1+1,j1+1,k1+1)-mVp(i1+1,j1+1,k1)), wgh);

      if (m_qf) {
        qp += (mQp(i1 + 1, j1 + 1, k1) + factor * (mQp(i1 + 1, j1 + 1, k1 + 1) -
                                                   mQp(i1 + 1, j1 + 1, k1))) *
              wgh;
        qs += (mQs(i1 + 1, j1 + 1, k1) + factor * (mQs(i1 + 1, j1 + 1, k1 + 1) -
                                                   mQs(i1 + 1, j1 + 1, k1))) *
              wgh;
      }
    }  // end for j1, i1

  // Now compute average properties by distance-weighted Gaussian average
  float_sw4 iw;
  // at this point, w holds the sum of the weigths
  if (w != 0.)
    iw = 1.0 / w;
  else {
    printf(
        "Error MaterialPfile::sample_latlon: weight w = 0 at lat=%e, lon=%e, "
        "depth=%e\n",
        lats, lons, zs);
    // tmp
    printf("ii=%i, ii2=%i, jj=%i, jj2=%i, lon_ppm[ii]=%e, lat_ppm[jj]=%e\n", ii,
           ii2, jj, jj2, m_lon[ii], m_lat[jj]);
    float_sw4 dist2 = (lons - m_lon[ii]) * (lons - m_lon[ii]) +
                      (lats - m_lat[jj]) * (lats - m_lat[jj]);
    float_sw4 exponent = -dist2 * appmi2;
    float_sw4 wgh = exp(exponent);
    printf("dist2=%e, appm=%e, appmi2=%e, exponent=%e, wgh=%e\n", dist2, appm,
           appmi2, exponent, wgh);

    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  vp *= iw;
  vs *= iw;
  rho *= iw;
  if (m_qf) {
    qp *= iw;
    qs *= iw;
  }
}

//-----------------------------------------------------------------------
void MaterialPfile::sample_cart(float_sw4 xs, float_sw4 ys, float_sw4 zs,
                                float_sw4 &vp, float_sw4 &vs, float_sw4 &rho,
                                float_sw4 &qp, float_sw4 &qs, bool debug)
//--------------------------------------------------------------------------
// return material properties (vp, vs, rho) at point (xs, ys, zs)
//--------------------------------------------------------------------------
{
  // tmp
  //  if (debug) printf("DEBUG::sample_cart: xs=%e, ys=%e, zs=%e, m_h=%e\n", xs,
  //  ys, zs, m_h);

  //  Check if xs and ys are out of range
  if (xs < m_xmin) {
    cerr << "MaterialPfile::sample xs out of range (min): " << xs << ", "
         << m_xmin << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (xs > m_xmax) {
    cerr << "MaterialPfile::sample xs out of range (max): " << xs << ", "
         << m_xmax << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (ys < m_ymin) {
    cerr << "MaterialPfile::sample ys out of range (min): " << ys << ", "
         << m_ymin << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (ys > m_ymax) {
    cerr << "MaterialPfile::sample ys out of range (max): " << ys << ", "
         << m_ymax << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int ii, jj;
  if (m_nstenc % 2 == 1) {
    // odd number of points in stencil
    int s = (m_nstenc - 1) / 2;
    ii = static_cast<int>(floor((xs - m_xmin) / m_h + 0.5)) - s;
    jj = static_cast<int>(floor((ys - m_ymin) / m_h + 0.5)) - s;
  } else {
    // even number of points in stencil
    int s = m_nstenc / 2;
    ii = static_cast<int>(floor((xs - m_xmin) / m_h)) - s + 1;
    jj = static_cast<int>(floor((ys - m_ymin) / m_h)) - s + 1;
  }

  // make sure we stay within array boundaries
  if (ii < 0) ii = 0;
  if (jj < 0) jj = 0;

  int ii2 = ii + m_nstenc - 1;
  int jj2 = jj + m_nstenc - 1;

  // m_x = new double[m_nx], m_nx=m_nlat
  if (ii2 >= m_nx) {
    ii2 = m_nx - 1;
    ii = ii2 - (m_nstenc - 1);
  }

  // m_y = new double[m_ny], m_ny=m_nlon
  if (jj2 >= m_ny) {
    jj2 = m_ny - 1;
    jj = jj2 - (m_nstenc - 1);
  }

  float_sw4 w = 0;
  vp = vs = rho = qp = qs = 0;
  float_sw4 appm = 0.5 * m_nstenc * m_h / sqrt(-log(1e-6));
  float_sw4 appmi2 = 1.0 / (appm * appm);

  for (int j1 = jj; j1 <= jj2; j1++)
    for (int i1 = ii; i1 <= ii2; i1++) {
      float_sw4 wgh = exp(
          -((xs - m_x[i1]) * (xs - m_x[i1]) + (ys - m_y[j1]) * (ys - m_y[j1])) *
          appmi2);
      w += wgh;

      // depth index
      int kk;
      for (kk = 1; kk < m_nmaxdepth; kk++)  // AP changed from kk <= m_nmaxdepth
      {
        if (mZ(i1 + 1, j1 + 1, kk) > zs) break;
      }
      // at this point we should have mZ(kk-1) <= zs < mZ(kk), kk < m_nmaxdepth

      int k1 = kk - 1;

      // Need to make sure that k1 is in range:
      if (k1 <= 0) k1 = 1;

      // now we should have mZ(k1) <= zs < mZ(k1+1)

      // linear interpolation factor
      float_sw4 dz = mZ(i1 + 1, j1 + 1, k1 + 1) - mZ(i1 + 1, j1 + 1, k1);
      float_sw4 factor = 0;
      //	  double factor =
      //(zs-mZ(i1+1,j1+1,k1))/(mZ(i1+1,j1+1,k1+1)-mZ(i1+1,j1+1,k1));
      //          if( factor < 0 )
      //	     factor = 0;
      if (dz != 0) {
        factor = (zs - mZ(i1 + 1, j1 + 1, k1)) / dz;
        if (factor < 0) factor = 0;
        if (factor > 1) factor = 1;
      }
      // new style
      vp += (mVp(i1 + 1, j1 + 1, k1) +
             factor * (mVp(i1 + 1, j1 + 1, k1 + 1) - mVp(i1 + 1, j1 + 1, k1))) *
            wgh;
      vs += (mVs(i1 + 1, j1 + 1, k1) +
             factor * (mVs(i1 + 1, j1 + 1, k1 + 1) - mVs(i1 + 1, j1 + 1, k1))) *
            wgh;
      rho +=
          (mRho(i1 + 1, j1 + 1, k1) +
           factor * (mRho(i1 + 1, j1 + 1, k1 + 1) - mRho(i1 + 1, j1 + 1, k1))) *
          wgh;
      // tmp
      // if (debug) printf("DEBUG: i1+1=%i, j1+1=%i, k1=%i, vp=%e, wgh=%e\n",
      // i1+1, j1+1, k1, 		    mVp(i1+1,j1+1,k1) +
      // factor*(mVp(i1+1,j1+1,k1+1)-mVp(i1+1,j1+1,k1)), wgh);

      if (m_qf) {
        // qp += (m_qp[m+k1] + factor*(m_qp[m+k1+1]-m_qp[m+k1]))*wgh;
        // qs += (m_qs[m+k1] + factor*(m_qs[m+k1+1]-m_qs[m+k1]))*wgh;
        qp += (mQp(i1 + 1, j1 + 1, k1) + factor * (mQp(i1 + 1, j1 + 1, k1 + 1) -
                                                   mQp(i1 + 1, j1 + 1, k1))) *
              wgh;
        qs += (mQs(i1 + 1, j1 + 1, k1) + factor * (mQs(i1 + 1, j1 + 1, k1 + 1) -
                                                   mQs(i1 + 1, j1 + 1, k1))) *
              wgh;
      }
    }

  // Normalize
  float_sw4 iw;
  if (w != 0.)
    iw = 1.0 / w;
  else {
    printf(
        "Error MaterialPfile::sample_cart: weight w = 0 at x=%e, y=%e, "
        "depth=%e\n",
        xs, ys, zs);
    // tmp
    printf("ii=%i, ii2=%i, jj=%i, jj2=%i, x[ii]=%e, y[jj]=%e\n", ii, ii2, jj,
           jj2, m_x[ii], m_y[jj]);
    float_sw4 dist2 =
        (xs - m_x[ii]) * (xs - m_x[ii]) + (ys - m_y[jj]) * (ys - m_y[jj]);
    float_sw4 exponent = -dist2 * appmi2;
    float_sw4 wgh = exp(exponent);
    printf("dist2=%e, appm=%e, appmi2=%e, exponent=%e, wgh=%e\n", dist2, appm,
           appmi2, exponent, wgh);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  vp *= iw;
  vs *= iw;
  rho *= iw;
  if (m_qf) {
    qp *= iw;
    qs *= iw;
  }
}

//-----------------------------------------------------------------------
MaterialPfile::~MaterialPfile() {
  // delete[] m_vp;
  // delete[] m_vs;
  // delete[] m_rho;
  // delete[] m_z;
  // if( m_qf )
  // {
  //    delete[] m_qp;
  //    delete[] m_qs;
  // }
  if (m_coords_geographic) {
    delete[] m_lat;
    delete[] m_lon;
  } else {
    delete[] m_x;
    delete[] m_y;
  }
}
