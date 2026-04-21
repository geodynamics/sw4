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

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "Byteswapper.h"
#include "EW.h"
#include "MaterialGMG.h"
#include "Require.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

using namespace std;

//-----------------------------------------------------------------------
MaterialGMG::MaterialGMG(EW* a_ew, const string a_file,
                         const string a_directory)
    : mEW(a_ew),
      m_model_file(a_file),
      m_model_dir(a_directory),
      m_use_attenuation(false),
      m_idx_rho(0),
      m_idx_vp(1),
      m_idx_vs(2),
      m_idx_qp(3),
      m_idx_qs(4) {
  mCoversAllPoints = false;
  // Check that the depths make sense
  if (a_ew != NULL) {
    m_use_attenuation = a_ew->usingAttenuation();
    read_gmg();
  }
}

//-----------------------------------------------------------------------
MaterialGMG::~MaterialGMG() {}

//-----------------------------------------------------------------------
void MaterialGMG::set_material_properties(std::vector<Sarray>& rho,
                                          std::vector<Sarray>& cs,
                                          std::vector<Sarray>& cp,
                                          std::vector<Sarray>& xis,
                                          std::vector<Sarray>& xip) {
  // Assume attenuation arrays defined on all grids if they are defined on grid
  // zero.
  bool use_q = m_use_attenuation && xis[0].is_defined() && xip[0].is_defined();
  size_t outside = 0, material = 0;

  const double yazimuthRad = m_Yaz * M_PI / 180.0;
  const double cosAz = cos(yazimuthRad);
  const double sinAz = sin(yazimuthRad);

  for (int g = 0; g < mEW->mNumberOfGrids; g++) {
    bool curvilinear =
        mEW->topographyExists() && g >= mEW->mNumberOfCartesianGrids;
    for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) {
      for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j) {
        float_sw4 x = (i - 1) * mEW->mGridSize[g];
        float_sw4 y = (j - 1) * mEW->mGridSize[g];

        int i0, j0;
        double sw4_lon, sw4_lat, gmg_x, gmg_y, gmg_x0, gmg_y0, top;
        mEW->computeGeographicCoord(x, y, sw4_lon, sw4_lat);
        /* printf("\ncomputeGeographicCoord: %f %f %f %f\n", x, y, sw4_lon,
         * sw4_lat); */

        // GMG x/y, lat/lon is switched from sw4 CRS
        mEW->computeCartesianCoordGMG(gmg_y0, gmg_x0, sw4_lon, sw4_lat, m_CRS);
        /* printf("computeCartesianCoordGMG : %f %f %f %f\n", gmg_x0, gmg_y0,
         * sw4_lon, sw4_lat); */

        const double xRel = gmg_x0 - m_Origin_x;
        const double yRel = gmg_y0 - m_Origin_y;
        gmg_x = xRel * cosAz - yRel * sinAz;
        gmg_y = xRel * sinAz + yRel * cosAz;

        const int top_i = static_cast<int>(floor(gmg_x / m_Top_hx));
        const int top_j = static_cast<int>(floor(gmg_y / m_Top_hy));
        if (top_i < 0 || top_i >= static_cast<int>(m_Top_dims[0]) ||
            top_j < 0 || top_j >= static_cast<int>(m_Top_dims[1])) {
          outside += static_cast<size_t>(mEW->m_kEnd[g] - mEW->m_kStart[g] + 1);
          continue;
        }

        top = -m_Top_surface[top_i * m_Top_dims[1] + top_j];

        for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) {
          float_sw4 z;
          if (curvilinear)
            z = mEW->mZ[g](i, j, k);
          else
            z = mEW->m_zmin[g] + (k - 1) * mEW->mGridSize[g];

          // Deal with some values on top grid that exceeds the topogrophy
          // interface
          if (g == mEW->mNumberOfGrids - 1 && z < m_Zmin) z = m_Zmin;

          // Find which block the current point belongs to
          int gr;
          for (gr = 1; gr < m_npatches; gr++)
            if (z <= top - m_ztop[gr]) break;

          gr--;
          double intf = top - m_ztop[gr];

          // When sw4 z exceeds gmg interface
          if (z < intf) z = intf;

          // i0, j0, k0 are the coordiates in GMG block
          i0 = static_cast<int>(floor(gmg_x / m_hh[gr]));
          j0 = static_cast<int>(floor(gmg_y / m_hh[gr]));
          int k0 = static_cast<int>(floor((z - intf) / m_hv[gr]));

          // (x, y, z) is the coordinate of current grid point
          if (m_Zmin <= z && z <= m_Zmax) {
            material++;

            // Extend the material value if simulation grid is larger than
            // material grid
            if (i0 >= m_ni[gr] - 1) i0 = m_ni[gr] - 2;
            if (j0 >= m_nj[gr] - 1) j0 = m_nj[gr] - 2;
            if (k0 >= m_nk[gr] - 1) k0 = m_nk[gr] - 2;
            if (k0 < 0) k0 = 0;

            // Use bilinear interpolation always:
            // Bias stencil near the boundary, need to communicate arrays
            // afterwards.
            float_sw4 wghx = (gmg_x - i0 * m_hh[gr]) / m_hh[gr];
            float_sw4 wghy = (gmg_y - j0 * m_hh[gr]) / m_hh[gr];
            float_sw4 wghz = (z - intf - k0 * m_hv[gr]) / m_hv[gr];

            /* if (x == 80000 && y == 9000) { */
            /*     printf("g=%d, ijk: %d %d %d, lalo: %f %f, converted gmg xyz:
             * %f %f %f, intf %f, gr %d, ijk %d %d %d, mat %f %f %f\n", */
            /*             g, i, j, k, sw4_lat, sw4_lon, gmg_x, gmg_y, z, intf,
             * gr, i0, j0, k0, mat(gr,0,i0,j0,k0), mat(gr,1,i0,j0,k0),
             * mat(gr,2,i0,j0,k0) ); */
            /* } */

            // weights should be within [0, 1]
            if (wghx > 1 || wghx < 0) {
#ifdef BZ_DEBUG
              printf("g=%d, sw4 (%d, %d, %d), gmg (%d, %d, %d) wghx = %.2f\n",
                     gr, i, j, k, i0, j0, k0, wghx);
#endif
              if (wghx > 1) wghx = 1;
              if (wghx < 0) wghx = 0;
            }

            if (wghy > 1 || wghy < 0) {
#ifdef BZ_DEBUG
              printf("g=%d, sw4 (%d, %d, %d), gmg (%d, %d, %d) wghy = %.2f\n",
                     gr, i, j, k, i0, j0, k0, wghy);
#endif
              if (wghy > 1) wghy = 1;
              if (wghy < 0) wghy = 0;
            }

            if (wghz > 1 || wghz < 0) {
#ifdef BZ_DEBUG
              printf("g=%d, sw4 (%d, %d, %d), gmg (%d, %d, %d) wghz = %.2f\n",
                     gr, i, j, k, i0, j0, k0, wghz);
#endif
              if (wghz > 1) wghz = 1;
              if (wghz < 0) wghz = 0;
            }

            rho[g](i, j, k) =
                (1 - wghz) *
                    ((1 - wghy) *
                         ((1 - wghx) * mat(gr, m_idx_rho, i0, j0, k0) +
                          wghx * mat(gr, m_idx_rho, i0 + 1, j0, k0)) +
                     wghy * ((1 - wghx) * mat(gr, m_idx_rho, i0, j0 + 1, k0) +
                             wghx *
                                 mat(gr, m_idx_rho, i0 + 1, j0 + 1, k0))) +
                wghz *
                    ((1 - wghy) *
                         ((1 - wghx) * mat(gr, m_idx_rho, i0, j0, k0 + 1) +
                          wghx * mat(gr, m_idx_rho, i0 + 1, j0, k0 + 1)) +
                     wghy * ((1 - wghx) *
                                 mat(gr, m_idx_rho, i0, j0 + 1, k0 + 1) +
                             wghx *
                                 mat(gr, m_idx_rho, i0 + 1, j0 + 1, k0 + 1)));

            /* if (x == 80000 && y == 9000) { */
            /*     printf("g=%d, ijk: %d %d %d, lalo: %f %f, converted gmg xyz:
             * %f %f %f, intf %f, gr %d, ijk %d %d %d, mat %f %f %f, rho=%f\n",
             */
            /*             g, i, j, k, sw4_lat, sw4_lon, gmg_x, gmg_y, z, intf,
             * gr, i0, j0, k0, mat(gr,0,i0,j0,k0), mat(gr,1,i0,j0,k0),
             * mat(gr,2,i0,j0,k0), rho[g](i, j, k)); */
            /* } */
            /* if (rho[g](i,j,k) < 1500) { */
            /*   printf("Rank %d, rho[%d](%d, %d, %d)=%.2f\n", mEW->getRank(),
             * g, i, j, k, rho[g](i,j,k)); */
            /*   ASSERT(0); */
            /* } */

            cp[g](i, j, k) =
                (1 - wghz) *
                    ((1 - wghy) *
                         ((1 - wghx) * mat(gr, m_idx_vp, i0, j0, k0) +
                          wghx * mat(gr, m_idx_vp, i0 + 1, j0, k0)) +
                     wghy * ((1 - wghx) * mat(gr, m_idx_vp, i0, j0 + 1, k0) +
                             wghx *
                                 mat(gr, m_idx_vp, i0 + 1, j0 + 1, k0))) +
                wghz *
                    ((1 - wghy) *
                         ((1 - wghx) * mat(gr, m_idx_vp, i0, j0, k0 + 1) +
                          wghx * mat(gr, m_idx_vp, i0 + 1, j0, k0 + 1)) +
                     wghy * ((1 - wghx) *
                                 mat(gr, m_idx_vp, i0, j0 + 1, k0 + 1) +
                             wghx *
                                 mat(gr, m_idx_vp, i0 + 1, j0 + 1, k0 + 1)));

            /* if (cp[g](i,j,k) < 700) { */
            /* printf("Rank %d, cp[%d](%d, %d, %d)=%.2f\n", mEW->getRank(), g,
             * i, j, k, cp[g](i,j,k)); */
            /* ASSERT(0); */
            /* } */

            cs[g](i, j, k) =
                (1 - wghz) *
                    ((1 - wghy) *
                         ((1 - wghx) * mat(gr, m_idx_vs, i0, j0, k0) +
                          wghx * mat(gr, m_idx_vs, i0 + 1, j0, k0)) +
                     wghy * ((1 - wghx) * mat(gr, m_idx_vs, i0, j0 + 1, k0) +
                             wghx *
                                 mat(gr, m_idx_vs, i0 + 1, j0 + 1, k0))) +
                wghz *
                    ((1 - wghy) *
                         ((1 - wghx) * mat(gr, m_idx_vs, i0, j0, k0 + 1) +
                          wghx * mat(gr, m_idx_vs, i0 + 1, j0, k0 + 1)) +
                     wghy * ((1 - wghx) *
                                 mat(gr, m_idx_vs, i0, j0 + 1, k0 + 1) +
                             wghx *
                                 mat(gr, m_idx_vs, i0 + 1, j0 + 1, k0 + 1)));

#ifdef BZ_DEBUG
            if (cs[g](i, j, k) < 0) {
              printf("Rank %d, cs[%d](%d, %d, %d)=%.2f\n", mEW->getRank(), g, i,
                     j, k, cs[g](i, j, k));
              ASSERT(0);
            }
#endif

            if (use_q) {
              xip[g](i, j, k) =
                  (1 - wghz) *
                      ((1 - wghy) *
                           ((1 - wghx) * mat(gr, m_idx_qp, i0, j0, k0) +
                            wghx * mat(gr, m_idx_qp, i0 + 1, j0, k0)) +
                       wghy * ((1 - wghx) *
                                   mat(gr, m_idx_qp, i0, j0 + 1, k0) +
                               wghx *
                                   mat(gr, m_idx_qp, i0 + 1, j0 + 1, k0))) +
                  wghz *
                      ((1 - wghy) *
                           ((1 - wghx) * mat(gr, m_idx_qp, i0, j0, k0 + 1) +
                            wghx * mat(gr, m_idx_qp, i0 + 1, j0, k0 + 1)) +
                       wghy * ((1 - wghx) *
                                   mat(gr, m_idx_qp, i0, j0 + 1, k0 + 1) +
                               wghx *
                                   mat(gr, m_idx_qp, i0 + 1, j0 + 1, k0 + 1)));

              xis[g](i, j, k) =
                  (1 - wghz) *
                      ((1 - wghy) *
                           ((1 - wghx) * mat(gr, m_idx_qs, i0, j0, k0) +
                            wghx * mat(gr, m_idx_qs, i0 + 1, j0, k0)) +
                       wghy * ((1 - wghx) *
                                   mat(gr, m_idx_qs, i0, j0 + 1, k0) +
                               wghx *
                                   mat(gr, m_idx_qs, i0 + 1, j0 + 1, k0))) +
                  wghz *
                      ((1 - wghy) *
                           ((1 - wghx) * mat(gr, m_idx_qs, i0, j0, k0 + 1) +
                            wghx * mat(gr, m_idx_qs, i0 + 1, j0, k0 + 1)) +
                       wghy * ((1 - wghx) *
                                   mat(gr, m_idx_qs, i0, j0 + 1, k0 + 1) +
                               wghx *
                                   mat(gr, m_idx_qs, i0 + 1, j0 + 1, k0 + 1)));
            }

          }  // End if inside
          else
            outside++;
        }  // End for i
      }    // End for j
    }      // End for k
  }        // end for g...

  free(m_CRS);
  for (int i = 0; i < m_npatches; i++) delete[] m_Material[i];
  delete[] m_Top_surface;

  mEW->communicate_arrays(rho);
  mEW->communicate_host_arrays(cs);
  mEW->communicate_host_arrays(cp);
  mEW->material_ic(rho);
  mEW->material_ic(cs);
  mEW->material_ic(cp);
  if (use_q) {
    mEW->communicate_host_arrays(xis);
    mEW->communicate_host_arrays(xip);
    mEW->material_ic(xis);
    mEW->material_ic(xip);
  }

  size_t materialSum, outsideSum;
  int mpisizelong, mpisizelonglong, mpisizeint;
  MPI_Type_size(MPI_LONG, &mpisizelong);
  MPI_Type_size(MPI_LONG_LONG, &mpisizelonglong);
  MPI_Type_size(MPI_INT, &mpisizeint);
  if (sizeof(size_t) == mpisizelong) {
    MPI_Reduce(&material, &materialSum, 1, MPI_LONG, MPI_SUM, 0,
               mEW->m_1d_communicator);
    MPI_Reduce(&outside, &outsideSum, 1, MPI_LONG, MPI_SUM, 0,
               mEW->m_1d_communicator);
  } else if (sizeof(size_t) == mpisizelonglong) {
    MPI_Reduce(&material, &materialSum, 1, MPI_LONG_LONG, MPI_SUM, 0,
               mEW->m_1d_communicator);
    MPI_Reduce(&outside, &outsideSum, 1, MPI_LONG_LONG, MPI_SUM, 0,
               mEW->m_1d_communicator);
  } else if (sizeof(size_t) == mpisizeint) {
    MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0,
               mEW->m_1d_communicator);
    MPI_Reduce(&outside, &outsideSum, 1, MPI_INT, MPI_SUM, 0,
               mEW->m_1d_communicator);
  } else {
    int materialsumi, outsidesumi, materiali = material, outsidei = outside;
    MPI_Reduce(&materiali, &materialsumi, 1, MPI_INT, MPI_SUM, 0,
               mEW->m_1d_communicator);
    MPI_Reduce(&outsidei, &outsidesumi, 1, MPI_INT, MPI_SUM, 0,
               mEW->m_1d_communicator);
    materialSum = materialsumi;
    outsideSum = outsidesumi;
  }
  if (mEW->getRank() == 0)
    //      cout << endl
    //           <<
    //           "--------------------------------------------------------------\n"
    //           << "GMG Initialized Node Types: " << endl
    //           << "   Material:        " << materialSum << endl
    //           << endl
    //           << "*Outside Domain:    " << outsideSum << endl
    //           << endl
    //           <<
    //           "--------------------------------------------------------------\n"
    //           << endl;
    cout << endl
         << "gmg command: outside = " << outsideSum
         << ", material = " << materialSum << endl;
}

#ifdef USE_HDF5
static void read_hdf5_attr(hid_t loc, hid_t dtype, const char* name,
                           void* data) {
  hid_t attr_id;
  int ierr;
  attr_id = H5Aopen(loc, name, H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  ierr = H5Aread(attr_id, dtype, data);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);
}

static char* read_hdf5_attr_str(hid_t loc, const char* name) {
  hid_t attr_id, dtype;
  int ierr;
  char* data = NULL;

  attr_id = H5Aopen(loc, name, H5P_DEFAULT);
  ASSERT(attr_id >= 0);

  dtype = H5Aget_type(attr_id);

  ierr = H5Aread(attr_id, dtype, &data);
  ASSERT(ierr >= 0);

  H5Tclose(dtype);
  H5Aclose(attr_id);

  /* fprintf(stderr, "Read data: [%s]\n", data); */
  return data;
}

static bool read_hdf5_attr_optional_f64(hid_t loc, const char* name,
                                        double& data) {
  if (H5Aexists(loc, name) <= 0) return false;
  read_hdf5_attr(loc, H5T_IEEE_F64LE, name, &data);
  return true;
}

static void read_gmg_surface_spacing(hid_t dataset_id, double& hx, double& hy) {
  double h = 0.0;
  const bool has_h = read_hdf5_attr_optional_f64(dataset_id, "resolution_horiz",
                                                 h);
  const bool has_hx =
      read_hdf5_attr_optional_f64(dataset_id, "x_resolution", hx);
  const bool has_hy =
      read_hdf5_attr_optional_f64(dataset_id, "y_resolution", hy);

  CHECK_INPUT(has_h || (has_hx && has_hy),
              "ERROR: GMG surface dataset must define resolution_horiz or "
              "both x_resolution and y_resolution");

  if (!has_hx) hx = h;
  if (!has_hy) hy = h;

  CHECK_INPUT(hx > 0 && hy > 0,
              "ERROR: GMG surface spacing must be positive, got hx="
                  << hx << " hy=" << hy);
}

static hid_t open_gmg_surface_dataset(hid_t group_id,
                                      const char** surface_name) {
  const char* candidates[] = {"top_surface", "topography_bathymetry"};
  const int ncandidates = sizeof(candidates) / sizeof(candidates[0]);

  for (int i = 0; i < ncandidates; i++) {
    if (H5Lexists(group_id, candidates[i], H5P_DEFAULT) > 0) {
      if (surface_name) *surface_name = candidates[i];
      return H5Dopen(group_id, candidates[i], H5P_DEFAULT);
    }
  }

  return -1;
}

static std::string trim_hdf5_string(const char* data, size_t len) {
  size_t end = len;
  while (end > 0 && (data[end - 1] == '\0' ||
                     isspace(static_cast<unsigned char>(data[end - 1]))))
    --end;
  return std::string(data, end);
}

static std::string normalize_gmg_component_name(const std::string& name) {
  std::string normalized;
  normalized.reserve(name.size());
  for (size_t i = 0; i < name.size(); i++) {
    unsigned char ch = static_cast<unsigned char>(name[i]);
    if (isalnum(ch)) normalized.push_back(static_cast<char>(tolower(ch)));
  }
  return normalized;
}

static std::vector<std::string> read_hdf5_attr_str_array_optional(hid_t loc,
                                                                  const char* name) {
  std::vector<std::string> values;
  if (H5Aexists(loc, name) <= 0) return values;

  hid_t attr_id = H5Aopen(loc, name, H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  hid_t dtype = H5Aget_type(attr_id);
  hid_t aspace = H5Aget_space(attr_id);
  ASSERT(dtype >= 0);
  ASSERT(aspace >= 0);

  CHECK_INPUT(H5Tget_class(dtype) == H5T_STRING,
              "ERROR: GMG attribute '" << name << "' must be a string array");

  hssize_t nvals = H5Sget_simple_extent_npoints(aspace);
  CHECK_INPUT(nvals >= 0, "ERROR: invalid size for GMG attribute '" << name << "'");
  values.reserve(static_cast<size_t>(nvals));

  if (H5Tis_variable_str(dtype)) {
    std::vector<char*> raw_values(static_cast<size_t>(nvals), NULL);
    int ierr = H5Aread(attr_id, dtype, raw_values.data());
    ASSERT(ierr >= 0);
    for (hssize_t i = 0; i < nvals; i++)
      values.push_back(raw_values[i] ? raw_values[i] : "");
#if H5_VERSION_GE(1, 12, 0)
    H5Treclaim(dtype, aspace, H5P_DEFAULT, raw_values.data());
#else
    H5Dvlen_reclaim(dtype, aspace, H5P_DEFAULT, raw_values.data());
#endif
  } else {
    const size_t elem_size = H5Tget_size(dtype);
    std::vector<char> raw_values(static_cast<size_t>(nvals) * elem_size);
    int ierr = H5Aread(attr_id, dtype, raw_values.data());
    ASSERT(ierr >= 0);
    for (hssize_t i = 0; i < nvals; i++)
      values.push_back(trim_hdf5_string(&raw_values[i * elem_size], elem_size));
  }

  H5Sclose(aspace);
  H5Tclose(dtype);
  H5Aclose(attr_id);

  return values;
}

static void resolve_gmg_component_indices(const std::vector<std::string>& components,
                                          bool use_attenuation, int& idx_rho,
                                          int& idx_vp, int& idx_vs, int& idx_qp,
                                          int& idx_qs) {
  idx_rho = 0;
  idx_vp = 1;
  idx_vs = 2;
  idx_qp = 3;
  idx_qs = 4;

  if (components.empty()) return;

  idx_rho = idx_vp = idx_vs = idx_qp = idx_qs = -1;
  for (size_t i = 0; i < components.size(); i++) {
    const std::string key = normalize_gmg_component_name(components[i]);
    if (key == "density" || key == "rho")
      idx_rho = static_cast<int>(i);
    else if (key == "vp" || key == "pvelocity")
      idx_vp = static_cast<int>(i);
    else if (key == "vs" || key == "svelocity")
      idx_vs = static_cast<int>(i);
    else if (key == "qp")
      idx_qp = static_cast<int>(i);
    else if (key == "qs")
      idx_qs = static_cast<int>(i);
  }

  CHECK_INPUT(idx_rho >= 0 && idx_vp >= 0 && idx_vs >= 0,
              "ERROR: GMG data_values must define density/rho, Vp, and Vs");
  if (use_attenuation)
    CHECK_INPUT(idx_qp >= 0 && idx_qs >= 0,
                "ERROR: GMG attenuation requires Qp and Qs in data_values");
}

struct GmgBlockInfo {
  std::string name;
  double ztop;
  double hh;
  double hv;
  hsize_t dims[4];
};

static herr_t collect_gmg_block_names(hid_t loc_id, const char* name,
                                      const H5L_info_t* /*info*/,
                                      void* operator_data) {
#if H5_VERSION_GE(1, 12, 0)
  H5O_info1_t object_info;
#else
  H5O_info_t object_info;
#endif
  std::vector<std::string>* block_names =
      static_cast<std::vector<std::string>*>(operator_data);

  ASSERT(operator_data != NULL);

#if H5_VERSION_GE(1, 12, 0)
  H5Oget_info_by_name1(loc_id, name, &object_info, H5P_DEFAULT);
#else
  H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
#endif

  if (object_info.type == H5O_TYPE_DATASET) block_names->push_back(name);

  return 0;
}
#endif

//-----------------------------------------------------------------------
void MaterialGMG::read_gmg() {
  // Timers
  double time_start, time_end;
  /* double intf_start, intf_end, mat_start, mat_end; */
  time_start = MPI_Wtime();

#ifdef USE_HDF5
  hid_t file_id, dataset_id, group_id, filespace_id, topo_grp;
  double alpha;
  herr_t ierr;
  hsize_t dims[4], top_dims[3];
  const char* surface_name = NULL;
  int str_len = 0, top_rank = 0;
  string fname = m_model_dir + "/" + m_model_file;
  std::vector<GmgBlockInfo> blocks;
  std::vector<std::string> components;
  m_npatches = 0;

  if (mEW->getRank() == 0) {
    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "Could not open hdf5 file: " << fname.c_str() << endl;
      MPI_Abort(MPI_COMM_WORLD, file_id);
    }

    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "origin_x", &m_Origin_x);
    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "origin_y", &m_Origin_y);
    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "y_azimuth", &m_Yaz);
    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "dim_z", &m_Zmax);

#ifdef BZ_DEBUG
    fprintf(stderr, "origin: %f %f, az %f, dim_z %f\n", m_Origin_x, m_Origin_y,
            m_Yaz, m_Zmax);
#endif

    m_CRS = read_hdf5_attr_str(file_id, "crs");
    str_len = (int)(strlen(m_CRS) + 1);
    components = read_hdf5_attr_str_array_optional(file_id, "data_values");
    resolve_gmg_component_indices(components, m_use_attenuation, m_idx_rho,
                                  m_idx_vp, m_idx_vs, m_idx_qp, m_idx_qs);

    group_id = H5Gopen(file_id, "blocks", H5P_DEFAULT);
    ASSERT(group_id >= 0);

    std::vector<std::string> block_names;
    ierr = H5Literate(group_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
                      collect_gmg_block_names, &block_names);
    ASSERT(ierr >= 0);

    m_npatches = (int)block_names.size();
    CHECK_INPUT(m_npatches > 0,
                "ERROR: GMG file contains no datasets in /blocks");

    blocks.resize(m_npatches);

    for (int p = 0; p < m_npatches; p++) {
      blocks[p].name = block_names[p];
      dataset_id = H5Dopen(group_id, block_names[p].c_str(), H5P_DEFAULT);
      ASSERT(dataset_id >= 0);

      filespace_id = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(filespace_id, dims, NULL);
      H5Sclose(filespace_id);

#ifdef BZ_DEBUG
      fprintf(stderr, "Rank %d, p=%d dims: %ld %ld %ld %ld\n", mEW->getRank(),
              p, dims[0], dims[1], dims[2], dims[3]);
#endif

      for (int d = 0; d < 4; d++) blocks[p].dims[d] = dims[d];

      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "z_top", &blocks[p].ztop);
      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "resolution_horiz",
                     &blocks[p].hh);
      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "resolution_vert",
                     &blocks[p].hv);

      H5Dclose(dataset_id);
    }

    std::sort(blocks.begin(), blocks.end(),
              [](const GmgBlockInfo& lhs, const GmgBlockInfo& rhs) {
                if (lhs.ztop != rhs.ztop) return lhs.ztop > rhs.ztop;
                return lhs.hv < rhs.hv;
              });

    m_hv.resize(m_npatches);
    m_hh.resize(m_npatches);
    m_ni.resize(m_npatches);
    m_nj.resize(m_npatches);
    m_nk.resize(m_npatches);
    m_nc.resize(m_npatches);
    m_ztop.resize(m_npatches);
    m_Material.resize(m_npatches);

    for (int p = 0; p < m_npatches; p++) {
      dataset_id = H5Dopen(group_id, blocks[p].name.c_str(), H5P_DEFAULT);
      ASSERT(dataset_id >= 0);

      filespace_id = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(filespace_id, dims, NULL);

      m_ni[p] = (int)blocks[p].dims[0];
      m_nj[p] = (int)blocks[p].dims[1];
      m_nk[p] = (int)blocks[p].dims[2];
      m_nc[p] = (int)blocks[p].dims[3];
      m_ztop[p] = blocks[p].ztop;
      m_hh[p] = blocks[p].hh;
      m_hv[p] = blocks[p].hv;

      int max_required_component =
          std::max(m_idx_rho, std::max(m_idx_vp, m_idx_vs));
      if (m_use_attenuation)
        max_required_component =
            std::max(max_required_component, std::max(m_idx_qp, m_idx_qs));
      CHECK_INPUT(dims[3] > static_cast<hsize_t>(max_required_component),
                  "ERROR: GMG block " << blocks[p].name
                                     << " does not contain required material "
                                        "components; nc="
                                     << dims[3]
                                     << " max required index="
                                     << max_required_component);

      m_Material[p] = new float[dims[0] * dims[1] * dims[2] * dims[3]]();
      ierr = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, filespace_id,
                     H5P_DEFAULT, &m_Material[p][0]);
      ASSERT(ierr >= 0);

      H5Sclose(filespace_id);
      H5Dclose(dataset_id);

      if (mEW->getVerbosity() >= 2) {
        printf("  GMG header block #%i (%s)\n", p, blocks[p].name.c_str());
        printf("    ztop=%f, hh=%f, hv=%f\n", m_ztop[p], m_hh[p], m_hv[p]);
        printf("    nc=%lld, ni=%lld, nj=%lld, nk=%lld\n", dims[3], dims[0],
               dims[1], dims[2]);
      }
    }  // End for each patch

    topo_grp = H5Gopen(file_id, "surfaces", H5P_DEFAULT);
    ASSERT(topo_grp >= 0);

    dataset_id = open_gmg_surface_dataset(topo_grp, &surface_name);
    CHECK_INPUT(dataset_id >= 0,
                "ERROR: GMG /surfaces must contain one of: top_surface, "
                "topography_bathymetry");

    filespace_id = H5Dget_space(dataset_id);
    top_rank = H5Sget_simple_extent_ndims(filespace_id);
    CHECK_INPUT(top_rank == 2 || top_rank == 3,
                "ERROR: GMG top_surface must be rank-2 or rank-3, got "
                    << top_rank);
    H5Sget_simple_extent_dims(filespace_id, top_dims, NULL);
    m_Top_dims[0] = top_dims[0];
    m_Top_dims[1] = top_dims[1];
    read_gmg_surface_spacing(dataset_id, m_Top_hx, m_Top_hy);
    if (top_rank == 3)
      CHECK_INPUT(top_dims[2] == 1,
                  "ERROR: GMG top_surface third dimension must be 1, got "
                      << top_dims[2]);

#ifdef BZ_DEBUG
    fprintf(stderr, "Top dims: %ld %ld\n", m_Top_dims[0], m_Top_dims[1]);
#endif

    m_Top_surface = new float[m_Top_dims[0] * m_Top_dims[1]]();
    ASSERT(m_Top_surface);

    ierr = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, filespace_id,
                   H5P_DEFAULT, m_Top_surface);
    ASSERT(ierr >= 0);

    H5Sclose(filespace_id);
    H5Dclose(dataset_id);
    H5Gclose(group_id);
    H5Gclose(topo_grp);
    H5Fclose(file_id);

    m_Zmin = 1e10;
    for (int i = 0; i < m_Top_dims[0] * m_Top_dims[1]; i++) {
      if (-m_Top_surface[i] < m_Zmin) m_Zmin = -m_Top_surface[i];
    }

  }  // End rank==0

  MPI_Bcast(&m_npatches, 1, MPI_INT, 0, mEW->m_1d_communicator);

  m_hv.resize(m_npatches);
  m_hh.resize(m_npatches);
  m_ni.resize(m_npatches);
  m_nj.resize(m_npatches);
  m_nk.resize(m_npatches);
  m_nc.resize(m_npatches);
  m_ztop.resize(m_npatches);
  m_Material.resize(m_npatches);

  MPI_Barrier(mEW->m_1d_communicator);

  MPI_Bcast(&m_Origin_x, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Origin_y, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Yaz, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Zmax, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Zmin, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Top_hx, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Top_hy, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_idx_rho, 1, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_idx_vp, 1, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_idx_vs, 1, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_idx_qp, 1, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_idx_qs, 1, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_hv[0], m_npatches, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_hh[0], m_npatches, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_ni[0], m_npatches, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_nj[0], m_npatches, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_nk[0], m_npatches, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_nc[0], m_npatches, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_ztop[0], m_npatches, MPI_DOUBLE, 0, mEW->m_1d_communicator);

  MPI_Bcast(&str_len, 1, MPI_INT, 0, mEW->m_1d_communicator);
  MPI_Bcast(m_Top_dims, 2, MPI_LONG_LONG, 0, mEW->m_1d_communicator);

  if (mEW->getRank() != 0) {
    /* fprintf(stderr, "Rank %d, strlen: %d, topo dims: %ld %ld\n",
     * mEW->getRank(), str_len, m_Top_dims[0], m_Top_dims[1]); */
    m_CRS = (char*)malloc(str_len * sizeof(char));
    m_Top_surface = new float[m_Top_dims[0] * m_Top_dims[1]];
    for (int p = 0; p < m_npatches; p++) {
      m_Material[p] = new float[m_ni[p] * m_nj[p] * m_nk[p] * m_nc[p]]();
      ASSERT(m_Material[p]);
    }
  }

  MPI_Bcast(m_CRS, str_len, MPI_CHAR, 0, mEW->m_1d_communicator);
  MPI_Bcast(m_Top_surface, m_Top_dims[0] * m_Top_dims[1], MPI_FLOAT, 0,
            mEW->m_1d_communicator);

  for (int p = 0; p < m_npatches; p++)
    MPI_Bcast(m_Material[p], m_ni[p] * m_nj[p] * m_nk[p] * m_nc[p], MPI_FLOAT,
              0, mEW->m_1d_communicator);

  CHECK_INPUT(m_Origin_x > 0 && m_Origin_y > 0,
              "ERROR: invalid GMG origin values origin_x="
                  << m_Origin_x << " origin_y=" << m_Origin_y);
  CHECK_INPUT(m_Yaz > 0, "ERROR: invalid GMG y_azimuth " << m_Yaz);
  CHECK_INPUT(m_Top_hx > 0 && m_Top_hy > 0,
              "ERROR: invalid GMG surface spacing hx="
                  << m_Top_hx << " hy=" << m_Top_hy);

  alpha = m_Yaz - 180.0;
  CHECK_INPUT(
      fabs(alpha - mEW->getGridAzimuth()) < 1e-6,
      "ERROR: gmg azimuth must be equal "
      "to coordinate system azimuth"
          << " azimuth on gmg = " << alpha
          << " azimuth of coordinate sytem = " << mEW->getGridAzimuth());

  if (mEW->getRank() == 0 && mEW->getVerbosity() >= 2) {
    printf("  GMG header: \n");
    printf("    y_azimuth=%e, origin_x=%f, origin_y=%f\n", m_Yaz, m_Origin_x,
           m_Origin_y);
    printf("    surface=%s, hx=%f, hy=%f\n",
           surface_name ? surface_name : "broadcast", m_Top_hx, m_Top_hy);
    printf("    components: rho=%d vp=%d vs=%d qp=%d qs=%d\n", m_idx_rho,
           m_idx_vp, m_idx_vs, m_idx_qp, m_idx_qs);
    printf("    nblocks=%d\n", m_npatches);
#ifdef BZ_DEBUG
    fprintf(stderr, "Rank %d, Done reading GMG data!\n", mEW->getRank());
    fprintf(stderr, "Rank %d, surface first last %f, %f\n", mEW->getRank(),
            m_Top_surface[0], m_Top_surface[m_Top_dims[0] * m_Top_dims[1] - 1]);

    for (int i = 0; i < m_npatches; i++) {
      fprintf(stderr, "Rank %d, p=%d, material first last %f, %f\n",
              mEW->getRank(), i, m_Material[i][0],
              m_Material[i][m_ni[i] * m_nj[i] * m_nk[i] * m_nc[i] - 1]);
    }
#endif
  }

  fill_in_fluids();

#ifdef BZ_DEBUG
  material_check(false);
#endif

  time_end = MPI_Wtime();
  if (mEW->getRank() == 0) {
    cout << "MaterialGMG::read_gmg, time to read material file: "
         << time_end - time_start << " seconds." << endl;
  }
#endif
}

//-----------------------------------------------------------------------
void MaterialGMG::fill_in_fluids() {
  // start from the last (bottom) block and progress upwards
  for (int p = m_npatches - 1; p >= 0; p--) {
    /* #pragma omp parallel for */
    for (int i = 0; i < m_ni[p]; i++) {
      for (int j = 0; j < m_nj[p]; j++) {
        int k0 = 0;
        while (mat(p, m_idx_vs, i, j, k0) < 0 && k0 < m_nk[p] - 1) k0++;

        // consider the case where the top block is all water. Then k0 = m_nk-1
        // and mat(Vs) <0
        if (k0 == m_nk[p] - 1 && mat(p, m_idx_vs, i, j, k0) < 0) {
          // get value from block p+1
          if (p < m_npatches - 1) {
            int pd = p + 1, id, jd, kd;  // index of donor block
            float_sw4 xm = i * m_hh[p];
            float_sw4 ym = j * m_hh[p];
            // get closest (id,jd) index on patch pd
            id = static_cast<int>(xm / m_hh[pd]);
            jd = static_cast<int>(ym / m_hh[pd]);
            kd = 0;  // get value from top of block pd

            if (!(id >= 0 && id < m_ni[pd] && jd >= 0 && jd < m_nj[pd])) {
              // out of bounds: find nearest interior point
              if (id > m_ni[pd] - 1) id = m_ni[pd] - 1;
              if (jd > m_nj[pd] - 1) jd = m_nj[pd] - 1;

              printf(
                  "WARNING: nearest grid point to (%e,%e) was outside local "
                  "part of block pd=%i\n"
                  " using id=%i, jd=%i, at (%e, %e)\n",
                  xm, ym, pd, id, jd, (id - 1) * m_hh[pd], (jd - 1) * m_hh[pd]);
            }

            // debug
            /* fprintf(stderr, "p=%d, ijk: %d %d %d, go to next block for valid
             * value, rho=%f\n", p, i, j, k0, mat(pd,0,id,jd,kd)); */

            // get values from block 'pd'
            mat_assign(p, m_idx_rho, i, j, k0,
                       mat(pd, m_idx_rho, id, jd, kd));
            mat_assign(p, m_idx_vp, i, j, k0, mat(pd, m_idx_vp, id, jd, kd));
            mat_assign(p, m_idx_vs, i, j, k0, mat(pd, m_idx_vs, id, jd, kd));
            if (m_use_attenuation) {
              mat_assign(p, m_idx_qp, i, j, k0,
                         mat(pd, m_idx_qp, id, jd, kd));
              mat_assign(p, m_idx_qs, i, j, k0,
                         mat(pd, m_idx_qs, id, jd, kd));
            }
          } else {
            printf(
                "ERROR: found undefined material properties in last material "
                "block\n"
                " patch p=%i, i=%i, j=%i, k0=%i\n",
                p, i, j, k0);
          }
        }

        // debug
        /* if (k0 > 0) { */
        /*    /1* fprintf(stderr, "p=%d, ijk: %d %d %d, s=%f\n", p, i, j, k0,
         * mat(p,2,i,j,k0-1)); *1/ */
        /*   fprintf(stderr, "p=%d, ijk: %d %d 0 to %d, assign rho=%f\n", p, i,
         * j, k0, mat(p,0,i,j,k0)); */
        /* } */

        for (int k = 0; k < k0; k++) {
          mat_assign(p, m_idx_rho, i, j, k, mat(p, m_idx_rho, i, j, k0));
          mat_assign(p, m_idx_vp, i, j, k, mat(p, m_idx_vp, i, j, k0));
          mat_assign(p, m_idx_vs, i, j, k, mat(p, m_idx_vs, i, j, k0));
          if (m_use_attenuation) {
            mat_assign(p, m_idx_qp, i, j, k, mat(p, m_idx_qp, i, j, k0));
            mat_assign(p, m_idx_qs, i, j, k, mat(p, m_idx_qs, i, j, k0));
          }
        }  // End for k

      }  // End for i
    }    // End for j
  }      // End for p
}

//-----------------------------------------------------------------------
void MaterialGMG::material_check(bool water) {
  bool printsmallcpcs = false;
  for (int p = 0; p < m_npatches; p++) {
    double csmin = 1e38, cpmin = 1e38, cratmin = 1e38, csmax = -1e38,
           cpmax = -1e38, cratmax = -1e38;
    double rhomin = 1e38, rhomax = -1e38;
    for (int i = 0; i < m_ni[p]; i++)
      for (int j = 0; j < m_nj[p]; j++)
        for (int k = 0; k < m_nk[p]; k++) {
          if (water || mat(p, m_idx_vs, i, j, k) > 0) {
            if (mat(p, m_idx_rho, i, j, k) < rhomin)
              rhomin = mat(p, m_idx_rho, i, j, k);
            if (mat(p, m_idx_rho, i, j, k) > rhomax)
              rhomax = mat(p, m_idx_rho, i, j, k);
            if (mat(p, m_idx_vp, i, j, k) < cpmin)
              cpmin = mat(p, m_idx_vp, i, j, k);
            if (mat(p, m_idx_vp, i, j, k) > cpmax)
              cpmax = mat(p, m_idx_vp, i, j, k);
            if (mat(p, m_idx_vs, i, j, k) < csmin)
              csmin = mat(p, m_idx_vs, i, j, k);
            if (mat(p, m_idx_vs, i, j, k) > csmax)
              csmax = mat(p, m_idx_vs, i, j, k);
            double crat =
                mat(p, m_idx_vp, i, j, k) / mat(p, m_idx_vs, i, j, k);
            if (crat < cratmin) {
              cratmin = crat;
              if (printsmallcpcs && crat < 1.41) {
                cout << "crat= " << crat << " at " << i << " " << j << " " << k
                     << endl;
                cout << " material is " << mat(p, m_idx_rho, i, j, k) << " "
                     << mat(p, m_idx_vp, i, j, k) << " "
                     << mat(p, m_idx_vs, i, j, k) << " ";
                if (m_idx_qp >= 0)
                  cout << mat(p, m_idx_qp, i, j, k) << " ";
                if (m_idx_qs >= 0)
                  cout << mat(p, m_idx_qs, i, j, k);
                cout << endl;
              }
            }
            if (crat > cratmax) cratmax = crat;
          }
        }
    double cmins[4] = {csmin, cpmin, cratmin, rhomin},
           cmaxs[4] = {csmax, cpmax, cratmax, rhomax};
    double cminstot[4], cmaxstot[4];
    MPI_Reduce(cmins, cminstot, 4, MPI_DOUBLE, MPI_MIN, 0,
               mEW->m_1d_communicator);
    MPI_Reduce(cmaxs, cmaxstot, 4, MPI_DOUBLE, MPI_MAX, 0,
               mEW->m_1d_communicator);
    int myid;
    MPI_Comm_rank(mEW->m_1d_communicator, &myid);
    if (myid == 0)
    //	 if( mEW->getRank()==0 )
    {
      if (p == 1 && !water)
        cout << "GMG-file limits, away from water: " << endl;
      else if (p == 1)
        cout << "GMG-file limits : " << endl;
      cout << "  Patch no " << p << " : " << endl;
      cout << "    cp    min and max " << cminstot[1] << " " << cmaxstot[1]
           << endl;
      cout << "    cs    min and max " << cminstot[0] << " " << cmaxstot[0]
           << endl;
      cout << "    cp/cs min and max " << cminstot[2] << " " << cmaxstot[2]
           << endl;
      cout << "    rho   min and max " << cminstot[3] << " " << cmaxstot[3]
           << endl;
    }
  }
}
