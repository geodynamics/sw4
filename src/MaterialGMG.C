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
      m_use_attenuation(false) {
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

        i0 = static_cast<int>(floor(gmg_x / m_hh[0]));
        j0 = static_cast<int>(floor(gmg_y / m_hh[0]));

        top = -m_Top_surface[i0 * m_Top_dims[1] + j0];

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
                    ((1 - wghy) * ((1 - wghx) * mat(gr, 0, i0, j0, k0) +
                                   wghx * mat(gr, 0, i0 + 1, j0, k0)) +
                     wghy * ((1 - wghx) * mat(gr, 0, i0, j0 + 1, k0) +
                             wghx * mat(gr, 0, i0 + 1, j0 + 1, k0))) +
                wghz * ((1 - wghy) * ((1 - wghx) * mat(gr, 0, i0, j0, k0 + 1) +
                                      wghx * mat(gr, 0, i0 + 1, j0, k0 + 1)) +
                        wghy * ((1 - wghx) * mat(gr, 0, i0, j0 + 1, k0 + 1) +
                                wghx * mat(gr, 0, i0 + 1, j0 + 1, k0 + 1)));

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
                    ((1 - wghy) * ((1 - wghx) * mat(gr, 1, i0, j0, k0) +
                                   wghx * mat(gr, 1, i0 + 1, j0, k0)) +
                     wghy * ((1 - wghx) * mat(gr, 1, i0, j0 + 1, k0) +
                             wghx * mat(gr, 1, i0 + 1, j0 + 1, k0))) +
                wghz * ((1 - wghy) * ((1 - wghx) * mat(gr, 1, i0, j0, k0 + 1) +
                                      wghx * mat(gr, 1, i0 + 1, j0, k0 + 1)) +
                        wghy * ((1 - wghx) * mat(gr, 1, i0, j0 + 1, k0 + 1) +
                                wghx * mat(gr, 1, i0 + 1, j0 + 1, k0 + 1)));

            /* if (cp[g](i,j,k) < 700) { */
            /* printf("Rank %d, cp[%d](%d, %d, %d)=%.2f\n", mEW->getRank(), g,
             * i, j, k, cp[g](i,j,k)); */
            /* ASSERT(0); */
            /* } */

            cs[g](i, j, k) =
                (1 - wghz) *
                    ((1 - wghy) * ((1 - wghx) * mat(gr, 2, i0, j0, k0) +
                                   wghx * mat(gr, 2, i0 + 1, j0, k0)) +
                     wghy * ((1 - wghx) * mat(gr, 2, i0, j0 + 1, k0) +
                             wghx * mat(gr, 2, i0 + 1, j0 + 1, k0))) +
                wghz * ((1 - wghy) * ((1 - wghx) * mat(gr, 2, i0, j0, k0 + 1) +
                                      wghx * mat(gr, 2, i0 + 1, j0, k0 + 1)) +
                        wghy * ((1 - wghx) * mat(gr, 2, i0, j0 + 1, k0 + 1) +
                                wghx * mat(gr, 2, i0 + 1, j0 + 1, k0 + 1)));

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
                      ((1 - wghy) * ((1 - wghx) * mat(gr, 3, i0, j0, k0) +
                                     wghx * mat(gr, 3, i0 + 1, j0, k0)) +
                       wghy * ((1 - wghx) * mat(gr, 3, i0, j0 + 1, k0) +
                               wghx * mat(gr, 3, i0 + 1, j0 + 1, k0))) +
                  wghz *
                      ((1 - wghy) * ((1 - wghx) * mat(gr, 3, i0, j0, k0 + 1) +
                                     wghx * mat(gr, 3, i0 + 1, j0, k0 + 1)) +
                       wghy * ((1 - wghx) * mat(gr, 3, i0, j0 + 1, k0 + 1) +
                               wghx * mat(gr, 3, i0 + 1, j0 + 1, k0 + 1)));

              xis[g](i, j, k) =
                  (1 - wghz) *
                      ((1 - wghy) * ((1 - wghx) * mat(gr, 4, i0, j0, k0) +
                                     wghx * mat(gr, 4, i0 + 1, j0, k0)) +
                       wghy * ((1 - wghx) * mat(gr, 4, i0, j0 + 1, k0) +
                               wghx * mat(gr, 4, i0 + 1, j0 + 1, k0))) +
                  wghz *
                      ((1 - wghy) * ((1 - wghx) * mat(gr, 4, i0, j0, k0 + 1) +
                                     wghx * mat(gr, 4, i0 + 1, j0, k0 + 1)) +
                       wghy * ((1 - wghx) * mat(gr, 4, i0, j0 + 1, k0 + 1) +
                               wghx * mat(gr, 4, i0 + 1, j0 + 1, k0 + 1)));
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
  hsize_t dims[4];
  char grid_name[128];
  int str_len, hv[5];
  string fname = m_model_dir + "/" + m_model_file;

  // Fixed for GMG grids
  m_npatches = 4;

  m_hv.resize(m_npatches);
  m_hh.resize(m_npatches);
  m_ni.resize(m_npatches);
  m_nj.resize(m_npatches);
  m_nk.resize(m_npatches);
  m_nc.resize(m_npatches);
  m_ztop.resize(m_npatches);
  m_Material.resize(m_npatches);

  hv[0] = 25;
  hv[1] = 50;
  hv[2] = 125;
  hv[3] = 250;

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

    // Origin_x is not correctly read sometimes
    if (m_Origin_x < 1) {
      m_Origin_x = 99286.2;
      if (mEW->getRank() == 0)
        printf("GMG origin_x read invalid value, correct to 99286.2 \n");
    }

    m_CRS = read_hdf5_attr_str(file_id, "crs");
    str_len = (int)(strlen(m_CRS) + 1);

    group_id = H5Gopen(file_id, "blocks", H5P_DEFAULT);
    ASSERT(group_id >= 0);

    for (int p = 0; p < m_npatches; p++) {
      sprintf(grid_name, "vres%dm", hv[p]);
      dataset_id = H5Dopen(group_id, grid_name, H5P_DEFAULT);
      ASSERT(dataset_id >= 0);

      filespace_id = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(filespace_id, dims, NULL);

#ifdef BZ_DEBUG
      fprintf(stderr, "Rank %d, p=%d dims: %ld %ld %ld %ld\n", mEW->getRank(),
              p, dims[0], dims[1], dims[2], dims[3]);
#endif

      m_ni[p] = (int)dims[0];
      m_nj[p] = (int)dims[1];
      m_nk[p] = (int)dims[2];
      m_nc[p] = (int)dims[3];

      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "z_top", &m_ztop[p]);
      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "resolution_horiz", &m_hh[p]);
      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "resolution_vert", &m_hv[p]);

      // Make assumption is correct with the data
      ASSERT(hv[p] == (int)m_hv[p]);
      ASSERT(dims[3] == 7);

      m_Material[p] = new float[dims[0] * dims[1] * dims[2] * dims[3]]();
      ierr = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, filespace_id,
                     H5P_DEFAULT, &m_Material[p][0]);
      ASSERT(ierr >= 0);

      H5Sclose(filespace_id);
      H5Dclose(dataset_id);

      if (mEW->getVerbosity() >= 2) {
        printf("  GMG header block #%i\n", p);
        printf("    hh=%f, hv=%f\n", m_hh[p], m_hv[p]);
        printf("    nc=%lld, ni=%lld, nj=%lld, nk=%lld\n", dims[3], dims[0],
               dims[1], dims[2]);
      }
    }  // End for each patch

    topo_grp = H5Gopen(file_id, "surfaces", H5P_DEFAULT);
    ASSERT(topo_grp >= 0);

    dataset_id = H5Dopen(topo_grp, "top_surface", H5P_DEFAULT);
    ASSERT(dataset_id >= 0);

    filespace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(filespace_id, &m_Top_dims[0], NULL);

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

  MPI_Barrier(mEW->m_1d_communicator);

  MPI_Bcast(&m_Origin_x, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Origin_y, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Yaz, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Zmax, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
  MPI_Bcast(&m_Zmin, 1, MPI_DOUBLE, 0, mEW->m_1d_communicator);
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
    m_CRS = new char[str_len]();
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

  ASSERT(m_Origin_x > 0);
  ASSERT(m_Origin_y > 0);
  ASSERT(m_Yaz > 0);

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
        // Vs is 2 in GMG model
        while (mat(p, 2, i, j, k0) < 0 && k0 < m_nk[p] - 1) k0++;

        // consider the case where the top block is all water. Then k0 = m_nk-1
        // and mat(Vs) <0
        if (k0 == m_nk[p] - 1 && mat(p, 2, i, j, k0) < 0) {
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
            mat_assign(p, 0, i, j, k0, mat(pd, 0, id, jd, kd));
            mat_assign(p, 1, i, j, k0, mat(pd, 1, id, jd, kd));
            mat_assign(p, 2, i, j, k0, mat(pd, 2, id, jd, kd));
            if (m_use_attenuation) {
              mat_assign(p, 3, i, j, k0, mat(pd, 3, id, jd, kd));
              mat_assign(p, 4, i, j, k0, mat(pd, 4, id, jd, kd));
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
          mat_assign(p, 0, i, j, k, mat(p, 0, i, j, k0));
          mat_assign(p, 1, i, j, k, mat(p, 1, i, j, k0));
          mat_assign(p, 2, i, j, k, mat(p, 2, i, j, k0));
          if (m_use_attenuation) {
            mat_assign(p, 3, i, j, k, mat(p, 3, i, j, k0));
            mat_assign(p, 4, i, j, k, mat(p, 4, i, j, k0));
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
          if (water || mat(p, 2, i, j, k) > 0) {
            if (mat(p, 0, i, j, k) < rhomin) rhomin = mat(p, 0, i, j, k);
            if (mat(p, 0, i, j, k) > rhomax) rhomax = mat(p, 0, i, j, k);
            if (mat(p, 1, i, j, k) < cpmin) cpmin = mat(p, 1, i, j, k);
            if (mat(p, 1, i, j, k) > cpmax) cpmax = mat(p, 1, i, j, k);
            if (mat(p, 2, i, j, k) < csmin) csmin = mat(p, 2, i, j, k);
            if (mat(p, 2, i, j, k) > csmax) csmax = mat(p, 2, i, j, k);
            double crat = mat(p, 1, i, j, k) / mat(p, 2, i, j, k);
            if (crat < cratmin) {
              cratmin = crat;
              if (printsmallcpcs && crat < 1.41) {
                cout << "crat= " << crat << " at " << i << " " << j << " " << k
                     << endl;
                cout << " material is " << mat(p, 0, i, j, k) << " "
                     << mat(p, 1, i, j, k) << " " << mat(p, 2, i, j, k) << " "
                     << mat(p, 3, i, j, k) << " " << mat(p, 4, i, j, k) << endl;
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
