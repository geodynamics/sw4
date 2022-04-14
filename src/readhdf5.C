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

#ifndef READHDF5_C
#define READHDF5_C

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <ctime>
#include <iostream>
#include <sstream>

#include "EW.h"
#include "Require.h"
#include "TimeSeries.h"
#include "readhdf5.h"

#ifdef USE_HDF5
#include "hdf5.h"

struct traverse_data_t {
  int myRank;
  EW *ew;
  /* bool cartCoordSet; */
  string inFileName;
  string outFileName;
  int writeEvery;
  int downSample;
  TimeSeries::receiverMode mode;
  int event;
  vector<vector<TimeSeries *> > *GlobalTimeSeries;
  float_sw4 m_global_xmax;
  float_sw4 m_global_ymax;
  bool is_obs;
  bool winlset;
  bool winrset;
  float_sw4 winl;
  float_sw4 winr;
  bool usex;
  bool usey;
  bool usez;
  float_sw4 t0;
  bool scalefactor_set;
  float_sw4 scalefactor;
} traverse_data_t;

struct traverse_data2_t {
  vector<string> *staname;
  vector<double> *x;
  vector<double> *y;
  vector<double> *z;
  vector<int> *is_nsew;
  int *n;
} traverse_data2_t;

struct srf_meta_t {
  float elon;
  float elat;
  int nstk;
  int ndip;
  float len;
  float wid;
  float stk;
  float dip;
  float dtop;
  float shyp;
  float dhyp;
} srf_meta_t;

struct srf_data_t {
  float lon;
  float lat;
  float dep;
  float stk;
  float dip;
  float area;
  float tinit;
  float dt;
  float vs;
  float den;
  float rake;
  float slip1;
  int nt1;
  float slip2;
  int nt2;
  float slip3;
  int nt3;
} srf_data_t;

static herr_t traverse_func(hid_t loc_id, const char *grp_name,
                            const H5L_info_t *info, void *operator_data) {
  hid_t grp, dset, attr;
#if H5_VERSION_GE(1, 12, 0)
  H5O_info1_t infobuf;
#else
  H5O_info_t infobuf;
#endif
  EW *a_ew;
  double data[3];
  double lon, lat, depth, x, y, z;
  bool geoCoordSet = true, topodepth = false, nsew = true;
  int isnsew, ret;

  ASSERT(operator_data != NULL);

  struct traverse_data_t *op_data = (struct traverse_data_t *)operator_data;
  a_ew = op_data->ew;
  ASSERT(a_ew != NULL);

#if H5_VERSION_GE(1, 12, 0)
  H5Oget_info_by_name1(loc_id, grp_name, &infobuf, H5P_DEFAULT);
#else
  H5Oget_info_by_name(loc_id, grp_name, &infobuf, H5P_DEFAULT);
#endif
  if (infobuf.type == H5O_TYPE_GROUP) {
    /* if (op_data->myRank == 0) */
    /*   printf ("Group: [%s] \n", grp_name); */

    // read x,y,z or ns,ew,up
    grp = H5Gopen(loc_id, grp_name, H5P_DEFAULT);
    if (grp < 0) {
      fprintf(stderr, "Error opening group [%s]\n", grp_name);
      return -1;
    }

    if (H5Lexists(grp, "ISNSEW", H5P_DEFAULT) > 0) {
      attr = H5Dopen(grp, "ISNSEW", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret =
          H5Dread(attr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &isnsew);
      ASSERT(ret >= 0);
      H5Dclose(attr);
      if (isnsew == 0) {
        nsew = false;
        geoCoordSet = false;
      }
    }

    if (H5Lexists(grp, "WindowL", H5P_DEFAULT) > 0) {
      op_data->winlset = true;
      attr = H5Dopen(grp, "WindowL", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret = H5Dread(attr, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &op_data->winl);
      ASSERT(ret >= 0);
      H5Dclose(attr);
    }

    if (H5Lexists(grp, "WindowR", H5P_DEFAULT) > 0) {
      op_data->winrset = true;
      attr = H5Dopen(grp, "WindowR", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret = H5Dread(attr, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &op_data->winr);
      ASSERT(ret >= 0);
      H5Dclose(attr);
    }

    if (H5Lexists(grp, "STX,STY,STZ", H5P_DEFAULT) > 0) {
      // X, Y, Z
      dset = H5Dopen(grp, "STX,STY,STZ", H5P_DEFAULT);
      if (dset < 0)
        fprintf(
            stderr,
            "Error reading from rechdf5 station %s, STX,STY,STZ open failed!\n",
            grp_name);
      ASSERT(attr > 0);
      ret =
          H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      ASSERT(ret >= 0);
      H5Dclose(dset);
      x = data[0];
      y = data[1];
      z = data[2];
      topodepth = true;
    } else if (H5Lexists(grp, "STLA,STLO,STDP", H5P_DEFAULT) > 0) {
      // STLA,STLO,STDP
      dset = H5Dopen(grp, "STLA,STLO,STDP", H5P_DEFAULT);
      if (dset < 0)
        fprintf(stderr,
                "Error reading from rechdf5 station %s, STLA,STLO,STDP open "
                "failed!\n",
                grp_name);
      ASSERT(attr > 0);
      ret =
          H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      ASSERT(ret >= 0);
      H5Dclose(dset);
      lat = data[0];
      lon = data[1];
      z = data[2];
      topodepth = true;
    } else {
      // Not a station group, ignore
      H5Gclose(grp);
      return 0;
    }

    depth = z;
    if (geoCoordSet) a_ew->computeCartesianCoord(x, y, lon, lat);

    bool inCurvilinear = false;
    // we are in or above the curvilinear grid
    if (a_ew->topographyExists() &&
        z < a_ew->m_zmin[a_ew->mNumberOfCartesianGrids - 1])
      inCurvilinear = true;

    // check if (x,y,z) is not in the global bounding box
    if (!((inCurvilinear || z >= 0) && x >= 0 && x <= op_data->m_global_xmax &&
          y >= 0 && y <= op_data->m_global_ymax)) {
      // The location of this station was outside the domain, so don't include
      // it in the global list
      if (op_data->myRank == 0 && a_ew->getVerbosity() > 0) {
        stringstream receivererr;

        receivererr << endl
                    << "***************************************************"
                    << endl
                    << " WARNING:  RECEIVER positioned outside grid!" << endl;
        receivererr << " No RECEIVER file will be generated for file = "
                    << op_data->outFileName << endl;
        if (geoCoordSet) {
          receivererr << " @ lon=" << lon << " lat=" << lat
                      << " depth=" << depth << endl
                      << endl;
        } else {
          receivererr << " @ x=" << x << " y=" << y << " z=" << z << endl
                      << endl;
        }

        receivererr << "***************************************************"
                    << endl;
        cerr << receivererr.str();
        cerr.flush();
      }
    } else {
      /* if (op_data->myRank == 0) */
      /*   cout << "x=" << x << ", y=" << y << ", z=" << z << ", writeEvery=" <<
       * op_data->writeEvery << endl; */

      TimeSeries *ts_ptr = new TimeSeries(
          a_ew, grp_name, grp_name, op_data->mode, false, false, true,
          op_data->outFileName, x, y, z, topodepth, op_data->writeEvery,
          op_data->downSample, !nsew, op_data->event);

      if ((*op_data->GlobalTimeSeries)[op_data->event].size() == 0) {
        ts_ptr->allocFid();
        ts_ptr->setTS0Ptr(ts_ptr);
      } else {
        ts_ptr->setFidPtr(
            (*op_data->GlobalTimeSeries)[op_data->event][0]->getFidPtr());
        ts_ptr->setTS0Ptr((*op_data->GlobalTimeSeries)[op_data->event][0]);
      }

      if (ts_ptr->myPoint()) {
        /* cout << "Rank " << op_data->myRank << "has this point x=" << x << "
         * y=" << y << " z=" << z << endl; */

        // Only for observation data
        if (op_data->is_obs) {
          // Read data
          bool ignore_utc = false;
          ts_ptr->readSACHDF5(op_data->ew, op_data->inFileName, ignore_utc);

          // Set reference UTC to simulation UTC, for easier plotting.
          ts_ptr->set_utc_to_simulation_utc();

          // Set window, in simulation time
          if (op_data->winlset || op_data->winrset) {
            if (op_data->winlset && !op_data->winrset) op_data->winr = 1e38;
            if (!op_data->winlset && op_data->winrset) op_data->winl = -1;
            ts_ptr->set_window(op_data->winl, op_data->winr);
          }

          // Exclude some components
          if (!op_data->usex || !op_data->usey || !op_data->usez)
            ts_ptr->exclude_component(op_data->usex, op_data->usey,
                                      op_data->usez);

          // Add extra shift from command line, use with care.
          if (op_data->t0 != 0) ts_ptr->add_shift(op_data->t0);

          // Set scale factor if given
          if (op_data->scalefactor_set)
            ts_ptr->set_scalefactor(op_data->scalefactor);
        }
      }

      // include the receiver in the global list
      (*op_data->GlobalTimeSeries)[op_data->event].push_back(ts_ptr);
    }

    H5Gclose(grp);
  }

  return 0;
}

void readStationHDF5(EW *ew, string inFileName, string outFileName,
                     int writeEvery, int downSample,
                     TimeSeries::receiverMode mode, int event,
                     vector<vector<TimeSeries *> > *GlobalTimeSeries,
                     float_sw4 m_global_xmax, float_sw4 m_global_ymax,
                     bool is_obs, bool winlset, bool winrset, float_sw4 winl,
                     float_sw4 winr, bool usex, bool usey, bool usez,
                     float_sw4 t0, bool scalefactor_set,
                     float_sw4 scalefactor) {
  hid_t fid, fapl;

  struct traverse_data_t tData;
  tData.myRank = ew->getRank();
  tData.ew = ew;
  tData.inFileName = inFileName;
  tData.outFileName = outFileName;
  tData.writeEvery = writeEvery;
  tData.downSample = downSample;
  tData.mode = mode;
  tData.event = event;
  tData.m_global_xmax = m_global_xmax;
  tData.m_global_ymax = m_global_ymax;
  tData.GlobalTimeSeries = GlobalTimeSeries;
  tData.is_obs = is_obs;
  tData.winlset = winlset;
  tData.winrset = winrset;
  tData.winl = winl;
  tData.winr = winr;
  tData.usex = usex;
  tData.usey = usey;
  tData.usez = usez;
  tData.t0 = t0;
  tData.scalefactor_set = scalefactor_set;
  tData.scalefactor = scalefactor;

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);

  fid = H5Fopen(inFileName.c_str(), H5F_ACC_RDONLY, fapl);
  if (fid < 0) {
    printf("%s Error opening file [%s]\n", __func__, inFileName.c_str());
    return;
  }

  if (inFileName == outFileName) {
    if (tData.myRank == 0) {
      printf("Warning: Same station input file and output file name [%s]\n",
             inFileName.c_str());
    }
  }

  H5Literate(fid, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, traverse_func, &tData);

  H5Pclose(fapl);
  H5Fclose(fid);
}

static herr_t traverse_func2(hid_t loc_id, const char *grp_name,
                             const H5L_info_t *info, void *operator_data) {
  hid_t grp, dset, attr;
#if H5_VERSION_GE(1, 12, 0)
  H5O_info1_t infobuf;
#else
  H5O_info_t infobuf;
#endif
  float data[3];
  int isnsew, ret;

  ASSERT(operator_data != NULL);

  struct traverse_data2_t *op_data = (struct traverse_data2_t *)operator_data;

#if H5_VERSION_GE(1, 12, 0)
  H5Oget_info_by_name1(loc_id, grp_name, &infobuf, H5P_DEFAULT);
#else
  H5Oget_info_by_name(loc_id, grp_name, &infobuf, H5P_DEFAULT);
#endif
  if (infobuf.type == H5O_TYPE_GROUP) {
    /* if (op_data->myRank == 0) */
    /*   printf ("Group: [%s] \n", grp_name); */

    // read x,y,z or ns,ew,up
    grp = H5Gopen(loc_id, grp_name, H5P_DEFAULT);
    if (grp < 0) {
      printf("Error opening group [%s]\n", grp_name);
      return -1;
    }

    if (H5Lexists(grp, "ISNSEW", H5P_DEFAULT) > 0) {
      attr = H5Dopen(grp, "ISNSEW", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret =
          H5Dread(attr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &isnsew);
      ASSERT(ret >= 0);
      H5Dclose(attr);
    }

    if (H5Lexists(grp, "STX,STY,STZ", H5P_DEFAULT) > 0) {
      // X, Y, Z
      dset = H5Dopen(grp, "STX,STY,STZ", H5P_DEFAULT);
      if (dset < 0)
        fprintf(
            stderr,
            "Error reading from rechdf5 station %s, STX,STY,STZ open failed!\n",
            grp_name);
      ASSERT(attr > 0);
      ret =
          H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      ASSERT(ret >= 0);
      H5Dclose(dset);
    } else if (H5Lexists(grp, "STLA,STLO,STDP", H5P_DEFAULT) > 0) {
      // STLA,STLO,STDP
      dset = H5Dopen(grp, "STLA,STLO,STDP", H5P_DEFAULT);
      if (dset < 0)
        fprintf(stderr,
                "Error reading from rechdf5 station %s, STLA,STLO,STDP open "
                "failed!\n",
                grp_name);
      ASSERT(attr > 0);
      ret =
          H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      ASSERT(ret >= 0);
      H5Dclose(dset);
    } else {
      // Not a station group, ignore
      H5Gclose(grp);
      return 0;
    }

    string staname = grp_name;
    (*op_data->staname).push_back(staname);
    (*op_data->x).push_back(data[0]);
    (*op_data->y).push_back(data[1]);
    (*op_data->z).push_back(data[2]);
    (*op_data->is_nsew).push_back(isnsew);
    (*op_data->n)++;

    H5Gclose(grp);
  }

  return 0;
}

void readStationInfoHDF5(string inFileName, vector<string> *staname,
                         vector<double> *x, vector<double> *y,
                         vector<double> *z, vector<int> *is_nsew, int *n) {
  hid_t fid, fapl;

  struct traverse_data2_t tData;
  tData.staname = staname;
  tData.x = x;
  tData.y = y;
  tData.z = z;
  tData.is_nsew = is_nsew;
  tData.n = n;

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL);

  fid = H5Fopen(inFileName.c_str(), H5F_ACC_RDONLY, fapl);
  if (fid < 0) {
    printf("%s Error opening file [%s]\n", __func__, inFileName.c_str());
    return;
  }

  H5Literate(fid, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, traverse_func2, &tData);

  H5Pclose(fapl);
  H5Fclose(fid);
}

void readRuptureHDF5(char *fname,
                     vector<vector<Source *> > &a_GlobalUniqueSources, EW *ew,
                     int event, float_sw4 m_global_xmax,
                     float_sw4 m_global_ymax, float_sw4 m_global_zmax,
                     float_sw4 mGeoAz, float_sw4 xmin, float_sw4 ymin,
                     float_sw4 zmin, int mVerbose, int nreader) {
  bool is_debug = true;
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (nreader <= 0) nreader = 1;

  int read_color = world_rank % (world_size / nreader) == 0 ? 0 : 1;
  int node_color = world_rank / (world_size / nreader);
  int read_rank, read_size;
  MPI_Comm read_comm, node_comm;
  MPI_Comm_split(MPI_COMM_WORLD, read_color, world_rank, &read_comm);
  MPI_Comm_split(MPI_COMM_WORLD, node_color, world_rank, &node_comm);
  MPI_Comm_rank(MPI_COMM_WORLD, &read_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &read_size);

  double stime, etime;
  hid_t fid, attr, ctype, dtype, dset, dspace, aspace, fapl;

  ctype = H5Tcreate(H5T_COMPOUND, 9 * sizeof(float) + 2 * sizeof(int));
  H5Tinsert(ctype, "ELON", 0, H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "ELAT", 1 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "NSTK", 2 * sizeof(float), H5T_NATIVE_INT);
  H5Tinsert(ctype, "NDIP", 2 * sizeof(float) + 1 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert(ctype, "LEN", 2 * sizeof(float) + 2 * sizeof(int),
            H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "WID", 3 * sizeof(float) + 2 * sizeof(int),
            H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "STK", 4 * sizeof(float) + 2 * sizeof(int),
            H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "DIP", 5 * sizeof(float) + 2 * sizeof(int),
            H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "DTOP", 6 * sizeof(float) + 2 * sizeof(int),
            H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "SHYP", 7 * sizeof(float) + 2 * sizeof(int),
            H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "DHYP", 8 * sizeof(float) + 2 * sizeof(int),
            H5T_NATIVE_FLOAT);

  dtype = H5Tcreate(H5T_COMPOUND, 14 * sizeof(float) + 3 * sizeof(int));
  H5Tinsert(dtype, "LON", 0, H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "LAT", 1 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "DEP", 2 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "STK", 3 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "DIP", 4 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "AREA", 5 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "TINIT", 6 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "DT", 7 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "VS", 8 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "DEN", 9 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "RAKE", 10 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "SLIP1", 11 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "NT1", 12 * sizeof(float), H5T_NATIVE_INT);
  H5Tinsert(dtype, "SLIP2", 12 * sizeof(float) + 1 * sizeof(int),
            H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "NT2", 13 * sizeof(float) + 1 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert(dtype, "SLIP3", 13 * sizeof(float) + 2 * sizeof(int),
            H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "NT3", 14 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_INT);

  struct srf_meta_t *srf_metadata;
  struct srf_data_t *point_data;
  float *sr_data;

  int npts = 0, nseg = 0, nsr1 = 0;
  hsize_t dims;
  double rVersion;
  int nSources = 0, nu1 = 0, nu2 = 0, nu3 = 0;

  stime = MPI_Wtime();
  // Only rank 0 reads data, then broadcast to all other processes
  if (read_color == 0) {
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, read_comm, MPI_INFO_NULL);
    fid = H5Fopen(fname, H5F_ACC_RDONLY, fapl);
    if (fid <= 0) cout << "Rupture HDF5 file " << fname << " not found" << endl;

    if (world_rank == 0) printf("Opened rupture file '%s'\n", fname);

    attr = H5Aopen(fid, "VERSION", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &rVersion);
    H5Aclose(attr);
    if (world_rank == 0) printf("Version = %.1f\n", rVersion);

    // read each header block
    attr = H5Aopen(fid, "PLANE", H5P_DEFAULT);
    aspace = H5Aget_space(attr);
    H5Sget_simple_extent_dims(aspace, &dims, NULL);
    nseg = (int)dims;
    srf_metadata =
        (struct srf_meta_t *)malloc(nseg * sizeof(struct srf_meta_t));
    H5Sclose(aspace);
    if (world_rank == 0)
      printf("Number of segments in header block: %i\n", nseg);
    H5Aread(attr, ctype, srf_metadata);
    H5Aclose(attr);

    if (world_rank == 0) {
      for (int seg = 0; seg < nseg; seg++) {
        printf("Seg #%i: elon=%g, elat=%g, nstk=%i, ndip=%i, len=%g, wid=%g\n",
               seg + 1, srf_metadata[seg].elon, srf_metadata[seg].elat,
               srf_metadata[seg].nstk, srf_metadata[seg].ndip,
               srf_metadata[seg].len, srf_metadata[seg].wid);
        printf("        stk=%g, dip=%g, dtop=%g, shyp=%g, dhyp=%g\n",
               srf_metadata[seg].stk, srf_metadata[seg].dip,
               srf_metadata[seg].dtop, srf_metadata[seg].shyp,
               srf_metadata[seg].dhyp);
      }
    }

    free(srf_metadata);
    dset = H5Dopen(fid, "POINTS", H5P_DEFAULT);
    if (dset < 0) {
      printf("Error with Rupture HDF5 file, no POINTS dataset found!\n");
      npts = -1;
    } else {
      dspace = H5Dget_space(dset);

      H5Sget_simple_extent_dims(dspace, &dims, NULL);
      npts = (int)dims;
      if (world_rank == 0)
        printf("Number of point sources in data block: %i\n", npts);

      point_data =
          (struct srf_data_t *)malloc(npts * sizeof(struct srf_data_t));
      H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, point_data);
      H5Sclose(dspace);
      H5Dclose(dset);
    }

    dset = H5Dopen(fid, "SR1", H5P_DEFAULT);
    if (dset < 0) {
      printf("Error with Rupture HDF5 file, no SR1 dataset found!\n");
      nsr1 = -1;
    } else {
      dspace = H5Dget_space(dset);

      H5Sget_simple_extent_dims(dspace, &dims, NULL);
      nsr1 = (int)dims;

      sr_data = (float *)malloc(nsr1 * sizeof(float));
      H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sr_data);
      H5Sclose(dspace);
      H5Dclose(dset);
    }

    H5Fclose(fid);
  }  // End read_color=0
  etime = MPI_Wtime();

  if (is_debug && world_rank == 0)
    printf("Read SRF-HDF5 takes %.2f seconds\n", etime - stime);

  MPI_Bcast(&npts, 1, MPI_INT, 0, node_comm);
  /* MPI_Bcast(&npts, 1, MPI_INT, 0, MPI_COMM_WORLD); */
  if (npts == -1) return;

  MPI_Bcast(&nsr1, 1, MPI_INT, 0, node_comm);
  /* MPI_Bcast(&nsr1, 1, MPI_INT, 0, MPI_COMM_WORLD); */
  if (nsr1 == -1) return;

  if (read_color != 0) {
    point_data = (struct srf_data_t *)malloc(npts * sizeof(struct srf_data_t));
    sr_data = (float *)malloc(nsr1 * sizeof(float));
  }

  MPI_Bcast(point_data, npts * sizeof(struct srf_data_t), MPI_CHAR, 0,
            node_comm);
  /* MPI_Bcast(point_data, npts*sizeof(struct srf_data_t), MPI_CHAR, 0,
   * MPI_COMM_WORLD); */

  MPI_Bcast(sr_data, nsr1, MPI_FLOAT, 0, node_comm);
  /* MPI_Bcast(sr_data, nsr1, MPI_FLOAT, 0, MPI_COMM_WORLD); */
  stime = MPI_Wtime();

  MPI_Comm_free(&node_comm);
  MPI_Comm_free(&read_comm);

  if (is_debug && world_rank == 0)
    printf("Bcast SRF-HDF5 takes %.2f seconds\n", stime - etime);

  Source *sourcePtr;
  timeDep tDep = iDiscrete;
  char formstring[100];
  strcpy(formstring, "Discrete");

  double x = 0.0, y = 0.0, z = 0.0;
  float_sw4 m0 = 1.0;
  float_sw4 t0 = 0.0, freq = 1.0;
  float_sw4 mxx = 0.0, mxy = 0.0, mxz = 0.0, myy = 0.0, myz = 0.0, mzz = 0.0;
  bool topodepth = true;
  // Discrete source time function
  float_sw4 *par = NULL;
  int *ipar = NULL;
  int npar = 0, nipar = 0, ncyc = 0, sr1pos = 0;
  // read all point sources
  for (int pts = 0; pts < npts; pts++) {
    double lon, lat, dep, stk, dip, area, tinit, dt, rake, slip1, slip2, slip3;
    int nt1, nt2, nt3;

    lon = (double)point_data[pts].lon;
    lat = (double)point_data[pts].lat;
    dep = (double)point_data[pts].dep;
    stk = (double)point_data[pts].stk;
    dip = (double)point_data[pts].dip;
    area = (double)point_data[pts].area;
    tinit = (double)point_data[pts].tinit;
    dt = (double)point_data[pts].dt;
    rake = (double)point_data[pts].rake;
    slip1 = (double)point_data[pts].slip1;
    slip2 = (double)point_data[pts].slip2;
    slip3 = (double)point_data[pts].slip3;
    nt1 = (int)point_data[pts].nt1;
    nt2 = (int)point_data[pts].nt2;
    nt3 = (int)point_data[pts].nt3;

    // nothing to do if nt1=nt2=nt3=0
    if (nt1 <= 0 && nt2 <= 0 && nt3 <= 0) continue;

    if (world_rank == 0 && mVerbose >= 2) {
      printf(
          "point #%i: lon=%g, lat=%g, dep=%g, stk=%g, dip=%g, area=%g, "
          "tinit=%g, dt=%g\n",
          pts + 1, lon, lat, dep, stk, dip, area, tinit, dt);
      printf(
          "          rake=%g, slip1=%g, nt1=%i, slip2=%g, nt2=%i, slip3=%g, "
          "nt3=%i\n",
          rake, slip1, nt1, slip2, nt2, slip3, nt3);
    }

    // read discrete time series for u1
    if (nt1 > 0) {
      nu1++;
      // note that the first data point is always zero, but the last is not
      // for this reason we always pad the time zeries with a '0'
      // also note that we need at least 7 data points, i.e. nt1>=6
      int nt1dim = max(6, nt1);
      par = new float_sw4[nt1dim + 2];
      par[0] = tinit;
      t0 = tinit;
      freq = 1 / dt;
      ipar = new int[1];
      ipar[0] = nt1dim + 1;  // add an extra point

      for (int i = 0; i < nt1; i++) par[i + 1] = sr_data[sr1pos++];

      // pad with 0
      if (nt1 < 6) {
        for (int j = nt1; j < 6; j++) par[j + 1] = 0.;
      }

      // last 0
      par[nt1dim + 1] = 0.0;

      // scale cm/s to m/s
      for (int i = 1; i <= nt1dim + 1; i++) {
        par[i] *= 1e-2;
      }

      // AP: Mar. 1, 2016: Additional scaling is needed to make the integral of
      // the time function = 1
      float_sw4 slip_m = slip1 * 1e-2;
      float_sw4 slip_sum = 0;
      for (int i = 1; i <= nt1dim + 1; i++) {
        slip_sum += par[i];
      }
      slip_sum *= dt;

      if (world_rank == 0 && mVerbose >= 2) {
        printf(
            "INFO: SRF file: dt*sum(slip_vel)=%e [m], total slip (from "
            "header)=%e [m]\n",
            slip_sum, slip_m);
      }
      // scale time series to sum to integrate to one
      for (int i = 1; i <= nt1dim + 1; i++) {
        par[i] /= slip_sum;
      }
      if (world_rank == 0 && mVerbose >= 2) {
        slip_sum = 0;
        for (int i = 1; i <= nt1dim + 1; i++) {
          slip_sum += par[i];
        }
        slip_sum *= dt;
        printf(
            "INFO: SRF file: After scaling time series: dt*sum(par)=%e [m]\n",
            slip_sum);
      }
      // done scaling

      npar = nt1dim + 2;
      nipar = 1;

      // printf("Read discrete time series: tinit=%g, dt=%g, nt1=%i\n", tinit,
      // dt, nt1); for (int i=0; i<nt1+1; i++)
      //   printf("Sv1[%i]=%g\n", i+1, par[i+1]);

      // convert lat, lon, depth to (x,y,z)
      ew->computeCartesianCoord(x, y, lon, lat);
      // convert depth in [km] to [m]
      z = dep * 1e3;

      // convert strike, dip, rake to Mij
      float_sw4 radconv = M_PI / 180.;
      float_sw4 S, D, R;
      stk -= mGeoAz;  // subtract off the grid azimuth
      S = stk * radconv;
      D = dip * radconv;
      R = rake * radconv;

      mxx = -1.0 * (sin(D) * cos(R) * sin(2 * S) +
                    sin(2 * D) * sin(R) * sin(S) * sin(S));
      myy = (sin(D) * cos(R) * sin(2 * S) -
             sin(2 * D) * sin(R) * cos(S) * cos(S));
      mzz = -1.0 * (mxx + myy);
      mxy = (sin(D) * cos(R) * cos(2 * S) +
             0.5 * sin(2 * D) * sin(R) * sin(2 * S));
      mxz = -1.0 * (cos(D) * cos(R) * cos(S) + cos(2 * D) * sin(R) * sin(S));
      myz = -1.0 * (cos(D) * cos(R) * sin(S) - cos(2 * D) * sin(R) * cos(S));

      // scale (note that the shear modulus is not yet available. Also note that
      // we convert [cm] to [m])
      m0 = area * 1e-4 * slip1 * 1e-2;

      mxx *= m0;
      mxy *= m0;
      mxz *= m0;
      myy *= m0;
      myz *= m0;
      mzz *= m0;

      // before creating the source, make sure (x,y,z) is inside the
      // computational domain

      // only check the z>zmin when we have topography. For a flat free surface,
      // we will remove sources too close or above the surface in the call to
      // mGlobalUniqueSources[i]->correct_Z_level()

      if (x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax ||
          z < zmin || z > m_global_zmax) {
        stringstream sourceposerr;
        sourceposerr << endl
                     << "***************************************************"
                     << endl
                     << " ERROR:  Source positioned outside grid!  " << endl
                     << endl
                     << " Source from rupture file @" << endl
                     << "  x=" << x << " y=" << y << " z=" << z << endl
                     << "  lat=" << lat << " lon=" << lon << " dep=" << dep
                     << endl
                     << endl;

        if (x < xmin)
          sourceposerr << " x is " << xmin - x << " meters away from min x ("
                       << xmin << ")" << endl;
        else if (x > m_global_xmax)
          sourceposerr << " x is " << x - m_global_xmax
                       << " meters away from max x (" << m_global_xmax << ")"
                       << endl;
        if (y < ymin)
          sourceposerr << " y is " << ymin - y << " meters away from min y ("
                       << ymin << ")" << endl;
        else if (y > m_global_ymax)
          sourceposerr << " y is " << y - m_global_ymax
                       << " meters away from max y (" << m_global_ymax << ")"
                       << endl;
        if (z < zmin)
          sourceposerr << " z is " << zmin - z << " meters away from min z ("
                       << zmin << ")" << endl;
        else if (z > m_global_zmax)
          sourceposerr << " z is " << z - m_global_zmax
                       << " meters away from max z (" << m_global_zmax << ")"
                       << endl;
        sourceposerr << "***************************************************"
                     << endl;
        if (world_rank == 0) cout << sourceposerr.str();
      } else {
        sourcePtr =
            new Source(ew, freq, t0, x, y, z, mxx, mxy, mxz, myy, myz, mzz,
                       tDep, formstring, topodepth, ncyc, par, npar, ipar,
                       nipar, true);  // true is correctStrengthForMu

        if (sourcePtr->ignore()) {
          delete sourcePtr;
        } else {
          a_GlobalUniqueSources[event].push_back(sourcePtr);
          nSources++;
        }
      }

      // deallocate temporary arrays...
      delete[] par;
      delete[] ipar;

    }  // end if nt1 >0

    // read past discrete time series for u2
    if (nt2 > 0) {
      nu2++;
      /* double dum; */
      if (world_rank == 0) printf("WARNING nt2=%i > 0 will be ignored\n", nt2);
    }  // end if nt2 > 0

    // read past discrete time series for u3
    if (nt3 > 0) {
      nu3++;
      /* double dum; */
      if (world_rank == 0) printf("WARNING nt3=%i > 0 will be ignored\n", nt3);
    }  // end if nt3 > 0

  }  // end for all sources
  if (world_rank == 0)
    printf(
        "Read npts=%i, made %i point moment tensor sources, nu1=%i, nu2=%i, "
        "nu3=%i\n",
        npts, nSources, nu1, nu2, nu3);

  etime = MPI_Wtime();
  if (is_debug && world_rank == 0)
    printf("Create source takes %.2f seconds\n", etime - stime);

  H5Tclose(ctype);
  H5Tclose(dtype);
  free(point_data);
  free(sr_data);
}

#endif  // USE_HDF5
#endif  // READHDF5_C
