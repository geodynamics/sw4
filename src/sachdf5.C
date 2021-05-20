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
#ifndef SACHDF5_C
#define SACHDF5_C

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

#define USE_DSET_ATTR 1

#ifdef USE_HDF5

#include "sachdf5.h"

int createAttr(hid_t loc, const char *name, hid_t type_id, hid_t space_id) {
  hid_t attr, dcpl;
  herr_t ret;

#ifdef USE_DSET_ATTR
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
  H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
  attr =
      H5Dcreate(loc, name, type_id, space_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  H5Pclose(dcpl);
#else
  attr = H5Acreate(loc, name, type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (attr < 0) {
    printf("%s: Error with H5Acreate [%s]\n", __func__, name);
    return -1;
  }
#ifdef USE_DSET_ATTR
  H5Dclose(attr);
#else
  H5Aclose(attr);
#endif
  return 1;
}

int createWriteAttr(hid_t loc, char const *name, hid_t type_id, hid_t space_id,
                    void *data) {
  hid_t attr, dcpl;
  herr_t ret;

#ifdef USE_DSET_ATTR
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
  H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
  attr =
      H5Dcreate(loc, name, type_id, space_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  H5Pclose(dcpl);
#else
  attr = H5Acreate(loc, name, type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (attr < 0) {
    printf("%s: Error with H5Acreate [%s]\n", __func__, name);
    return -1;
  }
#ifdef USE_DSET_ATTR
  ret = H5Dwrite(attr, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#else
  ret = H5Awrite(attr, type_id, data);
#endif
  if (ret < 0) {
    printf("%s: Error with H5Awrite [%s]\n", __func__, name);
    return -1;
  }
#ifdef USE_DSET_ATTR
  H5Dclose(attr);
#else
  H5Aclose(attr);
#endif
  return 1;
}

int openWriteAttr(hid_t loc, const char *name, hid_t type_id, void *data) {
  hid_t attr;
  herr_t ret;

  // debug
  /* printf("%s: start [%s]\n", __func__, name); */

#ifdef USE_DSET_ATTR
  attr = H5Dopen(loc, name, H5P_DEFAULT);
#else
  attr = H5Aopen(loc, name, H5P_DEFAULT);
#endif
  if (attr < 0) {
    printf("%s: Error with H5Aopen [%s]\n", __func__, name);
    return -1;
  }

#ifdef USE_DSET_ATTR
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  ret = H5Dwrite(attr, type_id, H5S_ALL, H5S_ALL, dxpl, data);
  H5Pclose(dxpl);
#else
  ret = H5Awrite(attr, type_id, data);
#endif
  if (ret < 0) {
    printf("%s: Error with H5Awrite [%s]\n", __func__, name);
    return -1;
  }

#ifdef USE_DSET_ATTR
  H5Dclose(attr);
#else
  H5Aclose(attr);
#endif

  // debug
  /* printf("%s: write [%s] success int=%d, float=%.2f!\n", __func__, name,
   * *((int*)data), *((float*)data)); */
  return 1;
}

int createWriteAttrStr(hid_t loc, const char *name, const char *str) {
  hid_t attr, atype, space;
  herr_t ret;

  ASSERT(name);
  ASSERT(str);

  space = H5Screate(H5S_SCALAR);
  atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, strlen(str) + 1);
  H5Tset_strpad(atype, H5T_STR_NULLTERM);
  attr = H5Acreate(loc, name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0) {
    printf("Error with H5Acreate [%s]\n", name);
    return -1;
  }

  ret = H5Awrite(attr, atype, str);
  if (ret < 0) {
    printf("Error with H5Awrite [%s]\n", name);
    return -1;
  }

  H5Tclose(atype);
  H5Sclose(space);
  H5Aclose(attr);

  return 1;
}

int openWriteData(hid_t loc, const char *name, hid_t type_id, void *data,
                  int ndim, hsize_t *start, hsize_t *count, int total_npts,
                  float btime, float cmpinc, float cmpaz, bool isIncAzWritten,
                  bool isLast) {
  bool is_debug = false;
  /* is_debug = true; */
  double stime, etime, etime1;
  hid_t dset, filespace, dxpl;
  herr_t ret;
  hsize_t dims[3];

  /* stime = MPI_Wtime(); */

  dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);

  // debug
  /* printf("%s: start write dset [%s]!\n", __func__, name); */

  dset = H5Dopen(loc, name, H5P_DEFAULT);
  if (dset < 0) {
    printf("%s: Error with H5Dopen [%s]\n", __func__, name);
    return -1;
  }

  filespace = H5Dget_space(dset);
  H5Sget_simple_extent_dims(filespace, dims, NULL);
  if (dims[0] < start[0] + count[0]) count[0] = dims[0] - start[0];
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

  ret = H5Dwrite(dset, type_id, H5S_ALL, filespace, dxpl, data);
  if (ret < 0) {
    printf("%s: Error with H5Dwrite [%s]\n", __func__, name);
    return -1;
  }

  /* etime = MPI_Wtime(); */

  if (isLast) openWriteAttr(loc, "NPTS", H5T_NATIVE_INT, &total_npts);

  /* etime1 = MPI_Wtime(); */

  /* int myRank; */
  /* MPI_Comm_rank(MPI_COMM_WORLD, &myRank); */
  /* printf("Rank %d: Attr write time %f, data write time %f\n", myRank,
   * etime1-etime, etime-stime); */
  /* fflush(stdout); */

  if (is_debug) {
    float *write_data = (float *)data;
    printf("npts=%d, first=%e, second=%e, middle=%e, second last=%e, last=%e\n",
           total_npts, write_data[0], write_data[1], write_data[total_npts / 2],
           write_data[total_npts - 2], write_data[total_npts - 1]);
    fflush(stdout);
  }

  /* H5Dflush(dset); */
  H5Pclose(dxpl);
  if (filespace != H5S_ALL) H5Sclose(filespace);
  H5Dclose(dset);

  // debug
  /* printf("%s: write [%s] success!\n", __func__, name); */
  /* fflush(stdout); */
  return 1;
}

int createTimeSeriesHDF5File(vector<TimeSeries *> &TimeSeries, int totalSteps,
                             float_sw4 delta, string suffix) {
  bool is_debug = false;

  hid_t fid, grp, attr, attr_space1, attr_space3, dset_space, dset, dcpl, fapl;
  herr_t ret;
  hsize_t dims1 = 1, dims3 = 3, total_dims;
  double start_time, elapsed_time;
  float stxyz[3], stlonlatdep[3], o, dt;
  std::string dset_names[10];
  TimeSeries::receiverMode mode;
  bool xyzcomponent;
  int ndset = 0, isnsew = 0;

  if (TimeSeries.size() == 0) return 0;

  /* // Do not create file for inverse */
  /* if (TimeSeries[0]->isInverse()) */
  /*     return 0; */

  start_time = MPI_Wtime();

  std::string path = TimeSeries[0]->getPath();
  std::string name = TimeSeries[0]->gethdf5FileName();
  std::string filename;

  char setstripe[4096], *env;
  int disablestripe = 0, stripecount = 128, stripesize = 512;

  env = getenv("DISABLE_LUSTRE_STRIPE");
  if (env != NULL) disablestripe = atoi(env);

  // Cori Lustre has cscratch in path
  if (path.find("cscratch") == std::string::npos) disablestripe = 1;

  // Set stripe parameters for time-series data
  if (disablestripe != 1) {
    env = getenv("SAC_LUSTRE_STRIPE_COUNT");
    if (env != NULL) stripecount = atoi(env);
    env = getenv("SAC_LUSTRE_STRIPE_SIZE");
    if (env != NULL) stripesize = atoi(env);

    if (stripecount < 1) stripecount = 1;
    if (stripesize < 128) stripesize = 512;

    fflush(stdout);
    sprintf(setstripe, "lfs setstripe -c %d -S %dk %s", stripecount, stripesize,
            path.c_str());
    if (system(setstripe) != 0) {
      disablestripe = 1;
      printf(
          "Failed to set Lustre stripe for SAC-HDF5 file, sw4 will continue to "
          "run (set DISABLE_LUSTRE_STRIPE=1 to disable this and the above lfs "
          "error messages, or set valid values to SAC_LUSTRE_STRIPE_COUNT and "
          "SAC_LUSTRE_STRIPE_SIZE)\n");
    } else
      printf("For SAC-HDF5 files Lustre stripe set to: %s\n", setstripe);
    fflush(stdout);
  }

  // Build the file name
  if (path != ".") filename = path;

  filename.append(name);
  filename.append(suffix);

  if (filename.find(".hdf5") == string::npos &&
      filename.find(".h5") == string::npos)
    filename.append(".hdf5");

  if (is_debug) {
    printf("Start create time-history HDF5 file [%s], %d steps\n",
           filename.c_str(), totalSteps);
    fflush(stdout);
  }

  if (access(filename.c_str(), F_OK) != -1) {
    // if the file exists, move it to a .bak before writing
    std::string bak;
    bak = filename + ".bak";
    ret = rename(filename.c_str(), bak.c_str());
    cout << "Rename existing file to [" << bak.c_str() << "]" << endl;
    if (ret == -1)
      cout << "ERROR: renaming SAC-HDF5 file to " << bak.c_str() << endl;
  }

  int alignment = 65536;
  /* char *env = getenv("HDF5_ALIGNMENT_SIZE"); */
  /* if (env != NULL) */
  /*     alignment = atoi(env); */
  /* if (alignment < 65536) */
  /*     alignment = 65536; */

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_alignment(fapl, 32767, alignment);
  fid = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  if (fid < 0) {
    printf("Error: H5Fcreate failed\n");
    return -1;
  }
  H5Pclose(fapl);

  // Set stripe parameters for files created after sac file (e.g. images)
  if (disablestripe != 1) {
    stripecount = 128, stripesize = 1024;
    env = getenv("IMAGE_LUSTRE_STRIPE_COUNT");
    if (env != NULL) stripecount = atoi(env);
    env = getenv("IMAGE_LUSTRE_STRIPE_SIZE");
    if (env != NULL) stripesize = atoi(env);

    if (stripecount < 1) stripecount = 1;
    if (stripesize < 1024) stripesize = 1024;

    fflush(stdout);
    sprintf(setstripe, "lfs setstripe -c %d -S %dk %s", stripecount, stripesize,
            path.c_str());
    if (system(setstripe) != 0)
      printf(
          "Failed to set Lustre stripe for files other than SAC-HDF5, sw4 will "
          "continue to run (set DISABLE_LUSTRE_STRIPE=1 to disable this and "
          "the above lfs error messages, or set valid values to "
          "IMAGE_LUSTRE_STRIPE_COUNT and IMAGE_LUSTRE_STRIPE_SIZE)\n");
    else
      printf("For files other than SAC-HDF5 Lustre stripe set to: %s\n",
             setstripe);
    fflush(stdout);
  }

  attr_space1 = H5Screate_simple(1, &dims1, NULL);
  attr_space3 = H5Screate_simple(1, &dims3, NULL);

  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
  H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

  // datetime, human readable string for date and time of start of SW4
  // calculation, e.g. 2019-07-23T01:02:03)
  char utcstr[32];
  sprintf(utcstr, "%d-%02d-%02dT%02d:%02d:%02d.%06d", TimeSeries[0]->getMUTC(0),
          TimeSeries[0]->getMUTC(1), TimeSeries[0]->getMUTC(2),
          TimeSeries[0]->getMUTC(3), TimeSeries[0]->getMUTC(4),
          TimeSeries[0]->getMUTC(5), TimeSeries[0]->getMUTC(6));
  createWriteAttrStr(fid, "DATETIME", utcstr);

  // delta, sample interval of time-series (seconds)
  int downsample = TimeSeries[0]->getDownSample();
  if (downsample < 1) downsample = 1;

  dt = (float)delta * downsample;
  createWriteAttr(fid, "DELTA", H5T_NATIVE_FLOAT, attr_space1, &dt);
  createWriteAttr(fid, "DOWNSAMPLE", H5T_NATIVE_INT, attr_space1, &downsample);

  // o, origin time (seconds, relative to start time of SW4 calculation and
  // seismogram, earliest source)
  createAttr(fid, "ORIGINTIME", H5T_NATIVE_FLOAT, attr_space1);

  // units, units for motion (m for displacement, m/s for velocity, m/s/s for
  // acceleration
  mode = TimeSeries[0]->getMode();
  if (mode == TimeSeries::Displacement)
    createWriteAttrStr(fid, "UNIT", "m");
  else if (mode == TimeSeries::Velocity)
    createWriteAttrStr(fid, "UNIT", "m/s");
  /* else if( mode == TimeSeries::??) */
  /* createWriteAttrStr(fid, "Unit", "m/s/s"); */

  float cmpazs[9] = {0}, cmpincs[9] = {0};

  for (int ts = 0; ts < TimeSeries.size(); ts++) {
    std::string stationname = TimeSeries[ts]->getStationName();
    grp = H5Gcreate(fid, stationname.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
    if (grp < 0) {
      printf("Error: H5Gcreate [%s] failed\n", stationname.c_str());
      return -1;
    }

    // Number of points
    createAttr(grp, "NPTS", H5T_NATIVE_INT, attr_space1);

    // x, y, z
    createAttr(grp, "STX,STY,STZ", H5T_NATIVE_DOUBLE, attr_space3);

    // Lon, lat, dep
    createAttr(grp, "STLA,STLO,STDP", H5T_NATIVE_DOUBLE, attr_space3);

    // TODO: Location, no value to write now
    createAttr(grp, "LOC", H5T_NATIVE_INT, attr_space1);

    xyzcomponent = TimeSeries[ts]->getXYZcomponent();
    if (!xyzcomponent) isnsew = 1;

    createWriteAttr(grp, "ISNSEW", H5T_NATIVE_INT, attr_space1, &isnsew);

    cmpazs[0] = TimeSeries[ts]->getXaz();
    cmpazs[1] = TimeSeries[ts]->getXaz() + 90.;
    cmpazs[2] = 0.;
    cmpincs[0] = 90.;
    cmpincs[1] = 90.;
    cmpincs[2] = 180.;
    mode = TimeSeries[ts]->getMode();
    // Datasets
    if (mode == TimeSeries::Displacement) {
      ndset = 3;
      if (xyzcomponent) {
        dset_names[0] = "X";
        dset_names[1] = "Y";
        dset_names[2] = "Z";
      } else {
        dset_names[0] = "EW";
        dset_names[1] = "NS";
        dset_names[2] = "UP";
        cmpincs[2] = 0.;
      }
    } else if (mode == TimeSeries::Velocity) {
      ndset = 3;
      if (xyzcomponent) {
        dset_names[0] = "Vx";
        dset_names[1] = "Vy";
        dset_names[2] = "Vz";
      } else {
        dset_names[0] = "Vew";
        dset_names[1] = "Vns";
        dset_names[2] = "Vup";
        cmpincs[2] = 0.;
      }
    } else if (mode == TimeSeries::Div) {
      ndset = 1;
      dset_names[0] = "Div";
    } else if (mode == TimeSeries::Curl) {
      ndset = 3;
      dset_names[0] = "Curlx";
      dset_names[1] = "Curly";
      dset_names[2] = "Curlz";
    } else if (mode == TimeSeries::Strains) {
      ndset = 6;
      dset_names[0] = "Uxx";
      dset_names[1] = "Uyy";
      dset_names[2] = "Uzz";
      dset_names[3] = "Uxy";
      dset_names[4] = "Uxz";
      dset_names[5] = "Uyz";
    } else if (mode == TimeSeries::DisplacementGradient) {
      ndset = 9;
      dset_names[0] = "DUXDX";
      dset_names[1] = "DUXDY";
      dset_names[2] = "DUXDZ";
      dset_names[3] = "DUYDX";
      dset_names[4] = "DUYDY";
      dset_names[5] = "DUYDZ";
      dset_names[6] = "DUZDX";
      dset_names[7] = "DUZDY";
      dset_names[8] = "DUZDZ";
    }

    for (int i = 0; i < ndset; i++) {
      total_dims = (hsize_t)(totalSteps / TimeSeries[ts]->getDownSample());
      if (total_dims == 0) {
        printf("%s: Error with dataset length 0\n", __func__);
        return -1;
      }
      dset_space = H5Screate_simple(1, &total_dims, NULL);
      dset = H5Dcreate(grp, dset_names[i].c_str(), H5T_NATIVE_FLOAT, dset_space,
                       H5P_DEFAULT, dcpl, H5P_DEFAULT);
      H5Sclose(dset_space);
      ASSERT(dset >= 0);
#ifdef USE_DSET_ATTR
      std::string incname = dset_names[i] + "CMPINC";
      std::string azname = dset_names[i] + "CMPAZ";
      /* createAttr(grp, incname.c_str(), H5T_NATIVE_FLOAT, attr_space1); */
      /* createAttr(grp, azname.c_str(), H5T_NATIVE_FLOAT, attr_space1); */
      createWriteAttr(grp, incname.c_str(), H5T_NATIVE_FLOAT, attr_space1,
                      &cmpazs[i]);
      createWriteAttr(grp, azname.c_str(), H5T_NATIVE_FLOAT, attr_space1,
                      &cmpincs[i]);
#else
      /* createAttr(dset, "CMPINC", H5T_NATIVE_FLOAT, attr_space1); */
      /* createAttr(dset, "CMPAZ", H5T_NATIVE_FLOAT, attr_space1); */
      createWriteAttr(grp, "CMPAZ", H5T_NATIVE_FLOAT, attr_space1, &cmpazs[i]);
      createWriteAttr(grp, "CMPINC", H5T_NATIVE_FLOAT, attr_space1,
                      &cmpincs[i]);
#endif
      H5Dclose(dset);
    }
    H5Gclose(grp);

    TimeSeries[ts]->setNsteps(totalSteps);
  }

  H5Pclose(dcpl);
  H5Sclose(attr_space1);
  H5Sclose(attr_space3);
  H5Fclose(fid);

  elapsed_time = MPI_Wtime() - start_time;
  printf("Created SAC HDF5 file [%s] time %e seconds\n", filename.c_str(),
         elapsed_time);
  fflush(stdout);
  return 1;
}

int readAttrStr(hid_t loc, const char *name, char *str) {
  herr_t ret;
  hid_t attr, atype;

  ASSERT(name);
  ASSERT(str);

  attr = H5Aopen(loc, name, H5P_DEFAULT);
  if (attr < 0) {
    printf("Error with H5Aopen [%s]\n", name);
    return -1;
  }
  atype = H5Aget_type(attr);

  ret = H5Aread(attr, atype, str);
  if (ret < 0) {
    printf("Error with H5Aread [%s]\n", name);
    return -1;
  }

  H5Tclose(atype);
  H5Aclose(attr);

  return 1;
}

int readAttrInt(hid_t loc, const char *name, int *data) {
  hid_t attr;
  herr_t ret;

#ifdef USE_DSET_ATTR
  attr = H5Dopen(loc, name, H5P_DEFAULT);
#else
  attr = H5Aopen(loc, name, H5P_DEFAULT);
#endif
  if (attr < 0) {
    printf("%s: Error with H5Aopen [%s]\n", __func__, name);
    return -1;
  }

#ifdef USE_DSET_ATTR
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  ret = H5Dread(attr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, dxpl, (void *)data);
  H5Pclose(dxpl);
#else
  ret = H5Aread(attr, H5T_NATIVE_INT, (void *)data);
#endif
  if (ret < 0) {
    printf("%s: Error with H5Aread [%s]\n", __func__, name);
    return -1;
  }

#ifdef USE_DSET_ATTR
  H5Dclose(attr);
#else
  H5Aclose(attr);
#endif

  return 1;
}

int readAttrFloat(hid_t loc, const char *name, float *data) {
  hid_t attr;
  herr_t ret;

#ifdef USE_DSET_ATTR
  attr = H5Dopen(loc, name, H5P_DEFAULT);
#else
  attr = H5Aopen(loc, name, H5P_DEFAULT);
#endif
  if (attr < 0) {
    printf("%s: Error with H5Aopen [%s]\n", __func__, name);
    return -1;
  }

#ifdef USE_DSET_ATTR
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  ret = H5Dread(attr, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, dxpl, (void *)data);
  H5Pclose(dxpl);
#else
  ret = H5Aread(attr, H5T_NATIVE_FLOAT, (void *)data);
#endif
  if (ret < 0) {
    printf("%s: Error with H5Aread [%s]\n", __func__, name);
    return -1;
  }

#ifdef USE_DSET_ATTR
  H5Dclose(attr);
#else
  H5Aclose(attr);
#endif

  return 1;
}

int readHDF5Data(hid_t loc, const char *name, int npts, void *data) {
  hid_t dset, filespace, dxpl;
  hsize_t start, count;
  herr_t ret;

  dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);

  dset = H5Dopen(loc, name, H5P_DEFAULT);
  if (dset < 0) {
    printf("%s: Error with H5Dopen [%s]\n", __func__, name);
    return -1;
  }

  start = 0;
  count = (hsize_t)npts;
  filespace = H5Dget_space(dset);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL);

  ret = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, filespace, dxpl, data);
  if (ret < 0) {
    printf("%s: Error with H5Dread [%s]\n", __func__, name);
    return -1;
  }

  H5Pclose(dxpl);
  H5Sclose(filespace);
  H5Dclose(dset);

  return 1;
}

#endif  // USE_HDF5

#endif  // SACHDF5_C
