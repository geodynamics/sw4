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
#ifndef SACHDF5_C
#define SACHDF5_C

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#include <ctime>

#include "Require.h"
#include "TimeSeries.h"

#ifdef USE_HDF5

#include "sachdf5.h"

#define USE_DSET_ATTR 1

int createAttr(hid_t loc, const char *name, hid_t type_id, hid_t space_id)
{
    hid_t attr, dcpl;
    herr_t ret;

#ifdef USE_DSET_ATTR
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
    attr = H5Dcreate(loc, name, type_id, space_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);
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

int createWriteAttr(hid_t loc, const char *name, hid_t type_id, hid_t space_id, void *data)
{
    hid_t attr, dcpl;
    herr_t ret;

#ifdef USE_DSET_ATTR
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
    attr = H5Dcreate(loc, name, type_id, space_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Pclose(dcpl);
#else
    attr = H5Acreate(loc, name, type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
#endif
    if (attr < 0) {
        printf("%s: Error with H5Acreate [%s]\n", __func__, name);
        return -1;
    }
#ifdef USE_DSET_ATTR
    ret  = H5Dwrite(attr, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#else
    ret  = H5Awrite(attr, type_id, data);
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

int openWriteAttr(hid_t loc, const char *name, hid_t type_id, void *data)
{
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
    ret  = H5Dwrite(attr, type_id, H5S_ALL, H5S_ALL, dxpl, data);
    H5Pclose(dxpl);
#else
    ret  = H5Awrite(attr, type_id, data);
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
    /* printf("%s: write [%s] success int=%d, float=%.2f!\n", __func__, name, *((int*)data), *((float*)data)); */
    return 1;
}

int createWriteAttrStr(hid_t loc, const char *name, const char* str)
{
    hid_t attr, atype, space;
    herr_t ret;

    ASSERT(name);
    ASSERT(str);

    space = H5Screate(H5S_SCALAR);
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(str)+1);
    H5Tset_strpad(atype,H5T_STR_NULLTERM);
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

int openWriteData(hid_t loc, const char *name, hid_t type_id, void *data, int ndim, hsize_t *start, hsize_t *count, int total_npts, 
                  float btime, float cmpinc, float cmpaz, bool isIncAzWritten, bool isLast)
{
    hid_t dset, filespace, dxpl;
    herr_t ret;

    dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);

    // debug
    /* printf("%s: start write dset [%s]!\n", __func__, name); */

    dset = H5Dopen(loc, name, H5P_DEFAULT);
    if (dset < 0) {
        printf("%s: Error with H5Dopen [%s]\n", __func__, name);
        return -1;
    }

    if (!isIncAzWritten) {
#ifdef USE_DSET_ATTR
      std::string newname = name;
      std::string incname = newname + "CMPINC";
      std::string azname  = newname + "CMPAZ";
      openWriteAttr(loc, incname.c_str(), H5T_NATIVE_FLOAT, &cmpinc);
      openWriteAttr(loc, azname.c_str(), H5T_NATIVE_FLOAT, &cmpaz);
#else
      openWriteAttr(dset, "CMPINC", H5T_NATIVE_FLOAT, &cmpinc);
      openWriteAttr(dset, "CMPAZ", H5T_NATIVE_FLOAT, &cmpaz);
#endif
    }

    if (start[0] == 0) {
        filespace = H5S_ALL;
    }
    else {
        filespace = H5Dget_space(dset);
        H5Sselect_hyperslab (filespace, H5S_SELECT_SET, start, NULL, count, NULL);
    }

    ret  = H5Dwrite(dset, type_id, H5S_ALL, filespace, dxpl, data);
    if (ret < 0) {
        printf("%s: Error with H5Dwrite [%s]\n", __func__, name);
        return -1;
    }

    if (isLast) {
        openWriteAttr(loc, "NPTS", H5T_NATIVE_INT, &total_npts);
    }

    H5Pclose(dxpl);
    if (filespace != H5S_ALL) 
        H5Sclose(filespace);
    H5Dclose(dset);

    // debug
    /* printf("%s: write [%s] success!\n", __func__, name); */
    /* fflush(stdout); */
    return 1;
}


int createTimeSeriesHDF5File(vector<TimeSeries*> & TimeSeries, int totalSteps, float_sw4 delta)
{
  bool is_debug = false;

  hid_t fid, grp, attr, attr_space1, attr_space3, dset_space, dset, dcpl;
  herr_t ret;
  hsize_t dims1 = 1, dims3 = 3, total_dims;
  double start_time, elapsed_time;
  float stxyz[3], stlonlatdep[3], o, dt;
  std::string dset_names[10];
  TimeSeries::receiverMode mode;
  bool xyzcomponent;
  int ndset = 0;

  if (TimeSeries.size() == 0) 
      return 0;

  start_time = MPI_Wtime();

  std::string path = TimeSeries[0]->getPath();
  std::string name = TimeSeries[0]->getFileName();
  std::string filename;

  // Build the file name
  if( path != "." )
    filename = path;

  filename.append(name);
  if (name.find(".hdf5") == string::npos && name.find(".h5") == string::npos) 
    filename.append(".hdf5");
 
  if (is_debug) {
      printf("Start createTimeSeriesHDF5File [%s], %d steps\n", filename.c_str(), totalSteps);
      fflush(stdout);
  }

  if( access( filename.c_str(), F_OK ) != -1) {
    // if the file exists, move it to a .bak before writing
    std::string bak;
    bak =  filename + ".bak";
    ret = rename(filename.c_str(), bak.c_str());
    if( ret == -1 )
      cout << "ERROR: renaming SAC HDF5 file to " << bak.c_str() <<  endl;
  }

  fid = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) {
     printf("Error: H5Fcreate failed\n");
     return -1;
  }

  attr_space1 = H5Screate_simple(1, &dims1, NULL);
  attr_space3 = H5Screate_simple(1, &dims3, NULL);

  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
  H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

  // datetime, human readable string for date and time of start of SW4 calculation, e.g. 2019-07-23T01:02:03)
  char utcstr[32];
  sprintf(utcstr, "%d-%02d-%02dT%02d:%02d:%02d.%06d", TimeSeries[0]->getMUTC(0), TimeSeries[0]->getMUTC(1), TimeSeries[0]->getMUTC(2), TimeSeries[0]->getMUTC(3), TimeSeries[0]->getMUTC(4), TimeSeries[0]->getMUTC(5), TimeSeries[0]->getMUTC(6));
  /* printf("DATATIME=%s\n", utcstr); */
  createWriteAttrStr(fid, "DATETIME", utcstr);
 
  // delta, sample interval of time-series (seconds)
  dt = (float)delta;
  dt *= TimeSeries[0]->getDownSample();
  createWriteAttr(fid, "DELTA", H5T_NATIVE_FLOAT, attr_space1, &dt);

  // o, origin time (seconds, relative to start time of SW4 calculation and seismogram, earliest source)
  /* o = (float)TimeSeries[0]->getEpiTimeOffset(); */
  createAttr(fid, "ORIGINTIME", H5T_NATIVE_FLOAT, attr_space1);

  /* int writesteps = totalSteps/TimeSeries[0]->getDownSample(); */
  /* createWriteAttr(fid, "NPTS", H5T_NATIVE_INT, attr_space1, &writesteps); */

  // units, units for motion (m for displacement, m/s for velocity, m/s/s for acceleration
  mode = TimeSeries[0]->getMode();
  if( mode == TimeSeries::Displacement )
      createWriteAttrStr(fid, "UNIT", "m");
  else if( mode == TimeSeries::Velocity)
      createWriteAttrStr(fid, "UNIT", "m/s");
  /* else if( mode == TimeSeries::??) */
      /* createWriteAttrStr(fid, "Unit", "m/s/s"); */


  for (int ts=0; ts<TimeSeries.size(); ts++)
  {
    std::string stationname = TimeSeries[ts]->getStationName();
    grp  = H5Gcreate(fid, stationname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (grp < 0) {
       printf("Error: H5Gcreate [%s] failed\n", stationname.c_str());
       return -1;
    }

    // Number of points
    createAttr(grp, "NPTS", H5T_NATIVE_INT, attr_space1);

    // x, y, z
    /* stxyz[0] = float(TimeSeries[ts]->getX()); */
    /* stxyz[1] = float(TimeSeries[ts]->getY()); */
    /* stxyz[2] = float(TimeSeries[ts]->getZ()); */
    /* createWriteAttr(grp, "STX,STY,STDP", H5T_NATIVE_FLOAT, attr_space3, stxyz); */
    createAttr(grp, "STX,STY,STZ", H5T_NATIVE_FLOAT, attr_space3);

    // Lon, lat, dep
    /* stlonlatdep[0] = float(TimeSeries[ts]->getLat()); */
    /* stlonlatdep[1] = float(TimeSeries[ts]->getLon()); */
    /* stlonlatdep[2] = float(TimeSeries[ts]->getZ()); */
    /* createWriteAttr(grp, "STLA,STLO,STDP", H5T_NATIVE_FLOAT, attr_space3, stlonlatdep); */
    createAttr(grp, "STLA,STLO,STDP", H5T_NATIVE_FLOAT, attr_space3);

    // TODO: Location, no value to write now
    createAttr(grp, "LOC", H5T_NATIVE_INT, attr_space1);

    xyzcomponent = TimeSeries[ts]->getXYZcomponent();
    mode         = TimeSeries[ts]->getMode();
    // Datasets
    if( mode == TimeSeries::Displacement )
    {
        ndset = 3;
       if( xyzcomponent )
       {
          dset_names[0] = "X";
          dset_names[1] = "Y";
          dset_names[2] = "Z";
          /* msg << "[x|y|z]" << endl; */
       }
       else
       {
          dset_names[0] = "EW";
          dset_names[1] = "NS";
          dset_names[2] = "UP";
          /* msg << "[e|n|u]" << endl; */

       }
    }
    else if( mode == TimeSeries::Velocity )
    {
        ndset = 3;
       if( xyzcomponent )
       {
          dset_names[0] = "Vx";
          dset_names[1] = "Vy";
          dset_names[2] = "Vz";
          /* msg << "[xv|yv|zv]" << endl; */
       }
       else
       {
          dset_names[0] = "Vew";
          dset_names[1] = "Vns";
          dset_names[2] = "Vup";
          /* msg << "[ev|nv|uv]" << endl; */
       }
    }
    else if( mode == TimeSeries::Div )
    {
        ndset = 1;
        dset_names[0] = "Div";
    	/* msg << "[div]" << endl; */
    }
    else if( mode == TimeSeries::Curl )
    {
        ndset = 3;
    	dset_names[0] = "Curlx";
    	dset_names[1] = "Curly";
    	dset_names[2] = "Curlz";
    	/* msg << "[curlx|curly|curlz]" << endl; */
    }
    else if( mode == TimeSeries::Strains )
    {
        ndset = 6;
    	dset_names[0] = "Uxx";
    	dset_names[1] = "Uyy";
    	dset_names[2] = "Uzz";
    	dset_names[3] = "Uxy";
    	dset_names[4] = "Uxz";
    	dset_names[5] = "Uyz";
    	/* msg << "[xx|yy|zz|xy|xz|yz]" << endl; */
    }
    else if( mode == TimeSeries::DisplacementGradient )
    {
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
    	/* msg << "[duxdx|duxdy|duxdz|duydx|duydy|duydz|duzdx|duzdy|duzdz]" << endl; */
    }

    for (int i = 0; i < ndset; i++) {
      total_dims = (hsize_t)(totalSteps/TimeSeries[i]->getDownSample());
      dset_space = H5Screate_simple(1, &total_dims, NULL);
      dset       = H5Dcreate(grp, dset_names[i].c_str(), H5T_NATIVE_FLOAT, dset_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
      H5Sclose(dset_space);
      ASSERT(dset >= 0);
#ifdef USE_DSET_ATTR
      std::string incname = dset_names[i] + "CMPINC";
      std::string azname  = dset_names[i] + "CMPAZ";
      createAttr(grp, incname.c_str(), H5T_NATIVE_FLOAT, attr_space1);
      createAttr(grp, azname.c_str(), H5T_NATIVE_FLOAT, attr_space1);
#else
      createAttr(dset, "CMPINC", H5T_NATIVE_FLOAT, attr_space1);
      createAttr(dset, "CMPAZ", H5T_NATIVE_FLOAT, attr_space1);
#endif
      H5Dclose(dset);
    }
    H5Gclose(grp);
  }

  H5Pclose(dcpl);
  H5Sclose(attr_space1);
  H5Sclose(attr_space3);
  H5Fclose(fid);

  elapsed_time = MPI_Wtime() - start_time;
  printf("Created SAC HDF5 file [%s] time %e seconds\n", filename.c_str(), elapsed_time);
  fflush(stdout);
  return 1;
}

int openHDF5file(vector<TimeSeries*> & TimeSeries)
{
  hid_t fapl, *fid_ptr;
  int myRank;

  if (TimeSeries.size() == 0) 
      return 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_sec2(fapl);

  /* H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL); */
  /* H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL); */
  /* H5Pset_coll_metadata_write(fapl, false); */
  /* H5Pset_all_coll_metadata_ops(fapl, false); */

  std::string path = TimeSeries[0]->getPath();
  std::string name = TimeSeries[0]->getFileName();
  std::string filename;

  // Build the file name
  if( path != "." )
    filename = path;

  filename.append(name);
  if (name.find(".hdf5") == string::npos && name.find(".h5") == string::npos) 
    filename.append(".hdf5");
 
  fid_ptr = TimeSeries[0]->getFidPtr();
  if (fid_ptr == NULL) {
    printf("%s fid_ptr is NULL, cannot open file [%s]\n", __func__, filename.c_str());
    return -1;
  }

  *fid_ptr = H5Fopen(filename.c_str(),  H5F_ACC_RDWR, fapl);
  if (*fid_ptr < 0) {
    printf("%s Error opening file [%s]\n", __func__, filename.c_str());
    return -1;
  }

  /* if (myRank == 0) { */
  /*     printf("Rank %d: HDF5 file [%s] successfully opened: %ld\n", myRank, filename.c_str(), *fid_ptr); */
  /*     fflush(stdout); */
  /* } */

  H5Pclose(fapl);

  return 0;
}


#endif // USE_HDF5

#endif // SACHDF5_C
