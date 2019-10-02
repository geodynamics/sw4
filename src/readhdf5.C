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

#ifndef READHDF5_C
#define READHDF5_C

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
#include "EW.h"
#include "TimeSeries.h"

#ifdef USE_HDF5

#include "hdf5.h"

struct traverse_data_t {
  int myRank;
  EW* ew;
  /* bool cartCoordSet; */
  string outFileName;
  int writeEvery;
  int downSample;
  TimeSeries::receiverMode mode;
  int event;
  vector<vector<TimeSeries*>> *GlobalTimeSeries;
  float_sw4 m_global_xmax;
  float_sw4 m_global_ymax;
} traverse_data_t;

static herr_t traverse_func (hid_t loc_id, const char *grp_name, const H5L_info_t *info, void *operator_data)
{
  hid_t grp, dset;
  herr_t status;
  H5O_info_t infobuf;
  EW *a_ew;
  float data[3];
  double lon, lat, depth, x, y, z;
  bool geoCoordSet = true, topodepth = false, nsew = true;

  ASSERT(operator_data != NULL);

  struct traverse_data_t *op_data = (struct traverse_data_t *)operator_data;
  a_ew = op_data->ew;
  ASSERT(a_ew != NULL);

  status = H5Oget_info_by_name (loc_id, grp_name, &infobuf, H5P_DEFAULT);
  if (infobuf.type == H5O_TYPE_GROUP) {
    /* if (op_data->myRank == 0) */
    /*   printf ("Group: [%s] \n", grp_name); */

    // TODO: read x,y,z or ns,ew,up
    grp = H5Gopen(loc_id, grp_name, H5P_DEFAULT);
    if (grp < 0) {
      printf("Error opening group [%s]\n", grp_name);
      return -1;
    }

    if (H5Lexists(grp, "STX,STY,STZ", H5P_DEFAULT) == true) {
      nsew = false;
      geoCoordSet = false;
    }

    if (nsew) {
      // STLA,STLO,STDP
      dset = H5Dopen(grp, "STLA,STLO,STDP", H5P_DEFAULT);
      status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      H5Dclose(dset);
      lat = data[0];
      lon = data[1];
      z   = data[2];
      topodepth = true;
    }
    else {
      // X, Y, Z
      dset = H5Dopen(grp, "STX,STY,STZ", H5P_DEFAULT);
      status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      H5Dclose(dset);
      x = data[0];
      y = data[1];
      z = data[2];
      topodepth = false;
    }

    depth = z;
    if (geoCoordSet)
      a_ew->computeCartesianCoord(x, y, lon, lat);

    bool inCurvilinear=false;
    // we are in or above the curvilinear grid 
    if ( a_ew->topographyExists() && z < a_ew->m_zmin[a_ew->mNumberOfCartesianGrids-1])
      inCurvilinear = true;

    // check if (x,y,z) is not in the global bounding box
    if ( !( (inCurvilinear || z >= 0) && x>=0 && x<=op_data->m_global_xmax && y>=0 && y<=op_data->m_global_ymax)) {
      // The location of this station was outside the domain, so don't include it in the global list
      if (op_data->myRank == 0 && a_ew->getVerbosity() > 0) {
        stringstream receivererr;
    
        receivererr << endl 
  		  << "***************************************************" << endl
  		  << " WARNING:  RECEIVER positioned outside grid!" << endl;
        receivererr << " No RECEIVER file will be generated for file = " << op_data->outFileName<< endl;
        if (geoCoordSet) {
  	  receivererr << " @ lon=" << lon << " lat=" << lat << " depth=" << depth << endl << endl;
        }
        else {
    	  receivererr << " @ x=" << x << " y=" << y << " z=" << z << endl << endl;
        }
        
        receivererr << "***************************************************" << endl;
        cerr << receivererr.str();
        cerr.flush();
      }
    }
    else
    {
      /* if (op_data->myRank == 0) */
      /*   cout << "x=" << x << ", y=" << y << ", z=" << z << ", writeEvery=" << op_data->writeEvery << endl; */

      TimeSeries *ts_ptr = new TimeSeries(a_ew, op_data->outFileName, grp_name, op_data->mode, false, false, true, x, y, z, 
  					topodepth, op_data->writeEvery, op_data->downSample, !nsew, op_data->event );
      if((*op_data->GlobalTimeSeries)[op_data->event].size() == 0) 
        ts_ptr->allocFid();
      else 
        ts_ptr->setFidPtr((*op_data->GlobalTimeSeries)[op_data->event][0]->getFidPtr());
  
      // include the receiver in the global list
      (*op_data->GlobalTimeSeries)[op_data->event].push_back(ts_ptr);
    }

    H5Gclose(grp);
  }

  return 0;
}


void readStationHDF5(EW *ew, string inFileName, string outFileName, int writeEvery, int downSample, TimeSeries::receiverMode mode, int event, vector< vector<TimeSeries*> > &GlobalTimeSeries, float_sw4 m_global_xmax, float_sw4 m_global_ymax)
{
  hid_t fid, fapl;
  struct traverse_data_t tData;

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);

  fid = H5Fopen(inFileName.c_str(),  H5F_ACC_RDONLY, fapl);
  if (fid < 0) {
    printf("%s Error opening file [%s]\n", __func__, inFileName.c_str());
    return;
  }

  if (inFileName == outFileName) {
    if (ew->getRank() == 0) {
      printf("Warning: Same station input file and output file name [%s]\n", inFileName.c_str());
    }
  }

  tData.myRank = ew->getRank();
  tData.ew = ew;
  tData.outFileName = outFileName;
  tData.writeEvery = writeEvery;
  tData.downSample = downSample;
  tData.mode = mode;
  tData.event = event;
  tData.m_global_xmax = m_global_xmax;
  tData.m_global_ymax = m_global_ymax;
  tData.GlobalTimeSeries = &GlobalTimeSeries;

  H5Literate (fid, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, traverse_func, &tData);

  H5Pclose(fapl);
  H5Fclose(fid);

}


#endif // USE_HDF5
#endif // READHDF5_C
