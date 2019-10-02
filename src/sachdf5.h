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
#ifndef _SACHDF5_H
#define _SACHDF5_H

#include <string>
#include <sstream>
#include <iostream> 

#include "TimeSeries.h"
#include "EW.h"

#ifdef USE_HDF5

#include "hdf5.h"


int openHDF5file(vector<TimeSeries*> & TimeSeries);
int createTimeSeriesHDF5File(vector<TimeSeries*> & TimeSeries, int totalSteps, float_sw4 delta);
int writeTimeSeriesHDF5File(vector<TimeSeries*> & TimeSeries, int npts, void *data);

int createAttr(hid_t loc, const char *name, hid_t type_id, hid_t space_id);
int createWriteAttr(hid_t loc, const char *name, hid_t type_id, hid_t space_id, void *data);
int openWriteAttr(hid_t loc, const char *name, hid_t type_id, void *data);
int createWriteAttrStr(hid_t loc, const char *name, const char* str);
int openWriteData(hid_t loc, const char *name, hid_t type_id, void *data, int ndim, hsize_t *start, hsize_t *count, int total_npts, float btime, float cmpinc, float cmpaz, bool isIncAzWritten, bool isLast);
int readAttrStr(hid_t loc, const char *name, char* str);
int readAttrInt(hid_t loc, const char *name, int *data);
int readAttrFloat(hid_t loc, const char *name, float *data);
int readData(hid_t loc, const char *name, int npts, void *data);

#endif // USE_HDF5

#endif // SACHDF5_H
