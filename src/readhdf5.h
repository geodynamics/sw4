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
#ifndef _READHDF5_H
#define _READHDF5_H

void readStationHDF5(EW *ew, string inFileName, string outFileName,
                     sw4_type writeEvery, sw4_type downSample,
                     TimeSeries::receiverMode mode, sw4_type event,
                     std::vector<std::vector<TimeSeries *> > *GlobalTimeSeries,
                     float_sw4 m_global_xmax, float_sw4 m_global_ymax,
                     bool is_obs, bool winlset, bool winrset, float_sw4 winl,
                     float_sw4 winr, bool usex, bool usey, bool usez,
                     float_sw4 t0, bool scalefactor_set, float_sw4 scalefactor);

void readRuptureHDF5(char *fname,
                     std::vector<std::vector<Source *> > &a_GlobalUniqueSource,
                     EW *ew, sw4_type event, float_sw4 m_global_xmax,
                     float_sw4 m_global_ymax, float_sw4 m_global_zmax,
                     float_sw4 mGeoAz, float_sw4 xmin, float_sw4 ymin,
                     float_sw4 zmin, sw4_type mVerboses, sw4_type nreader);

void readStationInfoHDF5(string inFileName, std::vector<string> *staname,
                         std::vector<double> *x, std::vector<double> *y,
                         std::vector<double> *z, std::vector<sw4_type> *is_nsew,
                         sw4_type *n);
#endif  // _READHDF5_H
