// -*-c++-*-
// #  WPP LICENSE
// # ----------------------------------------------------------------------
// # WPP - Wave propagation Program
// # ----------------------------------------------------------------------
// # Copyright (C) 2011, Lawrence Livermore National Security, LLC.  
// # Produced at the Lawrence Livermore National Laboratory
// # 
// # Written by:
// # 
// # Bjorn Sjogreen   (sjogreen2@llnl.gov)
// # Anders Petersson  (andersp@llnl.gov)
// # 
// # Alums:
// # Stefan Nilsson      
// # Daniel Appelo
// # Kathleen McCandless (mccandless2@llnl.gov)
// # Caroline Bono
// # 
// # CODE-227123 All rights reserved.
// # 
// # This file is part of WPP, v2.1
// # 
// # Please also read docs/GPLLICENSE.txt which contains 
// # "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License as published by
// # the Free Software Foundation; version 2, dated June 1991.
// # 
// # This program is distributed in the hope that it will be useful,
// # but WITHOUT ANY WARRANTY; without even the implied warranty of
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// # terms and conditions of the GNU General Public License for more details.
// # 
// # You should have received a copy of the GNU General Public License along with
// # this program; if not, write to the Free Software Foundation, Inc., 59 Temple
// # Place, Suite 330, Boston MA 02111-1307 USA.
// # ----------------------------------------------------------------------
#ifndef SAC_H
#define SAC_H

#include <string>
#include <vector>
#include <sys/time.h>

#include "GeographicCoord.h"
#include "boundaryConditionTypes.h"

class WPP2;

class SAC 
{

public:
   SAC();
   SAC(WPP2 *wpp,
       std::string fileName,
       std::string stationName,
       int freq, 
       int writeEvery,
       bool binaryMode = false,
       bool momentMode = false );

//    // Restart Constructor
//    SAC(bool restart,
//        const Grid* g,
//        const std::string& path,
//        std::vector<int>& intdata,
//        std::vector<float>& floatdata,
//        std::vector<std::string>& stringdata);

   void setBinaryMode(bool onoff);
  
   void set_coordinate( double x, double y, double z );

   void set_nsew( );

   void set_velocities( );

   void write_usgs_format( std::string fname );

  void initialize();
  
  //    void set_path( std::string path );
  
  void initializeEpicenter(const GeographicCoord& location,
                           double timeOffset);

void recordData(int cycle ); 
void writeFile();

void setEventDate(std::string date);
void setEventTime(std::string time);
void set_formats( int usgs, int sac );
void set_div();
void set_curl();
void set_strains();

static void initializeSystemTime(tm* tptr);
  
void set_z_relative_to_topography( bool tf ) { m_zRelativeToTopography = tf; };

private:
   void writeSACFile(int npts, char* fileName, float* data, float time, float dt,
                     char* var, float inclination, float azimuth);
   
   std::string mStationName;

   double mX, mY, mZ;
   std::string mHeaderName;
   bool mMomentMode;

   int mEventYear;
   int mEventMonth;
   int mEventDay;
   int mEventHour;
   int mEventMinute;
   float mEventSecond;
   bool mDateSet, mTimeSet;
   bool m_div, m_curl, m_strains;
  
   int mWriteEvery;

   WPP2 * mWPP;
   GeographicCoord mEpicenter;
   double mEpicenterTimeOffset; 

   bool mIOInitialized;
   bool mEpicenterInitialized;

   bool mBinaryMode;
   int mFrequency;
   int m_previous;

   std::vector<float> mRecordedUX;
   std::vector<float> mRecordedUZ;
   std::vector<float> mRecordedUY;
   std::vector<float> mRecordedUXY;
   std::vector<float> mRecordedUXZ;
   std::vector<float> mRecordedUYZ;
   
   static struct tm* mTimePtr;

   bool mRestarted;

// location relative to topography
   bool m_zRelativeToTopography;

// ignore this satation if it is above the topography
   bool mIgnore;

  /* This are the int we need to print the SAC */

  int m_iSAC   ;
  int m_jSAC   ;
  int m_kSAC   ;
  int m_gridSAC;
  bool m_writeSAC;

   bool m_xycomponent, m_velocities, m_usgsformat, m_sacformat;
   double m_calpha, m_salpha, m_thxnrm, m_thynrm, m_lat, m_lon, m_dthi;
   double m_dmx, m_dmy, m_dmz, m_d0x, m_d0y, m_d0z;
   double m_dmxy, m_dmxz, m_dmyz, m_d0xy, m_d0xz, m_d0yz;
};


#endif
