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
#ifndef ETREE_FILE
#define ETREE_FILE

#include <string>

#ifdef ENABLE_ETREE
#include "cencalvm/query/VMQuery.h"
#include "cencalvm/storage/ErrorHandler.h"
#include "cencalvm/storage/Geometry.h"
#endif

#include "MaterialData.h"


class EW;

class EtreeFile : public MaterialData
{
 public:
#ifdef ENABLE_ETREE
   EtreeFile(EW* ew,
             const std::string& access,
             const std::string& file,
             const std::string& xfile,
             const std::string& model,
             const std::string& log,
             const std::string& query,
             const double res );

   ~EtreeFile();

   void set_material_properties(std::vector<Sarray> & rho, 
                               std::vector<Sarray> & cs,
                               std::vector<Sarray> & cp, 
                               std::vector<Sarray> & xis, 
                               std::vector<Sarray> & xip);

   void getbox( double& latmin, double& latmax, double& lonmin, double& lonmax ) const;
   void getcorners( double& latse, double& lonse, double& latsw, double& lonsw,
                    double& latne, double& lonne, double& latnw, double& lonnw ) const;
   std::string getFileName() const;
   void setupEFile();
   const  cencalvm::storage::Geometry* getGeometry() const { return mQueryGeom; }
      
  protected:

// This is a very crude test! The actual box is rotated relative to the North and East directions
   inline bool inside( double lat, double lon, double elev )
      {
         return m_latmin <= lat && lat <= m_latmax && m_lonmin <= lon && lon <= m_lonmax 
            && m_elevmin <= elev && elev <= m_elevmax;
      }

   double max( double a, double b, double c, double d );

   double min( double a, double b, double c, double d );

   void initialize(const std::string& model);

   void readEFile(std::vector<Sarray> & rho, 
  	          std::vector<Sarray> & cs,
	          std::vector<Sarray> & cp, 
	          std::vector<Sarray> & qs, 
	          std::vector<Sarray> & qp,
	          int& outside, int& material);

   std::string mAccess;
   std::string mFileName;
   std::string mXFileName;
   std::string mLogName;
   std::string mQueryType;

   static const char* mQueryKeys[];
   int mPayloadSize;
   double* mPayload;
   
   double m_EtreeRes;

   cencalvm::query::VMQuery mQuery;
   cencalvm::storage::Geometry* mQueryGeom;

   double m_latse, m_latsw, m_latne, m_latnw;
   double m_lonse, m_lonsw, m_lonne, m_lonnw;
   
   double m_latmin, m_latmax, m_lonmin, m_lonmax, m_elevmin, m_elevmax;

   EW * mEw;

#endif
};


#endif
