#ifdef ENABLE_ETREE
#include "mpi.h"
#include "EtreeFile.h"
#include "Require.h"
#include "cencalvm/storage/Payload.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "nearlyEqual.h"

#include "EW.h"

// Cencalvm currently found at:  http://earthquake.usgs.gov/research/structure/3dgeologic/download.php
using namespace std;

const char* EtreeFile::mQueryKeys[] = \
   {"Density", "Vp", "Vs", "elevation", "Qp", "Qs","FaultBlock"};

//-----------------------------------------------------------------------
EtreeFile::EtreeFile(EW* ew,
          const string& accessmode,
          const string& file, 
          const string& xfile,
          const string& model,
          const string& logfile,
          const string& query, 
          const double res):
  mEw(ew),
  mAccess(accessmode),
  mFileName(file),
  mXFileName(xfile),
  mLogName(logfile),
  mQueryType(query),
  mQueryGeom(NULL),
  m_EtreeRes(res),
  mPayloadSize(7)
{
   mPayload = new double[mPayloadSize];
   initialize(model);
// does the efile cover all points?
// we could check the vertical extent, but the horizontal extent is currently not very precise (only uses the min and max (lat,lon)). 
// Therefore, we make the conservative guess that not all points are covered.
  mCoversAllPoints = false;
}

//-----------------------------------------------------------------------
EtreeFile::~EtreeFile()
{
   delete[] mPayload;
}
//-----------------------------------------------------------------------
void EtreeFile::initialize(const string& model)
{
   //---------------------------------------------------------------
   // 1. checks validity of setting query options
   // 2. determines if the etree files that the user specifies
   //    are valid
   // 3. Sets geographic coordinates of the efile
   // 4. Sets the geometry if appropriate
   //
   // End result:  
   //   GeoCoords set:  mGeoBox
   //   Geometry set:  mEfileGeom
   //---------------------------------------------------------------

   string err = "Efile Error: ";

   //----------------------------------------------
   // User specified, if it is not there, abort
   //----------------------------------------------
   VERIFY2(access(mFileName.c_str(), R_OK) == 0,
           err << "No read permission on etree file: " << mFileName);
   if (mXFileName != "NONE")
   {
      //----------------------------------------------
      // User specified, if it is not there, abort
      //----------------------------------------------
      VERIFY2(access(mXFileName.c_str(), R_OK) == 0,
              err << "No read permission on xefile: " << mXFileName);
   }
   //----------------------------------------------------------------
   // Set coordinates
   //----------------------------------------------------------------
   if (model == "SF")
   {
      //-------------------------------------------------------
      // These are hardcoded from the published website:
      //
      // http://www.sf06simulation.org/geology/blockmodel/
      //-------------------------------------------------------
      // lat lon depth
      
      if (mXFileName != "NONE")
      {
         // constrain on the extended model
         
         m_latse = 36.702176;
	 m_latsw = 35.009018;
         m_latne = 41.484869;
         m_latnw = 39.680558;
         m_lonse = -118.94452;
         m_lonsw = -121.930857;
         m_lonne = -123.273199;
         m_lonnw = -126.353173;
	 //         GeographicCoord SE(36.702176, -118.944524, 0.0);
	 //         GeographicCoord SW(35.009018, -121.930857, 0.0);
	 //         GeographicCoord NE(41.484869, -123.273199, 45000.0);
	 //         GeographicCoord NW(39.680558, -126.353173, 45000.0); 
	 //         mGeoBox = new GeographicBox(NW,NE,SW,SE);
// guessing these values
	 m_elevmax = 3000.;
	 m_elevmin = -45000.;
      }
      else
      {
         // constrain on detailed model

         m_latse = 37.050062;
	 m_latsw = 36.320331;
         m_latne = 39.174505;
         m_latnw = 38.424179;
         m_lonse = -120.644051;
         m_lonsw = -121.922036;
         m_lonne = -122.562365;
         m_lonnw = -123.858493;
	 //         GeographicCoord SE(37.050062, -120.644051, 0.0);
	 //         GeographicCoord SW(36.320331, -121.922036, 0.0);
	 //         GeographicCoord NW(38.424179, -123.858493, 45000.0);
	 //         GeographicCoord NE(39.174505, -122.562365, 45000.0);
	 //         mGeoBox = new GeographicBox(NW, NE, SW, SE);
// guessing these values
	 m_elevmax = 3000.;
	 m_elevmin = -45000.;
      }
      mQueryGeom = 0; // it is the default
   }
   else
      VERIFY2(0,
              "Efile model must be SF not: " << model << endl);

//   int myRank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//    if (myRank == 0) // proper verbosity?
//    {
//       cout << "Efile Model Bounds: " << endl;
//       mGeoBox->printRanges(); // currently does nothing
//    }

   // Compute bounding box
   m_latmin = min( m_latse, m_latsw, m_latne, m_latnw );
   m_latmax = max( m_latse, m_latsw, m_latne, m_latnw );
   m_lonmin = min( m_lonse, m_lonsw, m_lonne, m_lonnw );
   m_lonmax = max( m_lonse, m_lonsw, m_lonne, m_lonnw );
   //   GeographicCoord c1, c2, c3, c4;
   //   mGeoBox->getBounds( c1, c2, c3, c4 );
   //   m_latmin = min( c1.getLatitude(), c2.getLatitude(), c3.getLatitude(), c4.getLatitude() );
   //   m_latmax = max( c1.getLatitude(), c2.getLatitude(), c3.getLatitude(), c4.getLatitude() );
   //   m_lonmin = min( c1.getLongitude(), c2.getLongitude(), c3.getLongitude(), c4.getLongitude() );
   //   m_lonmax = max( c1.getLongitude(), c2.getLongitude(), c3.getLongitude(), c4.getLongitude() );
}

//-----------------------------------------------------------------------
void EtreeFile::getbox( double& latmin, double& latmax, double& lonmin, double& lonmax ) const
{ 
   latmin = m_latmin;
   latmax = m_latmax;
   lonmin = m_lonmin;
   lonmax = m_lonmax;
}

//-----------------------------------------------------------------------
void EtreeFile::getcorners( double& latse, double& lonse, double& latsw, double& lonsw,
			    double& latne, double& lonne, double& latnw, double& lonnw ) const
{
   latse = m_latse;
   lonse = m_lonse;
   latsw = m_latsw;
   lonsw = m_lonsw;
   latne = m_latne;
   lonne = m_lonne;
   latnw = m_latnw;
   lonnw = m_lonnw;
}

//-----------------------------------------------------------------------
//void EtreeFile::thresholdVp(double vpmin)
//{
//   mThresholdVP = true;
//   mMinVP = vpmin;
//}

//-----------------------------------------------------------------------
//void EtreeFile::thresholdVs(double vsmin)
//{
//   mThresholdVS = true;
//   mMinVS = vsmin;
//}

//-----------------------------------------------------------------------
std::string EtreeFile::getFileName() const { return mFileName; }

//-----------------------------------------------------------------------
void EtreeFile::setupEFile()
{
// set the database file name
   mQuery.filename(mFileName.c_str());
   if( mXFileName != "NONE" )
      mQuery.filenameExt(mXFileName.c_str());

// For anything other than the SF model
   if( mQueryGeom != 0 )
       mQuery.geometry(mQueryGeom);

   if( mQueryType == "FIXEDRES" )
   {
      mQuery.queryType(cencalvm::query::VMQuery::FIXEDRES);
      mQuery.queryRes(m_EtreeRes);
   }
   else
      mQuery.queryType(cencalvm::query::VMQuery::MAXRES);

// Set values to be returned in queries
  mQuery.queryVals(mQueryKeys, mPayloadSize);
}

//-----------------------------------------------------------------------
void EtreeFile::readEFile(std::vector<Sarray> & rho, 
             std::vector<Sarray> & cs,
             std::vector<Sarray> & cp, 
             std::vector<Sarray> & qs, 
             std::vector<Sarray> & qp,
             int& outside, int& material)             
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   float NODATAVAL = cencalvm::storage::Payload::NODATAVAL;

   //   bool surfaceThresholdOn = false;
   //   double surfaceThreshold = 0.0;
   bool use_attenuation = mEw->usingAttenuation();
   int ghost_points = mEw->getNumberOfGhostPoints();

   mQuery.open();

   // ----------------------------------------------------------
   // Note:  If we run into grid points which are outside of
   //        the etree database, we will flag them as outside,
   //        but not kill the run.  This way you can safely 
   //        embed an efile inside of a block 
   // ----------------------------------------------------------

   double x, y, z, lon, lat, elev, density, vp, vs, elevation, qup, qus, faultBlock;
   double zeta, elevDelta=25.0;
   int g, k, topLevel = mEw->mNumberOfGrids-1;

// start by figuring out the max elevation where the etree returns solid material
   for (int j = mEw->m_jStart[topLevel]; j <= mEw->m_jEnd[topLevel]; ++j)
      for (int i = mEw->m_iStart[topLevel]; i <= mEw->m_iEnd[topLevel]; ++i)
      {
	 x = (i-1)*mEw->mGridSize[topLevel];
	 y = (j-1)*mEw->mGridSize[topLevel];
	 mEw->computeGeographicCoord( x, y, lon, lat); 
	   
// initial query for elevation just below sealevel
	 elev = -25.0;

	 mQuery.query(&mPayload, mPayloadSize, lon, lat, elev);
                  
// Make sure the query didn't generated a warning or error
	 if (mQuery.errorHandler()->status() != cencalvm::storage::ErrorHandler::OK) 
	 {
// If query generated an error, then bail out, otherwise reset status
	    mQuery.errorHandler()->resetStatus();
	    cout << "WARNING: Etree query failed for initial elevation of topography at grid point (i,j)= (" << i << ", " << j 
	     << ") in grid g = " << topLevel << endl
	     << " lat= " << lat << " lon= " << lon << " query elevation= " << elev << endl;
	    mEw->mTopoMat(i,j,1) = 0.;
	    continue;
	 } // if

// get the max elevation where the etree will return solid material (for assigning material properties in the curvilinear grid)
	 elev = mPayload[3];
	 int maxIter=20;
	 for (int iter=0; iter<maxIter; iter++)
	 {
	    mQuery.query(&mPayload, mPayloadSize, lon, lat, elev);
// If query generated an error, we are most likely outside
	    if (mQuery.errorHandler()->status() != cencalvm::storage::ErrorHandler::OK) 
	    {
	       mQuery.errorHandler()->resetStatus();
	    }
	    if (mPayload[1] != NODATAVAL && mPayload[2] != NODATAVAL) // bail out if we have valid Vp and Vs values
	       break;
	    elev -= elevDelta; // otherwise reduce the elevation and try again
	 }
// make sure Vp and Vs are defined
	 CHECK_INPUT (mPayload[1] != NODATAVAL && mPayload[2] != NODATAVAL, "Vp and/or Vs are not defined for grid point (i,j)= (" 
		   << i << ", " << j << ") in curvilinear grid g = " << topLevel << endl
		   << " lat= " << lat << " lon= " << lon << " query elevation= " << elev 
		   << " topo/bathy elevation= " << mPayload[3] << " Vp= " << mPayload[1] << " Vs= " << mPayload[2]);
      
	 mEw->mTopoMat(i,j,1) = elev;
      } // end for i
// end filling in mTopoMat

   if(mEw->topographyExists())
   {
      g= mEw->mNumberOfGrids-1;
      for (int j = mEw->m_jStart[g]; j <= mEw->m_jEnd[g]; ++j)
	 for (int i = mEw->m_iStart[g]; i <= mEw->m_iEnd[g]; ++i)
	 {
	    for (k = 1; k <= mEw->m_kEnd[g]; ++k) // don't attempt querying the etree above the topography (start at k=1)
	    {
	       x = mEw->mX(i,j,k);
	       y = mEw->mY(i,j,k);
// zeta is a linear function of k which should be zeta=1 for k=1 and zeta=0 for k=m_kEnd - m_ghost_points
	       zeta = ((double) (mEw->m_kEnd[g] - ghost_points - k))/(mEw->m_kEnd[g] - ghost_points - 1.);
// modify the z coordinate to account for the difference between raw and smoothed topography. 
	       z = mEw->mZ(i,j,k) + ( mEw->mTopoGrid(i,j,1) - mEw->mTopoMat(i,j,1) )*zeta;
	       mEw->computeGeographicCoord(x, y, lon, lat);
	       elev = -z;
	       if( inside( lat, lon, elev )  )
	       {
//---------------------------------------------------------
// Query the location...
//---------------------------------------------------------
		  mQuery.query(&mPayload, mPayloadSize, lon, lat, elev);
                   
		  if (mQuery.errorHandler()->status() != cencalvm::storage::ErrorHandler::OK) 
		  {
		     mQuery.errorHandler()->resetStatus(); // C.B> If we have an error, we are outside.
		     outside++; 
		     continue;
		  }
                   
		  density    = mPayload[0];
		  vp         = mPayload[1];
		  vs         = mPayload[2];
		  elevation  = mPayload[3];
		  qup        = mPayload[4];
		  qus        = mPayload[5];
		  faultBlock = mPayload[6];
	     
                   
// Vp, Vs, Density = NODATAVAL (-999) means the material properties are undefined
		  if (vp == NODATAVAL || vs == NODATAVAL || density == NODATAVAL)
		  {
		     if (k==1) // game over...
		     {
			CHECK_INPUT (vp != NODATAVAL, "Vp undefined for grid point (i,j,k)= (" 
			 << i << ", " << j << ", " << k << ") in curvilinear grid g = " << g << endl
			 << " lat= " << lat << " lon= " << lon << " elevation= " << elev << endl
			 << "Density=" << mPayload[0] << " Vp=" << mPayload[1] << " Vs=" << mPayload[2] 
			 << " Elevation=" << mPayload[3] << " Qp=" << mPayload[4] << " Qs=" << mPayload[5]
			     << " Faultblock= " << mPayload[6] );
			CHECK_INPUT (vs != NODATAVAL, "Water material properties detected for grid point (i,j,k)= (" 
			 << i << ", " << j << ", " << k << ") in curvilinear grid g = " << g << endl
			 << " lat= " << lat << " lon= " << lon << " elevation= " << elev << endl
			 << "Density=" << mPayload[0] << " Vp=" << mPayload[1] << " Vs=" << mPayload[2]
			 << " Topo Elevation=" << mPayload[3] << " Qp=" << mPayload[4] << " Qs=" << mPayload[5]
			 << " Faultblock= " << mPayload[6] 
			 << " mZ(i,j,k)= " << mEw->mZ(i,j,k) << " topoGrid= " <<  mEw->mTopoGrid(i,j,1) 
			 << " topoMat= " << mEw->mTopoMat(i,j,1) );
		   
			CHECK_INPUT (density != NODATAVAL, "Density undefined for grid point (i,j,k)= (" 
			 << i << ", " << j << ", " << k << ") in curvilinear grid g = " << g << endl
			 << " lat= " << lat << " lon= " << lon << " elevation= " << elev << endl
			 << "Density=" << mPayload[0] << " Vp=" << mPayload[1] << " Vs=" << mPayload[2] 
			 << " Elevation=" << mPayload[3] << " Qp=" << mPayload[4] << " Qs=" << mPayload[5]
			 << " Faultblock= " << mPayload[6] );
		     }
		     else // copy values from the previous grid point in k
		     {
			density = rho[g](i,j,k-1);
			vp      = cp[g](i,j,k-1);
			vs      = cs[g](i,j,k-1);
			if( use_attenuation )
			{
			   qus = qs[g](i,j,k-1);
			   qup = qp[g](i,j,k-1);
			}
		     }
		  } // end if vp or vs or density values == NODATAVAL
	    
	    
//-----------------------------------------------------
// In Material
//-----------------------------------------------------
		  material++;
		  rho[g](i,j,k) = density;
		  cp[g](i,j,k) = vp;
		  cs[g](i,j,k) = vs;
		  if( use_attenuation )
		  {
		     qs[g](i,j,k) = qus;
		     qp[g](i,j,k) = qup;
		  }
		  //		  if (mThresholdVP && cp[g](i,j,k) < mMinVP)
		  //		  {
		  //		     cp[g](i,j,k) = mMinVP;
		  //		     thresholdedVp++;
		  //		  }
		  //		  if (mThresholdVS && cs[g](i,j,k) < mMinVS)
		  //		  {
		  //		     cs[g](i,j,k) = mMinVS;
		  //		     thresholdedVs++;
		  //		  }
	       } // end if inside
	       
	    } // end for k...
	 } // end for i...
     
   } // end if topogography

   for (g = 0; g < mEw->mNumberOfCartesianGrids; g++)
   {
      int kLow=mEw->m_kStart[g];
      if (g == topLevel) // no topography, so k=1 is at the top surface
      {
	 kLow = 1;
      }
     
      for (int j = mEw->m_jStart[g]; j <= mEw->m_jEnd[g]; ++j)
	 for (int i = mEw->m_iStart[g]; i <= mEw->m_iEnd[g]; ++i)
	 {
	    for (int k = kLow; k <= mEw->m_kEnd[g]; ++k)
	    {
	       x = (i-1)*mEw->mGridSize[g];
	       y = (j-1)*mEw->mGridSize[g];  

	       if (g == topLevel)
	       {
// zeta is a linear function of k which should be zeta=1 for k=1 and zeta=0 for k=m_kEnd - m_ghost_points
		  zeta = ((double) (mEw->m_kEnd[g] - ghost_points - k))/(mEw->m_kEnd[g] - ghost_points - 1.);
// modify the z coordinate to account for the difference between the z=0 flat grid and the top elevation where material data
// is available
		  z = mEw->m_zmin[g]+(k-1)*mEw->mGridSize[g]+ ( 0.0 - mEw->mTopoMat(i,j,1) )*zeta;
	       }
	       else
	       {
		  z = mEw->m_zmin[g]+(k-1)*mEw->mGridSize[g];
	       }
	       mEw->computeGeographicCoord(x, y, lon, lat);
	       elev = -z;

	       if( inside( lat, lon, elev )  )
	       {
	    //---------------------------------------------------------
	    // Query the location...
	    //---------------------------------------------------------
		  mQuery.query(&mPayload, mPayloadSize, lon, lat, elev);
		  if (cencalvm::storage::ErrorHandler::OK != mQuery.errorHandler()->status()) 
		  {
		     mQuery.errorHandler()->resetStatus(); // C.B> If we have an error, we are outside.
		     outside++; 
		     continue;
		  }
		  density   = mPayload[0];
		  vp        = mPayload[1];
		  vs        = mPayload[2];
		  elevation = mPayload[3];
		  qup       = mPayload[4];
		  qus       = mPayload[5];
                 
// NOTE: COPIED THE LOGIC FOR USING THE PREVIOUS k-VALUE IF WATER IS DETECTED (AS ABOVE)
// Vp, Vs, Density = NODATAVAL (-999) means the material properties are undefined
		  if (vp == NODATAVAL || vs == NODATAVAL || density == NODATAVAL)
		  {
		     if (k==kLow) // game over...
		     {
			CHECK_INPUT (vp != NODATAVAL, "Vp undefined for grid point (i,j,k)= (" 
			 << i << ", " << j << ", " << k << ") in Cartesian grid g = " << g << endl
			 << " lat= " << lat << " lon= " << lon << " elevation= " << elev << endl
			 << "Density=" << mPayload[0] << " Vp=" << mPayload[1] << " Vs=" << mPayload[2] 
			 << " Elevation=" << mPayload[3] << " Qp=" << mPayload[4] << " Qs=" << mPayload[5]
			     << " Faultblock= " << mPayload[6] );
			CHECK_INPUT (vs != NODATAVAL, "Water material properties detected for grid point (i,j,k)= (" 
			 << i << ", " << j << ", " << k << ") in Cartesian grid g = " << g << endl
			 << " lat= " << lat << " lon= " << lon << " elevation= " << elev << endl
			 << "Density=" << mPayload[0] << " Vp=" << mPayload[1] << " Vs=" << mPayload[2]
			 << " Topo Elevation=" << mPayload[3] << " Qp=" << mPayload[4] << " Qs=" << mPayload[5]
			 << " Faultblock= " << mPayload[6] 
			 << " mZ(i,j,k)= " << mEw->mZ(i,j,k) << " topoGrid= " <<  mEw->mTopoGrid(i,j,1) 
			 << " topoMat= " << mEw->mTopoMat(i,j,1) );
		   
			CHECK_INPUT (density != NODATAVAL, "Density undefined for grid point (i,j,k)= (" 
			 << i << ", " << j << ", " << k << ") in Cartesian grid g = " << g << endl
			 << " lat= " << lat << " lon= " << lon << " elevation= " << elev << endl
			 << "Density=" << mPayload[0] << " Vp=" << mPayload[1] << " Vs=" << mPayload[2] 
			 << " Elevation=" << mPayload[3] << " Qp=" << mPayload[4] << " Qs=" << mPayload[5]
			 << " Faultblock= " << mPayload[6] );
		     }
		     else // copy values from the previous grid point in k
		     {
			density = rho[g](i,j,k-1);
			vp      = cp[g](i,j,k-1);
			vs      = cs[g](i,j,k-1);
			if( use_attenuation )
			{
			   qus = qs[g](i,j,k-1);
			   qup = qp[g](i,j,k-1);
			}
		     }
		  } // end if vp or vs or density values == NODATAVAL

//-----------------------------------------------------
// In Material
//-----------------------------------------------------
		  material++;
		  rho[g](i,j,k) = density;
		  cp[g](i,j,k) = vp;
		  cs[g](i,j,k) = vs;
		  if( use_attenuation )
		  {
		     qs[g](i,j,k) = qus;
		     qp[g](i,j,k) = qup;
		  }
		  //		  if (mThresholdVP && cp[g](i,j,k) < mMinVP)
		  //		  {
		  //		     cp[g](i,j,k) = mMinVP;
		  //		     thresholdedVp++;
		  //		  }
		  //		  if (mThresholdVS && cs[g](i,j,k) < mMinVS)
		  //		  {
		  //		     cs[g](i,j,k) = mMinVS;
		  //		     thresholdedVs++;
		  //		  }
	       } // end if inside
	    } // end for k...
	 } // end for i...
   } // end for g...

   if (myRank ==0)
     cout << "Done reading etree properties..." << endl;
   mQuery.close();
}

//-----------------------------------------------------------------------
void EtreeFile::set_material_properties(std::vector<Sarray> & rho, 
					std::vector<Sarray> & cs,
					std::vector<Sarray> & cp, 
					std::vector<Sarray> & qs, 
					std::vector<Sarray> & qp)
{
// some checks and extra work + read etree data by readEFile
   int outside = 0, material=0;

   setupEFile();

   int myRank, comm_size;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   MPI_Comm_size(MPI_COMM_WORLD,&comm_size);

   //Report
   if (myRank == 0)
   {
      cout << "Querying at " << mQueryType;
      if (mQueryType != "MAXRES") cout << " of: " << m_EtreeRes;
      cout << endl;
   }
   
   // If the logfile was set, set it up on the error handler
   if (mLogName != "NONE")
   {
      stringstream logbyproc;
      logbyproc << mLogName;
      logbyproc << "_" << myRank;
      cout << "\tEtree log file: " <<  logbyproc.str() << endl;
      mQuery.errorHandler()->logFilename(logbyproc.str().c_str());
      logbyproc << ".max";
      //      mMaxQuery.errorHandler()->logFilename(logbyproc.str().c_str());
   }

   // -------------------------------------------------------
   //   If access is gapps, don't copy the file, but
   //   processors take turns
   // -------------------------------------------------------
   if ( mAccess == "serial")
   {
      for (int procID = 0; procID < comm_size; ++procID)
      {
         //-------------------------------------------------------
         // Loop through the filereaders one by one and read in
         // the data.  This is done in serial.
         //-------------------------------------------------------
         if (myRank == procID)
         {
	    readEFile(rho, cs, cp, qs, qp,
		      outside, material );
         } 
         MPI_Barrier (MPI_COMM_WORLD);
      }
   }
   else if ( mAccess == "parallel")
   {
      if (myRank == 0)
         cout << "Reading efile..." << endl;

      readEFile(rho, cs, cp, qs, qp,
		outside, material );
   }
   VERIFY2(mQuery.errorHandler()->status() !=	\
           cencalvm::storage::ErrorHandler::ERROR,
           "MCartError: Efile: " << mQuery.errorHandler()->message());

   // ----------------------------------------------------------------
   // sanity checks, make sure the data read in does not violate:
   //
   //       vp >= sqrt(2)*vs
   // DUE TO ROUNDOFF WHEN CONVERTING TO mu, lambda, the ratio needs to be slightly larger than sqrt(2).
   // sqrt(2)=1.4142135...   We use 1.42 to be on the safe side
   // ----------------------------------------------------------------
   // Enforce min Vp/Vs ratio
   int warningCount = 0;
   double sqrt2=1.42;
   
   for (int g = 0; g < mEw->mNumberOfGrids; g++)
      for (int i = mEw->m_iStart[g]; i <= mEw->m_iEnd[g]; ++i)
	 for (int j = mEw->m_jStart[g]; j <= mEw->m_jEnd[g]; ++j)
	    for (int k = mEw->m_kStart[g]; k <= mEw->m_kEnd[g]; ++k)
	    {
	       if (cp[g](i,j,k) < sqrt2*cs[g](i,j,k))
               {
		  warningCount++;
		  cp[g](i,j,k) = cs[g](i,j,k)*sqrt2; // better to increase cp than decrease cs
               }
	    }

   // Summarize the warnings
   int warningSum = 0;
   MPI_Reduce(&warningCount, &warningSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

   if (myRank == 0 && warningSum > 0)
      cout << endl
           << "------------------------------------------------------" 
           << endl
           << "***** EFILE Warning: " << endl
	   << "Detected " << warningSum << " grid points where vp < sqrt(2)*vs. "
	   << "Increased vp to satisfy vp = sqrt(2)*vs in all those points."
           << endl
           << "------------------------------------------------------"
           << endl;
   

   // Summarize the data
   int materialSum, outsideSum;
   MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD);
   MPI_Reduce(&outside, &outsideSum, 1, MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD);

   if (myRank == 0)
      cout << endl 
           << "--------------------------------------------------------------\n"
           << "EFile Initialized Node Types: " << endl
           << "   Material:        " << materialSum << endl
           << endl
           << "*Outside Domain:    " << outsideSum << endl
           << endl 
           << "--------------------------------------------------------------\n"
           << endl;
}

//-----------------------------------------------------------------------
double EtreeFile::max( double a, double b, double c, double d )
{
   double x = a>b?a:b;
   double y = c>d?c:d;
   double z = x>y?x:y;
   return z;
}

//-----------------------------------------------------------------------
double EtreeFile::min( double a, double b, double c, double d )
{
   double x = a<b?a:b;
   double y = c<d?c:d;
   double z = x<y?x:y;
   return z;
}

#endif

