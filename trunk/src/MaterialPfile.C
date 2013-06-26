// -*-c++-*-

#include "Require.h"

#include <cstring>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "EW.h"

#include "MaterialPfile.h"

using namespace std;

//-----------------------------------------------------------------------
MaterialPfile::MaterialPfile( EW * a_ew,
			     const std::string file,
			     const std::string directory,
			     const int nstenc,
			     const double vpmin,
			     const double vsmin,
			     const double rhomin,
			     const bool flatten,
			     const bool coords_geographic )
   :
   mEW(a_ew),
   m_model_file(file),
   m_model_dir(directory),
   m_nstenc(nstenc),
   m_vpmin(vpmin),
   m_vsmin(vsmin),
   m_rhomin(rhomin),
   m_flatten(flatten),
   m_coords_geographic(coords_geographic)
{
   read_pfile();

   mCoversAllPoints = false;
   if( m_coords_geographic )
   {
      double x1, x2, x3, x4, y1, y2, y3, y4;
      mEW->computeCartesianCoord( x1, y1, m_lonmin, m_latmin );
      mEW->computeCartesianCoord( x2, y2, m_lonmin, m_latmax );
      mEW->computeCartesianCoord( x3, y3, m_lonmax, m_latmin );
      mEW->computeCartesianCoord( x4, y4, m_lonmax, m_latmax );
      m_xmin = min(x1,min(x2,min(x3,x4)));
      m_xmax = max(x1,max(x2,max(x3,x4)));
      m_ymin = min(y1,min(y2,min(y3,y4)));
      m_ymax = max(y1,max(y2,max(y3,y4)));      
   }
   double bbox[6];
   mEW->getGlobalBoundingBox( bbox );

   if (m_xmin > bbox[0] || m_ymin > bbox[2] || m_depthmin > bbox[4] ||
       m_xmax < bbox[1] || m_ymax < bbox[3] || m_depthmax < bbox[5])
   {
      mCoversAllPoints = false;
//      cout << "This block does NOT cover all grid points" << endl;
   }
   else
   {
      mCoversAllPoints=true;
// tmp
//     cout << "This block COVERS all grid points" << endl;
   }
}

//-----------------------------------------------------------------------
void MaterialPfile::set_material_properties( std::vector<Sarray> & rho, 
					     std::vector<Sarray> & cs,
					     std::vector<Sarray> & cp, 
					     std::vector<Sarray> & qs, 
					     std::vector<Sarray> & qp  )
{
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
// the qs and qp arrays are always allocated to allow qs[g].is_defined() to be called
//  bool use_attenuation = !(qs.size()==0);
//  bool use_attenuation = m_ew->usingAttenuation();
  bool use_attenuation = false;

  int outside = 0; int material = 0;

  if (myRank == 0) cout << "Assigning material properties from pfile data..." << endl;

//tmp
//  printf("MPF: set_mat_prop: use_attenuation=%i\n", mWPP->usingAttenuation());
  
  double x, y, z, lon, lat, depth;
  double vp,vs,density,qup,qus,zsed,zmoho;
  bool foundcrust;
  
  int g, kLow, topLevel = mEW->mNumberOfGrids-1;
// first deal with the Cartesian grids
  for (g = 0; g < mEW->mNumberOfCartesianGrids; g++)
  {
    kLow=mEW->m_kStart[g];
    if (g == topLevel) // no topography, so k=1 is at the top surface
    {
      kLow = 1;
    }
    if( m_coords_geographic )
    {
//       for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k)
       for (int k = kLow; k <= mEW->m_kEnd[g]; ++k)
	  for (int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; ++j)
	     for (int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g]; ++i)
	     {
		x = (i-1)*mEW->mGridSize[g];
		y = (j-1)*mEW->mGridSize[g];
		z = mEW->m_zmin[g]+(k-1)*mEW->mGridSize[g];
                  
		mEW->computeGeographicCoord( x, y, lon, lat );
		mEW->getDepth(x,y,z,depth);

		if( inside( lat, lon, depth )  )
		{
		   //---------------------------------------------------------
		   // Query the location...
		   //---------------------------------------------------------
                   sample_latlon( lat, lon, depth, vp, vs, density, qup, qus, zsed, zmoho, foundcrust );
		   rho[g](i,j,k) = density;
		   cp[g](i,j,k) = vp;
		   cs[g](i,j,k) = vs;
		   if (m_qf) 
		   {
		      if (qp[g].is_defined())
			 qp[g](i,j,k) = qup;
		      if (qs[g].is_defined())
			 qs[g](i,j,k) = qus;
		   }
		   material++;
		}
		else
		{
		  if (mEW->getVerbosity() > 2)
		  {
		    printf("Point (i,j,k)=(%i, %i, %i) in grid g=%i\n"
			   "with (x,y,z)=(%e,%e,%e) and lat=%e, lon=%e, depth=%e\n"
			   "is outside the pfile domain: %e<= lat <= %e, %e <= lon <= %e, %e <= depth <= %e\n", 
			   i, j, k, g, 
			   x, y, z, lat, lon, depth, 
			   m_latmin, m_latmax, m_lonmin, m_lonmax, m_depthmin, m_depthmax);
		  }
		  outside++;
		}
		
	     } // end for i, j, k
    }
    else
    {
       // Cartesian p-file
       for (int k = kLow; k <= mEW->m_kEnd[g]; ++k)
	  for (int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; ++j)
	     for (int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g]; ++i)
	     {
		x = (i-1)*mEW->mGridSize[g];
		y = (j-1)*mEW->mGridSize[g];
		z = mEW->m_zmin[g]+(k-1)*mEW->mGridSize[g];

		mEW->getDepth(x,y,z,depth);

		if( inside_cart( x, y, depth )  )
		{
		   //---------------------------------------------------------
		   // Query the location...
		   //---------------------------------------------------------
		   sample_cart( x, y, z, vp, vs, density, qup, qus, zsed, zmoho, foundcrust );

		   rho[g](i,j,k) = density;
		   cp[g](i,j,k) = vp;
		   cs[g](i,j,k) = vs;

		   //		   cout << "x= " << x << " y= " << y << " depth= " << depth << " vp = " << vp << " vs = " << vs << " rho = " << density << endl;

		   if(m_qf) 
		   {
		      if (qp[g].is_defined())
			 qp[g](i,j,k) = qup;
		      if (qs[g].is_defined())
			 qs[g](i,j,k) = qus;
		   }
		   material++;
		}
		else
		   outside++;
	     } // end for i, j, k

    } // end cartesian pfile case
    
// communicate material properties to ghost points (necessary on refined meshes because ghost points don't have a well defined depth/topography)
    mEW->communicate_array( rho[g], g );
    mEW->communicate_array( cs[g], g );
    mEW->communicate_array( cp[g], g );
    
  } // end for all Cartesian grids
  

// Now, the curvilinear grid
  if (mEW->topographyExists())
  {
    g = mEW->mNumberOfGrids-1;

// the curvilinear grid is always on top
    kLow = 1;

    if( m_coords_geographic )
    {
//       for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) 
       for (int k = kLow; k <= mEW->m_kEnd[g]; ++k) // don't attempt querying the pfile above the topography (start at k=1)
	  for (int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; ++j)
	     for (int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g]; ++i)
	     {
		x = mEW->mX(i,j,k);
		y = mEW->mY(i,j,k);
		z = mEW->mZ(i,j,k);

		mEW->computeGeographicCoord( x, y, lon, lat );
		depth = z - mEW->mZ(i,j,1);

		if( inside( lat, lon, depth )  )
		{
		   //---------------------------------------------------------
		   // Query the location...
		   //---------------------------------------------------------
		   sample_latlon( lat, lon, depth, vp, vs, density, qup, qus, zsed, zmoho, foundcrust );
		   rho[g](i,j,k) = density;
		   cp[g](i,j,k) = vp;
		   cs[g](i,j,k) = vs;
		   if (m_qf) 
		   {
		      if (qp[g].is_defined())
			 qp[g](i,j,k) = qup;
		      if (qs[g].is_defined())
			 qs[g](i,j,k) = qus;
		   }
		   material++;
		}
		else
		{
		  if (mEW->getVerbosity() > 2)
		  {
		    printf("Point (i,j,k)=(%i, %i, %i) in grid g=%i\n"
			   "with (x,y,z)=(%e,%e,%e) and lat=%e, lon=%e, depth=%e\n"
			   "is outside the pfile domain: %e<= lat <= %e, %e <= lon <= %e, %e <= depth <= %e\n", 
			   i, j, k, g, 
			   x, y, z, lat, lon, depth, 
			   m_latmin, m_latmax, m_lonmin, m_lonmax, m_depthmin, m_depthmax);
		  }
		   outside++;
		}
		
	     } // end for i,j,k
    }
    else
    {
//       for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) // don't attempt querying the pfile above the topography (start at k=1)
       for (int k = kLow; k <= mEW->m_kEnd[g]; ++k) // don't attempt querying the pfile above the topography (start at k=1)
	  for (int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; ++j)
	     for (int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g]; ++i)
	     {
		x = mEW->mX(i,j,k);
		y = mEW->mY(i,j,k);
		z = mEW->mZ(i,j,k);
		depth = z - mEW->mZ(i,j,1);
		if( inside_cart( x, y, depth )  ) 
		{
		   //---------------------------------------------------------
		   // Query the location...
		   //---------------------------------------------------------
		   sample_cart( x, y, z, vp, vs, density, qup, qus, zsed, zmoho, foundcrust );
		   rho[g](i,j,k) = density;
		   cp[g](i,j,k) = vp;
		   cs[g](i,j,k) = vs;
		   if (m_qf) 
		   {
		      if (qp[g].is_defined())
			 qp[g](i,j,k) = qup;
		      if (qs[g].is_defined())
			 qs[g](i,j,k) = qus;
		   }
		   material++;
		}
		else
		   outside++;
	     } // end for i,j,k
    }
  } // end if topographyExists()
  
//  extrapolation is now done in WPP2:set_materials()

  int outsideSum, materialSum;
  MPI_Reduce(&outside, &outsideSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

  if (mEW->proc_zero()) 
    cout << "outside = " << outsideSum << ", " << "material = " << materialSum << endl;

}

//-----------------------------------------------------------------------
void MaterialPfile::read_pfile( )
{
   //   int kk, m;

   //   m_qf = false;

   //   m_model_dir  = ppdir;
   //   m_model_file = ppfile;
   string ppmfile = m_model_dir + "/" + m_model_file;

   //   m_vpmin   = vpmin;
   //   m_vsmin   = vsmin;
   //   m_rhomin  = rhomin;

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

// Open file
   FILE* fd=fopen(ppmfile.c_str(), "r" );
   if( fd == NULL )
   {
     if (myRank == 0) cerr << "Unable to open the pfile input file: '" << ppmfile << "'" << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   int bufsize = 1024;
   char* buf = new char[bufsize];
   
   // Read pfile header

   // Line 1
   CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error line 1 in pfile header not found\n");
   string tok0 = strtok( buf, " \t" );
// strip off any white space
   size_t nWhite = tok0.find_first_of(" \t\n");
   m_model_name = tok0.substr(0,nWhite);

   // Line 2
   CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error line 2 in pfile header not found\n");
   char* tok=strtok(buf," \t");
   CHECK_INPUT( tok != NULL, "Error on line 2 in pfile header, no grid spacing\n");
   m_h = atof(tok);

   // Line 3
   CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error line 3 in pfile header not found\n");
   tok = strtok(buf," \t");
   CHECK_INPUT( tok != NULL, "Error on line 3 in pfile header, no grid spacing\n");
   m_nlat = atoi(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 3 in pfile header, no min value\n");
   m_latmin = atof(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 3 in pfile header, no max value\n");
   m_latmax = atof(tok);

   // Line 4
   CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error line 4 in pfile header not found\n");
   tok = strtok(buf," \t");
   CHECK_INPUT( tok != NULL, "Error on line 4 in pfile header, no grid spacing\n");
   m_nlon = atoi(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 4 in pfile header, no min value\n");
   m_lonmin = atof(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 4 in pfile header, no max value\n");
   m_lonmax = atof(tok);

   // Line 5
   CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error line 5 in pfile header not found\n");
   tok = strtok(buf," \t");
   CHECK_INPUT( tok != NULL, "Error on line 5 in pfile header, no grid spacing\n");
   m_nmaxdepth = atoi(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 5 in pfile header, no min value\n");
   m_depthmin = atof(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 5 in pfile header, no max value\n");
   m_depthmax = atof(tok);

   // Line 6
   CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error line 6 in pfile header not found\n");
   tok = strtok(buf," \t");
   CHECK_INPUT( tok != NULL, "Error on line 6 in pfile header, no ksed\n");
   m_ksed = atoi(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 6 in pfile header, no kmoho\n");
   m_kmoho = atoi(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 6 in pfile header, no k410\n");
   m_k410 = atoi(tok);
   tok = strtok(NULL," \t");
   CHECK_INPUT( tok != NULL, "Error on line 6 in pfile header, no k660\n");
   m_k660 = atoi(tok);

   // Line 7
   CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error line 7 in pfile header not found\n");
   CHECK_INPUT( tok != NULL, "Error on line 7 in pfile header, no Q-available flag\n");      
   tok = strtok(buf," \t");
   CHECK_INPUT( tok != NULL, "Error on line 7 in pfile header, no Q-available flag\n");
   string cqf0 = tok;
// strip off any white space
   nWhite = cqf0.find_first_of(" \t\n");
   string cqf = cqf0.substr(0,nWhite);
   
// test
//   printf("Q-flag string '%s'\n", cqf.c_str());
   
   m_qf = ( cqf == "T") || (cqf == "t") || (cqf == ".TRUE.") || (cqf == ".true.");

   double km = 1000;
   if( m_coords_geographic )
   {
      m_depthmin = m_depthmin*km;
      m_depthmax = m_depthmax*km;
   }
   else
   {
      m_nx = m_nlat;
      m_ny = m_nlon;
      m_xmin = m_latmin;
      m_xmax = m_latmax;
      m_ymin = m_lonmin;
      m_ymax = m_lonmax;
   }

   if (myRank == 0)
   {
     cout << "Pfile model name (string): '" << m_model_name << "'" << endl;

     if( !m_coords_geographic )
     {
        cout << "Step size in x and y: " << m_h << endl;
	cout << "Number of x-direction points: " << m_nx << endl;
	cout << "Min x: " << m_xmin << " Max x: " << m_xmax << endl;
	cout << "Number of y-direction points: " << m_ny << endl;
	cout << "Min y: " << m_ymin << " Max y: " << m_ymax << endl;
     }
     else
     {
	cout << "Step size in lat and lon: " << m_h << endl;
	cout << "Number of latitude points: " << m_nlat << endl;
	cout << "Min Lat: " << m_latmin << " Max Lat: " << m_latmax << endl;
	cout << "Number of longitude points: " << m_nlon << endl;
	cout << "Min Lon: " << m_lonmin << " Max Lon: " << m_lonmax << endl;
     }
     cout << "Number of depth points: " << m_nmaxdepth << endl;
     cout << "Min depth: " << m_depthmin << " Max depth: " << m_depthmax << endl;
     cout << "Optional indices: Sediment: " << m_ksed << " MoHo: " << m_kmoho << " 410: " << m_k410 << " 660: " << m_k660 << endl;
     cout << "Attenuation Q-factors available: " << (m_qf? "yes":"no") << endl;
   }
   

   // Allocate arrays
   if( m_coords_geographic )
   {
      m_lat = new double[m_nlat];
      m_lon = new double[m_nlon];
   }
   else
   {
      m_x = new double[m_nx];
      m_y = new double[m_ny];
   }

   m_z   = new double[m_nlon*m_nlat*m_nmaxdepth];
   m_vp  = new double[m_nlon*m_nlat*m_nmaxdepth];
   m_vs  = new double[m_nlon*m_nlat*m_nmaxdepth];
   m_rho = new double[m_nlon*m_nlat*m_nmaxdepth];

   if(m_qf)
   {
      m_qp = new double[m_nlon*m_nlat*m_nmaxdepth];
      m_qs = new double[m_nlon*m_nlat*m_nmaxdepth];
   }
   else
   {
     if (myRank == 0) printf("ppmod: NOT allocating arrays for Qp and Qs\n");
   }
   
   m_st = new double[m_nlon*m_nlat];
   m_ct = new double[m_nlon*m_nlat];


   int m=0, kk, ndepth, line=7;

   
   if( !m_coords_geographic ) // cartesian
   {
      for(int j=0; j < m_ny; j++ )
	 for(int i=0; i < m_nx; i++ )
	 {
            CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error in pfile profile header at coordinate " 
			                              << i << " " << j << "\n" );
	    tok = strtok(buf," \t");
	    CHECK_INPUT( tok != NULL, "Error in pfile profile header, no x variable at " << i << " " << j << "\n");
            m_x[i] = atof(tok);

	    tok = strtok(NULL," \t");
	    CHECK_INPUT( tok != NULL, "Error in pfile profile header, no y variable at " << i << " " << j << "\n");
	    m_y[j] = atof(tok);

	    tok = strtok(NULL," \t");
	    CHECK_INPUT( tok != NULL, "Error in pfile profile header, no ndepth  at " << i << " " << j << "\n");
	    ndepth = atoi(tok);

	    line++;

// sanity check
	    double y = m_ymin + j*m_h;
	    double x = m_xmin + i*m_h;
	    if (fabs(y - m_y[j]) + fabs(x - m_x[i]) > 0.1*m_h)
	    {
	       if (myRank == 0)
	       {
		  cerr << "pfile reader error, ppmod file line=" << line << endl;
		  cerr << "read x[" << i << "]=" << m_x[i] << " but expected x=" << x << endl;
		  cerr << "read y[" << j << "]=" << m_y[j] << " but expected y=" << y << endl;
		  cerr << "CHECK THE PPMOD FILE." << endl;
		  cerr << "DEPTH PROFILES SHOULD BE ORDERED SUCH THAT X VARIES THE FASTEST" << endl;
	       }
	       MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	 
// sanity check 2
	    if (ndepth != m_nmaxdepth )
	    {
	       if (myRank == 0)
	       {
		  cerr << "pfile reader error, ppmod file line=" << line << endl;
		  cerr << "read ndepth=" << ndepth << " which is different from header nmaxdepth="
		       << m_nmaxdepth << endl;
	       }
	       MPI_Abort(MPI_COMM_WORLD, 1);
	    }

	    // Read depth profile       
	    for(int k=0; k < m_nmaxdepth; k++ )
	    {
	       CHECK_INPUT( fgets( buf, bufsize, fd ) != NULL, "Error in pfile profile at coordinate " 
			    << i << " " << j << " " << k << "\n" );

	       tok = strtok(buf," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading kk at " << i << " " << j << " " << k << "\n" );
	       kk = atoi( tok );

	       tok = strtok(NULL," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading z at " << i << " " << j << " " << k << "\n" );
	       m_z[m] = atof( tok );

	       tok = strtok(NULL," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading vp at " << i << " " << j << " " << k << "\n" );
	       m_vp[m] = atof( tok );

	       tok = strtok(NULL," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading vs at " << i << " " << j << " " << k << "\n" );
	       m_vs[m] = atof( tok );

	       tok = strtok(NULL," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading rho at " << i << " " << j << " " << k << "\n" );
	       m_rho[m] = atof( tok );

	       if( m_qf )
	       {
		  tok = strtok(NULL," \t");
		  CHECK_INPUT( tok != NULL, "Error in pfile reading qp at " << i << " " << j << " " << k << "\n" );
		  m_qp[m] = atof( tok );

		  tok = strtok(NULL," \t");
		  CHECK_INPUT( tok != NULL, "Error in pfile reading qs at " << i << " " << j << " " << k << "\n" );
		  m_qs[m] = atof( tok );
	       }
	       line++;
	       m_vp[m] = max(m_vp[m], m_vpmin );
	       m_vs[m] = max(m_vs[m], m_vsmin );
	       m_rho[m]= max(m_rho[m],m_rhomin);

	       if ( k == m_ksed-1 ) 
		  m_st[i+j*m_nx] = m_z[m];

	       if ( k == m_kmoho-1 ) 
		  m_ct[i+j*m_nx] = m_z[m];

// fundamental sanity checks
	       if (!(m_y[j] <= m_ymax && m_y[j] >= m_ymin && 
		     m_x[i] <= m_xmax && m_x[i] >= m_xmin) )
	       {
		  printf("Error reading pfile: x profile #%i, y profile #%i: x=%e or y=%e out of bounds!"
			 " min(x)=%e, max(x)=%e, min(y)=%e, max(y)=%e\n", i+1, j+1,
			 m_x[i], m_y[j], m_xmin, m_xmax, m_ymin, m_ymax );
		  MPI_Abort(MPI_COMM_WORLD, 1);
	       }
	       m++;
	    }
	 }
   }
   else // geographic coordinates
   {
      for(int j=0; j < m_nlat; j++ )
	 for(int i=0; i< m_nlon; i++ )
	 {
            CHECK_INPUT( fgets(buf,bufsize,fd) != NULL, "Error in pfile profile header at coordinate " 
			                              << i << " " << j << "\n" );
	    tok = strtok(buf," \t");
	    CHECK_INPUT( tok != NULL, "Error in pfile profile header, no lat. variable at " << i << " " << j << "\n");
            m_lat[j] = atof(tok);

	    tok = strtok(NULL," \t");
	    CHECK_INPUT( tok != NULL, "Error in pfile profile header, no lon. variable at " << i << " " << j << "\n");
	    m_lon[i] = atof(tok);

	    tok = strtok(NULL," \t");
	    CHECK_INPUT( tok != NULL, "Error in pfile profile header, no ndepth  at " << i << " " << j << "\n");
	    ndepth = atoi(tok);

	    line++;

// sanity check (have to relax this to allow for different step sizes in lat and lon)
	    double lat_j = m_latmin + j*m_h;
	    double lon_i = m_lonmin + i*m_h;
	    if (fabs(lat_j - m_lat[j]) + fabs(lon_i - m_lon[i]) > 0.1*m_h)
	    {
	       if (myRank == 0)
	       {
		  cerr << "pfile reader error, ppmod file line=" << line << endl;
		  cerr << "read lon_ppm[" << i << "]=" << m_lon[i] << " but expected lon=" << lon_i << endl;
		  cerr << "read lat_ppm[" << j << "]=" << m_lat[j] << " but expected lat=" << lat_j << endl;
		  cerr << "CHECK THE PPMOD FILE." << endl;
		  cerr << "DEPTH PROFILES SHOULD BE ORDERED SUCH THAT LONGITUDE VARIES THE FASTEST" << endl;
	       }
	       MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	 
// sanity check 2
	    if (ndepth != m_nmaxdepth )
	    {
	       if (myRank == 0)
	       {
		  cerr << "pfile reader error, ppmod file line=" << line << endl;
		  cerr << "read ndepth=" << ndepth << " which is different from header nmaxdepth="
		       << m_nmaxdepth << endl;
	       }
	       MPI_Abort(MPI_COMM_WORLD, 1);
	    }

	    // Read depth profile       
	    for(int k=0; k < m_nmaxdepth; k++ )
	    {
	       CHECK_INPUT( fgets( buf, bufsize, fd ) != NULL, "Error in pfile profile at coordinate " 
			    << i << " " << j << " " << k << "\n" );

	       tok = strtok(buf," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading kk at " << i << " " << j << " " << k << "\n" );
	       kk = atoi( tok );

	       tok = strtok(NULL," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading z at " << i << " " << j << " " << k << "\n" );
	       m_z[m] = km*atof( tok );

	       tok = strtok(NULL," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading vp at " << i << " " << j << " " << k << "\n" );
	       m_vp[m] = km*atof( tok );

	       tok = strtok(NULL," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading vs at " << i << " " << j << " " << k << "\n" );
	       m_vs[m] = km*atof( tok );

	       tok = strtok(NULL," \t");
	       CHECK_INPUT( tok != NULL, "Error in pfile reading rho at " << i << " " << j << " " << k << "\n" );
	       m_rho[m] = km*atof( tok );

	       if( m_qf )
	       {
		  tok = strtok(NULL," \t");
		  CHECK_INPUT( tok != NULL, "Error in pfile reading qp at " << i << " " << j << " " << k << "\n" );
		  m_qp[m] = atof( tok );

		  tok = strtok(NULL," \t");
		  CHECK_INPUT( tok != NULL, "Error in pfile reading qs at " << i << " " << j << " " << k << "\n" );
		  m_qs[m] = atof( tok );
	       }
	       line++;
	       m_vp[m] = max(m_vp[m], m_vpmin );
	       m_vs[m] = max(m_vs[m], m_vsmin );
	       m_rho[m]= max(m_rho[m],m_rhomin);

	       if ( k == m_ksed-1 ) 
		  m_st[i+j*m_nlon] = m_z[m];

	       if ( k == m_kmoho-1 ) 
		  m_ct[i+j*m_nlon] = m_z[m];

// fundamental sanity checks
	       if (!(m_lat[j] <= m_latmax && m_lat[j] >= m_latmin && 
		     m_lon[i] <= m_lonmax && m_lon[i] >= m_lonmin) )
	       {
		  printf("Error reading pfile: lat profile #%i, lon profile #%i: lat=%e or lon=%e out of bounds!"
			 " min(lat)=%e, max(lat)=%e, min(lon)=%e, max(lon)=%e\n", j+1, i+1,
			 m_lat[j], m_lon[i], m_latmin, m_latmax, m_lonmin, m_lonmax );
		  MPI_Abort(MPI_COMM_WORLD, 1);
	       }
	       m++;
	    }
	 }
   }
   delete[] buf;
}

//---------------------------------------------------------------------------------------
//int MaterialPfile::get_material_pt( double x, double y, double z, double& rho, double& cs, double& cp,
//				    double& qs, double& qp )
//  {
//   int retval = 0;
//   double zsed, zmoho;
//   bool foundcrust;
   
//   if( m_coords_geographic )
//   {
//     double lon, lat, depth;
     
//     mEW->computeGeographicCoord( x, y, lon, lat );
//     mEW->getDepth(x,y,z,depth);

//     if( inside( lat, lon, depth )  )
//     {
       //---------------------------------------------------------
       // Query the location...
       //---------------------------------------------------------
//       sample_latlon( lat, lon, depth, cp, cs, rho, qp, qs, zsed, zmoho, foundcrust );
//     }
//     else
//       retval = -1;
//   }
//   else
//   {
//     if( inside_cart( x, y, z )  ) // elev = -depth
//     {
       //---------------------------------------------------------
       // Query the location...
       //---------------------------------------------------------
//       sample_cart( x, y, z, cp, cs, rho, qp, qs, zsed, zmoho, foundcrust );
//     }
//     else
//       retval = -1;
//   }
//   return retval;
//  }


//-----------------------------------------------------------------------
void MaterialPfile::sample_latlon( double lats,double lons,double zs, double &vp, 
				    double &vs,double &rho, double &qp, double &qs,
				    double &zsed, double &zmoho, bool& foundcrust )
//--------------------------------------------------------------------------
// return material properties (vp, vs, rho) at point (lats, lons, zs)
//--------------------------------------------------------------------------
{
   //   double w;
   //   int i, j, k;
   //   int kk;
   //   int i1, j1, k1;
   //   int i1max, j1max;
   //   double  factor;
   //   double wt;
   //   int m;

   foundcrust=false;

// tmp
//   printf("ppmod::sample: lats=%e, lons=%e, zs=%e\n", lats, lons, zs);

   //  Check if lats and lons are out of range
   if ( lats < m_latmin )
   {
     cerr << "MaterialPfile::sample lats out of range (min): " <<  lats << ", " <<  m_latmin << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   if ( lats > m_latmax )
   {
     cerr << "MaterialPfile::sample lats out of range (max): "<< lats<< ", " << m_latmax << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
     
   if ( lons < m_lonmin)
   {
     cerr << "MaterialPfile::sample lons out of range (min): " << lons << ", " <<  m_lonmin << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
     
   if (lons > m_lonmax)
   {
     cerr << "MaterialPfile::sample lons out of range (max): " << lons << ", " << m_lonmax << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }

   int ii, jj;
   if( m_nstenc % 2 == 1 )
   {
      // odd number of points in stencil
      int s = (m_nstenc-1)/2;
      ii  = static_cast<int>( floor( (lons-m_lonmin+0.5*m_h)/m_h ) )-s;
      jj  = static_cast<int>( floor( (lats-m_latmin+0.5*m_h)/m_h ) )-s;
   }
   else
   {
      // even number of points in stencil
      int s = m_nstenc/2;
      ii  = static_cast<int>( floor( (lons-m_lonmin)/m_h ) )-s+1;
      jj  = static_cast<int>( floor( (lats-m_latmin)/m_h ) )-s+1;
   }

// make sure we stay within array boundaries
   if( ii < 0 )
      ii = 0;
   if( jj < 0 )
      jj = 0;

   int ii2 = ii + m_nstenc-1;
   int jj2 = jj + m_nstenc-1;

   if( ii2 >= m_nlon )
   {
     ii2 = m_nlon-1;
     ii = ii2 -(m_nstenc-1);
   }
   
   if( jj2 >= m_nlat )
   {
     jj2 = m_nlat-1;
     jj = jj2 - (m_nstenc-1);
   }
   
   double re=6371;
   if( m_flatten ) 
      zs = re*(1.0 - exp(-zs/re));


   double w=0;
   zsed=zmoho=vp=vs=rho=qp=qs=0; 
   double appm  = 0.5*m_nstenc*m_h/sqrt(-log(1e-6));
   double appmi2 = 1.0/(appm*appm);

   for( int j1 = 0 ; j1 < jj2-jj+1 ; j1++ )
      for( int i1 = 0 ; i1 < ii2-ii+1 ; i1++ )
      {
	  double wgh = exp(-( (lons-m_lon[ii+i1])*(lons-m_lon[ii+i1])
			     +(lats-m_lat[jj+j1])*(lats-m_lat[jj+j1]) )*appmi2 );
	  w += wgh;

	  int m = (ii+i1)*m_nmaxdepth + (jj+j1)*m_nmaxdepth*m_nlon;

          int kk;
	  for( kk=0; kk < m_nmaxdepth; kk++ )
	  {
	     if (m_z[m+kk] > zs) break;
	  }

	  if (kk+1 >= m_kmoho)
	     foundcrust = true;

	  int k1 = kk-1;
	  double factor = (zs-m_z[m+k1])/(m_z[m+k1+1]-m_z[m+k1]);

	  zsed  += m_z[m+m_ksed-1]*wgh;
	  zmoho += m_z[m+m_kmoho-1]*wgh;
	  vp    += (m_vp[m+k1]  + factor*(m_vp[m+k1+1]-m_vp[m+k1])  )*wgh;
	  vs    += (m_vs[m+k1]  + factor*(m_vs[m+k1+1]-m_vs[m+k1])  )*wgh;
	  rho   += (m_rho[m+k1] + factor*(m_rho[m+k1+1]-m_rho[m+k1]))*wgh;
	  if( m_qf )
	  {
	     qp += (m_qp[m+k1] + factor*(m_qp[m+k1+1]-m_qp[m+k1]))*wgh;
	     qs += (m_qs[m+k1] + factor*(m_qs[m+k1+1]-m_qs[m+k1]))*wgh;
	  }
      }

   // Now compute average properties by distance-weighted Gaussian average
   double iw;
   if (w != 0.)
     iw = 1.0/w;
   else
   {
     printf("Error MaterialPfile::sample_latlon: weight w = 0 at lat=%e, lon=%e, depth=%e\n", lats, lons, zs);
// tmp
     printf("ii=%i, ii2=%i, jj=%i, jj2=%i, lon_ppm[ii]=%e, lat_ppm[jj]=%e\n", ii, ii2, jj, jj2, m_lon[ii], m_lat[jj]);
     double dist2    = (lons-m_lon[ii])*(lons-m_lon[ii]) + (lats-m_lat[jj])*(lats-m_lat[jj]);
     double exponent = -dist2*appmi2;
     double wgh      =  exp( exponent );
     printf("dist2=%e, appm=%e, appmi2=%e, exponent=%e, wgh=%e\n", dist2, appm, appmi2, exponent, wgh);

     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   
   zsed  *= iw;
   zmoho *= iw;
   vp  *= iw;
   vs  *= iw;
   rho *= iw;
   if( m_qf ) 
   {
      qp *= iw;
      qs *= iw;
   }
   if( m_flatten )
   {
      zs  = re*log(re/(re-zs));
      vp  = vp*(re/(re-zs));
      vs  = vs*(re/(re-zs));
      rho = rho*pow(re/(re-zs),5.0);
   }
}

//-----------------------------------------------------------------------
void MaterialPfile::sample_cart( double xs,double ys,double zs, double &vp, 
				  double &vs,double &rho, double &qp, double &qs,
				  double &zsed, double &zmoho, bool& foundcrust )
//--------------------------------------------------------------------------
// return material properties (vp, vs, rho) at point (xs, ys, zs)
//--------------------------------------------------------------------------
{
   foundcrust=false;
   //  Check if xs and ys are out of range
   if ( xs < m_xmin )
   {
     cerr << "MaterialPfile::sample xs out of range (min): " <<  xs << ", " <<  m_xmin << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   if ( xs > m_xmax )
   {
     cerr << "MaterialPfile::sample xs out of range (max): " << xs << ", " << m_xmax << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
     
   if ( ys < m_ymin)
   {
     cerr << "MaterialPfile::sample ys out of range (min): " << ys << ", " <<  m_ymin << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   if( ys > m_ymax)
   {
     cerr << "MaterialPfile::sample ys out of range (max): " << ys << ", " << m_ymax << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   int ii, jj;
   if( m_nstenc % 2 == 1 )
   {
      // odd number of points in stencil
      int s = (m_nstenc-1)/2;
      ii  = static_cast<int>( floor( (xs-m_xmin+0.5*m_h)/m_h ) )-s;
      jj  = static_cast<int>( floor( (ys-m_ymin+0.5*m_h)/m_h ) )-s;
   }
   else
   {
      // even number of points in stencil
      int s = m_nstenc/2;
      ii  = static_cast<int>( floor( (xs-m_xmin)/m_h ) )-s+1;
      jj  = static_cast<int>( floor( (ys-m_ymin)/m_h ) )-s+1;
   }

// make sure we stay within array boundaries
   if( ii < 0 )
      ii = 0;
   if( jj < 0 )
      jj = 0;

   int ii2 = ii + m_nstenc-1;
   int jj2 = jj + m_nstenc-1;

// m_x = new double[m_nx], m_nx=m_nlat
   if( ii2 >= m_nx )
   {
     ii2 = m_nx-1;
     ii = ii2 -(m_nstenc-1);
   }
   
// m_y = new double[m_ny], m_ny=m_nlon
   if( jj2 >= m_ny )
   {
     jj2 = m_ny-1;
     jj = jj2 - (m_nstenc-1);
   }
   
   double re=6371e3;
   if( m_flatten ) 
      zs = re*(1.0 - exp(-zs/re));


   double w=0;
   zsed=zmoho=vp=vs=rho=qp=qs=0; 
   double appm  = 0.5*m_nstenc*m_h/sqrt(-log(1e-6));
   double appmi2 = 1.0/(appm*appm);

   for( int j1 = 0 ; j1 < jj2-jj+1 ; j1++ )
      for( int i1 = 0 ; i1 < ii2-ii+1 ; i1++ )
      {
	  double wgh = exp(-( (xs-m_x[ii+i1])*(xs-m_x[ii+i1])
			     +(ys-m_y[jj+j1])*(ys-m_y[jj+j1]) )*appmi2 );
	  w += wgh;
	  int m = (ii+i1)*m_nmaxdepth + (jj+j1)*m_nmaxdepth*m_nx;
          int kk;
	  for( kk=0; kk < m_nmaxdepth; kk++ )
	  {
	     if (m_z[m+kk] > zs) break;
	  }

	  if (kk+1 >= m_kmoho)
	     foundcrust = true;

	  int k1 = kk-1;
	  double factor = (zs-m_z[m+k1])/(m_z[m+k1+1]-m_z[m+k1]);
	  zsed  += m_z[m+m_ksed-1]*wgh;
	  zmoho += m_z[m+m_kmoho-1]*wgh;
	  vp    += (m_vp[m+k1]  + factor*(m_vp[m+k1+1]-m_vp[m+k1])  )*wgh;
	  vs    += (m_vs[m+k1]  + factor*(m_vs[m+k1+1]-m_vs[m+k1])  )*wgh;
	  rho   += (m_rho[m+k1] + factor*(m_rho[m+k1+1]-m_rho[m+k1]))*wgh;
	  if( m_qf )
	  {
	     qp += (m_qp[m+k1] + factor*(m_qp[m+k1+1]-m_qp[m+k1]))*wgh;
	     qs += (m_qs[m+k1] + factor*(m_qs[m+k1+1]-m_qs[m+k1]))*wgh;
	  }
      }

   // Now compute average properties by distance-weighted Gaussian average
   double iw;
   if (w != 0.)
     iw = 1.0/w;
   else
   {
     printf("Error MaterialPfile::sample_cart: weight w = 0 at x=%e, y=%e, depth=%e\n", xs, ys, zs);
// tmp
     printf("ii=%i, ii2=%i, jj=%i, jj2=%i, x[ii]=%e, y[jj]=%e\n", ii, ii2, jj, jj2, m_x[ii], m_y[jj]);
     double dist2    = (xs-m_x[ii])*(xs-m_x[ii]) +(ys-m_y[jj])*(ys-m_y[jj]);
     double exponent = -dist2*appmi2;
     double wgh      = exp(exponent);
     printf("dist2=%e, appm=%e, appmi2=%e, exponent=%e, wgh=%e\n", dist2, appm, appmi2, exponent, wgh);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   
   zsed  *= iw;
   zmoho *= iw;
   vp  *= iw;
   vs  *= iw;
   rho *= iw;
   if( m_qf ) 
   {
      qp *= iw;
      qs *= iw;
   }
   if( m_flatten )
   {
      zs  = re*log(re/(re-zs));
      vp  = vp*(re/(re-zs));
      vs  = vs*(re/(re-zs));
      rho = rho*pow(re/(re-zs),5.0);
   }
}

//-----------------------------------------------------------------------
MaterialPfile::~MaterialPfile()
{
   delete[] m_vp;
   delete[] m_vs;
   delete[] m_rho;
   delete[] m_z;
   delete[] m_st;
   delete[] m_ct;
   if( m_qf )
   {
      delete[] m_qp;
      delete[] m_qs;
   }
   if( m_coords_geographic )
   {
      delete[] m_lat;
      delete[] m_lon;
   }
   else
   {
      delete[] m_x;
      delete[] m_y;
   }
}
