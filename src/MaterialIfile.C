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
#include <mpi.h>

#include <iostream>
#include <cstring>
#include <unistd.h>

#include "MaterialIfile.h"
#include "EW.h"

using namespace std;

//-----------------------------------------------------------------------
MaterialIfile::MaterialIfile( EW * a_ew, std::string fileName, bool CartesianFormat )
{
   mEw   = a_ew;

// In general an ifile command only specifies properties down to a certain depth
   mCoversAllPoints = false;

// Is this file using Cartesian format?
   m_mat_Cartesian = CartesianFormat;
   
// read the grid surface
   if (CartesianFormat)
     extractSurfaceFromCartesianFile( fileName );
   else
     extractSurfaceFromGridFile( fileName );

}


//-----------------------------------------------------------------------
void MaterialIfile::set_material_properties( std::vector<Sarray> & rho, 
                                             std::vector<Sarray> & cs,
                                             std::vector<Sarray> & cp, 
                                             std::vector<Sarray> & qs, 
                                             std::vector<Sarray> & qp)
{
   if( mEw->m_materials.size() < m_number_material_surfaces )
   {
      cerr << "set_materials: ERROR: There are "<< m_number_material_surfaces << " material surfaces but only "
           << mEw->m_materials.size() << " defined materials" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   int totalPoints=0, definedPoints=0;
   int g;
   double lon, lat;
   float_sw4 x, y, z, depth;

   for( g = 0 ; g < mEw->mNumberOfCartesianGrids; g++) // Cartesian grids
   {
      int kLow = mEw->m_kStart[g];
      if (g == mEw->mNumberOfGrids-1 ) // no topography, so k=1 is at the top surface
	 kLow = 1;
#pragma omp parallel for reduction(+:totalPoints,definedPoints)
      for( int k = kLow ; k <= mEw->m_kEnd[g]; k++ )
      {

         
// what should the index boundaries be here to avoid parallel overlap points, but include real ghost points
	 for( int j = mEw->m_jStartInt[g]; j <= mEw->m_jEndInt[g]; j++ )
	 {
	    for( int i = mEw->m_iStartInt[g]; i <= mEw->m_iEndInt[g] ; i++ )
	    {
	       totalPoints += 1;
	       x = (i-1)*mEw->mGridSize[g];
	       y = (j-1)*mEw->mGridSize[g];
	       z = mEw->m_zmin[g]+(k-1)*mEw->mGridSize[g];
	  
	       mEw->getDepth(x, y, z, depth);
// (lon,lat) or (x,y) coordinates ?
	       if ( m_mat_Cartesian)
	       {
//		  if ( inside_cartesian_material_surfaces(x, y) )
// allow evaluation outside the surfaces (for ghost point values)
		  if ( true )
		  {
// interpolate material surfaces to find which material is at the specified depth 
		     int materialID = getCartesianMaterialID(x, y, depth );
		     if( 0 <= materialID && materialID < mEw->m_materials.size() )
		     {
			MaterialProperty* mprop = mEw->m_materials[materialID];
			definedPoints +=1;
			rho[g](i,j,k) = lookup_Rho(mprop, depth);
			cs[g](i,j,k)  = lookup_Vs(mprop, depth);
			cp[g](i,j,k)  = lookup_Vp(mprop, depth);
// attenuation model
			if( qp[g].is_defined())
			   qp[g](i,j,k) = lookup_Qp(mprop, depth);
			if( qs[g].is_defined())
			   qs[g](i,j,k) = lookup_Qs(mprop, depth);
		     } // end if knownMaterial
		  } // end if inside
	       } // end if (x,y)
	       else
	       {
// need to calculate (lon, lat) from (x, y) before calling getMaterialID
		  mEw->computeGeographicCoord(x, y, lon, lat);
		  if ( inside_material_surfaces(lat, lon))
		  {
// interpolate material surfaces to find which material is at the specified depth 
		     int materialID = getMaterialID(lat, lon, depth );
		     if( 0 <= materialID && materialID < mEw->m_materials.size() )
		     {
			MaterialProperty* mprop = mEw->m_materials[materialID];
			definedPoints +=1;
			rho[g](i,j,k) = lookup_Rho(mprop, depth);
			cs[g](i,j,k)  = lookup_Vs(mprop, depth);
			cp[g](i,j,k)  = lookup_Vp(mprop, depth);
// attenuation model
			if( qp[g].is_defined())
			   qp[g](i,j,k) = lookup_Qp(mprop, depth);
			if( qs[g].is_defined())
			   qs[g](i,j,k) = lookup_Qs(mprop, depth);
		       } // end if knownMaterial
		      } // end if inside
	       } // end if (lon,lat)
	    } // end for i
	 } // end for j

             cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> k=" << k << " z=" << z << " depth=" << depth << endl;

      } // end for k

// communicate material properties to ghost points (necessary on refined meshes because ghost points don't have a well defined depth/topography)
      mEw->communicate_array( rho[g], g );
      mEw->communicate_array( cs[g], g );
      mEw->communicate_array( cp[g], g );
// attenuation variables
      if( qp[g].is_defined())
	 mEw->communicate_array( qp[g], g );
      if( qs[g].is_defined())
	 mEw->communicate_array( qs[g], g );
   } // end for g (all Cartesian grids)
   if (mEw->topographyExists()) // curvilinear grid
   {
      int gTop = mEw->mNumberOfGrids-1;
      for (int g = mEw->mNumberOfCartesianGrids; g< mEw->mNumberOfGrids; g++)
      {
//      g = mEw->mNumberOfGrids-1; 
         int kLow = mEw->m_kStart[g];
         if (g == gTop)
            kLow = 1;

// Gradient and quadratic terms: Assume given (vp, vs, rho) constants are values at the previous depth surface.
//      double zsurf = 0; // set to zero for the first ifile command NEED TO GENERALIZE
#pragma omp parallel for reduction(+:totalPoints,definedPoints)
         for( int k = kLow ; k <= mEw->m_kEnd[g]; k++ )// don't attempt querying the ifile above the topography (start at k=1)
         {
            for( int j = mEw->m_jStart[g] ; j <= mEw->m_jEnd[g]; j++ )
            {
               for( int i = mEw->m_iStart[g] ; i <= mEw->m_iEnd[g] ; i++ )
               {
                  totalPoints += 1;
                  float_sw4 x = mEw->mX[g](i,j,k);
                  float_sw4 y = mEw->mY[g](i,j,k);
                  float_sw4 z = mEw->mZ[g](i,j,k);
                      
                  //printf("x ,y,z %f %f %f %f\n",x,y,z,mEw->m_zmin[g]);
                  float_sw4 depth;
                  mEw->getDepth(x, y, z, depth);
// (lon,lat) or (x,y) coordinates ?
                  if ( m_mat_Cartesian)
                  {
                     if ( inside_cartesian_material_surfaces(x, y))
                     {
// interpolate material surfaces to find which material is at the specified depth 
                        int materialID = getCartesianMaterialID(x, y, depth);
                        if( 0 <= materialID && materialID < mEw->m_materials.size() )
                        {
                           MaterialProperty* mprop = mEw->m_materials[materialID];
                           definedPoints +=1;
                           rho[g](i,j,k) = lookup_Rho(mprop, depth);
                           cs[g](i,j,k)  = lookup_Vs(mprop, depth);
                           cp[g](i,j,k)  = lookup_Vp(mprop, depth);
                           // attenuation model
                           if( qp[g].is_defined())
                              qp[g](i,j,k) = lookup_Qp(mprop, depth);
                           if( qs[g].is_defined())
                              qs[g](i,j,k) = lookup_Qs(mprop, depth);
                        } // end if knownMaterial
                     } // end if inside
                  } // end if (x,y)
                  else
                  {
// need to calculate (lon, lat) from (x, y) before calling getMaterialID
                     mEw->computeGeographicCoord(x, y, lon, lat );
                     if ( inside_material_surfaces(lat, lon))
                     {
// interpolate material surfaces to find which material is at the specified depth 
                        int materialID = getMaterialID(lat, lon, depth );
                        if( 0 <= materialID && materialID < mEw->m_materials.size() )
                        {
                           MaterialProperty* mprop = mEw->m_materials[materialID];
                           definedPoints +=1;
                           rho[g](i,j,k) = lookup_Rho(mprop, depth);
                           cs[g](i,j,k)  = lookup_Vs(mprop, depth);
                           cp[g](i,j,k)  = lookup_Vp(mprop, depth);
// attenuation model
                           if( qp[g].is_defined())
                              qp[g](i,j,k) = lookup_Qp(mprop, depth);
                           if( qs[g].is_defined())
                              qs[g](i,j,k) = lookup_Qs(mprop, depth);
                        } // end if knownMaterial
                     } // end if inside
                  } // end (lon, lat)
               } // end for i
            }// end for j
         } // end for k
      } // end for g (curvilinear)
      
   } // end if topographyExists
   int totalPointsSum, materialSum;
   MPI_Reduce(&totalPoints, &totalPointsSum, 1, MPI_INT, MPI_SUM, 0, mEw->m_cartesian_communicator);
   MPI_Reduce(&definedPoints, &materialSum, 1, MPI_INT, MPI_SUM, 0, mEw->m_cartesian_communicator);
   if (mEw->proc_zero())
      cout << "MaterialIfile:: set_material_properties: Total # points=" << totalPointsSum 
	   << ". Defined properties in # points=" << materialSum << endl;
} // end MaterialIfile::set_material_properties


//-----------------------------------------------------------------------
void MaterialIfile::extractSurfaceFromGridFile(string a_surfaceFileName)
{
   if (mEw->proc_zero())
      cout << "***inside MaterialIfile::extractSurfaceFromGridFile***"<< endl;

//----------------------------------------------
// Check user specified file name. Abort if they are not there or not readable
//----------------------------------------------
   VERIFY2(access(a_surfaceFileName.c_str(), R_OK) == 0,
	       "No read permission on ifile surface file: " << a_surfaceFileName);

   float_sw4 x, y, depth;
   double lat, lon;
   char buffer[256];
   char *token;

// 1. read the grid file
   int Nlon, Nlat, Nmat, i, j, k;
  
   ifstream gridFile(a_surfaceFileName.c_str());
   if (!gridFile.is_open())
   {    
      cerr << "# Error opening ifile surface file: '" << a_surfaceFileName << "'" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
  
//  gridFile >> Nlon >> Nlat >> Nmat;
   gridFile.getline(buffer, 256);
   token = strtok(buffer, " \t");
   if (token != NULL)
      Nlon = atoi(token);
   if ((token = strtok(NULL, " \t")) != NULL)
      Nlat = atoi(token);
   if ((token = strtok(NULL, " \t")) != NULL)
      Nmat = atoi(token);
  
   VERIFY2( Nlon>0 && Nlat>0 && Nmat>0, "ERROR: MaterialIfile:: extractSurfaceFromGridFile: Nlon, Nlat, Nmat= " << Nlon << " " << Nlat << " " << Nmat);

   m_number_material_surfaces = Nmat;
   m_Nlon = Nlon;
   m_Nlat = Nlat;
  
   m_materialDepth.define(Nmat,1,Nlon,1,Nlat,1,1);
   m_materialLat = new double[Nlat+1];
   m_materialLon = new double[Nlon+1];

   bool missingValues;
   for (j=1; j<=Nlat; j++)
      for (i=1; i<=Nlon; i++)
      {
// only works for one surface
//      gridFile >> m_materialLon[i] >> m_materialLat[j] >> m_materialDepth(1,i,j,1);
	 missingValues=false;
	 int nValues=0;
      
	 gridFile.getline(buffer, 256);
	 token = strtok(buffer, " \t");
	 if (token != NULL)
	 {
	    m_materialLon[i] = atof(token);
	    nValues++;
	 }
	 else
	    missingValues=true;
      
	 if (!missingValues && (token = strtok(NULL, " \t")) != NULL)
	 {
	    m_materialLat[j] = atof(token);
	    nValues++;
	 }
	 else
	    missingValues=true;

// read all depth values
	 if (!missingValues)
	 {
	    for (k=1; k<=Nmat; k++)
	    {
	       if ((token = strtok(NULL, " \t")) != NULL)
	       {
		  m_materialDepth(k,i,j,1) = atof(token);
		  nValues++;
	       }
	       else
	       {
		  missingValues=true;
		  break;
	       }
	    }// end for k
	 }
	 if (missingValues)
	 {
	    printf("extractSurfaceFromGridFile:: Fatal error: Expecting %i values on line %i, but only found: %i\n", 
	       2+Nmat, 1+(j-1)*Nlon+i, nValues);
	MPI_Abort(MPI_COMM_WORLD,1);
	 }
      
      } // end for i
   gridFile.close();
   if (mEw->proc_zero())
      printf("Nlon=%i Nlat=%i Nmat=%i\n", Nlon, Nlat, Nmat);

   m_materialLonMax=-1e30, m_materialLonMin=1e30, m_materialLatMax=-1e30, m_materialLatMin=1e30;
   float_sw4 depthMax=-1e10, depthMin=1e10;
   for (i=1; i<=Nlon; i++)
   {
      if (m_materialLon[i] < m_materialLonMin) m_materialLonMin=m_materialLon[i];
      if (m_materialLon[i] > m_materialLonMax) m_materialLonMax=m_materialLon[i];
   }
   for (i=1; i<=Nlat; i++)
   {
      if (m_materialLat[i] < m_materialLatMin) m_materialLatMin=m_materialLat[i];
      if (m_materialLat[i] > m_materialLatMax) m_materialLatMax=m_materialLat[i];
   }
  
   for (j=1; j<=Nlat; j++)
      for (i=1; i<=Nlon; i++)
	 for (k=1; k<=Nmat; k++)
	 {
	    if (m_materialDepth(k,i,j,1) < depthMin) depthMin=m_materialDepth(k,i,j,1);
	    if (m_materialDepth(k,i,j,1) > depthMax) depthMax=m_materialDepth(k,i,j,1);
	 }
   if (mEw->proc_zero())
      printf("Material interface surfaces: lonMin=%e, lonMax=%e\nlatMin=%e, latMax=%e\ndepthMin=%e, depthMax=%e\n", 
	     m_materialLonMin, m_materialLonMax, m_materialLatMin, m_materialLatMax, depthMin, depthMax);
  
// If the lat vector is not in increasing order, we need to reorder it
   if (m_materialLat[1] > m_materialLat[Nlat])
   {
      if (mEw->proc_zero()) 
	 printf("Reordering the latitude vector...\n");
      for (j=1; j<=Nlat/2; j++)
      {
	 lat=m_materialLat[Nlat+1-j];
	 m_materialLat[Nlat+1-j] = m_materialLat[j];
	 m_materialLat[j] = lat;
	 
	 for (i=1; i<=Nlon; i++)
	    for (k=1; k<=Nmat; k++)
	    {
	       depth = m_materialDepth(k,i,Nlat+1-j,1);
	       m_materialDepth(k,i,Nlat+1-j,1) = m_materialDepth(k,i,j,1);
	       m_materialDepth(k,i,j,1) = depth;
	    } // end for k,i
      
      }// end for j    
   } // end if m_materialLat[1] > m_materialLat[Nlat]
  
// If the lon vector is not in increasing order, we need to reorder it
   if (m_materialLon[1] > m_materialLon[Nlon])
   {
      if (mEw->proc_zero()) 
	 printf("Reordering the longitude vector...\n");
      for (i=1; i<=Nlon/2; i++)
      {
	 lon=m_materialLon[Nlon+1-i];
	 m_materialLon[Nlon+1-i] = m_materialLon[i];
	 m_materialLon[i] = lon;
      
	 for (j=1; j<=Nlat; j++)
	    for (k=1; k<=Nmat; k++)
	    {
	       depth = m_materialDepth(k,Nlon+1-i,j,1);
	       m_materialDepth(k,Nlon+1-i,j,1) = m_materialDepth(k,i,j,1);
	       m_materialDepth(k,i,j,1) = depth;
	    }// end for k,i
      }// end for i    
   } // end if m_materialLon[1] > m_materialLon[Nlon]
}

//-----------------------------------------------------------------------
void MaterialIfile::extractSurfaceFromCartesianFile(string a_surfaceFileName)
{
   if (mEw->proc_zero())
      cout << "***inside EW::extractSurfaceFromCartesianFile***"<< endl;

//----------------------------------------------
// Check user specified file name. Abort if they are not there or not readable
//----------------------------------------------
   VERIFY2(access(a_surfaceFileName.c_str(), R_OK) == 0,
	       "No read permission on ifile surface file: " << a_surfaceFileName);

   float_sw4 x, y;
   float_sw4 y0, x0, depth;
   char buffer[256];
   char *token;

// 1. read the grid file
   int Nx, Ny, Nmat, i, j, k;
  
   ifstream gridFile(a_surfaceFileName.c_str());
   if (!gridFile.is_open())
   {    
     cerr << "# Error opening ifile surface file: '" << a_surfaceFileName << "'" << endl;
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
  
//  gridFile >> Nx >> Ny >> Nmat;
   gridFile.getline(buffer, 256);
   token = strtok(buffer, " \t");
   if (token != NULL)
      Nx = atoi(token);
   if ((token = strtok(NULL, " \t")) != NULL)
      Ny = atoi(token);
   if ((token = strtok(NULL, " \t")) != NULL)
      Nmat = atoi(token);
  
   VERIFY2( Nx>0 && Ny>0 && Nmat>0, "ERROR: EW:: extractSurfaceFromCartesianFile: Nx, Ny, Nmat= " <<
	       Nx << " " << Ny << " " << Nmat);

   m_number_material_surfaces = Nmat;
   m_mat_Nx = Nx;
   m_mat_Ny = Ny;
  
   m_materialDepth.define(Nmat,1,Nx,1,Ny,1,1);
   m_mat_Yvec = new float_sw4[Ny+1];
   m_mat_Xvec = new float_sw4[Nx+1];

   bool missingValues;
   for (j=1; j<=Ny; j++)
      for (i=1; i<=Nx; i++)
      {
// only works for one surface
//      gridFile >> m_mat_Xvec[i] >> m_mat_Yvec[j] >> m_materialDepth(1,i,j,1);
	 missingValues=false;
	 int nValues=0;
      
	 gridFile.getline(buffer, 256);
	 token = strtok(buffer, " \t");
	 if (token != NULL)
	 {
	    m_mat_Xvec[i] = atof(token);
	    nValues++;
	 }
	 else
	    missingValues=true;
	 
	 if (!missingValues && (token = strtok(NULL, " \t")) != NULL)
	 {
	    m_mat_Yvec[j] = atof(token);
	    nValues++;
	 }
	 else
	    missingValues=true;

// read all depth values
	 if (!missingValues)
	 {
	    for (k=1; k<=Nmat; k++)
	    {
	       if ((token = strtok(NULL, " \t")) != NULL)
	       {
		  m_materialDepth(k,i,j,1) = atof(token);
		  nValues++;
	       }
	       else
	       {
		  missingValues=true;
		  break;
	       }
	    }// end for k
	 }
	 if (missingValues)
	 {
	    printf("extractSurfaceFromCartesianFile:: Fatal error: Expecting %i values on line %i, but only found: %i\n", 
		   2+Nmat, 1+(j-1)*Nx+i, nValues);
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
      
      } // end for i
   gridFile.close();
  
   if (mEw->proc_zero())
      printf("Nx=%i Ny=%i Nmat=%i\n", Nx, Ny, Nmat);

   m_mat_Xmax=-1e30, m_mat_Xmin=1e30, m_mat_Ymax=-1e30, m_mat_Ymin=1e30;
   float_sw4 depthMax=-1e10, depthMin=1e10;
   for (i=1; i<=Nx; i++)
   {
      if (m_mat_Xvec[i] < m_mat_Xmin) m_mat_Xmin=m_mat_Xvec[i];
      if (m_mat_Xvec[i] > m_mat_Xmax) m_mat_Xmax=m_mat_Xvec[i];
   }
   for (i=1; i<=Ny; i++)
   {
      if (m_mat_Yvec[i] < m_mat_Ymin) m_mat_Ymin=m_mat_Yvec[i];
      if (m_mat_Yvec[i] > m_mat_Ymax) m_mat_Ymax=m_mat_Yvec[i];
   }
  
   for (j=1; j<=Ny; j++)
      for (i=1; i<=Nx; i++)
	 for (k=1; k<=Nmat; k++)
	 {
	    if (m_materialDepth(k,i,j,1) < depthMin) depthMin=m_materialDepth(k,i,j,1);
	    if (m_materialDepth(k,i,j,1) > depthMax) depthMax=m_materialDepth(k,i,j,1);
	 }
   if (mEw->proc_zero())
      printf("Material interface surfaces: Xmin=%e, Xmax=%e\nYmin=%e, Ymax=%e\ndepthMin=%e, depthMax=%e\n", 
	     m_mat_Xmin, m_mat_Xmax, m_mat_Ymin, m_mat_Ymax, depthMin, depthMax);
  
// If the Y-vector is not in increasing order, we need to reorder it
   if (m_mat_Yvec[1] > m_mat_Yvec[Ny])
   {
      if (mEw->proc_zero())
	 printf("Reordering the Y-coordinate vector...\n");
      for (j=1; j<=Ny/2; j++)
      {
	 y0=m_mat_Yvec[Ny+1-j];
	 m_mat_Yvec[Ny+1-j] = m_mat_Yvec[j];
	 m_mat_Yvec[j] = y0;
      
	 for (i=1; i<=Nx; i++)
	    for (k=1; k<=Nmat; k++)
	    {
	       depth = m_materialDepth(k,i,Ny+1-j,1);
	       m_materialDepth(k,i,Ny+1-j,1) = m_materialDepth(k,i,j,1);
	       m_materialDepth(k,i,j,1) = depth;
	    } // end for k,i
      
      }// end for j    
   } // end if m_mat_Yvec[1] > m_mat_Yvec[Ny]
  
// If the X-vector is not in increasing order, we need to reorder it
   if (m_mat_Xvec[1] > m_mat_Xvec[Nx])
   {
      if (mEw->proc_zero())
	 printf("Reordering the X-coordinate vector...\n");
      for (i=1; i<=Nx/2; i++)
      {
	 x0=m_mat_Xvec[Nx+1-i];
	 m_mat_Xvec[Nx+1-i] = m_mat_Xvec[i];
	 m_mat_Xvec[i] = x0;
	 
	 for (j=1; j<=Ny; j++)
	    for (k=1; k<=Nmat; k++)
	    {
	       depth = m_materialDepth(k,Nx+1-i,j,1);
	       m_materialDepth(k,Nx+1-i,j,1) = m_materialDepth(k,i,j,1);
	       m_materialDepth(k,i,j,1) = depth;
	    }// end for k,i
      }// end for i    
   } // end if m_mat_Xvec[1] > m_mat_Xvec[Nx]
}

//-----------------------------------------------------------------------
int MaterialIfile::getMaterialID(double lat, double lon, float_sw4 depth )
{
// Interpolate the materialDepth surfaces to find which material is at the specified depth 

// 1. Interpolate in the grid file to get elevations on the computational grid
   double deltaLat = (m_materialLatMax-m_materialLatMin)/m_Nlat;
   double deltaLon = (m_materialLonMax-m_materialLonMin)/m_Nlon;
   double eInterp, xi, eta;
   int i0, j0;
   if (lat > m_materialLatMax || lat < m_materialLatMin || lon > m_materialLonMax || lon < m_materialLonMin)
   {
      printf("Warning: MaterialIfile::getMaterialID lon=%e, lat=%e is outside the material surface grid\n",
	     lon, lat);
      //      MPI_Abort(MPI_COMM_WORLD,1);
      return -1;
   }
   i0 = 1+(int)((lon-m_materialLonMin)/deltaLon);
   j0 = 1+(int)((lat-m_materialLatMin)/deltaLat);
   while ( lon < m_materialLon[i0] || m_materialLon[i0+1] < lon )
   {
      if (lon<m_materialLon[i0]) 
	 i0--;
      else if (lon>m_materialLon[i0+1])
	 i0++;
   }
   while (  lat < m_materialLat[j0] || m_materialLat[j0+1] < lat )
   {
      if (lat<m_materialLat[j0]) 
	 j0--;
      else if (lat>m_materialLat[j0+1])
	 j0++;
   }
   if (i0 > m_Nlon-1) i0 = m_Nlon-1;
   if (j0 > m_Nlat-1) j0 = m_Nlat-1;
      
// Test that we are inside the interval
   if (!(m_materialLon[i0] <= lon && lon < m_materialLon[i0+1]))
   {
      printf("lon=%e outside (m_materialLon[%i]=%e, m_materialLon[%i]=%e)\n",
	     lon, i0, m_materialLon[i0], i0+1, m_materialLon[i0+1]);
      return -1;
   }
   if (!(m_materialLat[j0] <= lat && lat < m_materialLat[j0+1]))
   {
      printf("lat=%e outside (m_materialLat[%i]=%e, m_materialLat[%i]=%e)\n",
	     lat, j0, m_materialLat[j0], j0+1, m_materialLat[j0+1]);
      return -1;
   }
// Local step sizes
   xi  = (lon - m_materialLon[i0])/(m_materialLon[i0+1]-m_materialLon[i0]);
   eta = (lat - m_materialLat[j0])/(m_materialLat[j0+1]-m_materialLat[j0]);
      
// Bi-linear interpolation for all surfaces
   int materialID=-1; // default material is unknown
   float_sw4 minDepth=0, maxDepth; 
   for (int q=1; q<=m_number_material_surfaces; q++) // surfaces are assumed to be in increasing depth order
   {
      maxDepth = (1.0-eta)*( (1.0-xi)*m_materialDepth(q,i0,j0,1) + xi*m_materialDepth(q,i0+1,j0,1) ) +
	 eta*( (1.0-xi)*m_materialDepth(q,i0,j0+1,1) + xi*m_materialDepth(q,i0+1,j0+1,1) );

      if (maxDepth > minDepth && depth <= maxDepth && depth >= minDepth)
	 // maxDepth>minDepth removes zero thickness layers
      {
	 materialID = q-1;
	 break;
      }
      minDepth=maxDepth; // next set of surfaces
   } // end for q
   
   return materialID;
}

//-----------------------------------------------------------------------
int MaterialIfile::getCartesianMaterialID(float_sw4 xP, float_sw4 yP, float_sw4 depth )
{
// interpolate the materialDepth surfaces to find which material is at the specified depth 
  
// Interpolate in the grid file to get elevations on the computational grid
// allow extrapolation (for ghost points)
   // if (yP > m_mat_Ymax || yP < m_mat_Ymin || xP > m_mat_Xmax || xP < m_mat_Xmin)
   // {
   //    printf("Warning: MaterialIfile::getCartesianMaterialID xP=%e, yP=%e is outside the material surface grid\n", xP, yP);
   //    //      MPI_Abort(MPI_COMM_WORLD,1);
   //    return -1;
   // }
   float_sw4 deltaY = (m_mat_Ymax-m_mat_Ymin)/m_mat_Ny;
   float_sw4 deltaX = (m_mat_Xmax-m_mat_Xmin)/m_mat_Nx;
   int i0 = 1+(int)((xP-m_mat_Xmin)/deltaX);
   int j0 = 1+(int)((yP-m_mat_Ymin)/deltaY);

// enforce the array bounds initially
   if (i0 < 1) i0 = 1;
   if (i0 > m_mat_Nx-1) i0 = m_mat_Nx-1;
   if (j0 < 1) j0 = 1;
   if (j0 > m_mat_Ny-1) j0 = m_mat_Ny-1;

// find the right interval
   while ( i0 >= 1 && i0 < m_mat_Nx && (xP < m_mat_Xvec[i0] || m_mat_Xvec[i0+1] < xP) )
   {
      if (xP<m_mat_Xvec[i0]) 
	 i0--;
      else if (xP>m_mat_Xvec[i0+1])
	 i0++;
   }
   while ( j0 >= 1 && j0 < m_mat_Ny && (yP < m_mat_Yvec[j0] || m_mat_Yvec[j0+1] < yP) )
   {
      if (yP<m_mat_Yvec[j0]) 
	 j0--;
      else if (yP>m_mat_Yvec[j0+1])
	 j0++;
   }
      
// enforce the array bounds again
   if (i0 < 1) i0 = 1;
   if (i0 > m_mat_Nx-1) i0 = m_mat_Nx-1;
   if (j0 < 1) j0 = 1;
   if (j0 > m_mat_Ny-1) j0 = m_mat_Ny-1;

// remove this test to allow extrapolation to ghost points
// test that we are inside the interval
   // if (!(m_mat_Xvec[i0] <= xP && xP <= m_mat_Xvec[i0+1]))
   // {
   //    printf("xP=%e outside (m_mat_Xvec[%i]=%e, m_mat_Xvec[%i]=%e)\n",
   // 	     xP, i0, m_mat_Xvec[i0], i0+1, m_mat_Xvec[i0+1]);
   //    return -1;
   // }
  
   // if (!(m_mat_Yvec[j0] <= yP && yP <= m_mat_Yvec[j0+1]))
   // {
   //    printf("yP=%e outside (m_mat_Yvec[%i]=%e, m_mat_Yvec[%i]=%e)\n",
   // 	     yP, j0, m_mat_Yvec[j0], j0+1, m_mat_Yvec[j0+1]);
   //    return -1;
   // }
  
// local step sizes
   float_sw4 xi, eta;
   xi  = (xP - m_mat_Xvec[i0])/(m_mat_Xvec[i0+1]-m_mat_Xvec[i0]);
   eta = (yP - m_mat_Yvec[j0])/(m_mat_Yvec[j0+1]-m_mat_Yvec[j0]);
   // bi-linear interpolation for all surfaces
   int materialID=-1; // default material is unknown
   float_sw4 minDepth=0, maxDepth; 
   for (int q=1; q<=m_number_material_surfaces; q++) // surfaces are assumed to be in increasing depth order
   {
      maxDepth = (1.0-eta)*( (1.0-xi)*m_materialDepth(q,i0,j0,1) + xi*m_materialDepth(q,i0+1,j0,1) ) +
	 eta*( (1.0-xi)*m_materialDepth(q,i0,j0+1,1) + xi*m_materialDepth(q,i0+1,j0+1,1) );

      if (maxDepth >= minDepth+.1 && depth <= maxDepth && depth >= minDepth) // maxDepth>=minDepth+.1 removes zero thickness layers
      {
	 materialID = q-1;
	 break;
      }
      minDepth=maxDepth; // next set of surfaces
   } // end for q
   return materialID;
}
