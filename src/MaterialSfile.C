// -*-c++-*-
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

#include "Require.h"

#include <cstring>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <fcntl.h>
#include "EW.h"
#include "MaterialSfile.h"
#include "Byteswapper.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

using namespace std;


//-----------------------------------------------------------------------
MaterialSfile::MaterialSfile( EW* a_ew, const string a_file, const string a_directory):
   mEW(a_ew),
   m_model_file(a_file),
   m_model_dir(a_directory),
   m_use_attenuation(false)
{
   mCoversAllPoints = false;
   // Check that the depths make sense
   if (a_ew != NULL) {
     m_use_attenuation = a_ew->usingAttenuation();
     read_sfile();
   }
}

//-----------------------------------------------------------------------
MaterialSfile::~MaterialSfile()
{
}

//-----------------------------------------------------------------------
void MaterialSfile::set_material_properties(std::vector<Sarray> & rho, 
                                             std::vector<Sarray> & cs,
                                             std::vector<Sarray> & cp, 
                                             std::vector<Sarray> & xis, 
                                             std::vector<Sarray> & xip )
{
// Assume attenuation arrays defined on all grids if they are defined on grid zero.
   bool use_q = m_use_attenuation && xis[0].is_defined() && xip[0].is_defined();
   size_t outside=0, material=0;
   float_sw4 z_min = m_zminloc;

   // Find the relative dimension size of upper and lower interface for each grid patch
   int* ist = new int[m_npatches];
   int* jst = new int[m_npatches];
   for( int g=0 ; g < m_npatches ; g++ ) {
     ist[g] = (int)ceil((double)mInterface[g].m_ni / mInterface[g+1].m_ni);
     jst[g] = (int)ceil((double)mInterface[g].m_nj / mInterface[g+1].m_nj);
   }

   for( int g=0 ; g < mEW->mNumberOfGrids ; g++ ) {
      bool curvilinear = mEW->topographyExists() && g >= mEW->mNumberOfCartesianGrids;
      size_t ni=mEW->m_iEnd[g]-mEW->m_iStart[g]+1;
      size_t nj=mEW->m_jEnd[g]-mEW->m_jStart[g]+1;
#pragma omp parallel for reduction(+:material,outside)
      for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) {
	 for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j) {
	    for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) {
                float_sw4 z0, hv;
		float_sw4 x = (i-1)*mEW->mGridSize[g];
		float_sw4 y = (j-1)*mEW->mGridSize[g];
		float_sw4 z;
		if( curvilinear )
		   z = mEW->mZ[g](i,j,k);
		else
		   z = mEW->m_zmin[g] + (k-1)*mEW->mGridSize[g];

                // Deal with some values on top grid that exceeds the topogrophy interface
                if (g == mEW->mNumberOfGrids - 1 && z < z_min) 
                  z = z_min;

                // (x, y, z) is the coordinate of current grid point
		/* if( inside( x, y, z ) ) { */
		if( m_zminloc <= z && z <= m_zmaxloc ) {
                   // Extend the material value if simulation grid is larger than material grid
                   if (x > m_xmaxloc)
                       x = m_xmaxloc;
                   if (y > m_ymaxloc)
                       y = m_ymaxloc;
                       
		   material++;
                   int i0, j0, i1, j1, k0, gr = m_npatches-1;
                   float_sw4 tmph, down_z;
                   // gr is the patch id that has the current sw4 grid point's data
                   // need to use the interface value to determine which patch the current point is in
		   while( gr >= 0 ) {
                     i0 = i1 = static_cast<int>( trunc( 1 + (x-m_x0)/m_hh[gr] ) );
                     j0 = j1 = static_cast<int>( trunc( 1 + (y-m_y0)/m_hh[gr] ) );
                     down_z = mInterface[gr+1](1, i0, j0, 1);

                     // Adjust the index if upper and lower interface have different dimension
                     if (ist[gr] > 1) { i1 = i1*ist[gr]-1; }
                     if (jst[gr] > 1) { j1 = j1*jst[gr]-1; }
                     z0 = mInterface[gr](1, i1, j1, 1);

                     if (gr == 0 || z > z0) break;
		     gr--;
                   }

                   tmph = down_z - z0;
                   if (z > down_z)
                       z = down_z;
                   if (z < z0)     
                       z = z0;

                   // Update the current vertical grid height and z-base with the sfile curvilinear grid
                   hv = tmph / (m_nk[gr]-1);

                   // we are using curvilinear grid in sfile 
                   k0 = static_cast<int>( trunc( 1 + (z-z0)/hv) );

		   // Use bilinear interpolation always:
        	   // Bias stencil near the boundary, need to communicate arrays afterwards.
		   if( i0 <= m_ifirst[gr] ) 
		      i0 = m_ifirst[gr];

		   if( i0 >= m_ilast[gr]-1 ) 
		      i0 = m_ilast[gr]-1;
		    
		   if( j0 <= m_jfirst[gr] ) 
		      j0 = m_jfirst[gr];
		    
		   if( j0 >= m_jlast[gr]-1 ) 
		      j0 = m_jlast[gr]-1;
		    
		   if( k0 <= m_kfirst[gr] ) 
		      k0 = m_kfirst[gr];
		    
		   if( k0 >= m_klast[gr]-1 ) 
		      k0 = m_klast[gr]-1;

   		   // bilinear intp.
                   float_sw4 wghx = (x-( (i0-1)*m_hh[gr]+m_x0) )/m_hh[gr];
                   float_sw4 wghy = (y-( (j0-1)*m_hh[gr]+m_y0) )/m_hh[gr];
                   float_sw4 wghz = (z-( (k0-1)*hv+z0) )/hv;

                   // Debug
                   // weights should be within [0, 1]
                   if (wghx > 1 || wghx < 0) { 
                       printf("g=%d, sw4 (%d, %d, %d), mat (%d, %d, %d) wghx = %.2f\n", gr, i, j, k, i0, j0, k0, wghx);
                       if (wghx > 1) wghx = 1;
                       if (wghx < 0) wghx = 0;
                   }

                   if (wghy > 1 || wghy < 0) { 
                       printf("g=%d, sw4 (%d, %d, %d), mat (%d, %d, %d) wghy = %.2f\n", gr, i, j, k, i0, j0, k0, wghy);
                       if (wghy > 1) wghy = 1;
                       if (wghy < 0) wghy = 0;
                   }

                   if (wghz > 1 || wghz < 0) { 
                       if (wghz > 1.001 || wghz < -0.001) 
                         printf("g=%d, sw4 (%d, %d, %d), mat (%d, %d, %d) wghz = %.2f, z=%.2f, z0=%.2f\n", 
                                 gr, i, j, k, i0, j0, k0, wghz, z, z0);
                       if (wghz > 1) wghz = 1;
                       if (wghz < 0) wghz = 0;
                   }

                   rho[g](i, j, k) = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_rho[gr](1,i0,j0,k0) + 
                                     wghx*mMaterial_rho[gr](1,i0+1,j0,k0) ) +
   		                        wghy*(     (1-wghx)*mMaterial_rho[gr](1,i0,j0+1,k0) +
                                     wghx*mMaterial_rho[gr](1,i0+1,j0+1,k0) ) ) + 
   		                     wghz*( (1-wghy)*( (1-wghx)*mMaterial_rho[gr](1,i0,j0,k0+1) +
                                     wghx*mMaterial_rho[gr](1,i0+1,j0,k0+1) ) +
   		     	             wghy*(    (1-wghx)*mMaterial_rho[gr](1,i0,j0+1,k0+1) + 
                                     wghx*mMaterial_rho[gr](1,i0+1,j0+1,k0+1) ) );

                   cp[g](i, j, k)  = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_cp[gr](1,i0,j0,k0) + 
                                     wghx*mMaterial_cp[gr](1,i0+1,j0,k0) ) +
                                     wghy*(     (1-wghx)*mMaterial_cp[gr](1,i0,j0+1,k0) + 
                                     wghx*mMaterial_cp[gr](1,i0+1,j0+1,k0) ) ) + 
                                     wghz*(  (1-wghy)*( (1-wghx)*mMaterial_cp[gr](1,i0,j0,k0+1) + 
                                     wghx*mMaterial_cp[gr](1,i0+1,j0,k0+1) ) +
                                     wghy*(     (1-wghx)*mMaterial_cp[gr](1,i0,j0+1,k0+1)+ 
                                     wghx*mMaterial_cp[gr](1,i0+1,j0+1,k0+1) ) );
       
                   cs[g](i, j, k)  = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_cs[gr](1,i0,j0,k0) + 
                                     wghx*mMaterial_cs[gr](1,i0+1,j0,k0) ) +
                                     wghy*(    (1-wghx)*mMaterial_cs[gr](1,i0,j0+1,k0) + 
                                     wghx*mMaterial_cs[gr](1,i0+1,j0+1,k0) ) ) + 
                                     wghz*(   (1-wghy)*( (1-wghx)*mMaterial_cs[gr](1,i0,j0,k0+1) + 
                                     wghx*mMaterial_cs[gr](1,i0+1,j0,k0+1) ) +
                                     wghy*(   (1-wghx)*mMaterial_cs[gr](1,i0,j0+1,k0+1)+ 
                                     wghx*mMaterial_cs[gr](1,i0+1,j0+1,k0+1) ) );

                   if( use_q ) {
                      xip[g](i, j, k) = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_qp[gr](1,i0,j0,k0) + 
                                        wghx*mMaterial_qp[gr](1,i0+1,j0,k0) ) +
                                        wghy*( (1-wghx)*mMaterial_qp[gr](1,i0,j0+1,k0) + 
                                        wghx*mMaterial_qp[gr](1,i0+1,j0+1,k0) ) ) + 
                                        wghz*( (1-wghy)*(   (1-wghx)*mMaterial_qp[gr](1,i0,j0,k0+1) + 
                                        wghx*mMaterial_qp[gr](1,i0+1,j0,k0+1) ) +
                                        wghy*( (1-wghx)*mMaterial_qp[gr](1,i0,j0+1,k0+1)+ 
                                        wghx*mMaterial_qp[gr](1,i0+1,j0+1,k0+1) ) );

                      xis[g](i, j, k) = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_qs[gr](1,i0,j0,k0) + 
                                        wghx*mMaterial_qs[gr](1,i0+1,j0,k0) ) +
                                        wghy*(      (1-wghx)*mMaterial_qs[gr](1,i0,j0+1,k0) + 
                                        wghx*mMaterial_qs[gr](1,i0+1,j0+1,k0) ) ) + 
                                        wghz*(   (1-wghy)*(   (1-wghx)*mMaterial_qs[gr](1,i0,j0,k0+1) + 
                                        wghx*mMaterial_qs[gr](1,i0+1,j0,k0+1) ) +
                                        wghy*(     (1-wghx)*mMaterial_qs[gr](1,i0,j0+1,k0+1)+ 
                                        wghx*mMaterial_qs[gr](1,i0+1,j0+1,k0+1) ) );
                   }

		} // End if inside
		else
		   outside++;
	    } // End for i
	  } // End for j
        } // End for k

    /* debug */ 
    /* if (mEW->getRank() == 0) */ 
    /*   printf("After interpolation from sfile\n"); */

    /* printf("Rank %d grid %d: rho min = %.2f, max = %.2f\n", mEW->getRank(), g, rho[g].minimum(), rho[g].maximum()); */
    /* printf("Rank %d grid %d: cp min = %.2f, max = %.2f\n", mEW->getRank(), g, cp[g].minimum(), cp[g].maximum()); */
    /* printf("Rank %d grid %d: cs min = %.2f, max = %.2f\n", mEW->getRank(), g, cs[g].minimum(), cs[g].maximum()); */

    /* if( use_q ) { */
    /*   printf("Rank %d grid %d: xip min = %.2f, max = %.2f\n", mEW->getRank(), g, xip[g].minimum(), xip[g].maximum()); */
    /*   printf("Rank %d grid %d: xis min = %.2f, max = %.2f\n", mEW->getRank(), g, xis[g].minimum(), xis[g].maximum()); */
    /* } */

   } // end for g...

   delete [] ist;
   delete [] jst;

   mEW->communicate_arrays( rho );
   mEW->communicate_arrays( cs );
   mEW->communicate_arrays( cp );
   mEW->material_ic( rho );
   mEW->material_ic( cs );
   mEW->material_ic( cp );
   if( use_q ) {
      mEW->communicate_arrays( xis );
      mEW->communicate_arrays( xip );
      mEW->material_ic( xis);
      mEW->material_ic( xip );
   }

   size_t materialSum, outsideSum;
   int mpisizelong, mpisizelonglong, mpisizeint;
   MPI_Type_size(MPI_LONG,&mpisizelong );
   MPI_Type_size(MPI_LONG_LONG,&mpisizelonglong );
   MPI_Type_size(MPI_INT,&mpisizeint );
   if( sizeof(size_t) == mpisizelong ) {
      MPI_Reduce(&material, &materialSum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else if( sizeof(size_t) == mpisizelonglong ) {
      MPI_Reduce(&material, &materialSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else if( sizeof(size_t) == mpisizeint ) {
      MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else {
      int materialsumi, outsidesumi, materiali=material, outsidei=outside;
      MPI_Reduce(&materiali, &materialsumi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outsidei,   &outsidesumi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      materialSum=materialsumi;
      outsideSum=outsidesumi;
   }
   if (mEW->getRank() == 0)
      //      cout << endl 
      //           << "--------------------------------------------------------------\n"
      //           << "Sfile Initialized Node Types: " << endl
      //           << "   Material:        " << materialSum << endl
      //           << endl
      //           << "*Outside Domain:    " << outsideSum << endl
      //           << endl 
      //           << "--------------------------------------------------------------\n"
      //           << endl;
      cout << endl
	   << "sfile command: outside = " << outsideSum << ", material = " << materialSum << endl;

}


//-----------------------------------------------------------------------
void MaterialSfile::read_sfile()
{
  // Timers
  double time_start, time_end;
  double intf_start, intf_end, mat_start, mat_end;
  time_start = MPI_Wtime();

  string rname = "MaterialSfile::read_sfile";

#ifdef USE_HDF5
  // Figure out bounding box in this processor
  float_sw4 xmin=1e38, xmax=-1e38, ymin=1e38, ymax=-1e38, zmin=1e38, zmax=-1e38;
  for( int g=0 ; g < mEW->mNumberOfGrids ; g++ ) {
     float_sw4 h=mEW->mGridSize[g];
     if( xmin > (mEW->m_iStart[g]-1)*h )
        xmin =  (mEW->m_iStart[g]-1)*h;
     if( xmax < (mEW->m_iEnd[g]-1)*h )
        xmax =  (mEW->m_iEnd[g]-1)*h;
     if( ymin > (mEW->m_jStart[g]-1)*h )
        ymin =  (mEW->m_jStart[g]-1)*h;
     if( ymax < (mEW->m_jEnd[g]-1)*h )
        ymax =  (mEW->m_jEnd[g]-1)*h;
     if( mEW->topographyExists() && g >= mEW->mNumberOfCartesianGrids ) {
        int kb=mEW->m_kStart[g], ke=mEW->m_kEnd[g];
        for( int j=mEW->m_jStart[g] ; j <= mEW->m_jEnd[g] ; j++ ) {
           for( int i=mEW->m_iStart[g] ; i <= mEW->m_iEnd[g] ; i++ ) {
              if( zmin > mEW->mZ[g](i,j,kb) )
                  zmin = mEW->mZ[g](i,j,kb);
              if( zmax < mEW->mZ[g](i,j,ke) )
                  zmax = mEW->mZ[g](i,j,ke);
           }
        }
     }
     else {
        if( zmin > (mEW->m_kStart[g]-1)*h + mEW->m_zmin[g] ) 
           zmin = (mEW->m_kStart[g]-1)*h + mEW->m_zmin[g];
        if( zmax < (mEW->m_kEnd[g]-1)*h + mEW->m_zmin[g] )
           zmax = (mEW->m_kEnd[g]-1)*h + mEW->m_zmin[g];
     }
  }
  m_xminloc = xmin;
  m_xmaxloc = xmax;
  m_yminloc = ymin;
  m_ymaxloc = ymax;
  m_zminloc = zmin;
  m_zmaxloc = zmax;

  string fname = m_model_dir + "/" + m_model_file;

  hid_t file_id, dataset_id, datatype_id, h5_dtype, group_id, grid_id, memspace_id, filespace_id, attr_id, plist_id, dxpl;
  int prec;
  double lonlataz[3];
  herr_t ierr;
  hsize_t dim[3];

  // Open file
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  dxpl = H5Pcreate(H5P_DATASET_XFER);
  /* H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE); */
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plist_id);
  if (file_id < 0) {
     cout << "Could not open hdf5 file: " << fname.c_str()<< endl;
     MPI_Abort(MPI_COMM_WORLD, file_id);
  }

  // Read sfile header. Translate each patch into SW4 Cartesian coordinate system
  // Origin longitude, latitude, azimuth
  attr_id = H5Aopen(file_id, "Origin longitude, latitude, azimuth", H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  ierr = H5Aread(attr_id, H5T_NATIVE_DOUBLE, lonlataz);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);

  // ---------- attenuation on file ?
  int att;
  attr_id = H5Aopen(file_id, "Attenuation", H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  ierr = H5Aread(attr_id, H5T_NATIVE_INT, &att);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);

  m_use_attenuation = (att==1);

  // ---------- azimuth on file
  double alpha = lonlataz[2], lon0 = lonlataz[0], lat0 = lonlataz[1];
  CHECK_INPUT( fabs(alpha-mEW->getGridAzimuth()) < 1e-6, "ERROR: sfile azimuth must be equal "
               "to coordinate system azimuth" << " azimuth on sfile = " << alpha << 
               " azimuth of coordinate sytem = " << mEW->getGridAzimuth() );

  // ---------- origin on file
  mEW->computeCartesianCoord( m_x0, m_y0, lon0, lat0 );

  // ---------- number of blocks on file
  attr_id = H5Aopen(file_id, "ngrids", H5P_DEFAULT);
  ierr = H5Aread(attr_id, H5T_NATIVE_INT, &m_npatches);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);

//test
  if (mEW->getRank()==0 && mEW->getVerbosity() >= 2) {
    printf("Sfile header: att=%i\n", att);
    printf("              azimuth=%e, lon0=%e, lat0=%e\n", alpha, lon0, lat0);
    printf("              nblocks=%i\n", m_npatches);
  }
     
  m_hh.resize(m_npatches);
  m_ni.resize(m_npatches);
  m_nj.resize(m_npatches);
  m_nk.resize(m_npatches);

  // ---------- read block headers
  vector<int> ncblock(m_npatches);
  char grid_name[16];

  group_id = H5Gopen(file_id, "Material_model", H5P_DEFAULT);
  ASSERT(group_id >= 0);
  for( int p=0 ; p < m_npatches ; p++ ) {

     sprintf(grid_name, "grid_%d", p);
     grid_id = H5Gopen(group_id, grid_name, H5P_DEFAULT);
     ASSERT(grid_id >= 0);

     // ---------- first part of block header
     double hs;
     attr_id = H5Aopen(grid_id, "Horizontal grid size", H5P_DEFAULT);
     ASSERT(attr_id >= 0);
     ierr = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &hs);
     ASSERT(ierr >= 0);
     H5Aclose(attr_id);

     m_hh[p] = static_cast<float_sw4>(hs);

     attr_id = H5Aopen(grid_id, "Number of components", H5P_DEFAULT);
     ASSERT(attr_id >= 0);
     ierr = H5Aread(attr_id, H5T_NATIVE_INT, &ncblock[p]);
     ASSERT(ierr >= 0);
     H5Aclose(attr_id);

     dataset_id = H5Dopen(grid_id, "Cp", H5P_DEFAULT);
     ASSERT(dataset_id >= 0);

     filespace_id = H5Dget_space(dataset_id);
     H5Sget_simple_extent_dims(filespace_id, dim, NULL);
     m_ni[p]    = (int)dim[0];
     m_nj[p]    = (int)dim[1];
     m_nk[p]    = (int)dim[2];

     H5Sclose(filespace_id);
     H5Dclose(dataset_id);
     H5Gclose(grid_id);

     if (mEW->getRank()==0 && mEW->getVerbosity() >= 2) {
       printf("  header block #%i\n", p);
       printf("  hh=%e\n", m_hh[p]);
       printf("  nc=%i, ni=%i, nj=%i, nk=%i\n", ncblock[p], m_ni[p], m_nj[p], m_nk[p]);
     }
  } // End for each patch

  double min_max_z[2];
  attr_id = H5Aopen(file_id, "Min, max depth", H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  ierr = H5Aread(attr_id, H5T_NATIVE_DOUBLE, min_max_z);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);

  // Intersect local grid with grid on sfile, assume all patches have same x- and y- extent. 
  float_sw4 xminrf = m_x0,    xmaxrf = m_x0+(m_ni[0]-1)*m_hh[0];
  float_sw4 yminrf = m_y0,    ymaxrf = m_y0+(m_nj[0]-1)*m_hh[0];
  float_sw4 zminrf = (float_sw4)min_max_z[0], zmaxrf = (float_sw4)min_max_z[1];

  if( xminrf > m_xminloc )
     m_xminloc = xminrf;
  if( xmaxrf < m_xmaxloc )
     m_xmaxloc = xmaxrf;
  if( yminrf > m_yminloc )
     m_yminloc = yminrf;
  if( ymaxrf < m_ymaxloc )
     m_ymaxloc = ymaxrf;
  if( zminrf > m_zminloc )
     m_zminloc = zminrf;
  if( zmaxrf < m_zmaxloc )
     m_zmaxloc = zmaxrf;
  
  
  mMaterial_rho.resize(m_npatches);
  mMaterial_cp.resize(m_npatches);
  mMaterial_cs.resize(m_npatches);
  if (m_use_attenuation) {
    mMaterial_qp.resize(m_npatches);
    mMaterial_qs.resize(m_npatches);
  }

  m_ifirst.resize(m_npatches);
  m_ilast.resize(m_npatches);
  m_jfirst.resize(m_npatches);
  m_jlast.resize(m_npatches);
  m_kfirst.resize(m_npatches);
  m_klast.resize(m_npatches);

  m_outside = m_xminloc >= m_xmaxloc || m_yminloc >= m_ymaxloc;
  m_isempty.resize(m_npatches);

  if( !m_outside ) {
     // each patch, figure out a box that encloses [xmin,xmax] x [ymin,ymax] x [zmin,zmax]
     for( int p=0 ; p < m_npatches ; p++ ) {
        m_ifirst[p] = static_cast<int>(floor( 1 + (m_xminloc-m_x0)/m_hh[p]));
        m_ilast[p]  = static_cast<int>( ceil(  1 + (m_xmaxloc-m_x0)/m_hh[p]));
        m_jfirst[p] = static_cast<int>(floor( 1 + (m_yminloc-m_y0)/m_hh[p]));
        m_jlast[p]  = static_cast<int>( ceil(  1 + (m_ymaxloc-m_y0)/m_hh[p]));
        m_kfirst[p] = 1;
        m_klast[p]  = m_nk[p];
        // Read all depth for simplicity
        // Limit index ranges to global size limits
        if( m_ifirst[p] < 1 )
           m_ifirst[p] = 1;
        if( m_ilast[p] > m_ni[p] )
           m_ilast[p] = m_ni[p];
        if( m_jfirst[p] < 1 )
           m_jfirst[p] = 1;
        if( m_jlast[p] > m_nj[p] )
           m_jlast[p] = m_nj[p];
        if( m_kfirst[p] < 1 )
           m_kfirst[p] = 1;
        if( m_klast[p] > m_nk[p] )
           m_klast[p] = m_nk[p];

        m_isempty[p] = false;
        if( m_klast[p] < m_kfirst[p] ) {
           m_isempty[p] = true;
           m_ilast[p] = 0;
           m_jlast[p] = 0;
           m_klast[p] = 0;
           m_ifirst[p] = 1;
           m_jfirst[p] = 1;
           m_kfirst[p] = 1;
        }
        /* if (mEW->getRank() == 0 && mEW->getVerbosity() >= 2) { */
        if (mEW->getVerbosity() >= 3) {
           cout << "myRank = " << mEW->getRank() << endl;
           cout << "patch nr " << p << " i " << m_ifirst[p] << " " << m_ilast[p] <<
    	  " j " << m_jfirst[p] << " " << m_jlast[p] << 
    	  " k " << m_kfirst[p] << " " << m_klast[p] << endl;
           cout << "nr components " << ncblock[p] << endl;
           cout << "patch nr global size " << m_ni[p] << " x " << m_nj[p] << " x " << m_nk[p] << endl;
        }
     }
  }
  else {
     for( int p=0 ; p < m_npatches ; p++ ) {
        m_isempty[p] = true;
        m_ilast[p] = 0;
        m_jlast[p] = 0;
        m_klast[p] = 0;
        m_ifirst[p] = 1;
        m_jfirst[p] = 1;
        m_kfirst[p] = 1;
     }
  }
  vector<int> isempty(m_npatches), isemptymin(m_npatches);
  for( int p=0 ; p < m_npatches ; p++ )
     isempty[p] = m_isempty[p];
  MPI_Allreduce( &isempty[0], &isemptymin[0], m_npatches, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
  for( int p=0 ; p < m_npatches ; p++ )
     m_isempty[p] = (isemptymin[p] == 1);

  //Allocate memory
  for( int p=0 ; p < m_npatches ; p++ ) {
     try {
        if( !m_isempty[p] ) {
           mMaterial_rho[p].define(1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
           mMaterial_cp[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
           mMaterial_cs[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
           if (m_use_attenuation) {
             mMaterial_qp[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
             mMaterial_qs[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
           }
        }
     }
     catch( bad_alloc& ba ) {
        cout << "Processor " << mEW->getRank() << " allocation of mMaterial failed." << endl;
        cout << "p= "<< p << " ncblock= " << ncblock[p] << " ifirst,ilast " << m_ifirst[p] << " " << m_ilast[p] <<
           " jfirst,jlast " << m_jfirst[p] << " " << m_jlast[p] <<
           " kfirst,klast " << m_kfirst[p] << " " << m_klast[p] << 
           " Exception= " << ba.what() << endl;
        MPI_Abort(MPI_COMM_WORLD,0);
     }
  }

  bool roworder = true;

  // Read interfaces
  hid_t topo_grp, topo_id;
  mInterface.resize(m_npatches+1);
  char intf_name[32];
  void  *in_data;

  topo_grp = H5Gopen(file_id, "Z_interfaces", H5P_DEFAULT);
  ASSERT(topo_grp >= 0);

  intf_start = MPI_Wtime();
  for (int p = 0; p < m_npatches+1; p++) {
    sprintf(intf_name, "z_values_%d", p);
    dataset_id = H5Dopen(topo_grp, intf_name, H5P_DEFAULT);
    ASSERT(dataset_id >= 0);

    filespace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(filespace_id, dim, NULL);
    H5Sclose(filespace_id);

    mInterface[p].define(1, 1, (int)dim[0], 1, (int)dim[1], 1, 1);

    float  *f_data = new  float[dim[0]*dim[1]];
    double *d_data = new double[dim[0]*dim[1]];

    if (p == 0) {
        // Get precision
        datatype_id = H5Dget_type(dataset_id);
        prec = (int)H5Tget_size(datatype_id);
        H5Tclose(datatype_id);
  
        if (prec == 4) 
            h5_dtype = H5T_NATIVE_FLOAT;
        else if (prec == 8) 
            h5_dtype = H5T_NATIVE_DOUBLE;
    }
    if (prec == 4) 
        in_data = (void*)f_data;
    else if (prec == 8) 
        in_data = (void*)d_data;
  
    if (mEW->getRank() == 0) {
      ierr = H5Dread(dataset_id, h5_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, in_data);
      ASSERT(ierr >= 0);
    }
    MPI_Bcast(in_data, dim[0] * dim[1] * prec, MPI_CHAR, 0, MPI_COMM_WORLD);
    H5Dclose(dataset_id);

    if (prec == 4) { mInterface[p].assign(f_data); }
    else {           mInterface[p].assign(d_data); }

    /* printf("Rank %d: interface %d min = %.2f, max = %.2f\n", */ 
    /*         mEW->getRank(), p,  mInterface[p].minimum(), mInterface[p].maximum()); */

    if (roworder)
      mInterface[p].transposeik();

    delete[] f_data;
    delete[] d_data;
  }
  H5Gclose(topo_grp);
  intf_end = MPI_Wtime();

  hid_t rho_id, cs_id, cp_id, qs_id, qp_id;
  for( int p = 0 ; p < m_npatches ; p++ ) {
    if( !m_isempty[p] ) {
      int global[3]={ m_ni[p], m_nj[p], m_nk[p] };
      int local[3] ={ m_ilast[p]-m_ifirst[p]+1, m_jlast[p]-m_jfirst[p]+1, m_klast[p]-m_kfirst[p]+1 };
      int start[3] ={ m_ifirst[p]-1, m_jfirst[p]-1, m_kfirst[p]-1 };

      /* // debug */
      /* printf("\nRank %d: start (%d, %d, %d), count (%d, %d, %d), global (%d, %d, %d), %d points\n\n", */
      /*         mEW->getRank(), start[0], start[1], start[2], local[0], local[1], local[2], */ 
      /*         global[0], global[1], global[2], mMaterial_rho[p].m_npts); */
      /* fflush(stdout); */

      float  *f_data = new float[mMaterial_rho[p].m_npts];
      double *d_data = new double[mMaterial_rho[p].m_npts];

      if (prec == 4) 
          in_data = (void*)f_data;
      else if (prec == 8) 
          in_data = (void*)d_data;


      sprintf(grid_name, "grid_%d", p);
      grid_id = H5Gopen(group_id, grid_name, H5P_DEFAULT);
      ASSERT(grid_id >= 0);

      rho_id = H5Dopen(grid_id, "Rho", H5P_DEFAULT);
      ASSERT(rho_id >= 0);
      cp_id  = H5Dopen(grid_id, "Cp", H5P_DEFAULT);
      ASSERT(cp_id >= 0);
      cs_id  = H5Dopen(grid_id, "Cs", H5P_DEFAULT);
      ASSERT(cs_id >= 0);
      if (m_use_attenuation) {
        qp_id  = H5Dopen(grid_id, "Qp", H5P_DEFAULT);
        ASSERT(qp_id >= 0);
        qs_id  = H5Dopen(grid_id, "Qs", H5P_DEFAULT);
        ASSERT(qs_id >= 0);
      }

      hsize_t h5_global[3], h5_count[3], h5_start[3];
      for (int i = 0; i < 3; i++) {
          h5_global[i] = (hsize_t)global[i];
          h5_count[i]  = (hsize_t)local[i];
          h5_start[i]  = (hsize_t)start[i];
      }

      memspace_id = H5Screate_simple(3, h5_count, NULL);
      ASSERT(memspace_id >= 0);

      filespace_id = H5Dget_space(rho_id);
      ASSERT(filespace_id >= 0);
      ierr = H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, h5_start, NULL, h5_count, NULL);
      ASSERT(ierr >= 0);

      // Read each var individually 
      ierr = H5Dread(rho_id, h5_dtype, memspace_id, filespace_id, dxpl, in_data);
      ASSERT(ierr >= 0);
      if( prec == 4 ) { mMaterial_rho[p].assign(f_data); }
      else {            mMaterial_rho[p].assign(d_data); }

      ierr = H5Dread(cp_id, h5_dtype, memspace_id, filespace_id, dxpl, in_data);
      ASSERT(ierr >= 0);
      if( prec == 4 ) { mMaterial_cp[p].assign(f_data); }
      else {            mMaterial_cp[p].assign(d_data); }

      ierr = H5Dread(cs_id, h5_dtype, memspace_id, filespace_id, dxpl, in_data);
      ASSERT(ierr >= 0);
      if( prec == 4 ) { mMaterial_cs[p].assign(f_data); }
      else {            mMaterial_cs[p].assign(d_data); }

      if (m_use_attenuation) {
          ierr = H5Dread(qp_id, h5_dtype, memspace_id, filespace_id, dxpl, in_data);
          ASSERT(ierr >= 0);
          if( prec == 4 ) { mMaterial_qp[p].assign(f_data); }
          else {            mMaterial_qp[p].assign(d_data); }

          ierr = H5Dread(qs_id, h5_dtype, memspace_id, filespace_id, dxpl, in_data);
          ASSERT(ierr >= 0);
          if( prec == 4 ) { mMaterial_qs[p].assign(f_data); }
          else {            mMaterial_qs[p].assign(d_data); }
      }

      H5Sclose(memspace_id);
      H5Sclose(filespace_id);
      H5Dclose(rho_id);
      H5Dclose(cp_id);
      H5Dclose(cs_id);
      if (m_use_attenuation) {
          H5Dclose(qp_id);
          H5Dclose(qs_id);
      }
      H5Gclose(grid_id);

      delete[] f_data;
      delete[] d_data;

      if( roworder ) {
         mMaterial_rho[p].transposeik();
         mMaterial_cp[p].transposeik();
         mMaterial_cs[p].transposeik();
         if (m_use_attenuation) {
           mMaterial_qp[p].transposeik();
           mMaterial_qs[p].transposeik();
         }
      }

    } // End if !m_isempty
  } // End for

  H5Pclose(plist_id);
  H5Pclose(dxpl);
  H5Gclose(group_id);
  H5Fclose(file_id);

  /* fill_in_fluids(); */
  // material_check(false);

  time_end = MPI_Wtime();
  if (mEW->getRank() == 0) {
     /* cout << "MaterialSfile::read_sfile, time to read interface: " << intf_end - intf_start << " seconds." << endl; */
     cout << "MaterialSfile::read_sfile, time to read material file: " << time_end - time_start << " seconds." << endl;
  }
  cout.flush();
#endif
}

//-----------------------------------------------------------------------
void MaterialSfile::fill_in_fluids()
{
// Start from p=0
// start from the last (bottom) block and progress upwards
  if( !m_outside ) {
    for( int p=m_npatches-1 ; p >= 0; p-- ) {
      if( !m_isempty[p] ) {
#pragma omp parallel for	 
        for( int j=mMaterial_cs[p].m_jb ; j <= mMaterial_cs[p].m_je ; j++ ) {
           for( int i=mMaterial_cs[p].m_ib ; i <= mMaterial_cs[p].m_ie ; i++ ) {
             int k0 = mMaterial_cs[p].m_kb;
             while( mMaterial_cs[p](1,i,j,k0) < 0 && k0 < mMaterial_cs[p].m_ke )
                k0++;
             // consider the case where the top block is all water. Then k0 = mMaterial[p].m_ke and mMaterial[p](3,i,j,k0)=-999
             // k0 is now the first k with cs > 0.
             if (mMaterial_cs[p](1,i,j,k0) < 0) {
                // get value from block p+1
                if (p<m_npatches-1) {
                   int pd=p+1, id, jd, kd; // index of donor block
                   float_sw4 xm=(i-1)*m_hh[p];
                   float_sw4 ym=(j-1)*m_hh[p];
                   // get closest (id,jd) index on patch pd
                   id = static_cast<int>( 1 + trunc(xm/m_hh[pd]) );
                   jd = static_cast<int>( 1 + trunc(ym/m_hh[pd]) );
                   kd = mMaterial_cs[pd].m_kb; // get value from top of block pd
              
                   if (! (id >= mMaterial_cs[pd].m_ib && id <= mMaterial_cs[pd].m_ie && 
                          jd >= mMaterial_cs[pd].m_jb && jd <= mMaterial_cs[pd].m_je )) {
                      // out of bounds: find nearest interior point
                      if (id < mMaterial_cs[pd].m_ib) id=mMaterial_cs[pd].m_ib;
                      if (id > mMaterial_cs[pd].m_ie) id=mMaterial_cs[pd].m_ie;
                      if (jd < mMaterial_cs[pd].m_jb) jd=mMaterial_cs[pd].m_jb;
                      if (jd > mMaterial_cs[pd].m_je) jd=mMaterial_cs[pd].m_je;
                      
                      printf("WARNING: nearest grid point to (%e,%e) was outside local part of block pd=%i\n"
                       " using id=%i, jd=%i, at (%e, %e)\n", xm, ym, pd, id, jd, (id-1)*m_hh[pd], (jd-1)*m_hh[pd]);

                   }
                   // get values from block 'pd'
                   mMaterial_rho[p](1,i,j,k0)= mMaterial_rho[pd](1,id,jd,kd);
                   mMaterial_cp[p](1,i,j,k0)= mMaterial_cp[pd](1,id,jd,kd);
                   mMaterial_cs[p](1,i,j,k0)= mMaterial_cs[pd](1,id,jd,kd);
                   if (m_use_attenuation) {
                      mMaterial_qp[p](1,i,j,k0)= mMaterial_qp[pd](1,id,jd,kd);
                      mMaterial_qs[p](1,i,j,k0)= mMaterial_qs[pd](1,id,jd,kd);
                   }
                }
                else {
                   printf("ERROR: found undefined material properties in last material block\n"
                          " patch p=%i, i=%i, j=%i, k0=%i\n", p, i, j, k0);
                }
             }

             for( int k=mMaterial_cs[p].m_kb ; k < k0 ; k++ ) {
                mMaterial_rho[p](1,i,j,k) = mMaterial_rho[p](1,i,j,k0);
                mMaterial_cp[p](1,i,j,k)  = mMaterial_cp[p](1,i,j,k0);
                mMaterial_cs[p](1,i,j,k)  = mMaterial_cs[p](1,i,j,k0);
                if( m_use_attenuation ) {
                   mMaterial_qp[p](1,i,j,k) = mMaterial_qp[p](1,i,j,k0);
                   mMaterial_qs[p](1,i,j,k) = mMaterial_qs[p](1,i,j,k0);
                }
             } // End for k

           } // End for i
         } // End for j
       } // End if !m_isempty
     } // End for p 
  } // End if !outside
}

//-----------------------------------------------------------------------
void MaterialSfile::material_check( bool water )
{
   bool printsmallcpcs=false;
   for( int p=1 ; p < m_npatches ; p++ )
   {
      double csmin=1e38,cpmin=1e38,cratmin=1e38,csmax=-1e38,cpmax=-1e38,cratmax=-1e38;
      double rhomin=1e38, rhomax=-1e38;
      for( int k=mMaterial_cs[p].m_kb ; k<= mMaterial_cs[p].m_ke ; k++ )
	 for( int j=mMaterial_cs[p].m_jb ; j<= mMaterial_cs[p].m_je ; j++ )
	    for( int i=mMaterial_cs[p].m_ib ; i<= mMaterial_cs[p].m_ie ; i++ )
	    {
	       if( water || mMaterial_cs[p](1,i,j,k) != -999 )
	       {
		  if( mMaterial_rho[p](1,i,j,k) < rhomin )
		     rhomin = mMaterial_rho[p](1,i,j,k);
		  if( mMaterial_rho[p](1,i,j,k) > rhomax )
		     rhomax = mMaterial_rho[p](1,i,j,k);
		  if( mMaterial_cs[p](1,i,j,k) < csmin )
		     csmin = mMaterial_cs[p](1,i,j,k);
		  if( mMaterial_cs[p](1,i,j,k) > csmax )
		     csmax = mMaterial_cs[p](1,i,j,k);
		  if( mMaterial_cp[p](1,i,j,k) < cpmin )
		     cpmin = mMaterial_cp[p](1,i,j,k);
		  if( mMaterial_cp[p](1,i,j,k) > cpmax )
		     cpmax = mMaterial_cp[p](1,i,j,k);
		  double crat = mMaterial_cp[p](1,i,j,k)/mMaterial_cs[p](1,i,j,k);
		  if( crat < cratmin ) {
		     cratmin = crat;
		     if( printsmallcpcs && crat < 1.41 ) {
			cout << "crat= " << crat << " at " << i << " " <<  j << " " << k << endl;
			cout << " material is " << mMaterial_rho[p](1,i,j,k) << " " << mMaterial_cp[p](1,i,j,k) << " " 
			     << mMaterial_cs[p](1,i,j,k) << " " << mMaterial_qp[p](1,i,j,k) << " " << mMaterial_qs[p](1,i,j,k) << endl;
		     }
		  }
		  if( crat > cratmax )
		     cratmax = crat;
	       }
	    }
      double cmins[4]={csmin,cpmin,cratmin,rhomin}, cmaxs[4]={csmax,cpmax,cratmax,rhomax};
      double cminstot[4], cmaxstot[4];
      MPI_Reduce(cmins, cminstot, 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
      MPI_Reduce(cmaxs, cmaxstot, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
      int myid;
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      if( myid == 0 )
	 //	 if( mEW->getRank()==0 )
      {
	 if( p== 1 && !water )
	    cout << "S-file limits, away from water: " << endl;
	 else if( p== 1 )
	    cout << "S-file limits : " << endl;
	 cout << "  Patch no " << p << " : " << endl;
	 cout << "    cp    min and max " << cminstot[1] << " " << cmaxstot[1] << endl;
	 cout << "    cs    min and max " << cminstot[0] << " " << cmaxstot[0] << endl;
	 cout << "    cp/cs min and max " << cminstot[2] << " " << cmaxstot[2] << endl;
	 cout << "    rho   min and max " << cminstot[3] << " " << cmaxstot[3] << endl;
      }
   }
}


