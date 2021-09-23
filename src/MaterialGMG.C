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
#include "MaterialGMG.h"
#include "Byteswapper.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

using namespace std;


//-----------------------------------------------------------------------
MaterialGMG::MaterialGMG( EW* a_ew, const string a_file, const string a_directory):
   mEW(a_ew),
   m_model_file(a_file),
   m_model_dir(a_directory),
   m_use_attenuation(false)
{
   mCoversAllPoints = false;
   // Check that the depths make sense
   if (a_ew != NULL) {
     m_use_attenuation = a_ew->usingAttenuation();
     read_gmg();
   }
}

//-----------------------------------------------------------------------
MaterialGMG::~MaterialGMG()
{
}

//-----------------------------------------------------------------------
void MaterialGMG::set_material_properties(std::vector<Sarray> & rho, 
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

                if (i == 6 && j == 1 && k == -1 && g == 1) {
                    int ttt=1;
                }

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
                   float_sw4 down_z;
                   // gr is the patch id that has the current sw4 grid point's data
                   // need to use the interface value to determine which patch the current point is in
		   while( gr >= 0 ) {
                     i0 = i1 = static_cast<int>( trunc( 1 + (x-m_x0)/m_hh[gr] ) );
                     j0 = j1 = static_cast<int>( trunc( 1 + (y-m_y0)/m_hh[gr] ) );
		     if( i0 <= m_ifirst[gr] )  i0 = m_ifirst[gr];
		     if( i0 >= m_ilast[gr]-1 ) i0 = m_ilast[gr]-1;
		     if( j0 <= m_jfirst[gr] )  j0 = m_jfirst[gr];
		     if( j0 >= m_jlast[gr]-1 ) j0 = m_jlast[gr]-1;
		     if( i1 <= m_ifirst[gr] )  i1 = m_ifirst[gr];
		     if( i1 >= m_ilast[gr]-1 ) i1 = m_ilast[gr]-1;
		     if( j1 <= m_jfirst[gr] )  j1 = m_jfirst[gr];
		     if( j1 >= m_jlast[gr]-1 ) j1 = m_jlast[gr]-1;

                     down_z = mInterface[gr+1](1, i0, j0, 1);

                     // Adjust the index if upper and lower interface have different dimension
                     if (ist[gr] > 1) { i1 = i1*ist[gr]-1; }
                     if (jst[gr] > 1) { j1 = j1*jst[gr]-1; }
                     z0 = mInterface[gr](1, i1, j1, 1);

                     if (gr == 0 || z > z0) break;
		     gr--;
                   }

                   if (z > down_z)
                       z = down_z;
                   if (z < z0)     
                       z = z0;

                   // Update the current vertical grid height and z-base with the GMG curvilinear grid
                   if (gr == m_npatches-1)
                     hv = (down_z-z0) / (m_nk[gr]-1);
                   else
                     hv = m_hv[gr];

                   // we are using curvilinear grid in GMG 
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
		      //                       printf("g=%d, sw4 (%d, %d, %d), mat (%d, %d, %d) wghx = %.2f\n", gr, i, j, k, i0, j0, k0, wghx);
                       if (wghx > 1) wghx = 1;
                       if (wghx < 0) wghx = 0;
                   }

                   if (wghy > 1 || wghy < 0) { 
		      //                       printf("g=%d, sw4 (%d, %d, %d), mat (%d, %d, %d) wghy = %.2f\n", gr, i, j, k, i0, j0, k0, wghy);
                       if (wghy > 1) wghy = 1;
                       if (wghy < 0) wghy = 0;
                   }

                   if (wghz > 1 || wghz < 0) { 
                       if (wghz > 1.001 || wghz < -0.001) 
			  //                         printf("g=%d, sw4 (%d, %d, %d), mat (%d, %d, %d) wghz = %.2f, z=%.2f, z0=%.2f\n", 
			  //                                 gr, i, j, k, i0, j0, k0, wghz, z, z0);
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

                   /* if (rho[g](i,j,k) < 1500) { */
                   /*   printf("Rank %d, rho[%d](%d, %d, %d)=%.2f\n", mEW->getRank(), g, i, j, k, rho[g](i,j,k)); */
                   /*   ASSERT(0); */
                   /* } */

                   cp[g](i, j, k)  = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_cp[gr](1,i0,j0,k0) + 
                                     wghx*mMaterial_cp[gr](1,i0+1,j0,k0) ) +
                                     wghy*(     (1-wghx)*mMaterial_cp[gr](1,i0,j0+1,k0) + 
                                     wghx*mMaterial_cp[gr](1,i0+1,j0+1,k0) ) ) + 
                                     wghz*(  (1-wghy)*( (1-wghx)*mMaterial_cp[gr](1,i0,j0,k0+1) + 
                                     wghx*mMaterial_cp[gr](1,i0+1,j0,k0+1) ) +
                                     wghy*(     (1-wghx)*mMaterial_cp[gr](1,i0,j0+1,k0+1)+ 
                                     wghx*mMaterial_cp[gr](1,i0+1,j0+1,k0+1) ) );
       
                   /* if (cp[g](i,j,k) < 700) { */
                     /* printf("Rank %d, cp[%d](%d, %d, %d)=%.2f\n", mEW->getRank(), g, i, j, k, cp[g](i,j,k)); */
                     /* ASSERT(0); */
                   /* } */

                   cs[g](i, j, k)  = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_cs[gr](1,i0,j0,k0) + 
                                     wghx*mMaterial_cs[gr](1,i0+1,j0,k0) ) +
                                     wghy*(    (1-wghx)*mMaterial_cs[gr](1,i0,j0+1,k0) + 
                                     wghx*mMaterial_cs[gr](1,i0+1,j0+1,k0) ) ) + 
                                     wghz*(   (1-wghy)*( (1-wghx)*mMaterial_cs[gr](1,i0,j0,k0+1) + 
                                     wghx*mMaterial_cs[gr](1,i0+1,j0,k0+1) ) +
                                     wghy*(   (1-wghx)*mMaterial_cs[gr](1,i0,j0+1,k0+1)+ 
                                     wghx*mMaterial_cs[gr](1,i0+1,j0+1,k0+1) ) );

                   /* if (cs[g](i,j,k) > 4000) { */
                   /*   printf("Rank %d, cs[%d](%d, %d, %d)=%.2f\n", mEW->getRank(), g, i, j, k, cs[g](i,j,k)); */
                     /* ASSERT(0); */
                   /* } */

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
      MPI_Reduce(&material, &materialSum, 1, MPI_LONG, MPI_SUM, 0, mEW->m_1d_communicator );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG, MPI_SUM, 0, mEW->m_1d_communicator );
   }
   else if( sizeof(size_t) == mpisizelonglong ) {
      MPI_Reduce(&material, &materialSum, 1, MPI_LONG_LONG, MPI_SUM, 0, mEW->m_1d_communicator );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG_LONG, MPI_SUM, 0, mEW->m_1d_communicator );
   }
   else if( sizeof(size_t) == mpisizeint ) {
      MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0, mEW->m_1d_communicator );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_INT, MPI_SUM, 0, mEW->m_1d_communicator );
   }
   else {
      int materialsumi, outsidesumi, materiali=material, outsidei=outside;
      MPI_Reduce(&materiali, &materialsumi, 1, MPI_INT, MPI_SUM, 0, mEW->m_1d_communicator );
      MPI_Reduce(&outsidei,   &outsidesumi, 1, MPI_INT, MPI_SUM, 0, mEW->m_1d_communicator );
      materialSum=materialsumi;
      outsideSum=outsidesumi;
   }
   if (mEW->getRank() == 0)
      //      cout << endl 
      //           << "--------------------------------------------------------------\n"
      //           << "GMG Initialized Node Types: " << endl
      //           << "   Material:        " << materialSum << endl
      //           << endl
      //           << "*Outside Domain:    " << outsideSum << endl
      //           << endl 
      //           << "--------------------------------------------------------------\n"
      //           << endl;
      cout << endl
	   << "gmg command: outside = " << outsideSum << ", material = " << materialSum << endl;

}

#ifdef USE_HDF5
static void read_hdf5_attr(hid_t loc, hid_t dtype, char *name, void* data)
{
  hid_t attr_id;
  int ierr;
  attr_id = H5Aopen(loc, name, H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  ierr = H5Aread(attr_id, dtype, data);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);
}
#endif

//-----------------------------------------------------------------------
void MaterialGMG::read_gmg()
{
  // Timers
  double time_start, time_end;
  double intf_start, intf_end, mat_start, mat_end;
  time_start = MPI_Wtime();

  string rname = "MaterialGMG::read_gmg";

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
  m_use_attenuation = true;

  string fname = m_model_dir + "/" + m_model_file;

  hid_t file_id, dataset_id, group_id, filespace_id, topo_grp, topo_id, fapl;
  double az, origin_x, origin_y, hh, hv, alpha, lon0, lat0, max_z, min_z, ztop[10];
  herr_t ierr;
  hsize_t dims[4];
  float *topo_data;
  char grid_name[32];

  // Fixed for GMG grids
  m_npatches = 4;

  m_hh.resize(m_npatches);
  m_hv.resize(m_npatches);
  m_ni.resize(m_npatches);
  m_nj.resize(m_npatches);
  m_nk.resize(m_npatches);
  vector<int> ncblock(m_npatches);

  m_hv[0] = 25;
  m_hv[1] = 50;
  m_hv[2] = 125;
  m_hv[3] = 250;

  if (mEW->getRank() == 0) {
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL);
    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, fapl);
    if (file_id < 0) {
       cout << "Could not open hdf5 file: " << fname.c_str()<< endl;
       MPI_Abort(MPI_COMM_WORLD, file_id);
    }
    H5Pclose(fapl);

    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "origin_x", &origin_x);
    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "origin_y", &origin_y);
    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "y_azimuth", &az);
    read_hdf5_attr(file_id, H5T_IEEE_F64LE, "dim_z", &max_z);

    group_id = H5Gopen(file_id, "blocks", H5P_DEFAULT);
    ASSERT(group_id >= 0);

    int factor = 1;
    for( int p = 0 ; p < m_npatches ; p++ ) {
      ncblock[p] = 7;
      sprintf(grid_name, "vres%dm", (int)m_hv[p]);
      dataset_id = H5Dopen(group_id, grid_name, H5P_DEFAULT);
      ASSERT(dataset_id >= 0);

      filespace_id = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(filespace_id, dims, NULL);
      m_nj[p]    = (int)dims[0];
      m_ni[p]    = (int)dims[1];
      m_nk[p]    = (int)dims[2];

      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "z_top", &ztop[p]);
      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "resolution_horiz", &m_hh[p]);
      read_hdf5_attr(dataset_id, H5T_IEEE_F64LE, "resolution_vert", &hv);
      ASSERT((int)hv == (int)m_hv[p]);

      H5Sclose(filespace_id);
      H5Dclose(dataset_id);

      if (mEW->getVerbosity() >= 2) {
        printf("  header block #%i\n", p);
        printf("  hh=%e\n", m_hh[p]);
        printf("  nc=%i, ni=%i, nj=%i, nk=%i\n", ncblock[p], m_ni[p], m_nj[p], m_nk[p]);
      }
      factor *= 2;
    } // End for each patch

    topo_grp = H5Gopen(file_id, "surfaces", H5P_DEFAULT);
    ASSERT(topo_grp >= 0);

    dataset_id = H5Dopen(topo_grp, "top_surface", H5P_DEFAULT);
    ASSERT(dataset_id >= 0);

    filespace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(filespace_id, dims, NULL);

    topo_data = new float[dims[0]*dims[1]]();
    ierr = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, filespace_id, H5P_DEFAULT, topo_data);
    ASSERT(ierr >= 0);
  
    H5Sclose(filespace_id);
    H5Dclose(dataset_id);
    H5Gclose(group_id);
    H5Gclose(topo_grp);
    H5Fclose(file_id);

    min_z = -1e10;
    for (int i = 0; i < dims[0]*dims[1]; i++)
      if (topo_data[i] > min_z)
        min_z = topo_data[i];

  } // End rank==0

  MPI_Bcast(&origin_x, 1,               MPI_DOUBLE,    0, mEW->m_1d_communicator);
  MPI_Bcast(&origin_y, 1,               MPI_DOUBLE,    0, mEW->m_1d_communicator);
  MPI_Bcast(&az,       1,               MPI_DOUBLE,    0, mEW->m_1d_communicator);
  MPI_Bcast(&max_z,    1,               MPI_DOUBLE,    0, mEW->m_1d_communicator);
  MPI_Bcast(&min_z,    1,               MPI_DOUBLE,    0, mEW->m_1d_communicator);
  MPI_Bcast(dims,      3,               MPI_LONG_LONG, 0, mEW->m_1d_communicator );
  MPI_Bcast(&m_hh[0],  m_npatches,      MPI_DOUBLE,    0, mEW->m_1d_communicator);
  MPI_Bcast(&m_ni[0],  m_npatches,      MPI_INT,       0, mEW->m_1d_communicator);
  MPI_Bcast(&m_nj[0],  m_npatches,      MPI_INT,       0, mEW->m_1d_communicator);
  MPI_Bcast(&m_nk[0],  m_npatches,      MPI_INT,       0, mEW->m_1d_communicator);
  MPI_Bcast(ztop,      m_npatches,      MPI_DOUBLE,    0, mEW->m_1d_communicator);

  if (mEW->getRank() != 0)
    topo_data = new float[dims[0]*dims[1]];

  MPI_Bcast(topo_data, dims[0]*dims[1], MPI_FLOAT,     0, mEW->m_1d_communicator);

  alpha = az - 180.0;

  // Tang: fix origin for GMG data
  lon0 = -122.562;
  lat0 = 39.1745;

  CHECK_INPUT( fabs(alpha-mEW->getGridAzimuth()) < 1e-6, "ERROR: gmg azimuth must be equal "
               "to coordinate system azimuth" << " azimuth on gmg = " << alpha << 
               " azimuth of coordinate sytem = " << mEW->getGridAzimuth() );

  // ---------- origin on file
  mEW->computeCartesianCoord( m_x0, m_y0, lon0, lat0 );
  // Convert to actual x0, y0 from GMG 
  m_x0 += 3600;
  m_y0 += 1300;

  if (mEW->getRank()==0 && mEW->getVerbosity() >= 2) {
    printf("GMG header: \n");
    printf("    azimuth=%e, lon0=%e, lat0=%e\n", alpha, lon0, lat0);
    printf("    nblocks=%i\n", m_npatches);
  }
     
  // Intersect local grid with grid on GMG, assume all patches have same x- and y- extent. 
  float_sw4 xminrf = m_x0,    xmaxrf = m_x0+(m_ni[0]-1)*m_hh[0];
  float_sw4 yminrf = m_y0,    ymaxrf = m_y0+(m_nj[0]-1)*m_hh[0];
  float_sw4 zminrf = (float_sw4)-min_z, zmaxrf = (float_sw4)max_z;

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
  MPI_Allreduce( &isempty[0], &isemptymin[0], m_npatches, MPI_INT, MPI_MIN, mEW->m_1d_communicator );
  for( int p=0 ; p < m_npatches ; p++ )
     m_isempty[p] = (isemptymin[p] == 1);

  // Read interfaces
  mInterface.resize(m_npatches+1);
  int factor = 1, nitop = dims[0], njtop = dims[1];
  for (int p = 0; p < m_npatches+1; p++) {
    float *in_data = new float[nitop * njtop]();

    // Convert to SW4 CRS
    for (int i = 0; i < nitop; i++) {
        for (int j = 0; j < njtop; j++) {
          if (p != m_npatches)
            in_data[(njtop-j-1)*nitop + nitop-i-1] = -topo_data[i*njtop*factor + j*factor] - ztop[p];
          else
            in_data[(njtop-j-1)*nitop + nitop-i-1] = max_z;
        }
    }

    mInterface[p].define(1, 1, njtop, 1, nitop, 1, 1);
    mInterface[p].assign(in_data);
    mInterface[p].transposeik();

    /* fprintf(stderr, "p=%d, factor=%d, mInterface %d %d\n", p, factor, njtop, nitop); */

    if (p > 0 && p < m_npatches) {
      factor = m_hh[p]/m_hh[0];
      nitop = (dims[0]-1)/factor + 1;
      njtop = (dims[1]-1)/factor + 1;
    }
    delete[] in_data;
  } // End for patch

  delete[] topo_data;

  const int nc = 7;
  for( int p = 0 ; p < m_npatches ; p++ ) {

    int global[4]={ m_ni[p], m_nj[p], m_nk[p], 0 };
    int local[4] ={ m_ilast[p]-m_ifirst[p]+1, m_jlast[p]-m_jfirst[p]+1, m_klast[p]-m_kfirst[p]+1, nc };
    int start[4] ={ m_ifirst[p]-1, m_jfirst[p]-1, m_kfirst[p]-1, 0 };

    /* fprintf(stderr, "\nRank %d: start (%d, %d, %d), count (%d, %d, %d), global (%d, %d, %d), hv=%d\n", */
    /*                 mEW->getRank(), start[0], start[1], start[2], local[0], local[1], local[2], */ 
    /*                 global[0], global[1], global[2], (int)m_hv[p]); */

    float *in_data;
    // Rank 0 read and broadcast data
    if (mEW->getRank() == 0) {
      if (p == 0) {
        fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL);
  
        file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, fapl);
        if (file_id < 0) {
           cout << "Could not open hdf5 file: " << fname.c_str()<< endl;
           MPI_Abort(MPI_COMM_WORLD, file_id);
        }
        H5Pclose(fapl);
  
        group_id = H5Gopen(file_id, "blocks", H5P_DEFAULT);
        ASSERT(group_id >= 0);
      }

      sprintf(grid_name, "vres%dm", (int)m_hv[p]);
      dataset_id = H5Dopen(group_id, grid_name, H5P_DEFAULT);
      ASSERT(dataset_id >= 0);
  
      filespace_id = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(filespace_id, dims, NULL);
      // Make sure number of components is as assumed
      ASSERT(dims[3] == nc);
  
      MPI_Bcast(dims, 4, MPI_LONG_LONG, 0, mEW->m_1d_communicator);
  
      in_data  = new float[dims[0]*dims[1]*dims[2]*dims[3]]();
      // Read all var
      ierr = H5Dread(dataset_id, H5T_IEEE_F32LE, H5S_ALL, filespace_id, H5P_DEFAULT, in_data);
      ASSERT(ierr >= 0);
      H5Sclose(filespace_id);
      H5Dclose(dataset_id);

      if (p == m_npatches-1) {
        H5Gclose(group_id);
        H5Fclose(file_id);
      }
  
      MPI_Bcast(in_data, dims[0]*dims[1]*dims[2]*dims[3], MPI_FLOAT, 0, mEW->m_1d_communicator );
    }
    else {
      MPI_Bcast(dims, 4, MPI_LONG_LONG, 0, mEW->m_1d_communicator);
      in_data  = new float[dims[0]*dims[1]*dims[2]*dims[3]]();
      MPI_Bcast(in_data, dims[0]*dims[1]*dims[2]*dims[3], MPI_FLOAT, 0, mEW->m_1d_communicator );
    }

    if( !m_isempty[p] ) {
      // Allocate memory
      try {
        mMaterial_rho[p].define(1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
        mMaterial_cp[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
        mMaterial_cs[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
        if (m_use_attenuation) {
          mMaterial_qp[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
          mMaterial_qs[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
        }
      }
      catch( bad_alloc& ba ) {
         cout << "Processor " << mEW->getRank() << " allocation of mMaterial failed." << endl;
         cout << "p= "<< p << " ncblock= " << ncblock[p] << " ifirst,ilast " << m_ifirst[p] << " " << m_ilast[p] <<
            " jfirst,jlast " << m_jfirst[p] << " " << m_jlast[p] <<
            " kfirst,klast " << m_kfirst[p] << " " << m_klast[p] << 
            " Exception= " << ba.what() << endl;
         MPI_Abort(MPI_COMM_WORLD, 0);
      }

      float *data   = new float[dims[0]*dims[1]]();
      float *f_data = new float[local[0]*local[1]]();
      Sarray hslice;
      hslice.define(1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],1,1);

      int idx1, idx2;
      for (int c = 0; c < 5; c++) {
        for (int k = 0; k < dims[2]; k++) {

          // Convert to SW4 CRS
          for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
              idx1 = (dims[1]-j-1)*dims[0] + dims[0]-i-1;
              idx2 = (i*dims[1]*dims[2] + j*dims[2] + k) * nc + c;
              if (idx1 > dims[0]*dims[1]) {
                  fprintf(stderr, "Rank %d, patch %d, 1 out %d / %d\n", mEW->getRank(), p, idx1, dims[0]*dims[1]);
                  ASSERT(0);
              }
              if (idx2 > dims[0]*dims[1]*dims[2]*dims[3]) {
                  fprintf(stderr, "Rank %d, patch %d, 2 out %d / %d\n", mEW->getRank(), p, idx2, dims[0]*dims[1]*dims[2]*dims[3]);
                  ASSERT(0);
              }
              data[idx1] = (float_sw4)in_data[idx2];

            } // end j
          } // end i

          for (int i = 0; i < local[0]; i++) {
            for (int j = 0; j < local[1]; j++) {
              idx1 = i*local[1] + j;
              idx2 = (i+m_ifirst[p])*dims[0] + (j+m_jfirst[p]);
#ifdef BZ_DEBUG
              if (idx1 > local[0]*local[1]) {
                  fprintf(stderr, "Rank %d, patch %d, 1 out %d / %d\n", mEW->getRank(), p, idx1, local[0]*local[1]);
                  ASSERT(0);
              }
              if (idx2 > dims[0]*dims[1]) {
                  fprintf(stderr, "Rank %d, patch %d, 2 out %d / %d\n", mEW->getRank(), p, idx2, dims[0]*dims[1]);
                  ASSERT(0);
              }
#endif
              f_data[idx1] = data[idx2];
            }
          }
          hslice.assign(f_data);
          hslice.transposeik();

#pragma omp parallel for	 
          for (int ii = m_ifirst[p]; ii <= m_ilast[p]; ii++) {
            for (int jj = m_jfirst[p]; jj <= m_jlast[p]; jj++) {
              switch (c) {
                case 0:
                  mMaterial_rho[p](1, ii, jj, m_kfirst[p]+k) = hslice(1, ii, jj, 1);
                  /* if (ii == m_ilast[p] && jj==m_jlast[p]) { */
                  /*   if (k == 0) { */
                  /*     fprintf(stderr, "p=%d: first %d, %d, %d,  last %d, %d, %d\n", */ 
                  /*             p, m_ifirst[p],  m_jfirst[p], m_kfirst[p], m_ilast[p], m_jlast[p], m_klast[p]); */
                  /*     fprintf(stderr, "p=%d: dims %d, %d, %d", p, dims[0], dims[1], dims[2]); */
                  /*   } */
                  /*   fprintf(stderr, "k=%d, rho=%.2f\n", k, mMaterial_rho[p](1, ii, jj, m_kfirst[p]+k)); */
                  /* } */
                  break;
                case 1:
                  mMaterial_cp[p](1, ii, jj, m_kfirst[p]+k) = hslice(1, ii, jj, 1);
                  /* if (ii == m_ilast[p] && jj==m_jlast[p]) { */
                  /*   if (k == 0) */
                  /*     fprintf(stderr, "p=%d: first %d, %d, %d,  last %d, %d, %d\n", */ 
                  /*             p, m_ifirst[p],  m_jfirst[p], m_kfirst[p], m_ilast[p], m_jlast[p], m_klast[p]); */
                  /*   fprintf(stderr, "k=%d, cp=%.2f\n", k, mMaterial_cp[p](1, ii, jj, m_kfirst[p]+k)); */
                  /* } */
                  break;
                case 2:
                  mMaterial_cs[p](1, ii, jj, m_kfirst[p]+k) = hslice(1, ii, jj, 1);
                  /* if (ii == m_ilast[p] && jj==m_jlast[p]) { */
                  /*   if (k == 0) */
                  /*     fprintf(stderr, "p=%d: first %d, %d, %d,  last %d, %d, %d\n", */ 
                  /*             p, m_ifirst[p],  m_jfirst[p], m_kfirst[p], m_ilast[p], m_jlast[p], m_klast[p]); */
                  /*   fprintf(stderr, "k=%d, cs=%.2f\n", k, mMaterial_cs[p](1, ii, jj, m_kfirst[p]+k)); */
                  /* } */
                  break;
                case 3:
                  mMaterial_qp[p](1, ii, jj, m_kfirst[p]+k) = hslice(1, ii, jj, 1);
                  /* if (ii == m_ilast[p] && jj==m_jlast[p]) { */
                  /*   if (k == 0) */
                  /*     fprintf(stderr, "p=%d: first %d, %d, %d,  last %d, %d, %d\n", */ 
                  /*             p, m_ifirst[p],  m_jfirst[p], m_kfirst[p], m_ilast[p], m_jlast[p], m_klast[p]); */
                  /*   fprintf(stderr, "k=%d, qp=%.2f\n", k, mMaterial_qp[p](1, ii, jj, m_kfirst[p]+k)); */
                  /* } */
                  break;
                case 4:
                  mMaterial_qs[p](1, ii, jj, m_kfirst[p]+k) = hslice(1, ii, jj, 1);
                  /* if (ii == m_ilast[p] && jj==m_jlast[p]) { */
                  /*   if (k == 0) */
                  /*     fprintf(stderr, "p=%d: first %d, %d, %d,  last %d, %d, %d\n", */ 
                  /*             p, m_ifirst[p],  m_jfirst[p], m_kfirst[p], m_ilast[p], m_jlast[p], m_klast[p]); */
                  /*   fprintf(stderr, "k=%d, qs=%.2f\n", k, mMaterial_qs[p](1, ii, jj, m_kfirst[p]+k)); */
                  /* } */
                  break;
                default:
                  ASSERT(0);
                  break;
              }
            }
          }

        } // End k
      } // End c

      delete[] in_data;
      delete[] f_data;
      delete[] data;

    } // End if !m_isempty
  } // End for

  fill_in_fluids();
#ifdef BZ_DEBUG
  material_check(false);
#endif

  time_end = MPI_Wtime();
  if (mEW->getRank() == 0) {
     cout << "MaterialGMG::read_gmg, time to read material file: " << time_end - time_start << " seconds." << endl;
  }
  cout.flush();
#endif
}

//-----------------------------------------------------------------------
void MaterialGMG::fill_in_fluids()
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
             // consider the case where the top block is all water. Then k0 = mMaterial[p].m_ke and mMaterial[p](3,i,j,k0)<0
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
void MaterialGMG::material_check( bool water )
{
   bool printsmallcpcs=false;
   for( int p=0 ; p < m_npatches ; p++ )
   {
      double csmin=1e38,cpmin=1e38,cratmin=1e38,csmax=-1e38,cpmax=-1e38,cratmax=-1e38;
      double rhomin=1e38, rhomax=-1e38;
      for( int k=mMaterial_cs[p].m_kb ; k<= mMaterial_cs[p].m_ke ; k++ )
	 for( int j=mMaterial_cs[p].m_jb ; j<= mMaterial_cs[p].m_je ; j++ )
	    for( int i=mMaterial_cs[p].m_ib ; i<= mMaterial_cs[p].m_ie ; i++ )
	    {
	       if( water || mMaterial_cs[p](1,i,j,k) > 0 )
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
      MPI_Reduce(cmins, cminstot, 4, MPI_DOUBLE, MPI_MIN, 0, mEW->m_1d_communicator );
      MPI_Reduce(cmaxs, cmaxstot, 4, MPI_DOUBLE, MPI_MAX, 0, mEW->m_1d_communicator );
      int myid;
      MPI_Comm_rank(mEW->m_1d_communicator,&myid);
      if( myid == 0 )
	 //	 if( mEW->getRank()==0 )
      {
	 if( p== 1 && !water )
	    cout << "GMG-file limits, away from water: " << endl;
	 else if( p== 1 )
	    cout << "GMG-file limits : " << endl;
	 cout << "  Patch no " << p << " : " << endl;
	 cout << "    cp    min and max " << cminstot[1] << " " << cmaxstot[1] << endl;
	 cout << "    cs    min and max " << cminstot[0] << " " << cmaxstot[0] << endl;
	 cout << "    cp/cs min and max " << cminstot[2] << " " << cmaxstot[2] << endl;
	 cout << "    rho   min and max " << cminstot[3] << " " << cmaxstot[3] << endl;
      }
   }
}


