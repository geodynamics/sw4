//-*-c++-*-
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
#include "mpi.h"

#include "EW.h"
#include "GridGenerator.h"

// #include "F77_FUNC.h"
// extern "C" {
//    void metric( int*, int*, int*, int*, int*, int*,
// 				 double*, double*, double*, double*, double*, int * );
//    void gridinfo( int*, int*, int*, int*, int*, int*,
// 				     double*, double*, double*, double* );
//    void metricexgh( int*, int*, int*, int*, int*, int*, int*, int*, int*,
// 					 double*, double*, double*, double*, double*,
// 					 int*, double*, double*, double*, double*,
// 					 double*, double*, double* );
//    void meterr4c( int*, int*, int*, int*, int*, int*,
// 		  double*, double*, double*, double*, double*, double*,
// 		  int*, int*, int*, int*, int*, int*, double* );
// }
#define SQR(x) ((x)*(x))


//-----------------------------------------------------------------------
void EW::smoothTopography(int maxIter)
{
//   int myRank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

// tmp
//   if (proc_zero())
//   {
//     cout << "***inside smoothTopography***  maxIter= " << maxIter << endl;
//   }
  
  float_sw4 rf=0.2; // rf<0.25 for stability
  //  int topLevel = mNumberOfGrids-1;
  int iter;

  copy_topo_to_topogridext();

  int imin = mTopoGridExt.m_ib;
  int imax = mTopoGridExt.m_ie;
  int jmin = mTopoGridExt.m_jb;
  int jmax = mTopoGridExt.m_je;
  
// temporary storage
  Sarray tmp;
  tmp.define(imin,imax,jmin,jmax,1,1);

// Laplacian filter
  for (iter=0; iter < maxIter; iter++)
  {
#pragma omp parallel for
    for (int i = imin+1; i <= imax-1; ++i)
      for (int j = jmin+1; j <= jmax-1; ++j)
      {
	tmp(i,j,1) = mTopoGridExt(i,j,1) + rf*(mTopoGridExt(i+1,j,1) + mTopoGridExt(i-1,j,1) + mTopoGridExt(i,j+1,1) + mTopoGridExt(i,j-1,1) - 4.*mTopoGridExt(i,j,1));
      }

// Neumann boundary conditions
#pragma omp parallel for
    for (int j = jmin+1; j <= jmax-1; ++j)
    {
      int i = imin;
      tmp(i,j,1) = tmp(i+1,j,1);
      i = imax;
      tmp(i,j,1) = tmp(i-1,j,1);
    }

#pragma omp parallel for
    for (int i = imin+1; i <= imax-1; ++i)
    {
      int j = jmin;
      tmp(i,j,1) = tmp(i,j+1,1);
      j = jmax;
      tmp(i,j,1) = tmp(i,j-1,1);
    }
// Corners
    int i = imin;
    int j = jmin;
    tmp(i,j,1) = tmp(i+1,j+1,1);

    i = imax;
    j = jmin;
    tmp(i,j,1) = tmp(i-1,j+1,1);

    i = imin;
    j = jmax;
    tmp(i,j,1) = tmp(i+1,j-1,1);
    
    i = imax;
    j = jmax;
    tmp(i,j,1) = tmp(i-1,j-1,1);

    communicate_array_2d_ext( tmp );

// update solution
#pragma omp parallel for
    for (int i = imin; i <= imax ; ++i)
      for (int j = jmin; j <= jmax ; ++j)
	 mTopoGridExt(i,j,1) = tmp(i,j,1);
  }// end for iter
}



//-----------------------------------------------------------------------
void EW::metric_derivatives_test()
{
   // Assumes mMetric and mJ have been computed by numerical differentiation
   // This function computes corresponding expressions by analytical differentiation

   int g=mNumberOfGrids-1;
   Sarray metex(mMetric[g]), jacex(mJ[g]); // TEMPORARY FIX

   int Bx=m_iStart[g];
   int By=m_jStart[g];
   int Bz=m_kStart[g];
   int Nx=m_iEnd[g];
   int Ny=m_jEnd[g];
   int Nz=m_kEnd[g];

   //   int nxg = m_global_nx[g];
   //   int nyg = m_global_ny[g];
   //   int nzg = m_global_nz[g];
   float_sw4 h= mGridSize[g];   
   //   float_sw4 zmax = m_zmin[g-1] - (nzg-1)*h*(1-m_zetaBreak);

//FTNC   if( m_croutines )
   m_gridGenerator->exact_metric( this, g, metex, jacex );
   //      metricexgh_ci( Bx, Nx, By, Ny, Bz, Nz, nzg, mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr(),
   //				    metex.c_ptr(), jacex.c_ptr(), m_grid_interpolation_order, m_zetaBreak, zmax, 
   //				    m_GaussianAmp, m_GaussianXc, m_GaussianYc, m_GaussianLx, m_GaussianLy ); 
//FTNC   else
//FTNC      metricexgh( &Bx, &Nx, &By, &Ny, &Bz, &Nz, &nxg, &nyg, &nzg, mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr(),
//FTNC				    metex.c_ptr(), jacex.c_ptr(), &m_grid_interpolation_order, &m_zetaBreak, &zmax, 
//FTNC				    &m_GaussianAmp, &m_GaussianXc, &m_GaussianYc, &m_GaussianLx, &m_GaussianLy ); 
   communicate_array( metex, mNumberOfGrids-1 );
   communicate_array( jacex, mNumberOfGrids-1 );

   float_sw4 li[5], l2[5];
   int imin = m_iStartInt[g];
   int imax = m_iEndInt[g];
   int jmin = m_jStartInt[g];
   int jmax = m_jEndInt[g];
   int kmin = m_kStartInt[g];
   int kmax = m_kEndInt[g];

//FTNC   if( m_croutines )
      meterr4c_ci( Bx, Nx, By, Ny, Bz, Nz, mMetric[g].c_ptr(), metex.c_ptr(), mJ[g].c_ptr(),
		   jacex.c_ptr(), li, l2, imin, imax, jmin, jmax, kmin, kmax, h );
//FTNC   else
//FTNC      meterr4c( &Bx, &Nx, &By, &Ny, &Bz, &Nz, mMetric[g].c_ptr(), metex.c_ptr(), mJ[g].c_ptr(),
//FTNC		jacex.c_ptr(), li, l2, &imin, &imax, &jmin, &jmax, &kmin, &kmax, &h );

   float_sw4 tmp[5];
   for( int c=0 ; c < 5 ;c++ )
      tmp[c] =li[c];
   MPI_Allreduce( tmp, li, 5, m_mpifloat, MPI_MAX, m_cartesian_communicator);
   for( int c=0 ; c < 5 ;c++ )
      tmp[c] =l2[c];
   MPI_Allreduce( tmp, l2, 5, m_mpifloat, MPI_SUM, m_cartesian_communicator);
   for( int c=0 ; c < 5 ;c++ )
      l2[c] = sqrt(l2[c]);
   if( proc_zero() )
   {
      cout << "Errors in metric, max norm and L2 norm \n";
      for( int c=0 ; c < 4 ; c++ )
	 cout << " " << li[c] << " " << l2[c] << endl;
      cout << "Error in Jacobian, max norm and L2 norm \n";
      cout << " " << li[4] << " " << l2[4] << endl;
   }
}

//-----------------------------------------------------------------------
void EW::copy_topo_to_topogridext()
{
   if( topographyExists() )
   {
      int gTop = mNumberOfGrids-1;
// copy raw topography
#pragma omp parallel for
      for (int i = m_iStart[gTop]; i <= m_iEnd[gTop]; ++i)
	 for (int j = m_jStart[gTop]; j <= m_jEnd[gTop]; ++j)
	    mTopoGridExt(i,j,1) = mTopo(i,j,1);

      int imin = mTopoGridExt.m_ib;
      int imax = mTopoGridExt.m_ie;
      int jmin = mTopoGridExt.m_jb;
      int jmax = mTopoGridExt.m_je;
// Number of extra ghost points = m_ext_ghost_points
      int egh = m_iStart[gTop]-imin;

// Update extra ghost points. Do not worry about processor boundaries,
// they will be overwritten with correct values in the communication update afterward.
#pragma omp parallel for
      for( int i=imin+egh ; i <= imax-egh ; i++ )
	 for( int q = 0 ; q < egh ; q++ )
	 {
	    mTopoGridExt(i,jmin+q,1)   = mTopoGridExt(i,jmin+egh,1);
	    mTopoGridExt(i,jmax-q,1)   = mTopoGridExt(i,jmax-egh,1);
	 }
#pragma omp parallel for
      for( int j=jmin ; j <= jmax ; j++ )
	 for( int q = 0 ; q < egh ; q++ )
	 {
	    mTopoGridExt(imin+q,j,1) = mTopoGridExt(imin+egh,j,1);
	    mTopoGridExt(imax-q,j,1) = mTopoGridExt(imax-egh,j,1);
	 }
      communicate_array_2d_ext( mTopoGridExt );
   }
}

//-----------------------------------------------------------------------
void EW::gettopowgh( float_sw4 ai, float_sw4 wgh[8] ) const
{
   float_sw4 pol = ai*ai*ai*ai*ai*ai*ai*(-251+135*ai+25*ai*ai-
                                      33*ai*ai*ai+6*ai*ai*ai*ai)/720;
   wgh[0] = -1.0/60*ai + 1.0/180*ai*ai + 1.0/48*ai*ai*ai + 23.0/144*ai*ai*ai*ai 
      - (17.0*ai + 223.0)*ai*ai*ai*ai*ai/720 - pol;
   wgh[1] = 3.0/20*ai -3.0/40*ai*ai -1.0/6*ai*ai*ai - 13.0/12*ai*ai*ai*ai + 
      97.0/45*ai*ai*ai*ai*ai + 1.0/6*ai*ai*ai*ai*ai*ai + 7*pol;
   wgh[2] = -0.75*ai +0.75*ai*ai+(13.0+155*ai)*ai*ai*ai/48 -103.0/16*ai*ai*ai*ai*ai
      - 121.0/240*ai*ai*ai*ai*ai*ai - 21*pol;
   wgh[3] = 1 - 49.0/36*ai*ai - 49.0/9*ai*ai*ai*ai+385.0/36*ai*ai*ai*ai*ai +
      61.0/72*ai*ai*ai*ai*ai*ai + 35*pol;
   wgh[4] = 0.75*ai + 0.75*ai*ai - 13.0/48*ai*ai*ai + 89.0/16*ai*ai*ai*ai - 
         1537.0/144*ai*ai*ai*ai*ai - 41.0/48*ai*ai*ai*ai*ai*ai - 35*pol;
   wgh[5] = -3.0/20*ai - 3.0/40*ai*ai + 1.0/6*ai*ai*ai - 41.0/12*ai*ai*ai*ai
      + 6.4*ai*ai*ai*ai*ai + 31.0/60*ai*ai*ai*ai*ai*ai + 21*pol;
   wgh[6] = 1.0/60*ai + 1.0/180*ai*ai - 1.0/48*ai*ai*ai + 167.0/144*ai*ai*ai*ai -
      1537.0/720*ai*ai*ai*ai*ai- 25.0/144*ai*ai*ai*ai*ai*ai - 7*pol;
   wgh[7] = -1.0/6*ai*ai*ai*ai + 11.0/36*ai*ai*ai*ai*ai + 1.0/40*ai*ai*ai*ai*ai*ai + pol;
}

//-----------------------------------------------------------------------
void EW::smooth_grid( int maxIter )
{
// Smooth the grid (only the Z component for now)
// NOTE: the current smoothing algorithm makes the error larger rather than smaller!
   float_sw4 rf=0.05; // rf<1/6 for stability
   if (mVerbose >= 1 && proc_zero() && maxIter>0)
      cout << "***smoothing the grid with " << maxIter << " Jacobi iterations and relaxation factor " << rf << " ***"<< endl;

   int topLevel = mNumberOfGrids-1;
   int g = mNumberOfGrids-1;
   int iter;
// temporary storage: How can I use mJ for temporary storage?
   Sarray tmp;
   tmp.define(m_iStart[topLevel],m_iEnd[topLevel],m_jStart[topLevel],m_jEnd[topLevel],m_kStart[topLevel],m_kEnd[topLevel]);

// initialize to make the Dirichlet boundary conditions work
#pragma omp parallel for
   for (int k = m_kStart[topLevel]; k <= m_kEnd[topLevel]; k++)
      for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
	 for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
	 {
	    tmp(i,j,k) = mZ[g](i,j,k);
	 }

// Laplacian filter
   for (iter=0; iter < maxIter; iter++)
   {
// loop over all interior points
#pragma omp parallel for
      for (int k = m_kStart[topLevel]+m_ghost_points+1; k <= m_kEnd[topLevel]-m_ghost_points-2; k++)
	 for (int j = m_jStart[topLevel]+1; j <= m_jEnd[topLevel]-1; j++)
	    for (int i = m_iStart[topLevel]+1; i <= m_iEnd[topLevel]-1; i++)
	    {
	       tmp(i,j,k) = mZ[g](i,j,k) + rf*(mZ[g](i+1,j,k) + mZ[g](i-1,j,k) + mZ[g](i,j+1,k) + mZ[g](i,j-1,k) + mZ[g](i,j,k+1) + mZ[g](i,j,k-1) - 6.*mZ[g](i,j,k));
	    }

// impose Neumann bc on the i and j sides
#pragma omp parallel for
      for (int k = m_kStart[topLevel]+m_ghost_points+1; k <= m_kEnd[topLevel]-m_ghost_points-2; k++)
      {
	 for (int j = m_jStart[topLevel]+1; j <= m_jEnd[topLevel]-1; ++j)
	 {
	    int i = m_iStart[topLevel];
	    tmp(i,j,k) = tmp(i+1,j,k);
	    i = m_iEnd[topLevel];
	    tmp(i,j,k) = tmp(i-1,j,k);
	 }

	 for (int i = m_iStart[topLevel]+1; i <= m_iEnd[topLevel]-1; ++i)
	 {
	    int j = m_jStart[topLevel];
	    tmp(i,j,k) = tmp(i,j+1,k);
	    j = m_jEnd[topLevel];
	    tmp(i,j,k) = tmp(i,j-1,k);
	 }
// Corners
	 int i = m_iStart[topLevel];
	 int j = m_jStart[topLevel];
	 tmp(i,j,k) = tmp(i+1,j+1,k);

	 i = m_iEnd[topLevel];
	 j = m_jStart[topLevel];
	 tmp(i,j,k) = tmp(i-1,j+1,k);

	 i = m_iStart[topLevel];
	 j = m_jEnd[topLevel];
	 tmp(i,j,k) = tmp(i+1,j-1,k);
    
	 i = m_iEnd[topLevel];
	 j = m_jEnd[topLevel];
	 tmp(i,j,k) = tmp(i-1,j-1,k);
      } // end Neumann loop

// communicate parallel ghost points
      communicate_array( tmp, topLevel );

// update solution (Dirichlet are imposed implicitly by never changing the tmp array along the top or bottom boundary)
#pragma omp parallel for
      for (int k = m_kStart[topLevel]; k <= m_kEnd[topLevel]; k++)
	 for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
	    for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
	       mZ[g](i,j,k) = tmp(i,j,k);

   }// end for iter (grid smoother)
}
