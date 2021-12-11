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
#include "EW.h"
#include "F77_FUNC.h"
#include "caliper.h"
#include "mpi.h"
extern "C" {
void metric(int*, int*, int*, int*, int*, int*, double*, double*, double*,
            double*, double*, int*);
void gridinfo(int*, int*, int*, int*, int*, int*, double*, double*, double*,
              double*);
void metricexgh(int*, int*, int*, int*, int*, int*, int*, int*, int*, double*,
                double*, double*, double*, double*, int*, double*, double*,
                double*, double*, double*, double*, double*);
void meterr4c(int*, int*, int*, int*, int*, int*, double*, double*, double*,
              double*, double*, double*, int*, int*, int*, int*, int*, int*,
              double*);
}
#define SQR(x) ((x) * (x))
//---------------------------------------------------------
// void EW::setup_metric() {
//   if (!m_topography_exists) return;
//   if (mVerbose >= 1 && proc_zero()) cout << "***inside setup_metric***" <<
//   endl; int g = mNumberOfGrids - 1; int Bx = m_iStart[g]; int By =
//   m_jStart[g]; int Bz = m_kStart[g]; int Nx = m_iEnd[g]; int Ny = m_jEnd[g];
//   int Nz = m_kEnd[g];

//   if (m_analytical_topo && m_use_analytical_metric) {
//     // Gaussina hill topography, analytical expressions for metric
//     derivatives. int nxg = m_global_nx[g]; int nyg = m_global_ny[g]; int nzg
//     = m_global_nz[g]; float_sw4 h = mGridSize[g]; float_sw4 zmax = m_zmin[g -
//     1] - (nzg - 1) * h * (1 - m_zetaBreak); if (m_croutines)
//       metricexgh_ci(Bx, Nx, By, Ny, Bz, Nz, nzg, mX.c_ptr(), mY.c_ptr(),
//                     mZ.c_ptr(), mMetric.c_ptr(), mJ.c_ptr(),
//                     m_grid_interpolation_order, m_zetaBreak, zmax,
//                     m_GaussianAmp, m_GaussianXc, m_GaussianYc, m_GaussianLx,
//                     m_GaussianLy);

//     else
//       metricexgh(&Bx, &Nx, &By, &Ny, &Bz, &Nz, &nxg, &nyg, &nzg, mX.c_ptr(),
//                  mY.c_ptr(), mZ.c_ptr(), mMetric.c_ptr(), mJ.c_ptr(),
//                  &m_grid_interpolation_order, &m_zetaBreak, &zmax,
//                  &m_GaussianAmp, &m_GaussianXc, &m_GaussianYc, &m_GaussianLx,
//                  &m_GaussianLy);
//   } else {
//     int ierr = 0;
//     if (m_croutines)
//       ierr = metric_ci(Bx, Nx, By, Ny, Bz, Nz, mX.c_ptr(), mY.c_ptr(),
//                        mZ.c_ptr(), mMetric.c_ptr(), mJ.c_ptr());
//     else
//       metric(&Bx, &Nx, &By, &Ny, &Bz, &Nz, mX.c_ptr(), mY.c_ptr(),
//       mZ.c_ptr(),
//              mMetric.c_ptr(), mJ.c_ptr(), &ierr);
//     CHECK_INPUT(ierr == 0, "Problems calculating the metric coefficients");
//   }

//   communicate_array(mMetric, mNumberOfGrids - 1);
//   communicate_array(mJ, mNumberOfGrids - 1);

//   if (m_analytical_topo && !m_use_analytical_metric && mVerbose > 3)
//     // Test metric derivatives if available
//     metric_derivatives_test();

//   float_sw4 minJ, maxJ;
//   if (m_croutines)
//     gridinfo_ci(Bx, Nx, By, Ny, Bz, Nz, mMetric.c_ptr(), mJ.c_ptr(), minJ,
//                 maxJ);
//   else
//     gridinfo(&Bx, &Nx, &By, &Ny, &Bz, &Nz, mMetric.c_ptr(), mJ.c_ptr(),
//     &minJ,
//              &maxJ);
//   float_sw4 minJglobal, maxJglobal;
//   MPI_Allreduce(&minJ, &minJglobal, 1, m_mpifloat, MPI_MIN,
//                 m_cartesian_communicator);
//   MPI_Allreduce(&maxJ, &maxJglobal, 1, m_mpifloat, MPI_MAX,
//                 m_cartesian_communicator);
//   if (mVerbose > 3 && proc_zero())
//     printf("*** Jacobian of metric: minJ = %e maxJ = %e\n", minJglobal,
//            maxJglobal);
//   // just save the results for now... do the sanity check later
//   m_minJacobian = minJglobal;
//   m_maxJacobian = maxJglobal;
// }

//-----------------------------------------------------------------------
// void EW::generate_grid() {
//   SW4_MARK_FUNCTION;
//   // Generate grid on domain: topography <= z <= zmax,
//   // The 2D grid on z=zmax, is given by ifirst <= i <= ilast, jfirst <= j <=
//   // jlast spacing h.
//   if (!m_topography_exists) return;

//   //  m_grid_interpolation_order = a_order;

//   if (mVerbose >= 1 && proc_zero())
//     cout << "***inside generate_grid***" << endl;

//   // get the size from the top Cartesian grid
//   int g = mNumberOfCartesianGrids - 1;
//   int ifirst = m_iStart[g];
//   int ilast = m_iEnd[g];
//   int jfirst = m_jStart[g];
//   int jlast = m_jEnd[g];

//   float_sw4 h = mGridSize[g];  // grid size must agree with top cartesian
//   grid float_sw4 zMaxCart = m_zmin[g];  // bottom z-level for curvilinear
//   grid

//   //  int i, j;
//   int gTop = mNumberOfGrids - 1;
//   int Nz = m_kEnd[gTop] - m_ghost_points;

//   if (mVerbose > 4 && proc_zero()) {
//     printf(
//         "generate_grid: Number of grid points in curvilinear grid = %i,
//         kStart "
//         "= %i, kEnd = %i\n",
//         Nz, m_kStart[gTop], m_kEnd[gTop]);
//   }

// // generate the grid by calling the curvilinear mapping function
// #pragma omp parallel for
//   for (int k = m_kStart[gTop]; k <= m_kEnd[gTop]; k++)
//     for (int j = m_jStart[gTop]; j <= m_jEnd[gTop]; j++)
//       for (int i = m_iStart[gTop]; i <= m_iEnd[gTop]; i++) {
//         float_sw4 X0, Y0, Z0;
//         curvilinear_grid_mapping((float_sw4)i, (float_sw4)j, (float_sw4)k,
//         X0,
//                                  Y0, Z0);
//         mX(i, j, k) = X0;
//         mY(i, j, k) = Y0;
//         mZ(i, j, k) = Z0;
//       }
//   // tmp
//   // test the inverse mapping
//   //  double q0, r0, s0, dist=0.;
//   //  for (k=m_kStart[gTop]; k<=m_kEnd[gTop]; k++)
//   //    for (j=m_jStart[gTop]; j<=m_jEnd[gTop]; j++)
//   //      for (i=m_iStart[gTop]; i<=m_iEnd[gTop]; i++)
//   //      {
//   // 	invert_curvilinear_grid_mapping(mX(i,j,k), mY(i,j,k), mZ(i,j,k),
//   q0, r0,
//   // s0); 	dist += SQR(q0 - (double) i) + SQR(r0 - (double) j) + SQR(s0 -
//   // (double) k);
//   //      }
//   //   double totalDist;
//   //   MPI_Allreduce( &dist, &totalDist, 1, MPI_DOUBLE, MPI_SUM,
//   //   m_cartesian_communicator ); if (m_myRank == 0)
//   //     printf("L2-error in inverse of curvilinear mapping = %e\n",
//   //     sqrt(totalDist));
//   // end test

//   // make sure all processors have made their grid before we continue
//   communicate_array(mZ, gTop);

//   // tmp
//   // calculate min and max((mZ(i,j,k)-mZ(i,j,k-1))/h) for k=Nz
//   int k = Nz;
//   // tmp
//   float_sw4 mZmin = 1.0e9, mZmax = 0;
// #pragma omp parallel for reduction(min : mZmin) reduction(max : mZmax)
//   for (int j = m_jStart[gTop]; j <= m_jEnd[gTop]; j++)
//     for (int i = m_iStart[gTop]; i <= m_iEnd[gTop]; i++) {
//       float_sw4 hRatio = (mZ(i, j, k) - mZ(i, j, k - 1)) / mGridSize[gTop];
//       if (hRatio < mZmin) mZmin = hRatio;
//       if (hRatio > mZmax) mZmax = hRatio;
//     }
//   float_sw4 zMinGlobal, zMaxGlobal;
//   MPI_Allreduce(&mZmin, &zMinGlobal, 1, m_mpifloat, MPI_MIN,
//                 m_cartesian_communicator);
//   MPI_Allreduce(&mZmax, &zMaxGlobal, 1, m_mpifloat, MPI_MAX,
//                 m_cartesian_communicator);
//   if (mVerbose > 3 && proc_zero()) {
//     printf(
//         "Curvilinear/Cartesian interface (k=Nz-1): Min grid size ratio - 1 =
//         "
//         "%e, max ratio z - 1 = %e, top grid # = %i\n",
//         zMinGlobal - 1., zMaxGlobal - 1., gTop);
//   }
//   // end tmp
// }

//-----------------------------------------------------------------------
bool EW::curvilinear_grid_mapping(float_sw4 q, float_sw4 r, float_sw4 s,
                                  float_sw4& X0, float_sw4& Y0, float_sw4& Z0) {
  // if (q,r) is on this processor (need a 2x2 interval in (i,j)-index space:
  // Return true and assign (X0,Y0,Z0) corresponding to (q,r,s)

  // Returns false if
  // 1) (q,r,s) is outside the global parameter domain (expanded by ghost
  // points) 2) There is no curvilinear grid. Still computes (X0, Y0) based on
  // the top Cartesian grid size 3) (q,r) is not on this processor. Still
  // computes (X0, Y0)

  // NOTE:
  // The parameters are normalized such that 1 <= q <= Nx is the full domain
  // (without ghost points),
  //  1 <= r <= Ny, 1 <= s <= Nz.

  int gCurv = mNumberOfGrids - 1;
  float_sw4 h = mGridSize[gCurv];
  // check global parameter space
  float_sw4 qMin = (float_sw4)(1 - m_ghost_points);
  float_sw4 qMax = (float_sw4)(m_global_nx[gCurv] + m_ghost_points);
  float_sw4 rMin = (float_sw4)(1 - m_ghost_points);
  float_sw4 rMax = (float_sw4)(m_global_ny[gCurv] + m_ghost_points);
  float_sw4 sMin = (float_sw4)m_kStart[gCurv];
  float_sw4 sMax = (float_sw4)m_kEnd[gCurv];

  if (!(q >= qMin && q <= qMax && r >= rMin && r <= rMax && s >= sMin &&
        s <= sMax)) {
    cout
        << "curvilinear_grid_mapping: input parameters out of bounds (q,r,s) = "
        << q << ", " << r << ", " << s << endl;
    cout << "limits are " << qMin << " " << qMax << " " << rMin << " " << rMax
         << " " << sMin << " " << sMax << endl;
    return false;
  }

  X0 = (q - 1.0) * h;
  Y0 = (r - 1.0) * h;

  if (!topographyExists()) return false;

  // bottom z-level for curvilinear grid = top z-level for highest Cartesian
  // grid
  float_sw4 zMaxCart = m_zmin[mNumberOfCartesianGrids - 1];

  // ************************
  // compute index interval based on (q,r)
  int iNear, jNear, kNear, g, i, j, k;

  Z0 = zMaxCart - h;  // to make computeNearestGridPoint think we are in the
                      // curvilinear grid
  computeNearestGridPoint(iNear, jNear, kNear, g, X0, Y0, Z0);

  if (g != gCurv) return false;

  float_sw4 tau;  // holds the elevation at (q,r). Recall that elevation=-z
  if (m_analytical_topo) {
    tau = m_GaussianAmp * exp(-SQR((X0 - m_GaussianXc) / m_GaussianLx) -
                              SQR((Y0 - m_GaussianYc) / m_GaussianLy));
  } else  // general case: interpolate mTopoGrid array
  {
    // if (X0, Y0) falls within roundoff of grid point (iNear,jNear), we only
    // need that grid point on this proc, otherwise we need the 2x2 area [i,i+1]
    // by [j,j+1]

    float_sw4 xPt = (iNear - 1) * h;
    float_sw4 yPt = (jNear - 1) * h;

    // first check if we are very close to a grid point
    bool smackOnTop =
        (fabs((xPt - X0) / h) < 1.e-9 && fabs((yPt - Y0) / h) < 1.e-9);

    if (smackOnTop) {
      if (!point_in_proc_ext(iNear, jNear, gCurv)) return false;
      tau = mTopoGridExt(iNear, jNear, 1);
    } else {
      computeNearestLowGridPoint(i, j, k, g, X0, Y0, Z0);
      // There are some subtle issues with the bi-cubic interpolation near
      // parallel processor boundaries, see invert_curvilinear_mapping (below)
      // for a discussion

      // bi-cubic interpolation for O(h^4) accuracy
      //       if ((point_in_proc(i-1,j-1,gCurv) && point_in_proc(i,j-1,gCurv)
      //       && point_in_proc(i+1,j-1,gCurv) && point_in_proc(i+2,j-1,gCurv)
      //       &&
      // 		point_in_proc(i-1,j,gCurv) && point_in_proc(i,j,gCurv)
      // && point_in_proc(i+1,j,gCurv) && point_in_proc(i+2,j,gCurv) &&
      // 		point_in_proc(i-1,j+1,gCurv) &&
      // point_in_proc(i,j+1,gCurv) && point_in_proc(i+1,j+1,gCurv) &&
      // point_in_proc(i+2,j+1,gCurv) &&
      // point_in_proc(i-1,j+2,gCurv) && point_in_proc(i,j+2,gCurv) &&
      // point_in_proc(i+1,j+2,gCurv) && point_in_proc(i+2,j+2,gCurv) ) )
      //       {
      // 	double Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj,
      // tjp1, tjp2; 	Qim1 = (q-i)*(q-i-1)*(q-i-2)/(-6.); 	Qi   =
      // (q-i+1)*(q-i-1)*(q-i-2)/(2.); 	Qip1 = (q-i+1)*(q-i)*(q-i-2)/(-2.);
      // Qip2 = (q-i+1)*(q-i)*(q-i-1)/(6.);

      // 	Rjm1 = (r-j)*(r-j-1)*(r-j-2)/(-6.);
      // 	Rj   = (r-j+1)*(r-j-1)*(r-j-2)/(2.);
      // 	Rjp1 = (r-j+1)*(r-j)*(r-j-2)/(-2.);
      // 	Rjp2 = (r-j+1)*(r-j)*(r-j-1)/(6.);

      // 	tjm1 = Qim1*mTopoGrid(i-1,j-1,1) + Qi*mTopoGrid(i,j-1,1) +
      // Qip1*mTopoGrid(i+1,j-1,1) +  Qip2*mTopoGrid(i+2,j-1,1); 	tj   =
      // Qim1*mTopoGrid(i-1,j,1) + Qi*mTopoGrid(i,j,1) + Qip1*mTopoGrid(i+1,j,1)
      // +  Qip2*mTopoGrid(i+2,j,1); 	tjp1 = Qim1*mTopoGrid(i-1,j+1,1) +
      // Qi*mTopoGrid(i,j+1,1) +  Qip1*mTopoGrid(i+1,j+1,1) +
      // Qip2*mTopoGrid(i+2,j+1,1); 	tjp2 = Qim1*mTopoGrid(i-1,j+2,1) +
      // Qi*mTopoGrid(i,j+2,1) +  Qip1*mTopoGrid(i+1,j+2,1) +
      // Qip2*mTopoGrid(i+2,j+2,1);

      // 	tau = Rjm1*tjm1 + Rj*tj + Rjp1*tjp1 + Rjp2*tjp2;
      //       }
      //      else if ( ( point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv)
      //      && point_in_proc(i,j+1,gCurv) &&
      //	     point_in_proc(i+1,j+1,gCurv) ) )
      //      if ( ( point_in_proc_ext(i,j,gCurv)   &&
      //      point_in_proc_ext(i+1,j,gCurv) &&
      //	     point_in_proc_ext(i,j+1,gCurv) &&
      // point_in_proc_ext(i+1,j+1,gCurv) ) )
      //      {
      //// linear interpolation to define topography between grid points
      //	double Qi, Qip1, Rj, Rjp1;
      //	Qi = (i+1 - q);
      //	Qip1 = (q - i);
      //	Rj = (j+1 - r);
      //	Rjp1 = (r - j);
      //	tau = mTopoGridExt(i,j,1)*Rj*Qi + mTopoGridExt(i,j+1,1)*Rjp1*Qi
      //+ mTopoGridExt(i+1,j,1)*Rj*Qip1 +
      // mTopoGridExt(i+1,j+1,1)*Rjp1*Qip1;
      //      }
      if (point_in_proc_ext(i - 3, j - 3, gCurv) &&
          point_in_proc_ext(i + 4, j + 4, gCurv)) {
        float_sw4 a6cofi[8], a6cofj[8];
        gettopowgh(q - i, a6cofi);
        gettopowgh(r - j, a6cofj);
        tau = 0;
        for (int l = j - 3; l <= j + 4; l++)
          for (int k = i - 3; k <= i + 4; k++)
            tau +=
                a6cofi[k - i + 3] * a6cofj[l - j + 3] * mTopoGridExt(k, l, 1);
      } else {
        return false;
      }
    }  // end else...(not smackOnTop)

  }  // end general case: interpolating mTopoGrid array

  // now we need to calculate Z0 = Z(q,r,s)

  // setup parameters for grid mapping
  int Nz = m_kEnd[gCurv] - m_ghost_points;

  // zeta  > zetaBreak gives constant grid size = h
  float_sw4 sBreak = 1. + m_zetaBreak * (Nz - 1);

  float_sw4 zeta, c1 = 0., c2 = 0., c3 = 0., c4 = 0.0, c5 = 0.0, c6 = 0.0, zMax;
  zMax = zMaxCart - (Nz - sBreak) * h;

  // quadratic term to make variation in grid size small at bottom boundary
  c1 = zMax + tau - mGridSize[gCurv] * (sBreak - 1);
  // cubic term to make 2nd derivative zero at zeta=1
  if (m_grid_interpolation_order >= 3)
    c2 = c1;
  else
    c2 = 0.;

  // 4th order term takes care of 3rd derivative, but can make grid warp itself
  // inside out
  if (m_grid_interpolation_order >= 4)
    c3 = c2;
  else
    c3 = 0.;

  // Added continuity of the 4th and 5th derivatives
  if (m_grid_interpolation_order >= 5)
    c4 = c3;
  else
    c4 = 0;
  if (m_grid_interpolation_order >= 6)
    c5 = c4;
  else
    c5 = 0;
  if (m_grid_interpolation_order >= 7)
    c6 = c5;
  else
    c6 = 0;

  // the forward mapping is ...
  if (s <= (float_sw4)sBreak) {
    zeta = (s - 1) / (sBreak - 1.);
    Z0 = (1. - zeta) * (-tau) +
         zeta * (zMax + c1 * (1. - zeta) + c2 * SQR(1. - zeta) +
                 c3 * (1. - zeta) * SQR(1. - zeta) +
                 (c4 + c5 * (1 - zeta)) * SQR(1 - zeta) * SQR(1 - zeta) +
                 c6 * SQR(1 - zeta) * SQR(1 - zeta) * SQR(1 - zeta));
  } else {
    Z0 = zMax + (s - sBreak) * h;
  }

  return true;
}

//-----------------------------------------------------------------------
bool EW::invert_curvilinear_grid_mapping(float_sw4 X0, float_sw4 Y0,
                                         float_sw4 Z0, float_sw4& q,
                                         float_sw4& r, float_sw4& s) {
  // If (X0, Y0, Z0) is on the curvilinear grid and (X0, Y0) is on this
  // processor: Return true and assigns (q,r,s) corresponding to point
  // (mX0,mY0,mZ0)

  // Returns false if
  // 1) There is no curvilinear grid. Computes (q,r) based on the top Cartesian
  // grid size 2) Z0 > m_zmin[topCartGrid]. Still computes (q,r) as if Z0 would
  // be on the curvilinear grid 3) (X0, Y0) is not on this processor. Still
  // computes (q,r)

  // NOTE:
  // Normalize parameters such that 1 <= q <= Nx is the full domain (without
  // ghost points),
  //  1 <= r <= Ny, 1 <= s <= Nz.

  int gCurv = mNumberOfGrids - 1;
  float_sw4 h = mGridSize[gCurv];
  q = X0 / h + 1.0;
  r = Y0 / h + 1.0;
  s = 0.;

  if (!topographyExists()) {
    if (mVerbose >= 3)
      printf(
          "invert_curvilinear_mapping returning false because there is no "
          "curvilinear mapping!\n");

    return false;
  }

  // compute index interval based on (q,r)
  int iNear, jNear, kNear, g, i, j, k;

  // for ghost points in the curvilinear grid, calling computeNearestGridPoint
  // with Z0 will pick up the top Cartesian grid

  // bottom z-level for curvilinear grid = top z-level for highest Cartesian
  // grid
  float_sw4 zMaxCart = m_zmin[mNumberOfCartesianGrids - 1];

  computeNearestGridPoint(iNear, jNear, kNear, g, X0, Y0, zMaxCart - h);

  if (g != gCurv) {
    cout << "invert_curvilinear_grid_mapping: computeNearestGridPoint returned "
            "g!=gCurv for (X0, Y0, Z0) = "
         << X0 << ", " << Y0 << ", " << Z0 << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  float_sw4 tau;  // holds the elevation at (q,r). Recall that elevation=-z

  if (m_analytical_topo) {
    tau = m_GaussianAmp * exp(-SQR((X0 - m_GaussianXc) / m_GaussianLx) -
                              SQR((Y0 - m_GaussianYc) / m_GaussianLy));
  } else  // general case: interpolate mTopoGrid array
  {
    // if (X0, Y0) falls within roundoff of grid point (iNear,jNear), we only
    // need that grid point on this proc, otherwise we need the 2x2 area [i,i+1]
    // by [j,j+1]

    float_sw4 xPt = (iNear - 1) * h;
    float_sw4 yPt = (jNear - 1) * h;

    // first check if we are very close to a grid point
    bool smackOnTop =
        (fabs((xPt - X0) / h) < 1.e-9 && fabs((yPt - Y0) / h) < 1.e-9);

    if (smackOnTop) {
      if (!point_in_proc_ext(iNear, jNear, gCurv)) {
        if (mVerbose >= 3)
          printf(
              "invert_curvilinear_mapping returning false because iNear=%i, "
              "jNear=%i does not live on this proc!\n",
              iNear, jNear);

        return false;
      }
      tau = mTopoGridExt(iNear, jNear, 1);
    } else {
      computeNearestLowGridPoint(i, j, k, g, X0, Y0, Z0);

      // float_sw4 Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1,
      // tjp2;

      // The bi-cubic interpolation has the following problem:
      // When the source is close to a processor boundary, the same source must
      // be discretized on two (or more processors). However, the topography
      // grid only has 2 points overlap, which sometimes is insufficient for the
      // bi-cubic interpolation formula. In practice, the difference between the
      // bi-linear and bi-cubic formula makes very little difference, so it is
      // more robust to always use the bi-linear formula.

      // // bi-cubic interpolation for O(h^4) accuracy
      //       if ( (point_in_proc(i-1,j-1,gCurv) && point_in_proc(i,j-1,gCurv)
      //       && point_in_proc(i+1,j-1,gCurv) && point_in_proc(i+2,j-1,gCurv)
      //       &&
      // 	    point_in_proc(i-1,j,gCurv) && point_in_proc(i,j,gCurv) &&
      // point_in_proc(i+1,j,gCurv) && point_in_proc(i+2,j,gCurv) &&
      // 	    point_in_proc(i-1,j+1,gCurv) && point_in_proc(i,j+1,gCurv)
      // && point_in_proc(i+1,j+1,gCurv) && point_in_proc(i+2,j+1,gCurv) &&
      // 	    point_in_proc(i-1,j+2,gCurv) && point_in_proc(i,j+2,gCurv)
      // && point_in_proc(i+1,j+2,gCurv) && point_in_proc(i+2,j+2,gCurv) ) )
      //       {
      // 	Qim1 = (q-i)*(q-i-1)*(q-i-2)/(-6.);
      // 	Qi   = (q-i+1)*(q-i-1)*(q-i-2)/(2.);
      // 	Qip1 = (q-i+1)*(q-i)*(q-i-2)/(-2.);
      // 	Qip2 = (q-i+1)*(q-i)*(q-i-1)/(6.);

      // 	Rjm1 = (r-j)*(r-j-1)*(r-j-2)/(-6.);
      // 	Rj   = (r-j+1)*(r-j-1)*(r-j-2)/(2.);
      // 	Rjp1 = (r-j+1)*(r-j)*(r-j-2)/(-2.);
      // 	Rjp2 = (r-j+1)*(r-j)*(r-j-1)/(6.);

      // 	tjm1 = Qim1*mTopoGrid(i-1,j-1,1) + Qi*mTopoGrid(i,j-1,1) +
      // Qip1*mTopoGrid(i+1,j-1,1) +  Qip2*mTopoGrid(i+2,j-1,1); 	tj   =
      // Qim1*mTopoGrid(i-1,j,1) + Qi*mTopoGrid(i,j,1) + Qip1*mTopoGrid(i+1,j,1)
      // +  Qip2*mTopoGrid(i+2,j,1); 	tjp1 = Qim1*mTopoGrid(i-1,j+1,1) +
      // Qi*mTopoGrid(i,j+1,1) +  Qip1*mTopoGrid(i+1,j+1,1) +
      // Qip2*mTopoGrid(i+2,j+1,1); 	tjp2 = Qim1*mTopoGrid(i-1,j+2,1) +
      // Qi*mTopoGrid(i,j+2,1) +  Qip1*mTopoGrid(i+1,j+2,1) +
      // Qip2*mTopoGrid(i+2,j+2,1);

      // 	tau = Rjm1*tjm1 + Rj*tj + Rjp1*tjp1 + Rjp2*tjp2;
      // // tmp
      // 	printf("invert_curvilinear_mapping: q=%e, r=%e, Cubic tau=%e\n",
      // q, r, tau);
      //       }
      //      else if ( ( point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv)
      //      && point_in_proc(i,j+1,gCurv) &&
      //		  point_in_proc(i+1,j+1,gCurv) ) )
      //      if ( ( point_in_proc_ext(i,j,gCurv)   &&
      //      point_in_proc_ext(i+1,j,gCurv) &&
      //	     point_in_proc_ext(i,j+1,gCurv) &&
      // point_in_proc_ext(i+1,j+1,gCurv) ) )
      //      {
      // use linear interpolation if there are not enough points for the
      // bi-cubic formula
      //	Qi = (i+1 - q);
      //	Qip1 = (q - i);
      //	Rj = (j+1 - r);
      //	Rjp1 = (r - j);
      //	tau = mTopoGridExt(i,j,1)*Rj*Qi + mTopoGridExt(i,j+1,1)*Rjp1*Qi
      //+ mTopoGridExt(i+1,j,1)*Rj*Qip1 +
      // mTopoGridExt(i+1,j+1,1)*Rjp1*Qip1;
      // tmp
      //	printf("invert_curvilinear_mapping: q=%e, r=%e, Linear
      // tau=%e\n", q, r, tau);
      //      }
      if (point_in_proc_ext(i - 3, j - 3, gCurv) &&
          point_in_proc_ext(i + 4, j + 4, gCurv)) {
        float_sw4 a6cofi[8], a6cofj[8];
        gettopowgh(q - i, a6cofi);
        gettopowgh(r - j, a6cofj);
        tau = 0;
        for (int l = j - 3; l <= j + 4; l++)
          for (int k = i - 3; k <= i + 4; k++)
            tau +=
                a6cofi[k - i + 3] * a6cofj[l - j + 3] * mTopoGridExt(k, l, 1);
      } else {
        if (mVerbose >= 3)
          printf(
              "invert_curvilinear_mapping returning false from bi-lin "
              "interpolation\n"
              "not all points around i=%i and j=%i are on this proc!\n",
              i, j);
        return false;
      }

    }  // end else... (not smackOnTop)

  }  // end else... general case

  // now we need to calculate s: Z(q,r,s) = Z0

  // setup parameters for grid mapping (same as curvilinear_grid_mapping)
  int Nz = m_kEnd[gCurv] - m_ghost_points;
  // zeta  > zetaBreak gives constant grid size = h
  float_sw4 sBreak = 1. + m_zetaBreak * (Nz - 1);

  float_sw4 zeta, c1 = 0., c2 = 0., c3 = 0., c4 = 0, c5 = 0, c6 = 0, zMax;
  zMax = zMaxCart - (Nz - sBreak) * h;

  // quadratic term to make variation in grid size small at bottom boundary
  c1 = zMax + tau - mGridSize[gCurv] * (sBreak - 1);
  // cubic term to make 2nd derivative zero at zeta=1
  if (m_grid_interpolation_order >= 3)
    c2 = c1;
  else
    c2 = 0.;

  // 4th order term takes care of 3rd derivative, but can make grid warp itself
  // inside out
  if (m_grid_interpolation_order >= 4)
    c3 = c2;
  else
    c3 = 0.;

  // Added continuity of the 4th and 5th derivatives
  if (m_grid_interpolation_order >= 5)
    c4 = c3;
  else
    c4 = 0;
  if (m_grid_interpolation_order >= 6)
    c5 = c4;
  else
    c5 = 0;
  if (m_grid_interpolation_order >= 7)
    c6 = c5;
  else
    c6 = 0;

  // the forward mapping is ...
  //  if (s <= (double) sBreak)
  //  {
  //    zeta = (s-1)/(sBreak-1.);
  //    Z0 = (1.-zeta)*(-tau)
  //      + zeta*(zMax + c1*(1.-zeta) + c2*SQR(1.-zeta) +
  //      c3*(1.-zeta)*SQR(1.-zeta));
  //  }
  //  else
  //  {
  //    Z0 = zMax + (s-sBreak)*h;
  //  }

  if (Z0 >= zMax) {
    s = (Z0 - zMax) / h + sBreak;
  } else {
    // Get initial guess by taking c1=c2=c3=0
    zeta = (Z0 + tau) / (zMax + tau);
    float_sw4 F0, Fp0, zeta0 = zeta;

    // Do a couple of Newton iterations
    int maxIter = 5;
    //      cout << "invert_curvilinear_grid_mapping: initial guess zeta0= " <<
    //      zeta0 << endl;
    for (int iter = 0; iter < maxIter; iter++) {
      F0 = (1. - zeta0) * (-tau) +
           zeta0 * (zMax + c1 * (1. - zeta0) + c2 * SQR(1. - zeta0) +
                    c3 * (1. - zeta0) * SQR(1. - zeta0) +
                    (c4 + c5 * (1 - zeta0) + c6 * SQR(1 - zeta0)) *
                        SQR(1 - zeta0) * SQR(1 - zeta0)) -
           Z0;
      Fp0 = tau + zMax + c1 * (1. - zeta0) + c2 * SQR(1. - zeta0) +
            c3 * (1. - zeta0) * SQR(1. - zeta0) +
            (c4 + (1 - zeta0) * c5 + c6 * SQR(1 - zeta0)) * SQR(1 - zeta0) *
                SQR(1 - zeta0) +
            zeta0 * (-c1 - 2. * c2 * (1. - zeta0) - 3. * c3 * SQR(1. - zeta0) -
                     4 * c4 * (1 - zeta0) * SQR(1 - zeta0) -
                     5 * c5 * SQR(1 - zeta0) * SQR(1 - zeta0) -
                     6 * c6 * (1 - zeta0) * SQR(1 - zeta0) * SQR(1 - zeta0));
      // tmp
      //        cout << "invert_curvilinear_grid_mapping: iter= " << iter << "
      //        zeta= " << zeta0 <<  " F0= " << F0 << " Fp0= " << Fp0 << endl;
      zeta0 -= F0 / Fp0;
    }
    // check convergence
    float_sw4 tol = 1.e-9;
    if (sizeof(float_sw4) == 4) tol = 1.e-5;
    if (fabs(F0) > tol) {
      cout << "invert_curvilinear_grid_mapping: poor convergence for X0, Y0, "
              "Z0 = "
           << X0 << ", " << Y0 << ", " << Z0 << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    zeta = zeta0;
    // tmp
    //    cout << "invert_curvilinear_grid_mapping: final zeta= " << zeta <<
    //    endl;
    s = 1. + (sBreak - 1.) * zeta;
  }  // end Z0 < zMax

  // tmp
  //    printf("invert_curvilinear_mapping: X0=%e, Y0=%e, Z0=%e, q=%e, r=%e,
  //    s=%e\n", X0, Y0, Z0, q, r, s);

  // if we got this far, the inversion was successful
  return true;
}

//-----------------------------------------------------------------------
void EW::smoothTopography(int maxIter) {
  //   int myRank;
  //   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // tmp
  //   if (proc_zero())
  //   {
  //     cout << "***inside smoothTopography***  maxIter= " << maxIter << endl;
  //   }

  float_sw4 rf = 0.2;  // rf<0.25 for stability
  //  int topLevel = mNumberOfGrids-1;
  int iter;

  copy_topo_to_topogridext();

  int imin = mTopoGridExt.m_ib;
  int imax = mTopoGridExt.m_ie;
  int jmin = mTopoGridExt.m_jb;
  int jmax = mTopoGridExt.m_je;

  // temporary storage
  Sarray tmp;
  tmp.define(imin, imax, jmin, jmax, 1, 1);

  // Laplacian filter
  for (iter = 0; iter < maxIter; iter++) {
#pragma omp parallel for
    for (int i = imin + 1; i <= imax - 1; ++i)
      for (int j = jmin + 1; j <= jmax - 1; ++j) {
        tmp(i, j, 1) =
            mTopoGridExt(i, j, 1) +
            rf * (mTopoGridExt(i + 1, j, 1) + mTopoGridExt(i - 1, j, 1) +
                  mTopoGridExt(i, j + 1, 1) + mTopoGridExt(i, j - 1, 1) -
                  4. * mTopoGridExt(i, j, 1));
      }

// Neumann boundary conditions
#pragma omp parallel for
    for (int j = jmin + 1; j <= jmax - 1; ++j) {
      int i = imin;
      tmp(i, j, 1) = tmp(i + 1, j, 1);
      i = imax;
      tmp(i, j, 1) = tmp(i - 1, j, 1);
    }

#pragma omp parallel for
    for (int i = imin + 1; i <= imax - 1; ++i) {
      int j = jmin;
      tmp(i, j, 1) = tmp(i, j + 1, 1);
      j = jmax;
      tmp(i, j, 1) = tmp(i, j - 1, 1);
    }
    // Corners
    int i = imin;
    int j = jmin;
    tmp(i, j, 1) = tmp(i + 1, j + 1, 1);

    i = imax;
    j = jmin;
    tmp(i, j, 1) = tmp(i - 1, j + 1, 1);

    i = imin;
    j = jmax;
    tmp(i, j, 1) = tmp(i + 1, j - 1, 1);

    i = imax;
    j = jmax;
    tmp(i, j, 1) = tmp(i - 1, j - 1, 1);

    communicate_array_2d_ext(tmp);

// update solution
#pragma omp parallel for
    for (int i = imin; i <= imax; ++i)
      for (int j = jmin; j <= jmax; ++j) mTopoGridExt(i, j, 1) = tmp(i, j, 1);
  }  // end for iter
}

//-----------------------------------------------------------------------
void EW::buildGaussianHillTopography(float_sw4 amp, float_sw4 Lx, float_sw4 Ly,
                                     float_sw4 x0, float_sw4 y0) {
  if (mVerbose >= 1 && proc_zero())
    cout << "***inside buildGaussianHillTopography***" << endl;

  int topLevel = mNumberOfGrids - 1;

  // copy data
  m_analytical_topo = true;
  //  m_analytical_topo = false;
  m_GaussianAmp = amp;
  m_GaussianLx = Lx;
  m_GaussianLy = Ly;
  m_GaussianXc = x0;
  m_GaussianYc = y0;

#pragma omp parallel for
  for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j) {
      float_sw4 x = (i - 1) * mGridSize[topLevel];
      float_sw4 y = (j - 1) * mGridSize[topLevel];
      // positive topography  is up (negative z)
      mTopo(i, j, 1) =
          m_GaussianAmp * exp(-SQR((x - m_GaussianXc) / m_GaussianLx) -
                              SQR((y - m_GaussianYc) / m_GaussianLy));
    }

#pragma omp parallel for
  for (int i = mTopoGridExt.m_ib; i <= mTopoGridExt.m_ie; ++i)
    for (int j = mTopoGridExt.m_jb; j <= mTopoGridExt.m_je; ++j) {
      float_sw4 x = (i - 1) * mGridSize[topLevel];
      float_sw4 y = (j - 1) * mGridSize[topLevel];
      // positive topography  is up (negative z)
      mTopoGridExt(i, j, 1) =
          m_GaussianAmp * exp(-SQR((x - m_GaussianXc) / m_GaussianLx) -
                              SQR((y - m_GaussianYc) / m_GaussianLy));
    }
}

//-----------------------------------------------------------------------
int EW::interpolate_topography(float_sw4 q, float_sw4 r, float_sw4& Z0,
                               bool smoothed) {
  // Interpolate the smoothed or raw topography
  //
  // if (q,r) is on this processor (need a 2x2 interval in (i,j)-index space:
  // Return true and assign Z0 corresponding to (q,r)

  // Returns false if
  // 1) (q,r) is outside the global parameter domain (expanded by ghost points)
  // 2) There is no curvilinear grid.
  // 3) (q,r) is not on this processor

  // NOTE:
  // The parameters are normalized such that 1 <= q <= Nx is the full domain
  // (without ghost points),
  //  1 <= r <= Ny.

  int gCurv = mNumberOfGrids - 1;
  float_sw4 h = mGridSize[gCurv];
  // check global parameter space

  float_sw4 qMin = (float_sw4)(1 - m_ghost_points);
  float_sw4 qMax = (float_sw4)(m_global_nx[gCurv] + m_ghost_points);
  float_sw4 rMin = (float_sw4)(1 - m_ghost_points);
  float_sw4 rMax = (float_sw4)(m_global_ny[gCurv] + m_ghost_points);

  // with mesh refinement, ghost points on coarser levels are further away
  //  int padding = m_ghost_points*((int) pow(2.0, mNumberOfCartesianGrids-1));

  float_sw4 sMin = (float_sw4)m_kStart[gCurv];
  float_sw4 sMax = (float_sw4)m_kEnd[gCurv];

  //  bool inside_extended_domain = (q >= qMin0 && q <= qMax0 && r >= rMin0 && r
  //  <= rMax0);
  bool inside_domain = (q >= qMin && q <= qMax && r >= rMin && r <= rMax);

  if (!inside_domain) {
    cout << "interpolate_topography: input parameters out of bounds (q,r) = "
         << q << ", " << r << endl;
    return -11;
  }

  float_sw4 X0 = (q - 1.0) * h;
  float_sw4 Y0 = (r - 1.0) * h;

  if (!topographyExists()) return -12;

  float_sw4 zMaxCart = m_zmin[mNumberOfCartesianGrids - 1];

  // ************************
  // compute index interval based on (q,r)
  int iNear, jNear, kNear, g, i, j, k;

  Z0 = zMaxCart - h;  // to make computeNearestGridPoint think we are in the
                      // curvilinear grid
  computeNearestGridPoint(iNear, jNear, kNear, g, X0, Y0, Z0);

  if (g != gCurv) {
    cout << "interpolate_topography: g = " << g << " gcurv = " << gCurv << endl;
    return -13;
  }
  float_sw4 tau;  // holds the elevation at (q,r). Recall that elevation=-z
  if (m_analytical_topo) {
    tau = m_GaussianAmp * exp(-SQR((X0 - m_GaussianXc) / m_GaussianLx) -
                              SQR((Y0 - m_GaussianYc) / m_GaussianLy));
  } else  // general case: interpolate mTopoGrid array
  {
    // if (X0, Y0) falls within roundoff of grid point (iNear,jNear), we only
    // need that grid point on this proc, otherwise we need the 2x2 area [i,i+1]
    // by [j,j+1]

    float_sw4 xPt = (iNear - 1) * h;
    float_sw4 yPt = (jNear - 1) * h;

    // first check if we are very close to a grid point
    bool smackOnTop =
        (fabs((xPt - X0) / h) < 1.e-9 && fabs((yPt - Y0) / h) < 1.e-9);

    // Note: Need to check point_in_proc, because with mesh refinement, (iNear,
    // jNear) can be right on top of a grid point outside the arrays on the
    // curvilinear grid
    if (smackOnTop && point_in_proc(iNear, jNear, gCurv)) {
      //       if (!point_in_proc(iNear,jNear,gCurv))
      //         return false;

      if (smoothed)
        tau = mTopoGridExt(iNear, jNear, 1);
      else
        tau = mTopo(iNear, jNear, 1);
    } else {
      computeNearestLowGridPoint(i, j, k, g, X0, Y0, Z0);

      //      if ( !( point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv) &&
      //      point_in_proc(i,j+1,gCurv) &&
      //              point_in_proc(i+1,j+1,gCurv) ) )
      //        return false;

      // linear interpolation to define topography between grid points
      //      double Qi, Qip1, Rj, Rjp1;
      //      Qi = (i+1 - q);
      //      Qip1 = (q - i);
      //      Rj = (j+1 - r);
      //      Rjp1 = (r - j);
      //      if (smoothed)
      //	tau = mTopoGrid(i,j,1)*Rj*Qi + mTopoGrid(i,j+1,1)*Rjp1*Qi +
      // mTopoGrid(i+1,j,1)*Rj*Qip1 + 	  mTopoGrid(i+1,j+1,1)*Rjp1*Qip1;
      //      else
      //	tau = mTopo(i,j,1)*Rj*Qi + mTopo(i,j+1,1)*Rjp1*Qi +
      // mTopo(i+1,j,1)*Rj*Qip1 + 	  mTopo(i+1,j+1,1)*Rjp1*Qip1;
      //    }
      if (point_in_proc_ext(i - 3, j - 3, gCurv) &&
          point_in_proc_ext(i + 4, j + 4, gCurv)) {
        float_sw4 a6cofi[8], a6cofj[8];
        gettopowgh(q - i, a6cofi);
        gettopowgh(r - j, a6cofj);
        tau = 0;
        for (int l = j - 3; l <= j + 4; l++)
          for (int k = i - 3; k <= i + 4; k++)
            tau +=
                a6cofi[k - i + 3] * a6cofj[l - j + 3] * mTopoGridExt(k, l, 1);
      } else {
        cout << "interpolate_topography: point not in ext domain i= " << i
             << " j = " << j << " limits " << m_iStart[gCurv] << " "
             << m_iEnd[gCurv] << " " << m_jStart[gCurv] << " " << m_jEnd[gCurv]
             << " nearest i = " << iNear << " j= " << jNear << endl;
        return -14;
      }
    }
  }  // end general case: interpolating mTopoGrid array
  Z0 = -tau;
  return 1;
}

//-----------------------------------------------------------------------
void EW::metric_derivatives_test() {
  // Assumes mMetric and mJ have been computed by numerical differentiation
  // This function computes corresponding expressions by analytical
  // differentiation

  int g = mNumberOfGrids - 1;
  Sarray metex(mMetric[g]), jacex(mJ[g]);  // TEMPORARY FIX

  int Bx = m_iStart[g];
  int By = m_jStart[g];
  int Bz = m_kStart[g];
  int Nx = m_iEnd[g];
  int Ny = m_jEnd[g];
  int Nz = m_kEnd[g];

  //   int nxg = m_global_nx[g];
  //   int nyg = m_global_ny[g];
  //   int nzg = m_global_nz[g];
  float_sw4 h = mGridSize[g];
  //   float_sw4 zmax = m_zmin[g-1] - (nzg-1)*h*(1-m_zetaBreak);

  // FTNC   if( m_croutines )
  m_gridGenerator->exact_metric(this, g, metex, jacex);
  //      metricexgh_ci( Bx, Nx, By, Ny, Bz, Nz, nzg, mX[g].c_ptr(),
  //      mY[g].c_ptr(), mZ[g].c_ptr(),
  //				    metex.c_ptr(), jacex.c_ptr(),
  // m_grid_interpolation_order, m_zetaBreak, zmax,
  // m_GaussianAmp, m_GaussianXc, m_GaussianYc, m_GaussianLx, m_GaussianLy );
  // FTNC   else
  // FTNC      metricexgh( &Bx, &Nx, &By, &Ny, &Bz, &Nz, &nxg, &nyg, &nzg,
  // mX[g].c_ptr(), mY[g].c_ptr(), mZ[g].c_ptr(), FTNC
  // metex.c_ptr(), jacex.c_ptr(), &m_grid_interpolation_order, &m_zetaBreak,
  // &zmax, FTNC				    &m_GaussianAmp,
  // &m_GaussianXc, &m_GaussianYc, &m_GaussianLx, &m_GaussianLy );
  communicate_array(metex, mNumberOfGrids - 1);
  communicate_array(jacex, mNumberOfGrids - 1);

  float_sw4 li[5], l2[5];
  int imin = m_iStartInt[g];
  int imax = m_iEndInt[g];
  int jmin = m_jStartInt[g];
  int jmax = m_jEndInt[g];
  int kmin = m_kStartInt[g];
  int kmax = m_kEndInt[g];

  // FTNC   if( m_croutines )
  meterr4c_ci(Bx, Nx, By, Ny, Bz, Nz, mMetric[g].c_ptr(), metex.c_ptr(),
              mJ[g].c_ptr(), jacex.c_ptr(), li, l2, imin, imax, jmin, jmax,
              kmin, kmax, h);
  // FTNC   else
  // FTNC      meterr4c( &Bx, &Nx, &By, &Ny, &Bz, &Nz, mMetric[g].c_ptr(),
  // metex.c_ptr(), mJ[g].c_ptr(), FTNC		jacex.c_ptr(), li, l2, &imin,
  // &imax, &jmin, &jmax, &kmin, &kmax, &h );

  float_sw4 tmp[5];
  for (int c = 0; c < 5; c++) tmp[c] = li[c];
  MPI_Allreduce(tmp, li, 5, m_mpifloat, MPI_MAX, m_cartesian_communicator);
  for (int c = 0; c < 5; c++) tmp[c] = l2[c];
  MPI_Allreduce(tmp, l2, 5, m_mpifloat, MPI_SUM, m_cartesian_communicator);
  for (int c = 0; c < 5; c++) l2[c] = sqrt(l2[c]);
  if (proc_zero()) {
    cout << "Errors in metric, max norm and L2 norm \n";
    for (int c = 0; c < 4; c++) cout << " " << li[c] << " " << l2[c] << endl;
    cout << "Error in Jacobian, max norm and L2 norm \n";
    cout << " " << li[4] << " " << l2[4] << endl;
  }
}

//-----------------------------------------------------------------------
void EW::copy_topo_to_topogridext() {
  if (topographyExists()) {
    int gTop = mNumberOfGrids - 1;
// copy raw topography
#pragma omp parallel for
    for (int i = m_iStart[gTop]; i <= m_iEnd[gTop]; ++i)
      for (int j = m_jStart[gTop]; j <= m_jEnd[gTop]; ++j)
        mTopoGridExt(i, j, 1) = mTopo(i, j, 1);

    int imin = mTopoGridExt.m_ib;
    int imax = mTopoGridExt.m_ie;
    int jmin = mTopoGridExt.m_jb;
    int jmax = mTopoGridExt.m_je;
    // Number of extra ghost points = m_ext_ghost_points
    int egh = m_iStart[gTop] - imin;

// Update extra ghost points. Do not worry about processor boundaries,
// they will be overwritten with correct values in the communication update
// afterward.
#pragma omp parallel for
    for (int i = imin + egh; i <= imax - egh; i++)
      for (int q = 0; q < egh; q++) {
        mTopoGridExt(i, jmin + q, 1) = mTopoGridExt(i, jmin + egh, 1);
        mTopoGridExt(i, jmax - q, 1) = mTopoGridExt(i, jmax - egh, 1);
      }
#pragma omp parallel for
    for (int j = jmin; j <= jmax; j++)
      for (int q = 0; q < egh; q++) {
        mTopoGridExt(imin + q, j, 1) = mTopoGridExt(imin + egh, j, 1);
        mTopoGridExt(imax - q, j, 1) = mTopoGridExt(imax - egh, j, 1);
      }
    communicate_array_2d_ext(mTopoGridExt);
  }
}

//-----------------------------------------------------------------------
void EW::gettopowgh(float_sw4 ai, float_sw4 wgh[8]) const {
  float_sw4 pol = ai * ai * ai * ai * ai * ai * ai *
                  (-251 + 135 * ai + 25 * ai * ai - 33 * ai * ai * ai +
                   6 * ai * ai * ai * ai) /
                  720;
  wgh[0] = -1.0 / 60 * ai + 1.0 / 180 * ai * ai + 1.0 / 48 * ai * ai * ai +
           23.0 / 144 * ai * ai * ai * ai -
           (17.0 * ai + 223.0) * ai * ai * ai * ai * ai / 720 - pol;
  wgh[1] = 3.0 / 20 * ai - 3.0 / 40 * ai * ai - 1.0 / 6 * ai * ai * ai -
           13.0 / 12 * ai * ai * ai * ai + 97.0 / 45 * ai * ai * ai * ai * ai +
           1.0 / 6 * ai * ai * ai * ai * ai * ai + 7 * pol;
  wgh[2] = -0.75 * ai + 0.75 * ai * ai + (13.0 + 155 * ai) * ai * ai * ai / 48 -
           103.0 / 16 * ai * ai * ai * ai * ai -
           121.0 / 240 * ai * ai * ai * ai * ai * ai - 21 * pol;
  wgh[3] = 1 - 49.0 / 36 * ai * ai - 49.0 / 9 * ai * ai * ai * ai +
           385.0 / 36 * ai * ai * ai * ai * ai +
           61.0 / 72 * ai * ai * ai * ai * ai * ai + 35 * pol;
  wgh[4] = 0.75 * ai + 0.75 * ai * ai - 13.0 / 48 * ai * ai * ai +
           89.0 / 16 * ai * ai * ai * ai -
           1537.0 / 144 * ai * ai * ai * ai * ai -
           41.0 / 48 * ai * ai * ai * ai * ai * ai - 35 * pol;
  wgh[5] = -3.0 / 20 * ai - 3.0 / 40 * ai * ai + 1.0 / 6 * ai * ai * ai -
           41.0 / 12 * ai * ai * ai * ai + 6.4 * ai * ai * ai * ai * ai +
           31.0 / 60 * ai * ai * ai * ai * ai * ai + 21 * pol;
  wgh[6] = 1.0 / 60 * ai + 1.0 / 180 * ai * ai - 1.0 / 48 * ai * ai * ai +
           167.0 / 144 * ai * ai * ai * ai -
           1537.0 / 720 * ai * ai * ai * ai * ai -
           25.0 / 144 * ai * ai * ai * ai * ai * ai - 7 * pol;
  wgh[7] = -1.0 / 6 * ai * ai * ai * ai + 11.0 / 36 * ai * ai * ai * ai * ai +
           1.0 / 40 * ai * ai * ai * ai * ai * ai + pol;
}

//-----------------------------------------------------------------------
void EW::smooth_grid(int maxIter) {
  // Smooth the grid (only the Z component for now)
  // NOTE: the current smoothing algorithm makes the error larger rather than
  // smaller!
  float_sw4 rf = 0.05;  // rf<1/6 for stability
  if (mVerbose >= 1 && proc_zero() && maxIter > 0)
    cout << "***smoothing the grid with " << maxIter
         << " Jacobi iterations and relaxation factor " << rf << " ***" << endl;

  int topLevel = mNumberOfGrids - 1;
  int g = mNumberOfGrids - 1;
  int iter;
  // temporary storage: How can I use mJ for temporary storage?
  Sarray tmp;
  tmp.define(m_iStart[topLevel], m_iEnd[topLevel], m_jStart[topLevel],
             m_jEnd[topLevel], m_kStart[topLevel], m_kEnd[topLevel]);

// initialize to make the Dirichlet boundary conditions work
#pragma omp parallel for
  for (int k = m_kStart[topLevel]; k <= m_kEnd[topLevel]; k++)
    for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
      for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i) {
        tmp(i, j, k) = mZ[g](i, j, k);
      }

  // Laplacian filter
  for (iter = 0; iter < maxIter; iter++) {
// loop over all interior points
#pragma omp parallel for
    for (int k = m_kStart[topLevel] + m_ghost_points + 1;
         k <= m_kEnd[topLevel] - m_ghost_points - 2; k++)
      for (int j = m_jStart[topLevel] + 1; j <= m_jEnd[topLevel] - 1; j++)
        for (int i = m_iStart[topLevel] + 1; i <= m_iEnd[topLevel] - 1; i++) {
          tmp(i, j, k) =
              mZ[g](i, j, k) + rf * (mZ[g](i + 1, j, k) + mZ[g](i - 1, j, k) +
                                     mZ[g](i, j + 1, k) + mZ[g](i, j - 1, k) +
                                     mZ[g](i, j, k + 1) + mZ[g](i, j, k - 1) -
                                     6. * mZ[g](i, j, k));
        }

// impose Neumann bc on the i and j sides
#pragma omp parallel for
    for (int k = m_kStart[topLevel] + m_ghost_points + 1;
         k <= m_kEnd[topLevel] - m_ghost_points - 2; k++) {
      for (int j = m_jStart[topLevel] + 1; j <= m_jEnd[topLevel] - 1; ++j) {
        int i = m_iStart[topLevel];
        tmp(i, j, k) = tmp(i + 1, j, k);
        i = m_iEnd[topLevel];
        tmp(i, j, k) = tmp(i - 1, j, k);
      }

      for (int i = m_iStart[topLevel] + 1; i <= m_iEnd[topLevel] - 1; ++i) {
        int j = m_jStart[topLevel];
        tmp(i, j, k) = tmp(i, j + 1, k);
        j = m_jEnd[topLevel];
        tmp(i, j, k) = tmp(i, j - 1, k);
      }
      // Corners
      int i = m_iStart[topLevel];
      int j = m_jStart[topLevel];
      tmp(i, j, k) = tmp(i + 1, j + 1, k);

      i = m_iEnd[topLevel];
      j = m_jStart[topLevel];
      tmp(i, j, k) = tmp(i - 1, j + 1, k);

      i = m_iStart[topLevel];
      j = m_jEnd[topLevel];
      tmp(i, j, k) = tmp(i + 1, j - 1, k);

      i = m_iEnd[topLevel];
      j = m_jEnd[topLevel];
      tmp(i, j, k) = tmp(i - 1, j - 1, k);
    }  // end Neumann loop

    // communicate parallel ghost points
    communicate_array(tmp, topLevel);

// update solution (Dirichlet are imposed implicitly by never changing the tmp
// array along the top or bottom boundary)
#pragma omp parallel for
    for (int k = m_kStart[topLevel]; k <= m_kEnd[topLevel]; k++)
      for (int j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
        for (int i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
          mZ[g](i, j, k) = tmp(i, j, k);

  }  // end for iter (grid smoother)
}
