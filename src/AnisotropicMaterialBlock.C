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
#include <iostream>

#include "AnisotropicMaterialBlock.h"
#include "EW.h"

using namespace std;
//-----------------------------------------------------------------------
AnisotropicMaterialBlock::AnisotropicMaterialBlock(
    EW* a_ew, float_sw4 rho, float_sw4 c[21], float_sw4 xmin, float_sw4 xmax,
    float_sw4 ymin, float_sw4 ymax, float_sw4 zmin, float_sw4 zmax) {
  m_rho = rho;
  for (int i = 0; i < 21; i++) {
    m_c[i] = c[i];
    m_cgrad[i] = 0;
  }
  m_xmin = xmin;
  m_xmax = xmax;
  m_ymin = ymin;
  m_ymax = ymax;
  m_zmin = zmin;
  m_zmax = zmax;
  m_tol = 1e-5;
  m_absoluteDepth = false;
  mEW = a_ew;
  float_sw4 bbox[6];
  mEW->getGlobalBoundingBox(bbox);

  // THE FOLLOWING ONLY WORKS IF m_absoluteDepth == true, or in the absence of
  // topography
  if (!mEW->topographyExists()) {
    // does this block cover all points?
    // note that the global_zmin can be non-zero when topography is present
    // global zmin is 0 in the absense of topography, and is assigned by the
    // grid generator when there is topography
    if (xmin > bbox[0] || ymin > bbox[2] || zmin > bbox[4] || xmax < bbox[1] ||
        ymax < bbox[3] || zmax < bbox[5]) {
      mCoversAllPoints = false;
      // tmp
      //     cout << "This block does NOT cover all grid points" << endl;
    } else {
      mCoversAllPoints = true;
      // tmp
      //     cout << "This block COVERS all grid points" << endl;
    }
  } else  // AP: with topography and m_absoluteDpeth == false, it is hard to say
          // if this block covers all points
  {
    mCoversAllPoints = false;
  }
}

//-----------------------------------------------------------------------
void AnisotropicMaterialBlock::set_absoluteDepth(bool absDepth) {
  m_absoluteDepth = absDepth;
}

//-----------------------------------------------------------------------
void AnisotropicMaterialBlock::set_gradients(float_sw4 rhograd,
                                             float_sw4 cgrad[21]) {
  m_rhograd = rhograd;
  for (int i = 0; i < 21; i++) m_cgrad[i] = cgrad[i];
}

//-----------------------------------------------------------------------
bool AnisotropicMaterialBlock::inside_block(float_sw4 x, float_sw4 y,
                                            float_sw4 z) {
  return m_xmin - m_tol <= x && x <= m_xmax + m_tol && m_ymin - m_tol <= y &&
         y <= m_ymax + m_tol && m_zmin - m_tol <= z && z <= m_zmax + m_tol;
}

//-----------------------------------------------------------------------
void AnisotropicMaterialBlock::set_material_properties(std::vector<Sarray>& rho,
                                                       std::vector<Sarray>& c) {
  //  int pc[4];
  // compute the number of parallel overlap points
  //  mEW->interiorPaddingCells( pc );
  int material = 0, outside = 0;

  for (int g = 0; g < mEW->mNumberOfCartesianGrids; g++)  // Cartesian grids
  {
#pragma omp parallel
    {
      // reference z-level for gradients is at z=0: AP changed this on 12/21/09
      float_sw4 zsurf = 0.;  // ?
#pragma omp for reduction(+ : material, outside)
      for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; k++) {
        for (int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; j++) {
          for (int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g]; i++) {
            float_sw4 x = (i - 1) * mEW->mGridSize[g];
            float_sw4 y = (j - 1) * mEW->mGridSize[g];
            float_sw4 z = mEW->m_zmin[g] + (k - 1) * mEW->mGridSize[g];

            // printf("x ,y,z %f %f %f %f\n",x,y,z,mEW->m_zmin[g]);

            float_sw4 depth;
            if (m_absoluteDepth) {
              depth = z;
            } else {
              mEW->getDepth(x, y, z, depth);
            }

            if (inside_block(x, y, depth)) {
              if (m_rho != -1)
                rho[g](i, j, k) = m_rho + m_rhograd * (depth - zsurf);
              for (int nr = 1; nr <= 21; nr++)
                c[g](nr, i, j, k) =
                    m_c[nr - 1] + m_cgrad[nr - 1] * (depth - zsurf);
              material++;
            } else {
              outside++;
              if (mEW->getVerbosity() > 5) {
                printf(
                    "Point (i,j,k)=(%i, %i, %i) in grid g=%i\n"
                    "with (x,y,z)=(%e,%e,%e) and depth=%e\n"
                    "is outside the block domain: %e<= x <= %e, %e <= y <= %e, "
                    "%e <= depth <= %e\n",
                    i, j, k, g, x, y, z, depth, m_xmin, m_xmax, m_ymin, m_ymax,
                    m_zmin, m_zmax);
              }
            }
          }
        }
      }
    }
    // communicate material properties to ghost points (necessary on refined
    // meshes because ghost points don't have a well defined depth/topography)
    mEW->communicate_array(rho[g], g);
    mEW->communicate_array(c[g], g);
  }  // end for all Cartesian grids

  if (mEW->topographyExists())  // curvilinear grid
  {
    int gTop = mEW->mNumberOfGrids - 1;
    for (int g = mEW->mNumberOfCartesianGrids; g < mEW->mNumberOfGrids; g++) {
//    int g = mEW->mNumberOfGrids-1;

// reference z-level for gradients is at z=0: AP changed this on 12/21/09
#pragma omp parallel
      {
        float_sw4 zsurf = 0.;

#pragma omp for reduction(+ : material, outside)
        for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; k++) {
          for (int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; j++) {
            for (int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g]; i++) {
              float_sw4 x = mEW->mX[g](i, j, k);
              float_sw4 y = mEW->mY[g](i, j, k);
              float_sw4 z = mEW->mZ[g](i, j, k);

              // printf("x ,y,z %f %f %f %f\n",x,y,z,mEW->m_zmin[g]);

              float_sw4 depth;
              if (m_absoluteDepth) {
                depth = z;
              } else {
                mEW->getDepth(x, y, z, depth);
              }

              if (inside_block(x, y, depth)) {
                if (m_rho != -1)
                  rho[g](i, j, k) = m_rho + m_rhograd * (depth - zsurf);
                for (int nr = 1; nr <= 21; nr++)
                  c[g](nr, i, j, k) =
                      m_c[nr - 1] + m_cgrad[nr - 1] * (depth - zsurf);
                material++;
              } else {
                if (mEW->getVerbosity() > 2) {
                  printf(
                      "Point (i,j,k)=(%i, %i, %i) in grid g=%i\n"
                      "with (x,y,z)=(%e,%e,%e) and depth=%e\n"
                      "is outside the block domain: %e<= x <= %e, %e <= y <= "
                      "%e, %e <= depth <= %e\n",
                      i, j, k, g, x, y, z, depth, m_xmin, m_xmax, m_ymin,
                      m_ymax, m_zmin, m_zmax);
                }
                outside++;
              }  // end else
            }    // end for i
          }      // end for j
        }        // end for k
      }          // end pragma omp
    }            // end for g (curvilinear)

  }  // end if topographyExists
  int outsideSum, materialSum;
  MPI_Reduce(&outside, &outsideSum, 1, MPI_INT, MPI_SUM, 0,
             mEW->m_cartesian_communicator);
  MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0,
             mEW->m_cartesian_communicator);

  if (mEW->proc_zero())
    cout << "block command: outside = " << outsideSum << ", "
         << "material = " << materialSum << endl;
}
