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
#include "Mspace.h"
#include "caliper.h"
#include "policies.h"
//-----------------------------------------------------------------------
void EW::bcfort_ci(int ib, int ie, int jb, int je, int kb, int ke, int wind[36],
                   int nx, int ny, int nz, float_sw4* u, float_sw4 h,
                   boundaryConditionType bccnd[6], float_sw4 sbop[5],
                   float_sw4* mu, float_sw4* la, float_sw4 t,
                   float_sw4* bforce1, float_sw4* bforce2, float_sw4* bforce3,
                   float_sw4* bforce4, float_sw4* bforce5, float_sw4* bforce6,
                   float_sw4 om, float_sw4 ph, float_sw4 cv, int curvilinear) {
  SW4_MARK_FUNCTION;
  const float_sw4 d4a = 2.0 / 3.0;
  const float_sw4 d4b = -1.0 / 12.0;
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  const size_t npts =
      static_cast<size_t>((ie - ib + 1)) * (je - jb + 1) * (ke - kb + 1);
  for (int s = 0; s < 6; s++) {
    if (bccnd[s] == bDirichlet || bccnd[s] == bSuperGrid) {
      size_t idel = 1 + wind[1 + 6 * s] - wind[6 * s];
      size_t ijdel = idel * (1 + wind[3 + 6 * s] - wind[2 + 6 * s]);
      if (s == 0) {
        //#pragma omp parallel for
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce1[3 * qq];
              u[ind + npts] = bforce1[1 + 3 * qq];
              u[ind + 2 * npts] = bforce1[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 1) {
        //#pragma omp parallel for
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce2[3 * qq];
              u[ind + npts] = bforce2[1 + 3 * qq];
              u[ind + 2 * npts] = bforce2[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 2) {
        //#pragma omp parallel for
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce3[3 * qq];
              u[ind + npts] = bforce3[1 + 3 * qq];
              u[ind + 2 * npts] = bforce3[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 3) {
        //#pragma omp parallel for
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce4[3 * qq];
              u[ind + npts] = bforce4[1 + 3 * qq];
              u[ind + 2 * npts] = bforce4[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 4) {
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          //#pragma omp parallel for
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce5[3 * qq];
              u[ind + npts] = bforce5[1 + 3 * qq];
              u[ind + 2 * npts] = bforce5[2 + 3 * qq];
              qq++;
            }
          }
        }
      } else if (s == 5) {
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++) {
          size_t qq = (k - wind[4 + 6 * s]) * ijdel;
          //#pragma omp parallel for
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++) {
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              u[ind] = bforce6[3 * qq];
              u[ind + npts] = bforce6[1 + 3 * qq];
              u[ind + 2 * npts] = bforce6[2 + 3 * qq];
              qq++;
            }
          }
        }
      }
    } else if (bccnd[s] == bPeriodic) {
      if (s == 0) {
#pragma omp parallel for
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + nx;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 1) {
#pragma omp parallel for
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - nx;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 2) {
#pragma omp parallel for
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + ni * ny;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 3) {
#pragma omp parallel for
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - ni * ny;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 4) {
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
#pragma omp parallel for
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + nij * nz;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      } else if (s == 5) {
        for (int k = wind[4 + 6 * s]; k <= wind[5 + 6 * s]; k++)
#pragma omp parallel for
          for (int j = wind[2 + 6 * s]; j <= wind[3 + 6 * s]; j++)
            for (int i = wind[6 * s]; i <= wind[1 + 6 * s]; i++) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - nij * nz;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            }
      }
    } else if (bccnd[s] == bStressFree) {
      REQUIRE2(s == 4 || s == 5, "EW::bcfort_ci,  ERROR: Free surface condition"
                                     << " not implemented for side " << s
                                     << endl);
      if (s == 4 && curvilinear == 0) {
        int k = 1, kl = 1;
#pragma omp parallel for
        for (int j = jb + 2; j <= je - 2; j++)
          for (int i = ib + 2; i <= ie - 2; i++) {
            size_t qq = i - ib + ni * (j - jb);
            size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
            float_sw4 wx =
                d4a * (u[2 * npts + ind + 1] - u[2 * npts + ind - 1]) +
                d4b * (u[2 * npts + ind + 2] - u[2 * npts + ind - 2]);
            float_sw4 ux = d4a * (u[ind + 1] - u[ind - 1]) +
                           d4b * (u[ind + 2] - u[ind - 2]);
            float_sw4 wy =
                d4a * (u[2 * npts + ind + ni] - u[2 * npts + ind - ni]) +
                d4b * (u[2 * npts + ind + 2 * ni] - u[2 * npts + ind - 2 * ni]);
            float_sw4 vy =
                d4a * (u[npts + ind + ni] - u[npts + ind - ni]) +
                d4b * (u[npts + ind + 2 * ni] - u[npts + ind - 2 * ni]);
            float_sw4 uz = 0, vz = 0, wz = 0;
            for (int w = 1; w <= 4; w++) {
              uz += sbop[w] * u[ind + nij * kl * (w - 1)];
              vz += sbop[w] * u[npts + ind + nij * kl * (w - 1)];
              wz += sbop[w] * u[2 * npts + ind + nij * kl * (w - 1)];
            }
            u[ind - nij * kl] =
                (-uz - kl * wx + kl * h * bforce5[3 * qq] / mu[ind]) / sbop[0];
            u[npts + ind - nij * kl] =
                (-vz - kl * wy + kl * h * bforce5[1 + 3 * qq] / mu[ind]) /
                sbop[0];
            u[2 * npts + ind - nij * kl] =
                (-wz +
                 (-kl * la[ind] * (ux + vy) + kl * h * bforce5[2 + 3 * qq]) /
                     (2 * mu[ind] + la[ind])) /
                sbop[0];
          }
      } else if (s == 5 && curvilinear == 0) {
        int k = nz, kl = -1;
#pragma omp parallel for
        for (int j = jb + 2; j <= je - 2; j++)
          for (int i = ib + 2; i <= ie - 2; i++) {
            size_t qq = i - ib + ni * (j - jb);
            size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
            float_sw4 wx =
                d4a * (u[2 * npts + ind + 1] - u[2 * npts + ind - 1]) +
                d4b * (u[2 * npts + ind + 2] - u[2 * npts + ind - 2]);
            float_sw4 ux = d4a * (u[ind + 1] - u[ind - 1]) +
                           d4b * (u[ind + 2] - u[ind - 2]);
            float_sw4 wy =
                d4a * (u[2 * npts + ind + ni] - u[2 * npts + ind - ni]) +
                d4b * (u[2 * npts + ind + 2 * ni] - u[2 * npts + ind - 2 * ni]);
            float_sw4 vy =
                d4a * (u[npts + ind + ni] - u[npts + ind - ni]) +
                d4b * (u[npts + ind + 2 * ni] - u[npts + ind - 2 * ni]);
            float_sw4 uz = 0, vz = 0, wz = 0;
            for (int w = 1; w <= 4; w++) {
              uz += sbop[w] * u[ind + nij * kl * (w - 1)];
              vz += sbop[w] * u[npts + ind + nij * kl * (w - 1)];
              wz += sbop[w] * u[2 * npts + ind + nij * kl * (w - 1)];
            }
            u[ind - nij * kl] =
                (-uz - kl * wx + kl * h * bforce6[3 * qq] / mu[ind]) / sbop[0];
            u[npts + ind - nij * kl] =
                (-vz - kl * wy + kl * h * bforce6[1 + 3 * qq] / mu[ind]) /
                sbop[0];
            u[2 * npts + ind - nij * kl] =
                (-wz +
                 (-kl * la[ind] * (ux + vy) + kl * h * bforce6[2 + 3 * qq]) /
                     (2 * mu[ind] + la[ind])) /
                sbop[0];
          }
      }
    }
  }
}

//-----------------------------------------------------------------------
void EW::bcfortsg_ci(int ib, int ie, int jb, int je, int kb, int ke,
                     int wind[36], int nx, int ny, int nz, float_sw4* u,
                     float_sw4 h, boundaryConditionType bccnd[6],
                     float_sw4 sbop[5], float_sw4* mu, float_sw4* la,
                     float_sw4 t, float_sw4* bforce1, float_sw4* bforce2,
                     float_sw4* bforce3, float_sw4* bforce4, float_sw4* bforce5,
                     float_sw4* bforce6, float_sw4 om, float_sw4 ph,
                     float_sw4 cv, float_sw4* strx, float_sw4* stry) {
  SW4_MARK_FUNCTION;
  const float_sw4 d4a = 2.0 / 3.0;
  const float_sw4 d4b = -1.0 / 12.0;
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  const size_t npts =
      static_cast<size_t>((ie - ib + 1)) * (je - jb + 1) * (ke - kb + 1);
#ifdef ENABLE_CUDA

#if SW4_RAJA_VERSION == 6
  using BCFORT_EXEC_POL2_ASYNC =
      RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::For<
          0, RAJA::cuda_threadblock_exec<4>,
          RAJA::statement::For<
              1, RAJA::cuda_threadblock_exec<4>,
              RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
                                   RAJA::statement::Lambda<0>>>>>>;

#elif SW4_RAJA_VERSION == 7

  using BCFORT_EXEC_POL2_ASYNC =
      RAJA::KernelPolicy<RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
          0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
          RAJA::statement::Tile<
              1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
              RAJA::statement::Tile<
                  2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
                  RAJA::statement::For<
                      0, RAJA::cuda_thread_y_direct,
                      RAJA::statement::For<
                          1, RAJA::cuda_thread_x_direct,
                          RAJA::statement::For<
                              2, RAJA::cuda_thread_z_direct,
                              RAJA::statement::Lambda<0>>>>>>>>>;

#endif

  // Policy below produces much lower GPU faults: 700 instead of 1000 but no
  // real difference in runtime. The faults could be coming from some code in
  // the RAJA nested loops
  using BCFORT_EXEC_POL2_X =
      RAJA::KernelPolicy<RAJA::statement::CudaKernel<RAJA::statement::For<
          2, RAJA::cuda_block_exec,
          RAJA::statement::For<
              1, RAJA::cuda_block_exec,
              RAJA::statement::For<0, RAJA::cuda_thread_exec,
                                   RAJA::statement::Lambda<0>>>>>>;
#else
  using BCFORT_EXEC_POL2_ASYNC = DEFAULT_LOOP3;
#endif
  for (int s = 0; s < 6; s++) {
    if (bccnd[s] == bDirichlet || bccnd[s] == bSuperGrid) {
      // std::cout<<"SET 1 "<<s<<"\n";
      // size_t idel = 1+wind[1+6*s]-wind[6*s];
      // size_t ijdel = idel * (1+wind[3+6*s]-wind[2+6*s]);
      if (s == 0) {
        // PREFETCH(bforce1);
        int lni = wind[1 + 6 * s] - wind[6 * s] + 1;
        int lnij = (wind[3 + 6 * s] - wind[2 + 6 * s] + 1) * lni;
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        //#pragma omp parallel for
        // for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ ) {
        //    size_t qq = (k-wind[4+6*s])*ijdel;
        //    for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ ) {
        // 	  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ ) {
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t qq =
                  (k - kstart) * lnij + (j - jstart) * lni + (i - istart);
              u[ind] = bforce1[3 * qq];
              u[ind + npts] = bforce1[1 + 3 * qq];
              u[ind + 2 * npts] = bforce1[2 + 3 * qq];
              // qq++;
            });  // SYNC_STREAM;
        //}
        // }
      } else if (s == 1) {
        // PREFETCH(bforce2);
        int lni = wind[1 + 6 * s] - wind[6 * s] + 1;
        int lnij = (wind[3 + 6 * s] - wind[2 + 6 * s] + 1) * lni;
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        //#pragma omp parallel for
        // for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ ) {
        //    size_t qq = (k-wind[4+6*s])*ijdel;
        //    for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ ) {
        // 	  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ ) {
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t qq =
                  (k - kstart) * lnij + (j - jstart) * lni + (i - istart);
              u[ind] = bforce2[3 * qq];
              u[ind + npts] = bforce2[1 + 3 * qq];
              u[ind + 2 * npts] = bforce2[2 + 3 * qq];
              // qq++;
            });  // SYNC_STREAM;
        // }
        // }
      } else if (s == 2) {
        // PREFETCH(bforce3);
        int lni = wind[1 + 6 * s] - wind[6 * s] + 1;
        int lnij = (wind[3 + 6 * s] - wind[2 + 6 * s] + 1) * lni;
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        //#pragma omp parallel for
        // for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ ) {
        //    size_t qq = (k-wind[4+6*s])*ijdel;
        //    for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ ) {
        // 	  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ ) {
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t qq =
                  (k - kstart) * lnij + (j - jstart) * lni + (i - istart);
              u[ind] = bforce3[3 * qq];
              u[ind + npts] = bforce3[1 + 3 * qq];
              u[ind + 2 * npts] = bforce3[2 + 3 * qq];
              // qq++;
            });  // SYNC_STREAM;
        // }
        //  }
      } else if (s == 3) {
        // PREFETCH(bforce4);
        int lni = wind[1 + 6 * s] - wind[6 * s] + 1;
        int lnij = (wind[3 + 6 * s] - wind[2 + 6 * s] + 1) * lni;
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        //#pragma omp parallel for
        // for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ ) {
        //    size_t qq = (k-wind[4+6*s])*ijdel;
        //    for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ ) {
        // 	  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ ) {
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t qq =
                  (k - kstart) * lnij + (j - jstart) * lni + (i - istart);
              u[ind] = bforce4[3 * qq];
              u[ind + npts] = bforce4[1 + 3 * qq];
              u[ind + 2 * npts] = bforce4[2 + 3 * qq];
              // qq++;
            });  // SYNC_STREAM;
        // }
        //   }
      } else if (s == 4) {
        PREFETCH(bforce5);
        int lni = wind[1 + 6 * s] - wind[6 * s] + 1;
        int lnij = (wind[3 + 6 * s] - wind[2 + 6 * s] + 1) * lni;
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        // for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ ) {
        //    size_t qq = (k-wind[4+6*s])*ijdel;
        //    //#pragma omp parallel for
        //    for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ ) {
        // 	  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ ) {
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t qq =
                  (k - kstart) * lnij + (j - jstart) * lni + (i - istart);
              u[ind] = bforce5[3 * qq];
              u[ind + npts] = bforce5[1 + 3 * qq];
              u[ind + 2 * npts] = bforce5[2 + 3 * qq];
              // qq++;
            });  // SYNC_STREAM;
        // }
        // }
      } else if (s == 5) {
        // PREFETCH(bforce6);
        int lni = wind[1 + 6 * s] - wind[6 * s] + 1;
        int lnij = (wind[3 + 6 * s] - wind[2 + 6 * s] + 1) * lni;
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        // for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ ) {
        //    size_t qq = (k-wind[4+6*s])*ijdel;
        //    //#pragma omp parallel for
        //    for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ ) {
        // 	  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ ) {
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t qq =
                  (k - kstart) * lnij + (j - jstart) * lni + (i - istart);
              u[ind] = bforce6[3 * qq];
              u[ind + npts] = bforce6[1 + 3 * qq];
              u[ind + 2 * npts] = bforce6[2 + 3 * qq];
              // qq++;
            });  // SYNC_STREAM;
        // }
        //  }
      }
    } else if (bccnd[s] == bPeriodic) {
      // std::cout<<"SET 2\n";
      if (s == 0) {
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        //#pragma omp parallel for
        // for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
        //    for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
        // 	  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
        // 	  {
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + nx;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            });  // SYNC_STREAM;
      } else if (s == 1) {
        // #pragma omp parallel for
        // 	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
        // 	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
        // 		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
        // 		  {
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - nx;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            });  // SYNC_STREAM;
      } else if (s == 2) {
        // #pragma omp parallel for
        // 	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
        // 	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
        // 		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
        // 		  {
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + ni * ny;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            });  // SYNC_STREAM;
      } else if (s == 3) {
        // #pragma omp parallel for
        // 	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
        // 	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
        // 		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
        // 		  {
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - ni * ny;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            });  // SYNC_STREAM;
      } else if (s == 4) {
        // 	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
        // #pragma omp parallel for
        // 	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
        // 		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
        // 		  {
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind + nij * nz;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            });  // SYNC_STREAM;
      } else if (s == 5) {
        // 	    for( int k=wind[4+6*s]; k <= wind[5+6*s] ; k++ )
        // #pragma omp parallel for
        // 	       for( int j=wind[2+6*s]; j <= wind[3+6*s] ; j++ )
        // 		  for( int i=wind[6*s]; i <= wind[1+6*s] ; i++ )
        // 		  {
        int istart = wind[6 * s];
        int jstart = wind[2 + 6 * s];
        int kstart = wind[4 + 6 * s];
        RAJA::RangeSegment k_range(kstart, wind[5 + 6 * s] + 1);
        RAJA::RangeSegment j_range(jstart, wind[3 + 6 * s] + 1);
        RAJA::RangeSegment i_range(istart, wind[1 + 6 * s] + 1);
        RAJA::kernel<BCFORT_EXEC_POL2_ASYNC>(
            RAJA::make_tuple(k_range, j_range, i_range),
            [=] RAJA_DEVICE(int k, int j, int i) {
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              size_t indp = ind - nij * nz;
              u[ind] = u[indp];
              u[ind + npts] = u[indp + npts];
              u[ind + 2 * npts] = u[indp + 2 * npts];
            });  // SYNC_STREAM;
      }
    } else if (bccnd[s] == bStressFree) {
      // std::cout<<"SET 3 \n";
      REQUIRE2(s == 4 || s == 5, "EW::bcfort_ci,  ERROR: Free surface condition"
                                     << " not implemented for side " << s
                                     << endl);
#ifdef ENABLE_CUDA

#if SW4_RAJA_VERSION == 6
      using BCFORT_EXEC_POL3_ASYNC = RAJA::KernelPolicy<
          RAJA::statement::CudaKernelAsync<RAJA::statement::For<
              0, RAJA::cuda_threadblock_exec<16>,
              RAJA::statement::For<1, RAJA::cuda_threadblock_exec<16>,
                                   RAJA::statement::Lambda<0>>>>>;

#elif SW4_RAJA_VERSION == 7
      using BCFORT_EXEC_POL3_ASYNC = RAJA::KernelPolicy<
          RAJA::statement::CudaKernelAsync<RAJA::statement::Tile<
              0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
              RAJA::statement::Tile<
                  1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
                  RAJA::statement::For<
                      0, RAJA::cuda_thread_x_direct,
                      RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
                                           RAJA::statement::Lambda<0>>>>>>>;

#endif

#else
      using BCFORT_EXEC_POL3_ASYNC = DEFAULT_LOOP2;
#endif
      if (s == 4) {
        // PREFETCH(bforce5);
        int k = 1, kl = 1;
        //#pragma omp parallel for
        RAJA::RangeSegment i_range(ib + 2, ie - 1);
        RAJA::RangeSegment j_range(jb + 2, je - 1);
        // for( int j=jb+2 ; j <= je-2 ; j++ )
        //    for( int i=ib+2 ; i <= ie-2 ; i++ )
        //    {
        RAJA::kernel<BCFORT_EXEC_POL3_ASYNC>(
            RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
              size_t qq = i - ib + ni * (j - jb);
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              float_sw4 wx =
                  strx[i - ib] *
                  (d4a * (u[2 * npts + ind + 1] - u[2 * npts + ind - 1]) +
                   d4b * (u[2 * npts + ind + 2] - u[2 * npts + ind - 2]));
              float_sw4 ux = strx[i - ib] * (d4a * (u[ind + 1] - u[ind - 1]) +
                                             d4b * (u[ind + 2] - u[ind - 2]));
              float_sw4 wy =
                  stry[j - jb] *
                  (d4a * (u[2 * npts + ind + ni] - u[2 * npts + ind - ni]) +
                   d4b * (u[2 * npts + ind + 2 * ni] -
                          u[2 * npts + ind - 2 * ni]));
              float_sw4 vy =
                  stry[j - jb] *
                  (d4a * (u[npts + ind + ni] - u[npts + ind - ni]) +
                   d4b * (u[npts + ind + 2 * ni] - u[npts + ind - 2 * ni]));
              float_sw4 uz = 0, vz = 0, wz = 0;
              for (int w = 1; w <= 4; w++) {
                uz += sbop[w] * u[ind + nij * kl * (w - 1)];
                vz += sbop[w] * u[npts + ind + nij * kl * (w - 1)];
                wz += sbop[w] * u[2 * npts + ind + nij * kl * (w - 1)];
              }
              u[ind - nij * kl] =
                  (-uz - kl * wx + kl * h * bforce5[3 * qq] / mu[ind]) /
                  sbop[0];
              u[npts + ind - nij * kl] =
                  (-vz - kl * wy + kl * h * bforce5[1 + 3 * qq] / mu[ind]) /
                  sbop[0];
              u[2 * npts + ind - nij * kl] =
                  (-wz +
                   (-kl * la[ind] * (ux + vy) + kl * h * bforce5[2 + 3 * qq]) /
                       (2 * mu[ind] + la[ind])) /
                  sbop[0];
            });  // SYNC_STREAM;
      } else {
        int k = nz, kl = -1;
        // #pragma omp parallel for
        // 	    for( int j=jb+2 ; j <= je-2 ; j++ )
        // 	       for( int i=ib+2 ; i <= ie-2 ; i++ )
        // 	       {
        // PREFETCH(bforce6);
        RAJA::RangeSegment i_range(ib + 2, ie - 1);
        RAJA::RangeSegment j_range(jb + 2, je - 1);
        RAJA::kernel<BCFORT_EXEC_POL3_ASYNC>(
            RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
              size_t qq = i - ib + ni * (j - jb);
              size_t ind = i - ib + ni * (j - jb) + nij * (k - kb);
              float_sw4 wx =
                  strx[i - ib] *
                  (d4a * (u[2 * npts + ind + 1] - u[2 * npts + ind - 1]) +
                   d4b * (u[2 * npts + ind + 2] - u[2 * npts + ind - 2]));
              float_sw4 ux = strx[i - ib] * (d4a * (u[ind + 1] - u[ind - 1]) +
                                             d4b * (u[ind + 2] - u[ind - 2]));
              float_sw4 wy =
                  stry[j - jb] *
                  (d4a * (u[2 * npts + ind + ni] - u[2 * npts + ind - ni]) +
                   d4b * (u[2 * npts + ind + 2 * ni] -
                          u[2 * npts + ind - 2 * ni]));
              float_sw4 vy =
                  stry[j - jb] *
                  (d4a * (u[npts + ind + ni] - u[npts + ind - ni]) +
                   d4b * (u[npts + ind + 2 * ni] - u[npts + ind - 2 * ni]));
              float_sw4 uz = 0, vz = 0, wz = 0;
              for (int w = 1; w <= 4; w++) {
                uz += sbop[w] * u[ind + nij * kl * (w - 1)];
                vz += sbop[w] * u[npts + ind + nij * kl * (w - 1)];
                wz += sbop[w] * u[2 * npts + ind + nij * kl * (w - 1)];
              }
              u[ind - nij * kl] =
                  (-uz - kl * wx + kl * h * bforce6[3 * qq] / mu[ind]) /
                  sbop[0];
              u[npts + ind - nij * kl] =
                  (-vz - kl * wy + kl * h * bforce6[1 + 3 * qq] / mu[ind]) /
                  sbop[0];
              u[2 * npts + ind - nij * kl] =
                  (-wz +
                   (-kl * la[ind] * (ux + vy) + kl * h * bforce6[2 + 3 * qq]) /
                       (2 * mu[ind] + la[ind])) /
                  sbop[0];
            });  // SYNC_STREAM;
      }
    }
  }
  // SYNC_STREAM;
}

//-----------------------------------------------------------------------
void EW::twdirbdry_ci(int wind[6], float_sw4 h, float_sw4 t, float_sw4 om,
                      float_sw4 cv, float_sw4 ph, float_sw4* bforce,
                      float_sw4 zmin) {
  SW4_MARK_FUNCTION;
  // size_t qq=0;
  int ni = wind[1] - wind[0] + 1;
  int nij = (wind[3] - wind[2] + 1) * ni;
  // int klen = wind[5]-wind[4]+1;
  int istart = wind[0];
  int jstart = wind[2];
  int kstart = wind[4];
  RAJA::RangeSegment k_range(wind[4], wind[5] + 1);
  RAJA::RangeSegment j_range(wind[2], wind[3] + 1);
  RAJA::RangeSegment i_range(wind[0], wind[1] + 1);
  RAJA::kernel<BCFORT_EXEC_POL1>(
      RAJA::make_tuple(k_range, j_range, i_range),
      [=] RAJA_DEVICE(int k, int j, int i) {
        // #pragma omp parallel for // qq variable will not work with OpenMP
        // for( int k=wind[4]; k <= wind[5] ; k++ ) {
        //    for( int j=wind[2]; j <= wind[3] ; j++ ) {
        // 	 for( int i=wind[0]; i <= wind[1] ; i++ ) {
        double x = (i - 1) * h, y = (j - 1) * h, z = zmin + (k - 1) * h;
        size_t qq = (k - kstart) * nij + (j - jstart) * ni + (i - istart);
        // if (lqq!=qq) std::cout<<i<<" "<<j<<" "<<k<<" "<<lqq<<" "<<qq;
        bforce[3 * qq] =
            sin(om * (x - cv * t)) * sin(om * y + ph) * sin(om * z + ph);
        bforce[1 + 3 * qq] =
            sin(om * x + ph) * sin(om * (y - cv * t)) * sin(om * z + ph);
        bforce[2 + 3 * qq] =
            sin(om * x + ph) * sin(om * y + ph) * sin(om * (z - cv * t));
        // qq++;
      });
  SYNC_STREAM;
}

//-----------------------------------------------------------------------
void EW::twdirbdryc_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                       int klast, int wind[6], float_sw4 t, float_sw4 om,
                       float_sw4 cv, float_sw4 ph, float_sw4* bforce,
                       float_sw4* x, float_sw4* y, float_sw4* z) {
  SW4_MARK_FUNCTION;
  size_t ni = ilast - ifirst + 1;
  size_t nij = ni * (jlast - jfirst + 1);
  size_t qq = 0;
  // #pragma omp parallel for qq variable will not work with OpenMP
  for (int k = wind[4]; k <= wind[5]; k++) {
    for (int j = wind[2]; j <= wind[3]; j++) {
      for (int i = wind[0]; i <= wind[1]; i++) {
        size_t ind = i - ifirst + ni * (j - jfirst) + nij * (k - kfirst);
	qq=(i-wind[0])+(j-wind[2])*(wind[1]-wind[0]+1)+(k-wind[4])*(wind[1]-wind[0]+1)*(wind[3]-wind[2]+1);
        bforce[3 * qq] = sin(om * (x[ind] - cv * t)) * sin(om * y[ind] + ph) *
                         sin(om * z[ind] + ph);
        bforce[1 + 3 * qq] = sin(om * x[ind] + ph) *
                             sin(om * (y[ind] - cv * t)) *
                             sin(om * z[ind] + ph);
        bforce[2 + 3 * qq] = sin(om * x[ind] + ph) * sin(om * y[ind] + ph) *
                             sin(om * (z[ind] - cv * t));
        //qq++;
      }
    }
  }
}

//-----------------------------------------------------------------------
void EW::twfrsurfz_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                      int klast, float_sw4 h, int kz, float_sw4 t,
                      float_sw4 omega, float_sw4 c, float_sw4 phase,
                      float_sw4* bforce, float_sw4* mu, float_sw4* lambda,
                      float_sw4 zmin) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  double z = (kz - 1) * h + zmin;
  RAJA::RangeSegment j_range(jfirst, jlast + 1);
  RAJA::RangeSegment i_range(ifirst, ilast + 1);
  RAJA::kernel<BCFORT_EXEC_POL2>(
      RAJA::make_tuple(j_range, i_range), [=] RAJA_DEVICE(int j, int i) {
        // #pragma omp parallel for
        //    for( int j=jfirst ; j<= jlast ; j++ )
        //       for( int i=ifirst ; i<= ilast ; i++ )
        //       {
        size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
        size_t qq = (i - ifirst) + ni * (j - jfirst);
        double x = (i - 1) * h, y = (j - 1) * h;
        double forces[3];
        double t13, t15, t16, t19, t20, t21, t23, t24, t28, t29, t32, t33, t34,
            t37, t38, t43, t44, t49, t60, t62, t65;

        t13 = mu[ind];
        t15 = omega * x + phase;
        t16 = cos(t15);
        t19 = omega * y + phase;
        t20 = sin(t19);
        t21 = c * t;
        t23 = omega * (z - t21);
        t24 = sin(t23);
        t28 = omega * (x - t21);
        t29 = sin(t28);
        t32 = omega * z + phase;
        t33 = cos(t32);
        t34 = t33 * omega;
        forces[0] = t13 * (t16 * omega * t20 * t24 + t29 * t20 * t34);
        t37 = sin(t15);
        t38 = cos(t19);
        t43 = omega * (y - t21);
        t44 = sin(t43);
        forces[1] = t13 * (t37 * t38 * omega * t24 + t37 * t44 * t34);
        t49 = cos(t23);
        t60 = cos(t28);
        t62 = sin(t32);
        t65 = cos(t43);
        forces[2] =
            2 * t13 * t37 * t20 * t49 * omega +
            lambda[ind] * (t60 * omega * t20 * t62 + t37 * t65 * omega * t62 +
                           t37 * t20 * t49 * omega);

        bforce[3 * qq] = forces[0];
        bforce[1 + 3 * qq] = forces[1];
        bforce[2 + 3 * qq] = forces[2];
      });
  SYNC_STREAM;
}

//-----------------------------------------------------------------------
void EW::twfrsurfzatt_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, float_sw4 h, int kz,
                         float_sw4 t, float_sw4 omega, float_sw4 c,
                         float_sw4 phase, float_sw4* bforce, float_sw4* mua,
                         float_sw4* lambdaa, float_sw4 zmin) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  double z = (kz - 1) * h + zmin;
#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++)
    for (int i = ifirst; i <= ilast; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      double x = (i - 1) * h, y = (j - 1) * h;
      double t2, t3, t6, t7, t8, t11, t12, t17, t18, t27, t30;
      double t20, t31, t35, t40, t45, t50, t52, t54;
      double t16, t23, t24, t34;
      double forces[3];

      t2 = omega * x + phase;
      t3 = sin(t2);
      t6 = omega * y + phase;
      t7 = cos(t6);
      t8 = c * t;
      t11 = -omega * (z - t8) - phase;
      t12 = sin(t11);
      t16 = omega * (x - t8);
      t17 = -t16 - phase;
      t18 = cos(t17);
      t20 = t12 * omega;
      forces[0] = mua[ind] * (t3 * omega * t7 * t12 + t18 * t3 * t20);
      t23 = cos(t2);
      t24 = sin(t6);
      t27 = sin(t16);
      t30 = -omega * (y - t8) - phase;
      t31 = cos(t30);
      t34 = omega * z + phase;
      t35 = sin(t34);
      forces[1] = mua[ind] * (t23 * t24 * t20 - t27 * t31 * t35 * omega);
      t40 = cos(t11);
      t45 = sin(t17);
      t50 = t40 * omega;
      t52 = sin(t30);
      t54 = cos(t34);
      forces[2] = 2.0 * mua[ind] * t23 * t7 * t40 * omega +
                  lambdaa[ind] * (t45 * omega * t3 * t40 + t18 * t23 * t50 +
                                  t27 * t52 * omega * t54 + t23 * t7 * t50);

      bforce[3 * qq] = forces[0];
      bforce[1 + 3 * qq] = forces[1];
      bforce[2 + 3 * qq] = forces[2];
    }
}

//-----------------------------------------------------------------------
void EW::twfrsurfzsgstr_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4 h, int kz,
                           float_sw4 t, float_sw4 om, float_sw4 c, float_sw4 ph,
                           float_sw4 omstrx, float_sw4 omstry,
                           float_sw4* bforce, float_sw4* mu, float_sw4* lambda,
                           float_sw4 zmin) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  double z = (kz - 1) * h + zmin;
#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++)
    for (int i = ifirst; i <= ilast; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      double x = (i - 1) * h, y = (j - 1) * h;
      double t1, t10, t11, t12, t15, t17, t19, t20, t22, t24, t25, t29, t3, t31,
          t32, t36, t39, t4, t40, t46, t51, t53, t56, t6, t7;
      double forces[3];

      t1 = c * t;
      t3 = om * (x - t1);
      t4 = sin(t3);
      t6 = om * y + ph;
      t7 = sin(t6);
      t10 = om * z + ph;
      t11 = cos(t10);
      t12 = t11 * om;
      t15 = sin(omstrx * x);
      t17 = 1 + t15 / 2;
      t19 = om * x + ph;
      t20 = cos(t19);
      t22 = om * t7;
      t24 = om * (z - t1);
      t25 = sin(t24);
      forces[0] = mu[ind] * (t4 * t7 * t12 + t17 * t20 * t22 * t25);
      t29 = sin(t19);
      t31 = om * (y - t1);
      t32 = sin(t31);
      t36 = sin(omstry * y);
      t39 = (1 + t36 / 2) * t29;
      t40 = cos(t6);
      forces[1] = mu[ind] * (t29 * t32 * t12 + t39 * t40 * om * t25);
      t46 = cos(t24);
      t51 = cos(t3);
      t53 = sin(t10);
      t56 = cos(t31);
      forces[2] = 2 * mu[ind] * t29 * t7 * t46 * om +
                  lambda[ind] * (t17 * t51 * t22 * t53 + t39 * t56 * om * t53 +
                                 t29 * t7 * t46 * om);

      bforce[3 * qq] = forces[0];
      bforce[1 + 3 * qq] = forces[1];
      bforce[2 + 3 * qq] = forces[2];
    }
}

//-----------------------------------------------------------------------
void EW::twfrsurfzsgstratt_ci(int ifirst, int ilast, int jfirst, int jlast,
                              int kfirst, int klast, float_sw4 h, int kz,
                              float_sw4 t, float_sw4 omega, float_sw4 c,
                              float_sw4 phase, float_sw4 omstrx,
                              float_sw4 omstry, float_sw4* bforce,
                              float_sw4* mua, float_sw4* lambdaa,
                              float_sw4 zmin) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  float_sw4 z = (kz - 1) * h + zmin;
#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++)
    for (int i = ifirst; i <= ilast; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      float_sw4 x = (i - 1) * h, y = (j - 1) * h;
      float_sw4 t1, t12, t13, t17, t19, t22, t23, t28, t3, t31, t32, t35, t36,
          t4, t40, t42, t43, t45, t5, t51, t56, t61, t66, t68, t7, t8;
      float_sw4 forces[3];

      t1 = c * t;
      t3 = omega * (x - t1);
      t4 = -t3 - phase;
      t5 = cos(t4);
      t7 = omega * x + phase;
      t8 = sin(t7);
      t12 = -omega * (z - t1) - phase;
      t13 = sin(t12);
      t17 = sin(omstrx * x);
      t19 = 1 + t17 / 2;
      t22 = omega * y + phase;
      t23 = cos(t22);
      forces[0] =
          mua[ind] * (t5 * t8 * t13 * omega + t19 * t8 * omega * t23 * t13);
      t28 = sin(t3);
      t31 = -omega * (y - t1) - phase;
      t32 = cos(t31);
      t35 = omega * z + phase;
      t36 = sin(t35);
      t40 = sin(omstry * y);
      t42 = 1 + t40 / 2;
      t43 = cos(t7);
      t45 = sin(t22);
      forces[1] =
          mua[ind] * (-t28 * t32 * t36 * omega + t42 * t43 * t45 * omega * t13);
      t51 = cos(t12);
      t56 = sin(t4);
      t61 = omega * t51;
      t66 = sin(t31);
      t68 = cos(t35);
      forces[2] =
          2 * mua[ind] * t43 * t23 * t51 * omega +
          lambdaa[ind] * (t19 * (t56 * omega * t8 * t51 + t5 * t43 * t61) +
                          t42 * t28 * t66 * omega * t68 + t43 * t23 * t61);

      bforce[3 * qq] -= forces[0];
      bforce[1 + 3 * qq] -= forces[1];
      bforce[2 + 3 * qq] -= forces[2];
    }
}

//-----------------------------------------------------------------------
void EW::twfrsurfz_wind_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4 h, int kz,
                           float_sw4 t, float_sw4 omega, float_sw4 c,
                           float_sw4 phase, float_sw4* __restrict__ bforce,
                           float_sw4* __restrict__ mu,
                           float_sw4* __restrict__ lambda, float_sw4 zmin,
                           int i1, int i2, int j1, int j2) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  float_sw4 z = (kz - 1) * h + zmin;
// the do loops should span j1,j2 and i1,i2
#pragma omp parallel for
  for (int j = j1; j <= j2; j++) {
    float_sw4 t13, t15, t16, t19, t20, t21, t23, t24, t28, t29, t32, t33, t34,
        t37, t38, t43, t44, t49, t60, t62, t65;
    float_sw4 y = (j - 1) * h;
#pragma ivdep
#pragma simd
    for (int i = i1; i <= i2; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      float_sw4 x = (i - 1) * h;
      float_sw4 forces[3];
      t13 = mu[ind];
      t15 = omega * x + phase;
      t16 = cos(t15);
      t19 = omega * y + phase;
      t20 = sin(t19);
      t21 = c * t;
      t23 = omega * (z - t21);
      t24 = sin(t23);
      t28 = omega * (x - t21);
      t29 = sin(t28);
      t32 = omega * z + phase;
      t33 = cos(t32);
      t34 = t33 * omega;
      forces[0] = t13 * (t16 * omega * t20 * t24 + t29 * t20 * t34);
      t37 = sin(t15);
      t38 = cos(t19);
      t43 = omega * (y - t21);
      t44 = sin(t43);
      forces[1] = t13 * (t37 * t38 * omega * t24 + t37 * t44 * t34);
      t49 = cos(t23);
      t60 = cos(t28);
      t62 = sin(t32);
      t65 = cos(t43);
      forces[2] =
          2 * t13 * t37 * t20 * t49 * omega +
          lambda[ind] * (t60 * omega * t20 * t62 + t37 * t65 * omega * t62 +
                         t37 * t20 * t49 * omega);
      bforce[qq] = forces[0];
      bforce[qq + nij] = forces[1];
      bforce[qq + 2 * nij] = forces[2];
    }
  }
}

//-----------------------------------------------------------------------
void EW::twfrsurfzsg_wind_ci(int ifirst, int ilast, int jfirst, int jlast,
                             int kfirst, int klast, float_sw4 h, int kz,
                             float_sw4 t, float_sw4 om, float_sw4 c,
                             float_sw4 ph, float_sw4 omstrx, float_sw4 omstry,
                             float_sw4* __restrict__ bforce,
                             float_sw4* __restrict__ mu,
                             float_sw4* __restrict__ lambda, float_sw4 zmin,
                             int i1, int i2, int j1, int j2) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  float_sw4 z = (kz - 1) * h + zmin;
#pragma omp parallel for
  for (int j = j1; j <= j2; j++) {
    float_sw4 t1, t10, t11, t12, t15, t17, t19, t20, t22, t24, t25, t29, t3,
        t31, t32, t36, t39, t4, t40, t46, t51, t53, t56, t6, t7;
    float_sw4 y = (j - 1) * h;
#pragma ivdep
#pragma simd
    for (int i = i1; i <= i2; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      float_sw4 x = (i - 1) * h;
      float_sw4 forces[3];
      t1 = c * t;
      t3 = om * (x - t1);
      t4 = sin(t3);
      t6 = om * y + ph;
      t7 = sin(t6);
      t10 = om * z + ph;
      t11 = cos(t10);
      t12 = t11 * om;
      t15 = sin(omstrx * x);
      t17 = 1 + t15 / 2;
      t19 = om * x + ph;
      t20 = cos(t19);
      t22 = om * t7;
      t24 = om * (z - t1);
      t25 = sin(t24);
      forces[0] = mu[ind] * (t4 * t7 * t12 + t17 * t20 * t22 * t25);
      t29 = sin(t19);
      t31 = om * (y - t1);
      t32 = sin(t31);
      t36 = sin(omstry * y);
      t39 = (1 + t36 / 2) * t29;
      t40 = cos(t6);
      forces[1] = mu[ind] * (t29 * t32 * t12 + t39 * t40 * om * t25);
      t46 = cos(t24);
      t51 = cos(t3);
      t53 = sin(t10);
      t56 = cos(t31);
      forces[2] = 2 * mu[ind] * t29 * t7 * t46 * om +
                  lambda[ind] * (t17 * t51 * t22 * t53 + t39 * t56 * om * t53 +
                                 t29 * t7 * t46 * om);
      bforce[qq] = forces[0];
      bforce[qq + nij] = forces[1];
      bforce[qq + 2 * nij] = forces[2];
    }
  }
}

//-----------------------------------------------------------------------
void EW::twfrsurfz_att_wind_ci(int ifirst, int ilast, int jfirst, int jlast,
                               int kfirst, int klast, float_sw4 h, int kz,
                               float_sw4 t, float_sw4 omega, float_sw4 c,
                               float_sw4 phase, float_sw4* __restrict__ bforce,
                               float_sw4* __restrict__ mua,
                               float_sw4* __restrict__ lambdaa, float_sw4 zmin,
                               int i1, int i2, int j1, int j2) {
  SW4_MARK_FUNCTION;
  // THIS ROUTINE ACCUMULATES CONTRIBUTIONS TO 'bforce'
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  float_sw4 z = (kz - 1) * h + zmin;
#pragma omp parallel for
  for (int j = j1; j <= j2; j++) {
    float_sw4 t2, t3, t6, t7, t8, t11, t12, t17, t18, t27, t30, t20, t31, t35;
    float_sw4 t40, t45, t50, t52, t54, t16, t23, t24, t34;
    float_sw4 y = (j - 1) * h;
#pragma ivdep
#pragma simd
    for (int i = i1; i <= i2; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      float_sw4 x = (i - 1) * h;
      float_sw4 forces[3];
      t2 = omega * x + phase;
      t3 = sin(t2);
      t6 = omega * y + phase;
      t7 = cos(t6);
      t8 = c * t;
      t11 = -omega * (z - t8) - phase;
      t12 = sin(t11);
      t16 = omega * (x - t8);
      t17 = -t16 - phase;
      t18 = cos(t17);
      t20 = t12 * omega;
      forces[0] = mua[ind] * (t3 * omega * t7 * t12 + t18 * t3 * t20);
      t23 = cos(t2);
      t24 = sin(t6);
      t27 = sin(t16);
      t30 = -omega * (y - t8) - phase;
      t31 = cos(t30);
      t34 = omega * z + phase;
      t35 = sin(t34);
      forces[1] = mua[ind] * (t23 * t24 * t20 - t27 * t31 * t35 * omega);
      t40 = cos(t11);
      t45 = sin(t17);
      t50 = t40 * omega;
      t52 = sin(t30);
      t54 = cos(t34);
      forces[2] = 2.0 * mua[ind] * t23 * t7 * t40 * omega +
                  lambdaa[ind] * (t45 * omega * t3 * t40 + t18 * t23 * t50 +
                                  t27 * t52 * omega * t54 + t23 * t7 * t50);
      bforce[qq] -= forces[0];
      bforce[qq + nij] -= forces[1];
      bforce[qq + 2 * nij] -= forces[2];
    }
  }
}

//-----------------------------------------------------------------------
void EW::twfrsurfzsg_att_wind_ci(
    int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
    float_sw4 h, int kz, float_sw4 t, float_sw4 omega, float_sw4 c,
    float_sw4 phase, float_sw4 omstrx, float_sw4 omstry,
    float_sw4* __restrict__ bforce, float_sw4* __restrict__ mua,
    float_sw4* __restrict__ lambdaa, float_sw4 zmin, int i1, int i2, int j1,
    int j2) {
  SW4_MARK_FUNCTION;
  // THIS ROUTINE ACCUMULATES CONTRIBUTIONS TO 'bforce'
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  float_sw4 z = (kz - 1) * h + zmin;
#pragma omp parallel for
  for (int j = j1; j <= j2; j++) {
    float_sw4 t1, t12, t13, t17, t19, t22, t23, t28, t3, t31, t32, t35, t36, t4,
        t40, t42, t43, t45, t5, t51, t56, t61, t66, t68, t7, t8;
    float_sw4 y = (j - 1) * h;
#pragma ivdep
#pragma simd
    for (int i = i1; i <= i2; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      float_sw4 x = (i - 1) * h;
      float_sw4 forces[3];
      t1 = c * t;
      t3 = omega * (x - t1);
      t4 = -t3 - phase;
      t5 = cos(t4);
      t7 = omega * x + phase;
      t8 = sin(t7);
      t12 = -omega * (z - t1) - phase;
      t13 = sin(t12);
      t17 = sin(omstrx * x);
      t19 = 1 + t17 / 2;
      t22 = omega * y + phase;
      t23 = cos(t22);
      forces[0] =
          mua[ind] * (t5 * t8 * t13 * omega + t19 * t8 * omega * t23 * t13);
      t28 = sin(t3);
      t31 = -omega * (y - t1) - phase;
      t32 = cos(t31);
      t35 = omega * z + phase;
      t36 = sin(t35);
      t40 = sin(omstry * y);
      t42 = 1 + t40 / 2;
      t43 = cos(t7);
      t45 = sin(t22);
      forces[1] =
          mua[ind] * (-t28 * t32 * t36 * omega + t42 * t43 * t45 * omega * t13);
      t51 = cos(t12);
      t56 = sin(t4);
      t61 = omega * t51;
      t66 = sin(t31);
      t68 = cos(t35);
      forces[2] =
          2 * mua[ind] * t43 * t23 * t51 * omega +
          lambdaa[ind] * (t19 * (t56 * omega * t8 * t51 + t5 * t43 * t61) +
                          t42 * t28 * t66 * omega * t68 + t43 * t23 * t61);
      bforce[qq] -= forces[0];
      bforce[qq + nij] -= forces[1];
      bforce[qq + 2 * nij] -= forces[2];
    }
  }
}

//-----------------------------------------------------------------------
void EW::twstensor_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                      int klast, int kz, float_sw4 t, float_sw4 om, float_sw4 c,
                      float_sw4 ph, float_sw4* xx, float_sw4* yy, float_sw4* zz,
                      float_sw4* tau, float_sw4* mu, float_sw4* lambda) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++)
    for (int i = ifirst; i <= ilast; i++) {
      double t1, t11, t12, t20, t21, t23, t24, t26, t3, t30, t31, t35, t36, t37,
          t38, t4, t41, t42, t46, t50, t51, t54, t7, t8;
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      double x = xx[ind], y = yy[ind], z = zz[ind];
      double muu = mu[ind], lambdaa = lambda[ind];

      t1 = c * t;
      t3 = om * (x - t1);
      t4 = cos(t3);
      t7 = om * y + ph;
      t8 = sin(t7);
      t11 = om * z + ph;
      t12 = sin(t11);
      t20 = om * x + ph;
      t21 = sin(t20);
      t23 = om * (y - t1);
      t24 = cos(t23);
      t26 = om * t12;
      t30 = om * (z - t1);
      t31 = cos(t30);
      t35 = lambdaa *
            (t4 * om * t8 * t12 + t21 * t24 * t26 + t21 * t8 * t31 * om);
      tau[qq] = 2 * muu * t4 * om * t8 * t12 + t35;
      t36 = cos(t20);
      t37 = t36 * om;
      t38 = sin(t23);
      t41 = sin(t3);
      t42 = cos(t7);
      tau[qq + nij] = muu * (t37 * t38 * t12 + t41 * t42 * t26);
      t46 = sin(t30);
      t50 = cos(t11);
      t51 = t50 * om;
      tau[qq + 2 * nij] = muu * (t37 * t8 * t46 + t41 * t8 * t51);
      t54 = muu * t21;
      tau[qq + 3 * nij] = 2 * t54 * t24 * om * t12 + t35;
      tau[qq + 4 * nij] = muu * (t21 * t42 * om * t46 + t21 * t38 * t51);
      tau[qq + 5 * nij] = 2 * t54 * t8 * t31 * om + t35;
    }
}

//-----------------------------------------------------------------------
void EW::twstensorsg_ci(int ifirst, int ilast, int jfirst, int jlast,
                        int kfirst, int klast, int kz, float_sw4 t,
                        float_sw4 om, float_sw4 c, float_sw4 ph, float_sw4* xx,
                        float_sw4* yy, float_sw4* zz, float_sw4* tau,
                        float_sw4* mu, float_sw4* lambda, float_sw4 omstrx,
                        float_sw4 omstry) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++)
    for (int i = ifirst; i <= ilast; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      double x = xx[ind], y = yy[ind], z = zz[ind];
      double muu = mu[ind], lambdaa = lambda[ind];

      double t12, t13, t14, t16, t17, t18, t2, t24, t26, t28, t29, t30, t32,
          t33, t35, t39, t4, t40, t44, t45;
      double t46, t47, t51, t53, t54, t59, t6, t60, t62, t8, t9;
      t2 = sin(omstrx * x);
      t4 = 1 + t2 / 2;
      t6 = c * t;
      t8 = om * (x - t6);
      t9 = cos(t8);
      t12 = om * y + ph;
      t13 = sin(t12);
      t14 = om * t13;
      t16 = om * z + ph;
      t17 = sin(t16);
      t18 = t14 * t17;
      t24 = sin(omstry * y);
      t26 = 1 + t24 / 2;
      t28 = om * x + ph;
      t29 = sin(t28);
      t30 = t26 * t29;
      t32 = om * (y - t6);
      t33 = cos(t32);
      t35 = t33 * om * t17;
      t39 = om * (z - t6);
      t40 = cos(t39);
      t44 = lambdaa * (t4 * t9 * t18 + t30 * t35 + t29 * t13 * t40 * om);
      tau[qq] = 2 * muu * t4 * t9 * t18 + t44;
      t45 = cos(t28);
      t46 = t4 * t45;
      t47 = sin(t32);
      t51 = sin(t8);
      t53 = cos(t12);
      t54 = t53 * om;
      tau[qq + nij] = muu * (t46 * om * t47 * t17 + t26 * t51 * t54 * t17);
      t59 = cos(t16);
      t60 = t59 * om;
      t62 = sin(t39);
      tau[qq + 2 * nij] = muu * (t51 * t13 * t60 + t46 * t14 * t62);
      tau[qq + 3 * nij] = 2 * muu * t26 * t29 * t35 + t44;
      tau[qq + 4 * nij] = muu * (t29 * t47 * t60 + t30 * t54 * t62);
      tau[qq + 5 * nij] = 2 * muu * t29 * t13 * t40 * om + t44;
    }
}

//-----------------------------------------------------------------------
void EW::twstensoratt_ci(int ifirst, int ilast, int jfirst, int jlast,
                         int kfirst, int klast, int kz, float_sw4 t,
                         float_sw4 omega, float_sw4 c, float_sw4 phase,
                         float_sw4* xx, float_sw4* yy, float_sw4* zz,
                         float_sw4* tau, float_sw4* mu, float_sw4* lambda) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++)
    for (int i = ifirst; i <= ilast; i++) {
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      double x = xx[ind], y = yy[ind], z = zz[ind];
      double muu = mu[ind], lambdaa = lambda[ind];
      double t1, t12, t13, t15, t16, t17, t19, t20, t24, t27, t28, t3, t31, t32,
          t36, t37, t4, t41, t42, t44, t48, t49, t5, t61, t64, t8, t9;
      double stensor[6];

      t1 = c * t;
      t3 = omega * (x - t1);
      t4 = -t3 - phase;
      t5 = sin(t4);
      t8 = omega * x + phase;
      t9 = sin(t8);
      t12 = -omega * (z - t1) - phase;
      t13 = cos(t12);
      t15 = t5 * omega * t9 * t13;
      t16 = cos(t4);
      t17 = cos(t8);
      t19 = omega * t13;
      t20 = t16 * t17 * t19;
      t24 = sin(t3);
      t27 = -omega * (y - t1) - phase;
      t28 = sin(t27);
      t31 = omega * z + phase;
      t32 = cos(t31);
      t36 = omega * y + phase;
      t37 = cos(t36);
      t41 = lambdaa * (t15 + t20 + t24 * t28 * omega * t32 + t17 * t37 * t19);
      stensor[0] = 2 * muu * (t15 + t20) + t41;
      t42 = cos(t3);
      t44 = cos(t27);
      stensor[1] = muu * t42 * omega * t44 * t32;
      t48 = sin(t12);
      t49 = t48 * omega;
      stensor[2] = muu * (t16 * t9 * t49 + t9 * omega * t37 * t48);
      stensor[3] = 2 * muu * t24 * t28 * omega * t32 + t41;
      t61 = sin(t31);
      t64 = sin(t36);
      stensor[4] = muu * (-t24 * t44 * t61 * omega + t17 * t64 * t49);
      stensor[5] = 2 * muu * t17 * t37 * t13 * omega + t41;

      tau[qq] = stensor[0];
      tau[qq + nij] = stensor[1];
      tau[qq + 2 * nij] = stensor[2];
      tau[qq + 3 * nij] = stensor[3];
      tau[qq + 4 * nij] = stensor[4];
      tau[qq + 5 * nij] = stensor[5];
    }
}

//-----------------------------------------------------------------------
void EW::twstensorsgatt_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, int kz, float_sw4 t,
                           float_sw4 omega, float_sw4 c, float_sw4 phase,
                           float_sw4* xx, float_sw4* yy, float_sw4* zz,
                           float_sw4* tau, float_sw4* mu, float_sw4* lambda,
                           float_sw4 omstrx, float_sw4 omstry) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
#pragma omp parallel for
  for (int j = jfirst; j <= jlast; j++)
    for (int i = ifirst; i <= ilast; i++) {
      double t10, t13, t14, t17, t18, t2, t21, t22, t24, t26, t31, t33, t34,
          t38, t39, t4, t42, t43, t44;
      double t47, t48, t5, t52, t53, t55, t59, t6, t72, t76, t8, t9;
      double stensor[6];
      size_t ind = (i - ifirst) + ni * (j - jfirst) + nij * (kz - kfirst);
      size_t qq = (i - ifirst) + ni * (j - jfirst);
      double x = xx[ind], y = yy[ind], z = zz[ind];
      double muu = mu[ind], lambdaa = lambda[ind];

      t2 = sin(omstrx * x);
      t4 = 1 + t2 / 2;
      t5 = muu * t4;
      t6 = c * t;
      t8 = omega * (x - t6);
      t9 = -t8 - phase;
      t10 = sin(t9);
      t13 = omega * x + phase;
      t14 = sin(t13);
      t17 = -omega * (z - t6) - phase;
      t18 = cos(t17);
      t21 = cos(t9);
      t22 = cos(t13);
      t24 = omega * t18;
      t26 = t10 * omega * t14 * t18 + t21 * t22 * t24;
      t31 = sin(omstry * y);
      t33 = 1 + t31 / 2;
      t34 = sin(t8);
      t38 = -omega * (y - t6) - phase;
      t39 = sin(t38);
      t42 = omega * z + phase;
      t43 = cos(t42);
      t44 = t39 * omega * t43;
      t47 = omega * y + phase;
      t48 = cos(t47);
      t52 = lambdaa * (t4 * t26 + t33 * t34 * t44 + t22 * t48 * t24);
      stensor[0] = 2 * t5 * t26 + t52;
      t53 = cos(t8);
      t55 = cos(t38);
      stensor[1] = t5 * t53 * omega * t55 * t43;
      t59 = sin(t17);
      stensor[2] =
          muu * (t21 * t14 * t59 * omega + t4 * t14 * omega * t48 * t59);
      stensor[3] = 2 * muu * t33 * t34 * t44 + t52;
      t72 = sin(t42);
      t76 = sin(t47);
      stensor[4] =
          muu * (-t34 * t55 * t72 * omega + t33 * t22 * t76 * omega * t59);
      stensor[5] = 2 * muu * t22 * t48 * t18 * omega + t52;

      tau[qq] = stensor[0];
      tau[qq + nij] = stensor[1];
      tau[qq + 2 * nij] = stensor[2];
      tau[qq + 3 * nij] = stensor[3];
      tau[qq + 4 * nij] = stensor[4];
      tau[qq + 5 * nij] = stensor[5];
    }
}
