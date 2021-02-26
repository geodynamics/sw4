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
#include "foralls.h"
#include "mpi.h"
#include "policies.h"
//-----------------------------------------------------------------------
bool EW::proc_decompose_2d(int ni, int nj, int nproc, int proc_max[2]) {
  // This routine determines a decomposition of nproc processors into
  // a 2D processor array  proc_max[0] x proc_max[1], which gives minimal
  // communication boundary for a grid with ni x nj points.

  float_sw4 fmin = ni + nj;
  bool first = true;
  int p1max = ni / m_ppadding;
  int p2max = nj / m_ppadding;
  for (int p1 = 1; p1 <= nproc; p1++)
    if (nproc % p1 == 0) {
      int p2 = nproc / p1;
      if (p1 <= p1max && p2 <= p2max) {
        // int w1 = p1==1?0:1;
        // int w2 = p2==1?0:1;
        // double f = w2*(double)(ni)/p1 + w1*(double)(nj)/p2;
        // try to make each subdomain as square as possible
        float_sw4 f = fabs((float_sw4)(ni) / p1 - (float_sw4)(nj) / p2);
        if (f < fmin || first) {
          fmin = f;
          proc_max[0] = p1;
          proc_max[1] = p2;
          first = false;
        }
      }
    }
  return !first;
}

//-----------------------------------------------------------------------
void EW::coarsen1d(int& n, int& ifirst, int& ilast, int periodic) {
  // n - total number of points 1<=i<=n,
  // Total index space is 1-ghosts <= i <= n + ghosts
  //
  // This routine coarsens the interval ifirst <= i <= ilast
  // by a factor two, and returns coarsened values of
  // all input parameters.
  //
  int nc;
  if (periodic)
    nc = n / 2;
  else
    nc = (n - 1) / 2 + 1;

  if (ilast == n + m_ghost_points)
    ilast = nc + m_ghost_points;
  else {
    ilast = ilast - m_ppadding;
    if (ilast % 2 == 0) ilast--;
    ilast = (ilast + 1) / 2;
    ilast = ilast + m_ppadding;
  }
  if (ifirst != 1 - m_ghost_points) {
    ifirst = ifirst + m_ppadding;
    if (ifirst % 2 == 0) ifirst++;
    ifirst = (ifirst + 1) / 2;
    ifirst = ifirst - m_ppadding;
  }
  n = nc;
}

//-----------------------------------------------------------------------
void EW::decomp1d(int nglobal, int myid, int nproc, int& s, int& e)
//
// Decompose index space 1 <= i <= nglobal into nproc blocks
// returns start and end indices for block nr. myid,
//          where 0 <= myid <= nproc-1
//
{
  int olap = 2 * m_ppadding;
  int nlocal = (nglobal + (nproc - 1) * olap) / nproc;
  int deficit = (nglobal + (nproc - 1) * olap) % nproc;

  if (myid < deficit)
    s = myid * (nlocal - olap) + myid + 1;
  else
    s = myid * (nlocal - olap) + deficit + 1;

  if (myid < deficit) nlocal = nlocal + 1;

  e = s + nlocal - 1;
}

// -----------------------------
void EW::setup2D_MPICommunications() {
  SW4_MARK_FUNCTION;
  if (mVerbose >= 2 && proc_zero())
    cout << "***inside setup2D_MPICommunications***" << endl;
  // // Define MPI datatypes for communication across processor boundaries
  // // For topography: finest grid (curvilinear) only, only one value per grid
  // // point (nc=1) get the size from the top Cartesian grid
  int g = mNumberOfGrids - 1;
  int ni = m_iEnd[g] - m_iStart[g] + 1, nj = m_jEnd[g] - m_jStart[g] + 1;
  // std::cout<<"GRID = "<<g<<" size "<<ni<<" "<<nj<<"\n";
  // MPI_Type_vector(nj, m_ppadding, ni, m_mpifloat, &m_send_type_2dfinest[0]);
  // MPI_Type_vector(1, m_ppadding * ni, ni * nj, m_mpifloat,
  //                 &m_send_type_2dfinest[1]);
  // MPI_Type_commit(&m_send_type_2dfinest[0]);
  // MPI_Type_commit(&m_send_type_2dfinest[1]);

  // // Extended number of padding points
  ni = ni + 2 * m_ext_ghost_points;
  nj = nj + 2 * m_ext_ghost_points;
  int extpadding = m_ppadding + m_ext_ghost_points;
  // MPI_Type_vector(nj, extpadding, ni, m_mpifloat,
  // &m_send_type_2dfinest_ext[0]); MPI_Type_vector(1, extpadding * ni, ni * nj,
  // m_mpifloat,
  //                 &m_send_type_2dfinest_ext[1]);

  // MPI_Type_commit(&m_send_type_2dfinest_ext[0]);
  // MPI_Type_commit(&m_send_type_2dfinest_ext[1]);

  send_type_2dfinest_ext.resize(2);
  bufs_type_2dfinest_ext.resize(2);

  make_type_2d(send_type_2dfinest_ext, bufs_type_2dfinest_ext, nj, extpadding,
               ni, 0);
  make_type_2d(send_type_2dfinest_ext, bufs_type_2dfinest_ext, 1,
               extpadding * ni, ni * nj, 1);

  // // For mesh refinement: 2D planes with three values per grid point (nc=3)
  // // Coarser grids
  // m_send_type_2dx.resize(mNumberOfCartesianGrids);
  // m_send_type_2dy.resize(mNumberOfCartesianGrids);
  // m_send_type_2dx1p.resize(mNumberOfCartesianGrids);  // padding=1
  // m_send_type_2dy1p.resize(mNumberOfCartesianGrids);
  // m_send_type_2dx3p.resize(mNumberOfCartesianGrids);  // padding=3
  // m_send_type_2dy3p.resize(mNumberOfCartesianGrids);

  // // Data for ASYNC_SEND_RECV
  send_type_2dx.resize(mNumberOfCartesianGrids);
  send_type_2dy.resize(mNumberOfCartesianGrids);

  bufs_type_2dx.resize(mNumberOfCartesianGrids);
  bufs_type_2dy.resize(mNumberOfCartesianGrids);

  for (int g = 0; g < mNumberOfCartesianGrids; g++) {
    ni = m_iEnd[g] - m_iStart[g] + 1;
    nj = m_jEnd[g] - m_jStart[g] + 1;
    //   if (m_croutines) {
    //     MPI_Type_vector(3 * nj, m_ppadding, ni, m_mpifloat,
    //     &m_send_type_2dx[g]); MPI_Type_vector(3, m_ppadding * ni, ni * nj,
    //     m_mpifloat,
    //                     &m_send_type_2dy[g]);

    make_type_2d(send_type_2dx, bufs_type_2dx, 3 * nj, m_ppadding, ni, g);
    make_type_2d(send_type_2dy, bufs_type_2dy, 3, m_ppadding * ni, ni * nj, g);

    //     MPI_Type_vector(3 * nj, 1, ni, m_mpifloat, &m_send_type_2dx1p[g]);
    //     MPI_Type_vector(3, ni, ni * nj, m_mpifloat, &m_send_type_2dy1p[g]);
    //     MPI_Type_vector(3 * nj, 3, ni, m_mpifloat, &m_send_type_2dx3p[g]);
    //     MPI_Type_vector(3, 3 * ni, ni * nj, m_mpifloat,
    //     &m_send_type_2dy3p[g]);
    //   } else {
    //     MPI_Type_vector(nj, 3 * m_ppadding, 3 * ni, m_mpifloat,
    //                     &m_send_type_2dx[g]);
    //     MPI_Type_vector(1, 3 * m_ppadding * ni, 3 * ni * nj, m_mpifloat,
    //                     &m_send_type_2dy[g]);
    //     MPI_Type_vector(nj, 3, 3 * ni, m_mpifloat, &m_send_type_2dx1p[g]);
    //     MPI_Type_vector(1, 3 * ni, 3 * ni * nj, m_mpifloat,
    //                     &m_send_type_2dy1p[g]);
    //     MPI_Type_vector(nj, 3 * 3, 3 * ni, m_mpifloat,
    //     &m_send_type_2dx3p[g]); MPI_Type_vector(1, 3 * 3 * ni, 3 * ni * nj,
    //     m_mpifloat,
    //                     &m_send_type_2dy3p[g]);
    //   }
    //   MPI_Type_commit(&m_send_type_2dx[g]);
    //   MPI_Type_commit(&m_send_type_2dy[g]);
    //   MPI_Type_commit(&m_send_type_2dx1p[g]);
    //   MPI_Type_commit(&m_send_type_2dy1p[g]);
    //   MPI_Type_commit(&m_send_type_2dx3p[g]);
    //   MPI_Type_commit(&m_send_type_2dy3p[g]);
  }

  // For topography: finest grid (curvilinear) only, only one value per grid
  // point (nc=1) get the size from the top curvlinear grid
  g = mNumberOfGrids - 1;
  ni = m_iEnd[g] - m_iStart[g] + 1;
  nj = m_jEnd[g] - m_jStart[g] + 1;
  MPI_Type_vector(nj, m_ppadding, ni, m_mpifloat, &m_send_type_2dfinest[0]);
  MPI_Type_vector(1, m_ppadding * ni, ni * nj, m_mpifloat,
                  &m_send_type_2dfinest[1]);
  MPI_Type_commit(&m_send_type_2dfinest[0]);
  MPI_Type_commit(&m_send_type_2dfinest[1]);

  // Extended number of padding points
  ni = ni + 2 * m_ext_ghost_points;
  nj = nj + 2 * m_ext_ghost_points;
  extpadding = m_ppadding + m_ext_ghost_points;
  MPI_Type_vector(nj, extpadding, ni, m_mpifloat, &m_send_type_2dfinest_ext[0]);
  MPI_Type_vector(1, extpadding * ni, ni * nj, m_mpifloat,
                  &m_send_type_2dfinest_ext[1]);

  MPI_Type_commit(&m_send_type_2dfinest_ext[0]);
  MPI_Type_commit(&m_send_type_2dfinest_ext[1]);

  // NEW: July-2019 communicators for interface surfaces
  int numSurfaces = mNumberOfGrids - mNumberOfCartesianGrids;
  m_send_type_isurfx.resize(numSurfaces);
  m_send_type_isurfy.resize(numSurfaces);

  for (int iSurf = 0; iSurf < numSurfaces; iSurf++) {
    int g = mNumberOfCartesianGrids + iSurf;
    int ni = m_iEnd[g] - m_iStart[g] + 1 + 2 * m_ext_ghost_points;
    int nj = m_jEnd[g] - m_jStart[g] + 1 + 2 * m_ext_ghost_points;

    MPI_Type_vector(nj, extpadding, ni, m_mpifloat, &m_send_type_isurfx[iSurf]);
    MPI_Type_vector(1, extpadding * ni, ni * nj, m_mpifloat,
                    &m_send_type_isurfy[iSurf]);

    MPI_Type_commit(&m_send_type_isurfx[iSurf]);
    MPI_Type_commit(&m_send_type_isurfy[iSurf]);
  }

  // For mesh refinement: 2D planes with three values per grid point (nc=3)
  // Coarser grids
  m_send_type_2dx.resize(mNumberOfGrids);
  m_send_type_2dy.resize(mNumberOfGrids);
  m_send_type_2dx1p.resize(mNumberOfGrids);  // padding=1
  m_send_type_2dy1p.resize(mNumberOfGrids);
  m_send_type_2dx3p.resize(mNumberOfGrids);  // padding=3
  m_send_type_2dy3p.resize(mNumberOfGrids);
  for (int g = 0; g < mNumberOfGrids; g++) {
    int ni = m_iEnd[g] - m_iStart[g] + 1, nj = m_jEnd[g] - m_jStart[g] + 1;

    MPI_Type_vector(3 * nj, m_ppadding, ni, m_mpifloat, &m_send_type_2dx[g]);
    MPI_Type_vector(3, m_ppadding * ni, ni * nj, m_mpifloat,
                    &m_send_type_2dy[g]);
    MPI_Type_vector(3 * nj, 1, ni, m_mpifloat, &m_send_type_2dx1p[g]);
    MPI_Type_vector(3, ni, ni * nj, m_mpifloat, &m_send_type_2dy1p[g]);
    MPI_Type_vector(3 * nj, 3, ni, m_mpifloat, &m_send_type_2dx3p[g]);
    MPI_Type_vector(3, 3 * ni, ni * nj, m_mpifloat, &m_send_type_2dy3p[g]);

    MPI_Type_commit(&m_send_type_2dx[g]);
    MPI_Type_commit(&m_send_type_2dy[g]);
    MPI_Type_commit(&m_send_type_2dx1p[g]);
    MPI_Type_commit(&m_send_type_2dy1p[g]);
    MPI_Type_commit(&m_send_type_2dx3p[g]);
    MPI_Type_commit(&m_send_type_2dy3p[g]);
  }

#ifdef SW4_STAGED_MPI_BUFFERS

#ifdef SW4_USE_UMPIRE

  global_variables.device_buffer =
      SW4_NEW(Space::Managed_temps, float_sw4[global_variables.buffer_size]);
#else
    void* ptr;
  if (cudaMalloc(&ptr, global_variables.buffer_size * 8) != cudaSuccess) {
    std::cerr << "cudaMalloc failed in line 387 of parallelStuff.C\n";
    abort();
  } else {
    std::cout << "Device buffer of size " << global_variables.buffer_size * 8
              << " bytes allocated in 2D\n";
    global_variables.device_buffer = (float_sw4*)ptr;
  }
#endif
#endif
}

// -----------------------------
void EW::setupMPICommunications() {
  if (mVerbose >= 2 && proc_zero())
    cout << "***inside setupMPICommunications***" << endl;
  // Define MPI datatypes for communication across processor boundaries
  m_send_type1.resize(2 * mNumberOfGrids);
  m_send_type3.resize(2 * mNumberOfGrids);
  m_send_type4.resize(2 * mNumberOfGrids);
  m_send_type21.resize(2 * mNumberOfGrids);

  // Data for ASYNC_SEND_RECV

  send_type1.resize(2 * mNumberOfGrids);
  send_type3.resize(2 * mNumberOfGrids);
  send_type4.resize(2 * mNumberOfGrids);
  send_type21.resize(2 * mNumberOfGrids);

  bufs_type1.resize(4 * mNumberOfGrids);
  bufs_type3.resize(4 * mNumberOfGrids);
  bufs_type4.resize(4 * mNumberOfGrids);
  bufs_type21.resize(4 * mNumberOfGrids);

  // End Data for ASYNC_SEND_RECV
  for (int g = 0; g < mNumberOfGrids; g++) {
    //      int ni = mU[g].m_ni, nj=mU[g].m_nj, nk=mU[g].m_nk;
    int ni = m_iEnd[g] - m_iStart[g] + 1;
    int nj = m_jEnd[g] - m_jStart[g] + 1;
    int nk = m_kEnd[g] - m_kStart[g] + 1;

    MPI_Type_vector(nj * nk, m_ppadding, ni, m_mpifloat, &m_send_type1[2 * g]);
    MPI_Type_vector(nk, m_ppadding * ni, ni * nj, m_mpifloat,
                    &m_send_type1[2 * g + 1]);

    // ASR START

    // send_type1[2*g]=std::make_tuple(nj*nk, m_ppadding, ni);
    // send_type1[2*g+1]=std::make_tuple(nk, m_ppadding*ni, ni*nj);

    // float_sw4* tbuf =
    // SW4_NEW(Managed,float_sw4[nj*nk*m_ppadding*4+nk*m_ppadding*ni*4]);
    // bufs_type1[4*g+0]=std::make_tuple(tbuf,tbuf+nj*nk*m_ppadding);
    // bufs_type1[4*g+1]=std::make_tuple(tbuf+2*nj*nk*m_ppadding,tbuf+3*nj*nk*m_ppadding);
    // tbuf+=4*nj*nk*m_ppadding;
    // bufs_type1[4*g+2]=std::make_tuple

    make_type(send_type1, bufs_type1, nj * nk, m_ppadding, ni, nk,
              m_ppadding * ni, ni * nj, g);
    // ASR END

    if (m_croutines) {
      MPI_Type_vector(3 * nj * nk, m_ppadding, ni, m_mpifloat,
                      &m_send_type3[2 * g]);
      MPI_Type_vector(3 * nk, m_ppadding * ni, ni * nj, m_mpifloat,
                      &m_send_type3[2 * g + 1]);

      make_type(send_type3, bufs_type3, 3 * nj * nk, m_ppadding, ni, 3 * nk,
                m_ppadding * ni, ni * nj, g);

      MPI_Type_vector(4 * nj * nk, m_ppadding, ni, m_mpifloat,
                      &m_send_type4[2 * g]);
      MPI_Type_vector(4 * nk, m_ppadding * ni, ni * nj, m_mpifloat,
                      &m_send_type4[2 * g + 1]);

      make_type(send_type4, bufs_type4, 4 * nj * nk, m_ppadding, ni, 4 * nk,
                m_ppadding * ni, ni * nj, g);

      MPI_Type_vector(21 * nj * nk, m_ppadding, ni, m_mpifloat,
                      &m_send_type21[2 * g]);
      MPI_Type_vector(21 * nk, m_ppadding * ni, ni * nj, m_mpifloat,
                      &m_send_type21[2 * g + 1]);

      make_type(send_type21, bufs_type21, 21 * nj * nk, m_ppadding, ni, 21 * nk,
                m_ppadding * ni, ni * nj, g);
    } else {
      MPI_Type_vector(nj * nk, 3 * m_ppadding, 3 * ni, m_mpifloat,
                      &m_send_type3[2 * g]);
      MPI_Type_vector(nk, 3 * m_ppadding * ni, 3 * ni * nj, m_mpifloat,
                      &m_send_type3[2 * g + 1]);

      MPI_Type_vector(nj * nk, 4 * m_ppadding, 4 * ni, m_mpifloat,
                      &m_send_type4[2 * g]);
      MPI_Type_vector(nk, 4 * m_ppadding * ni, 4 * ni * nj, m_mpifloat,
                      &m_send_type4[2 * g + 1]);

      MPI_Type_vector(nj * nk, 21 * m_ppadding, 21 * ni, m_mpifloat,
                      &m_send_type21[2 * g]);
      MPI_Type_vector(nk, 21 * m_ppadding * ni, 21 * ni * nj, m_mpifloat,
                      &m_send_type21[2 * g + 1]);
    }
    MPI_Type_commit(&m_send_type1[2 * g]);
    MPI_Type_commit(&m_send_type1[2 * g + 1]);

    MPI_Type_commit(&m_send_type3[2 * g]);
    MPI_Type_commit(&m_send_type3[2 * g + 1]);

    MPI_Type_commit(&m_send_type4[2 * g]);
    MPI_Type_commit(&m_send_type4[2 * g + 1]);

    MPI_Type_commit(&m_send_type21[2 * g]);
    MPI_Type_commit(&m_send_type21[2 * g + 1]);
  }

  // test call
  //   communicate_array( mRho[0], 0 );

  //   int g=mNumberOfGrids-1;
  //   mRho[g].set_to_zero();
  //   for( int j=mU[g].m_jb+m_ppadding ; j <= mU[g].m_je-m_ppadding ; j++ )
  //      for( int i=mU[g].m_ib+m_ppadding ; i <= mU[g].m_ie-m_ppadding ; i++ )
  //      {
  //	 mRho[g](i,j,1) = j;
  //      }
  //   communicate_array_2dfinest(mRho[g]);

  //   int myid;
  //   int seerank;
  //   MPI_Comm_rank( m_cartesian_communicator, &myid );
  //   do{
  //      if( myid == 0 )
  //      {
  //	 cout << "Give rank no (or -1 to exit) > " << endl;
  //	 cin >> seerank ;
  //      }
  //      MPI_Bcast( &seerank, 1, MPI_INT, 0, m_cartesian_communicator );
  //      if( myid == seerank )
  //      {
  //         cout << "Array bounds: " << mU[g].m_ib << " <= i <= " << mU[g].m_ie
  //         << " , "
  //	      << mU[g].m_jb << " <= j <= " << mU[g].m_je << endl;
  //	 for( int j=mU[g].m_je ; j >= mU[g].m_jb ; j-- )
  //	 {
  //	    for( int i=mU[g].m_ib ; i <= mU[g].m_ie ; i++ )
  //	       cout << " " << mRho[g](i,j,1);
  //	    cout << endl;
  //	 }
  //      }
  //   }
  //   while( seerank != -1 );

  //   MPI_Barrier(m_cartesian_communicator);
  //   REQUIRE2( 0==1,"dbg stop");

  //   if( dbg )
  //   {
  //      int tag = 3, dum;
  //      int myid;
  //      MPI_Comm_rank( m_cartesian_communicator, &myid );
  //      MPI_Status status;
  //      stringstream str;
  //      str << "dbg." << myid << ".dat";
  //      ofstream fileout(str.str().c_str());
  //      if( myid>0 )
  //	 MPI_Recv( &dum, 1, MPI_INT, myid-1, tag, m_cartesian_communicator,
  //&status );
  //      MPI_Comm_rank( m_cartesian_communicator, &myid );
  //      fileout << "Id in 1d " << myid << " Id in 2d " << my_proc_coords[0] <<
  //      " x " << my_proc_coords[1]
  //	      << " has Neighbors x-direction "  << neigh[0] << " and " <<
  // neigh[1] << " y-dir " << neigh[2]
  //	      << " and " << neigh[3] << endl;
  //      if( myid < nprocs-1 )
  //	 MPI_Send( &dum, 1, MPI_INT, myid+1, tag, m_cartesian_communicator );
  //      fileout.close();
  //   }

#ifdef SW4_STAGED_MPI_BUFFERS
#ifdef SW4_USE_UMPIRE
  if (global_variables.device_buffer != 0) {
    ::operator delete[](global_variables.device_buffer, Space::Managed_temps);
  }
  global_variables.device_buffer =
      SW4_NEW(Space::Managed_temps, float_sw4[global_variables.buffer_size]);
#else
  void* ptr;
  if (cudaMalloc(&ptr, global_variables.buffer_size * 8) != cudaSuccess) {
    std::cerr << "cudaMalloc failed in line 387 of parallelStuff.C\n";
    abort();
  } else {
    std::cout << "Device buffer of size " << global_variables.buffer_size * 8
              << " bytes allocated\n";
    global_variables.device_buffer = (float_sw4*)ptr;
  }
#endif
#endif
}

#if defined(ENABLE_GPU)
void EW::communicate_array(Sarray& u, int grid) {
  // The async version using either device or managed memory works
  // spectrum-mpi/2018.02.05 on Ray. And it is slower on the Hayward case:
  // baseline communicate_array ( without -gpu) : 25 minutes 20 secs
  // communicate_array_async with device buffers and -gpu 29 minutes 18 secs
  // communicate_array_async with UM buffers and -gpu 29 minutes 20 secs
  communicate_array_async(u, grid);
  return;
}
#else
//-----------------------------------------------------------------------
void EW::communicate_array(Sarray& u, int grid) {
  SW4_MARK_FUNCTION;
  // The async version using either device or managed memory works
  // spectrum-mpi/2018.02.05 on Ray. And it is slower on the Hayward case:
  // baseline communicate_array ( without -gpu) : 25 minutes 20 secs
  // communicate_array_async with device buffers and -gpu 29 minutes 18 secs
  // communicate_array_async with UM buffers and -gpu 29 minutes 20 secs
  // communicate_array_async(u,grid);
  // return;
  // REQUIRE2( 0 <= grid && grid < mU.size() ,
  // 	    " Error in communicate_array, grid = " << grid );

  REQUIRE2(u.m_nc == 21 || u.m_nc == 4 || u.m_nc == 3 || u.m_nc == 1,
           "Communicate array, only implemented for one-, three-, four-, and "
           "21-component arrays"
               << " nc = " << u.m_nc);
  int ie = u.m_ie, ib = u.m_ib, je = u.m_je, jb = u.m_jb, ke = u.m_ke,
      kb = u.m_kb;
  MPI_Status status;
  if (u.m_nc == 1) {
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    // X-direction communication
    MPI_Sendrecv(&u(ie - (2 * m_ppadding - 1), jb, kb), 1,
                 m_send_type1[2 * grid], m_neighbor[1], xtag1, &u(ib, jb, kb),
                 1, m_send_type1[2 * grid], m_neighbor[0], xtag1,
                 m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(ib + m_ppadding, jb, kb), 1, m_send_type1[2 * grid],
                 m_neighbor[0], xtag2, &u(ie - (m_ppadding - 1), jb, kb), 1,
                 m_send_type1[2 * grid], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);
    // Y-direction communication
    MPI_Sendrecv(&u(ib, je - (2 * m_ppadding - 1), kb), 1,
                 m_send_type1[2 * grid + 1], m_neighbor[3], ytag1,
                 &u(ib, jb, kb), 1, m_send_type1[2 * grid + 1], m_neighbor[2],
                 ytag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(ib, jb + m_ppadding, kb), 1, m_send_type1[2 * grid + 1],
                 m_neighbor[2], ytag2, &u(ib, je - (m_ppadding - 1), kb), 1,
                 m_send_type1[2 * grid + 1], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
  } else if (u.m_nc == 3) {
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    // X-direction communication
    MPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, kb), 1,
                 m_send_type3[2 * grid], m_neighbor[1], xtag1,
                 &u(1, ib, jb, kb), 1, m_send_type3[2 * grid], m_neighbor[0],
                 xtag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib + m_ppadding, jb, kb), 1, m_send_type3[2 * grid],
                 m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, kb), 1,
                 m_send_type3[2 * grid], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);
    // Y-direction communication
    MPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), kb), 1,
                 m_send_type3[2 * grid + 1], m_neighbor[3], ytag1,
                 &u(1, ib, jb, kb), 1, m_send_type3[2 * grid + 1],
                 m_neighbor[2], ytag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib, jb + m_ppadding, kb), 1, m_send_type3[2 * grid + 1],
                 m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), kb), 1,
                 m_send_type3[2 * grid + 1], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
  } else if (u.m_nc == 4) {
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    // X-direction communication
    MPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, kb), 1,
                 m_send_type4[2 * grid], m_neighbor[1], xtag1,
                 &u(1, ib, jb, kb), 1, m_send_type4[2 * grid], m_neighbor[0],
                 xtag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib + m_ppadding, jb, kb), 1, m_send_type4[2 * grid],
                 m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, kb), 1,
                 m_send_type4[2 * grid], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);
    // Y-direction communication
    MPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), kb), 1,
                 m_send_type4[2 * grid + 1], m_neighbor[3], ytag1,
                 &u(1, ib, jb, kb), 1, m_send_type4[2 * grid + 1],
                 m_neighbor[2], ytag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib, jb + m_ppadding, kb), 1, m_send_type4[2 * grid + 1],
                 m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), kb), 1,
                 m_send_type4[2 * grid + 1], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
  } else if (u.m_nc == 21) {
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    // X-direction communication
    MPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, kb), 1,
                 m_send_type21[2 * grid], m_neighbor[1], xtag1,
                 &u(1, ib, jb, kb), 1, m_send_type21[2 * grid], m_neighbor[0],
                 xtag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib + m_ppadding, jb, kb), 1, m_send_type21[2 * grid],
                 m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, kb), 1,
                 m_send_type21[2 * grid], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);
    // Y-direction communication
    MPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), kb), 1,
                 m_send_type21[2 * grid + 1], m_neighbor[3], ytag1,
                 &u(1, ib, jb, kb), 1, m_send_type21[2 * grid + 1],
                 m_neighbor[2], ytag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib, jb + m_ppadding, kb), 1, m_send_type21[2 * grid + 1],
                 m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), kb), 1,
                 m_send_type21[2 * grid + 1], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
  }
}
#endif
//-----------------------------------------------------------------------
void EW::communicate_arrays(vector<Sarray>& u) {
  SW4_MARK_FUNCTION;
  for (int g = 0; g < u.size(); g++) communicate_array(u[g], g);
}

#if defined(ENABLE_GPU)
void EW::communicate_array_2d(Sarray& u, int g, int k) {
  communicate_array_2d_async(u, g, k);
  return;
}
#else
//-----------------------------------------------------------------------
void EW::communicate_array_2d(Sarray& u, int g, int k) {
  SW4_MARK_FUNCTION;
  REQUIRE2(u.m_nc == 3,
           "Communicate array 2d, only implemented for three-component arrays");
  REQUIRE2(g < m_send_type_2dx.size(),
           "Communicate array 2d, only implemented for grid=0.."
               << m_send_type_2dx.size() - 1 << " but g= " << g);
  int ie = m_iEnd[g], ib = m_iStart[g];
  int je = m_jEnd[g], jb = m_jStart[g];

  MPI_Status status;
  int xtag1 = 345;
  int xtag2 = 346;
  int ytag1 = 347;
  int ytag2 = 348;
  // communicate_array_2d_async( u, g, k );
  // communicate_array_2d_async_memo( u, g, k ); // Memoized version to avoid
  // possible page fault return;

  if (m_croutines && u.m_ke - u.m_kb + 1 != 1) {
    Sarray u2d(3, u.m_ib, u.m_ie, u.m_jb, u.m_je, k, k, __FILE__, __LINE__);
    u2d.copy_kplane(u, k);
    SW4_MARK_BEGIN("comm_array_2d::MPI");
    // X-direction communication
    MPI_Sendrecv(&u2d(1, ie - (2 * m_ppadding - 1), jb, k), 1,
                 m_send_type_2dx[g], m_neighbor[1], xtag1, &u2d(1, ib, jb, k),
                 1, m_send_type_2dx[g], m_neighbor[0], xtag1,
                 m_cartesian_communicator, &status);
    MPI_Sendrecv(&u2d(1, ib + m_ppadding, jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[0], xtag2, &u2d(1, ie - (m_ppadding - 1), jb, k), 1,
                 m_send_type_2dx[g], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);
    // Y-direction communication
    MPI_Sendrecv(&u2d(1, ib, je - (2 * m_ppadding - 1), k), 1,
                 m_send_type_2dy[g], m_neighbor[3], ytag1, &u2d(1, ib, jb, k),
                 1, m_send_type_2dy[g], m_neighbor[2], ytag1,
                 m_cartesian_communicator, &status);
    MPI_Sendrecv(&u2d(1, ib, jb + m_ppadding, k), 1, m_send_type_2dy[g],
                 m_neighbor[2], ytag2, &u2d(1, ib, je - (m_ppadding - 1), k), 1,
                 m_send_type_2dy[g], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
    SW4_MARK_END("comm_array_2d::MPI");
    u.copy_kplane(u2d, k);
  } else {
    // X-direction communication
    MPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[1], xtag1, &u(1, ib, jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[0], xtag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib + m_ppadding, jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, k), 1,
                 m_send_type_2dx[g], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);

    // Y-direction communication
    MPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), k), 1, m_send_type_2dy[g],
                 m_neighbor[3], ytag1, &u(1, ib, jb, k), 1, m_send_type_2dy[g],
                 m_neighbor[2], ytag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib, jb + m_ppadding, k), 1, m_send_type_2dy[g],
                 m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), k), 1,
                 m_send_type_2dy[g], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
  }
}
#endif  // if defined(ENABLE_GPU)
#if defined(ENABLE_GPU)
void EW::communicate_array_2d_ext(Sarray& u) {
  communicate_array_2d_ext_async(u);
  return;
}
#else
//-----------------------------------------------------------------------
void EW::communicate_array_2d_ext(Sarray& u) {
  SW4_MARK_FUNCTION;
  // communicate_array_2d_ext_async(u);
  // return;
  REQUIRE2(
      u.m_nc == 1,
      "Communicate array 2d ext, only implemented for one-component arrays");
  int g = mNumberOfGrids - 1;
  int ie = m_iEnd[g] + m_ext_ghost_points,
      ib = m_iStart[g] - m_ext_ghost_points;
  int je = m_jEnd[g] + m_ext_ghost_points,
      jb = m_jStart[g] - m_ext_ghost_points;

  MPI_Status status;
  int xtag1 = 345;
  int xtag2 = 346;
  int ytag1 = 347;
  int ytag2 = 348;
  int k = 1;
  int extpadding = m_ppadding + m_ext_ghost_points;
  // X-direction communication
  MPI_Sendrecv(&u(1, ie - (2 * extpadding - 1), jb, k), 1,
               m_send_type_2dfinest_ext[0], m_neighbor[1], xtag1,
               &u(1, ib, jb, k), 1, m_send_type_2dfinest_ext[0], m_neighbor[0],
               xtag1, m_cartesian_communicator, &status);
  MPI_Sendrecv(&u(1, ib + extpadding, jb, k), 1, m_send_type_2dfinest_ext[0],
               m_neighbor[0], xtag2, &u(1, ie - (extpadding - 1), jb, k), 1,
               m_send_type_2dfinest_ext[0], m_neighbor[1], xtag2,
               m_cartesian_communicator, &status);

  // Y-direction communication
  MPI_Sendrecv(&u(1, ib, je - (2 * extpadding - 1), k), 1,
               m_send_type_2dfinest_ext[1], m_neighbor[3], ytag1,
               &u(1, ib, jb, k), 1, m_send_type_2dfinest_ext[1], m_neighbor[2],
               ytag1, m_cartesian_communicator, &status);
  MPI_Sendrecv(&u(1, ib, jb + extpadding, k), 1, m_send_type_2dfinest_ext[1],
               m_neighbor[2], ytag2, &u(1, ib, je - (extpadding - 1), k), 1,
               m_send_type_2dfinest_ext[1], m_neighbor[3], ytag2,
               m_cartesian_communicator, &status);
}
#endif
//-----------------------------------------------------------------------
void EW::communicate_array_2d_ext_async(Sarray& u) {
  SW4_MARK_FUNCTION;
  REQUIRE2(
      u.m_nc == 1,
      "Communicate array 2d ext, only implemented for one-component arrays");
  int g = mNumberOfGrids - 1;
  int ie = m_iEnd[g] + m_ext_ghost_points,
      ib = m_iStart[g] - m_ext_ghost_points;
  int je = m_jEnd[g] + m_ext_ghost_points,
      jb = m_jStart[g] - m_ext_ghost_points;

  MPI_Status status;
  int xtag1 = 345;
  int xtag2 = 346;
  int ytag1 = 347;
  int ytag2 = 348;
  int k = 1;
  int extpadding = m_ppadding + m_ext_ghost_points;
  // X-direction communication
  AMPI_Sendrecv(&u(1, ie - (2 * extpadding - 1), jb, k), 1,
                send_type_2dfinest_ext[0], m_neighbor[1], xtag1,
                &u(1, ib, jb, k), 1, send_type_2dfinest_ext[0], m_neighbor[0],
                xtag1, bufs_type_2dfinest_ext[0], m_cartesian_communicator,
                &status);
  AMPI_Sendrecv(&u(1, ib + extpadding, jb, k), 1, send_type_2dfinest_ext[0],
                m_neighbor[0], xtag2, &u(1, ie - (extpadding - 1), jb, k), 1,
                send_type_2dfinest_ext[0], m_neighbor[1], xtag2,
                bufs_type_2dfinest_ext[0], m_cartesian_communicator, &status);

  // Y-direction communication
  AMPI_Sendrecv(&u(1, ib, je - (2 * extpadding - 1), k), 1,
                send_type_2dfinest_ext[1], m_neighbor[3], ytag1,
                &u(1, ib, jb, k), 1, send_type_2dfinest_ext[1], m_neighbor[2],
                ytag1, bufs_type_2dfinest_ext[1], m_cartesian_communicator,
                &status);
  AMPI_Sendrecv(&u(1, ib, jb + extpadding, k), 1, send_type_2dfinest_ext[1],
                m_neighbor[2], ytag2, &u(1, ib, je - (extpadding - 1), k), 1,
                send_type_2dfinest_ext[1], m_neighbor[3], ytag2,
                bufs_type_2dfinest_ext[1], m_cartesian_communicator, &status);
}

//-----------------------------------------------------------------------
void EW::communicate_array_2d_asym(Sarray& u, int g, int k) {
  SW4_MARK_FUNCTION;
  REQUIRE2(
      u.m_nc == 3,
      "Communicate array 2d asym, only implemented for three-component arrays");
  REQUIRE2(g < m_send_type_2dx3p.size(),
           "Communicate array 2d asym, only implemented for grid=0.."
               << m_send_type_2dx3p.size() - 1 << " but g= " << g);
  int ie = m_iEnd[g], ib = m_iStart[g];
  int je = m_jEnd[g], jb = m_jStart[g];

  MPI_Status status;
  int xtag1 = 345;
  int xtag2 = 346;
  int ytag1 = 347;
  int ytag2 = 348;
  int extupper = ie % 2;
  int extlower = 1 - (ib % 2);

  // X-direction communication
  if (extupper && extlower) {
    MPI_Sendrecv(&u(1, ie - (3), jb, k), 1, m_send_type_2dx3p[g], m_neighbor[1],
                 xtag1, &u(1, ib, jb, k), 1, m_send_type_2dx3p[g],
                 m_neighbor[0], xtag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib + 3, jb, k), 1, m_send_type_2dx1p[g], m_neighbor[0],
                 xtag2, &u(1, ie, jb, k), 1, m_send_type_2dx1p[g],
                 m_neighbor[1], xtag2, m_cartesian_communicator, &status);
  } else if (extupper) {
    MPI_Sendrecv(&u(1, ie - 3, jb, k), 1, m_send_type_2dx3p[g], m_neighbor[1],
                 xtag1, &u(1, ib, jb, k), 1, m_send_type_2dx[g], m_neighbor[0],
                 xtag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib + m_ppadding, jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[0], xtag2, &u(1, ie, jb, k), 1,
                 m_send_type_2dx1p[g], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);
  } else if (extlower) {
    MPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[1], xtag1, &u(1, ib, jb, k), 1,
                 m_send_type_2dx3p[g], m_neighbor[0], xtag1,
                 m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib + 3, jb, k), 1, m_send_type_2dx1p[g], m_neighbor[0],
                 xtag2, &u(1, ie - (m_ppadding - 1), jb, k), 1,
                 m_send_type_2dx[g], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);
  } else {
    MPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[1], xtag1, &u(1, ib, jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[0], xtag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib + m_ppadding, jb, k), 1, m_send_type_2dx[g],
                 m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, k), 1,
                 m_send_type_2dx[g], m_neighbor[1], xtag2,
                 m_cartesian_communicator, &status);
  }

  extupper = je % 2;
  extlower = 1 - (jb % 2);
  // Y-direction communication

  if (extupper && extlower) {
    MPI_Sendrecv(&u(1, ib, je - (3), k), 1, m_send_type_2dy3p[g], m_neighbor[3],
                 ytag1, &u(1, ib, jb, k), 1, m_send_type_2dy3p[g],
                 m_neighbor[2], ytag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib, jb + 3, k), 1, m_send_type_2dy1p[g], m_neighbor[2],
                 ytag2, &u(1, ib, je, k), 1, m_send_type_2dy1p[g],
                 m_neighbor[3], ytag2, m_cartesian_communicator, &status);
  } else if (extupper) {
    MPI_Sendrecv(&u(1, ib, je - 3, k), 1, m_send_type_2dy3p[g], m_neighbor[3],
                 ytag1, &u(1, ib, jb, k), 1, m_send_type_2dy[g], m_neighbor[2],
                 ytag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib, jb + m_ppadding, k), 1, m_send_type_2dy[g],
                 m_neighbor[2], ytag2, &u(1, ib, je, k), 1,
                 m_send_type_2dy1p[g], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
  } else if (extlower) {
    MPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), k), 1, m_send_type_2dy[g],
                 m_neighbor[3], ytag1, &u(1, ib, jb, k), 1,
                 m_send_type_2dy3p[g], m_neighbor[2], ytag1,
                 m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib, jb + 3, k), 1, m_send_type_2dy1p[g], m_neighbor[2],
                 ytag2, &u(1, ib, je - (m_ppadding - 1), k), 1,
                 m_send_type_2dy[g], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
  } else {
    MPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), k), 1, m_send_type_2dy[g],
                 m_neighbor[3], ytag1, &u(1, ib, jb, k), 1, m_send_type_2dy[g],
                 m_neighbor[2], ytag1, m_cartesian_communicator, &status);
    MPI_Sendrecv(&u(1, ib, jb + m_ppadding, k), 1, m_send_type_2dy[g],
                 m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), k), 1,
                 m_send_type_2dy[g], m_neighbor[3], ytag2,
                 m_cartesian_communicator, &status);
  }
}

//-----------------------------------------------------------------------
void EW::communicate_array_2dfinest(Sarray& u) {
  SW4_MARK_FUNCTION;
  REQUIRE2(
      u.m_nc == 1,
      "Communicate array 2dfinest, only implemented for one-component arrays");

  int ie = u.m_ie, ib = u.m_ib, je = u.m_je, jb = u.m_jb;
  REQUIRE2(ib == m_iStart[mNumberOfGrids - 1] &&
               ie == m_iEnd[mNumberOfGrids - 1] &&
               jb == m_jStart[mNumberOfGrids - 1] &&
               je == m_jEnd[mNumberOfGrids - 1],
           "Communicate array 2d: Can only use it on the finest grid, grid "
           "sizes don't match");

  MPI_Status status;
  int xtag1 = 345;
  int xtag2 = 346;
  int ytag1 = 347;
  int ytag2 = 348;
  // X-direction communication
  MPI_Sendrecv(&u(ie - (2 * m_ppadding - 1), jb, 1), 1, m_send_type_2dfinest[0],
               m_neighbor[1], xtag1, &u(ib, jb, 1), 1, m_send_type_2dfinest[0],
               m_neighbor[0], xtag1, m_cartesian_communicator, &status);
  MPI_Sendrecv(&u(ib + m_ppadding, jb, 1), 1, m_send_type_2dfinest[0],
               m_neighbor[0], xtag2, &u(ie - (m_ppadding - 1), jb, 1), 1,
               m_send_type_2dfinest[0], m_neighbor[1], xtag2,
               m_cartesian_communicator, &status);
  // Y-direction communication
  MPI_Sendrecv(&u(ib, je - (2 * m_ppadding - 1), 1), 1, m_send_type_2dfinest[1],
               m_neighbor[3], ytag1, &u(ib, jb, 1), 1, m_send_type_2dfinest[1],
               m_neighbor[2], ytag1, m_cartesian_communicator, &status);
  MPI_Sendrecv(&u(ib, jb + m_ppadding, 1), 1, m_send_type_2dfinest[1],
               m_neighbor[2], ytag2, &u(ib, je - (m_ppadding - 1), 1), 1,
               m_send_type_2dfinest[1], m_neighbor[3], ytag2,
               m_cartesian_communicator, &status);
}
void EW::make_type(vector<std::tuple<int, int, int>>& send_type,
                   vector<std::tuple<float_sw4*, float_sw4*>>& bufs_type,
                   int i1, int j1, int k1, int i2, int j2, int k2, int g) {
  send_type[2 * g] = std::make_tuple(i1, j1, k1);
  send_type[2 * g + 1] = std::make_tuple(i2, j2, k2);
  // std::cout<<"MAKE_TYPE"<<i1<<" "<<j1<<" "<<i2<<" "<<j2<<"\n";
  float_sw4* tbuf =
      SW4_NEW(mpi_buffer_space, float_sw4[i1 * j1 * 4 + i2 * j2 * 4]);
  bufs_type[4 * g + 0] = std::make_tuple(tbuf, tbuf + i1 * j1);
  bufs_type[4 * g + 1] =
      std::make_tuple(tbuf + 2 * i1 * j1, tbuf + 3 * i1 * j1);
  tbuf += 4 * i1 * j1;
  bufs_type[4 * g + 2] = std::make_tuple(tbuf, tbuf + i2 * j2);
  bufs_type[4 * g + 3] =
      std::make_tuple(tbuf + 2 * i2 * j2, tbuf + 3 * i2 * j2);

  global_variables.buffer_size = std::max<size_t>(
      global_variables.buffer_size, std::max<size_t>(i1 * j1, i2 * j2));
}
void EW::make_type_2d(vector<std::tuple<int, int, int>>& send_type,
                      vector<std::tuple<float_sw4*, float_sw4*>>& bufs_type,
                      int i1, int j1, int k1, int g) {
  send_type[g] = std::make_tuple(i1, j1, k1);

  float_sw4* tbuf = SW4_NEW(Space::Pinned, float_sw4[i1 * j1 * 2]);
  bufs_type[g] = std::make_tuple(tbuf, tbuf + i1 * j1);

  global_variables.buffer_size =
      std::max<size_t>(global_variables.buffer_size, i1 * j1);
}
void EW::communicate_array_async(Sarray& u, int grid) {
  SW4_MARK_FUNCTION;
  // u.forceprefetch();
  // u.prefetch(cudaCpuDeviceId);
  REQUIRE2(u.m_nc == 1 || u.m_nc == 3 || u.m_nc == 4,
           "Communicate array, only implemented for nc=1,3, and 4 "
               << " nc = " << u.m_nc);
  int ie = u.m_ie, ib = u.m_ib, je = u.m_je, jb = u.m_jb,
      kb = u.m_kb;  //,ke=u.m_ke;
  MPI_Status status;
#ifdef THREADED_MPI
  const int threaded_mpi = 1;
#else
  //  const int threaded_mpi = 0;
#endif
  if (u.m_nc == 1) {
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    // X-direction communication
    // for (int i=0;i<4;i++)std::cout<<m_neighbor[i]<<" ";
    // std::cout<<"\n";
    //#pragma omp parallel default(shared) if (threaded_mpi)
    //#pragma omp sections
    {
      //#pragma omp section
      AMPI_Sendrecv(&u(ie - (2 * m_ppadding - 1), jb, kb), 1,
                    send_type1[2 * grid], m_neighbor[1], xtag1, &u(ib, jb, kb),
                    1, send_type1[2 * grid], m_neighbor[0], xtag1,
                    bufs_type1[4 * grid], m_cartesian_communicator, &status);
      //#pragma omp section
      AMPI_Sendrecv(&u(ib + m_ppadding, jb, kb), 1, send_type1[2 * grid],
                    m_neighbor[0], xtag2, &u(ie - (m_ppadding - 1), jb, kb), 1,
                    send_type1[2 * grid], m_neighbor[1], xtag2,
                    bufs_type1[4 * grid + 1], m_cartesian_communicator,
                    &status);
      // Y-direction communication
      //#pragma omp section
      AMPI_Sendrecv(&u(ib, je - (2 * m_ppadding - 1), kb), 1,
                    send_type1[2 * grid + 1], m_neighbor[3], ytag1,
                    &u(ib, jb, kb), 1, send_type1[2 * grid + 1], m_neighbor[2],
                    ytag1, bufs_type1[4 * grid + 2], m_cartesian_communicator,
                    &status);
      //#pragma omp section
      AMPI_Sendrecv(&u(ib, jb + m_ppadding, kb), 1, send_type1[2 * grid + 1],
                    m_neighbor[2], ytag2, &u(ib, je - (m_ppadding - 1), kb), 1,
                    send_type1[2 * grid + 1], m_neighbor[3], ytag2,
                    bufs_type1[4 * grid + 3], m_cartesian_communicator,
                    &status);
    }
  } else if (u.m_nc == 3) {
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    //#pragma omp parallel default(shared) if (threaded_mpi)
    //#pragma omp sections
    {
      // X-direction communication
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, kb), 1,
                    send_type3[2 * grid], m_neighbor[1], xtag1,
                    &u(1, ib, jb, kb), 1, send_type3[2 * grid], m_neighbor[0],
                    xtag1, bufs_type3[4 * grid], m_cartesian_communicator,
                    &status);
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ib + m_ppadding, jb, kb), 1, send_type3[2 * grid],
                    m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, kb),
                    1, send_type3[2 * grid], m_neighbor[1], xtag2,
                    bufs_type3[4 * grid + 1], m_cartesian_communicator,
                    &status);
      // Y-direction communication
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), kb), 1,
                    send_type3[2 * grid + 1], m_neighbor[3], ytag1,
                    &u(1, ib, jb, kb), 1, send_type3[2 * grid + 1],
                    m_neighbor[2], ytag1, bufs_type3[4 * grid + 2],
                    m_cartesian_communicator, &status);
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ib, jb + m_ppadding, kb), 1, send_type3[2 * grid + 1],
                    m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), kb),
                    1, send_type3[2 * grid + 1], m_neighbor[3], ytag2,
                    bufs_type3[4 * grid + 3], m_cartesian_communicator,
                    &status);
    }
  } else if (u.m_nc == 4) {
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    //#pragma omp parallel default(shared) if (threaded_mpi)
    //#pragma omp sections
    {
      // X-direction communication
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, kb), 1,
                    send_type4[2 * grid], m_neighbor[1], xtag1,
                    &u(1, ib, jb, kb), 1, send_type4[2 * grid], m_neighbor[0],
                    xtag1, bufs_type4[4 * grid], m_cartesian_communicator,
                    &status);
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ib + m_ppadding, jb, kb), 1, send_type4[2 * grid],
                    m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, kb),
                    1, send_type4[2 * grid], m_neighbor[1], xtag2,
                    bufs_type4[4 * grid + 1], m_cartesian_communicator,
                    &status);
      //#pragma omp section
      // Y-direction communication
      AMPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), kb), 1,
                    send_type4[2 * grid + 1], m_neighbor[3], ytag1,
                    &u(1, ib, jb, kb), 1, send_type4[2 * grid + 1],
                    m_neighbor[2], ytag1, bufs_type4[4 * grid + 2],
                    m_cartesian_communicator, &status);
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ib, jb + m_ppadding, kb), 1, send_type4[2 * grid + 1],
                    m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), kb),
                    1, send_type4[2 * grid + 1], m_neighbor[3], ytag2,
                    bufs_type4[4 * grid + 3], m_cartesian_communicator,
                    &status);
    }
  } else if (u.m_nc == 21) {
    std::cout << "This is happening\n";
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    //#pragma omp parallel default(shared) if (threaded_mpi)
    //#pragma omp sections
    {
      // X-direction communication
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ie - (2 * m_ppadding - 1), jb, kb), 1,
                    send_type21[2 * grid], m_neighbor[1], xtag1,
                    &u(1, ib, jb, kb), 1, send_type21[2 * grid], m_neighbor[0],
                    xtag1, bufs_type21[4 * grid], m_cartesian_communicator,
                    &status);
      //#pragma omp section
      AMPI_Sendrecv(&u(1, ib + m_ppadding, jb, kb), 1, send_type21[2 * grid],
                    m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, kb),
                    1, send_type21[2 * grid], m_neighbor[1], xtag2,
                    bufs_type21[4 * grid + 1], m_cartesian_communicator,
                    &status);
      //#pragma omp section
      // Y-direction communication
      AMPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), kb), 1,
                    send_type21[2 * grid + 1], m_neighbor[3], ytag1,
                    &u(1, ib, jb, kb), 1, send_type21[2 * grid + 1],
                    m_neighbor[2], ytag1, bufs_type21[4 * grid + 2],
                    m_cartesian_communicator, &status);
      //#pragma omp section
      AMPI_Sendrecv(
          &u(1, ib, jb + m_ppadding, kb), 1, send_type21[2 * grid + 1],
          m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), kb), 1,
          send_type21[2 * grid + 1], m_neighbor[3], ytag2,
          bufs_type21[4 * grid + 3], m_cartesian_communicator, &status);
    }
  }
  // u.prefetch();
}
void EW::AMPI_Sendrecv(float_sw4* a, int scount,
                       std::tuple<int, int, int>& sendt, int sendto, int stag,
                       float_sw4* b, int rcount,
                       std::tuple<int, int, int>& recvt, int recvfrom, int rtag,
                       std::tuple<float_sw4*, float_sw4*>& buf, MPI_Comm comm,
                       MPI_Status* status) {
  SW4_MARK_FUNCTION;
  SW4_MARK_BEGIN("THE REST");
  MPI_Request send_req = MPI_REQUEST_NULL, recv_req = MPI_REQUEST_NULL;

  int recv_count = std::get<0>(recvt) * std::get<1>(recvt);
  int send_count = std::get<0>(sendt) * std::get<1>(sendt);

  SW4_MARK_END("THE REST");
#if defined(ENABLE_MPI_TIMING_BARRIER)
#if defined(SW4_TRACK_MPI)
  std::chrono::high_resolution_clock::time_point t1, t2;
  t1 = SW4_CHRONO_NOW;
#endif
  MPI_Barrier(MPI_COMM_WORLD);
#if defined(SW4_TRACK_MPI)
  t2 = SW4_CHRONO_NOW;
  coll_sm.insert(0, SW4_CHRONO_DURATION_US(t1, t2));
#endif
#endif
  SW4_MARK_BEGIN("MPI_SENDRECV_ACTUAL");

  if (sendto != MPI_PROC_NULL) {
    getbuffer_device(a, std::get<0>(buf), sendt, true);
    // getbuffer_host(a, std::get<0>(buf), sendt);
  }
  // std::cout<<"send_count "<<send_count<<" recv_count "<<recv_count<<"\n";

#if defined(SW4_TRACK_MPI)
#if defined(ENABLE_MPI_TIMING_BARRIER)
  t1 = SW4_CHRONO_NOW;
  MPI_Barrier(MPI_COMM_WORLD);
  t2 = SW4_CHRONO_NOW;
  coll_sm.insert(1, SW4_CHRONO_DURATION_US(t1, t2));
#endif
  SYNC_STREAM;  // Avoid adding the buffering time to the MPI bandwdth
  t1 = SW4_CHRONO_NOW;
#endif
  if (recvfrom != MPI_PROC_NULL)
    if (MPI_Irecv(std::get<1>(buf), recv_count, MPI_DOUBLE, recvfrom, rtag,
                  comm, &recv_req) != MPI_SUCCESS)
      std::cerr << "MPI_Irecv failed in EW::AMPI_Sendrecv\n";

  if (sendto != MPI_PROC_NULL) {
    SYNC_STREAM;
    // getbuffer_device(a,std::get<0>(buf),sendt);
    // std::cout<<"SENDING :: "<<sendto<<" ";
    // for(int i=0;i<10;i++) std::cout<<std::get<0>(buf)[i]<<" ";
    // std::cout<<"\n";
    if (MPI_Isend(std::get<0>(buf), send_count, MPI_DOUBLE, sendto, stag, comm,
                  &send_req) != MPI_SUCCESS)
      std::cerr << "MPI_Isend failed in EW::AMPI_Sendrecv\n";
  }
  MPI_Status send_status, recv_status;

  SW4_MARK_BEGIN("MPI_RECV_WAIT");
  if (recvfrom != MPI_PROC_NULL) {
    if (MPI_Wait(&recv_req, &recv_status) != MPI_SUCCESS)
      std::cerr << "MPI_WAIT RECV FAILED IN AMPI_SENDrecv\n";

#if defined(SW4_TRACK_MPI)
    t2 = SW4_CHRONO_NOW;
#endif

    putbuffer_device(b, std::get<1>(buf), recvt, true);
    // putbuffer_host(b, std::get<1>(buf), recvt);
    // std::cout<<"RECEIVING :: "<<recvfrom<<" ";
    // for(int i=0;i<10;i++) std::cout<<std::get<1>(buf)[i]<<" ";
    // std::cout<<"\n";
  }
  SW4_MARK_END("MPI_RECV_WAIT");

  SW4_MARK_BEGIN("MPI_SEND_WAIT");
  if (sendto != MPI_PROC_NULL)
    if (MPI_Wait(&send_req, &send_status) != MPI_SUCCESS)
      std::cerr << "MPI_WAIT SEND FAILED IN AMPI_SENDrecv\n";

#if defined(SW4_TRACK_MPI)
  if (recvfrom == MPI_PROC_NULL) t2 = SW4_CHRONO_NOW;  // Timing only the send
#endif
  SW4_MARK_END("MPI_SEND_WAIT");

  SYNC_STREAM;  // For the putbuffer_device
  SW4_MARK_END("MPI_SENDRECV_ACTUAL");

#if defined(SW4_TRACK_MPI)
  size_t size = 0;
  if (sendto != MPI_PROC_NULL) size += send_count;
  if (recvfrom != MPI_PROC_NULL) size += recv_count;
  // auto got = mpi_times.find(size);
  // if (got==mpi_times.end()) {
  //   mpi_times[size]=SW4_CHRONO_DURATION_US(t1,t2);
  //   mpi_count[size]=1;
  //   //std::cout<<"INIT "<<SW4_CHRONO_DURATION_US(t1,t2)<<"\n";
  // } else {
  //   mpi_times[size]+=SW4_CHRONO_DURATION_US(t1,t2);
  //   mpi_count[size]++;
  //   //std::cout<<"UPDATE "<<SW4_CHRONO_DURATION_US(t1,t2)<<"\n";
  // }
  sm.insert(size, SW4_CHRONO_DURATION_US(t1, t2));
#endif
}
void EW::getbuffer_device(float_sw4* data, float_sw4* buf,
                          std::tuple<int, int, int>& mtype, bool async) {
  SW4_MARK_FUNCTION;
  int count = std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);
  // std::cout<<"getbuffer_device...";
  // PREFETCH(buf);
  // forall(0,count,[=] RAJA_DEVICE(int i){
  // RAJA::forall<DEFAULT_LOOP1> (0,count,[=] RAJA_DEVICE(int i){

  // Messages are greater than 2K for Hayward h=200 on 16 ranks.
  // for larger messages, host copy is slower.
  // if (bl*count*8>2048){

#ifdef SW4_STAGED_MPI_BUFFERS

  // Code for the staged option. Local device buffer + copy to pinned host
  // buffer A single large local buffer is used. Allocated in the 2D and 3D
  // setup routines
  float_sw4* lbuf = global_variables.device_buffer;
#ifndef UNRAJA
  RAJA::RangeSegment k_range(0, bl);
  RAJA::RangeSegment i_range(0, count);
  RAJA::kernel<BUFFER_POL>(RAJA::make_tuple(k_range, i_range),
                           [=] RAJA_DEVICE(int k, int i) {
                             lbuf[k + i * bl] = data[i * stride + k];
                           });
#else
  Range<16> k_range(0, bl);
  Range<16> i_range(0, count);
  forall2async(i_range, k_range, [=] RAJA_DEVICE(int i, int k) {
    lbuf[k + i * bl] = data[i * stride + k];
  });
#endif
  if (async)
    SW4_CheckDeviceError(
        cudaMemcpyAsync(buf, lbuf, count * bl * 8, cudaMemcpyDeviceToHost, 0));
  else
    SW4_CheckDeviceError(
        cudaMemcpy(buf, lbuf, count * bl * 8, cudaMemcpyDeviceToHost));
#else

  // Code for PINNED,DEVICE AND MANAGED BUFFERS
#ifndef UNRAJA
  RAJA::RangeSegment k_range(0, bl);
  RAJA::RangeSegment i_range(0, count);
  RAJA::kernel<BUFFER_POL>(RAJA::make_tuple(k_range, i_range),
                           [=] RAJA_DEVICE(int k, int i) {
                             buf[k + i * bl] = data[i * stride + k];
                           });
#else
  Range<16> k_range(0, bl);
  Range<16> i_range(0, count);
  forall2async(i_range, k_range, [=] RAJA_DEVICE(int i, int k) {
    buf[k + i * bl] = data[i * stride + k];
  });
#endif

  if (!async) {
    SYNC_STREAM;
  }

#endif

  // SW4_PEEK;
  // SYNC_STREAM;

  // } else {
  //   std::cout<<bl*count*8<<"\ bytes n";
  //   for(int i=0;i<count;i++) for(int k=0;k<bl;k++)
  //   buf[k+i*bl]=data[i*stride+k];
  // }
  // std::cout<<"Done\n";
}
void EW::getbuffer_host(float_sw4* data, float_sw4* buf,
                        std::tuple<int, int, int>& mtype) {
  SW4_MARK_FUNCTION;
  int count = std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);

  // std::cout<<bl*count*8<<"\ bytes n";
  for (int i = 0; i < count; i++)
    for (int k = 0; k < bl; k++) buf[k + i * bl] = data[i * stride + k];
  // std::cout<<"Done\n";
}

void EW::putbuffer_device(float_sw4* data, float_sw4* buf,
                          std::tuple<int, int, int>& mtype, bool async) {
  SW4_MARK_FUNCTION;
  int count = std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);
  // std::cout<<"putbuffer_device...";
  // PREFETCHFORCED(buf);

  // The STAGED option
#ifdef SW4_STAGED_MPI_BUFFERS
  float_sw4* lbuf = global_variables.device_buffer;
  SW4_CheckDeviceError(
      cudaMemcpyAsync(lbuf, buf, count * bl * 8, cudaMemcpyHostToDevice, 0));
#ifndef UNRAJA
  Range<16> k_range(0, bl);
  Range<16> i_range(0, count);
  forall2async(i_range, k_range, [=] RAJA_DEVICE(int i, int k) {
    data[i * stride + k] = lbuf[k + i * bl];
  });
#else
  RAJA::RangeSegment k_range(0, bl);
  RAJA::RangeSegment i_range(0, count);
  RAJA::kernel<BUFFER_POL>(RAJA::make_tuple(k_range, i_range),
                           [=] RAJA_DEVICE(int k, int i) {
                             // RAJA::forall<DEFAULT_LOOP1> (0,count,[=]
                             // RAJA_DEVICE(int i){
                             data[i * stride + k] = lbuf[k + i * bl];
                           });
#endif

  // The PINNED, DEVICE and MANAGED cases
#else

#ifndef UNRAJA
  Range<16> k_range(0, bl);
  Range<16> i_range(0, count);
  forall2async(i_range, k_range, [=] RAJA_DEVICE(int i, int k) {
    data[i * stride + k] = buf[k + i * bl];
  });
#else
  RAJA::RangeSegment k_range(0, bl);
  RAJA::RangeSegment i_range(0, count);
  RAJA::kernel<BUFFER_POL>(RAJA::make_tuple(k_range, i_range),
                           [=] RAJA_DEVICE(int k, int i) {
                             // RAJA::forall<DEFAULT_LOOP1> (0,count,[=]
                             // RAJA_DEVICE(int i){
                             data[i * stride + k] = buf[k + i * bl];
                           });
#endif

#endif
  // SW4_PEEK;
  // SYNC_STREAM;
  if (!async) {
    SYNC_STREAM;
  }
  // std::cout<<"Done\n";
}
void EW::putbuffer_host(float_sw4* data, float_sw4* buf,
                        std::tuple<int, int, int>& mtype) {
  SW4_MARK_FUNCTION;
  int count = std::get<0>(mtype);
  int bl = std::get<1>(mtype);
  int stride = std::get<2>(mtype);
  // std::cout<<"putbuffer_device...";

  for (int i = 0; i < count; i++)
    for (int k = 0; k < bl; k++) data[i * stride + k] = buf[k + i * bl];
}

void EW::communicate_array_2d_async(Sarray& u, int g, int k) {
  SW4_MARK_FUNCTION;
  REQUIRE2(u.m_nc == 3,
           "Communicate array 2d, only implemented for three-component arrays");
  REQUIRE2(g < m_send_type_2dx.size(),
           "Communicate array 2d, only implemented for grid=0.."
               << m_send_type_2dx.size() - 1 << " but g= " << g);
  int ie = m_iEnd[g], ib = m_iStart[g];
  int je = m_jEnd[g], jb = m_jStart[g];

  MPI_Status status;
  int xtag1 = 345;
  int xtag2 = 346;
  int ytag1 = 347;
  int ytag2 = 348;

  if (m_croutines && u.m_ke - u.m_kb + 1 != 1) {
    SW4_MARK_BEGIN("comm_array_2d_async::ALLOC");
    Sarray u2d(3, u.m_ib, u.m_ie, u.m_jb, u.m_je, k, k, __FILE__, __LINE__);
    // u2d.GetAtt(__FILE__,__LINE__);
    // u.GetAtt(__FILE__,__LINE__);
    SW4_MARK_END("comm_array_2d_async::ALLOC");
    u2d.copy_kplane(u, k);
    SW4_MARK_BEGIN("comm_array_2d_async::MPI-1");
    // X-direction communication
    AMPI_Sendrecv2(&u2d(1, ie - (2 * m_ppadding - 1), jb, k), 1,
                   send_type_2dx[g], m_neighbor[1], xtag1, &u2d(1, ib, jb, k),
                   1, send_type_2dx[g], m_neighbor[0], xtag1, bufs_type_2dx[g],
                   m_cartesian_communicator, &status);
    AMPI_Sendrecv2(&u2d(1, ib + m_ppadding, jb, k), 1, send_type_2dx[g],
                   m_neighbor[0], xtag2, &u2d(1, ie - (m_ppadding - 1), jb, k),
                   1, send_type_2dx[g], m_neighbor[1], xtag2, bufs_type_2dx[g],
                   m_cartesian_communicator, &status);
    // Y-direction communication
    AMPI_Sendrecv2(&u2d(1, ib, je - (2 * m_ppadding - 1), k), 1,
                   send_type_2dy[g], m_neighbor[3], ytag1, &u2d(1, ib, jb, k),
                   1, send_type_2dy[g], m_neighbor[2], ytag1, bufs_type_2dy[g],
                   m_cartesian_communicator, &status);
    AMPI_Sendrecv2(&u2d(1, ib, jb + m_ppadding, k), 1, send_type_2dy[g],
                   m_neighbor[2], ytag2, &u2d(1, ib, je - (m_ppadding - 1), k),
                   1, send_type_2dy[g], m_neighbor[3], ytag2, bufs_type_2dy[g],
                   m_cartesian_communicator, &status);
    SW4_MARK_END("comm_array_2d_async::MPI-1");
    u.copy_kplane(u2d, k);
  } else {
    // X-direction communication
    SYNC_STREAM;  // For the set above, the sync in copy_kplane takes care of
                  // races.
    SW4_MARK_BEGIN("comm_array_2d_async::MPI-2");
    AMPI_Sendrecv2(&u(1, ie - (2 * m_ppadding - 1), jb, k), 1, send_type_2dx[g],
                   m_neighbor[1], xtag1, &u(1, ib, jb, k), 1, send_type_2dx[g],
                   m_neighbor[0], xtag1, bufs_type_2dx[g],
                   m_cartesian_communicator, &status);
    AMPI_Sendrecv2(&u(1, ib + m_ppadding, jb, k), 1, send_type_2dx[g],
                   m_neighbor[0], xtag2, &u(1, ie - (m_ppadding - 1), jb, k), 1,
                   send_type_2dx[g], m_neighbor[1], xtag2, bufs_type_2dx[g],
                   m_cartesian_communicator, &status);

    // Y-direction communication
    AMPI_Sendrecv2(&u(1, ib, je - (2 * m_ppadding - 1), k), 1, send_type_2dy[g],
                   m_neighbor[3], ytag1, &u(1, ib, jb, k), 1, send_type_2dy[g],
                   m_neighbor[2], ytag1, bufs_type_2dy[g],
                   m_cartesian_communicator, &status);
    AMPI_Sendrecv2(&u(1, ib, jb + m_ppadding, k), 1, send_type_2dy[g],
                   m_neighbor[2], ytag2, &u(1, ib, je - (m_ppadding - 1), k), 1,
                   send_type_2dy[g], m_neighbor[3], ytag2, bufs_type_2dy[g],
                   m_cartesian_communicator, &status);
    SW4_MARK_END("comm_array_2d_async::MPI-2");
  }
}

void EW::communicate_array_2d_async_memo(Sarray& u, int g, int k) {
  SW4_MARK_FUNCTION;
  REQUIRE2(u.m_nc == 3,
           "Communicate array 2d, only implemented for three-component arrays");
  REQUIRE2(g < m_send_type_2dx.size(),
           "Communicate array 2d, only implemented for grid=0.."
               << m_send_type_2dx.size() - 1 << " but g= " << g);
  int ie = m_iEnd[g], ib = m_iStart[g];
  int je = m_jEnd[g], jb = m_jStart[g];

  MPI_Status status;
  int xtag1 = 345;
  int xtag2 = 346;
  int ytag1 = 347;
  int ytag2 = 348;

  if (m_croutines && u.m_ke - u.m_kb + 1 != 1) {
    SW4_MARK_BEGIN("comm_array_2d_async::ALLOC");
    Sarray u2d(3, u.m_ib, u.m_ie, u.m_jb, u.m_je, k, k, __FILE__, __LINE__);
    SW4_MARK_END("comm_array_2d_async::ALLOC");
    u2d.copy_kplane(u, k);
    SW4_MARK_BEGIN("comm_array_2d_async::MPI-1");
    // X-direction communication
    AMPI_Sendrecv(memoize(u2d, 1, ie - (2 * m_ppadding - 1), jb, k), 1,
                  send_type_2dx[g], m_neighbor[1], xtag1,
                  memoize(u2d, 1, ib, jb, k), 1, send_type_2dx[g],
                  m_neighbor[0], xtag1, bufs_type_2dx[g],
                  m_cartesian_communicator, &status);
    AMPI_Sendrecv(memoize(u2d, 1, ib + m_ppadding, jb, k), 1, send_type_2dx[g],
                  m_neighbor[0], xtag2,
                  memoize(u2d, 1, ie - (m_ppadding - 1), jb, k), 1,
                  send_type_2dx[g], m_neighbor[1], xtag2, bufs_type_2dx[g],
                  m_cartesian_communicator, &status);
    // Y-direction communication
    AMPI_Sendrecv(memoize(u2d, 1, ib, je - (2 * m_ppadding - 1), k), 1,
                  send_type_2dy[g], m_neighbor[3], ytag1,
                  memoize(u2d, 1, ib, jb, k), 1, send_type_2dy[g],
                  m_neighbor[2], ytag1, bufs_type_2dy[g],
                  m_cartesian_communicator, &status);
    AMPI_Sendrecv(memoize(u2d, 1, ib, jb + m_ppadding, k), 1, send_type_2dy[g],
                  m_neighbor[2], ytag2,
                  memoize(u2d, 1, ib, je - (m_ppadding - 1), k), 1,
                  send_type_2dy[g], m_neighbor[3], ytag2, bufs_type_2dy[g],
                  m_cartesian_communicator, &status);
    SW4_MARK_END("comm_array_2d_async::MPI-1");
    u.copy_kplane(u2d, k);
  } else {
    // X-direction communication
    SW4_MARK_BEGIN("comm_array_2d_async::MPI-2");
    AMPI_Sendrecv(memoize(u, 1, ie - (2 * m_ppadding - 1), jb, k), 1,
                  send_type_2dx[g], m_neighbor[1], xtag1,
                  memoize(u, 1, ib, jb, k), 1, send_type_2dx[g], m_neighbor[0],
                  xtag1, bufs_type_2dx[g], m_cartesian_communicator, &status);
    AMPI_Sendrecv(&u(1, ib + m_ppadding, jb, k), 1, send_type_2dx[g],
                  m_neighbor[0], xtag2,
                  memoize(u, 1, ie - (m_ppadding - 1), jb, k), 1,
                  send_type_2dx[g], m_neighbor[1], xtag2, bufs_type_2dx[g],
                  m_cartesian_communicator, &status);

    // Y-direction communication
    AMPI_Sendrecv(&u(1, ib, je - (2 * m_ppadding - 1), k), 1, send_type_2dy[g],
                  m_neighbor[3], ytag1, memoize(u, 1, ib, jb, k), 1,
                  send_type_2dy[g], m_neighbor[2], ytag1, bufs_type_2dy[g],
                  m_cartesian_communicator, &status);
    AMPI_Sendrecv(&u(1, ib, jb + m_ppadding, k), 1, send_type_2dy[g],
                  m_neighbor[2], ytag2,
                  memoize(u, 1, ib, je - (m_ppadding - 1), k), 1,
                  send_type_2dy[g], m_neighbor[3], ytag2, bufs_type_2dy[g],
                  m_cartesian_communicator, &status);
    SW4_MARK_END("comm_array_2d_async::MPI-2");
  }
}
void EW::AMPI_Sendrecv2(float_sw4* a, int scount,
                        std::tuple<int, int, int>& sendt, int sendto, int stag,
                        float_sw4* b, int rcount,
                        std::tuple<int, int, int>& recvt, int recvfrom,
                        int rtag, std::tuple<float_sw4*, float_sw4*>& buf,
                        MPI_Comm comm, MPI_Status* status) {
  SW4_MARK_FUNCTION;
  SW4_MARK_BEGIN("THE REST2");
  MPI_Request send_req = MPI_REQUEST_NULL, recv_req = MPI_REQUEST_NULL;

  int recv_count = std::get<0>(recvt) * std::get<1>(recvt);
  int send_count = std::get<0>(sendt) * std::get<1>(sendt);

  SW4_MARK_END("THE REST2");
#if defined(ENABLE_MPI_TIMING_BARRIER)
#if defined(SW4_TRACK_MPI)
  std::chrono::high_resolution_clock::time_point t1, t2;
  t1 = SW4_CHRONO_NOW;
#endif
  MPI_Barrier(MPI_COMM_WORLD);
#if defined(SW4_TRACK_MPI)
  t2 = SW4_CHRONO_NOW;
  coll_sm.insert(2, SW4_CHRONO_DURATION_US(t1, t2));
#endif
#endif
  SW4_MARK_BEGIN("MPI_SENDRECV_ACTUAL2");

  if (sendto != MPI_PROC_NULL) {
#if defined(SW4_TRACK_MPI)
    auto t1 = SW4_CHRONO_NOW;
    {
#endif
      // getbuffer_device(a,std::get<0>(buf),sendt,true);
      getbuffer_host(a, std::get<0>(buf), sendt);
#if defined(SW4_TRACK_MPI)
      auto t2 = SW4_CHRONO_NOW;
      size_t size = 0;
      if (sendto != MPI_PROC_NULL) size += send_count;
      if (recvfrom != MPI_PROC_NULL) size += recv_count;
      host_sm.insert(size, SW4_CHRONO_DURATION_US(t1, t2));
    }
#endif
  }
#if defined(SW4_TRACK_MPI)
  // SYNC_STREAM; // Avoid adding the buffering time to the MPI bandwdth
  t1 = SW4_CHRONO_NOW;
#endif
  // std::cout<<"send_count "<<send_count<<" recv_count "<<recv_count<<"\n";
  SW4_MARK_BEGIN("MPI_IRECV");
  if (recvfrom != MPI_PROC_NULL)
    if (MPI_Irecv(std::get<1>(buf), recv_count, MPI_DOUBLE, recvfrom, rtag,
                  comm, &recv_req) != MPI_SUCCESS)
      std::cerr << "MPI_Irecv failed in EW::AMPI_Sendrecv\n";
  SW4_MARK_END("MPI_IRECV");
  if (sendto != MPI_PROC_NULL) {
    // getbuffer_device(a,std::get<0>(buf),sendt);
    // std::cout<<"SENDING :: "<<sendto<<" ";
    // for(int i=0;i<10;i++) std::cout<<std::get<0>(buf)[i]<<" ";
    // std::cout<<"\n";
    // SYNC_STREAM;
    if (MPI_Isend(std::get<0>(buf), send_count, MPI_DOUBLE, sendto, stag, comm,
                  &send_req) != MPI_SUCCESS)
      std::cerr << "MPI_Isend failed in EW::AMPI_Sendrecv\n";
  }
  MPI_Status send_status, recv_status;

  SW4_MARK_BEGIN("MPI_RECV_WAIT2");
  if (recvfrom != MPI_PROC_NULL) {
    if (MPI_Wait(&recv_req, &recv_status) != MPI_SUCCESS)
      std::cerr << "MPI_WAIT RECV FAILED IN AMPI_SENDrecv\n";
      // putbuffer_device(b,std::get<1>(buf),recvt,true);

#if defined(SW4_TRACK_MPI)
    t2 = SW4_CHRONO_NOW;
#endif
    {
#if defined(SW4_TRACK_MPI)
      auto t1 = SW4_CHRONO_NOW;
#endif
      putbuffer_host(b, std::get<1>(buf), recvt);
#if defined(SW4_TRACK_MPI)
      auto t2 = SW4_CHRONO_NOW;
      size_t size = 0;
      if (sendto != MPI_PROC_NULL) size += send_count;
      if (recvfrom != MPI_PROC_NULL) size += recv_count;
      host_sm.insert(size, SW4_CHRONO_DURATION_US(t1, t2));
#endif
    }
    // std::cout<<"RECEIVING :: "<<recvfrom<<" ";
    // for(int i=0;i<10;i++) std::cout<<std::get<1>(buf)[i]<<" ";
    // std::cout<<"\n";
  }
  SW4_MARK_END("MPI_RECV_WAIT2");

  SW4_MARK_BEGIN("MPI_SEND_WAIT2");
  if (sendto != MPI_PROC_NULL)
    if (MPI_Wait(&send_req, &send_status) != MPI_SUCCESS)
      std::cerr << "MPI_WAIT SEND FAILED IN AMPI_SENDrecv\n";
  SW4_MARK_END("MPI_SEND_WAIT2");

#if defined(SW4_TRACK_MPI)
  if (recvfrom == MPI_PROC_NULL) t2 = SW4_CHRONO_NOW;  // Timing only the send
#endif

  // SYNC_STREAM;

  SW4_MARK_END("MPI_SENDRECV_ACTUAL2");
#if defined(SW4_TRACK_MPI)
  size_t size = 0;
  if (sendto != MPI_PROC_NULL) size += send_count;
  if (recvfrom != MPI_PROC_NULL) size += recv_count;
  // auto got = mpi_times.find(size);
  // if (got==mpi_times.end()) {
  //   mpi_times2[size]=SW4_CHRONO_DURATION_US(t1,t2);
  //   mpi_count2[size]=1;
  //   //std::cout<<"INIT "<<SW4_CHRONO_DURATION_US(t1,t2)<<"\n";
  // } else {
  //   mpi_times2[size]+=SW4_CHRONO_DURATION_US(t1,t2);
  //   mpi_count2[size]++;
  //   //std::cout<<"UPDATE "<<SW4_CHRONO_DURATION_US(t1,t2)<<"\n";
  // }
  sm2.insert(size, SW4_CHRONO_DURATION_US(t1, t2));
#endif
}
//-----------------------------------------------------------------------
void EW::communicate_array_2d_isurf(Sarray& u, int iSurf) {
  REQUIRE2(
      u.m_nc == 1,
      "Communicate array 2d isurf, only implemented for one-component arrays");
  SYNC_STREAM;
  int g = mNumberOfCartesianGrids + iSurf;
  int ie = m_iEnd[g] + m_ext_ghost_points,
      ib = m_iStart[g] - m_ext_ghost_points;
  int je = m_jEnd[g] + m_ext_ghost_points,
      jb = m_jStart[g] - m_ext_ghost_points;

  MPI_Status status;
  int xtag1 = 345;
  int xtag2 = 346;
  int ytag1 = 347;
  int ytag2 = 348;
  int k = 1;
  int extpadding = m_ppadding + m_ext_ghost_points;
  // X-direction communication
  // std::cout<<"COMM "<<&u(1, ie - (2 * extpadding - 1),jb,k)<<" "<<iSurf<<" "
  //	   <<&u(1, ib, jb, k)<<"\n"<<std::flush;
  // std::cout<<"NEIGHS"<<m_neighbor[0]<<" "<<m_neighbor[1]<<"\n"<<std::flush;

  MPI_Sendrecv(&u(1, ie - (2 * extpadding - 1), jb, k), 1,
               m_send_type_isurfx[iSurf], m_neighbor[1], xtag1,
               &u(1, ib, jb, k), 1, m_send_type_isurfx[iSurf], m_neighbor[0],
               xtag1, m_cartesian_communicator, &status);

  MPI_Sendrecv(&u(1, ib + extpadding, jb, k), 1, m_send_type_isurfx[iSurf],
               m_neighbor[0], xtag2, &u(1, ie - (extpadding - 1), jb, k), 1,
               m_send_type_isurfx[iSurf], m_neighbor[1], xtag2,
               m_cartesian_communicator, &status);

  // Y-direction communication

  MPI_Sendrecv(&u(1, ib, je - (2 * extpadding - 1), k), 1,
               m_send_type_isurfy[iSurf], m_neighbor[3], ytag1,
               &u(1, ib, jb, k), 1, m_send_type_isurfy[iSurf], m_neighbor[2],
               ytag1, m_cartesian_communicator, &status);

  MPI_Sendrecv(&u(1, ib, jb + extpadding, k), 1, m_send_type_isurfy[iSurf],
               m_neighbor[2], ytag2, &u(1, ib, je - (extpadding - 1), k), 1,
               m_send_type_isurfy[iSurf], m_neighbor[3], ytag2,
               m_cartesian_communicator, &status);
}
