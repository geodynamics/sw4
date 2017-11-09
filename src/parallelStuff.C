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

//-----------------------------------------------------------------------
bool EW::proc_decompose_2d( int ni, int nj, int nproc, int proc_max[2] )
{
   // This routine determines a decomposition of nproc processors into
   // a 2D processor array  proc_max[0] x proc_max[1], which gives minimal 
   // communication boundary for a grid with ni x nj points.

   float_sw4 fmin = ni+nj;
   bool first  = true;
   int p1max   = ni/m_ppadding;
   int p2max   = nj/m_ppadding;
   for( int p1 = 1 ; p1 <= nproc; p1++)
      if( nproc%p1 == 0 )
      {
        int p2 = nproc/p1;
        if( p1 <= p1max && p2 <= p2max )
        {
           // int w1 = p1==1?0:1;
           // int w2 = p2==1?0:1;
           // double f = w2*(double)(ni)/p1 + w1*(double)(nj)/p2;
// try to make each subdomain as square as possible
	  float_sw4 f = fabs((float_sw4)(ni)/p1 - (float_sw4)(nj)/p2);
           if( f < fmin || first )
           {
              fmin = f;
              proc_max[0]   = p1;
              proc_max[1]   = p2;
              first= false;
           }
        }
      }
   return !first;
}

//-----------------------------------------------------------------------
void EW::coarsen1d( int& n, int& ifirst, int& ilast, int periodic )
{
   // n - total number of points 1<=i<=n, 
   // Total index space is 1-ghosts <= i <= n + ghosts
   //
   // This routine coarsens the interval ifirst <= i <= ilast
   // by a factor two, and returns coarsened values of
   // all input parameters.
   //
   int nc;
   if( periodic )
      nc = n/2;
   else
      nc = (n-1)/2+1;
   
   if( ilast == n + m_ghost_points )
      ilast = nc + m_ghost_points;
   else
   {
      ilast = ilast - m_ppadding;
      if( ilast % 2 == 0 )
	 ilast--;
      ilast = (ilast+1)/2;
      ilast = ilast + m_ppadding;
   }
   if( ifirst != 1-m_ghost_points )
   {
      ifirst = ifirst + m_ppadding;
      if( ifirst % 2 == 0 )
	 ifirst++;
      ifirst = (ifirst+1)/2;
      ifirst = ifirst - m_ppadding;
   }
   n = nc;
}

//-----------------------------------------------------------------------
void EW::decomp1d( int nglobal, int myid, int nproc, int& s, int& e )
//
// Decompose index space 1 <= i <= nglobal into nproc blocks
// returns start and end indices for block nr. myid, 
//          where 0 <= myid <= nproc-1
//
{
   int olap    = 2*m_ppadding;
   int nlocal  = (nglobal + (nproc-1)*olap ) / nproc;
   int deficit = (nglobal + (nproc-1)*olap ) % nproc;

   if( myid < deficit )
      s = myid*(nlocal-olap) + myid+1;
   else
      s = myid*(nlocal-olap) + deficit+1;

   if (myid < deficit)
      nlocal = nlocal + 1;

   e = s + nlocal - 1;
}

// -----------------------------
void EW::setup2D_MPICommunications()
{
   if (mVerbose >= 2 && proc_zero())
      cout << "***inside setup2D_MPICommunications***"<< endl;
// Define MPI datatypes for communication across processor boundaries
// For topography: finest grid (curvilinear) only, only one value per grid point (nc=1)
// get the size from the top Cartesian grid
   int g= mNumberOfCartesianGrids-1;
   int ni = m_iEnd[g]-m_iStart[g]+1, nj=m_jEnd[g]-m_jStart[g]+1;
   MPI_Type_vector( nj, m_ppadding,    ni,    m_mpifloat, &m_send_type_2dfinest[0] );
   MPI_Type_vector( 1,  m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type_2dfinest[1] );
   MPI_Type_commit( &m_send_type_2dfinest[0] );
   MPI_Type_commit( &m_send_type_2dfinest[1] );

// Extended number of padding points
   ni = ni + 2*m_ext_ghost_points;
   nj = nj + 2*m_ext_ghost_points;
   int extpadding = m_ppadding + m_ext_ghost_points;
   MPI_Type_vector( nj, extpadding,    ni,    m_mpifloat, &m_send_type_2dfinest_ext[0] );
   MPI_Type_vector( 1,  extpadding*ni, ni*nj, m_mpifloat, &m_send_type_2dfinest_ext[1] );
   MPI_Type_commit( &m_send_type_2dfinest_ext[0] );
   MPI_Type_commit( &m_send_type_2dfinest_ext[1] );

// For mesh refinement: 2D planes with three values per grid point (nc=3)
// Coarser grids
   m_send_type_2dx.resize(mNumberOfCartesianGrids);
   m_send_type_2dy.resize(mNumberOfCartesianGrids);
   m_send_type_2dx1p.resize(mNumberOfCartesianGrids);//padding=1
   m_send_type_2dy1p.resize(mNumberOfCartesianGrids);
   m_send_type_2dx3p.resize(mNumberOfCartesianGrids);//padding=3
   m_send_type_2dy3p.resize(mNumberOfCartesianGrids);
   for( int g = 0 ; g < mNumberOfCartesianGrids ; g++ )
   {
      int ni = m_iEnd[g]-m_iStart[g]+1, nj=m_jEnd[g]-m_jStart[g]+1;
      if( m_croutines )
      {
	 MPI_Type_vector( 3*nj, m_ppadding,    ni,    m_mpifloat, &m_send_type_2dx[g] );
	 MPI_Type_vector( 3,    m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type_2dy[g] );
	 MPI_Type_vector( 3*nj, 1,    ni,    m_mpifloat, &m_send_type_2dx1p[g] );
	 MPI_Type_vector( 3,    ni, ni*nj, m_mpifloat, &m_send_type_2dy1p[g] );
	 MPI_Type_vector( 3*nj, 3,  ni,    m_mpifloat, &m_send_type_2dx3p[g] );
	 MPI_Type_vector( 3,    3*ni, ni*nj, m_mpifloat, &m_send_type_2dy3p[g] );
      }
      else
      {
	 MPI_Type_vector( nj, 3*m_ppadding,    3*ni,    m_mpifloat, &m_send_type_2dx[g] );
	 MPI_Type_vector( 1,  3*m_ppadding*ni, 3*ni*nj, m_mpifloat, &m_send_type_2dy[g] );
	 MPI_Type_vector( nj, 3,    3*ni,    m_mpifloat, &m_send_type_2dx1p[g] );
	 MPI_Type_vector( 1,  3*ni, 3*ni*nj, m_mpifloat, &m_send_type_2dy1p[g] );
	 MPI_Type_vector( nj, 3*3,    3*ni,    m_mpifloat, &m_send_type_2dx3p[g] );
	 MPI_Type_vector( 1,  3*3*ni, 3*ni*nj, m_mpifloat, &m_send_type_2dy3p[g] );
      }
      MPI_Type_commit( &m_send_type_2dx[g] );
      MPI_Type_commit( &m_send_type_2dy[g] );
      MPI_Type_commit( &m_send_type_2dx1p[g] );
      MPI_Type_commit( &m_send_type_2dy1p[g] );
      MPI_Type_commit( &m_send_type_2dx3p[g] );
      MPI_Type_commit( &m_send_type_2dy3p[g] );
   }      
}

// -----------------------------
void EW::setupMPICommunications()
{
   if (mVerbose >= 2 && proc_zero())
      cout << "***inside setupMPICommunications***"<< endl;
// Define MPI datatypes for communication across processor boundaries
   m_send_type1.resize(2*mNumberOfGrids);
   m_send_type3.resize(2*mNumberOfGrids);
   m_send_type4.resize(2*mNumberOfGrids);
   m_send_type21.resize(2*mNumberOfGrids);
   for( int g= 0 ; g < mNumberOfGrids ; g++ )
   {
//      int ni = mU[g].m_ni, nj=mU[g].m_nj, nk=mU[g].m_nk;
      int ni = m_iEnd[g] - m_iStart[g] + 1;
      int nj = m_jEnd[g] - m_jStart[g] + 1;
      int nk = m_kEnd[g] - m_kStart[g] + 1;

      MPI_Type_vector( nj*nk, m_ppadding, ni, m_mpifloat, &m_send_type1[2*g] );
      MPI_Type_vector( nk, m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type1[2*g+1] );

      if( m_croutines )
      {
	 MPI_Type_vector( 3*nj*nk, m_ppadding, ni, m_mpifloat, &m_send_type3[2*g] );
	 MPI_Type_vector( 3*nk, m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type3[2*g+1] );

	 MPI_Type_vector( 4*nj*nk, m_ppadding, ni, m_mpifloat, &m_send_type4[2*g] );
	 MPI_Type_vector( 4*nk, m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type4[2*g+1] );

	 MPI_Type_vector( 21*nj*nk, m_ppadding, ni, m_mpifloat, &m_send_type21[2*g] );
	 MPI_Type_vector( 21*nk, m_ppadding*ni, ni*nj, m_mpifloat, &m_send_type21[2*g+1] );
      }
      else
      {
	 MPI_Type_vector( nj*nk, 3*m_ppadding, 3*ni, m_mpifloat, &m_send_type3[2*g] );
	 MPI_Type_vector( nk, 3*m_ppadding*ni, 3*ni*nj, m_mpifloat, &m_send_type3[2*g+1] );

	 MPI_Type_vector( nj*nk, 4*m_ppadding, 4*ni, m_mpifloat, &m_send_type4[2*g] );
	 MPI_Type_vector( nk, 4*m_ppadding*ni, 4*ni*nj, m_mpifloat, &m_send_type4[2*g+1] );

	 MPI_Type_vector( nj*nk, 21*m_ppadding, 21*ni, m_mpifloat, &m_send_type21[2*g] );
	 MPI_Type_vector( nk, 21*m_ppadding*ni, 21*ni*nj, m_mpifloat, &m_send_type21[2*g+1] );
      }
      MPI_Type_commit( &m_send_type1[2*g] ); 
      MPI_Type_commit( &m_send_type1[2*g+1] ); 

      MPI_Type_commit( &m_send_type3[2*g] ); 
      MPI_Type_commit( &m_send_type3[2*g+1] ); 

      MPI_Type_commit( &m_send_type4[2*g] ); 
      MPI_Type_commit( &m_send_type4[2*g+1] ); 

      MPI_Type_commit( &m_send_type21[2*g] ); 
      MPI_Type_commit( &m_send_type21[2*g+1] ); 
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
//         cout << "Array bounds: " << mU[g].m_ib << " <= i <= " << mU[g].m_ie << " , "
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
//	 MPI_Recv( &dum, 1, MPI_INT, myid-1, tag, m_cartesian_communicator, &status );
//      MPI_Comm_rank( m_cartesian_communicator, &myid );
//      fileout << "Id in 1d " << myid << " Id in 2d " << my_proc_coords[0] << " x " << my_proc_coords[1]
//	      << " has Neighbors x-direction "  << neigh[0] << " and " << neigh[1] << " y-dir " << neigh[2]
//	      << " and " << neigh[3] << endl;
//      if( myid < nprocs-1 )
//	 MPI_Send( &dum, 1, MPI_INT, myid+1, tag, m_cartesian_communicator );
//      fileout.close();
//   }
}


//-----------------------------------------------------------------------
void EW::communicate_array( Sarray& u, int grid )
{
  // REQUIRE2( 0 <= grid && grid < mU.size() , 
  // 	    " Error in communicate_array, grid = " << grid );
   
   REQUIRE2( u.m_nc == 21 || u.m_nc == 4 || u.m_nc == 3 || u.m_nc == 1, "Communicate array, only implemented for one-, three-, four-, and 21-component arrays"
	     << " nc = " << u.m_nc );
   int ie = u.m_ie, ib=u.m_ib, je=u.m_je, jb=u.m_jb, ke=u.m_ke, kb=u.m_kb;
   MPI_Status status;
   if( u.m_nc == 1 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      // X-direction communication
      MPI_Sendrecv( &u(ie-(2*m_ppadding-1),jb,kb), 1, m_send_type1[2*grid], m_neighbor[1], xtag1,
		    &u(ib,jb,kb), 1, m_send_type1[2*grid], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(ib+m_ppadding,jb,kb), 1, m_send_type1[2*grid], m_neighbor[0], xtag2,
		    &u(ie-(m_ppadding-1),jb,kb), 1, m_send_type1[2*grid], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
      // Y-direction communication
      MPI_Sendrecv( &u(ib,je-(2*m_ppadding-1),kb), 1, m_send_type1[2*grid+1], m_neighbor[3], ytag1,
		    &u(ib,jb,kb), 1, m_send_type1[2*grid+1], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(ib,jb+m_ppadding,kb), 1, m_send_type1[2*grid+1], m_neighbor[2], ytag2,
		    &u(ib,je-(m_ppadding-1),kb), 1, m_send_type1[2*grid+1], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   else if( u.m_nc == 3 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      // X-direction communication
      MPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, m_send_type3[2*grid], m_neighbor[1], xtag1,
		    &u(1,ib,jb,kb), 1, m_send_type3[2*grid], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+m_ppadding,jb,kb), 1, m_send_type3[2*grid], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,kb), 1, m_send_type3[2*grid], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
      // Y-direction communication
      MPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),kb), 1, m_send_type3[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, m_send_type3[2*grid+1], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+m_ppadding,kb), 1, m_send_type3[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, m_send_type3[2*grid+1], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   else if( u.m_nc == 4 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      // X-direction communication
      MPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, m_send_type4[2*grid], m_neighbor[1], xtag1,
		    &u(1,ib,jb,kb), 1, m_send_type4[2*grid], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+m_ppadding,jb,kb), 1, m_send_type4[2*grid], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,kb), 1, m_send_type4[2*grid], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
      // Y-direction communication
      MPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),kb), 1, m_send_type4[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, m_send_type4[2*grid+1], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+m_ppadding,kb), 1, m_send_type4[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, m_send_type4[2*grid+1], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   else if( u.m_nc == 21 )
   {
      int xtag1 = 345;
      int xtag2 = 346;
      int ytag1 = 347;
      int ytag2 = 348;
      // X-direction communication
      MPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,kb), 1, m_send_type21[2*grid], m_neighbor[1], xtag1,
		    &u(1,ib,jb,kb), 1, m_send_type21[2*grid], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+m_ppadding,jb,kb), 1, m_send_type21[2*grid], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,kb), 1, m_send_type21[2*grid], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
      // Y-direction communication
      MPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),kb), 1, m_send_type21[2*grid+1], m_neighbor[3], ytag1,
		    &u(1,ib,jb,kb), 1, m_send_type21[2*grid+1], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+m_ppadding,kb), 1, m_send_type21[2*grid+1], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),kb), 1, m_send_type21[2*grid+1], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
}

//-----------------------------------------------------------------------
void EW::communicate_arrays( vector<Sarray>& u )
{
   for( int g= 0 ; g < u.size() ; g++ )
      communicate_array( u[g], g );
}


//-----------------------------------------------------------------------
void EW::communicate_array_2d( Sarray& u, int g, int k )
{
   REQUIRE2( u.m_nc == 3, "Communicate array 2d, only implemented for three-component arrays" );
   REQUIRE2( g < m_send_type_2dx.size(), "Communicate array 2d, only implemented for grid=0.." << m_send_type_2dx.size()-1
	     << " but g= " << g);
   int ie = m_iEnd[g], ib=m_iStart[g];
   int je = m_jEnd[g], jb=m_jStart[g];

   MPI_Status status;
   int xtag1 = 345;
   int xtag2 = 346;
   int ytag1 = 347;
   int ytag2 = 348;

      // X-direction communication
   MPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,k), 1, m_send_type_2dx[g], m_neighbor[1], xtag1,
		 &u(1,ib,jb,k), 1, m_send_type_2dx[g], m_neighbor[0], xtag1,
		 m_cartesian_communicator, &status );
   MPI_Sendrecv( &u(1,ib+m_ppadding,jb,k), 1, m_send_type_2dx[g], m_neighbor[0], xtag2,
		 &u(1,ie-(m_ppadding-1),jb,k), 1, m_send_type_2dx[g], m_neighbor[1], xtag2,
		 m_cartesian_communicator, &status );

      // Y-direction communication
   MPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),k), 1, m_send_type_2dy[g], m_neighbor[3], ytag1,
		 &u(1,ib,jb,k), 1, m_send_type_2dy[g], m_neighbor[2], ytag1,
		 m_cartesian_communicator, &status );
   MPI_Sendrecv( &u(1,ib,jb+m_ppadding,k), 1, m_send_type_2dy[g], m_neighbor[2], ytag2,
		 &u(1,ib,je-(m_ppadding-1),k), 1, m_send_type_2dy[g], m_neighbor[3], ytag2,
		 m_cartesian_communicator, &status );
}

//-----------------------------------------------------------------------
void EW::communicate_array_2d_ext( Sarray& u )
{
   REQUIRE2( u.m_nc == 1, "Communicate array 2d ext, only implemented for one-component arrays" );
   int g = mNumberOfGrids-1;
   int ie = m_iEnd[g]+m_ext_ghost_points, ib=m_iStart[g]-m_ext_ghost_points;
   int je = m_jEnd[g]+m_ext_ghost_points, jb=m_jStart[g]-m_ext_ghost_points;

   MPI_Status status;
   int xtag1 = 345;
   int xtag2 = 346;
   int ytag1 = 347;
   int ytag2 = 348;
   int k=1;
   int extpadding = m_ppadding+m_ext_ghost_points;
      // X-direction communication
   MPI_Sendrecv( &u(1,ie-(2*extpadding-1),jb,k), 1, m_send_type_2dfinest_ext[0], m_neighbor[1], xtag1,
		 &u(1,ib,jb,k), 1, m_send_type_2dfinest_ext[0], m_neighbor[0], xtag1,
		 m_cartesian_communicator, &status );
   MPI_Sendrecv( &u(1,ib+extpadding,jb,k), 1, m_send_type_2dfinest_ext[0], m_neighbor[0], xtag2,
		 &u(1,ie-(extpadding-1),jb,k), 1, m_send_type_2dfinest_ext[0], m_neighbor[1], xtag2,
		 m_cartesian_communicator, &status );

      // Y-direction communication
   MPI_Sendrecv( &u(1,ib,je-(2*extpadding-1),k), 1, m_send_type_2dfinest_ext[1], m_neighbor[3], ytag1,
		 &u(1,ib,jb,k), 1, m_send_type_2dfinest_ext[1], m_neighbor[2], ytag1,
		 m_cartesian_communicator, &status );
   MPI_Sendrecv( &u(1,ib,jb+extpadding,k), 1, m_send_type_2dfinest_ext[1], m_neighbor[2], ytag2,
		 &u(1,ib,je-(extpadding-1),k), 1, m_send_type_2dfinest_ext[1], m_neighbor[3], ytag2,
		 m_cartesian_communicator, &status );
}

//-----------------------------------------------------------------------
void EW::communicate_array_2d_asym( Sarray& u, int g, int k )
{
   REQUIRE2( u.m_nc == 3, "Communicate array 2d asym, only implemented for three-component arrays" );
   REQUIRE2( g < m_send_type_2dx3p.size(), "Communicate array 2d asym, only implemented for grid=0.." 
	     << m_send_type_2dx3p.size()-1 << " but g= " << g);
   int ie = m_iEnd[g], ib=m_iStart[g];
   int je = m_jEnd[g], jb=m_jStart[g];

   MPI_Status status;
   int xtag1 = 345;
   int xtag2 = 346;
   int ytag1 = 347;
   int ytag2 = 348;
   int extupper = ie % 2;
   int extlower = 1-(ib % 2);

      // X-direction communication
   if( extupper && extlower )
   {
      MPI_Sendrecv( &u(1,ie-(3),jb,k), 1, m_send_type_2dx3p[g], m_neighbor[1], xtag1,
		    &u(1,ib,jb,k), 1, m_send_type_2dx3p[g], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+3,jb,k), 1, m_send_type_2dx1p[g], m_neighbor[0], xtag2,
		    &u(1,ie,jb,k), 1, m_send_type_2dx1p[g], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
   }
   else if( extupper )
   {
      MPI_Sendrecv( &u(1,ie-3,jb,k), 1, m_send_type_2dx3p[g], m_neighbor[1], xtag1,
		    &u(1,ib,jb,k), 1, m_send_type_2dx[g], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+m_ppadding,jb,k), 1, m_send_type_2dx[g], m_neighbor[0], xtag2,
		    &u(1,ie,jb,k), 1, m_send_type_2dx1p[g], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
   }
   else if( extlower )
   {
      MPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,k), 1, m_send_type_2dx[g], m_neighbor[1], xtag1,
		    &u(1,ib,jb,k), 1, m_send_type_2dx3p[g], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+3,jb,k), 1, m_send_type_2dx1p[g], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,k), 1, m_send_type_2dx[g], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
   }
   else
   {
      MPI_Sendrecv( &u(1,ie-(2*m_ppadding-1),jb,k), 1, m_send_type_2dx[g], m_neighbor[1], xtag1,
		    &u(1,ib,jb,k), 1, m_send_type_2dx[g], m_neighbor[0], xtag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib+m_ppadding,jb,k), 1, m_send_type_2dx[g], m_neighbor[0], xtag2,
		    &u(1,ie-(m_ppadding-1),jb,k), 1, m_send_type_2dx[g], m_neighbor[1], xtag2,
		    m_cartesian_communicator, &status );
   }

   extupper = je % 2;
   extlower = 1-(jb % 2);
      // Y-direction communication

   if( extupper && extlower )
   {
      MPI_Sendrecv( &u(1,ib,je-(3),k), 1, m_send_type_2dy3p[g], m_neighbor[3], ytag1,
		    &u(1,ib,jb,k), 1, m_send_type_2dy3p[g], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+3,k), 1, m_send_type_2dy1p[g], m_neighbor[2], ytag2,
		    &u(1,ib,je,k), 1, m_send_type_2dy1p[g], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   else if( extupper )
   {
      MPI_Sendrecv( &u(1,ib,je-3,k), 1, m_send_type_2dy3p[g], m_neighbor[3], ytag1,
		    &u(1,ib,jb,k), 1, m_send_type_2dy[g], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+m_ppadding,k), 1, m_send_type_2dy[g], m_neighbor[2], ytag2,
		    &u(1,ib,je,k), 1, m_send_type_2dy1p[g], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   else if( extlower )
   {
      MPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),k), 1, m_send_type_2dy[g], m_neighbor[3], ytag1,
		    &u(1,ib,jb,k), 1, m_send_type_2dy3p[g], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+3,k), 1, m_send_type_2dy1p[g], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),k), 1, m_send_type_2dy[g], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
   else
   {
      MPI_Sendrecv( &u(1,ib,je-(2*m_ppadding-1),k), 1, m_send_type_2dy[g], m_neighbor[3], ytag1,
		    &u(1,ib,jb,k), 1, m_send_type_2dy[g], m_neighbor[2], ytag1,
		    m_cartesian_communicator, &status );
      MPI_Sendrecv( &u(1,ib,jb+m_ppadding,k), 1, m_send_type_2dy[g], m_neighbor[2], ytag2,
		    &u(1,ib,je-(m_ppadding-1),k), 1, m_send_type_2dy[g], m_neighbor[3], ytag2,
		    m_cartesian_communicator, &status );
   }
}

//-----------------------------------------------------------------------
void EW::communicate_array_2dfinest( Sarray& u )
{
   REQUIRE2( u.m_nc == 1, "Communicate array 2dfinest, only implemented for one-component arrays" );

   int ie = u.m_ie, ib=u.m_ib, je=u.m_je, jb=u.m_jb;
   REQUIRE2( ib == m_iStart[mNumberOfGrids-1] && ie == m_iEnd[mNumberOfGrids-1] &&
             jb == m_jStart[mNumberOfGrids-1] && je == m_jEnd[mNumberOfGrids-1] , 
             "Communicate array 2d: Can only use it on the finest grid, grid sizes don't match");

   MPI_Status status;
   int xtag1 = 345;
   int xtag2 = 346;
   int ytag1 = 347;
   int ytag2 = 348;
      // X-direction communication
   MPI_Sendrecv( &u(ie-(2*m_ppadding-1),jb,1), 1, m_send_type_2dfinest[0], m_neighbor[1], xtag1,
		 &u(ib,jb,1), 1, m_send_type_2dfinest[0], m_neighbor[0], xtag1,
		 m_cartesian_communicator, &status );
   MPI_Sendrecv( &u(ib+m_ppadding,jb,1), 1, m_send_type_2dfinest[0], m_neighbor[0], xtag2,
		 &u(ie-(m_ppadding-1),jb,1), 1, m_send_type_2dfinest[0], m_neighbor[1], xtag2,
		 m_cartesian_communicator, &status );
      // Y-direction communication
   MPI_Sendrecv( &u(ib,je-(2*m_ppadding-1),1), 1, m_send_type_2dfinest[1], m_neighbor[3], ytag1,
		 &u(ib,jb,1), 1, m_send_type_2dfinest[1], m_neighbor[2], ytag1,
		 m_cartesian_communicator, &status );
   MPI_Sendrecv( &u(ib,jb+m_ppadding,1), 1, m_send_type_2dfinest[1], m_neighbor[2], ytag2,
		 &u(ib,je-(m_ppadding-1),1), 1, m_send_type_2dfinest[1], m_neighbor[3], ytag2,
		 m_cartesian_communicator, &status );
}
