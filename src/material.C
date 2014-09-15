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
#include "EW.h"

//-----------------------------------------------------------------------
void EW::convert_material_to_mulambda( )
{
  for( int g = 0 ; g < mNumberOfGrids; g++)
  {
      // On input, we have stored cs in MU, cp in Lambda
      // use mu = rho*cs*cs and lambda = rho*cp*cp  - 2*mu
      
      for( int k = m_kStart[g] ; k <= m_kEnd[g]; k++ )
      {
          for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
            {
              for( int i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
                {
//                   if (k==32)
//                     printf("%i %i %i %i %f %f %f %f %f\n",i,j,k,g,mRho[g](i,j,k)*mMu[g](i,j,k)*mMu[g](i,j,k),mRho[g](i,j,k)*mLambda[g](i,j,k)*mLambda[g](i,j,k)-2*mMu[g](i,j,k),mMu[g](i,j,k),mLambda[g](i,j,k),mRho[g](i,j,k));

                  mMu[g](i,j,k)     = mRho[g](i,j,k)*mMu[g](i,j,k)*mMu[g](i,j,k);                  
                  mLambda[g](i,j,k) = mRho[g](i,j,k)*mLambda[g](i,j,k)*mLambda[g](i,j,k)-2*mMu[g](i,j,k);
                }
            }
      }
    } // end for all grids
} // end convert_material_to_mulambda

//-----------------------------------------------------------------------
void EW::check_materials()
{

  //---------------------------------------------------------------
  // Verify that the density is nonzero and positive in the 
  // internal grid points
  //---------------------------------------------------------------
  
  double mins[8],maxs[8];

// confusing with variables and functions that have names which only differ in capitalization 
  double lmin = localMin(mRho);
  MPI_Allreduce(&lmin,&mins[0],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  lmin = localMinVp();  
  MPI_Allreduce(&lmin,&mins[1],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  lmin = localMinVs();  
  MPI_Allreduce(&lmin,&mins[2],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  lmin = localMin(mMu);
  MPI_Allreduce(&lmin,&mins[3],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  lmin = localMin(mLambda);
  MPI_Allreduce(&lmin,&mins[4],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  
  CHECK_INPUT(mins[2] >= 0.0,
          "Error: the material data has s velocities that are negative.");

  lmin = localMinVpOverVs();  
  MPI_Allreduce(&lmin,&mins[5],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);

  double lmax = localMax(mRho);
  MPI_Allreduce(&lmax,&maxs[0],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVp();  
  MPI_Allreduce(&lmax,&maxs[1],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVs();  
  MPI_Allreduce(&lmax,&maxs[2],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  lmax = localMax(mMu);
  MPI_Allreduce(&lmax,&maxs[3],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  lmax = localMax(mLambda);
  MPI_Allreduce(&lmax,&maxs[4],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVpOverVs();  
  MPI_Allreduce(&lmax,&maxs[5],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);

  if( usingAttenuation() )
  {
      lmin = localMin(mQs);
      MPI_Allreduce(&lmin,&mins[6],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
      lmin = localMin(mQp);
      MPI_Allreduce(&lmin,&mins[7],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
      lmax = localMax(mQs);
      MPI_Allreduce(&lmax,&maxs[6],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
      lmax = localMax(mQp);
      MPI_Allreduce(&lmax,&maxs[7],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  }
  
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if( myRank == 0 )
  {
    string indent  = "\n       ";
    string indents =   "       ";
      
    cout << indent << "----------- Material properties ranges ---------------"
//            << indent << "  For grid [" << m_level << "][" << m_gridnr << "]" << endl
	 << indent << mins[0] << " kg/m^3 <=  Density <= " << maxs[0] << " kg/m^3"
	 << indent << mins[1] << " m/s    <=  Vp      <= " << maxs[1] << " m/s"
	 << indent << mins[2] << " m/s    <=  Vs      <= " << maxs[2] << " m/s" 
	 << indent << mins[5] << "        <=  Vp/Vs   <= " << maxs[5]
	 << indent << mins[3] << " Pa     <=  mu      <= " << maxs[3] << " Pa"
	 << indent << mins[4] << " Pa     <=  lambda  <= " << maxs[4] << " Pa" << endl;

    if( usingAttenuation() )
    {
      cout << indents << "Using attenuation "
           << indent << mins[6] << "        <=  Qs      <= " << maxs[6] << "  "
           << indent << mins[7] << "        <=  Qp      <= " << maxs[7] << "  " << endl;
    }
    cout  << indents << "------------------------------------------------------" << endl;
  }
   
  VERIFY2(mins[5] >= sqrt(2.),
	  "Error: vp/vs is smaller than sqrt(2).");

  if( mins[0] <= 0.0 )
  {
    for (int g = 0; g < mNumberOfGrids; g++)
      for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	  for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	  {
	    CHECK_INPUT( mRho[g](i,j,k) > 0., "Density= " << mRho[g](i,j,k)<< " in grid g= " << g << " at point " 
			  << " (" << i <<","<<j<<","<<k<<") ");
	  }
  }
   
  if( mins[3] < 0.0 )
  {
    for (int g = 0; g < mNumberOfGrids; g++)
      for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	  for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	  {
	    CHECK_INPUT( mMu[g](i,j,k) >= 0., "mu= " << mMu[g](i,j,k)<< " in grid g= " << g << " at point " 
			  << " (" << i <<","<<j<<","<<k<<") ");
	  }
  }
  if( mins[4] <= 0.0 )
  {
    for (int g = 0; g < mNumberOfGrids; g++)
      for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	  for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	  {
	    CHECK_INPUT( mLambda[g](i,j,k) >= 0., "lambda= " << mLambda[g](i,j,k)<< " in grid g= " << g << " at point " 
			 << " (" << i <<","<<j<<","<<k<<") ");
	  }
  }

   CHECK_INPUT(mins[0] > 0.0,
		"Error: the material data has density values less than or equal to zero.");
   CHECK_INPUT(mins[1] > 0.0,
           "Error: the material data has p velocities that are less than or equal to zero.");
   CHECK_INPUT(mins[3] >= 0.0,
           "Error: mu has values that are negative.");
   CHECK_INPUT(mins[4] > 0.0,
           "Error: lambda has values that are negative or zero.");

// check material ranges on each grid
   if (mVerbose >= 3)
   {
   double minRho, maxRho, minMu, maxMu, minLambda, maxLambda;
   
   for( int g = 0 ; g < mNumberOfGrids; g++)
   {
     minRho=1e100, minMu=1e100, minLambda=1e100;
     maxRho=0, maxMu=0, maxLambda=0;
     
     for( int k = m_kStart[g] ; k <= m_kEnd[g]; k++ )
       for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	 for( int i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	 {
	   if (mMu[g](i,j,k) < minMu) minMu=mMu[g](i,j,k);
	   if (mMu[g](i,j,k) > maxMu) maxMu=mMu[g](i,j,k);
	   if (mLambda[g](i,j,k) < minLambda) minLambda=mLambda[g](i,j,k);
	   if (mLambda[g](i,j,k) > maxLambda) maxLambda=mLambda[g](i,j,k);
	   if (mRho[g](i,j,k) < minRho) minRho=mRho[g](i,j,k);
	   if (mRho[g](i,j,k) > maxRho) maxRho=mRho[g](i,j,k);
	 }
// communicate min & max
     MPI_Allreduce(&minRho,&mins[0],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxRho,&maxs[0],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
     MPI_Allreduce(&minMu,&mins[1],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxMu,&maxs[1],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
     MPI_Allreduce(&minLambda,&mins[2],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxLambda,&maxs[2],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
// printout results
     if (proc_zero())
     {
       printf("Grid #%i:, %e <= Rho <= %e, %e <= Mu <= %e, %e <= Lambda <= %e\n", g, mins[0], maxs[0], mins[1], maxs[1], mins[2], maxs[2]);
     }
     
   } // end for all grids
   }
// end mVerbose >= 3   
}

//-----------------------------------------------------------------------
double EW::localMin(std::vector<Sarray> & a_field) 
{
  double lmin = a_field[0](m_iStart[0],m_jStart[0],m_kStart[0]);
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (a_field[g](i,j,k) < lmin)
                    {
                      lmin = a_field[g](i,j,k);
                    }
                }
            }
        }
    }

  return lmin; 
}

//-----------------------------------------------------------------------
double EW::localMax(std::vector<Sarray> & a_field) 
{
  double lmax = a_field[0](m_iStart[0],m_jStart[0],m_kStart[0]);
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (a_field[g](i,j,k) > lmax)
                    {
                      lmax = a_field[g](i,j,k);
                    }
                }
            }
        }
    }

  return lmax; 
}

//-----------------------------------------------------------------------
double EW::localMinVp() 
{
  double lmin = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k)) < lmin)
                    {
                      lmin = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return lmin; 
}

//-----------------------------------------------------------------------
double EW::localMaxVp() 
{
  double lmax = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k)) > lmax)
                    {
                      lmax = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return lmax; 
}

//-----------------------------------------------------------------------
double EW::localMinVs() 
{
  double lmin = sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])/mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) < lmin)
                    {
                      lmin = sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return lmin; 
}

//-----------------------------------------------------------------------
double EW::localMaxVs() 
{
  double lmax = sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])/mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if ( sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) > lmax)
                    {
                      lmax = sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return lmax; 
}

//-----------------------------------------------------------------------
double EW::localMinVpOverVs() 
{
  double lmin = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]))/sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])
									     /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
// Unneccessary to divided by rho in Vp and Vs because it cancels
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) < lmin)
                    {
                      lmin = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return lmin; 
}

//-----------------------------------------------------------------------
double EW::localMaxVpOverVs() 
{
  double lmax = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]))/sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])
									     /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
// Unneccessary to divided by rho in Vp and Vs because it cancels
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) > lmax)
                    {
                      lmax = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return lmax; 
}

// the same routine is defined in EtreeFile.C, but in the class EtreeFile
//-----------------------------------------------------------------------
//void EW::extrapolateInZ(Sarray& field, bool useThreshold, double thresHoldValue, bool linear) 
//{
//  // -------------------------------------------------------------
//  // We use linear extrapolation to update boundary points
//  //
//  // field(n+1) = 2*field(n) - field(n-1)
//  //
//  // Unless, 2*field(n) <= field(n-1), then we just set it to
//  // field(n).
//  // -------------------------------------------------------------
//   //  int numGhostPoints = getExternalGhostPointsPerBoundaryPoint();
//
//  int i, j, k;
//  double extField;

// tmp
//  if (proc_zero())
//  {
//    printf("extrapolateInZ: m_kb=%i, m_ghost_points=%i\n", field.m_kb, m_ghost_points);
//  }
  
// only extrapolate on the "low-k" side  
//  for( j = field.m_jb; j <= field.m_je; j++ )
//    for( i = field.m_ib; i <= field.m_ie ; i++ )
//      for( k = field.m_kb + m_ghost_points-1 ; k >= field.m_kb; k-- )
//      {
//	if (linear && 2.*field(i,j,k+1) > field(i,j,k+2)) // check if linear extrapolation will lead to a positive value
//	{
//	  extField = 2.*field(i,j,k+1)-field(i,j,k+2);
//	}
//	else // constant extrapolation
//	{
//	  extField = field(i,j,k+1);
//	}
//	if (useThreshold && extField<thresHoldValue)
//	  extField = thresHoldValue;
//	
//	field(i,j,k) = extField;
//      } // end for k...
//} // end extrapolateInZ

//-----------------------------------------------------------------------
void EW::extrapolateInXY( vector<Sarray>& field )
{
   for( int g= 0; g < mNumberOfGrids ; g++ )
   {
      if( m_iStartInt[g] == 1 )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iStart[g] ; i < 1 ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](1,j,k);
	       }
      if( m_iEndInt[g] == m_global_nx[g] )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iEndInt[g]+1 ; i <= m_iEnd[g] ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](m_iEndInt[g],j,k);
	       }
      if( m_jStartInt[g] == 1 )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j < 1 ; j++ )
	       for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](i,1,k);
	       }
      if( m_jEndInt[g] == m_global_ny[g] )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jEndInt[g]+1 ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](i,m_jEndInt[g],k);
	       }
// corners not necessary to treat explicitly???
      
   }
}

//-----------------------------------------------------------------------
void EW::extrapolateInZ( int g, Sarray& field, bool lowk, bool highk )
{
   if( lowk )
      for( int k=m_kStart[g] ; k < 1 ; k++ )
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       if( field(i,j,k) == -1 )
		  field(i,j,k) = field(i,j,1);
   if( highk )
      for( int k=m_kEndInt[g]+1 ; k <= m_kEnd[g] ; k++ )
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       if( field(i,j,k) == -1 )
		  field(i,j,k) = field(i,j,m_kEndInt[g]);
}

//-----------------------------------------------------------------------
void EW::extrapolateInXYvector( vector<Sarray>& field )
{

   for( int g= 0; g < mNumberOfGrids ; g++ )
   {
      int nc = field[g].m_nc;
      if( m_iStartInt[g] == 1 )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iStart[g] ; i < 1 ; i++ )
		  for( int m=1 ; m <= nc ; m++ )
		  {
		     if( field[g](m,i,j,k) == -1 )
			field[g](m,i,j,k) = field[g](m,1,j,k);
		  }
      if( m_iEndInt[g] == m_global_nx[g] )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iEndInt[g]+1 ; i <= m_iEnd[g] ; i++ )
		  for( int m=1 ; m <= nc ; m++ )
		  {
		     if( field[g](m,i,j,k) == -1 )
			field[g](m,i,j,k) = field[g](m,m_iEndInt[g],j,k);
		  }
      if( m_jStartInt[g] == 1 )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j < 1 ; j++ )
	       for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
		  for( int m=1 ; m <= nc ; m++ )
		  {
		     if( field[g](m,i,j,k) == -1 )
			field[g](m,i,j,k) = field[g](m,i,1,k);
		  }
      if( m_jEndInt[g] == m_global_ny[g] )
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jEndInt[g]+1 ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
		  for( int m=1 ; m <= nc ; m++ )
		  {
		     if( field[g](m,i,j,k) == -1 )
			field[g](m,i,j,k) = field[g](m,i,m_jEndInt[g],k);
		  }
// corners not necessary to treat explicitly???
      
   }
}

//-----------------------------------------------------------------------
void EW::extrapolateInZvector( int g, Sarray& field, bool lowk, bool highk )
{
   int nc = field.m_nc;
   if( lowk )
      for( int k=m_kStart[g] ; k < 1 ; k++ )
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       for( int m=1 ; m <= nc ; m++ )
	       {
		  if( field(m,i,j,k) == -1 )
		     field(m,i,j,k) = field(m,i,j,1);
	       }
   if( highk )
      for( int k=m_kEndInt[g]+1 ; k <= m_kEnd[g] ; k++ )
	 for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       for( int m=1 ; m <= nc ; m++ )
	       {
		  if( field(m,i,j,k) == -1 )
		     field(m,i,j,k) = field(m,i,j,m_kEndInt[g]);
	       }
}

