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
      
#pragma omp parallel for     
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


   // Minimum allowed  cp/cs, positive definite operator requires cp/cs > sqrt(4/3) = 1.155...
   //   lambda >0 requires cp/cs > sqrt(2)
  const float_sw4 mincpcsratio = 1.15;  // 1.2
  const float_sw4 la_min_fact = mincpcsratio*mincpcsratio-2;
  
  float_sw4 mins[8],maxs[8];

  float_sw4 lmin = localMin(mRho);
  MPI_Allreduce(&lmin,&mins[0],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  lmin = localMinVp();  
  MPI_Allreduce(&lmin,&mins[1],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  lmin = localMinVs();  
  MPI_Allreduce(&lmin,&mins[2],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  lmin = localMin(mMu);
  MPI_Allreduce(&lmin,&mins[3],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  lmin = localMin(mLambda);
  MPI_Allreduce(&lmin,&mins[4],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  
  CHECK_INPUT(mins[2] >= 0.0,
          "Error: the material data has s velocities that are negative.");

  lmin = localMinVpOverVs();  
  MPI_Allreduce(&lmin,&mins[5],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);

  float_sw4 lmax = localMax(mRho);
  MPI_Allreduce(&lmax,&maxs[0],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVp();  
  MPI_Allreduce(&lmax,&maxs[1],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVs();  
  MPI_Allreduce(&lmax,&maxs[2],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMax(mMu);
  MPI_Allreduce(&lmax,&maxs[3],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMax(mLambda);
  MPI_Allreduce(&lmax,&maxs[4],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVpOverVs();  
  MPI_Allreduce(&lmax,&maxs[5],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);

  if( usingAttenuation() && !m_twilight_forcing)
  {
      lmin = localMin(mQs);
      MPI_Allreduce(&lmin,&mins[6],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
      lmin = localMin(mQp);
      MPI_Allreduce(&lmin,&mins[7],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
      lmax = localMax(mQs);
      MPI_Allreduce(&lmax,&maxs[6],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
      lmax = localMax(mQp);
      MPI_Allreduce(&lmax,&maxs[7],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
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

    if( usingAttenuation() && !m_twilight_forcing)
    {
      cout << indents << "Using attenuation "
           << indent << mins[6] << "        <=  Qs      <= " << maxs[6] << "  "
           << indent << mins[7] << "        <=  Qp      <= " << maxs[7] << "  " << endl;
    }
    cout  << indents << "------------------------------------------------------" << endl;
  }
   

  if( mins[0] <= 0.0 )
  {
    for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
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
#pragma omp parallel for
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
#pragma omp parallel for
      for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	  for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	  {
	     CHECK_INPUT( mLambda[g](i,j,k) >= la_min_fact*mMu[g](i,j,k), "lambda= " << mLambda[g](i,j,k)<< " in grid g= " << g << " at point " 
			 << " (" << i <<","<<j<<","<<k<<") ");
	  }
  }
  if( m_use_attenuation && !m_twilight_forcing)
  {
     if( mins[6] <= 0.0 )
     {
	for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
	   for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	      for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
		 for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
		 {
		    CHECK_INPUT( mQs[g](i,j,k) >= 0., "Qs= " << mQs[g](i,j,k)<< " in grid g= " << g << " at point " 
				 << " (" << i <<","<<j<<","<<k<<") ");
		 }
     }
     if( mins[7] <= 0.0 )
     {
	for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
	   for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	      for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
		 for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
		 {
		    CHECK_INPUT( mQp[g](i,j,k) >= 0., "Qp= " << mQp[g](i,j,k)<< " in grid g= " << g << " at point " 
				 << " (" << i <<","<<j<<","<<k<<") ");
		 }
     }
  }

   CHECK_INPUT(mins[0] > 0.0,
		"Error: the material data has density values less than or equal to zero.");
   CHECK_INPUT(mins[1] > 0.0,
           "Error: the material data has p velocities that are less than or equal to zero.");
   CHECK_INPUT(mins[3] >= 0.0,
           "Error: mu has values that are negative.");
   //   CHECK_INPUT(mins[4] > 0.0,
   //           "Error: lambda has values that are negative or zero.");
  
  VERIFY2(mins[5] >= mincpcsratio,
	  "Error: vp/vs is smaller than set limit, " << mincpcsratio );

// check material ranges on each grid
   if (mVerbose >= 3)
   {
   float_sw4 minRho, maxRho, minMu, maxMu, minLambda, maxLambda;
   
   for( int g = 0 ; g < mNumberOfGrids; g++)
   {
     minRho=1e38, minMu=1e38, minLambda=1e38;
     maxRho=0, maxMu=0, maxLambda=0;
     
#pragma omp parallel for reduction(min:minMu,minLambda,minRho) reduction(max:maxMu,maxLambda,maxRho)
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
     MPI_Allreduce(&minRho,&mins[0],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxRho,&maxs[0],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
     MPI_Allreduce(&minMu,&mins[1],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxMu,&maxs[1],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
     MPI_Allreduce(&minLambda,&mins[2],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxLambda,&maxs[2],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
// printout results
     if (proc_zero())
     {
       printf("Grid #%i:, %e <= Rho <= %e, %e <= Mu <= %e, %e <= Lambda <= %e\n", g, mins[0], maxs[0], mins[1], maxs[1], mins[2], maxs[2]);
     }
     
   } // end for all grids
   }
// end mVerbose >= 3 

// evaluate min(Cs) and max(sqrt(Cp^2+2*Cs^2))for each grid  
   if (mVerbose >= 1)
   {
     double Cs, C_hat, minCs, maxCs, minC_hat, maxC_hat;
   
     for( int g = 0 ; g < mNumberOfGrids; g++)
     {
       minCs=1e100, minC_hat=1e100;
       maxCs=0, maxC_hat=0;
       
       for( int k = m_kStart[g] ; k <= m_kEnd[g]; k++ )
	 for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
	   for( int i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	   {
// take square root after computing the min and max
	     Cs = mMu[g](i,j,k)/mRho[g](i,j,k);
	     C_hat = (mLambda[g](i,j,k) + 4*mMu[g](i,j,k))/mRho[g](i,j,k);

	     
	     if (Cs < minCs) minCs=Cs;
	     if (Cs > maxCs) maxCs=Cs;

	     if (C_hat < minC_hat) minC_hat=C_hat;
	     if (C_hat > maxC_hat) maxC_hat=C_hat;
	   }
       minCs = sqrt(minCs);
       maxCs = sqrt(maxCs);
       minC_hat = sqrt(minC_hat);
       maxC_hat = sqrt(maxC_hat);
       
// communicate min & max
       MPI_Allreduce(&minCs,&mins[0],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
       MPI_Allreduce(&maxCs,&maxs[0],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
       MPI_Allreduce(&minC_hat,&mins[1],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
       MPI_Allreduce(&maxC_hat,&maxs[1],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
// printout results
       if (mVerbose >=2 && proc_zero())
       {
	 printf("Material model info, Grid g=%i: %e <= Cs <= %e, %e <= C-hat <= %e, h[g]/max(C-hat) = %e\n",
		g, mins[0], maxs[0], mins[1], maxs[1], mGridSize[g]/maxs[1]);
       }
     } // end for all grids
   }
// end mVerbose >= 1
}

//-----------------------------------------------------------------------
float_sw4 EW::localMin(std::vector<Sarray> & a_field) 
{
  float_sw4 lmin_all = a_field[0](m_iStart[0],m_jStart[0],m_kStart[0]);
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmin=1e38;
#pragma omp parallel for reduction(min:lmin)
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (a_field[g](i,j,k) < lmin)
                    {
                      lmin = a_field[g](i,j,k);
		      //		      cout << "lmin = " << lmin << " at " << i << " " << j << " " << k << endl;
                    }
                }
            }
        }
      lmin_all = lmin < lmin_all?lmin:lmin_all;
    }

  return lmin_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMax(std::vector<Sarray> & a_field) 
{
  float_sw4 lmax_all = a_field[0](m_iStart[0],m_jStart[0],m_kStart[0]);
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmax=-1e38;
#pragma omp parallel for reduction(max:lmax)
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
      lmax_all = lmax > lmax_all?lmax:lmax_all;
    }

  return lmax_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMinVp() 
{
  float_sw4 lmin_all = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmin=1e38;
#pragma omp parallel for reduction(min:lmin)
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
      lmin_all = lmin < lmin_all?lmin:lmin_all;
    }
  return lmin_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMaxVp() 
{
  float_sw4 lmax_all = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmax=-1e38;
#pragma omp parallel for reduction(max:lmax)
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
      lmax_all = lmax > lmax_all?lmax:lmax_all;
    }

  return lmax_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMinVs() 
{
  float_sw4 lmin_all = sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])/mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmin=1e38;
#pragma omp parallel for reduction(min:lmin)
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
      lmin_all = lmin < lmin_all?lmin:lmin_all;
    }

  return lmin_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMaxVs() 
{
  float_sw4 lmax_all = sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])/mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmax=-1e38;
#pragma omp parallel for reduction(max:lmax)
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
      lmax_all = lmax > lmax_all?lmax:lmax_all;
    }
  return lmax_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMinVpOverVs() 
{
  float_sw4 lmin_all = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]))/sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])
									     /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmin=1e38;
#pragma omp parallel for reduction(min:lmin)
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
      lmin_all = lmin < lmin_all?lmin:lmin_all;
    }
  return lmin_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMaxVpOverVs() 
{
  float_sw4 lmax_all = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]))/sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])
									     /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmax=-1e38;
#pragma omp parallel for reduction(max:lmax)
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
      lmax_all = lmax > lmax_all?lmax:lmax_all;
    }
  return lmax_all; 
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
#pragma omp parallel for
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iStart[g] ; i < 1 ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](1,j,k);
	       }
      if( m_iEndInt[g] == m_global_nx[g] )
#pragma omp parallel for
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iEndInt[g]+1 ; i <= m_iEnd[g] ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](m_iEndInt[g],j,k);
	       }
      if( m_jStartInt[g] == 1 )
#pragma omp parallel for
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j < 1 ; j++ )
	       for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
	       {
		  if( field[g](i,j,k) == -1 )
		     field[g](i,j,k) = field[g](i,1,k);
	       }
      if( m_jEndInt[g] == m_global_ny[g] )
#pragma omp parallel for
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
#pragma omp parallel for
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iStart[g] ; i < 1 ; i++ )
		  for( int m=1 ; m <= nc ; m++ )
		  {
		     if( field[g](m,i,j,k) == -1 )
			field[g](m,i,j,k) = field[g](m,1,j,k);
		  }
      if( m_iEndInt[g] == m_global_nx[g] )
#pragma omp parallel for
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
	       for( int i=m_iEndInt[g]+1 ; i <= m_iEnd[g] ; i++ )
		  for( int m=1 ; m <= nc ; m++ )
		  {
		     if( field[g](m,i,j,k) == -1 )
			field[g](m,i,j,k) = field[g](m,m_iEndInt[g],j,k);
		  }
      if( m_jStartInt[g] == 1 )
#pragma omp parallel for
         for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
	    for( int j=m_jStart[g] ; j < 1 ; j++ )
	       for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
		  for( int m=1 ; m <= nc ; m++ )
		  {
		     if( field[g](m,i,j,k) == -1 )
			field[g](m,i,j,k) = field[g](m,i,1,k);
		  }
      if( m_jEndInt[g] == m_global_ny[g] )
#pragma omp parallel for
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

//--------- Material properties for MR ---------------
void EW::setup_MR_coefficients()
{
// stretching on the fine side
#define str_x(i) m_sg_str_x[g][(i-m_iStart[g])]   
#define str_y(j) m_sg_str_y[g][(j-m_jStart[g])]   
// calculate m_Morc, etc for all Cartesian grids
   if( !m_anisotropic )
   {
      for(int g = 0; g < mNumberOfCartesianGrids; g++ )
      {
         int nk = m_global_nz[g];
         for (int c=1; c<=3; c++)
            for (int j=m_jStart[g]; j<=m_jEnd[g]; j++)
               for (int i=m_iStart[g]; i<=m_iEnd[g]; i++)
               {
                  float_sw4 irho=1/mRho[g](i,j,1);
                  m_Morc[g](i,j,1) = mMu[g](i,j,1)*irho; // mu/rho at k=1
                  m_Mlrc[g](i,j,1) = (2*mMu[g](i,j,1)+mLambda[g](i,j,1))*irho; // (2*mu+lambda)/rho at k=1
                  m_Mucs[g](i,j,1) = mMu[g](i,j,1)/(str_x(i)*str_y(j)); // mu/str at 1
                  m_Mlcs[g](i,j,1) = (2*mMu[g](i,j,1)+mLambda[g](i,j,1))/(str_x(i)*str_y(j)); //(2*mu + lambda)/str at 1

                  float_sw4 irhoN=1/mRho[g](i,j,nk);
                  m_Morf[g](i,j, nk) = mMu[g](i,j, nk)*irhoN; // mu/rho at nk
                  m_Mlrf[g](i,j, nk) = (2*mMu[g](i,j, nk)+mLambda[g](i,j, nk))*irhoN; // (2*mu+lambda)/rho at nk
                  m_Mufs[g](i,j,nk) = mMu[g](i,j,nk)/(str_x(i)*str_y(j)); // mu/str at nk
                  m_Mlfs[g](i,j,nk) = (2*mMu[g](i,j,nk)+mLambda[g](i,j,nk))/(str_x(i)*str_y(j)); //(2*mu + lambda)/str at nk
               }
      }
   }
#undef str_x
#undef str_y
}
