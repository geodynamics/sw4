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

// confusing with variables and functions which have names which only differ in capitalization 
  double localmin = localMin(mRho);
  MPI_Allreduce(&localmin,&mins[0],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  localmin = localMinVp();  
  MPI_Allreduce(&localmin,&mins[1],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  localmin = localMinVs();  
  MPI_Allreduce(&localmin,&mins[2],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  localmin = localMin(mMu);
  MPI_Allreduce(&localmin,&mins[3],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  localmin = localMin(mLambda);
  MPI_Allreduce(&localmin,&mins[4],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
  
  VERIFY2(mins[2] >= 0.0,
          "Error: the material data has s velocities that are negative.");
  localmin = localMinVpOverVs();  
  MPI_Allreduce(&localmin,&mins[5],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);

  double localmax = localMax(mRho);
  MPI_Allreduce(&localmax,&maxs[0],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  localmax = localMaxVp();  
  MPI_Allreduce(&localmax,&maxs[1],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  localmax = localMaxVs();  
  MPI_Allreduce(&localmax,&maxs[2],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  localmax = localMax(mMu);
  MPI_Allreduce(&localmax,&maxs[3],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  localmax = localMax(mLambda);
  MPI_Allreduce(&localmax,&maxs[4],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
  localmax = localMaxVpOverVs();  
  MPI_Allreduce(&localmax,&maxs[5],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);

  if( usingAttenuation() )
  {
      localmin = localMin(mQs);
      MPI_Allreduce(&localmin,&mins[6],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
      localmin = localMin(mQp);
      MPI_Allreduce(&localmin,&mins[7],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
      localmax = localMax(mQs);
      MPI_Allreduce(&localmax,&maxs[6],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
      localmax = localMax(mQp);
      MPI_Allreduce(&localmax,&maxs[7],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
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
  double localmin = a_field[0](m_iStart[0],m_jStart[0],m_kStart[0]);
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (a_field[g](i,j,k) < localmin)
                    {
                      localmin = a_field[g](i,j,k);
                    }
                }
            }
        }
    }

  return localmin; 
}

//-----------------------------------------------------------------------
double EW::localMax(std::vector<Sarray> & a_field) 
{
  double localmax = a_field[0](m_iStart[0],m_jStart[0],m_kStart[0]);
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (a_field[g](i,j,k) > localmax)
                    {
                      localmax = a_field[g](i,j,k);
                    }
                }
            }
        }
    }

  return localmax; 
}

//-----------------------------------------------------------------------
double EW::localMinVp() 
{
  double localmin = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k)) < localmin)
                    {
                      localmin = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return localmin; 
}

//-----------------------------------------------------------------------
double EW::localMaxVp() 
{
  double localmax = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k)) > localmax)
                    {
                      localmax = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return localmax; 
}

//-----------------------------------------------------------------------
double EW::localMinVs() 
{
  double localmin = sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])/mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if (sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) < localmin)
                    {
                      localmin = sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return localmin; 
}

//-----------------------------------------------------------------------
double EW::localMaxVs() 
{
  double localmax = sqrt(mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])/mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
      for (int k = m_kStart[g]; k <= m_kEnd[g]; k++ )
        {
          for (int j = m_jStart[g]; j <= m_jEnd[g]; j++ )
            {      
              for (int i = m_iStart[g]; i <= m_iEnd[g]; i++ )
                {
                  if ( sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) > localmax)
                    {
                      localmax = sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return localmax; 
}

//-----------------------------------------------------------------------
double EW::localMinVpOverVs() 
{
  double localmin = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
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
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) < localmin)
                    {
                      localmin = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return localmin; 
}

//-----------------------------------------------------------------------
double EW::localMaxVpOverVs() 
{
  double localmax = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
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
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) > localmax)
                    {
                      localmax = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
    }

  return localmax; 
}

// the same routine is defined in EtreeFile.C, but in the class EtreeFile
//-----------------------------------------------------------------------
void EW::extrapolateInZ(Sarray& field, bool useThreshold, double thresHoldValue, bool linear) 
{
  // -------------------------------------------------------------
  // We use linear extrapolation to update boundary points
  //
  // field(n+1) = 2*field(n) - field(n-1)
  //
  // Unless, 2*field(n) <= field(n-1), then we just set it to
  // field(n).
  // -------------------------------------------------------------
   //  int numGhostPoints = getExternalGhostPointsPerBoundaryPoint();

  int i, j, k;
  double extField;
  
// only extrapolate on the "low-k" side  
  for( j = field.m_jb; j <= field.m_je; j++ )
    for( i = field.m_ib; i <= field.m_ie ; i++ )
      for( k = field.m_kb + m_ghost_points-1 ; k >= field.m_kb; k-- )
      {
	if (linear && 2.*field(i,j,k+1) > field(i,j,k+2)) // check if extrapolation will lead to a negative value
	{
	  extField = 2.*field(i,j,k+1)-field(i,j,k+2);
	}
	else // constant extrapolation
	{
	  extField = field(i,j,k+1);
	}
	if (useThreshold && extField<thresHoldValue)
	  extField = thresHoldValue;
	
	field(i,j,k) = extField;
      } // end for k...
} // end extrapolateInZ

