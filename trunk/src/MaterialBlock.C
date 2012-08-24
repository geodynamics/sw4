// -*-c++-*-
#include "MaterialBlock.h"

#include <iostream>
#include "EW.h"

using namespace std;
//-----------------------------------------------------------------------
MaterialBlock::MaterialBlock( EW * a_ew, double rho, double vs, double vp, double xmin, 
                              double xmax, double ymin, double ymax, double zmin, double zmax,
			      double qs, double qp, double freq )
{
   m_rho = rho;
   m_vp  = vp;
   m_vs  = vs;
   m_xmin = xmin;
   m_xmax = xmax;
   m_ymin = ymin;
   m_ymax = ymax;
   m_zmin = zmin;
   m_zmax = zmax;
   m_tol = 1e-5;
   m_vpgrad  = 0;
   m_vsgrad  = 0;
   m_rhograd = 0;
   m_qs = qs;
   m_qp = qp;
   m_freq = freq;
   m_absoluteDepth = false;
   mEW = a_ew;

   double bbox[6];
   mEW->getGlobalBoundingBox( bbox );
  
// does this block cover all points?
// note that the global_zmin can be non-zero when topography is present
// global zmin is 0 in the absense of topography, and is assigned by the grid generator when there is topography
   if (xmin > bbox[0] || ymin > bbox[2] || zmin > bbox[4] ||
       xmax < bbox[1] || ymax < bbox[3] || zmax < bbox[5] )
   {
     mCoversAllPoints = false;
// tmp
//     cout << "This block does NOT cover all grid points" << endl;
   }
   else
   {
     mCoversAllPoints=true;
// tmp
//     cout << "This block COVERS all grid points" << endl;
   }
}

//-----------------------------------------------------------------------
void MaterialBlock::set_absoluteDepth( bool absDepth )
{
   m_absoluteDepth = absDepth;
}

//-----------------------------------------------------------------------
void MaterialBlock::set_gradients( double rhograd, double vsgrad, double vpgrad )
{
   m_rhograd = rhograd;
   m_vsgrad  = vsgrad;
   m_vpgrad  = vpgrad;
}

//-----------------------------------------------------------------------
bool MaterialBlock::inside_block( double x, double y, double z )
{
   return m_xmin-m_tol <= x && x <= m_xmax+m_tol && m_ymin-m_tol <= y && 
    y <= m_ymax+m_tol &&  m_zmin-m_tol <= z && z <= m_zmax+m_tol;
}

//-----------------------------------------------------------------------
void MaterialBlock::set_material_properties( std::vector<Sarray> & rho, 
                                             std::vector<Sarray> & cs,
                                             std::vector<Sarray> & cp, 
                                             std::vector<Sarray> & qs, 
                                             std::vector<Sarray> & qp)
{
  int pc[4];
// compute the number of parallel overlap points
//  mEW->interiorPaddingCells( pc );

  for( int g = 0 ; g < mEW->mNumberOfCartesianGrids; g++) // Cartesian grids
  {
// reference z-level for gradients is at z=0: AP changed this on 12/21/09
    double zsurf = 0.; // ?

    for( int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; k++ )
    {
      for( int j = mEW->m_jStart[g]; j <= mEW->m_jEnd[g]; j++ )
      {
	for( int i = mEW->m_iStart[g]; i <= mEW->m_iEnd[g] ; i++ )
	{
	  double x = (i-1)*mEW->mGridSize[g]                ;
	  double y = (j-1)*mEW->mGridSize[g]                ;
	  double z = mEW->m_zmin[g]+(k-1)*mEW->mGridSize[g];
                      
	  //printf("x ,y,z %f %f %f %f\n",x,y,z,mEW->m_zmin[g]);
                      
	  double depth;
	  if (m_absoluteDepth)
	  {
	    depth = z;
	  }
	  else
	  {
	    mEW->getDepth(x, y, z, depth);
	  }	  

	  if(inside_block(x,y,depth))
	  {
	    if( m_rho != -1 )
	      rho[g](i,j,k) = m_rho + m_rhograd*(depth-zsurf);
	    if( m_vs != -1 )
	      cs[g](i,j,k)  = m_vs + m_vsgrad*(depth-zsurf);
	    if( m_vp != -1 )
	      cp[g](i,j,k)  = m_vp + m_vpgrad*(depth-zsurf);
	    if( m_qp != -1 && qp[g].is_defined())
	      qp[g](i,j,k) = m_qp;
	    if( m_qs != -1 && qs[g].is_defined())
	      qs[g](i,j,k) = m_qs;
	  }
	}
      }
    }
    
// communicate material properties to ghost points (necessary on refined meshes because ghost points don't have a well defined depth/topography)
    mEW->communicate_array( rho[g], g );
    mEW->communicate_array( cs[g], g );
    mEW->communicate_array( cp[g], g );

    if (qs[g].is_defined())
      mEW->communicate_array( qs[g], g );
    if (qp[g].is_defined())
      mEW->communicate_array( qp[g], g );

  } // end for all Cartesian grids
  
  if (mEW->topographyExists()) // curvilinear grid
  {
    int g = mEW->mNumberOfGrids-1; 

// reference z-level for gradients is at z=0: AP changed this on 12/21/09
    double zsurf = 0.;

    for( int k = mEW->m_kStart[g] ; k <= mEW->m_kEnd[g]; k++ )
    {
      for( int j = mEW->m_jStart[g] ; j <= mEW->m_jEnd[g]; j++ )
      {
	for( int i = mEW->m_iStart[g] ; i <= mEW->m_iEnd[g] ; i++ )
	{
	  double x = mEW->mX(i,j,k);
	  double y = mEW->mY(i,j,k);
	  double z = mEW->mZ(i,j,k);
                      
	  //printf("x ,y,z %f %f %f %f\n",x,y,z,mEW->m_zmin[g]);
                      
	  double depth;
	  if (m_absoluteDepth)
	  {
	    depth = z;
	  }
	  else
	  {
	    mEW->getDepth(x, y, z, depth);
	  }	  

	  if(inside_block(x,y,depth))
	  {
	    if( m_rho != -1 )
	      rho[g](i,j,k) = m_rho + m_rhograd*(depth-zsurf);
	    if( m_vs != -1 )
	      cs[g](i,j,k)  = m_vs + m_vsgrad*(depth-zsurf);
	    if( m_vp != -1 )
	      cp[g](i,j,k)  = m_vp + m_vpgrad*(depth-zsurf);
	    if( m_qp != -1 && qp[g].is_defined())
	      qp[g](i,j,k) = m_qp;
	    if( m_qs != -1 && qs[g].is_defined())
	      qs[g](i,j,k) = m_qs;
	  }
	}
      }
    }
  } // end if topographyExists
} // end MaterialBlock::set_material_properties

//-----------------------------------------------------------------------
int MaterialBlock::set_material_pt( double x, double y, double z,
				    double& rho, double& cs, double& cp,
				    double& qs, double& qp )
{
   int retval = 0;
   if(inside_block(x,y,z))
   {
      double depth = z;
      double zsurf = 0;
      if( m_rho != -1 )
	 rho = m_rho + m_rhograd*(depth-zsurf);
      if( m_vs != -1 )
	 cs  = m_vs + m_vsgrad*(depth-zsurf);
      if( m_vp != -1 )
	 cp  = m_vp + m_vpgrad*(depth-zsurf);
      if( m_qp != -1 )
	 qp = m_qp;
      if( m_qs != -1 )
	 qs = m_qs;
   }
   else
      retval = -1;
   return retval;
}
