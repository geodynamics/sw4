
#include "EW.h"
#include "MaterialParAllpts.h"


//-----------------------------------------------------------------------
MaterialParAllpts::MaterialParAllpts( EW* a_ew, char* fname ) : MaterialParameterization( a_ew, fname )
{
   m_nms = 0;
   m_nmd = 0;
   m_nmd_global = 0;
   m_npts_per_grid.resize(m_ew->mNumberOfGrids);
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      if( m_ew->m_iEndAct[g]-m_ew->m_iStartAct[g]+1 > 0 && m_ew->m_jEndAct[g]-m_ew->m_jStartAct[g]+1 >0
	  && m_ew->m_kEndAct[g]-m_ew->m_kStartAct[g]+1 > 0 )
	 m_nmd += (m_ew->m_iEndAct[g]-m_ew->m_iStartAct[g]+1)*(m_ew->m_jEndAct[g]-m_ew->m_jStartAct[g]+1)*
	    (m_ew->m_kEndAct[g]-m_ew->m_kStartAct[g]+1)*3;
      size_t gpts = (m_ew->m_iEndActGlobal[g]-m_ew->m_iStartActGlobal[g]+1)*
	    (m_ew->m_jEndActGlobal[g]-m_ew->m_jStartActGlobal[g]+1)*
	    (m_ew->m_kEndActGlobal[g]-m_ew->m_kStartActGlobal[g]+1)*3;
      m_nmd_global += gpts;
      m_npts_per_grid[g] = gpts;
   }

}

//-----------------------------------------------------------------------
void MaterialParAllpts::get_material( int nmd, double* xmd, int nms, double* xms, vector<Sarray>& a_rho,
				      vector<Sarray>& a_mu, vector<Sarray>& a_lambda )
{
   size_t gp, ind;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      // Fill in base material, in case the parameterization do not cover all of the domain.
      a_rho[g].copy( m_ew->mRho[g] );
      a_mu[g].copy( m_ew->mMu[g] );
      a_lambda[g].copy( m_ew->mLambda[g] );

      if( g == 0 )
	 gp = 0;
      else
	 gp = gp + 3*ind;

      ind =0;
      for( int k=m_ew->m_kStartAct[g]; k <= m_ew->m_kEndAct[g]; k++ )
	 for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	    for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
	    {
	       a_rho[g](i,j,k)    = xmd[gp+ind*3];
	       a_mu[g](i,j,k)     = xmd[gp+ind*3+1];
	       a_lambda[g](i,j,k) = xmd[gp+ind*3+2];
	       ind++;
	    }
   }
}

//-----------------------------------------------------------------------
void MaterialParAllpts::get_parameters( int nmd, double* xmd, int nms, double* xms, vector<Sarray>& a_rho,
					vector<Sarray>& a_mu, vector<Sarray>& a_lambda )
{
   size_t gp, ind;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      if( g == 0 )
	 gp = 0;
      else
	 gp = gp + 3*ind;
      ind =0;
      for( int k=m_ew->m_kStartAct[g]; k <= m_ew->m_kEndAct[g]; k++ )
	 for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	    for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
	    {
	       xmd[gp+ind*3]   = a_rho[g](i,j,k);
	       xmd[gp+ind*3+1] = a_mu[g](i,j,k);
	       xmd[gp+ind*3+2] = a_lambda[g](i,j,k);
	       ind++;
	    }
   }
}

//-----------------------------------------------------------------------
void MaterialParAllpts::get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfms, double* dfmd,
				      vector<Sarray>& a_gradrho, vector<Sarray>& a_gradmu,
				      vector<Sarray>& a_gradlambda )
{
   size_t gp, ind;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      if( g == 0 )
	 gp = 0;
      else
	 gp = gp + 3*ind;
      ind =0;
      for( int k=m_ew->m_kStartAct[g]; k <= m_ew->m_kEndAct[g]; k++ )
	 for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	    for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
	    {
	       dfmd[gp+ind*3]   = a_gradrho[g](i,j,k);
	       dfmd[gp+ind*3+1] = a_gradmu[g](i,j,k);
	       dfmd[gp+ind*3+2] = a_gradlambda[g](i,j,k);
	       ind++;
	    }
   }

}

////-----------------------------------------------------------------------
//void MaterialParAllpts::perturb_material( int ip, int jp, int kp, int grid,
//					  int var, double h, double* xs, double* xm )
//// ignore xs
//{
//   if( h != 0 && m_ew->point_in_proc(ip,jp,grid) )
//   {
//      if( ip < m_ew->m_iStartAct[grid] || ip > m_ew->m_iEndAct[grid] )
//	 cout << "warning i-index outside active domain " << endl;
//      if( jp < m_ew->m_jStartAct[grid] || jp > m_ew->m_jEndAct[grid] )
//	 cout << "warning j-index outside active domain " << endl;
//      if( kp < m_ew->m_kStartAct[grid] || kp > m_ew->m_kEndAct[grid] )
//	 cout << "warning k-index outside active domain " << endl;
//      int ni = m_ew->m_iEndAct[grid]-m_ew->m_iStartAct[grid]+1;
//      int nj = m_ew->m_jEndAct[grid]-m_ew->m_jStartAct[grid]+1;
//      size_t ind = (ip-m_ew->m_iStartAct[grid]) + ni*(jp-m_ew->m_jStartAct[grid]) + ni*nj*(kp-m_ew->m_kStartAct[grid]); 
//      xm[3*ind+var] += h;
//   }
//}

//-----------------------------------------------------------------------
ssize_t MaterialParAllpts::parameter_index( int ip, int jp, int kp, int grid,
					   int var )
{
   if( m_ew->point_in_proc(ip,jp,grid) )
   {
      size_t ni = static_cast<ssize_t>(m_ew->m_iEndAct[grid]-m_ew->m_iStartAct[grid]+1);
      size_t nj = static_cast<ssize_t>(m_ew->m_jEndAct[grid]-m_ew->m_jStartAct[grid]+1);
      size_t ind = 3*((ip-m_ew->m_iStartAct[grid]) + ni*(jp-m_ew->m_jStartAct[grid]) +
		      ni*nj*(kp-m_ew->m_kStartAct[grid]) )+var; 
      return ind;
   }
   else
      return -1;
}

//-----------------------------------------------------------------------
ssize_t MaterialParAllpts::local_index( size_t ind_global )
{
   size_t ind=m_npts_per_grid[0];
   bool found = false;
   int g=0;
   while( g < m_npts_per_grid.size() && !found )
   {
      if( ind_global > ind )
      {
	 g++;
	 ind += m_npts_per_grid[g];
      }
      else
	 found = true;
   }
   if( found )
   {

      size_t nig = static_cast<ssize_t>(m_ew->m_iEndActGlobal[g]-m_ew->m_iStartActGlobal[g]+1);
      size_t njg = static_cast<ssize_t>(m_ew->m_jEndActGlobal[g]-m_ew->m_jStartActGlobal[g]+1);
      int r = ind_global % 3;
      ind = (ind_global-r)/3;
      int i = ind % nig;
      ind = (ind-i)/nig;
      int j = ind % njg;
      ind = (ind-j)/njg;
      int k = ind;
      i = i + m_ew->m_iStartActGlobal[g];
      j = j + m_ew->m_jStartActGlobal[g];
      k = k + m_ew->m_kStartActGlobal[g];
      if( m_ew->m_iStartAct[g] <= i && i <= m_ew->m_iEndAct[g] &&
	  m_ew->m_jStartAct[g] <= j && j <= m_ew->m_jEndAct[g] &&
	  m_ew->m_kStartAct[g] <= k && k <= m_ew->m_kEndAct[g] )
      {
	 size_t ni = static_cast<ssize_t>(m_ew->m_iEndAct[g]-m_ew->m_iStartAct[g]+1);
	 size_t nj = static_cast<ssize_t>(m_ew->m_jEndAct[g]-m_ew->m_jStartAct[g]+1);
	 return 3*((i-m_ew->m_iStartAct[g]) + ni*(j-m_ew->m_jStartAct[g]) +
		   ni*nj*(k-m_ew->m_kStartAct[g]) )+r; 
      }
      else
	 return -1;
   }
   else
      return -1;
}
