
#include "EW.h"
#include "MaterialParAllpts.h"

//-----------------------------------------------------------------------
MaterialParAllpts::MaterialParAllpts( EW* a_ew, char* fname, int variables )
  : MaterialParameterization( a_ew, fname )
{
   if( variables == 0 )
   {
      m_variables = RML;
      m_nc = 3;
   }
   else if( variables == 1 )
   {
      m_variables = RCSCP;
      m_nc = 3;
   }
   else if( variables == 2 )
   {
      m_variables = CSCP;
      m_nc = 2;
   }
   else if( variables == 3 )
   {
      m_variables = CP;
      m_nc = 1;
   }
   else
   {
      cout <<  "ERROR: MaterialParAllpts: constructor, variables = " << variables 
	  << " not defined" << endl;
   }
   m_ratio = 1.732;

   m_nms = 0;
   m_nmd = 0;
   m_nmd_global = 0;
   m_npts_per_grid.resize(m_ew->mNumberOfGrids);
   m_npts_per_grid_local.resize(m_ew->mNumberOfGrids);
   //
   // The limits of the active region as defined in EW are iStartAct[g] <= i <= iEndAct[g] (and similar for j,k)
   // Only interior points are in active region.
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      m_npts_per_grid_local[g]=0;
      if( m_ew->m_iEndAct[g]-m_ew->m_iStartAct[g]+1 > 0 && 
          m_ew->m_jEndAct[g]-m_ew->m_jStartAct[g]+1 > 0 &&
	  m_ew->m_kEndAct[g]-m_ew->m_kStartAct[g]+1 > 0 )
      {
         int npts = (m_ew->m_iEndAct[g]-m_ew->m_iStartAct[g]+1)*
                    (m_ew->m_jEndAct[g]-m_ew->m_jStartAct[g]+1)*
           	    (m_ew->m_kEndAct[g]-m_ew->m_kStartAct[g]+1)*m_nc;
	 m_nmd += npts;
         m_npts_per_grid_local[g] = npts;
      }
      size_t gpts = (m_ew->m_iEndActGlobal[g]-m_ew->m_iStartActGlobal[g]+1)*
	            (m_ew->m_jEndActGlobal[g]-m_ew->m_jStartActGlobal[g]+1)*
    	            (m_ew->m_kEndActGlobal[g]-m_ew->m_kStartActGlobal[g]+1)*m_nc;
      m_nmd_global += gpts;
      m_npts_per_grid[g] = gpts;
   }

   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      if( m_myrank == 0 )
	cout << "active region, index = " <<
	  m_ew->m_iStartActGlobal[g] << " " <<
	  m_ew->m_iEndActGlobal[g] << " " <<
	  m_ew->m_jStartActGlobal[g] << " " <<
	  m_ew->m_jEndActGlobal[g] << " " <<
	  m_ew->m_kStartActGlobal[g] << " " <<
	  m_ew->m_kEndActGlobal[g] << endl;
   }

}

//-----------------------------------------------------------------------
void MaterialParAllpts::get_material( int nmd, double* xmd, int nms, double* xms, vector<Sarray>& a_rho,
				      vector<Sarray>& a_mu, vector<Sarray>& a_lambda)
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
	 gp = gp + m_nc*ind;

      ind =0;
      for( int k=m_ew->m_kStartAct[g]; k <= m_ew->m_kEndAct[g]; k++ )
	 for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	    for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
	    {
	       if( m_variables == RML )
	       {
		  a_rho[g](i,j,k)    = xmd[gp+ind*3  ];
		  a_mu[g](i,j,k)     = xmd[gp+ind*3+1];
		  a_lambda[g](i,j,k) = xmd[gp+ind*3+2];
	       }
	       else if( m_variables == RCSCP )
	       {
		  a_rho[g](i,j,k)    = xmd[gp+ind*3  ];
		  // mu=rho*cs*cs
		  a_mu[g](i,j,k)     = xmd[gp+ind*3]*xmd[gp+ind*3+1]*xmd[gp+ind*3+1];
		  // lambda=rho*cp*cp-2*mu=rho*(cp*cp-2*cs*cs)
		  a_lambda[g](i,j,k) = xmd[gp+ind*3]*(  xmd[gp+ind*3+2]*xmd[gp+ind*3+2]-
						      2*xmd[gp+ind*3+1]*xmd[gp+ind*3+1]  );
	       }
	       else if( m_variables == CSCP )
	       {
		  // Leave rho unchanged
		  a_mu[g](i,j,k)     = a_rho[g](i,j,k)*xmd[gp+ind*2]*xmd[gp+ind*2];
		  a_lambda[g](i,j,k) = a_rho[g](i,j,k)*( xmd[gp+ind*2+1]*xmd[gp+ind*2+1]-
						       2*xmd[gp+ind*2  ]*xmd[gp+ind*2  ]  );
	       }
	       else if( m_variables == CP )
	       {
		  // cs  = cp/ratio, rho = fixed.
		  // Note: this is not consistent with values in sponge layer, fix later.
		  double cp = xmd[gp+ind];
		  //		  double rho = cp*m_gamma;
		  double rho = a_rho[g](i,j,k);
		  double cs  = cp/m_ratio;
		  //		  a_rho[g](i,j,k)    = rho;
		  a_mu[g](i,j,k)     = rho*cs*cs;
		  a_lambda[g](i,j,k) = rho*(cp*cp-2*cs*cs);
	       }
	       ind++;
	    }
      m_ew->communicate_array( a_rho[g], g );
      m_ew->communicate_array( a_mu[g], g );
      m_ew->communicate_array( a_lambda[g], g );
   }
}

//-----------------------------------------------------------------------
void MaterialParAllpts::get_parameters( int nmd, double* xmd, int nms, double* xms, vector<Sarray>& a_rho,
					vector<Sarray>& a_mu, vector<Sarray>& a_lambda, int nr )
{
   size_t gp, ind;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      if( g == 0 )
	 gp = 0;
      else
	 gp = gp + m_nc*ind;
      ind =0;
      for( int k=m_ew->m_kStartAct[g]; k <= m_ew->m_kEndAct[g]; k++ )
	 for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	    for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
	    {
	       if( m_variables == RML )
	       {
		  xmd[gp+ind*3  ] = a_rho[g](i,j,k);
		  xmd[gp+ind*3+1] = a_mu[g](i,j,k);
		  xmd[gp+ind*3+2] = a_lambda[g](i,j,k);
	       }
	       else if( m_variables == RCSCP )
	       {
		  xmd[gp+ind*3  ] = a_rho[g](i,j,k);
		  xmd[gp+ind*3+1] = sqrt(a_mu[g](i,j,k)/a_rho[g](i,j,k));
		  xmd[gp+ind*3+2] = sqrt((2*a_mu[g](i,j,k)+a_lambda[g](i,j,k))/a_rho[g](i,j,k));
	       }
	       else if( m_variables == CSCP )
	       {
		  xmd[gp+ind*2  ] = sqrt(a_mu[g](i,j,k)/a_rho[g](i,j,k));
		  xmd[gp+ind*2+1] = sqrt((2*a_mu[g](i,j,k)+a_lambda[g](i,j,k))/a_rho[g](i,j,k));
	       }
	       else if( m_variables == CP )
		  xmd[gp+ind] = sqrt((2*a_mu[g](i,j,k)+a_lambda[g](i,j,k))/a_rho[g](i,j,k));
	       ind++;
	    }
   }
}



//-----------------------------------------------------------------------
void MaterialParAllpts::get_gradient( int nmd, double* xmd, int nms, double* xms, double* dfms, double* dfmd,
				      std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
				      std::vector<Sarray>& a_lambda, 
				      vector<Sarray>& a_gradrho, vector<Sarray>& a_gradmu,
				      vector<Sarray>& a_gradlambda)
{
   double irat2=1/(m_ratio*m_ratio);

   size_t gp, ind;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      if( g == 0 )
	 gp = 0;
      else
	 gp = gp + m_nc*ind;
      ind =0;
      for( int k=m_ew->m_kStartAct[g]; k <= m_ew->m_kEndAct[g]; k++ )
	 for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	    for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
	    {
	       if( m_variables == RML )
	       {
		  dfmd[gp+ind*3]   = a_gradrho[g](i,j,k);
		  dfmd[gp+ind*3+1] = a_gradmu[g](i,j,k);
		  dfmd[gp+ind*3+2] = a_gradlambda[g](i,j,k);
	       }
	       else if( m_variables == RCSCP )
	       {
		  double cs = sqrt(a_mu[g](i,j,k)/a_rho[g](i,j,k));
		  double cp = sqrt((a_lambda[g](i,j,k)+2*a_mu[g](i,j,k))/a_rho[g](i,j,k));
		  dfmd[gp+ind*3]   = (    a_mu[g](i,j,k)*a_gradmu[g](i,j,k)+
				      a_lambda[g](i,j,k)*a_gradlambda[g](i,j,k)+
					  a_rho[g](i,j,k)*a_gradrho[g](i,j,k))/a_rho[g](i,j,k);
		  dfmd[gp+ind*3+1] = 2*a_rho[g](i,j,k)*cs*(a_gradmu[g](i,j,k)-2*a_gradlambda[g](i,j,k));
		  dfmd[gp+ind*3+2] = 2*a_rho[g](i,j,k)*cp*a_gradlambda[g](i,j,k);
	       }
	       else if( m_variables == CSCP )
	       {
		  double cs = sqrt(a_mu[g](i,j,k)/a_rho[g](i,j,k));
		  double cp = sqrt((a_lambda[g](i,j,k)+2*a_mu[g](i,j,k))/a_rho[g](i,j,k));
		  dfmd[gp+ind*2  ] = 2*a_rho[g](i,j,k)*cs*(a_gradmu[g](i,j,k)-2*a_gradlambda[g](i,j,k));
		  dfmd[gp+ind*2+1] = 2*a_rho[g](i,j,k)*cp*a_gradlambda[g](i,j,k);
	       }
	       else if( m_variables == CP )
	       {
		  double cp = sqrt((a_lambda[g](i,j,k)+2*a_mu[g](i,j,k))/a_rho[g](i,j,k));
		  dfmd[gp+ind] = 2*cp*a_rho[g](i,j,k)*(irat2*a_gradmu[g](i,j,k)+
						       (1-2*irat2)*a_gradlambda[g](i,j,k));
	       }
	       ind++;
	    }
   }

}

//-----------------------------------------------------------------------
void MaterialParAllpts::interpolate_pseudohessian( int nmpars, double* phs, 
                                                   int nmpard, double* phm, 
                                                   vector<Sarray>& phgrid )
{
   size_t gp, ind;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      if( g == 0 )
	 gp = 0;
      else
	 gp = gp + m_nc*ind;
      ind =0;
      for( int k=m_ew->m_kStartAct[g]; k <= m_ew->m_kEndAct[g]; k++ )
	 for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	    for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
            {
               if(  m_variables == RML || m_variables == RCSCP )
               {
                  phm[gp+ind*3  ] = phgrid[g](1,i,j,k);
                  phm[gp+ind*3+1] = phgrid[g](2,i,j,k);
                  phm[gp+ind*3+2] = phgrid[g](3,i,j,k);
               }
               else if( m_variables == CSCP )
               {
                  phm[gp+ind*2  ] = phgrid[g](2,i,j,k);
                  phm[gp+ind*2+1] = phgrid[g](3,i,j,k);
               }
               else if( m_variables == CP )
               {
                  phm[gp+ind] = phgrid[g](3,i,j,k);
               }
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
ssize_t MaterialParAllpts::parameter_index( int ip, int jp, int kp,
                                            int grid, int var )
{
   // Computes local index in xmd, from the global sw4-array index (ip,jp,kp)
   // var=0,1,or 2 
   size_t gp = 0;   // Local grid pointer
   for( int g=0; g < grid; g++)
      gp += m_npts_per_grid_local[g];
   ssize_t ind =-1;
   if( m_ew->point_in_proc(ip,jp,grid) )
   {
      if( m_ew->m_iStartAct[grid] <= ip && ip <= m_ew->m_iEndAct[grid] &&
	  m_ew->m_jStartAct[grid] <= jp && jp <= m_ew->m_jEndAct[grid] &&
	  m_ew->m_kStartAct[grid] <= kp && kp <= m_ew->m_kEndAct[grid] )
      {
	 size_t ni = static_cast<ssize_t>(m_ew->m_iEndAct[grid]-m_ew->m_iStartAct[grid]+1);
	 size_t nj = static_cast<ssize_t>(m_ew->m_jEndAct[grid]-m_ew->m_jStartAct[grid]+1);
	 ind = m_nc*((ip-m_ew->m_iStartAct[grid]) + ni*(jp-m_ew->m_jStartAct[grid]) +
		      ni*nj*(kp-m_ew->m_kStartAct[grid]) )+var+gp; 
      }
   }
   return ind;
}

//-----------------------------------------------------------------------
ssize_t MaterialParAllpts::local_index( size_t ind_global )
{
  // ind = c + nc*(i-ib) + nc*ni*(j-jb) + nc*ni*nj*(k-ib)
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

      // Transform to a global index (i,j,k)
      size_t nig = static_cast<ssize_t>(m_ew->m_iEndActGlobal[g]-m_ew->m_iStartActGlobal[g]+1);
      size_t njg = static_cast<ssize_t>(m_ew->m_jEndActGlobal[g]-m_ew->m_jStartActGlobal[g]+1);
      int r = ind_global % m_nc;
      ind = (ind_global-r)/m_nc;
      int i = ind % nig;
      ind = (ind-i)/nig;
      int j = ind % njg;
      ind = (ind-j)/njg;
      int k = ind;
      i = i + m_ew->m_iStartActGlobal[g];
      j = j + m_ew->m_jStartActGlobal[g];
      k = k + m_ew->m_kStartActGlobal[g];
      return parameter_index(i,j,k,g,r);

      //      if( m_ew->m_iStartAct[g] <= i && i <= m_ew->m_iEndAct[g] &&
      //	  m_ew->m_jStartAct[g] <= j && j <= m_ew->m_jEndAct[g] &&
      //	  m_ew->m_kStartAct[g] <= k && k <= m_ew->m_kEndAct[g] )
      //      {
      //	 size_t ni = static_cast<ssize_t>(m_ew->m_iEndAct[g]-m_ew->m_iStartAct[g]+1);
      //	 size_t nj = static_cast<ssize_t>(m_ew->m_jEndAct[g]-m_ew->m_jStartAct[g]+1);
      //	 return m_nc*((i-m_ew->m_iStartAct[g]) + ni*(j-m_ew->m_jStartAct[g]) +
      //		   ni*nj*(k-m_ew->m_kStartAct[g]) )+r; 
      //      }
      //      else
      //	 return -1;
   }
   else
      return -1;
}

//-----------------------------------------------------------------------
void MaterialParAllpts::set_scalefactors( int nmpars, double* sfs, 
                  double rho_ref, double mu_ref, double lambda_ref, double vs_ref, double vp_ref )
{
   if( m_variables == RML )
   {
      for( int i=0 ; i < nmpars ; i += 3 )
      {
	 sfs[i]   = rho_ref;
	 sfs[i+1] = mu_ref;
	 sfs[i+2] = lambda_ref;
      }
   }
   else if( m_variables == RCSCP )
   {
      for( int i=0 ; i < nmpars ; i += 3 )
      {
	 sfs[i]   = rho_ref;
	 sfs[i+1] = vs_ref;
	 sfs[i+2] = vp_ref;
      }
   }
   else if( m_variables == CSCP )
   {
      for( int i=0 ; i < nmpars ; i += 2 )
      {
	 sfs[i]   = vs_ref;
	 sfs[i+1] = vp_ref;
      }
   }
   else if( m_variables == CP )
   {
      for( int i=0 ; i < nmpars ; i++ )
      {
	 sfs[i] = vp_ref;
      }
   }
}


//-----------------------------------------------------------------------
void MaterialParAllpts::get_regularizer( int nmd, double* xmd, int nms, double* xms, 
					 double* xmd0, double* xms0, double regcoeff,
					 vector<Sarray>& a_rho, vector<Sarray>& a_mu,
					 vector<Sarray>& a_lambda, double& mf_reg,
					 double* sfd, double* sfs, 
					 bool compute_derivative, double* dmfd_reg,
					 double* dmfs_reg )
{
   if( regcoeff == 0 )
      return;

#define SQR(x) ((x)*(x))
   // Laplacian regularization
   Sarray cs, cp;
   int npts = 0;
   size_t ind0 = 0;
   int o=0;
   if( m_variables == RCSCP )
      o = 1;
   int myrank;
   MPI_Comm_rank(m_ew->m_1d_communicator,&myrank);

   mf_reg = 0;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      const int iba = m_ew->m_iStartAct[g];
      const int jba = m_ew->m_jStartAct[g];
      const int kba = m_ew->m_kStartAct[g];
      const int ni = ( m_ew->m_iEndAct[g]- m_ew->m_iStartAct[g]+1);
      const int nj = ( m_ew->m_jEndAct[g]- m_ew->m_jStartAct[g]+1);
      const int nk = ( m_ew->m_kEndAct[g]- m_ew->m_kStartAct[g]+1);
      const int nij = ni*nj;

      double rhoref=1.0, muref=1.0, laref=1.0;
      double csref=1.0,cpref=1.0;
      //      cout << myrank<<"act domain: " << iba << " " << ni << " " << jba << " " << nj
      //	   << " " << kba << " " << nk << endl;
      //      cout << myrank<<"nmd = " << nmd << " " << m_nmd << endl;
      //      cout << myrank<<"nms = " << nms << " " << m_nms << endl;
      if( m_variables == RML )
      {
	 if( m_nmd > 0 )
	 {
	    rhoref = sfd[ind0];
	    muref  = sfd[ind0+1];
	    laref  = sfd[ind0+2];
	 }
      }
      else if( m_variables == RCSCP || m_variables == CSCP )
      {
	 cs.define(a_mu[g]);
	 cp.define(a_lambda[g]);
	 double* csptr  = cs.c_ptr();
	 double* cpptr  = cp.c_ptr();
	 double* rhoptr = a_rho[g].c_ptr();
	 double* muptr  = a_mu[g].c_ptr();
	 double* laptr  = a_lambda[g].c_ptr();
	 for( int i = 0 ; i < cs.m_npts ; i++ )
	 {
	    csptr[i] = sqrt(muptr[i]/rhoptr[i]);
	    cpptr[i] = sqrt((2*muptr[i]+laptr[i])/rhoptr[i]);
	 }
	 csref = csptr[0];
	 cpref = cpptr[0];
	 rhoref=rhoptr[0];
      }
      else if( m_variables == CP )
      {
	 cp.define(a_lambda[g]);
	 double* cpptr  = cp.c_ptr();
	 double* rhoptr = a_rho[g].c_ptr();
	 double* muptr  = a_mu[g].c_ptr();
	 double* laptr  = a_lambda[g].c_ptr();
	 for( int i = 0 ; i < cp.m_npts ; i++ )
	 {
	    cpptr[i] = sqrt((2*muptr[i]+laptr[i])/rhoptr[i]);
	 }
	 cpref = cpptr[0];
      }
      // Scalar scale factors for now, turned out that using the full sfd array 
      // requires reworking of all formulas, since it enters as a variable coefficient.
      double irh2 = 1/(rhoref*rhoref);
      double imu2 = 1/(muref*muref);
      double ila2 = 1/(laref*laref);
      double ics2 = 1/(csref*csref);
      double icp2 = 1/(cpref*cpref);

      size_t ind;
      for( int k=m_ew->m_kStartAct[g]; k <= m_ew->m_kEndAct[g]; k++ )
	 for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	    for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
	    {
	       ind = ind0 + m_nc*(i-iba+ni*(j-jba)+nij*(k-kba));
	       if( m_variables == RCSCP || m_variables == RML )
	       {
		  mf_reg += (SQR(2*a_rho[g](i,j,k)-a_rho[g](i+1,j,k)-a_rho[g](i-1,j,k))+
			     SQR(2*a_rho[g](i,j,k)-a_rho[g](i,j+1,k)-a_rho[g](i,j-1,k)))*irh2;
		  if( k > 1 )
		     mf_reg += SQR(2*a_rho[g](i,j,k)-a_rho[g](i,j,k+1)-a_rho[g](i,j,k-1))*irh2;
		  else
		     mf_reg += SQR(2*(a_rho[g](i,j,k+1)-a_rho[g](i,j,k)))*irh2;
	       }
	       if( m_variables == RML )
	       {
		  mf_reg += (SQR(2*a_mu[g](i,j,k)-a_mu[g](i+1,j,k)-a_mu[g](i-1,j,k))+
			     SQR(2*a_mu[g](i,j,k)-a_mu[g](i,j+1,k)-a_mu[g](i,j-1,k)))*imu2;
		     mf_reg += ( k > 1 ) ? SQR(2*a_mu[g](i,j,k)-a_mu[g](i,j,k+1)-a_mu[g](i,j,k-1))*imu2:
			SQR(2*(a_mu[g](i,j,k+1)-a_mu[g](i,j,k)))*imu2;

		  mf_reg += (SQR(2*a_lambda[g](i,j,k)-a_lambda[g](i+1,j,k)-a_lambda[g](i-1,j,k))+
			     SQR(2*a_lambda[g](i,j,k)-a_lambda[g](i,j+1,k)-a_lambda[g](i,j-1,k))*ila2);
		  mf_reg += ( k > 1 ) ?
		     SQR(2*a_lambda[g](i,j,k)-a_lambda[g](i,j,k+1)-a_lambda[g](i,j,k-1))*ila2:
		     SQR(2*(a_lambda[g](i,j,k+1)-a_lambda[g](i,j,k)))*ila2;
	       }
	       else if( m_variables == RCSCP || m_variables == CSCP )
	       {
		  mf_reg += (SQR(2*cs(i,j,k)-cs(i+1,j,k)-cs(i-1,j,k))+
			     SQR(2*cs(i,j,k)-cs(i,j+1,k)-cs(i,j-1,k)) )*ics2;
		  mf_reg += ( k > 1 ) ? SQR(2*cs(i,j,k)-cs(i,j,k+1)-cs(i,j,k-1))*ics2:
		     SQR(2*(cs(i,j,k+1)-cs(i,j,k)))*ics2;
		  mf_reg += (SQR(2*cp(i,j,k)-cp(i+1,j,k)-cp(i-1,j,k)) +
			     SQR(2*cp(i,j,k)-cp(i,j+1,k)-cp(i,j-1,k)))*icp2;
		  mf_reg += ( k > 1 ) ? SQR(2*cp(i,j,k)-cp(i,j,k+1)-cp(i,j,k-1))*icp2:
		     SQR(2*(cp(i,j,k+1)-cp(i,j,k)))*icp2;
	       }
	       else if( m_variables == CP )
	       {
		  mf_reg += (SQR(2*cp(i,j,k)-cp(i+1,j,k)-cp(i-1,j,k))+
			     SQR(2*cp(i,j,k)-cp(i,j+1,k)-cp(i,j-1,k)))*icp2;
		  mf_reg += ( k > 1 ) ? SQR(2*cp(i,j,k)-cp(i,j,k+1)-cp(i,j,k-1))*icp2:
		     SQR(2*(cp(i,j,k+1)-cp(i,j,k)))*icp2;
	       }
	       npts++;
	    }
      if( compute_derivative )
      {

	 for( int i = 0 ; i < m_nc*nij*nk ; i++ )
	    dmfd_reg[ind0+i] = 0;

	 int kstart = m_ew->m_kStartAct[g];
	 if( kstart == 1 )
	 {
	    // Boundary modified regularizer in the k-direction
	    double rcof[8]={1.25,-1.5,0.25,0,-1.5,2.25,-1.0,0.25};
	    for( int k=1 ; k <= 2 ; k++ )
	       for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
		  for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
		  {
		     ind = ind0 + m_nc*(i-iba+ni*(j-jba)+nij*(k-kba));
		     if( m_variables == RML || m_variables == RCSCP )
		     {
			dmfd_reg[ind] += (3*a_rho[g](i,j,k)
			   -(a_rho[g](i+1,j,k)+a_rho[g](i-1,j,k) +
			     a_rho[g](i,j+1,k)+a_rho[g](i,j-1,k)) +
			   +0.25*(a_rho[g](i+2,j,k)+a_rho[g](i-2,j,k) +
				  a_rho[g](i,j+2,k)+a_rho[g](i,j-2,k)) +
			   rcof[4*(k-1)  ]*a_rho[g](i,j,1) + rcof[4*(k-1)+1]*a_rho[g](i,j,2)+
			   rcof[4*(k-1)+2]*a_rho[g](i,j,3) + rcof[4*(k-1)+3]*a_rho[g](i,j,4))*irh2;
		     }
		     if( m_variables == RML )
		     {
			dmfd_reg[ind+1] += (3*a_mu[g](i,j,k)
			   -(a_mu[g](i+1,j,k)+a_mu[g](i-1,j,k) +
			     a_mu[g](i,j+1,k)+a_mu[g](i,j-1,k)) +
			   +0.25*(a_mu[g](i+2,j,k)+a_mu[g](i-2,j,k) +
				  a_mu[g](i,j+2,k)+a_mu[g](i,j-2,k)) +
			   rcof[4*(k-1)  ]*a_mu[g](i,j,1) + rcof[4*(k-1)+1]*a_mu[g](i,j,2)+
			   rcof[4*(k-1)+2]*a_mu[g](i,j,3) + rcof[4*(k-1)+3]*a_mu[g](i,j,4))*imu2;
			dmfd_reg[ind+2] += (3*a_lambda[g](i,j,k)
			   -(a_lambda[g](i+1,j,k)+a_lambda[g](i-1,j,k) +
			     a_lambda[g](i,j+1,k)+a_lambda[g](i,j-1,k)) +
			   +0.25*(a_lambda[g](i+2,j,k)+a_lambda[g](i-2,j,k) +
				  a_lambda[g](i,j+2,k)+a_lambda[g](i,j-2,k)) +
			   rcof[4*(k-1)  ]*a_lambda[g](i,j,1) + rcof[4*(k-1)+1]*a_lambda[g](i,j,2)+
			   rcof[4*(k-1)+2]*a_lambda[g](i,j,3) + rcof[4*(k-1)+3]*a_lambda[g](i,j,4))*ila2;
		     }
		     else if( m_variables == RCSCP || m_variables == CSCP )
		     {
			dmfd_reg[ind+o] += (3*cs(i,j,k)
			   -(cs(i+1,j,k)+cs(i-1,j,k) +
			     cs(i,j+1,k)+cs(i,j-1,k)) +
			   +0.25*(cs(i+2,j,k)+cs(i-2,j,k) +
				  cs(i,j+2,k)+cs(i,j-2,k)) +
			   rcof[4*(k-1)  ]*cs(i,j,1) + rcof[4*(k-1)+1]*cs(i,j,2)+
			   rcof[4*(k-1)+2]*cs(i,j,3) + rcof[4*(k-1)+3]*cs(i,j,4))*ics2;
			dmfd_reg[ind+o+1] += (3*cp(i,j,k)
			   -(cp(i+1,j,k)+cp(i-1,j,k) +
			     cp(i,j+1,k)+cp(i,j-1,k)) +
			   +0.25*(cp(i+2,j,k)+cp(i-2,j,k) +
				  cp(i,j+2,k)+cp(i,j-2,k)) +
			   rcof[4*(k-1)  ]*cp(i,j,1) + rcof[4*(k-1)+1]*cp(i,j,2)+
		           rcof[4*(k-1)+2]*cp(i,j,3) + rcof[4*(k-1)+3]*cp(i,j,4))*icp2;
		     }
		     else if( m_variables == CP )
		     {
			dmfd_reg[ind] += (3*cp(i,j,k)
			   -(cp(i+1,j,k)+cp(i-1,j,k) +
			     cp(i,j+1,k)+cp(i,j-1,k)) +
			   +0.25*(cp(i+2,j,k)+cp(i-2,j,k) +
				  cp(i,j+2,k)+cp(i,j-2,k)) +
			   rcof[4*(k-1)  ]*cp(i,j,1) + rcof[4*(k-1)+1]*cp(i,j,2)+
			   rcof[4*(k-1)+2]*cp(i,j,3) + rcof[4*(k-1)+3]*cp(i,j,4))*icp2;
		     }
		  }
	    kstart = 3;
	 }
	 for( int k=kstart ; k <= m_ew->m_kEndAct[g]; k++ )
	    for( int j=m_ew->m_jStartAct[g]; j <= m_ew->m_jEndAct[g]; j++ )
	       for( int i=m_ew->m_iStartAct[g]; i <= m_ew->m_iEndAct[g]; i++ )
	       {
		  ind = ind0 + m_nc*(i-iba+ni*(j-jba)+nij*(k-kba));
		  if( m_variables == RML || m_variables == RCSCP )
		     dmfd_reg[ind] += (4.5*a_rho[g](i,j,k)
			-(a_rho[g](i+1,j,k)+a_rho[g](i-1,j,k) +
			  a_rho[g](i,j+1,k)+a_rho[g](i,j-1,k) +
			  a_rho[g](i,j,k+1)+a_rho[g](i,j,k-1) ) 
			+0.25*(a_rho[g](i+2,j,k)+a_rho[g](i-2,j,k) +
			       a_rho[g](i,j+2,k)+a_rho[g](i,j-2,k) +
			       a_rho[g](i,j,k+2)+a_rho[g](i,j,k-2) ))*irh2;
		  if( m_variables == RML )
		  {
		     dmfd_reg[ind+1] += (4.5*a_mu[g](i,j,k)
			-(a_mu[g](i+1,j,k)+a_mu[g](i-1,j,k) +
			  a_mu[g](i,j+1,k)+a_mu[g](i,j-1,k) +
			  a_mu[g](i,j,k+1)+a_mu[g](i,j,k-1) ) 
			+0.25*(a_mu[g](i+2,j,k)+a_mu[g](i-2,j,k) +
			       a_mu[g](i,j+2,k)+a_mu[g](i,j-2,k) +
			       a_mu[g](i,j,k+2)+a_mu[g](i,j,k-2) ))*imu2;
		     dmfd_reg[ind+2] += (4.5*a_lambda[g](i,j,k)
			-(a_lambda[g](i+1,j,k)+a_lambda[g](i-1,j,k) +
			  a_lambda[g](i,j+1,k)+a_lambda[g](i,j-1,k) +
			  a_lambda[g](i,j,k+1)+a_lambda[g](i,j,k-1) ) 
			+0.25*(a_lambda[g](i+2,j,k)+a_lambda[g](i-2,j,k) +
			       a_lambda[g](i,j+2,k)+a_lambda[g](i,j-2,k) +
			       a_lambda[g](i,j,k+2)+a_lambda[g](i,j,k-2) ))*ila2;
		  }
		  else if( m_variables == RCSCP || m_variables == CSCP )
		  {
		     dmfd_reg[ind+o] += (4.5*cs(i,j,k)
			-(cs(i+1,j,k)+cs(i-1,j,k) +
			  cs(i,j+1,k)+cs(i,j-1,k) +
			  cs(i,j,k+1)+cs(i,j,k-1) ) 
			+0.25*(cs(i+2,j,k)+cs(i-2,j,k) +
			       cs(i,j+2,k)+cs(i,j-2,k) +
			       cs(i,j,k+2)+cs(i,j,k-2) ))*ics2;
		     dmfd_reg[ind+o+1] += (4.5*cp(i,j,k)
			-(cp(i+1,j,k)+cp(i-1,j,k) +
			  cp(i,j+1,k)+cp(i,j-1,k) +
			  cp(i,j,k+1)+cp(i,j,k-1) ) 
			+0.25*(cp(i+2,j,k)+cp(i-2,j,k) +
			       cp(i,j+2,k)+cp(i,j-2,k) +
			       cp(i,j,k+2)+cp(i,j,k-2) ))*icp2;
		  }
		  else if( m_variables == CP )
		  {
		     dmfd_reg[ind] += (4.5*cp(i,j,k)
			-(cp(i+1,j,k)+cp(i-1,j,k) +
			  cp(i,j+1,k)+cp(i,j-1,k) +
			  cp(i,j,k+1)+cp(i,j,k-1) ) 
			+0.25*(cp(i+2,j,k)+cp(i-2,j,k) +
			       cp(i,j+2,k)+cp(i,j-2,k) +
			       cp(i,j,k+2)+cp(i,j,k-2) ))*icp2;
		  }
	       }
      }
      ind0 += m_npts_per_grid[g];
   }
   double mf_reg_tmp=mf_reg;
   MPI_Allreduce( &mf_reg_tmp, &mf_reg, 1, MPI_DOUBLE, MPI_SUM, m_ew->m_1d_communicator);
   int npts_tmp = npts;
   MPI_Allreduce( &npts_tmp, &npts, 1, MPI_INT, MPI_SUM, m_ew->m_1d_communicator );

   double inpts = 1.0/npts;
   mf_reg = 0.125*regcoeff*mf_reg*inpts;
   if( compute_derivative )
   {
      for( int i=0; i < m_nmd ; i++ )
	 dmfd_reg[i] *= regcoeff*inpts;
   }
#undef SQR
}

//-----------------------------------------------------------------------
int MaterialParAllpts::get_varcase()
{
   if( m_variables == RML )
      return 1;
   else if( m_variables == RCSCP )
      return 2;
   else if( m_variables == CSCP )
      return 3;
   else if( m_variables == CP )
      return 4;
   else
      return 0;
}


void MaterialParAllpts::get_material( int nmd, double* xmd, int nms,
					     double* xms, std::vector<Sarray>& a_rho,
					     std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda,
                    float_sw4 vp_min, float_sw4 vp_max, float_sw4 vs_min, float_sw4 vs_max, int wave_mode)
{

}
