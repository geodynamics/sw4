#include "MaterialParCartesian.h"
#include "EW.h"

MaterialParCartesian::MaterialParCartesian( EW* a_ew, int nx, int ny, int nz, int init, char* fname )
   : MaterialParameterization( a_ew, fname )
{
   //  VERIFY2( nx > 1 && ny > 1 && nz > 1, "MaterialParCartesian: The grid need at least two ponts in each direction")
     // Material represented on a coarse Cartesian grid, covering the 'active' domain.
     // points are x_0,..,x_{nx+1}, where x_0 and x_{nx+1} are fixed at zero.
   // the parameter vector represents offsets from a reference material, stored in (mRho,mMu,mLambda) in EW.

   m_init = init;
   m_nx = nx;
   m_ny = ny;
   m_nz = nz;

   double xmin, ymin, zmin, xmax, ymax, zmax;
   m_xmin = m_ymin = m_zmin =  1e38;
   m_xmax = m_ymax = m_zmax = -1e38;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      double hf = m_ew->mGridSize[g];
      xmin = (m_ew->m_iStartActGlobal[g]-1)*hf;
      ymin = (m_ew->m_jStartActGlobal[g]-1)*hf;
      xmax = (m_ew->m_iEndActGlobal[g]-1)*hf;
      ymax = (m_ew->m_jEndActGlobal[g]-1)*hf;
      zmax = m_ew->m_zmin[g] + (m_ew->m_kEndAct[g]-1)*hf;
      //      if( m_myrank == 0 )
      //	 cout << " xmin, xmax " << xmin << " " << xmax << endl;
      if( xmin < m_xmin )
	 m_xmin = xmin;
      if( ymin < m_ymin )
	 m_ymin = ymin;
      if( xmax > m_xmax )
	 m_xmax = xmax;
      if( ymax > m_ymax )
	 m_ymax = ymax;
      if( zmax > m_zmax )
         m_zmax = zmax;
  // z decreases when g increases, so zmin is always smallest on the last grid:
      m_zmin = m_ew->m_zmin[g];
      if( m_ew->topographyExists() && g == m_ew->mNumberOfGrids-1 )
      {
         zmin = 1e38;
	 for( int j= m_ew->m_jStartAct[g] ; j <= m_ew->m_jEndAct[g] ; j++ )
	    for( int i= m_ew->m_iStartAct[g] ; i <= m_ew->m_iEndAct[g] ; i++ )
	       if( m_ew->mZ(i,j,1) < zmin )
		  zmin = m_ew->mZ(i,j,1);
	 MPI_Allreduce( &zmin, &m_zmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
      }
   }
   //   double tmp[2] = {m_xmax,m_ymax}, vars[2];
   //   MPI_Allreduce( tmp, vars, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   //   m_xmax = vars[0];
   //   m_ymax = vars[1];

   //   tmp[0] = m_xmin;
   //   tmp[1] = m_ymin;
   //   MPI_Allreduce( tmp, vars, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
   //   m_xmin = vars[0];
   //   m_ymin = vars[1];
   //   cout << " xmin, xmax = " << m_xmin << " " << m_xmax << endl;
   //   cout << " ymin, ymax = " << m_ymin << " " << m_ymax << endl;

   // Determine h, such that x= i*h+xmin, i=0,..,nx+1

   m_hx = (m_xmax-m_xmin)/(nx+1);
   m_hy = (m_ymax-m_ymin)/(ny+1);
//   m_hz = (m_zmax-m_zmin)/(nz+1);
   m_hz = (m_zmax-m_zmin)/nz;
   m_zmin = m_zmin-m_hz;

   if( m_myrank == 0 )
   {
     cout << " xmin, xmax = " << m_xmin << " " << m_xmax << " hx = " << m_hx << endl;
     cout << " ymin, ymax = " << m_ymin << " " << m_ymax << " hy = " << m_hy << endl;
     cout << " zmin, zmax = " << m_zmin << " " << m_zmax << " hz = " << m_hz << endl;
   }
// Grid is determined.
   m_nms = 3*nx*ny*nz;
   m_nmd = 0;
   m_nmd_global = 0;
   m_rho.define(0,nx+1,0,ny+1,0,nz+1);
   m_mu.define(0,nx+1,0,ny+1,0,nz+1);
   m_lambda.define(0,nx+1,0,ny+1,0,nz+1);
   m_rho.set_to_zero();
   m_mu.set_to_zero();
   m_lambda.set_to_zero();
}

//-----------------------------------------------------------------------
void MaterialParCartesian::get_material( int nmd, double* xmd, int nms,
					 double* xms, vector<Sarray>& a_rho,
					 vector<Sarray>& a_mu, vector<Sarray>& a_lambda )
{
   double* rhop=m_rho.c_ptr();
   double* mup=m_mu.c_ptr();
   double* lambdap=m_lambda.c_ptr();
   size_t ind =0;
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
            rhop[indm]    = xms[3*ind];
	    mup[indm]     = xms[3*ind+1];
	    lambdap[indm] = xms[3*ind+2];
	    ind++;
	 }
   for( int g = 0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      m_ew->interpolate( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz, m_rho, m_mu,
			 m_lambda, g, a_rho[g], a_mu[g], a_lambda[g], true );
   }
}

//-----------------------------------------------------------------------
void MaterialParCartesian::interpolate_parameters( int nmd, double* xmd, int nms,
						   double* xms, std::vector<Sarray>& a_rho, 
			 std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )
{
   // Interpolates the difference a_rho-mRho, into local rho. i.e., m_rho = I(a_rho-mRho)
   // where a_rho,mRho are on the computational grid, rho on the parameter grid.
   //                             similarly for mu, lambda
   m_ew->interpolate_to_coarse( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz,
				m_rho, m_mu, m_lambda, a_rho, a_mu, a_lambda, true );
   double* rhop=m_rho.c_ptr();
   double* mup=m_mu.c_ptr();
   double* lambdap=m_lambda.c_ptr();
   size_t ind =0;
   double rhmin=100000,mumin=100000,lamin=100000;
   double rhmax=-100000,mumax=-100000,lamax=-100000; 
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
            xms[3*ind]   = rhop[indm];
	    xms[3*ind+1] = mup[indm];
	    xms[3*ind+2] = lambdap[indm];
	    if( rhop[indm]<rhmin )
	       rhmin = rhop[indm];
	    if( mup[indm]<mumin )
	       mumin = mup[indm];
	    if( lambdap[indm]<lamin )
	       lamin = lambdap[indm];
	    if( rhop[indm]>rhmax )
	    {
	       rhmax = rhop[indm];
       //               cout << " rhomax found " << m_myrank << " value " << rhop[ind] << " at " << i << " " << j << " " << k << endl;
	    }
	    if( mup[indm]>mumax )
	       mumax = mup[indm];
	    if( lambdap[indm]>lamax )
	       lamax = lambdap[indm];
	    ind++;
	 }
   //   cout << "proc " << m_myrank << " rho,mu,lambda min " << rhmin << " " <<mumin << " " << lamin << endl;
   //   cout << "proc " << m_myrank << " rho,mu,lambda max " << rhmax << " " <<mumax << " " << lamax << endl;
}

//-----------------------------------------------------------------------
void MaterialParCartesian::get_parameters( int nmd, double* xmd, int nms,
					   double* xms, std::vector<Sarray>& a_rho, 
					   std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )
{
   if( m_init == 0 )
   {
      for( int i=0 ; i < nms ; i++ )
	 xms[i] = 0;
   }
   else if( m_init == 1 )
   {
      //      cout << " hx, hy, hz " << m_hx  <<  " " << m_hy << " " << m_hz << endl;
      //      cout << " nx, ny, nz " << m_nx  <<  " " << m_ny << " " << m_nz << endl;
      //      cout << " xmin, xmax " << m_xmin  <<  " " << m_xmax << endl;
      //      cout << " ymin, ymax " << m_ymin  <<  " " << m_ymax << endl;
   // Test data for sine perturbed constant material
      double ep=0.01;
      double om = M_PI*2;
      size_t ind =0;
      for( int k=1 ; k <= m_nz ; k++ )
	 for( int j=1 ; j <= m_ny ; j++ )
	    for( int i=1 ; i <= m_nx ; i++ )
	    {
	       double x = i*m_hx + m_xmin;
	       double y = j*m_hy + m_ymin;
	       double z = k*m_hz + m_zmin;
	    
	       double rho = 1+ep*sin(om*x+0.13)*sin(om*y)*sin(om*z);
	       double cs  = 2+ep*cos(om*x)*sin(om*y)*cos(om*z+0.01);
	       double cp  = 4+ep*sin(om*x+0.4)*sin(om*y)*cos(om*z+0.1);
	       double mu = cs*cs*rho;
	       double lambda = rho*(cp*cp-2*cs*cs);
	       xms[3*ind]   = rho-1;
	       xms[3*ind+1] = mu-4;
	       xms[3*ind+2] = lambda-8;
	       //               if( m_myrank == 0 )
	       //		  cout << " xms " << xms[3*ind] << " " << xms[3*ind+1] << " " << xms[3*ind+2] << " x,y,z " << x << " " << y << " " << z << endl;
	       ind++;
	    }
   }
   else if( m_init == 2 )
   {
      read_parameters( nms, xms );
   }
   else if( m_init == 3 )
   {
      interpolate_parameters( nmd, xmd, nms, xms, a_rho, a_mu, a_lambda );
   }
}

//-----------------------------------------------------------------------
void MaterialParCartesian::get_gradient( int nmd, double* xmd, int nms, double* xms,
					 double* dfs, double* dfm,
					 std::vector<Sarray>& a_gradrho,
					 std::vector<Sarray>& a_gradmu,
					 std::vector<Sarray>& a_gradlambda )
{
   Sarray grho(0,m_nx+1,0,m_ny+1,0,m_nz+1), gmu(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray glambda(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   grho.set_to_zero();
   gmu.set_to_zero();
   glambda.set_to_zero();
   
   for( int g = 0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
   // Chain rule of interpolation relation
      m_ew->interpolation_gradient( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz,
			       grho, gmu, glambda, g, a_gradrho[g], a_gradmu[g], a_gradlambda[g] );
   }   
   double* grhop=grho.c_ptr();
   double* gmup=gmu.c_ptr();
   double* glambdap=glambda.c_ptr();
   int npts = (m_nx+2)*(m_ny+2)*(m_nz+2);

   double* tmp = new double[npts];
   for( int i=0 ; i < npts ; i++ )
      tmp[i] = grhop[i];
   MPI_Allreduce( tmp, grhop, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for( int i=0 ; i < npts ; i++ )
      tmp[i] = gmup[i];
   MPI_Allreduce( tmp, gmup, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for( int i=0 ; i < npts ; i++ )
      tmp[i] = glambdap[i];
   MPI_Allreduce( tmp, glambdap, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   delete[] tmp;

   size_t ind =0;
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
	    //            dfs[3*ind]   = grhop[ind];
	    //	    dfs[3*ind+1] = gmup[ind];
	    //	    dfs[3*ind+2] = glambdap[ind];
	    //	    ind++;
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
            dfs[3*ind]   = grhop[indm];
	    dfs[3*ind+1] = gmup[indm];
	    dfs[3*ind+2] = glambdap[indm];
	    ind++;
	 }
}

//-----------------------------------------------------------------------
//void MaterialParCartesian::perturb_material( int ip, int jp, int kp, int grid,
//					     int var, double h, double* xs, double* xm )
//// ignore grid, xm
//{
//   VERIFY2( 1 <= ip && ip <= m_nx && 1 <= jp && jp <= m_ny && 1 <= kp && kp <= m_nz,
//	    "ERROR in MaterialParCartesian::perturb_material, index (i,j,k) = " <<
//	    ip << " " << jp << " " << kp << " out of range " << endl );
//   VERIFY2( (var==0) || (var==1) || (var==2), 
//	    "ERROR in MaterialParCartesian::perturb_material, variable no. " << var
//	    << " out of range " << endl );
//   size_t ind = ip-1+m_nx*(jp-1)+m_nx*m_ny*(kp-1);
//   xs[3*ind+var] += h;
//}

//-----------------------------------------------------------------------
ssize_t MaterialParCartesian::parameter_index( int ip, int jp, int kp, int grid,
					    int var )
// Ignore grid.
{
   if( 1 <= ip && ip <= m_nx && 1 <= jp && jp <= m_ny && 1 <= kp && kp <= m_nz
       && (0 <= var) && (var <= 2 ) )
      return 3*(ip-1+static_cast<ssize_t>(m_nx)*(jp-1)+m_nx*m_ny*(kp-1))+var;
   else
      return -1;
}

//-----------------------------------------------------------------------
ssize_t MaterialParCartesian::local_index( size_t ind_global )
{
   return -1;
}
