#include "MaterialParCartesianVels.h"
#include "EW.h"

//-----------------------------------------------------------------------
//  Parameterize the material on a Cartesian coarse grid, with
//  (density, cs, cp)-update on a fixed material as the parameters. 
//  
//  rho=x[3*i]
//  cs =x[3*i+1]
//  cp =x[3*i+2]
//
//-----------------------------------------------------------------------
MaterialParCartesianVels::MaterialParCartesianVels( EW* a_ew, int nx, int ny, int nz, int init, char* fname )
   : MaterialParameterization( a_ew, fname )
{
   //  VERIFY2( nx > 1 && ny > 1 && nz > 1, "MaterialParCartesianVels: The grid need at least two ponts in each direction")
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
   m_cs.define(0,nx+1,0,ny+1,0,nz+1);
   m_cp.define(0,nx+1,0,ny+1,0,nz+1);

   m_rho.set_to_zero();
   m_mu.set_to_zero();
   m_lambda.set_to_zero();
   m_cs.set_to_zero();
   m_cp.set_to_zero();
}

//-----------------------------------------------------------------------
// Input parameter vector (xmd,xms) and output the material (a_rho,a_mu,a_lambda).
//
// mu = rho*cs*cs, 2*mu+lambda = rho*cp*cp --> lambda=rho*cp*cp-2*mu
//
//-----------------------------------------------------------------------
void MaterialParCartesianVels::get_material( int nmd, double* xmd, int nms,
					     double* xms, vector<Sarray>& a_rho,
					     vector<Sarray>& a_mu, vector<Sarray>& a_lambda )
{
   // 1.  (rho,cs,cp) := x
   // 2.  (a_rho,a_cs,a_cp) := I(rho,cs,cp)  where I(rho,cs,cp) is interpolation to f.d. grid.
   // 3.  (a_rho,a_mu,a_lambda) := T( mRho+a_rho,mCs+a_cs,mCp+a_cp), where T transforms from (rho,cs,cp) to (rho,mu,lambda)
   //          mRho,mCs,mCp is base material.
   double* rhop = m_rho.c_ptr();
   double* mup  = m_mu.c_ptr();
   double* lambdap=m_lambda.c_ptr();
   double* csp = m_cs.c_ptr();
   double* cpp = m_cp.c_ptr();
   size_t ind =0;
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
            rhop[indm] = xms[3*ind];
	    csp[indm]  = xms[3*ind+1];
	    cpp[indm]  = xms[3*ind+2];
	    //	    mup[indm]    = rhop[indm]*csp[indm]*csp[indm];
	    //	    lambdap[indm]= rhop[indm]*cpp[indm]*cpp[indm]-2*mup[indm];
	    ind++;
	 }
   for( int g = 0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      //      m_ew->interpolate( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz, m_rho, m_mu,
      //			 m_lambda, g, a_rho[g], a_mu[g], a_lambda[g], false );
      m_ew->interpolate( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz, m_rho, m_cs,
			 m_cp, g, a_rho[g], a_mu[g], a_lambda[g], false );
      m_ew->update_and_transform_material( g, a_rho[g], a_mu[g], a_lambda[g] );
   }
}

//-----------------------------------------------------------------------
// Input material on grid(s) (a_rho,a_mu,a_lambda), output corresponding
// parameter vector (xmd,xms).
//-----------------------------------------------------------------------
void MaterialParCartesianVels::interpolate_parameters( int nmd, double* xmd, int nms,
						   double* xms, std::vector<Sarray>& a_rho, 
			 std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )
{
   // Interpolates the difference a_rho-(m_ew->mRho), into local rho. i.e., m_rho = I(a_rho-(m_ew->mRho))
   // where a_rho,mRho are on the computational grid, rho on the parameter grid.
   //                             similarly for cs, cp
   // m_cs = I( a_cs-mCs)
   // m_cp = I( a_cp-mCp)
   //

   m_ew->interpolate_to_coarse_vel( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz,
				    m_rho, m_cs, m_cp, a_rho, a_mu, a_lambda );
   float_sw4* rhop=m_rho.c_ptr();
   float_sw4* csp =m_cs.c_ptr();
   float_sw4* cpp =m_cp.c_ptr();
   size_t ind =0;
   double rhmin=100000,csmin=100000,cpmin=100000;
   double rhmax=-100000,csmax=-100000,cpmax=-100000; 
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
            xms[3*ind]   = rhop[indm];
            xms[3*ind+1] = csp[indm];
	    xms[3*ind+2] = cpp[indm];
	    if( rhop[indm]<rhmin )
	       rhmin = rhop[indm];
	    if( csp[indm]<csmin )
	       csmin = csp[indm];
	    if( cpp[indm]<cpmin )
	       cpmin = cpp[indm];
	    if( rhop[indm]>rhmax )
	       rhmax = rhop[indm];
	    if( csp[indm]>csmax )
	       csmax = csp[indm];
	    if( cpp[indm]>cpmax )
	       cpmax = cpp[indm];
	    ind++;
	 }
}

//-----------------------------------------------------------------------
void MaterialParCartesianVels::get_parameters( int nmd, double* xmd, int nms,
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
	       xms[3*ind+1] = cs-2;
	       xms[3*ind+2] = cp-4;
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
void MaterialParCartesianVels::get_gradient( int nmd, double* xmd, int nms, double* xms,
					 double* dfs, double* dfm,
					 std::vector<Sarray>& a_gradrho,
					 std::vector<Sarray>& a_gradmu,
					 std::vector<Sarray>& a_gradlambda )
{
   // Computes gradient with respect to the material parameterization from given
   // gradients with respect to the material at grid points.
   // 
   Sarray grho(0,m_nx+1,0,m_ny+1,0,m_nz+1), gmu(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray glambda(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   grho.set_to_zero();
   gmu.set_to_zero();
   glambda.set_to_zero();
   //   cout << "getgrad nms=" << nms << " nx,ny,nz=" << m_nx <<","<<m_ny <<","<<m_nz <<endl;
   for( int g = 0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
// Chain rule of interpolation relation, multiplication by constant matrix, since interpolation
// is a linear operator.
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
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
            dfs[3*ind]   = grhop[indm];
	    dfs[3*ind+1] = gmup[indm];
	    dfs[3*ind+2] = glambdap[indm];
	    ind++;
	    //            float_sw4 rho=xms[3*ind];
	    //	    float_sw4 cs=xms[3*ind+1];
	    //	    float_sw4 cp=xms[3*ind+2];
	    //            dfs[3*ind]   = cs*cs*gmup[indm]+(cp*cp-2*cs*cs)*glambdap[indm]+grhop[indm];
	    //	    dfs[3*ind+1] = 2*rho*cs*gmup[indm]-4*rho*cs*glambdap[indm];
	    //	    dfs[3*ind+2] = 2*rho*cp*glambdap[indm];
	    //	    ind++;
	 }
}

//-----------------------------------------------------------------------
ssize_t MaterialParCartesianVels::parameter_index( int ip, int jp, int kp, int grid,
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
ssize_t MaterialParCartesianVels::local_index( size_t ind_global )
{
   return -1;
}

//-----------------------------------------------------------------------
void MaterialParCartesianVels::gradient_transformation( std::vector<Sarray>& a_rho,
							std::vector<Sarray>& a_mu,
							std::vector<Sarray>& a_lambda,
							std::vector<Sarray>& a_gradrho,
							std::vector<Sarray>& a_gradmu,
							std::vector<Sarray>& a_gradlambda )
{
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
      m_ew->transform_gradient( a_rho[g], a_mu[g], a_lambda[g], 
				a_gradrho[g], a_gradmu[g], a_gradlambda[g] );
}
//-----------------------------------------------------------------------
void MaterialParCartesianVels::set_scalefactors( int nmpars, double* sfs, double rho_ref,
	 double mu_ref, double lambda_ref, double vs_ref, double vp_ref )
{
   for( int i=0 ; i < nmpars ; i += 3 )
   {
      sfs[i]   = rho_ref;
      sfs[i+1] = vs_ref;
      sfs[i+2] = vp_ref;
   }
}
