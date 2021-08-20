#include "MaterialParCartesianVsVp.h"
#include "EW.h"
#include "MParGridFile.h"
//-----------------------------------------------------------------------
//  Parameterize the material on a Cartesian coarse grid, with
//  (cs, cp)-update on a fixed material as the parameters. 
//  
//  cs =x[2*i  ]
//  cp =x[2*i+1]
//
//-----------------------------------------------------------------------
MaterialParCartesianVsVp::MaterialParCartesianVsVp( EW* a_ew, int nx, int ny, int nz, int init, char* fname )
   : MaterialParameterization( a_ew, fname )
{
   //  VERIFY2( nx > 1 && ny > 1 && nz > 1, "MaterialParCartesianVsVp: The grid need at least two ponts in each direction")
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
      if( m_ew->topographyExists() && g >= m_ew->mNumberOfCartesianGrids )
      {
         zmin = 1e38;
	 for( int j= m_ew->m_jStartAct[g] ; j <= m_ew->m_jEndAct[g] ; j++ )
	    for( int i= m_ew->m_iStartAct[g] ; i <= m_ew->m_iEndAct[g] ; i++ )
	       if( m_ew->mZ[g](i,j,1) < zmin )
		  zmin = m_ew->mZ[g](i,j,1);
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
   m_nms = 2*nx*ny*nz;
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
void MaterialParCartesianVsVp::get_material( int nmd, double* xmd, int nms,
					     double* xms, vector<Sarray>& a_rho,
					     vector<Sarray>& a_mu, vector<Sarray>& a_lambda,
                    float_sw4 vp_min, float_sw4 vp_max, float_sw4 vs_min, float_sw4 vs_max,int wave_mode)
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

   //cout << "MaterialParCartesianVsVp::get_material: nx=" << m_nx << " ny=" << m_ny << " nz=" << m_nz << " nms=" << nms << " nmd=" << nmd << endl;

   for( int k=1 ; k <= m_nz ; k++ )
   for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
      size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
	    //            rhop[indm] = xms[3*ind];
       rhop[indm] = 0;
	    csp[indm]  = xms[2*ind ];
	    cpp[indm]  = xms[2*ind+1];
	    //	    mup[indm]    = rhop[indm]*csp[indm]*csp[indm];
	    //	    lambdap[indm]= rhop[indm]*cpp[indm]*cpp[indm]-2*mup[indm];
	    ind++;
	 }
   for( int g = 0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      //      m_ew->interpolate( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz, m_rho, m_mu,
      //			 m_lambda, g, a_rho[g], a_mu[g], a_lambda[g], false );
      // interpolate init model perturbation to the ref model      cs/cp->mu/lambda
      m_ew->interpolate( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz, m_rho, m_cs,
			 m_cp, g, a_rho[g], a_mu[g], a_lambda[g], false ); 
      // add mu/lambda update to the base      
      m_ew->update_and_transform_material( g, a_rho[g], a_mu[g], a_lambda[g], vp_min, vp_max, vs_min, vs_max, wave_mode);
   }
}

//-----------------------------------------------------------------------
// Input material on grid(s) (a_rho,a_mu,a_lambda), output corresponding
// parameter vector (xmd,xms).
//-----------------------------------------------------------------------
void MaterialParCartesianVsVp::interpolate_parameters( int nmd, double* xmd, int nms,
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
   //   double rhmin=100000,csmin=100000,cpmin=100000;
   //   double rhmax=-100000,csmax=-100000,cpmax=-100000; 
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
	    //            xms[3*ind]   = rhop[indm];
            xms[2*ind  ] = csp[indm];
	    xms[2*ind+1] = cpp[indm];
	    //	    if( rhop[indm]<rhmin )
	    //	       rhmin = rhop[indm];
	    //	    if( csp[indm]<csmin )
	    //	       csmin = csp[indm];
	    //	    if( cpp[indm]<cpmin )
	    //	       cpmin = cpp[indm];
	    //	    if( rhop[indm]>rhmax )
	    //	       rhmax = rhop[indm];
	    //	    if( csp[indm]>csmax )
	    //	       csmax = csp[indm];
	    //	    if( cpp[indm]>cpmax )
	    //	       cpmax = cpp[indm];
	    ind++;
	 }
}

//-----------------------------------------------------------------------
void MaterialParCartesianVsVp::get_parameters( int nmd, double* xmd, int nms,
					   double* xms, std::vector<Sarray>& a_rho, 
					   std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )
{

   if(m_myrank==0) cout << ">>>>>>>>> get_parameters: offset from reference model option m_init=" << m_init << endl;

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
	       //	       xms[3*ind]   = rho-1;
	       xms[2*ind  ] = cs-2;
	       xms[2*ind+1] = cp-4;
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
   else if( m_init == 4 )
   {
      MParGridFile mpfile( m_filename );
      mpfile.interpolate_to_other( xms, 4, m_nx, m_ny, m_nz, m_hx, m_hy, m_hz, m_xmin, m_ymin, m_zmin );
      //      for( int i=0 ; i < nms ; i++ )
      //	if( isnan(xms[i]) )
      //	  cout << "1:NAN at " << i << endl;
      if( !mpfile.is_update() )
	 subtract_base_mtrl( nms, xms );
      //      for( int i=0 ; i < nms ; i++ )
      //	if( isnan(xms[i]) )
      //	  cout << "2:NAN at " << i << endl;

   }

}


//-----------------------------------------------------------------------
void MaterialParCartesianVsVp::get_gradient( int nmd, double* xmd, int nms, double* xms,
					 double* dfs, double* dfm,
					 std::vector<Sarray>& a_rho,
					 std::vector<Sarray>& a_mu,
					 std::vector<Sarray>& a_lambda,
					 std::vector<Sarray>& a_gradrho,
					 std::vector<Sarray>& a_gradmu,
					 std::vector<Sarray>& a_gradlambda)
{
   // Computes gradient with respect to the material parameterization from given
   // gradients with respect to the material at grid points.
   // It is assumed that transform_gradient has been called before this routine.

  //if(rank==0) a_gradlambda[0].save_to_disk("gradlambda.say"); 

   // transform grads from lambda/mu to vp/vs
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
      m_ew->transform_gradient( a_rho[g], a_mu[g], a_lambda[g], 
				a_gradrho[g], a_gradmu[g], a_gradlambda[g] );

   Sarray grho(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray gmu(0,m_nx+1,0,m_ny+1,0,m_nz+1);
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

 
   //   double* grhop=grho.c_ptr();
   double* gmup=gmu.c_ptr();
   double* glambdap=glambda.c_ptr();
   int npts = (m_nx+2)*(m_ny+2)*(m_nz+2);

   double* tmp = new double[npts];
   //   for( int i=0 ; i < npts ; i++ )
   //      tmp[i] = grhop[i];
   //   MPI_Allreduce( tmp, grhop, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for( int i=0 ; i < npts ; i++ )
      tmp[i] = gmup[i];
   MPI_Allreduce( tmp, gmup, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for( int i=0 ; i < npts ; i++ )
      tmp[i] = glambdap[i];
   MPI_Allreduce( tmp, glambdap, npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   delete[] tmp;

   //MPI_Barrier(MPI_COMM_WORLD);
   //glambda.gaussian_smooth(31, 5.);
   //gmu.gaussian_smooth(21, 3.);
   //MPI_Barrier(MPI_COMM_WORLD);

   //glambda.save_to_disk("glambda.say");
   //gmu.save_to_disk("gmu.say");

   size_t ind =0;
   for( int k=1 ; k <= m_nz ; k++ )
    for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
	    //            dfs[3*ind]   = grhop[indm];
	    dfs[2*ind  ] = gmup[indm];
	    dfs[2*ind+1] = glambdap[indm];
	    ind++;
	    //            float_sw4 rho=xms[3*ind];
	    //	    float_sw4 cs=xms[3*ind+1];
	    //	    float_sw4 cp=xms[3*ind+2];
	    //            dfs[3*ind]   = cs*cs*gmup[indm]+(cp*cp-2*cs*cs)*glambdap[indm]+grhop[indm];
	    //	    dfs[3*ind+1] = 2*rho*cs*gmup[indm]-4*rho*cs*glambdap[indm];
	    //	    dfs[3*ind+2] = 2*rho*cp*glambdap[indm];
	    //	    ind++;
	 }
      //std::cout << " nx=" << m_nx << " ny=" << m_ny << " nz=" << m_nz << " ind=" << ind << std::endl;

}





//-----------------------------------------------------------------------
void MaterialParCartesianVsVp::interpolate_pseudohessian( int nmpars, double* phs,
                                                          int nmpard, double* phm,
                                                          vector<Sarray>& phgrid )
{
   int ig, jg, kg, g;
   size_t ind=0;
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            float_sw4 x = m_xmin + i*m_hx;
	    float_sw4 y = m_ymin + j*m_hy;
	    float_sw4 z = m_zmin + k*m_hz;
            m_ew->computeNearestLowGridPoint( ig, jg, kg, g, x, y, z );
            if( m_ew->interior_point_in_proc( ig, jg, g) )
	    {
               phs[ind*2  ] = phgrid[g](2,ig,jg,kg);
               phs[ind*2+1] = phgrid[g](3,ig,jg,kg);
            }
            else
            {
               phs[ind*2  ] = 0;
               phs[ind*2+1] = 0;
            }
            ind++;
         }
   int npts = m_nx*m_ny*m_nz;
   float_sw4* tmp =new float_sw4[2*npts];
   for( int i=0 ; i < 2*npts ; i++ )
      tmp[i] = phs[i];
   MPI_Allreduce( tmp, phs, 2*npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   delete[] tmp;
}

//-----------------------------------------------------------------------
ssize_t MaterialParCartesianVsVp::parameter_index( int ip, int jp, int kp, int grid,
					    int var )
// Ignore grid.
{
   if( 1 <= ip && ip <= m_nx && 1 <= jp && jp <= m_ny && 1 <= kp && kp <= m_nz
       && (0 <= var) && (var <= 1 ) )
      return 2*(ip-1+static_cast<ssize_t>(m_nx)*(jp-1)+m_nx*m_ny*(kp-1))+var;
   else
      return -1;
}

//-----------------------------------------------------------------------
ssize_t MaterialParCartesianVsVp::local_index( size_t ind_global )
{
   return -1;
}

//-----------------------------------------------------------------------
//void MaterialParCartesianVsVp::gradient_transformation( std::vector<Sarray>& a_rho,
//							std::vector<Sarray>& a_mu,
//							std::vector<Sarray>& a_lambda,
//							std::vector<Sarray>& a_gradrho,
//							std::vector<Sarray>& a_gradmu,
//							std::vector<Sarray>& a_gradlambda )
//{
//   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
//      m_ew->transform_gradient( a_rho[g], a_mu[g], a_lambda[g], 
//				a_gradrho[g], a_gradmu[g], a_gradlambda[g] );
//}
//-----------------------------------------------------------------------
void MaterialParCartesianVsVp::set_scalefactors( int nmpars, double* sfs, double rho_ref,
	 double mu_ref, double lambda_ref, double vs_ref, double vp_ref )
{
   for( int i=0 ; i < nmpars ; i += 2 )
   {
      sfs[i  ] = vs_ref;
      sfs[i+1] = vp_ref;
   }
}

//-----------------------------------------------------------------------
void MaterialParCartesianVsVp::subtract_base_mtrl( int nms, double* xms )
{
   // Assume xms are given as full material, interpolate and subtract the
   // base material to get xms as an update.

   Sarray rhobase(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray csbase(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray cpbase(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   m_ew->interpolate_base_to_coarse_vel( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz,
				    rhobase, csbase, cpbase );
   size_t ind = 0;
   double* csp  = csbase.c_ptr();
   double* cpp  = cpbase.c_ptr();
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
	    xms[2*ind  ] -= csp[indm];
	    xms[2*ind+1] -= cpp[indm];
	    ind++;
	 }
   
}

void MaterialParCartesianVsVp::get_base_parameters( int nmd, double* xmd, int nms,
					   double* xms, std::vector<Sarray>& a_rho, 
					   std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )
{
   std::cout << "VsVp::get_base_parameters" << std::endl;
   
   interpolate_base_parameters( nmd, xmd, nms, xms, a_rho, a_mu, a_lambda );
   
}

void MaterialParCartesianVsVp::interpolate_base_parameters( int nmd, double* xmd, int nms,
						   double* xms, std::vector<Sarray>& a_rho, 
			 std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda )
{
   // Interpolates the difference a_rho-(m_ew->mRho), into local rho. i.e., m_rho = I(a_rho-(m_ew->mRho))
   // where a_rho,mRho are on the computational grid, rho on the parameter grid.
   //                             similarly for cs, cp
   // m_cs = I( a_cs-mCs)
   // m_cp = I( a_cp-mCp)
   //
   std::cout << "nx=" << m_nx << " ny=" << m_ny << " nz=" << m_nz  << std::endl;
   std::cout << "rho min=" << a_rho[0].minimum() << " max=" << a_rho[0].maximum() << std::endl;
   std::cout << "lambda min=" << a_lambda[0].minimum() << " max=" << a_lambda[0].maximum() << std::endl;
   std::cout << "mu min=" << a_mu[0].minimum() << " max=" << a_mu[0].maximum() << std::endl;

   m_ew->interpolate_base_to_coarse_vel( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz,
				    m_rho, m_cs, m_cp);

   std::cout << "ew->interpolate_base_to_coarse_vel done" << std::endl;

   //float_sw4* rhop=m_rho.c_ptr();
   float_sw4* csp =m_cs.c_ptr();
   float_sw4* cpp =m_cp.c_ptr();
   size_t ind =0;
   double csmin=100000,cpmin=100000;
   double csmax=-100000,cpmax=-100000;

   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	      for( int i=1 ; i <= m_nx ; i++ )
	      {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
            xms[2*ind] = csp[indm];
	         xms[2*ind+1] = cpp[indm];
            
            if( csp[indm]<csmin )
               csmin = csp[indm];
            if( cpp[indm]<cpmin )
               cpmin = cpp[indm];
            
            if( csp[indm]>csmax )
               csmax = csp[indm];
            if( cpp[indm]>cpmax )
               cpmax = cpp[indm];
            ind++;
	    }

    std::cout << "cpmin=" << cpmin << " cpmax=" << cpmax << " csmin=" << csmin << " csmax=" << csmax << std::endl;

}


//-----------------------------------------------------------------------
void MaterialParCartesianVsVp::smooth_gradient(double* dfs, std::vector<Sarray>& a_Rho, std::vector<Sarray>& a_Mu, std::vector<Sarray>& a_Lambda, float_sw4 freq, float_sw4 sz)
{
   Sarray gmu(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray glambda(0,m_nx+1,0,m_ny+1,0,m_nz+1);

   double* gmup=gmu.c_ptr();
   double* glambdap=glambda.c_ptr();

   int count=0;
   float_sw4 cpavg=0.;
   float_sw4 csavg=0.;
   
   // find dominant wavelength
   int isz= (sz - m_zmin)/m_hz+0.5;
   std::cout << "smooth_gradient: sz=" << sz << " isz=" << isz << 
   " m_nx=" << m_nx << " m_ny=" << m_ny << " m_nz=" << m_nz << " freq=" << freq << std::endl;
   
   size_t ind=0;
   for( int k=1 ; k <= m_nz ; k++ )
    for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
        size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
	     gmup[indm] = dfs[2*ind];
	     glambdap[indm] = dfs[2*ind+1];
        ind++;
	 }

for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
   
      float_sw4* rhop = a_Rho[g].c_ptr();
      float_sw4* mup  = a_Mu[g].c_ptr();
      float_sw4* lap  = a_Lambda[g].c_ptr();

      for(int i=0; i< a_Rho[g].m_npts; i++) {
           cpavg += sqrt((lap[i]+2*mup[i])/rhop[i]);
           csavg += sqrt(mup[i]/rhop[i]);
           count++;
      }
      
   } // for all grids

      cpavg /= count;
      csavg /= count;

    std::cout << "solveTT cpavg=" <<  cpavg << " count=" << count << std::endl;
    std::cout << "solveTT csavg=" <<  csavg << " count=" << count << std::endl;

    
    // find averages from all procs
    double local[2]={cpavg, csavg};
    double global[2];

      MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      cpavg=global[0]/m_ew->no_of_procs();
      csavg=global[1]/m_ew->no_of_procs();
   
      //if( proc_zero() && verbose > 1 )
   
   int cp_len = cpavg / (freq*m_hx) /2 * 2+1;  // make it odd
   int cs_len = csavg / (freq*m_hx) /2 * 2+1;

    std::cout << "solveTT  reduced cpavg=" <<  cpavg << " cp_len=" << cp_len <<  std::endl;
    std::cout << "solveTT  reduced csavg=" <<  csavg << " cs_len=" << cs_len <<  std::endl;

  
   MPI_Barrier(MPI_COMM_WORLD);
   glambda.gaussian_smooth(cp_len, 5.);  // default 31 
   gmu.gaussian_smooth(cs_len, 3.);      // default 21
   MPI_Barrier(MPI_COMM_WORLD);

   //glambda.save_to_disk("glambda.say");
   //gmu.save_to_disk("gmu.say");
   ind =0;
   for( int k=1 ; k <= m_nz ; k++ )
    for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
      size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
	    dfs[2*ind  ] = gmup[indm];
	    dfs[2*ind+1] = glambdap[indm];
	    ind++;
	 }
std::cout << "reach end of smooth_gradient" << std::endl;

}