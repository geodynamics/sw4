#include <fcntl.h>
#include <unistd.h>

#include "MaterialParCartesian.h"
#include "EW.h"
#include "MParGridFile.h"

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
      if( m_myrank == 0 )
	cout << "active region, index = " <<
	  m_ew->m_iStartActGlobal[g] << " " <<
	  m_ew->m_iEndActGlobal[g] << " " <<
	  m_ew->m_jStartActGlobal[g] << " " <<
	  m_ew->m_jEndActGlobal[g] << " " <<
	  m_ew->m_kStartActGlobal[g] << " " <<
	  m_ew->m_kEndActGlobal[g] << endl;

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

   //cout << "MaterialParCartesian::get_material: nx=" << m_nx << " ny=" << m_ny << " nz=" << m_nz << " nms=" << nms << " nmd=" << nmd << endl;

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
   else if( m_init == 4 )
   {
      MParGridFile mpfile( m_filename );
      mpfile.interpolate_to_other( xms, 1, m_nx, m_ny, m_nz, m_hx, m_hy, m_hz, m_xmin, m_ymin, m_zmin );
      if( !mpfile.is_update() )
	 subtract_base_mtrl( nms, xms );
   }
}

//-----------------------------------------------------------------------
void MaterialParCartesian::get_gradient( int nmd, double* xmd, int nms, double* xms,
					 double* dfs, double* dfm,
					 std::vector<Sarray>& a_rho,
					 std::vector<Sarray>& a_mu,
					 std::vector<Sarray>& a_lambda,
					 std::vector<Sarray>& a_gradrho,
					 std::vector<Sarray>& a_gradmu,
					 std::vector<Sarray>& a_gradlambda)
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
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
            dfs[3*ind]   = grhop[indm];
	    dfs[3*ind+1] = gmup[indm];
	    dfs[3*ind+2] = glambdap[indm];
	    ind++;
	 }
}


//-----------------------------------------------------------------------
void MaterialParCartesian::interpolate_pseudohessian( int nmpars, double* phs,
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
               phs[ind*3  ] = phgrid[g](1,ig,jg,kg);
               phs[ind*3+1] = phgrid[g](2,ig,jg,kg);
               phs[ind*3+2] = phgrid[g](3,ig,jg,kg);
            }
            else
            {
               phs[ind*3  ] = 0;
               phs[ind*3+1] = 0;
               phs[ind*3+2] = 0;
            }
            ind++;
         }
   int npts = m_nx*m_ny*m_nz;
   float_sw4* tmp =new float_sw4[3*npts];
   for( int i=0 ; i < 3*npts ; i++ )
      tmp[i] = phs[i];
   MPI_Allreduce( tmp, phs, 3*npts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   delete[] tmp;
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

//-----------------------------------------------------------------------
extern "C" {
   void dgesv_( int*, int*, double*, int*, int*, double*, int*, int* );
}

//-----------------------------------------------------------------------
void MaterialParCartesian::project_and_write( std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
					      std::vector<Sarray>& a_lambda, std::string fname )
{
  float_sw4 mu_min= 1.02e7;
  float_sw4 rho_min=1593;
  float_sw4 lambda_min=7.67e8;
  
  float_sw4* projarray = new float_sw4[m_nx*m_ny*m_nz];
  float_sw4* xms = new float_sw4[m_nx*m_ny*m_nz*3];
  // Project rho
  projectl2( a_rho, projarray );
  for( int i=0 ; i < m_nx*m_ny*m_nz ; i++ )
    projarray[i] = projarray[i] > rho_min ? projarray[i] : rho_min;
  for( int i=0 ; i < m_nx*m_ny*m_nz ; i++ )
    xms[3*i] = projarray[i];

  // Project mu
  projectl2( a_mu, projarray );
  for( int i=0 ; i < m_nx*m_ny*m_nz ; i++ )
    projarray[i] = projarray[i] > mu_min ? projarray[i] : mu_min;
  for( int i=0 ; i < m_nx*m_ny*m_nz ; i++ )
    xms[3*i+1] = projarray[i];

  // Project lambda
  projectl2( a_lambda, projarray );
  for( int i=0 ; i < m_nx*m_ny*m_nz ; i++ )
    projarray[i] = projarray[i] > lambda_min ? projarray[i] : lambda_min;

  for( int i=0 ; i < m_nx*m_ny*m_nz ; i++ )
    xms[3*i+2] = projarray[i];
  MParGridFile mpfile( xms, 1, m_nx, m_ny, m_nz, m_hx, m_hy, m_hz, m_xmin, m_ymin, m_zmin, 0 );

  //  for( int i=0 ; i < m_nx*m_ny*m_nz ; i++ )
  //     if( xms[3*i] < 0 || xms[3*i+1] < 0 || xms[3*i+2]<0 )
  //       cout << "0:negval " << i << " " << xms[3*i] <<" " << xms[3*i+1] << " " << xms[3*i+2] << endl;

  mpfile.write_mparcartfile( fname );

  delete[] projarray;
  delete[] xms;
}
//-----------------------------------------------------------------------
void MaterialParCartesian::projectl2( std::vector<Sarray>& mtrl, float_sw4* rhs )
{
   // Project input mtrl array onto my parameter grid.
   //
   // parameter grid is x_0, x_1,..,x_{nx}, x_{nx+1}, where x_1,..,x_{nx} carry degrees of freedom.
   // The dimensions are such that x_0=xmin, x_{nx+1}=xmax.
   //
   // x_j = j*hx+xmin, j=0,..,nx+1  --> hx=(xmax-xmin)/(nx+1)
   //
   // In the z-direction, z_0 carries a degree of freedom, we set
   //   z_k = k*hz+zmin-hz, k=1,..,nz are still degrees of freedom, k=0 and k=nz+1 are 
   //  m_zmin is modified so that it really stores zmin-hz
   //
   //
   int nproc;
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
   CHECK_INPUT(nproc == 1,"ERROR: Projectl2 only implemented for single processor");
      
   // Initialize:
   int Ntot = m_nx*m_ny*m_nz;
   //   double* rhs = new double[Ntot];
   double* mat = new double[Ntot*Ntot];
   for( int i=0 ; i < Ntot ; i++ )
   {
      rhs[i] = 0;
      for( int j=0 ; j < Ntot ; j++ )
	 mat[j+Ntot*i]=0;
   }
   double ihx = 1.0/m_hx;
   double ihy = 1.0/m_hy;
   double ihz = 1.0/m_hz;
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      double h=m_ew->mGridSize[g];
      for( int k=mtrl[g].m_kb ; k<= mtrl[g].m_ke ; k++ )
	 for( int j=mtrl[g].m_jb ; j<= mtrl[g].m_je ; j++ )
	    for( int i=mtrl[g].m_ib ; i<= mtrl[g].m_ie ; i++ )
	    {
	       double x=(i-1)*h, y=(j-1)*h, z=(k-1)*h;
	       if( m_xmin+m_hx <= x && x <= m_xmin+m_nx*m_hx &&
		   m_ymin+m_hy <= y && y <= m_ymin+m_ny*m_hy &&
		   m_zmin+m_hz <= z && z <= m_zmin+m_nz*m_hz )
	       for( int m3=1 ; m3 <= m_nz ;m3++)
		  for( int m2=1 ; m2 <= m_ny ;m2++)
		     for( int m1=1 ; m1 <= m_nx ;m1++)
		     {
			double xpar1=m1*m_hx+m_xmin;
			double ypar1=m2*m_hy+m_ymin;
			double zpar1=m3*m_hz+m_zmin;
			if( std::abs(xpar1-x)< m_hx && std::abs(ypar1-y)<m_hy && std::abs(zpar1-z)<m_hz )
			{
			   int m=m1-1+m_nx*(m2-1)+m_nx*m_ny*(m3-1);
			   double phim;
			   if( xpar1-m_hx < x && x <= xpar1 )
			      phim = (x-(xpar1-m_hx))*ihx;
			   else
			      phim = (xpar1+m_hx-x)*ihx;
			   if( ypar1-m_hy < y && y <= ypar1 )
			      phim *= (y-(ypar1-m_hy))*ihy;
			   else
			      phim *= (ypar1+m_hy-y)*ihy;
			   if( zpar1-m_hz < z && z <= zpar1 )
			      phim *= (z-(zpar1-m_hz))*ihz;
			   else
			      phim *= (zpar1+m_hz-z)*ihz;
			   rhs[m] += mtrl[g](i,j,k)*phim;

			   int m3l=m3-1>=1   ? m3-1:1;
			   int m3u=m3+1<=m_nz? m3+1:m_nz;
			   int m2l=m2-1>=1   ? m2-1:1;
			   int m2u=m2+1<=m_ny? m2+1:m_ny;
			   int m1l=m1-1>=1   ? m1-1:1;
			   int m1u=m1+1<=m_nx? m1+1:m_nx;
			   for( int l3=m3l ; l3 <= m3u ; l3++ )
			      for( int l2=m2l ; l2 <= m2u ; l2++ )
				 for( int l1=m1l ; l1 <= m1u ; l1++ )
				 {
				    double xpar2=l1*m_hx+m_xmin;
				    double ypar2=l2*m_hy+m_ymin;
				    double zpar2=l3*m_hz+m_zmin;
				    if(  std::abs(xpar2-x)< m_hx && std::abs(ypar2-y)<m_hy && std::abs(zpar2-z)<m_hz )
				    {
				       int l=l1-1+m_nx*(l2-1)+m_nx*m_ny*(l3-1);
				       double phil;
				       if( xpar2-m_hx < x && x <= xpar2 )
					  phil = (x-(xpar2-m_hx))*ihx;
				       else
					  phil = (xpar2+m_hx-x)*ihx;
				       if( ypar2-m_hy < y && y <= ypar2 )
					  phil *= (y-(ypar2-m_hy))*ihy;
				       else
					  phil *= (ypar2+m_hy-y)*ihy;
				       if( zpar2-m_hz < z && z <= zpar2 )
					  phil *= (z-(zpar2-m_hz))*ihz;
				       else
					  phil *= (zpar2+m_hz-z)*ihz;
				       mat[m+Ntot*l] += phim*phil;
				    }
				 }
			}
		     }
	    }
   }	     
   //

   // TEST interpolate constant material:
   //   for( int i=0 ; i < Ntot ; i++ )
   //   {
   //      double rowsum=0;
   //      for( int j=0 ; j <Ntot ; j++ )
   //	 rowsum += mat[i+Ntot*j];
   //      if( std::abs(rhs[i]/rowsum-mtrl[0](1,1,1))>1e-3)
   //	  cout << "i= "<< i << " matsum = " << rowsum << " rhs= " << rhs[i] << endl;
   //   }

   // Solve system of linear equations
   int one=1, info=0;
   int* ipiv=new int[Ntot];
   dgesv_( &Ntot, &one, mat, &Ntot, ipiv, rhs, &Ntot, &info );
   if( info != 0 )
      cout <<"ERROR: in MaterialParCartesian::projectl2,  info = " << info << " from dgesv " << endl;
   delete[] ipiv;
   delete[] mat;

   //   // Output solution in 3D format:
   //   int fd=open(fname,O_WRONLY|O_TRUNC|O_CREAT, 0660);
   //
   //   // 1.header, nx,ny,nz,hx,hy,hz,xmin,ymin,zmin
   //   size_t nr=write(fd,&m_nx,sizeof(int));
   //   if( nr != sizeof(int))
   //      cout <<"ERROR writing m_nx, " << nr << " bytes written" << endl;
   //   nr=write(fd,&m_ny,sizeof(int));
   //   if( nr != sizeof(int))
   //      cout <<"ERROR writing m_ny, " << nr << " bytes written" << endl;
   //   nr=write(fd,&m_nz,sizeof(int));
   //   if( nr != sizeof(int))
   //      cout <<"ERROR writing m_nz, " << nr << " bytes written" << endl;
   //   nr=write(fd,&m_hx,sizeof(double));
   //   if( nr != sizeof(double))
   //      cout <<"ERROR writing m_hx, " << nr << " bytes written" << endl;
   //   nr=write(fd,&m_hy,sizeof(double));
   //   if( nr != sizeof(double))
   //      cout <<"ERROR writing m_hy, " << nr << " bytes written" << endl;
   //   nr=write(fd,&m_hz,sizeof(double));
   //   if( nr != sizeof(double))
   //      cout <<"ERROR writing m_hz, " << nr << " bytes written" << endl;
   //   nr=write(fd,&m_xmin,sizeof(double));
   //   if( nr != sizeof(double))
   //      cout <<"ERROR writing m_xmin, " << nr << " bytes written" << endl;
   //   nr=write(fd,&m_ymin,sizeof(double));
   //   if( nr != sizeof(double))
   //      cout <<"ERROR writing m_ymin, " << nr << " bytes written" << endl;
   //   nr=write(fd,&m_zmin,sizeof(double));
   //   if( nr != sizeof(double))
   //      cout <<"ERROR writing m_zmin, " << nr << " bytes written" << endl;
   //
   //   // 2.Projected material field
   //   nr=write(fd,rhs,Ntot*sizeof(double));
   //   if( nr != Ntot*sizeof(double))
   //      cout <<"ERROR writing rhs, " << nr << " bytes written" << endl;
   //
   //   close(fd);
   //   delete[] rhs;
}

//-----------------------------------------------------------------------
void MaterialParCartesian::subtract_base_mtrl( int nms, double* xms )
{
   // Assume xms are given as full material, interpolate and subtract the
   // base material to get xms as an update.

   Sarray rhobase(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray mubase(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray labase(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   m_ew->interpolate_base_to_coarse( m_nx, m_ny, m_nz, m_xmin, m_ymin, m_zmin, m_hx, m_hy, m_hz,
				    rhobase, mubase, labase );
   size_t ind = 0;
   double* rhop = rhobase.c_ptr();
   double* mup  = mubase.c_ptr();
   double* lap  = labase.c_ptr();
   for( int k=1 ; k <= m_nz ; k++ )
      for( int j=1 ; j <= m_ny ; j++ )
	 for( int i=1 ; i <= m_nx ; i++ )
	 {
            size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
	    xms[3*ind  ] -= rhop[indm];
	    xms[3*ind+1] -= mup[indm];
	    xms[3*ind+2] -= lap[indm];
	    ind++;
	 }
   
}

