#include "TestGrid.h"
#include "TestTwilight.h"
#include "TestEcons.h"
#include "Farray.h"
#include "sw4.h"
#include "F77_FUNC.h"
#include <iomanip>

#include "CurvilinearInterface.h"
#include "Farray.h"

extern "C" {
   void F77_FUNC(dgetrf,DGETRF)(int*,  int*, double*, int*, int*, int*);
   void F77_FUNC(dgetrs,DGETRS)(char*, int*, int*, double*, int*, int*, double*, int*, int*);
}

CurvilinearInterface::CurvilinearInterface( int a_gc, EW* a_ew )
{
   m_tw   = 0;
   m_etest= 0;
   nghost = 5;
   gc = a_gc;
   gf = a_gc+1;
   hc = a_ew->mGridSize[gc];
   hf = a_ew->mGridSize[gf];
   Rop.define(-4,2);
   P.define(-1,2);
   ghcof.define(1,6);
   acof.define(1,6,1,8,1,8);
   bof.define(1,4,1,6);
   ux_cof.define(-2,2);
   acof_no_gp.define(1,6,1,8,1,8);
   Sb.define(0,5);
   sbop_no_gp.define(0,5);

   a.dim  = 3;
   a.n1_c = a_ew->m_global_nx[gc];
   a.n2_c = a_ew->m_global_ny[gc];
   a.n3_c = a_ew->m_global_nz[gc];
   a.n1_f = a_ew->m_global_nx[gf];
   a.n2_f = a_ew->m_global_ny[gf];
   a.n3_f = a_ew->m_global_nz[gf];   
   a.nrg  = nghost;

   a.h1_c = 1.0/(a.n1_c-1);
   a.h2_c = 1.0/(a.n2_c-1);
   a.h3_c = 1.0/(a.n3_c-1);
   a.h1_f = 1.0/(a.n1_f-1);
   a.h2_f = 1.0/(a.n2_f-1);
   a.h3_f = 1.0/(a.n3_f-1);

   float_sw4 bbox[6];
   a_ew->getGlobalBoundingBox( bbox );
   a.l1 = bbox[1];
   a.l2 = bbox[3];
   a.l3 = bbox[5]-bbox[4];

   define_coeffs();

   // Sb used for k=1 boundary here, was k=nk boundary in the test code.
   // Easier to change sign of coefficients than changing sign of operator everywhere in the code.
   for( int i=0 ; i < 6 ;i++ )
      Sb(i)=-Sb(i);

   m_test_grid = a_ew->create_gaussianHill();
   m_tw    = a_ew->create_twilight();
   //   m_etest = a_ew->create_energytest();
   m_ew=a_ew;
}

void negate_z( Sarray& z )
{
   double* zp=z.c_ptr();
   for( int i=0; i < z.m_npts ; i++ )
      zp[i] = -zp[i];
}

void negate_zcomponent( Sarray& U )
{
   for( int k=U.m_kb ; k <= U.m_ke ;k++)
     for( int j=U.m_jb ; j <= U.m_je ;j++)
       for( int i=U.m_ib ; i <= U.m_ie ;i++)
	 U(3,i,j,k)=-U(3,i,j,k);
}
void swap( int& a, int& b )
{
  int tmp=a;
  a=b;
  b=tmp;
}
void swap( double& a, double& b )
{
  double tmp=a;
  a=b;
  b=tmp;
}
void zero_out_unknowns( Sarray& U_f, Sarray& U_c, int n3_f )
{

  for( int j=U_c.m_jb ; j <= U_c.m_je ;j++)
    for( int i=U_c.m_ib ; i <= U_c.m_ie ;i++)
    {
      U_c(i,j,0)=0;
    }
  for( int j=U_f.m_jb ; j <= U_f.m_je ;j++)
    for( int i=U_f.m_ib ; i <= U_f.m_ie ;i++)
    {
      U_f(i,j,n3_f)=0;
    }
}
//-----------------------------------------------------------------------
void CurvilinearInterface::init_arrays( std::vector<Sarray>& a_mu, std::vector<Sarray>& a_lambda, 
                                        std::vector<Sarray>& a_rho, std::vector<Sarray>& a_metric, 
                                        std::vector<Sarray>& a_jac )
{
   int n1_c=a.n1_c, n2_c=a.n2_c, n3_c=a.n3_c;
   int n1_f=a.n1_f, n2_f=a.n2_f, n3_f=a.n3_f;

   rho_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,1,1);
   rho_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f,n3_f);
   mu_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   lambda_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);               
   mu_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   lambda_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);

   Jacobian_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   Jacobian_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);

   //   Jacobian_c.insert_intersection(a_jac[gc]);
   //   Jacobian_f.insert_intersection(a_jac[gf]);

   Sarray jac_c(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   Sarray met_c(4,1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   XI13_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   XI23_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   XI33_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   //   compute_metric( XI13_c, XI23_c, XI33_c, a_metric[gc], a_jac[gc], hc, n1_c, n2_c, n3_c );

   Sarray jac_f(  1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   Sarray met_f(4,1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   XI13_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   XI23_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   XI33_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   //   compute_metric( XI13_f, XI23_f, XI33_f, a_metric[gf], a_jac[gf], hf, n1_f, n2_f, n3_f );
   
   x_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   y_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   z_c.define(1-nghost,n1_c+nghost,1-nghost,n2_c+nghost,0,8);
   //   regenerate_grid( a_ew, gc, x_c, y_c, z_c );

   x_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   y_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   z_f.define(1-nghost,n1_f+nghost,1-nghost,n2_f+nghost,n3_f-7,n3_f);
   //   regenerate_grid( a_ew, gf, x_f, y_f, z_f );

   m_test_grid->generate_grid_and_met( m_ew, gf, x_f, y_f, z_f, jac_f, met_f );
   m_test_grid->generate_grid_and_met( m_ew, gc, x_c, y_c, z_c, jac_c, met_c );   
   //   negate_z( z_c );
   //   negate_z( z_f );
   
   //   m_test_grid->generate_grid_and_met( m_ew, gf, y_f, x_f, z_f, jac_f, met_f );
   //   m_test_grid->generate_grid_and_met( m_ew, gc, y_c, x_c, z_c, jac_c, met_c );   
   convert_metric( jac_f, met_f, XI13_f, XI23_f, XI33_f, Jacobian_f, hf, n1_f, n2_f, n3_f );
   convert_metric( jac_c, met_c, XI13_c, XI23_c, XI33_c, Jacobian_c, hc, n1_c, n2_c, n3_c );   

   if( m_tw != 0 )
   {
      m_tw->get_rho(rho_c,x_c,y_c,z_c);
      m_tw->get_rho(rho_f,x_f,y_f,z_f);
      m_tw->get_mula(mu_c,lambda_c,x_c,y_c,z_c);
      m_tw->get_mula(mu_f,lambda_f,x_f,y_f,z_f);
   }
   if( m_etest != 0 )
   {
      int sides[6]={1,1,1,1,0,0};
      rho_c.insert_intersection(a_rho[gc]);
      rho_f.insert_intersection(a_rho[gf]);
      mu_c.insert_intersection(a_mu[gc]);
      mu_f.insert_intersection(a_mu[gf]);
      lambda_c.insert_intersection(a_lambda[gc]);
      lambda_f.insert_intersection(a_lambda[gf]);
      int extra_ghost = nghost-m_ew->getNumberOfGhostPoints();
      if( extra_ghost > 0 )
      {
         m_etest->get_rhobnd(rho_c,extra_ghost,sides);
         m_etest->get_rhobnd(rho_f,extra_ghost,sides);
         m_etest->get_mulabnd(mu_c,lambda_c,extra_ghost,sides);
         m_etest->get_mulabnd(mu_f,lambda_f,extra_ghost,sides);
      }
   }
   //   rho_c *= Jacobian_c;
   //   rho_f *= Jacobian_f;
   scale_rho( rho_c, Jacobian_c );
   scale_rho( rho_f, Jacobian_f );

   // Transpose x <--> y, rho, jacobian, mu, lambda, xi13, xi23, xi33
   rho_c.transposeij();
   rho_f.transposeij();
   mu_c.transposeij();
   mu_f.transposeij( );
   lambda_c.transposeij();
   lambda_f.transposeij();
   Jacobian_c.transposeij();
   Jacobian_f.transposeij();
   XI13_c.transposeij();
   XI13_f.transposeij();
   XI23_c.transposeij();
   XI23_f.transposeij();
   XI33_c.transposeij();
   XI33_f.transposeij();
   swap( a.n1_c, a.n2_c );
   swap( a.h1_c, a.h2_c );
   swap( a.n1_f, a.n2_f );
   swap( a.h1_f, a.h2_f );
   swap( a.l1, a.l2 );
   n1_c=a.n1_c;
   n1_f=a.n1_f;
   n2_c=a.n2_c;
   n2_f=a.n2_f;
   //Swap XI13 and XI23
   // moved to convert_metric
   //   double* tmp= XI13_c.c_ptr();
   //   XI13_c.reference( XI23_c.c_ptr() );
   //   XI23_c.reference( tmp );
   //   tmp= XI13_f.c_ptr();
   //   XI13_f.reference( XI23_f.c_ptr() );
   //   XI23_f.reference(tmp);

   Mass_block.define(1,3,1,3,1,n1_c,1,n2_c);
   interface_block( Rop, ghcof, Sb, rho_c, lambda_c, rho_f, Jacobian_c, mu_c, P, XI13_c, XI23_c, XI33_c, Mass_block, a );
   //   x_c.save_to_disk("xc.bin");
   //   y_c.save_to_disk("yc.bin");
   //   z_c.save_to_disk("zc.bin");
   //   x_f.save_to_disk("xf.bin");
   //   y_f.save_to_disk("yf.bin");
   //   z_f.save_to_disk("zf.bin");
   // std::cout<<"POST INTERFACE BLOCK  ";
   // for( int ii=1 ; ii<=3;ii++)
   //     for( int jj=1 ; jj<=3;jj++)
   //        std::cout << Mass_block(ii,jj,12,13) << " ";
   // std::cout << std::endl;
   // std::cout<<"Metric  ";
   // std::cout << XI13_c(12,13,1) << " " << XI23_c(12,13,1) << " " << XI33_c(12,13,1)<< " " << Jacobian_c(12,13,1);
   // std::cout << std::endl;
   // std::cout<<"Metric  ";
   // std::cout << XI13_f(23,25,n3_f) << " " << XI23_f(23,25,n3_f) << " " << XI33_f(23,25,n3_f)<< " " << Jacobian_f(23,25,n3_f);
   // std::cout << std::endl;
   // std::cout<<"Mtrlc  ";
   // std::cout << mu_c(12,13,1) << " " << lambda_c(12,13,1) << " " << rho_c(12,13,1);
   // std::cout << std::endl;
   // std::cout<<"Mtrlf  ";
   // std::cout << mu_f(23,25,n3_f) << " " << lambda_f(23,25,n3_f) << " " << rho_f(23,25,n3_f);
   // std::cout << std::endl;


   int three = 3;
   int INFO=0;
   IPIV_block = new int[n1_c*n2_c*3];
   for (int j = 1; j <= n2_c; j++) {
      for (int i = 1; i <= n1_c; i++) {
         F77_FUNC(dgetrf,DGETRF)(&three, &three, &Mass_block(1, 1, i, j), &three,
                     &IPIV_block[0 + (i - 1) * 3 + (j - 1) * 3 * n1_c], &INFO);
         if (INFO != 0) {
            std::cerr << "LU Fails at (i,j) equals" << i << "," << j
                      << " INFO = " << INFO << " " << Mass_block(INFO, INFO, i, j)
                      << "\n";
            for (int l = 1; l <= 3; l++) {
               for (int m = 1; m <= 3; m++)
                  std::cerr << Mass_block(l, m, i, j) << ",";
               std::cerr << "\n";
            }
         }
      }
   }
// Communicate arrays here
}

void CurvilinearInterface::convert_metric( Sarray& jac, Sarray& met, Sarray& XI13, Sarray& XI23, 
                                           Sarray& XI33, Sarray& Jacobian, float_sw4 h, int n1, int n2, int n3 )
{
   float_sw4 jfact = static_cast<float_sw4>(n1-1)*(n2-1)*(n3-1);
   float_sw4 h1f=1/(h*(n3-1)), h3f=1.0/(n3-1);
   for( int k=jac.m_kb ; k <= jac.m_ke ; k++ )
      for( int j=jac.m_jb ; j <= jac.m_je ; j++ )
         for( int i=jac.m_ib ; i <= jac.m_ie ; i++ )
         {
            float_sw4 im1=1/met(1,i,j,k);
            XI13(i,j,k) = h1f*met(3,i,j,k)*im1;
            XI23(i,j,k) = h1f*met(2,i,j,k)*im1;
	    //            XI13(i,j,k) = -h1f*met(2,i,j,k)*im1;
	    //            XI23(i,j,k) = -h1f*met(3,i,j,k)*im1;
            XI33(i,j,k) = h3f*im1*im1;
            Jacobian(i,j,k) = jac(i,j,k)*jfact;
         }
}

void CurvilinearInterface::scale_rho( Sarray& rho, Sarray& Jacobian )
{
   for( int k=rho.m_kb ; k <= rho.m_ke ; k++ )
      for( int j=rho.m_jb ; j <= rho.m_je ; j++ )
         for( int i=rho.m_ib ; i <= rho.m_ie ; i++ )
            rho(i,j,k) *= Jacobian(i,j,k);
}

void CurvilinearInterface::impose_ic( std::vector<Sarray>& a_U, float_sw4 t )
{
   int n1_f = a.n1_f;
   int n2_f = a.n2_f;
   int n3_f = a.n3_f;
   int n1_c = a.n1_c;
   int n2_c = a.n2_c;
   int n3_c = a.n3_c;
   int nrg  = a.nrg;
   bool force_dirichlet = false, check_stress_cont=false;
   int fg=0;
   if( force_dirichlet )
      fg = 1;

   //   Sarray U_f(3,-4,n1_f+5,-4,n2_f+5,n3_f-7,n3_f+fg), U_c(3,-4,n1_c+5,-4,n2_c+5,0,8);
   Sarray U_f(3,-4,n2_f+5,-4,n1_f+5,n3_f-7,n3_f+fg), U_c(3,-4,n2_c+5,-4,n1_c+5,0,8);

//  1. copy   a_U into U_f and U_c
   U_f.insert_intersection(a_U[gf]);
   U_c.insert_intersection(a_U[gc]);

// 2a. Impose dirichlet conditions at ghost points
   int sides[6]={1,1,1,1,0,0};
   if( m_tw != 0 )
   {
      m_tw->get_ubnd( U_f, x_f, y_f, z_f, t, nrg+1, sides );
      m_tw->get_ubnd( U_c, x_c, y_c, z_c, t, nrg+1, sides );
      if( force_dirichlet )
      {
      // Debug
         sides[0]=sides[1]=sides[2]=sides[3]=sides[4]=0;
         sides[5]=1;
         m_tw->get_ubnd( U_f, x_f, y_f, z_f, t, nrg+1, sides );
         sides[0]=sides[1]=sides[2]=sides[3] = sides[5]=0;
         sides[4]=1;
         m_tw->get_ubnd( U_c, x_c, y_c, z_c, t, nrg+1, sides );
         a_U[gc].copy_kplane2(U_c,0); // have computed U_c:s ghost points
         a_U[gc].copy_kplane2(U_c,1); // have computed U_c:s ghost points
         a_U[gf].copy_kplane2(U_f,n3_f);    // .. and U_f:s interface points
         a_U[gf].copy_kplane2(U_f,n3_f+1);    // .. and U_f:s ghost point
         return;
      // End debug
      }
   }
   if( m_etest != 0 )
   {
      m_etest->get_ubnd( U_f, nrg+1, sides );
      m_etest->get_ubnd( U_c, nrg+1, sides );
   }

// 2b. communicate U_f and U_c, all k=planes.

// 2c. Transform components
   U_f.transposeij();
   U_c.transposeij();
   U_f.swap12();
   U_c.swap12();
   negate_zcomponent( U_f );
   negate_zcomponent( U_c );
   
   //   zero_out_unknowns( U_f, U_c, n3_f );
//   U_f.transform_coordsystem();
//   U_c.transform_coordsystem();   

// 3. Inject U_f := U_c on interface
//   std::cout << "uc, uf before " << U_c(1,9,36,1) << " " << U_f(1,17,71,n3_f)<<std::endl;
   injection( U_f, U_c, P, a );
   //   std::cout << "uc, uf after " << U_c(1,9,36,1) << " " << U_f(1,17,71,n3_f)<<std::endl;

// 4. Solve equation for stress continuity
   Farray Vass(1,n1_c*n2_c*3), LHS(1,n1_c*n2_c*3), residual(1,n1_c*n2_c*3);
   Farray Mass_f1(-2,n1_f+3,-2,n2_f+3);
   Farray lh_c(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1, 1, 1, 3);
   Farray lh_f(1 - nrg, n1_f + nrg, 1 - nrg, n2_f + nrg, n3_f, n3_f, 1, 3);

   // 4.a Form right hand side of equation
   //   cout << "Before interface_rhs"<<endl;
   interface_rhs(Vass, lh_c, lh_f, Jacobian_c, Jacobian_f, mu_c, mu_f, lambda_c,
                lambda_f, rho_c, rho_f, XI13_c, XI23_c, XI33_c, XI13_f, XI23_f,
                XI33_f, P, Sb, Rop, sbop_no_gp, acof_no_gp, U_c, U_f, Mass_f1,
                ux_cof, ghcof, acof, bof, a);


   // 4.b Left hand side gives the initial residual
   //   cout << "Before interface_lhs"<<endl;
   interface_lhs( LHS, Jacobian_c, mu_c, lambda_c, rho_c, rho_f,  XI13_c, XI23_c, XI33_c,
                  P, Sb, Rop, U_c, ghcof, a );
   for (int i = 1; i <= n1_c * n2_c * 3; i++) 
      residual(i) = Vass(i) - LHS(i);

   // 4.c Jacobi iteration 
   float_sw4 tol=1e-10;
   int iter=0;
   float_sw4 res=0;
   int INFO=0, three=3, one=1;
   char trans='N';
   //   cout << "Before while"<<endl;
   while(((res = residual.maxabs())>tol) && iter <= 50 )
   {
      iter++;
    //      std::cout << "Iteration " << iter << " " << std::setprecision(17) << res << "\n";
      for (int j = 1; j <= n2_c; j++) {
         for (int i = 1; i <= n1_c; i++) {
            F77_FUNC(dgetrs,DGETRS)(&trans, &three, &one, &Mass_block(1, 1, i, j), &three,
                    &IPIV_block[0 + (i - 1) * 3 + (j - 1) * 3 * n1_c],
                    &residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 1), &three, &INFO);
            if (INFO != 0) {
               std::cerr << "SOLVE Fails at (i,j) equals" << i << "," << j
                         << " INFO = " << INFO << " " << Mass_block(INFO, INFO, i, j)
                         << "\n";
               abort();
            }
            U_c(1, i, j, 0) +=
               residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 1);
            U_c(2, i, j, 0) +=
               residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 2);
            U_c(3, i, j, 0) +=
               residual((j - 1) * 3 * n1_c + 3 * (i - 1) + 3);
         }
      }
  // 4.d Communicate U_c here, only k=0 plane
      interface_lhs( LHS, Jacobian_c, mu_c, lambda_c, rho_c, rho_f,  XI13_c, XI23_c, XI33_c,
                     P, Sb, Rop, U_c, ghcof, a );
      for (int i = 1; i <= n1_c * n2_c * 3; i++) 
         residual(i) = Vass(i) - LHS(i);
   }
   if( res > tol )
      std::cout << "WARNING, no convergence in curvilinear interface, res = " << res << " tol= " << tol << std::endl;
   
   //   std::cout << "uc, uf before transform " << U_c(1,9,36,1) << " " << U_f(1,17,71,n3_f)<<std::endl;
   // 4.e Inverse transform the components.
   U_f.transposeij();
   U_c.transposeij();
   U_f.swap12();
   U_c.swap12();
   negate_zcomponent( U_f );
   negate_zcomponent( U_c );
   //   U_f.transform_coordsystem();
   //   U_c.transform_coordsystem();   

   //   std::cout << "uc, uf after transform " << U_c(2,36,9,1) << " " << U_f(2,71,17,n3_f)<<std::endl;
// 5. Copy U_c and U_f back to a_U, only k=0 for U_c and k=n3f for U_f.
   a_U[gc].copy_kplane2(U_c,0); // have computed U_c:s ghost points
   a_U[gf].copy_kplane2(U_f,n3_f);    // .. and U_f:s interface points
   //   std::cout << "uc, uf after copyplane " << a_U[gc](2,36,9,1) << " " << a_U[gf](2,71,17,n3_f)<<std::endl;
// 6. Check stress continuity
   if( check_stress_cont )
   {
      int nk=m_ew->m_kEndInt[gf];

      Sarray Bc(3,a_U[gc].m_ib,a_U[gc].m_ie,a_U[gc].m_jb,a_U[gc].m_je,1,1);
      Sarray B (3,a_U[gf].m_ib,a_U[gf].m_ie,a_U[gf].m_jb,a_U[gf].m_je,nk,nk);

      double* strx=new double[a_U[gf].m_ie-a_U[gf].m_ib+1];
      double* stry=new double[a_U[gf].m_je-a_U[gf].m_jb+1];
      for( int i=0; i < a_U[gf].m_ie-a_U[gf].m_ib+1 ;i++ )
         strx[i]=1;
      for( int i=0; i < a_U[gf].m_je-a_U[gf].m_jb+1 ;i++ )
         strx[i]=1;
      double sbop[6]={-3.0/12,-10/12.0,18.0/12.0,-6.0/12.0,1.0/12.0,0};
      double sbopng[6]={0,-25.0/12.0,4.0,-3.0,4.0/3, -1.0/4};

      m_ew->compute_icstresses_curv( a_U[gc], Bc, 1, m_ew->mMetric[gc], m_ew->mMu[gc], m_ew->mLambda[gc], strx, stry, sbop );
      m_ew->compute_icstresses_curv( a_U[gf], B , nk, m_ew->mMetric[gf], m_ew->mMu[gf], m_ew->mLambda[gf], strx, stry, sbopng );

      double bdiff[3] = {0,0,0}, bnorm[3]={0,0,0}; 
      int imx[3], jmx[3];
      for( int j=Bc.m_jb+3 ; j <= Bc.m_je-3 ; j++ )
	for( int i=Bc.m_ib+3 ; i <= Bc.m_ie-3 ; i++ )
	  for( int c=1 ; c <= 3 ; c++ )
	  {
	    if( abs(4*B(c,2*i-1,2*j-1,nk)-Bc(c,i,j,1)) > bdiff[c-1] )
	      {
	      bdiff[c-1] = abs(4*B(c,2*i-1,2*j-1,nk)-Bc(c,i,j,1));
	      imx[c-1]=i;
	      jmx[c-1]=j;
	      }
	    if( abs(B(c,2*i-1,2*j-1,nk)) > bnorm[c-1] )
	      bnorm[c-1] = abs(B(c,2*i-1,2*j-1,nk));

	  }
      std::cout << "Bstress norm(diff) = " << bdiff[0] << " " << bdiff[1] << " " << bdiff[2] << " max diff at (i,j) " << imx[0] << ","<<jmx[0] << "   " << imx[1] << "," << jmx[1] <<
	"     " << imx[2] << "," << jmx[2] <<  std::endl;
      std::cout << "Bstress norm = " << bnorm[0] << " " << bnorm[1] << " " << bnorm[2] << std::endl;
      bool printout=true;
      if( printout)
      {
	 int ic=71, jc=17; //ff  (121,121) (i,j)-grid (middle grid, i.e. coarse curvilinear)
	 if( n2_c == 61 )
	 {
	    ic=36;
	    jc=9;
	 }
         std::cout << "checking at (i,j) = " << ic << " " << jc << std::endl;
      //      int ic=18, jc=46; //f 61x61 (i,j)-grid (middle grid, i.e. coarse curvilinear)
	 std::cout << "Boundary stresses Bf = " << 4*B(1,2*ic-1,2*jc-1,nk) << " " << 4*B(2,2*ic-1,2*jc-1,nk) << " " << 4*B(3,2*ic-1,2*jc-1,nk) << std::endl;
	 std::cout << "Boundary stresses Bc = " << Bc(1,ic,jc,1) << " " << Bc(2,ic,jc,1) << " " << Bc(3,ic,jc,1) << std::endl;
	 std::cout << " difference           = " << 4*B(1,2*ic-1,2*jc-1,nk)-Bc(1,ic,jc,1) << " "
             << 4*B(2,2*ic-1,2*jc-1,nk)-Bc(2,ic,jc,1) << " " 
             << 4*B(3,2*ic-1,2*jc-1,nk)-Bc(3,ic,jc,1) << std::endl;
	 std::cout << "Difference in displacements " << a_U[gc](1,ic,jc,1)-a_U[gf](1,2*ic-1,2*jc-1,nk) << " "  <<
         a_U[gc](2,ic,jc,1)-a_U[gf](2,2*ic-1,2*jc-1,nk) << " " << a_U[gc](3,ic,jc,1)-a_U[gf](3,2*ic-1,2*jc-1,nk) << std::endl;
      //      B.save_to_disk("B.bin");
      //      Bc.save_to_disk("Bc.bin");
	 delete[] strx, stry;
      //      exit(0);
      }
   }
}

//-----------------------------------------------------------------------
void CurvilinearInterface::interface_block(Farray &Rop, Farray &ghcof, Farray &Sb, Sarray &rho_c,
                     Sarray &lambda_c, Sarray &rho_f, Sarray &Jacobian_c,
                     Sarray &mu_c, Farray &P, Sarray &XI13_c, Sarray &XI23_c,
                     Sarray &XI33_c, Farray &Mass_block, PackArgs &a) {
  float_sw4 int_cof;
  int i, j, k, l;
  auto h3_c = a.h3_c;
  auto h3_f = a.h3_f;

  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;

  auto n3_f = a.n3_f;
  //
  int_cof = 17.0 / 48.0 * h3_f * ghcof(1) / pow(h3_c, 2);
  // SIGN CHANGE
  int_cof = -int_cof;
  //
  Mass_block = 0.0;
  for (l = 1; l <= n2_c; l++) {
    for (k = 1; k <= n1_c; k++) {
      //
      for (j = -4; j <= 2; j += 2) {
        for (i = -4; i <= 2; i += 2) {
          // first set equation w.r.t the first component
          Mass_block(1, 1, k, l) =
              Mass_block(1, 1, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                       pow(XI13_c(k, l, 1), 2) +
                   mu_c(k, l, 1) * (pow(XI23_c(k, l, 1), 2) +
                                       pow(XI33_c(k, l, 1), 2))) /
                  rho_c(k, l, 1) * int_cof;
          // first set equation w.r.t the second component
          Mass_block(1, 2, k, l) =
              Mass_block(1, 2, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                  XI13_c(k, l, 1) * XI23_c(k, l, 1) / rho_c(k, l, 1) *
                  int_cof;
          // first set equation w.r.t the third component
          Mass_block(1, 3, k, l) =
              Mass_block(1, 3, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                  XI13_c(k, l, 1) * XI33_c(k, l, 1) / rho_c(k, l, 1) *
                  int_cof;
          // second set equation w.r.t the first component
          Mass_block(2, 1, k, l) =
              Mass_block(2, 1, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                  XI13_c(k, l, 1) * XI23_c(k, l, 1) / rho_c(k, l, 1) *
                  int_cof;
          // second set equation w.r.t the second component
          Mass_block(2, 2, k, l) =
              Mass_block(2, 2, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                       pow(XI23_c(k, l, 1), 2) +
                   mu_c(k, l, 1) * (pow(XI13_c(k, l, 1), 2) +
                                       pow(XI33_c(k, l, 1), 2))) /
                  rho_c(k, l, 1) * int_cof;
          // second set equation w.r.t the third component
          Mass_block(2, 3, k, l) =
              Mass_block(2, 3, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                  XI23_c(k, l, 1) * XI33_c(k, l, 1) / rho_c(k, l, 1) *
                  int_cof;
          // third set equation w.r.t the first component
          Mass_block(3, 1, k, l) =
              Mass_block(3, 1, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                  XI13_c(k, l, 1) * XI33_c(k, l, 1) / rho_c(k, l, 1) *
                  int_cof;
          // third set equation w.r.t the second component
          Mass_block(3, 2, k, l) =
              Mass_block(3, 2, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                  XI23_c(k, l, 1) * XI33_c(k, l, 1) / rho_c(k, l, 1) *
                  int_cof;
          // third set equation w.r.t the third component
          Mass_block(3, 3, k, l) =
              Mass_block(3, 3, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, 1) *
                  ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                       pow(XI33_c(k, l, 1), 2) +
                   mu_c(k, l, 1) * (pow(XI13_c(k, l, 1), 2) +
                                       pow(XI23_c(k, l, 1), 2))) /
                  rho_c(k, l, 1) * int_cof;
        }
      }
      //
      for (j = -4; j <= 2; j += 2) {
        // first set equation w.r.t the first component
        Mass_block(1, 1, k, l) =
            Mass_block(1, 1, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                     pow(XI13_c(k, l, 1), 2) +
                 mu_c(k, l, 1) * (pow(XI23_c(k, l, 1), 2) +
                                     pow(XI33_c(k, l, 1), 2))) /
                rho_c(k, l, 1) * int_cof;
        // first set equation w.r.t the second component
        Mass_block(1, 2, k, l) =
            Mass_block(1, 2, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
                XI23_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // first set equation w.r.t the third component
        Mass_block(1, 3, k, l) =
            Mass_block(1, 3, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
                XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // second set equation w.r.t the first component
        Mass_block(2, 1, k, l) =
            Mass_block(2, 1, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
                XI23_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // second set equation w.r.t the second component
        Mass_block(2, 2, k, l) =
            Mass_block(2, 2, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                     pow(XI23_c(k, l, 1), 2) +
                 mu_c(k, l, 1) * (pow(XI13_c(k, l, 1), 2) +
                                     pow(XI33_c(k, l, 1), 2))) /
                rho_c(k, l, 1) * int_cof;
        // second set equation w.r.t the third component
        Mass_block(2, 3, k, l) =
            Mass_block(2, 3, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
                XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // third set equation w.r.t the first component
        Mass_block(3, 1, k, l) =
            Mass_block(3, 1, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
                XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // third set equation w.r.t the second component
        Mass_block(3, 2, k, l) =
            Mass_block(3, 2, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
                XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // third set equation w.r.t the third component
        Mass_block(3, 3, k, l) =
            Mass_block(3, 3, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(-j / 2) *
                Jacobian_c(k, l, 1) *
                ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                     pow(XI33_c(k, l, 1), 2) +
                 mu_c(k, l, 1) * (pow(XI13_c(k, l, 1), 2) +
                                     pow(XI23_c(k, l, 1), 2))) /
                rho_c(k, l, 1) * int_cof;
      }
      //
      for (i = -4; i <= 2; i += 2) {
        // first set equation w.r.t the first component
        Mass_block(1, 1, k, l) =
            Mass_block(1, 1, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                     pow(XI13_c(k, l, 1), 2) +
                 mu_c(k, l, 1) * (pow(XI23_c(k, l, 1), 2) +
                                     pow(XI33_c(k, l, 1), 2))) /
                rho_c(k, l, 1) * int_cof;
        // first set equation w.r.t the second component
        Mass_block(1, 2, k, l) =
            Mass_block(1, 2, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
                XI23_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // first set equation w.r.t the third component
        Mass_block(1, 3, k, l) =
            Mass_block(1, 3, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
                XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // second set equation w.r.t the first component
        Mass_block(2, 1, k, l) =
            Mass_block(2, 1, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
                XI23_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // second set equation w.r.t the second component
        Mass_block(2, 2, k, l) =
            Mass_block(2, 2, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                     pow(XI23_c(k, l, 1), 2) +
                 mu_c(k, l, 1) * (pow(XI13_c(k, l, 1), 2) +
                                     pow(XI33_c(k, l, 1), 2))) /
                rho_c(k, l, 1) * int_cof;
        // second set equation w.r.t the third component
        Mass_block(2, 3, k, l) =
            Mass_block(2, 3, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
                XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // third set equation w.r.t the first component
        Mass_block(3, 1, k, l) =
            Mass_block(3, 1, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
                XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // third set equation w.r.t the second component
        Mass_block(3, 2, k, l) =
            Mass_block(3, 2, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
                XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
        // third set equation w.r.t the third component
        Mass_block(3, 3, k, l) =
            Mass_block(3, 3, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(-i / 2) *
                Jacobian_c(k, l, 1) *
                ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                     pow(XI33_c(k, l, 1), 2) +
                 mu_c(k, l, 1) * (pow(XI13_c(k, l, 1), 2) +
                                     pow(XI23_c(k, l, 1), 2))) /
                rho_c(k, l, 1) * int_cof;
      }
      //
      // first set equation w.r.t the first component
      Mass_block(1, 1, k, l) =
          Mass_block(1, 1, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI13_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
               (pow(XI23_c(k, l, 1), 2) + pow(XI33_c(k, l, 1), 2))) /
              rho_c(k, l, 1) * int_cof;
      // first set equation w.r.t the second component
      Mass_block(1, 2, k, l) =
          Mass_block(1, 2, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI23_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
      // first set equation w.r.t the third component
      Mass_block(1, 3, k, l) =
          Mass_block(1, 3, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
      // second set equation w.r.t the first component
      Mass_block(2, 1, k, l) =
          Mass_block(2, 1, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI23_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
      // second set equation w.r.t the second component
      Mass_block(2, 2, k, l) =
          Mass_block(2, 2, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI23_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI13_c(k, l, 1), 2) + pow(XI33_c(k, l, 1), 2))) /
              rho_c(k, l, 1) * int_cof;
      // second set equation w.r.t the third component
      Mass_block(2, 3, k, l) =
          Mass_block(2, 3, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
              XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
      // third set equation w.r.t the first component
      Mass_block(3, 1, k, l) =
          Mass_block(3, 1, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
      // third set equation w.r.t the second component
      Mass_block(3, 2, k, l) =
          Mass_block(3, 2, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
              XI33_c(k, l, 1) / rho_c(k, l, 1) * int_cof;
      // third set equation w.r.t the third component
      Mass_block(3, 3, k, l) =
          Mass_block(3, 3, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI33_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI13_c(k, l, 1), 2) + pow(XI23_c(k, l, 1), 2))) /
              rho_c(k, l, 1) * int_cof;
      // from the norm derivative

      // first set equation w.r.t the first component
      Mass_block(1, 1, k, l) =
          Mass_block(1, 1, k, l) -
          Sb(0) * Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI13_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
               (pow(XI23_c(k, l, 1), 2) + pow(XI33_c(k, l, 1), 2))) /
              h3_c;
      // first set equation w.r.t the second component
      Mass_block(1, 2, k, l) = Mass_block(1, 2, k, l) -
                               Sb(0) * Jacobian_c(k, l, 1) *
                                   (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                                   XI13_c(k, l, 1) * XI23_c(k, l, 1) /
                                   h3_c;
      // first set equation w.r.t the third component
      Mass_block(1, 3, k, l) = Mass_block(1, 3, k, l) -
                               Sb(0) * Jacobian_c(k, l, 1) *
                                   (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                                   XI13_c(k, l, 1) * XI33_c(k, l, 1) /
                                   h3_c;
      // second set equation w.r.t the first component
      Mass_block(2, 1, k, l) = Mass_block(2, 1, k, l) -
                               Sb(0) * Jacobian_c(k, l, 1) *
                                   (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                                   XI13_c(k, l, 1) * XI23_c(k, l, 1) /
                                   h3_c;
      // second set equation w.r.t the second component
      Mass_block(2, 2, k, l) =
          Mass_block(2, 2, k, l) -
          Sb(0) * Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI23_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI13_c(k, l, 1), 2) + pow(XI33_c(k, l, 1), 2))) /
              h3_c;
      // second set equation w.r.t the third component
      Mass_block(2, 3, k, l) = Mass_block(2, 3, k, l) -
                               Sb(0) * Jacobian_c(k, l, 1) *
                                   (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                                   XI23_c(k, l, 1) * XI33_c(k, l, 1) /
                                   h3_c;
      // third set equation w.r.t the first component
      Mass_block(3, 1, k, l) = Mass_block(3, 1, k, l) -
                               Sb(0) * Jacobian_c(k, l, 1) *
                                   (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                                   XI13_c(k, l, 1) * XI33_c(k, l, 1) /
                                   h3_c;
      // third set equation w.r.t the second component
      Mass_block(3, 2, k, l) = Mass_block(3, 2, k, l) -
                               Sb(0) * Jacobian_c(k, l, 1) *
                                   (lambda_c(k, l, 1) + mu_c(k, l, 1)) *
                                   XI23_c(k, l, 1) * XI33_c(k, l, 1) /
                                   h3_c;
      // third set equation w.r.t the third component
      Mass_block(3, 3, k, l) =
          Mass_block(3, 3, k, l) -
          Sb(0) * Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI33_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
               (pow(XI13_c(k, l, 1), 2) + pow(XI23_c(k, l, 1), 2))) /
              h3_c;
    }
  }
}

//-----------------------------------------------------------------------
void CurvilinearInterface::injection(Sarray &u_f, Sarray &u_c, Farray &P, PackArgs &a ) {
  int i, j, k;
  // Injection at the interface

  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;

  auto n3_f = a.n3_f;
  auto nrg  = a.nrg;
  for (int l = 1; l <= a.dim; l++) {
    for (j = 1; j <= n2_c; j++) {
      for (i = 1; i <= n1_c; i++) {
        u_f(l, 2 * i - 1, 2 * j - 1, n3_f) = u_c(l, i, j, 1);
        u_f(l, 2 * i, 2 * j - 1, n3_f) =
            P(-1) * u_c(l, i - 1, j, 1) + P(0) * u_c(l, i, j, 1) +
            P(1) * u_c(l, i + 1, j, 1) + P(2) * u_c(l, i + 2, j, 1);
        u_f(l, 2 * i - 1, 2 * j, n3_f) =
            P(-1) * u_c(l, i, j - 1, 1) + P(0) * u_c(l, i, j, 1) +
            P(1) * u_c(l, i, j + 1, 1) + P(2) * u_c(l, i, j + 2, 1);
        u_f(l, 2 * i, 2 * j, n3_f) =
            P(-1) * (P(-1) * u_c(l, i - 1, j - 1, 1) +
                     P(0) * u_c(l, i, j - 1, 1) +
                     P(1) * u_c(l, i + 1, j - 1, 1) +
                     P(2) * u_c(l, i + 2, j - 1, 1)) +
            P(0) * (P(-1) * u_c(l, i - 1, j, 1) + P(0) * u_c(l, i, j, 1) +
                    P(1) * u_c(l, i + 1, j, 1) +
                    P(2) * u_c(l, i + 2, j, 1)) +
            P(1) * (P(-1) * u_c(l, i - 1, j + 1, 1) +
                    P(0) * u_c(l, i, j + 1, 1) +
                    P(1) * u_c(l, i + 1, j + 1, 1) +
                    P(2) * u_c(l, i + 2, j + 1, 1)) +
            P(2) * (P(-1) * u_c(l, i - 1, j + 2, 1) +
                    P(0) * u_c(l, i, j + 2, 1) +
                    P(1) * u_c(l, i + 1, j + 2, 1) +
                    P(2) * u_c(l, i + 2, j + 2, 1));
      }
    }
  }
}


//-----------------------------------------------------------------------
void CurvilinearInterface::interface_rhs(Farray &Vass, Farray &lh_c, Farray &lh_f, Sarray &Jacobian_c,
                   Sarray &Jacobian_f, Sarray &mu_c, Sarray &mu_f,
                   Sarray &lambda_c, Sarray &lambda_f, Sarray &rho_c,
                   Sarray &rho_f, Sarray &XI13_c, Sarray &XI23_c,
                   Sarray &XI33_c, Sarray &XI13_f, Sarray &XI23_f,
                   Sarray &XI33_f, Farray &P, Farray &Sb, Farray &Rop,
                   Farray &sbop_no_gp, Farray &acof_no_gp, Sarray &u_c,
                   Sarray &u_f, Farray &Mass_f1, Farray &ux_cof, Farray &ghcof,
                   Farray &acof, Farray &bof, PackArgs &a ) {
  int i, j, k, k1, l, m;
  //
  Vass = 0.0;
  lh_c = 0.0;
  lh_f = 0.0;

  float_sw4 l1 = a.l1;
  float_sw4 l2 = a.l2;
  //  float_sw4 l3 = a.l3;

  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;

  auto n1_f = a.n1_f;
  auto n2_f = a.n2_f;
  auto n3_f = a.n3_f;

  auto h1_c = a.h1_c;
  auto h2_c = a.h2_c;
  auto h3_c = a.h3_c;

  auto h1_f = a.h1_f;
  auto h2_f = a.h2_f;
  auto h3_f = a.h3_f;

  auto nrg = a.nrg;

  //
  // term 1
  // Vass := Vass - Bc*uc
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = 1; j <= 4; j++) {
        // 33
        // first set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) +
            Jacobian_c(i, k, 1) *
                ((2.0 * mu_c(i, k, 1) + lambda_c(i, k, 1)) *
                     pow(XI13_c(i, k, 1), 2) +
                 mu_c(i, k, 1) * (pow(XI23_c(i, k, 1), 2) +
                                     pow(XI33_c(i, k, 1), 2))) *
                Sb(j) * u_c(1, i, k, j) / h3_c +
            Jacobian_c(i, k, 1) * (lambda_c(i, k, 1) + mu_c(i, k, 1)) *
                XI13_c(i, k, 1) * XI23_c(i, k, 1) * Sb(j) *
                u_c(2, i, k, j) / h3_c +
            Jacobian_c(i, k, 1) * (lambda_c(i, k, 1) + mu_c(i, k, 1)) *
                XI13_c(i, k, 1) * XI33_c(i, k, 1) * Sb(j) *
                u_c(3, i, k, j) / h3_c;
        // second set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) +
            Jacobian_c(i, k, 1) * (lambda_c(i, k, 1) + mu_c(i, k, 1)) *
                XI13_c(i, k, 1) * XI23_c(i, k, 1) * Sb(j) *
                u_c(1, i, k, j) / h3_c +
            Jacobian_c(i, k, 1) *
                ((2.0 * mu_c(i, k, 1) + lambda_c(i, k, 1)) *
                     pow(XI23_c(i, k, 1), 2) +
                 mu_c(i, k, 1) * (pow(XI13_c(i, k, 1), 2) +
                                     pow(XI33_c(i, k, 1), 2))) *
                Sb(j) * u_c(2, i, k, j) / h3_c +
            Jacobian_c(i, k, 1) * (lambda_c(i, k, 1) + mu_c(i, k, 1)) *
                XI23_c(i, k, 1) * XI33_c(i, k, 1) * Sb(j) *
                u_c(3, i, k, j) / h3_c;
        // third set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) +
            Jacobian_c(i, k, 1) * (lambda_c(i, k, 1) + mu_c(i, k, 1)) *
                XI13_c(i, k, 1) * XI33_c(i, k, 1) * Sb(j) *
                u_c(1, i, k, j) / h3_c +
            Jacobian_c(i, k, 1) * (lambda_c(i, k, 1) + mu_c(i, k, 1)) *
                XI23_c(i, k, 1) * XI33_c(i, k, 1) * Sb(j) *
                u_c(2, i, k, j) / h3_c +
            Jacobian_c(i, k, 1) *
                ((2.0 * mu_c(i, k, 1) + lambda_c(i, k, 1)) *
                     pow(XI33_c(i, k, 1), 2) +
                 mu_c(i, k, 1) * (pow(XI13_c(i, k, 1), 2) +
                                     pow(XI23_c(i, k, 1), 2))) *
                Sb(j) * u_c(3, i, k, j) / h3_c;
      }
    }
  }

  //
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -2; j <= 2; j++) {
        // 31
        // first set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) -
            Jacobian_c(i, k, 1) *
                (2.0 * mu_c(i, k, 1) + lambda_c(i, k, 1)) / l1 *
                XI13_c(i, k, 1) * ux_cof(j) * u_c(1, i + j, k, 1) / h1_c -
            Jacobian_c(i, k, 1) * mu_c(i, k, 1) / l1 *
                XI23_c(i, k, 1) * ux_cof(j) * u_c(2, i + j, k, 1) / h1_c -
            Jacobian_c(i, k, 1) * mu_c(i, k, 1) / l1 *
                XI33_c(i, k, 1) * ux_cof(j) * u_c(3, i + j, k, 1) / h1_c;
        // second set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) -
            Jacobian_c(i, k, 1) * lambda_c(i, k, 1) / l1 *
                XI23_c(i, k, 1) * ux_cof(j) * u_c(1, i + j, k, 1) / h1_c -
            Jacobian_c(i, k, 1) * mu_c(i, k, 1) / l1 *
                XI13_c(i, k, 1) * ux_cof(j) * u_c(2, i + j, k, 1) / h1_c;
        // third set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) -
            Jacobian_c(i, k, 1) * lambda_c(i, k, 1) / l1 *
                XI33_c(i, k, 1) * ux_cof(j) * u_c(1, i + j, k, 1) / h1_c -
            Jacobian_c(i, k, 1) * mu_c(i, k, 1) / l1 *
                XI13_c(i, k, 1) * ux_cof(j) * u_c(3, i + j, k, 1) / h1_c;
        // 32
        // first set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) -
            Jacobian_c(i, k, 1) * mu_c(i, k, 1) / l2 *
                XI23_c(i, k, 1) * ux_cof(j) * u_c(1, i, k + j, 1) / h2_c -
            Jacobian_c(i, k, 1) * lambda_c(i, k, 1) / l2 *
                XI13_c(i, k, 1) * ux_cof(j) * u_c(2, i, k + j, 1) / h2_c;
        // second set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) -
            Jacobian_c(i, k, 1) * mu_c(i, k, 1) / l2 *
                XI13_c(i, k, 1) * ux_cof(j) * u_c(1, i, k + j, 1) / h2_c -
            Jacobian_c(i, k, 1) *
                (2.0 * mu_c(i, k, 1) + lambda_c(i, k, 1)) / l2 *
                XI23_c(i, k, 1) * ux_cof(j) * u_c(2, i, k + j, 1) / h2_c -
            Jacobian_c(i, k, 1) * mu_c(i, k, 1) / l2 *
                XI33_c(i, k, 1) * ux_cof(j) * u_c(3, i, k + j, 1) / h2_c;
        // third set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) -
            Jacobian_c(i, k, 1) * lambda_c(i, k, 1) / l2 *
                XI33_c(i, k, 1) * ux_cof(j) * u_c(2, i, k + j, 1) / h2_c -
            Jacobian_c(i, k, 1) * mu_c(i, k, 1) / l2 *
                XI23_c(i, k, 1) * ux_cof(j) * u_c(3, i, k + j, 1) / h2_c;
      }
    }
  }

  // term 2
  //  lh_c  := lh_c + J*L_c*uc
  // interior
  for (j = -2; j <= n2_c + 3; j++) {
    for (i = -2; i <= n1_c + 3; i++) {
      // second derivative 11 & 22 & 12 & 21
      // first set
      lh_c(i, j, 1, 1) =
          lh_c(i, j, 1, 1) +
          ((-Jacobian_c(i - 2, j, 1) *
                (2.0 * mu_c(i - 2, j, 1) + lambda_c(i - 2, j, 1)) / 8.0 +
            Jacobian_c(i - 1, j, 1) *
                (2.0 * mu_c(i - 1, j, 1) + lambda_c(i - 1, j, 1)) / 6.0 -
            Jacobian_c(i, j, 1) *
            (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) / 8.0) *
               u_c(1, i - 2, j, 1) +
           (Jacobian_c(i - 2, j, 1) *
                (2.0 * mu_c(i - 2, j, 1) + lambda_c(i - 2, j, 1)) / 6.0 +
            Jacobian_c(i - 1, j, 1) *
                (2.0 * mu_c(i - 1, j, 1) + lambda_c(i - 1, j, 1)) / 2.0 +
            Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) / 2.0 +
            Jacobian_c(i + 1, j, 1) *
                (2.0 * mu_c(i + 1, j, 1) + lambda_c(i + 1, j, 1)) / 6.0) *
               u_c(1, i - 1, j, 1) +
           (-Jacobian_c(i - 2, j, 1) *
                (2.0 * mu_c(i - 2, j, 1) + lambda_c(i - 2, j, 1)) / 24.0 -
            Jacobian_c(i - 1, j, 1) *
                (2.0 * mu_c(i - 1, j, 1) + lambda_c(i - 1, j, 1)) * 5.0 /
                6.0 -
            Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) * 3.0 / 4.0 -
            Jacobian_c(i + 1, j, 1) *
                (2.0 * mu_c(i + 1, j, 1) + lambda_c(i + 1, j, 1)) * 5.0 /
                6.0 -
            Jacobian_c(i + 2, j, 1) *
                (2.0 * mu_c(i + 2, j, 1) + lambda_c(i + 2, j, 1)) /
                24.0) *
               u_c(1, i - 0, j, 1) +
           (Jacobian_c(i - 1, j, 1) *
                (2.0 * mu_c(i - 1, j, 1) + lambda_c(i - 1, j, 1)) / 6.0 +
            Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) / 2.0 +
            Jacobian_c(i + 1, j, 1) *
                (2.0 * mu_c(i + 1, j, 1) + lambda_c(i + 1, j, 1)) / 2.0 +
            Jacobian_c(i + 2, j, 1) *
                (2.0 * mu_c(i + 2, j, 1) + lambda_c(i + 2, j, 1)) / 6.0) *
               u_c(1, i + 1, j, 1) +
           (-Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) / 8.0 +
            Jacobian_c(i + 1, j, 1) *
                (2.0 * mu_c(i + 1, j, 1) + lambda_c(i + 1, j, 1)) / 6.0 -
            Jacobian_c(i + 2, j, 1) *
                (2.0 * mu_c(i + 2, j, 1) + lambda_c(i + 2, j, 1)) / 8.0) *
               u_c(1, i + 2, j, 1)) /
              (h1_c * h1_c) / (l1 * l1) +
          ((-Jacobian_c(i, j - 2, 1) * mu_c(i, j - 2, 1) / 8.0 +
            Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) / 6.0 -
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 8.0) *
               u_c(1, i, j - 2, 1) +
           (Jacobian_c(i, j - 2, 1) * mu_c(i, j - 2, 1) / 6.0 +
            Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) / 2.0 +
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 2.0 +
            Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) / 6.0) *
               u_c(1, i, j - 1, 1) +
           (-Jacobian_c(i, j - 2, 1) * mu_c(i, j - 2, 1) / 24.0 -
            Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) * 5.0 / 6.0 -
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) * 3.0 / 4.0 -
            Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) * 5.0 / 6.0 -
            Jacobian_c(i, j + 2, 1) * mu_c(i, j + 2, 1) / 24.0) *
               u_c(1, i, j - 0, 1) +
           (Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) / 6.0 +
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 2.0 +
            Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) / 2.0 +
            Jacobian_c(i, j + 2, 1) * mu_c(i, j + 2, 1) / 6.0) *
               u_c(1, i, j + 1, 1) +
           (-Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 8.0 +
            Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) / 6.0 -
            Jacobian_c(i, j + 2, 1) * mu_c(i, j + 2, 1) / 8.0) *
               u_c(1, i, j + 2, 1)) /
              (h2_c * h2_c) / (l2 * l2) +
          (Jacobian_c(i - 2, j, 1) * lambda_c(i - 2, j, 1) *
               (u_c(2, i - 2, j - 2, 1) / 12.0 -
                u_c(2, i - 2, j - 1, 1) * 2.0 / 3.0 +
                u_c(2, i - 2, j + 1, 1) * 2.0 / 3.0 -
                u_c(2, i - 2, j + 2, 1) / 12.0) /
               12.0 -
           Jacobian_c(i - 1, j, 1) * lambda_c(i - 1, j, 1) *
               (u_c(2, i - 1, j - 2, 1) / 12.0 -
                u_c(2, i - 1, j - 1, 1) * 2.0 / 3.0 +
                u_c(2, i - 1, j + 1, 1) * 2.0 / 3.0 -
                u_c(2, i - 1, j + 2, 1) / 12.0) *
               2.0 / 3.0 +
           Jacobian_c(i + 1, j, 1) * lambda_c(i + 1, j, 1) *
               (u_c(2, i + 1, j - 2, 1) / 12.0 -
                u_c(2, i + 1, j - 1, 1) * 2.0 / 3.0 +
                u_c(2, i + 1, j + 1, 1) * 2.0 / 3.0 -
                u_c(2, i + 1, j + 2, 1) / 12.0) *
               2.0 / 3.0 -
           Jacobian_c(i + 2, j, 1) * lambda_c(i + 2, j, 1) *
               (u_c(2, i + 2, j - 2, 1) / 12.0 -
                u_c(2, i + 2, j - 1, 1) * 2.0 / 3.0 +
                u_c(2, i + 2, j + 1, 1) * 2.0 / 3.0 -
                u_c(2, i + 2, j + 2, 1) / 12.0) /
               12.0) /
              l1 / l2 / h1_c / h2_c +
          (Jacobian_c(i, j - 2, 1) * mu_c(i, j - 2, 1) *
               (u_c(2, i - 2, j - 2, 1) / 12.0 -
                u_c(2, i - 1, j - 2, 1) * 2.0 / 3.0 +
                u_c(2, i + 1, j - 2, 1) * 2.0 / 3.0 -
                u_c(2, i + 2, j - 2, 1) / 12.0) /
               12.0 -
           Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) *
               (u_c(2, i - 2, j - 1, 1) / 12.0 -
                u_c(2, i - 1, j - 1, 1) * 2.0 / 3.0 +
                u_c(2, i + 1, j - 1, 1) * 2.0 / 3.0 -
                u_c(2, i + 2, j - 1, 1) / 12.0) *
               2.0 / 3.0 +
           Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) *
               (u_c(2, i - 2, j + 1, 1) / 12.0 -
                u_c(2, i - 1, j + 1, 1) * 2.0 / 3.0 +
                u_c(2, i + 1, j + 1, 1) * 2.0 / 3.0 -
                u_c(2, i + 2, j + 1, 1) / 12.0) *
               2.0 / 3.0 -
           Jacobian_c(i, j + 2, 1) * mu_c(i, j + 2, 1) *
               (u_c(2, i - 2, j + 2, 1) / 12.0 -
                u_c(2, i - 1, j + 2, 1) * 2.0 / 3.0 +
                u_c(2, i + 1, j + 2, 1) * 2.0 / 3.0 -
                u_c(2, i + 2, j + 2, 1) / 12.0) /
               12.0) /
              l1 / l2 / h1_c / h2_c;
      //
      // second set
      lh_c(i, j, 1, 2) =
          lh_c(i, j, 1, 2) +
          ((-Jacobian_c(i - 2, j, 1) * mu_c(i - 2, j, 1) / 8.0 +
            Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) / 6.0 -
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 8.0) *
               u_c(2, i - 2, j, 1) +
           (Jacobian_c(i - 2, j, 1) * mu_c(i - 2, j, 1) / 6.0 +
            Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) / 2.0 +
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 2.0 +
            Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) / 6.0) *
               u_c(2, i - 1, j, 1) +
           (-Jacobian_c(i - 2, j, 1) * mu_c(i - 2, j, 1) / 24.0 -
            Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) * 5.0 / 6.0 -
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) * 3.0 / 4.0 -
            Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) * 5.0 / 6.0 -
            Jacobian_c(i + 2, j, 1) * mu_c(i + 2, j, 1) / 24.0) *
               u_c(2, i - 0, j, 1) +
           (Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) / 6.0 +
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 2.0 +
            Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) / 2.0 +
            Jacobian_c(i + 2, j, 1) * mu_c(i + 2, j, 1) / 6.0) *
               u_c(2, i + 1, j, 1) +
           (-Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 8.0 +
            Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) / 6.0 -
            Jacobian_c(i + 2, j, 1) * mu_c(i + 2, j, 1) / 8.0) *
               u_c(2, i + 2, j, 1)) /
              pow(h1_c, 2) / pow(l1, 2) +
          ((-Jacobian_c(i, j - 2, 1) *
                (2.0 * mu_c(i, j - 2, 1) + lambda_c(i, j - 2, 1)) / 8.0 +
            Jacobian_c(i, j - 1, 1) *
                (2.0 * mu_c(i, j - 1, 1) + lambda_c(i, j - 1, 1)) / 6.0 -
            Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) / 8.0) *
               u_c(2, i, j - 2, 1) +
           (Jacobian_c(i, j - 2, 1) *
                (2.0 * mu_c(i, j - 2, 1) + lambda_c(i, j - 2, 1)) / 6.0 +
            Jacobian_c(i, j - 1, 1) *
                (2.0 * mu_c(i, j - 1, 1) + lambda_c(i, j - 1, 1)) / 2.0 +
            Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) / 2.0 +
            Jacobian_c(i, j + 1, 1) *
                (2.0 * mu_c(i, j + 1, 1) + lambda_c(i, j + 1, 1)) / 6.0) *
               u_c(2, i, j - 1, 1) +
           (-Jacobian_c(i, j - 2, 1) *
                (2.0 * mu_c(i, j - 2, 1) + lambda_c(i, j - 2, 1)) / 24.0 -
            Jacobian_c(i, j - 1, 1) *
                (2.0 * mu_c(i, j - 1, 1) + lambda_c(i, j - 1, 1)) * 5.0 /
                6.0 -
            Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) * 3.0 / 4.0 -
            Jacobian_c(i, j + 1, 1) *
                (2.0 * mu_c(i, j + 1, 1) + lambda_c(i, j + 1, 1)) * 5.0 /
                6.0 -
            Jacobian_c(i, j + 2, 1) *
                (2.0 * mu_c(i, j + 2, 1) + lambda_c(i, j + 2, 1)) /
                24.0) *
               u_c(2, i, j - 0, 1) +
           (Jacobian_c(i, j - 1, 1) *
                (2.0 * mu_c(i, j - 1, 1) + lambda_c(i, j - 1, 1)) / 6.0 +
            Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) / 2.0 +
            Jacobian_c(i, j + 1, 1) *
                (2.0 * mu_c(i, j + 1, 1) + lambda_c(i, j + 1, 1)) / 2.0 +
            Jacobian_c(i, j + 2, 1) *
                (2.0 * mu_c(i, j + 2, 1) + lambda_c(i, j + 2, 1)) / 6.0) *
               u_c(2, i, j + 1, 1) +
           (-Jacobian_c(i, j, 1) *
                (2.0 * mu_c(i, j, 1) + lambda_c(i, j, 1)) / 8.0 +
            Jacobian_c(i, j + 1, 1) *
                (2.0 * mu_c(i, j + 1, 1) + lambda_c(i, j + 1, 1)) / 6.0 -
            Jacobian_c(i, j + 2, 1) *
                (2.0 * mu_c(i, j + 2, 1) + lambda_c(i, j + 2, 1)) / 8.0) *
               u_c(2, i, j + 2, 1)) /
              pow(h2_c, 2) / pow(l2, 2) +
          (Jacobian_c(i - 2, j, 1) * mu_c(i - 2, j, 1) *
               (u_c(1, i - 2, j - 2, 1) / 12.0 -
                u_c(1, i - 2, j - 1, 1) * 2.0 / 3.0 +
                u_c(1, i - 2, j + 1, 1) * 2.0 / 3.0 -
                u_c(1, i - 2, j + 2, 1) / 12.0) /
               12.0 -
           Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) *
               (u_c(1, i - 1, j - 2, 1) / 12.0 -
                u_c(1, i - 1, j - 1, 1) * 2.0 / 3.0 +
                u_c(1, i - 1, j + 1, 1) * 2.0 / 3.0 -
                u_c(1, i - 1, j + 2, 1) / 12.0) *
               2.0 / 3.0 +
           Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) *
               (u_c(1, i + 1, j - 2, 1) / 12.0 -
                u_c(1, i + 1, j - 1, 1) * 2.0 / 3.0 +
                u_c(1, i + 1, j + 1, 1) * 2.0 / 3.0 -
                u_c(1, i + 1, j + 2, 1) / 12.0) *
               2.0 / 3.0 -
           Jacobian_c(i + 2, j, 1) * mu_c(i + 2, j, 1) *
               (u_c(1, i + 2, j - 2, 1) / 12.0 -
                u_c(1, i + 2, j - 1, 1) * 2.0 / 3.0 +
                u_c(1, i + 2, j + 1, 1) * 2.0 / 3.0 -
                u_c(1, i + 2, j + 2, 1) / 12.0) /
               12.0) /
              l1 / l2 / h1_c / h2_c +
          (Jacobian_c(i, j - 2, 1) * lambda_c(i, j - 2, 1) *
               (u_c(1, i - 2, j - 2, 1) / 12.0 -
                u_c(1, i - 1, j - 2, 1) * 2.0 / 3.0 +
                u_c(1, i + 1, j - 2, 1) * 2.0 / 3.0 -
                u_c(1, i + 2, j - 2, 1) / 12.0) /
               12.0 -
           Jacobian_c(i, j - 1, 1) * lambda_c(i, j - 1, 1) *
               (u_c(1, i - 2, j - 1, 1) / 12.0 -
                u_c(1, i - 1, j - 1, 1) * 2.0 / 3.0 +
                u_c(1, i + 1, j - 1, 1) * 2.0 / 3.0 -
                u_c(1, i + 2, j - 1, 1) / 12.0) *
               2.0 / 3.0 +
           Jacobian_c(i, j + 1, 1) * lambda_c(i, j + 1, 1) *
               (u_c(1, i - 2, j + 1, 1) / 12.0 -
                u_c(1, i - 1, j + 1, 1) * 2.0 / 3.0 +
                u_c(1, i + 1, j + 1, 1) * 2.0 / 3.0 -
                u_c(1, i + 2, j + 1, 1) / 12.0) *
               2.0 / 3.0 -
           Jacobian_c(i, j + 2, 1) * lambda_c(i, j + 2, 1) *
               (u_c(1, i - 2, j + 2, 1) / 12.0 -
                u_c(1, i - 1, j + 2, 1) * 2.0 / 3.0 +
                u_c(1, i + 1, j + 2, 1) * 2.0 / 3.0 -
                u_c(1, i + 2, j + 2, 1) / 12.0) /
               12.0) /
              l1 / l2 / h1_c / h2_c;
      // third set
      lh_c(i, j, 1, 3) =
          lh_c(i, j, 1, 3) +
          ((-Jacobian_c(i - 2, j, 1) * mu_c(i - 2, j, 1) / 8.0 +
            Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) / 6.0 -
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 8.0) *
               u_c(3, i - 2, j, 1) +
           (Jacobian_c(i - 2, j, 1) * mu_c(i - 2, j, 1) / 6.0 +
            Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) / 2.0 +
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 2.0 +
            Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) / 6.0) *
               u_c(3, i - 1, j, 1) +
           (-Jacobian_c(i - 2, j, 1) * mu_c(i - 2, j, 1) / 24.0 -
            Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) * 5.0 / 6.0 -
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) * 3.0 / 4.0 -
            Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) * 5.0 / 6.0 -
            Jacobian_c(i + 2, j, 1) * mu_c(i + 2, j, 1) / 24.0) *
               u_c(3, i - 0, j, 1) +
           (Jacobian_c(i - 1, j, 1) * mu_c(i - 1, j, 1) / 6.0 +
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 2.0 +
            Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) / 2.0 +
            Jacobian_c(i + 2, j, 1) * mu_c(i + 2, j, 1) / 6.0) *
               u_c(3, i + 1, j, 1) +
           (-Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 8.0 +
            Jacobian_c(i + 1, j, 1) * mu_c(i + 1, j, 1) / 6.0 -
            Jacobian_c(i + 2, j, 1) * mu_c(i + 2, j, 1) / 8.0) *
               u_c(3, i + 2, j, 1)) /
              pow(h1_c, 2) / pow(l1, 2) +
          ((-Jacobian_c(i, j - 2, 1) * mu_c(i, j - 2, 1) / 8.0 +
            Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) / 6.0 -
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 8.0) *
               u_c(3, i, j - 2, 1) +
           (Jacobian_c(i, j - 2, 1) * mu_c(i, j - 2, 1) / 6.0 +
            Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) / 2.0 +
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 2.0 +
            Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) / 6.0) *
               u_c(3, i, j - 1, 1) +
           (-Jacobian_c(i, j - 2, 1) * mu_c(i, j - 2, 1) / 24.0 -
            Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) * 5.0 / 6.0 -
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) * 3.0 / 4.0 -
            Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) * 5.0 / 6.0 -
            Jacobian_c(i, j + 2, 1) * mu_c(i, j + 2, 1) / 24.0) *
               u_c(3, i, j - 0, 1) +
           (Jacobian_c(i, j - 1, 1) * mu_c(i, j - 1, 1) / 6.0 +
            Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 2.0 +
            Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) / 2.0 +
            Jacobian_c(i, j + 2, 1) * mu_c(i, j + 2, 1) / 6.0) *
               u_c(3, i, j + 1, 1) +
           (-Jacobian_c(i, j, 1) * mu_c(i, j, 1) / 8.0 +
            Jacobian_c(i, j + 1, 1) * mu_c(i, j + 1, 1) / 6.0 -
            Jacobian_c(i, j + 2, 1) * mu_c(i, j + 2, 1) / 8.0) *
               u_c(3, i, j + 2, 1)) /
              pow(h2_c, 2) / pow(l2, 2);
    }
  }
  //
  for (j = -2; j <= n2_c + 3; j++) {
    for (k = -2; k <= n1_c + 3; k++) {
      for (k1 = 1; k1 <= 8; k1++) {
        for (m = 1; m <= 8; m++) {
          // second derivative 33
          // first set equation
          lh_c(k, j, 1, 1) =
              lh_c(k, j, 1, 1) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               ((2.0 * mu_c(k, j,  m) +
                 lambda_c(k, j,  m)) *
                    pow(XI13_c(k, j,  m), 2) +
                mu_c(k, j,  m) *
                    (pow(XI23_c(k, j,  m), 2) +
                     pow(XI33_c(k, j,  m), 2))) *
               u_c(1, k, j,  k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               (mu_c(k, j,  m) + lambda_c(k, j,  m)) *
               XI13_c(k, j,  m) * XI23_c(k, j,  m) *
               u_c(2, k, j,  k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               (mu_c(k, j,  m) + lambda_c(k, j,  m)) *
               XI13_c(k, j,  m) * XI33_c(k, j,  m) *
               u_c(3, k, j,  k1)) /
                  pow(h3_c, 2);
          // second set equation PROBLEM FIXED HERE
          lh_c(k, j, 1, 2) =
              lh_c(k, j, 1, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               (mu_c(k, j,  m) + lambda_c(k, j,  m)) *
               XI13_c(k, j,  m) * XI23_c(k, j,  m) *
               u_c(1, k, j,  k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               ((2.0 * mu_c(k, j,  m) +
                 lambda_c(k, j,  m)) *
                    pow(XI23_c(k, j,  m), 2) +
                mu_c(k, j,  m) *
                    (pow(XI13_c(k, j,  m), 2) +
                     pow(XI33_c(k, j,  m), 2))) *
               u_c(2, k, j,  k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               (mu_c(k, j,  m) + lambda_c(k, j,  m)) *
               XI23_c(k, j,  m) * XI33_c(k, j,  m) *
               u_c(3, k, j,  k1)) /
                  pow(h3_c, 2);
          // third set equation
          lh_c(k, j, 1, 3) =
              lh_c(k, j, 1, 3) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               (mu_c(k, j,  m) + lambda_c(k, j,  m)) *
               XI13_c(k, j,  m) * XI33_c(k, j,  m) *
               u_c(1, k, j,  k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               (mu_c(k, j,  m) + lambda_c(k, j,  m)) *
               XI23_c(k, j,  m) * XI33_c(k, j,  m) *
               u_c(2, k, j,  k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j,  m) *
               ((2.0 * mu_c(k, j,  m) +
                 lambda_c(k, j,  m)) *
                    pow(XI33_c(k, j,  m), 2) +
                mu_c(k, j,  m) *
                    (pow(XI13_c(k, j,  m), 2) +
                     pow(XI23_c(k, j,  m), 2))) *
               u_c(3, k, j,  k1)) /
                  pow(h3_c, 2);
        }
      }
    }
  }
  // ghost points
  for (j = -2; j <= n2_c + 3; j++) {
    for (k = -2; k <= 0; k++) {
      // first set equation
      lh_c(k, j, 1, 1) =
          lh_c(k, j, 1, 1) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI13_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI23_c(k, j, 1), 2) + pow(XI33_c(k, j, 1), 2)))) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI23_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, 1, 2) =
          lh_c(k, j, 1, 2) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI23_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI23_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI13_c(k, j, 1), 2) + pow(XI33_c(k, j, 1), 2)))) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI23_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, 1, 3) =
          lh_c(k, j, 1, 3) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI23_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2) +
         (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI33_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI13_c(k, j, 1), 2) + pow(XI23_c(k, j, 1), 2)))) /
              pow(h3_c, 2);
    }
  }

  //
  for (j = -2; j <= n2_c + 3; j++) {
    for (k = n1_c + 1; k <= n1_c + 3; k++) {
      // first set equation
      lh_c(k, j, 1, 1) =
          lh_c(k, j, 1, 1) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI13_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI23_c(k, j, 1), 2) + pow(XI33_c(k, j, 1), 2)))) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI23_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, 1, 2) =
          lh_c(k, j, 1, 2) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI23_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI23_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI13_c(k, j, 1), 2) + pow(XI33_c(k, j, 1), 2)))) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI23_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, 1, 3) =
          lh_c(k, j, 1, 3) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI23_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI33_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI13_c(k, j, 1), 2) + pow(XI23_c(k, j, 1), 2)))) /
              pow(h3_c, 2);
    }
  }
  //

  for (j = -2; j <= 0; j++) {
    for (k = 1; k <= n1_c; k++) {
      // first set equation
      lh_c(k, j, 1, 1) =
          lh_c(k, j, 1, 1) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI13_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI23_c(k, j, 1), 2) + pow(XI33_c(k, j, 1), 2)))) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI23_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, 1, 2) =
          lh_c(k, j, 1, 2) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI23_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI23_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI13_c(k, j, 1), 2) + pow(XI33_c(k, j, 1), 2)))) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI23_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, 1, 3) =
          lh_c(k, j, 1, 3) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI23_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI33_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI13_c(k, j, 1), 2) + pow(XI23_c(k, j, 1), 2)))) /
              pow(h3_c, 2);
    }
  }
  //
  for (j = n2_c + 1; j <= n2_c + 3; j++) {
    for (k = 1; k <= n1_c; k++) {
      // first set equation
      lh_c(k, j, 1, 1) =
          lh_c(k, j, 1, 1) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI13_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI23_c(k, j, 1), 2) + pow(XI33_c(k, j, 1), 2)))) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI23_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, 1, 2) =
          lh_c(k, j, 1, 2) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI23_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI23_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI13_c(k, j, 1), 2) + pow(XI33_c(k, j, 1), 2)))) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI23_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, 1, 3) =
          lh_c(k, j, 1, 3) +
          (u_c(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI13_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           (mu_c(k, j, 1) + lambda_c(k, j, 1)) * XI23_c(k, j, 1) *
           XI33_c(k, j, 1)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
           ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                pow(XI33_c(k, j, 1), 2) +
            mu_c(k, j, 1) *
                (pow(XI13_c(k, j, 1), 2) + pow(XI23_c(k, j, 1), 2)))) /
              pow(h3_c, 2);
    }
    }
  //
  for (i = -2; i <= n2_c + 3; i++) {
    for (j = -2; j <= n1_c + 3; j++) {
      for (k1 = 1; k1 <= 6; k1++) {
        // mixed derivative 13  23  31  32
        // first set equation
        lh_c(j, i, 1, 1) =
            lh_c(j, i, 1, 1) +
            (-Jacobian_c(j - 2, i, 1) *
                 (2.0 * mu_c(j - 2, i, 1) + lambda_c(j - 2, i, 1)) *
                 XI13_c(j - 2, i, 1) * (-bof(1, k1)) *
                 u_c(1, j - 2, i,  k1) / 12.0 +
             Jacobian_c(j - 1, i, 1) *
                 (2.0 * mu_c(j - 1, i, 1) + lambda_c(j - 1, i, 1)) *
                 XI13_c(j - 1, i, 1) * (-bof(1, k1)) *
                 u_c(1, j - 1, i,  k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, 1) *
                 (2.0 * mu_c(j + 1, i, 1) + lambda_c(j + 1, i, 1)) *
                 XI13_c(j + 1, i, 1) * (-bof(1, k1)) *
                 u_c(1, j + 1, i,  k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, 1) *
                 (2.0 * mu_c(j + 2, i, 1) + lambda_c(j + 2, i, 1)) *
                 XI13_c(j + 2, i, 1) * (-bof(1, k1)) *
                 u_c(1, j + 2, i,  k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j - 2, i, 1) * lambda_c(j - 2, i, 1) *
                 XI23_c(j - 2, i, 1) * (-bof(1, k1)) *
                 u_c(2, j - 2, i,  k1) / 12.0 +
             Jacobian_c(j - 1, i, 1) * lambda_c(j - 1, i, 1) *
                 XI23_c(j - 1, i, 1) * (-bof(1, k1)) *
                 u_c(2, j - 1, i,  k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, 1) * lambda_c(j + 1, i, 1) *
                 XI23_c(j + 1, i, 1) * (-bof(1, k1)) *
                 u_c(2, j + 1, i,  k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, 1) * lambda_c(j + 2, i, 1) *
                 XI23_c(j + 2, i, 1) * (-bof(1, k1)) *
                 u_c(2, j + 2, i,  k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j - 2, i, 1) * lambda_c(j - 2, i, 1) *
                 XI33_c(j - 2, i, 1) * (-bof(1, k1)) *
                 u_c(3, j - 2, i,  k1) / 12.0 +
             Jacobian_c(j - 1, i, 1) * lambda_c(j - 1, i, 1) *
                 XI33_c(j - 1, i, 1) * (-bof(1, k1)) *
                 u_c(3, j - 1, i,  k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, 1) * lambda_c(j + 1, i, 1) *
                 XI33_c(j + 1, i, 1) * (-bof(1, k1)) *
                 u_c(3, j + 1, i,  k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, 1) * lambda_c(j + 2, i, 1) *
                 XI33_c(j + 2, i, 1) * (-bof(1, k1)) *
                 u_c(3, j + 2, i,  k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j, i - 2, 1) * mu_c(j, i - 2, 1) *
                 XI23_c(j, i - 2, 1) * (-bof(1, k1)) *
                 u_c(1, j, i - 2,  k1) / 12.0 +
             Jacobian_c(j, i - 1, 1) * mu_c(j, i - 1, 1) *
                 XI23_c(j, i - 1, 1) * (-bof(1, k1)) *
                 u_c(1, j, i - 1,  k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, 1) * mu_c(j, i + 1, 1) *
                 XI23_c(j, i + 1, 1) * (-bof(1, k1)) *
                 u_c(1, j, i + 1,  k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, 1) * mu_c(j, i + 2, 1) *
                 XI23_c(j, i + 2, 1) * (-bof(1, k1)) *
                 u_c(1, j, i + 2,  k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-Jacobian_c(j, i - 2, 1) * mu_c(j, i - 2, 1) *
                 XI13_c(j, i - 2, 1) * (-bof(1, k1)) *
                 u_c(2, j, i - 2,  k1) / 12.0 +
             Jacobian_c(j, i - 1, 1) * mu_c(j, i - 1, 1) *
                 XI13_c(j, i - 1, 1) * (-bof(1, k1)) *
                 u_c(2, j, i - 1,  k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, 1) * mu_c(j, i + 1, 1) *
                 XI13_c(j, i + 1, 1) * (-bof(1, k1)) *
                 u_c(2, j, i + 1,  k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, 1) * mu_c(j, i + 2, 1) *
                 XI13_c(j, i + 2, 1) * (-bof(1, k1)) *
                 u_c(2, j, i + 2,  k1) / 12.0) /
                l2 / h2_c / h3_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             (2.0 * mu_c(j, i,  k1) + lambda_c(j, i,  k1)) *
             XI13_c(j, i,  k1) *
             (u_c(1, j - 2, i,  k1) / 12.0 -
              u_c(1, j - 1, i,  k1) * 2.0 / 3.0 +
              u_c(1, j + 1, i,  k1) * 2.0 / 3.0 -
              u_c(1, j + 2, i,  k1) / 12.0)) /
                l1 / h3_c / h1_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             mu_c(j, i,  k1) * XI23_c(j, i,  k1) *
             (u_c(2, j - 2, i,  k1) / 12.0 -
              u_c(2, j - 1, i,  k1) * 2.0 / 3.0 +
              u_c(2, j + 1, i,  k1) * 2.0 / 3.0 -
              u_c(2, j + 2, i,  k1) / 12.0)) /
                l1 / h3_c / h1_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             mu_c(j, i,  k1) * XI33_c(j, i,  k1) *
             (u_c(3, j - 2, i,  k1) / 12.0 -
              u_c(3, j - 1, i,  k1) * 2.0 / 3.0 +
              u_c(3, j + 1, i,  k1) * 2.0 / 3.0 -
              u_c(3, j + 2, i,  k1) / 12.0)) /
                l1 / h1_c / h3_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             mu_c(j, i,  k1) * XI23_c(j, i,  k1) *
             (u_c(1, j, i - 2,  k1) / 12.0 -
              u_c(1, j, i - 1,  k1) * 2.0 / 3.0 +
              u_c(1, j, i + 1,  k1) * 2.0 / 3.0 -
              u_c(1, j, i + 2,  k1) / 12.0)) /
                l2 / h3_c / h2_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             lambda_c(j, i,  k1) * XI13_c(j, i,  k1) *
             (u_c(2, j, i - 2,  k1) / 12.0 -
              u_c(2, j, i - 1,  k1) * 2.0 / 3.0 +
              u_c(2, j, i + 1,  k1) * 2.0 / 3.0 -
              u_c(2, j, i + 2,  k1) / 12.0)) /
                l2 / h3_c / h2_c;
        // second set equation
        lh_c(j, i, 1, 2) =
            lh_c(j, i, 1, 2) +
            (-Jacobian_c(j - 2, i, 1) * mu_c(j - 2, i, 1) *
                 XI23_c(j - 2, i, 1) * (-bof(1, k1)) *
                 u_c(1, j - 2, i,  k1) / 12.0 +
             Jacobian_c(j - 1, i, 1) * mu_c(j - 1, i, 1) *
                 XI23_c(j - 1, i, 1) * (-bof(1, k1)) *
                 u_c(1, j - 1, i,  k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, 1) * mu_c(j + 1, i, 1) *
                 XI23_c(j + 1, i, 1) * (-bof(1, k1)) *
                 u_c(1, j + 1, i,  k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, 1) * mu_c(j + 2, i, 1) *
                 XI23_c(j + 2, i, 1) * (-bof(1, k1)) *
                 u_c(1, j + 2, i,  k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j - 2, i, 1) * mu_c(j - 2, i, 1) *
                 XI13_c(j - 2, i, 1) * (-bof(1, k1)) *
                 u_c(2, j - 2, i,  k1) / 12.0 +
             Jacobian_c(j - 1, i, 1) * mu_c(j - 1, i, 1) *
                 XI13_c(j - 1, i, 1) * (-bof(1, k1)) *
                 u_c(2, j - 1, i,  k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, 1) * mu_c(j + 1, i, 1) *
                 XI13_c(j + 1, i, 1) * (-bof(1, k1)) *
                 u_c(2, j + 1, i,  k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, 1) * mu_c(j + 2, i, 1) *
                 XI13_c(j + 2, i, 1) * (-bof(1, k1)) *
                 u_c(2, j + 2, i,  k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j, i - 2, 1) * lambda_c(j, i - 2, 1) *
                 XI13_c(j, i - 2, 1) * (-bof(1, k1)) *
                 u_c(1, j, i - 2,  k1) / 12.0 +
             Jacobian_c(j, i - 1, 1) * lambda_c(j, i - 1, 1) *
                 XI13_c(j, i - 1, 1) * (-bof(1, k1)) *
                 u_c(1, j, i - 1,  k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, 1) * lambda_c(j, i + 1, 1) *
                 XI13_c(j, i + 1, 1) * (-bof(1, k1)) *
                 u_c(1, j, i + 1,  k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, 1) * lambda_c(j, i + 2, 1) *
                 XI13_c(j, i + 2, 1) * (-bof(1, k1)) *
                 u_c(1, j, i + 2,  k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-Jacobian_c(j, i - 2, 1) *
                 (2.0 * mu_c(j, i - 2, 1) + lambda_c(j, i - 2, 1)) *
                 XI23_c(j, i - 2, 1) * (-bof(1, k1)) *
                 u_c(2, j, i - 2,  k1) / 12.0 +
             Jacobian_c(j, i - 1, 1) *
                 (2.0 * mu_c(j, i - 1, 1) + lambda_c(j, i - 1, 1)) *
                 XI23_c(j, i - 1, 1) * (-bof(1, k1)) *
                 u_c(2, j, i - 1,  k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, 1) *
                 (2.0 * mu_c(j, i + 1, 1) + lambda_c(j, i + 1, 1)) *
                 XI23_c(j, i + 1, 1) * (-bof(1, k1)) *
                 u_c(2, j, i + 1,  k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, 1) *
                 (2.0 * mu_c(j, i + 2, 1) + lambda_c(j, i + 2, 1)) *
                 XI23_c(j, i + 2, 1) * (-bof(1, k1)) *
                 u_c(2, j, i + 2,  k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-Jacobian_c(j, i - 2, 1) * lambda_c(j, i - 2, 1) *
                 XI33_c(j, i - 2, 1) * (-bof(1, k1)) *
                 u_c(3, j, i - 2,  k1) / 12.0 +
             Jacobian_c(j, i - 1, 1) * lambda_c(j, i - 1, 1) *
                 XI33_c(j, i - 1, 1) * (-bof(1, k1)) *
                 u_c(3, j, i - 1,  k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, 1) * lambda_c(j, i + 1, 1) *
                 XI33_c(j, i + 1, 1) * (-bof(1, k1)) *
                 u_c(3, j, i + 1,  k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, 1) * lambda_c(j, i + 2, 1) *
                 XI33_c(j, i + 2, 1) * (-bof(1, k1)) *
                 u_c(3, j, i + 2,  k1) / 12.0) /
                l2 / h2_c / h3_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             lambda_c(j, i,  k1) * XI23_c(j, i,  k1) *
             (u_c(1, j - 2, i,  k1) / 12.0 -
              u_c(1, j - 1, i,  k1) * 2.0 / 3.0 +
              u_c(1, j + 1, i,  k1) * 2.0 / 3.0 -
              u_c(1, j + 2, i,  k1) / 12.0)) /
                l1 / h1_c / h3_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             mu_c(j, i,  k1) * XI13_c(j, i,  k1) *
             (u_c(2, j - 2, i,  k1) / 12.0 -
              u_c(2, j - 1, i,  k1) * 2.0 / 3.0 +
              u_c(2, j + 1, i,  k1) * 2.0 / 3.0 -
              u_c(2, j + 2, i,  k1) / 12.0)) /
                l1 / h1_c / h3_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             mu_c(j, i,  k1) * XI13_c(j, i,  k1) *
             (u_c(1, j, i - 2,  k1) / 12.0 -
              u_c(1, j, i - 1,  k1) * 2.0 / 3.0 +
              u_c(1, j, i + 1,  k1) * 2.0 / 3.0 -
              u_c(1, j, i + 2,  k1) / 12.0)) /
                l2 / h3_c / h2_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             (2.0 * mu_c(j, i,  k1) + lambda_c(j, i,  k1)) *
             XI23_c(j, i,  k1) *
             (u_c(2, j, i - 2,  k1) / 12.0 -
              u_c(2, j, i - 1,  k1) * 2.0 / 3.0 +
              u_c(2, j, i + 1,  k1) * 2.0 / 3.0 -
              u_c(2, j, i + 2,  k1) / 12.0)) /
                l2 / h3_c / h2_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             mu_c(j, i,  k1) * XI33_c(j, i,  k1) *
             (u_c(3, j, i - 2,  k1) / 12.0 -
              u_c(3, j, i - 1,  k1) * 2.0 / 3.0 +
              u_c(3, j, i + 1,  k1) * 2.0 / 3.0 -
              u_c(3, j, i + 2,  k1) / 12.0)) /
                l2 / h3_c / h2_c;
        // third set equation
        lh_c(j, i, 1, 3) =
            lh_c(j, i, 1, 3) +
            (-Jacobian_c(j - 2, i, 1) * mu_c(j - 2, i, 1) *
                 XI33_c(j - 2, i, 1) * (-bof(1, k1)) *
                 u_c(1, j - 2, i,  k1) / 12.0 +
             Jacobian_c(j - 1, i, 1) * mu_c(j - 1, i, 1) *
                 XI33_c(j - 1, i, 1) * (-bof(1, k1)) *
                 u_c(1, j - 1, i,  k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, 1) * mu_c(j + 1, i, 1) *
                 XI33_c(j + 1, i, 1) * (-bof(1, k1)) *
                 u_c(1, j + 1, i,  k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, 1) * mu_c(j + 2, i, 1) *
                 XI33_c(j + 2, i, 1) * (-bof(1, k1)) *
                 u_c(1, j + 2, i,  k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j - 2, i, 1) * mu_c(j - 2, i, 1) *
                 XI13_c(j - 2, i, 1) * (-bof(1, k1)) *
                 u_c(3, j - 2, i,  k1) / 12.0 +
             Jacobian_c(j - 1, i, 1) * mu_c(j - 1, i, 1) *
                 XI13_c(j - 1, i, 1) * (-bof(1, k1)) *
                 u_c(3, j - 1, i,  k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, 1) * mu_c(j + 1, i, 1) *
                 XI13_c(j + 1, i, 1) * (-bof(1, k1)) *
                 u_c(3, j + 1, i,  k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, 1) * mu_c(j + 2, i, 1) *
                 XI13_c(j + 2, i, 1) * (-bof(1, k1)) *
                 u_c(3, j + 2, i,  k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j, i - 2, 1) * mu_c(j, i - 2, 1) *
                 XI33_c(j, i - 2, 1) * (-bof(1, k1)) *
                 u_c(2, j, i - 2,  k1) / 12.0 +
             Jacobian_c(j, i - 1, 1) * mu_c(j, i - 1, 1) *
                 XI33_c(j, i - 1, 1) * (-bof(1, k1)) *
                 u_c(2, j, i - 1,  k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, 1) * mu_c(j, i + 1, 1) *
                 XI33_c(j, i + 1, 1) * (-bof(1, k1)) *
                 u_c(2, j, i + 1,  k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, 1) * mu_c(j, i + 2, 1) *
                 XI33_c(j, i + 2, 1) * (-bof(1, k1)) *
                 u_c(2, j, i + 2,  k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-Jacobian_c(j, i - 2, 1) * mu_c(j, i - 2, 1) *
                 XI23_c(j, i - 2, 1) * (-bof(1, k1)) *
                 u_c(3, j, i - 2,  k1) / 12.0 +
             Jacobian_c(j, i - 1, 1) * mu_c(j, i - 1, 1) *
                 XI23_c(j, i - 1, 1) * (-bof(1, k1)) *
                 u_c(3, j, i - 1,  k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, 1) * mu_c(j, i + 1, 1) *
                 XI23_c(j, i + 1, 1) * (-bof(1, k1)) *
                 u_c(3, j, i + 1,  k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, 1) * mu_c(j, i + 2, 1) *
                 XI23_c(j, i + 2, 1) * (-bof(1, k1)) *
                 u_c(3, j, i + 2,  k1) / 12.0) /
                l2 / h2_c / h3_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             lambda_c(j, i,  k1) * XI33_c(j, i,  k1) *
             (u_c(1, j - 2, i,  k1) / 12.0 -
              u_c(1, j - 1, i,  k1) * 2.0 / 3.0 +
              u_c(1, j + 1, i,  k1) * 2.0 / 3.0 -
              u_c(1, j + 2, i,  k1) / 12.0)) /
                l1 / h3_c / h1_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             mu_c(j, i,  k1) * XI13_c(j, i,  k1) *
             (u_c(3, j - 2, i,  k1) / 12.0 -
              u_c(3, j - 1, i,  k1) * 2.0 / 3.0 +
              u_c(3, j + 1, i,  k1) * 2.0 / 3.0 -
              u_c(3, j + 2, i,  k1) / 12.0)) /
                l1 / h3_c / h1_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             lambda_c(j, i,  k1) * XI33_c(j, i,  k1) *
             (u_c(2, j, i - 2,  k1) / 12.0 -
              u_c(2, j, i - 1,  k1) * 2.0 / 3.0 +
              u_c(2, j, i + 1,  k1) * 2.0 / 3.0 -
              u_c(2, j, i + 2,  k1) / 12.0)) /
                l2 / h2_c / h3_c +
            ( bof(1, k1) * Jacobian_c(j, i,  k1) *
             mu_c(j, i,  k1) * XI23_c(j, i,  k1) *
             (u_c(3, j, i - 2,  k1) / 12.0 -
              u_c(3, j, i - 1,  k1) * 2.0 / 3.0 +
              u_c(3, j, i + 1,  k1) * 2.0 / 3.0 -
              u_c(3, j, i + 2,  k1) / 12.0)) /
                l2 / h2_c / h3_c;
      }
    }
  }

  // project
  // Mass_f1 := Mass_f1 + P(L_c/rho_c) (note lh_c = J*L_c, rho_c has been pre-multiplied by J)
  // first set
  Mass_f1 = 0.0;
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        for (l = -1; l <= 2; l++) {
          Mass_f1(2 * i, 2 * k) = Mass_f1(2 * i, 2 * k) +
                                  P(j) * (P(l) * lh_c(i + l, k + j, 1, 1) /
                                          rho_c(i + l, k + j, 1));
        }
      }
    }
  }
  //
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i - 1, 2 * k) =
            Mass_f1(2 * i - 1, 2 * k) +
            P(j) * lh_c(i, k + j, 1, 1) / rho_c(i, j + k, 1);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i, 2 * k - 1) =
            Mass_f1(2 * i, 2 * k - 1) +
            P(j) * lh_c(i + j, k, 1, 1) / rho_c(i + j, k, 1);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      Mass_f1(2 * i - 1, 2 * k - 1) = lh_c(i, k, 1, 1) / rho_c(i, k, 1);
    }
  }
  // restrict
  // first set 
  //     Vass := Vass -hf*w1*R(rho_f*J*P(L_c/(rho_c)))
  // SIGN CHANGE
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -4; j <= 2; j++) {
        for (l = -4; l <= 2; l++) {
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
	    //              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) -
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) +
              17.0 / 48.0 * h3_f * Rop(j) *
                  (Rop(l) * rho_f(2 * i + l, 2 * k + j, n3_f) *
                   Mass_f1(2 * i + l, 2 * k + j) * 1.0);
        }
      }
    }
  }
  // second set
  Mass_f1 = 0.0;
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        for (l = -1; l <= 2; l++) {
          Mass_f1(2 * i, 2 * k) = Mass_f1(2 * i, 2 * k) +
                                  P(j) * (P(l) * lh_c(i + l, k + j, 1, 2) /
                                          rho_c(i + l, k + j, 1));
        }
      }
    }
  }
  //
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i - 1, 2 * k) =
            Mass_f1(2 * i - 1, 2 * k) +
            P(j) * lh_c(i, k + j, 1, 2) / rho_c(i, j + k, 1);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i, 2 * k - 1) =
            Mass_f1(2 * i, 2 * k - 1) +
            P(j) * lh_c(i + j, k, 1, 2) / rho_c(i + j, k, 1);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      Mass_f1(2 * i - 1, 2 * k - 1) = lh_c(i, k, 1, 2) / rho_c(i, k, 1);
    }
  }
  // restriction
  // second set
  // SIGN CHANGE
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -4; j <= 2; j++) {
        for (l = -4; l <= 2; l++) {
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
	    //              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) -
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) +
              17.0 / 48.0 * h3_f * Rop(j) *
                  (Rop(l) * rho_f(2 * i + l, 2 * k + j, n3_f) *
                   Mass_f1(2 * i + l, 2 * k + j) * 1.0);
        }
      }
    }
  }
  // third set
  Mass_f1 = 0.0;
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        for (l = -1; l <= 2; l++) {
          Mass_f1(2 * i, 2 * k) = Mass_f1(2 * i, 2 * k) +
                                  P(j) * (P(l) * lh_c(i + l, k + j, 1, 3) /
                                          rho_c(i + l, k + j, 1));
        }
      }
    }
  }
  //
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i - 1, 2 * k) =
            Mass_f1(2 * i - 1, 2 * k) +
            P(j) * lh_c(i, k + j, 1, 3) / rho_c(i, j + k, 1);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i, 2 * k - 1) =
            Mass_f1(2 * i, 2 * k - 1) +
            P(j) * lh_c(i + j, k, 1, 3) / rho_c(i + j, k, 1);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      Mass_f1(2 * i - 1, 2 * k - 1) = lh_c(i, k, 1, 3) / rho_c(i, k, 1);
    }
  }
  // restrict
  // third set
  // SIGN CHANGE
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -4; j <= 2; j++) {
        for (l = -4; l <= 2; l++) {
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
	    //              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) -
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) +
              17.0 / 48.0 * h3_f * Rop(j) *
                  (Rop(l) * rho_f(2 * i + l, 2 * k + j, n3_f) *
                   Mass_f1(2 * i + l, 2 * k + j) * 1.0);
        }
      }
    }
  }
  // term 3
  // lh_f := lh_f + J*L_f 
  for (j = -2; j <= n2_f + 3; j++) {
    for (i = -2; i <= n1_f + 3; i++) {
      // second derivative 11  22  12  21
      // first set
      lh_f(i, j, n3_f, 1) =
          lh_f(i, j, n3_f, 1) +
          ((-Jacobian_f(i - 2, j, n3_f) *
                (2.0 * mu_f(i - 2, j, n3_f) + lambda_f(i - 2, j, n3_f)) / 8.0 +
            Jacobian_f(i - 1, j, n3_f) *
                (2.0 * mu_f(i - 1, j, n3_f) + lambda_f(i - 1, j, n3_f)) / 6.0 -
            Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) /
                8.0) *
               u_f(1, i - 2, j, n3_f) +
           (Jacobian_f(i - 2, j, n3_f) *
                (2.0 * mu_f(i - 2, j, n3_f) + lambda_f(i - 2, j, n3_f)) / 6.0 +
            Jacobian_f(i - 1, j, n3_f) *
                (2.0 * mu_f(i - 1, j, n3_f) + lambda_f(i - 1, j, n3_f)) / 2.0 +
            Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) /
                2.0 +
            Jacobian_f(i + 1, j, n3_f) *
                (2.0 * mu_f(i + 1, j, n3_f) + lambda_f(i + 1, j, n3_f)) / 6.0) *
               u_f(1, i - 1, j, n3_f) +
           (-Jacobian_f(i - 2, j, n3_f) *
                (2.0 * mu_f(i - 2, j, n3_f) + lambda_f(i - 2, j, n3_f)) / 24.0 -
            Jacobian_f(i - 1, j, n3_f) *
                (2.0 * mu_f(i - 1, j, n3_f) + lambda_f(i - 1, j, n3_f)) * 5.0 / 6.0 -
            Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                3.0 / 4.0 -
            Jacobian_f(i + 1, j, n3_f) *
                (2.0 * mu_f(i + 1, j, n3_f) + lambda_f(i + 1, j, n3_f)) * 5.0 / 6.0 -
            Jacobian_f(i + 2, j, n3_f) *
                (2.0 * mu_f(i + 2, j, n3_f) + lambda_f(i + 2, j, n3_f)) / 24.0) *
               u_f(1, i - 0, j, n3_f) +
           (Jacobian_f(i - 1, j, n3_f) *
                (2.0 * mu_f(i - 1, j, n3_f) + lambda_f(i - 1, j, n3_f)) / 6.0 +
            Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) /
                2.0 +
            Jacobian_f(i + 1, j, n3_f) *
                (2.0 * mu_f(i + 1, j, n3_f) + lambda_f(i + 1, j, n3_f)) / 2.0 +
            Jacobian_f(i + 2, j, n3_f) *
                (2.0 * mu_f(i + 2, j, n3_f) + lambda_f(i + 2, j, n3_f)) / 6.0) *
               u_f(1, i + 1, j, n3_f) +
           (-Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) /
                8.0 +
            Jacobian_f(i + 1, j, n3_f) *
                (2.0 * mu_f(i + 1, j, n3_f) + lambda_f(i + 1, j, n3_f)) / 6.0 -
            Jacobian_f(i + 2, j, n3_f) *
                (2.0 * mu_f(i + 2, j, n3_f) + lambda_f(i + 2, j, n3_f)) / 8.0) *
               u_f(1, i + 2, j, n3_f)) /
              pow(h1_f, 2) / pow(l1, 2) +
          ((-Jacobian_f(i, j - 2, n3_f) * mu_f(i, j - 2, n3_f) / 8.0 +
            Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) / 6.0 -
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 8.0) *
               u_f(1, i, j - 2, n3_f) +
           (Jacobian_f(i, j - 2, n3_f) * mu_f(i, j - 2, n3_f) / 6.0 +
            Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) / 2.0 +
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 2.0 +
            Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) / 6.0) *
               u_f(1, i, j - 1, n3_f) +
           (-Jacobian_f(i, j - 2, n3_f) * mu_f(i, j - 2, n3_f) / 24.0 -
            Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) * 5.0 / 6.0 -
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) * 3.0 / 4.0 -
            Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) * 5.0 / 6.0 -
            Jacobian_f(i, j + 2, n3_f) * mu_f(i, j + 2, n3_f) / 24.0) *
               u_f(1, i, j - 0, n3_f) +
           (Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) / 6.0 +
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 2.0 +
            Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) / 2.0 +
            Jacobian_f(i, j + 2, n3_f) * mu_f(i, j + 2, n3_f) / 6.0) *
               u_f(1, i, j + 1, n3_f) +
           (-Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 8.0 +
            Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) / 6.0 -
            Jacobian_f(i, j + 2, n3_f) * mu_f(i, j + 2, n3_f) / 8.0) *
               u_f(1, i, j + 2, n3_f)) /
              pow(h2_f, 2) / pow(l2, 2) +
          (Jacobian_f(i - 2, j, n3_f) * lambda_f(i - 2, j, n3_f) *
               (u_f(2, i - 2, j - 2, n3_f) / 12.0 -
                u_f(2, i - 2, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(2, i - 2, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(2, i - 2, j + 2, n3_f) / 12.0) /
               12.0 -
           Jacobian_f(i - 1, j, n3_f) * lambda_f(i - 1, j, n3_f) *
               (u_f(2, i - 1, j - 2, n3_f) / 12.0 -
                u_f(2, i - 1, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(2, i - 1, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(2, i - 1, j + 2, n3_f) / 12.0) *
               2.0 / 3.0 +
           Jacobian_f(i + 1, j, n3_f) * lambda_f(i + 1, j, n3_f) *
               (u_f(2, i + 1, j - 2, n3_f) / 12.0 -
                u_f(2, i + 1, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(2, i + 1, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(2, i + 1, j + 2, n3_f) / 12.0) *
               2.0 / 3.0 -
           Jacobian_f(i + 2, j, n3_f) * lambda_f(i + 2, j, n3_f) *
               (u_f(2, i + 2, j - 2, n3_f) / 12.0 -
                u_f(2, i + 2, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(2, i + 2, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(2, i + 2, j + 2, n3_f) / 12.0) /
               12.0) /
              l1 / l2 / h1_f / h2_f +
          (Jacobian_f(i, j - 2, n3_f) * mu_f(i, j - 2, n3_f) *
               (u_f(2, i - 2, j - 2, n3_f) / 12.0 -
                u_f(2, i - 1, j - 2, n3_f) * 2.0 / 3.0 +
                u_f(2, i + 1, j - 2, n3_f) * 2.0 / 3.0 -
                u_f(2, i + 2, j - 2, n3_f) / 12.0) /
               12.0 -
           Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) *
               (u_f(2, i - 2, j - 1, n3_f) / 12.0 -
                u_f(2, i - 1, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(2, i + 1, j - 1, n3_f) * 2.0 / 3.0 -
                u_f(2, i + 2, j - 1, n3_f) / 12.0) *
               2.0 / 3.0 +
           Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) *
               (u_f(2, i - 2, j + 1, n3_f) / 12.0 -
                u_f(2, i - 1, j + 1, n3_f) * 2.0 / 3.0 +
                u_f(2, i + 1, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(2, i + 2, j + 1, n3_f) / 12.0) *
               2.0 / 3.0 -
           Jacobian_f(i, j + 2, n3_f) * mu_f(i, j + 2, n3_f) *
               (u_f(2, i - 2, j + 2, n3_f) / 12.0 -
                u_f(2, i - 1, j + 2, n3_f) * 2.0 / 3.0 +
                u_f(2, i + 1, j + 2, n3_f) * 2.0 / 3.0 -
                u_f(2, i + 2, j + 2, n3_f) / 12.0) /
               12.0) /
              l1 / l2 / h1_f / h2_f;
      // second set
      lh_f(i, j, n3_f, 2) =
          lh_f(i, j, n3_f, 2) +
          ((-Jacobian_f(i - 2, j, n3_f) * mu_f(i - 2, j, n3_f) / 8.0 +
            Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) / 6.0 -
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 8.0) *
               u_f(2, i - 2, j, n3_f) +
           (Jacobian_f(i - 2, j, n3_f) * mu_f(i - 2, j, n3_f) / 6.0 +
            Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) / 2.0 +
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 2.0 +
            Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) / 6.0) *
               u_f(2, i - 1, j, n3_f) +
           (-Jacobian_f(i - 2, j, n3_f) * mu_f(i - 2, j, n3_f) / 24.0 -
            Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) * 5.0 / 6.0 -
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) * 3.0 / 4.0 -
            Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) * 5.0 / 6.0 -
            Jacobian_f(i + 2, j, n3_f) * mu_f(i + 2, j, n3_f) / 24.0) *
               u_f(2, i - 0, j, n3_f) +
           (Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) / 6.0 +
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 2.0 +
            Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) / 2.0 +
            Jacobian_f(i + 2, j, n3_f) * mu_f(i + 2, j, n3_f) / 6.0) *
               u_f(2, i + 1, j, n3_f) +
           (-Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 8.0 +
            Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) / 6.0 -
            Jacobian_f(i + 2, j, n3_f) * mu_f(i + 2, j, n3_f) / 8.0) *
               u_f(2, i + 2, j, n3_f)) /
              pow(h1_f, 2) / pow(l1, 2) +
          ((-Jacobian_f(i, j - 2, n3_f) *
                (2.0 * mu_f(i, j - 2, n3_f) + lambda_f(i, j - 2, n3_f)) / 8.0 +
            Jacobian_f(i, j - 1, n3_f) *
                (2.0 * mu_f(i, j - 1, n3_f) + lambda_f(i, j - 1, n3_f)) / 6.0 -
            Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) /
                8.0) *
               u_f(2, i, j - 2, n3_f) +
           (Jacobian_f(i, j - 2, n3_f) *
                (2.0 * mu_f(i, j - 2, n3_f) + lambda_f(i, j - 2, n3_f)) / 6.0 +
            Jacobian_f(i, j - 1, n3_f) *
                (2.0 * mu_f(i, j - 1, n3_f) + lambda_f(i, j - 1, n3_f)) / 2.0 +
            Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) /
                2.0 +
            Jacobian_f(i, j + 1, n3_f) *
                (2.0 * mu_f(i, j + 1, n3_f) + lambda_f(i, j + 1, n3_f)) / 6.0) *
               u_f(2, i, j - 1, n3_f) +
           (-Jacobian_f(i, j - 2, n3_f) *
                (2.0 * mu_f(i, j - 2, n3_f) + lambda_f(i, j - 2, n3_f)) / 24.0 -
            Jacobian_f(i, j - 1, n3_f) *
                (2.0 * mu_f(i, j - 1, n3_f) + lambda_f(i, j - 1, n3_f)) * 5.0 / 6.0 -
            Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                3.0 / 4.0 -
            Jacobian_f(i, j + 1, n3_f) *
                (2.0 * mu_f(i, j + 1, n3_f) + lambda_f(i, j + 1, n3_f)) * 5.0 / 6.0 -
            Jacobian_f(i, j + 2, n3_f) *
                (2.0 * mu_f(i, j + 2, n3_f) + lambda_f(i, j + 2, n3_f)) / 24.0) *
               u_f(2, i, j - 0, n3_f) +
           (Jacobian_f(i, j - 1, n3_f) *
                (2.0 * mu_f(i, j - 1, n3_f) + lambda_f(i, j - 1, n3_f)) / 6.0 +
            Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) /
                2.0 +
            Jacobian_f(i, j + 1, n3_f) *
                (2.0 * mu_f(i, j + 1, n3_f) + lambda_f(i, j + 1, n3_f)) / 2.0 +
            Jacobian_f(i, j + 2, n3_f) *
                (2.0 * mu_f(i, j + 2, n3_f) + lambda_f(i, j + 2, n3_f)) / 6.0) *
               u_f(2, i, j + 1, n3_f) +
           (-Jacobian_f(i, j, n3_f) * (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) /
                8.0 +
            Jacobian_f(i, j + 1, n3_f) *
                (2.0 * mu_f(i, j + 1, n3_f) + lambda_f(i, j + 1, n3_f)) / 6.0 -
            Jacobian_f(i, j + 2, n3_f) *
                (2.0 * mu_f(i, j + 2, n3_f) + lambda_f(i, j + 2, n3_f)) / 8.0) *
               u_f(2, i, j + 2, n3_f)) /
              pow(h2_f, 2) / pow(l2, 2) +
          (Jacobian_f(i - 2, j, n3_f) * mu_f(i - 2, j, n3_f) *
               (u_f(1, i - 2, j - 2, n3_f) / 12.0 -
                u_f(1, i - 2, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(1, i - 2, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(1, i - 2, j + 2, n3_f) / 12.0) /
               12.0 -
           Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) *
               (u_f(1, i - 1, j - 2, n3_f) / 12.0 -
                u_f(1, i - 1, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(1, i - 1, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(1, i - 1, j + 2, n3_f) / 12.0) *
               2.0 / 3.0 +
           Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) *
               (u_f(1, i + 1, j - 2, n3_f) / 12.0 -
                u_f(1, i + 1, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(1, i + 1, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(1, i + 1, j + 2, n3_f) / 12.0) *
               2.0 / 3.0 -
           Jacobian_f(i + 2, j, n3_f) * mu_f(i + 2, j, n3_f) *
               (u_f(1, i + 2, j - 2, n3_f) / 12.0 -
                u_f(1, i + 2, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(1, i + 2, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(1, i + 2, j + 2, n3_f) / 12.0) /
               12.0) /
              l1 / l2 / h1_f / h2_f +
          (Jacobian_f(i, j - 2, n3_f) * lambda_f(i, j - 2, n3_f) *
               (u_f(1, i - 2, j - 2, n3_f) / 12.0 -
                u_f(1, i - 1, j - 2, n3_f) * 2.0 / 3.0 +
                u_f(1, i + 1, j - 2, n3_f) * 2.0 / 3.0 -
                u_f(1, i + 2, j - 2, n3_f) / 12.0) /
               12.0 -
           Jacobian_f(i, j - 1, n3_f) * lambda_f(i, j - 1, n3_f) *
               (u_f(1, i - 2, j - 1, n3_f) / 12.0 -
                u_f(1, i - 1, j - 1, n3_f) * 2.0 / 3.0 +
                u_f(1, i + 1, j - 1, n3_f) * 2.0 / 3.0 -
                u_f(1, i + 2, j - 1, n3_f) / 12.0) *
               2.0 / 3.0 +
           Jacobian_f(i, j + 1, n3_f) * lambda_f(i, j + 1, n3_f) *
               (u_f(1, i - 2, j + 1, n3_f) / 12.0 -
                u_f(1, i - 1, j + 1, n3_f) * 2.0 / 3.0 +
                u_f(1, i + 1, j + 1, n3_f) * 2.0 / 3.0 -
                u_f(1, i + 2, j + 1, n3_f) / 12.0) *
               2.0 / 3.0 -
           Jacobian_f(i, j + 2, n3_f) * lambda_f(i, j + 2, n3_f) *
               (u_f(1, i - 2, j + 2, n3_f) / 12.0 -
                u_f(1, i - 1, j + 2, n3_f) * 2.0 / 3.0 +
                u_f(1, i + 1, j + 2, n3_f) * 2.0 / 3.0 -
                u_f(1, i + 2, j + 2, n3_f) / 12.0) /
               12.0) /
              l1 / l2 / h1_f / h2_f;
      // third set
      lh_f(i, j, n3_f, 3) =
          lh_f(i, j, n3_f, 3) +
          ((-Jacobian_f(i - 2, j, n3_f) * mu_f(i - 2, j, n3_f) / 8.0 +
            Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) / 6.0 -
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 8.0) *
               u_f(3, i - 2, j, n3_f) +
           (Jacobian_f(i - 2, j, n3_f) * mu_f(i - 2, j, n3_f) / 6.0 +
            Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) / 2.0 +
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 2.0 +
            Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) / 6.0) *
               u_f(3, i - 1, j, n3_f) +
           (-Jacobian_f(i - 2, j, n3_f) * mu_f(i - 2, j, n3_f) / 24.0 -
            Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) * 5.0 / 6.0 -
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) * 3.0 / 4.0 -
            Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) * 5.0 / 6.0 -
            Jacobian_f(i + 2, j, n3_f) * mu_f(i + 2, j, n3_f) / 24.0) *
               u_f(3, i - 0, j, n3_f) +
           (Jacobian_f(i - 1, j, n3_f) * mu_f(i - 1, j, n3_f) / 6.0 +
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 2.0 +
            Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) / 2.0 +
            Jacobian_f(i + 2, j, n3_f) * mu_f(i + 2, j, n3_f) / 6.0) *
               u_f(3, i + 1, j, n3_f) +
           (-Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 8.0 +
            Jacobian_f(i + 1, j, n3_f) * mu_f(i + 1, j, n3_f) / 6.0 -
            Jacobian_f(i + 2, j, n3_f) * mu_f(i + 2, j, n3_f) / 8.0) *
               u_f(3, i + 2, j, n3_f)) /
              pow(h1_f, 2) / pow(l1, 2) +
          ((-Jacobian_f(i, j - 2, n3_f) * mu_f(i, j - 2, n3_f) / 8.0 +
            Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) / 6.0 -
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 8.0) *
               u_f(3, i, j - 2, n3_f) +
           (Jacobian_f(i, j - 2, n3_f) * mu_f(i, j - 2, n3_f) / 6.0 +
            Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) / 2.0 +
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 2.0 +
            Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) / 6.0) *
               u_f(3, i, j - 1, n3_f) +
           (-Jacobian_f(i, j - 2, n3_f) * mu_f(i, j - 2, n3_f) / 24.0 -
            Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) * 5.0 / 6.0 -
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) * 3.0 / 4.0 -
            Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) * 5.0 / 6.0 -
            Jacobian_f(i, j + 2, n3_f) * mu_f(i, j + 2, n3_f) / 24.0) *
               u_f(3, i, j - 0, n3_f) +
           (Jacobian_f(i, j - 1, n3_f) * mu_f(i, j - 1, n3_f) / 6.0 +
            Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 2.0 +
            Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) / 2.0 +
            Jacobian_f(i, j + 2, n3_f) * mu_f(i, j + 2, n3_f) / 6.0) *
               u_f(3, i, j + 1, n3_f) +
           (-Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / 8.0 +
            Jacobian_f(i, j + 1, n3_f) * mu_f(i, j + 1, n3_f) / 6.0 -
            Jacobian_f(i, j + 2, n3_f) * mu_f(i, j + 2, n3_f) / 8.0) *
               u_f(3, i, j + 2, n3_f)) /
              pow(h2_f, 2) / pow(l2, 2);
    }
  }
  //
  for (j = -2; j <= n2_f + 3; j++) {
    for (k = -2; k <= n1_f + 3; k++) {
      for (k1 = n3_f-7; k1 <= n3_f; k1++) {
        for (m = n3_f-7; m <= n3_f; m++) {
          // second derivative 33
          // first set equation
          lh_f(k, j, n3_f, 1) =
              lh_f(k, j, n3_f, 1) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                    pow(XI13_f(k, j, m), 2) +
                mu_f(k, j, m) *
                    (pow(XI23_f(k, j, m), 2) + pow(XI33_f(k, j, m), 2))) *
               u_f(1, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
               XI23_f(k, j, m) * u_f(2, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
               XI33_f(k, j, m) * u_f(3, k, j, k1)) /
                  pow(h3_f, 2);
          // second set equation
          lh_f(k, j, n3_f, 2) =
              lh_f(k, j, n3_f, 2) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
               XI23_f(k, j, m) * u_f(1, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                    pow(XI23_f(k, j, m), 2) +
                mu_f(k, j, m) *
                    (pow(XI13_f(k, j, m), 2) + pow(XI33_f(k, j, m), 2))) *
               u_f(2, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI23_f(k, j, m) *
               XI33_f(k, j, m) * u_f(3, k, j, k1)) /
                  pow(h3_f, 2);
          // third set equation
          lh_f(k, j, n3_f, 3) =
              lh_f(k, j, n3_f, 3) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
               XI33_f(k, j, m) * u_f(1, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI23_f(k, j, m) *
               XI33_f(k, j, m) * u_f(2, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, n3_f-k1+1, n3_f-m+1) * Jacobian_f(k, j, m) *
               ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                    pow(XI33_f(k, j, m), 2) +
                mu_f(k, j, m) *
                    (pow(XI13_f(k, j, m), 2) + pow(XI23_f(k, j, m), 2))) *
               u_f(3, k, j, k1)) /
                  pow(h3_f, 2);
        }
      }
    }
  }
  //
  for (i = -2; i <= n2_f + 3; i++) {
    for (j = -2; j <= n1_f + 3; j++) {
      for (k1 = 1; k1 <= 6; k1++) {
        // mixed derivative 13  23  31  32
        // first set equation
        lh_f(j, i, n3_f, 1) =
            lh_f(j, i, n3_f, 1) +
            (Jacobian_f(j - 2, i, n3_f) *
                 (2.0 * mu_f(j - 2, i, n3_f) + lambda_f(j - 2, i, n3_f)) *
                 XI13_f(j - 2, i, n3_f) * (-bof(1, k1)) * u_f(1, j - 2, i, n3_f + 1 - k1) /
                 12.0 -
             Jacobian_f(j - 1, i, n3_f) *
                 (2.0 * mu_f(j - 1, i, n3_f) + lambda_f(j - 1, i, n3_f)) *
                 XI13_f(j - 1, i, n3_f) * (-bof(1, k1)) * u_f(1, j - 1, i, n3_f + 1 - k1) * 2.0 /
                 3.0 +
             Jacobian_f(j + 1, i, n3_f) *
                 (2.0 * mu_f(j + 1, i, n3_f) + lambda_f(j + 1, i, n3_f)) *
                 XI13_f(j + 1, i, n3_f) * (-bof(1, k1)) * u_f(1, j + 1, i, n3_f + 1 - k1) * 2.0 /
                 3.0 -
             Jacobian_f(j + 2, i, n3_f) *
                 (2.0 * mu_f(j + 2, i, n3_f) + lambda_f(j + 2, i, n3_f)) *
                 XI13_f(j + 2, i, n3_f) * (-bof(1, k1)) * u_f(1, j + 2, i, n3_f + 1 - k1) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, n3_f) * lambda_f(j - 2, i, n3_f) *
                 XI23_f(j - 2, i, n3_f) * (-bof(1, k1)) * u_f(2, j - 2, i, n3_f + 1 - k1) /
                 12.0 -
             Jacobian_f(j - 1, i, n3_f) * lambda_f(j - 1, i, n3_f) *
                 XI23_f(j - 1, i, n3_f) * (-bof(1, k1)) * u_f(2, j - 1, i, n3_f + 1 - k1) * 2.0 /
                 3.0 +
             Jacobian_f(j + 1, i, n3_f) * lambda_f(j + 1, i, n3_f) *
                 XI23_f(j + 1, i, n3_f) * (-bof(1, k1)) * u_f(2, j + 1, i, n3_f + 1 - k1) * 2.0 /
                 3.0 -
             Jacobian_f(j + 2, i, n3_f) * lambda_f(j + 2, i, n3_f) *
                 XI23_f(j + 2, i, n3_f) * (-bof(1, k1)) * u_f(2, j + 2, i, n3_f + 1 - k1) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, n3_f) * lambda_f(j - 2, i, n3_f) *
                 XI33_f(j - 2, i, n3_f) * (-bof(1, k1)) * u_f(3, j - 2, i, n3_f + 1 - k1) /
                 12.0 -
             Jacobian_f(j - 1, i, n3_f) * lambda_f(j - 1, i, n3_f) *
                 XI33_f(j - 1, i, n3_f) * (-bof(1, k1)) * u_f(3, j - 1, i, n3_f + 1 - k1) * 2.0 /
                 3.0 +
             Jacobian_f(j + 1, i, n3_f) * lambda_f(j + 1, i, n3_f) *
                 XI33_f(j + 1, i, n3_f) * (-bof(1, k1)) * u_f(3, j + 1, i, n3_f + 1 - k1) * 2.0 /
                 3.0 -
             Jacobian_f(j + 2, i, n3_f) * lambda_f(j + 2, i, n3_f) *
                 XI33_f(j + 2, i, n3_f) * (-bof(1, k1)) * u_f(3, j + 2, i, n3_f + 1 - k1) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, n3_f) * mu_f(j, i - 2, n3_f) * XI23_f(j, i - 2, n3_f) *
                 (-bof(1, k1)) * u_f(1, j, i - 2, n3_f + 1 - k1) / 12.0 -
             Jacobian_f(j, i - 1, n3_f) * mu_f(j, i - 1, n3_f) * XI23_f(j, i - 1, n3_f) *
                 (-bof(1, k1)) * u_f(1, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
             Jacobian_f(j, i + 1, n3_f) * mu_f(j, i + 1, n3_f) * XI23_f(j, i + 1, n3_f) *
                 (-bof(1, k1)) * u_f(1, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
             Jacobian_f(j, i + 2, n3_f) * mu_f(j, i + 2, n3_f) * XI23_f(j, i + 2, n3_f) *
                 (-bof(1, k1)) * u_f(1, j, i + 2, n3_f + 1 - k1) / 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, n3_f) * mu_f(j, i - 2, n3_f) * XI13_f(j, i - 2, n3_f) *
                 (-bof(1, k1)) * u_f(2, j, i - 2, n3_f + 1 - k1) / 12.0 -
             Jacobian_f(j, i - 1, n3_f) * mu_f(j, i - 1, n3_f) * XI13_f(j, i - 1, n3_f) *
                 (-bof(1, k1)) * u_f(2, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
             Jacobian_f(j, i + 1, n3_f) * mu_f(j, i + 1, n3_f) * XI13_f(j, i + 1, n3_f) *
                 (-bof(1, k1)) * u_f(2, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
             Jacobian_f(j, i + 2, n3_f) * mu_f(j, i + 2, n3_f) * XI13_f(j, i + 2, n3_f) *
                 (-bof(1, k1)) * u_f(2, j, i + 2, n3_f + 1 - k1) / 12.0) /
                l2 / h2_f / h3_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) *
             (2.0 * mu_f(j, i, n3_f + 1 - k1) + lambda_f(j, i, n3_f + 1 - k1)) * XI13_f(j, i, n3_f + 1 - k1) *
             (u_f(1, j - 2, i, n3_f + 1 - k1) / 12.0 - u_f(1, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(1, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(1, j + 2, i, n3_f + 1 - k1) / 12.0)) /
                l1 / h3_f / h1_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * mu_f(j, i, n3_f + 1 - k1) *
             XI23_f(j, i, n3_f + 1 - k1) *
             (u_f(2, j - 2, i, n3_f + 1 - k1) / 12.0 - u_f(2, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(2, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(2, j + 2, i, n3_f + 1 - k1) / 12.0)) /
                l1 / h3_f / h1_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * mu_f(j, i, n3_f + 1 - k1) *
             XI33_f(j, i, n3_f + 1 - k1) *
             (u_f(3, j - 2, i, n3_f + 1 - k1) / 12.0 - u_f(3, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(3, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(3, j + 2, i, n3_f + 1 - k1) / 12.0)) /
                l1 / h1_f / h3_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * mu_f(j, i, n3_f + 1 - k1) *
             XI23_f(j, i, n3_f + 1 - k1) *
             (u_f(1, j, i - 2, n3_f + 1 - k1) / 12.0 - u_f(1, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(1, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(1, j, i + 2, n3_f + 1 - k1) / 12.0)) /
                l2 / h3_f / h2_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * lambda_f(j, i, n3_f + 1 - k1) *
             XI13_f(j, i, n3_f + 1 - k1) *
             (u_f(2, j, i - 2, n3_f + 1 - k1) / 12.0 - u_f(2, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(2, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(2, j, i + 2, n3_f + 1 - k1) / 12.0)) /
                l2 / h3_f / h2_f;
     // second set equation
        lh_f(j, i, n3_f, 2) =
            lh_f(j, i, n3_f, 2) +
            (Jacobian_f(j - 2, i, n3_f) * mu_f(j - 2, i, n3_f) * XI23_f(j - 2, i, n3_f) *
                 (-bof(1, k1)) * u_f(1, j - 2, i, n3_f + 1 - k1) / 12.0 -
             Jacobian_f(j - 1, i, n3_f) * mu_f(j - 1, i, n3_f) * XI23_f(j - 1, i, n3_f) *
                 (-bof(1, k1)) * u_f(1, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
             Jacobian_f(j + 1, i, n3_f) * mu_f(j + 1, i, n3_f) * XI23_f(j + 1, i, n3_f) *
                 (-bof(1, k1)) * u_f(1, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
             Jacobian_f(j + 2, i, n3_f) * mu_f(j + 2, i, n3_f) * XI23_f(j + 2, i, n3_f) *
                 (-bof(1, k1)) * u_f(1, j + 2, i, n3_f + 1 - k1) / 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, n3_f) * mu_f(j - 2, i, n3_f) * XI13_f(j - 2, i, n3_f) *
                 (-bof(1, k1)) * u_f(2, j - 2, i, n3_f + 1 - k1) / 12.0 -
             Jacobian_f(j - 1, i, n3_f) * mu_f(j - 1, i, n3_f) * XI13_f(j - 1, i, n3_f) *
                 (-bof(1, k1)) * u_f(2, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
             Jacobian_f(j + 1, i, n3_f) * mu_f(j + 1, i, n3_f) * XI13_f(j + 1, i, n3_f) *
                 (-bof(1, k1)) * u_f(2, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
             Jacobian_f(j + 2, i, n3_f) * mu_f(j + 2, i, n3_f) * XI13_f(j + 2, i, n3_f) *
                 (-bof(1, k1)) * u_f(2, j + 2, i, n3_f + 1 - k1) / 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, n3_f) * lambda_f(j, i - 2, n3_f) *
                 XI13_f(j, i - 2, n3_f) * (-bof(1, k1)) * u_f(1, j, i - 2, n3_f + 1 - k1) /
                 12.0 -
             Jacobian_f(j, i - 1, n3_f) * lambda_f(j, i - 1, n3_f) *
                 XI13_f(j, i - 1, n3_f) * (-bof(1, k1)) * u_f(1, j, i - 1, n3_f + 1 - k1) * 2.0 /
                 3.0 +
             Jacobian_f(j, i + 1, n3_f) * lambda_f(j, i + 1, n3_f) *
                 XI13_f(j, i + 1, n3_f) * (-bof(1, k1)) * u_f(1, j, i + 1, n3_f + 1 - k1) * 2.0 /
                 3.0 -
             Jacobian_f(j, i + 2, n3_f) * lambda_f(j, i + 2, n3_f) *
                 XI13_f(j, i + 2, n3_f) * (-bof(1, k1)) * u_f(1, j, i + 2, n3_f + 1 - k1) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, n3_f) *
                 (2.0 * mu_f(j, i - 2, n3_f) + lambda_f(j, i - 2, n3_f)) *
                 XI23_f(j, i - 2, n3_f) * (-bof(1, k1)) * u_f(2, j, i - 2, n3_f + 1 - k1) /
                 12.0 -
             Jacobian_f(j, i - 1, n3_f) *
                 (2.0 * mu_f(j, i - 1, n3_f) + lambda_f(j, i - 1, n3_f)) *
                 XI23_f(j, i - 1, n3_f) * (-bof(1, k1)) * u_f(2, j, i - 1, n3_f + 1 - k1) * 2.0 /
                 3.0 +
             Jacobian_f(j, i + 1, n3_f) *
                 (2.0 * mu_f(j, i + 1, n3_f) + lambda_f(j, i + 1, n3_f)) *
                 XI23_f(j, i + 1, n3_f) * (-bof(1, k1)) * u_f(2, j, i + 1, n3_f + 1 - k1) * 2.0 /
                 3.0 -
             Jacobian_f(j, i + 2, n3_f) *
                 (2.0 * mu_f(j, i + 2, n3_f) + lambda_f(j, i + 2, n3_f)) *
                 XI23_f(j, i + 2, n3_f) * (-bof(1, k1)) * u_f(2, j, i + 2, n3_f + 1 - k1) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, n3_f) * lambda_f(j, i - 2, n3_f) *
                 XI33_f(j, i - 2, n3_f) * (-bof(1, k1)) * u_f(3, j, i - 2, n3_f + 1 - k1) /
                 12.0 -
             Jacobian_f(j, i - 1, n3_f) * lambda_f(j, i - 1, n3_f) *
                 XI33_f(j, i - 1, n3_f) * (-bof(1, k1)) * u_f(3, j, i - 1, n3_f + 1 - k1) * 2.0 /
                 3.0 +
             Jacobian_f(j, i + 1, n3_f) * lambda_f(j, i + 1, n3_f) *
                 XI33_f(j, i + 1, n3_f) * (-bof(1, k1)) * u_f(3, j, i + 1, n3_f + 1 - k1) * 2.0 /
                 3.0 -
             Jacobian_f(j, i + 2, n3_f) * lambda_f(j, i + 2, n3_f) *
                 XI33_f(j, i + 2, n3_f) * (-bof(1, k1)) * u_f(3, j, i + 2, n3_f + 1 - k1) /
                 12.0) /
                l2 / h2_f / h3_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * lambda_f(j, i, n3_f + 1 - k1) *
             XI23_f(j, i, n3_f + 1 - k1) *
             (u_f(1, j - 2, i, n3_f + 1 - k1) / 12.0 - u_f(1, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(1, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(1, j + 2, i, n3_f + 1 - k1) / 12.0)) /
                l1 / h1_f / h3_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * mu_f(j, i, n3_f + 1 - k1) *
             XI13_f(j, i, n3_f + 1 - k1) *
             (u_f(2, j - 2, i, n3_f + 1 - k1) / 12.0 - u_f(2, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(2, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(2, j + 2, i, n3_f + 1 - k1) / 12.0)) /
                l1 / h1_f / h3_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * mu_f(j, i, n3_f + 1 - k1) *
             XI13_f(j, i, n3_f + 1 - k1) *
             (u_f(1, j, i - 2, n3_f + 1 - k1) / 12.0 - u_f(1, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(1, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(1, j, i + 2, n3_f + 1 - k1) / 12.0)) /
                l2 / h3_f / h2_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) *
             (2.0 * mu_f(j, i, n3_f + 1 - k1) + lambda_f(j, i, n3_f + 1 - k1)) * XI23_f(j, i, n3_f + 1 - k1) *
             (u_f(2, j, i - 2, n3_f + 1 - k1) / 12.0 - u_f(2, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(2, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(2, j, i + 2, n3_f + 1 - k1) / 12.0)) /
                l2 / h3_f / h2_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * mu_f(j, i, n3_f + 1 - k1) *
             XI33_f(j, i, n3_f + 1 - k1) *
             (u_f(3, j, i - 2, n3_f + 1 - k1) / 12.0 - u_f(3, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(3, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(3, j, i + 2, n3_f + 1 - k1) / 12.0)) /
                l2 / h3_f / h2_f;
        // third set equation
        lh_f(j, i, n3_f, 3) =
            lh_f(j, i, n3_f, 3) +
            (Jacobian_f(j - 2, i, n3_f) * mu_f(j - 2, i, n3_f) * XI33_f(j - 2, i, n3_f) *
                 (-bof(1, k1)) * u_f(1, j - 2, i, n3_f + 1 - k1) / 12.0 -
             Jacobian_f(j - 1, i, n3_f) * mu_f(j - 1, i, n3_f) * XI33_f(j - 1, i, n3_f) *
                 (-bof(1, k1)) * u_f(1, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
             Jacobian_f(j + 1, i, n3_f) * mu_f(j + 1, i, n3_f) * XI33_f(j + 1, i, n3_f) *
                 (-bof(1, k1)) * u_f(1, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
             Jacobian_f(j + 2, i, n3_f) * mu_f(j + 2, i, n3_f) * XI33_f(j + 2, i, n3_f) *
                 (-bof(1, k1)) * u_f(1, j + 2, i, n3_f + 1 - k1) / 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, n3_f) * mu_f(j - 2, i, n3_f) * XI13_f(j - 2, i, n3_f) *
                 (-bof(1, k1)) * u_f(3, j - 2, i, n3_f + 1 - k1) / 12.0 -
             Jacobian_f(j - 1, i, n3_f) * mu_f(j - 1, i, n3_f) * XI13_f(j - 1, i, n3_f) *
                 (-bof(1, k1)) * u_f(3, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
             Jacobian_f(j + 1, i, n3_f) * mu_f(j + 1, i, n3_f) * XI13_f(j + 1, i, n3_f) *
                 (-bof(1, k1)) * u_f(3, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
             Jacobian_f(j + 2, i, n3_f) * mu_f(j + 2, i, n3_f) * XI13_f(j + 2, i, n3_f) *
                 (-bof(1, k1)) * u_f(3, j + 2, i, n3_f + 1 - k1) / 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, n3_f) * mu_f(j, i - 2, n3_f) * XI33_f(j, i - 2, n3_f) *
                 (-bof(1, k1)) * u_f(2, j, i - 2, n3_f + 1 - k1) / 12.0 -
             Jacobian_f(j, i - 1, n3_f) * mu_f(j, i - 1, n3_f) * XI33_f(j, i - 1, n3_f) *
                 (-bof(1, k1)) * u_f(2, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
             Jacobian_f(j, i + 1, n3_f) * mu_f(j, i + 1, n3_f) * XI33_f(j, i + 1, n3_f) *
                 (-bof(1, k1)) * u_f(2, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
             Jacobian_f(j, i + 2, n3_f) * mu_f(j, i + 2, n3_f) * XI33_f(j, i + 2, n3_f) *
                 (-bof(1, k1)) * u_f(2, j, i + 2, n3_f + 1 - k1) / 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, n3_f) * mu_f(j, i - 2, n3_f) * XI23_f(j, i - 2, n3_f) *
                 (-bof(1, k1)) * u_f(3, j, i - 2, n3_f + 1 - k1) / 12.0 -
             Jacobian_f(j, i - 1, n3_f) * mu_f(j, i - 1, n3_f) * XI23_f(j, i - 1, n3_f) *
                 (-bof(1, k1)) * u_f(3, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
             Jacobian_f(j, i + 1, n3_f) * mu_f(j, i + 1, n3_f) * XI23_f(j, i + 1, n3_f) *
                 (-bof(1, k1)) * u_f(3, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
             Jacobian_f(j, i + 2, n3_f) * mu_f(j, i + 2, n3_f) * XI23_f(j, i + 2, n3_f) *
                 (-bof(1, k1)) * u_f(3, j, i + 2, n3_f + 1 - k1) / 12.0) /
                l2 / h2_f / h3_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * lambda_f(j, i, n3_f + 1 - k1) *
             XI33_f(j, i, n3_f + 1 - k1) *
             (u_f(1, j - 2, i, n3_f + 1 - k1) / 12.0 - u_f(1, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(1, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(1, j + 2, i, n3_f + 1 - k1) / 12.0)) /
                l1 / h1_f / h3_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * mu_f(j, i, n3_f + 1 - k1) *
             XI13_f(j, i, n3_f + 1 - k1) *
             (u_f(3, j - 2, i, n3_f + 1 - k1) / 12.0 - u_f(3, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(3, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(3, j + 2, i, n3_f + 1 - k1) / 12.0)) /
                l1 / h3_f / h1_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * lambda_f(j, i, n3_f + 1 - k1) *
             XI33_f(j, i, n3_f + 1 - k1) *
             (u_f(2, j, i - 2, n3_f + 1 - k1) / 12.0 - u_f(2, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(2, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(2, j, i + 2, n3_f + 1 - k1) / 12.0)) /
                l2 / h2_f / h3_f +
            ((-bof(1, k1)) * Jacobian_f(j, i, n3_f + 1 - k1) * mu_f(j, i, n3_f + 1 - k1) *
             XI23_f(j, i, n3_f + 1 - k1) *
             (u_f(3, j, i - 2, n3_f + 1 - k1) / 12.0 - u_f(3, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
              u_f(3, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 - u_f(3, j, i + 2, n3_f + 1 - k1) / 12.0)) /
                l2 / h2_f / h3_f;
      }
    }
  }

  // scale L_f
  // lh_f := hf*w1*lh_f
  // SIGN CHANGE
  for (j = -2; j <= n2_f + 3; j++) {
    for (i = -2; i <= n1_f + 3; i++) {
      //      lh_f(i, j, n3_f, 1) = lh_f(i, j, n3_f, 1) * 17.0 / 48.0 * h3_f;
      //      lh_f(i, j, n3_f, 2) = lh_f(i, j, n3_f, 2) * 17.0 / 48.0 * h3_f;
      //      lh_f(i, j, n3_f, 3) = lh_f(i, j, n3_f, 3) * 17.0 / 48.0 * h3_f;
      lh_f(i, j, n3_f, 1) = -lh_f(i, j, n3_f, 1) * 17.0 / 48.0 * h3_f;
      lh_f(i, j, n3_f, 2) = -lh_f(i, j, n3_f, 2) * 17.0 / 48.0 * h3_f;
      lh_f(i, j, n3_f, 3) = -lh_f(i, j, n3_f, 3) * 17.0 / 48.0 * h3_f;
    }
  }
  //
  // lh_f := lh_f + B_f(uf)
  for (j = -2; j <= n2_f + 3; j++) {
    for (i = -2; i <= n1_f + 3; i++) {
      for (k = 1; k <= 5; k++) {
        // first set equation
        lh_f(i, j, n3_f, 1) =
            lh_f(i, j, n3_f, 1) +
            Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI13_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) *
                     (pow(XI23_f(i, j, n3_f), 2) + pow(XI33_f(i, j, n3_f), 2))) *
                (-sbop_no_gp(k)) * u_f(1, i, j, n3_f + 1 - k) / h3_f +
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI23_f(i, j, n3_f) * (-sbop_no_gp(k)) *
                u_f(2, i, j, n3_f + 1 - k) / h3_f +
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) * (-sbop_no_gp(k)) *
                u_f(3, i, j, n3_f + 1 - k) / h3_f;
        // second set equation
        lh_f(i, j, n3_f, 2) =
            lh_f(i, j, n3_f, 2) +
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI23_f(i, j, n3_f) * (-sbop_no_gp(k)) *
                u_f(1, i, j, n3_f + 1 - k) / h3_f +
            Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI23_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) *
                     (pow(XI13_f(i, j, n3_f), 2) + pow(XI33_f(i, j, n3_f), 2))) *
                (-sbop_no_gp(k)) * u_f(2, i, j, n3_f + 1 - k) / h3_f +
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI23_f(i, j, n3_f) * XI33_f(i, j, n3_f) * (-sbop_no_gp(k)) *
                u_f(3, i, j, n3_f + 1 - k) / h3_f;
        // third set equation
        lh_f(i, j, n3_f, 3) =
            lh_f(i, j, n3_f, 3) +
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) * (-sbop_no_gp(k)) *
                u_f(1, i, j, n3_f + 1 - k) / h3_f +
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI23_f(i, j, n3_f) * XI33_f(i, j, n3_f) * (-sbop_no_gp(k)) *
                u_f(2, i, j, n3_f + 1 - k) / h3_f +
            Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI33_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) *
                     (pow(XI13_f(i, j, n3_f), 2) + pow(XI23_f(i, j, n3_f), 2))) *
                (-sbop_no_gp(k)) * u_f(3, i, j, n3_f + 1 - k) / h3_f;
      }
    }
  }
  //
  for (k = -2; k <= n2_f + 3; k++) {
    for (i = -2; i <= n1_f + 3; i++) {
      for (j = -2; j <= 2; j++) {
        // 31  32
        // first set equation
        lh_f(i, k, n3_f, 1) =
            lh_f(i, k, n3_f, 1) +
            Jacobian_f(i, k, n3_f) * (2.0 * mu_f(i, k, n3_f) + lambda_f(i, k, n3_f)) /
                l1 * XI13_f(i, k, n3_f) * ux_cof(j) * u_f(1, i + j, k, n3_f) / h1_f +
            Jacobian_f(i, k, n3_f) * mu_f(i, k, n3_f) / l1 * XI23_f(i, k, n3_f) *
                ux_cof(j) * u_f(2, i + j, k, n3_f) / h1_f +
            Jacobian_f(i, k, n3_f) * mu_f(i, k, n3_f) / l1 * XI33_f(i, k, n3_f) *
                ux_cof(j) * u_f(3, i + j, k, n3_f) / h1_f +
            Jacobian_f(i, k, n3_f) * mu_f(i, k, n3_f) / l2 * XI23_f(i, k, n3_f) *
                ux_cof(j) * u_f(1, i, k + j, n3_f) / h2_f +
            Jacobian_f(i, k, n3_f) * lambda_f(i, k, n3_f) / l2 * XI13_f(i, k, n3_f) *
                ux_cof(j) * u_f(2, i, k + j, n3_f) / h2_f;
        // second set equation
        lh_f(i, k, n3_f, 2) =
            lh_f(i, k, n3_f, 2) +
            Jacobian_f(i, k, n3_f) * lambda_f(i, k, n3_f) / l1 * XI23_f(i, k, n3_f) *
                ux_cof(j) * u_f(1, i + j, k, n3_f) / h1_f +
            Jacobian_f(i, k, n3_f) * mu_f(i, k, n3_f) / l1 * XI13_f(i, k, n3_f) *
                ux_cof(j) * u_f(2, i + j, k, n3_f) / h1_f +
            Jacobian_f(i, k, n3_f) * mu_f(i, k, n3_f) / l2 * XI13_f(i, k, n3_f) *
                ux_cof(j) * u_f(1, i, k + j, n3_f) / h2_f +
            Jacobian_f(i, k, n3_f) * (2.0 * mu_f(i, k, n3_f) + lambda_f(i, k, n3_f)) /
                l2 * XI23_f(i, k, n3_f) * ux_cof(j) * u_f(2, i, k + j, n3_f) / h2_f +
            Jacobian_f(i, k, n3_f) * mu_f(i, k, n3_f) / l2 * XI33_f(i, k, n3_f) *
                ux_cof(j) * u_f(3, i, k + j, n3_f) / h2_f;
        // third set equation
        lh_f(i, k, n3_f, 3) =
            lh_f(i, k, n3_f, 3) +
            Jacobian_f(i, k, n3_f) * lambda_f(i, k, n3_f) / l1 * XI33_f(i, k, n3_f) *
                ux_cof(j) * u_f(1, i + j, k, n3_f) / h1_f +
            Jacobian_f(i, k, n3_f) * mu_f(i, k, n3_f) / l1 * XI13_f(i, k, n3_f) *
                ux_cof(j) * u_f(3, i + j, k, n3_f) / h1_f +
            Jacobian_f(i, k, n3_f) * lambda_f(i, k, n3_f) / l2 * XI33_f(i, k, n3_f) *
                ux_cof(j) * u_f(2, i, k + j, n3_f) / h2_f +
            Jacobian_f(i, k, n3_f) * mu_f(i, k, n3_f) / l2 * XI23_f(i, k, n3_f) *
                ux_cof(j) * u_f(3, i, k + j, n3_f) / h2_f;
      }
    }
  }
  // now restrict it to the coarse grid
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -4; j <= 2; j++) {
        for (l = -4; l <= 2; l++) {
          // first set
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) +
              Rop(j) * (Rop(l) * lh_f(2 * i + l, 2 * k + j, n3_f, 1));
          // second set
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) +
              Rop(j) * (Rop(l) * lh_f(2 * i + l, 2 * k + j, n3_f, 2));
          // third set
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) +
              Rop(j) * (Rop(l) * lh_f(2 * i + l, 2 * k + j, n3_f, 3));
        }
      }
    }
  }
} // END INTERFACE_RHS

//-----------------------------------------------------------------------
void CurvilinearInterface::interface_lhs(Farray &LHS, Sarray &Jacobian_c,
                   Sarray &mu_c, Sarray &lambda_c, Sarray &rho_c,
                   Sarray &rho_f, Sarray &XI13_c, Sarray &XI23_c,
                   Sarray &XI33_c, Farray &P, Farray &Sb, Farray &Rop,
                   Sarray &u_c, Farray &ghcof, PackArgs &a ) {
  float_sw4 int_cof;

  int i, j, k, i1, j1, k1, l;

  auto dim = a.dim;
  auto h3_c = a.h3_c;
  auto h3_f = a.h3_f;

  auto nrg = a.nrg;
  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;
  auto n3_f = a.n3_f;

  Farray u_temp(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1, dim);
  u_temp = 0.0;
  //  u_temp.owns_data=false;
  for (j = 1; j <= n2_c; j++) {
    for (i = 1; i <= n1_c; i++) {
      u_temp(i, j, 1) = u_c(1, i, j, 0);
      u_temp(i, j, 2) = u_c(2, i, j, 0);
      u_temp(i, j, 3) = u_c(3, i, j, 0);
    }
  }

  int_cof = 17.0 / 48.0 * h3_f * ghcof(1) / pow(h3_c, 2);
  // SIGN CHANGE
  int_cof = -int_cof;

  LHS = 0.0;
  for (l = 1; l <= n2_c; l++) {
    for (k = 1; k <= n1_c; k++) {
      for (j = -4; j <= 2; j += 2) {
        for (i = -4; i <= 2; i += 2) {
          for (j1 = -1; j1 <= 2; j1++) {
            for (i1 = -1; i1 <= 2; i1++) {
              // first set equation
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
                  LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      ((2.0 * mu_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                        lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                           pow(XI13_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                               2) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                           (pow(XI23_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                                2) +
                            pow(XI33_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                                2))) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 1) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                      XI13_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      XI23_c(k + i / 2 + i1, l + j / 2 + j1, 1) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 2) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                      XI13_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      XI33_c(k + i / 2 + i1, l + j / 2 + j1, 1) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 3) * int_cof;
              // second set equation
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
                  LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                      XI13_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      XI23_c(k + i / 2 + i1, l + j / 2 + j1, 1) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 1) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      ((2.0 * mu_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                        lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                           pow(XI23_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                               2) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                           (pow(XI13_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                                2) +
                            pow(XI33_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                                2))) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 2) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                      XI23_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      XI33_c(k + i / 2 + i1, l + j / 2 + j1, 1) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 3) * int_cof;
              // third set equation
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
                  LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                      XI13_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      XI33_c(k + i / 2 + i1, l + j / 2 + j1, 1) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 1) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                      XI23_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      XI33_c(k + i / 2 + i1, l + j / 2 + j1, 1) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 2) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, n3_f) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      ((2.0 * mu_c(k + i / 2 + i1, l + j / 2 + j1, 1) +
                        lambda_c(k + i / 2 + i1, l + j / 2 + j1, 1)) *
                           pow(XI33_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                               2) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                           (pow(XI13_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                                2) +
                            pow(XI23_c(k + i / 2 + i1, l + j / 2 + j1, 1),
                                2))) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, 1) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 3) * int_cof;
            }
          }
        }
      }
      //
      for (j = -4; j <= 2; j += 2) {
        for (j1 = -1; j1 <= 2; j1++) {
          // first set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  ((2.0 * mu_c(k, l + j / 2 + j1, 1) +
                    lambda_c(k, l + j / 2 + j1, 1)) *
                       pow(XI13_c(k, l + j / 2 + j1, 1), 2) +
                   mu_c(k, l + j / 2 + j1, 1) *
                       (pow(XI23_c(k, l + j / 2 + j1, 1), 2) +
                        pow(XI33_c(k, l + j / 2 + j1, 1), 2))) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 1) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  (lambda_c(k, l + j / 2 + j1, 1) +
                   mu_c(k, l + j / 2 + j1, 1)) *
                  XI13_c(k, l + j / 2 + j1, 1) *
                  XI23_c(k, l + j / 2 + j1, 1) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 2) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  (lambda_c(k, l + j / 2 + j1, 1) +
                   mu_c(k, l + j / 2 + j1, 1)) *
                  XI13_c(k, l + j / 2 + j1, 1) *
                  XI33_c(k, l + j / 2 + j1, 1) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 3) * int_cof;
          // second set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  (lambda_c(k, l + j / 2 + j1, 1) +
                   mu_c(k, l + j / 2 + j1, 1)) *
                  XI13_c(k, l + j / 2 + j1, 1) *
                  XI23_c(k, l + j / 2 + j1, 1) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 1) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  ((2.0 * mu_c(k, l + j / 2 + j1, 1) +
                    lambda_c(k, l + j / 2 + j1, 1)) *
                       pow(XI23_c(k, l + j / 2 + j1, 1), 2) +
                   mu_c(k, l + j / 2 + j1, 1) *
                       (pow(XI13_c(k, l + j / 2 + j1, 1), 2) +
                        pow(XI33_c(k, l + j / 2 + j1, 1), 2))) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 2) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  (lambda_c(k, l + j / 2 + j1, 1) +
                   mu_c(k, l + j / 2 + j1, 1)) *
                  XI23_c(k, l + j / 2 + j1, 1) *
                  XI33_c(k, l + j / 2 + j1, 1) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 3) * int_cof;
          // third set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  (lambda_c(k, l + j / 2 + j1, 1) +
                   mu_c(k, l + j / 2 + j1, 1)) *
                  XI13_c(k, l + j / 2 + j1, 1) *
                  XI33_c(k, l + j / 2 + j1, 1) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 1) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  (lambda_c(k, l + j / 2 + j1, 1) +
                   mu_c(k, l + j / 2 + j1, 1)) *
                  XI23_c(k, l + j / 2 + j1, 1) *
                  XI33_c(k, l + j / 2 + j1, 1) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 2) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, n3_f) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, 1) *
                  ((2.0 * mu_c(k, l + j / 2 + j1, 1) +
                    lambda_c(k, l + j / 2 + j1, 1)) *
                       pow(XI33_c(k, l + j / 2 + j1, 1), 2) +
                   mu_c(k, l + j / 2 + j1, 1) *
                       (pow(XI13_c(k, l + j / 2 + j1, 1), 2) +
                        pow(XI23_c(k, l + j / 2 + j1, 1), 2))) /
                  rho_c(k, l + j / 2 + j1, 1) *
                  u_temp(k, l + j / 2 + j1, 3) * int_cof;
        }
      }
      //
      for (i = -4; i <= 2; i += 2) {
        for (i1 = -1; i1 <= 2; i1++) {
          // first set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  ((2.0 * mu_c(k + i / 2 + i1, l, 1) +
                    lambda_c(k + i / 2 + i1, l, 1)) *
                       pow(XI13_c(k + i / 2 + i1, l, 1), 2) +
                   mu_c(k + i / 2 + i1, l, 1) *
                       (pow(XI23_c(k + i / 2 + i1, l, 1), 2) +
                        pow(XI33_c(k + i / 2 + i1, l, 1), 2))) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 1) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  (lambda_c(k + i / 2 + i1, l, 1) +
                   mu_c(k + i / 2 + i1, l, 1)) *
                  XI13_c(k + i / 2 + i1, l, 1) *
                  XI23_c(k + i / 2 + i1, l, 1) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 2) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  (lambda_c(k + i / 2 + i1, l, 1) +
                   mu_c(k + i / 2 + i1, l, 1)) *
                  XI13_c(k + i / 2 + i1, l, 1) *
                  XI33_c(k + i / 2 + i1, l, 1) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 3) * int_cof;
          // second set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  (lambda_c(k + i / 2 + i1, l, 1) +
                   mu_c(k + i / 2 + i1, l, 1)) *
                  XI13_c(k + i / 2 + i1, l, 1) *
                  XI23_c(k + i / 2 + i1, l, 1) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 1) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  ((2.0 * mu_c(k + i / 2 + i1, l, 1) +
                    lambda_c(k + i / 2 + i1, l, 1)) *
                       pow(XI23_c(k + i / 2 + i1, l, 1), 2) +
                   mu_c(k + i / 2 + i1, l, 1) *
                       (pow(XI13_c(k + i / 2 + i1, l, 1), 2) +
                        pow(XI33_c(k + i / 2 + i1, l, 1), 2))) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 2) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  (lambda_c(k + i / 2 + i1, l, 1) +
                   mu_c(k + i / 2 + i1, l, 1)) *
                  XI23_c(k + i / 2 + i1, l, 1) *
                  XI33_c(k + i / 2 + i1, l, 1) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 3) * int_cof;
          // third set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  (lambda_c(k + i / 2 + i1, l, 1) +
                   mu_c(k + i / 2 + i1, l, 1)) *
                  XI13_c(k + i / 2 + i1, l, 1) *
                  XI33_c(k + i / 2 + i1, l, 1) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 1) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  (lambda_c(k + i / 2 + i1, l, 1) +
                   mu_c(k + i / 2 + i1, l, 1)) *
                  XI23_c(k + i / 2 + i1, l, 1) *
                  XI33_c(k + i / 2 + i1, l, 1) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 2) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, n3_f) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, 1) *
                  ((2.0 * mu_c(k + i / 2 + i1, l, 1) +
                    lambda_c(k + i / 2 + i1, l, 1)) *
                       pow(XI33_c(k + i / 2 + i1, l, 1), 2) +
                   mu_c(k + i / 2 + i1, l, 1) *
                       (pow(XI13_c(k + i / 2 + i1, l, 1), 2) +
                        pow(XI23_c(k + i / 2 + i1, l, 1), 2))) /
                  rho_c(k + i / 2 + i1, l, 1) *
                  u_temp(k + i / 2 + i1, l, 3) * int_cof;
        }
      }
      //
      // first set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI13_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI23_c(k, l, 1), 2) + pow(XI33_c(k, l, 1), 2))) /
              rho_c(k, l, 1) * u_temp(k, l, 1) * int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI23_c(k, l, 1) / rho_c(k, l, 1) * u_temp(k, l, 2) *
              int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI33_c(k, l, 1) / rho_c(k, l, 1) * u_temp(k, l, 3) *
              int_cof;
      // second set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI23_c(k, l, 1) / rho_c(k, l, 1) * u_temp(k, l, 1) *
              int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI23_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI13_c(k, l, 1), 2) + pow(XI33_c(k, l, 1), 2))) /
              rho_c(k, l, 1) * u_temp(k, l, 2) * int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
              XI33_c(k, l, 1) / rho_c(k, l, 1) * u_temp(k, l, 3) *
              int_cof;
      // third set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI33_c(k, l, 1) / rho_c(k, l, 1) * u_temp(k, l, 1) *
              int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
              XI33_c(k, l, 1) / rho_c(k, l, 1) * u_temp(k, l, 2) *
              int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, n3_f) *
              Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI33_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI13_c(k, l, 1), 2) + pow(XI23_c(k, l, 1), 2))) /
              rho_c(k, l, 1) * u_temp(k, l, 3) * int_cof;
      //
      // first set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) -
          Sb(0) * Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI13_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI23_c(k, l, 1), 2) + pow(XI33_c(k, l, 1), 2))) /
              h3_c * u_temp(k, l, 1) -
          Sb(0) * Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI23_c(k, l, 1) / h3_c * u_temp(k, l, 2) -
          Sb(0) * Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI33_c(k, l, 1) / h3_c * u_temp(k, l, 3);
      // second set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) -
          Sb(0) * Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI23_c(k, l, 1) / h3_c * u_temp(k, l, 1) -
          Sb(0) * Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI23_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI13_c(k, l, 1), 2) + pow(XI33_c(k, l, 1), 2))) /
              h3_c * u_temp(k, l, 2) -
          Sb(0) * Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
              XI33_c(k, l, 1) / h3_c * u_temp(k, l, 3);
      // third set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) -
          Sb(0) * Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI13_c(k, l, 1) *
              XI33_c(k, l, 1) / h3_c * u_temp(k, l, 1) -
          Sb(0) * Jacobian_c(k, l, 1) *
              (lambda_c(k, l, 1) + mu_c(k, l, 1)) * XI23_c(k, l, 1) *
              XI33_c(k, l, 1) / h3_c * u_temp(k, l, 2) -
          Sb(0) * Jacobian_c(k, l, 1) *
              ((2.0 * mu_c(k, l, 1) + lambda_c(k, l, 1)) *
                   pow(XI33_c(k, l, 1), 2) +
               mu_c(k, l, 1) *
                   (pow(XI13_c(k, l, 1), 2) + pow(XI23_c(k, l, 1), 2))) /
              h3_c * u_temp(k, l, 3);
    }
  }
}  // END INTERFACE_LHS

//-----------------------------------------------------------------------
void CurvilinearInterface::define_coeffs()
{
  P(-1) = -1.0 / 16.0;
  P(0) = 9.0 / 16.0;
  P(1) = 9.0 / 16.0;
  P(2) = -1.0 / 16.0;

  Rop(-4) = -1.0 / 32.0;
  Rop(-3) = 0.0;
  Rop(-2) = 9.0 / 32.0;
  Rop(-1) = 1.0 / 2.0;
  Rop(0) = 9.0 / 32.0;
  Rop(1) = 0.0;
  Rop(2) = -1.0 / 32.0;

  ux_cof(-2) = 1.0 / 12.0;
  ux_cof(-1) = -2.0 / 3.0;
  ux_cof(0) = 0.0;
  ux_cof(1) = 2.0 / 3.0;
  ux_cof(2) = -1.0 / 12.0;

  varcoeffs4( acof, ghcof, Sb);
  varcoeff_noghost();
  dx_46(bof);  
}

//-----------------------------------------------------------------------
void CurvilinearInterface::varcoeffs4( Farray& a_acof, Farray& a_ghcof, Farray& a_Sb) {
  // acofs(i,j,k) is coefficient of a(k) in stencil coefficient (i,j);
  // ghcof is coefficient of ghost point, a(1)*ghcof*u(0) in stencil at i=1.;
  a_ghcof(1) = 12.0 / 17.0;
  a_ghcof(2) = 0;
  a_ghcof(3) = 0;
  a_ghcof(4) = 0;
  a_ghcof(5) = 0;
  a_ghcof(6) = 0;
  a_acof(1, 1, 1) = 104.0 / 289.0;
  a_acof(1, 1, 2) = -2476335.0 / 2435692.0;
  a_acof(1, 1, 3) = -16189.0 / 84966.0;
  a_acof(1, 1, 4) = -9.0 / 3332.0;
  a_acof(1, 1, 5) = 0;
  a_acof(1, 1, 6) = 0;
  a_acof(1, 1, 7) = 0;
  a_acof(1, 1, 8) = 0;
  a_acof(1, 2, 1) = -516.0 / 289.0;
  a_acof(1, 2, 2) = 544521.0 / 1217846.0;
  a_acof(1, 2, 3) = 2509879.0 / 3653538.0;
  a_acof(1, 2, 4) = 0;
  a_acof(1, 2, 5) = 0;
  a_acof(1, 2, 6) = 0;
  a_acof(1, 2, 7) = 0;
  a_acof(1, 2, 8) = 0;
  a_acof(1, 3, 1) = 312.0 / 289.0;
  a_acof(1, 3, 2) = 1024279.0 / 2435692.0;
  a_acof(1, 3, 3) = -687797.0 / 1217846.0;
  a_acof(1, 3, 4) = 177.0 / 3332.0;
  a_acof(1, 3, 5) = 0;
  a_acof(1, 3, 6) = 0;
  a_acof(1, 3, 7) = 0;
  a_acof(1, 3, 8) = 0;
  a_acof(1, 4, 1) = -104.0 / 289.0;
  a_acof(1, 4, 2) = 181507.0 / 1217846.0;
  a_acof(1, 4, 3) = 241309.0 / 3653538.0;
  a_acof(1, 4, 4) = 0;
  a_acof(1, 4, 5) = 0;
  a_acof(1, 4, 6) = 0;
  a_acof(1, 4, 7) = 0;
  a_acof(1, 4, 8) = 0;
  a_acof(1, 5, 1) = 0;
  a_acof(1, 5, 2) = 0;
  a_acof(1, 5, 3) = 5.0 / 2193.0;
  a_acof(1, 5, 4) = -48.0 / 833.0;
  a_acof(1, 5, 5) = 0;
  a_acof(1, 5, 6) = 0;
  a_acof(1, 5, 7) = 0;
  a_acof(1, 5, 8) = 0;
  a_acof(1, 6, 1) = 0;
  a_acof(1, 6, 2) = 0;
  a_acof(1, 6, 3) = 0;
  a_acof(1, 6, 4) = 6.0 / 833.0;
  a_acof(1, 6, 5) = 0;
  a_acof(1, 6, 6) = 0;
  a_acof(1, 6, 7) = 0;
  a_acof(1, 6, 8) = 0;
  a_acof(1, 7, 1) = 0;
  a_acof(1, 7, 2) = 0;
  a_acof(1, 7, 3) = 0;
  a_acof(1, 7, 4) = 0;
  a_acof(1, 7, 5) = 0;
  a_acof(1, 7, 6) = 0;
  a_acof(1, 7, 7) = 0;
  a_acof(1, 7, 8) = 0;
  a_acof(1, 8, 1) = 0;
  a_acof(1, 8, 2) = 0;
  a_acof(1, 8, 3) = 0;
  a_acof(1, 8, 4) = 0;
  a_acof(1, 8, 5) = 0;
  a_acof(1, 8, 6) = 0;
  a_acof(1, 8, 7) = 0;
  a_acof(1, 8, 8) = 0;
  a_acof(2, 1, 1) = 12.0 / 17.0;
  a_acof(2, 1, 2) = 544521.0 / 4226642.0;
  a_acof(2, 1, 3) = 2509879.0 / 12679926.0;
  a_acof(2, 1, 4) = 0;
  a_acof(2, 1, 5) = 0;
  a_acof(2, 1, 6) = 0;
  a_acof(2, 1, 7) = 0;
  a_acof(2, 1, 8) = 0;
  a_acof(2, 2, 1) = -59.0 / 68.0;
  a_acof(2, 2, 2) = -1633563.0 / 4226642.0;
  a_acof(2, 2, 3) = -21510077.0 / 25359852.0;
  a_acof(2, 2, 4) = -12655.0 / 372939.0;
  a_acof(2, 2, 5) = 0;
  a_acof(2, 2, 6) = 0;
  a_acof(2, 2, 7) = 0;
  a_acof(2, 2, 8) = 0;
  a_acof(2, 3, 1) = 2.0 / 17.0;
  a_acof(2, 3, 2) = 1633563.0 / 4226642.0;
  a_acof(2, 3, 3) = 2565299.0 / 4226642.0;
  a_acof(2, 3, 4) = 40072.0 / 372939.0;
  a_acof(2, 3, 5) = 0;
  a_acof(2, 3, 6) = 0;
  a_acof(2, 3, 7) = 0;
  a_acof(2, 3, 8) = 0;
  a_acof(2, 4, 1) = 3.0 / 68.0;
  a_acof(2, 4, 2) = -544521.0 / 4226642.0;
  a_acof(2, 4, 3) = 987685.0 / 25359852.0;
  a_acof(2, 4, 4) = -14762.0 / 124313.0;
  a_acof(2, 4, 5) = 0;
  a_acof(2, 4, 6) = 0;
  a_acof(2, 4, 7) = 0;
  a_acof(2, 4, 8) = 0;
  a_acof(2, 5, 1) = 0;
  a_acof(2, 5, 2) = 0;
  a_acof(2, 5, 3) = 1630.0 / 372939.0;
  a_acof(2, 5, 4) = 18976.0 / 372939.0;
  a_acof(2, 5, 5) = 0;
  a_acof(2, 5, 6) = 0;
  a_acof(2, 5, 7) = 0;
  a_acof(2, 5, 8) = 0;
  a_acof(2, 6, 1) = 0;
  a_acof(2, 6, 2) = 0;
  a_acof(2, 6, 3) = 0;
  a_acof(2, 6, 4) = -1.0 / 177.0;
  a_acof(2, 6, 5) = 0;
  a_acof(2, 6, 6) = 0;
  a_acof(2, 6, 7) = 0;
  a_acof(2, 6, 8) = 0;
  a_acof(2, 7, 1) = 0;
  a_acof(2, 7, 2) = 0;
  a_acof(2, 7, 3) = 0;
  a_acof(2, 7, 4) = 0;
  a_acof(2, 7, 5) = 0;
  a_acof(2, 7, 6) = 0;
  a_acof(2, 7, 7) = 0;
  a_acof(2, 7, 8) = 0;
  a_acof(2, 8, 1) = 0;
  a_acof(2, 8, 2) = 0;
  a_acof(2, 8, 3) = 0;
  a_acof(2, 8, 4) = 0;
  a_acof(2, 8, 5) = 0;
  a_acof(2, 8, 6) = 0;
  a_acof(2, 8, 7) = 0;
  a_acof(2, 8, 8) = 0;
  a_acof(3, 1, 1) = -96.0 / 731.0;
  a_acof(3, 1, 2) = 1024279.0 / 6160868.0;
  a_acof(3, 1, 3) = -687797.0 / 3080434.0;
  a_acof(3, 1, 4) = 177.0 / 8428.0;
  a_acof(3, 1, 5) = 0;
  a_acof(3, 1, 6) = 0;
  a_acof(3, 1, 7) = 0;
  a_acof(3, 1, 8) = 0;
  a_acof(3, 2, 1) = 118.0 / 731.0;
  a_acof(3, 2, 2) = 1633563.0 / 3080434.0;
  a_acof(3, 2, 3) = 2565299.0 / 3080434.0;
  a_acof(3, 2, 4) = 40072.0 / 271803.0;
  a_acof(3, 2, 5) = 0;
  a_acof(3, 2, 6) = 0;
  a_acof(3, 2, 7) = 0;
  a_acof(3, 2, 8) = 0;
  a_acof(3, 3, 1) = -16.0 / 731.0;
  a_acof(3, 3, 2) = -5380447.0 / 6160868.0;
  a_acof(3, 3, 3) = -3569115.0 / 3080434.0;
  a_acof(3, 3, 4) = -331815.0 / 362404.0;
  a_acof(3, 3, 5) = -283.0 / 6321.0;
  a_acof(3, 3, 6) = 0;
  a_acof(3, 3, 7) = 0;
  a_acof(3, 3, 8) = 0;
  a_acof(3, 4, 1) = -6.0 / 731.0;
  a_acof(3, 4, 2) = 544521.0 / 3080434.0;
  a_acof(3, 4, 3) = 2193521.0 / 3080434.0;
  a_acof(3, 4, 4) = 8065.0 / 12943.0;
  a_acof(3, 4, 5) = 381.0 / 2107.0;
  a_acof(3, 4, 6) = 0;
  a_acof(3, 4, 7) = 0;
  a_acof(3, 4, 8) = 0;
  a_acof(3, 5, 1) = 0;
  a_acof(3, 5, 2) = 0;
  a_acof(3, 5, 3) = -14762.0 / 90601.0;
  a_acof(3, 5, 4) = 32555.0 / 271803.0;
  a_acof(3, 5, 5) = -283.0 / 2107.0;
  a_acof(3, 5, 6) = 0;
  a_acof(3, 5, 7) = 0;
  a_acof(3, 5, 8) = 0;
  a_acof(3, 6, 1) = 0;
  a_acof(3, 6, 2) = 0;
  a_acof(3, 6, 3) = 0;
  a_acof(3, 6, 4) = 9.0 / 2107.0;
  a_acof(3, 6, 5) = -11.0 / 6321.0;
  a_acof(3, 6, 6) = 0;
  a_acof(3, 6, 7) = 0;
  a_acof(3, 6, 8) = 0;
  a_acof(3, 7, 1) = 0;
  a_acof(3, 7, 2) = 0;
  a_acof(3, 7, 3) = 0;
  a_acof(3, 7, 4) = 0;
  a_acof(3, 7, 5) = 0;
  a_acof(3, 7, 6) = 0;
  a_acof(3, 7, 7) = 0;
  a_acof(3, 7, 8) = 0;
  a_acof(3, 8, 1) = 0;
  a_acof(3, 8, 2) = 0;
  a_acof(3, 8, 3) = 0;
  a_acof(3, 8, 4) = 0;
  a_acof(3, 8, 5) = 0;
  a_acof(3, 8, 6) = 0;
  a_acof(3, 8, 7) = 0;
  a_acof(3, 8, 8) = 0;
  a_acof(4, 1, 1) = -36.0 / 833.0;
  a_acof(4, 1, 2) = 181507.0 / 3510262.0;
  a_acof(4, 1, 3) = 241309.0 / 10530786.0;
  a_acof(4, 1, 4) = 0;
  a_acof(4, 1, 5) = 0;
  a_acof(4, 1, 6) = 0;
  a_acof(4, 1, 7) = 0;
  a_acof(4, 1, 8) = 0;
  a_acof(4, 2, 1) = 177.0 / 3332.0;
  a_acof(4, 2, 2) = -544521.0 / 3510262.0;
  a_acof(4, 2, 3) = 987685.0 / 21061572.0;
  a_acof(4, 2, 4) = -14762.0 / 103243.0;
  a_acof(4, 2, 5) = 0;
  a_acof(4, 2, 6) = 0;
  a_acof(4, 2, 7) = 0;
  a_acof(4, 2, 8) = 0;
  a_acof(4, 3, 1) = -6.0 / 833.0;
  a_acof(4, 3, 2) = 544521.0 / 3510262.0;
  a_acof(4, 3, 3) = 2193521.0 / 3510262.0;
  a_acof(4, 3, 4) = 8065.0 / 14749.0;
  a_acof(4, 3, 5) = 381.0 / 2401.0;
  a_acof(4, 3, 6) = 0;
  a_acof(4, 3, 7) = 0;
  a_acof(4, 3, 8) = 0;
  a_acof(4, 4, 1) = -9.0 / 3332.0;
  a_acof(4, 4, 2) = -181507.0 / 3510262.0;
  a_acof(4, 4, 3) = -2647979.0 / 3008796.0;
  a_acof(4, 4, 4) = -80793.0 / 103243.0;
  a_acof(4, 4, 5) = -1927.0 / 2401.0;
  a_acof(4, 4, 6) = -2.0 / 49.0;
  a_acof(4, 4, 7) = 0;
  a_acof(4, 4, 8) = 0;
  a_acof(4, 5, 1) = 0;
  a_acof(4, 5, 2) = 0;
  a_acof(4, 5, 3) = 57418.0 / 309729.0;
  a_acof(4, 5, 4) = 51269.0 / 103243.0;
  a_acof(4, 5, 5) = 1143.0 / 2401.0;
  a_acof(4, 5, 6) = 8.0 / 49.0;
  a_acof(4, 5, 7) = 0;
  a_acof(4, 5, 8) = 0;
  a_acof(4, 6, 1) = 0;
  a_acof(4, 6, 2) = 0;
  a_acof(4, 6, 3) = 0;
  a_acof(4, 6, 4) = -283.0 / 2401.0;
  a_acof(4, 6, 5) = 403.0 / 2401.0;
  a_acof(4, 6, 6) = -6.0 / 49.0;
  a_acof(4, 6, 7) = 0;
  a_acof(4, 6, 8) = 0;
  a_acof(4, 7, 1) = 0;
  a_acof(4, 7, 2) = 0;
  a_acof(4, 7, 3) = 0;
  a_acof(4, 7, 4) = 0;
  a_acof(4, 7, 5) = 0;
  a_acof(4, 7, 6) = 0;
  a_acof(4, 7, 7) = 0;
  a_acof(4, 7, 8) = 0;
  a_acof(4, 8, 1) = 0;
  a_acof(4, 8, 2) = 0;
  a_acof(4, 8, 3) = 0;
  a_acof(4, 8, 4) = 0;
  a_acof(4, 8, 5) = 0;
  a_acof(4, 8, 6) = 0;
  a_acof(4, 8, 7) = 0;
  a_acof(4, 8, 8) = 0;
  a_acof(5, 1, 1) = 0;
  a_acof(5, 1, 2) = 0;
  a_acof(5, 1, 3) = 5.0 / 6192.0;
  a_acof(5, 1, 4) = -1.0 / 49.0;
  a_acof(5, 1, 5) = 0;
  a_acof(5, 1, 6) = 0;
  a_acof(5, 1, 7) = 0;
  a_acof(5, 1, 8) = 0;
  a_acof(5, 2, 1) = 0;
  a_acof(5, 2, 2) = 0;
  a_acof(5, 2, 3) = 815.0 / 151704.0;
  a_acof(5, 2, 4) = 1186.0 / 18963.0;
  a_acof(5, 2, 5) = 0;
  a_acof(5, 2, 6) = 0;
  a_acof(5, 2, 7) = 0;
  a_acof(5, 2, 8) = 0;
  a_acof(5, 3, 1) = 0;
  a_acof(5, 3, 2) = 0;
  a_acof(5, 3, 3) = -7381.0 / 50568.0;
  a_acof(5, 3, 4) = 32555.0 / 303408.0;
  a_acof(5, 3, 5) = -283.0 / 2352.0;
  a_acof(5, 3, 6) = 0;
  a_acof(5, 3, 7) = 0;
  a_acof(5, 3, 8) = 0;
  a_acof(5, 4, 1) = 0;
  a_acof(5, 4, 2) = 0;
  a_acof(5, 4, 3) = 28709.0 / 151704.0;
  a_acof(5, 4, 4) = 51269.0 / 101136.0;
  a_acof(5, 4, 5) = 381.0 / 784.0;
  a_acof(5, 4, 6) = 1.0 / 6.0;
  a_acof(5, 4, 7) = 0;
  a_acof(5, 4, 8) = 0;
  a_acof(5, 5, 1) = 0;
  a_acof(5, 5, 2) = 0;
  a_acof(5, 5, 3) = -349.0 / 7056.0;
  a_acof(5, 5, 4) = -247951.0 / 303408.0;
  a_acof(5, 5, 5) = -577.0 / 784.0;
  a_acof(5, 5, 6) = -5.0 / 6.0;
  a_acof(5, 5, 7) = -1.0 / 24.0;
  a_acof(5, 5, 8) = 0;
  a_acof(5, 6, 1) = 0;
  a_acof(5, 6, 2) = 0;
  a_acof(5, 6, 3) = 0;
  a_acof(5, 6, 4) = 1135.0 / 7056.0;
  a_acof(5, 6, 5) = 1165.0 / 2352.0;
  a_acof(5, 6, 6) = 1.0 / 2.0;
  a_acof(5, 6, 7) = 1.0 / 6.0;
  a_acof(5, 6, 8) = 0;
  a_acof(5, 7, 1) = 0;
  a_acof(5, 7, 2) = 0;
  a_acof(5, 7, 3) = 0;
  a_acof(5, 7, 4) = 0;
  a_acof(5, 7, 5) = -1.0 / 8.0;
  a_acof(5, 7, 6) = 1.0 / 6.0;
  a_acof(5, 7, 7) = -1.0 / 8.0;
  a_acof(5, 7, 8) = 0;
  a_acof(5, 8, 1) = 0;
  a_acof(5, 8, 2) = 0;
  a_acof(5, 8, 3) = 0;
  a_acof(5, 8, 4) = 0;
  a_acof(5, 8, 5) = 0;
  a_acof(5, 8, 6) = 0;
  a_acof(5, 8, 7) = 0;
  a_acof(5, 8, 8) = 0;
  a_acof(6, 1, 1) = 0;
  a_acof(6, 1, 2) = 0;
  a_acof(6, 1, 3) = 0;
  a_acof(6, 1, 4) = 1.0 / 392.0;
  a_acof(6, 1, 5) = 0;
  a_acof(6, 1, 6) = 0;
  a_acof(6, 1, 7) = 0;
  a_acof(6, 1, 8) = 0;
  a_acof(6, 2, 1) = 0;
  a_acof(6, 2, 2) = 0;
  a_acof(6, 2, 3) = 0;
  a_acof(6, 2, 4) = -1.0 / 144.0;
  a_acof(6, 2, 5) = 0;
  a_acof(6, 2, 6) = 0;
  a_acof(6, 2, 7) = 0;
  a_acof(6, 2, 8) = 0;
  a_acof(6, 3, 1) = 0;
  a_acof(6, 3, 2) = 0;
  a_acof(6, 3, 3) = 0;
  a_acof(6, 3, 4) = 3.0 / 784.0;
  a_acof(6, 3, 5) = -11.0 / 7056.0;
  a_acof(6, 3, 6) = 0;
  a_acof(6, 3, 7) = 0;
  a_acof(6, 3, 8) = 0;
  a_acof(6, 4, 1) = 0;
  a_acof(6, 4, 2) = 0;
  a_acof(6, 4, 3) = 0;
  a_acof(6, 4, 4) = -283.0 / 2352.0;
  a_acof(6, 4, 5) = 403.0 / 2352.0;
  a_acof(6, 4, 6) = -1.0 / 8.0;
  a_acof(6, 4, 7) = 0;
  a_acof(6, 4, 8) = 0;
  a_acof(6, 5, 1) = 0;
  a_acof(6, 5, 2) = 0;
  a_acof(6, 5, 3) = 0;
  a_acof(6, 5, 4) = 1135.0 / 7056.0;
  a_acof(6, 5, 5) = 1165.0 / 2352.0;
  a_acof(6, 5, 6) = 1.0 / 2.0;
  a_acof(6, 5, 7) = 1.0 / 6.0;
  a_acof(6, 5, 8) = 0;
  a_acof(6, 6, 1) = 0;
  a_acof(6, 6, 2) = 0;
  a_acof(6, 6, 3) = 0;
  a_acof(6, 6, 4) = -47.0 / 1176.0;
  a_acof(6, 6, 5) = -5869.0 / 7056.0;
  a_acof(6, 6, 6) = -3.0 / 4.0;
  a_acof(6, 6, 7) = -5.0 / 6.0;
  a_acof(6, 6, 8) = -1.0 / 24.0;
  a_acof(6, 7, 1) = 0;
  a_acof(6, 7, 2) = 0;
  a_acof(6, 7, 3) = 0;
  a_acof(6, 7, 4) = 0;
  a_acof(6, 7, 5) = 1.0 / 6.0;
  a_acof(6, 7, 6) = 1.0 / 2.0;
  a_acof(6, 7, 7) = 1.0 / 2.0;
  a_acof(6, 7, 8) = 1.0 / 6.0;
  a_acof(6, 8, 1) = 0;
  a_acof(6, 8, 2) = 0;
  a_acof(6, 8, 3) = 0;
  a_acof(6, 8, 4) = 0;
  a_acof(6, 8, 5) = 0;
  a_acof(6, 8, 6) = -1.0 / 8.0;
  a_acof(6, 8, 7) = 1.0 / 6.0;
  a_acof(6, 8, 8) = -1.0 / 8.0;
  // 129 non-zero out of 384.;
  ;
  a_Sb(0) = -3.0 / 12.0;
  a_Sb(1) = -10.0 / 12.0;
  a_Sb(2) = 18.0 / 12.0;
  a_Sb(3) = -6.0 / 12.0;
  a_Sb(4) = 1.0 / 12.0;
  a_Sb(5) = 0;
}

//-----------------------------------------------------------------------
void CurvilinearInterface::varcoeff_noghost() {
  Farray d5(0, 8);
  float_sw4 w0;
  Farray local_acof(1, 6, 1, 8, 1, 8), local_ghcof(1, 6);
  int i, j, k;
  varcoeffs4(local_acof, local_ghcof, sbop_no_gp);
  acof_no_gp = local_acof;
  d5 = 0.0;
  // Use 5th divided difference to cancel the ghost point contribution
  d5(0) = -1.00;
  d5(1) = 5.00;
  d5(2) = -10.00;
  d5(3) = 10.00;
  d5(4) = -5.00;
  d5(5) = 1.00;
  w0 = 17.0 / 48.0;
  i = 1;
  k = 1;  // only depends on the coefficient a(1)
  for (j = 1; j <= 8; j++)
    acof_no_gp(i, j, k) = local_acof(i, j, k) + d5(j) / (4.0 * w0);

  // the coeff for all ghost points are zero (don't divided by them!)
  //   ghcof_no_gp = 0.0;

  // boundary normal derivative, not using ghost points
  // sb = (-25*f(1)/12 + 4*f(2) - 3*f(3) + 4*f(4)/3 - f(5)/4)/h(q);
  sbop_no_gp(0) = 0.0;
  sbop_no_gp(1) = -25.0 / 12.0;
  sbop_no_gp(2) = 4.0;
  sbop_no_gp(3) = -3.0;
  sbop_no_gp(4) = 4.0 / 3.0;
  sbop_no_gp(5) = -1.0 / 4.0;
}

//-----------------------------------------------------------------------
void CurvilinearInterface::dx_46(Farray &bop) {
  bop(1, 1) = -24.0 / 17.0;
  bop(1, 2) = 59.0 / 34.0;
  bop(1, 3) = -4.0 / 17.0;
  bop(1, 4) = -3.0 / 34.0;
  bop(1, 5) = 0;
  bop(1, 6) = 0;
  bop(2, 1) = -1.0 / 2.0;
  bop(2, 2) = 0;
  bop(2, 3) = 1.0 / 2.0;
  bop(2, 4) = 0;
  bop(2, 5) = 0;
  bop(2, 6) = 0;
  bop(3, 1) = 4.0 / 43.0;
  bop(3, 2) = -59.0 / 86.0;
  bop(3, 3) = 0;
  bop(3, 4) = 59.0 / 86.0;
  bop(3, 5) = -4.0 / 43.0;
  bop(3, 6) = 0;
  bop(4, 1) = 3.0 / 98.0;
  bop(4, 2) = 0;
  bop(4, 3) = -59.0 / 98.0;
  bop(4, 4) = 0;
  bop(4, 5) = 32.0 / 49.0;
  bop(4, 6) = -4.0 / 49.0;
}
