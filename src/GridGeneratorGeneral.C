
#include "EW.h"
#include "GridGeneratorGeneral.h"

//-----------------------------------------------------------------------
GridGeneratorGeneral::GridGeneratorGeneral( float_sw4 topo_zmax, bool always_new, 
                                            int grid_interpolation_order, float_sw4 zetaBreak )
   : GridGenerator( topo_zmax, always_new, grid_interpolation_order, zetaBreak )
{
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::generate_grid_and_met( EW *a_ew, int g,
                                                  Sarray& a_x, Sarray& a_y, Sarray& a_z,
                                                  Sarray& a_jac, Sarray& a_met, bool a_comm )
{
   int ncurv = a_ew->mNumberOfGrids-a_ew->mNumberOfCartesianGrids;
   if( m_always_new || ncurv > 1 )
      generate_grid_and_met_new( a_ew, g, a_x, a_y, a_z, a_jac, a_met );
   else
      generate_grid_and_met_old( a_ew, a_x, a_y, a_z, a_jac, a_met );
   if( a_comm )
   {
      a_ew->communicate_array(a_jac,g);
      a_ew->communicate_array(a_met,g);
   }
}

//-----------------------------------------------------------------------
bool GridGeneratorGeneral::grid_mapping( EW* a_ew, float_sw4 q, float_sw4 r, float_sw4 s, int g,
                                         float_sw4& x, float_sw4& y, float_sw4& z )
{
   int ncurv = a_ew->mNumberOfGrids-a_ew->mNumberOfCartesianGrids;
   if( m_always_new || ncurv > 1 )
   {
   // New mapping.
      float_sw4 h = a_ew->mGridSize[g];
      int       Nz= a_ew->m_global_nz[g];
      return grid_mapping_new( a_ew, q, r, s, g, x, y, z, h, Nz );
   }
   else
   {
   // Old mapping, always only the topmost grid.
      float_sw4 h = a_ew->mGridSize[a_ew->mNumberOfGrids-1];
      int     Nz  = a_ew->m_global_nz[a_ew->mNumberOfGrids-1];
      return grid_mapping_old( q, r, s, g, x, y, z, a_ew->mTopoGridExt, h, Nz );
   }
}

//-----------------------------------------------------------------------
bool GridGeneratorGeneral::inverse_grid_mapping( EW* a_ew, float_sw4 x, float_sw4 y, float_sw4 z, int g,
                                                 float_sw4& q, float_sw4& r, float_sw4& s )
{
   int ncurv = a_ew->mNumberOfGrids-a_ew->mNumberOfCartesianGrids;
   if( m_always_new || ncurv > 1 )
   {
   // New mapping.
      float_sw4 h = a_ew->mGridSize[g];
      int     Nz  = a_ew->m_global_nz[g];
      return inverse_grid_mapping_new( a_ew, x, y, z, g, q, r, s, h, Nz );
   }
   else
   {
   // Old mapping, always only the topmost grid.
      float_sw4 h = a_ew->mGridSize[a_ew->mNumberOfGrids-1];
      int     Nz  = a_ew->m_global_nz[a_ew->mNumberOfGrids-1];

      float_sw4 bbox[6];
      a_ew->getGlobalBoundingBox( bbox );
      // 0. Check z
      if( !( bbox[4]-3*h < z && z < bbox[5]+3*h ) )
         return false;
      return inverse_grid_mapping_old( a_ew, x, y, z, g, q, r, s, a_ew->mTopoGridExt, h, Nz );
   }
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::grid_mapping_diff(
                             EW* a_ew, float_sw4 q, float_sw4 r, float_sw4 s, int g, 
                             int ic, int jc, int kc,
                             float_sw4& zq, float_sw4& zr, float_sw4& zs,
                             float_sw4& zqq, float_sw4& zqr, float_sw4& zqs,
                             float_sw4& zrr, float_sw4& zrs, float_sw4& zss )
                             
{
   int ncurv = a_ew->mNumberOfGrids-a_ew->mNumberOfCartesianGrids;
   if( m_always_new || ncurv > 1 )
   {
   // New mapping.
      float_sw4 h = a_ew->mGridSize[g];
      int      Nz = a_ew->m_global_nz[g];
      return grid_mapping_diff_new( a_ew, q, r, s, g, ic, jc, kc,
                                    zq, zr, zs, zqq, zqr, zqs, zrr, zrs, zss,
                                    h, Nz );
   }
   else
   {
   // Old mapping, always only the topmost grid.
      float_sw4 h = a_ew->mGridSize[a_ew->mNumberOfGrids-1];
      int     Nz  = a_ew->m_global_nz[a_ew->mNumberOfGrids-1];
      return grid_mapping_diff_old( a_ew, q, r, s, g, ic, jc, kc,
                                    zq, zr, zs, zqq, zqr, zqs, zrr, zrs, zss,
                                    a_ew->mTopoGridExt, h, Nz );
   }
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::generate_grid_and_met_old( EW *a_ew, Sarray& a_x, Sarray& a_y, Sarray& a_z,
                                           Sarray& a_jac, Sarray& a_met )
{
// Curvilinear grid that smoothly transitions to Cartesian at bottom.
// Single top curvilinear grid only.
   int g  = a_ew->mNumberOfGrids - 1;
   int nz = a_ew->m_global_nz[g];
   float_sw4 h=a_ew->mGridSize[g];

   float_sw4 izb = 1.0/(m_zetaBreak*(nz-1));
#pragma omp parallel for
   for (int k=a_x.m_kb; k<=a_x.m_ke; k++)
   {
      float_sw4 s  = (k-1)*izb;
      float_sw4 omsm = (1-s);
      for( int l=2 ; l <= m_grid_interpolation_order ; l++ )
         omsm *= (1-s);
      for (int j=a_x.m_jb; j<=a_x.m_je; j++)
         for (int i=a_x.m_ib; i<=a_x.m_ie; i++)
         {
            a_x(i,j,k) = (i-1)*h;
            a_y(i,j,k) = (j-1)*h;
            if( s >= 1 )
               a_z(i,j,k) = m_topo_zmax - (nz-k)*h;
            else
            {
               float_sw4 tau;
               tau = -m_curviInterface[0](i,j,1); 
               //               evaluate_topography(a_x(i,j,k),a_y(i,j,k),tau,m_topo);
               a_z(i,j,k) = m_topo_zmax - (nz-k)*h - omsm*(m_topo_zmax-(nz-1)*h+tau);
            }
         }
   }
   int ierr=0;
   ierr=metric_ci( a_x.m_ib, a_x.m_ie, a_x.m_jb, a_x.m_je, a_x.m_kb, a_x.m_ke, a_x.c_ptr(), 
                         a_y.c_ptr(), a_z.c_ptr(), a_met.c_ptr(), a_jac.c_ptr());
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::generate_grid_and_met_new( EW *a_ew, int g, Sarray& a_x, Sarray& a_y, Sarray& a_z,
                                           Sarray& a_jac, Sarray& a_met )
{
   int ncg=a_ew->mNumberOfCartesianGrids;
   float_sw4 h = a_ew->mGridSize[g]; 
   float_sw4 Nz_real = static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
   float_sw4 iNz_real = 1.0/Nz_real;
   float_sw4 scaleRatio=0;
   if( g > ncg )
      scaleRatio = (m_topo_zmax - a_ew->m_curviRefLev[g-1-ncg])/
                   (m_topo_zmax - a_ew->m_curviRefLev[g-ncg]);

#pragma omp parallel for
   for (int j=a_x.m_jb; j<=a_x.m_je; j++)
      for (int i=a_x.m_ib; i<=a_x.m_ie; i++)
      {
         float_sw4 X0  = (i-1)*h;
         float_sw4 Y0  = (j-1)*h;
         float_sw4 Ztop = m_curviInterface[g-ncg](i,j,1);
         float_sw4 Zbot = scaleRatio*Ztop+(1-scaleRatio)*m_topo_zmax;
#pragma omp simd
         for (int k=a_x.m_kb; k <= a_x.m_ke; k++)
         {
// Linear interpolation in the vertical direction
            float_sw4 zeta = static_cast<float_sw4>((k - a_ew->m_kStartInt[g])*iNz_real);
            a_x(i,j,k) = X0;
            a_y(i,j,k) = Y0;
            a_z(i,j,k) = (1.0- zeta)*Ztop + zeta*Zbot;
         }
      }

// Compute metric
   int ierr=0;
   ierr=metric_ci( a_x.m_ib, a_x.m_ie, a_x.m_jb, a_x.m_je, a_x.m_kb, a_x.m_ke, a_x.c_ptr(), 
                   a_y.c_ptr(), a_z.c_ptr(), a_met.c_ptr(), a_jac.c_ptr());
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::assignInterfaceSurfaces( EW* a_ew, Sarray& TopoGridExt )
{
   int ng=a_ew->mNumberOfGrids;
   int ncg=a_ew->mNumberOfCartesianGrids;
   m_curviInterface.resize(ng-ncg);
   int imin = TopoGridExt.m_ib;
   int imax = TopoGridExt.m_ie;
   int jmin = TopoGridExt.m_jb;
   int jmax = TopoGridExt.m_je;
   m_curviInterface[ng-1-ncg].define( imin, imax, jmin, jmax, 1, 1 );
   for (int i = imin; i <= imax; ++i)
      for (int j = jmin; j <= jmax; ++j)
      {
         m_curviInterface[ng-1-ncg](i,j,1) = -TopoGridExt(i,j,1);
      }
   int extragh = a_ew->m_iStart[ng-1]    - TopoGridExt.m_ib; // Number of extra ghost points
   int egh     = a_ew->m_iStartInt[ng-1] - TopoGridExt.m_ib; // Total number of ghost points (including extra points)
   int refFact=1;
   for( int g=ng-2 ; g>= ncg ; g-- )
   {
      int icib = a_ew->m_iStart[g]-extragh;
      int icie = a_ew->m_iEnd[g]  +extragh;
      int icjb = a_ew->m_jStart[g]-extragh;
      int icje = a_ew->m_jEnd[g]  +extragh;
      m_curviInterface[g-ncg].define(icib,icie,icjb,icje,1,1);
      refFact *= 2;
      float_sw4 scaleFact = (m_topo_zmax - a_ew->m_curviRefLev[g-ncg])/m_topo_zmax;

      // Interior only:
      for (int i=a_ew->m_iStartInt[g]; i<= a_ew->m_iEndInt[g]; i++)
         for (int j=a_ew->m_jStartInt[g]; j<= a_ew->m_jEndInt[g]; j++)
         {
            int iFine = 1 + (i-1)*refFact;
            int jFine = 1 + (j-1)*refFact;
            m_curviInterface[g-ncg](i,j,1) = scaleFact     * m_curviInterface[ng-1-ncg](iFine, jFine, 1) +
                                          (1.0 - scaleFact)* m_topo_zmax;
         }

      // Extrapolate to ghost points at domain boundaries:
      if (a_ew->m_jStartInt[g] == 1)
      {
         for( int i=icib+egh ; i <= icie-egh ; i++ )
            for( int q = 0 ; q < egh ; q++ )
               m_curviInterface[g-ncg](i,icjb+q,1)  = m_curviInterface[g-ncg](i,icjb+egh,1);
      }
      if (a_ew->m_jEndInt[g] == a_ew->m_global_ny[g])
      {
         for( int i=icib+egh ; i <= icie-egh ; i++ )
            for( int q = 0 ; q < egh ; q++ )
               m_curviInterface[g-ncg](i,icje-q,1)  = m_curviInterface[g-ncg](i,icje-egh,1);
      }
      if (a_ew->m_iStartInt[g] == 1)
      {
         for( int j=icjb ; j <= icje ; j++ )
            for( int q = 0 ; q < egh ; q++ )
               m_curviInterface[g-ncg](icib+q,j,1)  = m_curviInterface[g-ncg](icib+egh,j,1);
      }
      if (a_ew->m_iEndInt[g] == a_ew->m_global_nx[g])
      {
         for( int j=icjb ; j <= icje ; j++ )
            for( int q = 0 ; q < egh ; q++ )
               m_curviInterface[g-ncg](icie-q,j,1)  = m_curviInterface[g-ncg](icie-egh,j,1);
      }

      // Communicate padding points at processor boundaries
      a_ew->communicate_array_2d_isurf( m_curviInterface[g-ncg], g-ncg ); // Note: this routine adds ncg to its second argument.
   }
}

//-----------------------------------------------------------------------
bool GridGeneratorGeneral::grid_mapping_old( float_sw4 q, float_sw4 r, float_sw4 s, int g,
                                             float_sw4& x, float_sw4& y, float_sw4& z,
                                             Sarray& TopoGridExt, float_sw4 h, int Nz )
{
//
// Return (x,y) corresponding to (q,r).
// Return true and assign (X0,Y0,Z0) corresponding to (q,r,s) if (q,r) is inside processor.
//
// Returns false if 
// 1) (q,r,s) is outside the global parameter domain (expanded by ghost points)
// 2) There is no curvilinear grid. Still compute (X0, Y0) based on the top Cartesian grid size
// 3) (q,r) is not on this processor. Still compute (X0, Y0)
//
// The parameters are normalized such that 1 <= q <= Nx is the full domain (without ghost points),
//  1 <= r <= Ny, 1 <= s <= Nz.

   int nghost = 3; // 3 ghost points is maximum in solver.
// 0. Check if s is in range
   if( !( 1-nghost <= s && s <= Nz+nghost ) )
      return false;
   
// 1. Compute (x,y)
   x = (q-1.0)*h;
   y = (r-1.0)*h;


// 2. Compute z

// 2a. Find topography at (q,r), tau=tau(q,r)
   // Nearest grid point:
   int iNear = static_cast<int>(round(q));
   int jNear = static_cast<int>(round(r));
   float_sw4 tau;
   if ( fabs(iNear-q) < 1.e-9 && fabs(jNear-r) < 1.e-9 )
   {
// At a grid point, evaluate topography at that point
      if( TopoGridExt.in_range(1,iNear,jNear,1) )
         tau = TopoGridExt(iNear,jNear,1);
      else
         return false;
   }
   else
   {
  // Not at a grid  point, interpolate the topography
   // Nearest lower grid point
      iNear = static_cast<int>(floor(q));
      jNear = static_cast<int>(floor(r));
      if( TopoGridExt.in_range(1,iNear-3,jNear-3,1) &&  TopoGridExt.in_range(1,iNear+4,jNear+4,1) )
      {
	 float_sw4 a6cofi[8], a6cofj[8];
	 gettopowgh( q-iNear, a6cofi );
	 gettopowgh( r-jNear, a6cofj );
         tau = 0;
         for( int l=-3 ; l <= 4 ; l++ )
	    for( int k=-3 ; k <= 4 ; k++ )
	       tau += a6cofi[k+3]*a6cofj[l+3]*TopoGridExt(k+iNear,l+jNear,1);
      }
      else
      {
         return false;
      }
   } 
// 2b. Evaluate z-mapping
//   int Nz = m_global_nz[gFinest];
   z = m_topo_zmax - (Nz-s)*h;
   if( s-1 < m_zetaBreak*(Nz-1) )
   {
      float_sw4 omra = 1-(s-1)/(m_zetaBreak*(Nz-1));
      float_sw4 omsm = omra;
      for( int l=2 ; l <= m_grid_interpolation_order ; l++ )
         omsm *= omra;
      z -= omsm*(m_topo_zmax-(Nz-1)*h + tau);
   }
  return true;
}

//-----------------------------------------------------------------------
bool GridGeneratorGeneral::inverse_grid_mapping_old( EW* a_ew, 
                                                     float_sw4 x, float_sw4 y, float_sw4 z, int g,
                                                     float_sw4& q, float_sw4& r, float_sw4& s,
                                                     Sarray& TopoGridExt, float_sw4 h, int Nz )
{
//
// If (X0, Y0, Z0) is on the curvilinear grid and (X0, Y0) is on this processor:
// Return true and assigns (q,r,s) corresponding to point (mX0,mY0,mZ0)
//
// Normalize parameters such that 1 <= q <= Nx is the full domain (without ghost points),
//  1 <= r <= Ny, 1 <= s <= Nz.
//
//  int gCurv   = mNumberOfGrids - 1;
//  float_sw4 h = mGridSize[gCurv];


 // 1. Compute q and r
   q = x/h + 1.0;
   r = y/h + 1.0;
   int i= static_cast<int>(floor(q));
   int j= static_cast<int>(floor(r));   
   if( a_ew->interior_point_in_proc( i, j, g ) )
   {
      s = 0.;   
// 2. Compute s
      float_sw4 zlim = m_topo_zmax - (Nz-1)*(1-m_zetaBreak)*h;
      if( z >= zlim )
      {
// 2a. If z is in the Cartesian part of grid, this is the s value:
         s = (z-m_topo_zmax)/h + Nz;
      }
      else
      {
 // z is in curvilinear part of grid. 
// 2b. Find topography at (q,r), tau=tau(q,r)
   // Nearest grid point:
         int iNear = static_cast<int>(round(q));
         int jNear = static_cast<int>(round(r));
         float_sw4 tau;
         if ( fabs(iNear-q) < 1.e-9 && fabs(jNear-r) < 1.e-9 )
         {
// At a grid point, evaluate topography at that point
            if( TopoGridExt.in_range(1,iNear,jNear,1) )
               tau = TopoGridExt(iNear,jNear,1);
            else
               return false;
         }
         else
         {
            // Not at a grid  point, interpolate the topography
            // Nearest lower grid point
            iNear = static_cast<int>(floor(q));
            jNear = static_cast<int>(floor(r));
            if( TopoGridExt.in_range(1,iNear-3,jNear-3,1) &&  TopoGridExt.in_range(1,iNear+4,jNear+4,1) )
            {
               float_sw4 a6cofi[8], a6cofj[8];
               gettopowgh( q-iNear, a6cofi );
               gettopowgh( r-jNear, a6cofj );
               tau = 0;
               for( int l=-3 ; l <= 4 ; l++ )
                  for( int k=-3 ; k <= 4 ; k++ )
                     tau += a6cofi[k+3]*a6cofj[l+3]*TopoGridExt(k+iNear,l+jNear,1);
            }
            else
            {
               return false;
            }
         } 
 // 2c. Invert grid mapping to find s from z, i.e., solve Z(q,r,s)=z for s.
      // Invert polynomial, by Newton iteration
      // Use Cartesian value as initial guess.
         s = (z-m_topo_zmax)/h + Nz;
         float_sw4 z0  = m_topo_zmax - (Nz-1)*h + tau;
         float_sw4 izb = 1.0/(m_zetaBreak*(Nz-1));
         float_sw4 tol = 1e-12;
         float_sw4 er = tol+1;
         int maxit = 10;
         int it = 0;
         while( er > tol && it < maxit )
         {
            float_sw4 omra = 1-(s-1)*izb;
            float_sw4 omsm = omra;
            for( int l=2 ; l <= m_grid_interpolation_order-1 ; l++ )
               omsm *= omra;
            float_sw4 dfcn  = h + izb*m_grid_interpolation_order*omsm*z0;
            omsm *= omra;
            float_sw4 fcn = m_topo_zmax - (Nz-s)*h - omsm*z0 - z;
            float_sw4 sp    = s - fcn/dfcn;
            er    = abs(sp-s);
            s     = sp;
            it++;
         }
         if( er > tol )
         {
            cout << "invert_curvilinear_grid_mapping: poor convergence for X0, Y0, Z0 = " << x << ", " << y << ", " << z << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }
   // Successful if we get this far.
      return true;
   }
   else
      return false;
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::grid_mapping_diff_old(
                             EW* a_EW, float_sw4 q, float_sw4 r, float_sw4 s, int g, 
                             int ic, int jc, int kc,
                             float_sw4& zq, float_sw4& zr, float_sw4& zs,
                             float_sw4& zqq, float_sw4& zqr, float_sw4& zqs,
                             float_sw4& zrr, float_sw4& zrs, float_sw4& zss,
                             Sarray& TopoGridExt, float_sw4 h, int Nz )
{
   bool analytic_derivative = true;
   float_sw4 ai=q-ic;
   float_sw4 bi=r-jc;
   float_sw4 ci=s-kc;
   float_sw4 a6cofi[8], a6cofj[8], a6cofk[8];
   float_sw4 d6cofi[8], d6cofj[8], d6cofk[8];
   float_sw4 dd6cofi[8], dd6cofj[8], dd6cofk[8];
   float_sw4 ddd6cofi[8], ddd6cofj[8], ddd6cofk[8];
   getmetwgh( ai, a6cofi, d6cofi, dd6cofi, ddd6cofi );
   getmetwgh( bi, a6cofj, d6cofj, dd6cofj, ddd6cofj );
   getmetwgh( ci, a6cofk, d6cofk, dd6cofk, ddd6cofk );

   zq = zr = zs = 0;

   float_sw4 zpar = (s-1)/(m_zetaBreak*(Nz-1));
   float_sw4 kBreak = 1 + m_zetaBreak*(Nz-1);

   if( zpar >= 1 )
   {
      zq = 0;
      zr = 0;
      zs = h;
      zqq = zqr = zqs = zrr = zrs = zss = 0;
   }
   else 
   {
      int order=m_grid_interpolation_order;

      float_sw4 pp    = pow(1-zpar,order-1);
      float_sw4 powo  = (1-zpar)*pp;
      float_sw4 dpowo = -order*pp/m_zetaBreak;
      float_sw4 tauavg= 0;
      float_sw4 tauq=0, taur=0;
      float_sw4 tauqq=0, tauqr=0, taurr=0;
      for( int j=jc-3; j <= jc+4 ; j++ )
         for( int i=ic-3; i <= ic+4 ; i++ )
         {
            tauavg += a6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*TopoGridExt(i,j,1);
            tauq  +=  d6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*TopoGridExt(i,j,1);
            taur  +=  a6cofi[i-(ic-3)]* d6cofj[j-(jc-3)]*TopoGridExt(i,j,1);
            tauqq += dd6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*TopoGridExt(i,j,1);
            tauqr +=  d6cofi[i-(ic-3)]* d6cofj[j-(jc-3)]*TopoGridExt(i,j,1);
            taurr +=  a6cofi[i-(ic-3)]*dd6cofj[j-(jc-3)]*TopoGridExt(i,j,1);
         }
      zq  = (-tauq)*powo;
      zr  = (-taur)*powo;
      zqq = (-tauqq)*powo;
      zqr = (-tauqr)*powo;
      zrr = (-taurr)*powo;
      zqs = (-tauq)*dpowo/(Nz-1);
      zrs = (-taur)*dpowo/(Nz-1);

      float_sw4 zMax = a_EW->m_zmin[a_EW->mNumberOfCartesianGrids-1] - (Nz-kBreak)*h;
      float_sw4 c1   = zMax + tauavg - h*(kBreak-1);

      // Divide by Nz-1 to make consistent with undivided differences
      if( analytic_derivative )
      {
         zs  = h + c1*(-dpowo)/(Nz-1);
         zss = -c1*order*(order-1)*pow(1-zpar,order-2)/(m_zetaBreak*m_zetaBreak*(Nz-1)*(Nz-1));
      }
      else
      {
         zs = 0;
         zss= 0;
         float_sw4 z1d = 0;
         for( int k = kc-3 ; k <= kc+4; k++ ) 
         {
            zpar = (k-1)/(m_zetaBreak*(Nz-1));
            if( zpar >= 1 )
               z1d = zMax + (k-kBreak)*h;
            else
            {
               z1d = (1-zpar)*(-tauavg) + zpar*(zMax + c1*(1-zpar));
               for( int o=2 ; o < order ; o++ )
                  z1d += zpar*c1*pow(1-zpar,o);
            }
            zs  += d6cofk[k-(kc-3)]*z1d;
            zss += dd6cofk[k-(kc-3)]*z1d;
         }
      }
   }
}

//-----------------------------------------------------------------------
bool GridGeneratorGeneral::grid_mapping_new( EW* a_ew, float_sw4 q, float_sw4 r, float_sw4 s,
                                             int g,
                                             float_sw4& x, float_sw4& y, float_sw4& z,
                                             float_sw4 h, int Nz )
{
//
// Return (x,y) corresponding to (q,r).
// Return true and assign (x,y,z) corresponding to (q,r,s) if (q,r) is inside processor.
//
// Returns false if 
// 1) (q,r,s) is outside the global parameter domain (expanded by ghost points)
// 2) There is no curvilinear grid. Still compute (X0, Y0) based on the top Cartesian grid size
// 3) (q,r) is not on this processor. Still compute (X0, Y0)
//
// The parameters are normalized such that 1 <= q <= Nx is the full domain (without ghost points),
//  1 <= r <= Ny, 1 <= s <= Nz.

   int nghost = 3; // 3 ghost points is maximum in solver.
// 0. Check if s is in range
   if( !( 1-nghost <= s && s <= Nz+nghost ) )
      return false;
   
// 1. Compute (x,y)
   x = (q-1.0)*h;
   y = (r-1.0)*h;

// 2. Compute z
   int iSurfTop = g - a_ew->mNumberOfCartesianGrids;
   int iSurfBot = iSurfTop - 1;
   float_sw4 h0 = 2.0*h;
   float_sw4 Nz_real = static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
   float_sw4 iNz_real = 1.0/Nz_real;

   int i = static_cast<int>( floor(q) );
   int j = static_cast<int>( floor(r) );   

   float_sw4 Ztop = 0;
   float_sw4 a6cofi[8], a6cofj[8];
   gettopowgh( q-i, a6cofi );
   gettopowgh( r-j, a6cofj );
   for( int l=-3 ; l <= 4 ; l++ )
      for( int k=-3 ; k <= 4 ; k++ )
         Ztop += a6cofi[k+3]*a6cofj[l+3]*m_curviInterface[iSurfTop](k+i,l+j,1);

   float_sw4 Zbot;
   if (iSurfBot < 0)
   {
// Bottom interface of g=mNumberOfCartesianGrids is flat with z=m_topo_zmax
      Zbot = m_topo_zmax;
   }
   else
   {
// Bottom interface is non-planar (curvilinear)
      int iLow = static_cast<int>( floor(x/h0) )+1 ;
      int jLow = static_cast<int>( floor(y/h0) )+1;
                 
      float_sw4 xPt = (iLow-1)*h0;
      float_sw4 yPt = (jLow-1)*h0;

// First check if we are very close to a grid point
      if( fabs((xPt-x)/h0) < 1.e-9 && fabs((yPt-y)/h0) < 1.e-9 )
         Zbot = m_curviInterface[iSurfBot](iLow, jLow, 1);
      else
      {  // high order interpolation to get intermediate value of zBot
         if( true ) // point_in_proc_ext(i-3,j-3,gFinest) && point_in_proc_ext(i+4,j+4,gFinest)
         {
                  /* gettopowgh( q-i, a6cofi ); */
                  /* gettopowgh( r-j, a6cofj ); */
                  /* Zbot = 0; */
                  /* for( int l=j-3 ; l <= j+4 ; l++ ) */
                  /*    for( int k=i-3 ; k <= i+4 ; k++ ) */
                  /*       Zbot += a6cofi[k-i+3]*a6cofj[l-j+3]*m_curviInterface[iSurfBot](k,l,1); */
                  // for the purpose of plotting the grid, it suffices with linear interpolation
            float_sw4 xi  = (x - xPt)/h0;
            float_sw4 eta = (y - yPt)/h0;
            Zbot =
                     xi*eta       *(m_curviInterface[iSurfBot](iLow+1,jLow+1,1)) +
               (1.0-xi)*(1.0-eta) *(m_curviInterface[iSurfBot](iLow,jLow,1)) +
                     xi*(1.0-eta) *(m_curviInterface[iSurfBot](iLow+1,jLow,1)) +
               (1.0-xi)*eta       *(m_curviInterface[iSurfBot](iLow,jLow+1,1));
         }
      }
   }
// Linear interpolation in the vertical direction
   float_sw4 zeta = static_cast<float_sw4>((s - a_ew->m_kStartInt[g])*iNz_real);
   z = (1.0- zeta)*Ztop + zeta*Zbot;
   return true;
}

//-----------------------------------------------------------------------
bool GridGeneratorGeneral::inverse_grid_mapping_new( EW* a_ew, float_sw4 x,
                                                     float_sw4 y, float_sw4 z,
                                                     int g,
                                                     float_sw4& q, float_sw4& r,
                                                     float_sw4& s,
                                                     float_sw4 h, int Nz )
{
//
// If (X0, Y0, Z0) is in the curvilinear grid g and (X0, Y0) is on this processor:
// Return true and assigns (q,r,s) corresponding to point (mX0,mY0,mZ0)
//
// Normalize parameters such that 1 <= q <= Nx is the full domain (without ghost points),
//  1 <= r <= Ny, 1 <= s <= Nz.
//
//  int gCurv   = mNumberOfGrids - 1;
//  float_sw4 h = mGridSize[gCurv];
   bool retval = false;
 // 1. Compute q and r
   q = x/h + 1.0;
   r = y/h + 1.0;
   int i= static_cast<int>(floor(q));
   int j= static_cast<int>(floor(r));   
   if( a_ew->interior_point_in_proc( i, j, g ) )
   {
// 2. Compute s
      s = 0.;
      int grel =g-a_ew->mNumberOfCartesianGrids;
      float_sw4 ztop;
      // Find ztop at (x,y)
      if( fabs(x-(i-1)*h) < 1.e-9*h  && fabs(y-(j-1)*h) < 1.e-9*h )
      {
         ztop = m_curviInterface[grel](i, j, 1);
      }
      else
      {
         if( true )
         {
            // Use same interpolation order as for interpolate_topography.
            float_sw4 a6cofi[8], a6cofj[8];
            gettopowgh( q-i, a6cofi );
            gettopowgh( r-j, a6cofj );
            ztop = 0;
            for( int l=-3 ; l <= 4 ; l++ )
               for( int m=-3 ; m <= 4 ; m++ )
                  ztop += a6cofi[m+3]*a6cofj[l+3]*m_curviInterface[grel](m+i,l+j,1);
         }
         else
         {
            // Use bilinear interpolation for compatibility with lower interfaces
            float_sw4 xi  = (x - (i-1)*h)/(h);
            float_sw4 eta = (y - (j-1)*h)/(h);
            ztop =
               xi*eta             *(m_curviInterface[grel](i+1,j+1,1)) +
               (1.0-xi)*(1.0-eta) *(m_curviInterface[grel](i,j,1)) +
               xi*(1.0-eta)       *(m_curviInterface[grel](i+1,j,1)) +
               (1.0-xi)*eta       *(m_curviInterface[grel](i,j+1,1));

         }
      }

      // Find zbot at (x,y)
      float_sw4 zbot, hc=2*h;
      if( g==a_ew->mNumberOfCartesianGrids )
      {
         zbot = m_topo_zmax;
      }
      else
      {
         int ic= static_cast<int>(floor(x/(hc)+1));
         int jc= static_cast<int>(floor(y/(hc)+1));   
         if( fabs(x-(ic-1)*hc) < 1.e-9*hc  && fabs( y-(jc-1)*hc) < 1.e-9*hc )
            zbot = m_curviInterface[grel-1](ic, jc, 1);
         else
         {
            // Use same interpolation order as for interpolate_topography.
            float_sw4 xi  = (x - (ic-1)*hc)/(hc);
            float_sw4 eta = (y - (jc-1)*hc)/(hc);
            if( true )
            {
               float_sw4 a6cofi[8], a6cofj[8];
               gettopowgh( xi, a6cofi );
               gettopowgh( eta, a6cofj );
               zbot = 0;
               for( int l=-3 ; l <= 4 ; l++ )
                  for( int m=-3 ; m <= 4 ; m++ )
                     zbot += a6cofi[m+3]*a6cofj[l+3]*m_curviInterface[grel-1](m+ic,l+jc,1);
            }
            else
            {
            // Linear interpolation to get intermediate value of zbot
               zbot =
               xi*eta             *(m_curviInterface[grel-1](ic+1,jc+1,1)) +
               (1.0-xi)*(1.0-eta) *(m_curviInterface[grel-1](ic,jc,1)) +
               xi*(1.0-eta)       *(m_curviInterface[grel-1](ic+1,jc,1)) +
               (1.0-xi)*eta       *(m_curviInterface[grel-1](ic,jc+1,1));
            }
         }
      }
      // Allow coordinate a little bit above topography
      if( ztop-1e-6*h <= z && z <= zbot+1e-6*h || ( g == a_ew->mNumberOfGrids-1 && (ztop-h*0.5 <= z && z <= zbot)))
      {
            // Point is found on grid:
         s = (z-ztop)/(zbot-ztop)*(Nz-1)+1;
         //         if( s < 1 )
         //            s = 1;
         retval = true;
      }
   }
   return retval;
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::grid_mapping_diff_new(
                                 EW* a_ew, float_sw4 q, float_sw4 r, float_sw4 s, int g, 
                                 int ic, int jc, int kc,
                                 float_sw4& zq, float_sw4& zr, float_sw4& zs,
                                 float_sw4& zqq, float_sw4& zqr, float_sw4& zqs,
                                 float_sw4& zrr, float_sw4& zrs, float_sw4& zss,
                                 float_sw4 h, int Nz )
{
   // Computes derivatives of the grid mapping z=z(q,r,s) at the given location (q,r,s).
   //
   // Input: (q,r,s)    - Location in mapped space
   //        (ic,jc,kc) - Center stencils around this grid point 
   //                     It is assumed that (ic,jc) is in the interior of this processor.
   //        h          - Grid spacing
   //        Nz         - Number of (interior) grid points in the k-direction
   // Output: zq, zr, zs - Derivatives of the grid z-coordinate (grid mapping z=z(q,r,s) )
   //    zqq, zqr, zqs, zrr, zrs, zss - Second derivatives of z=z(q,r,s)
   //
   // The parameters (q,r,s) are normalized such that 1 <= q <= Nx is the
   // full domain (without ghost points). Similarly 1 <= r <= Ny, 1 <= s <= Nz.

   int nghost = 3; // 3 ghost points is maximum in solver.

   if( !( 1-nghost <= s && s <= Nz+nghost ) )
      return;
   
   float_sw4 ai=q-ic;
   float_sw4 bi=r-jc;
   float_sw4 a6cofi[8],   a6cofj[8];
   float_sw4 d6cofi[8],   d6cofj[8];
   float_sw4 dd6cofi[8],  dd6cofj[8];
   float_sw4 ddd6cofi[8], ddd6cofj[8];
   getmetwgh( ai, a6cofi, d6cofi, dd6cofi, ddd6cofi );
   getmetwgh( bi, a6cofj, d6cofj, dd6cofj, ddd6cofj );

   int iSurfTop = g - a_ew->mNumberOfCartesianGrids;
   int iSurfBot = iSurfTop - 1;
   float_sw4 x=(q-1)*h;
   float_sw4 y=(r-1)*h;   
   float_sw4 h0 = 2.0*h;
   float_sw4 Nz_real = static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
   float_sw4 iNz_real = 1.0/Nz_real;

   if( !(m_curviInterface[iSurfTop].in_range(1,ic-3,jc-3,1) && 
         m_curviInterface[iSurfTop].in_range(1,ic+4,jc+4,1)) )
      std::cout <<"ERROR in gridgen diff new, top " << ic << " " << jc << std::endl;
   //   std::cout << "in gridgen diff" << std::endl;

   float_sw4 Ztop = 0, Ztopq=0, Ztopr=0, Ztopqq=0, Ztopqr=0, Ztoprr=0;
   for( int j=jc-3; j <= jc+4 ; j++ )
      for( int i=ic-3; i <= ic+4 ; i++ )
      {
         Ztop   +=  a6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*m_curviInterface[iSurfTop](i,j,1);
         Ztopq  +=  d6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*m_curviInterface[iSurfTop](i,j,1);
         Ztopr  +=  a6cofi[i-(ic-3)]* d6cofj[j-(jc-3)]*m_curviInterface[iSurfTop](i,j,1);
         Ztopqq += dd6cofi[i-(ic-3)]* a6cofj[j-(jc-3)]*m_curviInterface[iSurfTop](i,j,1);
         Ztopqr +=  d6cofi[i-(ic-3)]* d6cofj[j-(jc-3)]*m_curviInterface[iSurfTop](i,j,1);
         Ztoprr +=  a6cofi[i-(ic-3)]*dd6cofj[j-(jc-3)]*m_curviInterface[iSurfTop](i,j,1);
      }

   float_sw4 Zbot=0, Zbotq=0, Zbotr=0, Zbotqq=0, Zbotqr=0, Zbotrr=0;
   if (iSurfBot < 0)
   {
// Bottom interface of g=mNumberOfCartesianGrids is flat with z=m_topo_zmax
      Zbot = m_topo_zmax;
   }
   else
   {
// Bottom interface is non-planar (curvilinear)
// Fix this formula for the coarser grid
//      int icc, jcc;
//      if( ic % 2 == 0 )
//         icc = ic/2;
//      else if( q-ic > 0 )
//         icc = (ic+1)/2;
//      else
//         icc = (ic-1)/2;
//      if( jc % 2 == 0 )
//         jcc = jc/2;
//      else if( r-jc > 0 )
//         jcc = (jc+1)/2;
//      else
//         jcc = (jc-1)/2;
      int icc = static_cast<int>( floor(x/h0) )+1 ;
      int jcc = static_cast<int>( floor(y/h0) )+1;
      //                 
      //      float_sw4 xPt = (iLow-1)*h0;
      //      float_sw4 yPt = (jLow-1)*h0;
      //
      float_sw4 xi  = x/h0 - icc+1;
      float_sw4 eta = y/h0 - jcc+1;

      getmetwgh( xi,  a6cofi, d6cofi, dd6cofi, ddd6cofi );
      getmetwgh( eta, a6cofj, d6cofj, dd6cofj, ddd6cofj );
      if( !(m_curviInterface[iSurfBot].in_range(1,icc-3,jcc-3,1) && 
            m_curviInterface[iSurfBot].in_range(1,icc+4,jcc+4,1)) )
         std::cout <<"ERROR in gridgen diff new, bot " << icc << " " << jcc << std::endl;
      for( int j=jcc-3; j <= jcc+4 ; j++ )
         for( int i=icc-3; i <= icc+4 ; i++ )
         {
            Zbot   +=  a6cofi[i-(icc-3)]* a6cofj[j-(jcc-3)]*m_curviInterface[iSurfBot](i,j,1);
            Zbotq  +=  d6cofi[i-(icc-3)]* a6cofj[j-(jcc-3)]*m_curviInterface[iSurfBot](i,j,1);
            Zbotr  +=  a6cofi[i-(icc-3)]* d6cofj[j-(jcc-3)]*m_curviInterface[iSurfBot](i,j,1);
            Zbotqq += dd6cofi[i-(icc-3)]* a6cofj[j-(jcc-3)]*m_curviInterface[iSurfBot](i,j,1);
            Zbotqr +=  d6cofi[i-(icc-3)]* d6cofj[j-(jcc-3)]*m_curviInterface[iSurfBot](i,j,1);
            Zbotrr +=  a6cofi[i-(icc-3)]*dd6cofj[j-(jcc-3)]*m_curviInterface[iSurfBot](i,j,1);
         }
      // Above derivatives are taken w.r.t. to the coarse grid parameter qc=(q+1)/2,.
      // Need to transform to derivatives w.r.t. q :
      Zbotq  *=0.5;
      Zbotr  *=0.5;      
      Zbotqq *=0.25;
      Zbotqr *=0.25;
      Zbotrr *=0.25;      
      //      int iLow = static_cast<int>( floor(x/h0) )+1 ;
      //      int jLow = static_cast<int>( floor(y/h0) )+1;
      //                 
      //      float_sw4 xPt = (iLow-1)*h0;
      //      float_sw4 yPt = (jLow-1)*h0;
      //
      //      float_sw4 xi  = (x - xPt)/h0;
      //      float_sw4 eta = (y - yPt)/h0;
      //      Zbot =
      //                     xi*eta       *(m_curviInterface[iSurfBot](iLow+1,jLow+1,1)) +
      //               (1.0-xi)*(1.0-eta) *(m_curviInterface[iSurfBot](iLow,jLow,1)) +
      //                     xi*(1.0-eta) *(m_curviInterface[iSurfBot](iLow+1,jLow,1)) +
      //               (1.0-xi)*eta       *(m_curviInterface[iSurfBot](iLow,jLow+1,1));
      //
      //      Zbotq = eta*( m_curviInterface[iSurfBot](iLow+1,jLow+1,1) -
      //                    m_curviInterface[iSurfBot](iLow,  jLow+1,1)  ) +
      //          (1-eta)*( m_curviInterface[iSurfBot](iLow+1,jLow,  1) -
      //                    m_curviInterface[iSurfBot](iLow,  jLow,  1)   );
      //      Zbotr = xi*( m_curviInterface[iSurfBot](iLow+1, jLow+1,1) -
      //                   m_curviInterface[iSurfBot](iLow+1, jLow,  1)  ) +
      //          (1-xi)*( m_curviInterface[iSurfBot](iLow,   jLow+1,1) -
      //                   m_curviInterface[iSurfBot](iLow,   jLow,  1)   );
      //      Zbotqr = m_curviInterface[iSurfBot](iLow+1,jLow+1,1) -
      //               m_curviInterface[iSurfBot](iLow,  jLow+1,1) - 
      //               m_curviInterface[iSurfBot](iLow+1,jLow,  1) +
      //               m_curviInterface[iSurfBot](iLow,  jLow,  1);
      //   }
   }
// Linear interpolation in the vertical direction
   float_sw4 zeta = static_cast<float_sw4>((s - a_ew->m_kStartInt[g])*iNz_real);
   zq = (1.0- zeta)*Ztopq + zeta*Zbotq;
   zr = (1.0- zeta)*Ztopr + zeta*Zbotr;   
   zs = (Zbot-Ztop)*iNz_real;
   zqq = (1.0- zeta)*Ztopqq + zeta*Zbotqq;
   zqr = (1.0- zeta)*Ztopqr + zeta*Zbotqr;   
   zrr = (1.0- zeta)*Ztoprr + zeta*Zbotrr;
   zqs = (Zbotq-Ztopq)*iNz_real;
   zrs = (Zbotr-Ztopr)*iNz_real;
   zss = 0;
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::getmetwgh( float_sw4 ai, float_sw4 wgh[8], float_sw4 dwgh[8],
			float_sw4 ddwgh[8], float_sw4 dddwgh[8] ) const
{
   float_sw4 pol = ai*ai*ai*ai*ai*ai*ai*(-251+135*ai+25*ai*ai-
                                      33*ai*ai*ai+6*ai*ai*ai*ai)/720;

   wgh[0] = -1.0/60*ai + 1.0/180*ai*ai + 1.0/48*ai*ai*ai + 23.0/144*ai*ai*ai*ai 
      - (17.0*ai + 223.0)*ai*ai*ai*ai*ai/720 - pol;
   wgh[1] = 3.0/20*ai -3.0/40*ai*ai -1.0/6*ai*ai*ai - 13.0/12*ai*ai*ai*ai + 
      97.0/45*ai*ai*ai*ai*ai + 1.0/6*ai*ai*ai*ai*ai*ai + 7*pol;
   wgh[2] = -0.75*ai +0.75*ai*ai+(13.0+155*ai)*ai*ai*ai/48 -103.0/16*ai*ai*ai*ai*ai
      - 121.0/240*ai*ai*ai*ai*ai*ai - 21*pol;
   wgh[3] = 1 - 49.0/36*ai*ai - 49.0/9*ai*ai*ai*ai+385.0/36*ai*ai*ai*ai*ai +
      61.0/72*ai*ai*ai*ai*ai*ai + 35*pol;
   wgh[4] = 0.75*ai + 0.75*ai*ai - 13.0/48*ai*ai*ai + 89.0/16*ai*ai*ai*ai - 
         1537.0/144*ai*ai*ai*ai*ai - 41.0/48*ai*ai*ai*ai*ai*ai - 35*pol;
   wgh[5] = -3.0/20*ai - 3.0/40*ai*ai + 1.0/6*ai*ai*ai - 41.0/12*ai*ai*ai*ai
      + 6.4*ai*ai*ai*ai*ai + 31.0/60*ai*ai*ai*ai*ai*ai + 21*pol;
   wgh[6] = 1.0/60*ai + 1.0/180*ai*ai - 1.0/48*ai*ai*ai + 167.0/144*ai*ai*ai*ai -
      1537.0/720*ai*ai*ai*ai*ai- 25.0/144*ai*ai*ai*ai*ai*ai - 7*pol;
   wgh[7] = -1.0/6*ai*ai*ai*ai + 11.0/36*ai*ai*ai*ai*ai + 1.0/40*ai*ai*ai*ai*ai*ai + pol;

   // Derivative wrt. ai
   pol = ai*ai*ai*ai*ai*ai*(-1757.0/720 + 1.5*ai + 0.31250*ai*ai - (1.375*ai*ai*ai-0.275*ai*ai*ai*ai)/3);
   dwgh[0] = -1.0/60 + 1.0/90*ai+ ai*ai/16 + 23.0/36*ai*ai*ai - 223.0/144*ai*ai*ai*ai -
      17.0/120*ai*ai*ai*ai*ai - pol;
   dwgh[1] = 3.0/20 - 3.0/20*ai - 0.5*ai*ai-13.0/3*ai*ai*ai + 97.0/9*ai*ai*ai*ai +
      ai*ai*ai*ai*ai + 7*pol;
   dwgh[2] = -0.75 + 1.5*ai + 13.0/16*ai*ai + 155.0*ai*ai*ai/12-103.0*5.0/16*ai*ai*ai*ai
      - 121.0/40*ai*ai*ai*ai*ai - 21*pol;
   dwgh[3] = -49.0/18*ai - 4*49.0/9.0*ai*ai*ai + 385.0*5.0/36*ai*ai*ai*ai +
      61.0/12*ai*ai*ai*ai*ai + 35*pol;
   dwgh[4] = 0.75 + 1.5*ai - 13.0/16*ai*ai + 89.0/4*ai*ai*ai - 1537.0*5/144.0*ai*ai*ai*ai -
      41.0/8*ai*ai*ai*ai*ai - 35*pol;
   dwgh[5] = -3.0/20 - 3.0/20*ai + 0.5*ai*ai-41.0/3*ai*ai*ai + 32*ai*ai*ai*ai +
      3.1*ai*ai*ai*ai*ai + 21*pol;
   dwgh[6] = 1.0/60 + 1.0/90*ai - 1.0/16*ai*ai + 167.0/36*ai*ai*ai - 1537.0/144*ai*ai*ai*ai -
      25.0/24*ai*ai*ai*ai*ai - 7*pol;
   dwgh[7] = -2.0/3*ai*ai*ai + 55.0/36*ai*ai*ai*ai + 3.0/20*ai*ai*ai*ai*ai + pol;

   // Second derivative wrt. ai
   pol = ai*ai*ai*ai*ai*(-1757.0/120 + 10.5*ai + 2.5*ai*ai - 4.125*ai*ai*ai + 11.0/12*ai*ai*ai*ai);
   ddwgh[0] = 1.0/90 + 0.125*ai + 23.0/12*ai*ai - 223.0/36*ai*ai*ai - 17.0/24*ai*ai*ai*ai - pol;
   ddwgh[1] = -3.0/20 - ai - 13.0*ai*ai + 4*97.0/9.0*ai*ai*ai + 5*ai*ai*ai*ai + 7*pol;
   ddwgh[2] = 1.5 + 13.0/8*ai + 155.0/4*ai*ai - 103.0*5.0/4*ai*ai*ai - 121.0/8*ai*ai*ai*ai - 21*pol;
   ddwgh[3] = -49.0/18 - 4*49.0/3.0*ai*ai + 385.0*5.0/9.0*ai*ai*ai + 5*61.0/12*ai*ai*ai*ai + 35*pol;
   ddwgh[4] = 1.5 -13.0/8*ai+89.0*3.0/4*ai*ai - 1537.0*5.0/36*ai*ai*ai - 205.0/8*ai*ai*ai*ai - 35*pol;
   ddwgh[5] = -3.0/20 + ai - 41.0*ai*ai + 128*ai*ai*ai + 15.5*ai*ai*ai*ai + 21*pol;
   ddwgh[6] = 1.0/90 - 0.125*ai + 167.0/12*ai*ai - 1537.0/36*ai*ai*ai - 125.0/24*ai*ai*ai*ai - 7*pol;
   ddwgh[7] = -2*ai*ai + 220.0/36*ai*ai*ai + 0.75*ai*ai*ai*ai + pol;

   // Third derivative wrt. ai
   pol = ai*ai*ai*ai*(-1757.0/24 + 63*ai + 17.5*ai*ai - 33*ai*ai*ai + 8.25*ai*ai*ai*ai);
   dddwgh[0] = 0.125 + 23.0/6*ai-223.0/12*ai*ai-17.0/6*ai*ai*ai - pol;
   dddwgh[1] = -1 - 26.0*ai + 4*97.0/3*ai*ai + 20*ai*ai*ai + 7*pol;
   dddwgh[2] =  1.625 + 77.5*ai - 386.25*ai*ai -60.5*ai*ai*ai - 21*pol;
   dddwgh[3] = -392.0/3*ai + 1925.0/3*ai*ai + 305.0/3*ai*ai*ai + 35*pol;
   dddwgh[4] = -1.625 + 133.5*ai-7685.0/12*ai*ai - 102.5*ai*ai*ai - 35*pol;
   dddwgh[5] = 1 - 82.0*ai + 384.0*ai*ai + 62.0*ai*ai*ai + 21*pol;
   dddwgh[6] = -0.125 + 167.0/6*ai - 1537.0/12*ai*ai - 125.0/6*ai*ai*ai - 7*pol;
   dddwgh[7] = -4*ai + 220.0/12*ai*ai + 3*ai*ai*ai + pol;
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::generate_z_and_j( EW* a_ew, int g, Sarray& z, Sarray& J )
{
   int ng=a_ew->mNumberOfGrids;
   int ncg=a_ew->mNumberOfCartesianGrids;
   int ref=1;
   for( int grid=ng-1 ; grid > g ; grid-- )
      ref *= 2;

   int iSurfTop = g - a_ew->mNumberOfCartesianGrids;
   int iSurfBot = iSurfTop - 1;
   float_sw4 h = a_ew->mGridSize[g]; 
   float_sw4 h0 = 2.0*h;
   float_sw4 Nz_real = static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
   float_sw4 iNz_real = 1.0/Nz_real;
   float_sw4 scaleFact=0, scaleRatio=0;
   if( g > ncg )
   {
      scaleFact = (m_topo_zmax - a_ew->m_curviRefLev[g-1-ncg])/m_topo_zmax;
      scaleRatio = (m_topo_zmax - a_ew->m_curviRefLev[g-1-ncg])/
         (m_topo_zmax - a_ew->m_curviRefLev[g-ncg]);
   }
   //   if( a_ew->getRank() == 4 )
   //   {
   //   std::cout << "zlims " << z.m_ib << " " << z.m_ie << " " << z.m_jb << " " << z.m_je << std::endl;
   //   std::cout << "curvii lims = " 
   //             << m_curviInterface[iSurfTop].m_ib << " "
   //             << m_curviInterface[iSurfTop].m_ie << " " 
   //             << m_curviInterface[iSurfTop].m_jb << " "
   //             << m_curviInterface[iSurfTop].m_je << std::endl;
   //   std::cout << "ref = " << ref << " curvitop lims = "
   //             << m_curviInterface[ng-1-ncg].m_ib << " "
   //             << m_curviInterface[ng-1-ncg].m_ie << " " 
   //             << m_curviInterface[ng-1-ncg].m_jb << " "
   //             << m_curviInterface[ng-1-ncg].m_je << std::endl;
   //   }
#pragma omp parallel for
   for (int j=z.m_jb; j<=z.m_je; j++)
      for (int i=z.m_ib; i<=z.m_ie; i++)
      {
         float_sw4 Ztop = m_curviInterface[iSurfTop](i,j,1);
         float_sw4 Zbot;
         if (iSurfBot < 0)
         {
// Bottom interface of g=mNumberOfCartesianGrids is flat with z=m_topo_zmax
            Zbot = m_topo_zmax;
         }
         else
         {
            //            Zbot = scaleFact*m_curviInterface[ng-1-ncg](ref*(i-1)+1,ref*(j-1)+1, 1) +
            //                          (1.0 - scaleFact)* m_topo_zmax;
            Zbot = scaleRatio*Ztop+(1-scaleRatio)*m_topo_zmax;
         }
#pragma omp parallel for
         for (int k=z.m_kb; k <= z.m_ke; k++)
         {
// Linear interpolation in the vertical direction
            float_sw4 zeta = static_cast<float_sw4>((k - a_ew->m_kStartInt[g])*iNz_real);
            z(i,j,k) = (1.0- zeta)*Ztop + zeta*Zbot;
            J(i,j,k) = h*h*iNz_real*(Zbot-Ztop); // (h*h*dz/dr)
            // note, mapping linear in k, exact and numerical derivatives are the same.
         }
      }
}
