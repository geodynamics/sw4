
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
   // Only old mapping supported this far. Always the topmost grid.
   float_sw4 h = a_ew->mGridSize[a_ew->mNumberOfGrids-1];
   int     Nz  = a_ew->m_global_nz[a_ew->mNumberOfGrids-1];
   return grid_mapping_old( q, r, s, g, x, y, z, a_ew->mTopoGridExt, h, Nz );
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
   ierr=a_ew->metric_ci( a_x.m_ib, a_x.m_ie, a_x.m_jb, a_x.m_je, a_x.m_kb, a_x.m_ke, a_x.c_ptr(), 
                         a_y.c_ptr(), a_z.c_ptr(), a_met.c_ptr(), a_jac.c_ptr());
}

//-----------------------------------------------------------------------
void GridGeneratorGeneral::generate_grid_and_met_new( EW *a_ew, int g, Sarray& a_x, Sarray& a_y, Sarray& a_z,
                                           Sarray& a_jac, Sarray& a_met )
{
   int iSurfTop = g - a_ew->mNumberOfCartesianGrids;
   int iSurfBot = iSurfTop - 1;
   float_sw4 h = a_ew->mGridSize[g]; 
   float_sw4 h0 = 2.0*h;
   float_sw4 Nz_real = static_cast<float_sw4>(a_ew->m_kEndInt[g] - a_ew->m_kStartInt[g]);
   float_sw4 iNz_real = 1.0/Nz_real;

#pragma omp parallel for
   for (int j=a_x.m_jb; j<=a_x.m_je; j++)
      for (int i=a_x.m_ib; i<=a_x.m_ie; i++)
      {
         float_sw4 X0  = (i-1)*h;
         float_sw4 Y0  = (j-1)*h;
         float_sw4 Ztop = m_curviInterface[iSurfTop](i,j,1);
         float_sw4 Zbot;
         if (iSurfBot < 0)
         {
// Bottom interface of g=mNumberOfCartesianGrids is flat with z=m_topo_zmax
            Zbot = m_topo_zmax;
         }
         else
         {
// Bottom interface is non-planar (curvilinear)
            int iLow = static_cast<int>( floor(X0/h0) )+1 ;
            int jLow = static_cast<int>( floor(Y0/h0) )+1;
                 
            float_sw4 xPt = (iLow-1)*h0;
            float_sw4 yPt = (jLow-1)*h0;

// First check if we are very close to a grid point
            if( fabs((xPt-X0)/h0) < 1.e-9 && fabs((yPt-Y0)/h0) < 1.e-9 )
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
                  float_sw4 xi  = (X0 - xPt)/h0;
                  float_sw4 eta = (Y0 - yPt)/h0;
                  Zbot =
                           xi*eta       *(m_curviInterface[iSurfBot](iLow+1,jLow+1,1)) +
                     (1.0-xi)*(1.0-eta) *(m_curviInterface[iSurfBot](iLow,jLow,1)) +
                           xi*(1.0-eta) *(m_curviInterface[iSurfBot](iLow+1,jLow,1)) +
                     (1.0-xi)*eta       *(m_curviInterface[iSurfBot](iLow,jLow+1,1));
               }
            }
         }
#pragma omp parallel for
         for (int k=a_x.m_kb; k <= a_x.m_ke; k++)
         {
// Linear interpolation in the vertical direction
            float_sw4 zeta = static_cast<float_sw4>((k - a_ew->m_kStartInt[g])*iNz_real);
            a_x(i,j,k) = X0;
            a_y(i,j,k) = Y0;
            a_z(i,j,k) = (1.0- zeta)*Ztop + zeta*Zbot;
         } 
      }
// make sure all processors have made their grid before we continue
//   a_ew->communicate_array( a_z, g ); 

// Compute metric
   int ierr=0;
   ierr=a_ew->metric_ci( a_x.m_ib, a_x.m_ie, a_x.m_jb, a_x.m_je, a_x.m_kb, a_x.m_ke, a_x.c_ptr(), 
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
   int i= static_cast<int>(round(q));
   int j= static_cast<int>(round(r));   
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
   int i= static_cast<int>(round(q));
   int j= static_cast<int>(round(r));   
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
         if( g == a_ew->mNumberOfGrids-1 )
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
         int ic= static_cast<int>(round(x/(hc)+1));
         int jc= static_cast<int>(round(y/(hc)+1));   
         if( fabs(x-(ic-1)*hc) < 1.e-9*hc  && fabs( y-(jc-1)*hc) < 1.e-9*hc )
            zbot = m_curviInterface[grel-1](ic, jc, 1);
         else
         {  // Linear interpolation to get intermediate value of zbot
            float_sw4 xi  = (x - (ic-1)*hc)/(hc);
            float_sw4 eta = (y - (jc-1)*hc)/(hc);
            zbot =
               xi*eta             *(m_curviInterface[grel-1](ic+1,jc+1,1)) +
               (1.0-xi)*(1.0-eta) *(m_curviInterface[grel-1](ic,jc,1)) +
               xi*(1.0-eta)       *(m_curviInterface[grel-1](ic+1,jc,1)) +
               (1.0-xi)*eta       *(m_curviInterface[grel-1](ic,jc+1,1));
         }
      }
      // Allow coordinate a little bit above topography
      if( ztop <= z && z <= zbot || ( g == a_ew->mNumberOfGrids-1 && (ztop-h*0.5 <= z && z <= zbot)))
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

