#include "mpi.h"

#include "EW.h"
#define SQR(x) ((x)*(x))
//---------------------------------------------------------
void EW::setup_metric()
{
  if (!m_topography_exists ) return;
  
  if (mVerbose >= 1 && proc_zero())
    cout << "***inside setup_metric***"<< endl;

  int g=mNumberOfGrids-1;
  
  int Bx=m_iStart[g];
  int By=m_jStart[g];
  int Bz=m_kStart[g];

  int Nx=m_iEnd[g];
  int Ny=m_jEnd[g];
  int Nz=m_kEnd[g];
  
  double xq,xr,xs,yq,yr,ys,zq,zr,zs;
  double minJ, maxJ;
  maxJ=0.;
  minJ=1.e18;
  for(int k=Bz; k<=Nz; k++)
    for(int j=By; j<=Ny; j++)
      for(int i=Bx; i<=Nx; i++)
      {
	if (i==Bx){
	  xq = -0.5*(3.*mX(Bx,j,k)-4.*mX(Bx+1,j,k)+mX(Bx+2,j,k));
	  yq = -0.5*(3.*mY(Bx,j,k)-4.*mY(Bx+1,j,k)+mY(Bx+2,j,k));
	  zq = -0.5*(3.*mZ(Bx,j,k)-4.*mZ(Bx+1,j,k)+mZ(Bx+2,j,k));
	} else if (i==Nx){
	  xq=0.5*(3.*mX(Nx,j,k)-4.*mX(Nx-1,j,k)+mX(Nx-2,j,k));
	  yq=0.5*(3.*mY(Nx,j,k)-4.*mY(Nx-1,j,k)+mY(Nx-2,j,k));
	  zq=0.5*(3.*mZ(Nx,j,k)-4.*mZ(Nx-1,j,k)+mZ(Nx-2,j,k));
	} else {
	  xq=0.5*(mX(i+1,j,k)-mX(i-1,j,k));
	  yq=0.5*(mY(i+1,j,k)-mY(i-1,j,k));
	  zq=0.5*(mZ(i+1,j,k)-mZ(i-1,j,k));
	}
	if (j==By){
	  xr = -0.5*(3.*mX(i,By,k)-4.*mX(i,By+1,k)+mX(i,By+2,k));
	  yr = -0.5*(3.*mY(i,By,k)-4.*mY(i,By+1,k)+mY(i,By+2,k));
	  zr = -0.5*(3.*mZ(i,By,k)-4.*mZ(i,By+1,k)+mZ(i,By+2,k));
	} else if (j==Ny){
	  xr=0.5*(3.*mX(i,Ny,k)-4.*mX(i,Ny-1,k)+mX(i,Ny-2,k));
	  yr=0.5*(3.*mY(i,Ny,k)-4.*mY(i,Ny-1,k)+mY(i,Ny-2,k));
	  zr=0.5*(3.*mZ(i,Ny,k)-4.*mZ(i,Ny-1,k)+mZ(i,Ny-2,k));
	} else {
	  xr=0.5*(mX(i,j+1,k)-mX(i,j-1,k));
	  yr=0.5*(mY(i,j+1,k)-mY(i,j-1,k));
	  zr=0.5*(mZ(i,j+1,k)-mZ(i,j-1,k));
	}
	if (k==Bz){
	  xs = -0.5*(3.*mX(i,j,Bz)-4.*mX(i,j,Bz+1)+mX(i,j,Bz+2));
	  ys = -0.5*(3.*mY(i,j,Bz)-4.*mY(i,j,Bz+1)+mY(i,j,Bz+2));
	  zs = -0.5*(3.*mZ(i,j,Bz)-4.*mZ(i,j,Bz+1)+mZ(i,j,Bz+2));
	} else if (k==Nz){
	  xs=0.5*(3.*mX(i,j,Nz)-4.*mX(i,j,Nz-1)+mX(i,j,Nz-2));
	  ys=0.5*(3.*mY(i,j,Nz)-4.*mY(i,j,Nz-1)+mY(i,j,Nz-2));
	  zs=0.5*(3.*mZ(i,j,Nz)-4.*mZ(i,j,Nz-1)+mZ(i,j,Nz-2));
	} else {
	  xs=0.5*(mX(i,j,k+1)-mX(i,j,k-1));
	  ys=0.5*(mY(i,j,k+1)-mY(i,j,k-1));
	  zs=0.5*(mZ(i,j,k+1)-mZ(i,j,k-1));
	}
	mJ(i,j,k) = xq*yr*zs+yq*zr*xs+zq*xr*ys-xq*zr*ys-yq*xr*zs-zq*yr*xs;
	mQ(1,i,j,k)= (1.0/mJ(i,j,k))*(yr*zs-zr*ys);
	mQ(2,i,j,k)= (1.0/mJ(i,j,k))*(zr*xs-xr*zs);
	mQ(3,i,j,k)= (1.0/mJ(i,j,k))*(xr*ys-yr*xs);
	mR(1,i,j,k)= (1.0/mJ(i,j,k))*(zq*ys-yq*zs);
	mR(2,i,j,k)= (1.0/mJ(i,j,k))*(xq*zs-zq*xs);
	mR(3,i,j,k)= (1.0/mJ(i,j,k))*(yq*xs-ys*xq);
	mS(1,i,j,k)= (1.0/mJ(i,j,k))*(yq*zr-zq*yr);
	mS(2,i,j,k)= (1.0/mJ(i,j,k))*(zq*xr-xq*zr);
	mS(3,i,j,k)= (1.0/mJ(i,j,k))*(xq*yr-xr*yq);
//
	if (mJ(i,j,k) < minJ) minJ=mJ(i,j,k);
	if (mJ(i,j,k) > maxJ) maxJ=mJ(i,j,k);
      } // end for i,j,k  
  double minJglobal, maxJglobal;
  MPI_Allreduce( &minJ, &minJglobal, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
  MPI_Allreduce( &maxJ, &maxJglobal, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator);
  if (mVerbose>3 && proc_zero())
    printf("*** Jacobian of metric: minJ = %e maxJ = %e\n", minJglobal, maxJglobal);

// just save the results for now... do the sanity check later
  m_minJacobian= minJglobal;
  m_maxJacobian= maxJglobal;
}

//-----------------------------------------------------------------------
void EW::generate_grid()
{
   // Generate grid on domain: topography <= z <= zmax, 
   // The 2D grid on z=zmax, is given by ifirst <= i <= ilast, jfirst <= j <= jlast
   // spacing h. 

  if (!m_topography_exists ) return;
  
//  m_grid_interpolation_order = a_order;

  if (mVerbose >= 1 && proc_zero())
    cout << "***inside generate_grid***"<< endl;

// get the size from the top Cartesian grid
  int g = mNumberOfCartesianGrids-1;
  int ifirst = m_iStart[g];
  int ilast  = m_iEnd[g];
  int jfirst = m_jStart[g];
  int jlast  = m_jEnd[g];

  double h = mGridSize[g]; // grid size must agree with top cartesian grid
  double zMaxCart = m_zmin[g]; // bottom z-level for curvilinear grid

  int i, j;
  int gTop = mNumberOfGrids-1;
  int Nz = m_kEnd[gTop] - m_ghost_points;

  if(mVerbose > 4 &&  proc_zero() )
  {
    printf("generate_grid: Number of grid points in curvilinear grid = %i, kStart = %i, kEnd = %i\n", 
	Nz, m_kStart[gTop], m_kEnd[gTop]);
  }

// generate the grid by calling the curvilinear mapping function
  double X0, Y0, Z0;
  int k;
  for (k=m_kStart[gTop]; k<=m_kEnd[gTop]; k++)
    for (j=m_jStart[gTop]; j<=m_jEnd[gTop]; j++)
      for (i=m_iStart[gTop]; i<=m_iEnd[gTop]; i++)
      {
	curvilinear_grid_mapping((double) i, (double) j, (double) k, X0, Y0, Z0);
	mX(i,j,k) = X0;
	mY(i,j,k) = Y0;
	mZ(i,j,k) = Z0;
      }

// tmp
// test the inverse mapping
//  double q0, r0, s0, dist=0.;
//  for (k=m_kStart[gTop]; k<=m_kEnd[gTop]; k++)
//    for (j=m_jStart[gTop]; j<=m_jEnd[gTop]; j++)
//      for (i=m_iStart[gTop]; i<=m_iEnd[gTop]; i++)
//      {
// 	invert_curvilinear_grid_mapping(mX(i,j,k), mY(i,j,k), mZ(i,j,k), q0, r0, s0);
// 	dist += SQR(q0 - (double) i) + SQR(r0 - (double) j) + SQR(s0 - (double) k);
//      }
//   double totalDist;
//   MPI_Allreduce( &dist, &totalDist, 1, MPI_DOUBLE, MPI_SUM, m_cartesian_communicator );
//   if (m_myRank == 0)
//     printf("L2-error in inverse of curvilinear mapping = %e\n", sqrt(totalDist));
// end test

// make sure all processors have made their grid before we continue
   MPI_Barrier(m_cartesian_communicator);

// Smooth the grid (only the Z component for now)
// NOTE: the current smoothing algorithm makes the error larger rather than smaller!
  int maxIter=0; // number of iterations
  double rf=0.05; // rf<1/6 for stability
  if (mVerbose >= 1 && proc_zero() && maxIter>0)
    cout << "***smoothing the grid with " << maxIter << " Jacobi iterations and relaxation factor " << rf << " ***"<< endl;

  int topLevel = mNumberOfGrids-1;
  int iter;
  
// temporary storage: How can I use mJ for temporary storage?
  Sarray tmp;
  tmp.define(m_iStart[topLevel],m_iEnd[topLevel],m_jStart[topLevel],m_jEnd[topLevel],m_kStart[topLevel],m_kEnd[topLevel]);

// initialize to make the Dirichlet boundary conditions work
    for (k = m_kStart[topLevel]; k <= m_kEnd[topLevel]; k++)
      for (j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
	for (i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
	{
	  tmp(i,j,k) = mZ(i,j,k);
	}

// Laplacian filter
  for (iter=0; iter < maxIter; iter++)
  {
// loop over all interior points
    for (k = m_kStart[topLevel]+m_ghost_points+1; k <= m_kEnd[topLevel]-m_ghost_points-2; k++)
      for (j = m_jStart[topLevel]+1; j <= m_jEnd[topLevel]-1; j++)
	for (i = m_iStart[topLevel]+1; i <= m_iEnd[topLevel]-1; i++)
	{
	  tmp(i,j,k) = mZ(i,j,k) + rf*(mZ(i+1,j,k) + mZ(i-1,j,k) + mZ(i,j+1,k) + mZ(i,j-1,k) + mZ(i,j,k+1) + mZ(i,j,k-1) - 6.*mZ(i,j,k));
	}

// impose Neumann bc on the i and j sides
    for (k = m_kStart[topLevel]+m_ghost_points+1; k <= m_kEnd[topLevel]-m_ghost_points-2; k++)
    {
      for (j = m_jStart[topLevel]+1; j <= m_jEnd[topLevel]-1; ++j)
      {
	i = m_iStart[topLevel];
	tmp(i,j,k) = tmp(i+1,j,k);
	i = m_iEnd[topLevel];
	tmp(i,j,k) = tmp(i-1,j,k);
      }

      for (i = m_iStart[topLevel]+1; i <= m_iEnd[topLevel]-1; ++i)
      {
	j = m_jStart[topLevel];
	tmp(i,j,k) = tmp(i,j+1,k);
	j = m_jEnd[topLevel];
	tmp(i,j,k) = tmp(i,j-1,k);
      }
// Corners
      i = m_iStart[topLevel];
      j = m_jStart[topLevel];
      tmp(i,j,k) = tmp(i+1,j+1,k);

      i = m_iEnd[topLevel];
      j = m_jStart[topLevel];
      tmp(i,j,k) = tmp(i-1,j+1,k);

      i = m_iStart[topLevel];
      j = m_jEnd[topLevel];
      tmp(i,j,k) = tmp(i+1,j-1,k);
    
      i = m_iEnd[topLevel];
      j = m_jEnd[topLevel];
      tmp(i,j,k) = tmp(i-1,j-1,k);
    } // end Neumann loop

// communicate parallel ghost points
    communicate_array( tmp, topLevel );

// update solution (Dirichlet are imposed implicitly by never changing the tmp array along the top or bottom boundary)
    for (k = m_kStart[topLevel]; k <= m_kEnd[topLevel]; k++)
      for (j = m_jStart[topLevel]; j <= m_jEnd[topLevel]; ++j)
	for (i = m_iStart[topLevel]; i <= m_iEnd[topLevel]; ++i)
	{
	  mZ(i,j,k) = tmp(i,j,k);
	}

  }// end for iter (grid smoother)
  
// tmp
// calculate min and max((mZ(i,j,k)-mZ(i,j,k-1))/h) for k=Nz
  k = Nz;
  double hRatio;
// tmp
  double mZmin, mZmax;
  mZmin = 1.e9;
  mZmax = 0.;
  for (j=m_jStart[gTop]; j<=m_jEnd[gTop]; j++)
    for (i=m_iStart[gTop]; i<=m_iEnd[gTop]; i++)
    {
      hRatio = (mZ(i,j,k)-mZ(i,j,k-1))/mGridSize[gTop];
	if (hRatio < mZmin) mZmin = hRatio;
	if (hRatio > mZmax) mZmax = hRatio;
    }
  double zMinGlobal, zMaxGlobal;
  MPI_Allreduce( &mZmin, &zMinGlobal, 1, MPI_DOUBLE, MPI_MIN, m_cartesian_communicator);
  MPI_Allreduce( &mZmax, &zMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator);
  if(mVerbose > 3 &&  proc_zero() )
  {
    printf("Curvilinear/Cartesian interface (k=Nz-1): Min grid size ratio - 1 = %e, max ratio z - 1 = %e, top grid # = %i\n", 
           zMinGlobal-1., zMaxGlobal-1., gTop);
  }
// end tmp

}

//-----------------------------------------------------------------------
bool EW::curvilinear_grid_mapping( double q, double r, double s, double & X0, double & Y0, double & Z0 )
{
// if (q,r) is on this processor (need a 2x2 interval in (i,j)-index space:
// Return true and assign (X0,Y0,Z0) corresponding to (q,r,s)

// Returns false if 
// 1) (q,r,s) is outside the global parameter domain (expanded by ghost points)
// 2) There is no curvilinear grid. Still computes (X0, Y0) based on the top Cartesian grid size
// 3) (q,r) is not on this processor. Still computes (X0, Y0)

// NOTE:
// The parameters are normalized such that 1 <= q <= Nx is the full domain (without ghost points),
//  1 <= r <= Ny, 1 <= s <= Nz.

  int gCurv = mNumberOfGrids - 1;
  double h = mGridSize[gCurv];
// check global parameter space
  double qMin = (double) (1- m_ghost_points);
  double qMax = (double) (m_global_nx[gCurv] + m_ghost_points);
  double rMin = (double) (1- m_ghost_points);
  double rMax = (double) (m_global_ny[gCurv] + m_ghost_points);
  double sMin = (double) m_kStart[gCurv];
  double sMax = (double) m_kEnd[gCurv];

  if (! (q >= qMin && q <= qMax && r >= rMin && r <= rMax && s >= sMin && s <= sMax))
  {
    cout << "curvilinear_grid_mapping: input parameters out of bounds (q,r,s) = " << q << ", " << r << ", " << s << endl;
    return false;
  }
  
  X0 = (q-1.0)*h;
  Y0 = (r-1.0)*h;

  if (!topographyExists())
    return false;

// bottom z-level for curvilinear grid = top z-level for highest Cartesian grid
  double zMaxCart = m_zmin[mNumberOfCartesianGrids-1]; 
  
// ************************
// compute index interval based on (q,r)
  int iNear, jNear, kNear, g, i, j, k;

  Z0 = zMaxCart - h; // to make computeNearestGridPoint think we are in the curvilinear grid
  computeNearestGridPoint(iNear, jNear, kNear, g, X0, Y0, Z0);

  if (g != gCurv)
    return false;

  double tau; // holds the elevation at (q,r). Recall that elevation=-z
  if (m_analytical_topo)
  {
    tau = m_GaussianAmp*exp(-SQR((X0-m_GaussianXc)/m_GaussianLx) 
                            -SQR((Y0-m_GaussianYc)/m_GaussianLy)); 
  }
  else // general case: interpolate mTopoGrid array
  {
// if (X0, Y0) falls within roundoff of grid point (iNear,jNear), we only need that grid point on this proc,
// otherwise we need the 2x2 area [i,i+1] by [j,j+1]

    double xPt = (iNear-1)*h;
    double yPt = (jNear-1)*h;

// first check if we are very close to a grid point
    bool smackOnTop = (fabs((xPt-X0)/h) < 1.e-9 && fabs((yPt-Y0)/h) < 1.e-9);

    if (smackOnTop)
    {
      if (!point_in_proc(iNear,jNear,gCurv))
        return false;
      tau = mTopoGrid(iNear,jNear,1);
    }
    else
    {
      computeNearestLowGridPoint(i, j, k, g, X0, Y0, Z0);
// There are some subtle issues with the bi-cubic interpolation near parallel processor boundaries, 
// see invert_curvilinear_mapping (below) for a discussion

// bi-cubic interpolation for O(h^4) accuracy
//       if ((point_in_proc(i-1,j-1,gCurv) && point_in_proc(i,j-1,gCurv) && point_in_proc(i+1,j-1,gCurv) && point_in_proc(i+2,j-1,gCurv) &&
// 		point_in_proc(i-1,j,gCurv) && point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv) && point_in_proc(i+2,j,gCurv) &&
// 		point_in_proc(i-1,j+1,gCurv) && point_in_proc(i,j+1,gCurv) && point_in_proc(i+1,j+1,gCurv) && point_in_proc(i+2,j+1,gCurv) &&
// 		point_in_proc(i-1,j+2,gCurv) && point_in_proc(i,j+2,gCurv) && point_in_proc(i+1,j+2,gCurv) && point_in_proc(i+2,j+2,gCurv) ) )
//       {
// 	double Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1, tjp2;
// 	Qim1 = (q-i)*(q-i-1)*(q-i-2)/(-6.);
// 	Qi   = (q-i+1)*(q-i-1)*(q-i-2)/(2.);
// 	Qip1 = (q-i+1)*(q-i)*(q-i-2)/(-2.);
// 	Qip2 = (q-i+1)*(q-i)*(q-i-1)/(6.);

// 	Rjm1 = (r-j)*(r-j-1)*(r-j-2)/(-6.);
// 	Rj   = (r-j+1)*(r-j-1)*(r-j-2)/(2.);
// 	Rjp1 = (r-j+1)*(r-j)*(r-j-2)/(-2.);
// 	Rjp2 = (r-j+1)*(r-j)*(r-j-1)/(6.);

// 	tjm1 = Qim1*mTopoGrid(i-1,j-1,1) + Qi*mTopoGrid(i,j-1,1) +  Qip1*mTopoGrid(i+1,j-1,1) +  Qip2*mTopoGrid(i+2,j-1,1);
// 	tj   = Qim1*mTopoGrid(i-1,j,1) + Qi*mTopoGrid(i,j,1) +  Qip1*mTopoGrid(i+1,j,1) +  Qip2*mTopoGrid(i+2,j,1);
// 	tjp1 = Qim1*mTopoGrid(i-1,j+1,1) + Qi*mTopoGrid(i,j+1,1) +  Qip1*mTopoGrid(i+1,j+1,1) +  Qip2*mTopoGrid(i+2,j+1,1);
// 	tjp2 = Qim1*mTopoGrid(i-1,j+2,1) + Qi*mTopoGrid(i,j+2,1) +  Qip1*mTopoGrid(i+1,j+2,1) +  Qip2*mTopoGrid(i+2,j+2,1);

// 	tau = Rjm1*tjm1 + Rj*tj + Rjp1*tjp1 + Rjp2*tjp2;
//       }
//      else if ( ( point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv) && point_in_proc(i,j+1,gCurv) && 
//	     point_in_proc(i+1,j+1,gCurv) ) )
      if ( ( point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv) && point_in_proc(i,j+1,gCurv) && 
	     point_in_proc(i+1,j+1,gCurv) ) )
      {
// linear interpolation to define topography between grid points
	double Qi, Qip1, Rj, Rjp1;
	Qi = (i+1 - q);
	Qip1 = (q - i);
	Rj = (j+1 - r);
	Rjp1 = (r - j);
	tau = mTopoGrid(i,j,1)*Rj*Qi + mTopoGrid(i,j+1,1)*Rjp1*Qi + mTopoGrid(i+1,j,1)*Rj*Qip1 + 
	  mTopoGrid(i+1,j+1,1)*Rjp1*Qip1;  
      }
      else
      {
        return false;
      }
    } // end else...(not smackOnTop)
    
  }// end general case: interpolating mTopoGrid array  

// now we need to calculate Z0 = Z(q,r,s)

// setup parameters for grid mapping
  int Nz = m_kEnd[gCurv] - m_ghost_points;
  double zetaBreak = 0.95; // zeta  > zetaBreak gives constant grid size = h
  double sBreak = 1. + zetaBreak*(Nz-1);

  double zeta, c1=0., c2=0., c3=0., zMax;
  zMax = zMaxCart - (Nz - sBreak)*h;

// quadratic term to make variation in grid size small at bottom boundary
  c1 = zMax  + tau - mGridSize[gCurv]*(sBreak-1);
// cubic term to make 2nd derivative zero at zeta=1
  if (m_grid_interpolation_order>=3) 
    c2 = c1;
  else
    c2 = 0.;
	
// 4th order term takes care of 3rd derivative, but can make grid warp itself inside out
  if (m_grid_interpolation_order>=4) 
    c3 = c2;
  else
    c3 = 0.;

// the forward mapping is ...
  if (s <= (double) sBreak)
  {
    zeta = (s-1)/(sBreak-1.);
    Z0 = (1.-zeta)*(-tau) 
      + zeta*(zMax + c1*(1.-zeta) + c2*SQR(1.-zeta) + c3*(1.-zeta)*SQR(1.-zeta));
  }
  else
  {
    Z0 = zMax + (s-sBreak)*h;
  }
  
  return true;
}

//-----------------------------------------------------------------------
bool EW::invert_curvilinear_grid_mapping( double X0, double Y0, double Z0, double& q, double& r, double& s )
{
// If (X0, Y0, Z0) is on the curvilinear grid and (X0, Y0) is on this processor:
// Return true and assigns (q,r,s) corresponding to point (mX0,mY0,mZ0)

// Returns false if 
// 1) There is no curvilinear grid. Computes (q,r) based on the top Cartesian grid size
// 2) Z0 > m_zmin[topCartGrid]. Still computes (q,r) as if Z0 would be on the curvilinear grid
// 3) (X0, Y0) is not on this processor. Still computes (q,r)

// NOTE:
// Normalize parameters such that 1 <= q <= Nx is the full domain (without ghost points),
//  1 <= r <= Ny, 1 <= s <= Nz.

  int gCurv = mNumberOfGrids - 1;
  double h = mGridSize[gCurv];
  q = X0/h + 1.0;
  r = Y0/h + 1.0;
  s = 0.;

  if (!topographyExists())
    return false;

// compute index interval based on (q,r)
  int iNear, jNear, kNear, g, i, j, k;

// for ghost points in the curvilinear grid, calling computeNearestGridPoint with Z0 will pick up the top Cartesian grid

// bottom z-level for curvilinear grid = top z-level for highest Cartesian grid
  double zMaxCart = m_zmin[mNumberOfCartesianGrids-1]; 

  computeNearestGridPoint(iNear, jNear, kNear, g, X0, Y0, zMaxCart - h);

  if (g != gCurv)
  {
    cout << "invert_curvilinear_grid_mapping: computeNearestGridPoint returned g!=gCurv for (X0, Y0, Z0) = " 
	 << X0 << ", " << Y0 << ", " << Z0 << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  
  double tau; // holds the elevation at (q,r). Recall that elevation=-z

  if (m_analytical_topo)
  {
    tau = m_GaussianAmp*exp(-SQR((X0-m_GaussianXc)/m_GaussianLx) 
                            -SQR((Y0-m_GaussianYc)/m_GaussianLy)); 
  }
  else // general case: interpolate mTopoGrid array
  {
// if (X0, Y0) falls within roundoff of grid point (iNear,jNear), we only need that grid point on this proc,
// otherwise we need the 2x2 area [i,i+1] by [j,j+1]

    double xPt = (iNear-1)*h;
    double yPt = (jNear-1)*h;

// first check if we are very close to a grid point
    bool smackOnTop = (fabs((xPt-X0)/h) < 1.e-9 && fabs((yPt-Y0)/h) < 1.e-9);

    if (smackOnTop)
    {
      if (!point_in_proc(iNear,jNear,gCurv))
        return false;
      tau = mTopoGrid(iNear,jNear,1);
    }
    else
    {
      computeNearestLowGridPoint(i, j, k, g, X0, Y0, Z0);

      double Qim1, Qi, Qip1, Qip2, Rjm1, Rj, Rjp1, Rjp2, tjm1, tj, tjp1, tjp2;

// The bi-cubic interpolation has the following problem:
// When the source is close to a processor boundary, the same source must be discretized on two (or more processors).
// However, the topography grid only has 2 points overlap, which sometimes
// is insufficient for the bi-cubic interpolation formula. 
// In practice, the difference between the bi-linear and bi-cubic formula makes very little difference, so it
// is more robust to always use the bi-linear formula.

// // bi-cubic interpolation for O(h^4) accuracy
//       if ( (point_in_proc(i-1,j-1,gCurv) && point_in_proc(i,j-1,gCurv) && point_in_proc(i+1,j-1,gCurv) && point_in_proc(i+2,j-1,gCurv) &&
// 	    point_in_proc(i-1,j,gCurv) && point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv) && point_in_proc(i+2,j,gCurv) &&
// 	    point_in_proc(i-1,j+1,gCurv) && point_in_proc(i,j+1,gCurv) && point_in_proc(i+1,j+1,gCurv) && point_in_proc(i+2,j+1,gCurv) &&
// 	    point_in_proc(i-1,j+2,gCurv) && point_in_proc(i,j+2,gCurv) && point_in_proc(i+1,j+2,gCurv) && point_in_proc(i+2,j+2,gCurv) ) )
//       {
// 	Qim1 = (q-i)*(q-i-1)*(q-i-2)/(-6.);
// 	Qi   = (q-i+1)*(q-i-1)*(q-i-2)/(2.);
// 	Qip1 = (q-i+1)*(q-i)*(q-i-2)/(-2.);
// 	Qip2 = (q-i+1)*(q-i)*(q-i-1)/(6.);

// 	Rjm1 = (r-j)*(r-j-1)*(r-j-2)/(-6.);
// 	Rj   = (r-j+1)*(r-j-1)*(r-j-2)/(2.);
// 	Rjp1 = (r-j+1)*(r-j)*(r-j-2)/(-2.);
// 	Rjp2 = (r-j+1)*(r-j)*(r-j-1)/(6.);

// 	tjm1 = Qim1*mTopoGrid(i-1,j-1,1) + Qi*mTopoGrid(i,j-1,1) +  Qip1*mTopoGrid(i+1,j-1,1) +  Qip2*mTopoGrid(i+2,j-1,1);
// 	tj   = Qim1*mTopoGrid(i-1,j,1) + Qi*mTopoGrid(i,j,1) +  Qip1*mTopoGrid(i+1,j,1) +  Qip2*mTopoGrid(i+2,j,1);
// 	tjp1 = Qim1*mTopoGrid(i-1,j+1,1) + Qi*mTopoGrid(i,j+1,1) +  Qip1*mTopoGrid(i+1,j+1,1) +  Qip2*mTopoGrid(i+2,j+1,1);
// 	tjp2 = Qim1*mTopoGrid(i-1,j+2,1) + Qi*mTopoGrid(i,j+2,1) +  Qip1*mTopoGrid(i+1,j+2,1) +  Qip2*mTopoGrid(i+2,j+2,1);

// 	tau = Rjm1*tjm1 + Rj*tj + Rjp1*tjp1 + Rjp2*tjp2;
// // tmp
// 	printf("invert_curvilinear_mapping: q=%e, r=%e, Cubic tau=%e\n", q, r, tau);
//       }
//      else if ( ( point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv) && point_in_proc(i,j+1,gCurv) && 
//		  point_in_proc(i+1,j+1,gCurv) ) )
      if ( ( point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv) && point_in_proc(i,j+1,gCurv) && 
	     point_in_proc(i+1,j+1,gCurv) ) )
      {	
// use linear interpolation if there are not enough points for the bi-cubic formula
	Qi = (i+1 - q);
	Qip1 = (q - i);
	Rj = (j+1 - r);
	Rjp1 = (r - j);
	tau = mTopoGrid(i,j,1)*Rj*Qi + mTopoGrid(i,j+1,1)*Rjp1*Qi + mTopoGrid(i+1,j,1)*Rj*Qip1 + 
	  mTopoGrid(i+1,j+1,1)*Rjp1*Qip1;  
// tmp
//	printf("invert_curvilinear_mapping: q=%e, r=%e, Linear tau=%e\n", q, r, tau);
      }
      else
      {
	return false;
      }
      
    } // end else... (not smackOnTop)
    
  } // end else... general case
  
// now we need to calculate s: Z(q,r,s) = Z0

// setup parameters for grid mapping (same as curvilinear_grid_mapping)
  int Nz = m_kEnd[gCurv] - m_ghost_points;
  double zetaBreak = 0.95; // zeta  > zetaBreak gives constant grid size = h
  double sBreak = 1. + zetaBreak*(Nz-1);

  double zeta, c1=0., c2=0., c3=0., zMax;
  zMax = zMaxCart - (Nz - sBreak)*h;

// quadratic term to make variation in grid size small at bottom boundary
  c1 = zMax  + tau - mGridSize[gCurv]*(sBreak-1);
// cubic term to make 2nd derivative zero at zeta=1
  if (m_grid_interpolation_order>=3) 
    c2 = c1;
  else
    c2 = 0.;
	
// 4th order term takes care of 3rd derivative, but can make grid warp itself inside out
  if (m_grid_interpolation_order>=4) 
    c3 = c2;
  else
    c3 = 0.;

// the forward mapping is ...
//  if (s <= (double) sBreak)
//  {
//    zeta = (s-1)/(sBreak-1.);
//    Z0 = (1.-zeta)*(-tau) 
//      + zeta*(zMax + c1*(1.-zeta) + c2*SQR(1.-zeta) + c3*(1.-zeta)*SQR(1.-zeta));
//  }
//  else
//  {
//    Z0 = zMax + (s-sBreak)*h;
//  }

    if (Z0 >= zMax)
    {
      s = (Z0 - zMax)/h + sBreak;
    }
    else
    {
// Get initial guess by taking c1=c2=c3=0
      zeta = (Z0 + tau)/(zMax + tau);
      double F0, Fp0, zeta0=zeta;
    
// Do a couple of Newton iterations
      int maxIter=5;
//      cout << "invert_curvilinear_grid_mapping: initial guess zeta0= " << zeta0 << endl;
      for (int iter=0; iter<maxIter; iter++)
      {
        F0  = (1.-zeta0)*(-tau) + zeta0*(zMax + c1*(1.-zeta0) + c2*SQR(1.-zeta0) + c3*(1.-zeta0)*SQR(1.-zeta0)) - Z0;
        Fp0 = tau 
             + zMax + c1*(1.-zeta0) + c2*SQR(1.-zeta0) + c3*(1.-zeta0)*SQR(1.-zeta0)
             + zeta0*(-c1 - 2.*c2*(1.-zeta0) - 3.*c3*SQR(1.-zeta0));
// tmp
//        cout << "invert_curvilinear_grid_mapping: iter= " << iter << " zeta= " << zeta0 <<  " F0= " << F0 << " Fp0= " << Fp0 << endl;
        zeta0 -= F0/Fp0;
      }
// check convergence
      if (fabs(F0) > 1.e-7)
      {
        cout << "invert_curvilinear_grid_mapping: poor convergence for X0, Y0, Z0 = " << X0 << ", " << Y0 << ", " << Z0 << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      zeta = zeta0;
// tmp
//    cout << "invert_curvilinear_grid_mapping: final zeta= " << zeta << endl;
      s = 1. + (sBreak-1.)*zeta;
    } // end Z0 < zMax

// tmp
//    printf("invert_curvilinear_mapping: X0=%e, Y0=%e, Z0=%e, q=%e, r=%e, s=%e\n", X0, Y0, Z0, q, r, s);

// if we got this far, the inversion was successful
    return true;
}

//-----------------------------------------------------------------------
bool EW::find_curvilinear_derivatives_at_point( double q, double r, double s,
						double qX[], double rX[], double sX[] )
{
   // return metric derivatives at parameter value (q,r,s), note 1 <= q <= Nq, 1 <= r <= Nr, 1 <= s <= Ns
   // Input: (q,r,s)
   // Output: (qx0,..,sz0)

// NOTE:
// The parameters are normalized such that 1 <= q <= Nx is the full domain (without ghost points),
//  1 <= r <= Ny, 1 <= s <= Nz.

  int gCurv = mNumberOfGrids - 1;
  double h = mGridSize[gCurv];
// check global parameter space
  double qMin = (double) (1- m_ghost_points);
  double qMax = (double) (m_global_nx[gCurv] + m_ghost_points);
  double rMin = (double) (1- m_ghost_points);
  double rMax = (double) (m_global_ny[gCurv] + m_ghost_points);
  double sMin = (double) m_kStart[gCurv];
  double sMax = (double) m_kEnd[gCurv];

  if (! (q >= qMin && q <= qMax && r >= rMin && r <= rMax && s >= sMin && s <= sMax))
  {
    cout << "find_curvilinear_derivatives_at_point: input parameters out of bounds (q,r,s) = " << q << ", " << r << ", " << s << endl;
    return false;
  }

  double X0 = (q-1.0)*h;
  double Y0 = (r-1.0)*h;

  if (!topographyExists())
    return false;

// bottom z-level for curvilinear grid = top z-level for highest Cartesian grid
  double zMaxCart = m_zmin[mNumberOfCartesianGrids-1]; 
  
// ************************
// compute index interval based on (q,r)
  int iNear, jNear, kNear, g, i, j, k;

  double Z0 = zMaxCart - h; // to make computeNearestGridPoint think we are in the curvilinear grid
  computeNearestGridPoint(iNear, jNear, kNear, g, X0, Y0, Z0);

  double xPt = (iNear-1)*h;
  double yPt = (jNear-1)*h;

// first check if we are very close to a grid point
  bool smackOnTop = (fabs((xPt-X0)/h) < 1.e-9 && fabs((yPt-Y0)/h) < 1.e-9);

  k = (int) floor(s);
  
  double Wi, Wip1, Wj, Wjp1, Wk, Wkp1;
  Wk = (k+1-s);
  Wkp1 = (s-k);
  
  if (smackOnTop)
  {
    if (!point_in_proc(iNear,jNear,gCurv))
      return false;
    for (int p=1; p<=3; p++)
    {
      qX[p-1] = Wk*mQ(p, iNear, jNear, k) + Wkp1*mQ(p, iNear, jNear, k+1);
      rX[p-1] = Wk*mR(p, iNear, jNear, k) + Wkp1*mR(p, iNear, jNear, k+1);
      sX[p-1] = Wk*mS(p, iNear, jNear, k) + Wkp1*mS(p, iNear, jNear, k+1);
    }
  }
  else
  {
    computeNearestLowGridPoint(i, j, k, g, X0, Y0, Z0);
    if ( !( point_in_proc(i,j,gCurv) && point_in_proc(i+1,j,gCurv) && point_in_proc(i,j+1,gCurv) && 
            point_in_proc(i+1,j+1,gCurv) ) )
      return false;

// linear interpolation to define topography between grid points
    Wi = (i+1 - q);
    Wip1 = (q - i);
    Wj = (j+1 - r);
    Wjp1 = (r - j);

    for (int p=1; p<=3; p++)
    {
      qX[p-1] = 
	Wj*Wi*    (Wk*mQ(p, i, j, k)     + Wkp1*mQ(p, i, j, k+1)) + 
	Wjp1*Wi*  (Wk*mQ(p, i, j+1, k)   + Wkp1*mQ(p, i, j+1, k+1)) + 
	Wj*Wip1*  (Wk*mQ(p, i+1, j, k)   + Wkp1*mQ(p, i+1, j, k+1)) + 
	Wjp1*Wip1*(Wk*mQ(p, i+1, j+1, k) + Wkp1*mQ(p, i+1, j+1, k+1));
      rX[p-1] = 
	Wj*Wi*    (Wk*mR(p, i, j, k)     + Wkp1*mR(p, i, j, k+1)) + 
	Wjp1*Wi*  (Wk*mR(p, i, j+1, k)   + Wkp1*mR(p, i, j+1, k+1)) + 
	Wj*Wip1*  (Wk*mR(p, i+1, j, k)   + Wkp1*mR(p, i+1, j, k+1)) + 
	Wjp1*Wip1*(Wk*mR(p, i+1, j+1, k) + Wkp1*mR(p, i+1, j+1, k+1));
      sX[p-1] = 
	Wj*Wi*    (Wk*mS(p, i, j, k)     + Wkp1*mS(p, i, j, k+1)) + 
	Wjp1*Wi*  (Wk*mS(p, i, j+1, k)   + Wkp1*mS(p, i, j+1, k+1)) + 
	Wj*Wip1*  (Wk*mS(p, i+1, j, k)   + Wkp1*mS(p, i+1, j, k+1)) + 
	Wjp1*Wip1*(Wk*mS(p, i+1, j+1, k) + Wkp1*mS(p, i+1, j+1, k+1));
    }
  } // !smackOnTop
  return true;
}
