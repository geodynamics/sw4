#include "EW.h"
#include "cf_interface.h"

#define SQR(x) ((x)*(x))

extern "C" 
{
   void twfrsurfz_wind( int *ifirst, int *ilast, int *jfirst, int *jlast, int *kfirst, int *klast,
                        float_sw4 *h, int* kz, float_sw4 *t, float_sw4 *omega, float_sw4 *c, float_sw4 *phase, float_sw4 *bforce,
                        float_sw4 *mu, float_sw4 *lambda, float_sw4 *zmin,
                        int *i1, int *i2, int *j1, int *j2 );

   void twfrsurfz_att_wind( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                          float_sw4 h, int kz, float_sw4 t, float_sw4 omega, float_sw4 c, float_sw4 phase,
                          float_sw4 *bforce, float_sw4 *mu, float_sw4 *lambda, float_sw4 zmin,
                          int i1, int i2, int j1, int j2 );
   void twfrsurfzsg_wind( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                          float_sw4 h, int kz, float_sw4 t, float_sw4 omega, float_sw4 c, float_sw4 phase, float_sw4 omstrx, float_sw4 omstry,
                          float_sw4 *bforce, float_sw4 *mu, float_sw4 *lambda, float_sw4 zmin,
                          int i1, int i2, int j1, int j2 );

   void twfrsurfzsg_att_wind( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
                          float_sw4 h, int kz, float_sw4 t, float_sw4 omega, float_sw4 c, float_sw4 phase, float_sw4 omstrx, float_sw4 omstry,
                          float_sw4 *bforce, float_sw4 *mu, float_sw4 *lambda, float_sw4 zmin,
                          int i1, int i2, int j1, int j2 );
}

//-----------------------------------------------------------------------
void EW::consintp( Sarray& Uf, Sarray& Unextf, Sarray& Bf, Sarray& Muf, Sarray& Lambdaf, Sarray& Rhof, float_sw4 hf,
		   Sarray& Uc, Sarray& Unextc, Sarray& Bc, Sarray& Muc, Sarray& Lambdac, Sarray& Rhoc, float_sw4 hc,
		   float_sw4 cof, int gc, int gf, int is_periodic[2])
{
   // At boundaries to the left and right, at least three ghost points are required
   // e.g., domain in i-direction:   i=-2,-1,0,1,2,...,Ni,Ni+1,Ni+2,Ni+3
   // we solve for ghost points at i=2,..,Ni-1, assuming Dirichlet conditions given on i=1,i=Ni and
   // at the ghost points i=-2,-1,0,Ni+1,Ni+2,Ni+3. 
   //
   // The arrays Bf,Unextf, etc are assumed to be computed with all z-ghost points zero,
   // i.e., Uf(i,Nkf+1) was zero for i=-2,-1,..,Ni+3 when Bf and Unextf were evaluated from Uf 
   // (similarly for Unextc and Bc).
   //
   // Before this routine is called, correct boundary values for
   // Uf(-2,Nkf+1),..Uf(1,Nkf+1) (and similarly at the upper i-boundary, and for Uc) must be imposed.
   // In this way the restriction and prolongation stencils can be computed without any special
   // treatment near the (i,j)-boundaries. 

   float_sw4 *a_strc_x = m_sg_str_x[gc];
   float_sw4 *a_strc_y = m_sg_str_y[gc];
   float_sw4 *a_strf_x  = m_sg_str_x[gf];
   float_sw4 *a_strf_y  = m_sg_str_y[gf];

// stretching on the coarse side
#define strc_x(i) a_strc_x[(i-m_iStart[gc])]   
#define strc_y(j) a_strc_y[(j-m_jStart[gc])]   

// stretching on the fine side
#define strf_x(i) a_strf_x[(i-m_iStart[gf])]   
#define strf_y(j) a_strf_y[(j-m_jStart[gf])]   

      
   const float_sw4 i16 = 1.0/16;
   const float_sw4 i256 = 1.0/256;
   const float_sw4 i1024 = 1.0/1024;

   int jcb, jce, icb, ice, jfb, jfe, ifb, ife, nkf;
   float_sw4 nuf = mDt*mDt/(cof*hf*hf); // cof=12 for the predictor, cof=1 for the corrector (argument to this routine)
   float_sw4 nuc = mDt*mDt/(cof*hc*hc);
   float_sw4 ihc = 1/hc, ihf=1/hf;
   float_sw4 jacerr = m_citol+1,jacerr0;
   float_sw4 relax;
   int it = 0;
   relax = m_cirelfact;
 
   icb = m_iStartInt[gc];
   ifb = m_iStartInt[gf];

   ice = m_iEndInt[gc];
   ife = m_iEndInt[gf];
   
   jcb = m_jStartInt[gc];
   jfb = m_jStartInt[gf];

   jce = m_jEndInt[gc];
   jfe = m_jEndInt[gf];

   nkf = m_global_nz[gf];
// material coefficients along the interface (fine grid)
   Sarray Mlfs(m_iStart[gf],m_iEnd[gf],m_jStart[gf],m_jEnd[gf],nkf,nkf);
// make a local copy of Muf to simplify the addition of stretching
   Sarray Mufs(m_iStart[gf],m_iEnd[gf],m_jStart[gf],m_jEnd[gf],nkf,nkf);

#pragma omp parallel for
   for( int j=m_jStart[gf] ; j<=m_jEnd[gf] ; j++ )
      for( int i=m_iStart[gf] ; i<=m_iEnd[gf] ; i++ )
      {
	 Mlfs(i,j,nkf) = (2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(strf_x(i)*strf_y(j)); // (2*mu + lambda)/stretching on the fine grid
         Mufs(i,j,nkf) = Muf(i,j,nkf)/(strf_x(i)*strf_y(j)); // mu/stretching on the fine grid
// include stretching terms in Bf
         for (int c=1; c<=3; c++)
         {
            Bf(c,i,j,nkf) = Bf(c,i,j,nkf)/(strf_x(i)*strf_y(j));
         }
      }

// material coefficients along the interface (coarse grid)
   Sarray Morc(m_iStart[gc],m_iEnd[gc],m_jStart[gc],m_jEnd[gc],1,1);
   Sarray Mlrc(m_iStart[gc],m_iEnd[gc],m_jStart[gc],m_jEnd[gc],1,1);
#pragma omp parallel for
   for( int jc=m_jStart[gc] ; jc<=m_jEnd[gc] ; jc++ )
      for( int ic=m_iStart[gc] ; ic<=m_iEnd[gc] ; ic++ )
      {
	 float_sw4 irho=1/Rhoc(ic,jc,1);
	 Morc(ic,jc,1) = Muc(ic,jc,1)*irho; // mu/rho on the coarse grid (no stretching)
	 Mlrc(ic,jc,1) = (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*irho; // (2*mu+lambda)/rho on the coarse grid
// scale normal stress by stretching
         for (int c=1; c<=3; c++)
            Bc(c,ic,jc,1) = Bc(c,ic,jc,1)/(strc_x(ic)*strc_y(jc));
      }
      
// start iteration
   while( jacerr > m_citol && it < m_cimaxiter )
   {
      float_sw4 rmax[6]={0,0,0,0,0,0};
//
// REMARK: check jump condition in the presence of stretching function;
// stretching function may be different in the fine and coarse grids!
//
// for i=2*ic-1 and j=2*jc-1: Enforce continuity of displacements and normal stresses along the interface

      if (m_croutines)
// optimized version of update for odd i and odd j
	oddIoddJinterpOpt(rmax, Uf.c_ptr(), Muf.c_ptr(), Lambdaf.c_ptr(), Rhof.c_ptr(),
			  Uc.c_ptr(), Muc.c_ptr(), Lambdac.c_ptr(), Rhoc.c_ptr(),
			  Mufs.c_ptr(), Mlfs.c_ptr(),
			  Unextf.c_ptr(), Bf.c_ptr(), Unextc.c_ptr(), Bc.c_ptr(),
			  m_iStart.data(), m_iEnd.data(), m_jStart.data(), m_jEnd.data(), m_kStart.data(), m_kEnd.data(), 
			  m_iStartInt.data(), m_iEndInt.data(), m_jStartInt.data(), m_jEndInt.data(),
			  gf, gc, nkf, mDt, hf, hc, cof, relax,
			  a_strf_x, a_strf_y, a_strc_x, a_strc_y, 
			  m_sbop, m_ghcof);
      else
	oddIoddJinterp(rmax, Uf, Muf, Lambdaf, Rhof, 
		     Uc, Muc, Lambdac, Rhoc,
		     Mufs, Mlfs,
		     Unextf, Bf, Unextc, Bc,
		     m_iStart.data(), m_jStart.data(), m_iStartInt.data(), m_iEndInt.data(), m_jStartInt.data(), m_jEndInt.data(),
		     gf, gc, nkf, mDt, hf, hc, cof, relax,
		     a_strf_x, a_strf_y, a_strc_x, a_strc_y, 
		     m_sbop, m_ghcof);
      
//      
// Enforce continuity of displacements along the interface (for fine ghost points in between coarse points)
//
      oddIevenJinterp(rmax, Uf, Muf, Lambdaf, Rhof, 
      		      Uc, Muc, Lambdac, Rhoc,
      		      Morc, Mlrc,
      		      Unextf, Bf, Unextc, Bc,
      		      m_iStart.data(), m_jStart.data(), m_iStartInt.data(), m_iEndInt.data(), m_jStartInt.data(), m_jEndInt.data(),
      		      gf, gc, nkf, mDt, hf, hc, cof, relax,
      		      a_strf_x, a_strf_y, a_strc_x, a_strc_y, 
      		      m_sbop, m_ghcof);

      evenIoddJinterp(rmax, Uf, Muf, Lambdaf, Rhof, 
      		      Uc, Muc, Lambdac, Rhoc,
      		      Morc, Mlrc,
      		      Unextf, Bf, Unextc, Bc,
      		      m_iStart.data(), m_jStart.data(), m_iStartInt.data(), m_iEndInt.data(), m_jStartInt.data(), m_jEndInt.data(),
      		      gf, gc, nkf, mDt, hf, hc, cof, relax,
      		      a_strf_x, a_strf_y, a_strc_x, a_strc_y, 
      		      m_sbop, m_ghcof);

      evenIevenJinterp(rmax, Uf, Muf, Lambdaf, Rhof, 
		       Uc, Muc, Lambdac, Rhoc,
		       Morc, Mlrc,
		       Unextf, Bf, Unextc, Bc,
		       m_iStart.data(), m_jStart.data(), m_iStartInt.data(), m_iEndInt.data(), m_jStartInt.data(), m_jEndInt.data(),
		       gf, gc, nkf, mDt, hf, hc, cof, relax,
		       a_strf_x, a_strf_y, a_strc_x, a_strc_y, 
		       m_sbop, m_ghcof);

      communicate_array_2d( Uf, gf, nkf+1 );
      communicate_array_2d( Uc, gc, 0 );
      float_sw4 jacerrtmp = 0;
      for (int q=0; q<6; q++)
         jacerrtmp += rmax[q];
      
      MPI_Allreduce( &jacerrtmp, &jacerr, 1, m_mpifloat, MPI_MAX, m_cartesian_communicator );
      if( it == 0 )
	 jacerr0 = jacerr;
      if( jacerr0 > 0 )
	 jacerr = jacerr/jacerr0;
      it++;

   } // end while jacerr > eps (Outer iteration)
   
   if( jacerr > m_citol && proc_zero() )
      cout << "EW::consintp, Warning, no convergence. err = " << jacerr << " tol= " << m_citol << endl;
      
   if( proc_zero() && mVerbose >= 4 )
      cout << "EW::consintp, no of iterations= " << it << " Jac iteration error= " << jacerr << endl;
#undef strc_x
#undef strc_y
#undef strf_x
#undef strf_y
}

//-----------------------------------------------------------------------
void EW::checkintp( Sarray& Uf, Sarray& Unextf, Sarray& Bf, Sarray& Muf, Sarray& Lambdaf, Sarray& Rhof, float_sw4 hf,
		   Sarray& Uc, Sarray& Unextc, Sarray& Bc, Sarray& Muc, Sarray& Lambdac, Sarray& Rhoc, float_sw4 hc,
                    float_sw4 cof, int gc, int gf, int is_periodic[2], float_sw4 time)
{
   //
   // This routine computes the residual in the interpolation conditions
   //
   // time: time corresponding to time level n+1
   //
   // Uf: Fine grid displacement, defined at interior and ghost points, at time level n+1
   // Unextf: Interior contributions (incl. forcing) to the fine grid displacement, at time level n+2
   //
   // Uc: Coarse grid displacement, defined at interior and ghost points, at time level n+1
   // Unextc: Interior contributions (incl. forcing) to the coarse grid displacement, at time level n+2
   //
   // Bf: Interior contributions (incl. forcing) to the fine grid boundary traction, at time level n+1, defined for k=nkf
   // Bc: Interior contributions (incl. forcing) to the coarse grid boundary traction, at time level n+1, defined for k=1
   //
   // At boundaries to the left and right, at least three ghost points are required
   // e.g., domain in i-direction:   i=-2,-1,0,1,2,...,Ni,Ni+1,Ni+2,Ni+3
   // we solve for ghost points at i=2,..,Ni-1, assuming Dirichlet conditions given on i=1,i=Ni and
   // at the ghost points i=-2,-1,0,Ni+1,Ni+2,Ni+3. 
   //
   // The arrays Bf,Unextf, etc are assumed to be computed with all z-ghost points zero,
   // i.e., Uf(i,Nkf+1) was zero for i=-2,-1,..,Ni+3 when Bf and Unextf were evaluated from Uf 
   // (similarly for Unextc and Bc).
   //
   // Before this routine is called, correct boundary values for
   // Uf(-2,Nkf+1),..Uf(1,Nkf+1) (and similarly at the upper i-boundary, and for Uc) must be imposed.
   // In this way the restriction and prolongation stencils can be computed without any special
   // treatment near the (i,j)-boundaries. 

   float_sw4 *a_strc_x = m_sg_str_x[gc];
   float_sw4 *a_strc_y = m_sg_str_y[gc];
   float_sw4 *a_strf_x  = m_sg_str_x[gf];
   float_sw4 *a_strf_y  = m_sg_str_y[gf];

// stretching on the coarse side
#define strc_x(i) a_strc_x[(i-m_iStart[gc])]   
#define strc_y(j) a_strc_y[(j-m_jStart[gc])]   

// stretching on the fine side
#define strf_x(i) a_strf_x[(i-m_iStart[gf])]   
#define strf_y(j) a_strf_y[(j-m_jStart[gf])]   

      
   const float_sw4 i16 = 1.0/16;
   const float_sw4 i256 = 1.0/256;
   const float_sw4 i1024 = 1.0/1024;

   int jcb, jce, icb, ice, jfb, jfe, ifb, ife, nkf;
   float_sw4 nuf = mDt*mDt/(cof*hf*hf); // cof=12 for the predictor, cof=1 for the corrector (argument to this routine)
   float_sw4 nuc = mDt*mDt/(cof*hc*hc);
   float_sw4 ihc = 1/hc, ihf=1/hf;
   float_sw4 jacerr = m_citol+1,jacerr0;
   float_sw4 relax;
   int it = 0;
   relax = m_cirelfact;
 
   icb = m_iStartInt[gc];
   ifb = m_iStartInt[gf];

   ice = m_iEndInt[gc];
   ife = m_iEndInt[gf];
   
   jcb = m_jStartInt[gc];
   jfb = m_jStartInt[gf];

   jce = m_jEndInt[gc];
   jfe = m_jEndInt[gf];

   nkf = m_global_nz[gf];
// material coefficients along the interface (fine grid)
   Sarray Mlfs(m_iStart[gf],m_iEnd[gf],m_jStart[gf],m_jEnd[gf],nkf,nkf);
// make a local copy of Muf to simplify the addition of stretching
   Sarray Mufs(m_iStart[gf],m_iEnd[gf],m_jStart[gf],m_jEnd[gf],nkf,nkf);

#pragma omp parallel for
   for( int j=m_jStart[gf] ; j<=m_jEnd[gf] ; j++ )
      for( int i=m_iStart[gf] ; i<=m_iEnd[gf] ; i++ )
      {
	 Mlfs(i,j,nkf) = (2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(strf_x(i)*strf_y(j)); // (2*mu + lambda)/stretching on the fine grid
         Mufs(i,j,nkf) = Muf(i,j,nkf)/(strf_x(i)*strf_y(j)); // mu/stretching on the fine grid
// include stretching terms in Bf
         for (int c=1; c<=3; c++)
         {
            Bf(c,i,j,nkf) = Bf(c,i,j,nkf)/(strf_x(i)*strf_y(j));
         }
      }

// material coefficients along the interface (coarse grid)
   Sarray Morc(m_iStart[gc],m_iEnd[gc],m_jStart[gc],m_jEnd[gc],1,1);
   Sarray Mlrc(m_iStart[gc],m_iEnd[gc],m_jStart[gc],m_jEnd[gc],1,1);
#pragma omp parallel for
   for( int jc=m_jStart[gc] ; jc<=m_jEnd[gc] ; jc++ )
      for( int ic=m_iStart[gc] ; ic<=m_iEnd[gc] ; ic++ )
      {
	 float_sw4 irho=1/Rhoc(ic,jc,1);
	 Morc(ic,jc,1) = Muc(ic,jc,1)*irho; // mu/rho on the coarse grid (no stretching)
	 Mlrc(ic,jc,1) = (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*irho; // (2*mu+lambda)/rho on the coarse grid
// scale normal stress by stretching
         for (int c=1; c<=3; c++)
            Bc(c,ic,jc,1) = Bc(c,ic,jc,1)/(strc_x(ic)*strc_y(jc));
      }
      
// get exact boundary stresses for twilight
   Sarray BcExact;
   BcExact.define(3,m_iStart[gc], m_iEnd[gc], m_jStart[gc], m_jEnd[gc],1,1);
   BcExact.set_to_zero();
   
// get array pointers for fortran
   float_sw4* mu_ptr    = mMu[gc].c_ptr();
   float_sw4* la_ptr    = mLambda[gc].c_ptr();
   float_sw4 om = m_twilight_forcing->m_omega;
   float_sw4 ph = m_twilight_forcing->m_phase;
   float_sw4 cv = m_twilight_forcing->m_c;
   float_sw4 h  = mGridSize[gc];
   float_sw4* b_ptr = BcExact.c_ptr();
   float_sw4 omstrx = m_supergrid_taper_x[gc].get_tw_omega();
   float_sw4 omstry = m_supergrid_taper_y[gc].get_tw_omega();

   float_sw4* mua_ptr = NULL;
   float_sw4* laa_ptr   = NULL;
   if (m_use_attenuation)
   {
      mua_ptr    = mMuVE[gc][0].c_ptr();
      laa_ptr    = mLambdaVE[gc][0].c_ptr();
   }

   // fill in the interior
   int i1 = icb, i2=ice;
   int j1 = jcb, j2=jce;
   int kic=1;

   if( usingSupergrid() )
   {
      if( m_croutines )
	 twfrsurfzsg_wind_ci( m_iStart[gc], m_iEnd[gc], m_jStart[gc], m_jEnd[gc], m_kStart[gc], m_kEnd[gc],
			      hc, kic, time, om, cv, ph, omstrx, omstry,
			      b_ptr, mu_ptr, la_ptr, m_zmin[gc], i1, i2, j1, j2 );
      else
	 twfrsurfzsg_wind( m_iStart[gc], m_iEnd[gc], m_jStart[gc], m_jEnd[gc], m_kStart[gc], m_kEnd[gc],
			   hc, kic, time, om, cv, ph, omstrx, omstry,
			   b_ptr, mu_ptr, la_ptr, m_zmin[gc], i1, i2, j1, j2 );
      if (m_use_attenuation) // only 1 mechanism with twilight forcing
      {
	 if( m_croutines )
	    twfrsurfzsg_att_wind_ci( m_iStart[gc], m_iEnd[gc], m_jStart[gc], m_jEnd[gc], m_kStart[gc], m_kEnd[gc],
				     hc, kic, time, om, cv, ph, omstrx, omstry,
				     b_ptr, mua_ptr, laa_ptr, m_zmin[gc], i1, i2, j1, j2 );
	 else
	    twfrsurfzsg_att_wind( m_iStart[gc], m_iEnd[gc], m_jStart[gc], m_jEnd[gc], m_kStart[gc], m_kEnd[gc],
				  hc, kic, time, om, cv, ph, omstrx, omstry,
				  b_ptr, mua_ptr, laa_ptr, m_zmin[gc], i1, i2, j1, j2 );
      }
   }
   else
   {
      if( m_croutines )
	 twfrsurfz_wind_ci( m_iStart[gc], m_iEnd[gc], m_jStart[gc], m_jEnd[gc], m_kStart[gc], m_kEnd[gc],
			    hc, kic, time, om, cv, ph,
			    b_ptr, mu_ptr, la_ptr, m_zmin[gc], i1, i2, j1, j2 );
      else
	 twfrsurfz_wind( &m_iStart[gc], &m_iEnd[gc], &m_jStart[gc], &m_jEnd[gc], &m_kStart[gc], &m_kEnd[gc],
			 &hc, &kic, &time, &om, &cv, &ph,
			 b_ptr, mu_ptr, la_ptr, &m_zmin[gc], &i1, &i2, &j1, &j2 );
      if (m_use_attenuation) // only 1 mechanism with twilight forcing
      {
	 if( m_croutines )
	    twfrsurfz_att_wind_ci( m_iStart[gc], m_iEnd[gc], m_jStart[gc], m_jEnd[gc], m_kStart[gc], m_kEnd[gc],
				   hc, kic, time, om, cv, ph, 
				   b_ptr, mua_ptr, laa_ptr, m_zmin[gc], i1, i2, j1, j2 );
	 else
	    twfrsurfz_att_wind( m_iStart[gc], m_iEnd[gc], m_jStart[gc], m_jEnd[gc], m_kStart[gc], m_kEnd[gc],
				hc, kic, time, om, cv, ph, 
				b_ptr, mua_ptr, laa_ptr, m_zmin[gc], i1, i2, j1, j2 );
      }
   }
      
// save displcements and boundary stresses on file (assume 1 proc)
   FILE *fpf=fopen("fine-stress.dat","w");
   FILE *fpc=fopen("coarse-stress.dat","w");
   FILE *fe=fopen("tw-stress.dat","w");

// compute residuals
   float_sw4 rmax[6]={0,0,0,0,0,0};
   float_sw4 err2_stress=0, errmax_stress=0;
   float_sw4 f2_stress=0, fmax_stress=0;
//
// REMARK: check jump condition in the presence of stretching function;
// stretching function may be different in the fine and coarse grids!
//
// for i=2*ic-1 and j=2*jc-1: Enforce continuity of displacements and normal stresses along the interface
   printf("checkintp: icb=%d, ice=%d, jcb=%d, jce=%d\n", icb, ice, jcb, jce);
   
   // Can not parallellize: file write in loop
   // #pragma omp parallel for reduction(+:err2_stress,f2_stress) reduction(max:errmax_stress,fmax_stress)
   for( int jc= jcb ; jc <= jce ; jc++ )
      for( int ic= icb ; ic <= ice ; ic++ )
      {
// i odd, j odd
         int i=2*ic-1, j=2*jc-1;
	 float_sw4 a11, a12, a21, a22, err_c, err_f;
         // setup 2x2 system matrix
         // unknowns: (Uf, Uc)
// eqn 1: continuity of normal stress: NEED stretching
         a11 = 0.25*Mufs(i,j,nkf)*m_sbop[0]*ihf; // ihf = 1/h on the fine grid; Mufs contains stretching
         a12 = Muc(ic,jc,1)*m_sbop[0]*ihc/(strc_x(ic)*strc_y(jc));  // ihc = 1/h on the coarse grid
// eqn 2: continuity of displacement
// nuf = dt^2/(cof * hf^2), nuc = dt^2/(cof * hc^2)
         a21 = nuf/Rhof(i,j,nkf)*Muf(i,j,nkf)*m_ghcof[0]; 
         a22 =-nuc/Rhoc(ic,jc,1)*Muc(ic,jc,1)*m_ghcof[0];
//
//       save results for c=1,2,3
         float_sw4 fstress[3], cstress[3];
         
         for( int c=1 ; c <= 2 ; c++ ) //  the 2 tangential components ? 
         {
// apply the restriction operator to the normal stress on the interface (Bf is on the fine grid)
// scale Bf by 1/strf ?
// fine grid stress
            fstress[c-1] = -0.25*Mufs(i,j,nkf)*m_sbop[0]*ihf*Uf(c,i,j,nkf+1); // -a11*Uf

// add in interpolation terms along the interface  (negative sign because this is the k=Nk boundary)
            fstress[c-1] += -i1024*m_sbop[0]*ihf*(
               Mufs(i-3,j-3,nkf)*Uf(c,i-3,j-3,nkf+1) - 9*Mufs(i-3,j-1,nkf)*Uf(c,i-3,j-1,nkf+1)
               -16*Mufs(i-3,j,  nkf)*Uf(c,i-3,  j,nkf+1) - 9*Mufs(i-3,j+1,nkf)*Uf(c,i-3,j+1,nkf+1)
               +Mufs(i-3,j+3,nkf)*Uf(c,i-3,j+3,nkf+1) +					    
               9*(  -Mufs(i-1,j-3,nkf)*Uf(c,i-1,j-3,nkf+1) + 9*Mufs(i-1,j-1,nkf)*Uf(c,i-1,j-1,nkf+1) 
                    +16*Mufs(i-1,j,  nkf)*Uf(c,i-1,j,  nkf+1) + 9*Mufs(i-1,j+1,nkf)*Uf(c,i-1,j+1,nkf+1) 
                    -Mufs(i-1,j+3,nkf)*Uf(c,i-1,j+3,nkf+1) ) +
               16*(  -Mufs(i,  j-3,nkf)*Uf(c,i,  j-3,nkf+1) + 9*Mufs(i,  j-1,nkf)*Uf(c,i,  j-1,nkf+1) // NOTE: the Uf(i,j) term is in a11
                     + 9*Mufs(i,  j+1,nkf)*Uf(c,i,  j+1,nkf+1)
                     -Mufs(i,  j+3,nkf)*Uf(c,i,  j+3,nkf+1) ) + 
               9*(  -Mufs(i+1,j-3,nkf)*Uf(c,i+1,j-3,nkf+1) + 9*Mufs(i+1,j-1,nkf)*Uf(c,i+1,j-1,nkf+1) 
                    +16*Mufs(i+1,j,  nkf)*Uf(c,i+1,j,  nkf+1) + 9*Mufs(i+1,j+1,nkf)*Uf(c,i+1,j+1,nkf+1)
                    -Mufs(i+1,j+3,nkf)*Uf(c,i+1,j+3,nkf+1) ) +
               Mufs(i+3,j-3,nkf)*Uf(c,i+3,j-3,nkf+1) - 9*Mufs(i+3,j-1,nkf)*Uf(c,i+3,j-1,nkf+1)
               -16*Mufs(i+3,j,  nkf)*Uf(c,i+3,j,  nkf+1) - 9*Mufs(i+3,j+1,nkf)*Uf(c,i+3,j+1,nkf+1)
               +Mufs(i+3,j+3,nkf)*Uf(c,i+3,j+3,nkf+1) );
            
//  add interpolated interior terms
            fstress[c-1]  += i1024*( 
               Bf(c,i-3,j-3,nkf)-9*Bf(c,i-3,j-1,nkf)-16*Bf(c,i-3,j,nkf)-9*Bf(c,i-3,j+1,nkf)+Bf(c,i-3,j+3,nkf)
               +9*(-Bf(c,i-1,j-3,nkf)+9*Bf(c,i-1,j-1,nkf)+16*Bf(c,i-1,j,nkf)+9*Bf(c,i-1,j+1,nkf)-Bf(c,i-1,j+3,nkf))
               +16*(-Bf(c,i,  j-3,nkf)+9*Bf(c,i,  j-1,nkf)+16*Bf(c,i,  j,nkf)+9*Bf(c,i,  j+1,nkf)-Bf(c,i,  j+3,nkf))  // NOTE: Includes Bf(i,j)
               +9*(-Bf(c,i+1,j-3,nkf)+9*Bf(c,i+1,j-1,nkf)+16*Bf(c,i+1,j,nkf)+9*Bf(c,i+1,j+1,nkf)-Bf(c,i+1,j+3,nkf)) +
               Bf(c,i+3,j-3,nkf)-9*Bf(c,i+3,j-1,nkf)-16*Bf(c,i+3,j,nkf)-9*Bf(c,i+3,j+1,nkf)+Bf(c,i+3,j+3,nkf)
               );

// fine stress without interpolation
//            fstress[c-1] = -Mufs(i,j,nkf)*m_sbop[0]*ihf*Uf(c,i,j,nkf+1) + Bf(c,i,  j,nkf);

// coarse stress: ghost point + interior contribution
            cstress[c-1] = Muc(ic,jc,1)*m_sbop[0]*ihc/(strc_x(ic)*strc_y(jc))*Uc(c,ic,jc,0) + Bc(c,ic,jc,1); // a12
            
         } // end for c=1,2

// setup the matrix for the 3rd component of the normal stress (different coefficients)
         // NEED stretching terms in a11 & a12
         a11 = 0.25*Mlfs(i,j,nkf)*m_sbop[0]*ihf; // Mlfs contains stretching
         a12 = (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*m_sbop[0]*ihc/(strc_x(ic)*strc_y(jc));

         a21 = nuf/Rhof(i,j,nkf)*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))*m_ghcof[0];
         a22 =-nuc/Rhoc(ic,jc,1)*(2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*m_ghcof[0];

// apply the restriction operator to the fine grid normal stress grid function (Bf)
// scale Bf for stretching?
//  (negative sign because this is the k=Nk boundary)         
         fstress[2] = -0.25*Mlfs(i,j,nkf)*m_sbop[0]*ihf*Uf(3,i,j,nkf+1); // -a11*Uf

// add in interpolation terms (negative sign because this is the k=Nk boundary)
         fstress[2] += -i1024*m_sbop[0]*ihf*(
            Mlfs(i-3,j-3,nkf)*Uf(3,i-3,j-3,nkf+1) - 9*Mlfs(i-3,j-1,nkf)*Uf(3,i-3,j-1,nkf+1)
            -16*Mlfs(i-3,j,  nkf)*Uf(3,i-3,  j,nkf+1) - 9*Mlfs(i-3,j+1,nkf)*Uf(3,i-3,j+1,nkf+1)
            +Mlfs(i-3,j+3,nkf)*Uf(3,i-3,j+3,nkf+1) +					    
            9*(  -Mlfs(i-1,j-3,nkf)*Uf(3,i-1,j-3,nkf+1) + 9*Mlfs(i-1,j-1,nkf)*Uf(3,i-1,j-1,nkf+1) 
                 +16*Mlfs(i-1,j,  nkf)*Uf(3,i-1,j,  nkf+1) + 9*Mlfs(i-1,j+1,nkf)*Uf(3,i-1,j+1,nkf+1) 
                 -Mlfs(i-1,j+3,nkf)*Uf(3,i-1,j+3,nkf+1) ) +
            16*(  -Mlfs(i,  j-3,nkf)*Uf(3,i,  j-3,nkf+1) + 9*Mlfs(i,  j-1,nkf)*Uf(3,i,  j-1,nkf+1)  // NOTE: excludes Uf(3,i,j)
                  + 9*Mlfs(i,  j+1,nkf)*Uf(3,i,  j+1,nkf+1)
                  -Mlfs(i,  j+3,nkf)*Uf(3,i,  j+3,nkf+1) ) + 
            9*(  -Mlfs(i+1,j-3,nkf)*Uf(3,i+1,j-3,nkf+1) + 9*Mlfs(i+1,j-1,nkf)*Uf(3,i+1,j-1,nkf+1)
                 +16*Mlfs(i+1,j,  nkf)*Uf(3,i+1,j,  nkf+1) + 9*Mlfs(i+1,j+1,nkf)*Uf(3,i+1,j+1,nkf+1)
                 -Mlfs(i+1,j+3,nkf)*Uf(3,i+1,j+3,nkf+1) ) +
            Mlfs(i+3,j-3,nkf)*Uf(3,i+3,j-3,nkf+1) - 9*Mlfs(i+3,j-1,nkf)*Uf(3,i+3,j-1,nkf+1)
            -16*Mlfs(i+3,j,  nkf)*Uf(3,i+3,j,  nkf+1) - 9*Mlfs(i+3,j+1,nkf)*Uf(3,i+3,j+1,nkf+1)
            +Mlfs(i+3,j+3,nkf)*Uf(3,i+3,j+3,nkf+1) );

// add in the interior contribution
         fstress[2] += i1024*( 
            Bf(3,i-3,j-3,nkf)-9*Bf(3,i-3,j-1,nkf)-16*Bf(3,i-3,j,nkf)-9*Bf(3,i-3,j+1,nkf)+Bf(3,i-3,j+3,nkf) +
            9*(-Bf(3,i-1,j-3,nkf)+9*Bf(3,i-1,j-1,nkf)+16*Bf(3,i-1,j,nkf)+9*Bf(3,i-1,j+1,nkf)-Bf(3,i-1,j+3,nkf)) +
            16*(-Bf(3,i,  j-3,nkf)+9*Bf(3,i,  j-1,nkf)+16*Bf(3,i,  j,nkf)+9*Bf(3,i,  j+1,nkf)-Bf(3,i,  j+3,nkf)) + // NOTE: includes Bf(3,i,j)
            9*(-Bf(3,i+1,j-3,nkf)+9*Bf(3,i+1,j-1,nkf)+16*Bf(3,i+1,j,nkf)+9*Bf(3,i+1,j+1,nkf)-Bf(3,i+1,j+3,nkf)) +
            Bf(3,i+3,j-3,nkf)-9*Bf(3,i+3,j-1,nkf)-16*Bf(3,i+3,j,nkf)-9*Bf(3,i+3,j+1,nkf)+Bf(3,i+3,j+3,nkf) );

// fine stress without interpolation
//            fstress[2] = -Mlfs(i,j,nkf)*m_sbop[0]*ihf*Uf(3,i,j,nkf+1) + Bf(3,i,  j,nkf);

// coarse stress: ghost point + interior contributions
         cstress[2] = a12*Uc(3,ic,jc,0) + Bc(3,ic,jc,1); 

         err_c = 0;
         for (int q=0; q<3; q++)
            err_c += SQR(cstress[q]-BcExact(q+1,ic,jc,1));
         err2_stress += err_c;
         err_c = sqrt(err_c);
         if (err_c > errmax_stress) errmax_stress = err_c;

// fine grid err
         err_f = 0;
         for (int q=0; q<3; q++)
            err_f += SQR(fstress[q]-BcExact(q+1,ic,jc,1));
         f2_stress += err_f;
         err_f = sqrt(err_f);
         if (err_f > fmax_stress) fmax_stress = err_f;

         
// save on file
         fprintf(fpf,"%d %d %e %e %e\n", ic, jc, fstress[0], fstress[1], fstress[2]);
         fprintf(fpc,"%d %d %e %e %e\n", ic, jc, cstress[0], cstress[1], cstress[2]);
         
         fprintf(fe,"%d %d %e %e %e\n", ic, jc, BcExact(1,ic,jc,1), BcExact(2,ic,jc,1), BcExact(3,ic,jc,1));
         
      }
   err2_stress = sqrt(err2_stress/((ice-icb+1)*(jce-jcb+1)));
   f2_stress = sqrt(f2_stress/((ice-icb+1)*(jce-jcb+1)));
   printf("checkintp: coinciding points: coarse grid interface stress error: max=%e, L2=%e\n", errmax_stress, err2_stress);
   printf("checkintp: coinciding points: fine grid interface stress error: max=%e, L2=%e\n", fmax_stress, f2_stress);
//      
// Enforce continuity of displacements along the interface (for fine ghost points in between coarse points)
//
// TODO: insert coarse and fine stretching functions below
//
   float_sw4 rmax1=0, rmax2=0, rmax3=0;
#pragma omp parallel for reduction(max:rmax1,rmax2,rmax3)
   for( int j=jfb ; j <= jfe ; j++ )
      for( int i=ifb ; i <= ife ; i++ )
      {
	 int ic, jc;
	 float_sw4 a11, b1, r3;
         if( !( (i % 2 == 1 && j % 2 == 1 ) ) ) // not both i and j are odd (handled above)
         {
// updated components 1,2 of the ghost point value of Uf
            for( int c=1 ; c <= 2 ; c++ )
            {
               if( (j % 2 == 0) && (i % 2 == 1) ) // j is even, i is odd
               {
                  ic = (i+1)/2;
                  jc = j/2;
// All Unextc terms
                  b1 = i16*(-Unextc(c,ic,jc-1,1)+9*(Unextc(c,ic,jc,1)+Unextc(c,ic,jc+1,1))-Unextc(c,ic,jc+2,1));
// All Uc terms
                  b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(c,ic,jc-1,0)*Morc(ic,jc-1,1) + 
                                                   9*Uc(c,ic,jc  ,0)*Morc(ic,jc  ,1) + 
                                                   9*Uc(c,ic,jc+1,0)*Morc(ic,jc+1,1)
                                                   -Uc(c,ic,jc+2,0)*Morc(ic,jc+2,1) );
               }
               if( (j % 2 == 1) && (i % 2 == 0) ) // j is odd, i is even
               {
                  ic = i/2;
                  jc = (j+1)/2;
// All Unextc terms
                  b1 = i16*(-Unextc(c,ic-1,jc,1)+9*(Unextc(c,ic,jc,1)+Unextc(c,ic+1,jc,1))-Unextc(c,ic+2,jc,1));
// All Uc terms
                  b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)+ 
                                                   9*Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+ 
                                                   9*Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1)
                                                   -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1));
               }
               if( (j % 2 == 0) && (i % 2 == 0) ) // i is even, j is even
               {
                  ic = i/2;
                  jc = j/2;
// All Unextc terms
                  b1 = i256*
                     ( Unextc(c,ic-1,jc-1,1)-9*(Unextc(c,ic,jc-1,1)+Unextc(c,ic+1,jc-1,1))+Unextc(c,ic+2,jc-1,1)
                       + 9*(-Unextc(c,ic-1,jc,  1)+9*(Unextc(c,ic,jc,  1)+Unextc(c,ic+1,jc,  1))-Unextc(c,ic+2,jc,  1)  
                            -Unextc(c,ic-1,jc+1,1)+9*(Unextc(c,ic,jc+1,1)+Unextc(c,ic+1,jc+1,1))-Unextc(c,ic+2,jc+1,1))
                       +Unextc(c,ic-1,jc+2,1)-9*(Unextc(c,ic,jc+2,1)+Unextc(c,ic+1,jc+2,1))+Unextc(c,ic+2,jc+2,1) );

// All Uc terms
                  b1 = b1 + nuc*m_ghcof[0]*i256*(
                     Uc(c,ic-1,jc-1,0)*Morc(ic-1,jc-1,1)-9*(Uc(c,ic,  jc-1,0)*Morc(ic,  jc-1,1)+Uc(c,ic+1,jc-1,0)*Morc(ic+1,jc-1,1)) +
                     Uc(c,ic+2,jc-1,0)*Morc(ic+2,jc-1,1)
                     + 9*(
                        -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)
                        +9*(Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1))
                        -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1)  
                        -Uc(c,ic-1,jc+1,0)*Morc(ic-1,jc+1,1)+
                        9*(Uc(c,ic,  jc+1,0)*Morc(ic,  jc+1,1)+Uc(c,ic+1,jc+1,0)*Morc(ic+1,jc+1,1))
                        -Uc(c,ic+2,jc+1,0)*Morc(ic+2,jc+1,1)
                        )
                     + Uc(c,ic-1,jc+2,0)*Morc(ic-1,jc+2,1)
                     -9*(Uc(c,ic,  jc+2,0)*Morc(ic,  jc+2,1) +Uc(c,ic+1,jc+2,0)*Morc(ic+1,jc+2,1))
                     +Uc(c,ic+2,jc+2,0)*Morc(ic+2,jc+2,1)
                     );
               }
               b1 = b1 - Unextf(c,i,j,nkf); 
               a11 = nuf*m_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
//                  a11 = nuf*m_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf))/(strf_x(i)*strf_y(j));
               r3 = Uf(c,i,j,nkf+1); // save old value for relaxation
// update ghost point value Uf(c,i,j,nkf+1)
               Uf(c,i,j,nkf+1) = b1/a11;
               Uf(c,i,j,nkf+1) = relax*Uf(c,i,j,nkf+1)+(1-relax)*r3;
               //		  if( i == 4 && j == 7 && c == 1)
               //		     cout << "in loop " << -a11*Uf(c,i,j,nkf+1) + b1  << endl;
// change in ghost point value
               r3 = r3 - Uf(c,i,j,nkf+1);
	       //               rmax[c-1+3] = rmax[c-1+3] > fabs(r3) ? rmax[c-1+3] : fabs(r3);
	       if( c== 1 )
		  rmax1 = rmax1 > fabs(r3) ? rmax1 : fabs(r3);
	       else
		  rmax2 = rmax2 > fabs(r3) ? rmax2 : fabs(r3);
            } // end for c=1,2
               
// work on componet 3 of the ghost point value of Uf
            if( (j % 2 == 0) && (i % 2 == 1) ) // j even, i odd
            {
               ic = (i+1)/2;
               jc = j/2;
// All Unextc terms
               b1 = i16*(-Unextc(3,ic,jc-1,1)+9*(Unextc(3,ic,jc,1)+Unextc(3,ic,jc+1,1))-Unextc(3,ic,jc+2,1));
// All Uc terms
               b1 = b1 + nuc*m_ghcof[0]*i16*(   - Uc(3,ic,jc-1,0)*Mlrc(ic,jc-1,1) + 
                                                9*Uc(3,ic,jc  ,0)*Mlrc(ic,jc  ,1) + 
                                                9*Uc(3,ic,jc+1,0)*Mlrc(ic,jc+1,1)
                                                -Uc(3,ic,jc+2,0)*Mlrc(ic,jc+2,1) );
            }
            if( (j % 2 == 1) && (i % 2 == 0) ) // j odd, i even
            {
               ic = i/2;
               jc = (j+1)/2;
// All Unextc terms
               b1 = i16*(-Unextc(3,ic-1,jc,1)+9*(Unextc(3,ic,jc,1)+Unextc(3,ic+1,jc,1))-Unextc(3,ic+2,jc,1));
// All Uc terms
               b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(3,ic-1,jc,0)*Mlrc(ic-1,jc,1)+ 
                                                9*Uc(3,ic,  jc,0)*Mlrc(ic,  jc,1)+ 
                                                9*Uc(3,ic+1,jc,0)*Mlrc(ic+1,jc,1)
                                                -Uc(3,ic+2,jc,0)*Mlrc(ic+2,jc,1));
            }
            if( (j % 2 == 0) && (i % 2 == 0) ) // j even, i even
            {
               ic = i/2;
               jc = j/2;
// All Unextc terms
               b1 = i256*
                  ( Unextc(3,ic-1,jc-1,1)-9*(Unextc(3,ic,jc-1,1)+Unextc(3,ic+1,jc-1,1))+Unextc(3,ic+2,jc-1,1)
                    + 9*(-Unextc(3,ic-1,jc,  1)+9*(Unextc(3,ic,jc,  1)+Unextc(3,ic+1,jc,  1))-Unextc(3,ic+2,jc,  1)  
                         -Unextc(3,ic-1,jc+1,1)+9*(Unextc(3,ic,jc+1,1)+Unextc(3,ic+1,jc+1,1))-Unextc(3,ic+2,jc+1,1))
                    +Unextc(3,ic-1,jc+2,1)-9*(Unextc(3,ic,jc+2,1)+Unextc(3,ic+1,jc+2,1))+Unextc(3,ic+2,jc+2,1) );

// All Uc terms
               b1 = b1 + nuc*m_ghcof[0]*i256*(
                  Uc(3,ic-1,jc-1,0)*Mlrc(ic-1,jc-1,1)
                  -9*(Uc(3,ic,  jc-1,0)*Mlrc(ic,  jc-1,1)+
                      Uc(3,ic+1,jc-1,0)*Mlrc(ic+1,jc-1,1))+
                  Uc(3,ic+2,jc-1,0)*Mlrc(ic+2,jc-1,1)
                  + 9*(
                     -Uc(3,ic-1,jc,0)*Mlrc(ic-1,jc,1)+
                     9*( Uc(3,ic,  jc,0)*Mlrc(ic,  jc,1)+
                         Uc(3,ic+1,jc,0)*Mlrc(ic+1,jc,1))
                     -Uc(3,ic+2,jc,0)*Mlrc(ic+2,jc,1)  
                     -Uc(3,ic-1,jc+1,0)*Mlrc(ic-1,jc+1,1)+
                     9*( Uc(3,ic,  jc+1,0)*Mlrc(ic,  jc+1,1)+
                         Uc(3,ic+1,jc+1,0)*Mlrc(ic+1,jc+1,1))
                     -Uc(3,ic+2,jc+1,0)*Mlrc(ic+2,jc+1,1)  )
                  + Uc(3,ic-1,jc+2,0)*Mlrc(ic-1,jc+2,1)
                  -9*(Uc(3,ic,  jc+2,0)*Mlrc(ic,  jc+2,1) +
                      Uc(3,ic+1,jc+2,0)*Mlrc(ic+1,jc+2,1)) +
                  Uc(3,ic+2,jc+2,0)*Mlrc(ic+2,jc+2,1)   );
            } // end  j even, i even
// right hand side is mismatch in displacement                
            b1 = b1 - Unextf(3,i,j,nkf); 
            a11 = nuf*m_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf)); // no str
//               a11 = nuf*m_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf))/(strf_x(i)*strf_y(j)); // with str

            r3 = Uf(3,i,j,nkf+1); // save previous value for relaxation below
// solve for the ghost point value Uf(3,i,j,nkf+1)
            Uf(3,i,j,nkf+1) = b1/a11;
            Uf(3,i,j,nkf+1) = relax*Uf(3,i,j,nkf+1)+(1-relax)*r3;
            r3 = r3 - Uf(3,i,j,nkf+1);
	    //            rmax[2+3] = rmax[2+3] > fabs(r3) ? rmax[2+3] : fabs(r3);
            rmax3 = rmax3 > fabs(r3) ? rmax3 : fabs(r3);

         } // end if not ( i%2=1 &&  j%2=1), i.e, not both i and j are odd

         // (i,j) both odd is handled by the first iteration
            
      } // end for all fine grid points on the interface
      
   rmax[3]=rmax1;
   rmax[4]=rmax2;
   rmax[5]=rmax3;
   communicate_array_2d( Uf, gf, nkf+1 );
   communicate_array_2d( Uc, gc, 0 );
   float_sw4 jacerrtmp = 0;
   for (int q=0; q<6; q++)
      jacerrtmp += rmax[q];
      
   fclose(fpc);
   fclose(fpf);
   fclose(fe);

#undef strc_x
#undef strc_y
#undef strf_x
#undef strf_y
}

