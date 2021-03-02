
#include "EW.h"
#include "MaterialParCart.h"
#include "Require.h"
#include <fstream>
#include <sstream>

#ifdef USE_HDF5
#include <hdf5.h>
#endif

MaterialParCart::MaterialParCart( EW* a_ew, int nx, int ny, int nz, int init, int varcase, 
                                  char* fname, float_sw4 amp, float_sw4 omega )
   : MaterialParameterization( a_ew, fname )
{
   //  VERIFY2( nx > 1 && ny > 1 && nz > 1, "MaterialParCartesian: The grid need at least two ponts in each direction")
     // Material represented on a coarse Cartesian grid, covering the 'active' domain.
     // points are x_0,..,x_{nx+1}, where x_0 and x_{nx+1} are fixed at zero.
   // the parameter vector represents offsets from a reference material, stored in (mRho,mMu,mLambda) in EW.

   m_variables = varcase;
   m_ratio = 1.732;   
   m_init = init;
   m_nx   = nx;
   m_ny   = ny;
   m_nz   = nz;

 // Amplitude and wave number for artificial material tests
   m_amplitude = amp;
   m_omega     = omega;

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
      if( m_myrank == 0 )
         std::cout << "active region, index = " <<
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

   // Determine h, such that x= i*h+xmin, i=0,..,nx+1
   // Note, z_k = k*h + zmin, k=0,..,nz+1, with zmin=-h at surface, 
   //  i.e., k=1 is upper bndry.
   m_hx = (m_xmax-m_xmin)/(nx+1);
   m_hy = (m_ymax-m_ymin)/(ny+1);
   m_hz = (m_zmax-m_zmin)/nz;
   m_zmin = m_zmin-m_hz;

   if( m_myrank == 0 )
   {
     cout << " xmin, xmax = " << m_xmin << " " << m_xmax << " hx = " << m_hx << endl;
     cout << " ymin, ymax = " << m_ymin << " " << m_ymax << " hy = " << m_hy << endl;
     cout << " zmin, zmax = " << m_zmin << " " << m_zmax << " hz = " << m_hz << endl;
   }

   int ncomp=3;
   if( varcase == 3 )
      ncomp=2;
   else if( varcase == 4 )
      ncomp=1;

// Global grid is determined.
   m_global = compute_overlap();

   bool dbg=false;
   if( dbg )
   {
      std::stringstream name;
      name << "dbginfo" << m_myrank << ".dat\0";

      std::ofstream utfil(name.str().c_str());
      utfil << "myid= " << m_myrank << " global=" << m_global << " index bounds, \nI: "
         << m_ib << " " << m_ibint << " " << m_ieint << " " << m_ie <<  "\nJ: "
         << m_jb << " " << m_jbint << " " << m_jeint << " " << m_je <<  "\nK: "
         << m_kb << " " << m_kbint << " " << m_keint << " " << m_ke <<  std::endl;
      utfil.close();
   }
   
   if( !m_global )
   {
      // Distributed arrays
      m_nms = 0;
      m_nmd = ncomp*(m_ieint-m_ibint+1)*(m_jeint-m_jbint+1)*(m_keint-m_kbint+1);
      m_nmd_global = ncomp*nx*ny*nz;
   }
   else
   {
      // Global arrays
      m_nms = ncomp*nx*ny*nz;
      m_nmd = 0;
      m_nmd_global = 0;
   }
}

//-----------------------------------------------------------------------
void MaterialParCart::get_material( int nmd, double* xmd, int nms, double* xms, 
                                    std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu, 
                                    std::vector<Sarray>& a_lambda )

//-----------------------------------------------------------------------
// Get the material corresponding to a parameter representation
//   
// Input: nmd, nms - Declared dimensions of xmd and xms resepectively.
//        xmd      - Distributed parameters.
//        xms      - Shared parameters.
// Output: a_rho, a_mu, a_lambda
//
//-----------------------------------------------------------------------
{
   REQUIRE2( nmd==m_nmd && nms==m_nms, "ERROR in MaterialParCart::get_material " << 
                                       " inconsistent dimensions\n" );

   double* xptr;
   if( m_nmd > 0 )
      xptr = xmd;
   else
      xptr = xms;

   Sarray matcart(3,m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   matcart.set_to_zero();
   if( m_variables == 1 || m_variables == 2 )
   {
      size_t ind =0;
      for( int k=m_kbint ; k <= m_keint ; k++ )
         for( int j=m_jbint ; j <= m_jeint ; j++ )
            for( int i=m_ibint ; i <= m_ieint ; i++ )
            {
               matcart(1,i,j,k) = xptr[3*ind];
               matcart(2,i,j,k) = xptr[3*ind+1];
               matcart(3,i,j,k) = xptr[3*ind+2];
               ind++;
            }
   }
   else if(  m_variables == 3 )
   {
      size_t ind =0;
      for( int k=m_kbint ; k <= m_keint ; k++ )
         for( int j=m_jbint ; j <= m_jeint ; j++ )
            for( int i=m_ibint ; i <= m_ieint ; i++ )
            {
               matcart(2,i,j,k) = xptr[2*ind];
               matcart(3,i,j,k) = xptr[2*ind+1];
               ind++;
            }
   }
   else
   {
      size_t ind =0;
      for( int k=m_kbint ; k <= m_keint ; k++ )
         for( int j=m_jbint ; j <= m_jeint ; j++ )
            for( int i=m_ibint ; i <= m_ieint ; i++ )
            {
               matcart(3,i,j,k) = xptr[ind];
               ind++;
            }
   }
   communicate(matcart);

   for( int g=0 ; g < m_ew->mNumberOfCartesianGrids; g++ )
   {
      a_rho[g].define(m_ew->mRho[g]);
      a_mu[g].define(m_ew->mMu[g]);
      a_lambda[g].define(m_ew->mLambda[g]);

      interpolate( matcart, g, a_rho[g], a_mu[g], a_lambda[g] );

      // Add base material
      float_sw4* mmup  = m_ew->mMu[g].c_ptr();
      float_sw4* mlap  = m_ew->mLambda[g].c_ptr();
      float_sw4* mrhop = m_ew->mRho[g].c_ptr();

      float_sw4* a_mup  = a_mu[g].c_ptr();
      float_sw4* a_lap  = a_lambda[g].c_ptr();
      float_sw4* a_rhop = a_rho[g].c_ptr();

      if( m_variables==1 )
      {
         for( size_t ind=0 ; ind < a_mu[g].npts() ; ind++ )
         {
            a_rhop[ind] = mrhop[ind] + a_rhop[ind];
            a_mup[ind]  = mmup[ind]  + a_mup[ind];
            a_lap[ind]  = mlap[ind]  + a_lap[ind];
         }
      }
      else if( m_variables==2 || m_variables == 3 )
      {
         for( size_t ind=0 ; ind < a_mu[g].npts() ; ind++ )
         {
            float_sw4 cs = sqrt(mmup[ind]/mrhop[ind]) + a_mup[ind];
            float_sw4 cp = sqrt((2*mmup[ind]+mlap[ind])/mrhop[ind])+ a_lap[ind];
            a_rhop[ind] = mrhop[ind] + a_rhop[ind];
            a_mup[ind]  = cs*cs*a_rhop[ind];
            a_lap[ind]  = (cp*cp-2*cs*cs)*a_rhop[ind];
         }
      }
      else
      {
         for( size_t ind=0 ; ind < a_mu[g].npts() ; ind++ )
         {
            float_sw4 cp = sqrt((2*mmup[ind]+mlap[ind])/mrhop[ind])+ a_lap[ind];
            float_sw4 cs = cp/m_ratio;
            a_rhop[ind]  = mrhop[ind];
            a_mup[ind]    = cs*cs*a_rhop[ind];
            a_lap[ind]   = (cp*cp-2*cs*cs)*a_rhop[ind];
         }
      }
      m_ew->communicate_array( a_rho[g], g );
      m_ew->communicate_array( a_mu[g], g );
      m_ew->communicate_array( a_lambda[g], g );
   }
}

//-----------------------------------------------------------------------
void MaterialParCart::find_lims( int ib, int ie, int iepm, int ibpp, 
                                 int& ibint, int& ieint )
{
   if( (ib+iepm) % 2 == 0 )
      ibint = (ib+iepm)/2;
   else
      ibint = (ib+iepm+1)/2;
//   iepmint = ibint-1;

   int ibppint; 
   if( (ibpp+ie) % 2 == 0 )
      ibppint = (ibpp+ie)/2;
   else
      ibppint = (ibpp+ie+1)/2;
   ieint = ibppint-1;
}

//-----------------------------------------------------------------------
bool MaterialParCart::compute_overlap()
{
// Material coarse grid limits in this processor
   int icb=10000000, ice=-100, jcb=10000000, jce=-100, kcb=10000000, kce=-100;

// 1. Interpolate from coarse to fine:
   // find stencil limits: sl,su
   int sl, su;
   float_sw4 wgh[2];
   getwgh( 0.0, wgh, sl, su );
   for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
      float_sw4 h=m_ew->mGridSize[g], zmin=m_ew->m_zmin[g];

    // Active region in fine grid in this processor. These
    // should be interpolated from material grid in this proc.
      int ibact = m_ew->m_iStartAct[g], ieact=m_ew->m_iEndAct[g];
      int jbact = m_ew->m_jStartAct[g], jeact=m_ew->m_jEndAct[g];
      int kbact = m_ew->m_kStartAct[g], keact=m_ew->m_kEndAct[g];   

   // I-direction
      float_sw4 x=(ibact-1)*h;
      float_sw4 q=(x-m_xmin)/m_hx;
      int ic =static_cast<int>(floor(q));
      //      icb = ic+sl;
      icb = std::min(ic+sl,icb);
      x = (ieact-1)*h;
      q=(x-m_xmin)/m_hx;
      ic =static_cast<int>(floor(q));
      //      ice = ic+su;
      ice = std::max(ic+su,ice);

   // J-direction
      float_sw4 y=(jbact-1)*h;
      q=(y-m_ymin)/m_hy;
      int jc =static_cast<int>(floor(q));
      //      jcb = jc+sl;
      jcb = std::min(jc+sl,jcb);
      y = (jeact-1)*h;
      q=(y-m_ymin)/m_hy;
      jc =static_cast<int>(floor(q));
      //      jce = jc+su;
      jce = std::max(jc+su,jce);

   // K-direction
      float_sw4 z=zmin+(kbact-1)*h;
      q=(z-m_zmin)/m_hz;
      int kc =static_cast<int>(floor(q));
      //      kcb = kc+sl;
      kcb = std::min(kc+sl,kcb);
      z=zmin+(keact-1)*h;
      q=(z-m_zmin)/m_hz;
      kc =static_cast<int>(floor(q));
      //      kce = kc+su;
      kce = std::max(kc+su,kce);
   }

   //   std::cout << m_myrank << " First sweep I " << icb << " " << ice << std::endl;
// 2. Interpolate from fine to coarse, adjust bounds if needed:
   int ngh=0;
   bool finetocoarse=false;
   if( finetocoarse )
   {
      bool done = false;
      int maxghost=10;
   //   if( (ice-icb)/2 < maxghost )
   //      maxghost=(ice-icb)/2;
   //   if( (jce-jcb)/2 < maxghost )
   //      maxghost=(jce-jcb)/2;

      while( !done && ngh < maxghost  )
      {
      // try with ngh ghost points, exit if interpolation stencil fits,
      // otherwise increase ngh
         done = true;
         int k=(kce+kcb)/2;
      //      for( int k=kcb+ngh ; k <= kce-ngh ; k++ )
         for( int j=jcb+ngh ; j <= jce-ngh ; j++ )
            for( int i=icb+ngh ; i <= ice-ngh ; i++ )
            {
               double x = m_xmin + i*m_hx;
               double y = m_ymin + j*m_hy;
               //               double z = m_zmin + k*m_hz;
               //               int g, ig, jg, kg;
               int g=0;
               float_sw4 h=m_ew->mGridSize[g], zmin=m_ew->m_zmin[g];
               int ig, jg;
               ig = static_cast<int>(floor(x/h+1));
               jg = static_cast<int>(floor(y/h+1));
               //               kg = static_cast<int>(floor((z-zmin)/h+1));
       //               m_ew->computeNearestLowGridPoint( ig, jg, kg, g, x, y, z );
         // Assume bilinear interpolation from fine grid,
         // only (i,i+1), (j,j+1), (and (k,k+1) ) involved:
               if(!(m_ew->point_in_proc( ig, jg, g) && m_ew->point_in_proc(ig+1,jg+1,g)) )
               {
                  done = false;
                  break;
               }
            }
         if( !done )
            ngh++;
      }
   }
   //   std::cout << m_myrank << " Second sweep ngh= " << ngh << std::endl;
   m_ib= icb-ngh;
   m_ie= ice+ngh;
   m_jb= jcb-ngh;
   m_je= jce+ngh;
   //   m_kb= kcb-ngh;
   //   m_ke= kce+ngh;
   if( m_ew->m_neighbor[0] == MPI_PROC_NULL )
      m_ib   = 0;
   if( m_ew->m_neighbor[1] == MPI_PROC_NULL )
      m_ie = m_nx+1;
   if( m_ew->m_neighbor[2] == MPI_PROC_NULL )
      m_jb   = 0;
   if( m_ew->m_neighbor[3] == MPI_PROC_NULL )
      m_je = m_ny+1;

   // Global index range in k:
   m_kb   = 0;
   m_ke   = m_nz+1;
   m_kbint=m_kb+1;
   m_keint=m_ke-1;

// 3. Compute overlap, and adjust if needed
   MPI_Status status;

   int iepm, ibpp, tag1=667, tag2=668;
   MPI_Sendrecv( &m_ib, 1, MPI_INT, m_ew->m_neighbor[0], tag1, 
                 &ibpp, 1, MPI_INT, m_ew->m_neighbor[1], tag1, 
                 m_ew->m_cartesian_communicator, &status );
   MPI_Sendrecv( &m_ie, 1, MPI_INT, m_ew->m_neighbor[1], tag2, 
                 &iepm, 1, MPI_INT, m_ew->m_neighbor[0], tag2, 
                 m_ew->m_cartesian_communicator, &status );

   if( m_ew->m_neighbor[0] == MPI_PROC_NULL )
      iepm = m_ib+1;
   if( m_ew->m_neighbor[1] == MPI_PROC_NULL )
      ibpp = m_ie;

   int jepm, jbpp;
   MPI_Sendrecv( &m_jb, 1, MPI_INT, m_ew->m_neighbor[2], tag1, 
                 &jbpp, 1, MPI_INT, m_ew->m_neighbor[3], tag1, 
                 m_ew->m_cartesian_communicator, &status );
   MPI_Sendrecv( &m_je, 1, MPI_INT, m_ew->m_neighbor[3], tag2, 
                 &jepm, 1, MPI_INT, m_ew->m_neighbor[2], tag2, 
                 m_ew->m_cartesian_communicator, &status );
   if( m_ew->m_neighbor[2] == MPI_PROC_NULL )
      jepm = m_jb+1;
   if( m_ew->m_neighbor[3] == MPI_PROC_NULL )
      jbpp = m_je;

   int go_global = 0;
   if( ibpp<iepm || jbpp < jepm )
   {
      //      std::cout << "making Material parameterization global  ib,ie = " << m_ib << " "<<m_ie
      //                << " iepm= " << iepm << " ibpp= " << ibpp << std::endl;
      //      std::cout << "              jb,je = " << m_jb << " "<<m_je
      //                << " jepm= " << jepm << " jbpp= " << jbpp << std::endl;
      go_global = 1;
   }
   int tmp=go_global;
   MPI_Allreduce( &tmp, &go_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
   //   // DEBUG
   //   go_global=1;
   if( go_global )
   {
      m_ib = m_jb = m_kb = 0;
      m_ie=m_nx+1;
      m_je=m_ny+1;
      m_ke=m_nz+1;
      m_ibint=m_jbint=m_kbint=1;
      m_ieint=m_nx;
      m_jeint=m_ny;
      m_keint=m_nz;
      if( m_myrank == 0 )
         std::cout << "making Material parameterization global  " << std::endl;
   }
   else
   {
      if( iepm < m_ib-1 || m_ie < ibpp-1 )
         std::cout << "ERROR, gap discovered ib,ie = " << m_ib << " "<<m_ie
                   << " iepm= " << iepm << " ibpp= " << ibpp << std::endl;
      if( jepm < m_jb-1 || m_je < jbpp-1 )
         std::cout << "ERROR, gap discovered jb,je = " << m_jb << " "<<m_je
                   << " jepm= " << iepm << " jbpp= " << ibpp << std::endl;
      find_lims( m_ib, m_ie, iepm, ibpp, m_ibint, m_ieint );
      m_iepm=iepm;
      m_ibpp=ibpp;
      find_lims( m_jb, m_je, jepm, jbpp, m_jbint, m_jeint );
      m_jepm=jepm;
      m_jbpp=jbpp;
   }
   return go_global;
}

//-----------------------------------------------------------------------
void MaterialParCart::getwgh( float_sw4 ai, float_sw4 wgh[2], int& sl, int& su )
{ // Linear interpolation
   wgh[0]=1-ai;
   wgh[1]=ai;
   sl =0;
   su =1;
}

//-----------------------------------------------------------------------
void MaterialParCart::interpolate( Sarray& matcart, int g, Sarray& rho, 
                                   Sarray& mu, Sarray& lambda )
{
   // Interpolate coarse grid material matcart to fine grid rho,mu,lambda.
   //
   // It is assumed that the dimensions of matcart has been determined
   // such that all grid points in the active domain of the computational grid 
   // in this processor can be interpolated from matcart.
   //
   // Input: matcart - Parameter grid material
   //        g       - Grid nr.
   //
   // Output: rho, mu, lambda - Interpolant of matcart on
   //                           computational grid nr. g
   //
   int ibact = m_ew->m_iStartAct[g], ieact=m_ew->m_iEndAct[g];
   int jbact = m_ew->m_jStartAct[g], jeact=m_ew->m_jEndAct[g];
   int kbact = m_ew->m_kStartAct[g], keact=m_ew->m_kEndAct[g];   
   float_sw4 h=m_ew->mGridSize[g], zmin=m_ew->m_zmin[g];
   int sl, su;
   float_sw4 wghdum[20];
   getwgh(0.0,wghdum,sl,su);// call to find sl, su
   int sw = su-sl+1;
   float_sw4* wghx = new float_sw4[sw*(ieact-ibact+1)];
   float_sw4* wghy = new float_sw4[sw*(jeact-jbact+1)];
   float_sw4* wghz = new float_sw4[sw*(keact-kbact+1)];

   for( int i=ibact ; i <= ieact ;i++ )
   {
      float_sw4 x=(i-1)*h;
      float_sw4 q=(x-m_xmin)/m_hx;
      int ic =static_cast<int>(floor(q));
      if( ic+su > m_nx+1 )
         ic = m_nx+1-su;
      if( ic+sl < 0 )
         ic = -sl;
      getwgh(q-ic,&wghx[sw*(i-ibact)],sl,su);
   }
   for( int j=jbact ; j <= jeact ;j++ )
   {
      float_sw4 y=(j-1)*h;
      float_sw4 q=(y-m_ymin)/m_hy;
      int jc =static_cast<int>(floor(q));
      if( jc+su > m_ny+1 )
         jc = m_ny+1-su;
      if( jc+sl < 0 )
         jc = -sl;
      getwgh(q-jc,&wghy[sw*(j-jbact)],sl,su);
   }
   for( int k=kbact ; k <= keact ;k++ )
   {
      float_sw4 z=zmin+(k-1)*h;
      float_sw4 q=(z-m_zmin)/m_hz;
      int kc =static_cast<int>(floor(q));
      if( kc+su > m_nz+1 )
         kc = m_nz+1-su;
      if( kc+sl < 0 )
         kc = -sl;
      getwgh(q-kc,&wghz[sw*(k-kbact)],sl,su);
   }

   rho.set_value(0.0);
   mu.set_value(0.0);
   lambda.set_value(0.0);
   for( int k=kbact; k <= keact ; k++ )
   {
      float_sw4 z=zmin+(k-1)*h;
      int kc =static_cast<int>(floor((z-m_zmin)/m_hz));
      if( kc+su > m_nz+1 )
         kc = m_nz+1-su;
      if( kc+sl < 0 )
         kc = -sl;
      for( int j=jbact; j <= jeact ; j++ )
      {
         float_sw4 y=(j-1)*h;
         int jc =static_cast<int>(floor((y-m_ymin)/m_hy));
         if( jc+su > m_ny+1 )
            jc = m_ny+1-su;
         if( jc+sl < 0 )
            jc = -sl;
         for( int i=ibact; i <= ieact ; i++ )
         {
            float_sw4 x=(i-1)*h;
            int ic =static_cast<int>(floor((x-m_xmin)/m_hx));
            if( ic+su > m_nx+1 )
               ic = m_nx+1-su;
            if( ic+sl < 0 )
               ic = -sl;
            for( int l=sl; l<= su ; l++ )
               for( int n=sl; n<= su ; n++ )
                  for( int m=sl; m<= su ; m++ )
                  {
                     float_sw4 wghtot=wghx[m-sl+sw*(i-ibact)]*wghy[n-sl+sw*(j-jbact)]*wghz[l-sl+sw*(k-kbact)];
                     rho(i,j,k)    += wghtot*matcart(1,ic+m,jc+n,kc+l);
                     mu(i,j,k)     += wghtot*matcart(2,ic+m,jc+n,kc+l);
                     lambda(i,j,k) += wghtot*matcart(3,ic+m,jc+n,kc+l);
                  }
         }
      }
   }
   delete[] wghx;
   delete[] wghy;
   delete[] wghz;
}

//-----------------------------------------------------------------------
void MaterialParCart::get_parameters( int nmd, double* xmd, 
                                      int nms, double* xms, 
                                      std::vector<Sarray>& a_rho, 
                                      std::vector<Sarray>& a_mu, 
                                      std::vector<Sarray>& a_lambda )
{
   if( m_init == 0 )
   {
      for( int i=0 ; i < nms ; i++ )
	 xms[i] = 0;
      for( int i=0 ; i < nmd ; i++ )
	 xmd[i] = 0;
   }
   else if( m_init == 1 )
   {
      //      cout << " hx, hy, hz " << m_hx  <<  " " << m_hy << " " << m_hz << endl;
      //      cout << " nx, ny, nz " << m_nx  <<  " " << m_ny << " " << m_nz << endl;
      //      cout << " xmin, xmax " << m_xmin  <<  " " << m_xmax << endl;
      //      cout << " ymin, ymax " << m_ymin  <<  " " << m_ymax << endl;

   // Test data for sine perturbed material

    // Get base material 
      if( m_variables == 1)
         interpolate_parameters( nmd, xmd, nms, xms, a_rho, a_mu, a_lambda, false );

      //      std::cout << "variables = " << m_variables << std::endl;
      float_sw4 om = m_omega;
      double* xptr = nms>0?xms:xmd;
      size_t ind =0;
      for( int k=m_kbint ; k <= m_keint ; k++ )
	 for( int j=m_jbint ; j <= m_jeint ; j++ )
	    for( int i=m_ibint ; i <= m_ieint ; i++ )
	    {
               bool dbg = i==1 && j==11 && k==5;
	       double x = i*m_hx + m_xmin;
	       double y = j*m_hy + m_ymin;
	       double z = k*m_hz + m_zmin;
               double rhopert = m_amplitude*sin(om*x+0.13)*sin(om*y)*sin(om*z);
               double cspert  = m_amplitude*cos(om*x)*sin(om*y)*cos(om*z+0.01);
               double cppert  = m_amplitude*sin(om*x+0.4)*sin(om*y)*cos(om*z+0.1);

               if( m_variables == 1 )
               {
                  double cs = sqrt( xptr[3*ind+1]/xptr[3*ind] )+cspert;
                  double cp = sqrt((xptr[3*ind+2]+2*xptr[3*ind+1])/xptr[3*ind])+cppert;
                  double mupert     =  xptr[3*ind]*cs*cs-xptr[3*ind+1];
                  double lambdapert =  xptr[3*ind]*(cp*cp-2*cs*cs)-xptr[3*ind+2];
                  xptr[3*ind]   = rhopert;
                  xptr[3*ind+1] = mupert;
                  xptr[3*ind+2] = lambdapert;
                  if( std::isnan(xptr[3*ind]) || std::isnan(xptr[3*ind+1]) || std::isnan(xptr[3*ind]+2) )
                  {
                     std::cout << "ind= " << ind << " " << xptr[3*ind] << " " << xptr[3*ind+1] <<
                        " " << xptr[3*ind+2] << std::endl;
                  }

               }
               else if( m_variables == 2 )
               {
                  xptr[3*ind]   = rhopert;
                  xptr[3*ind+1] = cspert;
                  xptr[3*ind+2] = cppert;
               }
               else if( m_variables == 3 )
               {
                  xptr[2*ind]   = cspert;
                  xptr[2*ind+1] = cppert;
               }
               else if( m_variables == 4 )
                  xptr[ind]   = cppert;
	       ind++;
	    }
   }
   else if( m_init == 2 )
   {
      if( nms > 0 )
         read_parameters( nms, xms );
      else
         read_parameters( nmd, xmd );
   }
   else if( m_init == 3 )
   {
      interpolate_parameters( nmd, xmd, nms, xms, a_rho, a_mu, a_lambda, true );
   }

}

//-----------------------------------------------------------------------
void MaterialParCart::interpolate_parameters( int nmd, double* xmd, 
                                              int nms, double* xms, 
                                              std::vector<Sarray>& a_rho, 
                                              std::vector<Sarray>& a_mu, 
                                              std::vector<Sarray>& a_lambda, 
                                              bool update )
{
   //-----------------------------------------------------------------------
   // Interpolate a material on the computational grid to the parameter grid:
   //
   //   x = I(a_rho-mRho)  if update =true
   //   x = I(a_rho)       if update=false
   //
   // Input: a_rho, a_mu, a_lambda - Material on the SW4 grids
   //
   // Output: xmd or xms - The material parameters in this processor
   //
   //-----------------------------------------------------------------------

   Sarray m_rho(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   Sarray m_cs(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   Sarray m_cp(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   if( m_variables == 1 )
   {
      interpolate_to_coarse( a_rho, a_mu, a_lambda, m_rho, m_cs, m_cp, update );
   }
   else 
   {
      interpolate_to_coarse_vel( a_rho, a_mu, a_lambda, m_rho, m_cs, m_cp );
   }

   double* xptr = m_nms>0? xms:xmd;
   size_t ind=0;
   for( int k=m_kbint ; k <= m_keint; k++ )
      for( int j=m_jbint ; j <= m_jeint ; j++ )
         for( int i=m_ibint ; i <= m_ieint ; i++ )
         {
            if( m_variables == 1 || m_variables == 2)
            {
               xptr[3*ind]   = m_rho(i,j,k);
               xptr[3*ind+1] = m_cs(i,j,k);
               xptr[3*ind+2] = m_cp(i,j,k);
            }
            else if( m_variables==3 )
            {
               xptr[2*ind]   = m_cs(i,j,k);
               xptr[2*ind+1] = m_cp(i,j,k);
            }
            else
               xptr[ind] = m_cp(i,j,k);
            ind++;
         }
}

//-----------------------------------------------------------------------
void MaterialParCart::communicate( Sarray& u )
{
   //--------------------------------------------------------------
   // General ghost point exchange at processor boundaries.
   //
   // Schematically: (o-ghost pt, x-interior pt):
   //
   //   (proc p-1)   iepm        ibpp  (proc p+1)
   //                  v           v
   // ...  x  x  x  x  o           o  o  x  x  x  x  x  ....
   //            o  o  x  x  x  x  x  x  o  o  o
   //            ^     ^              ^        ^
   // myproc(p): ib   ibint          ieint     ie
   //   
   //    myproc receives ib..ibint-1 from p-1, and sends ibint..iepm to p-1
   //                    ieint+1..ie from p+1, and sends ibpp..ieint to p+1
   //
   // Excplicit copy to buffers, not using fancy MPI-datatypes or sendrecv.
   //--------------------------------------------------------------

   if( m_global )
      return;

   const int ncomp=u.m_nc;
   const int ni = m_ie-m_ib+1;
   const int nj = m_je-m_jb+1;
   const int nk = m_ke-m_kb+1;
   float_sw4 *sbuf1, *sbuf2, *rbuf1, *rbuf2;

   MPI_Request req1, req2, req3, req4;
   MPI_Status status;
   int tag1=203, tag2=204;

   size_t npts1low = (m_ibint-m_ib)*nj*nk;
   size_t npts1up  = (m_ie-m_ieint)*nj*nk;

   size_t npts2low = ni*(m_jbint-m_jb)*nk;
   size_t npts2up  = ni*(m_je-m_jeint)*nk;
   size_t nptsmaxlow = max(npts1low,npts2low);
   size_t nptsmaxup  = max(npts1up, npts2up);
  // size: rbuf1 = (I dir) npts1up*ncomp   or (J-dir) npts2up*ncomp
  //       rbuf2 = (I dir) npts1low*ncomp  or (J-dir) npts2low*ncomp
  //       sbuf1 = (I dir) ncomp*ng1*nj*nk or (J-dir) ncomp*ng1*ni*nk
  //       subf2 = (I dir) ncomp*ng2*nj*nk or (J-dir) ncomp*ng2*ni*nk
   size_t nptss1=max((m_iepm-m_ibint+1)*nj*nk,ni*(m_jepm-m_jbint+1)*nk);
   sbuf1 = new float_sw4[ncomp*nptss1];
   rbuf1 = new float_sw4[ncomp*nptsmaxup];

   size_t nptss2=max((m_ieint-m_ibpp+1)*nj*nk,ni*(m_jeint-m_jbpp+1)*nk);
   sbuf2 = new float_sw4[ncomp*nptss2];
   rbuf2 = new float_sw4[ncomp*nptsmaxlow];

// I-direction communication  
   MPI_Irecv( rbuf1, ncomp*npts1up,  m_ew->m_mpifloat, m_ew->m_neighbor[1], tag1,
              m_ew->m_cartesian_communicator, &req1 );

   MPI_Irecv( rbuf2, ncomp*npts1low, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag2,
              m_ew->m_cartesian_communicator, &req2 );

   int ng1=m_iepm-m_ibint+1;
   int ng2=m_ieint-m_ibpp+1;
   if( m_ew->m_neighbor[0] != MPI_PROC_NULL )
   {
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_je; j++ )
               for(int i=m_ibint ; i <= m_iepm; i++ )
               {
                  size_t ind = i-(m_ibint)+ng1*(j-m_jb)+ng1*nj*(k-m_kb);
                  sbuf1[c-1+ncomp*ind]= u(c,i,j,k);
               }
   }
   MPI_Isend( sbuf1, ncomp*ng1*nj*nk, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag1,
              m_ew->m_cartesian_communicator, &req3 );
   if( m_ew->m_neighbor[1] != MPI_PROC_NULL )
   {
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_je; j++ )
               for(int i=m_ibpp; i <= m_ieint; i++ )
               {
                  size_t ind = i-m_ibpp+ng2*(j-m_jb)+ng2*nj*(k-m_kb);
                  sbuf2[c-1+ncomp*ind]= u(c,i,j,k);
               }
   }
   MPI_Isend( sbuf2, ncomp*ng2*nj*nk, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag2,
              m_ew->m_cartesian_communicator, &req4);
   MPI_Wait( &req1, &status );
   ng1 = m_ie-m_ieint;
   if( m_ew->m_neighbor[1] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_je; j++ )
               for(int i=m_ieint+1 ; i <= m_ie; i++ )
               {
                  size_t ind = i-(m_ieint+1)+ng1*(j-m_jb)+ng1*nj*(k-m_kb);
                  u(c,i,j,k) = rbuf1[c-1+ncomp*ind];
               }  
   MPI_Wait( &req2, &status );
   ng2 = m_ibint-m_ib;
   if( m_ew->m_neighbor[0] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_je; j++ )
               for(int i=m_ib ; i <= m_ibint-1; i++ )
               {
                  size_t ind = i-m_ib+ng2*(j-m_jb)+ng2*nj*(k-m_kb);
                  u(c,i,j,k) = rbuf2[c-1+ncomp*ind];
               }

   MPI_Wait( &req3, &status );
   MPI_Wait( &req4, &status );

   // J-direction communication  
   MPI_Irecv( rbuf1, ncomp*npts2up, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag1,
              m_ew->m_cartesian_communicator, &req1 );
   MPI_Irecv( rbuf2, ncomp*npts2low, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag2,
              m_ew->m_cartesian_communicator, &req2 );
   ng1=m_jepm-m_jbint+1;
   ng2=m_jeint-m_jbpp+1;
   if( m_ew->m_neighbor[2] != MPI_PROC_NULL )
      for( int c=1 ; c <= u.m_nc ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jbint ; j <= m_jepm; j++ )
               for(int i=m_ib ; i <= m_ie; i++ )
               {
                  size_t ind = i-m_ib+ni*(j-m_jbint)+ng1*ni*(k-m_kb);
                  sbuf1[c-1+ncomp*ind]= u(c,i,j,k);
               }
   MPI_Isend( sbuf1, ncomp*ng1*ni*nk, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag1,
              m_ew->m_cartesian_communicator, &req3 );
   if( m_ew->m_neighbor[3] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jbpp ; j <= m_jeint; j++ )
               for(int i=m_ib ; i <= m_ie; i++ )
               {
                  size_t ind = i-m_ib+ni*(j-m_jbpp)+ng2*ni*(k-m_kb);
                  sbuf2[c-1+ncomp*ind]= u(c,i,j,k);
               }
   MPI_Isend( sbuf2, ncomp*ni*ng2*nk, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag2,
              m_ew->m_cartesian_communicator, &req4);
   MPI_Wait( &req1, &status );
   ng1 = m_je-m_jeint;
   if( m_ew->m_neighbor[3] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jeint+1 ; j <= m_je; j++ )
               for(int i=m_ib ; i <= m_ie; i++ )
               {
                  size_t ind = i-m_ib + ni*(j-(m_jeint+1))+ng1*ni*(k-m_kb);
                  u(c,i,j,k) = rbuf1[c-1+ncomp*ind];
               }  
   MPI_Wait( &req2, &status );
   ng2 = m_jbint-m_jb;
   if( m_ew->m_neighbor[2] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_jbint-1; j++ )
               for(int i=m_ib ; i <= m_ie; i++ )
               {
                  size_t ind = i-m_ib+ni*(j-m_jb)+ng2*ni*(k-m_kb);
                  u(c,i,j,k) = rbuf2[c-1+ncomp*ind];
               }
    
   MPI_Wait( &req3, &status );
   MPI_Wait( &req4, &status );
   delete[] rbuf1, rbuf2, sbuf1, sbuf2;   
}

//-----------------------------------------------------------------------
void MaterialParCart::communicate_add( Sarray& u )
{
   //--------------------------------------------------------------
   // Get neighbor's ghost points and add to my interior points
   //
   // Schematically: (o-ghost pt, x-interior pt):
   //
   //   (proc p-1)   iepm        ibpp  (proc p+1)
   //                  v           v
   // ...  x  x  x  x  o           o  o  x  x  x  x  x  ....
   //            o  o  x  x  x  x  x  x  o  o  o
   //            ^     ^              ^        ^
   // myproc(p): ib   ibint          ieint     ie
   //   
   // myproc receives ibint..iepm from p-1 and adds the value to ibint..iepm in myproc
   //           sends ieint+1..ie to p+1
   //        receives ibpp..ieint from p+1 and adds the value to ibpp..ieint in myproc
   //           sends ib..ibint-1 to p-1
   //
   //--------------------------------------------------------------

   if( m_global )
      return;

   const int ncomp=u.m_nc;
   const int ni = m_ie-m_ib+1;
   const int nj = m_je-m_jb+1;
   const int nk = m_ke-m_kb+1;
   float_sw4 *sbuf1, *sbuf2, *rbuf1, *rbuf2;

   MPI_Request req1, req2, req3, req4;
   MPI_Status status;
   int tag1=203, tag2=204;

   size_t npts1low = (m_ibint-m_ib)*nj*nk;
   size_t npts1up  = (m_ie-m_ieint)*nj*nk;

   size_t npts2low = ni*(m_jbint-m_jb)*nk;
   size_t npts2up  = ni*(m_je-m_jeint)*nk;
   size_t nptsmaxlow = max(npts1low,npts2low);
   size_t nptsmaxup  = max(npts1up, npts2up);
  // size: rbuf1 = (I dir) npts1up*ncomp   or (J-dir) npts2up*ncomp
  //       rbuf2 = (I dir) npts1low*ncomp  or (J-dir) npts2low*ncomp
  //       sbuf1 = (I dir) ncomp*ng1*nj*nk or (J-dir) ncomp*ng1*ni*nk
  //       subf2 = (I dir) ncomp*ng2*nj*nk or (J-dir) ncomp*ng2*ni*nk
   size_t nptss1=max((m_iepm-m_ibint+1)*nj*nk,ni*(m_jepm-m_jbint+1)*nk);
   rbuf1 = new float_sw4[ncomp*nptss1];
   sbuf1 = new float_sw4[ncomp*nptsmaxup];

   size_t nptss2=max((m_ieint-m_ibpp+1)*nj*nk,ni*(m_jeint-m_jbpp+1)*nk);
   rbuf2 = new float_sw4[ncomp*nptss2];
   sbuf2 = new float_sw4[ncomp*nptsmaxlow];

// I-direction communication  
   MPI_Irecv( rbuf1, ncomp*(m_iepm-m_ibint+1)*nj*nk,  m_ew->m_mpifloat, m_ew->m_neighbor[0], tag1,
              m_ew->m_cartesian_communicator, &req1 );

   MPI_Irecv( rbuf2, ncomp*(m_ieint-m_ibpp+1)*nj*nk, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag2,
              m_ew->m_cartesian_communicator, &req2 );

   int ng1=m_ie-m_ieint;
   int ng2=m_ibint-m_ib;
   if( m_ew->m_neighbor[1] != MPI_PROC_NULL )
   {
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_je; j++ )
               for(int i=m_ieint+1 ; i <= m_ie; i++ )
               {
                  size_t ind = i-(m_ieint+1)+ng1*(j-m_jb)+ng1*nj*(k-m_kb);
                  sbuf1[c-1+ncomp*ind]= u(c,i,j,k);
               }
   }
   MPI_Isend( sbuf1, ncomp*ng1*nj*nk, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag1,
              m_ew->m_cartesian_communicator, &req3 );
   if( m_ew->m_neighbor[0] != MPI_PROC_NULL )
   {
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_je; j++ )
               for(int i=m_ib; i <= m_ibint-1; i++ )
               {
                  size_t ind = i-m_ib+ng2*(j-m_jb)+ng2*nj*(k-m_kb);
                  sbuf2[c-1+ncomp*ind]= u(c,i,j,k);
               }
   }
   MPI_Isend( sbuf2, ncomp*ng2*nj*nk, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag2,
              m_ew->m_cartesian_communicator, &req4);
   MPI_Wait( &req1, &status );
   ng1 = m_iepm-m_ibint+1;
   if( m_ew->m_neighbor[0] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_je; j++ )
               for(int i=m_ibint ; i <= m_iepm; i++ )
               {
                  size_t ind = i-m_ibint+ng1*(j-m_jb)+ng1*nj*(k-m_kb);
                  u(c,i,j,k) += rbuf1[c-1+ncomp*ind];
               }  
   MPI_Wait( &req2, &status );
   ng2 = m_ieint-m_ibpp+1;
   if( m_ew->m_neighbor[1] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_je; j++ )
               for(int i=m_ibpp ; i <= m_ieint; i++ )
               {
                  size_t ind = i-m_ibpp+ng2*(j-m_jb)+ng2*nj*(k-m_kb);
                  u(c,i,j,k) += rbuf2[c-1+ncomp*ind];
               }

   MPI_Wait( &req3, &status );
   MPI_Wait( &req4, &status );

   // J-direction communication  
   MPI_Irecv( rbuf1, ncomp*ni*(m_jepm-m_jbint+1)*nk, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag1,
              m_ew->m_cartesian_communicator, &req1 );
   MPI_Irecv( rbuf2, ncomp*ni*(m_jeint-m_jbpp+1)*nk, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag2,
              m_ew->m_cartesian_communicator, &req2 );
   ng1=m_je-m_jeint;
   ng2=m_jbint-m_jb;
   if( m_ew->m_neighbor[3] != MPI_PROC_NULL )
      for( int c=1 ; c <= u.m_nc ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jeint+1 ; j <= m_je; j++ )
               for(int i=m_ib ; i <= m_ie; i++ )
               {
                  size_t ind = i-m_ib+ni*(j-(m_jeint+1))+ng1*ni*(k-m_kb);
                  sbuf1[c-1+ncomp*ind]= u(c,i,j,k);
               }
   MPI_Isend( sbuf1, ncomp*ng1*ni*nk, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag1,
              m_ew->m_cartesian_communicator, &req3 );
   if( m_ew->m_neighbor[2] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jb ; j <= m_jbint-1; j++ )
               for(int i=m_ib ; i <= m_ie; i++ )
               {
                  size_t ind = i-m_ib+ni*(j-m_jb)+ng2*ni*(k-m_kb);
                  sbuf2[c-1+ncomp*ind]= u(c,i,j,k);
               }
   MPI_Isend( sbuf2, ncomp*ni*ng2*nk, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag2,
              m_ew->m_cartesian_communicator, &req4);
   MPI_Wait( &req1, &status );
   ng1 = m_jepm-m_jbint+1;
   if( m_ew->m_neighbor[2] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jbint ; j <= m_jepm; j++ )
               for(int i=m_ib ; i <= m_ie; i++ )
               {
                  size_t ind = i-m_ib + ni*(j-m_jbint)+ng1*ni*(k-m_kb);
                  u(c,i,j,k) += rbuf1[c-1+ncomp*ind];
               }  
   MPI_Wait( &req2, &status );
   ng2 = m_jeint-m_jbpp+1;
   if( m_ew->m_neighbor[3] != MPI_PROC_NULL )
      for( int c=1 ; c <= ncomp ; c++ )
         for( int k=m_kb ; k <= m_ke; k++ )
            for( int j=m_jbpp ; j <= m_jeint; j++ )
               for(int i=m_ib ; i <= m_ie; i++ )
               {
                  size_t ind = i-m_ib+ni*(j-m_jbpp)+ng2*ni*(k-m_kb);
                  u(c,i,j,k) += rbuf2[c-1+ncomp*ind];
               }
    
   MPI_Wait( &req3, &status );
   MPI_Wait( &req4, &status );
   delete[] rbuf1, rbuf2, sbuf1, sbuf2;   
}

//-----------------------------------------------------------------------
void MaterialParCart::interpolate_gradient( int g, Sarray& a_gradrho, Sarray& a_gradmu, 
                                            Sarray& a_gradlambda, Sarray& gradc )
{
   //-----------------------------------------------------------------------
   // Interpolate the gradient from the SW4 grid to the coarser parameter grid.
   //
   // Input: a_gradrho, a_gradmu, a_gradlambda - Three component gradient
   //                                            on the computational grid
   //
   // Output: gradc - The gradient on the parameter grid, three components as one
   //                 array, size of gradc is 3 x ni x nj x nk
   //
   // This is the transpose of the interpolation in MaterialParCart::interpolate,
   //-----------------------------------------------------------------------

   int ibact = m_ew->m_iStartAct[g], ieact=m_ew->m_iEndAct[g];
   int jbact = m_ew->m_jStartAct[g], jeact=m_ew->m_jEndAct[g];
   int kbact = m_ew->m_kStartAct[g], keact=m_ew->m_kEndAct[g];   
   float_sw4 h=m_ew->mGridSize[g], zmin=m_ew->m_zmin[g];
   int sl, su;
   float_sw4 wghdum[20];
   getwgh(0.0,wghdum,sl,su);// call to find sl, su
   int sw = su-sl+1;
   float_sw4* wghx = new float_sw4[sw*(ieact-ibact+1)];
   float_sw4* wghy = new float_sw4[sw*(jeact-jbact+1)];
   float_sw4* wghz = new float_sw4[sw*(keact-kbact+1)];
   for( int i=ibact ; i <= ieact ;i++ )
   {
      float_sw4 x=(i-1)*h;
      float_sw4 q=(x-m_xmin)/m_hx;
      int ic =static_cast<int>(floor(q));
      if( ic+su > m_nx+1 )
         ic = m_nx+1-su;
      if( ic-sl < 0 )
         ic = sl;
      getwgh(q-ic,&wghx[sw*(i-ibact)],sl,su);
   }
   for( int j=jbact ; j <= jeact ;j++ )
   {
      float_sw4 y=(j-1)*h;
      float_sw4 q=(y-m_ymin)/m_hy;
      int jc =static_cast<int>(floor(q));
      if( jc+su > m_ny+1 )
         jc = m_ny+1-su;
      if( jc-sl < 0 )
         jc = sl;
      getwgh(q-jc,&wghy[sw*(j-jbact)],sl,su);
   }
   for( int k=kbact ; k <= keact ;k++ )
   {
      float_sw4 z=zmin+(k-1)*h;
      float_sw4 q=(z-m_zmin)/m_hz;
      int kc =static_cast<int>(floor(q));
      if( kc+su > m_nz+1 )
         kc = m_nz+1-su;
      if( kc-sl < 0 )
         kc = sl;
      getwgh(q-kc,&wghz[sw*(k-kbact)],sl,su);
   }

   for( int k=kbact; k <= keact ; k++ )
   {
      float_sw4 z=zmin+(k-1)*h;
      int kc =static_cast<int>(floor((z-m_zmin)/m_hz));
      if( kc+su > m_nz+1 )
         kc = m_nz+1-su;
      if( kc-sl < 0 )
         kc = sl;
      for( int j=jbact; j <= jeact ; j++ )
      {
         float_sw4 y=(j-1)*h;
         int jc =static_cast<int>(floor((y-m_ymin)/m_hy));
         if( jc+su > m_ny+1 )
            jc = m_ny+1-su;
         if( jc-sl < 0 )
            jc = sl;
         for( int i=ibact; i <= ieact ; i++ )
         {
            float_sw4 x=(i-1)*h;
            int ic =static_cast<int>(floor((x-m_xmin)/m_hx));            
            if( ic+su > m_nx+1 )
               ic = m_nx+1-su;
            if( ic-sl < 0 )
               ic = sl;
            for( int l=sl; l<= su ; l++ )
               for( int n=sl; n<= su ; n++ )
                  for( int m=sl; m<= su ; m++ )
                  {
                     float_sw4 wghtot=wghx[m-sl+sw*(i-ibact)]*
                                      wghy[n-sl+sw*(j-jbact)]*
                                      wghz[l-sl+sw*(k-kbact)];
                     // Debug
                     //                     gradc(1,ic+m,jc+n,kc+l) += wghtot;
                     //                     gradc(2,ic+m,jc+n,kc+l) += wghtot;
                     //                     gradc(3,ic+m,jc+n,kc+l) += wghtot;
                     //End debug
                     gradc(1,ic+m,jc+n,kc+l) += wghtot*a_gradrho(i,j,k);
                     gradc(2,ic+m,jc+n,kc+l) += wghtot*a_gradmu(i,j,k);
                     gradc(3,ic+m,jc+n,kc+l) += wghtot*a_gradlambda(i,j,k);
                  }
         }
      }
   }
   delete[] wghx;
   delete[] wghy;
   delete[] wghz;
}

//-----------------------------------------------------------------------
void MaterialParCart::get_gradient( int nmd, double* xmd, int nms, double* xms,
                                    double* dfs, double* dfm,
                                    std::vector<Sarray>& a_rho,
                                    std::vector<Sarray>& a_mu,
                                    std::vector<Sarray>& a_lambda,
                                    std::vector<Sarray>& a_gradrho,
                                    std::vector<Sarray>& a_gradmu,
                                    std::vector<Sarray>& a_gradlambda )
{
   //-----------------------------------------------------------------------
   // Computes gradient with respect to the material parameterization from given
   // gradients with respect to the material at grid points.
   // 
   // Input: nmd, nms - Declared sizes of xmd and xms respectively
   //        xmd      - Distributed material parameters
   //        xms      - Shared material parameters
   //        a_rho, a_mu, a_lambda - Current material on the SW4 grid.
   //        a_gradrho, a_gradmu, a_gradlambda - Current material gradient with
   //           respect to (rho,mu,lambda) on the SW4 grid.
   // Output: dfs  - Gradient wrt. shared parameters (xms) 
   //         dfm  - Gradient wrt. distributed parameters (xmd)
   //
   //  Note, only one of dfs and dfm is used.
   //-----------------------------------------------------------------------

   if( m_global )
   {
      get_gradient_shared( nms, xms, dfs, a_rho, a_mu, a_lambda, 
                           a_gradrho, a_gradmu, a_gradlambda );
   }
   else
   {
   if( m_variables != 1 )
   {
   // 1. Transform gradient w.r.t. (rho,mu,lambda) to gradient w.r.t (rho,cs,cp):
      for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
         m_ew->transform_gradient( a_rho[g], a_mu[g], a_lambda[g], 
                                   a_gradrho[g], a_gradmu[g], a_gradlambda[g] );
   }

   Sarray gradc(3,m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
// 2. Gradient of interpolation operator
   gradc.set_to_zero();
   for( int g = 0 ; g < m_ew->mNumberOfGrids ; g++ )
   {
// Chain rule of interpolation relation, e.g., gradc(2,:)=gradmu*d(mu)/dx
      interpolate_gradient( g, a_gradrho[g], a_gradmu[g], a_gradlambda[g], gradc );
   }   
   communicate_add(gradc);
   if( m_variables == 1 || m_variables == 2 )
   {
      size_t ind =0;
      for( int k=m_kbint ; k <= m_keint ; k++ )
         for( int j=m_jbint ; j <= m_jeint ; j++ )
            for( int i=m_ibint ; i <= m_ieint ; i++ )
            {
               dfm[3*ind  ] = gradc(1,i,j,k);
               dfm[3*ind+1] = gradc(2,i,j,k);
               dfm[3*ind+2] = gradc(3,i,j,k);
               ind++;
            }
   }
   else if(  m_variables == 3 )
   {
      size_t ind =0;
      for( int k=m_kbint ; k <= m_keint ; k++ )
         for( int j=m_jbint ; j <= m_jeint ; j++ )
            for( int i=m_ibint ; i <= m_ieint ; i++ )
            {
               dfm[2*ind  ] = gradc(2,i,j,k);
               dfm[2*ind+1] = gradc(3,i,j,k);
               ind++;
            }
   }
   else
   {
      size_t ind =0;
      //      float_sw4 irat2=1/(m_ratio*m_ratio);
      float_sw4 irat=1/m_ratio;
      for( int k=m_kbint ; k <= m_keint ; k++ )
         for( int j=m_jbint ; j <= m_jeint ; j++ )
            for( int i=m_ibint ; i <= m_ieint ; i++ )
            {
               dfm[ind] = irat*gradc(2,i,j,k)+gradc(3,i,j,k);
               //               dfm[ind] = gradc(3,i,j,k);
               //               dfm[ind] = 2*(1-2*irat2)*gradc(3,i,j,k);
               ind++;
            }
   }
   }
}   

//-----------------------------------------------------------------------
void MaterialParCart::get_gradient_shared( int nms, double* xms, double* dfs,
					 std::vector<Sarray>& a_rho,
					 std::vector<Sarray>& a_mu,
					 std::vector<Sarray>& a_lambda,
					 std::vector<Sarray>& a_gradrho,
					 std::vector<Sarray>& a_gradmu,
					 std::vector<Sarray>& a_gradlambda )
{
   //-----------------------------------------------------------------------
   // Computes gradient with respect to the material parameterization from 
   // given gradients with respect to the material at grid points.
   //-----------------------------------------------------------------------

   // 1. Transform gradient w.r.t. (rho,mu,lambda) to gradient w.r.t (rho,cs,cp):
   if( m_variables != 1 )
   {
      for( int g=0 ; g < m_ew->mNumberOfGrids ; g++ )
         m_ew->transform_gradient( a_rho[g], a_mu[g], a_lambda[g], 
                                   a_gradrho[g], a_gradmu[g], a_gradlambda[g] );
   }

   Sarray grho(0,m_nx+1,0,m_ny+1,0,m_nz+1), gmu(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   Sarray glambda(0,m_nx+1,0,m_ny+1,0,m_nz+1);
   grho.set_to_zero();
   gmu.set_to_zero();
   glambda.set_to_zero();
   //   cout << "getgrad nms=" << nms << " nx,ny,nz=" << m_nx <<","<<m_ny <<","<<m_nz <<endl;
   // 2. Gradient of interpolation operator
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

   if( m_variables == 1 || m_variables == 2 )
   {
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
   else if( m_variables == 3 )
   {
      size_t ind =0;
      for( int k=1 ; k <= m_nz ; k++ )
         for( int j=1 ; j <= m_ny ; j++ )
            for( int i=1 ; i <= m_nx ; i++ )
            {
               size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
               dfs[2*ind]   = gmup[indm];
               dfs[2*ind+1] = glambdap[indm];
               ind++;
            }
   }
   else if( m_variables == 4 )
   {
      size_t ind =0;
      for( int k=1 ; k <= m_nz ; k++ )
         for( int j=1 ; j <= m_ny ; j++ )
            for( int i=1 ; i <= m_nx ; i++ )
            {
               size_t indm = i+(m_nx+2)*j + (m_nx+2)*(m_ny+2)*k;
               dfs[ind] = glambdap[indm];
               ind++;
            }
   }
}

//-----------------------------------------------------------------------
void MaterialParCart::set_scalefactors( int nmpars, double* sf, double rho_ref, 
                                        double mu_ref, double lambda_ref,  
                                        double vs_ref, double vp_ref )
{
   if( m_variables == 1 )
   {
      for( int i=0 ; i < nmpars ; i += 3 )
      {
         sf[i]   = rho_ref;
         sf[i+1] = mu_ref;
         sf[i+2] = lambda_ref;
      }
   }
   else if( m_variables == 2 )
   {
      for( int i=0 ; i < nmpars ; i += 3 )
      {
         sf[i]   = rho_ref;
         sf[i+1] = vs_ref;
         sf[i+2] = vp_ref;
      }
   }
   else if( m_variables == 3 )
   {
      for( int i=0 ; i < nmpars ; i += 2 )
      {
         sf[i]   = vs_ref;
         sf[i+1] = vp_ref;
      }
   }
   else if( m_variables == 4 )
   {
      for( int i=0 ; i < nmpars ; i++ )
         sf[i]   = vp_ref;
   }
}

//-----------------------------------------------------------------------
ssize_t MaterialParCart::parameter_index( int ip, int jp, int kp, int grid,
                                          int var )
{
   // Local 1D-index of parameter array at point (ip,jp,kp) where 
   // (ip,jp,kp) are global array indexes. If point not in this
   // processor, -1 is returned.

   int ncomp=3;
   if( m_variables == 3 )
      ncomp=2;
   else if( m_variables == 4 )
      ncomp=1;
   if( m_ibint <= ip && ip <= m_ieint && m_jbint <= jp && jp <= m_jeint
       && m_kbint <= kp && kp <= m_keint && 0 <= var && var <=ncomp-1 )
   {
      int ni=m_ieint-m_ibint+1;
      int nj=m_jeint-m_jbint+1;
      return var+ncomp*(ip-m_ibint+ni*(jp-m_jbint)+ni*nj*(kp-m_kbint) );
   }
   else
      return -1;
}

//-----------------------------------------------------------------------
ssize_t MaterialParCart::local_index( size_t ind_global )
{
   // 1.Transform to global (i,j,k)-space
   int ig, jg, kg, var;
   int ncomp=3;
   if( m_variables == 3 )
      ncomp=2;
   else if( m_variables == 4 )
      ncomp=1;
   var = ind_global % ncomp;
   ind_global = (ind_global-var)/ncomp;
   ig = ind_global % m_nx + 1;
   ind_global = (ind_global-(ig-1))/m_nx;
   jg = ind_global % m_ny + 1;
   ind_global = (ind_global-(jg-1))/m_ny;
   kg = ind_global+1;
   // 2. Reuse local 1D parameter index computation
   return parameter_index(ig,jg,kg,0,var);
}

//-----------------------------------------------------------------------
void MaterialParCart::interpolate_to_coarse( vector<Sarray>& rhogrid, 
                                             vector<Sarray>& mugrid,
                                             vector<Sarray>& lambdagrid,
                                             Sarray& rho, 
                                             Sarray& mu, 
                                             Sarray& lambda,
                                             bool update )
{
//-----------------------------------------------------------------------
// Compute material perturbation on parameter grid from a material that 
// is given on the computational grid.
//
// Computes:  rho := I(rhogrid-mRho)   (and similar for mu and lambda)
//            where I() interpolates from computational grid onto the parameter grid.
//   
//   Input: rhogrid, mugrid, lambdagrid - The material on the computational grids.
//          update - If true, compute rho  := I(rhogrid-mRho)
//                   if false, compute rho := I(rhogrid)
//
//   Output: rho, mu, lambda  - Material perturbation on parameter grid.
//
//-----------------------------------------------------------------------

   int a1=1;
   if( !update )
      a1 = 0;

   int ig, jg, kg, g;
   rho.set_to_zero();
   mu.set_to_zero();
   lambda.set_to_zero();
   //   for( int k=1 ; k <= nz ; k++ )
   //      for( int j=1 ; j <= ny ; j++ )
   //	 for( int i=1 ; i <= nx ; i++ )
   std::vector<Sarray>& mRho = m_ew->mRho;
   std::vector<Sarray>& mMu = m_ew->mMu;
   std::vector<Sarray>& mLambda = m_ew->mLambda;
   for( int k=m_kbint ; k <= m_keint ; k++ )
      for( int j=m_jbint ; j <= m_jeint ; j++ )
	 for( int i=m_ibint ; i <= m_ieint ; i++ )
	 {
            rho(i,j,k)=-1e38;
            mu(i,j,k)=-1e38;
            lambda(i,j,k)=-1e38;
            double x = m_xmin + i*m_hx;
	    double y = m_ymin + j*m_hy;
	    double z = m_zmin + k*m_hz;
            m_ew->computeNearestLowGridPoint( ig, jg, kg, g, x, y, z );
            //            if( interior_point_in_proc( ig, jg, g) )
            if( m_ew->m_iStart[g] <= ig && ig+1 <= m_ew->m_iEnd[g] &&
                m_ew->m_jStart[g] <= jg && jg+1 <= m_ew->m_jEnd[g] &&
                m_ew->m_kStart[g] <= kg && kg+1 <= m_ew->m_kEnd[g] )
	    {
	       double h = m_ew->mGridSize[g];
               double zming = m_ew->m_zmin[g];
	       double wghx = x/h-ig+1;
	       double wghy = y/h-jg+1;
	       double wghz = (z-zming)/h-kg+1;
               rho(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*rhogrid[g](ig,jg,kg)+wghx*rhogrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*rhogrid[g](ig,jg+1,kg)+wghx*rhogrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*rhogrid[g](ig,jg,kg+1)+wghx*rhogrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*rhogrid[g](ig,jg+1,kg+1)+wghx*rhogrid[g](ig+1,jg+1,kg+1)) 
		  - a1*((1-wghy)*(1-wghz)*(
			           (1-wghx)*mRho[g](ig,jg,kg)+wghx*mRho[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mRho[g](ig,jg+1,kg)+wghx*mRho[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mRho[g](ig,jg,kg+1)+wghx*mRho[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mRho[g](ig,jg+1,kg+1)+wghx*mRho[g](ig+1,jg+1,kg+1))  );
               mu(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*mugrid[g](ig,jg,kg)+wghx*mugrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mugrid[g](ig,jg+1,kg)+wghx*mugrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mugrid[g](ig,jg,kg+1)+wghx*mugrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mugrid[g](ig,jg+1,kg+1)+wghx*mugrid[g](ig+1,jg+1,kg+1))
		  - a1*( (1-wghy)*(1-wghz)*(
			           (1-wghx)*mMu[g](ig,jg,kg)+wghx*mMu[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mMu[g](ig,jg+1,kg)+wghx*mMu[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mMu[g](ig,jg,kg+1)+wghx*mMu[g](ig+1,jg,kg+1))+
		      (wghy)*(wghz)*(  (1-wghx)*mMu[g](ig,jg+1,kg+1)+wghx*mMu[g](ig+1,jg+1,kg+1)) );
               lambda(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*lambdagrid[g](ig,jg,kg)+wghx*lambdagrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*lambdagrid[g](ig,jg+1,kg)+wghx*lambdagrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*lambdagrid[g](ig,jg,kg+1)+wghx*lambdagrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*lambdagrid[g](ig,jg+1,kg+1)+wghx*lambdagrid[g](ig+1,jg+1,kg+1))
		  -a1*((1-wghy)*(1-wghz)*(
			           (1-wghx)*mLambda[g](ig,jg,kg)+wghx*mLambda[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mLambda[g](ig,jg+1,kg)+wghx*mLambda[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mLambda[g](ig,jg,kg+1)+wghx*mLambda[g](ig+1,jg,kg+1))+
		    (wghy)*(wghz)*(  (1-wghx)*mLambda[g](ig,jg+1,kg+1)+wghx*mLambda[g](ig+1,jg+1,kg+1)));
      // Could do trilinear intp through ig,ig+1,jg,jg+1,kg,kg+1 instead
	       //               rho(i,j,k)    =    rhogrid[g](ig,jg,kg)-mRho[g](ig,jg,kg);
	       //               mu(i,j,k)     =     mugrid[g](ig,jg,kg)-mMu[g](ig,jg,kg);
	       //               lambda(i,j,k) = lambdagrid[g](ig,jg,kg)-mLambda[g](ig,jg,kg);
	    }
	 }
   if( m_global )
   {
      Sarray tmp;
      tmp.copy(rho);
      MPI_Allreduce(tmp.c_ptr(),rho.c_ptr(),rho.npts(),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      tmp.copy(mu);
      MPI_Allreduce(tmp.c_ptr(),mu.c_ptr(),mu.npts(),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      tmp.copy(lambda);
      MPI_Allreduce(tmp.c_ptr(),lambda.c_ptr(),lambda.npts(),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   }
   else
   {
      communicate(rho);
      communicate(mu);
      communicate(lambda);
   }
}

//-----------------------------------------------------------------------
void MaterialParCart::interpolate_to_coarse_vel( vector<Sarray>& rhogrid, 
                                                 vector<Sarray>& mugrid,
                                                 vector<Sarray>& lambdagrid,
                                                 Sarray& rho, Sarray& cs, Sarray& cp )
{
// Compute material perturbation on parameter grid from a material that is given on the computational grid.
//
// Computes:  rho := I(rhogrid-mRho)   (and similar for mu and lambda)
//            cs  := I(csgrid-mCs)
//            cp  := I(cpgrid-mCp)
//            where I() interpolates from computational grid onto the parameter grid.
//            csgrid and cpgrid are obtained by transforming (rhogrid,mugrid,lambdagrid)
//            mCs and mCp are obtained by transforming (mRho, mMu, mLambda)
//   
//   Input: rhogrid, mugrid, lambdagrid - The material (rho,mu,lambda) on
//                                        the computational grids.
//
//   Output: rho, cs, cp  - Material velocity perturbation on parameter grid.
//

   int ig, jg, kg, g;
   rho.set_to_zero();
   cs.set_to_zero();
   cp.set_to_zero();
   int ngrids=m_ew->mNumberOfGrids;
   vector<bool> done(ngrids);
   for( int g=0 ; g < ngrids ; g++ )
      done[g] = false;

   std::vector<Sarray>& mRho    = m_ew->mRho;
   std::vector<Sarray>& mMu     = m_ew->mMu;
   std::vector<Sarray>& mLambda = m_ew->mLambda;

   vector<Sarray> csdiff(ngrids);
   vector<Sarray> cpdiff(ngrids);
   for( int k=m_kbint ; k <= m_keint ; k++ )
      for( int j=m_jbint ; j <= m_jeint ; j++ )
	 for( int i=m_ibint ; i <= m_ieint ; i++ )
	 {
            rho(i,j,k)= -1e38;
            cs(i,j,k) = -1e38;
            cp(i,j,k) = -1e38;
            double x = m_xmin + i*m_hx;
	    double y = m_ymin + j*m_hy;
	    double z = m_zmin + k*m_hz;
            m_ew->computeNearestLowGridPoint( ig, jg, kg, g, x, y, z );
               //            if( interior_point_in_proc( ig, jg, g) )
            if( m_ew->m_iStart[g] <= ig && ig+1 <= m_ew->m_iEnd[g] &&
                m_ew->m_jStart[g] <= jg && jg+1 <= m_ew->m_jEnd[g] &&
                m_ew->m_kStart[g] <= kg && kg+1 <= m_ew->m_kEnd[g] )
	    {
	       if( !done[g] )
	       {
		  csdiff[g].define(mugrid[g]);
		  cpdiff[g].define(mugrid[g]);
		  for( int k1=rhogrid[g].m_kb ; k1 <= rhogrid[g].m_ke ; k1++)
		     for( int j1=rhogrid[g].m_jb ; j1 <= rhogrid[g].m_je ; j1++)
			for( int i1=rhogrid[g].m_ib ; i1 <= rhogrid[g].m_ie ; i1++)
			{
			   csdiff[g](i1,j1,k1)=sqrt(mugrid[g](i1,j1,k1)/rhogrid[g](i1,j1,k1))-
			                       sqrt(mMu[g](i1,j1,k1)/mRho[g](i1,j1,k1));
			   cpdiff[g](i1,j1,k1)=
			     sqrt((2*mugrid[g](i1,j1,k1)+lambdagrid[g](i1,j1,k1))/rhogrid[g](i1,j1,k1)) -
		             sqrt((2*mMu[g](i1,j1,k1)   +mLambda[g](i1,j1,k1)   )/mRho[g](i1,j1,k1));
			}
		  done[g] = true;
	       }
	       double h = m_ew->mGridSize[g];
	       double wghx = x/h-ig+1;
	       double wghy = y/h-jg+1;
	       double wghz = z/h-kg+1;
               rho(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*rhogrid[g](ig,jg,kg)+wghx*rhogrid[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*rhogrid[g](ig,jg+1,kg)+wghx*rhogrid[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*rhogrid[g](ig,jg,kg+1)+wghx*rhogrid[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*rhogrid[g](ig,jg+1,kg+1)+wghx*rhogrid[g](ig+1,jg+1,kg+1)) 
		  -      ((1-wghy)*(1-wghz)*(
			           (1-wghx)*mRho[g](ig,jg,kg)+wghx*mRho[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*mRho[g](ig,jg+1,kg)+wghx*mRho[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*mRho[g](ig,jg,kg+1)+wghx*mRho[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*mRho[g](ig,jg+1,kg+1)+wghx*mRho[g](ig+1,jg+1,kg+1))  );
               cs(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*csdiff[g](ig,jg,kg)+wghx*csdiff[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*csdiff[g](ig,jg+1,kg)+wghx*csdiff[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*csdiff[g](ig,jg,kg+1)+wghx*csdiff[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*csdiff[g](ig,jg+1,kg+1)+wghx*csdiff[g](ig+1,jg+1,kg+1));
               cp(i,j,k) = (1-wghy)*(1-wghz)*(
			           (1-wghx)*cpdiff[g](ig,jg,kg)+wghx*cpdiff[g](ig+1,jg,kg))+
		  (wghy)*(1-wghz)*((1-wghx)*cpdiff[g](ig,jg+1,kg)+wghx*cpdiff[g](ig+1,jg+1,kg))+
		  (1-wghy)*(wghz)*((1-wghx)*cpdiff[g](ig,jg,kg+1)+wghx*cpdiff[g](ig+1,jg,kg+1))+
		  (wghy)*(wghz)*(  (1-wghx)*cpdiff[g](ig,jg+1,kg+1)+wghx*cpdiff[g](ig+1,jg+1,kg+1));
	    }
	 }
   if( m_global )
   {
      Sarray tmp;
      tmp.copy(rho);
      MPI_Allreduce(tmp.c_ptr(),rho.c_ptr(),rho.npts(),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      tmp.copy(cs);
      MPI_Allreduce(tmp.c_ptr(),cs.c_ptr(),cs.npts(),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      tmp.copy(cp);
      MPI_Allreduce(tmp.c_ptr(),cp.c_ptr(),cp.npts(),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   }
   else
   {
      communicate(rho);
      communicate(cs);
      communicate(cp);
   }
}

void MaterialParCart::write_dfm_hdf5(double *dfm, std::string fname, MPI_Comm comm)
{
#ifdef USE_HDF5
    hid_t fid, dset, fapl, dxpl, dspace, mspace;
    hsize_t dims[3], start[3], count[3];
    int ret;

    int ncomp=3;
    if(m_variables == 3)
        ncomp=2;
    else if(m_variables == 4)
        ncomp=1;

    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);

    dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

    fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    if (fid < 0) {
        std::cout << "Error with dfm file createion!" << endl;
        return;
    }

    dims[0] = m_nz;
    dims[1] = m_ny;
    dims[2] = m_nx*ncomp;

    start[0] = m_kbint - 1;
    start[1] = m_jbint - 1;
    start[2] = (m_ibint - 1)*ncomp;

    count[0] = m_keint - m_kbint + 1;
    count[1] = m_jeint - m_jbint + 1;
    count[2] = (m_ieint - m_ibint + 1) * ncomp;

    hsize_t total = count[0]*count[1]*count[2];
    printf("Rank %d, Global: %llu %llu %llu, ncomp %d\n", m_myrank, dims[0], dims[1], dims[2], ncomp);
    printf("Rank %d, Local start: %llu %llu %llu count: %llu %llu %llu, total %llu\n", 
            m_myrank, start[0], start[1], start[2], count[0], count[1], count[2], total);
    printf("Rank %d, %f %f, %f %f\n", m_myrank, dfm[0], dfm[1], dfm[total-2], dfm[total-1]);
    fflush(stdout);

    dspace = H5Screate_simple(3, dims, NULL);
    mspace = H5Screate_simple(3, count, NULL);

    dset = H5Dcreate(fid, "dfm", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        std::cout << "Error with dfm dset createion!" << endl;
        return;
    }

    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL, count, NULL);

    ret = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, dspace, dxpl, dfm);
    if (dset < 0) {
        std::cout << "Error with dfm dset createion!" << endl;
    }

    H5Dclose(dset);
    H5Sclose(dspace);
    H5Sclose(mspace);
    H5Pclose(fapl);
    H5Pclose(dxpl);
    H5Fclose(fid);
#else
   if( m_myrank == 0 )
       std::cout << "Warning: calling write_dfm_hdf5 without compiling SW4 with HDF5!" << endl;
#endif
}
