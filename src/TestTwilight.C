#include "TestTwilight.h"

#if defined(SW4_USE_RAJA) // SW4 backend
#include "caliper.h"

#else // SW4 backend
#endif // SW4 backend
TestTwilight::TestTwilight( float_sw4 omega, float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase,
                            float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda ) :
   m_omega(omega),
   m_c(c),
   m_phase(phase),
   m_momega(momega),
   m_mphase(mphase),
   m_amprho(amprho),
   m_ampmu(ampmu),
   m_amplambda(amplambda),
   m_sw4twilight(true)
{}

void TestTwilight::get_rho( Sarray& rho, Sarray& x, Sarray& y, Sarray& z )
{

#if defined(SW4_USE_RAJA) // SW4 backend
SW4_MARK_FUNCTION;

#else // SW4 backend
#endif // SW4 backend
if( m_sw4twilight )
   {
      for( int k=rho.m_kb ; k <= rho.m_ke ; k++ )
         for( int j=rho.m_jb ; j <= rho.m_je ; j++ )
            for( int i=rho.m_ib ; i <= rho.m_ie ; i++ )
               rho(i,j,k) = m_amprho*(2 +
                                   sin(m_momega*x(i,j,k)+m_mphase)*
                                   cos(m_momega*y(i,j,k)+m_mphase)*
                                   sin(m_momega*z(i,j,k)+m_mphase) );
   }
   else
   {
      // Test code material
      for( int k=rho.m_kb ; k <= rho.m_ke ; k++ )
         for( int j=rho.m_jb ; j <= rho.m_je ; j++ )
            for( int i=rho.m_ib ; i <= rho.m_ie ; i++ )
               rho(i,j,k) = 2 + sin(x(i,j,k)+0.3)*sin(y(i,j,k)+0.3)*sin(z(i,j,k)-0.2);
   }
}

void TestTwilight::get_mula( Sarray& mu, Sarray& lambda, Sarray& x, Sarray& y, Sarray& z )
{

#if defined(SW4_USE_RAJA) // SW4 backend
SW4_MARK_FUNCTION;

#else // SW4 backend
#endif // SW4 backend
if( m_sw4twilight )
   {
      for( int k=mu.m_kb ; k <= mu.m_ke ; k++ )
         for( int j=mu.m_jb ; j <= mu.m_je ; j++ )
            for( int i=mu.m_ib ; i <= mu.m_ie ; i++ )
            {
               mu(i,j,k) = m_ampmu*(3 +
                                 cos(m_momega*x(i,j,k)+m_mphase)*
                                 sin(m_momega*y(i,j,k)+m_mphase)*
                                 sin(m_momega*z(i,j,k)+m_mphase) );
               lambda(i,j,k)= m_amplambda*(2 +
                                        sin(m_momega*x(i,j,k)+m_mphase)*
                                        sin(m_momega*y(i,j,k)+m_mphase)*
                                        cos(m_momega*z(i,j,k)+m_mphase) );
            }
   }
   else
   {
      // Test code material
      for( int k=mu.m_kb ; k <= mu.m_ke ; k++ )
         for( int j=mu.m_jb ; j <= mu.m_je ; j++ )
            for( int i=mu.m_ib ; i <= mu.m_ie ; i++ )
            {
               mu(i,j,k) = 3.0 +sin(3*x(i,j,k)+0.1)*sin(3*y(i,j,k)+0.1)*sin(z(i,j,k));
               lambda(i,j,k) = 21.0 + cos(x(i,j,k)+0.1)*cos(y(i,j,k)+0.1)*pow(sin(3*z(i,j,k)),2);
            }
   }
}

#if defined(SW4_USE_RAJA) // SW4 backend
void TestTwilight::get_ubnd(Sarray& u_in, Sarray& x_in, Sarray& y_in, Sarray& z_in,
                            float_sw4 t, int npts, int sides[6]) {
  SW4_MARK_FUNCTION;

  SView& u = u_in.getview();
  SView& x = x_in.getview();
  SView& y = y_in.getview();
  SView& z = z_in.getview();

#else // SW4 backend
void TestTwilight::get_ubnd( Sarray& u, Sarray& x, Sarray& y, Sarray& z, float_sw4 t, int npts, int sides[6] )
{

#endif // SW4 backend
for( int s=0 ; s < 6 ; s++ )
      if( sides[s]==1 )
      {

#if defined(SW4_USE_RAJA) // SW4 backend
int kb = u_in.m_kb, ke = u_in.m_ke, jb = u_in.m_jb, je = u_in.m_je, ib = u_in.m_ib,
          ie = u_in.m_ie;

#else // SW4 backend
int kb=u.m_kb, ke=u.m_ke, jb=u.m_jb, je=u.m_je, ib=u.m_ib, ie=u.m_ie;

#endif // SW4 backend
if( s == 0 )
            ie = ib+npts-1;
         if( s == 1 )
            ib = ie-npts+1;
         if( s == 2 )
            je = jb+npts-1;
         if( s == 3 )
            jb = je-npts+1;
         if( s == 4 )
         {
            ke = kb+npts-1;

#if defined(SW4_USE_RAJA) // SW4 backend
if (ke > u_in.m_ke) ke = u_in.m_ke;

#else // SW4 backend
if( ke > u.m_ke )
               ke = u.m_ke;

#endif // SW4 backend
}
         if( s == 5 )
         {
            kb = ke-npts+1;

#if defined(SW4_USE_RAJA) // SW4 backend
if (kb < u_in.m_kb) kb = u_in.m_kb;

#else // SW4 backend
if( kb < u.m_kb )
               kb = u.m_kb;

#endif // SW4 backend
}

#if defined(SW4_USE_RAJA) // SW4 backend
auto lm_omega = m_omega;
      auto lm_phase = m_phase;
      auto lm_c = m_c;
      RAJA::RangeSegment k_range(kb, ke + 1);
      RAJA::RangeSegment j_range(jb, je + 1);
      RAJA::RangeSegment i_range(ib, ie + 1);
      RAJA::kernel<TGU_POL_ASYNC>(RAJA::make_tuple(k_range, j_range, i_range),
                                  [=] RAJA_DEVICE(int k, int j, int i) {
                                    //for (int k = kb; k <= ke; k++)
                                    //for (int j = jb; j <= je; j++)
                                    //for (int i = ib; i <= ie; i++) {
            u(1, i, j, k) = sin(lm_omega * (x(i, j, k) - lm_c * t)) *
                            sin(lm_omega * y(i, j, k) + lm_phase) *
                            sin(lm_omega * z(i, j, k) + lm_phase);
            u(2, i, j, k) = sin(lm_omega * x(i, j, k) + lm_phase) *
                            sin(lm_omega * (y(i, j, k) - lm_c * t)) *
                            sin(lm_omega * z(i, j, k) + lm_phase);
            u(3, i, j, k) = sin(lm_omega * x(i, j, k) + lm_phase) *
                            sin(lm_omega * y(i, j, k) + lm_phase) *
                            sin(lm_omega * (z(i, j, k) - lm_c * t));
                                  });

#else // SW4 backend
for( int k=kb ; k <= ke ; k++ )
            for( int j=jb ; j <= je ; j++ )
               for( int i=ib ; i <= ie ; i++ )
               {
                  u(1,i,j,k) = sin(m_omega*(x(i,j,k)-m_c*t))*sin(m_omega*y(i,j,k)+m_phase)*sin(m_omega*z(i,j,k)+m_phase);
                  u(2,i,j,k) = sin(m_omega*x(i,j,k)+m_phase)*sin(m_omega*(y(i,j,k)-m_c*t))*sin(m_omega*z(i,j,k)+m_phase);
                  u(3,i,j,k) = sin(m_omega*x(i,j,k)+m_phase)*sin(m_omega*y(i,j,k)+m_phase)*sin(m_omega*(z(i,j,k)-m_c*t));

#endif // SW4 backend
}

#if defined(SW4_USE_RAJA) // SW4 backend
SYNC_STREAM;
#else // SW4 backend
}
#endif // SW4 backend
}

void TestTwilight::get_ubnd( Sarray& u, float_sw4 h, float_sw4 zmin, float_sw4 t, int npts, int sides[6] )
{

#if defined(SW4_USE_RAJA) // SW4 backend
SW4_MARK_FUNCTION;
  SYNC_STREAM;

#else // SW4 backend
#endif // SW4 backend
for( int s=0 ; s < 6 ; s++ )
      if( sides[s]==1 )
      {
         int kb=u.m_kb, ke=u.m_ke, jb=u.m_jb, je=u.m_je, ib=u.m_ib, ie=u.m_ie;
         if( s == 0 )
            ie = ib+npts-1;
         if( s == 1 )
            ib = ie-npts+1;
         if( s == 2 )
            je = jb+npts-1;
         if( s == 3 )
            jb = je-npts+1;
         if( s == 4 )
         {
            ke = kb+npts-1;
            if( ke > u.m_ke )
               ke = u.m_ke;
         }
         if( s == 5 )
         {
            kb = ke-npts+1;
            if( kb < u.m_kb )
               kb = u.m_kb;
         }
         for( int k=kb ; k <= ke ; k++ )
            for( int j=jb ; j <= je ; j++ )
               for( int i=ib ; i <= ie ; i++ )
               {
                  float_sw4 x = h*(i-1), y=h*(j-1), z=h*(k-1)+zmin;
                  u(1,i,j,k) = sin(m_omega*(x-m_c*t))*sin(m_omega*y+m_phase)*sin(m_omega*z+m_phase);
                  u(2,i,j,k) = sin(m_omega*x+m_phase)*sin(m_omega*(y-m_c*t))*sin(m_omega*z+m_phase);
                  u(3,i,j,k) = sin(m_omega*x+m_phase)*sin(m_omega*y+m_phase)*sin(m_omega*(z-m_c*t));
               }
      }
}

void TestTwilight::get_mula_att( Sarray& muve, Sarray& lambdave, Sarray& x, Sarray& y, Sarray& z )
{

#if defined(SW4_USE_RAJA) // SW4 backend
SW4_MARK_FUNCTION;

#else // SW4 backend
#endif // SW4 backend
if( m_sw4twilight )
   {
      for( int k=muve.m_kb ; k <= muve.m_ke ; k++ )
         for( int j=muve.m_jb ; j <= muve.m_je ; j++ )
            for( int i=muve.m_ib ; i <= muve.m_ie ; i++ )
            {
               muve(i,j,k) = m_ampmu*(1.5 + 0.5*
                                 cos(m_momega*x(i,j,k)+m_mphase)*
                                 cos(m_momega*y(i,j,k)+m_mphase)*
                                 sin(m_momega*z(i,j,k)+m_mphase) );
               lambdave(i,j,k)= m_amplambda*(0.5 + 0.25*
                                        sin(m_momega*x(i,j,k)+m_mphase)*
                                        cos(m_momega*y(i,j,k)+m_mphase)*
                                        sin(m_momega*z(i,j,k)+m_mphase) );
            }
   }
}

#if defined(SW4_USE_RAJA) // SW4 backend
void TestTwilight::get_bnd_att(Sarray& AlphaVE_in, Sarray& x_in, Sarray& y_in, Sarray& z_in,
                               float_sw4 t, int npts, int sides[6]) {
  SW4_MARK_FUNCTION;
  SView& AlphaVE = AlphaVE_in.getview();
  SView& x = x_in.getview();
  SView& y = y_in.getview();
  SView& z = z_in.getview();
  //std::cout << "WARNING TestTwilight::get_bnd_att running on CPU\n"
  //          << std::flush;

#else // SW4 backend
void TestTwilight::get_bnd_att( Sarray& AlphaVE, Sarray& x, Sarray& y, Sarray& z, float_sw4 t, int npts, int sides[6] )
{

#endif // SW4 backend
for( int s=0 ; s < 6 ; s++ )
      if( sides[s]==1 )
      {

#if defined(SW4_USE_RAJA) // SW4 backend
int kb = AlphaVE_in.m_kb, ke = AlphaVE_in.m_ke, jb = AlphaVE_in.m_jb,
          je = AlphaVE_in.m_je, ib = AlphaVE_in.m_ib, ie = AlphaVE_in.m_ie;

#else // SW4 backend
int kb=AlphaVE.m_kb, ke=AlphaVE.m_ke, jb=AlphaVE.m_jb, je=AlphaVE.m_je, ib=AlphaVE.m_ib, ie=AlphaVE.m_ie;

#endif // SW4 backend
if( s == 0 )
            ie = ib+npts-1;
         if( s == 1 )
            ib = ie-npts+1;
         if( s == 2 )
            je = jb+npts-1;
         if( s == 3 )
            jb = je-npts+1;
         if( s == 4 )
         {
            ke = kb+npts-1;

#if defined(SW4_USE_RAJA) // SW4 backend
if (ke > AlphaVE_in.m_ke) ke = AlphaVE_in.m_ke;

#else // SW4 backend
if( ke > AlphaVE.m_ke )
               ke = AlphaVE.m_ke;

#endif // SW4 backend
}
         if( s == 5 )
         {
            kb = ke-npts+1;

#if defined(SW4_USE_RAJA) // SW4 backend
if (kb < AlphaVE_in.m_kb) kb = AlphaVE_in.m_kb;

#else // SW4 backend
if( kb < AlphaVE.m_kb )
               kb = AlphaVE.m_kb;

#endif // SW4 backend
}

#if defined(SW4_USE_RAJA) // SW4 backend
auto lm_omega = m_omega;
      auto lm_phase = m_phase;
      auto lm_c = m_c;
      RAJA::RangeSegment k_range(kb, ke + 1);
      RAJA::RangeSegment j_range(jb, je + 1);
      RAJA::RangeSegment i_range(ib, ie + 1);
      RAJA::kernel<TGU_POL_ASYNC>(RAJA::make_tuple(k_range, j_range, i_range),
                                  [=] RAJA_DEVICE(int k, int j, int i) {
      // for (int k = kb; k <= ke; k++)
       //  for (int j = jb; j <= je; j++)
       //    for (int i = ib; i <= ie; i++) {

#else // SW4 backend
for( int k=kb ; k <= ke ; k++ )
            for( int j=jb ; j <= je ; j++ )
               for( int i=ib ; i <= ie ; i++ )
               {

#endif // SW4 backend
AlphaVE(1,i,j,k) =
#if defined(SW4_USE_RAJA) // SW4 backend
cos(lm_omega * (x(i, j, k) - lm_c * t) + lm_phase) *
                sin(lm_omega * x(i, j, k) + lm_phase) *
                cos(lm_omega * (z(i, j, k) - lm_c * t) + lm_phase);


#else // SW4 backend
cos(m_omega*(x(i,j,k)-m_c*t)+m_phase)*
                                     sin(m_omega*x(i,j,k)        +m_phase)*
                                     cos(m_omega*(z(i,j,k)-m_c*t)+m_phase);


#endif // SW4 backend
AlphaVE(2,i,j,k) =
#if defined(SW4_USE_RAJA) // SW4 backend
sin(lm_omega * (x(i, j, k) - lm_c * t)) *
                cos(lm_omega * (y(i, j, k) - lm_c * t) + lm_phase) *
                cos(lm_omega * z(i, j, k) + lm_phase);


#else // SW4 backend
sin(m_omega*(x(i,j,k)-m_c*t)        )*
                                     cos(m_omega*(y(i,j,k)-m_c*t)+m_phase)*
                                     cos(m_omega*z(i,j,k)        +m_phase);


#endif // SW4 backend
AlphaVE(3,i,j,k) =
#if defined(SW4_USE_RAJA) // SW4 backend
cos(lm_omega * x(i, j, k) + lm_phase) *
                cos(lm_omega * y(i, j, k) + lm_phase) *
                sin(lm_omega * (z(i, j, k) - lm_c * t) + lm_phase);
                                  });

#else // SW4 backend
cos(m_omega*x(i,j,k)        +m_phase)*
                                     cos(m_omega*y(i,j,k)        +m_phase)*
                                     sin(m_omega*(z(i,j,k)-m_c*t)+m_phase);

#endif // SW4 backend
}
      }
#if defined(SW4_USE_RAJA) // SW4 backend
#else // SW4 backend
}
#endif // SW4 backend
