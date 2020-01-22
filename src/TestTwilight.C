#include "TestTwilight.h"

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

void TestTwilight::get_ubnd( Sarray& u, Sarray& x, Sarray& y, Sarray& z, float_sw4 t, int npts, int sides[6] )
{
   for( int s=0 ; s < 5 ; s++ )
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
            ke = kb+npts-1;
         if( s == 5 )
            kb = ke-npts+1;
         for( int k=kb ; k <= ke ; k++ )
            for( int j=jb ; j <= je ; j++ )
               for( int i=ib ; i <= ie ; i++ )
               {
                  u(1,i,j,k) = sin(m_omega*(x(i,j,k)-m_c*t))*sin(m_omega*y(i,j,k)+m_phase)*sin(m_omega*z(i,j,k)+m_phase);
                  u(2,i,j,k) = sin(m_omega*x(i,j,k)+m_phase)*sin(m_omega*(y(i,j,k)-m_c*t))*sin(m_omega*z(i,j,k)+m_phase);
                  u(3,i,j,k) = sin(m_omega*x(i,j,k)+m_phase)*sin(m_omega*y(i,j,k)+m_phase)*sin(m_omega*(z(i,j,k)-m_c*t));
               }
      }
}
