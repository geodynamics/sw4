#include "TestTwilight.h"
#include "caliper.h"

TestTwilight::TestTwilight(float_sw4 omega, float_sw4 c, float_sw4 phase,
                           float_sw4 momega, float_sw4 mphase, float_sw4 amprho,
                           float_sw4 ampmu, float_sw4 amplambda)
    : m_omega(omega),
      m_c(c),
      m_phase(phase),
      m_momega(momega),
      m_mphase(mphase),
      m_amprho(amprho),
      m_ampmu(ampmu),
      m_amplambda(amplambda),
      m_sw4twilight(true) {}

void TestTwilight::get_rho(Sarray& rho, Sarray& x, Sarray& y, Sarray& z) {
  SW4_MARK_FUNCTION;
  if (m_sw4twilight) {
    for (int k = rho.m_kb; k <= rho.m_ke; k++)
      for (int j = rho.m_jb; j <= rho.m_je; j++)
        for (int i = rho.m_ib; i <= rho.m_ie; i++)
          rho(i, j, k) =
              m_amprho * (2 + sin(m_momega * x(i, j, k) + m_mphase) *
                                  cos(m_momega * y(i, j, k) + m_mphase) *
                                  sin(m_momega * z(i, j, k) + m_mphase));
  } else {
    // Test code material
    for (int k = rho.m_kb; k <= rho.m_ke; k++)
      for (int j = rho.m_jb; j <= rho.m_je; j++)
        for (int i = rho.m_ib; i <= rho.m_ie; i++)
          rho(i, j, k) = 2 + sin(x(i, j, k) + 0.3) * sin(y(i, j, k) + 0.3) *
                                 sin(z(i, j, k) - 0.2);
  }
}

void TestTwilight::get_mula(Sarray& mu, Sarray& lambda, Sarray& x, Sarray& y,
                            Sarray& z) {
  SW4_MARK_FUNCTION;
  if (m_sw4twilight) {
    for (int k = mu.m_kb; k <= mu.m_ke; k++)
      for (int j = mu.m_jb; j <= mu.m_je; j++)
        for (int i = mu.m_ib; i <= mu.m_ie; i++) {
          mu(i, j, k) =
              m_ampmu * (3 + cos(m_momega * x(i, j, k) + m_mphase) *
                                 sin(m_momega * y(i, j, k) + m_mphase) *
                                 sin(m_momega * z(i, j, k) + m_mphase));
          lambda(i, j, k) =
              m_amplambda * (2 + sin(m_momega * x(i, j, k) + m_mphase) *
                                     sin(m_momega * y(i, j, k) + m_mphase) *
                                     cos(m_momega * z(i, j, k) + m_mphase));
        }
  } else {
    // Test code material
    for (int k = mu.m_kb; k <= mu.m_ke; k++)
      for (int j = mu.m_jb; j <= mu.m_je; j++)
        for (int i = mu.m_ib; i <= mu.m_ie; i++) {
          mu(i, j, k) = 3.0 + sin(3 * x(i, j, k) + 0.1) *
                                  sin(3 * y(i, j, k) + 0.1) * sin(z(i, j, k));
          lambda(i, j, k) = 21.0 + cos(x(i, j, k) + 0.1) *
                                       cos(y(i, j, k) + 0.1) *
                                       pow(sin(3 * z(i, j, k)), 2);
        }
  }
}

void TestTwilight::get_ubnd(Sarray& u_in, Sarray& x_in, Sarray& y_in, Sarray& z_in,
                            float_sw4 t, int npts, int sides[6]) {
  SW4_MARK_FUNCTION;
  
  SView& u = u_in.getview();
  SView& x = x_in.getview();
  SView& y = y_in.getview();
  SView& z = z_in.getview();
  for (int s = 0; s < 6; s++)
    if (sides[s] == 1) {
      int kb = u_in.m_kb, ke = u_in.m_ke, jb = u_in.m_jb, je = u_in.m_je, ib = u_in.m_ib,
          ie = u_in.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;
      if (s == 4) {
        ke = kb + npts - 1;
        if (ke > u_in.m_ke) ke = u_in.m_ke;
      }
      if (s == 5) {
        kb = ke - npts + 1;
        if (kb < u_in.m_kb) kb = u_in.m_kb;
      }
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
    }
  SYNC_STREAM;
}

void TestTwilight::get_ubnd(Sarray& u, float_sw4 h, float_sw4 zmin, float_sw4 t,
                            int npts, int sides[6]) {
  SW4_MARK_FUNCTION;
  SYNC_STREAM;
  for (int s = 0; s < 6; s++)
    if (sides[s] == 1) {
      int kb = u.m_kb, ke = u.m_ke, jb = u.m_jb, je = u.m_je, ib = u.m_ib,
          ie = u.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;
      if (s == 4) {
        ke = kb + npts - 1;
        if (ke > u.m_ke) ke = u.m_ke;
      }
      if (s == 5) {
        kb = ke - npts + 1;
        if (kb < u.m_kb) kb = u.m_kb;
      }
      for (int k = kb; k <= ke; k++)
        for (int j = jb; j <= je; j++)
          for (int i = ib; i <= ie; i++) {
            float_sw4 x = h * (i - 1), y = h * (j - 1), z = h * (k - 1) + zmin;
            u(1, i, j, k) = sin(m_omega * (x - m_c * t)) *
                            sin(m_omega * y + m_phase) *
                            sin(m_omega * z + m_phase);
            u(2, i, j, k) = sin(m_omega * x + m_phase) *
                            sin(m_omega * (y - m_c * t)) *
                            sin(m_omega * z + m_phase);
            u(3, i, j, k) = sin(m_omega * x + m_phase) *
                            sin(m_omega * y + m_phase) *
                            sin(m_omega * (z - m_c * t));
          }
    }
}

void TestTwilight::get_mula_att(Sarray& muve, Sarray& lambdave, Sarray& x,
                                Sarray& y, Sarray& z) {
  SW4_MARK_FUNCTION;
  if (m_sw4twilight) {
    for (int k = muve.m_kb; k <= muve.m_ke; k++)
      for (int j = muve.m_jb; j <= muve.m_je; j++)
        for (int i = muve.m_ib; i <= muve.m_ie; i++) {
          muve(i, j, k) =
              m_ampmu * (1.5 + 0.5 * cos(m_momega * x(i, j, k) + m_mphase) *
                                   cos(m_momega * y(i, j, k) + m_mphase) *
                                   sin(m_momega * z(i, j, k) + m_mphase));
          lambdave(i, j, k) =
              m_amplambda *
              (0.5 + 0.25 * sin(m_momega * x(i, j, k) + m_mphase) *
                         cos(m_momega * y(i, j, k) + m_mphase) *
                         sin(m_momega * z(i, j, k) + m_mphase));
        }
  }
}

//void TestTwilight::get_bnd_att_device(Sarray& AlphaVE_i, Sarray& x_i, Sarray& y_i, Sarray& z_i,
//                                float_sw4 t, int npts, int sides[6]) {
//   SW4_MARK_FUNCTION;
//   //std::cout << "WARNING TestTwilight::get_bnd_att running on CPU\n"
//   //          << std::flush;
//   auto& AlphaVE = AlphaVE_i.getview();
//   auto& x = x_i.getview();
//   auto& y = y_i.getview();
//   auto& z = z_i.getview();
//   for (int s = 0; s < 6; s++)
//     if (sides[s] == 1) {
//       int kb = AlphaVE_i.m_kb, ke = AlphaVE_i.m_ke, jb = AlphaVE_i.m_jb,
//           je = AlphaVE_i.m_je, ib = AlphaVE_i.m_ib, ie = AlphaVE_i.m_ie;
//       if (s == 0) ie = ib + npts - 1;
//       if (s == 1) ib = ie - npts + 1;
//       if (s == 2) je = jb + npts - 1;
//       if (s == 3) jb = je - npts + 1;
//       if (s == 4) {
//         ke = kb + npts - 1;
//         if (ke > AlphaVE_i.m_ke) ke = AlphaVE_i.m_ke;
//       }
//       if (s == 5) {
//         kb = ke - npts + 1;
//         if (kb < AlphaVE_i.m_kb) kb = AlphaVE_i.m_kb;
//       }

//       auto l_omega = m_omega;
//       auto l_phase = m_phase;
      
//       RAJA::RangeSegment k_range(kb, ke + 1);
//       RAJA::RangeSegment j_range(jb, je + 1);
//       RAJA::RangeSegment i_range(ib, ie + 1);
//       RAJA::kernel<DEFAULT_LOOP3>(
//       RAJA::make_tuple(k_range, j_range, i_range),
//       [=] RAJA_DEVICE(int k, int j, int i) {
//       // for (int k = kb; k <= ke; k++)
//       //   for (int j = jb; j <= je; j++)
//       //     for (int i = ib; i <= ie; i++) {
//             AlphaVE(1, i, j, k) =
//                 cos(l_omega * (x(i, j, k) - m_c * t) + l_phase) *
//                 sin(l_omega * x(i, j, k) + l_phase) *
//                 cos(l_omega * (z(i, j, k) - m_c * t) + l_phase);

//             AlphaVE(2, i, j, k) =
//                 sin(l_omega * (x(i, j, k) - m_c * t)) *
//                 cos(l_omega * (y(i, j, k) - m_c * t) + l_phase) *
//                 cos(l_omega * z(i, j, k) + l_phase);

//             AlphaVE(3, i, j, k) =
//                 cos(l_omega * x(i, j, k) + l_phase) *
//                 cos(l_omega * y(i, j, k) + l_phase) *
//                 sin(l_omega * (z(i, j, k) - m_c * t) + l_phase);
//       });
//     }
// }
void TestTwilight::get_bnd_att(Sarray& AlphaVE_in, Sarray& x_in, Sarray& y_in, Sarray& z_in,
                               float_sw4 t, int npts, int sides[6]) {
  SW4_MARK_FUNCTION;
  SView& AlphaVE = AlphaVE_in.getview();
  SView& x = x_in.getview();
  SView& y = y_in.getview();
  SView& z = z_in.getview();
  //std::cout << "WARNING TestTwilight::get_bnd_att running on CPU\n"
  //          << std::flush;
  for (int s = 0; s < 6; s++)
    if (sides[s] == 1) {
      int kb = AlphaVE_in.m_kb, ke = AlphaVE_in.m_ke, jb = AlphaVE_in.m_jb,
          je = AlphaVE_in.m_je, ib = AlphaVE_in.m_ib, ie = AlphaVE_in.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;
      if (s == 4) {
        ke = kb + npts - 1;
        if (ke > AlphaVE_in.m_ke) ke = AlphaVE_in.m_ke;
      }
      if (s == 5) {
        kb = ke - npts + 1;
        if (kb < AlphaVE_in.m_kb) kb = AlphaVE_in.m_kb;
      }
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
            AlphaVE(1, i, j, k) =
                cos(lm_omega * (x(i, j, k) - lm_c * t) + lm_phase) *
                sin(lm_omega * x(i, j, k) + lm_phase) *
                cos(lm_omega * (z(i, j, k) - lm_c * t) + lm_phase);

            AlphaVE(2, i, j, k) =
                sin(lm_omega * (x(i, j, k) - lm_c * t)) *
                cos(lm_omega * (y(i, j, k) - lm_c * t) + lm_phase) *
                cos(lm_omega * z(i, j, k) + lm_phase);

            AlphaVE(3, i, j, k) =
                cos(lm_omega * x(i, j, k) + lm_phase) *
                cos(lm_omega * y(i, j, k) + lm_phase) *
                sin(lm_omega * (z(i, j, k) - lm_c * t) + lm_phase);
				  });
    }
}
