#include <stdlib.h>

#include "TestEcons.h"
#include "caliper.h"
TestEcons::TestEcons(float_sw4 amp, float_sw4 cpocs)
    : m_amp(amp), m_cpocs(cpocs) {}

void TestEcons::get_rhobnd(Sarray& rho, int npts, int sides[6]) {
  SW4_MARK_FUNCTION;
  // Fillin random numbers at ghost points.
  for (int s = 0; s < 5; s++)
    if (sides[s] == 1) {
      int kb = rho.m_kb, ke = rho.m_ke, jb = rho.m_jb, je = rho.m_je,
          ib = rho.m_ib, ie = rho.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;
      if (s == 4) ke = kb + npts - 1;
      if (s == 5) kb = ke - npts + 1;
      for (int k = kb; k <= ke; k++)
        for (int j = jb; j <= je; j++)
          for (int i = ib; i <= ie; i++) {
            rho(i, j, k) = m_amp * drand48() + 2;
          }
    }
}

void TestEcons::get_mulabnd(Sarray& mu, Sarray& lambda, int npts,
                            int sides[6]) {
  SW4_MARK_FUNCTION;
  // Fillin random numbers at ghost points.
  for (int s = 0; s < 5; s++)
    if (sides[s] == 1) {
      int kb = mu.m_kb, ke = mu.m_ke, jb = mu.m_jb, je = mu.m_je, ib = mu.m_ib,
          ie = mu.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;
      if (s == 4) ke = kb + npts - 1;
      if (s == 5) kb = ke - npts + 1;
      for (int k = kb; k <= ke; k++)
        for (int j = jb; j <= je; j++)
          for (int i = ib; i <= ie; i++) {
            mu(i, j, k) = m_amp * drand48() + 2;
            lambda(i, j, k) =
                mu(i, j, k) * (m_cpocs * m_cpocs - 2) + m_amp * drand48();
          }
    }
}

void TestEcons::get_ubnd(Sarray& u, int npts, int sides[6]) {
  SW4_MARK_FUNCTION;
  // Homogeneous Dirichet at boundaries

  SView& uV = u.getview();
  for (int s = 0; s < 6; s++)
    if (sides[s] == 1) {
      int kb = u.m_kb, ke = u.m_ke, jb = u.m_jb, je = u.m_je, ib = u.m_ib,
          ie = u.m_ie;
      if (s == 0) ie = ib + npts - 1;
      if (s == 1) ib = ie - npts + 1;
      if (s == 2) je = jb + npts - 1;
      if (s == 3) jb = je - npts + 1;
      if (s == 4) ke = kb + npts - 1;
      if (s == 5) kb = ke - npts + 1;
      // for (int k = kb; k <= ke; k++)
      //   for (int j = jb; j <= je; j++)
      //     for (int i = ib; i <= ie; i++) {

      RAJA::RangeSegment k_range(kb, ke + 1);
      RAJA::RangeSegment j_range(jb, je + 1);
      RAJA::RangeSegment i_range(ib, ie + 1);
      RAJA::kernel<TGU_POL_ASYNC>(RAJA::make_tuple(k_range, j_range, i_range),
                                  [=] RAJA_DEVICE(int k, int j, int i) {
                                    uV(1, i, j, k) = 0;
                                    uV(2, i, j, k) = 0;
                                    uV(3, i, j, k) = 0;
                                  });
    }
  SYNC_STREAM;
}
