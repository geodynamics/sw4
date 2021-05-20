#include <cmath>
#include <iostream>
#include "TestPointSource.h"

#include "Source.h"

void TestPointSource::ubnd(Sarray& u, Sarray& xx, Sarray& yy, Sarray& zz,
                           float_sw4 t, float_sw4 h, int npts, int sides[6]) {
  Source& source = *m_source_ptr;
  timeDep tD;
  if (!(source.getName() == "SmoothWave" ||
        source.getName() == "VerySmoothBump" ||
        source.getName() == "C6SmoothBump" || source.getName() == "Gaussian")) {
    std::cout << "TestPointSource::ubnd: Error, time dependency must be "
                 "SmoothWave, VerySmoothBump, C6SmoothBump, or Gaussian, not "
              << source.getName() << std::endl;
    return;
  } else if (source.getName() == "SmoothWave")
    tD = iSmoothWave;
  else if (source.getName() == "VerySmoothBump")
    tD = iVerySmoothBump;
  else if (source.getName() == "C6SmoothBump")
    tD = iC6SmoothBump;
  else
    tD = iGaussian;

  float_sw4 alpha = m_cp;
  float_sw4 beta = m_cs;
  float_sw4 rho = m_rho;

  float_sw4 x0 = source.getX0();
  float_sw4 y0 = source.getY0();
  float_sw4 z0 = source.getZ0();
  float_sw4 fr = source.getFrequency();
  float_sw4 time = (t - source.getOffset()) * source.getFrequency();
  if (tD == iGaussian) {
    fr = 1 / fr;
    time = time * fr;
  }
  bool ismomentsource = source.isMomentSource();
  float_sw4 fx, fy, fz;
  float_sw4 mxx, myy, mzz, mxy, mxz, myz, m0;

  if (!ismomentsource) {
    source.getForces(fx, fy, fz);
  } else {
    source.getMoments(mxx, mxy, mxz, myy, myz, mzz);
    m0 = 1;
  }
  //   float_sw4 h   = mGridSize[g];
  float_sw4 eps = 1e-3 * h;

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
            float_sw4 x = xx(i, j, k), y = yy(i, j, k), z = zz(i, j, k);
            if (!ismomentsource) {
              float_sw4 R = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) +
                                 (z - z0) * (z - z0));
              if (R < eps)
                u(1, i, j, k) = u(2, i, j, k) = u(3, i, j, k) = 0;
              else {
                float_sw4 A, B;
                if (tD == iSmoothWave) {
                  A = (1 / pow(alpha, 2) * SmoothWave(time, fr * R, alpha) -
                       1 / pow(beta, 2) * SmoothWave(time, fr * R, beta) +
                       3 / pow(fr * R, 2) *
                           SmoothWave_x_T_Integral(time, fr * R, alpha, beta)) /
                      (4 * M_PI * rho * R * R * R);

                  B = (1 / pow(beta, 2) * SmoothWave(time, fr * R, beta) -
                       1 / pow(fr * R, 2) *
                           SmoothWave_x_T_Integral(time, fr * R, alpha, beta)) /
                      (4 * M_PI * rho * R);
                } else if (tD == iVerySmoothBump) {
                  A = (1 / pow(alpha, 2) * VerySmoothBump(time, fr * R, alpha) -
                       1 / pow(beta, 2) * VerySmoothBump(time, fr * R, beta) +
                       3 / pow(fr * R, 2) *
                           VerySmoothBump_x_T_Integral(time, fr * R, alpha,
                                                       beta)) /
                      (4 * M_PI * rho * R * R * R);

                  B = (1 / pow(beta, 2) * VerySmoothBump(time, fr * R, beta) -
                       1 / pow(fr * R, 2) *
                           VerySmoothBump_x_T_Integral(time, fr * R, alpha,
                                                       beta)) /
                      (4 * M_PI * rho * R);
                } else if (tD == iC6SmoothBump) {
                  A = (1 / pow(alpha, 2) * C6SmoothBump(time, fr * R, alpha) -
                       1 / pow(beta, 2) * C6SmoothBump(time, fr * R, beta) +
                       3 / pow(fr * R, 2) *
                           C6SmoothBump_x_T_Integral(time, fr * R, alpha,
                                                     beta)) /
                      (4 * M_PI * rho * R * R * R);

                  B = (1 / pow(beta, 2) * C6SmoothBump(time, fr * R, beta) -
                       1 / pow(fr * R, 2) *
                           C6SmoothBump_x_T_Integral(time, fr * R, alpha,
                                                     beta)) /
                      (4 * M_PI * rho * R);
                } else if (tD == iGaussian) {
                  A = (1 / pow(alpha, 2) * Gaussian(time, R, alpha, fr) -
                       1 / pow(beta, 2) * Gaussian(time, R, beta, fr) +
                       3 / pow(R, 2) *
                           Gaussian_x_T_Integral(time, R, fr, alpha, beta)) /
                      (4 * M_PI * rho * R * R * R);

                  B = (1 / pow(beta, 2) * Gaussian(time, R, beta, fr) -
                       1 / pow(R, 2) *
                           Gaussian_x_T_Integral(time, R, fr, alpha, beta)) /
                      (4 * M_PI * rho * R);
                }
                u(1, i, j, k) =
                    ((x - x0) * (x - x0) * fx + (x - x0) * (y - y0) * fy +
                     (x - x0) * (z - z0) * fz) *
                        A +
                    fx * B;
                u(2, i, j, k) =
                    ((y - y0) * (x - x0) * fx + (y - y0) * (y - y0) * fy +
                     (y - y0) * (z - z0) * fz) *
                        A +
                    fy * B;
                u(3, i, j, k) =
                    ((z - z0) * (x - x0) * fx + (z - z0) * (y - y0) * fy +
                     (z - z0) * (z - z0) * fz) *
                        A +
                    fz * B;
              }
            } else {
              u(1, i, j, k) = u(2, i, j, k) = u(3, i, j, k) = 0;
              // Here, ismomentsource == true
              float_sw4 R = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) +
                                 (z - z0) * (z - z0));
              if (R < eps) {
                u(1, i, j, k) = u(2, i, j, k) = u(3, i, j, k) = 0;
              } else {
                float_sw4 A, B, C, D, E;
                if (tD == iSmoothWave) {
                  A = SmoothWave(time, R, alpha);
                  B = SmoothWave(time, R, beta);
                  C = SmoothWave_x_T_Integral(time, R, alpha, beta);
                  D = d_SmoothWave_dt(time, R, alpha) / pow(alpha, 3) / R;
                  E = d_SmoothWave_dt(time, R, beta) / pow(beta, 3) / R;
                } else if (tD == iVerySmoothBump) {
                  A = VerySmoothBump(time, R, alpha);
                  B = VerySmoothBump(time, R, beta);
                  C = VerySmoothBump_x_T_Integral(time, R, alpha, beta);
                  D = d_VerySmoothBump_dt(time, R, alpha) / pow(alpha, 3) / R;
                  E = d_VerySmoothBump_dt(time, R, beta) / pow(beta, 3) / R;
                } else if (tD == iC6SmoothBump) {
                  A = C6SmoothBump(time, R, alpha);
                  B = C6SmoothBump(time, R, beta);
                  C = C6SmoothBump_x_T_Integral(time, R, alpha, beta);
                  D = d_C6SmoothBump_dt(time, R, alpha) / pow(alpha, 3) / R;
                  E = d_C6SmoothBump_dt(time, R, beta) / pow(beta, 3) / R;
                } else if (tD == iGaussian) {
                  A = Gaussian(time, R, alpha, fr);
                  B = Gaussian(time, R, beta, fr);
                  C = Gaussian_x_T_Integral(time, R, fr, alpha, beta);
                  D = d_Gaussian_dt(time, R, alpha, fr) / pow(alpha, 3) / R;
                  E = d_Gaussian_dt(time, R, beta, fr) / pow(beta, 3) / R;
                }
                u(1, i, j, k) +=
                    // m_xx*G_xx,x
                    +m0 * mxx / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (x - x0) * (x - x0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - 2 * (x - x0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (x - x0) * (x - x0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (x - x0) * (x - x0) / pow(R, 7) -
                        6 * (x - x0) / pow(R, 5)) *
                           C

                     + (x - x0) * (x - x0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     - 1 / pow(R, 3) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     - 3 * (x - x0) / pow(R, 5) * C

                     + (x - x0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (x - x0) * E);
                u(1, i, j, k) +=
                    // m_yy*G_xy,y
                    +m0 * myy / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (x - x0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (y - y0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     + 3 * (x - x0) * (y - y0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (y - y0) * (y - y0) / pow(R, 7) -
                        3 * (x - x0) / pow(R, 5)) *
                           C);
                u(1, i, j, k) +=
                    // m_zz*G_xz,z
                    +m0 * mzz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (z - z0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (x - x0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (z - z0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     + 3 * (x - x0) * (z - z0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (z - z0) * (z - z0) / pow(R, 7) -
                        3 * (x - x0) / pow(R, 5)) *
                           C);
                u(1, i, j, k) +=
                    // m_xy*G_xy,x
                    +m0 * mxy / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (x - x0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (y - y0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (y - y0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     + 3 * (x - x0) * (y - y0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (x - x0) * (y - y0) / pow(R, 7) -
                        3 * (y - y0) / pow(R, 5)) *
                           C);
                u(1, i, j, k) +=
                    // m_xy*G_xx,y
                    +m0 * mxy / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (x - x0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (x - x0) * (x - x0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (x - x0) * (y - y0) / pow(R, 7) * C

                     + (x - x0) * (x - x0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     - 1 / pow(R, 3) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     - 3 * (y - y0) / pow(R, 5) * C

                     + (y - y0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (y - y0) * E);
                u(1, i, j, k) +=
                    // m_xz*G_xz,x
                    +m0 * mxz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (x - x0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (z - z0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (z - z0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     + 3 * (x - x0) * (z - z0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (x - x0) * (z - z0) / pow(R, 7) -
                        3 * (z - z0) / pow(R, 5)) *
                           C);
                u(1, i, j, k) +=
                    // m_yz*G_xz,y
                    +m0 * myz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (z - z0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     + 3 * (x - x0) * (z - z0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (y - y0) * (z - z0) / pow(R, 7) * C);
                u(1, i, j, k) +=
                    // m_xz*G_xx,z
                    +m0 * mxz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (x - x0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (x - x0) * (x - x0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (x - x0) * (z - z0) / pow(R, 7) * C

                     + (x - x0) * (x - x0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     - 1 / pow(R, 3) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     - 3 * (z - z0) / pow(R, 5) * C

                     + (z - z0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (z - z0) * E);
                u(1, i, j, k) +=
                    // m_yz*G_yx,z
                    +m0 * myz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (y - y0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     + 3 * (x - x0) * (y - y0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (y - y0) * (z - z0) / pow(R, 7) * C);
                //------------------------------------------------------------
                u(2, i, j, k) +=
                    // m_xx*G_xy,x
                    m0 * mxx / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (x - x0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (y - y0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (y - y0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     + 3 * (x - x0) * (y - y0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (x - x0) * (y - y0) / pow(R, 7) -
                        3 * (y - y0) / pow(R, 5)) *
                           C);
                u(2, i, j, k) +=
                    // m_yy**G_yy,y
                    +m0 * myy / (4 * M_PI * rho) *
                    (+3 * (y - y0) * (y - y0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - 2 * (y - y0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (y - y0) * (y - y0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + (15 * (y - y0) * (y - y0) * (y - y0) / pow(R, 7) -
                        6 * (y - y0) / pow(R, 5)) *
                           C

                     + (y - y0) * (y - y0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     - 1 / pow(R, 3) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     - 3 * (y - y0) / pow(R, 5) * C

                     + (y - y0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (y - y0) * E);
                u(2, i, j, k) +=
                    // m_zz*G_zy,z
                    +m0 * mzz / (4 * M_PI * rho) *
                    (+3 * (z - z0) * (z - z0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (y - y0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (z - z0) * (y - y0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     + 3 * (z - z0) * (y - y0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + (15 * (z - z0) * (z - z0) * (y - y0) / pow(R, 7) -
                        3 * (y - y0) / pow(R, 5)) *
                           C);
                u(2, i, j, k) +=
                    // m_xy*G_yy,x
                    +m0 * mxy / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (y - y0) * (y - y0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (y - y0) * (y - y0) / pow(R, 7) * C

                     + (y - y0) * (y - y0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     - 1 / pow(R, 3) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     - 3 * (x - x0) / pow(R, 5) * C

                     + (x - x0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (x - x0) * E);
                u(2, i, j, k) +=
                    // m_xz*G_zy,x
                    +m0 * mxz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (y - y0) * (z - z0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     + 3 * (y - y0) * (z - z0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (y - y0) * (z - z0) / pow(R, 7) * C);
                u(2, i, j, k) +=
                    // m_xy*G_xy,y
                    +m0 * mxy / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (x - x0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (y - y0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     + 3 * (x - x0) * (y - y0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (y - y0) * (y - y0) / pow(R, 7) -
                        3 * (x - x0) / pow(R, 5)) *
                           C);
                u(2, i, j, k) +=
                    // m_yz*G_zy,y
                    +m0 * myz / (4 * M_PI * rho) *
                    (+3 * (z - z0) * (y - y0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (z - z0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (z - z0) * (y - y0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     + 3 * (z - z0) * (y - y0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + (15 * (z - z0) * (y - y0) * (y - y0) / pow(R, 7) -
                        3 * (z - z0) / pow(R, 5)) *
                           C);
                u(2, i, j, k) +=
                    // m_xz*G_xy,z
                    +m0 * mxz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (y - y0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     + 3 * (x - x0) * (y - y0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (y - y0) * (z - z0) / pow(R, 7) * C);
                u(2, i, j, k) +=
                    // m_yz*G_yy,z
                    +m0 * myz / (4 * M_PI * rho) *
                    (+3 * (z - z0) * (y - y0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (y - y0) * (y - y0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + 15 * (z - z0) * (y - y0) * (y - y0) / pow(R, 7) * C

                     + (y - y0) * (y - y0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     - 1 / pow(R, 3) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     - 3 * (z - z0) / pow(R, 5) * C

                     + (z - z0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (z - z0) * E);
                //------------------------------------------------------------
                u(3, i, j, k) +=
                    // m_xx*G_zx,x
                    +m0 * mxx / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (x - x0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (z - z0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (z - z0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     + 3 * (x - x0) * (z - z0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (x - x0) * (z - z0) / pow(R, 7) -
                        3 * (z - z0) / pow(R, 5)) *
                           C);
                u(3, i, j, k) +=
                    // m_yy*G_zy,y
                    +m0 * myy / (4 * M_PI * rho) *
                    (+3 * (y - y0) * (y - y0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (z - z0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (y - y0) * (z - z0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     + 3 * (y - y0) * (z - z0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + (15 * (y - y0) * (y - y0) * (z - z0) / pow(R, 7) -
                        3 * (z - z0) / pow(R, 5)) *
                           C);
                u(3, i, j, k) +=
                    // m_zz**G_zz,z
                    +m0 * mzz / (4 * M_PI * rho) *
                    (+3 * (z - z0) * (z - z0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - 2 * (z - z0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (z - z0) * (z - z0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + (15 * (z - z0) * (z - z0) * (z - z0) / pow(R, 7) -
                        6 * (z - z0) / pow(R, 5)) *
                           C

                     + (z - z0) * (z - z0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     - 1 / pow(R, 3) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     - 3 * (z - z0) / pow(R, 5) * C

                     + (z - z0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (z - z0) * E);
                u(3, i, j, k) +=
                    // m_xy*G_zy,x
                    +m0 * mxy / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (y - y0) * (z - z0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     + 3 * (y - y0) * (z - z0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (y - y0) * (z - z0) / pow(R, 7) * C);
                u(3, i, j, k) +=
                    // m_xz**G_zz,x
                    +m0 * mxz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (z - z0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (z - z0) * (z - z0) / pow(R, 5) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (z - z0) * (z - z0) / pow(R, 7) * C

                     + (z - z0) * (z - z0) / pow(R, 3) *
                           ((x - x0) * D - (x - x0) * E)

                     - 1 / pow(R, 3) *
                           ((x - x0) * A / pow(alpha, 2) -
                            (x - x0) * B / pow(beta, 2))

                     - 3 * (x - x0) / pow(R, 5) * C

                     + (x - x0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (x - x0) * E);
                u(3, i, j, k) +=
                    // m_xy*G_xz,y
                    +m0 * mxy / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (y - y0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (z - z0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     + 3 * (x - x0) * (z - z0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + 15 * (x - x0) * (y - y0) * (z - z0) / pow(R, 7) * C);
                u(3, i, j, k) +=
                    // m_yz*G_zz,y
                    +m0 * myz / (4 * M_PI * rho) *
                    (+3 * (y - y0) * (z - z0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     + 3 * (z - z0) * (z - z0) / pow(R, 5) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     + 15 * (y - y0) * (z - z0) * (z - z0) / pow(R, 7) * C

                     + (z - z0) * (z - z0) / pow(R, 3) *
                           ((y - y0) * D - (y - y0) * E)

                     - 1 / pow(R, 3) *
                           ((y - y0) * A / pow(alpha, 2) -
                            (y - y0) * B / pow(beta, 2))

                     - 3 * (y - y0) / pow(R, 5) * C

                     + (y - y0) / (pow(R, 3) * pow(beta, 2)) * B

                     + 1 / R * (y - y0) * E);
                u(3, i, j, k) +=
                    // m_xz*G_xz,z
                    +m0 * mxz / (4 * M_PI * rho) *
                    (+3 * (x - x0) * (z - z0) * (z - z0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (x - x0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (x - x0) * (z - z0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     + 3 * (x - x0) * (z - z0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + (15 * (x - x0) * (z - z0) * (z - z0) / pow(R, 7) -
                        3 * (x - x0) / pow(R, 5)) *
                           C);
                u(3, i, j, k) +=
                    // m_yz*G_yz,z
                    +m0 * myz / (4 * M_PI * rho) *
                    (+3 * (z - z0) * (z - z0) * (y - y0) / pow(R, 5) *
                         (A / pow(alpha, 2) - B / pow(beta, 2))

                     - (y - y0) / pow(R, 3) *
                           (A / pow(alpha, 2) - B / pow(beta, 2))

                     + (z - z0) * (y - y0) / pow(R, 3) *
                           ((z - z0) * D - (z - z0) * E)

                     + 3 * (z - z0) * (y - y0) / pow(R, 5) *
                           ((z - z0) * A / pow(alpha, 2) -
                            (z - z0) * B / pow(beta, 2))

                     + (15 * (z - z0) * (z - z0) * (y - y0) / pow(R, 7) -
                        3 * (y - y0) / pow(R, 5)) *
                           C);
              }
            }
          }
    }
}

//-----------------------------------------------------------------------
// smooth wave for time dependence to test point force term with
float_sw4 TestPointSource::SmoothWave(float_sw4 t, float_sw4 R, float_sw4 c) {
  float_sw4 temp = R;
  float_sw4 c0 = 2187. / 8., c1 = -10935. / 8., c2 = 19683. / 8.,
            c3 = -15309. / 8., c4 = 2187. / 4.;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1,
  //  (c0*pow(t-R/c,3)+c1*pow(t-R/c,4)+c2*pow(t-R/c,5)+c3*pow(t-R/c,6)+c4*pow(t-R/c,7)),
  //  0);
  if ((t - R / c) > 0 && (t - R / c) < 1)
    temp = (c0 * pow(t - R / c, 3) + c1 * pow(t - R / c, 4) +
            c2 * pow(t - R / c, 5) + c3 * pow(t - R / c, 6) +
            c4 * pow(t - R / c, 7));
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// very smooth bump for time dependence for further testing of point force
float_sw4 TestPointSource::VerySmoothBump(float_sw4 t, float_sw4 R,
                                          float_sw4 c) {
  float_sw4 temp = R;
  float_sw4 c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120,
            c5 = -1024;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1,
  //  (c0*pow(t-R/c,5)+c1*pow(t-R/c,6)+c2*pow(t-R/c,7)+c3*pow(t-R/c,8)+c4*pow(t-R/c,9)+c5*pow(t-R/c,10)),
  //  0);
  if ((t - R / c) > 0 && (t - R / c) < 1)
    temp = (c0 * pow(t - R / c, 5) + c1 * pow(t - R / c, 6) +
            c2 * pow(t - R / c, 7) + c3 * pow(t - R / c, 8) +
            c4 * pow(t - R / c, 9) + c5 * pow(t - R / c, 10));
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// C6 smooth bump for time dependence for further testing of point force
float_sw4 TestPointSource::C6SmoothBump(float_sw4 t, float_sw4 R, float_sw4 c) {
  float_sw4 retval = 0;
  if ((t - R / c) > 0 && (t - R / c) < 1)
    retval = 51480.0 * pow((t - R / c) * (1 - t + R / c), 7);
  return retval;
}

//-----------------------------------------------------------------------
// derivative of smooth wave
float_sw4 TestPointSource::d_SmoothWave_dt(float_sw4 t, float_sw4 R,
                                           float_sw4 c) {
  float_sw4 temp = R;
  float_sw4 c0 = 2187. / 8., c1 = -10935. / 8., c2 = 19683. / 8.,
            c3 = -15309. / 8., c4 = 2187. / 4.;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1,
  //  (3*c0*pow(t-R/c,2)+4*c1*pow(t-R/c,3)+5*c2*pow(t-R/c,4)+6*c3*pow(t-R/c,5)+7*c4*pow(t-R/c,6)),
  //  0);
  if ((t - R / c) > 0 && (t - R / c) < 1)
    temp = (3 * c0 * pow(t - R / c, 2) + 4 * c1 * pow(t - R / c, 3) +
            5 * c2 * pow(t - R / c, 4) + 6 * c3 * pow(t - R / c, 5) +
            7 * c4 * pow(t - R / c, 6));
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// very smooth bump for time dependence to further testing of point force
float_sw4 TestPointSource::d_VerySmoothBump_dt(float_sw4 t, float_sw4 R,
                                               float_sw4 c) {
  float_sw4 temp = R;
  float_sw4 c0 = 1024, c1 = -5120, c2 = 10240, c3 = -10240, c4 = 5120,
            c5 = -1024;

  //  temp = where ( (t-R/c) > 0 && (t-R/c) < 1,
  //  (5*c0*pow(t-R/c,4)+6*c1*pow(t-R/c,5)+7*c2*pow(t-R/c,6)+8*c3*pow(t-R/c,7)+9*c4*pow(t-R/c,8))+10*c5*pow(t-R/c,9),
  //  0);
  if ((t - R / c) > 0 && (t - R / c) < 1)
    temp = (5 * c0 * pow(t - R / c, 4) + 6 * c1 * pow(t - R / c, 5) +
            7 * c2 * pow(t - R / c, 6) + 8 * c3 * pow(t - R / c, 7) +
            9 * c4 * pow(t - R / c, 8)) +
           10 * c5 * pow(t - R / c, 9);
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// C6 smooth bump for time dependence to further testing of point force
float_sw4 TestPointSource::d_C6SmoothBump_dt(float_sw4 t, float_sw4 R,
                                             float_sw4 c) {
  float_sw4 retval = 0;
  if ((t - R / c) > 0 && (t - R / c) < 1)
    retval = 51480.0 * 7 * (1 - 2 * (t - R / c)) *
             pow((t - R / c) * (1 - t + R / c), 6);
  return retval;
}

//-----------------------------------------------------------------------
// Primitive function (for T) of SmoothWave(t-T)*T
float_sw4 TestPointSource::SWTP(float_sw4 Lim, float_sw4 t) {
  float_sw4 temp = Lim;

  float_sw4 c0 = 2187. / 8., c1 = -10935. / 8., c2 = 19683. / 8.,
            c3 = -15309. / 8., c4 = 2187. / 4.;

  temp = (pow(t, 3) *
          (c0 + c1 * t + c2 * pow(t, 2) + c3 * pow(t, 3) + c4 * pow(t, 4)) *
          pow(Lim, 2)) /
             2. -
         (pow(t, 2) *
          (3 * c0 + 4 * c1 * t + 5 * c2 * pow(t, 2) + 6 * c3 * pow(t, 3) +
           7 * c4 * pow(t, 4)) *
          pow(Lim, 3)) /
             3. +
         (t *
          (3 * c0 + 6 * c1 * t + 10 * c2 * pow(t, 2) + 15 * c3 * pow(t, 3) +
           21 * c4 * pow(t, 4)) *
          pow(Lim, 4)) /
             4. +
         ((-c0 - 4 * c1 * t - 10 * c2 * pow(t, 2) - 20 * c3 * pow(t, 3) -
           35 * c4 * pow(t, 4)) *
          pow(Lim, 5)) /
             5. +
         ((c1 + 5 * c2 * t + 15 * c3 * pow(t, 2) + 35 * c4 * pow(t, 3)) *
          pow(Lim, 6)) /
             6. +
         ((-c2 - 6 * c3 * t - 21 * c4 * pow(t, 2)) * pow(Lim, 7)) / 7. +
         ((c3 + 7 * c4 * t) * pow(Lim, 8)) / 8. - (c4 * pow(Lim, 9)) / 9.;

  return temp;
}

//-----------------------------------------------------------------------
// Primitive function (for T) of VerySmoothBump(t-T)*T
float_sw4 TestPointSource::VSBTP(float_sw4 Lim, float_sw4 t) {
  float_sw4 temp = Lim;
  float_sw4 f = 1024., g = -5120., h = 10240., i = -10240., j = 5120.,
            k = -1024.;

  temp =
      (pow(Lim, 11) * (-25200 * k * t - 2520 * j) + 2310 * k * pow(Lim, 12) +
       (124740 * k * pow(t, 2) + 24948 * j * t + 2772 * i) * pow(Lim, 10) +
       (-369600 * k * pow(t, 3) - 110880 * j * pow(t, 2) - 24640 * i * t -
        3080 * h) *
           pow(Lim, 9) +
       (727650 * k * pow(t, 4) + 291060 * j * pow(t, 3) +
        97020 * i * pow(t, 2) + 24255 * h * t + 3465 * g) *
           pow(Lim, 8) +
       (-997920 * k * pow(t, 5) - 498960 * j * pow(t, 4) -
        221760 * i * pow(t, 3) - 83160 * h * pow(t, 2) - 23760 * g * t -
        3960 * f) *
           pow(Lim, 7) +
       (970200 * k * pow(t, 6) + 582120 * j * pow(t, 5) +
        323400 * i * pow(t, 4) + 161700 * h * pow(t, 3) +
        69300 * g * pow(t, 2) + 23100 * f * t) *
           pow(Lim, 6) +
       (-665280 * k * pow(t, 7) - 465696 * j * pow(t, 6) -
        310464 * i * pow(t, 5) - 194040 * h * pow(t, 4) -
        110880 * g * pow(t, 3) - 55440 * f * pow(t, 2)) *
           pow(Lim, 5) +
       (311850 * k * pow(t, 8) + 249480 * j * pow(t, 7) +
        194040 * i * pow(t, 6) + 145530 * h * pow(t, 5) +
        103950 * g * pow(t, 4) + 69300 * f * pow(t, 3)) *
           pow(Lim, 4) +
       (-92400 * k * pow(t, 9) - 83160 * j * pow(t, 8) - 73920 * i * pow(t, 7) -
        64680 * h * pow(t, 6) - 55440 * g * pow(t, 5) - 46200 * f * pow(t, 4)) *
           pow(Lim, 3) +
       (13860 * k * pow(t, 10) + 13860 * j * pow(t, 9) + 13860 * i * pow(t, 8) +
        13860 * h * pow(t, 7) + 13860 * g * pow(t, 6) + 13860 * f * pow(t, 5)) *
           pow(Lim, 2)) /
      27720.0;

  return temp;
}
//-----------------------------------------------------------------------
// Primitive function (for T) of C6SmoothBump(t-T)*T
float_sw4 TestPointSource::C6SBTP(float_sw4 Lim, float_sw4 t) {
  float_sw4 x = t - Lim;
  return pow(x, 8) *
         (-3217.5 * pow(x, 8) + 3432.0 * (7 + t) * pow(x, 7) -
          25740.0 * (3 + t) * pow(x, 6) + 27720.0 * (5 + 3 * t) * pow(x, 5) -
          150150.0 * (t + 1) * x * x * x * x +
          32760.0 * (3 + 5 * t) * x * x * x - 36036.0 * (1 + 3 * t) * x * x +
          5720.0 * (1 + 7 * t) * x - 6435.0 * t);
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*SmoothWave(t-T)*T from R/alpha to R/beta
float_sw4 TestPointSource::SmoothWave_x_T_Integral(float_sw4 t, float_sw4 R,
                                                   float_sw4 alpha,
                                                   float_sw4 beta) {
  float_sw4 temp = R;

  float_sw4 lowL, hiL;

  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t,
  //  R / beta, t);
  if ((R / alpha > t - 1))
    lowL = R / alpha;
  else
    lowL = t - 1;
  if (R / beta < t)
    hiL = R / beta;
  else
    hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, SWTP(hiL, t) - SWTP(lowL, t), 0.0);
  if (lowL < t && hiL > t - 1)
    temp = SWTP(hiL, t) - SWTP(lowL, t);
  else
    temp = 0;

  return temp;
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*VerySmoothBump(t-T)*T from R/alpha to R/beta
float_sw4 TestPointSource::VerySmoothBump_x_T_Integral(float_sw4 t, float_sw4 R,
                                                       float_sw4 alpha,
                                                       float_sw4 beta) {
  float_sw4 temp = R;

  float_sw4 lowL, hiL;

  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t,
  //  R / beta, t);
  if (R / alpha > t - 1)
    lowL = R / alpha;
  else
    lowL = t - 1;
  if (R / beta < t)
    hiL = R / beta;
  else
    hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t),
  //  0.0);
  if (lowL < t && hiL > t - 1)
    temp = VSBTP(hiL, t) - VSBTP(lowL, t);
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
// Integral of H(t-T)*H(1-t+T)*C6SmoothBump(t-T)*T from R/alpha to R/beta
float_sw4 TestPointSource::C6SmoothBump_x_T_Integral(float_sw4 t, float_sw4 R,
                                                     float_sw4 alpha,
                                                     float_sw4 beta) {
  float_sw4 temp = R;

  float_sw4 lowL, hiL;

  //  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t,
  //  R / beta, t);
  if (R / alpha > t - 1)
    lowL = R / alpha;
  else
    lowL = t - 1;
  if (R / beta < t)
    hiL = R / beta;
  else
    hiL = t;

  //  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t),
  //  0.0);
  if (lowL < t && hiL > t - 1)
    temp = C6SBTP(hiL, t) - C6SBTP(lowL, t);
  else
    temp = 0;
  return temp;
}

//-----------------------------------------------------------------------
float_sw4 TestPointSource::Gaussian(float_sw4 t, float_sw4 R, float_sw4 c,
                                    float_sw4 f) {
  float_sw4 temp = R;
  temp = 1 / (f * sqrt(2 * M_PI)) * exp(-pow(t - R / c, 2) / (2 * f * f));
  return temp;
}

//-----------------------------------------------------------------------
float_sw4 TestPointSource::d_Gaussian_dt(float_sw4 t, float_sw4 R, float_sw4 c,
                                         float_sw4 f) {
  float_sw4 temp = R;
  temp = 1 / (f * sqrt(2 * M_PI)) *
         (-exp(-pow(t - R / c, 2) / (2 * f * f)) * (t - R / c)) / pow(f, 2);
  return temp;
}

//-----------------------------------------------------------------------
float_sw4 TestPointSource::Gaussian_x_T_Integral(float_sw4 t, float_sw4 R,
                                                 float_sw4 f, float_sw4 alpha,
                                                 float_sw4 beta) {
  float_sw4 temp = R;
  temp = -0.5 * t *
             (erf((t - R / beta) / (sqrt(2.0) * f)) -
              erf((t - R / alpha) / (sqrt(2.0) * f))) -
         f / sqrt(2 * M_PI) *
             (exp(-pow(t - R / beta, 2) / (2 * f * f)) -
              exp(-pow(t - R / alpha, 2) / (2 * f * f)));
  //  temp = 1/(f*sqrt(2*M_PI))*(
  //  f*f*(-exp(-pow(t-R/beta,2)/(2*f*f))+exp(-pow(t-R/alpha,2)/(2*f*f)) ) +
  //	     t*0.5*sqrt(M_PI*2)*f*( erf((t-R/alpha)/(sqrt(2.0)*f)) -
  //erf((t-R/beta)/(sqrt(2.0)*f)) ) );
  //  temp = 1 /(f*sqrt(2*M_PI))*(f*( (-exp(-pow(t-R / alpha,2)/pow(f,2)) +
  //  exp(-pow(t-R / beta,2)/pow(f,2)) )*f + sqrt(M_PI)*t*(-erf((t-R / alpha) /
  //  f) + erf(R / beta / f))))/2.;
  return temp;
}
