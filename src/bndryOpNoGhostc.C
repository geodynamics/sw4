#include "sw4.h"
void varcoeffs4(float_sw4* m_acof, float_sw4* m_ghcof);
void bndryOpNoGhostc(float_sw4* acof_no_gp, float_sw4* ghcof_no_gp,
                     float_sw4* sbop_no_gp) {
  // acofs(i,j,k) is coefficient of a(k) in stencil coefficient (i,j)
  // ghcof is coefficient of ghost point, a(1)*ghcof*u(0) in stencil at i=1.
#define acof(q, k, m) acof_no_gp[q - 1 + 6 * (k - 1) + 48 * (m - 1)]
  varcoeffs4(acof_no_gp, ghcof_no_gp);
  float_sw4 d5[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  // modified coefficients for d(a(x) du): NOT using the ghost point
  //  d5 = 0;
  // Use 5th divided difference to cancel the ghost point contribution
  d5[0] = -1.0;
  d5[1] = 5.0;
  d5[2] = -10.0;
  d5[3] = 10.0;
  d5[4] = -5.0;
  d5[5] = 1.0;
  float_sw4 w0 = 17.0 / 48.0;
  for (int j = 1; j <= 8; j++) acof(1, j, 1) = acof(1, j, 1) + d5[j] / (4 * w0);

  // the coeff for all ghost points are zero (don't divided by them!)
  for (int k = 0; k < 6; k++) ghcof_no_gp[k] = 0;

  // boundary normal derivative, not using ghost points
  // sb = (-25*f(1)/12 + 4*f(2) - 3*f(3) + 4*f(4)/3 - f(5)/4)/h(q);
  sbop_no_gp[0] = 0;
  sbop_no_gp[1] = -25.0 / 12.0;
  sbop_no_gp[2] = 4.0;
  sbop_no_gp[3] = -3.0;
  sbop_no_gp[4] = 4.0 / 3.0;
  sbop_no_gp[5] = -1.0 / 4.0;

#undef acof
}
